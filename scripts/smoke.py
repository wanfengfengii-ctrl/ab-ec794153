"""End-to-end decode smoke test against the running HTTP API.

Exercises the happy path with a read deliberately carrying one
substitution, one insertion and one deletion, plus the unique, ambiguous
and constraint-failure branches.  Exits non-zero on any problem.
"""

from __future__ import annotations

import json
import os
import sys
import urllib.error
import urllib.request

API_URL = os.environ.get("API_URL", "http://127.0.0.1:8080").rstrip("/")

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from app.solver import replay  # noqa: E402


def post(payload: dict):
    req = urllib.request.Request(
        f"{API_URL}/api/concatemers/decode",
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            return resp.status, json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        return exc.code, json.loads(exc.read().decode("utf-8"))


def get(path: str):
    with urllib.request.urlopen(f"{API_URL}{path}", timeout=15) as resp:
        return resp.status, json.loads(resp.read().decode("utf-8"))


def check(cond: bool, message: str) -> None:
    if not cond:
        raise AssertionError(message)
    print(f"  ok - {message}")


def main() -> int:
    print(f"smoke target: {API_URL}")

    status, body = get("/healthz")
    check(status == 200 and body.get("status") == "ok", "health endpoint is ok")

    # --- read carrying a substitution, an insertion and a deletion ------
    ref = "ACGTAGCCTA"  # 10 bases
    rotated = ref[3:] + ref[:3]
    c1 = rotated[:2] + next(b for b in "ACGT" if b != rotated[2]) + rotated[3:]
    c2 = rotated[:6] + "T" + rotated[6:]
    c3 = rotated[:4] + rotated[5:]
    read = c1 + c2 + c3
    check(len(read) == 30, f"crafted read has length 30 (got {len(read)})")

    status, body = post({
        "reference": ref,
        "read": read,
        "copies": 3,
        "max_edits_per_copy": 3,
    })
    check(status == 200, f"indel/sub decode returns 200 (got {status})")
    check(body["status"] in ("unique", "ambiguous"),
          f"decoding status is unique or ambiguous (got {body['status']})")
    witness = body if body["status"] == "unique" else body["witnesses"][0]

    # every emitted CIGAR must replay and reconstruct the whole read
    kinds = set()
    prev = 0
    for seg in witness["segments"]:
        piece = read[prev:seg["end"]]
        trace = replay(piece, witness["rotated_reference"], seg["cigar"])
        check(trace["edits"] == seg["edits"],
              f"segment {prev}:{seg['end']} CIGAR replays with same edits")
        kinds.update(o["op"] for o in trace["operations"] if o["op"] != "=")
        prev = seg["end"]
    check(prev == len(read), "segments cover the whole read")
    check(kinds == {"X", "I", "D"},
          f"optimal alignment contains substitution, insertion and deletion (got {sorted(kinds)})")
    check(witness["total_edits"] <= 9, "total edits within budget")

    # --- unique zero-edit repeat ----------------------------------------
    ref2 = "ACGTGGTCAA"
    rot2 = ref2[4:] + ref2[:4]
    status, body = post({
        "reference": ref2,
        "read": rot2 * 3,
        "copies": 3,
        "max_edits_per_copy": 0,
    })
    check(status == 200 and body["status"] == "unique",
          "clean triple repeat decodes as unique")
    check(body["shift"] == 4 and body["total_edits"] == 0,
          "unique witness reports shift=4 and zero edits")
    check(all(s["cigar"] == "10=" for s in body["segments"]),
          "unique witness CIGARs are three 10= runs")

    # --- ambiguity over a homopolymer -----------------------------------
    status, body = post({
        "reference": "AAAAAAAA",
        "read": "A" * 30,
        "copies": 3,
        "max_edits_per_copy": 3,
    })
    check(status == 200 and body["status"] == "ambiguous",
          "homopolymer alignment is reported ambiguous")
    check(len(body["witnesses"]) == 2, "ambiguous response carries two witnesses")
    keys = [
        (w["shift"], tuple(s["end"] for s in w["segments"]),
         tuple(s["cigar"] for s in w["segments"]))
        for w in body["witnesses"]
    ]
    check(keys == sorted(keys) and keys[0] < keys[1],
          "witnesses follow the stable (shift, boundaries, CIGAR) order")

    # --- locatable constraint failure -----------------------------------
    status, body = post({
        "reference": "ACGTACGTA",
        "read": "A" * 30,
        "copies": 3,
        "max_edits_per_copy": 0,
    })
    check(status == 422 and body.get("error") == "constraint_failure",
          "infeasible length window returns 422 constraint_failure")
    check(body.get("field") == "read"
          and body.get("reason") == "read_length_outside_window",
          "constraint failure is locatable (field=read)")

    # --- invalid request -------------------------------------------------
    status, body = post({
        "reference": "ACGT",          # too short
        "read": "A" * 30,
        "copies": 3,
        "max_edits_per_copy": 1,
    })
    check(status == 422 and body.get("field") == "reference",
          "invalid reference is rejected with a locatable field")

    print("SMOKE OK")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:  # noqa: BLE001 - report any failure via exit code
        print(f"SMOKE FAILED: {exc}", file=sys.stderr)
        sys.exit(1)
