"""One-shot verification service entrypoint.

It (1) waits for the API health endpoint, (2) checks the build by
byte-compiling every source file, (3) runs the unit tests, and (4) runs
the decode smoke test (substitution/insertion/deletion + unique +
ambiguous + constraint-failure paths).  The process exit code reports the
overall result: 0 only if every stage passed.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
import urllib.request

API_URL = os.environ.get("API_URL", "http://127.0.0.1:8080").rstrip("/")
HEALTH_TIMEOUT = float(os.environ.get("HEALTH_TIMEOUT", "60"))
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def wait_for_health() -> bool:
    deadline = time.time() + HEALTH_TIMEOUT
    attempt = 0
    while time.time() < deadline:
        attempt += 1
        try:
            with urllib.request.urlopen(f"{API_URL}/healthz", timeout=3) as r:
                if r.status == 200 and json.load(r).get("status") == "ok":
                    print(f"[verify] API healthy after {attempt} attempt(s)")
                    return True
        except Exception as exc:  # noqa: BLE001 - any failure means retry
            print(f"[verify] waiting for API ({attempt}): {exc}")
        time.sleep(1.0)
    return False


def run_step(name: str, cmd: list[str]) -> bool:
    print(f"\n[verify] === {name} ===")
    print(f"[verify] $ {' '.join(cmd)}")
    completed = subprocess.run(cmd, cwd=ROOT)
    if completed.returncode != 0:
        print(f"[verify] {name} FAILED with exit code {completed.returncode}")
        return False
    print(f"[verify] {name} passed")
    return True


def main() -> int:
    print(f"[verify] root={ROOT} api={API_URL}")
    if not wait_for_health():
        print("[verify] API did not become healthy in time", file=sys.stderr)
        return 1

    steps = [
        ("build check (byte-compile)",
         [sys.executable, "-m", "compileall", "-q", "app", "tests", "scripts"]),
        ("unit tests",
         [sys.executable, "-m", "unittest", "discover", "-s", "tests", "-v"]),
        ("decode smoke (sub/ins/del + unique/ambiguous/failure)",
         [sys.executable, "scripts/smoke.py"]),
    ]
    failed = [name for name, cmd in steps if not run_step(name, cmd)]
    if failed:
        print(f"\n[verify] FAILED stages: {', '.join(failed)}", file=sys.stderr)
        return 1
    print("\n[verify] ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
