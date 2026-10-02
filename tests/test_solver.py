"""Tests for the joint concatemer decoder.

The strongest checks brute-force every distinct rotation and every
segmentation on small instances, including every optimal per-segment
CIGAR (up to the same cap the solver keeps), and compare against the
dynamic-programming solver.
"""

from __future__ import annotations

import itertools
import random
import unittest

from app.solver import (
    COUNT_CAP,
    ConstraintFailure,
    InvalidRequest,
    decode,
    replay,
    align_all,
    _encode_cigar,
    _needleman,
    _rotations,
)


def levenshtein(a: str, b: str) -> int:
    m, n = len(a), len(b)
    prev = list(range(n + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1,
                           prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[n]


def mutate(seq: str, rng: random.Random, edits: int) -> str:
    out = list(seq)
    budget = edits
    while budget:
        kind = rng.choice(["S", "I", "D"] if len(out) > 1 else ["S", "I"])
        if kind == "S":
            i = rng.randrange(len(out))
            out[i] = rng.choice([b for b in "ACGT" if b != out[i]])
        elif kind == "I":
            i = rng.randrange(len(out) + 1)
            out.insert(i, rng.choice("ACGT"))
        else:
            i = rng.randrange(len(out))
            del out[i]
        budget -= 1
    return "".join(out)


def apply_edit(seq: str, kind: str, where: int, base: str = "T") -> str:
    """Deterministic single edit: sub/ins/del at a fixed position."""
    out = list(seq)
    if kind == "S":
        out[where] = next(b for b in "ACGT" if b != out[where])
    elif kind == "I":
        out.insert(where, base)
    else:
        del out[where]
    return "".join(out)


def brute_force(read: str, ref: str, copies: int, per: int):
    """Enumerate every interpretation the same way the solver defines one.

    Returns ``(best_unconstrained, optimal)`` where ``optimal`` is a list
    of ``(shift, bounds, [cigar_tuple, ...])`` sorted by
    ``(shift, bounds, cigar_tuple)``; at most two CIGAR tuples are stored
    per boundary tuple (sufficient for the first two global witnesses).
    Boundary tuples with several optimal CIGARs are flagged via a count.
    """
    length = len(ref)
    lo, hi = length - per, length + per
    n = len(read)
    best_unconstrained = None
    feasible = []
    if not (copies * lo <= n <= copies * hi):
        return best_unconstrained, []
    for shift, rotated in _rotations(ref):
        for cuts in itertools.combinations(range(1, n), copies - 1):
            bounds = cuts + (n,)
            prev = 0
            pair = (0, 0)
            per_seg = []
            length_ok = True
            for end in bounds:
                if not (lo <= end - prev <= hi):
                    length_ok = False
                    break
                cost, cigars, count = align_all(
                    read[prev:end], rotated, top_k=3
                )
                per_seg.append((cost, cigars, count))
                pair = (pair[0] + cost, max(pair[1], cost))
                prev = end
            if not length_ok:
                continue
            if best_unconstrained is None or pair < best_unconstrained:
                best_unconstrained = pair
            if pair[1] <= per and pair[0] <= copies * per:
                feasible.append((shift, bounds, per_seg, pair))
    if not feasible:
        return best_unconstrained, []
    target = min(w[3] for w in feasible)
    optimal = []
    for shift, bounds, per_seg, pair in feasible:
        if pair != target:
            continue
        lists = [list(ps[1]) for ps in per_seg]
        # first two CIGAR tuples of the sorted cartesian product
        combos = []
        if all(lists):
            combos.append(tuple(lst[0] for lst in lists))
            if any(ps[2] >= COUNT_CAP for ps in per_seg):
                best_second = None
                for d, ps in enumerate(per_seg):
                    if ps[2] >= COUNT_CAP and len(lists[d]) >= 2:
                        cand = tuple(
                            lists[k][1] if k == d else lists[k][0]
                            for k in range(len(lists))
                        )
                        if best_second is None or cand < best_second:
                            best_second = cand
                if best_second is not None:
                    combos.append(best_second)
        optimal.append((shift, bounds, combos, pair))
    optimal.sort(key=lambda w: (w[0], w[1],
                                w[2][0] if w[2] else None))
    return best_unconstrained, optimal


def first_two_witnesses(optimal):
    """Flatten the brute-force enumeration into the first two witnesses."""
    out = []
    for shift, bounds, combos, _pair in optimal:
        for combo in combos:
            if len(out) < 2:
                out.append((shift, bounds, combo))
            else:
                return out, True
    return out, False


class AlignmentTests(unittest.TestCase):
    def test_edit_distance_and_canonical_cigar(self):
        cases = [
            ("ACGT", "ACGT"),
            ("ACGT", "AGGT"),
            ("ACGT", "ACG"),
            ("ACG", "ACGT"),
            ("AAAA", "AA"),
            ("", "ACGT"),
            ("ACGTACGT", "TGCATGCA"),
        ]
        for seg, ref in cases:
            with self.subTest(seg=seg, ref=ref):
                cost, cigar = _needleman(seg, ref)
                self.assertEqual(cost, levenshtein(seg, ref))
                self.assertEqual(replay(seg, ref, cigar)["edits"], cost)

    def test_cigar_replay_rejects_bad_inputs(self):
        self.assertEqual(_needleman("ACGTA", "ACGTA")[1], "5=")
        with self.assertRaises(ValueError):
            replay("ACGTA", "ACGTC", "5=")
        with self.assertRaises(ValueError):
            replay("ACG", "ACGT", "3=")
        with self.assertRaises(ValueError):
            replay("ACGTA", "ACGTA", "3=2X")

    def test_alignment_exercising_all_operation_classes(self):
        # Search guarantees an *optimal* alignment whose CIGAR uses = X I D.
        ref = "ACGTACGTAA"
        found = None
        rng = random.Random(9)
        for _ in range(2000):
            seg = mutate(ref, rng, 3)
            if not (8 <= len(seg) <= 13):
                continue
            cost, cigars, _ = align_all(seg, ref, top_k=3)
            for cigar in cigars:
                ops = {o["op"] for o in replay(seg, ref, cigar)["operations"]}
                if ops == {"=", "X", "I", "D"}:
                    found = (seg, cost, cigar)
                    break
            if found:
                break
        self.assertIsNotNone(found, "random search should hit =,X,I,D")
        seg, cost, cigar = found
        self.assertEqual(cost, levenshtein(seg, ref))

    def test_multiple_optimal_cigars_over_homopolymer(self):
        # 10 A's aligned to 8 A's needs two deletions placed in >=2 spots.
        cost, cigars, count = align_all("AAAAAAAAAA", "AAAAAAAA", top_k=3)
        self.assertEqual(cost, 2)
        self.assertGreaterEqual(count, 2)
        self.assertGreater(len(cigars), 1)
        self.assertEqual(list(cigars), sorted(cigars))
        for cigar in cigars:
            self.assertEqual(
                replay("AAAAAAAAAA", "AAAAAAAA", cigar)["edits"], 2
            )


class SolverBruteForceTests(unittest.TestCase):
    def assert_matches_bruteforce(self, read, ref, copies, per):
        best_un, optimal = brute_force(read, ref, copies, per)
        length = len(ref)
        length_ok = copies * (length - per) <= len(read) <= copies * (length + per)
        try:
            result = decode(ref, read, copies, per)
        except ConstraintFailure as exc:
            self.assertFalse(
                optimal,
                f"solver failed but brute force found {optimal[:1]}",
            )
            payload = exc.payload
            if not length_ok:
                self.assertEqual(payload["reason"], "read_length_outside_window")
            else:
                self.assertEqual(payload["reason"], "edit_budget_exceeded")
                self.assertEqual(
                    payload["best_achieved"]["total_edits"], best_un[0]
                )
                self.assertEqual(
                    payload["best_achieved"]["max_segment_edits"], best_un[1]
                )
            return

        self.assertTrue(optimal, "solver succeeded but brute force found nothing")
        target = optimal[0][3]
        self.assertEqual(result["total_edits"], target[0])
        self.assertEqual(result["max_segment_edits"], target[1])

        expected_wits, expected_more = first_two_witnesses(optimal)
        self.assertGreaterEqual(len(expected_wits), 1)

        if len(expected_wits) == 1 and not expected_more:
            self.assertEqual(result["status"], "unique")
            shift, bounds, combo = expected_wits[0]
            self.assertEqual(result["shift"], shift)
            self.assertEqual([s["end"] for s in result["segments"]],
                             list(bounds))
            self.assertEqual([s["cigar"] for s in result["segments"]],
                             list(combo))
            self._validate_witness(result, read)
        else:
            self.assertEqual(result["status"], "ambiguous")
            self.assertEqual(len(result["witnesses"]), 2)
            for got, (shift, bounds, combo) in zip(
                    result["witnesses"], expected_wits):
                self.assertEqual(got["shift"], shift)
                self.assertEqual(
                    [s["end"] for s in got["segments"]], list(bounds)
                )
                self.assertEqual(
                    [s["cigar"] for s in got["segments"]], list(combo)
                )
                self._validate_witness(got, read)
            keys = [
                (w["shift"],
                 tuple(s["end"] for s in w["segments"]),
                 tuple(s["cigar"] for s in w["segments"]))
                for w in result["witnesses"]
            ]
            self.assertEqual(keys, sorted(keys))
            self.assertLess(keys[0], keys[1])

    def _validate_witness(self, w, read):
        prev = 0
        for seg in w["segments"]:
            piece = read[prev:seg["end"]]
            trace = replay(piece, w["rotated_reference"], seg["cigar"])
            self.assertEqual(seg["edits"], trace["edits"])
            prev = seg["end"]
        self.assertEqual(prev, len(read))
        self.assertEqual(w["total_edits"],
                         sum(s["edits"] for s in w["segments"]))
        self.assertEqual(w["max_segment_edits"],
                         max(s["edits"] for s in w["segments"]))

    def test_designed_read_with_substitution_insertion_deletion(self):
        ref = "ACGTAGCCTA"  # 10 bases
        rotated = ref[3:] + ref[:3]
        c1 = apply_edit(rotated, "S", 2)
        c2 = apply_edit(rotated, "I", 6, base="T")
        c3 = apply_edit(rotated, "D", 4)
        read = c1 + c2 + c3
        self.assertEqual(len(read), 30)
        result = decode(ref, read, 3, 3)
        self.assertIn(result["status"], ("unique", "ambiguous"))
        wits = [result] if result["status"] == "unique" else result["witnesses"]
        w = wits[0]
        kinds = set()
        prev = 0
        for seg in w["segments"]:
            for o in replay(read[prev:seg["end"]], w["rotated_reference"],
                            seg["cigar"])["operations"]:
                if o["op"] != "=":
                    kinds.add(o["op"])
            prev = seg["end"]
        self.assertEqual(kinds, {"X", "I", "D"})
        self.assert_matches_bruteforce(read, ref, 3, 3)

    def test_battery(self):
        rng = random.Random(20261002)
        ref = "ACGTAGCCTA"
        cases = []
        for shift in range(10):
            for per in (0, 1, 3):
                read = (ref[shift:] + ref[:shift]) * 3
                cases.append((read, ref, 3, per))
        for seed in range(6):
            r = random.Random(100 + seed)
            rotated = ref[r.randrange(10):]
            shift = r.randrange(10)
            rotated = ref[shift:] + ref[:shift]
            read = "".join(
                mutate(rotated, r, r.randrange(4)) for _ in range(3)
            )
            if 30 <= len(read) <= 160:
                cases.append((read, ref, 3, 3))
        ref8 = "ACGTAGGC"
        for seed in range(4):
            r = random.Random(300 + seed)
            shift = r.randrange(8)
            rotated = ref8[shift:] + ref8[:shift]
            read = "".join(mutate(rotated, r, r.randrange(3)) for _ in range(4))
            if 30 <= len(read) <= 160:
                cases.append((read, ref8, 4, 2))
        for length in (30, 32, 40):
            cases.append(
                ("".join(rng.choice("ACGT") for _ in range(length)),
                 ref, 3, 3)
            )
        for read, r, copies, per in cases:
            with self.subTest(read=read, ref=r, copies=copies, per=per):
                self.assert_matches_bruteforce(read, r, copies, per)


class FixedScenariosTests(unittest.TestCase):
    def test_homopolymer_cigar_ambiguity_reports_two_witnesses(self):
        ref = "AAAAAAAA"          # one distinct rotation
        read = "A" * 30           # unique boundary optimum: 10/10/10 (6, 2)
        result = decode(ref, read, 3, 3)
        self.assertEqual(result["status"], "ambiguous")
        self.assertEqual(result["total_edits"], 6)
        self.assertEqual(result["max_segment_edits"], 2)
        w1, w2 = result["witnesses"]
        for w in (w1, w2):
            self.assertEqual(w["shift"], 0)
            self.assertEqual([s["end"] for s in w["segments"]], [10, 20, 30])
            self.assertEqual([s["edits"] for s in w["segments"]], [2, 2, 2])
        c1 = [s["cigar"] for s in w1["segments"]]
        c2 = [s["cigar"] for s in w2["segments"]]
        self.assertNotEqual(c1, c2)
        # global stable ordering: same shift/boundaries, so CIGARs order it
        self.assertLess(tuple(c1), tuple(c2))

    def test_zero_edit_cap_exact_repeat_is_unique(self):
        ref = "ACGTGGTCAA"  # 10 bases, no inner period
        shift = 4
        rotated = ref[shift:] + ref[:shift]
        read = rotated * 3
        self.assertEqual(len(read), 30)
        result = decode(ref, read, 3, 0)
        self.assertEqual(result["status"], "unique")
        self.assertEqual(result["shift"], shift)
        self.assertEqual(result["total_edits"], 0)
        self.assertEqual(result["boundaries"], [10, 20])
        self.assertEqual([s["length"] for s in result["segments"]],
                         [10, 10, 10])
        self.assertTrue(all(s["cigar"] == "10=" for s in result["segments"]))

    def test_length_window_constraint_failure_is_locatable(self):
        ref = "ACGTACGTA"  # length 9
        read = "A" * 30    # copies=3, cap 0 forces total length exactly 27
        with self.assertRaises(ConstraintFailure) as ctx:
            decode(ref, read, 3, 0)
        p = ctx.exception.payload
        self.assertEqual(p["error"], "constraint_failure")
        self.assertEqual(p["field"], "read")
        self.assertEqual(p["reason"], "read_length_outside_window")
        self.assertEqual(p["constraints"]["total_length_window"], [27, 27])
        self.assertIn("30", p["message"])

    def test_edit_budget_constraint_failure_reports_best_achieved(self):
        ref = "ACGTACGTAC"  # length 10
        read = "G" * 30
        with self.assertRaises(ConstraintFailure) as ctx:
            decode(ref, read, 3, 1)
        p = ctx.exception.payload
        self.assertEqual(p["reason"], "edit_budget_exceeded")
        self.assertEqual(p["field"], "max_edits_per_copy")
        self.assertIsNotNone(p["best_achieved"])
        self.assertGreater(p["best_achieved"]["total_edits"], 3)
        self.assertEqual(len(p["per_shift"]),
                         len({(ref[s:] + ref[:s]) for s in range(10)}))

    def test_validation_errors_are_locatable(self):
        with self.assertRaises(InvalidRequest) as ctx:
            decode("ACGT", "A" * 30, 3, 1)
        self.assertEqual(ctx.exception.field, "reference")
        with self.assertRaises(InvalidRequest) as ctx:
            decode("ACGTACGTAA", "A" * 20, 3, 1)
        self.assertEqual(ctx.exception.field, "read")
        with self.assertRaises(InvalidRequest) as ctx:
            decode("ACGTACGTAA", "A" * 30, 9, 1)
        self.assertEqual(ctx.exception.field, "copies")
        with self.assertRaises(InvalidRequest) as ctx:
            decode("ACGTACGTAA", "A" * 30, 3, 5)
        self.assertEqual(ctx.exception.field, "max_edits_per_copy")
        with self.assertRaises(InvalidRequest) as ctx:
            decode("ACGTACNTAA", "A" * 30, 3, 1)
        self.assertEqual(ctx.exception.field, "reference")

    def test_inner_period_de_duplicates_rotations(self):
        ref = "ACGTACGT"  # period 4
        rotated = ref[1:] + ref[:1]
        read = rotated * 4  # length 32
        result = decode(ref, read, 4, 0)
        self.assertEqual(result["status"], "unique")
        self.assertEqual(result["shift"], 1)

    def test_four_copies_objective_total_then_worst(self):
        ref = "ACGTACGT"
        rng = random.Random(55)
        rotated = ref[2:] + ref[:2]

        def sub_one(seq):
            out = list(seq)
            i = rng.randrange(len(out))
            out[i] = rng.choice([b for b in "ACGT" if b != out[i]])
            return "".join(out)

        read = "".join(sub_one(rotated) for _ in range(4))
        self.assertEqual(len(read), 32)
        best, optimal = brute_force(read, ref, 4, 2)
        self.assertTrue(optimal)
        result = decode(ref, read, 4, 2)
        # derive expected pair from the brute force witness CIGAR edits
        shift, bounds, combos, pair = optimal[0]
        combo = combos[0]
        edits = [
            replay(read[(0 if k == 0 else bounds[k - 1]):bounds[k]],
                   ref[shift:] + ref[:shift], c)["edits"]
            for k, c in enumerate(combo)
        ]
        self.assertEqual(result["total_edits"], pair[0])
        self.assertEqual(result["max_segment_edits"], pair[1])
        self.assertEqual(best, pair)


if __name__ == "__main__":
    unittest.main(verbosity=2)
