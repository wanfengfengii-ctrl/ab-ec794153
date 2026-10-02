"""Joint decoding of rolling-circle concatemer reads.

For a circular reference we jointly choose:

1. a cyclic shift (rotation) of the reference,
2. exactly ``copies`` consecutive, non-empty segmentation boundaries,
3. a global (Needleman-Wunsch) alignment of every segment against the
   rotated reference.

The objective is lexicographic: first minimise the *total* edit count,
then the maximum edit count of any single segment.  When several optimal
alignments exist (common over homopolymer runs), every distinct optimal
CIGAR is itself a distinct witness.  Witnesses are ordered stably by
``(shift, boundaries, CIGARs)``.

Rotations producing the identical rotated sequence (references with an
inner period) are de-duplicated: the explanations differ only by an
unobservable phase coordinate.

Only the Python standard library is used.
"""

from __future__ import annotations

import heapq
import re
from typing import Dict, List, Optional, Sequence, Tuple

ALPHABET = frozenset("ACGT")

Pair = Tuple[int, int]  # (total edits, worst single-segment edits)

CIGAR_RE = re.compile(r"(\d+)([=XID])")

# CIGAR strings kept per (cell, diagonal-count) group.  Two witnesses
# already prove ambiguity; a margin of one is carried for safety.
TOP_K = 3
# Counts saturate here; the solver only needs to distinguish 1 from "many".
COUNT_CAP = 2


class InvalidRequest(ValueError):
    """Request fields fail validation. Carries a locatable field path."""

    def __init__(self, message: str, field: str):
        super().__init__(message)
        self.field = field
        self.message = message


class ConstraintFailure(RuntimeError):
    """No interpretation exists inside the edit / length budget."""

    def __init__(self, payload: dict):
        super().__init__(payload["message"])
        self.payload = payload


# ---------------------------------------------------------------------------
# CIGAR helpers
# ---------------------------------------------------------------------------

def _encode_cigar(ops: Sequence[str]) -> str:
    out: List[str] = []
    i = 0
    while i < len(ops):
        j = i + 1
        while j < len(ops) and ops[j] == ops[i]:
            j += 1
        out.append(f"{j - i}{ops[i]}")
        i = j
    return "".join(out)


def replay(seg: str, ref: str, cigar: str) -> dict:
    """Replay a CIGAR against the segment and reference.

    Verifies that the CIGAR reconstructs exactly both strings and returns
    the per-operation trace plus the edit count.  Raises ``ValueError`` on
    any inconsistency, which makes every emitted CIGAR directly auditable.
    """
    operations: List[dict] = []
    si = ri = 0
    pos = 0
    for match in CIGAR_RE.finditer(cigar):
        if match.start() != pos:
            raise ValueError(f"malformed CIGAR near offset {pos}: {cigar!r}")
        pos = match.end()
        length = int(match.group(1))
        op = match.group(2)
        if length == 0:
            raise ValueError("zero-length CIGAR operation")
        if op == "=":
            if (len(seg) - si < length or len(ref) - ri < length or
                    seg[si:si + length] != ref[ri:ri + length]):
                raise ValueError(f"={length} run does not match both strings")
            si += length
            ri += length
            edits = 0
        elif op == "X":
            if (len(seg) - si < length or len(ref) - ri < length or
                    any(seg[si + k] == ref[ri + k] for k in range(length))):
                raise ValueError("X run is invalid")
            si += length
            ri += length
            edits = length
        elif op == "I":
            if len(seg) - si < length:
                raise ValueError("I run exceeds the segment")
            si += length
            edits = length
        else:  # D
            if len(ref) - ri < length:
                raise ValueError("D run exceeds the reference")
            ri += length
            edits = length
        operations.append({"op": op, "length": length, "edits": edits})
    if pos != len(cigar):
        raise ValueError(f"trailing garbage in CIGAR: {cigar!r}")
    if si != len(seg) or ri != len(ref):
        raise ValueError("CIGAR does not reconstruct both strings fully")
    return {
        "edits": sum(o["edits"] for o in operations),
        "operations": operations,
    }


# ---------------------------------------------------------------------------
# Needleman-Wunsch with grouped retention of the best CIGAR strings
# ---------------------------------------------------------------------------

def _cell_fill(piece: str, ref: str, top_k: int):
    """Fill the NW grid.

    Returns
    -------
    f : edit-cost table
    cf : table with the exact number of minimum-cost tracebacks per cell,
         saturated at :data:`COUNT_CAP`
    best : ``best[i][j]`` maps the diagonal-move count ``q`` to up to
         ``top_k`` best *raw* op strings (one character per operation).

    Why group by ``q``: a traceback reaching ``(i, j)`` with ``q`` diagonal
    moves contains exactly ``q`` match/substitution ops, ``i - q``
    insertions and ``j - q`` deletions; within one q-group all op strings
    are equal-length anagrams.  Pruning such a group to its ``top_k``
    lexicographically smallest strings is exact — equal-length strings
    can never prefix one another, and appending a common continuation
    preserves serialized-CIGAR order (this monotonicity was verified
    exhaustively).  Different q groups trade substitutions against
    insertion/deletion pairs at equal edit cost and must not prune each
    other.

    Why the counts are exact: predecessor kinds (diagonal / down / right)
    append different final operations, so their traceback sets are
    disjoint and the counts simply add.  Move sequences and op strings are
    in bijection, and run-length encoding is injective.
    """
    m, n = len(piece), len(ref)
    f = [[0] * (n + 1) for _ in range(m + 1)]
    cf = [[0] * (n + 1) for _ in range(m + 1)]
    best: List[List[Dict[int, List[str]]]] = [
        [{} for _ in range(n + 1)] for _ in range(m + 1)
    ]
    best[0][0] = {0: [""]}
    cf[0][0] = 1
    for i in range(1, m + 1):
        f[i][0] = i
        cf[i][0] = 1
        best[i][0] = {0: ["I" * i]}
    for j in range(1, n + 1):
        f[0][j] = j
        cf[0][j] = 1
        best[0][j] = {0: ["D" * j]}

    for i in range(1, m + 1):
        rb = piece[i - 1]
        for j in range(1, n + 1):
            diag_cost = f[i - 1][j - 1] + (rb != ref[j - 1])
            down_cost = f[i - 1][j] + 1   # read insertion -> I
            right_cost = f[i][j - 1] + 1  # reference deletion -> D
            cost = min(diag_cost, down_cost, right_cost)
            f[i][j] = cost
            cell: Dict[int, List[str]] = {}
            count = 0

            def absorb(q: int, additions: List[str]) -> None:
                merged = list(cell.get(q, ()))
                merged.extend(additions)
                merged.sort(key=_encode_cigar)
                uniq: List[str] = []
                for c in merged:
                    if not uniq or c != uniq[-1]:
                        uniq.append(c)
                cell[q] = uniq[:top_k]

            if diag_cost == cost:
                op = "=" if rb == ref[j - 1] else "X"
                for q, strings in best[i - 1][j - 1].items():
                    absorb(q + 1, [s + op for s in strings])
                count += cf[i - 1][j - 1]
            if down_cost == cost:
                for q, strings in best[i - 1][j].items():
                    absorb(q, [s + "I" for s in strings])
                count += cf[i - 1][j]
            if right_cost == cost:
                for q, strings in best[i][j - 1].items():
                    absorb(q, [s + "D" for s in strings])
                count += cf[i][j - 1]
            best[i][j] = cell
            cf[i][j] = min(COUNT_CAP, count)
    return f, cf, best


def _final_cigars(best_cell: Dict[int, List[str]], top_k: int) -> Tuple[str, ...]:
    """Merge diagonal-count groups at the sink into the best RLE CIGARs."""
    merged: List[str] = []
    for strings in best_cell.values():
        merged.extend(_encode_cigar(s) for s in strings)
    merged.sort()
    picked: List[str] = []
    for c in merged:
        if not picked or c != picked[-1]:
            picked.append(c)
        if len(picked) == top_k:
            break
    return tuple(picked)


def align_all(seg: str, ref: str,
              top_k: int = TOP_K) -> Tuple[int, Tuple[str, ...], int]:
    """Return ``(distance, best CIGARs, saturated optimal-CIGAR count)``."""
    f, cf, best = _cell_fill(seg, ref, top_k)
    m, n = len(seg), len(ref)
    return f[m][n], _final_cigars(best[m][n], top_k), cf[m][n]


def _int_fill(piece: str, ref: str):
    """Scalar NW fill: edit distances and saturated traceback counts.

    This is the cheap phase used for every candidate segment; it keeps no
    strings and only the sink column of each row is read later (segments
    share the start position, so one fill serves every allowed length).
    """
    m, n = len(piece), len(ref)
    prev = list(range(n + 1))
    prev_cf = [1] * (n + 1)
    dist_col = [0] * (m + 1)
    count_col = [1] * (m + 1)
    for i in range(1, m + 1):
        cur = [i] + [0] * n
        cur_cf = [1] + [0] * n
        rb = piece[i - 1]
        for j in range(1, n + 1):
            diag_cost = prev[j - 1] + (rb != ref[j - 1])
            down_cost = prev[j] + 1
            right_cost = cur[j - 1] + 1
            cost = min(diag_cost, down_cost, right_cost)
            cur[j] = cost
            count = 0
            if diag_cost == cost:
                count += prev_cf[j - 1]
            if down_cost == cost:
                count += prev_cf[j]
            if right_cost == cost:
                count += cur_cf[j - 1]
            cur_cf[j] = min(COUNT_CAP, count)
        dist_col[i] = cur[n]
        count_col[i] = cur_cf[n]
        prev, prev_cf = cur, cur_cf
    return dist_col, count_col


def _needleman(seg: str, ref: str) -> Tuple[int, str]:
    """Distance and the single canonical (lexicographically first) CIGAR."""
    cost, cigars, _ = align_all(seg, ref, top_k=1)
    return cost, cigars[0]


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def _validate_sequence(value: object, field: str, lo: int, hi: int) -> str:
    if not isinstance(value, str):
        raise InvalidRequest(f"{field} must be a string", field)
    seq = value.strip().upper()
    if not (lo <= len(seq) <= hi):
        raise InvalidRequest(
            f"{field} length must be between {lo} and {hi} bases", field
        )
    bad = sorted({ch for ch in seq if ch not in ALPHABET})
    if bad:
        raise InvalidRequest(
            f"{field} contains unsupported bases: {''.join(bad)}; allowed: ACGT",
            field,
        )
    return seq


def _validate_int(value: object, field: str, lo: int, hi: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise InvalidRequest(f"{field} must be an integer", field)
    if not (lo <= value <= hi):
        raise InvalidRequest(f"{field} must be between {lo} and {hi}", field)
    return value


# ---------------------------------------------------------------------------
# Partition DP over segmentation boundaries (fixed rotation)
# ---------------------------------------------------------------------------

def _add_cost(p: Pair, cost: int) -> Pair:
    return (p[0] + cost, max(p[1], cost))


def _merge(p: Pair, q: Pair) -> Pair:
    return (p[0] + q[0], max(p[1], q[1]))


def _fill_costs(read: str, ref: str, lo: int, hi: int):
    """costs[a][b] = ``(edit distance, saturated optimal-alignment count)``
    for read[a:b] against ``ref``.

    One scalar NW fill per start position serves every allowed length via
    the sink column.  No strings are produced in this phase.
    """
    n = len(read)
    costs: List[Optional[Dict[int, Tuple[int, int]]]] = [None] * (n + 1)
    for a in range(0, n - lo + 1):
        max_len = min(hi, n - a)
        if max_len < lo:
            continue
        dist_col, count_col = _int_fill(read[a:a + max_len], ref)
        costs[a] = {
            a + length: (dist_col[length], count_col[length])
            for length in range(lo, max_len + 1)
        }
    return costs


def _edge_cigars(piece: str, ref: str) -> Tuple[str, ...]:
    """Best CIGAR strings for one segment (only used for emitted witnesses)."""
    _cost, cigars, _count = align_all(piece, ref, top_k=TOP_K)
    return cigars


def _forward(costs, n: int, k_total: int, lo: int, hi: int):
    """Forward partition DP.

    Layer ``k`` maps end position to ``(best pair, saturated witness
    count)``; witness counts multiply along edges (each edge carries
    ``align_count`` optimal global alignments) and add over incoming
    edges, saturated at :data:`COUNT_CAP`.
    """
    dp: List[Dict[int, Tuple[Pair, int]]] = [{0: ((0, 0), 1)}]
    for k in range(1, k_total + 1):
        prev, cur = dp[k - 1], {}
        remaining = k_total - k
        for a, (pair, ways) in prev.items():
            ends = costs[a]
            if not ends:
                continue
            for b, (cost, align_ways) in ends.items():
                left = n - b
                if left < remaining * lo or left > remaining * hi:
                    continue  # the remaining segments could never fit
                cand = _add_cost(pair, cost)
                extra = min(COUNT_CAP, ways * align_ways)
                old = cur.get(b)
                if old is None or cand < old[0]:
                    cur[b] = (cand, extra)
                elif cand == old[0]:
                    cur[b] = (cand, min(COUNT_CAP, old[1] + extra))
        dp.append(cur)
    return dp


def _reverse(costs, n: int, k_total: int, lo: int, hi: int):
    rdp: List[Dict[int, Pair]] = [{n: (0, 0)}]
    for k in range(1, k_total + 1):
        prev, cur = rdp[k - 1], {}
        for a in range(0, n - k * lo + 1):
            ends = costs[a]
            if not ends:
                continue
            best_pair: Optional[Pair] = None
            for b, (cost, _ways) in ends.items():
                tail = prev.get(b)
                if tail is None:
                    continue
                cand = _merge((cost, cost), tail)
                if best_pair is None or cand < best_pair:
                    best_pair = cand
            if best_pair is not None:
                cur[a] = best_pair
        rdp.append(cur)
    return rdp


def _enumerate_boundaries(costs, rdp, n, k_total, lo, hi, target: Pair,
                          limit: int):
    """Yield up to ``limit`` optimal boundary tuples in lexicographic order."""
    results: List[Tuple[int, ...]] = []

    def rec(k: int, start: int, bounds: List[int], pair: Pair) -> None:
        if len(results) >= limit:
            return
        if k == k_total:
            if start == n and pair == target:
                results.append(tuple(bounds))
            return
        remaining_after = k_total - k - 1
        for b in sorted(costs[start] or {}):
            left = n - b
            if left < remaining_after * lo or left > remaining_after * hi:
                continue
            cand = _add_cost(pair, costs[start][b][0])
            tail = rdp[remaining_after].get(b)
            if tail is None or _merge(cand, tail) != target:
                continue
            bounds.append(b)
            rec(k + 1, b, bounds, cand)
            bounds.pop()

    rec(0, 0, [], (0, 0))
    return results


def _combo_iter(per_segment_lists: Sequence[Tuple[str, ...]]):
    """Yield tuples of the cartesian product in lexicographic order.

    Per-segment lists are already sorted; a heap over index vectors
    produces the product lazily without materialising it.
    """
    lists = [list(seq) for seq in per_segment_lists]
    if any(not seq for seq in lists):
        return
    width = len(lists)
    start = (0,) * width
    heap: List[Tuple[Tuple[str, ...], Tuple[int, ...]]] = [
        (tuple(lists[i][0] for i in range(width)), start)
    ]
    seen = {start}
    while heap:
        key, indices = heapq.heappop(heap)
        yield key
        for i in range(width):
            if indices[i] + 1 < len(lists[i]):
                nxt = list(indices)
                nxt[i] += 1
                nxt_t = tuple(nxt)
                if nxt_t not in seen:
                    seen.add(nxt_t)
                    heapq.heappush(
                        heap,
                        (tuple(lists[d][nxt_t[d]] for d in range(width)),
                         nxt_t),
                    )


def _rotations(ref: str) -> List[Tuple[int, str]]:
    """Distinct rotations as (smallest shift, rotated string), by shift."""
    first: Dict[str, int] = {}
    for shift in range(len(ref)):
        rotated = ref[shift:] + ref[:shift]
        if rotated not in first:
            first[rotated] = shift
    return sorted(((shift, rotated) for rotated, shift in first.items()),
                  key=lambda item: item[0])


def _build_witness(shift: int, rotated: str, read: str,
                   bounds: Tuple[int, ...], cigar_tuple: Tuple[str, ...],
                   costs) -> dict:
    segments = []
    total = 0
    worst = 0
    prev = 0
    for end, cigar in zip(bounds, cigar_tuple):
        piece = read[prev:end]
        edits = replay(piece, rotated, cigar)["edits"]  # validates the CIGAR
        if costs[prev][end][0] != edits:
            raise AssertionError("CIGAR cost disagrees with partition DP")
        segments.append({
            "start": prev,
            "end": end,
            "length": end - prev,
            "edits": edits,
            "cigar": cigar,
        })
        total += edits
        worst = max(worst, edits)
        prev = end
    return {
        "shift": shift,
        "rotated_reference": rotated,
        "boundaries": list(bounds[:-1]),
        "total_edits": total,
        "max_segment_edits": worst,
        "segments": segments,
    }


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def decode(reference: object, read: object, copies: object,
           max_edits_per_copy: object) -> dict:
    """Decode one concatemer request.

    Returns a JSON-serialisable response (status ``unique`` or
    ``ambiguous``) or raises :class:`InvalidRequest` /
    :class:`ConstraintFailure`.
    """
    ref = _validate_sequence(reference, "reference", 8, 20)
    seq = _validate_sequence(read, "read", 30, 160)
    k_total = _validate_int(copies, "copies", 3, 8)
    per_copy = _validate_int(max_edits_per_copy, "max_edits_per_copy", 0, 3)

    length = len(ref)
    n = len(seq)
    seg_lo, seg_hi = length - per_copy, length + per_copy
    total_lo, total_hi = k_total * seg_lo, k_total * seg_hi
    budget = k_total * per_copy

    constraints = {
        "copies": k_total,
        "max_edits_per_copy": per_copy,
        "total_edit_budget": budget,
        "segment_length_window": [seg_lo, seg_hi],
        "total_length_window": [total_lo, total_hi],
        "read_length": n,
    }

    if not (total_lo <= n <= total_hi):
        raise ConstraintFailure({
            "error": "constraint_failure",
            "message": (
                f"Read length {n} cannot be split into {k_total} non-empty "
                f"segments of length {seg_lo}..{seg_hi} "
                f"(required total length {total_lo}..{total_hi})."
            ),
            "reason": "read_length_outside_window",
            "field": "read",
            "constraints": constraints,
            "best_achieved": None,
        })

    shifts = _rotations(ref)
    solved = []  # (shift, rotated, best_pair, ways, costs, rdp)
    global_best: Optional[Pair] = None
    for shift, rotated in shifts:
        costs = _fill_costs(seq, rotated, seg_lo, seg_hi)
        dp = _forward(costs, n, k_total, seg_lo, seg_hi)
        entry = dp[k_total].get(n)
        if entry is None:
            continue
        best, ways = entry
        rdp = _reverse(costs, n, k_total, seg_lo, seg_hi)
        solved.append((shift, rotated, best, ways, costs, rdp))
        if global_best is None or best < global_best:
            global_best = best

    if global_best is None:
        raise ConstraintFailure({
            "error": "constraint_failure",
            "message": "No segmentation exists within the segment length window.",
            "reason": "no_segmentation",
            "field": "read",
            "constraints": constraints,
            "best_achieved": None,
        })

    if global_best[0] > budget or global_best[1] > per_copy:
        raise ConstraintFailure({
            "error": "constraint_failure",
            "message": (
                f"Best decoding needs {global_best[0]} edits "
                f"(budget {budget}) with a worst segment of "
                f"{global_best[1]} edits (cap {per_copy})."
            ),
            "reason": "edit_budget_exceeded",
            "field": "max_edits_per_copy",
            "constraints": constraints,
            "best_achieved": {
                "shift": next(s[0] for s in solved if s[2] == global_best),
                "total_edits": global_best[0],
                "max_segment_edits": global_best[1],
            },
            "per_shift": [
                {
                    "shift": shift,
                    "total_edits": pair[0],
                    "max_segment_edits": pair[1],
                }
                for shift, _rot, pair, _w, _c, _r in solved
            ],
        })

    # Global witness multiplicity follows from the forward DP's saturated
    # counts.  Witnesses are enumerated in stable order: rotations by
    # shift, boundary tuples lexicographically, CIGAR tuples per boundary
    # from a sorted cartesian product.  The second witness can only be the
    # second CIGAR tuple of the first boundary tuple or a tuple of the
    # second boundary tuple, so two of each suffice.
    total_ways = min(
        COUNT_CAP,
        sum(s[3] for s in solved if s[2] == global_best),
    )

    witnesses: List[dict] = []
    more_exist = total_ways >= COUNT_CAP
    for shift, rotated, pair, _ways, costs, rdp in solved:
        if pair != global_best:
            continue
        for bounds in _enumerate_boundaries(
            costs, rdp, n, k_total, seg_lo, seg_hi, global_best, limit=2
        ):
            prev = 0
            cigar_lists = []
            for end in bounds:
                cigar_lists.append(_edge_cigars(seq[prev:end], rotated))
                prev = end
            for cigar_tuple in _combo_iter(cigar_lists):
                if len(witnesses) < 2:
                    witnesses.append(
                        _build_witness(shift, rotated, seq, bounds,
                                       cigar_tuple, costs)
                    )
                else:
                    break
            if len(witnesses) >= 2:
                break
        if len(witnesses) >= 2:
            break

    if len(witnesses) == 1 and not more_exist:
        result = dict(witnesses[0])
        result["status"] = "unique"
        return result

    return {
        "status": "ambiguous",
        "total_edits": global_best[0],
        "max_segment_edits": global_best[1],
        "witnesses": witnesses[:2],
        "message": (
            "Multiple optimal decodings exist; the first two witnesses "
            "(ordered by shift, boundaries, CIGAR) are returned."
        ),
    }
