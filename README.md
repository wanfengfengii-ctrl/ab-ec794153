# Concatemer Decoder

Joint decoding of rolling-circle amplification (RCA) concatemer reads.
Given a circular reference, a noisy read, an exact copy count and an
edit cap per copy, the service jointly selects:

1. one cyclic shift (rotation) of the reference,
2. exactly `copies` consecutive, non-empty segmentation boundaries,
3. one global (Needleman–Wunsch) alignment per segment,

optimising lexicographically by **(total edits, maximum edits of any
single segment)**.  Witnesses are stably ordered by
**(shift, boundaries, per-segment CIGAR)**.  No third-party packages are
required — only the Python 3.11 standard library.

## API

`POST /api/concatemers/decode`

```json
{
  "reference": "ACGTGGTCAA",
  "read": "GGTCAAACGTGGTCAAACGTGGTCAA",
  "copies": 3,
  "max_edits_per_copy": 0
}
```

| field | range |
| --- | --- |
| `reference` | 8–20 bases, `ACGT` |
| `read` | 30–160 bases, `ACGT` |
| `copies` | 3–8 |
| `max_edits_per_copy` | 0–3 |

Responses:

- **200 `unique`** — one optimal interpretation: `shift`,
  `rotated_reference`, internal `boundaries`, `total_edits`,
  `max_segment_edits` and each segment (`start`, `end`, `length`,
  `edits`, replayable `cigar`).
- **200 `ambiguous`** — several optimal interpretations; returns the
  objective values and the first two witnesses in stable order.
  Distinct optimal global alignments (e.g. where a deletion falls inside
  a homopolymer) count as distinct witnesses.
- **422 `constraint_failure`** — no feasible interpretation.  The body is
  locatable (`field`), explains the `reason`
  (`read_length_outside_window` / `edit_budget_exceeded` /
  `no_segmentation`), echoes the active constraints, and — for edit
  failures — reports the best unconstrained result and the value at every
  distinct rotation.
- **422 `invalid_request`** — malformed/out-of-range fields, identified by
  `field`.

CIGARs use `=` (match), `X` (substitution), `I` (read insertion) and
`D` (read deletion), run-length encoded (e.g. `8=1X1I`).  Every emitted
CIGAR is replayed against both strings before it leaves the solver, so it
always reconstructs the segment and the rotated reference exactly.

`GET /healthz` returns `{"status":"ok"}`.

## Algorithm notes

- Segments of length `L ± e` are aligned to the rotated reference with
  unit-cost global alignment; a partition DP over boundary positions
  optimises (total edits, worst segment).
- The forward partition DP also carries a saturated count of optimal
  witnesses (per-edge optimal-alignment counts multiply, incoming edges
  add), giving an exact unique/ambiguous decision without enumeration.
- Per-cell alignment counts are exact: diagonal/down/right predecessors
  append different final operations, so their traceback sets are
  disjoint.  Candidate CIGAR strings are grouped by diagonal-move count
  before lexicographic k-best pruning, which keeps the pruning provably
  exact even when substitutions tie with insertion/deletion pairs.
- Rotations producing the same sequence (references with an inner
  period) are de-duplicated to the smallest representing shift.

## Run with Docker Compose

```bash
# build and start the API on the default host port 8080
docker compose up --build api

# custom host port
PORT=9090 docker compose up --build api
```

The API container includes a `/healthz` healthcheck.

### One-shot verification service

```bash
docker compose run --rm verify
```

`verify` waits for the API to become healthy, then runs:

1. a build check (byte-compilation of all sources),
2. the unit-test suite (including brute-force cross-validation of the
   dynamic program on many small instances),
3. an end-to-end decode smoke test covering substitution, insertion and
   deletion plus the unique, ambiguous and constraint-failure paths.

Its process exit code is 0 only when every stage passes.

## Local development

```bash
python -m unittest discover -s tests -v          # tests
python -m app.main                                # serve on $PORT (8080)
API_URL=http://127.0.0.1:8080 python scripts/smoke.py
```
