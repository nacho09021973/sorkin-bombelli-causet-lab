# Foundation Benchmarks

This directory holds the frozen *protocol* of the foundation
validation benchmark for causal-set diagnostics and historical
annealer recoverability probes.

The benchmark is a grid of uniform Poisson sprinklings of the
unit Minkowski causal diamond. Each cell is identified by
`(d_spacetime, n, seed)` and produced deterministically by
`validation_suite.sprinkle_minkowski_diamond(n, seed, d_spacetime)`.
Because the sprinkler is reproducible, we freeze only the
precomputed order-theoretic invariants and let the actual
causets be regenerated on demand.

## Cells

### d_spacetime = 2

- sizes: 16, 32, 64
- seeds: 1959, 1962, 1987, 2009, 2026
- expected Myrheim-Meyer dimension: 2 (closed form for unit diamond)

### d_spacetime = 3

- sizes: 16, 32, 64
- seeds: 1959, 1962, 1987, 2009, 2026
- expected Myrheim-Meyer dimension: 3 (closed form for unit diamond)

### d_spacetime = 4

- sizes: 16, 32, 64
- seeds: 1959, 1962, 1987, 2009, 2026
- expected Myrheim-Meyer dimension: 4 (closed form for unit diamond)

## Frozen invariants

`invariants.json` records the full fingerprint of each cell.
A single entry looks like this:

```json
{
  "d_spacetime": 2,
  "n": 16,
  "seed": 1959,
  "fingerprint": {
    "n": 16,
    "relation_count": 70,
    "link_count": 21,
    "ordering_fraction": 0.5833333333333334,
    "myrheim_meyer_dim": 1.7919472299399786,
    "height": 7,
    "antichain_profile": [
      3,
      1,
      2,
      3,
      4,
      2,
      1
    ],
    "chain_counts": {
      "2": 70,
      "3": 125,
      "4": 106
    }
  }
}
```

The fields are described in `causet_invariants.py`.

## Regenerating

Run `python3 tools/build_foundation_benchmarks.py` from the
repository root. The script overwrites `invariants.json` and
this README. The accompanying regression test
`tests/test_foundation_benchmarks.py` then verifies that the
current code reproduces the frozen invariants bit for bit.

## Why not freeze the matrices?

Two reasons. First, the sprinkler is deterministic under a
fixed seed, so the matrix is fully determined by the cell
key; storing it again would be redundant and would invite
drift. Second, freezing the *invariants* directly tests the
scientifically meaningful claim: that the order-theoretic
diagnostics of these reference causets are stable. A change
to the sprinkler or the invariant computation is detected
immediately, while incidental serialization changes do not
perturb the benchmark.

## SORKIN-2 interpretation guardrails

The v1.0.0 Bombelli revival remains complete as a historical
reconstruction. SORKIN-2 uses the foundation artifacts as diagnostics
for algorithmic recoverability in the Bombelli annealer.

The distinction to preserve is:

- causal realization / embedding existence;
- annealer accessibility / algorithmic recoverability.

Annealer failure must not be interpreted as non-embeddability, and low
final energy must not be interpreted as manifoldlikeness. The current
diagnostic question is whether the historical energy, move set,
schedule, and acceptance rule access known-truth zero-energy
realizations.

Existing Phase 1 and Phase 2 provenance should be preserved. If new
wording conflicts with old provenance, keep the old provenance and add a
TODO rather than deleting historical context.

Diagnostic artifact status:

- Phase2E/2F: warmup diagnostics.
- Phase4C/4D: optimizer-seed and robustness diagnostics.
- Phase3/4A/4B/5: exploratory material unless later formalized as
  provenance-frozen known-truth diagnostics.

KAN, PySR, and GAM are deferred until a clean multi-family known-truth
dataset exists.
