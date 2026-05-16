# Foundation Benchmarks

This directory holds the frozen *protocol* of the foundation
validation benchmark for causet embedding methods.

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

## Scaling atlases

`phase1_atlas.{csv,md}` records the fixed-size order-theoretic
atlas comparing Minkowski sprinklings with Kleitman-Rothschild
controls.

`phase1b_scaling_atlas.{csv,md}` extends that comparison across
`n in {32, 64, 128, 256}`.

`phase1c_scaling_atlas.{csv,md}` keeps the Phase 1B grid and
adds suspended corona/crown posets as a second non-manifoldlike
control.

`phase1d_structural_atlas.{csv,md}` keeps the Phase 1C grid and
adds raw 2-chain/3-chain counts plus normalized 3-chain
abundance as a third order-theoretic diagnostic.

`phase2_embedding_bridge.{csv,md}` runs a small fixed probe of
the historical Bombelli-Sorkin annealer on selected Phase 1D
families and records both pre-embedding diagnostics and
post-embedding metrics.

`phase2b_annealer_schedule_probe.{csv,md}` varies only the
historical annealer's iteration budget (``warmup_limit``,
``anneal_limit``, ``max_data``) over a Minkowski-only grid in
``d in {2, 3, 4}``, ``n in {32, 64}`` and three seeds, and
records ``energy_gap``, ``interval_rmse`` and a conservative
``success_flag`` per run. Non-manifoldlike controls are
deliberately excluded; the probe is *not* a manifoldness
classifier.

`phase2c_oracle_embedding_audit.{csv,md}` is a ground-truth
consistency check: for each Minkowski case in the Phase 2/2B
grid the oracle evaluates the Bombelli energy, reconstructs
the causal matrix, and computes the Lorentz-invariant interval
residual — all directly at the ground-truth coordinates,
without any annealing. The verdict determines whether the
failure in Phase 2/2B is a convention inconsistency (oracle
fails) or an optimizer problem (oracle passes).

`phase2d_initialization_basin_audit.{csv,md}` audits the basin
structure around the truth by injecting four initialization
strategies (truth, truth+small noise, truth+medium noise,
random_init) into ConesSimulator and measuring energy, interval
RMSE, and coordinate distance before and after warmup+anneal.
Determines whether the failure is in the warmup dynamics, the
move set, or the basin structure.

## Regenerating

Run `python3 tools/build_foundation_benchmarks.py` from the
repository root. The script overwrites `invariants.json` and
this README. The accompanying regression test
`tests/test_foundation_benchmarks.py` then verifies that the
current code reproduces the frozen invariants bit for bit.

Run `make regen-phase1`, `make regen-phase1b`,
`make regen-phase1c`, or `make regen-phase1d` to regenerate the
corresponding atlas and its focused regression test.

Run `make regen-phase2` to regenerate the embedding bridge and
its focused regression test, `make regen-phase2b` to regenerate
the annealer schedule probe, `make regen-phase2c` to regenerate
the oracle embedding audit, or `make regen-phase2d` to regenerate
the initialization/basin audit.

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
