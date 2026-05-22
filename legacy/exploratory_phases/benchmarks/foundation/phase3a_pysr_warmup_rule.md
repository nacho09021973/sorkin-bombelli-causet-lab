# Phase 3A — PySR Symbolic Regression on Causal Set Pipeline Features

**Frontier note**: this is the first application of PySR / symbolic
regression to data derived from Causal Set Theory.  No prior published
work has trained PySR on causal matrix observables.

Two experiments search for closed-form rules over features already
produced by the Bombelli-Sorkin annealing pipeline.

## Experiment A — Order-theoretic predictor

**Question**: can a formula built from *pure* combinatorial invariants
of the causal matrix (computable without any embedding) predict how
much the energy drifts during warmup+anneal from a near-truth start?

If yes, that formula encodes a relationship between discrete causal
order and the current optimizer-response target. This is an exploratory optimizer-response result, not a standalone physical claim.

- Samples: 180
- Features (X): mm_dim, midpoint_dim, abs_discrepancy_mm_midpoint, chain2_count, chain3_count, chain3_abundance, n, target_dim, noise_epsilon, initial_energy
- Target (y): `log1p(max(0, delta_energy))`
- Rows: `truth_plus_small_noise` + `truth_plus_medium_noise`, all warmup modes
- PySR iterations: 100

### Discovered equations (Pareto front)

| complexity | loss | equation | best |
| ---: | ---: | --- | :---: |
| 1 | 7.829 | `2.8777664` |  |
| 2 | 3.632 | `sqrt(initial_energy)` |  |
| 3 | 0.9332 | `noise_epsilon * 110.168755` |  |
| 5 | 0.7916 | `sqrt(sqrt(initial_energy * 91.53967))` |  |
| 6 | 0.6576 | `log((29.975407 * initial_energy) + 1.0806619)` |  |
| 8 | 0.386 | `abs(log(abs((n * initial_energy) - 1.2976338)))` |  |
| 10 | 0.2463 | `abs(log(abs((initial_energy * n) + -1.2938378)) * -0.9135673)` |  |
| 12 | 0.2127 | `abs(log(abs((n * initial_energy) - 1.2984009))) - (initial_energy * 0.043432042)` |  |
| 13 | 0.2126 | `abs(log(abs((initial_energy * n) - 1.2980788))) - abs(initial_energy * -0.042988848)` |  |
| 14 | 0.2082 | `(abs(log(abs(1.2749882 - (initial_energy * n)))) - (0.044887207 * initial_energy)) - 0.06988208` |  |
| 15 | 0.1958 | `abs(log(abs((initial_energy * n) - 1.2962997))) - abs(-0.14997852 - (-0.05244016 * initial_energy))` |  |
| 16 | 0.1958 | `abs(abs((initial_energy * -0.05244016) - -0.14997852) - abs(log(abs((n * initial_energy) - 1.2962997))))` |  |
| 17 | 0.1949 | `abs(log(abs((initial_energy * n) - 1.2968028))) - (abs(-0.14381775 - (initial_energy * -0.05200797)) + 0.022145815)` |  |
| 18 | 0.1916 | `abs(abs((-0.047207233 * initial_energy) - -0.13998367) - abs(log(abs(1.2964337 - (n * initial_energy))))) - 0.08523293` |  |
| 19 | 0.1916 | `abs(abs(abs(log(abs((n * initial_energy) - 1.2964337))) - abs((initial_energy * -0.047207233) - -0.13998367)) - 0.08523293)` | **★** |

## Experiment B — Warmup dynamics predictor

**Question**: can a formula built from warmup move statistics
(acceptance rate, energy change during warmup) and order-theoretic
invariants predict whether a near-truth initialization is preserved?

- Samples: 108
- Features (X): mm_dim, midpoint_dim, abs_discrepancy_mm_midpoint, chain3_abundance, n, target_dim, noise_epsilon, initial_energy, warmup_acceptance_rate, warmup_delta_energy
- Target (y): `preserved_near_truth` ∈ {0, 1}
- Rows: `truth_plus_small_noise` + `truth_plus_medium_noise`,
  warmup_mode ∈ {legacy_warmup, guarded_warmup}
- PySR iterations: 100
- `warmup_acceptance_rate` = accepted_moves / attempted_moves

### Discovered equations (Pareto front)

| complexity | loss | equation | best |
| ---: | ---: | --- | :---: |
| 1 | 0.2486 | `0.4629699` |  |
| 3 | 0.03446 | `0.00092555606 / noise_epsilon` |  |
| 5 | 0.02871 | `0.041143913 / (initial_energy + 0.040237892)` |  |
| 6 | 0.02851 | `square(0.08370819 / (initial_energy + 0.08257621))` |  |
| 7 | 0.0203 | `0.030685078 / ((initial_energy * abs_discrepancy_mm_midpoint) + 0.02979767)` |  |
| 8 | 0.01938 | `square(0.06068767 / ((abs_discrepancy_mm_midpoint * initial_energy) + 0.05938534))` |  |
| 9 | 0.01574 | `0.031480845 / ((abs_discrepancy_mm_midpoint * (initial_energy + warmup_delta_energy)) + 0.030225035)` |  |
| 10 | 0.01517 | `square(0.067162216 / ((abs_discrepancy_mm_midpoint * (initial_energy + warmup_delta_energy)) + 0.06592686))` |  |
| 12 | 0.01393 | `0.048731018 / ((((initial_energy + warmup_delta_energy) * abs_discrepancy_mm_midpoint) / log(midpoint_dim)) + 0.046959013)` |  |
| 13 | 0.01276 | `0.08785749 / ((((warmup_delta_energy + initial_energy) * abs_discrepancy_mm_midpoint) / square(log(midpoint_dim))) + 0.08496517)` |  |
| 15 | 0.01052 | `square(0.6297455 / (((warmup_delta_energy + initial_energy) * (abs_discrepancy_mm_midpoint / square(square(log(midpoint_dim))))) + 0.6219816))` |  |
| 16 | 0.01052 | `square(abs(0.6271797 / (((abs_discrepancy_mm_midpoint * (initial_energy + warmup_delta_energy)) / square(square(log(midpoint_dim)))) + 0.61946356)))` |  |
| 17 | 0.009885 | `square(0.36135474 / ((((initial_energy / midpoint_dim) + warmup_delta_energy) * (abs_discrepancy_mm_midpoint / square(square(log(midpoint_dim))))) + 0.35860553))` |  |
| 18 | 0.009599 | `square(0.5492915 / square(((warmup_delta_energy + (initial_energy / midpoint_dim)) * (abs_discrepancy_mm_midpoint / square(square(log(midpoint_dim))))) + 0.7379694))` |  |
| 19 | 0.009343 | `square(square(0.72691613 / (((((initial_energy / sqrt(target_dim)) + warmup_delta_energy) * abs_discrepancy_mm_midpoint) / square(square(log(midpoint_dim)))) + 0.7235194)))` |  |
| 20 | 0.007542 | `square(square(0.55567336 / ((((((initial_energy / midpoint_dim) + warmup_delta_energy) * initial_energy) / chain3_abundance) / square(square(log(midpoint_dim)))) + 0.5534699)))` | **★** |

## Data provenance

- Invariants: `phase1d_structural_atlas.csv`
- Dynamics: `phase2e_warmup_skip_probe.csv` + `phase2f_guarded_warmup_probe.csv`
- Join key: `(family, target_dim, n, seed)` — minkowski only.

Regenerate via `make regen-phase3a`.
Source: `tools/build_phase3a_pysr_warmup_rule.py`.
