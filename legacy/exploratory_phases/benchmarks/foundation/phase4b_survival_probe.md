# Phase 4B — Exploratory survival probe

**Status:** exploratory survival test of the Phase 4A post-hoc pattern. No PySR is run here.

## Objective

Phase 4B asks whether a V-like aggregate optimizer-response morphology is concentrated in cells with high `dim_discrepancy_rel_midpoint` and without floor saturation. This is not a confirmatory physical claim.

## Grid design

- Grid mode: `pilot`
- Sizes: 32, 48, 64
- Target spacetime dims: 2, 3, 4
- Epsilons: 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2
- Seeds: 1900, 1916, 1923, 1939, 1953, 1973, 1981, 1995, 2003, 2020
- Threshold theta: 0.35
- Floor tolerance: 1e-06
- Wall-clock: 79s

The `pilot` grid is the default Makefile target. The `full` grid is available explicitly via `--grid full`.

## Provenance note

`phase4b_survival_probe.csv` is the one-row-per-cell aggregate used for the global exploratory outcome. `phase4b_survival_probe_per_epsilon.csv` persists the per-epsilon curve values (`mean_loss = mean(abs_relative_drift)`) used to audit curve morphology. `phase4b_survival_probe_per_seed.csv` persists one row per seed/epsilon run so later visual audits can select representative, lowest-loss and highest-loss runs before drawing Hasse diagrams. These files add provenance only; no thresholds, classifications, or physical conclusions are changed.

## Curve-shape summary

| n | target_dim | curve_shape | survival_label | dim_disc_rel | min_val | epsilon_at_min | rise_frac | ordering_fraction | chain3_abundance |
| ---: | :---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 32 | 2 | monotone_decay | control_negative | 0.135 | 0.001637 | 0.06 | 0.06411 | 0.4847 | 0.1554 |
| 32 | 3 | v_shape | supporting | 0.5711 | 0.02538 | 0.08 | 0.1219 | 0.2177 | 0.02347 |
| 32 | 4 | v_shape | supporting | 1.51 | 0.06894 | 0.1 | 0.1388 | 0.1305 | 0.004637 |
| 48 | 2 | monotone_decay | control_negative | 0.06299 | 0 | 0.06 | 0 | 0.5087 | 0.1704 |
| 48 | 3 | v_shape | counterexample | 0.2218 | 0.0002744 | 0.1 | 0.05781 | 0.2284 | 0.02377 |
| 48 | 4 | monotone_decay | counterexample | 0.9227 | 0.01772 | 0.06 | 0.07284 | 0.1214 | 0.004169 |
| 64 | 2 | noisy | control_negative | 0.07348 | 0 | 0.04 | 0 | 0.4916 | 0.1585 |
| 64 | 3 | monotone_decay | ambiguous | 0.1526 | 0 | 0.08 | 0.1227 | 0.2214 | 0.02296 |
| 64 | 4 | v_shape | supporting | 0.6799 | 0.002943 | 0.06 | 0.05508 | 0.116 | 0.003947 |

## Tail-cleanliness / borderline audit

This post-hoc audit records whether a curve has an interior minimum, how cleanly the right tail rises, and whether the V-like behavior is marginal. It does not change `curve_shape`, `survival_label`, or the global outcome.

| n | target_dim | curve_shape | survival_label | tail_pattern | rise_frac | rise_frac_margin | borderline_v_like |
| ---: | :---: | --- | --- | --- | ---: | ---: | :---: |
| 32 | 2 | monotone_decay | control_negative | positive,negative,negative,positive | 0.06411 | 0.01411 | true |
| 32 | 3 | v_shape | supporting | positive,positive,positive | 0.1219 | 0.07189 | false |
| 32 | 4 | v_shape | supporting | positive,positive | 0.1388 | 0.0888 | false |
| 48 | 2 | monotone_decay | control_negative | positive,negative,positive,negative | 0 | -0.05 | false |
| 48 | 3 | v_shape | counterexample | positive,positive | 0.05781 | 0.00781 | true |
| 48 | 4 | monotone_decay | counterexample | positive,negative,positive,negative | 0.07284 | 0.02284 | true |
| 64 | 2 | noisy | control_negative | zero,zero,zero,zero,zero | 0 | -0.05 | false |
| 64 | 3 | monotone_decay | ambiguous | positive,negative,positive | 0.1227 | 0.0727 | false |
| 64 | 4 | v_shape | supporting | positive,positive,positive,positive | 0.05508 | 0.005084 | true |

## Survival test of Phase 4A hypothesis

- V-shape curves: 4
- Supporting cells: 3
- Counterexamples: 2
- Censored floor cases: 0

`survival_label` is an exploratory interpretation, not a per-cell physical verdict. Floor-saturated high-discrepancy non-V curves are classified as `censored_floor` and are not counted as strong negative evidence.

## Negative controls target_dim=2

| n | curve_shape | survival_label | min_val | dim_disc_rel |
| ---: | --- | --- | ---: | ---: |
| 32 | monotone_decay | control_negative | 0.001637 | 0.135 |
| 48 | monotone_decay | control_negative | 0 | 0.06299 |
| 64 | noisy | control_negative | 0 | 0.07348 |

## Floor-saturated / censored cases

| n | target_dim | curve_shape | survival_label | min_val | epsilon_at_min |
| ---: | :---: | --- | --- | ---: | ---: |
| 48 | 2 | monotone_decay | control_negative | 0 | 0.06 |
| 64 | 2 | noisy | control_negative | 0 | 0.04 |
| 64 | 3 | monotone_decay | ambiguous | 0 | 0.08 |

## Loss semantics

In Phase 4B, `loss` is inherited directly from the Phase 4A `abs_relative_drift` column: `loss = |warmup_delta_energy / initial_energy|` for valid runs. Here `initial_energy` is the simulator energy after initializing the embedding from the sprinkled coordinates plus epsilon-scaled coordinate noise, and `warmup_delta_energy` is the change produced by the guarded warmup stage before the later anneal result is used.

`loss` is therefore an optimizer/embedding response diagnostic under this specific energy, initialization, epsilon, seed, and coordinate parametrization. It is not an intrinsic observable of the partial order, not the dimensional discrepancy, not a Lorentz-invariant distance to Minkowski space, and not an absolute physical quality score for the causet.

The visual audit roles `best`, `near_mean`, and `worst` rank seeds only by this Phase 4B loss within the same cell and epsilon. They are provenance labels for selecting examples, not physical labels for good or bad causal sets.

## Global exploratory outcome

**MIXED**

Outcome definitions:

- `PASS_EXPLORATORY_SURVIVAL`: most V-shapes occur in high-discrepancy, non-floor cells; target_dim=2 controls do not produce strong false positives; floor cases are censored rather than counted as clean negatives.
- `MIXED`: partial signal with enough counterexamples or false positives to require refinement.
- `FAIL`: the high-discrepancy plus no-floor separation does not survive the expanded grid.

## Conservative conclusion

Phase 4B shows partial survival of the Phase 4A pattern, but counterexamples or ambiguous cells require refinement before any stronger claim.

Regenerate via `make regen-phase4b` for the pilot grid, or run `python3 tools/build_phase4b_survival_probe.py --grid full` explicitly for the full grid.
Source: `tools/build_phase4b_survival_probe.py`.
