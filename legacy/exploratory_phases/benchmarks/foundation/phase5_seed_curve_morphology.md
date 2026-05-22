# Phase 5 — Seed-level curve morphology audit

**Status:** exploratory seed-level audit using existing Phase 4B CSV provenance only. No new simulations, PySR, BDG action, SVG generation, or Phase 4B relabeling are performed.

## Objective

Phase 5 asks whether V-like or interior-minimum morphology observed in aggregate Phase 4B cell curves also appears in individual seed curves, or whether it is an artifact of averaging heterogeneous seeds.

## Inputs

- `phase4b_survival_probe_per_seed.csv`: reconstructs `loss(epsilon)` by `(n, target_dim, seed)`.
- `phase4b_survival_probe.csv`: supplies aggregate Phase 4B cell labels for comparison.
- `phase4b_survival_probe_per_epsilon.csv`: documents aggregate per-epsilon provenance.

## Seed-level classification rule

For each `(n, target_dim, seed)`, Phase 5 sorts valid rows by epsilon and classifies the individual `loss(epsilon)` curve. If fewer than three valid epsilons are available, the seed is marked `seed_insufficient_valid_points`.

Floor saturation uses the Phase 4B floor tolerance `1e-06`. A floor-saturated seed is censored and is not counted as strong negative evidence against V-like behavior.

The shape taxonomy is `seed_v_shape`, `seed_monotone_decay`, `seed_floor_saturated`, `seed_interior_min_noisy_tail`, `seed_flat_noisy`, and `seed_insufficient_valid_points`.

## Cell-level summary

| n | target_dim | seeds | seed_v_shape | seed_monotone_decay | seed_interior_min_noisy_tail | seed_floor_saturated | seed_insufficient_valid_points | majority_seed_shape | cell_curve_shape_phase4b | aggregate_represents_majority |
| ---: | :---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | :---: |
| 32 | 2 | 10 | 0 | 0 | 0 | 10 | 0 | seed_floor_saturated | monotone_decay | false |
| 32 | 3 | 10 | 0 | 0 | 0 | 10 | 0 | seed_floor_saturated | v_shape | false |
| 32 | 4 | 10 | 1 | 0 | 2 | 5 | 2 | seed_floor_saturated | v_shape | false |
| 48 | 2 | 10 | 0 | 0 | 0 | 10 | 0 | seed_floor_saturated | monotone_decay | false |
| 48 | 3 | 10 | 0 | 0 | 0 | 10 | 0 | seed_floor_saturated | v_shape | false |
| 48 | 4 | 10 | 0 | 0 | 0 | 9 | 1 | seed_floor_saturated | monotone_decay | false |
| 64 | 2 | 10 | 0 | 0 | 0 | 10 | 0 | seed_floor_saturated | noisy | false |
| 64 | 3 | 10 | 0 | 0 | 0 | 10 | 0 | seed_floor_saturated | monotone_decay | false |
| 64 | 4 | 10 | 0 | 0 | 0 | 9 | 1 | seed_floor_saturated | v_shape | false |

## Seed-level questions

- Aggregate V-shapes with at least one seed-level V-shape: 1/4.
- `(48,3)` remains a seed-level mixed/borderline audit target: 0 seed-level V-shapes, 0 interior-minimum noisy-tail seeds, majority `seed_floor_saturated`.
- `(48,4)` remains a counterexample/noisy-tail audit target: 0 seed-level V-shapes, 0 interior-minimum noisy-tail seeds, majority `seed_floor_saturated`.
- Cells with floor-saturated seeds: 9. These seeds are censored and are not counted as strong negative evidence against V-like behavior.

## Global Phase 5 outcome

**INSUFFICIENT**

Outcome definitions:

- `SEED_LEVEL_SUPPORT`: aggregate V-like morphology is broadly visible in individual seed curves.
- `SEED_LEVEL_MIXED`: seed-level morphology is present but heterogeneous across cells or seeds.
- `AGGREGATE_ARTIFACT`: aggregate V-like morphology is not visible at seed level.
- `INSUFFICIENT`: available seed curves are not sufficient for the audit.

## Floor-censoring audit

| n | target_dim | n_seeds_total | n_min_loss_eq_0 | n_le_1e-15 | n_le_1e-12 | n_le_1e-9 | n_le_1e-6 | n_le_1e-4 | min(min_loss) | median(min_loss) | max(min_loss) |
| ---: | :---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 32 | 2 | 10 | 10 | 10 | 10 | 10 | 10 | 10 | 0 | 0 | 0 |
| 32 | 3 | 10 | 10 | 10 | 10 | 10 | 10 | 10 | 0 | 0 | 0 |
| 32 | 4 | 10 | 5 | 5 | 5 | 5 | 5 | 5 | 0 | 0 | 0.07596641182 |
| 48 | 2 | 10 | 10 | 10 | 10 | 10 | 10 | 10 | 0 | 0 | 0 |
| 48 | 3 | 10 | 10 | 10 | 10 | 10 | 10 | 10 | 0 | 0 | 0 |
| 48 | 4 | 10 | 9 | 9 | 9 | 9 | 9 | 9 | 0 | 0 | 0 |
| 64 | 2 | 10 | 10 | 10 | 10 | 10 | 10 | 10 | 0 | 0 | 0 |
| 64 | 3 | 10 | 10 | 10 | 10 | 10 | 10 | 10 | 0 | 0 | 0 |
| 64 | 4 | 10 | 9 | 9 | 9 | 9 | 9 | 9 | 0 | 0 | 0 |

- Global valid seed curves with finite `min_loss`: 86.
- Exact-zero minima: 83/86.
- `min_loss <= 1e-15`: 83/86.
- `min_loss <= 1e-12`: 83/86.
- `min_loss <= 1e-9`: 83/86.
- `min_loss <= 1e-6`: 83/86.
- `min_loss <= 1e-4`: 83/86.

Seed-level morphology is currently censored by optimizer-floor saturation; aggregate Phase 4B curves should not be interpreted as direct evidence of seed-level V-shapes.

In the current pilot grid this censoring is dominated by exact zeros rather than by sensitivity to the chosen `floor_tolerance` threshold.

## Conservative conclusion

Phase 5 distinguishes seed-level morphology from aggregate morphology under the Phase 4B optimizer-response loss. It does not establish a physical law, validate Phase 4B as a physical claim, or replace the existing Phase 4B `MIXED` outcome.

Rows written: 90 seed curves across 9 pilot-grid cells.

Source: `tools/build_phase5_seed_curve_morphology.py`.
