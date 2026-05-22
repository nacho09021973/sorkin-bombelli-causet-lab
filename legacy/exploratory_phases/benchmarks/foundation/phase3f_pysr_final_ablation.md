# Phase 3F — PySR final ablation on expanded dataset

**Status:** decisive test of the Phase 3E candidate signal on the
Phase 2G expanded dataset.  Not a physics law.  An exploratory
diagnostic with sufficient samples to either replicate or reject
the 3E `abs_discrepancy_mm_midpoint` finding.

## Verdict (automatic)

**INTERMEDIATE**

Panel D's improvement is +6.2% — between null (5%) and strong (10%) thresholds.  Weak; not interpretable as positive.

## Hypotheses tested

- **H1**: Panel D crosses +10% over constant baseline.
- **H2**: `abs_discrepancy_mm_midpoint` appears in D's (or C's) best equation.
- **H3**: Panel A (design-only sanity) is below +10% — residualization is clean.

## Dataset

- Source: `phase2g_extended_guarded_warmup_probe.csv` (15 seeds, 3 sizes, 3 dims, 2 near-truth init labels, guarded warmup only).
- Invariants: `phase1e_extended_structural_atlas.csv`.
- Rows after E₀ filter (|E₀| ≥ 0.0001): 158
- Strata: 12, singletons dropped: 0
- Rows kept for regression: 158
- Residual std: 0.1521
- Constant-predictor loss (variance of y): **0.0231248**

## Stratum sizes

| noise | n | target_dim | count | mean raw | min | max | kept? |
| :---: | ---: | :---: | ---: | ---: | ---: | ---: | :---: |
| small | 32 | 2 | 8 | -0.6968 | -0.875 | -0.4118 | ✓ |
| small | 32 | 3 | 13 | -0.6532 | -1 | -0.2857 | ✓ |
| small | 32 | 4 | 5 | -0.6607 | -0.8182 | -0.4474 | ✓ |
| small | 64 | 2 | 14 | -0.5483 | -0.8202 | -0.1188 | ✓ |
| small | 64 | 3 | 14 | -0.48 | -0.9765 | -0.05405 | ✓ |
| small | 64 | 4 | 14 | -0.3993 | -0.7529 | -0.05714 | ✓ |
| medium | 32 | 2 | 15 | -0.03728 | -0.1317 | +0 | ✓ |
| medium | 32 | 3 | 15 | -0.0735 | -0.1899 | +0 | ✓ |
| medium | 32 | 4 | 15 | -0.1329 | -0.4298 | +0 | ✓ |
| medium | 64 | 2 | 15 | +0 | +0 | +0 | ✓ |
| medium | 64 | 3 | 15 | -0.002407 | -0.01389 | +0 | ✓ |
| medium | 64 | 4 | 15 | -0.0171 | -0.1006 | +0 | ✓ |

## Panel definitions

| panel | description | features |
| --- | --- | --- |
| A | design-only sanity (noise_level + n + target_dim) | `noise_level, n, target_dim` |
| B | order + design (noise_level + n + target_dim + order features) | `noise_level, n, target_dim, mm_dim, midpoint_dim, abs_discrepancy_mm_mp, chain2_count, chain3_count, chain3_abundance, chain4_count, link_count, link_density, relation_count, ordering_fraction, height` |
| C | order + known-d (n + target_dim + order features, no noise) | `n, target_dim, mm_dim, midpoint_dim, abs_discrepancy_mm_mp, chain2_count, chain3_count, chain3_abundance, chain4_count, link_count, link_density, relation_count, ordering_fraction, height` |
| D | order-only no-oracle (n + order features only) | `n, mm_dim, midpoint_dim, abs_discrepancy_mm_mp, chain2_count, chain3_count, chain3_abundance, chain4_count, link_count, link_density, relation_count, ordering_fraction, height` |

**Panel D leakage guard** (verified at build time):

  `delta_energy, final_distance_to_truth_rms, final_energy, final_interval_rmse, improved_energy, improved_interval_rmse, initial_distance_to_truth_rms, initial_energy, initial_interval_rmse, noise_level, preserved_near_truth, raw_relative_drift, residual_relative_drift, stratum_mean, target_dim, warmup_accepted_moves, warmup_attempted_moves, warmup_delta_energy, warmup_energy_after, warmup_energy_before, warmup_mode, warmup_mode_code, warmup_rejected_moves`

## Summary: best loss per panel

| panel | description | best loss | Δ vs constant | rel. Δ |
| --- | --- | ---: | ---: | ---: |
| A | design-only sanity (noise_level + n + target_dim) | 0.02312 | +1.298e-08 | +0.00% |
| B | order + design (noise_level + n + target_dim + order features) | 0.02136 | +0.001768 | +7.64% |
| C | order + known-d (n + target_dim + order features, no noise) | 0.02204 | +0.001086 | +4.70% |
| D | order-only no-oracle (n + order features only) | 0.02168 | +0.001441 | +6.23% |

## Pareto fronts (per panel, complexity ≤ 10)

### Panel A — design-only sanity (noise_level + n + target_dim)

| complexity | loss | equation | best | order? | design? | abs_disc? |
| ---: | ---: | --- | :---: | :---: | :---: | :---: |
| 1 | 0.02312 | `-1.6880264e-05` | **★** | — | — | — |

### Panel B — order + design (noise_level + n + target_dim + order features)

| complexity | loss | equation | best | order? | design? | abs_disc? |
| ---: | ---: | --- | :---: | :---: | :---: | :---: |
| 1 | 0.02312 | `-1.6923079e-5` |  | — | — | — |
| 3 | 0.0231 | `-0.000692426 / ordering_fraction` |  | ✓ | — | — |
| 4 | 0.02305 | `-0.00012320875 / square(ordering_fraction)` |  | ✓ | — | — |
| 5 | 0.0225 | `-1.4022067 / (relation_count - chain4_count)` |  | ✓ | — | — |
| 7 | 0.02209 | `(noise_level - 1.9634908) / (chain2_count - chain4_count)` |  | ✓ | ✓ | — |
| 9 | 0.02136 | `((midpoint_dim * noise_level) - link_density) / (chain2_count - chain4_count)` | **★** | ✓ | ✓ | — |

### Panel C — order + known-d (n + target_dim + order features, no noise)

| complexity | loss | equation | best | order? | design? | abs_disc? |
| ---: | ---: | --- | :---: | :---: | :---: | :---: |
| 1 | 0.02312 | `-1.6891574e-5` |  | — | — | — |
| 3 | 0.0231 | `-0.00069250586 / ordering_fraction` |  | ✓ | — | — |
| 4 | 0.02305 | `-0.16303729 / square(height)` |  | ✓ | — | — |
| 5 | 0.0225 | `-1.4022048 / (relation_count - chain4_count)` |  | ✓ | — | — |
| 6 | 0.02234 | `log(link_density) / (chain4_count - relation_count)` |  | ✓ | — | — |
| 7 | 0.02227 | `link_density * (0.5388285 / (chain4_count - chain2_count))` |  | ✓ | — | — |
| 8 | 0.02216 | `-1.5957824 / (chain2_count - square(chain4_count + -0.71434915))` |  | ✓ | — | — |
| 9 | 0.02211 | `-1.595925 / (chain2_count - square(square(chain4_count + -0.7092614)))` |  | ✓ | — | — |
| 10 | 0.02204 | `link_density / (-3.3345299 - ((chain2_count - square(chain4_count)) + chain4_count))` | **★** | ✓ | — | — |

### Panel D — order-only no-oracle (n + order features only)

| complexity | loss | equation | best | order? | design? | abs_disc? |
| ---: | ---: | --- | :---: | :---: | :---: | :---: |
| 1 | 0.02312 | `-2.4020537e-6` |  | — | — | — |
| 3 | 0.0231 | `-0.00069250556 / ordering_fraction` |  | ✓ | — | — |
| 4 | 0.02282 | `0.006839401 / log(link_density)` |  | ✓ | — | — |
| 5 | 0.02265 | `0.066558406 / (chain3_count - 0.26373217)` |  | ✓ | — | — |
| 6 | 0.0226 | `0.06995433 / (square(chain3_count) - 0.271868)` |  | ✓ | — | — |
| 7 | 0.02208 | `0.004561078 / ((link_density * 1.9845618) - height)` |  | ✓ | — | — |
| 9 | 0.02202 | `(abs_discrepancy_mm_mp * -0.0015563451) / log(mm_dim / square(link_density))` |  | ✓ | — | ✓ |
| 10 | 0.02168 | `(square(abs_discrepancy_mm_mp) * -0.0014540266) / log(mm_dim / square(link_density))` | **★** | ✓ | — | ✓ |

## Decision rule

Thresholds: strong = 10%, null = 5% relative improvement over constant.

1. **RESIDUALIZATION_FAILED_OR_LEAKAGE**: Panel A above strong → broken.
2. **POSSIBLE_ORDER_SIGNAL_WITH_ABS_DISC**: D above strong, uses
   abs_discrepancy_mm_midpoint or C's best does.  → 3E candidate replicates.
3. **ORDER_SIGNAL_WITHOUT_ABS_DISC**: D above strong but abs_disc absent.
4. **NULL_SIGNAL**: D below null threshold.
5. **SIGNAL_CONDITIONED_ON_DESIGN**: B strong, D not.
6. **INTERMEDIATE**: D between thresholds.

## Reproducibility

- PySR iterations per panel: 100
- maxsize cap: 10
- random_state=1959, parallelism=serial, initial_energy_floor=0.0001

Regenerate via `make regen-phase3f`.
Source: `tools/build_phase3f_pysr_final_ablation.py`.
