# Phase 3D — PySR on the stratum-residual target

**Status:** exploratory diagnostic.  Not a physics result.

Phase 3C showed that with the binary target `preserved_near_truth`,
order-only panels could not improve over the majority-class constant.
Phase 3D keeps the dataset unchanged but switches to a continuous
target — the within-stratum residual of warmup-relative drift —
designed to expose whatever variability stratification cannot absorb.

## Verdict (automatic)

**SIGNAL_CONDITIONED_ON_DESIGN_VARIABLES**

Panel C improves but Panel D does not.  The signal needs target_dim and/or warmup_mode_code (design variables) even after removing noise_level.  Order features by themselves cannot explain the residual.

## Target definition

```
raw_relative_drift     = warmup_delta_energy / initial_energy
stratum_key            = (noise_level, n, target_dim)
stratum_mean           = mean(raw_relative_drift | stratum)
residual_relative_drift= raw_relative_drift - stratum_mean   ← target y
```

Rows with |initial_energy| < 0.0001 are excluded to avoid division blow-ups (documented; no row is silently kept).

## Sample sizes per stratum

| noise_level | n | target_dim | count | mean raw | min | max |
| :---: | ---: | :---: | ---: | ---: | ---: | ---: |
| small | 32 | 2 | 4 | +0.252 | -0.4724 | +1.48 |
| small | 32 | 3 | 4 | +0.2522 | -0.4828 | +1.125 |
| small | 32 | 4 | 2 | +0.8289 | -0.4474 | +2.105 |
| small | 64 | 2 | 6 | -0.2453 | -0.9667 | +1.731 |
| small | 64 | 3 | 4 | +11.94 | -0.6923 | +49.27 |
| small | 64 | 4 | 6 | -0.02227 | -0.68 | +1.868 |
| medium | 32 | 2 | 6 | +30.4 | -0.1317 | +76.76 |
| medium | 32 | 3 | 6 | +22.46 | -0.02678 | +48.39 |
| medium | 32 | 4 | 6 | +34.42 | -0.0722 | +89.97 |
| medium | 64 | 2 | 6 | +40.46 | +0 | +118.7 |
| medium | 64 | 3 | 6 | +48.27 | -0.008668 | +112.9 |
| medium | 64 | 4 | 6 | +130.3 | -0.01466 | +366.1 |

Total rows: 62.  Source: phase2f_guarded_warmup_probe.csv
(phase2e excluded because it lacks the `warmup_delta_energy` column).

Filter: init_label ∈ {truth_plus_small_noise, truth_plus_medium_noise}, warmup_mode ∈ {legacy_warmup, guarded_warmup}.

## Panel definitions

| panel | description | features |
| --- | --- | --- |
| A | noise-only baseline | `noise_level` |
| B | order + noise + warmup_mode + target_dim | `noise_level, warmup_mode_code, n, target_dim, mm_dim, midpoint_dim, abs_discrepancy_mm_mp, chain2_count, chain3_count, chain3_abundance, link_count, link_density, relation_count, ordering_fraction, height` |
| C | order + warmup_mode + target_dim (no noise_level) | `warmup_mode_code, n, target_dim, mm_dim, midpoint_dim, abs_discrepancy_mm_mp, chain2_count, chain3_count, chain3_abundance, link_count, link_density, relation_count, ordering_fraction, height` |
| D | order-only no-oracle (only n is a coarse design var) | `n, mm_dim, midpoint_dim, abs_discrepancy_mm_mp, chain2_count, chain3_count, chain3_abundance, link_count, link_density, relation_count, ordering_fraction, height` |

**Panel D leakage guard:** the following columns are excluded from
Panel D and verified at build time. If any reappears, the run fails:

  `delta_energy, final_distance_to_truth_rms, final_energy, final_interval_rmse, improved_energy, improved_interval_rmse, initial_distance_to_truth_rms, initial_energy, initial_interval_rmse, noise_level, preserved_near_truth, raw_relative_drift, residual_relative_drift, stratum_mean, target_dim, warmup_accepted_moves, warmup_attempted_moves, warmup_delta_energy, warmup_energy_after, warmup_energy_before, warmup_mode_code, warmup_rejected_moves`

## Summary: best loss per panel

Constant-predictor loss (variance of y): **2690.86**

| panel | description | best loss | Δ vs constant | rel. Δ |
| --- | --- | ---: | ---: | ---: |
| A | noise-only baseline | 2691 | +0.0005377 | +0.00% |
| B | order + noise + warmup_mode + target_dim | 462.8 | +2228 | +82.80% |
| C | order + warmup_mode + target_dim (no noise_level) | 1367 | +1324 | +49.19% |
| D | order-only no-oracle (only n is a coarse design var) | 2612 | +78.44 | +2.92% |

## Pareto fronts (per panel, complexity ≤ 12)

### Panel A — noise-only baseline

| complexity | loss | equation | best | order? | design? |
| ---: | ---: | --- | :---: | :---: | :---: |
| 1 | 2691 | `0.0087040765` |  | — | — |
| 5 | 2691 | `(noise_level - 0.31507555) * 0.011615529` | **★** | — | ✓ |

### Panel B — order + noise + warmup_mode + target_dim

| complexity | loss | equation | best | order? | design? |
| ---: | ---: | --- | :---: | :---: | :---: |
| 1 | 2691 | `chain3_abundance` |  | ✓ | — |
| 2 | 2690 | `log(midpoint_dim)` |  | ✓ | — |
| 3 | 2597 | `warmup_mode_code * -6.1343303` |  | — | ✓ |
| 4 | 2222 | `log(warmup_mode_code) * -46.182617` |  | — | ✓ |
| 5 | 1751 | `(122.643166 / warmup_mode_code) + -91.98233` |  | — | ✓ |
| 7 | 1176 | `noise_level * ((204.29231 / warmup_mode_code) + -153.21942)` |  | — | ✓ |
| 9 | 842.6 | `noise_level * (((109.75351 / warmup_mode_code) + -81.39457) * midpoint_dim)` |  | ✓ | ✓ |
| 10 | 609.2 | `(noise_level * square(midpoint_dim)) * ((51.701862 / warmup_mode_code) + -37.944366)` |  | ✓ | ✓ |
| 11 | 486.3 | `square(square(noise_level * midpoint_dim)) * ((8.443689 / warmup_mode_code) + -6.0887847)` |  | ✓ | ✓ |
| 12 | 462.8 | `(((15.883984 / warmup_mode_code) + -11.629286) * square(noise_level * midpoint_dim)) * target_dim` | **★** | ✓ | ✓ |

### Panel C — order + warmup_mode + target_dim (no noise_level)

| complexity | loss | equation | best | order? | design? |
| ---: | ---: | --- | :---: | :---: | :---: |
| 1 | 2691 | `ordering_fraction` |  | ✓ | — |
| 2 | 2690 | `log(midpoint_dim)` |  | ✓ | — |
| 3 | 2597 | `warmup_mode_code * -6.129034` |  | — | ✓ |
| 4 | 2222 | `log(warmup_mode_code) * -46.03808` |  | — | ✓ |
| 5 | 1680 | `(1.5000092 - warmup_mode_code) * n` |  | — | ✓ |
| 7 | 1581 | `(-47.617744 - (-64.11832 / warmup_mode_code)) * midpoint_dim` |  | ✓ | ✓ |
| 8 | 1412 | `square(midpoint_dim) * (target_dim / (1.4520153 - warmup_mode_code))` |  | ✓ | ✓ |
| 9 | 1404 | `(link_density * (2.5437698 - square(warmup_mode_code))) * square(target_dim)` |  | ✓ | ✓ |
| 10 | 1371 | `target_dim * ((square(midpoint_dim) / (1.3275819 - warmup_mode_code)) - target_dim)` |  | ✓ | ✓ |
| 11 | 1367 | `target_dim * ((square(midpoint_dim) / (1.3275819 - warmup_mode_code)) - log(n))` | **★** | ✓ | ✓ |

### Panel D — order-only no-oracle (only n is a coarse design var)

| complexity | loss | equation | best | order? | design? |
| ---: | ---: | --- | :---: | :---: | :---: |
| 1 | 2691 | `ordering_fraction` |  | ✓ | — |
| 2 | 2690 | `log(midpoint_dim)` |  | ✓ | — |
| 3 | 2672 | `-0.012344133 / chain3_abundance` |  | ✓ | — |
| 5 | 2648 | `-0.67332983 / (midpoint_dim + -3.0266519)` |  | ✓ | — |
| 7 | 2645 | `square(square(square(square(midpoint_dim + -1.7721493))))` |  | ✓ | — |
| 8 | 2638 | `((square(midpoint_dim) * 0.0132613415) - 0.057432607) / chain3_abundance` |  | ✓ | — |
| 9 | 2630 | `square(square(midpoint_dim) + -4.139362) - (0.018383298 / chain3_abundance)` |  | ✓ | — |
| 10 | 2617 | `square(square(square(midpoint_dim + -1.4778565))) - (0.016311547 / chain3_abundance)` |  | ✓ | — |
| 11 | 2614 | `square(square(square(square(midpoint_dim + -1.7624222)))) - (0.015887719 / chain3_abundance)` |  | ✓ | — |
| 12 | 2612 | `square(square(square(midpoint_dim + -1.3671405)) + -1.5998611) - (0.018984566 / chain3_abundance)` | **★** | ✓ | — |

## Decision rule (applied above to produce verdict)

Improvement threshold: a panel "improves clearly" if best loss is
at least 10% below the constant-predictor loss.

- **INTRA_STRATUM_ORDER_SIGNAL**: D improves AND best equation in D uses
  order features.  → autonomous order-theoretic signal.
- **SIGNAL_CONDITIONED_ON_DESIGN_VARIABLES**: C improves but D does not.
  → signal needs target_dim and/or warmup_mode_code.
- **SIGNAL_CONDITIONED_ON_NOISE_LEVEL**: B improves but C and D do not.
  → signal needs noise_level specifically.
- **NO_DETECTABLE_SIGNAL**: no panel improves meaningfully.
- **MIXED**: heterogeneous, no single conclusion.

## Reproducibility

- PySR niterations per panel: 100
- maxsize cap: 12
- random_state=1959, parallelism=serial (deterministic).
- initial_energy_floor: 0.0001

Regenerate via `make regen-phase3d`.
Source: `tools/build_phase3d_pysr_residual_target.py`.
