# Phase 3E — PySR on warmup-mode-aware stratum residuals

**Status:** methodological cleanup test, not a physics result.

Phase 3D residualized by `(noise_level, n, target_dim)` but the
residuals were dominated by `warmup_mode_code` — legacy vs guarded
warmup differ by construction and should have been stratified out.

Phase 3E uses the corrected stratification
`(noise_level, n, target_dim, warmup_mode)`, so the residual
captures only seed-to-seed variation within an otherwise fixed
protocol.  Different seeds give different sprinklings with different
order-theoretic invariants — this is the variance order features
should explain if any signal exists at this dataset scale.

## Verdict (automatic)

**POSSIBLE_INTRA_PROTOCOL_ORDER_SIGNAL**

Panel D improves by +15.0% over constant and uses order features in its best equation.  After controlling for all design variables via stratification, order-theoretic features explain part of the seed-to-seed residual.  Exploratory positive; not a physical law and needs replication with more data.

## Target definition

```
raw_relative_drift     = warmup_delta_energy / initial_energy
stratum_key            = (noise_level, n, target_dim, warmup_mode)
stratum_mean           = mean(raw_relative_drift | stratum)
residual_relative_drift= raw_relative_drift - stratum_mean   ← target y
```

Rows with `|initial_energy| < 0.0001` are excluded to avoid division blow-ups.

Singleton strata (count = 1) are dropped: the residual would be identically zero and carry no within-stratum information.

## Strata

- Rows after E0 filter: 62
- Strata: 24
- Singleton strata dropped: 2
- Rows kept for regression: 60
- Target std (kept rows): 20.32
- Constant-predictor loss (variance of y): **412.711**

### Per-stratum counts

| noise | n | target_dim | warmup_mode | count | mean raw | min | max | kept? |
| :---: | ---: | :---: | --- | ---: | ---: | ---: | ---: | :---: |
| small | 32 | 2 | guarded_warmup | 2 | -0.4421 | -0.4724 | -0.4118 | ✓ |
| small | 32 | 2 | legacy_warmup | 2 | +0.946 | +0.4118 | +1.48 | ✓ |
| small | 32 | 3 | guarded_warmup | 2 | -0.4289 | -0.4828 | -0.375 | ✓ |
| small | 32 | 3 | legacy_warmup | 2 | +0.9332 | +0.7414 | +1.125 | ✓ |
| small | 32 | 4 | guarded_warmup | 1 | — | — | — | — |
| small | 32 | 4 | legacy_warmup | 1 | — | — | — | — |
| small | 64 | 2 | guarded_warmup | 3 | -0.6306 | -0.75 | -0.4615 | ✓ |
| small | 64 | 2 | legacy_warmup | 3 | +0.1399 | -0.9667 | +1.731 | ✓ |
| small | 64 | 3 | guarded_warmup | 2 | -0.4045 | -0.4359 | -0.3731 | ✓ |
| small | 64 | 3 | legacy_warmup | 2 | +24.29 | -0.6923 | +49.27 | ✓ |
| small | 64 | 4 | guarded_warmup | 3 | -0.4528 | -0.68 | -0.2632 | ✓ |
| small | 64 | 4 | legacy_warmup | 3 | +0.4082 | -0.6038 | +1.868 | ✓ |
| medium | 32 | 2 | guarded_warmup | 3 | -0.04694 | -0.1317 | +0 | ✓ |
| medium | 32 | 2 | legacy_warmup | 3 | +60.85 | +47.44 | +76.76 | ✓ |
| medium | 32 | 3 | guarded_warmup | 3 | -0.008928 | -0.02678 | +0 | ✓ |
| medium | 32 | 3 | legacy_warmup | 3 | +44.93 | +41.09 | +48.39 | ✓ |
| medium | 32 | 4 | guarded_warmup | 3 | -0.0287 | -0.0722 | +0 | ✓ |
| medium | 32 | 4 | legacy_warmup | 3 | +68.86 | +38.11 | +89.97 | ✓ |
| medium | 64 | 2 | guarded_warmup | 3 | +0 | +0 | +0 | ✓ |
| medium | 64 | 2 | legacy_warmup | 3 | +80.92 | +53.83 | +118.7 | ✓ |
| medium | 64 | 3 | guarded_warmup | 3 | -0.002889 | -0.008668 | +0 | ✓ |
| medium | 64 | 3 | legacy_warmup | 3 | +96.55 | +80.78 | +112.9 | ✓ |
| medium | 64 | 4 | guarded_warmup | 3 | -0.004885 | -0.01466 | +0 | ✓ |
| medium | 64 | 4 | legacy_warmup | 3 | +260.7 | +175.6 | +366.1 | ✓ |

## Panel definitions

| panel | description | features |
| --- | --- | --- |
| A | design-only sanity (no order features) | `noise_level, warmup_mode_code, n, target_dim` |
| B | order + design (noise, warmup, n, target_dim) | `noise_level, warmup_mode_code, n, target_dim, mm_dim, midpoint_dim, abs_discrepancy_mm_mp, chain2_count, chain3_count, chain3_abundance, link_count, link_density, relation_count, ordering_fraction, height` |
| C | order + n + target_dim (no noise, no warmup) | `n, target_dim, mm_dim, midpoint_dim, abs_discrepancy_mm_mp, chain2_count, chain3_count, chain3_abundance, link_count, link_density, relation_count, ordering_fraction, height` |
| D | order-only no-oracle (only n is a coarse design var) | `n, mm_dim, midpoint_dim, abs_discrepancy_mm_mp, chain2_count, chain3_count, chain3_abundance, link_count, link_density, relation_count, ordering_fraction, height` |

**Panel D leakage guard** (verified at build time):

  `delta_energy, final_distance_to_truth_rms, final_energy, final_interval_rmse, improved_energy, improved_interval_rmse, initial_distance_to_truth_rms, initial_energy, initial_interval_rmse, noise_level, preserved_near_truth, raw_relative_drift, residual_relative_drift, stratum_mean, target_dim, warmup_accepted_moves, warmup_attempted_moves, warmup_delta_energy, warmup_energy_after, warmup_energy_before, warmup_mode, warmup_mode_code, warmup_rejected_moves`

## Summary: best loss per panel

| panel | description | best loss | Δ vs constant | rel. Δ |
| --- | --- | ---: | ---: | ---: |
| A | design-only sanity (no order features) | 412.7 | +5.912e-05 | +0.00% |
| B | order + design (noise, warmup, n, target_dim) | 253.5 | +159.3 | +38.59% |
| C | order + n + target_dim (no noise, no warmup) | 354.6 | +58.07 | +14.07% |
| D | order-only no-oracle (only n is a coarse design var) | 351 | +61.73 | +14.96% |

## Pareto fronts (per panel, complexity ≤ 10)

### Panel A — design-only sanity (no order features)

| complexity | loss | equation | best | order? | design? |
| ---: | ---: | --- | :---: | :---: | :---: |
| 1 | 412.7 | `0.0029067164` | **★** | — | — |

### Panel B — order + design (noise, warmup, n, target_dim)

| complexity | loss | equation | best | order? | design? |
| ---: | ---: | --- | :---: | :---: | :---: |
| 1 | 412.6 | `ordering_fraction` |  | ✓ | — |
| 2 | 411.3 | `log(midpoint_dim)` |  | ✓ | — |
| 3 | 392.4 | `-0.012524564 / chain3_abundance` |  | ✓ | — |
| 4 | 381.2 | `-1.3747413e-5 / square(chain3_abundance)` |  | ✓ | — |
| 5 | 345.2 | `noise_level / (abs_discrepancy_mm_mp - 0.7246527)` |  | ✓ | ✓ |
| 7 | 256.1 | `noise_level / (abs_discrepancy_mm_mp - (0.7341591 * warmup_mode_code))` |  | ✓ | ✓ |
| 9 | 254 | `(noise_level / (abs_discrepancy_mm_mp - (0.7341591 * warmup_mode_code))) / warmup_mode_code` |  | ✓ | ✓ |
| 10 | 253.5 | `(noise_level / (abs_discrepancy_mm_mp - (0.7341591 * warmup_mode_code))) / square(warmup_mode_code)` | **★** | ✓ | ✓ |

### Panel C — order + n + target_dim (no noise, no warmup)

| complexity | loss | equation | best | order? | design? |
| ---: | ---: | --- | :---: | :---: | :---: |
| 1 | 412.6 | `ordering_fraction` |  | ✓ | — |
| 2 | 411.3 | `log(midpoint_dim)` |  | ✓ | — |
| 3 | 392.4 | `-0.01253002 / chain3_abundance` |  | ✓ | — |
| 4 | 381.2 | `-1.3747293e-5 / square(chain3_abundance)` |  | ✓ | — |
| 5 | 381.2 | `0.29902807 / (link_density + -1.7936398)` |  | ✓ | — |
| 6 | 374.9 | `(-1.3747293e-5 / square(chain3_abundance)) + midpoint_dim` |  | ✓ | — |
| 7 | 368.4 | `square(square(square(square(square(log(midpoint_dim))))))` |  | ✓ | — |
| 9 | 362.1 | `square(square(square(square(square(log(midpoint_dim)))))) - abs_discrepancy_mm_mp` |  | ✓ | — |
| 10 | 354.6 | `square(square(square(square(square(log(midpoint_dim)))))) - square(abs_discrepancy_mm_mp)` | **★** | ✓ | — |

### Panel D — order-only no-oracle (only n is a coarse design var)

| complexity | loss | equation | best | order? | design? |
| ---: | ---: | --- | :---: | :---: | :---: |
| 1 | 412.6 | `ordering_fraction` |  | ✓ | — |
| 2 | 411.3 | `log(midpoint_dim)` |  | ✓ | — |
| 3 | 392.4 | `-0.012524527 / chain3_abundance` |  | ✓ | — |
| 5 | 381.7 | `(-0.019458313 / chain3_abundance) + mm_dim` |  | ✓ | — |
| 6 | 369.9 | `square(square(square(midpoint_dim - 1.5082731)))` |  | ✓ | — |
| 7 | 365.2 | `square(square(square(square(midpoint_dim - 1.7721802))))` |  | ✓ | — |
| 8 | 362.5 | `square(square(square(1.5010043 - midpoint_dim))) - abs_discrepancy_mm_mp` |  | ✓ | — |
| 9 | 355.1 | `square(square(square(1.5010043 - midpoint_dim))) - square(abs_discrepancy_mm_mp)` |  | ✓ | — |
| 10 | 351 | `square(square(square(square(1.7700238 - midpoint_dim)))) - square(abs_discrepancy_mm_mp)` | **★** | ✓ | — |

## Decision rule

Improvement thresholds:
- *strong* = best loss is at least **10%** below constant baseline
- *null*   = best loss is less than **5%** below constant baseline

Applied in this order:

1. **RESIDUALIZATION_FAILED_OR_LEAKAGE**: A improves strongly.
   → residualization or feature set is broken; do not interpret.
2. **POSSIBLE_INTRA_PROTOCOL_ORDER_SIGNAL**: D improves strongly AND
   uses order features.  → exploratory positive.
3. **NULL_SIGNAL**: D's improvement below null threshold.
   → no detectable order signal at this dataset scale.
4. **SIGNAL_CONDITIONED_ON_DESIGN_VARIABLES**: B strong, D not strong.
5. **INTERMEDIATE**: D between null and strong thresholds.

## Reproducibility

- PySR iterations per panel: 100
- maxsize cap: 10
- random_state=1959, parallelism=serial (deterministic).
- initial_energy_floor: 0.0001

Regenerate via `make regen-phase3e`.
Source: `tools/build_phase3e_pysr_residual_by_warmup_mode.py`.
