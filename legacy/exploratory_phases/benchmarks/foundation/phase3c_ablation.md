# Phase 3C — Ablation: Four Feature Panels

Controlled ablation to determine whether order-theoretic features
carry signal independent of experimental design variables.

## Setup

- Samples: 180 (84 preserved, 96 destroyed)
- PySR iterations per panel: 100
- maxsize: 12  ← complexity cap to suppress opportunistic fitting

## Panel definitions

| panel | features | noise_level | target_dim | order features |
| --- | --- | :---: | :---: | :---: |
| A | noise-only baseline | ✓ | — | — |
| B | order + noise | ✓ | ✓ | ✓ |
| C | order-only known-d | — | ✓ | ✓ |
| D | order-only no-oracle | — | — | ✓ |

## Summary: best loss per panel

| panel | description | best loss | Δ vs constant |
| --- | --- | ---: | ---: |
| A | noise-only baseline — no order features, no target_dim | 0.03111 | +0.4667 |
| B | order + noise — all order features plus noise_level and warmup_mode | 0.01741 | +0.4804 |
| C | order-only known-d — order features + target_dim, no noise_level | 0.2453 | +0.2525 |
| D | order-only no-oracle — fully autonomous (n + order features only) | 0.247 | +0.2508 |

Constant-predictor baseline loss (majority class): ~0.4978

**Decision rule:**
- If D's Δ vs constant is meaningful and order features appear in D's best equation:
  → genuine order-theoretic signal, warrants more data.
- If only A/B achieve large Δ: signal depends on experimental design variables.
- If all panels are near baseline: dataset too small or noise dominates completely.

## Panel A — noise-only baseline — no order features, no target_dim

Features: `noise_level`

| complexity | loss | equation | best | order? | design? |
| ---: | ---: | --- | :---: | :---: | :---: |
| 1 | 0.2489 | `0.46665585` |  | — | — |
| 3 | 0.03222 | `0.9666716 - noise_level` |  | — | ✓ |
| 4 | 0.03111 | `square(noise_level - 0.96610177)` |  | — | ✓ |
| 5 | 0.03111 | `(noise_level * -0.9333492) + 0.9333378` | **★** | — | ✓ |

## Panel B — order + noise — all order features plus noise_level and warmup_mode

Features: `noise_level, warmup_mode_code, mm_dim, midpoint_dim, abs_discrepancy_mm_mp, chain2_count, chain3_count, chain3_abundance, link_count, link_density, relation_count, ordering_fraction, height`

| complexity | loss | equation | best | order? | design? |
| ---: | ---: | --- | :---: | :---: | :---: |
| 1 | 0.2489 | `0.4666556` |  | — | — |
| 3 | 0.03222 | `0.9666648 - noise_level` |  | — | ✓ |
| 4 | 0.03111 | `square(0.9661015 - noise_level)` |  | — | ✓ |
| 5 | 0.03111 | `0.9333341 - (noise_level * 0.93333405)` |  | — | ✓ |
| 7 | 0.03035 | `(0.9980096 - noise_level) * (chain3_abundance + 0.8787318)` |  | ✓ | ✓ |
| 8 | 0.02965 | `square(((abs_discrepancy_mm_mp / -32.10843) + 1.0005388) - noise_level)` |  | ✓ | ✓ |
| 9 | 0.02516 | `((-3.0833972 / (chain2_count - chain3_count)) - noise_level) + 1.0226815` |  | ✓ | ✓ |
| 10 | 0.01835 | `square(((-3.0601213 / (relation_count - chain3_count)) + 1.0200946) - noise_level)` |  | ✓ | ✓ |
| 11 | 0.01783 | `square((-3.2800174 / abs(relation_count - chain3_count)) + (1.0315465 - noise_level))` |  | ✓ | ✓ |
| 12 | 0.01741 | `square(((midpoint_dim - height) / (chain2_count - chain3_count)) + (1.0204413 - noise_level))` | **★** | ✓ | ✓ |

## Panel C — order-only known-d — order features + target_dim, no noise_level

Features: `n, target_dim, mm_dim, midpoint_dim, abs_discrepancy_mm_mp, chain2_count, chain3_count, chain3_abundance, link_count, link_density, relation_count, ordering_fraction, height`

| complexity | loss | equation | best | order? | design? |
| ---: | ---: | --- | :---: | :---: | :---: |
| 1 | 0.2489 | `0.46665585` |  | — | — |
| 4 | 0.2485 | `square(chain3_abundance) + 0.45833024` |  | ✓ | — |
| 5 | 0.2477 | `log(target_dim) * sqrt(ordering_fraction)` |  | ✓ | ✓ |
| 6 | 0.2465 | `(ordering_fraction + 0.22811761) * log(target_dim)` |  | ✓ | ✓ |
| 7 | 0.2461 | `(log(sqrt(target_dim)) + ordering_fraction) - 0.3212177` |  | ✓ | ✓ |
| 8 | 0.2459 | `log(ordering_fraction + ((target_dim * 0.19099434) + 0.76419413))` |  | ✓ | ✓ |
| 9 | 0.2456 | `log(log((target_dim - (ordering_fraction / -0.18843183)) / 0.88146645))` |  | ✓ | ✓ |
| 10 | 0.2454 | `sqrt(log(log(target_dim - (ordering_fraction / -0.1942946))) + -0.15338896)` |  | ✓ | ✓ |
| 12 | 0.2453 | `sqrt(log(log((target_dim - (ordering_fraction / -0.18680938)) + -0.40594018)) / 1.4314072)` | **★** | ✓ | ✓ |

## Panel D — order-only no-oracle — fully autonomous (n + order features only)

Features: `n, mm_dim, midpoint_dim, abs_discrepancy_mm_mp, chain2_count, chain3_count, chain3_abundance, link_count, link_density, relation_count, ordering_fraction, height`

| complexity | loss | equation | best | order? | design? |
| ---: | ---: | --- | :---: | :---: | :---: |
| 1 | 0.2489 | `0.466627` |  | — | — |
| 4 | 0.2485 | `square(chain3_abundance) + 0.45833045` |  | ✓ | — |
| 5 | 0.248 | `sqrt(sqrt(link_count) / n)` |  | ✓ | — |
| 7 | 0.2479 | `sqrt(sqrt(link_density / (abs_discrepancy_mm_mp + n)))` |  | ✓ | — |
| 8 | 0.2479 | `sqrt(sqrt(sqrt(midpoint_dim * link_density) / n))` |  | ✓ | — |
| 9 | 0.2477 | `sqrt(sqrt(link_density / ((-0.020058928 / chain3_abundance) + n)))` |  | ✓ | — |
| 11 | 0.2473 | `sqrt(sqrt((link_density + -0.37415913) / ((-0.03155158 / chain3_abundance) + n)))` |  | ✓ | — |
| 12 | 0.247 | `sqrt(sqrt(sqrt(link_density / (n + (-0.024405533 / chain3_abundance))) + -0.18019497))` | **★** | ✓ | — |

## Data provenance

- `phase1d_structural_atlas.csv`, `invariants.json`
- `phase2e_warmup_skip_probe.csv`, `phase2f_guarded_warmup_probe.csv`
- Join key: `(target_dim, n, seed)` — minkowski only.

Regenerate via `make regen-phase3c`.
Source: `tools/build_phase3c_ablation.py`.
