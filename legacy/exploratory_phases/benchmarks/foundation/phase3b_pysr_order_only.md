# Phase 3B — PySR: Strictly Order-Theoretic Features

Critical test following Phase 3A: does the signal survive when
`noise_epsilon`, `initial_energy`, and all warmup-dynamic features
are removed?  Only causal-matrix observables and coarse design
variables are allowed as predictors.

## Setup

- Samples: 180 (84 preserved, 96 destroyed)
- PySR iterations: 100

**Input features (X):**

| feature | type |
| --- | --- |
| `n` | design |
| `target_dim` | design |
| `noise_level` | design (0=small, 1=medium) |
| `warmup_mode_code` | design (0=skip, 1=legacy, 2=guarded) |
| `mm_dim` | order-theoretic |
| `midpoint_dim` | order-theoretic |
| `abs_discrepancy_mm_mp` | order-theoretic |
| `chain2_count`, `chain3_count`, `chain3_abundance` | order-theoretic |
| `link_count`, `link_density` | order-theoretic |
| `relation_count`, `ordering_fraction` | order-theoretic |
| `height` | order-theoretic |

**Excluded vs Phase 3A:** `noise_epsilon`, `initial_energy`,
`warmup_delta_energy`, `warmup_acceptance_rate`, and all
embedding-dependent features.

**Target:** `preserved_near_truth` ∈ {0, 1}

**Key diagnostic:** equations flagged `uses_order_features=true` contain
at least one purely combinatorial causet observable.  If such equations
achieve lower loss than design-variable-only equations of comparable
complexity, the order-theoretic signal is real.

## Discovered equations (Pareto front)

| complexity | loss | equation | best | order features? |
| ---: | ---: | --- | :---: | :---: |
| 1 | 0.2489 | `0.46665585` |  | — |
| 3 | 0.03222 | `0.9666648 - noise_level` |  | — |
| 4 | 0.03111 | `square(0.9661014 - noise_level)` |  | — |
| 5 | 0.03111 | `(1.0000014 - noise_level) * 0.93332934` |  | — |
| 6 | 0.0299 | `sqrt(target_dim / mm_dim) - noise_level` |  | ✓ |
| 7 | 0.0268 | `square(noise_level - sqrt(target_dim / mm_dim))` |  | ✓ |
| 9 | 0.02613 | `square(noise_level - sqrt((target_dim / mm_dim) - 0.036516163))` |  | ✓ |
| 11 | 0.02435 | `square((noise_level + abs(0.980196 - (mm_dim / target_dim))) - 1.0550709)` |  | ✓ |
| 12 | 0.02414 | `square(square((noise_level + square(0.91809475 - (mm_dim / target_dim))) - 1.0061866))` |  | ✓ |
| 13 | 0.0213 | `square((noise_level + square((0.97227865 - (mm_dim / target_dim)) * link_density)) - 1.0204235)` |  | ✓ |
| 14 | 0.02004 | `square(square(square(link_density * (0.9934859 - (mm_dim / target_dim))) + (noise_level - 1.0232702)))` |  | ✓ |
| 15 | 0.01982 | `square((noise_level + square((link_density * (0.991036 - (mm_dim / target_dim))) / 0.7511943)) - 1.0422794)` |  | ✓ |
| 16 | 0.01606 | `square(noise_level + (square(square((link_density / 0.5411445) * (0.9811886 - (mm_dim / target_dim)))) - 1.0198423))` |  | ✓ |
| 17 | 0.01346 | `square(1.0076478 - (noise_level + square(square(square(link_density * ((0.9695584 - (mm_dim / target_dim)) / 0.48843363))))))` |  | ✓ |
| 18 | 0.01285 | `square(square(noise_level + (square(square(square(link_density * (((mm_dim / target_dim) - 0.96768034) / 0.524281)))) - 1.0046176)))` |  | ✓ |
| 20 | 0.01259 | `square((target_dim / target_dim) - (noise_level + square(square(square(square(((mm_dim / target_dim) - 0.97010595) * (link_density * 2.1930192)))))))` | **★** | ✓ |

## Interpretation guide

- If the best equation uses only `noise_level` and/or `warmup_mode_code`:
  the order-theoretic signal does **not** survive — the dataset size or
  the noise-level dominance is too strong.
- If `abs_discrepancy_mm_mp`, `ordering_fraction`, or any other
  order-theoretic feature appears in competitive equations: the signal
  **survives** leakage removal and warrants deeper study.

## Data provenance

- `phase1d_structural_atlas.csv` (mm_dim, midpoint_dim, chain stats)
- `invariants.json` (link_count, relation_count, ordering_fraction, height)
- `phase2e_warmup_skip_probe.csv` + `phase2f_guarded_warmup_probe.csv`
- Join key: `(target_dim, n, seed)` — minkowski only.

Regenerate via `make regen-phase3b`.
Source: `tools/build_phase3b_pysr_order_only.py`.
