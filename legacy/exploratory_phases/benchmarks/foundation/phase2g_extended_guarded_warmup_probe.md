# Phase 2G — Extended Guarded-Warmup Probe

Data foundation for the Phase 3F PySR ablation.

- Total runs: 180
- Preserved (delta_energy ≤ 1e-6): 82
- Seeds (15): 1959, 1962, 1987, 2009, 2026, 1812, 1848, 1871, 1905, 1929, 1945, 1968, 1989, 2001, 2017
- Sizes: 32, 64
- Spacetime dims: 2, 3, 4
- Init labels: truth_plus_small_noise, truth_plus_medium_noise
- Warmup mode: guarded_warmup only (guard threshold = 0.0, warmup_limit = 10)
- Anneal: anneal_limit=10, max_data=4, T₀=100.0, γ=0.9, optimizer_seed=1987

## Preservation by (d, n, init)

| d | n | init | runs | preserved | mean ΔE | mean ΔE/E₀ |
| :---: | ---: | --- | ---: | ---: | ---: | ---: |
| 2 | 32 | truth_plus_medium_noise | 15 | 0/15 | 119.3253 | 20.4112 |
| 2 | 32 | truth_plus_small_noise | 15 | 15/15 | -0.0020 | -0.9050 |
| 2 | 64 | truth_plus_medium_noise | 15 | 0/15 | 484.6693 | 20.2724 |
| 2 | 64 | truth_plus_small_noise | 15 | 9/15 | 46.1239 | 4948.2508 |
| 3 | 32 | truth_plus_medium_noise | 15 | 0/15 | 98.9847 | 22.6197 |
| 3 | 32 | truth_plus_small_noise | 15 | 15/15 | -0.0019 | -0.8493 |
| 3 | 64 | truth_plus_medium_noise | 15 | 0/15 | 363.9072 | 20.6118 |
| 3 | 64 | truth_plus_small_noise | 15 | 13/15 | 0.6722 | 65.4745 |
| 4 | 32 | truth_plus_medium_noise | 15 | 0/15 | 82.7324 | 28.0480 |
| 4 | 32 | truth_plus_small_noise | 15 | 15/15 | -0.0012 | -0.8774 |
| 4 | 64 | truth_plus_medium_noise | 15 | 0/15 | 274.6512 | 24.1318 |
| 4 | 64 | truth_plus_small_noise | 15 | 15/15 | -0.0036 | -0.7396 |

## Warmup statistics (guarded-warmup only)

| d | n | init | mean attempted | mean accepted | mean rejected | mean W ΔE |
| :---: | ---: | --- | ---: | ---: | ---: | ---: |
| 2 | 32 | truth_plus_medium_noise | 10.00 | 2.20 | 7.80 | -0.2220 |
| 2 | 32 | truth_plus_small_noise | 5.33 | 2.60 | 2.73 | -0.0015 |
| 2 | 64 | truth_plus_medium_noise | 10.00 | 1.53 | 8.47 | 0.0000 |
| 2 | 64 | truth_plus_small_noise | 9.33 | 4.47 | 4.87 | -0.0046 |
| 3 | 32 | truth_plus_medium_noise | 10.00 | 2.40 | 7.60 | -0.3587 |
| 3 | 32 | truth_plus_small_noise | 9.20 | 4.47 | 4.73 | -0.0015 |
| 3 | 64 | truth_plus_medium_noise | 10.00 | 1.47 | 8.53 | -0.0407 |
| 3 | 64 | truth_plus_small_noise | 9.33 | 4.53 | 4.80 | -0.0042 |
| 4 | 32 | truth_plus_medium_noise | 10.00 | 3.53 | 6.47 | -0.3842 |
| 4 | 32 | truth_plus_small_noise | 3.33 | 1.53 | 1.80 | -0.0009 |
| 4 | 64 | truth_plus_medium_noise | 10.00 | 2.13 | 7.87 | -0.2068 |
| 4 | 64 | truth_plus_small_noise | 9.33 | 3.80 | 5.53 | -0.0022 |

Regenerate via `make regen-phase2g`.
Source: `tools/build_phase2g_extended_guarded_warmup_probe.py`.
