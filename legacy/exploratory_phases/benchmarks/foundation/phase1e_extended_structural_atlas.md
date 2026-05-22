# Phase 1E — Extended Structural Atlas

Self-contained invariant table for the Phase 3F PySR ablation.
Replaces the join across `phase1d_structural_atlas.csv` +
`invariants.json` with a single CSV that covers the expanded
grid (15 seeds, n ∈ {32, 64, 128}, d ∈ {2, 3, 4}) and includes
every invariant the Phase 3 panels need.

Minkowski family only — controls (KR, corona) belong to Phase 1D.

- Total rows: 135
- Seeds (15): 1959, 1962, 1987, 2009, 2026, 1812, 1848, 1871, 1905, 1929, 1945, 1968, 1989, 2001, 2017
- Sizes: 32, 64, 128
- Spacetime dims: 2, 3, 4

## Ensemble statistics by (d, n)

| d | n | seeds | mean MM | std MM | mean midpoint | std midpoint | mean \|disc\| | std \|disc\| | mean link_density | mean ordering_frac |
| :---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2 | 32 | 15 | 2.0318 | 0.1781 | 1.8003 | 0.2299 | 0.2527 | 0.1764 | 2.1792 | 0.4919 |
| 2 | 64 | 15 | 2.0139 | 0.1330 | 1.9281 | 0.1920 | 0.1396 | 0.0755 | 2.8156 | 0.4969 |
| 2 | 128 | 15 | 2.0375 | 0.0761 | 2.0065 | 0.1364 | 0.1057 | 0.0610 | 3.4734 | 0.4868 |
| 3 | 32 | 15 | 3.1819 | 0.3633 | 1.9415 | 0.4109 | 1.2404 | 0.6650 | 2.0167 | 0.2051 |
| 3 | 64 | 15 | 3.0942 | 0.1875 | 2.3843 | 0.4607 | 0.7764 | 0.4941 | 3.2635 | 0.2141 |
| 3 | 128 | 15 | 3.0476 | 0.1497 | 2.7390 | 0.2606 | 0.3393 | 0.3544 | 5.1297 | 0.2214 |
| 4 | 32 | 15 | 3.8674 | 0.4623 | 1.6440 | 0.4718 | 2.0849 | 0.6585 | 1.5167 | 0.1195 |
| 4 | 64 | 15 | 3.8924 | 0.2870 | 2.3186 | 0.4681 | 1.5738 | 0.6576 | 2.6396 | 0.1125 |
| 4 | 128 | 15 | 4.0013 | 0.1940 | 3.0656 | 0.3746 | 0.9357 | 0.4655 | 4.1750 | 0.1012 |

Regenerate via `make regen-phase1e`.
Source: `tools/build_phase1e_extended_structural_atlas.py`.
