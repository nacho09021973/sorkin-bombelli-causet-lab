# Phase 1D Structural Atlas

Ensemble statistics for the Phase 1C finite-size grid with
one additional order-theoretic observable: 3-chain abundance.
The raw CSV records `chain2_count`, `chain3_count`, and
`chain3_abundance = chain3_count / binom(n, 3)`.

This is not a calibrated dimension estimator. It is a minimal
structural statistic asking whether higher-order chain
abundance helps separate Minkowski sprinklings from
Kleitman-Rothschild and suspended corona controls without
coordinates, energy, embedding, fitting, or optimization.

| family | target d | n | seeds | mean MM | std MM | mean midpoint | std midpoint | mean \|disc\| | mean C3 abundance | std C3 abundance |
| --- | :---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| corona_poset | - | 32 | 5 | 1.8809 | 0.0000 | 4.0000 | 0.0000 | 2.1191 | 0.0907 | 0.0000 |
| corona_poset | - | 64 | 5 | 1.9390 | 0.0000 | 5.0000 | 0.0000 | 3.0610 | 0.0461 | 0.0000 |
| corona_poset | - | 128 | 5 | 1.9691 | 0.0000 | 6.0000 | 0.0000 | 4.0309 | 0.0233 | 0.0000 |
| corona_poset | - | 256 | 5 | 1.9845 | 0.0000 | 7.0000 | 0.0000 | 5.0155 | 0.0117 | 0.0000 |
| kleitman_rothschild | - | 32 | 5 | 2.3431 | 0.0491 | 2.2582 | 0.1166 | 0.1684 | 0.0521 | 0.0058 |
| kleitman_rothschild | - | 64 | 5 | 2.3611 | 0.0211 | 3.0846 | 0.0906 | 0.7235 | 0.0489 | 0.0023 |
| kleitman_rothschild | - | 128 | 5 | 2.3651 | 0.0157 | 3.8549 | 0.0936 | 1.4898 | 0.0480 | 0.0018 |
| kleitman_rothschild | - | 256 | 5 | 2.3719 | 0.0096 | 4.7094 | 0.0783 | 2.3374 | 0.0474 | 0.0011 |
| minkowski | 2 | 32 | 5 | 1.9301 | 0.1567 | 1.6546 | 0.1729 | 0.2754 | 0.1843 | 0.0501 |
| minkowski | 2 | 64 | 5 | 2.0144 | 0.1253 | 1.9074 | 0.1951 | 0.1574 | 0.1634 | 0.0364 |
| minkowski | 2 | 128 | 5 | 2.0364 | 0.0915 | 1.9553 | 0.1581 | 0.1315 | 0.1585 | 0.0251 |
| minkowski | 2 | 256 | 5 | 2.0206 | 0.0763 | 2.0580 | 0.1456 | 0.0704 | 0.1611 | 0.0202 |
| minkowski | 3 | 32 | 5 | 3.3164 | 0.4448 | 1.7428 | 0.2602 | 1.5735 | 0.0187 | 0.0160 |
| minkowski | 3 | 64 | 5 | 3.2151 | 0.1692 | 2.3550 | 0.4130 | 0.8855 | 0.0169 | 0.0071 |
| minkowski | 3 | 128 | 5 | 3.1101 | 0.1463 | 2.7175 | 0.2712 | 0.3979 | 0.0195 | 0.0056 |
| minkowski | 3 | 256 | 5 | 3.0642 | 0.1119 | 2.9978 | 0.3060 | 0.2830 | 0.0206 | 0.0045 |
| minkowski | 4 | 32 | 5 | 3.7065 | 0.3431 | 1.5546 | 0.4623 | 2.1519 | 0.0042 | 0.0026 |
| minkowski | 4 | 64 | 5 | 3.9174 | 0.3646 | 2.4853 | 0.3930 | 1.4321 | 0.0032 | 0.0015 |
| minkowski | 4 | 128 | 5 | 3.9532 | 0.2287 | 3.3465 | 0.3875 | 0.6067 | 0.0030 | 0.0013 |
| minkowski | 4 | 256 | 5 | 4.0673 | 0.1537 | 3.8041 | 0.2794 | 0.4229 | 0.0023 | 0.0007 |

Reading by family:

- Minkowski d=2: high relation density and high C3 abundance;
  the observable is expected to overlap corona controls at
  larger n because both sit near effective dimension two by
  ordering fraction.
- Minkowski d=3,d=4: C3 abundance falls rapidly with dimension
  and gives an independent structural scale beside MM and
  midpoint.
- Kleitman-Rothschild: intermediate C3 abundance with growing
  MM/midpoint discrepancy.
- Corona: high C3 abundance but a sharply growing midpoint
  discrepancy, so C3 alone is not a manifoldness classifier.

Regenerate via `make regen-phase1d`. Source tool:
`tools/build_phase1d_structural_atlas.py`.

For interpretation see the *Phase 1D* section of
`results_note_2026.md`.
