# Phase 1C Two-Control Scaling Atlas

Ensemble statistics of two order-theoretic dimension
estimators across `n` in {32, 64, 128, 256} for canonical
Minkowski-diamond sprinklings (d_spacetime in {2, 3, 4}) and
two non-manifoldlike controls: Kleitman-Rothschild three-level
posets and suspended corona/crown posets. Each cell aggregates
five seeds: 1959, 1962, 1987, 2009, 2026.

Columns:

- `family`, `d`, `n`, `seeds`: the cell descriptor.
- `mean_mm`, `std_mm`: Myrheim-Meyer dim, ensemble mean and
  population standard deviation.
- `mean_midpoint`, `std_midpoint`: Meyer's midpoint scaling
  dim, same convention.
- `mean |disc|`: ensemble mean of `|mm - midpoint|`.

What to look for:

- Manifoldlike (Minkowski): `mean_mm` should converge toward
  `d_spacetime` as `n` grows, and the two estimators should
  move toward agreement at larger `n`.
- Non-manifoldlike controls: `mean_mm` should be comparatively
  flat within a family, while `mean_midpoint` grows with
  finite-size scale rather than converging to the same target.
- Agreement of KR and corona trends would support the claim that
  the Phase 1B signature is not KR-specific.

| family | d | n | seeds | mean MM | std MM | mean midpoint | std midpoint | mean \|disc\| |
| --- | :---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| corona_poset | - | 32 | 5 | 1.8809 | 0.0000 | 4.0000 | 0.0000 | 2.1191 |
| corona_poset | - | 64 | 5 | 1.9390 | 0.0000 | 5.0000 | 0.0000 | 3.0610 |
| corona_poset | - | 128 | 5 | 1.9691 | 0.0000 | 6.0000 | 0.0000 | 4.0309 |
| corona_poset | - | 256 | 5 | 1.9845 | 0.0000 | 7.0000 | 0.0000 | 5.0155 |
| kleitman_rothschild | - | 32 | 5 | 2.3431 | 0.0491 | 2.2582 | 0.1166 | 0.1684 |
| kleitman_rothschild | - | 64 | 5 | 2.3611 | 0.0211 | 3.0846 | 0.0906 | 0.7235 |
| kleitman_rothschild | - | 128 | 5 | 2.3651 | 0.0157 | 3.8549 | 0.0936 | 1.4898 |
| kleitman_rothschild | - | 256 | 5 | 2.3719 | 0.0096 | 4.7094 | 0.0783 | 2.3374 |
| minkowski | 2 | 32 | 5 | 1.9301 | 0.1567 | 1.6546 | 0.1729 | 0.2754 |
| minkowski | 2 | 64 | 5 | 2.0144 | 0.1253 | 1.9074 | 0.1951 | 0.1574 |
| minkowski | 2 | 128 | 5 | 2.0364 | 0.0915 | 1.9553 | 0.1581 | 0.1315 |
| minkowski | 2 | 256 | 5 | 2.0206 | 0.0763 | 2.0580 | 0.1456 | 0.0704 |
| minkowski | 3 | 32 | 5 | 3.3164 | 0.4448 | 1.7428 | 0.2602 | 1.5735 |
| minkowski | 3 | 64 | 5 | 3.2151 | 0.1692 | 2.3550 | 0.4130 | 0.8855 |
| minkowski | 3 | 128 | 5 | 3.1101 | 0.1463 | 2.7175 | 0.2712 | 0.3979 |
| minkowski | 3 | 256 | 5 | 3.0642 | 0.1119 | 2.9978 | 0.3060 | 0.2830 |
| minkowski | 4 | 32 | 5 | 3.7065 | 0.3431 | 1.5546 | 0.4623 | 2.1519 |
| minkowski | 4 | 64 | 5 | 3.9174 | 0.3646 | 2.4853 | 0.3930 | 1.4321 |
| minkowski | 4 | 128 | 5 | 3.9532 | 0.2287 | 3.3465 | 0.3875 | 0.6067 |
| minkowski | 4 | 256 | 5 | 4.0673 | 0.1537 | 3.8041 | 0.2794 | 0.4229 |

Regenerate via `make regen-phase1c`. Source tool:
`tools/build_phase1c_scaling_atlas.py`.

For interpretation see the *Phase 1C* section of
`results_note_2026.md`.
