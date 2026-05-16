# Phase 2 Embedding Bridge

A minimal bridge between Phase 1D order-theoretic diagnostics
and the historical Bombelli-Sorkin annealing code. This is a
fixed, small probe rather than an optimizer search.

Protocol:

- cases: n=64, case seed=1959; Minkowski d=2,3,4,
  Kleitman-Rothschild, and suspended corona.
- optimizer seed: 1987.
- schedule: warmup_limit=10, anneal_limit=10,
  max_data=4, initial_temp=100.0,
  cooling_factor=0.9.
- controls have no ground-truth coordinates, so `truth_energy`,
  `energy_gap`, and `interval_rmse` are recorded as NA.
- controls are embedded at spatial dim=2;
  Minkowski cases use spatial dim=d_spacetime-1.

| family | d | embed dim | MM | midpoint | \|disc\| | C3 abundance | final E | truth E | gap | RMSE | status |
| --- | :---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| minkowski | 2 | 1 | 1.9934 | 1.8981 | 0.0953 | 0.1659 | 821.8675 | 0.0000 | 821.8675 | 77567276.6828 | completed |
| minkowski | 3 | 2 | 3.2596 | 2.5025 | 0.7571 | 0.0141 | 637.2047 | 0.0000 | 637.2047 | 40941762.2346 | completed |
| minkowski | 4 | 3 | 3.7626 | 2.7004 | 1.0622 | 0.0044 | 661.6395 | 0.0000 | 661.6395 | 17146618.9655 | completed |
| kleitman_rothschild | - | 2 | 2.3362 | 3.0875 | 0.7513 | 0.0521 | 866.4288 | NA | NA | NA | completed |
| corona_poset | - | 2 | 1.9390 | 5.0000 | 3.0610 | 0.0461 | 952.6428 | NA | NA | NA | completed |

Interpretation:

- This probe is designed to expose qualitative alignment or
  tension between structural diagnostics and annealing outcomes.
- A low final energy alone is not treated as proof of a faithful
  embedding; for Minkowski cases the useful checks are also
  truth energy, energy gap, and interval RMSE.
- If a Minkowski case has good structural diagnostics but a large
  gap or RMSE, that points to optimizer/schedule failure rather
  than immediate non-manifoldlikeness.
- If a control reaches a deceptively low final energy, that is a
  warning that the energy can reward non-geometric artifacts.

Regenerate via `make regen-phase2`. Source tool:
`tools/build_phase2_embedding_bridge.py`.
