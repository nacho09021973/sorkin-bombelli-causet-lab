# N=24 Gamma-Grid Stress Test

This folder contains a pure exploratory gamma-grid stress test for
SORKIN-2/Bombelli annealer accessibility. It is not a confirmation run.

Protocol:

- `N = 24`
- `T0 = 100`
- budget: short Phase 2B schedule
- `warmup_limit = 10`
- `anneal_limit = 10`
- `max_data = 4`
- completed gamma-grid points: `193`
- grid families tested: `8`

Current summary:

- `final_energy` and `interval_rmse` agree in only `1` of `8` grid
  families.
- All `success_flag` values are `false`.
- `final_energy` prefers low gamma values, approximately `0.45` to
  `0.56`.
- `interval_rmse` prefers higher gamma values, approximately `0.72` to
  `0.82`.
- There is no robust gamma optimum in this exploratory stress test.

This result makes no physical resonance claim, no embeddability claim,
and no recovery claim.

## Current Interpretation

The result suggests metric instability rather than a gamma resonance.
Energy minimization and interval-structure preservation are not aligned
under this budget.

## Next Methodological Question

Before increasing `N` or refining gamma, decide which diagnostic should
be primary for SORKIN-2 recoverability: final energy, interval RMSE, MM
dimension recovery, success flag, or a composite score.

## Files

- `run_gamma_grid_stress_n24.py`: permanent exploratory runner.
- `gamma_grid_stress_n24.csv`: numerical stress-test results.
- `gamma_grid_stress_n24.md`: generated run summary.
- `gamma_grid_stress_n24_summary.svg`: generated SVG summary figure.
