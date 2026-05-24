# N=24 Gamma-Grid Stress Test

**Status:** exploratory only; not confirmation.

This stress test asks whether the apparent gamma preference in the N=24 probe is stable under changes of gamma grid, or whether it follows arbitrary grid nodes. It is an accessibility diagnostic for the historical Bombelli annealer only.

## Provenance

- Output directory: `explore/gamma_grid_stress_n24`
- Runner: `explore/gamma_grid_stress_n24/run_gamma_grid_stress_n24.py`
- Command: `python3 explore/gamma_grid_stress_n24/run_gamma_grid_stress_n24.py`
- Generated at UTC: `2026-05-24T08:25:02+00:00`
- CSV: `explore/gamma_grid_stress_n24/gamma_grid_stress_n24.csv`
- Figure: `explore/gamma_grid_stress_n24/gamma_grid_stress_n24_summary.svg`
- N: `24`
- T0: `100.0`
- warmup_limit: `10`
- anneal_limit: `10`
- max_data: `4`
- completed grid points: `193`
- total runtime seconds: `2.563`
- stopped early: `false`
- stop reason: `completed all requested grid families`

## Per-Family Winners

| grid_family | points | best final_energy gamma | final_energy | best interval_rmse gamma | interval_rmse | agree | successes |
| --- | ---: | ---: | ---: | ---: | ---: | :---: | ---: |
| log_0p5_0p9_21 | 21 | 0.546085796 | 130.446696 | 0.546085796 | 55649.304629 | yes | 0 |
| lin_0p5_0p9_21 | 21 | 0.560000000 | 130.465017 | 0.820000000 | 113319.222524 | no | 0 |
| log_0p45_0p95_21 | 21 | 0.450000000 | 125.577796 | 0.731382767 | 18585.571053 | no | 0 |
| lin_0p45_0p95_21 | 21 | 0.450000000 | 125.577796 | 0.725000000 | 18585.571053 | no | 0 |
| log_0p5_0p9_31 | 31 | 0.551461785 | 130.446696 | 0.725508644 | 18585.571053 | no | 0 |
| lin_0p5_0p9_31 | 31 | 0.553333333 | 130.465017 | 0.726666667 | 18585.571053 | no | 0 |
| jitter_log_0p5_0p9_21_seed24051959 | 21 | 0.560749555 | 130.465017 | 0.817363397 | 113319.222524 | no | 0 |
| custom_dense_0p75_0p9 | 26 | 0.540000000 | 131.092389 | 0.815625000 | 113319.222524 | no | 0 |

## Required Questions

1. Which gamma wins by `final_energy` for each grid family? See the per-family table above.
2. Which gamma wins by `interval_rmse` for each grid family? See the per-family table above.
3. Do winners cluster in a stable gamma region or follow arbitrary grid nodes?
   `final_energy` winners span `0.450000000` to `0.560749555`; `interval_rmse` winners span `0.546085796` to `0.820000000`. The winners are grid- and metric-dependent in this exploratory run.
4. Do `final_energy` and `interval_rmse` agree or disagree?
   They agree in `1` of `8` completed grid families.
5. Are any `success_flag` values true?
   `no`.
6. Is there any evidence for a robust gamma optimum?
   No robust optimum is established by this exploratory stress test.
7. Warning: no physical resonance claim, no embeddability claim, and no recovery claim.

## Interpretation Guardrails

- Exploration only.
- N=24 only.
- Fixed T0=100.
- Short Phase 2B budget only.
- Accessibility diagnostic only.
