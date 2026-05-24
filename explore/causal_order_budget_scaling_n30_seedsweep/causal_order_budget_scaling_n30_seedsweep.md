# N=30 Optimizer-Seed Budget/Misalignment Probe

exploratory only; not confirmation.

At N=30, fixed family/case/schedule/budgets, does the short -> medium energy/causal-metric misalignment observed at optimizer seed 1987 persist across optimizer seeds, or was it seed-specific?

## Provenance

- Output directory: `explore/causal_order_budget_scaling_n30_seedsweep`
- Runner: `explore/causal_order_budget_scaling_n30_seedsweep/run_causal_order_budget_scaling_n30_seedsweep.py`
- Command: `python3 explore/causal_order_budget_scaling_n30_seedsweep/run_causal_order_budget_scaling_n30_seedsweep.py`
- Generated at UTC: `2026-05-24T16:50:27+00:00`
- CSV: `explore/causal_order_budget_scaling_n30_seedsweep/causal_order_budget_scaling_n30_seedsweep.csv`
- Figure: `explore/causal_order_budget_scaling_n30_seedsweep/causal_order_budget_scaling_n30_seedsweep.svg`
- Log: `explore/causal_order_budget_scaling_n30_seedsweep/causal_order_budget_scaling_n30_seedsweep.log`
- N: `30`
- family: `minkowski`
- d_spacetime: `2`
- case seed: `1959`
- optimizer seeds: `1959, 1962, 1987, 2001, 2026`
- T0: `100.0`
- gamma: `0.5`
- h: `1`
- budgets: `short_10_10_4, medium_25_25_8`
- reused rows (from `explore/n_ladder_budget_misalignment_probe/n_ladder_budget_misalignment_probe.csv`): `seed 1987`
- newly run optimizer seeds: `1959, 1962, 2001, 2026`
- skipped rows: `none`
- per-run timeout seconds: `900.0`
- total wall runtime seconds: `798.823`
- Final causal-order diagnostics use `sim.rold`/`sim.xold`, not `rnew`/`xnew`.

## Result table

| seed | budget | status | final E | interval RMSE | F1 | recall | precision | missing | extra | exact | success | runtime s | note |
| ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: | --- |
| 1959 | short_10_10_4 | completed | 175.540541 | 14868.3868 | 0.240000 | 0.148760 | 0.620690 | 206 | 22 | no | no | 0.031 |  |
| 1959 | medium_25_25_8 | completed | 162.613029 | 7768.9424 | 0.534653 | 0.446281 | 0.666667 | 134 | 54 | no | no | 164.693 |  |
| 1962 | short_10_10_4 | completed | 163.598711 | 14589933.0884 | 0.390244 | 0.297521 | 0.566929 | 170 | 55 | no | no | 0.030 |  |
| 1962 | medium_25_25_8 | completed | 152.196622 | 211080.4171 | 0.457895 | 0.359504 | 0.630435 | 155 | 51 | no | no | 142.528 |  |
| 1987 | short_10_10_4 | completed | 192.900753 | 422233.2087 | 0.466488 | 0.359504 | 0.664122 | 155 | 44 | no | no | 0.025 | reused from explore/n_ladder_budget_misalignment_probe/n_ladder_budget_misalignment_probe.csv |
| 1987 | medium_25_25_8 | completed | 134.963611 | 5973.8217 | 0.385269 | 0.280992 | 0.612613 | 174 | 43 | no | no | 260.110 | reused from explore/n_ladder_budget_misalignment_probe/n_ladder_budget_misalignment_probe.csv |
| 2001 | short_10_10_4 | completed | 168.559313 | 22446.3107 | 0.386167 | 0.276860 | 0.638095 | 175 | 38 | no | no | 0.088 |  |
| 2001 | medium_25_25_8 | completed | 140.934924 | 989.3000 | 0.766316 | 0.752066 | 0.781116 | 60 | 51 | no | no | 267.045 |  |
| 2026 | short_10_10_4 | completed | 189.295404 | 9923801.7088 | 0.338192 | 0.239669 | 0.574257 | 184 | 43 | no | no | 0.031 |  |
| 2026 | medium_25_25_8 | completed | 144.752144 | 1810.9720 | 0.619647 | 0.508264 | 0.793548 | 119 | 32 | no | no | 224.258 |  |

## Seed-level comparison (short -> medium)

| seed | E_short | E_medium | delta_E | F1_short | F1_medium | delta_F1 | recall_short | recall_medium | missing_short | missing_medium | label |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 1959 | 175.540541 | 162.613029 | -12.927511 | 0.240000 | 0.534653 | 0.294653 | 0.148760 | 0.446281 | 206 | 134 | aligned_improvement |
| 1962 | 163.598711 | 152.196622 | -11.402089 | 0.390244 | 0.457895 | 0.067651 | 0.297521 | 0.359504 | 170 | 155 | aligned_improvement |
| 1987 | 192.900753 | 134.963611 | -57.937142 | 0.466488 | 0.385269 | -0.081219 | 0.359504 | 0.280992 | 155 | 174 | energy_causality_misalignment |
| 2001 | 168.559313 | 140.934924 | -27.624389 | 0.386167 | 0.766316 | 0.380149 | 0.276860 | 0.752066 | 175 | 60 | aligned_improvement |
| 2026 | 189.295404 | 144.752144 | -44.543260 | 0.338192 | 0.619647 | 0.281455 | 0.239669 | 0.508264 | 184 | 119 | aligned_improvement |

Interpretation labels (thresholds: |delta_E| > 1.0, |delta_F1| > 0.02):
- aligned_improvement: delta_E < -1.0 and delta_F1 > +0.02.
- energy_causality_misalignment: delta_E < -1.0 and delta_F1 < -0.02.
- causal_improves_without_energy: delta_E > +1.0 and delta_F1 > +0.02.
- no_clear_change: otherwise.

## Decision rule

- seeds with energy_causality_misalignment: `1/5` (1987).
- seeds with aligned_improvement: `4/5` (1959, 1962, 2001, 2026).
- seeds with causal_improves_without_energy: `0/5` (none).
- seeds with no_clear_change: `0/5` (none).
- rule: >=4/5 misalignment -> persists; <=2/5 -> seed-specific; otherwise mixed/inconclusive.
- outcome: `misalignment_appears_seed_specific`.

## Conservative readout

- Does the N=30 short->medium misalignment persist across optimizer seeds at this case seed/family/schedule? `no`.
- Is the N=24 -> N=30 bracket strengthened or weakened by this probe? `weakened (N=30 misalignment looks seed-specific in this case seed)`.
- Does this justify a topology / case-seed follow-up at N=30? `low priority for topology probe before broadening optimizer seeds at the closer N values`.
- This probe does not establish: a transition in N, a physical gamma claim, an embeddability claim, a theorem, or a general N-scaling statement.

## Guardrails

- one N only (`30`).
- one case seed only (`1959`).
- one family only (`minkowski`).
- one schedule only (`T0=100`, `gamma=0.5`, `h=1`).
- two budgets only (`short_10_10_4`, `medium_25_25_8`).
- exploratory only.
- no embeddability claim.
- no physical gamma claim.
- no theorem.
- no general N-scaling claim.
- optimizer seed 1987 rows are reused from `explore/n_ladder_budget_misalignment_probe/n_ladder_budget_misalignment_probe.csv` and were not regenerated here.
