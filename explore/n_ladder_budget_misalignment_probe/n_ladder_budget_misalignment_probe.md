# N-Ladder Budget/Misalignment Probe

exploratory only; not confirmation.

Across a small ladder of N values at fixed family/seed/h/gamma/T0, when does increasing the historical annealer budget from short to medium improve causal-order recovery, and when does it lower final energy while worsening causal-order recovery?

## Provenance

- Output directory: `explore/n_ladder_budget_misalignment_probe`
- Runner: `explore/n_ladder_budget_misalignment_probe/run_n_ladder_budget_misalignment_probe.py`
- Command: `python3 explore/n_ladder_budget_misalignment_probe/run_n_ladder_budget_misalignment_probe.py`
- Generated at UTC: `2026-05-24T16:28:13+00:00`
- CSV: `explore/n_ladder_budget_misalignment_probe/n_ladder_budget_misalignment_probe.csv`
- Figure: `explore/n_ladder_budget_misalignment_probe/n_ladder_budget_misalignment_probe.svg`
- N values requested: `12, 18, 30`
- Completed N values: `12, 18, 30`
- Skipped N values: `none`
- Completed short budget rows: `N=12, N=18, N=30`
- Completed medium budget rows: `N=12, N=18, N=30`
- T0: `100.0`
- gamma: `0.5`
- h: `1`
- case seed: `1959`
- optimizer seed: `1987`
- family: `minkowski`
- d_spacetime: `2`
- budgets attempted: `short_10_10_4, medium_25_25_8`
- per-run timeout seconds: `900.0`
- total wall runtime seconds: `260.749`
- N=24 reference CSV: `explore/causal_order_budget_scaling_n24/causal_order_budget_scaling_n24.csv` (rows found for optimizer seed 1987: long_50_50_16, medium_25_25_8, short_10_10_4)
- N=36 reference CSV: `explore/causal_order_budget_scaling_n36/causal_order_budget_scaling_n36.csv` (rows found for optimizer seed 1987: medium_25_25_8, short_10_10_4)
- Final causal-order diagnostics use `sim.rold`/`sim.xold`, not `rnew`/`xnew`.

## Result table by N and budget

| N | source | budget | status | final E | interval RMSE | F1 | recall | precision | missing | extra | exact | success | runtime s |
| ---: | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: |
| 12 | probe | short_10_10_4 | completed | 25.979490 | 12080.6904 | 0.192308 | 0.121951 | 0.454545 | 36 | 6 | no | no | 0.005 |
| 12 | probe | medium_25_25_8 | completed | 18.678446 | 60.9942 | 0.693333 | 0.634146 | 0.764706 | 15 | 8 | no | no | 0.033 |
| 18 | probe | short_10_10_4 | completed | 65.250277 | 12807.0993 | 0.449275 | 0.356322 | 0.607843 | 56 | 20 | no | no | 0.009 |
| 18 | probe | medium_25_25_8 | completed | 56.433963 | 33.1356 | 0.643836 | 0.540230 | 0.796610 | 40 | 12 | no | no | 0.541 |
| 24 | reference | short_10_10_4 | completed | 131.822227 | 71268273.4800 | 0.354978 | 0.259494 | 0.561644 | 117 | 32 | no | no | 0.015 |
| 24 | reference | medium_25_25_8 | completed | 80.989247 | 61.1911 | 0.712871 | 0.683544 | 0.744828 | 50 | 37 | no | no | 25.936 |
| 30 | probe | short_10_10_4 | completed | 192.900753 | 422233.2087 | 0.466488 | 0.359504 | 0.664122 | 155 | 44 | no | no | 0.025 |
| 30 | probe | medium_25_25_8 | completed | 134.963611 | 5973.8217 | 0.385269 | 0.280992 | 0.612613 | 174 | 43 | no | no | 260.110 |
| 36 | reference | short_10_10_4 | completed | 251.424456 | 413612.8960 | 0.488971 | 0.387755 | 0.661692 | 210 | 68 | no | no | 0.059 |
| 36 | reference | medium_25_25_8 | completed | 215.628980 | 664950.5619 | 0.373225 | 0.268222 | 0.613333 | 251 | 58 | no | no | 688.715 |

## Derived comparison table (short -> medium)

| N | source | E_short | E_medium | delta_E | F1_short | F1_medium | delta_F1 | recall_short | recall_medium | missing_short | missing_medium | interpretation |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 12 | probe | 25.979490 | 18.678446 | -7.301044 | 0.192308 | 0.693333 | 0.501026 | 0.121951 | 0.634146 | 36 | 15 | aligned_improvement |
| 18 | probe | 65.250277 | 56.433963 | -8.816313 | 0.449275 | 0.643836 | 0.194560 | 0.356322 | 0.540230 | 56 | 40 | aligned_improvement |
| 24 | reference | 131.822227 | 80.989247 | -50.832981 | 0.354978 | 0.712871 | 0.357893 | 0.259494 | 0.683544 | 117 | 50 | aligned_improvement |
| 30 | probe | 192.900753 | 134.963611 | -57.937142 | 0.466488 | 0.385269 | -0.081219 | 0.359504 | 0.280992 | 155 | 174 | energy_causality_misalignment |
| 36 | reference | 251.424456 | 215.628980 | -35.795477 | 0.488971 | 0.373225 | -0.115745 | 0.387755 | 0.268222 | 210 | 251 | energy_causality_misalignment |

Interpretation labels:
- aligned_improvement: energy decreases and F1 increases (threshold: |delta_E| > 1.0, |delta_F1| > 0.01).
- energy_causality_misalignment: energy decreases and F1 decreases.
- causal_improves_without_energy: energy increases and F1 increases.
- no_clear_change: otherwise or change too small/ambiguous.

## Conservative readout

- N values with aligned improvement: `12, 18, 24`.
- N values with energy/causality misalignment: `30, 36`.
- N values with causal-improves-without-energy: `none`.
- N values with no clear change: `none`.
- N=24 interpretation label: `aligned_improvement`.
- N=36 interpretation label: `energy_causality_misalignment`.
- Apparent transition between N=24 and N=36 (aligned at 24, misaligned at 36)? `yes`.
- Does this exploratory ladder support looking for geometry-dependent thermal mobility? exploratory hint only; one case seed and one optimizer seed per N is not enough to establish geometry dependence, but the table provides a working bracket worth examining with more seeds before any claim.
- Does this exploratory ladder justify more seeds or topologies? yes as an exploratory next step; broadening the seed pool and adding non-Minkowski topologies is appropriate before any promotion.

## Guardrails

- exploratory only.
- one family only (`minkowski`).
- one case seed only (`1959`).
- one optimizer seed only (`1987`).
- two budgets only (`short_10_10_4`, `medium_25_25_8`).
- no embeddability claim.
- no physical gamma claim.
- no theorem.
- no general N-scaling claim.
- N=24 and N=36 rows in the comparison table are read from existing exploratory CSVs and were not regenerated here.
