# N=36 Optimizer-Seed Budget/Misalignment Probe

exploratory only; not confirmation.

At N=36, fixed family/case/schedule/budgets, does the short -> medium energy/causal-metric misalignment observed at optimizer seed 1987 persist across optimizer seeds, or was it seed-specific?

## Provenance

- Output directory: `explore/causal_order_budget_scaling_n36_seedsweep`
- Runner: `explore/causal_order_budget_scaling_n36_seedsweep/run_causal_order_budget_scaling_n36_seedsweep.py`
- Command: `python3 explore/causal_order_budget_scaling_n36_seedsweep/run_causal_order_budget_scaling_n36_seedsweep.py`
- Generated at UTC: `2026-05-24T17:57:33+00:00`
- CSV: `explore/causal_order_budget_scaling_n36_seedsweep/causal_order_budget_scaling_n36_seedsweep.csv`
- Figure: `explore/causal_order_budget_scaling_n36_seedsweep/causal_order_budget_scaling_n36_seedsweep.svg`
- Log: `explore/causal_order_budget_scaling_n36_seedsweep/causal_order_budget_scaling_n36_seedsweep.log`
- N: `36`
- family: `minkowski`
- d_spacetime: `2`
- case seed: `1959`
- optimizer seeds: `1959, 1962, 1987, 2001, 2026`
- T0: `100.0`
- gamma: `0.5`
- h: `1`
- budgets: `short_10_10_4, medium_25_25_8`
- reused rows (from `explore/causal_order_budget_scaling_n36/causal_order_budget_scaling_n36.csv`): `none`
- newly run optimizer seeds: `none`
- skipped rows: `none`
- per-run timeout seconds: `900.0`
- total wall runtime seconds: `0.001`
- Final causal-order diagnostics use `sim.rold`/`sim.xold`, not `rnew`/`xnew`.

## Result table

| seed | budget | status | final E | interval RMSE | F1 | recall | precision | missing | extra | exact | success | runtime s | note |
| ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: | --- |
| 1959 | short_10_10_4 | completed | 281.044498 | 583447.3943 | 0.429119 | 0.326531 | 0.625698 | 231 | 67 | no | no | 0.044 |  |
| 1959 | medium_25_25_8 | completed | 225.887112 | 7224.7552 | 0.421456 | 0.320700 | 0.614525 | 233 | 69 | no | no | 275.778 |  |
| 1962 | short_10_10_4 | completed | 274.458006 | 2669540.8180 | 0.404669 | 0.303207 | 0.608187 | 239 | 67 | no | no | 0.214 |  |
| 1962 | medium_25_25_8 | completed | 220.347317 | 7590.6575 | 0.609555 | 0.539359 | 0.700758 | 158 | 79 | no | no | 328.023 |  |
| 1987 | short_10_10_4 | completed | 251.424456 | 413612.8960 | 0.488971 | 0.387755 | 0.661692 | 210 | 68 | no | no | 0.059 | reused from explore/causal_order_budget_scaling_n36/causal_order_budget_scaling_n36.csv |
| 1987 | medium_25_25_8 | completed | 215.628980 | 664950.5619 | 0.373225 | 0.268222 | 0.613333 | 251 | 58 | no | no | 688.715 | reused from explore/causal_order_budget_scaling_n36/causal_order_budget_scaling_n36.csv |
| 2001 | short_10_10_4 | completed | 216.837536 | 113689.0693 | 0.385686 | 0.282799 | 0.606250 | 246 | 63 | no | no | 0.158 |  |
| 2001 | medium_25_25_8 | completed | 234.068370 | 335124.0292 | 0.375000 | 0.271137 | 0.607843 | 250 | 60 | no | no | 622.363 |  |
| 2026 | short_10_10_4 | completed | 253.549550 | 664918.7426 | 0.451376 | 0.358601 | 0.608911 | 220 | 79 | no | no | 0.064 |  |
| 2026 | medium_25_25_8 | timeout | NA | NA | NA | NA | NA | 0 | 0 | no | no | 900.005 | timeout after 900.0s |

## Seed-level comparison (short -> medium)

| seed | E_short | E_medium | delta_E | F1_short | F1_medium | delta_F1 | recall_short | recall_medium | missing_short | missing_medium | label |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 1959 | 281.044498 | 225.887112 | -55.157386 | 0.429119 | 0.421456 | -0.007663 | 0.326531 | 0.320700 | 231 | 233 | no_clear_change |
| 1962 | 274.458006 | 220.347317 | -54.110689 | 0.404669 | 0.609555 | 0.204886 | 0.303207 | 0.539359 | 239 | 158 | aligned_improvement |
| 1987 | 251.424456 | 215.628980 | -35.795477 | 0.488971 | 0.373225 | -0.115745 | 0.387755 | 0.268222 | 210 | 251 | energy_causality_misalignment |
| 2001 | 216.837536 | 234.068370 | 17.230834 | 0.385686 | 0.375000 | -0.010686 | 0.282799 | 0.271137 | 246 | 250 | no_clear_change |

Interpretation labels (thresholds: |delta_E| > 1.0, |delta_F1| > 0.02):
- aligned_improvement: delta_E < -1.0 and delta_F1 > +0.02.
- energy_causality_misalignment: delta_E < -1.0 and delta_F1 < -0.02.
- causal_improves_without_energy: delta_E > +1.0 and delta_F1 > +0.02.
- no_clear_change: otherwise.

## Decision rule

- seeds with energy_causality_misalignment: `1/4` (1987).
- seeds with aligned_improvement: `1/4` (1962).
- seeds with causal_improves_without_energy: `0/4` (none).
- seeds with no_clear_change: `2/4` (1959, 2001).
- rule: >=4/5 misalignment -> persists; <=2/5 -> seed-specific; otherwise mixed/inconclusive.
- outcome: `misalignment_appears_seed_specific`.

## Conservative readout

- Does the N=36 short->medium misalignment persist across optimizer seeds at this case seed/family/schedule? `no`.
- Is the N=30 seed-specific result repeated or overturned at N=36? `repeated (N=36 misalignment also looks seed-specific, mirroring the N=30 finding) -- weakens any N-driven transition story`.
- Does this justify a topology / case-seed follow-up at N=36? `low priority for topology probe; the multi-seed pattern aligns with N=30 and suggests basin-selection drives the seed-1987 misalignment rather than topology`.
- This probe does not establish: a transition in N, a physical gamma claim, an embeddability claim, a theorem, or a general N-scaling statement.

## Guardrails

- one N only (`36`).
- one case seed only (`1959`).
- one family only (`minkowski`).
- one schedule only (`T0=100`, `gamma=0.5`, `h=1`).
- two budgets only (`short_10_10_4`, `medium_25_25_8`).
- exploratory only.
- no embeddability claim.
- no physical gamma claim.
- no theorem.
- no general N-scaling claim.
- optimizer seed 1987 rows are reused from `explore/causal_order_budget_scaling_n36/causal_order_budget_scaling_n36.csv` and were not regenerated here.
