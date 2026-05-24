# N=36 Causal-Order Budget Scaling Pilot

**Status:** exploratory only; not confirmation.

Does increasing annealer budget improve N=36 causal-order recoverability at fixed T0 and gamma?

## Provenance

- Output directory: `explore/causal_order_budget_scaling_n36`
- Runner: `explore/causal_order_budget_scaling_n36/run_causal_order_budget_scaling_n36.py`
- Command: `python3 explore/causal_order_budget_scaling_n36/run_causal_order_budget_scaling_n36.py`
- Generated at UTC: `2026-05-24T10:27:56+00:00`
- CSV: `explore/causal_order_budget_scaling_n36/causal_order_budget_scaling_n36.csv`
- Figure: `explore/causal_order_budget_scaling_n36/causal_order_budget_scaling_n36.svg`
- Log: `explore/causal_order_budget_scaling_n36/causal_order_budget_scaling_n36.log`
- N: `36`
- T0: `100.0`
- gamma: `0.5`
- case seed: `1959`
- optimizer seed(s): `1987`
- budgets attempted: `short_10_10_4, medium_25_25_8`
- completed rows: `2`
- completed budgets: `short_10_10_4, medium_25_25_8`
- deferred rows: `long_50_50_16`
- timeout rows: `none`
- timeout policy: `900.0` seconds per budget row
- long-budget deferral policy: defer long if medium runtime exceeds `300.0` seconds
- total wall runtime seconds: `688.792`
- Final causal-order diagnostics use `sim.rold`/`sim.xold`, not `rnew`/`xnew`.

## Aggregate Table

| budget | runs | mean F1 | mean recall | mean precision | mean final E | mean RMSE | missing | extra | mean missing | mean extra | exact | success | mean runtime s |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| short_10_10_4 | 1 | 0.488971 | 0.387755 | 0.661692 | 251.424457 | 413612.895982 | 210 | 68 | 210.000 | 68.000 | 0 | 0 | 0.059 |
| medium_25_25_8 | 1 | 0.373225 | 0.268222 | 0.613333 | 215.628980 | 664950.561916 | 251 | 58 | 251.000 | 58.000 | 0 | 0 | 688.715 |

## Conservative Answers

- Does F1 improve with budget? `no`.
- Does recall improve with budget? `no`.
- Are missing relations reduced? `no`.
- Do extra relations increase? `no`.
- Any exact match? `no`.
- Any success flag? `no`.
- Is improvement monotonic over attempted budgets? `no`.
- Does this look like a budget/accessibility effect? `no`.

## Row Status

| budget | status | runtime s | note |
| --- | --- | ---: | --- |
| short_10_10_4 | completed | 0.059 |  |
| medium_25_25_8 | completed | 688.715 |  |
| long_50_50_16 | deferred | 0.000 | deferred because medium runtime exceeded long-budget gate |

## Guardrails

- One N=36 case only.
- One/few seeds only.
- Exploratory only.
- No embeddability claim.
- No physical gamma claim.
- No theorem.
- No general N=36 conclusion.
- Further seed/topology sweeps required before promotion.
