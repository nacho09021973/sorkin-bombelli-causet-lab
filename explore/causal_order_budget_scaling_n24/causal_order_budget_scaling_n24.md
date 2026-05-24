# N=24 Causal-Order Budget Scaling Probe

**Status:** exploratory only; not confirmation.

This probe tests whether the dominant missing-relation failure mode improves when the historical annealer budget is increased at fixed gamma 0.50.

## Provenance

- Output directory: `explore/causal_order_budget_scaling_n24`
- Runner: `explore/causal_order_budget_scaling_n24/run_causal_order_budget_scaling_n24.py`
- Command: `python3 explore/causal_order_budget_scaling_n24/run_causal_order_budget_scaling_n24.py`
- Generated at UTC: `2026-05-24T09:02:44+00:00`
- CSV: `explore/causal_order_budget_scaling_n24/causal_order_budget_scaling_n24.csv`
- Figure: `explore/causal_order_budget_scaling_n24/causal_order_budget_scaling_n24.svg`
- N: `24`
- T0: `100.0`
- gamma: `0.5`
- optimizer seeds: `1959, 1962, 1987, 2001, 2026`
- runs completed: `15`
- new rows completed in this resume: `5`
- completed short budget seeds: `5`
- completed medium budget seeds: `5`
- completed long budget seeds: `5`
- long budget status: `deferred / not attempted in this resume`
- planned active resume rows: `5`
- total runtime seconds: `533.222`
- stopped early: `false`
- stop reason: `completed active medium-budget resume rows; long budget deferred`
- per-run timeout seconds: `900.0`
- timeout note: timeout was raised from 60 seconds to 180 seconds for this medium-budget resume.
- Final causal-order diagnostics use `sim.rold`/`sim.xold`, not `rnew`/`xnew`.

## Answers

1. Does medium mean `causal_f1` remain above short mean `causal_f1`? `yes`.
2. Does medium mean recall remain above short mean recall? `yes`.
3. Are missing relations still reduced at medium budget? `yes`.
4. Does extra relation count increase significantly? `no`.
5. Pearson correlation between `final_energy` and `causal_f1`: `-0.958042`.
6. Any `exact_match` true? `yes`. Any `success_flag` true? `yes`.
7. Is the medium improvement stable across seeds or driven by one outlier? `not one outlier`.
8. Conservative conclusion only: accessibility/recoverability diagnostic, no embeddability claim.

## Aggregate Table

| budget | runs | mean F1 | mean recall | mean precision | mean final E | missing | extra | mean missing | mean extra | exact | success | mean runtime s |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| short_10_10_4 | 5 | 0.397550 | 0.288608 | 0.640612 | 114.345228 | 562 | 129 | 112.400 | 25.800 | 0 | 0 | 0.019 |
| medium_25_25_8 | 5 | 0.686196 | 0.627848 | 0.762365 | 84.710796 | 294 | 155 | 58.800 | 31.000 | 0 | 0 | 36.828 |
| long_50_50_16 | 5 | 0.993610 | 0.992405 | 0.994935 | 0.137267 | 6 | 4 | 1.200 | 0.800 | 2 | 5 | 106.640 |

## Guardrails

- Exploration only.
- One N=24 known-truth case.
- No physical gamma claim.
- No embeddability claim.
- No recovery claim.
