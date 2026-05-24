# N=24 Causal-Order Metric Probe

**Status:** exploratory only; not confirmation.

This probe compares direct causal-order preservation with final energy and interval RMSE for three gamma values. It is an accessibility/recoverability diagnostic only.

## Provenance

- Output directory: `explore/causal_order_metric_probe_n24`
- Runner: `explore/causal_order_metric_probe_n24/run_causal_order_metric_probe_n24.py`
- Command: `python3 explore/causal_order_metric_probe_n24/run_causal_order_metric_probe_n24.py`
- Generated at UTC: `2026-05-24T08:32:40+00:00`
- CSV: `explore/causal_order_metric_probe_n24/causal_order_metric_probe_n24.csv`
- Figure: `explore/causal_order_metric_probe_n24/causal_order_metric_probe_n24.svg`
- N: `24`
- T0: `100.0`
- warmup_limit: `10`
- anneal_limit: `10`
- max_data: `4`
- total runtime seconds: `0.039`
- Final causal-order diagnostics use `sim.rold`/`sim.xold`, not `rnew`/`xnew`.

## Answers

1. Best gamma by `final_energy`: `0.55`.
2. Best gamma by `interval_rmse`: `0.73`.
3. Best gamma by `causal_f1`: `0.55`.
4. `causal_f1` aligns with `final_energy` in this three-point probe.
5. Dominant causal-order failure mode by total count: `missing` (`missing=352`, `extra=85`).
6. Any `exact_match` true? `no`. Any `success_flag` true? `no`.
7. Conservative conclusion: accessibility/recoverability diagnostic only; no embeddability claim.

## Key Table

| gamma | final_energy | interval_rmse | causal_f1 | precision | recall | missing | extra | exact | success |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :---: | :---: |
| 0.55 | 130.446696 | 55649.304629 | 0.473029 | 0.686747 | 0.360759 | 101 | 26 | no | no |
| 0.73 | 145.005191 | 18585.571053 | 0.405405 | 0.703125 | 0.284810 | 113 | 19 | no | no |
| 0.82 | 142.986146 | 113319.222524 | 0.183486 | 0.333333 | 0.126582 | 138 | 40 | no | no |

## Guardrails

- Exploration only.
- One N=24 known-truth case.
- No physical gamma claim.
- No embeddability claim.
- No recovery claim unless exact-match/success criteria pass.
