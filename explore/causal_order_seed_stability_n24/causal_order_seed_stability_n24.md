# N=24 Causal-Order Seed Stability Probe

**Status:** exploratory only; not confirmation.

This probe tests whether the apparent alignment between low gamma, final energy, and causal F1 is stable across optimizer seeds. It is an accessibility/recoverability diagnostic only.

## Provenance

- Output directory: `explore/causal_order_seed_stability_n24`
- Runner: `explore/causal_order_seed_stability_n24/run_causal_order_seed_stability_n24.py`
- Command: `python3 explore/causal_order_seed_stability_n24/run_causal_order_seed_stability_n24.py`
- Generated at UTC: `2026-05-24T08:36:02+00:00`
- CSV: `explore/causal_order_seed_stability_n24/causal_order_seed_stability_n24.csv`
- Figure: `explore/causal_order_seed_stability_n24/causal_order_seed_stability_n24.svg`
- N: `24`
- T0: `100.0`
- warmup_limit: `10`
- anneal_limit: `10`
- max_data: `4`
- gammas: `0.50, 0.53, 0.55, 0.57, 0.60`
- optimizer seeds: `1959, 1962, 1987, 2001, 2026`
- runs completed: `25`
- total runtime seconds: `0.386`
- Final causal-order diagnostics use `sim.rold`/`sim.xold`, not `rnew`/`xnew`.

## Answers

1. Best mean `causal_f1`: gamma `0.50`.
2. Best median `causal_f1`: gamma `0.50`.
3. Best worst-case `causal_f1`: gamma `0.50`.
4. Pearson correlation between `final_energy` and `causal_f1` across seeds: `-0.132610`.
5. Missing relations dominate? `yes` (`missing=2968`, `extra=633`).
6. Any `exact_match` true? `no`. Any `success_flag` true? `no`.
7. Conservative conclusion only: no recovery claim, no embeddability claim.

## Aggregate Table

| gamma | runs | mean F1 | median F1 | worst F1 | best F1 | mean final E | mean RMSE | missing | extra | exact | success |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.50 | 5 | 0.397550 | 0.407240 | 0.354978 | 0.425532 | 114.345228 | 14268220.446929 | 562 | 129 | 0 | 0 |
| 0.53 | 5 | 0.360681 | 0.387387 | 0.282927 | 0.393162 | 117.603788 | 108649.561516 | 589 | 118 | 0 | 0 |
| 0.55 | 5 | 0.366021 | 0.377193 | 0.238095 | 0.473029 | 119.470835 | 434025.348386 | 581 | 132 | 0 | 0 |
| 0.57 | 5 | 0.321146 | 0.324324 | 0.243902 | 0.377193 | 120.146490 | 3410464.297193 | 612 | 136 | 0 | 0 |
| 0.60 | 5 | 0.308833 | 0.309735 | 0.283019 | 0.325581 | 121.980397 | 1316986.536152 | 624 | 118 | 0 | 0 |

## Guardrails

- Exploration only.
- One N=24 known-truth case.
- No physical gamma claim.
- No embeddability claim.
- No recovery claim.
