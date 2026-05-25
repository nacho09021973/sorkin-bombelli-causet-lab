# Schedule seed stability N=36

Exploratory SORKIN-2 diagnostic testing whether the schedule signal persists across optimizer seeds for one fixed known-truth causal set.

## Configuration

- Command: `python3 explore/schedule_seed_stability_n36/run_schedule_seed_stability_n36.py`
- Generated at UTC: `2026-05-25T08:17:13+00:00`
- Output directory: `explore/schedule_seed_stability_n36`
- CSV: `explore/schedule_seed_stability_n36/schedule_seed_stability_n36.csv`
- Summary CSV: `explore/schedule_seed_stability_n36/schedule_seed_stability_n36_summary.csv`
- SVG: `explore/schedule_seed_stability_n36/schedule_seed_stability_n36.svg`
- family: `minkowski`
- N: `36`
- d_spacetime: `2`
- target spatial dim: `1`
- case_seed: `1959`
- optimizer_seeds: `1959, 1962, 1987, 2001`
- T0: `100.0`
- budget: `medium_25_25_8`
- warmup_limit: `25`
- anneal_limit: `25`
- max_data: `8`
- Schedules vary only native `cooling_factor`: `0.5`, `0.8`, `0.9`, `0.95`.
- block_callback reads `sim.rold`/`sim.xold` only; historical energy, move set, acceptance rule, and schedule mechanism are unchanged.

## Summary matrix

| seed | schedule | gamma | final F1 | best F1 seen | block best F1 | final recall | final missing | min Eave | block min E | Edown F1down | Edown recalldown | Edown missingup |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1959 | gamma_0p5 | 0.5 | 0.421456 | 0.548822 | 6 | 0.3207 | 233 | 225.887 | 8 | 2 | 2 | 2 |
| 1959 | gamma_0p8 | 0.8 | 0.374245 | 0.461818 | 7 | 0.271137 | 250 | 318.792 | 1 | 2 | 2 | 2 |
| 1959 | gamma_0p9 | 0.9 | 0.297593 | 0.371257 | 1 | 0.198251 | 275 | 318.792 | 1 | 1 | 2 | 2 |
| 1959 | gamma_0p95 | 0.95 | 0.235294 | 0.371257 | 1 | 0.151603 | 291 | 318.792 | 1 | 2 | 2 | 2 |
| 1962 | gamma_0p5 | 0.5 | 0.609555 | 0.609555 | 8 | 0.539359 | 158 | 220.347 | 8 | 2 | 2 | 2 |
| 1962 | gamma_0p8 | 0.8 | 0.464419 | 0.464419 | 8 | 0.361516 | 219 | 332.13 | 1 | 2 | 2 | 2 |
| 1962 | gamma_0p9 | 0.9 | 0.268657 | 0.390057 | 5 | 0.183673 | 280 | 332.13 | 1 | 1 | 1 | 1 |
| 1962 | gamma_0p95 | 0.95 | 0.322851 | 0.403846 | 3 | 0.22449 | 266 | 332.13 | 1 | 1 | 1 | 1 |
| 1987 | gamma_0p5 | 0.5 | 0.373225 | 0.548701 | 6 | 0.268222 | 251 | 215.629 | 8 | 3 | 3 | 3 |
| 1987 | gamma_0p8 | 0.8 | 0.387226 | 0.429688 | 5 | 0.282799 | 246 | 322.55 | 1 | 2 | 3 | 3 |
| 1987 | gamma_0p9 | 0.9 | 0.3762 | 0.447273 | 2 | 0.285714 | 245 | 322.55 | 1 | 4 | 4 | 4 |
| 1987 | gamma_0p95 | 0.95 | 0.354331 | 0.447273 | 2 | 0.262391 | 253 | 322.55 | 1 | 2 | 2 | 2 |
| 2001 | gamma_0p5 | 0.5 | 0.375 | 0.584229 | 7 | 0.271137 | 250 | 234.068 | 8 | 2 | 3 | 3 |
| 2001 | gamma_0p8 | 0.8 | 0.350305 | 0.422481 | 7 | 0.250729 | 257 | 281.665 | 1 | 1 | 1 | 1 |
| 2001 | gamma_0p9 | 0.9 | 0.257511 | 0.385069 | 7 | 0.174927 | 283 | 281.665 | 1 | 2 | 2 | 2 |
| 2001 | gamma_0p95 | 0.95 | 0.26087 | 0.375969 | 1 | 0.174927 | 283 | 281.665 | 1 | 2 | 2 | 2 |

## Per-seed winners

| seed | final F1 winner | best F1 winner | lowest tension winner | min energy winner | min energy = final F1? | min energy = best F1? |
| ---: | --- | --- | --- | --- | --- | --- |
| 1959 | gamma_0p5 | gamma_0p5 | gamma_0p9 | gamma_0p5 | true | true |
| 1962 | gamma_0p5 | gamma_0p5 | gamma_0p9 | gamma_0p5 | true | true |
| 1987 | gamma_0p8 | gamma_0p5 | gamma_0p8 | gamma_0p5 | false | true |
| 2001 | gamma_0p5 | gamma_0p5 | gamma_0p8 | gamma_0p5 | true | true |

## Aggregate winner counts

| schedule | final F1 wins | best F1 wins | lowest tension wins |
| --- | ---: | ---: | ---: |
| gamma_0p5 | 3 | 4 | 0 |
| gamma_0p8 | 1 | 0 | 2 |
| gamma_0p9 | 0 | 0 | 2 |
| gamma_0p95 | 0 | 0 | 0 |

## Readout

- Minimum energy coincides with the final-F1 winner in `3` of `4` optimizer seeds.
- Minimum energy coincides with the best-trajectory-F1 winner in `4` of `4` optimizer seeds.
- At least one slower gamma reduces `energy_down_f1_down` relative to gamma 0.5 in `4` of `4` optimizer seeds.
One schedule dominates one causal-F1 readout across all optimizer seeds in this small matrix; that would support cooling as a robust algorithmic factor for this fixed causal set.
Slower cooling reduces at least one tension count in part of the matrix, so thermal mobility remains a plausible diagnostic lever.

## Guardrails

This is exploratory: one N, one family, one case seed, four optimizer seeds, one budget, and the unchanged historical objective function.
It is not evidence of embeddability or non-embeddability, does not identify a physical gamma, does not establish an N transition, and does not by itself justify changing the objective function.
