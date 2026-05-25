# Oracle checkpoint ceiling N=36

Post-run SORKIN-2 oracle diagnostic over the existing schedule/optimizer-seed trajectory matrix.

## Configuration

- Command: `python3 explore/oracle_checkpoint_ceiling_n36/run_oracle_checkpoint_ceiling_n36.py`
- Generated at UTC: `2026-05-25T09:12:10+00:00`
- Source CSV: `explore/schedule_seed_stability_n36/schedule_seed_stability_n36.csv`
- Output CSV: `explore/oracle_checkpoint_ceiling_n36/oracle_checkpoint_ceiling_n36.csv`
- Summary CSV: `explore/oracle_checkpoint_ceiling_n36/oracle_checkpoint_ceiling_n36_summary.csv`
- SVG: `explore/oracle_checkpoint_ceiling_n36/oracle_checkpoint_ceiling_n36.svg`
- This script only reads trajectory CSV rows; it does not run the annealer.
- Grouping key: `(optimizer_seed, schedule_label, cooling_factor, budget_label)`.
- Oracle checkpoint maximizes `causal_f1`, then `causal_recall`, then minimizes missing, then extra, then chooses the earliest block.

## Per group

| seed | schedule | gamma | budget | final F1 | best F1 | delta F1 | final recall | best recall | final missing | best missing | final E | best E | min E | final block | best block | minE block | best before final | best=minE |
| ---: | --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |
| 1959 | gamma_0p5 | 0.5 | medium_25_25_8 | 0.421456 | 0.548822 | 0.127366 | 0.3207 | 0.475219 | 233 | 180 | 225.887 | 318.668 | 225.887 | 8 | 6 | 8 | true | false |
| 1962 | gamma_0p5 | 0.5 | medium_25_25_8 | 0.609555 | 0.609555 | 0 | 0.539359 | 0.539359 | 158 | 158 | 220.347 | 220.347 | 220.347 | 8 | 8 | 8 | false | true |
| 1987 | gamma_0p5 | 0.5 | medium_25_25_8 | 0.373225 | 0.548701 | 0.175476 | 0.268222 | 0.492711 | 251 | 174 | 215.629 | 314.577 | 215.629 | 8 | 6 | 8 | true | false |
| 2001 | gamma_0p5 | 0.5 | medium_25_25_8 | 0.375 | 0.584229 | 0.209229 | 0.271137 | 0.475219 | 250 | 180 | 234.068 | 267.223 | 234.068 | 8 | 7 | 8 | true | false |
| 1959 | gamma_0p8 | 0.8 | medium_25_25_8 | 0.374245 | 0.461818 | 0.0875727 | 0.271137 | 0.370262 | 250 | 216 | 438.358 | 493.874 | 318.792 | 8 | 7 | 1 | true | false |
| 1962 | gamma_0p8 | 0.8 | medium_25_25_8 | 0.464419 | 0.464419 | 0 | 0.361516 | 0.361516 | 219 | 219 | 482.384 | 482.384 | 332.13 | 8 | 8 | 1 | false | false |
| 1987 | gamma_0p8 | 0.8 | medium_25_25_8 | 0.387226 | 0.429688 | 0.042462 | 0.282799 | 0.3207 | 246 | 233 | 447.107 | 548.521 | 322.55 | 8 | 5 | 1 | true | false |
| 2001 | gamma_0p8 | 0.8 | medium_25_25_8 | 0.350305 | 0.422481 | 0.0721751 | 0.250729 | 0.317784 | 257 | 234 | 488.843 | 509.627 | 281.665 | 8 | 7 | 1 | true | false |
| 1959 | gamma_0p9 | 0.9 | medium_25_25_8 | 0.297593 | 0.371257 | 0.0736645 | 0.198251 | 0.271137 | 275 | 250 | 509.643 | 318.792 | 318.792 | 8 | 1 | 1 | true | true |
| 1962 | gamma_0p9 | 0.9 | medium_25_25_8 | 0.268657 | 0.390057 | 0.121401 | 0.183673 | 0.297376 | 280 | 241 | 580.434 | 583.853 | 332.13 | 8 | 5 | 1 | true | false |
| 1987 | gamma_0p9 | 0.9 | medium_25_25_8 | 0.3762 | 0.447273 | 0.0710731 | 0.285714 | 0.358601 | 245 | 220 | 504.925 | 473.264 | 322.55 | 8 | 2 | 1 | true | false |
| 2001 | gamma_0p9 | 0.9 | medium_25_25_8 | 0.257511 | 0.385069 | 0.127558 | 0.174927 | 0.285714 | 283 | 245 | 521.377 | 538.599 | 281.665 | 8 | 7 | 1 | true | false |
| 1959 | gamma_0p95 | 0.95 | medium_25_25_8 | 0.235294 | 0.371257 | 0.135963 | 0.151603 | 0.271137 | 291 | 250 | 547.206 | 318.792 | 318.792 | 8 | 1 | 1 | true | true |
| 1962 | gamma_0p95 | 0.95 | medium_25_25_8 | 0.322851 | 0.403846 | 0.080995 | 0.22449 | 0.306122 | 266 | 238 | 592.986 | 632.192 | 332.13 | 8 | 3 | 1 | true | false |
| 1987 | gamma_0p95 | 0.95 | medium_25_25_8 | 0.354331 | 0.447273 | 0.092942 | 0.262391 | 0.358601 | 253 | 220 | 598.076 | 473.264 | 322.55 | 8 | 2 | 1 | true | false |
| 2001 | gamma_0p95 | 0.95 | medium_25_25_8 | 0.26087 | 0.375969 | 0.115099 | 0.174927 | 0.282799 | 283 | 246 | 563.107 | 281.665 | 281.665 | 8 | 1 | 1 | true | true |

## Schedule aggregates

| schedule | gamma | budget | groups | avg final F1 | avg oracle best F1 | avg delta F1 |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
| gamma_0p5 | 0.5 | medium_25_25_8 | 4 | 0.444809 | 0.572827 | 0.128018 |
| gamma_0p8 | 0.8 | medium_25_25_8 | 4 | 0.394049 | 0.444601 | 0.0505524 |
| gamma_0p9 | 0.9 | medium_25_25_8 | 4 | 0.29999 | 0.398414 | 0.0984241 |
| gamma_0p95 | 0.95 | medium_25_25_8 | 4 | 0.293336 | 0.399586 | 0.10625 |

## Global summary

- Groups: `16`.
- Average final causal F1: `0.358046`.
- Average oracle best checkpoint causal F1: `0.453857`.
- Average oracle gain: `0.0958111`.
- Median oracle gain: `0.0902574`.
- Max oracle gain: `0.209229`.
- Best checkpoint > final: `14` of `16`.
- Best checkpoint > final by 0.02: `14` of `16`.
- Best checkpoint before final: `14` of `16`.
- Best checkpoint matches minimum energy block: `4` of `16`.
- Best avg oracle schedule/budget/cooling: `gamma_0p5` / `medium_25_25_8` / `0.5` with avg best F1 `0.572827`.
- Best avg oracle gain schedule/budget/cooling: `gamma_0p5` / `medium_25_25_8` / `0.5` with avg gain `0.128018`.

## Conservative interpretation

This is an oracle diagnostic: it selects checkpoints using causal F1 against the known-truth target order.
It is not a causally deployable checkpoint-selection criterion for truth-free cases.
The oracle ceiling substantially improves the endpoint in most groups, so the N=36 limitation in this matrix is not only search access or budget: the annealer often visits better causal checkpoints that the endpoint selector discards.
Because the minimum-energy block coincides with the best-F1 block in only a minority of groups, historical energy should not be used alone as a recoverability selector in this diagnostic setting.

## Guardrails

This is a post-run diagnostic only, using benchmark cases with known truth.
It is not an embeddability claim, not a physical gamma claim, not an N-transition claim, and not proof of general annealer failure.
It is not a deployable criterion for truth-free cases yet.
