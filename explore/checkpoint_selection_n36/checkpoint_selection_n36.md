# Checkpoint selection N=36

Post-run SORKIN-2 diagnostic over the existing schedule/optimizer-seed trajectory matrix.

## Configuration

- Command: `python3 explore/checkpoint_selection_n36/run_checkpoint_selection_n36.py`
- Generated at UTC: `2026-05-25T08:27:52+00:00`
- Source CSV: `explore/schedule_seed_stability_n36/schedule_seed_stability_n36.csv`
- Output CSV: `explore/checkpoint_selection_n36/checkpoint_selection_n36.csv`
- Summary CSV: `explore/checkpoint_selection_n36/checkpoint_selection_n36_summary.csv`
- SVG: `explore/checkpoint_selection_n36/checkpoint_selection_n36.svg`
- This script only reads trajectory CSV rows; it does not run the annealer.
- Best checkpoint maximizes `causal_f1`, then `causal_recall`, then minimizes missing, then extra, then chooses the earliest block.

## Per group

| schedule | gamma | seed | final block | final F1 | best block | best F1 | delta F1 | final recall | best recall | final missing | best missing | min E block | best=minE |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| gamma_0p5 | 0.5 | 1959 | 8 | 0.421456 | 6 | 0.548822 | 0.127366 | 0.3207 | 0.475219 | 233 | 180 | 8 | false |
| gamma_0p5 | 0.5 | 1962 | 8 | 0.609555 | 8 | 0.609555 | 0 | 0.539359 | 0.539359 | 158 | 158 | 8 | true |
| gamma_0p5 | 0.5 | 1987 | 8 | 0.373225 | 6 | 0.548701 | 0.175476 | 0.268222 | 0.492711 | 251 | 174 | 8 | false |
| gamma_0p5 | 0.5 | 2001 | 8 | 0.375 | 7 | 0.584229 | 0.209229 | 0.271137 | 0.475219 | 250 | 180 | 8 | false |
| gamma_0p8 | 0.8 | 1959 | 8 | 0.374245 | 7 | 0.461818 | 0.0875727 | 0.271137 | 0.370262 | 250 | 216 | 1 | false |
| gamma_0p8 | 0.8 | 1962 | 8 | 0.464419 | 8 | 0.464419 | 0 | 0.361516 | 0.361516 | 219 | 219 | 1 | false |
| gamma_0p8 | 0.8 | 1987 | 8 | 0.387226 | 5 | 0.429688 | 0.042462 | 0.282799 | 0.3207 | 246 | 233 | 1 | false |
| gamma_0p8 | 0.8 | 2001 | 8 | 0.350305 | 7 | 0.422481 | 0.0721751 | 0.250729 | 0.317784 | 257 | 234 | 1 | false |
| gamma_0p9 | 0.9 | 1959 | 8 | 0.297593 | 1 | 0.371257 | 0.0736645 | 0.198251 | 0.271137 | 275 | 250 | 1 | true |
| gamma_0p9 | 0.9 | 1962 | 8 | 0.268657 | 5 | 0.390057 | 0.121401 | 0.183673 | 0.297376 | 280 | 241 | 1 | false |
| gamma_0p9 | 0.9 | 1987 | 8 | 0.3762 | 2 | 0.447273 | 0.0710731 | 0.285714 | 0.358601 | 245 | 220 | 1 | false |
| gamma_0p9 | 0.9 | 2001 | 8 | 0.257511 | 7 | 0.385069 | 0.127558 | 0.174927 | 0.285714 | 283 | 245 | 1 | false |
| gamma_0p95 | 0.95 | 1959 | 8 | 0.235294 | 1 | 0.371257 | 0.135963 | 0.151603 | 0.271137 | 291 | 250 | 1 | true |
| gamma_0p95 | 0.95 | 1962 | 8 | 0.322851 | 3 | 0.403846 | 0.080995 | 0.22449 | 0.306122 | 266 | 238 | 1 | false |
| gamma_0p95 | 0.95 | 1987 | 8 | 0.354331 | 2 | 0.447273 | 0.092942 | 0.262391 | 0.358601 | 253 | 220 | 1 | false |
| gamma_0p95 | 0.95 | 2001 | 8 | 0.26087 | 1 | 0.375969 | 0.115099 | 0.174927 | 0.282799 | 283 | 246 | 1 | true |

## Summary by schedule

| schedule | gamma | groups | avg final F1 | avg best F1 | avg delta F1 | better | better >0.02 | before final | best is final | best=minE | avg final recall | avg best recall | avg final missing | avg best missing |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| gamma_0p5 | 0.5 | 4 | 0.444809 | 0.572827 | 0.128018 | 3 | 3 | 3 | 1 | 1 | 0.349854 | 0.495627 | 223 | 173 |
| gamma_0p8 | 0.8 | 4 | 0.394049 | 0.444601 | 0.0505524 | 3 | 3 | 3 | 1 | 0 | 0.291545 | 0.342566 | 243 | 225.5 |
| gamma_0p9 | 0.9 | 4 | 0.29999 | 0.398414 | 0.0984241 | 4 | 4 | 4 | 0 | 1 | 0.210641 | 0.303207 | 270.75 | 239 |
| gamma_0p95 | 0.95 | 4 | 0.293336 | 0.399586 | 0.10625 | 4 | 4 | 4 | 0 | 2 | 0.203353 | 0.304665 | 273.25 | 238.5 |

## Readout

- Total groups: `16`.
- Best checkpoint better than final: `14` of `16`.
- Best checkpoint better than final by more than 0.02 F1: `14` of `16`.
- Best checkpoint before final: `14` of `16`.
- Best avg checkpoint F1 schedule: `gamma_0p5` with avg best F1 `0.572827`.
- Best-F1 block equals minimum-energy block in `4` of `16` groups.
Checkpoint selection appears more promising here than simply lowering gamma more slowly: the best causal checkpoint often occurs before the final block, while slower schedules did not robustly dominate best causal F1 in the seed-stability matrix.

## Guardrails

This is a post-run diagnostic only, using benchmark cases with known truth.
It is not an embeddability claim, not a physical gamma claim, not an N-transition claim, and not proof of general annealer failure.
It is also not a deployable criterion for truth-free cases yet, because selecting by causal F1 uses known-truth labels.
