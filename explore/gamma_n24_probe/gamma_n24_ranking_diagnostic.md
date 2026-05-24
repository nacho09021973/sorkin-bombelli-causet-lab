# N=24 Gamma Ranking Diagnostic

**Status:** post-analysis of the completed exploratory N=24 CSV only. No annealer run is performed here.

This note is an accessibility diagnostic for the historical annealer. It makes no embeddability claim, no physical gamma claim, and no success claim.

## Provenance

- Input CSV: `explore/gamma_n24_probe/gamma_n24_probe.csv`
- Analysis script: `explore/gamma_n24_probe/build_gamma_n24_ranking_diagnostic.py`
- Derived CSV: `explore/gamma_n24_probe/gamma_n24_ranking_diagnostic.csv`
- Figure: `explore/gamma_n24_probe/gamma_n24_energy_vs_rmse.svg`
- Exact command: `python3 explore/gamma_n24_probe/build_gamma_n24_ranking_diagnostic.py`

## Metric Equivalence and Limits

- `truth_energy` is zero for every row, so `final_energy` and `energy_gap` are equivalent in this run.
- Every `success_flag` value is false, so the success ranking cannot distinguish gamma values and cannot establish recovery.

## Answers

- Best gamma by `final_energy`: `0.562373` with final energy `130.465017`.
- Best gamma by `interval_rmse`: `0.848624` with interval RMSE `113319.222500`.
- Do the rankings agree? no.
- Apparent optimum: metric-dependent. The final-energy optimum sits in a broad 10% final-energy band, while interval RMSE favors a different gamma.
- Next exploratory priority: stabilize the diagnostic metric first, before refining gamma or changing T0.

## Ranking Table

| gamma | final_energy | final rank | energy_gap | gap rank | interval_rmse | RMSE rank | success | within 10% final |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | :---: | :---: |
| 0.500000 | 131.822227 | 3 | 131.822227 | 3 | 71268273.480000 | 9 | no | yes |
| 0.530270 | 131.092389 | 2 | 131.092389 | 2 | 340916.400700 | 3 | no | yes |
| 0.562373 | 130.465017 | 1 | 130.465017 | 1 | 9137960.630000 | 8 | no | yes |
| 0.596419 | 135.645205 | 4 | 135.645205 | 4 | 6106828.456000 | 6 | no | yes |
| 0.632527 | 135.645205 | 4 | 135.645205 | 4 | 6106828.456000 | 6 | no | yes |
| 0.670820 | 135.645205 | 4 | 135.645205 | 4 | 6106828.456000 | 6 | no | yes |
| 0.711432 | 151.237389 | 9 | 151.237389 | 9 | 344640.349000 | 4 | no | no |
| 0.754503 | 139.465256 | 6 | 139.465256 | 6 | 268366.496900 | 2 | no | yes |
| 0.800181 | 139.321485 | 5 | 139.321485 | 5 | 7210040.598000 | 7 | no | yes |
| 0.848624 | 142.986146 | 7 | 142.986146 | 7 | 113319.222500 | 1 | no | yes |
| 0.900000 | 143.563028 | 8 | 143.563028 | 8 | 729055.731900 | 5 | no | no |

## Conservative Readout

- `final_energy` and `energy_gap` select the same gamma because the ground-truth energy is zero in every row.
- `interval_rmse` selects a different gamma, so the apparent gamma optimum is not robust across diagnostics.
- Since all success flags are false, this post-analysis does not establish recovery.
