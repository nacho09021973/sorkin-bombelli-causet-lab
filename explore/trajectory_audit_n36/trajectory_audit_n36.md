# Trajectory audit N=36

Exploratory SORKIN-2 diagnostic for one known-truth case and one optimizer seed.

## Configuration

- Command: `python3 explore/trajectory_audit_n36/run_trajectory_audit_n36.py`
- Generated at UTC: `2026-05-25T07:09:09+00:00`
- Output directory: `explore/trajectory_audit_n36`
- CSV: `explore/trajectory_audit_n36/trajectory_audit_n36.csv`
- SVG: `explore/trajectory_audit_n36/trajectory_audit_n36.svg`
- family: `minkowski`
- N: `36`
- d_spacetime: `2`
- target spatial dim: `1`
- case_seed: `1959`
- optimizer_seed: `1987`
- T0: `100.0`
- gamma: `0.5`
- h: `1`
- backend: `cpu`
- This probe uses `ConesSimulator.block_callback` in read-only mode.
- Checkpoint causal metrics use `sim.rold`/`sim.xold`, the accepted state documented by `validation_suite.verify_recovery` and reported by `ConesSimulator.writeout`.
- The probe does not modify the historical annealer, selection rule, acceptance rule, energy function, temperature schedule, coordinates, or causal order.

## Blocks

| budget | block | temp | energy_eave | F1 | recall | precision | missing | extra | induced | correct | exact |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| medium_25_25_8 | 1 | 100 | 322.55 | 0.292035 | 0.19242 | 0.605505 | 277 | 43 | 109 | 66 | false |
| medium_25_25_8 | 2 | 50 | 444.849 | 0.341176 | 0.253644 | 0.520958 | 256 | 80 | 167 | 87 | false |
| medium_25_25_8 | 3 | 25 | 548.037 | 0.396825 | 0.291545 | 0.621118 | 243 | 61 | 161 | 100 | false |
| medium_25_25_8 | 4 | 12.5 | 477.652 | 0.285714 | 0.201166 | 0.492857 | 274 | 71 | 140 | 69 | false |
| medium_25_25_8 | 5 | 6.25 | 389.694 | 0.476703 | 0.387755 | 0.618605 | 210 | 82 | 215 | 133 | false |
| medium_25_25_8 | 6 | 3.125 | 314.577 | 0.548701 | 0.492711 | 0.619048 | 174 | 104 | 273 | 169 | false |
| medium_25_25_8 | 7 | 1.5625 | 263.518 | 0.53833 | 0.460641 | 0.647541 | 185 | 86 | 244 | 158 | false |
| medium_25_25_8 | 8 | 0.78125 | 215.629 | 0.373225 | 0.268222 | 0.613333 | 251 | 58 | 150 | 92 | false |
| short_10_10_4 | 1 | 100 | 110.557 | 0.202586 | 0.137026 | 0.38843 | 296 | 74 | 121 | 47 | false |
| short_10_10_4 | 2 | 50 | 172.755 | 0.408644 | 0.303207 | 0.626506 | 239 | 62 | 166 | 104 | false |
| short_10_10_4 | 3 | 25 | 214.807 | 0.239278 | 0.154519 | 0.53 | 290 | 47 | 100 | 53 | false |
| short_10_10_4 | 4 | 12.5 | 251.424 | 0.488971 | 0.387755 | 0.661692 | 210 | 68 | 201 | 133 | false |

## Consecutive deltas

| budget | from | to | delta_energy_eave | delta_causal_f1 | delta_recall | delta_missing | delta_extra | energy_down_f1_down | energy_down_recall_down | energy_down_missing_up | energy_down_extra_up |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | --- | --- |
| short_10_10_4 | 1 | 2 | 62.1973 | 0.206058 | 0.166181 | -57 | -12 | false | false | false | false |
| short_10_10_4 | 2 | 3 | 42.0527 | -0.169367 | -0.148688 | 51 | -15 | false | false | false | false |
| short_10_10_4 | 3 | 4 | 36.6173 | 0.249693 | 0.233236 | -80 | 21 | false | false | false | false |
| medium_25_25_8 | 1 | 2 | 122.299 | 0.0491411 | 0.0612245 | -21 | 37 | false | false | false | false |
| medium_25_25_8 | 2 | 3 | 103.188 | 0.0556489 | 0.0379009 | -13 | -19 | false | false | false | false |
| medium_25_25_8 | 3 | 4 | -70.3848 | -0.111111 | -0.090379 | 31 | 10 | true | true | true | true |
| medium_25_25_8 | 4 | 5 | -87.9579 | 0.190988 | 0.186589 | -64 | 11 | false | false | false | true |
| medium_25_25_8 | 5 | 6 | -75.1175 | 0.0719988 | 0.104956 | -36 | 22 | false | false | false | true |
| medium_25_25_8 | 6 | 7 | -51.0589 | -0.0103708 | -0.03207 | 11 | -18 | true | true | true | false |
| medium_25_25_8 | 7 | 8 | -47.8889 | -0.165105 | -0.19242 | 66 | -28 | true | true | true | false |

## Pattern counts

- energy_down_f1_down: `3`
- energy_down_recall_down: `3`
- energy_down_missing_up: `3`
- energy_down_extra_up: `3`

## Conservative interpretation

At least one consecutive block shows within-trajectory energy/causal tension: the historical energy decreases while one audited causal metric worsens.
This is a diagnostic for one N=36 Minkowski case, one case seed, one optimizer seed, and two short budgets.
It is not evidence of embeddability or non-embeddability, does not identify a physical gamma, does not establish a transition in N, and does not by itself justify changing the objective function.
