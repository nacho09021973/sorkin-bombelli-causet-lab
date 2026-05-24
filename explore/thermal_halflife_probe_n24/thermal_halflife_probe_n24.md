# N=24 Thermal Half-Life Probe

## Status

Exploratory only; not confirmation.

This probe scans h where gamma = 2^(-1/h), so h is the number of annealing steps required to halve the temperature.

## Run Configuration

- Command: `python3 explore/thermal_halflife_probe_n24/run_thermal_halflife_probe_n24.py`
- Generated UTC: `2026-05-24T11:37:46+00:00`
- N: `24`
- T0: `100.0`
- case seed: `1959`
- optimizer seed: `1987`
- budget label: `medium_25_25_8`
- budget: warmup `25`, anneal `25`, max data `8`

## h/gamma Table

| h | gamma |
|---:|---:|
| 1 | 0.5000000000 |
| 2 | 0.7071067812 |
| 3 | 0.7937005260 |
| 4 | 0.8408964153 |
| 5 | 0.8705505633 |
| 6 | 0.8908987181 |
| 8 | 0.9170040432 |
| 10 | 0.9330329915 |
| 12 | 0.9438743127 |

## Results

| h | gamma | final_energy | interval_rmse | precision | recall | F1 | missing | extra | exact_match | success | runtime_s |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|:---:|---:|
| 1 | 0.5000000000 | 80.9892 | 61.1911 | 0.744828 | 0.683544 | 0.712871 | 50 | 37 | false | false | 24.597 |
| 2 | 0.7071067812 | 181.051 | 7.69883e+06 | 0.698795 | 0.367089 | 0.481328 | 100 | 25 | false | false | 0.103 |
| 3 | 0.7937005260 | 214.557 | 587571 | 0.642857 | 0.227848 | 0.336449 | 122 | 20 | false | false | 0.056 |
| 4 | 0.8408964153 | 241.943 | 183430 | 0.541667 | 0.246835 | 0.339130 | 119 | 33 | false | false | 0.050 |
| 5 | 0.8705505633 | 251.974 | 4.92638e+07 | 0.654545 | 0.227848 | 0.338028 | 122 | 19 | false | false | 0.047 |
| 6 | 0.8908987181 | 258.914 | 353401 | 0.546667 | 0.259494 | 0.351931 | 117 | 34 | false | false | 0.049 |
| 8 | 0.9170040432 | 250.698 | 259112 | 0.605263 | 0.145570 | 0.234694 | 135 | 15 | false | false | 0.053 |
| 10 | 0.9330329915 | 271.393 | 18992 | 0.759259 | 0.259494 | 0.386792 | 117 | 13 | false | false | 0.047 |
| 12 | 0.9438743127 | 281.894 | 1.26761e+06 | 0.447368 | 0.107595 | 0.173469 | 141 | 21 | false | false | 0.046 |

## Conservative Readout

- Best causal F1: h `1` (gamma `0.5000000000`, F1 `0.712871`).
- Best recall: h `1` (gamma `0.5000000000`, recall `0.683544`).
- Fewest missing relations: h `1` (gamma `0.5000000000`, missing `50`).
- Do metrics agree? Yes under this restricted probe.
- Broad window or sharp optimum? Within this one-run probe, the near-best causal-F1 set is narrow, but the single case, single optimizer seed, and single budget are not enough to separate a sharp optimum from noise.
- Does this organize the previous gamma ambiguity better? This half-life view makes the schedule scale explicit, but with one N=24 case and one optimizer seed it does not yet resolve the previous gamma ambiguity.

## Guardrails

- One N=24 case only.
- One optimizer seed only.
- One budget only.
- Exploratory only.
- No physical constant claim.
- No embeddability claim.
- No theorem.
- Future N=36 use requires careful budget control because medium N=36 is expensive.
