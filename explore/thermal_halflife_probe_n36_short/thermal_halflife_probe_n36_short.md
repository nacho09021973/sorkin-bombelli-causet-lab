# N=36 Short-Budget Thermal Half-Life Mini-Probe

## Status

Exploratory only; not confirmation.

This probe scans h where gamma = 2^(-1/h), so h is the number of annealing steps required to halve the temperature.

## Run Configuration

- Command: `python3 explore/thermal_halflife_probe_n36_short/run_thermal_halflife_probe_n36_short.py`
- Generated UTC: `2026-05-24T11:46:44+00:00`
- Family: `minkowski`
- N: `36`
- spacetime dimension: `2`
- T0: `100.0`
- case seed: `1959`
- optimizer seed: `1987`
- budget label: `short_10_10_4`
- budget: warmup `10`, anneal `10`, max data `4`
- timeout policy: `900.0` seconds per h

## h/gamma Table

| h | gamma |
|---:|---:|
| 1 | 0.5000000000 |
| 3 | 0.7937005260 |
| 6 | 0.8908987181 |
| 10 | 0.9330329915 |

## Results

| h | gamma | final_energy | interval_rmse | precision | recall | F1 | missing | extra | exact_match | success | runtime_s |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|:---:|---:|
| 1 | 0.5000000000 | 251.424 | 413613 | 0.661692 | 0.387755 | 0.488971 | 210 | 68 | false | false | 0.062 |
| 3 | 0.7937005260 | 263.937 | 1.72351e+07 | 0.516340 | 0.230321 | 0.318548 | 264 | 74 | false | false | 0.029 |
| 6 | 0.8908987181 | 285.842 | 549583 | 0.520000 | 0.113703 | 0.186603 | 304 | 36 | false | false | 0.027 |
| 10 | 0.9330329915 | 276.526 | 8.01363e+06 | 0.513514 | 0.166181 | 0.251101 | 286 | 54 | false | false | 0.025 |

## Required Readout

- Which h gives best causal F1? h `1` (gamma `0.5000000000`, F1 `0.488971`).
- Which h gives best recall? h `1` (gamma `0.5000000000`, recall `0.387755`).
- Which h minimizes missing relations? h `1` (gamma `0.5000000000`, missing `210`).
- Does any h improve over h=1? `no` under these causal-order criteria.
- Is there evidence that N=36 prefers a longer thermal half-life than N=24? In this restricted N=36 mini-probe, the primary causal-order criteria do not prefer h > 1. This is only weak exploratory evidence because it is one N=36 case, one optimizer seed, and short budget only.

## Guardrails

- Exploratory only.
- One N=36 case only.
- One optimizer seed only.
- Short budget only.
- No embeddability claim.
- No physical constant claim.
- No inference from annealer failure to non-existence.
- No inference from low final energy to manifoldlikeness.
