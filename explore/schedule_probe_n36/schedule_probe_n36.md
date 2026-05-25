# Schedule probe N=36

Exploratory SORKIN-2 diagnostic comparing native geometric cooling schedules for one known-truth case and one optimizer seed.

## Configuration

- Command: `python3 explore/schedule_probe_n36/run_schedule_probe_n36.py`
- Generated at UTC: `2026-05-25T07:30:30+00:00`
- Output directory: `explore/schedule_probe_n36`
- CSV: `explore/schedule_probe_n36/schedule_probe_n36.csv`
- SVG: `explore/schedule_probe_n36/schedule_probe_n36.svg`
- family: `minkowski`
- N: `36`
- d_spacetime: `2`
- target spatial dim: `1`
- case_seed: `1959`
- optimizer_seed: `1987`
- T0: `100.0`
- budget: `medium_25_25_8`
- warmup_limit: `25`
- anneal_limit: `25`
- max_data: `8`
- backend: `cpu`
- Schedules vary only `cooling_factor`: `0.5`, `0.8`, `0.9`, `0.95`.
- This probe uses `ConesSimulator.block_callback` in read-only mode.
- Checkpoint causal metrics use `sim.rold`/`sim.xold`, the accepted state documented by `validation_suite.verify_recovery` and reported by `ConesSimulator.writeout`.
- The historical Bombelli energy, acceptance rule, move set, and internal annealer dynamics are unchanged.

## Summary by schedule

| schedule | gamma | final_energy_eave | final F1 | final recall | final missing | final extra | best F1 seen | block best F1 | min energy seen | block min energy | Edown F1down | Edown recalldown | Edown missingup | Edown extraup |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| gamma_0p5 | 0.5 | 215.629 | 0.373225 | 0.268222 | 251 | 58 | 0.548701 | 6 | 215.629 | 8 | 3 | 3 | 3 | 3 |
| gamma_0p8 | 0.8 | 447.107 | 0.387226 | 0.282799 | 246 | 61 | 0.429688 | 5 | 322.55 | 1 | 2 | 3 | 3 | 3 |
| gamma_0p9 | 0.9 | 504.925 | 0.3762 | 0.285714 | 245 | 80 | 0.447273 | 2 | 322.55 | 1 | 4 | 4 | 4 | 3 |
| gamma_0p95 | 0.95 | 598.076 | 0.354331 | 0.262391 | 253 | 75 | 0.447273 | 2 | 322.55 | 1 | 2 | 2 | 2 | 2 |

## Readout

- Best final causal F1: `gamma_0p8` with F1 `0.387226`.
- Best causal F1 seen anywhere in the trajectory: `gamma_0p5` with F1 `0.548701` at block `6`.
- Minimum energy seen: `gamma_0p5` with Eave `215.629` at block `8`.
- In this run, the schedule with the lowest observed energy matches the best F1 seen along the trajectory, but not the best final F1 schedule.
- Fewest `energy_down_f1_down` steps: `gamma_0p8`, `gamma_0p95` with `2` such steps.
Higher gamma improves at least one audited causal-F1 readout in this run, which is exploratory evidence that cooling rate and thermal mobility matter for this case.
A slower schedule reduces `energy_down_f1_down` counts relative to gamma 0.5 in this run.

## Guardrails

This is exploratory: one N, one family, one case seed, one optimizer seed, one budget, and the unchanged historical objective function.
It is not evidence of embeddability or non-embeddability, does not identify a physical gamma, does not establish an N transition, and does not by itself justify changing the objective function.
