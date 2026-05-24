# Exploratory N=24 Gamma Probe

**Status:** exploratory result only. This is not confirmation and not publication-ready.

This is one N=24 known-truth case generated through the existing SORKIN-2 Minkowski sprinkling machinery, using one case seed and one optimizer seed. It uses the documented short Phase 2B budget only. It is not evidence for a physical constant, not evidence for causal-set embeddability, and only a probe of algorithmic accessibility in the historical Bombelli annealer.

## Provenance

- Output directory: `explore/gamma_n24_probe`
- Runner: `explore/gamma_n24_probe/run_gamma_n24_probe.py`
- Command: `python3 explore/gamma_n24_probe/run_gamma_n24_probe.py`
- Generated at UTC: `2026-05-24T08:14:00+00:00`
- CSV: `explore/gamma_n24_probe/gamma_n24_probe.csv`
- Plot: `explore/gamma_n24_probe/gamma_n24_probe.svg`
- Family: `minkowski`
- N: `24`
- spacetime dimension: `2`
- target embedding dimension: `1`
- case seed: `1959`
- optimizer seed: `1987`
- T0: `100.0`
- warmup_limit: `10`
- anneal_limit: `10`
- max_data: `4`
- backend: `cpu`
- total runtime seconds: `0.138`
- budget source: Phase 2B short schedule.
- N=36 note: the earlier N=36 scan was too expensive for this exploratory pass and produced no CSV/markdown/SVG/log artifact before interruption or disappearance from the process table.

## Result

- Best gamma by final energy: `0.562373`
- Best final energy: `130.465017`
- Best energy gap: `130.465017`
- Best interval RMSE: `9137960.630239`
- Window readout: broad by the exploratory 10% final-energy criterion.
- Thermal half-life at best gamma: `1.204` annealing temperature steps.

## Gamma Table

| gamma | final_energy | energy_gap | interval_rmse | success | runtime_s |
| ---: | ---: | ---: | ---: | :---: | ---: |
| 0.500000 | 131.822227 | 131.822227 | 71268273.480146 | no | 0.014 |
| 0.530270 | 131.092389 | 131.092389 | 340916.400716 | no | 0.015 |
| 0.562373 | 130.465017 | 130.465017 | 9137960.630239 | no | 0.014 |
| 0.596419 | 135.645205 | 135.645205 | 6106828.456190 | no | 0.012 |
| 0.632527 | 135.645205 | 135.645205 | 6106828.456190 | no | 0.012 |
| 0.670820 | 135.645205 | 135.645205 | 6106828.456190 | no | 0.012 |
| 0.711432 | 151.237389 | 151.237389 | 344640.348984 | no | 0.012 |
| 0.754503 | 139.465256 | 139.465256 | 268366.496860 | no | 0.012 |
| 0.800181 | 139.321485 | 139.321485 | 7210040.597804 | no | 0.011 |
| 0.848624 | 142.986146 | 142.986146 | 113319.222524 | no | 0.011 |
| 0.900000 | 143.563028 | 143.563028 | 729055.731927 | no | 0.011 |

## Required Questions

- Is the best gamma near 0.8? In this one-run probe, best gamma is `0.562373`.
- Is the optimum broad or narrow? broad by the exploratory 10% final-energy criterion.
- Does the result suggest a thermal half-life scale? The best-gamma half-life is `1.204` steps, but this single case and single optimizer seed cannot establish a scale.
- Does the result justify a follow-up sweep over seeds/topologies? It can motivate one only if treated as an algorithmic accessibility follow-up, not as physical evidence.

## Interpretation Guardrails

- Exploratory result only.
- One N=24 case.
- One case seed and one optimizer seed.
- Not evidence for a physical constant.
- Not evidence for causal-set embeddability.
- Only a probe of algorithmic accessibility in the historical Bombelli annealer.
