# S4-KERR-K18C residual boundary-sensitivity diagnostic 001

## Status

Completed (single-candidate boundary-sensitivity diagnostic). Overall classification: **INTERIOR** (b-cut: INTERIOR, lambda-cut: INTERIOR).

## Fixed candidate

- N = 24
- seed = 1961
- spin_a = 0.5
- event_A = 15
- event_B = 4
- direction = ingoing (probe direction held fixed at K18B best_direction)
- K17d best_b = 1.0 (upper edge of old PROBE_B_GRID)
- K17d best_lambda = 0.5 (lower edge of old PROBE_LAMBDA_GRID)
- sector m: optimized internally by k17._eval_trial exactly as in K17d (recorded as best_sector_m).

## Input artifacts / code actually used

- `explore/sorkin4_kerr_benchmark/audit_kerr_k17_controlled_candidate_pair_sandbox_001.py` (`k17._eval_trial` — the residual evaluation, reused unchanged)
- `explore/sorkin4_kerr_benchmark/audit_kerr_k17d_cloud_size_seed_scan_001.py` (`_probe_best` call pattern and PROBE grids reproduced)
- `explore/sorkin4_kerr_benchmark/run_kerr_minimal_benchmark.py` (`kerr_horizon_radius`, `generate_exterior_events`, `Event`)
- `explore/sorkin4_schwarzschild_benchmark/run_schwarzschild_minimal_benchmark.py` (`EXTERIOR_MARGIN`)
- `explore/sorkin4_kerr_benchmark/kerr_k18b_single_candidate_geometric_audit_001.md` (candidate provenance / reference values)

cones.py and production geometry were NOT modified or imported beyond the existing K17 audit stack.

## Exact commands used

```bash
python3 -m explore.sorkin4_kerr_benchmark.audit_kerr_k18c_residual_boundary_sensitivity_001
```

## Provenance / reproduction guard

- regenerated cloud: N=24, seed=1961, spin_a=0.5, n_events=24, r_plus=1.86602540378
- event_A coords (t,r,phi) = (2.91035674173, 4.91768037017, 4.95722245691) — matches K18B: True
- event_B coords (t,r,phi) = (3.89832985076, 4.55058878914, 5.10652257915) — matches K18B: True
- recomputed residual at (b=1.0, lambda=0.5, ingoing) = 0.125469416355 (K18B reference 0.125469416355); reproduced = True

## Probe grid used for b and lambda

- B_CUT_GRID = (-0.5, 0.0, 0.5, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0) (lambda fixed at 0.5)
- LAMBDA_CUT_GRID = (0.03125, 0.0625, 0.125, 0.25, 0.5, 1.0, 2.0) (b fixed at 1.0)
- CORNER_B_GRID = (0.5, 1.0, 1.5, 2.0)
- CORNER_LAMBDA_GRID = (0.125, 0.25, 0.5, 1.0)
- (reference) K17d PROBE_B_GRID = (-1.0, -0.5, 0.0, 0.5, 1.0)
- (reference) K17d PROBE_LAMBDA_GRID = (0.5, 1.0, 2.0)

## 1D b-cut table

Fixed: lambda = 0.5, direction = ingoing. Extension direction: increasing b beyond old upper edge 1.0.

| b | weighted_residual | best_sector_m | t_residual | r_residual | phi_residual_sector_adj | note |
|---:|---:|---:|---:|---:|---:|---|
| -0.5 | 0.152064657187 | 0 | -0.110625405502 | -0.137836879073 | -0.152064657187 | in K17d grid |
| 0 | 0.140777070426 | 0 | -0.114942394894 | -0.137008915454 | -0.140777070426 | in K17d grid |
| 0.5 | 0.132908418971 | 0 | -0.119485382039 | -0.132908418971 | -0.129512926228 | in K17d grid |
| 1 | 0.125469416355 | 0 | -0.124239835393 | -0.125469416355 | -0.11829642045 | old best_b (upper edge); in K17d grid |
| 1.25 | 0.126692588759 | 0 | -0.126692588759 | -0.120460405047 | -0.112713845497 | extension |
| 1.5 | 0.129194058286 | 0 | -0.129194058286 | -0.11456280398 | -0.107152740234 | extension |
| 2 | 0.134339543208 | 0 | -0.134339543208 | -0.0999828745228 | -0.0961088055523 | extension |
| 2.5 | 0.139671755012 | 0 | -0.139671755012 | -0.0814232878119 | -0.0851943030434 | extension |
| 3 | 0.145191612238 | 0 | -0.145191612238 | -0.0584351183212 | -0.0744432660713 | extension |
| 4 | 0.156844649011 | 0 | -0.156844649011 | 0.00386133469726 | -0.053607511498 | extension |

b-cut classification: **INTERIOR**

## 1D lambda-cut table

Fixed: b = 1.0, direction = ingoing. Extension direction: decreasing lambda below old lower edge 0.5.

| lambda | weighted_residual | best_sector_m | t_residual | r_residual | phi_residual_sector_adj | note |
|---:|---:|---:|---:|---:|---:|---|
| 0.03125 | 0.935780055206 | 0 | -0.935780055206 | 0.336265777143 | -0.147581908817 | extension |
| 0.0625 | 0.883371269397 | 0 | -0.883371269397 | 0.305444965046 | -0.145837612179 | extension |
| 0.125 | 0.777888649986 | 0 | -0.777888649986 | 0.243818667618 | -0.142268303763 | extension |
| 0.25 | 0.564130267405 | 0 | -0.564130267405 | 0.120629955801 | -0.134788415312 | extension |
| 0.5 | 0.125469416355 | 0 | -0.124239835393 | -0.125469416355 | -0.11829642045 | old best_lambda (lower edge); in K17d grid |
| 1 | 0.817590889728 | 0 | 0.817590889728 | -0.616338892994 | -0.077421219285 | in K17d grid |
| 2 | 3.14968093396 | 0 | 3.14968093396 | -1.59012979095 | 0.0642904939157 | in K17d grid |

lambda-cut classification: **INTERIOR**

## 2D corner table

Direction = ingoing. Small 2D grid around the corner (b=1.0, lambda=0.5). Cell = weighted_residual.

| b \ lambda | 0.125 | 0.25 | 0.5 | 1 |
|---:|---:|---:|---:|---:|
| 0.5 | 0.776951190596 | 0.562106450398 | 0.132908418971 | 0.831244273734 |
| 1 | 0.777888649986 | 0.564130267405 | 0.125469416355 | 0.817590889728 |
| 1.5 | 0.778836510477 | 0.566198231774 | 0.129194058286 | 0.802894538159 |
| 2 | 0.779794587693 | 0.568309230433 | 0.134339543208 | 0.787250131769 |

## Classification: INTERIOR / RUNAWAY / INCONCLUSIVE

- b-cut: **INTERIOR**
- lambda-cut: **INTERIOR**
- overall: **INTERIOR**

Pre-registered rule:
- RUNAWAY if the minimum residual on a cut sits at the extreme grid point in the extension direction (residual keeps improving off the new edge).
- INTERIOR if the minimum is strictly interior / on the non-extension side (extending past the old K17d edge does not keep lowering the residual).
- INCONCLUSIVE if a cut has fewer than three finite residuals, or if the provenance/reproduction guard fails.
- overall = INCONCLUSIVE if any cut is INCONCLUSIVE; else RUNAWAY if any cut is RUNAWAY; else INTERIOR.

## Interpretation

The residual attains an interior minimum within the extended b/lambda range (it stops improving once the old K17d edge is passed). Under the pre-registered reading this indicates the old K17d probe grid was too narrow, so K17d/K18A residual rankings are NOT final and would need re-evaluation on a corrected grid. This is an objective-grid statement only; no causal/reachability claim follows (see guardrails).

## Guardrails

- residual profile != causal_true
- interior minimum != reachability
- runaway residual != spacelike separation
- K18C is a numerical objective diagnostic only
- no global causal claim is made
- boundary-sensitivity of a residual is a statement about the objective's parametrization, not about Kerr causal structure
- this audit is local to ONE K18A/K18B candidate family

## Next operational recommendation

Single next step: treat K17d/K18A residual rankings as non-final and re-evaluate the top candidates on a corrected (interior-containing) b/lambda grid before any further interpretation. Still no causal/reachability claim.

