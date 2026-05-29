# S4-KERR-K18D residual calibration audit 001 (Phase 0 + Phase 1)

## Status

PARTIAL — Phase 0 ok; Phase 1 complete; floor_off_grid = 0.0416334976665 >= 1e-2 → trending OBJECTIVE_UNRESOLVED

This is a PARTIAL artifact. Phase 2 (background distribution) and Phase 3 (pre-registered verdict) are not yet run. Review floor numbers before proceeding.

## Fixed point

- N = 24
- seed = 1961
- spin_a = 0.5
- event_A = 15
- event_B = 4  (B_real; used in Phase 2 only, not used here)
- direction = ingoing (fixed throughout)
- Reference residual (b=1.0, lambda=0.5, ingoing) = 0.12546941635537312

## Definitions

- **objective(A, target, b, lambda, direction)** = K17 weighted residual `max(|dt|, |dr|, |dphi_sector_adjusted|)` of the integrated null-geodesic endpoint from A versus target, with sector m optimized internally. Exact reuse of `k17._eval_trial(...)['endpoint_weighted_residual']`.
- **probe_best(A, target, grid, direction)** = `min` of `objective` over the CAL (b, lambda) grid at fixed direction.
- **floor_on_grid** = min `probe_best` over on-grid constructions (b_star, lambda_star coincide with CAL grid nodes). Measures integrator/float-noise.
- **floor_off_grid** = min `probe_best` over off-grid constructions (b_star, lambda_star do not coincide with CAL grid nodes). Measures the realistic floor a generic target faces given CAL spacing.

## Code reused

- `explore/sorkin4_kerr_benchmark/run_kerr_minimal_benchmark.py` (`kerr_horizon_radius`, `generate_exterior_events`, `Event`)
- `explore/sorkin4_kerr_benchmark/audit_kerr_k17_controlled_candidate_pair_sandbox_001.py` (`_eval_trial`, `integrate_to_lambda`, `W_TOL`, `MASS`)
- `explore/sorkin4_schwarzschild_benchmark/run_schwarzschild_minimal_benchmark.py` (`EXTERIOR_MARGIN`)
- K18C provenance constants reproduced (same fixed point, same reference values).

`cones.py` and production geometry were NOT modified or imported.

## Exact commands

```bash
python3 -m explore.sorkin4_kerr_benchmark.audit_kerr_k18d_residual_calibration_001
```

## Calibration grid

- CAL_B_GRID = (0.6, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4)
- CAL_LAMBDA_GRID = (0.25, 0.375, 0.5, 0.625, 0.75, 1.0)
- direction = ingoing (fixed)
- Edge guard: if best (b, lambda) sits on a CAL grid boundary, result is flagged EDGE. Do NOT auto-expand the grid; surface it for the next decision.

## Phase 0: Provenance guard

- r_plus = 1.86602540378
- n_events regenerated = 24
- A coords (t,r,phi) = (2.91035674173, 4.91768037017, 4.95722245691) — matches K18C: True
- B_real coords (t,r,phi) = (3.89832985076, 4.55058878914, 5.10652257915) — matches K18C: True
- recomputed residual at (b=1.0, lambda=0.5, ingoing) = 0.125469416355 (reference: 0.125469416355); reproduced = True
- **provenance_ok = True**

## Phase 1: Floor (synthetic reachable targets)

Each construction integrates from A at (b_star, lambda_star, ingoing) to get B_syn, then runs probe_best over the full CAL grid (ingoing only). The floor is the minimum probe_best residual across constructions of each type.

### On-grid constructions

| b_star | lambda_star | status | probe_best_residual | probe_best_b | probe_best_lambda | edge |
|---:|---:|---|---:|---:|---:|---|
| 1 | 0.5 | ok | 0 | 1 | 0.5 | False |

### Off-grid constructions

| b_star | lambda_star | status | probe_best_residual | probe_best_b | probe_best_lambda | edge |
|---:|---:|---|---:|---:|---:|---|
| 0.95 | 0.45 | ok | 0.0850341770677 | 1.4 | 0.5 | True |
| 1.05 | 0.55 | ok | 0.0857954523347 | 0.6 | 0.5 | True |
| 1.1 | 0.6 | ok | 0.0416334976665 | 1.4 | 0.625 | True |

**floor_on_grid  = 0**
**floor_off_grid = 0.0416334976665**
(near-hit threshold = 1e-2; W_TOL = 0.001)

### Phase 1 interpretation

floor_off_grid = 0.0416334976665 >= 1e-2 (near-hit threshold). The probe cannot recover even a known-reachable off-grid target below the near-hit scale. This indicates OBJECTIVE_UNRESOLVED: the residual does not measure at the required scale, and the 0.125 value carries no information. Phase 2+3 would be moot — recommend stopping here.
floor_on_grid = 0 (integrator/float-noise floor; expected near machine epsilon for on-grid construction).

## Guardrails

- residual profile != causal_true
- interior minimum != reachability
- low residual != proof
- synthetic reachable target != causal_true (it is a numerical self-consistency anchor)
- candidate_hit != reachability; candidate_miss != spacelike separation
- K18D is a numerical objective calibration only; it makes no causal claim
- floor_off_grid measures the floor the probe faces for a generic target at CAL spacing
- floor_on_grid measures the integrator/float-noise floor (expected near machine epsilon)

## Next operational recommendation

floor_off_grid = 0.0416334976665 >= 1e-2. The objective is unresolved at the required scale. Recommended action: close or redefine the residual cloud-cloud line (OBJECTIVE_UNRESOLVED verdict; Phase 2+3 would be moot). Do NOT re-rank candidates.

