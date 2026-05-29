# Plan: S4-KERR-K18D residual calibration audit

Status: PLAN ONLY (no implementation in this document).
Date: 2026-05-29
Author context: SORKIN-4, Kerr causal-set benchmark continuation.

---

## 0. Purpose and the one question

K18C established (for a single point) that the K17d endpoint-weighted residual has a
**genuine interior local minimum** in `(b, lambda)` for the strongest K18A family —
it is not a boundary/grid artifact. But the interior minimum value is still
`~0.125`, which is large: K17's own tolerances are `W_TOL = 1e-3` and
near-hit `= 10*W_TOL = 1e-2`, so `0.125` is `~125x W_TOL` and `~12.5x` the near-hit
threshold.

The binding open question is therefore no longer *well-posedness* of the objective
but its **scale / meaning**:

> Is `residual ~= 0.125` close to the smallest value this objective can attain
> under its own search procedure, or is it effectively background-level (i.e. the
> "best candidate" is just an ordinary cloud-pair residual)?

K18D answers exactly this, for exactly one fixed point, before any re-ranking of more
families (K18E) or any decision to close/pivot the residual cloud-cloud line.

K18D is a **numerical objective calibration only**. It makes no causal claim.

### Mandatory boundaries (carry into the artifact)

- residual profile != causal_true
- interior minimum != reachability
- low residual != proof
- synthetic reachable target != causal_true (it is a numerical self-consistency anchor)
- candidate_hit != reachability; candidate_miss != spacelike separation
- do NOT touch `cones.py`
- do NOT run any new random cloud-cloud search
- do NOT add another simple filter (no K17F-style step)
- do NOT expand to production geometry
- do NOT sweep spin (stay at the single audited spin)
- do NOT add the outgoing branch (stay ingoing, matching K18C)
- do NOT sell anything on small N

---

## 1. Fixed point (do not change)

Identical to the K18C point:

- `N = 24`
- `seed = 1961`
- `spin_a = 0.50`
- `event_A = index 15`
- `event_B = index 4`  (the real K18A/K18B target; used in Phase 2 only)
- `direction = ingoing` (`-1.0`)
- Reference: `residual(A -> B_real)` at `(b=1.0, lambda=0.5, ingoing)` = `0.12546941635537312`
- Reference A coords (t,r,phi) = (2.9103567417252916, 4.917680370166319, 4.9572224569140015)
- Reference B coords (t,r,phi) = (3.8983298507599984, 4.550588789137428, 5.106522579152854)

---

## 2. Definitions (write these into the artifact verbatim)

- **objective(A, target, b, lambda, direction)** = the K17 weighted residual
  `max(|dt|, |dr|, |dphi_sector_adjusted|)` of the integrated null-geodesic endpoint
  from A versus `target`, with the sector `m` optimized internally. This is exactly
  `k17._eval_trial(...)["endpoint_weighted_residual"]`. No re-definition is introduced.
- **probe_best(A, target, grid, direction)** = `min` of `objective` over a `(b, lambda)`
  grid at fixed `direction`. This is the K17d `_probe_best` pattern, restricted to one
  direction and run on the calibration grid.
- **floor** = the smallest residual `probe_best` can attain for a target that is *known
  to be reachable* from A within the same null-geodesic family (Phase 1).
- **background** = the distribution of `probe_best` over arbitrary already-existing
  cloud events used as targets (Phase 2b). NOT a new random search.

---

## 3. Code to reuse (no reimplementation, no geometry edits)

All already present in the repo:

- `explore.sorkin4_kerr_benchmark.run_kerr_minimal_benchmark` as `kerr`
  - `kerr.kerr_horizon_radius(MASS, spin)`
  - `kerr.generate_exterior_events(N, seed, r_min, equatorial=True)` -> list of `Event`
  - `Event` exposes `.t, .r, .phi, .index`
- `explore.sorkin4_kerr_benchmark.audit_kerr_k17_controlled_candidate_pair_sandbox_001` as `k17`
  - `k17._eval_trial(*, spin, A, B, b, lam, direction)` -> dict with
    `endpoint_weighted_residual`, `endpoint_t_residual`, `endpoint_r_residual`,
    `endpoint_phi_residual_sector_adjusted`, `best_sector_m`, `direction_best`
  - `k17.integrate_to_lambda(spin, b, direction, state0=(t,r,phi), lambda_end)` ->
    dict with `states` (list of `(t,r,phi)`), `rhs`, `failed_reason`
  - `k17.W_TOL`, `k17.MASS`
- `explore.sorkin4_schwarzschild_benchmark.run_schwarzschild_minimal_benchmark` as `schwarz`
  - `schwarz.EXTERIOR_MARGIN`
- `explore.sorkin4_kerr_benchmark.audit_kerr_k18c_residual_boundary_sensitivity_001` as `k18c`
  - reuse its provenance constants / guard logic (regenerate A,B; assert coords + residual)

Precedent to follow (do not import for results, only for construction pattern):
`audit_kerr_k15_synthetic_causal_sandbox_001` and
`audit_kerr_k16_semi_synthetic_pair_sandbox_001` already build *synthetic reachable*
pairs by integrating a geodesic and taking its endpoint. Phase 1 below uses the same
construction, but for **scale calibration**, not for hit-hunting.

### Key reuse boundary (important, prevents scope creep)

`k17._eval_trial` reads only `B.t`, `B.r`, `B.phi` from its target. A synthetic target
therefore needs only an object exposing `.t, .r, .phi` (e.g. `types.SimpleNamespace`).
No `Event` constructor, no geometry module, is required to build a target.

---

## 4. Calibration grid (corrected, interior-containing, from K18C)

K18C found the interior minimum at `(b=1.0, lambda=0.5)` with residual rising on both
sides. Refine *around* that point only (do not widen into new territory):

- `CAL_B_GRID = (0.6, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4)`
- `CAL_LAMBDA_GRID = (0.25, 0.375, 0.5, 0.625, 0.75, 1.0)`
- `direction = ingoing` (fixed)

Edge guard (bounded, in-scope): if any reported best `(b, lambda)` lands on a CAL grid
edge, record it and mark that sub-result `EDGE` — do NOT auto-expand the grid in this
audit; surface it for the next decision instead. (This is what K18C did right.)

---

## 5. Phases

### Phase 0 — Provenance / reproduction guard

- Regenerate the cloud: `r_plus = kerr.kerr_horizon_radius(k17.MASS, 0.50)`,
  `events = kerr.generate_exterior_events(24, 1961, r_plus + schwarz.EXTERIOR_MARGIN, equatorial=True)`.
  `A = {e.index: e}[15]`, `B_real = {e.index: e}[4]`.
- Assert A,B coords match the Section 1 references within `1e-9`.
- Assert `objective(A, B_real, b=1.0, lambda=0.5, ingoing)` matches `0.12546941635537312`
  within `1e-9`.
- **If any guard fails -> overall verdict `INCONCLUSIVE`, stop, write artifact.** Do not
  proceed to Phases 1-3.

Verification: artifact `Provenance` section shows `provenance_ok = True` and the
recomputed residual equal to the reference.

### Phase 1 — Floor (synthetic reachable target recovery)

Goal: the smallest residual the probe can reach for a *known-reachable* target.

Construction set (fixed, ingoing):
- on-grid controls: `(b*, lambda*)` in `{(1.0, 0.5)}`
- off-grid controls: `(b*, lambda*)` in `{(0.95, 0.45), (1.05, 0.55), (1.1, 0.6)}`

For each `(b*, lambda*)`:
1. `run = k17.integrate_to_lambda(spin=0.50, b=b*, direction=-1.0, state0=(A.t, A.r, A.phi), lambda_end=lambda*)`.
   If `run["failed_reason"] is not None`, record `failed` and skip this control.
2. `t*, r*, phi* = run["states"][-1]`; build `B_syn = SimpleNamespace(t=t*, r=r*, phi=phi*, index=-1)`.
3. `floor_value = probe_best(A, B_syn, CAL grid, ingoing)`; record the best `(b, lambda)`
   and whether it is `EDGE`.

Outputs to record:
- `floor_on_grid` = min `floor_value` over on-grid constructions
  (expected ~ integrator/float noise; the objective's intrinsic floor).
- `floor_off_grid` = min `floor_value` over off-grid constructions
  (the realistic floor a generic target faces given CAL spacing).

Verification: artifact `Phase 1` table lists each `(b*, lambda*)`, recovered residual,
recovered `(b, lambda)`, and `EDGE`/`failed` flags. `floor_on_grid` and `floor_off_grid`
are reported as single numbers.

### Phase 2 — Scale placement (real target + background)

- 2a. `residual_real = probe_best(A, B_real, CAL grid, ingoing)`. (Should land near
  `0.125`; if CAL lowers it, record the corrected value and the `(b, lambda)`.)
- 2b. Background, using ONLY the already-generated N=24 cloud (no new sampling):
  deterministically select the target set = every cloud event `T` with `T.index != 15`
  and `T.t > A.t` (forward in time) and `T.r > r_plus + schwarz.EXTERIOR_MARGIN`
  (exterior). For each such `T`, compute `probe_best(A, T, CAL grid, ingoing)`.
  Report `{count, min, median, max}` and where `B_real` falls in this set (rank /
  percentile).

Verification: artifact `Phase 2` section reports `residual_real`, the background
`{count, min, median, max}`, and `B_real`'s rank within the background.

### Phase 3 — Pre-registered verdict

Let `F = floor_off_grid`, `R = residual_real`, `Bmed = background median`,
`Bmin = background min`.

Classify into exactly one:

1. `OBJECTIVE_UNRESOLVED` if `F >= 1e-2`
   (the probe cannot recover even a known-reachable target below the near-hit scale;
   `R ~ 0.125` is then within objective noise and carries no information).
   -> Consequence: the residual lacks resolution at this scale; supports closing /
   redefining the cloud-cloud residual line. Do NOT re-rank.

2. `BACKGROUND_LEVEL` if `F < 1e-2` AND `R >= 0.5 * Bmed`
   (the floor is fine, but the "best candidate" residual is comparable to the typical
   residual of arbitrary cloud events; the candidate is not distinguished).
   -> Consequence: the candidate is ordinary; supports E (redefine what relation the
   residual should detect) and/or A (close cloud-cloud). Do NOT re-rank.

3. `DISTINGUISHED` if `F < 1e-2` AND `R < 0.5 * Bmed` AND `R` is at/near `Bmin`
   (the floor is fine and the candidate is markedly below background).
   -> Consequence: the candidate is genuinely distinguished from background on a
   corrected grid; this is the ONLY outcome that justifies a disciplined K18E
   re-rank of <= 3 deduplicated families on the CAL grid. Still no causal claim.

4. `INCONCLUSIVE` if Phase 0 failed, or all Phase 1 constructions `failed`, or any
   required number is non-finite.

(The `0.5 * Bmed` split is a pre-registered, deliberately conservative threshold; record
the raw numbers so the verdict can be re-derived by hand.)

---

## 6. Stop criteria (hard)

- Fixed scope: one `A` (index 15), one spin (`0.50`), `ingoing` only, the four Phase-1
  constructions, the one deterministic background set, the one CAL grid.
- Do not expand CAL beyond Section 4 (only flag `EDGE`).
- Do not add spin, outgoing, new clouds, new seeds, or new filters.
- If integration fails, record and stop the affected sub-step; do NOT tune the integrator,
  tolerances, or `cones.py`.
- Exactly one artifact pair is produced; no CSV/PNG sweep.

---

## 7. What each verdict authorizes next (decision table, not action now)

- `DISTINGUISHED` -> next step may be **K18E**: re-rank <= 3 deduplicated K18A families on
  the CAL grid (selection criterion: lowest K18A residual after dedup, distinct
  `(seed, event_A, event_B)` signatures, all at `spin_a = 0.50`). Stop K18E if no family
  beats background by the same `0.5 * Bmed` margin.
- `BACKGROUND_LEVEL` -> next step is **E**: reframe what geometric relation the residual
  should encode in Kerr (a new conceptual target), or **A**: close the cloud-cloud line.
- `OBJECTIVE_UNRESOLVED` -> next step is **A/E**: close or redefine; the current objective
  is not measuring at the required scale.
- `INCONCLUSIVE` -> fix the failing guard/construction first; no close/continue decision.

K18D does NOT itself decide causality, reachability, or spacelike separation under any
verdict. Those words appear only as negative caveats.

---

## 8. Deliverable artifacts (when implemented)

- Script: `explore/sorkin4_kerr_benchmark/audit_kerr_k18d_residual_calibration_001.py`
- Output: `explore/sorkin4_kerr_benchmark/kerr_k18d_residual_calibration_001.md`
- Output: `explore/sorkin4_kerr_benchmark/kerr_k18d_residual_calibration_001.json`

Required artifact sections: Status; Fixed point; Definitions; Code reused; Exact commands;
Calibration grid; Provenance guard; Phase 1 floor table (on-grid / off-grid);
Phase 2 real + background table; Verdict
(`DISTINGUISHED` / `BACKGROUND_LEVEL` / `OBJECTIVE_UNRESOLVED` / `INCONCLUSIVE`);
Interpretation; Guardrails; Next operational recommendation.

---

## 9. Single first operational step (when you say go)

Write `audit_kerr_k18d_residual_calibration_001.py` implementing Phase 0 + Phase 1 only
(provenance guard + synthetic-reachable floor), emit the partial artifact, and review the
`floor_on_grid` / `floor_off_grid` numbers **before** running Phases 2-3. If the floor is
already `>= 1e-2` (verdict trending `OBJECTIVE_UNRESOLVED`), stop and report — Phases 2-3
would be moot.
