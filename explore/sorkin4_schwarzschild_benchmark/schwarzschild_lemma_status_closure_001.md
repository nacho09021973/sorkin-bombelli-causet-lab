# S4-SCHW-LEMMA-STATUS-CLOSURE-001

Generated: 2026-05-28 — documentation-only synthesis, no new code.

## Purpose

This document records the current known-truth status of the Schwarzschild
causal-boundary / fastest-geodesic lemma front as of the K1–K6 Kerr pause.
It does not add new claims and does not implement anything.

---

## What the lemma is

Every artifact in this directory uses the phrase "not a proof of the
fastest-geodesic lemma" or "not a proof of the causal-boundary lemma" but
none defines the lemma in one place.  From the accumulated context, the
implicit lemma is:

> **Fastest-geodesic lemma (exterior, He-Rideout Appendix B scope):**
> For any exterior pair (E1, E2) with t1 ≤ t2 and both events strictly
> outside the Schwarzschild horizon (r > 2M), the direct monotone-u null
> geodesic from u1=1/r1 to u2=1/r2 on the no-root c² branch is the
> fastest future-directed null geodesic connecting E1 to the stationary
> worldline through E2.  In particular, no one-turn or multi-winding
> outgoing competitor arrives earlier.

This is implied by He & Rideout's Appendix A ("fastest relevant geodesic is
the one with |Δφ|=φ₂, the direct branch") but that appendix proves optimality
among geodesics with the same winding number, not across winding sectors.
The multi-winding suppression is a separate claim.

---

## A. Verified / tested facts

These are items covered by frozen test artifacts and passing regression tests.

### A1. He-Rideout radial/sufficient-bound layer (S4-2)
- Radial exact criteria implemented and regression-tested.
- Exterior spacelike lower bounds (radial and angular f(r)-based) implemented.
- Composed null-curve timelike sufficient bounds implemented.
- All order checks (antisymmetry, transitivity, contradiction) pass across
  N=12, N=16, seeds 1959–1961, margin=0.25/0.35/0.50 (S4-SCHW-STABILITY-001,
  `schwarzschild_stability_001.csv`).

### A2. Direct He-Rideout generic shooting branch (S4-3a/b)
- Direct monotone-u Simpson integrator implemented and validated against
  He & Rideout Table 2 (pair `1 6`, exterior generic example).
- `target_above_direct_max` documented: 5 residual undecided pairs across
  seeds 1959–1968 at N=12; the direct branch cannot reach the required angle.
  These are reported as undecided, not decided false.

### A3. Horizon-crossing plunging branch (direct plunging audit)
- Direct plunging-branch audit: 16 horizon-crossing links across seeds 1..40.
- Angular shooting error ≤ 9.009e-10 rad; Simpson relative drift ≤ 2.155e-12.
- Arrives no later than target event time for all 16 audited links.
- No positive-root obstruction on the accepted c² for any audited link.
- Source: `schwarzschild_horizon_shooting_branch_audit.{json,md}`.
- Regression: `tests/test_sorkin4_schwarzschild_horizon.py`.

### A4. Turning-branch sample audit
- 64 outgoing direct-branch pairs checked (N=12, seeds 1959–1968).
- 0 one-turn branches found that hit the same angular separation.
- Outcome: all 64 pairs are `no_turning_solution`.
- Source: `schwarzschild_exterior_turning_branch_audit.{json,md}`.
- Regression: `tests/test_sorkin4_exterior_turning_branch_audit.py`.

### A5. Aligned horizon-area sweep
- Mean horizon-link count is monotone non-decreasing over M=0.75..2.00.
- Covers N_exterior=16, N_interior=8, 40 seeds, aligned (radial-only) mode.
- Source: `schwarzschild_horizon_area_sweep.{json,md}`.
- Regression: `tests/test_sorkin4_schwarzschild_area_sweep.py`.

### A6. S4-SCHW-STABILITY-001 sweep
- Turning-branch gap positive for all 18 cells (N=12/16, seeds 1959–1961,
  margin=0.25/0.35/0.50).  Min gap across all cells: 0.1402.
- Source: `schwarzschild_stability_001.{csv,json,md}`.
- Regression: `tests/test_sorkin4_schw_stability_001.py`.

---

## B. Analytic / asymptotic arguments

These are pen-and-paper or closed-form arguments; not executable tests.

### B1. Ingoing branch turning-point exclusion
The genuine one-turn competitor exists only for outgoing pairs
(u2 < u1 < u_turn < 1/(3M)).  For ingoing pairs (u2 > u1), the turning
point is the maximum u, so a turned branch cannot later reach a larger u.
This eliminates ingoing one-turn competitors by construction.
Source: README S4 Exterior Turning-Branch Competitor Audit section; code
comment in `audit_exterior_turning_branch.py`.

### B2. Asymptotic gap in the weak-field nearby-radius corner
For the outgoing case with r1=R, r2=R+ε, R >> M, ε/R << 1:

```
phi_direct_max  ~ 3√3 M ε / R²        [O(M ε / R²)]
phi_turning_min ~ √(2 ε / R)           [O(√(ε/R))]
gap             ~ √(2 ε / R) > 0       [leading term, positive]
```

The ratio of direct-max to turning-min is O(M√ε / R^(3/2)), small in the
stated regime.  The gap falls as R^(-1/2) for fixed small ε.
Source: `schwarzschild_exterior_turning_asymptotic_note.md`.

This is a local asymptotic argument, not a global theorem.

---

## C. Numerical exploratory evidence

These are finite-grid diagnostics, not proofs.

### C1. Phase-space angular disjointness audit
- Grid: r1=3.1..12.0, r2≤20.0, step=0.1 (11,205 points, r2 > r1 > 3M).
- All 11,205 points: `disjoint_ranges`.  Min gap: 0.1446.
- Conclusion: on this finite grid, the direct branch (c² > c²_crit) and the
  one-turn branch (c² < c²_crit) serve disjoint angular intervals; they
  compete for no target angle.
- Source: `schwarzschild_exterior_turning_phase_space_audit.{json,md}`.
- Regression: `tests/test_sorkin4_exterior_turning_phase_space_audit.py`.

### C2. Asymptotic corner audit
- Grid: R ∈ {6,8,10,12,16,20,30,50,80,120}, ε ∈ {1e-3,3e-3,...,1} (70 pts).
- All 70 points: `disjoint_ranges`.  Min gap: 0.00429 at R=120, ε=0.001.
- Source: `schwarzschild_exterior_turning_asymptotic_audit.{json,md}`.
- Regression: `tests/test_sorkin4_exterior_turning_asymptotic_audit.py`.

---

## D. Not yet shown

These are the genuine open items, stated explicitly.

### D1. Proof of angular disjointness for all exterior r1, r2 > 2M

The phase-space audit (C1) covers r1 ∈ [3.1, 12.0] and r2 ≤ 20.0.  It does
not cover:
- r1 close to 2M (near-horizon exterior, 2M < r1 < 3.1M).
- r2 > 20M (large-separation pairs).
- r1 close to r2 for r1 > 12M (asymptotic near-equal-radius).

The asymptotic note (B2) covers the weak-field nearby-radius corner in the
R >> M regime.  Together C1 and B2 leave a gap for near-photon-sphere and
moderate-far-exterior pairs.

### D2. Earlier-arrival proof across winding sectors

The phase-space audit shows disjoint angular ranges, which implies the two
branches compete for different angular targets — but it does not directly
compare arrival times for a pair that has the same target angle and could be
reached by either branch.  In practice, if angular ranges are always disjoint,
no pair can ever be reached by both, so an arrival-time comparison is vacuous.
Confirming this vacuousness is itself a logical step that is not yet recorded
as a standalone audit.

### D3. Multi-winding geodesics (|Δφ| > π)

The turning-branch audit and phase-space audit cover the one-turn outgoing
competitor only.  Multi-winding geodesics (winding number ≥ 2) are not audited.
He & Rideout's Appendix A excludes them by arguing the fastest geodesic has the
minimum |Δφ|, but this is not verified numerically in this project.

### D4. Near-horizon exterior pairs (2M < r < 3.1M)

The phase-space audit's lower bound r1 > 3.1M was set to stay above the photon
sphere r = 3M with margin.  Near-horizon exterior pairs with 2M < r1 < 3M are
geometrically distinct: the photon sphere affects turning-point structure in
ways that the current grid does not probe.

### D5. Generic-pair coverage for `target_above_direct_max` residuals

As documented in S4-3b (README), the 5 residual undecided pairs have
`phi_target > phi_direct_max`.  These are left as `None` (undecided), which is
correct for the He-Rideout direct-branch procedure.  No additional branch or
bound has been implemented to decide them.

---

## a/b/c route options: status

**The a/b/c route options mentioned in conversation are NOT documented in the
repository.**  No file in `explore/`, `docs/`, or test code records a labelled
a/b/c or option-A/B/C choice for the next step on the lemma front.  This note
cannot recover them from the repo and will not fabricate them.

The following are **reconstructed routes**, not Claude's original options:

**Route 1 — Extend the numerical phase-space grid.**
Extend C1 to cover r1 ∈ (2M, 3M) (near-photon-sphere exterior) and r2 up to
50M or 100M.  Straightforward: the `audit_exterior_turning_phase_space.py`
script is already parameterized.  Expected outcome: min gap decreases near
r1→3M (photon sphere) but remains positive.  This addresses D1 for the missing
near-horizon and large-separation corners, but still leaves a finite-grid
caveat.  Cost: low (one script run + artifact).

**Route 2 — Analytic proof of angular disjointness for all r > 2M.**
Derive a closed-form lower bound on `phi_turning_min - phi_direct_max` that is
positive for all r1, r2 > 2M.  The asymptotic note (B2) gives the weak-field
piece; the photon-sphere piece would need a separate local expansion near r=3M.
Combining them would give a global argument.  This directly addresses D1 and
D2.  Cost: analytical work, no code; result would be a documentation file.

**Route 3 — Accept the current state and advance to the Dou-Sorkin L0 horizon
program.**
Per `docs/sorkin4_level0_scope.md`, the active scientific front is the
horizon-molecule count (N_links ∝ A = 16πM²), not the turning-branch lemma.
The current exterior solver already reports residual `target_above_direct_max`
pairs as undecided (D5), which is the correct output.  The shooting step and
turning-branch audit have sufficient coverage for the exterior pairs actually
produced by the benchmark.  Moving forward to a full interior+exterior solver
and the Dou-Sorkin link count is independent of closing D1–D4 completely.
Cost: deferral of lemma closure; next step is the interior IEF solver.

---

## Recommendation

**Route 3 is the appropriate next move**, for the following reasons:

1. The active scientific front, per `docs/sorkin4_level0_scope.md`, is L0
   (Dou-Sorkin horizon molecules), not turning-branch completeness.
2. The turning-branch gap (D1–D4) does not block the interior solver or the
   horizon-link count.  The exterior solver correctly reports `None` for
   cases it cannot decide; those are a small minority (≤1/N*(N-1)/2 per seed).
3. Route 1 would extend numerical evidence but not close the analytic gap.
   Route 2 would close D1–D2 analytically but is a detour from L0.
4. The turning-branch evidence (A4, C1, C2, B2) is sufficient to mark the
   sample behavior as "no turning competitor observed; disjoint ranges confirmed
   on large exterior grid; asymptotic gap positive in dangerous weak-field
   corner."  That is an honest diagnostic description, not a proof gap that
   blocks progress.

If the analytic closure of D1 is desired before L0, Route 2 can be done as a
documentation-only file (no new script), extending the asymptotic note to cover
the near-photon-sphere limit.  That is a bounded analytic exercise.

---

## Caveats

- This document does not claim a global proof of the fastest-geodesic lemma.
- It does not claim Kerr causal validation of any kind.
- It does not claim Hawking/Bekenstein discrete rediscovery.
- The Schwarzschild scope is restricted to the exterior region r > 2M unless
  files explicitly say otherwise.
- "Verified" (section A) means: frozen artifact + passing regression test.
  It does not mean mathematically proved.
- "Analytic argument" (section B) means: pen-and-paper asymptotic expansion
  recorded in a documentation file.  It has not been peer-reviewed.
- "Numerical evidence" (section C) means: finite-grid audit with positive
  result.  It is falsifiable by extending the grid.

---

## Proposed next implementation (if Route 3 is accepted)

No new script yet — this file is documentation only.

The implementation proposal for Route 3 is:

1. Implement a full interior IEF event generator and radial causal criterion
   for r < 2M, as specified in `docs/sorkin4_level0_scope.md` items 1–3.
2. Implement a horizon-link counter (covering relations with one end in
   r > 2M, one end in r < 2M) as specified in item 4.
3. Run the Dou-Sorkin sweep: vary M across a range, record mean horizon-link
   count, verify N_links ∝ 16πM².

This is the Level 0 causal-set black-hole diagnostic.  It does not require
closing D1–D4.
