# S4-KERR-K3-LOCAL-CONE-001: Kerr Equatorial Local Null-Cone Diagnostic

Generated: 2026-05-27T12:35:11.963410+00:00

## What this is

This is a **local metric-sign diagnostic**, not a global Kerr causal solver.

It computes Boyer-Lindquist equatorial metric coefficients and evaluates
the quadratic interval `ds²` at the midpoint radius for small-displacement
equatorial pairs. Each evaluated pair is labelled as one of:

- `timelike_local_candidate`  — `ds² < -tol`
- `nullish_local_candidate`   — `|ds²| ≤ tol`
- `spacelike_local_candidate` — `ds² > tol`

**It does NOT:**

- Establish null geodesic connectivity between the two events.
- Integrate Kerr geodesics of any kind.
- Decide prograde or retrograde causal relations.
- Constitute a Kerr causal solver of any kind.

This is only the first local consistency check before any Kerr causal inference.

## Parameters

- M = 1.0, theta = pi/2 (equatorial), spins = [0.0, 0.25, 0.5, 0.75]
- N = 12, seed = 1959, margin = 0.35
- Local-pair filter: |dr| ≤ 1.0 (BL units), |dphi| ≤ 0.5 rad (~28.6°)

## Summary

| Check | Result |
|-------|--------|
| **all_checks_pass** | **True** |
| positive_spin_cases_all_undecided | True |

## Per-Spin Results

| a | r_+ | min_Δ | min_g_rr | min_g_φφ | evaluated | timelike | nullish | spacelike | global_true | undecided | checks |
|---|-----|-------|----------|----------|-----------|----------|---------|-----------|-------------|-----------|--------|
| 0.00 | 2.0000 | 1.4280 | 1.5364 | 6.5444 | 4 | 1 | 0 | 3 | 1 | 5 | True |
| 0.25 | 1.9682 | 1.4095 | 1.5320 | 6.5228 | 4 | 1 | 0 | 3 | 0 | 66 | True |
| 0.50 | 1.8660 | 1.3491 | 1.5190 | 6.4510 | 4 | 1 | 0 | 3 | 0 | 66 | True |
| 0.75 | 1.6614 | 1.2238 | 1.4978 | 6.2931 | 4 | 1 | 0 | 3 | 0 | 66 | True |

## Interpretation

- `evaluated` = pairs passing both displacement filters (|dr|≤1.0, |dphi|≤0.5 rad).
- `global_true` = global causal assertions. For a=0 these are the Schwarzschild
  control counts. For a>0 this is always 0 (all pairs remain undecided globally).
- `undecided` = global causal pairs not decided by this diagnostic.
- **Local labels are not global causal decisions.** A `timelike_local_candidate`
  label means `ds² < 0` at the midpoint radius under the local BL metric.
  It does not imply global causal reachability.
