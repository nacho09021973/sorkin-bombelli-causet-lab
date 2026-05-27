# S4-KERR-K3-LOCAL-CONE-001 Design

Date: 2026-05-27

## What it is

A conservative Kerr equatorial **local metric-sign diagnostic**. K3 computes
the Boyer-Lindquist quadratic interval `ds²` at a midpoint radius for
small-displacement equatorial pairs, classifies each as
`timelike_local_candidate / nullish_local_candidate / spacelike_local_candidate`,
and preserves the K1/K2 invariant that for a>0 **no global causal relations are
decided**.

## Line boundary

```
local_timelike_candidate  ≠  causal_relation_true
local_spacelike_candidate ≠  causal_relation_false
```

This must be stated explicitly in the script, the MD artifact, and the test.

## Physics setup

- Boyer-Lindquist equatorial plane, theta = pi/2, M = 1
- Spins a = 0.0, 0.25, 0.5, 0.75  (all sub-extremal)
- Horizon: r_plus = M + sqrt(M² - a²)
- Ergosphere equatorial: r_erg_eq = 2M
- Delta = r² - 2Mr + a²    (computed explicitly before kerr_metric_components)
- Sigma_eq = r²             (theta=pi/2)

Equatorial metric coefficients (signature -+++):

  g_tt     = -(1 - 2M/r)
  g_tphi   = -2Ma/r
  g_rr     = r²/Delta
  g_phiphi = r² + a² + 2Ma²/r

dtheta = 0 by construction (all events at theta=pi/2).

## Local-pair filter (Adjustment 1)

Only "small-displacement" pairs are classified as local cone samples.
Threshold: the Euclidean coordinate displacement must satisfy

  |dr| < DR_LOCAL_THRESHOLD   and   |dphi| < DPHI_LOCAL_THRESHOLD

where thresholds are set empirically to keep the sample physically local.
Pairs outside the threshold are counted but not classified; their count is
exposed as `local_skipped_pair_count`.

Fields added:
  local_evaluated_pair_count
  local_skipped_pair_count
  local_max_abs_dr
  local_max_abs_dphi
  local_max_abs_dt

Test invariant: local_timelike_count + local_nullish_count + local_spacelike_count
                = local_evaluated_pair_count  (NOT N*(N-1)/2)

## Delta pre-computation (Adjustment 2)

For each event, compute Delta = r² - 2Mr + a² before calling
kerr_metric_components.  Record min_Delta in the row.  This makes the
exterior-check explicit and traceable in the artifact even if the sampling
changes.

## Causal accounting (Adjustment 3)

a=0: preserve Schwarzschild control counts (reproduce K1/K2 a=0 row:
     true=1, false=64, undecided=1 for N=12, seed=1959, default margin).
a>0: global_true_relations=0, global_false_relations=0,
     global_undecided_pairs=N*(N-1)/2.

## CSV row fields (one row per spin)

spin_a, M, N, seed, margin, r_plus, r_erg_eq, r_min_observed,
all_points_exterior, min_Delta, min_g_rr, min_g_phiphi, max_abs_g_tphi,
local_evaluated_pair_count, local_skipped_pair_count,
local_max_abs_dr, local_max_abs_dphi, local_max_abs_dt,
local_timelike_count, local_nullish_count, local_spacelike_count,
global_true_relations, global_false_relations, global_undecided_pairs,
schwarzschild_reduction_pass, frame_dragging_sign_pass, all_checks_pass

## Files created

- explore/sorkin4_kerr_benchmark/audit_kerr_k3_local_cone_001.py
- explore/sorkin4_kerr_benchmark/kerr_k3_local_cone_001_n12_seed1959.{csv,json,md}
- tests/test_sorkin4_kerr_k3_local_cone.py
- explore/sorkin4_kerr_benchmark/README.md  (K3 section appended)

## Files NOT touched

- run_schwarzschild_minimal_benchmark.py
- run_kerr_minimal_benchmark.py  (no changes; imports only)
- cones.py
- Any existing benchmark CSVs
