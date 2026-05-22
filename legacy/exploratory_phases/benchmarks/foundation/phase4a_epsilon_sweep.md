# Phase 4A — Epsilon sweep + morphology diagnostic

**Status:** experimental-design probe.  No PySR is run here.
Phase 4B is gated on this verdict.

## Verdict (automatic)

**AMBIGUOUS**

3/6 curves show monotonic ε dependence — between the strong and null thresholds. Inspect the per-curve table below and decide manually whether a V-like or interior-minimum morphology concentrates in a particular (n, target_dim) cell worth probing further.

## Configuration

- Warmup mode: guarded_warmup only (threshold=0.0, limit=10)
- Sizes: 32, 64
- Spacetime dims: 2, 3, 4
- Epsilons: 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2
- Seeds (10, disjoint from Phase 2G): 1900, 1916, 1923, 1939, 1953, 1973, 1981, 1995, 2003, 2020
- Optimizer seed: 1987 (anneal_limit=10, max_data=4, T₀=100.0, γ=0.9)
- Initial-energy floor for validity: 0.0001

## Sample bookkeeping

- Expected rows: 480
- Generated rows: 480
- Invalid (|initial_energy| below floor): 0
- Invalid (midpoint_dim degenerate): 24
- Rows used for diagnostic: 456
- Wall-clock: 67s

## New adimensional features

Per-row, computed once from the causal matrix:

```
dim_discrepancy_abs           = |mm_dim - midpoint_dim|
dim_discrepancy_rel_midpoint  = |mm_dim - midpoint_dim| / midpoint_dim
dim_ratio_mm_midpoint         = mm_dim / midpoint_dim
```

Rows where `midpoint_dim` is zero, NaN, or non-finite are excluded.

## Targets

```
relative_drift          = warmup_delta_energy / initial_energy
abs_relative_drift      = |relative_drift|
log_abs_relative_drift  = log1p(|relative_drift|)
```

## Per-curve optimizer-response (mean across seeds)

| n | target_dim | metric | ε=0.01 | ε=0.02 | ε=0.04 | ε=0.06 | ε=0.08 | ε=0.1 | ε=0.15 | ε=0.2 |
| ---: | :---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 32 | 2 | mean abs rel drift | 0.361 | 0.143 | 0.0444 | 0.00164 | 0.0461 | 0.0107 | 0.00205 | 0.0247 |
| 32 | 2 | mean log1p(abs)    | 0.304 | 0.132 | 0.0425 | 0.00163 | 0.0434 | 0.0105 | 0.00203 | 0.0237 |
| 32 | 3 | mean abs rel drift | 0.301 | 0.196 | 0.124 | 0.0447 | 0.0254 | 0.0326 | 0.0453 | 0.059 |
| 32 | 3 | mean log1p(abs)    | 0.252 | 0.175 | 0.112 | 0.0415 | 0.0246 | 0.0314 | 0.0436 | 0.056 |
| 32 | 4 | mean abs rel drift | 0.424 | 0.284 | 0.174 | 0.0995 | 0.101 | 0.0689 | 0.108 | 0.118 |
| 32 | 4 | mean log1p(abs)    | 0.344 | 0.245 | 0.157 | 0.0929 | 0.0921 | 0.0659 | 0.096 | 0.104 |
| 64 | 2 | mean abs rel drift | 0.046 | 0.0138 | 0 | 0 | 0 | 0 | 0 | 0 |
| 64 | 2 | mean log1p(abs)    | 0.0443 | 0.0135 | 0 | 0 | 0 | 0 | 0 | 0 |
| 64 | 3 | mean abs rel drift | 0.0996 | 0.0286 | 0.00491 | 0.000509 | 0 | 0.00409 | 0.00199 | 0.0122 |
| 64 | 3 | mean log1p(abs)    | 0.0934 | 0.028 | 0.00482 | 0.000507 | 0 | 0.00401 | 0.00197 | 0.0119 |
| 64 | 4 | mean abs rel drift | 0.157 | 0.0662 | 0.0158 | 0.00294 | 0.00438 | 0.00554 | 0.00977 | 0.0114 |
| 64 | 4 | mean log1p(abs)    | 0.141 | 0.0622 | 0.0156 | 0.00292 | 0.00434 | 0.00548 | 0.00937 | 0.011 |

## Monotonicity counts

- Curves (n × target_dim): 6
- Curves (quasi-)monotonic in mean abs_relative_drift vs ε: 2
- Curves (quasi-)monotonic in mean log1p(abs) vs ε: 3
- Curves monotonic in either: 3

Quasi-monotonic = at least 5 of 7 successive ε-differences move in the dominant direction.

Per-curve flags:

| n×target_dim | mono in abs | mono in log1p(abs) |
| --- | :---: | :---: |
| 32|2 | ✓ | ✓ |
| 32|3 | — | — |
| 32|4 | — | ✓ |
| 64|2 | — | — |
| 64|3 | ✓ | ✓ |
| 64|4 | — | — |

## Curve-shape diagnostic (post-hoc, exploratory)

This section is derived **after** observing the AMBIGUOUS verdict. The pre-registered monotonic criterion penalises V-shapes because a V-shaped curve is not monotone, even though a minimum followed by a clear rise is a potentially informative optimizer-response pattern.

**V-shape detection rule** (applied to mean `abs_relative_drift` vs ε):

- Interior minimum: 0 < argmin < len−1
- Fall fraction: (vals[0] − vals[argmin]) / vals[0] ≥ 30%
- Rise fraction: (vals[−1] − vals[argmin]) / fall ≥ 5%
- Left wing: ≥ 60% of consecutive diffs are negative
- Right wing: 0 negative diffs if n_pairs < 4; ≤ 1 if n_pairs ≥ 4

| n×target_dim | shape (abs_rd) | shape (log1p_abs_rd) | consensus shape |
| --- | --- | --- | --- |
| 32|2 | monotone_decay | monotone_decay | **monotone_decay** |
| 32|3 | v_shape | v_shape | **v_shape** |
| 32|4 | v_shape | v_shape | **v_shape** |
| 64|2 | noisy | noisy | **noisy** |
| 64|3 | monotone_decay | monotone_decay | **monotone_decay** |
| 64|4 | v_shape | v_shape | **v_shape** |

**Summary:**

- `verdict_original` = **AMBIGUOUS**
- `shape_diagnostic` = **PARTIAL_V_LIKE_AGGREGATE_MORPHOLOGY**
- `n_v_shape` = 3/6
- `v_shape_cells` = [(32,3), (32,4), (64,4)]

The pre-registered monotonic criterion yields **AMBIGUOUS**. A post-hoc curve-shape diagnostic detects V-shapes in 3/6 curves ((32,3), (32,4), (64,4)). These findings are consistent: the monotonic criterion does not recognise V-shapes by construction. This supports an exploratory (not confirmatory) Phase 4B.

## Physical interpretation (exploratory, N=6)

The following observations are based on six (n, target_dim) cells. With N=6 a clean empirical separation does not yet constitute a physical law; it motivates a hypothesis to test in Phase 4B.

**Observed pattern:**
V-like aggregate morphology (interior minimum or decline-then-rise) appears only when two conditions hold jointly:

```
V-shape candidate if:
    dim_discrepancy_rel_midpoint > θ   (θ ≈ 0.35–0.40, exploratory, not calibrated)
    and min_val > floor_tolerance       (curve does not saturate to zero)
```

- Cells with low `dim_discrepancy_rel` (< 0.16) show no V-shape regardless of other invariants — the embedding is well-conditioned and the optimizer is insensitive to the ε regime.
- `(64,3)` is the informative boundary case: its `dim_discrepancy_rel` and `ordering_fraction` are comparable to the V-shape cells, but its minimum reaches the numerical floor (`min_val = 0`). The potential V-signal is censored by floor saturation, not absent.
- `ordering_fraction` and `chain3_abundance` vary consistently with `dim_discrepancy_rel` across the six cells, but it is not yet possible to determine which invariant is the primary predictor and which is a correlated proxy.

**Phase 4B question priority:**

1. Does `dim_discrepancy_rel` separate V/non-V cells when the (n, d) grid is expanded? (survival test — must pass before any stronger claim)
2. Does `min_val ≈ 0` explain apparent false negatives such as `(64,3)`?
3. Do `ordering_fraction` and `chain3_abundance` contribute independent predictive information once conditioned on `dim_discrepancy_rel`?

Only if questions 1–3 survive should Phase 4B attempt symbolic regression of `rise_frac` or `ε_at_min` as continuous targets.

**Conservative conclusion:**

Phase 4A does not establish a morphological regime shift, but it identifies a clean exploratory pattern: V-like aggregate morphology occurs only in cells with large dimension-estimator discrepancy and without numerical floor saturation. The primary candidate control variable is `dim_discrepancy_rel_midpoint`; floor saturation is a secondary censoring mechanism. Phase 4B should test whether this separation survives a larger grid before attempting symbolic regression of `rise_frac` or `ε_at_min`.

## Full per-stratum table

| n | target_dim | ε | count | mean rel | std rel | mean abs | mean log1p(abs) | mean dim_disc_rel |
| ---: | :---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 32 | 2 | 0.01 | 10 | -0.361 | 0.1233 | 0.361 | 0.3042 | 0.135 |
| 32 | 2 | 0.02 | 10 | -0.1428 | 0.06062 | 0.1428 | 0.1321 | 0.135 |
| 32 | 2 | 0.04 | 10 | -0.04438 | 0.04519 | 0.04438 | 0.04251 | 0.135 |
| 32 | 2 | 0.06 | 10 | -0.001637 | 0.003198 | 0.001637 | 0.001631 | 0.135 |
| 32 | 2 | 0.08 | 10 | -0.04614 | 0.06209 | 0.04614 | 0.04341 | 0.135 |
| 32 | 2 | 0.1 | 10 | -0.01074 | 0.01703 | 0.01074 | 0.01054 | 0.135 |
| 32 | 2 | 0.15 | 10 | -0.002047 | 0.00614 | 0.002047 | 0.002026 | 0.135 |
| 32 | 2 | 0.2 | 10 | -0.02468 | 0.03904 | 0.02468 | 0.02367 | 0.135 |
| 32 | 3 | 0.01 | 10 | -0.3009 | 0.1943 | 0.3009 | 0.2521 | 0.5711 |
| 32 | 3 | 0.02 | 10 | -0.1965 | 0.1132 | 0.1965 | 0.1748 | 0.5711 |
| 32 | 3 | 0.04 | 10 | -0.1241 | 0.1177 | 0.1241 | 0.1117 | 0.5711 |
| 32 | 3 | 0.06 | 10 | -0.04471 | 0.07231 | 0.04471 | 0.04154 | 0.5711 |
| 32 | 3 | 0.08 | 10 | -0.02538 | 0.03004 | 0.02538 | 0.02464 | 0.5711 |
| 32 | 3 | 0.1 | 10 | -0.03255 | 0.03687 | 0.03255 | 0.03141 | 0.5711 |
| 32 | 3 | 0.15 | 10 | -0.04529 | 0.04015 | 0.04529 | 0.04356 | 0.5711 |
| 32 | 3 | 0.2 | 10 | -0.05897 | 0.05469 | 0.05897 | 0.05601 | 0.5711 |
| 32 | 4 | 0.01 | 8 | -0.4243 | 0.1965 | 0.4243 | 0.3435 | 1.51 |
| 32 | 4 | 0.02 | 8 | -0.2837 | 0.1279 | 0.2837 | 0.2447 | 1.51 |
| 32 | 4 | 0.04 | 8 | -0.1739 | 0.09811 | 0.1739 | 0.1568 | 1.51 |
| 32 | 4 | 0.06 | 8 | -0.09949 | 0.069 | 0.09949 | 0.09288 | 1.51 |
| 32 | 4 | 0.08 | 8 | -0.1014 | 0.1061 | 0.1014 | 0.09206 | 1.51 |
| 32 | 4 | 0.1 | 8 | -0.06894 | 0.04152 | 0.06894 | 0.06591 | 1.51 |
| 32 | 4 | 0.15 | 8 | -0.1076 | 0.1279 | 0.1076 | 0.09597 | 1.51 |
| 32 | 4 | 0.2 | 8 | -0.1183 | 0.143 | 0.1183 | 0.1043 | 1.51 |
| 64 | 2 | 0.01 | 10 | -0.04605 | 0.03893 | 0.04605 | 0.04433 | 0.07348 |
| 64 | 2 | 0.02 | 10 | -0.01375 | 0.02011 | 0.01375 | 0.01346 | 0.07348 |
| 64 | 2 | 0.04 | 10 | +0 | 0 | 0 | 0 | 0.07348 |
| 64 | 2 | 0.06 | 10 | +0 | 0 | 0 | 0 | 0.07348 |
| 64 | 2 | 0.08 | 10 | +0 | 0 | 0 | 0 | 0.07348 |
| 64 | 2 | 0.1 | 10 | +0 | 0 | 0 | 0 | 0.07348 |
| 64 | 2 | 0.15 | 10 | +0 | 0 | 0 | 0 | 0.07348 |
| 64 | 2 | 0.2 | 10 | +0 | 0 | 0 | 0 | 0.07348 |
| 64 | 3 | 0.01 | 10 | -0.09963 | 0.06225 | 0.09963 | 0.09338 | 0.1526 |
| 64 | 3 | 0.02 | 10 | -0.02862 | 0.02326 | 0.02862 | 0.02797 | 0.1526 |
| 64 | 3 | 0.04 | 10 | -0.004909 | 0.01284 | 0.004909 | 0.004817 | 0.1526 |
| 64 | 3 | 0.06 | 10 | -0.0005086 | 0.001526 | 0.0005086 | 0.0005073 | 0.1526 |
| 64 | 3 | 0.08 | 10 | +0 | 0 | 0 | 0 | 0.1526 |
| 64 | 3 | 0.1 | 10 | -0.004095 | 0.01228 | 0.004095 | 0.004013 | 0.1526 |
| 64 | 3 | 0.15 | 10 | -0.001987 | 0.00596 | 0.001987 | 0.001967 | 0.1526 |
| 64 | 3 | 0.2 | 10 | -0.01222 | 0.02157 | 0.01222 | 0.01193 | 0.1526 |
| 64 | 4 | 0.01 | 9 | -0.1569 | 0.1165 | 0.1569 | 0.1408 | 0.6799 |
| 64 | 4 | 0.02 | 9 | -0.06625 | 0.06668 | 0.06625 | 0.06225 | 0.6799 |
| 64 | 4 | 0.04 | 9 | -0.01583 | 0.01492 | 0.01583 | 0.0156 | 0.6799 |
| 64 | 4 | 0.06 | 9 | -0.002943 | 0.005708 | 0.002943 | 0.002923 | 0.6799 |
| 64 | 4 | 0.08 | 9 | -0.004378 | 0.007411 | 0.004378 | 0.004342 | 0.6799 |
| 64 | 4 | 0.1 | 9 | -0.00554 | 0.009924 | 0.00554 | 0.005477 | 0.6799 |
| 64 | 4 | 0.15 | 9 | -0.009774 | 0.02765 | 0.009774 | 0.009368 | 0.6799 |
| 64 | 4 | 0.2 | 9 | -0.01142 | 0.02563 | 0.01142 | 0.01105 | 0.6799 |

## Recommendation

Phase 4B is justified as an **exploratory** (not confirmatory) run. See the Physical interpretation section above for the question priority and the conditions under which symbolic regression of `rise_frac` / `ε_at_min` becomes appropriate.

Design notes:

- Include all (n, d) panels; d=2 serves as a negative control.
- Report V-shape cells ((32,3), (32,4), (64,4)) as a secondary analysis, not as the primary stratification.
- Do not restrict the run to V-shape cells only — that would condition on a post-hoc observation.
- Label all Phase 4B results as exploratory until the `dim_discrepancy_rel` separation survives a larger (n, d) grid.

## Reproducibility

Regenerate via `make regen-phase4a`.
Source: `tools/build_phase4a_epsilon_sweep.py`.

No PySR is invoked by this script. No data from Phase 2G or earlier phases is reused.
