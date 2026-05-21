# Phase 4D — Robustness-vs-invariants audit (no-PySR)

**Status:** descriptive audit.  No PySR.  No new simulations.  No modification to `cones.py`, Phase 4A/4B/4C/5 CSVs, thresholds, or verdicts.

## Semantic caveats (read first)

- `loss = |warmup_delta_energy / initial_energy|` is an optimizer/embedding response diagnostic of the (causet × pipeline) pair.  It is **not** a physical observable of the partial order, not a Lorentz-invariant residual, not a manifoldness score, and not an embeddability witness.
- `iqr_loss_K`, `per_seed_label_distinct_shapes_K`, and `floor_saturated_fraction_K` are properties of (causet × ε × schedule × energy × optimizer family × K optimizer seeds).  They are **not** properties of the causet alone.
- A non-zero correlation between an order-theoretic invariant and a robustness target means: "this invariant coexists with seed instability under the current pipeline".  It does **not** establish an embedding rule, a manifoldness signature, or any pipeline-independent geometric claim.

## Objective

Phase 4C showed that the Phase 4B `MIXED` / Phase 5 `INSUFFICIENT` picture is `OPTIMIZER_SEED_LIMITED` (5/9 cells flip curve-shape label across K=3 optimizer seeds; per-cell IQR/loss ratio ≈ 0.92; floor saturation invariant under seed at 0.474 ≈ 0.473).  Phase 4D asks whether any order-theoretic invariant **coexists** with this residual variance or with the floor pathology, under the current pipeline.

## Inputs

- `phase4c_optimizer_seed_probe_per_run.csv`: per-run optimizer-seed multi-start from Phase 4C.
- `phase4c_optimizer_seed_probe_per_cell_epsilon.csv`: per-cell-epsilon aggregate from Phase 4C (label_stability_cell, mean/iqr/floor over K×seeds).
- Order-theoretic invariants: recomputed per (n, target_dim, causet_seed) via `p4a.compute_invariants` (identical to Phase 1/4A).

## Method

- K = 3 optimizer seeds: 1987, 1990, 1993
- Causet seeds (per cell): 10
- Per-seed level (N = 90): one row per (n, target_dim, causet_seed).  For each epsilon, compute IQR/mean/median/min/floor across K runs; then average over the 8 epsilons.  Shape per K optimizer seed (per causet) classified individually; the number of distinct shape labels across K is `per_seed_label_distinct_shapes_K` ∈ {1, 2, 3}.
- Per-cell level (N = 9): one row per (n, target_dim).  Invariants averaged across the 10 causet_seeds.  `label_stability_cell` and `curve_shape_per_optimizer_seed` propagated from Phase 4C.
- Correlation: Spearman ρ and Pearson r (pairwise complete on finite values).  Spearman ρ is the primary statistic.

## Targets

Robustness targets (per-seed level):

- `per_seed_iqr_loss_K_mean_eps` — seed-dispersion robustness
- `per_seed_label_distinct_shapes_K` — morphology robustness (1, 2 or 3)
- `per_seed_floor_saturated_fraction_K_mean_eps` — floor pathology

`min_over_K` is explicitly excluded as a target: trivially decreases with K (optimistic lucky-seed selection).

## Verdict

**ORDER_THEORETIC_CORRELATE_DETECTED**

Decision rule (per-seed Spearman ρ, N = 90):

- `ORDER_THEORETIC_CORRELATE_DETECTED` if max |ρ| ≥ 0.6 AND a second invariant reaches |ρ| ≥ 0.6 against the same target with the same sign.
- `WEAK_CORRELATE` if 0.3 ≤ max |ρ| < 0.6.
- `NO_ROBUST_ORDER_THEORETIC_CORRELATE` if max |ρ| < 0.3.

Verdict inputs:

- Top invariant: `relation_count`
- Top target: `per_seed_floor_saturated_fraction_K_mean_eps`
- Top |ρ_spearman|: 0.9405720378
- Top signed ρ_spearman: 0.9405720378
- Pairs with |ρ| ≥ 0.3: 42
- Pairs with |ρ| ≥ 0.6: 25
- Detected pair: `('relation_count', 'chain2_count')`

## Correlation matrix (per-seed level, signed Spearman ρ)

| invariant | ρ vs per_seed_iqr_loss_K_mean_eps | ρ vs per_seed_label_distinct_shapes_K | ρ vs per_seed_floor_saturated_fraction_K_mean_eps |
| --- | ---: | ---: | ---: |
| `mm_dim` | 0.5378014299 | 0.4462158695 | -0.6981440588 |
| `midpoint_dim` | -0.2743752927 | 0.02888267893 | 0.1616886578 |
| `abs_discrepancy_mm_midpoint` | 0.6075097882 | 0.435905438 | -0.7593890765 |
| `dim_discrepancy_abs` | 0.6075097882 | 0.435905438 | -0.7593890765 |
| `dim_discrepancy_rel_midpoint` | 0.6142648238 | 0.4194790218 | -0.746921015 |
| `dim_ratio_mm_midpoint` | 0.6243407708 | 0.3963555997 | -0.7490384113 |
| `ordering_fraction` | -0.5378014299 | -0.4462158695 | 0.6981440588 |
| `chain2_count` | -0.8406560618 | -0.5537940859 | 0.9405720378 |
| `chain3_count` | -0.7679396791 | -0.5456042698 | 0.9105405838 |
| `chain3_abundance` | -0.5427156775 | -0.4545274859 | 0.7062581297 |
| `chain4_count` | -0.7237279512 | -0.5327763105 | 0.8773179646 |
| `link_count` | -0.7710013241 | -0.3723655419 | 0.7503845083 |
| `link_density` | -0.7120516796 | -0.3733931904 | 0.7195354252 |
| `height` | -0.6850714637 | -0.4991042736 | 0.8362909837 |
| `relation_count` | -0.8406560618 | -0.5537940859 | 0.9405720378 |

## n-control (interpretive layer)

The n-control is interpretive auditing, not mechanical reclassification. The Phase 4D verdict (`ORDER_THEORETIC_CORRELATE_DETECTED` or otherwise) is computed from raw Spearman as before. The partial-n and stratified-by-n values are reported here so the reader can judge whether the raw correlation reflects n-scaling rather than an order-theoretic correlate.

### Target: `per_seed_iqr_loss_K_mean_eps`

| invariant | raw ρ | partial ρ (ctrl n) | min_abs stratified | ρ_n32 | ρ_n48 | ρ_n64 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `mm_dim` | 0.5378014299 | 0.701682152 | 0.6269885454 | 0.7186561397 | 0.7174638487 | 0.6269885454 |
| `midpoint_dim` | -0.2743752927 | -0.02368030145 | 0.06094260261 | -0.1914694548 | -0.06094260261 | 0.1976586943 |
| `abs_discrepancy_mm_midpoint` | 0.6075097882 | 0.6123267924 | 0.3625615764 | 0.6967706623 | 0.645320197 | 0.3625615764 |
| `dim_discrepancy_abs` | 0.6075097882 | 0.6123267924 | 0.3625615764 | 0.6967706623 | 0.645320197 | 0.3625615764 |
| `dim_discrepancy_rel_midpoint` | 0.6142648238 | 0.5928798851 | 0.3246305419 | 0.6787082649 | 0.5975369458 | 0.3246305419 |
| `dim_ratio_mm_midpoint` | 0.6243407708 | 0.6253882864 | 0.3802955665 | 0.7000547345 | 0.6379310345 | 0.3802955665 |
| `ordering_fraction` | -0.5378014299 | -0.701682152 | 0.6269885454 | -0.7186561397 | -0.7174638487 | -0.6269885454 |
| `chain2_count` | -0.8406560618 | -0.7292487733 | 0.6269885454 | -0.7186561397 | -0.7174638487 | -0.6269885454 |
| `chain3_count` | -0.7679396791 | -0.7117556524 | 0.6266963293 | -0.6870619689 | -0.7104238558 | -0.6266963293 |
| `chain3_abundance` | -0.5427156775 | -0.6994328807 | 0.6266963293 | -0.6870619689 | -0.7104238558 | -0.6266963293 |
| `chain4_count` | -0.7237279512 | -0.6948500524 | 0.64434309 | -0.6999038163 | -0.686700162 | -0.64434309 |
| `link_count` | -0.7710013241 | -0.5429525103 | 0.1974182257 | -0.5966596127 | -0.4541407958 | -0.1974182257 |
| `link_density` | -0.7120516796 | -0.510079032 | 0.1974182257 | -0.5966596127 | -0.4541407958 | -0.1974182257 |
| `height` | -0.6850714637 | -0.6821886165 | 0.6317516704 | -0.7202922479 | -0.6317516704 | -0.6798894831 |
| `relation_count` | -0.8406560618 | -0.7292487733 | 0.6269885454 | -0.7186561397 | -0.7174638487 | -0.6269885454 |

> Top raw correlate is an extensive size-like invariant. The raw correlation is expected to track `n` independently of any order-theoretic content. Refer to the partial-n and stratified-n columns to judge whether the correlation survives after the finite-size effect is removed.

### Target: `per_seed_label_distinct_shapes_K`

| invariant | raw ρ | partial ρ (ctrl n) | min_abs stratified | ρ_n32 | ρ_n48 | ρ_n64 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `mm_dim` | 0.4462158695 | 0.4619409296 | 0.247882157 | 0.247882157 | 0.5784975848 | 0.5554829749 |
| `midpoint_dim` | 0.02888267893 | 0.1566764839 | 0.0002801120448 | 0.2529010654 | 0.0002801120448 | 0.3080137613 |
| `abs_discrepancy_mm_midpoint` | 0.435905438 | 0.3920300488 | 0.2601329909 | 0.2824625763 | 0.5477925762 | 0.2601329909 |
| `dim_discrepancy_abs` | 0.435905438 | 0.3920300488 | 0.2601329909 | 0.2824625763 | 0.5477925762 | 0.2601329909 |
| `dim_discrepancy_rel_midpoint` | 0.4194790218 | 0.3667226341 | 0.238982382 | 0.238982382 | 0.5243037674 | 0.2427907915 |
| `dim_ratio_mm_midpoint` | 0.3963555997 | 0.3464445656 | 0.1977238034 | 0.1977238034 | 0.5164741645 | 0.2601329909 |
| `ordering_fraction` | -0.4462158695 | -0.4619409296 | 0.247882157 | -0.247882157 | -0.5784975848 | -0.5554829749 |
| `chain2_count` | -0.5537940859 | -0.5007820079 | 0.247882157 | -0.247882157 | -0.5784975848 | -0.5554829749 |
| `chain3_count` | -0.5456042698 | -0.4877566527 | 0.2588991418 | -0.2588991418 | -0.6192736897 | -0.5234544957 |
| `chain3_abundance` | -0.4545274859 | -0.4676398916 | 0.2588991418 | -0.2588991418 | -0.6192736897 | -0.5234544957 |
| `chain4_count` | -0.5327763105 | -0.4807821995 | 0.2310421584 | -0.2310421584 | -0.6232176927 | -0.475557367 |
| `link_count` | -0.3723655419 | -0.2673069696 | 0.1764485288 | -0.1764485288 | -0.2296054588 | -0.2078759766 |
| `link_density` | -0.3733931904 | -0.2585054248 | 0.1764485288 | -0.1764485288 | -0.2296054588 | -0.2078759766 |
| `height` | -0.4991042736 | -0.4527308965 | 0.2600753589 | -0.2600753589 | -0.5299915719 | -0.4740435825 |
| `relation_count` | -0.5537940859 | -0.5007820079 | 0.247882157 | -0.247882157 | -0.5784975848 | -0.5554829749 |

> Top raw correlate is an extensive size-like invariant. The raw correlation is expected to track `n` independently of any order-theoretic content. Refer to the partial-n and stratified-n columns to judge whether the correlation survives after the finite-size effect is removed.

### Target: `per_seed_floor_saturated_fraction_K_mean_eps`

| invariant | raw ρ | partial ρ (ctrl n) | min_abs stratified | ρ_n32 | ρ_n48 | ρ_n64 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `mm_dim` | -0.6981440588 | -0.8770419631 | 0.7967446317 | -0.8959108044 | -0.9443904166 | -0.7967446317 |
| `midpoint_dim` | 0.1616886578 | -0.11331165 | 0.0890997568 | 0.0890997568 | -0.2362467711 | -0.2943680902 |
| `abs_discrepancy_mm_midpoint` | -0.7593890765 | -0.8011535016 | 0.6341661868 | -0.877051793 | -0.8438868381 | -0.6341661868 |
| `dim_discrepancy_abs` | -0.7593890765 | -0.8011535016 | 0.6341661868 | -0.877051793 | -0.8438868381 | -0.6341661868 |
| `dim_discrepancy_rel_midpoint` | -0.746921015 | -0.7610606074 | 0.5903537687 | -0.8284647249 | -0.7871176153 | -0.5903537687 |
| `dim_ratio_mm_midpoint` | -0.7490384113 | -0.7786943302 | 0.6079282415 | -0.8174222094 | -0.8209323263 | -0.6079282415 |
| `ordering_fraction` | 0.6981440588 | 0.8770419631 | 0.7967446317 | 0.8959108044 | 0.9443904166 | 0.7967446317 |
| `chain2_count` | 0.9405720378 | 0.9030016703 | 0.7967446317 | 0.8959108044 | 0.9443904166 | 0.7967446317 |
| `chain3_count` | 0.9105405838 | 0.9073105038 | 0.808052763 | 0.9047645639 | 0.9112809854 | 0.808052763 |
| `chain3_abundance` | 0.7062581297 | 0.8795065169 | 0.808052763 | 0.9047645639 | 0.9112809854 | 0.808052763 |
| `chain4_count` | 0.8773179646 | 0.8941353212 | 0.8119420047 | 0.9056861281 | 0.8760043525 | 0.8119420047 |
| `link_count` | 0.7503845083 | 0.5615020511 | 0.2612361863 | 0.6128682754 | 0.5574889887 | 0.2612361863 |
| `link_density` | 0.7195354252 | 0.5417292858 | 0.2612361863 | 0.6128682754 | 0.5574889887 | 0.2612361863 |
| `height` | 0.8362909837 | 0.8704864658 | 0.8075580758 | 0.9165848899 | 0.8491943683 | 0.8075580758 |
| `relation_count` | 0.9405720378 | 0.9030016703 | 0.7967446317 | 0.8959108044 | 0.9443904166 | 0.7967446317 |

> Top raw correlate is an extensive size-like invariant. The raw correlation is expected to track `n` independently of any order-theoretic content. Refer to the partial-n and stratified-n columns to judge whether the correlation survives after the finite-size effect is removed.

## Correlation matrix (per-cell level, signed Spearman ρ)

No n-control is computed at per-cell level (N=9 is too small for partial or stratified statistics; per-cell is reported as descriptive sanity only).

Per-cell N = 9; treat as sanity / consistency check.

| invariant | ρ vs cell_mean_iqr_loss_K_over_eps | ρ vs label_stability_cell | ρ vs cell_mean_floor_saturated_fraction_K_over_eps |
| --- | ---: | ---: | ---: |
| `mm_dim` | 0.45 | 0.1732050808 | -0.8 |
| `midpoint_dim` | -0.55 | -0.5196152423 | 0.2166666667 |
| `abs_discrepancy_mm_midpoint` | 0.6 | 0.08660254038 | -0.8833333333 |
| `dim_discrepancy_abs` | 0.6 | 0.08660254038 | -0.8833333333 |
| `dim_discrepancy_rel_midpoint` | 0.6 | 0.08660254038 | -0.8833333333 |
| `dim_ratio_mm_midpoint` | 0.6 | 0.08660254038 | -0.8833333333 |
| `ordering_fraction` | -0.4166666667 | -0.08660254038 | 0.7666666667 |
| `chain2_count` | -0.8666666667 | -0.08660254038 | 1 |
| `chain3_count` | -0.8 | 0 | 0.9666666667 |
| `chain3_abundance` | -0.3166666667 | -0.08660254038 | 0.6833333333 |
| `chain4_count` | -0.6666666667 | 0 | 0.9 |
| `link_count` | -0.95 | -0.08660254038 | 0.7666666667 |
| `link_density` | -0.9333333333 | -0.2598076211 | 0.85 |
| `height` | -0.6666666667 | 0 | 0.9 |
| `relation_count` | -0.8666666667 | -0.08660254038 | 1 |

## Cells with label_stability_cell = 0 (Phase 4C)

| n | target_dim | curve_shape_per_optimizer_seed | dim_disc_rel | ordering_fraction | chain3_abundance | link_density | height | mm_dim |
| ---: | :---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 32 | 3 | `v_shape|monotone_decay|v_shape` | 0.5711174782 | 0.2177419355 | 0.02346774194 | 2.03125 | 4.8 | 3.121316808 |
| 48 | 2 | `monotone_decay|noisy|noisy` | 0.06298588276 | 0.5086879433 | 0.1704151249 | 2.658333333 | 10.8 | 1.978913857 |
| 48 | 3 | `v_shape|v_shape|monotone_decay` | 0.2218230296 | 0.2283687943 | 0.02376850139 | 2.872916667 | 5.5 | 3.03556888 |
| 48 | 4 | `noisy|monotone_decay|noisy` | 0.9226915362 | 0.1148049645 | 0.003763876041 | 2.052083333 | 4.1 | 3.892211667 |
| 64 | 3 | `monotone_decay|v_shape|noisy` | 0.1526208251 | 0.2213789683 | 0.02295986943 | 3.259375 | 6.4 | 3.072535506 |

## Interpretation

A non-trivial coexistence pattern was detected between at least two order-theoretic invariants and one robustness target, with consistent sign.  This is descriptive: it does NOT establish a physical claim, an embedding rule, or a manifoldness signature.  It does establish that the optimizer-seed instability observed in Phase 4C has structural coexistence with order-theoretic features at this pipeline scale.  The natural follow-up is Phase 4E or a re-scoped Phase 3G that uses the aggregated target with an explicit caveat on the multiple-comparisons risk at N = 90.

## Scope

- No PySR.  No new optimizer.  No reheating.  No parallel tempering.
- No modification to Phase 4A/4B/4C/5 CSVs, thresholds, labels, or verdicts.
- No physical claim, no manifoldness claim, no embeddability claim.
- Targets are properties of (causet × pipeline), not of the causet alone.

## Reproducibility

Regenerate via `make regen-phase4d`.
Source: `tools/build_phase4d_robustness_audit.py`.
