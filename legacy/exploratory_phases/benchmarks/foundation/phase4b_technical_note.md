# Phase 4B/5 Technical Note: Optimizer-Response Morphology in Finite Causal Sets

## Scope

This is an internal methodological note. It does not claim a new physical law, and it does not establish a continuum Hauptvermutung result.

The purpose is to document a reproducible exploratory pipeline and its limitations. Phase 4B fixes the aggregate-level state of the Phase 4A survival probe, and Phase 5 closes that audit at seed level. The combined result is methodological, not a standalone physical claim.

## Definition of loss

The Phase 4B `loss` column is inherited from Phase 4A. In `tools/build_phase4b_survival_probe.py`, `build_rows()` calls `p4a._run_one(d, n, seed, epsilon, matrix, points)` and then `p4a._row_from_sim(...)`. The per-seed writer persists:

```text
loss = row["abs_relative_drift"]
```

In `tools/build_phase4a_epsilon_sweep.py`, `_row_from_sim(...)` defines:

```text
relative_drift = warmup_delta_energy / initial_energy
abs_relative_drift = |relative_drift|
```

Therefore:

```text
loss = |warmup_delta_energy / initial_energy|
```

The inputs are the causal matrix, sprinkled coordinates `points`, `target_dim`, `epsilon`, `seed`, the epsilon-noised initialization, and the guarded warmup behavior of `ConesSimulator`. `initial_energy` is the simulator energy after initializing the embedding from the sprinkled coordinates plus epsilon-scaled coordinate noise. `warmup_delta_energy` is the energy change during guarded warmup, before the later anneal result is used.

`loss` is an optimizer/embedding response diagnostic under this specific energy, initialization, epsilon, seed, and coordinate parametrization. It is not an intrinsic observable of the partial order, not a Lorentz-invariant distance to Minkowski space, not the dimensional discrepancy, and not an absolute physical quality score for the causet.

The visual audit roles `best`, `near_mean`, and `worst` rank seeds only by this Phase 4B `loss` within the same cell and epsilon. They are representatives under that metric, not absolute physical labels for good or bad causal sets.

## Phase 4B design

Phase 4B uses the `pilot` grid:

- Sizes: `n = 32, 48, 64`
- Target spacetime dimensions: `target_dim = 2, 3, 4`
- Epsilons: `0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2`
- Seeds: `1900, 1916, 1923, 1939, 1953, 1973, 1981, 1995, 2003, 2020`
- Exploratory threshold: `theta = 0.35`
- Floor tolerance: `1e-06`

Persisted outputs:

- `benchmarks/foundation/phase4b_survival_probe.csv`: one row per `(n, target_dim)` cell.
- `benchmarks/foundation/phase4b_survival_probe_per_epsilon.csv`: one row per `(n, target_dim, epsilon)`.
- `benchmarks/foundation/phase4b_survival_probe_per_seed.csv`: one row per `(n, target_dim, epsilon, seed)`.
- `benchmarks/foundation/phase4b_causet_visual_index.csv`: selected visual audit seeds and SVG paths.
- `benchmarks/foundation/phase4b_visual_link_audit.csv`: selected visual audit seeds with Hasse link counts.
- `benchmarks/foundation/phase4b_intrinsic_poset_audit.csv`: selected visual audit seeds with intrinsic poset observables.
- `benchmarks/foundation/phase4b_survival_probe.md`: generated exploratory report.

## Main result

Phase 4B finds a real but fragile optimizer-response morphology: interior minima and V-like curves appear in structured regions of the pilot grid, but the binary V/non-V classification and best/worst seed rankings do not reduce to a robust scalar intrinsic observable of the causal set. The result is therefore methodological rather than a standalone physical claim.

Phase 4B identifies a real but fragile aggregate optimizer-response morphology, while Phase 5 shows that seed-level morphology is largely censored by exact optimizer-floor saturation. The combined result is methodological: the observed V-like structure is reproducible at aggregate level but cannot yet be promoted to a seed-level or poset-intrinsic physical claim.

## Outcome

Global outcome: **MIXED**.

The `target_dim=2` controls remain clean in the sense that they are not promoted to strong positive evidence. Some V-like or interior-minimum behavior survives beyond Phase 4A, especially in structured regions of the pilot grid. However, counterexamples and borderline cases prevent a `PASS_EXPLORATORY_SURVIVAL` claim.

The current result should be read as partial survival of the Phase 4A morphology under stronger provenance, not as confirmation of a physical transition.

## Per-epsilon audit

The per-epsilon CSV was needed because the aggregate CSV records one row per cell and cannot show whether a counterexample comes from a clean curve, a marginal tail rise, floor effects, or alternating/noisy tail behavior.

The per-epsilon audit showed:

- `(48,3)` is a marginal V-like case. It is classified as `v_shape`, but its `rise_frac_margin` is small and the morphology is close to the threshold.
- `(48,4)` has an interior minimum and an alternating/noisy tail. It is not a clean monotone-decay story, but it also does not become a clean V-shape under the existing classifier.

This audit did not change thresholds, `curve_shape`, `survival_label`, or the global `MIXED` outcome.

## Per-seed audit

The per-seed CSV was needed because the causet depends on `(n, target_dim, seed)`, while `epsilon` selects an optimizer/noise row for that fixed causet. Without per-seed provenance, drawing a Hasse diagram for an arbitrary seed would be visually interesting but weakly traceable.

The visual audit therefore selected, at each priority cell's `epsilon_at_min`:

- `near_mean`: valid seed with `loss` closest to the per-epsilon mean.
- `best`: valid seed with minimum `loss`.
- `worst`: valid seed with maximum `loss`.

These roles are selection/provenance labels under the Phase 4B loss definition. They are not physical quality labels.

## Visual audit

The visual audit uses Hasse/projection SVGs generated from the transitive reduction links, not all transitive relations. This avoids drawing the full transitive closure.

The SVG projection uses only the `(t, x1)` plane. For `target_dim=3` and `target_dim=4`, the images are therefore shadows of the causet geometry rather than full spatial representations.

The visual audit supports seed-by-seed heterogeneity, but visual regularity does not track `loss` robustly. The diagrams are useful for audit and intuition, not for confirmatory evidence.

## Intrinsic poset audit

No robust scalar predictor is established: there is no stable monotonic relationship between loss and any single observable among ordering_fraction, n_links_Hasse, or chain3_abundance.

`dim_discrepancy_rel_midpoint` shows partial alignment in some cells, but it is not a complete explanation of the Phase 4B behavior. In the inspected `target_dim=4` cases, seed `2003` appears intrinsically sparse in the tabulated observables, not merely visually sparse in the SVG projection.

The intrinsic audit reinforces the methodological reading: the aggregate morphology is not explained by a single scalar observable currently persisted in the Phase 4B artifacts.

## Phase 5 seed-level audit

### Objective

Phase 5 asks whether the V-like and interior-minimum morphology observed in aggregate Phase 4B curves also appears in `loss(epsilon)` curves reconstructed for individual seeds.

### Method

Phase 5 reconstructs one curve per `(n, target_dim, seed)` directly from `benchmarks/foundation/phase4b_survival_probe_per_seed.csv` by grouping rows by seed, sorting by epsilon, and using valid `loss(epsilon)` values only. No new simulations are run, and no Phase 4B thresholds, labels, or outcomes are changed.

As in Phase 4B, `loss` is an optimizer-response quantity, not physical quality. The roles `best`, `near_mean`, and `worst` remain representatives under `loss`, not physical labels for individual causets.

### Outcome

Phase 5 global outcome: **INSUFFICIENT**.

This does not mean that Phase 4B failed. It means the current seed-level audit cannot confirm aggregate morphology cleanly because most seed curves are censored by the optimizer floor before a seed-level shape can be read out robustly.

### Floor-censoring result

The dominant seed-level effect is not a loose threshold choice. It is exact floor saturation.

- Phase 5 reconstructed seed-level `loss(epsilon)` curves from `phase4b_survival_probe_per_seed.csv`.
- `83/86` seed-curves with finite `min_loss` reached `min_loss == 0.0` exactly.
- The same `83/86` count persists for `min_loss <= 1e-15`, `1e-12`, `1e-9`, `1e-6`, and `1e-4`.
- Therefore the main censoring mechanism is exact optimizer-floor saturation, not sensitivity to the chosen `floor_tolerance`.

### Interpretation

Phase 4B outcome remains **MIXED** and Phase 5 outcome remains **INSUFFICIENT**.

Phase 4B found a real but fragile aggregate optimizer-response morphology. Phase 5 then showed that the corresponding seed-level morphology is largely censored by exact floor saturation. Therefore the aggregate V-shapes from Phase 4B should not be interpreted as direct evidence of V-shapes at the level of individual causets.

This also sharpens the reading of the aggregate counterexamples. They are not yet clean seed-level physical negatives; they are aggregate optimizer-response patterns observed under strong seed-level censoring.

### Limitation

Phase 5 does not establish that seed-level morphology is absent. It establishes that, in the current pilot grid and under the current optimizer-response metric, seed-level morphology is mostly unresolvable because the curves saturate at exact zero.

## Limitations

- The result is based on the pilot grid only.
- All observations are finite-`n` observations.
- `loss` is optimizer- and pipeline-specific.
- Phase 5 seed-level morphology is dominated by exact floor censoring in the current pilot grid.
- SVGs are `(t, x1)` projections and do not show the full geometry for `target_dim=3` or `target_dim=4`.
- Selected visual seeds are representative under `loss`, not physical quality labels.
- No BDG action is included.
- There is no full distributional audit over all seeds beyond the available per-seed CSV.
- No symbolic regression or PySR is used.

## Reproducibility

Existing commands:

```bash
make regen-phase4b
make regen-phase5
python3 -m unittest tests.test_phase4b_survival_probe -v
python3 -m unittest tests.test_phase4b_visualize_causets -v
python3 -m unittest tests.test_phase5_seed_curve_morphology -v
```

These commands regenerate the Phase 4B and Phase 5 reports from the current pipeline and run the Phase 4B visual and Phase 5 seed-level tests.

## Recommended next step

Document loss semantics and floor-censoring behavior before any Phase 4C, PySR, BDG action, or grid expansion.
