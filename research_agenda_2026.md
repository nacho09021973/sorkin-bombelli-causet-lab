# Research Agenda 2026

This note frames the revived Bombelli 1987 annealing program as a research instrument rather than only a historical port.

Current empirical results are summarized in [`results_note_2026.md`](/home/adnac/sorkin/results_note_2026.md).

## Central Question

Which finite causal sets are naturally embeddable into low-dimensional Minkowski space, and what controls the transition between easy, hard, and apparently non-embeddable cases?

The old program asks whether one causet can be fitted. The modern tool should ask how whole families behave.

## Working Hypotheses

- Embeddability is not controlled only by `n`; density, dimension, and order structure matter.
- The annealing landscape has basins of attraction that can be mapped statistically.
- Schedule sensitivity is itself physical information if it correlates with causal-set structure.
- Failure to embed cleanly may separate numerical difficulty from genuinely non-manifoldlike order.

## Five-Step Program

1. Establish canonical benchmarks.

Use a small benchmark family instead of isolated examples:

- easy: small sprinkled causets where energy often reaches zero
- medium: cases like `tesis_like_12.in` where schedule choice matters strongly
- hard: larger or denser causets where most seeds fail
- controls: chains, antichains, layered posets, and random non-manifoldlike orders

The goal is not to win on one example. The goal is to make failure modes comparable.

2. Build phase maps.

For each benchmark family, scan:

- number of elements `n`
- target embedding dimension `dim`
- relation density
- initial temperature
- cooling factor
- random seed

Primary output:

- success probability
- mean and median final energy
- variance across seeds
- best seed per region
- heatmaps over dimension and schedule

This turns the old annealing run into a map of where the method works.

3. Study structure versus difficulty.

For every input causet, compute structural descriptors before annealing:

- number of relations
- relation density
- height
- width
- layer profile
- number of links in the transitive reduction
- interval-size distribution

Then compare those descriptors against:

- final energy
- success rate
- number of annealing blocks
- sensitivity to seed
- sensitivity to schedule

The scientific target is to identify which order-theoretic features predict geometric embeddability.

4. Compare energy models.

Keep the thesis energy as the historical baseline, then test variants:

- original Bombelli energy
- regularized radius penalties
- penalties that weight links differently from transitive relations
- interval-aware penalties
- multiobjective versions that separate false positives from false negatives

The question is not only whether the old energy works. It is which energy best detects manifoldlike causal structure.

5. Separate algorithmic failure from physical failure.

For hard cases, rerun with stronger methods:

- more seeds
- broader schedules
- local refinements near the best embedding
- higher target dimension
- alternative energy functions

If all methods fail, treat the causet as a candidate non-manifoldlike order.
If only the original schedule fails, the lesson is algorithmic rather than physical.

## Near-Term Experiments

Start with experiments that the current code can already support:

1. Run ensemble scans on `tesis_like_6.in` and `tesis_like_12.in`.
2. Generate sprinkled benchmarks for `n = 16, 24, 32, 48, 64`.
3. For each `n`, scan `dim = 1, 2, 3, 4`.
4. For each `(n, dim)`, estimate success probability over many seeds.
5. Produce a table of the largest `n` where the success rate remains above `50%`.

That last number is a clean first answer to the historical question: how far can the revived experiment go now?

## Metrics To Record

Every experiment should preserve:

- benchmark name
- generation seed
- number of elements
- target dimension
- number of relations
- number of links
- schedule parameters
- annealing seed
- initial energy
- warmup energy
- final energy
- success threshold
- runtime

Without these fields, the run is only anecdotal.

## Interpretation Rules

- Energy near zero means a good embedding was found, not that the causet is uniquely geometric.
- High final energy after one run is not evidence of non-embeddability.
- Repeated high energy across many schedules and seeds is evidence of difficulty.
- If higher dimension fixes the problem, the original target dimension was too restrictive.
- If schedule tuning fixes the problem, the original historical run was computationally limited.

## First Milestone

The first serious milestone is a phase diagram:

- x-axis: number of elements `n`
- y-axis: target dimension `dim`
- cell value: success probability

A second layer should show median final energy.

This is the experiment that Bombelli and Sorkin could not practically do in 1987: not a single optimized causet, but a statistical map of embeddability.

The first implementation is [`phase_diagram.py`](/home/adnac/sorkin/phase_diagram.py).

Example:

```bash
python3 phase_diagram.py --n-values 6,12 --dim-min 1 --dim-max 4 \
  --seed-start 1959 --seed-count 4 \
  --gpu-first --backend auto \
  --runs-csv /tmp/phase_runs.csv \
  --summary-csv /tmp/phase_summary.csv \
  --report-md /tmp/phase.md \
  --heatmap-svg /tmp/phase.svg
```

The command prints progress by default so long runs show which `(n, dim, seed)` is being evaluated.
It also writes partial outputs after each completed cell. For frontier scans, set `--max-run-seconds` so one pathological run cannot dominate the whole experiment.
Use `--fast-frontier` for exploratory runs near the computational boundary; it reduces the annealing budget and stops cells early when the available evidence is already decisive.

## Current Tools

- [`cones.py`](/home/adnac/sorkin/cones.py): revived annealing program
- [`phase_diagram.py`](/home/adnac/sorkin/phase_diagram.py): n-vs-dim embeddability maps
- [`ensemble_scan.py`](/home/adnac/sorkin/ensemble_scan.py): seed and schedule ensembles
- [`dimension_sweep.py`](/home/adnac/sorkin/dimension_sweep.py): dimension comparison
- [`schedule_sweep.py`](/home/adnac/sorkin/schedule_sweep.py): schedule comparison
- [`analyze_sweep.py`](/home/adnac/sorkin/analyze_sweep.py): sweep summaries
- [`causet_invariants.py`](/home/adnac/sorkin/causet_invariants.py): order-theoretic invariants (v2)
- [`validation_suite.py`](/home/adnac/sorkin/validation_suite.py): sprinkle-then-recover protocol (v2)
- [`results_note_2026.md`](/home/adnac/sorkin/results_note_2026.md): current empirical findings
- [`benchmarks/README.md`](/home/adnac/sorkin/benchmarks/README.md): benchmark documentation
- [`benchmarks/foundation/README.md`](/home/adnac/sorkin/benchmarks/foundation/README.md): foundation benchmark (v2)

## v2: Diagnostic Foundation

Before scaling the original annealing program to larger ``n`` or
substituting it with new optimization strategies, the v2 line of
work is laying down a diagnostic foundation that lets us tell
*algorithmic* failure from *physical* non-embeddability. Without
this layer, every failure of the optimizer is uninterpretable.

The foundation has four components.

### Order-theoretic invariants

[`causet_invariants.py`](/home/adnac/sorkin/causet_invariants.py)
provides a battery of structural descriptors that depend only on the
causal matrix, not on any embedding or random seed. They include the
ordering fraction, height, antichain profile, chain counts, link
count, and the **Myrheim-Meyer dimension** (Myrheim 1978, Meyer
1988), obtained by numerically inverting the closed-form expected
ordering fraction for a uniform Poisson sprinkling of a Minkowski
diamond,

```
f(d) = Gamma(d + 1) Gamma(d / 2) / (2 Gamma(3 d / 2)).
```

When the annealing optimizer fails on a causet, the Myrheim-Meyer
dimension is an *independent* witness: if it indicates a
well-defined ``d`` near the target embedding dimension, the failure
is most likely algorithmic; if it diverges or sits far from the
target, the causet is plausibly not manifoldlike.

### Canonical Minkowski sprinkler

The legacy ``cones.generate_sprinkled_causet`` samples the unit
causal diamond *uniformly* only for ``d_spacetime = 2``. For
``d_spacetime >= 3`` it distributes the time coordinate uniformly on
``[0, 1]``, whereas the canonical uniform measure on the diamond
has marginal proportional to the spatial cross-section volume,
``p(t) ~ min(t, 1 - t)^{d_spacetime - 1}``.

[`validation_suite.sprinkle_minkowski_diamond`](/home/adnac/sorkin/validation_suite.py)
adds a canonical sprinkler that samples the diamond exactly
uniformly in any dimension by rejection against the bounding box
``[0, 1] x B^{d_spatial}``. This is the function used by every
diagnostic and benchmark below; the legacy sprinkler is kept
unchanged to preserve reproducibility of the historical benchmarks.

### Sprinkle-then-recover validation

[`validation_suite.py`](/home/adnac/sorkin/validation_suite.py)
implements the validation protocol:

1. Sprinkle a known Minkowski diamond.
2. Run the embedding optimizer on the resulting causet, exposing
   only the causal matrix.
3. Compare the recovered embedding to the truth using the
   **Lorentz-invariant interval residual**

   ```
   RMSE_{ij} (s^2_recovered[i,j] - s^2_truth[i,j]),
   ```

   where ``s^2 = -dt^2 + |dx|^2`` is the squared Minkowski interval.
   Two embeddings related by a proper Poincare transform produce
   identical interval matrices, so this residual is the natural
   embedding-quality metric and avoids the explicit Lorentzian
   Procrustes problem.
4. Independently compute the Myrheim-Meyer dimension from the
   causal matrix as an embedding-free sanity check.
5. Compute the Bombelli energy at the *ground-truth coordinates*
   as well as at the recovered configuration. If the optimizer's
   final energy is far above the ground-truth energy, the
   optimizer is stuck; if both are similar, the energy function
   itself is the bottleneck.

### Foundation benchmark

[`benchmarks/foundation/`](/home/adnac/sorkin/benchmarks/foundation)
contains the frozen *protocol* of the v2 benchmark: a 45-cell grid
of canonical sprinklings at ``d_spacetime in {2, 3, 4}``, ``n in
{16, 32, 64}``, ``seed in {1959, 1962, 1987, 2009, 2026}``, together
with the precomputed order-theoretic invariants for each cell. The
matrices themselves are *not* stored as bytes; they are regenerated
on demand from the deterministic sprinkler, and a regression test
verifies that the recomputed invariants exactly match the frozen
JSON. Any change to the sprinkler or the invariants is therefore
caught immediately.

Regenerate the benchmark with

```
python3 tools/build_foundation_benchmarks.py
```

and verify integrity with

```
python3 -m unittest tests.test_foundation_benchmarks
```

### Invariance tests

[`tests/test_invariances.py`](/home/adnac/sorkin/tests/test_invariances.py)
checks the structural symmetries the diagnostic foundation must
respect: relabeling the events of a causet (composed with a fresh
topological sort) does not change any invariant, and a composite
Poincare transform (translation + boost + rotation + uniform
rescale) leaves the interval residual at machine zero. Existing
optimizer-seed reproducibility is already covered by the original
regression tests on ``cones.py``.

### Next steps

The v2 foundation is now in place. The follow-on work that it
unblocks, in order:

1. Run the foundation benchmark through ``cones.ConesSimulator``
   and record, per cell, the ratio
   ``final_energy / truth_energy`` and the interval residual.
   This is the first calibrated map of where the historical
   algorithm reaches the global minimum and where it gets stuck.
2. Build ``embedding_lab.py``: a harness that runs multiple
   engines (legacy annealing, smart initializations, basin
   hopping, parallel tempering) on the same foundation grid and
   scores them on the metrics above.
3. Investigate which order-theoretic features correlate with
   embedding difficulty by joining the invariant fingerprints
   against the recovery metrics.

## Phase 1: order-theoretic negative controls

A minimum Phase 1 has been completed before any embedding,
energy, or annealing is invoked. The physical question is

> Can internal invariants of the partial order alone separate
> manifoldlike sprinklings from non-manifoldlike controls?

The artifacts:

- **Kleitman-Rothschild control generator** in
  [`validation_suite.generate_kleitman_rothschild`](/home/adnac/sorkin/validation_suite.py).
  Following Kleitman and Rothschild (Trans. AMS 205, 1975), the
  typical labeled poset on ``n`` elements has three antichain
  levels with sizes near ``(n/4, n/2, n/4)``. This generator
  samples such a structure with deterministic levels and
  Bernoulli inter-level edges, then takes transitive closure. It
  is a **diagnostic control**, not a physical model.
- **Second order-theoretic estimator** in
  [`causet_invariants.midpoint_scaling_dimension`](/home/adnac/sorkin/causet_invariants.py),
  implementing Meyer's midpoint interval scaler. It is
  algorithmically independent of the Myrheim-Meyer ordering
  fraction: it picks the largest causal interval, finds the
  order-theoretic midpoint by maximizing the smaller
  sub-interval cardinality, and returns the log-ratio of full
  to balanced sub-interval. Agreement with MM is one signal of
  manifoldlike behavior; disagreement is informative.
- **Atlas comparator** in
  [`tools/build_phase1_atlas.py`](/home/adnac/sorkin/tools/build_phase1_atlas.py).
  For each foundation cell ``(n, seed)`` and each ``d_spacetime
  in {2, 3, 4}`` it builds a Minkowski-diamond sprinkling and a
  size-matched KR control, computes both dimension estimators,
  and writes a CSV plus markdown summary at
  ``benchmarks/foundation/phase1_atlas.{csv,md}``.

Regenerate via ``make regen-phase1``. The regression test
``tests/test_phase1_atlas.py`` verifies that recomputation
matches the recorded artifact byte for byte.

This phase does **not** evaluate any embedding algorithm. The
question it answers is upstream of optimization: are the
invariants we already know expressive enough to recognize
non-manifoldlike structure on their own?

## Phase 1B: finite-size scaling of the order-theoretic atlas

Phase 1 fixed ``n in {16, 32, 64}`` and left open whether the
diagnostic improves systematically with ``n``. Phase 1B answers
that question by extending the same atlas across
``n in {32, 64, 128, 256}`` while leaving the rest of the
protocol untouched. No embedding, energy, or optimizer is
invoked anywhere in this phase.

The single artifact is
[`benchmarks/foundation/phase1b_scaling_atlas.{csv,md}`](/home/adnac/sorkin/benchmarks/foundation/phase1b_scaling_atlas.md),
regenerable via ``make regen-phase1b``. It uses the same five
seeds as Phase 1 (so every Phase 1 cell at ``n in {32, 64}``
appears as a subset of the Phase 1B grid), the same Minkowski
sprinklings and Kleitman-Rothschild controls, and the same two
estimators.

[`tools/build_phase1b_scaling_atlas.py`](/home/adnac/sorkin/tools/build_phase1b_scaling_atlas.py)
reuses the Phase 1 helpers ``_estimate_dimensions``,
``_discrepancy``, and ``_format_field`` verbatim, so the two
atlases share their numerical conventions and a divergence
between them would point unambiguously at the new size axis.

The regression test ``tests/test_phase1b_scaling_atlas.py``
enforces schema, full ``(family, d, n, seed)`` coverage, and a
finiteness check on every cell — so a silent failure of
generator or estimator at any cell is caught rather than
absorbed into an ensemble mean.

Empirical interpretation of the resulting table lives in the
*Phase 1B* section of
[`results_note_2026.md`](/home/adnac/sorkin/results_note_2026.md).
The headline finding — that the *finite-size trajectory* of the
inter-estimator discrepancy ``|disc|`` separates manifoldlike
from KR cells with opposite sign of ``d|disc|/dn`` across the
entire grid — is the strongest single diagnostic the
order-theoretic foundation produces so far, but it is still
limited to a single non-manifoldlike control family.

## Phase 1C: second non-manifoldlike control

Phase 1C tests whether the Phase 1B signature survives a
structurally different negative control by adding suspended
corona/crown posets beside Kleitman-Rothschild controls. The
artifact is
[`benchmarks/foundation/phase1c_scaling_atlas.{csv,md}`](/home/adnac/sorkin/benchmarks/foundation/phase1c_scaling_atlas.md),
regenerable via ``make regen-phase1c``. It keeps the Phase 1B
size grid, seeds, Minkowski sprinklings, and dimension
estimators, and adds one new control family generated by
[`validation_suite.generate_corona_poset`](/home/adnac/sorkin/validation_suite.py).

[`tools/build_phase1c_scaling_atlas.py`](/home/adnac/sorkin/tools/build_phase1c_scaling_atlas.py)
reuses the Phase 1/1B numerical helpers, and the regression
test ``tests/test_phase1c_scaling_atlas.py`` enforces schema,
full grid coverage, and finite ``mm_dim``/``midpoint_dim``
values for every cell.

The empirical result is that corona controls show the same
qualitative non-manifoldlike scaling as KR: ``mean_mm`` is
comparatively flat within the family, while ``mean_midpoint``
grows with finite-size scale and the inter-estimator
discrepancy increases rather than shrinks. This does not yet
prove per-causet classification, but it removes the most direct
"KR-only artifact" objection to Phase 1B.

## Phase 1D: third order-theoretic invariant

Phase 1D adds one independent structural observable to the
Phase 1C grid: 3-chain abundance. The raw count ``C3`` is the
number of triples ``i prec j prec k``; the normalized observable
recorded in the atlas is ``C3 / binom(n, 3)``. This is a
minimal chain-abundance statistic, not a calibrated effective
dimension, because no analytic inversion formula is being
assumed.

The artifact is
[`benchmarks/foundation/phase1d_structural_atlas.{csv,md}`](/home/adnac/sorkin/benchmarks/foundation/phase1d_structural_atlas.md),
regenerable via ``make regen-phase1d``. The CSV records
``chain2_count``, ``chain3_count``, and ``chain3_abundance`` in
addition to MM, midpoint, and their absolute discrepancy.
[`tools/build_phase1d_structural_atlas.py`](/home/adnac/sorkin/tools/build_phase1d_structural_atlas.py)
keeps the same families, sizes, seeds, and no-embedding protocol
as Phase 1C.

The associated regression test
``tests/test_phase1d_structural_atlas.py`` enforces schema,
full grid coverage, family coverage, finite numeric values, and
``0 <= chain3_abundance <= 1``. The invariant-level tests in
``tests/test_causet_invariants.py`` check normalization on
chains/antichains and finiteness on Minkowski, KR, and corona
cells.

Interpretation is conservative: C3 abundance helps characterize
the ensemble structure and is particularly informative for
separating high-dimensional Minkowski sprinklings from dense
low-dimensional or layered controls, but it is not a standalone
manifoldness classifier. The right object remains the joint
trajectory of MM, midpoint, and structural chain abundance.
