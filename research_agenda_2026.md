# Research Agenda 2026

This note frames the revived Bombelli 1987 annealing program as a research instrument rather than only a historical port.

SORKIN-2 update: the v1.0.0 revival is complete as a historical reconstruction.
The active continuation is algorithmic recoverability in the Bombelli annealer:
the gap between known zero-energy causal realization and accessibility by the
historical energy, move set, schedule, and acceptance rule. This is not a
causal-set physics classification program, and ML is deferred until a clean
multi-family known-truth dataset exists.

Current empirical results are summarized in [`results_note_2026.md`](/home/adnac/sorkin/results_note_2026.md).

## Central Question

Which known-truth causal sets are recoverable with zero energy or verified exact relation recovery under the current historical annealing pipeline, and what controls the transition between accessible and inaccessible cases?

The old program asks whether one causet can be fitted by a particular constructive procedure. The modern diagnostic should ask how known-truth families behave under that procedure.

## Working Hypotheses

- Optimizer recoverability is not controlled only by `n`; density, dimension, and order structure matter.
- The annealing landscape has basins of attraction that can be mapped statistically.
- Schedule sensitivity may become diagnostically informative if it correlates reproducibly with causal-set structure, but it is not by itself a physical observable.
- Failure to reach zero energy separates algorithmic accessibility from existence only when the case has known truth; it does not establish non-embeddability.

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
- lowest-final-energy seed per region under the current optimizer metric
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

The scientific target is to identify which order-theoretic features predict optimizer recoverability under the current embedding pipeline.

4. Compare energy and move-set diagnostics.

Keep the thesis energy as the historical baseline, then test variants:

- original Bombelli energy
- regularized radius penalties
- penalties that weight links differently from transitive relations
- interval-aware penalties
- multiobjective versions that separate false positives from false negatives

The question is not only whether the old energy works. It is which historical energy/move-set choices control optimizer success or failure on known-truth reference cases.

5. Separate algorithmic inaccessibility from existence.

For hard cases, rerun with stronger methods:

- more seeds
- broader schedules
- local refinements near the lowest-energy configuration found under the current metric
- higher target dimension
- alternative energy functions

If all methods fail, treat the case as a candidate hard instance for the current pipeline, not as evidence of non-embeddability.
If only the original schedule fails, the lesson is algorithmic.

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

- Energy near zero means a low-energy configuration was found under the current objective, not that the causet is uniquely geometric.
- High final energy after one run is not evidence of non-embeddability.
- Repeated high energy across many schedules and seeds is evidence of algorithmic difficulty.
- If higher dimension fixes the problem, the original target dimension was too restrictive.
- If schedule tuning fixes the problem, the original historical run was computationally limited.

## First Milestone

The first serious milestone is a phase diagram:

- x-axis: number of elements `n`
- y-axis: target dimension `dim`
- cell value: success probability

A second layer should show median final energy.

This is the experiment that Bombelli and Sorkin could not practically do in 1987: not a single optimized causet, but a statistical map of optimizer recoverability under the current objective.

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
- [`phase_diagram.py`](/home/adnac/sorkin/phase_diagram.py): n-vs-dim optimizer-recoverability maps
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
*algorithmic* failure from *candidate structural obstruction*. Without
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
dimension is an *independent* diagnostic: if it indicates a
well-defined ``d`` near the target embedding dimension, the failure
is plausibly algorithmic. If it diverges or sits far from the target,
that is a structural warning, not a verdict from the annealer.

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

## Phase 2: embedding bridge

Phase 2 begins the connection back to the historical
Bombelli-Sorkin embedding code without turning the project into
an optimizer search. The artifact is
[`benchmarks/foundation/phase2_embedding_bridge.{csv,md}`](/home/adnac/sorkin/benchmarks/foundation/phase2_embedding_bridge.md),
regenerable via ``make regen-phase2``.

The first bridge uses a deliberately small fixed probe:
``n = 64``, case seed ``1959``, optimizer seed ``1987``, one
case each for Minkowski ``d = 2, 3, 4``, Kleitman-Rothschild,
and suspended corona. The schedule is short
(``warmup_limit = 10``, ``anneal_limit = 10``, ``max_data = 4``)
so the bridge is cheap enough to keep under regression. It is
therefore a diagnostic probe, not a production annealing run.

The CSV records the Phase 1D pre-embedding observables
(``mm_dim``, ``midpoint_dim``, discrepancy, raw chain counts,
and ``chain3_abundance``) beside embedding outputs available
from the existing code: ``initial_energy``, ``warmup_energy``,
``final_energy``, ``truth_energy`` and ``interval_rmse`` for
Minkowski sprinklings with known coordinates, and ``NA`` for
truth-dependent columns on controls.

The immediate diagnostic use is triage. If a Minkowski case has
good structural diagnostics but a large energy gap or interval
RMSE, the failure is evidence about the annealing schedule or
energy landscape, not evidence that the causet is
non-manifoldlike. If a control has bad structural diagnostics
but a deceptively modest final energy, the energy functional
needs auditing before it is used for any interpretation.

## Phase 2B: annealer schedule probe

Phase 2B is a minimal diagnostic of where the Phase 2 failure
mode sits on the budget/quality curve. The Phase 2 short
schedule leaves Minkowski cases with large ``energy_gap`` and
large ``interval_rmse`` despite a known ``truth_energy = 0``.
The natural follow-up question is whether widening the
historical annealer's budget alone closes that gap.

The artifact is
[`benchmarks/foundation/phase2b_annealer_schedule_probe.{csv,md}`](/home/adnac/sorkin/benchmarks/foundation/phase2b_annealer_schedule_probe.md),
regenerable via ``make regen-phase2b``. The probe varies only
``warmup_limit``, ``anneal_limit`` and ``max_data`` over a
small Minkowski-only grid in ``d in {2, 3, 4}``, ``n in
{32, 64}``, three case seeds drawn from the Phase 1B atlas
(the first three of ``SEEDS``), and one optimizer seed
(``1987``). The temperature schedule (``initial_temp = 100``,
``cooling_factor = 0.9``) is held fixed so the comparison is a
pure budget probe rather than a temperature search. Three
schedule labels — ``short`` (matching Phase 2 exactly),
``medium`` and ``long`` — span a reconfigure-budget range of
50 -> 140 -> 330 (roughly 7x).

By construction Phase 2B does **not** include
Kleitman-Rothschild or suspended corona controls. They have
no ground-truth coordinates, so ``truth_energy``,
``energy_gap`` and ``interval_rmse`` are undefined for them.
This probe is therefore not a manifoldness classifier and
should never be cited as one. The conservative interpretation
rule is fixed: a budget-induced drop in ``energy_gap`` is
evidence of schedule/optimizer failure, and a flat gap across
budgets is evidence about the energy or implementation, not
about the manifoldlikeness of the underlying sprinkling.

Phase 2B does not introduce any new optimizer (no basin
hopping, parallel tempering, or ML-assisted initialization).
Those remain explicitly outside the v2 foundation scope.

## Phase 2C: oracle embedding audit

Phase 2B narrowed the failure candidates to energy definition,
parametrization, or move set. Phase 2C eliminates the energy
and convention candidates by evaluating the Bombelli energy
*directly at the ground-truth coordinates*, with no annealing.
The artifact is
[`benchmarks/foundation/phase2c_oracle_embedding_audit.{csv,md}`](/home/adnac/sorkin/benchmarks/foundation/phase2c_oracle_embedding_audit.md),
regenerable via ``make regen-phase2c``.

Three oracle checks are run for each Minkowski case in the
Phase 2/2B grid (d ∈ {2,3,4}, n ∈ {32,64}, three seeds):

1. **oracle_pass_energy** — Bombelli energy at ground-truth
   coordinates is numerically zero (threshold |E| ≤ 1e-9).
2. **oracle_pass_causal_matrix** — causal matrix reconstructed
   from the stored coordinates matches the stored matrix
   bit-for-bit (zero discordant pairs).
3. **oracle_pass_interval_rmse** — Lorentz-invariant RMSE of
   the truth embedding against itself is numerically zero.

Verdict: **ORACLE PASSES** on all 18/18 cases.

Conservative interpretation:

- The Bombelli energy formula returns exactly 0.0 at the
  ground-truth coordinates in every case. The energy objective
  correctly identifies the truth as the global minimum.
- The causal matrices reconstructed from stored coordinates
  match the originals with zero discordant pairs in every case.
  No floating-point drift or sign-convention mismatch.
- The Lorentz-invariant interval residual is identically zero
  when comparing the truth to itself.
- All three checks pass; there is **no convention or formula
  inconsistency** to fix.
- The failure in Phase 2/2B is therefore localized to the
  **optimizer**: move set, initialization, or annealing
  landscape. More budget (Phase 2B) did not help; the optimizer
  is not finding the geometry the energy already recognises.
- The next diagnostic step is a move-set or initialization
  audit — not an energy redesign, and not more budget.

Phase 2C does not introduce a new optimizer or touch the energy
definition.

## Phase 2D: initialization / basin audit

Phase 2C localized the failure to the optimizer: move set,
initialization, or annealing landscape. Phase 2D discriminates
between those candidates by injecting four controlled
initializations into ConesSimulator and measuring what happens
during warmup and annealing. The artifact is
[`benchmarks/foundation/phase2d_initialization_basin_audit.{csv,md}`](/home/adnac/sorkin/benchmarks/foundation/phase2d_initialization_basin_audit.md),
regenerable via ``make regen-phase2d``.

Four initialization strategies are run on the same Minkowski
grid as Phase 2B (d ∈ {2,3,4}, n ∈ {32,64}, seeds
1959/1962/1987) with the Phase 2 short schedule
(warmup_limit=10, anneal_limit=10, max_data=4):

1. **truth** — exact ground-truth coordinates (ε = 0).
2. **truth_plus_small_noise** — truth + Gaussian noise at
   ε = 1e-3 per coordinate.
3. **truth_plus_medium_noise** — truth + Gaussian noise at
   ε = 5e-2 per coordinate.
4. **random_init** — the historical ConesSimulator default
   (rnew[i] = i+2, xnew[i] = 0, then energy()+update()).

Custom initialization bypasses ``startup()`` by setting
``sim.rnew``/``sim.xnew`` directly and calling
``sim.energy()`` + ``sim.update()``. No changes are made to
``cones.py``.

Metrics recorded at two checkpoints (last-accepted positions
``rold``/``xold``):

- ``initial_energy`` / ``final_energy``: Bombelli energy
  before and after warmup+anneal.
- ``initial_interval_rmse`` / ``final_interval_rmse``:
  Lorentz-invariant RMSE relative to scaled truth.
- ``initial_distance_to_truth_rms`` /
  ``final_distance_to_truth_rms``: RMS coordinate distance
  to ground-truth positions.
- ``preserved_near_truth``: ``delta_energy ≤ 1e-6``.

Verdict: **NARROW_BASIN**.

Truth initialization (energy = 0) is preserved exactly in
18/18 cells. The warmup exits immediately because
``energies[0] ≤ 0`` and no moves are made. Any
positive-energy perturbation — even ε = 1e-3, which gives
initial energy ≈ 0.005 — activates the warmup, which makes
unconditional accepts for 10 steps and systematically
scrambles the configuration. The mean final energy for
small-noise starts is 18.9 vs. an initial mean of 0.005.
Medium-noise and random starts both worsen substantially.

Conservative interpretation:

- The warmup loop ``while count < warmup_limit and
  energies[0] > 0`` exits immediately at energy = 0.
  For any positive energy it makes 10 unconditional moves
  with no Metropolis criterion, deliberately exploring
  high-energy regions before annealing.
- The effective basin of attraction is extremely narrow in this grid: only
  the exact zero-energy configuration is preserved.
- This is a **warmup-dynamics failure**, not a move-set
  failure in isolation, and not an energy failure (Phase 2C).
- Recommended next step: skip warmup or replace it with a
  conditioned equilibration when starting near a known
  low-energy configuration. No new optimizer is assumed.

Phase 2D does not introduce a new optimizer or modify the
energy definition or the move set.

## Phase 2E: warmup-skip probe

Phase 2D attributed the near-truth destruction to the unconditional
warmup loop. Phase 2E tests this directly with a paired comparison:
the same grid and the same four initialization strategies are run
with ``with_warmup`` (Phase 2D baseline) and ``skip_warmup``
(anneal-only) side by side. Each physical case is identified by a
``paired_key`` so the two modes are compared row-for-row. The
artifact is
[`benchmarks/foundation/phase2e_warmup_skip_probe.{csv,md}`](/home/adnac/sorkin/benchmarks/foundation/phase2e_warmup_skip_probe.md),
regenerable via ``make regen-phase2e``.

The diagnostic question: is the warmup the primary cause of
near-truth destruction, or does the annealing phase fail
independently?

Verdict: **WARMUP_IS_PRIMARY_FAILURE**.

Skipping the warmup improves small-noise preservation (17/18
vs 16/18 with warmup) and reduces mean final energy for small-
noise starts from 18.9 to 12.1. Medium-noise starts are also
lower without warmup (mean 286 vs 396), although both remain
non-preserved (0/18). Random-init unexpectedly also benefits
from skipping warmup (11/18 vs 8/18 preserved).

Paired deltas (skip − with, mean over all cells):

| init | mean Δ final E | skip pres | with pres |
| --- | ---: | ---: | ---: |
| truth | 0.0000 | 18/18 | 18/18 |
| truth_plus_small_noise | −6.8 | 17/18 | 16/18 |
| truth_plus_medium_noise | −109.8 | 0/18 | 0/18 |
| random_init | −97.3 | 11/18 | 8/18 |

Conservative interpretation:

- The warmup phase is confirmed as the primary cause of near-truth
  destruction. Removing it strictly improves or preserves all
  label categories.
- The one small-noise row that fails even without warmup indicates
  the annealing phase has some residual instability near the
  truth minimum. Skipping warmup is necessary but not universally
  sufficient for small perturbations.
- Medium-noise starts are not recovered by either mode. The basin
  that the anneal-only phase can hold is narrower than ε = 5e-2.
- Random-init improvement is a secondary finding: the 10
  unconditional warmup steps are not helping random starts
  either; they raise the energy from the default linear ladder
  rather than exploring productively.
- Next diagnostic step: a conditioned equilibration or
  energy-gated warmup that makes unconditional accepts only when
  the proposed move does not increase energy by more than a
  threshold.

Phase 2E does not introduce a new optimizer or modify the energy,
move set, or cooling schedule.

## Phase 2F: guarded-warmup probe

Phase 2E confirmed that the warmup is the primary failure. Phase 2F
tests whether a non-destructive alternative — the guarded warmup —
can preserve the exploratory intent of warmup for random starts
while not destroying near-truth configurations. The artifact is
[`benchmarks/foundation/phase2f_guarded_warmup_probe.{csv,md}`](/home/adnac/sorkin/benchmarks/foundation/phase2f_guarded_warmup_probe.md),
regenerable via ``make regen-phase2f``.

Three warmup modes are compared in a three-way paired comparison on
the same Phase 2D/2E grid (d∈{2,3,4}, n∈{32,64}, seeds
1959/1962/1987, 4 init labels) = 216 rows:

1. **legacy_warmup** — ``sim.warmup(buf)`` + ``sim.anneal(buf)``
   (Phase 2E baseline).
2. **skip_warmup** — ``sim.anneal(buf)`` only (Phase 2E baseline).
3. **guarded_warmup** — energy-gated warmup + ``sim.anneal(buf)``.
   External wrapper; no changes to ``cones.py``.

Guarded-warmup implementation:

- ``GUARD_THRESHOLD = 0.0`` (strictly non-worsening).
- Accept a proposed move iff ``sim.deltae ≤ 0`` (pre-normalization).
- On rejection: clear ``change[i]`` flags, restore ``sim.rave``
  to ``sim.r``. ``rnew``/``xnew`` are not explicitly reset because
  ``reconfigure()`` overwrites them from ``rold``/``xold`` at the
  start of each new call.
- Records ``warmup_attempted_moves``, ``warmup_accepted_moves``,
  ``warmup_rejected_moves``, ``warmup_energy_before``,
  ``warmup_energy_after``, ``warmup_delta_energy`` per row.

Verdict: **GUARDED_WARMUP_FIXES_PRIMARY_FAILURE**.

Per-label aggregate (mean over 18 cells per label):

| init | legacy final E | skip final E | guarded final E | guarded pres |
| --- | ---: | ---: | ---: | ---: |
| truth | 0.0000 | 0.0000 | 0.0000 | 18/18 |
| truth_plus_small_noise | 18.92 | 12.12 | **0.0013** | **18/18** |
| truth_plus_medium_noise | 395.80 | 286.03 | 255.32 | 0/18 |
| random_init | 405.16 | 307.89 | **271.42** | **12/18** |

Conservative interpretation:

- **guarded_warmup dominates on every label.** Small-noise: 18/18
  preserved, mean final E = 0.0013 (vs 12.12 skip, 18.92 legacy).
  The energy-gated warmup removes the observed small-noise failure
  on this Phase2F grid.
- **random_init also improves** (12/18 vs 11/18 skip, 8/18
  legacy). The greedy warmup moves do help random starts explore
  without the unconditional-accept damage.
- **Medium-noise not recovered** by any mode (0/18 all). The
  anneal-only basin does not extend to ε = 5e-2 regardless of
  warmup strategy. This is a move-set or cooling-schedule limit,
  not a warmup limit.
- **Normalization note.** The guard is applied to pre-normalization
  ``sim.deltae``; post-normalization ``warmup_energy_after`` may
  differ from ``warmup_energy_before`` by a normalization factor.
  This is documented and does not invalidate the guard.
- This is a diagnostic finding, not a new optimizer. The guarded
  warmup is an external wrapper around the same
  ``ConesSimulator`` internals. No new energy, no new move set,
  no changes to ``cones.py``.

Phase 2F does not introduce a new optimizer or modify the energy,
move set, or cooling schedule.
