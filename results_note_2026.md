# Results Note 2026

This note records what has been learned so far from reviving Bombelli's 1987 annealing program.

## Historical Baseline

The revived Python code preserves the structure of the thesis program:

- Pascal-style upper-triangular incidence input
- thesis-style `ran2` and `gasdev` random routines
- warmup phase
- annealing phase
- acceptance proportional to `4 * exp(-deltaE / T)`
- multiplicative cooling

The port therefore acts as a faithful historical baseline, not a new optimizer disguised as the old one.

## Schedule Sensitivity

On the reproducible benchmark `benchmarks/tesis_like_12.in`, the historical defaults give:

- initial temperature: `100.0`
- cooling factor: `0.9`
- mean final energy: `20.021376`
- zero rate: `0.00`

A refined empirical schedule search found:

- initial temperature: `180.0`
- cooling factor: `0.8`
- mean final energy: `0.166158`
- zero rate: `0.00`

This is about a `99.17%` reduction in mean final energy for that benchmark.
The main lesson is that the thesis algorithm is very schedule-sensitive.

## Statistical Interpretation

The revived project changes the question from a single embedding attempt to an ensemble measurement.

The useful observables are now:

- success probability
- median final energy
- mean final energy
- best seed
- timeout rate
- sensitivity to target dimension
- sensitivity to schedule

A single failed run is not evidence of non-embeddability.
Repeated failure across seeds, schedules, and dimensions is evidence of computational or structural difficulty.

## First Phase Diagram

The first frontier scan used:

- `n = 6, 12, 16, 24`
- `dim = 1, 2, 3, 4`
- `8` annealing seeds per cell
- `120s` timeout per run
- timeout energy marker: `9999`

Results:

| n | dim | runs | success | median final | note |
| ---: | ---: | ---: | ---: | ---: | --- |
| 6 | 4 | 8 | 0.50 | 0.000010 | best cell in the deeper frontier scan |
| 12 | 2 | 8 | 0.125 | 0.002364 | best median for n=12 |
| 16 | 2 | 8 | 0.125 | 0.392161 | frontier of useful search |
| 24 | 1-4 | 8 each | 0.00 | 9999.000000 | all runs timed out |

The `9999` values are not physical energies. They indicate timed-out runs.

The important result is that `n=24` is beyond the current brute-force annealing budget with these settings.

## Fast Frontier Scan

An adaptive fast scan was added to avoid spending the full budget on obviously poor cells.

It used:

- reduced annealing budget
- partial output after each cell
- early cell stopping
- timeout and status fields

Results:

| n | dim | runs | success | median final | interpretation |
| ---: | ---: | ---: | ---: | ---: | --- |
| 6 | 3 | 8 | 0.25 | 0.027485 | best fast cell by success |
| 6 | 4 | 8 | 0.125 | 0.074269 | still finds exact successes |
| 12 | 1 | 8 | 0.00 | 46.222369 | poor under fast budget |
| 12 | 2 | 8 | 0.00 | 62.739481 | poor under fast budget |
| 16 | 1 | 8 | 0.00 | 84.846101 | completes but poor |
| 16 | 2-4 | 3 each | 0.00 | 100+ | stopped early |
| 24 | 1-4 | 3 each | 0.00 | 170+ | stopped early |

The fast scan is not meant to produce final physics numbers. It is a triage tool.
It identifies where deeper runs are worth spending time.

## Main Discoveries So Far

1. The historical algorithm works, but it is fragile.

2. Schedule choice can change outcomes by orders of magnitude.

3. Embeddability should be treated statistically, not as a one-run yes/no question.

4. The first computational frontier appears between `n=16` and `n=24` for the current implementation and settings.

5. Higher target dimension does not monotonically improve the result. In several scans, higher `dim` increased cost and often worsened median energy.

6. Timeouts are themselves informative: they mark the boundary where the original method stops being a practical instrument.

## Working Interpretation

The revived thesis program now supports a stronger scientific question:

Which regions of causal-set space are embeddable with high probability, in which target dimension, and under which annealing schedules?

That is the experiment that was out of reach in 1987.

## v2 Diagnostic Foundation (in place)

A new layer of tooling has been added so that, going forward, the
question "did the optimizer fail because the causet is not
manifoldlike, or because annealing got stuck?" can be answered
honestly. See the *v2: Diagnostic Foundation* section of
[`research_agenda_2026.md`](/home/adnac/sorkin/research_agenda_2026.md)
for full discussion. The headline additions are:

- [`causet_invariants.py`](/home/adnac/sorkin/causet_invariants.py):
  order-theoretic diagnostics including a working Myrheim-Meyer
  dimension estimator. Verified against canonical sprinklings at
  ``d in {2, 3, 4}`` to within Monte Carlo tolerance.
- [`validation_suite.py`](/home/adnac/sorkin/validation_suite.py):
  a canonical uniform sprinkler for the Minkowski diamond (filling
  in for the legacy sprinkler, which is not uniform for
  ``d_spacetime >= 3``), plus the sprinkle-then-recover protocol
  with a Lorentz-invariant interval-matrix residual as the
  embedding quality metric.
- [`benchmarks/foundation/`](/home/adnac/sorkin/benchmarks/foundation):
  a 45-cell frozen grid of canonical sprinklings with precomputed
  invariants and an integrity test.
- Invariance tests on label permutation and composite Poincare
  transformations.

One subtlety surfaced during this work and is worth recording: the
Myrheim-Meyer ordering-fraction formula is commonly quoted with a
factor of ``4 Gamma(3d/2)`` in the denominator, which is the
*ordered-pair* convention ``R / [N (N - 1)]``. The project counts
unordered pairs everywhere (``R / C(N, 2)``), so the matching
formula is ``Gamma(d + 1) Gamma(d / 2) / (2 Gamma(3 d / 2))``.
Both conventions yield the same recovered dimension when each is
inverted against the matching empirical fraction.

## Phase 1: pure order-theoretic diagnostic (no embedding)

A minimal Phase 1 has been added that asks one question:

> Do invariants of the partial order alone separate manifoldlike
> sprinklings from non-manifoldlike controls?

This phase does **not** invoke any embedding algorithm, any
energy function, or any annealing. It compares two independent
order-theoretic dimension estimators on Minkowski sprinklings
and on Kleitman-Rothschild three-level posets at matched
``(n, seed)`` cells. See the section "Phase 1: order-theoretic
negative controls" in
[`research_agenda_2026.md`](/home/adnac/sorkin/research_agenda_2026.md)
for the components.

The reproducible artifact lives at
[`benchmarks/foundation/phase1_atlas.{csv,md}`](/home/adnac/sorkin/benchmarks/foundation/phase1_atlas.csv)
and is regenerated with ``make regen-phase1``. First-pass
empirical observation (ensembles of 5 seeds at ``n in {16, 32,
64}``):

- For Minkowski sprinklings the Myrheim-Meyer estimator
  converges toward the true ``d_spacetime`` as ``n`` grows
  (``MM(d=2) -> 2.01``, ``MM(d=3) -> 3.22``, ``MM(d=4) -> 3.92``
  at ``n = 64``).
- For Kleitman-Rothschild the Myrheim-Meyer estimator is
  essentially flat in ``n`` (``MM ~ 2.36`` at every ``n``). The
  ``n``-dependence of MM is therefore itself a discriminator:
  manifoldlike causets *converge*, KR controls do not.
- The midpoint scaling estimator has substantial finite-size
  bias at these ``n``. It undershoots in higher dimensions
  (``midpoint(d=4, n=64) ~ 2.49``) and overshoots on KR
  (``midpoint(KR, n=64) ~ 3.08``).
- At ``n = 64``, ``d_spacetime = 2`` the two estimators *together*
  separate Minkowski (MM ~ 2.01, midpoint ~ 1.91) from KR
  (MM ~ 2.36, midpoint ~ 3.08) cleanly. For higher target
  dimensions, the separation is muddier and a richer invariant
  set (or larger ``n``) will be needed.

The conservative reading: *current* order-theoretic invariants
distinguish KR from ``1+1`` Minkowski sprinklings at ``n = 64``,
but they do not yet provide robust per-causet separation across
all ``d_spacetime``. Phase 1 makes this honest gap visible,
which is itself the point: we now know what the foundation can
and cannot do before any annealing is invoked.

## Phase 1B: finite-size scaling of the order-theoretic atlas

Phase 1 found that at fixed ``n = 64`` the two estimators
jointly separate Minkowski ``d = 2`` from Kleitman-Rothschild,
but the separation degrades for higher target dimensions.
Phase 1B asks the next question: does the separation improve
systematically as ``n`` grows?

The reproducible artifact is
[`benchmarks/foundation/phase1b_scaling_atlas.{csv,md}`](/home/adnac/sorkin/benchmarks/foundation/phase1b_scaling_atlas.md),
ensemble-aggregated over the same five seeds as Phase 1 at
``n in {32, 64, 128, 256}`` for the same families. Regenerable
via ``make regen-phase1b``. Headline numbers (ensemble means):

| family | n | mean MM | mean midpoint | mean \|disc\| |
| --- | ---: | ---: | ---: | ---: |
| Minkowski d=2 | 32 | 1.93 | 1.65 | 0.28 |
| Minkowski d=2 | 256 | 2.02 | 2.06 | 0.07 |
| Minkowski d=3 | 32 | 3.32 | 1.74 | 1.57 |
| Minkowski d=3 | 256 | 3.06 | 3.00 | 0.28 |
| Minkowski d=4 | 32 | 3.71 | 1.55 | 2.15 |
| Minkowski d=4 | 256 | 4.07 | 3.80 | 0.42 |
| Kleitman-Rothschild | 32 | 2.34 | 2.26 | 0.17 |
| Kleitman-Rothschild | 256 | 2.37 | 4.71 | 2.34 |

Three observations are robust across the grid:

1. **Manifoldlike cells converge.** For every Minkowski
   ``d_spacetime in {2, 3, 4}``, ``mean_mm`` lands within ~0.1
   of the truth at ``n = 256`` and ``std_mm`` shrinks by roughly
   a factor of two between ``n = 32`` and ``n = 256``. The
   midpoint estimator, which had substantial finite-size bias
   in Phase 1, also converges toward the truth at the same
   ``n``.

2. **Kleitman-Rothschild does not converge to any common
   dimension.** ``mean_mm`` is essentially flat (~2.36) in
   ``n``; the midpoint estimator grows approximately linearly
   in ``log n`` (2.26, 3.08, 3.85, 4.71 at ``n = 32, 64, 128,
   256``). The two estimators disagree by a value that
   *increases* with ``n`` (``|disc| = 0.17, 0.72, 1.49, 2.34``).

3. **The scaling trajectory itself is the diagnostic.** The
   single-cell separation that Phase 1 found insufficient at
   higher ``d`` becomes clean once the trajectory across ``n``
   is read: Minkowski cells exhibit shrinking inter-estimator
   discrepancy, KR cells exhibit growing discrepancy. The sign
   of ``d|disc|/dn`` is opposite for the two families across
   the entire grid.

The conservative reading: order-theoretic invariants now
provide a working *ensemble-level* diagnostic for manifoldlike
vs. non-manifoldlike structure across ``d_spacetime in {2, 3,
4}`` and ``n in [32, 256]``, provided one reads the
finite-size trajectory rather than any single ``n``. One honest
limit remains:

- **Ensemble-level, not per-causet.** At any fixed
  ``(d_spacetime, n)``, the seed-to-seed spread in ``mm_dim``
  and ``midpoint_dim`` is wide enough that one cannot, from a
  single causet's invariants alone, certify manifoldness.
  Ensemble-level separation does not yet imply per-causet
  classification.

## Phase 1C: second non-manifoldlike control

Phase 1C addresses the main limitation of Phase 1B by adding a
second cheap negative control: suspended corona/crown posets.
The reproducible artifact is
[`benchmarks/foundation/phase1c_scaling_atlas.{csv,md}`](/home/adnac/sorkin/benchmarks/foundation/phase1c_scaling_atlas.md),
regenerable via ``make regen-phase1c``. It keeps the same
``n in {32, 64, 128, 256}`` grid and five seeds, and compares
Minkowski sprinklings against both Kleitman-Rothschild and
corona controls.

Headline control rows (ensemble means):

| family | n | mean MM | mean midpoint | mean \|disc\| |
| --- | ---: | ---: | ---: | ---: |
| Kleitman-Rothschild | 32 | 2.34 | 2.26 | 0.17 |
| Kleitman-Rothschild | 256 | 2.37 | 4.71 | 2.34 |
| corona poset | 32 | 1.88 | 4.00 | 2.12 |
| corona poset | 256 | 1.98 | 7.00 | 5.02 |

The result supports the Phase 1B reading. The two controls do
not share the same Myrheim-Meyer plateau (KR sits near 2.36,
the corona construction approaches 2.0), but both exhibit the
same qualitative non-manifoldlike trajectory: ``mean_mm`` is
comparatively flat within the family, while ``mean_midpoint``
grows with the finite-size scale instead of converging toward a
common dimension. The growing inter-estimator discrepancy is
therefore not specific to the KR three-level random generator.

## Phase 1D: third order-theoretic invariant

Phase 1D adds a third embedding-free observable: 3-chain
abundance. The raw count ``C3`` is the number of triples
``i prec j prec k``; the normalized observable is
``C3 / binom(n, 3)``. This is not inverted to a dimension
estimator, because no analytic calibration is being assumed.
The reproducible artifact is
[`benchmarks/foundation/phase1d_structural_atlas.{csv,md}`](/home/adnac/sorkin/benchmarks/foundation/phase1d_structural_atlas.md),
regenerable via ``make regen-phase1d``.

Headline rows at ``n = 256`` (ensemble means):

| family | d | mean MM | mean midpoint | mean \|disc\| | mean C3 abundance |
| --- | :---: | ---: | ---: | ---: | ---: |
| Minkowski | 2 | 2.02 | 2.06 | 0.07 | 0.1611 |
| Minkowski | 3 | 3.06 | 3.00 | 0.28 | 0.0206 |
| Minkowski | 4 | 4.07 | 3.80 | 0.42 | 0.0023 |
| Kleitman-Rothschild | - | 2.37 | 4.71 | 2.34 | 0.0474 |
| corona poset | - | 1.98 | 7.00 | 5.02 | 0.0117 |

The observable improves the ensemble-level structural picture,
especially for the higher-dimensional sprinklings: Minkowski
``d = 3`` and ``d = 4`` have much lower 3-chain abundance than
Minkowski ``d = 2`` and sit on different trajectories from KR.
It also shows why no single raw invariant should be treated as a
manifoldness test: corona controls have a C3 abundance comparable
to the ``d = 3`` scale at large ``n`` while their midpoint
discrepancy grows sharply. Thus Phase 1D strengthens the
multi-observable ensemble diagnostic, but it still does not give
robust per-causet classification at ``n <= 256``.

## Next Plan

The next step is now a two-stage search:

1. Use `phase_diagram.py --fast-frontier` to identify promising cells.
2. Use `phase_refine.py` to rerun only those cells with deeper annealing and more seeds.

For now, the manual refinement queue is:

- `n=6, dim=3`
- `n=6, dim=4`
- `n=12, dim=2`
- possibly `n=16, dim=2`, but only with careful schedule tuning

Cells like `n=24` should not be attacked with uniform brute force until the algorithm is improved.

The exact command for this second phase is:

```bash
python3 phase_refine.py --cells 6:3,6:4,12:2,16:2 \
  --seed-start 1959 --seed-count 16 \
  --gpu-first --backend auto \
  --partial cell \
  --runs-csv /tmp/phase_refine_runs.csv \
  --summary-csv /tmp/phase_refine_summary.csv \
  --report-md /tmp/phase_refine.md \
  --heatmap-svg /tmp/phase_refine.svg
```

This keeps the scientific protocol honest: the fast scan decides where to look, but the deeper rerun is what should be used for interpretation.
