# Phase 4D — extension: control by n for the order-theoretic correlate audit

## Scope statement

This is a **descriptive audit extension** to an existing script. It is **not**
a new phase, **not** a new script, **not** a change to the algorithm or to
`cones.py`, and **not** a new mechanical verdict.

Its single purpose is to answer one question:

> Does the existing per-seed Spearman correlation between causal-set
> invariants and optimizer-robustness metrics in Phase 4D survive
> controlling for `n`, or is it driven by finite-size scaling?

The answer is reported as numbers in the existing CSV and as a new section in
the existing Markdown audit. **No new DETECTED / NO_DETECTED label is
introduced.** The original Phase 4D verdict logic is left untouched.

---

## Context (what already exists)

Phase 4D is an existing script that:

1. Reads outputs of Phase 4C:
   - `phase4c_optimizer_seed_probe_per_run.csv`
   - `phase4c_optimizer_seed_probe_per_cell_epsilon.csv`
2. Recomputes causal-set invariants per `(n, target_dim, causet_seed)`.
3. Cross-correlates those invariants with optimizer-robustness metrics:
   - `IQR_loss_K`
   - `mean_loss_K`
   - `floor_saturated_fraction_K`
   - `label_distinct_shapes_K` (label stability)
4. Reports Spearman/Pearson at two levels:
   - per-seed, N=90
   - per-cell, N=9
5. Produces:
   - `benchmarks/foundation/phase4d_robustness_per_seed.csv`
   - `benchmarks/foundation/phase4d_robustness_per_cell.csv`
   - `benchmarks/foundation/phase4d_robustness_audit.md`

Current mechanical verdict on the full run:
`ORDER_THEORETIC_CORRELATE_DETECTED`.

The verdict is driven mostly by `relation_count` and `chain2_count` vs
`floor_saturated_fraction_K`, with Spearman `|ρ| ≈ 0.94`.

These are extensive counts. The grid uses only `n ∈ {32, 48, 64}`, so the
correlation may be finite-size scaling, not an order-theoretic correlate.
For contrast, `label_distinct_shapes_K`, which is the metric that actually
measures morphological instability under `optimizer_seed`, reaches only
`max |ρ| ≈ 0.55`.

---

## Hard constraints

1. **Do not create Phase 4E.** Extend the existing Phase 4D script in place.
2. **Do not create new scripts, modules, or top-level files** other than
   tests for the new helpers (tests go wherever existing Phase 4D tests live).
3. **Do not modify `cones.py`.**
4. **Do not hand-edit any CSV.** All CSV changes happen by re-running
   Phase 4D after the code change.
5. **Do not add new invariants or new dimensionless ratios.** Use only the
   invariants Phase 4D already computes.
6. **Do not introduce a new automatic verdict** based on partial or
   stratified correlations. The Phase 4D verdict logic stays exactly as is.
7. **n-control is computed only at per-seed level (N=90).** Per-cell (N=9)
   is too small for partial or stratified statistics; it stays as
   descriptive sanity, unchanged.

---

## Step 0 — Discover the existing code (mandatory before any edit)

Before writing any code:

- Locate the Phase 4D script in the repository (search by output filename
  `phase4d_robustness_audit.md`, by the input filename
  `phase4c_optimizer_seed_probe_per_run.csv`, or by the phase name).
- Read it completely. Identify and report back:
  - The function that computes per-seed correlations.
  - The function that writes `phase4d_robustness_per_seed.csv`.
  - The function that writes `phase4d_robustness_audit.md`.
  - The exact column name in the per-seed dataframe that holds the value
    of `n` (it may be `n`, `N`, `n_elements`, `causet_n`, etc.).
  - The exact column names already used for the raw Spearman output, so
    new columns follow the same naming convention and dtype.
  - Where Phase 4D tests live, and the test naming convention.

Do not start editing until that map is on the table.

---

## Step 1 — Add partial Spearman controlling for n (per-seed only)

For each `(invariant, target)` pair already correlated at per-seed level,
add a partial Spearman correlation controlling for `n`.

Implementation:

- Method: rank-transform invariant, target, and `n` independently, then
  compute partial Pearson by residualizing rank(invariant) and rank(target)
  on rank(n) and correlating the residuals. This is the standard
  rank-partial-correlation method.
- Add to the per-seed output (next to the existing raw Spearman columns,
  matching their dtype and the existing naming convention discovered in
  Step 0):
  - one column for the partial Spearman coefficient,
  - one column for the partial Spearman p-value **if and only if** the
    existing raw Spearman output also carries a p-value column,
  - one column documenting the method (a short string such as
    `rank-partial-pearson_residualize_n`).

Do not introduce a threshold or "significant / not significant" decision
based on these columns.

---

## Step 2 — Add stratified Spearman by n (per-seed only)

For each `(invariant, target)` pair, compute regular Spearman within each
`n` stratum. With `n ∈ {32, 48, 64}` and N=90 total, each stratum has N=30.

Add to the per-seed output:

- One column per stratum holding the Spearman coefficient within that
  stratum. Naming: follow whatever convention is already used in
  Phase 4D, e.g. `spearman_n32`, `spearman_n48`, `spearman_n64`.
- One summary column holding the minimum of `|ρ|` across the three
  strata. Rationale: if the raw correlation is driven by `n` alone, then
  within any single `n` the correlation collapses, so the minimum
  absolute value across strata is a stricter check than the average.

Do not introduce a threshold or verdict on the summary column.

---

## Step 3 — Update the audit Markdown

In `phase4d_robustness_audit.md`, add a new section after the existing
top-level correlation tables. Title it:

> `n-control (interpretive layer)`

Inside it:

1. A short paragraph, included **literally**:

   > The n-control is interpretive auditing, not mechanical
   > reclassification. The Phase 4D verdict
   > (`ORDER_THEORETIC_CORRELATE_DETECTED` or otherwise) is computed from
   > raw Spearman as before. The partial-n and stratified-by-n values are
   > reported here so the reader can judge whether the raw correlation
   > reflects n-scaling rather than an order-theoretic correlate.

2. For each `(invariant, target)` pair appearing in the existing raw
   top-correlations table at per-seed level, show three values side by
   side:
   - raw Spearman,
   - partial Spearman controlling for n,
   - the stratified summary `min_abs` across the three n-strata
     (and, if it fits the existing table style, the three per-stratum
     values too).

3. A **size-like-invariants caveat block**. List explicitly the
   invariants that are extensive in `n` and therefore expected to
   correlate with anything that also grows with `n`:
   - `relation_count`
   - `chain2_count`
   - `chain3_count`
   - `chain4_count`
   - `height`

   Whenever any of these appears as the top raw correlate for a given
   target, the caveat block must be rendered in the prose immediately
   under the table where that target's top correlate is shown. Wording
   (literal):

   > Top raw correlate is an extensive size-like invariant. The raw
   > correlation is expected to track `n` independently of any
   > order-theoretic content. Refer to the partial-n and stratified-n
   > columns to judge whether the correlation survives after the
   > finite-size effect is removed.

The per-cell section of the audit is left unchanged. Add one sentence at
the top of that section: *"No n-control is computed at per-cell level
(N=9 is too small for partial or stratified statistics; per-cell is
reported as descriptive sanity only)."*

---

## Step 4 — Tests

Add unit tests for the new partial-Spearman and stratified-Spearman
helpers. Tests go wherever existing Phase 4D tests live, following the
existing naming convention.

Minimum coverage:

1. Pure-`n` synthetic: construct a synthetic dataset where `y` is a
   monotone function of `n` and `x` is also a monotone function of `n`,
   but `x` and `y` are conditionally independent given `n`. Assert:
   - raw Spearman `|ρ|` is large (close to 1),
   - partial Spearman controlling for `n` is small (close to 0),
   - the stratified `min_abs` is small (close to 0).

2. Genuine signal synthetic: construct a synthetic dataset where `y`
   depends on `x` directly and not on `n`. Assert:
   - raw Spearman `|ρ|` is large,
   - partial Spearman controlling for `n` is still large,
   - the stratified `min_abs` is still large.

3. Existing Phase 4D tests must keep passing unchanged.

---

## Step 5 — Re-run Phase 4D end-to-end

Re-run Phase 4D with the existing inputs from Phase 4C.

Verify:

- `phase4d_robustness_per_seed.csv` has the new columns from Steps 1 and
  2 populated. All existing columns and rows are unchanged (same number
  of rows, same values in the original columns).
- `phase4d_robustness_per_cell.csv` is byte-equivalent to the previous
  run (per-cell did not get n-control).
- `phase4d_robustness_audit.md` has the new section, with all three
  numbers side by side per `(invariant, target)`, the literal
  interpretive disclaimer, and the size-like caveat block where
  applicable.
- The Phase 4D verdict variable, whatever it is internally, has the same
  value as before. The n-control did not alter it.
- `pytest` exit code is 0. Report the test summary line.

---

## Acceptance criteria (single checklist)

- [ ] Step 0 map reported (script path, function names, n column name,
      existing column naming convention, tests location).
- [ ] Partial Spearman columns added in per-seed CSV.
- [ ] Stratified Spearman columns + `min_abs` summary added in
      per-seed CSV.
- [ ] Per-cell CSV unchanged.
- [ ] Audit Markdown has the new section with literal disclaimer and
      caveat block.
- [ ] New tests pass; existing tests still pass.
- [ ] Phase 4D verdict logic untouched (no new DETECTED / NO_DETECTED
      branch was added).

---

## Out of scope (explicit non-goals)

- Modifying the optimizer or `cones.py`.
- Adding new invariants or new dimensionless ratios.
- Introducing coarse-graining.
- Comparing or fusing the per-seed and per-cell verdicts.
- Re-running Phase 4C or any earlier phase.

---

## What the user will decide after this extension runs

If the partial and stratified Spearman values collapse on the size-like
invariants (`relation_count`, `chain2_count`, etc.), and
`label_distinct_shapes_K` correlations stay weak, the conservative
reading is that Phase 4D's raw `ORDER_THEORETIC_CORRELATE_DETECTED` was
finite-size scaling. The order-theoretic-correlate branch closes
without modifying the algorithm.

If, instead, partial and stratified Spearman values stay strong on at
least one non-extensive invariant, or on `label_distinct_shapes_K`,
that is a candidate signal worth a separate planning step (which is
**not** part of this extension).
