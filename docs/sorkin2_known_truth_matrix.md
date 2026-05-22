# SORKIN-2 Known-Truth Matrix

## Purpose

This document defines the design of the minimal known-truth case matrix for
SORKIN-2. It is a design document, not an experimental result. No experiments
are run here.

SORKIN-2 studies algorithmic recoverability in the historical Bombelli annealer.
It does not study physical embeddability in the general sense.

### Central Distinction

Two questions must be kept separate throughout:

**Existence question (out of scope for SORKIN-2):**
Does a zero-energy realization of a given causal order exist in Minkowski
spacetime of a given dimension? This is a statement about the structure of the
causal set and the target manifold.

**Accessibility question (in scope for SORKIN-2):**
Given a causal order for which a zero-energy realization is known or expected,
can the historical Bombelli annealer — with its specific energy function, move
set, temperature schedule, and acceptance rule — find a zero-energy
configuration? This is a statement about a particular constructive procedure.

Failure of the annealer to reach zero energy is evidence about the algorithmic
setup. It is not evidence that no zero-energy realization exists.

---

## Definition of a Known-Truth Case

A case may be included in the known-truth matrix as a positive recoverability
case only if all five criteria below are satisfied:

1. **Defined target order.** The causal order matrix `z[i][j]` is fully
   specified and reproducible from a fixed input or construction.

2. **Independent reason for realizability.** There exists a reason, independent
   of the annealer's energy, to believe that a zero-energy realization exists.
   Acceptable reasons include: mathematical proof (chains, antichains),
   explicit construction of ground-truth coordinates (Minkowski sprinklings
   with known seed), or a theorem guaranteeing embeddability.
   The annealer's own output is not an acceptable reason.

3. **Energy-independent truth source.** The truth source must not depend on the
   annealer finding zero energy. Ground-truth coordinates from the sprinkler
   satisfy this: they exist and have zero energy by construction, before any
   annealing.

4. **Delimited allowed interpretation.** Before any run, the allowed
   interpretation of a result must be stated. What the result can and cannot
   be taken as evidence for must be written explicitly.

5. **Explicit forbidden interpretation.** The forbidden interpretation must be
   stated in the row. In particular: failure to recover must not be interpreted
   as non-embeddability, and success must not be generalized beyond the
   specific instance.

Cases that do not satisfy all five criteria may be included as structural
contrast cases, with that role stated explicitly and no positive recoverability
claim attached.

---

## Schema

| column | type | description |
|---|---|---|
| `case_id` | string | unique identifier for the case |
| `family` | string | structural family: chain, antichain, minkowski\_sprinkling, corona |
| `n` | int | number of events |
| `target_dimension` | int | spacetime dimension of the target Minkowski space |
| `truth_status` | string | one of: guaranteed\_embeddable, constructed\_truth\_exists, unknown\_embedding\_status |
| `why_truth_known` | string | independent reason for realizability; must not reference annealer output |
| `source_in_repo` | string | file(s) where the case or its construction exists |
| `annealer_mode` | string | which annealer configuration is used; must be labeled before any run |
| `expected_role` | string | positive\_sanity\_check, positive\_known\_recovered, documented\_failure, structural\_contrast |
| `allowed_target` | string | what the result may be interpreted as evidence for |
| `forbidden_interpretation` | string | what the result must not be taken as evidence for |
| `missing_piece_before_run` | string | infrastructure or decision gaps that must be resolved before running |

---

## Proposed Cases

The rows below are proposals for the first matrix. They are not experimental
results. Each row must be independently validated before it can produce a
result.

### Row A

| field | value |
|---|---|
| `case_id` | chain\_4\_d2 |
| `family` | chain |
| `n` | 4 |
| `target_dimension` | 2 |
| `truth_status` | guaranteed\_embeddable |
| `why_truth_known` | Any finite chain is realizable as a sequence of timelike-ordered events in 1+1 Minkowski spacetime. This is provable without reference to the annealer. |
| `source_in_repo` | `tests/test_causet_invariants.py` (`make_chain` helper); no `.in` file exists |
| `annealer_mode` | historical/default; must be specified explicitly before run |
| `expected_role` | minimal positive sanity check |
| `allowed_target` | `reaches_zero_energy`; later `induced_order == target` once a standalone verifier exists |
| `forbidden_interpretation` | Success does not generalize to other families or dimensions. Failure does not indicate non-embeddability of chains in Minkowski spacetime. No claim about manifoldlikeness. |
| `missing_piece_before_run` | No `.in` file in Pascal format; `make_chain` is a test helper, not an annealer input. An adapter from causal matrix to `.in` format is required. The annealer mode must be recorded. |

### Row B

| field | value |
|---|---|
| `case_id` | antichain\_4\_d2 |
| `family` | antichain |
| `n` | 4 |
| `target_dimension` | 2 |
| `truth_status` | guaranteed\_embeddable |
| `why_truth_known` | Any finite antichain is realizable as a set of spacelike-separated events. In d=2 (1+1D), four events can be placed at equal coordinate time with distinct spatial positions; no two are causally related. |
| `source_in_repo` | `tests/test_causet_invariants.py` (`make_antichain` helper); no `.in` file exists |
| `annealer_mode` | historical/default; must be specified explicitly before run |
| `expected_role` | minimal complementary sanity check |
| `allowed_target` | `reaches_zero_energy`; later `induced_order == target` once a standalone verifier exists |
| `forbidden_interpretation` | Success does not generalize to other families. Failure does not indicate non-embeddability of antichains. No claim about manifoldlikeness. |
| `missing_piece_before_run` | Same as chain\_4\_d2: no `.in` file; adapter required. Also verify that the energy function correctly handles all-spacelike configurations at initialization. |

### Row C

| field | value |
|---|---|
| `case_id` | minkowski\_6\_s1959\_d2 |
| `family` | Minkowski sprinkling |
| `n` | 6 |
| `target_dimension` | 2 |
| `truth_status` | constructed\_truth\_exists |
| `why_truth_known` | Ground-truth coordinates are generated by the sprinkler in `cones.py` with seed 1959 in a 1+1D Minkowski causal diamond. The causal order is induced by those coordinates. The energy is zero at those coordinates by construction, before any annealing. `test_regression.py` confirms the annealer reaches final energy ≈ 1.6e-11 on `benchmarks/tesis_like_6.in` with seed 1959 under the historical default schedule. |
| `source_in_repo` | `benchmarks/tesis_like_6.in`; `test_regression.py` (seed 1959 regression); sprinkler in `cones.py` |
| `annealer_mode` | historical/default (T₀=100, γ=0.9, warmup\_limit=100, anneal\_limit=100, max\_data=35) |
| `expected_role` | positive known-recovered case; baseline for the annealer succeeding |
| `allowed_target` | `reaches_zero_energy` (confirmed by regression); `induced_order == target` once a standalone verifier exists |
| `forbidden_interpretation` | This single case does not establish annealer completeness. Success here is not evidence that the annealer recovers all Minkowski sprinklings. No claim about manifoldlikeness. No generalization to n=12 or other seeds. |
| `missing_piece_before_run` | A standalone binary verifier of `induced_order == target` does not yet exist independently of the energy function. Ground-truth coordinates are regenerated each run, not persisted; reproducibility depends on fixed seed. |

### Row D

| field | value |
|---|---|
| `case_id` | minkowski\_12\_s1962\_d2\_hist |
| `family` | Minkowski sprinkling |
| `n` | 12 |
| `target_dimension` | 2 |
| `truth_status` | constructed\_truth\_exists |
| `why_truth_known` | Ground-truth coordinates are generated by the sprinkler with seed 1962. The causal order and zero-energy realization exist by construction. The question here is whether the historical annealer accesses that realization. `test_regression.py` documents that the historical default schedule does not converge: final energy ≈ 16.8, 35 cooling blocks exhausted. |
| `source_in_repo` | `benchmarks/tesis_like_12.in`; `test_regression.py` (seed 1962 regression, historical schedule) |
| `annealer_mode` | historical/default (T₀=100, γ=0.9) — the schedule under which failure is documented |
| `expected_role` | documented algorithmic failure case; characterizes the gap between truth existence and annealer accessibility |
| `allowed_target` | `characterize_failure_of_historical_schedule`: record final energy, number of cooling blocks, distance to truth if measurable. The target is characterization of the failure, not recovery. |
| `forbidden_interpretation` | Failure to converge is not evidence of non-embeddability. The zero-energy realization exists; the annealer did not access it under this schedule. Do not report this as a negative result about causal set physics. Do not conflate with Row E (tuned schedule). Results from this row and Row E must not be mixed. |
| `missing_piece_before_run` | The exact historical default parameters used in the original 1987 thesis must be explicitly recorded for this run. The schedule must match `test_regression.py` precisely. A decision is needed on whether `warmup_limit` is 10 (diagnostic grid) or 100 (full historical). |

### Row E

| field | value |
|---|---|
| `case_id` | minkowski\_12\_s1962\_d2\_tuned |
| `family` | Minkowski sprinkling |
| `n` | 12 |
| `target_dimension` | 2 |
| `truth_status` | constructed\_truth\_exists |
| `why_truth_known` | Same ground truth as Row D. Same seed, same causal order. The question here is whether a non-historical, empirically adjusted schedule can access the known zero-energy realization. `test_regression.py` documents that a schedule with T₀=180, γ=0.8 achieves final energy ≈ 0.17, a large reduction from the default schedule's ≈ 16.8. |
| `source_in_repo` | `benchmarks/tesis_like_12.in`; `test_regression.py` (tuned schedule test) |
| `annealer_mode` | tuned/non-historical schedule (T₀=180, γ=0.8); this is not the Bombelli 1987 default |
| `expected_role` | improved algorithmic recoverability case under a non-historical schedule |
| `allowed_target` | `near_zero_energy` under the tuned schedule; the threshold for "near zero" must be explicitly defined before the run (suggested: final\_energy ≤ 1.0) |
| `forbidden_interpretation` | This is not the historical Bombelli annealer. Results from this row must not be attributed to the historical schedule. Do not conflate with Row D. Near-zero energy under a tuned schedule is not a claim about the historical annealer's capabilities. No manifoldlikeness claim. |
| `missing_piece_before_run` | The provenance of T₀=180, γ=0.8 must be documented (source: `test_regression.py`). A threshold for "near zero" must be agreed before the run. The exact tuned parameters must be recorded in the row before any result is written. |

### Row F

| field | value |
|---|---|
| `case_id` | corona\_8\_d2 |
| `family` | corona |
| `n` | 8 |
| `target_dimension` | 2 |
| `truth_status` | unknown\_embedding\_status |
| `why_truth_known` | Not applicable. The corona poset is a manifoldlikeness negative control. Whether a zero-energy Minkowski realization exists for this specific corona is not established. This case is included only as a structural contrast. |
| `source_in_repo` | `validation_suite.py` (`generate_corona_poset`) |
| `annealer_mode` | any; must be labeled before run |
| `expected_role` | structural contrast only; marks a case where annealer behavior cannot be attributed to known realizability or non-realizability |
| `allowed_target` | diagnostic contrast only; the result may describe what the annealer does on a non-Minkowski-like structure |
| `forbidden_interpretation` | Failure does not mean non-embeddable. Success does not mean manifoldlike. This case must not appear as a positive recoverability result. It must not be used to claim that the annealer detects manifoldlikeness or its absence. It must not be cited alongside Rows A–E as if it has the same known-truth status. |
| `missing_piece_before_run` | The embedding status of this specific corona must remain marked as unknown. The case must be labeled as contrast, not known-truth positive, in any output file header or column. A decision is needed on whether to include this case in Phase 6 at all given the absence of ground truth. |

---

## Annealer Modes

Three annealer configurations appear in this matrix. They are not interchangeable.
Results from one configuration must not be reported as if they came from another.

**historical/default**

The configuration implemented in `cones.py` as ported from the Bombelli 1987
thesis (Pascal source). Parameters: T₀=100.0, γ=0.9, warmup\_limit=100 (or 10
in diagnostic grids), anneal\_limit=100, max\_data=35, acceptance\_scale=4.0.
The warmup makes unconditional accepts; Phase 2D documents that this destroys
near-truth configurations.

**guarded\_warmup**

The non-destructive warmup variant introduced in Phase 2F. Accept-only-if-ΔE ≤ 0
during warmup. All other parameters unchanged. Phase 2F documents 18/18
small-noise preservation on the tested grid versus 16/18 with the historical
warmup. This is not the historical Bombelli annealer; it is a diagnostic variant.
Any row using this mode must label it explicitly.

**tuned\_schedule**

A non-historical schedule with adjusted T₀ and γ, sourced from
`test_regression.py`. Parameters: T₀=180, γ=0.8. Reduces final energy on
`tesis_like_12.in` (seed 1962) from ≈ 16.8 to ≈ 0.17. This is not the
historical Bombelli annealer. Results under this mode characterize what an
adjusted optimizer can do; they say nothing about the historical annealer.

**Mixing modes invalidates comparison.** A row must record its annealer mode
before the run. Rows from different modes must be kept in separate rows in the
matrix and must not be aggregated into a single result. In particular, Rows D
and E in this matrix use different modes and must be read separately.

---

## Do Not Use Yet

The following artifacts and approaches must not feed SORKIN-2 in the current
phase.

| artifact or approach | reason |
|---|---|
| `legacy/` (all contents) | Provenance archive from exploratory branches. Requires a fresh audit before any claim can be drawn from it. |
| Phase 3 / PySR artifacts | Symbolic regression requires a clean multi-family known-truth dataset. That dataset does not yet exist. |
| Phase 4A / Phase 4B | Archived exploratory phases. Not active. Located in `legacy/`. |
| Phase 5 (seed curve morphology) | Exploratory. Not reactivated without an explicit decision. |
| `validation_suite.generate_kleitman_rothschild` as a positive case | The Kleitman–Rothschild poset is a manifoldlikeness negative control with unknown embedding status. It must not appear in the matrix as a known-truth positive. |
| `validation_suite.generate_corona_poset` as an embeddability claim | The corona has unknown embedding status. It may appear only as structural contrast (see Row F). Do not use it to claim annealer success or failure on an embeddable structure. |
| tuned schedule attributed to the historical annealer | The tuned schedule (T₀=180, γ=0.8) is not the Bombelli 1987 default. Results obtained with it must not be labeled or described as results of the historical annealer. |
| `myrheim_meyer_dimension` as the primary recoverability target | Myrheim-Meyer dimension is an order-theoretic invariant. It does not verify that the annealer recovered the target causal order. It may be a secondary diagnostic but must not replace `reaches_zero_energy` or `induced_order == target` as the primary target. |
| low final energy alone as a recoverability claim | Final energy near zero is a necessary condition for exact relation recovery, not a sufficient one. A standalone verifier of `induced_order == target` does not yet exist. Until it does, `reaches_zero_energy` is the primary target, and its limitations must be stated. |
| KAN, GAM, or any symbolic model | Deferred until a clean multi-family known-truth dataset exists with stable provenance. The current data are not sufficient. |

---

## Graphical outputs policy

Every accepted SORKIN-2 known-truth run must produce permanent graphical
outputs, not only CSV, JSON, or text artifacts. These figures must be stored in
a traceable path belonging to the run itself, alongside or directly linked to
the machine-readable result files for that run.

Figures do not replace the binary order verifier. `exact_match` remains the
strong criterion for exact recovery. Low final energy is not sufficient by
itself, and plots are diagnostic and interpretive artifacts, not proof of
recovery by themselves.

The minimum figure set for each run is:

1. `target_order_matrix.png`: visualization of the target causal order matrix.

2. `induced_order_matrix.png`: visualization of the causal order induced by
   the recovered final configuration.

3. `order_difference_matrix.png`: discrete difference between target and
   induced order, marking missing relations and extra relations.

4. `recovered_coordinates.png`: scatter or geometry plot of the recovered final
   configuration, when the dimension permits a meaningful visualization.

5. `energy_trace.png`: energy evolution during annealing, if the run records an
   energy trajectory.

Each figure must be traceable to the same provenance record as the numerical
result: `case_id`, `annealer_mode`, seed, input `.in` file, code used,
timestamp or `run_id`, and the JSON or CSV result artifacts from the same run.

Interpretation is constrained. A visually clean figure does not imply general
embeddability. A visually poor figure does not imply non-embeddability. Any
visual mismatch must be confirmed with `compare_causal_orders`. The figures
exist for human audit and diagnosis of failure modes.

---

## Missing Infrastructure

The following gaps must be resolved before any row in this matrix can produce a
verified SORKIN-2 result.

**Standalone induced-order verifier: IMPLEMENTED.** `validation_suite.py`
now provides an energy-independent binary verifier for recovered coordinates
against a target causal order: `OrderComparison`,
`induced_order_from_coords`, `compare_causal_orders`, and
`verify_recovery`. This closes the verifier gap; remaining rows still require
their own inputs and run-mode decisions before any known-truth case is run.

**Adapter from causal matrix to .in format.** Rows A and B (chain, antichain)
require inputs in Pascal `.in` format. The helpers `make_chain` and
`make_antichain` in `tests/test_causet_invariants.py` produce causal matrices
but not `.in` files. An adapter is needed before those rows can be run.

**Explicit annealer-mode record.** Each row must have its `annealer_mode` field
filled before any run. The Phase 2F finding (guarded warmup is best on the
tested grid) does not make guarded warmup the default; `cones.py` still uses
the historical warmup. Any run that departs from the historical default must
document the departure.

**Near-zero threshold definition.** Row E requires a stated threshold for "near
zero" before the run. Using a post-hoc threshold introduces selection bias.

**Ground-truth coordinate persistence.** For Rows C, D, and E, ground-truth
coordinates are regenerated from the sprinkler at each run. They are not stored
as a file. This is reproducible given a fixed seed but means that exact
`induced_order == target` verification requires re-running the sprinkler at the
same seed before comparing.

**Permanent plotting/reporting layer for SORKIN-2 known-truth runs: PARTIALLY
IMPLEMENTED.** Order-matrix diagnostic plots are implemented for
`target_order_matrix.png`, `induced_order_matrix.png`, and
`order_difference_matrix.png`. Coordinate plots and energy traces remain
pending until the first real run-output schema is fixed.

---

## Confirmation

- no files changed except `docs/sorkin2_known_truth_matrix.md`
- no experiments run
- no data generated
- no git add, commit, or push
