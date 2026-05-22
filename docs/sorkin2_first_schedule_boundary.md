# SORKIN-2 first schedule boundary

## Physical question

Is the historical failure of the Bombelli annealer on the N=12 constructed causal
set a geometric or energetic barrier — intrinsic to the causal set or the Bombelli
energy function — or is it a barrier of the historical annealing protocol?

The experiment below tests this by running the same input, same seed, same energy
function, and same verifier under two different protocols: the historical Bombelli
schedule and a non-historical tuned schedule. Only the schedule parameters differ.

---

## Results

| case_id | N | family | annealer_mode | T₀ | γ | target | induced | final_energy | exact_match |
|---|---|---|---|---|---|---|---|---|---|
| chain_4_d2 | 4 | chain | historical/default | 100.0 | 0.9 | 6 | 6 | 0.0 | true |
| antichain_4_d2 | 4 | antichain | historical/default | 100.0 | 0.9 | 0 | 0 | 4.757e-02 | true |
| minkowski_6\_s1959\_d2 | 6 | Minkowski sprinkling | historical/default | 100.0 | 0.9 | 4 | 4 | 1.645e-11 | true |
| minkowski_12\_s1962\_d2\_hist | 12 | Minkowski sprinkling | historical/default | 100.0 | 0.9 | 19 | 1 | 1.681e+01 | **false** |
| minkowski_12\_s1962\_d2\_tuned | 12 | Minkowski sprinkling | tuned/non-historical | 180.0 | 0.8 | 19 | 19 | 6.919e-06 | **true** |

Run paths:

- `results/sorkin2_known_truth/chain_4_d2/20260522T083444Z_chain_4_d2_seed1959/`
- `results/sorkin2_known_truth/antichain_4_d2/20260522T084014Z_antichain_4_d2_seed1959/`
- `results/sorkin2_known_truth/minkowski_6_s1959_d2/20260522T084455Z_minkowski_6_s1959_d2_seed1959/`
- `results/sorkin2_known_truth/minkowski_12_s1962_d2_hist/20260522T084722Z_minkowski_12_s1962_d2_hist_seed1962/`
- `results/sorkin2_known_truth/minkowski_12_s1962_d2_tuned/20260522T085340Z_minkowski_12_s1962_d2_tuned_seed1962/`

The hist and tuned rows for N=12 are two separate measurements. They must not be
aggregated or compared without explicit labeling of which protocol produced which result.

---

## Main interpretation

**The historical schedule fails on N=12.** Under T₀=100, γ=0.9, the annealer
induces only 1 of 19 target causal relations, with 19 missing and 1 spurious,
exhausting all 35 cooling blocks at final energy ≈ 16.81. `exact_match = false`.

**The tuned schedule recovers the same N=12.** Under T₀=180, γ=0.8, with the
same input file, same seed 1962, same Bombelli energy function, and same verifier,
the annealer induces all 19 of 19 target causal relations, with zero missing and
zero extra, at final energy ≈ 6.9e-6. `exact_match = true`.

**The boundary belongs to the historical protocol, not to the causal set or the
energy.** The only variable that differs between these two runs is the annealing
schedule. The Bombelli energy function is identical in both cases. The ground-truth
realization — the sprinkled coordinates at seed 1962 — is the same; it has zero
energy by construction and is accessible to the tuned protocol. The historical
protocol fails to access it.

**The energy function permits recovery.** The fact that the tuned schedule reaches
`exact_match = true` using the same Bombelli energy is direct evidence that the
energy landscape contains a path to the ground-truth configuration. The historical
protocol does not find that path. The failure of the historical protocol is not a
failure of the energy function.

**The causal set is not non-embeddable.** It was constructed by sprinkling
ground-truth coordinates in 1+1D Minkowski spacetime. The zero-energy realization
exists before any annealing, independent of the protocol. The historical failure
says nothing about embeddability.

---

## Cautions

- **One N=12, one seed.** This comparison is between two protocols on a single
  instance: `tesis_like_12.in`, seed 1962. No statement is made about recovery
  rates across seeds, input sizes, or causal set families.

- **No generalization to other sprinklings.** The tuned schedule recovers this
  specific constructed case. Whether it recovers other N=12 sprinklings, other
  seeds, or larger sprinklings is not established here.

- **No manifoldlikeness claim.** Success or failure at any case in this table is
  not interpretable as manifoldlikeness or its absence. The Bombelli annealer is
  not a manifoldlikeness detector.

- **Rows D and E are separate results.** `minkowski_12_s1962_d2_hist` and
  `minkowski_12_s1962_d2_tuned` correspond to Rows D and E in
  `docs/sorkin2_known_truth_matrix.md`. They must not be merged, averaged, or
  summarized as if they were variants of the same run. They answer different
  questions about the same causal set.

- **`tuned/non-historical` is not the Bombelli 1987 annealer.** Results from the
  tuned row must not be attributed to the historical protocol or cited as
  performance of the historical annealer. The label `tuned/non-historical` is
  explicit and must be preserved in any downstream reference.

- **`exact_match` is the primary recovery observable.** Low final energy is
  necessary but not sufficient for exact causal order recovery. The antichain row
  illustrates this: `final_energy = 4.757e-02 > 0` yet `exact_match = true`
  because the induced order (empty) matches the target (empty) exactly.

---

## Figures

Each run directory contains three diagnostic figures. Not embedded here; listed
for traceability.

**chain_4_d2** (`results/sorkin2_known_truth/chain_4_d2/20260522T083444Z_chain_4_d2_seed1959/`):
- `target_order_matrix.png`
- `induced_order_matrix.png`
- `order_difference_matrix.png`

**antichain_4_d2** (`results/sorkin2_known_truth/antichain_4_d2/20260522T084014Z_antichain_4_d2_seed1959/`):
- `target_order_matrix.png`
- `induced_order_matrix.png`
- `order_difference_matrix.png`

**minkowski_6_s1959_d2** (`results/sorkin2_known_truth/minkowski_6_s1959_d2/20260522T084455Z_minkowski_6_s1959_d2_seed1959/`):
- `target_order_matrix.png`
- `induced_order_matrix.png`
- `order_difference_matrix.png`

**minkowski_12_s1962_d2_hist** (`results/sorkin2_known_truth/minkowski_12_s1962_d2_hist/20260522T084722Z_minkowski_12_s1962_d2_hist_seed1962/`):
- `target_order_matrix.png`
- `induced_order_matrix.png`
- `order_difference_matrix.png`

**minkowski_12_s1962_d2_tuned** (`results/sorkin2_known_truth/minkowski_12_s1962_d2_tuned/20260522T085340Z_minkowski_12_s1962_d2_tuned_seed1962/`):
- `target_order_matrix.png`
- `induced_order_matrix.png`
- `order_difference_matrix.png`

The `order_difference_matrix.png` figures for the hist and tuned N=12 runs are the
primary diagnostic artifacts for the next phase: they encode which causal relations
were missing or spurious under each protocol and may reveal the geometric signature
of the historical failure.

---

## Next physical question (not yet executed)

Given that the hist/tuned comparison localizes the failure to the schedule, the
question becomes:

> What is the mechanism of the historical failure?

Several non-exclusive hypotheses can be formulated before any new run:

1. **Insufficient initial temperature.** T₀=100 may not provide enough thermal
   energy to escape the initial basin. The annealer may be exploring a local region
   from the start, never reaching configurations near the ground truth.

2. **Overly aggressive cooling.** γ=0.9 with max_data=35 blocks corresponds to a
   fixed cooling budget. If the energy landscape at N=12 requires finer temperature
   resolution in a particular range, the historical schedule may cool through that
   range too fast to find the right basin.

3. **Destructive warmup.** Phase 2D/2F established that the historical unconditional
   warmup can destroy near-truth configurations. If the initial configuration happens
   to be close to ground truth, the warmup may move it away before the annealing
   phase begins.

4. **Trapping in a near-antichain basin.** The induced order under hist is 1 of 19
   relations — almost an antichain. This pattern may indicate that the annealer
   settles into a nearly-spacelike configuration early and cannot escape it under
   the historical cooling rate.

5. **Move-set locality at scale.** The local transposition move set of the
   Bombelli annealer may require many sequential moves to reorder N=12 events
   correctly, and the cooling budget may run out before enough moves accumulate in
   the right direction.

Testing these hypotheses requires controlled comparisons that isolate each variable:
guarded warmup (fixes hypothesis 3), varied T₀ at fixed γ (tests hypothesis 1),
varied γ at fixed T₀ (tests hypothesis 2), energy-trace analysis at N=12 (tests
hypothesis 4). None of those runs are made here.

---

## Note on scope

This document marks the transition from historical reproduction to SORKIN-2
original research. The first four rows of the results table (through N=6 hist) are
reproductions or confirmations of known behavior. The N=12 hist/tuned comparison is
the first controlled experiment designed and executed as part of SORKIN-2: a
variable-isolation test to localize the failure of a historical scientific procedure
to a specific component of that procedure.

The result — that the barrier belongs to the protocol, not the causal set or the
energy — is the first SORKIN-2 finding.

---

*No new runs were executed to produce this document.*
*All numerical values are read from existing result.json files.*
*Date: 2026-05-22.*
