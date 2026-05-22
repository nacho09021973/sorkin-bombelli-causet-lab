# SORKIN-2 N=12 topology panel plan

## Physical question

Does historical annealer accessibility at N=12 depend on the topology of the
target poset?

This is an algorithmic recoverability question. It is not an embeddability test,
not a manifoldlikeness detector, and not a claim about physical realizability
beyond the specified known-truth construction of each case.

## Motivation from `minkowski_12_s1962_d2`

The studied N=12 case shows an endpoint/topology failure mode:

- HIST gamma=0.9 fails even though it visits relation-rich configurations.
- The failure is not absence of relations and not near-null ambiguity.
- Key spurious relations are deep timelike in the recovered geometry.
- HIST develops spurious hubs on structurally central endpoints.
- GAMMA CONTROL T100 gamma=0.8 reaches all 19 correct relations before the
  final exact match, then purges the last two extras.

The next question is whether this endpoint-hub failure is specific to the
`tesis_like_12.in` topology, specific to seed 1962, or a broader failure mode of
the historical schedule on N=12 known-truth cases.

## Existing reusable infrastructure

Existing N=12 inputs:

| input | status |
|---|---|
| `benchmarks/tesis_like_12.in` | only existing N=12 `.in` benchmark found; already used by `minkowski_12_s1962_d2_*` harness cases |

Other existing inputs are not N=12:

| input | N | role |
|---|---:|---|
| `benchmarks/tesis_like_6.in` | 6 | small Minkowski benchmark |
| `benchmarks/known_truth/chain_4_d2.in` | 4 | trivial known-truth sanity check |
| `benchmarks/known_truth/antichain_4_d2.in` | 4 | trivial known-truth sanity check |

Reusable generators and metrics already present:

- `validation_suite.sprinkle_minkowski_diamond`: canonical 1+1D and higher-D
  Minkowski diamond sprinkler with ground-truth coordinates.
- `cones.generate_sprinkled_causet`: legacy sprinkler used by the historical
  benchmark path.
- `cones.transitive_reduction`: cover/link extraction.
- `causet_invariants.relation_count`, `link_count`, `height`,
  `antichain_profile`, `chain_counts`, `ordering_fraction`.
- `validation_suite.write_order_matrix_plots`: target order matrix PNGs.

Missing before the panel can be run:

- Additional permanent N=12 `.in` files for the screened sprinkling cases.
- Ground-truth coordinate files for screened sprinkling cases.
- Per-case `topology_summary.json` for screened sprinkling cases.
- Per-case `target_order_matrix.png` for screened sprinkling cases.
- Harness registration of new case IDs after the inputs exist.

## Proposed panel

The panel is deliberately small: six N=12 cases, all interpretable first in
1+1D. The first four are fixed structural anchors. The last two are screened
1+1D sprinklings selected to bracket the observed `s1962` topology without
turning this into a sweep.

| label | proposed case_id | input_file | family | N | target_relations | covers | transitives | max_out_degree | max_in_degree | width_estimate | height_estimate | hub_nodes | reason_for_inclusion |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---|---:|---|---|
| A | `minkowski_12_s1962_d2` | `benchmarks/tesis_like_12.in` | 1+1D Minkowski sprinkling | 12 | 19 | 17 | 2 | 4 | 7 | `max(antichain_profile)` pending recompute in summary | 3 | 1, 9, 11 by high endpoint degree | Baseline studied case with endpoint-hub failure. |
| B | `chain_12_d2` | `benchmarks/known_truth/n12_topology_panel/chain_12_d2.in` | explicit 1+1D construction | 12 | 66 | 11 | 55 | 11 | 11 | 1 | 12 | broad high degree, no branching ambiguity | Input created; topology summary created; target matrix created; not run. Dense total order extreme; tests whether difficulty is caused by relation count alone. |
| C | `antichain_12_d2` | `benchmarks/known_truth/n12_topology_panel/antichain_12_d2.in` | explicit 1+1D construction | 12 | 0 | 0 | 0 | 0 | 0 | 12 | 1 | none | Input created; topology summary created; target matrix created; not run. Sparse extreme; checks all-spacelike recovery and guards against relation-count-only narratives. |
| D | `layered_4_4_4_d2` | `benchmarks/known_truth/n12_topology_panel/layered_4_4_4_d2.in` | explicit 1+1D layered construction | 12 | 48 | 32 | 16 | 8 | 8 | 4 | 3 | all layer nodes have high total degree | Input created; topology summary created; target matrix created; not run. Controlled hub-rich topology with many transitives and fixed width/depth. |
| E | `minkowski_12_sparse_d2_seedTBD` | `benchmarks/known_truth/n12_topology_panel/<case_id>.in` | screened 1+1D Minkowski sprinkling | 12 | pending, target below 19 | pending | pending | pending | pending | pending | pending | pending | Sprinkling with lower relation count and weaker hub structure than s1962. |
| F | `minkowski_12_hub_d2_seedTBD` | `benchmarks/known_truth/n12_topology_panel/<case_id>.in` | screened 1+1D Minkowski sprinkling | 12 | pending, target near or above 19 | pending | pending | pending | pending | pending | pending | pending | Sprinkling with strong endpoint hubs to test whether hub topology predicts historical failure. |

Notes:

- `width_estimate` is initially the maximum canonical level size from
  `causet_invariants.antichain_profile`, not an exact maximum antichain unless
  an exact width helper is added later.
- The hand-constructed cases B-D should include explicit coordinates in their
  provenance so their known-truth status is independent of the annealer.
- Cases E-F require a small selection pass over candidate 1+1D seeds, but that
  selection is a future input-construction step, not an annealer run.

## Topological metrics

Each case should receive a `topology_summary.json` before any annealer run:

| field | definition |
|---|---|
| `case_id` | permanent identifier used by the harness |
| `input_file` | Pascal-format `.in` file |
| `family` | explicit construction, 1+1D Minkowski sprinkling, or structural contrast |
| `N` | number of elements; fixed to 12 for this panel |
| `target_relations` | number of true upper-triangular target entries |
| `covers` | number of transitive-reduction links |
| `transitives` | `target_relations - covers` |
| `max_out_degree` | maximum number of target future relations from any node |
| `max_in_degree` | maximum number of target past relations into any node |
| `degree_profile` | sorted or per-node in/out degree list |
| `width_estimate` | `max(causet_invariants.antichain_profile(z))` initially; exact width pending |
| `height_estimate` | `causet_invariants.height(z)` |
| `hub_nodes` | nodes with high total endpoint degree under a stated threshold, e.g. max-degree nodes or degree >= 75th percentile |
| `chain_counts` | `causet_invariants.chain_counts(z, k_max=4)` |
| `ordering_fraction` | relation density among all unordered pairs |
| `reason_for_inclusion` | the topology contrast the case is meant to test |

Exact width is not required for the first panel. If an exact maximum-antichain
helper is added later, it should be recorded separately from the level-profile
estimate to avoid silently changing the meaning of old summaries.

## Input construction rule

Prefer existing inputs where available. Case A and the explicit anchor cases
B-D now exist as N=12 inputs. The screened sprinkling cases E-F still need new
permanent inputs before any run:

`benchmarks/known_truth/n12_topology_panel/<case_id>.in`

For generated or constructed inputs, also store:

- ground-truth coordinates when they exist;
- `topology_summary.json`;
- `target_order_matrix.png`.

No annealer output belongs in the input directory.

If a generator is added later, it should be a minimal input-construction tool,
not a run harness. Its job is to freeze known-truth N=12 inputs and topology
summaries. It should not run `ConesSimulator`.

## Comparison rule

For each panel case, compare schedules in this order:

1. Run `historical/default` first.
2. Only if historical/default fails exact recovery, run the gamma control
   `T0=100, gamma=0.8`.

The comparison must keep fixed:

- input causal matrix;
- optimizer seed for the paired comparison;
- Bombelli energy;
- move set;
- acceptance rule;
- verifier.

Do not introduce the tuned T180 gamma=0.8 schedule in the first topology panel.
Do not use guarded warmup in this panel.

## Cautions

- This is a small diagnostic panel, not a sweep.
- Do not mix 1+1D and 2+1D cases in the first pass.
- Do not infer manifoldlikeness from success or failure.
- Do not infer non-embeddability from annealer failure.
- Each case is a constructed-truth or explicitly labeled structural diagnostic
  for algorithmic accessibility.
- The screened sprinkling cases should be selected by topology before annealer
  runs, otherwise the panel becomes outcome-selected.

## Next step

Create or locate the missing panel inputs and topology summaries. Do not run the
annealer until the panel inputs, case IDs, and topology summaries are frozen.

## Historical/default results for explicit anchors

The three explicit N=12 topology anchors were run once under
`historical/default` with optimizer seed 1962. No gamma-control, tuned schedule,
guarded warmup, screened sprinkling, multi-seed run, or sweep was run.

Figure:

`docs/figures/sorkin2_n12_topology_panel_historical_results.png`

| case_id | role | run path | target_relations | covers | transitives | width | height | max_out | max_in | final_energy | exact_match | induced_final | correct_final | missing_final | extra_final |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|
| `chain_12_d2` | dense total-order extreme | `results/sorkin2_known_truth/chain_12_d2/20260522T120950Z_chain_12_d2_seed1962/` | 66 | 11 | 55 | 1 | 12 | 11 | 11 | 0.000000 | true | 66 | 66 | 0 | 0 |
| `antichain_12_d2` | sparse all-spacelike extreme | `results/sorkin2_known_truth/antichain_12_d2/20260522T120954Z_antichain_12_d2_seed1962/` | 0 | 0 | 0 | 12 | 1 | 0 | 0 | 0.976436 | true | 0 | 0 | 0 | 0 |
| `layered_4_4_4_d2` | controlled layered hub-rich topology | `results/sorkin2_known_truth/layered_4_4_4_d2/20260522T120957Z_layered_4_4_4_d2_seed1962/` | 48 | 32 | 16 | 4 | 3 | 8 | 8 | 31.915135 | false | 20 | 14 | 34 | 6 |

Conservative reading:

- Historical/default recovers the N=12 chain exactly, despite maximum relation
  density. The historical failure is therefore not caused by N=12 or relation
  count alone.
- Historical/default recovers the N=12 antichain exactly at the order level.
  Its final energy is not zero, so exact order recovery and low energy remain
  distinct diagnostics.
- Historical/default fails the layered 4-4-4 case. The final extras are all
  within-layer relations: 4 extras inside L1 and 2 extras inside L2. The missing
  relations are all required inter-layer relations: 10 missing L0->L1, 11
  missing L0->L2, and 13 missing L1->L2.
- This points away from simple density, width, or height explanations. The
  diagnostic failure is concentrated in a layered topology where the annealer
  must preserve antichains within layers while simultaneously realizing many
  complete inter-layer causal relations.

Next decision:

- If analyzing the layered failure, inspect pair-level anatomy for within-layer
  extras and inter-layer missing relations.
- If all explicit anchors had passed, the next step would have been screened
  sparse/hub sprinkling selection. Because layered fails, the immediate decision
  is whether to analyze that failure before selecting the screened sprinklings.

## Anatomy of historical failure in layered_4_4_4_d2

Run analyzed:

`results/sorkin2_known_truth/layered_4_4_4_d2/20260522T120957Z_layered_4_4_4_d2_seed1962/`

Figure:

`docs/figures/sorkin2_layered_4_4_4_historical_failure_anatomy.png`

Final block summary:

| block | induced | correct | missing | extra | exact_match | final_energy |
|---:|---:|---:|---:|---:|---|---:|
| 35 | 20 | 14 | 34 | 6 | false | 31.915135 |

Final correct pairs:

`(0,5), (0,7), (0,10), (1,7), (1,10), (2,5), (2,7), (2,9), (2,10), (3,7), (3,10), (4,10), (5,10), (6,10)`

Final missing pairs by required layer relation:

| type | count | pairs |
|---|---:|---|
| missing L0->L1 | 10 | `(0,4), (0,6), (1,4), (1,5), (1,6), (2,4), (2,6), (3,4), (3,5), (3,6)` |
| missing L1->L2 | 13 | `(4,8), (4,9), (4,11), (5,8), (5,9), (5,11), (6,8), (6,9), (6,11), (7,8), (7,9), (7,10), (7,11)` |
| missing L0->L2 | 11 | `(0,8), (0,9), (0,11), (1,8), (1,9), (1,11), (2,8), (2,11), (3,8), (3,9), (3,11)` |

Final extras by forbidden within-layer relation:

| type | count | pairs |
|---|---:|---|
| extra within L0 | 0 | none |
| extra within L1 | 4 | `(4,5), (4,7), (5,7), (6,7)` |
| extra within L2 | 2 | `(8,10), (9,10)` |
| extra cross-layer wrong | 0 | none |

Geometry of final extras:

| pair | type | Delta R | |Delta X| | interval | classification |
|---|---|---:|---:|---:|---|
| `(4,5)` | L1->L1 | 66.353 | 54.241 | 1460.644 | deep_timelike |
| `(4,7)` | L1->L1 | 612.160 | 308.739 | 279419.816 | deep_timelike |
| `(5,7)` | L1->L1 | 545.807 | 362.800 | 166281.041 | deep_timelike |
| `(6,7)` | L1->L1 | 610.552 | 250.095 | 310226.660 | deep_timelike |
| `(8,10)` | L2->L2 | 292.987 | 157.004 | 61190.894 | deep_timelike |
| `(9,10)` | L2->L2 | 256.305 | 101.404 | 55409.650 | deep_timelike |

Missing relation causes:

| required layer type | MISSING_SPACELIKE | MISSING_WRONG_ORDER |
|---|---:|---:|
| L0->L1 | 6 | 4 |
| L1->L2 | 5 | 8 |
| L0->L2 | 6 | 5 |

Layer geometry at final block:

| layer | R span | X[0] span | X[1] span | induced internal timelike pairs |
|---|---:|---:|---:|---|
| L0 | 88.995 | 105.536 | 80.324 | none |
| L1 | 612.160 | 80.555 | 353.744 | `(4,5), (4,7), (5,7), (6,7)` |
| L2 | 292.987 | 269.636 | 130.610 | `(8,10), (9,10)` |

Dynamic key blocks:

| diagnostic | block | induced | correct | missing | extra | exact_match | score=correct-extra |
|---|---:|---:|---:|---:|---:|---|---:|
| maximum correct | 14 | 21 | 21 | 27 | 0 | false | 21 |
| minimum missing | 14 | 21 | 21 | 27 | 0 | false | 21 |
| minimum extra | 14 | 21 | 21 | 27 | 0 | false | 21 |
| best score | 14 | 21 | 21 | 27 | 0 | false | 21 |

Interpretation:

- The final failure is a combination of missing required inter-layer relations
  and forbidden intra-layer timelike relations.
- The six final extras are not cross-layer orientation errors. They are all
  forbidden within-layer relations, concentrated in L1 and L2.
- The missing relations are distributed across all required layer connections,
  with the largest count in L1->L2.
- The final geometry strongly violates the intended layered structure: L1 and
  L2 become internally timelike, while many required inter-layer pairs are
  either spacelike or temporally inverted.
- The best block is block 14, not the final block. At block 14 the run has
  21 correct relations and 0 extras, but still misses 27 required relations.
  Later cooling introduces forbidden within-layer extras and loses correct
  inter-layer structure.

Comparison with `minkowski_12_s1962_d2`:

- Both failures involve deep timelike wrong relations rather than near-null
  ambiguity.
- The Minkowski case showed endpoint-hub competition in a sparse sprinkling
  topology. The layered case shows a different but related failure: the annealer
  does not preserve antichain layers while also realizing dense inter-layer
  order.
- The layered failure is therefore not simply the same endpoint-hub mechanism.
  It is a layer-separation failure: within-layer antichains collapse into
  timelike substructure while required cross-layer relations remain incomplete.

## Gamma-control test for layered_4_4_4_d2

A single gamma-control run was executed for `layered_4_4_4_d2`:

`results/sorkin2_known_truth/layered_4_4_4_d2_T100_g08/20260522T121920Z_layered_4_4_4_d2_T100_g08_seed1962/`

Parameters: same input, optimizer seed 1962, same energy, same move set, same
acceptance rule, same verifier, `T0=100`, `gamma=0.8`. No chain/antichain rerun,
T180 run, tuned schedule, guarded warmup, screened sprinkling, multi-seed run,
or sweep was run.

Figure:

`docs/figures/sorkin2_layered_4_4_4_hist_vs_g08.png`

| run | blocks | final_energy | exact_match | induced_final | correct_final | missing_final | extra_final | first exact |
|---|---:|---:|---|---:|---:|---:|---:|---:|
| HIST gamma=0.9 | 35 | 31.915135 | false | 20 | 14 | 34 | 6 | never |
| GAMMA CONTROL gamma=0.8 | 24 | 0.050869 | true | 48 | 48 | 0 | 0 | 24 |

Dynamic comparison:

| run | diagnostic block | induced | correct | missing | extra | exact_match | score=correct-extra |
|---|---:|---:|---:|---:|---:|---|---:|
| HIST best block | 14 | 21 | 21 | 27 | 0 | false | 21 |
| HIST final | 35 | 20 | 14 | 34 | 6 | false | 8 |
| GAMMA CONTROL best/final | 24 | 48 | 48 | 0 | 0 | true | 48 |

Gamma-control final anatomy:

- `missing_pairs`: none.
- `extra_pairs`: none.
- missing by layer type: L0->L1 = 0, L1->L2 = 0, L0->L2 = 0.
- extra by layer type: within L0 = 0, within L1 = 0, within L2 = 0,
  wrong cross-layer = 0.

Interpretation:

- `gamma=0.8` also crosses the layered barrier for this fixed N=12 case and
  seed. The same gamma-control that recovered `minkowski_12_s1962_d2` also
  recovers the controlled layered topology.
- For historical/default, the best block had no extras but still missed 27
  required inter-layer relations; the final state then degraded into both
  missing inter-layer relations and forbidden within-layer extras.
- For gamma-control, exact match appears at the final block 24. This supports
  the interpretation that deeper/faster cooling can complete the late
  organization step: preserve all required inter-layer relations while purging
  within-layer extras.

Caution: this is still one topology and one optimizer seed. It does not prove a
universal gamma effect, nor does it generalize to other N or screened
sprinklings.

## Cross-topology gamma-control result

The first cross-topology N=12 comparison now contains two nontrivial target
orders where the historical/default protocol fails and the explicit
`T100_g08` gamma-control succeeds. Both comparisons use the same input within
each case, optimizer seed 1962, historical energy, move set, acceptance rule,
and exact-order verifier.

| case_id | topology | annealer_mode | T0 | gamma | target_relations | correct_final | missing_final | extra_final | final_energy | exact_match | run_path |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|---|
| `minkowski_12_s1962_d2_hist` | 1+1D sprinkling, endpoint-hub sensitive | historical/default | 100 | 0.9 | 19 | 0 | 19 | 1 | 16.808901 | false | `results/sorkin2_known_truth/minkowski_12_s1962_d2_hist/20260522T111350Z_minkowski_12_s1962_d2_hist_seed1962/` |
| `minkowski_12_s1962_d2_T100_g08` | 1+1D sprinkling, endpoint-hub sensitive | mechanism/T100_g08 | 100 | 0.8 | 19 | 19 | 0 | 0 | 0.008072 | true | `results/sorkin2_known_truth/minkowski_12_s1962_d2_T100_g08/20260522T111357Z_minkowski_12_s1962_d2_T100_g08_seed1962/` |
| `layered_4_4_4_d2` | controlled three-layer poset | historical/default | 100 | 0.9 | 48 | 14 | 34 | 6 | 31.915135 | false | `results/sorkin2_known_truth/layered_4_4_4_d2/20260522T120957Z_layered_4_4_4_d2_seed1962/` |
| `layered_4_4_4_d2_T100_g08` | controlled three-layer poset | mechanism/T100_g08 | 100 | 0.8 | 48 | 48 | 0 | 0 | 0.050869 | true | `results/sorkin2_known_truth/layered_4_4_4_d2_T100_g08/20260522T121920Z_layered_4_4_4_d2_T100_g08_seed1962/` |

Physical interpretation:

- Historical/default fails in two qualitatively different N=12 topologies.
- `gamma=0.8` recovers both cases with the same optimizer seed, 1962.
- This favors the working hypothesis that the observed N=12 boundary is a
  historical cooling-protocol accessibility boundary, not simply a consequence
  of relation density, width, height, or number of target relations.
- The concrete failure mechanism is topology-dependent. In
  `minkowski_12_s1962_d2`, the historical run fails through endpoint-hub
  competition and deep-timelike spurious relations. In `layered_4_4_4_d2`, it
  fails through layer-separation breakdown: missing inter-layer relations plus
  forbidden intra-layer extras.
- The shared feature is not the same final error pattern, but the fact that the
  historical cooling path does not settle into the exact target order while the
  explicit gamma-control does.

Cautions:

- This comparison contains only two nontrivial N=12 topologies.
- It uses one optimizer seed.
- It should not be generalized to larger N yet.
- It is not a manifoldlikeness claim.
- It is not a universal claim about `gamma=0.8`.
- It is not a broad tuned-schedule result; it is only the explicit
  gamma-control comparison `T0=100`, `gamma=0.8`.

Next physical question:

Does the historical fail / `gamma=0.8` pass pattern persist in N=12 sprinklings
selected by topology, especially sparse and hub-rich cases?
