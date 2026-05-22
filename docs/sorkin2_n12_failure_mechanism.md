# SORKIN-2 N=12 failure mechanism: historical schedule

## Physical question

How does the historical Bombelli schedule fail on the N=12 constructed causal set?
Given that the same input, same energy, and same verifier produce `exact_match = true`
under the tuned schedule, the failure is localized to the historical protocol. This
note characterizes the geometric and causal signature of that failure using the
final coordinate artifacts already generated — no new annealer runs.

---

## Data sources

All analysis in this note is derived from existing result artifacts. No new runs
were executed.

| run | path |
|---|---|
| HIST | `results/sorkin2_known_truth/minkowski_12_s1962_d2_hist/20260522T084722Z_minkowski_12_s1962_d2_hist_seed1962/` |
| TUNED | `results/sorkin2_known_truth/minkowski_12_s1962_d2_tuned/20260522T085340Z_minkowski_12_s1962_d2_tuned_seed1962/` |

Both runs share: `benchmarks/tesis_like_12.in`, seed 1962, `warmup_limit=100`,
`anneal_limit=100`, `max_data=35`, `acceptance_scale=4.0`. Only T₀ and γ differ.

---

## Causal diagnostic

### Summary

| | HIST | TUNED |
|---|---|---|
| target\_relations | 19 | 19 |
| induced\_relations | 1 | 19 |
| missing\_relations | 19 | 0 |
| extra\_relations | 1 | 0 |
| exact\_match | false | **true** |

### HIST: structure of missing relations

All 19 target relations are missing in HIST. The single induced relation, `(10,11)`,
is spurious — it is not in the target order.

The 19 missing pairs decompose by failure mode:

**16 of 19: SPACELIKE.** Temporal ordering is correct (ΔR > 0) but spatial
separation dominates: |ΔX|² >> ΔR², so the pair falls outside each other's
lightcone. Examples:

- `(0,5)`: ΔR = +0.51, |ΔX| = 13.1 — events nearly simultaneous, far apart
- `(0,8)`: ΔR = +2.06, |ΔX| = 22.6 — small time gap, large spatial gap
- `(6,9)`: ΔR = +11.66, |ΔX| = 19.7 — even a moderate ΔR is overwhelmed by |ΔX|

**3 of 19: WRONG ORDER.** The event that should be in the future has a smaller R
than its predecessor:

- `(0,10)`: ΔR = −2.66 (event 10 is earlier than event 0 in R)
- `(1,10)`: ΔR = −0.33 (event 10 nearly simultaneous with event 1, but reversed)
- `(2,9)`:  ΔR = −1.90 (event 9 is earlier than event 2 in R)

Events 9 and 10 are misplaced in the temporal ordering. They should be in the
future of several earlier events but end up at R=12.66 and R=3.36 respectively —
both too early relative to some of their required predecessors.

---

## Geometric diagnostic

### Coordinate statistics

| | HIST | TUNED |
|---|---|---|
| R min | 1.000 | 1.000 |
| R max | 30.181 | 91.235 |
| **R span** | **29.2** | **90.2** |
| X[0] span | 31.8 | 96.7 |
| R/X[0] aspect | ≈ 0.92 | ≈ 0.93 |

Both runs have roughly equal R and X[0] spans in ratio, but TUNED is 3× larger in
absolute scale. This matters for causality: a causal relation requires ΔR² ≥ |ΔX|².
In HIST, the absolute R-spread is insufficient to place most target pairs within
the lightcone, even when their temporal ordering is correct.

### R ordering (sorted by R)

**HIST (R ascending):**  
`6(1.00), 7(1.61), 10(3.36), 4(3.63), 1(3.69), 3(5.20), 0(6.03), 5(6.54), 8(8.08), 9(12.66), 2(14.56), 11(30.18)`

**TUNED (R ascending):**  
`2(1.00), 1(2.29), 0(3.75), 3(5.51), 6(5.98), 7(11.48), 4(15.54), 5(33.62), 9(37.75), 10(44.65), 11(51.82), 8(91.24)`

In HIST, event 10 sits at R=3.36, while events 0 (R=6.03) and 1 (R=3.69) need
to precede it — they don't. Events cluster in the R range 1–6 (seven events in
a span of 5 R-units), making spacelike relations unavoidable for short ΔR pairs
with non-negligible |ΔX|.

In TUNED, no such compression occurs. Events are spread to R≈91, and their
temporal ordering is consistent with the target causal structure.

### Warmup observation (both runs)

Both HIST and TUNED share the same warmup behavior:

- Initial energy: 14.407 (same initial configuration at seed 1962)
- Post-warmup energy: 30.501 (same unconditional warmup, 100 steps)

The historical warmup increases the energy by +16.1 from the initial state —
the destructive warmup effect documented in Phase 2D/2F. The annealing phase then
starts from this elevated E ≈ 30.5.

- **HIST** starts annealing at T₀=100. In the Metropolis criterion, T₀ scales
  the acceptance probability for energy-raising moves; lower T₀ means fewer
  uphill moves accepted, limiting reorganization from the post-warmup state.
- **TUNED** starts at T₀=180, giving higher acceptance probability for uphill
  moves and broader exploration before the schedule cools. Note on cooling
  rates: T_k = T₀·γ^k, so γ=0.8 reaches lower temperatures faster than γ=0.9
  within the same 35 cooling blocks. TUNED both starts higher in T₀ and cools
  to a lower final temperature (≈0.33) than HIST (≈2.24) after 35 blocks.
  T₀ and E are in different units; a comparison "T₀ vs E_warmup" is a
  category error and is not made here.

### Coordinate figure

`docs/figures/sorkin2_n12_hist_vs_tuned_coordinates.png`

Side-by-side plot of final event positions in the (X[0], R) plane. Blue arrows:
target relations that are induced. Red dashed arrows: target relations that are
missing. Orange arrows: extra (spurious) relations. The visual compression of the
HIST R-axis relative to TUNED is the geometric signature of the failure.

---

## Mechanistic hypothesis

The data support the following picture, in order of evidential strength:

**1. Low initial temperature as primary candidate.** After warmup, both runs begin
annealing from the same post-warmup state (E≈30.5). HIST uses T₀=100; TUNED uses
T₀=180. T₀ controls the Metropolis acceptance probability for uphill moves: higher
T₀ means more uphill moves accepted, enabling broader exploration of the energy
landscape. With T₀=100, the annealer accepts fewer uphill moves from the outset,
limiting its ability to redistribute 12 events into a configuration that satisfies
all 19 target causal relations. The schedule then exhausts 35 cooling blocks
without reaching the correct basin. The result is a spatially overcrowded,
temporally compressed configuration — near-antichain in induced order.
Note: HIST uses γ=0.9 (slower cooling) while TUNED uses γ=0.8 (faster cooling);
the role of the cooling rate is not separated by the HIST/TUNED comparison alone.

**2. Destructive warmup as a compounding factor.** The warmup is identical in
both runs and is not the differentiating factor — it affects both equally. However,
it sets a difficult starting point. The warmup increases energy 14.4→30.5,
destroying whatever partial structure the initial configuration had. The historical
warmup (unconditional accept) is already documented in Phase 2D/2F as potentially
destroying near-truth configurations. The annealing phase in HIST must recover from
30.5 with T₀=100; in TUNED, from 30.5 with T₀=180.

**3. Near-antichain trap (consequence, not cause).** The HIST final configuration
is effectively an antichain in induced order: 1/19 causal relations, spatially
overcrowded relative to its R-spread. This is the outcome, not the mechanism. The
trap is entered because T₀=100 gives lower acceptance probability for uphill
moves compared to T₀=180, limiting exploration from the post-warmup state. The
schedule then exhausts available temperatures before the events reorganize into
the wider, temporally separated geometry needed to satisfy all 19 target relations.
Note: HIST (γ=0.9) cools more slowly than TUNED (γ=0.8); the near-antichain trap
is not caused by cooling too fast — HIST's schedule is the slower one.

The dominant geometric signal: 16 of 19 failures are SPACELIKE (not WRONG_ORDER),
meaning the annealer preserved approximate temporal ordering in most cases but
failed to spread events far enough in R relative to their spatial separation.
Three WRONG_ORDER failures (for events 9 and 10) indicate additional temporal
displacement errors.

---

## Cautions

- **One case, one seed.** This analysis covers a single run pair: `tesis_like_12.in`,
  seed 1962. No statement is made about the mechanism of failure in other seeds,
  other N, or other causal set families.

- **Mechanism not demonstrated.** The hypothesis is supported by the geometric
  comparison of two runs but is not confirmed by a controlled single-variable
  experiment. Confirming that T₀ is the primary cause would require holding γ
  fixed and varying only T₀, or vice versa. That experiment has not been done.

- **The X[1] component is not plotted.** The coordinate figure shows (X[0], R)
  only. The lightcone condition involves both X[0] and X[1]. Some pairs that
  appear near-causal in the figure may be spacelike due to the X[1] component.
  Quantitative conclusions use the full |ΔX|² = ΔX[0]² + ΔX[1]².

- **No manifoldlikeness claim.** The analysis of final coordinates does not
  constitute a claim about manifoldlikeness, embeddability, or any physical
  property of the causal set beyond what is directly measured.

- **The tuned result is not generalizable.** TUNED's success (exact_match = true)
  means the energy function permits recovery in this case; it does not imply
  general recovery across all N=12 sprinklings.

---

## Block-trace comparison: gamma=0.9 vs gamma=0.8

This comparison isolates the cooling factor at fixed input, seed, warmup, anneal
block length, initial temperature, acceptance scale, and verifier:

| run | path |
|---|---|
| HIST: T0=100, gamma=0.9 | `results/sorkin2_known_truth/minkowski_12_s1962_d2_hist/20260522T103058Z_minkowski_12_s1962_d2_hist_seed1962/` |
| GAMMA CONTROL: T0=100, gamma=0.8 | `results/sorkin2_known_truth/minkowski_12_s1962_d2_T100_g08/20260522T103101Z_minkowski_12_s1962_d2_T100_g08_seed1962/` |

Both runs write real per-block coordinate-derived traces:
`trace.csv`, `trace_energy.png`, and `trace_relations.png`.

| run | blocks | T first | T final | E first | E final | first induced change | first induced=19 | first exact |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| HIST, gamma=0.9 | 35 | 100.000 | 2.781 | 35.544 | 16.809 | 3 | never | never |
| GAMMA CONTROL, gamma=0.8 | 28 | 100.000 | 0.242 | 35.544 | 0.008 | 3 | 28 | 28 |

Induced relation counts by block:

- HIST, gamma=0.9:
  `4, 4, 5, 2, 8, 6, 9, 10, 4, 1, 10, 3, 3, 10, 8, 7, 15, 9, 6, 5, 13, 11, 7, 11, 22, 2, 9, 6, 18, 7, 25, 21, 11, 13, 1`
- GAMMA CONTROL, gamma=0.8:
  `4, 4, 5, 2, 8, 7, 9, 10, 6, 20, 1, 17, 2, 7, 8, 13, 13, 8, 11, 7, 9, 15, 16, 18, 21, 20, 21, 19`

Physical reading:

- gamma=0.9 does not simply fail because it never builds enough induced
  relations. It transiently builds many relations: 15 at block 17, 18 at block
  29, 22 at block 25, 25 at block 31, and 21 at block 32. The failure is that
  those relations are not the target order and are not retained. The run ends
  with only 1 induced relation, 19 missing target relations, and 1 extra relation.
- gamma=0.8 reaches exact recovery only at the final block. It approaches the
  target gradually in the last part of the run: 15, 16, 18, 21, 20, 21, 19 over
  blocks 22-28, with missing relations dropping to 1 at blocks 25-26 and to 0
  at blocks 27-28. Exact match occurs at block 28, when extra relations also
  drop to 0.
- The gamma=0.8 success therefore looks like late commitment after a gradual
  relaxation, not early freezing. The gamma=0.9 failure looks like unstable
  basin wandering: it creates relation-rich configurations, including too many
  induced relations, but does not settle into the exact target order before the
  run terminates.

Caution: this is still one fixed N=12 case and one optimizer seed. It diagnoses
the block-level mechanism for this known-truth run pair only; it is not a claim
about other seeds, N, or causal set families.

---

## Correct vs spurious relations along the trace

The per-block trace separates induced relations into:

- `correct_relations = 19 - missing_relations_count`
- `extra_relations = extra_relations_count`

The derived comparison figure is:

`docs/figures/sorkin2_n12_trace_correct_vs_extra.png`

| run | max induced block | max induced | correct there | extra there | max correct block | max correct | extra there | final induced | final correct | final extra | final exact |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| HIST, gamma=0.9 | 31 | 25 | 9 | 16 | 25 | 12 | 10 | 1 | 0 | 1 | false |
| GAMMA CONTROL, gamma=0.8 | 25 | 21 | 18 | 3 | 27 | 19 | 2 | 19 | 19 | 0 | true |

Physical reading:

- HIST does generate many induced relations, but relation count alone is
  misleading. Its maximum induced block has 25 induced relations, of which only
  9 are target relations and 16 are spurious. Even its best block by correct
  relations reaches only 12/19 and still carries 10 extras. The failure is not
  an inability to make causal links; it is an inability to align those links
  with the target order and then retain them.
- The gamma=0.8 control succeeds through both mechanisms: it first raises the
  correct relation count to 18/19 by block 25 and 19/19 by block 27, then removes
  the remaining extras by block 28. Block 27 already has all target relations,
  but exact recovery is still false because 2 spurious relations remain.
- Therefore the final success is not only "maximizing correct relations" and not
  only "eliminating extras." It requires both: all 19 target relations present
  and all spurious relations absent. In this run, the last step is primarily the
  cleanup of extras, from 2 to 0, while preserving all 19 correct relations.

Caution: this remains a single known-truth N=12 case at optimizer seed 1962. It
is a diagnostic of recoverability dynamics for this run pair, not a general
statement about the annealer, embeddability, or manifoldlikeness.

---

## Anatomy of spurious relations

The trace block index is 1-based: the first data row in `trace.csv` is
`block=1`. The current traced artifacts contain per-block counts but not
per-block coordinate snapshots or pair lists. Therefore concrete `extra_pairs`
can be reconstructed only for final states, where `annealer_output.txt`,
`result.json`, and the final order-matrix figures record the recovered geometry
and final comparison.

Derived figure:

`docs/figures/sorkin2_n12_spurious_pairs_anatomy.png`

| run/block | source available | induced | correct | extra | exact | concrete extra pairs |
|---|---|---:|---:|---:|---|---|
| HIST block 31 | counts only | 25 | 9 | 16 | false | not reconstructible from saved artifacts |
| HIST block 25 | counts only | 22 | 12 | 10 | false | not reconstructible from saved artifacts |
| HIST block 35 | final coordinates + final comparison | 1 | 0 | 1 | false | `(10, 11)` |
| GAMMA CONTROL block 27 | counts only | 21 | 19 | 2 | false | not reconstructible from saved artifacts |
| GAMMA CONTROL block 28 | final coordinates + final comparison | 19 | 19 | 0 | true | none |

For the reconstructible HIST final extra, using the final recovered coordinates:

| pair | Delta R | |Delta X| | interval DeltaR^2 - |DeltaX|^2 | classification |
|---|---:|---:|---:|---|
| `(10, 11)` | 26.818 | 23.550 | 164.592 | `EXTRA_DEEP_TIMELIKE` |

The diagnostic near-null threshold used here is `|interval| < 1e-6`. This is a
numerical diagnostic threshold only, not a universal physical definition.

Physical reading:

- The final HIST extra `(10, 11)` is not near-null. It is a deep timelike
  relation in the recovered geometry, but it is absent from the target order.
  The final spurious relation is therefore geometrically decisive, not an
  ambiguous boundary relation.
- The existing artifacts do not identify which 16 extra pairs occur at HIST
  block 31, which 10 occur at HIST block 25, or which 2 occur at GAMMA CONTROL
  block 27. Consequently, the current saved data cannot decide whether the
  intermediate spurious relations concentrate on a few hub events or form a
  broad overconnection pattern.
- The GAMMA CONTROL final state has no extras. Since block 27 had 19 correct
  and 2 extra relations, and block 28 has 19 correct and 0 extras, the final
  transition is a cleanup of two spurious relations while preserving every
  target relation. The identities and geometry of those two purged relations
  were not saved in the current artifacts.

To answer the hub/topology question at intermediate blocks without rerunning,
future traced artifacts would need to persist either per-block `rold`/`xold`
snapshots or per-block `correct_pairs`, `missing_pairs`, and `extra_pairs`.

### Required next observable: pair-level block snapshots

The current N=12 traced runs cannot support pair-level anatomy of blocks 25,
31, or 27 because they predate pair-level block snapshots. Future traced reruns
must include two additional JSONL artifacts:

- `trace_pairs.jsonl`: one JSON object per block with `correct_pairs`,
  `missing_pairs`, and `extra_pairs`, plus the corresponding relation counts.
- `trace_coordinates.jsonl`: one JSON object per block with the recovered
  `R` and `X` coordinates for every event.

These artifacts are sufficient to reconstruct the induced order and geometry at
each block without using PNGs and without rerunning the annealer.

### Pair-level block snapshots

The two discriminating N=12, seed 1962 runs were repeated with pair-level block
snapshots:

| run | path |
|---|---|
| HIST, gamma=0.9 | `results/sorkin2_known_truth/minkowski_12_s1962_d2_hist/20260522T111350Z_minkowski_12_s1962_d2_hist_seed1962/` |
| GAMMA CONTROL, gamma=0.8 | `results/sorkin2_known_truth/minkowski_12_s1962_d2_T100_g08/20260522T111357Z_minkowski_12_s1962_d2_T100_g08_seed1962/` |

Both runs now include `trace_pairs.jsonl` and `trace_coordinates.jsonl`, with one
line per cooling block.

Derived figure:

`docs/figures/sorkin2_n12_spurious_pairs_anatomy_snapshots.png`

Key blocks:

| run/block | role | induced | correct | missing | extra | exact |
|---|---|---:|---:|---:|---:|---|
| HIST block 31 | maximum induced relations | 25 | 9 | 10 | 16 | false |
| HIST block 25 | maximum correct relations | 22 | 12 | 7 | 10 | false |
| HIST block 35 | final | 1 | 0 | 19 | 1 | false |
| GAMMA CONTROL block 27 | previous block before exact match | 21 | 19 | 0 | 2 | false |
| GAMMA CONTROL block 28 | first exact match and final | 19 | 19 | 0 | 0 | true |

Concrete extra pairs:

- HIST block 31:
  `(0,1), (0,9), (2,8), (2,10), (3,8), (3,10), (4,8), (4,10), (5,10), (5,11), (6,8), (6,10), (7,8), (7,10), (8,10), (8,11)`
- HIST block 25:
  `(0,9), (1,2), (1,6), (1,7), (1,9), (3,5), (3,8), (4,6), (5,9), (8,9)`
- HIST block 35:
  `(10,11)`
- GAMMA CONTROL block 27:
  `(1,9), (6,10)`
- GAMMA CONTROL block 28:
  none

Geometric classification of extras uses the diagnostic threshold
`near_null` if `|DeltaR^2 - |DeltaX|^2| < 1e-6`; all extras in these key blocks
are `deep_timelike`.

| run/block | pair | Delta R | |Delta X| | interval | class |
|---|---|---:|---:|---:|---|
| HIST 31 | `(0,1)` | 41.889 | 29.169 | 903.861 | deep_timelike |
| HIST 31 | `(0,9)` | 57.941 | 57.732 | 24.214 | deep_timelike |
| HIST 31 | `(2,8)` | 21.311 | 9.613 | 361.731 | deep_timelike |
| HIST 31 | `(2,10)` | 39.305 | 1.683 | 1542.081 | deep_timelike |
| HIST 31 | `(3,8)` | 23.307 | 16.538 | 269.682 | deep_timelike |
| HIST 31 | `(3,10)` | 41.301 | 10.339 | 1598.902 | deep_timelike |
| HIST 31 | `(4,8)` | 19.530 | 19.254 | 10.703 | deep_timelike |
| HIST 31 | `(4,10)` | 37.524 | 15.086 | 1180.506 | deep_timelike |
| HIST 31 | `(5,10)` | 36.378 | 6.762 | 1277.633 | deep_timelike |
| HIST 31 | `(5,11)` | 43.011 | 17.327 | 1549.750 | deep_timelike |
| HIST 31 | `(6,8)` | 21.284 | 20.189 | 45.402 | deep_timelike |
| HIST 31 | `(6,10)` | 39.278 | 11.890 | 1401.418 | deep_timelike |
| HIST 31 | `(7,8)` | 20.027 | 3.912 | 385.759 | deep_timelike |
| HIST 31 | `(7,10)` | 38.021 | 12.155 | 1297.885 | deep_timelike |
| HIST 31 | `(8,10)` | 17.995 | 8.361 | 253.906 | deep_timelike |
| HIST 31 | `(8,11)` | 24.628 | 24.461 | 8.214 | deep_timelike |
| HIST 25 | `(0,9)` | 32.078 | 10.854 | 911.191 | deep_timelike |
| HIST 25 | `(1,2)` | 9.399 | 3.157 | 78.371 | deep_timelike |
| HIST 25 | `(1,6)` | 10.998 | 3.863 | 106.035 | deep_timelike |
| HIST 25 | `(1,7)` | 8.180 | 7.235 | 14.580 | deep_timelike |
| HIST 25 | `(1,9)` | 33.624 | 19.848 | 736.645 | deep_timelike |
| HIST 25 | `(3,5)` | 6.115 | 5.860 | 3.057 | deep_timelike |
| HIST 25 | `(3,8)` | 8.562 | 7.154 | 22.123 | deep_timelike |
| HIST 25 | `(4,6)` | 9.347 | 2.865 | 79.167 | deep_timelike |
| HIST 25 | `(5,9)` | 24.293 | 13.349 | 411.937 | deep_timelike |
| HIST 25 | `(8,9)` | 21.846 | 14.343 | 271.520 | deep_timelike |
| HIST 35 | `(10,11)` | 26.818 | 23.550 | 164.596 | deep_timelike |
| GAMMA CONTROL 27 | `(1,9)` | 17.810 | 17.641 | 5.992 | deep_timelike |
| GAMMA CONTROL 27 | `(6,10)` | 13.555 | 13.541 | 0.358 | deep_timelike |

Physical reading:

- HIST block 31 is not globally random overconnection. The extras concentrate
  strongly on event 8 and event 10, each appearing in 7 of the 16 extra pairs.
  This is a hub-like spurious structure.
- HIST block 25 is less dominated by a single late hub, but event 1 and event 9
  each appear in 4 of 10 extras. The best-correct block still carries substantial
  spurious hub structure.
- The extras are not near-null ambiguities. Even the smallest key-block interval
  in GAMMA CONTROL block 27, `(6,10)`, has interval about 0.358, well above the
  diagnostic `1e-6` threshold.
- GAMMA CONTROL reaches all 19 correct relations at block 27 but still has two
  deep-timelike extras, `(1,9)` and `(6,10)`. Exact match at block 28 is achieved
  by eliminating those extras while preserving all 19 correct relations.

---

## Next physical test (not yet executed)

To discriminate between "insufficient T₀" and "aggressive cooling" as the primary
bottleneck, the minimal test is:

> **Run the same N=12, seed 1962 case at T₀=180, γ=0.9 (historical cooling rate).**

This isolates T₀ from γ. If T₀=180 + γ=0.9 also recovers the case, the bottleneck
was T₀ alone. If it fails but T₀=180 + γ=0.8 succeeds, both T₀ and γ contribute.
If T₀=100 + γ=0.8 is also tested, the γ contribution can be isolated.

This is a 2×2 schedule grid (T₀ ∈ {100, 180} × γ ∈ {0.8, 0.9}) on a single fixed
case — not a sweep. The result determines whether the historical failure is
primarily a temperature-scale problem, a cooling-rate problem, or a joint one.

The two variable-isolation runs are now being executed as separate SORKIN-2
harness cases:

- `minkowski_12_s1962_d2_T180_g09` — T₀=180, γ=0.9 (high start, historical rate)
- `minkowski_12_s1962_d2_T100_g08` — T₀=100, γ=0.8 (historical start, faster rate)

A separate test: **guarded warmup at historical T₀=100.** Phase 2F shows that the
guarded warmup (accept only ΔE ≤ 0 during warmup) preserves near-truth
configurations on small grids. If the historical schedule with guarded warmup
recovers the N=12 case, the warmup destruction is the primary bottleneck, not T₀.
This test remains pending.

---

## Next needed observable: annealing trace

The 2×2 grid localizes the dominant bottleneck to γ (cooling depth), but it does
not characterize **when** the annealer commits to the wrong basin. Three hypotheses
remain open:

1. **Early freezing.** The annealer enters a near-antichain configuration in the
   first few cooling blocks and the schedule never offers enough thermal fluctuation
   to escape. Energy drops sharply after block 1–5 and stays flat.

2. **Late commitment.** Energy remains high or oscillating through most of the
   schedule, then drops to a local minimum in the final blocks — too late to
   reorganize into the correct geometry.

3. **Gradual drift into a wrong basin.** Energy decreases monotonically but
   converges to a non-zero local minimum, with no sharp transition.

To distinguish these, the required observable is the **per-block energy trajectory**
for each of the four 2×2 cases. The SORKIN-2 harness now writes this as
`trace.csv` (block, temperature, energy) and two PNG figures:

- `trace_energy.png`: energy per cooling block (plain energy trace)
- `trace_relations.png`: energy trace (top panel) + final-block relation bar chart
  (bottom panel; induced/missing/extra vs target)

These are written automatically for every new run. Existing runs (pre-trace) do not
have these files.

**Current trace status.** `cones.py` now exposes a no-op-by-default block callback
after each `self.data.append()` call in `ConesSimulator.anneal()`. The SORKIN-2
known-truth harness uses it to snapshot `rold`/`xold` and populate
`induced_relations`, `missing_relations_count`, `extra_relations_count`, and
`exact_match` at every block. This instrumentation does not change energy,
move set, acceptance rule, schedules, or verifier.

---

*Coordinate analysis derived from existing annealer\_output.txt files; no new annealer runs for this analysis.*
*Two mechanism runs (T180\_g09, T100\_g08) added to harness after initial drafting.*
*Trace infrastructure section added 2026-05-22.*
*Date: 2026-05-22.*

---

## Mechanism summary after pair-level snapshots

The central result is that HIST gamma=0.9 does not fail by absence of induced
relations. It fails by forming deep timelike relations that are geometrically
real in the recovered coordinates but topologically wrong relative to the target
order.

Evidence:

- HIST block 31 has 25 induced relations: 9 correct and 16 extra. The spurious
  relations concentrate on events 8 and 10, each appearing in 7 of the 16 extras.
- HIST block 25 has 22 induced relations: 12 correct and 10 extra. The spurious
  relations concentrate on events 1 and 9, each appearing in 4 of the 10 extras.
- All diagnosed extras in the key blocks are `deep_timelike`, not `near_null`,
  under the diagnostic threshold `|DeltaR^2 - |DeltaX|^2| < 1e-6`.

The gamma=0.8 contrast isolates the final recovery step. At block 27, the run
already has all 19 target relations but still has two extras, `(1,9)` and
`(6,10)`. At block 28, it still has 19 correct relations, has 0 extras, and
`exact_match=true`. The last step is therefore a purge of spurious relations,
not a gain of missing correct relations.

Physical interpretation: the historical schedule visits causally dense
geometries, but those geometries realize the wrong suborder. The gamma=0.8 run
permits a final cleanup of spurious causal links while preserving the target
order. In this case, the observed mechanism is selection and purge of the
correct suborder, not mere creation of causal relations.

Cautions: this is one known-truth N=12 case and one optimizer seed. It should
not be generalized to larger N, other seeds, or other families. It is not a
claim of manifoldlikeness, and it is not yet a universal mechanism for the
annealer.

Next question, Option B: classify the 19 target relations by structural role
(covers, transitive relations, long relations, short relations) to determine
which kinds stabilize first and which kinds compete with spurious extras.

---

## Target relation classes and recovery dynamics

The N=12 target order contains 19 relations:

- 17 `COVER` relations.
- 2 `TRANSITIVE` relations: `(0,8)` and `(1,8)`, both mediated by event 5.

The complete per-relation table is stored in:

`docs/tables/sorkin2_n12_target_relation_recovery.csv`

The diagnostic target-class figure is:

`docs/figures/sorkin2_n12_target_relation_classes.png`

Recovery by structural class, using `K=5` final blocks for stability:

| class | count | final HIST | final GAMMA CONTROL | stable last5 HIST | stable last5 GAMMA CONTROL |
|---|---:|---:|---:|---:|---:|
| COVER | 17 | 0 | 17 | 0 | 10 |
| TRANSITIVE | 2 | 0 | 2 | 0 | 0 |

By diagnostic index-distance class:

| index-distance class | count | final HIST | final GAMMA CONTROL | stable last5 HIST | stable last5 GAMMA CONTROL |
|---|---:|---:|---:|---:|---:|
| short (`j-i <= 4`) | 5 | 0 | 5 | 0 | 2 |
| medium (`5 <= j-i <= 8`) | 9 | 0 | 9 | 0 | 6 |
| long (`j-i >= 9`) | 5 | 0 | 5 | 0 | 2 |

The index-distance classes are only label-diagnostic; they are not invariant
physical observables.

Physical reading:

- HIST does not retain either covers or transitive target relations in the final
  state. Its final failure is not selective by cover/transitive type; all 19
  target relations are absent.
- GAMMA CONTROL recovers all covers and both transitive relations in the final
  state. However, only 10 of 17 covers are stable across the final 5 blocks, and
  neither transitive relation is stable across all final 5 blocks. The final
  exact match is therefore a late stabilization event, not a long-standing
  fixed target order throughout the end of the run.
- The stable final relations in GAMMA CONTROL are mostly cover relations and are
  concentrated in the medium index-distance class.

Competition with critical extras:

- The GAMMA CONTROL block-27 extra `(1,9)` connects two structurally central
  endpoints: event 1 has target out-degree 4, and event 9 has target in-degree
  5. The endpoints touch 9 target relations total: 8 covers and 1 transitive
  relation. This extra competes with a high-degree source/sink pattern, not with
  an isolated pair.
- The GAMMA CONTROL block-27 extra `(6,10)` connects lower-degree endpoints:
  event 6 has target out-degree 2, and event 10 has target in-degree 2. Its
  endpoints touch 4 target relations, all covers.
- HIST block-31 hub event 8 is the target of 3 relations, including both
  transitive target relations `(0,8)` and `(1,8)`. HIST block-31 hub event 10
  is the target of 2 cover relations. The spurious hub structure therefore
  overlaps both transitive structure around event 8 and cover structure around
  event 10.
- HIST block-25 hubs event 1 and event 9 are both high-degree target nodes.
  Event 1 is a source for 4 target relations, and event 9 is a sink for 5 cover
  relations. The best-correct HIST block is still contaminated by extras attached
  to structurally central endpoints.

Interpretation: the historical failure is not primarily a failure to create
covers, transitive relations, short relations, or long relations as separate
classes. It is a failure to select the correct endpoint topology among
structurally central nodes. The gamma=0.8 run succeeds only after preserving all
target relations while purging extras that compete with high-degree target
endpoints.

---

## N=12 mechanism synthesis

Minimal experimental result:

- HIST gamma=0.9 fails.
- GAMMA CONTROL gamma=0.8 passes.
- Both use the same input, same optimizer seed, same verifier, same energy, same
  move set, and same acceptance rule. The discriminating schedule parameter is
  the cooling factor.

What the failure is not:

- It is not absence of induced relations.
- It is not near-null ambiguity: all diagnosed extras in the key blocks are
  `deep_timelike`.
- It is not a simple cover/transitive preference. The target has 17 covers and
  2 transitives; HIST retains neither class in the final state, while GAMMA
  CONTROL recovers both classes in the final state.
- It is not evidence of non-embeddability. A failed annealer run is an
  accessibility/recoverability statement, not an existence statement.

What the failure appears to be:

- HIST visits causally dense geometries.
- Many induced relations are deep timelike in the recovered geometry but
  topologically wrong relative to the target order.
- Spurious hubs appear on structurally central endpoints.
- GAMMA CONTROL conserves the correct relations and purges the final extras at
  the end of the run.

Quantitative evidence:

- HIST block 31: 25 induced, 9 correct, 16 extra.
- HIST block 25: 22 induced, 12 correct, 10 extra.
- HIST final: 1 induced, 0 correct, 1 extra.
- GAMMA CONTROL block 27: 21 induced, 19 correct, 2 extra.
- GAMMA CONTROL block 28: 19 induced, 19 correct, 0 extra.
- All diagnosed extras are `deep_timelike`, not `near_null`.
- Target structure: 17 cover relations and 2 transitive relations.

Endpoint topology:

- HIST block 31 hubs: events 8 and 10.
- HIST block 25 hubs: events 1 and 9.
- GAMMA CONTROL final purged extras: `(1,9)` and `(6,10)`.
- `(1,9)` is structurally central: its endpoints touch 9 target relations.

Physical interpretation: the observed barrier is a problem of endpoint and
topological selection inside the causal-order landscape, not merely a problem of
creating causal relations.

Cautions:

- One causal set.
- One optimizer seed.
- No universal mechanism is claimed.
- No manifoldlikeness claim is made.
- No extrapolation to larger N is justified yet.

Next physical question: now that the N=12 mechanism is understood for this seed,
the next real question is whether this endpoint-hub failure is seed-specific,
N=12-specific, or a generic failure mode of the historical schedule.
