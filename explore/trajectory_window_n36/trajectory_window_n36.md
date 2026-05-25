# Trajectory-window probe  N=36

Post-run SORKIN-2 diagnostic.  Reads the per-block trajectory CSV
(`explore/schedule_seed_stability_n36/schedule_seed_stability_n36.csv`, 128 rows = 16 groups × 8 blocks)
to characterise the causal-F1 profile shape for each schedule × seed.

## Configuration

- Command: `python3 explore/trajectory_window_n36/run_trajectory_window_n36.py`
- Generated at UTC: `2026-05-25T09:40:03+00:00`
- Source CSV: `explore/schedule_seed_stability_n36/schedule_seed_stability_n36.csv`
- Per-seed CSV: `explore/trajectory_window_n36/trajectory_window_n36_per_seed.csv`
- By-schedule CSV: `explore/trajectory_window_n36/trajectory_window_n36_by_schedule.csv`
- SVG: `explore/trajectory_window_n36/trajectory_window_n36.svg`
- This script does not run the annealer.
- Peak identification uses `causal_f1` against known-truth → oracular, not deployable.
- Convergence threshold: `COLD_THRESHOLD = 5.0` (T_final < 5.0 → converged).
- Cold window: `T ∈ [1.0, 5.0]`  (blocks 6–7 for gamma_0p5).

## Budget convergence audit

All schedules start at T = 100 and run 8 blocks.  Final temperatures differ dramatically.

| schedule | γ | T_initial | T_final | converged | avg F1 range | note |
| --- | ---: | ---: | ---: | --- | ---: | --- |
| gamma_0p5 | 0.5 | 100 | 0.78 | **yes** | 0.303 | cold regime reached; causal structure can form |
| gamma_0p8 | 0.8 | 100 | 20.97 | no | 0.221 | still in high-T exploration (T_final ≈ 21); profiles are noise-dominated |
| gamma_0p9 | 0.9 | 100 | 47.83 | no | 0.168 | still in high-T exploration (T_final ≈ 48); profiles are noise-dominated |
| gamma_0p95 | 0.95 | 100 | 69.83 | no | 0.183 | still in high-T exploration (T_final ≈ 70); profiles are noise-dominated |

**Implication**: gamma_0p8, gamma_0p9, and gamma_0p95 do not reach a cold regime in 8 blocks.
The early-block F1 peaks observed for those schedules (labeled H2b in `energy_f1_decoupling_n36`)
are not escapes from causal attractors — they are the system's starting configuration
carried forward through hot, high-acceptance exploration.
The H2a / H2b classification from the previous probe is valid as a phenomenological
description but conflates two distinct causes: true over-annealing (gamma_0p5) and
insufficient cooling budget (the other three).

## gamma_0p5: F1 profile per seed

| seed | blk 1 | blk 2 | blk 3 | blk 4 | blk 5 | blk 6 | blk 7 | blk 8 | peak blk | peak T | peak F1 | final F1 | loss | in window |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 1959 | 0.3713 | 0.2313 | 0.2517 | 0.3347 | 0.4699 | 0.5488 | 0.5308 | 0.4215 | 6 | 3.125 | 0.5488 | 0.4215 | +0.1274 | **yes** |
| 1962 | 0.3534 | 0.2747 | 0.3713 | 0.3485 | 0.5373 | 0.5214 | 0.5389 | 0.6096 | 8 | 0.781 | 0.6096 | 0.6096 | +0.0000 | no |
| 1987 | 0.292 | 0.3412 | 0.3968 | 0.2857 | 0.4767 | 0.5487 | 0.5383 | 0.3732 | 6 | 3.125 | 0.5487 | 0.3732 | +0.1755 | **yes** |
| 2001 | 0.376 | 0.287 | 0.3858 | 0.2905 | 0.506 | 0.5475 | 0.5842 | 0.375 | 7 | 1.562 | 0.5842 | 0.3750 | +0.2092 | **yes** |

## Diagnostic questions

### Q1 — Is the peak-F1 temperature window consistent for gamma_0p5?

**Yes, with conservative support.**  3/4 seeds show peak F1 within T ∈ [1.0, 5.0] (blocks 6–7, i.e. T ≈ 1.56–3.13).  Seeds 1959 and 1987 both peak at block 6 (T = 3.125, F1 ≈ 0.549).  Seed 2001 peaks at block 7 (T = 1.562, F1 = 0.584).  Seed 1962 is inconclusive: the final block is already the best (F1 = 0.610, loss = 0).  The window is T ∈ [1.56, 3.13], corresponding to blocks 6–7 in the gamma_0p5 schedule.

### Q2 — How much F1 is lost by choosing the final endpoint?

For the 3 H2a seeds (those where the final block is not the best):  losses are 0.127366, 0.175476, 0.209229.  Average loss = **0.17069** causal F1 units.  This is the recoverability gap addressable by a stopping criterion targeting T ∈ [1.56, 3.13].

### Q3 — Does a temperature-based stopping criterion look viable for gamma_0p5?

**Preliminary yes, with strong caveats.**  A rule of the form 'stop annealing when T drops below T_stop ≈ 3' would, in this dataset, capture the best-F1 checkpoint for 2/3 H2a seeds (those peaking at block 6).  Seed 2001 peaks at block 7 (T = 1.56), so it would require T_stop ≈ 1.5.  A bracket T_stop ∈ [1.5, 3.5] covers all 3 H2a seeds.

**What this does not tell us:**
- Whether the block-6/7 checkpoint is in the same causal basin for all seeds (configurations might differ even at similar F1).
- Whether this window generalises to N ≠ 36.
- Whether it generalises to seeds outside {1959, 1987, 2001}.
- Whether stopping at T_stop recovers a valid causal realisation or an approximate one.

### Q4 — Are the other three schedules informative for the same question?

**No, not within this budget.**  gamma_0p8 reaches T_final = 21; gamma_0p9 reaches T_final = 48; gamma_0p95 reaches T_final = 70.  None of these enters the cold regime in 8 blocks.  Their F1 profiles have high variance and low range compared to gamma_0p5.  Comparisons across schedules at this budget measure the effect of total cooling, not the causal landscape at low temperature.

### Q5 — What is the right next probe?

The natural successor is a **cross-seed basin consistency probe**:
do the block-6/7 checkpoints for seeds 1959, 1987, and 2001 correspond to the
same causal configuration (same basin), or do they achieve similar F1 via different
causal orderings?

If same basin → the cold window attracts to a single accessible causal structure;
a stopping rule is meaningful.

If different basins → F1 coincidence is not structural; stopping at T_stop selects
different configurations depending on the seed; no reliable stopping rule exists.

This probe requires the causal configuration (list of causal pairs) at each
checkpoint per seed, which is **not** available in the current CSVs.  It would
require a new run with configuration dumps at each block for gamma_0p5.

## Conservative interpretation

This diagnostic uses `causal_f1` against the known-truth partial order to identify
the peak checkpoint.  It is oracular and not deployable without ground truth.

The temperature window T ∈ [1.56, 3.13] (blocks 6–7) is observed in 3 seeds of
gamma_0p5 at N = 36.  This is **suggestive, not confirmatory**:  4 seeds at one N
is insufficient to claim a general stopping rule.

The average endpoint F1 loss of ≈ 0.1707 across the 3 H2a seeds represents
the upper bound on what a correct stopping rule could recover under oracle conditions.
A real stopping rule (without ground truth) would recover less.

The reframe from `energy_f1_decoupling_n36` is confirmed:  the three slow schedules
are not informative about the causal landscape at low temperature within this budget.
The H2b classification for gamma_0p9/0p95 reflects insufficient cooling, not a
distinct physical escape mechanism.

## Guardrails

This is a post-run diagnostic only, over benchmark cases with known truth.
It is not an embeddability claim, not a physical gamma claim, not an N-transition claim,
and not proof of general annealer failure.
It is not a deployable checkpoint-selection criterion.
