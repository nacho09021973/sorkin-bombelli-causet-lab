# gamma0p5_basin_consistency_n36

Post-run SORKIN-2 probe.  Part of the N=36 checkpoint analysis chain.

## What this probe does

Runs the Bombelli annealer with schedule `gamma_0p5` / `medium_25_25_8`
for N=36, case_seed=1959, four optimizer seeds {1959, 1962, 1987, 2001}.

At each annealing block the full induced causal pair set is captured and
serialized as a `"|"`-separated `"i:j"` string.  After the run, pairwise
Jaccard similarity is computed between:

- **Candidates** — the H2a peak checkpoints identified by `trajectory_window_n36`
  (blocks 6–7 for seeds 1959, 1987, 2001) plus the inconclusive seed 1962 at
  its final block (block 8).
- **Controls** — block 1 (T=100) for all four seeds.

The central question: do the H2a peak checkpoints (which share similar
causal_f1 scalars) also share a common causal configuration (high Jaccard)?

## Why this question matters

`trajectory_window_n36` found that 3/4 seeds under gamma_0p5 achieve their
best causal_f1 at blocks 6–7 (T ≈ 1.56–3.13), with an average endpoint loss
of ≈ 0.171 causal F1 units.  A temperature-based stopping criterion targeting
that window is only meaningful if the configurations are structurally similar
across seeds — i.e., they are in the same causal basin.  If different seeds
peak at similar F1 but different configurations, stopping at T ≈ 3 selects
seed-dependent structures and the criterion is unreliable.

## Reproducibility guardrail

Before any Jaccard interpretation, the script validates F1 and energy values
against `explore/schedule_seed_stability_n36/schedule_seed_stability_n36.csv`
block-by-block (tolerance: |Δ causal_f1| ≤ 1e-6).  If the check fails, the
script writes a warning-only markdown and exits with code 1.

## Checkpoints

| checkpoint | type | seed | block | expected T |
|---|---|---:|---:|---:|
| s1959_b6 | H2a_peak | 1959 | 6 | 3.125 |
| s1987_b6 | H2a_peak | 1987 | 6 | 3.125 |
| s2001_b7 | H2a_peak | 2001 | 7 | 1.562 |
| s1962_b8 | inconclusive | 1962 | 8 | 0.781 |
| s1959_b1 | control_blk1 | 1959 | 1 | 50.000 |
| s1987_b1 | control_blk1 | 1987 | 1 | 50.000 |
| s2001_b1 | control_blk1 | 2001 | 1 | 50.000 |
| s1962_b1 | control_blk1 | 1962 | 1 | 50.000 |

## Jaccard thresholds

| threshold | value | interpretation |
|---|---:|---|
| JACCARD_HIGH | 0.75 | above → basin compartido candidate |
| JACCARD_LOW | 0.50 | below → basins distintos |

## Artefacts

| File | Content |
|---|---|
| `run_gamma0p5_basin_consistency_n36.py` | Script (runs annealer; 4 seeds) |
| `gamma0p5_basin_consistency_n36_trajectory.csv` | 32 rows: 4 seeds × 8 blocks, includes `induced_pairs_str` |
| `gamma0p5_basin_consistency_n36_target_pairs.csv` | 343 rows: ground-truth causal pairs (i, j) |
| `gamma0p5_basin_consistency_n36_jaccard.csv` | 36 rows: upper triangle of 8×8 pairwise Jaccard matrix |
| `gamma0p5_basin_consistency_n36.md` | Human-readable verdict (6 diagnostic questions) |
| `gamma0p5_basin_consistency_n36.svg` | 8×8 Jaccard heatmap (matplotlib imshow) |

## Interpretation limits

- Classification uses `causal_f1` against known truth → oracular, not deployable.
- Sample: 3 H2a seeds, 1 case seed (1959), 1 N (36), 1 schedule (gamma_0p5).
- Jaccard measures induced-pair overlap, not topological or metrical proximity.
- High Jaccard is necessary but not sufficient for a "shared causal basin" claim.

## Dependency chain

```
schedule_seed_stability_n36/schedule_seed_stability_n36.csv   ← reference for validation
    └─ trajectory_window_n36/trajectory_window_n36.md         ← identified H2a peak checkpoints
            └─ gamma0p5_basin_consistency_n36/run_gamma0p5_basin_consistency_n36.py  ← this script
```

## Run

```bash
python3 explore/gamma0p5_basin_consistency_n36/run_gamma0p5_basin_consistency_n36.py
```
