# energy_f1_decoupling_n36

Post-run SORKIN-2 diagnostic.  Part of the N=36 checkpoint analysis chain.

## What this probe does

Reads the already-generated oracle checkpoint ceiling CSV
(`explore/oracle_checkpoint_ceiling_n36/oracle_checkpoint_ceiling_n36.csv`)
and asks whether the mismatch between energy minimisation and causal F1
maximisation can be split into two qualitatively different failure modes.

## Source CSV choice

This probe reads `oracle_checkpoint_ceiling_n36.csv`, not
`checkpoint_selection_n36.csv`.  Reason: the oracle CSV includes `budget_label`
in every row, which is required for the per-group identity key and the
by-schedule aggregation.  The checkpoint-selection CSV does not expose
`budget_label` at the group level.  The oracle CSV contains all columns needed
(`final_energy_eave`, `best_checkpoint_energy_eave`, `min_energy_eave`,
`block_index_final`, `block_index_best_causal_f1`, `block_index_min_energy`,
`best_matches_min_energy`) and is verified to have 16 rows matching the 16
groups in the 4×4 schedule × seed matrix.

## Failure-mode hypotheses

### H2a — over-annealing candidate

`delta_energy_best_minus_final > 0` and `delta_best_minus_final > 0`

The oracle best-F1 checkpoint has **higher energy** than the final endpoint.
The annealer cooled through the good causal region and froze into a
lower-energy but causally worse local minimum.

### H2b — escape / non-convergence candidate

`delta_energy_best_minus_final < 0` and `delta_best_minus_final > 0`

The oracle best-F1 checkpoint has **lower energy** than the final endpoint.
The system visited the good causal region early (at low energy), then moved
away to states with higher energy and worse causal quality.

### inconclusive

All other cases: best == final, or ΔF1 == 0, or ΔE == 0 with ΔF1 > 0.

## Artefacts

| File | Content |
|---|---|
| `run_energy_f1_decoupling_n36.py` | Script (reads oracle CSV; no annealer runs) |
| `energy_f1_decoupling_n36.csv` | Per-group classification table |
| `energy_f1_decoupling_n36_by_schedule.csv` | Aggregates per schedule |
| `energy_f1_decoupling_n36.md` | Human-readable verdict |
| `energy_f1_decoupling_n36.svg` | Scatter: ΔE vs ΔF1, coloured by schedule |

## Interpretation limits

- The oracle best checkpoint is selected using `causal_f1` against the
  known-truth partial order.  The classification is therefore oracular and
  not deployable without ground truth.
- N = 16 groups, 4 seeds per schedule.  The two-regime hypothesis is
  suggestive, not confirmed.
- This diagnostic does not distinguish an energy function that is
  intrinsically misaligned from schedule-specific convergence failure.

## Dependency chain

```
schedule_seed_stability_n36/schedule_seed_stability_n36.csv
    └─ oracle_checkpoint_ceiling_n36/run_oracle_checkpoint_ceiling_n36.py
            └─ oracle_checkpoint_ceiling_n36/oracle_checkpoint_ceiling_n36.csv
                    └─ energy_f1_decoupling_n36/run_energy_f1_decoupling_n36.py  ← this script
```

## Run

```bash
python3 explore/energy_f1_decoupling_n36/run_energy_f1_decoupling_n36.py
```
