# Checkpoint Selection N=36

Post-run SORKIN-2 diagnostic over the existing trajectory matrix:

`explore/schedule_seed_stability_n36/schedule_seed_stability_n36.csv`

This probe does not run the annealer. It reads the per-block trajectory CSV
and compares each trajectory's final checkpoint against the best causal-F1
checkpoint seen during that trajectory.

Best checkpoint tie-break:

1. Maximize `causal_f1`.
2. Maximize `causal_recall`.
3. Minimize `missing_relations_count`.
4. Minimize `extra_relations_count`.
5. Choose the earliest block.

Run from the repository root:

```bash
python3 explore/checkpoint_selection_n36/run_checkpoint_selection_n36.py
```

Outputs:

- `checkpoint_selection_n36.csv`
- `checkpoint_selection_n36_summary.csv`
- `checkpoint_selection_n36.md`
- `checkpoint_selection_n36.svg`

This is a benchmark diagnostic with known truth. It is not an embeddability
claim, not a physical gamma claim, not an N-transition claim, and not a
truth-free deployable selection criterion yet.
