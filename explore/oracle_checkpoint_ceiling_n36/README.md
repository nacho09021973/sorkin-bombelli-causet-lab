# Oracle Checkpoint Ceiling N=36

Post-run SORKIN-2 diagnostic that reads:

`explore/schedule_seed_stability_n36/schedule_seed_stability_n36.csv`

It asks how much causal recoverability is available if an oracle selects the
best observed checkpoint by known-truth causal F1 instead of taking the final
annealer endpoint.

This script does not run the annealer. It groups by:

`optimizer_seed, schedule_label, cooling_factor, budget_label`

Oracle checkpoint tie-break:

1. Maximize `causal_f1`.
2. Maximize `causal_recall`.
3. Minimize `missing_relations_count`.
4. Minimize `extra_relations_count`.
5. Choose the earliest block.

Run from the repository root:

```bash
python3 explore/oracle_checkpoint_ceiling_n36/run_oracle_checkpoint_ceiling_n36.py
```

Outputs:

- `oracle_checkpoint_ceiling_n36.csv`
- `oracle_checkpoint_ceiling_n36_summary.csv`
- `oracle_checkpoint_ceiling_n36.md`
- `oracle_checkpoint_ceiling_n36.svg`

This is an oracle diagnostic with known truth, not a deployable selector for
truth-free cases, not an embeddability claim, and not a physical gamma claim.
