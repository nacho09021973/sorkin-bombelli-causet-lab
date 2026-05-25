# Trajectory Audit N=36

Exploratory SORKIN-2 diagnostic for one fixed known-truth Minkowski case.
It records per-block historical energy and causal-order metrics for:

- `N = 36`
- `family = minkowski`
- `d_spacetime = 2`
- `case_seed = 1959`
- `optimizer_seed = 1987`
- `T0 = 100`
- `gamma = 0.5`
- `h = 1`
- budgets `short_10_10_4` and `medium_25_25_8`

The probe uses `ConesSimulator.block_callback` in read-only mode. At each
annealing block it copies `sim.rold`/`sim.xold`, computes the induced causal
order with `validation_suite.induced_order_from_coords`, and compares it to
the known target order with `validation_suite.compare_causal_orders`.

This does not modify the historical annealer, the objective function, the
selection or acceptance rule, the temperature schedule, coordinates, or causal
order. It is a trajectory diagnostic only, not an embeddability claim and not
a physics claim about `gamma` or an `N` transition.

Run from the repository root:

```bash
python3 explore/trajectory_audit_n36/run_trajectory_audit_n36.py
```

Outputs:

- `trajectory_audit_n36.csv`
- `trajectory_audit_n36.md`
- `trajectory_audit_n36.svg`
