# Schedule Seed Stability N=36

Exploratory SORKIN-2 diagnostic for whether the geometric schedule signal
seen in `explore/schedule_probe_n36/` persists across optimizer seeds.

Fixed configuration:

- `N = 36`
- `family = minkowski`
- `d_spacetime = 2`
- `case_seed = 1959`
- `T0 = 100`
- budget `medium_25_25_8`
- unchanged historical Bombelli energy

Matrix:

- `optimizer_seeds = 1959, 1962, 1987, 2001`
- `cooling_factor = 0.5, 0.8, 0.9, 0.95`

The probe uses `ConesSimulator.block_callback` in read-only mode. At each
annealing block it copies `sim.rold`/`sim.xold`, computes the induced causal
order with `validation_suite.induced_order_from_coords`, and compares it to
the known target order with `validation_suite.compare_causal_orders`.

This does not modify `cones.py`, `validation_suite.py`, the historical
annealer, the objective function, the move/acceptance rule, or the internal
dynamics. It is a schedule/seed diagnostic only, not an embeddability claim
and not a physical claim about `gamma` or an `N` transition.

Run from the repository root:

```bash
python3 explore/schedule_seed_stability_n36/run_schedule_seed_stability_n36.py
```

Outputs:

- `schedule_seed_stability_n36.csv`
- `schedule_seed_stability_n36_summary.csv`
- `schedule_seed_stability_n36.md`
- `schedule_seed_stability_n36.svg`
