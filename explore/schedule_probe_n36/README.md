# Schedule Probe N=36

Exploratory SORKIN-2 diagnostic comparing native geometric cooling
schedules for one fixed known-truth Minkowski case:

- `N = 36`
- `family = minkowski`
- `d_spacetime = 2`
- `case_seed = 1959`
- `optimizer_seed = 1987`
- `T0 = 100`
- budget `medium_25_25_8`
- unchanged historical Bombelli energy

The only varied parameter is the existing `ConesSimulator.cooling_factor`:

- `gamma_0p5`: `0.5`
- `gamma_0p8`: `0.8`
- `gamma_0p9`: `0.9`
- `gamma_0p95`: `0.95`

The probe uses `ConesSimulator.block_callback` in read-only mode. At each
annealing block it copies `sim.rold`/`sim.xold`, computes the induced causal
order with `validation_suite.induced_order_from_coords`, and compares it to
the known target order with `validation_suite.compare_causal_orders`.

This does not modify `cones.py`, `validation_suite.py`, the historical
annealer, the objective function, the move/acceptance rule, or the internal
dynamics. It is a schedule diagnostic only, not an embeddability claim and
not a physical claim about `gamma` or an `N` transition.

Run from the repository root:

```bash
python3 explore/schedule_probe_n36/run_schedule_probe_n36.py
```

Outputs:

- `schedule_probe_n36.csv`
- `schedule_probe_n36.md`
- `schedule_probe_n36.svg`
