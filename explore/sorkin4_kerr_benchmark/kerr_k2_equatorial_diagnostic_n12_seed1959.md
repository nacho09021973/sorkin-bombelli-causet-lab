# S4-K2 Kerr equatorial diagnostic scaffold

Scope: equatorial kinematic diagnostic only. This does not decide Kerr causal relations.

- `M`: `1.0`
- `theta`: `1.5707963267948966`
- `spins`: `[0.0, 0.25, 0.5, 0.75]`
- `N`: `12`
- `seed`: `1959`
- `all_checks_pass`: `True`

| a | r_plus | r_ergosphere_eq | r_min_observed | inside_ergo | outside_ergo | prograde | retrograde | undecided |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.00 | 2 | 2 | 2.55819746857 | 0 | 12 | 0 | 0 | 5 |
| 0.25 | 1.96824583655 | 2 | 2.53198444031 | 0 | 12 | 34 | 32 | 66 |
| 0.50 | 1.86602540378 | 2 | 2.44883542862 | 0 | 12 | 34 | 32 | 66 |
| 0.75 | 1.66143782777 | 2 | 2.28890021611 | 0 | 12 | 34 | 32 | 66 |

Interpretation: these are scaffold diagnostics only. Undecided pairs are unknown, not non-relations.
