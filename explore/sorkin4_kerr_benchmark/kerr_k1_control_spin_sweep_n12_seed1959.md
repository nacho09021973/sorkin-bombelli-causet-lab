# S4-K1 Kerr scaffold invariant control audit

Scope: conservative scaffold invariant/control audit only. This does not implement Kerr causal physics.

- `M`: `1.0`
- `spins`: `[0.0, 0.25, 0.5, 0.75]`
- `N`: `12`
- `seed`: `1959`
- `possible_pairs`: `66`
- `all_checks_pass`: `True`

Controls:

- `a=0.0` must exactly match the frozen Schwarzschild exterior benchmark.
- `a>0` uses only exterior points with `r > r_+ + margin`.
- `a>0` leaves all unordered pairs undecided: true `0`, false `0`, undecided `N*(N-1)/2`.

Per-spin counts:

| a | r_plus | r_min | true | false | undecided | checks |
|---:|---:|---:|---:|---:|---:|:---:|
| 0.00 | 2 | 2.35 | 1 | 64 | 1 | True |
| 0.25 | 1.96824583655 | 2.31824583655 | 0 | 0 | 66 | True |
| 0.50 | 1.86602540378 | 2.21602540378 | 0 | 0 | 66 | True |
| 0.75 | 1.66143782777 | 2.01143782777 | 0 | 0 | 66 | True |

Interpretation: undecided means not decided by this scaffold. It is not a Kerr non-relation.
