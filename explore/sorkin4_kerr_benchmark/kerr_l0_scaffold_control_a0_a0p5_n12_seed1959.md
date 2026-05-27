# S4-K0 Kerr L0 scaffold control audit

Scope: this is a bookkeeping scaffold, not a Kerr causal-relation solver.

Checks frozen by this artifact:

- `a=0`: exact control against the Schwarzschild exterior subset for N=12, seed=1959.
- `a=0.5`: scaffold mode leaves all 66 unordered pairs undecided.
- Undecided pairs are unknown pairs, not decided non-relations.

Result:

- `a0_checks_pass`: `True`
- `spin_checks_pass`: `True`
- `all_checks_pass`: `True`

Counts:

- `a=0`: true `1`, false `64`, undecided `1`
- `a=0.5`: true `0`, false `0`, undecided `66`

Next admissible Kerr step: add only diagnostics with an explicit `a=0` regression gate before making any physical claim.
