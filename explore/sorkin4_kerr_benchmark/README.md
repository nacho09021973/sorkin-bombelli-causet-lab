# SORKIN-4 Kerr Scaffold

This directory is currently a Kerr exterior scaffold, not a Kerr causal
reconstruction.

The active L0 check is `audit_kerr_l0_scaffold.py`.  It freezes only two
bookkeeping guarantees:

- `a=0` reproduces the existing Schwarzschild exterior benchmark subset with
  the same generated events and relation states.
- `a=0.5` in scaffold mode leaves every pair undecided.

The corresponding artifact is
`kerr_l0_scaffold_control_a0_a0p5_n12_seed1959.{csv,json,md}`.

The K1 control audit is `audit_kerr_k1_control.py`.  It is still only a
scaffold invariant audit, not Kerr causal validation.  It fixes `M=1`, sweeps
`a = 0.0, 0.25, 0.5, 0.75`, enforces `|a| < M`, computes
`r_+ = M + sqrt(M^2-a^2)`, and samples only exterior points with
`r > r_+ + margin`.

K1 freezes these controls:

- `a=0.0` exactly matches the frozen Schwarzschild exterior benchmark.
- `a>0` does not decide causal relations; all unordered pairs remain
  undecided with true `0`, false `0`, undecided `N*(N-1)/2`.
- `all_checks_pass` is explicit in the CSV/JSON/MD artifacts.

The K1 artifact is
`kerr_k1_control_spin_sweep_n12_seed1959.{csv,json,md}`.

The K2 equatorial diagnostic is `audit_kerr_k2_equatorial_diagnostic.py`.
It keeps the same conservative boundary: no Kerr causal relations are decided.
It fixes `theta = pi/2`, records the equatorial ergosphere radius `2M`, and
counts only scaffold diagnostics such as exterior/ergosphere membership and
signed prograde/retrograde azimuthal pairs.

K2 freezes these controls:

- `M=1`, `a = 0.0, 0.25, 0.5, 0.75`, and `r > r_+ + margin`.
- `a=0.0` is the Schwarzschild control subset on equatorial points.
- `a>0` leaves all unordered pairs undecided with true `0`, false `0`,
  undecided `N*(N-1)/2`.
- Prograde/retrograde counts are kinematic bookkeeping, not a frame-dragging
  claim and not a causal-relation criterion.

The K2 artifact is
`kerr_k2_equatorial_diagnostic_n12_seed1959.{csv,json,md}`.

Interpretation:

- The `a=0` branch is a regression/control gate.
- The `a!=0` branch is coordinate/event scaffolding only.
- Undecided pairs are unknown pairs, not decided non-relations.
- Local cone and equatorial files in this directory are diagnostics only; they
  are not global Kerr null-geodesic causal decisions.

Before any Kerr physical claim, the next implementation must pass through an
explicit `a=0` regression gate and state exactly which diagnostic is being
added.
