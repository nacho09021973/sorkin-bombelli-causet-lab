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

The K3 local cone diagnostic is `audit_kerr_k3_local_cone_001.py`.
It is still conservative: no Kerr causal relations are decided.
It fixes `theta = pi/2` and introduces the first genuine Kerr-geometry
computation: the Boyer-Lindquist equatorial metric coefficients.

K3 computes for each sub-extremal spin in `(0.0, 0.25, 0.5, 0.75)`:
- `Delta = r² - 2Mr + a²`, `g_tt`, `g_tphi`, `g_rr`, `g_phiphi` at each
  event point (Delta is computed explicitly before the metric call).
- For small-displacement equatorial pairs satisfying `|dr| ≤ 1.0` and
  `|dphi| ≤ 0.5 rad`, the local quadratic interval `ds²` at the midpoint
  radius: `ds² = g_tt dt² + 2 g_tphi dt dphi + g_rr dr² + g_phiphi dphi²`.
- Classification of each evaluated pair as `timelike_local_candidate`,
  `nullish_local_candidate`, or `spacelike_local_candidate`.

K3 freezes these controls:

- `a=0.0` Schwarzschild reduction: `g_tphi = 0`, `g_rr = 1/(1-2M/r)`,
  `g_phiphi = r²`, all exact to floating-point tolerance.
- `a>0` frame-dragging sign: `g_tphi = -2Ma/r < 0` at all event points.
- `a>0` global causal accounting: all unordered pairs remain undecided
  with true `0`, false `0`, undecided `N*(N-1)/2`.
- Local interval labels are **metric-sign diagnostics only**.  They do not
  imply null geodesic connectivity, do not integrate Kerr geodesics, and do
  not constitute a Kerr causal solver of any kind.

The K3 artifact is
`kerr_k3_local_cone_001_n12_seed1959.{csv,json,md}`.

The K4 Schwarzschild-limit audit is `audit_kerr_k4_schwarzschild_limit_001.py`.
It is a **known-truth perturbative metric audit**, not a Kerr causal solver.

It verifies the analytic `a -> 0` Schwarzschild limit of the Boyer-Lindquist
equatorial metric, at a **fixed radial grid** `r = [2.5, 3.0, 4.0, 6.0, 10.0]`
(invariant to spin, safely outside `r_+(a=0.1) + margin ≈ 2.345`).

K4 sweeps `a = 0.0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1` and checks:

1. **Schwarzschild limit at a=0**: `g_tphi = 0`, `g_rr = 1/(1-2M/r)`,
   `g_phiphi = r²`, all within `1e-12`.
2. **Frame-dragging linearity** (a>0): `|g_tphi/a - (-2M/r)| ≤ 1e-12`
   (no catastrophic cancellation; both sides are O(M/r)).
3. **Azimuthal quadratic** (a>0): `|g_phiphi - r² - a²(1+2M/r)| ≤ 1e-12`
   (absolute residual; the ratio form `(g_phiphi-r²)/a²` has
   catastrophic cancellation for small `a`).
4. **g_rr formula** (all a): `g_rr = r²/Δ` (identity check, near machine epsilon).
5. **Horizon quadratic shift** (0 < a ≤ 0.01): `|(2M-r_+)/a² - 1/(2M)| ≤ 1e-3`
   (genuine O(a²) perturbative check; error ≈ a²/8 ≤ 1.25e-5 at a=0.01).

K4 causal accounting follows K1/K2/K3 invariants:

- `a=0`: Schwarzschild/Kerr scaffold control counts preserved.
- `a>0`: all unordered pairs remain undecided
  (true=0, false=0, undecided=N*(N-1)/2).

K4 does NOT:

- Implement Kerr causal inference of any kind.
- Integrate null geodesics.
- Claim global causal reachability.
- Create causal true/false relations for `a != 0`.

The K4 artifact is
`kerr_k4_schwarzschild_limit_001_n12_seed1959.{csv,json,md,png}`.

The PNG diagnostic figure shows perturbative scaling residuals vs. spin `a`
on a log-scale x-axis (excluding `a=0`), in a 2×2 panel.

---

Interpretation:

- The `a=0` branch is a regression/control gate.
- The `a!=0` branch is coordinate/event scaffolding only.
- Undecided pairs are unknown pairs, not decided non-relations.
- Local cone and equatorial files in this directory are diagnostics only; they
  are not global Kerr null-geodesic causal decisions.

Before any Kerr physical claim, the next implementation must pass through an
explicit `a=0` regression gate and state exactly which diagnostic is being
added.
