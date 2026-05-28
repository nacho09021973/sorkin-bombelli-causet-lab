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

The K5 prograde/retrograde local cone diagnostic is
`audit_kerr_k5_prograde_retrograde_local_cone_001.py`.
It is a **local null-slope diagnostic**, not a Kerr causal solver.

It measures the local prograde/retrograde asymmetry of the equatorial Kerr
light cone caused by the frame-dragging term g_tphi.  At fixed r with
dr = dtheta = 0, imposing ds² = 0 gives:

    omega_± = (-g_tphi ± sqrt(g_tphi² - g_tt g_phiphi)) / g_phiphi
    omega_center = (omega_+ + omega_-)/2 = -g_tphi / g_phiphi
    omega_width  = (omega_+ - omega_-)/2 = sqrt(disc) / g_phiphi

K5 sweeps `a = 0.0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 0.25, 0.5, 0.75`
at the same fixed radial grid `r = [2.5, 3.0, 4.0, 6.0, 10.0]` as K4 and checks:

1. **Discriminant positivity** (all a): `disc > 0` outside the horizon.
2. **Cone width positivity** (all a): `omega_width > 0`.
3. **Exact center identity** (all a): `(omega_+ + omega_-)/2 = -g_tphi/g_phiphi`
   (algebraic identity; residual ≤ 1e-12).
4. **Schwarzschild symmetry** (a=0): `omega_+ = -omega_-`, `omega_center = 0`.
5. **Frame-dragging sign** (a>0): `omega_center > 0`
   (g_tphi < 0 tilts the cone in the prograde direction).
6. **Linear scaling** (0<a≤0.01): `omega_center/a → 2M/r³`
   (leading-order frame-dragging; residual ≤ 1e-4).

K5 does NOT:

- Implement Kerr causal inference of any kind.
- Integrate null geodesics.
- Claim global causal reachability.
- Assert global causal relations.  Local null slopes are not global causal
  relations.

K5 causal accounting follows K1–K4 invariants:

- `a=0`: Schwarzschild/Kerr scaffold control counts preserved.
- `a>0`: all unordered pairs remain undecided
  (true=0, false=0, undecided=N*(N-1)/2).

The K5 artifact is
`kerr_k5_prograde_retrograde_local_cone_001_n12_seed1959.{csv,json,md,png}`.

The PNG shows a 2×2 diagnostic panel: mean cone center vs spin, small-a
linear scaling residuals, omega_± vs r for representative spins, and
discriminant/width positivity vs spin.

The K6 near-horizon convergence audit is `audit_kerr_k6_zamo_omega_horizon_001.py`.
It is a **near-horizon convergence diagnostic**, not a Kerr causal solver.

It verifies that the local ZAMO angular velocity
`omega_ZAMO = -g_tphi/g_phiphi` (same as omega_center from K5) converges to
the known-truth horizon angular velocity `Omega_H = a/(r_+^2 + a^2)` as
`r -> r_+` from outside the horizon.  The convergence is linear in
`delta = r - r_+`.

K6 fixes `M=1`, sweeps `a = 0.0, 0.25, 0.5, 0.75, 0.9`, and evaluates
`omega_ZAMO` at `r_+ + delta` for `delta = (1e-1, 3e-2, 1e-2, 3e-3, 1e-3)`.

K6 freezes these controls:

- `a=0`: `Omega_H = omega_ZAMO = 0` everywhere; trivially True by convention.
- `a>0`: `omega_ZAMO < Omega_H` at all evaluation points (monotone from below).
- `a>0`: residuals `|omega_ZAMO - Omega_H|` decrease strictly as `delta -> 0`,
  consistent with `O(delta)` linear convergence.
- Causal accounting: `a>0` => all global pairs undecided
  (true=0, false=0, undecided=N*(N-1)/2).

K6 does NOT:

- Cross the horizon.
- Implement Kerr causal inference of any kind.
- Integrate null geodesics.
- Claim global causal reachability.
- Convergence of `omega_ZAMO` to `Omega_H` is a **local metric identity**,
  not a causal relation.  It satisfies the level-A criterion from the Hawking
  consistency guardrail (AGENTS.md): a closed-form identity check, not a
  discrete pipeline rediscovery.

Connection to the K-sequence: K5 measured `omega_ZAMO` asymmetry at fixed `r`;
K6 measures `omega_ZAMO` convergence to `Omega_H` near the horizon.  Together
they bridge local frame-dragging -> horizon angular velocity without crossing
the Hawking/Bekenstein thermodynamic guardrail.

The K6 artifact is
`kerr_k6_zamo_omega_horizon_001_n12_seed1959.{csv,json,md,png}`.

The PNG shows a 2×2 diagnostic panel: `omega_ZAMO` vs delta per spin (with
`Omega_H` dashed), residual `|omega_ZAMO - Omega_H|` vs delta on a log-log
scale with a slope-1 reference, `Omega_H` vs `a` (analytic curve + test
points), and residual/delta vs delta (linear convergence rate plateau).

The K7 equatorial null-potential audit is
`audit_kerr_k7_equatorial_null_potential_001.py`.
It is an **analytic radial-potential known-truth audit**, not a Kerr causal
solver.

It verifies the Kerr equatorial null-geodesic radial potential

```text
R(r; a, b) = [r^2 + a^2 - a*b]^2 - Delta*(b - a)^2
Delta = r^2 - 2*M*r + a^2
```

and its analytic derivative `dR/dr` at the circular photon orbit radii.

K7 fixes `M=1`, sweeps `a = 0.0, 1e-4, 1e-3, 1e-2, 0.1, 0.25, 0.5, 0.75, 0.9`,
and computes for each spin:
- Prograde circular photon orbit radius:
  `r_ph_pro = 2M[1 + cos((2/3)*arccos(-a/M))]`
- Retrograde circular photon orbit radius:
  `r_ph_retro = 2M[1 + cos((2/3)*arccos(+a/M))]`
- Impact parameters:
  `b_ph_pro = a + 2*r_ph_pro*sqrt(Delta(r_ph_pro))/(r_ph_pro - M)` (prograde, b>0)
  `b_ph_retro = a - 2*r_ph_retro*sqrt(Delta(r_ph_retro))/(r_ph_retro - M)` (retrograde, b<0)
- `R(r_ph_pro; b_ph_pro)`, `dR/dr` at `r_ph_pro`
- `R(r_ph_retro; b_ph_retro)`, `dR/dr` at `r_ph_retro`

K7 freezes these controls:

- `a=0`: Schwarzschild photon sphere at `r=3M`, `b_pro=+3√3M`, `b_retro=-3√3M`.
- All spins: `|R| <= 1e-9` and `|dR/dr| <= 1e-9` at circular orbit radii.
- `a>0`: prograde orbit inside `3M`, retrograde orbit outside `3M`.
- `a>0`: `b_pro > 0`, `b_retro < 0`.
- Causal accounting: `a>0` => all global pairs undecided.

K7 does NOT:

- Integrate null geodesics.
- Decide causal reachability between sprinkled events.
- Create Kerr causal relations between any events.
- Cross the Hawking/Bekenstein thermodynamic guardrail.

It is a **preflight check** before any future Kerr geodesic integration.
It satisfies the level-A criterion from the Hawking consistency guardrail
(AGENTS.md): a closed-form identity check, not a discrete pipeline rediscovery.

Connection to the K-sequence:
- K5 measured local null slopes `dphi/dt` at fixed `r`.
- K6 measured `omega_ZAMO` convergence to `Omega_H` near the horizon.
- K7 introduces the radial potential `R(r; a, b)`, verifying the circular
  photon orbit structure.

The K7 artifact is
`kerr_k7_equatorial_null_potential_001_n12_seed1959.{csv,json,md,png}`.

The PNG shows a 2×2 diagnostic panel: circular photon orbit radii vs spin,
impact parameters vs spin, residuals `|R|` and `|dR/dr|` on a log-log scale
with the tolerance line, and photon orbit clearance above the horizon vs spin.

---

## S4 Kerr K8: Equatorial Null Radial-Flow Preflight

The K8 equatorial null radial-flow preflight is
`audit_kerr_k8_equatorial_null_radial_flow_001.py`.
It is a **first numerical ODE preflight audit**, not a Kerr causal solver.

It integrates the equatorial Kerr null radial-flow equation

```text
dr/dlambda = s * sqrt(R(r; a, b)) / r^2     [Sigma = r^2 for theta=pi/2]
s = +1 outgoing, -1 ingoing
```

using a local RK4 integrator (no new dependencies) and verifies
known-truth trajectory properties.

K8 fixes `M=1`, sweeps `a = 0.0, 0.25, 0.5, 0.75`, and for each spin
runs three cases:

- `outgoing_b0`: `r0=5M`, `b=0`, `s=+1`
- `ingoing_b0`: `r0=10M`, `b=0`, `s=-1`
- `circular_pro`: `r0=r_ph_pro`, `b=b_ph_pro`, `s=+1`

The safe choice `b=0` gives `R(r; a, 0) = r²(r²+a²) + 2Ma²r ≥ 0` for all
`r > 0`; no turning point for any `r > 0`.

K8 freezes these controls:

- All trajectory points remain exterior to the outer horizon.
- `R(r) >= 0` along the trajectory (no forbidden-region excursion).
- RHS consistency: max central-difference error ≤ 0.01 for all cases.
- `b=0` outgoing/ingoing trajectories are monotonically increasing/decreasing.
- Schwarzschild limit (`a=0`, `b=0`): `dr/dlambda = ±1` everywhere (constant);
  `r(lambda) = r0 ± lambda`; RK4 error ≤ 1e-10 (machine precision).
- Circular orbit drift: `|r_final - r0| < 1e-6` (advisory; prograde circular
  orbit has `R ≈ 0`, so the RHS ≈ 0 and the trajectory barely moves).
- Causal accounting: `a>0` => all global pairs undecided (K1–K8 invariant).

K8 does NOT:

- Decide causal reachability between sprinkled events.
- Create Kerr causal relations between any pair.
- Constitute a global Kerr causal solver of any kind.
- Cross the Hawking/Bekenstein thermodynamic guardrail.

Connection to the K-sequence:
- K7 verified `R(r_ph; a, b_ph) = 0` and `dR/dr = 0` at circular photon orbit radii.
- K8 takes the first numerical step: integrating `dr/dlambda = s*sqrt(R)/r^2`
  with `b=0` (safe) trajectories and verifying the Schwarzschild limit exactly.

The K8 artifact is
`kerr_k8_equatorial_null_radial_flow_001_n12_seed1959.{csv,json,md,png}`.

The PNG shows a 2×2 diagnostic panel: `b=0` outgoing trajectories `r(λ)` for
all spins + analytic line (`a=0`); `b=0` ingoing trajectories; RHS consistency
error vs step index (semilog y); circular orbit radial drift vs spin (log y).

---

## K1–K8 known-truth status (as of K8)

What the K-sequence has verified by known-truth checks:

- `a -> 0` Schwarzschild limit: Boyer-Lindquist equatorial metric reduces to
  Schwarzschild exact to floating-point tolerance (K3, K4).
- Frame-dragging sign: `g_tphi = -2Ma/r < 0` for all `a > 0` (K3).
- Perturbative metric scaling: `g_tphi/a`, `g_phiphi - r²`, `g_rr = r²/Δ`,
  and horizon shift `r_+ - 2M = O(a²)` all confirmed to within stated
  tolerances across `a = 0` to `a = 0.1` (K4).
- Local light cone structure: equatorial null slopes `omega_±` are real,
  well-separated, and symmetric at `a=0`; discriminant positive at all
  exterior points (K3, K5).
- Prograde/retrograde asymmetry: `omega_center = -g_tphi/g_phiphi > 0`
  for `a > 0`; linear scaling `omega_center/a -> 2M/r³` for small `a` (K5).
- Horizon angular velocity recovered from outside: `omega_ZAMO(r_+ + δ) -> Ω_H`
  with `O(δ)` linear convergence for `a = 0.25, 0.5, 0.75, 0.9` (K6).
- Circular photon orbit structure: `R(r_ph; a, b_ph) = 0` and `dR/dr = 0`
  at circular orbit radii for all spins to within 1e-9; prograde orbit inside
  `3M`, retrograde outside `3M` for `a > 0`; exact Schwarzschild limit
  `r_ph=3M`, `b=±3√3M` at `a=0` (K7).
- Equatorial null radial-flow integration: `dr/dlambda = s*sqrt(R)/r^2`
  with `b=0` integrates correctly via RK4; Schwarzschild limit error ≤ 1e-10
  (machine precision); circular orbit drift < 1e-6; RHS consistency ≤ 2e-6
  for all tested spins (K8).

**Addendum after K7–K8:**

- K7 verifies the equatorial Kerr null radial potential `R(r; a, b)` and the
  circular photon orbit conditions `R = 0` and `dR/dr = 0`, recovering the
  Schwarzschild photon sphere (`r=3M`, `b=±3√3M`) exactly at `a=0`.
- K8 integrates only the radial equation `dr/dλ = ±sqrt(R)/r²` for safe
  control cases (`b=0`, circular-orbit starts) and verifies radial RHS
  consistency, exterior support, circular-orbit drift control, and the
  Schwarzschild radial-flow limit `r(λ) = r0 ± λ` to machine precision.

What has NOT been shown:

- No full Kerr geodesic integration: `t(λ)` and `φ(λ)` are not integrated;
  K8 advances only the `r`-equation.
- No point-to-point shooting between arbitrary events.
- No causal reachability between sprinkled events; all global pairs remain
  undecided for `a != 0` (K1–K8 invariant).
- `b=0` in K8 is a safe radial-flow control case, not a generic Kerr
  null-geodesic family; generic `b` trajectories and angular evolution
  `dphi/dlambda` are not yet implemented.
- No global causal reachability has been claimed.
- No Hawking/Bekenstein thermodynamic quantity has been reconstructed from
  the discrete pipeline output (level B of the Hawking guardrail, AGENTS.md).
- The circular photon orbit identity and the near-horizon omega_ZAMO convergence
  are closed-form metric identity checks (level A), not discrete pipeline
  rediscoveries.

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
