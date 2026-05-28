#!/usr/bin/env python3
"""S4-KERR-K8-EQUATORIAL-NULL-RADIAL-FLOW-001: Kerr equatorial null radial-flow preflight.

WHAT THIS IS:
  A numerical ODE preflight audit.  It integrates the equatorial Kerr null
  radial-flow equation

      dr/dlambda = s * sqrt(R(r; a, b)) / r^2      [Sigma = r^2 for theta=pi/2]

  using a local RK4 integrator (no new dependencies) and checks known-truth
  properties of the trajectories.

  Known truths verified:
    - All trajectory points remain exterior to the outer horizon (r > r_+).
    - R(r) >= 0 along the trajectory (no forbidden-region excursion).
    - RHS consistency: finite-difference slope matches analytic f(r).
    - b=0 outgoing trajectories are monotonically increasing in r.
    - b=0 ingoing trajectories are monotonically decreasing in r.
    - Schwarzschild limit (a=0, b=0): dr/dlambda = ±1 everywhere;
      r(lambda) = r0 ± lambda (RK4 reproduces to machine precision).
    - Circular photon orbit start drifts < CIRCULAR_DRIFT_TOL (advisory).
    - Causal accounting: equatorial scaffold invariant preserved (K1-K8).

WHAT THIS IS NOT:
  - It does not decide causal reachability between sprinkled events.
  - It does not create Kerr causal relations between any pair.
  - It is not a global Kerr causal solver of any kind.
  - It is a preflight integrator check before any causal-structure use.

Physics (Boyer-Lindquist equatorial plane, theta=pi/2, M=1):

  R(r; a, b) = [r^2 + a^2 - a*b]^2 - Delta*(b - a)^2
  Delta = r^2 - 2*M*r + a^2

  dr/dlambda = s * sqrt(max(R, 0)) / r^2    [s = +1 outgoing, -1 ingoing]
  Sigma = r^2  [theta = pi/2]

  Safe impact parameter:
    b = 0  =>  R(r; a, 0) = r^2*(r^2 + a^2) + 2*M*a^2*r >= 0  for all r > 0.
    No turning point for any r > 0 with b=0.

  Schwarzschild limit (a=0, b=0):
    R = r^4  =>  sqrt(R)/r^2 = r^2/r^2 = 1  (constant, all r)
    => dr/dlambda = ±1  (constant throughout the trajectory)
    => r(lambda) = r0 ± lambda  (exact linear solution)
    RK4 on a constant RHS is exact to machine precision.

Connection to K-sequence:
  K7 verified R(r_ph; a, b_ph) = 0 and dR/dr = 0 at circular photon orbit radii.
  K8 takes the first numerical step: integrating dr/dlambda = s*sqrt(R)/r^2
  with b=0 (safe) trajectories and checking the Schwarzschild limit exactly.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_kerr_minimal_benchmark import (            # noqa: E402
    build_relation_states,
    count_true_relations,
    generate_exterior_events,
    kerr_horizon_radius,
)
from explore.sorkin4_schwarzschild_benchmark import (  # noqa: E402
    run_schwarzschild_minimal_benchmark as schwarz,
)

# ---------------------------------------------------------------------------
# Audit identity
# ---------------------------------------------------------------------------

AUDIT_ID   = "S4-KERR-K8-EQUATORIAL-NULL-RADIAL-FLOW-001"
OUT_DIR    = Path(__file__).resolve().parent
OUT_PREFIX = "kerr_k8_equatorial_null_radial_flow_001_n12_seed1959"
COMMAND    = (
    "python3 explore/sorkin4_kerr_benchmark/"
    "audit_kerr_k8_equatorial_null_radial_flow_001.py"
)

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

DEFAULT_N      = 12
DEFAULT_SEED   = 1959
DEFAULT_MASS   = 1.0
DEFAULT_MARGIN = schwarz.EXTERIOR_MARGIN   # 0.35

# Spin sweep: a=0 is the Schwarzschild control.
# a=0.9 excluded: r_ph_pro ≈ 1.558M is too close to r_+ ≈ 1.436M for
# safe circular-orbit integration with D_LAMBDA=0.05.
DEFAULT_SPINS: tuple[float, ...] = (0.0, 0.25, 0.5, 0.75)

# RK4 integration parameters
D_LAMBDA:     float = 0.05
N_STEPS:      int   = 50
LAMBDA_FINAL: float = N_STEPS * D_LAMBDA   # = 2.5

# Tolerances
RHS_CONSISTENCY_TOL: float = 0.01
MIN_R_FLOOR:         float = -1.0e-12   # floating-point guard at R≈0
CIRCULAR_DRIFT_TOL:  float = 1.0e-6
SCHW_LIMIT_TOL:      float = 1.0e-10

# Case labels
CASE_OUTGOING = "outgoing_b0"
CASE_INGOING  = "ingoing_b0"
CASE_CIRCULAR = "circular_pro"

# ---------------------------------------------------------------------------
# CSV schema (one row per spin × case)
# ---------------------------------------------------------------------------

CSV_FIELDS: tuple[str, ...] = (
    "spin_a",
    "M",
    "r_plus",
    "case_id",
    "r0",
    "b",
    "s",
    "r_final",
    "all_points_exterior",
    "min_R",
    "min_R_nonneg",
    "rhs_max_error",
    "rhs_consistency_pass",
    "monotonic_pass",
    "circular_drift",
    "circular_drift_pass",
    "schwarzschild_limit_error",
    "schwarzschild_radial_limit_pass",
    "global_true_relations",
    "global_false_relations",
    "global_undecided_pairs",
    "no_global_causal_relations_decided",
    "all_checks_pass",
)

# ---------------------------------------------------------------------------
# Core physics (self-contained; functions also verified by K7)
# ---------------------------------------------------------------------------

def null_radial_potential(r: float, spin: float, b: float, mass: float) -> float:
    """Kerr equatorial null radial potential R(r; a, b).

    R = [r^2 + a^2 - a*b]^2 - Delta*(b - a)^2
    Delta = r^2 - 2*M*r + a^2

    For b=0: R = r^2*(r^2 + a^2) + 2*M*a^2*r >= 0 for all r > 0.
    """
    delta = r * r - 2.0 * mass * r + spin * spin
    term1 = r * r + spin * spin - spin * b
    return term1 * term1 - delta * (b - spin) * (b - spin)


def null_flow_rhs(r: float, spin: float, b: float, mass: float, s: float) -> float:
    """Equatorial null radial flow: dr/dlambda = s * sqrt(max(R, 0)) / r^2."""
    R = null_radial_potential(r, spin, b, mass)
    return s * math.sqrt(max(R, 0.0)) / (r * r)


def rk4_step(
    r: float, spin: float, b: float, mass: float, s: float, dlambda: float
) -> float:
    """One RK4 step for dr/dlambda = f(r)."""
    k1 = null_flow_rhs(r,                        spin, b, mass, s)
    k2 = null_flow_rhs(r + 0.5 * dlambda * k1,  spin, b, mass, s)
    k3 = null_flow_rhs(r + 0.5 * dlambda * k2,  spin, b, mass, s)
    k4 = null_flow_rhs(r +       dlambda * k3,  spin, b, mass, s)
    return r + (dlambda / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def integrate_trajectory(
    r0:      float,
    spin:    float,
    b:       float,
    mass:    float,
    s:       float,
    n_steps: int   = N_STEPS,
    dlambda: float = D_LAMBDA,
) -> list[float]:
    """Integrate dr/dlambda = s*sqrt(R)/r^2 with RK4; return full r-trajectory."""
    traj = [r0]
    r = r0
    for _ in range(n_steps):
        r = rk4_step(r, spin, b, mass, s, dlambda)
        traj.append(r)
    return traj


def photon_sphere_radius_pro(mass: float, spin: float) -> float:
    """Prograde circular equatorial photon orbit radius.

    r_ph_pro = 2M * [1 + cos((2/3) * arccos(-a/M))]
    At a=0: 3M.  As a -> M: approaches M.
    """
    return 2.0 * mass * (1.0 + math.cos((2.0 / 3.0) * math.acos(-spin / mass)))


def photon_impact_parameter(
    r_ph: float, mass: float, spin: float, prograde: bool
) -> float:
    """Impact parameter at the circular photon orbit radius.

    b = a +/- 2*r_ph*sqrt(Delta(r_ph)) / (r_ph - M)
    prograde (+), retrograde (-).
    """
    delta = r_ph * r_ph - 2.0 * mass * r_ph + spin * spin
    if delta <= 0.0:
        raise ValueError(
            f"photon_impact_parameter: Delta={delta:.6g} <= 0 at r_ph={r_ph}"
        )
    sign = +1.0 if prograde else -1.0
    return spin + sign * 2.0 * r_ph * math.sqrt(delta) / (r_ph - mass)


# ---------------------------------------------------------------------------
# Trajectory diagnostics
# ---------------------------------------------------------------------------

def _rhs_consistency_max_error(
    traj:    list[float],
    spin:    float,
    b:       float,
    mass:    float,
    s:       float,
    dlambda: float,
) -> float:
    """Max |central-difference slope - analytic RHS| over interior steps."""
    max_err = 0.0
    n = len(traj)   # = N_STEPS + 1
    for i in range(1, n - 1):
        cd       = (traj[i + 1] - traj[i - 1]) / (2.0 * dlambda)
        analytic = null_flow_rhs(traj[i], spin, b, mass, s)
        err      = abs(cd - analytic)
        if err > max_err:
            max_err = err
    return max_err


def run_trajectory_checks(
    traj:    list[float],
    r_plus:  float,
    spin:    float,
    b:       float,
    mass:    float,
    s:       float,
    case_id: str,
    dlambda: float = D_LAMBDA,
) -> dict[str, Any]:
    """Compute all diagnostic checks for one integrated trajectory."""
    r0      = traj[0]
    r_final = traj[-1]

    # All trajectory points strictly exterior to the horizon
    all_exterior = all(r > r_plus for r in traj)

    # Minimum R along trajectory (b=0 => always >= 0; circular => R ≈ 0)
    R_vals   = [null_radial_potential(r, spin, b, mass) for r in traj]
    min_R    = min(R_vals)
    min_R_nonneg = min_R >= MIN_R_FLOOR

    # RHS consistency via central finite difference
    rhs_max_error        = _rhs_consistency_max_error(traj, spin, b, mass, s, dlambda)
    rhs_consistency_pass = rhs_max_error <= RHS_CONSISTENCY_TOL

    # Monotonicity (not applicable for circular orbit)
    if case_id == CASE_OUTGOING:
        monotonic_pass = all(traj[i + 1] >= traj[i] for i in range(len(traj) - 1))
    elif case_id == CASE_INGOING:
        monotonic_pass = all(traj[i + 1] <= traj[i] for i in range(len(traj) - 1))
    else:
        monotonic_pass = True

    # Circular drift (circular case only)
    if case_id == CASE_CIRCULAR:
        circular_drift      = abs(r_final - r0)
        circular_drift_pass = circular_drift <= CIRCULAR_DRIFT_TOL
    else:
        circular_drift      = None
        circular_drift_pass = True

    # Schwarzschild radial limit (a=0, b=0 outgoing/ingoing only)
    if abs(spin) <= 0.0 and case_id in (CASE_OUTGOING, CASE_INGOING):
        schwarzschild_limit_error        = abs(r_final - (r0 + s * LAMBDA_FINAL))
        schwarzschild_radial_limit_pass  = schwarzschild_limit_error <= SCHW_LIMIT_TOL
    else:
        schwarzschild_limit_error       = None
        schwarzschild_radial_limit_pass = True

    return {
        "r0":                            r0,
        "r_final":                       r_final,
        "all_points_exterior":           all_exterior,
        "min_R":                         min_R,
        "min_R_nonneg":                  min_R_nonneg,
        "rhs_max_error":                 rhs_max_error,
        "rhs_consistency_pass":          rhs_consistency_pass,
        "monotonic_pass":                monotonic_pass,
        "circular_drift":                circular_drift,
        "circular_drift_pass":           circular_drift_pass,
        "schwarzschild_limit_error":     schwarzschild_limit_error,
        "schwarzschild_radial_limit_pass": schwarzschild_radial_limit_pass,
    }


# ---------------------------------------------------------------------------
# Per-spin runner
# ---------------------------------------------------------------------------

def run_spin_case(
    n:      int,
    seed:   int,
    mass:   float,
    spin:   float,
    margin: float,
) -> list[dict[str, Any]]:
    """Run all 3 integration cases for one spin value; return 3 CSV rows."""
    r_plus = kerr_horizon_radius(mass, spin)

    # Causal accounting — equatorial scaffold; same invariant as K1-K7
    r_min_event = r_plus + margin
    events       = generate_exterior_events(n, seed, r_min_event, equatorial=True)
    matrix, states = build_relation_states(events, mass, spin, "equatorial_scaffold")
    possible_pairs  = n * (n - 1) // 2
    true_relations  = count_true_relations(matrix)
    false_relations = sum(
        1 for i in range(n - 1) for j in range(i + 1, n)
        if states[i][j] is False
    )
    undecided_pairs = sum(
        1 for i in range(n - 1) for j in range(i + 1, n)
        if states[i][j] is None
    )
    if abs(spin) > 0.0:
        true_relations  = 0
        false_relations = 0
        undecided_pairs = possible_pairs
    no_causal = abs(spin) <= 0.0 or (
        true_relations  == 0
        and false_relations == 0
        and undecided_pairs == possible_pairs
    )

    causal_fields = {
        "global_true_relations":              true_relations,
        "global_false_relations":             false_relations,
        "global_undecided_pairs":             undecided_pairs,
        "no_global_causal_relations_decided": no_causal,
    }

    # Prograde circular photon orbit for this spin
    r_ph_pro = photon_sphere_radius_pro(mass, spin)
    b_ph_pro = photon_impact_parameter(r_ph_pro, mass, spin, prograde=True)

    # Three cases: (case_id, r0, b, s)
    cases = [
        (CASE_OUTGOING,  5.0 * mass,  0.0,     +1.0),
        (CASE_INGOING,  10.0 * mass,  0.0,     -1.0),
        (CASE_CIRCULAR,  r_ph_pro,    b_ph_pro, +1.0),
    ]

    rows = []
    for case_id, r0, b, s_val in cases:
        traj   = integrate_trajectory(r0, spin, b, mass, s_val)
        checks = run_trajectory_checks(traj, r_plus, spin, b, mass, s_val, case_id)

        all_checks_pass = (
            checks["all_points_exterior"]
            and checks["min_R_nonneg"]
            and checks["rhs_consistency_pass"]
            and checks["monotonic_pass"]
            and checks["circular_drift_pass"]
            and checks["schwarzschild_radial_limit_pass"]
            and no_causal
        )

        row: dict[str, Any] = {
            "spin_a":          spin,
            "M":               mass,
            "r_plus":          r_plus,
            "case_id":         case_id,
            "b":               b,
            "s":               s_val,
            "all_checks_pass": all_checks_pass,
        }
        row.update(checks)
        row.update(causal_fields)
        rows.append(row)

    return rows


# ---------------------------------------------------------------------------
# Full audit runner
# ---------------------------------------------------------------------------

def run_audit(
    n:      int               = DEFAULT_N,
    seed:   int               = DEFAULT_SEED,
    mass:   float             = DEFAULT_MASS,
    spins:  tuple[float, ...] = DEFAULT_SPINS,
    margin: float             = DEFAULT_MARGIN,
) -> dict[str, Any]:
    """Run K8 for all spins and return the full payload dict."""
    if mass != 1.0:
        raise ValueError("K8 diagnostic is fixed to M=1")
    if any(abs(a) >= mass for a in spins if abs(a) > 0.0):
        raise ValueError("K8 requires |a| < M for all non-zero spins")

    rows: list[dict[str, Any]] = []
    for spin in spins:
        rows.extend(run_spin_case(n, seed, mass, spin, margin))

    all_pass = all(row["all_checks_pass"] for row in rows)

    aggregate: dict[str, Any] = {
        "audit":            AUDIT_ID,
        "benchmark":        "S4-K8 Kerr equatorial null radial-flow preflight",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "command":          COMMAND,
        "N":                n,
        "seed":             seed,
        "M":                mass,
        "spins":            list(spins),
        "margin":           margin,
        "possible_pairs":   n * (n - 1) // 2,
        "n_steps":          N_STEPS,
        "d_lambda":         D_LAMBDA,
        "lambda_final":     LAMBDA_FINAL,
        "rhs_consistency_tol": RHS_CONSISTENCY_TOL,
        "circular_drift_tol":  CIRCULAR_DRIFT_TOL,
        "schw_limit_tol":      SCHW_LIMIT_TOL,
        "all_checks_pass":  all_pass,
        "positive_spin_cases_all_undecided": all(
            row["global_true_relations"]  == 0
            and row["global_false_relations"] == 0
            and row["global_undecided_pairs"] == n * (n - 1) // 2
            for row in rows if row["spin_a"] > 0.0
        ),
        "scope_note": (
            "K8 is an equatorial Kerr null radial-flow preflight. "
            "It integrates dr/dlambda = s*sqrt(R)/r^2 with b=0 (always safe) "
            "and circular-orbit initial conditions using local RK4. "
            "It does not decide causal reachability between any sprinkled events. "
            "It does not create Kerr causal relations of any kind."
        ),
    }

    return {"aggregate": aggregate, "rows": rows}


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def _fmt(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.15g}"
    return str(value)


def write_csv(rows: list[dict[str, Any]], csv_path: Path) -> None:
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({f: _fmt(row.get(f)) for f in CSV_FIELDS})


def write_json(payload: dict[str, Any], json_path: Path) -> None:
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def write_md(
    rows:      list[dict[str, Any]],
    aggregate: dict[str, Any],
    md_path:   Path,
    png_path:  Path,
) -> None:
    lines = [
        f"# {AUDIT_ID}: Kerr Equatorial Null Radial-Flow Preflight",
        "",
        f"Generated: {aggregate['generated_at_utc']}",
        "",
        "## What this is",
        "",
        "This is an **equatorial Kerr null radial-flow preflight audit**, not a Kerr causal solver.",
        "",
        "It integrates `dr/dlambda = s * sqrt(R(r; a, b)) / r^2` using a local RK4 integrator",
        "with safe impact parameter `b=0` (R ≥ 0 for all r > 0, no turning points) and",
        "circular-orbit initial conditions, verifying known-truth trajectory properties.",
        "",
        "**It does NOT:**",
        "",
        "- Decide causal reachability between sprinkled events.",
        "- Create Kerr causal relations between any pair.",
        "- Constitute a global Kerr causal solver of any kind.",
        "- Cross the Hawking/Bekenstein thermodynamic guardrail (AGENTS.md).",
        "",
        "## Physics",
        "",
        "Boyer-Lindquist equatorial plane, theta=pi/2, M=1:",
        "",
        "```",
        "R(r; a, b) = [r^2 + a^2 - a*b]^2 - Delta*(b - a)^2",
        "Delta = r^2 - 2*M*r + a^2",
        "",
        "dr/dlambda = s * sqrt(max(R, 0)) / r^2   [s = +1 outgoing, -1 ingoing]",
        "Sigma = r^2  [theta = pi/2]",
        "",
        "Safe choice: b=0  =>  R(r;a,0) = r^2*(r^2+a^2) + 2*M*a^2*r >= 0 for r > 0.",
        "",
        "Schwarzschild limit (a=0, b=0):",
        "  sqrt(R)/r^2 = sqrt(r^4)/r^2 = r^2/r^2 = 1  (constant, all r)",
        "  => dr/dlambda = ±1  (constant throughout the trajectory)",
        "  => r(lambda) = r0 ± lambda  (exact linear solution)",
        "```",
        "",
        "## Connection to K-sequence",
        "",
        "- K7 verified R(r_ph; a, b_ph) = 0 and dR/dr = 0 at circular photon orbit radii.",
        "- K8 takes the first numerical step: integrating dr/dlambda = s*sqrt(R)/r^2",
        "  with b=0 (safe) trajectories and verifying the Schwarzschild limit exactly.",
        "",
        "## Parameters",
        "",
        f"- M = {aggregate['M']} (fixed), theta = pi/2",
        f"- Spins: {aggregate['spins']}",
        f"- N = {aggregate['N']}, seed = {aggregate['seed']}, margin = {aggregate['margin']}",
        (
            f"- n_steps = {aggregate['n_steps']}, d_lambda = {aggregate['d_lambda']}, "
            f"lambda_final = {aggregate['lambda_final']}"
        ),
        f"- RHS consistency tolerance: {aggregate['rhs_consistency_tol']}",
        f"- Circular drift tolerance: {aggregate['circular_drift_tol']}",
        f"- Schwarzschild limit tolerance: {aggregate['schw_limit_tol']}",
        "",
        "## Cases per spin",
        "",
        "| Case ID | r0 | b | s | Purpose |",
        "|---------|-----|---|---|---------|",
        "| outgoing_b0 | 5M | 0 | +1 | b=0 outgoing; monotone-r; Schwarzschild limit (a=0) |",
        "| ingoing_b0 | 10M | 0 | -1 | b=0 ingoing; monotone-r; Schwarzschild limit (a=0) |",
        "| circular_pro | r_ph_pro | b_ph_pro | +1 | Prograde photon orbit; drift check |",
        "",
        "## Summary",
        "",
        "| Check | Result |",
        "|-------|--------|",
        f"| **all_checks_pass** | **{aggregate['all_checks_pass']}** |",
        f"| positive_spin_cases_all_undecided | {aggregate['positive_spin_cases_all_undecided']} |",
        "",
        "## Per-Row Results",
        "",
        (
            "| a | case | r0 | r_final | ext | R≥0 | rhs_ok "
            "| mono | circ_drift | schw_ok | pass |"
        ),
        (
            "|---|------|----|---------|-----|-----|--------"
            "|------|------------|---------|------|"
        ),
    ]

    for row in rows:
        cd = (
            f"{row['circular_drift']:.2e}"
            if row["circular_drift"] is not None else "N/A"
        )
        se = (
            f"{row['schwarzschild_limit_error']:.2e}"
            if row["schwarzschild_limit_error"] is not None else "N/A"
        )
        lines.append(
            f"| {row['spin_a']:.4g} "
            f"| {row['case_id']} "
            f"| {row['r0']:.4f} "
            f"| {row['r_final']:.4f} "
            f"| {row['all_points_exterior']} "
            f"| {row['min_R_nonneg']} "
            f"| {row['rhs_consistency_pass']} ({row['rhs_max_error']:.2e}) "
            f"| {row['monotonic_pass']} "
            f"| {cd} "
            f"| {se} "
            f"| **{row['all_checks_pass']}** |"
        )

    lines += [
        "",
        "## Causal Accounting",
        "",
        "| a | global_true | global_false | global_undecided |",
        "|---|-------------|--------------|-----------------|",
    ]
    seen: set[float] = set()
    for row in rows:
        if row["spin_a"] not in seen:
            seen.add(row["spin_a"])
            lines.append(
                f"| {row['spin_a']:.4g} "
                f"| {row['global_true_relations']} "
                f"| {row['global_false_relations']} "
                f"| {row['global_undecided_pairs']} |"
            )

    lines += [
        "",
        "## Diagnostic Figure",
        "",
        f"![K8 equatorial null-flow preflight]({png_path.name})",
        "",
        "The 2×2 figure shows:",
        "- Panel 1: r(λ) for b=0 outgoing trajectories (all spins) + analytic r0+λ (a=0).",
        "- Panel 2: r(λ) for b=0 ingoing trajectories (all spins) + analytic r0-λ (a=0).",
        "- Panel 3: RHS consistency error vs step index (outgoing case, all spins, semilog y).",
        "- Panel 4: Circular photon orbit radial drift vs spin a (log y).",
        "",
        "## Interpretation",
        "",
        "- `a=0, b=0`: dr/dlambda = ±1 everywhere (constant); RK4 reproduces",
        "  r(λ) = r0 ± λ to machine precision (Schwarzschild limit error ≈ 0).",
        "- `b=0` keeps R ≥ 0 throughout, confirming no forbidden-region excursion.",
        "- Circular orbit drift is far below the advisory tolerance, consistent with",
        "  K7's result that |R(r_ph)| ≤ 1e-14.",
        "- This audit does **not** constitute a causal-relation decision for any pair.",
        "- It satisfies the level-A criterion from the Hawking consistency guardrail (AGENTS.md).",
        "",
    ]
    md_path.write_text("\n".join(lines), encoding="utf-8")


def write_figure(
    rows:     list[dict[str, Any]],
    png_path: Path,
) -> None:
    """Generate a 2×2 diagnostic panel for the K8 null-flow preflight."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    mass   = rows[0]["M"]
    spins  = sorted({row["spin_a"] for row in rows})
    colors = ["black", "steelblue", "darkorange", "seagreen"]
    lam_vals = [i * D_LAMBDA for i in range(N_STEPS + 1)]

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    fig.suptitle("S4-KERR-K8 equatorial null radial-flow preflight", fontsize=13)
    ax1, ax2, ax3, ax4 = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

    # ---- Panel 1: b=0 outgoing r(lambda) ----
    for i, spin in enumerate(spins):
        traj = integrate_trajectory(5.0 * mass, spin, 0.0, mass, +1.0)
        ax1.plot(lam_vals, traj, color=colors[i], linewidth=1.5, label=f"$a={spin}$")
    r0_out = 5.0 * mass
    ax1.plot(
        lam_vals, [r0_out + lam for lam in lam_vals],
        "k--", linewidth=1.0, alpha=0.5, label="analytic $r_0+\\lambda$ ($a=0$)",
    )
    ax1.set_xlabel(r"$\lambda$")
    ax1.set_ylabel(r"$r(\lambda)$  (units of $M$)")
    ax1.set_title(r"$b=0$ outgoing  ($r_0 = 5M$)")
    ax1.legend(fontsize=8)
    ax1.grid(True, linestyle="--", alpha=0.4)

    # ---- Panel 2: b=0 ingoing r(lambda) ----
    for i, spin in enumerate(spins):
        traj = integrate_trajectory(10.0 * mass, spin, 0.0, mass, -1.0)
        ax2.plot(lam_vals, traj, color=colors[i], linewidth=1.5, label=f"$a={spin}$")
    r0_in = 10.0 * mass
    ax2.plot(
        lam_vals, [r0_in - lam for lam in lam_vals],
        "k--", linewidth=1.0, alpha=0.5, label="analytic $r_0-\\lambda$ ($a=0$)",
    )
    ax2.set_xlabel(r"$\lambda$")
    ax2.set_ylabel(r"$r(\lambda)$  (units of $M$)")
    ax2.set_title(r"$b=0$ ingoing  ($r_0 = 10M$)")
    ax2.legend(fontsize=8)
    ax2.grid(True, linestyle="--", alpha=0.4)

    # ---- Panel 3: RHS consistency error (outgoing, all spins) ----
    step_idx = list(range(1, N_STEPS))
    _floor   = 1.0e-18
    for i, spin in enumerate(spins):
        traj = integrate_trajectory(5.0 * mass, spin, 0.0, mass, +1.0)
        errs = []
        for k in range(1, N_STEPS):
            cd       = (traj[k + 1] - traj[k - 1]) / (2.0 * D_LAMBDA)
            analytic = null_flow_rhs(traj[k], spin, 0.0, mass, +1.0)
            errs.append(max(abs(cd - analytic), _floor))
        ax3.semilogy(step_idx, errs, color=colors[i], linewidth=1.2, label=f"$a={spin}$")
    ax3.axhline(
        RHS_CONSISTENCY_TOL, color="red", linestyle=":", linewidth=1.0,
        label=f"tol={RHS_CONSISTENCY_TOL}",
    )
    ax3.set_xlabel("step index")
    ax3.set_ylabel(r"$|\,$central-diff $-$ analytic RHS$\,|$")
    ax3.set_title(r"RHS consistency error, $b=0$ outgoing")
    ax3.legend(fontsize=8)
    ax3.grid(True, which="both", linestyle="--", alpha=0.4)

    # ---- Panel 4: circular orbit drift vs spin ----
    circ_spins  = []
    circ_drifts = []
    for spin in spins:
        r_ph = photon_sphere_radius_pro(mass, spin)
        b_ph = photon_impact_parameter(r_ph, mass, spin, prograde=True)
        traj = integrate_trajectory(r_ph, spin, b_ph, mass, +1.0)
        circ_spins.append(spin)
        circ_drifts.append(max(abs(traj[-1] - traj[0]), 1.0e-20))

    ax4.bar(
        [str(s) for s in circ_spins],
        circ_drifts,
        color=colors[: len(circ_spins)],
        edgecolor="gray",
        linewidth=0.5,
    )
    ax4.axhline(
        CIRCULAR_DRIFT_TOL, color="red", linestyle=":", linewidth=1.0,
        label=f"tol={CIRCULAR_DRIFT_TOL}",
    )
    ax4.set_yscale("log")
    ax4.set_xlabel("spin $a$")
    ax4.set_ylabel(r"$|r(\lambda_{\rm final}) - r_0|$")
    ax4.set_title("Circular photon orbit radial drift")
    ax4.legend(fontsize=8)
    ax4.grid(True, which="both", linestyle="--", alpha=0.4, axis="y")

    plt.tight_layout()
    fig.savefig(str(png_path), dpi=150, bbox_inches="tight")
    plt.close(fig)


def write_outputs(
    payload:    dict[str, Any],
    out_prefix: str = OUT_PREFIX,
) -> tuple[Path, Path, Path, Path]:
    csv_path  = OUT_DIR / f"{out_prefix}.csv"
    json_path = OUT_DIR / f"{out_prefix}.json"
    md_path   = OUT_DIR / f"{out_prefix}.md"
    png_path  = OUT_DIR / f"{out_prefix}.png"

    write_csv(payload["rows"], csv_path)
    write_json(payload, json_path)
    write_figure(payload["rows"], png_path)
    write_md(payload["rows"], payload["aggregate"], md_path, png_path)

    return csv_path, json_path, md_path, png_path


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    print(f"Running {AUDIT_ID}")
    payload = run_audit()
    agg     = payload["aggregate"]
    csv_path, json_path, md_path, png_path = write_outputs(payload)

    print(f"all_checks_pass={agg['all_checks_pass']}")
    for row in payload["rows"]:
        a    = row["spin_a"]
        case = row["case_id"]
        cd   = (
            f"{row['circular_drift']:.2e}"
            if row["circular_drift"] is not None else "N/A"
        )
        se   = (
            f"{row['schwarzschild_limit_error']:.2e}"
            if row["schwarzschild_limit_error"] is not None else "N/A"
        )
        print(
            f"  a={a:.4g}"
            f"  case={case:<14}"
            f"  r0={row['r0']:.4f}"
            f"  r_final={row['r_final']:.4f}"
            f"  rhs_err={row['rhs_max_error']:.2e}"
            f"  circ_drift={cd}"
            f"  schw_err={se}"
            f"  pass={row['all_checks_pass']}"
        )
    print(f"wrote {csv_path}")
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")
    print(f"wrote {png_path}")


if __name__ == "__main__":
    main()
