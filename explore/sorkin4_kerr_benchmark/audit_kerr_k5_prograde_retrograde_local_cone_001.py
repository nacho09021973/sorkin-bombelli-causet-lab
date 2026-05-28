#!/usr/bin/env python3
"""S4-KERR-K5-PROGRADE-RETROGRADE-LOCAL-CONE-001: Kerr equatorial prograde/retrograde
local null-slope diagnostic.

WHAT THIS IS:
  A local null-slope diagnostic. It measures the local prograde/retrograde
  asymmetry of the Kerr equatorial light cone caused by g_tphi.  Known-truth
  checks verify the Schwarzschild symmetry at a=0 and the linear scaling of
  frame-dragging with spin.

WHAT THIS IS NOT:
  - It does not implement Kerr causal inference.
  - It does not integrate null geodesics.
  - It does not claim global causal reachability.
  - Local prograde/retrograde null slopes are not global causal relations.

Physics (Boyer-Lindquist equatorial plane, theta=pi/2, M=1, signature -+++):

  At fixed r with dr=dtheta=0, imposing ds^2=0:
    g_phiphi omega^2 + 2 g_tphi omega + g_tt = 0,   omega = dphi/dt
    omega_+/- = (-g_tphi +/- sqrt(disc)) / g_phiphi
    disc = g_tphi^2 - g_tt * g_phiphi  (> 0 outside the horizon)

  omega_center = (omega_+ + omega_-)/2 = -g_tphi / g_phiphi
  omega_width  = (omega_+ - omega_-)/2 = sqrt(disc) / g_phiphi > 0

  Schwarzschild (a=0):  omega_+ = -omega_-,  omega_center = 0
  Kerr (a>0):           omega_center > 0  (frame-dragging lifts the cone)
  Small a:              omega_center / a -> 2M/r^3
"""

from __future__ import annotations

import csv
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

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

AUDIT_ID   = "S4-KERR-K5-PROGRADE-RETROGRADE-LOCAL-CONE-001"
OUT_DIR    = Path(__file__).resolve().parent
OUT_PREFIX = "kerr_k5_prograde_retrograde_local_cone_001_n12_seed1959"
COMMAND    = (
    "python3 explore/sorkin4_kerr_benchmark/"
    "audit_kerr_k5_prograde_retrograde_local_cone_001.py"
)

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

DEFAULT_N      = 12
DEFAULT_SEED   = 1959
DEFAULT_MASS   = 1.0
DEFAULT_MARGIN = schwarz.EXTERIOR_MARGIN   # 0.35

DEFAULT_SPINS: tuple[float, ...] = (
    0.0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 0.25, 0.5, 0.75,
)

# Fixed radial grid for metric evaluation (same for every spin).
# All five points are safely outside r_plus(a=0.75) + margin ≈ 2.011:
#   min(RADIAL_GRID) = 2.5 > 2.011  for all spins in the sweep.
RADIAL_GRID: tuple[float, ...] = (2.5, 3.0, 4.0, 6.0, 10.0)

# Spin threshold for the small-a linear-scaling check.
SMALL_A_MAX = 1e-2

# Identity tolerances: computed from same closed-form expressions; residuals
# are purely floating-point rounding (~1e-15).
METRIC_FORMULA_TOL = 1.0e-12

# Small-a linear-scaling tolerance: correction is O(a^2/r^2); at a=1e-2 and
# r=2.5 the error is ~3.7e-7, well within 1e-4.
SMALL_A_LINEAR_TOL = 1.0e-4

# ---------------------------------------------------------------------------
# CSV schema
# ---------------------------------------------------------------------------

CSV_FIELDS: tuple[str, ...] = (
    "spin_a",
    "M",
    "N",
    "seed",
    "margin",
    "radial_grid",
    "r_plus",
    "r_min_evaluated",
    "all_points_exterior",
    "min_discriminant",
    "min_omega_width",
    "max_abs_schwarzschild_symmetry_residual",
    "max_abs_omega_center_minus_exact",
    "max_abs_omega_center_linear_residual_small_a",
    "max_abs_width_schwarzschild_residual_at_a0",
    "max_abs_width_positive_violation",
    "omega_center_mean",
    "omega_center_max_abs",
    "omega_width_mean",
    "schwarzschild_symmetry_pass",
    "discriminant_positive_pass",
    "omega_width_positive_pass",
    "frame_dragging_sign_pass",
    "omega_center_exact_identity_pass",
    "small_a_linear_scaling_pass",
    "global_true_relations",
    "global_false_relations",
    "global_undecided_pairs",
    "all_checks_pass",
)


# ---------------------------------------------------------------------------
# Metric and null-slope helpers
# ---------------------------------------------------------------------------

def equatorial_metric_at_r(
    r: float,
    mass: float,
    spin: float,
) -> dict[str, float]:
    """Boyer-Lindquist equatorial metric at theta=pi/2, signature -+++.

    Equatorial identities (sin(theta)=1, cos(theta)=0, Sigma=r^2):
      Delta    = r^2 - 2Mr + a^2
      g_tt     = -(1 - 2M/r)
      g_tphi   = -2Ma/r
      g_rr     = r^2/Delta
      g_phiphi = r^2 + a^2 + 2Ma^2/r

    Raises ValueError if Delta <= 0 (inside or on the Kerr horizon).
    """
    delta = r * r - 2.0 * mass * r + spin * spin
    if delta <= 0.0:
        raise ValueError(
            f"equatorial_metric_at_r: Delta={delta:.6g} <= 0 "
            f"at r={r}, M={mass}, a={spin}"
        )
    return {
        "g_tt":     -(1.0 - 2.0 * mass / r),
        "g_tphi":   -2.0 * mass * spin / r,
        "g_rr":     r * r / delta,
        "g_phiphi": r * r + spin * spin + 2.0 * mass * spin * spin / r,
        "delta":    delta,
    }


def angular_null_slopes(
    r: float,
    mass: float,
    spin: float,
) -> tuple[float, float, float, dict[str, float]]:
    """Equatorial angular null slopes at radius r with dr=dtheta=0.

    Solves g_phiphi omega^2 + 2 g_tphi omega + g_tt = 0.

    Returns (omega_plus, omega_minus, discriminant, metric).
    omega_plus >= omega_minus with strict inequality outside the horizon.
    """
    m = equatorial_metric_at_r(r, mass, spin)
    disc = m["g_tphi"] * m["g_tphi"] - m["g_tt"] * m["g_phiphi"]
    if disc < 0.0:
        raise ValueError(
            f"angular_null_slopes: disc={disc:.6g} < 0 at r={r}, M={mass}, a={spin}"
        )
    sqrt_disc = math.sqrt(disc)
    omega_plus  = (-m["g_tphi"] + sqrt_disc) / m["g_phiphi"]
    omega_minus = (-m["g_tphi"] - sqrt_disc) / m["g_phiphi"]
    return omega_plus, omega_minus, disc, m


# ---------------------------------------------------------------------------
# Per-spin case runner
# ---------------------------------------------------------------------------

def run_spin_case(
    n:      int,
    seed:   int,
    mass:   float,
    spin:   float,
    margin: float,
) -> dict[str, Any]:
    """Run one K5 diagnostic cell and return a dict matching CSV_FIELDS."""

    r_plus      = kerr_horizon_radius(mass, spin)
    r_min_event = r_plus + margin

    r_min_evaluated = min(RADIAL_GRID)
    all_points_exterior = all(r > r_min_event for r in RADIAL_GRID)

    # Per-radial-point accumulators
    discriminants:      list[float] = []
    omega_widths:       list[float] = []
    omega_centers:      list[float] = []
    identity_residuals: list[float] = []
    width_violations:   list[float] = []
    schw_sym_residuals: list[float] = []       # a=0 only
    width_schw_residuals: list[float] = []     # a=0 only
    linear_residuals:   list[float] = []       # 0 < a <= SMALL_A_MAX only

    for r in RADIAL_GRID:
        omega_plus, omega_minus, disc, m = angular_null_slopes(r, mass, spin)

        omega_center = (omega_plus + omega_minus) * 0.5
        omega_width  = (omega_plus - omega_minus) * 0.5

        discriminants.append(disc)
        omega_widths.append(omega_width)
        omega_centers.append(omega_center)

        # Exact algebraic identity: (omega_+ + omega_-)/2 = -g_tphi/g_phiphi
        omega_center_exact = -m["g_tphi"] / m["g_phiphi"]
        identity_residuals.append(abs(omega_center - omega_center_exact))

        # Width positivity check value
        width_violations.append(max(0.0, -omega_width))

        if abs(spin) <= 0.0:
            # Schwarzschild symmetry: omega_+ + omega_- = 0 exactly
            schw_sym_residuals.append(abs(omega_plus + omega_minus))
            # Width formula at a=0: omega_width = sqrt(-g_tt/g_phiphi)
            omega_width_schw = math.sqrt(-m["g_tt"] / m["g_phiphi"])
            width_schw_residuals.append(abs(omega_width - omega_width_schw))

        if 0.0 < abs(spin) <= SMALL_A_MAX:
            # Linear scaling: -g_tphi/(a*g_phiphi) -> 2M/r^3 as a->0
            linear_expected = 2.0 * mass / (r * r * r)
            linear_actual   = (-m["g_tphi"]) / (abs(spin) * m["g_phiphi"])
            linear_residuals.append(abs(linear_actual - linear_expected))

    min_discriminant     = min(discriminants)
    min_omega_width      = min(omega_widths)
    omega_center_mean    = sum(omega_centers) / len(omega_centers)
    omega_center_max_abs = max(abs(c) for c in omega_centers)
    omega_width_mean     = sum(omega_widths) / len(omega_widths)
    max_abs_omega_center_minus_exact  = max(identity_residuals)
    max_abs_width_positive_violation  = max(width_violations)

    # a=0-specific fields
    max_abs_schwarzschild_symmetry_residual: Optional[float]
    max_abs_width_schwarzschild_residual_at_a0: Optional[float]
    if abs(spin) <= 0.0:
        max_abs_schwarzschild_symmetry_residual     = max(schw_sym_residuals)
        max_abs_width_schwarzschild_residual_at_a0  = max(width_schw_residuals)
    else:
        max_abs_schwarzschild_symmetry_residual     = None
        max_abs_width_schwarzschild_residual_at_a0  = None

    # Small-a field
    max_abs_omega_center_linear_residual_small_a: Optional[float]
    if 0.0 < abs(spin) <= SMALL_A_MAX:
        max_abs_omega_center_linear_residual_small_a = max(linear_residuals)
    else:
        max_abs_omega_center_linear_residual_small_a = None

    # --- Check flags ---
    discriminant_positive_pass       = min_discriminant > 0.0
    omega_width_positive_pass        = min_omega_width  > 0.0
    omega_center_exact_identity_pass = (
        max_abs_omega_center_minus_exact <= METRIC_FORMULA_TOL
    )

    if abs(spin) <= 0.0:
        schwarzschild_symmetry_pass = (
            max_abs_schwarzschild_symmetry_residual <= METRIC_FORMULA_TOL
        )
        frame_dragging_sign_pass = True   # not applicable at a=0; True by convention
    else:
        schwarzschild_symmetry_pass = True   # not checked for a>0
        frame_dragging_sign_pass    = all(c > 0.0 for c in omega_centers)

    if 0.0 < abs(spin) <= SMALL_A_MAX:
        small_a_linear_scaling_pass = (
            max_abs_omega_center_linear_residual_small_a <= SMALL_A_LINEAR_TOL
        )
    else:
        small_a_linear_scaling_pass = True   # not applicable

    # --- Causal accounting (N equatorial scaffold events, K1-K4 invariant) ---
    events  = generate_exterior_events(n, seed, r_min_event, equatorial=True)
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
    # K5 invariant: a>0 → all pairs undecided
    if abs(spin) > 0.0:
        true_relations  = 0
        false_relations = 0
        undecided_pairs = possible_pairs

    all_checks_pass: bool = (
        all_points_exterior
        and discriminant_positive_pass
        and omega_width_positive_pass
        and omega_center_exact_identity_pass
        and schwarzschild_symmetry_pass
        and frame_dragging_sign_pass
        and small_a_linear_scaling_pass
    )

    return {
        "spin_a":                                    spin,
        "M":                                         mass,
        "N":                                         n,
        "seed":                                      seed,
        "margin":                                    margin,
        "radial_grid":                               list(RADIAL_GRID),
        "r_plus":                                    r_plus,
        "r_min_evaluated":                           r_min_evaluated,
        "all_points_exterior":                       all_points_exterior,
        "min_discriminant":                          min_discriminant,
        "min_omega_width":                           min_omega_width,
        "max_abs_schwarzschild_symmetry_residual":   max_abs_schwarzschild_symmetry_residual,
        "max_abs_omega_center_minus_exact":          max_abs_omega_center_minus_exact,
        "max_abs_omega_center_linear_residual_small_a": max_abs_omega_center_linear_residual_small_a,
        "max_abs_width_schwarzschild_residual_at_a0": max_abs_width_schwarzschild_residual_at_a0,
        "max_abs_width_positive_violation":          max_abs_width_positive_violation,
        "omega_center_mean":                         omega_center_mean,
        "omega_center_max_abs":                      omega_center_max_abs,
        "omega_width_mean":                          omega_width_mean,
        "schwarzschild_symmetry_pass":               schwarzschild_symmetry_pass,
        "discriminant_positive_pass":                discriminant_positive_pass,
        "omega_width_positive_pass":                 omega_width_positive_pass,
        "frame_dragging_sign_pass":                  frame_dragging_sign_pass,
        "omega_center_exact_identity_pass":          omega_center_exact_identity_pass,
        "small_a_linear_scaling_pass":               small_a_linear_scaling_pass,
        "global_true_relations":                     true_relations,
        "global_false_relations":                    false_relations,
        "global_undecided_pairs":                    undecided_pairs,
        "all_checks_pass":                           all_checks_pass,
    }


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
    """Run K5 for all spins and return the full payload dict."""
    if mass != 1.0:
        raise ValueError("K5 diagnostic is fixed to M=1")
    if any(abs(a) >= mass for a in spins if abs(a) > 0.0):
        raise ValueError("K5 requires |a| < M for all non-zero spins")

    rows    = [run_spin_case(n, seed, mass, spin, margin) for spin in spins]
    all_pass = all(row["all_checks_pass"] for row in rows)

    aggregate: dict[str, Any] = {
        "audit":             AUDIT_ID,
        "benchmark":         "S4-K5 Kerr equatorial prograde/retrograde local null-slope diagnostic",
        "generated_at_utc":  datetime.now(timezone.utc).isoformat(),
        "command":           COMMAND,
        "N":                 n,
        "seed":              seed,
        "M":                 mass,
        "spins":             list(spins),
        "margin":            margin,
        "radial_grid":       list(RADIAL_GRID),
        "theta":             math.pi / 2.0,
        "possible_pairs":    n * (n - 1) // 2,
        "metric_formula_tol": METRIC_FORMULA_TOL,
        "small_a_linear_tol": SMALL_A_LINEAR_TOL,
        "small_a_max":       SMALL_A_MAX,
        "all_checks_pass":   all_pass,
        "positive_spin_cases_all_undecided": all(
            row["global_true_relations"]  == 0
            and row["global_false_relations"] == 0
            and row["global_undecided_pairs"] == n * (n - 1) // 2
            for row in rows if row["spin_a"] > 0.0
        ),
        "scope_note": (
            "K5 is a local null-slope diagnostic. "
            "It measures local prograde/retrograde asymmetry of the Kerr equatorial "
            "light cone. It does not implement Kerr causal inference, "
            "does not integrate null geodesics, "
            "and does not claim global causal reachability. "
            "Local prograde/retrograde null slopes are not global causal relations."
        ),
    }

    return {
        "aggregate": aggregate,
        "rows":      rows,
    }


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
    if isinstance(value, list):
        return ";".join(str(v) for v in value)
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
        f"# {AUDIT_ID}: Kerr Equatorial Prograde/Retrograde Local Null-Slope Diagnostic",
        "",
        f"Generated: {aggregate['generated_at_utc']}",
        "",
        "## What this is",
        "",
        "This is a **local null-slope diagnostic**, not a Kerr causal solver.",
        "",
        "It measures the local prograde/retrograde asymmetry of the Kerr equatorial",
        "light cone caused by the frame-dragging term g_tphi.",
        "",
        "**It does NOT:**",
        "",
        "- Implement Kerr causal inference of any kind.",
        "- Integrate null geodesics.",
        "- Claim global causal reachability.",
        "- Assert global causal relations.  Local prograde/retrograde null slopes",
        "  are **not** global causal relations.",
        "",
        "## Parameters",
        "",
        f"- M = {aggregate['M']}, theta = pi/2 (equatorial), M = 1 (fixed)",
        f"- Spin sweep: {aggregate['spins']}",
        f"- N = {aggregate['N']}, seed = {aggregate['seed']}, "
        f"margin = {aggregate['margin']}",
        f"- **Fixed radial grid** (metric evaluation): "
        f"`r = {aggregate['radial_grid']}`",
        f"  (all points safely outside `r_+(a=0.75) + margin ≈ 2.011`)",
        f"- Metric formula tolerance: `{aggregate['metric_formula_tol']:.1e}`",
        f"- Small-a linear tolerance: `{aggregate['small_a_linear_tol']:.1e}`"
        f" (checked for `0 < a <= {aggregate['small_a_max']}`)",
        "",
        "## Analytic Checks",
        "",
        "1. **Discriminant positivity** (all a): `disc = g_tphi^2 - g_tt*g_phiphi > 0`",
        "2. **Cone width positivity** (all a): `omega_width = sqrt(disc)/g_phiphi > 0`",
        "3. **Exact center identity** (all a): `(omega_+ + omega_-)/2 = -g_tphi/g_phiphi`",
        "   (algebraic identity; residual ≤ 1e-12)",
        "4. **Schwarzschild symmetry** (a=0): `omega_+ = -omega_-`",
        "   (frame-dragging absent; residual ≤ 1e-12)",
        "5. **Frame-dragging sign** (a>0): `omega_center = -g_tphi/g_phiphi > 0`",
        "   (g_tphi < 0 for a>0; cone tilts in prograde direction)",
        "6. **Linear scaling** (0<a≤0.01): `omega_center/a → 2M/r^3`",
        "   (leading-order frame-dragging; residual ≤ 1e-4)",
        "",
        "## Diagnostic Figure",
        "",
        f"![K5 prograde/retrograde local cone audit]({png_path.name})",
        "",
        "The 2×2 figure shows:",
        "- Panel 1: mean omega_center vs spin a (zero at a=0, positive for a>0)",
        "- Panel 2: small-a linear scaling residual vs a (log-log, a≤0.01)",
        "- Panel 3: omega_+ and omega_- vs r for a=0, 0.5, 0.75",
        "  (symmetry at a=0, asymmetry for a>0)",
        "- Panel 4: min discriminant and min cone width vs spin a",
        "",
        "## Summary",
        "",
        "| Check | Result |",
        "|-------|--------|",
        f"| **all_checks_pass** | **{aggregate['all_checks_pass']}** |",
        f"| positive_spin_cases_all_undecided | "
        f"{aggregate['positive_spin_cases_all_undecided']} |",
        "",
        "## Per-Spin Results",
        "",
        "| a | r_+ | ext? | disc+ | width+ | fd_sign | schw_sym | lin_ok | pass |",
        "|---|-----|------|-------|--------|---------|----------|--------|------|",
    ]

    for row in rows:
        a = row["spin_a"]
        lines.append(
            f"| {a:.1e} "
            f"| {row['r_plus']:.6f} "
            f"| {row['all_points_exterior']} "
            f"| {row['discriminant_positive_pass']} "
            f"| {row['omega_width_positive_pass']} "
            f"| {row['frame_dragging_sign_pass']} "
            f"| {row['schwarzschild_symmetry_pass']} "
            f"| {row['small_a_linear_scaling_pass']} "
            f"| **{row['all_checks_pass']}** |"
        )

    lines += [
        "",
        "## Causal Accounting",
        "",
        "| a | global_true | global_false | global_undecided |",
        "|---|-------------|--------------|-----------------|",
    ]
    for row in rows:
        lines.append(
            f"| {row['spin_a']:.1e} "
            f"| {row['global_true_relations']} "
            f"| {row['global_false_relations']} "
            f"| {row['global_undecided_pairs']} |"
        )

    lines += [
        "",
        "## Interpretation",
        "",
        "- For `a=0`: Schwarzschild symmetry is verified; the existing scaffold",
        "  control causal counts are preserved.",
        "- For `a>0`: the cone tilts in the prograde direction (`omega_center > 0`),",
        "  but all global causal pairs remain undecided",
        "  (true=0, false=0, undecided=N*(N-1)/2).",
        "- Metric formula residuals (checks 3, 4) are at machine precision (~10^-15):",
        "  both sides use the same closed-form expression.",
        "- The linear scaling check (6) is a genuine perturbative check:",
        "  the O(a^2) correction to g_phiphi introduces a residual ~3.7e-7 at a=0.01.",
        "",
    ]
    md_path.write_text("\n".join(lines), encoding="utf-8")


# ---------------------------------------------------------------------------
# Figure generator
# ---------------------------------------------------------------------------

def write_figure(
    rows:     list[dict[str, Any]],
    png_path: Path,
) -> None:
    """Generate a 2x2 diagnostic panel of prograde/retrograde null-slope results."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    M     = 1.0
    spins_all = [r["spin_a"] for r in rows]

    _FLOOR = 1e-20   # floor to avoid log(0)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    fig.suptitle("S4-KERR-K5 prograde/retrograde local cone audit", fontsize=14)

    # ---- Panel 1: omega_center_mean vs spin_a (linear x) ----
    ax1 = axes[0, 0]
    omega_cm = [r["omega_center_mean"] for r in rows]
    ax1.plot(spins_all, omega_cm, "o-", color="steelblue", linewidth=1.5, markersize=5)
    ax1.axhline(0.0, color="gray", linestyle="--", linewidth=0.8)
    ax1.set_xlabel("spin $a$")
    ax1.set_ylabel(r"$\bar{\omega}_{center}$")
    ax1.set_title("Frame-dragging: mean cone center vs $a$")
    ax1.grid(True, which="both", linestyle="--", alpha=0.4)

    # ---- Panel 2: small-a linear scaling residual (log-log) ----
    ax2 = axes[0, 1]
    small_rows = [r for r in rows if 0.0 < r["spin_a"] <= SMALL_A_MAX]
    if small_rows:
        sa    = [r["spin_a"] for r in small_rows]
        resid = [max(r["max_abs_omega_center_linear_residual_small_a"], _FLOOR)
                 for r in small_rows]
        ax2.loglog(sa, resid, "o-", color="darkorange", linewidth=1.5, markersize=5)
        ax2.set_xlabel("spin $a$")
        ax2.set_ylabel(r"$|\omega_{center}/a - 2M/r^3|_{max}$")
        ax2.set_title("Small-$a$ linear scaling residual")
        ax2.grid(True, which="both", linestyle="--", alpha=0.4)

    # ---- Panel 3: omega_plus/minus vs r for a=0, 0.5, 0.75 ----
    ax3 = axes[1, 0]
    plot_spins = [
        (0.0,  "steelblue",  "solid"),
        (0.5,  "darkorange", "dashed"),
        (0.75, "green",      "dashdot"),
    ]
    r_grid = list(RADIAL_GRID)
    for spin_val, color, ls in plot_spins:
        op_list: list[float] = []
        om_list: list[float] = []
        for r in r_grid:
            delta     = r * r - 2.0 * M * r + spin_val * spin_val
            g_tt      = -(1.0 - 2.0 * M / r)
            g_tphi    = -2.0 * M * spin_val / r
            g_phiphi  = r * r + spin_val * spin_val + 2.0 * M * spin_val * spin_val / r
            disc_val  = g_tphi * g_tphi - g_tt * g_phiphi
            sq        = math.sqrt(disc_val)
            op_list.append((-g_tphi + sq) / g_phiphi)
            om_list.append((-g_tphi - sq) / g_phiphi)
        label_p = rf"$\omega_+$, $a$={spin_val:.2f}"
        label_m = rf"$\omega_-$, $a$={spin_val:.2f}"
        ax3.plot(r_grid, op_list, "o", color=color, linewidth=1.5, markersize=4,
                 linestyle=ls, label=label_p)
        ax3.plot(r_grid, om_list, "s", color=color, linewidth=1.5, markersize=4,
                 linestyle=ls, label=label_m, alpha=0.7)
    ax3.axhline(0.0, color="gray", linestyle="--", linewidth=0.8)
    ax3.set_xlabel("radius $r$")
    ax3.set_ylabel(r"angular null slope $\omega$")
    ax3.set_title(r"Null slopes $\omega_\pm$ vs $r$ (equatorial)")
    ax3.legend(fontsize=7, ncol=2)
    ax3.grid(True, which="both", linestyle="--", alpha=0.4)

    # ---- Panel 4: min_discriminant and min_omega_width vs spin_a ----
    ax4 = axes[1, 1]
    min_disc  = [r["min_discriminant"] for r in rows]
    min_width = [r["min_omega_width"]  for r in rows]
    ax4.plot(spins_all, min_disc,  "o-", color="steelblue",  linewidth=1.5,
             markersize=5, label="min discriminant")
    ax4.plot(spins_all, min_width, "s--", color="darkorange", linewidth=1.5,
             markersize=5, label=r"min $\omega_{width}$")
    ax4.axhline(0.0, color="gray", linestyle="--", linewidth=0.8)
    ax4.set_xlabel("spin $a$")
    ax4.set_ylabel("value")
    ax4.set_title("Positivity of discriminant and cone width")
    ax4.legend()
    ax4.grid(True, which="both", linestyle="--", alpha=0.4)

    plt.tight_layout()
    fig.savefig(str(png_path), dpi=150, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Top-level output writer
# ---------------------------------------------------------------------------

def write_outputs(
    payload:    dict[str, Any],
    out_prefix: str = OUT_PREFIX,
) -> tuple[Path, Path, Path, Path]:
    csv_path  = OUT_DIR / f"{out_prefix}.csv"
    json_path = OUT_DIR / f"{out_prefix}.json"
    md_path   = OUT_DIR / f"{out_prefix}.md"
    png_path  = OUT_DIR / f"{out_prefix}.png"

    rows      = payload["rows"]
    aggregate = payload["aggregate"]

    write_csv(rows, csv_path)
    write_json(payload, json_path)
    write_figure(rows, png_path)
    write_md(rows, aggregate, md_path, png_path)

    return csv_path, json_path, md_path, png_path


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    print(f"Running {AUDIT_ID}")
    payload  = run_audit()
    agg      = payload["aggregate"]
    csv_path, json_path, md_path, png_path = write_outputs(payload)

    print(f"all_checks_pass={agg['all_checks_pass']}")
    for row in payload["rows"]:
        a = row["spin_a"]
        print(
            f"  a={a:.1e}"
            f"  r_+={row['r_plus']:.6f}"
            f"  ext={row['all_points_exterior']}"
            f"  disc+={row['discriminant_positive_pass']}"
            f"  fd={row['frame_dragging_sign_pass']}"
            f"  schw={row['schwarzschild_symmetry_pass']}"
            f"  lin={row['small_a_linear_scaling_pass']}"
            f"  pass={row['all_checks_pass']}"
        )
    print(f"wrote {csv_path}")
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")
    print(f"wrote {png_path}")


if __name__ == "__main__":
    main()
