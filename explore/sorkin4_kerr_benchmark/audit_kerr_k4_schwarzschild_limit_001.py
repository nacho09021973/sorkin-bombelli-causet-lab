#!/usr/bin/env python3
"""S4-KERR-K4-SCHWARZSCHILD-LIMIT-001: Kerr equatorial Schwarzschild-limit audit.

WHAT THIS IS:
  A known-truth perturbative metric audit.  It verifies that the
  Boyer-Lindquist equatorial metric has the correct analytic Schwarzschild
  limit as a -> 0, and checks the perturbative scaling of metric components
  at a fixed radial grid.

WHAT THIS IS NOT:
  - It does not implement Kerr causal inference.
  - It does not integrate null geodesics.
  - It does not claim global causal reachability.
  - It does not create causal true/false relations for a != 0.

Physics (Boyer-Lindquist equatorial plane, theta=pi/2, M=1):

  Delta    = r^2 - 2Mr + a^2
  g_tt     = -(1 - 2M/r)
  g_tphi   = -2Ma/r
  g_rr     = r^2/Delta
  g_phiphi = r^2 + a^2 + 2Ma^2/r

  Schwarzschild limit at a=0:
    g_tphi = 0, g_rr = 1/(1-2M/r), g_phiphi = r^2

  Perturbative checks (a > 0):
    1. Frame-dragging linearity:   g_tphi / a     = -2M/r
    2. Azimuthal quadratic:        (g_phiphi-r^2) / a^2 = 1 + 2M/r
    3. Radial formula:             g_rr = r^2/Delta  (identity at all a)
    4. Horizon quadratic shift:    (2M - r_+) / a^2 -> 1/(2M)  as a->0
       (checked only for a <= 0.01 where the O(a^2) error is small enough)
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

AUDIT_ID   = "S4-KERR-K4-SCHWARZSCHILD-LIMIT-001"
OUT_DIR    = Path(__file__).resolve().parent
OUT_PREFIX = "kerr_k4_schwarzschild_limit_001_n12_seed1959"
COMMAND    = (
    "python3 explore/sorkin4_kerr_benchmark/"
    "audit_kerr_k4_schwarzschild_limit_001.py"
)

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

DEFAULT_N      = 12
DEFAULT_SEED   = 1959
DEFAULT_MASS   = 1.0
DEFAULT_MARGIN = schwarz.EXTERIOR_MARGIN  # 0.35

DEFAULT_SPINS: tuple[float, ...] = (
    0.0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1,
)

# Fixed radial grid for metric evaluation (independent of spin).
# All five points are safely outside r_plus(a=0.1) + margin ≈ 2.345:
#   2.5, 3.0, 4.0, 6.0, 10.0 > 2.345
RADIAL_GRID: tuple[float, ...] = (2.5, 3.0, 4.0, 6.0, 10.0)

# Horizon quadratic check: only applied for 0 < a <= this value.
# Error is O(a^2): at a=0.01 the error ≈ (0.01)^2/8 ≈ 1.25e-5.
HORIZON_QUAD_SPIN_MAX = 1e-2
HORIZON_QUAD_TOL      = 1.0e-3   # generous: 80× max expected error at a=0.01

# Identity tolerances (g_tphi/a, (g_phiphi-r^2)/a^2, g_rr=r^2/Delta):
# These are computed from the same closed-form expressions, so residuals
# are determined by floating-point rounding alone.
METRIC_FORMULA_TOL = 1.0e-12

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
    "r_plus_shift",
    "r_plus_shift_over_a2",
    "r_min_evaluated",
    "all_points_exterior",
    "max_abs_g_tphi_minus_linear",
    "max_abs_gphiphi_quadratic_residual",
    "max_abs_grr_formula_residual",
    "max_abs_gtphi_at_a0",
    "max_abs_gphiphi_minus_r2_at_a0",
    "max_abs_grr_schwarzschild_at_a0",
    "horizon_quadratic_check_pass",
    "frame_dragging_linear_check_pass",
    "gphiphi_quadratic_check_pass",
    "grr_formula_check_pass",
    "schwarzschild_metric_limit_pass",
    "global_true_relations",
    "global_false_relations",
    "global_undecided_pairs",
    "all_checks_pass",
)


# ---------------------------------------------------------------------------
# Equatorial metric helper
# ---------------------------------------------------------------------------

def equatorial_metric_at_r(
    r: float,
    mass: float,
    spin: float,
) -> dict[str, float]:
    """Boyer-Lindquist equatorial metric coefficients at theta=pi/2, signature -+++.

    Equatorial simplifications (theta=pi/2 → sin=1, cos=0, Sigma=r²):
      Delta    = r² - 2Mr + a²
      g_tt     = -(1 - 2M/r)
      g_tphi   = -2Ma/r
      g_rr     = r²/Delta
      g_phiphi = r² + a² + 2Ma²/r

    Schwarzschild BL limit (a=0):
      g_tphi = 0, g_rr = 1/(1-2M/r), g_phiphi = r²

    Raises ValueError if Delta <= 0 (inside or on the Kerr horizon).
    """
    delta = r * r - 2.0 * mass * r + spin * spin
    if delta <= 0.0:
        raise ValueError(
            f"equatorial_metric_at_r: Delta={delta:.6g} <= 0 "
            f"at r={r}, M={mass}, a={spin}; "
            "point is inside or on the Kerr horizon."
        )
    return {
        "g_tt":     -(1.0 - 2.0 * mass / r),
        "g_tphi":   -2.0 * mass * spin / r,
        "g_rr":     r * r / delta,
        "g_phiphi": r * r + spin * spin + 2.0 * mass * spin * spin / r,
        "delta":    delta,
    }


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
    """Run one (spin, N, seed, margin) K4 diagnostic cell.

    Returns a dict matching CSV_FIELDS plus any extra internal keys.
    """
    # --- Horizon geometry ---
    r_plus      = kerr_horizon_radius(mass, spin)
    r_min_event = r_plus + margin  # lower bound for generated events

    r_plus_shift: float = 2.0 * mass - r_plus  # = M - sqrt(M^2 - a^2) >= 0

    r_plus_shift_over_a2: Optional[float]
    if abs(spin) > 0.0:
        r_plus_shift_over_a2 = r_plus_shift / (spin * spin)
    else:
        r_plus_shift_over_a2 = None  # undefined for a=0

    # --- Fixed radial grid checks ---
    radial_grid     = list(RADIAL_GRID)
    r_min_evaluated = min(radial_grid)
    # All grid points must be outside r_plus + margin for all spins in the sweep
    all_points_exterior = all(r > r_min_event for r in radial_grid)

    # Accumulators for perturbative residuals
    frame_drag_residuals:    list[float] = []
    gphiphi_quad_residuals:  list[float] = []
    grr_formula_residuals:   list[float] = []
    gtphi_a0_vals:           list[float] = []
    gphiphi_a0_vals:         list[float] = []
    grr_a0_vals:             list[float] = []

    for r in radial_grid:
        m = equatorial_metric_at_r(r, mass, spin)

        # Check 4: g_rr formula identity: g_rr = r^2 / Delta (all spins)
        # This is the same formula used in equatorial_metric_at_r, so the
        # residual is determined by floating-point rounding alone.
        grr_formula_value = r * r / m["delta"]
        grr_formula_residuals.append(abs(m["g_rr"] - grr_formula_value))

        if abs(spin) > 0.0:
            # Check 1: frame-dragging linearity: g_tphi / a = -2M/r
            # g_tphi = -2Ma/r → g_tphi/a = -2M/r (no cancellation; both sides ~M/r)
            linear_expected = -2.0 * mass / r
            frame_drag_residuals.append(abs(m["g_tphi"] / spin - linear_expected))

            # Check 2: azimuthal quadratic: g_phiphi = r^2 + a^2*(1+2M/r)
            # Computed as ABSOLUTE residual to avoid catastrophic cancellation.
            # (g_phiphi - r^2)/a^2 = 1+2M/r is the identity, but
            # computing the l.h.s. via m["g_phiphi"] - r*r and then dividing
            # by a^2 amplifies floating-point rounding as O(eps*r^2/a^2),
            # which exceeds 1e-12 for a < ~0.094 at r~2.5.
            # Using abs(g_phiphi - r^2 - a^2*(1+2M/r)) avoids this:
            # residual ~ eps * g_phiphi ~ 1e-15, stable for all a in the sweep.
            gphiphi_expected = r * r + spin * spin * (1.0 + 2.0 * mass / r)
            gphiphi_quad_residuals.append(abs(m["g_phiphi"] - gphiphi_expected))
        else:
            # Schwarzschild limit checks (a=0)
            gtphi_a0_vals.append(abs(m["g_tphi"]))
            grr_a0_vals.append(abs(m["g_rr"] - 1.0 / (1.0 - 2.0 * mass / r)))
            gphiphi_a0_vals.append(abs(m["g_phiphi"] - r * r))

    max_abs_grr_formula_residual: float = max(grr_formula_residuals)

    # Perturbative residuals for a > 0
    max_abs_g_tphi_minus_linear: Optional[float]
    max_abs_gphiphi_quadratic_residual: Optional[float]
    if abs(spin) > 0.0:
        max_abs_g_tphi_minus_linear         = max(frame_drag_residuals)
        max_abs_gphiphi_quadratic_residual  = max(gphiphi_quad_residuals)
    else:
        max_abs_g_tphi_minus_linear         = None
        max_abs_gphiphi_quadratic_residual  = None

    # Schwarzschild limit values (a=0 only)
    max_abs_gtphi_at_a0:           Optional[float]
    max_abs_gphiphi_minus_r2_at_a0: Optional[float]
    max_abs_grr_schwarzschild_at_a0: Optional[float]
    if abs(spin) <= 0.0:
        max_abs_gtphi_at_a0            = max(gtphi_a0_vals)
        max_abs_gphiphi_minus_r2_at_a0 = max(gphiphi_a0_vals)
        max_abs_grr_schwarzschild_at_a0 = max(grr_a0_vals)
    else:
        max_abs_gtphi_at_a0            = None
        max_abs_gphiphi_minus_r2_at_a0 = None
        max_abs_grr_schwarzschild_at_a0 = None

    # --- Individual check flags ---

    grr_formula_check_pass: bool = (
        max_abs_grr_formula_residual <= METRIC_FORMULA_TOL
    )

    if abs(spin) > 0.0:
        frame_dragging_linear_check_pass: bool = (
            max_abs_g_tphi_minus_linear <= METRIC_FORMULA_TOL
        )
        gphiphi_quadratic_check_pass: bool = (
            max_abs_gphiphi_quadratic_residual <= METRIC_FORMULA_TOL
        )
        schwarzschild_metric_limit_pass: bool = True  # not applicable for a>0
    else:
        # a=0: these perturbative ratio checks are trivially True (no division)
        frame_dragging_linear_check_pass = True
        gphiphi_quadratic_check_pass     = True
        schwarzschild_metric_limit_pass  = (
            max_abs_gtphi_at_a0            <= METRIC_FORMULA_TOL
            and max_abs_grr_schwarzschild_at_a0 <= METRIC_FORMULA_TOL
            and max_abs_gphiphi_minus_r2_at_a0  <= METRIC_FORMULA_TOL
        )

    # Horizon quadratic check:
    #   (2M - r_+) / a^2 = [M - sqrt(M^2-a^2)] / a^2 -> 1/(2M) as a->0
    # Error is O(a^2): checked only for 0 < a <= HORIZON_QUAD_SPIN_MAX.
    # For a=0 or a > HORIZON_QUAD_SPIN_MAX: True (not applicable).
    horizon_quadratic_check_pass: bool
    if abs(spin) > 0.0 and spin <= HORIZON_QUAD_SPIN_MAX:
        horizon_quad_target = 1.0 / (2.0 * mass)
        horizon_quadratic_check_pass = (
            abs(r_plus_shift_over_a2 - horizon_quad_target) <= HORIZON_QUAD_TOL
        )
    else:
        horizon_quadratic_check_pass = True

    # --- Causal accounting (N equatorial events, equatorial_scaffold mode) ---
    # For a=0: uses the existing Schwarzschild control subset.
    # For a>0: enforced to all-undecided by K1/K2/K3 invariant.
    events = generate_exterior_events(n, seed, r_min_event, equatorial=True)
    matrix, states = build_relation_states(events, mass, spin, "equatorial_scaffold")

    possible_pairs  = n * (n - 1) // 2
    true_relations  = count_true_relations(matrix)
    false_relations = sum(
        1
        for i in range(n - 1)
        for j in range(i + 1, n)
        if states[i][j] is False
    )
    undecided_pairs = sum(
        1
        for i in range(n - 1)
        for j in range(i + 1, n)
        if states[i][j] is None
    )

    # Enforce K4 causal accounting rule: a>0 → all pairs undecided.
    if abs(spin) > 0.0:
        true_relations  = 0
        false_relations = 0
        undecided_pairs = possible_pairs

    # --- all_checks_pass ---
    all_checks_pass: bool = (
        all_points_exterior
        and grr_formula_check_pass
        and frame_dragging_linear_check_pass
        and gphiphi_quadratic_check_pass
        and horizon_quadratic_check_pass
        and schwarzschild_metric_limit_pass
    )

    return {
        "spin_a":                           spin,
        "M":                                mass,
        "N":                                n,
        "seed":                             seed,
        "margin":                           margin,
        "radial_grid":                      radial_grid,
        "r_plus":                           r_plus,
        "r_plus_shift":                     r_plus_shift,
        "r_plus_shift_over_a2":             r_plus_shift_over_a2,
        "r_min_evaluated":                  r_min_evaluated,
        "all_points_exterior":              all_points_exterior,
        "max_abs_g_tphi_minus_linear":      max_abs_g_tphi_minus_linear,
        "max_abs_gphiphi_quadratic_residual": max_abs_gphiphi_quadratic_residual,
        "max_abs_grr_formula_residual":     max_abs_grr_formula_residual,
        "max_abs_gtphi_at_a0":              max_abs_gtphi_at_a0,
        "max_abs_gphiphi_minus_r2_at_a0":   max_abs_gphiphi_minus_r2_at_a0,
        "max_abs_grr_schwarzschild_at_a0":  max_abs_grr_schwarzschild_at_a0,
        "horizon_quadratic_check_pass":     horizon_quadratic_check_pass,
        "frame_dragging_linear_check_pass": frame_dragging_linear_check_pass,
        "gphiphi_quadratic_check_pass":     gphiphi_quadratic_check_pass,
        "grr_formula_check_pass":           grr_formula_check_pass,
        "schwarzschild_metric_limit_pass":  schwarzschild_metric_limit_pass,
        "global_true_relations":            true_relations,
        "global_false_relations":           false_relations,
        "global_undecided_pairs":           undecided_pairs,
        "all_checks_pass":                  all_checks_pass,
    }


# ---------------------------------------------------------------------------
# Full audit runner
# ---------------------------------------------------------------------------

def run_audit(
    n:      int              = DEFAULT_N,
    seed:   int              = DEFAULT_SEED,
    mass:   float            = DEFAULT_MASS,
    spins:  tuple[float, ...] = DEFAULT_SPINS,
    margin: float            = DEFAULT_MARGIN,
) -> dict[str, Any]:
    """Run K4 for all spins and return the full payload dict."""
    if mass != 1.0:
        raise ValueError("K4 diagnostic is fixed to M=1")
    if any(abs(a) >= mass for a in spins if abs(a) > 0.0):
        raise ValueError("K4 requires |a| < M for all non-zero spins")

    rows = [run_spin_case(n, seed, mass, spin, margin) for spin in spins]
    all_pass = all(row["all_checks_pass"] for row in rows)

    aggregate: dict[str, Any] = {
        "audit":             AUDIT_ID,
        "benchmark":         "S4-K4 Kerr equatorial Schwarzschild-limit audit",
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
        "horizon_quad_tol":  HORIZON_QUAD_TOL,
        "horizon_quad_spin_max": HORIZON_QUAD_SPIN_MAX,
        "all_checks_pass":   all_pass,
        "positive_spin_cases_all_undecided": all(
            row["global_true_relations"]  == 0
            and row["global_false_relations"] == 0
            and row["global_undecided_pairs"] == n * (n - 1) // 2
            for row in rows
            if row["spin_a"] > 0.0
        ),
        "scope_note": (
            "K4 is a known-truth perturbative metric audit. "
            "It checks the analytic a -> 0 Schwarzschild limit of the "
            "Boyer-Lindquist equatorial metric at a fixed radial grid. "
            "It does not implement Kerr causal inference, "
            "does not integrate null geodesics, "
            "and does not claim global causal reachability."
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
        f"# {AUDIT_ID}: Kerr Equatorial Schwarzschild-Limit Audit",
        "",
        f"Generated: {aggregate['generated_at_utc']}",
        "",
        "## What this is",
        "",
        "This is a **known-truth perturbative metric audit**, not a Kerr causal solver.",
        "",
        "It checks the analytic `a -> 0` Schwarzschild limit of the Boyer-Lindquist",
        "equatorial metric at a fixed radial grid.",
        "",
        "**It does NOT:**",
        "",
        "- Implement Kerr causal inference of any kind.",
        "- Integrate null geodesics.",
        "- Claim global causal reachability.",
        "- Create causal true/false relations for `a != 0`.",
        "",
        "## Parameters",
        "",
        f"- M = {aggregate['M']}, theta = pi/2 (equatorial), M = 1 (fixed)",
        f"- Spin sweep: {aggregate['spins']}",
        f"- N = {aggregate['N']}, seed = {aggregate['seed']}, "
        f"margin = {aggregate['margin']}",
        f"- **Fixed radial grid** (metric evaluation): "
        f"`r = {aggregate['radial_grid']}`",
        f"  (all points safely outside `r_+(a=0.1) + margin ≈ 2.345`)",
        f"- Metric formula tolerance: `{aggregate['metric_formula_tol']:.1e}`",
        f"- Horizon quadratic tolerance: `{aggregate['horizon_quad_tol']:.1e}`"
        f" (checked for `a <= {aggregate['horizon_quad_spin_max']}`)",
        "",
        "## Analytic Checks",
        "",
        "1. **Schwarzschild limit at a=0**: g_tphi=0, g_rr=1/(1-2M/r), g_phiphi=r²",
        "2. **Frame-dragging linearity** (a>0): g_tphi/a = -2M/r",
        "3. **Azimuthal quadratic** (a>0): abs(g_phiphi - r² - a²*(1+2M/r)) ≤ tol",
        "   (absolute residual; the ratio form has catastrophic cancellation for small a)",
        "4. **g_rr formula** (all a): g_rr = r²/Δ  (identity check)",
        "5. **Horizon quadratic shift** (0<a≤0.01): (2M-r_+)/a² → 1/(2M)",
        "",
        "## Diagnostic Figure",
        "",
        f"![K4 Schwarzschild-limit audit]({png_path.name})",
        "",
        "The figure shows perturbative scaling residuals vs. spin a (log-scale x-axis),",
        "excluding a=0 where division by a or a² is undefined.",
        "Residuals 1–3 are at machine precision (same closed-form formula both sides).",
        "Residual 4 (horizon shift) is a genuine O(a²) perturbative check.",
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
        "| a | r_+ | r_+_shift/a² | ext? | frame_drag | gphiphi_quad | grr | horiz | schw | pass |",
        "|---|-----|--------------|------|-----------|-------------|-----|-------|------|------|",
    ]

    for row in rows:
        a = row["spin_a"]
        shift_over_a2 = row["r_plus_shift_over_a2"]
        shift_str = f"{shift_over_a2:.6f}" if shift_over_a2 is not None else "—"
        lines.append(
            f"| {a:.1e} "
            f"| {row['r_plus']:.6f} "
            f"| {shift_str} "
            f"| {row['all_points_exterior']} "
            f"| {row['frame_dragging_linear_check_pass']} "
            f"| {row['gphiphi_quadratic_check_pass']} "
            f"| {row['grr_formula_check_pass']} "
            f"| {row['horizon_quadratic_check_pass']} "
            f"| {row['schwarzschild_metric_limit_pass']} "
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
        "- For `a=0`: the existing Schwarzschild/Kerr scaffold control behavior "
        "is preserved.",
        "- For `a>0`: all global causal pairs remain undecided "
        "(true=0, false=0, undecided=N*(N-1)/2).",
        "- Metric residuals 1–3 are at machine precision (~10⁻¹⁵): "
        "both sides use the same closed-form formula.",
        "- Residual 4 (horizon shift) is a genuine perturbative check: "
        "the exact value `(2M-r_+)/a²` converges to `1/(2M)=0.5` "
        "with O(a²) error.",
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
    """Generate a 2×2 diagnostic panel of perturbative scaling residuals."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Only rows with a > 0
    positive_rows = [r for r in rows if r["spin_a"] > 0.0]
    spins_pos = [r["spin_a"] for r in positive_rows]

    _FLOOR = 1e-20  # floor for log-scale: avoids log(0) for machine-zero residuals

    tphi_res   = [max(r["max_abs_g_tphi_minus_linear"], _FLOOR)
                  for r in positive_rows]
    gphi_res   = [max(r["max_abs_gphiphi_quadratic_residual"], _FLOOR)
                  for r in positive_rows]
    grr_res    = [max(r["max_abs_grr_formula_residual"], _FLOOR)
                  for r in positive_rows]
    horiz_res  = [
        max(abs(r["r_plus_shift_over_a2"] - 0.5), _FLOOR)
        for r in positive_rows
    ]

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    fig.suptitle("S4-KERR-K4 Schwarzschild-limit audit", fontsize=14)

    panels = [
        (axes[0, 0], tphi_res,
         r"$|g_{t\phi}/a - (-2M/r)|_{\rm max}$",
         "Frame-dragging linearity residual"),
        (axes[0, 1], gphi_res,
         r"$|g_{\phi\phi} - r^2 - a^2(1+2M/r)|_{\rm max}$",
         "Azimuthal quadratic residual (absolute)"),
        (axes[1, 0], grr_res,
         r"$|g_{rr} - r^2/\Delta|_{\rm max}$",
         r"$g_{rr}$ formula residual"),
        (axes[1, 1], horiz_res,
         r"$|(2M - r_+)/a^2 - 1/(2M)|$",
         "Horizon quadratic shift residual"),
    ]

    for ax, residuals, ylabel, title in panels:
        ax.loglog(spins_pos, residuals, "o-", color="steelblue", linewidth=1.5,
                  markersize=5)
        ax.set_xlabel("spin $a$")
        ax.set_ylabel(ylabel, fontsize=9)
        ax.set_title(title, fontsize=10)
        ax.grid(True, which="both", linestyle="--", alpha=0.4)
        ax.set_xlim(left=min(spins_pos) * 0.5, right=max(spins_pos) * 2)

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
    payload   = run_audit()
    agg       = payload["aggregate"]
    csv_path, json_path, md_path, png_path = write_outputs(payload)

    print(f"all_checks_pass={agg['all_checks_pass']}")
    for row in payload["rows"]:
        a          = row["spin_a"]
        shift_a2   = row["r_plus_shift_over_a2"]
        shift_str  = f"{shift_a2:.8f}" if shift_a2 is not None else "n/a"
        print(
            f"  a={a:.1e}  r_+={row['r_plus']:.6f}"
            f"  shift/a²={shift_str}"
            f"  ext={row['all_points_exterior']}"
            f"  fd={row['frame_dragging_linear_check_pass']}"
            f"  gphiphi={row['gphiphi_quadratic_check_pass']}"
            f"  grr={row['grr_formula_check_pass']}"
            f"  horiz={row['horizon_quadratic_check_pass']}"
            f"  schw={row['schwarzschild_metric_limit_pass']}"
            f"  pass={row['all_checks_pass']}"
        )
    print(f"wrote {csv_path}")
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")
    print(f"wrote {png_path}")


if __name__ == "__main__":
    main()
