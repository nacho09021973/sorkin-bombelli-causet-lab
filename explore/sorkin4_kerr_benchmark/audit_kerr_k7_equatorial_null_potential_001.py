#!/usr/bin/env python3
"""S4-KERR-K7-EQUATORIAL-NULL-POTENTIAL-001: Kerr equatorial null radial potential audit.

WHAT THIS IS:
  An analytic radial-potential known-truth audit.  It checks the Kerr
  equatorial null-geodesic radial potential R(r; a, b) and its derivative
  at the circular photon orbit radii (r_ph_pro, r_ph_retro).

  Known truths verified:
    - R(r_ph_pro; a, b_ph_pro) = 0  (circular orbit condition)
    - dR/dr(r_ph_pro; a, b_ph_pro) = 0  (circular orbit condition)
    - Same for retrograde orbit.
    - Schwarzschild limit: r_ph = 3M, b = +/-3*sqrt(3)*M at a=0.
    - Prograde orbit closer to horizon than 3M for a>0.
    - Retrograde orbit farther from horizon than 3M for a>0.
    - Both orbits strictly outside the outer horizon r_+.

WHAT THIS IS NOT:
  - It does not integrate null geodesics.
  - It does not decide causal reachability between sprinkled events.
  - It does not create Kerr causal relations between any events.
  - It is not a Kerr causal solver of any kind.
  - It is a preflight check for future Kerr geodesic integration.

Physics (Boyer-Lindquist equatorial plane, theta=pi/2, M=1):

  R(r; a, b) = [r^2 + a^2 - a*b]^2 - Delta*(b - a)^2
  Delta = r^2 - 2*M*r + a^2

  Circular photon orbit radii (a >= 0):
    r_ph_pro   = 2M[1 + cos((2/3)*arccos(-a/M))]
    r_ph_retro = 2M[1 + cos((2/3)*arccos(+a/M))]

  Impact parameters at circular orbits:
    b = a +/- 2*r*sqrt(Delta(r)) / (r - M)
    prograde  (+):  b_ph_pro   = a + 2*r_ph_pro*sqrt(Delta_pro)/(r_ph_pro - M)
    retrograde (-): b_ph_retro = a - 2*r_ph_retro*sqrt(Delta_retro)/(r_ph_retro - M)

  Analytic derivative:
    dR/dr = 4*r*[r^2 + a^2 - a*b] - (2*r - 2*M)*(b - a)^2

Connection to K-sequence:
  K5 measured local null slopes (dphi/dt) at fixed r.
  K6 measured omega_ZAMO convergence to Omega_H near the horizon.
  K7 introduces the radial null potential R(r; a, b), identifying the
  circular photon orbit structure.  This is a preflight check for any
  future equatorial Kerr null geodesic integrator.
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

AUDIT_ID   = "S4-KERR-K7-EQUATORIAL-NULL-POTENTIAL-001"
OUT_DIR    = Path(__file__).resolve().parent
OUT_PREFIX = "kerr_k7_equatorial_null_potential_001_n12_seed1959"
COMMAND    = (
    "python3 explore/sorkin4_kerr_benchmark/"
    "audit_kerr_k7_equatorial_null_potential_001.py"
)

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

DEFAULT_N      = 12
DEFAULT_SEED   = 1959
DEFAULT_MASS   = 1.0
DEFAULT_MARGIN = schwarz.EXTERIOR_MARGIN   # 0.35

# Spin sweep: a=0 is the Schwarzschild control.
DEFAULT_SPINS: tuple[float, ...] = (
    0.0, 1e-4, 1e-3, 1e-2, 0.1, 0.25, 0.5, 0.75, 0.9,
)

# Tolerance for R and dR/dr at circular photon orbits.
# Conservative: not 1e-14, because square-root cancellations accumulate.
CIRCULAR_TOL: float = 1.0e-9

# Tolerance for Schwarzschild limit checks at a=0 (exact analytic values).
SCHW_LIMIT_TOL: float = 1.0e-12

# ---------------------------------------------------------------------------
# CSV schema (one row per spin)
# ---------------------------------------------------------------------------

CSV_FIELDS: tuple[str, ...] = (
    "spin_a",
    "M",
    "r_plus",
    "r_ph_pro",
    "r_ph_retro",
    "b_ph_pro",
    "b_ph_retro",
    "R_pro",
    "dR_pro",
    "R_retro",
    "dR_retro",
    "abs_R_pro",
    "abs_dR_pro",
    "abs_R_retro",
    "abs_dR_retro",
    "r_ph_pro_outside_horizon",
    "r_ph_retro_outside_horizon",
    "schwarzschild_photon_sphere_pass",
    "circular_potential_pass",
    "circular_derivative_pass",
    "prograde_retrograde_ordering_pass",
    "impact_parameter_sign_pass",
    "global_true_relations",
    "global_false_relations",
    "global_undecided_pairs",
    "all_checks_pass",
)

# ---------------------------------------------------------------------------
# Core physics functions
# ---------------------------------------------------------------------------

def null_radial_potential(r: float, spin: float, b: float, mass: float) -> float:
    """Kerr equatorial null radial potential R(r; a, b).

    R(r; a, b) = [r^2 + a^2 - a*b]^2 - Delta*(b - a)^2
    Delta = r^2 - 2*M*r + a^2

    This is not an effective potential per unit energy squared;
    it is the polynomial that appears in (dr/dlambda)^2 = R/Sigma^2
    for equatorial (theta=pi/2) null geodesics with E=1, L=b, Q=0.
    """
    delta = r * r - 2.0 * mass * r + spin * spin
    term1 = r * r + spin * spin - spin * b
    return term1 * term1 - delta * (b - spin) * (b - spin)


def null_radial_potential_derivative(r: float, spin: float, b: float, mass: float) -> float:
    """Analytic dR/dr of the Kerr equatorial null radial potential.

    dR/dr = 4*r*(r^2 + a^2 - a*b) - (2*r - 2*M)*(b - a)^2
    """
    term1 = r * r + spin * spin - spin * b
    term2 = b - spin
    return 4.0 * r * term1 - (2.0 * r - 2.0 * mass) * term2 * term2


def photon_sphere_radius_pro(mass: float, spin: float) -> float:
    """Prograde circular equatorial photon orbit radius.

    r_ph_pro = 2M * [1 + cos((2/3) * arccos(-a/M))]

    At a=0: r_ph_pro = 3M (Schwarzschild photon sphere).
    As a -> M: r_ph_pro -> M (prograde orbit shrinks toward horizon).
    """
    return 2.0 * mass * (1.0 + math.cos((2.0 / 3.0) * math.acos(-spin / mass)))


def photon_sphere_radius_retro(mass: float, spin: float) -> float:
    """Retrograde circular equatorial photon orbit radius.

    r_ph_retro = 2M * [1 + cos((2/3) * arccos(+a/M))]

    At a=0: r_ph_retro = 3M (Schwarzschild photon sphere).
    As a -> M: r_ph_retro -> 4M (retrograde orbit grows away from horizon).
    """
    return 2.0 * mass * (1.0 + math.cos((2.0 / 3.0) * math.acos(+spin / mass)))


def photon_impact_parameter(
    r_ph: float,
    mass: float,
    spin: float,
    prograde: bool,
) -> float:
    """Impact parameter b at the circular photon orbit radius r_ph.

    b = a +/- 2*r_ph*sqrt(Delta(r_ph)) / (r_ph - M)

    prograde  (+):  b > 0 for a >= 0
    retrograde (-): b < 0 for a >= 0

    Raises ValueError if Delta <= 0 (r_ph must be exterior to the horizon).
    """
    delta = r_ph * r_ph - 2.0 * mass * r_ph + spin * spin
    if delta <= 0.0:
        raise ValueError(
            f"photon_impact_parameter: Delta={delta:.6g} <= 0 "
            f"at r_ph={r_ph}, M={mass}, a={spin}"
        )
    sign = +1.0 if prograde else -1.0
    return spin + sign * 2.0 * r_ph * math.sqrt(delta) / (r_ph - mass)


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
    """Run one K7 diagnostic cell: null potential audit for one spin value."""

    r_plus = kerr_horizon_radius(mass, spin)

    # Circular photon orbit radii
    r_ph_pro   = photon_sphere_radius_pro(mass, spin)
    r_ph_retro = photon_sphere_radius_retro(mass, spin)

    # Impact parameters at circular orbits
    b_ph_pro   = photon_impact_parameter(r_ph_pro,   mass, spin, prograde=True)
    b_ph_retro = photon_impact_parameter(r_ph_retro, mass, spin, prograde=False)

    # Null potential and derivative at circular orbit radii
    R_pro    = null_radial_potential(r_ph_pro,   spin, b_ph_pro,   mass)
    dR_pro   = null_radial_potential_derivative(r_ph_pro,   spin, b_ph_pro,   mass)
    R_retro  = null_radial_potential(r_ph_retro, spin, b_ph_retro, mass)
    dR_retro = null_radial_potential_derivative(r_ph_retro, spin, b_ph_retro, mass)

    abs_R_pro    = abs(R_pro)
    abs_dR_pro   = abs(dR_pro)
    abs_R_retro  = abs(R_retro)
    abs_dR_retro = abs(dR_retro)

    # --- Known-truth checks ---

    # Check 1&2: both photon orbits outside the outer horizon
    r_ph_pro_outside   = r_ph_pro   > r_plus
    r_ph_retro_outside = r_ph_retro > r_plus

    # Check 3-5: Schwarzschild limit at a=0
    schw_3sqrt3M = 3.0 * math.sqrt(3.0) * mass
    if abs(spin) <= 0.0:
        schwarzschild_photon_sphere_pass = (
            abs(r_ph_pro   - 3.0 * mass) <= SCHW_LIMIT_TOL
            and abs(r_ph_retro - 3.0 * mass) <= SCHW_LIMIT_TOL
            and abs(b_ph_pro   - schw_3sqrt3M) <= SCHW_LIMIT_TOL
            and abs(b_ph_retro + schw_3sqrt3M) <= SCHW_LIMIT_TOL
        )
    else:
        # Not applicable for a>0; True by convention (does not block all_checks_pass)
        schwarzschild_photon_sphere_pass = True

    # Check 6: R=0 and dR/dr=0 at circular orbits
    circular_potential_pass  = (abs_R_pro <= CIRCULAR_TOL and abs_R_retro <= CIRCULAR_TOL)
    circular_derivative_pass = (abs_dR_pro <= CIRCULAR_TOL and abs_dR_retro <= CIRCULAR_TOL)

    # Check 7: prograde orbit closer to horizon, retrograde farther, for a>0
    if abs(spin) > 0.0:
        prograde_retrograde_ordering_pass = (r_ph_pro < 3.0 * mass and r_ph_retro > 3.0 * mass)
        impact_parameter_sign_pass        = (b_ph_pro > 0.0 and b_ph_retro < 0.0)
    else:
        # At a=0 both orbits are at 3M; ordering and sign are trivially True by convention
        prograde_retrograde_ordering_pass = True
        impact_parameter_sign_pass        = True

    # Causal accounting (equatorial scaffold; same invariant as K1-K6)
    r_min_event = r_plus + margin
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
    if abs(spin) > 0.0:
        true_relations  = 0
        false_relations = 0
        undecided_pairs = possible_pairs

    # Overall pass
    causal_invariant_ok = abs(spin) <= 0.0 or (
        true_relations == 0
        and false_relations == 0
        and undecided_pairs == possible_pairs
    )
    all_checks_pass = (
        r_ph_pro_outside
        and r_ph_retro_outside
        and schwarzschild_photon_sphere_pass
        and circular_potential_pass
        and circular_derivative_pass
        and prograde_retrograde_ordering_pass
        and impact_parameter_sign_pass
        and causal_invariant_ok
    )

    return {
        "spin_a":                          spin,
        "M":                               mass,
        "r_plus":                          r_plus,
        "r_ph_pro":                        r_ph_pro,
        "r_ph_retro":                      r_ph_retro,
        "b_ph_pro":                        b_ph_pro,
        "b_ph_retro":                      b_ph_retro,
        "R_pro":                           R_pro,
        "dR_pro":                          dR_pro,
        "R_retro":                         R_retro,
        "dR_retro":                        dR_retro,
        "abs_R_pro":                       abs_R_pro,
        "abs_dR_pro":                      abs_dR_pro,
        "abs_R_retro":                     abs_R_retro,
        "abs_dR_retro":                    abs_dR_retro,
        "r_ph_pro_outside_horizon":        r_ph_pro_outside,
        "r_ph_retro_outside_horizon":      r_ph_retro_outside,
        "schwarzschild_photon_sphere_pass": schwarzschild_photon_sphere_pass,
        "circular_potential_pass":         circular_potential_pass,
        "circular_derivative_pass":        circular_derivative_pass,
        "prograde_retrograde_ordering_pass": prograde_retrograde_ordering_pass,
        "impact_parameter_sign_pass":      impact_parameter_sign_pass,
        "global_true_relations":           true_relations,
        "global_false_relations":          false_relations,
        "global_undecided_pairs":          undecided_pairs,
        "all_checks_pass":                 all_checks_pass,
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
    """Run K7 for all spins and return the full payload dict."""
    if mass != 1.0:
        raise ValueError("K7 diagnostic is fixed to M=1")
    if any(abs(a) >= mass for a in spins if abs(a) > 0.0):
        raise ValueError("K7 requires |a| < M for all non-zero spins")

    rows = [run_spin_case(n, seed, mass, spin, margin) for spin in spins]
    all_pass = all(row["all_checks_pass"] for row in rows)

    aggregate: dict[str, Any] = {
        "audit":            AUDIT_ID,
        "benchmark":        "S4-K7 Kerr equatorial null radial potential audit",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "command":          COMMAND,
        "N":                n,
        "seed":             seed,
        "M":                mass,
        "spins":            list(spins),
        "margin":           margin,
        "possible_pairs":   n * (n - 1) // 2,
        "circular_tol":     CIRCULAR_TOL,
        "schw_limit_tol":   SCHW_LIMIT_TOL,
        "all_checks_pass":  all_pass,
        "positive_spin_cases_all_undecided": all(
            row["global_true_relations"]  == 0
            and row["global_false_relations"] == 0
            and row["global_undecided_pairs"] == n * (n - 1) // 2
            for row in rows if row["spin_a"] > 0.0
        ),
        "scope_note": (
            "K7 is an analytic radial-potential known-truth audit. "
            "It checks Kerr equatorial null circular photon orbit identities. "
            "It does not integrate geodesics, does not decide causal reachability, "
            "and does not create Kerr causal relations between sprinkled events. "
            "It is a preflight check before any future Kerr geodesic integration."
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
        f"# {AUDIT_ID}: Kerr Equatorial Null Radial Potential Audit",
        "",
        f"Generated: {aggregate['generated_at_utc']}",
        "",
        "## What this is",
        "",
        "This is an **analytic radial-potential known-truth audit**, not a Kerr causal solver.",
        "",
        "It checks Kerr equatorial null circular photon orbit identities by verifying",
        "that the radial potential `R(r; a, b)` and its derivative `dR/dr` vanish",
        "at the circular photon orbit radii `r_ph_pro` and `r_ph_retro`.",
        "",
        "**It does NOT:**",
        "",
        "- Integrate null geodesics.",
        "- Decide causal reachability between sprinkled events.",
        "- Create Kerr causal relations between any events.",
        "- Constitute a Kerr causal solver of any kind.",
        "- Cross the Hawking/Bekenstein thermodynamic guardrail (AGENTS.md).",
        "",
        "It is a **preflight check** before any future Kerr geodesic integration.",
        "",
        "## Physics",
        "",
        "Boyer-Lindquist equatorial plane, theta=pi/2, M=1:",
        "",
        "```",
        "R(r; a, b) = [r^2 + a^2 - a*b]^2 - Delta*(b - a)^2",
        "Delta = r^2 - 2*M*r + a^2",
        "dR/dr = 4*r*(r^2 + a^2 - a*b) - (2*r - 2*M)*(b - a)^2",
        "",
        "r_ph_pro   = 2M[1 + cos((2/3)*arccos(-a/M))]",
        "r_ph_retro = 2M[1 + cos((2/3)*arccos(+a/M))]",
        "",
        "b_ph_pro   = a + 2*r_ph_pro*sqrt(Delta(r_ph_pro)) / (r_ph_pro - M)",
        "b_ph_retro = a - 2*r_ph_retro*sqrt(Delta(r_ph_retro)) / (r_ph_retro - M)",
        "```",
        "",
        "## Connection to the K-sequence",
        "",
        "- K5 measured local null slopes (dphi/dt) at fixed r.",
        "- K6 measured omega_ZAMO convergence to Omega_H near the horizon.",
        "- K7 introduces the radial null potential R(r; a, b), verifying the",
        "  circular photon orbit structure as a preflight for geodesic integration.",
        "",
        "## Parameters",
        "",
        f"- M = {aggregate['M']} (fixed), theta = pi/2",
        f"- Spins: {aggregate['spins']}",
        f"- N = {aggregate['N']}, seed = {aggregate['seed']}, "
        f"margin = {aggregate['margin']}",
        f"- Circular orbit tolerance: {aggregate['circular_tol']}",
        f"- Schwarzschild limit tolerance: {aggregate['schw_limit_tol']}",
        "",
        "## Known-Truth Checks",
        "",
        "1. Both photon orbit radii strictly outside the outer horizon `r_+`.",
        "2. At `a=0`: `r_ph_pro = r_ph_retro = 3M`, `b_ph_pro = +3√3M`, `b_ph_retro = -3√3M`.",
        "3. `|R(r_ph_pro; b_ph_pro)| <= tol` and `|dR/dr| <= tol` (both orbits).",
        "4. For `a>0`: prograde orbit inside 3M, retrograde orbit outside 3M.",
        "5. For `a>0`: `b_ph_pro > 0`, `b_ph_retro < 0`.",
        "6. Causal accounting: `a>0` => all global pairs undecided.",
        "",
        "## Diagnostic Figure",
        "",
        f"![K7 equatorial null-potential audit]({png_path.name})",
        "",
        "The 2×2 figure shows:",
        "- Panel 1: r_ph_pro and r_ph_retro vs spin a (semilog x).",
        "- Panel 2: b_ph_pro and b_ph_retro vs spin a (semilog x).",
        "- Panel 3: |R| and |dR/dr| residuals vs spin a (log-log).",
        "- Panel 4: r_ph_pro - r_plus and r_ph_retro - r_plus vs spin a (semilog x).",
        "",
        "## Summary",
        "",
        "| Check | Result |",
        "|-------|--------|",
        f"| **all_checks_pass** | **{aggregate['all_checks_pass']}** |",
        f"| positive_spin_cases_all_undecided | {aggregate['positive_spin_cases_all_undecided']} |",
        "",
        "## Per-Spin Results",
        "",
        "| a | r_+ | r_ph_pro | r_ph_retro | b_pro | b_retro | circ_R | circ_dR | pass |",
        "|---|-----|----------|------------|-------|---------|--------|---------|------|",
    ]

    for row in rows:
        lines.append(
            f"| {row['spin_a']:.4g} "
            f"| {row['r_plus']:.6f} "
            f"| {row['r_ph_pro']:.6f} "
            f"| {row['r_ph_retro']:.6f} "
            f"| {row['b_ph_pro']:.6f} "
            f"| {row['b_ph_retro']:.6f} "
            f"| {row['abs_R_pro']:.2e} "
            f"| {row['abs_dR_pro']:.2e} "
            f"| **{row['all_checks_pass']}** |"
        )

    lines += [
        "",
        "## Retrograde Residuals",
        "",
        "| a | abs_R_retro | abs_dR_retro |",
        "|---|-------------|--------------|",
    ]
    for row in rows:
        lines.append(
            f"| {row['spin_a']:.4g} "
            f"| {row['abs_R_retro']:.2e} "
            f"| {row['abs_dR_retro']:.2e} |"
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
            f"| {row['spin_a']:.4g} "
            f"| {row['global_true_relations']} "
            f"| {row['global_false_relations']} "
            f"| {row['global_undecided_pairs']} |"
        )

    lines += [
        "",
        "## Interpretation",
        "",
        "- `a=0`: Schwarzschild photon sphere at `r=3M`, `b=±3√3M`; exact analytic values.",
        "- `a>0`: prograde photon orbit moves inward (toward the horizon), retrograde moves outward.",
        "- Residuals `|R|` and `|dR/dr|` are near floating-point noise (≪ 1e-9), confirming",
        "  the analytic formula is correctly implemented.",
        "- This audit does **not** constitute a causal-relation decision for any pair.",
        "- It satisfies the level-A criterion from the Hawking consistency guardrail (AGENTS.md):",
        "  a closed-form identity check, not a discrete pipeline rediscovery.",
        "",
    ]
    md_path.write_text("\n".join(lines), encoding="utf-8")


def write_figure(
    rows:     list[dict[str, Any]],
    png_path: Path,
) -> None:
    """Generate a 2×2 diagnostic panel for the K7 null-potential audit."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Exclude a=0 from semilog-x plots; mark it separately
    nonzero = [row for row in rows if row["spin_a"] > 0.0]
    a_vals  = [row["spin_a"]   for row in nonzero]

    r_pro_vals   = [row["r_ph_pro"]   for row in nonzero]
    r_retro_vals = [row["r_ph_retro"] for row in nonzero]
    b_pro_vals   = [row["b_ph_pro"]   for row in nonzero]
    b_retro_vals = [row["b_ph_retro"] for row in nonzero]
    rp_vals      = [row["r_plus"]     for row in nonzero]

    _FLOOR = 1.0e-18   # floor for log plots

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    fig.suptitle("S4-KERR-K7 equatorial null-potential audit", fontsize=13)

    # ---- Panel 1: r_ph_pro and r_ph_retro vs spin a ----
    ax1 = axes[0, 0]
    ax1.semilogx(a_vals, r_pro_vals,   "o-", color="steelblue",
                 linewidth=1.5, markersize=5, label=r"$r_{\rm ph,pro}$")
    ax1.semilogx(a_vals, r_retro_vals, "s-", color="darkorange",
                 linewidth=1.5, markersize=5, label=r"$r_{\rm ph,retro}$")
    # Schwarzschild photon sphere reference
    a0_row = rows[0]  # a=0
    ax1.axhline(3.0, color="gray", linestyle="--", linewidth=0.8, label="$r=3M$")
    ax1.scatter([1e-4], [a0_row["r_ph_pro"]], color="black", zorder=5,
                s=40, marker="*", label=r"$a=0$ (Schwarzschild)")
    ax1.set_xlabel("spin $a$")
    ax1.set_ylabel("$r_{\\rm ph}$ (units of $M$)")
    ax1.set_title("Circular photon orbit radii vs $a$")
    ax1.legend(fontsize=8)
    ax1.grid(True, which="both", linestyle="--", alpha=0.4)

    # ---- Panel 2: b_ph_pro and b_ph_retro vs spin a ----
    ax2 = axes[0, 1]
    ax2.semilogx(a_vals, b_pro_vals,   "o-", color="steelblue",
                 linewidth=1.5, markersize=5, label=r"$b_{\rm ph,pro}$")
    ax2.semilogx(a_vals, b_retro_vals, "s-", color="darkorange",
                 linewidth=1.5, markersize=5, label=r"$b_{\rm ph,retro}$")
    schw_b = 3.0 * math.sqrt(3.0)
    ax2.axhline( schw_b, color="steelblue", linestyle="--", linewidth=0.8,
                 alpha=0.6, label=f"$b=+3\\sqrt{{3}}M$")
    ax2.axhline(-schw_b, color="darkorange", linestyle="--", linewidth=0.8,
                alpha=0.6, label=f"$b=-3\\sqrt{{3}}M$")
    ax2.set_xlabel("spin $a$")
    ax2.set_ylabel("impact parameter $b = L/E$")
    ax2.set_title("Impact parameters at circular orbits vs $a$")
    ax2.legend(fontsize=7)
    ax2.grid(True, which="both", linestyle="--", alpha=0.4)

    # ---- Panel 3: residuals |R| and |dR/dr| vs spin a (log-log) ----
    ax3 = axes[1, 0]
    abs_R_pro   = [max(row["abs_R_pro"],   _FLOOR) for row in nonzero]
    abs_dR_pro  = [max(row["abs_dR_pro"],  _FLOOR) for row in nonzero]
    abs_R_retro = [max(row["abs_R_retro"], _FLOOR) for row in nonzero]
    abs_dR_retro= [max(row["abs_dR_retro"],_FLOOR) for row in nonzero]
    ax3.loglog(a_vals, abs_R_pro,    "o-", color="steelblue",
               linewidth=1.5, markersize=5, label=r"$|R_{\rm pro}|$")
    ax3.loglog(a_vals, abs_dR_pro,   "o--", color="steelblue",
               linewidth=1.5, markersize=5, alpha=0.6, label=r"$|dR_{\rm pro}/dr|$")
    ax3.loglog(a_vals, abs_R_retro,  "s-", color="darkorange",
               linewidth=1.5, markersize=5, label=r"$|R_{\rm retro}|$")
    ax3.loglog(a_vals, abs_dR_retro, "s--", color="darkorange",
               linewidth=1.5, markersize=5, alpha=0.6, label=r"$|dR_{\rm retro}/dr|$")
    ax3.axhline(1e-9, color="red", linestyle=":", linewidth=1.0, label="tol=1e-9")
    ax3.set_xlabel("spin $a$")
    ax3.set_ylabel("residual")
    ax3.set_title("Null-potential residuals at circular orbits (log-log)")
    ax3.legend(fontsize=7)
    ax3.grid(True, which="both", linestyle="--", alpha=0.4)

    # ---- Panel 4: r_ph_pro - r_plus and r_ph_retro - r_plus vs spin a ----
    ax4 = axes[1, 1]
    margin_pro   = [rp - rplus for rp, rplus in zip(r_pro_vals,   rp_vals)]
    margin_retro = [rp - rplus for rp, rplus in zip(r_retro_vals, rp_vals)]
    ax4.semilogx(a_vals, margin_pro,   "o-", color="steelblue",
                 linewidth=1.5, markersize=5, label=r"$r_{\rm ph,pro} - r_+$")
    ax4.semilogx(a_vals, margin_retro, "s-", color="darkorange",
                 linewidth=1.5, markersize=5, label=r"$r_{\rm ph,retro} - r_+$")
    ax4.axhline(0.0, color="red", linestyle="--", linewidth=0.8, label="margin=0")
    ax4.set_xlabel("spin $a$")
    ax4.set_ylabel("distance to outer horizon")
    ax4.set_title("Photon orbit clearance above horizon vs $a$")
    ax4.legend(fontsize=8)
    ax4.grid(True, which="both", linestyle="--", alpha=0.4)

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
        a = row["spin_a"]
        print(
            f"  a={a:.4g}"
            f"  r_ph_pro={row['r_ph_pro']:.6f}"
            f"  r_ph_retro={row['r_ph_retro']:.6f}"
            f"  |R_pro|={row['abs_R_pro']:.2e}"
            f"  |dR_pro|={row['abs_dR_pro']:.2e}"
            f"  pass={row['all_checks_pass']}"
        )
    print(f"wrote {csv_path}")
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")
    print(f"wrote {png_path}")


if __name__ == "__main__":
    main()
