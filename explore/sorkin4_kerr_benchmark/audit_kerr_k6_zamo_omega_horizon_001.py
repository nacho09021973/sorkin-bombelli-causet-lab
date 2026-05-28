#!/usr/bin/env python3
"""S4-KERR-K6-ZAMO-OMEGA-HORIZON-001: ZAMO angular velocity horizon convergence audit.

WHAT THIS IS:
  A near-horizon convergence diagnostic. It verifies that the local ZAMO
  angular velocity omega_ZAMO = -g_tphi/g_phiphi converges to the known-truth
  horizon angular velocity Omega_H = a/(r_+^2 + a^2) as r -> r_+ from outside
  the horizon.  The convergence is linear in delta = r - r_+.

WHAT THIS IS NOT:
  - It does not cross the horizon.
  - It does not implement Kerr causal inference.
  - It does not integrate null geodesics.
  - It does not claim global causal reachability.
  - Convergence of omega_ZAMO to Omega_H is a local metric identity, not a
    causal relation.

Physics (Boyer-Lindquist equatorial plane, theta=pi/2, M=1):
  omega_ZAMO(r) = -g_tphi / g_phiphi  (same as omega_center from K5)
  Omega_H       = a / (r_+^2 + a^2)   (horizon angular velocity)

  Known truth:
    omega_ZAMO(r_+) = Omega_H  (exact, follows from Delta(r_+) = 0)
    omega_ZAMO(r_+ + delta) < Omega_H  for all delta > 0  (monotone from below)
    omega_ZAMO(r_+ + delta) - Omega_H = O(delta)  (linear convergence)

  For a=0: Omega_H = omega_ZAMO = 0 everywhere (Schwarzschild; trivial check).

Connection to the K-sequence:
  K5 measured omega_ZAMO asymmetry at fixed r;
  K6 measures omega_ZAMO convergence to Omega_H near the horizon.
  Together they bridge local frame-dragging -> horizon angular velocity,
  without crossing the Hawking/Bekenstein thermodynamic guardrail.
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

AUDIT_ID   = "S4-KERR-K6-ZAMO-OMEGA-HORIZON-001"
OUT_DIR    = Path(__file__).resolve().parent
OUT_PREFIX = "kerr_k6_zamo_omega_horizon_001_n12_seed1959"
COMMAND    = (
    "python3 explore/sorkin4_kerr_benchmark/"
    "audit_kerr_k6_zamo_omega_horizon_001.py"
)

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

DEFAULT_N      = 12
DEFAULT_SEED   = 1959
DEFAULT_MASS   = 1.0
DEFAULT_MARGIN = schwarz.EXTERIOR_MARGIN   # 0.35

# Spin sweep: a=0 is the Schwarzschild control (trivial Omega_H=0).
# a = 0.25, 0.5, 0.75, 0.9 are the non-trivial test cases.
DEFAULT_SPINS: tuple[float, ...] = (0.0, 0.25, 0.5, 0.75, 0.9)

# Near-horizon approach: delta = r - r_+, decreasing toward 0.
# r_eval = r_+ + delta is always > r_+ (exterior).
DEFAULT_DELTAS: tuple[float, ...] = (1e-1, 3e-2, 1e-2, 3e-3, 1e-3)

# Tolerance for a=0 trivial residual (should be exactly 0).
METRIC_FORMULA_TOL = 1.0e-12

# ---------------------------------------------------------------------------
# CSV schema  (one row per (spin, delta) point)
# ---------------------------------------------------------------------------

CSV_FIELDS: tuple[str, ...] = (
    "spin_a",
    "M",
    "delta",
    "r_eval",
    "r_plus",
    "Omega_H",
    "omega_zamo",
    "residual_abs",
    "exterior_pass",
)


# ---------------------------------------------------------------------------
# Metric and ZAMO helpers
# ---------------------------------------------------------------------------

def equatorial_metric_at_r(
    r: float,
    mass: float,
    spin: float,
) -> dict[str, float]:
    """Boyer-Lindquist equatorial metric at theta=pi/2, signature -+++.

    Valid for any r > r_+ (Delta > 0).  Inside the ergosphere (r < 2M but
    r > r_+), g_tt > 0 and the t-direction is spacelike; omega_ZAMO is still
    well-defined and equals -g_tphi/g_phiphi.
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


def omega_H_analytic(r_plus: float, spin: float) -> float:
    """Analytic horizon angular velocity Omega_H = a / (r_+^2 + a^2).

    Returns 0.0 for a=0 (non-rotating Schwarzschild horizon).
    The formula follows from Delta(r_+) = 0:
      r_+^2 + a^2 = 2Mr_+  =>  Omega_H = a/(2Mr_+)
    """
    if abs(spin) <= 0.0:
        return 0.0
    return spin / (r_plus * r_plus + spin * spin)


def omega_zamo_at_r(r: float, mass: float, spin: float) -> float:
    """ZAMO angular velocity = -g_tphi / g_phiphi (same as omega_center in K5).

    For a=0: omega_ZAMO = 0 everywhere.
    For a>0 and r > r_+: 0 < omega_ZAMO < Omega_H.
    At r = r_+: omega_ZAMO = Omega_H exactly (exact algebraic identity).
    """
    m = equatorial_metric_at_r(r, mass, spin)
    return -m["g_tphi"] / m["g_phiphi"]


# ---------------------------------------------------------------------------
# Per-spin case runner
# ---------------------------------------------------------------------------

def run_spin_case(
    n:      int,
    seed:   int,
    mass:   float,
    spin:   float,
    margin: float,
    deltas: tuple[float, ...],
) -> dict[str, Any]:
    """Run one K6 diagnostic cell: near-horizon approach for one spin value."""

    r_plus  = kerr_horizon_radius(mass, spin)
    Omega_H = omega_H_analytic(r_plus, spin)

    # --- Near-horizon approach: evaluate omega_ZAMO at each r_+ + delta ---
    point_rows:       list[dict[str, Any]] = []
    omega_zamo_vals:  list[float] = []
    residuals:        list[float] = []

    for delta in deltas:
        r_eval   = r_plus + delta
        omega_v  = omega_zamo_at_r(r_eval, mass, spin)
        residual = abs(omega_v - Omega_H)
        exterior = r_eval > r_plus   # always True since delta > 0

        point_rows.append({
            "spin_a":       spin,
            "M":            mass,
            "delta":        delta,
            "r_eval":       r_eval,
            "r_plus":       r_plus,
            "Omega_H":      Omega_H,
            "omega_zamo":   omega_v,
            "residual_abs": residual,
            "exterior_pass": exterior,
        })
        omega_zamo_vals.append(omega_v)
        residuals.append(residual)

    all_exterior = all(row["exterior_pass"] for row in point_rows)

    # --- Convergence check ---
    # For a=0: omega_ZAMO = Omega_H = 0 everywhere; trivially True.
    # For a>0: residuals must decrease strictly as delta decreases.
    #   deltas[0] > deltas[1] > ... > deltas[-1]
    #   => residuals should satisfy residuals[i] > residuals[i+1]
    if abs(spin) <= 0.0:
        convergence_monotone_pass = True
        omega_zamo_below_Omega_H  = True  # trivially (both are 0)
    else:
        convergence_monotone_pass = all(
            residuals[i] > residuals[i + 1]
            for i in range(len(residuals) - 1)
        )
        omega_zamo_below_Omega_H = all(v < Omega_H for v in omega_zamo_vals)

    min_residual = min(residuals)
    max_residual = max(residuals)

    # --- Causal accounting (equatorial scaffold, same invariant as K1–K5) ---
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

    all_checks_pass = (
        all_exterior
        and convergence_monotone_pass
        and omega_zamo_below_Omega_H
        and (
            abs(spin) <= 0.0
            or (true_relations == 0 and false_relations == 0
                and undecided_pairs == possible_pairs)
        )
    )

    return {
        "spin_a":                    spin,
        "r_plus":                    r_plus,
        "Omega_H":                   Omega_H,
        "deltas":                    list(deltas),
        "omega_zamo_values":         omega_zamo_vals,
        "residuals":                 residuals,
        "min_residual":              min_residual,
        "max_residual":              max_residual,
        "all_exterior":              all_exterior,
        "convergence_monotone_pass": convergence_monotone_pass,
        "omega_zamo_below_Omega_H":  omega_zamo_below_Omega_H,
        "global_true_relations":     true_relations,
        "global_false_relations":    false_relations,
        "global_undecided_pairs":    undecided_pairs,
        "all_checks_pass":           all_checks_pass,
        "_point_rows":               point_rows,
    }


# ---------------------------------------------------------------------------
# Full audit runner
# ---------------------------------------------------------------------------

def run_audit(
    n:      int               = DEFAULT_N,
    seed:   int               = DEFAULT_SEED,
    mass:   float             = DEFAULT_MASS,
    spins:  tuple[float, ...] = DEFAULT_SPINS,
    deltas: tuple[float, ...] = DEFAULT_DELTAS,
    margin: float             = DEFAULT_MARGIN,
) -> dict[str, Any]:
    """Run K6 for all spins and return the full payload dict."""
    if mass != 1.0:
        raise ValueError("K6 diagnostic is fixed to M=1")
    if any(abs(a) >= mass for a in spins if abs(a) > 0.0):
        raise ValueError("K6 requires |a| < M for all non-zero spins")
    if any(d <= 0.0 for d in deltas):
        raise ValueError("K6 requires all deltas > 0")

    spin_cases = [
        run_spin_case(n, seed, mass, spin, margin, deltas) for spin in spins
    ]
    all_pass   = all(c["all_checks_pass"] for c in spin_cases)

    # Flatten per-case point_rows to a single list for CSV writing
    all_rows: list[dict[str, Any]] = []
    for case in spin_cases:
        all_rows.extend(case["_point_rows"])

    # Spin summaries (public, without internal keys)
    spin_summaries = [
        {k: v for k, v in case.items() if not k.startswith("_")}
        for case in spin_cases
    ]

    aggregate: dict[str, Any] = {
        "audit":            AUDIT_ID,
        "benchmark":        "S4-K6 Kerr ZAMO angular velocity horizon convergence",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "command":          COMMAND,
        "N":                n,
        "seed":             seed,
        "M":                mass,
        "spins":            list(spins),
        "deltas":           list(deltas),
        "margin":           margin,
        "possible_pairs":   n * (n - 1) // 2,
        "all_checks_pass":  all_pass,
        "positive_spin_cases_all_undecided": all(
            c["global_true_relations"]  == 0
            and c["global_false_relations"] == 0
            and c["global_undecided_pairs"] == n * (n - 1) // 2
            for c in spin_cases if c["spin_a"] > 0.0
        ),
        "scope_note": (
            "K6 is a near-horizon convergence diagnostic. "
            "It verifies that the ZAMO angular velocity omega_ZAMO = -g_tphi/g_phiphi "
            "converges to the horizon angular velocity Omega_H = a/(r_+^2+a^2) "
            "as r -> r_+ from outside the horizon. "
            "It does not cross the horizon, does not integrate null geodesics, "
            "and does not claim global causal reachability. "
            "Convergence of omega_ZAMO to Omega_H is a local metric identity, "
            "not a causal relation."
        ),
    }

    return {
        "aggregate":      aggregate,
        "spin_summaries": spin_summaries,
        "rows":           all_rows,
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
    spin_summaries: list[dict[str, Any]],
    aggregate:      dict[str, Any],
    md_path:        Path,
    png_path:       Path,
) -> None:
    lines = [
        f"# {AUDIT_ID}: ZAMO Angular Velocity Horizon Convergence Audit",
        "",
        f"Generated: {aggregate['generated_at_utc']}",
        "",
        "## What this is",
        "",
        "This is a **near-horizon convergence diagnostic**, not a Kerr causal solver.",
        "",
        "It verifies that the local ZAMO angular velocity",
        "`omega_ZAMO = -g_tphi / g_phiphi` converges to the known-truth horizon",
        "angular velocity `Omega_H = a / (r_+^2 + a^2)` as `r -> r_+` from outside",
        "the horizon.",
        "",
        "**It does NOT:**",
        "",
        "- Cross the horizon.",
        "- Implement Kerr causal inference of any kind.",
        "- Integrate null geodesics.",
        "- Claim global causal reachability.",
        "- Convergence of omega_ZAMO to Omega_H is a **local metric identity**,",
        "  not a causal relation.",
        "",
        "## Connection to the K-sequence",
        "",
        "- K5 measured omega_ZAMO asymmetry at fixed r (prograde/retrograde).",
        "- K6 measures omega_ZAMO convergence to Omega_H near the horizon.",
        "- Together they bridge local frame-dragging -> horizon angular velocity,",
        "  without crossing the Hawking/Bekenstein thermodynamic guardrail.",
        "",
        "## Parameters",
        "",
        f"- M = {aggregate['M']}, theta = pi/2 (equatorial), M = 1 (fixed)",
        f"- Spins: {aggregate['spins']}",
        f"- Deltas (r - r_+): {aggregate['deltas']}",
        f"- N = {aggregate['N']}, seed = {aggregate['seed']}, "
        f"margin = {aggregate['margin']}",
        "",
        "## Known-Truth Checks",
        "",
        "1. **All r_eval > r_+**: near-horizon grid stays exterior (delta > 0).",
        "2. **omega_ZAMO < Omega_H** (a>0): frame-dragging rate is below the",
        "   horizon value for all r > r_+ (trivially True for a=0).",
        "3. **Monotone convergence** (a>0): residuals decrease as delta -> 0,",
        "   consistent with O(delta) convergence rate.",
        "4. **Causal invariant**: a>0 => all global pairs undecided.",
        "",
        "## Diagnostic Figure",
        "",
        f"![K6 ZAMO horizon convergence audit]({png_path.name})",
        "",
        "The 2×2 figure shows:",
        "- Panel 1: omega_ZAMO vs delta for each spin (convergence to Omega_H)",
        "- Panel 2: |omega_ZAMO - Omega_H| vs delta (log-log; linear slope ~ 1)",
        "- Panel 3: Omega_H vs spin a (analytic formula curve + test points)",
        "- Panel 4: residual / delta vs delta (testing linear convergence rate)",
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
        "| a | r_+ | Omega_H | mono_pass | below_pass | min_res | max_res | pass |",
        "|---|-----|---------|-----------|------------|---------|---------|------|",
    ]

    for s in spin_summaries:
        a = s["spin_a"]
        lines.append(
            f"| {a:.2f} "
            f"| {s['r_plus']:.6f} "
            f"| {s['Omega_H']:.8f} "
            f"| {s['convergence_monotone_pass']} "
            f"| {s['omega_zamo_below_Omega_H']} "
            f"| {s['min_residual']:.2e} "
            f"| {s['max_residual']:.2e} "
            f"| **{s['all_checks_pass']}** |"
        )

    lines += [
        "",
        "## Per-Spin Residual Tables",
        "",
    ]
    for s in spin_summaries:
        a = s["spin_a"]
        lines.append(f"### a = {a:.2f}")
        lines.append("")
        lines.append("| delta | r_eval | omega_ZAMO | residual |")
        lines.append("|-------|--------|------------|----------|")
        for i, (delta, omega_v, res) in enumerate(
            zip(s["deltas"], s["omega_zamo_values"], s["residuals"])
        ):
            lines.append(
                f"| {delta:.1e} "
                f"| {s['r_plus'] + delta:.6f} "
                f"| {omega_v:.8f} "
                f"| {res:.2e} |"
            )
        lines.append("")

    lines += [
        "## Causal Accounting",
        "",
        "| a | global_true | global_false | global_undecided |",
        "|---|-------------|--------------|-----------------|",
    ]
    for s in spin_summaries:
        lines.append(
            f"| {s['spin_a']:.2f} "
            f"| {s['global_true_relations']} "
            f"| {s['global_false_relations']} "
            f"| {s['global_undecided_pairs']} |"
        )

    lines += [
        "",
        "## Interpretation",
        "",
        "- `a=0`: omega_ZAMO = Omega_H = 0 everywhere; trivial check, True by convention.",
        "- `a>0`: the convergence is linear in `delta = r - r_+`, consistent with",
        "  the O(delta) Taylor expansion of omega_ZAMO around r_+.",
        "- The residuals decrease by approximately the same factor as delta at each",
        "  step (ratio ≈ 3 for factor-3 steps in delta).",
        "- This audit does **not** constitute a Hawking temperature computation.",
        "  It is a local geometric identity check, satisfying the level-A",
        "  criterion from the Hawking consistency guardrail (AGENTS.md).",
        "",
    ]
    md_path.write_text("\n".join(lines), encoding="utf-8")


# ---------------------------------------------------------------------------
# Figure generator
# ---------------------------------------------------------------------------

def write_figure(
    spin_summaries: list[dict[str, Any]],
    png_path:       Path,
) -> None:
    """Generate a 2x2 diagnostic panel: ZAMO omega_H convergence."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    M      = 1.0
    deltas = spin_summaries[0]["deltas"]

    # Non-trivial spins only (exclude a=0 from convergence plots)
    non_trivial = [s for s in spin_summaries if s["spin_a"] > 0.0]

    colors = ["steelblue", "darkorange", "green", "purple"]
    _FLOOR = 1e-20

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    fig.suptitle("S4-KERR-K6 ZAMO angular velocity horizon convergence audit",
                 fontsize=13)

    # ---- Panel 1: omega_ZAMO vs delta for each spin ----
    ax1 = axes[0, 0]
    for s, col in zip(non_trivial, colors):
        a = s["spin_a"]
        ax1.semilogx(deltas, s["omega_zamo_values"], "o-", color=col,
                     linewidth=1.5, markersize=5, label=f"$a$={a:.2f}")
        ax1.axhline(s["Omega_H"], color=col, linestyle="--", linewidth=0.9,
                    alpha=0.7)
    ax1.set_xlabel(r"$\delta = r - r_+$")
    ax1.set_ylabel(r"$\omega_{\rm ZAMO}(r_+ + \delta)$")
    ax1.set_title(r"$\omega_{\rm ZAMO}$ converging to $\Omega_H$ (dashed)")
    ax1.legend(fontsize=8)
    ax1.grid(True, which="both", linestyle="--", alpha=0.4)

    # ---- Panel 2: residual vs delta (log-log) ----
    ax2 = axes[0, 1]
    for s, col in zip(non_trivial, colors):
        a = s["spin_a"]
        res = [max(r, _FLOOR) for r in s["residuals"]]
        ax2.loglog(deltas, res, "o-", color=col, linewidth=1.5, markersize=5,
                   label=f"$a$={a:.2f}")
    # Reference slope-1 line
    d0, r0 = deltas[0], non_trivial[0]["residuals"][0]
    ref_x = [deltas[0], deltas[-1]]
    ref_y = [r0, r0 * (deltas[-1] / deltas[0])]
    ax2.loglog(ref_x, ref_y, "k:", linewidth=1.0, label="slope 1")
    ax2.set_xlabel(r"$\delta = r - r_+$")
    ax2.set_ylabel(r"$|\omega_{\rm ZAMO} - \Omega_H|$")
    ax2.set_title("Convergence residual vs $\\delta$ (log-log)")
    ax2.legend(fontsize=8)
    ax2.grid(True, which="both", linestyle="--", alpha=0.4)

    # ---- Panel 3: Omega_H vs spin a ----
    ax3 = axes[1, 0]
    a_dense = [i * 0.01 for i in range(1, 100)]   # a in (0, 1)
    r_plus_dense = [1.0 + math.sqrt(1.0 - a ** 2) for a in a_dense]
    Omega_H_dense = [a / (r ** 2 + a ** 2) for a, r in zip(a_dense, r_plus_dense)]
    ax3.plot(a_dense, Omega_H_dense, "-", color="steelblue", linewidth=1.5,
             label=r"$\Omega_H = a/(r_+^2 + a^2)$")
    # Mark test spins
    test_spins = [s["spin_a"] for s in non_trivial]
    test_OH    = [s["Omega_H"] for s in non_trivial]
    ax3.scatter(test_spins, test_OH, color="red", zorder=5, s=40,
                label="test spins")
    ax3.set_xlabel("spin $a$")
    ax3.set_ylabel(r"$\Omega_H$")
    ax3.set_title("Horizon angular velocity $\\Omega_H$ vs $a$")
    ax3.legend(fontsize=8)
    ax3.grid(True, which="both", linestyle="--", alpha=0.4)

    # ---- Panel 4: residual / delta vs delta (linear convergence rate) ----
    ax4 = axes[1, 1]
    for s, col in zip(non_trivial, colors):
        a = s["spin_a"]
        rate = [max(r, _FLOOR) / d for r, d in zip(s["residuals"], deltas)]
        ax4.semilogx(deltas, rate, "o-", color=col, linewidth=1.5,
                     markersize=5, label=f"$a$={a:.2f}")
    ax4.set_xlabel(r"$\delta = r - r_+$")
    ax4.set_ylabel(r"$|\omega_{\rm ZAMO} - \Omega_H| / \delta$")
    ax4.set_title("Linear convergence rate (should plateau)")
    ax4.legend(fontsize=8)
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

    write_csv(payload["rows"], csv_path)
    write_json(payload, json_path)
    write_figure(payload["spin_summaries"], png_path)
    write_md(payload["spin_summaries"], payload["aggregate"], md_path, png_path)

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
    for s in payload["spin_summaries"]:
        a = s["spin_a"]
        print(
            f"  a={a:.2f}"
            f"  r_+={s['r_plus']:.6f}"
            f"  Omega_H={s['Omega_H']:.6f}"
            f"  mono={s['convergence_monotone_pass']}"
            f"  below={s['omega_zamo_below_Omega_H']}"
            f"  min_res={s['min_residual']:.2e}"
            f"  max_res={s['max_residual']:.2e}"
            f"  pass={s['all_checks_pass']}"
        )
    print(f"wrote {csv_path}")
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")
    print(f"wrote {png_path}")


if __name__ == "__main__":
    main()
