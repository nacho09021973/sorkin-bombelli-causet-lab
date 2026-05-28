#!/usr/bin/env python3
"""S4-KERR-K9-EQUATORIAL-FULL-RHS-PREFLIGHT-001 (Step 1/4).

Scope in this step:
- Equatorial full RHS core for Kerr null geodesics.
- Local fixed-step RK4 integrator.
- Schwarzschild radial gate only: a=0, b=0, outgoing/ingoing.
- Stdout summary only (no artifacts).
"""

from __future__ import annotations

import csv
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_kerr_minimal_benchmark import kerr_horizon_radius  # noqa: E402
from explore.sorkin4_schwarzschild_benchmark import (  # noqa: E402
    run_schwarzschild_minimal_benchmark as schwarz,
)


NULL_RESIDUAL_TOL = 1.0e-7
EL_RESIDUAL_TOL = 1.0e-7
RADIAL_RHS_TOL = 1.0e-6
R_MIN_TOL = 1.0e-10
N_STEPS = 80
DLAMBDA = 0.01
CIRCULAR_DRIFT_TOL = 1.0e-6
N_EVENTS = 12
ARTIFACT_DIR = Path(__file__).resolve().parent
OUT_PREFIX = "kerr_k9_equatorial_full_rhs_preflight_001_n12_seed1959"
COMMAND = (
    "python3 explore/sorkin4_kerr_benchmark/"
    "audit_kerr_k9_equatorial_full_rhs_preflight_001.py"
)


def delta(r: float, mass: float, spin: float) -> float:
    return r * r - 2.0 * mass * r + spin * spin


def radial_potential(
    r: float, mass: float, spin: float, impact_b: float, energy: float = 1.0
) -> float:
    p_r = (r * r + spin * spin) * energy - spin * impact_b
    l_minus_ae = impact_b - spin * energy
    return p_r * p_r - delta(r, mass, spin) * l_minus_ae * l_minus_ae


def photon_sphere_radius_pro(mass: float, spin: float) -> float:
    return 2.0 * mass * (1.0 + math.cos((2.0 / 3.0) * math.acos(-spin / mass)))


def photon_sphere_radius_retro(mass: float, spin: float) -> float:
    return 2.0 * mass * (1.0 + math.cos((2.0 / 3.0) * math.acos(+spin / mass)))


def photon_impact_parameter(
    r_ph: float, mass: float, spin: float, prograde: bool
) -> float:
    dlt = delta(r_ph, mass, spin)
    if dlt <= 0.0:
        raise ValueError(f"Delta<=0 at circular orbit r={r_ph:.9g}")
    sign = +1.0 if prograde else -1.0
    return spin + sign * 2.0 * r_ph * math.sqrt(dlt) / (r_ph - mass)


def kerr_equatorial_rhs(
    state: tuple[float, float, float],
    mass: float,
    spin: float,
    impact_b: float,
    direction: float,
    energy: float = 1.0,
) -> tuple[float, float, float]:
    _, r, _ = state
    sigma = r * r
    dlt = delta(r, mass, spin)
    if dlt <= 0.0:
        raise ValueError(f"Delta<=0 at r={r:.9g}")
    p_r = (r * r + spin * spin) * energy - spin * impact_b
    l_minus_ae = impact_b - spin * energy
    r_pot = radial_potential(r, mass, spin, impact_b, energy=energy)
    tdot = (((r * r + spin * spin) / dlt) * p_r + spin * l_minus_ae) / sigma
    phidot = ((spin / dlt) * p_r + l_minus_ae) / sigma
    rdot = direction * math.sqrt(max(r_pot, 0.0)) / sigma
    return (tdot, rdot, phidot)


def rk4_step(
    state: tuple[float, float, float],
    h: float,
    mass: float,
    spin: float,
    impact_b: float,
    direction: float,
    energy: float = 1.0,
) -> tuple[float, float, float]:
    def add(s: tuple[float, float, float], k: tuple[float, float, float], c: float) -> tuple[float, float, float]:
        return (s[0] + c * k[0], s[1] + c * k[1], s[2] + c * k[2])

    k1 = kerr_equatorial_rhs(state, mass, spin, impact_b, direction, energy=energy)
    k2 = kerr_equatorial_rhs(add(state, k1, 0.5 * h), mass, spin, impact_b, direction, energy=energy)
    k3 = kerr_equatorial_rhs(add(state, k2, 0.5 * h), mass, spin, impact_b, direction, energy=energy)
    k4 = kerr_equatorial_rhs(add(state, k3, h), mass, spin, impact_b, direction, energy=energy)
    return (
        state[0] + (h / 6.0) * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]),
        state[1] + (h / 6.0) * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]),
        state[2] + (h / 6.0) * (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]),
    )


def integrate_trajectory(
    state0: tuple[float, float, float],
    mass: float,
    spin: float,
    impact_b: float,
    direction: float,
    energy: float = 1.0,
    n_steps: int = N_STEPS,
    h: float = DLAMBDA,
    r_horizon_safety: float = 1.0e-6,
    r_potential_tol: float = R_MIN_TOL,
) -> dict[str, Any]:
    r_plus = kerr_horizon_radius(mass, spin)
    states = [state0]
    rhs_vals = []
    failed_reason = None
    for _ in range(n_steps):
        st = states[-1]
        rr = st[1]
        rpot = radial_potential(rr, mass, spin, impact_b, energy=energy)
        if rpot < -r_potential_tol:
            failed_reason = f"R< -tol at r={rr:.9g}, R={rpot:.9g}"
            break
        if rr <= r_plus + r_horizon_safety:
            failed_reason = f"r reached safety boundary at r={rr:.9g}"
            break
        rhs = kerr_equatorial_rhs(st, mass, spin, impact_b, direction, energy=energy)
        rhs_vals.append(rhs)
        if not all(math.isfinite(x) for x in rhs):
            failed_reason = "non-finite RHS"
            break
        nxt = rk4_step(st, h, mass, spin, impact_b, direction, energy=energy)
        if not all(math.isfinite(x) for x in nxt):
            failed_reason = "non-finite solution"
            break
        if nxt[1] <= r_plus + r_horizon_safety:
            failed_reason = f"next step approaches horizon safety boundary r={nxt[1]:.9g}"
            break
        states.append(nxt)
    return {"states": states, "rhs": rhs_vals, "failed_reason": failed_reason}


def metric_equatorial(r: float, mass: float, spin: float) -> tuple[float, float, float, float]:
    g_tt = -(1.0 - (2.0 * mass / r))
    g_tphi = -(2.0 * mass * spin / r)
    g_rr = (r * r) / delta(r, mass, spin)
    g_phiphi = r * r + spin * spin + (2.0 * mass * spin * spin / r)
    return (g_tt, g_tphi, g_rr, g_phiphi)


def null_residual(
    r: float, tdot: float, rdot: float, phidot: float, mass: float, spin: float
) -> float:
    g_tt, g_tphi, g_rr, g_phiphi = metric_equatorial(r, mass, spin)
    return (
        g_tt * tdot * tdot
        + 2.0 * g_tphi * tdot * phidot
        + g_rr * rdot * rdot
        + g_phiphi * phidot * phidot
    )


def conserved_energy_angular_momentum(
    r: float, tdot: float, phidot: float, mass: float, spin: float
) -> tuple[float, float]:
    g_tt, g_tphi, _, g_phiphi = metric_equatorial(r, mass, spin)
    e_calc = -g_tt * tdot - g_tphi * phidot
    l_calc = g_tphi * tdot + g_phiphi * phidot
    return (e_calc, l_calc)


def check_trajectory(
    states: list[tuple[float, float, float]],
    rhs_vals: list[tuple[float, float, float]],
    mass: float,
    spin: float,
    impact_b: float,
    direction: float,
    energy: float,
    r_plus: float,
    apply_schwarzschild_limit: bool,
    circular_target_r: float | None = None,
    circular_advisory: bool = False,
) -> dict[str, Any]:
    rs = [s[1] for s in states]
    ts = [s[0] for s in states]
    all_points_exterior = all(r > r_plus for r in rs)
    min_delta = min(delta(r, mass, spin) for r in rs)
    min_rpot = min(radial_potential(r, mass, spin, impact_b, energy=energy) for r in rs)
    finite_solution = all(all(math.isfinite(v) for v in s) for s in states)
    finite_rhs = all(all(math.isfinite(v) for v in rr) for rr in rhs_vals)
    t_monotonic = all(ts[i + 1] >= ts[i] for i in range(len(ts) - 1))
    radial_sign = all(
        (states[i + 1][1] - states[i][1]) >= -1.0e-12 if direction > 0 else (states[i + 1][1] - states[i][1]) <= 1.0e-12
        for i in range(len(states) - 1)
    )

    max_radial_rhs_res = 0.0
    max_null_res = 0.0
    max_e_res = 0.0
    max_l_res = 0.0
    max_dr_exact_res = 0.0
    max_dphi_abs = 0.0
    max_dt_exact_res = 0.0
    max_circular_drift = 0.0

    for i, st in enumerate(states[:-1]):
        tdot, rdot, phidot = rhs_vals[i]
        r = st[1]
        rpot = radial_potential(r, mass, spin, impact_b, energy=energy)
        radial_rhs = direction * math.sqrt(max(rpot, 0.0)) / (r * r)
        max_radial_rhs_res = max(max_radial_rhs_res, abs(rdot - radial_rhs))
        max_null_res = max(max_null_res, abs(null_residual(r, tdot, rdot, phidot, mass, spin)))
        e_calc, l_calc = conserved_energy_angular_momentum(r, tdot, phidot, mass, spin)
        max_e_res = max(max_e_res, abs(e_calc - energy))
        max_l_res = max(max_l_res, abs(l_calc - impact_b))
        max_dr_exact_res = max(max_dr_exact_res, abs(rdot - direction * energy))
        max_dphi_abs = max(max_dphi_abs, abs(phidot))
        max_dt_exact_res = max(max_dt_exact_res, abs(tdot - (energy / (1.0 - (2.0 * mass / r)))))
        if circular_target_r is not None:
            max_circular_drift = max(max_circular_drift, abs(r - circular_target_r))

    schwarzschild_radial_limit_pass = (
        max_dphi_abs <= EL_RESIDUAL_TOL
        and max_dr_exact_res <= RADIAL_RHS_TOL
        and max_dt_exact_res <= EL_RESIDUAL_TOL
    )
    circular_orbit_radial_drift_pass = (
        max_circular_drift <= CIRCULAR_DRIFT_TOL if circular_target_r is not None else True
    )

    checks = {
        "all_points_exterior": all_points_exterior,
        "min_Delta_gt_zero": min_delta > 0.0,
        "min_R_ge_neg_tol": min_rpot >= -R_MIN_TOL,
        "finite_rhs_all_steps": finite_rhs,
        "finite_solution_all_steps": finite_solution,
        "t_monotonic_future_pass": t_monotonic,
        "radial_sign_consistency_pass": radial_sign,
        "radial_rhs_consistency_pass": max_radial_rhs_res <= RADIAL_RHS_TOL,
        "null_condition_pass": max_null_res <= NULL_RESIDUAL_TOL,
        "constants_consistency_pass": (max_e_res <= EL_RESIDUAL_TOL and max_l_res <= EL_RESIDUAL_TOL),
        "schwarzschild_radial_limit_pass": (schwarzschild_radial_limit_pass if apply_schwarzschild_limit else True),
        "circular_orbit_radial_drift_pass": (
            circular_orbit_radial_drift_pass if (circular_target_r is not None and not circular_advisory) else True
        ),
    }
    all_checks_pass = all(checks.values())
    return {
        "checks": checks,
        "all_checks_pass": all_checks_pass,
        "circular_advisory": circular_advisory,
        "metrics": {
            "min_Delta": min_delta,
            "min_R": min_rpot,
            "max_abs_null_residual": max_null_res,
            "max_abs_E_residual": max_e_res,
            "max_abs_L_residual": max_l_res,
            "max_abs_radial_rhs_residual": max_radial_rhs_res,
            "max_abs_dphi": max_dphi_abs,
            "max_abs_dr_pm1_residual": max_dr_exact_res,
            "max_abs_dt_schw_residual": max_dt_exact_res,
            "max_abs_circular_drift": max_circular_drift,
        },
    }


def _run_case(
    *,
    label: str,
    mass: float,
    spin: float,
    impact_b: float,
    direction: float,
    r0: float,
    n_steps: int,
    h: float,
    apply_schwarzschild_limit: bool,
    circular_target_r: float | None = None,
    circular_advisory: bool = False,
) -> dict[str, Any]:
    energy = 1.0
    r_plus = kerr_horizon_radius(mass, spin)
    integ = integrate_trajectory(
        state0=(0.0, r0, 0.0),
        mass=mass,
        spin=spin,
        impact_b=impact_b,
        direction=direction,
        energy=energy,
        n_steps=n_steps,
        h=h,
        r_horizon_safety=1.0e-6,
        r_potential_tol=R_MIN_TOL,
    )
    fail_reason = integ["failed_reason"]
    states = integ["states"]
    rhs_vals = integ["rhs"]
    checked = check_trajectory(
        states,
        rhs_vals,
        mass,
        spin,
        impact_b,
        direction,
        energy,
        r_plus,
        apply_schwarzschild_limit=apply_schwarzschild_limit,
        circular_target_r=circular_target_r,
        circular_advisory=circular_advisory,
    )
    checks = checked["checks"]
    metrics = checked["metrics"]
    all_checks_pass = checked["all_checks_pass"]
    case_pass = (fail_reason is None) and all_checks_pass
    print(
        f"[{label}] pass={case_pass} steps={len(states)-1}/{n_steps} "
        f"a={spin:.2f} b={impact_b:.6f} dir={int(direction):+d} "
        f"min_Delta={metrics['min_Delta']:.3e} min_R={metrics['min_R']:.3e} "
        f"null_max={metrics['max_abs_null_residual']:.3e} "
        f"Eres_max={metrics['max_abs_E_residual']:.3e} Lres_max={metrics['max_abs_L_residual']:.3e} "
        f"dr_rhs_max={metrics['max_abs_radial_rhs_residual']:.3e}"
    )
    if apply_schwarzschild_limit:
        print(
            f"  schw: dr_pm1_max={metrics['max_abs_dr_pm1_residual']:.3e} "
            f"dt_schw_max={metrics['max_abs_dt_schw_residual']:.3e} dphi_max={metrics['max_abs_dphi']:.3e}"
        )
    if circular_target_r is not None:
        print(
            f"  circular: R_start={radial_potential(r0, mass, spin, impact_b):.3e} "
            f"drift_max={metrics['max_abs_circular_drift']:.3e} "
            f"pass={checks['circular_orbit_radial_drift_pass']} advisory={circular_advisory}"
        )
    if fail_reason is not None:
        print(f"  fail_reason={fail_reason}")
    failed_checks = [k for k, v in checks.items() if not v]
    if failed_checks:
        print(f"  failed_checks={','.join(failed_checks)}")
    return {
        "label": label,
        "mass": mass,
        "spin": spin,
        "energy": energy,
        "impact_b": impact_b,
        "direction": direction,
        "r_plus": r_plus,
        "r0": r0,
        "n_steps": n_steps,
        "h": h,
        "states": states,
        "rhs_vals": rhs_vals,
        "pass": case_pass,
        "fail_reason": fail_reason,
        "checks": checks,
        "metrics": metrics,
        "is_advisory": circular_target_r is not None and circular_advisory,
        "is_spin_positive": spin > 0.0,
        "case_type": (
            "safe_b0_radial_flow_control"
            if impact_b == 0.0
            else "circular_photon_orbit_hold_diagnostic"
        ),
        "advisory_note": (
            "advisory circular hold; not evidence of orbital stability"
            if (circular_target_r is not None and circular_advisory)
            else ""
        ),
    }


def build_preflight_results() -> list[dict[str, Any]]:
    mass = 1.0
    results = []
    print("S4-KERR-K9-EQUATORIAL-FULL-RHS-PREFLIGHT-001 | Step 2/4")
    print("Family 1+2: safe b=0 radial-flow control cases (not generic Kerr null-geodesic families).")

    for spin in (0.0, 0.25, 0.5, 0.75):
        r_plus = kerr_horizon_radius(mass, spin)
        r0 = r_plus + schwarz.EXTERIOR_MARGIN + 1.0
        for label, direction in (("outgoing", +1.0), ("ingoing", -1.0)):
            results.append(
                _run_case(
                    label=f"safe_b0_spin_{spin:.2f}_{label}",
                    mass=mass,
                    spin=spin,
                    impact_b=0.0,
                    direction=direction,
                    r0=r0,
                    n_steps=N_STEPS,
                    h=DLAMBDA,
                    apply_schwarzschild_limit=(spin == 0.0),
                )
            )

    print("Family 2 optional: a=0.90 skipped in Step 2 (near-horizon margin sensitivity).")
    print("Family 3: circular photon orbit hold (short-interval, advisory where needed; no stability claim).")
    for spin in (0.0, 0.25, 0.5, 0.75):
        r_plus = kerr_horizon_radius(mass, spin)
        for suffix, prograde in (("pro", True), ("retro", False)):
            r_ph = photon_sphere_radius_pro(mass, spin) if prograde else photon_sphere_radius_retro(mass, spin)
            b_ph = photon_impact_parameter(r_ph, mass, spin, prograde=prograde)
            advisory = spin >= 0.75
            if r_ph <= r_plus + 5.0e-4:
                print(f"[circular_spin_{spin:.2f}_{suffix}] skipped: not clearly exterior.")
                continue
            results.append(
                _run_case(
                    label=f"circular_spin_{spin:.2f}_{suffix}",
                    mass=mass,
                    spin=spin,
                    impact_b=b_ph,
                    direction=+1.0,
                    r0=r_ph,
                    n_steps=40,
                    h=0.002,
                    apply_schwarzschild_limit=False,
                    circular_target_r=r_ph,
                    circular_advisory=advisory,
                )
            )

    print("Family 4 optional: skipped (no additional non-circular +/-b cases added in Step 2).")
    return results


def _case_record(case: dict[str, Any]) -> dict[str, Any]:
    states = case["states"]
    metrics = case["metrics"]
    checks = case["checks"]
    direction_label = "outgoing" if case["direction"] > 0 else "ingoing"
    min_r = min(st[1] for st in states)
    max_r = max(st[1] for st in states)
    radial_drift_abs = abs(states[-1][1] - case["r0"])
    return {
        "case_id": case["label"],
        "spin_a": case["spin"],
        "M": case["mass"],
        "E": case["energy"],
        "b": case["impact_b"],
        "direction": direction_label,
        "case_type": case["case_type"],
        "r_plus": case["r_plus"],
        "r0": case["r0"],
        "lambda_end": case["n_steps"] * case["h"],
        "n_steps": case["n_steps"],
        "min_r": min_r,
        "max_r": max_r,
        "min_Delta": metrics["min_Delta"],
        "min_R": metrics["min_R"],
        "max_abs_null_residual": metrics["max_abs_null_residual"],
        "max_abs_E_residual": metrics["max_abs_E_residual"],
        "max_abs_L_residual": metrics["max_abs_L_residual"],
        "max_abs_radial_rhs_residual": metrics["max_abs_radial_rhs_residual"],
        "max_abs_schwarzschild_radial_residual": metrics["max_abs_dr_pm1_residual"],
        "radial_drift_abs": radial_drift_abs,
        "all_points_exterior": checks["all_points_exterior"],
        "finite_rhs_all_steps": checks["finite_rhs_all_steps"],
        "finite_solution_all_steps": checks["finite_solution_all_steps"],
        "t_monotonic_future_pass": checks["t_monotonic_future_pass"],
        "radial_sign_consistency_pass": checks["radial_sign_consistency_pass"],
        "radial_rhs_consistency_pass": checks["radial_rhs_consistency_pass"],
        "null_condition_pass": checks["null_condition_pass"],
        "constants_consistency_pass": checks["constants_consistency_pass"],
        "schwarzschild_radial_limit_pass": checks["schwarzschild_radial_limit_pass"],
        "circular_orbit_radial_drift_pass": checks["circular_orbit_radial_drift_pass"],
        "advisory_only": case["is_advisory"],
        "all_checks_pass": case["pass"],
    }


def write_csv(records: list[dict[str, Any]], path: Path) -> None:
    fields = [
        "case_id", "spin_a", "M", "E", "b", "direction", "case_type", "r_plus", "r0",
        "lambda_end", "n_steps", "min_r", "max_r", "min_Delta", "min_R",
        "max_abs_null_residual", "max_abs_E_residual", "max_abs_L_residual",
        "max_abs_radial_rhs_residual", "max_abs_schwarzschild_radial_residual",
        "radial_drift_abs", "all_points_exterior", "finite_rhs_all_steps",
        "finite_solution_all_steps", "t_monotonic_future_pass",
        "radial_sign_consistency_pass", "radial_rhs_consistency_pass",
        "null_condition_pass", "constants_consistency_pass",
        "schwarzschild_radial_limit_pass", "circular_orbit_radial_drift_pass",
        "advisory_only", "all_checks_pass",
    ]
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for rec in records:
            writer.writerow(rec)


def write_json(records: list[dict[str, Any]], summary: dict[str, Any], path: Path) -> None:
    payload = {
        "benchmark": "S4-KERR-K9 equatorial full RHS preflight",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "command": COMMAND,
        "scope_note": (
            "K9 integrates equatorial Kerr null RHS t(lambda), r(lambda), phi(lambda) "
            "as a preflight diagnostic only. No shooting, no causal reachability, "
            "no event-pair classification."
        ),
        "cases": records,
        "global_summary": summary,
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def write_md(records: list[dict[str, Any]], summary: dict[str, Any], path: Path) -> None:
    advisory = [r for r in records if r["advisory_only"]]
    lines = [
        "# S4-KERR-K9 equatorial full RHS preflight",
        "",
        "K9 integrates the full equatorial RHS `t(lambda)`, `r(lambda)`, and `phi(lambda)`.",
        "K9 is not point-to-point shooting.",
        "K9 does not decide causal reachability.",
        "K9 does not classify sprinkled event pairs.",
        "K9 is a preflight for a future geodesic integrator/shooter.",
        "The `b=0` cases are safe radial-flow control cases, not generic Kerr null-geodesic families.",
        "Circular photon orbit cases are drift/hold diagnostics, not stability claims.",
        "Advisory circular cases, when present, are not used as evidence of orbital stability.",
        "",
        f"- Total cases: {summary['total_cases']}",
        f"- Passed cases: {summary['passed_cases']}",
        f"- Failed cases: {summary['failed_cases']}",
        f"- Advisory cases: {summary['advisory_cases']}",
        f"- Global undecided pairs (a>0 control accounting): {summary['global_undecided_pairs']}",
        "",
        "## Artifact Set",
        "",
        f"- `{OUT_PREFIX}.csv`",
        f"- `{OUT_PREFIX}.json`",
        f"- `{OUT_PREFIX}.md`",
        f"- `{OUT_PREFIX}.png`",
        "",
        "## Advisory Cases",
        "",
    ]
    if advisory:
        for row in advisory:
            lines.append(f"- `{row['case_id']}`")
    else:
        lines.append("- None")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_png(results: list[dict[str, Any]], records: list[dict[str, Any]], path: Path) -> None:
    fig, axs = plt.subplots(2, 2, figsize=(11, 8))
    fig.suptitle("S4-KERR-K9 equatorial full RHS preflight")

    ax = axs[0, 0]
    for row, raw in zip(records, results):
        if row["case_type"] == "safe_b0_radial_flow_control":
            lam = [i * raw["h"] for i in range(len(raw["states"]))]
            rr = [s[1] for s in raw["states"]]
            ax.plot(lam, rr, label=row["case_id"])
    ax.set_xlabel("lambda")
    ax.set_ylabel("r(lambda)")
    ax.set_title("Safe b=0 radial-flow controls")
    ax.legend(fontsize=6, ncol=2)

    ax = axs[0, 1]
    for row, raw in zip(records, results):
        if row["case_type"] == "safe_b0_radial_flow_control":
            lam = [i * raw["h"] for i in range(len(raw["states"]))]
            ph = [s[2] for s in raw["states"]]
            ax.plot(lam, ph, label=row["case_id"])
    ax.set_xlabel("lambda")
    ax.set_ylabel("phi(lambda)")
    ax.set_title("Integrated coordinate output only")

    ax = axs[1, 0]
    x = list(range(len(records)))
    nulls = [r["max_abs_null_residual"] for r in records]
    ers = [r["max_abs_E_residual"] for r in records]
    lrs = [r["max_abs_L_residual"] for r in records]
    ax.semilogy(x, nulls, "o-", label="null residual")
    ax.semilogy(x, ers, "s-", label="E residual")
    ax.semilogy(x, lrs, "^-", label="L residual")
    ax.set_xlabel("case index")
    ax.set_title("Residual diagnostics by case")
    ax.legend(fontsize=7)

    ax = axs[1, 1]
    circ = [r for r in records if "circular" in r["case_id"]]
    spins = [r["spin_a"] for r in circ]
    drift = [r["radial_drift_abs"] for r in circ]
    ax.semilogy(spins, drift, "o")
    ax.set_xlabel("spin a")
    ax.set_ylabel("|r_end-r0|")
    ax.set_title("Circular hold radial drift by spin")

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def run_preflight_step3() -> int:
    results = build_preflight_results()

    non_advisory_failures = [
        r for r in results if (not r["pass"]) and (not r["is_advisory"])
    ]
    if non_advisory_failures:
        print("STOP: non-advisory case failure detected.")
        for rr in non_advisory_failures:
            print(f"  failed_case={rr['label']}")
        return 1

    records = [_case_record(r) for r in results]
    total_cases = len(records)
    passed_cases = sum(1 for r in records if r["all_checks_pass"])
    failed_cases = total_cases - passed_cases
    advisory_cases = sum(1 for r in records if r["advisory_only"])
    a_positive_cases = [r for r in results if r["is_spin_positive"]]
    undecided_pairs = N_EVENTS * (N_EVENTS - 1) // 2
    summary = {
        "total_cases": total_cases,
        "passed_cases": passed_cases,
        "failed_cases": failed_cases,
        "advisory_cases": advisory_cases,
        "global_true_relations": 0,
        "global_false_relations": 0,
        "global_undecided_pairs": undecided_pairs,
        "all_checks_pass": failed_cases == 0,
    }

    csv_path = ARTIFACT_DIR / f"{OUT_PREFIX}.csv"
    json_path = ARTIFACT_DIR / f"{OUT_PREFIX}.json"
    md_path = ARTIFACT_DIR / f"{OUT_PREFIX}.md"
    png_path = ARTIFACT_DIR / f"{OUT_PREFIX}.png"
    write_csv(records, csv_path)
    write_json(records, summary, json_path)
    write_md(records, summary, md_path)
    write_png(results, records, png_path)

    print("Global summary (stdout only):")
    print(f"  total_cases={total_cases}")
    print(f"  passed_cases={passed_cases}")
    print(f"  failed_cases={failed_cases}")
    print(f"  advisory_cases={advisory_cases}")
    print(f"  global_true_relations=0")
    print(f"  global_false_relations=0")
    print(f"  global_undecided_pairs={undecided_pairs}")
    print(f"  a_gt_0_case_count={len(a_positive_cases)}")
    print(f"  all_checks_pass={failed_cases == 0}")
    print(f"  artifacts_written={csv_path.name},{json_path.name},{md_path.name},{png_path.name}")
    return 0


def main() -> None:
    raise SystemExit(run_preflight_step3())


if __name__ == "__main__":
    main()
