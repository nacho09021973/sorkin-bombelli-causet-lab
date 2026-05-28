#!/usr/bin/env python3
"""S4-KERR-K10-EQUATORIAL-GEODESIC-SEGMENT-AUDIT-001.

Initial-value equatorial Kerr null-geodesic segment audit.
This is not point-to-point shooting, not a boundary-value solver,
and not a causal reachability classifier.
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
from audit_kerr_k9_equatorial_full_rhs_preflight_001 import (  # noqa: E402
    EL_RESIDUAL_TOL,
    NULL_RESIDUAL_TOL,
    RADIAL_RHS_TOL,
    R_MIN_TOL,
    conserved_energy_angular_momentum,
    delta,
    kerr_equatorial_rhs,
    null_residual,
    photon_impact_parameter,
    photon_sphere_radius_pro,
    photon_sphere_radius_retro,
    radial_potential,
    rk4_step,
)


N_EVENTS = 12
ARTIFACT_DIR = Path(__file__).resolve().parent
OUT_PREFIX = "kerr_k10_equatorial_geodesic_segment_audit_001_n12_seed1959"
COMMAND = (
    "python3 explore/sorkin4_kerr_benchmark/"
    "audit_kerr_k10_equatorial_geodesic_segment_audit_001.py"
)

N_STEPS = 100
DLAMBDA = 0.01
MASS = 1.0
ENERGY = 1.0
HORIZON_SAFETY = 1.0e-6
NEAR_PHOTON_STEPS = 60
NEAR_PHOTON_H = 0.002


def integrate_segment(
    *,
    state0: tuple[float, float, float],
    mass: float,
    spin: float,
    impact_b: float,
    direction: float,
    energy: float,
    n_steps: int,
    h: float,
    r_horizon_safety: float = HORIZON_SAFETY,
    r_potential_tol: float = R_MIN_TOL,
) -> dict[str, Any]:
    r_plus = kerr_horizon_radius(mass, spin)
    states = [state0]
    rhs_vals: list[tuple[float, float, float]] = []
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


def check_segment(
    *,
    states: list[tuple[float, float, float]],
    rhs_vals: list[tuple[float, float, float]],
    mass: float,
    spin: float,
    impact_b: float,
    direction: float,
    energy: float,
    r_plus: float,
    enforce_radial_sign: bool,
    apply_schwarzschild_limit: bool,
) -> dict[str, Any]:
    rs = [s[1] for s in states]
    ts = [s[0] for s in states]
    phis = [s[2] for s in states]
    all_points_exterior = all(r > r_plus for r in rs)
    min_delta = min(delta(r, mass, spin) for r in rs)
    min_rpot = min(radial_potential(r, mass, spin, impact_b, energy=energy) for r in rs)
    finite_solution = all(all(math.isfinite(v) for v in s) for s in states)
    finite_rhs = all(all(math.isfinite(v) for v in rr) for rr in rhs_vals)
    t_monotonic = all(ts[i + 1] >= ts[i] for i in range(len(ts) - 1))

    radial_sign_pass = True
    if enforce_radial_sign:
        radial_sign_pass = all(
            (states[i + 1][1] - states[i][1]) >= -1.0e-12
            if direction > 0
            else (states[i + 1][1] - states[i][1]) <= 1.0e-12
            for i in range(len(states) - 1)
        )

    max_radial_rhs_res = 0.0
    max_null_res = 0.0
    max_e_res = 0.0
    max_l_res = 0.0
    max_dr_exact_res = 0.0
    max_dphi_abs = 0.0
    max_dt_exact_res = 0.0

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

    schwarzschild_radial_limit_pass = (
        max_dphi_abs <= EL_RESIDUAL_TOL
        and max_dr_exact_res <= RADIAL_RHS_TOL
        and max_dt_exact_res <= EL_RESIDUAL_TOL
    )

    delta_t = ts[-1] - ts[0]
    delta_phi = phis[-1] - phis[0]

    checks = {
        "all_points_exterior": all_points_exterior,
        "min_Delta_gt_zero": min_delta > 0.0,
        "min_R_ge_neg_tol": min_rpot >= -R_MIN_TOL,
        "finite_rhs_all_steps": finite_rhs,
        "finite_solution_all_steps": finite_solution,
        "t_monotonic_future_pass": t_monotonic,
        "radial_rhs_consistency_pass": max_radial_rhs_res <= RADIAL_RHS_TOL,
        "null_condition_pass": max_null_res <= NULL_RESIDUAL_TOL,
        "constants_consistency_pass": (max_e_res <= EL_RESIDUAL_TOL and max_l_res <= EL_RESIDUAL_TOL),
        "schwarzschild_radial_limit_pass": schwarzschild_radial_limit_pass if apply_schwarzschild_limit else True,
        "radial_sign_consistency_pass": radial_sign_pass,
        "angular_accumulation_finite_pass": math.isfinite(delta_t) and math.isfinite(delta_phi),
        "no_endpoint_targeting_pass": True,
        "no_global_causal_relations_decided": True,
    }
    return {
        "checks": checks,
        "all_checks_pass": all(checks.values()),
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
            "delta_t": delta_t,
            "delta_phi": delta_phi,
            "radial_drift": rs[-1] - rs[0],
            "min_r": min(rs),
            "max_r": max(rs),
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
    case_type: str,
    enforce_radial_sign: bool,
    apply_schwarzschild_limit: bool,
    advisory_only: bool = False,
) -> dict[str, Any]:
    r_plus = kerr_horizon_radius(mass, spin)
    integ = integrate_segment(
        state0=(0.0, r0, 0.0),
        mass=mass,
        spin=spin,
        impact_b=impact_b,
        direction=direction,
        energy=ENERGY,
        n_steps=n_steps,
        h=h,
    )
    states = integ["states"]
    rhs_vals = integ["rhs"]
    fail_reason = integ["failed_reason"]
    checked = check_segment(
        states=states,
        rhs_vals=rhs_vals,
        mass=mass,
        spin=spin,
        impact_b=impact_b,
        direction=direction,
        energy=ENERGY,
        r_plus=r_plus,
        enforce_radial_sign=enforce_radial_sign,
        apply_schwarzschild_limit=apply_schwarzschild_limit,
    )
    checks = checked["checks"]
    metrics = checked["metrics"]
    case_pass = (fail_reason is None) and checked["all_checks_pass"]
    return {
        "label": label,
        "mass": mass,
        "spin": spin,
        "energy": ENERGY,
        "impact_b": impact_b,
        "direction": direction,
        "r_plus": r_plus,
        "r0": r0,
        "n_steps": n_steps,
        "h": h,
        "states": states,
        "rhs_vals": rhs_vals,
        "checks": checks,
        "metrics": metrics,
        "pass": case_pass,
        "fail_reason": fail_reason,
        "case_type": case_type,
        "advisory_only": advisory_only,
    }


def _segment_record(case: dict[str, Any]) -> dict[str, Any]:
    c = case["checks"]
    m = case["metrics"]
    return {
        "case_id": case["label"],
        "spin_a": case["spin"],
        "M": case["mass"],
        "E": case["energy"],
        "b": case["impact_b"],
        "direction": "outgoing" if case["direction"] > 0 else "ingoing",
        "case_type": case["case_type"],
        "r_plus": case["r_plus"],
        "r0": case["r0"],
        "lambda_end": case["n_steps"] * case["h"],
        "n_steps": case["n_steps"],
        "delta_t": m["delta_t"],
        "delta_phi": m["delta_phi"],
        "radial_drift": m["radial_drift"],
        "min_r": m["min_r"],
        "max_r": m["max_r"],
        "min_Delta": m["min_Delta"],
        "min_R": m["min_R"],
        "max_abs_null_residual": m["max_abs_null_residual"],
        "max_abs_E_residual": m["max_abs_E_residual"],
        "max_abs_L_residual": m["max_abs_L_residual"],
        "max_abs_radial_rhs_residual": m["max_abs_radial_rhs_residual"],
        "all_points_exterior": c["all_points_exterior"],
        "finite_rhs_all_steps": c["finite_rhs_all_steps"],
        "finite_solution_all_steps": c["finite_solution_all_steps"],
        "t_monotonic_future_pass": c["t_monotonic_future_pass"],
        "radial_rhs_consistency_pass": c["radial_rhs_consistency_pass"],
        "null_condition_pass": c["null_condition_pass"],
        "constants_consistency_pass": c["constants_consistency_pass"],
        "schwarzschild_radial_limit_pass": c["schwarzschild_radial_limit_pass"],
        "radial_sign_consistency_pass": c["radial_sign_consistency_pass"],
        "angular_accumulation_finite_pass": c["angular_accumulation_finite_pass"],
        "no_endpoint_targeting_pass": c["no_endpoint_targeting_pass"],
        "all_checks_pass": case["pass"],
        "advisory_only": case["advisory_only"],
    }


def _is_safe_trial(*, spin: float, b: float, direction: float, r0: float) -> bool:
    trial = integrate_segment(
        state0=(0.0, r0, 0.0),
        mass=MASS,
        spin=spin,
        impact_b=b,
        direction=direction,
        energy=ENERGY,
        n_steps=30,
        h=0.01,
    )
    if trial["failed_reason"] is not None:
        return False
    rs = [s[1] for s in trial["states"]]
    min_r = min(rs)
    min_rpot = min(radial_potential(r, MASS, spin, b, energy=ENERGY) for r in rs)
    return min_r > kerr_horizon_radius(MASS, spin) + 1.0e-4 and min_rpot >= -R_MIN_TOL


def _pick_safe_non_circular_b(spin: float, prograde: bool) -> float | None:
    r_plus = kerr_horizon_radius(MASS, spin)
    r0 = r_plus + schwarz.EXTERIOR_MARGIN + 1.2
    candidates = [0.5, 1.0, 1.5, 2.0, 2.5] if prograde else [-0.5, -1.0, -1.5, -2.0, -2.5]
    for b in candidates:
        if _is_safe_trial(spin=spin, b=b, direction=+1.0, r0=r0):
            return b
    return None


def build_results() -> list[dict[str, Any]]:
    results: list[dict[str, Any]] = []

    # Family 1 and 2: radial control segments
    for spin in (0.0, 0.25, 0.5, 0.75):
        r_plus = kerr_horizon_radius(MASS, spin)
        r0 = r_plus + schwarz.EXTERIOR_MARGIN + 1.0
        for lbl, direction in (("outgoing", +1.0), ("ingoing", -1.0)):
            results.append(
                _run_case(
                    label=f"segment_safe_b0_spin_{spin:.2f}_{lbl}",
                    mass=MASS,
                    spin=spin,
                    impact_b=0.0,
                    direction=direction,
                    r0=r0,
                    n_steps=N_STEPS,
                    h=DLAMBDA,
                    case_type="safe_b0_control_segment",
                    enforce_radial_sign=True,
                    apply_schwarzschild_limit=(spin == 0.0),
                )
            )

    # Family 3: non-circular safe +/-b segments for a>0
    for spin in (0.25, 0.5, 0.75):
        r_plus = kerr_horizon_radius(MASS, spin)
        r0 = r_plus + schwarz.EXTERIOR_MARGIN + 1.2
        b_pro = _pick_safe_non_circular_b(spin, prograde=True)
        b_retro = _pick_safe_non_circular_b(spin, prograde=False)
        if b_pro is not None:
            results.append(
                _run_case(
                    label=f"segment_non_circular_spin_{spin:.2f}_pro_b{b_pro:+.2f}",
                    mass=MASS,
                    spin=spin,
                    impact_b=b_pro,
                    direction=+1.0,
                    r0=r0,
                    n_steps=N_STEPS,
                    h=DLAMBDA,
                    case_type="non_circular_prograde_safe_segment",
                    enforce_radial_sign=True,
                    apply_schwarzschild_limit=False,
                )
            )
        if b_retro is not None:
            results.append(
                _run_case(
                    label=f"segment_non_circular_spin_{spin:.2f}_retro_b{b_retro:+.2f}",
                    mass=MASS,
                    spin=spin,
                    impact_b=b_retro,
                    direction=+1.0,
                    r0=r0,
                    n_steps=N_STEPS,
                    h=DLAMBDA,
                    case_type="non_circular_retrograde_safe_segment",
                    enforce_radial_sign=True,
                    apply_schwarzschild_limit=False,
                )
            )

    # Family 4: near-photon-sphere diagnostic segments (a>0)
    for spin in (0.25, 0.5, 0.75):
        r_ph_pro = photon_sphere_radius_pro(MASS, spin)
        b_ph_pro = photon_impact_parameter(r_ph_pro, MASS, spin, prograde=True)
        r0_pro = r_ph_pro + 0.02
        b_near = 0.98 * b_ph_pro
        results.append(
            _run_case(
                label=f"segment_near_photon_spin_{spin:.2f}_pro",
                mass=MASS,
                spin=spin,
                impact_b=b_near,
                direction=+1.0,
                r0=r0_pro,
                n_steps=NEAR_PHOTON_STEPS,
                h=NEAR_PHOTON_H,
                case_type="near_photon_sphere_diagnostic_segment",
                enforce_radial_sign=False,
                apply_schwarzschild_limit=False,
            )
        )

        r_ph_retro = photon_sphere_radius_retro(MASS, spin)
        b_ph_retro = photon_impact_parameter(r_ph_retro, MASS, spin, prograde=False)
        r0_retro = r_ph_retro + 0.02
        b_near_retro = 0.98 * b_ph_retro
        results.append(
            _run_case(
                label=f"segment_near_photon_spin_{spin:.2f}_retro",
                mass=MASS,
                spin=spin,
                impact_b=b_near_retro,
                direction=+1.0,
                r0=r0_retro,
                n_steps=NEAR_PHOTON_STEPS,
                h=NEAR_PHOTON_H,
                case_type="near_photon_sphere_diagnostic_segment",
                enforce_radial_sign=False,
                apply_schwarzschild_limit=False,
            )
        )
    return results


def write_csv(records: list[dict[str, Any]], path: Path) -> None:
    fields = [
        "case_id", "spin_a", "M", "E", "b", "direction", "case_type", "r_plus", "r0",
        "lambda_end", "n_steps", "delta_t", "delta_phi", "radial_drift", "min_r", "max_r",
        "min_Delta", "min_R", "max_abs_null_residual", "max_abs_E_residual",
        "max_abs_L_residual", "max_abs_radial_rhs_residual", "all_points_exterior",
        "finite_rhs_all_steps", "finite_solution_all_steps", "t_monotonic_future_pass",
        "radial_rhs_consistency_pass", "null_condition_pass", "constants_consistency_pass",
        "schwarzschild_radial_limit_pass", "radial_sign_consistency_pass",
        "angular_accumulation_finite_pass", "no_endpoint_targeting_pass", "all_checks_pass",
    ]
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for rec in records:
            writer.writerow({k: rec[k] for k in fields})


def write_json(records: list[dict[str, Any]], summary: dict[str, Any], path: Path) -> None:
    payload = {
        "benchmark": "S4-KERR-K10 equatorial geodesic segment audit",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "command": COMMAND,
        "scope_note": (
            "K10 integrates initial-value equatorial Kerr null-geodesic segments only. "
            "No point-to-point shooting, no boundary-value solving, no causal reachability, "
            "and no sprinkled event-pair classification."
        ),
        "segments": records,
        "global_summary": summary,
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def write_md(records: list[dict[str, Any]], summary: dict[str, Any], path: Path) -> None:
    lines = [
        "# S4-KERR-K10 equatorial geodesic segment audit",
        "",
        "K10 integrates initial-value equatorial Kerr null-geodesic segments.",
        "K10 does not do point-to-point shooting.",
        "K10 does not solve boundary-value problems.",
        "K10 does not decide causal reachability.",
        "K10 does not classify sprinkled event pairs.",
        "Delta_phi and Delta_t are diagnostics only, not evidence of endpoint connection.",
        "K10 is the final preflight before a future K11 shooting sandbox.",
        "",
        f"- Total segments: {summary['total_segments']}",
        f"- Passed segments: {summary['passed_segments']}",
        f"- Failed segments: {summary['failed_segments']}",
        f"- Advisory segments: {summary['advisory_segments']}",
        f"- Global undecided pairs (a>0 control accounting): {summary['global_undecided_pairs']}",
        "",
        "## Artifact Set",
        "",
        f"- `{OUT_PREFIX}.csv`",
        f"- `{OUT_PREFIX}.json`",
        f"- `{OUT_PREFIX}.md`",
        f"- `{OUT_PREFIX}.png`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_png(results: list[dict[str, Any]], records: list[dict[str, Any]], path: Path) -> None:
    fig, axs = plt.subplots(2, 2, figsize=(11, 8))
    fig.suptitle("S4-KERR-K10 equatorial geodesic segment audit")

    ax = axs[0, 0]
    for row, raw in zip(records, results):
        if row["case_type"] in {
            "safe_b0_control_segment",
            "non_circular_prograde_safe_segment",
            "non_circular_retrograde_safe_segment",
        }:
            lam = [i * raw["h"] for i in range(len(raw["states"]))]
            rr = [s[1] for s in raw["states"]]
            ax.plot(lam, rr, label=row["case_id"])
    ax.set_xlabel("lambda")
    ax.set_ylabel("r(lambda)")
    ax.set_title("Control/non-circular segments")
    ax.legend(fontsize=6, ncol=2)

    ax = axs[0, 1]
    for row, raw in zip(records, results):
        if "prograde" in row["case_type"] or "retrograde" in row["case_type"]:
            lam = [i * raw["h"] for i in range(len(raw["states"]))]
            ph = [s[2] for s in raw["states"]]
            ax.plot(lam, ph, label=row["case_id"])
    ax.set_xlabel("lambda")
    ax.set_ylabel("phi(lambda)")
    ax.set_title("Angular accumulation diagnostics")

    ax = axs[1, 0]
    x = list(range(len(records)))
    dphi = [r["delta_phi"] for r in records]
    ax.plot(x, dphi, "o-")
    ax.set_xlabel("case index")
    ax.set_ylabel("Delta_phi")
    ax.set_title("Delta_phi by segment")

    ax = axs[1, 1]
    x = list(range(len(records)))
    nulls = [r["max_abs_null_residual"] for r in records]
    ers = [r["max_abs_E_residual"] for r in records]
    lrs = [r["max_abs_L_residual"] for r in records]
    ax.semilogy(x, nulls, "o-", label="null residual")
    ax.semilogy(x, ers, "s-", label="E residual")
    ax.semilogy(x, lrs, "^-", label="L residual")
    ax.set_xlabel("case index")
    ax.set_title("Residual diagnostics by segment")
    ax.legend(fontsize=7)

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def run_audit() -> int:
    results = build_results()

    non_advisory_failures = [r for r in results if (not r["pass"]) and (not r["advisory_only"])]
    if non_advisory_failures:
        print("STOP: non-advisory segment failure detected.")
        for rr in non_advisory_failures:
            print(f"  failed_segment={rr['label']} reason={rr['fail_reason']}")
        return 1

    records = [_segment_record(r) for r in results]
    total_segments = len(records)
    passed_segments = sum(1 for r in records if r["all_checks_pass"])
    failed_segments = total_segments - passed_segments
    advisory_segments = sum(1 for r in results if r["advisory_only"])
    summary = {
        "total_segments": total_segments,
        "passed_segments": passed_segments,
        "failed_segments": failed_segments,
        "advisory_segments": advisory_segments,
        "global_true_relations": 0,
        "global_false_relations": 0,
        "global_undecided_pairs": N_EVENTS * (N_EVENTS - 1) // 2,
        "all_checks_pass": failed_segments == 0,
    }

    csv_path = ARTIFACT_DIR / f"{OUT_PREFIX}.csv"
    json_path = ARTIFACT_DIR / f"{OUT_PREFIX}.json"
    md_path = ARTIFACT_DIR / f"{OUT_PREFIX}.md"
    png_path = ARTIFACT_DIR / f"{OUT_PREFIX}.png"
    write_csv(records, csv_path)
    write_json(records, summary, json_path)
    write_md(records, summary, md_path)
    write_png(results, records, png_path)
    print(f"segments={total_segments} passed={passed_segments} failed={failed_segments}")
    return 0


def main() -> None:
    raise SystemExit(run_audit())


if __name__ == "__main__":
    main()
