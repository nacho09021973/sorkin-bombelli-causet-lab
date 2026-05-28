#!/usr/bin/env python3
"""S4-KERR-K13B-NEAR-PHOTON-WHIRLING-PROBE-001.

Direction-aware near-photon initial-value whirling probe in equatorial Kerr.
This is diagnostic only; it does not perform causal classification.
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
from audit_kerr_k10_equatorial_geodesic_segment_audit_001 import check_segment  # noqa: E402
from audit_kerr_k11_equatorial_shooting_sandbox_001 import ENERGY, MASS, HORIZON_SAFETY  # noqa: E402
from audit_kerr_k9_equatorial_full_rhs_preflight_001 import (  # noqa: E402
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


ARTIFACT_DIR = Path(__file__).resolve().parent
OUT_PREFIX = "kerr_k13b_near_photon_whirling_probe_001_n12_seed1959"
SPINS = (0.0, 0.25, 0.5)
BRANCHES = ("prograde", "retrograde")
DIRECTIONS: tuple[tuple[str, float], ...] = (("outgoing", +1.0), ("ingoing", -1.0))
DELTA_RS = (0.02, 0.05, 0.10)
EPSILONS = (-1.0e-2, -3.0e-3, -1.0e-3, -3.0e-4, 3.0e-4, 1.0e-3, 3.0e-3, 1.0e-2)
LAMBDA_ENDS = (2.4, 5.0, 8.0)
STEP_H = 0.01
N_EVENTS = 12
NEAR_PHOTON_BAND = 0.10
HARD_NULL_TOL = 1.0e-7
HARD_EL_TOL = 1.0e-7
HARD_RADIAL_RHS_TOL = 1.0e-6
R_NEG_TOL = R_MIN_TOL


def _branch_params(spin: float, branch: str) -> tuple[float, float]:
    if branch == "prograde":
        r_ph = photon_sphere_radius_pro(MASS, spin)
        b_ph = photon_impact_parameter(r_ph, MASS, spin, prograde=True)
    else:
        r_ph = photon_sphere_radius_retro(MASS, spin)
        b_ph = photon_impact_parameter(r_ph, MASS, spin, prograde=False)
    return r_ph, b_ph


def integrate_segment(
    *,
    spin: float,
    b: float,
    r0: float,
    direction: float,
    lambda_end: float,
) -> dict[str, Any]:
    r_plus = kerr_horizon_radius(MASS, spin)
    n_steps = int(round(lambda_end / STEP_H))
    state = (0.0, r0, 0.0)
    states = [state]
    rhs_vals: list[tuple[float, float, float]] = []
    stop_reason = "normal"

    if r0 <= r_plus + HORIZON_SAFETY:
        return {"states": states, "rhs_vals": rhs_vals, "stop_reason": "safety_stop"}

    r0_potential = radial_potential(r0, MASS, spin, b, energy=ENERGY)
    if r0_potential < -R_NEG_TOL:
        return {"states": states, "rhs_vals": rhs_vals, "stop_reason": "forbidden_initial_R"}

    for _ in range(n_steps):
        _, r, _ = state
        dval = delta(r, MASS, spin)
        if dval <= 0.0:
            stop_reason = "safety_stop"
            break
        if r <= r_plus + HORIZON_SAFETY:
            stop_reason = "r_crossed_horizon_margin"
            break

        rpot = radial_potential(r, MASS, spin, b, energy=ENERGY)
        if rpot < -R_NEG_TOL:
            stop_reason = "R_negative"
            break

        try:
            rhs = kerr_equatorial_rhs(state, MASS, spin, b, direction, energy=ENERGY)
        except ValueError as exc:
            msg = str(exc)
            if "Delta<=0" in msg:
                stop_reason = "safety_stop"
            elif "Negative R" in msg:
                stop_reason = "R_negative"
            else:
                stop_reason = "nonfinite_rhs"
            break
        if not all(math.isfinite(x) for x in rhs):
            stop_reason = "nonfinite_rhs"
            break
        rhs_vals.append(rhs)

        tdot, rdot, phidot = rhs
        radial_rhs = direction * math.sqrt(max(rpot, 0.0)) / (r * r)
        nres = abs(null_residual(r, tdot, rdot, phidot, MASS, spin))
        e_calc, l_calc = conserved_energy_angular_momentum(r, tdot, phidot, MASS, spin)
        if (
            nres > HARD_NULL_TOL
            or abs(e_calc - ENERGY) > HARD_EL_TOL
            or abs(l_calc - b) > HARD_EL_TOL
            or abs(rdot - radial_rhs) > HARD_RADIAL_RHS_TOL
        ):
            stop_reason = "invariant_failure"
            break

        try:
            nxt = rk4_step(state, STEP_H, MASS, spin, b, direction, energy=ENERGY)
        except ValueError as exc:
            msg = str(exc)
            if "Delta<=0" in msg:
                stop_reason = "safety_stop"
            elif "Negative R" in msg:
                stop_reason = "R_negative"
            else:
                stop_reason = "nonfinite_rhs"
            break
        if not all(math.isfinite(x) for x in nxt):
            stop_reason = "nonfinite_rhs"
            break
        state = nxt
        states.append(state)

    return {"states": states, "rhs_vals": rhs_vals, "stop_reason": stop_reason}


def _time_near_photon(rs: list[float], r_ph: float) -> float:
    if not rs:
        return 0.0
    count = sum(1 for r in rs if abs(r - r_ph) < NEAR_PHOTON_BAND)
    return count / float(len(rs))


def _row(
    *,
    case_id: str,
    spin: float,
    branch: str,
    direction_name: str,
    direction_sign: float,
    b: float,
    b_ph: float,
    epsilon_b: float,
    r_ph: float,
    r0: float,
    delta_r: float,
    lambda_end: float,
    run: dict[str, Any],
) -> dict[str, Any]:
    states = run["states"]
    rhs_vals = run["rhs_vals"]
    stop_reason = run["stop_reason"]
    r_plus = kerr_horizon_radius(MASS, spin)
    rs = [s[1] for s in states]

    checked = check_segment(
        states=states,
        rhs_vals=rhs_vals,
        mass=MASS,
        spin=spin,
        impact_b=b,
        direction=direction_sign,
        energy=ENERGY,
        r_plus=r_plus,
        enforce_radial_sign=False,
        apply_schwarzschild_limit=False,
    )
    checks = checked["checks"]
    metrics = checked["metrics"]

    t0, _, p0 = states[0]
    tf, _, pf = states[-1]
    delta_t = tf - t0
    delta_phi_raw = pf - p0
    abs_delta_phi = abs(delta_phi_raw)
    winding_m_estimate = int(round(delta_phi_raw / (2.0 * math.pi)))
    min_abs = min(abs(r - r_ph) for r in rs) if rs else float("inf")
    r_closest = min(rs, key=lambda r: abs(r - r_ph)) if rs else r0

    weak_winding_pass = abs_delta_phi > (0.5 * math.pi)
    strong_winding_pass = abs_delta_phi > math.pi
    full_turn_probe_pass = abs_delta_phi > (2.0 * math.pi)
    angular_accumulation_finite_pass = math.isfinite(delta_t) and math.isfinite(delta_phi_raw)
    time_near_photon = _time_near_photon(rs, r_ph)

    invariants_ok = (
        checks["all_points_exterior"]
        and checks["finite_rhs_all_steps"]
        and checks["finite_solution_all_steps"]
        and checks["t_monotonic_future_pass"]
        and checks["radial_rhs_consistency_pass"]
        and checks["null_condition_pass"]
        and checks["constants_consistency_pass"]
        and metrics["min_Delta"] > 0.0
        and metrics["min_R"] >= -R_MIN_TOL
        and angular_accumulation_finite_pass
        and metrics["max_abs_null_residual"] <= HARD_NULL_TOL
        and metrics["max_abs_E_residual"] <= HARD_EL_TOL
        and metrics["max_abs_L_residual"] <= HARD_EL_TOL
        and metrics["max_abs_radial_rhs_residual"] <= HARD_RADIAL_RHS_TOL
    )

    unresolved = stop_reason in ("nonfinite_rhs", "invariant_failure")
    advisory_only = stop_reason != "normal"

    no_sprinkling_pair_used = True
    no_global_causal_relations_decided = True
    no_causal_classifier_introduced = True

    all_checks_pass = (
        (not unresolved)
        and (not advisory_only)
        and invariants_ok
        and no_sprinkling_pair_used
        and no_global_causal_relations_decided
        and no_causal_classifier_introduced
    )

    return {
        "case_id": case_id,
        "spin_a": spin,
        "branch": branch,
        "direction": direction_name,
        "M": MASS,
        "E": ENERGY,
        "b": b,
        "b_ph": b_ph,
        "epsilon_b": epsilon_b,
        "r_plus": r_plus,
        "r_ph": r_ph,
        "r0": r0,
        "delta_r": delta_r,
        "lambda_end": lambda_end,
        "lambda_until_stop": (len(states) - 1) * STEP_H,
        "n_steps": len(states) - 1,
        "stop_reason": stop_reason,
        "delta_t": delta_t,
        "delta_phi_raw": delta_phi_raw,
        "abs_delta_phi": abs_delta_phi,
        "winding_m_estimate": winding_m_estimate,
        "r_closest_to_r_ph": r_closest,
        "min_abs_r_minus_r_ph": min_abs,
        "radial_drift": metrics["radial_drift"],
        "time_near_photon": time_near_photon,
        "min_r": metrics["min_r"],
        "max_r": metrics["max_r"],
        "min_Delta": metrics["min_Delta"],
        "min_R": metrics["min_R"],
        "max_abs_null_residual": metrics["max_abs_null_residual"],
        "max_abs_E_residual": metrics["max_abs_E_residual"],
        "max_abs_L_residual": metrics["max_abs_L_residual"],
        "all_points_exterior": checks["all_points_exterior"],
        "finite_rhs_all_steps": checks["finite_rhs_all_steps"],
        "finite_solution_all_steps": checks["finite_solution_all_steps"],
        "t_monotonic_future_pass": checks["t_monotonic_future_pass"],
        "radial_rhs_consistency_pass": checks["radial_rhs_consistency_pass"],
        "null_condition_pass": checks["null_condition_pass"],
        "constants_consistency_pass": checks["constants_consistency_pass"],
        "angular_accumulation_finite_pass": angular_accumulation_finite_pass,
        "weak_winding_pass": weak_winding_pass,
        "strong_winding_pass": strong_winding_pass,
        "full_turn_probe_pass": full_turn_probe_pass,
        "no_sprinkling_pair_used": no_sprinkling_pair_used,
        "no_global_causal_relations_decided": no_global_causal_relations_decided,
        "no_causal_classifier_introduced": no_causal_classifier_introduced,
        "advisory_only": advisory_only,
        "unresolved": unresolved,
        "all_checks_pass": all_checks_pass,
        "_states": states,
    }


def scan_cases() -> tuple[list[dict[str, Any]], int]:
    rows: list[dict[str, Any]] = []
    scanned = 0
    for spin in SPINS:
        for branch in BRANCHES:
            r_ph, b_ph = _branch_params(spin, branch)
            for delta_r in DELTA_RS:
                r0 = r_ph + delta_r
                for eps in EPSILONS:
                    b = b_ph * (1.0 + eps)
                    for direction_name, direction_sign in DIRECTIONS:
                        for lambda_end in LAMBDA_ENDS:
                            scanned += 1
                            run = integrate_segment(
                                spin=spin,
                                b=b,
                                r0=r0,
                                direction=direction_sign,
                                lambda_end=lambda_end,
                            )
                            cid = (
                                f"k13b_a{spin:.2f}_{branch}_{direction_name}"
                                f"_dr{delta_r:.2f}_eps{eps:+.4f}_le{lambda_end:.1f}"
                            )
                            rows.append(
                                _row(
                                    case_id=cid,
                                    spin=spin,
                                    branch=branch,
                                    direction_name=direction_name,
                                    direction_sign=direction_sign,
                                    b=b,
                                    b_ph=b_ph,
                                    epsilon_b=eps,
                                    r_ph=r_ph,
                                    r0=r0,
                                    delta_r=delta_r,
                                    lambda_end=lambda_end,
                                    run=run,
                                )
                            )
    return rows, scanned


def write_csv(rows: list[dict[str, Any]], path: Path) -> None:
    fields = [
        "case_id", "spin_a", "branch", "direction", "M", "E", "b", "b_ph", "epsilon_b",
        "r_plus", "r_ph", "r0", "delta_r", "lambda_end", "lambda_until_stop", "n_steps",
        "stop_reason", "delta_t", "delta_phi_raw", "abs_delta_phi", "winding_m_estimate",
        "r_closest_to_r_ph", "min_abs_r_minus_r_ph", "radial_drift", "time_near_photon",
        "min_r", "max_r", "min_Delta", "min_R", "max_abs_null_residual", "max_abs_E_residual",
        "max_abs_L_residual", "all_points_exterior", "finite_rhs_all_steps",
        "finite_solution_all_steps", "t_monotonic_future_pass", "radial_rhs_consistency_pass",
        "null_condition_pass", "constants_consistency_pass", "angular_accumulation_finite_pass",
        "weak_winding_pass", "strong_winding_pass", "full_turn_probe_pass",
        "no_sprinkling_pair_used", "no_global_causal_relations_decided",
        "no_causal_classifier_introduced", "advisory_only", "unresolved", "all_checks_pass",
    ]
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k) for k in fields})


def write_json(rows: list[dict[str, Any]], scanned: int, path: Path) -> None:
    accepted = [r for r in rows if (not r["advisory_only"]) and (not r["unresolved"])]
    best = max(rows, key=lambda r: r["abs_delta_phi"]) if rows else None
    summary = {
        "total_candidates_scanned": scanned,
        "total_cases_recorded": len(rows),
        "accepted_cases": len(accepted),
        "advisory_cases": sum(1 for r in rows if r["advisory_only"]),
        "unresolved_cases": sum(1 for r in rows if r["unresolved"]),
        "max_abs_delta_phi": (best["abs_delta_phi"] if best else 0.0),
        "best_case_id": (best["case_id"] if best else ""),
        "best_case_spin_a": (best["spin_a"] if best else None),
        "best_case_branch": (best["branch"] if best else ""),
        "best_case_direction": (best["direction"] if best else ""),
        "best_case_epsilon_b": (best["epsilon_b"] if best else None),
        "best_case_delta_r": (best["delta_r"] if best else None),
        "best_case_stop_reason": (best["stop_reason"] if best else ""),
        "best_case_min_abs_r_minus_r_ph": (best["min_abs_r_minus_r_ph"] if best else None),
        "best_case_time_near_photon": (best["time_near_photon"] if best else None),
        "any_weak_winding_pass": any(r["weak_winding_pass"] for r in rows),
        "any_strong_winding_pass": any(r["strong_winding_pass"] for r in rows),
        "any_full_turn_probe_pass": any(r["full_turn_probe_pass"] for r in rows),
        "global_true_relations": 0,
        "global_false_relations": 0,
        "global_undecided_pairs": N_EVENTS * (N_EVENTS - 1) // 2,
        "all_checks_pass": all(r["all_checks_pass"] for r in accepted),
    }
    payload = {
        "benchmark": "S4-KERR-K13b near-photon whirling probe",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "cases": [{k: v for k, v in r.items() if not k.startswith("_")} for r in rows],
        "global_summary": summary,
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def write_md(summary: dict[str, Any], path: Path) -> None:
    lines = [
        "# S4-KERR-K13b near-photon whirling probe",
        "",
        "K13b tests physical near-photon whirling with both ingoing/outgoing directions.",
        "K13b is motivated by the K13 diagnosis: K13 fired outward from outside r_ph, so accepted rays escaped instead of whirling.",
        "K13b does not classify causal reachability.",
        "K13b does not use sprinkling event pairs.",
        "K13b does not implement a production Kerr causal classifier.",
        "Large Delta_phi here is angular accumulation along an initial-value segment, not endpoint reachability.",
        "If no strong winding is found, that is reported as a probe result and not hidden.",
        "",
        f"- Total candidates scanned: {summary['total_candidates_scanned']}",
        f"- Total cases recorded: {summary['total_cases_recorded']}",
        f"- Accepted cases: {summary['accepted_cases']}",
        f"- Advisory cases: {summary['advisory_cases']}",
        f"- Unresolved cases: {summary['unresolved_cases']}",
        f"- Max |Delta_phi|: {summary['max_abs_delta_phi']:.6g}",
        f"- Any weak_winding_pass: {summary['any_weak_winding_pass']}",
        f"- Any strong_winding_pass: {summary['any_strong_winding_pass']}",
        f"- Any full_turn_probe_pass: {summary['any_full_turn_probe_pass']}",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_png(rows: list[dict[str, Any]], path: Path) -> None:
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle("S4-KERR-K13b near-photon whirling probe")

    accepted = [r for r in rows if (not r["advisory_only"]) and (not r["unresolved"])]
    draw = accepted if accepted else rows

    ax = axs[0, 0]
    labels = [f"{r['branch'][:3]}-{r['direction'][:3]}" for r in draw]
    ax.scatter(range(len(draw)), [r["abs_delta_phi"] for r in draw], c="tab:blue", s=14)
    ax.set_xticks(range(len(draw)))
    ax.set_xticklabels(labels, rotation=90, fontsize=6)
    ax.set_title("abs(delta_phi) by direction/branch")
    ax.set_ylabel("abs(delta_phi)")

    top = sorted(draw, key=lambda r: r["abs_delta_phi"], reverse=True)[:4]
    ax = axs[0, 1]
    for r in top:
        rs = [s[1] for s in r["_states"]]
        lam = [i * STEP_H for i in range(len(rs))]
        ax.plot(lam, [rv - r["r_ph"] for rv in rs], label=r["case_id"])
    ax.set_title("r(lambda)-r_ph top winding cases")
    ax.set_xlabel("lambda")
    ax.legend(fontsize=6)

    ax = axs[1, 0]
    for r in top:
        ph = [s[2] for s in r["_states"]]
        lam = [i * STEP_H for i in range(len(ph))]
        ax.plot(lam, ph, label=r["case_id"])
    ax.set_title("phi(lambda) top winding cases")
    ax.set_xlabel("lambda")

    ax = axs[1, 1]
    xs = list(range(len(top)))
    ax.semilogy(xs, [max(r["max_abs_null_residual"], 1e-18) for r in top], "o-", label="null")
    ax.semilogy(xs, [max(r["max_abs_E_residual"], 1e-18) for r in top], "s-", label="E")
    ax.semilogy(xs, [max(r["max_abs_L_residual"], 1e-18) for r in top], "^-", label="L")
    ax.set_xticks(xs)
    ax.set_xticklabels([r["case_id"] for r in top], rotation=90, fontsize=6)
    ax.set_title("Invariant residuals (top accepted)")
    ax.legend(fontsize=7)

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main() -> None:
    rows, scanned = scan_cases()
    csv_path = ARTIFACT_DIR / f"{OUT_PREFIX}.csv"
    json_path = ARTIFACT_DIR / f"{OUT_PREFIX}.json"
    md_path = ARTIFACT_DIR / f"{OUT_PREFIX}.md"
    png_path = ARTIFACT_DIR / f"{OUT_PREFIX}.png"

    write_csv(rows, csv_path)
    write_json(rows, scanned, json_path)
    summary = json.loads(json_path.read_text(encoding="utf-8"))["global_summary"]
    write_md(summary, md_path)
    write_png(rows, png_path)
    print(
        f"scanned={summary['total_candidates_scanned']} recorded={summary['total_cases_recorded']} "
        f"accepted={summary['accepted_cases']} max_abs_delta_phi={summary['max_abs_delta_phi']:.6g}"
    )


if __name__ == "__main__":
    main()
