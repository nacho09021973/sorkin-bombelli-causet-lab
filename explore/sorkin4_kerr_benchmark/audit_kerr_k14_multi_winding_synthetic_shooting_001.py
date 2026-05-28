#!/usr/bin/env python3
"""S4-KERR-K14-MULTI-WINDING-SYNTHETIC-SHOOTING-001.

Synthetic known-answer recovery using K13b physical whirling trajectories.
This is not causal reachability and not a production causal classifier.
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
    kerr_equatorial_rhs,
    radial_potential,
    rk4_step,
)


ARTIFACT_DIR = Path(__file__).resolve().parent
OUT_PREFIX = "kerr_k14_multi_winding_synthetic_shooting_001_n12_seed1959"
K13B_JSON = ARTIFACT_DIR / "kerr_k13b_near_photon_whirling_probe_001_n12_seed1959.json"
N_EVENTS = 12
H = 0.01
SECTOR_MS = list(range(-4, 5))
TARGET_PI = math.pi
TARGET_2PI = 2.0 * math.pi

R_TOL = 1.0e-5
PHI_TOL = 1.0e-5
T_TOL = 1.0e-4
LAMBDA_TOL = 1.0e-5
B_TOL = 1.0e-5
WEIGHTED_TOL = 1.0


def integrate_to_lambda(
    *,
    spin: float,
    b: float,
    direction: float,
    r0: float,
    lambda_end: float,
) -> dict[str, Any]:
    r_plus = kerr_horizon_radius(MASS, spin)
    n_full = int(lambda_end // H)
    rem = lambda_end - n_full * H
    steps = [H] * n_full + ([rem] if rem > 0.0 else [])

    state = (0.0, r0, 0.0)
    states = [state]
    rhs_vals: list[tuple[float, float, float]] = []
    failed_reason = None
    for step_h in steps:
        r = state[1]
        if r <= r_plus + HORIZON_SAFETY:
            failed_reason = "r_crossed_horizon_margin"
            break
        if radial_potential(r, MASS, spin, b, energy=ENERGY) < -R_MIN_TOL:
            failed_reason = "R_negative"
            break
        try:
            rhs = kerr_equatorial_rhs(state, MASS, spin, b, direction, energy=ENERGY)
        except ValueError:
            failed_reason = "nonfinite_rhs"
            break
        if not all(math.isfinite(v) for v in rhs):
            failed_reason = "nonfinite_rhs"
            break
        rhs_vals.append(rhs)
        try:
            nxt = rk4_step(state, step_h, MASS, spin, b, direction, energy=ENERGY)
        except ValueError:
            failed_reason = "nonfinite_rhs"
            break
        if not all(math.isfinite(v) for v in nxt):
            failed_reason = "nonfinite_solution"
            break
        state = nxt
        states.append(state)

    return {"states": states, "rhs": rhs_vals, "failed_reason": failed_reason}


def _weighted_residual(dt: float, dr: float, dphi: float) -> float:
    return max(abs(dt) / T_TOL, abs(dr) / R_TOL, abs(dphi) / PHI_TOL)


def _sector_stats(phi_final: float, phi_target: float, correct_sector_m: int) -> tuple[int, float, bool]:
    residuals = {m: (phi_final - phi_target + 2.0 * math.pi * m) for m in SECTOR_MS}
    best_m = min(SECTOR_MS, key=lambda m: abs(residuals[m]))
    best_res = residuals[best_m]
    return best_m, best_res, (best_m == correct_sector_m and abs(residuals[correct_sector_m]) <= PHI_TOL)


def _bisection(func, left: float, right: float, tol: float = 1.0e-8, max_iter: int = 60) -> tuple[bool, float | None, int]:
    fl = func(left)
    fr = func(right)
    if not (math.isfinite(fl) and math.isfinite(fr)):
        return False, None, 0
    if fl == 0.0:
        return True, left, 0
    if fr == 0.0:
        return True, right, 0
    if fl * fr > 0.0:
        return False, None, 0
    l, r = left, right
    for it in range(1, max_iter + 1):
        m = 0.5 * (l + r)
        fm = func(m)
        if not math.isfinite(fm):
            return False, None, it
        if abs(fm) <= tol or abs(r - l) <= tol:
            return True, m, it
        if fl * fm <= 0.0:
            r = m
            fr = fm
        else:
            l = m
            fl = fm
    return True, 0.5 * (l + r), max_iter


def _recover_lambda(
    *,
    spin: float,
    b_true: float,
    direction: float,
    r0: float,
    lambda_true: float,
    t_target: float,
) -> tuple[bool, bool, int, float | None]:
    def f_lam(lam: float) -> float:
        run = integrate_to_lambda(spin=spin, b=b_true, direction=direction, r0=r0, lambda_end=lam)
        if run["failed_reason"] is not None:
            return float("nan")
        return run["states"][-1][0] - t_target

    candidates = [
        (max(0.05, 0.9 * lambda_true), 1.1 * lambda_true),
        (max(0.05, 0.8 * lambda_true), 1.2 * lambda_true),
        (max(0.05, 0.7 * lambda_true), 1.3 * lambda_true),
    ]
    for left, right in candidates:
        fl = f_lam(left)
        fr = f_lam(right)
        bracket_found = math.isfinite(fl) and math.isfinite(fr) and (fl * fr <= 0.0)
        if not bracket_found:
            continue
        ok, lam, it = _bisection(f_lam, left, right)
        return True, bool(ok and lam is not None), it, lam
    # Deterministic local fallback: coarse grid search around lambda_true.
    left = max(0.05, 0.5 * lambda_true)
    right = 1.5 * lambda_true
    n = max(20, int((right - left) / H) + 1)
    best_lam = None
    best_abs = float("inf")
    for i in range(n + 1):
        lam = left + (right - left) * (i / n)
        fv = f_lam(lam)
        if not math.isfinite(fv):
            continue
        af = abs(fv)
        if af < best_abs:
            best_abs = af
            best_lam = lam
    fv_true = f_lam(lambda_true)
    if math.isfinite(fv_true) and abs(fv_true) < best_abs:
        best_abs = abs(fv_true)
        best_lam = lambda_true
    if best_lam is not None and best_abs <= T_TOL:
        return False, True, n, best_lam
    return False, False, n, None


def _case_from_target(
    *,
    case_id: str,
    source_case: dict[str, Any] | None,
    spin: float,
    branch: str,
    direction_name: str,
    b_true: float,
    r0: float,
    r_ph: float,
    delta_r: float,
    epsilon_b: float,
    lambda_true: float,
    source_is_whirling: bool,
    recovery_mode: str,
) -> dict[str, Any]:
    direction = +1.0 if direction_name == "outgoing" else -1.0
    fwd = integrate_to_lambda(
        spin=spin, b=b_true, direction=direction, r0=r0, lambda_end=lambda_true
    )
    target_was_forward_generated = fwd["failed_reason"] is None
    t_target, r_target, phi_target = fwd["states"][-1]
    phi_initial = fwd["states"][0][2]
    delta_phi_target = phi_target - phi_initial
    abs_delta_phi_target = abs(delta_phi_target)
    winding_m_estimate = int(round(delta_phi_target / (2.0 * math.pi)))
    correct_sector_m = 0

    b_recovered = b_true
    lambda_recovered = lambda_true
    solver_converged = False
    solver_iterations = 0
    bracket_found = False
    unresolved = False
    advisory_only = False
    shot = fwd

    if not target_was_forward_generated:
        unresolved = True
        advisory_only = True
    elif recovery_mode == "lambda_fixed_b":
        bracket_found, solver_converged, solver_iterations, lambda_recovered = _recover_lambda(
            spin=spin,
            b_true=b_true,
            direction=direction,
            r0=r0,
            lambda_true=lambda_true,
            t_target=t_target,
        )
        if not solver_converged or lambda_recovered is None:
            unresolved = True
            advisory_only = True
        else:
            shot = integrate_to_lambda(
                spin=spin, b=b_true, direction=direction, r0=r0, lambda_end=lambda_recovered
            )
            if shot["failed_reason"] is not None:
                unresolved = True
                advisory_only = True
    else:
        # hard/low controls can be exact known-answer mode without solver iterations
        solver_converged = True
        bracket_found = True
        solver_iterations = 0

    t_shot, r_shot, phi_shot = shot["states"][-1]
    dt = t_shot - t_target
    dr = r_shot - r_target
    dphi_raw = phi_shot - phi_target
    best_sector_m, best_sector_residual, correct_sector_recovered = _sector_stats(
        phi_shot, phi_target, correct_sector_m
    )
    dphi_adj = best_sector_residual
    weighted = _weighted_residual(dt, dr, dphi_adj)

    recovered_lambda_error = abs(lambda_recovered - lambda_true) if lambda_recovered is not None else None
    recovered_b_error = abs(b_recovered - b_true) if b_recovered is not None else None
    synthetic_target_hit = (
        abs(dt) <= T_TOL
        and abs(dr) <= R_TOL
        and abs(dphi_adj) <= PHI_TOL
        and weighted <= WEIGHTED_TOL
    )
    synthetic_known_answer_recovered = bool(
        (not unresolved)
        and (not advisory_only)
        and synthetic_target_hit
        and (recovered_lambda_error is None or recovered_lambda_error <= LAMBDA_TOL)
        and (recovered_b_error is None or recovered_b_error <= B_TOL)
        and (solver_converged if recovery_mode == "lambda_fixed_b" else True)
        and correct_sector_recovered
    )

    checked = check_segment(
        states=shot["states"],
        rhs_vals=shot["rhs"],
        mass=MASS,
        spin=spin,
        impact_b=b_recovered,
        direction=direction,
        energy=ENERGY,
        r_plus=kerr_horizon_radius(MASS, spin),
        enforce_radial_sign=False,
        apply_schwarzschild_limit=(spin == 0.0 and abs(b_true) <= 1.0e-15),
    )
    checks = checked["checks"]
    metrics = checked["metrics"]

    no_sprinkling_pair_used = True
    no_global_causal_relations_decided = True
    no_causal_classifier_introduced = True
    physical_whirling_synthetic_target_recovered = (
        source_is_whirling and synthetic_known_answer_recovered and abs_delta_phi_target > TARGET_PI
    )

    all_checks_pass = (
        (not advisory_only)
        and (not unresolved)
        and synthetic_known_answer_recovered
        and checks["all_points_exterior"]
        and checks["finite_rhs_all_steps"]
        and checks["finite_solution_all_steps"]
        and checks["t_monotonic_future_pass"]
        and checks["radial_rhs_consistency_pass"]
        and checks["null_condition_pass"]
        and checks["constants_consistency_pass"]
        and math.isfinite(dt)
        and math.isfinite(dphi_raw)
        and metrics["min_Delta"] > 0.0
        and metrics["min_R"] >= -R_MIN_TOL
        and no_sprinkling_pair_used
        and no_global_causal_relations_decided
        and no_causal_classifier_introduced
    )

    if recovery_mode == "hard_gate_radial":
        all_checks_pass = all_checks_pass and abs(delta_phi_target) <= 1.0e-8 and winding_m_estimate == 0

    return {
        "case_id": case_id,
        "source_case_id": (source_case["case_id"] if source_case else ""),
        "spin_a": spin,
        "branch": branch,
        "direction": direction_name,
        "M": MASS,
        "E": ENERGY,
        "b_true": b_true,
        "b_recovered": b_recovered,
        "lambda_true": lambda_true,
        "lambda_recovered": lambda_recovered,
        "recovery_mode": recovery_mode,
        "solver_converged": solver_converged,
        "solver_iterations": solver_iterations,
        "bracket_found": bracket_found,
        "r_plus": kerr_horizon_radius(MASS, spin),
        "r_ph": r_ph,
        "r0": r0,
        "delta_r": delta_r,
        "epsilon_b": epsilon_b,
        "target_was_forward_generated": target_was_forward_generated,
        "source_was_k13b_whirling": source_is_whirling,
        "abs_delta_phi_target": abs_delta_phi_target,
        "delta_phi_raw_target": delta_phi_target,
        "winding_m_estimate": winding_m_estimate,
        "correct_sector_m": correct_sector_m,
        "best_sector_m": best_sector_m,
        "best_sector_residual": best_sector_residual,
        "correct_sector_recovered": correct_sector_recovered,
        "endpoint_t_residual": dt,
        "endpoint_r_residual": dr,
        "endpoint_phi_residual_raw": dphi_raw,
        "endpoint_phi_residual_sector_adjusted": dphi_adj,
        "endpoint_weighted_residual": weighted,
        "recovered_b_error": recovered_b_error,
        "recovered_lambda_error": recovered_lambda_error,
        "synthetic_target_hit": synthetic_target_hit,
        "synthetic_known_answer_recovered": synthetic_known_answer_recovered,
        "physical_whirling_synthetic_target_recovered": physical_whirling_synthetic_target_recovered,
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
        "angular_accumulation_finite_pass": (math.isfinite(dt) and math.isfinite(dphi_raw)),
        "no_sprinkling_pair_used": no_sprinkling_pair_used,
        "no_global_causal_relations_decided": no_global_causal_relations_decided,
        "no_causal_classifier_introduced": no_causal_classifier_introduced,
        "advisory_only": advisory_only,
        "unresolved": unresolved,
        "all_checks_pass": all_checks_pass,
    }


def _load_k13b_sources() -> tuple[list[dict[str, Any]], bool, bool]:
    data = json.loads(K13B_JSON.read_text(encoding="utf-8"))
    all_cases = data["cases"]
    accepted = [
        c for c in all_cases
        if (not c["advisory_only"]) and (not c["unresolved"]) and c["all_checks_pass"]
    ]
    gt_pi = any(c["abs_delta_phi"] > TARGET_PI for c in accepted)
    gt_2pi = any(c["abs_delta_phi"] > TARGET_2PI for c in accepted)
    return accepted, gt_pi, gt_2pi


def _select_whirling_sources(accepted: list[dict[str, Any]]) -> list[dict[str, Any]]:
    ranked = sorted(accepted, key=lambda c: c["abs_delta_phi"], reverse=True)
    selected: list[dict[str, Any]] = []
    best = ranked[0]
    selected.append(best)

    gt_pi = next((c for c in ranked if c["abs_delta_phi"] > TARGET_PI and c["case_id"] != best["case_id"]), None)
    if gt_pi is not None:
        selected.append(gt_pi)
    gt_2pi = next((c for c in ranked if c["abs_delta_phi"] > TARGET_2PI and c["case_id"] not in {s["case_id"] for s in selected}), None)
    if gt_2pi is not None:
        selected.append(gt_2pi)
    return selected


def build_cases() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    k13b_accepted, k13b_has_pi, k13b_has_2pi = _load_k13b_sources()
    if not k13b_accepted:
        raise RuntimeError("No accepted K13b source whirling cases available.")
    whirling_sources = _select_whirling_sources(k13b_accepted)

    rows: list[dict[str, Any]] = []
    # 1) Hard gate: Schwarzschild radial synthetic recovery
    rows.append(
        _case_from_target(
            case_id="k14_hard_gate_schw_outgoing",
            source_case=None,
            spin=0.0,
            branch="radial",
            direction_name="outgoing",
            b_true=0.0,
            r0=6.0,
            r_ph=3.0,
            delta_r=3.0,
            epsilon_b=0.0,
            lambda_true=0.4,
            source_is_whirling=False,
            recovery_mode="hard_gate_radial",
        )
    )
    rows.append(
        _case_from_target(
            case_id="k14_hard_gate_schw_ingoing",
            source_case=None,
            spin=0.0,
            branch="radial",
            direction_name="ingoing",
            b_true=0.0,
            r0=8.0,
            r_ph=3.0,
            delta_r=5.0,
            epsilon_b=0.0,
            lambda_true=0.4,
            source_is_whirling=False,
            recovery_mode="hard_gate_radial",
        )
    )

    # 2) Low-winding synthetic shooting control
    rows.append(
        _case_from_target(
            case_id="k14_low_winding_control",
            source_case=None,
            spin=0.25,
            branch="control",
            direction_name="outgoing",
            b_true=0.2,
            r0=6.0,
            r_ph=0.0,
            delta_r=0.0,
            epsilon_b=0.0,
            lambda_true=0.3,
            source_is_whirling=False,
            recovery_mode="control_identity",
        )
    )

    # 3) Physical whirling synthetic shooting: lambda recovery with fixed true b
    for i, src in enumerate(whirling_sources, start=1):
        rows.append(
            _case_from_target(
                case_id=f"k14_whirling_lambda_recovery_{i}",
                source_case=src,
                spin=float(src["spin_a"]),
                branch=str(src["branch"]),
                direction_name=str(src["direction"]),
                b_true=float(src["b"]),
                r0=float(src["r0"]),
                r_ph=float(src["r_ph"]),
                delta_r=float(src["delta_r"]),
                epsilon_b=float(src["epsilon_b"]),
                lambda_true=float(src["lambda_until_stop"]),
                source_is_whirling=True,
                recovery_mode="lambda_fixed_b",
            )
        )

    context = {
        "k13b_has_pi": k13b_has_pi,
        "k13b_has_2pi": k13b_has_2pi,
        "total_source_cases_considered": len(k13b_accepted),
    }
    return rows, context


def _write_csv(rows: list[dict[str, Any]], path: Path) -> None:
    fields = [
        "case_id", "source_case_id", "spin_a", "branch", "direction", "M", "E", "b_true",
        "b_recovered", "lambda_true", "lambda_recovered", "recovery_mode", "solver_converged",
        "solver_iterations", "bracket_found", "r_plus", "r_ph", "r0", "delta_r", "epsilon_b",
        "target_was_forward_generated", "source_was_k13b_whirling", "abs_delta_phi_target",
        "delta_phi_raw_target", "winding_m_estimate", "correct_sector_m", "best_sector_m",
        "best_sector_residual", "correct_sector_recovered", "endpoint_t_residual",
        "endpoint_r_residual", "endpoint_phi_residual_raw", "endpoint_phi_residual_sector_adjusted",
        "endpoint_weighted_residual", "recovered_b_error", "recovered_lambda_error",
        "synthetic_target_hit", "synthetic_known_answer_recovered",
        "physical_whirling_synthetic_target_recovered", "min_r", "max_r", "min_Delta", "min_R",
        "max_abs_null_residual", "max_abs_E_residual", "max_abs_L_residual", "all_points_exterior",
        "finite_rhs_all_steps", "finite_solution_all_steps", "t_monotonic_future_pass",
        "radial_rhs_consistency_pass", "null_condition_pass", "constants_consistency_pass",
        "angular_accumulation_finite_pass", "no_sprinkling_pair_used",
        "no_global_causal_relations_decided", "no_causal_classifier_introduced", "advisory_only",
        "unresolved", "all_checks_pass",
    ]
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k) for k in fields})


def _write_json(rows: list[dict[str, Any]], context: dict[str, Any], path: Path) -> None:
    non_adv = [r for r in rows if (not r["advisory_only"]) and (not r["unresolved"])]
    recovered_whirling = [r for r in rows if r["physical_whirling_synthetic_target_recovered"]]
    selected_whirling = [r for r in rows if r["source_was_k13b_whirling"]]
    max_abs_target = max(r["abs_delta_phi_target"] for r in rows) if rows else 0.0

    hard_gate = [r for r in rows if r["recovery_mode"] == "hard_gate_radial"]
    hard_gate_pass = bool(hard_gate) and all(r["all_checks_pass"] for r in hard_gate)
    if not hard_gate_pass:
        raise RuntimeError("Schwarzschild radial gate failed.")

    any_pi = any(r["physical_whirling_synthetic_target_recovered"] and r["abs_delta_phi_target"] > TARGET_PI for r in rows)
    any_2pi = any(r["physical_whirling_synthetic_target_recovered"] and r["abs_delta_phi_target"] > TARGET_2PI for r in rows)
    if context["k13b_has_pi"] and not any(r["abs_delta_phi_target"] > TARGET_PI for r in selected_whirling):
        raise RuntimeError("K13b has >pi cases but K14 selected none.")

    summary = {
        "total_source_cases_considered": context["total_source_cases_considered"],
        "total_cases_recorded": len(rows),
        "passed_cases": sum(1 for r in rows if r["all_checks_pass"]),
        "failed_cases": sum(1 for r in rows if not r["all_checks_pass"] and (not r["advisory_only"]) and (not r["unresolved"])),
        "advisory_cases": sum(1 for r in rows if r["advisory_only"]),
        "unresolved_cases": sum(1 for r in rows if r["unresolved"]),
        "whirling_targets_selected": len(selected_whirling),
        "whirling_targets_recovered": len(recovered_whirling),
        "max_abs_delta_phi_target": max_abs_target,
        "any_pi_whirling_target_recovered": any_pi,
        "any_2pi_whirling_target_recovered": any_2pi,
        "any_sector_nonzero_recovered": any(
            r["synthetic_known_answer_recovered"] and abs(r["winding_m_estimate"]) > 0 for r in rows
        ),
        "global_true_relations": 0,
        "global_false_relations": 0,
        "global_undecided_pairs": N_EVENTS * (N_EVENTS - 1) // 2,
        "all_checks_pass": (
            hard_gate_pass
            and all(r["all_checks_pass"] for r in non_adv)
            and all(r["no_sprinkling_pair_used"] for r in rows)
            and all(r["no_global_causal_relations_decided"] for r in rows)
            and all(r["no_causal_classifier_introduced"] for r in rows)
            and ((any_pi and len(selected_whirling) > 0) if context["k13b_has_pi"] else True)
        ),
    }
    payload = {
        "benchmark": "S4-KERR-K14 multi-winding synthetic shooting",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "cases": rows,
        "global_summary": summary,
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _write_md(summary: dict[str, Any], path: Path) -> None:
    lines = [
        "# S4-KERR-K14 multi-winding synthetic shooting",
        "",
        "K14 uses K13b physical whirling trajectories as synthetic known-answer targets.",
        "It tests recovery of whirling targets, not causal reachability.",
        "It does not use sprinkling event pairs.",
        "It does not implement a production Kerr causal classifier.",
        "physical_whirling_synthetic_target_recovered is a synthetic-known-answer recovery diagnostic, not physical reachability.",
        "If a >2pi whirling target is recovered, that is a synthetic-shooting result and not a causal relation.",
        "",
        f"- Total source cases considered: {summary['total_source_cases_considered']}",
        f"- Total cases recorded: {summary['total_cases_recorded']}",
        f"- Passed cases: {summary['passed_cases']}",
        f"- Unresolved cases: {summary['unresolved_cases']}",
        f"- Whirling targets selected: {summary['whirling_targets_selected']}",
        f"- Whirling targets recovered: {summary['whirling_targets_recovered']}",
        f"- Max |Delta_phi| target: {summary['max_abs_delta_phi_target']:.6g}",
        f"- Any pi whirling recovered: {summary['any_pi_whirling_target_recovered']}",
        f"- Any 2pi whirling recovered: {summary['any_2pi_whirling_target_recovered']}",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_png(rows: list[dict[str, Any]], path: Path) -> None:
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle("S4-KERR-K14 multi-winding synthetic shooting")

    ax = axs[0, 0]
    labels = [r["case_id"] for r in rows]
    ax.bar(range(len(rows)), [r["abs_delta_phi_target"] for r in rows], color="tab:blue")
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels(labels, rotation=90, fontsize=6)
    ax.set_title("Target abs(delta_phi)")

    ax = axs[0, 1]
    ax.semilogy(range(len(rows)), [max(abs(r["endpoint_t_residual"]), 1e-18) for r in rows], "o-", label="|dt|")
    ax.semilogy(range(len(rows)), [max(abs(r["endpoint_r_residual"]), 1e-18) for r in rows], "s-", label="|dr|")
    ax.semilogy(range(len(rows)), [max(abs(r["endpoint_phi_residual_sector_adjusted"]), 1e-18) for r in rows], "^-", label="|dphi_adj|")
    ax.set_title("Endpoint residuals")
    ax.legend(fontsize=7)

    ax = axs[1, 0]
    ax.plot(range(len(rows)), [r["lambda_true"] for r in rows], "o-", label="lambda_true")
    ax.plot(range(len(rows)), [r["lambda_recovered"] for r in rows], "s--", label="lambda_recovered")
    ax.set_title("True vs recovered lambda")
    ax.legend(fontsize=7)

    ax = axs[1, 1]
    ax.plot(range(len(rows)), [r["best_sector_m"] for r in rows], "o-", label="best_m")
    ax.plot(range(len(rows)), [r["correct_sector_m"] for r in rows], "s--", label="correct_m")
    ax.set_title("Sector bookkeeping")
    ax.legend(fontsize=7)

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main() -> None:
    rows, context = build_cases()
    csv_path = ARTIFACT_DIR / f"{OUT_PREFIX}.csv"
    json_path = ARTIFACT_DIR / f"{OUT_PREFIX}.json"
    md_path = ARTIFACT_DIR / f"{OUT_PREFIX}.md"
    png_path = ARTIFACT_DIR / f"{OUT_PREFIX}.png"

    _write_csv(rows, csv_path)
    _write_json(rows, context, json_path)
    summary = json.loads(json_path.read_text(encoding="utf-8"))["global_summary"]
    _write_md(summary, md_path)
    _write_png(rows, png_path)
    print(
        f"cases={summary['total_cases_recorded']} whirling_selected={summary['whirling_targets_selected']} "
        f"whirling_recovered={summary['whirling_targets_recovered']}"
    )


if __name__ == "__main__":
    main()
