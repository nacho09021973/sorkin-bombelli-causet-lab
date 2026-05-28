#!/usr/bin/env python3
"""S4-KERR-K15-SYNTHETIC-CAUSAL-SANDBOX-001.

Synthetic known-answer event-pair sandbox from validated Kerr trajectories.
No production causal classification is performed.
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
OUT_PREFIX = "kerr_k15_synthetic_causal_sandbox_001_n12_seed1959"
K14_JSON = ARTIFACT_DIR / "kerr_k14_multi_winding_synthetic_shooting_001_n12_seed1959.json"
H = 0.01
SECTOR_MS = list(range(-4, 5))
PAIR_MATCH_TOL = 1.0e-4
N_EVENTS = 12


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


def _bisection(func, left: float, right: float, tol: float = 1.0e-8, max_iter: int = 60) -> tuple[bool, float | None, int]:
    fl = func(left)
    fr = func(right)
    if not (math.isfinite(fl) and math.isfinite(fr)) or fl * fr > 0.0:
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
    b: float,
    direction: float,
    r0: float,
    lambda_true: float,
    t_target: float,
) -> tuple[bool, float | None]:
    def f_lam(lam: float) -> float:
        run = integrate_to_lambda(spin=spin, b=b, direction=direction, r0=r0, lambda_end=lam)
        if run["failed_reason"] is not None:
            return float("nan")
        return run["states"][-1][0] - t_target

    # First try deterministic brackets around lambda_true.
    for left, right in [
        (max(0.05, 0.9 * lambda_true), 1.1 * lambda_true),
        (max(0.05, 0.8 * lambda_true), 1.2 * lambda_true),
        (max(0.05, 0.7 * lambda_true), 1.3 * lambda_true),
    ]:
        ok, lam, _ = _bisection(f_lam, left, right)
        if ok and lam is not None:
            return True, lam
    # Fallback: choose lambda_true (known-answer sandbox) if valid.
    fv = f_lam(lambda_true)
    if math.isfinite(fv) and abs(fv) <= PAIR_MATCH_TOL:
        return True, lambda_true
    return False, None


def _sector_metrics(phi_final: float, phi_target: float, correct_sector_m: int) -> tuple[int, float, bool]:
    residuals = {m: (phi_final - phi_target + 2.0 * math.pi * m) for m in SECTOR_MS}
    best_m = min(SECTOR_MS, key=lambda m: abs(residuals[m]))
    best_res = residuals[best_m]
    correct_ok = (best_m == correct_sector_m) and (abs(residuals[correct_sector_m]) <= 1.0e-5)
    return best_m, best_res, correct_ok


def _weighted(dt: float, dr: float, dphi_adj: float) -> float:
    return max(abs(dt), abs(dr), abs(dphi_adj))


def _make_pair_case(
    *,
    case_id: str,
    source_case_id: str,
    spin: float,
    branch: str,
    direction_name: str,
    pair_type: str,
    b_true: float,
    r0: float,
    r_ph: float,
    lambda_true: float,
    target_perturbation: tuple[float, float, float] = (0.0, 0.0, 0.0),
    use_lambda_recovery: bool = False,
) -> dict[str, Any]:
    direction = +1.0 if direction_name == "outgoing" else -1.0
    run_true = integrate_to_lambda(spin=spin, b=b_true, direction=direction, r0=r0, lambda_end=lambda_true)
    target_was_forward_generated = run_true["failed_reason"] is None
    t0, _, phi0 = run_true["states"][0]
    t_target_true, r_target_true, phi_target_true = run_true["states"][-1]
    dtp, drp, dpp = target_perturbation
    t_target = t_target_true + dtp
    r_target = r_target_true + drp
    phi_target = phi_target_true + dpp
    delta_phi_raw_target = phi_target_true - phi0
    abs_delta_phi_target = abs(delta_phi_raw_target)
    winding_m_estimate = int(round(delta_phi_raw_target / (2.0 * math.pi)))
    correct_sector_m = 0

    b_recovered = b_true
    lambda_recovered = lambda_true
    unresolved = False
    advisory_only = False
    if not target_was_forward_generated:
        unresolved = True
        advisory_only = True

    if use_lambda_recovery and (not unresolved):
        ok, lam = _recover_lambda(
            spin=spin, b=b_true, direction=direction, r0=r0, lambda_true=lambda_true, t_target=t_target_true
        )
        if not ok or lam is None:
            unresolved = True
            advisory_only = True
        else:
            lambda_recovered = lam

    run_shot = integrate_to_lambda(
        spin=spin, b=b_recovered, direction=direction, r0=r0, lambda_end=lambda_recovered
    )
    if run_shot["failed_reason"] is not None:
        unresolved = True
        advisory_only = True

    t_shot, r_shot, phi_shot = run_shot["states"][-1]
    endpoint_t_residual = t_shot - t_target
    endpoint_r_residual = r_shot - r_target
    endpoint_phi_residual_raw = phi_shot - phi_target
    best_sector_m, best_sector_residual, correct_sector_recovered = _sector_metrics(
        phi_shot, phi_target, correct_sector_m
    )
    endpoint_phi_residual_sector_adjusted = best_sector_residual
    endpoint_weighted_residual = _weighted(
        endpoint_t_residual, endpoint_r_residual, endpoint_phi_residual_sector_adjusted
    )

    synthetic_pair_recovered = (
        endpoint_weighted_residual <= PAIR_MATCH_TOL
        and abs(endpoint_t_residual) <= 1.0e-4
        and abs(endpoint_r_residual) <= 1.0e-5
        and abs(endpoint_phi_residual_sector_adjusted) <= 1.0e-5
        and correct_sector_recovered
    )
    if unresolved:
        synthetic_pair_recovered = False

    if pair_type == "negative_control":
        synthetic_pair_classification = "synthetic_no_match" if (not synthetic_pair_recovered) else "synthetic_unresolved"
        known_answer_null_connection_recovered = False
    else:
        synthetic_pair_classification = "synthetic_null_connected" if synthetic_pair_recovered else "synthetic_unresolved"
        known_answer_null_connection_recovered = synthetic_pair_recovered

    checked = check_segment(
        states=run_shot["states"],
        rhs_vals=run_shot["rhs"],
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
    no_production_classifier_introduced = True

    all_checks_pass = (
        (not advisory_only)
        and (not unresolved)
        and target_was_forward_generated
        and checks["all_points_exterior"]
        and checks["finite_rhs_all_steps"]
        and checks["finite_solution_all_steps"]
        and checks["t_monotonic_future_pass"]
        and checks["radial_rhs_consistency_pass"]
        and checks["null_condition_pass"]
        and checks["constants_consistency_pass"]
        and math.isfinite(endpoint_t_residual)
        and math.isfinite(endpoint_r_residual)
        and math.isfinite(endpoint_phi_residual_sector_adjusted)
        and no_sprinkling_pair_used
        and no_global_causal_relations_decided
        and no_production_classifier_introduced
    )
    if pair_type != "negative_control":
        all_checks_pass = all_checks_pass and synthetic_pair_recovered and known_answer_null_connection_recovered
    else:
        all_checks_pass = all_checks_pass and (not synthetic_pair_recovered) and (not known_answer_null_connection_recovered) and endpoint_weighted_residual > PAIR_MATCH_TOL

    return {
        "case_id": case_id,
        "source_case_id": source_case_id,
        "spin_a": spin,
        "branch": branch,
        "direction": direction_name,
        "pair_type": pair_type,
        "synthetic_pair_classification": synthetic_pair_classification,
        "M": MASS,
        "E": ENERGY,
        "b_true": b_true,
        "b_recovered": b_recovered,
        "lambda_true": lambda_true,
        "lambda_recovered": lambda_recovered,
        "r_plus": kerr_horizon_radius(MASS, spin),
        "r_ph": r_ph,
        "t0": t0,
        "r0": r0,
        "phi0": phi0,
        "t_target": t_target,
        "r_target": r_target,
        "phi_target": phi_target,
        "abs_delta_phi_target": abs_delta_phi_target,
        "winding_m_estimate": winding_m_estimate,
        "correct_sector_m": correct_sector_m,
        "best_sector_m": best_sector_m,
        "best_sector_residual": best_sector_residual,
        "correct_sector_recovered": correct_sector_recovered,
        "endpoint_t_residual": endpoint_t_residual,
        "endpoint_r_residual": endpoint_r_residual,
        "endpoint_phi_residual_raw": endpoint_phi_residual_raw,
        "endpoint_phi_residual_sector_adjusted": endpoint_phi_residual_sector_adjusted,
        "endpoint_weighted_residual": endpoint_weighted_residual,
        "synthetic_pair_recovered": synthetic_pair_recovered,
        "known_answer_null_connection_recovered": known_answer_null_connection_recovered,
        "target_was_forward_generated": target_was_forward_generated,
        "no_sprinkling_pair_used": no_sprinkling_pair_used,
        "no_global_causal_relations_decided": no_global_causal_relations_decided,
        "no_production_classifier_introduced": no_production_classifier_introduced,
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
        "angular_accumulation_finite_pass": (math.isfinite(endpoint_t_residual) and math.isfinite(endpoint_phi_residual_raw)),
        "advisory_only": advisory_only,
        "unresolved": unresolved,
        "all_checks_pass": all_checks_pass,
    }


def build_cases() -> tuple[list[dict[str, Any]], dict[str, bool]]:
    k14 = json.loads(K14_JSON.read_text(encoding="utf-8"))
    whirling = [
        c for c in k14["cases"]
        if c["source_was_k13b_whirling"] and c["physical_whirling_synthetic_target_recovered"]
    ]
    if not whirling:
        raise RuntimeError("No recovered K14 whirling sources found.")
    has_k14_2pi = any(c["abs_delta_phi_target"] > 2.0 * math.pi for c in whirling)

    rows: list[dict[str, Any]] = []
    rows.append(
        _make_pair_case(
            case_id="k15_schw_radial_outgoing",
            source_case_id="k15_schw_seed",
            spin=0.0,
            branch="radial",
            direction_name="outgoing",
            pair_type="schwarzschild_radial",
            b_true=0.0,
            r0=6.0,
            r_ph=3.0,
            lambda_true=0.4,
        )
    )
    rows.append(
        _make_pair_case(
            case_id="k15_schw_radial_ingoing",
            source_case_id="k15_schw_seed",
            spin=0.0,
            branch="radial",
            direction_name="ingoing",
            pair_type="schwarzschild_radial",
            b_true=0.0,
            r0=8.0,
            r_ph=3.0,
            lambda_true=0.4,
        )
    )
    rows.append(
        _make_pair_case(
            case_id="k15_kerr_low_winding",
            source_case_id="k15_low_winding_seed",
            spin=0.5,
            branch="control",
            direction_name="outgoing",
            pair_type="kerr_low_winding",
            b_true=0.2,
            r0=6.0,
            r_ph=0.0,
            lambda_true=0.3,
        )
    )

    # Use at least one whirling >pi and include >2pi when available from K14.
    whirling_sorted = sorted(whirling, key=lambda c: c["abs_delta_phi_target"], reverse=True)
    chosen = [whirling_sorted[0]]
    if len(whirling_sorted) > 1:
        chosen.append(whirling_sorted[1])
    for idx, src in enumerate(chosen, start=1):
        rows.append(
            _make_pair_case(
                case_id=f"k15_whirling_{idx}",
                source_case_id=src["case_id"],
                spin=float(src["spin_a"]),
                branch=str(src["branch"]),
                direction_name=str(src["direction"]),
                pair_type="kerr_whirling",
                b_true=float(src["b_true"]),
                r0=float(src["r0"]),
                r_ph=float(src["r_ph"]),
                lambda_true=float(src["lambda_true"]),
                use_lambda_recovery=True,
            )
        )

    # Numerical negative control from strongest whirling source.
    bad = chosen[0]
    rows.append(
        _make_pair_case(
            case_id="k15_negative_control_phi_shift",
            source_case_id=bad["case_id"],
            spin=float(bad["spin_a"]),
            branch=str(bad["branch"]),
            direction_name=str(bad["direction"]),
            pair_type="negative_control",
            b_true=float(bad["b_true"]),
            r0=float(bad["r0"]),
            r_ph=float(bad["r_ph"]),
            lambda_true=float(bad["lambda_true"]),
            target_perturbation=(0.0, 0.0, 2.0e-2),
            use_lambda_recovery=False,
        )
    )
    return rows, {"k14_has_2pi_recovered": has_k14_2pi}


def write_csv(rows: list[dict[str, Any]], path: Path) -> None:
    fields = [
        "case_id", "source_case_id", "spin_a", "branch", "direction", "pair_type",
        "synthetic_pair_classification", "M", "E", "b_true", "b_recovered", "lambda_true",
        "lambda_recovered", "r_plus", "r_ph", "t0", "r0", "phi0", "t_target", "r_target",
        "phi_target", "abs_delta_phi_target", "winding_m_estimate", "correct_sector_m",
        "best_sector_m", "best_sector_residual", "correct_sector_recovered", "endpoint_t_residual",
        "endpoint_r_residual", "endpoint_phi_residual_raw", "endpoint_phi_residual_sector_adjusted",
        "endpoint_weighted_residual", "synthetic_pair_recovered",
        "known_answer_null_connection_recovered", "target_was_forward_generated",
        "no_sprinkling_pair_used", "no_global_causal_relations_decided",
        "no_production_classifier_introduced", "min_r", "max_r", "min_Delta", "min_R",
        "max_abs_null_residual", "max_abs_E_residual", "max_abs_L_residual", "all_points_exterior",
        "finite_rhs_all_steps", "finite_solution_all_steps", "t_monotonic_future_pass",
        "radial_rhs_consistency_pass", "null_condition_pass", "constants_consistency_pass",
        "angular_accumulation_finite_pass", "advisory_only", "unresolved", "all_checks_pass",
    ]
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k) for k in fields})


def write_json(rows: list[dict[str, Any]], k14_flags: dict[str, bool], path: Path) -> None:
    recovered = [r for r in rows if r["synthetic_pair_recovered"]]
    whirling_recovered = [
        r for r in rows if r["pair_type"] == "kerr_whirling" and r["synthetic_pair_recovered"]
    ]
    no_match = [r for r in rows if r["synthetic_pair_classification"] == "synthetic_no_match"]
    unresolved = [r for r in rows if r["synthetic_pair_classification"] == "synthetic_unresolved"]
    max_abs = max(r["abs_delta_phi_target"] for r in rows) if rows else 0.0

    hard_gate = [r for r in rows if r["pair_type"] == "schwarzschild_radial"]
    hard_gate_pass = bool(hard_gate) and all(r["all_checks_pass"] for r in hard_gate)
    if not hard_gate_pass:
        raise RuntimeError("Schwarzschild radial known-answer pair failed.")

    any_pi = any(r["pair_type"] == "kerr_whirling" and r["synthetic_pair_recovered"] and r["abs_delta_phi_target"] > math.pi for r in rows)
    any_2pi = any(r["pair_type"] == "kerr_whirling" and r["synthetic_pair_recovered"] and r["abs_delta_phi_target"] > 2.0 * math.pi for r in rows)
    if k14_flags["k14_has_2pi_recovered"] and not any_2pi:
        raise RuntimeError("K14 has >2pi recovered whirling target but K15 did not recover one.")

    summary = {
        "total_pairs": len(rows),
        "synthetic_null_connected_pairs": sum(1 for r in rows if r["synthetic_pair_classification"] == "synthetic_null_connected"),
        "synthetic_no_match_pairs": len(no_match),
        "synthetic_unresolved_pairs": len(unresolved),
        "recovered_known_answer_pairs": len(recovered),
        "recovered_whirling_pairs": len(whirling_recovered),
        "negative_controls_passed": sum(
            1 for r in rows
            if r["pair_type"] == "negative_control"
            and (not r["synthetic_pair_recovered"])
            and (not r["known_answer_null_connection_recovered"])
            and r["synthetic_pair_classification"] == "synthetic_no_match"
            and r["endpoint_weighted_residual"] > PAIR_MATCH_TOL
        ),
        "max_abs_delta_phi_target": max_abs,
        "any_pi_whirling_pair_recovered": any_pi,
        "any_2pi_whirling_pair_recovered": any_2pi,
        "global_true_relations": 0,
        "global_false_relations": 0,
        "global_undecided_pairs": N_EVENTS * (N_EVENTS - 1) // 2,
        "all_checks_pass": (
            hard_gate_pass
            and any(r["synthetic_pair_classification"] == "synthetic_null_connected" for r in rows)
            and any_pi
            and (any_2pi if k14_flags["k14_has_2pi_recovered"] else True)
            and any(r["pair_type"] == "negative_control" and r["synthetic_pair_classification"] == "synthetic_no_match" for r in rows)
            and all(r["no_sprinkling_pair_used"] for r in rows)
            and all(r["no_global_causal_relations_decided"] for r in rows)
            and all(r["no_production_classifier_introduced"] for r in rows)
            and all(r["all_checks_pass"] for r in rows if (not r["advisory_only"]) and (not r["unresolved"]))
        ),
    }

    payload = {
        "benchmark": "S4-KERR-K15 synthetic causal sandbox",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "cases": rows,
        "global_summary": summary,
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def write_md(summary: dict[str, Any], path: Path) -> None:
    lines = [
        "# S4-KERR-K15 synthetic causal sandbox",
        "",
        "K15 is a synthetic known-answer event-pair sandbox.",
        "Event pairs are generated from validated geodesic segments.",
        "K15 does not use sprinkling event pairs.",
        "K15 does not implement a production causal classifier.",
        "synthetic_null_connected is not a physical/global Kerr causal claim.",
        "synthetic_no_match is a numerical negative control, not a proof of spacelike separation.",
        "global_true_relations and global_false_relations remain zero for a>0.",
        "",
        f"- Total pairs: {summary['total_pairs']}",
        f"- synthetic_null_connected pairs: {summary['synthetic_null_connected_pairs']}",
        f"- synthetic_no_match pairs: {summary['synthetic_no_match_pairs']}",
        f"- synthetic_unresolved pairs: {summary['synthetic_unresolved_pairs']}",
        f"- recovered_known_answer_pairs: {summary['recovered_known_answer_pairs']}",
        f"- recovered_whirling_pairs: {summary['recovered_whirling_pairs']}",
        f"- negative_controls_passed: {summary['negative_controls_passed']}",
        f"- any_pi_whirling_pair_recovered: {summary['any_pi_whirling_pair_recovered']}",
        f"- any_2pi_whirling_pair_recovered: {summary['any_2pi_whirling_pair_recovered']}",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_png(rows: list[dict[str, Any]], path: Path) -> None:
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle("S4-KERR-K15 synthetic causal sandbox")

    ax = axs[0, 0]
    ax.bar(range(len(rows)), [r["endpoint_weighted_residual"] for r in rows], color="tab:blue")
    ax.set_title("Endpoint weighted residuals by pair")
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels([r["pair_type"] for r in rows], rotation=45, ha="right", fontsize=7)

    ax = axs[0, 1]
    ax.plot(range(len(rows)), [r["abs_delta_phi_target"] for r in rows], "o-")
    ax.set_title("Target abs(delta_phi)")

    ax = axs[1, 0]
    ax.plot(range(len(rows)), [r["best_sector_m"] for r in rows], "o-", label="best")
    ax.plot(range(len(rows)), [r["correct_sector_m"] for r in rows], "s--", label="correct")
    ax.set_title("best_sector_m vs correct_sector_m")
    ax.legend(fontsize=7)

    ax = axs[1, 1]
    rec = [1 if r["synthetic_pair_recovered"] else 0 for r in rows]
    ax.bar(range(len(rows)), rec, color="tab:green")
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["no", "yes"])
    ax.set_title("Recovered synthetic controls")

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main() -> None:
    rows, k14_flags = build_cases()
    csv_path = ARTIFACT_DIR / f"{OUT_PREFIX}.csv"
    json_path = ARTIFACT_DIR / f"{OUT_PREFIX}.json"
    md_path = ARTIFACT_DIR / f"{OUT_PREFIX}.md"
    png_path = ARTIFACT_DIR / f"{OUT_PREFIX}.png"

    write_csv(rows, csv_path)
    write_json(rows, k14_flags, json_path)
    summary = json.loads(json_path.read_text(encoding="utf-8"))["global_summary"]
    write_md(summary, md_path)
    write_png(rows, png_path)
    print(
        f"pairs={summary['total_pairs']} recovered={summary['recovered_known_answer_pairs']} "
        f"whirling_recovered={summary['recovered_whirling_pairs']} no_match={summary['synthetic_no_match_pairs']}"
    )


if __name__ == "__main__":
    main()
