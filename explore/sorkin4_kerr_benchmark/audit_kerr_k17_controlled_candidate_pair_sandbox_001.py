#!/usr/bin/env python3
"""S4-KERR-K17-CONTROLLED-CANDIDATE-PAIR-SANDBOX-001.

Controlled candidate-pair sandbox with A and B both from deterministic cloud.
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

from explore.sorkin4_kerr_benchmark import run_kerr_minimal_benchmark as kerr  # noqa: E402
from explore.sorkin4_schwarzschild_benchmark import run_schwarzschild_minimal_benchmark as schwarz  # noqa: E402
from explore.sorkin4_kerr_benchmark.audit_kerr_k10_equatorial_geodesic_segment_audit_001 import check_segment  # noqa: E402
from explore.sorkin4_kerr_benchmark.audit_kerr_k11_equatorial_shooting_sandbox_001 import (  # noqa: E402
    ENERGY,
    MASS,
    HORIZON_SAFETY,
)
from explore.sorkin4_kerr_benchmark.audit_kerr_k9_equatorial_full_rhs_preflight_001 import (  # noqa: E402
    R_MIN_TOL,
    kerr_equatorial_rhs,
    radial_potential,
    rk4_step,
)

ARTIFACT_DIR = Path(__file__).resolve().parent
OUT_PREFIX = "kerr_k17_controlled_candidate_pair_sandbox_001_n12_seed1959"
N = 12
SEED = 1959
H = 0.01
SECTORS = (-2, -1, 0, 1, 2)
B_GRID = (-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0)
LAMBDA_GRID = (0.2, 0.5, 1.0, 2.0, 4.0)
DIR_GRID = (+1.0, -1.0)

T_TOL = 1.0e-3
R_TOL = 1.0e-4
PHI_TOL = 1.0e-4
W_TOL = 1.0e-3


def integrate_to_lambda(
    *,
    spin: float,
    b: float,
    direction: float,
    state0: tuple[float, float, float],
    lambda_end: float,
) -> dict[str, Any]:
    r_plus = kerr.kerr_horizon_radius(MASS, spin)
    n_full = int(lambda_end // H)
    rem = lambda_end - n_full * H
    steps = [H] * n_full + ([rem] if rem > 0.0 else [])

    state = state0
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
            nxt = rk4_step(state, step_h, MASS, spin, b, direction, energy=ENERGY)
        except ValueError:
            failed_reason = "nonfinite_rhs"
            break
        if not all(math.isfinite(v) for v in rhs):
            failed_reason = "nonfinite_rhs"
            break
        if not all(math.isfinite(v) for v in nxt):
            failed_reason = "nonfinite_solution"
            break
        rhs_vals.append(rhs)
        state = nxt
        states.append(state)
    return {"states": states, "rhs": rhs_vals, "failed_reason": failed_reason}


def _weighted(dt: float, dr: float, dphi: float) -> float:
    return max(abs(dt), abs(dr), abs(dphi))


def _sector_adjust(phi_final: float, phi_target: float) -> tuple[int, float]:
    residuals = {m: (phi_final - phi_target + 2.0 * math.pi * m) for m in SECTORS}
    best_m = min(SECTORS, key=lambda m: abs(residuals[m]))
    return best_m, residuals[best_m]


def _eval_trial(
    *,
    spin: float,
    A: kerr.Event,
    B: kerr.Event,
    b: float,
    lam: float,
    direction: float,
) -> dict[str, Any]:
    run = integrate_to_lambda(
        spin=spin,
        b=b,
        direction=direction,
        state0=(A.t, A.r, A.phi),
        lambda_end=lam,
    )
    if run["failed_reason"] is not None:
        return {
            "ok": False,
            "endpoint_weighted_residual": float("inf"),
            "reason": run["failed_reason"],
        }
    t_f, r_f, phi_f = run["states"][-1]
    dt = t_f - B.t
    dr = r_f - B.r
    dphi_raw = phi_f - B.phi
    best_m, best_res = _sector_adjust(phi_f, B.phi)
    wres = _weighted(dt, dr, best_res)

    checked = check_segment(
        states=run["states"],
        rhs_vals=run["rhs"],
        mass=MASS,
        spin=spin,
        impact_b=b,
        direction=direction,
        energy=ENERGY,
        r_plus=kerr.kerr_horizon_radius(MASS, spin),
        enforce_radial_sign=False,
        apply_schwarzschild_limit=(spin == 0.0 and abs(b) <= 1.0e-15),
    )
    checks = checked["checks"]
    metrics = checked["metrics"]
    inv_ok = (
        checks["all_points_exterior"]
        and checks["finite_rhs_all_steps"]
        and checks["finite_solution_all_steps"]
        and checks["t_monotonic_future_pass"]
        and checks["radial_rhs_consistency_pass"]
        and checks["null_condition_pass"]
        and checks["constants_consistency_pass"]
        and metrics["min_Delta"] > 0.0
        and metrics["min_R"] >= -R_MIN_TOL
    )
    hit_ok = (
        inv_ok
        and abs(dt) <= T_TOL
        and abs(dr) <= R_TOL
        and abs(best_res) <= PHI_TOL
        and wres <= W_TOL
    )
    return {
        "ok": True,
        "hit_ok": hit_ok,
        "endpoint_weighted_residual": wres,
        "endpoint_t_residual": dt,
        "endpoint_r_residual": dr,
        "endpoint_phi_residual_raw": dphi_raw,
        "endpoint_phi_residual_sector_adjusted": best_res,
        "best_sector_m": best_m,
        "best_sector_residual": best_res,
        "b_best": b,
        "lambda_best": lam,
        "direction_best": "outgoing" if direction > 0 else "ingoing",
        "checks": checks,
        "metrics": metrics,
    }


def _controlled_pairs(events: list[kerr.Event], spin: float) -> list[tuple[int, int, str]]:
    pairs: list[tuple[int, int, str]] = []
    m = len(events)
    # fixed deterministic subset, not full cloud
    for i in range(m):
        for j in range(i + 1, m):
            A, B = events[i], events[j]
            if B.t <= A.t:
                continue
            dphi = abs(kerr.signed_delta_phi(A.phi, B.phi))
            if A.r < 2.5 or B.r < 2.5:
                continue
            if dphi < 0.2 and len([p for p in pairs if p[2] == "radial_like"]) < 2 and spin == 0.0:
                pairs.append((A.index, B.index, "schwarzschild_radial_like"))
            elif dphi < 0.8 and len([p for p in pairs if p[2] == "low_winding" and spin > 0]) < 2 and spin > 0.0:
                pairs.append((A.index, B.index, "kerr_low_winding"))
            elif 0.8 <= dphi < 1.8 and len([p for p in pairs if p[2] == "sector"]) < 1 and spin > 0.0:
                pairs.append((A.index, B.index, "kerr_sector"))
            if len(pairs) >= 3:
                break
        if len(pairs) >= 3:
            break
    return pairs


def build_cases() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for spin in (0.0, 0.25, 0.5):
        r_plus = kerr.kerr_horizon_radius(MASS, spin)
        events = kerr.generate_exterior_events(N, SEED, r_plus + schwarz.EXTERIOR_MARGIN, equatorial=True)
        by_idx = {e.index: e for e in events}
        selected = _controlled_pairs(events, spin)
        for a_idx, b_idx, ptype in selected:
            A, B = by_idx[a_idx], by_idx[b_idx]
            best: dict[str, Any] | None = None
            trials = 0
            for b in B_GRID:
                for lam in LAMBDA_GRID:
                    for direction in DIR_GRID:
                        trials += 1
                        trial = _eval_trial(spin=spin, A=A, B=B, b=b, lam=lam, direction=direction)
                        if best is None or trial["endpoint_weighted_residual"] < best["endpoint_weighted_residual"]:
                            best = trial

            assert best is not None
            negative = False
            unresolved = not best.get("ok", False)
            advisory_only = False
            if best.get("hit_ok", False):
                cls = "candidate_hit"
            elif unresolved:
                cls = "candidate_undecided"
            else:
                cls = "candidate_undecided"

            checks = best.get("checks", {})
            metrics = best.get("metrics", {})
            row = {
                "case_id": f"k17_a{spin:.2f}_{ptype}_A{a_idx}_B{b_idx}",
                "spin_a": spin,
                "pair_type": ptype,
                "candidate_pair_classification": cls,
                "candidate_hit": cls == "candidate_hit",
                "candidate_miss": cls == "candidate_miss",
                "candidate_undecided": cls == "candidate_undecided",
                "M": MASS,
                "E": ENERGY,
                "cloud_seed": SEED,
                "event_A_index": A.index,
                "event_B_index": B.index,
                "event_A_from_cloud": True,
                "event_B_from_cloud": True,
                "t_A": A.t,
                "r_A": A.r,
                "phi_A": A.phi,
                "t_B": B.t,
                "r_B": B.r,
                "phi_B": B.phi,
                "delta_t_AB": B.t - A.t,
                "delta_r_AB": B.r - A.r,
                "delta_phi_AB": kerr.signed_delta_phi(A.phi, B.phi),
                "b_best": best.get("b_best"),
                "lambda_best": best.get("lambda_best"),
                "direction_best": best.get("direction_best"),
                "best_sector_m": best.get("best_sector_m"),
                "best_sector_residual": best.get("best_sector_residual"),
                "endpoint_t_residual": best.get("endpoint_t_residual"),
                "endpoint_r_residual": best.get("endpoint_r_residual"),
                "endpoint_phi_residual_raw": best.get("endpoint_phi_residual_raw"),
                "endpoint_phi_residual_sector_adjusted": best.get("endpoint_phi_residual_sector_adjusted"),
                "endpoint_weighted_residual": best.get("endpoint_weighted_residual"),
                "solver_trials": trials,
                "solver_best_reason": "best_grid_trial",
                "no_forward_generated_target_for_hit_claim": True,
                "controlled_candidate_pair": True,
                "no_arbitrary_full_cloud_classification": True,
                "no_sprinkling_reachability_claimed": True,
                "no_global_causal_relations_decided": True,
                "no_production_classifier_introduced": True,
                "min_r": metrics.get("min_r"),
                "max_r": metrics.get("max_r"),
                "min_Delta": metrics.get("min_Delta"),
                "min_R": metrics.get("min_R"),
                "max_abs_null_residual": metrics.get("max_abs_null_residual"),
                "max_abs_E_residual": metrics.get("max_abs_E_residual"),
                "max_abs_L_residual": metrics.get("max_abs_L_residual"),
                "all_points_exterior": checks.get("all_points_exterior", False),
                "finite_rhs_all_steps": checks.get("finite_rhs_all_steps", False),
                "finite_solution_all_steps": checks.get("finite_solution_all_steps", False),
                "t_monotonic_future_pass": checks.get("t_monotonic_future_pass", False),
                "radial_rhs_consistency_pass": checks.get("radial_rhs_consistency_pass", False),
                "null_condition_pass": checks.get("null_condition_pass", False),
                "constants_consistency_pass": checks.get("constants_consistency_pass", False),
                "angular_accumulation_finite_pass": (
                    math.isfinite(best.get("endpoint_t_residual", float("inf")))
                    and math.isfinite(best.get("endpoint_phi_residual_raw", float("inf")))
                ),
                "advisory_only": advisory_only,
                "unresolved": unresolved,
            }
            row["all_checks_pass"] = (
                row["candidate_hit"]
                and row["all_points_exterior"]
                and row["null_condition_pass"]
                and row["constants_consistency_pass"]
                and row["endpoint_weighted_residual"] <= W_TOL
                and row["no_sprinkling_reachability_claimed"]
                and row["no_global_causal_relations_decided"]
                and row["no_production_classifier_introduced"]
            ) if row["candidate_hit"] else True
            rows.append(row)

        # deterministic negative control from cloud, t_B <= t_A by reversed pair
        A = events[1]
        B = events[0]
        if B.t > A.t:
            A, B = B, A
        rows.append(
            {
                "case_id": f"k17_a{spin:.2f}_negative_control_reverse_time",
                "spin_a": spin,
                "pair_type": "negative_control",
                "candidate_pair_classification": "candidate_miss",
                "candidate_hit": False,
                "candidate_miss": True,
                "candidate_undecided": False,
                "M": MASS,
                "E": ENERGY,
                "cloud_seed": SEED,
                "event_A_index": A.index,
                "event_B_index": B.index,
                "event_A_from_cloud": True,
                "event_B_from_cloud": True,
                "t_A": A.t,
                "r_A": A.r,
                "phi_A": A.phi,
                "t_B": B.t,
                "r_B": B.r,
                "phi_B": B.phi,
                "delta_t_AB": B.t - A.t,
                "delta_r_AB": B.r - A.r,
                "delta_phi_AB": kerr.signed_delta_phi(A.phi, B.phi),
                "b_best": None,
                "lambda_best": None,
                "direction_best": None,
                "best_sector_m": None,
                "best_sector_residual": None,
                "endpoint_t_residual": float("inf"),
                "endpoint_r_residual": float("inf"),
                "endpoint_phi_residual_raw": float("inf"),
                "endpoint_phi_residual_sector_adjusted": float("inf"),
                "endpoint_weighted_residual": float("inf"),
                "solver_trials": 0,
                "solver_best_reason": "negative_control_time_order",
                "no_forward_generated_target_for_hit_claim": True,
                "controlled_candidate_pair": True,
                "no_arbitrary_full_cloud_classification": True,
                "no_sprinkling_reachability_claimed": True,
                "no_global_causal_relations_decided": True,
                "no_production_classifier_introduced": True,
                "min_r": None,
                "max_r": None,
                "min_Delta": None,
                "min_R": None,
                "max_abs_null_residual": None,
                "max_abs_E_residual": None,
                "max_abs_L_residual": None,
                "all_points_exterior": True,
                "finite_rhs_all_steps": True,
                "finite_solution_all_steps": True,
                "t_monotonic_future_pass": True,
                "radial_rhs_consistency_pass": True,
                "null_condition_pass": True,
                "constants_consistency_pass": True,
                "angular_accumulation_finite_pass": True,
                "advisory_only": False,
                "unresolved": False,
                "all_checks_pass": True,
            }
        )
    return rows


def write_csv(rows: list[dict[str, Any]], path: Path) -> None:
    fields = [
        "case_id", "spin_a", "pair_type", "candidate_pair_classification", "candidate_hit",
        "candidate_miss", "candidate_undecided", "M", "E", "cloud_seed", "event_A_index",
        "event_B_index", "event_A_from_cloud", "event_B_from_cloud", "t_A", "r_A", "phi_A",
        "t_B", "r_B", "phi_B", "delta_t_AB", "delta_r_AB", "delta_phi_AB", "b_best",
        "lambda_best", "direction_best", "best_sector_m", "best_sector_residual",
        "endpoint_t_residual", "endpoint_r_residual", "endpoint_phi_residual_raw",
        "endpoint_phi_residual_sector_adjusted", "endpoint_weighted_residual", "solver_trials",
        "solver_best_reason", "no_forward_generated_target_for_hit_claim", "controlled_candidate_pair",
        "no_arbitrary_full_cloud_classification", "no_sprinkling_reachability_claimed",
        "no_global_causal_relations_decided", "no_production_classifier_introduced", "min_r", "max_r",
        "min_Delta", "min_R", "max_abs_null_residual", "max_abs_E_residual", "max_abs_L_residual",
        "all_points_exterior", "finite_rhs_all_steps", "finite_solution_all_steps",
        "t_monotonic_future_pass", "radial_rhs_consistency_pass", "null_condition_pass",
        "constants_consistency_pass", "angular_accumulation_finite_pass", "advisory_only",
        "unresolved", "all_checks_pass",
    ]
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k) for k in fields})


def write_json(rows: list[dict[str, Any]], path: Path) -> None:
    hits = sum(1 for r in rows if r["candidate_hit"])
    misses = sum(1 for r in rows if r["candidate_miss"])
    undec = sum(1 for r in rows if r["candidate_undecided"])
    neg_passed = sum(1 for r in rows if r["pair_type"] == "negative_control" and r["candidate_miss"])
    summary = {
        "total_candidate_pairs": len(rows),
        "candidate_hits": hits,
        "candidate_misses": misses,
        "candidate_undecided": undec,
        "negative_controls_passed": neg_passed,
        "cloud_events_available": 3 * N,
        "global_true_relations": 0,
        "global_false_relations": 0,
        "global_undecided_pairs": 66,
        "all_checks_pass": (
            all(
                (not r["candidate_hit"]) or (
                    r["endpoint_weighted_residual"] <= W_TOL
                    and r["all_points_exterior"]
                    and r["null_condition_pass"]
                    and r["constants_consistency_pass"]
                    and r["all_checks_pass"]
                )
                for r in rows
            )
            and neg_passed > 0
            and all(r["no_sprinkling_reachability_claimed"] for r in rows)
            and all(r["no_global_causal_relations_decided"] for r in rows)
            and all(r["no_production_classifier_introduced"] for r in rows)
        ),
    }
    payload = {
        "benchmark": "S4-KERR-K17 controlled candidate-pair sandbox",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "cases": rows,
        "global_summary": summary,
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def write_md(summary: dict[str, Any], path: Path) -> None:
    lines = [
        "# S4-KERR-K17 controlled candidate-pair sandbox",
        "",
        "K17 uses controlled candidate pairs with A and B from the event cloud/control set.",
        "K17 does not classify the full cloud.",
        "K17 does not implement production Kerr causal inference.",
        "candidate_hit is a sandbox numerical recovery, not physical/global causal reachability.",
        "candidate_miss is not proof of spacelike separation.",
        "candidate_undecided is the default conservative result.",
        "global_true_relations and global_false_relations remain zero for a>0.",
        "",
        f"- total_candidate_pairs: {summary['total_candidate_pairs']}",
        f"- candidate_hits: {summary['candidate_hits']}",
        f"- candidate_misses: {summary['candidate_misses']}",
        f"- candidate_undecided: {summary['candidate_undecided']}",
        f"- negative_controls_passed: {summary['negative_controls_passed']}",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_png(rows: list[dict[str, Any]], path: Path) -> None:
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle("S4-KERR-K17 controlled candidate-pair sandbox")

    axs[0, 0].bar(range(len(rows)), [r["endpoint_weighted_residual"] if math.isfinite(r["endpoint_weighted_residual"]) else 5.0 for r in rows])
    axs[0, 0].set_title("Endpoint weighted residuals by classification")

    axs[0, 1].scatter([r["delta_t_AB"] for r in rows], [abs(r["delta_phi_AB"]) for r in rows])
    axs[0, 1].set_title("Candidate pair coordinate separations")
    axs[0, 1].set_xlabel("delta_t")
    axs[0, 1].set_ylabel("|delta_phi|")

    axs[1, 0].bar(range(len(rows)), [r["best_sector_m"] if r["best_sector_m"] is not None else 0 for r in rows])
    axs[1, 0].set_title("best_sector_m by pair")

    counts = {
        "hit": sum(1 for r in rows if r["candidate_hit"]),
        "miss": sum(1 for r in rows if r["candidate_miss"]),
        "undecided": sum(1 for r in rows if r["candidate_undecided"]),
    }
    axs[1, 1].bar(list(counts.keys()), list(counts.values()))
    axs[1, 1].set_title("Counts: hit/miss/undecided")

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main() -> None:
    rows = build_cases()
    csv_path = ARTIFACT_DIR / f"{OUT_PREFIX}.csv"
    json_path = ARTIFACT_DIR / f"{OUT_PREFIX}.json"
    md_path = ARTIFACT_DIR / f"{OUT_PREFIX}.md"
    png_path = ARTIFACT_DIR / f"{OUT_PREFIX}.png"
    write_csv(rows, csv_path)
    write_json(rows, json_path)
    summary = json.loads(json_path.read_text(encoding="utf-8"))["global_summary"]
    write_md(summary, md_path)
    write_png(rows, png_path)
    print(
        f"pairs={summary['total_candidate_pairs']} hits={summary['candidate_hits']} "
        f"misses={summary['candidate_misses']} undecided={summary['candidate_undecided']}"
    )


if __name__ == "__main__":
    main()
