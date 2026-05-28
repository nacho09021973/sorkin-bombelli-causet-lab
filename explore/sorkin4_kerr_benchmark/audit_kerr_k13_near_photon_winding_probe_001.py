#!/usr/bin/env python3
"""S4-KERR-K13-NEAR-PHOTON-WINDING-PROBE-001.

Cautious physical angular-accumulation probe near equatorial photon spheres.
No causal classification or sprinkling reachability is performed.
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
    delta,
    kerr_equatorial_rhs,
    photon_impact_parameter,
    photon_sphere_radius_pro,
    photon_sphere_radius_retro,
    radial_potential,
    rk4_step,
)


ARTIFACT_DIR = Path(__file__).resolve().parent
OUT_PREFIX = "kerr_k13_near_photon_winding_probe_001_n12_seed1959"
SPINS = (0.0, 0.25, 0.5)
DELTA_RS = (0.05, 0.10, 0.20)
EPSILONS = (-0.02, -0.01, 0.01, 0.02)
MAX_STEPS = 240
STEP_H = 0.01
N_EVENTS = 12
R_MIN_SAFE = 1.0e-12
HARD_NULL_TOL = 1.0e-7
HARD_EL_TOL = 1.0e-7
HARD_RHS_TOL = 1.0e-6


def integrate_segment(
    *,
    spin: float,
    b: float,
    r0: float,
    n_steps: int = MAX_STEPS,
    h: float = STEP_H,
    direction: float = +1.0,
) -> dict[str, Any]:
    r_plus = kerr_horizon_radius(MASS, spin)
    state = (0.0, r0, 0.0)
    states = [state]
    rhs_vals: list[tuple[float, float, float]] = []
    stop_reason = "normal"

    for _ in range(n_steps):
        t, r, _ = state
        _ = t
        if r <= r_plus + HORIZON_SAFETY:
            stop_reason = "safety_stop"
            break
        if delta(r, MASS, spin) <= 0.0:
            stop_reason = "safety_stop"
            break
        rpot = radial_potential(r, MASS, spin, b, energy=ENERGY)
        if rpot < -R_MIN_TOL:
            stop_reason = "invariant_failure"
            break
        if rpot <= R_MIN_SAFE:
            stop_reason = "safety_stop"
            break

        rhs = kerr_equatorial_rhs(state, MASS, spin, b, direction, energy=ENERGY)
        if not all(math.isfinite(x) for x in rhs):
            stop_reason = "invariant_failure"
            break
        rhs_vals.append(rhs)
        nxt = rk4_step(state, h, MASS, spin, b, direction, energy=ENERGY)
        if not all(math.isfinite(x) for x in nxt):
            stop_reason = "invariant_failure"
            break
        state = nxt
        states.append(state)

    return {
        "states": states,
        "rhs_vals": rhs_vals,
        "stop_reason": stop_reason,
    }


def _row(
    *,
    case_id: str,
    spin: float,
    branch: str,
    b: float,
    b_ph: float,
    epsilon_b: float,
    r_ph: float,
    r0: float,
    delta_r: float,
    run: dict[str, Any],
    advisory_only: bool,
    unresolved: bool,
) -> dict[str, Any]:
    states = run["states"]
    rhs_vals = run["rhs_vals"]
    stop_reason = run["stop_reason"]
    r_plus = kerr_horizon_radius(MASS, spin)

    checked = check_segment(
        states=states,
        rhs_vals=rhs_vals,
        mass=MASS,
        spin=spin,
        impact_b=b,
        direction=+1.0,
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

    weak_winding_pass = abs_delta_phi > (0.5 * math.pi)
    strong_winding_pass = abs_delta_phi > math.pi
    full_turn_probe_pass = abs_delta_phi > (2.0 * math.pi)
    angular_accumulation_finite_pass = math.isfinite(delta_t) and math.isfinite(delta_phi_raw)

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
        and metrics["max_abs_radial_rhs_residual"] <= HARD_RHS_TOL
    )

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
        "M": MASS,
        "E": ENERGY,
        "b": b,
        "b_ph": b_ph,
        "epsilon_b": epsilon_b,
        "r_plus": r_plus,
        "r_ph": r_ph,
        "r0": r0,
        "delta_r": delta_r,
        "lambda_end": (len(states) - 1) * STEP_H,
        "n_steps": len(states) - 1,
        "stop_reason": stop_reason,
        "delta_t": delta_t,
        "delta_phi_raw": delta_phi_raw,
        "abs_delta_phi": abs_delta_phi,
        "winding_m_estimate": winding_m_estimate,
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
        r_plus = kerr_horizon_radius(MASS, spin)
        branches = (
            ("prograde", photon_sphere_radius_pro(MASS, spin), photon_impact_parameter(photon_sphere_radius_pro(MASS, spin), MASS, spin, prograde=True)),
            ("retrograde", photon_sphere_radius_retro(MASS, spin), photon_impact_parameter(photon_sphere_radius_retro(MASS, spin), MASS, spin, prograde=False)),
        )
        for branch, r_ph, b_ph in branches:
            for delta_r in DELTA_RS:
                for eps in EPSILONS:
                    scanned += 1
                    r0 = r_ph + delta_r
                    b = b_ph * (1.0 + eps)
                    unresolved = False
                    advisory_only = True
                    stop_reason = "normal"

                    if r0 <= r_plus + HORIZON_SAFETY:
                        unresolved = True
                        stop_reason = "safety_stop"
                        run = {"states": [(0.0, r0, 0.0)], "rhs_vals": [], "stop_reason": stop_reason}
                    else:
                        run = integrate_segment(spin=spin, b=b, r0=r0)
                        stop_reason = run["stop_reason"]
                        unresolved = stop_reason != "normal"
                        # Accept as non-advisory if fully controlled invariants.
                        advisory_only = not (stop_reason == "normal")

                    case_id = (
                        f"k13_a{spin:.2f}_{branch}_dr{delta_r:.2f}_eps{eps:+.2f}"
                    )
                    row = _row(
                        case_id=case_id,
                        spin=spin,
                        branch=branch,
                        b=b,
                        b_ph=b_ph,
                        epsilon_b=eps,
                        r_ph=r_ph,
                        r0=r0,
                        delta_r=delta_r,
                        run=run,
                        advisory_only=advisory_only,
                        unresolved=unresolved,
                    )
                    if stop_reason != "normal":
                        row["stop_reason"] = stop_reason
                    rows.append(row)
    return rows, scanned


def write_csv(rows: list[dict[str, Any]], path: Path) -> None:
    fields = [
        "case_id", "spin_a", "branch", "M", "E", "b", "b_ph", "epsilon_b", "r_plus", "r_ph",
        "r0", "delta_r", "lambda_end", "n_steps", "stop_reason", "delta_t", "delta_phi_raw",
        "abs_delta_phi", "winding_m_estimate", "min_r", "max_r", "min_Delta", "min_R",
        "max_abs_null_residual", "max_abs_E_residual", "max_abs_L_residual",
        "all_points_exterior", "finite_rhs_all_steps", "finite_solution_all_steps",
        "t_monotonic_future_pass", "radial_rhs_consistency_pass", "null_condition_pass",
        "constants_consistency_pass", "angular_accumulation_finite_pass", "weak_winding_pass",
        "strong_winding_pass", "full_turn_probe_pass", "no_sprinkling_pair_used",
        "no_global_causal_relations_decided", "no_causal_classifier_introduced",
        "advisory_only", "unresolved", "all_checks_pass",
    ]
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k) for k in fields})


def write_json(rows: list[dict[str, Any]], scanned: int, path: Path) -> None:
    non_adv_non_unres = [r for r in rows if (not r["advisory_only"]) and (not r["unresolved"])]
    summary = {
        "total_candidates_scanned": scanned,
        "total_cases_recorded": len(rows),
        "accepted_cases": len(non_adv_non_unres),
        "advisory_cases": sum(1 for r in rows if r["advisory_only"]),
        "unresolved_cases": sum(1 for r in rows if r["unresolved"]),
        "max_abs_delta_phi": max(r["abs_delta_phi"] for r in rows) if rows else 0.0,
        "any_weak_winding_pass": any(r["weak_winding_pass"] for r in rows),
        "any_strong_winding_pass": any(r["strong_winding_pass"] for r in rows),
        "any_full_turn_probe_pass": any(r["full_turn_probe_pass"] for r in rows),
        "global_true_relations": 0,
        "global_false_relations": 0,
        "global_undecided_pairs": N_EVENTS * (N_EVENTS - 1) // 2,
        "all_checks_pass": all(r["all_checks_pass"] for r in non_adv_non_unres),
    }
    payload = {
        "benchmark": "S4-KERR-K13 near-photon physical winding probe",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "cases": [{k: v for k, v in r.items() if not k.startswith("_")} for r in rows],
        "global_summary": summary,
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def write_md(summary: dict[str, Any], path: Path) -> None:
    lines = [
        "# S4-KERR-K13 near-photon physical winding probe",
        "",
        "K13 is a physical angular-accumulation probe near the photon sphere.",
        "It does not classify causal reachability.",
        "It does not use sprinkling event pairs.",
        "It does not implement a production Kerr causal classifier.",
        "strong_winding_pass means large angular accumulation in an initial-value segment, not endpoint reachability.",
        "Failure to reach |Delta_phi| > pi is not a physics failure if invariants remain controlled; it only means the probe did not find a strong-winding safe case.",
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
    fig, axs = plt.subplots(2, 2, figsize=(11, 8))
    fig.suptitle("S4-KERR-K13 near-photon physical winding probe")

    x = list(range(len(rows)))
    axs[0, 0].plot(x, [r["abs_delta_phi"] for r in rows], "o-")
    axs[0, 0].set_title("abs(delta_phi) by case")
    axs[0, 0].set_xlabel("case index")

    top = sorted(rows, key=lambda r: r["abs_delta_phi"], reverse=True)[:4]
    ax = axs[0, 1]
    for r in top:
        states = r["_states"]
        lam = [i * STEP_H for i in range(len(states))]
        rr = [s[1] for s in states]
        ax.plot(lam, rr, label=r["case_id"])
    ax.set_title("r(lambda) top winding cases")
    ax.set_xlabel("lambda")
    ax.legend(fontsize=6)

    ax = axs[1, 0]
    for r in top:
        states = r["_states"]
        lam = [i * STEP_H for i in range(len(states))]
        pp = [s[2] for s in states]
        ax.plot(lam, pp, label=r["case_id"])
    ax.set_title("phi(lambda) top winding cases")
    ax.set_xlabel("lambda")

    ax = axs[1, 1]
    nulls = [r["max_abs_null_residual"] for r in rows]
    ers = [r["max_abs_E_residual"] for r in rows]
    lrs = [r["max_abs_L_residual"] for r in rows]
    ax.semilogy(x, nulls, "o-", label="null")
    ax.semilogy(x, ers, "s-", label="E")
    ax.semilogy(x, lrs, "^-", label="L")
    ax.set_title("Invariant residuals by case")
    ax.set_xlabel("case index")
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
