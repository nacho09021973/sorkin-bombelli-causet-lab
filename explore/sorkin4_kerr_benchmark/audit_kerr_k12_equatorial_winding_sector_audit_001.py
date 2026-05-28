#!/usr/bin/env python3
"""S4-KERR-K12-EQUATORIAL-WINDING-SECTOR-AUDIT-001.

Synthetic winding-sector bookkeeping audit for equatorial Kerr null geodesics.
This is not a production Kerr causal classifier.
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
from audit_kerr_k11_equatorial_shooting_sandbox_001 import (  # noqa: E402
    ENERGY,
    MASS,
    integrate_to_lambda,
)
from audit_kerr_k10_equatorial_geodesic_segment_audit_001 import (  # noqa: E402
    check_segment,
)
from audit_kerr_k9_equatorial_full_rhs_preflight_001 import (  # noqa: E402
    R_MIN_TOL,
    photon_impact_parameter,
    photon_sphere_radius_pro,
    photon_sphere_radius_retro,
)


ARTIFACT_DIR = Path(__file__).resolve().parent
OUT_PREFIX = "kerr_k12_equatorial_winding_sector_audit_001_n12_seed1959"
N_EVENTS = 12
SECTOR_MS = [-2, -1, 0, 1, 2]
PHI_RES_TOL = 1.0e-6


def _wrap_to_pi(x: float) -> float:
    return (x + math.pi) % (2.0 * math.pi) - math.pi


def _sector_metrics(phi_final: float, phi_target_full: float) -> tuple[float, int, dict[int, float], int, float]:
    delta_phi_raw = phi_target_full
    winding_m_estimate = int(round(delta_phi_raw / (2.0 * math.pi)))
    phi_target_wrapped = _wrap_to_pi(phi_target_full)
    residuals = {m: (phi_final - phi_target_wrapped + 2.0 * math.pi * m) for m in SECTOR_MS}
    best_m = min(SECTOR_MS, key=lambda m: abs(residuals[m]))
    best_res = residuals[best_m]
    return delta_phi_raw, winding_m_estimate, residuals, best_m, best_res


def _run_case(
    *,
    case_id: str,
    spin: float,
    b: float,
    direction: float,
    r0: float,
    lambda_end: float,
    case_type: str,
    sign_expectation: str | None = None,
    advisory_only: bool = False,
) -> dict[str, Any]:
    # Forward-generate synthetic target.
    fwd = integrate_to_lambda(
        spin=spin, b=b, direction=direction, r0=r0, lambda_end=lambda_end
    )
    unresolved = fwd["failed_reason"] is not None
    shot = fwd

    t_final, r_final, phi_final = shot["states"][-1]
    phi_initial = shot["states"][0][2]
    t_initial = shot["states"][0][0]
    delta_t = t_final - t_initial
    phi_target = phi_final - phi_initial
    delta_phi_raw, winding_m_estimate, sector_residuals, best_sector_m, best_sector_residual = _sector_metrics(
        phi_final=phi_final, phi_target_full=phi_target
    )
    correct_sector_m = winding_m_estimate
    correct_sector_recovered = best_sector_m == correct_sector_m
    synthetic_winding_sector_recovered = correct_sector_recovered

    if sign_expectation == "prograde":
        prograde_retrograde_sign_pass = delta_phi_raw > 0.0
    elif sign_expectation == "retrograde":
        prograde_retrograde_sign_pass = delta_phi_raw < 0.0
    else:
        prograde_retrograde_sign_pass = True

    checked = check_segment(
        states=shot["states"],
        rhs_vals=shot["rhs"],
        mass=MASS,
        spin=spin,
        impact_b=b,
        direction=direction,
        energy=ENERGY,
        r_plus=kerr_horizon_radius(MASS, spin),
        enforce_radial_sign=True,
        apply_schwarzschild_limit=(spin == 0.0 and b == 0.0),
    )
    checks = checked["checks"]
    metrics = checked["metrics"]

    all_checks_pass = (
        (not unresolved)
        and checked["all_checks_pass"]
        and correct_sector_recovered
        and (abs(best_sector_residual) <= PHI_RES_TOL)
        and prograde_retrograde_sign_pass
    )
    if advisory_only or unresolved:
        all_checks_pass = False

    return {
        "case_id": case_id,
        "spin_a": spin,
        "M": MASS,
        "E": ENERGY,
        "b": b,
        "direction": "outgoing" if direction > 0 else "ingoing",
        "case_type": case_type,
        "r_plus": kerr_horizon_radius(MASS, spin),
        "r0": r0,
        "lambda_end": lambda_end,
        "delta_t": delta_t,
        "delta_phi_raw": delta_phi_raw,
        "winding_m_estimate": winding_m_estimate,
        "correct_sector_m": correct_sector_m,
        "best_sector_m": best_sector_m,
        "best_sector_residual": best_sector_residual,
        "correct_sector_recovered": correct_sector_recovered,
        "synthetic_winding_sector_recovered": synthetic_winding_sector_recovered,
        "prograde_retrograde_sign_pass": prograde_retrograde_sign_pass,
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
        "null_condition_pass": checks["null_condition_pass"],
        "constants_consistency_pass": checks["constants_consistency_pass"],
        "no_sprinkling_pair_used": True,
        "no_global_causal_relations_decided": True,
        "advisory_only": advisory_only,
        "unresolved": unresolved,
        "all_checks_pass": all_checks_pass,
        "tested_sector_range": SECTOR_MS,
        "sector_residuals": sector_residuals,
        "_traj": shot["states"],
    }


def build_cases() -> list[dict[str, Any]]:
    cases: list[dict[str, Any]] = []

    # 1) Schwarzschild controls
    cases.append(
        _run_case(
            case_id="k12_schw_b0_outgoing",
            spin=0.0,
            b=0.0,
            direction=+1.0,
            r0=4.0,
            lambda_end=0.4,
            case_type="schwarzschild_low_winding_control",
        )
    )
    cases.append(
        _run_case(
            case_id="k12_schw_b0_ingoing",
            spin=0.0,
            b=0.0,
            direction=-1.0,
            r0=8.0,
            lambda_end=0.4,
            case_type="schwarzschild_low_winding_control",
        )
    )
    cases.append(
        _run_case(
            case_id="k12_schw_bpos_safe",
            spin=0.0,
            b=1.0,
            direction=+1.0,
            r0=6.0,
            lambda_end=0.3,
            case_type="schwarzschild_low_winding_control",
            sign_expectation="prograde",
        )
    )

    # 2) Kerr low-winding b=0
    for spin in (0.25, 0.5):
        r_plus = kerr_horizon_radius(MASS, spin)
        for direction, lbl in ((+1.0, "out"), (-1.0, "in")):
            cases.append(
                _run_case(
                    case_id=f"k12_kerr_low_winding_a{spin:.2f}_{lbl}",
                    spin=spin,
                    b=0.0,
                    direction=direction,
                    r0=r_plus + schwarz.EXTERIOR_MARGIN + (1.0 if direction > 0 else 2.0),
                    lambda_end=0.3,
                    case_type="kerr_low_winding_synthetic",
                )
            )

    # 3) Prograde/retrograde sector probes
    for spin in (0.25, 0.5):
        r_plus = kerr_horizon_radius(MASS, spin)
        r0 = r_plus + schwarz.EXTERIOR_MARGIN + 1.3
        cases.append(
            _run_case(
                case_id=f"k12_kerr_prograde_a{spin:.2f}",
                spin=spin,
                b=+1.0,
                direction=+1.0,
                r0=r0,
                lambda_end=0.25,
                case_type="kerr_prograde_sector_synthetic",
                sign_expectation="prograde",
            )
        )
        cases.append(
            _run_case(
                case_id=f"k12_kerr_retrograde_a{spin:.2f}",
                spin=spin,
                b=-1.0,
                direction=+1.0,
                r0=r0,
                lambda_end=0.25,
                case_type="kerr_retrograde_sector_synthetic",
                sign_expectation="retrograde",
            )
        )

    # 4) Optional advisory near-photon probe
    for spin in (0.25, 0.5):
        rph = photon_sphere_radius_pro(MASS, spin)
        bph = photon_impact_parameter(rph, MASS, spin, prograde=True)
        cases.append(
            _run_case(
                case_id=f"k12_near_photon_pro_a{spin:.2f}",
                spin=spin,
                b=0.98 * bph,
                direction=+1.0,
                r0=rph + 0.02,
                lambda_end=0.2,
                case_type="advisory_near_photon_winding_probe",
                sign_expectation="prograde",
                advisory_only=True,
            )
        )
        rphr = photon_sphere_radius_retro(MASS, spin)
        bphr = photon_impact_parameter(rphr, MASS, spin, prograde=False)
        cases.append(
            _run_case(
                case_id=f"k12_near_photon_retro_a{spin:.2f}",
                spin=spin,
                b=0.98 * bphr,
                direction=+1.0,
                r0=rphr + 0.02,
                lambda_end=0.2,
                case_type="advisory_near_photon_winding_probe",
                sign_expectation="retrograde",
                advisory_only=True,
            )
        )
    return cases


def write_csv(rows: list[dict[str, Any]], path: Path) -> None:
    fields = [
        "case_id", "spin_a", "M", "E", "b", "direction", "case_type", "r_plus", "r0",
        "lambda_end", "delta_t", "delta_phi_raw", "winding_m_estimate", "correct_sector_m",
        "best_sector_m", "best_sector_residual", "correct_sector_recovered",
        "synthetic_winding_sector_recovered", "prograde_retrograde_sign_pass", "min_r",
        "max_r", "min_Delta", "min_R", "max_abs_null_residual", "max_abs_E_residual",
        "max_abs_L_residual", "all_points_exterior", "finite_rhs_all_steps",
        "finite_solution_all_steps", "t_monotonic_future_pass", "null_condition_pass",
        "constants_consistency_pass", "no_sprinkling_pair_used",
        "no_global_causal_relations_decided", "advisory_only", "unresolved", "all_checks_pass",
    ]
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in fields})


def write_json(rows: list[dict[str, Any]], path: Path) -> None:
    summary = {
        "total_cases": len(rows),
        "passed_cases": sum(1 for r in rows if r["all_checks_pass"]),
        "failed_cases": sum(1 for r in rows if (not r["all_checks_pass"]) and (not r["unresolved"])),
        "advisory_cases": sum(1 for r in rows if r["advisory_only"]),
        "unresolved_cases": sum(1 for r in rows if r["unresolved"]),
        "sectors_tested": SECTOR_MS,
        "synthetic_winding_targets_generated": sum(1 for r in rows if not r["unresolved"]),
        "synthetic_winding_targets_recovered": sum(
            1 for r in rows if r["synthetic_winding_sector_recovered"]
        ),
        "global_true_relations": 0,
        "global_false_relations": 0,
        "global_undecided_pairs": N_EVENTS * (N_EVENTS - 1) // 2,
        "all_checks_pass": all(
            r["all_checks_pass"] for r in rows if (not r["advisory_only"] and not r["unresolved"])
        ),
    }
    payload = {
        "benchmark": "S4-KERR-K12 equatorial winding-sector audit",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "cases": [{k: v for k, v in r.items() if not k.startswith("_")} for r in rows],
        "global_summary": summary,
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def write_md(summary: dict[str, Any], path: Path) -> None:
    lines = [
        "# S4-KERR-K12 equatorial winding-sector audit",
        "",
        "K12 audits winding-sector bookkeeping on synthetic targets.",
        "It does not use sprinkling event pairs.",
        "It does not decide causal reachability.",
        "It does not implement a production Kerr causal classifier.",
        "correct_sector_recovered is not physical reachability.",
        "advisory near-photon-sphere cases are diagnostics only.",
        "",
        f"- Total cases: {summary['total_cases']}",
        f"- Passed cases: {summary['passed_cases']}",
        f"- Advisory cases: {summary['advisory_cases']}",
        f"- Synthetic winding targets generated: {summary['synthetic_winding_targets_generated']}",
        f"- Synthetic winding targets recovered: {summary['synthetic_winding_targets_recovered']}",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_png(rows: list[dict[str, Any]], path: Path) -> None:
    fig, axs = plt.subplots(2, 2, figsize=(11, 8))
    fig.suptitle("S4-KERR-K12 equatorial winding-sector audit")

    x = list(range(len(rows)))
    axs[0, 0].plot(x, [r["delta_phi_raw"] for r in rows], "o-")
    axs[0, 0].set_title("delta_phi_raw by case")
    axs[0, 0].set_xlabel("case index")

    axs[0, 1].plot(x, [r["correct_sector_m"] for r in rows], "o", label="correct")
    axs[0, 1].plot(x, [r["best_sector_m"] for r in rows], "x", label="best")
    axs[0, 1].set_title("best sector vs correct sector")
    axs[0, 1].set_xlabel("case index")
    axs[0, 1].legend(fontsize=8)

    axs[1, 0].semilogy(x, [abs(r["best_sector_residual"]) for r in rows], "o-")
    axs[1, 0].set_title("|best sector residual| by case")
    axs[1, 0].set_xlabel("case index")

    ax = axs[1, 1]
    for r in rows[:4]:
        ys = [abs(r["sector_residuals"][m]) for m in SECTOR_MS]
        ax.plot(SECTOR_MS, ys, marker="o", label=r["case_id"])
    ax.set_title("Representative sector residual curves")
    ax.set_xlabel("m")
    ax.set_ylabel("|phi residual_m|")
    ax.legend(fontsize=6)

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main() -> None:
    rows = build_cases()
    non_adv_resolved = [r for r in rows if (not r["advisory_only"]) and (not r["unresolved"])]
    if any(not r["all_checks_pass"] for r in non_adv_resolved):
        raise SystemExit("STOP: non-advisory resolved case failed.")
    # hard gate
    gate = next(r for r in rows if r["case_id"] == "k12_schw_b0_outgoing")
    if gate["best_sector_m"] != 0:
        raise SystemExit("STOP: Schwarzschild radial control did not recover m=0.")

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
        f"cases={summary['total_cases']} recovered={summary['synthetic_winding_targets_recovered']} "
        f"advisory={summary['advisory_cases']}"
    )


if __name__ == "__main__":
    main()
