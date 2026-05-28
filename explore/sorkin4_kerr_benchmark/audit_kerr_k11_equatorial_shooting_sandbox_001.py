#!/usr/bin/env python3
"""S4-KERR-K11-EQUATORIAL-SHOOTING-SANDBOX-001.

Synthetic known-answer shooting sandbox only.
Targets are forward generated from controlled initial data.
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
    kerr_equatorial_rhs,
    radial_potential,
    rk4_step,
)
from audit_kerr_k10_equatorial_geodesic_segment_audit_001 import (  # noqa: E402
    check_segment,
)


N_EVENTS = 12
ARTIFACT_DIR = Path(__file__).resolve().parent
OUT_PREFIX = "kerr_k11_equatorial_shooting_sandbox_001_n12_seed1959"
COMMAND = (
    "python3 explore/sorkin4_kerr_benchmark/"
    "audit_kerr_k11_equatorial_shooting_sandbox_001.py"
)

MASS = 1.0
ENERGY = 1.0
H = 0.01
HORIZON_SAFETY = 1.0e-6

T_TOL = 1.0e-5
R_TOL = 1.0e-6
PHI_TOL = 1.0e-6
B_TOL = 1.0e-5
WEIGHTED_TOL = 1.0e-5


def integrate_to_lambda(
    *,
    spin: float,
    b: float,
    direction: float,
    r0: float,
    lambda_end: float,
    h: float = H,
) -> dict[str, Any]:
    r_plus = kerr_horizon_radius(MASS, spin)
    state = (0.0, r0, 0.0)
    states = [state]
    rhs_vals = []
    failed_reason = None
    n_full = int(lambda_end // h)
    rem = lambda_end - n_full * h

    for step_h in [h] * n_full + ([rem] if rem > 0 else []):
        rr = state[1]
        rpot = radial_potential(rr, MASS, spin, b, energy=ENERGY)
        if rpot < -R_MIN_TOL:
            failed_reason = f"R< -tol at r={rr:.9g}, R={rpot:.9g}"
            break
        if rr <= r_plus + HORIZON_SAFETY:
            failed_reason = f"r reached safety boundary at r={rr:.9g}"
            break
        rhs = kerr_equatorial_rhs(state, MASS, spin, b, direction, energy=ENERGY)
        rhs_vals.append(rhs)
        if not all(math.isfinite(x) for x in rhs):
            failed_reason = "non-finite RHS"
            break
        nxt = rk4_step(state, step_h, MASS, spin, b, direction, energy=ENERGY)
        if not all(math.isfinite(x) for x in nxt):
            failed_reason = "non-finite solution"
            break
        if nxt[1] <= r_plus + HORIZON_SAFETY:
            failed_reason = f"next step approaches horizon safety boundary r={nxt[1]:.9g}"
            break
        state = nxt
        states.append(state)
    return {"states": states, "rhs": rhs_vals, "failed_reason": failed_reason}


def _endpoint(states: list[tuple[float, float, float]]) -> tuple[float, float, float]:
    t, r, phi = states[-1]
    return t, r, phi


def _weighted_residual(dt: float, dr: float, dphi: float) -> float:
    return max(abs(dt) / T_TOL, abs(dr) / R_TOL, abs(dphi) / PHI_TOL)


def _bisection(
    func,
    left: float,
    right: float,
    tol: float = 1.0e-7,
    max_iter: int = 40,
) -> tuple[bool, float | None, int]:
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


def _case(
    *,
    case_id: str,
    spin: float,
    direction: float,
    r0: float,
    lambda_true: float,
    b_true: float,
    case_type: str,
    recover_b: bool,
    recover_lambda: bool,
    b_bracket: tuple[float, float] | None = None,
) -> dict[str, Any]:
    # forward-generated synthetic target
    fwd = integrate_to_lambda(
        spin=spin, b=b_true, direction=direction, r0=r0, lambda_end=lambda_true
    )
    fwd_fail = fwd["failed_reason"]
    target_generated = fwd_fail is None
    t_target, r_target, phi_target = _endpoint(fwd["states"])

    b_recovered = None
    lambda_recovered = None
    bracket_found = False
    solver_converged = False
    solver_iterations = 0
    unresolved = False
    advisory_only = False
    shot = fwd

    if not target_generated:
        unresolved = True
        advisory_only = True
    elif recover_lambda:
        left, right = 0.2 * lambda_true, 2.0 * lambda_true

        def fr(lam: float) -> float:
            run = integrate_to_lambda(spin=spin, b=b_true, direction=direction, r0=r0, lambda_end=lam)
            if run["failed_reason"] is not None:
                return float("nan")
            return run["states"][-1][1] - r_target

        fl = fr(left)
        frv = fr(right)
        bracket_found = math.isfinite(fl) and math.isfinite(frv) and (fl * frv <= 0.0)
        if bracket_found:
            ok, lam_rec, it = _bisection(fr, left, right, tol=1.0e-8, max_iter=40)
            solver_converged = ok and lam_rec is not None
            solver_iterations = it
            lambda_recovered = lam_rec
            if solver_converged:
                shot = integrate_to_lambda(
                    spin=spin,
                    b=b_true,
                    direction=direction,
                    r0=r0,
                    lambda_end=lambda_recovered,
                )
            else:
                unresolved = True
                advisory_only = True
        else:
            unresolved = True
            advisory_only = True
    elif recover_b:
        assert b_bracket is not None
        left, right = b_bracket

        def fb(bv: float) -> float:
            run = integrate_to_lambda(
                spin=spin, b=bv, direction=direction, r0=r0, lambda_end=lambda_true
            )
            if run["failed_reason"] is not None:
                return float("nan")
            return run["states"][-1][2] - phi_target

        fl = fb(left)
        frv = fb(right)
        bracket_found = math.isfinite(fl) and math.isfinite(frv) and (fl * frv <= 0.0)
        if bracket_found:
            ok, b_rec, it = _bisection(fb, left, right, tol=1.0e-8, max_iter=40)
            solver_converged = ok and b_rec is not None
            solver_iterations = it
            b_recovered = b_rec
            if solver_converged:
                shot = integrate_to_lambda(
                    spin=spin, b=b_recovered, direction=direction, r0=r0, lambda_end=lambda_true
                )
            else:
                unresolved = True
                advisory_only = True
        else:
            unresolved = True
            advisory_only = True

    t_shot, r_shot, phi_shot = _endpoint(shot["states"])
    dt = t_shot - t_target
    dr = r_shot - r_target
    dphi = phi_shot - phi_target
    weighted = _weighted_residual(dt, dr, dphi)

    b_err = abs((b_recovered if b_recovered is not None else b_true) - b_true)
    lam_err = abs((lambda_recovered if lambda_recovered is not None else lambda_true) - lambda_true)
    synthetic_target_hit = (
        abs(dt) <= T_TOL and abs(dr) <= R_TOL and abs(dphi) <= PHI_TOL and weighted <= 1.0
    )
    if not solver_converged and (recover_b or recover_lambda):
        synthetic_target_hit = False

    checked = check_segment(
        states=shot["states"],
        rhs_vals=shot["rhs"],
        mass=MASS,
        spin=spin,
        impact_b=(b_recovered if b_recovered is not None else b_true),
        direction=direction,
        energy=ENERGY,
        r_plus=kerr_horizon_radius(MASS, spin),
        enforce_radial_sign=True,
        apply_schwarzschild_limit=(spin == 0.0 and b_true == 0.0),
    )
    checks = checked["checks"]
    metrics = checked["metrics"]

    checks_extra = {
        "target_was_forward_generated": target_generated,
        "no_sprinkling_pair_used": True,
        "no_global_causal_relations_decided": True,
        "bracket_found": (bracket_found if (recover_b or recover_lambda) else True),
        "solver_converged": (solver_converged if (recover_b or recover_lambda) else True),
        "synthetic_target_hit": synthetic_target_hit,
    }
    all_pass = (
        shot["failed_reason"] is None
        and checked["all_checks_pass"]
        and all(checks_extra.values())
        and (b_err <= B_TOL if recover_b and solver_converged else True)
        and (lam_err <= 1.0e-5 if recover_lambda and solver_converged else True)
    )
    if unresolved:
        all_pass = False
    synthetic_known_answer_recovered = synthetic_target_hit

    return {
        "case_id": case_id,
        "spin_a": spin,
        "M": MASS,
        "E": ENERGY,
        "direction": "outgoing" if direction > 0 else "ingoing",
        "case_type": case_type,
        "r_plus": kerr_horizon_radius(MASS, spin),
        "r0": r0,
        "lambda_true": lambda_true,
        "b_true": b_true,
        "b_recovered": b_recovered,
        "lambda_recovered": lambda_recovered,
        "bracket_found": checks_extra["bracket_found"],
        "solver_converged": checks_extra["solver_converged"],
        "solver_iterations": solver_iterations,
        "target_was_forward_generated": checks_extra["target_was_forward_generated"],
        "no_sprinkling_pair_used": checks_extra["no_sprinkling_pair_used"],
        "no_global_causal_relations_decided": checks_extra["no_global_causal_relations_decided"],
        "synthetic_target_hit": checks_extra["synthetic_target_hit"],
        "synthetic_known_answer_recovered": synthetic_known_answer_recovered,
        "endpoint_t_residual": dt,
        "endpoint_r_residual": dr,
        "endpoint_phi_residual": dphi,
        "endpoint_weighted_residual": weighted,
        "recovered_b_error": b_err if recover_b and solver_converged else None,
        "recovered_lambda_error": lam_err if recover_lambda and solver_converged else None,
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
        "schwarzschild_radial_limit_pass": checks["schwarzschild_radial_limit_pass"],
        "radial_sign_consistency_pass": checks["radial_sign_consistency_pass"],
        "advisory_only": advisory_only,
        "unresolved": unresolved,
        "all_checks_pass": all_pass,
        "_trajectory": shot["states"],
    }


def build_cases() -> list[dict[str, Any]]:
    cases = []
    # 1) Schwarzschild radial hard gate
    for direction, label, r0 in ((+1.0, "outgoing", 4.0), (-1.0, "ingoing", 8.0)):
        cases.append(
            _case(
                case_id=f"k11_schw_radial_{label}",
                spin=0.0,
                direction=direction,
                r0=r0,
                lambda_true=0.4,
                b_true=0.0,
                case_type="schwarzschild_radial_known_answer_shooting",
                recover_b=False,
                recover_lambda=True,
            )
        )
    # 2) Kerr b=0 control shooting
    for spin in (0.25, 0.5):
        r_plus = kerr_horizon_radius(MASS, spin)
        for direction, label in ((+1.0, "outgoing"), (-1.0, "ingoing")):
            cases.append(
                _case(
                    case_id=f"k11_safe_b0_spin_{spin:.2f}_{label}",
                    spin=spin,
                    direction=direction,
                    r0=r_plus + schwarz.EXTERIOR_MARGIN + (1.0 if direction > 0 else 2.0),
                    lambda_true=0.35,
                    b_true=0.0,
                    case_type="safe_b0_control_shooting",
                    recover_b=False,
                    recover_lambda=True,
                )
            )
    # 3) one-parameter b-shooting synthetic test
    for spin in (0.25, 0.5):
        r_plus = kerr_horizon_radius(MASS, spin)
        r0 = r_plus + schwarz.EXTERIOR_MARGIN + 1.3
        cases.append(
            _case(
                case_id=f"k11_bshoot_spin_{spin:.2f}_pro",
                spin=spin,
                direction=+1.0,
                r0=r0,
                lambda_true=0.25,
                b_true=+1.0,
                case_type="synthetic_b_shooting_prograde_like",
                recover_b=True,
                recover_lambda=False,
                b_bracket=(0.2, 2.2),
            )
        )
        cases.append(
            _case(
                case_id=f"k11_bshoot_spin_{spin:.2f}_retro",
                spin=spin,
                direction=+1.0,
                r0=r0,
                lambda_true=0.25,
                b_true=-1.0,
                case_type="synthetic_b_shooting_retrograde_like",
                recover_b=True,
                recover_lambda=False,
                b_bracket=(-2.2, -0.2),
            )
        )
    return cases


def write_csv(rows: list[dict[str, Any]], path: Path) -> None:
    fields = [
        "case_id", "spin_a", "M", "E", "direction", "case_type", "r_plus", "r0",
        "lambda_true", "b_true", "b_recovered", "lambda_recovered", "bracket_found",
        "solver_converged", "solver_iterations", "target_was_forward_generated",
        "no_sprinkling_pair_used", "no_global_causal_relations_decided",
        "synthetic_target_hit", "synthetic_known_answer_recovered", "endpoint_t_residual",
        "endpoint_r_residual", "endpoint_phi_residual", "endpoint_weighted_residual",
        "recovered_b_error", "recovered_lambda_error", "min_r", "max_r", "min_Delta",
        "min_R", "max_abs_null_residual", "max_abs_E_residual", "max_abs_L_residual",
        "all_points_exterior", "finite_rhs_all_steps", "finite_solution_all_steps",
        "t_monotonic_future_pass", "radial_rhs_consistency_pass", "null_condition_pass",
        "constants_consistency_pass", "schwarzschild_radial_limit_pass",
        "radial_sign_consistency_pass", "advisory_only", "unresolved", "all_checks_pass",
    ]
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k) for k in fields})


def write_json(rows: list[dict[str, Any]], summary: dict[str, Any], path: Path) -> None:
    payload = {
        "benchmark": "S4-KERR-K11 equatorial shooting sandbox",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "command": COMMAND,
        "scope_note": (
            "Synthetic known-answer shooting only. Targets are forward generated. "
            "No sprinkling event-pair usage and no production Kerr causal classifier."
        ),
        "cases": [{k: v for k, v in r.items() if not k.startswith("_")} for r in rows],
        "global_summary": summary,
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def write_md(summary: dict[str, Any], path: Path) -> None:
    lines = [
        "# S4-KERR-K11 equatorial shooting sandbox",
        "",
        "K11 is a synthetic known-answer shooting sandbox.",
        "Targets are generated by forward integration.",
        "K11 does not use sprinkling event pairs.",
        "K11 does not decide causal reachability.",
        "K11 does not implement a production Kerr causal classifier.",
        "synthetic_target_hit is not the same as physical causal reachability.",
        "S4-THERMO-001 remains a Level-A horizon/thermo guardrail, not Level-B discrete rediscovery.",
        "Unresolved or unbracketed cases are allowed and are reported explicitly.",
        "",
        f"- Total cases: {summary['total_cases']}",
        f"- Passed cases: {summary['passed_cases']}",
        f"- Failed cases: {summary['failed_cases']}",
        f"- Advisory cases: {summary['advisory_cases']}",
        f"- Unresolved cases: {summary['unresolved_cases']}",
        f"- Synthetic targets generated: {summary['synthetic_targets_generated']}",
        f"- Synthetic targets hit: {summary['synthetic_targets_hit']}",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_png(rows: list[dict[str, Any]], path: Path) -> None:
    fig, axs = plt.subplots(2, 2, figsize=(11, 8))
    fig.suptitle("S4-KERR-K11 equatorial shooting sandbox")

    ax = axs[0, 0]
    x = list(range(len(rows)))
    wr = [r["endpoint_weighted_residual"] for r in rows]
    ax.semilogy(x, wr, "o-")
    ax.set_title("Endpoint weighted residuals")
    ax.set_xlabel("case index")

    ax = axs[0, 1]
    b_rows = [r for r in rows if "b_shooting" in r["case_type"] and r["solver_converged"]]
    if b_rows:
        ax.plot([r["b_true"] for r in b_rows], [r["b_recovered"] for r in b_rows], "o")
    ax.set_xlabel("b_true")
    ax.set_ylabel("b_recovered")
    ax.set_title("True b vs recovered b")

    ax = axs[1, 0]
    phi_res = [abs(r["endpoint_phi_residual"]) for r in rows]
    ax.semilogy(x, phi_res, "o-")
    ax.set_title("|endpoint phi residual|")
    ax.set_xlabel("case index")

    ax = axs[1, 1]
    for r in rows[:6]:
        lam = [i * H for i in range(len(r["_trajectory"]))]
        rr = [s[1] for s in r["_trajectory"]]
        ax.plot(lam, rr, label=r["case_id"])
    ax.set_xlabel("lambda")
    ax.set_ylabel("r(lambda)")
    ax.set_title("Representative recovered trajectories")
    ax.legend(fontsize=6)

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def run() -> int:
    rows = build_cases()
    non_advisory_non_unresolved = [r for r in rows if not r["advisory_only"] and not r["unresolved"]]
    if any(not r["all_checks_pass"] for r in non_advisory_non_unresolved):
        print("STOP: non-advisory resolved case failed.")
        return 1
    summary = {
        "total_cases": len(rows),
        "passed_cases": sum(1 for r in rows if r["all_checks_pass"]),
        "failed_cases": sum(1 for r in rows if (not r["all_checks_pass"]) and (not r["unresolved"])),
        "advisory_cases": sum(1 for r in rows if r["advisory_only"]),
        "unresolved_cases": sum(1 for r in rows if r["unresolved"]),
        "synthetic_targets_generated": sum(1 for r in rows if r["target_was_forward_generated"]),
        "synthetic_targets_hit": sum(1 for r in rows if r["synthetic_target_hit"]),
        "global_true_relations": 0,
        "global_false_relations": 0,
        "global_undecided_pairs": N_EVENTS * (N_EVENTS - 1) // 2,
        "all_checks_pass": all(
            r["all_checks_pass"] for r in rows if (not r["advisory_only"] and not r["unresolved"])
        ),
    }

    csv_path = ARTIFACT_DIR / f"{OUT_PREFIX}.csv"
    json_path = ARTIFACT_DIR / f"{OUT_PREFIX}.json"
    md_path = ARTIFACT_DIR / f"{OUT_PREFIX}.md"
    png_path = ARTIFACT_DIR / f"{OUT_PREFIX}.png"
    write_csv(rows, csv_path)
    write_json(rows, summary, json_path)
    write_md(summary, md_path)
    write_png(rows, png_path)
    print(
        f"cases={summary['total_cases']} passed={summary['passed_cases']} "
        f"unresolved={summary['unresolved_cases']} hits={summary['synthetic_targets_hit']}"
    )
    return 0


def main() -> None:
    raise SystemExit(run())


if __name__ == "__main__":
    main()
