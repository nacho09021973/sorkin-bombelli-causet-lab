#!/usr/bin/env python3
"""Kerr K1 scaffold invariant/control audit.

This audit deliberately does not implement Kerr causal physics.  It extends
the L0 bookkeeping scaffold to a small spin sweep:

* a=0 must exactly match the frozen Schwarzschild exterior benchmark subset;
* a>0 samples only outside r_+=M+sqrt(M^2-a^2) plus the existing margin;
* a>0 leaves every pair undecided until a justified Kerr causal solver exists.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from explore.sorkin4_kerr_benchmark import audit_kerr_l0_scaffold as l0  # noqa: E402
from explore.sorkin4_kerr_benchmark import run_kerr_minimal_benchmark as kerr  # noqa: E402
from explore.sorkin4_schwarzschild_benchmark import run_schwarzschild_minimal_benchmark as schwarz  # noqa: E402


OUT_DIR = Path(__file__).resolve().parent
OUT_PREFIX = "kerr_k1_control_spin_sweep_n12_seed1959"
COMMAND = "python3 explore/sorkin4_kerr_benchmark/audit_kerr_k1_control.py"

DEFAULT_N = 12
DEFAULT_SEED = 1959
DEFAULT_MASS = 1.0
DEFAULT_MARGIN = schwarz.EXTERIOR_MARGIN
DEFAULT_SPINS = (0.0, 0.25, 0.5, 0.75)
FLOAT_TOL = 1.0e-12


def _events_to_rows(events: list[kerr.Event]) -> list[dict[str, float | int]]:
    return [asdict(event) for event in events]


def _run_kerr_case(
    n: int,
    seed: int,
    mass: float,
    spin: float,
    margin: float,
) -> tuple[list[kerr.Event], list[list[bool]], list[list[Optional[bool]]], dict[str, object]]:
    r_plus = kerr.kerr_horizon_radius(mass, spin)
    r_min = r_plus + margin
    events = kerr.generate_exterior_events(n, seed, r_min)
    matrix, states = kerr.build_relation_states(events, mass, spin, "scaffold")
    summary = kerr.summarize_case(
        events,
        matrix,
        states,
        n=n,
        seed=seed,
        mass=mass,
        spin=spin,
        r_plus=r_plus,
        r_min=r_min,
        margin=margin,
        kerr_mode="scaffold",
    )
    return events, matrix, states, summary


def _upper_triangle_all_undecided(states: list[list[Optional[bool]]]) -> bool:
    n = len(states)
    return all(states[i][j] is None for i in range(n - 1) for j in range(i + 1, n))


def _strictly_exterior(events: list[kerr.Event], r_min: float) -> bool:
    return all(event.r > r_min for event in events)


def _case_row(
    summary: dict[str, object],
    checks: dict[str, bool],
    all_checks_pass: bool,
) -> dict[str, object]:
    spin = float(summary["a"])
    return {
        "benchmark": "S4-K1 Kerr scaffold invariant control audit",
        "status": "pass" if all_checks_pass else "fail",
        "generated_at_utc": summary["generated_at_utc"],
        "command": COMMAND,
        "N": summary["N"],
        "seed": summary["seed"],
        "M": summary["M"],
        "a": spin,
        "r_plus": summary["r_plus"],
        "r_min_margin": summary["r_min_margin"],
        "r_min": summary["r_min"],
        "possible_pairs": summary["possible_pairs"],
        "true_relations": summary["true_relations"],
        "false_relations": summary["false_relations"],
        "undecided_pairs": summary["undecided_pairs"],
        "decided_pairs": summary["decided_pairs"],
        "case_checks_pass": all(checks.values()),
        "all_checks_pass": all_checks_pass,
        "case_role": "schwarzschild_control" if spin == 0.0 else "kerr_scaffold_undecided_control",
        "warning": (
            "a=0 Schwarzschild exterior regression/control"
            if spin == 0.0
            else "Kerr a>0 remains scaffold only; undecided pairs are not non-relations."
        ),
    }


def run_audit(
    n: int = DEFAULT_N,
    seed: int = DEFAULT_SEED,
    mass: float = DEFAULT_MASS,
    margin: float = DEFAULT_MARGIN,
    spins: tuple[float, ...] = DEFAULT_SPINS,
) -> dict[str, object]:
    if mass != 1.0:
        raise ValueError("K1 control audit is fixed to M=1")
    if margin != schwarz.EXTERIOR_MARGIN:
        raise ValueError("K1 a=0 control uses the frozen Schwarzschild exterior margin")
    if spins != DEFAULT_SPINS:
        raise ValueError("K1 control audit uses the frozen spin sweep (0.0, 0.25, 0.5, 0.75)")
    if any(abs(spin) >= mass for spin in spins):
        raise ValueError("K1 spin sweep requires |a| < M")

    sch_events, sch_matrix, sch_states, sch_summary, _debug_rows = schwarz.run_case(
        n,
        seed,
        enable_shooting=False,
    )
    sch_event_rows = _events_to_rows(sch_events)
    possible_pairs = n * (n - 1) // 2

    case_payloads: list[dict[str, object]] = []
    rows: list[dict[str, object]] = []

    for spin in spins:
        events, matrix, states, summary = _run_kerr_case(n, seed, mass, spin, margin)
        event_rows = _events_to_rows(events)
        r_plus_formula = mass + math.sqrt(mass * mass - spin * spin)
        r_min = float(summary["r_min"])

        common_checks = {
            "spin_subextremal": abs(spin) < mass,
            "r_plus_matches_formula": abs(float(summary["r_plus"]) - r_plus_formula) <= FLOAT_TOL,
            "r_min_is_r_plus_plus_margin": abs(r_min - (r_plus_formula + margin)) <= FLOAT_TOL,
            "events_strictly_exterior": _strictly_exterior(events, r_min),
            "possible_pairs_match": summary["possible_pairs"] == possible_pairs,
        }

        if spin == 0.0:
            checks = {
                **common_checks,
                "r_plus_equals_2M": abs(float(summary["r_plus"]) - 2.0 * mass) <= FLOAT_TOL,
                "events_match_schwarzschild": l0._event_rows_match(event_rows, sch_event_rows),
                "causal_matrix_matches_schwarzschild": matrix == sch_matrix,
                "relation_states_match_schwarzschild": states == sch_states,
                "true_relations_match_schwarzschild": summary["true_relations"] == sch_summary["true_relations"],
                "false_relations_match_schwarzschild": summary["false_relations"] == sch_summary["false_relations"],
                "undecided_pairs_match_schwarzschild": summary["undecided_pairs"] == sch_summary["undecided_pairs"],
            }
        else:
            checks = {
                **common_checks,
                "status_is_scaffold_undecided": summary["status"] == "kerr_scaffold_pairs_undecided",
                "true_relations_zero": summary["true_relations"] == 0,
                "false_relations_zero": summary["false_relations"] == 0,
                "undecided_pairs_all": summary["undecided_pairs"] == possible_pairs,
                "decided_pairs_zero": summary["decided_pairs"] == 0,
                "causal_matrix_empty": not any(any(row) for row in matrix),
                "relation_states_upper_triangle_undecided": _upper_triangle_all_undecided(states),
            }

        case_payloads.append(
            {
                "spin": spin,
                "summary": summary,
                "checks": checks,
                "events": event_rows,
                "relation_states": states,
            }
        )
        rows.append(_case_row(summary, checks, all_checks_pass=False))

    all_checks_pass = all(all(case["checks"].values()) for case in case_payloads if isinstance(case["checks"], dict))
    for row in rows:
        row["status"] = "pass" if all_checks_pass else "fail"
        row["all_checks_pass"] = all_checks_pass

    aggregate = {
        "benchmark": "S4-K1 Kerr scaffold invariant control audit",
        "status": "pass" if all_checks_pass else "fail",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "command": COMMAND,
        "N": n,
        "seed": seed,
        "M": mass,
        "spins": list(spins),
        "r_min_margin": margin,
        "possible_pairs": possible_pairs,
        "all_checks_pass": all_checks_pass,
        "a0_exact_schwarzschild_control": bool(case_payloads[0]["checks"]["relation_states_match_schwarzschild"]),
        "positive_spin_cases_all_undecided": all(
            case["summary"]["undecided_pairs"] == possible_pairs
            and case["summary"]["true_relations"] == 0
            and case["summary"]["false_relations"] == 0
            for case in case_payloads
            if float(case["spin"]) > 0.0
        ),
        "warning": "K1 is a scaffold invariant/control audit only; it does not validate Kerr causal relations.",
    }

    return {
        "aggregate": aggregate,
        "rows": rows,
        "cases": case_payloads,
        "schwarzschild_reference": {
            "summary": sch_summary,
            "events": sch_event_rows,
            "relation_states": sch_states,
        },
        "notes": [
            "This is a conservative scaffold invariant/control audit, not a Kerr causal solver.",
            "M is fixed to 1 and all swept spins satisfy |a| < M.",
            "Each case samples only r > r_+ + margin.",
            "The a=0 case is an exact regression against the frozen Schwarzschild exterior benchmark.",
            "For a>0, every unordered pair remains undecided by construction.",
            "Undecided pairs are unknown pairs, not decided non-relations.",
        ],
    }


def _fmt(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return f"{value:.12g}"
    if isinstance(value, list):
        return ";".join(str(item) for item in value)
    return str(value)


def write_outputs(payload: dict[str, object], out_prefix: str = OUT_PREFIX) -> tuple[Path, Path, Path]:
    csv_path = OUT_DIR / f"{out_prefix}.csv"
    json_path = OUT_DIR / f"{out_prefix}.json"
    md_path = OUT_DIR / f"{out_prefix}.md"

    rows = payload["rows"]
    aggregate = payload["aggregate"]
    if not isinstance(rows, list) or not rows:
        raise TypeError("rows payload must be a non-empty list")
    if not isinstance(aggregate, dict):
        raise TypeError("aggregate payload must be a dict")

    headers = list(rows[0])
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(headers)
        for row in rows:
            writer.writerow([_fmt(row[header]) for header in headers])

    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    md_lines = [
        "# S4-K1 Kerr scaffold invariant control audit",
        "",
        "Scope: conservative scaffold invariant/control audit only. This does not implement Kerr causal physics.",
        "",
        f"- `M`: `{aggregate['M']}`",
        f"- `spins`: `{aggregate['spins']}`",
        f"- `N`: `{aggregate['N']}`",
        f"- `seed`: `{aggregate['seed']}`",
        f"- `possible_pairs`: `{aggregate['possible_pairs']}`",
        f"- `all_checks_pass`: `{aggregate['all_checks_pass']}`",
        "",
        "Controls:",
        "",
        "- `a=0.0` must exactly match the frozen Schwarzschild exterior benchmark.",
        "- `a>0` uses only exterior points with `r > r_+ + margin`.",
        "- `a>0` leaves all unordered pairs undecided: true `0`, false `0`, undecided `N*(N-1)/2`.",
        "",
        "Per-spin counts:",
        "",
        "| a | r_plus | r_min | true | false | undecided | checks |",
        "|---:|---:|---:|---:|---:|---:|:---:|",
    ]
    for row in rows:
        md_lines.append(
            "| {a:.2f} | {r_plus:.12g} | {r_min:.12g} | {true_relations} | "
            "{false_relations} | {undecided_pairs} | {case_checks_pass} |".format(**row)
        )
    md_lines.extend(
        [
            "",
            "Interpretation: undecided means not decided by this scaffold. It is not a Kerr non-relation.",
            "",
        ]
    )
    md_path.write_text("\n".join(md_lines), encoding="utf-8")
    return csv_path, json_path, md_path


def main() -> None:
    payload = run_audit()
    csv_path, json_path, md_path = write_outputs(payload)
    aggregate = payload["aggregate"]
    if not isinstance(aggregate, dict):
        raise TypeError("aggregate payload must be a dict")
    print("S4-K1 Kerr scaffold invariant control audit")
    print(f"all_checks_pass={aggregate['all_checks_pass']}")
    for row in payload["rows"]:
        print(
            f"a={row['a']} r_plus={row['r_plus']:.12g} "
            f"true={row['true_relations']} false={row['false_relations']} "
            f"undecided={row['undecided_pairs']} checks={row['case_checks_pass']}"
        )
    print(f"wrote {csv_path}")
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")


if __name__ == "__main__":
    main()
