#!/usr/bin/env python3
"""Kerr L0 scaffold audit with an a=0 Schwarzschild control.

This audit deliberately does not implement Kerr causal decisions.  It freezes
two bookkeeping guarantees:

* a=0 uses the same generated events and relation states as the Schwarzschild
  exterior benchmark subset;
* a!=0 scaffold mode leaves every pair undecided.
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

from explore.sorkin4_kerr_benchmark import run_kerr_minimal_benchmark as kerr  # noqa: E402
from explore.sorkin4_schwarzschild_benchmark import run_schwarzschild_minimal_benchmark as schwarz  # noqa: E402


OUT_DIR = Path(__file__).resolve().parent
OUT_PREFIX = "kerr_l0_scaffold_control_a0_a0p5_n12_seed1959"
COMMAND = "python3 explore/sorkin4_kerr_benchmark/audit_kerr_l0_scaffold.py"

DEFAULT_N = 12
DEFAULT_SEED = 1959
DEFAULT_MASS = 1.0
DEFAULT_SPIN = 0.5
DEFAULT_MARGIN = 0.35
FLOAT_TOL = 1.0e-12


def _events_to_rows(events: list[kerr.Event]) -> list[dict[str, float | int]]:
    return [asdict(event) for event in events]


def _event_rows_match(
    left: list[dict[str, float | int]],
    right: list[dict[str, float | int]],
    tol: float = FLOAT_TOL,
) -> bool:
    if len(left) != len(right):
        return False
    for p, q in zip(left, right):
        if p["index"] != q["index"]:
            return False
        for key in ("t", "r", "theta", "phi"):
            if abs(float(p[key]) - float(q[key])) > tol:
                return False
    return True


def _states_match(
    left: list[list[Optional[bool]]],
    right: list[list[Optional[bool]]],
) -> bool:
    return left == right


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


def run_audit(
    n: int = DEFAULT_N,
    seed: int = DEFAULT_SEED,
    mass: float = DEFAULT_MASS,
    spin: float = DEFAULT_SPIN,
    margin: float = DEFAULT_MARGIN,
) -> dict[str, object]:
    if mass != 1.0 or margin != schwarz.EXTERIOR_MARGIN:
        raise ValueError("L0 a=0 control is defined against the default Schwarzschild benchmark")
    if not (0.0 < abs(spin) < mass):
        raise ValueError("L0 spin case requires 0 < |a| < M")

    a0_events, a0_matrix, a0_states, a0_summary = _run_kerr_case(n, seed, mass, 0.0, margin)
    sch_events, sch_matrix, sch_states, sch_summary, _debug_rows = schwarz.run_case(
        n,
        seed,
        enable_shooting=False,
    )
    spin_events, spin_matrix, spin_states, spin_summary = _run_kerr_case(n, seed, mass, spin, margin)

    a0_event_rows = _events_to_rows(a0_events)
    sch_event_rows = _events_to_rows(sch_events)
    spin_event_rows = _events_to_rows(spin_events)
    possible_pairs = n * (n - 1) // 2

    a0_checks = {
        "r_plus_equals_2M": abs(float(a0_summary["r_plus"]) - 2.0 * mass) <= FLOAT_TOL,
        "events_match_schwarzschild": _event_rows_match(a0_event_rows, sch_event_rows),
        "causal_matrix_matches_schwarzschild": a0_matrix == sch_matrix,
        "relation_states_match_schwarzschild": _states_match(a0_states, sch_states),
        "true_relations_match": a0_summary["true_relations"] == sch_summary["true_relations"],
        "false_relations_match": a0_summary["false_relations"] == sch_summary["false_relations"],
        "undecided_pairs_match": a0_summary["undecided_pairs"] == sch_summary["undecided_pairs"],
    }
    spin_checks = {
        "r_plus_matches_formula": abs(float(spin_summary["r_plus"]) - (mass + math.sqrt(mass * mass - spin * spin)))
        <= FLOAT_TOL,
        "all_pairs_undecided": spin_summary["undecided_pairs"] == possible_pairs,
        "no_true_relations": spin_summary["true_relations"] == 0,
        "no_false_relations": spin_summary["false_relations"] == 0,
        "no_decided_pairs": spin_summary["decided_pairs"] == 0,
        "matrix_empty": not any(any(row) for row in spin_matrix),
        "relation_states_upper_triangle_undecided": all(
            spin_states[i][j] is None for i in range(n - 1) for j in range(i + 1, n)
        ),
    }

    a0_checks_pass = all(a0_checks.values())
    spin_checks_pass = all(spin_checks.values())
    summary = {
        "benchmark": "S4-K0 Kerr L0 scaffold control audit",
        "status": "pass" if a0_checks_pass and spin_checks_pass else "fail",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "command": COMMAND,
        "N": n,
        "seed": seed,
        "M": mass,
        "a0": 0.0,
        "a_scaffold": spin,
        "r_min_margin": margin,
        "possible_pairs": possible_pairs,
        "a0_status": a0_summary["status"],
        "a0_true_relations": a0_summary["true_relations"],
        "a0_false_relations": a0_summary["false_relations"],
        "a0_undecided_pairs": a0_summary["undecided_pairs"],
        "a0_decided_pairs": a0_summary["decided_pairs"],
        "a0_events_match_schwarzschild": a0_checks["events_match_schwarzschild"],
        "a0_relation_states_match_schwarzschild": a0_checks["relation_states_match_schwarzschild"],
        "a0_checks_pass": a0_checks_pass,
        "spin_status": spin_summary["status"],
        "spin_r_plus": spin_summary["r_plus"],
        "spin_r_min": spin_summary["r_min"],
        "spin_true_relations": spin_summary["true_relations"],
        "spin_false_relations": spin_summary["false_relations"],
        "spin_undecided_pairs": spin_summary["undecided_pairs"],
        "spin_decided_pairs": spin_summary["decided_pairs"],
        "spin_checks_pass": spin_checks_pass,
        "all_checks_pass": a0_checks_pass and spin_checks_pass,
        "warning": "Kerr a!=0 remains scaffold only; undecided pairs are not non-relations.",
    }

    return {
        "summary": summary,
        "a0_checks": a0_checks,
        "spin_checks": spin_checks,
        "cases": {
            "kerr_a0_control": {
                "summary": a0_summary,
                "events": a0_event_rows,
                "relation_states": a0_states,
            },
            "schwarzschild_reference": {
                "summary": sch_summary,
                "events": sch_event_rows,
                "relation_states": sch_states,
            },
            "kerr_spin_scaffold": {
                "summary": spin_summary,
                "events": spin_event_rows,
                "relation_states": spin_states,
            },
        },
        "notes": [
            "This audit is a scaffold/control check, not a Kerr causal solver.",
            "The a=0 case must remain exactly aligned with the existing Schwarzschild exterior subset.",
            "The a!=0 scaffold intentionally leaves every pair undecided.",
            "Undecided relation states are unknown pairs, not decided non-relations.",
        ],
    }


def _fmt(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return f"{value:.12g}"
    return str(value)


def write_outputs(payload: dict[str, object], out_prefix: str = OUT_PREFIX) -> tuple[Path, Path, Path]:
    csv_path = OUT_DIR / f"{out_prefix}.csv"
    json_path = OUT_DIR / f"{out_prefix}.json"
    md_path = OUT_DIR / f"{out_prefix}.md"
    summary = payload["summary"]
    if not isinstance(summary, dict):
        raise TypeError("summary payload must be a dict")

    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(list(summary))
        writer.writerow([_fmt(summary[key]) for key in summary])

    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    md_lines = [
        "# S4-K0 Kerr L0 scaffold control audit",
        "",
        "Scope: this is a bookkeeping scaffold, not a Kerr causal-relation solver.",
        "",
        "Checks frozen by this artifact:",
        "",
        f"- `a=0`: exact control against the Schwarzschild exterior subset for N={summary['N']}, seed={summary['seed']}.",
        f"- `a={summary['a_scaffold']}`: scaffold mode leaves all {summary['possible_pairs']} unordered pairs undecided.",
        "- Undecided pairs are unknown pairs, not decided non-relations.",
        "",
        "Result:",
        "",
        f"- `a0_checks_pass`: `{summary['a0_checks_pass']}`",
        f"- `spin_checks_pass`: `{summary['spin_checks_pass']}`",
        f"- `all_checks_pass`: `{summary['all_checks_pass']}`",
        "",
        "Counts:",
        "",
        (
            f"- `a=0`: true `{summary['a0_true_relations']}`, false `{summary['a0_false_relations']}`, "
            f"undecided `{summary['a0_undecided_pairs']}`"
        ),
        (
            f"- `a={summary['a_scaffold']}`: true `{summary['spin_true_relations']}`, "
            f"false `{summary['spin_false_relations']}`, undecided `{summary['spin_undecided_pairs']}`"
        ),
        "",
        "Next admissible Kerr step: add only diagnostics with an explicit `a=0` regression gate before making any physical claim.",
        "",
    ]
    md_path.write_text("\n".join(md_lines), encoding="utf-8")
    return csv_path, json_path, md_path


def main() -> None:
    payload = run_audit()
    csv_path, json_path, md_path = write_outputs(payload)
    summary = payload["summary"]
    if not isinstance(summary, dict):
        raise TypeError("summary payload must be a dict")
    print("S4-K0 Kerr L0 scaffold control audit")
    print(f"all_checks_pass={summary['all_checks_pass']}")
    print(
        f"a0 true={summary['a0_true_relations']} false={summary['a0_false_relations']} "
        f"undecided={summary['a0_undecided_pairs']}"
    )
    print(
        f"a={summary['a_scaffold']} true={summary['spin_true_relations']} "
        f"false={summary['spin_false_relations']} undecided={summary['spin_undecided_pairs']}"
    )
    print(f"wrote {csv_path}")
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")


if __name__ == "__main__":
    main()
