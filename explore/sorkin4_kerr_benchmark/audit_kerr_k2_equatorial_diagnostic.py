#!/usr/bin/env python3
"""Kerr K2 equatorial diagnostic scaffold.

This audit does not implement Kerr causal decisions.  It only freezes a small
equatorial diagnostic sweep that can later inform a justified solver:

* M=1, theta=pi/2, a in (0.0, 0.25, 0.5, 0.75);
* points are sampled only outside r_+ plus the existing margin;
* a=0 remains the Schwarzschild control;
* a>0 leaves all causal relation states undecided.
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
OUT_PREFIX = "kerr_k2_equatorial_diagnostic_n12_seed1959"
COMMAND = "python3 explore/sorkin4_kerr_benchmark/audit_kerr_k2_equatorial_diagnostic.py"

DEFAULT_N = 12
DEFAULT_SEED = 1959
DEFAULT_MASS = 1.0
DEFAULT_MARGIN = schwarz.EXTERIOR_MARGIN
DEFAULT_SPINS = (0.0, 0.25, 0.5, 0.75)
FLOAT_TOL = 1.0e-12


def _events_to_rows(events: list[kerr.Event]) -> list[dict[str, float | int]]:
    return [asdict(event) for event in events]


def _run_equatorial_case(
    n: int,
    seed: int,
    mass: float,
    spin: float,
    margin: float,
) -> tuple[list[kerr.Event], list[list[bool]], list[list[Optional[bool]]], dict[str, object]]:
    r_plus = kerr.kerr_horizon_radius(mass, spin)
    r_min = r_plus + margin
    events = kerr.generate_exterior_events(n, seed, r_min, equatorial=True)
    matrix, states = kerr.build_relation_states(events, mass, spin, "equatorial_scaffold")
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
        kerr_mode="equatorial_scaffold",
    )
    return events, matrix, states, summary


def _positive_spin_pair_diagnostics(events: list[kerr.Event], spin: float) -> dict[str, float | int]:
    prograde: list[float] = []
    retrograde: list[float] = []
    if spin == 0.0:
        return {
            "prograde_pair_count": 0,
            "retrograde_pair_count": 0,
            "mean_delta_phi_prograde": 0.0,
            "mean_delta_phi_retrograde": 0.0,
        }

    for i in range(len(events) - 1):
        for j in range(i + 1, len(events)):
            delta_phi = kerr.signed_delta_phi(events[i].phi, events[j].phi)
            if abs(delta_phi) <= kerr.RADIAL_PHI_EPS:
                continue
            if delta_phi * spin > 0.0:
                prograde.append(delta_phi)
            else:
                retrograde.append(delta_phi)

    return {
        "prograde_pair_count": len(prograde),
        "retrograde_pair_count": len(retrograde),
        "mean_delta_phi_prograde": sum(prograde) / len(prograde) if prograde else 0.0,
        "mean_delta_phi_retrograde": sum(retrograde) / len(retrograde) if retrograde else 0.0,
    }


def _all_pairs_undecided(states: list[list[Optional[bool]]]) -> bool:
    n = len(states)
    return all(states[i][j] is None for i in range(n - 1) for j in range(i + 1, n))


def _row_for_case(
    summary: dict[str, object],
    events: list[kerr.Event],
    states: list[list[Optional[bool]]],
    matrix: list[list[bool]],
    mass: float,
    checks: dict[str, bool],
    all_checks_pass: bool,
) -> dict[str, object]:
    spin = float(summary["a"])
    r_ergosphere_equatorial = 2.0 * mass
    pair_diagnostics = _positive_spin_pair_diagnostics(events, spin)
    r_min_observed = min(event.r for event in events)
    return {
        "benchmark": "S4-K2 Kerr equatorial diagnostic scaffold",
        "status": "pass" if all_checks_pass else "fail",
        "generated_at_utc": summary["generated_at_utc"],
        "command": COMMAND,
        "N": summary["N"],
        "seed": summary["seed"],
        "M": mass,
        "a": spin,
        "theta": math.pi / 2.0,
        "r_plus": summary["r_plus"],
        "r_ergosphere_equatorial": r_ergosphere_equatorial,
        "r_min_margin": summary["r_min_margin"],
        "r_min": summary["r_min"],
        "r_min_observed": r_min_observed,
        "r_min_observed_gt_r_plus_plus_margin": r_min_observed > float(summary["r_min"]),
        "outside_equatorial_ergosphere_count": sum(1 for event in events if event.r > r_ergosphere_equatorial),
        "inside_equatorial_ergosphere_count": sum(1 for event in events if event.r <= r_ergosphere_equatorial),
        "prograde_pair_count": pair_diagnostics["prograde_pair_count"],
        "retrograde_pair_count": pair_diagnostics["retrograde_pair_count"],
        "mean_delta_phi_prograde": pair_diagnostics["mean_delta_phi_prograde"],
        "mean_delta_phi_retrograde": pair_diagnostics["mean_delta_phi_retrograde"],
        "possible_pairs": summary["possible_pairs"],
        "true_relations": summary["true_relations"],
        "false_relations": summary["false_relations"],
        "undecided_pairs": summary["undecided_pairs"],
        "decided_pairs": summary["decided_pairs"],
        "all_pairs_undecided": _all_pairs_undecided(states),
        "causal_matrix_empty": not any(any(row) for row in matrix),
        "case_checks_pass": all(checks.values()),
        "all_checks_pass": all_checks_pass,
        "warning": "Equatorial kinematic diagnostic only; no Kerr causal relations are decided.",
    }


def run_audit(
    n: int = DEFAULT_N,
    seed: int = DEFAULT_SEED,
    mass: float = DEFAULT_MASS,
    margin: float = DEFAULT_MARGIN,
    spins: tuple[float, ...] = DEFAULT_SPINS,
) -> dict[str, object]:
    if mass != 1.0:
        raise ValueError("K2 equatorial diagnostic is fixed to M=1")
    if margin != schwarz.EXTERIOR_MARGIN:
        raise ValueError("K2 uses the frozen Schwarzschild exterior margin")
    if spins != DEFAULT_SPINS:
        raise ValueError("K2 uses the frozen spin sweep (0.0, 0.25, 0.5, 0.75)")
    if any(abs(spin) >= mass for spin in spins):
        raise ValueError("K2 spin sweep requires |a| < M")

    possible_pairs = n * (n - 1) // 2
    case_payloads: list[dict[str, object]] = []
    rows: list[dict[str, object]] = []

    for spin in spins:
        events, matrix, states, summary = _run_equatorial_case(n, seed, mass, spin, margin)
        r_plus_formula = mass + math.sqrt(mass * mass - spin * spin)
        r_min = float(summary["r_min"])
        common_checks = {
            "spin_subextremal": abs(spin) < mass,
            "theta_all_equatorial": all(abs(event.theta - math.pi / 2.0) <= FLOAT_TOL for event in events),
            "r_plus_matches_formula": abs(float(summary["r_plus"]) - r_plus_formula) <= FLOAT_TOL,
            "r_min_is_r_plus_plus_margin": abs(r_min - (r_plus_formula + margin)) <= FLOAT_TOL,
            "r_min_observed_exterior": min(event.r for event in events) > r_min,
            "possible_pairs_match": summary["possible_pairs"] == possible_pairs,
        }
        if spin == 0.0:
            checks = {
                **common_checks,
                "status_is_schwarzschild_control": summary["status"] == "a0_schwarzschild_regression",
                "a0_uses_schwarzschild_subset": summary["decided_pairs"] > 0,
            }
        else:
            checks = {
                **common_checks,
                "status_is_equatorial_scaffold": summary["status"] == "kerr_equatorial_scaffold",
                "true_relations_zero": summary["true_relations"] == 0,
                "false_relations_zero": summary["false_relations"] == 0,
                "undecided_pairs_all": summary["undecided_pairs"] == possible_pairs,
                "decided_pairs_zero": summary["decided_pairs"] == 0,
                "all_pairs_undecided": _all_pairs_undecided(states),
                "causal_matrix_empty": not any(any(row) for row in matrix),
            }
        case_payloads.append(
            {
                "spin": spin,
                "summary": summary,
                "checks": checks,
                "events": _events_to_rows(events),
                "relation_states": states,
            }
        )
        rows.append(
            _row_for_case(
                summary,
                events,
                states,
                matrix,
                mass,
                checks,
                all_checks_pass=False,
            )
        )

    all_checks_pass = all(all(case["checks"].values()) for case in case_payloads if isinstance(case["checks"], dict))
    for row in rows:
        row["status"] = "pass" if all_checks_pass else "fail"
        row["all_checks_pass"] = all_checks_pass

    aggregate = {
        "benchmark": "S4-K2 Kerr equatorial diagnostic scaffold",
        "status": "pass" if all_checks_pass else "fail",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "command": COMMAND,
        "N": n,
        "seed": seed,
        "M": mass,
        "spins": list(spins),
        "theta": math.pi / 2.0,
        "r_min_margin": margin,
        "possible_pairs": possible_pairs,
        "all_checks_pass": all_checks_pass,
        "positive_spin_cases_all_undecided": all(
            row["true_relations"] == 0
            and row["false_relations"] == 0
            and row["undecided_pairs"] == possible_pairs
            for row in rows
            if float(row["a"]) > 0.0
        ),
        "warning": "K2 is equatorial kinematic scaffolding only; it does not validate Kerr causality.",
    }

    return {
        "aggregate": aggregate,
        "rows": rows,
        "cases": case_payloads,
        "notes": [
            "This audit is a Kerr equatorial diagnostic scaffold, not a causal solver.",
            "All events are generated at theta=pi/2.",
            "The equatorial ergosphere radius is recorded as r=2M.",
            "Prograde/retrograde counts use the signed shortest delta_phi relative to spin sign.",
            "For a=0, relation decisions are the Schwarzschild control subset.",
            "For a>0, all relation states remain undecided.",
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
        "# S4-K2 Kerr equatorial diagnostic scaffold",
        "",
        "Scope: equatorial kinematic diagnostic only. This does not decide Kerr causal relations.",
        "",
        f"- `M`: `{aggregate['M']}`",
        f"- `theta`: `{aggregate['theta']}`",
        f"- `spins`: `{aggregate['spins']}`",
        f"- `N`: `{aggregate['N']}`",
        f"- `seed`: `{aggregate['seed']}`",
        f"- `all_checks_pass`: `{aggregate['all_checks_pass']}`",
        "",
        "| a | r_plus | r_ergosphere_eq | r_min_observed | inside_ergo | outside_ergo | prograde | retrograde | undecided |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        md_lines.append(
            "| {a:.2f} | {r_plus:.12g} | {r_ergosphere_equatorial:.12g} | "
            "{r_min_observed:.12g} | {inside_equatorial_ergosphere_count} | "
            "{outside_equatorial_ergosphere_count} | {prograde_pair_count} | "
            "{retrograde_pair_count} | {undecided_pairs} |".format(**row)
        )
    md_lines.extend(
        [
            "",
            "Interpretation: these are scaffold diagnostics only. Undecided pairs are unknown, not non-relations.",
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
    print("S4-K2 Kerr equatorial diagnostic scaffold")
    print(f"all_checks_pass={aggregate['all_checks_pass']}")
    for row in payload["rows"]:
        print(
            f"a={row['a']} r_plus={row['r_plus']:.12g} "
            f"inside_ergo={row['inside_equatorial_ergosphere_count']} "
            f"outside_ergo={row['outside_equatorial_ergosphere_count']} "
            f"prograde={row['prograde_pair_count']} retrograde={row['retrograde_pair_count']} "
            f"undecided={row['undecided_pairs']}"
        )
    print(f"wrote {csv_path}")
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")


if __name__ == "__main__":
    main()
