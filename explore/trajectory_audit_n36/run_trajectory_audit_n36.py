#!/usr/bin/env python3
"""N=36 within-trajectory energy/causal diagnostic for SORKIN-2.

Exploratory only: fixed known-truth Minkowski case and optimizer seed,
recording per-block energy and causal-order metrics via the annealer's
read-only block_callback hook.
"""

from __future__ import annotations

import contextlib
import csv
import io
import math
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import cones  # noqa: E402
import validation_suite as vs  # noqa: E402


OUT_DIR = Path(__file__).resolve().parent
CSV_PATH = OUT_DIR / "trajectory_audit_n36.csv"
MD_PATH = OUT_DIR / "trajectory_audit_n36.md"
SVG_PATH = OUT_DIR / "trajectory_audit_n36.svg"
COMMAND = "python3 explore/trajectory_audit_n36/run_trajectory_audit_n36.py"

FAMILY = "minkowski"
D_SPACETIME = 2
N = 36
CASE_SEED = 1959
OPTIMIZER_SEED = 1987
TARGET_DIM = D_SPACETIME - 1
INITIAL_TEMP = 100.0
GAMMA = 0.50
H = 1
BACKEND = "cpu"

BUDGETS = (
    {"budget_label": "short_10_10_4", "warmup_limit": 10, "anneal_limit": 10, "max_data": 4},
    {"budget_label": "medium_25_25_8", "warmup_limit": 25, "anneal_limit": 25, "max_data": 8},
)

CSV_HEADERS = (
    "budget_label",
    "block_index",
    "temperature",
    "energy_eave",
    "causal_precision",
    "causal_recall",
    "causal_f1",
    "missing_relations_count",
    "extra_relations_count",
    "total_relations_target",
    "total_relations_induced",
    "correct_relations",
    "exact_match",
    "success_flag",
    "phase",
    "accepted_or_current",
    "note",
)


def _fmt(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if math.isfinite(value):
            return f"{value:.10g}"
        return "NA"
    return str(value)


def _fmt_f(value: float) -> str:
    if math.isfinite(value):
        return f"{value:.6g}"
    return "NA"


def _safe_div(numerator: float, denominator: float) -> float:
    if denominator == 0.0:
        return 0.0
    return numerator / denominator


def _make_case() -> vs.SprinkleCase:
    matrix, points = vs.sprinkle_minkowski_diamond(
        n=N,
        seed=CASE_SEED,
        d_spacetime=D_SPACETIME,
    )
    return vs.SprinkleCase(
        d_spacetime=D_SPACETIME,
        n=N,
        seed=CASE_SEED,
        matrix=matrix,
        points=points,
    )


def _target_pairs(matrix: list[list[bool]]) -> set[tuple[int, int]]:
    return {
        (i, j)
        for i in range(len(matrix) - 1)
        for j in range(i + 1, len(matrix))
        if matrix[i][j]
    }


def _induced_pairs(matrix: list[list[bool]]) -> set[tuple[int, int]]:
    return {
        (i, j)
        for i in range(len(matrix) - 1)
        for j in range(i + 1, len(matrix))
        if matrix[i][j]
    }


def _causal_metrics(
    target: list[list[bool]],
    coords: list[tuple[float, ...]],
) -> dict[str, object]:
    induced = vs.induced_order_from_coords(coords)
    comparison = vs.compare_causal_orders(target, induced)
    correct_relations = len(_target_pairs(target) & _induced_pairs(induced))

    precision = _safe_div(correct_relations, comparison.total_relations_induced)
    recall = _safe_div(correct_relations, comparison.total_relations_target)
    f1 = _safe_div(2.0 * precision * recall, precision + recall)

    return {
        "causal_precision": precision,
        "causal_recall": recall,
        "causal_f1": f1,
        "missing_relations_count": len(comparison.missing_relations),
        "extra_relations_count": len(comparison.extra_relations),
        "total_relations_target": comparison.total_relations_target,
        "total_relations_induced": comparison.total_relations_induced,
        "correct_relations": correct_relations,
        "exact_match": comparison.exact_match,
        "success_flag": comparison.exact_match,
    }


def _run_one(case: vs.SprinkleCase, budget: dict[str, object]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    budget_label = str(budget["budget_label"])

    def _block_callback(sim: Any, block_idx: int, temp: float, eave: float) -> None:
        # cones.update() writes accepted states into rold/xold, and
        # validation_suite.verify_recovery documents rold/xold as the final
        # accepted state. Copy values here so later blocks cannot mutate rows.
        coords = [
            (float(sim.rold[i]), *[float(value) for value in sim.xold[i]])
            for i in range(sim.n)
        ]
        metrics = _causal_metrics(case.matrix, coords)
        rows.append({
            "budget_label": budget_label,
            "block_index": block_idx,
            "temperature": temp,
            "energy_eave": eave,
            **metrics,
            "phase": "anneal",
            "accepted_or_current": "accepted_rold_xold",
            "note": "captured by read-only block_callback after statistics/data append",
        })

    start = time.perf_counter()
    with tempfile.TemporaryDirectory() as tmpdir, contextlib.redirect_stdout(io.StringIO()):
        sim = cones.ConesSimulator(
            z=case.matrix,
            dim=TARGET_DIM,
            seed=OPTIMIZER_SEED,
            interactive=False,
            max_data=int(budget["max_data"]),
            plot_path=None,
            warmup_limit=int(budget["warmup_limit"]),
            anneal_limit=int(budget["anneal_limit"]),
            initial_temp=INITIAL_TEMP,
            cooling_factor=GAMMA,
            backend=BACKEND,
            block_callback=_block_callback,
        )
        sim.run(Path(tmpdir) / "annealer_output.txt")
    runtime = time.perf_counter() - start

    for row in rows:
        row["runtime_seconds"] = runtime
        row["warmup_limit"] = int(budget["warmup_limit"])
        row["anneal_limit"] = int(budget["anneal_limit"])
        row["max_data"] = int(budget["max_data"])

    return rows


def _write_csv(rows: list[dict[str, object]]) -> None:
    with CSV_PATH.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_HEADERS)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _fmt(row[key]) for key in CSV_HEADERS})


def _pair_deltas(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    grouped: dict[str, list[dict[str, object]]] = {}
    for row in rows:
        grouped.setdefault(str(row["budget_label"]), []).append(row)

    for budget_label, budget_rows in grouped.items():
        ordered = sorted(budget_rows, key=lambda row: int(row["block_index"]))
        for prev, curr in zip(ordered, ordered[1:]):
            delta_energy = float(curr["energy_eave"]) - float(prev["energy_eave"])
            delta_f1 = float(curr["causal_f1"]) - float(prev["causal_f1"])
            delta_recall = float(curr["causal_recall"]) - float(prev["causal_recall"])
            delta_missing = int(curr["missing_relations_count"]) - int(prev["missing_relations_count"])
            delta_extra = int(curr["extra_relations_count"]) - int(prev["extra_relations_count"])
            out.append({
                "budget_label": budget_label,
                "from_block": int(prev["block_index"]),
                "to_block": int(curr["block_index"]),
                "delta_energy_eave": delta_energy,
                "delta_causal_f1": delta_f1,
                "delta_recall": delta_recall,
                "delta_missing": delta_missing,
                "delta_extra": delta_extra,
                "energy_down_f1_down": delta_energy < 0.0 and delta_f1 < 0.0,
                "energy_down_recall_down": delta_energy < 0.0 and delta_recall < 0.0,
                "energy_down_missing_up": delta_energy < 0.0 and delta_missing > 0,
                "energy_down_extra_up": delta_energy < 0.0 and delta_extra > 0,
            })
    return out


def _pattern_counts(deltas: list[dict[str, object]]) -> dict[str, int]:
    keys = (
        "energy_down_f1_down",
        "energy_down_recall_down",
        "energy_down_missing_up",
        "energy_down_extra_up",
    )
    return {key: sum(1 for row in deltas if bool(row[key])) for key in keys}


def _write_svg(rows: list[dict[str, object]]) -> None:
    width = 980
    height = 620
    margin_left = 72
    margin_right = 28
    panel_h = 220
    top_energy = 58
    top_f1 = 346
    plot_w = width - margin_left - margin_right
    colors = {
        "short_10_10_4": "#1f77b4",
        "medium_25_25_8": "#d1495b",
    }
    grouped: dict[str, list[dict[str, object]]] = {}
    for row in rows:
        grouped.setdefault(str(row["budget_label"]), []).append(row)

    max_block = max(int(row["block_index"]) for row in rows)
    energies = [float(row["energy_eave"]) for row in rows]
    f1s = [float(row["causal_f1"]) for row in rows]

    def sx(block: int) -> float:
        if max_block <= 1:
            return margin_left
        return margin_left + (block - 1) * plot_w / (max_block - 1)

    def sy(value: float, values: list[float], top: float) -> float:
        lo = min(values)
        hi = max(values)
        if lo == hi:
            return top + panel_h / 2.0
        return top + panel_h - (value - lo) * panel_h / (hi - lo)

    parts = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "<rect width='100%' height='100%' fill='#f7f4ed'/>",
        f"<text x='{margin_left}' y='30' font-family='monospace' font-size='18'>trajectory_audit_n36: energy and causal F1 by block</text>",
    ]

    for top, title, ylabel in (
        (top_energy, "Energy Eave", "energy_eave"),
        (top_f1, "Causal F1", "causal_f1"),
    ):
        parts.extend([
            f"<text x='{margin_left}' y='{top - 16}' font-family='monospace' font-size='14'>{title}</text>",
            f"<line x1='{margin_left}' y1='{top + panel_h}' x2='{width - margin_right}' y2='{top + panel_h}' stroke='#333' stroke-width='1.5'/>",
            f"<line x1='{margin_left}' y1='{top}' x2='{margin_left}' y2='{top + panel_h}' stroke='#333' stroke-width='1.5'/>",
        ])
        for block in range(1, max_block + 1):
            x = sx(block)
            parts.append(f"<line x1='{x:.2f}' y1='{top + panel_h}' x2='{x:.2f}' y2='{top + panel_h + 5}' stroke='#333'/>")
            parts.append(f"<text x='{x:.2f}' y='{top + panel_h + 20}' text-anchor='middle' font-family='monospace' font-size='11'>{block}</text>")

        for budget_label, budget_rows in grouped.items():
            ordered = sorted(budget_rows, key=lambda row: int(row["block_index"]))
            values = energies if ylabel == "energy_eave" else f1s
            points = " ".join(
                f"{sx(int(row['block_index'])):.2f},{sy(float(row[ylabel]), values, top):.2f}"
                for row in ordered
            )
            color = colors.get(budget_label, "#555")
            parts.append(f"<polyline fill='none' stroke='{color}' stroke-width='2.5' points='{points}'/>")
            for row in ordered:
                x = sx(int(row["block_index"]))
                y = sy(float(row[ylabel]), values, top)
                parts.append(f"<circle cx='{x:.2f}' cy='{y:.2f}' r='4' fill='{color}'/>")

    legend_y = 592
    legend_x = margin_left
    for budget_label, color in colors.items():
        parts.append(f"<line x1='{legend_x}' y1='{legend_y}' x2='{legend_x + 26}' y2='{legend_y}' stroke='{color}' stroke-width='3'/>")
        parts.append(f"<text x='{legend_x + 34}' y='{legend_y + 4}' font-family='monospace' font-size='12'>{budget_label}</text>")
        legend_x += 245
    parts.append("</svg>")
    SVG_PATH.write_text("\n".join(parts) + "\n", encoding="utf-8")


def _write_markdown(rows: list[dict[str, object]], deltas: list[dict[str, object]]) -> None:
    generated_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
    counts = _pattern_counts(deltas)
    any_tension = any(
        counts[key] > 0
        for key in (
            "energy_down_f1_down",
            "energy_down_recall_down",
            "energy_down_missing_up",
            "energy_down_extra_up",
        )
    )

    lines = [
        "# Trajectory audit N=36",
        "",
        "Exploratory SORKIN-2 diagnostic for one known-truth case and one optimizer seed.",
        "",
        "## Configuration",
        "",
        f"- Command: `{COMMAND}`",
        f"- Generated at UTC: `{generated_at}`",
        f"- Output directory: `{OUT_DIR.relative_to(ROOT)}`",
        f"- CSV: `{CSV_PATH.relative_to(ROOT)}`",
        f"- SVG: `{SVG_PATH.relative_to(ROOT)}`",
        f"- family: `{FAMILY}`",
        f"- N: `{N}`",
        f"- d_spacetime: `{D_SPACETIME}`",
        f"- target spatial dim: `{TARGET_DIM}`",
        f"- case_seed: `{CASE_SEED}`",
        f"- optimizer_seed: `{OPTIMIZER_SEED}`",
        f"- T0: `{INITIAL_TEMP}`",
        f"- gamma: `{GAMMA}`",
        f"- h: `{H}`",
        f"- backend: `{BACKEND}`",
        "- This probe uses `ConesSimulator.block_callback` in read-only mode.",
        "- Checkpoint causal metrics use `sim.rold`/`sim.xold`, the accepted state documented by `validation_suite.verify_recovery` and reported by `ConesSimulator.writeout`.",
        "- The probe does not modify the historical annealer, selection rule, acceptance rule, energy function, temperature schedule, coordinates, or causal order.",
        "",
        "## Blocks",
        "",
        "| budget | block | temp | energy_eave | F1 | recall | precision | missing | extra | induced | correct | exact |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]

    for row in sorted(rows, key=lambda item: (str(item["budget_label"]), int(item["block_index"]))):
        lines.append(
            "| {budget} | {block} | {temp} | {energy} | {f1} | {recall} | {precision} | {missing} | {extra} | {induced} | {correct} | {exact} |".format(
                budget=row["budget_label"],
                block=row["block_index"],
                temp=_fmt_f(float(row["temperature"])),
                energy=_fmt_f(float(row["energy_eave"])),
                f1=_fmt_f(float(row["causal_f1"])),
                recall=_fmt_f(float(row["causal_recall"])),
                precision=_fmt_f(float(row["causal_precision"])),
                missing=row["missing_relations_count"],
                extra=row["extra_relations_count"],
                induced=row["total_relations_induced"],
                correct=row["correct_relations"],
                exact=_fmt(row["exact_match"]),
            )
        )

    lines.extend([
        "",
        "## Consecutive deltas",
        "",
        "| budget | from | to | delta_energy_eave | delta_causal_f1 | delta_recall | delta_missing | delta_extra | energy_down_f1_down | energy_down_recall_down | energy_down_missing_up | energy_down_extra_up |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | --- | --- |",
    ])
    for row in deltas:
        lines.append(
            "| {budget} | {from_block} | {to_block} | {de} | {df1} | {dr} | {dm} | {dx} | {p1} | {p2} | {p3} | {p4} |".format(
                budget=row["budget_label"],
                from_block=row["from_block"],
                to_block=row["to_block"],
                de=_fmt_f(float(row["delta_energy_eave"])),
                df1=_fmt_f(float(row["delta_causal_f1"])),
                dr=_fmt_f(float(row["delta_recall"])),
                dm=row["delta_missing"],
                dx=row["delta_extra"],
                p1=_fmt(row["energy_down_f1_down"]),
                p2=_fmt(row["energy_down_recall_down"]),
                p3=_fmt(row["energy_down_missing_up"]),
                p4=_fmt(row["energy_down_extra_up"]),
            )
        )

    lines.extend([
        "",
        "## Pattern counts",
        "",
        f"- energy_down_f1_down: `{counts['energy_down_f1_down']}`",
        f"- energy_down_recall_down: `{counts['energy_down_recall_down']}`",
        f"- energy_down_missing_up: `{counts['energy_down_missing_up']}`",
        f"- energy_down_extra_up: `{counts['energy_down_extra_up']}`",
        "",
        "## Conservative interpretation",
        "",
    ])
    if any_tension:
        lines.append(
            "At least one consecutive block shows within-trajectory energy/causal tension: the historical energy decreases while one audited causal metric worsens."
        )
    else:
        lines.append(
            "No consecutive block in this run shows the audited energy-down/causal-worse patterns."
        )
    lines.extend([
        "This is a diagnostic for one N=36 Minkowski case, one case seed, one optimizer seed, and two short budgets.",
        "It is not evidence of embeddability or non-embeddability, does not identify a physical gamma, does not establish a transition in N, and does not by itself justify changing the objective function.",
        "",
    ])
    MD_PATH.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    case = _make_case()
    rows: list[dict[str, object]] = []
    for budget in BUDGETS:
        rows.extend(_run_one(case, budget))
    deltas = _pair_deltas(rows)
    _write_csv(rows)
    _write_svg(rows)
    _write_markdown(rows, deltas)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
