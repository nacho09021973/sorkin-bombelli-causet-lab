#!/usr/bin/env python3
"""N=36 geometric schedule probe for SORKIN-2.

Exploratory only: fixed known-truth Minkowski case and optimizer seed,
varying only the annealer's native cooling_factor while auditing each
annealing block through the read-only block_callback hook.
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
CSV_PATH = OUT_DIR / "schedule_probe_n36.csv"
MD_PATH = OUT_DIR / "schedule_probe_n36.md"
SVG_PATH = OUT_DIR / "schedule_probe_n36.svg"
COMMAND = "python3 explore/schedule_probe_n36/run_schedule_probe_n36.py"

FAMILY = "minkowski"
D_SPACETIME = 2
N = 36
CASE_SEED = 1959
OPTIMIZER_SEED = 1987
TARGET_DIM = D_SPACETIME - 1
INITIAL_TEMP = 100.0
BACKEND = "cpu"

BUDGET = {
    "budget_label": "medium_25_25_8",
    "warmup_limit": 25,
    "anneal_limit": 25,
    "max_data": 8,
}

SCHEDULES = (
    {"schedule_label": "gamma_0p5", "cooling_factor": 0.50},
    {"schedule_label": "gamma_0p8", "cooling_factor": 0.80},
    {"schedule_label": "gamma_0p9", "cooling_factor": 0.90},
    {"schedule_label": "gamma_0p95", "cooling_factor": 0.95},
)

CSV_HEADERS = (
    "schedule_label",
    "cooling_factor",
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
    "delta_energy_eave",
    "delta_causal_f1",
    "delta_causal_recall",
    "delta_missing_relations",
    "delta_extra_relations",
    "energy_down_f1_down",
    "energy_down_recall_down",
    "energy_down_missing_up",
    "energy_down_extra_up",
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


def _run_one(
    case: vs.SprinkleCase,
    schedule: dict[str, object],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    schedule_label = str(schedule["schedule_label"])
    cooling_factor = float(schedule["cooling_factor"])

    def _block_callback(sim: Any, block_idx: int, temp: float, eave: float) -> None:
        # cones.update() stores accepted states in rold/xold; validation_suite
        # verify_recovery uses the same state. Copy here to avoid later mutation.
        coords = [
            (float(sim.rold[i]), *[float(value) for value in sim.xold[i]])
            for i in range(sim.n)
        ]
        metrics = _causal_metrics(case.matrix, coords)
        rows.append({
            "schedule_label": schedule_label,
            "cooling_factor": cooling_factor,
            "budget_label": str(BUDGET["budget_label"]),
            "block_index": block_idx,
            "temperature": temp,
            "energy_eave": eave,
            **metrics,
        })

    start = time.perf_counter()
    with tempfile.TemporaryDirectory() as tmpdir, contextlib.redirect_stdout(io.StringIO()):
        sim = cones.ConesSimulator(
            z=case.matrix,
            dim=TARGET_DIM,
            seed=OPTIMIZER_SEED,
            interactive=False,
            max_data=int(BUDGET["max_data"]),
            plot_path=None,
            warmup_limit=int(BUDGET["warmup_limit"]),
            anneal_limit=int(BUDGET["anneal_limit"]),
            initial_temp=INITIAL_TEMP,
            cooling_factor=cooling_factor,
            backend=BACKEND,
            block_callback=_block_callback,
        )
        sim.run(Path(tmpdir) / "annealer_output.txt")
    runtime = time.perf_counter() - start

    for row in rows:
        row["runtime_seconds"] = runtime
    return rows


def _with_deltas(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    grouped: dict[str, list[dict[str, object]]] = {}
    for row in rows:
        grouped.setdefault(str(row["schedule_label"]), []).append(row)

    for _schedule_label, schedule_rows in grouped.items():
        prev: dict[str, object] | None = None
        for row in sorted(schedule_rows, key=lambda item: int(item["block_index"])):
            item = dict(row)
            if prev is None:
                item.update({
                    "delta_energy_eave": float("nan"),
                    "delta_causal_f1": float("nan"),
                    "delta_causal_recall": float("nan"),
                    "delta_missing_relations": float("nan"),
                    "delta_extra_relations": float("nan"),
                    "energy_down_f1_down": False,
                    "energy_down_recall_down": False,
                    "energy_down_missing_up": False,
                    "energy_down_extra_up": False,
                })
            else:
                delta_energy = float(row["energy_eave"]) - float(prev["energy_eave"])
                delta_f1 = float(row["causal_f1"]) - float(prev["causal_f1"])
                delta_recall = float(row["causal_recall"]) - float(prev["causal_recall"])
                delta_missing = int(row["missing_relations_count"]) - int(prev["missing_relations_count"])
                delta_extra = int(row["extra_relations_count"]) - int(prev["extra_relations_count"])
                item.update({
                    "delta_energy_eave": delta_energy,
                    "delta_causal_f1": delta_f1,
                    "delta_causal_recall": delta_recall,
                    "delta_missing_relations": delta_missing,
                    "delta_extra_relations": delta_extra,
                    "energy_down_f1_down": delta_energy < 0.0 and delta_f1 < 0.0,
                    "energy_down_recall_down": delta_energy < 0.0 and delta_recall < 0.0,
                    "energy_down_missing_up": delta_energy < 0.0 and delta_missing > 0,
                    "energy_down_extra_up": delta_energy < 0.0 and delta_extra > 0,
                })
            out.append(item)
            prev = row

    return sorted(out, key=lambda item: (str(item["schedule_label"]), int(item["block_index"])))


def _write_csv(rows: list[dict[str, object]]) -> None:
    with CSV_PATH.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_HEADERS)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _fmt(row[key]) for key in CSV_HEADERS})


def _summaries(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    summaries: list[dict[str, object]] = []
    grouped: dict[str, list[dict[str, object]]] = {}
    for row in rows:
        grouped.setdefault(str(row["schedule_label"]), []).append(row)

    for schedule_label, schedule_rows in grouped.items():
        ordered = sorted(schedule_rows, key=lambda item: int(item["block_index"]))
        final = ordered[-1]
        best_f1 = max(ordered, key=lambda item: float(item["causal_f1"]))
        min_energy = min(ordered, key=lambda item: float(item["energy_eave"]))
        summaries.append({
            "schedule_label": schedule_label,
            "cooling_factor": float(final["cooling_factor"]),
            "final_energy_eave": float(final["energy_eave"]),
            "final_causal_f1": float(final["causal_f1"]),
            "final_recall": float(final["causal_recall"]),
            "final_missing": int(final["missing_relations_count"]),
            "final_extra": int(final["extra_relations_count"]),
            "best_causal_f1_seen": float(best_f1["causal_f1"]),
            "block_of_best_causal_f1": int(best_f1["block_index"]),
            "minimum_energy_seen": float(min_energy["energy_eave"]),
            "block_of_minimum_energy": int(min_energy["block_index"]),
            "count_energy_down_f1_down": sum(1 for row in ordered if bool(row["energy_down_f1_down"])),
            "count_energy_down_recall_down": sum(1 for row in ordered if bool(row["energy_down_recall_down"])),
            "count_energy_down_missing_up": sum(1 for row in ordered if bool(row["energy_down_missing_up"])),
            "count_energy_down_extra_up": sum(1 for row in ordered if bool(row["energy_down_extra_up"])),
        })
    return sorted(summaries, key=lambda item: float(item["cooling_factor"]))


def _write_svg(rows: list[dict[str, object]]) -> None:
    width = 1080
    height = 660
    margin_left = 76
    margin_right = 34
    panel_h = 230
    top_energy = 62
    top_f1 = 368
    plot_w = width - margin_left - margin_right
    colors = {
        "gamma_0p5": "#1f77b4",
        "gamma_0p8": "#d1495b",
        "gamma_0p9": "#2a9d8f",
        "gamma_0p95": "#8a5a44",
    }
    grouped: dict[str, list[dict[str, object]]] = {}
    for row in rows:
        grouped.setdefault(str(row["schedule_label"]), []).append(row)

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
        f"<text x='{margin_left}' y='32' font-family='monospace' font-size='18'>schedule_probe_n36: geometric cooling comparison</text>",
    ]

    for top, title, metric, values in (
        (top_energy, "Energy Eave", "energy_eave", energies),
        (top_f1, "Causal F1", "causal_f1", f1s),
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

        for schedule_label, schedule_rows in grouped.items():
            ordered = sorted(schedule_rows, key=lambda item: int(item["block_index"]))
            points = " ".join(
                f"{sx(int(row['block_index'])):.2f},{sy(float(row[metric]), values, top):.2f}"
                for row in ordered
            )
            color = colors.get(schedule_label, "#555")
            parts.append(f"<polyline fill='none' stroke='{color}' stroke-width='2.4' points='{points}'/>")
            for row in ordered:
                x = sx(int(row["block_index"]))
                y = sy(float(row[metric]), values, top)
                parts.append(f"<circle cx='{x:.2f}' cy='{y:.2f}' r='3.8' fill='{color}'/>")

    legend_y = 632
    legend_x = margin_left
    for schedule_label, color in colors.items():
        parts.append(f"<line x1='{legend_x}' y1='{legend_y}' x2='{legend_x + 26}' y2='{legend_y}' stroke='{color}' stroke-width='3'/>")
        parts.append(f"<text x='{legend_x + 34}' y='{legend_y + 4}' font-family='monospace' font-size='12'>{schedule_label}</text>")
        legend_x += 220
    parts.append("</svg>")
    SVG_PATH.write_text("\n".join(parts) + "\n", encoding="utf-8")


def _write_markdown(rows: list[dict[str, object]]) -> None:
    generated_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
    summaries = _summaries(rows)
    best_final = max(summaries, key=lambda item: float(item["final_causal_f1"]))
    best_seen = max(summaries, key=lambda item: float(item["best_causal_f1_seen"]))
    min_energy = min(summaries, key=lambda item: float(item["minimum_energy_seen"]))
    min_tension_count = min(int(item["count_energy_down_f1_down"]) for item in summaries)
    tension_best_labels = [
        str(item["schedule_label"])
        for item in summaries
        if int(item["count_energy_down_f1_down"]) == min_tension_count
    ]

    energy_matches_best_final = min_energy["schedule_label"] == best_final["schedule_label"]
    energy_matches_best_seen = min_energy["schedule_label"] == best_seen["schedule_label"]

    lines = [
        "# Schedule probe N=36",
        "",
        "Exploratory SORKIN-2 diagnostic comparing native geometric cooling schedules for one known-truth case and one optimizer seed.",
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
        f"- budget: `{BUDGET['budget_label']}`",
        f"- warmup_limit: `{BUDGET['warmup_limit']}`",
        f"- anneal_limit: `{BUDGET['anneal_limit']}`",
        f"- max_data: `{BUDGET['max_data']}`",
        f"- backend: `{BACKEND}`",
        "- Schedules vary only `cooling_factor`: `0.5`, `0.8`, `0.9`, `0.95`.",
        "- This probe uses `ConesSimulator.block_callback` in read-only mode.",
        "- Checkpoint causal metrics use `sim.rold`/`sim.xold`, the accepted state documented by `validation_suite.verify_recovery` and reported by `ConesSimulator.writeout`.",
        "- The historical Bombelli energy, acceptance rule, move set, and internal annealer dynamics are unchanged.",
        "",
        "## Summary by schedule",
        "",
        "| schedule | gamma | final_energy_eave | final F1 | final recall | final missing | final extra | best F1 seen | block best F1 | min energy seen | block min energy | Edown F1down | Edown recalldown | Edown missingup | Edown extraup |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]

    for row in summaries:
        lines.append(
            "| {schedule} | {gamma} | {final_energy} | {final_f1} | {final_recall} | {final_missing} | {final_extra} | {best_f1} | {block_best} | {min_energy} | {block_min} | {c1} | {c2} | {c3} | {c4} |".format(
                schedule=row["schedule_label"],
                gamma=_fmt_f(float(row["cooling_factor"])),
                final_energy=_fmt_f(float(row["final_energy_eave"])),
                final_f1=_fmt_f(float(row["final_causal_f1"])),
                final_recall=_fmt_f(float(row["final_recall"])),
                final_missing=row["final_missing"],
                final_extra=row["final_extra"],
                best_f1=_fmt_f(float(row["best_causal_f1_seen"])),
                block_best=row["block_of_best_causal_f1"],
                min_energy=_fmt_f(float(row["minimum_energy_seen"])),
                block_min=row["block_of_minimum_energy"],
                c1=row["count_energy_down_f1_down"],
                c2=row["count_energy_down_recall_down"],
                c3=row["count_energy_down_missing_up"],
                c4=row["count_energy_down_extra_up"],
            )
        )

    lines.extend([
        "",
        "## Readout",
        "",
        f"- Best final causal F1: `{best_final['schedule_label']}` with F1 `{_fmt_f(float(best_final['final_causal_f1']))}`.",
        f"- Best causal F1 seen anywhere in the trajectory: `{best_seen['schedule_label']}` with F1 `{_fmt_f(float(best_seen['best_causal_f1_seen']))}` at block `{best_seen['block_of_best_causal_f1']}`.",
        f"- Minimum energy seen: `{min_energy['schedule_label']}` with Eave `{_fmt_f(float(min_energy['minimum_energy_seen']))}` at block `{min_energy['block_of_minimum_energy']}`.",
    ])
    if energy_matches_best_final and energy_matches_best_seen:
        lines.append("- In this run, the schedule with the lowest observed energy also has the best final and best seen causal F1.")
    elif energy_matches_best_final:
        lines.append("- In this run, the schedule with the lowest observed energy matches the best final F1 schedule, but not the best F1 seen along the trajectory.")
    elif energy_matches_best_seen:
        lines.append("- In this run, the schedule with the lowest observed energy matches the best F1 seen along the trajectory, but not the best final F1 schedule.")
    else:
        lines.append("- In this run, the lowest observed energy does not coincide with the best final F1 schedule or the best trajectory F1 schedule.")

    lines.append(
        "- Fewest `energy_down_f1_down` steps: `{labels}` with `{count}` such steps.".format(
            labels="`, `".join(tension_best_labels),
            count=min_tension_count,
        )
    )
    if float(best_final["cooling_factor"]) > 0.5 or float(best_seen["cooling_factor"]) > 0.5:
        lines.append(
            "Higher gamma improves at least one audited causal-F1 readout in this run, which is exploratory evidence that cooling rate and thermal mobility matter for this case."
        )
    else:
        lines.append(
            "Higher gamma does not improve the causal-F1 readouts in this run; this is compatible with the problem sitting deeper in the energy surrogate or basin selection for this case."
        )

    if min_tension_count < int(
        next(row for row in summaries if row["schedule_label"] == "gamma_0p5")["count_energy_down_f1_down"]
    ):
        lines.append(
            "A slower schedule reduces `energy_down_f1_down` counts relative to gamma 0.5 in this run."
        )
    else:
        lines.append(
            "The slower schedules do not reduce `energy_down_f1_down` counts relative to gamma 0.5 in this run."
        )

    lines.extend([
        "",
        "## Guardrails",
        "",
        "This is exploratory: one N, one family, one case seed, one optimizer seed, one budget, and the unchanged historical objective function.",
        "It is not evidence of embeddability or non-embeddability, does not identify a physical gamma, does not establish an N transition, and does not by itself justify changing the objective function.",
        "",
    ])
    MD_PATH.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    case = _make_case()
    rows: list[dict[str, object]] = []
    for schedule in SCHEDULES:
        rows.extend(_run_one(case, schedule))
    rows_with_deltas = _with_deltas(rows)
    _write_csv(rows_with_deltas)
    _write_svg(rows_with_deltas)
    _write_markdown(rows_with_deltas)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
