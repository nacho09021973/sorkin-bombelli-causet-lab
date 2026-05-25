#!/usr/bin/env python3
"""N=36 schedule/optimizer-seed stability probe for SORKIN-2.

Exploratory only: fixed known-truth Minkowski case and medium budget,
varying optimizer_seed and the annealer's native geometric
cooling_factor while auditing each annealing block with block_callback.
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
CSV_PATH = OUT_DIR / "schedule_seed_stability_n36.csv"
SUMMARY_CSV_PATH = OUT_DIR / "schedule_seed_stability_n36_summary.csv"
MD_PATH = OUT_DIR / "schedule_seed_stability_n36.md"
SVG_PATH = OUT_DIR / "schedule_seed_stability_n36.svg"
COMMAND = "python3 explore/schedule_seed_stability_n36/run_schedule_seed_stability_n36.py"

FAMILY = "minkowski"
D_SPACETIME = 2
N = 36
CASE_SEED = 1959
OPTIMIZER_SEEDS = (1959, 1962, 1987, 2001)
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
    "optimizer_seed",
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

SUMMARY_HEADERS = (
    "optimizer_seed",
    "schedule_label",
    "cooling_factor",
    "final_energy_eave",
    "final_causal_f1",
    "best_causal_f1_seen",
    "block_of_best_f1",
    "final_recall",
    "final_missing",
    "final_extra",
    "minimum_energy_seen",
    "block_of_minimum_energy",
    "energy_down_f1_down_count",
    "energy_down_recall_down_count",
    "energy_down_missing_up_count",
    "energy_down_extra_up_count",
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


def _pairs(matrix: list[list[bool]]) -> set[tuple[int, int]]:
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
    correct_relations = len(_pairs(target) & _pairs(induced))
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
    optimizer_seed: int,
    schedule: dict[str, object],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    schedule_label = str(schedule["schedule_label"])
    cooling_factor = float(schedule["cooling_factor"])

    def _block_callback(sim: Any, block_idx: int, temp: float, eave: float) -> None:
        # rold/xold are the accepted state used by validation_suite.verify_recovery.
        # Copy values inside the callback so later annealing blocks cannot mutate rows.
        coords = [
            (float(sim.rold[i]), *[float(value) for value in sim.xold[i]])
            for i in range(sim.n)
        ]
        rows.append({
            "optimizer_seed": optimizer_seed,
            "schedule_label": schedule_label,
            "cooling_factor": cooling_factor,
            "budget_label": str(BUDGET["budget_label"]),
            "block_index": block_idx,
            "temperature": temp,
            "energy_eave": eave,
            **_causal_metrics(case.matrix, coords),
        })

    start = time.perf_counter()
    with tempfile.TemporaryDirectory() as tmpdir, contextlib.redirect_stdout(io.StringIO()):
        sim = cones.ConesSimulator(
            z=case.matrix,
            dim=TARGET_DIM,
            seed=optimizer_seed,
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
    grouped: dict[tuple[int, str], list[dict[str, object]]] = {}
    for row in rows:
        grouped.setdefault((int(row["optimizer_seed"]), str(row["schedule_label"])), []).append(row)

    out: list[dict[str, object]] = []
    for _key, group_rows in grouped.items():
        prev: dict[str, object] | None = None
        for row in sorted(group_rows, key=lambda item: int(item["block_index"])):
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
    return sorted(out, key=lambda item: (int(item["optimizer_seed"]), str(item["schedule_label"]), int(item["block_index"])))


def _summary_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[int, str], list[dict[str, object]]] = {}
    for row in rows:
        grouped.setdefault((int(row["optimizer_seed"]), str(row["schedule_label"])), []).append(row)

    summaries: list[dict[str, object]] = []
    for (optimizer_seed, schedule_label), group_rows in grouped.items():
        ordered = sorted(group_rows, key=lambda item: int(item["block_index"]))
        final = ordered[-1]
        best_f1 = max(ordered, key=lambda item: float(item["causal_f1"]))
        min_energy = min(ordered, key=lambda item: float(item["energy_eave"]))
        summaries.append({
            "optimizer_seed": optimizer_seed,
            "schedule_label": schedule_label,
            "cooling_factor": float(final["cooling_factor"]),
            "final_energy_eave": float(final["energy_eave"]),
            "final_causal_f1": float(final["causal_f1"]),
            "best_causal_f1_seen": float(best_f1["causal_f1"]),
            "block_of_best_f1": int(best_f1["block_index"]),
            "final_recall": float(final["causal_recall"]),
            "final_missing": int(final["missing_relations_count"]),
            "final_extra": int(final["extra_relations_count"]),
            "minimum_energy_seen": float(min_energy["energy_eave"]),
            "block_of_minimum_energy": int(min_energy["block_index"]),
            "energy_down_f1_down_count": sum(1 for row in ordered if bool(row["energy_down_f1_down"])),
            "energy_down_recall_down_count": sum(1 for row in ordered if bool(row["energy_down_recall_down"])),
            "energy_down_missing_up_count": sum(1 for row in ordered if bool(row["energy_down_missing_up"])),
            "energy_down_extra_up_count": sum(1 for row in ordered if bool(row["energy_down_extra_up"])),
        })
    return sorted(summaries, key=lambda item: (int(item["optimizer_seed"]), float(item["cooling_factor"])))


def _write_csv(path: Path, headers: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=headers)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _fmt(row[key]) for key in headers})


def _winner_counts(summaries: list[dict[str, object]], metric: str, *, lower_is_better: bool = False) -> dict[str, int]:
    counts = {str(schedule["schedule_label"]): 0 for schedule in SCHEDULES}
    for seed in OPTIMIZER_SEEDS:
        seed_rows = [row for row in summaries if int(row["optimizer_seed"]) == seed]
        winner = min(seed_rows, key=lambda row: float(row[metric])) if lower_is_better else max(seed_rows, key=lambda row: float(row[metric]))
        counts[str(winner["schedule_label"])] += 1
    return counts


def _write_svg(summaries: list[dict[str, object]]) -> None:
    width = 1120
    height = 700
    margin_left = 78
    margin_right = 34
    top = 72
    panel_h = 160
    panel_gap = 70
    plot_w = width - margin_left - margin_right
    colors = {
        "gamma_0p5": "#1f77b4",
        "gamma_0p8": "#d1495b",
        "gamma_0p9": "#2a9d8f",
        "gamma_0p95": "#8a5a44",
    }
    x_positions = {seed: margin_left + idx * plot_w / (len(OPTIMIZER_SEEDS) - 1) for idx, seed in enumerate(OPTIMIZER_SEEDS)}

    def y_scale(value: float, values: list[float], panel_top: float, invert: bool = False) -> float:
        lo = min(values)
        hi = max(values)
        if lo == hi:
            return panel_top + panel_h / 2.0
        frac = (value - lo) / (hi - lo)
        if invert:
            frac = 1.0 - frac
        return panel_top + panel_h - frac * panel_h

    panels = (
        ("Final causal F1", "final_causal_f1", False),
        ("Best causal F1 seen", "best_causal_f1_seen", False),
        ("energy_down_f1_down count", "energy_down_f1_down_count", True),
    )
    parts = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "<rect width='100%' height='100%' fill='#f7f4ed'/>",
        f"<text x='{margin_left}' y='32' font-family='monospace' font-size='18'>schedule_seed_stability_n36: schedule x optimizer_seed</text>",
    ]
    for panel_idx, (title, metric, invert) in enumerate(panels):
        panel_top = top + panel_idx * (panel_h + panel_gap)
        values = [float(row[metric]) for row in summaries]
        parts.extend([
            f"<text x='{margin_left}' y='{panel_top - 16}' font-family='monospace' font-size='14'>{title}</text>",
            f"<line x1='{margin_left}' y1='{panel_top + panel_h}' x2='{width - margin_right}' y2='{panel_top + panel_h}' stroke='#333' stroke-width='1.5'/>",
            f"<line x1='{margin_left}' y1='{panel_top}' x2='{margin_left}' y2='{panel_top + panel_h}' stroke='#333' stroke-width='1.5'/>",
        ])
        for seed, x in x_positions.items():
            parts.append(f"<line x1='{x:.2f}' y1='{panel_top + panel_h}' x2='{x:.2f}' y2='{panel_top + panel_h + 5}' stroke='#333'/>")
            parts.append(f"<text x='{x:.2f}' y='{panel_top + panel_h + 20}' text-anchor='middle' font-family='monospace' font-size='11'>{seed}</text>")
        for schedule in SCHEDULES:
            label = str(schedule["schedule_label"])
            schedule_rows = [row for row in summaries if str(row["schedule_label"]) == label]
            points = " ".join(
                f"{x_positions[int(row['optimizer_seed'])]:.2f},{y_scale(float(row[metric]), values, panel_top, invert):.2f}"
                for row in schedule_rows
            )
            color = colors[label]
            parts.append(f"<polyline fill='none' stroke='{color}' stroke-width='2.4' points='{points}'/>")
            for row in schedule_rows:
                x = x_positions[int(row["optimizer_seed"])]
                y = y_scale(float(row[metric]), values, panel_top, invert)
                parts.append(f"<circle cx='{x:.2f}' cy='{y:.2f}' r='3.8' fill='{color}'/>")

    legend_y = 676
    legend_x = margin_left
    for label, color in colors.items():
        parts.append(f"<line x1='{legend_x}' y1='{legend_y}' x2='{legend_x + 26}' y2='{legend_y}' stroke='{color}' stroke-width='3'/>")
        parts.append(f"<text x='{legend_x + 34}' y='{legend_y + 4}' font-family='monospace' font-size='12'>{label}</text>")
        legend_x += 220
    parts.append("</svg>")
    SVG_PATH.write_text("\n".join(parts) + "\n", encoding="utf-8")


def _per_seed_winners(summaries: list[dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for seed in OPTIMIZER_SEEDS:
        seed_rows = [row for row in summaries if int(row["optimizer_seed"]) == seed]
        final = max(seed_rows, key=lambda row: float(row["final_causal_f1"]))
        best = max(seed_rows, key=lambda row: float(row["best_causal_f1_seen"]))
        tension = min(seed_rows, key=lambda row: int(row["energy_down_f1_down_count"]))
        energy = min(seed_rows, key=lambda row: float(row["minimum_energy_seen"]))
        rows.append({
            "optimizer_seed": seed,
            "winner_final_f1": final["schedule_label"],
            "winner_best_f1": best["schedule_label"],
            "winner_lowest_tension": tension["schedule_label"],
            "winner_min_energy": energy["schedule_label"],
            "min_energy_matches_final_f1": energy["schedule_label"] == final["schedule_label"],
            "min_energy_matches_best_f1": energy["schedule_label"] == best["schedule_label"],
        })
    return rows


def _write_markdown(rows: list[dict[str, object]], summaries: list[dict[str, object]]) -> None:
    generated_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
    final_counts = _winner_counts(summaries, "final_causal_f1")
    best_counts = _winner_counts(summaries, "best_causal_f1_seen")
    tension_counts = _winner_counts(summaries, "energy_down_f1_down_count", lower_is_better=True)
    per_seed = _per_seed_winners(summaries)
    gamma_0p5_tension = {
        int(row["optimizer_seed"]): int(row["energy_down_f1_down_count"])
        for row in summaries
        if row["schedule_label"] == "gamma_0p5"
    }
    slower_reductions = 0
    for seed in OPTIMIZER_SEEDS:
        seed_rows = [
            row for row in summaries
            if int(row["optimizer_seed"]) == seed and row["schedule_label"] != "gamma_0p5"
        ]
        if min(int(row["energy_down_f1_down_count"]) for row in seed_rows) < gamma_0p5_tension[seed]:
            slower_reductions += 1

    lines = [
        "# Schedule seed stability N=36",
        "",
        "Exploratory SORKIN-2 diagnostic testing whether the schedule signal persists across optimizer seeds for one fixed known-truth causal set.",
        "",
        "## Configuration",
        "",
        f"- Command: `{COMMAND}`",
        f"- Generated at UTC: `{generated_at}`",
        f"- Output directory: `{OUT_DIR.relative_to(ROOT)}`",
        f"- CSV: `{CSV_PATH.relative_to(ROOT)}`",
        f"- Summary CSV: `{SUMMARY_CSV_PATH.relative_to(ROOT)}`",
        f"- SVG: `{SVG_PATH.relative_to(ROOT)}`",
        f"- family: `{FAMILY}`",
        f"- N: `{N}`",
        f"- d_spacetime: `{D_SPACETIME}`",
        f"- target spatial dim: `{TARGET_DIM}`",
        f"- case_seed: `{CASE_SEED}`",
        f"- optimizer_seeds: `{', '.join(str(seed) for seed in OPTIMIZER_SEEDS)}`",
        f"- T0: `{INITIAL_TEMP}`",
        f"- budget: `{BUDGET['budget_label']}`",
        f"- warmup_limit: `{BUDGET['warmup_limit']}`",
        f"- anneal_limit: `{BUDGET['anneal_limit']}`",
        f"- max_data: `{BUDGET['max_data']}`",
        "- Schedules vary only native `cooling_factor`: `0.5`, `0.8`, `0.9`, `0.95`.",
        "- block_callback reads `sim.rold`/`sim.xold` only; historical energy, move set, acceptance rule, and schedule mechanism are unchanged.",
        "",
        "## Summary matrix",
        "",
        "| seed | schedule | gamma | final F1 | best F1 seen | block best F1 | final recall | final missing | min Eave | block min E | Edown F1down | Edown recalldown | Edown missingup |",
        "| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in summaries:
        lines.append(
            "| {seed} | {schedule} | {gamma} | {final_f1} | {best_f1} | {block_best} | {recall} | {missing} | {min_energy} | {block_min} | {c1} | {c2} | {c3} |".format(
                seed=row["optimizer_seed"],
                schedule=row["schedule_label"],
                gamma=_fmt_f(float(row["cooling_factor"])),
                final_f1=_fmt_f(float(row["final_causal_f1"])),
                best_f1=_fmt_f(float(row["best_causal_f1_seen"])),
                block_best=row["block_of_best_f1"],
                recall=_fmt_f(float(row["final_recall"])),
                missing=row["final_missing"],
                min_energy=_fmt_f(float(row["minimum_energy_seen"])),
                block_min=row["block_of_minimum_energy"],
                c1=row["energy_down_f1_down_count"],
                c2=row["energy_down_recall_down_count"],
                c3=row["energy_down_missing_up_count"],
            )
        )

    lines.extend([
        "",
        "## Per-seed winners",
        "",
        "| seed | final F1 winner | best F1 winner | lowest tension winner | min energy winner | min energy = final F1? | min energy = best F1? |",
        "| ---: | --- | --- | --- | --- | --- | --- |",
    ])
    for row in per_seed:
        lines.append(
            "| {seed} | {final} | {best} | {tension} | {energy} | {ef} | {eb} |".format(
                seed=row["optimizer_seed"],
                final=row["winner_final_f1"],
                best=row["winner_best_f1"],
                tension=row["winner_lowest_tension"],
                energy=row["winner_min_energy"],
                ef=_fmt(row["min_energy_matches_final_f1"]),
                eb=_fmt(row["min_energy_matches_best_f1"]),
            )
        )

    lines.extend([
        "",
        "## Aggregate winner counts",
        "",
        "| schedule | final F1 wins | best F1 wins | lowest tension wins |",
        "| --- | ---: | ---: | ---: |",
    ])
    for schedule in SCHEDULES:
        label = str(schedule["schedule_label"])
        lines.append(f"| {label} | {final_counts[label]} | {best_counts[label]} | {tension_counts[label]} |")

    total_min_energy_matches_final = sum(1 for row in per_seed if bool(row["min_energy_matches_final_f1"]))
    total_min_energy_matches_best = sum(1 for row in per_seed if bool(row["min_energy_matches_best_f1"]))
    lines.extend([
        "",
        "## Readout",
        "",
        f"- Minimum energy coincides with the final-F1 winner in `{total_min_energy_matches_final}` of `{len(OPTIMIZER_SEEDS)}` optimizer seeds.",
        f"- Minimum energy coincides with the best-trajectory-F1 winner in `{total_min_energy_matches_best}` of `{len(OPTIMIZER_SEEDS)}` optimizer seeds.",
        f"- At least one slower gamma reduces `energy_down_f1_down` relative to gamma 0.5 in `{slower_reductions}` of `{len(OPTIMIZER_SEEDS)}` optimizer seeds.",
    ])
    if max(final_counts.values()) == len(OPTIMIZER_SEEDS) or max(best_counts.values()) == len(OPTIMIZER_SEEDS):
        lines.append(
            "One schedule dominates one causal-F1 readout across all optimizer seeds in this small matrix; that would support cooling as a robust algorithmic factor for this fixed causal set."
        )
    elif max(final_counts.values()) >= 2 or max(best_counts.values()) >= 2:
        lines.append(
            "The winners are not fully seed-invariant, but at least one schedule wins multiple seeds; this is mixed exploratory evidence for schedule sensitivity plus optimizer-seed dependence."
        )
    else:
        lines.append(
            "Winners are dispersed across optimizer seeds; this favors a basin-accessibility interpretation over a single robust geometric gamma."
        )
    if slower_reductions > 0:
        lines.append(
            "Slower cooling reduces at least one tension count in part of the matrix, so thermal mobility remains a plausible diagnostic lever."
        )
    else:
        lines.append(
            "Slower cooling does not reduce the audited tension count in this matrix; that would point more strongly toward basin selection or the surrogate energy."
        )

    lines.extend([
        "",
        "## Guardrails",
        "",
        "This is exploratory: one N, one family, one case seed, four optimizer seeds, one budget, and the unchanged historical objective function.",
        "It is not evidence of embeddability or non-embeddability, does not identify a physical gamma, does not establish an N transition, and does not by itself justify changing the objective function.",
        "",
    ])
    MD_PATH.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    case = _make_case()
    rows: list[dict[str, object]] = []
    for optimizer_seed in OPTIMIZER_SEEDS:
        for schedule in SCHEDULES:
            rows.extend(_run_one(case, optimizer_seed, schedule))
    rows_with_deltas = _with_deltas(rows)
    summaries = _summary_rows(rows_with_deltas)
    _write_csv(CSV_PATH, CSV_HEADERS, rows_with_deltas)
    _write_csv(SUMMARY_CSV_PATH, SUMMARY_HEADERS, summaries)
    _write_svg(summaries)
    _write_markdown(rows_with_deltas, summaries)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
