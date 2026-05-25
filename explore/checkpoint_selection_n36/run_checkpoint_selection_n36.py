#!/usr/bin/env python3
"""Post-run checkpoint-selection diagnostic for SORKIN-2 N=36.

Reads the already generated schedule/seed trajectory CSV and compares
the final checkpoint against the best causal-F1 checkpoint seen during
each trajectory. This script does not run the annealer.
"""

from __future__ import annotations

import csv
import math
import sys
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


OUT_DIR = Path(__file__).resolve().parent
SOURCE_CSV = ROOT / "explore" / "schedule_seed_stability_n36" / "schedule_seed_stability_n36.csv"
CSV_PATH = OUT_DIR / "checkpoint_selection_n36.csv"
SUMMARY_CSV_PATH = OUT_DIR / "checkpoint_selection_n36_summary.csv"
MD_PATH = OUT_DIR / "checkpoint_selection_n36.md"
SVG_PATH = OUT_DIR / "checkpoint_selection_n36.svg"
COMMAND = "python3 explore/checkpoint_selection_n36/run_checkpoint_selection_n36.py"

CSV_HEADERS = (
    "schedule_label",
    "cooling_factor",
    "optimizer_seed",
    "final_block",
    "final_causal_f1",
    "final_recall",
    "final_missing_relations_count",
    "final_extra_relations_count",
    "final_energy_eave",
    "best_block",
    "best_causal_f1",
    "best_recall",
    "best_missing_relations_count",
    "best_extra_relations_count",
    "energy_at_best_causal_f1",
    "delta_best_minus_final_f1",
    "delta_best_minus_final_recall",
    "delta_missing_best_minus_final",
    "delta_extra_best_minus_final",
    "best_checkpoint_before_final",
    "best_checkpoint_is_final",
    "min_energy_block",
    "min_energy_eave",
    "best_f1_block_equals_min_energy_block",
)

SUMMARY_HEADERS = (
    "schedule_label",
    "cooling_factor",
    "groups",
    "avg_final_causal_f1",
    "avg_best_causal_f1",
    "avg_delta_best_minus_final_f1",
    "count_best_better_than_final",
    "count_best_better_than_final_gt_0p02",
    "count_best_before_final",
    "count_best_is_final",
    "count_best_f1_block_equals_min_energy_block",
    "avg_final_recall",
    "avg_best_recall",
    "avg_final_missing",
    "avg_best_missing",
    "avg_final_extra",
    "avg_best_extra",
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


def _read_source() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with SOURCE_CSV.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for raw in reader:
            rows.append({
                "optimizer_seed": int(raw["optimizer_seed"]),
                "schedule_label": raw["schedule_label"],
                "cooling_factor": float(raw["cooling_factor"]),
                "block_index": int(raw["block_index"]),
                "energy_eave": float(raw["energy_eave"]),
                "causal_recall": float(raw["causal_recall"]),
                "causal_f1": float(raw["causal_f1"]),
                "missing_relations_count": int(raw["missing_relations_count"]),
                "extra_relations_count": int(raw["extra_relations_count"]),
            })
    return rows


def _best_checkpoint(rows: list[dict[str, object]]) -> dict[str, object]:
    return sorted(
        rows,
        key=lambda row: (
            -float(row["causal_f1"]),
            -float(row["causal_recall"]),
            int(row["missing_relations_count"]),
            int(row["extra_relations_count"]),
            int(row["block_index"]),
        ),
    )[0]


def _group_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, int], list[dict[str, object]]] = {}
    for row in rows:
        grouped.setdefault((str(row["schedule_label"]), int(row["optimizer_seed"])), []).append(row)

    out: list[dict[str, object]] = []
    for (schedule_label, optimizer_seed), group in sorted(grouped.items()):
        ordered = sorted(group, key=lambda row: int(row["block_index"]))
        final = ordered[-1]
        best = _best_checkpoint(ordered)
        min_energy = min(ordered, key=lambda row: (float(row["energy_eave"]), int(row["block_index"])))

        final_f1 = float(final["causal_f1"])
        best_f1 = float(best["causal_f1"])
        final_recall = float(final["causal_recall"])
        best_recall = float(best["causal_recall"])
        final_missing = int(final["missing_relations_count"])
        best_missing = int(best["missing_relations_count"])
        final_extra = int(final["extra_relations_count"])
        best_extra = int(best["extra_relations_count"])
        best_block = int(best["block_index"])
        final_block = int(final["block_index"])
        min_energy_block = int(min_energy["block_index"])

        out.append({
            "schedule_label": schedule_label,
            "cooling_factor": float(final["cooling_factor"]),
            "optimizer_seed": optimizer_seed,
            "final_block": final_block,
            "final_causal_f1": final_f1,
            "final_recall": final_recall,
            "final_missing_relations_count": final_missing,
            "final_extra_relations_count": final_extra,
            "final_energy_eave": float(final["energy_eave"]),
            "best_block": best_block,
            "best_causal_f1": best_f1,
            "best_recall": best_recall,
            "best_missing_relations_count": best_missing,
            "best_extra_relations_count": best_extra,
            "energy_at_best_causal_f1": float(best["energy_eave"]),
            "delta_best_minus_final_f1": best_f1 - final_f1,
            "delta_best_minus_final_recall": best_recall - final_recall,
            "delta_missing_best_minus_final": best_missing - final_missing,
            "delta_extra_best_minus_final": best_extra - final_extra,
            "best_checkpoint_before_final": best_block < final_block,
            "best_checkpoint_is_final": best_block == final_block,
            "min_energy_block": min_energy_block,
            "min_energy_eave": float(min_energy["energy_eave"]),
            "best_f1_block_equals_min_energy_block": best_block == min_energy_block,
        })
    return sorted(out, key=lambda row: (float(row["cooling_factor"]), int(row["optimizer_seed"])))


def _mean(values: list[float]) -> float:
    return sum(values) / len(values)


def _summary_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[str, list[dict[str, object]]] = {}
    for row in rows:
        grouped.setdefault(str(row["schedule_label"]), []).append(row)

    summaries: list[dict[str, object]] = []
    for schedule_label, group in grouped.items():
        summaries.append({
            "schedule_label": schedule_label,
            "cooling_factor": float(group[0]["cooling_factor"]),
            "groups": len(group),
            "avg_final_causal_f1": _mean([float(row["final_causal_f1"]) for row in group]),
            "avg_best_causal_f1": _mean([float(row["best_causal_f1"]) for row in group]),
            "avg_delta_best_minus_final_f1": _mean([float(row["delta_best_minus_final_f1"]) for row in group]),
            "count_best_better_than_final": sum(1 for row in group if float(row["delta_best_minus_final_f1"]) > 0.0),
            "count_best_better_than_final_gt_0p02": sum(1 for row in group if float(row["delta_best_minus_final_f1"]) > 0.02),
            "count_best_before_final": sum(1 for row in group if bool(row["best_checkpoint_before_final"])),
            "count_best_is_final": sum(1 for row in group if bool(row["best_checkpoint_is_final"])),
            "count_best_f1_block_equals_min_energy_block": sum(1 for row in group if bool(row["best_f1_block_equals_min_energy_block"])),
            "avg_final_recall": _mean([float(row["final_recall"]) for row in group]),
            "avg_best_recall": _mean([float(row["best_recall"]) for row in group]),
            "avg_final_missing": _mean([float(row["final_missing_relations_count"]) for row in group]),
            "avg_best_missing": _mean([float(row["best_missing_relations_count"]) for row in group]),
            "avg_final_extra": _mean([float(row["final_extra_relations_count"]) for row in group]),
            "avg_best_extra": _mean([float(row["best_extra_relations_count"]) for row in group]),
        })
    return sorted(summaries, key=lambda row: float(row["cooling_factor"]))


def _write_csv(path: Path, headers: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=headers)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _fmt(row[key]) for key in headers})


def _write_svg(summaries: list[dict[str, object]]) -> None:
    width = 980
    height = 620
    margin_left = 76
    margin_right = 34
    top_f1 = 72
    top_delta = 365
    panel_h = 200
    plot_w = width - margin_left - margin_right
    schedules = [str(row["schedule_label"]) for row in summaries]

    def sx(index: int) -> float:
        if len(schedules) == 1:
            return margin_left + plot_w / 2.0
        return margin_left + index * plot_w / (len(schedules) - 1)

    def sy(value: float, values: list[float], top: float) -> float:
        lo = min(values)
        hi = max(values)
        if lo == hi:
            return top + panel_h / 2.0
        return top + panel_h - (value - lo) * panel_h / (hi - lo)

    final_values = [float(row["avg_final_causal_f1"]) for row in summaries]
    best_values = [float(row["avg_best_causal_f1"]) for row in summaries]
    all_f1_values = final_values + best_values
    delta_values = [float(row["avg_delta_best_minus_final_f1"]) for row in summaries]

    parts = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "<rect width='100%' height='100%' fill='#f7f4ed'/>",
        f"<text x='{margin_left}' y='32' font-family='monospace' font-size='18'>checkpoint_selection_n36: final vs best checkpoint</text>",
    ]
    for top, title in ((top_f1, "Average causal F1"), (top_delta, "Average best-minus-final F1")):
        parts.extend([
            f"<text x='{margin_left}' y='{top - 16}' font-family='monospace' font-size='14'>{title}</text>",
            f"<line x1='{margin_left}' y1='{top + panel_h}' x2='{width - margin_right}' y2='{top + panel_h}' stroke='#333' stroke-width='1.5'/>",
            f"<line x1='{margin_left}' y1='{top}' x2='{margin_left}' y2='{top + panel_h}' stroke='#333' stroke-width='1.5'/>",
        ])
        for idx, label in enumerate(schedules):
            x = sx(idx)
            parts.append(f"<line x1='{x:.2f}' y1='{top + panel_h}' x2='{x:.2f}' y2='{top + panel_h + 5}' stroke='#333'/>")
            parts.append(f"<text x='{x:.2f}' y='{top + panel_h + 20}' text-anchor='middle' font-family='monospace' font-size='11'>{label}</text>")

    final_points = " ".join(f"{sx(i):.2f},{sy(value, all_f1_values, top_f1):.2f}" for i, value in enumerate(final_values))
    best_points = " ".join(f"{sx(i):.2f},{sy(value, all_f1_values, top_f1):.2f}" for i, value in enumerate(best_values))
    delta_points = " ".join(f"{sx(i):.2f},{sy(value, delta_values, top_delta):.2f}" for i, value in enumerate(delta_values))
    parts.append(f"<polyline fill='none' stroke='#d1495b' stroke-width='2.8' points='{final_points}'/>")
    parts.append(f"<polyline fill='none' stroke='#1f77b4' stroke-width='2.8' points='{best_points}'/>")
    parts.append(f"<polyline fill='none' stroke='#2a9d8f' stroke-width='2.8' points='{delta_points}'/>")
    for i, value in enumerate(final_values):
        parts.append(f"<circle cx='{sx(i):.2f}' cy='{sy(value, all_f1_values, top_f1):.2f}' r='4' fill='#d1495b'/>")
    for i, value in enumerate(best_values):
        parts.append(f"<circle cx='{sx(i):.2f}' cy='{sy(value, all_f1_values, top_f1):.2f}' r='4' fill='#1f77b4'/>")
    for i, value in enumerate(delta_values):
        parts.append(f"<circle cx='{sx(i):.2f}' cy='{sy(value, delta_values, top_delta):.2f}' r='4' fill='#2a9d8f'/>")
    parts.extend([
        "<line x1='76' y1='598' x2='102' y2='598' stroke='#d1495b' stroke-width='3'/>",
        "<text x='110' y='602' font-family='monospace' font-size='12'>avg final F1</text>",
        "<line x1='260' y1='598' x2='286' y2='598' stroke='#1f77b4' stroke-width='3'/>",
        "<text x='294' y='602' font-family='monospace' font-size='12'>avg best checkpoint F1</text>",
        "<line x1='520' y1='598' x2='546' y2='598' stroke='#2a9d8f' stroke-width='3'/>",
        "<text x='554' y='602' font-family='monospace' font-size='12'>avg delta</text>",
        "</svg>",
    ])
    SVG_PATH.write_text("\n".join(parts) + "\n", encoding="utf-8")


def _write_markdown(rows: list[dict[str, object]], summaries: list[dict[str, object]]) -> None:
    generated_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
    groups_total = len(rows)
    total_better = sum(1 for row in rows if float(row["delta_best_minus_final_f1"]) > 0.0)
    total_better_gt = sum(1 for row in rows if float(row["delta_best_minus_final_f1"]) > 0.02)
    total_before = sum(1 for row in rows if bool(row["best_checkpoint_before_final"]))
    total_energy_match = sum(1 for row in rows if bool(row["best_f1_block_equals_min_energy_block"]))
    best_avg = max(summaries, key=lambda row: float(row["avg_best_causal_f1"]))

    lines = [
        "# Checkpoint selection N=36",
        "",
        "Post-run SORKIN-2 diagnostic over the existing schedule/optimizer-seed trajectory matrix.",
        "",
        "## Configuration",
        "",
        f"- Command: `{COMMAND}`",
        f"- Generated at UTC: `{generated_at}`",
        f"- Source CSV: `{SOURCE_CSV.relative_to(ROOT)}`",
        f"- Output CSV: `{CSV_PATH.relative_to(ROOT)}`",
        f"- Summary CSV: `{SUMMARY_CSV_PATH.relative_to(ROOT)}`",
        f"- SVG: `{SVG_PATH.relative_to(ROOT)}`",
        "- This script only reads trajectory CSV rows; it does not run the annealer.",
        "- Best checkpoint maximizes `causal_f1`, then `causal_recall`, then minimizes missing, then extra, then chooses the earliest block.",
        "",
        "## Per group",
        "",
        "| schedule | gamma | seed | final block | final F1 | best block | best F1 | delta F1 | final recall | best recall | final missing | best missing | min E block | best=minE |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in rows:
        lines.append(
            "| {schedule} | {gamma} | {seed} | {final_block} | {final_f1} | {best_block} | {best_f1} | {delta} | {final_recall} | {best_recall} | {final_missing} | {best_missing} | {min_block} | {match} |".format(
                schedule=row["schedule_label"],
                gamma=_fmt_f(float(row["cooling_factor"])),
                seed=row["optimizer_seed"],
                final_block=row["final_block"],
                final_f1=_fmt_f(float(row["final_causal_f1"])),
                best_block=row["best_block"],
                best_f1=_fmt_f(float(row["best_causal_f1"])),
                delta=_fmt_f(float(row["delta_best_minus_final_f1"])),
                final_recall=_fmt_f(float(row["final_recall"])),
                best_recall=_fmt_f(float(row["best_recall"])),
                final_missing=row["final_missing_relations_count"],
                best_missing=row["best_missing_relations_count"],
                min_block=row["min_energy_block"],
                match=_fmt(row["best_f1_block_equals_min_energy_block"]),
            )
        )

    lines.extend([
        "",
        "## Summary by schedule",
        "",
        "| schedule | gamma | groups | avg final F1 | avg best F1 | avg delta F1 | better | better >0.02 | before final | best is final | best=minE | avg final recall | avg best recall | avg final missing | avg best missing |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ])
    for row in summaries:
        lines.append(
            "| {schedule} | {gamma} | {groups} | {avg_final} | {avg_best} | {avg_delta} | {better} | {better_gt} | {before} | {is_final} | {match} | {final_recall} | {best_recall} | {final_missing} | {best_missing} |".format(
                schedule=row["schedule_label"],
                gamma=_fmt_f(float(row["cooling_factor"])),
                groups=row["groups"],
                avg_final=_fmt_f(float(row["avg_final_causal_f1"])),
                avg_best=_fmt_f(float(row["avg_best_causal_f1"])),
                avg_delta=_fmt_f(float(row["avg_delta_best_minus_final_f1"])),
                better=row["count_best_better_than_final"],
                better_gt=row["count_best_better_than_final_gt_0p02"],
                before=row["count_best_before_final"],
                is_final=row["count_best_is_final"],
                match=row["count_best_f1_block_equals_min_energy_block"],
                final_recall=_fmt_f(float(row["avg_final_recall"])),
                best_recall=_fmt_f(float(row["avg_best_recall"])),
                final_missing=_fmt_f(float(row["avg_final_missing"])),
                best_missing=_fmt_f(float(row["avg_best_missing"])),
            )
        )

    lines.extend([
        "",
        "## Readout",
        "",
        f"- Total groups: `{groups_total}`.",
        f"- Best checkpoint better than final: `{total_better}` of `{groups_total}`.",
        f"- Best checkpoint better than final by more than 0.02 F1: `{total_better_gt}` of `{groups_total}`.",
        f"- Best checkpoint before final: `{total_before}` of `{groups_total}`.",
        f"- Best avg checkpoint F1 schedule: `{best_avg['schedule_label']}` with avg best F1 `{_fmt_f(float(best_avg['avg_best_causal_f1']))}`.",
        f"- Best-F1 block equals minimum-energy block in `{total_energy_match}` of `{groups_total}` groups.",
        "Checkpoint selection appears more promising here than simply lowering gamma more slowly: the best causal checkpoint often occurs before the final block, while slower schedules did not robustly dominate best causal F1 in the seed-stability matrix.",
        "",
        "## Guardrails",
        "",
        "This is a post-run diagnostic only, using benchmark cases with known truth.",
        "It is not an embeddability claim, not a physical gamma claim, not an N-transition claim, and not proof of general annealer failure.",
        "It is also not a deployable criterion for truth-free cases yet, because selecting by causal F1 uses known-truth labels.",
        "",
    ])
    MD_PATH.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    source_rows = _read_source()
    rows = _group_rows(source_rows)
    summaries = _summary_rows(rows)
    _write_csv(CSV_PATH, CSV_HEADERS, rows)
    _write_csv(SUMMARY_CSV_PATH, SUMMARY_HEADERS, summaries)
    _write_svg(summaries)
    _write_markdown(rows, summaries)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
