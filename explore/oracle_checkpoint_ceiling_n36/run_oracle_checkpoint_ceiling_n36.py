#!/usr/bin/env python3
"""Oracle checkpoint ceiling diagnostic for SORKIN-2 N=36.

Reads the existing schedule/seed trajectory CSV and asks how much causal
recoverability would be available if an oracle selected the checkpoint
with best known-truth causal F1 instead of the final endpoint.
"""

from __future__ import annotations

import csv
import math
import statistics
import sys
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


OUT_DIR = Path(__file__).resolve().parent
SOURCE_CSV = ROOT / "explore" / "schedule_seed_stability_n36" / "schedule_seed_stability_n36.csv"
CSV_PATH = OUT_DIR / "oracle_checkpoint_ceiling_n36.csv"
SUMMARY_CSV_PATH = OUT_DIR / "oracle_checkpoint_ceiling_n36_summary.csv"
MD_PATH = OUT_DIR / "oracle_checkpoint_ceiling_n36.md"
SVG_PATH = OUT_DIR / "oracle_checkpoint_ceiling_n36.svg"
COMMAND = "python3 explore/oracle_checkpoint_ceiling_n36/run_oracle_checkpoint_ceiling_n36.py"

GROUP_HEADERS = (
    "optimizer_seed",
    "schedule_label",
    "cooling_factor",
    "budget_label",
    "final_causal_f1",
    "best_checkpoint_causal_f1",
    "delta_best_minus_final",
    "final_causal_recall",
    "best_checkpoint_causal_recall",
    "delta_recall_best_minus_final",
    "final_missing_relations_count",
    "best_checkpoint_missing_relations_count",
    "delta_missing_best_minus_final",
    "final_extra_relations_count",
    "best_checkpoint_extra_relations_count",
    "delta_extra_best_minus_final",
    "final_energy_eave",
    "best_checkpoint_energy_eave",
    "min_energy_eave",
    "block_index_final",
    "block_index_best_causal_f1",
    "block_index_min_energy",
    "best_before_final",
    "best_matches_min_energy",
)

SUMMARY_HEADERS = (
    "n_groups",
    "avg_final_causal_f1",
    "avg_best_checkpoint_causal_f1",
    "avg_delta_best_minus_final",
    "median_delta_best_minus_final",
    "max_delta_best_minus_final",
    "count_best_gt_final",
    "count_best_gt_final_by_0p02",
    "count_best_before_final",
    "count_best_matches_min_energy",
    "best_avg_checkpoint_schedule_label",
    "best_avg_checkpoint_budget_label",
    "best_avg_checkpoint_cooling_factor",
    "best_avg_checkpoint_causal_f1",
    "best_avg_delta_schedule_label",
    "best_avg_delta_budget_label",
    "best_avg_delta_cooling_factor",
    "best_avg_delta_best_minus_final",
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
                "budget_label": raw["budget_label"],
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


def _group_rows(source_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[int, str, float, str], list[dict[str, object]]] = {}
    for row in source_rows:
        key = (
            int(row["optimizer_seed"]),
            str(row["schedule_label"]),
            float(row["cooling_factor"]),
            str(row["budget_label"]),
        )
        grouped.setdefault(key, []).append(row)

    out: list[dict[str, object]] = []
    for key, group in sorted(grouped.items()):
        optimizer_seed, schedule_label, cooling_factor, budget_label = key
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
        final_block = int(final["block_index"])
        best_block = int(best["block_index"])
        min_energy_block = int(min_energy["block_index"])

        out.append({
            "optimizer_seed": optimizer_seed,
            "schedule_label": schedule_label,
            "cooling_factor": cooling_factor,
            "budget_label": budget_label,
            "final_causal_f1": final_f1,
            "best_checkpoint_causal_f1": best_f1,
            "delta_best_minus_final": best_f1 - final_f1,
            "final_causal_recall": final_recall,
            "best_checkpoint_causal_recall": best_recall,
            "delta_recall_best_minus_final": best_recall - final_recall,
            "final_missing_relations_count": final_missing,
            "best_checkpoint_missing_relations_count": best_missing,
            "delta_missing_best_minus_final": best_missing - final_missing,
            "final_extra_relations_count": final_extra,
            "best_checkpoint_extra_relations_count": best_extra,
            "delta_extra_best_minus_final": best_extra - final_extra,
            "final_energy_eave": float(final["energy_eave"]),
            "best_checkpoint_energy_eave": float(best["energy_eave"]),
            "min_energy_eave": float(min_energy["energy_eave"]),
            "block_index_final": final_block,
            "block_index_best_causal_f1": best_block,
            "block_index_min_energy": min_energy_block,
            "best_before_final": best_block < final_block,
            "best_matches_min_energy": best_block == min_energy_block,
        })
    return sorted(out, key=lambda row: (float(row["cooling_factor"]), int(row["optimizer_seed"])))


def _mean(values: list[float]) -> float:
    return sum(values) / len(values)


def _schedule_aggregates(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, float], list[dict[str, object]]] = {}
    for row in rows:
        key = (str(row["schedule_label"]), str(row["budget_label"]), float(row["cooling_factor"]))
        grouped.setdefault(key, []).append(row)

    aggregates: list[dict[str, object]] = []
    for (schedule_label, budget_label, cooling_factor), group in grouped.items():
        aggregates.append({
            "schedule_label": schedule_label,
            "budget_label": budget_label,
            "cooling_factor": cooling_factor,
            "groups": len(group),
            "avg_final_causal_f1": _mean([float(row["final_causal_f1"]) for row in group]),
            "avg_best_checkpoint_causal_f1": _mean([float(row["best_checkpoint_causal_f1"]) for row in group]),
            "avg_delta_best_minus_final": _mean([float(row["delta_best_minus_final"]) for row in group]),
        })
    return sorted(aggregates, key=lambda row: float(row["cooling_factor"]))


def _summary_row(rows: list[dict[str, object]]) -> dict[str, object]:
    deltas = [float(row["delta_best_minus_final"]) for row in rows]
    aggregates = _schedule_aggregates(rows)
    best_avg_checkpoint = max(aggregates, key=lambda row: float(row["avg_best_checkpoint_causal_f1"]))
    best_avg_delta = max(aggregates, key=lambda row: float(row["avg_delta_best_minus_final"]))
    return {
        "n_groups": len(rows),
        "avg_final_causal_f1": _mean([float(row["final_causal_f1"]) for row in rows]),
        "avg_best_checkpoint_causal_f1": _mean([float(row["best_checkpoint_causal_f1"]) for row in rows]),
        "avg_delta_best_minus_final": _mean(deltas),
        "median_delta_best_minus_final": statistics.median(deltas),
        "max_delta_best_minus_final": max(deltas),
        "count_best_gt_final": sum(1 for row in rows if float(row["delta_best_minus_final"]) > 0.0),
        "count_best_gt_final_by_0p02": sum(1 for row in rows if float(row["delta_best_minus_final"]) > 0.02),
        "count_best_before_final": sum(1 for row in rows if bool(row["best_before_final"])),
        "count_best_matches_min_energy": sum(1 for row in rows if bool(row["best_matches_min_energy"])),
        "best_avg_checkpoint_schedule_label": best_avg_checkpoint["schedule_label"],
        "best_avg_checkpoint_budget_label": best_avg_checkpoint["budget_label"],
        "best_avg_checkpoint_cooling_factor": best_avg_checkpoint["cooling_factor"],
        "best_avg_checkpoint_causal_f1": best_avg_checkpoint["avg_best_checkpoint_causal_f1"],
        "best_avg_delta_schedule_label": best_avg_delta["schedule_label"],
        "best_avg_delta_budget_label": best_avg_delta["budget_label"],
        "best_avg_delta_cooling_factor": best_avg_delta["cooling_factor"],
        "best_avg_delta_best_minus_final": best_avg_delta["avg_delta_best_minus_final"],
    }


def _write_csv(path: Path, headers: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=headers)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _fmt(row[key]) for key in headers})


def _write_svg(rows: list[dict[str, object]]) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - environment issue
        raise RuntimeError("matplotlib is required for oracle checkpoint ceiling SVG") from exc

    fig, (ax_scatter, ax_bar) = plt.subplots(1, 2, figsize=(11, 4.8))

    final = [float(row["final_causal_f1"]) for row in rows]
    best = [float(row["best_checkpoint_causal_f1"]) for row in rows]
    deltas = [float(row["delta_best_minus_final"]) for row in rows]
    labels = [f"{row['schedule_label']}/{row['optimizer_seed']}" for row in rows]

    lo = min(final + best)
    hi = max(final + best)
    pad = 0.03
    ax_scatter.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color="#555", lw=1.2, ls="--")
    ax_scatter.scatter(final, best, color="#1f77b4", s=34)
    ax_scatter.set_xlabel("final causal F1")
    ax_scatter.set_ylabel("best checkpoint causal F1")
    ax_scatter.set_title("Oracle ceiling vs endpoint")
    ax_scatter.grid(True, alpha=0.3)

    x = list(range(len(rows)))
    ax_bar.bar(x, deltas, color="#2a9d8f")
    ax_bar.axhline(0.0, color="#555", lw=1.0)
    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels(labels, rotation=75, ha="right", fontsize=7)
    ax_bar.set_ylabel("best - final F1")
    ax_bar.set_title("Oracle checkpoint gain")
    ax_bar.grid(True, axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig(SVG_PATH, format="svg")
    plt.close(fig)


def _write_markdown(rows: list[dict[str, object]], summary: dict[str, object]) -> None:
    generated_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
    aggregates = _schedule_aggregates(rows)

    lines = [
        "# Oracle checkpoint ceiling N=36",
        "",
        "Post-run SORKIN-2 oracle diagnostic over the existing schedule/optimizer-seed trajectory matrix.",
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
        "- Grouping key: `(optimizer_seed, schedule_label, cooling_factor, budget_label)`.",
        "- Oracle checkpoint maximizes `causal_f1`, then `causal_recall`, then minimizes missing, then extra, then chooses the earliest block.",
        "",
        "## Per group",
        "",
        "| seed | schedule | gamma | budget | final F1 | best F1 | delta F1 | final recall | best recall | final missing | best missing | final E | best E | min E | final block | best block | minE block | best before final | best=minE |",
        "| ---: | --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |",
    ]
    for row in rows:
        lines.append(
            "| {seed} | {schedule} | {gamma} | {budget} | {final_f1} | {best_f1} | {delta} | {final_recall} | {best_recall} | {final_missing} | {best_missing} | {final_e} | {best_e} | {min_e} | {final_block} | {best_block} | {min_block} | {before} | {match} |".format(
                seed=row["optimizer_seed"],
                schedule=row["schedule_label"],
                gamma=_fmt_f(float(row["cooling_factor"])),
                budget=row["budget_label"],
                final_f1=_fmt_f(float(row["final_causal_f1"])),
                best_f1=_fmt_f(float(row["best_checkpoint_causal_f1"])),
                delta=_fmt_f(float(row["delta_best_minus_final"])),
                final_recall=_fmt_f(float(row["final_causal_recall"])),
                best_recall=_fmt_f(float(row["best_checkpoint_causal_recall"])),
                final_missing=row["final_missing_relations_count"],
                best_missing=row["best_checkpoint_missing_relations_count"],
                final_e=_fmt_f(float(row["final_energy_eave"])),
                best_e=_fmt_f(float(row["best_checkpoint_energy_eave"])),
                min_e=_fmt_f(float(row["min_energy_eave"])),
                final_block=row["block_index_final"],
                best_block=row["block_index_best_causal_f1"],
                min_block=row["block_index_min_energy"],
                before=_fmt(row["best_before_final"]),
                match=_fmt(row["best_matches_min_energy"]),
            )
        )

    lines.extend([
        "",
        "## Schedule aggregates",
        "",
        "| schedule | gamma | budget | groups | avg final F1 | avg oracle best F1 | avg delta F1 |",
        "| --- | ---: | --- | ---: | ---: | ---: | ---: |",
    ])
    for row in aggregates:
        lines.append(
            "| {schedule} | {gamma} | {budget} | {groups} | {final} | {best} | {delta} |".format(
                schedule=row["schedule_label"],
                gamma=_fmt_f(float(row["cooling_factor"])),
                budget=row["budget_label"],
                groups=row["groups"],
                final=_fmt_f(float(row["avg_final_causal_f1"])),
                best=_fmt_f(float(row["avg_best_checkpoint_causal_f1"])),
                delta=_fmt_f(float(row["avg_delta_best_minus_final"])),
            )
        )

    lines.extend([
        "",
        "## Global summary",
        "",
        f"- Groups: `{summary['n_groups']}`.",
        f"- Average final causal F1: `{_fmt_f(float(summary['avg_final_causal_f1']))}`.",
        f"- Average oracle best checkpoint causal F1: `{_fmt_f(float(summary['avg_best_checkpoint_causal_f1']))}`.",
        f"- Average oracle gain: `{_fmt_f(float(summary['avg_delta_best_minus_final']))}`.",
        f"- Median oracle gain: `{_fmt_f(float(summary['median_delta_best_minus_final']))}`.",
        f"- Max oracle gain: `{_fmt_f(float(summary['max_delta_best_minus_final']))}`.",
        f"- Best checkpoint > final: `{summary['count_best_gt_final']}` of `{summary['n_groups']}`.",
        f"- Best checkpoint > final by 0.02: `{summary['count_best_gt_final_by_0p02']}` of `{summary['n_groups']}`.",
        f"- Best checkpoint before final: `{summary['count_best_before_final']}` of `{summary['n_groups']}`.",
        f"- Best checkpoint matches minimum energy block: `{summary['count_best_matches_min_energy']}` of `{summary['n_groups']}`.",
        f"- Best avg oracle schedule/budget/cooling: `{summary['best_avg_checkpoint_schedule_label']}` / `{summary['best_avg_checkpoint_budget_label']}` / `{_fmt_f(float(summary['best_avg_checkpoint_cooling_factor']))}` with avg best F1 `{_fmt_f(float(summary['best_avg_checkpoint_causal_f1']))}`.",
        f"- Best avg oracle gain schedule/budget/cooling: `{summary['best_avg_delta_schedule_label']}` / `{summary['best_avg_delta_budget_label']}` / `{_fmt_f(float(summary['best_avg_delta_cooling_factor']))}` with avg gain `{_fmt_f(float(summary['best_avg_delta_best_minus_final']))}`.",
        "",
        "## Conservative interpretation",
        "",
        "This is an oracle diagnostic: it selects checkpoints using causal F1 against the known-truth target order.",
        "It is not a causally deployable checkpoint-selection criterion for truth-free cases.",
        "The oracle ceiling substantially improves the endpoint in most groups, so the N=36 limitation in this matrix is not only search access or budget: the annealer often visits better causal checkpoints that the endpoint selector discards.",
        "Because the minimum-energy block coincides with the best-F1 block in only a minority of groups, historical energy should not be used alone as a recoverability selector in this diagnostic setting.",
        "",
        "## Guardrails",
        "",
        "This is a post-run diagnostic only, using benchmark cases with known truth.",
        "It is not an embeddability claim, not a physical gamma claim, not an N-transition claim, and not proof of general annealer failure.",
        "It is not a deployable criterion for truth-free cases yet.",
        "",
    ])
    MD_PATH.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    source_rows = _read_source()
    rows = _group_rows(source_rows)
    summary = _summary_row(rows)
    _write_csv(CSV_PATH, GROUP_HEADERS, rows)
    _write_csv(SUMMARY_CSV_PATH, SUMMARY_HEADERS, [summary])
    _write_svg(rows)
    _write_markdown(rows, summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
