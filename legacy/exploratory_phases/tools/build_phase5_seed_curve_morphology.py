#!/usr/bin/env python3
"""Phase 5 — seed-level curve morphology audit.

Phase 5 does not run new simulations.  It consumes the Phase 4B per-seed
provenance CSV and asks whether aggregate cell-level V-like/interior-minimum
morphology also exists in individual seed curves.
"""

from __future__ import annotations

import csv
import math
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools import build_phase4a_epsilon_sweep as p4a  # noqa: E402
from tools import build_phase4b_survival_probe as p4b  # noqa: E402


FOUNDATION = ROOT / "benchmarks" / "foundation"
PHASE4B_PER_SEED_CSV = FOUNDATION / "phase4b_survival_probe_per_seed.csv"
PHASE4B_CSV = FOUNDATION / "phase4b_survival_probe.csv"
PHASE4B_PER_EPSILON_CSV = FOUNDATION / "phase4b_survival_probe_per_epsilon.csv"
PHASE5_CSV = FOUNDATION / "phase5_seed_curve_morphology.csv"
PHASE5_MD = FOUNDATION / "phase5_seed_curve_morphology.md"

PHASE5_CSV_HEADERS = (
    "phase",
    "source_phase",
    "grid",
    "family",
    "n",
    "target_dim",
    "seed",
    "n_valid_epsilon",
    "n_invalid_epsilon",
    "complete_epsilon_coverage",
    "seed_curve_shape",
    "seed_censoring_label",
    "seed_floor_saturated",
    "epsilon_at_min",
    "min_loss",
    "rise_frac",
    "fall_frac",
    "tail_positive_count",
    "tail_negative_count",
    "tail_zero_count",
    "tail_n_pairs",
    "tail_positive_fraction",
    "tail_pattern",
    "ordering_fraction",
    "chain3_abundance",
    "dim_discrepancy_rel_midpoint",
    "curve_shape_cell_phase4b",
    "survival_label_cell_phase4b",
    "borderline_v_like_cell_phase4b",
)


CELL_SUMMARY_HEADERS = (
    "n",
    "target_dim",
    "n_seeds",
    "count_seed_v_shape",
    "count_seed_monotone_decay",
    "count_seed_interior_min_noisy_tail",
    "count_seed_floor_saturated",
    "count_seed_insufficient_valid_points",
    "majority_seed_shape",
    "cell_curve_shape_phase4b",
    "aggregate_representative_of_seed_majority",
)


def _parse_bool(value: str) -> bool:
    return str(value).strip().lower() == "true"


def _parse_float(value: str) -> float:
    if value in ("", "NA", None):
        return float("nan")
    return float(value)


def _format(value) -> str:
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, str):
        return value
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if math.isnan(value):
            return "NA"
        if math.isinf(value):
            return "inf" if value > 0 else "-inf"
        return f"{value:.10g}"
    return str(value)


def load_phase4b_per_seed(path: Path = PHASE4B_PER_SEED_CSV) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def load_phase4b_summary(path: Path = PHASE4B_CSV) -> list[dict]:
    return p4b.load_summary_csv(path)


def group_seed_curves(rows: list[dict]) -> dict[tuple[int, int, int], list[dict]]:
    grouped: dict[tuple[int, int, int], list[dict]] = {}
    for row in rows:
        key = (int(row["n"]), int(row["target_dim"]), int(row["seed"]))
        grouped.setdefault(key, []).append(row)
    for key in grouped:
        grouped[key] = sorted(grouped[key], key=lambda r: _parse_float(r["epsilon"]))
    return grouped


def _loss_curve(valid_rows: list[dict]) -> list[tuple[float, float, float]]:
    return [
        (_parse_float(r["epsilon"]), _parse_float(r["loss"]), _parse_float(r["loss"]))
        for r in sorted(valid_rows, key=lambda row: _parse_float(row["epsilon"]))
    ]


def classify_seed_curve(
    curve: list[tuple[float, float, float]],
    floor_tolerance: float = p4b.DEFAULT_FLOOR_TOLERANCE,
) -> tuple[str, str, bool, dict]:
    """Return (seed_curve_shape, censoring_label, floor_saturated, tail_audit)."""
    if len(curve) < 3:
        empty_tail = {
            "tail_positive_count": 0,
            "tail_negative_count": 0,
            "tail_zero_count": 0,
            "tail_n_pairs": 0,
            "tail_positive_fraction": float("nan"),
            "tail_pattern": "",
        }
        return (
            "seed_insufficient_valid_points",
            "insufficient_valid_points",
            False,
            empty_tail,
        )

    stats = p4b._curve_stats(curve)
    floor_saturated = (
        math.isfinite(float(stats["min_val"]))
        and float(stats["min_val"]) <= floor_tolerance
    )
    if floor_saturated:
        censoring_label = "floor_saturated"
        seed_shape = "seed_floor_saturated"
    else:
        censoring_label = "none"
        raw_shape = p4a._classify_curve_shape(curve, idx=1)
        tail = p4b._tail_audit(curve, raw_shape, floor_saturated)
        has_interior = bool(tail["has_interior_minimum"])
        rise_frac = float(stats["rise_frac"])
        if raw_shape == "v_shape":
            seed_shape = "seed_v_shape"
        elif (
            has_interior
            and math.isfinite(rise_frac)
            and rise_frac > 0.0
            and int(tail["tail_negative_count"]) > 0
        ):
            seed_shape = "seed_interior_min_noisy_tail"
        elif raw_shape == "monotone_decay":
            seed_shape = "seed_monotone_decay"
        elif has_interior and math.isfinite(rise_frac) and rise_frac > 0.0:
            seed_shape = "seed_interior_min_noisy_tail"
        else:
            seed_shape = "seed_flat_noisy"
        return seed_shape, censoring_label, floor_saturated, tail

    tail = p4b._tail_audit(curve, "noisy", floor_saturated)
    return seed_shape, censoring_label, floor_saturated, tail


def summarize_seed_group(
    group_rows: list[dict],
    expected_epsilon_count: int,
    floor_tolerance: float = p4b.DEFAULT_FLOOR_TOLERANCE,
) -> dict:
    first = group_rows[0]
    valid_rows = [r for r in group_rows if _parse_bool(r["valid"])]
    invalid_rows = [r for r in group_rows if not _parse_bool(r["valid"])]
    curve = _loss_curve(valid_rows)
    stats = p4b._curve_stats(curve)
    seed_shape, censoring_label, floor_saturated, tail = classify_seed_curve(
        curve, floor_tolerance=floor_tolerance
    )
    complete = (
        len(valid_rows) == expected_epsilon_count
        and len(group_rows) == expected_epsilon_count
    )
    return {
        "phase": "phase5_seed_curve_morphology",
        "source_phase": "phase4b_exploratory",
        "grid": first["grid"],
        "family": first["family"],
        "n": int(first["n"]),
        "target_dim": int(first["target_dim"]),
        "seed": int(first["seed"]),
        "n_valid_epsilon": len(valid_rows),
        "n_invalid_epsilon": len(invalid_rows),
        "complete_epsilon_coverage": complete,
        "seed_curve_shape": seed_shape,
        "seed_censoring_label": censoring_label,
        "seed_floor_saturated": floor_saturated,
        "epsilon_at_min": stats["epsilon_at_min"],
        "min_loss": stats["min_val"],
        "rise_frac": stats["rise_frac"],
        "fall_frac": stats["fall_frac"],
        "tail_positive_count": tail["tail_positive_count"],
        "tail_negative_count": tail["tail_negative_count"],
        "tail_zero_count": tail["tail_zero_count"],
        "tail_n_pairs": tail["tail_n_pairs"],
        "tail_positive_fraction": tail["tail_positive_fraction"],
        "tail_pattern": tail["tail_pattern"],
        "ordering_fraction": _parse_float(first["ordering_fraction"]),
        "chain3_abundance": _parse_float(first["chain3_abundance"]),
        "dim_discrepancy_rel_midpoint": _parse_float(
            first["dim_discrepancy_rel_midpoint"]
        ),
        "curve_shape_cell_phase4b": first["curve_shape_cell"],
        "survival_label_cell_phase4b": first["survival_label_cell"],
        "borderline_v_like_cell_phase4b": _parse_bool(
            first["borderline_v_like_cell"]
        ),
    }


def build_seed_morphology_rows(
    per_seed_rows: list[dict],
    floor_tolerance: float = p4b.DEFAULT_FLOOR_TOLERANCE,
) -> list[dict]:
    grouped = group_seed_curves(per_seed_rows)
    expected_epsilon_count = len({
        _parse_float(r["epsilon"]) for r in per_seed_rows
    })
    out = [
        summarize_seed_group(
            rows,
            expected_epsilon_count=expected_epsilon_count,
            floor_tolerance=floor_tolerance,
        )
        for _, rows in sorted(grouped.items())
    ]
    return out


def _shape_count(rows: list[dict], shape: str) -> int:
    return sum(1 for r in rows if r["seed_curve_shape"] == shape)


def _majority_shape(rows: list[dict]) -> str:
    counts = Counter(r["seed_curve_shape"] for r in rows)
    if not counts:
        return "NA"
    max_count = max(counts.values())
    winners = sorted(shape for shape, count in counts.items() if count == max_count)
    return winners[0] if len(winners) == 1 else "tie:" + "|".join(winners)


def _phase4b_shape_represented(phase4b_shape: str, majority_seed_shape: str) -> bool:
    if majority_seed_shape == "seed_v_shape":
        return phase4b_shape == "v_shape"
    if majority_seed_shape == "seed_monotone_decay":
        return phase4b_shape == "monotone_decay"
    if majority_seed_shape == "seed_flat_noisy":
        return phase4b_shape in {"flat", "noisy"}
    if majority_seed_shape == "seed_floor_saturated":
        return False
    if majority_seed_shape == "seed_interior_min_noisy_tail":
        return False
    return False


def build_cell_summary(seed_rows: list[dict]) -> list[dict]:
    by_cell: dict[tuple[int, int], list[dict]] = {}
    for row in seed_rows:
        by_cell.setdefault((row["n"], row["target_dim"]), []).append(row)
    out: list[dict] = []
    for (n, d), rows in sorted(by_cell.items()):
        majority = _majority_shape(rows)
        phase4b_shape = rows[0]["curve_shape_cell_phase4b"]
        out.append({
            "n": n,
            "target_dim": d,
            "n_seeds": len(rows),
            "count_seed_v_shape": _shape_count(rows, "seed_v_shape"),
            "count_seed_monotone_decay": _shape_count(rows, "seed_monotone_decay"),
            "count_seed_interior_min_noisy_tail": _shape_count(
                rows, "seed_interior_min_noisy_tail"
            ),
            "count_seed_floor_saturated": _shape_count(rows, "seed_floor_saturated"),
            "count_seed_insufficient_valid_points": _shape_count(
                rows, "seed_insufficient_valid_points"
            ),
            "majority_seed_shape": majority,
            "cell_curve_shape_phase4b": phase4b_shape,
            "aggregate_representative_of_seed_majority": _phase4b_shape_represented(
                phase4b_shape, majority
            ),
        })
    return out


def phase5_outcome(cell_summary: list[dict]) -> str:
    if not cell_summary:
        return "INSUFFICIENT"
    if all(r["count_seed_insufficient_valid_points"] == r["n_seeds"] for r in cell_summary):
        return "INSUFFICIENT"
    n_seeds = sum(r["n_seeds"] for r in cell_summary)
    n_censored = sum(
        r["count_seed_floor_saturated"] + r["count_seed_insufficient_valid_points"]
        for r in cell_summary
    )
    if n_seeds and n_censored / n_seeds >= 0.5:
        return "INSUFFICIENT"
    aggregate_v = [r for r in cell_summary if r["cell_curve_shape_phase4b"] == "v_shape"]
    if aggregate_v:
        v_supported = sum(1 for r in aggregate_v if r["count_seed_v_shape"] > 0)
        v_majority = sum(1 for r in aggregate_v if r["majority_seed_shape"] == "seed_v_shape")
        if v_supported == len(aggregate_v) and v_majority >= max(1, len(aggregate_v) - 1):
            return "SEED_LEVEL_SUPPORT"
    artifact_like = [
        r for r in aggregate_v
        if r["count_seed_v_shape"] == 0
        and r["count_seed_interior_min_noisy_tail"] == 0
    ]
    if aggregate_v and len(artifact_like) == len(aggregate_v):
        return "AGGREGATE_ARTIFACT"
    return "SEED_LEVEL_MIXED"


def write_csv(rows: list[dict], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=PHASE5_CSV_HEADERS)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _format(row[key]) for key in PHASE5_CSV_HEADERS})


def _md_cell_table(rows: list[dict]) -> list[str]:
    lines = [
        "| n | target_dim | seeds | seed_v_shape | seed_monotone_decay | "
        "seed_interior_min_noisy_tail | seed_floor_saturated | "
        "seed_insufficient_valid_points | majority_seed_shape | "
        "cell_curve_shape_phase4b | aggregate_represents_majority |",
        "| ---: | :---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | :---: |",
    ]
    for row in rows:
        lines.append(
            f"| {row['n']} | {row['target_dim']} | {row['n_seeds']} "
            f"| {row['count_seed_v_shape']} "
            f"| {row['count_seed_monotone_decay']} "
            f"| {row['count_seed_interior_min_noisy_tail']} "
            f"| {row['count_seed_floor_saturated']} "
            f"| {row['count_seed_insufficient_valid_points']} "
            f"| {row['majority_seed_shape']} "
            f"| {row['cell_curve_shape_phase4b']} "
            f"| {'true' if row['aggregate_representative_of_seed_majority'] else 'false'} |"
        )
    return lines


def _answer_lines(cell_summary: list[dict]) -> list[str]:
    by_cell = {(r["n"], r["target_dim"]): r for r in cell_summary}
    aggregate_v = [r for r in cell_summary if r["cell_curve_shape_phase4b"] == "v_shape"]
    v_with_seed_v = sum(1 for r in aggregate_v if r["count_seed_v_shape"] > 0)
    c483 = by_cell.get((48, 3))
    c484 = by_cell.get((48, 4))
    floor_cells = [
        r for r in cell_summary
        if r["count_seed_floor_saturated"] > 0
    ]
    lines = [
        "## Seed-level questions",
        "",
        f"- Aggregate V-shapes with at least one seed-level V-shape: {v_with_seed_v}/{len(aggregate_v)}.",
    ]
    if c483:
        lines.append(
            "- `(48,3)` remains a seed-level mixed/borderline audit target: "
            f"{c483['count_seed_v_shape']} seed-level V-shapes, "
            f"{c483['count_seed_interior_min_noisy_tail']} interior-minimum noisy-tail seeds, "
            f"majority `{c483['majority_seed_shape']}`."
        )
    if c484:
        lines.append(
            "- `(48,4)` remains a counterexample/noisy-tail audit target: "
            f"{c484['count_seed_v_shape']} seed-level V-shapes, "
            f"{c484['count_seed_interior_min_noisy_tail']} interior-minimum noisy-tail seeds, "
            f"majority `{c484['majority_seed_shape']}`."
        )
    lines.append(
        f"- Cells with floor-saturated seeds: {len(floor_cells)}. These seeds are censored and are not counted as strong negative evidence against V-like behavior."
    )
    return lines


def _count_threshold(rows: list[dict], threshold: float) -> int:
    return sum(
        1 for r in rows
        if math.isfinite(float(r["min_loss"])) and float(r["min_loss"]) <= threshold
    )


def _floor_audit_table(seed_rows: list[dict]) -> list[str]:
    by_cell: dict[tuple[int, int], list[dict]] = {}
    for row in seed_rows:
        by_cell.setdefault((row["n"], row["target_dim"]), []).append(row)
    lines = [
        "| n | target_dim | n_seeds_total | n_min_loss_eq_0 | n_le_1e-15 | n_le_1e-12 | n_le_1e-9 | n_le_1e-6 | n_le_1e-4 | min(min_loss) | median(min_loss) | max(min_loss) |",
        "| ---: | :---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for cell in sorted(by_cell):
        rows = by_cell[cell]
        vals = [
            float(r["min_loss"])
            for r in rows
            if math.isfinite(float(r["min_loss"]))
        ]
        vals.sort()
        median = vals[len(vals) // 2] if len(vals) % 2 == 1 else (
            (vals[len(vals) // 2 - 1] + vals[len(vals) // 2]) / 2.0
            if vals else float("nan")
        )
        lines.append(
            f"| {cell[0]} | {cell[1]} | {len(rows)} "
            f"| {sum(1 for v in vals if v == 0.0)} "
            f"| {_count_threshold(rows, 1e-15)} "
            f"| {_count_threshold(rows, 1e-12)} "
            f"| {_count_threshold(rows, 1e-9)} "
            f"| {_count_threshold(rows, 1e-6)} "
            f"| {_count_threshold(rows, 1e-4)} "
            f"| {_format(min(vals) if vals else float('nan'))} "
            f"| {_format(median)} "
            f"| {_format(max(vals) if vals else float('nan'))} |"
        )
    return lines


def _floor_audit_lines(seed_rows: list[dict]) -> list[str]:
    valid_with_min = [
        r for r in seed_rows if math.isfinite(float(r["min_loss"]))
    ]
    n_valid = len(valid_with_min)
    n_eq0 = sum(1 for r in valid_with_min if float(r["min_loss"]) == 0.0)
    lines = [
        "## Floor-censoring audit",
        "",
        *_floor_audit_table(seed_rows),
        "",
        f"- Global valid seed curves with finite `min_loss`: {n_valid}.",
        f"- Exact-zero minima: {n_eq0}/{n_valid}.",
        f"- `min_loss <= 1e-15`: {_count_threshold(valid_with_min, 1e-15)}/{n_valid}.",
        f"- `min_loss <= 1e-12`: {_count_threshold(valid_with_min, 1e-12)}/{n_valid}.",
        f"- `min_loss <= 1e-9`: {_count_threshold(valid_with_min, 1e-9)}/{n_valid}.",
        f"- `min_loss <= 1e-6`: {_count_threshold(valid_with_min, 1e-6)}/{n_valid}.",
        f"- `min_loss <= 1e-4`: {_count_threshold(valid_with_min, 1e-4)}/{n_valid}.",
        "",
        "Seed-level morphology is currently censored by optimizer-floor saturation; aggregate Phase 4B curves should not be interpreted as direct evidence of seed-level V-shapes.",
        "",
        "In the current pilot grid this censoring is dominated by exact zeros rather than by sensitivity to the chosen `floor_tolerance` threshold.",
        "",
    ]
    return lines


def write_markdown(
    seed_rows: list[dict],
    cell_summary: list[dict],
    outcome: str,
    path: Path,
) -> None:
    lines = [
        "# Phase 5 — Seed-level curve morphology audit",
        "",
        "**Status:** exploratory seed-level audit using existing Phase 4B CSV provenance only. No new simulations, PySR, BDG action, SVG generation, or Phase 4B relabeling are performed.",
        "",
        "## Objective",
        "",
        "Phase 5 asks whether V-like or interior-minimum morphology observed in aggregate Phase 4B cell curves also appears in individual seed curves, or whether it is an artifact of averaging heterogeneous seeds.",
        "",
        "## Inputs",
        "",
        "- `phase4b_survival_probe_per_seed.csv`: reconstructs `loss(epsilon)` by `(n, target_dim, seed)`.",
        "- `phase4b_survival_probe.csv`: supplies aggregate Phase 4B cell labels for comparison.",
        "- `phase4b_survival_probe_per_epsilon.csv`: documents aggregate per-epsilon provenance.",
        "",
        "## Seed-level classification rule",
        "",
        "For each `(n, target_dim, seed)`, Phase 5 sorts valid rows by epsilon and classifies the individual `loss(epsilon)` curve. If fewer than three valid epsilons are available, the seed is marked `seed_insufficient_valid_points`.",
        "",
        f"Floor saturation uses the Phase 4B floor tolerance `{p4b.DEFAULT_FLOOR_TOLERANCE:g}`. A floor-saturated seed is censored and is not counted as strong negative evidence against V-like behavior.",
        "",
        "The shape taxonomy is `seed_v_shape`, `seed_monotone_decay`, `seed_floor_saturated`, `seed_interior_min_noisy_tail`, `seed_flat_noisy`, and `seed_insufficient_valid_points`.",
        "",
        "## Cell-level summary",
        "",
        *_md_cell_table(cell_summary),
        "",
        *_answer_lines(cell_summary),
        "",
        "## Global Phase 5 outcome",
        "",
        f"**{outcome}**",
        "",
        "Outcome definitions:",
        "",
        "- `SEED_LEVEL_SUPPORT`: aggregate V-like morphology is broadly visible in individual seed curves.",
        "- `SEED_LEVEL_MIXED`: seed-level morphology is present but heterogeneous across cells or seeds.",
        "- `AGGREGATE_ARTIFACT`: aggregate V-like morphology is not visible at seed level.",
        "- `INSUFFICIENT`: available seed curves are not sufficient for the audit.",
        "",
        *_floor_audit_lines(seed_rows),
        "## Conservative conclusion",
        "",
        "Phase 5 distinguishes seed-level morphology from aggregate morphology under the Phase 4B optimizer-response loss. It does not establish a physical law, validate Phase 4B as a physical claim, or replace the existing Phase 4B `MIXED` outcome.",
        "",
        f"Rows written: {len(seed_rows)} seed curves across {len(cell_summary)} pilot-grid cells.",
        "",
        "Source: `tools/build_phase5_seed_curve_morphology.py`.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    if not PHASE4B_PER_SEED_CSV.exists():
        sys.exit(f"Missing input: {PHASE4B_PER_SEED_CSV}")
    if not PHASE4B_CSV.exists():
        sys.exit(f"Missing input: {PHASE4B_CSV}")
    if not PHASE4B_PER_EPSILON_CSV.exists():
        sys.exit(f"Missing input: {PHASE4B_PER_EPSILON_CSV}")
    per_seed_rows = load_phase4b_per_seed(PHASE4B_PER_SEED_CSV)
    seed_rows = build_seed_morphology_rows(per_seed_rows)
    cell_summary = build_cell_summary(seed_rows)
    outcome = phase5_outcome(cell_summary)
    write_csv(seed_rows, PHASE5_CSV)
    write_markdown(seed_rows, cell_summary, outcome, PHASE5_MD)
    print(f"Phase 5 CSV: {PHASE5_CSV}")
    print(f"Phase 5 Markdown: {PHASE5_MD}")
    print(f"Seed curves: {len(seed_rows)}")
    print(f"Global outcome: {outcome}")


if __name__ == "__main__":
    main()
