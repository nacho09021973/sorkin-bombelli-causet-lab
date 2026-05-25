#!/usr/bin/env python3
"""N=20 manufactured known-truth validation for the historical annealer.

This probe runs the unchanged Bombelli annealer on one deterministic
Minkowski sprinkling and audits each block through block_callback.
"""

from __future__ import annotations

import contextlib
import csv
import io
import math
import os
import statistics
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
TRAJECTORY_CSV = OUT_DIR / "n20_known_truth_validation_trajectory.csv"
SUMMARY_CSV = OUT_DIR / "n20_known_truth_validation_summary.csv"
MD_PATH = OUT_DIR / "n20_known_truth_validation.md"
SVG_PATH = OUT_DIR / "n20_known_truth_validation.svg"
COMMAND = "python3 explore/n20_known_truth_validation/run_n20_known_truth_validation.py"

KNOWN_TRUTH_DIR = ROOT / "explore" / "known_truth_n20"

FAMILY = "minkowski"
D_SPACETIME = 2
N = 20
CASE_SEED = 1959
TARGET_DIM = D_SPACETIME - 1
INITIAL_TEMP = 100.0
SCHEDULE_LABEL = "gamma_0p5"
COOLING_FACTOR = 0.5
BUDGET_LABEL = "medium_25_25_8"
WARMUP_LIMIT = 25
ANNEAL_LIMIT = 25
MAX_DATA = 8
BACKEND = "cpu"
OPTIMIZER_SEEDS = (1959, 1962, 1987, 2001)

TRAJECTORY_HEADERS = (
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
    "induced_pairs_str",
)

SUMMARY_HEADERS = (
    "n_runs",
    "avg_final_causal_f1",
    "avg_best_checkpoint_causal_f1",
    "count_exact_match_final",
    "count_exact_match_any_checkpoint",
    "count_best_gt_final",
    "count_best_matches_min_energy",
    "avg_delta_best_minus_final",
    "max_delta_best_minus_final",
)


def _fmt(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return f"{value:.10g}" if math.isfinite(value) else "NA"
    return str(value)


def _fmt_f(value: float) -> str:
    return f"{value:.6g}" if math.isfinite(value) else "NA"


def _safe_div(num: float, den: float) -> float:
    return 0.0 if den == 0.0 else num / den


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


def _causal_metrics_with_pairs(
    target: list[list[bool]],
    coords: list[tuple[float, ...]],
) -> dict[str, object]:
    induced = vs.induced_order_from_coords(coords)
    induced_pairs = _pairs(induced)
    target_pairs = _pairs(target)
    comparison = vs.compare_causal_orders(target, induced)
    correct = len(target_pairs & induced_pairs)
    precision = _safe_div(correct, comparison.total_relations_induced)
    recall = _safe_div(correct, comparison.total_relations_target)
    f1 = _safe_div(2.0 * precision * recall, precision + recall)
    return {
        "causal_precision": precision,
        "causal_recall": recall,
        "causal_f1": f1,
        "missing_relations_count": len(comparison.missing_relations),
        "extra_relations_count": len(comparison.extra_relations),
        "total_relations_target": comparison.total_relations_target,
        "total_relations_induced": comparison.total_relations_induced,
        "correct_relations": correct,
        "exact_match": comparison.exact_match,
        "success_flag": comparison.exact_match,
        "_induced_pairs": induced_pairs,
    }


def _run_one(case: vs.SprinkleCase, optimizer_seed: int) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []

    def _block_callback(sim: Any, block_idx: int, temp: float, eave: float) -> None:
        # sim.rold/sim.xold are the accepted current state used by the existing
        # known-truth validation path; copy them before later blocks mutate state.
        coords = [
            (float(sim.rold[i]), *[float(value) for value in sim.xold[i]])
            for i in range(sim.n)
        ]
        metrics = _causal_metrics_with_pairs(case.matrix, coords)
        induced_pairs_str = "|".join(
            f"{i}:{j}" for i, j in sorted(metrics["_induced_pairs"])
        )
        rows.append({
            "optimizer_seed": optimizer_seed,
            "schedule_label": SCHEDULE_LABEL,
            "cooling_factor": COOLING_FACTOR,
            "budget_label": BUDGET_LABEL,
            "block_index": block_idx,
            "temperature": temp,
            "energy_eave": eave,
            "causal_precision": metrics["causal_precision"],
            "causal_recall": metrics["causal_recall"],
            "causal_f1": metrics["causal_f1"],
            "missing_relations_count": metrics["missing_relations_count"],
            "extra_relations_count": metrics["extra_relations_count"],
            "total_relations_target": metrics["total_relations_target"],
            "total_relations_induced": metrics["total_relations_induced"],
            "correct_relations": metrics["correct_relations"],
            "exact_match": metrics["exact_match"],
            "success_flag": metrics["success_flag"],
            "induced_pairs_str": induced_pairs_str,
        })

    with tempfile.TemporaryDirectory() as tmpdir, contextlib.redirect_stdout(io.StringIO()):
        sim = cones.ConesSimulator(
            z=case.matrix,
            dim=TARGET_DIM,
            seed=optimizer_seed,
            interactive=False,
            max_data=MAX_DATA,
            plot_path=None,
            warmup_limit=WARMUP_LIMIT,
            anneal_limit=ANNEAL_LIMIT,
            initial_temp=INITIAL_TEMP,
            cooling_factor=COOLING_FACTOR,
            backend=BACKEND,
            block_callback=_block_callback,
        )
        sim.run(Path(tmpdir) / "annealer_output.txt")
    return rows


def _best_checkpoint(rows: list[dict[str, object]]) -> dict[str, object]:
    return max(
        rows,
        key=lambda row: (
            float(row["causal_f1"]),
            float(row["causal_recall"]),
            -int(row["missing_relations_count"]),
            -int(row["extra_relations_count"]),
            -int(row["block_index"]),
        ),
    )


def _per_seed_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[int, list[dict[str, object]]] = {}
    for row in rows:
        grouped.setdefault(int(row["optimizer_seed"]), []).append(row)

    out: list[dict[str, object]] = []
    for seed, seed_rows in sorted(grouped.items()):
        ordered = sorted(seed_rows, key=lambda row: int(row["block_index"]))
        final = ordered[-1]
        best = _best_checkpoint(ordered)
        min_energy = min(ordered, key=lambda row: float(row["energy_eave"]))
        out.append({
            "optimizer_seed": seed,
            "final_block": int(final["block_index"]),
            "final_temperature": float(final["temperature"]),
            "final_causal_f1": float(final["causal_f1"]),
            "final_recall": float(final["causal_recall"]),
            "final_missing": int(final["missing_relations_count"]),
            "final_extra": int(final["extra_relations_count"]),
            "final_energy_eave": float(final["energy_eave"]),
            "final_exact_match": bool(final["exact_match"]),
            "best_block": int(best["block_index"]),
            "best_temperature": float(best["temperature"]),
            "best_causal_f1": float(best["causal_f1"]),
            "best_recall": float(best["causal_recall"]),
            "best_missing": int(best["missing_relations_count"]),
            "best_extra": int(best["extra_relations_count"]),
            "best_energy_eave": float(best["energy_eave"]),
            "best_exact_match": bool(best["exact_match"]),
            "delta_best_minus_final": float(best["causal_f1"]) - float(final["causal_f1"]),
            "any_exact_match": any(bool(row["exact_match"]) for row in ordered),
            "best_matches_min_energy": int(best["block_index"]) == int(min_energy["block_index"]),
            "min_energy_block": int(min_energy["block_index"]),
            "min_energy_eave": float(min_energy["energy_eave"]),
        })
    return out


def _global_summary(per_seed: list[dict[str, object]]) -> dict[str, object]:
    deltas = [float(row["delta_best_minus_final"]) for row in per_seed]
    return {
        "n_runs": len(per_seed),
        "avg_final_causal_f1": statistics.fmean(float(row["final_causal_f1"]) for row in per_seed),
        "avg_best_checkpoint_causal_f1": statistics.fmean(float(row["best_causal_f1"]) for row in per_seed),
        "count_exact_match_final": sum(1 for row in per_seed if bool(row["final_exact_match"])),
        "count_exact_match_any_checkpoint": sum(1 for row in per_seed if bool(row["any_exact_match"])),
        "count_best_gt_final": sum(1 for row in per_seed if float(row["delta_best_minus_final"]) > 0.0),
        "count_best_matches_min_energy": sum(1 for row in per_seed if bool(row["best_matches_min_energy"])),
        "avg_delta_best_minus_final": statistics.fmean(deltas),
        "max_delta_best_minus_final": max(deltas),
    }


def _write_csv(path: Path, headers: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=headers)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _fmt(row[key]) for key in headers})


def _write_svg(per_seed: list[dict[str, object]]) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-sorkin")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.0, 5.5))
    final_values = [float(row["final_causal_f1"]) for row in per_seed]
    best_values = [float(row["best_causal_f1"]) for row in per_seed]
    lo = min(final_values + best_values + [0.0])
    hi = max(final_values + best_values + [1.0])
    ax.plot([lo, hi], [lo, hi], color="#777777", linewidth=1.2, linestyle="--", label="endpoint = best")
    ax.scatter(final_values, best_values, s=58, color="#2a9d8f", edgecolor="#222222", zorder=3)
    for row in per_seed:
        ax.annotate(
            str(row["optimizer_seed"]),
            (float(row["final_causal_f1"]), float(row["best_causal_f1"])),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8,
        )
    ax.set_xlabel("final causal F1")
    ax.set_ylabel("best checkpoint causal F1")
    ax.set_title("N=20 known-truth validation: endpoint vs best checkpoint")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.grid(True, alpha=0.25)
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(SVG_PATH, format="svg")
    plt.close(fig)


def _write_markdown(
    rows: list[dict[str, object]],
    per_seed: list[dict[str, object]],
    summary: dict[str, object],
    runtime_seconds: float,
) -> None:
    generated_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
    target_relations = int(rows[0]["total_relations_target"]) if rows else 0
    any_exact = int(summary["count_exact_match_any_checkpoint"])
    final_exact = int(summary["count_exact_match_final"])
    best_gt_final = int(summary["count_best_gt_final"])
    best_matches_min = int(summary["count_best_matches_min_energy"])

    if any_exact > 0:
        recoverability = "At least one optimizer seed reaches exact causal recovery in the audited trajectory."
    elif float(summary["avg_best_checkpoint_causal_f1"]) >= 0.9:
        recoverability = "The case is recovered with high average causal fidelity, but no exact match appears in this run matrix."
    else:
        recoverability = "This budget does not show high-fidelity recovery for this N=20 case."

    if best_gt_final > 0:
        endpoint_readout = "The endpoint does not always preserve the best causal checkpoint."
    else:
        endpoint_readout = "The endpoint preserves the best causal checkpoint in this four-seed matrix."

    if best_matches_min == len(per_seed):
        energy_readout = "Minimum energy coincides with best causal F1 for every optimizer seed in this matrix."
    else:
        energy_readout = "Minimum energy does not consistently select the best causal F1 checkpoint."

    n36_note = (
        "Compared with the N=36 probes, this N=20 validation is smaller and uses a single gamma, "
        "so qualitative comparison is limited to endpoint/checkpoint behavior rather than a schedule matrix."
    )

    lines = [
        "# N=20 known-truth validation",
        "",
        "This is a manufactured known-truth SORKIN-2 validation case. It does not use a community golden table.",
        "",
        "## Configuration",
        "",
        f"- Command: `{COMMAND}`",
        f"- Generated at UTC: `{generated_at}`",
        f"- Runtime seconds: `{runtime_seconds:.3f}`",
        f"- Known-truth directory: `{KNOWN_TRUTH_DIR.relative_to(ROOT)}`",
        f"- family: `{FAMILY}`",
        f"- N: `{N}`",
        f"- d_spacetime: `{D_SPACETIME}`",
        f"- case_seed: `{CASE_SEED}`",
        f"- optimizer_seeds: `{', '.join(str(seed) for seed in OPTIMIZER_SEEDS)}`",
        f"- schedule_label: `{SCHEDULE_LABEL}`",
        f"- cooling_factor: `{COOLING_FACTOR}`",
        f"- T0: `{INITIAL_TEMP}`",
        f"- budget_label: `{BUDGET_LABEL}`",
        f"- warmup_limit: `{WARMUP_LIMIT}`",
        f"- anneal_limit: `{ANNEAL_LIMIT}`",
        f"- max_data: `{MAX_DATA}`",
        f"- total_relations_target: `{target_relations}`",
        "- Instrumentation: read-only `block_callback` over `sim.rold`/`sim.xold`.",
        "",
        "## Per-seed checkpoint analysis",
        "",
        "| optimizer_seed | final block | final F1 | final recall | final missing | final extra | final exact | best block | best temp | best F1 | best recall | best missing | best extra | any exact | delta best-final | min energy block | best=min energy |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | --- |",
    ]
    for row in per_seed:
        lines.append(
            "| {seed} | {fb} | {ff1} | {fr} | {fm} | {fe} | {fx} | {bb} | {bt} | {bf1} | {br} | {bm} | {be} | {ae} | {delta} | {meb} | {bme} |".format(
                seed=row["optimizer_seed"],
                fb=row["final_block"],
                ff1=_fmt_f(float(row["final_causal_f1"])),
                fr=_fmt_f(float(row["final_recall"])),
                fm=row["final_missing"],
                fe=row["final_extra"],
                fx=_fmt(row["final_exact_match"]),
                bb=row["best_block"],
                bt=_fmt_f(float(row["best_temperature"])),
                bf1=_fmt_f(float(row["best_causal_f1"])),
                br=_fmt_f(float(row["best_recall"])),
                bm=row["best_missing"],
                be=row["best_extra"],
                ae=_fmt(row["any_exact_match"]),
                delta=_fmt_f(float(row["delta_best_minus_final"])),
                meb=row["min_energy_block"],
                bme=_fmt(row["best_matches_min_energy"]),
            )
        )

    lines.extend([
        "",
        "## Global summary",
        "",
        "| metric | value |",
        "| --- | ---: |",
        f"| n_runs | {summary['n_runs']} |",
        f"| avg_final_causal_f1 | {_fmt_f(float(summary['avg_final_causal_f1']))} |",
        f"| avg_best_checkpoint_causal_f1 | {_fmt_f(float(summary['avg_best_checkpoint_causal_f1']))} |",
        f"| count_exact_match_final | {summary['count_exact_match_final']} |",
        f"| count_exact_match_any_checkpoint | {summary['count_exact_match_any_checkpoint']} |",
        f"| count_best_gt_final | {summary['count_best_gt_final']} |",
        f"| count_best_matches_min_energy | {summary['count_best_matches_min_energy']} |",
        f"| avg_delta_best_minus_final | {_fmt_f(float(summary['avg_delta_best_minus_final']))} |",
        f"| max_delta_best_minus_final | {_fmt_f(float(summary['max_delta_best_minus_final']))} |",
        "",
        "## Readout",
        "",
        f"1. N=20 recoverability: {recoverability}",
        f"2. Exact match: endpoint exact matches occur in `{final_exact}/{len(per_seed)}` runs; any-checkpoint exact matches occur in `{any_exact}/{len(per_seed)}` runs.",
        f"3. Endpoint selection: {endpoint_readout}",
        f"4. Energy selection: {energy_readout}",
        f"5. N=20 vs N=36: {n36_note}",
        "6. Implication for N=36: if N=20 preserves best checkpoints more often than N=36, the N=36 selection/parada hypothesis is strengthened; if not, checkpoint selection remains a broader annealer diagnostic rather than an N=36-only effect.",
        "",
        "## Guardrails",
        "",
        "This is one manufactured known-truth causal set, one family, one case_seed, one schedule, four optimizer seeds, and the unchanged historical Bombelli objective.",
        "It is not a community benchmark, not an embeddability claim, not a physical gamma claim, not an N-transition claim, and not proof of general annealer failure.",
        "",
    ])
    MD_PATH.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    start = time.perf_counter()
    case = _make_case()
    rows: list[dict[str, object]] = []
    for optimizer_seed in OPTIMIZER_SEEDS:
        rows.extend(_run_one(case, optimizer_seed))
    per_seed = _per_seed_rows(rows)
    summary = _global_summary(per_seed)
    _write_csv(TRAJECTORY_CSV, TRAJECTORY_HEADERS, rows)
    _write_csv(SUMMARY_CSV, SUMMARY_HEADERS, [summary])
    _write_svg(per_seed)
    _write_markdown(rows, per_seed, summary, time.perf_counter() - start)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
