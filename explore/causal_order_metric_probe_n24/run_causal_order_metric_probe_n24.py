#!/usr/bin/env python3
"""Minimal causal-order metric probe for SORKIN-2.

Exploratory only: compares final energy, interval RMSE, and direct
causal-order preservation for three N=24 gamma values.
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


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import causet_invariants  # noqa: E402
import cones  # noqa: E402
import validation_suite as vs  # noqa: E402


OUT_DIR = Path(__file__).resolve().parent
CSV_PATH = OUT_DIR / "causal_order_metric_probe_n24.csv"
MD_PATH = OUT_DIR / "causal_order_metric_probe_n24.md"
SVG_PATH = OUT_DIR / "causal_order_metric_probe_n24.svg"
COMMAND = "python3 explore/causal_order_metric_probe_n24/run_causal_order_metric_probe_n24.py"

FAMILY = "minkowski"
D_SPACETIME = 2
N = 24
CASE_SEED = 1959
OPTIMIZER_SEED = 1987
TARGET_DIM = D_SPACETIME - 1
INITIAL_TEMP = 100.0
WARMUP_LIMIT = 10
ANNEAL_LIMIT = 10
MAX_DATA = 4
BACKEND = "cpu"
SUCCESS_GAP_THRESHOLD = 1.0
GAMMAS = (0.55, 0.73, 0.82)

CSV_HEADERS = (
    "gamma",
    "n",
    "initial_temp",
    "warmup_limit",
    "anneal_limit",
    "max_data",
    "final_energy",
    "energy_gap",
    "interval_rmse",
    "mm_dim_truth",
    "mm_dim_recovered",
    "success_flag",
    "total_relations_target",
    "total_relations_induced",
    "correct_relations",
    "missing_relations_count",
    "extra_relations_count",
    "exact_match",
    "causal_precision",
    "causal_recall",
    "causal_f1",
    "false_positive_relation_rate",
    "false_negative_relation_rate",
    "comparability_fraction_error",
    "runtime_seconds",
)


def _fmt(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if math.isfinite(value):
            return f"{value:.10g}"
        return "NA"
    return str(value)


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


def _safe_div(numerator: float, denominator: float) -> float:
    if denominator == 0.0:
        return 0.0
    return numerator / denominator


def _run_one(case: vs.SprinkleCase, gamma: float) -> dict[str, object]:
    start = time.perf_counter()
    with tempfile.TemporaryDirectory() as tmpdir, contextlib.redirect_stdout(io.StringIO()):
        sim = cones.ConesSimulator(
            z=case.matrix,
            dim=TARGET_DIM,
            seed=OPTIMIZER_SEED,
            interactive=False,
            max_data=MAX_DATA,
            plot_path=None,
            warmup_limit=WARMUP_LIMIT,
            anneal_limit=ANNEAL_LIMIT,
            initial_temp=INITIAL_TEMP,
            cooling_factor=gamma,
            backend=BACKEND,
        )
        sim.run(Path(tmpdir) / "annealer_output.txt")
    runtime = time.perf_counter() - start

    recovered = [(sim.rnew[i], *sim.xnew[i]) for i in range(sim.n)]
    rave_truth = sum(p[0] for p in case.points) / case.n
    rave_recovered = sum(p[0] for p in recovered) / sim.n
    if rave_truth > 0.0 and rave_recovered > 0.0:
        scale = rave_recovered / rave_truth
        truth_scaled = [tuple(c * scale for c in p) for p in case.points]
    else:
        truth_scaled = list(case.points)

    interval_residual = vs.interval_rmse(recovered, truth_scaled)
    truth_energy = vs.bombelli_energy_at(case.matrix, truth_scaled, d_spatial=TARGET_DIM)
    final_energy = sim.data[-1][1] if sim.data else float("nan")
    energy_gap = final_energy - truth_energy
    success = (
        math.isfinite(energy_gap)
        and math.isfinite(interval_residual)
        and energy_gap <= SUCCESS_GAP_THRESHOLD
    )

    coords = [(sim.rold[i], *sim.xold[i]) for i in range(sim.n)]
    induced = vs.induced_order_from_coords(coords)
    comparison = vs.compare_causal_orders(case.matrix, induced)
    target_pairs = {
        (i, j)
        for i in range(comparison.n - 1)
        for j in range(i + 1, comparison.n)
        if case.matrix[i][j]
    }
    induced_pairs = {
        (i, j)
        for i in range(comparison.n - 1)
        for j in range(i + 1, comparison.n)
        if induced[i][j]
    }
    correct_relations = len(target_pairs & induced_pairs)
    missing_count = len(comparison.missing_relations)
    extra_count = len(comparison.extra_relations)
    total_pairs = comparison.n * (comparison.n - 1) // 2
    non_target_pairs = total_pairs - comparison.total_relations_target

    precision = _safe_div(correct_relations, comparison.total_relations_induced)
    recall = _safe_div(correct_relations, comparison.total_relations_target)
    f1 = _safe_div(2.0 * precision * recall, precision + recall)
    fp_rate = _safe_div(extra_count, non_target_pairs)
    fn_rate = _safe_div(missing_count, comparison.total_relations_target)
    comparability_error = (
        (comparison.total_relations_induced - comparison.total_relations_target)
        / total_pairs
    )

    return {
        "gamma": gamma,
        "n": N,
        "initial_temp": INITIAL_TEMP,
        "warmup_limit": WARMUP_LIMIT,
        "anneal_limit": ANNEAL_LIMIT,
        "max_data": MAX_DATA,
        "final_energy": final_energy,
        "energy_gap": energy_gap,
        "interval_rmse": interval_residual,
        "mm_dim_truth": float(case.d_spacetime),
        "mm_dim_recovered": causet_invariants.myrheim_meyer_dimension(case.matrix),
        "success_flag": success,
        "total_relations_target": comparison.total_relations_target,
        "total_relations_induced": comparison.total_relations_induced,
        "correct_relations": correct_relations,
        "missing_relations_count": missing_count,
        "extra_relations_count": extra_count,
        "exact_match": comparison.exact_match,
        "causal_precision": precision,
        "causal_recall": recall,
        "causal_f1": f1,
        "false_positive_relation_rate": fp_rate,
        "false_negative_relation_rate": fn_rate,
        "comparability_fraction_error": comparability_error,
        "runtime_seconds": runtime,
    }


def write_csv(rows: list[dict[str, object]]) -> None:
    with CSV_PATH.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(CSV_HEADERS)
        for row in rows:
            writer.writerow([_fmt(row[header]) for header in CSV_HEADERS])


def _best(rows: list[dict[str, object]], metric: str, *, maximize: bool = False) -> dict[str, object]:
    return (max if maximize else min)(rows, key=lambda row: float(row[metric]))


def write_svg(rows: list[dict[str, object]]) -> None:
    width = 980
    height = 620
    margin = 66
    panel_gap = 38
    panel_width = (width - 2 * margin - 2 * panel_gap) / 3
    panel_height = height - 2 * margin
    gammas = [float(row["gamma"]) for row in rows]

    panels = [
        ("final energy", "final_energy", False, margin),
        ("log10 interval RMSE", "interval_rmse", False, margin + panel_width + panel_gap),
        ("causal F1", "causal_f1", True, margin + 2 * (panel_width + panel_gap)),
    ]

    def sx(gamma: float, x0: float) -> float:
        if min(gammas) == max(gammas):
            return x0 + panel_width / 2
        return x0 + (gamma - min(gammas)) * panel_width / (max(gammas) - min(gammas))

    def sy(value: float, values: list[float], maximize: bool) -> float:
        low = min(values)
        high = max(values)
        if low == high:
            low -= 1.0
            high += 1.0
        y = height - margin - (value - low) * panel_height / (high - low)
        return y

    chunks = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "  <rect width='100%' height='100%' fill='#f7f4ed'/>",
        f"  <text x='{width / 2:.0f}' y='30' text-anchor='middle' font-family='serif' font-size='22' fill='#222'>N=24 causal-order metric probe</text>",
    ]
    for title, metric, maximize, x0 in panels:
        values = [
            math.log10(float(row[metric])) if metric == "interval_rmse" else float(row[metric])
            for row in rows
        ]
        best_row = _best(rows, metric, maximize=maximize)
        points = " ".join(
            f"{sx(float(row['gamma']), x0):.2f},{sy(math.log10(float(row[metric])) if metric == 'interval_rmse' else float(row[metric]), values, maximize):.2f}"
            for row in rows
        )
        chunks.append(f"  <text x='{x0 + panel_width / 2:.1f}' y='58' text-anchor='middle' font-family='monospace' font-size='13'>{title}</text>")
        chunks.append(f"  <line x1='{x0:.1f}' y1='{height - margin}' x2='{x0 + panel_width:.1f}' y2='{height - margin}' stroke='#333' stroke-width='2'/>")
        chunks.append(f"  <line x1='{x0:.1f}' y1='{margin}' x2='{x0:.1f}' y2='{height - margin}' stroke='#333' stroke-width='2'/>")
        chunks.append(f"  <polyline fill='none' stroke='#1f77b4' stroke-width='2.5' points='{points}'/>")
        for row in rows:
            raw_value = math.log10(float(row[metric])) if metric == "interval_rmse" else float(row[metric])
            fill = "#d1495b" if row is best_row else "#2a9d8f"
            chunks.append(
                f"  <circle cx='{sx(float(row['gamma']), x0):.2f}' cy='{sy(raw_value, values, maximize):.2f}' r='5' fill='{fill}'/>"
            )
            chunks.append(
                f"  <text x='{sx(float(row['gamma']), x0):.2f}' y='{height - margin + 22}' text-anchor='middle' font-family='monospace' font-size='10'>{float(row['gamma']):.2f}</text>"
            )
    chunks.append("  <text x='66' y='600' font-family='monospace' font-size='12' fill='#444'>red marks best by panel metric; lower is better except causal F1</text>")
    chunks.append("</svg>")
    SVG_PATH.write_text("\n".join(chunks) + "\n", encoding="utf-8")


def write_markdown(rows: list[dict[str, object]], generated_at_utc: str, total_runtime: float) -> None:
    best_energy = _best(rows, "final_energy")
    best_rmse = _best(rows, "interval_rmse")
    best_f1 = _best(rows, "causal_f1", maximize=True)
    f1_aligns = (
        "final_energy"
        if float(best_f1["gamma"]) == float(best_energy["gamma"])
        else "interval_rmse"
        if float(best_f1["gamma"]) == float(best_rmse["gamma"])
        else "neither"
    )
    total_missing = sum(int(row["missing_relations_count"]) for row in rows)
    total_extra = sum(int(row["extra_relations_count"]) for row in rows)
    dominant = "missing" if total_missing > total_extra else "extra" if total_extra > total_missing else "balanced"
    any_exact = any(bool(row["exact_match"]) for row in rows)
    any_success = any(bool(row["success_flag"]) for row in rows)

    lines = [
        "# N=24 Causal-Order Metric Probe",
        "",
        "**Status:** exploratory only; not confirmation.",
        "",
        "This probe compares direct causal-order preservation with final energy and interval RMSE for three gamma values. It is an accessibility/recoverability diagnostic only.",
        "",
        "## Provenance",
        "",
        f"- Output directory: `{OUT_DIR.relative_to(ROOT)}`",
        f"- Runner: `{Path(__file__).relative_to(ROOT)}`",
        f"- Command: `{COMMAND}`",
        f"- Generated at UTC: `{generated_at_utc}`",
        f"- CSV: `{CSV_PATH.relative_to(ROOT)}`",
        f"- Figure: `{SVG_PATH.relative_to(ROOT)}`",
        f"- N: `{N}`",
        f"- T0: `{INITIAL_TEMP}`",
        f"- warmup_limit: `{WARMUP_LIMIT}`",
        f"- anneal_limit: `{ANNEAL_LIMIT}`",
        f"- max_data: `{MAX_DATA}`",
        f"- total runtime seconds: `{total_runtime:.3f}`",
        "- Final causal-order diagnostics use `sim.rold`/`sim.xold`, not `rnew`/`xnew`.",
        "",
        "## Answers",
        "",
        f"1. Best gamma by `final_energy`: `{float(best_energy['gamma']):.2f}`.",
        f"2. Best gamma by `interval_rmse`: `{float(best_rmse['gamma']):.2f}`.",
        f"3. Best gamma by `causal_f1`: `{float(best_f1['gamma']):.2f}`.",
        f"4. `causal_f1` aligns with `{f1_aligns}` in this three-point probe.",
        f"5. Dominant causal-order failure mode by total count: `{dominant}` (`missing={total_missing}`, `extra={total_extra}`).",
        f"6. Any `exact_match` true? `{'yes' if any_exact else 'no'}`. Any `success_flag` true? `{'yes' if any_success else 'no'}`.",
        "7. Conservative conclusion: accessibility/recoverability diagnostic only; no embeddability claim.",
        "",
        "## Key Table",
        "",
        "| gamma | final_energy | interval_rmse | causal_f1 | precision | recall | missing | extra | exact | success |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :---: | :---: |",
    ]
    for row in rows:
        lines.append(
            "| {gamma:.2f} | {energy:.6f} | {rmse:.6f} | {f1:.6f} | {precision:.6f} | {recall:.6f} | {missing} | {extra} | {exact} | {success} |".format(
                gamma=float(row["gamma"]),
                energy=float(row["final_energy"]),
                rmse=float(row["interval_rmse"]),
                f1=float(row["causal_f1"]),
                precision=float(row["causal_precision"]),
                recall=float(row["causal_recall"]),
                missing=int(row["missing_relations_count"]),
                extra=int(row["extra_relations_count"]),
                exact="yes" if bool(row["exact_match"]) else "no",
                success="yes" if bool(row["success_flag"]) else "no",
            )
        )
    lines += [
        "",
        "## Guardrails",
        "",
        "- Exploration only.",
        "- One N=24 known-truth case.",
        "- No physical gamma claim.",
        "- No embeddability claim.",
        "- No recovery claim unless exact-match/success criteria pass.",
    ]
    MD_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    generated_at_utc = datetime.now(timezone.utc).isoformat(timespec="seconds")
    start = time.perf_counter()
    case = _make_case()
    rows = [_run_one(case, gamma) for gamma in GAMMAS]
    total_runtime = time.perf_counter() - start
    write_csv(rows)
    write_svg(rows)
    write_markdown(rows, generated_at_utc, total_runtime)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
