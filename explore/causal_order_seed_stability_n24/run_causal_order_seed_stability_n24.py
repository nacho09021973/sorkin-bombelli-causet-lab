#!/usr/bin/env python3
"""Multi-seed causal-order stability probe for SORKIN-2.

Exploratory only: tests whether low gamma, final energy, and causal F1
alignment is stable across optimizer seeds for one N=24 known-truth case.
"""

from __future__ import annotations

import contextlib
import csv
import io
import math
import statistics
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
CSV_PATH = OUT_DIR / "causal_order_seed_stability_n24.csv"
MD_PATH = OUT_DIR / "causal_order_seed_stability_n24.md"
SVG_PATH = OUT_DIR / "causal_order_seed_stability_n24.svg"
COMMAND = "python3 explore/causal_order_seed_stability_n24/run_causal_order_seed_stability_n24.py"

D_SPACETIME = 2
N = 24
CASE_SEED = 1959
TARGET_DIM = D_SPACETIME - 1
INITIAL_TEMP = 100.0
WARMUP_LIMIT = 10
ANNEAL_LIMIT = 10
MAX_DATA = 4
BACKEND = "cpu"
SUCCESS_GAP_THRESHOLD = 1.0
GAMMAS = (0.50, 0.53, 0.55, 0.57, 0.60)
OPTIMIZER_SEEDS = (1959, 1962, 1987, 2001, 2026)

CSV_HEADERS = (
    "gamma",
    "optimizer_seed",
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


def _run_one(case: vs.SprinkleCase, gamma: float, optimizer_seed: int) -> dict[str, object]:
    start = time.perf_counter()
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
        "optimizer_seed": optimizer_seed,
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


def _by_gamma(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped = []
    for gamma in GAMMAS:
        g_rows = [row for row in rows if float(row["gamma"]) == gamma]
        f1s = [float(row["causal_f1"]) for row in g_rows]
        energies = [float(row["final_energy"]) for row in g_rows]
        rmses = [float(row["interval_rmse"]) for row in g_rows]
        grouped.append(
            {
                "gamma": gamma,
                "runs": len(g_rows),
                "mean_causal_f1": statistics.fmean(f1s),
                "median_causal_f1": statistics.median(f1s),
                "min_causal_f1": min(f1s),
                "max_causal_f1": max(f1s),
                "mean_final_energy": statistics.fmean(energies),
                "mean_interval_rmse": statistics.fmean(rmses),
                "missing_total": sum(int(row["missing_relations_count"]) for row in g_rows),
                "extra_total": sum(int(row["extra_relations_count"]) for row in g_rows),
                "exact_matches": sum(1 for row in g_rows if bool(row["exact_match"])),
                "successes": sum(1 for row in g_rows if bool(row["success_flag"])),
            }
        )
    return grouped


def _pearson(xs: list[float], ys: list[float]) -> float:
    if len(xs) != len(ys) or len(xs) < 2:
        return float("nan")
    mean_x = statistics.fmean(xs)
    mean_y = statistics.fmean(ys)
    dx = [x - mean_x for x in xs]
    dy = [y - mean_y for y in ys]
    denom = math.sqrt(sum(x * x for x in dx) * sum(y * y for y in dy))
    if denom == 0.0:
        return float("nan")
    return sum(x * y for x, y in zip(dx, dy)) / denom


def write_svg(rows: list[dict[str, object]], aggregates: list[dict[str, object]]) -> None:
    width = 1040
    height = 620
    margin = 68
    gap = 54
    panel_width = (width - 2 * margin - gap) / 2
    panel_height = height - 2 * margin
    gammas = [float(row["gamma"]) for row in rows]
    f1s = [float(row["causal_f1"]) for row in rows]
    energies = [float(row["final_energy"]) for row in rows]

    def sx_gamma(gamma: float) -> float:
        return margin + (gamma - min(gammas)) * panel_width / (max(gammas) - min(gammas))

    def sy(value: float, values: list[float]) -> float:
        low = min(values)
        high = max(values)
        if low == high:
            low -= 1.0
            high += 1.0
        return height - margin - (value - low) * panel_height / (high - low)

    right_x = margin + panel_width + gap

    chunks = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "  <rect width='100%' height='100%' fill='#f7f4ed'/>",
        f"  <text x='{width / 2:.0f}' y='30' text-anchor='middle' font-family='serif' font-size='22' fill='#222'>N=24 causal-order seed stability</text>",
        f"  <text x='{margin + panel_width / 2:.1f}' y='58' text-anchor='middle' font-family='monospace' font-size='13'>causal F1 by gamma</text>",
        f"  <text x='{right_x + panel_width / 2:.1f}' y='58' text-anchor='middle' font-family='monospace' font-size='13'>final energy vs causal F1</text>",
        f"  <line x1='{margin}' y1='{height - margin}' x2='{margin + panel_width}' y2='{height - margin}' stroke='#333' stroke-width='2'/>",
        f"  <line x1='{margin}' y1='{margin}' x2='{margin}' y2='{height - margin}' stroke='#333' stroke-width='2'/>",
        f"  <line x1='{right_x}' y1='{height - margin}' x2='{right_x + panel_width}' y2='{height - margin}' stroke='#333' stroke-width='2'/>",
        f"  <line x1='{right_x}' y1='{margin}' x2='{right_x}' y2='{height - margin}' stroke='#333' stroke-width='2'/>",
    ]
    for row in rows:
        gamma = float(row["gamma"])
        seed = int(row["optimizer_seed"])
        jitter = ((seed % 11) - 5) * 0.0012
        chunks.append(
            f"  <circle cx='{sx_gamma(gamma + jitter):.2f}' cy='{sy(float(row['causal_f1']), f1s):.2f}' r='4.5' fill='#2a9d8f' opacity='0.78'/>"
        )
        chunks.append(
            f"  <circle cx='{right_x + (float(row['final_energy']) - min(energies)) * panel_width / (max(energies) - min(energies)):.2f}' cy='{sy(float(row['causal_f1']), f1s):.2f}' r='4.5' fill='#1f77b4' opacity='0.78'/>"
        )
    mean_points = " ".join(
        f"{sx_gamma(float(row['gamma'])):.2f},{sy(float(row['mean_causal_f1']), f1s):.2f}"
        for row in aggregates
    )
    chunks.append(f"  <polyline fill='none' stroke='#d1495b' stroke-width='2.5' points='{mean_points}'/>")
    for gamma in GAMMAS:
        chunks.append(
            f"  <text x='{sx_gamma(gamma):.2f}' y='{height - margin + 22}' text-anchor='middle' font-family='monospace' font-size='10'>{gamma:.2f}</text>"
        )
    chunks.append(
        f"  <text x='{right_x}' y='{height - margin + 22}' font-family='monospace' font-size='10'>final energy low -> high</text>"
    )
    chunks.append("  <text x='68' y='600' font-family='monospace' font-size='12' fill='#444'>red line = mean causal F1 by gamma; points = optimizer seeds</text>")
    chunks.append("</svg>")
    SVG_PATH.write_text("\n".join(chunks) + "\n", encoding="utf-8")


def write_markdown(
    rows: list[dict[str, object]],
    aggregates: list[dict[str, object]],
    generated_at_utc: str,
    total_runtime: float,
) -> None:
    best_mean = max(aggregates, key=lambda row: float(row["mean_causal_f1"]))
    best_median = max(aggregates, key=lambda row: float(row["median_causal_f1"]))
    best_worst = max(aggregates, key=lambda row: float(row["min_causal_f1"]))
    corr = _pearson(
        [float(row["final_energy"]) for row in rows],
        [float(row["causal_f1"]) for row in rows],
    )
    total_missing = sum(int(row["missing_relations_count"]) for row in rows)
    total_extra = sum(int(row["extra_relations_count"]) for row in rows)
    any_exact = any(bool(row["exact_match"]) for row in rows)
    any_success = any(bool(row["success_flag"]) for row in rows)

    lines = [
        "# N=24 Causal-Order Seed Stability Probe",
        "",
        "**Status:** exploratory only; not confirmation.",
        "",
        "This probe tests whether the apparent alignment between low gamma, final energy, and causal F1 is stable across optimizer seeds. It is an accessibility/recoverability diagnostic only.",
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
        f"- gammas: `{', '.join(f'{g:.2f}' for g in GAMMAS)}`",
        f"- optimizer seeds: `{', '.join(str(s) for s in OPTIMIZER_SEEDS)}`",
        f"- runs completed: `{len(rows)}`",
        f"- total runtime seconds: `{total_runtime:.3f}`",
        "- Final causal-order diagnostics use `sim.rold`/`sim.xold`, not `rnew`/`xnew`.",
        "",
        "## Answers",
        "",
        f"1. Best mean `causal_f1`: gamma `{float(best_mean['gamma']):.2f}`.",
        f"2. Best median `causal_f1`: gamma `{float(best_median['gamma']):.2f}`.",
        f"3. Best worst-case `causal_f1`: gamma `{float(best_worst['gamma']):.2f}`.",
        f"4. Pearson correlation between `final_energy` and `causal_f1` across seeds: `{corr:.6f}`.",
        f"5. Missing relations dominate? `{'yes' if total_missing > total_extra else 'no'}` (`missing={total_missing}`, `extra={total_extra}`).",
        f"6. Any `exact_match` true? `{'yes' if any_exact else 'no'}`. Any `success_flag` true? `{'yes' if any_success else 'no'}`.",
        "7. Conservative conclusion only: no recovery claim, no embeddability claim.",
        "",
        "## Aggregate Table",
        "",
        "| gamma | runs | mean F1 | median F1 | worst F1 | best F1 | mean final E | mean RMSE | missing | extra | exact | success |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in aggregates:
        lines.append(
            "| {gamma:.2f} | {runs} | {mean:.6f} | {median:.6f} | {worst:.6f} | {best:.6f} | {energy:.6f} | {rmse:.6f} | {missing} | {extra} | {exact} | {success} |".format(
                gamma=float(row["gamma"]),
                runs=int(row["runs"]),
                mean=float(row["mean_causal_f1"]),
                median=float(row["median_causal_f1"]),
                worst=float(row["min_causal_f1"]),
                best=float(row["max_causal_f1"]),
                energy=float(row["mean_final_energy"]),
                rmse=float(row["mean_interval_rmse"]),
                missing=int(row["missing_total"]),
                extra=int(row["extra_total"]),
                exact=int(row["exact_matches"]),
                success=int(row["successes"]),
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
        "- No recovery claim.",
    ]
    MD_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    generated_at_utc = datetime.now(timezone.utc).isoformat(timespec="seconds")
    start = time.perf_counter()
    case = _make_case()
    rows = [
        _run_one(case, gamma, optimizer_seed)
        for gamma in GAMMAS
        for optimizer_seed in OPTIMIZER_SEEDS
    ]
    total_runtime = time.perf_counter() - start
    aggregates = _by_gamma(rows)
    write_csv(rows)
    write_svg(rows, aggregates)
    write_markdown(rows, aggregates, generated_at_utc, total_runtime)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
