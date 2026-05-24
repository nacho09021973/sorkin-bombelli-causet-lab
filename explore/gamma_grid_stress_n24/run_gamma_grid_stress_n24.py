#!/usr/bin/env python3
"""Exploratory N=24 gamma-grid stress test for SORKIN-2.

This script tests whether a gamma preference is stable under changes of
grid construction. It is exploratory only and does not run N=32 or N=36.
"""

from __future__ import annotations

import contextlib
import csv
import io
import math
import random
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import validation_suite as vs  # noqa: E402


OUT_DIR = Path(__file__).resolve().parent
CSV_PATH = OUT_DIR / "gamma_grid_stress_n24.csv"
MD_PATH = OUT_DIR / "gamma_grid_stress_n24.md"
SVG_PATH = OUT_DIR / "gamma_grid_stress_n24_summary.svg"
COMMAND = "python3 explore/gamma_grid_stress_n24/run_gamma_grid_stress_n24.py"

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
JITTER_SEED = 24051959
STOP_AFTER_FAMILY_SECONDS = 120.0

CSV_HEADERS = (
    "grid_family",
    "grid_index",
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
    "runtime_seconds",
)


def _linspace(start: float, stop: float, count: int) -> list[float]:
    if count == 1:
        return [start]
    return [start + (stop - start) * index / (count - 1) for index in range(count)]


def _logspace(start: float, stop: float, count: int) -> list[float]:
    if count == 1:
        return [start]
    ratio = stop / start
    return [start * ratio ** (index / (count - 1)) for index in range(count)]


def _jittered_logspace(start: float, stop: float, count: int, seed: int) -> list[float]:
    base = _logspace(start, stop, count)
    rng = random.Random(seed)
    values = [base[0]]
    for index in range(1, count - 1):
        left = base[index - 1]
        right = base[index + 1]
        step = min(base[index] - left, right - base[index])
        values.append(base[index] + rng.uniform(-0.35 * step, 0.35 * step))
    values.append(base[-1])
    return sorted(values)


def _custom_dense_grid() -> list[float]:
    coarse = _linspace(0.5, 0.74, 7)
    dense = _linspace(0.75, 0.9, 17)
    high = [0.925, 0.95]
    return sorted(set(round(value, 12) for value in [*coarse, *dense, *high]))


def _grid_families() -> list[tuple[str, list[float]]]:
    return [
        ("log_0p5_0p9_21", _logspace(0.5, 0.9, 21)),
        ("lin_0p5_0p9_21", _linspace(0.5, 0.9, 21)),
        ("log_0p45_0p95_21", _logspace(0.45, 0.95, 21)),
        ("lin_0p45_0p95_21", _linspace(0.45, 0.95, 21)),
        ("log_0p5_0p9_31", _logspace(0.5, 0.9, 31)),
        ("lin_0p5_0p9_31", _linspace(0.5, 0.9, 31)),
        ("jitter_log_0p5_0p9_21_seed24051959", _jittered_logspace(0.5, 0.9, 21, JITTER_SEED)),
        ("custom_dense_0p75_0p9", _custom_dense_grid()),
    ]


def _fmt(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if math.isfinite(value):
            return f"{value:.10g}"
        return "NA"
    return str(value)


def _case() -> vs.SprinkleCase:
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


def _run_one(case: vs.SprinkleCase, grid_family: str, grid_index: int, gamma: float) -> dict[str, object]:
    start = time.perf_counter()
    with contextlib.redirect_stdout(io.StringIO()):
        result = vs.run_recovery(
            case,
            optimizer_seed=OPTIMIZER_SEED,
            target_dim=TARGET_DIM,
            warmup_limit=WARMUP_LIMIT,
            anneal_limit=ANNEAL_LIMIT,
            max_data=MAX_DATA,
            initial_temp=INITIAL_TEMP,
            cooling_factor=gamma,
            backend=BACKEND,
        )
    runtime = time.perf_counter() - start
    energy_gap = result.final_energy - result.truth_energy
    success = (
        math.isfinite(energy_gap)
        and math.isfinite(result.interval_rmse)
        and energy_gap <= SUCCESS_GAP_THRESHOLD
    )
    return {
        "grid_family": grid_family,
        "grid_index": grid_index,
        "gamma": gamma,
        "n": N,
        "initial_temp": INITIAL_TEMP,
        "warmup_limit": WARMUP_LIMIT,
        "anneal_limit": ANNEAL_LIMIT,
        "max_data": MAX_DATA,
        "final_energy": result.final_energy,
        "energy_gap": energy_gap,
        "interval_rmse": result.interval_rmse,
        "mm_dim_truth": result.mm_dim_truth,
        "mm_dim_recovered": result.mm_dim_recovered,
        "success_flag": success,
        "runtime_seconds": runtime,
    }


def write_csv(rows: list[dict[str, object]]) -> None:
    with CSV_PATH.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(CSV_HEADERS)
        for row in rows:
            writer.writerow([_fmt(row[header]) for header in CSV_HEADERS])


def _best_by_family(rows: list[dict[str, object]], metric: str) -> list[dict[str, object]]:
    best_rows = []
    families = []
    for row in rows:
        family = str(row["grid_family"])
        if family not in families:
            families.append(family)
    for family in families:
        family_rows = [row for row in rows if row["grid_family"] == family]
        best_rows.append(min(family_rows, key=lambda row: float(row[metric])))
    return best_rows


def write_svg(rows: list[dict[str, object]]) -> None:
    width = 1180
    height = 760
    margin_left = 190
    margin_right = 42
    margin_top = 54
    margin_bottom = 64
    row_gap = 78
    plot_width = width - margin_left - margin_right
    families = []
    for row in rows:
        family = str(row["grid_family"])
        if family not in families:
            families.append(family)
    best_energy = {(row["grid_family"], row["grid_index"]) for row in _best_by_family(rows, "final_energy")}
    best_rmse = {(row["grid_family"], row["grid_index"]) for row in _best_by_family(rows, "interval_rmse")}

    def sx(gamma: float) -> float:
        return margin_left + (gamma - 0.45) * plot_width / (0.95 - 0.45)

    items = []
    axis = []
    for family_index, family in enumerate(families):
        y = margin_top + family_index * row_gap
        axis.append(
            f"<text x='18' y='{y + 5:.1f}' font-family='monospace' font-size='12' fill='#222'>{family}</text>"
        )
        axis.append(
            f"<line x1='{margin_left}' y1='{y}' x2='{width - margin_right}' y2='{y}' stroke='#bbb' stroke-width='1'/>"
        )
        family_rows = [row for row in rows if row["grid_family"] == family]
        for row in family_rows:
            key = (row["grid_family"], row["grid_index"])
            if key in best_energy and key in best_rmse:
                fill = "#7b3294"
                radius = 6.5
            elif key in best_energy:
                fill = "#d1495b"
                radius = 5.5
            elif key in best_rmse:
                fill = "#1f77b4"
                radius = 5.5
            else:
                fill = "#9c9c9c"
                radius = 3.0
            items.append(
                f"<circle cx='{sx(float(row['gamma'])):.2f}' cy='{y:.2f}' r='{radius}' fill='{fill}' opacity='0.88'/>"
            )

    ticks = []
    for tick in [0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]:
        x = sx(tick)
        ticks.append(f"<line x1='{x:.2f}' y1='{margin_top - 16}' x2='{x:.2f}' y2='{height - margin_bottom + 8}' stroke='#e2dfd8' stroke-width='1'/>")
        ticks.append(f"<text x='{x:.2f}' y='{height - 28}' text-anchor='middle' font-family='monospace' font-size='12'>{tick:.2f}</text>")

    svg = f"""<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>
  <rect width='100%' height='100%' fill='#f7f4ed'/>
  <text x='{width / 2:.0f}' y='28' text-anchor='middle' font-family='serif' font-size='22' fill='#222'>N=24 gamma-grid stress test</text>
  {''.join(ticks)}
  {''.join(axis)}
  {''.join(items)}
  <text x='{margin_left}' y='{height - 8}' font-family='monospace' font-size='12' fill='#444'>red = best final energy; blue = best interval RMSE; purple = both</text>
</svg>
"""
    SVG_PATH.write_text(svg, encoding="utf-8")


def write_markdown(
    rows: list[dict[str, object]],
    *,
    generated_at_utc: str,
    total_runtime_seconds: float,
    stopped_early: bool,
    stop_reason: str,
) -> None:
    best_energy = _best_by_family(rows, "final_energy")
    best_rmse = _best_by_family(rows, "interval_rmse")
    rmse_by_family = {row["grid_family"]: row for row in best_rmse}
    any_success = any(bool(row["success_flag"]) for row in rows)
    agreements = [
        float(row["gamma"]) == float(rmse_by_family[row["grid_family"]]["gamma"])
        for row in best_energy
    ]
    energy_gammas = [float(row["gamma"]) for row in best_energy]
    rmse_gammas = [float(row["gamma"]) for row in best_rmse]

    lines = [
        "# N=24 Gamma-Grid Stress Test",
        "",
        "**Status:** exploratory only; not confirmation.",
        "",
        "This stress test asks whether the apparent gamma preference in the N=24 probe is stable under changes of gamma grid, or whether it follows arbitrary grid nodes. It is an accessibility diagnostic for the historical Bombelli annealer only.",
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
        f"- completed grid points: `{len(rows)}`",
        f"- total runtime seconds: `{total_runtime_seconds:.3f}`",
        f"- stopped early: `{'true' if stopped_early else 'false'}`",
        f"- stop reason: `{stop_reason}`",
        "",
        "## Per-Family Winners",
        "",
        "| grid_family | points | best final_energy gamma | final_energy | best interval_rmse gamma | interval_rmse | agree | successes |",
        "| --- | ---: | ---: | ---: | ---: | ---: | :---: | ---: |",
    ]
    for energy_row in best_energy:
        family = str(energy_row["grid_family"])
        rmse_row = rmse_by_family[family]
        family_rows = [row for row in rows if row["grid_family"] == family]
        successes = sum(1 for row in family_rows if bool(row["success_flag"]))
        agree = float(energy_row["gamma"]) == float(rmse_row["gamma"])
        lines.append(
            "| {family} | {points} | {energy_gamma:.9f} | {energy:.6f} | {rmse_gamma:.9f} | {rmse:.6f} | {agree} | {successes} |".format(
                family=family,
                points=len(family_rows),
                energy_gamma=float(energy_row["gamma"]),
                energy=float(energy_row["final_energy"]),
                rmse_gamma=float(rmse_row["gamma"]),
                rmse=float(rmse_row["interval_rmse"]),
                agree="yes" if agree else "no",
                successes=successes,
            )
        )

    lines += [
        "",
        "## Required Questions",
        "",
        "1. Which gamma wins by `final_energy` for each grid family? See the per-family table above.",
        "2. Which gamma wins by `interval_rmse` for each grid family? See the per-family table above.",
        "3. Do winners cluster in a stable gamma region or follow arbitrary grid nodes?",
    ]
    lines.append(
        f"   `final_energy` winners span `{min(energy_gammas):.9f}` to `{max(energy_gammas):.9f}`; `interval_rmse` winners span `{min(rmse_gammas):.9f}` to `{max(rmse_gammas):.9f}`. The winners are grid- and metric-dependent in this exploratory run."
    )
    lines += [
        "4. Do `final_energy` and `interval_rmse` agree or disagree?",
        f"   They agree in `{sum(agreements)}` of `{len(agreements)}` completed grid families.",
        "5. Are any `success_flag` values true?",
        f"   `{'yes' if any_success else 'no'}`.",
        "6. Is there any evidence for a robust gamma optimum?",
        "   No robust optimum is established by this exploratory stress test.",
        "7. Warning: no physical resonance claim, no embeddability claim, and no recovery claim.",
        "",
        "## Interpretation Guardrails",
        "",
        "- Exploration only.",
        "- N=24 only.",
        "- Fixed T0=100.",
        "- Short Phase 2B budget only.",
        "- Accessibility diagnostic only.",
    ]
    MD_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_artifacts(
    rows: list[dict[str, object]],
    generated_at_utc: str,
    total_runtime_seconds: float,
    stopped_early: bool,
    stop_reason: str,
) -> None:
    write_csv(rows)
    write_svg(rows)
    write_markdown(
        rows,
        generated_at_utc=generated_at_utc,
        total_runtime_seconds=total_runtime_seconds,
        stopped_early=stopped_early,
        stop_reason=stop_reason,
    )


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    generated_at_utc = datetime.now(timezone.utc).isoformat(timespec="seconds")
    total_start = time.perf_counter()
    rows: list[dict[str, object]] = []
    case = _case()
    stopped_early = False
    stop_reason = "completed all requested grid families"

    for grid_family, gammas in _grid_families():
        family_start = time.perf_counter()
        for grid_index, gamma in enumerate(gammas):
            rows.append(_run_one(case, grid_family, grid_index, gamma))
        total_runtime_seconds = time.perf_counter() - total_start
        if time.perf_counter() - family_start > STOP_AFTER_FAMILY_SECONDS:
            stopped_early = True
            stop_reason = (
                f"stopped after completing {grid_family}; family runtime exceeded "
                f"{STOP_AFTER_FAMILY_SECONDS:.1f} seconds"
            )
            _write_artifacts(
                rows,
                generated_at_utc,
                total_runtime_seconds,
                stopped_early,
                stop_reason,
            )
            return 0
        _write_artifacts(
            rows,
            generated_at_utc,
            total_runtime_seconds,
            stopped_early,
            stop_reason,
        )

    total_runtime_seconds = time.perf_counter() - total_start
    _write_artifacts(
        rows,
        generated_at_utc,
        total_runtime_seconds,
        stopped_early,
        stop_reason,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
