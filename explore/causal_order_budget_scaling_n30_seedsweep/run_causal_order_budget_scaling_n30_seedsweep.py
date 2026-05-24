#!/usr/bin/env python3
"""N=30 optimizer-seed budget/misalignment probe for SORKIN-2.

Exploratory only: fixed family/case/schedule, sweep the canonical
optimizer-seed set at N=30 with short and medium budgets to test
whether the short->medium energy/causal-metric misalignment seen at
optimizer seed 1987 persists across seeds or was seed-specific.

The optimizer seed 1987 short/medium rows already exist in
``explore/n_ladder_budget_misalignment_probe/`` and are reused here.
"""

from __future__ import annotations

import contextlib
import csv
import io
import math
import multiprocessing as mp
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import causet_invariants  # noqa: E402, F401
import cones  # noqa: E402
import validation_suite as vs  # noqa: E402


OUT_DIR = Path(__file__).resolve().parent
CSV_PATH = OUT_DIR / "causal_order_budget_scaling_n30_seedsweep.csv"
MD_PATH = OUT_DIR / "causal_order_budget_scaling_n30_seedsweep.md"
SVG_PATH = OUT_DIR / "causal_order_budget_scaling_n30_seedsweep.svg"
LOG_PATH = OUT_DIR / "causal_order_budget_scaling_n30_seedsweep.log"
COMMAND = (
    "python3 explore/causal_order_budget_scaling_n30_seedsweep/"
    "run_causal_order_budget_scaling_n30_seedsweep.py"
)

FAMILY = "minkowski"
D_SPACETIME = 2
N = 30
CASE_SEED = 1959
OPTIMIZER_SEEDS = (1959, 1962, 1987, 2001, 2026)
TARGET_DIM = D_SPACETIME - 1
INITIAL_TEMP = 100.0
GAMMA = 0.50
H = 1
BACKEND = "cpu"
SUCCESS_GAP_THRESHOLD = 1.0
PER_RUN_TIMEOUT_SECONDS = 900.0

BUDGETS = (
    {"budget_label": "short_10_10_4", "warmup_limit": 10, "anneal_limit": 10, "max_data": 4},
    {"budget_label": "medium_25_25_8", "warmup_limit": 25, "anneal_limit": 25, "max_data": 8},
)

REUSE_CSV = (
    ROOT / "explore" / "n_ladder_budget_misalignment_probe" /
    "n_ladder_budget_misalignment_probe.csv"
)
REUSE_OPTIMIZER_SEED = 1987

CSV_HEADERS = (
    "n",
    "family",
    "d_spacetime",
    "case_seed",
    "optimizer_seed",
    "t0",
    "gamma",
    "h",
    "budget_label",
    "row_status",
    "final_energy",
    "interval_rmse",
    "causal_precision",
    "causal_recall",
    "causal_f1",
    "missing_relations_count",
    "extra_relations_count",
    "exact_match",
    "success_flag",
    "runtime_seconds",
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


def _safe_div(numerator: float, denominator: float) -> float:
    if denominator == 0.0:
        return 0.0
    return numerator / denominator


def _parse_bool(value: object) -> bool:
    return str(value).strip().lower() == "true"


def _parse_float(value: str) -> float:
    if value == "NA":
        return float("nan")
    return float(value)


def _log(message: str) -> None:
    stamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
    with LOG_PATH.open("a", encoding="utf-8") as fh:
        fh.write(f"{stamp} {message}\n")


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


def _blank_row(
    optimizer_seed: int,
    budget: dict[str, object],
    row_status: str,
    runtime: float,
    note: str,
) -> dict[str, object]:
    return {
        "n": N,
        "family": FAMILY,
        "d_spacetime": D_SPACETIME,
        "case_seed": CASE_SEED,
        "optimizer_seed": optimizer_seed,
        "t0": INITIAL_TEMP,
        "gamma": GAMMA,
        "h": H,
        "budget_label": str(budget["budget_label"]),
        "row_status": row_status,
        "final_energy": float("nan"),
        "interval_rmse": float("nan"),
        "causal_precision": float("nan"),
        "causal_recall": float("nan"),
        "causal_f1": float("nan"),
        "missing_relations_count": 0,
        "extra_relations_count": 0,
        "exact_match": False,
        "success_flag": False,
        "runtime_seconds": runtime,
        "note": note,
    }


def _run_one(
    case: vs.SprinkleCase,
    budget: dict[str, object],
    optimizer_seed: int,
) -> dict[str, object]:
    start = time.perf_counter()
    with tempfile.TemporaryDirectory() as tmpdir, contextlib.redirect_stdout(io.StringIO()):
        sim = cones.ConesSimulator(
            z=case.matrix,
            dim=TARGET_DIM,
            seed=optimizer_seed,
            interactive=False,
            max_data=int(budget["max_data"]),
            plot_path=None,
            warmup_limit=int(budget["warmup_limit"]),
            anneal_limit=int(budget["anneal_limit"]),
            initial_temp=INITIAL_TEMP,
            cooling_factor=GAMMA,
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

    precision = _safe_div(correct_relations, comparison.total_relations_induced)
    recall = _safe_div(correct_relations, comparison.total_relations_target)
    f1 = _safe_div(2.0 * precision * recall, precision + recall)

    return {
        "n": N,
        "family": FAMILY,
        "d_spacetime": D_SPACETIME,
        "case_seed": CASE_SEED,
        "optimizer_seed": optimizer_seed,
        "t0": INITIAL_TEMP,
        "gamma": GAMMA,
        "h": H,
        "budget_label": str(budget["budget_label"]),
        "row_status": "completed",
        "final_energy": final_energy,
        "interval_rmse": interval_residual,
        "causal_precision": precision,
        "causal_recall": recall,
        "causal_f1": f1,
        "missing_relations_count": missing_count,
        "extra_relations_count": extra_count,
        "exact_match": comparison.exact_match,
        "success_flag": success,
        "runtime_seconds": runtime,
        "note": "",
    }


def _run_one_worker(
    queue: mp.Queue,
    case: vs.SprinkleCase,
    budget: dict[str, object],
    optimizer_seed: int,
) -> None:
    try:
        queue.put(("ok", _run_one(case, budget, optimizer_seed)))
    except BaseException as exc:  # pragma: no cover - defensive child process path
        queue.put(("error", repr(exc)))


def _run_one_with_timeout(
    case: vs.SprinkleCase,
    budget: dict[str, object],
    optimizer_seed: int,
) -> tuple[dict[str, object] | None, str | None]:
    queue: mp.Queue = mp.Queue()
    process = mp.Process(
        target=_run_one_worker,
        args=(queue, case, budget, optimizer_seed),
    )
    start = time.perf_counter()
    process.start()
    process.join(PER_RUN_TIMEOUT_SECONDS)
    runtime = time.perf_counter() - start
    if process.is_alive():
        process.terminate()
        process.join(5.0)
        if process.is_alive():
            process.kill()
            process.join()
        return _blank_row(
            optimizer_seed,
            budget,
            "timeout",
            runtime,
            f"timeout after {PER_RUN_TIMEOUT_SECONDS:.1f}s",
        ), None
    if queue.empty():
        return None, (
            f"child exited without result at "
            f"{budget['budget_label']} seed {optimizer_seed}"
        )
    status, payload = queue.get()
    if status == "ok":
        return payload, None
    return None, str(payload)


def append_csv_row(row: dict[str, object]) -> None:
    needs_header = not CSV_PATH.exists() or CSV_PATH.stat().st_size == 0
    with CSV_PATH.open("a", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        if needs_header:
            writer.writerow(CSV_HEADERS)
        writer.writerow([_fmt(row[header]) for header in CSV_HEADERS])


def read_existing_rows() -> list[dict[str, object]]:
    if not CSV_PATH.exists():
        return []
    rows: list[dict[str, object]] = []
    with CSV_PATH.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            parsed: dict[str, object] = {}
            for header in CSV_HEADERS:
                value = row[header]
                if header in {"budget_label", "row_status", "note", "family"}:
                    parsed[header] = value
                elif header in {"success_flag", "exact_match"}:
                    parsed[header] = _parse_bool(value)
                elif header in {
                    "n",
                    "d_spacetime",
                    "case_seed",
                    "optimizer_seed",
                    "h",
                    "missing_relations_count",
                    "extra_relations_count",
                }:
                    parsed[header] = int(value)
                else:
                    parsed[header] = _parse_float(value)
            rows.append(parsed)
    return rows


def _load_reuse_rows() -> list[dict[str, object]]:
    """Pull seed=1987 N=30 short/medium rows from the n_ladder probe CSV."""
    if not REUSE_CSV.exists():
        return []
    reuse_rows: list[dict[str, object]] = []
    with REUSE_CSV.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            try:
                row_n = int(row.get("n", "0"))
                row_seed = int(row.get("optimizer_seed", "0"))
            except ValueError:
                continue
            if row_n != N or row_seed != REUSE_OPTIMIZER_SEED:
                continue
            status = row.get("row_status", "")
            if status != "completed":
                continue
            label = row.get("budget_label", "")
            if label not in {b["budget_label"] for b in BUDGETS}:
                continue
            parsed: dict[str, object] = {
                "n": N,
                "family": FAMILY,
                "d_spacetime": D_SPACETIME,
                "case_seed": CASE_SEED,
                "optimizer_seed": REUSE_OPTIMIZER_SEED,
                "t0": INITIAL_TEMP,
                "gamma": GAMMA,
                "h": H,
                "budget_label": label,
                "row_status": "completed",
                "final_energy": _parse_float(row.get("final_energy", "NA")),
                "interval_rmse": _parse_float(row.get("interval_rmse", "NA")),
                "causal_precision": _parse_float(row.get("causal_precision", "NA")),
                "causal_recall": _parse_float(row.get("causal_recall", "NA")),
                "causal_f1": _parse_float(row.get("causal_f1", "NA")),
                "missing_relations_count": int(row.get("missing_relations_count", "0")),
                "extra_relations_count": int(row.get("extra_relations_count", "0")),
                "exact_match": _parse_bool(row.get("exact_match", "false")),
                "success_flag": _parse_bool(row.get("success_flag", "false")),
                "runtime_seconds": _parse_float(row.get("runtime_seconds", "NA")),
                "note": f"reused from {REUSE_CSV.relative_to(ROOT)}",
            }
            reuse_rows.append(parsed)
    return reuse_rows


def _interpretation_label(
    e_short: float,
    e_medium: float,
    f1_short: float,
    f1_medium: float,
    *,
    e_eps: float = 1.0,
    f1_eps: float = 0.02,
) -> str:
    if not all(map(math.isfinite, (e_short, e_medium, f1_short, f1_medium))):
        return "no_clear_change"
    delta_e = e_medium - e_short
    delta_f1 = f1_medium - f1_short
    if delta_e < -e_eps and delta_f1 > f1_eps:
        return "aligned_improvement"
    if delta_e < -e_eps and delta_f1 < -f1_eps:
        return "energy_causality_misalignment"
    if delta_e > e_eps and delta_f1 > f1_eps:
        return "causal_improves_without_energy"
    return "no_clear_change"


def _seed_summary(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    short_by_seed: dict[int, dict[str, object]] = {}
    medium_by_seed: dict[int, dict[str, object]] = {}
    for row in rows:
        if str(row["row_status"]) != "completed":
            continue
        seed = int(row["optimizer_seed"])
        label = str(row["budget_label"])
        if label == "short_10_10_4":
            short_by_seed[seed] = row
        elif label == "medium_25_25_8":
            medium_by_seed[seed] = row

    out: list[dict[str, object]] = []
    for seed in OPTIMIZER_SEEDS:
        short = short_by_seed.get(seed)
        medium = medium_by_seed.get(seed)
        if short is None or medium is None:
            continue
        e_short = float(short["final_energy"])
        e_medium = float(medium["final_energy"])
        f1_short = float(short["causal_f1"])
        f1_medium = float(medium["causal_f1"])
        out.append(
            {
                "optimizer_seed": seed,
                "E_short": e_short,
                "E_medium": e_medium,
                "delta_E": e_medium - e_short,
                "F1_short": f1_short,
                "F1_medium": f1_medium,
                "delta_F1": f1_medium - f1_short,
                "recall_short": float(short["causal_recall"]),
                "recall_medium": float(medium["causal_recall"]),
                "missing_short": int(short["missing_relations_count"]),
                "missing_medium": int(medium["missing_relations_count"]),
                "label": _interpretation_label(
                    e_short, e_medium, f1_short, f1_medium
                ),
            }
        )
    return out


def _decision_outcome(summary: list[dict[str, object]]) -> tuple[str, int, int]:
    total = len(summary)
    misalign = sum(
        1 for row in summary if row["label"] == "energy_causality_misalignment"
    )
    if total == 0:
        return ("no_seeds_completed", misalign, total)
    if misalign >= 4:
        return ("misalignment_persists_across_optimizer_seeds", misalign, total)
    if misalign <= 2:
        return ("misalignment_appears_seed_specific", misalign, total)
    return ("mixed_inconclusive", misalign, total)


def _fmt_f(value: float, digits: int = 6) -> str:
    if not math.isfinite(value):
        return "NA"
    return f"{value:.{digits}f}"


def write_svg(summary: list[dict[str, object]]) -> None:
    width = 1120
    height = 640
    margin_left = 96
    margin_right = 96
    margin_top = 80
    margin_bottom = 100
    plot_w = width - margin_left - margin_right
    plot_h = height - margin_top - margin_bottom

    chunks = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' "
        f"viewBox='0 0 {width} {height}'>",
        "  <rect width='100%' height='100%' fill='#f7f4ed'/>",
        f"  <text x='{width / 2:.0f}' y='34' text-anchor='middle' "
        "font-family='serif' font-size='22' fill='#222'>N=30 optimizer-seed "
        "budget/misalignment probe: delta_F1 and delta_E by seed "
        "(short -> medium)</text>",
    ]

    if not summary:
        chunks.append(
            "  <text x='560' y='320' text-anchor='middle' "
            "font-family='monospace' font-size='16' fill='#444'>"
            "no completed seeds</text>"
        )
        chunks.append("</svg>")
        SVG_PATH.write_text("\n".join(chunks) + "\n", encoding="utf-8")
        return

    seeds = [int(row["optimizer_seed"]) for row in summary]
    delta_f1 = [float(row["delta_F1"]) for row in summary]
    delta_e = [float(row["delta_E"]) for row in summary]
    labels = [str(row["label"]) for row in summary]

    n_seeds = len(seeds)
    if n_seeds == 1:
        def x_of(i: int) -> float:
            return margin_left + plot_w / 2
    else:
        def x_of(i: int) -> float:
            return margin_left + i * plot_w / (n_seeds - 1)

    f1_min = min(delta_f1 + [0.0])
    f1_max = max(delta_f1 + [0.0])
    if f1_min == f1_max:
        f1_min -= 0.1
        f1_max += 0.1

    def y_of_f1(v: float) -> float:
        return margin_top + plot_h - (v - f1_min) * plot_h / (f1_max - f1_min)

    e_min = min(delta_e + [0.0])
    e_max = max(delta_e + [0.0])
    if e_min == e_max:
        e_min -= 1.0
        e_max += 1.0

    def y_of_e(v: float) -> float:
        return margin_top + plot_h - (v - e_min) * plot_h / (e_max - e_min)

    chunks.append(
        f"  <rect x='{margin_left}' y='{margin_top}' width='{plot_w}' "
        f"height='{plot_h}' fill='none' stroke='#333' stroke-width='2'/>"
    )

    y_zero_f1 = y_of_f1(0.0)
    chunks.append(
        f"  <line x1='{margin_left}' y1='{y_zero_f1:.2f}' "
        f"x2='{margin_left + plot_w}' y2='{y_zero_f1:.2f}' "
        "stroke='#999' stroke-width='1' stroke-dasharray='4 4'/>"
    )

    for i, seed in enumerate(seeds):
        x = x_of(i)
        chunks.append(
            f"  <line x1='{x:.2f}' y1='{margin_top + plot_h}' "
            f"x2='{x:.2f}' y2='{margin_top + plot_h + 6}' "
            "stroke='#333' stroke-width='1.5'/>"
        )
        chunks.append(
            f"  <text x='{x:.2f}' y='{margin_top + plot_h + 22}' "
            f"text-anchor='middle' font-family='monospace' font-size='12' "
            f"fill='#222'>seed {seed}</text>"
        )

    chunks.append(
        f"  <text x='{margin_left + plot_w / 2:.2f}' "
        f"y='{margin_top + plot_h + 58}' text-anchor='middle' "
        "font-family='monospace' font-size='13' fill='#222'>"
        "optimizer seed</text>"
    )

    for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
        v = f1_min + frac * (f1_max - f1_min)
        y = y_of_f1(v)
        chunks.append(
            f"  <line x1='{margin_left - 6}' y1='{y:.2f}' "
            f"x2='{margin_left}' y2='{y:.2f}' stroke='#333' stroke-width='1'/>"
        )
        chunks.append(
            f"  <text x='{margin_left - 10}' y='{y + 4:.2f}' text-anchor='end' "
            f"font-family='monospace' font-size='11' fill='#1f77b4'>"
            f"{v:.3f}</text>"
        )

    for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
        v = e_min + frac * (e_max - e_min)
        y = y_of_e(v)
        chunks.append(
            f"  <line x1='{margin_left + plot_w}' y1='{y:.2f}' "
            f"x2='{margin_left + plot_w + 6}' y2='{y:.2f}' "
            "stroke='#333' stroke-width='1'/>"
        )
        chunks.append(
            f"  <text x='{margin_left + plot_w + 10}' y='{y + 4:.2f}' "
            f"text-anchor='start' font-family='monospace' font-size='11' "
            f"fill='#d1495b'>{v:.2f}</text>"
        )

    chunks.append(
        f"  <text x='28' y='{margin_top + plot_h / 2:.2f}' text-anchor='middle' "
        "font-family='monospace' font-size='13' fill='#1f77b4' "
        f"transform='rotate(-90 28 {margin_top + plot_h / 2:.2f})'>"
        "delta_F1 = F1_medium - F1_short</text>"
    )
    chunks.append(
        f"  <text x='{width - 28}' y='{margin_top + plot_h / 2:.2f}' "
        "text-anchor='middle' font-family='monospace' font-size='13' "
        f"fill='#d1495b' transform='rotate(90 {width - 28} "
        f"{margin_top + plot_h / 2:.2f})'>delta_E = E_medium - E_short</text>"
    )

    bar_half = max(min(plot_w / (max(n_seeds, 1) * 4.0), 28.0), 6.0)
    for i, (vf1, ve, lbl) in enumerate(zip(delta_f1, delta_e, labels)):
        x = x_of(i)
        yf = y_of_f1(vf1)
        # delta_F1 bar from zero line.
        fill_f1 = "#1f77b4" if vf1 >= 0.0 else "#d35a85"
        chunks.append(
            f"  <rect x='{x - bar_half:.2f}' y='{min(yf, y_zero_f1):.2f}' "
            f"width='{2 * bar_half:.2f}' height='{abs(yf - y_zero_f1):.2f}' "
            f"fill='{fill_f1}' fill-opacity='0.55' stroke='#1f77b4' "
            "stroke-width='1'/>"
        )
        # delta_E marker (red ring).
        ye = y_of_e(ve)
        chunks.append(
            f"  <circle cx='{x:.2f}' cy='{ye:.2f}' r='7' fill='none' "
            "stroke='#d1495b' stroke-width='2.5'/>"
        )
        label_y = y_zero_f1 + (14 if yf >= y_zero_f1 else -8)
        chunks.append(
            f"  <text x='{x:.2f}' y='{label_y:.2f}' text-anchor='middle' "
            f"font-family='monospace' font-size='10' fill='#222'>{lbl}</text>"
        )

    # delta_E polyline through the rings.
    if n_seeds >= 2:
        e_points = " ".join(
            f"{x_of(i):.2f},{y_of_e(ve):.2f}" for i, ve in enumerate(delta_e)
        )
        chunks.append(
            f"  <polyline fill='none' stroke='#d1495b' stroke-width='2' "
            f"stroke-dasharray='6 5' points='{e_points}'/>"
        )

    # Legend.
    lg_x = margin_left + 16
    lg_y = margin_top + 16
    chunks.append(
        f"  <rect x='{lg_x - 8}' y='{lg_y - 14}' width='320' height='70' "
        "fill='#ffffff' fill-opacity='0.85' stroke='#aaa' stroke-width='1'/>"
    )
    chunks.append(
        f"  <rect x='{lg_x}' y='{lg_y - 8}' width='28' height='14' "
        "fill='#1f77b4' fill-opacity='0.55' stroke='#1f77b4'/>"
    )
    chunks.append(
        f"  <text x='{lg_x + 36}' y='{lg_y + 4}' font-family='monospace' "
        "font-size='12' fill='#222'>delta_F1 bar (left axis, blue)</text>"
    )
    chunks.append(
        f"  <circle cx='{lg_x + 14}' cy='{lg_y + 22}' r='7' fill='none' "
        "stroke='#d1495b' stroke-width='2.5'/>"
    )
    chunks.append(
        f"  <text x='{lg_x + 36}' y='{lg_y + 26}' font-family='monospace' "
        "font-size='12' fill='#222'>delta_E ring (right axis, red)</text>"
    )
    chunks.append(
        f"  <text x='{lg_x}' y='{lg_y + 48}' font-family='monospace' "
        "font-size='10' fill='#444'>label below zero line: "
        "interpretation per seed</text>"
    )

    chunks.append(
        f"  <text x='{margin_left}' y='{height - 18}' font-family='monospace' "
        "font-size='11' fill='#444'>Exploratory only; family minkowski; "
        f"d_spacetime {D_SPACETIME}; case seed {CASE_SEED}; N {N}; "
        f"T0 {INITIAL_TEMP:.0f}; gamma {GAMMA}; h {H}.</text>"
    )
    chunks.append("</svg>")
    SVG_PATH.write_text("\n".join(chunks) + "\n", encoding="utf-8")


def write_markdown(
    rows: list[dict[str, object]],
    summary: list[dict[str, object]],
    reused_seeds: list[int],
    newly_run_seeds: list[int],
    skipped: list[tuple[int, str, str]],
    generated_at_utc: str,
    total_runtime: float,
) -> None:
    outcome, misalign, total = _decision_outcome(summary)
    aligned = [int(row["optimizer_seed"]) for row in summary
               if row["label"] == "aligned_improvement"]
    misaligned = [int(row["optimizer_seed"]) for row in summary
                  if row["label"] == "energy_causality_misalignment"]
    causal_only = [int(row["optimizer_seed"]) for row in summary
                   if row["label"] == "causal_improves_without_energy"]
    unclear = [int(row["optimizer_seed"]) for row in summary
               if row["label"] == "no_clear_change"]

    lines = [
        "# N=30 Optimizer-Seed Budget/Misalignment Probe",
        "",
        "exploratory only; not confirmation.",
        "",
        "At N=30, fixed family/case/schedule/budgets, does the short -> medium "
        "energy/causal-metric misalignment observed at optimizer seed 1987 "
        "persist across optimizer seeds, or was it seed-specific?",
        "",
        "## Provenance",
        "",
        f"- Output directory: `{OUT_DIR.relative_to(ROOT)}`",
        f"- Runner: `{Path(__file__).relative_to(ROOT)}`",
        f"- Command: `{COMMAND}`",
        f"- Generated at UTC: `{generated_at_utc}`",
        f"- CSV: `{CSV_PATH.relative_to(ROOT)}`",
        f"- Figure: `{SVG_PATH.relative_to(ROOT)}`",
        f"- Log: `{LOG_PATH.relative_to(ROOT)}`",
        f"- N: `{N}`",
        f"- family: `{FAMILY}`",
        f"- d_spacetime: `{D_SPACETIME}`",
        f"- case seed: `{CASE_SEED}`",
        f"- optimizer seeds: "
        f"`{', '.join(str(s) for s in OPTIMIZER_SEEDS)}`",
        f"- T0: `{INITIAL_TEMP}`",
        f"- gamma: `{GAMMA}`",
        f"- h: `{H}`",
        f"- budgets: "
        f"`{', '.join(str(b['budget_label']) for b in BUDGETS)}`",
        f"- reused rows (from `{REUSE_CSV.relative_to(ROOT)}`): "
        f"`{', '.join(f'seed {s}' for s in reused_seeds) if reused_seeds else 'none'}`",
        f"- newly run optimizer seeds: "
        f"`{', '.join(str(s) for s in newly_run_seeds) if newly_run_seeds else 'none'}`",
        f"- skipped rows: "
        f"`{', '.join(f'{s}/{lbl}: {note}' for s, lbl, note in skipped) if skipped else 'none'}`",
        f"- per-run timeout seconds: `{PER_RUN_TIMEOUT_SECONDS:.1f}`",
        f"- total wall runtime seconds: `{total_runtime:.3f}`",
        "- Final causal-order diagnostics use `sim.rold`/`sim.xold`, not "
        "`rnew`/`xnew`.",
        "",
        "## Result table",
        "",
        "| seed | budget | status | final E | interval RMSE | F1 | recall | precision | missing | extra | exact | success | runtime s | note |",
        "| ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: | --- |",
    ]
    ordered = sorted(
        rows,
        key=lambda r: (
            int(r["optimizer_seed"]),
            0 if str(r["budget_label"]) == "short_10_10_4" else 1,
        ),
    )
    for row in ordered:
        lines.append(
            "| {seed} | {label} | {status} | {energy} | {rmse} | {f1} | "
            "{recall} | {precision} | {missing} | {extra} | {exact} | "
            "{success} | {runtime} | {note} |".format(
                seed=int(row["optimizer_seed"]),
                label=str(row["budget_label"]),
                status=str(row["row_status"]),
                energy=_fmt_f(float(row["final_energy"])),
                rmse=_fmt_f(float(row["interval_rmse"]), 4),
                f1=_fmt_f(float(row["causal_f1"])),
                recall=_fmt_f(float(row["causal_recall"])),
                precision=_fmt_f(float(row["causal_precision"])),
                missing=int(row["missing_relations_count"]),
                extra=int(row["extra_relations_count"]),
                exact="yes" if bool(row["exact_match"]) else "no",
                success="yes" if bool(row["success_flag"]) else "no",
                runtime=_fmt_f(float(row["runtime_seconds"]), 3),
                note=str(row["note"]) if str(row["note"]) else "",
            )
        )

    lines += [
        "",
        "## Seed-level comparison (short -> medium)",
        "",
        "| seed | E_short | E_medium | delta_E | F1_short | F1_medium | delta_F1 | recall_short | recall_medium | missing_short | missing_medium | label |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in summary:
        lines.append(
            "| {seed} | {e_short} | {e_medium} | {de} | {f1_short} | "
            "{f1_medium} | {df1} | {r_short} | {r_medium} | {m_short} | "
            "{m_medium} | {label} |".format(
                seed=int(row["optimizer_seed"]),
                e_short=_fmt_f(float(row["E_short"])),
                e_medium=_fmt_f(float(row["E_medium"])),
                de=_fmt_f(float(row["delta_E"])),
                f1_short=_fmt_f(float(row["F1_short"])),
                f1_medium=_fmt_f(float(row["F1_medium"])),
                df1=_fmt_f(float(row["delta_F1"])),
                r_short=_fmt_f(float(row["recall_short"])),
                r_medium=_fmt_f(float(row["recall_medium"])),
                m_short=int(row["missing_short"]),
                m_medium=int(row["missing_medium"]),
                label=str(row["label"]),
            )
        )

    lines += [
        "",
        "Interpretation labels (thresholds: |delta_E| > 1.0, "
        "|delta_F1| > 0.02):",
        "- aligned_improvement: delta_E < -1.0 and delta_F1 > +0.02.",
        "- energy_causality_misalignment: delta_E < -1.0 and "
        "delta_F1 < -0.02.",
        "- causal_improves_without_energy: delta_E > +1.0 and "
        "delta_F1 > +0.02.",
        "- no_clear_change: otherwise.",
        "",
        "## Decision rule",
        "",
        f"- seeds with energy_causality_misalignment: `{misalign}/{total}` "
        f"({', '.join(str(s) for s in misaligned) if misaligned else 'none'}).",
        f"- seeds with aligned_improvement: "
        f"`{len(aligned)}/{total}` "
        f"({', '.join(str(s) for s in aligned) if aligned else 'none'}).",
        f"- seeds with causal_improves_without_energy: "
        f"`{len(causal_only)}/{total}` "
        f"({', '.join(str(s) for s in causal_only) if causal_only else 'none'}).",
        f"- seeds with no_clear_change: `{len(unclear)}/{total}` "
        f"({', '.join(str(s) for s in unclear) if unclear else 'none'}).",
        "- rule: >=4/5 misalignment -> persists; <=2/5 -> seed-specific; "
        "otherwise mixed/inconclusive.",
        f"- outcome: `{outcome}`.",
        "",
        "## Conservative readout",
        "",
        f"- Does the N=30 short->medium misalignment persist across "
        f"optimizer seeds at this case seed/family/schedule? "
        f"`{'yes' if outcome == 'misalignment_persists_across_optimizer_seeds' else 'no' if outcome == 'misalignment_appears_seed_specific' else 'inconclusive'}`.",
        "- Is the N=24 -> N=30 bracket strengthened or weakened by this "
        "probe? "
        f"`{'strengthened (multi-seed N=30 reproduces the misalignment)' if outcome == 'misalignment_persists_across_optimizer_seeds' else 'weakened (N=30 misalignment looks seed-specific in this case seed)' if outcome == 'misalignment_appears_seed_specific' else 'unchanged; pattern is mixed'}`.",
        "- Does this justify a topology / case-seed follow-up at N=30? "
        f"`{'yes; a multi-case-seed probe at N=30 short and medium would test topology specificity next' if outcome == 'misalignment_persists_across_optimizer_seeds' else 'low priority for topology probe before broadening optimizer seeds at the closer N values' if outcome == 'misalignment_appears_seed_specific' else 'yes; mixed outcome calls for both a larger seed pool and a case-seed variation before any claim'}`.",
        "- This probe does not establish: a transition in N, a physical "
        "gamma claim, an embeddability claim, a theorem, or a general "
        "N-scaling statement.",
        "",
        "## Guardrails",
        "",
        "- one N only (`30`).",
        "- one case seed only (`1959`).",
        "- one family only (`minkowski`).",
        "- one schedule only (`T0=100`, `gamma=0.5`, `h=1`).",
        "- two budgets only (`short_10_10_4`, `medium_25_25_8`).",
        "- exploratory only.",
        "- no embeddability claim.",
        "- no physical gamma claim.",
        "- no theorem.",
        "- no general N-scaling claim.",
        f"- optimizer seed {REUSE_OPTIMIZER_SEED} rows are reused from "
        f"`{REUSE_CSV.relative_to(ROOT)}` and were not regenerated here.",
    ]
    MD_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    generated_at_utc = datetime.now(timezone.utc).isoformat(timespec="seconds")
    start = time.perf_counter()
    _log("starting N=30 optimizer-seed budget/misalignment probe")

    rows = read_existing_rows()
    have_rows = {
        (int(row["optimizer_seed"]), str(row["budget_label"])) for row in rows
    }

    # Inject reused seed-1987 rows from the n_ladder probe if not yet
    # present in this probe's CSV.
    reused_rows: list[dict[str, object]] = []
    for reuse_row in _load_reuse_rows():
        key = (int(reuse_row["optimizer_seed"]), str(reuse_row["budget_label"]))
        if key in have_rows:
            continue
        rows.append(reuse_row)
        append_csv_row(reuse_row)
        have_rows.add(key)
        reused_rows.append(reuse_row)
        _log(f"reused {key[1]} optimizer_seed={key[0]} from {REUSE_CSV.name}")

    reused_seeds_sorted = sorted({int(r["optimizer_seed"]) for r in reused_rows})

    case = _make_case()
    newly_run_seeds: list[int] = []
    skipped: list[tuple[int, str, str]] = []

    for optimizer_seed in OPTIMIZER_SEEDS:
        if optimizer_seed == REUSE_OPTIMIZER_SEED and all(
            (optimizer_seed, str(b["budget_label"])) in have_rows for b in BUDGETS
        ):
            continue
        seed_started = False
        for budget in BUDGETS:
            label = str(budget["budget_label"])
            if (optimizer_seed, label) in have_rows:
                continue
            _log(f"attempting {label} optimizer_seed={optimizer_seed}")
            row, error = _run_one_with_timeout(case, budget, optimizer_seed)
            if error is not None:
                _log(f"error at {label} seed {optimizer_seed}: {error}")
                skipped.append((optimizer_seed, label, error))
                break
            assert row is not None
            rows.append(row)
            have_rows.add((optimizer_seed, label))
            append_csv_row(row)
            _log(
                f"{row['row_status']} {label} optimizer_seed={optimizer_seed} "
                f"runtime={float(row['runtime_seconds']):.3f}s"
            )
            seed_started = True
            summary = _seed_summary(rows)
            write_svg(summary)
            write_markdown(
                rows,
                summary,
                reused_seeds_sorted,
                newly_run_seeds,
                skipped,
                generated_at_utc,
                time.perf_counter() - start,
            )
            if str(row["row_status"]) == "timeout":
                break
        if seed_started:
            newly_run_seeds.append(optimizer_seed)

    summary = _seed_summary(rows)
    write_svg(summary)
    write_markdown(
        rows,
        summary,
        reused_seeds_sorted,
        newly_run_seeds,
        skipped,
        generated_at_utc,
        time.perf_counter() - start,
    )
    _log("finished N=30 optimizer-seed budget/misalignment probe")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
