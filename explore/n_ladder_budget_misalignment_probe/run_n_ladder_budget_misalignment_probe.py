#!/usr/bin/env python3
"""N-ladder budget/misalignment probe for SORKIN-2.

Exploratory only: at fixed family/seed/h/gamma/T0, walk a small ladder of
N values and measure whether increasing the historical annealer budget
from short to medium improves causal-order recovery or only lowers
energy. Brackets the gap between the existing N=24 and N=36 probes.

This probe does not rerun N=24 or N=36; their values are read from the
existing CSV artifacts and surfaced only in the Markdown comparison
table.
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

import causet_invariants  # noqa: E402
import cones  # noqa: E402
import validation_suite as vs  # noqa: E402


OUT_DIR = Path(__file__).resolve().parent
CSV_PATH = OUT_DIR / "n_ladder_budget_misalignment_probe.csv"
MD_PATH = OUT_DIR / "n_ladder_budget_misalignment_probe.md"
SVG_PATH = OUT_DIR / "n_ladder_budget_misalignment_probe.svg"
COMMAND = (
    "python3 explore/n_ladder_budget_misalignment_probe/"
    "run_n_ladder_budget_misalignment_probe.py"
)

FAMILY = "minkowski"
D_SPACETIME = 2
CASE_SEED = 1959
OPTIMIZER_SEED = 1987
TARGET_DIM = D_SPACETIME - 1
INITIAL_TEMP = 100.0
GAMMA = 0.50
H = 1
BACKEND = "cpu"
SUCCESS_GAP_THRESHOLD = 1.0
PER_RUN_TIMEOUT_SECONDS = 900.0

REQUESTED_N_VALUES = (12, 18, 30)

BUDGETS = (
    {"budget_label": "short_10_10_4", "warmup_limit": 10, "anneal_limit": 10, "max_data": 4},
    {"budget_label": "medium_25_25_8", "warmup_limit": 25, "anneal_limit": 25, "max_data": 8},
)

REFERENCE_N24_CSV = (
    ROOT / "explore" / "causal_order_budget_scaling_n24" /
    "causal_order_budget_scaling_n24.csv"
)
REFERENCE_N36_CSV = (
    ROOT / "explore" / "causal_order_budget_scaling_n36" /
    "causal_order_budget_scaling_n36.csv"
)

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


def _make_case(n: int) -> vs.SprinkleCase:
    matrix, points = vs.sprinkle_minkowski_diamond(
        n=n,
        seed=CASE_SEED,
        d_spacetime=D_SPACETIME,
    )
    return vs.SprinkleCase(
        d_spacetime=D_SPACETIME,
        n=n,
        seed=CASE_SEED,
        matrix=matrix,
        points=points,
    )


def _blank_row(
    n: int,
    budget: dict[str, object],
    row_status: str,
    runtime: float,
    note: str,
) -> dict[str, object]:
    return {
        "n": n,
        "family": FAMILY,
        "d_spacetime": D_SPACETIME,
        "case_seed": CASE_SEED,
        "optimizer_seed": OPTIMIZER_SEED,
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


def _run_one(case: vs.SprinkleCase, budget: dict[str, object]) -> dict[str, object]:
    start = time.perf_counter()
    with tempfile.TemporaryDirectory() as tmpdir, contextlib.redirect_stdout(io.StringIO()):
        sim = cones.ConesSimulator(
            z=case.matrix,
            dim=TARGET_DIM,
            seed=OPTIMIZER_SEED,
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
        "n": case.n,
        "family": FAMILY,
        "d_spacetime": D_SPACETIME,
        "case_seed": CASE_SEED,
        "optimizer_seed": OPTIMIZER_SEED,
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
) -> None:
    try:
        queue.put(("ok", _run_one(case, budget)))
    except BaseException as exc:  # pragma: no cover - defensive child process path
        queue.put(("error", repr(exc)))


def _run_one_with_timeout(
    case: vs.SprinkleCase,
    budget: dict[str, object],
) -> tuple[dict[str, object] | None, str | None]:
    queue: mp.Queue = mp.Queue()
    process = mp.Process(target=_run_one_worker, args=(queue, case, budget))
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
            case.n,
            budget,
            "timeout",
            runtime,
            f"timeout after {PER_RUN_TIMEOUT_SECONDS:.1f}s",
        ), None
    if queue.empty():
        return None, f"child exited without result at N={case.n} {budget['budget_label']}"
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


def _load_reference(csv_path: Path, n: int, opt_seed: int) -> dict[str, dict[str, float]]:
    """Read existing N=24/N=36 CSV and pick the row matching the optimizer seed.

    Returns a mapping from budget_label -> dict of relevant metrics, or an
    empty dict if the file is absent or the rows are not present.
    """
    out: dict[str, dict[str, float]] = {}
    if not csv_path.exists():
        return out
    with csv_path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            try:
                row_n = int(row.get("n", "0"))
                row_seed = int(row.get("optimizer_seed", "0"))
            except ValueError:
                continue
            if row_n != n or row_seed != opt_seed:
                continue
            label = row.get("budget_label", "")
            if not label:
                continue
            status = row.get("row_status", "completed")
            if status != "completed":
                continue
            try:
                out[label] = {
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
                    "source": str(csv_path.relative_to(ROOT)),
                }
            except ValueError:
                continue
    return out


def _interpretation_label(
    e_short: float,
    e_medium: float,
    f1_short: float,
    f1_medium: float,
    *,
    e_eps: float = 1.0,
    f1_eps: float = 0.01,
) -> str:
    if not (math.isfinite(e_short) and math.isfinite(e_medium)
            and math.isfinite(f1_short) and math.isfinite(f1_medium)):
        return "no_clear_change"
    delta_e = e_medium - e_short
    delta_f1 = f1_medium - f1_short
    e_down = delta_e < -e_eps
    e_up = delta_e > e_eps
    f1_up = delta_f1 > f1_eps
    f1_down = delta_f1 < -f1_eps
    if e_down and f1_up:
        return "aligned_improvement"
    if e_down and f1_down:
        return "energy_causality_misalignment"
    if e_up and f1_up:
        return "causal_improves_without_energy"
    return "no_clear_change"


def _comparison_table(
    rows: list[dict[str, object]],
    ref_n24: dict[str, dict[str, float]],
    ref_n36: dict[str, dict[str, float]],
) -> list[dict[str, object]]:
    """One row per N with short/medium metrics and an interpretation label."""

    def pick_completed(n: int, label: str) -> dict[str, object] | None:
        for row in rows:
            if (
                int(row["n"]) == n
                and str(row["budget_label"]) == label
                and str(row["row_status"]) == "completed"
            ):
                return row
        return None

    def pick_reference(
        ref: dict[str, dict[str, float]], label: str
    ) -> dict[str, float] | None:
        return ref.get(label)

    summary: list[dict[str, object]] = []
    # The ladder: 12, 18, 24 (reference), 30, 36 (reference).
    n_plan: list[tuple[int, str]] = [
        (12, "probe"),
        (18, "probe"),
        (24, "reference"),
        (30, "probe"),
        (36, "reference"),
    ]
    for n, source in n_plan:
        if source == "probe":
            short = pick_completed(n, "short_10_10_4")
            medium = pick_completed(n, "medium_25_25_8")
        else:
            ref = ref_n24 if n == 24 else ref_n36
            short = pick_reference(ref, "short_10_10_4")
            medium = pick_reference(ref, "medium_25_25_8")
        if short is None or medium is None:
            continue
        e_short = float(short["final_energy"])
        e_medium = float(medium["final_energy"])
        f1_short = float(short["causal_f1"])
        f1_medium = float(medium["causal_f1"])
        rec_short = float(short["causal_recall"])
        rec_medium = float(medium["causal_recall"])
        miss_short = int(short["missing_relations_count"])
        miss_medium = int(medium["missing_relations_count"])
        summary.append(
            {
                "n": n,
                "source": source,
                "E_short": e_short,
                "E_medium": e_medium,
                "delta_E": e_medium - e_short,
                "F1_short": f1_short,
                "F1_medium": f1_medium,
                "delta_F1": f1_medium - f1_short,
                "recall_short": rec_short,
                "recall_medium": rec_medium,
                "missing_short": miss_short,
                "missing_medium": miss_medium,
                "interpretation": _interpretation_label(
                    e_short, e_medium, f1_short, f1_medium
                ),
            }
        )
    return summary


def write_svg(comparison: list[dict[str, object]]) -> None:
    width = 1120
    height = 640
    margin_left = 84
    margin_right = 84
    margin_top = 78
    margin_bottom = 96
    plot_w = width - margin_left - margin_right
    plot_h = height - margin_top - margin_bottom

    chunks = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' "
        f"viewBox='0 0 {width} {height}'>",
        "  <rect width='100%' height='100%' fill='#f7f4ed'/>",
        f"  <text x='{width / 2:.0f}' y='34' text-anchor='middle' font-family='serif' "
        "font-size='22' fill='#222'>N-ladder budget/misalignment probe: "
        "delta_F1 and delta_E vs N (short -> medium)</text>",
    ]

    if not comparison:
        chunks.append(
            "  <text x='560' y='320' text-anchor='middle' font-family='monospace' "
            "font-size='16' fill='#444'>no completed comparison rows</text>"
        )
        chunks.append("</svg>")
        SVG_PATH.write_text("\n".join(chunks) + "\n", encoding="utf-8")
        return

    ns = [int(row["n"]) for row in comparison]
    delta_f1 = [float(row["delta_F1"]) for row in comparison]
    delta_e = [float(row["delta_E"]) for row in comparison]
    sources = [str(row["source"]) for row in comparison]

    n_min = min(ns)
    n_max = max(ns)
    n_span = max(n_max - n_min, 1)

    def x_of(n: int) -> float:
        return margin_left + (n - n_min) * plot_w / n_span

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

    # Axes box.
    chunks.append(
        f"  <rect x='{margin_left}' y='{margin_top}' width='{plot_w}' "
        f"height='{plot_h}' fill='none' stroke='#333' stroke-width='2'/>"
    )

    # Zero line for delta values.
    y_zero_f1 = y_of_f1(0.0)
    chunks.append(
        f"  <line x1='{margin_left}' y1='{y_zero_f1:.2f}' "
        f"x2='{margin_left + plot_w}' y2='{y_zero_f1:.2f}' "
        "stroke='#999' stroke-width='1' stroke-dasharray='4 4'/>"
    )

    # X-axis tick marks at each N.
    for n in ns:
        x = x_of(n)
        chunks.append(
            f"  <line x1='{x:.2f}' y1='{margin_top + plot_h}' "
            f"x2='{x:.2f}' y2='{margin_top + plot_h + 6}' "
            "stroke='#333' stroke-width='1.5'/>"
        )
        chunks.append(
            f"  <text x='{x:.2f}' y='{margin_top + plot_h + 22}' text-anchor='middle' "
            f"font-family='monospace' font-size='12' fill='#222'>N={n}</text>"
        )

    # X-axis label.
    chunks.append(
        f"  <text x='{margin_left + plot_w / 2:.2f}' "
        f"y='{margin_top + plot_h + 56}' text-anchor='middle' "
        "font-family='monospace' font-size='13' fill='#222'>N</text>"
    )

    # Left y-axis ticks for delta_F1 (5 levels).
    for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
        v = f1_min + frac * (f1_max - f1_min)
        y = y_of_f1(v)
        chunks.append(
            f"  <line x1='{margin_left - 6}' y1='{y:.2f}' "
            f"x2='{margin_left}' y2='{y:.2f}' stroke='#333' stroke-width='1'/>"
        )
        chunks.append(
            f"  <text x='{margin_left - 10}' y='{y + 4:.2f}' text-anchor='end' "
            f"font-family='monospace' font-size='11' fill='#1f77b4'>{v:.3f}</text>"
        )

    # Right y-axis ticks for delta_E (5 levels).
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

    # Y-axis labels.
    chunks.append(
        f"  <text x='24' y='{margin_top + plot_h / 2:.2f}' text-anchor='middle' "
        "font-family='monospace' font-size='13' fill='#1f77b4' "
        f"transform='rotate(-90 24 {margin_top + plot_h / 2:.2f})'>"
        "delta_F1 = F1_medium - F1_short</text>"
    )
    chunks.append(
        f"  <text x='{width - 24}' y='{margin_top + plot_h / 2:.2f}' "
        "text-anchor='middle' font-family='monospace' font-size='13' fill='#d1495b' "
        f"transform='rotate(90 {width - 24} {margin_top + plot_h / 2:.2f})'>"
        "delta_E = E_medium - E_short</text>"
    )

    # Plot delta_F1 polyline (blue).
    f1_points = " ".join(
        f"{x_of(n):.2f},{y_of_f1(v):.2f}" for n, v in zip(ns, delta_f1)
    )
    chunks.append(
        f"  <polyline fill='none' stroke='#1f77b4' stroke-width='2.5' "
        f"points='{f1_points}'/>"
    )

    # Plot delta_E polyline (red, dashed).
    e_points = " ".join(
        f"{x_of(n):.2f},{y_of_e(v):.2f}" for n, v in zip(ns, delta_e)
    )
    chunks.append(
        f"  <polyline fill='none' stroke='#d1495b' stroke-width='2' "
        f"stroke-dasharray='6 5' points='{e_points}'/>"
    )

    # Markers: filled circle for delta_F1, ring for delta_E. Reference points hollow.
    for n, vf1, ve, src in zip(ns, delta_f1, delta_e, sources):
        x = x_of(n)
        yf = y_of_f1(vf1)
        ye = y_of_e(ve)
        if src == "probe":
            chunks.append(
                f"  <circle cx='{x:.2f}' cy='{yf:.2f}' r='6' fill='#1f77b4' "
                "stroke='#1f77b4' stroke-width='1'/>"
            )
            chunks.append(
                f"  <circle cx='{x:.2f}' cy='{ye:.2f}' r='6' fill='none' "
                "stroke='#d1495b' stroke-width='2'/>"
            )
        else:
            chunks.append(
                f"  <circle cx='{x:.2f}' cy='{yf:.2f}' r='6' fill='#f7f4ed' "
                "stroke='#1f77b4' stroke-width='2'/>"
            )
            chunks.append(
                f"  <circle cx='{x:.2f}' cy='{ye:.2f}' r='6' fill='#f7f4ed' "
                "stroke='#d1495b' stroke-width='2' stroke-dasharray='3 2'/>"
            )

    # Legend.
    lg_x = margin_left + 16
    lg_y = margin_top + 16
    chunks.append(
        f"  <rect x='{lg_x - 8}' y='{lg_y - 14}' width='320' height='66' "
        "fill='#ffffff' fill-opacity='0.82' stroke='#aaa' stroke-width='1'/>"
    )
    chunks.append(
        f"  <line x1='{lg_x}' y1='{lg_y}' x2='{lg_x + 28}' y2='{lg_y}' "
        "stroke='#1f77b4' stroke-width='2.5'/>"
    )
    chunks.append(
        f"  <text x='{lg_x + 36}' y='{lg_y + 4}' font-family='monospace' "
        "font-size='12' fill='#222'>delta_F1 (left axis, blue)</text>"
    )
    chunks.append(
        f"  <line x1='{lg_x}' y1='{lg_y + 20}' x2='{lg_x + 28}' y2='{lg_y + 20}' "
        "stroke='#d1495b' stroke-width='2' stroke-dasharray='6 5'/>"
    )
    chunks.append(
        f"  <text x='{lg_x + 36}' y='{lg_y + 24}' font-family='monospace' "
        "font-size='12' fill='#222'>delta_E (right axis, red)</text>"
    )
    chunks.append(
        f"  <text x='{lg_x}' y='{lg_y + 44}' font-family='monospace' "
        "font-size='11' fill='#444'>filled = this probe; hollow = N=24/N=36 reference"
        "</text>"
    )

    # Footnote.
    chunks.append(
        f"  <text x='{margin_left}' y='{height - 18}' font-family='monospace' "
        "font-size='11' fill='#444'>Exploratory only; family minkowski; "
        f"d_spacetime {D_SPACETIME}; case seed {CASE_SEED}; optimizer seed "
        f"{OPTIMIZER_SEED}; T0 {INITIAL_TEMP:.0f}; gamma {GAMMA}; h {H}.</text>"
    )
    chunks.append("</svg>")
    SVG_PATH.write_text("\n".join(chunks) + "\n", encoding="utf-8")


def _fmt_f(value: float, digits: int = 6) -> str:
    if not math.isfinite(value):
        return "NA"
    return f"{value:.{digits}f}"


def write_markdown(
    rows: list[dict[str, object]],
    comparison: list[dict[str, object]],
    completed_ns: list[int],
    skipped_ns: list[tuple[int, str]],
    ref_n24: dict[str, dict[str, float]],
    ref_n36: dict[str, dict[str, float]],
    generated_at_utc: str,
    total_runtime: float,
) -> None:
    completed_short = sorted(
        int(row["n"]) for row in rows
        if str(row["row_status"]) == "completed"
        and str(row["budget_label"]) == "short_10_10_4"
    )
    completed_medium = sorted(
        int(row["n"]) for row in rows
        if str(row["row_status"]) == "completed"
        and str(row["budget_label"]) == "medium_25_25_8"
    )
    aligned = [int(row["n"]) for row in comparison
               if row["interpretation"] == "aligned_improvement"]
    misaligned = [int(row["n"]) for row in comparison
                  if row["interpretation"] == "energy_causality_misalignment"]
    causal_only = [int(row["n"]) for row in comparison
                   if row["interpretation"] == "causal_improves_without_energy"]
    unclear = [int(row["n"]) for row in comparison
               if row["interpretation"] == "no_clear_change"]
    n24_label = next(
        (str(row["interpretation"]) for row in comparison if int(row["n"]) == 24),
        "not in comparison",
    )
    n36_label = next(
        (str(row["interpretation"]) for row in comparison if int(row["n"]) == 36),
        "not in comparison",
    )
    transition = (
        n24_label == "aligned_improvement"
        and n36_label == "energy_causality_misalignment"
    )

    lines = [
        "# N-Ladder Budget/Misalignment Probe",
        "",
        "exploratory only; not confirmation.",
        "",
        "Across a small ladder of N values at fixed family/seed/h/gamma/T0, "
        "when does increasing the historical annealer budget from short to "
        "medium improve causal-order recovery, and when does it lower final "
        "energy while worsening causal-order recovery?",
        "",
        "## Provenance",
        "",
        f"- Output directory: `{OUT_DIR.relative_to(ROOT)}`",
        f"- Runner: `{Path(__file__).relative_to(ROOT)}`",
        f"- Command: `{COMMAND}`",
        f"- Generated at UTC: `{generated_at_utc}`",
        f"- CSV: `{CSV_PATH.relative_to(ROOT)}`",
        f"- Figure: `{SVG_PATH.relative_to(ROOT)}`",
        f"- N values requested: `{', '.join(str(n) for n in REQUESTED_N_VALUES)}`",
        f"- Completed N values: "
        f"`{', '.join(str(n) for n in completed_ns) if completed_ns else 'none'}`",
        f"- Skipped N values: "
        f"`{', '.join(f'{n} ({reason})' for n, reason in skipped_ns) if skipped_ns else 'none'}`",
        f"- Completed short budget rows: "
        f"`{', '.join(f'N={n}' for n in completed_short) if completed_short else 'none'}`",
        f"- Completed medium budget rows: "
        f"`{', '.join(f'N={n}' for n in completed_medium) if completed_medium else 'none'}`",
        f"- T0: `{INITIAL_TEMP}`",
        f"- gamma: `{GAMMA}`",
        f"- h: `{H}`",
        f"- case seed: `{CASE_SEED}`",
        f"- optimizer seed: `{OPTIMIZER_SEED}`",
        f"- family: `{FAMILY}`",
        f"- d_spacetime: `{D_SPACETIME}`",
        f"- budgets attempted: "
        f"`{', '.join(str(b['budget_label']) for b in BUDGETS)}`",
        f"- per-run timeout seconds: `{PER_RUN_TIMEOUT_SECONDS:.1f}`",
        f"- total wall runtime seconds: `{total_runtime:.3f}`",
        f"- N=24 reference CSV: `{REFERENCE_N24_CSV.relative_to(ROOT)}` "
        f"(rows found for optimizer seed {OPTIMIZER_SEED}: "
        f"{', '.join(sorted(ref_n24.keys())) if ref_n24 else 'none'})",
        f"- N=36 reference CSV: `{REFERENCE_N36_CSV.relative_to(ROOT)}` "
        f"(rows found for optimizer seed {OPTIMIZER_SEED}: "
        f"{', '.join(sorted(ref_n36.keys())) if ref_n36 else 'none'})",
        "- Final causal-order diagnostics use `sim.rold`/`sim.xold`, not `rnew`/`xnew`.",
        "",
        "## Result table by N and budget",
        "",
        "| N | source | budget | status | final E | interval RMSE | F1 | recall | precision | missing | extra | exact | success | runtime s |",
        "| ---: | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: |",
    ]

    def row_line(
        n: int,
        source: str,
        label: str,
        status: str,
        final_energy: float,
        rmse: float,
        f1: float,
        recall: float,
        precision: float,
        missing: int,
        extra: int,
        exact: bool,
        success: bool,
        runtime: float,
    ) -> str:
        return (
            f"| {n} | {source} | {label} | {status} | {_fmt_f(final_energy)} | "
            f"{_fmt_f(rmse, 4)} | {_fmt_f(f1)} | {_fmt_f(recall)} | "
            f"{_fmt_f(precision)} | {missing} | {extra} | "
            f"{'yes' if exact else 'no'} | {'yes' if success else 'no'} | "
            f"{_fmt_f(runtime, 3)} |"
        )

    # Probe rows (N=12, 18, 30) in the ladder order, then references for context.
    ladder_order = [12, 18, 24, 30, 36]
    for n in ladder_order:
        if n in REQUESTED_N_VALUES:
            for budget in BUDGETS:
                label = str(budget["budget_label"])
                row = next(
                    (r for r in rows
                     if int(r["n"]) == n and str(r["budget_label"]) == label),
                    None,
                )
                if row is None:
                    lines.append(
                        f"| {n} | probe | {label} | not_run | NA | NA | NA | NA | "
                        "NA | 0 | 0 | no | no | NA |"
                    )
                    continue
                lines.append(
                    row_line(
                        n,
                        "probe",
                        label,
                        str(row["row_status"]),
                        float(row["final_energy"]),
                        float(row["interval_rmse"]),
                        float(row["causal_f1"]),
                        float(row["causal_recall"]),
                        float(row["causal_precision"]),
                        int(row["missing_relations_count"]),
                        int(row["extra_relations_count"]),
                        bool(row["exact_match"]),
                        bool(row["success_flag"]),
                        float(row["runtime_seconds"]),
                    )
                )
        elif n == 24:
            for label in ("short_10_10_4", "medium_25_25_8"):
                ref = ref_n24.get(label)
                if ref is None:
                    lines.append(
                        f"| 24 | reference | {label} | unavailable | NA | NA | NA | NA | "
                        "NA | 0 | 0 | no | no | NA |"
                    )
                    continue
                lines.append(
                    row_line(
                        24,
                        "reference",
                        label,
                        "completed",
                        float(ref["final_energy"]),
                        float(ref["interval_rmse"]),
                        float(ref["causal_f1"]),
                        float(ref["causal_recall"]),
                        float(ref["causal_precision"]),
                        int(ref["missing_relations_count"]),
                        int(ref["extra_relations_count"]),
                        bool(ref["exact_match"]),
                        bool(ref["success_flag"]),
                        float(ref["runtime_seconds"]),
                    )
                )
        elif n == 36:
            for label in ("short_10_10_4", "medium_25_25_8"):
                ref = ref_n36.get(label)
                if ref is None:
                    lines.append(
                        f"| 36 | reference | {label} | unavailable | NA | NA | NA | NA | "
                        "NA | 0 | 0 | no | no | NA |"
                    )
                    continue
                lines.append(
                    row_line(
                        36,
                        "reference",
                        label,
                        "completed",
                        float(ref["final_energy"]),
                        float(ref["interval_rmse"]),
                        float(ref["causal_f1"]),
                        float(ref["causal_recall"]),
                        float(ref["causal_precision"]),
                        int(ref["missing_relations_count"]),
                        int(ref["extra_relations_count"]),
                        bool(ref["exact_match"]),
                        bool(ref["success_flag"]),
                        float(ref["runtime_seconds"]),
                    )
                )

    lines += [
        "",
        "## Derived comparison table (short -> medium)",
        "",
        "| N | source | E_short | E_medium | delta_E | F1_short | F1_medium | delta_F1 | recall_short | recall_medium | missing_short | missing_medium | interpretation |",
        "| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in comparison:
        lines.append(
            "| {n} | {source} | {e_short} | {e_medium} | {de} | {f1_short} | "
            "{f1_medium} | {df1} | {r_short} | {r_medium} | {m_short} | "
            "{m_medium} | {interp} |".format(
                n=int(row["n"]),
                source=str(row["source"]),
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
                interp=str(row["interpretation"]),
            )
        )

    lines += [
        "",
        "Interpretation labels:",
        "- aligned_improvement: energy decreases and F1 increases (threshold: |delta_E| > 1.0, |delta_F1| > 0.01).",
        "- energy_causality_misalignment: energy decreases and F1 decreases.",
        "- causal_improves_without_energy: energy increases and F1 increases.",
        "- no_clear_change: otherwise or change too small/ambiguous.",
        "",
        "## Conservative readout",
        "",
        f"- N values with aligned improvement: "
        f"`{', '.join(str(n) for n in aligned) if aligned else 'none'}`.",
        f"- N values with energy/causality misalignment: "
        f"`{', '.join(str(n) for n in misaligned) if misaligned else 'none'}`.",
        f"- N values with causal-improves-without-energy: "
        f"`{', '.join(str(n) for n in causal_only) if causal_only else 'none'}`.",
        f"- N values with no clear change: "
        f"`{', '.join(str(n) for n in unclear) if unclear else 'none'}`.",
        f"- N=24 interpretation label: `{n24_label}`.",
        f"- N=36 interpretation label: `{n36_label}`.",
        f"- Apparent transition between N=24 and N=36 "
        f"(aligned at 24, misaligned at 36)? `{'yes' if transition else 'no'}`.",
        "- Does this exploratory ladder support looking for geometry-dependent "
        "thermal mobility? exploratory hint only; one case seed and one "
        "optimizer seed per N is not enough to establish geometry dependence, "
        "but the table provides a working bracket worth examining with more "
        "seeds before any claim.",
        "- Does this exploratory ladder justify more seeds or topologies? "
        "yes as an exploratory next step; broadening the seed pool and adding "
        "non-Minkowski topologies is appropriate before any promotion.",
        "",
        "## Guardrails",
        "",
        "- exploratory only.",
        "- one family only (`minkowski`).",
        "- one case seed only (`1959`).",
        "- one optimizer seed only (`1987`).",
        "- two budgets only (`short_10_10_4`, `medium_25_25_8`).",
        "- no embeddability claim.",
        "- no physical gamma claim.",
        "- no theorem.",
        "- no general N-scaling claim.",
        "- N=24 and N=36 rows in the comparison table are read from existing "
        "exploratory CSVs and were not regenerated here.",
    ]
    MD_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    generated_at_utc = datetime.now(timezone.utc).isoformat(timespec="seconds")
    start = time.perf_counter()
    rows = read_existing_rows()
    seen = {
        (int(row["n"]), str(row["budget_label"])): row for row in rows
    }
    completed_ns: list[int] = []
    skipped_ns: list[tuple[int, str]] = []

    ref_n24 = _load_reference(REFERENCE_N24_CSV, n=24, opt_seed=OPTIMIZER_SEED)
    ref_n36 = _load_reference(REFERENCE_N36_CSV, n=36, opt_seed=OPTIMIZER_SEED)

    for n in REQUESTED_N_VALUES:
        try:
            case = _make_case(n)
        except Exception as exc:  # noqa: BLE001
            skipped_ns.append((n, f"case construction failed: {exc!r}"))
            continue
        n_ok = True
        for budget in BUDGETS:
            label = str(budget["budget_label"])
            if (n, label) in seen:
                continue
            row, error = _run_one_with_timeout(case, budget)
            if error is not None:
                n_ok = False
                skipped_ns.append(
                    (n, f"runner error at {label}: {error}")
                )
                break
            assert row is not None
            rows.append(row)
            seen[(n, label)] = row
            append_csv_row(row)
            comparison = _comparison_table(rows, ref_n24, ref_n36)
            write_svg(comparison)
            write_markdown(
                rows,
                comparison,
                completed_ns,
                skipped_ns,
                ref_n24,
                ref_n36,
                generated_at_utc,
                time.perf_counter() - start,
            )
        if n_ok:
            n_completed = all(
                (n, str(b["budget_label"])) in seen
                and str(seen[(n, str(b["budget_label"]))]["row_status"]) == "completed"
                for b in BUDGETS
            )
            if n_completed:
                completed_ns.append(n)

    comparison = _comparison_table(rows, ref_n24, ref_n36)
    write_svg(comparison)
    write_markdown(
        rows,
        comparison,
        completed_ns,
        skipped_ns,
        ref_n24,
        ref_n36,
        generated_at_utc,
        time.perf_counter() - start,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
