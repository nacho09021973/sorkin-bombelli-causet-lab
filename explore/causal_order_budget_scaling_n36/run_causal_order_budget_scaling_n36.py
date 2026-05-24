#!/usr/bin/env python3
"""Budget-aware N=36 causal-order recoverability pilot for SORKIN-2.

Exploratory only: fixed N=36, one known-truth Minkowski case, one
optimizer seed, T0=100, gamma=0.50, varying only annealer budget.
"""

from __future__ import annotations

import contextlib
import csv
import io
import math
import multiprocessing as mp
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
CSV_PATH = OUT_DIR / "causal_order_budget_scaling_n36.csv"
MD_PATH = OUT_DIR / "causal_order_budget_scaling_n36.md"
SVG_PATH = OUT_DIR / "causal_order_budget_scaling_n36.svg"
LOG_PATH = OUT_DIR / "causal_order_budget_scaling_n36.log"
COMMAND = "python3 explore/causal_order_budget_scaling_n36/run_causal_order_budget_scaling_n36.py"

D_SPACETIME = 2
N = 36
CASE_SEED = 1959
OPTIMIZER_SEED = 1987
TARGET_DIM = D_SPACETIME - 1
INITIAL_TEMP = 100.0
GAMMA = 0.50
BACKEND = "cpu"
SUCCESS_GAP_THRESHOLD = 1.0
PER_RUN_TIMEOUT_SECONDS = 900.0
DEFER_LONG_IF_MEDIUM_RUNTIME_EXCEEDS_SECONDS = 300.0
BUDGETS = (
    {"budget_label": "short_10_10_4", "warmup_limit": 10, "anneal_limit": 10, "max_data": 4},
    {"budget_label": "medium_25_25_8", "warmup_limit": 25, "anneal_limit": 25, "max_data": 8},
    {"budget_label": "long_50_50_16", "warmup_limit": 50, "anneal_limit": 50, "max_data": 16},
)

CSV_HEADERS = (
    "budget_label",
    "row_status",
    "gamma",
    "case_seed",
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


def _parse_bool(value: object) -> bool:
    return str(value).strip().lower() == "true"


def _parse_float(value: str) -> float:
    if value == "NA":
        return float("nan")
    return float(value)


def _safe_div(numerator: float, denominator: float) -> float:
    if denominator == 0.0:
        return 0.0
    return numerator / denominator


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


def _blank_row(budget: dict[str, object], row_status: str, runtime: float, note: str) -> dict[str, object]:
    return {
        "budget_label": str(budget["budget_label"]),
        "row_status": row_status,
        "gamma": GAMMA,
        "case_seed": CASE_SEED,
        "optimizer_seed": OPTIMIZER_SEED,
        "n": N,
        "initial_temp": INITIAL_TEMP,
        "warmup_limit": int(budget["warmup_limit"]),
        "anneal_limit": int(budget["anneal_limit"]),
        "max_data": int(budget["max_data"]),
        "final_energy": float("nan"),
        "energy_gap": float("nan"),
        "interval_rmse": float("nan"),
        "mm_dim_truth": float("nan"),
        "mm_dim_recovered": float("nan"),
        "success_flag": False,
        "total_relations_target": 0,
        "total_relations_induced": 0,
        "correct_relations": 0,
        "missing_relations_count": 0,
        "extra_relations_count": 0,
        "exact_match": False,
        "causal_precision": float("nan"),
        "causal_recall": float("nan"),
        "causal_f1": float("nan"),
        "false_positive_relation_rate": float("nan"),
        "false_negative_relation_rate": float("nan"),
        "comparability_fraction_error": float("nan"),
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
        "budget_label": str(budget["budget_label"]),
        "row_status": "completed",
        "gamma": GAMMA,
        "case_seed": CASE_SEED,
        "optimizer_seed": OPTIMIZER_SEED,
        "n": N,
        "initial_temp": INITIAL_TEMP,
        "warmup_limit": int(budget["warmup_limit"]),
        "anneal_limit": int(budget["anneal_limit"]),
        "max_data": int(budget["max_data"]),
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
        "note": "",
    }


def _run_one_worker(queue: mp.Queue, case: vs.SprinkleCase, budget: dict[str, object]) -> None:
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
            budget,
            "timeout",
            runtime,
            f"timeout after {PER_RUN_TIMEOUT_SECONDS:.1f}s",
        ), None
    if queue.empty():
        return None, f"child exited without result at {budget['budget_label']}"
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
                if header in {"budget_label", "row_status", "note"}:
                    parsed[header] = value
                elif header in {"success_flag", "exact_match"}:
                    parsed[header] = _parse_bool(value)
                elif header in {
                    "case_seed",
                    "optimizer_seed",
                    "n",
                    "warmup_limit",
                    "anneal_limit",
                    "max_data",
                    "total_relations_target",
                    "total_relations_induced",
                    "correct_relations",
                    "missing_relations_count",
                    "extra_relations_count",
                }:
                    parsed[header] = int(value)
                else:
                    parsed[header] = _parse_float(value)
            rows.append(parsed)
    return rows


def _completed_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    return [row for row in rows if row["row_status"] == "completed"]


def _aggregates(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    out = []
    completed = _completed_rows(rows)
    for budget in BUDGETS:
        label = str(budget["budget_label"])
        b_rows = [row for row in completed if row["budget_label"] == label]
        if not b_rows:
            continue
        out.append(
            {
                "budget_label": label,
                "runs": len(b_rows),
                "mean_causal_f1": statistics.fmean(float(row["causal_f1"]) for row in b_rows),
                "mean_causal_recall": statistics.fmean(float(row["causal_recall"]) for row in b_rows),
                "mean_causal_precision": statistics.fmean(float(row["causal_precision"]) for row in b_rows),
                "mean_final_energy": statistics.fmean(float(row["final_energy"]) for row in b_rows),
                "mean_interval_rmse": statistics.fmean(float(row["interval_rmse"]) for row in b_rows),
                "missing_total": sum(int(row["missing_relations_count"]) for row in b_rows),
                "extra_total": sum(int(row["extra_relations_count"]) for row in b_rows),
                "mean_missing": statistics.fmean(int(row["missing_relations_count"]) for row in b_rows),
                "mean_extra": statistics.fmean(int(row["extra_relations_count"]) for row in b_rows),
                "exact_matches": sum(1 for row in b_rows if bool(row["exact_match"])),
                "successes": sum(1 for row in b_rows if bool(row["success_flag"])),
                "mean_runtime_seconds": statistics.fmean(float(row["runtime_seconds"]) for row in b_rows),
            }
        )
    return out


def _is_monotonic(values: list[float], *, increasing: bool) -> bool:
    if len(values) < 2:
        return False
    pairs = zip(values, values[1:])
    if increasing:
        return all(right >= left for left, right in pairs)
    return all(right <= left for left, right in pairs)


def _answer(value: bool | None) -> str:
    if value is None:
        return "not assessable"
    return "yes" if value else "no"


def write_svg(rows: list[dict[str, object]], aggregates: list[dict[str, object]]) -> None:
    width = 1120
    height = 720
    margin = 72
    panel_gap = 38
    panel_width = (width - 2 * margin - panel_gap) / 2
    panel_height = (height - 2 * margin - panel_gap) / 2
    labels = [row["budget_label"] for row in aggregates]
    metrics = [
        ("mean causal F1", "mean_causal_f1", False),
        ("mean recall / precision", "mean_causal_recall", False),
        ("mean missing count", "mean_missing", True),
        ("mean final energy", "mean_final_energy", True),
    ]

    def sx(index: int, x0: float) -> float:
        if len(labels) == 1:
            return x0 + panel_width / 2
        return x0 + index * panel_width / (len(labels) - 1)

    def sy(value: float, values: list[float]) -> float:
        low = min(values)
        high = max(values)
        if low == high:
            low -= 1.0
            high += 1.0
        return y0 + panel_height - (value - low) * panel_height / (high - low)

    chunks = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "  <rect width='100%' height='100%' fill='#f7f4ed'/>",
        f"  <text x='{width / 2:.0f}' y='30' text-anchor='middle' font-family='serif' font-size='22' fill='#222'>N=36 causal-order budget scaling pilot</text>",
    ]
    if not aggregates:
        chunks.append("  <text x='560' y='360' text-anchor='middle' font-family='monospace' font-size='16' fill='#444'>no completed budget rows</text>")
        chunks.append("</svg>")
        SVG_PATH.write_text("\n".join(chunks) + "\n", encoding="utf-8")
        return

    for panel_index, (title, metric, lower_better) in enumerate(metrics):
        col = panel_index % 2
        row = panel_index // 2
        x0 = margin + col * (panel_width + panel_gap)
        y0 = margin + row * (panel_height + panel_gap)
        values = [float(item[metric]) for item in aggregates]
        points = " ".join(
            f"{sx(i, x0):.2f},{sy(float(item[metric]), values):.2f}"
            for i, item in enumerate(aggregates)
        )
        chunks.append(f"  <text x='{x0 + panel_width / 2:.1f}' y='{y0 - 12:.1f}' text-anchor='middle' font-family='monospace' font-size='13'>{title}</text>")
        chunks.append(f"  <line x1='{x0}' y1='{y0 + panel_height}' x2='{x0 + panel_width}' y2='{y0 + panel_height}' stroke='#333' stroke-width='2'/>")
        chunks.append(f"  <line x1='{x0}' y1='{y0}' x2='{x0}' y2='{y0 + panel_height}' stroke='#333' stroke-width='2'/>")
        chunks.append(f"  <polyline fill='none' stroke='#1f77b4' stroke-width='2.5' points='{points}'/>")
        best_value = min(values) if lower_better else max(values)
        for i, item in enumerate(aggregates):
            value = float(item[metric])
            fill = "#d1495b" if value == best_value else "#2a9d8f"
            chunks.append(f"  <circle cx='{sx(i, x0):.2f}' cy='{sy(value, values):.2f}' r='5' fill='{fill}'/>")
            chunks.append(f"  <text x='{sx(i, x0):.2f}' y='{y0 + panel_height + 20:.1f}' text-anchor='middle' font-family='monospace' font-size='9'>{i + 1}</text>")
        if metric == "mean_causal_recall":
            precision_values = [float(item["mean_causal_precision"]) for item in aggregates]
            precision_points = " ".join(
                f"{sx(i, x0):.2f},{sy(float(item['mean_causal_precision']), values + precision_values):.2f}"
                for i, item in enumerate(aggregates)
            )
            chunks.append(f"  <polyline fill='none' stroke='#d1495b' stroke-width='2' stroke-dasharray='5 4' points='{precision_points}'/>")
    chunks.append("  <text x='72' y='698' font-family='monospace' font-size='12' fill='#444'>budget labels: 1=10/10/4, 2=25/25/8, 3=50/50/16; red point marks best panel value</text>")
    chunks.append("</svg>")
    SVG_PATH.write_text("\n".join(chunks) + "\n", encoding="utf-8")


def write_markdown(
    rows: list[dict[str, object]],
    aggregates: list[dict[str, object]],
    generated_at_utc: str,
    total_runtime: float,
) -> None:
    completed = _completed_rows(rows)
    status_by_budget = {str(row["budget_label"]): str(row["row_status"]) for row in rows}
    completed_labels = [str(row["budget_label"]) for row in completed]
    deferred_labels = [label for label, status in status_by_budget.items() if status == "deferred"]
    timeout_labels = [label for label, status in status_by_budget.items() if status == "timeout"]
    attempted_labels = [str(row["budget_label"]) for row in rows if row["row_status"] in {"completed", "timeout"}]

    first = aggregates[0] if aggregates else None
    last = aggregates[-1] if aggregates else None
    f1_improves = None if not first or not last or first is last else float(last["mean_causal_f1"]) > float(first["mean_causal_f1"])
    recall_improves = None if not first or not last or first is last else float(last["mean_causal_recall"]) > float(first["mean_causal_recall"])
    missing_reduced = None if not first or not last or first is last else int(last["missing_total"]) < int(first["missing_total"])
    extra_increase = None if not first or not last or first is last else int(last["extra_total"]) > int(first["extra_total"])
    any_exact = any(bool(row["exact_match"]) for row in completed)
    any_success = any(bool(row["success_flag"]) for row in completed)
    f1_values = [float(row["mean_causal_f1"]) for row in aggregates]
    recall_values = [float(row["mean_causal_recall"]) for row in aggregates]
    missing_values = [float(row["mean_missing"]) for row in aggregates]
    monotonic = (
        _is_monotonic(f1_values, increasing=True)
        and _is_monotonic(recall_values, increasing=True)
        and _is_monotonic(missing_values, increasing=False)
        if len(aggregates) >= 2
        else None
    )
    accessibility_effect = (
        bool(f1_improves or recall_improves or missing_reduced)
        if len(aggregates) >= 2
        else None
    )

    lines = [
        "# N=36 Causal-Order Budget Scaling Pilot",
        "",
        "**Status:** exploratory only; not confirmation.",
        "",
        "Does increasing annealer budget improve N=36 causal-order recoverability at fixed T0 and gamma?",
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
        f"- T0: `{INITIAL_TEMP}`",
        f"- gamma: `{GAMMA}`",
        f"- case seed: `{CASE_SEED}`",
        f"- optimizer seed(s): `{OPTIMIZER_SEED}`",
        f"- budgets attempted: `{', '.join(attempted_labels) if attempted_labels else 'none'}`",
        f"- completed rows: `{len(completed)}`",
        f"- completed budgets: `{', '.join(completed_labels) if completed_labels else 'none'}`",
        f"- deferred rows: `{', '.join(deferred_labels) if deferred_labels else 'none'}`",
        f"- timeout rows: `{', '.join(timeout_labels) if timeout_labels else 'none'}`",
        f"- timeout policy: `{PER_RUN_TIMEOUT_SECONDS:.1f}` seconds per budget row",
        f"- long-budget deferral policy: defer long if medium runtime exceeds `{DEFER_LONG_IF_MEDIUM_RUNTIME_EXCEEDS_SECONDS:.1f}` seconds",
        f"- total wall runtime seconds: `{total_runtime:.3f}`",
        "- Final causal-order diagnostics use `sim.rold`/`sim.xold`, not `rnew`/`xnew`.",
        "",
        "## Aggregate Table",
        "",
        "| budget | runs | mean F1 | mean recall | mean precision | mean final E | mean RMSE | missing | extra | mean missing | mean extra | exact | success | mean runtime s |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in aggregates:
        lines.append(
            "| {label} | {runs} | {f1:.6f} | {recall:.6f} | {precision:.6f} | {energy:.6f} | {rmse:.6f} | {missing} | {extra} | {mean_missing:.3f} | {mean_extra:.3f} | {exact} | {success} | {runtime:.3f} |".format(
                label=row["budget_label"],
                runs=int(row["runs"]),
                f1=float(row["mean_causal_f1"]),
                recall=float(row["mean_causal_recall"]),
                precision=float(row["mean_causal_precision"]),
                energy=float(row["mean_final_energy"]),
                rmse=float(row["mean_interval_rmse"]),
                missing=int(row["missing_total"]),
                extra=int(row["extra_total"]),
                mean_missing=float(row["mean_missing"]),
                mean_extra=float(row["mean_extra"]),
                exact=int(row["exact_matches"]),
                success=int(row["successes"]),
                runtime=float(row["mean_runtime_seconds"]),
            )
        )
    lines += [
        "",
        "## Conservative Answers",
        "",
        f"- Does F1 improve with budget? `{_answer(f1_improves)}`.",
        f"- Does recall improve with budget? `{_answer(recall_improves)}`.",
        f"- Are missing relations reduced? `{_answer(missing_reduced)}`.",
        f"- Do extra relations increase? `{_answer(extra_increase)}`.",
        f"- Any exact match? `{'yes' if any_exact else 'no'}`.",
        f"- Any success flag? `{'yes' if any_success else 'no'}`.",
        f"- Is improvement monotonic over attempted budgets? `{_answer(monotonic)}`.",
        f"- Does this look like a budget/accessibility effect? `{_answer(accessibility_effect)}`.",
        "",
        "## Row Status",
        "",
        "| budget | status | runtime s | note |",
        "| --- | --- | ---: | --- |",
    ]
    for row in rows:
        lines.append(
            "| {label} | {status} | {runtime:.3f} | {note} |".format(
                label=row["budget_label"],
                status=row["row_status"],
                runtime=float(row["runtime_seconds"]),
                note=str(row["note"]) if str(row["note"]) else "",
            )
        )
    lines += [
        "",
        "## Guardrails",
        "",
        "- One N=36 case only.",
        "- One/few seeds only.",
        "- Exploratory only.",
        "- No embeddability claim.",
        "- No physical gamma claim.",
        "- No theorem.",
        "- No general N=36 conclusion.",
        "- Further seed/topology sweeps required before promotion.",
    ]
    MD_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    generated_at_utc = datetime.now(timezone.utc).isoformat(timespec="seconds")
    start = time.perf_counter()
    _log("starting N=36 causal-order budget scaling pilot")
    case = _make_case()
    rows = read_existing_rows()
    seen = {str(row["budget_label"]) for row in rows}

    for budget in BUDGETS:
        label = str(budget["budget_label"])
        if label in seen:
            continue
        if label == "long_50_50_16":
            medium = next((row for row in rows if row["budget_label"] == "medium_25_25_8"), None)
            if (
                medium is not None
                and medium["row_status"] == "completed"
                and float(medium["runtime_seconds"]) > DEFER_LONG_IF_MEDIUM_RUNTIME_EXCEEDS_SECONDS
            ):
                row = _blank_row(
                    budget,
                    "deferred",
                    0.0,
                    "deferred because medium runtime exceeded long-budget gate",
                )
                rows.append(row)
                append_csv_row(row)
                seen.add(label)
                _log(f"deferred {label}")
                break
        _log(f"attempting {label}")
        row, error = _run_one_with_timeout(case, budget)
        if error is not None:
            _log(f"error at {label}: {error}")
            return 1
        assert row is not None
        rows.append(row)
        append_csv_row(row)
        seen.add(label)
        _log(f"{row['row_status']} {label} runtime={float(row['runtime_seconds']):.3f}s")
        aggregates = _aggregates(rows)
        write_svg(rows, aggregates)
        write_markdown(rows, aggregates, generated_at_utc, time.perf_counter() - start)
        if row["row_status"] == "timeout":
            break

    aggregates = _aggregates(rows)
    write_svg(rows, aggregates)
    write_markdown(rows, aggregates, generated_at_utc, time.perf_counter() - start)
    _log("finished N=36 causal-order budget scaling pilot")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
