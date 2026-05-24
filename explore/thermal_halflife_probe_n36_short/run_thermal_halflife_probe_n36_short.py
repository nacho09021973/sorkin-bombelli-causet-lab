#!/usr/bin/env python3
"""Exploratory N=36 short-budget thermal half-life mini-probe.

Scans h where gamma = 2^(-1/h) for one known-truth N=36 case, one
optimizer seed, and the short_10_10_4 budget. This is a thermal
schedule diagnostic only, not an embeddability or physics claim.
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
CSV_PATH = OUT_DIR / "thermal_halflife_probe_n36_short.csv"
MD_PATH = OUT_DIR / "thermal_halflife_probe_n36_short.md"
SVG_PATH = OUT_DIR / "thermal_halflife_probe_n36_short.svg"
COMMAND = "python3 explore/thermal_halflife_probe_n36_short/run_thermal_halflife_probe_n36_short.py"

FAMILY = "minkowski"
D_SPACETIME = 2
N = 36
CASE_SEED = 1959
OPTIMIZER_SEED = 1987
TARGET_DIM = D_SPACETIME - 1
INITIAL_TEMP = 100.0
BUDGET_LABEL = "short_10_10_4"
WARMUP_LIMIT = 10
ANNEAL_LIMIT = 10
MAX_DATA = 4
BACKEND = "cpu"
SUCCESS_GAP_THRESHOLD = 1.0
PER_RUN_TIMEOUT_SECONDS = 900.0
HALF_LIVES = (1, 3, 6, 10)

CSV_HEADERS = (
    "h",
    "gamma",
    "family",
    "n",
    "d_spacetime",
    "t0",
    "case_seed",
    "optimizer_seed",
    "budget_label",
    "warmup_limit",
    "anneal_limit",
    "max_data",
    "final_energy",
    "energy_gap",
    "interval_rmse",
    "mm_dim_truth",
    "mm_dim_recovered",
    "causal_precision",
    "causal_recall",
    "causal_f1",
    "missing_relations_count",
    "extra_relations_count",
    "exact_match",
    "success_flag",
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


def gamma_from_half_life(h: int) -> float:
    return 2.0 ** (-1.0 / h)


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


def _run_one(case: vs.SprinkleCase, h: int) -> dict[str, object]:
    gamma = gamma_from_half_life(h)
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

    precision = _safe_div(correct_relations, comparison.total_relations_induced)
    recall = _safe_div(correct_relations, comparison.total_relations_target)
    f1 = _safe_div(2.0 * precision * recall, precision + recall)

    return {
        "h": h,
        "gamma": gamma,
        "family": FAMILY,
        "n": N,
        "d_spacetime": D_SPACETIME,
        "t0": INITIAL_TEMP,
        "case_seed": CASE_SEED,
        "optimizer_seed": OPTIMIZER_SEED,
        "budget_label": BUDGET_LABEL,
        "warmup_limit": WARMUP_LIMIT,
        "anneal_limit": ANNEAL_LIMIT,
        "max_data": MAX_DATA,
        "final_energy": final_energy,
        "energy_gap": energy_gap,
        "interval_rmse": interval_residual,
        "mm_dim_truth": float(case.d_spacetime),
        "mm_dim_recovered": causet_invariants.myrheim_meyer_dimension(case.matrix),
        "causal_precision": precision,
        "causal_recall": recall,
        "causal_f1": f1,
        "missing_relations_count": missing_count,
        "extra_relations_count": extra_count,
        "exact_match": comparison.exact_match,
        "success_flag": success,
        "runtime_seconds": runtime,
    }


def _run_one_worker(queue: mp.Queue, case: vs.SprinkleCase, h: int) -> None:
    try:
        queue.put(("ok", _run_one(case, h)))
    except BaseException as exc:  # pragma: no cover - defensive child process path
        queue.put(("error", repr(exc)))


def _run_one_with_timeout(case: vs.SprinkleCase, h: int) -> tuple[dict[str, object] | None, str | None]:
    queue: mp.Queue = mp.Queue()
    process = mp.Process(target=_run_one_worker, args=(queue, case, h))
    process.start()
    process.join(PER_RUN_TIMEOUT_SECONDS)
    if process.is_alive():
        process.terminate()
        process.join(5.0)
        if process.is_alive():
            process.kill()
            process.join()
        return None, f"timeout after {PER_RUN_TIMEOUT_SECONDS:.1f}s at h={h}"
    if queue.empty():
        return None, f"child exited without result at h={h}"
    status, payload = queue.get()
    if status == "ok":
        return payload, None
    return None, str(payload)


def write_csv(rows: list[dict[str, object]]) -> None:
    with CSV_PATH.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(CSV_HEADERS)
        for row in rows:
            writer.writerow([_fmt(row[header]) for header in CSV_HEADERS])


def _best(rows: list[dict[str, object]], metric: str, *, maximize: bool) -> dict[str, object]:
    return (max if maximize else min)(rows, key=lambda row: float(row[metric]))


def write_svg(rows: list[dict[str, object]]) -> None:
    width = 1120
    height = 620
    margin_left = 76
    margin_right = 34
    margin_top = 62
    margin_bottom = 78
    panel_gap = 48
    panel_width = (width - margin_left - margin_right - 2 * panel_gap) / 3
    panel_height = height - margin_top - margin_bottom
    hs = [float(row["h"]) for row in rows]
    panels = (
        ("causal F1", "causal_f1", "#1f77b4"),
        ("causal recall", "causal_recall", "#2a9d8f"),
        ("missing relations", "missing_relations_count", "#d1495b"),
    )

    def sx(h: float, x0: float) -> float:
        if min(hs) == max(hs):
            return x0 + panel_width / 2.0
        return x0 + (h - min(hs)) * panel_width / (max(hs) - min(hs))

    def sy(value: float, values: list[float]) -> float:
        low = min(values)
        high = max(values)
        if low == high:
            low -= 1.0
            high += 1.0
        return margin_top + (high - value) * panel_height / (high - low)

    chunks = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "  <rect width='100%' height='100%' fill='#f7f4ed'/>",
        f"  <text x='{width / 2:.0f}' y='30' text-anchor='middle' font-family='serif' font-size='22' fill='#222'>N=36 short-budget thermal half-life mini-probe</text>",
    ]
    for index, (title, metric, color) in enumerate(panels):
        x0 = margin_left + index * (panel_width + panel_gap)
        values = [float(row[metric]) for row in rows]
        points = " ".join(
            f"{sx(float(row['h']), x0):.2f},{sy(float(row[metric]), values):.2f}"
            for row in rows
        )
        best = _best(rows, metric, maximize=(metric != "missing_relations_count"))
        chunks.extend(
            [
                f"  <line x1='{x0:.2f}' y1='{height - margin_bottom}' x2='{x0 + panel_width:.2f}' y2='{height - margin_bottom}' stroke='#333' stroke-width='1.5'/>",
                f"  <line x1='{x0:.2f}' y1='{margin_top}' x2='{x0:.2f}' y2='{height - margin_bottom}' stroke='#333' stroke-width='1.5'/>",
                f"  <text x='{x0 + panel_width / 2:.2f}' y='{margin_top - 16}' text-anchor='middle' font-family='monospace' font-size='14' fill='#222'>{title}</text>",
                f"  <polyline fill='none' stroke='{color}' stroke-width='2.5' points='{points}'/>",
            ]
        )
        for row in rows:
            h = float(row["h"])
            value = float(row[metric])
            fill = "#111" if row is best else color
            chunks.append(
                f"  <circle cx='{sx(h, x0):.2f}' cy='{sy(value, values):.2f}' r='4.8' fill='{fill}' opacity='0.9'/>"
            )
        for h in HALF_LIVES:
            chunks.append(
                f"  <text x='{sx(float(h), x0):.2f}' y='{height - margin_bottom + 22}' text-anchor='middle' font-family='monospace' font-size='10' fill='#222'>{h}</text>"
            )
        chunks.append(
            f"  <text x='{x0 + panel_width / 2:.2f}' y='{height - 18}' text-anchor='middle' font-family='monospace' font-size='13' fill='#222'>half-life h</text>"
        )
    chunks.append("</svg>\n")
    SVG_PATH.write_text("\n".join(chunks), encoding="utf-8")


def _md_bool(value: object) -> str:
    return "true" if bool(value) else "false"


def write_markdown(rows: list[dict[str, object]], errors: list[str]) -> None:
    best_f1 = _best(rows, "causal_f1", maximize=True)
    best_recall = _best(rows, "causal_recall", maximize=True)
    best_missing = _best(rows, "missing_relations_count", maximize=False)
    h1 = next((row for row in rows if int(row["h"]) == 1), None)
    improves_over_h1 = (
        any(
            int(row["h"]) != 1
            and (
                float(row["causal_f1"]) > float(h1["causal_f1"])
                or float(row["causal_recall"]) > float(h1["causal_recall"])
                or int(row["missing_relations_count"]) < int(h1["missing_relations_count"])
            )
            for row in rows
        )
        if h1 is not None
        else False
    )
    best_h_gt_1 = (
        int(best_f1["h"]) > 1
        or int(best_recall["h"]) > 1
        or int(best_missing["h"]) > 1
    )
    longer_half_life_text = (
        "In this restricted N=36 mini-probe, at least one primary causal-order criterion is best at h > 1."
        if best_h_gt_1
        else "In this restricted N=36 mini-probe, the primary causal-order criteria do not prefer h > 1."
    )

    lines = [
        "# N=36 Short-Budget Thermal Half-Life Mini-Probe",
        "",
        "## Status",
        "",
        "Exploratory only; not confirmation.",
        "",
        "This probe scans h where gamma = 2^(-1/h), so h is the number of annealing steps required to halve the temperature.",
        "",
        "## Run Configuration",
        "",
        f"- Command: `{COMMAND}`",
        f"- Generated UTC: `{datetime.now(timezone.utc).isoformat(timespec='seconds')}`",
        f"- Family: `{FAMILY}`",
        f"- N: `{N}`",
        f"- spacetime dimension: `{D_SPACETIME}`",
        f"- T0: `{INITIAL_TEMP}`",
        f"- case seed: `{CASE_SEED}`",
        f"- optimizer seed: `{OPTIMIZER_SEED}`",
        f"- budget label: `{BUDGET_LABEL}`",
        f"- budget: warmup `{WARMUP_LIMIT}`, anneal `{ANNEAL_LIMIT}`, max data `{MAX_DATA}`",
        f"- timeout policy: `{PER_RUN_TIMEOUT_SECONDS:.1f}` seconds per h",
        "",
        "## h/gamma Table",
        "",
        "| h | gamma |",
        "|---:|---:|",
    ]
    for h in HALF_LIVES:
        lines.append(f"| {h} | {gamma_from_half_life(h):.10f} |")
    lines.extend(
        [
            "",
            "## Results",
            "",
            "| h | gamma | final_energy | interval_rmse | precision | recall | F1 | missing | extra | exact_match | success | runtime_s |",
            "|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|:---:|---:|",
        ]
    )
    for row in rows:
        lines.append(
            "| {h} | {gamma:.10f} | {final:.6g} | {rmse:.6g} | {precision:.6f} | {recall:.6f} | {f1:.6f} | {missing} | {extra} | {exact} | {success} | {runtime:.3f} |".format(
                h=int(row["h"]),
                gamma=float(row["gamma"]),
                final=float(row["final_energy"]),
                rmse=float(row["interval_rmse"]),
                precision=float(row["causal_precision"]),
                recall=float(row["causal_recall"]),
                f1=float(row["causal_f1"]),
                missing=int(row["missing_relations_count"]),
                extra=int(row["extra_relations_count"]),
                exact=_md_bool(row["exact_match"]),
                success=_md_bool(row["success_flag"]),
                runtime=float(row["runtime_seconds"]),
            )
        )
    lines.extend(
        [
            "",
            "## Required Readout",
            "",
            f"- Which h gives best causal F1? h `{int(best_f1['h'])}` (gamma `{float(best_f1['gamma']):.10f}`, F1 `{float(best_f1['causal_f1']):.6f}`).",
            f"- Which h gives best recall? h `{int(best_recall['h'])}` (gamma `{float(best_recall['gamma']):.10f}`, recall `{float(best_recall['causal_recall']):.6f}`).",
            f"- Which h minimizes missing relations? h `{int(best_missing['h'])}` (gamma `{float(best_missing['gamma']):.10f}`, missing `{int(best_missing['missing_relations_count'])}`).",
            f"- Does any h improve over h=1? `{'yes' if improves_over_h1 else 'no'}` under these causal-order criteria.",
            f"- Is there evidence that N=36 prefers a longer thermal half-life than N=24? {longer_half_life_text} This is only weak exploratory evidence because it is one N=36 case, one optimizer seed, and short budget only.",
            "",
            "## Guardrails",
            "",
            "- Exploratory only.",
            "- One N=36 case only.",
            "- One optimizer seed only.",
            "- Short budget only.",
            "- No embeddability claim.",
            "- No physical constant claim.",
            "- No inference from annealer failure to non-existence.",
            "- No inference from low final energy to manifoldlikeness.",
        ]
    )
    if errors:
        lines.extend(["", "## Run Errors", ""])
        lines.extend(f"- {error}" for error in errors)
    lines.append("")
    MD_PATH.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    case = _make_case()
    rows: list[dict[str, object]] = []
    errors: list[str] = []
    for h in HALF_LIVES:
        row, error = _run_one_with_timeout(case, h)
        if error is not None:
            errors.append(error)
            print(error, file=sys.stderr)
            break
        assert row is not None
        rows.append(row)
        rows.sort(key=lambda item: int(item["h"]))
        write_csv(rows)
        write_svg(rows)
        write_markdown(rows, errors)
        print(
            "h={h} gamma={gamma:.10f} f1={f1:.6f} recall={recall:.6f} missing={missing} runtime={runtime:.3f}s".format(
                h=int(row["h"]),
                gamma=float(row["gamma"]),
                f1=float(row["causal_f1"]),
                recall=float(row["causal_recall"]),
                missing=int(row["missing_relations_count"]),
                runtime=float(row["runtime_seconds"]),
            ),
            flush=True,
        )
    if rows:
        write_csv(rows)
        write_svg(rows)
        write_markdown(rows, errors)
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
