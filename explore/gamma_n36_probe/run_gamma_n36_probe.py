#!/usr/bin/env python3
"""Exploratory N=36 gamma probe for the Bombelli annealer.

Provenance:
- Requested as an exploratory SORKIN-2 gamma-resonance probe.
- Reuses the existing known-truth sprinkling and recovery machinery in
  validation_suite.py.
- Writes only to explore/gamma_n36_probe/ and does not modify benchmark
  or stable documentation artifacts.
"""

from __future__ import annotations

import contextlib
import csv
import io
import math
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import validation_suite as vs  # noqa: E402


OUT_DIR = Path(__file__).resolve().parent
CSV_PATH = OUT_DIR / "gamma_n36_probe.csv"
MD_PATH = OUT_DIR / "gamma_n36_probe.md"
SVG_PATH = OUT_DIR / "gamma_n36_probe.svg"

GAMMAS = (
    0.500000,
    0.530270,
    0.562373,
    0.596419,
    0.632527,
    0.670820,
    0.711432,
    0.754503,
    0.800181,
    0.848624,
    0.900000,
)

# Existing SORKIN-2 diagnostic conventions.
FAMILY = "minkowski"
D_SPACETIME = 2
N = 36
CASE_SEED = 1959
OPTIMIZER_SEED = 1987
TARGET_DIM = D_SPACETIME - 1
INITIAL_TEMP = 100.0
WARMUP_LIMIT = 100
ANNEAL_LIMIT = 100
MAX_DATA = 35
BACKEND = "cpu"
SUCCESS_GAP_THRESHOLD = 1.0


CSV_HEADERS = (
    "family",
    "d_spacetime",
    "n",
    "case_seed",
    "optimizer_seed",
    "target_dim",
    "initial_temp",
    "gamma",
    "warmup_limit",
    "anneal_limit",
    "max_data",
    "initial_energy",
    "warmup_energy",
    "truth_energy",
    "final_energy",
    "energy_gap",
    "interval_rmse",
    "mm_dim_truth",
    "mm_dim_recovered",
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


def _run_one(gamma: float) -> dict[str, object]:
    matrix, points = vs.sprinkle_minkowski_diamond(
        n=N,
        seed=CASE_SEED,
        d_spacetime=D_SPACETIME,
    )
    case = vs.SprinkleCase(
        d_spacetime=D_SPACETIME,
        n=N,
        seed=CASE_SEED,
        matrix=matrix,
        points=points,
    )
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
        "family": FAMILY,
        "d_spacetime": D_SPACETIME,
        "n": N,
        "case_seed": CASE_SEED,
        "optimizer_seed": OPTIMIZER_SEED,
        "target_dim": TARGET_DIM,
        "initial_temp": INITIAL_TEMP,
        "gamma": gamma,
        "warmup_limit": WARMUP_LIMIT,
        "anneal_limit": ANNEAL_LIMIT,
        "max_data": MAX_DATA,
        "initial_energy": result.initial_energy,
        "warmup_energy": result.warmup_energy,
        "truth_energy": result.truth_energy,
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
            writer.writerow([_fmt(row[h]) for h in CSV_HEADERS])


def write_svg(rows: list[dict[str, object]]) -> None:
    width = 920
    height = 540
    margin = 70
    gammas = [float(r["gamma"]) for r in rows]
    finals = [float(r["final_energy"]) for r in rows]
    low = min(finals)
    high = max(finals)
    if low == high:
        low -= 1.0
        high += 1.0

    def sx(x: float) -> float:
        return margin + (x - min(gammas)) * (width - 2 * margin) / (max(gammas) - min(gammas))

    def sy(y: float) -> float:
        return height - margin - (y - low) * (height - 2 * margin) / (high - low)

    points = " ".join(f"{sx(float(r['gamma'])):.2f},{sy(float(r['final_energy'])):.2f}" for r in rows)
    circles = []
    labels = []
    best = min(rows, key=lambda r: float(r["final_energy"]))
    for row in rows:
        gamma = float(row["gamma"])
        final = float(row["final_energy"])
        fill = "#d1495b" if row is best else "#2a9d8f"
        circles.append(
            f"<circle cx='{sx(gamma):.2f}' cy='{sy(final):.2f}' r='5' fill='{fill}' />"
        )
        labels.append(
            f"<text x='{sx(gamma):.2f}' y='{height - margin + 24}' text-anchor='middle' "
            f"font-family='monospace' font-size='10' fill='#222'>{gamma:.3f}</text>"
        )

    svg = f"""<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>
  <rect width='100%' height='100%' fill='#f7f4ed'/>
  <line x1='{margin}' y1='{height - margin}' x2='{width - margin}' y2='{height - margin}' stroke='#333' stroke-width='2'/>
  <line x1='{margin}' y1='{margin}' x2='{margin}' y2='{height - margin}' stroke='#333' stroke-width='2'/>
  <text x='{width/2:.0f}' y='34' text-anchor='middle' font-family='serif' font-size='22' fill='#222'>N=36 gamma probe</text>
  <text x='{width/2:.0f}' y='{height - 18}' text-anchor='middle' font-family='monospace' font-size='15' fill='#222'>gamma</text>
  <text x='22' y='{height/2:.0f}' text-anchor='middle' font-family='monospace' font-size='15' fill='#222' transform='rotate(-90 22 {height/2:.0f})'>final energy</text>
  <text x='{margin}' y='{margin - 20}' font-family='monospace' font-size='12' fill='#444'>lower is better for this diagnostic</text>
  <polyline fill='none' stroke='#1f77b4' stroke-width='2.5' points='{points}' />
  {''.join(circles)}
  {''.join(labels)}
</svg>
"""
    SVG_PATH.write_text(svg, encoding="utf-8")


def _window_label(rows: list[dict[str, object]]) -> str:
    finals = [float(r["final_energy"]) for r in rows]
    best = min(finals)
    worst = max(finals)
    if worst <= 0:
        return "not classifiable from final energy"
    near = [r for r in rows if float(r["final_energy"]) <= best * 1.10 + 1e-9]
    if len(near) >= 3:
        return "broad by the exploratory 10% final-energy criterion"
    return "narrow by the exploratory 10% final-energy criterion"


def write_markdown(rows: list[dict[str, object]]) -> None:
    best = min(rows, key=lambda r: float(r["final_energy"]))
    best_gamma = float(best["gamma"])
    near_label = _window_label(rows)
    half_life = math.log(0.5) / math.log(best_gamma)

    lines = [
        "# Exploratory N=36 Gamma Probe",
        "",
        "**Status:** exploratory result only. This is not confirmation and not publication-ready.",
        "",
        "This is one N=36 known-truth case generated through the existing SORKIN-2 Minkowski sprinkling machinery, using one case seed and one optimizer seed. It is not evidence for a physical constant, not evidence for causal-set embeddability, and only a probe of algorithmic accessibility in the historical Bombelli annealer.",
        "",
        "## Provenance",
        "",
        f"- Output directory: `{OUT_DIR.relative_to(ROOT)}`",
        f"- Runner: `{Path(__file__).relative_to(ROOT)}`",
        f"- CSV: `{CSV_PATH.relative_to(ROOT)}`",
        f"- Plot: `{SVG_PATH.relative_to(ROOT)}`",
        f"- Family: `{FAMILY}`",
        f"- N: `{N}`",
        f"- spacetime dimension: `{D_SPACETIME}`",
        f"- target embedding dimension: `{TARGET_DIM}`",
        f"- case seed: `{CASE_SEED}`",
        f"- optimizer seed: `{OPTIMIZER_SEED}`",
        f"- T0: `{INITIAL_TEMP}`",
        f"- warmup_limit: `{WARMUP_LIMIT}`",
        f"- anneal_limit: `{ANNEAL_LIMIT}`",
        f"- max_data: `{MAX_DATA}`",
        f"- backend: `{BACKEND}`",
        "",
        "## Result",
        "",
        f"- Best gamma by final energy: `{best_gamma:.6f}`",
        f"- Best final energy: `{float(best['final_energy']):.6f}`",
        f"- Best energy gap: `{float(best['energy_gap']):.6f}`",
        f"- Best interval RMSE: `{float(best['interval_rmse']):.6f}`",
        f"- Window readout: {near_label}.",
        f"- Thermal half-life at best gamma: `{half_life:.3f}` annealing temperature steps.",
        "",
        "## Gamma Table",
        "",
        "| gamma | final_energy | energy_gap | interval_rmse | success | runtime_s |",
        "| ---: | ---: | ---: | ---: | :---: | ---: |",
    ]
    for row in rows:
        lines.append(
            "| {gamma:.6f} | {final:.6f} | {gap:.6f} | {rmse:.6f} | {success} | {runtime:.3f} |".format(
                gamma=float(row["gamma"]),
                final=float(row["final_energy"]),
                gap=float(row["energy_gap"]),
                rmse=float(row["interval_rmse"]),
                success="yes" if bool(row["success_flag"]) else "no",
                runtime=float(row["runtime_seconds"]),
            )
        )

    lines += [
        "",
        "## Required Questions",
        "",
        f"- Is the best gamma near 0.8? In this one-run probe, best gamma is `{best_gamma:.6f}`.",
        f"- Is the optimum broad or narrow? {near_label}.",
        f"- Does the result suggest a thermal half-life scale? The best-gamma half-life is `{half_life:.3f}` steps, but this single case and single optimizer seed cannot establish a scale.",
        "- Does the result justify a follow-up sweep over seeds/topologies? It can motivate one only if treated as an algorithmic accessibility follow-up, not as physical evidence.",
        "",
        "## Interpretation Guardrails",
        "",
        "- Exploratory result only.",
        "- One N=36 case.",
        "- One case seed and one optimizer seed.",
        "- Not evidence for a physical constant.",
        "- Not evidence for causal-set embeddability.",
        "- Only a probe of algorithmic accessibility in the historical Bombelli annealer.",
    ]
    MD_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = [_run_one(gamma) for gamma in GAMMAS]
    write_csv(rows)
    write_svg(rows)
    write_markdown(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
