#!/usr/bin/env python3
"""Sweep annealing schedules for a single benchmark causet.

This script compares cooling schedules directly instead of comparing
different seeds only. It keeps the thesis defaults as the baseline and
lets us ask which schedule best reproduces low final energies.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import io
import statistics
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

from cones import ConesSimulator, parse_cones_input


@dataclass
class RunRow:
    initial_temp: float
    cooling_factor: float
    seed: int
    initial_energy: float
    final_energy: float
    warmup_energy: float
    points: int
    output_file: str


@dataclass
class SummaryRow:
    initial_temp: float
    cooling_factor: float
    runs: int
    mean_final: float
    median_final: float
    stdev_final: float
    min_final: float
    max_final: float
    zero_rate: float
    best_seed: int


def linspace(start: float, stop: float, count: int) -> List[float]:
    if count <= 1:
        return [start]
    step = (stop - start) / (count - 1)
    return [start + step * i for i in range(count)]


def summarize(rows: Sequence[RunRow], zero_threshold: float) -> List[SummaryRow]:
    grouped: Dict[Tuple[float, float], List[RunRow]] = {}
    for row in rows:
        grouped.setdefault((row.initial_temp, row.cooling_factor), []).append(row)

    summary: List[SummaryRow] = []
    for (initial_temp, cooling_factor), items in sorted(grouped.items()):
        finals = [r.final_energy for r in items]
        best = min(items, key=lambda r: r.final_energy)
        zero_rate = sum(1 for value in finals if value <= zero_threshold) / len(finals)
        summary.append(
            SummaryRow(
                initial_temp=initial_temp,
                cooling_factor=cooling_factor,
                runs=len(items),
                mean_final=statistics.fmean(finals),
                median_final=statistics.median(finals),
                stdev_final=statistics.pstdev(finals) if len(items) > 1 else 0.0,
                min_final=min(finals),
                max_final=max(finals),
                zero_rate=zero_rate,
                best_seed=best.seed,
            )
        )
    return summary


def color_for(value: float, low: float, high: float) -> str:
    if high <= low:
        return "#2a9d8f"
    t = max(0.0, min(1.0, (value - low) / (high - low)))
    r = int(40 + 190 * t)
    g = int(175 - 110 * t)
    b = int(120 - 70 * t)
    return f"rgb({r},{g},{b})"


def write_summary_csv(summary: Sequence[SummaryRow], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "initial_temp",
                "cooling_factor",
                "runs",
                "mean_final",
                "median_final",
                "stdev_final",
                "min_final",
                "max_final",
                "zero_rate",
                "best_seed",
            ]
        )
        for row in summary:
            writer.writerow(
                [
                    f"{row.initial_temp:.6f}",
                    f"{row.cooling_factor:.6f}",
                    row.runs,
                    f"{row.mean_final:.6f}",
                    f"{row.median_final:.6f}",
                    f"{row.stdev_final:.6f}",
                    f"{row.min_final:.6f}",
                    f"{row.max_final:.6f}",
                    f"{row.zero_rate:.6f}",
                    row.best_seed,
                ]
            )


def write_run_csv(rows: Sequence[RunRow], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "initial_temp",
                "cooling_factor",
                "seed",
                "initial_energy",
                "final_energy",
                "warmup_energy",
                "points",
                "output_file",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    f"{row.initial_temp:.6f}",
                    f"{row.cooling_factor:.6f}",
                    row.seed,
                    f"{row.initial_energy:.6f}",
                    f"{row.final_energy:.6f}",
                    f"{row.warmup_energy:.6f}",
                    row.points,
                    row.output_file,
                ]
            )


def write_heatmap(summary: Sequence[SummaryRow], path: Path) -> None:
    if not summary:
        path.write_text(
            "<svg xmlns='http://www.w3.org/2000/svg' width='800' height='300'></svg>",
            encoding="utf-8",
        )
        return

    temps = sorted({row.initial_temp for row in summary})
    factors = sorted({row.cooling_factor for row in summary})
    low = min(row.mean_final for row in summary)
    high = max(row.mean_final for row in summary)

    width = 120 + 110 * len(factors)
    height = 120 + 78 * len(temps)
    cell_w = 110
    cell_h = 78
    margin_x = 90
    margin_y = 60
    lookup = {(row.initial_temp, row.cooling_factor): row for row in summary}

    cells: List[str] = []
    for yi, temp in enumerate(temps):
        cells.append(
            f"<text x='18' y='{margin_y + yi * cell_h + 34}' font-family='monospace' font-size='13' fill='#222'>T={temp:.1f}</text>"
        )
        for xi, factor in enumerate(factors):
            row = lookup[(temp, factor)]
            x = margin_x + xi * cell_w
            y = margin_y + yi * cell_h
            fill = color_for(row.mean_final, low, high)
            cells.append(
                f"<rect x='{x}' y='{y}' width='{cell_w - 8}' height='{cell_h - 8}' rx='8' fill='{fill}' />"
            )
            cells.append(
                f"<text x='{x + (cell_w - 8)/2:.1f}' y='{y + 28}' text-anchor='middle' font-family='monospace' font-size='12' fill='#111'>{row.mean_final:.3f}</text>"
            )
            cells.append(
                f"<text x='{x + (cell_w - 8)/2:.1f}' y='{y + 48}' text-anchor='middle' font-family='monospace' font-size='11' fill='#111'>z={row.zero_rate:.2f}</text>"
            )

    factor_labels = []
    for xi, factor in enumerate(factors):
        factor_labels.append(
            f"<text x='{margin_x + xi * cell_w + (cell_w - 8)/2:.1f}' y='{margin_y - 18}' text-anchor='middle' font-family='monospace' font-size='13' fill='#222'>{factor:.3f}</text>"
        )

    svg = f"""<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>
  <rect width='100%' height='100%' fill='#f7f3eb'/>
  <text x='{width/2:.0f}' y='32' text-anchor='middle' font-family='serif' font-size='22' fill='#222'>Schedule sweep mean final energy</text>
  <text x='{width/2:.0f}' y='{height - 18}' text-anchor='middle' font-family='monospace' font-size='16' fill='#222'>cooling factor</text>
  {''.join(factor_labels)}
  {''.join(cells)}
  <text x='{width - 220}' y='40' font-family='monospace' font-size='12' fill='#444'>green = lower mean energy</text>
</svg>
"""
    path.write_text(svg, encoding="utf-8")


def write_markdown(
    input_file: Path,
    rows: Sequence[RunRow],
    summary: Sequence[SummaryRow],
    heatmap: Path,
    path: Path,
    zero_threshold: float,
) -> None:
    lines: List[str] = ["# Schedule Sweep", ""]
    lines.append(f"- Input: `{input_file}`")
    lines.append(f"- Heatmap: `{heatmap}`")
    lines.append(f"- Zero threshold: `{zero_threshold}`")
    lines.append("")

    if not rows:
        lines.append("No rows found.")
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return

    best = min(summary, key=lambda row: row.mean_final)
    worst = max(summary, key=lambda row: row.mean_final)

    lines.append("## Summary")
    lines.append("")
    lines.append(
        f"- Best schedule: T0={best.initial_temp:.3f}, cooling={best.cooling_factor:.6f}, mean final energy={best.mean_final:.6f}, zero rate={best.zero_rate:.2f}"
    )
    lines.append(
        f"- Worst schedule: T0={worst.initial_temp:.3f}, cooling={worst.cooling_factor:.6f}, mean final energy={worst.mean_final:.6f}, zero rate={worst.zero_rate:.2f}"
    )
    lines.append("")

    lines.append("## Schedule Summary")
    lines.append("")
    lines.append("| T0 | cooling | runs | mean_final | median_final | stdev_final | min_final | max_final | zero_rate | best_seed |")
    lines.append("| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for row in sorted(summary, key=lambda r: (r.mean_final, r.initial_temp, r.cooling_factor)):
        lines.append(
            f"| {row.initial_temp:.3f} | {row.cooling_factor:.6f} | {row.runs} | {row.mean_final:.6f} | {row.median_final:.6f} | {row.stdev_final:.6f} | {row.min_final:.6f} | {row.max_final:.6f} | {row.zero_rate:.2f} | {row.best_seed} |"
        )

    lines.append("")
    lines.append("## Runs")
    lines.append("")
    lines.append("| T0 | cooling | seed | final_energy | warmup_energy | points | output_file |")
    lines.append("| ---: | ---: | ---: | ---: | ---: | ---: | --- |")
    for row in sorted(rows, key=lambda r: (r.initial_temp, r.cooling_factor, r.final_energy, r.seed)):
        lines.append(
            f"| {row.initial_temp:.3f} | {row.cooling_factor:.6f} | {row.seed} | {row.final_energy:.6f} | {row.warmup_energy:.6f} | {row.points} | {row.output_file} |"
        )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Sweep annealing schedules for cones.py.")
    parser.add_argument("input_file", type=Path, help="benchmark incidence file")
    parser.add_argument("--dim", type=int, required=True, help="space dimension for the embedding")
    parser.add_argument("--seed-start", type=int, default=1959)
    parser.add_argument("--seed-count", type=int, default=12)
    parser.add_argument("--seed-step", type=int, default=1)
    parser.add_argument("--temp-min", type=float, default=50.0)
    parser.add_argument("--temp-max", type=float, default=150.0)
    parser.add_argument("--temp-count", type=int, default=5)
    parser.add_argument("--cooling-min", type=float, default=0.80)
    parser.add_argument("--cooling-max", type=float, default=0.99)
    parser.add_argument("--cooling-count", type=int, default=6)
    parser.add_argument("--warmup-limit", type=int, default=100)
    parser.add_argument("--anneal-limit", type=int, default=100)
    parser.add_argument("--acceptance-scale", type=float, default=4.0)
    parser.add_argument("--backend", choices=["cpu", "cuda", "auto"], default="cpu")
    parser.add_argument("--max-data", type=int, default=35)
    parser.add_argument("--zero-threshold", type=float, default=1e-6)
    parser.add_argument("--run-csv", type=Path, default=Path("schedule_sweep_runs.csv"))
    parser.add_argument("--summary-csv", type=Path, default=Path("schedule_sweep_summary.csv"))
    parser.add_argument("--report-md", type=Path, default=Path("schedule_sweep.md"))
    parser.add_argument("--heatmap-svg", type=Path, default=Path("schedule_sweep.svg"))
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    z = parse_cones_input(args.input_file)
    seeds = [args.seed_start + i * args.seed_step for i in range(args.seed_count)]
    temps = linspace(args.temp_min, args.temp_max, args.temp_count)
    factors = linspace(args.cooling_min, args.cooling_max, args.cooling_count)

    rows: List[RunRow] = []
    with tempfile.TemporaryDirectory(prefix="schedule_sweep_") as tmpdir:
        tmp = Path(tmpdir)
        for initial_temp in temps:
            for cooling_factor in factors:
                for seed in seeds:
                    out_path = tmp / f"T{initial_temp:.3f}_C{cooling_factor:.3f}_S{seed}.out"
                    sim = ConesSimulator(
                        z=z,
                        dim=args.dim,
                        seed=seed,
                        interactive=False,
                        max_data=args.max_data,
                        plot_path=None,
                        backend=args.backend,
                        warmup_limit=args.warmup_limit,
                        anneal_limit=args.anneal_limit,
                        initial_temp=initial_temp,
                        cooling_factor=cooling_factor,
                        acceptance_scale=args.acceptance_scale,
                    )
                    with contextlib.redirect_stdout(io.StringIO()):
                        sim.run(out_path)
                    rows.append(
                        RunRow(
                            initial_temp=initial_temp,
                            cooling_factor=cooling_factor,
                            seed=seed,
                            initial_energy=sim.initial_energy,
                            final_energy=sim.data[-1][1] if sim.data else sim.eave,
                            warmup_energy=sim.warmup_energy,
                            points=len(sim.data),
                            output_file=str(out_path),
                        )
                    )

    summary = summarize(rows, args.zero_threshold)
    write_run_csv(rows, args.run_csv)
    write_summary_csv(summary, args.summary_csv)
    write_heatmap(summary, args.heatmap_svg)
    write_markdown(args.input_file, rows, summary, args.heatmap_svg, args.report_md, args.zero_threshold)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
