#!/usr/bin/env python3
"""Analyze sweep results from ``cones.py``.

This script stays dependency-free and turns a sweep CSV into:
- a concise markdown report
- a simple SVG plot of final energy by run
- an optional CSV with derived statistics is intentionally omitted for now

The focus is on quickly comparing many seeds without needing a notebook.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence


@dataclass
class SweepRow:
    seed: int
    n: int
    dim: int
    initial_energy: float
    final_energy: float
    warmup_energy: float
    points: int
    input_file: str
    output_file: str
    plot_file: str
    causet_plot_file: str


def parse_float(value: str) -> float:
    return float(value) if value else 0.0


def parse_int(value: str) -> int:
    return int(value) if value else 0


def read_sweep_csv(path: Path) -> List[SweepRow]:
    rows: List[SweepRow] = []
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for raw in reader:
            rows.append(
                SweepRow(
                    seed=parse_int(raw.get("seed", "")),
                    n=parse_int(raw.get("n", "")),
                    dim=parse_int(raw.get("dim", "")),
                    initial_energy=parse_float(raw.get("initial_energy", "")),
                    final_energy=parse_float(raw.get("final_energy", "")),
                    warmup_energy=parse_float(raw.get("warmup_energy", "")),
                    points=parse_int(raw.get("points", "")),
                    input_file=raw.get("input_file", "") or "",
                    output_file=raw.get("output_file", "") or "",
                    plot_file=raw.get("plot_file", "") or "",
                    causet_plot_file=raw.get("causet_plot_file", "") or "",
                )
            )
    return rows


def summarize(rows: Sequence[SweepRow]) -> dict:
    if not rows:
        return {}

    final_energies = [r.final_energy for r in rows]
    initial_energies = [r.initial_energy for r in rows]
    warmup_energies = [r.warmup_energy for r in rows]
    points = [r.points for r in rows]

    best = min(rows, key=lambda r: r.final_energy)
    worst = max(rows, key=lambda r: r.final_energy)

    return {
        "count": len(rows),
        "mean_final": statistics.fmean(final_energies),
        "median_final": statistics.median(final_energies),
        "stdev_final": statistics.pstdev(final_energies) if len(rows) > 1 else 0.0,
        "min_final": best.final_energy,
        "max_final": worst.final_energy,
        "mean_initial": statistics.fmean(initial_energies),
        "mean_warmup": statistics.fmean(warmup_energies),
        "mean_points": statistics.fmean(points),
        "best_seed": best.seed,
        "worst_seed": worst.seed,
    }


def rank_rows(rows: Sequence[SweepRow]) -> List[SweepRow]:
    return sorted(rows, key=lambda r: (r.final_energy, r.seed))


def svg_line_plot(rows: Sequence[SweepRow], path: Path) -> None:
    width = 1000
    height = 600
    margin = 70
    if not rows:
        path.write_text(
            "<svg xmlns='http://www.w3.org/2000/svg' width='1000' height='600'>"
            "<text x='40' y='50' font-family='monospace' font-size='20'>No data</text>"
            "</svg>",
            encoding="utf-8",
        )
        return

    xs = list(range(len(rows)))
    ys = [r.final_energy for r in rows]
    x_min, x_max = 0, max(xs) if xs else 1
    y_min, y_max = min(ys), max(ys)
    if y_min == y_max:
        y_min -= 1.0
        y_max += 1.0
    if x_max == x_min:
        x_max += 1

    def sx(x: float) -> float:
        return margin + (x - x_min) * (width - 2 * margin) / (x_max - x_min)

    def sy(y: float) -> float:
        return height - margin - (y - y_min) * (height - 2 * margin) / (y_max - y_min)

    polyline = " ".join(f"{sx(i):.2f},{sy(y):.2f}" for i, y in zip(xs, ys))
    circles = "\n".join(
        f"<circle cx='{sx(i):.2f}' cy='{sy(row.final_energy):.2f}' r='4' fill='#d1495b' />"
        for i, row in enumerate(rows)
    )
    labels = "\n".join(
        f"<text x='{sx(i):.2f}' y='{height - 40}' text-anchor='middle' font-family='monospace' font-size='12' fill='#444'>{row.seed}</text>"
        for i, row in enumerate(rows)
    )

    svg = f"""<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>
  <rect width='100%' height='100%' fill='#fbf7f0'/>
  <line x1='{margin}' y1='{height - margin}' x2='{width - margin}' y2='{height - margin}' stroke='#222' stroke-width='2'/>
  <line x1='{margin}' y1='{margin}' x2='{margin}' y2='{height - margin}' stroke='#222' stroke-width='2'/>
  <text x='{width/2:.0f}' y='34' text-anchor='middle' font-family='serif' font-size='22' fill='#222'>Final energy by sweep run</text>
  <text x='{width/2:.0f}' y='{height - 12}' text-anchor='middle' font-family='monospace' font-size='16' fill='#222'>run index</text>
  <text x='18' y='{height/2:.0f}' text-anchor='middle' font-family='monospace' font-size='16' fill='#222' transform='rotate(-90 18 {height/2:.0f})'>final energy</text>
  <polyline fill='none' stroke='#1d4ed8' stroke-width='2.5' points='{polyline}' />
  {circles}
  {labels}
</svg>
"""
    path.write_text(svg, encoding="utf-8")


def write_markdown(
    rows: Sequence[SweepRow],
    summary: dict,
    path: Path,
    top_k: int,
    plot_path: Path | None,
) -> None:
    lines: List[str] = ["# Sweep analysis", ""]
    if not rows:
        lines.append("No rows found.")
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return

    lines.append("## Summary")
    lines.append("")
    lines.append(f"- Runs: {summary['count']}")
    lines.append(f"- Mean final energy: {summary['mean_final']:.6f}")
    lines.append(f"- Median final energy: {summary['median_final']:.6f}")
    lines.append(f"- Std dev final energy: {summary['stdev_final']:.6f}")
    lines.append(f"- Min final energy: {summary['min_final']:.6f} (seed {summary['best_seed']})")
    lines.append(f"- Max final energy: {summary['max_final']:.6f} (seed {summary['worst_seed']})")
    lines.append(f"- Mean initial energy: {summary['mean_initial']:.6f}")
    lines.append(f"- Mean warmup energy: {summary['mean_warmup']:.6f}")
    lines.append(f"- Mean number of recorded points: {summary['mean_points']:.2f}")
    if plot_path is not None:
        lines.append(f"- Plot: `{plot_path}`")

    lines.append("")
    lines.append("## Best Runs")
    lines.append("")
    lines.append("| seed | final_energy | initial_energy | warmup_energy | points | output_file |")
    lines.append("| --- | ---: | ---: | ---: | ---: | --- |")
    for row in rank_rows(rows)[:top_k]:
        lines.append(
            f"| {row.seed} | {row.final_energy:.6f} | {row.initial_energy:.6f} | {row.warmup_energy:.6f} | {row.points} | {row.output_file} |"
        )

    lines.append("")
    lines.append("## Full Run List")
    lines.append("")
    lines.append("| run | seed | final_energy | plot | causet |")
    lines.append("| ---: | ---: | ---: | --- | --- |")
    for idx, row in enumerate(rows):
        lines.append(
            f"| {idx} | {row.seed} | {row.final_energy:.6f} | {row.plot_file or '-'} | {row.causet_plot_file or '-'} |"
        )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Analyze cones.py sweep outputs.")
    parser.add_argument("csv", type=Path, help="sweep CSV produced by cones.py")
    parser.add_argument("--report-md", type=Path, default=Path("sweep_analysis.md"), help="markdown report path")
    parser.add_argument("--plot-svg", type=Path, default=Path("sweep_analysis.svg"), help="SVG summary plot path")
    parser.add_argument("--top-k", type=int, default=5, help="how many best runs to show in the report")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    rows = read_sweep_csv(args.csv)
    summary = summarize(rows)
    svg_line_plot(rows, args.plot_svg)
    write_markdown(rows, summary, args.report_md, args.top_k, args.plot_svg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
