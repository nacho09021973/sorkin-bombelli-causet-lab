#!/usr/bin/env python3
"""Sweep dimensions and seeds for a single causet benchmark.

This is the next layer above the per-seed sweep:
- compare `dim` values directly
- find which dimension gives the lowest energies for a given input
- produce a compact markdown summary and a heatmap SVG
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
from typing import Dict, List, Sequence

from cones import ConesSimulator, parse_cones_input


@dataclass
class DimRow:
    dim: int
    seed: int
    initial_energy: float
    final_energy: float
    warmup_energy: float
    points: int
    output_file: str


def read_rows(path: Path) -> List[DimRow]:
    rows: List[DimRow] = []
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for raw in reader:
            rows.append(
                DimRow(
                    dim=int(raw["dim"]),
                    seed=int(raw["seed"]),
                    initial_energy=float(raw["initial_energy"]),
                    final_energy=float(raw["final_energy"]),
                    warmup_energy=float(raw["warmup_energy"]),
                    points=int(raw["points"]),
                    output_file=raw["output_file"],
                )
            )
    return rows


def color_for(value: float, low: float, high: float) -> str:
    if high <= low:
        return "#2a9d8f"
    t = max(0.0, min(1.0, (value - low) / (high - low)))
    # Green for low, red for high.
    r = int(45 + 190 * t)
    g = int(170 - 120 * t)
    b = int(120 - 60 * t)
    return f"rgb({r},{g},{b})"


def summarize(rows: Sequence[DimRow]) -> Dict[int, Dict[str, float]]:
    grouped: Dict[int, List[DimRow]] = {}
    for row in rows:
        grouped.setdefault(row.dim, []).append(row)

    summary: Dict[int, Dict[str, float]] = {}
    for dim, items in grouped.items():
        finals = [r.final_energy for r in items]
        summary[dim] = {
            "count": float(len(items)),
            "mean_final": statistics.fmean(finals),
            "median_final": statistics.median(finals),
            "min_final": min(finals),
            "max_final": max(finals),
            "best_seed": float(min(items, key=lambda r: r.final_energy).seed),
        }
    return summary


def write_heatmap(rows: Sequence[DimRow], path: Path) -> None:
    if not rows:
        path.write_text(
            "<svg xmlns='http://www.w3.org/2000/svg' width='800' height='300'></svg>",
            encoding="utf-8",
        )
        return

    dims = sorted({r.dim for r in rows})
    seeds = sorted({r.seed for r in rows})
    finals = [r.final_energy for r in rows]
    low, high = min(finals), max(finals)

    width = 120 + 90 * len(seeds)
    height = 120 + 70 * len(dims)
    cell_w = 90
    cell_h = 70
    margin_x = 80
    margin_y = 60
    cells = []
    lookup = {(r.dim, r.seed): r for r in rows}
    for yi, dim in enumerate(dims):
        for xi, seed in enumerate(seeds):
            r = lookup[(dim, seed)]
            fill = color_for(r.final_energy, low, high)
            x = margin_x + xi * cell_w
            y = margin_y + yi * cell_h
            cells.append(
                f"<rect x='{x}' y='{y}' width='{cell_w - 6}' height='{cell_h - 6}' rx='8' fill='{fill}' />"
            )
            cells.append(
                f"<text x='{x + (cell_w - 6)/2:.1f}' y='{y + 28}' text-anchor='middle' font-family='monospace' font-size='13' fill='#111'>{r.final_energy:.3f}</text>"
            )
            cells.append(
                f"<text x='{x + (cell_w - 6)/2:.1f}' y='{y + 50}' text-anchor='middle' font-family='monospace' font-size='11' fill='#111'>s{seed}</text>"
            )
        cells.append(
            f"<text x='20' y='{margin_y + yi*cell_h + 30}' font-family='monospace' font-size='13' fill='#222'>dim {dim}</text>"
        )

    seed_labels = []
    for xi, seed in enumerate(seeds):
        seed_labels.append(
            f"<text x='{margin_x + xi*cell_w + (cell_w - 6)/2:.1f}' y='{margin_y - 18}' text-anchor='middle' font-family='monospace' font-size='13' fill='#222'>{seed}</text>"
        )

    svg = f"""<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>
  <rect width='100%' height='100%' fill='#f7f3eb'/>
  <text x='{width/2:.0f}' y='32' text-anchor='middle' font-family='serif' font-size='22' fill='#222'>Dimension sweep final energy</text>
  <text x='{width/2:.0f}' y='{height - 18}' text-anchor='middle' font-family='monospace' font-size='16' fill='#222'>seed</text>
  {''.join(seed_labels)}
  {''.join(cells)}
  <text x='{width - 180}' y='40' font-family='monospace' font-size='12' fill='#444'>green = lower energy</text>
</svg>
"""
    path.write_text(svg, encoding="utf-8")


def write_markdown(rows: Sequence[DimRow], summary: Dict[int, Dict[str, float]], path: Path, heatmap: Path) -> None:
    lines = ["# Dimension Sweep", ""]
    if not rows:
        lines.append("No rows found.")
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return

    lines.append(f"- Heatmap: `{heatmap}`")
    lines.append("")
    lines.append("## Per-dimension summary")
    lines.append("")
    lines.append("| dim | runs | mean_final | median_final | min_final | max_final | best_seed |")
    lines.append("| ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for dim in sorted(summary):
        s = summary[dim]
        lines.append(
            f"| {dim} | {int(s['count'])} | {s['mean_final']:.6f} | {s['median_final']:.6f} | {s['min_final']:.6f} | {s['max_final']:.6f} | {int(s['best_seed'])} |"
        )

    lines.append("")
    lines.append("## All runs")
    lines.append("")
    lines.append("| dim | seed | final_energy | warmup_energy | points | output_file |")
    lines.append("| ---: | ---: | ---: | ---: | ---: | --- |")
    for row in sorted(rows, key=lambda r: (r.dim, r.final_energy, r.seed)):
        lines.append(
            f"| {row.dim} | {row.seed} | {row.final_energy:.6f} | {row.warmup_energy:.6f} | {row.points} | {row.output_file} |"
        )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Sweep dimensions for a benchmark causet.")
    parser.add_argument("input_file", type=Path, help="benchmark incidence file")
    parser.add_argument("--dim-min", type=int, default=1)
    parser.add_argument("--dim-max", type=int, default=4)
    parser.add_argument("--seed-start", type=int, default=1959)
    parser.add_argument("--seed-count", type=int, default=12)
    parser.add_argument("--seed-step", type=int, default=1)
    parser.add_argument("--output-csv", type=Path, default=Path("dimension_sweep.csv"))
    parser.add_argument("--report-md", type=Path, default=Path("dimension_sweep.md"))
    parser.add_argument("--heatmap-svg", type=Path, default=Path("dimension_sweep.svg"))
    parser.add_argument("--backend", choices=["cpu", "cuda", "auto"], default="cpu")
    parser.add_argument("--max-data", type=int, default=35)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    z = parse_cones_input(args.input_file)
    seeds = [args.seed_start + i * args.seed_step for i in range(args.seed_count)]
    rows: List[DimRow] = []

    with tempfile.TemporaryDirectory(prefix="dimension_sweep_") as tmpdir:
        tmp = Path(tmpdir)
        for dim in range(args.dim_min, args.dim_max + 1):
            for seed in seeds:
                out_path = tmp / f"dim{dim}_seed{seed}.out"
                sim = ConesSimulator(
                    z=z,
                    dim=dim,
                    seed=seed,
                    interactive=False,
                    max_data=args.max_data,
                    plot_path=None,
                    backend=args.backend,
                )
                with contextlib.redirect_stdout(io.StringIO()):
                    sim.run(out_path)
                rows.append(
                    DimRow(
                        dim=dim,
                        seed=seed,
                        initial_energy=sim.initial_energy,
                        final_energy=sim.data[-1][1] if sim.data else sim.eave,
                        warmup_energy=sim.warmup_energy,
                        points=len(sim.data),
                        output_file=str(out_path),
                    )
                )

    rows_csv = args.output_csv
    with rows_csv.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["dim", "seed", "initial_energy", "final_energy", "warmup_energy", "points", "output_file"])
        for row in rows:
            writer.writerow(
                [
                    row.dim,
                    row.seed,
                    f"{row.initial_energy:.6f}",
                    f"{row.final_energy:.6f}",
                    f"{row.warmup_energy:.6f}",
                    row.points,
                    row.output_file,
                ]
            )

    summary = summarize(rows)
    write_heatmap(rows, args.heatmap_svg)
    write_markdown(rows, summary, args.report_md, args.heatmap_svg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
