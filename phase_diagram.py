#!/usr/bin/env python3
"""Generate n-vs-dim phase diagrams for the revived cones experiment."""

from __future__ import annotations

import argparse
import contextlib
import csv
import io
import signal
import statistics
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

from cones import ConesSimulator, generate_sprinkled_causet, transitive_reduction

try:
    from cuda_backend import CUDAEnergyBackend, cuda_available
except Exception:  # pragma: no cover - optional accelerator
    CUDAEnergyBackend = None

    def cuda_available() -> bool:  # type: ignore
        return False


@dataclass
class PhaseRun:
    n: int
    dim: int
    generation_seed: int
    anneal_seed: int
    relations: int
    links: int
    relation_density: float
    initial_energy: float
    warmup_energy: float
    final_energy: float
    points: int
    output_file: str
    status: str = "ok"


@dataclass
class PhaseCell:
    n: int
    dim: int
    runs: int
    success_rate: float
    mean_final: float
    median_final: float
    min_final: float
    max_final: float
    best_seed: int
    relations: int
    links: int
    relation_density: float


def parse_int_list(value: str) -> List[int]:
    items: List[int] = []
    for raw in value.split(","):
        raw = raw.strip()
        if not raw:
            continue
        items.append(int(raw))
    if not items:
        raise argparse.ArgumentTypeError("expected at least one integer")
    return items


def count_relations(z: Sequence[Sequence[bool]]) -> int:
    return sum(1 for i in range(len(z)) for j in range(i + 1, len(z)) if z[i][j])


def relation_density(relations: int, n: int) -> float:
    total = n * (n - 1) // 2
    return relations / total if total else 0.0


def log_progress(args: argparse.Namespace, message: str) -> None:
    if not args.quiet:
        print(message, flush=True)


class RunTimeoutError(TimeoutError):
    pass


@contextlib.contextmanager
def run_timeout(seconds: float):
    if seconds <= 0:
        yield
        return

    def handler(_signum, _frame):
        raise RunTimeoutError(f"run exceeded {seconds:.1f}s")

    old_handler = signal.getsignal(signal.SIGALRM)
    old_timer = signal.setitimer(signal.ITIMER_REAL, seconds)
    signal.signal(signal.SIGALRM, handler)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, old_timer[0], old_timer[1])
        signal.signal(signal.SIGALRM, old_handler)


def summarize(rows: Sequence[PhaseRun], success_threshold: float) -> List[PhaseCell]:
    grouped: Dict[Tuple[int, int], List[PhaseRun]] = {}
    for row in rows:
        grouped.setdefault((row.n, row.dim), []).append(row)

    cells: List[PhaseCell] = []
    for (n, dim), items in sorted(grouped.items()):
        finals = [row.final_energy for row in items]
        best = min(items, key=lambda row: row.final_energy)
        cells.append(
            PhaseCell(
                n=n,
                dim=dim,
                runs=len(items),
                success_rate=sum(1 for value in finals if value <= success_threshold) / len(finals),
                mean_final=statistics.fmean(finals),
                median_final=statistics.median(finals),
                min_final=min(finals),
                max_final=max(finals),
                best_seed=best.anneal_seed,
                relations=best.relations,
                links=best.links,
                relation_density=best.relation_density,
            )
        )
    return cells


def color_for_success(value: float) -> str:
    t = max(0.0, min(1.0, value))
    r = int(220 - 165 * t)
    g = int(70 + 120 * t)
    b = int(65 + 55 * t)
    return f"rgb({r},{g},{b})"


def write_runs_csv(rows: Sequence[PhaseRun], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "n",
                "dim",
                "generation_seed",
                "anneal_seed",
                "relations",
                "links",
                "relation_density",
                "initial_energy",
                "warmup_energy",
                "final_energy",
                "points",
                "output_file",
                "status",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row.n,
                    row.dim,
                    row.generation_seed,
                    row.anneal_seed,
                    row.relations,
                    row.links,
                    f"{row.relation_density:.6f}",
                    f"{row.initial_energy:.6f}",
                    f"{row.warmup_energy:.6f}",
                    f"{row.final_energy:.6f}",
                    row.points,
                    row.output_file,
                    row.status,
                ]
            )


def write_summary_csv(cells: Sequence[PhaseCell], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "n",
                "dim",
                "runs",
                "success_rate",
                "mean_final",
                "median_final",
                "min_final",
                "max_final",
                "best_seed",
                "relations",
                "links",
                "relation_density",
            ]
        )
        for cell in cells:
            writer.writerow(
                [
                    cell.n,
                    cell.dim,
                    cell.runs,
                    f"{cell.success_rate:.6f}",
                    f"{cell.mean_final:.6f}",
                    f"{cell.median_final:.6f}",
                    f"{cell.min_final:.6f}",
                    f"{cell.max_final:.6f}",
                    cell.best_seed,
                    cell.relations,
                    cell.links,
                    f"{cell.relation_density:.6f}",
                ]
            )


def write_heatmap(cells: Sequence[PhaseCell], path: Path) -> None:
    if not cells:
        path.write_text("<svg xmlns='http://www.w3.org/2000/svg' width='800' height='300'></svg>", encoding="utf-8")
        return

    ns = sorted({cell.n for cell in cells})
    dims = sorted({cell.dim for cell in cells})
    cell_w = 112
    cell_h = 78
    margin_x = 90
    margin_y = 62
    width = margin_x + cell_w * len(ns) + 40
    height = margin_y + cell_h * len(dims) + 70
    lookup = {(cell.n, cell.dim): cell for cell in cells}

    pieces: List[str] = []
    for xi, n in enumerate(ns):
        pieces.append(
            f"<text x='{margin_x + xi * cell_w + 52:.1f}' y='44' text-anchor='middle' font-family='monospace' font-size='13' fill='#222'>n={n}</text>"
        )
    for yi, dim in enumerate(dims):
        pieces.append(
            f"<text x='20' y='{margin_y + yi * cell_h + 32}' font-family='monospace' font-size='13' fill='#222'>dim {dim}</text>"
        )
        for xi, n in enumerate(ns):
            x = margin_x + xi * cell_w
            y = margin_y + yi * cell_h
            cell = lookup.get((n, dim))
            if cell is None:
                pieces.append(
                    f"<rect x='{x}' y='{y}' width='{cell_w - 8}' height='{cell_h - 8}' rx='8' fill='#e4ded4' />"
                )
                pieces.append(
                    f"<text x='{x + 52:.1f}' y='{y + 40}' text-anchor='middle' font-family='monospace' font-size='12' fill='#777'>pending</text>"
                )
                continue
            pieces.append(
                f"<rect x='{x}' y='{y}' width='{cell_w - 8}' height='{cell_h - 8}' rx='8' fill='{color_for_success(cell.success_rate)}' />"
            )
            pieces.append(
                f"<text x='{x + 52:.1f}' y='{y + 28}' text-anchor='middle' font-family='monospace' font-size='13' fill='#111'>p={cell.success_rate:.2f}</text>"
            )
            pieces.append(
                f"<text x='{x + 52:.1f}' y='{y + 50}' text-anchor='middle' font-family='monospace' font-size='11' fill='#111'>m={cell.median_final:.3f}</text>"
            )

    svg = f"""<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>
  <rect width='100%' height='100%' fill='#f7f3eb'/>
  <text x='{width/2:.0f}' y='28' text-anchor='middle' font-family='serif' font-size='22' fill='#222'>Causet embeddability phase diagram</text>
  {''.join(pieces)}
  <text x='{width - 260}' y='{height - 22}' font-family='monospace' font-size='12' fill='#444'>green = higher success probability</text>
</svg>
"""
    path.write_text(svg, encoding="utf-8")


def write_markdown(cells: Sequence[PhaseCell], path: Path, heatmap: Path, success_threshold: float) -> None:
    lines = ["# Phase Diagram", ""]
    lines.append(f"- Heatmap: `{heatmap}`")
    lines.append(f"- Success threshold: `{success_threshold}`")
    lines.append("")
    if not cells:
        lines.append("No cells found.")
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return

    best = max(cells, key=lambda cell: (cell.success_rate, -cell.median_final))
    lines.append("## Summary")
    lines.append("")
    lines.append(
        f"- Best cell: n={best.n}, dim={best.dim}, success={best.success_rate:.2f}, median_final={best.median_final:.6f}"
    )
    lines.append("")
    lines.append("## Cells")
    lines.append("")
    lines.append("| n | dim | runs | success_rate | median_final | mean_final | min_final | max_final | relations | links | density | best_seed |")
    lines.append("| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for cell in sorted(cells, key=lambda c: (c.n, c.dim)):
        lines.append(
            f"| {cell.n} | {cell.dim} | {cell.runs} | {cell.success_rate:.2f} | {cell.median_final:.6f} | {cell.mean_final:.6f} | {cell.min_final:.6f} | {cell.max_final:.6f} | {cell.relations} | {cell.links} | {cell.relation_density:.3f} | {cell.best_seed} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_outputs(args: argparse.Namespace, rows: Sequence[PhaseRun]) -> None:
    cells = summarize(rows, args.success_threshold) if rows else []
    write_runs_csv(rows, args.runs_csv)
    write_summary_csv(cells, args.summary_csv)
    write_heatmap(cells, args.heatmap_svg)
    write_markdown(cells, args.report_md, args.heatmap_svg, args.success_threshold)


def maybe_write_partial(args: argparse.Namespace, rows: Sequence[PhaseRun], scope: str) -> None:
    if args.partial == scope and rows:
        write_outputs(args, rows)
        log_progress(args, f"phase_diagram: partial {scope} results written")


def cell_stop_reason(args: argparse.Namespace, cell_rows: Sequence[PhaseRun]) -> str:
    if not cell_rows:
        return ""
    if len(cell_rows) < args.min_cell_runs:
        return ""

    timeouts = sum(1 for row in cell_rows if row.status == "timeout")
    if args.stop_cell_timeouts > 0 and timeouts >= args.stop_cell_timeouts:
        return f"{timeouts} timeouts"

    timeout_rate = timeouts / len(cell_rows)
    if args.stop_cell_timeout_rate <= 1.0 and timeout_rate >= args.stop_cell_timeout_rate:
        return f"timeout_rate={timeout_rate:.2f}"

    if args.stop_cell_median_above > 0.0:
        finals = [row.final_energy for row in cell_rows]
        median = statistics.median(finals)
        successes = sum(1 for value in finals if value <= args.success_threshold)
        if successes == 0 and median >= args.stop_cell_median_above:
            return f"median={median:.6f}"

    return ""


def run_phase_diagram(args: argparse.Namespace) -> Tuple[List[PhaseRun], List[PhaseCell]]:
    rows: List[PhaseRun] = []
    dims = list(range(args.dim_min, args.dim_max + 1))
    anneal_seeds = [args.seed_start + i * args.seed_step for i in range(args.seed_count)]
    total_cells = len(args.n_values) * len(dims)
    total_runs = total_cells * len(anneal_seeds)
    cell_index = 0
    run_index = 0
    started = time.monotonic()

    log_progress(
        args,
        (
            f"phase_diagram: starting {total_runs} runs "
            f"({len(args.n_values)} n-values x {len(dims)} dims x {len(anneal_seeds)} seeds), "
            f"backend={args.backend}"
        ),
    )

    with tempfile.TemporaryDirectory(prefix="phase_diagram_") as tmpdir:
        tmp = Path(tmpdir)
        for n in args.n_values:
            generation_seed = args.generation_seed + n * args.generation_seed_step
            z, _ = generate_sprinkled_causet(n, generation_seed, args.sprinkle_spacetime_dim)
            relations = count_relations(z)
            links = len(transitive_reduction(z))
            density = relation_density(relations, n)
            log_progress(
                args,
                (
                    f"phase_diagram: generated n={n} seed={generation_seed} "
                    f"relations={relations} links={links} density={density:.3f}"
                ),
            )

            for dim in dims:
                cell_index += 1
                cell_started = time.monotonic()
                log_progress(args, f"phase_diagram: cell {cell_index}/{total_cells} n={n} dim={dim} started")
                shared_cuda_backend: Any | None = None
                try:
                    if (
                        args.gpu_first
                        and args.backend in ("cuda", "auto")
                        and CUDAEnergyBackend is not None
                        and cuda_available()
                    ):
                        shared_cuda_backend = CUDAEnergyBackend(n, dim)
                        shared_cuda_backend.set_z(z)
                        log_progress(args, f"phase_diagram: cell {cell_index}/{total_cells} using shared CUDA backend")

                    for seed in anneal_seeds:
                        current_cell_rows = [row for row in rows if row.n == n and row.dim == dim]
                        reason = cell_stop_reason(args, current_cell_rows)
                        if reason:
                            log_progress(
                                args,
                                (
                                    f"phase_diagram: cell {cell_index}/{total_cells} n={n} dim={dim} "
                                    f"stopped early after {len(current_cell_rows)} runs ({reason})"
                                ),
                            )
                            break

                        run_index += 1
                        run_started = time.monotonic()
                        log_progress(
                            args,
                            (
                                f"phase_diagram: run {run_index}/{total_runs} "
                                f"n={n} dim={dim} anneal_seed={seed} started"
                            ),
                        )
                        out_path = tmp / f"n{n}_dim{dim}_seed{seed}.out"
                        sim = ConesSimulator(
                            z=z,
                            dim=dim,
                            seed=seed,
                            interactive=False,
                            max_data=args.max_data,
                            plot_path=None,
                            backend=args.backend,
                            shared_cuda_backend=shared_cuda_backend,
                            warmup_limit=args.warmup_limit,
                            anneal_limit=args.anneal_limit,
                            initial_temp=args.initial_temp,
                            cooling_factor=args.cooling_factor,
                            acceptance_scale=args.acceptance_scale,
                        )
                        status = "ok"
                        try:
                            with run_timeout(args.max_run_seconds), contextlib.redirect_stdout(io.StringIO()):
                                sim.run(out_path)
                            final_energy = sim.data[-1][1] if sim.data else sim.eave
                            initial_energy = sim.initial_energy
                            warmup_energy = sim.warmup_energy
                            points = len(sim.data)
                        except RunTimeoutError:
                            status = "timeout"
                            final_energy = args.timeout_energy
                            initial_energy = sim.initial_energy
                            warmup_energy = sim.warmup_energy
                            points = len(sim.data)
                            log_progress(
                                args,
                                (
                                    f"phase_diagram: run {run_index}/{total_runs} "
                                    f"n={n} dim={dim} anneal_seed={seed} timed out after "
                                    f"{args.max_run_seconds:.1f}s"
                                ),
                            )
                            if shared_cuda_backend is not None:
                                shared_cuda_backend.close()
                                shared_cuda_backend = None
                        log_progress(
                            args,
                            (
                                f"phase_diagram: run {run_index}/{total_runs} "
                                f"n={n} dim={dim} anneal_seed={seed} final={final_energy:.6f} status={status} "
                                f"elapsed={time.monotonic() - run_started:.2f}s"
                            ),
                        )
                        rows.append(
                            PhaseRun(
                                n=n,
                                dim=dim,
                                generation_seed=generation_seed,
                                anneal_seed=seed,
                                relations=relations,
                                links=links,
                                relation_density=density,
                                initial_energy=initial_energy,
                                warmup_energy=warmup_energy,
                                final_energy=final_energy,
                                points=points,
                                output_file=str(out_path),
                                status=status,
                            )
                        )
                        maybe_write_partial(args, rows, "run")
                finally:
                    if shared_cuda_backend is not None:
                        shared_cuda_backend.close()
                cell_rows = [row for row in rows if row.n == n and row.dim == dim]
                if cell_rows:
                    finals = [row.final_energy for row in cell_rows]
                    success_rate = sum(1 for value in finals if value <= args.success_threshold) / len(finals)
                    log_progress(
                        args,
                        (
                            f"phase_diagram: cell {cell_index}/{total_cells} n={n} dim={dim} done "
                            f"success={success_rate:.2f} median={statistics.median(finals):.6f} "
                            f"elapsed={time.monotonic() - cell_started:.2f}s"
                        ),
                    )
                maybe_write_partial(args, rows, "cell")

    cells = summarize(rows, args.success_threshold)
    log_progress(args, f"phase_diagram: completed in {time.monotonic() - started:.2f}s")
    return rows, cells


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build a causet embeddability phase diagram.")
    parser.add_argument("--n-values", type=parse_int_list, default=[6, 12], help="comma-separated causet sizes")
    parser.add_argument("--dim-min", type=int, default=1)
    parser.add_argument("--dim-max", type=int, default=4)
    parser.add_argument("--sprinkle-spacetime-dim", type=int, default=1)
    parser.add_argument("--generation-seed", type=int, default=1987)
    parser.add_argument("--generation-seed-step", type=int, default=17)
    parser.add_argument("--seed-start", type=int, default=1959)
    parser.add_argument("--seed-count", type=int, default=4)
    parser.add_argument("--seed-step", type=int, default=1)
    parser.add_argument("--backend", choices=["cpu", "cuda", "auto"], default="cpu")
    parser.add_argument("--gpu-first", action="store_true")
    parser.add_argument("--warmup-limit", type=int, default=100)
    parser.add_argument("--anneal-limit", type=int, default=100)
    parser.add_argument("--initial-temp", type=float, default=180.0)
    parser.add_argument("--cooling-factor", type=float, default=0.8)
    parser.add_argument("--acceptance-scale", type=float, default=4.0)
    parser.add_argument("--max-data", type=int, default=35)
    parser.add_argument("--success-threshold", type=float, default=1e-6)
    parser.add_argument("--runs-csv", type=Path, default=Path("phase_diagram_runs.csv"))
    parser.add_argument("--summary-csv", type=Path, default=Path("phase_diagram_summary.csv"))
    parser.add_argument("--report-md", type=Path, default=Path("phase_diagram.md"))
    parser.add_argument("--heatmap-svg", type=Path, default=Path("phase_diagram.svg"))
    parser.add_argument("--quiet", action="store_true", help="suppress progress messages")
    parser.add_argument(
        "--partial",
        choices=["cell", "run", "none"],
        default="cell",
        help="write partial CSV/Markdown/SVG outputs after each completed cell or run",
    )
    parser.add_argument("--max-run-seconds", type=float, default=0.0, help="timeout for a single annealing run; 0 disables it")
    parser.add_argument("--timeout-energy", type=float, default=1000000.0, help="final energy assigned to timed-out runs")
    parser.add_argument("--fast-frontier", action="store_true", help="apply conservative settings for large frontier scans")
    parser.add_argument("--min-cell-runs", type=int, default=0, help="minimum runs before early cell stopping is allowed")
    parser.add_argument("--stop-cell-timeouts", type=int, default=0, help="stop a cell after this many timed-out runs; 0 disables it")
    parser.add_argument(
        "--stop-cell-timeout-rate",
        type=float,
        default=2.0,
        help="stop a cell when timeout rate reaches this value after min-cell-runs; >1 disables it",
    )
    parser.add_argument(
        "--stop-cell-median-above",
        type=float,
        default=0.0,
        help="stop a cell with no successes when median final energy reaches this value; 0 disables it",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.dim_min <= 0 or args.dim_max < args.dim_min:
        parser.error("invalid dimension range")
    if args.seed_count <= 0:
        parser.error("seed-count must be positive")
    if args.gpu_first and args.backend == "cpu":
        args.backend = "auto"
    if args.fast_frontier:
        if args.warmup_limit == 100:
            args.warmup_limit = 40
        if args.anneal_limit == 100:
            args.anneal_limit = 40
        if args.max_data == 35:
            args.max_data = 12
        if args.max_run_seconds == 0.0:
            args.max_run_seconds = 45.0
        if args.min_cell_runs == 0:
            args.min_cell_runs = 3
        if args.stop_cell_timeouts == 0:
            args.stop_cell_timeouts = 2
        if args.stop_cell_timeout_rate > 1.0:
            args.stop_cell_timeout_rate = 0.67
        if args.stop_cell_median_above == 0.0:
            args.stop_cell_median_above = 100.0

    try:
        rows, cells = run_phase_diagram(args)
    except KeyboardInterrupt:
        log_progress(args, "phase_diagram: interrupted; keeping last partial outputs on disk")
        raise
    log_progress(args, f"phase_diagram: writing {args.runs_csv}")
    write_runs_csv(rows, args.runs_csv)
    log_progress(args, f"phase_diagram: writing {args.summary_csv}")
    write_summary_csv(cells, args.summary_csv)
    log_progress(args, f"phase_diagram: writing {args.heatmap_svg}")
    write_heatmap(cells, args.heatmap_svg)
    log_progress(args, f"phase_diagram: writing {args.report_md}")
    write_markdown(cells, args.report_md, args.heatmap_svg, args.success_threshold)
    log_progress(args, "phase_diagram: done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
