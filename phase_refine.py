#!/usr/bin/env python3
"""Refine selected cells from a fast causet phase-diagram scan."""

from __future__ import annotations

import argparse
import contextlib
import csv
import io
import statistics
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, List, Sequence, Set, Tuple

from cones import ConesSimulator, generate_sprinkled_causet, transitive_reduction
from phase_diagram import (
    CUDAEnergyBackend,
    PhaseRun,
    RunTimeoutError,
    count_relations,
    cuda_available,
    log_progress,
    relation_density,
    run_timeout,
    summarize,
    write_heatmap,
    write_markdown,
    write_runs_csv,
    write_summary_csv,
)


@dataclass(frozen=True, order=True)
class RefineCell:
    n: int
    dim: int


def parse_cells(value: str) -> List[RefineCell]:
    cells: List[RefineCell] = []
    seen: Set[Tuple[int, int]] = set()
    for raw in value.split(","):
        item = raw.strip()
        if not item:
            continue
        if ":" not in item:
            raise argparse.ArgumentTypeError(f"cell must be formatted as n:dim, got {item!r}")
        n_raw, dim_raw = item.split(":", 1)
        n = int(n_raw)
        dim = int(dim_raw)
        if n <= 0 or dim <= 0:
            raise argparse.ArgumentTypeError("n and dim must be positive")
        key = (n, dim)
        if key not in seen:
            cells.append(RefineCell(n=n, dim=dim))
            seen.add(key)
    if not cells:
        raise argparse.ArgumentTypeError("expected at least one cell")
    return sorted(cells)


def select_cells_from_summary(args: argparse.Namespace) -> List[RefineCell]:
    if args.cells:
        return args.cells
    if args.source_summary_csv is None:
        raise ValueError("provide --cells or --source-summary-csv")

    cells: List[RefineCell] = []
    seen: Set[Tuple[int, int]] = set()
    with args.source_summary_csv.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            n = int(row["n"])
            dim = int(row["dim"])
            runs = int(row["runs"])
            success_rate = float(row["success_rate"])
            median_final = float(row["median_final"])
            min_final = float(row["min_final"])
            if args.max_n > 0 and n > args.max_n:
                continue
            if args.max_dim > 0 and dim > args.max_dim:
                continue
            if runs < args.min_source_runs:
                continue
            selected = (
                success_rate >= args.min_source_success_rate
                or median_final <= args.max_source_median
                or min_final <= args.max_source_min
            )
            if selected and (n, dim) not in seen:
                cells.append(RefineCell(n=n, dim=dim))
                seen.add((n, dim))
    return sorted(cells)


def write_outputs(args: argparse.Namespace, rows: Sequence[PhaseRun]) -> None:
    cells = summarize(rows, args.success_threshold) if rows else []
    write_runs_csv(rows, args.runs_csv)
    write_summary_csv(cells, args.summary_csv)
    write_heatmap(cells, args.heatmap_svg)
    write_markdown(cells, args.report_md, args.heatmap_svg, args.success_threshold)


def maybe_write_partial(args: argparse.Namespace, rows: Sequence[PhaseRun], scope: str) -> None:
    if args.partial == scope and rows:
        write_outputs(args, rows)
        log_progress(args, f"phase_refine: partial {scope} results written")


def run_refinement(args: argparse.Namespace) -> List[PhaseRun]:
    selected_cells = select_cells_from_summary(args)
    if not selected_cells:
        log_progress(args, "phase_refine: no cells selected")
        return []

    anneal_seeds = [args.seed_start + i * args.seed_step for i in range(args.seed_count)]
    total_runs = len(selected_cells) * len(anneal_seeds)
    rows: List[PhaseRun] = []
    run_index = 0
    started = time.monotonic()

    log_progress(
        args,
        (
            f"phase_refine: starting {total_runs} runs "
            f"({len(selected_cells)} cells x {len(anneal_seeds)} seeds), backend={args.backend}"
        ),
    )
    log_progress(args, "phase_refine: selected cells " + ", ".join(f"n={c.n} dim={c.dim}" for c in selected_cells))

    with tempfile.TemporaryDirectory(prefix="phase_refine_") as tmpdir:
        tmp = Path(tmpdir)
        for cell_index, cell in enumerate(selected_cells, start=1):
            cell_started = time.monotonic()
            generation_seed = args.generation_seed + cell.n * args.generation_seed_step
            z, _ = generate_sprinkled_causet(cell.n, generation_seed, args.sprinkle_spacetime_dim)
            relations = count_relations(z)
            links = len(transitive_reduction(z))
            density = relation_density(relations, cell.n)
            log_progress(
                args,
                (
                    f"phase_refine: cell {cell_index}/{len(selected_cells)} n={cell.n} dim={cell.dim} "
                    f"generated seed={generation_seed} relations={relations} links={links} density={density:.3f}"
                ),
            )

            shared_cuda_backend: Any | None = None
            try:
                if (
                    args.gpu_first
                    and args.backend in ("cuda", "auto")
                    and CUDAEnergyBackend is not None
                    and cuda_available()
                ):
                    shared_cuda_backend = CUDAEnergyBackend(cell.n, cell.dim)
                    shared_cuda_backend.set_z(z)
                    log_progress(args, f"phase_refine: cell {cell_index}/{len(selected_cells)} using shared CUDA backend")

                for seed in anneal_seeds:
                    run_index += 1
                    run_started = time.monotonic()
                    log_progress(
                        args,
                        (
                            f"phase_refine: run {run_index}/{total_runs} "
                            f"n={cell.n} dim={cell.dim} anneal_seed={seed} started"
                        ),
                    )
                    out_path = tmp / f"n{cell.n}_dim{cell.dim}_seed{seed}.out"
                    sim = ConesSimulator(
                        z=z,
                        dim=cell.dim,
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
                                f"phase_refine: run {run_index}/{total_runs} n={cell.n} dim={cell.dim} "
                                f"anneal_seed={seed} timed out after {args.max_run_seconds:.1f}s"
                            ),
                        )
                        if shared_cuda_backend is not None:
                            shared_cuda_backend.close()
                            shared_cuda_backend = None

                    rows.append(
                        PhaseRun(
                            n=cell.n,
                            dim=cell.dim,
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
                    log_progress(
                        args,
                        (
                            f"phase_refine: run {run_index}/{total_runs} n={cell.n} dim={cell.dim} "
                            f"anneal_seed={seed} final={final_energy:.6f} status={status} "
                            f"elapsed={time.monotonic() - run_started:.2f}s"
                        ),
                    )
                    maybe_write_partial(args, rows, "run")
            finally:
                if shared_cuda_backend is not None:
                    shared_cuda_backend.close()

            cell_rows = [row for row in rows if row.n == cell.n and row.dim == cell.dim]
            finals = [row.final_energy for row in cell_rows]
            success_rate = sum(1 for value in finals if value <= args.success_threshold) / len(finals)
            log_progress(
                args,
                (
                    f"phase_refine: cell {cell_index}/{len(selected_cells)} n={cell.n} dim={cell.dim} done "
                    f"success={success_rate:.2f} median={statistics.median(finals):.6f} "
                    f"elapsed={time.monotonic() - cell_started:.2f}s"
                ),
            )
            maybe_write_partial(args, rows, "cell")

    log_progress(args, f"phase_refine: completed in {time.monotonic() - started:.2f}s")
    return rows


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Deep-rerun selected cells from a causet phase diagram.")
    parser.add_argument("--source-summary-csv", type=Path, default=None, help="summary CSV produced by phase_diagram.py")
    parser.add_argument("--cells", type=parse_cells, default=None, help="explicit cells as n:dim,n:dim; overrides CSV selection")
    parser.add_argument("--min-source-runs", type=int, default=1)
    parser.add_argument("--min-source-success-rate", type=float, default=0.10)
    parser.add_argument("--max-source-median", type=float, default=0.05)
    parser.add_argument("--max-source-min", type=float, default=1e-6)
    parser.add_argument("--max-n", type=int, default=0, help="ignore source cells above this n; 0 disables")
    parser.add_argument("--max-dim", type=int, default=0, help="ignore source cells above this dim; 0 disables")
    parser.add_argument("--sprinkle-spacetime-dim", type=int, default=1)
    parser.add_argument("--generation-seed", type=int, default=1987)
    parser.add_argument("--generation-seed-step", type=int, default=17)
    parser.add_argument("--seed-start", type=int, default=1959)
    parser.add_argument("--seed-count", type=int, default=16)
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
    parser.add_argument("--max-run-seconds", type=float, default=180.0, help="timeout for one deep run; 0 disables it")
    parser.add_argument("--timeout-energy", type=float, default=1000000.0)
    parser.add_argument("--runs-csv", type=Path, default=Path("phase_refine_runs.csv"))
    parser.add_argument("--summary-csv", type=Path, default=Path("phase_refine_summary.csv"))
    parser.add_argument("--report-md", type=Path, default=Path("phase_refine.md"))
    parser.add_argument("--heatmap-svg", type=Path, default=Path("phase_refine.svg"))
    parser.add_argument("--quiet", action="store_true", help="suppress progress messages")
    parser.add_argument(
        "--partial",
        choices=["cell", "run", "none"],
        default="cell",
        help="write partial outputs after each completed cell or run",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.seed_count <= 0:
        parser.error("seed-count must be positive")
    if args.gpu_first and args.backend == "cpu":
        args.backend = "auto"

    try:
        rows = run_refinement(args)
    except ValueError as exc:
        parser.error(str(exc))
    except KeyboardInterrupt:
        log_progress(args, "phase_refine: interrupted; keeping last partial outputs on disk")
        raise

    cells = summarize(rows, args.success_threshold) if rows else []
    log_progress(args, f"phase_refine: writing {args.runs_csv}")
    write_runs_csv(rows, args.runs_csv)
    log_progress(args, f"phase_refine: writing {args.summary_csv}")
    write_summary_csv(cells, args.summary_csv)
    log_progress(args, f"phase_refine: writing {args.heatmap_svg}")
    write_heatmap(cells, args.heatmap_svg)
    log_progress(args, f"phase_refine: writing {args.report_md}")
    write_markdown(cells, args.report_md, args.heatmap_svg, args.success_threshold)
    log_progress(args, "phase_refine: done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
