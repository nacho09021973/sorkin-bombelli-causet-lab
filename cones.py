#!/usr/bin/env python3
"""Reimplementation of the Pascal thesis program ``cones``.

The original listing appears in ``Pascal.pdf`` as appendix A.2,
``An application of simulated annealing``.

This module keeps the structure of the Pascal program but exposes a
small command-line interface so the simulation can be rerun on modern
systems.
"""

from __future__ import annotations

import argparse
import math
import random
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Sequence, TextIO, Tuple

try:
    from cuda_backend import CUDAEnergyBackend
except Exception:  # pragma: no cover - optional accelerator
    CUDAEnergyBackend = None


def parse_cones_input(path: Path) -> List[List[bool]]:
    """Parse the Pascal input file format.

    The file contains N, followed by the upper-triangular entries of the
    incidence matrix in row-major order: z[i][j] for 1 <= i < j <= N.
    """

    tokens = path.read_text().split()
    if not tokens:
        raise ValueError(f"{path} is empty")

    values = [int(tok) for tok in tokens]
    n = values[0]
    if n <= 0:
        raise ValueError("N must be positive")

    expected = 1 + n * (n - 1) // 2
    if len(values) < expected:
        raise ValueError(
            f"expected at least {expected} integers, found {len(values)}"
        )

    z = [[False for _ in range(n)] for _ in range(n)]
    idx = 1
    for i in range(n - 1):
        for j in range(i + 1, n):
            zij = values[idx]
            idx += 1
            if zij not in (0, 1):
                raise ValueError(f"invalid matrix entry {zij} at position {(i, j)}")
            z[i][j] = bool(zij)
    return z


def format_matrix_bool(matrix: Sequence[Sequence[bool]]) -> str:
    lines = []
    for row in matrix:
        lines.append(" ".join("1" if cell else "0" for cell in row))
    return "\n".join(lines)


def format_matrix_float(matrix: Sequence[Sequence[float]]) -> str:
    lines = []
    for row in matrix:
        lines.append(" ".join(f"{val:12.6f}" for val in row))
    return "\n".join(lines)


@dataclass
class PascalRNG:
    """Port of the thesis ``ran2`` and ``gasdev`` routines.

    The generator keeps the original shuffle table and Box-Muller cache.
    It initializes on the first use, and also honors a negative seed for
    explicit reinitialization, matching the Pascal routine's convention.
    """

    seed: int
    m: int = 714025
    ia: int = 1366
    ic: int = 150889
    rm: float = 1.400512e-6
    gliy: int = 0
    glir: List[int] = field(default_factory=lambda: [0] * 97)
    gliset: int = 0
    glgset: float = 0.0
    initialized: bool = False

    def _initialize(self) -> None:
        idum = (-self.seed if self.seed > 0 else self.seed) or -1
        idum = (self.ic - idum) % self.m
        for j in range(97):
            idum = (self.ia * idum + self.ic) % self.m
            self.glir[j] = idum
        idum = (self.ia * idum + self.ic) % self.m
        self.gliy = idum
        self.seed = idum
        self.initialized = True

    def ran2(self) -> float:
        if self.seed < 0 or not self.initialized:
            self._initialize()
        j = 1 + (97 * self.gliy) // self.m
        if j < 1 or j > 97:
            raise RuntimeError("pause in routine RAN2")
        self.gliy = self.glir[j - 1]
        value = self.gliy * self.rm
        self.seed = (self.ia * self.seed + self.ic) % self.m
        self.glir[j - 1] = self.seed
        return value

    def gasdev(self) -> float:
        if self.gliset == 0:
            while True:
                v1 = 2.0 * self.ran2() - 1.0
                v2 = 2.0 * self.ran2() - 1.0
                r = v1 * v1 + v2 * v2
                if r < 1.0 and r > 0.0:
                    break
            fac = math.sqrt(-2.0 * math.log(r) / r)
            self.glgset = v1 * fac
            self.gliset = 1
            return v2 * fac
        self.gliset = 0
        return self.glgset


def write_svg_plot(points: Sequence[Tuple[float, float]], path: Path) -> None:
    """Write a simple SVG scatter/line plot for temperature vs average energy."""

    width = 900
    height = 600
    margin = 70
    x_vals = [p[0] for p in points]
    y_vals = [p[1] for p in points]
    if not x_vals:
        path.write_text(
            "<svg xmlns='http://www.w3.org/2000/svg' width='900' height='600'>"
            "<text x='40' y='50' font-family='monospace' font-size='20'>No data</text>"
            "</svg>",
            encoding="utf-8",
        )
        return

    x_min = min(x_vals)
    x_max = max(x_vals)
    y_min = min(y_vals)
    y_max = max(y_vals)
    if x_min == x_max:
        x_min -= 1.0
        x_max += 1.0
    if y_min == y_max:
        y_min -= 1.0
        y_max += 1.0

    def sx(x: float) -> float:
        return margin + (x - x_min) * (width - 2 * margin) / (x_max - x_min)

    def sy(y: float) -> float:
        return height - margin - (y - y_min) * (height - 2 * margin) / (y_max - y_min)

    polyline = " ".join(f"{sx(x):.2f},{sy(y):.2f}" for x, y in points)
    circles = "\n".join(
        f"<circle cx='{sx(x):.2f}' cy='{sy(y):.2f}' r='3.5' fill='#e4572e' />"
        for x, y in points
    )
    svg = f"""<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>
  <rect width='100%' height='100%' fill='#f7f4ed'/>
  <line x1='{margin}' y1='{height - margin}' x2='{width - margin}' y2='{height - margin}' stroke='#333' stroke-width='2'/>
  <line x1='{margin}' y1='{margin}' x2='{margin}' y2='{height - margin}' stroke='#333' stroke-width='2'/>
  <text x='{width/2:.0f}' y='34' text-anchor='middle' font-family='serif' font-size='22' fill='#222'>T vs Eave</text>
  <text x='{width/2:.0f}' y='{height - 18}' text-anchor='middle' font-family='monospace' font-size='16' fill='#222'>T</text>
  <text x='20' y='{height/2:.0f}' text-anchor='middle' font-family='monospace' font-size='16' fill='#222' transform='rotate(-90 20 {height/2:.0f})'>Eave</text>
  <polyline fill='none' stroke='#1f77b4' stroke-width='2.5' points='{polyline}' />
  {circles}
  <text x='{margin}' y='{margin - 18}' font-family='monospace' font-size='14' fill='#444'>n={len(points)}</text>
</svg>
"""
    path.write_text(svg, encoding="utf-8")


def generate_sprinkled_causet(
    n: int,
    seed: int,
    spacetime_dim: int = 1,
) -> Tuple[List[List[bool]], List[Tuple[float, ...]]]:
    """Generate a random causet by sprinkling in a unit causal diamond.

    The construction is 1+spacetime_dim dimensional Minkowski space.
    Events are uniformly distributed in the causal diamond defined by
    0 <= u <= 1 and 0 <= v <= 1 with light-cone coordinates
    t = (u + v)/2, x = (u - v)/2 for the 1+1 case.

    For higher dimensions we use a simple rejection sampler inside the
    unit ball in the spatial directions. The causal relation is defined
    by Minkowski timelike or null separation.
    """

    if n <= 0:
        raise ValueError("sprinkle count must be positive")
    if spacetime_dim <= 0:
        raise ValueError("spacetime_dim must be positive")

    rng = random.Random(seed)

    points: List[Tuple[float, ...]] = []
    if spacetime_dim == 1:
        for _ in range(n):
            u = rng.random()
            v = rng.random()
            t = 0.5 * (u + v)
            x = 0.5 * (u - v)
            points.append((t, x))
    else:
        for _ in range(n):
            t = rng.random()
            while True:
                coords = [2.0 * rng.random() - 1.0 for _ in range(spacetime_dim)]
                norm_sq = sum(c * c for c in coords)
                if norm_sq <= 1.0:
                    break
            # Shrink spatial radius so the point remains inside the diamond.
            scale = min(t, 1.0 - t)
            spatial = tuple(scale * c for c in coords)
            points.append((t, *spatial))

    # Sort by time to get a stable labelling.
    points = sorted(points, key=lambda p: (p[0],) + p[1:])

    matrix = [[False for _ in range(n)] for _ in range(n)]
    for i in range(n - 1):
        pi = points[i]
        for j in range(i + 1, n):
            pj = points[j]
            dt = pj[0] - pi[0]
            if spacetime_dim == 1:
                dx = pj[1] - pi[1]
                if dt >= abs(dx):
                    matrix[i][j] = True
            else:
                spatial_sq = 0.0
                for k in range(1, spacetime_dim + 1):
                    diff = pj[k] - pi[k]
                    spatial_sq += diff * diff
                if dt * dt >= spatial_sq:
                    matrix[i][j] = True

    return matrix, points


def transitive_reduction(matrix: Sequence[Sequence[bool]]) -> List[Tuple[int, int]]:
    """Return the covering relations of an upper-triangular causal matrix."""

    n = len(matrix)
    links: List[Tuple[int, int]] = []
    for i in range(n - 1):
        for j in range(i + 1, n):
            if not matrix[i][j]:
                continue
            covered = False
            for k in range(i + 1, j):
                if matrix[i][k] and matrix[k][j]:
                    covered = True
                    break
            if not covered:
                links.append((i, j))
    return links


def write_causet_svg(
    points: Sequence[Tuple[float, ...]],
    path: Path,
    matrix: Sequence[Sequence[bool]] | None = None,
) -> None:
    """Plot a sprinkled causet in a simple t-x plane SVG."""

    width = 900
    height = 600
    margin = 70
    if not points:
        path.write_text(
            "<svg xmlns='http://www.w3.org/2000/svg' width='900' height='600'>"
            "<text x='40' y='50' font-family='monospace' font-size='20'>No points</text>"
            "</svg>",
            encoding="utf-8",
        )
        return

    # For display, use time on x-axis and the first spatial coordinate on y-axis.
    xs = [p[0] for p in points]
    ys = [p[1] if len(p) > 1 else 0.0 for p in points]
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    if x_min == x_max:
        x_min -= 1.0
        x_max += 1.0
    if y_min == y_max:
        y_min -= 1.0
        y_max += 1.0

    def sx(x: float) -> float:
        return margin + (x - x_min) * (width - 2 * margin) / (x_max - x_min)

    def sy(y: float) -> float:
        return height - margin - (y - y_min) * (height - 2 * margin) / (y_max - y_min)

    edges = ""
    if matrix is not None:
        edge_items = []
        for i, j in transitive_reduction(matrix):
            xi = sx(points[i][0])
            yi = sy(points[i][1] if len(points[i]) > 1 else 0.0)
            xj = sx(points[j][0])
            yj = sy(points[j][1] if len(points[j]) > 1 else 0.0)
            edge_items.append(
                f"<line x1='{xi:.2f}' y1='{yi:.2f}' x2='{xj:.2f}' y2='{yj:.2f}' "
                f"stroke='#9b9b9b' stroke-width='1.2' opacity='0.75' />"
            )
        edges = "\n  ".join(edge_items)

    circles = "\n".join(
        f"<circle cx='{sx(p[0]):.2f}' cy='{sy(p[1] if len(p) > 1 else 0.0):.2f}' "
        f"r='4' fill='#2a9d8f' />"
        for p in points
    )
    svg = f"""<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>
  <rect width='100%' height='100%' fill='#f4f1ea'/>
  <line x1='{margin}' y1='{height - margin}' x2='{width - margin}' y2='{height - margin}' stroke='#333' stroke-width='2'/>
  <line x1='{margin}' y1='{margin}' x2='{margin}' y2='{height - margin}' stroke='#333' stroke-width='2'/>
  <text x='{width/2:.0f}' y='34' text-anchor='middle' font-family='serif' font-size='22' fill='#222'>Sprinkled causet</text>
  <text x='{width/2:.0f}' y='{height - 18}' text-anchor='middle' font-family='monospace' font-size='16' fill='#222'>t</text>
  <text x='20' y='{height/2:.0f}' text-anchor='middle' font-family='monospace' font-size='16' fill='#222' transform='rotate(-90 20 {height/2:.0f})'>x</text>
  {edges}
  {circles}
  <text x='{margin}' y='{margin - 18}' font-family='monospace' font-size='14' fill='#444'>n={len(points)}</text>
</svg>
"""
    path.write_text(svg, encoding="utf-8")


def write_sweep_report(results: Sequence[Dict[str, Any]], path: Path) -> None:
    """Write a compact markdown report for a batch of runs."""

    if not results:
        path.write_text("# Sweep report\n\nNo runs.\n", encoding="utf-8")
        return

    lines = ["# Sweep report", ""]
    lines.append(
        "| seed | n | dim | initial_energy | final_energy | warmup_energy | points | plot_file |"
    )
    lines.append("| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |")
    for row in results:
        lines.append(
            "| {seed} | {n} | {dim} | {initial_energy:.6f} | {final_energy:.6f} | {warmup_energy:.6f} | {points} | {plot_file} |".format(
                **row
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_sweep_csv(results: Sequence[Dict[str, Any]], path: Path) -> None:
    """Write a CSV summary for batch runs."""

    headers = [
        "seed",
        "n",
        "dim",
        "initial_energy",
        "final_energy",
        "warmup_energy",
        "points",
        "input_file",
        "output_file",
        "plot_file",
        "causet_plot_file",
    ]
    lines = [",".join(headers)]
    for row in results:
        lines.append(
            ",".join(
                [
                    str(row.get("seed", "")),
                    str(row.get("n", "")),
                    str(row.get("dim", "")),
                    f"{row.get('initial_energy', 0.0):.6f}",
                    f"{row.get('final_energy', 0.0):.6f}",
                    f"{row.get('warmup_energy', 0.0):.6f}",
                    str(row.get("points", "")),
                    str(row.get("input_file", "")),
                    str(row.get("output_file", "")),
                    str(row.get("plot_file", "")),
                    str(row.get("causet_plot_file", "")),
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_single_case(
    z: List[List[bool]],
    dim: int,
    seed: int,
    output_path: Path,
    plot_path: Path | None,
    interactive: bool,
    max_data: int,
    backend: str,
) -> ConesSimulator:
    sim = ConesSimulator(
        z=z,
        dim=dim,
        seed=seed,
        interactive=interactive,
        max_data=max_data,
        plot_path=plot_path,
        backend=backend,
    )
    sim.run(output_path)
    return sim


def run_sweep(args: argparse.Namespace) -> int:
    if args.sprinkle is not None:
        base_z, base_points = generate_sprinkled_causet(
            args.sprinkle,
            seed=args.seed,
            spacetime_dim=args.sprinkle_spacetime_dim,
        )
        if args.causet_plot is not None:
            write_causet_svg(base_points, args.causet_plot, base_z)
        if args.save_generated_input is not None:
            n = len(base_z)
            entries = []
            for i in range(n - 1):
                for j in range(i + 1, n):
                    entries.append("1" if base_z[i][j] else "0")
            args.save_generated_input.write_text(
                str(n) + "\n" + " ".join(entries) + "\n",
                encoding="utf-8",
            )
    else:
        if args.input_file is None:
            raise ValueError("either input_file or --sprinkle must be provided")
        base_z = parse_cones_input(args.input_file)
        base_points = None
        if args.causet_plot is not None:
            # Reconstruct a simple t-x display from the matrix only if the user
            # wants a plot of the supplied causal structure.
            points = [(float(i), 0.0) for i in range(len(base_z))]
            write_causet_svg(points, args.causet_plot, base_z)

    results: List[Dict[str, Any]] = []
    sweep_dir = args.sweep_dir
    sweep_dir.mkdir(parents=True, exist_ok=True)
    for idx in range(args.sweep):
        seed = args.seed + idx * args.sweep_step
        prefix = sweep_dir / f"run_{idx:03d}"
        output_path = prefix.with_suffix(".out")
        plot_path = None if args.no_plot else prefix.with_suffix(".svg")
        sim = run_single_case(
            base_z,
            dim=args.dim,
            seed=seed,
            output_path=output_path,
            plot_path=plot_path,
            interactive=args.interactive,
            max_data=args.max_data,
            backend=args.backend,
        )
        row = sim.summary()
        row.update(
            {
                "output_file": str(output_path),
                "plot_file": str(plot_path) if plot_path is not None else "",
                "causet_plot_file": str(args.causet_plot or ""),
                "input_file": str(args.save_generated_input or args.input_file or ""),
            }
        )
        results.append(row)

    write_sweep_csv(results, args.sweep_report_csv)
    write_sweep_report(results, args.sweep_report_md)
    if args.sweep_summary_plot is not None:
        summary_points = [(float(r["seed"]), float(r["final_energy"])) for r in results]
        write_svg_plot(summary_points, args.sweep_summary_plot)
    return 0


@dataclass
class ConesSimulator:
    z: List[List[bool]]
    dim: int
    seed: int = 1959
    interactive: bool = False
    max_data: int = 35
    plot_path: Path | None = None
    backend: str = "cpu"
    shared_cuda_backend: Any | None = None
    warmup_limit: int = 100
    anneal_limit: int = 100
    initial_temp: float = 100.0
    cooling_factor: float = 0.9
    acceptance_scale: float = 4.0
    block_callback: Any | None = field(default=None, repr=False)
    rng: PascalRNG = field(init=False)
    cuda_backend: Any | None = field(init=False, default=None)

    def __post_init__(self) -> None:
        self.n = len(self.z)
        if self.n == 0:
            raise ValueError("empty incidence matrix")
        if self.dim <= 0:
            raise ValueError("dim must be positive")
        if self.warmup_limit <= 0:
            raise ValueError("warmup_limit must be positive")
        if self.anneal_limit <= 0:
            raise ValueError("anneal_limit must be positive")
        if self.initial_temp <= 0.0:
            raise ValueError("initial_temp must be positive")
        if not (0.0 < self.cooling_factor <= 1.0):
            raise ValueError("cooling_factor must be in (0, 1]")
        if self.acceptance_scale < 0.0:
            raise ValueError("acceptance_scale must be nonnegative")
        if self.backend not in ("cpu", "cuda", "auto"):
            raise ValueError("backend must be cpu, cuda, or auto")
        self.original_seed = self.seed
        self.rng = PascalRNG(self.seed)

        self.change = [True] * self.n
        self.xnew = [[0.0] * self.dim for _ in range(self.n)]
        self.xold = [[0.0] * self.dim for _ in range(self.n)]
        self.rnew = [0.0] * self.n
        self.rold = [0.0] * self.n
        self.energies = [0.0] * 100
        self.enew = [[0.0] * self.n for _ in range(self.n)]
        self.eold = [[0.0] * self.n for _ in range(self.n)]
        self.data: List[Tuple[float, float]] = []

        self.count = 0
        self.ndata = 0
        self.totalcount = 0
        self.deltae = 0.0
        self.rave = 0.0
        self.r = 0.0
        self.rmin = 0.0
        self.t = self.initial_temp
        self.eave = 0.0
        self.evar = 0.0
        self.spheat = 0.0
        self.initial_energy = 0.0
        self.warmup_energy = 0.0

        if self.shared_cuda_backend is not None:
            self.cuda_backend = self.shared_cuda_backend
        elif self.backend in ("cuda", "auto"):
            self._init_cuda_backend()

    @staticmethod
    def _sort(value: float) -> float:
        return math.sqrt(max(value, 0.0))

    def ran2(self) -> float:
        return self.rng.ran2()

    def gasdev(self) -> float:
        return self.rng.gasdev()

    def _init_cuda_backend(self) -> None:
        if CUDAEnergyBackend is None:
            if self.backend == "cuda":
                raise RuntimeError("CUDA backend requested but cuda_backend.py is unavailable")
            self.backend = "cpu"
            return
        try:
            self.cuda_backend = CUDAEnergyBackend(self.n, self.dim)
            self.cuda_backend.set_z(self.z)
        except Exception:
            if self.backend == "cuda":
                raise
            self.cuda_backend = None
            self.backend = "cpu"

    def startup(self, out: TextIO) -> None:
        for i in range(self.n):
            self.change[i] = True
            self.xnew[i] = [0.0] * self.dim
            self.rnew[i] = float(i + 2)
        self.rave = (self.n + 4) / 2.0
        self.energy()
        self.update()
        self.initial_energy = self.energies[0]
        self._writeln(out, f"Energy of initial configuration: {self.energies[0]:9.3f}")

    def energy(self) -> None:
        if self.cuda_backend is not None:
            flat = self.cuda_backend.compute(self.xnew, self.rnew, self.rave)
            self.deltae = 0.0
            for i in range(self.n - 1):
                base = i * self.n
                for j in range(i + 1, self.n):
                    value = flat[base + j]
                    self.enew[i][j] = value
                    self.deltae += value - self.eold[i][j]
            return

        roottwo = math.sqrt(2.0)
        self.deltae = 0.0

        for i in range(self.n - 1):
            for j in range(i + 1, self.n):
                if self.change[i] or self.change[j]:
                    rij = self.rnew[i] - self.rnew[j]
                    xij_sq = 0.0
                    for k in range(self.dim):
                        diff = self.xnew[i][k] - self.xnew[j][k]
                        xij_sq += diff * diff
                    s2 = -(rij ** 2) + xij_sq
                    xij = self._sort(xij_sq)
                    if self.z[i][j]:
                        if s2 > 0.0:
                            self.enew[i][j] = (xij + rij) / (roottwo * self.rave)
                        elif rij > 0.0:
                            self.enew[i][j] = math.sqrt(s2 + 2.0 * (rij ** 2)) / self.rave
                        else:
                            self.enew[i][j] = 0.0
                    else:
                        if s2 > 0.0:
                            self.enew[i][j] = 0.0
                        else:
                            self.enew[i][j] = (abs(rij) - xij) / (roottwo * self.rave)
                    self.deltae += self.enew[i][j] - self.eold[i][j]

    def reconfigure(self) -> None:
        self.rave = self.r
        for i in range(self.n):
            efraction = 2.0 * self.eold[i][i] / self.energies[0] if self.energies[0] else 0.0
            if self.ran2() < efraction:
                self.change[i] = True

            while True:
                while True:
                    norm = 0.0
                    displacement = [0.0] * (self.dim + 1)
                    for k in range(self.dim + 1):
                        displacement[k] = 2.0 * self.ran2() - 1.0
                        norm += displacement[k] * displacement[k]
                    if norm < 1.0:
                        break
                norm = self.eold[i][i] * self.r * self.gasdev() / math.sqrt(norm) if norm else 0.0
                for k in range(self.dim):
                    self.xnew[i][k] = self.xold[i][k] + displacement[k] * norm
                self.rnew[i] = self.rold[i] + displacement[self.dim] * norm
                if self.rnew[i] > 0.0:
                    break

            self.rave += (self.rnew[i] - self.rold[i]) / self.n

    def update(self) -> None:
        for idx in range(99, 0, -1):
            self.energies[idx] = self.energies[idx - 1]
        self.energies[0] = 0.0

        self.rmin = self.rnew[0]
        for i in range(self.n):
            if self.rnew[i] < self.rmin:
                self.rmin = self.rnew[i]

            if self.change[i]:
                for k in range(self.dim):
                    self.xold[i][k] = self.xnew[i][k]
                self.rold[i] = self.rnew[i]

            for j in range(i + 1, self.n):
                if self.change[i] or self.change[j]:
                    self.eold[i][j] = self.enew[i][j]

            self.eold[i][i] = 0.0
            for j in range(i):
                self.eold[i][i] += self.eold[j][i] / 2.0
            for j in range(i + 1, self.n):
                self.eold[i][i] += self.eold[i][j] / 2.0
            self.energies[0] += self.eold[i][i]

        for i in range(self.n):
            for k in range(self.dim):
                self.xold[i][k] /= self.rmin
            self.rold[i] /= self.rmin

        self.rave /= self.rmin
        self.r = self.rave

    def decide(self) -> bool:
        if self.deltae < 0.0:
            self.update()
            self.count += 1
            return True
        if self.t > 0.0 and self.ran2() < self.acceptance_scale * math.exp(-self.deltae / self.t):
            self.update()
            self.count += 1
            return True
        return False

    def statistics(self) -> None:
        self.eave = 0.0
        e2 = 0.0
        for value in self.energies:
            self.eave += value
            e2 += value * value
        self.eave /= 100.0
        e2 /= 100.0
        self.evar = math.sqrt(max(e2 - (self.eave ** 2), 0.0))
        if self.totalcount == 0:
            self.t = self.initial_temp
        if self.t > 0.0:
            self.spheat = (self.evar / self.t) ** 2
        self.totalcount += self.warmup_limit

    def warmup(self, out: TextIO) -> None:
        self.count = 0
        while self.count < self.warmup_limit and self.energies[0] > 0.0:
            self.reconfigure()
            self.energy()
            self.update()
            self.count += 1
            print(self.count)
        self.statistics()
        self._writeln(out, "Count Temp.   E[1].   Eave.   Evar.  spheat")
        self._writeln(
            out,
            f"warmup ------- {self.energies[0]:8.3f} {self.eave:8.3f} {self.evar:8.3f} {self.spheat:8.3f}",
        )
        self.warmup_energy = self.energies[0]
        if self.interactive:
            input("waiting for <cr> ... ")

    def anneal(self, out: TextIO) -> None:
        self.ndata = 0
        continue_flag = True
        while continue_flag:
            self.count = 0
            while self.count < self.anneal_limit and self.energies[0] > 0.0:
                self.reconfigure()
                self.energy()
                self.decide()

            self.statistics()
            self.ndata += 1
            self.data.append((self.t, self.eave))
            if self.block_callback is not None:
                self.block_callback(self, self.ndata, self.t, self.eave)

            if self.interactive:
                reply = input("y/t/n: ").strip().lower()[:1] or "n"
                if reply == "t":
                    self.t = float(input("new T: "))
                    continue_flag = True
                elif reply == "y":
                    self.t *= self.cooling_factor
                    continue_flag = True
                else:
                    continue_flag = False
            else:
                if self.ndata >= self.max_data or self.energies[0] <= 0.0:
                    continue_flag = False
                else:
                    self.t *= self.cooling_factor

    def writeout(self, where: TextIO) -> None:
        self._writeln(where, "")
        self._writeln(where, "Incidence matrix:")
        self._writeln(where, "")
        self._writeln(where, format_matrix_bool(self.z))
        self._writeln(where, "")
        self._writeln(where, f"dimension: {self.dim:2d}")
        self._writeln(where, "")
        for i in range(self.n):
            coords = " ".join(f"{v:8.3f}" for v in self.xold[i])
            self._writeln(where, f"i={i + 1:2d} X= {coords}  R={self.rold[i]:8.3f}")
        self._writeln(where, "")
        self._writeln(where, "Energy matrix:")
        self._writeln(where, "")
        self._writeln(where, format_matrix_float(self.eold))
        self._writeln(where, "")
        self._writeln(where, f"Total energy: {self.energies[0]:.6f}")
        self._writeln(
            where,
            (
                f"Schedule: warmup_limit={self.warmup_limit} "
                f"anneal_limit={self.anneal_limit} "
                f"initial_temp={self.initial_temp:.3f} "
                f"cooling_factor={self.cooling_factor:.6f} "
                f"acceptance_scale={self.acceptance_scale:.3f}"
            ),
        )

    @staticmethod
    def _writeln(where: TextIO, text: str = "") -> None:
        where.write(text + "\n")

    def run(self, output_path: Path) -> None:
        try:
            with output_path.open("w", encoding="utf-8") as out:
                self.startup(out)
                self.warmup(out)
                self.anneal(out)
                self.writeout(sys.stdout)
                self.writeout(out)
            if self.plot_path is not None:
                write_svg_plot(self.data, self.plot_path)
        finally:
            if self.cuda_backend is not None and self.shared_cuda_backend is None:
                self.cuda_backend.close()

    def summary(self) -> Dict[str, Any]:
        final_energy = self.data[-1][1] if self.data else self.eave
        return {
            "seed": self.original_seed,
            "n": self.n,
            "dim": self.dim,
            "initial_energy": self.initial_energy,
            "final_energy": final_energy,
            "warmup_energy": self.warmup_energy,
            "points": len(self.data),
        }


def build_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run the revived cones simulation.")
    parser.add_argument("input_file", type=Path, nargs="?", help="Pascal-style incidence input file")
    parser.add_argument("--sprinkle", type=int, help="generate a random causet with this many events")
    parser.add_argument(
        "--sprinkle-spacetime-dim",
        type=int,
        default=1,
        help="spatial dimension used by the sprinkling generator",
    )
    parser.add_argument("--dim", type=int, help="space dimension")
    parser.add_argument("--output", type=Path, default=Path("cone.out"), help="output file")
    parser.add_argument("--seed", type=int, default=1959, help="random seed")
    parser.add_argument("--interactive", action="store_true", help="use the old prompt-driven annealing loop")
    parser.add_argument("--max-data", type=int, default=35, help="maximum annealing data points in noninteractive mode")
    parser.add_argument("--backend", choices=["cpu", "cuda", "auto"], default="cpu", help="energy backend")
    parser.add_argument("--warmup-limit", type=int, default=100, help="maximum reconfigure steps during warmup")
    parser.add_argument("--anneal-limit", type=int, default=100, help="maximum reconfigure steps per anneal block")
    parser.add_argument("--initial-temp", type=float, default=100.0, help="initial temperature")
    parser.add_argument("--cooling-factor", type=float, default=0.9, help="temperature multiplier between anneal blocks")
    parser.add_argument("--acceptance-scale", type=float, default=4.0, help="scale factor in the Metropolis test")
    parser.add_argument("--plot", type=Path, default=Path("cone.svg"), help="write a simple SVG plot of T vs Eave")
    parser.add_argument("--no-plot", action="store_true", help="disable SVG plot generation")
    parser.add_argument("--causet-plot", type=Path, help="write an SVG plot of the generated causet points")
    parser.add_argument("--save-generated-input", type=Path, help="write the generated incidence matrix to a file")
    parser.add_argument("--sweep", type=int, default=0, help="run a batch of simulations with different seeds")
    parser.add_argument("--sweep-step", type=int, default=1, help="seed increment between sweep runs")
    parser.add_argument("--sweep-dir", type=Path, default=Path("sweep_runs"), help="directory for per-run outputs")
    parser.add_argument("--sweep-report-csv", type=Path, default=Path("sweep_report.csv"), help="CSV summary for sweeps")
    parser.add_argument("--sweep-report-md", type=Path, default=Path("sweep_report.md"), help="markdown summary for sweeps")
    parser.add_argument("--sweep-summary-plot", type=Path, default=Path("sweep_summary.svg"), help="SVG summary plot for sweeps")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_argparser()
    args = parser.parse_args(argv)

    if args.sweep > 0:
        if args.dim is None:
            parser.error("--dim is required when using --sweep")
        return run_sweep(args)

    points: List[Tuple[float, ...]] | None = None
    if args.sprinkle is not None:
        z, points = generate_sprinkled_causet(
            args.sprinkle,
            seed=args.seed,
            spacetime_dim=args.sprinkle_spacetime_dim,
        )
        if args.save_generated_input is not None:
            n = len(z)
            entries = []
            for i in range(n - 1):
                for j in range(i + 1, n):
                    entries.append("1" if z[i][j] else "0")
            args.save_generated_input.write_text(
                str(n) + "\n" + " ".join(entries) + "\n",
                encoding="utf-8",
            )
        if args.causet_plot is not None and points is not None:
            write_causet_svg(points, args.causet_plot, z)
    else:
        if args.input_file is None:
            parser.error("either input_file or --sprinkle must be provided")
        z = parse_cones_input(args.input_file)

    dim = args.dim
    if dim is None:
        dim = int(input("how many space dimensions? "))

    sim = ConesSimulator(
        z=z,
        dim=dim,
        seed=args.seed,
        interactive=args.interactive,
        max_data=args.max_data,
        plot_path=None if args.no_plot else args.plot,
        backend=args.backend,
        warmup_limit=args.warmup_limit,
        anneal_limit=args.anneal_limit,
        initial_temp=args.initial_temp,
        cooling_factor=args.cooling_factor,
        acceptance_scale=args.acceptance_scale,
    )
    sim.run(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
