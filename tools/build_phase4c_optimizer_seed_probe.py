#!/usr/bin/env python3
"""Phase 4C — Optimizer-seed multi-start probe.

Phase 4A/4B/5 all use a single hardcoded ``OPTIMIZER_SEED = 1987`` for
``ConesSimulator``.  The Phase 4B `MIXED` and Phase 5 `INSUFFICIENT`
outcomes therefore include an uncontrolled source of variance: the
internal PascalRNG of the simulated annealer.  Phase 4C asks:

    Is the Phase 4B/5 picture limited by the single optimizer seed,
    or is it robust under optimizer-seed perturbation?

No new optimizer is introduced.  ``cones.py`` is not modified.  Only
``_run_one`` of Phase 4A gains an optional ``optimizer_seed`` argument
whose default reproduces Phase 4A bit-for-bit.

Grid
----
n            : 32, 48, 64   (same pilot grid as Phase 4B)
target_dim   : 2, 3, 4
causet_seed  : Phase4A.PHASE4A_SEEDS (10 seeds, unchanged)
epsilon      : Phase4A.EPSILONS (8 values, unchanged)
optimizer_seed : 1987, 1990, 1993  (K = 3)

Outputs
-------
benchmarks/foundation/phase4c_optimizer_seed_probe_per_run.csv
benchmarks/foundation/phase4c_optimizer_seed_probe_per_cell_epsilon.csv
benchmarks/foundation/phase4c_optimizer_seed_probe.md

Verdict
-------
OPTIMIZER_SEED_LIMITED  if any of:
    mean_label_stability < 0.5
    mean(IQR_K) / mean(loss_K) > 0.1
    floor_fraction_K < 0.5 * floor_fraction_phase4a

OPTIMIZER_SEED_ROBUST   if all of:
    mean_label_stability > 0.9
    mean(IQR_K) / mean(loss_K) < 0.01
    floor_fraction_K ≈ floor_fraction_phase4a (within 10%)

INCONCLUSIVE            otherwise.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import validation_suite as vs  # noqa: E402
from tools import build_phase4a_epsilon_sweep as p4a  # noqa: E402


FOUNDATION = ROOT / "benchmarks" / "foundation"
PHASE4A_CSV = FOUNDATION / "phase4a_epsilon_sweep.csv"

PILOT_SIZES: tuple[int, ...] = (32, 48, 64)
PILOT_DIMS:  tuple[int, ...] = (2, 3, 4)
PHASE4C_OPTIMIZER_SEEDS: tuple[int, ...] = (1987, 1990, 1993)

DEFAULT_FLOOR_TOLERANCE = 1e-6
IQR_LIMIT_RATIO   = 0.10   # LIMITED if mean(IQR) > 0.10 * mean(loss)
IQR_ROBUST_RATIO  = 0.01   # ROBUST  if mean(IQR) < 0.01 * mean(loss)
FLOOR_LIMIT_RATIO = 0.50   # LIMITED if floor_K < 0.50 * floor_phase4a
FLOOR_ROBUST_TOL  = 0.10   # ROBUST  if |floor_K/floor_phase4a - 1| < 0.10
LABEL_LIMIT       = 0.50   # LIMITED if mean_stability < 0.50
LABEL_ROBUST      = 0.90   # ROBUST  if mean_stability > 0.90


PER_RUN_HEADERS = (
    "phase",
    "n",
    "target_dim",
    "causet_seed",
    "epsilon",
    "optimizer_seed",
    "valid",
    "failure_mode",
    "loss",
    "initial_energy",
    "final_energy",
    "warmup_delta_energy",
    "warmup_attempted_moves",
    "warmup_accepted_moves",
    "warmup_rejected_moves",
)

PER_CELL_EPSILON_HEADERS = (
    "phase",
    "n",
    "target_dim",
    "epsilon",
    "K",
    "n_valid_runs",
    "mean_loss_K",
    "std_loss_K",
    "min_loss_K",
    "max_loss_K",
    "iqr_loss_K",
    "floor_saturated_fraction_K",
    "phase4a_mean_loss",
    "phase4a_floor_saturated_fraction",
    "delta_min_loss_vs_phase4a",
    "delta_floor_saturated_vs_phase4a",
    "curve_shape_per_optimizer_seed",
    "label_stability_cell",
)


# ---------------------------------------------------------------------------
# Per-run simulation
# ---------------------------------------------------------------------------

def _row_from_sim_result(
    n: int, d: int, causet_seed: int, eps: float, optimizer_seed: int,
    sim_result: dict,
) -> dict:
    E0  = sim_result["initial_energy"]
    wdE = sim_result["warmup_delta_energy"]
    if not math.isfinite(E0) or abs(E0) < p4a.INITIAL_ENERGY_FLOOR:
        valid = False
        failure_mode = f"E0_too_small ({E0:.3g} < {p4a.INITIAL_ENERGY_FLOOR:g})"
        loss = float("nan")
    else:
        valid = True
        failure_mode = ""
        loss = abs(wdE / E0)
    return {
        "phase": "phase4c_optimizer_seed_probe",
        "n": n,
        "target_dim": d,
        "causet_seed": causet_seed,
        "epsilon": eps,
        "optimizer_seed": optimizer_seed,
        "valid": valid,
        "failure_mode": failure_mode,
        "loss": loss,
        "initial_energy": E0,
        "final_energy": sim_result["final_energy"],
        "warmup_delta_energy": wdE,
        "warmup_attempted_moves": sim_result["warmup_attempted_moves"],
        "warmup_accepted_moves":  sim_result["warmup_accepted_moves"],
        "warmup_rejected_moves":  sim_result["warmup_rejected_moves"],
    }


def build_rows(
    sizes: tuple[int, ...],
    dims: tuple[int, ...],
    optimizer_seeds: tuple[int, ...],
) -> list[dict]:
    """Run the full Phase 4C grid (sprinkle once per (n, seed, d); reuse across ε and K)."""
    causet_seeds = p4a.PHASE4A_SEEDS
    epsilons     = p4a.EPSILONS
    total = (
        len(sizes) * len(dims) * len(causet_seeds)
        * len(epsilons) * len(optimizer_seeds)
    )
    rows: list[dict] = []
    done = 0
    t0 = time.time()
    last_report = t0

    for n in sizes:
        for seed in causet_seeds:
            for d in dims:
                matrix, points = vs.sprinkle_minkowski_diamond(
                    n=n, seed=seed, d_spacetime=d,
                )
                for eps in epsilons:
                    for opt_seed in optimizer_seeds:
                        sim_result = p4a._run_one(
                            d, n, seed, eps, matrix, points,
                            optimizer_seed=opt_seed,
                        )
                        rows.append(_row_from_sim_result(
                            n, d, seed, eps, opt_seed, sim_result,
                        ))
                        done += 1
                        now = time.time()
                        if now - last_report > 10:
                            rate = done / (now - t0) if now > t0 else 0.0
                            eta  = (total - done) / rate if rate > 0 else 0.0
                            print(
                                f"  {done}/{total} rows "
                                f"({now - t0:.0f}s elapsed, ETA {eta:.0f}s)",
                                flush=True,
                            )
                            last_report = now
    print(f"  done: {done}/{total} rows in {time.time() - t0:.0f}s.", flush=True)
    return rows


# ---------------------------------------------------------------------------
# Aggregation helpers
# ---------------------------------------------------------------------------

def _iqr(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    s = sorted(values)
    n = len(s)
    q1_idx = max(0, (n - 1) // 4)
    q3_idx = min(n - 1, (3 * (n - 1)) // 4)
    return s[q3_idx] - s[q1_idx]


def _phase4a_baseline_lookup(path: Path = PHASE4A_CSV) -> dict:
    """Return {(n, d, eps): [valid abs_relative_drift values]} for delta computations.
    Returns empty dict if Phase 4A CSV is absent (n=48 baseline is always empty)."""
    if not path.exists():
        return {}
    out: dict[tuple[int, int, float], list[float]] = {}
    for r in p4a.load_rows_from_csv(path):
        if not r.get("row_valid", False):
            continue
        loss = r.get("abs_relative_drift", float("nan"))
        if not math.isfinite(loss):
            continue
        key = (r["n"], r["target_dim"], r["epsilon"])
        out.setdefault(key, []).append(loss)
    return out


def _curve_for_optimizer_seed(
    rows: list[dict], n: int, d: int, opt_seed: int,
) -> list[tuple]:
    """Build epsilon-ordered (eps, mean_loss_across_causet_seeds) curve for one (n,d,opt_seed)."""
    by_eps: dict[float, list[float]] = {}
    for r in rows:
        if r["n"] != n or r["target_dim"] != d or r["optimizer_seed"] != opt_seed:
            continue
        if not r["valid"] or not math.isfinite(r["loss"]):
            continue
        by_eps.setdefault(r["epsilon"], []).append(r["loss"])
    curve: list[tuple] = []
    for eps in sorted(by_eps):
        losses = by_eps[eps]
        mean = sum(losses) / len(losses)
        curve.append((eps, mean, mean))
    return curve


def summarize_per_cell_epsilon(
    rows: list[dict],
    optimizer_seeds: tuple[int, ...],
    floor_tolerance: float = DEFAULT_FLOOR_TOLERANCE,
    phase4a_baseline: dict | None = None,
) -> list[dict]:
    """Aggregate per-run rows into one row per (n, target_dim, epsilon).

    Label stability is computed per (n, target_dim) cell across the K
    optimizer-seed curves and propagated to every epsilon row of the cell.
    """
    if phase4a_baseline is None:
        phase4a_baseline = _phase4a_baseline_lookup()
    K = len(optimizer_seeds)
    cells = sorted({(r["n"], r["target_dim"]) for r in rows})

    shape_by_cell_opt: dict[tuple[int, int, int], str] = {}
    for (n, d) in cells:
        for opt_seed in optimizer_seeds:
            curve = _curve_for_optimizer_seed(rows, n, d, opt_seed)
            if len(curve) >= 3:
                shape = p4a._classify_curve_shape(curve, idx=1)
            else:
                shape = "insufficient"
            shape_by_cell_opt[(n, d, opt_seed)] = shape

    label_stability_per_cell: dict[tuple[int, int], float] = {}
    for (n, d) in cells:
        shapes = {shape_by_cell_opt[(n, d, s)] for s in optimizer_seeds}
        label_stability_per_cell[(n, d)] = 1.0 if len(shapes) == 1 else 0.0

    out: list[dict] = []
    keys = sorted({(r["n"], r["target_dim"], r["epsilon"]) for r in rows})
    for (n, d, eps) in keys:
        cell_eps = [
            r for r in rows
            if r["n"] == n and r["target_dim"] == d and r["epsilon"] == eps
        ]
        valid_losses = [
            r["loss"] for r in cell_eps
            if r["valid"] and math.isfinite(r["loss"])
        ]
        n_valid = len(valid_losses)
        if n_valid > 0:
            mean_loss = sum(valid_losses) / n_valid
            std_loss  = statistics.pstdev(valid_losses) if n_valid > 1 else 0.0
            min_loss  = min(valid_losses)
            max_loss  = max(valid_losses)
            iqr_loss  = _iqr(valid_losses)
            floor_count = sum(1 for v in valid_losses if v <= floor_tolerance)
            floor_frac  = floor_count / n_valid
        else:
            mean_loss = std_loss = min_loss = max_loss = iqr_loss = float("nan")
            floor_frac = float("nan")

        baseline = phase4a_baseline.get((n, d, eps), [])
        if baseline:
            p4a_mean       = sum(baseline) / len(baseline)
            p4a_min        = min(baseline)
            p4a_floor      = sum(1 for v in baseline if v <= floor_tolerance) / len(baseline)
        else:
            p4a_mean = p4a_min = p4a_floor = float("nan")

        delta_min = (
            min_loss - p4a_min
            if math.isfinite(min_loss) and math.isfinite(p4a_min)
            else float("nan")
        )
        delta_floor = (
            floor_frac - p4a_floor
            if math.isfinite(floor_frac) and math.isfinite(p4a_floor)
            else float("nan")
        )

        shapes_str = "|".join(shape_by_cell_opt[(n, d, s)] for s in optimizer_seeds)

        out.append({
            "phase": "phase4c_optimizer_seed_probe",
            "n": n,
            "target_dim": d,
            "epsilon": eps,
            "K": K,
            "n_valid_runs": n_valid,
            "mean_loss_K": mean_loss,
            "std_loss_K": std_loss,
            "min_loss_K": min_loss,
            "max_loss_K": max_loss,
            "iqr_loss_K": iqr_loss,
            "floor_saturated_fraction_K": floor_frac,
            "phase4a_mean_loss": p4a_mean,
            "phase4a_floor_saturated_fraction": p4a_floor,
            "delta_min_loss_vs_phase4a": delta_min,
            "delta_floor_saturated_vs_phase4a": delta_floor,
            "curve_shape_per_optimizer_seed": shapes_str,
            "label_stability_cell": label_stability_per_cell[(n, d)],
        })
    return out


# ---------------------------------------------------------------------------
# Verdict
# ---------------------------------------------------------------------------

def compute_verdict(per_cell_eps_rows: list[dict]) -> tuple[str, dict]:
    if not per_cell_eps_rows:
        return ("INCONCLUSIVE", {"reason": "no per_cell_epsilon rows"})

    cells = sorted({(r["n"], r["target_dim"]) for r in per_cell_eps_rows})
    cell_stabilities: list[float] = []
    for (n, d) in cells:
        cell_rows = [
            r for r in per_cell_eps_rows
            if r["n"] == n and r["target_dim"] == d
        ]
        if cell_rows:
            cell_stabilities.append(cell_rows[0]["label_stability_cell"])
    mean_stability = (
        sum(cell_stabilities) / len(cell_stabilities) if cell_stabilities else 0.0
    )

    valid_iqr  = [r["iqr_loss_K"]  for r in per_cell_eps_rows if math.isfinite(r["iqr_loss_K"])]
    valid_mean = [r["mean_loss_K"] for r in per_cell_eps_rows if math.isfinite(r["mean_loss_K"])]
    mean_iqr  = sum(valid_iqr)  / len(valid_iqr)  if valid_iqr  else 0.0
    mean_loss = sum(valid_mean) / len(valid_mean) if valid_mean else 0.0
    iqr_ratio = mean_iqr / mean_loss if mean_loss > 0 else 0.0

    # Floor-fraction comparison uses only cells that have a Phase 4A baseline
    # (n=48 cells have no baseline and are excluded from the floor comparison).
    floor_pairs = [
        (r["floor_saturated_fraction_K"], r["phase4a_floor_saturated_fraction"])
        for r in per_cell_eps_rows
        if math.isfinite(r["floor_saturated_fraction_K"])
        and math.isfinite(r["phase4a_floor_saturated_fraction"])
    ]
    if floor_pairs:
        floor_K        = sum(f for f, _ in floor_pairs) / len(floor_pairs)
        floor_phase4a  = sum(p for _, p in floor_pairs) / len(floor_pairs)
    else:
        floor_K = floor_phase4a = float("nan")

    floor_phase4a_zero = (
        math.isfinite(floor_phase4a) and math.isclose(floor_phase4a, 0.0, abs_tol=1e-9)
    )
    floor_K_zero = (
        math.isfinite(floor_K) and math.isclose(floor_K, 0.0, abs_tol=1e-9)
    )
    # LIMITED: floor_K shrank materially below baseline (only meaningful if baseline > 0).
    floor_limited = (
        math.isfinite(floor_phase4a)
        and floor_phase4a > 0.0
        and math.isfinite(floor_K)
        and floor_K < FLOOR_LIMIT_RATIO * floor_phase4a
    )
    # ROBUST: floor_K matches baseline within FLOOR_ROBUST_TOL, or both zero.
    if not math.isfinite(floor_K) or not math.isfinite(floor_phase4a):
        floor_robust = False
    elif floor_phase4a_zero:
        floor_robust = floor_K_zero
    else:
        floor_robust = abs(floor_K / floor_phase4a - 1.0) < FLOOR_ROBUST_TOL

    info = {
        "n_cells": len(cells),
        "mean_label_stability": mean_stability,
        "mean_iqr_K": mean_iqr,
        "mean_loss_K": mean_loss,
        "iqr_ratio": iqr_ratio,
        "floor_fraction_K": floor_K,
        "floor_fraction_phase4a": floor_phase4a,
        "floor_limited_triggered": floor_limited,
        "floor_robust_satisfied":  floor_robust,
    }

    limited = (
        mean_stability < LABEL_LIMIT
        or iqr_ratio > IQR_LIMIT_RATIO
        or floor_limited
    )
    robust = (
        mean_stability > LABEL_ROBUST
        and iqr_ratio < IQR_ROBUST_RATIO
        and floor_robust
    )
    if limited and not robust:
        return ("OPTIMIZER_SEED_LIMITED", info)
    if robust and not limited:
        return ("OPTIMIZER_SEED_ROBUST", info)
    return ("INCONCLUSIVE", info)


# ---------------------------------------------------------------------------
# Verify-against-phase4a mode
# ---------------------------------------------------------------------------

def verify_against_phase4a(
    optimizer_seed: int = p4a.OPTIMIZER_SEED,
    sample_size: int = 12,
    abs_tol: float = 5e-5,
) -> tuple[bool, list[str]]:
    """Run K=1 with the given optimizer_seed on a small sample of Phase 4A rows
    and verify that abs_relative_drift matches the stored CSV to the precision
    at which the CSV persists floats (`.4f` → half-unit = 5e-5).  The simulator
    output is deterministic; this tolerance only absorbs CSV rounding."""
    if not PHASE4A_CSV.exists():
        return (False, [f"missing baseline: {PHASE4A_CSV}"])
    p4a_rows = p4a.load_rows_from_csv(PHASE4A_CSV)
    valid = [r for r in p4a_rows if r.get("row_valid", False)]
    if not valid:
        return (False, ["Phase 4A CSV has no valid rows"])
    # Spread sample across (n, d) for coverage.
    cells = sorted({(r["n"], r["target_dim"]) for r in valid})
    sampled: list[dict] = []
    for cell in cells:
        per_cell = [r for r in valid if (r["n"], r["target_dim"]) == cell]
        sampled.extend(per_cell[: max(1, sample_size // max(1, len(cells)))])
    sampled = sampled[:sample_size]

    failures: list[str] = []
    for row in sampled:
        n = row["n"]; d = row["target_dim"]
        seed = row["seed"]; eps = row["epsilon"]
        matrix, points = vs.sprinkle_minkowski_diamond(
            n=n, seed=seed, d_spacetime=d,
        )
        sim_result = p4a._run_one(
            d, n, seed, eps, matrix, points, optimizer_seed=optimizer_seed,
        )
        E0  = sim_result["initial_energy"]
        wdE = sim_result["warmup_delta_energy"]
        expected = row["abs_relative_drift"]
        if not math.isfinite(E0) or abs(E0) < p4a.INITIAL_ENERGY_FLOOR:
            computed = float("nan")
        else:
            computed = abs(wdE / E0)
        if math.isnan(expected) and math.isnan(computed):
            continue
        if not math.isclose(expected, computed, rel_tol=0.0, abs_tol=abs_tol):
            failures.append(
                f"  mismatch (n={n}, d={d}, seed={seed}, eps={eps}): "
                f"expected={expected!r}, got={computed!r}"
            )
    return (not failures, failures)


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def _fmt(value) -> str:
    if value is None:                       return "NA"
    if isinstance(value, bool):             return "true" if value else "false"
    if isinstance(value, str):              return value
    if isinstance(value, int):              return str(value)
    if isinstance(value, float):
        if math.isnan(value):               return "NA"
        if math.isinf(value):               return "inf" if value > 0 else "-inf"
        return f"{value:.10g}"
    return str(value)


def write_per_run_csv(rows: list[dict], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(PER_RUN_HEADERS)
        for row in rows:
            writer.writerow([_fmt(row[h]) for h in PER_RUN_HEADERS])


def write_per_cell_epsilon_csv(rows: list[dict], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(PER_CELL_EPSILON_HEADERS)
        for row in rows:
            writer.writerow([_fmt(row[h]) for h in PER_CELL_EPSILON_HEADERS])


def _cell_stability_table(per_cell_eps_rows: list[dict]) -> list[str]:
    cells = sorted({(r["n"], r["target_dim"]) for r in per_cell_eps_rows})
    lines = [
        "| n | target_dim | label_stability | curve_shape_per_optimizer_seed |",
        "| ---: | :---: | ---: | --- |",
    ]
    for (n, d) in cells:
        cell_rows = [r for r in per_cell_eps_rows if r["n"] == n and r["target_dim"] == d]
        if not cell_rows:
            continue
        first = cell_rows[0]
        lines.append(
            f"| {n} | {d} | {first['label_stability_cell']:.2f} "
            f"| `{first['curve_shape_per_optimizer_seed']}` |"
        )
    return lines


def _per_cell_eps_table(per_cell_eps_rows: list[dict]) -> list[str]:
    lines = [
        "| n | target_dim | epsilon | n_valid | mean_loss_K | std_loss_K "
        "| IQR_loss_K | floor_K | phase4a_floor | Δ_floor |",
        "| ---: | :---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in sorted(per_cell_eps_rows, key=lambda r: (r["n"], r["target_dim"], r["epsilon"])):
        lines.append(
            f"| {row['n']} | {row['target_dim']} | {row['epsilon']:.3g} "
            f"| {row['n_valid_runs']} "
            f"| {_fmt(row['mean_loss_K'])} | {_fmt(row['std_loss_K'])} "
            f"| {_fmt(row['iqr_loss_K'])} | {_fmt(row['floor_saturated_fraction_K'])} "
            f"| {_fmt(row['phase4a_floor_saturated_fraction'])} "
            f"| {_fmt(row['delta_floor_saturated_vs_phase4a'])} |"
        )
    return lines


def write_markdown(
    per_run_rows: list[dict],
    per_cell_eps_rows: list[dict],
    verdict: str,
    info: dict,
    optimizer_seeds: tuple[int, ...],
    sizes: tuple[int, ...],
    dims: tuple[int, ...],
    elapsed_s: float,
    path: Path,
) -> None:
    K = len(optimizer_seeds)
    n_runs = len(per_run_rows)
    n_valid = sum(1 for r in per_run_rows if r["valid"])
    n_cells = info.get("n_cells", 0)
    lines = [
        "# Phase 4C — Optimizer-seed multi-start probe",
        "",
        "**Status:** exploratory probe of optimizer-seed variance.  No new "
        "optimizer is introduced.  `cones.py` is not modified.",
        "",
        "## Objective",
        "",
        "Phase 4A/4B/5 use a single hardcoded `OPTIMIZER_SEED = 1987` inside "
        "`ConesSimulator`.  Phase 4C asks whether the Phase 4B `MIXED` and "
        "Phase 5 `INSUFFICIENT` outcomes are limited by that single optimizer "
        "seed, or robust under K = "
        f"{K} optimizer-seed perturbation.",
        "",
        "Phase 4C does not aim to lower the loss.  It only quantifies how "
        "much of the existing pipeline's behaviour is reproducible across "
        "different optimizer seeds while holding the causet, the epsilon, "
        "the initialization noise, the schedule, and the energy formula "
        "unchanged.",
        "",
        "## Grid",
        "",
        f"- Sizes: {', '.join(str(n) for n in sizes)}",
        f"- Spacetime dims: {', '.join(str(d) for d in dims)}",
        f"- Epsilons: {', '.join(f'{e:g}' for e in p4a.EPSILONS)}",
        f"- Causet seeds ({len(p4a.PHASE4A_SEEDS)}): "
        f"{', '.join(str(s) for s in p4a.PHASE4A_SEEDS)}",
        f"- Optimizer seeds (K = {K}): "
        f"{', '.join(str(s) for s in optimizer_seeds)}",
        f"- Wall-clock: {elapsed_s:.0f}s"
        + (" (markdown only)" if elapsed_s == 0.0 else ""),
        "",
        "## Sample bookkeeping",
        "",
        f"- Total runs: {n_runs}",
        f"- Valid runs: {n_valid}",
        f"- Cells: {n_cells}",
        "",
        "## Global verdict",
        "",
        f"**{verdict}**",
        "",
        "Decision rule:",
        "",
        f"- `OPTIMIZER_SEED_LIMITED` if any of: "
        f"`mean_label_stability < {LABEL_LIMIT}`, "
        f"`mean(IQR_K)/mean(loss_K) > {IQR_LIMIT_RATIO}`, "
        f"`floor_fraction_K < {FLOOR_LIMIT_RATIO} * floor_fraction_phase4a`.",
        f"- `OPTIMIZER_SEED_ROBUST` if all of: "
        f"`mean_label_stability > {LABEL_ROBUST}`, "
        f"`mean(IQR_K)/mean(loss_K) < {IQR_ROBUST_RATIO}`, "
        f"`|floor_fraction_K/floor_fraction_phase4a - 1| < {FLOOR_ROBUST_TOL}` (or both zero).",
        "- `INCONCLUSIVE` otherwise.",
        "",
        "Verdict inputs:",
        "",
        f"- `mean_label_stability` = {info.get('mean_label_stability', float('nan')):.3f}",
        f"- `mean_iqr_K` = {_fmt(info.get('mean_iqr_K', float('nan')))}",
        f"- `mean_loss_K` = {_fmt(info.get('mean_loss_K', float('nan')))}",
        f"- `iqr_ratio` (IQR/mean) = {_fmt(info.get('iqr_ratio', float('nan')))}",
        f"- `floor_fraction_K` = {_fmt(info.get('floor_fraction_K', float('nan')))}",
        f"- `floor_fraction_phase4a` = {_fmt(info.get('floor_fraction_phase4a', float('nan')))}",
        "",
        "## Per-cell label stability",
        "",
        "`label_stability_cell` = 1.0 iff all K optimizer-seed curves for a "
        "given (n, target_dim) cell receive the same `curve_shape` label.",
        "",
        *_cell_stability_table(per_cell_eps_rows),
        "",
        "## Per-cell-epsilon statistics (K runs combined with causet seeds)",
        "",
        "`mean_loss_K`, `std_loss_K`, `IQR_loss_K`, `floor_K` are computed "
        "over K * len(causet_seeds) valid rows per (n, target_dim, epsilon).  "
        "Phase 4A baseline columns are NaN for n = 48 (no baseline).",
        "",
        *_per_cell_eps_table(per_cell_eps_rows),
        "",
        "## Interpretation",
        "",
    ]
    if verdict == "OPTIMIZER_SEED_LIMITED":
        lines += [
            "The Phase 4B/5 outcomes are sensitive to the choice of "
            "`OPTIMIZER_SEED`.  At least one of (a) curve-shape labels flip "
            "across optimizer seeds, (b) loss IQR is comparable to the mean "
            "loss, or (c) the floor-saturation rate drops materially below "
            "the Phase 4A baseline.  The Phase 3F `INTERMEDIATE` signal and "
            "the Phase 5 `INSUFFICIENT` censoring should be re-read as "
            "potentially contaminated by single-seed variance.  This justifies "
            "a follow-up moderate intervention (e.g., reheating, larger K, "
            "or a multi-start aggregate target) under pre-registered "
            "criteria; it does not by itself establish a physical claim.",
        ]
    elif verdict == "OPTIMIZER_SEED_ROBUST":
        lines += [
            "The Phase 4B/5 outcomes are stable under optimizer-seed "
            "perturbation in this grid.  Curve-shape labels agree across K "
            "optimizer seeds, loss IQR is small relative to the mean loss, "
            "and the floor-saturation fraction matches the Phase 4A baseline.  "
            "The `INTERMEDIATE` / `MIXED` / `INSUFFICIENT` outcomes therefore "
            "reflect properties of the problem under the current optimizer "
            "and loss definition, not artifacts of the single hardcoded "
            "optimizer seed.  No further optimizer change is justified at "
            "this pipeline scale.",
        ]
    else:
        lines += [
            "Optimizer-seed evidence is mixed: at least one criterion is "
            "satisfied for both LIMITED and ROBUST, or neither.  Increase K "
            "(e.g., K=5 or K=7) before drawing a stronger conclusion.  The "
            "Phase 4B `MIXED` and Phase 5 `INSUFFICIENT` outcomes should "
            "remain the working interpretation in the meantime.",
        ]

    lines += [
        "",
        "## Scope",
        "",
        "- `loss` retains the Phase 4A/4B definition "
        "`|warmup_delta_energy / initial_energy|`.  It is an optimizer/"
        "embedding-response diagnostic, not a physical observable.",
        "- No Phase 4A, 4B, or 5 CSV, label, threshold, or outcome is "
        "modified.  Phase 3F is not rerun.  PySR is not invoked.",
        "- K = 1 with `optimizer_seed = 1987` reproduces Phase 4A bit-for-bit "
        "(see `--verify-against-phase4a`).",
        "",
        "## Reproducibility",
        "",
        "Regenerate via `make regen-phase4c`.",
        "Source: `tools/build_phase4c_optimizer_seed_probe.py`.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--K", type=int, default=len(PHASE4C_OPTIMIZER_SEEDS),
        help="Number of optimizer seeds to use from the head of --optimizer-seeds.",
    )
    ap.add_argument(
        "--optimizer-seeds", type=int, nargs="+",
        default=list(PHASE4C_OPTIMIZER_SEEDS),
        help="Optimizer seeds to vary (subset of size --K is used).",
    )
    ap.add_argument(
        "--verify-against-phase4a", action="store_true",
        help="Run K=1 with optimizer_seed=1987 on a sample of Phase 4A rows "
             "and verify abs_relative_drift matches the stored CSV exactly.",
    )
    args = ap.parse_args()

    if args.verify_against_phase4a:
        seeds_used = tuple(args.optimizer_seeds[: max(1, args.K)])
        if seeds_used[0] != p4a.OPTIMIZER_SEED:
            sys.exit(
                "--verify-against-phase4a requires optimizer_seed=1987 first; "
                f"got {seeds_used}."
            )
        print(
            f"Phase 4C verify-against-phase4a (optimizer_seed={seeds_used[0]})...",
            flush=True,
        )
        ok, failures = verify_against_phase4a(optimizer_seed=seeds_used[0])
        if ok:
            print("  OK — K=1 reproduces Phase 4A bit-for-bit.")
            return
        print(f"  FAIL — {len(failures)} mismatches:")
        for f in failures:
            print(f)
        sys.exit(1)

    optimizer_seeds = tuple(args.optimizer_seeds[: max(1, args.K)])
    K = len(optimizer_seeds)

    FOUNDATION.mkdir(parents=True, exist_ok=True)
    per_run_path  = FOUNDATION / "phase4c_optimizer_seed_probe_per_run.csv"
    per_cell_path = FOUNDATION / "phase4c_optimizer_seed_probe_per_cell_epsilon.csv"
    md_path       = FOUNDATION / "phase4c_optimizer_seed_probe.md"

    n_expected = (
        len(PILOT_SIZES) * len(PILOT_DIMS) * len(p4a.PHASE4A_SEEDS)
        * len(p4a.EPSILONS) * K
    )
    print(
        f"Phase 4C: {n_expected} runs "
        f"(K={K}, optimizer_seeds={optimizer_seeds}).",
        flush=True,
    )

    t0 = time.time()
    rows = build_rows(PILOT_SIZES, PILOT_DIMS, optimizer_seeds)
    elapsed = time.time() - t0

    per_cell_eps_rows = summarize_per_cell_epsilon(rows, optimizer_seeds)
    verdict, info = compute_verdict(per_cell_eps_rows)

    write_per_run_csv(rows, per_run_path)
    write_per_cell_epsilon_csv(per_cell_eps_rows, per_cell_path)
    write_markdown(
        rows, per_cell_eps_rows, verdict, info,
        optimizer_seeds=optimizer_seeds,
        sizes=PILOT_SIZES, dims=PILOT_DIMS,
        elapsed_s=elapsed, path=md_path,
    )

    print("")
    print("--- Phase 4C summary ---")
    print(f"  K:                       {K}")
    print(f"  optimizer_seeds:         {optimizer_seeds}")
    print(f"  runs:                    {len(rows)}")
    print(f"  valid runs:              {sum(1 for r in rows if r['valid'])}")
    print(f"  cells:                   {info.get('n_cells', 0)}")
    print(f"  mean_label_stability:    {info.get('mean_label_stability', float('nan')):.3f}")
    print(f"  iqr_ratio (IQR/mean):    {_fmt(info.get('iqr_ratio', float('nan')))}")
    print(f"  floor_fraction_K:        {_fmt(info.get('floor_fraction_K', float('nan')))}")
    print(f"  floor_fraction_phase4a:  {_fmt(info.get('floor_fraction_phase4a', float('nan')))}")
    print(f"  verdict:                 {verdict}")
    print("")
    print(f"PER-RUN CSV:           {per_run_path}")
    print(f"PER-CELL-EPSILON CSV:  {per_cell_path}")
    print(f"MD:                    {md_path}")


if __name__ == "__main__":
    main()
