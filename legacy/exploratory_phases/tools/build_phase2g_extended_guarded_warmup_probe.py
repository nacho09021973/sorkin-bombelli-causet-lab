#!/usr/bin/env python3
"""Phase 2G — extended guarded-warmup probe for Phase 3F.

Phase 2F established that guarded warmup preserves near-truth starts
without destroying the random_init exploration that legacy warmup
provides.  Phase 3E then showed that within (noise, n, target_dim,
warmup_mode) strata of the Phase 2F data, the order-only PySR panel
hovers around 15% loss reduction — borderline, with midpoint_dim
acting partly as a target_dim proxy.

Phase 2G is the data-expansion experiment that the 3E diagnostic
called for:

  - 15 seeds (instead of 3): five times more samples per stratum.
  - Sizes (32, 64, 128): broader range in n.
  - Spacetime dims (2, 3, 4): unchanged.
  - Init labels: truth_plus_small_noise, truth_plus_medium_noise only.
    truth is excluded (energy is identically zero, no drift to learn);
    random_init is excluded (the residual analysis is about near-truth
    preservation, not landscape exploration).
  - Warmup mode: guarded_warmup only.  Eliminates the legacy/guarded
    heteroscedasticity that contaminated Phase 3D.

Grid: 15 × 3 × 3 × 2 = 270 rows.

This is the same physical experiment as Phase 2F's guarded_warmup
branch, scaled up and reduced to the cases the Phase 3F residual
ablation cares about.  No changes to cones.py.

Output
------
benchmarks/foundation/phase2g_extended_guarded_warmup_probe.csv
benchmarks/foundation/phase2g_extended_guarded_warmup_probe.md
"""

from __future__ import annotations

import contextlib
import io
import math
import random
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import cones  # noqa: E402
from tools.build_phase1_atlas import _format_field  # noqa: E402
from tools.build_phase1e_extended_structural_atlas import (  # noqa: E402
    EXTENDED_SEEDS,
    SPACETIME_DIMS,
    FOUNDATION,
)
import validation_suite as vs  # noqa: E402


# Probe sizes are a subset of Phase 1E's atlas sizes.  n=128 was attempted
# but excluded from the probe because the pure-Python ConesSimulator
# update() is O(n^2) per accepted move and a single n=128 d=4 run exceeds
# 60 s in this codebase.  Phase 1E's atlas does cover n=128 invariants so
# downstream join logic in Phase 3F still resolves cleanly; only the
# dynamics rows are restricted to {32, 64} here.
PROBE_SIZES: tuple[int, ...] = (32, 64)


OPTIMIZER_SEED  = 1987
WARMUP_LIMIT    = 10
ANNEAL_LIMIT    = 10
MAX_DATA        = 4
INITIAL_TEMP    = 100.0
COOLING_FACTOR  = 0.9

NOISE_SMALL  = 1e-3
NOISE_MEDIUM = 5e-2
GUARD_THRESHOLD = 0.0

INIT_LABELS: tuple[tuple[str, float], ...] = (
    ("truth_plus_small_noise",  NOISE_SMALL),
    ("truth_plus_medium_noise", NOISE_MEDIUM),
)

CSV_HEADERS = (
    "family",
    "target_dim",
    "n",
    "seed",
    "init_label",
    "warmup_mode",
    "noise_epsilon",
    "paired_key",
    "initial_energy",
    "final_energy",
    "delta_energy",
    "initial_interval_rmse",
    "final_interval_rmse",
    "initial_distance_to_truth_rms",
    "final_distance_to_truth_rms",
    "improved_energy",
    "improved_interval_rmse",
    "preserved_near_truth",
    "warmup_attempted_moves",
    "warmup_accepted_moves",
    "warmup_rejected_moves",
    "warmup_energy_before",
    "warmup_energy_after",
    "warmup_delta_energy",
    "notes",
)


# ------------------------------------------------------------------
# Coordinate helpers (Phase 2D/2E/2F conventions, unchanged)
# ------------------------------------------------------------------

def _coord_distance_rms(
    coords: list[vs.Coord],
    truth_scaled: list[vs.Coord],
    n: int,
    d_spatial: int,
) -> float:
    sq_sum = 0.0
    for i in range(n):
        ci, ti = coords[i], truth_scaled[i]
        dr = ci[0] - ti[0]
        sq_sum += dr * dr
        for k in range(1, d_spatial + 1):
            dx = (ci[k] if k < len(ci) else 0.0) - (ti[k] if k < len(ti) else 0.0)
            sq_sum += dx * dx
    return math.sqrt(sq_sum / n)


def _coords_from_sim(sim: cones.ConesSimulator) -> list[vs.Coord]:
    return [(sim.rold[i], *sim.xold[i]) for i in range(sim.n)]


def _scale_truth(points, rave_target, n):
    rave_truth = sum(p[0] for p in points) / n
    if rave_truth <= 0.0:
        return list(points)
    scale = rave_target / rave_truth
    return [tuple(c * scale for c in p) for p in points]


def _paired_key(d, n, seed, label, eps) -> str:
    return f"{d}|{n}|{seed}|{label}|{eps}"


def _custom_startup(sim, points, noise_epsilon, noise_rng):
    n = sim.n
    for i in range(n):
        sim.change[i] = True
        t_i = points[i][0]
        if noise_epsilon > 0.0:
            t_i = max(1e-12, t_i + noise_epsilon * noise_rng.gauss(0, 1))
        sim.rnew[i] = t_i
        for k in range(sim.dim):
            x_ik = points[i][k + 1] if k + 1 < len(points[i]) else 0.0
            if noise_epsilon > 0.0:
                x_ik += noise_epsilon * noise_rng.gauss(0, 1)
            sim.xnew[i][k] = x_ik
    sim.rave = sum(sim.rnew) / n
    sim.energy()
    sim.update()
    sim.initial_energy = sim.energies[0]


def _guarded_warmup(sim):
    energy_before = sim.energies[0]
    attempted = accepted = rejected = 0
    for _ in range(WARMUP_LIMIT):
        if sim.energies[0] <= 0.0:
            break
        for i in range(sim.n):
            sim.change[i] = False
        rave_saved = sim.r
        attempted += 1
        sim.reconfigure()
        sim.energy()
        if sim.deltae <= GUARD_THRESHOLD:
            sim.update()
            accepted += 1
        else:
            for i in range(sim.n):
                sim.change[i] = False
            sim.rave = rave_saved
            sim.deltae = 0.0
            rejected += 1
    sim.statistics()
    sim.warmup_energy = sim.energies[0]
    return {
        "warmup_attempted_moves": attempted,
        "warmup_accepted_moves":  accepted,
        "warmup_rejected_moves":  rejected,
        "warmup_energy_before":   energy_before,
        "warmup_energy_after":    sim.energies[0],
        "warmup_delta_energy":    sim.energies[0] - energy_before,
    }


# ------------------------------------------------------------------
# Per-row computation
# ------------------------------------------------------------------

def _probe_row(d, n, seed, init_label, noise_epsilon) -> dict:
    d_spatial = d - 1
    matrix, points = vs.sprinkle_minkowski_diamond(n=n, seed=seed, d_spacetime=d)

    noise_seed = seed * 10007 + d * 1009 + n * 97
    noise_rng = random.Random(noise_seed)

    with contextlib.redirect_stdout(io.StringIO()):
        sim = cones.ConesSimulator(
            z=matrix, dim=d_spatial,
            seed=OPTIMIZER_SEED, interactive=False,
            max_data=MAX_DATA, plot_path=None,
            warmup_limit=WARMUP_LIMIT, anneal_limit=ANNEAL_LIMIT,
            initial_temp=INITIAL_TEMP, cooling_factor=COOLING_FACTOR,
            backend="cpu",
        )
        buf = io.StringIO()

        _custom_startup(sim, points, noise_epsilon, noise_rng)

        initial_energy = sim.initial_energy
        initial_coords = _coords_from_sim(sim)
        rave_initial = sum(c[0] for c in initial_coords) / n
        truth_scaled_initial = _scale_truth(points, rave_initial, n)
        initial_interval_rmse = vs.interval_rmse(initial_coords, truth_scaled_initial)
        initial_dist = _coord_distance_rms(initial_coords, truth_scaled_initial, n, d_spatial)

        wmeta = _guarded_warmup(sim)
        sim.anneal(buf)

        final_energy = sim.data[-1][1] if sim.data else sim.eave
        final_coords = _coords_from_sim(sim)
        rave_final = sum(c[0] for c in final_coords) / n
        truth_scaled_final = _scale_truth(points, rave_final, n)
        final_interval_rmse = vs.interval_rmse(final_coords, truth_scaled_final)
        final_dist = _coord_distance_rms(final_coords, truth_scaled_final, n, d_spatial)

    delta_energy = final_energy - initial_energy
    improved_energy = delta_energy < -1e-9
    improved_interval_rmse = final_interval_rmse < initial_interval_rmse - 1e-9
    preserved_near_truth = delta_energy <= 1e-6

    if improved_energy:
        note = f"guarded+anneal improved energy {initial_energy:.4g} -> {final_energy:.4g}"
    elif preserved_near_truth:
        note = f"guarded+anneal preserved near-truth: energy {initial_energy:.4g} -> {final_energy:.4g}"
    else:
        note = f"guarded+anneal destroyed near-truth: energy {initial_energy:.4g} -> {final_energy:.4g}"

    return {
        "family": "minkowski",
        "target_dim": d,
        "n": n,
        "seed": seed,
        "init_label": init_label,
        "warmup_mode": "guarded_warmup",
        "noise_epsilon": noise_epsilon,
        "paired_key": _paired_key(d, n, seed, init_label, noise_epsilon),
        "initial_energy": initial_energy,
        "final_energy": final_energy,
        "delta_energy": delta_energy,
        "initial_interval_rmse": initial_interval_rmse,
        "final_interval_rmse": final_interval_rmse,
        "initial_distance_to_truth_rms": initial_dist,
        "final_distance_to_truth_rms": final_dist,
        "improved_energy": improved_energy,
        "improved_interval_rmse": improved_interval_rmse,
        "preserved_near_truth": preserved_near_truth,
        **wmeta,
        "notes": note,
    }


def build_rows() -> list[dict]:
    rows: list[dict] = []
    total = len(SPACETIME_DIMS) * len(PROBE_SIZES) * len(EXTENDED_SEEDS) * len(INIT_LABELS)
    done = 0
    t_start = time.time()
    last_report = t_start
    for d in SPACETIME_DIMS:
        for n in PROBE_SIZES:
            for seed in EXTENDED_SEEDS:
                for label, eps in INIT_LABELS:
                    rows.append(_probe_row(d, n, seed, label, eps))
                    done += 1
                    now = time.time()
                    if now - last_report > 10:
                        elapsed = now - t_start
                        rate = done / elapsed
                        eta = (total - done) / rate if rate > 0 else 0
                        print(f"  {done}/{total} rows  "
                              f"({elapsed:.0f}s elapsed, ETA {eta:.0f}s)",
                              flush=True)
                        last_report = now
    print(f"  done: {done}/{total} in {time.time()-t_start:.0f}s.")
    return rows


# ------------------------------------------------------------------
# Output
# ------------------------------------------------------------------

def _fmt(value) -> str:
    if value is None: return "NA"
    if isinstance(value, bool): return "true" if value else "false"
    if isinstance(value, str): return value
    if isinstance(value, int): return str(value)
    return _format_field(value)


def write_csv(rows: list[dict], path: Path) -> None:
    lines = [",".join(CSV_HEADERS)]
    for row in rows:
        lines.append(",".join(_fmt(row[h]) for h in CSV_HEADERS))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_markdown(rows: list[dict], path: Path) -> None:
    by_d_n_eps: dict = {}
    for r in rows:
        key = (r["target_dim"], r["n"], r["init_label"])
        by_d_n_eps.setdefault(key, []).append(r)

    n_total = len(rows)
    n_preserved = sum(1 for r in rows if r["preserved_near_truth"])

    lines = [
        "# Phase 2G — Extended Guarded-Warmup Probe",
        "",
        "Data foundation for the Phase 3F PySR ablation.",
        "",
        f"- Total runs: {n_total}",
        f"- Preserved (delta_energy ≤ 1e-6): {n_preserved}",
        f"- Seeds ({len(EXTENDED_SEEDS)}): {', '.join(str(s) for s in EXTENDED_SEEDS)}",
        f"- Sizes: {', '.join(str(s) for s in PROBE_SIZES)}",
        f"- Spacetime dims: {', '.join(str(s) for s in SPACETIME_DIMS)}",
        f"- Init labels: {', '.join(l for l, _ in INIT_LABELS)}",
        "- Warmup mode: guarded_warmup only "
        f"(guard threshold = {GUARD_THRESHOLD}, warmup_limit = {WARMUP_LIMIT})",
        f"- Anneal: anneal_limit={ANNEAL_LIMIT}, max_data={MAX_DATA}, "
        f"T₀={INITIAL_TEMP}, γ={COOLING_FACTOR}, optimizer_seed={OPTIMIZER_SEED}",
        "",
        "## Preservation by (d, n, init)",
        "",
        "| d | n | init | runs | preserved | mean ΔE | mean ΔE/E₀ |",
        "| :---: | ---: | --- | ---: | ---: | ---: | ---: |",
    ]
    for key in sorted(by_d_n_eps.keys()):
        d, n, init = key
        rs = by_d_n_eps[key]
        n_p = sum(1 for r in rs if r["preserved_near_truth"])
        mean_de = sum(r["delta_energy"] for r in rs) / len(rs)
        # mean |ΔE/E₀| skipping tiny-E0 rows
        rs_ratio = [r["delta_energy"]/r["initial_energy"] for r in rs
                    if abs(r["initial_energy"]) > 1e-12]
        mean_ratio = sum(rs_ratio)/len(rs_ratio) if rs_ratio else float("nan")
        lines.append(
            f"| {d} | {n} | {init} | {len(rs)} | {n_p}/{len(rs)} | "
            f"{_format_field(mean_de)} | {_format_field(mean_ratio)} |"
        )

    lines += [
        "",
        "## Warmup statistics (guarded-warmup only)",
        "",
        "| d | n | init | mean attempted | mean accepted | mean rejected | mean W ΔE |",
        "| :---: | ---: | --- | ---: | ---: | ---: | ---: |",
    ]
    for key in sorted(by_d_n_eps.keys()):
        d, n, init = key
        rs = by_d_n_eps[key]
        att = sum(r["warmup_attempted_moves"] for r in rs) / len(rs)
        acc = sum(r["warmup_accepted_moves"]  for r in rs) / len(rs)
        rej = sum(r["warmup_rejected_moves"]  for r in rs) / len(rs)
        wde = sum(r["warmup_delta_energy"]    for r in rs) / len(rs)
        lines.append(
            f"| {d} | {n} | {init} | {att:.2f} | {acc:.2f} | {rej:.2f} | "
            f"{_format_field(wde)} |"
        )

    lines += [
        "",
        "Regenerate via `make regen-phase2g`.",
        "Source: `tools/build_phase2g_extended_guarded_warmup_probe.py`.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    n_cells = len(EXTENDED_SEEDS) * len(PROBE_SIZES) * len(SPACETIME_DIMS) * len(INIT_LABELS)
    print(f"Phase 2G: {n_cells} guarded-warmup probe runs.")
    rows = build_rows()

    FOUNDATION.mkdir(parents=True, exist_ok=True)
    write_csv(rows, FOUNDATION / "phase2g_extended_guarded_warmup_probe.csv")
    write_markdown(rows, FOUNDATION / "phase2g_extended_guarded_warmup_probe.md")

    n_preserved = sum(1 for r in rows if r["preserved_near_truth"])
    print(f"Wrote {len(rows)} rows ({n_preserved} preserved) to {FOUNDATION}.")


if __name__ == "__main__":
    main()
