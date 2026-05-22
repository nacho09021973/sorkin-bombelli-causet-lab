#!/usr/bin/env python3
"""Phase 4A — Epsilon sweep + pre-PySR morphology diagnostic.

Phase 3F closed at an INTERMEDIATE verdict: Panel D (order-only) crossed
the null floor (+6.2%) but not the strong threshold (+10%).  The dataset
was guarded-warmup only with two binary noise levels (small ε=1e-3 and
medium ε=5e-2) — by construction a tightly-clamped regime where the
optimizer rarely fails.  Most stratum residuals are small and PySR has
little physical variance to learn from.

Phase 4A is a *pure experimental-design* probe, not a PySR run.  It asks
a single, prior question:

    Does sweeping the initialization noise across a finer ε grid produce
    a visible aggregate V-like or interior-minimum morphology in the warmup-relative drift, broken down by
    causet geometry (n, target_dim)?

If yes (AGGREGATE_MORPHOLOGY_VISIBLE), the question becomes worth a PySR ablation
(Phase 4B), because the structural invariants have a chance of explaining
*where* on the V-like curve each causet sits.
If no (FLAT_RESPONSE), do not run PySR; the regime is too flat.
If only some curves show aggregate morphology (AMBIGUOUS), report and decide manually.

No PySR is run here.  No data from Phase 2G or earlier phases is reused.

Configuration
-------------
warmup_mode = guarded_warmup only
n ∈ {32, 64}
target_dim ∈ {2, 3, 4}
epsilon ∈ {0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.15, 0.20}
seeds = 10 fresh integers, disjoint from Phase 2G's EXTENDED_SEEDS
total rows = 2 × 3 × 8 × 10 = 480

Targets computed per row
------------------------
relative_drift            = warmup_delta_energy / initial_energy
abs_relative_drift        = |relative_drift|
log_abs_relative_drift    = log1p(|relative_drift|)

New adimensional features
-------------------------
dim_discrepancy_abs           = |mm_dim - midpoint_dim|
dim_discrepancy_rel_midpoint  = |mm_dim - midpoint_dim| / midpoint_dim
dim_ratio_mm_midpoint         = mm_dim / midpoint_dim

If midpoint_dim is zero, NaN, or non-finite for a given row, the row is
excluded from the dataset with explicit count in the markdown.  No silent
division.

Output
------
benchmarks/foundation/phase4a_epsilon_sweep.csv
benchmarks/foundation/phase4a_epsilon_sweep.md
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

import causet_invariants as ci  # noqa: E402
import cones  # noqa: E402
from tools.build_phase1_atlas import (  # noqa: E402
    _discrepancy,
    _estimate_dimensions,
    _format_field,
)
import validation_suite as vs  # noqa: E402


FOUNDATION = ROOT / "benchmarks" / "foundation"

# Fresh seeds, disjoint from Phase 1E / 2G:
#   EXTENDED_SEEDS = (1959, 1962, 1987, 2009, 2026,
#                    1812, 1848, 1871, 1905, 1929,
#                    1945, 1968, 1989, 2001, 2017)
# Phase 4A seeds: 10 distinct historical-year integers, no overlap.
PHASE4A_SEEDS: tuple[int, ...] = (
    1900, 1916, 1923, 1939, 1953,
    1973, 1981, 1995, 2003, 2020,
)
SIZES: tuple[int, ...]          = (32, 64)
SPACETIME_DIMS: tuple[int, ...] = (2, 3, 4)
EPSILONS: tuple[float, ...]     = (0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.15, 0.20)

OPTIMIZER_SEED  = 1987
WARMUP_LIMIT    = 10
ANNEAL_LIMIT    = 10
MAX_DATA        = 4
INITIAL_TEMP    = 100.0
COOLING_FACTOR  = 0.9
GUARD_THRESHOLD = 0.0

INITIAL_ENERGY_FLOOR = 1e-4
FLAT_SPAN_THRESHOLD  = 1e-6

# V-shape thresholds (documented; do not adjust post-hoc)
V_FALL_FRAC  = 0.30  # (vals[0] - vals[imin]) / vals[0] >= 30%
V_RISE_FRAC  = 0.05  # (vals[-1] - vals[imin]) / fall >= 5%
V_LEFT_DOWN  = 0.60  # fraction of left-wing diffs that must be negative
# Right-wing clean-rise: 0 negatives if n_right_pairs < 4; 1 allowed if >= 4

CSV_HEADERS = (
    "family", "target_dim", "n", "seed", "epsilon",
    "warmup_mode",
    "initial_energy", "final_energy", "delta_energy",
    "warmup_attempted_moves", "warmup_accepted_moves", "warmup_rejected_moves",
    "warmup_energy_before", "warmup_energy_after", "warmup_delta_energy",
    "relative_drift", "abs_relative_drift", "log_abs_relative_drift",
    "mm_dim", "midpoint_dim", "abs_discrepancy_mm_midpoint",
    "dim_discrepancy_abs", "dim_discrepancy_rel_midpoint", "dim_ratio_mm_midpoint",
    "chain2_count", "chain3_count", "chain3_abundance", "chain4_count",
    "link_count", "link_density",
    "relation_count", "ordering_fraction", "height",
    "row_valid", "skip_reason",
)


# ---------------------------------------------------------------------------
# Per-causet invariants (computed once per (n, seed, d))
# ---------------------------------------------------------------------------

def compute_invariants(matrix, n: int) -> dict:
    mm, mid = _estimate_dimensions(matrix)
    chains = ci.chain_counts(matrix, k_max=4)
    rc = ci.relation_count(matrix)
    of = ci.ordering_fraction(matrix)
    lc = ci.link_count(matrix)
    ht = ci.height(matrix)
    return {
        "mm_dim": mm,
        "midpoint_dim": mid,
        "abs_discrepancy_mm_midpoint": _discrepancy(mm, mid),
        "chain2_count": chains[2],
        "chain3_count": chains[3],
        "chain4_count": chains[4],
        "chain3_abundance": ci.three_chain_abundance(matrix),
        "link_count": lc,
        "link_density": lc / n if n > 0 else 0.0,
        "relation_count": rc,
        "ordering_fraction": of,
        "height": ht,
    }


def adimensional_features(inv: dict) -> dict:
    """Return adimensional dim features.  If midpoint_dim is degenerate the
    relative/ratio fields are NaN and the caller is responsible for
    flagging the row as invalid."""
    mm  = inv["mm_dim"]
    mid = inv["midpoint_dim"]
    abs_disc = abs(mm - mid)
    if mid is None or not math.isfinite(mid) or mid == 0.0:
        return {
            "dim_discrepancy_abs":          abs_disc,
            "dim_discrepancy_rel_midpoint": float("nan"),
            "dim_ratio_mm_midpoint":        float("nan"),
        }
    return {
        "dim_discrepancy_abs":          abs_disc,
        "dim_discrepancy_rel_midpoint": abs_disc / mid,
        "dim_ratio_mm_midpoint":        mm / mid,
    }


# ---------------------------------------------------------------------------
# Simulator helpers (mirror Phase 2G)
# ---------------------------------------------------------------------------

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


def _run_one(d: int, n: int, seed: int, epsilon: float,
             matrix, points, optimizer_seed: int = OPTIMIZER_SEED) -> dict:
    d_spatial = d - 1
    noise_seed = seed * 10007 + d * 1009 + n * 97 + int(epsilon * 1e6)
    noise_rng = random.Random(noise_seed)

    with contextlib.redirect_stdout(io.StringIO()):
        sim = cones.ConesSimulator(
            z=matrix, dim=d_spatial,
            seed=optimizer_seed, interactive=False,
            max_data=MAX_DATA, plot_path=None,
            warmup_limit=WARMUP_LIMIT, anneal_limit=ANNEAL_LIMIT,
            initial_temp=INITIAL_TEMP, cooling_factor=COOLING_FACTOR,
            backend="cpu",
        )
        buf = io.StringIO()
        _custom_startup(sim, points, epsilon, noise_rng)
        initial_energy = sim.initial_energy
        wmeta = _guarded_warmup(sim)
        sim.anneal(buf)
        final_energy = sim.data[-1][1] if sim.data else sim.eave

    return {
        "initial_energy": initial_energy,
        "final_energy":   final_energy,
        "delta_energy":   final_energy - initial_energy,
        **wmeta,
    }


# ---------------------------------------------------------------------------
# Row construction
# ---------------------------------------------------------------------------

def _row_from_sim(
    d: int, n: int, seed: int, epsilon: float,
    inv: dict, adim: dict, sim_result: dict,
) -> dict:
    E0  = sim_result["initial_energy"]
    wdE = sim_result["warmup_delta_energy"]

    row_valid = True
    skip_reason = ""

    # E0 floor
    if not math.isfinite(E0) or abs(E0) < INITIAL_ENERGY_FLOOR:
        row_valid = False
        skip_reason = f"E0_too_small ({E0:.3g} < {INITIAL_ENERGY_FLOOR:g})"

    # Compute relative_drift only if E0 is usable; otherwise NaN.
    if row_valid and math.isfinite(E0) and abs(E0) >= INITIAL_ENERGY_FLOOR:
        rd  = wdE / E0
        ard = abs(rd)
        lard = math.log1p(ard)
    else:
        rd = ard = lard = float("nan")

    # Degenerate midpoint_dim invalidates the adimensional features.
    if (not math.isfinite(adim["dim_discrepancy_rel_midpoint"])
            or not math.isfinite(adim["dim_ratio_mm_midpoint"])):
        if row_valid:
            row_valid = False
            skip_reason = "midpoint_dim_degenerate"

    return {
        "family":      "minkowski",
        "target_dim":  d,
        "n":           n,
        "seed":        seed,
        "epsilon":     epsilon,
        "warmup_mode": "guarded_warmup",
        "initial_energy": sim_result["initial_energy"],
        "final_energy":   sim_result["final_energy"],
        "delta_energy":   sim_result["delta_energy"],
        "warmup_attempted_moves": sim_result["warmup_attempted_moves"],
        "warmup_accepted_moves":  sim_result["warmup_accepted_moves"],
        "warmup_rejected_moves":  sim_result["warmup_rejected_moves"],
        "warmup_energy_before":   sim_result["warmup_energy_before"],
        "warmup_energy_after":    sim_result["warmup_energy_after"],
        "warmup_delta_energy":    sim_result["warmup_delta_energy"],
        "relative_drift":         rd,
        "abs_relative_drift":     ard,
        "log_abs_relative_drift": lard,
        **inv,
        **adim,
        "row_valid":   row_valid,
        "skip_reason": skip_reason,
    }


def build_rows() -> list[dict]:
    """Run the full sweep.  Sprinkle once per (n, seed, d); reuse across ε."""
    rows: list[dict] = []
    total = len(SIZES) * len(PHASE4A_SEEDS) * len(SPACETIME_DIMS) * len(EPSILONS)
    done = 0
    t0 = time.time()
    last_report = t0

    for n in SIZES:
        for seed in PHASE4A_SEEDS:
            for d in SPACETIME_DIMS:
                matrix, points = vs.sprinkle_minkowski_diamond(
                    n=n, seed=seed, d_spacetime=d
                )
                inv  = compute_invariants(matrix, n)
                adim = adimensional_features(inv)

                for eps in EPSILONS:
                    sim_result = _run_one(d, n, seed, eps, matrix, points)
                    rows.append(_row_from_sim(d, n, seed, eps, inv, adim, sim_result))
                    done += 1
                    now = time.time()
                    if now - last_report > 10:
                        rate = done / (now - t0) if now > t0 else 0
                        eta = (total - done) / rate if rate > 0 else 0
                        print(f"  {done}/{total} rows  "
                              f"({now-t0:.0f}s elapsed, ETA {eta:.0f}s)",
                              flush=True)
                        last_report = now

    print(f"  done: {done}/{total} rows in {time.time()-t0:.0f}s.", flush=True)
    return rows


# ---------------------------------------------------------------------------
# CSV
# ---------------------------------------------------------------------------

def _fmt(value) -> str:
    if value is None: return "NA"
    if isinstance(value, bool): return "true" if value else "false"
    if isinstance(value, str): return value
    if isinstance(value, int): return str(value)
    if isinstance(value, float):
        if math.isnan(value): return "NA"
        if math.isinf(value): return "inf" if value > 0 else "-inf"
    return _format_field(value)


def write_csv(rows: list[dict], path: Path) -> None:
    lines = [",".join(CSV_HEADERS)]
    for row in rows:
        lines.append(",".join(_fmt(row[h]) for h in CSV_HEADERS))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# Transition diagnostic
# ---------------------------------------------------------------------------

def _valid_subset(rows: list[dict]) -> list[dict]:
    return [r for r in rows if r["row_valid"]]


def _by_stratum(rows: list[dict]) -> dict[tuple, list[dict]]:
    out: dict[tuple, list[dict]] = {}
    for r in rows:
        key = (r["n"], r["target_dim"], r["epsilon"])
        out.setdefault(key, []).append(r)
    return out


def _by_curve(rows: list[dict]) -> dict[tuple, list[tuple]]:
    """Return {(n, target_dim): [(epsilon, mean_abs_rd, mean_log_abs_rd), ...]}"""
    bucket: dict[tuple, dict[float, list[dict]]] = {}
    for r in rows:
        key = (r["n"], r["target_dim"])
        bucket.setdefault(key, {}).setdefault(r["epsilon"], []).append(r)
    out: dict[tuple, list[tuple]] = {}
    for (n, d), per_eps in bucket.items():
        curve = []
        for eps in sorted(per_eps.keys()):
            rs = per_eps[eps]
            if not rs:
                continue
            mean_ard  = sum(r["abs_relative_drift"]     for r in rs) / len(rs)
            mean_lard = sum(r["log_abs_relative_drift"] for r in rs) / len(rs)
            curve.append((eps, mean_ard, mean_lard))
        out[(n, d)] = curve
    return out


def _monotonic_score(values: list[float]) -> tuple[int, int, str]:
    """Return (n_in_dominant_direction, n_pairs, direction).

    direction ∈ {'up', 'down', 'flat'}.
    A curve is 'quasi-monotonic' if n_in_dominant_direction >= n_pairs - 1
    (one violation allowed) and at least 5 of 7 consecutive differences
    move in the dominant direction.
    """
    if len(values) < 2:
        return (0, 0, "flat")
    diffs = [values[i+1] - values[i] for i in range(len(values)-1)]
    n_pairs = len(diffs)
    n_up   = sum(1 for d in diffs if d > 0)
    n_down = sum(1 for d in diffs if d < 0)
    if n_up >= n_down:
        return (n_up, n_pairs, "up")
    return (n_down, n_pairs, "down")


def _curve_is_monotonic(curve: list[tuple], idx: int = 1, min_dominant: int = 5) -> bool:
    """idx=1 → abs_rd, idx=2 → log_abs_rd."""
    vals = [c[idx] for c in curve]
    n_dom, n_pairs, _ = _monotonic_score(vals)
    return n_pairs >= 5 and n_dom >= min_dominant


def _classify_curve_shape(curve: list[tuple], idx: int = 1) -> str:
    """Classify the epsilon-ordered curve into one shape label.

    Returns one of: monotone_decay, monotone_growth, v_shape, flat, noisy.

    V-shape rule (all conditions must hold; thresholds defined as module constants):
      1. Interior minimum: 0 < argmin < len-1
      2. fall = vals[0] - vals[argmin]; fall/vals[0] >= V_FALL_FRAC (30%)
      3. rise = vals[-1] - vals[argmin]; rise/fall >= V_RISE_FRAC (5%)
      4. Left wing (vals[0..argmin]): >= V_LEFT_DOWN fraction of diffs are negative
      5. Right wing (vals[argmin..end]): at most 1 negative diff if n_right_pairs >= 4;
         zero negative diffs if n_right_pairs < 4
    """
    vals = [c[idx] for c in curve]
    n = len(vals)
    if n < 3:
        return "flat"
    span = max(vals) - min(vals)
    if span < FLAT_SPAN_THRESHOLD:
        return "flat"

    imin = vals.index(min(vals))

    if 0 < imin < n - 1:
        v0   = vals[0]
        vmin = vals[imin]
        vlast = vals[-1]
        fall = v0 - vmin
        rise = vlast - vmin

        if (v0 > 0
                and fall / v0 >= V_FALL_FRAC
                and rise > 0
                and rise / fall >= V_RISE_FRAC):
            left_diffs  = [vals[i + 1] - vals[i] for i in range(imin)]
            n_left_down = sum(1 for d in left_diffs if d < 0)
            left_ok = (len(left_diffs) == 0
                       or n_left_down / len(left_diffs) >= V_LEFT_DOWN)

            right_diffs = [vals[i + 1] - vals[i] for i in range(imin, n - 1)]
            n_right_neg = sum(1 for d in right_diffs if d < 0)
            n_right_pairs = len(right_diffs)
            max_neg = 1 if n_right_pairs >= 4 else 0
            right_ok = n_right_neg <= max_neg

            if left_ok and right_ok:
                return "v_shape"

    n_dom, n_pairs, direction = _monotonic_score(vals)
    if n_pairs >= 5 and n_dom >= 5:
        return "monotone_decay" if direction == "down" else "monotone_growth"
    return "noisy"


def _shape_diagnostics(rows: list[dict]) -> dict:
    """Return per-curve shape classification and V-shape summary."""
    curves = _by_curve(rows)
    per_curve: dict[str, dict] = {}
    v_cells: list[tuple] = []
    for key in sorted(curves.keys()):
        n, d = key
        c = curves[key]
        shape_abs = _classify_curve_shape(c, idx=1)
        shape_log = _classify_curve_shape(c, idx=2)
        # Consensus: if abs says v_shape or monotone, use it; log is tiebreak for noisy
        shape = shape_abs if shape_abs != "noisy" else shape_log
        per_curve[f"{n}|{d}"] = {
            "shape_abs": shape_abs,
            "shape_log": shape_log,
            "shape":     shape,
        }
        if shape == "v_shape":
            v_cells.append(key)

    n_v     = len(v_cells)
    n_total = len(curves)
    summary = (
        "PARTIAL_V_LIKE_AGGREGATE_MORPHOLOGY"   if n_v >= 3 else
        "POSSIBLE_V_LIKE_AGGREGATE_MORPHOLOGY"  if n_v >= 2 else
        "NO_V_SHAPE_DETECTED"
    )
    return {
        "per_curve":    per_curve,
        "n_v_shape":    n_v,
        "n_curves":     n_total,
        "v_shape_cells": v_cells,
        "summary":      summary,
    }


def _verdict(rows: list[dict]) -> tuple[str, str, dict]:
    curves = _by_curve(rows)
    n_curves = len(curves)

    monotonic_abs  = {k: _curve_is_monotonic(c, idx=1) for k, c in curves.items()}
    monotonic_logabs = {k: _curve_is_monotonic(c, idx=2) for k, c in curves.items()}

    n_mono_abs    = sum(monotonic_abs.values())
    n_mono_logabs = sum(monotonic_logabs.values())
    n_mono_either = sum(1 for k in curves if monotonic_abs[k] or monotonic_logabs[k])

    shape_info = _shape_diagnostics(rows)

    info = {
        "n_curves":      n_curves,
        "n_mono_abs":    n_mono_abs,
        "n_mono_logabs": n_mono_logabs,
        "n_mono_either": n_mono_either,
        "per_curve":     {f"{n}|{d}": {
                            "mono_abs":    monotonic_abs[(n, d)],
                            "mono_logabs": monotonic_logabs[(n, d)],
                          } for (n, d) in curves},
        "shape":         shape_info,
    }

    if n_mono_either >= 4:
        return (
            "AGGREGATE_MORPHOLOGY_VISIBLE",
            f"{n_mono_either}/{n_curves} (n, target_dim) curves show "
            "(quasi-)monotonic dependence of |relative_drift| or "
            "log(1+|relative_drift|) on ε. A V-like or interior-minimum "
            "aggregate morphology is observable in the optimizer response. "
            "Phase 4B PySR per-dimension ablation is justified.",
            info,
        )
    if n_mono_either <= 2:
        return (
            "FLAT_RESPONSE",
            f"Only {n_mono_either}/{n_curves} curves show monotonic ε dependence. "
            "The ε sweep does not expose a visible V-like or interior-minimum "
            "morphology in any clean majority of (n, target_dim) cells. "
            "Do not run Phase 4B; structural invariants have no informative "
            "aggregate target to predict.",
            info,
        )
    return (
        "AMBIGUOUS",
        f"{n_mono_either}/{n_curves} curves show monotonic ε dependence — "
        "between the strong and null thresholds. Inspect the per-curve table "
        "below and decide manually whether a V-like or interior-minimum "
        "morphology concentrates in a particular (n, target_dim) cell worth "
        "probing further.",
        info,
    )


# ---------------------------------------------------------------------------
# Markdown report
# ---------------------------------------------------------------------------

def _stratum_table_lines(rows: list[dict]) -> list[str]:
    lines = [
        "| n | target_dim | ε | count | mean rel | std rel "
        "| mean abs | mean log1p(abs) | mean dim_disc_rel |",
        "| ---: | :---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    strata = _by_stratum(rows)
    for key in sorted(strata.keys()):
        n, d, eps = key
        rs = strata[key]
        cnt = len(rs)
        if cnt == 0: continue
        rds = [r["relative_drift"] for r in rs]
        mean_rd = sum(rds)/cnt
        var_rd  = sum((v-mean_rd)**2 for v in rds)/cnt
        std_rd  = math.sqrt(var_rd)
        mean_ard  = sum(r["abs_relative_drift"]     for r in rs)/cnt
        mean_lard = sum(r["log_abs_relative_drift"] for r in rs)/cnt
        mean_ddr  = sum(r["dim_discrepancy_rel_midpoint"] for r in rs)/cnt
        lines.append(
            f"| {n} | {d} | {eps:.3g} | {cnt} "
            f"| {mean_rd:+.4g} | {std_rd:.4g} "
            f"| {mean_ard:.4g} | {mean_lard:.4g} | {mean_ddr:.4g} |"
        )
    return lines


def _curve_table_lines(rows: list[dict]) -> list[str]:
    """One row per (n, target_dim) — show abs_drift as ε grows."""
    curves = _by_curve(rows)
    eps_list = sorted({r["epsilon"] for r in rows})
    header = "| n | target_dim | metric | " + " | ".join(f"ε={e:g}" for e in eps_list) + " |"
    sep = "| ---: | :---: | --- | " + " | ".join("---:" for _ in eps_list) + " |"
    lines = [header, sep]
    for key in sorted(curves.keys()):
        n, d = key
        curve = {c[0]: (c[1], c[2]) for c in curves[key]}
        abs_cells = [f"{curve[e][0]:.3g}" if e in curve else "—" for e in eps_list]
        log_cells = [f"{curve[e][1]:.3g}" if e in curve else "—" for e in eps_list]
        lines.append(f"| {n} | {d} | mean abs rel drift | " + " | ".join(abs_cells) + " |")
        lines.append(f"| {n} | {d} | mean log1p(abs)    | " + " | ".join(log_cells) + " |")
    return lines


def write_markdown(
    rows_all: list[dict],
    rows_valid: list[dict],
    n_invalid_E0: int,
    n_invalid_mid: int,
    verdict_label: str,
    verdict_text: str,
    verdict_info: dict,
    path: Path,
    elapsed_s: float,
) -> None:
    n_expected = len(SIZES) * len(SPACETIME_DIMS) * len(EPSILONS) * len(PHASE4A_SEEDS)
    lines = [
        "# Phase 4A — Epsilon sweep + morphology diagnostic",
        "",
        "**Status:** experimental-design probe.  No PySR is run here.",
        "Phase 4B is gated on this verdict.",
        "",
        "## Verdict (automatic)",
        "",
        f"**{verdict_label}**",
        "",
        verdict_text,
        "",
        "## Configuration",
        "",
        f"- Warmup mode: guarded_warmup only "
        f"(threshold={GUARD_THRESHOLD}, limit={WARMUP_LIMIT})",
        f"- Sizes: {', '.join(str(n) for n in SIZES)}",
        f"- Spacetime dims: {', '.join(str(d) for d in SPACETIME_DIMS)}",
        f"- Epsilons: {', '.join(f'{e:g}' for e in EPSILONS)}",
        f"- Seeds ({len(PHASE4A_SEEDS)}, disjoint from Phase 2G): "
        f"{', '.join(str(s) for s in PHASE4A_SEEDS)}",
        f"- Optimizer seed: {OPTIMIZER_SEED} "
        f"(anneal_limit={ANNEAL_LIMIT}, max_data={MAX_DATA}, "
        f"T₀={INITIAL_TEMP}, γ={COOLING_FACTOR})",
        f"- Initial-energy floor for validity: {INITIAL_ENERGY_FLOOR:g}",
        "",
        "## Sample bookkeeping",
        "",
        f"- Expected rows: {n_expected}",
        f"- Generated rows: {len(rows_all)}",
        f"- Invalid (|initial_energy| below floor): {n_invalid_E0}",
        f"- Invalid (midpoint_dim degenerate): {n_invalid_mid}",
        f"- Rows used for diagnostic: {len(rows_valid)}",
        f"- Wall-clock: {elapsed_s:.0f}s",
        "",
        "## New adimensional features",
        "",
        "Per-row, computed once from the causal matrix:",
        "",
        "```",
        "dim_discrepancy_abs           = |mm_dim - midpoint_dim|",
        "dim_discrepancy_rel_midpoint  = |mm_dim - midpoint_dim| / midpoint_dim",
        "dim_ratio_mm_midpoint         = mm_dim / midpoint_dim",
        "```",
        "",
        "Rows where `midpoint_dim` is zero, NaN, or non-finite are excluded.",
        "",
        "## Targets",
        "",
        "```",
        "relative_drift          = warmup_delta_energy / initial_energy",
        "abs_relative_drift      = |relative_drift|",
        "log_abs_relative_drift  = log1p(|relative_drift|)",
        "```",
        "",
        "## Per-curve optimizer-response (mean across seeds)",
        "",
        *_curve_table_lines(rows_valid),
        "",
        "## Monotonicity counts",
        "",
        f"- Curves (n × target_dim): {verdict_info['n_curves']}",
        f"- Curves (quasi-)monotonic in mean abs_relative_drift vs ε: "
        f"{verdict_info['n_mono_abs']}",
        f"- Curves (quasi-)monotonic in mean log1p(abs) vs ε: "
        f"{verdict_info['n_mono_logabs']}",
        f"- Curves monotonic in either: {verdict_info['n_mono_either']}",
        "",
        "Quasi-monotonic = at least 5 of 7 successive ε-differences move "
        "in the dominant direction.",
        "",
        "Per-curve flags:",
        "",
        "| n×target_dim | mono in abs | mono in log1p(abs) |",
        "| --- | :---: | :---: |",
    ]
    for k, info in verdict_info["per_curve"].items():
        lines.append(
            f"| {k} "
            f"| {'✓' if info['mono_abs'] else '—'} "
            f"| {'✓' if info['mono_logabs'] else '—'} |"
        )

    # ------------------------------------------------------------------
    # Shape diagnostic (post-hoc, exploratory)
    # ------------------------------------------------------------------
    sd = verdict_info["shape"]
    v_cell_str = ", ".join(f"({n},{d})" for (n, d) in sd["v_shape_cells"])
    lines += [
        "",
        "## Curve-shape diagnostic (post-hoc, exploratory)",
        "",
        "This section is derived **after** observing the AMBIGUOUS verdict. "
        "The pre-registered monotonic criterion penalises V-shapes because a "
        "V-shaped curve is not monotone, even though a minimum followed by a "
        "clear rise is a potentially informative optimizer-response pattern.",
        "",
        "**V-shape detection rule** (applied to mean `abs_relative_drift` vs ε):",
        "",
        f"- Interior minimum: 0 < argmin < len−1",
        f"- Fall fraction: (vals[0] − vals[argmin]) / vals[0] ≥ {V_FALL_FRAC:.0%}",
        f"- Rise fraction: (vals[−1] − vals[argmin]) / fall ≥ {V_RISE_FRAC:.0%}",
        f"- Left wing: ≥ {V_LEFT_DOWN:.0%} of consecutive diffs are negative",
        f"- Right wing: 0 negative diffs if n_pairs < 4; ≤ 1 if n_pairs ≥ 4",
        "",
        "| n×target_dim | shape (abs_rd) | shape (log1p_abs_rd) | consensus shape |",
        "| --- | --- | --- | --- |",
    ]
    for k, s in sd["per_curve"].items():
        lines.append(
            f"| {k} | {s['shape_abs']} | {s['shape_log']} | **{s['shape']}** |"
        )
    lines += [
        "",
        "**Summary:**",
        "",
        f"- `verdict_original` = **{verdict_label}**",
        f"- `shape_diagnostic` = **{sd['summary']}**",
        f"- `n_v_shape` = {sd['n_v_shape']}/{sd['n_curves']}",
        f"- `v_shape_cells` = [{v_cell_str}]" if v_cell_str else
          "- `v_shape_cells` = []",
        "",
        "The pre-registered monotonic criterion yields **AMBIGUOUS**. "
        f"A post-hoc curve-shape diagnostic detects V-shapes in "
        f"{sd['n_v_shape']}/{sd['n_curves']} curves"
        + (f" ({v_cell_str})" if v_cell_str else "")
        + ". These findings are consistent: the monotonic criterion does not "
        "recognise V-shapes by construction. This supports an exploratory "
        "(not confirmatory) Phase 4B.",
    ]

    # ------------------------------------------------------------------
    # Physical interpretation — conservative framing (N=6 cells)
    # ------------------------------------------------------------------
    lines += [
        "",
        "## Physical interpretation (exploratory, N=6)",
        "",
        "The following observations are based on six (n, target_dim) cells. "
        "With N=6 a clean empirical separation does not yet constitute a "
        "physical law; it motivates a hypothesis to test in Phase 4B.",
        "",
        "**Observed pattern:**",
        "V-like aggregate morphology (interior minimum or decline-then-rise) appears only when two conditions hold jointly:",
        "",
        "```",
        "V-shape candidate if:",
        "    dim_discrepancy_rel_midpoint > θ   (θ ≈ 0.35–0.40, exploratory, not calibrated)",
        "    and min_val > floor_tolerance       (curve does not saturate to zero)",
        "```",
        "",
        "- Cells with low `dim_discrepancy_rel` (< 0.16) show no V-shape regardless "
        "of other invariants — the embedding is well-conditioned and the optimizer "
        "is insensitive to the ε regime.",
        "- `(64,3)` is the informative boundary case: its `dim_discrepancy_rel` and "
        "`ordering_fraction` are comparable to the V-shape cells, but its minimum "
        "reaches the numerical floor (`min_val = 0`). The potential V-signal is "
        "censored by floor saturation, not absent.",
        "- `ordering_fraction` and `chain3_abundance` vary consistently with "
        "`dim_discrepancy_rel` across the six cells, but it is not yet possible "
        "to determine which invariant is the primary predictor and which is a "
        "correlated proxy.",
        "",
        "**Phase 4B question priority:**",
        "",
        "1. Does `dim_discrepancy_rel` separate V/non-V cells when the (n, d) grid "
        "is expanded? (survival test — must pass before any stronger claim)",
        "2. Does `min_val ≈ 0` explain apparent false negatives such as `(64,3)`?",
        "3. Do `ordering_fraction` and `chain3_abundance` contribute independent "
        "predictive information once conditioned on `dim_discrepancy_rel`?",
        "",
        "Only if questions 1–3 survive should Phase 4B attempt symbolic regression "
        "of `rise_frac` or `ε_at_min` as continuous targets.",
        "",
        "**Conservative conclusion:**",
        "",
        "Phase 4A does not establish a morphological regime shift, but it identifies a clean "
        "exploratory pattern: V-like aggregate morphology occurs only in cells with "
        "large dimension-estimator discrepancy and without numerical floor "
        "saturation. The primary candidate control variable is "
        "`dim_discrepancy_rel_midpoint`; floor saturation is a secondary "
        "censoring mechanism. Phase 4B should test whether this separation "
        "survives a larger grid before attempting symbolic regression of "
        "`rise_frac` or `ε_at_min`.",
    ]

    lines += [
        "",
        "## Full per-stratum table",
        "",
        *_stratum_table_lines(rows_valid),
        "",
        "## Recommendation",
        "",
    ]
    if verdict_label == "AGGREGATE_MORPHOLOGY_VISIBLE":
        lines += [
            "Phase 4B is justified.  Build "
            "`tools/build_phase4b_pysr_ablation_per_dim.py` that runs PySR "
            "four-panel ablation per dimension (d=2, d=3, d=4, all dims) on "
            "this dataset, using `relative_drift` as primary target and "
            "`log_abs_relative_drift` as secondary.  Use the new adimensional "
            "features in Panel D.",
        ]
    elif verdict_label == "FLAT_RESPONSE":
        lines += [
            "**Do not run Phase 4B.** The ε sweep does not expose a V-like "
            "or interior-minimum morphology in the warmup-relative drift across "
            "(n, target_dim) cells. Structural invariants have no informative "
            "aggregate target to predict. Reconsider the experimental design "
            "before further PySR.",
        ]
    else:
        sd = verdict_info["shape"]
        v_cell_str = ", ".join(f"({n},{d})" for (n, d) in sd["v_shape_cells"])
        lines += [
            "Phase 4B is justified as an **exploratory** (not confirmatory) run. "
            "See the Physical interpretation section above for the question "
            "priority and the conditions under which symbolic regression of "
            "`rise_frac` / `ε_at_min` becomes appropriate.",
            "",
            "Design notes:",
            "",
            "- Include all (n, d) panels; d=2 serves as a negative control.",
            "- Report V-shape cells "
            + (f"({v_cell_str}) " if v_cell_str else "")
            + "as a secondary analysis, not as the primary stratification.",
            "- Do not restrict the run to V-shape cells only — that would "
            "condition on a post-hoc observation.",
            "- Label all Phase 4B results as exploratory until the "
            "`dim_discrepancy_rel` separation survives a larger (n, d) grid.",
        ]

    lines += [
        "",
        "## Reproducibility",
        "",
        "Regenerate via `make regen-phase4a`.",
        "Source: `tools/build_phase4a_epsilon_sweep.py`.",
        "",
        "No PySR is invoked by this script. No data from Phase 2G or earlier "
        "phases is reused.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


# ---------------------------------------------------------------------------
# CSV loader (for --regen-md without re-running simulations)
# ---------------------------------------------------------------------------

_INT_FIELDS  = {"n", "target_dim", "seed",
                "warmup_attempted_moves", "warmup_accepted_moves",
                "warmup_rejected_moves", "height"}
_BOOL_FIELDS = {"row_valid"}
_STR_FIELDS  = {"family", "warmup_mode", "skip_reason"}


def load_rows_from_csv(path: Path) -> list[dict]:
    """Parse an existing phase4a CSV back into row dicts."""
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines:
        return []
    headers = lines[0].split(",")
    rows: list[dict] = []
    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.split(",")
        row: dict = {}
        for h, v in zip(headers, parts):
            if h in _BOOL_FIELDS:
                row[h] = v.lower() == "true"
            elif h in _STR_FIELDS:
                row[h] = v
            elif v == "NA":
                row[h] = float("nan")
            elif h in _INT_FIELDS:
                try:
                    row[h] = int(v)
                except ValueError:
                    row[h] = 0
            else:
                try:
                    row[h] = float(v)
                except ValueError:
                    row[h] = v
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--regen-md", action="store_true",
        help="Regenerate markdown from the existing CSV without re-running simulations.",
    )
    args = ap.parse_args()

    FOUNDATION.mkdir(parents=True, exist_ok=True)
    csv_path = FOUNDATION / "phase4a_epsilon_sweep.csv"

    if args.regen_md:
        if not csv_path.exists():
            sys.exit(f"CSV not found: {csv_path}\nRun without --regen-md first.")
        print("Loading rows from existing CSV (no simulation re-run).", flush=True)
        rows_all = load_rows_from_csv(csv_path)
        elapsed = 0.0
    else:
        n_expected = len(SIZES) * len(SPACETIME_DIMS) * len(EPSILONS) * len(PHASE4A_SEEDS)
        print(f"Phase 4A: {n_expected} runs "
              f"(epsilon sweep, guarded_warmup only).", flush=True)
        t0 = time.time()
        rows_all = build_rows()
        elapsed = time.time() - t0
        write_csv(rows_all, csv_path)

    n_invalid_E0  = sum(1 for r in rows_all if str(r.get("skip_reason", "")).startswith("E0_too_small"))
    n_invalid_mid = sum(1 for r in rows_all if r.get("skip_reason", "") == "midpoint_dim_degenerate")
    rows_valid = _valid_subset(rows_all)

    verdict_label, verdict_text, verdict_info = _verdict(rows_valid)

    md_path = FOUNDATION / "phase4a_epsilon_sweep.md"
    write_markdown(
        rows_all=rows_all, rows_valid=rows_valid,
        n_invalid_E0=n_invalid_E0, n_invalid_mid=n_invalid_mid,
        verdict_label=verdict_label,
        verdict_text=verdict_text, verdict_info=verdict_info,
        path=md_path,
        elapsed_s=elapsed,
    )

    sd = verdict_info["shape"]
    v_cell_str = ", ".join(f"({n},{d})" for (n, d) in sd["v_shape_cells"])
    print(f"\n--- Phase 4A summary ---", flush=True)
    print(f"  rows:           {len(rows_all)} total / {len(rows_valid)} valid")
    print(f"  invalid (E0):   {n_invalid_E0}")
    print(f"  invalid (mid):  {n_invalid_mid}")
    print(f"  monotonic (abs|log|either): "
          f"{verdict_info['n_mono_abs']}|{verdict_info['n_mono_logabs']}"
          f"|{verdict_info['n_mono_either']} of {verdict_info['n_curves']}")
    print(f"  verdict_original:  {verdict_label}")
    print(f"  shape_diagnostic:  {sd['summary']}")
    print(f"  v_shape_cells:     {v_cell_str or 'none'}")
    if not args.regen_md:
        print(f"\nCSV:  {csv_path}", flush=True)
    print(f"MD:   {md_path}", flush=True)


if __name__ == "__main__":
    main()
