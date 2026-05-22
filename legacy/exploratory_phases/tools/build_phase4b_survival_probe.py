#!/usr/bin/env python3
"""Phase 4B — Exploratory survival probe for the Phase 4A pattern.

Phase 4B does not run PySR and does not promote the Phase 4A post-hoc
pattern into a law.  It asks whether the empirical separation observed in
Phase 4A survives a larger (n, target_dim) grid:

    a V-like aggregate morphology appeared in Phase 4A/4B only when
    dim_discrepancy_rel_midpoint was high and the curve was not censored by
    floor saturation.

Outputs are aggregated one row per (n, target_dim) curve:

    benchmarks/foundation/phase4b_survival_probe.csv
    benchmarks/foundation/phase4b_survival_probe_per_epsilon.csv
    benchmarks/foundation/phase4b_survival_probe_per_seed.csv
    benchmarks/foundation/phase4b_survival_probe.md

Implementation note
-------------------
This script imports Phase 4A simulation and curve-shape helpers rather
than refactoring Phase 4A.  Phase 4A is a historical artifact with its own
fixed narrative, constants, and output paths; keeping Phase 4B separate
preserves that boundary while reusing the validated experimental protocol.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import validation_suite as vs  # noqa: E402
from tools import build_phase4a_epsilon_sweep as p4a  # noqa: E402


FOUNDATION = ROOT / "benchmarks" / "foundation"

PILOT_SIZES: tuple[int, ...] = (32, 48, 64)
PILOT_DIMS: tuple[int, ...] = (2, 3, 4)
FULL_SIZES: tuple[int, ...] = (32, 48, 64, 96)
FULL_DIMS: tuple[int, ...] = (2, 3, 4, 5)

DEFAULT_THETA = 0.35
DEFAULT_FLOOR_TOLERANCE = 1e-6

CSV_HEADERS = (
    "phase",
    "grid",
    "family",
    "target_dim",
    "n",
    "curve_shape",
    "survival_label",
    "phase4a_reference_cell",
    "high_dim_discrepancy",
    "floor_saturated",
    "theta",
    "floor_tolerance",
    "epsilon_at_min",
    "min_val",
    "rise_frac",
    "fall_frac",
    "dim_discrepancy_rel_midpoint",
    "ordering_fraction",
    "chain3_abundance",
    "has_interior_minimum",
    "tail_positive_count",
    "tail_negative_count",
    "tail_zero_count",
    "tail_n_pairs",
    "tail_positive_fraction",
    "tail_pattern",
    "rise_frac_margin",
    "borderline_v_like",
    "valid_rows",
    "invalid_rows",
    "valid_seed_count",
    "epsilon_count",
)

PER_EPSILON_CSV_HEADERS = (
    "phase",
    "grid",
    "family",
    "target_dim",
    "n",
    "epsilon",
    "mean_loss",
    "std_loss",
    "n_valid",
    "n_invalid",
    "delta_from_prev_epsilon",
    "delta_sign",
    "is_min_epsilon",
    "near_floor",
    "floor_saturated_cell",
    "has_interior_minimum_cell",
    "tail_positive_count_cell",
    "tail_negative_count_cell",
    "tail_zero_count_cell",
    "tail_n_pairs_cell",
    "tail_positive_fraction_cell",
    "tail_pattern_cell",
    "rise_frac_margin_cell",
    "borderline_v_like_cell",
    "curve_shape",
    "survival_label",
)

PER_SEED_CSV_HEADERS = (
    "phase",
    "grid",
    "family",
    "target_dim",
    "n",
    "epsilon",
    "seed",
    "loss",
    "valid",
    "failure_mode",
    "ordering_fraction",
    "chain3_abundance",
    "dim_discrepancy_rel_midpoint",
    "curve_shape_cell",
    "survival_label_cell",
    "mean_loss_cell_epsilon",
    "std_loss_cell_epsilon",
    "delta_from_prev_epsilon_cell",
    "is_min_epsilon_cell",
    "borderline_v_like_cell",
)

PHASE4A_V_CELLS = {(32, 3), (32, 4), (64, 4)}
PHASE4A_NON_V_CELLS = {(32, 2), (64, 2), (64, 3)}


def _grid_values(grid: str) -> tuple[tuple[int, ...], tuple[int, ...]]:
    if grid == "pilot":
        return PILOT_SIZES, PILOT_DIMS
    if grid == "full":
        return FULL_SIZES, FULL_DIMS
    raise ValueError(f"unknown grid: {grid}")


def build_rows(grid: str) -> list[dict]:
    """Run the selected Phase 4B row-level simulation grid."""
    sizes, dims = _grid_values(grid)
    total = len(sizes) * len(dims) * len(p4a.PHASE4A_SEEDS) * len(p4a.EPSILONS)
    rows: list[dict] = []
    done = 0
    t0 = time.time()
    last_report = t0

    for n in sizes:
        for seed in p4a.PHASE4A_SEEDS:
            for d in dims:
                matrix, points = vs.sprinkle_minkowski_diamond(
                    n=n, seed=seed, d_spacetime=d
                )
                inv = p4a.compute_invariants(matrix, n)
                adim = p4a.adimensional_features(inv)
                for eps in p4a.EPSILONS:
                    sim_result = p4a._run_one(d, n, seed, eps, matrix, points)
                    rows.append(
                        p4a._row_from_sim(d, n, seed, eps, inv, adim, sim_result)
                    )
                    done += 1
                    now = time.time()
                    if now - last_report > 10:
                        rate = done / (now - t0) if now > t0 else 0.0
                        eta = (total - done) / rate if rate > 0 else 0.0
                        print(
                            f"  {done}/{total} rows "
                            f"({now - t0:.0f}s elapsed, ETA {eta:.0f}s)",
                            flush=True,
                        )
                        last_report = now
    print(f"  done: {done}/{total} rows in {time.time() - t0:.0f}s.", flush=True)
    return rows


def _mean(rows: list[dict], key: str) -> float:
    vals = [float(r[key]) for r in rows if math.isfinite(float(r[key]))]
    return sum(vals) / len(vals) if vals else float("nan")


def _curve_stats(curve: list[tuple]) -> dict:
    vals = [float(c[1]) for c in curve]
    eps = [float(c[0]) for c in curve]
    if not vals:
        return {
            "epsilon_at_min": float("nan"),
            "min_val": float("nan"),
            "rise_frac": float("nan"),
            "fall_frac": float("nan"),
        }
    imin = vals.index(min(vals))
    v0 = vals[0]
    vmin = vals[imin]
    vlast = vals[-1]
    fall = v0 - vmin
    rise = vlast - vmin
    return {
        "epsilon_at_min": eps[imin],
        "min_val": vmin,
        "rise_frac": rise / fall if fall > 0 else 0.0,
        "fall_frac": fall / v0 if v0 > 0 else 0.0,
    }


def _tail_audit(curve: list[tuple], curve_shape: str, floor_saturated: bool) -> dict:
    vals = [float(c[1]) for c in curve]
    n = len(vals)
    if n < 3:
        return {
            "has_interior_minimum": False,
            "tail_positive_count": 0,
            "tail_negative_count": 0,
            "tail_zero_count": 0,
            "tail_n_pairs": 0,
            "tail_positive_fraction": float("nan"),
            "tail_pattern": "",
            "rise_frac_margin": float("nan"),
            "borderline_v_like": False,
        }

    imin = vals.index(min(vals))
    has_interior = 0 < imin < n - 1
    if has_interior:
        tail_diffs = [vals[i + 1] - vals[i] for i in range(imin, n - 1)]
    else:
        tail_diffs = []

    signs = []
    for delta in tail_diffs:
        if delta > 0:
            signs.append("positive")
        elif delta < 0:
            signs.append("negative")
        else:
            signs.append("zero")

    tail_n = len(signs)
    tail_pos = signs.count("positive")
    tail_neg = signs.count("negative")
    tail_zero = signs.count("zero")
    tail_positive_fraction = tail_pos / tail_n if tail_n else float("nan")
    stats = _curve_stats(curve)
    rise_frac = float(stats["rise_frac"])
    rise_margin = rise_frac - p4a.V_RISE_FRAC if math.isfinite(rise_frac) else float("nan")
    max_tail_neg = 1 if tail_n >= 4 else 0
    tail_clean = tail_neg <= max_tail_neg
    margin_small = math.isfinite(rise_margin) and 0.0 <= rise_margin <= p4a.V_RISE_FRAC
    borderline = (
        has_interior
        and not floor_saturated
        and math.isfinite(rise_frac)
        and rise_frac > 0.0
        and (curve_shape != "v_shape" or not tail_clean or margin_small)
    )
    return {
        "has_interior_minimum": has_interior,
        "tail_positive_count": tail_pos,
        "tail_negative_count": tail_neg,
        "tail_zero_count": tail_zero,
        "tail_n_pairs": tail_n,
        "tail_positive_fraction": tail_positive_fraction,
        "tail_pattern": ",".join(signs),
        "rise_frac_margin": rise_margin,
        "borderline_v_like": borderline,
    }


def classify_survival_cell(
    target_dim: int,
    curve_shape: str,
    dim_discrepancy_rel_midpoint: float,
    min_val: float,
    theta: float = DEFAULT_THETA,
    floor_tolerance: float = DEFAULT_FLOOR_TOLERANCE,
) -> tuple[str, bool, bool]:
    """Return (survival_label, high_dim_discrepancy, floor_saturated).

    The label is explicitly exploratory.  Floor-saturated non-V curves are
    censored before they can count as strong negative evidence.
    """
    high_disc = (
        math.isfinite(dim_discrepancy_rel_midpoint)
        and dim_discrepancy_rel_midpoint > theta
    )
    floor_saturated = math.isfinite(min_val) and min_val <= floor_tolerance
    is_v = curve_shape == "v_shape"

    if is_v and high_disc and not floor_saturated:
        return "supporting", high_disc, floor_saturated
    if is_v:
        return "counterexample", high_disc, floor_saturated
    if target_dim == 2:
        return "control_negative", high_disc, floor_saturated
    if high_disc and floor_saturated:
        return "censored_floor", high_disc, floor_saturated
    if high_disc and not floor_saturated:
        return "counterexample", high_disc, floor_saturated
    return "ambiguous", high_disc, floor_saturated


def summarize_curves(
    rows_all: list[dict],
    grid: str,
    theta: float = DEFAULT_THETA,
    floor_tolerance: float = DEFAULT_FLOOR_TOLERANCE,
) -> list[dict]:
    rows_valid = p4a._valid_subset(rows_all)
    curves = p4a._by_curve(rows_valid)
    summaries: list[dict] = []

    for key in sorted(curves):
        n, d = key
        curve = curves[key]
        shape_abs = p4a._classify_curve_shape(curve, idx=1)
        shape_log = p4a._classify_curve_shape(curve, idx=2)
        curve_shape = shape_abs if shape_abs != "noisy" else shape_log
        stats = _curve_stats(curve)
        cell_valid = [r for r in rows_valid if r["n"] == n and r["target_dim"] == d]
        cell_all = [r for r in rows_all if r["n"] == n and r["target_dim"] == d]
        dim_disc = _mean(cell_valid, "dim_discrepancy_rel_midpoint")
        label, high_disc, floor_sat = classify_survival_cell(
            target_dim=d,
            curve_shape=curve_shape,
            dim_discrepancy_rel_midpoint=dim_disc,
            min_val=stats["min_val"],
            theta=theta,
            floor_tolerance=floor_tolerance,
        )
        tail_audit = _tail_audit(curve, curve_shape, floor_sat)
        summaries.append({
            "phase": "phase4b_exploratory",
            "grid": grid,
            "family": "minkowski",
            "target_dim": d,
            "n": n,
            "curve_shape": curve_shape,
            "survival_label": label,
            "phase4a_reference_cell": (n, d) in PHASE4A_V_CELLS or (n, d) in PHASE4A_NON_V_CELLS,
            "high_dim_discrepancy": high_disc,
            "floor_saturated": floor_sat,
            "theta": theta,
            "floor_tolerance": floor_tolerance,
            "epsilon_at_min": stats["epsilon_at_min"],
            "min_val": stats["min_val"],
            "rise_frac": stats["rise_frac"],
            "fall_frac": stats["fall_frac"],
            "dim_discrepancy_rel_midpoint": dim_disc,
            "ordering_fraction": _mean(cell_valid, "ordering_fraction"),
            "chain3_abundance": _mean(cell_valid, "chain3_abundance"),
            **tail_audit,
            "valid_rows": len(cell_valid),
            "invalid_rows": len(cell_all) - len(cell_valid),
            "valid_seed_count": len({r["seed"] for r in cell_valid}),
            "epsilon_count": len({r["epsilon"] for r in cell_valid}),
        })
    return summaries


def phase4b_outcome(summary_rows: list[dict]) -> str:
    v_rows = [r for r in summary_rows if r["curve_shape"] == "v_shape"]
    v_supporting = [r for r in v_rows if r["survival_label"] == "supporting"]
    strong_false_positive_controls = [
        r for r in summary_rows
        if r["target_dim"] == 2 and r["curve_shape"] == "v_shape"
    ]
    counterexamples = [
        r for r in summary_rows
        if r["survival_label"] == "counterexample"
    ]

    majority_v_high_no_floor = bool(v_rows) and len(v_supporting) > len(v_rows) / 2
    controls_clean = not strong_false_positive_controls

    if majority_v_high_no_floor and controls_clean and not counterexamples:
        return "PASS_EXPLORATORY_SURVIVAL"
    if not majority_v_high_no_floor or not controls_clean:
        return "FAIL"
    return "MIXED"


def _fmt(value) -> str:
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, str):
        return value
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if math.isnan(value):
            return "NA"
        if math.isinf(value):
            return "inf" if value > 0 else "-inf"
        return f"{value:.10g}"
    return str(value)


def write_summary_csv(summary_rows: list[dict], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(CSV_HEADERS)
        for row in summary_rows:
            writer.writerow([_fmt(row[h]) for h in CSV_HEADERS])


def _std_population(values: list[float]) -> float:
    if not values:
        return float("nan")
    mean = sum(values) / len(values)
    var = sum((v - mean) ** 2 for v in values) / len(values)
    return math.sqrt(var)


def _delta_sign(delta: float) -> str:
    if not math.isfinite(delta):
        return "NA"
    if delta > 0:
        return "positive"
    if delta < 0:
        return "negative"
    return "zero"


def build_per_epsilon_rows(
    rows_all: list[dict],
    summary_rows: list[dict],
    grid: str,
) -> list[dict]:
    """Persist the per-epsilon curve values used by the classifier.

    mean_loss is mean(abs_relative_drift) across valid seeds for a given
    (n, target_dim, epsilon), matching Phase 4A's curve-shape classifier.
    std_loss is the population standard deviation across those valid rows.
    """
    summary_by_cell = {(r["n"], r["target_dim"]): r for r in summary_rows}
    out: list[dict] = []

    for key in sorted(summary_by_cell):
        n, d = key
        summary = summary_by_cell[key]
        previous = float("nan")
        epsilons = sorted({
            float(r["epsilon"])
            for r in rows_all
            if r["n"] == n and r["target_dim"] == d
        })
        for eps in epsilons:
            cell_eps_all = [
                r for r in rows_all
                if r["n"] == n and r["target_dim"] == d and float(r["epsilon"]) == eps
            ]
            cell_eps_valid = [r for r in cell_eps_all if r["row_valid"]]
            losses = [float(r["abs_relative_drift"]) for r in cell_eps_valid]
            mean_loss = sum(losses) / len(losses) if losses else float("nan")
            std_loss = _std_population(losses)
            delta = (
                mean_loss - previous
                if math.isfinite(mean_loss) and math.isfinite(previous)
                else float("nan")
            )
            out.append({
                "phase": "phase4b_exploratory",
                "grid": grid,
                "family": "minkowski",
                "target_dim": d,
                "n": n,
                "epsilon": eps,
                "mean_loss": mean_loss,
                "std_loss": std_loss,
                "n_valid": len(cell_eps_valid),
                "n_invalid": len(cell_eps_all) - len(cell_eps_valid),
                "delta_from_prev_epsilon": delta,
                "delta_sign": _delta_sign(delta),
                "is_min_epsilon": (
                    math.isfinite(mean_loss)
                    and math.isclose(eps, float(summary["epsilon_at_min"]), rel_tol=0.0, abs_tol=1e-12)
                ),
                "near_floor": (
                    math.isfinite(mean_loss)
                    and mean_loss <= float(summary["floor_tolerance"])
                ),
                "floor_saturated_cell": summary["floor_saturated"],
                "has_interior_minimum_cell": summary["has_interior_minimum"],
                "tail_positive_count_cell": summary["tail_positive_count"],
                "tail_negative_count_cell": summary["tail_negative_count"],
                "tail_zero_count_cell": summary["tail_zero_count"],
                "tail_n_pairs_cell": summary["tail_n_pairs"],
                "tail_positive_fraction_cell": summary["tail_positive_fraction"],
                "tail_pattern_cell": summary["tail_pattern"],
                "rise_frac_margin_cell": summary["rise_frac_margin"],
                "borderline_v_like_cell": summary["borderline_v_like"],
                "curve_shape": summary["curve_shape"],
                "survival_label": summary["survival_label"],
            })
            if math.isfinite(mean_loss):
                previous = mean_loss
    return out


def write_per_epsilon_csv(per_epsilon_rows: list[dict], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(PER_EPSILON_CSV_HEADERS)
        for row in per_epsilon_rows:
            writer.writerow([_fmt(row[h]) for h in PER_EPSILON_CSV_HEADERS])


def build_per_seed_rows(
    rows_all: list[dict],
    summary_rows: list[dict],
    per_epsilon_rows: list[dict],
    grid: str,
) -> list[dict]:
    """Persist one row per simulated Phase 4B seed/epsilon run.

    Per-cell and per-cell+epsilon columns are repeated deliberately so a
    later visual audit can select representative, lowest-loss and highest-loss runs
    without re-deriving joins from multiple files.
    """
    summary_by_cell = {(r["n"], r["target_dim"]): r for r in summary_rows}
    per_eps_by_cell = {
        (r["n"], r["target_dim"], r["epsilon"]): r
        for r in per_epsilon_rows
    }
    out: list[dict] = []
    for row in sorted(
        rows_all,
        key=lambda r: (r["n"], r["target_dim"], r["epsilon"], r["seed"]),
    ):
        n = row["n"]
        d = row["target_dim"]
        eps = row["epsilon"]
        summary = summary_by_cell[(n, d)]
        per_eps = per_eps_by_cell[(n, d, eps)]
        out.append({
            "phase": "phase4b_exploratory",
            "grid": grid,
            "family": row["family"],
            "target_dim": d,
            "n": n,
            "epsilon": eps,
            "seed": row["seed"],
            "loss": row["abs_relative_drift"],
            "valid": row["row_valid"],
            "failure_mode": row["skip_reason"],
            "ordering_fraction": row["ordering_fraction"],
            "chain3_abundance": row["chain3_abundance"],
            "dim_discrepancy_rel_midpoint": row["dim_discrepancy_rel_midpoint"],
            "curve_shape_cell": summary["curve_shape"],
            "survival_label_cell": summary["survival_label"],
            "mean_loss_cell_epsilon": per_eps["mean_loss"],
            "std_loss_cell_epsilon": per_eps["std_loss"],
            "delta_from_prev_epsilon_cell": per_eps["delta_from_prev_epsilon"],
            "is_min_epsilon_cell": per_eps["is_min_epsilon"],
            "borderline_v_like_cell": summary["borderline_v_like"],
        })
    return out


def write_per_seed_csv(per_seed_rows: list[dict], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(PER_SEED_CSV_HEADERS)
        for row in per_seed_rows:
            writer.writerow([_fmt(row[h]) for h in PER_SEED_CSV_HEADERS])


_INT_FIELDS = {
    "target_dim", "n", "valid_rows", "invalid_rows", "valid_seed_count",
    "epsilon_count", "tail_positive_count", "tail_negative_count",
    "tail_zero_count", "tail_n_pairs",
}
_BOOL_FIELDS = {
    "phase4a_reference_cell", "high_dim_discrepancy", "floor_saturated",
    "has_interior_minimum", "borderline_v_like",
}
_STR_FIELDS = {
    "phase", "grid", "family", "curve_shape", "survival_label", "tail_pattern",
}


def load_summary_csv(path: Path) -> list[dict]:
    rows: list[dict] = []
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            return []
        for raw in reader:
            row: dict = {}
            for h, v in raw.items():
                if v is None:
                    continue
                if h in _BOOL_FIELDS:
                    row[h] = v.lower() == "true"
                elif h in _INT_FIELDS:
                    row[h] = int(v)
                elif h in _STR_FIELDS:
                    row[h] = v
                elif v == "NA":
                    row[h] = float("nan")
                else:
                    row[h] = float(v)
            rows.append(row)
    return rows


def _md_table(summary_rows: list[dict]) -> list[str]:
    lines = [
        "| n | target_dim | curve_shape | survival_label | dim_disc_rel | min_val | epsilon_at_min | rise_frac | ordering_fraction | chain3_abundance |",
        "| ---: | :---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for r in sorted(summary_rows, key=lambda x: (x["n"], x["target_dim"])):
        lines.append(
            f"| {r['n']} | {r['target_dim']} | {r['curve_shape']} "
            f"| {r['survival_label']} "
            f"| {r['dim_discrepancy_rel_midpoint']:.4g} "
            f"| {r['min_val']:.4g} "
            f"| {r['epsilon_at_min']:.4g} "
            f"| {r['rise_frac']:.4g} "
            f"| {r['ordering_fraction']:.4g} "
            f"| {r['chain3_abundance']:.4g} |"
        )
    return lines


def _tail_audit_table(summary_rows: list[dict]) -> list[str]:
    interior = [
        r for r in sorted(summary_rows, key=lambda x: (x["n"], x["target_dim"]))
        if r["has_interior_minimum"]
    ]
    if not interior:
        return ["No cells have an interior minimum."]
    lines = [
        "| n | target_dim | curve_shape | survival_label | tail_pattern | rise_frac | rise_frac_margin | borderline_v_like |",
        "| ---: | :---: | --- | --- | --- | ---: | ---: | :---: |",
    ]
    for r in interior:
        lines.append(
            f"| {r['n']} | {r['target_dim']} | {r['curve_shape']} "
            f"| {r['survival_label']} | {r['tail_pattern']} "
            f"| {r['rise_frac']:.4g} | {r['rise_frac_margin']:.4g} "
            f"| {'true' if r['borderline_v_like'] else 'false'} |"
        )
    return lines


def write_markdown(
    summary_rows: list[dict],
    path: Path,
    grid: str,
    elapsed_s: float,
) -> None:
    sizes, dims = _grid_values(grid)
    outcome = phase4b_outcome(summary_rows)
    n_v = sum(1 for r in summary_rows if r["curve_shape"] == "v_shape")
    n_support = sum(1 for r in summary_rows if r["survival_label"] == "supporting")
    n_counter = sum(1 for r in summary_rows if r["survival_label"] == "counterexample")
    n_censored = sum(1 for r in summary_rows if r["survival_label"] == "censored_floor")
    controls = [r for r in summary_rows if r["target_dim"] == 2]
    censored = [r for r in summary_rows if r["floor_saturated"]]

    lines = [
        "# Phase 4B — Exploratory survival probe",
        "",
        "**Status:** exploratory survival test of the Phase 4A post-hoc pattern. No PySR is run here.",
        "",
        "## Objective",
        "",
        "Phase 4B asks whether a V-like aggregate optimizer-response morphology is concentrated in cells with high `dim_discrepancy_rel_midpoint` and without floor saturation. This is not a confirmatory physical claim.",
        "",
        "## Grid design",
        "",
        f"- Grid mode: `{grid}`",
        f"- Sizes: {', '.join(str(n) for n in sizes)}",
        f"- Target spacetime dims: {', '.join(str(d) for d in dims)}",
        f"- Epsilons: {', '.join(f'{e:g}' for e in p4a.EPSILONS)}",
        f"- Seeds: {', '.join(str(s) for s in p4a.PHASE4A_SEEDS)}",
        f"- Threshold theta: {summary_rows[0]['theta'] if summary_rows else DEFAULT_THETA:g}",
        f"- Floor tolerance: {summary_rows[0]['floor_tolerance'] if summary_rows else DEFAULT_FLOOR_TOLERANCE:g}",
        f"- Wall-clock: {elapsed_s:.0f}s"
        + (" (markdown regenerated from CSV; simulations not rerun)" if elapsed_s == 0.0 else ""),
        "",
        "The `pilot` grid is the default Makefile target. The `full` grid is available explicitly via `--grid full`.",
        "",
        "## Provenance note",
        "",
        "`phase4b_survival_probe.csv` is the one-row-per-cell aggregate used for the global exploratory outcome. `phase4b_survival_probe_per_epsilon.csv` persists the per-epsilon curve values (`mean_loss = mean(abs_relative_drift)`) used to audit curve morphology. `phase4b_survival_probe_per_seed.csv` persists one row per seed/epsilon run so later visual audits can select representative, lowest-loss and highest-loss runs before drawing Hasse diagrams. These files add provenance only; no thresholds, classifications, or physical conclusions are changed.",
        "",
        "## Curve-shape summary",
        "",
        *_md_table(summary_rows),
        "",
        "## Tail-cleanliness / borderline audit",
        "",
        "This post-hoc audit records whether a curve has an interior minimum, how cleanly the right tail rises, and whether the V-like behavior is marginal. It does not change `curve_shape`, `survival_label`, or the global outcome.",
        "",
        *_tail_audit_table(summary_rows),
        "",
        "## Survival test of Phase 4A hypothesis",
        "",
        f"- V-shape curves: {n_v}",
        f"- Supporting cells: {n_support}",
        f"- Counterexamples: {n_counter}",
        f"- Censored floor cases: {n_censored}",
        "",
        "`survival_label` is an exploratory interpretation, not a per-cell physical verdict. Floor-saturated high-discrepancy non-V curves are classified as `censored_floor` and are not counted as strong negative evidence.",
        "",
        "## Negative controls target_dim=2",
        "",
    ]
    if controls:
        lines += [
            "| n | curve_shape | survival_label | min_val | dim_disc_rel |",
            "| ---: | --- | --- | ---: | ---: |",
        ]
        for r in sorted(controls, key=lambda x: x["n"]):
            lines.append(
                f"| {r['n']} | {r['curve_shape']} | {r['survival_label']} "
                f"| {r['min_val']:.4g} | {r['dim_discrepancy_rel_midpoint']:.4g} |"
            )
    else:
        lines.append("No target_dim=2 controls are present in this grid.")

    lines += [
        "",
        "## Floor-saturated / censored cases",
        "",
    ]
    if censored:
        lines += [
            "| n | target_dim | curve_shape | survival_label | min_val | epsilon_at_min |",
            "| ---: | :---: | --- | --- | ---: | ---: |",
        ]
        for r in sorted(censored, key=lambda x: (x["n"], x["target_dim"])):
            lines.append(
                f"| {r['n']} | {r['target_dim']} | {r['curve_shape']} "
                f"| {r['survival_label']} | {r['min_val']:.4g} "
                f"| {r['epsilon_at_min']:.4g} |"
            )
    else:
        lines.append("No curves reached the configured floor tolerance.")

    lines += [
        "",
        "## Loss semantics",
        "",
        "In Phase 4B, `loss` is inherited directly from the Phase 4A `abs_relative_drift` column: `loss = |warmup_delta_energy / initial_energy|` for valid runs. Here `initial_energy` is the simulator energy after initializing the embedding from the sprinkled coordinates plus epsilon-scaled coordinate noise, and `warmup_delta_energy` is the change produced by the guarded warmup stage before the later anneal result is used.",
        "",
        "`loss` is therefore an optimizer/embedding response diagnostic under this specific energy, initialization, epsilon, seed, and coordinate parametrization. It is not an intrinsic observable of the partial order, not the dimensional discrepancy, not a Lorentz-invariant distance to Minkowski space, and not an absolute physical quality score for the causet.",
        "",
        "The visual audit roles `best`, `near_mean`, and `worst` rank seeds only by this Phase 4B loss within the same cell and epsilon. They are provenance labels for selecting examples, not physical labels for good or bad causal sets.",
        "",
        "## Global exploratory outcome",
        "",
        f"**{outcome}**",
        "",
        "Outcome definitions:",
        "",
        "- `PASS_EXPLORATORY_SURVIVAL`: most V-shapes occur in high-discrepancy, non-floor cells; target_dim=2 controls do not produce strong false positives; floor cases are censored rather than counted as clean negatives.",
        "- `MIXED`: partial signal with enough counterexamples or false positives to require refinement.",
        "- `FAIL`: the high-discrepancy plus no-floor separation does not survive the expanded grid.",
        "",
        "## Conservative conclusion",
        "",
    ]
    if outcome == "PASS_EXPLORATORY_SURVIVAL":
        lines.append(
            "Phase 4B supports survival of the exploratory Phase 4A pattern under this larger grid, but does not establish a physical law."
        )
    elif outcome == "MIXED":
        lines.append(
            "Phase 4B shows partial survival of the Phase 4A pattern, but counterexamples or ambiguous cells require refinement before any stronger claim."
        )
    else:
        lines.append(
            "Phase 4B does not support clean survival of the Phase 4A separation on this grid."
        )
    lines += [
        "",
        "Regenerate via `make regen-phase4b` for the pilot grid, or run `python3 tools/build_phase4b_survival_probe.py --grid full` explicitly for the full grid.",
        "Source: `tools/build_phase4b_survival_probe.py`.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", choices=("pilot", "full"), default="pilot")
    ap.add_argument("--regen-md", action="store_true")
    ap.add_argument("--theta", type=float, default=DEFAULT_THETA)
    ap.add_argument("--floor-tolerance", type=float, default=DEFAULT_FLOOR_TOLERANCE)
    args = ap.parse_args()

    FOUNDATION.mkdir(parents=True, exist_ok=True)
    csv_path = FOUNDATION / "phase4b_survival_probe.csv"
    per_epsilon_csv_path = FOUNDATION / "phase4b_survival_probe_per_epsilon.csv"
    per_seed_csv_path = FOUNDATION / "phase4b_survival_probe_per_seed.csv"
    md_path = FOUNDATION / "phase4b_survival_probe.md"

    if args.regen_md:
        if not csv_path.exists():
            sys.exit(f"CSV not found: {csv_path}\nRun without --regen-md first.")
        print("Loading Phase 4B summary from existing CSV (no simulation re-run).")
        summary_rows = load_summary_csv(csv_path)
        grid = summary_rows[0]["grid"] if summary_rows else args.grid
        elapsed = 0.0
    else:
        sizes, dims = _grid_values(args.grid)
        n_expected = len(sizes) * len(dims) * len(p4a.PHASE4A_SEEDS) * len(p4a.EPSILONS)
        print(f"Phase 4B: {n_expected} runs ({args.grid} grid).", flush=True)
        t0 = time.time()
        rows_all = build_rows(args.grid)
        elapsed = time.time() - t0
        summary_rows = summarize_curves(
            rows_all,
            grid=args.grid,
            theta=args.theta,
            floor_tolerance=args.floor_tolerance,
        )
        per_epsilon_rows = build_per_epsilon_rows(rows_all, summary_rows, args.grid)
        per_seed_rows = build_per_seed_rows(
            rows_all, summary_rows, per_epsilon_rows, args.grid
        )
        write_summary_csv(summary_rows, csv_path)
        write_per_epsilon_csv(per_epsilon_rows, per_epsilon_csv_path)
        write_per_seed_csv(per_seed_rows, per_seed_csv_path)
        grid = args.grid

    write_markdown(summary_rows, md_path, grid=grid, elapsed_s=elapsed)
    outcome = phase4b_outcome(summary_rows)
    print("\n--- Phase 4B summary ---")
    print(f"  grid:      {grid}")
    print(f"  cells:     {len(summary_rows)}")
    print(f"  outcome:   {outcome}")
    if not args.regen_md:
        print(f"CSV:  {csv_path}")
        print(f"Per-epsilon CSV:  {per_epsilon_csv_path}")
        print(f"Per-seed CSV:  {per_seed_csv_path}")
    print(f"MD:   {md_path}")


if __name__ == "__main__":
    main()
