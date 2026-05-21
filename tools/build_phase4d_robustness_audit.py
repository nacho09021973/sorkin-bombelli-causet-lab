#!/usr/bin/env python3
"""Phase 4D — Robustness-vs-invariants audit (no-PySR).

Phase 4C established that the Phase 4B `MIXED` / Phase 5 `INSUFFICIENT`
picture is `OPTIMIZER_SEED_LIMITED`: 5/9 (n, target_dim) cells flip their
curve-shape label across K=3 optimizer seeds, the per-cell IQR/loss ratio
is 0.92, and the floor-saturation fraction is invariant under optimizer
seed (0.474 ≈ 0.473).

Phase 4D asks the orthogonal question:

    When the target is marginalized over optimizer_seed, does any
    order-theoretic invariant *coexist* with the residual variance
    (instability, IQR_loss_K) or with the floor pathology?

This is a pre-PySR descriptive audit.  No PySR is run.  No new optimizer
is introduced.  No new simulations are run (Phase 4C per-run and
per-cell-epsilon CSVs are reused; only order-theoretic invariants are
recomputed via `p4a.compute_invariants` for each (n, d, causet_seed)).

Important semantic caveats
--------------------------
- `loss = |warmup_delta_energy / initial_energy|` is an optimizer/embedding
  response diagnostic of the (causet × pipeline) pair, NOT a physical
  observable of the partial order.
- `iqr_loss_K`, `per_seed_label_distinct_shapes_K`, and
  `floor_saturated_fraction_K` are properties of (causet × ε × schedule ×
  energy × optimizer family × K optimizer seeds), not of the causet alone.
- A non-zero Spearman correlation between an invariant and a robustness
  target means: "this invariant coexists with seed-instability under the
  current pipeline".  It does NOT mean: "this invariant predicts
  embeddability", "manifoldlike evidence", or "physical transition".

Targets of robustness
---------------------
- `per_seed_iqr_loss_K_mean_eps`     → robustness target (seed dispersion)
- `per_seed_label_distinct_shapes_K` → morphology robustness (1, 2 or 3)
- `per_seed_floor_saturated_fraction_K_mean_eps` → floor pathology target

`min_over_K` is explicitly excluded as a target: it is an optimistic
lucky-seed selector, not a robustness measure.

Verdict
-------
ORDER_THEORETIC_CORRELATE_DETECTED  if max |ρ_spearman| ≥ 0.6 at the
    per-seed level (N=90), AND at least one further invariant reaches
    |ρ_spearman| ≥ 0.6 against the same target with consistent sign.
WEAK_CORRELATE                      if 0.3 ≤ max |ρ_spearman| < 0.6 at
    the per-seed level.
NO_ROBUST_ORDER_THEORETIC_CORRELATE if max |ρ_spearman| < 0.3 at the
    per-seed level (across all 13 invariants × 3 targets).

If `NO_ROBUST_ORDER_THEORETIC_CORRELATE` triggers, the conservative
conclusion to record is:

    "No autonomous robust order-theoretic rule detected at this pipeline
    scale; observed Phase 3F/4B morphology is dominated by optimizer-seed
    variance and metric floor pathology."

Outputs
-------
benchmarks/foundation/phase4d_robustness_per_seed.csv      (90 rows)
benchmarks/foundation/phase4d_robustness_per_cell.csv      (9 rows)
benchmarks/foundation/phase4d_robustness_audit.md
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
PHASE4C_PER_RUN_CSV = FOUNDATION / "phase4c_optimizer_seed_probe_per_run.csv"
PHASE4C_PER_CELL_EPS_CSV = (
    FOUNDATION / "phase4c_optimizer_seed_probe_per_cell_epsilon.csv"
)
PHASE4D_PER_SEED_CSV = FOUNDATION / "phase4d_robustness_per_seed.csv"
PHASE4D_PER_CELL_CSV = FOUNDATION / "phase4d_robustness_per_cell.csv"
PHASE4D_MD          = FOUNDATION / "phase4d_robustness_audit.md"

DEFAULT_FLOOR_TOLERANCE = 1e-6

THRESHOLD_WEAK     = 0.30
THRESHOLD_DETECTED = 0.60

INVARIANTS = (
    "mm_dim",
    "midpoint_dim",
    "abs_discrepancy_mm_midpoint",
    "dim_discrepancy_abs",
    "dim_discrepancy_rel_midpoint",
    "dim_ratio_mm_midpoint",
    "ordering_fraction",
    "chain2_count",
    "chain3_count",
    "chain3_abundance",
    "chain4_count",
    "link_count",
    "link_density",
    "height",
    "relation_count",
)

PER_SEED_TARGETS = (
    "per_seed_iqr_loss_K_mean_eps",
    "per_seed_label_distinct_shapes_K",
    "per_seed_floor_saturated_fraction_K_mean_eps",
)

PER_CELL_TARGETS = (
    "cell_mean_iqr_loss_K_over_eps",
    "label_stability_cell",
    "cell_mean_floor_saturated_fraction_K_over_eps",
)


PER_CELL_HEADERS = (
    "phase",
    "n",
    "target_dim",
    "n_causet_seeds",
    "n_optimizer_seeds",
    "label_stability_cell",
    "curve_shape_per_optimizer_seed",
    "cell_mean_loss_K_over_eps",
    "cell_mean_iqr_loss_K_over_eps",
    "cell_max_iqr_loss_K_over_eps",
    "cell_iqr_ratio",
    "cell_mean_floor_saturated_fraction_K_over_eps",
    *(f"cell_mean_{inv}" for inv in INVARIANTS),
)


# ---------------------------------------------------------------------------
# n-control constants
# ---------------------------------------------------------------------------

N_CONTROL_METHOD = "rank-partial-pearson_residualize_n"
SIZE_LIKE_INVARIANTS = frozenset({
    "relation_count", "chain2_count", "chain3_count", "chain4_count", "height",
})
NCONTROL_N_STRATA = (32, 48, 64)


def _ncontrol_col_partial(inv: str, tgt: str) -> str:
    return f"nctrl_partial__{inv}__{tgt}"


def _ncontrol_col_stratum(n_val: int, inv: str, tgt: str) -> str:
    return f"nctrl_sn{n_val}__{inv}__{tgt}"


def _ncontrol_col_minabs(inv: str, tgt: str) -> str:
    return f"nctrl_minabs__{inv}__{tgt}"


def _ncontrol_pair_headers(inv: str, tgt: str) -> tuple[str, ...]:
    return (
        _ncontrol_col_partial(inv, tgt),
        *(_ncontrol_col_stratum(n, inv, tgt) for n in NCONTROL_N_STRATA),
        _ncontrol_col_minabs(inv, tgt),
    )


NCONTROL_PAIR_HEADERS: tuple[str, ...] = tuple(
    col
    for inv in INVARIANTS
    for tgt in PER_SEED_TARGETS
    for col in _ncontrol_pair_headers(inv, tgt)
)

PER_SEED_HEADERS = (
    "phase",
    "n",
    "target_dim",
    "causet_seed",
    "per_seed_iqr_loss_K_mean_eps",
    "per_seed_iqr_loss_K_max_eps",
    "per_seed_mean_loss_K_mean_eps",
    "per_seed_median_loss_K_mean_eps",
    "per_seed_min_loss_K_mean_eps",
    "per_seed_floor_saturated_fraction_K_mean_eps",
    "per_seed_label_distinct_shapes_K",
    "per_seed_curve_shape_per_optimizer_seed",
    *INVARIANTS,
    "nctrl_method",
    *NCONTROL_PAIR_HEADERS,
)


# ---------------------------------------------------------------------------
# Spearman / Pearson  (manual, no scipy dependency)
# ---------------------------------------------------------------------------

def _ranks_with_average_ties(values: list[float]) -> list[float]:
    n = len(values)
    sorted_pairs = sorted(enumerate(values), key=lambda p: p[1])
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and sorted_pairs[j + 1][1] == sorted_pairs[i][1]:
            j += 1
        avg_rank = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[sorted_pairs[k][0]] = avg_rank
        i = j + 1
    return ranks


def _pairwise_complete(xs: list, ys: list) -> tuple[list[float], list[float]]:
    out_x: list[float] = []
    out_y: list[float] = []
    for x, y in zip(xs, ys):
        if (
            x is None or y is None
            or (isinstance(x, float) and not math.isfinite(x))
            or (isinstance(y, float) and not math.isfinite(y))
        ):
            continue
        out_x.append(float(x))
        out_y.append(float(y))
    return out_x, out_y


def pearson_r(xs_raw: list, ys_raw: list) -> tuple[float, int]:
    xs, ys = _pairwise_complete(xs_raw, ys_raw)
    n = len(xs)
    if n < 2:
        return (float("nan"), n)
    mx = sum(xs) / n
    my = sum(ys) / n
    num = sum((xs[i] - mx) * (ys[i] - my) for i in range(n))
    den_x = math.sqrt(sum((xs[i] - mx) ** 2 for i in range(n)))
    den_y = math.sqrt(sum((ys[i] - my) ** 2 for i in range(n)))
    if den_x == 0.0 or den_y == 0.0:
        return (0.0, n)
    return (num / (den_x * den_y), n)


def spearman_rho(xs_raw: list, ys_raw: list) -> tuple[float, int]:
    xs, ys = _pairwise_complete(xs_raw, ys_raw)
    n = len(xs)
    if n < 2:
        return (float("nan"), n)
    rx = _ranks_with_average_ties(xs)
    ry = _ranks_with_average_ties(ys)
    return pearson_r(rx, ry)


# ---------------------------------------------------------------------------
# n-control: partial and stratified Spearman
# ---------------------------------------------------------------------------

def _residualize(values: list[float], covariate: list[float]) -> list[float]:
    """OLS-residualize values on covariate; returns residuals."""
    n = len(values)
    if n == 0:
        return []
    mean_c = sum(covariate) / n
    mean_v = sum(values) / n
    cov = sum((covariate[i] - mean_c) * (values[i] - mean_v) for i in range(n))
    var_c = sum((covariate[i] - mean_c) ** 2 for i in range(n))
    if var_c == 0.0:
        return [v - mean_v for v in values]
    beta = cov / var_c
    return [values[i] - mean_v - beta * (covariate[i] - mean_c) for i in range(n)]


def partial_spearman_rho_n(
    rows: list[dict],
    xs_col: str,
    ys_col: str,
    n_col: str = "n",
) -> tuple[float, str]:
    """Rank-partial Spearman of xs_col on ys_col controlling for n_col.

    Method: rank-transform x, y, n independently; residualize rank(x) and
    rank(y) on rank(n); return Pearson of the residuals.  This is the
    standard rank-partial-correlation approach.
    """
    xs_raw = [r[xs_col] for r in rows]
    ys_raw = [r[ys_col] for r in rows]
    ns_raw = [r[n_col] for r in rows]

    xs: list[float] = []
    ys: list[float] = []
    ns: list[float] = []
    for x, y, n_val in zip(xs_raw, ys_raw, ns_raw):
        if (
            x is None or y is None
            or (isinstance(x, float) and not math.isfinite(x))
            or (isinstance(y, float) and not math.isfinite(y))
        ):
            continue
        xs.append(float(x))
        ys.append(float(y))
        ns.append(float(n_val))

    if len(xs) < 3:
        return (float("nan"), N_CONTROL_METHOD)

    rx = _ranks_with_average_ties(xs)
    ry = _ranks_with_average_ties(ys)
    rn = _ranks_with_average_ties(ns)

    ex = _residualize(rx, rn)
    ey = _residualize(ry, rn)

    rho, _ = pearson_r(ex, ey)
    return (rho, N_CONTROL_METHOD)


def stratified_spearman_by_n(
    rows: list[dict],
    xs_col: str,
    ys_col: str,
    n_strata: tuple[int, ...],
    n_col: str = "n",
) -> tuple[dict[int, float], float]:
    """Spearman within each n-stratum.

    Returns ({n_val: rho}, min_abs_across_strata).
    """
    buckets: dict[int, tuple[list[float], list[float]]] = {
        n_val: ([], []) for n_val in n_strata
    }
    for row in rows:
        x = row[xs_col]
        y = row[ys_col]
        n_val = int(row[n_col])
        if n_val not in buckets:
            continue
        if (
            x is None or y is None
            or (isinstance(x, float) and not math.isfinite(x))
            or (isinstance(y, float) and not math.isfinite(y))
        ):
            continue
        buckets[n_val][0].append(float(x))
        buckets[n_val][1].append(float(y))

    result: dict[int, float] = {}
    for n_val in n_strata:
        xs, ys = buckets[n_val]
        if len(xs) < 2:
            result[n_val] = float("nan")
        else:
            rho, _ = spearman_rho(xs, ys)
            result[n_val] = rho

    finite_abs = [abs(v) for v in result.values() if math.isfinite(v)]
    min_abs = min(finite_abs) if finite_abs else float("nan")
    return (result, min_abs)


def compute_ncontrol_matrix(
    rows: list[dict],
    invariants: tuple[str, ...],
    targets: tuple[str, ...],
    n_strata: tuple[int, ...],
    n_col: str = "n",
) -> dict[tuple[str, str], dict]:
    """Compute partial and stratified Spearman for each (invariant, target) pair."""
    out: dict[tuple[str, str], dict] = {}
    for inv in invariants:
        for tgt in targets:
            partial_rho, method = partial_spearman_rho_n(rows, inv, tgt, n_col)
            strat_dict, min_abs = stratified_spearman_by_n(rows, inv, tgt, n_strata, n_col)
            entry: dict = {
                "partial_spearman_n": partial_rho,
                "method": method,
                "min_abs_stratified": min_abs,
            }
            for n_val in n_strata:
                entry[f"spearman_n{n_val}"] = strat_dict.get(n_val, float("nan"))
            out[(inv, tgt)] = entry
    return out


def inject_ncontrol_into_rows(
    rows: list[dict],
    ncontrol_corr: dict[tuple[str, str], dict],
) -> None:
    """Stamp dataset-level n-control statistics into each per-seed row."""
    for row in rows:
        row["nctrl_method"] = N_CONTROL_METHOD
        for (inv, tgt), stats in ncontrol_corr.items():
            row[_ncontrol_col_partial(inv, tgt)] = stats["partial_spearman_n"]
            for n_val in NCONTROL_N_STRATA:
                row[_ncontrol_col_stratum(n_val, inv, tgt)] = stats.get(
                    f"spearman_n{n_val}", float("nan")
                )
            row[_ncontrol_col_minabs(inv, tgt)] = stats["min_abs_stratified"]


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def _parse_bool(s: str) -> bool:
    return s.strip().lower() == "true"


def _parse_float_or_nan(s: str) -> float:
    if s == "" or s == "NA":
        return float("nan")
    try:
        return float(s)
    except ValueError:
        return float("nan")


def load_phase4c_per_run(path: Path = PHASE4C_PER_RUN_CSV) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def load_phase4c_per_cell_eps(path: Path = PHASE4C_PER_CELL_EPS_CSV) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


# ---------------------------------------------------------------------------
# Invariants per (n, target_dim, causet_seed)
# ---------------------------------------------------------------------------

def compute_invariants_for_cells(
    sizes: tuple[int, ...], dims: tuple[int, ...],
    causet_seeds: tuple[int, ...],
) -> dict[tuple[int, int, int], dict]:
    out: dict[tuple[int, int, int], dict] = {}
    for n in sizes:
        for d in dims:
            for seed in causet_seeds:
                matrix, _ = vs.sprinkle_minkowski_diamond(
                    n=n, seed=seed, d_spacetime=d,
                )
                inv = p4a.compute_invariants(matrix, n)
                adim = p4a.adimensional_features(inv)
                out[(n, d, seed)] = {**inv, **adim}
    return out


# ---------------------------------------------------------------------------
# Per-seed aggregation
# ---------------------------------------------------------------------------

def _safe_iqr(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    s = sorted(values)
    n = len(s)
    q1_idx = max(0, (n - 1) // 4)
    q3_idx = min(n - 1, (3 * (n - 1)) // 4)
    return s[q3_idx] - s[q1_idx]


def _median(values: list[float]) -> float:
    if not values:
        return float("nan")
    s = sorted(values)
    n = len(s)
    if n % 2 == 1:
        return s[n // 2]
    return (s[n // 2 - 1] + s[n // 2]) / 2.0


def build_per_seed_rows(
    per_run_rows: list[dict],
    invariants_lookup: dict[tuple[int, int, int], dict],
    optimizer_seeds: tuple[int, ...],
    floor_tolerance: float = DEFAULT_FLOOR_TOLERANCE,
) -> list[dict]:
    """One row per (n, target_dim, causet_seed).  Aggregates K runs at each
    epsilon into (mean/median/min/iqr/floor) then averages over the 8 epsilons.
    Also classifies the K-curve shapes for this individual causet."""
    keyed: dict[tuple[int, int, int, float, int], dict] = {}
    for r in per_run_rows:
        if not _parse_bool(r["valid"]):
            continue
        key = (
            int(r["n"]), int(r["target_dim"]), int(r["causet_seed"]),
            float(r["epsilon"]), int(r["optimizer_seed"]),
        )
        keyed[key] = r

    seed_cells: dict[tuple[int, int, int], dict[float, list[float]]] = {}
    for (n, d, cseed, eps, _opt), r in keyed.items():
        loss = _parse_float_or_nan(r["loss"])
        if not math.isfinite(loss):
            continue
        seed_cells.setdefault((n, d, cseed), {}).setdefault(eps, []).append(loss)

    shape_seeds: dict[tuple[int, int, int], dict[int, list[tuple[float, float]]]] = {}
    for (n, d, cseed, eps, opt), r in keyed.items():
        loss = _parse_float_or_nan(r["loss"])
        if not math.isfinite(loss):
            continue
        shape_seeds.setdefault((n, d, cseed), {}).setdefault(opt, []).append((eps, loss))

    out: list[dict] = []
    cell_keys = sorted(seed_cells.keys())
    for (n, d, cseed) in cell_keys:
        per_eps_losses = seed_cells[(n, d, cseed)]
        eps_iqrs: list[float] = []
        eps_means: list[float] = []
        eps_medians: list[float] = []
        eps_mins: list[float] = []
        eps_floors: list[float] = []
        for eps in sorted(per_eps_losses):
            losses = per_eps_losses[eps]
            eps_iqrs.append(_safe_iqr(losses))
            eps_means.append(sum(losses) / len(losses))
            eps_medians.append(_median(losses))
            eps_mins.append(min(losses))
            n_floor = sum(1 for v in losses if v <= floor_tolerance)
            eps_floors.append(n_floor / len(losses))

        opt_curves = shape_seeds.get((n, d, cseed), {})
        shapes: list[str] = []
        for opt in optimizer_seeds:
            pts = sorted(opt_curves.get(opt, []), key=lambda x: x[0])
            curve = [(e, v, v) for (e, v) in pts]
            if len(curve) >= 3:
                shapes.append(p4a._classify_curve_shape(curve, idx=1))
            else:
                shapes.append("insufficient")
        distinct = len(set(shapes))
        shapes_str = "|".join(shapes)

        inv = invariants_lookup.get((n, d, cseed), {})
        row = {
            "phase": "phase4d_robustness_audit",
            "n": n,
            "target_dim": d,
            "causet_seed": cseed,
            "per_seed_iqr_loss_K_mean_eps":
                sum(eps_iqrs) / len(eps_iqrs) if eps_iqrs else float("nan"),
            "per_seed_iqr_loss_K_max_eps":
                max(eps_iqrs) if eps_iqrs else float("nan"),
            "per_seed_mean_loss_K_mean_eps":
                sum(eps_means) / len(eps_means) if eps_means else float("nan"),
            "per_seed_median_loss_K_mean_eps":
                sum(eps_medians) / len(eps_medians) if eps_medians else float("nan"),
            "per_seed_min_loss_K_mean_eps":
                sum(eps_mins) / len(eps_mins) if eps_mins else float("nan"),
            "per_seed_floor_saturated_fraction_K_mean_eps":
                sum(eps_floors) / len(eps_floors) if eps_floors else float("nan"),
            "per_seed_label_distinct_shapes_K": distinct,
            "per_seed_curve_shape_per_optimizer_seed": shapes_str,
        }
        for inv_name in INVARIANTS:
            row[inv_name] = inv.get(inv_name, float("nan"))
        out.append(row)
    return out


# ---------------------------------------------------------------------------
# Per-cell aggregation
# ---------------------------------------------------------------------------

def build_per_cell_rows(
    per_seed_rows: list[dict],
    per_cell_eps_phase4c: list[dict],
    optimizer_seeds: tuple[int, ...],
) -> list[dict]:
    p4c_cell_eps: dict[tuple[int, int, float], dict] = {}
    for r in per_cell_eps_phase4c:
        key = (int(r["n"]), int(r["target_dim"]), float(r["epsilon"]))
        p4c_cell_eps[key] = r

    cells: dict[tuple[int, int], list[dict]] = {}
    for r in per_seed_rows:
        cells.setdefault((r["n"], r["target_dim"]), []).append(r)

    out: list[dict] = []
    for (n, d) in sorted(cells.keys()):
        seed_rows = cells[(n, d)]
        eps_rows_p4c = [
            r for r in per_cell_eps_phase4c
            if int(r["n"]) == n and int(r["target_dim"]) == d
        ]
        if eps_rows_p4c:
            label_stability_cell = _parse_float_or_nan(
                eps_rows_p4c[0]["label_stability_cell"]
            )
            curve_shape_per_opt = eps_rows_p4c[0]["curve_shape_per_optimizer_seed"]
            mean_loss_K_values = [
                _parse_float_or_nan(r["mean_loss_K"]) for r in eps_rows_p4c
            ]
            iqr_K_values = [
                _parse_float_or_nan(r["iqr_loss_K"]) for r in eps_rows_p4c
            ]
            floor_K_values = [
                _parse_float_or_nan(r["floor_saturated_fraction_K"]) for r in eps_rows_p4c
            ]
            mean_loss_K_finite = [v for v in mean_loss_K_values if math.isfinite(v)]
            iqr_finite        = [v for v in iqr_K_values if math.isfinite(v)]
            floor_finite      = [v for v in floor_K_values if math.isfinite(v)]
            cell_mean_loss = (
                sum(mean_loss_K_finite) / len(mean_loss_K_finite)
                if mean_loss_K_finite else float("nan")
            )
            cell_mean_iqr = (
                sum(iqr_finite) / len(iqr_finite) if iqr_finite else float("nan")
            )
            cell_max_iqr = max(iqr_finite) if iqr_finite else float("nan")
            cell_iqr_ratio = (
                cell_mean_iqr / cell_mean_loss
                if (math.isfinite(cell_mean_iqr) and math.isfinite(cell_mean_loss)
                    and cell_mean_loss > 0.0)
                else float("nan")
            )
            cell_mean_floor = (
                sum(floor_finite) / len(floor_finite)
                if floor_finite else float("nan")
            )
        else:
            label_stability_cell = float("nan")
            curve_shape_per_opt = ""
            cell_mean_loss = cell_mean_iqr = cell_max_iqr = float("nan")
            cell_iqr_ratio = cell_mean_floor = float("nan")

        row = {
            "phase": "phase4d_robustness_audit",
            "n": n,
            "target_dim": d,
            "n_causet_seeds": len(seed_rows),
            "n_optimizer_seeds": len(optimizer_seeds),
            "label_stability_cell": label_stability_cell,
            "curve_shape_per_optimizer_seed": curve_shape_per_opt,
            "cell_mean_loss_K_over_eps": cell_mean_loss,
            "cell_mean_iqr_loss_K_over_eps": cell_mean_iqr,
            "cell_max_iqr_loss_K_over_eps": cell_max_iqr,
            "cell_iqr_ratio": cell_iqr_ratio,
            "cell_mean_floor_saturated_fraction_K_over_eps": cell_mean_floor,
        }
        for inv_name in INVARIANTS:
            vals = [seed_row[inv_name] for seed_row in seed_rows]
            finite = [v for v in vals if isinstance(v, (int, float)) and math.isfinite(v)]
            row[f"cell_mean_{inv_name}"] = (
                sum(finite) / len(finite) if finite else float("nan")
            )
        out.append(row)
    return out


# ---------------------------------------------------------------------------
# Correlation matrix and verdict
# ---------------------------------------------------------------------------

def compute_correlation_matrix(
    rows: list[dict],
    invariants: tuple[str, ...],
    targets: tuple[str, ...],
) -> dict[tuple[str, str], dict]:
    """Return {(invariant, target): {'spearman': r, 'pearson': r, 'n': n}}."""
    out: dict[tuple[str, str], dict] = {}
    for inv in invariants:
        for tgt in targets:
            xs = [r[inv] for r in rows]
            ys = [r[tgt] for r in rows]
            rho, n = spearman_rho(xs, ys)
            r_p, _ = pearson_r(xs, ys)
            out[(inv, tgt)] = {"spearman": rho, "pearson": r_p, "n": n}
    return out


def compute_verdict(
    per_seed_corr: dict[tuple[str, str], dict],
    invariants: tuple[str, ...],
    targets: tuple[str, ...],
) -> tuple[str, dict]:
    """Apply the 0.30 / 0.60 rule on per-seed Spearman ρ.

    DETECTED:   max|ρ| ≥ 0.60  AND  ≥ 2 invariants reach |ρ| ≥ 0.60 against the
                same target with consistent sign.
    WEAK:       0.30 ≤ max|ρ| < 0.60.
    NO_ROBUST:  max|ρ| < 0.30.
    """
    abs_rhos: list[tuple[float, str, str, float]] = []  # (|rho|, inv, tgt, signed_rho)
    for (inv, tgt), entry in per_seed_corr.items():
        rho = entry["spearman"]
        if math.isfinite(rho):
            abs_rhos.append((abs(rho), inv, tgt, rho))
    if not abs_rhos:
        return ("NO_ROBUST_ORDER_THEORETIC_CORRELATE", {"reason": "no finite correlations"})

    abs_rhos.sort(reverse=True)
    top_abs, top_inv, top_tgt, top_rho = abs_rhos[0]

    detected = False
    detected_pair: tuple[str, str] | None = None
    if top_abs >= THRESHOLD_DETECTED:
        for (a, inv2, tgt2, rho2) in abs_rhos[1:]:
            if (
                tgt2 == top_tgt
                and inv2 != top_inv
                and a >= THRESHOLD_DETECTED
                and (rho2 * top_rho) > 0.0
            ):
                detected = True
                detected_pair = (top_inv, inv2)
                break

    info = {
        "top_invariant": top_inv,
        "top_target": top_tgt,
        "top_abs_spearman": top_abs,
        "top_signed_spearman": top_rho,
        "n_above_weak":     sum(1 for (a, *_rest) in abs_rhos if a >= THRESHOLD_WEAK),
        "n_above_detected": sum(1 for (a, *_rest) in abs_rhos if a >= THRESHOLD_DETECTED),
        "detected_pair": detected_pair,
    }

    if detected:
        return ("ORDER_THEORETIC_CORRELATE_DETECTED", info)
    if top_abs >= THRESHOLD_WEAK:
        return ("WEAK_CORRELATE", info)
    return ("NO_ROBUST_ORDER_THEORETIC_CORRELATE", info)


# ---------------------------------------------------------------------------
# Unstable-cell identification
# ---------------------------------------------------------------------------

def identify_unstable_cells(per_cell_rows: list[dict]) -> list[dict]:
    return [r for r in per_cell_rows if r["label_stability_cell"] == 0.0]


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


def write_per_seed_csv(rows: list[dict], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(PER_SEED_HEADERS)
        for row in rows:
            writer.writerow([_fmt(row[h]) for h in PER_SEED_HEADERS])


def write_per_cell_csv(rows: list[dict], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(PER_CELL_HEADERS)
        for row in rows:
            writer.writerow([_fmt(row[h]) for h in PER_CELL_HEADERS])


def _corr_table_lines(
    corr: dict[tuple[str, str], dict],
    invariants: tuple[str, ...],
    targets: tuple[str, ...],
) -> list[str]:
    header = "| invariant | " + " | ".join(
        f"ρ vs {t}" for t in targets
    ) + " |"
    sep = "| --- " + "| ---: " * len(targets) + "|"
    lines = [header, sep]
    for inv in invariants:
        cells = []
        for tgt in targets:
            entry = corr.get((inv, tgt), {})
            rho = entry.get("spearman", float("nan"))
            cells.append(_fmt(rho))
        lines.append(f"| `{inv}` | " + " | ".join(cells) + " |")
    return lines


def _unstable_cells_table(unstable: list[dict]) -> list[str]:
    if not unstable:
        return [
            "No (n, target_dim) cells have label_stability_cell == 0 in the "
            "current Phase 4C K=3 grid.",
        ]
    header = ("| n | target_dim | curve_shape_per_optimizer_seed "
              "| dim_disc_rel | ordering_fraction | chain3_abundance "
              "| link_density | height | mm_dim |")
    sep = "| ---: | :---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |"
    lines = [header, sep]
    for r in sorted(unstable, key=lambda x: (x["n"], x["target_dim"])):
        lines.append(
            f"| {r['n']} | {r['target_dim']} "
            f"| `{r['curve_shape_per_optimizer_seed']}` "
            f"| {_fmt(r['cell_mean_dim_discrepancy_rel_midpoint'])} "
            f"| {_fmt(r['cell_mean_ordering_fraction'])} "
            f"| {_fmt(r['cell_mean_chain3_abundance'])} "
            f"| {_fmt(r['cell_mean_link_density'])} "
            f"| {_fmt(r['cell_mean_height'])} "
            f"| {_fmt(r['cell_mean_mm_dim'])} |"
        )
    return lines


def _ncontrol_section_lines(
    per_seed_corr: dict,
    ncontrol_corr: dict,
    invariants: tuple[str, ...],
    targets: tuple[str, ...],
) -> list[str]:
    """Render the n-control interpretive section for the markdown."""
    lines: list[str] = [
        "## n-control (interpretive layer)",
        "",
        "The n-control is interpretive auditing, not mechanical "
        "reclassification. The Phase 4D verdict "
        "(`ORDER_THEORETIC_CORRELATE_DETECTED` or otherwise) is computed from "
        "raw Spearman as before. The partial-n and stratified-by-n values are "
        "reported here so the reader can judge whether the raw correlation "
        "reflects n-scaling rather than an order-theoretic correlate.",
        "",
    ]
    for tgt in targets:
        lines.append(f"### Target: `{tgt}`")
        lines.append("")
        # Identify top raw correlate for this target
        top_inv = max(
            invariants,
            key=lambda inv: abs(per_seed_corr.get((inv, tgt), {}).get("spearman", 0.0)
                                if math.isfinite(per_seed_corr.get((inv, tgt), {}).get("spearman", float("nan")))
                                else 0.0),
        )
        # Table header
        strata_headers = " | ".join(f"ρ_n{n}" for n in NCONTROL_N_STRATA)
        lines.append(
            f"| invariant | raw ρ | partial ρ (ctrl n) | min_abs stratified "
            f"| {strata_headers} |"
        )
        lines.append(
            "| --- | ---: | ---: | ---: | "
            + " | ".join("---:" for _ in NCONTROL_N_STRATA)
            + " |"
        )
        for inv in invariants:
            raw_entry = per_seed_corr.get((inv, tgt), {})
            raw_rho = raw_entry.get("spearman", float("nan"))
            nc = ncontrol_corr.get((inv, tgt), {})
            partial = nc.get("partial_spearman_n", float("nan"))
            min_abs = nc.get("min_abs_stratified", float("nan"))
            strata_vals = " | ".join(
                _fmt(nc.get(f"spearman_n{n}", float("nan")))
                for n in NCONTROL_N_STRATA
            )
            lines.append(
                f"| `{inv}` | {_fmt(raw_rho)} | {_fmt(partial)} "
                f"| {_fmt(min_abs)} | {strata_vals} |"
            )
        lines.append("")
        # Size-like caveat if top raw correlate for this target is extensive
        if top_inv in SIZE_LIKE_INVARIANTS:
            lines.append(
                "> Top raw correlate is an extensive size-like invariant. The raw "
                "correlation is expected to track `n` independently of any "
                "order-theoretic content. Refer to the partial-n and stratified-n "
                "columns to judge whether the correlation survives after the "
                "finite-size effect is removed."
            )
            lines.append("")
    return lines


def write_markdown(
    per_seed_rows: list[dict],
    per_cell_rows: list[dict],
    per_seed_corr: dict,
    per_cell_corr: dict,
    verdict: str,
    info: dict,
    optimizer_seeds: tuple[int, ...],
    path: Path,
    ncontrol_corr: dict | None = None,
) -> None:
    K = len(optimizer_seeds)
    n_seed_rows = len(per_seed_rows)
    n_cell_rows = len(per_cell_rows)
    unstable = identify_unstable_cells(per_cell_rows)
    lines = [
        "# Phase 4D — Robustness-vs-invariants audit (no-PySR)",
        "",
        "**Status:** descriptive audit.  No PySR.  No new simulations.  No "
        "modification to `cones.py`, Phase 4A/4B/4C/5 CSVs, thresholds, or "
        "verdicts.",
        "",
        "## Semantic caveats (read first)",
        "",
        "- `loss = |warmup_delta_energy / initial_energy|` is an optimizer/"
        "embedding response diagnostic of the (causet × pipeline) pair.  It "
        "is **not** a physical observable of the partial order, not a "
        "Lorentz-invariant residual, not a manifoldness score, and not an "
        "embeddability witness.",
        "- `iqr_loss_K`, `per_seed_label_distinct_shapes_K`, and "
        "`floor_saturated_fraction_K` are properties of (causet × ε × "
        "schedule × energy × optimizer family × K optimizer seeds).  They are "
        "**not** properties of the causet alone.",
        "- A non-zero correlation between an order-theoretic invariant and a "
        "robustness target means: \"this invariant coexists with seed "
        "instability under the current pipeline\".  It does **not** establish "
        "an embedding rule, a manifoldness signature, or any pipeline-"
        "independent geometric claim.",
        "",
        "## Objective",
        "",
        "Phase 4C showed that the Phase 4B `MIXED` / Phase 5 `INSUFFICIENT` "
        "picture is `OPTIMIZER_SEED_LIMITED` (5/9 cells flip curve-shape "
        "label across K=3 optimizer seeds; per-cell IQR/loss ratio ≈ 0.92; "
        "floor saturation invariant under seed at 0.474 ≈ 0.473).  Phase 4D "
        "asks whether any order-theoretic invariant **coexists** with this "
        "residual variance or with the floor pathology, under the current "
        "pipeline.",
        "",
        "## Inputs",
        "",
        f"- `{PHASE4C_PER_RUN_CSV.name}`: per-run optimizer-seed multi-start "
        "from Phase 4C.",
        f"- `{PHASE4C_PER_CELL_EPS_CSV.name}`: per-cell-epsilon aggregate "
        "from Phase 4C (label_stability_cell, mean/iqr/floor over K×seeds).",
        "- Order-theoretic invariants: recomputed per (n, target_dim, "
        "causet_seed) via `p4a.compute_invariants` (identical to Phase 1/4A).",
        "",
        "## Method",
        "",
        f"- K = {K} optimizer seeds: {', '.join(str(s) for s in optimizer_seeds)}",
        f"- Causet seeds (per cell): {len(p4a.PHASE4A_SEEDS)}",
        f"- Per-seed level (N = {n_seed_rows}): one row per (n, target_dim, "
        "causet_seed).  For each epsilon, compute IQR/mean/median/min/floor "
        "across K runs; then average over the 8 epsilons.  Shape per K "
        "optimizer seed (per causet) classified individually; the number of "
        "distinct shape labels across K is `per_seed_label_distinct_shapes_K` "
        "∈ {1, 2, 3}.",
        f"- Per-cell level (N = {n_cell_rows}): one row per (n, target_dim).  "
        "Invariants averaged across the 10 causet_seeds.  `label_stability_"
        "cell` and `curve_shape_per_optimizer_seed` propagated from Phase 4C.",
        "- Correlation: Spearman ρ and Pearson r (pairwise complete on "
        "finite values).  Spearman ρ is the primary statistic.",
        "",
        "## Targets",
        "",
        "Robustness targets (per-seed level):",
        "",
        "- `per_seed_iqr_loss_K_mean_eps` — seed-dispersion robustness",
        "- `per_seed_label_distinct_shapes_K` — morphology robustness "
        "(1, 2 or 3)",
        "- `per_seed_floor_saturated_fraction_K_mean_eps` — floor pathology",
        "",
        "`min_over_K` is explicitly excluded as a target: trivially "
        "decreases with K (optimistic lucky-seed selection).",
        "",
        "## Verdict",
        "",
        f"**{verdict}**",
        "",
        "Decision rule (per-seed Spearman ρ, N = "
        f"{n_seed_rows}):",
        "",
        f"- `ORDER_THEORETIC_CORRELATE_DETECTED` if max |ρ| ≥ "
        f"{THRESHOLD_DETECTED} AND a second invariant reaches |ρ| ≥ "
        f"{THRESHOLD_DETECTED} against the same target with the same sign.",
        f"- `WEAK_CORRELATE` if {THRESHOLD_WEAK} ≤ max |ρ| < "
        f"{THRESHOLD_DETECTED}.",
        f"- `NO_ROBUST_ORDER_THEORETIC_CORRELATE` if max |ρ| < "
        f"{THRESHOLD_WEAK}.",
        "",
        "Verdict inputs:",
        "",
        f"- Top invariant: `{info.get('top_invariant', 'NA')}`",
        f"- Top target: `{info.get('top_target', 'NA')}`",
        f"- Top |ρ_spearman|: {_fmt(info.get('top_abs_spearman', float('nan')))}",
        f"- Top signed ρ_spearman: {_fmt(info.get('top_signed_spearman', float('nan')))}",
        f"- Pairs with |ρ| ≥ {THRESHOLD_WEAK}: "
        f"{info.get('n_above_weak', 0)}",
        f"- Pairs with |ρ| ≥ {THRESHOLD_DETECTED}: "
        f"{info.get('n_above_detected', 0)}",
        f"- Detected pair: `{info.get('detected_pair', None)}`",
        "",
        "## Correlation matrix (per-seed level, signed Spearman ρ)",
        "",
        *_corr_table_lines(per_seed_corr, INVARIANTS, PER_SEED_TARGETS),
        "",
        *(_ncontrol_section_lines(per_seed_corr, ncontrol_corr, INVARIANTS, PER_SEED_TARGETS)
          if ncontrol_corr is not None else []),
        "## Correlation matrix (per-cell level, signed Spearman ρ)",
        "",
        "No n-control is computed at per-cell level (N=9 is too small for partial "
        "or stratified statistics; per-cell is reported as descriptive sanity only).",
        "",
        "Per-cell N = 9; treat as sanity / consistency check.",
        "",
        *_corr_table_lines(per_cell_corr, INVARIANTS, PER_CELL_TARGETS),
        "",
        "## Cells with label_stability_cell = 0 (Phase 4C)",
        "",
        *_unstable_cells_table(unstable),
        "",
        "## Interpretation",
        "",
    ]
    if verdict == "ORDER_THEORETIC_CORRELATE_DETECTED":
        lines += [
            "A non-trivial coexistence pattern was detected between at least "
            "two order-theoretic invariants and one robustness target, with "
            "consistent sign.  This is descriptive: it does NOT establish a "
            "physical claim, an embedding rule, or a manifoldness signature.  "
            "It does establish that the optimizer-seed instability observed "
            "in Phase 4C has structural coexistence with order-theoretic "
            "features at this pipeline scale.  The natural follow-up is "
            "Phase 4E or a re-scoped Phase 3G that uses the aggregated "
            "target with an explicit caveat on the multiple-comparisons "
            "risk at N = "
            f"{n_seed_rows}.",
        ]
    elif verdict == "WEAK_CORRELATE":
        lines += [
            "A single invariant reaches |ρ| ≥ "
            f"{THRESHOLD_WEAK} but no second invariant confirms the pattern.  "
            "This is consistent with marginal signal or with multiple-"
            "comparisons noise across 13 invariants × 3 targets.  The "
            "conservative reading is that no robust order-theoretic "
            "correlate has been established; before any follow-up, the K "
            "should be increased (Phase 4C extension to K = 5 or 7) so the "
            "robustness targets are less coarse.",
        ]
    else:
        lines += [
            "No autonomous robust order-theoretic rule detected at this "
            "pipeline scale; observed Phase 3F/4B morphology is dominated "
            "by optimizer-seed variance and metric floor pathology.",
            "",
            "Concretely: across all "
            f"{len(INVARIANTS)} order-theoretic invariants × "
            f"{len(PER_SEED_TARGETS)} robustness targets, the maximum "
            f"|Spearman ρ| at N = {n_seed_rows} is "
            f"{_fmt(info.get('top_abs_spearman', float('nan')))}, below the "
            f"{THRESHOLD_WEAK} weak-correlate threshold.",
            "",
            "Recommended next step: do not introduce new optimizer "
            "machinery (no reheating, no parallel tempering).  The "
            "next experimental lever is the target / metric "
            "definition itself (e.g., a Lorentz-invariant interval "
            "RMSE under guarded warmup), not the optimizer.  PySR on "
            "the current target family is not justified at this scale.",
        ]
    lines += [
        "",
        "## Scope",
        "",
        "- No PySR.  No new optimizer.  No reheating.  No parallel tempering.",
        "- No modification to Phase 4A/4B/4C/5 CSVs, thresholds, labels, or "
        "verdicts.",
        "- No physical claim, no manifoldness claim, no embeddability claim.",
        "- Targets are properties of (causet × pipeline), not of the causet "
        "alone.",
        "",
        "## Reproducibility",
        "",
        "Regenerate via `make regen-phase4d`.",
        "Source: `tools/build_phase4d_robustness_audit.py`.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--floor-tolerance", type=float, default=DEFAULT_FLOOR_TOLERANCE,
        help="Floor tolerance for per_seed_floor_saturated_fraction_K.",
    )
    args = ap.parse_args()

    for path in (PHASE4C_PER_RUN_CSV, PHASE4C_PER_CELL_EPS_CSV):
        if not path.exists():
            sys.exit(f"Missing Phase 4C input: {path}\nRun make regen-phase4c first.")

    print("Loading Phase 4C inputs...", flush=True)
    per_run_rows = load_phase4c_per_run()
    per_cell_eps_rows = load_phase4c_per_cell_eps()

    sizes = tuple(sorted({int(r["n"]) for r in per_run_rows}))
    dims  = tuple(sorted({int(r["target_dim"]) for r in per_run_rows}))
    causet_seeds = tuple(sorted({int(r["causet_seed"]) for r in per_run_rows}))
    optimizer_seeds = tuple(sorted({int(r["optimizer_seed"]) for r in per_run_rows}))

    print(
        f"  Phase 4C grid: sizes={sizes}, dims={dims}, "
        f"causet_seeds={len(causet_seeds)}, optimizer_seeds={optimizer_seeds}",
        flush=True,
    )

    print("Computing invariants for "
          f"{len(sizes) * len(dims) * len(causet_seeds)} (n, d, causet_seed) "
          "triples...", flush=True)
    t0 = time.time()
    invariants_lookup = compute_invariants_for_cells(sizes, dims, causet_seeds)
    print(f"  done in {time.time() - t0:.2f}s.", flush=True)

    per_seed_rows = build_per_seed_rows(
        per_run_rows, invariants_lookup, optimizer_seeds,
        floor_tolerance=args.floor_tolerance,
    )
    per_cell_rows = build_per_cell_rows(
        per_seed_rows, per_cell_eps_rows, optimizer_seeds,
    )

    per_seed_corr = compute_correlation_matrix(
        per_seed_rows, INVARIANTS, PER_SEED_TARGETS,
    )
    per_cell_corr = compute_correlation_matrix(
        per_cell_rows,
        tuple(f"cell_mean_{inv}" for inv in INVARIANTS),
        PER_CELL_TARGETS,
    )
    per_cell_corr_renamed = {
        (inv.replace("cell_mean_", ""), tgt): v
        for (inv, tgt), v in per_cell_corr.items()
    }

    verdict, info = compute_verdict(
        per_seed_corr, INVARIANTS, PER_SEED_TARGETS,
    )

    print("Computing n-control (partial + stratified Spearman)...", flush=True)
    ncontrol_corr = compute_ncontrol_matrix(
        per_seed_rows, INVARIANTS, PER_SEED_TARGETS, NCONTROL_N_STRATA,
    )
    inject_ncontrol_into_rows(per_seed_rows, ncontrol_corr)

    write_per_seed_csv(per_seed_rows, PHASE4D_PER_SEED_CSV)
    write_per_cell_csv(per_cell_rows, PHASE4D_PER_CELL_CSV)
    write_markdown(
        per_seed_rows, per_cell_rows,
        per_seed_corr, per_cell_corr_renamed,
        verdict, info, optimizer_seeds, PHASE4D_MD,
        ncontrol_corr=ncontrol_corr,
    )

    print("")
    print("--- Phase 4D summary ---")
    print(f"  per-seed rows:           {len(per_seed_rows)}")
    print(f"  per-cell rows:           {len(per_cell_rows)}")
    print(f"  top invariant:           {info.get('top_invariant', 'NA')}")
    print(f"  top target:              {info.get('top_target', 'NA')}")
    print(f"  top |ρ_spearman|:        "
          f"{_fmt(info.get('top_abs_spearman', float('nan')))}")
    print(f"  top signed ρ_spearman:   "
          f"{_fmt(info.get('top_signed_spearman', float('nan')))}")
    print(f"  pairs |ρ| ≥ {THRESHOLD_WEAK}:      "
          f"{info.get('n_above_weak', 0)}")
    print(f"  pairs |ρ| ≥ {THRESHOLD_DETECTED}:      "
          f"{info.get('n_above_detected', 0)}")
    print(f"  verdict:                 {verdict}")
    print("")
    print(f"PER-SEED CSV:  {PHASE4D_PER_SEED_CSV}")
    print(f"PER-CELL CSV:  {PHASE4D_PER_CELL_CSV}")
    print(f"MD:            {PHASE4D_MD}")


if __name__ == "__main__":
    main()
