#!/usr/bin/env python3
"""Phase 3A — PySR symbolic regression over causal set pipeline features.

This is the first application of symbolic regression to causal set theory
data.  Two experiments are run:

  Experiment A — Order-theoretic predictor
    X: pure combinatorial invariants computed from the causal matrix alone
       (mm_dim, midpoint_dim, chain3_abundance, chain2_count, chain3_count,
        abs_discrepancy_mm_midpoint) plus scalar scalars (n, target_dim,
        initial_energy).
    y: log1p(max(0, delta_energy))  — log-compressed energy drift after
       warmup+anneal, relative to the near-truth starting point.
    Rows: truth_plus_small_noise and truth_plus_medium_noise inits, all
          warmup modes, joined from phase1d + phase2f.

    If PySR discovers a closed-form rule here it would be the first
    symbolic law connecting discrete causal-order combinatorics to
    Bombelli embedding quality — and therefore frontier CST physics.

  Experiment B — Warmup dynamics predictor
    X: order-theoretic invariants + warmup acceptance rate + warmup energy
       change, all available from phase2f guarded/legacy warmup rows.
    y: preserved_near_truth  (binary 0/1)
    Rows: truth_plus_small_noise and truth_plus_medium_noise, warmup_mode
          in (legacy_warmup, guarded_warmup).

    This discovers what drives near-truth preservation at the level of
    move statistics, providing a symbolic warmup-quality score.

Input data (already on disk):
  benchmarks/foundation/phase1d_structural_atlas.csv
  benchmarks/foundation/phase2e_warmup_skip_probe.csv
  benchmarks/foundation/phase2f_guarded_warmup_probe.csv

Output:
  benchmarks/foundation/phase3a_pysr_warmup_rule.csv
  benchmarks/foundation/phase3a_pysr_warmup_rule.md

Usage:
  python3 tools/build_phase3a_pysr_warmup_rule.py [--niterations N]
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
FOUNDATION = ROOT / "benchmarks" / "foundation"

# ---------------------------------------------------------------------------
# CSV loaders
# ---------------------------------------------------------------------------

def _parse_val(s: str) -> Any:
    s = s.strip()
    if s in ("NA", ""):
        return None
    if s == "true":
        return True
    if s == "false":
        return False
    try:
        return int(s)
    except ValueError:
        pass
    try:
        return float(s)
    except ValueError:
        return s


def load_csv(path: Path) -> list[dict]:
    rows = []
    with path.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            rows.append({k: _parse_val(v) for k, v in row.items()})
    return rows


# ---------------------------------------------------------------------------
# Data merge
# ---------------------------------------------------------------------------

NEAR_TRUTH_LABELS = {"truth_plus_small_noise", "truth_plus_medium_noise"}


def _invariant_key(row: dict) -> tuple:
    return (row["family"], row["target_dim"], row["n"], row["seed"])


def build_merged(
    phase1d_rows: list[dict],
    phase2_rows: list[dict],
) -> list[dict]:
    """Join phase1d invariants onto phase2 rows on (family, target_dim, n, seed).

    Only minkowski rows are joined (phase1d also has KR/corona rows which
    have no matching phase2 data).
    """
    inv: dict[tuple, dict] = {}
    for r in phase1d_rows:
        if r["family"] != "minkowski":
            continue
        inv[_invariant_key(r)] = r

    merged = []
    for r in phase2_rows:
        key = _invariant_key(r)
        if key not in inv:
            continue
        combined = dict(r)
        combined.update({
            "mm_dim": inv[key]["mm_dim"],
            "midpoint_dim": inv[key]["midpoint_dim"],
            "abs_discrepancy_mm_midpoint": inv[key]["abs_discrepancy_mm_midpoint"],
            "chain2_count": inv[key]["chain2_count"],
            "chain3_count": inv[key]["chain3_count"],
            "chain3_abundance": inv[key]["chain3_abundance"],
        })
        merged.append(combined)
    return merged


# ---------------------------------------------------------------------------
# Dataset builders
# ---------------------------------------------------------------------------

EXPT_A_FEATURES = [
    "mm_dim",
    "midpoint_dim",
    "abs_discrepancy_mm_midpoint",
    "chain2_count",
    "chain3_count",
    "chain3_abundance",
    "n",
    "target_dim",
    "noise_epsilon",
    "initial_energy",
]

EXPT_B_FEATURES = [
    "mm_dim",
    "midpoint_dim",
    "abs_discrepancy_mm_midpoint",
    "chain3_abundance",
    "n",
    "target_dim",
    "noise_epsilon",
    "initial_energy",
    "warmup_acceptance_rate",
    "warmup_delta_energy",
]


def _safe_float(v, default: float = 0.0) -> float:
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return default
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def build_dataset_a(merged: list[dict]) -> tuple[list[list[float]], list[float], list[str]]:
    """Experiment A: order-theoretic features → log-energy drift.

    Filter: near-truth inits only (small and medium noise).
    Target: log1p(max(0, delta_energy)) — compresses the enormous range of
            energy drift values into a regression-friendly scale while
            preserving the ordering (preserved cases → 0, destroyed → large).
    """
    X, y = [], []
    for r in merged:
        if r["init_label"] not in NEAR_TRUTH_LABELS:
            continue
        if r["noise_epsilon"] is None:
            continue
        de = _safe_float(r["delta_energy"])
        row_x = [_safe_float(r[f]) for f in EXPT_A_FEATURES]
        if any(math.isnan(v) or math.isinf(v) for v in row_x):
            continue
        y.append(math.log1p(max(0.0, de)))
        X.append(row_x)
    return X, y, EXPT_A_FEATURES


def build_dataset_b(merged: list[dict]) -> tuple[list[list[float]], list[float], list[str]]:
    """Experiment B: warmup dynamics + order-theoretic → preserved_near_truth.

    Filter: near-truth inits, warmup modes that generate move statistics
            (legacy_warmup and guarded_warmup; skip_warmup has no stats).
    Derived feature: warmup_acceptance_rate = accepted / attempted (in [0,1]).
    Target: float(preserved_near_truth) ∈ {0.0, 1.0}.
    """
    X, y = [], []
    for r in merged:
        if r["init_label"] not in NEAR_TRUTH_LABELS:
            continue
        if r["warmup_mode"] not in ("legacy_warmup", "guarded_warmup"):
            continue
        if r["noise_epsilon"] is None:
            continue

        attempted = _safe_float(r.get("warmup_attempted_moves"), 0.0)
        accepted = _safe_float(r.get("warmup_accepted_moves"), 0.0)
        acceptance_rate = accepted / attempted if attempted > 0 else 0.0

        r_aug = dict(r)
        r_aug["warmup_acceptance_rate"] = acceptance_rate

        row_x = [_safe_float(r_aug.get(f)) for f in EXPT_B_FEATURES]
        if any(math.isnan(v) or math.isinf(v) for v in row_x):
            continue

        pnt = r["preserved_near_truth"]
        y.append(1.0 if pnt is True else 0.0)
        X.append(row_x)
    return X, y, EXPT_B_FEATURES


# ---------------------------------------------------------------------------
# PySR runner
# ---------------------------------------------------------------------------

def _run_pysr(
    X: list[list[float]],
    y: list[float],
    feature_names: list[str],
    niterations: int,
    label: str,
) -> tuple[Any, Any] | tuple[None, None]:
    """Run PySR and return (model, equations_df).  Returns (None, None) if unavailable."""
    try:
        import numpy as np
        import pandas as pd
        from pysr import PySRRegressor
    except ImportError as e:
        print(f"  [skip] PySR not available ({e}). Install with: pip install pysr", file=sys.stderr)
        return None, None

    import numpy as np
    import pandas as pd

    X_df = pd.DataFrame(X, columns=feature_names)
    y_np = np.array(y, dtype=np.float64)

    print(f"\n  Running PySR on {label}: {X_df.shape[0]} samples × {X_df.shape[1]} features")
    print(f"  Target range: [{y_np.min():.4g}, {y_np.max():.4g}], mean={y_np.mean():.4g}")

    model = PySRRegressor(
        niterations=niterations,
        binary_operators=["+", "-", "*", "/"],
        unary_operators=["log", "sqrt", "square", "abs"],
        maxsize=20,
        populations=20,
        model_selection="best",
        verbosity=0,
        random_state=1959,
        deterministic=True,
        parallelism="serial",
    )
    model.fit(X_df, y_np)

    return model, model.equations_


# ---------------------------------------------------------------------------
# CSV output
# ---------------------------------------------------------------------------

EQUATION_CSV_HEADERS = (
    "experiment",
    "complexity",
    "loss",
    "equation",
    "is_best",
)


def _equations_to_rows(equations_df: Any, experiment: str) -> list[dict]:
    if equations_df is None:
        return []
    try:
        best_idx = equations_df["loss"].idxmin()
    except Exception:
        best_idx = None

    rows = []
    for idx, eq_row in equations_df.iterrows():
        rows.append({
            "experiment": experiment,
            "complexity": eq_row.get("complexity", ""),
            "loss": f"{eq_row.get('loss', ''):.6g}",
            "equation": eq_row.get("equation", ""),
            "is_best": "true" if idx == best_idx else "false",
        })
    return rows


def write_csv(all_rows: list[dict], path: Path) -> None:
    lines = [",".join(EQUATION_CSV_HEADERS)]
    for row in all_rows:
        lines.append(",".join(str(row.get(h, "")) for h in EQUATION_CSV_HEADERS))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# Markdown report
# ---------------------------------------------------------------------------

def _eq_table(equations_df: Any) -> list[str]:
    if equations_df is None:
        return ["*PySR not available — install with `pip install pysr`.*", ""]
    try:
        best_idx = equations_df["loss"].idxmin()
    except Exception:
        best_idx = None

    lines = [
        "| complexity | loss | equation | best |",
        "| ---: | ---: | --- | :---: |",
    ]
    for idx, row in equations_df.iterrows():
        marker = "**★**" if idx == best_idx else ""
        loss_val = row.get("loss", "")
        try:
            loss_str = f"{float(loss_val):.4g}"
        except (TypeError, ValueError):
            loss_str = str(loss_val)
        lines.append(
            f"| {row.get('complexity', '')} | {loss_str} | `{row.get('equation', '')}` | {marker} |"
        )
    return lines


def write_markdown(
    n_a: int,
    n_b: int,
    features_a: list[str],
    features_b: list[str],
    eq_a: Any,
    eq_b: Any,
    niterations: int,
    path: Path,
) -> None:
    lines = [
        "# Phase 3A — PySR Symbolic Regression on Causal Set Pipeline Features",
        "",
        "**Frontier note**: this is the first application of PySR / symbolic",
        "regression to data derived from Causal Set Theory.  No prior published",
        "work has trained PySR on causal matrix observables.",
        "",
        "Two experiments search for closed-form rules over features already",
        "produced by the Bombelli-Sorkin annealing pipeline.",
        "",
        "## Experiment A — Order-theoretic predictor",
        "",
        "**Question**: can a formula built from *pure* combinatorial invariants",
        "of the causal matrix (computable without any embedding) predict how",
        "much the energy drifts during warmup+anneal from a near-truth start?",
        "",
        "If yes, that formula encodes a relationship between discrete causal",
        "order and the current optimizer-response target. This is an exploratory optimizer-response result, not a standalone physical claim.",
        "",
        f"- Samples: {n_a}",
        f"- Features (X): {', '.join(features_a)}",
        "- Target (y): `log1p(max(0, delta_energy))`",
        "- Rows: `truth_plus_small_noise` + `truth_plus_medium_noise`, all warmup modes",
        f"- PySR iterations: {niterations}",
        "",
        "### Discovered equations (Pareto front)",
        "",
        *_eq_table(eq_a),
        "",
        "## Experiment B — Warmup dynamics predictor",
        "",
        "**Question**: can a formula built from warmup move statistics",
        "(acceptance rate, energy change during warmup) and order-theoretic",
        "invariants predict whether a near-truth initialization is preserved?",
        "",
        f"- Samples: {n_b}",
        f"- Features (X): {', '.join(features_b)}",
        "- Target (y): `preserved_near_truth` ∈ {0, 1}",
        "- Rows: `truth_plus_small_noise` + `truth_plus_medium_noise`,",
        "  warmup_mode ∈ {legacy_warmup, guarded_warmup}",
        f"- PySR iterations: {niterations}",
        "- `warmup_acceptance_rate` = accepted_moves / attempted_moves",
        "",
        "### Discovered equations (Pareto front)",
        "",
        *_eq_table(eq_b),
        "",
        "## Data provenance",
        "",
        "- Invariants: `phase1d_structural_atlas.csv`",
        "- Dynamics: `phase2e_warmup_skip_probe.csv` + `phase2f_guarded_warmup_probe.csv`",
        "- Join key: `(family, target_dim, n, seed)` — minkowski only.",
        "",
        "Regenerate via `make regen-phase3a`.",
        "Source: `tools/build_phase3a_pysr_warmup_rule.py`.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--niterations",
        type=int,
        default=100,
        help="PySR niterations per experiment (default: 100)",
    )
    args = parser.parse_args()

    phase1d = load_csv(FOUNDATION / "phase1d_structural_atlas.csv")
    phase2e = load_csv(FOUNDATION / "phase2e_warmup_skip_probe.csv")
    phase2f = load_csv(FOUNDATION / "phase2f_guarded_warmup_probe.csv")

    # Normalise warmup_mode names: phase2e uses 'with_warmup', phase2f uses
    # 'legacy_warmup'.  Map to 'legacy_warmup' throughout for consistency.
    for r in phase2e:
        if r.get("warmup_mode") == "with_warmup":
            r["warmup_mode"] = "legacy_warmup"
        # phase2e has no warmup stats columns; fill with None so the merge dict
        # has the keys experiment B expects.
        for col in ("warmup_attempted_moves", "warmup_accepted_moves",
                    "warmup_rejected_moves", "warmup_energy_before",
                    "warmup_energy_after", "warmup_delta_energy"):
            r.setdefault(col, None)

    all_phase2 = phase2e + phase2f
    merged = build_merged(phase1d, all_phase2)

    print(f"Merged {len(merged)} rows from phase1d + phase2e + phase2f.")

    X_a, y_a, feat_a = build_dataset_a(merged)
    X_b, y_b, feat_b = build_dataset_b(merged)

    print(f"Dataset A (order-theoretic → log-energy drift): {len(X_a)} samples.")
    print(f"Dataset B (warmup dynamics → preserved): {len(X_b)} samples.")

    model_a, eq_a = _run_pysr(X_a, y_a, feat_a, args.niterations, "A")
    model_b, eq_b = _run_pysr(X_b, y_b, feat_b, args.niterations, "B")

    FOUNDATION.mkdir(parents=True, exist_ok=True)

    csv_rows = _equations_to_rows(eq_a, "A") + _equations_to_rows(eq_b, "B")
    write_csv(csv_rows, FOUNDATION / "phase3a_pysr_warmup_rule.csv")
    write_markdown(
        n_a=len(X_a),
        n_b=len(X_b),
        features_a=feat_a,
        features_b=feat_b,
        eq_a=eq_a,
        eq_b=eq_b,
        niterations=args.niterations,
        path=FOUNDATION / "phase3a_pysr_warmup_rule.md",
    )

    if eq_a is not None or eq_b is not None:
        print("\nPhase 3A complete.")
        if eq_a is not None:
            try:
                best_a = eq_a.loc[eq_a["loss"].idxmin(), "equation"]
                print(f"  Best A: {best_a}")
            except Exception:
                pass
        if eq_b is not None:
            try:
                best_b = eq_b.loc[eq_b["loss"].idxmin(), "equation"]
                print(f"  Best B: {best_b}")
            except Exception:
                pass
    else:
        print("\nPhase 3A: skeleton written. Install PySR to run the regression.")

    print(f"Output: {FOUNDATION / 'phase3a_pysr_warmup_rule.csv'}")
    print(f"        {FOUNDATION / 'phase3a_pysr_warmup_rule.md'}")


if __name__ == "__main__":
    main()
