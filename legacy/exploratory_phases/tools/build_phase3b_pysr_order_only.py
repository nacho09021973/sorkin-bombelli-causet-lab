#!/usr/bin/env python3
"""Phase 3B — PySR with strictly order-theoretic features.

Phase 3A showed that PySR recovers noise_epsilon and initial_energy as dominant
predictors — both are experimental design variables or embedding-dependent
quantities, not properties of the causal order.  Phase 3A therefore diagnosed
the pipeline correctly but did not find a CST law.

Phase 3B asks the harder question:

    Does the signal survive when we remove every feature that carries
    information about the embedding or the perturbation scale?

Feature set (X):
  Order-theoretic (computable from causal matrix alone, no embedding):
    mm_dim                   Myrheim-Meyer dimension estimator
    midpoint_dim             Midpoint-scaling dimension estimator
    abs_discrepancy_mm_mp    |mm_dim - midpoint_dim|
    chain2_count             # 2-chains (2-element chains)
    chain3_count             # 3-chains
    chain3_abundance         chain3 / chain4 fraction
    link_count               Hasse-diagram link count
    link_density             link_count / n  (scale-invariant)
    relation_count           total causal pairs
    ordering_fraction        relation_count / C(n,2)
    height                   longest chain length

  Coarse design variables (known a priori, not causet-internal):
    n                        system size
    target_dim               spacetime dimension
    noise_level              binary: 0 = small noise (ε=0.001),
                                     1 = medium noise (ε=0.05)
    warmup_mode_code         0 = skip_warmup
                             1 = legacy_warmup
                             2 = guarded_warmup

Excluded vs Phase 3A:
    noise_epsilon            continuous perturbation scale (design variable)
    initial_energy           embedding-dependent
    warmup_delta_energy      dynamic information
    warmup_acceptance_rate   dynamic information
    initial_interval_rmse    embedding-dependent
    initial_distance_to_*    embedding-dependent

Target (y):
    preserved_near_truth  (binary 0/1)

Rows:
    truth_plus_small_noise + truth_plus_medium_noise,
    all three warmup modes,
    joined across phase2e + phase2f.

The critical test: does abs_discrepancy_mm_mp (or any other pure
order-theoretic invariant) appear in the Pareto-front equations when
competing only with n, target_dim, noise_level, and warmup_mode_code?
If yes, that is a genuine causal set signal.

Output:
  benchmarks/foundation/phase3b_pysr_order_only.csv
  benchmarks/foundation/phase3b_pysr_order_only.md
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
FOUNDATION = ROOT / "benchmarks" / "foundation"

# ---------------------------------------------------------------------------
# Loaders
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


def load_invariants_json(path: Path) -> dict[tuple, dict]:
    """Return dict keyed by (d_spacetime, n, seed) → fingerprint dict."""
    data = json.loads(path.read_text(encoding="utf-8"))
    index: dict[tuple, dict] = {}
    for cell in data["cells"]:
        key = (cell["d_spacetime"], cell["n"], cell["seed"])
        fp = cell["fingerprint"]
        n = cell["n"]
        index[key] = {
            "link_count":        fp.get("link_count", 0),
            "link_density":      fp.get("link_count", 0) / n if n > 0 else 0.0,
            "relation_count":    fp.get("relation_count", 0),
            "ordering_fraction": fp.get("ordering_fraction", 0.0),
            "height":            fp.get("height", 0),
        }
    return index


def load_phase1d(path: Path) -> dict[tuple, dict]:
    """Return dict keyed by (family, target_dim, n, seed) → invariant cols."""
    index: dict[tuple, dict] = {}
    for row in load_csv(path):
        if row["family"] != "minkowski":
            continue
        key = (row["family"], row["target_dim"], row["n"], row["seed"])
        index[key] = {
            "mm_dim":                  row["mm_dim"],
            "midpoint_dim":            row["midpoint_dim"],
            "abs_discrepancy_mm_mp":   row["abs_discrepancy_mm_midpoint"],
            "chain2_count":            row["chain2_count"],
            "chain3_count":            row["chain3_count"],
            "chain3_abundance":        row["chain3_abundance"],
        }
    return index


# ---------------------------------------------------------------------------
# Encoding helpers
# ---------------------------------------------------------------------------

WARMUP_CODE = {"skip_warmup": 0, "legacy_warmup": 1, "guarded_warmup": 2}
NOISE_CODE  = {"truth_plus_small_noise": 0, "truth_plus_medium_noise": 1}

NEAR_TRUTH_LABELS = set(NOISE_CODE)

FEATURE_NAMES = [
    "n",
    "target_dim",
    "noise_level",
    "warmup_mode_code",
    "mm_dim",
    "midpoint_dim",
    "abs_discrepancy_mm_mp",
    "chain2_count",
    "chain3_count",
    "chain3_abundance",
    "link_count",
    "link_density",
    "relation_count",
    "ordering_fraction",
    "height",
]


def _safe(v, default: float = 0.0) -> float:
    if v is None:
        return default
    try:
        f = float(v)
        return default if (math.isnan(f) or math.isinf(f)) else f
    except (TypeError, ValueError):
        return default


# ---------------------------------------------------------------------------
# Dataset builder
# ---------------------------------------------------------------------------

def build_dataset(
    phase2_rows: list[dict],
    phase1d_idx: dict[tuple, dict],
    inv_idx:     dict[tuple, dict],
) -> tuple[list[list[float]], list[float]]:
    X, y = [], []

    for r in phase2_rows:
        if r["init_label"] not in NEAR_TRUTH_LABELS:
            continue
        if r["warmup_mode"] not in WARMUP_CODE:
            continue

        p1key = ("minkowski", r["target_dim"], r["n"], r["seed"])
        if p1key not in phase1d_idx:
            continue

        inv_key = (r["target_dim"], r["n"], r["seed"])
        if inv_key not in inv_idx:
            continue

        p1 = phase1d_idx[p1key]
        inv = inv_idx[inv_key]

        row_x = [
            _safe(r["n"]),
            _safe(r["target_dim"]),
            float(NOISE_CODE[r["init_label"]]),
            float(WARMUP_CODE[r["warmup_mode"]]),
            _safe(p1["mm_dim"]),
            _safe(p1["midpoint_dim"]),
            _safe(p1["abs_discrepancy_mm_mp"]),
            _safe(p1["chain2_count"]),
            _safe(p1["chain3_count"]),
            _safe(p1["chain3_abundance"]),
            _safe(inv["link_count"]),
            _safe(inv["link_density"]),
            _safe(inv["relation_count"]),
            _safe(inv["ordering_fraction"]),
            _safe(inv["height"]),
        ]

        if any(math.isnan(v) or math.isinf(v) for v in row_x):
            continue

        pnt = r["preserved_near_truth"]
        y.append(1.0 if pnt is True else 0.0)
        X.append(row_x)

    return X, y


# ---------------------------------------------------------------------------
# PySR runner
# ---------------------------------------------------------------------------

def run_pysr(
    X: list[list[float]],
    y: list[float],
    niterations: int,
) -> tuple[Any, Any] | tuple[None, None]:
    try:
        import numpy as np
        import pandas as pd
        from pysr import PySRRegressor
    except ImportError as e:
        print(f"  [skip] PySR unavailable ({e}). pip install pysr", file=sys.stderr)
        return None, None

    import numpy as np
    import pandas as pd

    X_df = pd.DataFrame(X, columns=FEATURE_NAMES)
    y_np = np.array(y, dtype=np.float64)

    pos = int(y_np.sum())
    print(f"\n  Samples: {len(y)}, preserved={pos}, destroyed={len(y)-pos}")
    print(f"  Features: {FEATURE_NAMES}")

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
# Output
# ---------------------------------------------------------------------------

CSV_HEADERS = ("complexity", "loss", "equation", "is_best", "uses_order_features")

ORDER_FEATURES = {
    "mm_dim", "midpoint_dim", "abs_discrepancy_mm_mp",
    "chain2_count", "chain3_count", "chain3_abundance",
    "link_count", "link_density", "relation_count",
    "ordering_fraction", "height",
}


def _uses_order(eq: str) -> bool:
    return any(f in eq for f in ORDER_FEATURES)


def write_csv(equations_df: Any, path: Path) -> None:
    if equations_df is None:
        path.write_text(",".join(CSV_HEADERS) + "\n", encoding="utf-8")
        return
    try:
        best_idx = equations_df["loss"].idxmin()
    except Exception:
        best_idx = None
    lines = [",".join(CSV_HEADERS)]
    for idx, row in equations_df.iterrows():
        eq = str(row.get("equation", ""))
        lines.append(",".join([
            str(row.get("complexity", "")),
            f"{row.get('loss', ''):.6g}",
            eq,
            "true" if idx == best_idx else "false",
            "true" if _uses_order(eq) else "false",
        ]))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_markdown(
    n_samples: int,
    n_pos: int,
    equations_df: Any,
    niterations: int,
    path: Path,
) -> None:
    lines = [
        "# Phase 3B — PySR: Strictly Order-Theoretic Features",
        "",
        "Critical test following Phase 3A: does the signal survive when",
        "`noise_epsilon`, `initial_energy`, and all warmup-dynamic features",
        "are removed?  Only causal-matrix observables and coarse design",
        "variables are allowed as predictors.",
        "",
        "## Setup",
        "",
        f"- Samples: {n_samples} ({n_pos} preserved, {n_samples - n_pos} destroyed)",
        f"- PySR iterations: {niterations}",
        "",
        "**Input features (X):**",
        "",
        "| feature | type |",
        "| --- | --- |",
        "| `n` | design |",
        "| `target_dim` | design |",
        "| `noise_level` | design (0=small, 1=medium) |",
        "| `warmup_mode_code` | design (0=skip, 1=legacy, 2=guarded) |",
        "| `mm_dim` | order-theoretic |",
        "| `midpoint_dim` | order-theoretic |",
        "| `abs_discrepancy_mm_mp` | order-theoretic |",
        "| `chain2_count`, `chain3_count`, `chain3_abundance` | order-theoretic |",
        "| `link_count`, `link_density` | order-theoretic |",
        "| `relation_count`, `ordering_fraction` | order-theoretic |",
        "| `height` | order-theoretic |",
        "",
        "**Excluded vs Phase 3A:** `noise_epsilon`, `initial_energy`,",
        "`warmup_delta_energy`, `warmup_acceptance_rate`, and all",
        "embedding-dependent features.",
        "",
        "**Target:** `preserved_near_truth` ∈ {0, 1}",
        "",
        "**Key diagnostic:** equations flagged `uses_order_features=true` contain",
        "at least one purely combinatorial causet observable.  If such equations",
        "achieve lower loss than design-variable-only equations of comparable",
        "complexity, the order-theoretic signal is real.",
        "",
        "## Discovered equations (Pareto front)",
        "",
        "| complexity | loss | equation | best | order features? |",
        "| ---: | ---: | --- | :---: | :---: |",
    ]

    if equations_df is None:
        lines.append("| — | — | *PySR not available* | | |")
    else:
        try:
            best_idx = equations_df["loss"].idxmin()
        except Exception:
            best_idx = None
        for idx, row in equations_df.iterrows():
            eq = str(row.get("equation", ""))
            marker = "**★**" if idx == best_idx else ""
            order_flag = "✓" if _uses_order(eq) else "—"
            try:
                loss_str = f"{float(row.get('loss', '')):.4g}"
            except (TypeError, ValueError):
                loss_str = str(row.get("loss", ""))
            lines.append(
                f"| {row.get('complexity', '')} | {loss_str} | `{eq}` | {marker} | {order_flag} |"
            )

    lines += [
        "",
        "## Interpretation guide",
        "",
        "- If the best equation uses only `noise_level` and/or `warmup_mode_code`:",
        "  the order-theoretic signal does **not** survive — the dataset size or",
        "  the noise-level dominance is too strong.",
        "- If `abs_discrepancy_mm_mp`, `ordering_fraction`, or any other",
        "  order-theoretic feature appears in competitive equations: the signal",
        "  **survives** leakage removal and warrants deeper study.",
        "",
        "## Data provenance",
        "",
        "- `phase1d_structural_atlas.csv` (mm_dim, midpoint_dim, chain stats)",
        "- `invariants.json` (link_count, relation_count, ordering_fraction, height)",
        "- `phase2e_warmup_skip_probe.csv` + `phase2f_guarded_warmup_probe.csv`",
        "- Join key: `(target_dim, n, seed)` — minkowski only.",
        "",
        "Regenerate via `make regen-phase3b`.",
        "Source: `tools/build_phase3b_pysr_order_only.py`.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--niterations", type=int, default=100)
    args = parser.parse_args()

    phase1d_idx = load_phase1d(FOUNDATION / "phase1d_structural_atlas.csv")
    inv_idx     = load_invariants_json(FOUNDATION / "invariants.json")

    phase2e = load_csv(FOUNDATION / "phase2e_warmup_skip_probe.csv")
    phase2f = load_csv(FOUNDATION / "phase2f_guarded_warmup_probe.csv")
    for r in phase2e:
        if r.get("warmup_mode") == "with_warmup":
            r["warmup_mode"] = "legacy_warmup"

    all_rows = phase2e + phase2f
    X, y = build_dataset(all_rows, phase1d_idx, inv_idx)

    print(f"Phase 3B: {len(X)} samples after joining and filtering.")
    if not X:
        print("ERROR: empty dataset. Check join keys.", file=sys.stderr)
        sys.exit(1)

    n_pos = int(sum(y))
    print(f"  preserved={n_pos}, destroyed={len(y)-n_pos}")

    _, equations_df = run_pysr(X, y, args.niterations)

    FOUNDATION.mkdir(parents=True, exist_ok=True)
    write_csv(equations_df, FOUNDATION / "phase3b_pysr_order_only.csv")
    write_markdown(len(X), n_pos, equations_df, args.niterations,
                   FOUNDATION / "phase3b_pysr_order_only.md")

    if equations_df is not None:
        try:
            best = equations_df.loc[equations_df["loss"].idxmin(), "equation"]
            uses = _uses_order(best)
            verdict = "ORDER SIGNAL PRESENT" if uses else "NO ORDER SIGNAL — design variables dominate"
            print(f"\nBest equation: {best}")
            print(f"Verdict: {verdict}")
        except Exception:
            pass

    print(f"\nOutput: {FOUNDATION / 'phase3b_pysr_order_only.csv'}")
    print(f"        {FOUNDATION / 'phase3b_pysr_order_only.md'}")


if __name__ == "__main__":
    main()
