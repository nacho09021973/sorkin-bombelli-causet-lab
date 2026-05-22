#!/usr/bin/env python3
"""Phase 3C — Ablation: four feature panels, maxsize=12.

Phase 3B showed order-theoretic features entering the Pareto front, but
noise_level (a design variable) was still present in every competitive
equation.  Phase 3C closes the leakage question with a controlled ablation
across four panels:

  A  noise-only baseline
       X = [noise_level]
       Establishes the ceiling a design variable alone gives.

  B  order + noise
       X = [noise_level, warmup_mode_code] + all order features
       Measures the incremental lift of order features over the noise baseline.

  C  order-only known-d
       X = [n, target_dim] + all order features
       noise_level excluded.  target_dim accepted as generation metadata.

  D  order-only no-oracle
       X = [n] + all order features
       noise_level AND target_dim excluded.
       This is the only fully autonomous version: everything computable
       from the causal matrix alone, plus the cardinality n.

Order features used in B/C/D:
  mm_dim, midpoint_dim, abs_discrepancy_mm_mp,
  chain2_count, chain3_count, chain3_abundance,
  link_count, link_density,
  relation_count, ordering_fraction, height

Complexity cap: maxsize=12.
  With 180 samples and 3 seeds, complexity-20 equations risk opportunistic
  overfitting.  The cap forces parsimony.

Decision criterion:
  If D's best loss beats the constant baseline (~0.249) by more than A's
  margin AND order features appear in D's equations → genuine order signal.
  If D cannot beat the constant baseline, order features alone are
  insufficient at this dataset scale.

Target: preserved_near_truth ∈ {0, 1}
Rows:   truth_plus_small_noise + truth_plus_medium_noise,
        all warmup modes, phase2e + phase2f.

Output:
  benchmarks/foundation/phase3c_ablation.csv
  benchmarks/foundation/phase3c_ablation.md
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
# Loaders  (identical pattern to Phase 3B)
# ---------------------------------------------------------------------------

def _parse_val(s: str) -> Any:
    s = s.strip()
    if s in ("NA", ""):
        return None
    if s == "true":  return True
    if s == "false": return False
    try:    return int(s)
    except ValueError: pass
    try:    return float(s)
    except ValueError: return s


def load_csv(path: Path) -> list[dict]:
    rows = []
    with path.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            rows.append({k: _parse_val(v) for k, v in row.items()})
    return rows


def load_invariants_json(path: Path) -> dict[tuple, dict]:
    data = json.loads(path.read_text(encoding="utf-8"))
    index: dict[tuple, dict] = {}
    for cell in data["cells"]:
        key = (cell["d_spacetime"], cell["n"], cell["seed"])
        fp  = cell["fingerprint"]
        n   = cell["n"]
        index[key] = {
            "link_count":        fp.get("link_count", 0),
            "link_density":      fp.get("link_count", 0) / n if n > 0 else 0.0,
            "relation_count":    fp.get("relation_count", 0),
            "ordering_fraction": fp.get("ordering_fraction", 0.0),
            "height":            fp.get("height", 0),
        }
    return index


def load_phase1d(path: Path) -> dict[tuple, dict]:
    index: dict[tuple, dict] = {}
    for row in load_csv(path):
        if row["family"] != "minkowski":
            continue
        key = (row["family"], row["target_dim"], row["n"], row["seed"])
        index[key] = {
            "mm_dim":               row["mm_dim"],
            "midpoint_dim":         row["midpoint_dim"],
            "abs_discrepancy_mm_mp": row["abs_discrepancy_mm_midpoint"],
            "chain2_count":         row["chain2_count"],
            "chain3_count":         row["chain3_count"],
            "chain3_abundance":     row["chain3_abundance"],
        }
    return index


# ---------------------------------------------------------------------------
# Panel definitions
# ---------------------------------------------------------------------------

NEAR_TRUTH_LABELS = {"truth_plus_small_noise", "truth_plus_medium_noise"}
WARMUP_CODE = {"skip_warmup": 0, "legacy_warmup": 1, "guarded_warmup": 2}
NOISE_CODE  = {"truth_plus_small_noise": 0, "truth_plus_medium_noise": 1}

ORDER_FEATURE_NAMES = [
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

ORDER_FEATURE_SET = set(ORDER_FEATURE_NAMES)

PANELS: dict[str, list[str]] = {
    "A": ["noise_level"],
    "B": ["noise_level", "warmup_mode_code"] + ORDER_FEATURE_NAMES,
    "C": ["n", "target_dim"] + ORDER_FEATURE_NAMES,
    "D": ["n"] + ORDER_FEATURE_NAMES,
}

PANEL_DESCRIPTIONS: dict[str, str] = {
    "A": "noise-only baseline — no order features, no target_dim",
    "B": "order + noise — all order features plus noise_level and warmup_mode",
    "C": "order-only known-d — order features + target_dim, no noise_level",
    "D": "order-only no-oracle — fully autonomous (n + order features only)",
}


# ---------------------------------------------------------------------------
# Dataset assembly
# ---------------------------------------------------------------------------

def _safe(v, default: float = 0.0) -> float:
    if v is None: return default
    try:
        f = float(v)
        return default if (math.isnan(f) or math.isinf(f)) else f
    except (TypeError, ValueError):
        return default


def build_base_rows(
    phase2_rows: list[dict],
    phase1d_idx: dict[tuple, dict],
    inv_idx:     dict[tuple, dict],
) -> list[dict]:
    """Build one enriched dict per usable row (all features + target)."""
    out = []
    for r in phase2_rows:
        if r["init_label"] not in NEAR_TRUTH_LABELS:
            continue
        if r["warmup_mode"] not in WARMUP_CODE:
            continue

        p1key  = ("minkowski", r["target_dim"], r["n"], r["seed"])
        invkey = (r["target_dim"], r["n"], r["seed"])
        if p1key not in phase1d_idx or invkey not in inv_idx:
            continue

        p1  = phase1d_idx[p1key]
        inv = inv_idx[invkey]

        row = {
            "noise_level":      float(NOISE_CODE[r["init_label"]]),
            "warmup_mode_code": float(WARMUP_CODE[r["warmup_mode"]]),
            "n":                _safe(r["n"]),
            "target_dim":       _safe(r["target_dim"]),
            "mm_dim":                p1["mm_dim"],
            "midpoint_dim":          p1["midpoint_dim"],
            "abs_discrepancy_mm_mp": p1["abs_discrepancy_mm_mp"],
            "chain2_count":          p1["chain2_count"],
            "chain3_count":          p1["chain3_count"],
            "chain3_abundance":      p1["chain3_abundance"],
            "link_count":            inv["link_count"],
            "link_density":          inv["link_density"],
            "relation_count":        inv["relation_count"],
            "ordering_fraction":     inv["ordering_fraction"],
            "height":                inv["height"],
            "_preserved":       1.0 if r["preserved_near_truth"] is True else 0.0,
        }
        if any(math.isnan(float(v)) or math.isinf(float(v))
               for v in row.values() if isinstance(v, (int, float))):
            continue
        out.append(row)
    return out


def panel_arrays(
    base_rows: list[dict],
    feature_names: list[str],
) -> tuple[list[list[float]], list[float]]:
    X = [[_safe(r[f]) for f in feature_names] for r in base_rows]
    y = [r["_preserved"] for r in base_rows]
    return X, y


# ---------------------------------------------------------------------------
# PySR runner
# ---------------------------------------------------------------------------

def run_pysr(
    X: list[list[float]],
    y: list[float],
    feature_names: list[str],
    niterations: int,
    maxsize: int,
    label: str,
) -> Any | None:
    try:
        import numpy as np
        import pandas as pd
        from pysr import PySRRegressor
    except ImportError as e:
        print(f"  [skip] PySR unavailable: {e}", file=sys.stderr)
        return None

    import numpy as np
    import pandas as pd

    X_df = pd.DataFrame(X, columns=feature_names)
    y_np = np.array(y, dtype=np.float64)

    print(f"\n  Panel {label} ({len(feature_names)} features): "
          f"{len(y)} samples, preserved={int(y_np.sum())}/{len(y)}")

    model = PySRRegressor(
        niterations=niterations,
        binary_operators=["+", "-", "*", "/"],
        unary_operators=["log", "sqrt", "square", "abs"],
        maxsize=maxsize,
        populations=20,
        model_selection="best",
        verbosity=0,
        random_state=1959,
        deterministic=True,
        parallelism="serial",
    )
    model.fit(X_df, y_np)
    return model.equations_


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def _uses_order(eq: str) -> bool:
    return any(f in eq for f in ORDER_FEATURE_SET)


def _uses_design(eq: str) -> bool:
    return any(f in eq for f in ("noise_level", "warmup_mode_code", "target_dim"))


CSV_HEADERS = (
    "panel", "complexity", "loss", "equation",
    "is_best", "uses_order", "uses_design",
)


def equations_to_rows(eqs: Any, panel: str) -> list[dict]:
    if eqs is None:
        return []
    try:
        best_idx = eqs["loss"].idxmin()
    except Exception:
        best_idx = None
    rows = []
    for idx, row in eqs.iterrows():
        eq = str(row.get("equation", ""))
        rows.append({
            "panel":       panel,
            "complexity":  row.get("complexity", ""),
            "loss":        f"{row.get('loss', ''):.6g}",
            "equation":    eq,
            "is_best":     "true" if idx == best_idx else "false",
            "uses_order":  "true" if _uses_order(eq) else "false",
            "uses_design": "true" if _uses_design(eq) else "false",
        })
    return rows


def write_csv(all_rows: list[dict], path: Path) -> None:
    lines = [",".join(CSV_HEADERS)]
    for r in all_rows:
        lines.append(",".join(str(r.get(h, "")) for h in CSV_HEADERS))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# Markdown report
# ---------------------------------------------------------------------------

def _panel_table(eqs: Any, panel: str) -> list[str]:
    lines = [
        f"| complexity | loss | equation | best | order? | design? |",
        "| ---: | ---: | --- | :---: | :---: | :---: |",
    ]
    if eqs is None:
        lines.append("| — | — | *PySR unavailable* | | | |")
        return lines
    try:
        best_idx = eqs["loss"].idxmin()
    except Exception:
        best_idx = None
    for idx, row in eqs.iterrows():
        eq = str(row.get("equation", ""))
        marker = "**★**" if idx == best_idx else ""
        try:
            loss_str = f"{float(row.get('loss', '')):.4g}"
        except (TypeError, ValueError):
            loss_str = str(row.get("loss", ""))
        lines.append(
            f"| {row.get('complexity','')} | {loss_str} | `{eq}` | {marker}"
            f" | {'✓' if _uses_order(eq) else '—'}"
            f" | {'✓' if _uses_design(eq) else '—'} |"
        )
    return lines


def _best_loss(eqs: Any) -> float | None:
    if eqs is None:
        return None
    try:
        return float(eqs["loss"].min())
    except Exception:
        return None


def write_markdown(
    n_samples: int,
    n_pos: int,
    panel_eqs: dict[str, Any],
    niterations: int,
    maxsize: int,
    path: Path,
) -> None:
    constant_loss = n_pos / n_samples * (1 - n_pos / n_samples) * 2
    lines = [
        "# Phase 3C — Ablation: Four Feature Panels",
        "",
        "Controlled ablation to determine whether order-theoretic features",
        "carry signal independent of experimental design variables.",
        "",
        "## Setup",
        "",
        f"- Samples: {n_samples} ({n_pos} preserved, {n_samples-n_pos} destroyed)",
        f"- PySR iterations per panel: {niterations}",
        f"- maxsize: {maxsize}  ← complexity cap to suppress opportunistic fitting",
        "",
        "## Panel definitions",
        "",
        "| panel | features | noise_level | target_dim | order features |",
        "| --- | --- | :---: | :---: | :---: |",
        "| A | noise-only baseline | ✓ | — | — |",
        "| B | order + noise | ✓ | ✓ | ✓ |",
        "| C | order-only known-d | — | ✓ | ✓ |",
        "| D | order-only no-oracle | — | — | ✓ |",
        "",
        "## Summary: best loss per panel",
        "",
        "| panel | description | best loss | Δ vs constant |",
        "| --- | --- | ---: | ---: |",
    ]

    for p in "ABCD":
        bl = _best_loss(panel_eqs.get(p))
        if bl is None:
            lines.append(f"| {p} | {PANEL_DESCRIPTIONS[p]} | — | — |")
        else:
            delta = constant_loss - bl
            lines.append(
                f"| {p} | {PANEL_DESCRIPTIONS[p]} | {bl:.4g} | {delta:+.4g} |"
            )

    lines += [
        "",
        f"Constant-predictor baseline loss (majority class): ~{constant_loss:.4g}",
        "",
        "**Decision rule:**",
        "- If D's Δ vs constant is meaningful and order features appear in D's best equation:",
        "  → genuine order-theoretic signal, warrants more data.",
        "- If only A/B achieve large Δ: signal depends on experimental design variables.",
        "- If all panels are near baseline: dataset too small or noise dominates completely.",
        "",
    ]

    for p in "ABCD":
        eqs = panel_eqs.get(p)
        lines += [
            f"## Panel {p} — {PANEL_DESCRIPTIONS[p]}",
            "",
            f"Features: `{', '.join(PANELS[p])}`",
            "",
            *_panel_table(eqs, p),
            "",
        ]

    lines += [
        "## Data provenance",
        "",
        "- `phase1d_structural_atlas.csv`, `invariants.json`",
        "- `phase2e_warmup_skip_probe.csv`, `phase2f_guarded_warmup_probe.csv`",
        "- Join key: `(target_dim, n, seed)` — minkowski only.",
        "",
        "Regenerate via `make regen-phase3c`.",
        "Source: `tools/build_phase3c_ablation.py`.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--niterations", type=int, default=100)
    parser.add_argument("--maxsize",     type=int, default=12)
    args = parser.parse_args()

    phase1d_idx = load_phase1d(FOUNDATION / "phase1d_structural_atlas.csv")
    inv_idx     = load_invariants_json(FOUNDATION / "invariants.json")

    phase2e = load_csv(FOUNDATION / "phase2e_warmup_skip_probe.csv")
    phase2f = load_csv(FOUNDATION / "phase2f_guarded_warmup_probe.csv")
    for r in phase2e:
        if r.get("warmup_mode") == "with_warmup":
            r["warmup_mode"] = "legacy_warmup"

    base_rows = build_base_rows(phase2e + phase2f, phase1d_idx, inv_idx)
    print(f"Phase 3C: {len(base_rows)} base rows.")

    n_pos = int(sum(r["_preserved"] for r in base_rows))

    panel_eqs: dict[str, Any] = {}
    all_csv_rows: list[dict] = []

    for panel_id, feature_names in PANELS.items():
        X, y = panel_arrays(base_rows, feature_names)
        eqs  = run_pysr(X, y, feature_names, args.niterations, args.maxsize, panel_id)
        panel_eqs[panel_id] = eqs
        all_csv_rows.extend(equations_to_rows(eqs, panel_id))

    FOUNDATION.mkdir(parents=True, exist_ok=True)
    write_csv(all_csv_rows, FOUNDATION / "phase3c_ablation.csv")
    write_markdown(len(base_rows), n_pos, panel_eqs, args.niterations, args.maxsize,
                   FOUNDATION / "phase3c_ablation.md")

    print("\n--- Phase 3C summary ---")
    for p in "ABCD":
        bl = _best_loss(panel_eqs.get(p))
        if bl is not None:
            eqs = panel_eqs[p]
            try:
                best_eq = eqs.loc[eqs["loss"].idxmin(), "equation"]
            except Exception:
                best_eq = "—"
            order_flag = "ORDER" if _uses_order(best_eq) else "design-only"
            print(f"  Panel {p}: best loss={bl:.4g}  [{order_flag}]  {best_eq[:80]}")

    print(f"\nOutput: {FOUNDATION / 'phase3c_ablation.csv'}")
    print(f"        {FOUNDATION / 'phase3c_ablation.md'}")


if __name__ == "__main__":
    main()
