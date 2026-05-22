#!/usr/bin/env python3
"""Phase 3E — Methodological cleanup of Phase 3D.

Phase 3D residualized the warmup-relative drift by
  groupby = (noise_level, n, target_dim)
but the residuals were dominated by warmup_mode_code: PySR learned the
legacy-vs-guarded contrast, not order structure.  The stratification
should have absorbed the warmup-mode signal.

Phase 3E repeats Phase 3D with the corrected stratification:
  groupby = (noise_level, n, target_dim, warmup_mode)

This residual contains only seed-to-seed variation within an otherwise
fixed protocol.  Since different seeds produce different sprinklings
with different order-theoretic invariants, this is the variance that
genuine order-theoretic features should explain — if any signal exists
at this dataset scale.

Target
------
raw_relative_drift     = warmup_delta_energy / initial_energy
stratum_key            = (noise_level, n, target_dim, warmup_mode)
stratum_mean           = mean(raw_relative_drift | stratum)
residual_relative_drift= raw_relative_drift - stratum_mean   ← target y

Singleton strata
----------------
A stratum with exactly one row produces a deterministic residual of zero.
That row carries no within-stratum information and is dropped.
The number dropped is reported in the markdown.

Panels
------
  A  design-only sanity
       X = [noise_level, warmup_mode_code, n, target_dim]
       Expected: little or no signal if residualization is correct.
  B  order + design
       X = order features + [noise_level, warmup_mode_code, n, target_dim]
  C  order + known-d, no noise / no warmup
       X = order features + [n, target_dim]
  D  order-only no-oracle
       X = [n] + order features
       Excludes noise_level, warmup_mode_code, target_dim, and every
       target-derived or embedding-dependent column.

Reproducibility
---------------
maxsize=10  (cap reduced from Phase 3D's 12 to suppress opportunistic fits).
niterations=100, deterministic seed.

Output
------
benchmarks/foundation/phase3e_pysr_residual_by_warmup_mode.csv
benchmarks/foundation/phase3e_pysr_residual_by_warmup_mode.md

This is a methodological cleanliness test.  Not a physics result.
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
    if s == "true":  return True
    if s == "false": return False
    try:    return int(s)
    except ValueError: pass
    try:    return float(s)
    except ValueError: return s


def load_csv(path: Path) -> tuple[list[str], list[dict]]:
    with path.open(encoding="utf-8") as f:
        reader = csv.DictReader(f)
        header = list(reader.fieldnames or [])
        rows   = [{k: _parse_val(v) for k, v in row.items()} for row in reader]
    return header, rows


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
    _, rows = load_csv(path)
    index: dict[tuple, dict] = {}
    for row in rows:
        if row["family"] != "minkowski":
            continue
        key = (row["family"], row["target_dim"], row["n"], row["seed"])
        index[key] = {
            "mm_dim":                row["mm_dim"],
            "midpoint_dim":          row["midpoint_dim"],
            "abs_discrepancy_mm_mp": row["abs_discrepancy_mm_midpoint"],
            "chain2_count":          row["chain2_count"],
            "chain3_count":          row["chain3_count"],
            "chain3_abundance":      row["chain3_abundance"],
        }
    return index


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

NEAR_TRUTH_LABELS = {"truth_plus_small_noise", "truth_plus_medium_noise"}
NOISE_CODE  = {"truth_plus_small_noise": 0, "truth_plus_medium_noise": 1}
WARMUP_CODE = {"legacy_warmup": 1, "guarded_warmup": 2}

REQUIRED_PHASE2F_COLUMNS = {
    "family", "target_dim", "n", "seed", "init_label", "warmup_mode",
    "initial_energy", "warmup_delta_energy",
}

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
    "A": ["noise_level", "warmup_mode_code", "n", "target_dim"],
    "B": ["noise_level", "warmup_mode_code", "n", "target_dim"] + ORDER_FEATURE_NAMES,
    "C": ["n", "target_dim"] + ORDER_FEATURE_NAMES,
    "D": ["n"] + ORDER_FEATURE_NAMES,
}

PANEL_DESCRIPTIONS: dict[str, str] = {
    "A": "design-only sanity (no order features)",
    "B": "order + design (noise, warmup, n, target_dim)",
    "C": "order + n + target_dim (no noise, no warmup)",
    "D": "order-only no-oracle (only n is a coarse design var)",
}

LEAKAGE_COLUMNS = {
    "noise_level",
    "warmup_mode",
    "warmup_mode_code",
    "target_dim",
    "initial_energy",
    "warmup_delta_energy",
    "warmup_energy_before",
    "warmup_energy_after",
    "warmup_accepted_moves",
    "warmup_rejected_moves",
    "warmup_attempted_moves",
    "final_energy",
    "delta_energy",
    "preserved_near_truth",
    "improved_energy",
    "improved_interval_rmse",
    "initial_interval_rmse",
    "final_interval_rmse",
    "initial_distance_to_truth_rms",
    "final_distance_to_truth_rms",
    "raw_relative_drift",
    "residual_relative_drift",
    "stratum_mean",
}


# ---------------------------------------------------------------------------
# Target construction
# ---------------------------------------------------------------------------

def _safe(v, default: float = 0.0) -> float:
    if v is None: return default
    try:
        f = float(v)
        return default if (math.isnan(f) or math.isinf(f)) else f
    except (TypeError, ValueError):
        return default


def validate_panel_d_features(panel_d_features: list[str]) -> None:
    leaks = [f for f in panel_d_features if f in LEAKAGE_COLUMNS]
    if leaks:
        raise RuntimeError(
            f"Panel D feature set contains leakage columns: {leaks}.  "
            "Refusing to run."
        )


def assemble_rows(
    phase2f_rows: list[dict],
    phase1d_idx: dict[tuple, dict],
    inv_idx:     dict[tuple, dict],
    initial_energy_floor: float,
) -> list[dict]:
    out = []
    skipped_small_E0 = 0
    skipped_missing_join = 0

    for r in phase2f_rows:
        if r["init_label"] not in NEAR_TRUTH_LABELS:
            continue
        if r["warmup_mode"] not in WARMUP_CODE:
            continue

        E0  = r.get("initial_energy")
        wdE = r.get("warmup_delta_energy")
        if E0 is None or wdE is None:
            continue
        E0_f, wdE_f = float(E0), float(wdE)
        if not math.isfinite(E0_f) or not math.isfinite(wdE_f):
            continue
        if abs(E0_f) < initial_energy_floor:
            skipped_small_E0 += 1
            continue

        p1key  = ("minkowski", r["target_dim"], r["n"], r["seed"])
        invkey = (r["target_dim"], r["n"], r["seed"])
        if p1key not in phase1d_idx or invkey not in inv_idx:
            skipped_missing_join += 1
            continue

        p1  = phase1d_idx[p1key]
        inv = inv_idx[invkey]

        row = {
            "noise_level":      float(NOISE_CODE[r["init_label"]]),
            "warmup_mode_code": float(WARMUP_CODE[r["warmup_mode"]]),
            "warmup_mode":      r["warmup_mode"],
            "n":                _safe(r["n"]),
            "target_dim":       _safe(r["target_dim"]),
            "seed":             r["seed"],
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
            "initial_energy":        E0_f,
            "warmup_delta_energy":   wdE_f,
            "raw_relative_drift":    wdE_f / E0_f,
        }
        out.append(row)

    print(f"  assembled {len(out)} rows "
          f"(skipped {skipped_small_E0} for |initial_energy|<{initial_energy_floor:g}, "
          f"{skipped_missing_join} missing invariants join).")
    return out


def stratify_and_residualize(rows: list[dict]) -> tuple[list[dict], dict[tuple, dict], int]:
    """Group by (noise_level, n, target_dim, warmup_mode); compute residuals.

    Drop singleton strata (residual identically zero, no within-stratum
    information).  Return (kept_rows, stratum_info, n_singletons_dropped).
    """
    groups: dict[tuple, list[dict]] = {}
    for r in rows:
        key = (r["noise_level"], r["n"], r["target_dim"], r["warmup_mode"])
        groups.setdefault(key, []).append(r)

    kept_rows: list[dict] = []
    stratum_info: dict[tuple, dict] = {}
    n_singletons = 0

    for key, members in groups.items():
        if len(members) < 2:
            n_singletons += 1
            stratum_info[key] = {
                "count": len(members),
                "mean":  float("nan"),
                "min":   float("nan"),
                "max":   float("nan"),
                "dropped": True,
            }
            continue
        vals = [m["raw_relative_drift"] for m in members]
        mean = sum(vals) / len(vals)
        for m in members:
            m["stratum_mean"] = mean
            m["residual_relative_drift"] = m["raw_relative_drift"] - mean
            kept_rows.append(m)
        stratum_info[key] = {
            "count": len(members),
            "mean":  mean,
            "min":   min(vals),
            "max":   max(vals),
            "dropped": False,
        }

    return kept_rows, stratum_info, n_singletons


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
          f"{len(y)} samples, target std={y_np.std():.4g}")

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
# Output
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
    out = []
    for idx, row in eqs.iterrows():
        eq = str(row.get("equation", ""))
        out.append({
            "panel":       panel,
            "complexity":  row.get("complexity", ""),
            "loss":        f"{row.get('loss', ''):.6g}",
            "equation":    eq,
            "is_best":     "true" if idx == best_idx else "false",
            "uses_order":  "true" if _uses_order(eq) else "false",
            "uses_design": "true" if _uses_design(eq) else "false",
        })
    return out


def write_csv(all_rows: list[dict], path: Path) -> None:
    lines = [",".join(CSV_HEADERS)]
    for r in all_rows:
        lines.append(",".join(str(r.get(h, "")) for h in CSV_HEADERS))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _best_loss(eqs: Any) -> float | None:
    if eqs is None: return None
    try:    return float(eqs["loss"].min())
    except Exception: return None


def _panel_table(eqs: Any) -> list[str]:
    lines = [
        "| complexity | loss | equation | best | order? | design? |",
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


def _verdict(
    constant_loss: float,
    panel_losses: dict[str, float | None],
    panel_eqs: dict[str, Any],
    strong_threshold: float = 0.10,
    null_threshold:   float = 0.05,
) -> tuple[str, str]:
    """Apply the decision rule the user specified."""
    def rel(p: str) -> float:
        bl = panel_losses.get(p)
        if bl is None or constant_loss <= 0: return 0.0
        return (constant_loss - bl) / constant_loss

    def uses_order_in_best(p: str) -> bool:
        eqs = panel_eqs.get(p)
        if eqs is None: return False
        try:
            best_eq = str(eqs.loc[eqs["loss"].idxmin(), "equation"])
        except Exception:
            return False
        return _uses_order(best_eq)

    rA, rB, rC, rD = (rel(p) for p in "ABCD")

    if rA > strong_threshold:
        return (
            "RESIDUALIZATION_FAILED_OR_LEAKAGE",
            f"Panel A (design-only) improves by {rA:+.1%} over constant. "
            "If residualization absorbed the design variables, A should be "
            "near zero.  Either the stratification did not fully absorb the "
            "design signal, or a target-derived column is leaking into A's "
            "features.  Investigate before interpreting B/C/D."
        )
    if rD > strong_threshold and uses_order_in_best("D"):
        return (
            "POSSIBLE_INTRA_PROTOCOL_ORDER_SIGNAL",
            f"Panel D improves by {rD:+.1%} over constant and uses order "
            "features in its best equation.  After controlling for all design "
            "variables via stratification, order-theoretic features explain "
            "part of the seed-to-seed residual.  Exploratory positive; not a "
            "physical law and needs replication with more data."
        )
    if rD < null_threshold:
        return (
            "NULL_SIGNAL",
            f"Panel D's relative improvement is {rD:+.1%}, below the {null_threshold:.0%} "
            "null threshold.  After warmup-mode-aware residualization, order-only "
            "features cannot explain the within-stratum variance.  Likely "
            "interpretation: the dataset (3 seeds × 3 dims × 2 sizes × 2 noise "
            "× 2 warmup modes) lacks sufficient intra-stratum variability for "
            "PySR to find a signal at this scale."
        )
    if rB > strong_threshold and rD <= strong_threshold:
        return (
            "SIGNAL_CONDITIONED_ON_DESIGN_VARIABLES",
            f"Panel B improves ({rB:+.1%}) but Panel D ({rD:+.1%}) does not "
            "cross the strong-signal threshold.  Order features help only "
            "in combination with design variables."
        )
    return (
        "INTERMEDIATE",
        f"Panel D's relative improvement is {rD:+.1%} (between {null_threshold:.0%} "
        f"and {strong_threshold:.0%}).  Weak suggestion of order-theoretic signal; "
        "not strong enough to interpret as positive."
    )


def write_markdown(
    n_rows_after_singletons: int,
    n_singletons_dropped: int,
    n_rows_total: int,
    initial_energy_floor: float,
    stratum_info: dict[tuple, dict],
    constant_loss: float,
    y_std: float,
    panel_eqs: dict[str, Any],
    niterations: int,
    maxsize: int,
    path: Path,
) -> None:
    panel_losses = {p: _best_loss(panel_eqs.get(p)) for p in "ABCD"}
    verdict_label, verdict_text = _verdict(constant_loss, panel_losses, panel_eqs)

    lines = [
        "# Phase 3E — PySR on warmup-mode-aware stratum residuals",
        "",
        "**Status:** methodological cleanup test, not a physics result.",
        "",
        "Phase 3D residualized by `(noise_level, n, target_dim)` but the",
        "residuals were dominated by `warmup_mode_code` — legacy vs guarded",
        "warmup differ by construction and should have been stratified out.",
        "",
        "Phase 3E uses the corrected stratification",
        "`(noise_level, n, target_dim, warmup_mode)`, so the residual",
        "captures only seed-to-seed variation within an otherwise fixed",
        "protocol.  Different seeds give different sprinklings with different",
        "order-theoretic invariants — this is the variance order features",
        "should explain if any signal exists at this dataset scale.",
        "",
        "## Verdict (automatic)",
        "",
        f"**{verdict_label}**",
        "",
        verdict_text,
        "",
        "## Target definition",
        "",
        "```",
        "raw_relative_drift     = warmup_delta_energy / initial_energy",
        "stratum_key            = (noise_level, n, target_dim, warmup_mode)",
        "stratum_mean           = mean(raw_relative_drift | stratum)",
        "residual_relative_drift= raw_relative_drift - stratum_mean   ← target y",
        "```",
        "",
        f"Rows with `|initial_energy| < {initial_energy_floor:g}` are excluded "
        "to avoid division blow-ups.",
        "",
        "Singleton strata (count = 1) are dropped: the residual would be "
        "identically zero and carry no within-stratum information.",
        "",
        "## Strata",
        "",
        f"- Rows after E0 filter: {n_rows_total}",
        f"- Strata: {len(stratum_info)}",
        f"- Singleton strata dropped: {n_singletons_dropped}",
        f"- Rows kept for regression: {n_rows_after_singletons}",
        f"- Target std (kept rows): {y_std:.4g}",
        f"- Constant-predictor loss (variance of y): **{constant_loss:.6g}**",
        "",
        "### Per-stratum counts",
        "",
        "| noise | n | target_dim | warmup_mode | count | mean raw | min | max | kept? |",
        "| :---: | ---: | :---: | --- | ---: | ---: | ---: | ---: | :---: |",
    ]
    for key in sorted(stratum_info.keys()):
        nl, n, td, wm = key
        info = stratum_info[key]
        nl_label = "small" if nl == 0 else "medium"
        kept = "—" if info["dropped"] else "✓"
        if info["dropped"]:
            lines.append(
                f"| {nl_label} | {int(n)} | {int(td)} | {wm} | {info['count']}"
                f" | — | — | — | {kept} |"
            )
        else:
            lines.append(
                f"| {nl_label} | {int(n)} | {int(td)} | {wm} | {info['count']}"
                f" | {info['mean']:+.4g} | {info['min']:+.4g} | {info['max']:+.4g} | {kept} |"
            )

    lines += [
        "",
        "## Panel definitions",
        "",
        "| panel | description | features |",
        "| --- | --- | --- |",
    ]
    for p in "ABCD":
        lines.append(
            f"| {p} | {PANEL_DESCRIPTIONS[p]} | `{', '.join(PANELS[p])}` |"
        )

    lines += [
        "",
        "**Panel D leakage guard** (verified at build time):",
        "",
        f"  `{', '.join(sorted(LEAKAGE_COLUMNS))}`",
        "",
        "## Summary: best loss per panel",
        "",
        "| panel | description | best loss | Δ vs constant | rel. Δ |",
        "| --- | --- | ---: | ---: | ---: |",
    ]
    for p in "ABCD":
        bl = panel_losses[p]
        if bl is None:
            lines.append(f"| {p} | {PANEL_DESCRIPTIONS[p]} | — | — | — |")
        else:
            delta = constant_loss - bl
            rel = delta / constant_loss if constant_loss > 0 else 0.0
            lines.append(
                f"| {p} | {PANEL_DESCRIPTIONS[p]} | {bl:.4g} | "
                f"{delta:+.4g} | {rel:+.2%} |"
            )

    lines += [
        "",
        f"## Pareto fronts (per panel, complexity ≤ {maxsize})",
        "",
    ]
    for p in "ABCD":
        lines += [
            f"### Panel {p} — {PANEL_DESCRIPTIONS[p]}",
            "",
            *_panel_table(panel_eqs.get(p)),
            "",
        ]

    lines += [
        "## Decision rule",
        "",
        "Improvement thresholds:",
        "- *strong* = best loss is at least **10%** below constant baseline",
        "- *null*   = best loss is less than **5%** below constant baseline",
        "",
        "Applied in this order:",
        "",
        "1. **RESIDUALIZATION_FAILED_OR_LEAKAGE**: A improves strongly.",
        "   → residualization or feature set is broken; do not interpret.",
        "2. **POSSIBLE_INTRA_PROTOCOL_ORDER_SIGNAL**: D improves strongly AND",
        "   uses order features.  → exploratory positive.",
        "3. **NULL_SIGNAL**: D's improvement below null threshold.",
        "   → no detectable order signal at this dataset scale.",
        "4. **SIGNAL_CONDITIONED_ON_DESIGN_VARIABLES**: B strong, D not strong.",
        "5. **INTERMEDIATE**: D between null and strong thresholds.",
        "",
        "## Reproducibility",
        "",
        f"- PySR iterations per panel: {niterations}",
        f"- maxsize cap: {maxsize}",
        "- random_state=1959, parallelism=serial (deterministic).",
        f"- initial_energy_floor: {initial_energy_floor:g}",
        "",
        "Regenerate via `make regen-phase3e`.",
        "Source: `tools/build_phase3e_pysr_residual_by_warmup_mode.py`.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--niterations", type=int, default=100)
    parser.add_argument("--maxsize",     type=int, default=10)
    parser.add_argument("--initial-energy-floor", type=float, default=1e-4)
    args = parser.parse_args()

    # ---- Required-column check (fail fast) ----
    p2f_header, p2f_rows = load_csv(FOUNDATION / "phase2f_guarded_warmup_probe.csv")
    missing = REQUIRED_PHASE2F_COLUMNS - set(p2f_header)
    if missing:
        print(
            f"ERROR: phase2f missing required columns: {sorted(missing)}\n"
            f"Available columns: {p2f_header}",
            file=sys.stderr,
        )
        sys.exit(1)

    validate_panel_d_features(PANELS["D"])

    phase1d_idx = load_phase1d(FOUNDATION / "phase1d_structural_atlas.csv")
    inv_idx     = load_invariants_json(FOUNDATION / "invariants.json")

    rows_all = assemble_rows(p2f_rows, phase1d_idx, inv_idx, args.initial_energy_floor)
    if not rows_all:
        print("ERROR: zero usable rows after filtering.", file=sys.stderr)
        sys.exit(2)

    rows_kept, stratum_info, n_singletons = stratify_and_residualize(rows_all)
    print(f"  strata={len(stratum_info)} "
          f"(singletons dropped={n_singletons}), kept={len(rows_kept)} rows.")

    if not rows_kept:
        print("ERROR: all strata are singletons. Cannot compute residuals.",
              file=sys.stderr)
        sys.exit(3)

    y = [r["residual_relative_drift"] for r in rows_kept]
    y_mean = sum(y) / len(y)
    y_std  = math.sqrt(sum((yi - y_mean)**2 for yi in y) / len(y))
    constant_loss = sum((yi - y_mean)**2 for yi in y) / len(y)
    print(f"  target std={y_std:.4g}, constant-predictor loss={constant_loss:.6g}")

    panel_eqs: dict[str, Any] = {}
    all_csv_rows: list[dict] = []

    for p, feats in PANELS.items():
        X = [[_safe(r[f]) for f in feats] for r in rows_kept]
        eqs = run_pysr(X, y, feats, args.niterations, args.maxsize, p)
        panel_eqs[p] = eqs
        all_csv_rows.extend(equations_to_rows(eqs, p))

    FOUNDATION.mkdir(parents=True, exist_ok=True)
    write_csv(all_csv_rows, FOUNDATION / "phase3e_pysr_residual_by_warmup_mode.csv")
    write_markdown(
        n_rows_after_singletons=len(rows_kept),
        n_singletons_dropped=n_singletons,
        n_rows_total=len(rows_all),
        initial_energy_floor=args.initial_energy_floor,
        stratum_info=stratum_info,
        constant_loss=constant_loss,
        y_std=y_std,
        panel_eqs=panel_eqs,
        niterations=args.niterations,
        maxsize=args.maxsize,
        path=FOUNDATION / "phase3e_pysr_residual_by_warmup_mode.md",
    )

    print("\n--- Phase 3E summary ---")
    for p in "ABCD":
        bl = _best_loss(panel_eqs.get(p))
        if bl is None:
            print(f"  Panel {p}: (no result)")
            continue
        rel = (constant_loss - bl) / constant_loss if constant_loss > 0 else 0.0
        eqs = panel_eqs[p]
        try:
            best_eq = str(eqs.loc[eqs["loss"].idxmin(), "equation"])
        except Exception:
            best_eq = "—"
        flag = "ORDER" if _uses_order(best_eq) else "design/none"
        print(f"  Panel {p}: loss={bl:.4g}  (rel Δ={rel:+.2%})  [{flag}]")
        print(f"           {best_eq[:90]}")

    print(f"\nOutput: {FOUNDATION / 'phase3e_pysr_residual_by_warmup_mode.csv'}")
    print(f"        {FOUNDATION / 'phase3e_pysr_residual_by_warmup_mode.md'}")


if __name__ == "__main__":
    main()
