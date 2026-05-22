#!/usr/bin/env python3
"""Phase 3D — PySR on the stratum-residual target.

Phase 3C showed that ``preserved_near_truth`` is dominated by ``noise_level``:
order-only panels could not improve over the constant baseline.  Phase 3D
keeps the dataset unchanged but switches to a continuous target designed
to expose *intra-stratum* variability — the part of the warmup drift that
``noise_level`` and ``n`` and ``target_dim`` cannot explain by themselves.

Target construction
-------------------
1. raw_relative_drift = warmup_delta_energy / initial_energy
   (computed only on rows where both columns are present and the
    denominator is large enough — fail-fast if columns are missing)

2. Group rows by stratum (noise_level, n, target_dim) and compute the
   per-stratum mean of raw_relative_drift.

3. residual_relative_drift = raw_relative_drift - stratum_mean
   This is the within-stratum deviation: the part of the drift that
   stratification (i.e. design variables) cannot explain.

If order-theoretic features explain the residual, then there is signal
above and beyond the experimental design.  Phase 3C's verdict was
"order-only no-oracle is indistinguishable from constant"; Phase 3D
asks the same question with a target that gives intra-stratum variance
a chance to surface.

Panels (parallel to Phase 3C):
  A  noise-only         X = [noise_level]
  B  order + noise      X = [noise_level, warmup_mode_code, n, target_dim]
                            + order features
  C  order + known-d    X = [n, target_dim, warmup_mode_code] + order
                            (noise_level excluded)
  D  order-only         X = [n] + order features
                            (no noise_level, no target_dim, no warmup_mode_code,
                             no initial_energy, no warmup_delta_energy,
                             no preserved_near_truth)

Source rows
-----------
phase2f only.  phase2e lacks the warmup_delta_energy column.
Filter: init_label in {truth_plus_small_noise, truth_plus_medium_noise}
        warmup_mode in {legacy_warmup, guarded_warmup}
        (skip_warmup is excluded because warmup_delta_energy is identically
         zero by construction in that branch and would contaminate the
         stratum mean with a zero point that has no warmup dynamics.)

Reproducibility
---------------
maxsize=12, niterations=100, deterministic seed.  Same caps as Phase 3C.

Output
------
benchmarks/foundation/phase3d_pysr_residual_target.csv
benchmarks/foundation/phase3d_pysr_residual_target.md

This is an exploratory diagnostic, not a physics result.  The .md report
labels it as such.
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
# Loaders (identical pattern to Phase 3B/3C)
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
    """Return (header, rows) so that missing-column failures can be precise."""
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
    "A": ["noise_level"],
    "B": ["noise_level", "warmup_mode_code", "n", "target_dim"] + ORDER_FEATURE_NAMES,
    "C": ["warmup_mode_code", "n", "target_dim"] + ORDER_FEATURE_NAMES,
    "D": ["n"] + ORDER_FEATURE_NAMES,
}

PANEL_DESCRIPTIONS: dict[str, str] = {
    "A": "noise-only baseline",
    "B": "order + noise + warmup_mode + target_dim",
    "C": "order + warmup_mode + target_dim (no noise_level)",
    "D": "order-only no-oracle (only n is a coarse design var)",
}

# Columns that are forbidden in Panel D: leak the target or its components.
LEAKAGE_COLUMNS = {
    "noise_level",
    "target_dim",
    "warmup_mode_code",
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


def assemble_rows(
    phase2f_rows: list[dict],
    phase1d_idx: dict[tuple, dict],
    inv_idx:     dict[tuple, dict],
    initial_energy_floor: float,
) -> list[dict]:
    """Build enriched rows with raw_relative_drift.  Apply filters."""
    out = []
    skipped_small_E0 = 0
    skipped_missing_join = 0

    for r in phase2f_rows:
        if r["init_label"] not in NEAR_TRUTH_LABELS:
            continue
        if r["warmup_mode"] not in WARMUP_CODE:
            continue

        E0     = r.get("initial_energy")
        wdE    = r.get("warmup_delta_energy")
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
            "initial_energy":        E0_f,
            "warmup_delta_energy":   wdE_f,
            "raw_relative_drift":    wdE_f / E0_f,
        }
        out.append(row)

    print(f"  assembled {len(out)} rows "
          f"(skipped {skipped_small_E0} for |initial_energy|<{initial_energy_floor:g}, "
          f"{skipped_missing_join} missing invariants join).")
    return out


def add_stratum_residuals(rows: list[dict]) -> dict[tuple, dict]:
    """Compute per-stratum mean of raw_relative_drift, attach residuals.

    Stratum key: (noise_level, n, target_dim).
    Returns a stratum_info dict for the markdown summary.
    """
    strata: dict[tuple, list[float]] = {}
    for r in rows:
        key = (r["noise_level"], r["n"], r["target_dim"])
        strata.setdefault(key, []).append(r["raw_relative_drift"])

    stratum_info: dict[tuple, dict] = {}
    for key, vals in strata.items():
        m = sum(vals) / len(vals)
        stratum_info[key] = {
            "count": len(vals),
            "mean":  m,
            "min":   min(vals),
            "max":   max(vals),
        }

    for r in rows:
        key = (r["noise_level"], r["n"], r["target_dim"])
        r["stratum_mean"] = stratum_info[key]["mean"]
        r["residual_relative_drift"] = r["raw_relative_drift"] - stratum_info[key]["mean"]

    return stratum_info


# ---------------------------------------------------------------------------
# Panel feature integrity check
# ---------------------------------------------------------------------------

def validate_panel_d_features(panel_d_features: list[str]) -> None:
    """Panel D must not contain any leakage column. Fail fast if it does."""
    leaks = [f for f in panel_d_features if f in LEAKAGE_COLUMNS]
    if leaks:
        raise RuntimeError(
            f"Panel D feature set contains leakage columns: {leaks}.  "
            "Refusing to run."
        )


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
          f"{len(y)} samples, target mean={y_np.mean():+.4g}, std={y_np.std():.4g}")

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
    threshold_rel: float = 0.10,
) -> tuple[str, str]:
    """Decide which automatic verdict applies based on relative improvement.

    threshold_rel = 0.10: a panel "improves clearly" if best_loss is at least
    10% below the constant-predictor loss.
    """
    def improves(p: str) -> bool:
        bl = panel_losses.get(p)
        return bl is not None and bl < constant_loss * (1.0 - threshold_rel)

    def uses_order_in_best(p: str) -> bool:
        eqs = panel_eqs.get(p)
        if eqs is None: return False
        try:
            best_eq = str(eqs.loc[eqs["loss"].idxmin(), "equation"])
        except Exception:
            return False
        return _uses_order(best_eq)

    a, b, c, d = (improves(p) for p in "ABCD")
    d_uses_order = uses_order_in_best("D")

    if d and d_uses_order:
        return (
            "INTRA_STRATUM_ORDER_SIGNAL",
            "Panel D improves clearly over the constant baseline and its best "
            "equation uses order-theoretic features.  Intra-stratum variability "
            "is partially explained by combinatorial invariants.  This is a "
            "diagnostic positive, not a physical law."
        )
    if not (a or b or c or d):
        return (
            "NO_DETECTABLE_SIGNAL",
            "No panel improves meaningfully over the constant baseline.  "
            "With this dataset, residualizing the target did not surface "
            "intra-stratum signal in any panel."
        )
    if c and not d:
        return (
            "SIGNAL_CONDITIONED_ON_DESIGN_VARIABLES",
            "Panel C improves but Panel D does not.  The signal needs "
            "target_dim and/or warmup_mode_code (design variables) even after "
            "removing noise_level.  Order features by themselves cannot "
            "explain the residual."
        )
    if b and not c and not d:
        return (
            "SIGNAL_CONDITIONED_ON_NOISE_LEVEL",
            "Panel B improves but Panels C and D do not.  Order features help "
            "only in combination with noise_level.  Same status as Phase 3C: "
            "no autonomous rule."
        )
    return (
        "MIXED",
        "Heterogeneous panel improvements; no single clean verdict."
    )


def write_markdown(
    n_rows: int,
    initial_energy_floor: float,
    stratum_info: dict[tuple, dict],
    constant_loss: float,
    panel_eqs: dict[str, Any],
    niterations: int,
    maxsize: int,
    path: Path,
) -> None:
    panel_losses = {p: _best_loss(panel_eqs.get(p)) for p in "ABCD"}
    verdict_label, verdict_text = _verdict(constant_loss, panel_losses, panel_eqs)

    lines = [
        "# Phase 3D — PySR on the stratum-residual target",
        "",
        "**Status:** exploratory diagnostic.  Not a physics result.",
        "",
        "Phase 3C showed that with the binary target `preserved_near_truth`,",
        "order-only panels could not improve over the majority-class constant.",
        "Phase 3D keeps the dataset unchanged but switches to a continuous",
        "target — the within-stratum residual of warmup-relative drift —",
        "designed to expose whatever variability stratification cannot absorb.",
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
        "stratum_key            = (noise_level, n, target_dim)",
        "stratum_mean           = mean(raw_relative_drift | stratum)",
        "residual_relative_drift= raw_relative_drift - stratum_mean   ← target y",
        "```",
        "",
        f"Rows with |initial_energy| < {initial_energy_floor:g} are excluded to "
        "avoid division blow-ups (documented; no row is silently kept).",
        "",
        "## Sample sizes per stratum",
        "",
        "| noise_level | n | target_dim | count | mean raw | min | max |",
        "| :---: | ---: | :---: | ---: | ---: | ---: | ---: |",
    ]
    for key in sorted(stratum_info.keys()):
        nl, n, td = key
        info = stratum_info[key]
        nl_label = "small" if nl == 0 else "medium"
        lines.append(
            f"| {nl_label} | {int(n)} | {int(td)} | {info['count']} "
            f"| {info['mean']:+.4g} | {info['min']:+.4g} | {info['max']:+.4g} |"
        )

    lines += [
        "",
        f"Total rows: {n_rows}.  Source: phase2f_guarded_warmup_probe.csv",
        "(phase2e excluded because it lacks the `warmup_delta_energy` column).",
        "",
        f"Filter: init_label ∈ {{truth_plus_small_noise, truth_plus_medium_noise}}, "
        "warmup_mode ∈ {legacy_warmup, guarded_warmup}.",
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
        "**Panel D leakage guard:** the following columns are excluded from",
        "Panel D and verified at build time. If any reappears, the run fails:",
        "",
        f"  `{', '.join(sorted(LEAKAGE_COLUMNS))}`",
        "",
        "## Summary: best loss per panel",
        "",
        f"Constant-predictor loss (variance of y): **{constant_loss:.6g}**",
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
        "## Pareto fronts (per panel, complexity ≤ {maxsize})".format(maxsize=maxsize),
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
        "## Decision rule (applied above to produce verdict)",
        "",
        "Improvement threshold: a panel \"improves clearly\" if best loss is",
        "at least 10% below the constant-predictor loss.",
        "",
        "- **INTRA_STRATUM_ORDER_SIGNAL**: D improves AND best equation in D uses",
        "  order features.  → autonomous order-theoretic signal.",
        "- **SIGNAL_CONDITIONED_ON_DESIGN_VARIABLES**: C improves but D does not.",
        "  → signal needs target_dim and/or warmup_mode_code.",
        "- **SIGNAL_CONDITIONED_ON_NOISE_LEVEL**: B improves but C and D do not.",
        "  → signal needs noise_level specifically.",
        "- **NO_DETECTABLE_SIGNAL**: no panel improves meaningfully.",
        "- **MIXED**: heterogeneous, no single conclusion.",
        "",
        "## Reproducibility",
        "",
        f"- PySR niterations per panel: {niterations}",
        f"- maxsize cap: {maxsize}",
        "- random_state=1959, parallelism=serial (deterministic).",
        f"- initial_energy_floor: {initial_energy_floor:g}",
        "",
        "Regenerate via `make regen-phase3d`.",
        "Source: `tools/build_phase3d_pysr_residual_target.py`.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--niterations", type=int, default=100,
                        help="PySR iterations per panel (default 100, same as Phase 3C)")
    parser.add_argument("--maxsize",     type=int, default=12,
                        help="PySR complexity cap (default 12, same as Phase 3C)")
    parser.add_argument("--initial-energy-floor", type=float, default=1e-4,
                        help="Skip rows with |initial_energy| below this; default 1e-4")
    args = parser.parse_args()

    # ---- Required-column check (fail fast) ----
    p2f_header, p2f_rows = load_csv(FOUNDATION / "phase2f_guarded_warmup_probe.csv")
    missing = REQUIRED_PHASE2F_COLUMNS - set(p2f_header)
    if missing:
        print(
            f"ERROR: phase2f is missing required columns: {sorted(missing)}\n"
            f"Available columns: {p2f_header}",
            file=sys.stderr,
        )
        sys.exit(1)

    # ---- Panel-D leakage guard (fail fast) ----
    validate_panel_d_features(PANELS["D"])

    # ---- Load auxiliary indices ----
    phase1d_idx = load_phase1d(FOUNDATION / "phase1d_structural_atlas.csv")
    inv_idx     = load_invariants_json(FOUNDATION / "invariants.json")

    # ---- Assemble rows + residualize target ----
    rows = assemble_rows(p2f_rows, phase1d_idx, inv_idx, args.initial_energy_floor)
    if not rows:
        print("ERROR: zero usable rows after filtering. Cannot run PySR.",
              file=sys.stderr)
        sys.exit(2)

    stratum_info = add_stratum_residuals(rows)
    print(f"  {len(stratum_info)} strata, {len(rows)} rows total.")

    # ---- Build target ----
    y = [r["residual_relative_drift"] for r in rows]
    constant_loss = sum((yi - sum(y)/len(y))**2 for yi in y) / len(y)
    print(f"  Constant-predictor loss (variance of y): {constant_loss:.6g}")

    # ---- Run all four panels ----
    panel_eqs: dict[str, Any] = {}
    all_csv_rows: list[dict] = []

    for p, feats in PANELS.items():
        X = [[_safe(r[f]) for f in feats] for r in rows]
        eqs = run_pysr(X, y, feats, args.niterations, args.maxsize, p)
        panel_eqs[p] = eqs
        all_csv_rows.extend(equations_to_rows(eqs, p))

    # ---- Write outputs ----
    FOUNDATION.mkdir(parents=True, exist_ok=True)
    write_csv(all_csv_rows, FOUNDATION / "phase3d_pysr_residual_target.csv")
    write_markdown(
        n_rows=len(rows),
        initial_energy_floor=args.initial_energy_floor,
        stratum_info=stratum_info,
        constant_loss=constant_loss,
        panel_eqs=panel_eqs,
        niterations=args.niterations,
        maxsize=args.maxsize,
        path=FOUNDATION / "phase3d_pysr_residual_target.md",
    )

    print("\n--- Phase 3D summary ---")
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

    print(f"\nOutput: {FOUNDATION / 'phase3d_pysr_residual_target.csv'}")
    print(f"        {FOUNDATION / 'phase3d_pysr_residual_target.md'}")


if __name__ == "__main__":
    main()
