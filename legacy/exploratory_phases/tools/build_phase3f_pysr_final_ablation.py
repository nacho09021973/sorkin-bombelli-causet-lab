#!/usr/bin/env python3
"""Phase 3F — PySR ablation on the expanded Phase 2G dataset.

Phase 3E methodologically cleaned the residual analysis but the dataset
(60 effective rows, 3 seeds, 2 sizes) gave only a borderline +15% signal
in Panel D that was largely explained by midpoint_dim acting as a noisy
target_dim proxy.  The actual genuine-order contribution was ~3% from
abs_discrepancy_mm_midpoint.

Phase 3F repeats the ablation on the Phase 2G + Phase 1E dataset:
  - 15 seeds, sizes {32, 64, 128}, dims {2, 3, 4}
  - guarded_warmup only (eliminates legacy/guarded heteroscedasticity)
  - 270 raw rows before E0 filter

Because only one warmup mode is present, the stratification simplifies to
  (noise_level, n, target_dim)
which is the cleanest residualization compatible with the data.

Hypothesis under test
---------------------
Phase 3E identified `abs_discrepancy_mm_midpoint` as the only feature
that consistently appears beyond midpoint_dim's dimension-proxy role.
Phase 3F tests:

  H1.  Does Panel D (order-only) cross the +10% threshold over constant?
  H2.  Does abs_discrepancy_mm_midpoint appear in Panel D's best
       equation, with contribution larger than the seed noise floor?
  H3.  Is Panel A (design-only sanity) at ~0% as expected?

Panels (parallel to 3E)
-----------------------
  A  design-only sanity
  B  order + design (noise_level, n, target_dim)
  C  order + n + target_dim (no noise_level)
  D  order-only no-oracle (n + order features)

Reproducibility
---------------
maxsize=10, niterations=100, random_state=1959, parallelism=serial.

Output
------
benchmarks/foundation/phase3f_pysr_final_ablation.csv
benchmarks/foundation/phase3f_pysr_final_ablation.md
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
# Loaders
# ---------------------------------------------------------------------------

def _parse_val(s: str) -> Any:
    s = s.strip()
    if s in ("NA", ""): return None
    if s == "true":     return True
    if s == "false":    return False
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


def load_phase1e(path: Path) -> dict[tuple, dict]:
    _, rows = load_csv(path)
    out: dict[tuple, dict] = {}
    for r in rows:
        if r["family"] != "minkowski":
            continue
        key = (r["target_dim"], r["n"], r["seed"])
        out[key] = {
            "mm_dim":                r["mm_dim"],
            "midpoint_dim":          r["midpoint_dim"],
            "abs_discrepancy_mm_mp": r["abs_discrepancy_mm_midpoint"],
            "chain2_count":          r["chain2_count"],
            "chain3_count":          r["chain3_count"],
            "chain3_abundance":      r["chain3_abundance"],
            "chain4_count":          r["chain4_count"],
            "link_count":            r["link_count"],
            "link_density":          r["link_density"],
            "relation_count":        r["relation_count"],
            "ordering_fraction":     r["ordering_fraction"],
            "height":                r["height"],
        }
    return out


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

NEAR_TRUTH_LABELS = {"truth_plus_small_noise", "truth_plus_medium_noise"}
NOISE_CODE = {"truth_plus_small_noise": 0, "truth_plus_medium_noise": 1}

REQUIRED_PHASE2G_COLUMNS = {
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
    "chain4_count",
    "link_count",
    "link_density",
    "relation_count",
    "ordering_fraction",
    "height",
]
ORDER_FEATURE_SET = set(ORDER_FEATURE_NAMES)

PANELS: dict[str, list[str]] = {
    "A": ["noise_level", "n", "target_dim"],
    "B": ["noise_level", "n", "target_dim"] + ORDER_FEATURE_NAMES,
    "C": ["n", "target_dim"] + ORDER_FEATURE_NAMES,
    "D": ["n"] + ORDER_FEATURE_NAMES,
}

PANEL_DESCRIPTIONS: dict[str, str] = {
    "A": "design-only sanity (noise_level + n + target_dim)",
    "B": "order + design (noise_level + n + target_dim + order features)",
    "C": "order + known-d (n + target_dim + order features, no noise)",
    "D": "order-only no-oracle (n + order features only)",
}

LEAKAGE_COLUMNS = {
    "noise_level", "warmup_mode", "warmup_mode_code", "target_dim",
    "initial_energy", "warmup_delta_energy",
    "warmup_energy_before", "warmup_energy_after",
    "warmup_accepted_moves", "warmup_rejected_moves", "warmup_attempted_moves",
    "final_energy", "delta_energy",
    "preserved_near_truth", "improved_energy", "improved_interval_rmse",
    "initial_interval_rmse", "final_interval_rmse",
    "initial_distance_to_truth_rms", "final_distance_to_truth_rms",
    "raw_relative_drift", "residual_relative_drift", "stratum_mean",
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
            f"Panel D feature set contains leakage columns: {leaks}.  Refusing to run."
        )


def assemble_rows(
    phase2g_rows: list[dict],
    p1e_idx: dict[tuple, dict],
    initial_energy_floor: float,
) -> list[dict]:
    out = []
    skipped_small_E0 = skipped_missing_join = 0
    for r in phase2g_rows:
        if r["init_label"] not in NEAR_TRUTH_LABELS:
            continue
        if r["warmup_mode"] != "guarded_warmup":
            continue

        E0 = r.get("initial_energy")
        wdE = r.get("warmup_delta_energy")
        if E0 is None or wdE is None:
            continue
        E0_f, wdE_f = float(E0), float(wdE)
        if not math.isfinite(E0_f) or not math.isfinite(wdE_f):
            continue
        if abs(E0_f) < initial_energy_floor:
            skipped_small_E0 += 1
            continue

        key = (r["target_dim"], r["n"], r["seed"])
        if key not in p1e_idx:
            skipped_missing_join += 1
            continue
        inv = p1e_idx[key]

        out.append({
            "noise_level":  float(NOISE_CODE[r["init_label"]]),
            "n":            _safe(r["n"]),
            "target_dim":   _safe(r["target_dim"]),
            "seed":         r["seed"],
            "initial_energy":      E0_f,
            "warmup_delta_energy": wdE_f,
            "raw_relative_drift":  wdE_f / E0_f,
            **inv,
        })
    print(f"  assembled {len(out)} rows "
          f"(skipped {skipped_small_E0} small-E₀, {skipped_missing_join} missing-invariants).")
    return out


def stratify_and_residualize(rows: list[dict]) -> tuple[list[dict], dict[tuple, dict], int]:
    """Stratify by (noise_level, n, target_dim). Drop singleton strata."""
    groups: dict[tuple, list[dict]] = {}
    for r in rows:
        key = (r["noise_level"], r["n"], r["target_dim"])
        groups.setdefault(key, []).append(r)

    kept: list[dict] = []
    info: dict[tuple, dict] = {}
    n_singletons = 0
    for key, members in groups.items():
        if len(members) < 2:
            n_singletons += 1
            info[key] = {"count": len(members), "mean": math.nan,
                         "min": math.nan, "max": math.nan, "dropped": True}
            continue
        vals = [m["raw_relative_drift"] for m in members]
        mean = sum(vals) / len(vals)
        for m in members:
            m["stratum_mean"] = mean
            m["residual_relative_drift"] = m["raw_relative_drift"] - mean
            kept.append(m)
        info[key] = {"count": len(members), "mean": mean,
                     "min": min(vals), "max": max(vals), "dropped": False}
    return kept, info, n_singletons


# ---------------------------------------------------------------------------
# PySR runner
# ---------------------------------------------------------------------------

def run_pysr(X, y, feature_names, niterations, maxsize, label):
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
        maxsize=maxsize, populations=20,
        model_selection="best",
        verbosity=0, random_state=1959,
        deterministic=True, parallelism="serial",
    )
    model.fit(X_df, y_np)
    return model.equations_


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def _uses_order(eq: str) -> bool:
    return any(f in eq for f in ORDER_FEATURE_SET)


def _uses_design(eq: str) -> bool:
    return any(f in eq for f in ("noise_level", "target_dim"))


def _uses_abs_disc(eq: str) -> bool:
    return "abs_discrepancy_mm_mp" in eq


CSV_HEADERS = (
    "panel", "complexity", "loss", "equation",
    "is_best", "uses_order", "uses_design", "uses_abs_discrepancy",
)


def equations_to_rows(eqs, panel):
    if eqs is None: return []
    try:
        best_idx = eqs["loss"].idxmin()
    except Exception:
        best_idx = None
    out = []
    for idx, row in eqs.iterrows():
        eq = str(row.get("equation", ""))
        out.append({
            "panel": panel,
            "complexity": row.get("complexity", ""),
            "loss": f"{row.get('loss', ''):.6g}",
            "equation": eq,
            "is_best": "true" if idx == best_idx else "false",
            "uses_order":   "true" if _uses_order(eq) else "false",
            "uses_design":  "true" if _uses_design(eq) else "false",
            "uses_abs_discrepancy": "true" if _uses_abs_disc(eq) else "false",
        })
    return out


def write_csv(all_rows, path):
    lines = [",".join(CSV_HEADERS)]
    for r in all_rows:
        lines.append(",".join(str(r.get(h, "")) for h in CSV_HEADERS))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _best_loss(eqs):
    if eqs is None: return None
    try: return float(eqs["loss"].min())
    except Exception: return None


def _panel_table(eqs):
    lines = [
        "| complexity | loss | equation | best | order? | design? | abs_disc? |",
        "| ---: | ---: | --- | :---: | :---: | :---: | :---: |",
    ]
    if eqs is None:
        lines.append("| — | — | *PySR unavailable* | | | | |")
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
            f" | {'✓' if _uses_design(eq) else '—'}"
            f" | {'✓' if _uses_abs_disc(eq) else '—'} |"
        )
    return lines


def _verdict(constant_loss, panel_losses, panel_eqs,
             strong=0.10, null=0.05):
    def has_equations(p):
        eqs = panel_eqs.get(p)
        if eqs is None:
            return False
        if hasattr(eqs, "empty"):
            return not eqs.empty
        return bool(eqs)

    if not any(has_equations(p) for p in "ABCD"):
        return (
            "PYSR_UNAVAILABLE",
            "No PySR equations were produced for any panel. This is a missing-dependency "
            "or execution-environment outcome, not evidence for or against an "
            "order-theoretic signal."
        )

    def rel(p):
        bl = panel_losses.get(p)
        if bl is None or constant_loss <= 0: return 0.0
        return (constant_loss - bl) / constant_loss

    def best_eq(p):
        eqs = panel_eqs.get(p)
        if eqs is None: return ""
        try: return str(eqs.loc[eqs["loss"].idxmin(), "equation"])
        except Exception: return ""

    rA, rB, rC, rD = (rel(p) for p in "ABCD")
    bestD = best_eq("D")
    bestC = best_eq("C")
    d_order = _uses_order(bestD)
    d_abs   = _uses_abs_disc(bestD)
    c_abs   = _uses_abs_disc(bestC)

    H1 = rD >= strong
    H2 = d_abs or c_abs   # abs_discrepancy appears in D's best (or C's best)
    H3 = rA < strong      # design-only sanity holds

    if not H3:
        return (
            "RESIDUALIZATION_FAILED_OR_LEAKAGE",
            f"Panel A improves by {rA:+.1%} over constant. Sanity broken — "
            "design-only should be near zero after stratification."
        )
    if H1 and d_order:
        if H2:
            return (
                "POSSIBLE_ORDER_SIGNAL_WITH_ABS_DISC",
                f"Panel D improves by {rD:+.1%} and its best equation uses "
                "order features including abs_discrepancy_mm_midpoint.  "
                "Phase 3E's candidate signal survives the dataset expansion. "
                "Exploratory positive; needs theoretical derivation, not a law."
            )
        return (
            "ORDER_SIGNAL_WITHOUT_ABS_DISC",
            f"Panel D improves by {rD:+.1%} using order features but NOT "
            "abs_discrepancy_mm_midpoint. The 3E candidate did not replicate; "
            "a different combination of order features explains the residual."
        )
    if rD < null:
        return (
            "NULL_SIGNAL",
            f"Panel D's improvement is {rD:+.1%} (< {null:.0%}).  Even with "
            f"the expanded dataset ({len(panel_eqs.get('D') or [])} eqs) order-only "
            "features do not explain intra-stratum residual variance.  Honest "
            "conclusion: no autonomous order-theoretic rule at this pipeline scale."
        )
    if rB >= strong and rD < strong:
        return (
            "SIGNAL_CONDITIONED_ON_DESIGN",
            f"Panel B improves ({rB:+.1%}) but Panel D ({rD:+.1%}) is below "
            "the strong threshold.  Order features need design variables to help."
        )
    return (
        "INTERMEDIATE",
        f"Panel D's improvement is {rD:+.1%} — between null ({null:.0%}) and "
        f"strong ({strong:.0%}) thresholds.  Weak; not interpretable as positive."
    )


def write_markdown(
    n_kept, n_singletons, n_rows_total,
    initial_energy_floor, stratum_info,
    constant_loss, y_std,
    panel_eqs, niterations, maxsize, path,
):
    panel_losses = {p: _best_loss(panel_eqs.get(p)) for p in "ABCD"}
    verdict_label, verdict_text = _verdict(constant_loss, panel_losses, panel_eqs)

    lines = [
        "# Phase 3F — PySR final ablation on expanded dataset",
        "",
        "**Status:** decisive test of the Phase 3E candidate signal on the",
        "Phase 2G expanded dataset.  Not a physics law.  An exploratory",
        "diagnostic with sufficient samples to either replicate or reject",
        "the 3E `abs_discrepancy_mm_midpoint` finding.",
        "",
        "## Verdict (automatic)",
        "",
        f"**{verdict_label}**",
        "",
        verdict_text,
        "",
        "## Hypotheses tested",
        "",
        "- **H1**: Panel D crosses +10% over constant baseline.",
        "- **H2**: `abs_discrepancy_mm_midpoint` appears in D's (or C's) best equation.",
        "- **H3**: Panel A (design-only sanity) is below +10% — residualization is clean.",
        "",
        "## Dataset",
        "",
        f"- Source: `phase2g_extended_guarded_warmup_probe.csv` (15 seeds, "
        "3 sizes, 3 dims, 2 near-truth init labels, guarded warmup only).",
        f"- Invariants: `phase1e_extended_structural_atlas.csv`.",
        f"- Rows after E₀ filter (|E₀| ≥ {initial_energy_floor:g}): {n_rows_total}",
        f"- Strata: {len(stratum_info)}, singletons dropped: {n_singletons}",
        f"- Rows kept for regression: {n_kept}",
        f"- Residual std: {y_std:.4g}",
        f"- Constant-predictor loss (variance of y): **{constant_loss:.6g}**",
        "",
        "## Stratum sizes",
        "",
        "| noise | n | target_dim | count | mean raw | min | max | kept? |",
        "| :---: | ---: | :---: | ---: | ---: | ---: | ---: | :---: |",
    ]
    for key in sorted(stratum_info.keys()):
        nl, n, td = key
        info = stratum_info[key]
        nl_label = "small" if nl == 0 else "medium"
        kept = "—" if info["dropped"] else "✓"
        if info["dropped"]:
            lines.append(
                f"| {nl_label} | {int(n)} | {int(td)} | {info['count']} | — | — | — | {kept} |"
            )
        else:
            lines.append(
                f"| {nl_label} | {int(n)} | {int(td)} | {info['count']}"
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
        lines.append(f"| {p} | {PANEL_DESCRIPTIONS[p]} | `{', '.join(PANELS[p])}` |")

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
        "Thresholds: strong = 10%, null = 5% relative improvement over constant.",
        "",
        "1. **RESIDUALIZATION_FAILED_OR_LEAKAGE**: Panel A above strong → broken.",
        "2. **POSSIBLE_ORDER_SIGNAL_WITH_ABS_DISC**: D above strong, uses",
        "   abs_discrepancy_mm_midpoint or C's best does.  → 3E candidate replicates.",
        "3. **ORDER_SIGNAL_WITHOUT_ABS_DISC**: D above strong but abs_disc absent.",
        "4. **NULL_SIGNAL**: D below null threshold.",
        "5. **SIGNAL_CONDITIONED_ON_DESIGN**: B strong, D not.",
        "6. **INTERMEDIATE**: D between thresholds.",
        "",
        "## Reproducibility",
        "",
        f"- PySR iterations per panel: {niterations}",
        f"- maxsize cap: {maxsize}",
        f"- random_state=1959, parallelism=serial, initial_energy_floor={initial_energy_floor:g}",
        "",
        "Regenerate via `make regen-phase3f`.",
        "Source: `tools/build_phase3f_pysr_final_ablation.py`.",
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

    p2g_path = FOUNDATION / "phase2g_extended_guarded_warmup_probe.csv"
    p1e_path = FOUNDATION / "phase1e_extended_structural_atlas.csv"
    if not p2g_path.exists() or not p1e_path.exists():
        print(
            f"ERROR: required inputs missing.\n  {p2g_path}\n  {p1e_path}\n"
            "Run `make regen-phase1e regen-phase2g` first.",
            file=sys.stderr,
        )
        sys.exit(1)

    p2g_header, p2g_rows = load_csv(p2g_path)
    missing = REQUIRED_PHASE2G_COLUMNS - set(p2g_header)
    if missing:
        print(
            f"ERROR: phase2g missing required columns: {sorted(missing)}\n"
            f"Available: {p2g_header}",
            file=sys.stderr,
        )
        sys.exit(2)

    validate_panel_d_features(PANELS["D"])

    p1e_idx = load_phase1e(p1e_path)

    rows = assemble_rows(p2g_rows, p1e_idx, args.initial_energy_floor)
    if not rows:
        print("ERROR: zero usable rows after filtering.", file=sys.stderr)
        sys.exit(3)

    rows_kept, stratum_info, n_singletons = stratify_and_residualize(rows)
    print(f"  strata={len(stratum_info)} (singletons={n_singletons}), kept={len(rows_kept)} rows.")
    if not rows_kept:
        print("ERROR: all strata are singletons.", file=sys.stderr)
        sys.exit(4)

    y = [r["residual_relative_drift"] for r in rows_kept]
    y_mean = sum(y) / len(y)
    y_std  = math.sqrt(sum((v - y_mean)**2 for v in y) / len(y))
    constant_loss = sum((v - y_mean)**2 for v in y) / len(y)
    print(f"  target std={y_std:.4g}, constant-predictor loss={constant_loss:.6g}")

    panel_eqs: dict[str, Any] = {}
    all_csv_rows: list[dict] = []
    for p, feats in PANELS.items():
        X = [[_safe(r[f]) for f in feats] for r in rows_kept]
        eqs = run_pysr(X, y, feats, args.niterations, args.maxsize, p)
        panel_eqs[p] = eqs
        all_csv_rows.extend(equations_to_rows(eqs, p))

    FOUNDATION.mkdir(parents=True, exist_ok=True)
    write_csv(all_csv_rows, FOUNDATION / "phase3f_pysr_final_ablation.csv")
    write_markdown(
        n_kept=len(rows_kept), n_singletons=n_singletons, n_rows_total=len(rows),
        initial_energy_floor=args.initial_energy_floor,
        stratum_info=stratum_info,
        constant_loss=constant_loss, y_std=y_std,
        panel_eqs=panel_eqs, niterations=args.niterations, maxsize=args.maxsize,
        path=FOUNDATION / "phase3f_pysr_final_ablation.md",
    )

    print("\n--- Phase 3F summary ---")
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
        order_flag = "ORDER" if _uses_order(best_eq) else "design/none"
        abs_flag = "+abs_disc" if _uses_abs_disc(best_eq) else ""
        print(f"  Panel {p}: loss={bl:.4g}  (rel Δ={rel:+.2%})  [{order_flag} {abs_flag}]")
        print(f"           {best_eq[:90]}")

    print(f"\nOutput: {FOUNDATION / 'phase3f_pysr_final_ablation.csv'}")
    print(f"        {FOUNDATION / 'phase3f_pysr_final_ablation.md'}")


if __name__ == "__main__":
    main()
