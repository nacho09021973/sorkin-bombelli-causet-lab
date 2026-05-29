#!/usr/bin/env python3
"""S4-KERR-K17E-NECESSARY-CAUSAL-FILTER-AUDIT-001.

Apply necessary / pre-shooting filters to the K17d selected pairs.

This is not K18, not a causal classifier, not reachability. A pair that
survives the filter is "not_excluded" — it has not been ruled out by the
necessary conditions encoded here. A pair that fails is "rejected_by_filter"
— this is not proof of spacelike separation.

Filters (all "necessary, not sufficient"):
  1. time_order_pass:                           delta_t_AB > 0
  2. schwarzschild_radial_time_bound_pass:      hard only for spin_a == 0
  3. angular_sector_admissibility_pass:         min_m |delta_phi + 2 pi m| <= PHI_SECTOR_MAX
  4. radial_potential_admissibility_pass:       exists b with min_r R(r;M,a,b,E=1) >= -R_TOL
  5. combined_not_excluded:                     conjunction of the above

Critical constraints respected:
  * Reuses radial_potential from K9 unmodified.
  * Reuses _principal_angle / _best_sector_abs / SECTORS / TIME_TOLERANCE_BAND
    / PHI_SECTOR_MAX from K17c unmodified.
  * Reuses B_GRID and W_TOL from K17 unmodified.
  * Reuses outgoing_radial_trip from run_kerr_minimal_benchmark unmodified, and
    ONLY for spin_a == 0 (Schwarzschild). For spin_a > 0 the radial-time filter
    is recorded as not_applied_kerr — no Kerr radial-time lower bound invented.
  * Does not modify any K17d artifact; reads CSV/JSON only.
  * Does not modify cones.py or the Schwarzschild benchmark code.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from explore.sorkin4_kerr_benchmark import run_kerr_minimal_benchmark as kerr  # noqa: E402
from explore.sorkin4_kerr_benchmark import (  # noqa: E402
    audit_kerr_k17_controlled_candidate_pair_sandbox_001 as k17,
)
from explore.sorkin4_kerr_benchmark import (  # noqa: E402
    audit_kerr_k17c_candidate_pair_selection_audit_001 as k17c,
)
from explore.sorkin4_kerr_benchmark.audit_kerr_k9_equatorial_full_rhs_preflight_001 import (  # noqa: E402
    R_MIN_TOL,
    radial_potential,
)

ARTIFACT_DIR = Path(__file__).resolve().parent
INPUT_CSV = ARTIFACT_DIR / "kerr_k17d_cloud_size_seed_scan_001.csv"
INPUT_JSON = ARTIFACT_DIR / "kerr_k17d_cloud_size_seed_scan_001.json"
OUT_PREFIX = "kerr_k17e_necessary_causal_filter_audit_001"

MASS = 1.0
ENERGY = 1.0
TIME_TOLERANCE_BAND = k17c.TIME_TOLERANCE_BAND  # 0.05
PHI_SECTOR_MAX = k17c.PHI_SECTOR_MAX  # 1.6
SECTORS = k17c.SECTORS  # (-2, -1, 0, 1, 2)
B_GRID = k17.B_GRID  # (-2,-1,-0.5,0,0.5,1,2) — broader than K17d PROBE_B_GRID
W_TOL = k17.W_TOL  # 1.0e-3
R_TOL = R_MIN_TOL  # 1.0e-10
N_R_SAMPLES = 32  # radial sampling resolution for R(r) scan

CAVEATS = [
    "not_excluded is not reachability.",
    "rejected_by_filter is not proof of spacelike separation.",
    "no causal_true/false relations decided.",
    "no production classifier introduced.",
    "no global Kerr causal claim introduced.",
    "no Level-B Hawking/Bekenstein claim introduced.",
    "Schwarzschild radial-time bound applied as hard filter only for spin_a == 0.",
    "For spin_a > 0 the radial-time filter is recorded as not_applied_kerr; no Kerr lower bound invented.",
    "The radial_potential filter is necessary, not sufficient.",
]


# ---------------------------------------------------------------------------
# Filter primitives — each returns enough metadata to populate the row
# ---------------------------------------------------------------------------


def _time_order_pass(dt: float) -> bool:
    return dt > 0.0


def _schwarzschild_radial_time_pass(
    *,
    spin: float,
    dt: float,
    r_a: float,
    r_b: float,
    mass: float,
    tol: float,
) -> tuple[bool, str, float | None]:
    """Hard Schwarzschild radial-time filter; only applied when spin == 0.

    Returns (survives, mode, bound). For spin != 0 the bound is None and
    survives is True (filter not applied — no Kerr lower bound invented).
    """
    if spin != 0.0:
        return True, "not_applied_kerr", None
    r1 = min(r_a, r_b)
    r2 = max(r_a, r_b)
    try:
        bound = kerr.outgoing_radial_trip(r1, r2, mass)
    except (ValueError, ZeroDivisionError):
        return True, "schwarzschild_bound_undefined", None
    return (dt >= bound - tol), "schwarzschild_hard", bound


def _angular_sector_pass(dphi_raw: float) -> tuple[bool, int, float]:
    dphi_mod = k17c._principal_angle(dphi_raw)
    sector_abs = k17c._best_sector_abs(dphi_mod)
    best_m = min(SECTORS, key=lambda m: abs(dphi_mod + 2.0 * math.pi * m))
    return (sector_abs <= PHI_SECTOR_MAX), best_m, sector_abs


def _radial_potential_pass(
    *,
    spin: float,
    r_a: float,
    r_b: float,
    mass: float,
    b_grid: tuple[float, ...],
    r_tol: float,
    n_r: int,
) -> tuple[bool, float | None, float | None]:
    """Necessary radial-potential filter.

    Pass if exists b in b_grid with min over sampled r in [min(r_a,r_b), max(r_a,r_b)]
    of R(r; M, a, b, E=1) >= -r_tol. Returns (survives, best_b, min_R_at_best_b).
    "Best b" = the b with the largest sampled min_R; not necessarily the b that
    a real geodesic would take.
    """
    r1 = min(r_a, r_b)
    r2 = max(r_a, r_b)
    if not math.isfinite(r1) or not math.isfinite(r2):
        return False, None, None
    if r2 <= r1:
        rs = [r1]
    else:
        step = (r2 - r1) / max(n_r - 1, 1)
        rs = [r1 + i * step for i in range(n_r)]

    best_b: float | None = None
    best_min_R: float = -math.inf
    for b in b_grid:
        try:
            min_R = min(radial_potential(r, mass, spin, b, energy=ENERGY) for r in rs)
        except (ValueError, ZeroDivisionError):
            continue
        if not math.isfinite(min_R):
            continue
        if min_R > best_min_R:
            best_min_R = min_R
            best_b = b
    if best_b is None:
        return False, None, None
    return (best_min_R >= -r_tol), best_b, best_min_R


# ---------------------------------------------------------------------------
# CSV ingestion + per-row processing
# ---------------------------------------------------------------------------


def _to_float(x: Any) -> float:
    if x is None or x == "":
        return float("nan")
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def _to_bool(x: Any) -> bool:
    if isinstance(x, bool):
        return x
    return str(x).strip().lower() in {"true", "1", "yes"}


def _load_k17d_rows() -> tuple[list[dict[str, str]], int, int]:
    """Return (rows_selected, n_total_rows, n_selected_rows).

    Sanity-check: K17d only emits rows for selected_* labels, but we filter
    explicitly to be robust against schema drift.
    """
    if not INPUT_CSV.exists():
        raise FileNotFoundError(f"Missing K17d CSV: {INPUT_CSV}")
    rows_all: list[dict[str, str]] = []
    rows_selected: list[dict[str, str]] = []
    with INPUT_CSV.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            rows_all.append(row)
            label = row.get("selection_label", "")
            if label.startswith("selected_"):
                rows_selected.append(row)
    return rows_selected, len(rows_all), len(rows_selected)


def _process_row(row: dict[str, str]) -> dict[str, Any]:
    spin = _to_float(row["spin_a"])
    n_val = int(row["N"])
    seed = int(row["seed"])
    r_a = _to_float(row["r_A"])
    r_b = _to_float(row["r_B"])
    dt = _to_float(row["delta_t_AB"])
    dphi_raw = _to_float(row["delta_phi_AB"])
    best_residual_k17d = _to_float(row.get("best_residual", ""))
    probe_succeeded = _to_bool(row.get("probe_succeeded", "False"))

    f1 = _time_order_pass(dt)
    f2_pass, f2_mode, f2_bound = _schwarzschild_radial_time_pass(
        spin=spin, dt=dt, r_a=r_a, r_b=r_b, mass=MASS, tol=TIME_TOLERANCE_BAND
    )
    f3_pass, best_m, sector_abs = _angular_sector_pass(dphi_raw)
    f4_pass, best_b, min_R = _radial_potential_pass(
        spin=spin,
        r_a=r_a,
        r_b=r_b,
        mass=MASS,
        b_grid=B_GRID,
        r_tol=R_TOL,
        n_r=N_R_SAMPLES,
    )

    combined = bool(f1 and f2_pass and f3_pass and f4_pass)

    return {
        "case_id": row.get("case_id"),
        "N": n_val,
        "seed": seed,
        "spin_a": spin,
        "event_A_index": int(row["event_A_index"]),
        "event_B_index": int(row["event_B_index"]),
        "t_A": _to_float(row["t_A"]),
        "r_A": r_a,
        "phi_A": _to_float(row["phi_A"]),
        "t_B": _to_float(row["t_B"]),
        "r_B": r_b,
        "phi_B": _to_float(row["phi_B"]),
        "delta_t_AB": dt,
        "delta_r_AB": _to_float(row["delta_r_AB"]),
        "delta_phi_AB": dphi_raw,
        "selection_label_k17d": row.get("selection_label"),
        "probe_succeeded_k17d": probe_succeeded,
        "best_residual_k17d": best_residual_k17d,
        "best_b_k17d": row.get("best_b"),
        "best_lambda_k17d": row.get("best_lambda"),
        "best_direction_k17d": row.get("best_direction"),
        "best_sector_m_k17d": row.get("best_sector_m"),
        # K17E filter outputs
        "survives_time_order_filter": f1,
        "schwarzschild_radial_time_bound_pass": f2_pass,
        "radial_time_filter_mode": f2_mode,
        "schwarzschild_radial_time_bound": f2_bound,
        "survives_sector_filter": f3_pass,
        "best_sector_m_k17e": best_m,
        "best_sector_abs_k17e": sector_abs,
        "survives_radial_potential_filter": f4_pass,
        "best_b_k17e": best_b,
        "min_R_at_best_b_k17e": min_R,
        "combined_not_excluded": combined,
        # Framing flags (selection-only — not causal verdicts)
        "not_excluded_is_not_reachability": True,
        "no_causal_claim_introduced": True,
    }


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------


def _finite_residuals(
    rows: list[dict[str, Any]], *, only_combined: bool
) -> list[float]:
    out: list[float] = []
    for r in rows:
        if only_combined and not r["combined_not_excluded"]:
            continue
        v = r["best_residual_k17d"]
        if isinstance(v, float) and math.isfinite(v):
            out.append(v)
    return out


def _dist(vals: list[float]) -> dict[str, Any]:
    if not vals:
        return {"n": 0, "min": None, "median": None, "max": None}
    s = sorted(vals)
    n = len(s)
    median = s[n // 2] if n % 2 == 1 else 0.5 * (s[n // 2 - 1] + s[n // 2])
    return {"n": n, "min": s[0], "median": median, "max": s[-1]}


def _summarize(
    rows: list[dict[str, Any]], n_csv_total: int, n_csv_selected: int
) -> dict[str, Any]:
    survivors = {
        "time_order": sum(1 for r in rows if r["survives_time_order_filter"]),
        "schwarzschild_radial_time_bound": sum(
            1 for r in rows if r["schwarzschild_radial_time_bound_pass"]
        ),
        "angular_sector": sum(1 for r in rows if r["survives_sector_filter"]),
        "radial_potential": sum(
            1 for r in rows if r["survives_radial_potential_filter"]
        ),
        "combined_not_excluded": sum(1 for r in rows if r["combined_not_excluded"]),
    }
    # Independent per-filter rejection counts. A pair can be counted under
    # multiple reasons; the sum can exceed (n_total - combined_survivors).
    rejection_reasons = {
        "rejected_time_order": sum(
            1 for r in rows if not r["survives_time_order_filter"]
        ),
        "rejected_schwarzschild_radial_time": sum(
            1
            for r in rows
            if r["radial_time_filter_mode"] == "schwarzschild_hard"
            and not r["schwarzschild_radial_time_bound_pass"]
        ),
        "rejected_angular_sector": sum(
            1 for r in rows if not r["survives_sector_filter"]
        ),
        "rejected_radial_potential": sum(
            1 for r in rows if not r["survives_radial_potential_filter"]
        ),
    }
    radial_time_mode_histogram: dict[str, int] = {}
    for r in rows:
        m = r["radial_time_filter_mode"]
        radial_time_mode_histogram[m] = radial_time_mode_histogram.get(m, 0) + 1

    before_vals = _finite_residuals(rows, only_combined=False)
    after_vals = _finite_residuals(rows, only_combined=True)

    near_hits_before = sum(1 for v in before_vals if W_TOL < v <= 10.0 * W_TOL)
    near_hits_after = sum(1 for v in after_vals if W_TOL < v <= 10.0 * W_TOL)
    hits_before = sum(1 for v in before_vals if v <= W_TOL)
    hits_after = sum(1 for v in after_vals if v <= W_TOL)

    per_cell: dict[tuple[int, int, float], dict[str, int]] = {}
    for r in rows:
        key = (r["N"], r["seed"], r["spin_a"])
        rec = per_cell.setdefault(key, {"n_in": 0, "n_out": 0})
        rec["n_in"] += 1
        if r["combined_not_excluded"]:
            rec["n_out"] += 1

    per_cell_rows = [
        {
            "N": k[0],
            "seed": k[1],
            "spin_a": k[2],
            "n_in": v["n_in"],
            "n_out": v["n_out"],
            "survival_rate": (v["n_out"] / v["n_in"]) if v["n_in"] else 0.0,
        }
        for k, v in sorted(per_cell.items())
    ]

    return {
        "total_k17d_csv_rows": n_csv_total,
        "total_k17d_pairs_read": n_csv_selected,
        "survivors": survivors,
        "rejection_reasons": rejection_reasons,
        "radial_time_mode_histogram": radial_time_mode_histogram,
        "residual_distribution_before": _dist(before_vals),
        "residual_distribution_after": _dist(after_vals),
        "best_residual_before": (min(before_vals) if before_vals else None),
        "best_residual_after": (min(after_vals) if after_vals else None),
        "hits_le_w_tol_before": hits_before,
        "hits_le_w_tol_after": hits_after,
        "near_hits_before": near_hits_before,
        "near_hits_after": near_hits_after,
        "per_cell_survival": per_cell_rows,
        "constants": {
            "MASS": MASS,
            "ENERGY": ENERGY,
            "W_TOL": W_TOL,
            "R_TOL": R_TOL,
            "TIME_TOLERANCE_BAND": TIME_TOLERANCE_BAND,
            "PHI_SECTOR_MAX": PHI_SECTOR_MAX,
            "B_GRID": list(B_GRID),
            "SECTORS": list(SECTORS),
            "N_R_SAMPLES": N_R_SAMPLES,
        },
        "caveats": CAVEATS,
    }


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

CSV_FIELDS = (
    "case_id",
    "N",
    "seed",
    "spin_a",
    "event_A_index",
    "event_B_index",
    "t_A",
    "r_A",
    "phi_A",
    "t_B",
    "r_B",
    "phi_B",
    "delta_t_AB",
    "delta_r_AB",
    "delta_phi_AB",
    "selection_label_k17d",
    "probe_succeeded_k17d",
    "best_residual_k17d",
    "best_b_k17d",
    "best_lambda_k17d",
    "best_direction_k17d",
    "best_sector_m_k17d",
    "survives_time_order_filter",
    "schwarzschild_radial_time_bound_pass",
    "radial_time_filter_mode",
    "schwarzschild_radial_time_bound",
    "survives_sector_filter",
    "best_sector_m_k17e",
    "best_sector_abs_k17e",
    "survives_radial_potential_filter",
    "best_b_k17e",
    "min_R_at_best_b_k17e",
    "combined_not_excluded",
    "not_excluded_is_not_reachability",
    "no_causal_claim_introduced",
)


def write_csv(rows: list[dict[str, Any]], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in CSV_FIELDS})


def _json_default(o: Any) -> Any:
    if isinstance(o, float) and not math.isfinite(o):
        return None
    return str(o)


def write_json(rows: list[dict[str, Any]], summary: dict[str, Any], path: Path) -> None:
    payload = {
        "schema": "S4-KERR-K17E-NECESSARY-CAUSAL-FILTER-AUDIT-001",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "input_csv": INPUT_CSV.name,
        "input_json": INPUT_JSON.name,
        "summary": summary,
        "rows": rows,
    }
    path.write_text(
        json.dumps(payload, indent=2, default=_json_default), encoding="utf-8"
    )


def write_md(summary: dict[str, Any], path: Path) -> None:
    lines: list[str] = []
    lines.append("# S4-KERR-K17E necessary causal filter audit")
    lines.append("")
    lines.append("Necessary / pre-shooting filters applied to K17d selected pairs.")
    lines.append(
        "not_excluded is not reachability. rejected_by_filter is not proof of spacelike separation."
    )
    lines.append("")
    lines.append(f"K17d CSV rows total: {summary['total_k17d_csv_rows']}")
    lines.append(f"K17d selected pairs read: {summary['total_k17d_pairs_read']}")
    lines.append("")
    lines.append("## Per-filter survivor counts")
    lines.append("")
    for key, val in summary["survivors"].items():
        lines.append(f"- {key}: {val}")
    lines.append("")
    lines.append(
        "## Independent rejection counts (a pair may be counted in multiple reasons)"
    )
    lines.append("")
    for key, val in summary["rejection_reasons"].items():
        lines.append(f"- {key}: {val}")
    lines.append("")
    lines.append("## radial_time_filter_mode histogram")
    lines.append("")
    for k, v in sorted(summary["radial_time_mode_histogram"].items()):
        lines.append(f"- {k}: {v}")
    lines.append("")
    lines.append("## Residual distribution before / after filtering")
    lines.append("")
    rb = summary["residual_distribution_before"]
    ra = summary["residual_distribution_after"]
    lines.append(
        f"- before: n={rb['n']}, min={rb['min']}, median={rb['median']}, max={rb['max']}"
    )
    lines.append(
        f"- after:  n={ra['n']}, min={ra['min']}, median={ra['median']}, max={ra['max']}"
    )
    lines.append(f"- best_residual_before = {summary['best_residual_before']}")
    lines.append(f"- best_residual_after  = {summary['best_residual_after']}")
    lines.append(
        f"- hits_le_W_TOL before / after = {summary['hits_le_w_tol_before']} / {summary['hits_le_w_tol_after']}"
    )
    lines.append(
        f"- near_hits before / after = {summary['near_hits_before']} / {summary['near_hits_after']}"
    )
    lines.append("")
    lines.append("## Per (N, seed, spin) survival")
    lines.append("")
    lines.append("| N | seed | spin | n_in | n_out | survival_rate |")
    lines.append("|---|---|---|---|---|---|")
    for c in summary["per_cell_survival"]:
        lines.append(
            f"| {c['N']} | {c['seed']} | {c['spin_a']:.2f} | "
            f"{c['n_in']} | {c['n_out']} | {c['survival_rate']:.3f} |"
        )
    lines.append("")
    lines.append("## Constants used")
    lines.append("")
    for k, v in summary["constants"].items():
        lines.append(f"- {k} = {v}")
    lines.append("")
    lines.append("## Defaults documented")
    lines.append("")
    lines.append(
        "- B_GRID is the 7-value K17.B_GRID, not the 5-value K17d.PROBE_B_GRID — chosen so the necessary radial-potential filter does not over-reject."
    )
    lines.append(
        "- TIME_TOLERANCE_BAND reuses K17c (0.05) — same band that entered the K17d selection."
    )
    lines.append(
        "- R_TOL reuses K9.R_MIN_TOL (1e-10) — the same tolerance K17 uses internally for R >= 0."
    )
    lines.append("- N_R_SAMPLES = 32 is a sampling resolution decision, not physics.")
    lines.append(
        "- For spin_a > 0, radial_time_filter_mode = not_applied_kerr. No Kerr lower bound invented."
    )
    lines.append("")
    lines.append("## Caveats")
    lines.append("")
    for c in summary["caveats"]:
        lines.append(f"- {c}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_png(rows: list[dict[str, Any]], summary: dict[str, Any], path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    before = _finite_residuals(rows, only_combined=False)
    after = _finite_residuals(rows, only_combined=True)

    ax = axes[0, 0]
    if before:
        ax.hist(before, bins=30, alpha=0.55, label=f"K17d selected (n={len(before)})")
    if after:
        ax.hist(after, bins=30, alpha=0.75, label=f"K17E not_excluded (n={len(after)})")
    ax.axvline(W_TOL, linestyle="--", linewidth=0.8, label="W_TOL")
    ax.axvline(10.0 * W_TOL, linestyle=":", linewidth=0.8, label="10 W_TOL")
    ax.set_xlabel("best_residual (from K17d)")
    ax.set_ylabel("count")
    ax.set_xscale("log") if before and min(before) > 0 else None
    ax.legend(fontsize=8)
    ax.set_title("Residual distribution before / after K17E filter")

    ax = axes[0, 1]
    surv = summary["survivors"]
    ax.bar(list(surv.keys()), list(surv.values()))
    ax.set_title("Per-filter survivor counts")
    ax.tick_params(axis="x", rotation=25, labelsize=8)
    for label in ax.get_xticklabels():
        label.set_horizontalalignment("right")

    ax = axes[1, 0]
    cells = summary["per_cell_survival"]
    labels = [f"N={c['N']},s{c['seed']},a{c['spin_a']:.2f}" for c in cells]
    rates = [c["survival_rate"] for c in cells]
    y = list(range(len(labels)))
    ax.barh(y, rates)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=6)
    ax.set_xlim(0.0, 1.0)
    ax.set_xlabel("survival_rate")
    ax.set_title("Per (N, seed, spin) survival rate")

    ax = axes[1, 1]
    rej = summary["rejection_reasons"]
    nonzero = [(k, v) for k, v in rej.items() if v > 0]
    if nonzero:
        ax.pie(
            [v for _, v in nonzero], labels=[k for k, _ in nonzero], autopct="%1.1f%%"
        )
        ax.set_title("Independent rejection-reason counts")
    else:
        ax.text(0.5, 0.5, "no rejections", ha="center", va="center")
        ax.set_title("Independent rejection-reason counts")
        ax.set_axis_off()

    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Entrypoint
# ---------------------------------------------------------------------------


def main() -> None:
    rows_selected, n_csv_total, n_csv_selected = _load_k17d_rows()
    processed = [_process_row(r) for r in rows_selected]
    summary = _summarize(processed, n_csv_total, n_csv_selected)

    csv_path = ARTIFACT_DIR / f"{OUT_PREFIX}.csv"
    json_path = ARTIFACT_DIR / f"{OUT_PREFIX}.json"
    md_path = ARTIFACT_DIR / f"{OUT_PREFIX}.md"
    png_path = ARTIFACT_DIR / f"{OUT_PREFIX}.png"

    write_csv(processed, csv_path)
    write_json(processed, summary, json_path)
    write_md(summary, md_path)
    write_png(processed, summary, png_path)

    print(
        f"K17E: read {summary['total_k17d_pairs_read']} K17d selected pairs; "
        f"not_excluded = {summary['survivors']['combined_not_excluded']}; "
        f"best_residual_before = {summary['best_residual_before']}; "
        f"best_residual_after = {summary['best_residual_after']}; "
        f"near_hits_before / after = "
        f"{summary['near_hits_before']} / {summary['near_hits_after']}"
    )


if __name__ == "__main__":
    main()
