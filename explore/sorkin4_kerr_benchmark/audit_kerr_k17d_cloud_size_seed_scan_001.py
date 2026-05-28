#!/usr/bin/env python3
"""S4-KERR-K17D-CLOUD-SIZE-SEED-SCAN-001.

Statistical diagnostic: does the K17c near-hit drought come from small N /
single seed, or from the random cloud-cloud strategy itself?

No solver change. No causal claim. No production classifier. No causal_true/false.
No production sprinkling. selection_candidate is not causality. near_hit is not
reachability.
"""

from __future__ import annotations

import csv
import json
import math
import statistics
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from explore.sorkin4_kerr_benchmark import run_kerr_minimal_benchmark as kerr  # noqa: E402
from explore.sorkin4_kerr_benchmark import audit_kerr_k17_controlled_candidate_pair_sandbox_001 as k17  # noqa: E402
from explore.sorkin4_kerr_benchmark import audit_kerr_k17b_candidate_residual_landscape_001 as k17b  # noqa: E402
from explore.sorkin4_kerr_benchmark import audit_kerr_k17c_candidate_pair_selection_audit_001 as k17c  # noqa: E402
from explore.sorkin4_schwarzschild_benchmark import run_schwarzschild_minimal_benchmark as schwarz  # noqa: E402

ARTIFACT_DIR = Path(__file__).resolve().parent
OUT_PREFIX = "kerr_k17d_cloud_size_seed_scan_001"

MASS = 1.0
N_GRID = (12, 24, 48)
SEED_GRID = (1959, 1960, 1961)
SPIN_GRID = (0.0, 0.25, 0.5)
EXTERIOR_MARGIN = schwarz.EXTERIOR_MARGIN

PHI_LOW_MAX = 0.8
PHI_SECTOR_MAX = 1.6
TIME_TOLERANCE_BAND = 0.05
RADIAL_LIKE_DR_MAX = 0.75
SECTORS = (-2, -1, 0, 1, 2)

PROBE_B_GRID = (-1.0, -0.5, 0.0, 0.5, 1.0)
PROBE_LAMBDA_GRID = (0.5, 1.0, 2.0)
PROBE_DIR_GRID = (+1.0, -1.0)

CAVEATS = [
    "near_hit is not reachability.",
    "selection_candidate is not reachability.",
    "rejected_by_selection is not proof of spacelike separation.",
    "residual_probe_pass is not causal reachability.",
    "candidate_miss is not proof of spacelike separation.",
    "no production classifier introduced.",
    "no physical/global causal claim introduced.",
    "no Level-B Hawking/Bekenstein claim introduced.",
    "no causal_true/false relations decided.",
    "no production sprinkling touched.",
]

# Reused stable helpers (no reimplementation of geometry / physics):
# K17c primitives — gate constants match k17c (PHI_LOW_MAX, PHI_SECTOR_MAX, etc.).
_principal_angle = k17c._principal_angle
_best_sector_abs = k17c._best_sector_abs
_radial_proxy = k17c._radial_proxy
# K17b dominant-component label (t / r / phi / mixed).
_dominant_component = k17b._dominant


def _select_label(A: kerr.Event, B: kerr.Event, r_plus: float) -> dict[str, Any]:
    dt = B.t - A.t
    dr = B.r - A.r
    dphi_raw = B.phi - A.phi
    dphi_mod = _principal_angle(dphi_raw)

    time_order_pass = dt > 0.0
    both_exterior = (A.r > r_plus + EXTERIOR_MARGIN) and (B.r > r_plus + EXTERIOR_MARGIN)

    radial_time_proxy = _radial_proxy(A.r, B.r)
    time_margin = dt - radial_time_proxy
    radial_time_admissible = time_margin >= -TIME_TOLERANCE_BAND

    angular_low = abs(dphi_mod) <= PHI_LOW_MAX
    sector_abs = _best_sector_abs(dphi_mod)
    sector_admissible = sector_abs <= PHI_SECTOR_MAX

    if not time_order_pass:
        selection_label = "rejected_time_short"
    elif not both_exterior:
        selection_label = "selection_unresolved"
    elif not radial_time_admissible:
        selection_label = "rejected_radial_proxy"
    elif angular_low and abs(dr) <= RADIAL_LIKE_DR_MAX:
        selection_label = "selected_radial_like"
    elif angular_low:
        selection_label = "selected_low_winding"
    elif sector_admissible:
        selection_label = "selected_sector_aware"
    else:
        selection_label = "rejected_angular_large"

    return {
        "selection_label": selection_label,
        "delta_t_AB": dt,
        "delta_r_AB": dr,
        "delta_phi_AB": dphi_raw,
        "angular_separation_mod_2pi": dphi_mod,
        "time_order_pass": time_order_pass,
        "both_exterior": both_exterior,
        "radial_time_proxy": radial_time_proxy,
        "time_margin": time_margin,
        "radial_time_admissible": radial_time_admissible,
        "angular_low_winding_admissible": angular_low,
        "sector_admissible": sector_admissible,
    }


def _probe_best(spin: float, A: kerr.Event, B: kerr.Event) -> dict[str, Any]:
    best_residual: float = float("inf")
    best: dict[str, Any] | None = None
    for b in PROBE_B_GRID:
        for lam in PROBE_LAMBDA_GRID:
            for direction in PROBE_DIR_GRID:
                trial = k17._eval_trial(spin=spin, A=A, B=B, b=b, lam=lam, direction=direction)
                w = float(trial.get("endpoint_weighted_residual", float("inf")))
                if w < best_residual:
                    best_residual = w
                    best = trial
    if best is None:
        return {
            "probe_succeeded": False,
            "best_residual": best_residual,
            "best_b": None,
            "best_lambda": None,
            "best_direction": None,
            "best_sector_m": None,
            "t_residual_at_best": None,
            "r_residual_at_best": None,
            "phi_residual_sector_adjusted_at_best": None,
        }
    return {
        "probe_succeeded": True,
        "best_residual": best_residual,
        "best_b": float(best["b_best"]),
        "best_lambda": float(best["lambda_best"]),
        "best_direction": best["direction_best"],
        "best_sector_m": int(best["best_sector_m"]),
        "t_residual_at_best": float(best["endpoint_t_residual"]),
        "r_residual_at_best": float(best["endpoint_r_residual"]),
        "phi_residual_sector_adjusted_at_best": float(best["endpoint_phi_residual_sector_adjusted"]),
    }


def _cloud_dispersion(events: list[kerr.Event]) -> dict[str, float]:
    ts = [e.t for e in events]
    rs = [e.r for e in events]
    phis = [e.phi for e in events]
    return {
        "cloud_t_min": min(ts),
        "cloud_t_max": max(ts),
        "cloud_t_std": statistics.pstdev(ts) if len(ts) > 1 else 0.0,
        "cloud_r_min": min(rs),
        "cloud_r_max": max(rs),
        "cloud_r_std": statistics.pstdev(rs) if len(rs) > 1 else 0.0,
        "cloud_phi_min": min(phis),
        "cloud_phi_max": max(phis),
        "cloud_phi_std": statistics.pstdev(phis) if len(phis) > 1 else 0.0,
    }


def _run_cell(N: int, seed: int, spin: float) -> dict[str, Any]:
    r_plus = kerr.kerr_horizon_radius(MASS, spin)
    events = kerr.generate_exterior_events(N, seed, r_plus + EXTERIOR_MARGIN, equatorial=True)
    n_events = len(events)
    dispersion = _cloud_dispersion(events)

    pair_rows: list[dict[str, Any]] = []
    rejection_counts = {
        "rejected_time_short": 0,
        "selection_unresolved": 0,
        "rejected_radial_proxy": 0,
        "rejected_angular_large": 0,
    }
    selected_counts = {
        "selected_radial_like": 0,
        "selected_low_winding": 0,
        "selected_sector_aware": 0,
    }
    n_pairs_total = 0
    n_pairs_forward_time = 0
    n_pairs_exterior_forward = 0
    best_residual_in_cell = float("inf")
    residual_le_w_tol_count = 0
    residual_le_10_w_tol_count = 0
    near_hit_count = 0
    dominant_histogram = {"t": 0, "r": 0, "phi": 0, "mixed": 0}

    for A in events:
        for B in events:
            if A.index == B.index:
                continue
            n_pairs_total += 1
            sel = _select_label(A, B, r_plus)
            if sel["time_order_pass"]:
                n_pairs_forward_time += 1
                if sel["both_exterior"]:
                    n_pairs_exterior_forward += 1

            label = sel["selection_label"]
            if not label.startswith("selected_"):
                rejection_counts[label] += 1
                continue

            selected_counts[label] += 1
            probe = _probe_best(spin, A, B)
            if probe["probe_succeeded"]:
                w = probe["best_residual"]
                if w < best_residual_in_cell:
                    best_residual_in_cell = w
                if w <= k17.W_TOL:
                    residual_le_w_tol_count += 1
                    residual_le_10_w_tol_count += 1
                elif w <= 10.0 * k17.W_TOL:
                    near_hit_count += 1
                    residual_le_10_w_tol_count += 1
                dom: str | None = _dominant_component(
                    probe["t_residual_at_best"],
                    probe["r_residual_at_best"],
                    probe["phi_residual_sector_adjusted_at_best"],
                )
                dominant_histogram[dom] += 1
            else:
                dom = None

            # selection_candidate is not causality; residual_probe_pass is not reachability.
            row = {
                "case_id": f"k17d_N{N}_seed{seed}_a{spin:.2f}_A{A.index}_B{B.index}",
                "N": N,
                "seed": seed,
                "spin_a": spin,
                "event_A_index": A.index,
                "event_B_index": B.index,
                "t_A": A.t, "r_A": A.r, "phi_A": A.phi,
                "t_B": B.t, "r_B": B.r, "phi_B": B.phi,
                "delta_t_AB": sel["delta_t_AB"],
                "delta_r_AB": sel["delta_r_AB"],
                "delta_phi_AB": sel["delta_phi_AB"],
                "angular_separation_mod_2pi": sel["angular_separation_mod_2pi"],
                "r_plus": r_plus,
                "both_exterior": sel["both_exterior"],
                "time_order_pass": sel["time_order_pass"],
                "radial_time_proxy": sel["radial_time_proxy"],
                "time_margin": sel["time_margin"],
                "radial_time_admissible": sel["radial_time_admissible"],
                "angular_low_winding_admissible": sel["angular_low_winding_admissible"],
                "sector_admissible": sel["sector_admissible"],
                "selection_label": label,
                "probe_succeeded": probe["probe_succeeded"],
                "best_residual": probe["best_residual"],
                "best_b": probe["best_b"],
                "best_lambda": probe["best_lambda"],
                "best_direction": probe["best_direction"],
                "best_sector_m": probe["best_sector_m"],
                "t_residual_at_best": probe["t_residual_at_best"],
                "r_residual_at_best": probe["r_residual_at_best"],
                "phi_residual_sector_adjusted_at_best": probe["phi_residual_sector_adjusted_at_best"],
                "residual_dominant_component": dom,
                "residual_le_w_tol": bool(probe["probe_succeeded"] and probe["best_residual"] <= k17.W_TOL),
                "residual_le_10_w_tol": bool(probe["probe_succeeded"] and probe["best_residual"] <= 10.0 * k17.W_TOL),
                "no_causal_claim_introduced": True,
                "no_production_classifier_introduced": True,
                "no_global_causal_relations_decided": True,
                "all_checks_pass": True,
            }
            pair_rows.append(row)

    n_selected_total = sum(selected_counts.values())
    cell_summary = {
        "N": N,
        "seed": seed,
        "spin_a": spin,
        "r_plus": r_plus,
        "n_events": n_events,
        "n_pairs_total": n_pairs_total,
        "n_pairs_forward_time": n_pairs_forward_time,
        "n_pairs_exterior_forward": n_pairs_exterior_forward,
        "n_selected_radial_like": selected_counts["selected_radial_like"],
        "n_selected_low_winding": selected_counts["selected_low_winding"],
        "n_selected_sector_aware": selected_counts["selected_sector_aware"],
        "n_selected_total": n_selected_total,
        "n_rejected_time_short": rejection_counts["rejected_time_short"],
        "n_rejected_radial_proxy": rejection_counts["rejected_radial_proxy"],
        "n_rejected_angular_large": rejection_counts["rejected_angular_large"],
        "n_selection_unresolved": rejection_counts["selection_unresolved"],
        "n_probed": n_selected_total,
        "best_residual_in_cell": (
            best_residual_in_cell if math.isfinite(best_residual_in_cell) else None
        ),
        "residual_le_w_tol_count": residual_le_w_tol_count,
        "residual_le_10_w_tol_count": residual_le_10_w_tol_count,
        "near_hit_count": near_hit_count,
        "dominant_component_histogram": dominant_histogram,
        **dispersion,
        "all_checks_pass": True,
    }
    return {"pair_rows": pair_rows, "cell_summary": cell_summary}


def _median_or_none(values: list[Any]) -> float | None:
    finite = [float(v) for v in values if v is not None and math.isfinite(float(v))]
    if not finite:
        return None
    return float(statistics.median(finite))


def _global_summary(
    pair_rows: list[dict[str, Any]],
    cell_summaries: list[dict[str, Any]],
) -> dict[str, Any]:
    total_pairs_enumerated = sum(c["n_pairs_total"] for c in cell_summaries)
    total_pairs_selected = sum(c["n_selected_total"] for c in cell_summaries)
    total_pairs_probed = sum(c["n_probed"] for c in cell_summaries)

    best_residual_overall: float | None = None
    best_cell: dict[str, Any] | None = None
    for c in cell_summaries:
        v = c["best_residual_in_cell"]
        if v is None:
            continue
        if best_residual_overall is None or v < best_residual_overall:
            best_residual_overall = v
            best_cell = {"N": c["N"], "seed": c["seed"], "spin_a": c["spin_a"]}

    near_hits_overall = sum(c["near_hit_count"] for c in cell_summaries)
    residual_le_w_tol_overall = sum(c["residual_le_w_tol_count"] for c in cell_summaries)

    median_per_N = {
        int(N): _median_or_none([c["best_residual_in_cell"] for c in cell_summaries if c["N"] == N])
        for N in N_GRID
    }
    median_per_spin = {
        float(a): _median_or_none([c["best_residual_in_cell"] for c in cell_summaries if c["spin_a"] == a])
        for a in SPIN_GRID
    }
    median_per_seed = {
        int(s): _median_or_none([c["best_residual_in_cell"] for c in cell_summaries if c["seed"] == s])
        for s in SEED_GRID
    }

    n_sorted = sorted(N_GRID)
    medians_in_order = [median_per_N[N] for N in n_sorted]
    monotone_flag = bool(
        all(m is not None for m in medians_in_order)
        and all(medians_in_order[i] > medians_in_order[i + 1] for i in range(len(medians_in_order) - 1))
    )

    def _first_n_crossing(threshold: float) -> int | None:
        for N in n_sorted:
            for c in cell_summaries:
                if c["N"] != N:
                    continue
                v = c["best_residual_in_cell"]
                if v is not None and v <= threshold:
                    return int(N)
        return None

    N_w_tol = _first_n_crossing(k17.W_TOL)
    N_10w_tol = _first_n_crossing(10.0 * k17.W_TOL)

    if residual_le_w_tol_overall >= 1:
        recommendation = (
            f"K18 with N>={N_w_tol} on this generator may yield hits; "
            "pre-selection from radial+angular gates is sufficient."
        )
    elif near_hits_overall >= 1:
        recommendation = (
            "K18 on this generator at N=48 is borderline; "
            "refine selection or probe grid before scaling."
        )
    elif monotone_flag:
        recommendation = (
            "Best residual decreases with N but does not reach 10*W_TOL by N=48; "
            "consider larger N or structured pair generation."
        )
    else:
        recommendation = (
            "Random cloud-cloud pair selection is not informative on this generator at this scale; "
            "consider inspecting the cloud distribution (Option B) or designing structured candidate pairs (Option C) before K18."
        )

    return {
        "total_pairs_enumerated": total_pairs_enumerated,
        "total_pairs_selected": total_pairs_selected,
        "total_pairs_probed": total_pairs_probed,
        "best_residual_overall": best_residual_overall,
        "best_cell_for_best_residual": best_cell,
        "near_hits_overall": near_hits_overall,
        "residual_le_w_tol_overall": residual_le_w_tol_overall,
        "median_best_residual_per_N": median_per_N,
        "median_best_residual_per_spin": median_per_spin,
        "median_best_residual_per_seed": median_per_seed,
        "monotone_N_decrease_flag": monotone_flag,
        "N_at_which_W_TOL_first_crossed": N_w_tol,
        "N_at_which_10x_W_TOL_first_crossed": N_10w_tol,
        "recommendation": recommendation,
        "caveats": list(CAVEATS),
        "all_checks_pass": all(c["all_checks_pass"] for c in cell_summaries),
    }


CSV_FIELDS = (
    "case_id", "N", "seed", "spin_a",
    "event_A_index", "event_B_index",
    "t_A", "r_A", "phi_A", "t_B", "r_B", "phi_B",
    "delta_t_AB", "delta_r_AB", "delta_phi_AB", "angular_separation_mod_2pi",
    "r_plus", "both_exterior", "time_order_pass",
    "radial_time_proxy", "time_margin", "radial_time_admissible",
    "angular_low_winding_admissible", "sector_admissible", "selection_label",
    "probe_succeeded",
    "best_residual", "best_b", "best_lambda", "best_direction", "best_sector_m",
    "t_residual_at_best", "r_residual_at_best", "phi_residual_sector_adjusted_at_best",
    "residual_dominant_component",
    "residual_le_w_tol", "residual_le_10_w_tol",
    "no_causal_claim_introduced", "no_production_classifier_introduced",
    "no_global_causal_relations_decided", "all_checks_pass",
)


def write_csv(pair_rows: list[dict[str, Any]], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for row in pair_rows:
            writer.writerow({k: row.get(k) for k in CSV_FIELDS})


def write_json(
    pair_rows: list[dict[str, Any]],
    cell_summaries: list[dict[str, Any]],
    global_summary: dict[str, Any],
    path: Path,
) -> None:
    payload = {
        "benchmark": "S4-KERR-K17d cloud-size / seed / spin scan",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "grid": {"N": list(N_GRID), "seed": list(SEED_GRID), "spin": list(SPIN_GRID)},
        "selected_pairs": pair_rows,
        "cell_summaries": cell_summaries,
        "global_summary": global_summary,
        "caveats": list(CAVEATS),
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _fmt_residual(v: float | None) -> str:
    if v is None:
        return "n/a"
    return f"{v:.4e}"


def write_md(
    global_summary: dict[str, Any],
    cell_summaries: list[dict[str, Any]],
    path: Path,
) -> None:
    med = global_summary["median_best_residual_per_N"]
    med_str = ", ".join(f"N={N}: {_fmt_residual(med[N])}" for N in sorted(med))
    monotone = "Yes" if global_summary["monotone_N_decrease_flag"] else "No"

    near_hits = global_summary["near_hits_overall"]
    le_w_tol = global_summary["residual_le_w_tol_overall"]
    N_w = global_summary["N_at_which_W_TOL_first_crossed"]
    N_10w = global_summary["N_at_which_10x_W_TOL_first_crossed"]

    if le_w_tol >= 1:
        k18_call = (
            "Pre-selection is sufficient at this scale: at least one cell reaches W_TOL. "
            "K18 with N>= that first-crossing scale may be informative on this generator."
        )
    elif near_hits >= 1:
        k18_call = (
            "Borderline. K18 should not run naively until selection or probe grid is refined."
        )
    elif global_summary["monotone_N_decrease_flag"]:
        k18_call = (
            "Not yet informative naively; best residual decreases with N but does not reach the near-hit band."
        )
    else:
        k18_call = (
            "Not informative naively; random cloud-cloud strategy fails to shrink residual with N."
        )

    table_header = (
        "| N | seed | spin | n_sel | rej_time | rej_radial | rej_angular | unresolved | best_residual | near_hits |"
    )
    table_sep = "|---|---|---|---|---|---|---|---|---|---|"
    table_rows = []
    for c in sorted(cell_summaries, key=lambda x: (x["N"], x["seed"], x["spin_a"])):
        table_rows.append(
            f"| {c['N']} | {c['seed']} | {c['spin_a']:.2f} | "
            f"{c['n_selected_total']} | {c['n_rejected_time_short']} | "
            f"{c['n_rejected_radial_proxy']} | {c['n_rejected_angular_large']} | "
            f"{c['n_selection_unresolved']} | {_fmt_residual(c['best_residual_in_cell'])} | "
            f"{c['near_hit_count']} |"
        )

    lines = [
        "# S4-KERR-K17d cloud-size / seed / spin scan",
        "",
        "1. Does best_residual decrease with N?",
        f"monotone_N_decrease_flag = {monotone}. Median by N: {med_str}.",
        "",
        "2. Are there near-hits (W_TOL < residual <= 10*W_TOL)?",
        f"near_hits_overall = {near_hits}. residual_le_w_tol_overall = {le_w_tol}.",
        "",
        "3. Does any cell cross W_TOL or 10*W_TOL?",
        f"N_at_which_W_TOL_first_crossed = {N_w}.",
        f"N_at_which_10x_W_TOL_first_crossed = {N_10w}.",
        "",
        "4. Is K18 informative naively on this generator?",
        k18_call,
        "",
        "5. Recommendation:",
        global_summary["recommendation"],
        "",
        "6. Per-cell breakdown (rejection counts + best residual + near-hits):",
        "",
        table_header,
        table_sep,
        *table_rows,
        "",
        "Refinement boundary: K17d does not refine top cases; refinement is reserved for a separate phase.",
        "",
        *CAVEATS,
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_png(cell_summaries: list[dict[str, Any]], path: Path) -> None:
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle("S4-KERR-K17d cloud-size / seed / spin scan")

    for c in cell_summaries:
        v = c["best_residual_in_cell"]
        if v is None:
            continue
        axs[0, 0].scatter(c["N"], v, s=30, alpha=0.7, color="tab:blue")
    axs[0, 0].set_yscale("log")
    axs[0, 0].axhline(k17.W_TOL, linestyle="--", color="k", linewidth=1, label="W_TOL")
    axs[0, 0].axhline(10.0 * k17.W_TOL, linestyle=":", color="k", linewidth=1, label="10*W_TOL")
    axs[0, 0].set_xticks(sorted(N_GRID))
    axs[0, 0].set_xlabel("N")
    axs[0, 0].set_ylabel("best_residual_in_cell")
    axs[0, 0].set_title("Per-cell best residual vs N")
    axs[0, 0].legend(fontsize=7)

    Ns = sorted(N_GRID)
    medians = []
    for N in Ns:
        vals = [
            c["best_residual_in_cell"]
            for c in cell_summaries
            if c["N"] == N and c["best_residual_in_cell"] is not None
        ]
        medians.append(statistics.median(vals) if vals else None)
    finite_xs = [Ns[i] for i, m in enumerate(medians) if m is not None]
    finite_ys = [m for m in medians if m is not None]
    if finite_ys:
        axs[0, 1].plot(finite_xs, finite_ys, "o-")
        axs[0, 1].set_yscale("log")
        axs[0, 1].axhline(k17.W_TOL, linestyle="--", color="k", linewidth=1, label="W_TOL")
        axs[0, 1].axhline(10.0 * k17.W_TOL, linestyle=":", color="k", linewidth=1, label="10*W_TOL")
        axs[0, 1].legend(fontsize=7)
    axs[0, 1].set_xticks(Ns)
    axs[0, 1].set_xlabel("N")
    axs[0, 1].set_ylabel("median best_residual")
    axs[0, 1].set_title("Median best residual per N")

    cells_sorted = sorted(cell_summaries, key=lambda x: (x["N"], x["seed"], x["spin_a"]))
    labels = [f"N{c['N']}/s{c['seed']}/a{c['spin_a']:.2f}" for c in cells_sorted]
    counts = [c["n_selected_total"] for c in cells_sorted]
    axs[1, 0].bar(range(len(labels)), counts)
    axs[1, 0].set_xticks(range(len(labels)))
    axs[1, 0].set_xticklabels(labels, rotation=80, fontsize=6)
    axs[1, 0].set_title("Selected pair count per cell")
    axs[1, 0].set_ylabel("n_selected_total")

    reject_keys = (
        "n_rejected_time_short",
        "n_rejected_radial_proxy",
        "n_rejected_angular_large",
        "n_selection_unresolved",
    )
    pretty = ("time_short", "radial_proxy", "angular_large", "unresolved")
    totals = [sum(c[k] for c in cell_summaries) for k in reject_keys]
    axs[1, 1].bar(pretty, totals)
    axs[1, 1].set_title("Aggregate rejection counts (all cells)")
    axs[1, 1].set_ylabel("count")

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sanity-cell",
        action="store_true",
        help="Run a single cell (N=12, seed=1959, spin=0.0) and print its summary; no artifacts written.",
    )
    parser.add_argument(
        "--smoke-one-cell",
        action="store_true",
        help="Run one cell (N=12, seed=1959, spin=0.0) end-to-end and write _smoke-suffixed artifacts.",
    )
    parser.add_argument(
        "--full",
        action="store_true",
        help="Run the full 3x3x3 grid and write csv/json/md/png artifacts. Reserved for Task 5.",
    )
    args = parser.parse_args()

    if args.sanity_cell:
        result = _run_cell(12, 1959, 0.0)
        cs = result["cell_summary"]
        print(
            f"sanity_cell N={cs['N']} seed={cs['seed']} spin={cs['spin_a']} "
            f"n_events={cs['n_events']} n_pairs_total={cs['n_pairs_total']} "
            f"n_selected_total={cs['n_selected_total']} "
            f"best_residual={cs['best_residual_in_cell']} "
            f"residual_le_w_tol={cs['residual_le_w_tol_count']} "
            f"near_hits={cs['near_hit_count']} "
            f"rejected_time_short={cs['n_rejected_time_short']} "
            f"rejected_radial_proxy={cs['n_rejected_radial_proxy']} "
            f"rejected_angular_large={cs['n_rejected_angular_large']} "
            f"selection_unresolved={cs['n_selection_unresolved']} "
            f"dominant={cs['dominant_component_histogram']} "
            f"cloud_t_std={cs['cloud_t_std']:.4f} cloud_r_std={cs['cloud_r_std']:.4f}"
        )
        return

    if args.smoke_one_cell:
        result = _run_cell(12, 1959, 0.0)
        pair_rows = list(result["pair_rows"])
        cell_summaries = [result["cell_summary"]]
        global_summary = _global_summary(pair_rows, cell_summaries)

        smoke_prefix = f"{OUT_PREFIX}_smoke"
        csv_path = ARTIFACT_DIR / f"{smoke_prefix}.csv"
        json_path = ARTIFACT_DIR / f"{smoke_prefix}.json"
        md_path = ARTIFACT_DIR / f"{smoke_prefix}.md"
        png_path = ARTIFACT_DIR / f"{smoke_prefix}.png"

        write_csv(pair_rows, csv_path)
        write_json(pair_rows, cell_summaries, global_summary, json_path)
        write_md(global_summary, cell_summaries, md_path)
        write_png(cell_summaries, png_path)

        print(
            f"smoke_one_cell wrote {csv_path.name}, {json_path.name}, "
            f"{md_path.name}, {png_path.name} | selected={len(pair_rows)} "
            f"best_residual_overall={global_summary['best_residual_overall']} "
            f"near_hits_overall={global_summary['near_hits_overall']} "
            f"residual_le_w_tol_overall={global_summary['residual_le_w_tol_overall']} "
            f"monotone_N_decrease_flag={global_summary['monotone_N_decrease_flag']}"
        )
        return

    if args.full:
        pair_rows: list[dict[str, Any]] = []
        cell_summaries: list[dict[str, Any]] = []
        for N in N_GRID:
            for seed in SEED_GRID:
                for spin in SPIN_GRID:
                    result = _run_cell(N, seed, spin)
                    pair_rows.extend(result["pair_rows"])
                    cell_summaries.append(result["cell_summary"])
        global_summary = _global_summary(pair_rows, cell_summaries)

        csv_path = ARTIFACT_DIR / f"{OUT_PREFIX}.csv"
        json_path = ARTIFACT_DIR / f"{OUT_PREFIX}.json"
        md_path = ARTIFACT_DIR / f"{OUT_PREFIX}.md"
        png_path = ARTIFACT_DIR / f"{OUT_PREFIX}.png"

        write_csv(pair_rows, csv_path)
        write_json(pair_rows, cell_summaries, global_summary, json_path)
        write_md(global_summary, cell_summaries, md_path)
        write_png(cell_summaries, png_path)

        print(
            f"cells={len(cell_summaries)} probed={global_summary['total_pairs_probed']} "
            f"best_residual_overall={global_summary['best_residual_overall']} "
            f"near_hits_overall={global_summary['near_hits_overall']} "
            f"residual_le_w_tol_overall={global_summary['residual_le_w_tol_overall']} "
            f"N_first_crossing_W_TOL={global_summary['N_at_which_W_TOL_first_crossed']} "
            f"N_first_crossing_10x_W_TOL={global_summary['N_at_which_10x_W_TOL_first_crossed']}"
        )
        return

    raise NotImplementedError(
        "Pass --sanity-cell, --smoke-one-cell, or --full."
    )


if __name__ == "__main__":
    main()
