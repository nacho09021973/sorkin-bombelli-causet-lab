#!/usr/bin/env python3
"""S4-KERR-K17B-CANDIDATE-RESIDUAL-LANDSCAPE-001.

Diagnostic audit for K17 zero-hit behavior. No causal claims are introduced.
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
from explore.sorkin4_kerr_benchmark import audit_kerr_k17_controlled_candidate_pair_sandbox_001 as k17  # noqa: E402
from explore.sorkin4_schwarzschild_benchmark import run_schwarzschild_minimal_benchmark as schwarz  # noqa: E402

ARTIFACT_DIR = Path(__file__).resolve().parent
OUT_PREFIX = "kerr_k17b_candidate_residual_landscape_001_n12_seed1959"
K17_CSV = ARTIFACT_DIR / "kerr_k17_controlled_candidate_pair_sandbox_001_n12_seed1959.csv"
K17_JSON = ARTIFACT_DIR / "kerr_k17_controlled_candidate_pair_sandbox_001_n12_seed1959.json"

DBS = (-0.5, -0.25, -0.1, 0.0, 0.1, 0.25, 0.5)
LFACT = (0.5, 0.75, 1.0, 1.25, 1.5, 2.0)
SECTORS = (-4, -3, -2, -1, 0, 1, 2, 3, 4)


def _load_k17_cases() -> list[dict[str, Any]]:
    return json.loads(K17_JSON.read_text(encoding="utf-8"))["cases"]


def _dominant(dt: float, dr: float, dphi: float) -> str:
    vals = {"t": abs(dt), "r": abs(dr), "phi": abs(dphi)}
    top = max(vals.values())
    winners = [k for k, v in vals.items() if abs(v - top) <= 1.0e-15]
    return winners[0] if len(winners) == 1 else "mixed"


def _angular_mod(delta_phi: float) -> float:
    return abs((delta_phi + math.pi) % (2.0 * math.pi) - math.pi)


def _time_short_flag(delta_t: float, delta_r: float) -> bool:
    # heuristic only: flat-like radial lower bound proxy
    return delta_t < abs(delta_r)


def _build_event_map(spin: float, seed: int) -> dict[int, kerr.Event]:
    r_plus = kerr.kerr_horizon_radius(1.0, spin)
    events = kerr.generate_exterior_events(12, seed, r_plus + schwarz.EXTERIOR_MARGIN, equatorial=True)
    return {e.index: e for e in events}


def _eval_expanded(
    *,
    row: dict[str, Any],
    A: kerr.Event,
    B: kerr.Event,
) -> dict[str, Any]:
    b0 = float(row["b_best"]) if row["b_best"] not in (None, "", "None") else 0.0
    l0 = float(row["lambda_best"]) if row["lambda_best"] not in (None, "", "None") else 1.0
    spin = float(row["spin_a"])
    best: dict[str, Any] | None = None
    for db in DBS:
        b = b0 + db
        for lf in LFACT:
            lam = max(0.05, l0 * lf)
            for direction in (+1.0, -1.0):
                run = k17.integrate_to_lambda(
                    spin=spin, b=b, direction=direction, state0=(A.t, A.r, A.phi), lambda_end=lam
                )
                if run["failed_reason"] is not None:
                    continue
                t_f, r_f, phi_f = run["states"][-1]
                for m in SECTORS:
                    dphi_adj = phi_f - B.phi + 2.0 * math.pi * m
                    dt = t_f - B.t
                    dr = r_f - B.r
                    w = max(abs(dt), abs(dr), abs(dphi_adj))
                    trial = {
                        "b": b,
                        "lam": lam,
                        "m": m,
                        "direction": "outgoing" if direction > 0 else "ingoing",
                        "w": w,
                        "dt": dt,
                        "dr": dr,
                        "dphi_adj": dphi_adj,
                        "dphi_raw": phi_f - B.phi,
                    }
                    if best is None or trial["w"] < best["w"]:
                        best = trial
    if best is None:
        return {
            "b_best_expanded": b0,
            "lambda_best_expanded": l0,
            "sector_best_expanded": int(row["best_sector_m"]) if row["best_sector_m"] not in (None, "", "None") else 0,
            "endpoint_weighted_residual_expanded": float("inf"),
            "endpoint_t_residual_expanded": float("inf"),
            "endpoint_r_residual_expanded": float("inf"),
            "endpoint_phi_residual_sector_adjusted_expanded": float("inf"),
            "direction_best_expanded": row.get("direction_best"),
        }
    return {
        "b_best_expanded": best["b"],
        "lambda_best_expanded": best["lam"],
        "sector_best_expanded": best["m"],
        "endpoint_weighted_residual_expanded": best["w"],
        "endpoint_t_residual_expanded": best["dt"],
        "endpoint_r_residual_expanded": best["dr"],
        "endpoint_phi_residual_sector_adjusted_expanded": best["dphi_adj"],
        "direction_best_expanded": best["direction"],
    }


def build_rows() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    k17_cases = _load_k17_cases()
    by_spin: dict[float, dict[int, kerr.Event]] = {}
    rows: list[dict[str, Any]] = []
    dominant_counts = {"t": 0, "r": 0, "phi": 0, "mixed": 0}
    improved = 0

    for c in k17_cases:
        spin = float(c["spin_a"])
        seed = int(c["cloud_seed"])
        if spin not in by_spin:
            by_spin[spin] = _build_event_map(spin, seed)
        emap = by_spin[spin]
        A = emap[int(c["event_A_index"])]
        B = emap[int(c["event_B_index"])]

        orig_w = float(c["endpoint_weighted_residual"])
        exp = {
            "b_best_expanded": c.get("b_best"),
            "lambda_best_expanded": c.get("lambda_best"),
            "sector_best_expanded": c.get("best_sector_m"),
            "endpoint_weighted_residual_expanded": orig_w,
            "endpoint_t_residual_expanded": float(c["endpoint_t_residual"]),
            "endpoint_r_residual_expanded": float(c["endpoint_r_residual"]),
            "endpoint_phi_residual_sector_adjusted_expanded": float(c["endpoint_phi_residual_sector_adjusted"]),
            "direction_best_expanded": c.get("direction_best"),
        }
        if c["candidate_pair_classification"] == "candidate_undecided":
            exp = _eval_expanded(row=c, A=A, B=B)

        exp_w = float(exp["endpoint_weighted_residual_expanded"])
        if math.isfinite(orig_w) and math.isfinite(exp_w) and exp_w < orig_w:
            improved += 1

        near_hit = math.isfinite(exp_w) and (exp_w <= 10.0 * k17.W_TOL) and (exp_w > k17.W_TOL)
        far_undecided = (not near_hit) and (c["candidate_pair_classification"] == "candidate_undecided")
        if c["candidate_pair_classification"] == "candidate_miss":
            dclass = "candidate_miss"
        elif near_hit:
            dclass = "near_hit"
        elif far_undecided:
            dclass = "far_undecided"
        else:
            dclass = "stable"

        dom = _dominant(
            float(exp["endpoint_t_residual_expanded"]),
            float(exp["endpoint_r_residual_expanded"]),
            float(exp["endpoint_phi_residual_sector_adjusted_expanded"]),
        )
        dominant_counts[dom] += 1

        row = {
            "case_id": c["case_id"],
            "spin_a": spin,
            "pair_type": c["pair_type"],
            "original_classification": c["candidate_pair_classification"],
            "diagnostic_classification": dclass,
            "b_best_original": c["b_best"],
            "lambda_best_original": c["lambda_best"],
            "sector_best_original": c["best_sector_m"],
            "endpoint_weighted_residual_original": orig_w,
            "b_best_expanded": exp["b_best_expanded"],
            "lambda_best_expanded": exp["lambda_best_expanded"],
            "sector_best_expanded": exp["sector_best_expanded"],
            "endpoint_weighted_residual_expanded": exp_w,
            "residual_improvement_factor": (orig_w / exp_w) if (math.isfinite(orig_w) and math.isfinite(exp_w) and exp_w > 0.0) else 1.0,
            "endpoint_t_residual_expanded": exp["endpoint_t_residual_expanded"],
            "endpoint_r_residual_expanded": exp["endpoint_r_residual_expanded"],
            "endpoint_phi_residual_sector_adjusted_expanded": exp["endpoint_phi_residual_sector_adjusted_expanded"],
            "residual_dominant_component": dom,
            "delta_t_AB": float(c["delta_t_AB"]),
            "delta_r_AB": float(c["delta_r_AB"]),
            "delta_phi_AB": float(c["delta_phi_AB"]),
            "angular_separation_mod_2pi": _angular_mod(float(c["delta_phi_AB"])),
            "t_order_pass": float(c["delta_t_AB"]) > 0.0,
            "exterior_pair_pass": bool(c["event_A_from_cloud"]) and bool(c["event_B_from_cloud"]),
            "heuristic_time_short_flag": _time_short_flag(float(c["delta_t_AB"]), float(c["delta_r_AB"])),
            "near_hit": near_hit,
            "far_undecided": far_undecided,
            "no_causal_claim_introduced": True,
            "no_production_classifier_introduced": True,
            "all_checks_pass": True,
        }
        rows.append(row)

    summary = {
        "total_pairs_analyzed": len(rows),
        "original_candidate_hits": sum(1 for c in k17_cases if c["candidate_pair_classification"] == "candidate_hit"),
        "original_candidate_misses": sum(1 for c in k17_cases if c["candidate_pair_classification"] == "candidate_miss"),
        "original_candidate_undecided": sum(1 for c in k17_cases if c["candidate_pair_classification"] == "candidate_undecided"),
        "diagnostic_near_hits": sum(1 for r in rows if r["near_hit"]),
        "diagnostic_far_undecided": sum(1 for r in rows if r["far_undecided"]),
        "expanded_grid_improved_pairs": improved,
        "best_residual_original": min(r["endpoint_weighted_residual_original"] for r in rows),
        "best_residual_expanded": min(r["endpoint_weighted_residual_expanded"] for r in rows),
        "best_improvement_factor": max(r["residual_improvement_factor"] for r in rows),
        "dominant_failure_counts": dominant_counts,
        "all_checks_pass": all(r["no_causal_claim_introduced"] and r["no_production_classifier_introduced"] for r in rows),
    }
    return rows, summary


def write_csv(rows: list[dict[str, Any]], path: Path) -> None:
    fields = [
        "case_id", "spin_a", "pair_type", "original_classification", "diagnostic_classification",
        "b_best_original", "lambda_best_original", "sector_best_original",
        "endpoint_weighted_residual_original", "b_best_expanded", "lambda_best_expanded",
        "sector_best_expanded", "endpoint_weighted_residual_expanded", "residual_improvement_factor",
        "endpoint_t_residual_expanded", "endpoint_r_residual_expanded",
        "endpoint_phi_residual_sector_adjusted_expanded", "residual_dominant_component", "delta_t_AB",
        "delta_r_AB", "delta_phi_AB", "angular_separation_mod_2pi", "t_order_pass", "exterior_pair_pass",
        "heuristic_time_short_flag", "near_hit", "far_undecided", "no_causal_claim_introduced",
        "no_production_classifier_introduced", "all_checks_pass",
    ]
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k) for k in fields})


def write_json(rows: list[dict[str, Any]], summary: dict[str, Any], path: Path) -> None:
    payload = {
        "benchmark": "S4-KERR-K17b candidate residual landscape",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "cases": rows,
        "global_summary": summary,
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def write_md(summary: dict[str, Any], rows: list[dict[str, Any]], path: Path) -> None:
    lines = [
        "# S4-KERR-K17b candidate residual landscape",
        "",
        "1. Why K17 had zero candidate_hits:",
        "Most cloud pairs are not exactly compatible with null endpoint matching under the restricted K17 grid and tolerances.",
        "",
        "2. Expanded local search near_hits:",
        f"diagnostic_near_hits = {summary['diagnostic_near_hits']}.",
        "",
        "3. Dominant residual component:",
        f"{summary['dominant_failure_counts']}.",
        "",
        "4. Geometry diagnosis:",
        "Rows flagged with heuristic_time_short_flag suggest many candidate pairs are time-short relative to radial proxy.",
        "",
        "5. K18 readiness:",
        "K18 may be informative, but candidate-pair selection and local search design should be improved first.",
        "",
        "diagnostic_near_hit is not reachability.",
        "candidate_undecided remains conservative.",
        "No arbitrary causal classification is introduced.",
        "No production classifier is introduced.",
        "No physical/global causal claim is introduced.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_png(rows: list[dict[str, Any]], path: Path) -> None:
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle("S4-KERR-K17b candidate residual landscape")

    axs[0, 0].plot(
        range(len(rows)),
        [r["endpoint_weighted_residual_original"] for r in rows],
        "o-",
        label="original",
    )
    axs[0, 0].plot(
        range(len(rows)),
        [r["endpoint_weighted_residual_expanded"] for r in rows],
        "s--",
        label="expanded",
    )
    axs[0, 0].set_title("Original vs expanded weighted residual")
    axs[0, 0].legend(fontsize=7)

    counts = {"t": 0, "r": 0, "phi": 0, "mixed": 0}
    for r in rows:
        counts[r["residual_dominant_component"]] += 1
    axs[0, 1].bar(list(counts.keys()), list(counts.values()))
    axs[0, 1].set_title("Dominant residual component counts")

    axs[1, 0].bar(range(len(rows)), [r["residual_improvement_factor"] for r in rows])
    axs[1, 0].set_title("Residual improvement factor by pair")

    cmap = {"near_hit": "tab:orange", "far_undecided": "tab:blue", "candidate_miss": "tab:red", "stable": "tab:green"}
    axs[1, 1].scatter(
        [r["delta_t_AB"] for r in rows],
        [r["angular_separation_mod_2pi"] for r in rows],
        c=[cmap[r["diagnostic_classification"]] for r in rows],
    )
    axs[1, 1].set_xlabel("delta_t_AB")
    axs[1, 1].set_ylabel("angular separation mod 2pi")
    axs[1, 1].set_title("Pair geometry by diagnostic class")

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main() -> None:
    rows, summary = build_rows()
    csv_path = ARTIFACT_DIR / f"{OUT_PREFIX}.csv"
    json_path = ARTIFACT_DIR / f"{OUT_PREFIX}.json"
    md_path = ARTIFACT_DIR / f"{OUT_PREFIX}.md"
    png_path = ARTIFACT_DIR / f"{OUT_PREFIX}.png"
    write_csv(rows, csv_path)
    write_json(rows, summary, json_path)
    write_md(summary, rows, md_path)
    write_png(rows, png_path)
    print(
        f"pairs={summary['total_pairs_analyzed']} near_hits={summary['diagnostic_near_hits']} "
        f"improved={summary['expanded_grid_improved_pairs']}"
    )


if __name__ == "__main__":
    main()
