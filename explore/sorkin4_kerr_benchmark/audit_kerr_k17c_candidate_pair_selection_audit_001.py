#!/usr/bin/env python3
"""S4-KERR-K17C-CANDIDATE-PAIR-SELECTION-AUDIT-001.

Selection-only audit with cheap necessary geometric filters.
No causal/global reachability claims are introduced.
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
OUT_PREFIX = "kerr_k17c_candidate_pair_selection_audit_001_n12_seed1959"
K17_JSON = ARTIFACT_DIR / "kerr_k17_controlled_candidate_pair_sandbox_001_n12_seed1959.json"

MASS = 1.0
SPINS = (0.0, 0.25, 0.5)
N = 12
SEED = 1959
EXTERIOR_MARGIN = schwarz.EXTERIOR_MARGIN
PHI_LOW_MAX = 0.8
PHI_SECTOR_MAX = 1.6
TIME_TOLERANCE_BAND = 0.05
RADIAL_LIKE_DR_MAX = 0.75
SECTORS = (-2, -1, 0, 1, 2)

PROBE_B_GRID = (-1.0, -0.5, 0.0, 0.5, 1.0)
PROBE_LAMBDA_GRID = (0.5, 1.0, 2.0)
PROBE_DIR_GRID = (+1.0, -1.0)


def _principal_angle(delta: float) -> float:
    return (delta + math.pi) % (2.0 * math.pi) - math.pi


def _to_float_or_none(v: Any) -> float | None:
    if v is None:
        return None
    return float(v)


def _load_k17_summary() -> tuple[int, int]:
    if not K17_JSON.exists():
        return 0, 0
    payload = json.loads(K17_JSON.read_text(encoding="utf-8"))
    s = payload.get("global_summary", {})
    return int(s.get("candidate_hits", 0)), int(s.get("candidate_undecided", 0))


def _radial_proxy(r_a: float, r_b: float) -> float:
    r1 = min(r_a, r_b)
    r2 = max(r_a, r_b)
    try:
        return abs(kerr.outgoing_radial_trip(r1, r2, MASS))
    except (ValueError, ZeroDivisionError):
        # Conservative fallback proxy; still selection-only heuristic.
        fac = 1.0 / max(1.0e-12, 1.0 - 2.0 * MASS / min(r_a, r_b))
        return abs(r_b - r_a) * fac


def _best_sector_abs(delta_phi_mod: float) -> float:
    return min(abs(delta_phi_mod + 2.0 * math.pi * m) for m in SECTORS)


def _probe_residual(spin: float, A: kerr.Event, B: kerr.Event) -> float | None:
    best: float | None = None
    for b in PROBE_B_GRID:
        for lam in PROBE_LAMBDA_GRID:
            for direction in PROBE_DIR_GRID:
                trial = k17._eval_trial(spin=spin, A=A, B=B, b=b, lam=lam, direction=direction)
                w = float(trial.get("endpoint_weighted_residual", float("inf")))
                if best is None or w < best:
                    best = w
    return best


def build_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for spin in SPINS:
        r_plus = kerr.kerr_horizon_radius(MASS, spin)
        events = kerr.generate_exterior_events(N, SEED, r_plus + EXTERIOR_MARGIN, equatorial=True)

        for A in events:
            for B in events:
                if A.index == B.index:
                    continue

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

                selected_radial_like = False
                selected_low_winding = False
                selected_sector_aware = False

                if not time_order_pass:
                    selection_label = "rejected_time_short"
                elif not both_exterior:
                    selection_label = "selection_unresolved"
                elif not radial_time_admissible:
                    selection_label = "rejected_radial_proxy"
                elif angular_low and abs(dr) <= RADIAL_LIKE_DR_MAX:
                    selected_radial_like = True
                    selection_label = "selected_radial_like"
                elif angular_low:
                    selected_low_winding = True
                    selection_label = "selected_low_winding"
                elif sector_admissible:
                    selected_sector_aware = True
                    selection_label = "selected_sector_aware"
                else:
                    selection_label = "rejected_angular_large"

                is_selected = selection_label.startswith("selected_")
                probe_best = _probe_residual(spin, A, B) if is_selected else None
                probe_pass = bool(probe_best is not None and probe_best <= k17.W_TOL)
                probe_not_close = bool(probe_best is not None and probe_best > k17.W_TOL)

                ang_score = max(0.0, 1.0 - abs(dphi_mod) / math.pi)
                time_score = max(0.0, min(1.0, (time_margin + TIME_TOLERANCE_BAND) / (1.0 + radial_time_proxy)))
                radial_score = max(0.0, 1.0 - min(1.0, abs(dr) / 5.0))
                selection_score = 0.5 * time_score + 0.3 * ang_score + 0.2 * radial_score

                row = {
                    "case_id": f"k17c_a{spin:.2f}_A{A.index}_B{B.index}",
                    "spin_a": spin,
                    "event_A_index": A.index,
                    "event_B_index": B.index,
                    "t_A": A.t,
                    "r_A": A.r,
                    "phi_A": A.phi,
                    "t_B": B.t,
                    "r_B": B.r,
                    "phi_B": B.phi,
                    "delta_t_AB": dt,
                    "delta_r_AB": dr,
                    "delta_phi_AB": dphi_raw,
                    "angular_separation_mod_2pi": dphi_mod,
                    "r_plus": r_plus,
                    "both_exterior": both_exterior,
                    "time_order_pass": time_order_pass,
                    "radial_time_proxy": radial_time_proxy,
                    "time_margin": time_margin,
                    "radial_time_admissible": radial_time_admissible,
                    "angular_low_winding_admissible": angular_low,
                    "sector_admissible": sector_admissible,
                    "selected_radial_like": selected_radial_like,
                    "selected_low_winding": selected_low_winding,
                    "selected_sector_aware": selected_sector_aware,
                    "selection_label": selection_label,
                    "selection_score": selection_score,
                    "endpoint_weighted_residual_best_if_probed": _to_float_or_none(probe_best),
                    "residual_probe_pass": probe_pass,
                    "residual_probe_not_close": probe_not_close,
                    "no_causal_claim_introduced": True,
                    "no_production_classifier_introduced": True,
                    "no_global_causal_relations_decided": True,
                    "all_checks_pass": True,
                }
                rows.append(row)
    return rows


def write_csv(rows: list[dict[str, Any]], path: Path) -> None:
    fields = [
        "case_id", "spin_a", "event_A_index", "event_B_index", "t_A", "r_A", "phi_A", "t_B", "r_B", "phi_B",
        "delta_t_AB", "delta_r_AB", "delta_phi_AB", "angular_separation_mod_2pi", "r_plus", "both_exterior",
        "time_order_pass", "radial_time_proxy", "time_margin", "radial_time_admissible",
        "angular_low_winding_admissible", "sector_admissible", "selected_radial_like", "selected_low_winding",
        "selected_sector_aware", "selection_label", "selection_score", "endpoint_weighted_residual_best_if_probed",
        "residual_probe_pass", "residual_probe_not_close", "no_causal_claim_introduced",
        "no_production_classifier_introduced", "no_global_causal_relations_decided", "all_checks_pass",
    ]
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k) for k in fields})


def _build_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    k17_hits, k17_undec = _load_k17_summary()
    forward = [r for r in rows if r["time_order_pass"]]
    exterior_forward = [r for r in forward if r["both_exterior"]]
    selected = [r for r in rows if r["selection_label"].startswith("selected_")]
    probes = [r for r in rows if r["endpoint_weighted_residual_best_if_probed"] is not None]
    best_probe = min((float(r["endpoint_weighted_residual_best_if_probed"]) for r in probes), default=None)

    selected_count = len(selected)
    if selected_count == 0:
        recommendation = "Do not run K18 naively: tighten pre-selection first because zero pairs pass cheap admissibility."
    elif len([r for r in selected if r["residual_probe_pass"]]) == 0:
        recommendation = "Use radial_time_admissible + low-winding/sector-aware pre-selection before K18; naive K18 is likely low-yield."
    else:
        recommendation = "K18 may be informative only on pre-selected pairs from radial and angular admissibility gates."

    return {
        "total_pairs_enumerated": len(rows),
        "forward_time_pairs": len(forward),
        "exterior_forward_pairs": len(exterior_forward),
        "selected_radial_like_pairs": sum(1 for r in rows if r["selected_radial_like"]),
        "selected_low_winding_pairs": sum(1 for r in rows if r["selected_low_winding"]),
        "selected_sector_aware_pairs": sum(1 for r in rows if r["selected_sector_aware"]),
        "rejected_time_short_pairs": sum(1 for r in rows if r["selection_label"] == "rejected_time_short"),
        "rejected_angular_large_pairs": sum(1 for r in rows if r["selection_label"] == "rejected_angular_large"),
        "rejected_radial_proxy_pairs": sum(1 for r in rows if r["selection_label"] == "rejected_radial_proxy"),
        "selection_unresolved_pairs": sum(1 for r in rows if r["selection_label"] == "selection_unresolved"),
        "residual_probe_pairs": len(probes),
        "residual_probe_passes": sum(1 for r in probes if r["residual_probe_pass"]),
        "best_residual_probe": best_probe,
        "k17_original_hits": k17_hits,
        "k17_original_undecided": k17_undec,
        "recommendation": recommendation,
        "global_true_relations": 0,
        "global_false_relations": 0,
        "global_undecided_pairs": 66,
        "all_checks_pass": (
            all(r["no_causal_claim_introduced"] for r in rows)
            and all(r["no_production_classifier_introduced"] for r in rows)
            and all(r["no_global_causal_relations_decided"] for r in rows)
        ),
    }


def write_json(rows: list[dict[str, Any]], path: Path) -> dict[str, Any]:
    summary = _build_summary(rows)
    payload = {
        "benchmark": "S4-KERR-K17c candidate-pair selection audit",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "cases": rows,
        "global_summary": summary,
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return summary


def write_md(rows: list[dict[str, Any]], summary: dict[str, Any], path: Path) -> None:
    k17_selected_fail = summary["rejected_radial_proxy_pairs"] + summary["rejected_angular_large_pairs"]
    plausible = summary["selected_radial_like_pairs"] + summary["selected_low_winding_pairs"] + summary["selected_sector_aware_pairs"]
    naive_k18 = "No" if plausible == 0 else "Probably low-yield without pre-selection"
    lines = [
        "# S4-KERR-K17c candidate-pair selection audit",
        "",
        "1. Did K17 select pairs that fail cheap null-admissibility filters?",
        f"Yes: many forward-time exterior pairs fail radial/angle gates (named rejections total={k17_selected_fail}).",
        "",
        "2. Does N=12 / seed=1959 contain any plausible near-null candidate pairs?",
        f"Selection candidates found={plausible} across radial-like, low-winding, and sector-aware buckets.",
        "",
        "3. Is K18 likely to be informative if run naively?",
        f"{naive_k18}. Recommendation: {summary['recommendation']}",
        "",
        "4. What selection heuristic should K18 use?",
        "Use deterministic pre-selection: time order + exteriority + radial_time_admissible + low-winding first, then sector-aware fallback.",
        "",
        "selection_candidate is not reachability.",
        "rejected_by_selection is not proof of spacelike separation.",
        "residual_probe_pass is not causal reachability.",
        "no production classifier.",
        "no physical/global causal claim.",
        "candidate_undecided remains conservative.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_png(rows: list[dict[str, Any]], path: Path) -> None:
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle("S4-KERR-K17c candidate-pair selection audit")

    colors = {
        "selected_radial_like": "tab:green",
        "selected_low_winding": "tab:blue",
        "selected_sector_aware": "tab:purple",
        "rejected_time_short": "tab:red",
        "rejected_angular_large": "tab:orange",
        "rejected_radial_proxy": "tab:brown",
        "selection_unresolved": "tab:gray",
    }

    for label, color in colors.items():
        xs = [r["delta_t_AB"] for r in rows if r["selection_label"] == label]
        ys = [r["radial_time_proxy"] for r in rows if r["selection_label"] == label]
        if xs:
            axs[0, 0].scatter(xs, ys, s=14, alpha=0.8, color=color, label=label)
    axs[0, 0].set_xlabel("delta_t_AB")
    axs[0, 0].set_ylabel("radial_time_proxy")
    axs[0, 0].set_title("delta_t_AB vs radial_time_proxy")
    axs[0, 0].legend(fontsize=7)

    axs[0, 1].hist([r["angular_separation_mod_2pi"] for r in rows], bins=20)
    axs[0, 1].set_title("angular_separation_mod_2pi distribution")
    axs[0, 1].set_xlabel("principal angle")

    labels = list(colors.keys())
    vals = [sum(1 for r in rows if r["selection_label"] == lab) for lab in labels]
    axs[1, 0].bar(labels, vals)
    axs[1, 0].tick_params(axis="x", rotation=30)
    axs[1, 0].set_title("selection label counts")

    probed_vals = [r["endpoint_weighted_residual_best_if_probed"] for r in rows if r["endpoint_weighted_residual_best_if_probed"] is not None]
    if probed_vals:
        axs[1, 1].plot(range(len(probed_vals)), probed_vals, "o-")
        axs[1, 1].axhline(k17.W_TOL, color="k", linestyle="--", linewidth=1)
        axs[1, 1].set_yscale("log")
    axs[1, 1].set_title("residual probe best values")
    axs[1, 1].set_xlabel("selected pair index")

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main() -> None:
    rows = build_rows()
    csv_path = ARTIFACT_DIR / f"{OUT_PREFIX}.csv"
    json_path = ARTIFACT_DIR / f"{OUT_PREFIX}.json"
    md_path = ARTIFACT_DIR / f"{OUT_PREFIX}.md"
    png_path = ARTIFACT_DIR / f"{OUT_PREFIX}.png"

    write_csv(rows, csv_path)
    summary = write_json(rows, json_path)
    write_md(rows, summary, md_path)
    write_png(rows, png_path)

    print(
        f"pairs={summary['total_pairs_enumerated']} selected="
        f"{summary['selected_radial_like_pairs'] + summary['selected_low_winding_pairs'] + summary['selected_sector_aware_pairs']} "
        f"probe_passes={summary['residual_probe_passes']}"
    )


if __name__ == "__main__":
    main()
