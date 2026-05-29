#!/usr/bin/env python3
"""S4-KERR-K18C-RESIDUAL-BOUNDARY-SENSITIVITY-001.

Local boundary-sensitivity diagnostic for exactly ONE K18B candidate family.

Scientific question (the only one this script answers):
    For the fixed candidate
        N=24, seed=1961, spin_a=0.50, event_A=15, event_B=4, direction=ingoing,
    does the K17d endpoint-weighted residual have an INTERIOR minimum in the
    probe controls (b, lambda), or does it keep decreasing monotonically toward
    / outside the previous K17d grid boundary (best_b=1.0 upper edge,
    best_lambda=0.5 lower edge)?

This is NOT:
    - a causal classifier;
    - a reachability claim;
    - a new random cloud-cloud search;
    - a new necessary/sufficient filter;
    - a production geometry change;
    - an attempt to improve the candidate.

It reuses, without reimplementation, the *exact* K17d residual evaluation:
    k17._eval_trial(spin=, A=, B=, b=, lam=, direction=)
called precisely as K17d._probe_best calls it. The probe controls b and lambda
are the only things extended beyond the K17d grid; direction is held at the
K18B best_direction (ingoing); the sector m is optimized internally by
_eval_trial exactly as in K17d (its returned best_sector_m is recorded).

Guardrails (also written into the artifact):
    - residual profile != causal_true
    - interior minimum != reachability
    - runaway residual != spacelike separation
    - K18C is a numerical objective diagnostic only
    - no global causal claim is made

Run:
    python3 -m explore.sorkin4_kerr_benchmark.audit_kerr_k18c_residual_boundary_sensitivity_001

It writes (only):
    explore/sorkin4_kerr_benchmark/kerr_k18c_residual_boundary_sensitivity_001.md
    explore/sorkin4_kerr_benchmark/kerr_k18c_residual_boundary_sensitivity_001.json
No CSV/PNG are emitted; this is a small local diagnostic, not a sweep.
"""

from __future__ import annotations

import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

# Reuse the same modules and the same residual engine as K17d. No geometry is
# reimplemented here; cones.py / production geometry are untouched.
from explore.sorkin4_kerr_benchmark import run_kerr_minimal_benchmark as kerr  # noqa: E402
from explore.sorkin4_kerr_benchmark import (  # noqa: E402
    audit_kerr_k17_controlled_candidate_pair_sandbox_001 as k17,
)
from explore.sorkin4_schwarzschild_benchmark import (  # noqa: E402
    run_schwarzschild_minimal_benchmark as schwarz,
)

ARTIFACT_DIR = Path(__file__).resolve().parent
OUT_PREFIX = "kerr_k18c_residual_boundary_sensitivity_001"

# ----------------------------------------------------------------------------
# Fixed candidate (from K18B / K18A rank-1 family). Do not change.
# ----------------------------------------------------------------------------
MASS = 1.0
EXTERIOR_MARGIN = schwarz.EXTERIOR_MARGIN

CAND_N = 24
CAND_SEED = 1961
CAND_SPIN = 0.50
CAND_A_INDEX = 15
CAND_B_INDEX = 4
CAND_DIRECTION = -1.0  # ingoing  (K17b convention: outgoing if dir>0 else ingoing)
CAND_DIRECTION_LABEL = "ingoing"

# Previous K17d probe lattice (for reference / boundary identification only):
K17D_PROBE_B_GRID = (-1.0, -0.5, 0.0, 0.5, 1.0)  # best_b = 1.0  -> upper edge
K17D_PROBE_LAMBDA_GRID = (0.5, 1.0, 2.0)  # best_lambda = 0.5 -> lower edge
K17D_BEST_B = 1.0
K17D_BEST_LAMBDA = 0.5

# Provenance reference values recorded by K18B for this exact family at a=0.50.
# Used only as a sanity guard that we regenerated the identical cloud/pair.
REF_RESIDUAL_AT_BEST = 0.12546941635537312
REF_A = {"t": 2.9103567417252916, "r": 4.917680370166319, "phi": 4.9572224569140015}
REF_B = {"t": 3.8983298507599984, "r": 4.550588789137428, "phi": 5.106522579152854}
PROVENANCE_TOL = 1.0e-9

# ----------------------------------------------------------------------------
# Extended probe grids. These EXTEND beyond the K17d boundaries on purpose.
#   - b-cut (lambda fixed at the old lower edge 0.5): extend b ABOVE 1.0.
#   - lambda-cut (b fixed at the old upper edge 1.0): extend lambda BELOW 0.5.
#   - small 2D corner around (b=1.0, lambda=0.5).
# Old grid points are retained so the cut reproduces the K17d optimum exactly.
# ----------------------------------------------------------------------------
B_CUT_GRID = (
    -0.5,
    0.0,
    0.5,
    1.0,
    1.25,
    1.5,
    2.0,
    2.5,
    3.0,
    4.0,
)  # extension: increasing b
B_CUT_FIXED_LAMBDA = 0.5

LAMBDA_CUT_GRID = (
    0.03125,
    0.0625,
    0.125,
    0.25,
    0.5,
    1.0,
    2.0,
)  # extension: decreasing lambda
LAMBDA_CUT_FIXED_B = 1.0

CORNER_B_GRID = (0.5, 1.0, 1.5, 2.0)
CORNER_LAMBDA_GRID = (0.125, 0.25, 0.5, 1.0)

CAVEATS = [
    "residual profile != causal_true",
    "interior minimum != reachability",
    "runaway residual != spacelike separation",
    "K18C is a numerical objective diagnostic only",
    "no global causal claim is made",
    "boundary-sensitivity of a residual is a statement about the objective's "
    "parametrization, not about Kerr causal structure",
    "this audit is local to ONE K18A/K18B candidate family",
]


# ----------------------------------------------------------------------------
# Reused residual evaluation: identical call site to K17d._probe_best.
# ----------------------------------------------------------------------------
def _eval_residual(
    spin: float,
    A: "kerr.Event",
    B: "kerr.Event",
    b: float,
    lam: float,
    direction: float,
) -> dict[str, Any]:
    """Single residual evaluation, mirroring K17d._probe_best exactly.

    Returns the K17d weighted residual plus the internally chosen sector and the
    per-axis residuals. No re-definition of the residual is introduced.
    """
    trial = k17._eval_trial(spin=spin, A=A, B=B, b=b, lam=lam, direction=direction)
    w = float(trial.get("endpoint_weighted_residual", float("inf")))
    return {
        "b": b,
        "lambda": lam,
        "weighted_residual": w,
        "best_sector_m": trial.get("best_sector_m"),
        "t_residual": trial.get("endpoint_t_residual"),
        "r_residual": trial.get("endpoint_r_residual"),
        "phi_residual_sector_adjusted": trial.get(
            "endpoint_phi_residual_sector_adjusted"
        ),
        "direction_best": trial.get("direction_best"),
    }


def _build_fixed_pair() -> dict[str, Any]:
    """Regenerate the exact K17d cloud for the fixed candidate and pull A, B."""
    r_plus = kerr.kerr_horizon_radius(MASS, CAND_SPIN)
    events = kerr.generate_exterior_events(
        CAND_N, CAND_SEED, r_plus + EXTERIOR_MARGIN, equatorial=True
    )
    emap = {e.index: e for e in events}
    A = emap[CAND_A_INDEX]
    B = emap[CAND_B_INDEX]
    return {"r_plus": r_plus, "A": A, "B": B, "n_events": len(events)}


def _provenance_check(A: "kerr.Event", B: "kerr.Event") -> dict[str, Any]:
    def _close(x: float, y: float) -> bool:
        return math.isfinite(x) and math.isfinite(y) and abs(x - y) <= PROVENANCE_TOL

    a_ok = (
        _close(A.t, REF_A["t"])
        and _close(A.r, REF_A["r"])
        and _close(A.phi, REF_A["phi"])
    )
    b_ok = (
        _close(B.t, REF_B["t"])
        and _close(B.r, REF_B["r"])
        and _close(B.phi, REF_B["phi"])
    )
    return {
        "event_A_matches_k18b": a_ok,
        "event_B_matches_k18b": b_ok,
        "A": {"t": A.t, "r": A.r, "phi": A.phi},
        "B": {"t": B.t, "r": B.r, "phi": B.phi},
    }


# ----------------------------------------------------------------------------
# Pre-registered classification of a 1D cut.
# ----------------------------------------------------------------------------
def _classify_cut(
    points: list[dict[str, Any]], param_key: str, extension_is_increasing: bool
) -> str:
    """INTERIOR / RUNAWAY / INCONCLUSIVE for one 1D residual cut.

    RUNAWAY  : the minimum residual sits at the extreme grid point in the
               extension direction (residual keeps improving off the new edge).
    INTERIOR : the minimum residual is strictly interior (or on the
               non-extension side), i.e. extending past the old K17d edge does
               not keep lowering the residual.
    INCONCLUSIVE: fewer than 3 finite residuals -> profile not safely evaluable.
    """
    finite = [
        (p[param_key], p["weighted_residual"])
        for p in points
        if math.isfinite(p["weighted_residual"])
    ]
    if len(finite) < 3:
        return "INCONCLUSIVE"
    p_at_min = min(finite, key=lambda t: t[1])[0]
    if extension_is_increasing:
        extreme_p = max(p for p, _ in finite)
    else:
        extreme_p = min(p for p, _ in finite)
    if abs(p_at_min - extreme_p) <= 1.0e-12:
        return "RUNAWAY"
    return "INTERIOR"


def _combine(class_b: str, class_lambda: str) -> str:
    if "INCONCLUSIVE" in (class_b, class_lambda):
        return "INCONCLUSIVE"
    if "RUNAWAY" in (class_b, class_lambda):
        return "RUNAWAY"
    return "INTERIOR"


# ----------------------------------------------------------------------------
# Run the diagnostic.
# ----------------------------------------------------------------------------
def run() -> dict[str, Any]:
    pair = _build_fixed_pair()
    A, B = pair["A"], pair["B"]
    prov = _provenance_check(A, B)

    b_cut = [
        _eval_residual(
            CAND_SPIN, A, B, b=b, lam=B_CUT_FIXED_LAMBDA, direction=CAND_DIRECTION
        )
        for b in B_CUT_GRID
    ]
    lambda_cut = [
        _eval_residual(
            CAND_SPIN, A, B, b=LAMBDA_CUT_FIXED_B, lam=lam, direction=CAND_DIRECTION
        )
        for lam in LAMBDA_CUT_GRID
    ]
    corner = [
        _eval_residual(CAND_SPIN, A, B, b=b, lam=lam, direction=CAND_DIRECTION)
        for b in CORNER_B_GRID
        for lam in CORNER_LAMBDA_GRID
    ]

    # Reproduction of the K18B optimum at (b=1.0, lambda=0.5, ingoing).
    repro = _eval_residual(
        CAND_SPIN, A, B, b=K17D_BEST_B, lam=K17D_BEST_LAMBDA, direction=CAND_DIRECTION
    )
    repro_ok = (
        math.isfinite(repro["weighted_residual"])
        and abs(repro["weighted_residual"] - REF_RESIDUAL_AT_BEST) <= 1.0e-9
    )

    class_b = _classify_cut(b_cut, "b", extension_is_increasing=True)
    class_lambda = _classify_cut(lambda_cut, "lambda", extension_is_increasing=False)
    overall = _combine(class_b, class_lambda)

    provenance_ok = (
        prov["event_A_matches_k18b"] and prov["event_B_matches_k18b"] and repro_ok
    )
    if not provenance_ok:
        # Cannot safely trust the profile if we did not reproduce the exact pair
        # and the exact K18B optimum residual.
        overall = "INCONCLUSIVE"

    return {
        "pair": pair,
        "provenance": prov,
        "reproduction": {
            "recomputed_residual_at_best": repro["weighted_residual"],
            "reference_residual_at_best": REF_RESIDUAL_AT_BEST,
            "reproduced": repro_ok,
        },
        "b_cut": b_cut,
        "lambda_cut": lambda_cut,
        "corner": corner,
        "classification_b_cut": class_b,
        "classification_lambda_cut": class_lambda,
        "classification_overall": overall,
        "provenance_ok": provenance_ok,
    }


# ----------------------------------------------------------------------------
# Markdown / JSON artifact emission.
# ----------------------------------------------------------------------------
def _fmt(v: Any, nd: int = 12) -> str:
    if v is None:
        return "None"
    if isinstance(v, float):
        if not math.isfinite(v):
            return "inf"
        return f"{v:.{nd}g}"
    return str(v)


def _b_cut_table(points: list[dict[str, Any]]) -> list[str]:
    lines = [
        f"Fixed: lambda = {B_CUT_FIXED_LAMBDA}, direction = {CAND_DIRECTION_LABEL}. "
        f"Extension direction: increasing b beyond old upper edge {K17D_BEST_B}.",
        "",
        "| b | weighted_residual | best_sector_m | t_residual | r_residual | phi_residual_sector_adj | note |",
        "|---:|---:|---:|---:|---:|---:|---|",
    ]
    for p in points:
        note = []
        if abs(p["b"] - K17D_BEST_B) <= 1e-12:
            note.append("old best_b (upper edge)")
        if p["b"] in K17D_PROBE_B_GRID:
            note.append("in K17d grid")
        else:
            note.append("extension")
        lines.append(
            f"| {_fmt(p['b'])} | {_fmt(p['weighted_residual'])} | {_fmt(p['best_sector_m'])} | "
            f"{_fmt(p['t_residual'])} | {_fmt(p['r_residual'])} | "
            f"{_fmt(p['phi_residual_sector_adjusted'])} | {'; '.join(note)} |"
        )
    return lines


def _lambda_cut_table(points: list[dict[str, Any]]) -> list[str]:
    lines = [
        f"Fixed: b = {LAMBDA_CUT_FIXED_B}, direction = {CAND_DIRECTION_LABEL}. "
        f"Extension direction: decreasing lambda below old lower edge {K17D_BEST_LAMBDA}.",
        "",
        "| lambda | weighted_residual | best_sector_m | t_residual | r_residual | phi_residual_sector_adj | note |",
        "|---:|---:|---:|---:|---:|---:|---|",
    ]
    for p in points:
        note = []
        if abs(p["lambda"] - K17D_BEST_LAMBDA) <= 1e-12:
            note.append("old best_lambda (lower edge)")
        if p["lambda"] in K17D_PROBE_LAMBDA_GRID:
            note.append("in K17d grid")
        else:
            note.append("extension")
        lines.append(
            f"| {_fmt(p['lambda'])} | {_fmt(p['weighted_residual'])} | {_fmt(p['best_sector_m'])} | "
            f"{_fmt(p['t_residual'])} | {_fmt(p['r_residual'])} | "
            f"{_fmt(p['phi_residual_sector_adjusted'])} | {'; '.join(note)} |"
        )
    return lines


def _corner_table(points: list[dict[str, Any]]) -> list[str]:
    lines = [
        f"Direction = {CAND_DIRECTION_LABEL}. Small 2D grid around the corner "
        f"(b={K17D_BEST_B}, lambda={K17D_BEST_LAMBDA}). Cell = weighted_residual.",
        "",
    ]
    header = (
        "| b \\ lambda | " + " | ".join(_fmt(la) for la in CORNER_LAMBDA_GRID) + " |"
    )
    sep = "|---:" * (len(CORNER_LAMBDA_GRID) + 1) + "|"
    lines.append(header)
    lines.append(sep)
    cell = {(p["b"], p["lambda"]): p["weighted_residual"] for p in points}
    for b in CORNER_B_GRID:
        row = [f"{_fmt(b)}"]
        for la in CORNER_LAMBDA_GRID:
            row.append(_fmt(cell.get((b, la), float("inf"))))
        lines.append("| " + " | ".join(row) + " |")
    return lines


def _interpretation(result: dict[str, Any]) -> list[str]:
    overall = result["classification_overall"]
    if not result["provenance_ok"]:
        return [
            "Provenance/reproduction failed: the regenerated pair or the recomputed "
            "residual at (b=1.0, lambda=0.5, ingoing) did not match the recorded K18B "
            "values within tolerance. The boundary profile is therefore not trusted and "
            "the result is INCONCLUSIVE. No closing/continuation decision is supported.",
        ]
    if overall == "RUNAWAY":
        return [
            "The residual keeps decreasing monotonically toward the extended grid edge in "
            "at least one probe control; no interior minimum is attained within the explored "
            "b/lambda range. Under the pre-registered reading this indicates the K18B "
            "candidate is consistent with a boundary/parametrization artifact of the K17d "
            "residual objective rather than an interior geometric optimum. This SUPPORTS "
            "closing the residual cloud-cloud path on objective-wellposedness grounds. It "
            "says nothing about Kerr causal structure (see guardrails).",
        ]
    if overall == "INTERIOR":
        return [
            "The residual attains an interior minimum within the extended b/lambda range "
            "(it stops improving once the old K17d edge is passed). Under the pre-registered "
            "reading this indicates the old K17d probe grid was too narrow, so K17d/K18A "
            "residual rankings are NOT final and would need re-evaluation on a corrected "
            "grid. This is an objective-grid statement only; no causal/reachability claim "
            "follows (see guardrails).",
        ]
    return [
        "Result INCONCLUSIVE: fewer than three finite residuals on a cut, so the profile "
        "shape cannot be classified safely without further (still non-production, "
        "non-geometry) work.",
    ]


def write_artifacts(result: dict[str, Any]) -> tuple[Path, Path]:
    md_path = ARTIFACT_DIR / f"{OUT_PREFIX}.md"
    json_path = ARTIFACT_DIR / f"{OUT_PREFIX}.json"

    prov = result["provenance"]
    repro = result["reproduction"]

    lines: list[str] = []
    lines += ["# S4-KERR-K18C residual boundary-sensitivity diagnostic 001", ""]

    lines += ["## Status", ""]
    lines += [
        f"Completed (single-candidate boundary-sensitivity diagnostic). "
        f"Overall classification: **{result['classification_overall']}** "
        f"(b-cut: {result['classification_b_cut']}, lambda-cut: {result['classification_lambda_cut']}).",
        "",
    ]

    lines += ["## Fixed candidate", ""]
    lines += [
        f"- N = {CAND_N}",
        f"- seed = {CAND_SEED}",
        f"- spin_a = {CAND_SPIN}",
        f"- event_A = {CAND_A_INDEX}",
        f"- event_B = {CAND_B_INDEX}",
        f"- direction = {CAND_DIRECTION_LABEL} (probe direction held fixed at K18B best_direction)",
        f"- K17d best_b = {K17D_BEST_B} (upper edge of old PROBE_B_GRID)",
        f"- K17d best_lambda = {K17D_BEST_LAMBDA} (lower edge of old PROBE_LAMBDA_GRID)",
        "- sector m: optimized internally by k17._eval_trial exactly as in K17d (recorded as best_sector_m).",
        "",
    ]

    lines += ["## Input artifacts / code actually used", ""]
    lines += [
        "- `explore/sorkin4_kerr_benchmark/audit_kerr_k17_controlled_candidate_pair_sandbox_001.py` "
        "(`k17._eval_trial` — the residual evaluation, reused unchanged)",
        "- `explore/sorkin4_kerr_benchmark/audit_kerr_k17d_cloud_size_seed_scan_001.py` "
        "(`_probe_best` call pattern and PROBE grids reproduced)",
        "- `explore/sorkin4_kerr_benchmark/run_kerr_minimal_benchmark.py` "
        "(`kerr_horizon_radius`, `generate_exterior_events`, `Event`)",
        "- `explore/sorkin4_schwarzschild_benchmark/run_schwarzschild_minimal_benchmark.py` "
        "(`EXTERIOR_MARGIN`)",
        "- `explore/sorkin4_kerr_benchmark/kerr_k18b_single_candidate_geometric_audit_001.md` "
        "(candidate provenance / reference values)",
        "",
        "cones.py and production geometry were NOT modified or imported beyond the existing "
        "K17 audit stack.",
        "",
    ]

    lines += ["## Exact commands used", ""]
    lines += [
        "```bash",
        "python3 -m explore.sorkin4_kerr_benchmark.audit_kerr_k18c_residual_boundary_sensitivity_001",
        "```",
        "",
    ]

    lines += ["## Provenance / reproduction guard", ""]
    lines += [
        f"- regenerated cloud: N={CAND_N}, seed={CAND_SEED}, spin_a={CAND_SPIN}, "
        f"n_events={result['pair']['n_events']}, r_plus={_fmt(result['pair']['r_plus'])}",
        f"- event_A coords (t,r,phi) = "
        f"({_fmt(prov['A']['t'])}, {_fmt(prov['A']['r'])}, {_fmt(prov['A']['phi'])}) "
        f"— matches K18B: {prov['event_A_matches_k18b']}",
        f"- event_B coords (t,r,phi) = "
        f"({_fmt(prov['B']['t'])}, {_fmt(prov['B']['r'])}, {_fmt(prov['B']['phi'])}) "
        f"— matches K18B: {prov['event_B_matches_k18b']}",
        f"- recomputed residual at (b={K17D_BEST_B}, lambda={K17D_BEST_LAMBDA}, "
        f"{CAND_DIRECTION_LABEL}) = {_fmt(repro['recomputed_residual_at_best'])} "
        f"(K18B reference {_fmt(repro['reference_residual_at_best'])}); "
        f"reproduced = {repro['reproduced']}",
        "",
    ]

    lines += ["## Probe grid used for b and lambda", ""]
    lines += [
        f"- B_CUT_GRID = {B_CUT_GRID} (lambda fixed at {B_CUT_FIXED_LAMBDA})",
        f"- LAMBDA_CUT_GRID = {LAMBDA_CUT_GRID} (b fixed at {LAMBDA_CUT_FIXED_B})",
        f"- CORNER_B_GRID = {CORNER_B_GRID}",
        f"- CORNER_LAMBDA_GRID = {CORNER_LAMBDA_GRID}",
        f"- (reference) K17d PROBE_B_GRID = {K17D_PROBE_B_GRID}",
        f"- (reference) K17d PROBE_LAMBDA_GRID = {K17D_PROBE_LAMBDA_GRID}",
        "",
    ]

    lines += ["## 1D b-cut table", ""]
    lines += _b_cut_table(result["b_cut"])
    lines += ["", f"b-cut classification: **{result['classification_b_cut']}**", ""]

    lines += ["## 1D lambda-cut table", ""]
    lines += _lambda_cut_table(result["lambda_cut"])
    lines += [
        "",
        f"lambda-cut classification: **{result['classification_lambda_cut']}**",
        "",
    ]

    lines += ["## 2D corner table", ""]
    lines += _corner_table(result["corner"])
    lines += [""]

    lines += ["## Classification: INTERIOR / RUNAWAY / INCONCLUSIVE", ""]
    lines += [
        f"- b-cut: **{result['classification_b_cut']}**",
        f"- lambda-cut: **{result['classification_lambda_cut']}**",
        f"- overall: **{result['classification_overall']}**",
        "",
        "Pre-registered rule:",
        "- RUNAWAY if the minimum residual on a cut sits at the extreme grid point in the "
        "extension direction (residual keeps improving off the new edge).",
        "- INTERIOR if the minimum is strictly interior / on the non-extension side "
        "(extending past the old K17d edge does not keep lowering the residual).",
        "- INCONCLUSIVE if a cut has fewer than three finite residuals, or if the "
        "provenance/reproduction guard fails.",
        "- overall = INCONCLUSIVE if any cut is INCONCLUSIVE; else RUNAWAY if any cut is "
        "RUNAWAY; else INTERIOR.",
        "",
    ]

    lines += ["## Interpretation", ""]
    lines += _interpretation(result)
    lines += [""]

    lines += ["## Guardrails", ""]
    lines += [f"- {c}" for c in CAVEATS]
    lines += [""]

    lines += ["## Next operational recommendation", ""]
    overall = result["classification_overall"]
    if not result["provenance_ok"]:
        lines += [
            "Resolve the provenance/reproduction mismatch first (it indicates the regenerated "
            "pair or residual does not match K18B). No close/continue decision until the exact "
            "candidate is reproduced.",
        ]
    elif overall == "RUNAWAY":
        lines += [
            "Single next step: close the residual cloud-cloud path as objective-ill-posed in its "
            "own controls, and (separately, later) reframe what relation the residual should "
            "encode in Kerr. Do NOT add another simple filter (K17F-style) or another random "
            "cloud-cloud search.",
        ]
    elif overall == "INTERIOR":
        lines += [
            "Single next step: treat K17d/K18A residual rankings as non-final and re-evaluate the "
            "top candidates on a corrected (interior-containing) b/lambda grid before any further "
            "interpretation. Still no causal/reachability claim.",
        ]
    else:
        lines += [
            "Single next step: widen only the finite-residual coverage of the failing cut "
            "(still local, still non-production) so the profile becomes classifiable.",
        ]
    lines += [""]

    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    payload = {
        "benchmark": "S4-KERR-K18C residual boundary sensitivity (single candidate)",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "fixed_candidate": {
            "N": CAND_N,
            "seed": CAND_SEED,
            "spin_a": CAND_SPIN,
            "event_A": CAND_A_INDEX,
            "event_B": CAND_B_INDEX,
            "direction": CAND_DIRECTION_LABEL,
            "k17d_best_b": K17D_BEST_B,
            "k17d_best_lambda": K17D_BEST_LAMBDA,
        },
        "grids": {
            "b_cut_grid": list(B_CUT_GRID),
            "b_cut_fixed_lambda": B_CUT_FIXED_LAMBDA,
            "lambda_cut_grid": list(LAMBDA_CUT_GRID),
            "lambda_cut_fixed_b": LAMBDA_CUT_FIXED_B,
            "corner_b_grid": list(CORNER_B_GRID),
            "corner_lambda_grid": list(CORNER_LAMBDA_GRID),
        },
        "provenance": result["provenance"],
        "reproduction": result["reproduction"],
        "b_cut": result["b_cut"],
        "lambda_cut": result["lambda_cut"],
        "corner": result["corner"],
        "classification_b_cut": result["classification_b_cut"],
        "classification_lambda_cut": result["classification_lambda_cut"],
        "classification_overall": result["classification_overall"],
        "provenance_ok": result["provenance_ok"],
        "guardrails": CAVEATS,
    }
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return md_path, json_path


def main() -> None:
    result = run()
    md_path, json_path = write_artifacts(result)
    print(
        f"K18C overall={result['classification_overall']} "
        f"b_cut={result['classification_b_cut']} lambda_cut={result['classification_lambda_cut']} "
        f"provenance_ok={result['provenance_ok']} -> {md_path.name}, {json_path.name}"
    )


if __name__ == "__main__":
    main()
