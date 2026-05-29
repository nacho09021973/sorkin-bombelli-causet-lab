#!/usr/bin/env python3
"""S4-KERR-K18D-RESIDUAL-CALIBRATION-001.

Residual scale calibration for the K18A/K18C best candidate.

Scientific question (the only one this script answers):
    Is endpoint_weighted_residual ~ 0.125 close to the smallest value the
    objective can attain for a known-reachable synthetic target (the floor),
    or is 0.125 effectively background-level?

This script implements Phases 0 and 1 only:
    Phase 0: provenance guard (reproduce exact cloud + residual).
    Phase 1: floor measurement via synthetic reachable targets.

Phase 2 (background distribution) and Phase 3 (verdict) are NOT run here.
The floor numbers must be reviewed before proceeding.

This is NOT:
    - a causal classifier;
    - a reachability claim;
    - a new random cloud-cloud search;
    - a geometry change of any kind.

Guardrails (carried into the artifact):
    - residual profile != causal_true
    - interior minimum != reachability
    - low residual != proof
    - synthetic reachable target != causal_true (it is a numerical self-consistency anchor)
    - candidate_hit != reachability; candidate_miss != spacelike separation
    - do NOT touch cones.py

Run:
    python3 -m explore.sorkin4_kerr_benchmark.audit_kerr_k18d_residual_calibration_001

Writes (partial, Phase 0+1 only):
    explore/sorkin4_kerr_benchmark/kerr_k18d_residual_calibration_001.md
    explore/sorkin4_kerr_benchmark/kerr_k18d_residual_calibration_001.json
"""

from __future__ import annotations

import json
import math
import sys
import types
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from explore.sorkin4_kerr_benchmark import run_kerr_minimal_benchmark as kerr  # noqa: E402
from explore.sorkin4_kerr_benchmark import (  # noqa: E402
    audit_kerr_k17_controlled_candidate_pair_sandbox_001 as k17,
)
from explore.sorkin4_schwarzschild_benchmark import (  # noqa: E402
    run_schwarzschild_minimal_benchmark as schwarz,
)

ARTIFACT_DIR = Path(__file__).resolve().parent
OUT_PREFIX = "kerr_k18d_residual_calibration_001"

# ----------------------------------------------------------------------------
# Fixed point (identical to K18C; do not change).
# ----------------------------------------------------------------------------
MASS = k17.MASS
EXTERIOR_MARGIN = schwarz.EXTERIOR_MARGIN

CAND_N = 24
CAND_SEED = 1961
CAND_SPIN = 0.50
CAND_A_INDEX = 15
CAND_B_INDEX = 4
CAND_DIRECTION = -1.0  # ingoing
CAND_DIRECTION_LABEL = "ingoing"

REF_RESIDUAL_AT_BEST = 0.12546941635537312
REF_A = {"t": 2.9103567417252916, "r": 4.917680370166319, "phi": 4.9572224569140015}
REF_B = {"t": 3.8983298507599984, "r": 4.550588789137428, "phi": 5.106522579152854}
PROVENANCE_TOL = 1.0e-9

# Previous K17d best (b, lambda) — provenance reproduction anchor.
K17D_BEST_B = 1.0
K17D_BEST_LAMBDA = 0.5

# ----------------------------------------------------------------------------
# Calibration grid (K18C interior-containing; from plan Section 4).
# ----------------------------------------------------------------------------
CAL_B_GRID = (0.6, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4)
CAL_LAMBDA_GRID = (0.25, 0.375, 0.5, 0.625, 0.75, 1.0)
# direction is fixed at ingoing for all CAL evaluations.

# ----------------------------------------------------------------------------
# Phase 1 construction set (plan Section 5, Phase 1).
# on-grid: the synthetic target is the integrated endpoint of the exact
#          K17d best-fit (b*, lambda*) pair — expected to recover ~float noise.
# off-grid: plausible interior points not coinciding with CAL grid nodes.
# ----------------------------------------------------------------------------
CONSTRUCTIONS_ON_GRID = [
    {"b_star": 1.0, "lambda_star": 0.5},
]
CONSTRUCTIONS_OFF_GRID = [
    {"b_star": 0.95, "lambda_star": 0.45},
    {"b_star": 1.05, "lambda_star": 0.55},
    {"b_star": 1.1, "lambda_star": 0.6},
]

CAVEATS = [
    "residual profile != causal_true",
    "interior minimum != reachability",
    "low residual != proof",
    "synthetic reachable target != causal_true (it is a numerical self-consistency anchor)",
    "candidate_hit != reachability; candidate_miss != spacelike separation",
    "K18D is a numerical objective calibration only; it makes no causal claim",
    "floor_off_grid measures the floor the probe faces for a generic target at CAL spacing",
    "floor_on_grid measures the integrator/float-noise floor (expected near machine epsilon)",
]


# ----------------------------------------------------------------------------
# Helpers.
# ----------------------------------------------------------------------------
def _is_edge(b: float, lam: float) -> bool:
    """Return True if (b, lam) sits on a CAL grid boundary."""
    b_min, b_max = min(CAL_B_GRID), max(CAL_B_GRID)
    l_min, l_max = min(CAL_LAMBDA_GRID), max(CAL_LAMBDA_GRID)
    return (
        abs(b - b_min) < 1e-12
        or abs(b - b_max) < 1e-12
        or abs(lam - l_min) < 1e-12
        or abs(lam - l_max) < 1e-12
    )


def _probe_best_cal(
    A: "kerr.Event",
    B_target: Any,
    direction: float,
) -> dict[str, Any]:
    """probe_best over the K18D CAL grid at fixed direction.

    Returns the grid minimum residual, the best (b, lambda), and an EDGE flag.
    B_target only needs .t, .r, .phi (may be SimpleNamespace).
    """
    best_residual = float("inf")
    best_b: float | None = None
    best_lam: float | None = None
    for b in CAL_B_GRID:
        for lam in CAL_LAMBDA_GRID:
            trial = k17._eval_trial(
                spin=CAND_SPIN, A=A, B=B_target, b=b, lam=lam, direction=direction
            )
            w = float(trial.get("endpoint_weighted_residual", float("inf")))
            if w < best_residual:
                best_residual = w
                best_b = b
                best_lam = lam
    edge = (
        _is_edge(best_b, best_lam)
        if (best_b is not None and best_lam is not None)
        else False
    )
    return {
        "best_residual": best_residual,
        "best_b": best_b,
        "best_lambda": best_lam,
        "edge": edge,
    }


# ----------------------------------------------------------------------------
# Phase 0: provenance guard.
# ----------------------------------------------------------------------------
def phase0() -> dict[str, Any]:
    r_plus = kerr.kerr_horizon_radius(MASS, CAND_SPIN)
    events = kerr.generate_exterior_events(
        CAND_N, CAND_SEED, r_plus + EXTERIOR_MARGIN, equatorial=True
    )
    emap = {e.index: e for e in events}
    A = emap[CAND_A_INDEX]
    B_real = emap[CAND_B_INDEX]

    def _close(x: float, y: float) -> bool:
        return math.isfinite(x) and math.isfinite(y) and abs(x - y) <= PROVENANCE_TOL

    a_ok = (
        _close(A.t, REF_A["t"])
        and _close(A.r, REF_A["r"])
        and _close(A.phi, REF_A["phi"])
    )
    b_ok = (
        _close(B_real.t, REF_B["t"])
        and _close(B_real.r, REF_B["r"])
        and _close(B_real.phi, REF_B["phi"])
    )

    repro_trial = k17._eval_trial(
        spin=CAND_SPIN,
        A=A,
        B=B_real,
        b=K17D_BEST_B,
        lam=K17D_BEST_LAMBDA,
        direction=CAND_DIRECTION,
    )
    repro_res = float(repro_trial.get("endpoint_weighted_residual", float("inf")))
    repro_ok = (
        math.isfinite(repro_res)
        and abs(repro_res - REF_RESIDUAL_AT_BEST) <= PROVENANCE_TOL
    )

    provenance_ok = a_ok and b_ok and repro_ok

    return {
        "r_plus": r_plus,
        "n_events": len(events),
        "A": {"t": A.t, "r": A.r, "phi": A.phi},
        "B_real": {"t": B_real.t, "r": B_real.r, "phi": B_real.phi},
        "event_A_matches_k18c": a_ok,
        "event_B_matches_k18c": b_ok,
        "recomputed_residual": repro_res,
        "reference_residual": REF_RESIDUAL_AT_BEST,
        "residual_reproduced": repro_ok,
        "provenance_ok": provenance_ok,
        "_events": events,
        "_A": A,
        "_B_real": B_real,
    }


# ----------------------------------------------------------------------------
# Phase 1: floor (synthetic reachable targets).
# ----------------------------------------------------------------------------
def _run_construction(
    A: "kerr.Event",
    b_star: float,
    lambda_star: float,
) -> dict[str, Any]:
    run = k17.integrate_to_lambda(
        spin=CAND_SPIN,
        b=b_star,
        direction=CAND_DIRECTION,
        state0=(A.t, A.r, A.phi),
        lambda_end=lambda_star,
    )
    if run["failed_reason"] is not None:
        return {
            "b_star": b_star,
            "lambda_star": lambda_star,
            "status": "failed",
            "failed_reason": run["failed_reason"],
            "B_syn": None,
            "probe_best_residual": None,
            "probe_best_b": None,
            "probe_best_lambda": None,
            "edge": None,
        }

    t_f, r_f, phi_f = run["states"][-1]
    B_syn = types.SimpleNamespace(t=t_f, r=r_f, phi=phi_f, index=-1)

    probe = _probe_best_cal(A, B_syn, direction=CAND_DIRECTION)

    return {
        "b_star": b_star,
        "lambda_star": lambda_star,
        "status": "ok",
        "failed_reason": None,
        "B_syn": {"t": t_f, "r": r_f, "phi": phi_f},
        "probe_best_residual": probe["best_residual"],
        "probe_best_b": probe["best_b"],
        "probe_best_lambda": probe["best_lambda"],
        "edge": probe["edge"],
    }


def phase1(A: "kerr.Event") -> dict[str, Any]:
    on_grid = [
        _run_construction(A, c["b_star"], c["lambda_star"])
        for c in CONSTRUCTIONS_ON_GRID
    ]
    off_grid = [
        _run_construction(A, c["b_star"], c["lambda_star"])
        for c in CONSTRUCTIONS_OFF_GRID
    ]

    def _min_residual(rows: list[dict[str, Any]]) -> float | None:
        vals = [
            r["probe_best_residual"]
            for r in rows
            if r["status"] == "ok" and r["probe_best_residual"] is not None
        ]
        return min(vals) if vals else None

    floor_on_grid = _min_residual(on_grid)
    floor_off_grid = _min_residual(off_grid)

    return {
        "on_grid_constructions": on_grid,
        "off_grid_constructions": off_grid,
        "floor_on_grid": floor_on_grid,
        "floor_off_grid": floor_off_grid,
    }


# ----------------------------------------------------------------------------
# Artifact emission.
# ----------------------------------------------------------------------------
def _fmt(v: Any, nd: int = 12) -> str:
    if v is None:
        return "None"
    if isinstance(v, float):
        if not math.isfinite(v):
            return "inf"
        return f"{v:.{nd}g}"
    return str(v)


def _construction_table(rows: list[dict[str, Any]], label: str) -> list[str]:
    lines = [
        f"### {label}",
        "",
        "| b_star | lambda_star | status | probe_best_residual | probe_best_b | probe_best_lambda | edge |",
        "|---:|---:|---|---:|---:|---:|---|",
    ]
    for r in rows:
        lines.append(
            f"| {_fmt(r['b_star'])} | {_fmt(r['lambda_star'])} | {r['status']} | "
            f"{_fmt(r['probe_best_residual'])} | {_fmt(r['probe_best_b'])} | "
            f"{_fmt(r['probe_best_lambda'])} | {_fmt(r['edge'])} |"
        )
    return lines


def _phase1_interpretation(p1: dict[str, Any]) -> list[str]:
    floor_on = p1["floor_on_grid"]
    floor_off = p1["floor_off_grid"]
    threshold = 1.0e-2
    lines = []
    if floor_off is None:
        lines.append(
            "All off-grid constructions failed: floor_off_grid is undefined. "
            "Cannot classify. Verdict trending INCONCLUSIVE."
        )
    elif floor_off >= threshold:
        lines.append(
            f"floor_off_grid = {_fmt(floor_off)} >= 1e-2 (near-hit threshold). "
            "The probe cannot recover even a known-reachable off-grid target below the "
            "near-hit scale. This indicates OBJECTIVE_UNRESOLVED: the residual does not "
            "measure at the required scale, and the 0.125 value carries no information. "
            "Phase 2+3 would be moot — recommend stopping here."
        )
    else:
        lines.append(
            f"floor_off_grid = {_fmt(floor_off)} < 1e-2. "
            "The objective can recover off-grid synthetic targets below the near-hit scale. "
            "The floor is resolved. Phase 2 (background comparison) is warranted: "
            "the question of whether B_real is distinguished from background remains open."
        )
    if floor_on is not None:
        lines.append(
            f"floor_on_grid = {_fmt(floor_on)} (integrator/float-noise floor; "
            "expected near machine epsilon for on-grid construction)."
        )
    return lines


def write_artifacts(p0: dict[str, Any], p1: dict[str, Any]) -> tuple[Path, Path]:
    md_path = ARTIFACT_DIR / f"{OUT_PREFIX}.md"
    json_path = ARTIFACT_DIR / f"{OUT_PREFIX}.json"

    provenance_ok = p0["provenance_ok"]
    floor_on = p1["floor_on_grid"]
    floor_off = p1["floor_off_grid"]
    threshold = 1.0e-2

    if not provenance_ok:
        status_str = "PARTIAL — Phase 0 FAILED (INCONCLUSIVE; provenance guard failed)"
    elif floor_off is None:
        status_str = "PARTIAL — Phase 0 ok; Phase 1 INCONCLUSIVE (all off-grid constructions failed)"
    elif floor_off >= threshold:
        status_str = (
            f"PARTIAL — Phase 0 ok; Phase 1 complete; "
            f"floor_off_grid = {_fmt(floor_off)} >= 1e-2 → trending OBJECTIVE_UNRESOLVED"
        )
    else:
        status_str = (
            f"PARTIAL — Phase 0 ok; Phase 1 complete; "
            f"floor_off_grid = {_fmt(floor_off)} < 1e-2 → Phase 2 warranted"
        )

    lines: list[str] = []
    lines += ["# S4-KERR-K18D residual calibration audit 001 (Phase 0 + Phase 1)", ""]

    lines += ["## Status", ""]
    lines += [status_str, ""]
    lines += [
        "This is a PARTIAL artifact. Phase 2 (background distribution) and Phase 3 "
        "(pre-registered verdict) are not yet run. Review floor numbers before proceeding.",
        "",
    ]

    lines += ["## Fixed point", ""]
    lines += [
        f"- N = {CAND_N}",
        f"- seed = {CAND_SEED}",
        f"- spin_a = {CAND_SPIN}",
        f"- event_A = {CAND_A_INDEX}",
        f"- event_B = {CAND_B_INDEX}  (B_real; used in Phase 2 only, not used here)",
        f"- direction = {CAND_DIRECTION_LABEL} (fixed throughout)",
        f"- Reference residual (b=1.0, lambda=0.5, ingoing) = {REF_RESIDUAL_AT_BEST}",
        "",
    ]

    lines += ["## Definitions", ""]
    lines += [
        "- **objective(A, target, b, lambda, direction)** = K17 weighted residual "
        "`max(|dt|, |dr|, |dphi_sector_adjusted|)` of the integrated null-geodesic "
        "endpoint from A versus target, with sector m optimized internally. "
        "Exact reuse of `k17._eval_trial(...)['endpoint_weighted_residual']`.",
        "- **probe_best(A, target, grid, direction)** = `min` of `objective` over "
        "the CAL (b, lambda) grid at fixed direction.",
        "- **floor_on_grid** = min `probe_best` over on-grid constructions "
        "(b_star, lambda_star coincide with CAL grid nodes). Measures integrator/float-noise.",
        "- **floor_off_grid** = min `probe_best` over off-grid constructions "
        "(b_star, lambda_star do not coincide with CAL grid nodes). Measures "
        "the realistic floor a generic target faces given CAL spacing.",
        "",
    ]

    lines += ["## Code reused", ""]
    lines += [
        "- `explore/sorkin4_kerr_benchmark/run_kerr_minimal_benchmark.py` "
        "(`kerr_horizon_radius`, `generate_exterior_events`, `Event`)",
        "- `explore/sorkin4_kerr_benchmark/audit_kerr_k17_controlled_candidate_pair_sandbox_001.py` "
        "(`_eval_trial`, `integrate_to_lambda`, `W_TOL`, `MASS`)",
        "- `explore/sorkin4_schwarzschild_benchmark/run_schwarzschild_minimal_benchmark.py` "
        "(`EXTERIOR_MARGIN`)",
        "- K18C provenance constants reproduced (same fixed point, same reference values).",
        "",
        "`cones.py` and production geometry were NOT modified or imported.",
        "",
    ]

    lines += ["## Exact commands", ""]
    lines += [
        "```bash",
        "python3 -m explore.sorkin4_kerr_benchmark.audit_kerr_k18d_residual_calibration_001",
        "```",
        "",
    ]

    lines += ["## Calibration grid", ""]
    lines += [
        f"- CAL_B_GRID = {CAL_B_GRID}",
        f"- CAL_LAMBDA_GRID = {CAL_LAMBDA_GRID}",
        f"- direction = {CAND_DIRECTION_LABEL} (fixed)",
        "- Edge guard: if best (b, lambda) sits on a CAL grid boundary, result is flagged EDGE. "
        "Do NOT auto-expand the grid; surface it for the next decision.",
        "",
    ]

    lines += ["## Phase 0: Provenance guard", ""]
    lines += [
        f"- r_plus = {_fmt(p0['r_plus'])}",
        f"- n_events regenerated = {p0['n_events']}",
        f"- A coords (t,r,phi) = ({_fmt(p0['A']['t'])}, {_fmt(p0['A']['r'])}, {_fmt(p0['A']['phi'])}) "
        f"— matches K18C: {p0['event_A_matches_k18c']}",
        f"- B_real coords (t,r,phi) = ({_fmt(p0['B_real']['t'])}, {_fmt(p0['B_real']['r'])}, {_fmt(p0['B_real']['phi'])}) "
        f"— matches K18C: {p0['event_B_matches_k18c']}",
        f"- recomputed residual at (b={K17D_BEST_B}, lambda={K17D_BEST_LAMBDA}, {CAND_DIRECTION_LABEL}) = "
        f"{_fmt(p0['recomputed_residual'])} (reference: {_fmt(p0['reference_residual'])}); "
        f"reproduced = {p0['residual_reproduced']}",
        f"- **provenance_ok = {p0['provenance_ok']}**",
        "",
    ]
    if not provenance_ok:
        lines += [
            "STOP: Provenance guard failed. The regenerated pair or recomputed residual does not match "
            "the K18C/K18B reference within tolerance. Overall verdict = INCONCLUSIVE. "
            "Phase 1 was not run.",
            "",
        ]

    if provenance_ok:
        lines += ["## Phase 1: Floor (synthetic reachable targets)", ""]
        lines += [
            "Each construction integrates from A at (b_star, lambda_star, ingoing) to get B_syn, "
            "then runs probe_best over the full CAL grid (ingoing only). "
            "The floor is the minimum probe_best residual across constructions of each type.",
            "",
        ]
        lines += _construction_table(
            p1["on_grid_constructions"], "On-grid constructions"
        )
        lines += [""]
        lines += _construction_table(
            p1["off_grid_constructions"], "Off-grid constructions"
        )
        lines += [
            "",
            f"**floor_on_grid  = {_fmt(floor_on)}**",
            f"**floor_off_grid = {_fmt(floor_off)}**",
            f"(near-hit threshold = 1e-2; W_TOL = {k17.W_TOL})",
            "",
        ]
        lines += ["### Phase 1 interpretation", ""]
        lines += _phase1_interpretation(p1)
        lines += [""]

    lines += ["## Guardrails", ""]
    lines += [f"- {c}" for c in CAVEATS]
    lines += [""]

    lines += ["## Next operational recommendation", ""]
    if not provenance_ok:
        lines += [
            "Resolve the provenance mismatch before any other step. "
            "No floor/background/verdict decision until Phase 0 passes.",
        ]
    elif floor_off is None:
        lines += [
            "All off-grid Phase 1 constructions failed: investigate whether integration "
            "is failing for this spin/direction at these (b*, lambda*) values before "
            "proceeding. No Phase 2 until at least one off-grid floor value is available.",
        ]
    elif floor_off >= threshold:
        lines += [
            f"floor_off_grid = {_fmt(floor_off)} >= 1e-2. The objective is unresolved at "
            "the required scale. Recommended action: close or redefine the residual cloud-cloud "
            "line (OBJECTIVE_UNRESOLVED verdict; Phase 2+3 would be moot). "
            "Do NOT re-rank candidates.",
        ]
    else:
        lines += [
            f"floor_off_grid = {_fmt(floor_off)} < 1e-2. The floor is resolved. "
            "Proceed to Phase 2: compare B_real against background from the same N=24 cloud "
            "(all cloud events T with T.t > A.t and T exterior, using probe_best on the CAL grid). "
            "Review the B_real rank in background before any verdict.",
        ]
    lines += [""]

    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    payload: dict[str, Any] = {
        "benchmark": "S4-KERR-K18D residual calibration (Phase 0 + Phase 1)",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "phases_run": ["phase0", "phase1"] if provenance_ok else ["phase0"],
        "status": status_str,
        "fixed_point": {
            "N": CAND_N,
            "seed": CAND_SEED,
            "spin_a": CAND_SPIN,
            "event_A": CAND_A_INDEX,
            "event_B": CAND_B_INDEX,
            "direction": CAND_DIRECTION_LABEL,
            "reference_residual": REF_RESIDUAL_AT_BEST,
        },
        "calibration_grid": {
            "CAL_B_GRID": list(CAL_B_GRID),
            "CAL_LAMBDA_GRID": list(CAL_LAMBDA_GRID),
            "direction": CAND_DIRECTION_LABEL,
        },
        "phase0": {
            "r_plus": p0["r_plus"],
            "n_events": p0["n_events"],
            "A": p0["A"],
            "B_real": p0["B_real"],
            "event_A_matches_k18c": p0["event_A_matches_k18c"],
            "event_B_matches_k18c": p0["event_B_matches_k18c"],
            "recomputed_residual": p0["recomputed_residual"],
            "reference_residual": p0["reference_residual"],
            "residual_reproduced": p0["residual_reproduced"],
            "provenance_ok": p0["provenance_ok"],
        },
        "phase1": {
            "on_grid_constructions": p1["on_grid_constructions"],
            "off_grid_constructions": p1["off_grid_constructions"],
            "floor_on_grid": p1["floor_on_grid"],
            "floor_off_grid": p1["floor_off_grid"],
        }
        if provenance_ok
        else None,
        "guardrails": CAVEATS,
    }
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return md_path, json_path


def main() -> None:
    p0 = phase0()

    if not p0["provenance_ok"]:
        p1: dict[str, Any] = {
            "on_grid_constructions": [],
            "off_grid_constructions": [],
            "floor_on_grid": None,
            "floor_off_grid": None,
        }
        md_path, json_path = write_artifacts(p0, p1)
        print(
            f"K18D Phase0 FAILED (provenance_ok=False) -> {md_path.name}, {json_path.name}"
        )
        return

    A = p0["_A"]
    p1 = phase1(A)

    md_path, json_path = write_artifacts(p0, p1)

    floor_on = p1["floor_on_grid"]
    floor_off = p1["floor_off_grid"]
    print(
        f"K18D Phase0+1 done: provenance_ok={p0['provenance_ok']} "
        f"floor_on_grid={_fmt(floor_on)} floor_off_grid={_fmt(floor_off)} "
        f"-> {md_path.name}, {json_path.name}"
    )


if __name__ == "__main__":
    main()
