#!/usr/bin/env python3
"""Phase-space audit for the exterior outgoing turning branch.

This is the C' diagnostic: before comparing arrival times, check whether the
outgoing one-turn branch can share any angular target with the direct no-root
branch.  For each grid point with r2 > r1 > 3M, it records

    phi_direct_max = direct angular reach as c2 -> 1/(27M^2)+
    phi_turning_min = turning angular reach as u_turn -> u1+
    gap = phi_turning_min - phi_direct_max

If the gap is positive, the two branch angular ranges do not overlap at that
grid point.  This is a numerical diagnostic, not a theorem.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

sys.path.insert(0, str(Path(__file__).resolve().parent))

from audit_exterior_turning_branch import branch_angular_ranges  # noqa: E402
from run_schwarzschild_minimal_benchmark import MASS  # noqa: E402


OUT_DIR = Path(__file__).resolve().parent
DEFAULT_OUT_PREFIX = "schwarzschild_exterior_turning_phase_space_audit"
DEFAULT_R1_MIN = 3.1
DEFAULT_R1_MAX = 12.0
DEFAULT_R2_MAX = 20.0
DEFAULT_STEP = 0.1

CSV_FIELDS = (
    "r1",
    "r2",
    "u1",
    "u2",
    "phi_direct_max",
    "phi_turning_min",
    "turning_minus_direct_phi_gap",
    "status",
)


def _jsonable(value: Any) -> Any:
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return None
    return value


def _grid_values(start: float, stop: float, step: float) -> list[float]:
    scale = 10 ** max(0, -math.floor(math.log10(step)) + 2)
    start_i = math.ceil(start * scale - 1.0e-9)
    stop_i = math.floor(stop * scale + 1.0e-9)
    step_i = max(1, round(step * scale))
    return [value / scale for value in range(start_i, stop_i + 1, step_i)]


def audit_grid(
    r1_min: float = DEFAULT_R1_MIN,
    r1_max: float = DEFAULT_R1_MAX,
    r2_max: float = DEFAULT_R2_MAX,
    step: float = DEFAULT_STEP,
    mass: float = MASS,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    photon_radius = 3.0 * mass
    for r1 in _grid_values(max(r1_min, photon_radius + step), r1_max, step):
        for r2 in _grid_values(r1 + step, r2_max, step):
            if not (r2 > r1 > photon_radius):
                continue
            u1 = 1.0 / r1
            u2 = 1.0 / r2
            direct_max, turning_min, gap = branch_angular_ranges(u1, u2, mass)
            if direct_max is None or turning_min is None or gap is None:
                status = "range_failed"
            elif gap > 0.0:
                status = "disjoint_ranges"
            else:
                status = "range_overlap"
            rows.append(
                {
                    "r1": r1,
                    "r2": r2,
                    "u1": u1,
                    "u2": u2,
                    "phi_direct_max": direct_max,
                    "phi_turning_min": turning_min,
                    "turning_minus_direct_phi_gap": gap,
                    "status": status,
                }
            )
    return rows


def write_outputs(
    rows: list[dict[str, Any]],
    out_prefix: str,
    r1_min: float,
    r1_max: float,
    r2_max: float,
    step: float,
) -> tuple[Path, Path, Path]:
    csv_path = OUT_DIR / f"{out_prefix}.csv"
    json_path = OUT_DIR / f"{out_prefix}.json"
    md_path = OUT_DIR / f"{out_prefix}.md"

    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in CSV_FIELDS})

    status_counts = {status: sum(1 for row in rows if row["status"] == status) for status in sorted({row["status"] for row in rows})}
    gaps = [
        row["turning_minus_direct_phi_gap"]
        for row in rows
        if row["turning_minus_direct_phi_gap"] is not None
    ]
    min_gap_row = min(
        (row for row in rows if row["turning_minus_direct_phi_gap"] is not None),
        key=lambda row: row["turning_minus_direct_phi_gap"],
        default=None,
    )
    summary = {
        "audit": "S4 Schwarzschild exterior outgoing turning-branch phase space",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "r1_min": r1_min,
        "r1_max": r1_max,
        "r2_max": r2_max,
        "step": step,
        "grid_points": len(rows),
        "status_counts": status_counts,
        "min_turning_minus_direct_phi_gap": min(gaps) if gaps else None,
        "max_turning_minus_direct_phi_gap": max(gaps) if gaps else None,
        "min_gap_row": {key: _jsonable(value) for key, value in min_gap_row.items()} if min_gap_row else None,
        "all_ranges_disjoint": bool(rows) and status_counts == {"disjoint_ranges": len(rows)},
        "scope": (
            "Numerical phase-space audit only. It checks whether the direct "
            "no-root branch and outgoing one-turn branch have overlapping "
            "angular ranges on a finite r-grid with r2 > r1 > 3M. It is not "
            "a proof of the global fastest-geodesic lemma."
        ),
    }
    json_path.write_text(
        json.dumps({"summary": summary, "rows": [{key: _jsonable(value) for key, value in row.items()} for row in rows]}, indent=2),
        encoding="utf-8",
    )

    md_lines = [
        "# Schwarzschild Exterior Turning-Branch Phase-Space Audit",
        "",
        "This is the C' numerical diagnostic for the outgoing one-turn competitor.",
        "It checks angular-range overlap before any arrival-time comparison.",
        "",
        f"- r1 range: {r1_min}..{r1_max}",
        f"- r2 max: {r2_max}",
        f"- step: {step}",
        f"- Grid points: {len(rows)}",
        f"- Status counts: {status_counts}",
        f"- All ranges disjoint: {summary['all_ranges_disjoint']}",
        f"- Min angular gap: {summary['min_turning_minus_direct_phi_gap']}",
        f"- Min-gap row: {summary['min_gap_row']}",
        "",
        "Definitions:",
        "",
        "- `phi_direct_max`: direct branch angular reach as `c2 -> 1/(27M^2)+`.",
        "- `phi_turning_min`: one-turn branch angular reach as `u_turn -> u1+`.",
        "- Positive gap means the two angular ranges do not overlap at that grid point.",
    ]
    md_path.write_text("\n".join(md_lines) + "\n", encoding="utf-8")
    return csv_path, json_path, md_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit outgoing turning-branch angular range over r-space.")
    parser.add_argument("--r1-min", type=float, default=DEFAULT_R1_MIN)
    parser.add_argument("--r1-max", type=float, default=DEFAULT_R1_MAX)
    parser.add_argument("--r2-max", type=float, default=DEFAULT_R2_MAX)
    parser.add_argument("--step", type=float, default=DEFAULT_STEP)
    parser.add_argument("--out-prefix", default=DEFAULT_OUT_PREFIX)
    args = parser.parse_args()
    if args.step <= 0.0:
        raise SystemExit("--step must be positive")
    if args.r1_max <= args.r1_min:
        raise SystemExit("--r1-max must be greater than --r1-min")
    if args.r2_max <= args.r1_min:
        raise SystemExit("--r2-max must be greater than --r1-min")
    return args


def main() -> None:
    args = parse_args()
    rows = audit_grid(args.r1_min, args.r1_max, args.r2_max, args.step)
    csv_path, json_path, md_path = write_outputs(
        rows,
        args.out_prefix,
        args.r1_min,
        args.r1_max,
        args.r2_max,
        args.step,
    )
    print(f"wrote {csv_path}")
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")
    print(f"grid_points={len(rows)}")


if __name__ == "__main__":
    main()
