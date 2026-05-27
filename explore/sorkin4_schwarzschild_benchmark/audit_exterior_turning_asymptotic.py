#!/usr/bin/env python3
"""Asymptotic audit for the exterior outgoing turning-branch angular gap.

This probes the numerically dangerous corner suggested by the C' phase-space
grid: weak field and nearby radii,

    r1 = R,  r2 = R + eps,  R >> M,  eps << R.

The diagnostic records the angular-range gap

    phi_turning_min - phi_direct_max

for a small product grid in R and eps.  It is a numerical asymptotic audit, not
an analytic proof.
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


OUT_DIR = Path(__file__).resolve().parent
DEFAULT_OUT_PREFIX = "schwarzschild_exterior_turning_asymptotic_audit"
DEFAULT_R_VALUES = (6.0, 8.0, 10.0, 12.0, 16.0, 20.0, 30.0, 50.0, 80.0, 120.0)
DEFAULT_EPS_VALUES = (1.0e-3, 3.0e-3, 1.0e-2, 3.0e-2, 1.0e-1, 3.0e-1, 1.0)

CSV_FIELDS = (
    "R",
    "eps",
    "eps_over_R",
    "r1",
    "r2",
    "u1",
    "u2",
    "phi_direct_max",
    "phi_turning_min",
    "turning_minus_direct_phi_gap",
    "gap_times_sqrt_R",
    "gap_times_R",
    "status",
)


def _jsonable(value: Any) -> Any:
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return None
    return value


def parse_float_list(raw: str) -> list[float]:
    values = [float(part.strip()) for part in raw.split(",") if part.strip()]
    if not values:
        raise argparse.ArgumentTypeError("expected at least one comma-separated float")
    return values


def audit_asymptotic_grid(
    r_values: list[float] | tuple[float, ...] = DEFAULT_R_VALUES,
    eps_values: list[float] | tuple[float, ...] = DEFAULT_EPS_VALUES,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for radius in r_values:
        for eps in eps_values:
            r1 = float(radius)
            r2 = r1 + float(eps)
            if not (r2 > r1 > 3.0):
                continue
            u1 = 1.0 / r1
            u2 = 1.0 / r2
            direct_max, turning_min, gap = branch_angular_ranges(u1, u2)
            if direct_max is None or turning_min is None or gap is None:
                status = "range_failed"
            elif gap > 0.0:
                status = "disjoint_ranges"
            else:
                status = "range_overlap"
            rows.append(
                {
                    "R": r1,
                    "eps": eps,
                    "eps_over_R": eps / r1,
                    "r1": r1,
                    "r2": r2,
                    "u1": u1,
                    "u2": u2,
                    "phi_direct_max": direct_max,
                    "phi_turning_min": turning_min,
                    "turning_minus_direct_phi_gap": gap,
                    "gap_times_sqrt_R": gap * math.sqrt(r1) if gap is not None else None,
                    "gap_times_R": gap * r1 if gap is not None else None,
                    "status": status,
                }
            )
    return rows


def write_outputs(
    rows: list[dict[str, Any]],
    out_prefix: str,
    r_values: list[float] | tuple[float, ...],
    eps_values: list[float] | tuple[float, ...],
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
    valid_rows = [row for row in rows if row["turning_minus_direct_phi_gap"] is not None]
    min_gap_row = min(valid_rows, key=lambda row: row["turning_minus_direct_phi_gap"], default=None)
    max_r_min_eps_rows = [
        row
        for row in rows
        if row["R"] == max(r_values) and row["eps"] == min(eps_values)
    ]
    summary = {
        "audit": "S4 Schwarzschild exterior outgoing turning-branch asymptotic corner",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "R_values": list(r_values),
        "eps_values": list(eps_values),
        "grid_points": len(rows),
        "status_counts": status_counts,
        "all_ranges_disjoint": bool(rows) and status_counts == {"disjoint_ranges": len(rows)},
        "min_turning_minus_direct_phi_gap": min_gap_row["turning_minus_direct_phi_gap"] if min_gap_row else None,
        "min_gap_row": {key: _jsonable(value) for key, value in min_gap_row.items()} if min_gap_row else None,
        "largest_R_smallest_eps_row": (
            {key: _jsonable(value) for key, value in max_r_min_eps_rows[0].items()}
            if max_r_min_eps_rows
            else None
        ),
        "scope": (
            "Numerical asymptotic audit only. It probes R large and eps/R small "
            "for r1=R, r2=R+eps, checking whether the angular-range gap remains "
            "positive. It is not an analytic proof."
        ),
    }
    json_path.write_text(
        json.dumps({"summary": summary, "rows": [{key: _jsonable(value) for key, value in row.items()} for row in rows]}, indent=2),
        encoding="utf-8",
    )

    md_lines = [
        "# Schwarzschild Exterior Turning-Branch Asymptotic Audit",
        "",
        "This probes the weak-field, nearby-radius corner suggested by the phase-space audit.",
        "",
        f"- R values: {list(r_values)}",
        f"- eps values: {list(eps_values)}",
        f"- Grid points: {len(rows)}",
        f"- Status counts: {status_counts}",
        f"- All ranges disjoint: {summary['all_ranges_disjoint']}",
        f"- Min angular gap: {summary['min_turning_minus_direct_phi_gap']}",
        f"- Min-gap row: {summary['min_gap_row']}",
        f"- Largest-R / smallest-eps row: {summary['largest_R_smallest_eps_row']}",
        "",
        "Scaling columns:",
        "",
        "- `gap_times_sqrt_R` checks whether the gap is roughly `O(R^-1/2)`.",
        "- `gap_times_R` checks whether the gap is roughly `O(R^-1)`.",
    ]
    md_path.write_text("\n".join(md_lines) + "\n", encoding="utf-8")
    return csv_path, json_path, md_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit the weak-field nearby-radius angular-gap corner.")
    parser.add_argument("--R-values", type=parse_float_list, default=list(DEFAULT_R_VALUES))
    parser.add_argument("--eps-values", type=parse_float_list, default=list(DEFAULT_EPS_VALUES))
    parser.add_argument("--out-prefix", default=DEFAULT_OUT_PREFIX)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = audit_asymptotic_grid(args.R_values, args.eps_values)
    csv_path, json_path, md_path = write_outputs(rows, args.out_prefix, args.R_values, args.eps_values)
    print(f"wrote {csv_path}")
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")
    print(f"grid_points={len(rows)}")


if __name__ == "__main__":
    main()
