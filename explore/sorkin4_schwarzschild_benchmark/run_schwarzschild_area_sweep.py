#!/usr/bin/env python3
"""Small SORKIN-4 Schwarzschild horizon-area sweep.

This is an aligned radial-strand diagnostic, not a 4D area-law measurement.
It varies Schwarzschild mass, counts exterior-to-interior links in the
transitive reduction, and checks whether the mean count grows with the formal
horizon area A = 16*pi*M^2 under the current IEF radial criterion.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from statistics import mean, pstdev
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_schwarzschild_horizon_benchmark import run_horizon_case  # noqa: E402


OUT_DIR = Path(__file__).resolve().parent
DEFAULT_OUT_PREFIX = "schwarzschild_horizon_area_sweep"
DEFAULT_N_EXTERIOR = 16
DEFAULT_N_INTERIOR = 8
DEFAULT_SEED_START = 1
DEFAULT_SEED_STOP = 40
DEFAULT_MASSES = (0.75, 1.0, 1.25, 1.5, 1.75, 2.0)

CSV_FIELDS = (
    "mass",
    "horizon_radius",
    "formal_area_16pi_M2",
    "seed_start",
    "seed_stop",
    "num_seeds",
    "N_exterior",
    "N_interior",
    "mean_horizon_crossing_links",
    "std_horizon_crossing_links",
    "min_horizon_crossing_links",
    "max_horizon_crossing_links",
    "nonzero_seed_count",
    "mean_total_links",
    "mean_solver_coverage",
    "failed_order_check_count",
)


def parse_float_list(raw: str) -> list[float]:
    values = [float(part.strip()) for part in raw.split(",") if part.strip()]
    if not values:
        raise argparse.ArgumentTypeError("expected at least one comma-separated float")
    return values


def _jsonable(value: Any) -> Any:
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return None
    return value


def run_area_sweep(
    masses: list[float] | tuple[float, ...] = DEFAULT_MASSES,
    seed_start: int = DEFAULT_SEED_START,
    seed_stop: int = DEFAULT_SEED_STOP,
    n_exterior: int = DEFAULT_N_EXTERIOR,
    n_interior: int = DEFAULT_N_INTERIOR,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    aggregate_rows: list[dict[str, Any]] = []
    per_seed_rows: list[dict[str, Any]] = []

    seeds = list(range(seed_start, seed_stop + 1))
    for mass in masses:
        link_counts: list[int] = []
        total_links: list[int] = []
        solver_coverages: list[float] = []
        failed_checks = 0
        for seed in seeds:
            _events, _matrix, _states, summary = run_horizon_case(
                n_exterior=n_exterior,
                n_interior=n_interior,
                seed=seed,
                mass=mass,
                aligned=True,
            )
            horizon_links = int(summary["horizon_crossing_links"])
            links = int(summary["links"])
            coverage = float(summary["solver_coverage"])
            ok = bool(summary["antisymmetric"]) and bool(summary["transitive"])
            if not ok:
                failed_checks += 1
            link_counts.append(horizon_links)
            total_links.append(links)
            solver_coverages.append(coverage)
            per_seed_rows.append(
                {
                    "mass": mass,
                    "horizon_radius": 2.0 * mass,
                    "formal_area_16pi_M2": 16.0 * math.pi * mass * mass,
                    "seed": seed,
                    "N_exterior": n_exterior,
                    "N_interior": n_interior,
                    "horizon_crossing_links": horizon_links,
                    "total_links": links,
                    "solver_coverage": coverage,
                    "antisymmetric": bool(summary["antisymmetric"]),
                    "transitive": bool(summary["transitive"]),
                }
            )

        aggregate_rows.append(
            {
                "mass": mass,
                "horizon_radius": 2.0 * mass,
                "formal_area_16pi_M2": 16.0 * math.pi * mass * mass,
                "seed_start": seed_start,
                "seed_stop": seed_stop,
                "num_seeds": len(seeds),
                "N_exterior": n_exterior,
                "N_interior": n_interior,
                "mean_horizon_crossing_links": mean(link_counts),
                "std_horizon_crossing_links": pstdev(link_counts),
                "min_horizon_crossing_links": min(link_counts),
                "max_horizon_crossing_links": max(link_counts),
                "nonzero_seed_count": sum(1 for value in link_counts if value > 0),
                "mean_total_links": mean(total_links),
                "mean_solver_coverage": mean(solver_coverages),
                "failed_order_check_count": failed_checks,
            }
        )
    return aggregate_rows, per_seed_rows


def write_outputs(
    aggregate_rows: list[dict[str, Any]],
    per_seed_rows: list[dict[str, Any]],
    out_prefix: str,
) -> tuple[Path, Path, Path, Path]:
    csv_path = OUT_DIR / f"{out_prefix}.csv"
    per_seed_csv_path = OUT_DIR / f"{out_prefix}_per_seed.csv"
    json_path = OUT_DIR / f"{out_prefix}.json"
    md_path = OUT_DIR / f"{out_prefix}.md"

    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for row in aggregate_rows:
            writer.writerow({field: row.get(field, "") for field in CSV_FIELDS})

    per_seed_fields = list(per_seed_rows[0]) if per_seed_rows else []
    with per_seed_csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=per_seed_fields)
        writer.writeheader()
        for row in per_seed_rows:
            writer.writerow(row)

    means = [row["mean_horizon_crossing_links"] for row in aggregate_rows]
    monotone_non_decreasing = all(left <= right for left, right in zip(means, means[1:]))
    failed_checks = sum(int(row["failed_order_check_count"]) for row in aggregate_rows)
    summary = {
        "audit": "S4 Schwarzschild aligned radial horizon-area sweep",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "num_mass_points": len(aggregate_rows),
        "mass_values": [row["mass"] for row in aggregate_rows],
        "seed_start": aggregate_rows[0]["seed_start"] if aggregate_rows else None,
        "seed_stop": aggregate_rows[0]["seed_stop"] if aggregate_rows else None,
        "N_exterior": aggregate_rows[0]["N_exterior"] if aggregate_rows else None,
        "N_interior": aggregate_rows[0]["N_interior"] if aggregate_rows else None,
        "monotone_non_decreasing_mean_horizon_links": monotone_non_decreasing,
        "failed_order_check_count": failed_checks,
        "scope": (
            "Aligned radial-strand diagnostic only. This varies M and the formal "
            "Schwarzschild area A=16*pi*M^2, but it is not a 4D Dou-Sorkin area-law "
            "measurement and not a proof of horizon-molecule scaling."
        ),
    }
    json_path.write_text(
        json.dumps(
            {
                "summary": summary,
                "aggregate_rows": [
                    {key: _jsonable(value) for key, value in row.items()}
                    for row in aggregate_rows
                ],
                "per_seed_rows": [
                    {key: _jsonable(value) for key, value in row.items()}
                    for row in per_seed_rows
                ],
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    md_lines = [
        "# Schwarzschild Horizon Area Sweep",
        "",
        "This is an aligned radial-strand SORKIN-4 diagnostic.",
        "It is not a 4D area-law measurement.",
        "",
        f"- Mass values: {summary['mass_values']}",
        f"- Seeds: {summary['seed_start']}..{summary['seed_stop']}",
        f"- N exterior/interior: {summary['N_exterior']}/{summary['N_interior']}",
        f"- Mean horizon-link counts monotone non-decreasing: {monotone_non_decreasing}",
        f"- Failed order checks: {failed_checks}",
        "",
        "| M | A=16piM^2 | mean horizon links | std | min | max | nonzero seeds |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in aggregate_rows:
        md_lines.append(
            "| {mass:.2f} | {area:.6g} | {mean_links:.6g} | {std:.6g} | {min_links} | {max_links} | {nonzero} |".format(
                mass=row["mass"],
                area=row["formal_area_16pi_M2"],
                mean_links=row["mean_horizon_crossing_links"],
                std=row["std_horizon_crossing_links"],
                min_links=row["min_horizon_crossing_links"],
                max_links=row["max_horizon_crossing_links"],
                nonzero=row["nonzero_seed_count"],
            )
        )
    md_lines.extend(
        [
            "",
            "Interpretation:",
            "",
            "- `horizon_crossing_links` counts exterior-to-interior links in the transitive reduction.",
            "- The aligned radial setup gives a controlled 1+1-like diagnostic with no non-radial undecided pairs.",
            "- The result checks monotonic response to increasing formal horizon area, not the 4D area coefficient.",
        ]
    )
    md_path.write_text("\n".join(md_lines) + "\n", encoding="utf-8")
    return csv_path, per_seed_csv_path, json_path, md_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run a small aligned Schwarzschild horizon-area sweep.")
    parser.add_argument("--masses", type=parse_float_list, default=list(DEFAULT_MASSES))
    parser.add_argument("--seed-start", type=int, default=DEFAULT_SEED_START)
    parser.add_argument("--seed-stop", type=int, default=DEFAULT_SEED_STOP)
    parser.add_argument("--n-exterior", type=int, default=DEFAULT_N_EXTERIOR)
    parser.add_argument("--n-interior", type=int, default=DEFAULT_N_INTERIOR)
    parser.add_argument("--out-prefix", default=DEFAULT_OUT_PREFIX)
    args = parser.parse_args()
    if args.seed_stop < args.seed_start:
        raise SystemExit("--seed-stop must be >= --seed-start")
    if args.n_exterior <= 0 or args.n_interior <= 0:
        raise SystemExit("--n-exterior and --n-interior must be positive")
    return args


def main() -> None:
    args = parse_args()
    aggregate_rows, per_seed_rows = run_area_sweep(
        masses=args.masses,
        seed_start=args.seed_start,
        seed_stop=args.seed_stop,
        n_exterior=args.n_exterior,
        n_interior=args.n_interior,
    )
    paths = write_outputs(aggregate_rows, per_seed_rows, args.out_prefix)
    for path in paths:
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
