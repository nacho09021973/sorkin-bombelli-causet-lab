#!/usr/bin/env python3
"""Numerical audit for the horizon-crossing direct plunging branch.

This is not a proof of the exterior-to-interior fastest-geodesic lemma.  It is
an explicit reproducibility check for the branch used by
``--horizon-shooting``: for horizon-crossing links in a small seed sweep, the
script verifies that the selected direct branch is a regular plunging null,
has a locally unique shooting parameter, and is numerically stable under
Simpson refinement.
"""

from __future__ import annotations

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

import run_schwarzschild_horizon_benchmark as horizon  # noqa: E402
from run_schwarzschild_minimal_benchmark import (  # noqa: E402
    CUBIC_ROOT_EPS,
    MASS,
    direct_phi_integral,
    direct_time_integral,
    find_direct_shooting_c2,
    positive_real_roots_cubic,
)


OUT_DIR = Path(__file__).resolve().parent
DEFAULT_SEED_START = 1
DEFAULT_SEED_STOP = 40
DEFAULT_OUT_PREFIX = "schwarzschild_horizon_shooting_branch_audit"

CSV_FIELDS = (
    "seed",
    "i",
    "j",
    "r_i",
    "r_j",
    "phi_target",
    "c2",
    "c2_over_critical",
    "positive_root_count",
    "phi_value",
    "phi_error",
    "event_dt",
    "null_dt_512",
    "null_dt_1024",
    "null_dt_2048",
    "dt_refine_abs_1024_2048",
    "dt_refine_rel_1024_2048",
    "horizon_integrand",
    "exterior_regular_minus_raw_abs",
    "phi_lower_c2",
    "phi_upper_c2",
    "local_phi_monotone",
    "related_margin",
)


def _regular_integrand_at_horizon(c2: float, mass: float = MASS) -> float:
    """Closed-form value of the rationalized t-tilde integrand at r=2M."""

    return 4.0 * mass * mass + 1.0 / (2.0 * c2)


def _regular_integrand(u: float, c2: float, mass: float = MASS) -> float:
    c = math.sqrt(c2)
    f = horizon.cubic_f(u, c2, mass)
    sf = math.sqrt(f)
    num = c2 * (1.0 + 2.0 * mass * u) + 4.0 * mass * mass * u**4
    den = u * u * sf * (c + 2.0 * mass * u * sf)
    return num / den


def _exterior_rationalization_error(u1: float, c2: float, mass: float = MASS) -> float | None:
    """Compare regular and raw He-Rideout t integrals on an exterior-only segment."""

    horizon_u = 1.0 / (2.0 * mass)
    if u1 >= horizon_u:
        return None
    u_mid = 0.5 * (u1 + horizon_u)
    raw = direct_time_integral(u1, u_mid, c2, mass, intervals=1024)
    regular = horizon._direct_ttilde_integral_horizon_regular(u1, u_mid, c2, mass, intervals=1024)
    if raw is None or regular is None:
        return None
    return abs(raw - regular)


def _local_phi_monotonicity(
    u1: float,
    u2: float,
    c2: float,
    phi_value: float,
    mass: float = MASS,
) -> tuple[float | None, float | None, bool]:
    critical = 1.0 / (27.0 * mass * mass)
    lower_room = max(0.0, (c2 - critical * (1.0 + 1.0e-9)) / c2)
    delta = min(0.01, 0.5 * lower_room)
    if delta <= 0.0:
        return None, None, False

    phi_lower = direct_phi_integral(u1, u2, c2 * (1.0 - delta), mass, intervals=1024)
    phi_upper = direct_phi_integral(u1, u2, c2 * (1.0 + delta), mass, intervals=1024)
    if phi_lower is None or phi_upper is None:
        return phi_lower, phi_upper, False
    return phi_lower, phi_upper, phi_lower >= phi_value >= phi_upper


def audit_crossing_link(seed: int, i: int, j: int, events: list[horizon.Event]) -> dict[str, Any]:
    p = events[i]
    q = events[j]
    u1 = 1.0 / p.r
    u2 = 1.0 / q.r
    phi_target = horizon.angular_separation(p, q)
    found = find_direct_shooting_c2(u1, u2, phi_target, MASS)
    if found is None:
        raise AssertionError(f"seed={seed} link {(i, j)} has no direct shooting c2")

    c2, phi_value = found
    critical = 1.0 / (27.0 * MASS * MASS)
    roots = positive_real_roots_cubic(c2, mass=MASS)
    null_dt_512 = horizon._direct_ttilde_integral_horizon_regular(u1, u2, c2, MASS, intervals=512)
    null_dt_1024 = horizon._direct_ttilde_integral_horizon_regular(u1, u2, c2, MASS, intervals=1024)
    null_dt_2048 = horizon._direct_ttilde_integral_horizon_regular(u1, u2, c2, MASS, intervals=2048)
    if null_dt_512 is None or null_dt_1024 is None or null_dt_2048 is None:
        raise AssertionError(f"seed={seed} link {(i, j)} has non-finite t-tilde integral")

    exterior_error = _exterior_rationalization_error(u1, c2, MASS)
    phi_lower, phi_upper, local_monotone = _local_phi_monotonicity(u1, u2, c2, phi_value, MASS)
    event_dt = q.t - p.t
    refine_abs = abs(null_dt_2048 - null_dt_1024)
    refine_rel = refine_abs / max(abs(null_dt_2048), CUBIC_ROOT_EPS)

    return {
        "seed": seed,
        "i": i,
        "j": j,
        "r_i": p.r,
        "r_j": q.r,
        "phi_target": phi_target,
        "c2": c2,
        "c2_over_critical": c2 / critical,
        "positive_root_count": len(roots),
        "phi_value": phi_value,
        "phi_error": abs(phi_value - phi_target),
        "event_dt": event_dt,
        "null_dt_512": null_dt_512,
        "null_dt_1024": null_dt_1024,
        "null_dt_2048": null_dt_2048,
        "dt_refine_abs_1024_2048": refine_abs,
        "dt_refine_rel_1024_2048": refine_rel,
        "horizon_integrand": _regular_integrand_at_horizon(c2, MASS),
        "regular_integrand_sample": _regular_integrand(0.5, c2, MASS),
        "exterior_regular_minus_raw_abs": exterior_error,
        "phi_lower_c2": phi_lower,
        "phi_upper_c2": phi_upper,
        "local_phi_monotone": local_monotone,
        "related_margin": event_dt - null_dt_2048,
    }


def collect_audit_rows(seed_start: int, seed_stop: int) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for seed in range(seed_start, seed_stop + 1):
        events, matrix, _states, summary = horizon.run_horizon_case(
            8,
            4,
            seed,
            enable_horizon_shooting=True,
        )
        if not summary["antisymmetric"] or not summary["transitive"]:
            raise AssertionError(f"order check failed for seed={seed}")
        for i, j in horizon.transitive_reduction_links(matrix):
            if (
                horizon.region(events[i].r) == horizon.EXTERIOR
                and horizon.region(events[j].r) == horizon.INTERIOR
                and horizon.angular_separation(events[i], events[j]) > horizon.ANGLE_EPS
            ):
                rows.append(audit_crossing_link(seed, i, j, events))
    return rows


def _jsonable(value: Any) -> Any:
    if isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            return None
    return value


def write_outputs(rows: list[dict[str, Any]], out_prefix: str) -> tuple[Path, Path, Path]:
    csv_path = OUT_DIR / f"{out_prefix}.csv"
    json_path = OUT_DIR / f"{out_prefix}.json"
    md_path = OUT_DIR / f"{out_prefix}.md"

    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in CSV_FIELDS})

    max_phi_error = max((row["phi_error"] for row in rows), default=0.0)
    max_refine_rel = max((row["dt_refine_rel_1024_2048"] for row in rows), default=0.0)
    max_exterior_error = max(
        (
            row["exterior_regular_minus_raw_abs"]
            for row in rows
            if row["exterior_regular_minus_raw_abs"] is not None
        ),
        default=0.0,
    )
    all_checks_pass = all(
        row["c2_over_critical"] > 1.0
        and row["positive_root_count"] == 0
        and row["phi_error"] <= 1.0e-8
        and row["dt_refine_rel_1024_2048"] <= 1.0e-8
        and row["local_phi_monotone"]
        and row["related_margin"] >= -horizon.TIME_EPS
        and row["exterior_regular_minus_raw_abs"] is not None
        and row["exterior_regular_minus_raw_abs"] <= 1.0e-12
        for row in rows
    )

    summary = {
        "audit": "S4 Schwarzschild horizon-crossing direct plunging branch",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "seed_start": DEFAULT_SEED_START,
        "seed_stop": DEFAULT_SEED_STOP,
        "audited_crossing_links": len(rows),
        "seeds_with_crossing_links": sorted({row["seed"] for row in rows}),
        "all_checks_pass": all_checks_pass,
        "max_phi_error": max_phi_error,
        "max_dt_refine_rel_1024_2048": max_refine_rel,
        "max_exterior_regular_minus_raw_abs": max_exterior_error,
        "scope": (
            "Numerical branch audit only.  It checks regularity, uniqueness, "
            "and stability of the direct plunging branch used by --horizon-shooting; "
            "it is not a mathematical proof that no non-direct crossing geodesic "
            "can arrive earlier."
        ),
        "rows": [{key: _jsonable(value) for key, value in row.items()} for row in rows],
    }
    json_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    md_lines = [
        "# Schwarzschild Horizon-Shooting Branch Audit",
        "",
        "This is a numerical audit of the direct plunging null branch used by `--horizon-shooting`.",
        "It is evidence for the implementation branch, not a proof of the fastest-geodesic lemma.",
        "",
        f"- Seed range: {DEFAULT_SEED_START}..{DEFAULT_SEED_STOP}",
        f"- Audited horizon-crossing links: {len(rows)}",
        f"- Seeds with crossing links: {summary['seeds_with_crossing_links']}",
        f"- All checks pass: {all_checks_pass}",
        f"- Max angular shooting error: {max_phi_error:.3e}",
        f"- Max Simpson relative drift, 1024 vs 2048: {max_refine_rel:.3e}",
        f"- Max exterior raw-vs-regular time-integral error: {max_exterior_error:.3e}",
        "",
        "Checks per audited link:",
        "",
        "- `c2 > 1/(27M^2)`, so the cubic has no positive root obstruction.",
        "- The rationalized EF-time integrand is finite at the horizon.",
        "- The direct shooting solution reproduces the target angular separation.",
        "- Simpson refinement from 1024 to 2048 intervals is stable.",
        "- Local `phi(c2)` is monotone around the selected shooting parameter.",
        "- The selected null arrives no later than the target event time.",
        "",
        "Residual caveat: this does not prove that every possible non-direct or turning branch is slower.",
    ]
    md_path.write_text("\n".join(md_lines) + "\n", encoding="utf-8")
    return csv_path, json_path, md_path


def main() -> None:
    rows = collect_audit_rows(DEFAULT_SEED_START, DEFAULT_SEED_STOP)
    csv_path, json_path, md_path = write_outputs(rows, DEFAULT_OUT_PREFIX)
    print(f"audited_crossing_links={len(rows)}")
    print(f"wrote {csv_path.name}, {json_path.name}, {md_path.name}")


if __name__ == "__main__":
    main()
