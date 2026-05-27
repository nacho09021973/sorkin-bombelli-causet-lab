#!/usr/bin/env python3
"""Numerical audit for the exterior outgoing turning-branch competitor.

This is a falsifiable diagnostic, not a proof of the fastest-geodesic lemma.
It checks exterior outgoing pairs for which the implemented direct shooting
branch exists, then asks whether a one-turn competitor can hit the same angular
separation and, when it can, whether its EF arrival time is later.

The only competitor considered here is:

    u1 -> u_turn -> u2, with u2 < u1 < u_turn < 1/(3M)

so both endpoints are outside the photon sphere.  Ingoing pairs are skipped:
the turning branch cannot reach a larger u after turning because u_turn is the
maximum u on that geodesic.
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

from run_schwarzschild_minimal_benchmark import (  # noqa: E402
    ANGLE_EPS,
    CUBIC_ROOT_EPS,
    MASS,
    SHOOTING_PHI_TOL,
    TIME_EPS,
    angular_separation,
    cubic_f,
    direct_time_integral,
    find_direct_shooting_c2,
    generate_exterior_events,
    positive_real_roots_cubic,
)


OUT_DIR = Path(__file__).resolve().parent
DEFAULT_N = 12
DEFAULT_SEED_START = 1959
DEFAULT_SEED_STOP = 1968
DEFAULT_OUT_PREFIX = "schwarzschild_exterior_turning_branch_audit"
SIMPSON_INTERVALS = 2048
TURNING_PHI_TOL = 1.0e-6

CSV_FIELDS = (
    "seed",
    "i",
    "j",
    "r_i",
    "r_j",
    "u1",
    "u2",
    "phi_target",
    "direct_c2",
    "direct_phi",
    "direct_dt",
    "direct_phi_max",
    "turning_phi_min",
    "turning_range_gap",
    "turning_c2",
    "turning_u",
    "turning_phi",
    "turning_dt",
    "turning_minus_direct_dt",
    "turning_c2_over_direct_c2",
    "outcome",
)


def _jsonable(value: Any) -> Any:
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return None
    return value


def _simpson_integral_unit(func, intervals: int = SIMPSON_INTERVALS) -> float | None:
    if intervals <= 0:
        raise ValueError("intervals must be positive")
    if intervals % 2:
        intervals += 1
    h = 1.0 / intervals
    total = 0.0
    for idx in range(intervals + 1):
        x = idx * h
        value = func(x)
        if not math.isfinite(value):
            return None
        coeff = 1 if idx == 0 or idx == intervals else 4 if idx % 2 else 2
        total += coeff * value
    result = total * h / 3.0
    return result if math.isfinite(result) else None


def _lower_turning_root(c2: float, mass: float = MASS) -> float | None:
    roots = positive_real_roots_cubic(c2, mass=mass)
    if len(roots) < 2:
        return None
    return roots[0]


def _integral_to_turn(
    u_start: float,
    u_turn: float,
    c2: float,
    regular_sign: float | None,
    mass: float = MASS,
    intervals: int = SIMPSON_INTERVALS,
) -> float | None:
    """Integrate from u_start to the lower root with endpoint regularisation.

    If regular_sign is None, integrate only dphi = du/sqrt(f).  Otherwise
    integrate the EF-time segment after reversing any outgoing leg so the
    singular term is positive.  regular_sign=-1 is the inward leg
    u_start -> u_turn; regular_sign=+1 is the reversed outward leg
    u_start -> u_turn corresponding to u_turn -> u_start.
    """

    if not u_start < u_turn:
        return None
    length = u_turn - u_start
    c = math.sqrt(c2)
    fprime_turn = 6.0 * mass * u_turn * u_turn - 2.0 * u_turn
    if fprime_turn >= -CUBIC_ROOT_EPS:
        return None
    endpoint_singular_limit = 2.0 * math.sqrt(length / abs(fprime_turn))
    denom_turn = u_turn * u_turn - 2.0 * mass * u_turn**3
    if regular_sign is not None and abs(denom_turn) <= CUBIC_ROOT_EPS:
        return None

    def transformed(x: float) -> float:
        one_minus_x = 1.0 - x
        u = u_start + length * (2.0 * x - x * x)
        du_dx = 2.0 * length * one_minus_x

        if abs(one_minus_x) <= 1.0e-14:
            if regular_sign is None:
                return endpoint_singular_limit
            return c * endpoint_singular_limit / denom_turn

        f_value = cubic_f(u, c2, mass)
        if f_value <= 0.0:
            return math.nan
        singular = du_dx / math.sqrt(f_value)
        if regular_sign is None:
            return singular

        denom = u * u - 2.0 * mass * u**3
        if abs(denom) <= CUBIC_ROOT_EPS:
            return math.nan
        regular = regular_sign * 2.0 * mass * u * du_dx
        return (c * singular + regular) / denom

    return _simpson_integral_unit(transformed, intervals)


def turning_phi_integral(
    u1: float,
    u2: float,
    c2: float,
    mass: float = MASS,
    intervals: int = SIMPSON_INTERVALS,
) -> tuple[float, float] | None:
    u_turn = _lower_turning_root(c2, mass)
    if u_turn is None or not (u2 < u1 < u_turn < 1.0 / (3.0 * mass)):
        return None
    first = _integral_to_turn(u1, u_turn, c2, None, mass, intervals)
    second = _integral_to_turn(u2, u_turn, c2, None, mass, intervals)
    if first is None or second is None:
        return None
    return first + second, u_turn


def turning_time_integral(
    u1: float,
    u2: float,
    c2: float,
    mass: float = MASS,
    intervals: int = SIMPSON_INTERVALS,
) -> tuple[float, float] | None:
    u_turn = _lower_turning_root(c2, mass)
    if u_turn is None or not (u2 < u1 < u_turn < 1.0 / (3.0 * mass)):
        return None
    inward = _integral_to_turn(u1, u_turn, c2, -1.0, mass, intervals)
    outward_reversed = _integral_to_turn(u2, u_turn, c2, +1.0, mass, intervals)
    if inward is None or outward_reversed is None:
        return None
    return inward + outward_reversed, u_turn


def find_turning_c2(
    u1: float,
    u2: float,
    phi_target: float,
    mass: float = MASS,
) -> tuple[float, float, float] | None:
    """Find the c2<critical one-turn branch hitting phi_target."""

    if not (u2 < u1 < 1.0 / (3.0 * mass)) or phi_target <= ANGLE_EPS:
        return None

    critical = 1.0 / (27.0 * mass * mass)
    # The lower root equals u1 when c2 = u1^2 - 2M u1^3.  The turn exists
    # only above this value and below the critical photon-sphere value.
    low = (u1 * u1 - 2.0 * mass * u1**3) * (1.0 + 1.0e-8)
    high = critical * (1.0 - 1.0e-10)
    if not (0.0 < low < high):
        return None

    phi_low_turn = turning_phi_integral(u1, u2, low, mass)
    for _ in range(12):
        if phi_low_turn is not None:
            break
        low = low + 0.25 * (high - low)
        phi_low_turn = turning_phi_integral(u1, u2, low, mass)
    phi_high_turn = turning_phi_integral(u1, u2, high, mass)
    if phi_low_turn is None or phi_high_turn is None:
        return None
    phi_low, _ = phi_low_turn
    phi_high, _ = phi_high_turn

    if not (phi_low <= phi_target <= phi_high):
        return None

    best_c2 = low
    best_phi = phi_low
    best_turn = phi_low_turn[1]
    for _ in range(100):
        mid = 0.5 * (low + high)
        value = turning_phi_integral(u1, u2, mid, mass)
        if value is None:
            low = mid
            continue
        phi_mid, u_turn = value
        best_c2 = mid
        best_phi = phi_mid
        best_turn = u_turn
        if abs(phi_mid - phi_target) <= TURNING_PHI_TOL:
            return best_c2, best_phi, best_turn
        if phi_mid < phi_target:
            low = mid
        else:
            high = mid

    if abs(best_phi - phi_target) <= 10.0 * TURNING_PHI_TOL:
        return best_c2, best_phi, best_turn
    return None


def branch_angular_ranges(
    u1: float,
    u2: float,
    mass: float = MASS,
) -> tuple[float | None, float | None, float | None]:
    """Return direct max, turning min, and turning_min-direct_max gap."""

    critical = 1.0 / (27.0 * mass * mass)
    direct_phi_max = None
    if u2 < u1:
        from run_schwarzschild_minimal_benchmark import direct_phi_integral

        direct_phi_max = direct_phi_integral(
            u1,
            u2,
            critical * (1.0 + 1.0e-9),
            mass,
            intervals=SIMPSON_INTERVALS,
        )

    turning_phi_min = None
    if u2 < u1 < 1.0 / (3.0 * mass):
        # Estimate the lower angular endpoint by placing the turning root
        # just above u1.  This is more stable than stepping in c2, especially
        # in the weak-field nearby-radius corner where u1-u2 is tiny.
        room = 1.0 / (3.0 * mass) - u1
        radial_gap = u1 - u2
        for scale in (
            1.0e-8,
            3.0e-8,
            1.0e-7,
            3.0e-7,
            1.0e-6,
            3.0e-6,
            1.0e-5,
            3.0e-5,
            1.0e-4,
            3.0e-4,
            1.0e-3,
            3.0e-3,
            1.0e-2,
        ):
            delta = max(scale * radial_gap, 1.0e-14 * max(1.0, u1))
            if delta >= room:
                continue
            u_turn = u1 + delta
            c2 = u_turn * u_turn - 2.0 * mass * u_turn**3
            if not (0.0 < c2 < critical):
                continue
            value = turning_phi_integral(u1, u2, c2, mass)
            if value is not None:
                turning_phi_min = value[0]
                break

    gap = (
        turning_phi_min - direct_phi_max
        if turning_phi_min is not None and direct_phi_max is not None
        else None
    )
    return direct_phi_max, turning_phi_min, gap


def audit_pair(seed: int, i: int, j: int, p, q) -> dict[str, Any] | None:
    """Return an audit row for an outgoing exterior pair, or None if skipped."""

    if q.r <= p.r:
        return None
    u1 = 1.0 / p.r
    u2 = 1.0 / q.r
    phi_target = angular_separation(p, q)
    if phi_target <= ANGLE_EPS:
        return None

    direct = find_direct_shooting_c2(u1, u2, phi_target, MASS)
    if direct is None:
        return None
    direct_c2, direct_phi = direct
    direct_dt = direct_time_integral(u1, u2, direct_c2, MASS, intervals=SIMPSON_INTERVALS)
    if direct_dt is None or not math.isfinite(direct_dt):
        return None
    direct_phi_max, turning_phi_min, turning_range_gap = branch_angular_ranges(u1, u2, MASS)

    turning = find_turning_c2(u1, u2, phi_target, MASS)
    if turning is None:
        return {
            "seed": seed,
            "i": i,
            "j": j,
            "r_i": p.r,
            "r_j": q.r,
            "u1": u1,
            "u2": u2,
            "phi_target": phi_target,
            "direct_c2": direct_c2,
            "direct_phi": direct_phi,
            "direct_dt": direct_dt,
            "direct_phi_max": direct_phi_max,
            "turning_phi_min": turning_phi_min,
            "turning_range_gap": turning_range_gap,
            "turning_c2": None,
            "turning_u": None,
            "turning_phi": None,
            "turning_dt": None,
            "turning_minus_direct_dt": None,
            "turning_c2_over_direct_c2": None,
            "outcome": "no_turning_solution",
        }

    turning_c2, turning_phi, turning_u = turning
    turning_time = turning_time_integral(u1, u2, turning_c2, MASS, intervals=SIMPSON_INTERVALS)
    if turning_time is None:
        outcome = "turning_time_failed"
        turning_dt = None
        margin = None
    else:
        turning_dt, turning_u_time = turning_time
        if abs(turning_u_time - turning_u) > 1.0e-8:
            outcome = "turning_root_mismatch"
        else:
            margin = turning_dt - direct_dt
            outcome = "turning_later" if margin > TIME_EPS else "turning_not_later"

    return {
        "seed": seed,
        "i": i,
        "j": j,
        "r_i": p.r,
        "r_j": q.r,
        "u1": u1,
        "u2": u2,
        "phi_target": phi_target,
        "direct_c2": direct_c2,
        "direct_phi": direct_phi,
        "direct_dt": direct_dt,
        "direct_phi_max": direct_phi_max,
        "turning_phi_min": turning_phi_min,
        "turning_range_gap": turning_range_gap,
        "turning_c2": turning_c2,
        "turning_u": turning_u,
        "turning_phi": turning_phi,
        "turning_dt": turning_dt,
        "turning_minus_direct_dt": margin,
        "turning_c2_over_direct_c2": turning_c2 / direct_c2,
        "outcome": outcome,
    }


def collect_rows(n: int, seed_start: int, seed_stop: int) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for seed in range(seed_start, seed_stop + 1):
        events = generate_exterior_events(n, seed)
        for i in range(n - 1):
            for j in range(i + 1, n):
                row = audit_pair(seed, i, j, events[i], events[j])
                if row is not None:
                    rows.append(row)
    return rows


def write_outputs(rows: list[dict[str, Any]], out_prefix: str, n: int, seed_start: int, seed_stop: int) -> tuple[Path, Path, Path]:
    csv_path = OUT_DIR / f"{out_prefix}.csv"
    json_path = OUT_DIR / f"{out_prefix}.json"
    md_path = OUT_DIR / f"{out_prefix}.md"

    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in CSV_FIELDS})

    outcome_counts = {outcome: sum(1 for row in rows if row["outcome"] == outcome) for outcome in sorted({row["outcome"] for row in rows})}
    margins = [
        row["turning_minus_direct_dt"]
        for row in rows
        if row["turning_minus_direct_dt"] is not None
    ]
    range_gaps = [
        row["turning_range_gap"]
        for row in rows
        if row["turning_range_gap"] is not None
    ]
    checked_turning_branches = len(margins)
    all_checked_later = (
        all(row["outcome"] != "turning_not_later" for row in rows)
        if checked_turning_branches
        else None
    )
    summary = {
        "audit": "S4 Schwarzschild exterior outgoing turning-branch competitor",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "N": n,
        "seed_start": seed_start,
        "seed_stop": seed_stop,
        "audited_outgoing_direct_pairs": len(rows),
        "checked_turning_branches": checked_turning_branches,
        "outcome_counts": outcome_counts,
        "all_checked_turning_branches_later": all_checked_later,
        "min_turning_minus_direct_dt": min(margins) if margins else None,
        "max_turning_minus_direct_dt": max(margins) if margins else None,
        "min_turning_range_gap": min(range_gaps) if range_gaps else None,
        "max_turning_range_gap": max(range_gaps) if range_gaps else None,
        "scope": (
            "Numerical competitor audit only. Ingoing pairs are skipped because "
            "a one-turn branch cannot reach larger u after turning. Outgoing "
            "turning competitors are checked only when both endpoints are outside "
            "the photon sphere and the one-turn branch can hit the same angular "
            "separation. This is not a proof of the fastest-geodesic lemma."
        ),
        "rows": [{key: _jsonable(value) for key, value in row.items()} for row in rows],
    }
    json_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    md_lines = [
        "# Schwarzschild Exterior Turning-Branch Audit",
        "",
        "This is a numerical audit of the outgoing one-turn competitor to the direct exterior shooting branch.",
        "It is a falsifiable diagnostic, not a proof of the fastest-geodesic lemma.",
        "",
        f"- Seed range: {seed_start}..{seed_stop}",
        f"- N: {n}",
        f"- Audited outgoing direct pairs: {len(rows)}",
        f"- Checked turning branches: {checked_turning_branches}",
        f"- Outcome counts: {outcome_counts}",
        f"- All checked turning branches later: {all_checked_later}",
        f"- Min turning minus direct EF-time: {summary['min_turning_minus_direct_dt']}",
        f"- Min turning angular range gap: {summary['min_turning_range_gap']}",
        "",
        "Scope notes:",
        "",
        "- Ingoing pairs are skipped: the turning point is the maximum u, so a turned branch cannot later reach a larger u.",
        "- The genuine competitor is outgoing and requires u2 < u1 < u_turn < 1/(3M), i.e. both endpoints outside the photon sphere.",
        "- The code uses c2=(E/L)^2. Direct no-root branches are on c2 > 1/(27M^2); one-turn branches are on c2 < 1/(27M^2).",
    ]
    md_path.write_text("\n".join(md_lines) + "\n", encoding="utf-8")
    return csv_path, json_path, md_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit exterior outgoing turning-branch competitors.")
    parser.add_argument("--N", type=int, default=DEFAULT_N)
    parser.add_argument("--seed-start", type=int, default=DEFAULT_SEED_START)
    parser.add_argument("--seed-stop", type=int, default=DEFAULT_SEED_STOP)
    parser.add_argument("--out-prefix", default=DEFAULT_OUT_PREFIX)
    args = parser.parse_args()
    if args.N <= 0:
        raise SystemExit("--N must be positive")
    if args.seed_stop < args.seed_start:
        raise SystemExit("--seed-stop must be >= --seed-start")
    return args


def main() -> None:
    args = parse_args()
    rows = collect_rows(args.N, args.seed_start, args.seed_stop)
    csv_path, json_path, md_path = write_outputs(rows, args.out_prefix, args.N, args.seed_start, args.seed_stop)
    print(f"wrote {csv_path}")
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")
    print(f"audited_outgoing_direct_pairs={len(rows)}")


if __name__ == "__main__":
    main()
