#!/usr/bin/env python3
"""S4-SCHW-STABILITY-001: stability sweep for the exterior Schwarzschild benchmark.

Sweeps N, seed, and margin to answer four diagnostic questions:

1. Do basic order checks (antisymmetric, transitive, decided-transitivity) hold
   across the whole grid?
2. Does the outgoing one-turn turning branch stay outside the direct-branch
   angular range on all sampled pairs?
3. Is the angular-range gap consistent and positive across the parameter grid
   (compatible with the single frozen benchmark at margin=0.35, seed=1959)?
4. Do any True assertions appear suspiciously close to the inner margin or outer
   domain boundary?

This is an exploratory diagnostic, not a claim about embeddability or physics.
A pass on all checks means the benchmark behaves stably across the swept
parameters, not that any physical causal structure has been confirmed.
"""

from __future__ import annotations

import csv
import json
import math
import random
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_schwarzschild_minimal_benchmark import (  # noqa: E402
    Event,
    MASS,
    R_MAX,
    SCHWARZSCHILD_RADIUS,
    T_MIN,
    T_MAX,
    build_causal_matrix,
    check_antisymmetric,
    check_decided_transitivity,
    check_transitive,
)
from audit_exterior_turning_branch import branch_angular_ranges  # noqa: E402


# ---------------------------------------------------------------------------
# Audit identity
# ---------------------------------------------------------------------------

AUDIT_ID = "S4-SCHW-STABILITY-001"
OUT_DIR = Path(__file__).resolve().parent
DEFAULT_OUT_PREFIX = "schwarzschild_stability_001"

# ---------------------------------------------------------------------------
# Sweep parameters (fixed for this audit)
# ---------------------------------------------------------------------------

SWEEP_N = (12, 16)
SWEEP_SEED = (1959, 1960, 1961)
SWEEP_MARGIN = (0.25, 0.35, 0.50)

# ---------------------------------------------------------------------------
# Near-margin diagnostic thresholds
# ---------------------------------------------------------------------------

# A True assertion is flagged as "near-inner" if either event has
#   r < r_min + NEAR_INNER_FRAC * (R_MAX - r_min)
# and "near-outer" if either event has
#   r > R_MAX - NEAR_OUTER_FRAC * (R_MAX - r_min).
# These thresholds are diagnostic only; they do not gate all_checks_pass.
NEAR_INNER_FRAC = 0.10
NEAR_OUTER_FRAC = 0.10

# ---------------------------------------------------------------------------
# CSV schema
# ---------------------------------------------------------------------------

CSV_FIELDS = (
    "N",
    "seed",
    "margin",
    "r_min",
    "true_relations",
    "false_relations",
    "undecided_pairs",
    "decided_pairs",
    "antisymmetric",
    "transitive",
    "decided_transitivity",
    "outgoing_pairs_with_gap",
    "outgoing_gaps_all_positive",
    "min_gap",
    "near_inner_true_count",
    "near_outer_true_count",
    "local_checks_pass",
)


# ---------------------------------------------------------------------------
# Event generation with explicit margin
# ---------------------------------------------------------------------------

def generate_exterior_events_with_margin(
    n: int,
    seed: int,
    margin: float,
) -> list[Event]:
    """Generate n exterior events with r > schwarzschild_radius + margin.

    Uses the same volume-weighted r^2 sampling and coordinate generation as
    run_schwarzschild_minimal_benchmark.generate_exterior_events.  The only
    difference is that r_min is computed from the explicit margin argument
    instead of the module-level EXTERIOR_MARGIN constant, so we can sweep it.
    """
    r_min = SCHWARZSCHILD_RADIUS + margin
    rng = random.Random(seed)
    events: list[Event] = []
    for index in range(n):
        t = rng.uniform(T_MIN, T_MAX)
        u_frac = rng.random()
        r = (r_min ** 3 + u_frac * (R_MAX ** 3 - r_min ** 3)) ** (1.0 / 3.0)
        cos_theta = rng.uniform(-1.0, 1.0)
        theta = math.acos(cos_theta)
        phi = rng.uniform(0.0, 2.0 * math.pi)
        events.append(Event(index=index, t=t, r=r, theta=theta, phi=phi))
    return sorted(events, key=lambda e: (e.t, e.r, e.theta, e.phi))


# ---------------------------------------------------------------------------
# Per-pair turning-branch gap diagnostic
# ---------------------------------------------------------------------------

def _collect_turning_gaps(events: list[Event]) -> list[float]:
    """For each outgoing pair (r_j > r_i), collect the angular-range gap.

    The gap is phi_turning_min - phi_direct_max from branch_angular_ranges.
    Only pairs where branch_angular_ranges returns a valid gap are included.
    A positive gap means the turning branch cannot share an angular target
    with the direct branch for that pair geometry.
    """
    gaps: list[float] = []
    n = len(events)
    for i in range(n - 1):
        for j in range(i + 1, n):
            if events[j].r <= events[i].r:
                # Skip ingoing pairs; turning branch constraint only applies
                # to the outgoing direction.
                continue
            u1 = 1.0 / events[i].r
            u2 = 1.0 / events[j].r
            # u2 < u1 since r_j > r_i, i.e. this is outgoing in u-space.
            _, _, gap = branch_angular_ranges(u1, u2, MASS)
            if gap is not None:
                gaps.append(gap)
    return gaps


# ---------------------------------------------------------------------------
# Per-pair near-margin True-assertion diagnostic
# ---------------------------------------------------------------------------

def _near_margin_true_counts(
    events: list[Event],
    matrix: list[list[bool]],
    r_min: float,
) -> tuple[int, int]:
    """Count True assertions involving events near the inner or outer boundary.

    Returns (near_inner_count, near_outer_count).  An assertion is
    "near-inner" if either endpoint satisfies r < r_min + NEAR_INNER_FRAC *
    domain_width, and "near-outer" if either endpoint satisfies
    r > R_MAX - NEAR_OUTER_FRAC * domain_width.
    These counts are diagnostic; a nonzero value is not an error.
    """
    domain_width = R_MAX - r_min
    near_inner_threshold = r_min + NEAR_INNER_FRAC * domain_width
    near_outer_threshold = R_MAX - NEAR_OUTER_FRAC * domain_width
    n = len(events)
    near_inner_count = 0
    near_outer_count = 0
    for i in range(n - 1):
        for j in range(i + 1, n):
            if not matrix[i][j]:
                continue
            ri, rj = events[i].r, events[j].r
            if ri < near_inner_threshold or rj < near_inner_threshold:
                near_inner_count += 1
            if ri > near_outer_threshold or rj > near_outer_threshold:
                near_outer_count += 1
    return near_inner_count, near_outer_count


# ---------------------------------------------------------------------------
# Single-cell runner
# ---------------------------------------------------------------------------

def run_stability_cell(
    n: int,
    seed: int,
    margin: float,
) -> dict[str, Any]:
    """Run one (N, seed, margin) stability cell and return a result dict."""
    r_min = SCHWARZSCHILD_RADIUS + margin
    events = generate_exterior_events_with_margin(n, seed, margin)
    matrix, states, _ = build_causal_matrix(events, seed=seed, enable_shooting=True)

    n_events = len(events)
    possible = n_events * (n_events - 1) // 2
    true_rel = sum(1 for i in range(n_events - 1) for j in range(i + 1, n_events) if matrix[i][j])
    false_rel = sum(1 for i in range(n_events - 1) for j in range(i + 1, n_events) if states[i][j] is False)
    undecided = sum(1 for i in range(n_events - 1) for j in range(i + 1, n_events) if states[i][j] is None)

    antisym = check_antisymmetric(matrix)
    transitive = check_transitive(matrix)
    decided_trans = check_decided_transitivity(states)

    gaps = _collect_turning_gaps(events)
    gaps_all_positive = all(g > 0.0 for g in gaps) if gaps else True
    min_gap = min(gaps) if gaps else None

    near_inner, near_outer = _near_margin_true_counts(events, matrix, r_min)

    local_pass = antisym and transitive and decided_trans and gaps_all_positive

    return {
        "N": n,
        "seed": seed,
        "margin": margin,
        "r_min": r_min,
        "true_relations": true_rel,
        "false_relations": false_rel,
        "undecided_pairs": undecided,
        "decided_pairs": true_rel + false_rel,
        "antisymmetric": antisym,
        "transitive": transitive,
        "decided_transitivity": decided_trans,
        "outgoing_pairs_with_gap": len(gaps),
        "outgoing_gaps_all_positive": gaps_all_positive,
        "min_gap": min_gap,
        "near_inner_true_count": near_inner,
        "near_outer_true_count": near_outer,
        "local_checks_pass": local_pass,
    }


# ---------------------------------------------------------------------------
# Full sweep
# ---------------------------------------------------------------------------

def run_full_sweep(
    n_values: tuple[int, ...] = SWEEP_N,
    seeds: tuple[int, ...] = SWEEP_SEED,
    margins: tuple[float, ...] = SWEEP_MARGIN,
) -> list[dict[str, Any]]:
    """Run the full stability sweep and return one row per (N, seed, margin)."""
    rows: list[dict[str, Any]] = []
    for n in n_values:
        for seed in seeds:
            for margin in margins:
                rows.append(run_stability_cell(n, seed, margin))
    return rows


# ---------------------------------------------------------------------------
# Aggregate
# ---------------------------------------------------------------------------

def compute_aggregate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Summarise the sweep rows into a single aggregate with all_checks_pass."""
    n_rows = len(rows)
    order_ok = all(r["local_checks_pass"] for r in rows)
    gaps_ok = all(r["outgoing_gaps_all_positive"] for r in rows)
    all_pass = order_ok and gaps_ok

    valid_gaps = [r["min_gap"] for r in rows if r["min_gap"] is not None]
    global_min_gap = min(valid_gaps) if valid_gaps else None
    global_max_gap = max(valid_gaps) if valid_gaps else None

    return {
        "audit": AUDIT_ID,
        "description": (
            "Stability sweep: does the exterior Schwarzschild benchmark stay "
            "well-behaved across seed, N, and margin?"
        ),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "sweep_N": list(SWEEP_N),
        "sweep_seed": list(SWEEP_SEED),
        "sweep_margin": list(SWEEP_MARGIN),
        "n_cells": n_rows,
        "Q1_order_checks_all_pass": order_ok,
        "Q2_turning_branch_no_invasion": gaps_ok,
        "Q3_min_turning_gap_all_cells": global_min_gap,
        "Q3_max_turning_gap_all_cells": global_max_gap,
        "Q3_gap_consistent_positive": (global_min_gap is not None and global_min_gap > 0.0),
        "Q4_near_inner_true_assertions_total": sum(r["near_inner_true_count"] for r in rows),
        "Q4_near_outer_true_assertions_total": sum(r["near_outer_true_count"] for r in rows),
        "all_checks_pass": all_pass,
        "notes": [
            "Q1: antisymmetric + transitive + decided-transitivity must hold for every cell.",
            "Q2: every outgoing pair with a valid turning-branch gap must have gap > 0.",
            "Q3: the gap min across all cells is reported; a positive value is consistent with the frozen single case.",
            "Q4: near-margin True assertion counts are informational only; they do not gate all_checks_pass.",
            "all_checks_pass = Q1 and Q2 only (Q3 positivity check is folded into Q2; Q4 is advisory).",
        ],
    }


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def _fmt(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.12g}"
    return str(value)


def write_outputs(
    rows: list[dict[str, Any]],
    aggregate: dict[str, Any],
    out_prefix: str = DEFAULT_OUT_PREFIX,
) -> tuple[Path, Path, Path]:
    csv_path = OUT_DIR / f"{out_prefix}.csv"
    json_path = OUT_DIR / f"{out_prefix}.json"
    md_path = OUT_DIR / f"{out_prefix}.md"

    # CSV
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({f: _fmt(row[f]) for f in CSV_FIELDS})

    # JSON
    payload = {
        "aggregate": aggregate,
        "rows": rows,
    }
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    # Markdown
    _write_md(md_path, rows, aggregate)

    return csv_path, json_path, md_path


def _write_md(
    md_path: Path,
    rows: list[dict[str, Any]],
    aggregate: dict[str, Any],
) -> None:
    lines: list[str] = [
        f"# {AUDIT_ID}: Exterior Schwarzschild Stability Sweep",
        "",
        f"Generated: {aggregate['generated_at_utc']}",
        "",
        "## Purpose",
        "",
        "Exploratory diagnostic. Sweeps `N`, `seed`, and `margin` to check whether",
        "the exterior Schwarzschild benchmark remains well-behaved outside the single",
        "frozen reference case (N=12, seed=1959, margin=0.35).",
        "",
        "## Sweep Parameters",
        "",
        f"- N: {aggregate['sweep_N']}",
        f"- seed: {aggregate['sweep_seed']}",
        f"- margin: {aggregate['sweep_margin']}",
        f"- Total cells: {aggregate['n_cells']}",
        "",
        "## Aggregate Results",
        "",
        f"| Check | Result |",
        f"|-------|--------|",
        f"| Q1: all order checks pass | {aggregate['Q1_order_checks_all_pass']} |",
        f"| Q2: turning branch no invasion | {aggregate['Q2_turning_branch_no_invasion']} |",
        f"| Q3: min turning gap (all cells) | {aggregate['Q3_min_turning_gap_all_cells']} |",
        f"| Q3: max turning gap (all cells) | {aggregate['Q3_max_turning_gap_all_cells']} |",
        f"| Q3: gap consistent positive | {aggregate['Q3_gap_consistent_positive']} |",
        f"| Q4: near-inner True assertions (total) | {aggregate['Q4_near_inner_true_assertions_total']} |",
        f"| Q4: near-outer True assertions (total) | {aggregate['Q4_near_outer_true_assertions_total']} |",
        f"| **all_checks_pass** | **{aggregate['all_checks_pass']}** |",
        "",
        "## Per-Cell Results",
        "",
        "| N | seed | margin | true | undecided | antisym | transitive | gaps | min_gap | local_pass |",
        "|---|------|--------|------|-----------|---------|------------|------|---------|------------|",
    ]
    for r in rows:
        min_gap_str = f"{r['min_gap']:.4f}" if r["min_gap"] is not None else "—"
        lines.append(
            f"| {r['N']} | {r['seed']} | {r['margin']:.2f} "
            f"| {r['true_relations']} | {r['undecided_pairs']} "
            f"| {r['antisymmetric']} | {r['transitive']} "
            f"| {r['outgoing_pairs_with_gap']} | {min_gap_str} "
            f"| {r['local_checks_pass']} |"
        )

    lines += [
        "",
        "## Notes",
        "",
        "- `gaps` = number of outgoing pairs for which a valid turning-branch gap was computed.",
        "- `min_gap` = minimum of (phi_turning_min − phi_direct_max) over those pairs; a positive",
        "  value means the turning branch cannot reach any target that the direct branch can reach.",
        "- Q4 counts near-margin True assertions using a 10 % domain-width threshold; these are",
        "  informational only and do not affect `all_checks_pass`.",
        "- This diagnostic does not claim embeddability, manifoldlikeness, or physical correctness.",
        "  It only reports algorithmic stability of the partial He-Rideout exterior model across",
        "  the swept parameters.",
        "",
    ]
    md_path.write_text("\n".join(lines), encoding="utf-8")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    print(f"Running {AUDIT_ID}: exterior Schwarzschild stability sweep")
    print(f"  N in {SWEEP_N}, seed in {SWEEP_SEED}, margin in {SWEEP_MARGIN}")
    rows = run_full_sweep()
    aggregate = compute_aggregate(rows)
    csv_path, json_path, md_path = write_outputs(rows, aggregate)
    print(f"wrote {csv_path}")
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")
    print(f"cells={aggregate['n_cells']}")
    print(f"Q1_order_checks_all_pass={aggregate['Q1_order_checks_all_pass']}")
    print(f"Q2_turning_branch_no_invasion={aggregate['Q2_turning_branch_no_invasion']}")
    print(f"Q3_gap_consistent_positive={aggregate['Q3_gap_consistent_positive']}")
    print(f"Q3_min_turning_gap={aggregate['Q3_min_turning_gap_all_cells']}")
    print(f"Q4_near_inner_true_total={aggregate['Q4_near_inner_true_assertions_total']}")
    print(f"Q4_near_outer_true_total={aggregate['Q4_near_outer_true_assertions_total']}")
    print(f"all_checks_pass={aggregate['all_checks_pass']}")


if __name__ == "__main__":
    main()
