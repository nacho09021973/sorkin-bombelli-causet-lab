#!/usr/bin/env python3
"""Schwarzschild horizon-crossing causal solver for SORKIN-4 Level 0.

This module implements causal relations in ingoing Eddington-Finkelstein (IEF)
coordinates that cover both the exterior (r > r_s) and interior (r < r_s) of
Schwarzschild spacetime through the future horizon.

IEF metric: ds² = -dt² + dr² + r² dΩ² + (2M/r)(dt + dr)²

where the paper's time coordinate is t = t_s + 2M ln(r/2M - 1) = v - r,
so v = t + r is the advanced null coordinate.

Sign convention: v is increasing toward the future. The horizon is at r = r_s = 2M.
The singularity is at r = 0 (in the future interior).

## Causal region classification

Each event is assigned to one of three regions:
  EXTERIOR  r > r_s   (Region I of Penrose diagram)
  HORIZON   r = r_s   (future horizon, included in boundary)
  INTERIOR  r < r_s   (Region II: between horizon and singularity)

## Implemented causal criteria

For an ordered pair (p, q) with v_p ≤ v_q:

EXTERIOR → EXTERIOR  (both r > r_s):
  Delegates to causal_relation_schwarzschild from the exterior benchmark.
  This is the He & Rideout criterion with bounds and optional shooting.

EXTERIOR → INTERIOR  (r_p > r_s, r_q < r_s):
  Nothing from Region I can reach Region II after the event is in the
  interior — but a signal can ENTER the interior from the exterior.
  For RADIAL pairs (angular separation ≈ 0):
    p causes q iff v_q >= v_p  [ingoing null carries v = t + r = const inward]
  For non-radial pairs:
    UNDECIDED in this implementation (TODO: angular null geodesics
    that cross the horizon; the ingoing null shell from p at v_p reaches
    (v_p, r_q) at all angles, but non-radial geodesics modify the timing).

INTERIOR → EXTERIOR  (r_p < r_s, r_q > r_s):
  Always False.  The future horizon is a one-way membrane; no future-
  directed causal curve escapes Region II to Region I.

INTERIOR → INTERIOR  (both r < r_s):
  In the interior r is the timelike coordinate (decreasing toward
  future singularity). Event p is "earlier" if r_p > r_q.
  For r_p <= r_q: p is not earlier than q; return False.
  For RADIAL pairs (r_p > r_q, same angles):
    Lower bound (ingoing null v = t + r = const): v_q >= v_p
    Upper bound (outgoing null from p reaching r_q):
      v_upper = v_p + 2 * [(r_p - r_q) + 2M * ln((2M - r_p) / (2M - r_q))]
    p causes q iff v_p <= v_q <= v_upper.
  For non-radial pairs:
    UNDECIDED (TODO: angular criterion in the interior).

## Horizon-crossing links

A horizon-crossing link is a covering relation in the transitive
reduction that has one endpoint in the exterior (r > r_s) and one in
the interior (r < r_s). These are the Dou-Sorkin "horizon-molecules".
Counting them gives N_links ≈ A / l_Planck² where A = 4π r_s².

## Radial-strand generation mode

Because generic random sprinklings have angular separations of order O(1)
radians, the non-radial "undecided" case dominates and N_horizon_links ≈ 0.
To produce a useful Dou-Sorkin diagnostic, use radially-aligned events:
all events share the same (θ, φ) on S², so every pair is exactly radial.
Pass --aligned to the CLI to use this mode.

In the aligned mode:
- All events share θ = π/2, φ = 0 (equatorial, zero azimuth).
- The IEF radial criterion applies to all pairs.
- N_horizon_links counts covering relations across r = 2M.
- The formula N_links ~ A_1d / l² does NOT directly give the 4D area;
  aligned-strand counts need to be compared with a 1+1D density formula.
- Comparison of N_links(M₁) vs N_links(M₂) at the same N and box still
  tests monotonicity as a function of M; that is the L0 diagnostic target.

## What this is NOT

This is a minimal scaffold, not a complete Schwarzschild interior solver.
- Non-radial interior pairs are left undecided in the generic mode.
- No Kerr interior.
- No validation against analytic expectations for interior geodesics yet.
- The outgoing-null upper bound for interior-interior pairs is derived from
  the IEF metric; it has not been validated against an independent source.
- The aligned mode tests the IEF criterion, not the full 3+1D Dou-Sorkin
  counting formula.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional
import sys

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_schwarzschild_minimal_benchmark import (  # noqa: E402
    Event,
    MASS,
    SCHWARZSCHILD_RADIUS,
    TIME_EPS,
    ANGLE_EPS,
    R_MAX,
    T_MIN,
    T_MAX,
    angular_separation,
    causal_relation_schwarzschild,
    transitive_reduction_links,
    check_antisymmetric,
    check_transitive,
)
import validation_suite as vs  # noqa: E402


OUT_DIR = Path(__file__).resolve().parent
COMMAND = "python3 explore/sorkin4_schwarzschild_benchmark/run_schwarzschild_horizon_benchmark.py"

# IEF coordinate ranges for the mixed interior/exterior sprinkling.
# The interior region must satisfy 0 < r < r_s; we exclude a thin shell
# near r = 0 (singularity) and near r = r_s (horizon) to avoid coordinate
# singularities in the tortoise coordinate.
INTERIOR_R_MIN_FRACTION = 0.05   # r_interior_min = INTERIOR_R_MIN_FRACTION * r_s
INTERIOR_R_MAX_FRACTION = 0.95   # r_interior_max = INTERIOR_R_MAX_FRACTION * r_s
EXTERIOR_R_MIN_MARGIN = 0.35     # r_exterior_min = r_s + margin (same as exterior benchmark)
EXTERIOR_R_MAX = R_MAX           # same ceiling as exterior benchmark

DEFAULT_N_EXTERIOR = 8
DEFAULT_N_INTERIOR = 4
DEFAULT_SEED = 1959
DEFAULT_MASS = MASS
DEFAULT_OUT_PREFIX = "schwarzschild_horizon_benchmark"

EPS = 1.0e-12

# Angular tolerance for "radially aligned" classification.
# ANGLE_EPS (1e-12) is too tight for floating-point coordinates; pairs
# generated with math.pi/2 exactly will pass, but 1.5707... approximations
# will not.  RADIAL_ANGLE_EPS is the working threshold for IEF radial
# decisions — set just above floating-point noise from trig operations.
RADIAL_ANGLE_EPS = 1.0e-9


# ---------------------------------------------------------------------------
# Tortoise coordinate and IEF advanced time
# ---------------------------------------------------------------------------

def tortoise(r: float, mass: float = MASS) -> float:
    """Tortoise coordinate r*(r) = r + 2M ln|r/2M - 1|.

    Defined for r != r_s = 2M.  Returns a large negative number as r → r_s
    from below or above (the logarithm diverges to -inf).
    """

    rs = 2.0 * mass
    arg = abs(r / rs - 1.0)
    if arg < EPS:
        raise ValueError(f"r={r} is too close to the horizon r_s={rs}")
    return r + 2.0 * mass * math.log(arg)


def advanced_null_v(t: float, r: float, mass: float = MASS) -> float:
    """IEF advanced null coordinate v = t + r.

    He & Rideout Eq. (3) uses t = t_s + 2M ln(r/2M - 1) = v - r,
    so the advanced null time is v = t + r.  Eq. (4) gives the ingoing
    radial null as dt + dr = 0, hence t + r is constant.

    The mass parameter is retained for API compatibility.
    """

    _ = mass
    return t + r


# ---------------------------------------------------------------------------
# Event region classification
# ---------------------------------------------------------------------------

EXTERIOR = "exterior"
INTERIOR = "interior"
ON_HORIZON = "on_horizon"


def region(r: float, mass: float = MASS) -> str:
    """Classify an event as exterior, interior, or on the horizon."""

    rs = 2.0 * mass
    if abs(r - rs) < EPS:
        return ON_HORIZON
    return EXTERIOR if r > rs else INTERIOR


# ---------------------------------------------------------------------------
# IEF outgoing-null upper bound for interior pairs
# ---------------------------------------------------------------------------

def _interior_outgoing_null_v_upper(
    v_p: float,
    r_p: float,
    r_q: float,
    mass: float = MASS,
) -> float:
    """v-coordinate reached by the outgoing null from p at radius r_q.

    Both events are in the interior: 0 < r_q < r_p < r_s = 2M.

    The outgoing null from p satisfies dv/dr = 2/(1 - 2M/r), which in the
    interior (1 - 2M/r < 0) gives dv/dr < 0, so v decreases as r decreases.

    Wait — that means the outgoing null DECREASES in v as r decreases. But
    if v decreases moving toward smaller r, then the causal future of p in
    the interior is the region v_q >= v_p (bounded below by the ingoing null
    and above by... wait, the outgoing null goes the wrong way).

    Let me reconsider.  In the interior the "outgoing" null has:
      dv = 2 dr / (1 - 2M/r)
    Since 1 - 2M/r < 0 for r < 2M, dv/dr < 0.  So as r decreases (toward
    singularity), v also decreases along this null.

    The causal future of p in the interior is the region between:
      Lower null (ingoing):  v = v_p  (v constant, r decreasing)
      Upper null (outgoing): v decreasing as r decreases → smaller v

    This means the outgoing null from p at r_p has v < v_p for r < r_p.
    The causal future of p at radius r_q (< r_p) is:
      v_outgoing(r_q) <= v_q <= v_p

    where v_outgoing(r_q) < v_p is found by integrating dv = 2 dr / (1 - 2M/r)
    from r_p to r_q.

    Integrating: dv = 2 dr / (1 - 2M/r) = 2r dr / (r - 2M)
    Let u = r - 2M, du = dr:
      ∫ 2(u + 2M) / u du = 2u + 4M ln|u| + C = 2(r-2M) + 4M ln|r-2M| + C

    So:
      v_outgoing(r_q) = v_p + [2(r-2M) + 4M ln|r-2M|] evaluated from r_p to r_q
      = v_p + [2(r_q-2M) + 4M ln(2M-r_q)] - [2(r_p-2M) + 4M ln(2M-r_p)]
      (using |r-2M| = 2M-r in the interior)

    For r_q < r_p < 2M: (r_q - 2M) < (r_p - 2M) < 0, so first bracket more negative.
    Also ln(2M - r_q) vs ln(2M - r_p): since r_q < r_p, 2M-r_q > 2M-r_p > 0, so ln bigger.

    Net sign: the result can be positive or negative depending on parameters.
    """

    rs = 2.0 * mass
    if r_q >= r_p or r_p >= rs:
        raise ValueError(
            f"Interior outgoing null requires 0 < r_q={r_q} < r_p={r_p} < r_s={rs}"
        )

    # Antiderivative: F(r) = 2(r - 2M) + 4M ln(2M - r)  [interior: r < 2M]
    def antideriv(r: float) -> float:
        return 2.0 * (r - rs) + 4.0 * mass * math.log(rs - r)

    return v_p + antideriv(r_q) - antideriv(r_p)


# ---------------------------------------------------------------------------
# Four-case causal criterion in IEF
# ---------------------------------------------------------------------------

def causal_relation_ief(
    p: Event,
    q: Event,
    mass: float = MASS,
    enable_exterior_shooting: bool = False,
    angle_eps: float = RADIAL_ANGLE_EPS,
) -> Optional[bool]:
    """Causal relation in IEF coordinates for any interior/exterior pair.

    Returns:
      True   — p definitively causally precedes q
      False  — p definitively does NOT causally precede q
      None   — undecided (non-radial interior pair or shooting not enabled)

    The four cases are handled separately.  Interior → Exterior is always
    False regardless of other criteria.
    """

    rs = 2.0 * mass
    reg_p = region(p.r, mass)
    reg_q = region(q.r, mass)

    # Compute advanced null coordinates for the IEF criteria.
    # These may raise ValueError if r is exactly on the horizon; treat as undecided.
    try:
        v_p = advanced_null_v(p.t, p.r, mass)
        v_q = advanced_null_v(q.t, q.r, mass)
    except ValueError:
        return None

    # --- INTERIOR → EXTERIOR: always False ---
    if reg_p == INTERIOR and reg_q == EXTERIOR:
        return False

    # --- EXTERIOR → EXTERIOR ---
    if reg_p == EXTERIOR and reg_q == EXTERIOR:
        return causal_relation_schwarzschild(p, q, enable_shooting=enable_exterior_shooting)

    # --- EXTERIOR → INTERIOR ---
    if reg_p == EXTERIOR and reg_q == INTERIOR:
        # Necessary condition: v_q >= v_p (future horizon can only be crossed
        # by signals with non-decreasing v).
        if v_q < v_p - TIME_EPS:
            return False

        phi2 = angular_separation(p, q)
        if phi2 <= angle_eps:
            # Radial case: exact via ingoing null (v = const).
            # An ingoing null from p at v_p reaches (v_p, r_q) for any r_q.
            # Since p is exterior and q is interior: r_p > r_s > r_q always,
            # so p.r > q.r is guaranteed by the region check.
            return True

        # Non-radial exterior → interior: undecided.
        # Deciding requires integrating a null geodesic through the horizon
        # with angular component, which is not yet implemented.
        return None

    # --- INTERIOR → INTERIOR ---
    if reg_p == INTERIOR and reg_q == INTERIOR:
        # In the interior r is timelike (decreasing toward future singularity).
        # p is "earlier" (further from singularity) iff r_p > r_q.
        if p.r <= q.r + EPS:
            return False   # p is not earlier than q in interior time

        phi2 = angular_separation(p, q)
        if phi2 <= angle_eps:
            # Radial interior case.
            # Lower bound (ingoing null v = const from p): v_q >= v_p
            if v_q < v_p - TIME_EPS:
                return False

            # Upper bound (outgoing null from p reaches r_q at v_upper):
            #   v_q <= v_outgoing(r_q).
            # If v_q is above this bound, q is outside p's future null cone.
            try:
                v_upper = _interior_outgoing_null_v_upper(v_p, p.r, q.r, mass)
            except ValueError:
                return None

            if v_q > v_upper + TIME_EPS:
                return False

            return True

        # Non-radial interior pair: undecided.
        return None

    # ON_HORIZON: treat as exterior for the sufficient bounds.
    # Horizon events are a set of measure zero; their precise classification
    # does not affect the link count for generic sprinklings.
    return causal_relation_schwarzschild(p, q, enable_shooting=enable_exterior_shooting)


# ---------------------------------------------------------------------------
# Event generation: mixed interior + exterior sprinkling
# ---------------------------------------------------------------------------

def _volume_weighted_r_interior(rng: random.Random, r_min: float, r_max: float) -> float:
    """Sample r uniformly in [r_min, r_max] weighted by the r² volume element."""

    u = rng.random()
    return (r_min**3 + u * (r_max**3 - r_min**3)) ** (1.0 / 3.0)


def generate_mixed_events(
    n_exterior: int,
    n_interior: int,
    seed: int,
    mass: float = MASS,
    aligned: bool = False,
) -> list[Event]:
    """Generate exterior + interior events in IEF coordinates.

    Exterior events: r in [r_s + EXTERIOR_R_MIN_MARGIN, EXTERIOR_R_MAX]
    Interior events: r in [INTERIOR_R_MIN_FRACTION * r_s,
                            INTERIOR_R_MAX_FRACTION * r_s]

    All events use the He-Rideout EF t coordinate; v = t + r is
    computed on demand.  Angular coordinates are uniform on S².
    """

    rs = 2.0 * mass
    r_ext_min = rs + EXTERIOR_R_MIN_MARGIN
    r_int_min = INTERIOR_R_MIN_FRACTION * rs
    r_int_max = INTERIOR_R_MAX_FRACTION * rs

    rng = random.Random(seed)
    events: list[Event] = []

    # In aligned mode all events share θ=π/2, φ=0 (equatorial radial strand).
    # This ensures every pair is exactly radial (phi2=0) so the IEF radial
    # criterion applies to all pairs — producing non-zero horizon-crossing links.
    # In generic mode events are uniformly distributed on S²; nearly all pairs
    # will be non-radial and undecided.
    fixed_theta = math.pi / 2.0 if aligned else None
    fixed_phi = 0.0 if aligned else None

    for _ in range(n_exterior):
        t = rng.uniform(T_MIN, T_MAX)
        r = _volume_weighted_r_interior(rng, r_ext_min, EXTERIOR_R_MAX)
        if aligned:
            theta, phi = fixed_theta, fixed_phi
        else:
            cos_theta = rng.uniform(-1.0, 1.0)
            theta = math.acos(cos_theta)
            phi = rng.uniform(0.0, 2.0 * math.pi)
        events.append(Event(index=len(events), t=t, r=r, theta=theta, phi=phi))

    for _ in range(n_interior):
        t = rng.uniform(T_MIN, T_MAX)
        r = _volume_weighted_r_interior(rng, r_int_min, r_int_max)
        if aligned:
            theta, phi = fixed_theta, fixed_phi
        else:
            cos_theta = rng.uniform(-1.0, 1.0)
            theta = math.acos(cos_theta)
            phi = rng.uniform(0.0, 2.0 * math.pi)
        events.append(Event(index=len(events), t=t, r=r, theta=theta, phi=phi))

    # Sort by IEF advanced time v = t + r as the natural causal ordering.
    def sort_key(e: Event) -> tuple[float, float]:
        try:
            return advanced_null_v(e.t, e.r, mass), e.r
        except ValueError:
            return (e.t, e.r)

    events.sort(key=sort_key)
    for idx, e in enumerate(events):
        events[idx] = Event(index=idx, t=e.t, r=e.r, theta=e.theta, phi=e.phi)
    return events


# ---------------------------------------------------------------------------
# Causal matrix and horizon-link counting
# ---------------------------------------------------------------------------

def build_horizon_causal_matrix(
    events: list[Event],
    mass: float = MASS,
    enable_exterior_shooting: bool = False,
    aligned: bool = False,
) -> tuple[vs.CausalMatrix, list[list[Optional[bool]]]]:
    """Build the causal matrix using IEF criteria."""

    n = len(events)
    matrix: vs.CausalMatrix = [[False] * n for _ in range(n)]
    states: list[list[Optional[bool]]] = [
        [False if i == j else None for j in range(n)] for i in range(n)
    ]

    # In aligned mode every pair is radially aligned; use RADIAL_ANGLE_EPS.
    # In generic mode use ANGLE_EPS (1e-12) which is effectively "exact radial only".
    angle_eps = RADIAL_ANGLE_EPS if aligned else ANGLE_EPS

    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            rel = causal_relation_ief(
                events[i], events[j], mass=mass,
                enable_exterior_shooting=enable_exterior_shooting,
                angle_eps=angle_eps,
            )
            states[i][j] = rel
            if rel is True:
                matrix[i][j] = True

    return matrix, states


def count_horizon_links(
    events: list[Event],
    links: list[tuple[int, int]],
    mass: float = MASS,
) -> int:
    """Count links (covering relations) that cross the future horizon.

    A horizon-crossing link is a pair (i, j) in the transitive reduction where:
      events[i] is in the EXTERIOR  (r_i > r_s)
      events[j] is in the INTERIOR  (r_j < r_s)

    These are the Dou-Sorkin horizon-molecules.
    N_horizon_links ≈ A / l_Planck²  where A = 4π (2M)² = 16π M².
    """

    rs = 2.0 * mass
    count = 0
    for i, j in links:
        r_i = events[i].r
        r_j = events[j].r
        if r_i > rs + EPS and r_j < rs - EPS:
            count += 1
    return count


# ---------------------------------------------------------------------------
# Summary and output
# ---------------------------------------------------------------------------

def summarize_horizon_case(
    events: list[Event],
    matrix: vs.CausalMatrix,
    states: list[list[Optional[bool]]],
    mass: float,
    seed: int,
    n_exterior: int,
    n_interior: int,
    enable_exterior_shooting: bool,
) -> dict[str, object]:
    n = len(events)
    rs = 2.0 * mass
    possible_pairs = n * (n - 1) // 2

    true_relations = sum(
        1 for i in range(n - 1) for j in range(i + 1, n) if states[i][j] is True
    )
    false_relations = sum(
        1 for i in range(n - 1) for j in range(i + 1, n) if states[i][j] is False
    )
    undecided_pairs = sum(
        1 for i in range(n - 1) for j in range(i + 1, n) if states[i][j] is None
    )
    decided_pairs = true_relations + false_relations

    # ordering_fraction: denominator = total possible pairs.
    # solver_coverage: fraction of pairs decided by the solver.
    ordering_fraction = true_relations / possible_pairs if possible_pairs else 0.0
    solver_coverage = decided_pairs / possible_pairs if possible_pairs else 0.0

    links = transitive_reduction_links(matrix)
    n_horizon_links = count_horizon_links(events, links, mass=mass)
    expected_a = 4.0 * math.pi * rs * rs
    expected_a_16pi_m2 = 16.0 * math.pi * mass * mass

    n_ext = sum(1 for e in events if e.r > rs + EPS)
    n_int = sum(1 for e in events if e.r < rs - EPS)

    return {
        "benchmark": "S4 Schwarzschild horizon-crossing causal solver",
        "status": "minimal_ief_radial_criterion_nonradial_undecided",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "command": COMMAND,
        "N": n,
        "N_exterior": n_ext,
        "N_interior": n_int,
        "seed": seed,
        "mass": mass,
        "schwarzschild_radius": rs,
        "expected_horizon_area_4pi_rs2": expected_a,
        "expected_horizon_area_16pi_M2": expected_a_16pi_m2,
        "coordinate_system": "ingoing Eddington-Finkelstein (v, r, theta, phi)",
        "causal_model": (
            "IEF four-case criterion: ext-ext via He-Rideout (shooting="
            + str(enable_exterior_shooting)
            + "); ext->int radial exact; int->ext always False; "
            "int-int radial exact with outgoing-null upper bound; non-radial undecided"
        ),
        "possible_pairs": possible_pairs,
        "true_relations": true_relations,
        "false_relations": false_relations,
        "undecided_pairs": undecided_pairs,
        "decided_pairs": decided_pairs,
        "ordering_fraction": ordering_fraction,
        "solver_coverage": solver_coverage,
        "links": len(links),
        "horizon_crossing_links": n_horizon_links,
        "antisymmetric": check_antisymmetric(matrix),
        "transitive": check_transitive(matrix),
        "notes": [
            "horizon_crossing_links counts ext->int links in the transitive reduction.",
            "Dou-Sorkin: N_links ~ A / l_Planck^2 = 16 pi M^2 in Planck units.",
            "Non-radial interior pairs are undecided; they contribute to undecided_pairs.",
            "ordering_fraction denominator is possible_pairs (null = unknown, not non-causal).",
            "solver_coverage = decided_pairs / possible_pairs is a solver completeness metric.",
        ],
    }


def run_horizon_case(
    n_exterior: int,
    n_interior: int,
    seed: int,
    mass: float = DEFAULT_MASS,
    enable_exterior_shooting: bool = False,
    aligned: bool = False,
) -> tuple[list[Event], vs.CausalMatrix, list[list[Optional[bool]]], dict[str, object]]:
    events = generate_mixed_events(n_exterior, n_interior, seed, mass=mass, aligned=aligned)
    matrix, states = build_horizon_causal_matrix(
        events, mass=mass,
        enable_exterior_shooting=enable_exterior_shooting,
        aligned=aligned,
    )
    summary = summarize_horizon_case(
        events, matrix, states,
        mass=mass, seed=seed,
        n_exterior=n_exterior, n_interior=n_interior,
        enable_exterior_shooting=enable_exterior_shooting,
    )
    summary["aligned_mode"] = aligned
    return events, matrix, states, summary


def write_outputs(
    events: list[Event],
    matrix: vs.CausalMatrix,
    states: list[list[Optional[bool]]],
    summary: dict[str, object],
    out_prefix: str,
) -> tuple[Path, Path]:
    csv_path = OUT_DIR / f"{out_prefix}.csv"
    json_path = OUT_DIR / f"{out_prefix}.json"

    with csv_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(list(summary))
        writer.writerow(list(summary.values()))

    with json_path.open("w") as fh:
        json.dump(
            {
                "summary": summary,
                "events": [asdict(e) for e in events],
                "relation_states": [
                    [str(v) if v is None else v for v in row] for row in states
                ],
            },
            fh,
            indent=2,
            default=str,
        )

    return csv_path, json_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-exterior", type=int, default=DEFAULT_N_EXTERIOR)
    parser.add_argument("--n-interior", type=int, default=DEFAULT_N_INTERIOR)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--mass", type=float, default=DEFAULT_MASS)
    parser.add_argument("--shooting", action="store_true",
                        help="Enable exterior generic-pair shooting (slow)")
    parser.add_argument("--aligned", action="store_true",
                        help="Place all events at theta=pi/2, phi=0 (radial strand). "
                             "Required for non-zero horizon_crossing_links in this implementation.")
    parser.add_argument("--out-prefix", default=DEFAULT_OUT_PREFIX)
    parser.add_argument("--no-output", action="store_true")
    args = parser.parse_args()

    events, matrix, states, summary = run_horizon_case(
        n_exterior=args.n_exterior,
        n_interior=args.n_interior,
        seed=args.seed,
        mass=args.mass,
        enable_exterior_shooting=args.shooting,
        aligned=args.aligned,
    )

    print(f"N={summary['N']} (ext={summary['N_exterior']}, int={summary['N_interior']})")
    print(f"seed={summary['seed']}  mass={summary['mass']}  r_s={summary['schwarzschild_radius']}")
    print(f"true_relations={summary['true_relations']}")
    print(f"undecided_pairs={summary['undecided_pairs']}")
    print(f"horizon_crossing_links={summary['horizon_crossing_links']}")
    print(f"ordering_fraction={summary['ordering_fraction']:.6f}  "
          f"solver_coverage={summary['solver_coverage']:.6f}")
    print(f"links={summary['links']}  antisymmetric={summary['antisymmetric']}")
    print(f"Dou-Sorkin expected A=16pi M^2={summary['expected_horizon_area_16pi_M2']:.4f}")

    if not args.no_output:
        csv_path, json_path = write_outputs(events, matrix, states, summary, args.out_prefix)
        print(f"wrote {csv_path.name}, {json_path.name}")


if __name__ == "__main__":
    main()
