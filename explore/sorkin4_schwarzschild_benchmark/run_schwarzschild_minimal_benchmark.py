#!/usr/bin/env python3
"""Minimal Schwarzschild exterior benchmark for SORKIN-4.

This is deliberately not a Kerr benchmark and deliberately not a fake
Schwarzschild causal-relation solver.  It implements only the He & Rideout
radial exact tests and sufficient exterior bounds that are short and directly
traceable to their paper.  Generic pairs remain undecided until the numerical
null-geodesic shooting procedure is implemented.
"""

from __future__ import annotations

import csv
import json
import math
import random
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import validation_suite as vs  # noqa: E402


OUT_DIR = Path(__file__).resolve().parent
CSV_PATH = OUT_DIR / "schwarzschild_minimal_benchmark.csv"
JSON_PATH = OUT_DIR / "schwarzschild_minimal_benchmark.json"
COMMAND = "python3 explore/sorkin4_schwarzschild_benchmark/run_schwarzschild_minimal_benchmark.py"

N = 12
SEED = 1959
SCHWARZSCHILD_RADIUS = 2.0
MASS = SCHWARZSCHILD_RADIUS / 2.0
EXTERIOR_MARGIN = 0.35
R_MIN = SCHWARZSCHILD_RADIUS + EXTERIOR_MARGIN
R_MAX = 6.0
T_MIN = 0.0
T_MAX = 4.0
ANGLE_EPS = 1.0e-12
TIME_EPS = 1.0e-12


@dataclass(frozen=True)
class Event:
    """Known-coordinate event in exterior Eddington-Finkelstein form."""

    index: int
    t: float
    r: float
    theta: float
    phi: float


def _sample_radius_volume_weighted(rng: random.Random) -> float:
    """Sample r with density proportional to r^2 on [R_MIN, R_MAX]."""

    u = rng.random()
    return (R_MIN**3 + u * (R_MAX**3 - R_MIN**3)) ** (1.0 / 3.0)


def generate_exterior_events(n: int, seed: int) -> list[Event]:
    """Generate a reproducible bounded exterior Schwarzschild point set.

    The coordinate patch is outside the horizon only: r > r_s + margin.
    The radial/angular sampling follows the Schwarzschild volume element
    factor r^2 sin(theta) in a bounded coordinate box.  This is only the
    event-generation scaffold; pairwise causal relations are not inferred
    until the He & Rideout null-geodesic decision procedure is implemented.
    """

    rng = random.Random(seed)
    events: list[Event] = []
    for index in range(n):
        t = rng.uniform(T_MIN, T_MAX)
        r = _sample_radius_volume_weighted(rng)
        cos_theta = rng.uniform(-1.0, 1.0)
        theta = math.acos(cos_theta)
        phi = rng.uniform(0.0, 2.0 * math.pi)
        events.append(Event(index=index, t=t, r=r, theta=theta, phi=phi))
    return sorted(events, key=lambda event: (event.t, event.r, event.theta, event.phi))


def angular_separation(p: Event, q: Event) -> float:
    """Great-circle angular separation after the He-Rideout rotation step."""

    cos_angle = (
        math.cos(p.theta) * math.cos(q.theta)
        + math.sin(p.theta) * math.sin(q.theta) * math.cos(p.phi - q.phi)
    )
    return math.acos(max(-1.0, min(1.0, cos_angle)))


def _outgoing_radial_trip(r1: float, r2: float) -> float:
    """EF-time for an outgoing radial null ray outside the horizon."""

    return r2 - r1 + 4.0 * MASS * math.log((r2 - SCHWARZSCHILD_RADIUS) / (r1 - SCHWARZSCHILD_RADIUS))


def _angular_f(r: float) -> float:
    """He-Rideout f(r)=r/sqrt(1-2M/r), valid only outside the horizon."""

    return r / math.sqrt(1.0 - SCHWARZSCHILD_RADIUS / r)


def _angular_spacelike_r0(r1: float, r2: float) -> float:
    """Minimizer of f(r) on the outgoing exterior interval [r1, r2]."""

    photon_radius = 3.0 * MASS
    if photon_radius <= r1 <= r2:
        return r1
    if r1 < photon_radius < r2:
        return photon_radius
    return r2


def _timelike_bound_r0(r1: float, r2: float) -> float:
    """He-Rideout r0 for the composed null-curve sufficient timelike bound."""

    photon_radius = 3.0 * MASS
    if r1 >= photon_radius and r2 >= photon_radius:
        return min(r1, r2)
    if (r1 > photon_radius > r2) or (r2 > photon_radius > r1):
        return photon_radius
    return max(r1, r2)


def causal_relation_schwarzschild(p: Event, q: Event) -> Optional[bool]:
    """Return whether p causally precedes q in the implemented S4-2 subset.

    Implemented from He & Rideout (2009):
    - exact radial null-geodesic criteria for angular separation zero;
    - radial and angular sufficient spacelike bounds;
    - composed radial/angular null-curve sufficient timelike bounds.

    TODO(S4-2): implement the generic-pair numerical shooting procedure:
    solve for c^2=(E/L)^2 so the integral of dphi/du reaches the rotated
    angular separation, then integrate dt/du and compare arrival time with
    q.t.  Until then, None means "generic pair not decided"; callers must not
    treat undecided pairs as physical non-relations.
    """

    dt = q.t - p.t
    if dt < -TIME_EPS:
        return False

    r1 = p.r
    r2 = q.r
    phi2 = angular_separation(p, q)
    horizon = SCHWARZSCHILD_RADIUS

    if r2 > r1 and r1 < horizon:
        return False

    if phi2 <= ANGLE_EPS:
        if r1 >= r2:
            return dt + TIME_EPS >= r1 - r2
        if r2 >= r1 > horizon:
            return dt + TIME_EPS >= _outgoing_radial_trip(r1, r2)
        return False

    if r1 >= r2:
        radial_lower_bound = r1 - r2
        if dt < radial_lower_bound - TIME_EPS:
            return False
        angular_lower_bound = r2 * phi2
        if dt < angular_lower_bound - TIME_EPS:
            return False
    elif r2 >= r1 > horizon:
        radial_lower_bound = _outgoing_radial_trip(r1, r2)
        if dt < radial_lower_bound - TIME_EPS:
            return False
        r0_spacelike = _angular_spacelike_r0(r1, r2)
        angular_lower_bound = _angular_f(r0_spacelike) * phi2
        if dt < angular_lower_bound - TIME_EPS:
            return False

    r0_timelike = _timelike_bound_r0(r1, r2)
    angular_trip = _angular_f(r0_timelike) * phi2
    if r1 >= r2 and r1 > horizon:
        composed_trip = r1 - r2 + angular_trip
        if dt + TIME_EPS >= composed_trip:
            return True
    elif r2 >= r1 > horizon:
        composed_trip = _outgoing_radial_trip(r1, r2) + angular_trip
        if dt + TIME_EPS >= composed_trip:
            return True

    return None


def build_causal_matrix(events: list[Event]) -> tuple[vs.CausalMatrix, list[list[Optional[bool]]]]:
    """Build the asserted causal matrix from implemented Schwarzschild tests."""

    n = len(events)
    matrix: vs.CausalMatrix = [[False] * n for _ in range(n)]
    states: list[list[Optional[bool]]] = [[False if i == j else None for j in range(n)] for i in range(n)]
    for i in range(n - 1):
        for j in range(i + 1, n):
            relation = causal_relation_schwarzschild(events[i], events[j])
            states[i][j] = relation
            if relation:
                matrix[i][j] = True
    return matrix, states


def count_relations(matrix: vs.CausalMatrix) -> int:
    return sum(1 for i, row in enumerate(matrix) for j, value in enumerate(row) if i < j and value)


def check_antisymmetric(matrix: vs.CausalMatrix) -> bool:
    n = len(matrix)
    for i in range(n):
        if matrix[i][i]:
            return False
        for j in range(i + 1, n):
            if matrix[i][j] and matrix[j][i]:
                return False
    return True


def check_transitive(matrix: vs.CausalMatrix) -> bool:
    n = len(matrix)
    for i in range(n):
        for j in range(n):
            if not matrix[i][j]:
                continue
            for k in range(n):
                if matrix[j][k] and not matrix[i][k]:
                    return False
    return True


def check_decided_transitivity(states: list[list[Optional[bool]]]) -> bool:
    """Return False only when decided true chains contradict a decided false pair."""

    n = len(states)
    for i in range(n):
        for j in range(n):
            if states[i][j] is not True:
                continue
            for k in range(n):
                if states[j][k] is True and states[i][k] is False:
                    return False
    return True


def transitive_reduction_links(matrix: vs.CausalMatrix) -> list[tuple[int, int]]:
    """Return links for an already-transitive finite order matrix."""

    n = len(matrix)
    links: list[tuple[int, int]] = []
    for i in range(n):
        for j in range(n):
            if not matrix[i][j]:
                continue
            has_intermediate = any(matrix[i][k] and matrix[k][j] for k in range(n))
            if not has_intermediate:
                links.append((i, j))
    return links


def _fmt(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return f"{value:.12g}"
    return str(value)


def write_outputs(
    events: list[Event],
    matrix: vs.CausalMatrix,
    states: list[list[Optional[bool]]],
) -> dict[str, object]:
    n = len(events)
    possible_pairs = n * (n - 1) // 2
    true_relations = count_relations(matrix)
    false_relations = sum(1 for i in range(n - 1) for j in range(i + 1, n) if states[i][j] is False)
    undecided_pairs = sum(1 for i in range(n - 1) for j in range(i + 1, n) if states[i][j] is None)
    decided_pairs = true_relations + false_relations
    ordering_fraction = true_relations / decided_pairs if decided_pairs else 0.0
    links = transitive_reduction_links(matrix)
    self_comparison = vs.compare_causal_orders(matrix, matrix)

    summary: dict[str, object] = {
        "benchmark": "Schwarzschild benchmark before Kerr ordinal diagnostics.",
        "status": "partial_he_rideout_exterior_bounds_generic_pairs_undecided",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "command": COMMAND,
        "N": n,
        "seed": SEED,
        "schwarzschild_radius": SCHWARZSCHILD_RADIUS,
        "mass": MASS,
        "exterior_margin": EXTERIOR_MARGIN,
        "r_min": R_MIN,
        "r_max": R_MAX,
        "t_min": T_MIN,
        "t_max": T_MAX,
        "coordinate_patch": "exterior Eddington-Finkelstein coordinates",
        "causal_relation_model": "He & Rideout radial exact tests plus sufficient bounds; generic geodesic shooting TODO",
        "true_relations": true_relations,
        "false_relations": false_relations,
        "decided_pairs": decided_pairs,
        "relations": true_relations,
        "ordering_fraction": ordering_fraction,
        "ordering_fraction_denominator": "decided_pairs",
        "links": len(links),
        "transitive_reduction_implemented": True,
        "antisymmetric": check_antisymmetric(matrix),
        "transitive_true_matrix": check_transitive(matrix),
        "decided_transitivity_no_false_contradictions": check_decided_transitivity(states),
        "possible_pairs": possible_pairs,
        "undecided_pairs": undecided_pairs,
        "warning": "undecided generic pairs remain" if undecided_pairs else "",
        "self_compare_exact_match": self_comparison.exact_match,
    }

    headers = list(summary)
    with CSV_PATH.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(headers)
        writer.writerow([_fmt(summary[header]) for header in headers])

    payload = {
        "summary": summary,
        "events": [asdict(event) for event in events],
        "causal_matrix": matrix,
        "relation_states": states,
        "links": links,
        "notes": [
            "S4-2 asserts only relations decided by He & Rideout radial exact tests or sufficient bounds.",
            "False relation_states entries are decided non-relations; null entries are generic undecided pairs.",
            "False causal_matrix entries include both decided false and undecided pairs; inspect relation_states for the distinction.",
            "The next required step is He & Rideout generic null-geodesic shooting.",
        ],
    }
    JSON_PATH.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return summary


def main() -> None:
    events = generate_exterior_events(N, SEED)
    matrix, states = build_causal_matrix(events)
    summary = write_outputs(events, matrix, states)

    print("S4-2 Schwarzschild minimal benchmark")
    print("Regime: exterior only, r > r_s + margin")
    print(f"N={summary['N']} seed={summary['seed']}")
    print(f"M={summary['mass']} r_s={summary['schwarzschild_radius']} margin={summary['exterior_margin']}")
    print(
        f"true_relations={summary['true_relations']} "
        f"false_relations={summary['false_relations']} "
        f"undecided_pairs={summary['undecided_pairs']}"
    )
    print(f"ordering_fraction_decided={summary['ordering_fraction']:.6g}")
    print(f"links={summary['links']} undecided_pairs={summary['undecided_pairs']}")
    print(
        f"antisymmetric={summary['antisymmetric']} "
        f"transitive_true_matrix={summary['transitive_true_matrix']} "
        "decided_transitivity_no_false_contradictions="
        f"{summary['decided_transitivity_no_false_contradictions']}"
    )
    if summary["warning"]:
        print(f"warning={summary['warning']}")
    print(f"causal_relation_model={summary['causal_relation_model']}")
    print(f"wrote {CSV_PATH}")
    print(f"wrote {JSON_PATH}")


if __name__ == "__main__":
    main()
