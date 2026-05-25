#!/usr/bin/env python3
"""Minimal Schwarzschild exterior benchmark scaffold for SORKIN-4.

This is deliberately not a Kerr benchmark and deliberately not a fake
Schwarzschild causal-relation solver.  It creates a small reproducible
known-coordinate exterior sprinkling and writes the same output structure
that a later He & Rideout-style causal-relation implementation should fill.
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
EXTERIOR_MARGIN = 0.35
R_MIN = SCHWARZSCHILD_RADIUS + EXTERIOR_MARGIN
R_MAX = 6.0
V_MIN = 0.0
V_MAX = 4.0


@dataclass(frozen=True)
class Event:
    """Known-coordinate event in exterior ingoing Eddington-Finkelstein form."""

    index: int
    v: float
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
        v = rng.uniform(V_MIN, V_MAX)
        r = _sample_radius_volume_weighted(rng)
        cos_theta = rng.uniform(-1.0, 1.0)
        theta = math.acos(cos_theta)
        phi = rng.uniform(0.0, 2.0 * math.pi)
        events.append(Event(index=index, v=v, r=r, theta=theta, phi=phi))
    return sorted(events, key=lambda event: (event.v, event.r, event.theta, event.phi))


def causal_relation_schwarzschild(p: Event, q: Event) -> Optional[bool]:
    """Return whether p causally precedes q in Schwarzschild, once implemented.

    TODO(S4-1): implement the He & Rideout (2009) Schwarzschild pairwise
    causal-relation algorithm: cheap sufficient spacelike/timelike tests,
    then numerical null-geodesic arrival checks for generic exterior pairs.
    Until then, returning None means "not decided"; callers must not treat
    undecided pairs as physical non-relations.
    """

    _ = (p, q)
    return None


def build_causal_matrix(events: list[Event]) -> tuple[vs.CausalMatrix, int]:
    """Build the asserted causal matrix from implemented Schwarzschild tests."""

    n = len(events)
    matrix: vs.CausalMatrix = [[False] * n for _ in range(n)]
    undecided_pairs = 0
    for i in range(n - 1):
        for j in range(i + 1, n):
            relation = causal_relation_schwarzschild(events[i], events[j])
            if relation is None:
                undecided_pairs += 1
            elif relation:
                matrix[i][j] = True
    return matrix, undecided_pairs


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


def write_outputs(events: list[Event], matrix: vs.CausalMatrix, undecided_pairs: int) -> dict[str, object]:
    n = len(events)
    possible_pairs = n * (n - 1) // 2
    relation_count = count_relations(matrix)
    ordering_fraction = relation_count / possible_pairs if possible_pairs else 0.0
    links = transitive_reduction_links(matrix)
    self_comparison = vs.compare_causal_orders(matrix, matrix)

    summary: dict[str, object] = {
        "benchmark": "Schwarzschild benchmark before Kerr ordinal diagnostics.",
        "status": "scaffold_no_physical_causal_relations_asserted",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "command": COMMAND,
        "N": n,
        "seed": SEED,
        "schwarzschild_radius": SCHWARZSCHILD_RADIUS,
        "exterior_margin": EXTERIOR_MARGIN,
        "r_min": R_MIN,
        "r_max": R_MAX,
        "v_min": V_MIN,
        "v_max": V_MAX,
        "coordinate_patch": "exterior ingoing Eddington-Finkelstein-like coordinates",
        "causal_relation_model": "TODO He & Rideout 2009 null-geodesic decision procedure",
        "relations": relation_count,
        "ordering_fraction": ordering_fraction,
        "links": len(links),
        "transitive_reduction_implemented": True,
        "antisymmetric": check_antisymmetric(matrix),
        "transitive": check_transitive(matrix),
        "possible_pairs": possible_pairs,
        "undecided_pairs": undecided_pairs,
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
        "links": links,
        "notes": [
            "No pairwise Schwarzschild causal relation is asserted in S4-1.",
            "False matrix entries mean unasserted in this scaffold, not physically spacelike.",
            "The next required step is the He & Rideout Schwarzschild relation algorithm.",
        ],
    }
    JSON_PATH.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return summary


def main() -> None:
    events = generate_exterior_events(N, SEED)
    matrix, undecided_pairs = build_causal_matrix(events)
    summary = write_outputs(events, matrix, undecided_pairs)

    print("S4-1 Schwarzschild minimal benchmark")
    print("Regime: exterior only, r > r_s + margin")
    print(f"N={summary['N']} seed={summary['seed']}")
    print(f"r_s={summary['schwarzschild_radius']} margin={summary['exterior_margin']}")
    print(f"relations={summary['relations']} ordering_fraction={summary['ordering_fraction']:.6g}")
    print(f"links={summary['links']} undecided_pairs={summary['undecided_pairs']}")
    print(f"antisymmetric={summary['antisymmetric']} transitive={summary['transitive']}")
    print(f"causal_relation_model={summary['causal_relation_model']}")
    print(f"wrote {CSV_PATH}")
    print(f"wrote {JSON_PATH}")


if __name__ == "__main__":
    main()
