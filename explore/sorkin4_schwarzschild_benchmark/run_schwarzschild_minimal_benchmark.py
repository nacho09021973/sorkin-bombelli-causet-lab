#!/usr/bin/env python3
"""Minimal Schwarzschild exterior benchmark for SORKIN-4.

This is deliberately not a Kerr benchmark and deliberately not a fake
Schwarzschild causal-relation solver.  It implements only the He & Rideout
radial exact tests and sufficient exterior bounds that are short and directly
traceable to their paper.  Generic pairs remain undecided until the numerical
null-geodesic shooting procedure is implemented.
"""

from __future__ import annotations

import argparse
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
COMMAND = "python3 explore/sorkin4_schwarzschild_benchmark/run_schwarzschild_minimal_benchmark.py"
DEFAULT_OUT_PREFIX = "schwarzschild_minimal_benchmark"

DEFAULT_N = 12
DEFAULT_SEED = 1959
SCHWARZSCHILD_RADIUS = 2.0
MASS = SCHWARZSCHILD_RADIUS / 2.0
EXTERIOR_MARGIN = 0.35
R_MIN = SCHWARZSCHILD_RADIUS + EXTERIOR_MARGIN
R_MAX = 6.0
T_MIN = 0.0
T_MAX = 4.0
ANGLE_EPS = 1.0e-12
TIME_EPS = 1.0e-12
CUBIC_ROOT_EPS = 1.0e-12
SHOOTING_PHI_TOL = 1.0e-9
SHOOTING_MAX_ITER = 80
SIMPSON_INTERVALS = 512
ROW_HEADERS = (
    "N",
    "seed",
    "true_relations",
    "false_relations",
    "undecided_pairs",
    "decided_pairs",
    "ordering_fraction_decided",
    "links",
    "antisymmetric",
    "transitive_true_matrix",
    "decided_transitivity_no_false_contradictions",
)
DEBUG_HEADERS = (
    "seed",
    "i",
    "j",
    "t_i",
    "r_i",
    "theta_i",
    "phi_i",
    "t_j",
    "r_j",
    "theta_j",
    "phi_j",
    "delta_t",
    "angular_separation",
    "u1",
    "u2",
    "direct_branch_sign",
    "reason_code",
    "phi_target",
    "phi_min_sampled",
    "phi_max_sampled",
    "c2_at_phi_max",
    "c2_valid_sample_count",
    "c2_invalid_root_count",
    "bracket_failure_type",
)


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


def cubic_f(u: float, c2: float, mass: float = MASS) -> float:
    """He-Rideout cubic f(u)=2*M*u^3-u^2+c^2."""

    return 2.0 * mass * u * u * u - u * u + c2


def _bisect_root(
    left: float,
    right: float,
    c2: float,
    mass: float,
    eps: float,
    max_iter: int = 200,
) -> float:
    """Find a root of cubic_f in a bracket with opposite signs."""

    f_left = cubic_f(left, c2, mass)
    f_right = cubic_f(right, c2, mass)
    if abs(f_left) <= eps:
        return left
    if abs(f_right) <= eps:
        return right
    if f_left * f_right > 0.0:
        raise ValueError("root bracket does not change sign")

    lo = left
    hi = right
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        f_mid = cubic_f(mid, c2, mass)
        if abs(f_mid) <= eps or abs(hi - lo) <= eps:
            return mid
        if f_left * f_mid <= 0.0:
            hi = mid
            f_right = f_mid
            _ = f_right
        else:
            lo = mid
            f_left = f_mid
    return 0.5 * (lo + hi)


def positive_real_roots_cubic(c2: float, mass: float = MASS, eps: float = CUBIC_ROOT_EPS) -> list[float]:
    """Positive real roots of 2*M*u^3-u^2+c2 on u>0.

    This uses the Schwarzschild-specific cubic shape from He & Rideout:
    for c2 < 1/(27*M^2) there are two positive roots, at equality there
    is one double positive root, and above it there are no positive roots.
    The c2=0 root at u=0 is nonnegative but not positive, so it is omitted.
    """

    if mass <= 0.0:
        raise ValueError("mass must be positive")
    if c2 < -eps:
        raise ValueError("c2 must be nonnegative for the He-Rideout shooting parameter")

    c2 = max(0.0, c2)
    critical = 1.0 / (27.0 * mass * mass)
    u_min = 1.0 / (3.0 * mass)
    u_horizon = 1.0 / (2.0 * mass)

    if c2 > critical + eps:
        return []
    if abs(c2 - critical) <= eps:
        return [u_min]
    if abs(c2) <= eps:
        return [u_horizon]

    lower = _bisect_root(0.0, u_min, c2, mass, eps)
    upper = _bisect_root(u_min, u_horizon, c2, mass, eps)
    roots = [root for root in (lower, upper) if root > eps]
    roots.sort()
    return roots


def interval_has_root_between(
    u1: float,
    u2: float,
    roots: list[float],
    eps: float = CUBIC_ROOT_EPS,
) -> bool:
    """Return True if a root lies strictly inside the direct u interval."""

    lo = min(u1, u2)
    hi = max(u1, u2)
    return any(lo + eps < root < hi - eps for root in roots)


def direct_branch_sign(u1: float, u2: float) -> Optional[float]:
    """He-Rideout direct-branch sign for dphi/du and dt/du."""

    if abs(u2 - u1) <= CUBIC_ROOT_EPS:
        return None
    return 1.0 if u2 > u1 else -1.0


def shooting_candidate_valid(u1: float, u2: float, c2: float, mass: float = MASS) -> bool:
    """Whether c2 gives a real finite direct integrand on u1 -> u2."""

    roots = positive_real_roots_cubic(c2, mass=mass)
    if interval_has_root_between(u1, u2, roots):
        return False
    for u in (u1, u2, 0.5 * (u1 + u2)):
        if cubic_f(u, c2, mass) <= CUBIC_ROOT_EPS:
            return False
    return True


def _simpson_integral(func, a: float, b: float, intervals: int = SIMPSON_INTERVALS) -> Optional[float]:
    """Composite Simpson integration with finite-value checks."""

    if intervals <= 0:
        raise ValueError("intervals must be positive")
    if intervals % 2:
        intervals += 1
    if abs(b - a) <= CUBIC_ROOT_EPS:
        return 0.0

    h = (b - a) / intervals
    total = 0.0
    for idx in range(intervals + 1):
        u = a + idx * h
        value = func(u)
        if not math.isfinite(value):
            return None
        coeff = 1 if idx == 0 or idx == intervals else 4 if idx % 2 else 2
        total += coeff * value
    result = total * h / 3.0
    if not math.isfinite(result):
        return None
    return result


def direct_phi_integral(
    u1: float,
    u2: float,
    c2: float,
    mass: float = MASS,
    intervals: int = SIMPSON_INTERVALS,
) -> Optional[float]:
    """Integrate He-Rideout dphi/du on the direct branch."""

    sign = direct_branch_sign(u1, u2)
    if sign is None or not shooting_candidate_valid(u1, u2, c2, mass):
        return None

    def integrand(u: float) -> float:
        f_value = cubic_f(u, c2, mass)
        if f_value <= CUBIC_ROOT_EPS:
            return math.nan
        return sign / math.sqrt(f_value)

    value = _simpson_integral(integrand, u1, u2, intervals)
    if value is None:
        return None
    return abs(value)


def direct_time_integral(
    u1: float,
    u2: float,
    c2: float,
    mass: float = MASS,
    intervals: int = SIMPSON_INTERVALS,
) -> Optional[float]:
    """Integrate He-Rideout dt/du on the same direct branch."""

    sign = direct_branch_sign(u1, u2)
    if sign is None or not shooting_candidate_valid(u1, u2, c2, mass):
        return None
    c = math.sqrt(max(0.0, c2))

    def integrand(u: float) -> float:
        f_value = cubic_f(u, c2, mass)
        denom = u * u - 2.0 * mass * u * u * u
        if f_value <= CUBIC_ROOT_EPS or abs(denom) <= CUBIC_ROOT_EPS:
            return math.nan
        return (sign * c / math.sqrt(f_value) - 2.0 * mass * u) / denom

    value = _simpson_integral(integrand, u1, u2, intervals)
    if value is None:
        return None
    return value


def find_direct_shooting_c2(
    u1: float,
    u2: float,
    phi_target: float,
    mass: float = MASS,
) -> Optional[tuple[float, float]]:
    """Find c2 whose direct angular integral reaches phi_target."""

    if phi_target <= ANGLE_EPS or direct_branch_sign(u1, u2) is None:
        return None

    critical = 1.0 / (27.0 * mass * mass)
    low = critical * (1.0 + 1.0e-9)
    phi_low = direct_phi_integral(u1, u2, low, mass)
    if phi_low is None:
        return None

    if phi_low < phi_target:
        return None

    high = max(2.0 * low, low + 1.0e-6)
    phi_high = direct_phi_integral(u1, u2, high, mass)
    expand_count = 0
    while phi_high is not None and phi_high > phi_target and expand_count < 80:
        high *= 2.0
        phi_high = direct_phi_integral(u1, u2, high, mass)
        expand_count += 1
    if phi_high is None or phi_high > phi_target:
        return None

    best_c2 = low
    best_phi = phi_low
    for _ in range(SHOOTING_MAX_ITER):
        mid = 0.5 * (low + high)
        phi_mid = direct_phi_integral(u1, u2, mid, mass)
        if phi_mid is None:
            low = mid
            continue
        best_c2 = mid
        best_phi = phi_mid
        if abs(phi_mid - phi_target) <= SHOOTING_PHI_TOL:
            return best_c2, best_phi
        if phi_mid > phi_target:
            low = mid
        else:
            high = mid
    if abs(best_phi - phi_target) <= 10.0 * SHOOTING_PHI_TOL:
        return best_c2, best_phi
    return None


def find_direct_shooting_c2_with_reason(
    u1: float,
    u2: float,
    phi_target: float,
    mass: float = MASS,
) -> tuple[Optional[tuple[float, float]], str]:
    """Debug wrapper for c2 shooting with a coarse failure reason."""

    if phi_target <= ANGLE_EPS or direct_branch_sign(u1, u2) is None:
        return None, "same_u_generic"

    critical = 1.0 / (27.0 * mass * mass)
    low = critical * (1.0 + 1.0e-9)
    roots_low = positive_real_roots_cubic(low, mass=mass)
    if interval_has_root_between(u1, u2, roots_low):
        return None, "root_obstruction"

    phi_low = direct_phi_integral(u1, u2, low, mass)
    if phi_low is None:
        return None, "phi_integral_failed"
    if phi_low < phi_target:
        return None, "no_valid_c2_bracket"

    high = max(2.0 * low, low + 1.0e-6)
    phi_high = direct_phi_integral(u1, u2, high, mass)
    expand_count = 0
    while phi_high is not None and phi_high > phi_target and expand_count < 80:
        high *= 2.0
        phi_high = direct_phi_integral(u1, u2, high, mass)
        expand_count += 1
    if phi_high is None:
        return None, "phi_integral_failed"
    if phi_high > phi_target:
        return None, "no_valid_c2_bracket"

    best_c2 = low
    best_phi = phi_low
    saw_integral_failure = False
    for _ in range(SHOOTING_MAX_ITER):
        mid = 0.5 * (low + high)
        phi_mid = direct_phi_integral(u1, u2, mid, mass)
        if phi_mid is None:
            saw_integral_failure = True
            low = mid
            continue
        best_c2 = mid
        best_phi = phi_mid
        if abs(phi_mid - phi_target) <= SHOOTING_PHI_TOL:
            return (best_c2, best_phi), ""
        if phi_mid > phi_target:
            low = mid
        else:
            high = mid
    if abs(best_phi - phi_target) <= 10.0 * SHOOTING_PHI_TOL:
        return (best_c2, best_phi), ""
    if saw_integral_failure:
        return None, "phi_integral_failed"
    return None, "shooting_no_convergence"


def sample_direct_phi_range(
    u1: float,
    u2: float,
    phi_target: float,
    mass: float = MASS,
    sample_count: int = 96,
) -> dict[str, object]:
    """Sample the direct branch angular range over valid c2 candidates."""

    if direct_branch_sign(u1, u2) is None:
        return {
            "phi_target": phi_target,
            "phi_min_sampled": "NA",
            "phi_max_sampled": "NA",
            "c2_at_phi_max": "NA",
            "c2_valid_sample_count": 0,
            "c2_invalid_root_count": 0,
            "bracket_failure_type": "no_valid_samples",
        }

    critical = 1.0 / (27.0 * mass * mass)
    c2_min = critical * (1.0 + 1.0e-9)
    c2_max = max(1.0, 64.0 * c2_min)
    valid_values: list[tuple[float, float]] = []
    invalid_root_count = 0
    nan_count = 0

    for index in range(sample_count):
        if sample_count == 1:
            c2 = c2_min
        else:
            frac = index / (sample_count - 1)
            c2 = c2_min * (c2_max / c2_min) ** frac
        roots = positive_real_roots_cubic(c2, mass=mass)
        if interval_has_root_between(u1, u2, roots):
            invalid_root_count += 1
            continue
        if not shooting_candidate_valid(u1, u2, c2, mass):
            invalid_root_count += 1
            continue
        phi_value = direct_phi_integral(u1, u2, c2, mass)
        if phi_value is None or not math.isfinite(phi_value):
            nan_count += 1
            continue
        valid_values.append((c2, phi_value))

    if not valid_values:
        failure_type = "no_valid_samples" if nan_count == 0 else "nonmonotonic_or_nan"
        return {
            "phi_target": phi_target,
            "phi_min_sampled": "NA",
            "phi_max_sampled": "NA",
            "c2_at_phi_max": "NA",
            "c2_valid_sample_count": 0,
            "c2_invalid_root_count": invalid_root_count,
            "bracket_failure_type": failure_type,
        }

    phi_values = [value for _c2, value in valid_values]
    c2_at_phi_max, phi_max = max(valid_values, key=lambda item: item[1])
    phi_min = min(phi_values)
    if phi_target > phi_max + SHOOTING_PHI_TOL:
        failure_type = "target_above_direct_max"
    elif phi_target < phi_min - SHOOTING_PHI_TOL:
        failure_type = "target_below_direct_min"
    elif nan_count:
        failure_type = "nonmonotonic_or_nan"
    else:
        failure_type = "other"
    return {
        "phi_target": phi_target,
        "phi_min_sampled": phi_min,
        "phi_max_sampled": phi_max,
        "c2_at_phi_max": c2_at_phi_max,
        "c2_valid_sample_count": len(valid_values),
        "c2_invalid_root_count": invalid_root_count,
        "bracket_failure_type": failure_type,
    }


def causal_relation_schwarzschild_direct_shooting(p: Event, q: Event) -> Optional[bool]:
    """Try the direct generic He-Rideout shooting step for an undecided pair."""

    phi_target = angular_separation(p, q)
    if phi_target <= ANGLE_EPS:
        return None
    u1 = 1.0 / p.r
    u2 = 1.0 / q.r
    c2_phi = find_direct_shooting_c2(u1, u2, phi_target)
    if c2_phi is None:
        return None
    c2, _phi_value = c2_phi
    null_dt = direct_time_integral(u1, u2, c2)
    if null_dt is None or not math.isfinite(null_dt):
        return None
    event_dt = q.t - p.t
    return null_dt <= event_dt + TIME_EPS


def direct_shooting_failure_reason(p: Event, q: Event) -> str:
    """Coarse reason code for a generic pair left undecided by direct shooting."""

    phi_target = angular_separation(p, q)
    u1 = 1.0 / p.r
    u2 = 1.0 / q.r
    c2_phi, reason = find_direct_shooting_c2_with_reason(u1, u2, phi_target)
    if c2_phi is None:
        return reason or "other"
    c2, _phi_value = c2_phi
    null_dt = direct_time_integral(u1, u2, c2)
    if null_dt is None or not math.isfinite(null_dt):
        return "time_integral_failed"
    return "other"


def causal_relation_schwarzschild(p: Event, q: Event, enable_shooting: bool = False) -> Optional[bool]:
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

    if enable_shooting:
        return causal_relation_schwarzschild_direct_shooting(p, q)
    return None


def debug_undecided_row(seed: int, i: int, j: int, p: Event, q: Event, reason: str) -> dict[str, object]:
    """Build one debug row for a pair left undecided."""

    u1 = 1.0 / p.r
    u2 = 1.0 / q.r
    sign = direct_branch_sign(u1, u2)
    phi_target = angular_separation(p, q)
    phi_range = sample_direct_phi_range(u1, u2, phi_target)
    return {
        "seed": seed,
        "i": i,
        "j": j,
        "t_i": p.t,
        "r_i": p.r,
        "theta_i": p.theta,
        "phi_i": p.phi,
        "t_j": q.t,
        "r_j": q.r,
        "theta_j": q.theta,
        "phi_j": q.phi,
        "delta_t": q.t - p.t,
        "angular_separation": phi_target,
        "u1": u1,
        "u2": u2,
        "direct_branch_sign": "None" if sign is None else sign,
        "reason_code": reason,
        **phi_range,
    }


def build_causal_matrix(
    events: list[Event],
    seed: int,
    enable_shooting: bool = False,
    debug_undecided: bool = False,
) -> tuple[vs.CausalMatrix, list[list[Optional[bool]]], list[dict[str, object]]]:
    """Build the asserted causal matrix from implemented Schwarzschild tests."""

    n = len(events)
    matrix: vs.CausalMatrix = [[False] * n for _ in range(n)]
    states: list[list[Optional[bool]]] = [[False if i == j else None for j in range(n)] for i in range(n)]
    debug_rows: list[dict[str, object]] = []
    for i in range(n - 1):
        for j in range(i + 1, n):
            relation = causal_relation_schwarzschild(events[i], events[j], enable_shooting=enable_shooting)
            states[i][j] = relation
            if relation:
                matrix[i][j] = True
            elif relation is None and debug_undecided:
                reason = direct_shooting_failure_reason(events[i], events[j]) if enable_shooting else "other"
                debug_rows.append(debug_undecided_row(seed, i, j, events[i], events[j], reason))
    return matrix, states, debug_rows


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


def summarize_case(
    events: list[Event],
    matrix: vs.CausalMatrix,
    states: list[list[Optional[bool]]],
    n: int,
    seed: int,
    enable_shooting: bool,
) -> dict[str, object]:
    _ = events
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
        "seed": seed,
        "schwarzschild_radius": SCHWARZSCHILD_RADIUS,
        "mass": MASS,
        "exterior_margin": EXTERIOR_MARGIN,
        "r_min": R_MIN,
        "r_max": R_MAX,
        "t_min": T_MIN,
        "t_max": T_MAX,
        "coordinate_patch": "exterior Eddington-Finkelstein coordinates",
        "causal_relation_model": (
            "He & Rideout radial exact tests plus sufficient bounds and direct generic shooting"
            if enable_shooting
            else "He & Rideout radial exact tests plus sufficient bounds; generic geodesic shooting TODO"
        ),
        "direct_shooting_enabled": enable_shooting,
        "true_relations": true_relations,
        "false_relations": false_relations,
        "decided_pairs": decided_pairs,
        "relations": true_relations,
        "ordering_fraction": ordering_fraction,
        "ordering_fraction_decided": ordering_fraction,
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
    return summary


def run_case(
    n: int,
    seed: int,
    enable_shooting: bool = False,
    debug_undecided: bool = False,
) -> tuple[list[Event], vs.CausalMatrix, list[list[Optional[bool]]], dict[str, object], list[dict[str, object]]]:
    events = generate_exterior_events(n, seed)
    matrix, states, debug_rows = build_causal_matrix(
        events,
        seed=seed,
        enable_shooting=enable_shooting,
        debug_undecided=debug_undecided,
    )
    summary = summarize_case(events, matrix, states, n=n, seed=seed, enable_shooting=enable_shooting)
    return events, matrix, states, summary, debug_rows


def write_single_outputs(
    events: list[Event],
    matrix: vs.CausalMatrix,
    states: list[list[Optional[bool]]],
    summary: dict[str, object],
    out_prefix: str,
    debug_rows: list[dict[str, object]] | None = None,
) -> tuple[Path, Path]:
    csv_path = OUT_DIR / f"{out_prefix}.csv"
    json_path = OUT_DIR / f"{out_prefix}.json"
    links = transitive_reduction_links(matrix)

    headers = list(summary)
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(headers)
        writer.writerow([_fmt(summary[header]) for header in headers])

    payload = {
        "summary": summary,
        "events": [asdict(event) for event in events],
        "causal_matrix": matrix,
        "relation_states": states,
        "links": links,
        "undecided_debug": debug_rows or [],
        "notes": [
            "S4-2 asserts only relations decided by He & Rideout radial exact tests or sufficient bounds.",
            "False relation_states entries are decided non-relations; null entries are generic undecided pairs.",
            "False causal_matrix entries include both decided false and undecided pairs; inspect relation_states for the distinction.",
            "The next required step is He & Rideout generic null-geodesic shooting.",
        ],
    }
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return csv_path, json_path


def write_debug_csv(debug_rows: list[dict[str, object]], out_prefix: str) -> Optional[Path]:
    """Write undecided-pair debug rows when requested."""

    if not debug_rows:
        return None
    debug_path = OUT_DIR / f"{out_prefix}_undecided_debug.csv"
    with debug_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(DEBUG_HEADERS)
        for row in debug_rows:
            writer.writerow([_fmt(row[header]) for header in DEBUG_HEADERS])
    return debug_path


def _min_mean_max(rows: list[dict[str, object]], key: str) -> dict[str, float]:
    values = [float(row[key]) for row in rows]
    return {
        "min": min(values),
        "mean": sum(values) / len(values),
        "max": max(values),
    }


def failed_check_count(rows: list[dict[str, object]]) -> int:
    return sum(
        1
        for row in rows
        if not (
            bool(row["antisymmetric"])
            and bool(row["transitive_true_matrix"])
            and bool(row["decided_transitivity_no_false_contradictions"])
        )
    )


def write_sweep_outputs(
    rows: list[dict[str, object]],
    out_prefix: str,
    debug_rows: list[dict[str, object]] | None = None,
) -> tuple[Path, Path, dict[str, object], Optional[Path]]:
    csv_path = OUT_DIR / f"{out_prefix}.csv"
    json_path = OUT_DIR / f"{out_prefix}.json"

    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(ROW_HEADERS)
        for row in rows:
            writer.writerow([_fmt(row[header]) for header in ROW_HEADERS])

    aggregate = {
        "benchmark": "Schwarzschild benchmark before Kerr ordinal diagnostics.",
        "status": "partial_model_seed_sweep_generic_pairs_undecided",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "command": COMMAND,
        "num_rows": len(rows),
        "N_values": sorted({int(row["N"]) for row in rows}),
        "seed_min": min(int(row["seed"]) for row in rows),
        "seed_max": max(int(row["seed"]) for row in rows),
        "undecided_pairs": _min_mean_max(rows, "undecided_pairs"),
        "true_relations": _min_mean_max(rows, "true_relations"),
        "ordering_fraction_decided": _min_mean_max(rows, "ordering_fraction_decided"),
        "failed_checks": failed_check_count(rows),
        "notes": [
            "This sweep exercises only the partial He & Rideout radial/bound model.",
            "It does not validate the full Schwarzschild causal relation.",
            "Undecided pairs remain explicit and require generic null-geodesic shooting.",
        ],
    }
    json_path.write_text(
        json.dumps({"aggregate": aggregate, "rows": rows, "undecided_debug": debug_rows or []}, indent=2),
        encoding="utf-8",
    )
    debug_path = write_debug_csv(debug_rows or [], out_prefix)
    return csv_path, json_path, aggregate, debug_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the partial He-Rideout Schwarzschild exterior benchmark."
    )
    parser.add_argument("--N", type=int, default=DEFAULT_N, help="number of exterior events")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help="initial random seed")
    parser.add_argument("--num-seeds", type=int, default=1, help="number of consecutive seeds to run")
    parser.add_argument(
        "--out-prefix",
        default=DEFAULT_OUT_PREFIX,
        help="output filename prefix under explore/sorkin4_schwarzschild_benchmark/",
    )
    parser.add_argument(
        "--self-test-roots",
        action="store_true",
        help="run internal checks for the He-Rideout cubic root helpers",
    )
    parser.add_argument(
        "--enable-shooting",
        action="store_true",
        help="try direct generic He-Rideout shooting for pairs left undecided by bounds",
    )
    parser.add_argument(
        "--debug-undecided",
        action="store_true",
        help="write and print diagnostic rows for pairs left undecided",
    )
    args = parser.parse_args()
    if args.N <= 0:
        raise SystemExit("--N must be positive")
    if args.num_seeds <= 0:
        raise SystemExit("--num-seeds must be positive")
    return args


def _assert_close(actual: float, expected: float, tol: float, label: str) -> None:
    if abs(actual - expected) > tol:
        raise AssertionError(f"{label}: got {actual:.17g}, expected {expected:.17g}")


def self_test_roots() -> None:
    """Small internal checks for the He-Rideout cubic root helpers."""

    mass = 1.0
    critical = 1.0 / (27.0 * mass * mass)

    roots_zero = positive_real_roots_cubic(0.0, mass=mass)
    if len(roots_zero) != 1:
        raise AssertionError(f"c2=0 expected one positive root, got {roots_zero}")
    _assert_close(roots_zero[0], 0.5, 1.0e-10, "c2=0 horizon root")
    _assert_close(cubic_f(roots_zero[0], 0.0, mass), 0.0, 1.0e-10, "c2=0 f(root)")

    c2_small = 0.01
    roots_small = positive_real_roots_cubic(c2_small, mass=mass)
    if len(roots_small) != 2:
        raise AssertionError(f"small c2 expected two positive roots, got {roots_small}")
    if not roots_small[0] < 1.0 / 3.0 < roots_small[1] < 0.5:
        raise AssertionError(f"small c2 roots in unexpected locations: {roots_small}")
    for index, root in enumerate(roots_small):
        _assert_close(cubic_f(root, c2_small, mass), 0.0, 1.0e-10, f"small c2 f(root {index})")

    roots_critical = positive_real_roots_cubic(critical, mass=mass)
    if len(roots_critical) != 1:
        raise AssertionError(f"critical c2 expected one double positive root, got {roots_critical}")
    _assert_close(roots_critical[0], 1.0 / 3.0, 1.0e-10, "critical double root")

    roots_large = positive_real_roots_cubic(0.1, mass=mass)
    if roots_large:
        raise AssertionError(f"large c2 expected no positive roots, got {roots_large}")
    if interval_has_root_between(0.1, 0.4, roots_large):
        raise AssertionError("large c2 interval unexpectedly reported a root")

    forward = interval_has_root_between(0.05, 0.2, roots_small)
    backward = interval_has_root_between(0.2, 0.05, roots_small)
    if not (forward and backward):
        raise AssertionError("interval root check should be symmetric under endpoint order")

    edge = interval_has_root_between(roots_small[0], 0.2, roots_small)
    if edge:
        raise AssertionError("root exactly at interval endpoint should not count as strictly between")

    print("S4 root self-test")
    print(f"c2=0 positive_roots={[_fmt(root) for root in roots_zero]}")
    print(f"c2=0.01 positive_roots={[_fmt(root) for root in roots_small]}")
    print(f"c2=1/27 positive_roots={[_fmt(root) for root in roots_critical]}")
    print(f"c2=0.1 positive_roots={[_fmt(root) for root in roots_large]}")
    print(f"interval_symmetry={forward == backward == True}")
    print("status=pass")


def print_debug_rows(debug_rows: list[dict[str, object]], debug_path: Optional[Path]) -> None:
    """Print a compact undecided-pair debug table."""

    if not debug_rows:
        print("undecided_debug_rows=0")
        return
    print(f"undecided_debug_rows={len(debug_rows)}")
    if debug_path is not None:
        print(f"wrote {debug_path}")
    print(",".join(DEBUG_HEADERS))
    for row in debug_rows:
        print(",".join(_fmt(row[header]) for header in DEBUG_HEADERS))


def print_single_summary(
    summary: dict[str, object],
    csv_path: Path,
    json_path: Path,
    debug_rows: list[dict[str, object]] | None = None,
    debug_path: Optional[Path] = None,
) -> None:
    print("S4-2b Schwarzschild minimal benchmark")
    print("Regime: exterior only, r > r_s + margin")
    print(f"N={summary['N']} seed={summary['seed']}")
    print(f"M={summary['mass']} r_s={summary['schwarzschild_radius']} margin={summary['exterior_margin']}")
    print(
        f"true_relations={summary['true_relations']} "
        f"false_relations={summary['false_relations']} "
        f"undecided_pairs={summary['undecided_pairs']}"
    )
    print(f"ordering_fraction_decided={summary['ordering_fraction_decided']:.6g}")
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
    print(f"wrote {csv_path}")
    print(f"wrote {json_path}")
    if debug_rows is not None:
        print_debug_rows(debug_rows, debug_path)


def print_sweep_summary(
    rows: list[dict[str, object]],
    aggregate: dict[str, object],
    csv_path: Path,
    json_path: Path,
    debug_rows: list[dict[str, object]] | None = None,
    debug_path: Optional[Path] = None,
) -> None:
    print("S4-2b Schwarzschild partial-model robustness sweep")
    print("Regime: exterior only, r > r_s + margin")
    print(
        f"N={rows[0]['N']} seeds={aggregate['seed_min']}..{aggregate['seed_max']} "
        f"num_seeds={aggregate['num_rows']}"
    )
    print(
        "undecided_pairs="
        f"min/mean/max {aggregate['undecided_pairs']['min']:.6g}/"
        f"{aggregate['undecided_pairs']['mean']:.6g}/"
        f"{aggregate['undecided_pairs']['max']:.6g}"
    )
    print(
        "true_relations="
        f"min/mean/max {aggregate['true_relations']['min']:.6g}/"
        f"{aggregate['true_relations']['mean']:.6g}/"
        f"{aggregate['true_relations']['max']:.6g}"
    )
    print(
        "ordering_fraction_decided="
        f"min/mean/max {aggregate['ordering_fraction_decided']['min']:.6g}/"
        f"{aggregate['ordering_fraction_decided']['mean']:.6g}/"
        f"{aggregate['ordering_fraction_decided']['max']:.6g}"
    )
    print(f"failed_checks={aggregate['failed_checks']}")
    if aggregate["undecided_pairs"]["max"] > 0:
        print("warning=undecided generic pairs remain; full He-Rideout shooting is still required")
    print(f"wrote {csv_path}")
    print(f"wrote {json_path}")
    if debug_rows is not None:
        print_debug_rows(debug_rows, debug_path)


def main() -> None:
    args = parse_args()
    if args.self_test_roots:
        self_test_roots()
        return

    if args.num_seeds == 1:
        events, matrix, states, summary, debug_rows = run_case(
            args.N,
            args.seed,
            enable_shooting=args.enable_shooting,
            debug_undecided=args.debug_undecided,
        )
        csv_path, json_path = write_single_outputs(
            events,
            matrix,
            states,
            summary,
            args.out_prefix,
            debug_rows=debug_rows if args.debug_undecided else None,
        )
        debug_path = write_debug_csv(debug_rows, args.out_prefix) if args.debug_undecided else None
        print_single_summary(
            summary,
            csv_path,
            json_path,
            debug_rows=debug_rows if args.debug_undecided else None,
            debug_path=debug_path,
        )
        return

    rows: list[dict[str, object]] = []
    debug_rows_all: list[dict[str, object]] = []
    for seed in range(args.seed, args.seed + args.num_seeds):
        _events, _matrix, _states, summary, debug_rows = run_case(
            args.N,
            seed,
            enable_shooting=args.enable_shooting,
            debug_undecided=args.debug_undecided,
        )
        rows.append({header: summary[header] for header in ROW_HEADERS})
        debug_rows_all.extend(debug_rows)
    csv_path, json_path, aggregate, debug_path = write_sweep_outputs(
        rows,
        args.out_prefix,
        debug_rows=debug_rows_all if args.debug_undecided else None,
    )
    print_sweep_summary(
        rows,
        aggregate,
        csv_path,
        json_path,
        debug_rows=debug_rows_all if args.debug_undecided else None,
        debug_path=debug_path,
    )


if __name__ == "__main__":
    main()
