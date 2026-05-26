#!/usr/bin/env python3
"""Minimal Kerr exterior scaffold for SORKIN-4.

This is not a Kerr causal-relation solver.  For a=0 it regresses to the
implemented Schwarzschild exterior benchmark subset.  For a!=0 it only
generates exterior Boyer-Lindquist-like coordinates and marks all pairs
undecided until a justified Kerr causal criterion is added.
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


OUT_DIR = Path(__file__).resolve().parent
COMMAND = "python3 explore/sorkin4_kerr_benchmark/run_kerr_minimal_benchmark.py"

DEFAULT_N = 12
DEFAULT_SEED = 1959
DEFAULT_M = 1.0
DEFAULT_A = 0.0
DEFAULT_R_MIN_MARGIN = 0.35
DEFAULT_OUT_PREFIX = "kerr_minimal_benchmark"
DEFAULT_KERR_MODE = "scaffold"
DEFAULT_NUM_SEEDS = 1

R_MAX = 6.0
T_MIN = 0.0
T_MAX = 4.0
ANGLE_EPS = 1.0e-12
TIME_EPS = 1.0e-12
LOCAL_CONE_MARGIN = 0.75
LOCAL_CONE_EPS = 1.0e-9
RADIAL_PHI_EPS = 1.0e-9
KERR_MODES = ("scaffold", "local_cone_diagnostic", "equatorial_scaffold")
CALIBRATION_HEADERS = (
    "seed",
    "i",
    "j",
    "schwarzschild_relation",
    "local_relation",
    "agree",
    "delta_t",
    "r_i",
    "r_j",
    "theta_i",
    "theta_j",
    "phi_i",
    "phi_j",
    "ds2_midpoint",
)


@dataclass(frozen=True)
class Event:
    """Known-coordinate event in an exterior Boyer-Lindquist-like chart."""

    index: int
    t: float
    r: float
    theta: float
    phi: float


def kerr_horizon_radius(mass: float, spin: float) -> float:
    """Outer Kerr horizon radius r_+=M+sqrt(M^2-a^2)."""

    if mass <= 0.0:
        raise ValueError("M must be positive")
    if abs(spin) >= mass:
        raise ValueError("|a| must be < M for this exterior scaffold")
    return mass + math.sqrt(mass * mass - spin * spin)


def sample_radius_volume_weighted(rng: random.Random, r_min: float) -> float:
    """Sample r with density proportional to r^2 on [r_min, R_MAX]."""

    if r_min >= R_MAX:
        raise ValueError("r_min must be smaller than R_MAX")
    u = rng.random()
    return (r_min**3 + u * (R_MAX**3 - r_min**3)) ** (1.0 / 3.0)


def generate_exterior_events(n: int, seed: int, r_min: float, equatorial: bool = False) -> list[Event]:
    """Generate a reproducible bounded exterior point set."""

    rng = random.Random(seed)
    events: list[Event] = []
    for index in range(n):
        t = rng.uniform(T_MIN, T_MAX)
        r = sample_radius_volume_weighted(rng, r_min)
        cos_theta = rng.uniform(-1.0, 1.0)
        theta = math.pi / 2.0 if equatorial else math.acos(cos_theta)
        phi = rng.uniform(0.0, 2.0 * math.pi)
        events.append(Event(index=index, t=t, r=r, theta=theta, phi=phi))
    return sorted(events, key=lambda event: (event.t, event.r, event.theta, event.phi))


def angular_separation(p: Event, q: Event) -> float:
    """Great-circle angular separation used by the Schwarzschild regression path."""

    cos_angle = (
        math.cos(p.theta) * math.cos(q.theta)
        + math.sin(p.theta) * math.sin(q.theta) * math.cos(p.phi - q.phi)
    )
    return math.acos(max(-1.0, min(1.0, cos_angle)))


def angular_delta(phi1: float, phi2: float) -> float:
    """Shortest signed azimuthal coordinate difference."""

    return signed_delta_phi(phi1, phi2)


def signed_delta_phi(phi_i: float, phi_j: float) -> float:
    """Shortest signed azimuthal difference in (-pi, pi]."""

    delta = (phi_j - phi_i + math.pi) % (2.0 * math.pi) - math.pi
    if delta <= -math.pi + ANGLE_EPS:
        return math.pi
    return delta


def outgoing_radial_trip(r1: float, r2: float, mass: float) -> float:
    """EF-time for an outgoing radial Schwarzschild null ray outside the horizon."""

    horizon = 2.0 * mass
    return r2 - r1 + 4.0 * mass * math.log((r2 - horizon) / (r1 - horizon))


def angular_f(r: float, mass: float) -> float:
    """He-Rideout f(r)=r/sqrt(1-2M/r), valid outside the Schwarzschild horizon."""

    return r / math.sqrt(1.0 - (2.0 * mass) / r)


def angular_spacelike_r0(r1: float, r2: float, mass: float) -> float:
    """Minimizer of f(r) on the outgoing exterior interval [r1, r2]."""

    photon_radius = 3.0 * mass
    if photon_radius <= r1 <= r2:
        return r1
    if r1 < photon_radius < r2:
        return photon_radius
    return r2


def timelike_bound_r0(r1: float, r2: float, mass: float) -> float:
    """He-Rideout r0 for the composed null-curve sufficient timelike bound."""

    photon_radius = 3.0 * mass
    if r1 >= photon_radius and r2 >= photon_radius:
        return min(r1, r2)
    if (r1 > photon_radius > r2) or (r2 > photon_radius > r1):
        return photon_radius
    return max(r1, r2)


def causal_relation_schwarzschild(p: Event, q: Event, mass: float) -> Optional[bool]:
    """Return the implemented Schwarzschild benchmark subset for p before q.

    This is the same default subset as the Schwarzschild benchmark: radial
    exact tests plus sufficient bounds.  Generic shooting is still TODO.
    """

    dt = q.t - p.t
    if dt < -TIME_EPS:
        return False

    r1 = p.r
    r2 = q.r
    phi2 = angular_separation(p, q)
    horizon = 2.0 * mass

    if r2 > r1 and r1 < horizon:
        return False

    if phi2 <= ANGLE_EPS:
        if r1 >= r2:
            return dt + TIME_EPS >= r1 - r2
        if r2 >= r1 > horizon:
            return dt + TIME_EPS >= outgoing_radial_trip(r1, r2, mass)
        return False

    if r1 >= r2:
        radial_lower_bound = r1 - r2
        if dt < radial_lower_bound - TIME_EPS:
            return False
        angular_lower_bound = r2 * phi2
        if dt < angular_lower_bound - TIME_EPS:
            return False
    elif r2 >= r1 > horizon:
        radial_lower_bound = outgoing_radial_trip(r1, r2, mass)
        if dt < radial_lower_bound - TIME_EPS:
            return False
        r0_spacelike = angular_spacelike_r0(r1, r2, mass)
        angular_lower_bound = angular_f(r0_spacelike, mass) * phi2
        if dt < angular_lower_bound - TIME_EPS:
            return False

    r0_timelike = timelike_bound_r0(r1, r2, mass)
    angular_trip = angular_f(r0_timelike, mass) * phi2
    if r1 >= r2 and r1 > horizon:
        composed_trip = r1 - r2 + angular_trip
        if dt + TIME_EPS >= composed_trip:
            return True
    elif r2 >= r1 > horizon:
        composed_trip = outgoing_radial_trip(r1, r2, mass) + angular_trip
        if dt + TIME_EPS >= composed_trip:
            return True

    return None


def kerr_metric_components(r: float, theta: float, mass: float, spin: float) -> dict[str, float]:
    """Kerr metric in Boyer-Lindquist coordinates, signature -+++."""

    sin_theta = math.sin(theta)
    cos_theta = math.cos(theta)
    sigma = r * r + spin * spin * cos_theta * cos_theta
    delta = r * r - 2.0 * mass * r + spin * spin
    if sigma <= 0.0 or delta <= 0.0:
        raise ValueError("metric point must be outside the Kerr horizon")
    return {
        "tt": -(1.0 - (2.0 * mass * r) / sigma),
        "tphi": -(2.0 * mass * spin * r * sin_theta * sin_theta) / sigma,
        "rr": sigma / delta,
        "thetatheta": sigma,
        "phiphi": (
            r * r
            + spin * spin
            + (2.0 * mass * spin * spin * r * sin_theta * sin_theta) / sigma
        )
        * sin_theta
        * sin_theta,
    }


def local_cone_ds2_midpoint(p: Event, q: Event, mass: float, spin: float) -> tuple[float, float]:
    """Return midpoint ds^2 and local spatial distance in BL coordinates."""

    dt = q.t - p.t
    r_mid = 0.5 * (p.r + q.r)
    theta_mid = 0.5 * (p.theta + q.theta)
    dr = q.r - p.r
    dtheta = q.theta - p.theta
    dphi = angular_delta(p.phi, q.phi)
    metric = kerr_metric_components(r_mid, theta_mid, mass, spin)

    spatial_sq = (
        metric["rr"] * dr * dr
        + metric["thetatheta"] * dtheta * dtheta
        + metric["phiphi"] * dphi * dphi
    )
    if spatial_sq < 0.0:
        raise ValueError("local spatial distance became negative")
    spatial_distance = math.sqrt(spatial_sq)

    ds2 = (
        metric["tt"] * dt * dt
        + 2.0 * metric["tphi"] * dt * dphi
        + metric["rr"] * dr * dr
        + metric["thetatheta"] * dtheta * dtheta
        + metric["phiphi"] * dphi * dphi
    )
    return ds2, spatial_distance


def local_cone_relation_kerr(p: Event, q: Event, mass: float, spin: float) -> Optional[bool]:
    """Local Kerr cone diagnostic from the BL metric at the pair midpoint.

    This is only an infinitesimal-cone proxy for nearby pairs.  It is not a
    global Kerr geodesic causal relation.
    """

    dt = q.t - p.t
    if dt <= TIME_EPS:
        return None

    r_mid = 0.5 * (p.r + q.r)
    ds2, spatial_distance = local_cone_ds2_midpoint(p, q, mass, spin)
    if spatial_distance > LOCAL_CONE_MARGIN * max(r_mid, mass):
        return None

    metric = kerr_metric_components(r_mid, 0.5 * (p.theta + q.theta), mass, spin)
    dphi = angular_delta(p.phi, q.phi)
    spatial_sq = spatial_distance * spatial_distance
    scale = max(1.0, abs(metric["tt"] * dt * dt), spatial_sq, abs(2.0 * metric["tphi"] * dt * dphi))
    eps = LOCAL_CONE_EPS * scale
    if ds2 < -eps:
        return True
    if ds2 > eps:
        return False
    return None


def causal_relation_model(
    p: Event,
    q: Event,
    mass: float,
    spin: float,
    kerr_mode: str,
) -> Optional[bool]:
    """Dispatch the current causal relation scaffold."""

    if abs(spin) <= 0.0:
        return causal_relation_schwarzschild(p, q, mass)
    if kerr_mode == "local_cone_diagnostic":
        return None
    if kerr_mode == "equatorial_scaffold":
        return None
    # TODO(S4-K): add a justified Kerr exterior causal criterion.
    return None


def count_local_cone_candidates(events: list[Event], mass: float, spin: float) -> tuple[int, int]:
    """Count local cone signs without using them as causal decisions."""

    timelike = 0
    spacelike = 0
    for i in range(len(events) - 1):
        for j in range(i + 1, len(events)):
            relation = local_cone_relation_kerr(events[i], events[j], mass, spin)
            if relation is True:
                timelike += 1
            elif relation is False:
                spacelike += 1
    return timelike, spacelike


def equatorial_pair_diagnostics(events: list[Event], spin: float) -> dict[str, float | int]:
    """Classify equatorial pair directions without deciding causality."""

    prograde_pairs = 0
    retrograde_pairs = 0
    radial_ish_pairs = 0
    abs_delta_phi_total = 0.0
    delta_t_total = 0.0
    pair_count = 0

    for i in range(len(events) - 1):
        for j in range(i + 1, len(events)):
            delta_phi = signed_delta_phi(events[i].phi, events[j].phi)
            abs_delta_phi = abs(delta_phi)
            abs_delta_phi_total += abs_delta_phi
            delta_t_total += events[j].t - events[i].t
            pair_count += 1
            if abs_delta_phi <= RADIAL_PHI_EPS:
                radial_ish_pairs += 1
            elif abs(spin) > 0.0 and delta_phi * spin > 0.0:
                prograde_pairs += 1
            elif abs(spin) > 0.0:
                retrograde_pairs += 1

    return {
        "prograde_pairs": prograde_pairs,
        "retrograde_pairs": retrograde_pairs,
        "radial_ish_pairs": radial_ish_pairs,
        "mean_abs_delta_phi": abs_delta_phi_total / pair_count if pair_count else 0.0,
        "mean_delta_t": delta_t_total / pair_count if pair_count else 0.0,
        "r_min_observed": min((event.r for event in events), default=0.0),
        "r_max_observed": max((event.r for event in events), default=0.0),
    }


def build_relation_states(
    events: list[Event],
    mass: float,
    spin: float,
    kerr_mode: str,
) -> tuple[list[list[bool]], list[list[Optional[bool]]]]:
    """Build asserted true matrix and explicit relation states."""

    n = len(events)
    matrix = [[False] * n for _ in range(n)]
    states: list[list[Optional[bool]]] = [[False if i == j else None for j in range(n)] for i in range(n)]
    for i in range(n - 1):
        for j in range(i + 1, n):
            relation = causal_relation_model(events[i], events[j], mass, spin, kerr_mode)
            states[i][j] = relation
            if relation is True:
                matrix[i][j] = True
    return matrix, states


def count_true_relations(matrix: list[list[bool]]) -> int:
    return sum(1 for i, row in enumerate(matrix) for j, value in enumerate(row) if i < j and value)


def check_antisymmetric(matrix: list[list[bool]]) -> bool:
    n = len(matrix)
    for i in range(n):
        if matrix[i][i]:
            return False
        for j in range(i + 1, n):
            if matrix[i][j] and matrix[j][i]:
                return False
    return True


def check_transitive(matrix: list[list[bool]]) -> bool:
    n = len(matrix)
    for i in range(n):
        for j in range(n):
            if not matrix[i][j]:
                continue
            for k in range(n):
                if matrix[j][k] and not matrix[i][k]:
                    return False
    return True


def failed_checks_for(matrix: list[list[bool]]) -> list[str]:
    failed: list[str] = []
    if not check_antisymmetric(matrix):
        failed.append("antisymmetric")
    if not check_transitive(matrix):
        failed.append("transitive_true_matrix")
    return failed


def model_label(spin: float, kerr_mode: str) -> str:
    if abs(spin) <= 0.0:
        return "Schwarzschild a=0 regression: He & Rideout radial exact tests plus sufficient bounds; generic geodesic shooting TODO"
    if kerr_mode == "local_cone_diagnostic":
        return "Kerr local cone diagnostic, not global geodesic causal relation"
    if kerr_mode == "equatorial_scaffold":
        return "Kerr equatorial scaffold: prograde/retrograde diagnostics only; no global geodesic causal relation"
    return "Kerr exterior scaffold only: causal relations undecided; Kerr null geodesic criterion TODO"


def summarize_case(
    events: list[Event],
    matrix: list[list[bool]],
    states: list[list[Optional[bool]]],
    n: int,
    seed: int,
    mass: float,
    spin: float,
    r_plus: float,
    r_min: float,
    margin: float,
    kerr_mode: str,
) -> dict[str, object]:
    possible_pairs = n * (n - 1) // 2
    true_relations = count_true_relations(matrix)
    false_relations = sum(1 for i in range(n - 1) for j in range(i + 1, n) if states[i][j] is False)
    undecided_pairs = sum(1 for i in range(n - 1) for j in range(i + 1, n) if states[i][j] is None)
    decided_pairs = true_relations + false_relations
    ordering_fraction = true_relations / decided_pairs if decided_pairs else 0.0
    failed_checks = failed_checks_for(matrix)
    local_diagnostic_active = abs(spin) > 0.0 and kerr_mode == "local_cone_diagnostic"
    equatorial_active = kerr_mode == "equatorial_scaffold"
    local_timelike_candidates = 0
    local_spacelike_candidates = 0
    if local_diagnostic_active:
        local_timelike_candidates, local_spacelike_candidates = count_local_cone_candidates(events, mass, spin)
    equatorial = equatorial_pair_diagnostics(events, spin) if equatorial_active else {
        "prograde_pairs": 0,
        "retrograde_pairs": 0,
        "radial_ish_pairs": 0,
        "mean_abs_delta_phi": 0.0,
        "mean_delta_t": 0.0,
        "r_min_observed": min((event.r for event in events), default=0.0),
        "r_max_observed": max((event.r for event in events), default=0.0),
    }

    return {
        "benchmark": "S4-K3 Kerr exterior minimal diagnostic scaffold",
        "status": (
            "a0_schwarzschild_regression"
            if abs(spin) <= 0.0
            else "kerr_equatorial_scaffold"
            if equatorial_active
            else "kerr_local_cone_diagnostic"
            if local_diagnostic_active
            else "kerr_scaffold_pairs_undecided"
        ),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "command": COMMAND,
        "N": n,
        "seed": seed,
        "M": mass,
        "a": spin,
        "kerr_mode": kerr_mode,
        "r_plus": r_plus,
        "r_min_margin": margin,
        "r_min": r_min,
        "r_max": R_MAX,
        "t_min": T_MIN,
        "t_max": T_MAX,
        "coordinate_patch": "Boyer-Lindquist-like exterior coordinates for generation only",
        "true_relations": true_relations,
        "false_relations": false_relations,
        "local_true_relations": true_relations if local_diagnostic_active else 0,
        "local_false_relations": false_relations if local_diagnostic_active else 0,
        "local_timelike_candidates": local_timelike_candidates,
        "local_spacelike_candidates": local_spacelike_candidates,
        "local_diagnostic_only": local_diagnostic_active,
        "prograde_pairs": equatorial["prograde_pairs"],
        "retrograde_pairs": equatorial["retrograde_pairs"],
        "radial_ish_pairs": equatorial["radial_ish_pairs"],
        "mean_abs_delta_phi": equatorial["mean_abs_delta_phi"],
        "mean_delta_t": equatorial["mean_delta_t"],
        "r_min_observed": equatorial["r_min_observed"],
        "r_max_observed": equatorial["r_max_observed"],
        "undecided_pairs": undecided_pairs,
        "decided_pairs": decided_pairs,
        "ordering_fraction_decided": ordering_fraction,
        "local_cone_margin_used": LOCAL_CONE_MARGIN if local_diagnostic_active else 0.0,
        "antisymmetric": check_antisymmetric(matrix),
        "transitive_true_matrix": check_transitive(matrix),
        "failed_checks": failed_checks,
        "possible_pairs": possible_pairs,
        "warning": (
            "local cone diagnostic is not used as causal decision; failed a=0 calibration as False filter"
            if local_diagnostic_active
            else "equatorial scaffold only; no Kerr geodesic causal decisions implemented"
            if equatorial_active and abs(spin) > 0.0
            else ""
        ),
        "causal_relation_model": model_label(spin, kerr_mode),
    }


def _fmt(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return f"{value:.12g}"
    if isinstance(value, list):
        return "[]" if not value else ";".join(str(item) for item in value)
    return str(value)


def write_outputs(
    events: list[Event],
    matrix: list[list[bool]],
    states: list[list[Optional[bool]]],
    summary: dict[str, object],
    out_prefix: str,
) -> tuple[Path, Path]:
    csv_path = OUT_DIR / f"{out_prefix}.csv"
    json_path = OUT_DIR / f"{out_prefix}.json"

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
        "notes": [
            "S4-K2 is a diagnostic scaffold, not a Kerr causal reconstruction.",
            "For a=0, the relation model is the Schwarzschild benchmark subset.",
            "For scaffold mode, a!=0 null relation states mean undecided pairs, not non-relations.",
            "For local_cone_diagnostic mode, local metric-cone signs are counted but not used as relation decisions.",
            "For equatorial_scaffold mode, theta is fixed to pi/2 and prograde/retrograde counts are diagnostic only.",
            "TODO(S4-K): justify and implement global Kerr exterior causal decisions.",
        ],
    }
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return csv_path, json_path


def relation_label(value: Optional[bool]) -> str:
    if value is True:
        return "true"
    if value is False:
        return "false"
    return "undecided"


def calibrate_local_cone_seed(n: int, seed: int, mass: float, margin: float) -> tuple[dict[str, object], list[dict[str, object]]]:
    r_plus = kerr_horizon_radius(mass, 0.0)
    events = generate_exterior_events(n, seed, r_plus + margin)
    possible_pairs = n * (n - 1) // 2
    rows: list[dict[str, object]] = []
    summary = {
        "local_decided_pairs": 0,
        "local_true_agree": 0,
        "local_false_agree": 0,
        "local_true_disagree": 0,
        "local_false_disagree": 0,
        "local_false_positive_vs_schwarzschild": 0,
        "local_false_negative_vs_schwarzschild": 0,
        "local_undecided_pairs": 0,
    }

    for i in range(n - 1):
        for j in range(i + 1, n):
            p = events[i]
            q = events[j]
            schwarzschild_relation = causal_relation_schwarzschild(p, q, mass)
            local_relation = local_cone_relation_kerr(p, q, mass, 0.0)
            if local_relation is None:
                summary["local_undecided_pairs"] += 1
                continue

            ds2, _spatial_distance = local_cone_ds2_midpoint(p, q, mass, 0.0)
            summary["local_decided_pairs"] += 1
            if local_relation is True and schwarzschild_relation is True:
                summary["local_true_agree"] += 1
            elif local_relation is False and schwarzschild_relation is False:
                summary["local_false_agree"] += 1
            elif local_relation is True and schwarzschild_relation is False:
                summary["local_true_disagree"] += 1
                summary["local_false_positive_vs_schwarzschild"] += 1
            elif local_relation is False and schwarzschild_relation is True:
                summary["local_false_disagree"] += 1
                summary["local_false_negative_vs_schwarzschild"] += 1

            if schwarzschild_relation is None:
                agree: object = "reference_undecided"
            else:
                agree = local_relation is schwarzschild_relation
            rows.append(
                {
                    "seed": seed,
                    "i": i,
                    "j": j,
                    "schwarzschild_relation": relation_label(schwarzschild_relation),
                    "local_relation": relation_label(local_relation),
                    "agree": agree,
                    "delta_t": q.t - p.t,
                    "r_i": p.r,
                    "r_j": q.r,
                    "theta_i": p.theta,
                    "theta_j": q.theta,
                    "phi_i": p.phi,
                    "phi_j": q.phi,
                    "ds2_midpoint": ds2,
                }
            )

    if summary["local_decided_pairs"] + summary["local_undecided_pairs"] != possible_pairs:
        raise AssertionError("calibration pair accounting mismatch")
    return summary, rows


def write_calibration_outputs(
    summary: dict[str, object],
    rows: list[dict[str, object]],
    out_prefix: str,
) -> tuple[Path, Path]:
    csv_path = OUT_DIR / f"{out_prefix}.csv"
    json_path = OUT_DIR / f"{out_prefix}.json"

    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(CALIBRATION_HEADERS)
        for row in rows:
            writer.writerow([_fmt(row[header]) for header in CALIBRATION_HEADERS])

    payload = {
        "summary": summary,
        "local_decided_pair_rows": rows,
        "notes": [
            "Calibration compares local BL midpoint-cone decisions at a=0 against the existing Schwarzschild benchmark subset.",
            "Schwarzschild undecided reference pairs are reported but are not counted as calibration contradictions.",
        ],
    }
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return csv_path, json_path


def run_local_cone_a0_calibration(args: argparse.Namespace) -> tuple[dict[str, object], list[dict[str, object]]]:
    if abs(args.a) > 0.0:
        raise SystemExit("--calibrate-local-cone-a0 requires --a 0.0")

    aggregate = {
        "benchmark": "S4-K2b local cone a=0 calibration",
        "status": "local_cone_a0_calibrated_against_schwarzschild_subset",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "command": COMMAND,
        "N": args.N,
        "seed": args.seed,
        "num_seeds": args.num_seeds,
        "M": args.M,
        "a": args.a,
        "r_min_margin": args.r_min_margin,
        "local_cone_margin_used": LOCAL_CONE_MARGIN,
        "local_decided_pairs": 0,
        "local_true_agree": 0,
        "local_false_agree": 0,
        "local_true_disagree": 0,
        "local_false_disagree": 0,
        "local_false_positive_vs_schwarzschild": 0,
        "local_false_negative_vs_schwarzschild": 0,
        "local_undecided_pairs": 0,
        "calibration_pass": True,
        "warning": "",
    }
    all_rows: list[dict[str, object]] = []
    count_keys = (
        "local_decided_pairs",
        "local_true_agree",
        "local_false_agree",
        "local_true_disagree",
        "local_false_disagree",
        "local_false_positive_vs_schwarzschild",
        "local_false_negative_vs_schwarzschild",
        "local_undecided_pairs",
    )

    for seed in range(args.seed, args.seed + args.num_seeds):
        summary, rows = calibrate_local_cone_seed(args.N, seed, args.M, args.r_min_margin)
        for key in count_keys:
            aggregate[key] += int(summary[key])
        all_rows.extend(rows)

    aggregate["calibration_pass"] = (
        aggregate["local_false_positive_vs_schwarzschild"] == 0
        and aggregate["local_false_negative_vs_schwarzschild"] == 0
    )
    if not aggregate["calibration_pass"]:
        aggregate["warning"] = "local cone diagnostic disagrees with decided Schwarzschild benchmark pairs"
    return aggregate, all_rows


def print_calibration_summary(summary: dict[str, object], csv_path: Path, json_path: Path) -> None:
    print("S4-K2b local cone a=0 calibration")
    print(f"N={summary['N']} seed={summary['seed']} num_seeds={summary['num_seeds']}")
    print(f"M={summary['M']} a={summary['a']} local_cone_margin_used={summary['local_cone_margin_used']:.6g}")
    print(
        f"local_decided_pairs={summary['local_decided_pairs']} "
        f"local_undecided_pairs={summary['local_undecided_pairs']}"
    )
    print(
        f"local_true_agree={summary['local_true_agree']} "
        f"local_false_agree={summary['local_false_agree']} "
        f"local_true_disagree={summary['local_true_disagree']} "
        f"local_false_disagree={summary['local_false_disagree']}"
    )
    print(
        "local_false_positive_vs_schwarzschild="
        f"{summary['local_false_positive_vs_schwarzschild']} "
        "local_false_negative_vs_schwarzschild="
        f"{summary['local_false_negative_vs_schwarzschild']}"
    )
    print(f"calibration_pass={summary['calibration_pass']}")
    if summary["warning"]:
        print(f"warning={summary['warning']}")
    print(f"wrote {csv_path}")
    print(f"wrote {json_path}")


def print_summary(summary: dict[str, object], csv_path: Path, json_path: Path) -> None:
    print("S4-K3 Kerr minimal benchmark")
    print("Regime: exterior only, r > r_plus + margin")
    print(f"N={summary['N']} seed={summary['seed']}")
    print(
        f"M={summary['M']} a={summary['a']} "
        f"r_plus={summary['r_plus']:.12g} margin={summary['r_min_margin']}"
    )
    print(f"kerr_mode={summary['kerr_mode']}")
    print(
        f"true_relations={summary['true_relations']} "
        f"false_relations={summary['false_relations']} "
        f"undecided_pairs={summary['undecided_pairs']}"
    )
    print(
        f"local_true_relations={summary['local_true_relations']} "
        f"local_false_relations={summary['local_false_relations']} "
        f"local_cone_margin_used={summary['local_cone_margin_used']:.6g}"
    )
    print(
        f"local_timelike_candidates={summary['local_timelike_candidates']} "
        f"local_spacelike_candidates={summary['local_spacelike_candidates']} "
        f"local_diagnostic_only={summary['local_diagnostic_only']}"
    )
    print(
        f"prograde_pairs={summary['prograde_pairs']} "
        f"retrograde_pairs={summary['retrograde_pairs']} "
        f"mean_abs_delta_phi={summary['mean_abs_delta_phi']:.6g} "
        f"mean_delta_t={summary['mean_delta_t']:.6g}"
    )
    print(
        f"r_min_observed={summary['r_min_observed']:.12g} "
        f"r_max_observed={summary['r_max_observed']:.12g}"
    )
    print(f"decided_pairs={summary['decided_pairs']}")
    print(f"ordering_fraction_decided={summary['ordering_fraction_decided']:.6g}")
    print(
        f"antisymmetric={summary['antisymmetric']} "
        f"transitive_true_matrix={summary['transitive_true_matrix']} "
        f"failed_checks={_fmt(summary['failed_checks'])}"
    )
    if summary["warning"]:
        print(f"warning={summary['warning']}")
    print(f"causal_relation_model={summary['causal_relation_model']}")
    print(f"wrote {csv_path}")
    print(f"wrote {json_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the minimal Kerr exterior scaffold benchmark.")
    parser.add_argument("--N", type=int, default=DEFAULT_N, help="number of exterior events")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help="random seed")
    parser.add_argument("--num-seeds", type=int, default=DEFAULT_NUM_SEEDS, help="number of consecutive seeds")
    parser.add_argument("--M", type=float, default=DEFAULT_M, help="Kerr mass parameter")
    parser.add_argument("--a", type=float, default=DEFAULT_A, help="Kerr spin parameter")
    parser.add_argument(
        "--r-min-margin",
        type=float,
        default=DEFAULT_R_MIN_MARGIN,
        help="minimum radius margin above r_plus",
    )
    parser.add_argument(
        "--out-prefix",
        default=DEFAULT_OUT_PREFIX,
        help="output filename prefix under explore/sorkin4_kerr_benchmark/",
    )
    parser.add_argument(
        "--kerr-mode",
        choices=KERR_MODES,
        default=DEFAULT_KERR_MODE,
        help="Kerr a!=0 relation mode",
    )
    parser.add_argument(
        "--calibrate-local-cone-a0",
        action="store_true",
        help="compare local_cone_diagnostic decisions at a=0 against the Schwarzschild benchmark subset",
    )
    args = parser.parse_args()
    if args.N <= 0:
        raise SystemExit("--N must be positive")
    if args.num_seeds <= 0:
        raise SystemExit("--num-seeds must be positive")
    if args.r_min_margin <= 0.0:
        raise SystemExit("--r-min-margin must be positive")
    return args


def main() -> None:
    args = parse_args()
    if args.calibrate_local_cone_a0:
        summary, rows = run_local_cone_a0_calibration(args)
        csv_path, json_path = write_calibration_outputs(summary, rows, args.out_prefix)
        print_calibration_summary(summary, csv_path, json_path)
        return
    if args.num_seeds != 1:
        raise SystemExit("--num-seeds is only supported with --calibrate-local-cone-a0")

    r_plus = kerr_horizon_radius(args.M, args.a)
    r_min = r_plus + args.r_min_margin
    events = generate_exterior_events(
        args.N,
        args.seed,
        r_min,
        equatorial=args.kerr_mode == "equatorial_scaffold",
    )
    matrix, states = build_relation_states(events, args.M, args.a, args.kerr_mode)
    summary = summarize_case(
        events,
        matrix,
        states,
        n=args.N,
        seed=args.seed,
        mass=args.M,
        spin=args.a,
        r_plus=r_plus,
        r_min=r_min,
        margin=args.r_min_margin,
        kerr_mode=args.kerr_mode,
    )
    csv_path, json_path = write_outputs(events, matrix, states, summary, args.out_prefix)
    print_summary(summary, csv_path, json_path)


if __name__ == "__main__":
    main()
