#!/usr/bin/env python3
"""Validate causal_relation_schwarzschild against Table 2 of He & Rideout (2009).

Reference: He & Rideout, "A Causal Set Black Hole" (2009).

Table 1 events and Table 2 intermediate values are transcribed from the paper.
This script checks selected exterior pairs against the paper's known decisions
and intermediate values before accepting the implementation as ready for Kerr.

Selected pairs tested:
  0-1  outgoing, unrelated by radial spacelike bound
  0-6  ingoing,  unrelated by angular spacelike bound
  1-6  ingoing,  generic exterior, related (requires shooting)

Intermediate values checked:
  phi2       angular separation (Eq. 10 angular rotation to equatorial plane)
  dt         EF time separation
  rad_bound  radial spacelike lower bound (Eq. 6/7)
  ang_bound  angular spacelike lower bound (Eq. 9)
  tot_trip   composed radial+angular timelike upper bound (Sec. 2.3.2)
  c2         shooting parameter (E/L)^2 for pair 1-6
  null_dt    elapsed EF time along fastest null geodesic for pair 1-6

Tolerances are generous (1e-3) because Table 2 values are printed with
limited decimal places; tighter tolerances test internal consistency.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

# Make the runner importable from this directory.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

# Import internal helpers from the benchmark runner.
sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_schwarzschild_minimal_benchmark import (  # noqa: E402
    Event,
    MASS,
    SCHWARZSCHILD_RADIUS,
    TIME_EPS,
    angular_separation,
    _outgoing_radial_trip,
    _angular_spacelike_r0,
    _angular_f,
    _timelike_bound_r0,
    find_direct_shooting_c2,
    direct_time_integral,
    causal_relation_schwarzschild,
)

# ---------------------------------------------------------------------------
# Table 1 events (exact paper coordinates, Eddington-Finkelstein)
# He & Rideout Table 1: 9 events, M=1 (horizon at r=2)
# ---------------------------------------------------------------------------
TABLE1: dict[int, Event] = {
    0: Event(index=0, t=0.410895, r=2.36161,  theta=1.80295, phi=0.57951),
    1: Event(index=1, t=1.109415, r=2.89891,  theta=1.04335, phi=4.25531),
    2: Event(index=2, t=1.133105, r=1.36083,  theta=1.89919, phi=1.06482),  # interior
    3: Event(index=3, t=2.743428, r=2.74093,  theta=2.97906, phi=4.22204),
    4: Event(index=4, t=3.235970, r=0.65462,  theta=0.11664, phi=5.06884),  # interior
    5: Event(index=5, t=3.972871, r=0.96354,  theta=2.33727, phi=1.38169),  # interior
    6: Event(index=6, t=5.230757, r=2.34476,  theta=1.11855, phi=3.47242),
    7: Event(index=7, t=6.014261, r=0.664739, theta=2.82235, phi=0.95459),  # interior
    8: Event(index=8, t=6.193089, r=0.429636, theta=2.20122, phi=1.99644),  # interior
}

# ---------------------------------------------------------------------------
# Expected intermediate values from Table 2
# All quantities are from the paper (M=1, horizon=2).
# ---------------------------------------------------------------------------

# Pair 0-1: outgoing (r1=2.36161, r2=2.89891), unrelated by radial bound.
#   dir=out, r0=2.898906, phi2=2.567258, dt=0.698520,
#   rad_trip=4.179694, ang_bnd=13.36484, tot_trip=—
#   result: unrelated: either bound
PAIR_01_PHI2        = 2.567258
PAIR_01_DT          = 0.698520
PAIR_01_RAD_BOUND   = 4.179694   # (r2-r1) + 4M*ln((r2-2M)/(r1-2M))
PAIR_01_RESULT      = False

# Pair 0-6: ingoing (r1=2.36161, r2=2.34476), unrelated by angular bound.
#   dir=in, r0=2.361614, phi2=2.820685, dt=4.819862,
#   rad_trip=0.016859, ang_bnd=6.613817, tot_trip=—
#   result: unrelated: angular bound
PAIR_06_PHI2        = 2.820685
PAIR_06_DT          = 4.819862
PAIR_06_RAD_BOUND   = 0.016859   # r1 - r2 (ingoing, simple)
PAIR_06_ANG_BOUND   = 6.613817   # r2 * phi2
PAIR_06_RESULT      = False

# Pair 1-6: ingoing (r1=2.89891, r2=2.34476), generic exterior, related.
#   dir=in, r0=2.898906, phi2=0.690536, dt=4.121342,
#   rad_trip=0.554150, ang_bnd=1.619139, tot_trip=4.148999,
#   c2=0.0476468, null_dt=2.60973
#   result: generic, exterior, related
PAIR_16_PHI2        = 0.690536
PAIR_16_DT          = 4.121342
PAIR_16_RAD_BOUND   = 0.554150   # r1 - r2 (ingoing, simple)
PAIR_16_ANG_BOUND   = 1.619139   # r2 * phi2
PAIR_16_TOT_TRIP    = 4.148999   # composed bound
PAIR_16_C2          = 0.0476468
PAIR_16_NULL_DT     = 2.60973
PAIR_16_RESULT      = True

# ---------------------------------------------------------------------------
# Tolerance
# ---------------------------------------------------------------------------
TOL_ANGLE    = 2e-4   # angular values printed to 6 decimal places in table
TOL_TRIP     = 2e-4   # bound values printed to 6 decimal places
TOL_C2       = 1e-3   # c2 printed to 6-7 significant figures
TOL_NULL_DT  = 5e-3   # null_dt printed to 5 significant figures


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _close(actual: float, expected: float, tol: float, label: str) -> bool:
    ok = abs(actual - expected) <= tol
    status = "PASS" if ok else "FAIL"
    marker = "" if ok else f"  <-- expected {expected:.8g}, got {actual:.8g}, diff {actual-expected:+.4g}"
    print(f"  [{status}] {label}: {actual:.8g}{marker}")
    return ok


def _result_ok(actual, expected, label: str) -> bool:
    ok = actual == expected
    status = "PASS" if ok else "FAIL"
    marker = f"  <-- expected {expected!r}, got {actual!r}" if not ok else ""
    print(f"  [{status}] {label}: {actual!r}{marker}")
    return ok


# ---------------------------------------------------------------------------
# Test implementations
# ---------------------------------------------------------------------------

def test_pair_01() -> bool:
    """Pair 0-1: outgoing exterior, unrelated by radial spacelike bound."""
    print("=== Pair 0-1: outgoing, unrelated by radial bound ===")
    p = TABLE1[0]
    q = TABLE1[1]

    # p.t < q.t → E1=p, E2=q; r1=p.r=2.36161, r2=q.r=2.89891 (outgoing)
    r1, r2 = p.r, q.r
    assert r2 >= r1, f"expected outgoing: r2={r2} < r1={r1}"

    phi2 = angular_separation(p, q)
    dt = q.t - p.t

    # Outgoing radial bound: (r2-r1) + 4M*ln((r2-2M)/(r1-2M))
    # Our code: _outgoing_radial_trip(r1, r2) = (r2-r1) + 2M*ln(...)  <- CHECK
    rad_bound_code = _outgoing_radial_trip(r1, r2)
    rad_bound_paper = (r2 - r1) + 4.0 * MASS * math.log(
        (r2 - SCHWARZSCHILD_RADIUS) / (r1 - SCHWARZSCHILD_RADIUS)
    )

    passes: list[bool] = []
    passes.append(_close(phi2, PAIR_01_PHI2, TOL_ANGLE, "phi2"))
    passes.append(_close(dt, PAIR_01_DT, TOL_TRIP, "dt"))
    passes.append(_close(rad_bound_paper, PAIR_01_RAD_BOUND, TOL_TRIP,
                         "rad_bound 4M*ln (paper formula)"))
    passes.append(_close(rad_bound_code, PAIR_01_RAD_BOUND, TOL_TRIP,
                         "rad_bound _outgoing_radial_trip (our code)"))

    result = causal_relation_schwarzschild(p, q, enable_shooting=False)
    passes.append(_result_ok(result, PAIR_01_RESULT, "causal_relation"))

    ok = all(passes)
    print(f"  --> pair 0-1: {'PASS' if ok else 'FAIL'}\n")
    return ok


def test_pair_06() -> bool:
    """Pair 0-6: ingoing exterior, unrelated by angular spacelike bound."""
    print("=== Pair 0-6: ingoing, unrelated by angular bound ===")
    p = TABLE1[0]
    q = TABLE1[6]

    # p.t=0.410895 < q.t=5.230757 → E1=p, E2=q; r1=2.36161, r2=2.34476 (ingoing)
    r1, r2 = p.r, q.r
    assert r1 >= r2, f"expected ingoing: r1={r1} < r2={r2}"

    phi2 = angular_separation(p, q)
    dt = q.t - p.t
    rad_bound = r1 - r2            # ingoing simple bound
    ang_bound = r2 * phi2          # He & Rideout ingoing angular spacelike bound

    passes: list[bool] = []
    passes.append(_close(phi2, PAIR_06_PHI2, TOL_ANGLE, "phi2"))
    passes.append(_close(dt, PAIR_06_DT, TOL_TRIP, "dt"))
    passes.append(_close(rad_bound, PAIR_06_RAD_BOUND, TOL_TRIP, "rad_bound r1-r2"))
    passes.append(_close(ang_bound, PAIR_06_ANG_BOUND, TOL_TRIP, "ang_bound r2*phi2"))

    result = causal_relation_schwarzschild(p, q, enable_shooting=False)
    passes.append(_result_ok(result, PAIR_06_RESULT, "causal_relation"))

    ok = all(passes)
    print(f"  --> pair 0-6: {'PASS' if ok else 'FAIL'}\n")
    return ok


def test_pair_16() -> bool:
    """Pair 1-6: ingoing exterior, generic, related (requires shooting)."""
    print("=== Pair 1-6: ingoing, generic exterior, related (requires shooting) ===")
    p = TABLE1[1]
    q = TABLE1[6]

    # p.t=1.109415 < q.t=5.230757; r1=2.89891, r2=2.34476 (ingoing)
    r1, r2 = p.r, q.r
    assert r1 >= r2, f"expected ingoing: r1={r1} < r2={r2}"

    phi2 = angular_separation(p, q)
    dt = q.t - p.t
    rad_bound = r1 - r2
    ang_bound = r2 * phi2

    # Composed timelike bound: radial_trip + angular_trip
    r0_timelike = _timelike_bound_r0(r1, r2)   # should be r1 for both below photon orbit
    angular_trip = _angular_f(r0_timelike) * phi2
    tot_trip = rad_bound + angular_trip

    passes: list[bool] = []
    passes.append(_close(phi2, PAIR_16_PHI2, TOL_ANGLE, "phi2"))
    passes.append(_close(dt, PAIR_16_DT, TOL_TRIP, "dt"))
    passes.append(_close(rad_bound, PAIR_16_RAD_BOUND, TOL_TRIP, "rad_bound r1-r2"))
    passes.append(_close(ang_bound, PAIR_16_ANG_BOUND, TOL_TRIP, "ang_bound r2*phi2"))
    passes.append(_close(tot_trip, PAIR_16_TOT_TRIP, TOL_TRIP, "tot_trip (composed bound)"))

    # Without shooting: pair should be undecided (tot_trip > dt, so timelike bound doesn't fire)
    result_no_shoot = causal_relation_schwarzschild(p, q, enable_shooting=False)
    passes.append(_result_ok(result_no_shoot, None, "causal_relation (no shooting)"))

    # Shooting: find c2
    u1, u2 = 1.0 / r1, 1.0 / r2
    c2_result = find_direct_shooting_c2(u1, u2, phi2)
    if c2_result is None:
        passes.append(_result_ok(None, "c2_found", "find_direct_shooting_c2"))
    else:
        c2, phi_achieved = c2_result
        passes.append(_close(c2, PAIR_16_C2, TOL_C2, "c2 (E/L)^2"))
        null_dt = direct_time_integral(u1, u2, c2)
        if null_dt is None:
            passes.append(_result_ok(None, "null_dt_found", "direct_time_integral"))
        else:
            passes.append(_close(null_dt, PAIR_16_NULL_DT, TOL_NULL_DT, "null_dt (elapsed EF time)"))
            # Causal criterion: null_dt <= dt → related
            criterion = null_dt <= dt + TIME_EPS
            passes.append(_result_ok(criterion, True, "null_dt <= dt (causal criterion)"))

    # With shooting: should return True
    result_shoot = causal_relation_schwarzschild(p, q, enable_shooting=True)
    passes.append(_result_ok(result_shoot, PAIR_16_RESULT, "causal_relation (with shooting)"))

    ok = all(passes)
    print(f"  --> pair 1-6: {'PASS' if ok else 'FAIL'}\n")
    return ok


def test_outgoing_radial_formula() -> bool:
    """Explicit check of the outgoing radial trip formula against He & Rideout Eq. 7.

    The paper's recipe (Section 2.2) states:
      for r2 >= r1 > 2M, related iff t2 >= t1 + r2 - r1 + 4M*ln((r2-2M)/(r1-2M))

    So _outgoing_radial_trip(r1, r2) must return (r2-r1) + 4M*ln((r2-2M)/(r1-2M)).
    """
    print("=== _outgoing_radial_trip formula check ===")
    r1 = TABLE1[0].r   # 2.36161
    r2 = TABLE1[1].r   # 2.89891
    expected = PAIR_01_RAD_BOUND

    code_result = _outgoing_radial_trip(r1, r2)
    paper_result = (r2 - r1) + 4.0 * MASS * math.log(
        (r2 - SCHWARZSCHILD_RADIUS) / (r1 - SCHWARZSCHILD_RADIUS)
    )

    passes: list[bool] = []
    passes.append(_close(paper_result, expected, TOL_TRIP, "paper 4M formula matches Table 2"))
    passes.append(_close(code_result, expected, TOL_TRIP, "_outgoing_radial_trip matches Table 2"))

    ok = all(passes)
    print(f"  --> formula check: {'PASS' if ok else 'FAIL'}\n")
    return ok


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print("He & Rideout (2009) Table 2 validation")
    print(f"M={MASS} horizon={SCHWARZSCHILD_RADIUS}")
    print()

    results = [
        test_outgoing_radial_formula(),
        test_pair_01(),
        test_pair_06(),
        test_pair_16(),
    ]

    passed = sum(results)
    total = len(results)
    print(f"=== Summary: {passed}/{total} test groups passed ===")
    if passed < total:
        print("VALIDATION FAILED — do not proceed to Kerr")
        sys.exit(1)
    else:
        print("VALIDATION PASSED — Schwarzschild implementation verified against Table 2")
        sys.exit(0)


if __name__ == "__main__":
    main()
