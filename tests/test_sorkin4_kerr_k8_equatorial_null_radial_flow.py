"""Regression tests for S4-KERR-K8-EQUATORIAL-NULL-RADIAL-FLOW-001.

Three test classes:
  KerrK8ArtifactFileTests  — artifact files exist and have the right shape.
  KerrK8ArtifactRowTests   — CSV row invariants and all_checks_pass flags.
  KerrK8ComputationTests   — unit tests on physics and integrator functions.
"""

from __future__ import annotations

import csv
import json
import math
import sys
import unittest
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

KERR_DIR = ROOT / "explore" / "sorkin4_kerr_benchmark"
sys.path.insert(0, str(KERR_DIR))

from audit_kerr_k8_equatorial_null_radial_flow_001 import (
    CASE_CIRCULAR,
    CASE_INGOING,
    CASE_OUTGOING,
    CIRCULAR_DRIFT_TOL,
    D_LAMBDA,
    DEFAULT_MASS,
    DEFAULT_SEED,
    DEFAULT_SPINS,
    LAMBDA_FINAL,
    N_STEPS,
    OUT_PREFIX,
    RHS_CONSISTENCY_TOL,
    SCHW_LIMIT_TOL,
    _rhs_consistency_max_error,
    integrate_trajectory,
    null_flow_rhs,
    null_radial_potential,
    photon_impact_parameter,
    photon_sphere_radius_pro,
    rk4_step,
    run_audit,
    run_trajectory_checks,
)

ARTIFACT_DIR = KERR_DIR
CSV_PATH  = ARTIFACT_DIR / f"{OUT_PREFIX}.csv"
JSON_PATH = ARTIFACT_DIR / f"{OUT_PREFIX}.json"
MD_PATH   = ARTIFACT_DIR / f"{OUT_PREFIX}.md"
PNG_PATH  = ARTIFACT_DIR / f"{OUT_PREFIX}.png"

MASS = DEFAULT_MASS


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_csv() -> list[dict[str, str]]:
    with CSV_PATH.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _load_json() -> dict:
    return json.loads(JSON_PATH.read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------
# KerrK8ArtifactFileTests
# ---------------------------------------------------------------------------

class KerrK8ArtifactFileTests(unittest.TestCase):

    def test_csv_exists(self):
        assert CSV_PATH.exists(), f"Missing {CSV_PATH}"

    def test_json_exists(self):
        assert JSON_PATH.exists(), f"Missing {JSON_PATH}"

    def test_md_exists(self):
        assert MD_PATH.exists(), f"Missing {MD_PATH}"

    def test_png_exists(self):
        assert PNG_PATH.exists(), f"Missing {PNG_PATH}"

    def test_csv_row_count(self):
        rows = _load_csv()
        expected = len(DEFAULT_SPINS) * 3  # 3 cases per spin
        assert len(rows) == expected, f"Expected {expected} rows, got {len(rows)}"

    def test_json_audit_id(self):
        payload = _load_json()
        assert payload["aggregate"]["audit"] == "S4-KERR-K8-EQUATORIAL-NULL-RADIAL-FLOW-001"

    def test_json_all_checks_pass(self):
        payload = _load_json()
        assert payload["aggregate"]["all_checks_pass"] is True

    def test_json_spins_match_defaults(self):
        payload = _load_json()
        assert payload["aggregate"]["spins"] == list(DEFAULT_SPINS)

    def test_json_n_steps(self):
        payload = _load_json()
        assert payload["aggregate"]["n_steps"] == N_STEPS

    def test_json_d_lambda(self):
        payload = _load_json()
        assert abs(payload["aggregate"]["d_lambda"] - D_LAMBDA) < 1e-15

    def test_json_positive_spin_all_undecided(self):
        payload = _load_json()
        assert payload["aggregate"]["positive_spin_cases_all_undecided"] is True

    def test_md_contains_audit_id(self):
        text = MD_PATH.read_text(encoding="utf-8")
        assert "S4-KERR-K8-EQUATORIAL-NULL-RADIAL-FLOW-001" in text

    def test_md_not_causal_solver(self):
        text = MD_PATH.read_text(encoding="utf-8").lower()
        assert "not a kerr causal solver" in text or "not a global kerr causal" in text

    def test_png_nonzero_size(self):
        assert PNG_PATH.stat().st_size > 1000


# ---------------------------------------------------------------------------
# KerrK8ArtifactRowTests
# ---------------------------------------------------------------------------

class KerrK8ArtifactRowTests(unittest.TestCase):

    def test_all_rows_pass(self):
        for row in _load_csv():
            spin   = float(row["spin_a"])
            case   = row["case_id"]
            passed = row["all_checks_pass"].strip().lower() == "true"
            assert passed, f"all_checks_pass=False for spin={spin}, case={case}"

    def test_case_ids_present(self):
        rows     = _load_csv()
        case_ids = {row["case_id"] for row in rows}
        assert CASE_OUTGOING in case_ids
        assert CASE_INGOING  in case_ids
        assert CASE_CIRCULAR in case_ids

    def test_each_spin_has_three_cases(self):
        rows = _load_csv()
        from collections import Counter
        counts = Counter(row["spin_a"] for row in rows)
        for spin, count in counts.items():
            assert count == 3, f"spin={spin} has {count} rows, expected 3"

    def test_all_points_exterior(self):
        for row in _load_csv():
            assert row["all_points_exterior"].strip().lower() == "true", (
                f"all_points_exterior=False for spin={row['spin_a']}, case={row['case_id']}"
            )

    def test_min_R_nonneg(self):
        for row in _load_csv():
            assert row["min_R_nonneg"].strip().lower() == "true", (
                f"min_R_nonneg=False for spin={row['spin_a']}, case={row['case_id']}"
            )

    def test_rhs_consistency_pass(self):
        for row in _load_csv():
            assert row["rhs_consistency_pass"].strip().lower() == "true", (
                f"rhs_consistency_pass=False for spin={row['spin_a']}, case={row['case_id']}"
            )

    def test_rhs_max_error_below_tol(self):
        for row in _load_csv():
            err = float(row["rhs_max_error"])
            assert err <= RHS_CONSISTENCY_TOL, (
                f"rhs_max_error={err:.2e} > {RHS_CONSISTENCY_TOL} "
                f"for spin={row['spin_a']}, case={row['case_id']}"
            )

    def test_monotonic_pass(self):
        for row in _load_csv():
            assert row["monotonic_pass"].strip().lower() == "true", (
                f"monotonic_pass=False for spin={row['spin_a']}, case={row['case_id']}"
            )

    def test_circular_drift_pass(self):
        for row in _load_csv():
            assert row["circular_drift_pass"].strip().lower() == "true", (
                f"circular_drift_pass=False for spin={row['spin_a']}, case={row['case_id']}"
            )

    def test_circular_drift_present_only_for_circular(self):
        for row in _load_csv():
            if row["case_id"] == CASE_CIRCULAR:
                assert row["circular_drift"].strip() != "", (
                    f"circular_drift empty for circular case, spin={row['spin_a']}"
                )
            else:
                assert row["circular_drift"].strip() == "", (
                    f"circular_drift non-empty for case={row['case_id']}"
                )

    def test_schwarzschild_limit_pass(self):
        for row in _load_csv():
            assert row["schwarzschild_radial_limit_pass"].strip().lower() == "true", (
                f"schwarzschild_radial_limit_pass=False for "
                f"spin={row['spin_a']}, case={row['case_id']}"
            )

    def test_schwarzschild_limit_error_near_zero_at_a0_b0(self):
        rows = _load_csv()
        for row in rows:
            spin = float(row["spin_a"])
            case = row["case_id"]
            if abs(spin) <= 0.0 and case in (CASE_OUTGOING, CASE_INGOING):
                err_str = row["schwarzschild_limit_error"].strip()
                assert err_str != "", f"schwarzschild_limit_error empty for a=0, case={case}"
                err = float(err_str)
                assert err <= SCHW_LIMIT_TOL, (
                    f"Schwarzschild limit error {err:.2e} > {SCHW_LIMIT_TOL} "
                    f"for a=0, case={case}"
                )

    def test_outgoing_r_final_greater_than_r0(self):
        for row in _load_csv():
            if row["case_id"] == CASE_OUTGOING:
                assert float(row["r_final"]) > float(row["r0"]), (
                    f"outgoing r_final <= r0 for spin={row['spin_a']}"
                )

    def test_ingoing_r_final_less_than_r0(self):
        for row in _load_csv():
            if row["case_id"] == CASE_INGOING:
                assert float(row["r_final"]) < float(row["r0"]), (
                    f"ingoing r_final >= r0 for spin={row['spin_a']}"
                )

    def test_causal_accounting_positive_spins(self):
        payload       = _load_json()
        possible_pairs = payload["aggregate"]["possible_pairs"]
        for row in payload["rows"]:
            if row["spin_a"] > 0.0:
                assert row["global_true_relations"]  == 0,  f"spin={row['spin_a']}"
                assert row["global_false_relations"] == 0,  f"spin={row['spin_a']}"
                assert row["global_undecided_pairs"] == possible_pairs, (
                    f"spin={row['spin_a']}: undecided={row['global_undecided_pairs']} "
                    f"!= {possible_pairs}"
                )

    def test_no_global_causal_relations_decided(self):
        for row in _load_json()["rows"]:
            if row["spin_a"] > 0.0:
                assert row["no_global_causal_relations_decided"] is True, (
                    f"no_global_causal_relations_decided=False for spin={row['spin_a']}"
                )


# ---------------------------------------------------------------------------
# KerrK8ComputationTests
# ---------------------------------------------------------------------------

class KerrK8ComputationTests(unittest.TestCase):

    # --- null_radial_potential ---

    def test_null_potential_b0_always_positive(self):
        """R(r; a, b=0) = r^2*(r^2+a^2) + 2*M*a^2*r > 0 for r > 0."""
        mass = 1.0
        for spin in (0.0, 0.25, 0.5, 0.75):
            for r in (1.5, 2.0, 3.0, 5.0, 10.0):
                R = null_radial_potential(r, spin, 0.0, mass)
                assert R > 0.0, f"R <= 0 for a={spin}, r={r}: R={R}"

    def test_null_potential_schw_b0(self):
        """At a=0, b=0: R = r^4 exactly."""
        mass = 1.0
        spin = 0.0
        for r in (2.1, 3.0, 5.0, 10.0):
            R = null_radial_potential(r, spin, 0.0, mass)
            assert abs(R - r ** 4) < 1e-12 * r ** 4, f"R={R}, expected r^4={r**4}"

    def test_null_potential_circular_orbit_a0(self):
        """At a=0, r=3M, b=3sqrt(3)M: R=0 (Schwarzschild photon sphere)."""
        mass = 1.0
        spin = 0.0
        r_ph = 3.0 * mass
        b_ph = 3.0 * math.sqrt(3.0) * mass
        R = null_radial_potential(r_ph, spin, b_ph, mass)
        assert abs(R) < 1e-12, f"|R|={abs(R):.2e} at Schwarzschild photon sphere"

    def test_null_potential_circular_orbit_nonzero_spin(self):
        """R(r_ph_pro; a, b_ph_pro) ≈ 0 for a > 0."""
        mass = 1.0
        for spin in (0.25, 0.5, 0.75):
            r_ph = photon_sphere_radius_pro(mass, spin)
            b_ph = photon_impact_parameter(r_ph, mass, spin, prograde=True)
            R = null_radial_potential(r_ph, spin, b_ph, mass)
            assert abs(R) < 1e-9, (
                f"|R|={abs(R):.2e} at prograde circular orbit for a={spin}"
            )

    # --- null_flow_rhs ---

    def test_null_flow_rhs_schw_b0_constant(self):
        """At a=0, b=0: dr/dlambda = s * r^2/r^2 = s (constant, all r)."""
        mass = 1.0
        spin = 0.0
        for r in (2.1, 3.0, 5.0, 10.0, 20.0):
            for s in (+1.0, -1.0):
                rhs = null_flow_rhs(r, spin, 0.0, mass, s)
                assert abs(rhs - s) < 1e-14, (
                    f"RHS={rhs} != {s} at a=0, b=0, r={r}"
                )

    def test_null_flow_rhs_outgoing_positive_b0(self):
        """b=0 outgoing RHS > 0 for all r > 0."""
        mass = 1.0
        for spin in (0.0, 0.25, 0.5, 0.75):
            for r in (2.1, 3.0, 5.0, 10.0):
                rhs = null_flow_rhs(r, spin, 0.0, mass, +1.0)
                assert rhs > 0.0, f"outgoing RHS <= 0 for a={spin}, r={r}"

    def test_null_flow_rhs_circular_near_zero(self):
        """At circular photon orbit, RHS ≈ 0 (R ≈ 0)."""
        mass = 1.0
        for spin in (0.0, 0.25, 0.5, 0.75):
            r_ph = photon_sphere_radius_pro(mass, spin)
            b_ph = photon_impact_parameter(r_ph, mass, spin, prograde=True)
            rhs  = null_flow_rhs(r_ph, spin, b_ph, mass, +1.0)
            assert abs(rhs) < 1e-4, f"|RHS|={abs(rhs):.2e} at circular orbit for a={spin}"

    # --- rk4_step ---

    def test_rk4_step_schw_b0_exact(self):
        """a=0, b=0: RHS=+1 (constant), so rk4_step = r + dlambda exactly."""
        mass = 1.0
        r    = 5.0
        for dlambda in (0.05, 0.1, 0.2):
            r_new = rk4_step(r, 0.0, 0.0, mass, +1.0, dlambda)
            assert abs(r_new - (r + dlambda)) < 1e-14 * r, (
                f"rk4_step error={abs(r_new - (r+dlambda)):.2e} for dlambda={dlambda}"
            )

    def test_rk4_step_outgoing_increases_r(self):
        """b=0 outgoing: each RK4 step increases r."""
        mass = 1.0
        r    = 5.0
        for spin in (0.0, 0.25, 0.5, 0.75):
            r_new = rk4_step(r, spin, 0.0, mass, +1.0, D_LAMBDA)
            assert r_new > r, f"rk4 outgoing step did not increase r for a={spin}"

    def test_rk4_step_ingoing_decreases_r(self):
        """b=0 ingoing: each RK4 step decreases r."""
        mass = 1.0
        r    = 10.0
        for spin in (0.0, 0.25, 0.5, 0.75):
            r_new = rk4_step(r, spin, 0.0, mass, -1.0, D_LAMBDA)
            assert r_new < r, f"rk4 ingoing step did not decrease r for a={spin}"

    # --- integrate_trajectory ---

    def test_integrate_trajectory_length(self):
        """integrate_trajectory returns N_STEPS+1 points."""
        traj = integrate_trajectory(5.0, 0.0, 0.0, MASS, +1.0)
        assert len(traj) == N_STEPS + 1

    def test_integrate_trajectory_schw_b0_exact(self):
        """a=0, b=0, outgoing: r_final = r0 + LAMBDA_FINAL to machine precision."""
        r0   = 5.0
        traj = integrate_trajectory(r0, 0.0, 0.0, MASS, +1.0)
        err  = abs(traj[-1] - (r0 + LAMBDA_FINAL))
        assert err <= SCHW_LIMIT_TOL, (
            f"Schwarzschild outgoing limit error {err:.2e} > {SCHW_LIMIT_TOL}"
        )

    def test_integrate_trajectory_schw_b0_ingoing_exact(self):
        """a=0, b=0, ingoing: r_final = r0 - LAMBDA_FINAL to machine precision."""
        r0   = 10.0
        traj = integrate_trajectory(r0, 0.0, 0.0, MASS, -1.0)
        err  = abs(traj[-1] - (r0 - LAMBDA_FINAL))
        assert err <= SCHW_LIMIT_TOL, (
            f"Schwarzschild ingoing limit error {err:.2e} > {SCHW_LIMIT_TOL}"
        )

    def test_integrate_trajectory_b0_monotone_outgoing(self):
        """b=0 outgoing trajectory is non-decreasing for all spins."""
        for spin in DEFAULT_SPINS:
            traj = integrate_trajectory(5.0, spin, 0.0, MASS, +1.0)
            for i in range(len(traj) - 1):
                assert traj[i + 1] >= traj[i], (
                    f"Non-monotone outgoing at step {i} for a={spin}"
                )

    def test_integrate_trajectory_b0_monotone_ingoing(self):
        """b=0 ingoing trajectory is non-increasing for all spins."""
        for spin in DEFAULT_SPINS:
            traj = integrate_trajectory(10.0, spin, 0.0, MASS, -1.0)
            for i in range(len(traj) - 1):
                assert traj[i + 1] <= traj[i], (
                    f"Non-monotone ingoing at step {i} for a={spin}"
                )

    # --- RHS consistency ---

    def test_rhs_consistency_schw_b0_machine_precision(self):
        """a=0, b=0: RHS is constant; central-diff error ≈ 0."""
        traj = integrate_trajectory(5.0, 0.0, 0.0, MASS, +1.0)
        err  = _rhs_consistency_max_error(traj, 0.0, 0.0, MASS, +1.0, D_LAMBDA)
        assert err < 1e-12, f"RHS consistency error {err:.2e} for a=0 (expected ≈0)"

    def test_rhs_consistency_all_spins_b0_within_tol(self):
        """All spins, b=0 outgoing: RHS consistency ≤ RHS_CONSISTENCY_TOL."""
        for spin in DEFAULT_SPINS:
            traj = integrate_trajectory(5.0, spin, 0.0, MASS, +1.0)
            err  = _rhs_consistency_max_error(traj, spin, 0.0, MASS, +1.0, D_LAMBDA)
            assert err <= RHS_CONSISTENCY_TOL, (
                f"a={spin}: rhs_max_error={err:.2e} > {RHS_CONSISTENCY_TOL}"
            )

    # --- circular orbit drift ---

    def test_circular_orbit_drift_below_tol(self):
        """Prograde circular orbit drift < CIRCULAR_DRIFT_TOL for all spins."""
        for spin in DEFAULT_SPINS:
            r_ph = photon_sphere_radius_pro(MASS, spin)
            b_ph = photon_impact_parameter(r_ph, MASS, spin, prograde=True)
            traj = integrate_trajectory(r_ph, spin, b_ph, MASS, +1.0)
            drift = abs(traj[-1] - traj[0])
            assert drift < CIRCULAR_DRIFT_TOL, (
                f"a={spin}: circular drift={drift:.2e} >= {CIRCULAR_DRIFT_TOL}"
            )

    def test_circular_orbit_schw_drift_exactly_zero(self):
        """a=0, circular orbit (r=3M, b=3sqrt(3)M): RHS=0 everywhere => zero drift."""
        spin = 0.0
        r_ph = photon_sphere_radius_pro(MASS, spin)
        b_ph = photon_impact_parameter(r_ph, MASS, spin, prograde=True)
        traj = integrate_trajectory(r_ph, spin, b_ph, MASS, +1.0)
        assert traj[-1] == traj[0], f"a=0 circular drift non-zero: {traj[-1] - traj[0]}"

    # --- photon sphere helpers ---

    def test_photon_sphere_a0_is_3M(self):
        assert abs(photon_sphere_radius_pro(MASS, 0.0) - 3.0 * MASS) < 1e-12

    def test_photon_sphere_decreases_with_spin(self):
        spins = [0.0, 0.25, 0.5, 0.75]
        radii = [photon_sphere_radius_pro(MASS, a) for a in spins]
        for i in range(len(radii) - 1):
            assert radii[i + 1] < radii[i], (
                f"r_ph_pro did not decrease from a={spins[i]} to a={spins[i+1]}"
            )

    def test_photon_impact_parameter_a0(self):
        """At a=0: b_ph_pro = +3sqrt(3)M."""
        r_ph = photon_sphere_radius_pro(MASS, 0.0)
        b_ph = photon_impact_parameter(r_ph, MASS, 0.0, prograde=True)
        assert abs(b_ph - 3.0 * math.sqrt(3.0) * MASS) < 1e-12

    def test_photon_impact_raises_inside_horizon(self):
        """photon_impact_parameter raises if Delta <= 0."""
        with pytest.raises(ValueError):
            photon_impact_parameter(0.5, MASS, 0.0, prograde=True)

    # --- run_audit smoke test ---

    def test_run_audit_all_pass(self):
        """run_audit with default parameters returns all_checks_pass=True."""
        payload = run_audit()
        assert payload["aggregate"]["all_checks_pass"] is True

    def test_run_audit_row_count(self):
        payload = run_audit()
        assert len(payload["rows"]) == len(DEFAULT_SPINS) * 3

    def test_run_audit_invalid_mass_raises(self):
        with pytest.raises(ValueError, match="M=1"):
            run_audit(mass=2.0)

    def test_run_audit_invalid_spin_raises(self):
        with pytest.raises(ValueError, match=r"\|a\| < M"):
            run_audit(spins=(0.0, 1.0, 1.5))
