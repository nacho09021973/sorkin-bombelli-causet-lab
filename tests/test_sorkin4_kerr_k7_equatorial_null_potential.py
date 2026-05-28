"""Tests for S4-KERR-K7-EQUATORIAL-NULL-POTENTIAL-001.

Test requirements:
 1.  Artifact CSV/JSON/MD/PNG files exist.
 2.  CSV is parseable and has exactly 9 rows (one per spin).
 3.  JSON is parseable; has "aggregate" and "rows".
 4.  MD exists.
 5.  PNG exists and has nonzero size.
 6.  JSON aggregate all_checks_pass=True.
 7.  Spin list is exactly
     [0.0, 1e-4, 1e-3, 1e-2, 0.1, 0.25, 0.5, 0.75, 0.9].
 8.  At a=0: r_ph_pro = r_ph_retro = 3M, b_ph_pro = +3√3M, b_ph_retro = -3√3M.
 9.  For all spins: r_ph_pro > r_plus and r_ph_retro > r_plus.
10.  For all spins: circular_potential_pass=True, circular_derivative_pass=True.
11.  For all a>0: r_ph_pro < 3M, r_ph_retro > 3M, b_ph_pro > 0, b_ph_retro < 0.
12.  For all a>0: global causal accounting invariant
         (global_true=0, global_false=0, global_undecided=N*(N-1)/2).
13.  Existing K1–K6 and Schwarzschild stability tests still pass
     (covered by the full validation command; not duplicated here).
"""
from __future__ import annotations

import csv
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "explore" / "sorkin4_kerr_benchmark"))

import audit_kerr_k7_equatorial_null_potential_001 as audit

ARTIFACT_DIR  = ROOT / "explore" / "sorkin4_kerr_benchmark"
ARTIFACT_CSV  = ARTIFACT_DIR / "kerr_k7_equatorial_null_potential_001_n12_seed1959.csv"
ARTIFACT_JSON = ARTIFACT_DIR / "kerr_k7_equatorial_null_potential_001_n12_seed1959.json"
ARTIFACT_MD   = ARTIFACT_DIR / "kerr_k7_equatorial_null_potential_001_n12_seed1959.md"
ARTIFACT_PNG  = ARTIFACT_DIR / "kerr_k7_equatorial_null_potential_001_n12_seed1959.png"

EXPECTED_SPINS = [0.0, 1e-4, 1e-3, 1e-2, 0.1, 0.25, 0.5, 0.75, 0.9]
CIRCULAR_TOL   = 1.0e-9
SCHW_LIMIT_TOL = 1.0e-12

_MISSING_MSG = (
    "Missing K7 artifact. Run: python3 explore/sorkin4_kerr_benchmark/"
    "audit_kerr_k7_equatorial_null_potential_001.py"
)


class KerrK7ArtifactFileTests(unittest.TestCase):
    """Tests 1–7: artifact files exist, parseable, correct metadata."""

    def test_csv_exists_and_parseable(self) -> None:
        self.assertTrue(ARTIFACT_CSV.exists(), msg=_MISSING_MSG)
        with ARTIFACT_CSV.open(encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        self.assertEqual(len(rows), len(EXPECTED_SPINS))

    def test_json_exists_and_parseable(self) -> None:
        self.assertTrue(ARTIFACT_JSON.exists(), msg=_MISSING_MSG)
        with ARTIFACT_JSON.open(encoding="utf-8") as fh:
            payload = json.load(fh)
        self.assertIn("aggregate", payload)
        self.assertIn("rows", payload)

    def test_md_exists(self) -> None:
        self.assertTrue(ARTIFACT_MD.exists(), msg=_MISSING_MSG)

    def test_png_exists_and_nonempty(self) -> None:
        self.assertTrue(ARTIFACT_PNG.exists(), msg=_MISSING_MSG)
        self.assertGreater(ARTIFACT_PNG.stat().st_size, 0)

    def test_json_all_checks_pass(self) -> None:
        self.assertTrue(ARTIFACT_JSON.exists(), msg=_MISSING_MSG)
        with ARTIFACT_JSON.open(encoding="utf-8") as fh:
            payload = json.load(fh)
        self.assertTrue(
            payload["aggregate"]["all_checks_pass"],
            msg="aggregate all_checks_pass is False in frozen artifact",
        )

    def test_json_spin_list(self) -> None:
        self.assertTrue(ARTIFACT_JSON.exists(), msg=_MISSING_MSG)
        with ARTIFACT_JSON.open(encoding="utf-8") as fh:
            payload = json.load(fh)
        self.assertEqual(
            payload["aggregate"]["spins"],
            EXPECTED_SPINS,
            msg="Spin list in frozen artifact does not match expected sweep",
        )

    def test_json_row_count(self) -> None:
        self.assertTrue(ARTIFACT_JSON.exists(), msg=_MISSING_MSG)
        with ARTIFACT_JSON.open(encoding="utf-8") as fh:
            payload = json.load(fh)
        self.assertEqual(len(payload["rows"]), len(EXPECTED_SPINS))


class KerrK7ArtifactRowTests(unittest.TestCase):
    """Tests 8–12: per-row checks on the frozen JSON artifact."""

    def setUp(self) -> None:
        self.assertTrue(ARTIFACT_JSON.exists(), msg=_MISSING_MSG)
        with ARTIFACT_JSON.open(encoding="utf-8") as fh:
            self.payload = json.load(fh)
        self.rows = self.payload["rows"]
        self.n    = self.payload["aggregate"]["N"]

    # Test 8
    def test_a0_schwarzschild_photon_sphere(self) -> None:
        a0 = next(r for r in self.rows if r["spin_a"] == 0.0)
        M = a0["M"]
        self.assertAlmostEqual(a0["r_ph_pro"],   3.0 * M,              delta=SCHW_LIMIT_TOL)
        self.assertAlmostEqual(a0["r_ph_retro"], 3.0 * M,              delta=SCHW_LIMIT_TOL)
        self.assertAlmostEqual(a0["b_ph_pro"],   3.0 * math.sqrt(3.0) * M, delta=SCHW_LIMIT_TOL)
        self.assertAlmostEqual(a0["b_ph_retro"], -3.0 * math.sqrt(3.0) * M, delta=SCHW_LIMIT_TOL)

    # Test 9
    def test_photon_orbits_outside_horizon(self) -> None:
        for row in self.rows:
            self.assertTrue(
                row["r_ph_pro_outside_horizon"],
                msg=f"r_ph_pro not outside horizon for a={row['spin_a']}",
            )
            self.assertTrue(
                row["r_ph_retro_outside_horizon"],
                msg=f"r_ph_retro not outside horizon for a={row['spin_a']}",
            )
            self.assertGreater(row["r_ph_pro"],   row["r_plus"], msg=f"a={row['spin_a']}")
            self.assertGreater(row["r_ph_retro"], row["r_plus"], msg=f"a={row['spin_a']}")

    # Test 10
    def test_circular_potential_pass(self) -> None:
        for row in self.rows:
            self.assertTrue(
                row["circular_potential_pass"],
                msg=f"circular_potential_pass=False for a={row['spin_a']}",
            )

    def test_circular_derivative_pass(self) -> None:
        for row in self.rows:
            self.assertTrue(
                row["circular_derivative_pass"],
                msg=f"circular_derivative_pass=False for a={row['spin_a']}",
            )

    # Test 11
    def test_positive_spins_orbit_ordering(self) -> None:
        for row in self.rows:
            if row["spin_a"] > 0.0:
                M = row["M"]
                self.assertLess(
                    row["r_ph_pro"], 3.0 * M,
                    msg=f"r_ph_pro >= 3M for a={row['spin_a']}",
                )
                self.assertGreater(
                    row["r_ph_retro"], 3.0 * M,
                    msg=f"r_ph_retro <= 3M for a={row['spin_a']}",
                )

    def test_positive_spins_impact_parameter_signs(self) -> None:
        for row in self.rows:
            if row["spin_a"] > 0.0:
                self.assertGreater(
                    row["b_ph_pro"], 0.0,
                    msg=f"b_ph_pro <= 0 for a={row['spin_a']}",
                )
                self.assertLess(
                    row["b_ph_retro"], 0.0,
                    msg=f"b_ph_retro >= 0 for a={row['spin_a']}",
                )

    # Test 12
    def test_global_causal_accounting_positive_spins(self) -> None:
        possible_pairs = self.n * (self.n - 1) // 2
        for row in self.rows:
            if row["spin_a"] > 0.0:
                self.assertEqual(row["global_true_relations"],  0,
                                 msg=f"a={row['spin_a']}")
                self.assertEqual(row["global_false_relations"], 0,
                                 msg=f"a={row['spin_a']}")
                self.assertEqual(row["global_undecided_pairs"], possible_pairs,
                                 msg=f"a={row['spin_a']}")

    def test_positive_spin_cases_all_undecided_aggregate(self) -> None:
        self.assertTrue(self.payload["aggregate"]["positive_spin_cases_all_undecided"])


class KerrK7ComputationTests(unittest.TestCase):
    """In-process tests: run the computation without relying on frozen artifacts."""

    def setUp(self) -> None:
        self.payload   = audit.run_audit()
        self.rows      = self.payload["rows"]
        self.aggregate = self.payload["aggregate"]

    def test_all_checks_pass(self) -> None:
        self.assertTrue(self.aggregate["all_checks_pass"])

    def test_spin_list(self) -> None:
        self.assertEqual(self.aggregate["spins"], EXPECTED_SPINS)

    def test_row_count(self) -> None:
        self.assertEqual(len(self.rows), len(EXPECTED_SPINS))

    def test_n_and_mass(self) -> None:
        self.assertEqual(self.aggregate["N"], 12)
        self.assertEqual(self.aggregate["M"], 1.0)

    def test_photon_orbits_outside_horizon(self) -> None:
        for row in self.rows:
            self.assertGreater(row["r_ph_pro"],   row["r_plus"], msg=f"a={row['spin_a']}")
            self.assertGreater(row["r_ph_retro"], row["r_plus"], msg=f"a={row['spin_a']}")

    def test_circular_potential_pass_all_spins(self) -> None:
        for row in self.rows:
            self.assertTrue(row["circular_potential_pass"],  msg=f"a={row['spin_a']}")
            self.assertTrue(row["circular_derivative_pass"], msg=f"a={row['spin_a']}")

    def test_a0_schwarzschild_photon_sphere(self) -> None:
        a0 = next(r for r in self.rows if r["spin_a"] == 0.0)
        M = a0["M"]
        self.assertAlmostEqual(a0["r_ph_pro"],   3.0 * M,                  delta=SCHW_LIMIT_TOL)
        self.assertAlmostEqual(a0["r_ph_retro"], 3.0 * M,                  delta=SCHW_LIMIT_TOL)
        self.assertAlmostEqual(a0["b_ph_pro"],   3.0 * math.sqrt(3.0) * M, delta=SCHW_LIMIT_TOL)
        self.assertAlmostEqual(a0["b_ph_retro"],-3.0 * math.sqrt(3.0) * M, delta=SCHW_LIMIT_TOL)
        self.assertTrue(a0["schwarzschild_photon_sphere_pass"])
        self.assertTrue(a0["all_checks_pass"])

    def test_positive_spins_orbit_ordering(self) -> None:
        for row in self.rows:
            if row["spin_a"] > 0.0:
                self.assertLess(row["r_ph_pro"],    3.0 * row["M"], msg=f"a={row['spin_a']}")
                self.assertGreater(row["r_ph_retro"], 3.0 * row["M"], msg=f"a={row['spin_a']}")

    def test_positive_spins_impact_signs(self) -> None:
        for row in self.rows:
            if row["spin_a"] > 0.0:
                self.assertGreater(row["b_ph_pro"],   0.0, msg=f"a={row['spin_a']}")
                self.assertLess(row["b_ph_retro"],    0.0, msg=f"a={row['spin_a']}")

    def test_global_causal_accounting_positive_spins(self) -> None:
        n = self.aggregate["N"]
        possible = n * (n - 1) // 2
        for row in self.rows:
            if row["spin_a"] > 0.0:
                self.assertEqual(row["global_true_relations"],  0,       msg=f"a={row['spin_a']}")
                self.assertEqual(row["global_false_relations"], 0,       msg=f"a={row['spin_a']}")
                self.assertEqual(row["global_undecided_pairs"], possible, msg=f"a={row['spin_a']}")

    def test_positive_spin_cases_all_undecided(self) -> None:
        self.assertTrue(self.aggregate["positive_spin_cases_all_undecided"])

    # --- Unit tests on core physics helpers ---

    def test_null_radial_potential_schw_a0(self) -> None:
        """At a=0, r=3M, b=3√3M: R=0 exactly."""
        M = 1.0
        b = 3.0 * math.sqrt(3.0) * M
        R = audit.null_radial_potential(3.0 * M, 0.0, b, M)
        self.assertAlmostEqual(R, 0.0, delta=1.0e-14, msg="R != 0 at Schwarzschild photon sphere")

    def test_null_radial_potential_derivative_schw_a0(self) -> None:
        """At a=0, r=3M, b=3√3M: dR/dr=0 exactly."""
        M = 1.0
        b = 3.0 * math.sqrt(3.0) * M
        dR = audit.null_radial_potential_derivative(3.0 * M, 0.0, b, M)
        self.assertAlmostEqual(dR, 0.0, delta=1.0e-14, msg="dR/dr != 0 at Schwarzschild photon sphere")

    def test_photon_sphere_radius_pro_a0(self) -> None:
        """Prograde photon orbit at a=0 is exactly 3M."""
        M = 1.0
        r = audit.photon_sphere_radius_pro(M, 0.0)
        self.assertAlmostEqual(r, 3.0 * M, delta=SCHW_LIMIT_TOL)

    def test_photon_sphere_radius_retro_a0(self) -> None:
        """Retrograde photon orbit at a=0 is exactly 3M."""
        M = 1.0
        r = audit.photon_sphere_radius_retro(M, 0.0)
        self.assertAlmostEqual(r, 3.0 * M, delta=SCHW_LIMIT_TOL)

    def test_impact_parameter_pro_a0(self) -> None:
        """Prograde impact parameter at a=0, r=3M: b = +3√3M."""
        M = 1.0
        b = audit.photon_impact_parameter(3.0 * M, M, 0.0, prograde=True)
        self.assertAlmostEqual(b, 3.0 * math.sqrt(3.0) * M, delta=SCHW_LIMIT_TOL)

    def test_impact_parameter_retro_a0(self) -> None:
        """Retrograde impact parameter at a=0, r=3M: b = -3√3M."""
        M = 1.0
        b = audit.photon_impact_parameter(3.0 * M, M, 0.0, prograde=False)
        self.assertAlmostEqual(b, -3.0 * math.sqrt(3.0) * M, delta=SCHW_LIMIT_TOL)

    def test_r_and_dR_at_circular_orbits_all_spins(self) -> None:
        """R and dR/dr vanish at circular orbit radii for all spins."""
        M = 1.0
        for a in EXPECTED_SPINS:
            r_ph_pro   = audit.photon_sphere_radius_pro(M, a)
            r_ph_retro = audit.photon_sphere_radius_retro(M, a)
            b_pro   = audit.photon_impact_parameter(r_ph_pro,   M, a, prograde=True)
            b_retro = audit.photon_impact_parameter(r_ph_retro, M, a, prograde=False)

            R_pro    = audit.null_radial_potential(r_ph_pro,   a, b_pro,   M)
            dR_pro   = audit.null_radial_potential_derivative(r_ph_pro,   a, b_pro,   M)
            R_retro  = audit.null_radial_potential(r_ph_retro, a, b_retro, M)
            dR_retro = audit.null_radial_potential_derivative(r_ph_retro, a, b_retro, M)

            self.assertLessEqual(abs(R_pro),    CIRCULAR_TOL, msg=f"|R_pro| > tol at a={a}")
            self.assertLessEqual(abs(dR_pro),   CIRCULAR_TOL, msg=f"|dR_pro| > tol at a={a}")
            self.assertLessEqual(abs(R_retro),  CIRCULAR_TOL, msg=f"|R_retro| > tol at a={a}")
            self.assertLessEqual(abs(dR_retro), CIRCULAR_TOL, msg=f"|dR_retro| > tol at a={a}")

    def test_prograde_orbit_decreases_with_spin(self) -> None:
        """r_ph_pro decreases monotonically with a > 0."""
        M = 1.0
        positive_spins = [a for a in EXPECTED_SPINS if a > 0.0]
        r_values = [audit.photon_sphere_radius_pro(M, a) for a in positive_spins]
        for i in range(len(r_values) - 1):
            self.assertGreater(
                r_values[i], r_values[i + 1],
                msg=f"r_ph_pro not decreasing at a={positive_spins[i+1]}",
            )

    def test_retrograde_orbit_increases_with_spin(self) -> None:
        """r_ph_retro increases monotonically with a > 0."""
        M = 1.0
        positive_spins = [a for a in EXPECTED_SPINS if a > 0.0]
        r_values = [audit.photon_sphere_radius_retro(M, a) for a in positive_spins]
        for i in range(len(r_values) - 1):
            self.assertLess(
                r_values[i], r_values[i + 1],
                msg=f"r_ph_retro not increasing at a={positive_spins[i+1]}",
            )

    def test_photon_impact_parameter_raises_inside_horizon(self) -> None:
        """photon_impact_parameter raises ValueError when Delta <= 0."""
        M = 1.0
        a = 0.5
        r_plus = M + math.sqrt(M * M - a * a)
        r_inside = r_plus - 0.1
        with self.assertRaises(ValueError):
            audit.photon_impact_parameter(r_inside, M, a, prograde=True)

    def test_null_potential_positive_far_exterior(self) -> None:
        """R(r; a, b=0) > 0 for r >> M, since [r^2+a^2]^2 > Delta*a^2."""
        M = 1.0
        for a in [0.0, 0.5, 0.9]:
            r = 100.0 * M
            R = audit.null_radial_potential(r, a, 0.0, M)
            self.assertGreater(R, 0.0, msg=f"R not positive at r=100M, a={a}, b=0")


if __name__ == "__main__":
    unittest.main()
