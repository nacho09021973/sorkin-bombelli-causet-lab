"""Tests for S4-KERR-K5-PROGRADE-RETROGRADE-LOCAL-CONE-001.

Test requirements (see spec):
 1.  Artifact CSV/JSON/MD/PNG files exist.
 2.  CSV and JSON are parseable.
 3.  PNG exists and has nonzero file size.
 4.  JSON summary has all_checks_pass=True.
 5.  The spin list is exactly
     [0.0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 0.25, 0.5, 0.75].
 6.  The radial grid is exactly [2.5, 3.0, 4.0, 6.0, 10.0].
 7.  For all spins, all_points_exterior=True.
 8.  For all spins, discriminant_positive_pass=True.
 9.  For all spins, omega_width_positive_pass=True.
10.  For a=0, schwarzschild_symmetry_pass=True and
     max_abs_schwarzschild_symmetry_residual is near zero.
11.  For a=0, frame_dragging_sign_pass=True (True by convention, must not
     fail all_checks_pass).
12.  For all a>0, frame_dragging_sign_pass=True and omega_center_mean > 0.
13.  For small positive a <= 1e-2, small_a_linear_scaling_pass=True.
14.  For all a>0:
         global_true_relations=0
         global_false_relations=0
         global_undecided_pairs=N*(N-1)/2
15.  Existing K4, K3, K2, K1, L0, and Schwarzschild stability tests still pass
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

import audit_kerr_k5_prograde_retrograde_local_cone_001 as audit

ARTIFACT_DIR  = ROOT / "explore" / "sorkin4_kerr_benchmark"
ARTIFACT_CSV  = ARTIFACT_DIR / "kerr_k5_prograde_retrograde_local_cone_001_n12_seed1959.csv"
ARTIFACT_JSON = ARTIFACT_DIR / "kerr_k5_prograde_retrograde_local_cone_001_n12_seed1959.json"
ARTIFACT_MD   = ARTIFACT_DIR / "kerr_k5_prograde_retrograde_local_cone_001_n12_seed1959.md"
ARTIFACT_PNG  = ARTIFACT_DIR / "kerr_k5_prograde_retrograde_local_cone_001_n12_seed1959.png"

EXPECTED_SPINS = [0.0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 0.25, 0.5, 0.75]
EXPECTED_GRID  = [2.5, 3.0, 4.0, 6.0, 10.0]
SMALL_A_MAX    = 1e-2
METRIC_TOL     = 1.0e-12
SMALL_A_TOL    = 1.0e-4

_MISSING_MSG = (
    "Missing K5 artifact. Run: python3 explore/sorkin4_kerr_benchmark/"
    "audit_kerr_k5_prograde_retrograde_local_cone_001.py"
)


class KerrK5ArtifactFileTests(unittest.TestCase):
    """Tests 1–6: artifact files exist, are parseable, and have correct metadata."""

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

    def test_json_radial_grid(self) -> None:
        self.assertTrue(ARTIFACT_JSON.exists(), msg=_MISSING_MSG)
        with ARTIFACT_JSON.open(encoding="utf-8") as fh:
            payload = json.load(fh)
        self.assertEqual(
            payload["aggregate"]["radial_grid"],
            EXPECTED_GRID,
            msg="Radial grid in frozen artifact does not match expected grid",
        )


class KerrK5ArtifactRowTests(unittest.TestCase):
    """Tests 7–14: per-row checks on the frozen JSON artifact."""

    def setUp(self) -> None:
        self.assertTrue(ARTIFACT_JSON.exists(), msg=_MISSING_MSG)
        with ARTIFACT_JSON.open(encoding="utf-8") as fh:
            self.payload = json.load(fh)
        self.rows = self.payload["rows"]
        self.n    = self.payload["aggregate"]["N"]

    def test_row_count(self) -> None:
        self.assertEqual(len(self.rows), len(EXPECTED_SPINS))

    # Test 7
    def test_all_points_exterior(self) -> None:
        for row in self.rows:
            self.assertTrue(
                row["all_points_exterior"],
                msg=f"all_points_exterior=False for a={row['spin_a']}",
            )

    # Test 8
    def test_discriminant_positive_pass(self) -> None:
        for row in self.rows:
            self.assertTrue(
                row["discriminant_positive_pass"],
                msg=f"discriminant_positive_pass=False for a={row['spin_a']}",
            )

    # Test 9
    def test_omega_width_positive_pass(self) -> None:
        for row in self.rows:
            self.assertTrue(
                row["omega_width_positive_pass"],
                msg=f"omega_width_positive_pass=False for a={row['spin_a']}",
            )

    # Test 10
    def test_a0_schwarzschild_symmetry(self) -> None:
        a0_rows = [r for r in self.rows if r["spin_a"] == 0.0]
        self.assertEqual(len(a0_rows), 1)
        row = a0_rows[0]
        self.assertTrue(
            row["schwarzschild_symmetry_pass"],
            msg="schwarzschild_symmetry_pass=False at a=0",
        )
        self.assertLessEqual(
            row["max_abs_schwarzschild_symmetry_residual"],
            METRIC_TOL,
            msg=f"symmetry residual={row['max_abs_schwarzschild_symmetry_residual']} > {METRIC_TOL}",
        )

    # Test 11
    def test_a0_frame_dragging_sign_pass_by_convention(self) -> None:
        a0_rows = [r for r in self.rows if r["spin_a"] == 0.0]
        self.assertEqual(len(a0_rows), 1)
        row = a0_rows[0]
        self.assertTrue(
            row["frame_dragging_sign_pass"],
            msg="frame_dragging_sign_pass must be True at a=0 (by convention)",
        )
        self.assertTrue(
            row["all_checks_pass"],
            msg="all_checks_pass must be True at a=0",
        )

    # Test 12
    def test_positive_spins_frame_dragging_sign(self) -> None:
        for row in self.rows:
            if row["spin_a"] > 0.0:
                self.assertTrue(
                    row["frame_dragging_sign_pass"],
                    msg=f"frame_dragging_sign_pass=False for a={row['spin_a']}",
                )

    def test_positive_spins_omega_center_mean_positive(self) -> None:
        for row in self.rows:
            if row["spin_a"] > 0.0:
                self.assertGreater(
                    row["omega_center_mean"],
                    0.0,
                    msg=f"omega_center_mean <= 0 for a={row['spin_a']}",
                )

    # Test 13
    def test_small_a_linear_scaling_pass(self) -> None:
        for row in self.rows:
            a = row["spin_a"]
            if 0.0 < a <= SMALL_A_MAX:
                self.assertTrue(
                    row["small_a_linear_scaling_pass"],
                    msg=f"small_a_linear_scaling_pass=False for a={a}",
                )

    # Test 14
    def test_global_causal_accounting_positive_spins(self) -> None:
        possible_pairs = self.n * (self.n - 1) // 2
        for row in self.rows:
            if row["spin_a"] > 0.0:
                self.assertEqual(
                    row["global_true_relations"], 0,
                    msg=f"global_true_relations != 0 for a={row['spin_a']}",
                )
                self.assertEqual(
                    row["global_false_relations"], 0,
                    msg=f"global_false_relations != 0 for a={row['spin_a']}",
                )
                self.assertEqual(
                    row["global_undecided_pairs"], possible_pairs,
                    msg=f"global_undecided_pairs != {possible_pairs} for a={row['spin_a']}",
                )


class KerrK5ComputationTests(unittest.TestCase):
    """In-process tests: run the computation without relying on frozen artifacts."""

    def setUp(self) -> None:
        self.payload   = audit.run_audit()
        self.rows      = self.payload["rows"]
        self.aggregate = self.payload["aggregate"]

    def test_all_checks_pass(self) -> None:
        self.assertTrue(self.aggregate["all_checks_pass"])

    def test_spin_list(self) -> None:
        self.assertEqual(self.aggregate["spins"], EXPECTED_SPINS)

    def test_radial_grid(self) -> None:
        self.assertEqual(self.aggregate["radial_grid"], EXPECTED_GRID)

    def test_row_count(self) -> None:
        self.assertEqual(len(self.rows), len(EXPECTED_SPINS))

    def test_n_and_mass(self) -> None:
        self.assertEqual(self.aggregate["N"], 12)
        self.assertEqual(self.aggregate["M"], 1.0)

    def test_all_points_exterior(self) -> None:
        for row in self.rows:
            self.assertTrue(
                row["all_points_exterior"],
                msg=f"a={row['spin_a']}",
            )

    def test_discriminant_positive_all_spins(self) -> None:
        for row in self.rows:
            self.assertTrue(
                row["discriminant_positive_pass"],
                msg=f"a={row['spin_a']}",
            )

    def test_omega_width_positive_all_spins(self) -> None:
        for row in self.rows:
            self.assertTrue(
                row["omega_width_positive_pass"],
                msg=f"a={row['spin_a']}",
            )

    def test_omega_center_exact_identity_all_spins(self) -> None:
        for row in self.rows:
            self.assertTrue(
                row["omega_center_exact_identity_pass"],
                msg=f"a={row['spin_a']}",
            )

    def test_a0_schwarzschild_symmetry(self) -> None:
        a0 = next(r for r in self.rows if r["spin_a"] == 0.0)
        self.assertTrue(a0["schwarzschild_symmetry_pass"])
        self.assertLessEqual(a0["max_abs_schwarzschild_symmetry_residual"], METRIC_TOL)

    def test_a0_omega_center_is_zero(self) -> None:
        a0 = next(r for r in self.rows if r["spin_a"] == 0.0)
        self.assertAlmostEqual(a0["omega_center_mean"], 0.0, places=13)
        self.assertAlmostEqual(a0["omega_center_max_abs"], 0.0, places=13)

    def test_a0_frame_dragging_sign_pass_by_convention(self) -> None:
        a0 = next(r for r in self.rows if r["spin_a"] == 0.0)
        self.assertTrue(a0["frame_dragging_sign_pass"])
        self.assertTrue(a0["all_checks_pass"])

    def test_positive_spins_frame_dragging_sign(self) -> None:
        for row in self.rows:
            if row["spin_a"] > 0.0:
                self.assertTrue(
                    row["frame_dragging_sign_pass"],
                    msg=f"a={row['spin_a']}",
                )

    def test_positive_spins_omega_center_positive(self) -> None:
        for row in self.rows:
            if row["spin_a"] > 0.0:
                self.assertGreater(
                    row["omega_center_mean"], 0.0,
                    msg=f"a={row['spin_a']}",
                )

    def test_small_a_linear_scaling(self) -> None:
        for row in self.rows:
            a = row["spin_a"]
            if 0.0 < a <= SMALL_A_MAX:
                self.assertTrue(
                    row["small_a_linear_scaling_pass"],
                    msg=f"a={a}",
                )
                self.assertLessEqual(
                    row["max_abs_omega_center_linear_residual_small_a"],
                    SMALL_A_TOL,
                    msg=f"a={a}",
                )

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

    def test_a0_schw_fields_null_for_positive_spins(self) -> None:
        for row in self.rows:
            if row["spin_a"] > 0.0:
                self.assertIsNone(
                    row["max_abs_schwarzschild_symmetry_residual"],
                    msg=f"a={row['spin_a']}",
                )
                self.assertIsNone(
                    row["max_abs_width_schwarzschild_residual_at_a0"],
                    msg=f"a={row['spin_a']}",
                )

    def test_small_a_field_null_for_large_spins(self) -> None:
        for row in self.rows:
            a = row["spin_a"]
            if a == 0.0 or a > SMALL_A_MAX:
                self.assertIsNone(
                    row["max_abs_omega_center_linear_residual_small_a"],
                    msg=f"a={a}",
                )

    # --- Unit tests on angular_null_slopes helper ---

    def test_angular_null_slopes_a0_schwarzschild(self) -> None:
        """At a=0, omega_+ = sqrt(1-2M/r)/r and omega_- = -omega_+."""
        M = 1.0
        for r in [2.5, 3.0, 4.0, 6.0, 10.0]:
            op, om, disc, _ = audit.angular_null_slopes(r, M, 0.0)
            expected_width = math.sqrt(1.0 - 2.0 * M / r) / r
            self.assertAlmostEqual(op,  expected_width, places=13, msg=f"r={r}")
            self.assertAlmostEqual(om, -expected_width, places=13, msg=f"r={r}")

    def test_angular_null_slopes_symmetry_a0(self) -> None:
        """omega_+ + omega_- = 0 exactly at a=0."""
        M = 1.0
        for r in [2.5, 3.0, 4.0, 6.0, 10.0]:
            op, om, _, _ = audit.angular_null_slopes(r, M, 0.0)
            self.assertAlmostEqual(op + om, 0.0, places=13, msg=f"r={r}")

    def test_angular_null_slopes_center_identity(self) -> None:
        """(omega_+ + omega_-)/2 = -g_tphi/g_phiphi for all test cases."""
        M = 1.0
        for a in [0.0, 1e-4, 1e-2, 0.25, 0.5, 0.75]:
            for r in [2.5, 3.0, 4.0, 6.0, 10.0]:
                op, om, _, m = audit.angular_null_slopes(r, M, a)
                center_formula = (op + om) * 0.5
                center_exact   = -m["g_tphi"] / m["g_phiphi"]
                self.assertAlmostEqual(
                    center_formula, center_exact, places=12,
                    msg=f"center identity fails at r={r}, a={a}",
                )

    def test_angular_null_slopes_positive_spin_prograde(self) -> None:
        """For a>0, omega_center > 0 at all exterior radii."""
        M = 1.0
        for a in [1e-4, 1e-2, 0.25, 0.5, 0.75]:
            for r in [2.5, 3.0, 4.0, 6.0, 10.0]:
                op, om, _, _ = audit.angular_null_slopes(r, M, a)
                center = (op + om) * 0.5
                self.assertGreater(
                    center, 0.0,
                    msg=f"omega_center <= 0 for a={a}, r={r}",
                )

    def test_discriminant_always_positive(self) -> None:
        """Discriminant must be positive for all exterior points."""
        M = 1.0
        for a in [0.0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 0.25, 0.5, 0.75]:
            for r in [2.5, 3.0, 4.0, 6.0, 10.0]:
                _, _, disc, _ = audit.angular_null_slopes(r, M, a)
                self.assertGreater(disc, 0.0, msg=f"a={a}, r={r}")

    def test_omega_plus_gt_omega_minus(self) -> None:
        """omega_+ > omega_- for all test cases."""
        M = 1.0
        for a in [0.0, 1e-4, 1e-2, 0.25, 0.5, 0.75]:
            for r in [2.5, 3.0, 4.0, 6.0, 10.0]:
                op, om, _, _ = audit.angular_null_slopes(r, M, a)
                self.assertGreater(op, om, msg=f"a={a}, r={r}")


if __name__ == "__main__":
    unittest.main()
