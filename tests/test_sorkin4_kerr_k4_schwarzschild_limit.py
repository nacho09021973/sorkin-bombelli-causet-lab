"""Tests for S4-KERR-K4-SCHWARZSCHILD-LIMIT-001.

Test requirements (see spec):
 1.  Artifact CSV/JSON/MD/PNG files exist.
 2.  CSV and JSON are parseable.
 3.  PNG exists and has nonzero file size.
 4.  JSON summary has all_checks_pass=True.
 5.  The spin list is exactly [0.0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1].
 6.  a=0 passes Schwarzschild metric limit:
         g_tphi = 0, g_rr = 1/(1-2M/r), g_phiphi = r^2.
 7.  For all a>0, frame_dragging_linear_check_pass=True.
 8.  For all a>0, gphiphi_quadratic_check_pass=True.
 9.  For all spins, grr_formula_check_pass=True.
10.  For small a <= 1e-2, horizon_quadratic_check_pass=True.
11.  For all spins, all_points_exterior=True.
12.  For all a>0:
         global_true_relations=0
         global_false_relations=0
         global_undecided_pairs=N*(N-1)/2
13.  Existing K3, K2, K1, and Schwarzschild stability tests still pass
     (covered by the full validation command; not duplicated here).
"""
from __future__ import annotations

import csv
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "explore" / "sorkin4_kerr_benchmark"))

import audit_kerr_k4_schwarzschild_limit_001 as audit

ARTIFACT_DIR = ROOT / "explore" / "sorkin4_kerr_benchmark"
ARTIFACT_CSV  = ARTIFACT_DIR / "kerr_k4_schwarzschild_limit_001_n12_seed1959.csv"
ARTIFACT_JSON = ARTIFACT_DIR / "kerr_k4_schwarzschild_limit_001_n12_seed1959.json"
ARTIFACT_MD   = ARTIFACT_DIR / "kerr_k4_schwarzschild_limit_001_n12_seed1959.md"
ARTIFACT_PNG  = ARTIFACT_DIR / "kerr_k4_schwarzschild_limit_001_n12_seed1959.png"

EXPECTED_SPINS = [0.0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1]
HORIZON_QUAD_SPIN_MAX = 1e-2
METRIC_TOL = 1.0e-12


class KerrK4ArtifactFileTests(unittest.TestCase):
    """Tests 1–5: artifact files exist and are parseable."""

    def _missing_msg(self, path: Path) -> str:
        return (
            f"Missing K4 artifact {path}. Run: "
            "python3 explore/sorkin4_kerr_benchmark/"
            "audit_kerr_k4_schwarzschild_limit_001.py"
        )

    # Test 1 & 2: files exist
    def test_csv_exists_and_parseable(self) -> None:
        self.assertTrue(ARTIFACT_CSV.exists(), msg=self._missing_msg(ARTIFACT_CSV))
        with ARTIFACT_CSV.open(encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        self.assertEqual(len(rows), len(EXPECTED_SPINS))

    def test_json_exists_and_parseable(self) -> None:
        self.assertTrue(ARTIFACT_JSON.exists(), msg=self._missing_msg(ARTIFACT_JSON))
        with ARTIFACT_JSON.open(encoding="utf-8") as fh:
            payload = json.load(fh)
        self.assertIn("aggregate", payload)
        self.assertIn("rows", payload)

    def test_md_exists(self) -> None:
        self.assertTrue(ARTIFACT_MD.exists(), msg=self._missing_msg(ARTIFACT_MD))

    # Test 3: PNG exists and has nonzero size
    def test_png_exists_and_nonempty(self) -> None:
        self.assertTrue(ARTIFACT_PNG.exists(), msg=self._missing_msg(ARTIFACT_PNG))
        self.assertGreater(ARTIFACT_PNG.stat().st_size, 0)

    # Test 4: JSON all_checks_pass=True
    def test_json_all_checks_pass(self) -> None:
        with ARTIFACT_JSON.open(encoding="utf-8") as fh:
            payload = json.load(fh)
        self.assertTrue(
            payload["aggregate"]["all_checks_pass"],
            msg="aggregate all_checks_pass is False in frozen artifact",
        )

    # Test 5: correct spin list
    def test_json_spin_list(self) -> None:
        with ARTIFACT_JSON.open(encoding="utf-8") as fh:
            payload = json.load(fh)
        self.assertEqual(
            payload["aggregate"]["spins"],
            EXPECTED_SPINS,
            msg="Spin list in frozen artifact does not match expected sweep",
        )


class KerrK4ArtifactRowTests(unittest.TestCase):
    """Tests 6–12: per-row checks on the frozen JSON artifact."""

    def setUp(self) -> None:
        self.assertTrue(
            ARTIFACT_JSON.exists(),
            msg=(
                "Missing K4 artifact. Run: "
                "python3 explore/sorkin4_kerr_benchmark/"
                "audit_kerr_k4_schwarzschild_limit_001.py"
            ),
        )
        with ARTIFACT_JSON.open(encoding="utf-8") as fh:
            self.payload = json.load(fh)
        self.rows = self.payload["rows"]
        self.n    = self.payload["aggregate"]["N"]

    def test_row_count(self) -> None:
        self.assertEqual(len(self.rows), len(EXPECTED_SPINS))

    # Test 11: all_points_exterior for all spins
    def test_all_points_exterior(self) -> None:
        for row in self.rows:
            self.assertTrue(
                row["all_points_exterior"],
                msg=f"all_points_exterior=False for a={row['spin_a']}",
            )

    # Test 9: grr_formula_check_pass for all spins
    def test_grr_formula_check_pass_all_spins(self) -> None:
        for row in self.rows:
            self.assertTrue(
                row["grr_formula_check_pass"],
                msg=f"grr_formula_check_pass=False for a={row['spin_a']}",
            )

    # Test 6: a=0 Schwarzschild metric limit
    def test_a0_schwarzschild_metric_limit(self) -> None:
        a0_rows = [r for r in self.rows if r["spin_a"] == 0.0]
        self.assertEqual(len(a0_rows), 1, msg="Expected exactly one a=0 row")
        row = a0_rows[0]
        self.assertTrue(
            row["schwarzschild_metric_limit_pass"],
            msg="schwarzschild_metric_limit_pass=False at a=0",
        )
        self.assertLessEqual(
            row["max_abs_gtphi_at_a0"],
            METRIC_TOL,
            msg=f"max_abs_gtphi_at_a0={row['max_abs_gtphi_at_a0']} > {METRIC_TOL}",
        )
        self.assertLessEqual(
            row["max_abs_grr_schwarzschild_at_a0"],
            METRIC_TOL,
            msg=(
                f"max_abs_grr_schwarzschild_at_a0="
                f"{row['max_abs_grr_schwarzschild_at_a0']} > {METRIC_TOL}"
            ),
        )
        self.assertLessEqual(
            row["max_abs_gphiphi_minus_r2_at_a0"],
            METRIC_TOL,
            msg=(
                f"max_abs_gphiphi_minus_r2_at_a0="
                f"{row['max_abs_gphiphi_minus_r2_at_a0']} > {METRIC_TOL}"
            ),
        )

    # Test 7: frame_dragging_linear_check_pass for all a>0
    def test_frame_dragging_linear_check_positive_spins(self) -> None:
        for row in self.rows:
            if row["spin_a"] > 0.0:
                self.assertTrue(
                    row["frame_dragging_linear_check_pass"],
                    msg=f"frame_dragging_linear_check_pass=False for a={row['spin_a']}",
                )

    # Test 8: gphiphi_quadratic_check_pass for all a>0
    def test_gphiphi_quadratic_check_positive_spins(self) -> None:
        for row in self.rows:
            if row["spin_a"] > 0.0:
                self.assertTrue(
                    row["gphiphi_quadratic_check_pass"],
                    msg=f"gphiphi_quadratic_check_pass=False for a={row['spin_a']}",
                )

    # Test 10: horizon_quadratic_check_pass for small a <= 1e-2
    def test_horizon_quadratic_check_small_spins(self) -> None:
        for row in self.rows:
            a = row["spin_a"]
            if 0.0 < a <= HORIZON_QUAD_SPIN_MAX:
                self.assertTrue(
                    row["horizon_quadratic_check_pass"],
                    msg=f"horizon_quadratic_check_pass=False for a={a}",
                )

    # Test 12: global causal accounting for a>0
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


class KerrK4ComputationTests(unittest.TestCase):
    """In-process tests: run the computation without frozen artifacts."""

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

    def test_radial_grid(self) -> None:
        self.assertEqual(
            self.aggregate["radial_grid"],
            list(audit.RADIAL_GRID),
        )

    def test_all_points_exterior(self) -> None:
        for row in self.rows:
            self.assertTrue(
                row["all_points_exterior"],
                msg=f"all_points_exterior=False for a={row['spin_a']}",
            )

    def test_grr_formula_check_all_spins(self) -> None:
        for row in self.rows:
            self.assertTrue(
                row["grr_formula_check_pass"],
                msg=f"grr_formula_check_pass=False for a={row['spin_a']}",
            )

    def test_a0_schwarzschild_limit(self) -> None:
        a0 = next(r for r in self.rows if r["spin_a"] == 0.0)
        self.assertTrue(a0["schwarzschild_metric_limit_pass"])
        self.assertLessEqual(a0["max_abs_gtphi_at_a0"],           METRIC_TOL)
        self.assertLessEqual(a0["max_abs_grr_schwarzschild_at_a0"], METRIC_TOL)
        self.assertLessEqual(a0["max_abs_gphiphi_minus_r2_at_a0"], METRIC_TOL)

    def test_frame_dragging_linear_positive_spins(self) -> None:
        for row in self.rows:
            if row["spin_a"] > 0.0:
                self.assertTrue(
                    row["frame_dragging_linear_check_pass"],
                    msg=f"a={row['spin_a']}",
                )

    def test_gphiphi_quadratic_positive_spins(self) -> None:
        for row in self.rows:
            if row["spin_a"] > 0.0:
                self.assertTrue(
                    row["gphiphi_quadratic_check_pass"],
                    msg=f"a={row['spin_a']}",
                )

    def test_horizon_quadratic_small_spins(self) -> None:
        for row in self.rows:
            a = row["spin_a"]
            if 0.0 < a <= HORIZON_QUAD_SPIN_MAX:
                self.assertTrue(
                    row["horizon_quadratic_check_pass"],
                    msg=f"a={a}",
                )

    def test_global_causal_accounting_positive_spins(self) -> None:
        n = self.aggregate["N"]
        possible = n * (n - 1) // 2
        for row in self.rows:
            if row["spin_a"] > 0.0:
                self.assertEqual(row["global_true_relations"],  0,        msg=f"a={row['spin_a']}")
                self.assertEqual(row["global_false_relations"], 0,        msg=f"a={row['spin_a']}")
                self.assertEqual(row["global_undecided_pairs"], possible, msg=f"a={row['spin_a']}")

    def test_equatorial_metric_schwarzschild_a0(self) -> None:
        """equatorial_metric_at_r at a=0 returns exact Schwarzschild coefficients."""
        M = 1.0
        for r in [2.5, 3.0, 4.0, 6.0, 10.0]:
            m = audit.equatorial_metric_at_r(r, M, 0.0)
            self.assertAlmostEqual(
                m["g_tt"], -(1.0 - 2*M/r), places=14,
                msg=f"g_tt mismatch at r={r}",
            )
            self.assertAlmostEqual(
                m["g_tphi"], 0.0, places=14,
                msg=f"g_tphi != 0 at r={r}, a=0",
            )
            self.assertAlmostEqual(
                m["g_rr"], 1.0/(1.0 - 2*M/r), places=14,
                msg=f"g_rr mismatch at r={r}",
            )
            self.assertAlmostEqual(
                m["g_phiphi"], r*r, places=14,
                msg=f"g_phiphi != r^2 at r={r}, a=0",
            )

    def test_equatorial_metric_frame_dragging_sign(self) -> None:
        """g_tphi < 0 for a > 0 at all exterior radii."""
        M = 1.0
        for a in [1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1]:
            for r in [2.5, 3.0, 4.0, 6.0, 10.0]:
                m = audit.equatorial_metric_at_r(r, M, a)
                self.assertLess(
                    m["g_tphi"], 0.0,
                    msg=f"g_tphi >= 0 for a={a}, r={r}",
                )

    def test_positive_spin_cases_all_undecided(self) -> None:
        self.assertTrue(self.aggregate["positive_spin_cases_all_undecided"])

    def test_a0_schwarzschild_fields_null_for_positive_spins(self) -> None:
        """a>0 rows must have None for a=0-specific fields."""
        for row in self.rows:
            if row["spin_a"] > 0.0:
                self.assertIsNone(row["max_abs_gtphi_at_a0"],             msg=f"a={row['spin_a']}")
                self.assertIsNone(row["max_abs_gphiphi_minus_r2_at_a0"],  msg=f"a={row['spin_a']}")
                self.assertIsNone(row["max_abs_grr_schwarzschild_at_a0"], msg=f"a={row['spin_a']}")

    def test_r_plus_shift_over_a2_is_none_at_a0(self) -> None:
        a0 = next(r for r in self.rows if r["spin_a"] == 0.0)
        self.assertIsNone(a0["r_plus_shift_over_a2"])

    def test_horizon_shift_converges_to_half(self) -> None:
        """r_plus_shift_over_a2 → 1/2 as a → 0 (for a > 0)."""
        for row in self.rows:
            a = row["spin_a"]
            if 0.0 < a <= HORIZON_QUAD_SPIN_MAX:
                self.assertAlmostEqual(
                    row["r_plus_shift_over_a2"],
                    0.5,
                    delta=audit.HORIZON_QUAD_TOL,
                    msg=f"horizon shift not near 0.5 for a={a}",
                )


if __name__ == "__main__":
    unittest.main()
