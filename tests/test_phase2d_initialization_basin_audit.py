"""Regression tests for the Phase 2D initialization / basin audit."""

import csv
import math
import os
import unittest

CSV_PATH = os.path.join(
    os.path.dirname(__file__),
    "..",
    "benchmarks",
    "foundation",
    "phase2d_initialization_basin_audit.csv",
)

EXPECTED_COLUMNS = {
    "family",
    "target_dim",
    "n",
    "seed",
    "init_label",
    "noise_epsilon",
    "initial_energy",
    "final_energy",
    "delta_energy",
    "initial_interval_rmse",
    "final_interval_rmse",
    "initial_distance_to_truth_rms",
    "final_distance_to_truth_rms",
    "improved_energy",
    "improved_interval_rmse",
    "preserved_near_truth",
    "notes",
}

EXPECTED_INIT_LABELS = {"truth", "truth_plus_small_noise", "truth_plus_medium_noise", "random_init"}


def _load_rows():
    with open(CSV_PATH, newline="") as fh:
        return list(csv.DictReader(fh))


class TestPhase2DOutputExists(unittest.TestCase):
    def test_csv_exists(self):
        self.assertTrue(os.path.exists(CSV_PATH), f"CSV not found at {CSV_PATH}")

    def test_md_exists(self):
        md_path = CSV_PATH.replace(".csv", ".md")
        self.assertTrue(os.path.exists(md_path), f"Markdown not found at {md_path}")


class TestPhase2DSchema(unittest.TestCase):
    def setUp(self):
        self.rows = _load_rows()

    def test_columns_match_expected(self):
        actual = set(self.rows[0].keys())
        self.assertEqual(actual, EXPECTED_COLUMNS)

    def test_row_count(self):
        # 3d x 2n x 3seeds x 4init_labels = 72
        self.assertEqual(len(self.rows), 72)


class TestPhase2DFamilyFilter(unittest.TestCase):
    def setUp(self):
        self.rows = _load_rows()

    def test_only_minkowski(self):
        families = {r["family"] for r in self.rows}
        self.assertEqual(families, {"minkowski"})

    def test_no_kr_or_corona(self):
        for r in self.rows:
            self.assertNotIn("kleitman", r["family"].lower())
            self.assertNotIn("corona", r["family"].lower())

    def test_dims_cover_2_3_4(self):
        dims = {int(r["target_dim"]) for r in self.rows}
        self.assertEqual(dims, {2, 3, 4})

    def test_sizes_cover_32_and_64(self):
        sizes = {int(r["n"]) for r in self.rows}
        self.assertEqual(sizes, {32, 64})


class TestPhase2DInitLabels(unittest.TestCase):
    def setUp(self):
        self.rows = _load_rows()

    def test_all_init_labels_present(self):
        labels = {r["init_label"] for r in self.rows}
        self.assertEqual(labels, EXPECTED_INIT_LABELS)

    def test_truth_rows_present(self):
        truth_rows = [r for r in self.rows if r["init_label"] == "truth"]
        self.assertGreater(len(truth_rows), 0)

    def test_each_label_has_18_rows(self):
        for label in EXPECTED_INIT_LABELS:
            count = sum(1 for r in self.rows if r["init_label"] == label)
            self.assertEqual(count, 18, f"Expected 18 rows for {label}, got {count}")


class TestPhase2DTruthInit(unittest.TestCase):
    def setUp(self):
        self.truth_rows = [r for r in _load_rows() if r["init_label"] == "truth"]

    def test_truth_initial_energy_near_zero(self):
        for r in self.truth_rows:
            self.assertAlmostEqual(float(r["initial_energy"]), 0.0, places=6,
                                   msg=f"truth initial_energy not near zero: {r}")

    def test_truth_final_energy_near_zero(self):
        for r in self.truth_rows:
            self.assertAlmostEqual(float(r["final_energy"]), 0.0, places=6,
                                   msg=f"truth final_energy not near zero: {r}")

    def test_truth_preserved_near_truth(self):
        for r in self.truth_rows:
            self.assertEqual(r["preserved_near_truth"], "true",
                             msg=f"truth row not preserved: {r}")

    def test_truth_initial_distance_zero(self):
        for r in self.truth_rows:
            self.assertAlmostEqual(float(r["initial_distance_to_truth_rms"]), 0.0, places=6,
                                   msg=f"truth initial distance not zero: {r}")


class TestPhase2DNumericSanity(unittest.TestCase):
    def setUp(self):
        self.rows = _load_rows()

    def _is_finite(self, val_str):
        try:
            v = float(val_str)
            return math.isfinite(v)
        except ValueError:
            return val_str == "NA"

    def test_initial_energy_finite_or_na(self):
        for r in self.rows:
            self.assertTrue(self._is_finite(r["initial_energy"]),
                            f"initial_energy not finite: {r['initial_energy']}")

    def test_final_energy_finite_or_na(self):
        for r in self.rows:
            self.assertTrue(self._is_finite(r["final_energy"]),
                            f"final_energy not finite: {r['final_energy']}")

    def test_delta_energy_finite_or_na(self):
        for r in self.rows:
            self.assertTrue(self._is_finite(r["delta_energy"]),
                            f"delta_energy not finite: {r['delta_energy']}")

    def test_boolean_columns_valid(self):
        for r in self.rows:
            for col in ("improved_energy", "improved_interval_rmse", "preserved_near_truth"):
                self.assertIn(r[col], ("true", "false"),
                              f"{col} not boolean: {r[col]}")

    def test_noise_epsilon_na_only_for_random_init(self):
        for r in self.rows:
            if r["init_label"] == "random_init":
                self.assertEqual(r["noise_epsilon"], "NA",
                                 f"random_init row should have noise_epsilon=NA: {r}")
            else:
                val = float(r["noise_epsilon"])
                self.assertGreaterEqual(val, 0.0)


if __name__ == "__main__":
    unittest.main()
