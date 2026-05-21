"""Regression tests for the Phase 2E warmup-skip probe."""

import csv
import math
import os
import unittest

CSV_PATH = os.path.join(
    os.path.dirname(__file__),
    "..",
    "benchmarks",
    "foundation",
    "phase2e_warmup_skip_probe.csv",
)

EXPECTED_COLUMNS = {
    "family",
    "target_dim",
    "n",
    "seed",
    "init_label",
    "warmup_mode",
    "noise_epsilon",
    "paired_key",
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

EXPECTED_INIT_LABELS = {
    "truth",
    "truth_plus_small_noise",
    "truth_plus_medium_noise",
    "random_init",
}


def _load_rows():
    with open(CSV_PATH, newline="") as fh:
        return list(csv.DictReader(fh))


class TestPhase2EOutputExists(unittest.TestCase):
    def test_csv_exists(self):
        self.assertTrue(os.path.exists(CSV_PATH), f"CSV not found: {CSV_PATH}")

    def test_md_exists(self):
        md = CSV_PATH.replace(".csv", ".md")
        self.assertTrue(os.path.exists(md), f"Markdown not found: {md}")


class TestPhase2ESchema(unittest.TestCase):
    def setUp(self):
        self.rows = _load_rows()

    def test_columns_match_expected(self):
        actual = set(self.rows[0].keys())
        self.assertEqual(actual, EXPECTED_COLUMNS)

    def test_row_count(self):
        # 3d × 2n × 3seeds × 4labels × 2modes = 144
        self.assertEqual(len(self.rows), 144)


class TestPhase2EFamilyFilter(unittest.TestCase):
    def setUp(self):
        self.rows = _load_rows()

    def test_only_minkowski(self):
        self.assertEqual({r["family"] for r in self.rows}, {"minkowski"})

    def test_no_kr_or_corona(self):
        for r in self.rows:
            self.assertNotIn("kleitman", r["family"].lower())
            self.assertNotIn("corona", r["family"].lower())


class TestPhase2EWarmupModes(unittest.TestCase):
    def setUp(self):
        self.rows = _load_rows()

    def test_both_warmup_modes_present(self):
        modes = {r["warmup_mode"] for r in self.rows}
        self.assertEqual(modes, {"with_warmup", "skip_warmup"})

    def test_equal_row_count_per_mode(self):
        n_with = sum(1 for r in self.rows if r["warmup_mode"] == "with_warmup")
        n_skip = sum(1 for r in self.rows if r["warmup_mode"] == "skip_warmup")
        self.assertEqual(n_with, n_skip)


class TestPhase2EPairedKey(unittest.TestCase):
    def setUp(self):
        self.rows = _load_rows()

    def test_each_paired_key_has_exactly_two_rows(self):
        from collections import Counter
        counts = Counter(r["paired_key"] for r in self.rows)
        for key, count in counts.items():
            self.assertEqual(count, 2, f"paired_key {key!r} has {count} rows, expected 2")

    def test_paired_key_has_both_modes(self):
        from collections import defaultdict
        by_key: dict[str, set] = defaultdict(set)
        for r in self.rows:
            by_key[r["paired_key"]].add(r["warmup_mode"])
        for key, modes in by_key.items():
            self.assertEqual(modes, {"with_warmup", "skip_warmup"},
                             f"paired_key {key!r} missing a warmup mode")

    def test_paired_key_format(self):
        for r in self.rows:
            parts = r["paired_key"].split("|")
            self.assertEqual(len(parts), 5,
                             f"paired_key should have 5 parts: {r['paired_key']}")


class TestPhase2ETruthInit(unittest.TestCase):
    def setUp(self):
        self.truth_rows = [r for r in _load_rows() if r["init_label"] == "truth"]

    def test_truth_initial_energy_near_zero(self):
        for r in self.truth_rows:
            self.assertAlmostEqual(float(r["initial_energy"]), 0.0, places=6,
                                   msg=f"truth initial_energy not near zero: {r}")

    def test_truth_preserved_in_both_modes(self):
        for r in self.truth_rows:
            self.assertEqual(r["preserved_near_truth"], "true",
                             msg=f"truth row not preserved: {r['warmup_mode']}")


class TestPhase2ENumericSanity(unittest.TestCase):
    def setUp(self):
        self.rows = _load_rows()

    def _is_finite_or_na(self, val_str):
        if val_str == "NA":
            return True
        try:
            return math.isfinite(float(val_str))
        except ValueError:
            return False

    def test_initial_energy_finite_or_na(self):
        for r in self.rows:
            self.assertTrue(self._is_finite_or_na(r["initial_energy"]),
                            f"initial_energy not finite: {r['initial_energy']}")

    def test_final_energy_finite_or_na(self):
        for r in self.rows:
            self.assertTrue(self._is_finite_or_na(r["final_energy"]),
                            f"final_energy not finite: {r['final_energy']}")

    def test_delta_energy_finite_or_na(self):
        for r in self.rows:
            self.assertTrue(self._is_finite_or_na(r["delta_energy"]),
                            f"delta_energy not finite: {r['delta_energy']}")

    def test_boolean_columns_valid(self):
        for r in self.rows:
            for col in ("improved_energy", "improved_interval_rmse", "preserved_near_truth"):
                self.assertIn(r[col], ("true", "false"),
                              f"{col} not boolean: {r[col]}")

    def test_noise_epsilon_na_only_for_random_init(self):
        for r in self.rows:
            if r["init_label"] == "random_init":
                self.assertEqual(r["noise_epsilon"], "NA")
            else:
                self.assertGreaterEqual(float(r["noise_epsilon"]), 0.0)


if __name__ == "__main__":
    unittest.main()
