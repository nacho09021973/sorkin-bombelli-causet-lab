"""Regression tests for the Phase 2F guarded-warmup probe."""

import csv
import math
import os
import unittest

CSV_PATH = os.path.join(
    os.path.dirname(__file__),
    "..",
    "benchmarks",
    "foundation",
    "phase2f_guarded_warmup_probe.csv",
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
    "warmup_attempted_moves",
    "warmup_accepted_moves",
    "warmup_rejected_moves",
    "warmup_energy_before",
    "warmup_energy_after",
    "warmup_delta_energy",
    "notes",
}

EXPECTED_WARMUP_MODES = {"legacy_warmup", "skip_warmup", "guarded_warmup"}
EXPECTED_INIT_LABELS = {
    "truth",
    "truth_plus_small_noise",
    "truth_plus_medium_noise",
    "random_init",
}

# Tolerance for warmup energy guard test (normalization effects).
NORMALIZATION_TOLERANCE = 0.5


def _load_rows():
    with open(CSV_PATH, newline="") as fh:
        return list(csv.DictReader(fh))


class TestPhase2FOutputExists(unittest.TestCase):
    def test_csv_exists(self):
        self.assertTrue(os.path.exists(CSV_PATH), f"CSV not found: {CSV_PATH}")

    def test_md_exists(self):
        md = CSV_PATH.replace(".csv", ".md")
        self.assertTrue(os.path.exists(md), f"Markdown not found: {md}")


class TestPhase2FSchema(unittest.TestCase):
    def setUp(self):
        self.rows = _load_rows()

    def test_columns_match_expected(self):
        self.assertEqual(set(self.rows[0].keys()), EXPECTED_COLUMNS)

    def test_row_count(self):
        # 3d × 2n × 3seeds × 4labels × 3modes = 216
        self.assertEqual(len(self.rows), 216)


class TestPhase2FFamilyFilter(unittest.TestCase):
    def setUp(self):
        self.rows = _load_rows()

    def test_only_minkowski(self):
        self.assertEqual({r["family"] for r in self.rows}, {"minkowski"})

    def test_no_kr_or_corona(self):
        for r in self.rows:
            self.assertNotIn("kleitman", r["family"].lower())
            self.assertNotIn("corona", r["family"].lower())


class TestPhase2FWarmupModes(unittest.TestCase):
    def setUp(self):
        self.rows = _load_rows()

    def test_all_warmup_modes_present(self):
        modes = {r["warmup_mode"] for r in self.rows}
        self.assertEqual(modes, EXPECTED_WARMUP_MODES)

    def test_equal_count_per_mode(self):
        for mode in EXPECTED_WARMUP_MODES:
            count = sum(1 for r in self.rows if r["warmup_mode"] == mode)
            self.assertEqual(count, 72, f"Expected 72 rows for {mode}")


class TestPhase2FPairedKey(unittest.TestCase):
    def setUp(self):
        self.rows = _load_rows()

    def test_each_paired_key_has_exactly_three_rows(self):
        from collections import Counter
        counts = Counter(r["paired_key"] for r in self.rows)
        for key, count in counts.items():
            self.assertEqual(count, 3, f"paired_key {key!r} has {count} rows, expected 3")

    def test_paired_key_has_all_three_modes(self):
        from collections import defaultdict
        by_key: dict[str, set] = defaultdict(set)
        for r in self.rows:
            by_key[r["paired_key"]].add(r["warmup_mode"])
        for key, modes in by_key.items():
            self.assertEqual(modes, EXPECTED_WARMUP_MODES,
                             f"paired_key {key!r} missing warmup modes: {modes}")


class TestPhase2FTruthInit(unittest.TestCase):
    def setUp(self):
        self.truth_rows = [r for r in _load_rows() if r["init_label"] == "truth"]

    def test_truth_initial_energy_near_zero(self):
        for r in self.truth_rows:
            self.assertAlmostEqual(float(r["initial_energy"]), 0.0, places=6,
                                   msg=f"truth initial_energy not near zero: {r}")

    def test_truth_preserved_in_all_modes(self):
        for r in self.truth_rows:
            self.assertEqual(r["preserved_near_truth"], "true",
                             msg=f"truth not preserved in mode {r['warmup_mode']}")


class TestPhase2FNumericSanity(unittest.TestCase):
    def setUp(self):
        self.rows = _load_rows()

    def _is_finite_or_na(self, val_str):
        if val_str == "NA":
            return True
        try:
            return math.isfinite(float(val_str))
        except ValueError:
            return False

    def test_energies_finite(self):
        for r in self.rows:
            for col in ("initial_energy", "final_energy", "delta_energy"):
                self.assertTrue(self._is_finite_or_na(r[col]),
                                f"{col} not finite: {r[col]}")

    def test_interval_rmse_finite(self):
        for r in self.rows:
            for col in ("initial_interval_rmse", "final_interval_rmse"):
                self.assertTrue(self._is_finite_or_na(r[col]),
                                f"{col} not finite: {r[col]}")

    def test_distances_finite(self):
        for r in self.rows:
            for col in ("initial_distance_to_truth_rms", "final_distance_to_truth_rms"):
                self.assertTrue(self._is_finite_or_na(r[col]),
                                f"{col} not finite: {r[col]}")

    def test_boolean_columns_valid(self):
        for r in self.rows:
            for col in ("improved_energy", "improved_interval_rmse", "preserved_near_truth"):
                self.assertIn(r[col], ("true", "false"),
                              f"{col} not boolean: {r[col]}")

    def test_warmup_move_counts_non_negative(self):
        for r in self.rows:
            self.assertGreaterEqual(int(r["warmup_attempted_moves"]), 0)
            self.assertGreaterEqual(int(r["warmup_accepted_moves"]), 0)
            self.assertGreaterEqual(int(r["warmup_rejected_moves"]), 0)

    def test_warmup_move_counts_consistent(self):
        for r in self.rows:
            att = int(r["warmup_attempted_moves"])
            acc = int(r["warmup_accepted_moves"])
            rej = int(r["warmup_rejected_moves"])
            self.assertEqual(att, acc + rej,
                             f"attempted != accepted + rejected: {att} != {acc} + {rej}")


class TestPhase2FGuardedWarmupProperty(unittest.TestCase):
    """Test 8: guarded_warmup should not worsen warmup energy beyond tolerance."""

    def setUp(self):
        self.guarded_rows = [
            r for r in _load_rows() if r["warmup_mode"] == "guarded_warmup"
        ]

    def test_guarded_warmup_energy_not_worsened(self):
        for r in self.guarded_rows:
            e_before = float(r["warmup_energy_before"])
            e_after = float(r["warmup_energy_after"])
            # Guard is applied pre-normalization; allow normalization tolerance.
            tolerance = max(NORMALIZATION_TOLERANCE * e_before, 1e-6)
            self.assertLessEqual(
                e_after,
                e_before + tolerance,
                msg=(
                    f"guarded_warmup increased energy beyond tolerance "
                    f"({e_after:.6g} > {e_before:.6g} + {tolerance:.6g}) "
                    f"for {r['init_label']} d={r['target_dim']} n={r['n']} seed={r['seed']}"
                ),
            )

    def test_skip_warmup_has_zero_moves(self):
        skip_rows = [r for r in _load_rows() if r["warmup_mode"] == "skip_warmup"]
        for r in skip_rows:
            self.assertEqual(int(r["warmup_attempted_moves"]), 0)
            self.assertEqual(int(r["warmup_accepted_moves"]), 0)
            self.assertEqual(int(r["warmup_rejected_moves"]), 0)

    def test_legacy_warmup_has_zero_rejected(self):
        legacy_rows = [r for r in _load_rows() if r["warmup_mode"] == "legacy_warmup"]
        for r in legacy_rows:
            self.assertEqual(int(r["warmup_rejected_moves"]), 0,
                             f"legacy_warmup should have 0 rejected: {r}")


if __name__ == "__main__":
    unittest.main()
