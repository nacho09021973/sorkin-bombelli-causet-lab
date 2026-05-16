"""Integrity tests for the Phase 2B annealer schedule probe artifact."""

from __future__ import annotations

import csv
import math
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROBE_CSV = ROOT / "benchmarks" / "foundation" / "phase2b_annealer_schedule_probe.csv"
PROBE_MD = ROOT / "benchmarks" / "foundation" / "phase2b_annealer_schedule_probe.md"

EXPECTED_HEADERS = (
    "family",
    "target_dim",
    "n",
    "seed",
    "schedule_label",
    "warmup_limit",
    "anneal_limit",
    "max_data",
    "initial_temp",
    "cooling_factor",
    "optimizer_seed",
    "embedding_dim",
    "initial_energy",
    "warmup_energy",
    "truth_energy",
    "final_energy",
    "energy_gap",
    "interval_rmse",
    "success_flag",
    "runtime_seconds",
)

# An exact Minkowski sprinkling has Bombelli energy zero at its true
# coordinates; we still allow a tiny floating-point slack to keep the
# test from being fragile under harmless numerical noise.
TRUTH_ENERGY_TOL = 1e-9


class Phase2BAnnealerScheduleProbeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.assertTrue(
            PROBE_CSV.exists(),
            msg=(
                f"missing Phase 2B CSV at {PROBE_CSV}; "
                "run `make regen-phase2b`"
            ),
        )
        with PROBE_CSV.open(newline="", encoding="utf-8") as fh:
            self.rows = list(csv.DictReader(fh))

    def test_markdown_exists(self) -> None:
        self.assertTrue(
            PROBE_MD.exists(),
            msg=(
                f"missing Phase 2B markdown at {PROBE_MD}; "
                "run `make regen-phase2b`"
            ),
        )

    def test_columns_match_expected(self) -> None:
        with PROBE_CSV.open(newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh)
            header = tuple(next(reader))
        self.assertEqual(header, EXPECTED_HEADERS)

    def test_has_at_least_one_row(self) -> None:
        self.assertGreater(len(self.rows), 0)

    def test_only_minkowski_family(self) -> None:
        families = {row["family"] for row in self.rows}
        self.assertEqual(
            families,
            {"minkowski"},
            msg=(
                "Phase 2B must not include KR or corona controls; "
                f"got families={families}"
            ),
        )

    def test_truth_energy_compatible_with_zero(self) -> None:
        for row in self.rows:
            raw = row["truth_energy"]
            self.assertNotIn(raw.lower(), {"nan", "inf", "-inf", "na"})
            value = float(raw)
            self.assertTrue(math.isfinite(value), msg=f"truth_energy={raw}")
            self.assertLessEqual(
                abs(value),
                TRUTH_ENERGY_TOL,
                msg=(
                    f"truth_energy={value} should be numerically zero "
                    f"for a Minkowski sprinkling (row={row})"
                ),
            )

    def test_energy_gap_and_rmse_are_finite(self) -> None:
        for row in self.rows:
            for col in ("energy_gap", "interval_rmse"):
                raw = row[col]
                self.assertNotIn(raw.lower(), {"nan", "inf", "-inf", "na"})
                value = float(raw)
                self.assertTrue(math.isfinite(value), msg=f"{col}={raw}")

    def test_schedule_labels_cover_short_medium_long(self) -> None:
        labels = {row["schedule_label"] for row in self.rows}
        self.assertEqual(labels, {"short", "medium", "long"})

    def test_target_dims_cover_2_3_4(self) -> None:
        dims = {int(row["target_dim"]) for row in self.rows}
        self.assertEqual(dims, {2, 3, 4})

    def test_success_flag_is_boolean(self) -> None:
        for row in self.rows:
            self.assertIn(row["success_flag"], {"true", "false"})

    def test_embedding_dim_is_d_minus_one(self) -> None:
        for row in self.rows:
            d = int(row["target_dim"])
            edim = int(row["embedding_dim"])
            self.assertEqual(edim, d - 1)


if __name__ == "__main__":
    unittest.main()
