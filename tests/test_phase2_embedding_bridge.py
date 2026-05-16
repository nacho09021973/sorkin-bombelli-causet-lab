"""Integrity tests for the Phase 2 embedding bridge artifact."""

from __future__ import annotations

import csv
import math
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BRIDGE_CSV = ROOT / "benchmarks" / "foundation" / "phase2_embedding_bridge.csv"
BRIDGE_MD = ROOT / "benchmarks" / "foundation" / "phase2_embedding_bridge.md"

EXPECTED_HEADERS = (
    "family",
    "target_dim",
    "n",
    "seed",
    "optimizer_seed",
    "embedding_dim",
    "mm_dim",
    "midpoint_dim",
    "abs_discrepancy_mm_midpoint",
    "chain2_count",
    "chain3_count",
    "chain3_abundance",
    "initial_energy",
    "warmup_energy",
    "final_energy",
    "truth_energy",
    "energy_gap",
    "interval_rmse",
    "optimizer_status",
    "failure_mode",
    "runtime_seconds",
    "optimizer_steps",
)

STRUCTURAL_NUMERIC_COLUMNS = (
    "n",
    "seed",
    "optimizer_seed",
    "embedding_dim",
    "mm_dim",
    "midpoint_dim",
    "abs_discrepancy_mm_midpoint",
    "chain2_count",
    "chain3_count",
    "chain3_abundance",
)

POST_EMBEDDING_COLUMNS = (
    "initial_energy",
    "warmup_energy",
    "final_energy",
    "truth_energy",
    "energy_gap",
    "interval_rmse",
    "optimizer_status",
    "failure_mode",
    "runtime_seconds",
    "optimizer_steps",
)


class Phase2EmbeddingBridgeIntegrityTests(unittest.TestCase):
    def setUp(self) -> None:
        self.assertTrue(
            BRIDGE_CSV.exists(),
            msg=(
                f"missing Phase 2 CSV at {BRIDGE_CSV}; "
                "run `make regen-phase2`"
            ),
        )
        with BRIDGE_CSV.open(newline="", encoding="utf-8") as fh:
            self.rows = list(csv.DictReader(fh))

    def test_markdown_exists(self) -> None:
        self.assertTrue(
            BRIDGE_MD.exists(),
            msg=(
                f"missing Phase 2 markdown at {BRIDGE_MD}; "
                "run `make regen-phase2`"
            ),
        )

    def test_columns_match_expected(self) -> None:
        with BRIDGE_CSV.open(newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh)
            header = tuple(next(reader))
        self.assertEqual(header, EXPECTED_HEADERS)

    def test_has_minkowski_and_control_rows(self) -> None:
        families = {row["family"] for row in self.rows}
        self.assertIn("minkowski", families)
        self.assertTrue(
            {"kleitman_rothschild", "corona_poset"} & families,
            msg=f"expected at least one negative control, got {families}",
        )

    def test_structural_columns_are_finite(self) -> None:
        for row in self.rows:
            for col in STRUCTURAL_NUMERIC_COLUMNS:
                raw = row[col]
                self.assertNotIn(raw.lower(), {"nan", "inf", "-inf"})
                value = float(raw)
                self.assertTrue(math.isfinite(value), msg=f"{col}={raw}")

    def test_post_embedding_columns_exist(self) -> None:
        for row in self.rows:
            for col in POST_EMBEDDING_COLUMNS:
                self.assertIn(col, row)

    def test_completed_rows_have_finite_core_embedding_metrics(self) -> None:
        for row in self.rows:
            if row["optimizer_status"] != "completed":
                continue
            for col in ("initial_energy", "warmup_energy", "final_energy"):
                value = float(row[col])
                self.assertTrue(math.isfinite(value), msg=f"{col}={value}")
            self.assertGreater(int(row["optimizer_steps"]), 0)


if __name__ == "__main__":
    unittest.main()
