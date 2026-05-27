"""Tests for S4-SCHW-STABILITY-001: exterior Schwarzschild stability sweep.

Covers two goals:
1. Artifact guard: the generated CSV/JSON/MD files exist, are current, and
   report all_checks_pass=True.
2. In-process functional test: a minimal single-cell sweep (small N) runs
   without errors and returns correct structural output.
"""

from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "explore" / "sorkin4_schwarzschild_benchmark"))

import audit_schw_stability_001 as audit


ARTIFACT_JSON = (
    ROOT
    / "explore"
    / "sorkin4_schwarzschild_benchmark"
    / "schwarzschild_stability_001.json"
)
ARTIFACT_CSV = ARTIFACT_JSON.with_suffix(".csv")
ARTIFACT_MD = ARTIFACT_JSON.with_suffix(".md")


class StabilityArtifactTests(unittest.TestCase):
    """Guard the committed stability-sweep artefacts."""

    def setUp(self) -> None:
        for path in (ARTIFACT_JSON, ARTIFACT_CSV, ARTIFACT_MD):
            self.assertTrue(
                path.exists(),
                msg=(
                    f"Missing S4-SCHW-STABILITY-001 artifact at {path}. "
                    "Run: python explore/sorkin4_schwarzschild_benchmark/"
                    "audit_schw_stability_001.py"
                ),
            )
        with ARTIFACT_JSON.open(encoding="utf-8") as fh:
            self.payload = json.load(fh)

    def test_aggregate_sweep_parameters(self) -> None:
        agg = self.payload["aggregate"]
        self.assertEqual(agg["sweep_N"], [12, 16])
        self.assertEqual(agg["sweep_seed"], [1959, 1960, 1961])
        self.assertAlmostEqual(agg["sweep_margin"][0], 0.25)
        self.assertAlmostEqual(agg["sweep_margin"][1], 0.35)
        self.assertAlmostEqual(agg["sweep_margin"][2], 0.50)
        self.assertEqual(agg["n_cells"], 18)

    def test_all_checks_pass(self) -> None:
        agg = self.payload["aggregate"]
        self.assertTrue(agg["all_checks_pass"])
        self.assertTrue(agg["Q1_order_checks_all_pass"])
        self.assertTrue(agg["Q2_turning_branch_no_invasion"])

    def test_gap_positive(self) -> None:
        agg = self.payload["aggregate"]
        self.assertTrue(agg["Q3_gap_consistent_positive"])
        self.assertGreater(agg["Q3_min_turning_gap_all_cells"], 0.0)

    def test_per_cell_order_checks(self) -> None:
        for row in self.payload["rows"]:
            self.assertTrue(
                row["antisymmetric"],
                msg=f"antisymmetric failed for N={row['N']} seed={row['seed']} margin={row['margin']}",
            )
            self.assertTrue(
                row["transitive"],
                msg=f"transitive failed for N={row['N']} seed={row['seed']} margin={row['margin']}",
            )
            self.assertTrue(
                row["decided_transitivity"],
                msg=f"decided_transitivity failed for N={row['N']} seed={row['seed']} margin={row['margin']}",
            )
            self.assertTrue(
                row["local_checks_pass"],
                msg=f"local_checks_pass failed for N={row['N']} seed={row['seed']} margin={row['margin']}",
            )

    def test_per_cell_turning_gaps_positive(self) -> None:
        for row in self.payload["rows"]:
            self.assertTrue(
                row["outgoing_gaps_all_positive"],
                msg=(
                    f"Turning-branch gap ≤ 0 for N={row['N']} seed={row['seed']} "
                    f"margin={row['margin']}: min_gap={row['min_gap']}"
                ),
            )


class StabilityFunctionalTests(unittest.TestCase):
    """In-process functional tests: a tiny single-cell run exercises the code paths."""

    def test_single_cell_structure(self) -> None:
        row = audit.run_stability_cell(n=8, seed=1959, margin=0.35)
        for field in audit.CSV_FIELDS:
            self.assertIn(field, row, msg=f"field {field!r} missing from run_stability_cell output")

    def test_single_cell_order_checks(self) -> None:
        row = audit.run_stability_cell(n=8, seed=1959, margin=0.35)
        self.assertTrue(row["antisymmetric"])
        self.assertTrue(row["transitive"])
        self.assertTrue(row["decided_transitivity"])

    def test_single_cell_gap_positive_or_none(self) -> None:
        row = audit.run_stability_cell(n=8, seed=1959, margin=0.35)
        if row["min_gap"] is not None:
            self.assertGreater(row["min_gap"], 0.0)

    def test_event_generation_with_margin(self) -> None:
        """Events must stay outside r_s + margin."""
        for margin in (0.25, 0.35, 0.50):
            r_min = 2.0 + margin
            events = audit.generate_exterior_events_with_margin(12, 1959, margin)
            self.assertEqual(len(events), 12)
            for event in events:
                self.assertGreaterEqual(
                    event.r,
                    r_min - 1.0e-12,
                    msg=f"Event r={event.r} below r_min={r_min} for margin={margin}",
                )
                self.assertLessEqual(event.r, audit.R_MAX + 1.0e-12)

    def test_mini_sweep_aggregate(self) -> None:
        """A 2×2×2 mini-sweep should return 8 cells and aggregate cleanly."""
        rows = audit.run_full_sweep(
            n_values=(8, 10),
            seeds=(1959, 1960),
            margins=(0.35, 0.50),
        )
        self.assertEqual(len(rows), 8)
        agg = audit.compute_aggregate(rows)
        self.assertIn("all_checks_pass", agg)
        self.assertTrue(agg["Q1_order_checks_all_pass"])


if __name__ == "__main__":
    unittest.main()
