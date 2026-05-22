"""Regression tests for the Phase 5 seed-level morphology audit."""

from __future__ import annotations

import csv
import unittest
from pathlib import Path

from tools import build_phase4b_survival_probe as p4b
from tools import build_phase5_seed_curve_morphology as p5


ROOT = Path(__file__).resolve().parents[1]
PHASE4B_CSV = ROOT / "benchmarks" / "foundation" / "phase4b_survival_probe.csv"
PHASE4B_PER_SEED_CSV = (
    ROOT / "benchmarks" / "foundation" / "phase4b_survival_probe_per_seed.csv"
)
PHASE5_CSV = ROOT / "benchmarks" / "foundation" / "phase5_seed_curve_morphology.csv"
PHASE5_MD = ROOT / "benchmarks" / "foundation" / "phase5_seed_curve_morphology.md"


def _curve(vals: list[float]) -> list[tuple[float, float, float]]:
    eps = [0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.15, 0.20]
    return [(e, v, v) for e, v in zip(eps, vals)]


class Phase5ClassifierTests(unittest.TestCase):
    def test_synthetic_v_shape_detected(self) -> None:
        shape, censoring, floor, _ = p5.classify_seed_curve(
            _curve([0.30, 0.20, 0.10, 0.03, 0.04, 0.05, 0.06, 0.07])
        )
        self.assertEqual(shape, "seed_v_shape")
        self.assertEqual(censoring, "none")
        self.assertFalse(floor)

    def test_synthetic_monotone_decay_detected(self) -> None:
        shape, censoring, floor, _ = p5.classify_seed_curve(
            _curve([0.30, 0.20, 0.10, 0.05, 0.03, 0.02, 0.01, 0.005])
        )
        self.assertEqual(shape, "seed_monotone_decay")
        self.assertEqual(censoring, "none")
        self.assertFalse(floor)

    def test_synthetic_interior_min_noisy_tail_detected(self) -> None:
        shape, censoring, floor, tail = p5.classify_seed_curve(
            _curve([0.30, 0.12, 0.04, 0.03, 0.05, 0.04, 0.06, 0.05])
        )
        self.assertEqual(shape, "seed_interior_min_noisy_tail")
        self.assertEqual(censoring, "none")
        self.assertFalse(floor)
        self.assertGreater(tail["tail_negative_count"], 0)

    def test_floor_saturated_seed_is_censored(self) -> None:
        shape, censoring, floor, _ = p5.classify_seed_curve(
            _curve([0.20, 0.08, 0.0, 0.01, 0.02]),
            floor_tolerance=1e-6,
        )
        self.assertEqual(shape, "seed_floor_saturated")
        self.assertEqual(censoring, "floor_saturated")
        self.assertTrue(floor)

    def test_insufficient_valid_points(self) -> None:
        shape, censoring, floor, _ = p5.classify_seed_curve(_curve([0.2, 0.1]))
        self.assertEqual(shape, "seed_insufficient_valid_points")
        self.assertEqual(censoring, "insufficient_valid_points")
        self.assertFalse(floor)


class Phase5FixtureTests(unittest.TestCase):
    def test_grouping_reconstructs_ordered_seed_curves_if_generated(self) -> None:
        if not PHASE4B_PER_SEED_CSV.exists():
            self.skipTest("Phase 4B per-seed CSV not generated yet")
        rows = p5.load_phase4b_per_seed(PHASE4B_PER_SEED_CSV)
        grouped = p5.group_seed_curves(rows)
        cell_seed = grouped[(32, 2, 1900)]
        eps = [float(r["epsilon"]) for r in cell_seed]
        self.assertEqual(eps, sorted(eps))
        self.assertEqual(len(eps), 8)

    def test_pilot_grid_contains_nine_cells_if_generated(self) -> None:
        if not PHASE4B_PER_SEED_CSV.exists():
            self.skipTest("Phase 4B per-seed CSV not generated yet")
        rows = p5.load_phase4b_per_seed(PHASE4B_PER_SEED_CSV)
        cells = {(int(r["n"]), int(r["target_dim"])) for r in rows}
        self.assertEqual(len(cells), 9)

    def test_phase4b_outcome_remains_mixed(self) -> None:
        if not PHASE4B_CSV.exists():
            self.skipTest("Phase 4B CSV not generated yet")
        summary = p4b.load_summary_csv(PHASE4B_CSV)
        self.assertEqual(p4b.phase4b_outcome(summary), "MIXED")

    def test_phase5_csv_schema_if_generated(self) -> None:
        if not PHASE5_CSV.exists():
            self.skipTest("Phase 5 CSV not generated yet")
        with PHASE5_CSV.open(newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh)
            header = tuple(next(reader))
        self.assertEqual(header, p5.PHASE5_CSV_HEADERS)

    def test_phase5_markdown_mentions_conservative_scope_if_generated(self) -> None:
        if not PHASE5_MD.exists():
            self.skipTest("Phase 5 Markdown not generated yet")
        text = PHASE5_MD.read_text(encoding="utf-8")
        self.assertIn("No new simulations", text)
        self.assertIn("does not establish a physical law", text)


if __name__ == "__main__":
    unittest.main()
