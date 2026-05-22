"""Regression tests for the Phase 4B exploratory survival probe."""

from __future__ import annotations

import csv
import unittest
from pathlib import Path

from tools import build_phase4a_epsilon_sweep as p4a
from tools import build_phase4b_survival_probe as p4b


ROOT = Path(__file__).resolve().parents[1]
PHASE4A_CSV = ROOT / "benchmarks" / "foundation" / "phase4a_epsilon_sweep.csv"
PHASE4B_CSV = ROOT / "benchmarks" / "foundation" / "phase4b_survival_probe.csv"
PHASE4B_PER_EPSILON_CSV = (
    ROOT / "benchmarks" / "foundation" / "phase4b_survival_probe_per_epsilon.csv"
)
PHASE4B_PER_SEED_CSV = (
    ROOT / "benchmarks" / "foundation" / "phase4b_survival_probe_per_seed.csv"
)
PHASE4B_MD = ROOT / "benchmarks" / "foundation" / "phase4b_survival_probe.md"

EXPECTED_HEADERS = p4b.CSV_HEADERS
EXPECTED_PER_EPSILON_HEADERS = p4b.PER_EPSILON_CSV_HEADERS
EXPECTED_PER_SEED_HEADERS = p4b.PER_SEED_CSV_HEADERS


def _curve(vals: list[float]) -> list[tuple[float, float, float]]:
    eps = list(p4a.EPSILONS[:len(vals)])
    return [(e, v, v) for e, v in zip(eps, vals)]


class Phase4BClassifierTests(unittest.TestCase):
    def test_synthetic_v_shape_detected(self) -> None:
        curve = _curve([0.30, 0.20, 0.10, 0.03, 0.04, 0.05, 0.06, 0.07])
        self.assertEqual(p4a._classify_curve_shape(curve, idx=1), "v_shape")

    def test_synthetic_monotone_decay_is_not_v_shape(self) -> None:
        curve = _curve([0.30, 0.20, 0.10, 0.05, 0.03, 0.02, 0.01, 0.005])
        self.assertEqual(p4a._classify_curve_shape(curve, idx=1), "monotone_decay")

    def test_floor_saturation_censors_high_discrepancy_non_v(self) -> None:
        label, high_disc, floor_sat = p4b.classify_survival_cell(
            target_dim=3,
            curve_shape="monotone_decay",
            dim_discrepancy_rel_midpoint=0.50,
            min_val=0.0,
        )
        self.assertTrue(high_disc)
        self.assertTrue(floor_sat)
        self.assertEqual(label, "censored_floor")

    def test_non_floor_high_discrepancy_non_v_is_counterexample(self) -> None:
        label, high_disc, floor_sat = p4b.classify_survival_cell(
            target_dim=3,
            curve_shape="monotone_decay",
            dim_discrepancy_rel_midpoint=0.50,
            min_val=0.02,
        )
        self.assertTrue(high_disc)
        self.assertFalse(floor_sat)
        self.assertEqual(label, "counterexample")


class Phase4BPhase4ACompatibilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.assert_path_exists(PHASE4A_CSV)
        rows = p4a.load_rows_from_csv(PHASE4A_CSV)
        cls.summary = {
            (r["n"], r["target_dim"]): r
            for r in p4b.summarize_curves(rows, grid="phase4a_reference")
        }

    @staticmethod
    def assert_path_exists(path: Path) -> None:
        if not path.exists():
            raise AssertionError(f"missing fixture: {path}")

    def test_phase4a_v_shape_cells_remain_v_shape(self) -> None:
        for cell in ((32, 3), (32, 4), (64, 4)):
            self.assertEqual(
                self.summary[cell]["curve_shape"],
                "v_shape",
                msg=f"Phase 4A reference cell {cell} lost v_shape classification",
            )

    def test_phase4a_negative_controls_do_not_become_false_positive_v_shapes(self) -> None:
        for cell in ((32, 2), (64, 2)):
            self.assertNotEqual(
                self.summary[cell]["curve_shape"],
                "v_shape",
                msg=f"Phase 4A d=2 control {cell} became a V-shape false positive",
            )

    def test_floor_case_is_not_strong_negative_when_at_floor(self) -> None:
        row = self.summary[(64, 3)]
        if row["min_val"] <= row["floor_tolerance"]:
            self.assertNotEqual(row["survival_label"], "counterexample")


class Phase4BOutputFixtureTests(unittest.TestCase):
    def test_csv_schema_if_generated(self) -> None:
        if not PHASE4B_CSV.exists():
            self.skipTest("Phase 4B CSV not generated yet; run make regen-phase4b")
        with PHASE4B_CSV.open(newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh)
            header = tuple(next(reader))
        if "borderline_v_like" not in header:
            self.skipTest("Phase 4B CSV has not been regenerated with tail audit columns")
        self.assertEqual(header, EXPECTED_HEADERS)

    def test_markdown_contains_required_sections_if_generated(self) -> None:
        if not PHASE4B_MD.exists():
            self.skipTest("Phase 4B markdown not generated yet; run make regen-phase4b")
        text = PHASE4B_MD.read_text(encoding="utf-8")
        for section in (
            "## Objective",
            "## Grid design",
            "## Curve-shape summary",
            "## Survival test of Phase 4A hypothesis",
            "## Negative controls target_dim=2",
            "## Floor-saturated / censored cases",
            "## Global exploratory outcome",
            "## Conservative conclusion",
        ):
            self.assertIn(section, text)
        if PHASE4B_PER_EPSILON_CSV.exists():
            self.assertIn("## Provenance note", text)
        if PHASE4B_CSV.exists():
            with PHASE4B_CSV.open(newline="", encoding="utf-8") as fh:
                header = next(csv.reader(fh))
            if "borderline_v_like" in header:
                self.assertIn("## Tail-cleanliness / borderline audit", text)
        if PHASE4B_PER_SEED_CSV.exists():
            self.assertIn("phase4b_survival_probe_per_seed.csv", text)

    def test_per_epsilon_csv_schema_if_generated(self) -> None:
        if not PHASE4B_PER_EPSILON_CSV.exists():
            self.skipTest(
                "Phase 4B per-epsilon CSV not generated yet; run make regen-phase4b"
            )
        with PHASE4B_PER_EPSILON_CSV.open(newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh)
            header = tuple(next(reader))
        if "borderline_v_like_cell" not in header:
            self.skipTest(
                "Phase 4B per-epsilon CSV has not been regenerated with tail audit columns"
            )
        self.assertEqual(header, EXPECTED_PER_EPSILON_HEADERS)

    def test_per_epsilon_csv_contains_counterexample_cells_if_generated(self) -> None:
        if not PHASE4B_PER_EPSILON_CSV.exists():
            self.skipTest(
                "Phase 4B per-epsilon CSV not generated yet; run make regen-phase4b"
            )
        with PHASE4B_PER_EPSILON_CSV.open(newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        cells = {(int(r["n"]), int(r["target_dim"])) for r in rows}
        self.assertIn((48, 3), cells)
        self.assertIn((48, 4), cells)

    def test_first_delta_is_na_per_cell_if_generated(self) -> None:
        if not PHASE4B_PER_EPSILON_CSV.exists():
            self.skipTest(
                "Phase 4B per-epsilon CSV not generated yet; run make regen-phase4b"
            )
        with PHASE4B_PER_EPSILON_CSV.open(newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        by_cell: dict[tuple[int, int], list[dict]] = {}
        for row in rows:
            by_cell.setdefault((int(row["n"]), int(row["target_dim"])), []).append(row)
        for cell_rows in by_cell.values():
            first = sorted(cell_rows, key=lambda r: float(r["epsilon"]))[0]
            self.assertEqual(first["delta_from_prev_epsilon"], "NA")
            self.assertEqual(first["delta_sign"], "NA")

    def test_is_min_epsilon_matches_aggregate_if_generated(self) -> None:
        if not PHASE4B_CSV.exists() or not PHASE4B_PER_EPSILON_CSV.exists():
            self.skipTest("Phase 4B CSVs not generated yet; run make regen-phase4b")
        with PHASE4B_CSV.open(newline="", encoding="utf-8") as fh:
            aggregate_rows = list(csv.DictReader(fh))
        with PHASE4B_PER_EPSILON_CSV.open(newline="", encoding="utf-8") as fh:
            per_eps_rows = list(csv.DictReader(fh))
        expected = {
            (int(r["n"]), int(r["target_dim"])): float(r["epsilon_at_min"])
            for r in aggregate_rows
        }
        by_cell: dict[tuple[int, int], list[dict]] = {}
        for row in per_eps_rows:
            by_cell.setdefault((int(row["n"]), int(row["target_dim"])), []).append(row)
        for cell, cell_rows in by_cell.items():
            marked = [r for r in cell_rows if r["is_min_epsilon"] == "true"]
            self.assertEqual(len(marked), 1, msg=f"{cell} has {len(marked)} minima")
            self.assertAlmostEqual(float(marked[0]["epsilon"]), expected[cell])

    def test_pilot_outcome_remains_mixed_if_generated(self) -> None:
        if not PHASE4B_CSV.exists():
            self.skipTest("Phase 4B CSV not generated yet; run make regen-phase4b")
        summary = p4b.load_summary_csv(PHASE4B_CSV)
        self.assertEqual(p4b.phase4b_outcome(summary), "MIXED")

    def test_per_seed_csv_schema_if_generated(self) -> None:
        if not PHASE4B_PER_SEED_CSV.exists():
            self.skipTest(
                "Phase 4B per-seed CSV not generated yet; run make regen-phase4b"
            )
        with PHASE4B_PER_SEED_CSV.open(newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh)
            header = tuple(next(reader))
        self.assertEqual(header, EXPECTED_PER_SEED_HEADERS)

    def test_per_seed_csv_contains_priority_cells_if_generated(self) -> None:
        if not PHASE4B_PER_SEED_CSV.exists():
            self.skipTest(
                "Phase 4B per-seed CSV not generated yet; run make regen-phase4b"
            )
        with PHASE4B_PER_SEED_CSV.open(newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        cells = {(int(r["n"]), int(r["target_dim"])) for r in rows}
        for cell in ((32, 4), (48, 3), (48, 4), (64, 4)):
            self.assertIn(cell, cells)

    def test_per_seed_csv_covers_pilot_epsilons_and_seeds_if_generated(self) -> None:
        if not PHASE4B_PER_SEED_CSV.exists():
            self.skipTest(
                "Phase 4B per-seed CSV not generated yet; run make regen-phase4b"
            )
        with PHASE4B_PER_SEED_CSV.open(newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        expected_eps = set(p4a.EPSILONS)
        expected_seeds = set(p4a.PHASE4A_SEEDS)
        for cell in ((32, 4), (48, 3), (48, 4), (64, 4)):
            cell_rows = [
                r for r in rows
                if (int(r["n"]), int(r["target_dim"])) == cell
            ]
            self.assertEqual({float(r["epsilon"]) for r in cell_rows}, expected_eps)
            self.assertEqual({int(r["seed"]) for r in cell_rows}, expected_seeds)

    def test_aggregate_and_per_epsilon_csvs_exist_if_per_seed_generated(self) -> None:
        if not PHASE4B_PER_SEED_CSV.exists():
            self.skipTest(
                "Phase 4B per-seed CSV not generated yet; run make regen-phase4b"
            )
        self.assertTrue(PHASE4B_CSV.exists())
        self.assertTrue(PHASE4B_PER_EPSILON_CSV.exists())

    def test_tail_audit_marks_phase4b_counterexamples_without_reclassification(self) -> None:
        if not PHASE4B_CSV.exists():
            self.skipTest("Phase 4B CSV not generated yet; run make regen-phase4b")
        with PHASE4B_CSV.open(newline="", encoding="utf-8") as fh:
            header = next(csv.reader(fh))
        if "borderline_v_like" not in header:
            self.skipTest("Phase 4B CSV has no tail audit columns yet")
        rows = {
            (r["n"], r["target_dim"]): r
            for r in p4b.load_summary_csv(PHASE4B_CSV)
        }

        row_483 = rows[(48, 3)]
        self.assertEqual(row_483["curve_shape"], "v_shape")
        self.assertEqual(row_483["survival_label"], "counterexample")
        self.assertTrue(row_483["borderline_v_like"])
        self.assertGreaterEqual(row_483["rise_frac_margin"], 0.0)
        self.assertLessEqual(row_483["rise_frac_margin"], p4a.V_RISE_FRAC)

        row_484 = rows[(48, 4)]
        self.assertEqual(row_484["curve_shape"], "monotone_decay")
        self.assertTrue(row_484["has_interior_minimum"])
        self.assertEqual(row_484["tail_pattern"], "positive,negative,positive,negative")
        self.assertTrue(row_484["borderline_v_like"])

    def test_d4_supporting_cells_keep_clean_positive_tails(self) -> None:
        if not PHASE4B_CSV.exists():
            self.skipTest("Phase 4B CSV not generated yet; run make regen-phase4b")
        with PHASE4B_CSV.open(newline="", encoding="utf-8") as fh:
            header = next(csv.reader(fh))
        if "tail_pattern" not in header:
            self.skipTest("Phase 4B CSV has no tail audit columns yet")
        rows = {
            (r["n"], r["target_dim"]): r
            for r in p4b.load_summary_csv(PHASE4B_CSV)
        }
        self.assertEqual(rows[(32, 4)]["tail_pattern"], "positive,positive")
        self.assertEqual(
            rows[(64, 4)]["tail_pattern"],
            "positive,positive,positive,positive",
        )


if __name__ == "__main__":
    unittest.main()
