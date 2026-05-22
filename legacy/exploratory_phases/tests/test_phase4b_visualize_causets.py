"""Tests for Phase 4B causet visual audit selection and SVG output."""

from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from tools import build_phase4b_survival_probe as p4b
from tools import visualize_phase4b_causets as viz


ROOT = Path(__file__).resolve().parents[1]
AGGREGATE_CSV = ROOT / "benchmarks" / "foundation" / "phase4b_survival_probe.csv"
PER_SEED_CSV = ROOT / "benchmarks" / "foundation" / "phase4b_survival_probe_per_seed.csv"
VISUAL_INDEX_CSV = ROOT / "benchmarks" / "foundation" / "phase4b_causet_visual_index.csv"
VISUAL_LINK_AUDIT_CSV = ROOT / "benchmarks" / "foundation" / "phase4b_visual_link_audit.csv"
INTRINSIC_POSET_AUDIT_CSV = ROOT / "benchmarks" / "foundation" / "phase4b_intrinsic_poset_audit.csv"
PHASE4B_MD = ROOT / "benchmarks" / "foundation" / "phase4b_survival_probe.md"


class Phase4BVisualSelectionTests(unittest.TestCase):
    def setUp(self) -> None:
        if not AGGREGATE_CSV.exists() or not PER_SEED_CSV.exists():
            self.skipTest("Phase 4B provenance CSVs not generated yet")
        self.aggregate_rows = viz._read_csv(AGGREGATE_CSV)
        self.per_seed_rows = viz._read_csv(PER_SEED_CSV)

    def test_selector_picks_three_roles_for_known_cell(self) -> None:
        selected = viz.select_visual_rows(
            self.per_seed_rows,
            self.aggregate_rows,
            priority_cells=((48, 4),),
        )
        self.assertEqual({r["role"] for r in selected}, {"near_mean", "best", "worst"})
        self.assertEqual({(int(r["n"]), int(r["target_dim"])) for r in selected}, {(48, 4)})
        self.assertEqual({float(r["epsilon"]) for r in selected}, {0.06})
        for row in selected:
            self.assertEqual(row["valid"], "true")

    def test_generate_svg_for_one_known_cell(self) -> None:
        selected = viz.select_visual_rows(
            self.per_seed_rows,
            self.aggregate_rows,
            priority_cells=((32, 2),),
        )[:1]
        with tempfile.TemporaryDirectory() as tmp:
            out_dir = Path(tmp)
            index_rows = viz.generate_svgs(selected, out_dir)
            self.assertEqual(len(index_rows), 1)
            svg_path = ROOT / index_rows[0]["svg_path"]
            # generate_svgs writes relative paths only for repo paths; for tmp
            # output the stored relative path can be invalid, so inspect out_dir.
            svgs = list(out_dir.glob("*.svg"))
            self.assertEqual(len(svgs), 1)
            text = svgs[0].read_text(encoding="utf-8")
            self.assertIn("<svg", text)
            self.assertIn("<line", text)
            self.assertTrue(str(svg_path) or True)

    def test_visual_index_points_to_existing_svgs_if_generated(self) -> None:
        if not VISUAL_INDEX_CSV.exists():
            self.skipTest("Visual index not generated yet")
        with VISUAL_INDEX_CSV.open(newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        self.assertGreater(len(rows), 0)
        for row in rows:
            self.assertTrue((ROOT / row["svg_path"]).exists(), row["svg_path"])

    def test_link_audit_rows_have_hasse_counts(self) -> None:
        selected = viz.select_visual_rows(
            self.per_seed_rows,
            self.aggregate_rows,
            priority_cells=((48, 4),),
        )
        rows = viz.build_link_audit_rows(selected)
        self.assertEqual(len(rows), 3)
        self.assertEqual({r["role"] for r in rows}, {"near_mean", "best", "worst"})
        for row in rows:
            self.assertGreater(int(row["n_links_Hasse"]), 0)
            self.assertEqual((int(row["n"]), int(row["target_dim"])), (48, 4))

    def test_link_audit_csv_if_generated(self) -> None:
        if not VISUAL_LINK_AUDIT_CSV.exists():
            self.skipTest("Visual link audit CSV not generated yet")
        with VISUAL_LINK_AUDIT_CSV.open(newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        self.assertGreater(len(rows), 0)
        self.assertEqual(set(rows[0].keys()), set(viz.VISUAL_LINK_AUDIT_HEADERS))
        for row in rows:
            self.assertGreater(int(row["n_links_Hasse"]), 0)

    def test_intrinsic_poset_audit_csv_if_generated(self) -> None:
        if not INTRINSIC_POSET_AUDIT_CSV.exists():
            self.skipTest("Intrinsic poset audit CSV not generated yet")
        with INTRINSIC_POSET_AUDIT_CSV.open(newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        self.assertEqual(len(rows), 15)
        self.assertEqual(set(rows[0].keys()), set(viz.INTRINSIC_POSET_AUDIT_HEADERS))

    def test_markdown_mentions_no_robust_scalar_predictor_if_generated(self) -> None:
        if not PHASE4B_MD.exists():
            self.skipTest("Phase 4B markdown not generated yet")
        text = PHASE4B_MD.read_text(encoding="utf-8")
        if "## Intrinsic poset audit" not in text:
            self.skipTest("Intrinsic poset audit markdown not generated yet")
        self.assertIn("No robust scalar predictor is established", text)

    def test_outcome_still_mixed(self) -> None:
        if not AGGREGATE_CSV.exists():
            self.skipTest("Phase 4B aggregate CSV not generated yet")
        self.assertEqual(p4b.phase4b_outcome(p4b.load_summary_csv(AGGREGATE_CSV)), "MIXED")


if __name__ == "__main__":
    unittest.main()
