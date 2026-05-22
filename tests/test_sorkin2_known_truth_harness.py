"""Tests for the minimal SORKIN-2 known-truth run harness.

These tests exercise result directory construction, JSON/report writing, and
case gating without running the Bombelli annealer.
"""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tools import run_sorkin2_known_truth_case as harness


def make_chain(n: int) -> list[list[bool]]:
    return [[i < j for j in range(n)] for i in range(n)]


def make_antichain(n: int) -> list[list[bool]]:
    return [[False] * n for _ in range(n)]


class FakeChainSim:
    def __init__(self) -> None:
        self.z = make_chain(4)
        self.n = 4
        self.rold = [1.0, 1.5, 2.0, 2.5]
        self.xold = [[0.0, 0.0] for _ in range(4)]
        self.data = [(100.0, 0.0)]


class FakeAntichainSim:
    def __init__(self) -> None:
        self.z = make_antichain(4)
        self.n = 4
        self.rold = [1.0, 1.0, 1.0, 1.0]
        self.xold = [[0.0, 0.0], [2.0, 0.0], [4.0, 0.0], [6.0, 0.0]]
        self.data = [(100.0, 0.0)]


class TestSorkin2KnownTruthHarness(unittest.TestCase):
    def assert_png_exists_and_nonempty(self, path: Path) -> None:
        self.assertTrue(path.exists(), msg=f"{path} should exist")
        self.assertGreater(path.stat().st_size, 0, msg=f"{path} should be nonempty")
        self.assertEqual(path.read_bytes()[:8], b"\x89PNG\r\n\x1a\n")

    def test_run_id_and_result_directory_construction(self) -> None:
        now = datetime(2026, 5, 22, 10, 11, 12, tzinfo=timezone.utc)
        run_id = harness.make_run_id("chain_4_d2", 1959, now)
        self.assertEqual(run_id, "20260522T101112Z_chain_4_d2_seed1959")

        root = Path("/tmp/sorkin2-test-root")
        run_dir = harness.build_run_dir(root, "chain_4_d2", run_id)
        self.assertEqual(
            run_dir,
            root / "chain_4_d2" / "20260522T101112Z_chain_4_d2_seed1959",
        )

    def test_enabled_case_ids_allowed(self) -> None:
        self.assertEqual(harness.get_enabled_case("chain_4_d2").case_id, "chain_4_d2")
        self.assertEqual(
            harness.get_enabled_case("antichain_4_d2").case_id,
            "antichain_4_d2",
        )
        chain12 = harness.get_enabled_case("chain_12_d2")
        self.assertEqual(chain12.case_id, "chain_12_d2")
        self.assertEqual(chain12.annealer_mode, "historical/default")
        self.assertTrue(str(chain12.input_file).endswith("chain_12_d2.in"))
        antichain12 = harness.get_enabled_case("antichain_12_d2")
        self.assertEqual(antichain12.case_id, "antichain_12_d2")
        self.assertEqual(antichain12.annealer_mode, "historical/default")
        self.assertTrue(str(antichain12.input_file).endswith("antichain_12_d2.in"))
        layered12 = harness.get_enabled_case("layered_4_4_4_d2")
        self.assertEqual(layered12.case_id, "layered_4_4_4_d2")
        self.assertEqual(layered12.annealer_mode, "historical/default")
        self.assertTrue(str(layered12.input_file).endswith("layered_4_4_4_d2.in"))
        layered12_g08 = harness.get_enabled_case("layered_4_4_4_d2_T100_g08")
        self.assertEqual(layered12_g08.case_id, "layered_4_4_4_d2_T100_g08")
        self.assertEqual(layered12_g08.annealer_mode, "mechanism/T100_g08")
        self.assertTrue(str(layered12_g08.input_file).endswith("layered_4_4_4_d2.in"))
        self.assertEqual(layered12_g08.initial_temp, 100.0)
        self.assertEqual(layered12_g08.cooling_factor, 0.8)
        case6 = harness.get_enabled_case("minkowski_6_s1959_d2")
        self.assertEqual(case6.case_id, "minkowski_6_s1959_d2")
        self.assertEqual(case6.annealer_mode, "historical/default")
        case12h = harness.get_enabled_case("minkowski_12_s1962_d2_hist")
        self.assertEqual(case12h.case_id, "minkowski_12_s1962_d2_hist")
        self.assertEqual(case12h.annealer_mode, "historical/default")
        case12t = harness.get_enabled_case("minkowski_12_s1962_d2_tuned")
        self.assertEqual(case12t.case_id, "minkowski_12_s1962_d2_tuned")
        self.assertEqual(case12t.annealer_mode, "tuned/non-historical")
        self.assertEqual(case12t.initial_temp, 180.0)
        self.assertEqual(case12t.cooling_factor, 0.8)
        case_t180 = harness.get_enabled_case("minkowski_12_s1962_d2_T180_g09")
        self.assertEqual(case_t180.case_id, "minkowski_12_s1962_d2_T180_g09")
        self.assertEqual(case_t180.annealer_mode, "mechanism/T180_g09")
        self.assertEqual(case_t180.initial_temp, 180.0)
        self.assertEqual(case_t180.cooling_factor, 0.9)
        case_t100 = harness.get_enabled_case("minkowski_12_s1962_d2_T100_g08")
        self.assertEqual(case_t100.case_id, "minkowski_12_s1962_d2_T100_g08")
        self.assertEqual(case_t100.annealer_mode, "mechanism/T100_g08")
        self.assertEqual(case_t100.initial_temp, 100.0)
        self.assertEqual(case_t100.cooling_factor, 0.8)

    def test_disabled_case_id_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "unsupported"):
            harness.get_enabled_case("minkowski_6")
        with self.assertRaisesRegex(ValueError, "unsupported"):
            harness.get_enabled_case("minkowski_12_s1962_d2")
        with self.assertRaisesRegex(ValueError, "unsupported"):
            harness.get_enabled_case("minkowski_12_sparse_d2_seedTBD")
        with self.assertRaisesRegex(ValueError, "unsupported"):
            harness.get_enabled_case("minkowski_12_hub_d2_seedTBD")
        with self.assertRaisesRegex(ValueError, "unsupported"):
            harness.get_enabled_case("chain_12_d2_T100_g08")
        with self.assertRaisesRegex(ValueError, "unsupported"):
            harness.get_enabled_case("layered_4_4_4_d2_tuned")
        with self.assertRaisesRegex(ValueError, "unsupported"):
            harness.get_enabled_case("layered_4_4_4_d2_T180_g09")
        with self.assertRaisesRegex(ValueError, "unsupported"):
            harness.get_enabled_case("unknown_case")

    def test_writes_result_json_and_required_figures(self) -> None:
        cases = [
            (harness.get_enabled_case("chain_4_d2"), FakeChainSim()),
            (harness.get_enabled_case("antichain_4_d2"), FakeAntichainSim()),
        ]
        for case, sim in cases:
            with self.subTest(case_id=case.case_id), tempfile.TemporaryDirectory() as tmpdir:
                run_dir = Path(tmpdir) / "results" / case.case_id / "unit"
                run_dir.mkdir(parents=True)

                result = harness.write_result_artifacts(
                    case=case,
                    seed=1959,
                    run_id="unit",
                    run_dir=run_dir,
                    sim=sim,
                )

                result_path = run_dir / "result.json"
                self.assertTrue(result_path.exists())
                loaded = json.loads(result_path.read_text(encoding="utf-8"))
                self.assertEqual(loaded, result)

                required = {
                    "case_id",
                    "input_file",
                    "annealer_mode",
                    "seed",
                    "n",
                    "final_energy",
                    "exact_match",
                    "total_relations_target",
                    "total_relations_induced",
                    "missing_relations",
                    "extra_relations",
                    "generated_figures",
                    "primary_recovery_criterion_warning",
                }
                self.assertTrue(required.issubset(loaded))
                self.assertEqual(loaded["case_id"], case.case_id)
                self.assertTrue(loaded["exact_match"])
                self.assertEqual(loaded["missing_relations"], [])
                self.assertEqual(loaded["extra_relations"], [])
                self.assertIn("exact_match", loaded["primary_recovery_criterion_warning"])

                figures = loaded["generated_figures"]
                self.assertEqual(figures["target_order_matrix"], "target_order_matrix.png")
                self.assertEqual(figures["induced_order_matrix"], "induced_order_matrix.png")
                self.assertEqual(figures["order_difference_matrix"], "order_difference_matrix.png")
                self.assert_png_exists_and_nonempty(run_dir / figures["target_order_matrix"])
                self.assert_png_exists_and_nonempty(run_dir / figures["induced_order_matrix"])
                self.assert_png_exists_and_nonempty(run_dir / figures["order_difference_matrix"])

                trace = loaded.get("trace_artifacts", {})
                self.assertIn("trace_csv", trace)
                self.assertIn("trace_pairs", trace)
                self.assertIn("trace_coordinates", trace)
                self.assertIn("trace_energy", trace)
                self.assertIn("trace_relations", trace)
                self.assertTrue((run_dir / trace["trace_pairs"]).exists())
                self.assertTrue((run_dir / trace["trace_coordinates"]).exists())
                self.assert_png_exists_and_nonempty(run_dir / trace["trace_energy"])
                self.assert_png_exists_and_nonempty(run_dir / trace["trace_relations"])

    def test_write_trace_artifacts_csv_columns_and_pngs(self) -> None:
        import csv as csv_mod
        import validation_suite as vs

        sim = FakeChainSim()
        comparison = vs.verify_recovery(sim)
        with tempfile.TemporaryDirectory() as tmpdir:
            run_dir = Path(tmpdir)
            artifacts = harness.write_trace_artifacts(
                sim=sim, comparison=comparison, run_dir=run_dir
            )

            # CSV columns
            trace_csv = run_dir / artifacts["trace_csv"]
            self.assertTrue(trace_csv.exists())
            with trace_csv.open(newline="", encoding="utf-8") as fh:
                reader = csv_mod.DictReader(fh)
                self.assertEqual(reader.fieldnames, harness._TRACE_FIELDS)
                rows = list(reader)
            self.assertEqual(len(rows), len(sim.data))
            last = rows[-1]
            self.assertIn(last["exact_match"], ("true", "false"))
            self.assertNotEqual(last["induced_relations"], "")
            self.assertNotEqual(last["missing_relations_count"], "")
            self.assertNotEqual(last["extra_relations_count"], "")

            # PNG figures
            self.assert_png_exists_and_nonempty(run_dir / artifacts["trace_energy"])
            self.assert_png_exists_and_nonempty(run_dir / artifacts["trace_relations"])
            self.assertTrue((run_dir / artifacts["trace_pairs"]).exists())
            self.assertTrue((run_dir / artifacts["trace_coordinates"]).exists())

    def test_write_trace_artifacts_multiblock(self) -> None:
        import csv as csv_mod
        import validation_suite as vs

        sim = FakeChainSim()
        sim.data = [(100.0, 5.0), (50.0, 2.5), (25.0, 0.5), (12.5, 0.0)]
        comparison = vs.verify_recovery(sim)
        with tempfile.TemporaryDirectory() as tmpdir:
            run_dir = Path(tmpdir)
            artifacts = harness.write_trace_artifacts(
                sim=sim, comparison=comparison, run_dir=run_dir
            )
            trace_csv = run_dir / artifacts["trace_csv"]
            with trace_csv.open(newline="", encoding="utf-8") as fh:
                rows = list(csv_mod.DictReader(fh))
            self.assertEqual(len(rows), 4)
            # Intermediate blocks have empty causal-order columns
            for row in rows[:-1]:
                self.assertEqual(row["induced_relations"], "")
                self.assertEqual(row["missing_relations_count"], "")
                self.assertEqual(row["extra_relations_count"], "")
                self.assertEqual(row["exact_match"], "")
            # Final block is populated
            last = rows[-1]
            self.assertNotEqual(last["induced_relations"], "")
            self.assertIn(last["exact_match"], ("true", "false"))


class TestBlockCallback(unittest.TestCase):
    def test_callback_invoked_per_block(self) -> None:
        import contextlib
        import io as io_mod
        import cones

        z = make_chain(4)
        calls: list[tuple[int, float, float, int]] = []

        def cb(sim: Any, block_idx: int, temp: float, eave: float) -> None:
            calls.append((block_idx, temp, eave, len(sim.rold)))

        with tempfile.TemporaryDirectory() as tmpdir:
            sim = cones.ConesSimulator(
                z=z, dim=2, seed=1959, interactive=False,
                max_data=3, warmup_limit=10, anneal_limit=10,
                block_callback=cb,
            )
            with contextlib.redirect_stdout(io_mod.StringIO()):
                sim.run(Path(tmpdir) / "out.txt")

        self.assertEqual(len(calls), len(sim.data))
        for i, (block_idx, temp, eave, n_rold) in enumerate(calls):
            self.assertEqual(block_idx, i + 1)
            self.assertEqual(n_rold, 4)

    def test_callback_none_does_not_change_result(self) -> None:
        import contextlib
        import io as io_mod
        import cones

        z = make_chain(4)
        with tempfile.TemporaryDirectory() as tmpdir:
            sim_no_cb = cones.ConesSimulator(
                z=z, dim=2, seed=1959, interactive=False,
                max_data=5, warmup_limit=10, anneal_limit=10,
            )
            with contextlib.redirect_stdout(io_mod.StringIO()):
                sim_no_cb.run(Path(tmpdir) / "out_no_cb.txt")

            sim_with_cb = cones.ConesSimulator(
                z=z, dim=2, seed=1959, interactive=False,
                max_data=5, warmup_limit=10, anneal_limit=10,
                block_callback=lambda s, b, t, e: None,
            )
            with contextlib.redirect_stdout(io_mod.StringIO()):
                sim_with_cb.run(Path(tmpdir) / "out_with_cb.txt")

        self.assertEqual(sim_no_cb.data, sim_with_cb.data)
        self.assertAlmostEqual(sim_no_cb.eave, sim_with_cb.eave)

    def test_write_trace_artifacts_with_full_block_records(self) -> None:
        import csv as csv_mod
        import validation_suite as vs
        from typing import Any

        sim = FakeChainSim()
        comparison = vs.verify_recovery(sim)
        block_records: list[dict[str, Any]] = [
            {
                "block": 1, "temperature": 100.0, "energy": 3.0,
                "induced_relations": 2, "correct_relations": 2,
                "missing_relations_count": 4, "extra_relations_count": 0,
                "exact_match": False,
                "correct_pairs": [[0, 1], [1, 2]],
                "missing_pairs": [[0, 2], [0, 3], [1, 3], [2, 3]],
                "extra_pairs": [],
                "coordinates": [
                    {"i": 0, "R": 1.0, "X": [0.0, 0.0]},
                    {"i": 1, "R": 1.5, "X": [0.0, 0.0]},
                    {"i": 2, "R": 2.0, "X": [0.0, 0.0]},
                    {"i": 3, "R": 2.5, "X": [0.0, 0.0]},
                ],
            },
            {
                "block": 2, "temperature": 90.0, "energy": 1.0,
                "induced_relations": 5, "correct_relations": 5,
                "missing_relations_count": 1, "extra_relations_count": 0,
                "exact_match": False,
                "correct_pairs": [[0, 1], [0, 2], [1, 2], [1, 3], [2, 3]],
                "missing_pairs": [[0, 3]],
                "extra_pairs": [],
                "coordinates": [
                    {"i": 0, "R": 1.0, "X": [0.0, 0.0]},
                    {"i": 1, "R": 1.5, "X": [0.0, 0.0]},
                    {"i": 2, "R": 2.0, "X": [0.0, 0.0]},
                    {"i": 3, "R": 2.5, "X": [0.0, 0.0]},
                ],
            },
            {
                "block": 3, "temperature": 81.0, "energy": 0.0,
                "induced_relations": 6, "correct_relations": 6,
                "missing_relations_count": 0, "extra_relations_count": 0,
                "exact_match": True,
                "correct_pairs": [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]],
                "missing_pairs": [],
                "extra_pairs": [],
                "coordinates": [
                    {"i": 0, "R": 1.0, "X": [0.0, 0.0]},
                    {"i": 1, "R": 1.5, "X": [0.0, 0.0]},
                    {"i": 2, "R": 2.0, "X": [0.0, 0.0]},
                    {"i": 3, "R": 2.5, "X": [0.0, 0.0]},
                ],
            },
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            run_dir = Path(tmpdir)
            artifacts = harness.write_trace_artifacts(
                sim=sim, comparison=comparison, run_dir=run_dir,
                block_records=block_records,
            )
            trace_csv = run_dir / artifacts["trace_csv"]
            with trace_csv.open(newline="", encoding="utf-8") as fh:
                rows = list(csv_mod.DictReader(fh))
            trace_pairs = run_dir / artifacts["trace_pairs"]
            trace_coordinates = run_dir / artifacts["trace_coordinates"]
            pair_rows = [
                json.loads(line)
                for line in trace_pairs.read_text(encoding="utf-8").splitlines()
            ]
            coordinate_rows = [
                json.loads(line)
                for line in trace_coordinates.read_text(encoding="utf-8").splitlines()
            ]

        self.assertEqual(len(rows), 3)
        self.assertEqual(len(pair_rows), 3)
        self.assertEqual(len(coordinate_rows), 3)
        for row in rows:
            self.assertNotEqual(row["induced_relations"], "")
            self.assertNotEqual(row["missing_relations_count"], "")
            self.assertNotEqual(row["exact_match"], "")
        for pair_row in pair_rows:
            self.assertIn("correct_pairs", pair_row)
            self.assertIn("missing_pairs", pair_row)
            self.assertIn("extra_pairs", pair_row)
            self.assertIn("correct_relations", pair_row)
        for coordinate_row in coordinate_rows:
            self.assertIn("coordinates", coordinate_row)
            self.assertEqual(len(coordinate_row["coordinates"]), 4)
            self.assertIn("R", coordinate_row["coordinates"][0])
            self.assertIn("X", coordinate_row["coordinates"][0])
        self.assertEqual(rows[0]["exact_match"], "false")
        self.assertEqual(rows[-1]["exact_match"], "true")
        self.assertEqual(pair_rows[0]["correct_pairs"], [[0, 1], [1, 2]])
        self.assertEqual(pair_rows[-1]["missing_pairs"], [])


if __name__ == "__main__":
    unittest.main()
