from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "explore" / "sorkin4_schwarzschild_benchmark"))

from explore.sorkin4_schwarzschild_benchmark import run_schwarzschild_area_sweep as area


AREA_SWEEP_JSON = (
    ROOT
    / "explore"
    / "sorkin4_schwarzschild_benchmark"
    / "schwarzschild_horizon_area_sweep.json"
)


class SchwarzschildAreaSweepTests(unittest.TestCase):
    def test_default_artifact_records_monotone_aligned_sweep(self) -> None:
        self.assertTrue(AREA_SWEEP_JSON.exists())
        payload = json.loads(AREA_SWEEP_JSON.read_text(encoding="utf-8"))
        summary = payload["summary"]

        self.assertEqual(summary["mass_values"], [0.75, 1.0, 1.25, 1.5, 1.75, 2.0])
        self.assertEqual(summary["seed_start"], 1)
        self.assertEqual(summary["seed_stop"], 40)
        self.assertEqual(summary["N_exterior"], 16)
        self.assertEqual(summary["N_interior"], 8)
        self.assertTrue(summary["monotone_non_decreasing_mean_horizon_links"])
        self.assertEqual(summary["failed_order_check_count"], 0)

        means = [row["mean_horizon_crossing_links"] for row in payload["aggregate_rows"]]
        self.assertEqual(means, [2.775, 3.375, 3.8, 4.575, 4.9, 5.05])

    def test_tiny_sweep_preserves_order_checks(self) -> None:
        rows, per_seed = area.run_area_sweep(
            masses=[1.0, 1.5],
            seed_start=1,
            seed_stop=3,
            n_exterior=8,
            n_interior=4,
        )

        self.assertEqual(len(rows), 2)
        self.assertEqual(len(per_seed), 6)
        self.assertTrue(all(row["failed_order_check_count"] == 0 for row in rows))
        self.assertTrue(all(seed_row["antisymmetric"] for seed_row in per_seed))
        self.assertTrue(all(seed_row["transitive"] for seed_row in per_seed))


if __name__ == "__main__":
    unittest.main()
