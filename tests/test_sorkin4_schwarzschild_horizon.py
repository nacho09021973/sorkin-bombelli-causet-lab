from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from explore.sorkin4_schwarzschild_benchmark import run_schwarzschild_horizon_benchmark as horizon


class SchwarzschildHorizonShootingTests(unittest.TestCase):
    def test_seed1959_generic_horizon_shooting_control(self) -> None:
        _, _, _, summary = horizon.run_horizon_case(
            8,
            4,
            1959,
            enable_horizon_shooting=True,
        )

        self.assertEqual(summary["status"], "minimal_ief_radial_criterion_horizon_shooting")
        self.assertEqual(summary["true_relations"], 0)
        self.assertEqual(summary["undecided_pairs"], 3)
        self.assertEqual(summary["horizon_crossing_links"], 0)
        self.assertTrue(summary["antisymmetric"])
        self.assertTrue(summary["transitive"])

    def test_seed31_generic_horizon_shooting_adds_crossing_links(self) -> None:
        _, _, _, baseline = horizon.run_horizon_case(8, 4, 31)
        _, _, _, shooting = horizon.run_horizon_case(
            8,
            4,
            31,
            enable_horizon_shooting=True,
        )

        self.assertEqual(baseline["status"], "minimal_ief_radial_criterion_nonradial_undecided")
        self.assertEqual(baseline["true_relations"], 0)
        self.assertEqual(baseline["undecided_pairs"], 14)
        self.assertEqual(baseline["horizon_crossing_links"], 0)

        self.assertEqual(shooting["status"], "minimal_ief_radial_criterion_horizon_shooting")
        self.assertIn("non-radial plunging-null shooting", str(shooting["causal_model"]))
        self.assertEqual(shooting["true_relations"], 3)
        self.assertEqual(shooting["undecided_pairs"], 3)
        self.assertEqual(shooting["horizon_crossing_links"], 3)
        self.assertTrue(shooting["antisymmetric"])
        self.assertTrue(shooting["transitive"])


if __name__ == "__main__":
    unittest.main()
