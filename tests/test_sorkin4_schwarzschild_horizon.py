from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

AUDIT_JSON = (
    ROOT
    / "explore"
    / "sorkin4_schwarzschild_benchmark"
    / "schwarzschild_horizon_shooting_branch_audit.json"
)

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


class SchwarzschildHorizonShootingBranchAuditTests(unittest.TestCase):
    def setUp(self) -> None:
        self.assertTrue(
            AUDIT_JSON.exists(),
            msg=(
                f"missing S4 horizon shooting branch audit at {AUDIT_JSON}; "
                "run `python explore/sorkin4_schwarzschild_benchmark/"
                "audit_horizon_shooting_branch.py`"
            ),
        )
        with AUDIT_JSON.open(encoding="utf-8") as fh:
            self.audit = json.load(fh)

    def test_seed_range_and_crossing_link_count(self) -> None:
        self.assertEqual(self.audit["seed_start"], 1)
        self.assertEqual(self.audit["seed_stop"], 40)
        self.assertEqual(self.audit["audited_crossing_links"], 16)
        self.assertEqual(len(self.audit["rows"]), 16)

    def test_summary_checks_pass_with_audited_tolerances(self) -> None:
        self.assertTrue(self.audit["all_checks_pass"])
        self.assertLessEqual(self.audit["max_phi_error"], 1.0e-8)
        self.assertLessEqual(self.audit["max_dt_refine_rel_1024_2048"], 1.0e-8)
        self.assertLessEqual(self.audit["max_exterior_regular_minus_raw_abs"], 1.0e-12)

    def test_each_crossing_link_passes_branch_audit_checks(self) -> None:
        for row in self.audit["rows"]:
            self.assertGreater(row["c2_over_critical"], 1.0)
            self.assertEqual(row["positive_root_count"], 0)
            self.assertLessEqual(row["phi_error"], 1.0e-8)
            self.assertLessEqual(row["dt_refine_rel_1024_2048"], 1.0e-8)
            self.assertTrue(row["local_phi_monotone"])
            self.assertGreaterEqual(row["related_margin"], -horizon.TIME_EPS)
            self.assertIsNotNone(row["exterior_regular_minus_raw_abs"])
            self.assertLessEqual(row["exterior_regular_minus_raw_abs"], 1.0e-12)


if __name__ == "__main__":
    unittest.main()
