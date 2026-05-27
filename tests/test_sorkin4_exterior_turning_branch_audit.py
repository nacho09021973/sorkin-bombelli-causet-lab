from __future__ import annotations

import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "explore" / "sorkin4_schwarzschild_benchmark"))

from explore.sorkin4_schwarzschild_benchmark import audit_exterior_turning_branch as audit


class SchwarzschildExteriorTurningBranchAuditTests(unittest.TestCase):
    def test_current_short_sweep_has_no_turning_competitor(self) -> None:
        rows = audit.collect_rows(12, 1959, 1968)

        self.assertEqual(len(rows), 64)
        self.assertEqual({row["outcome"] for row in rows}, {"no_turning_solution"})

    def test_synthetic_turning_solver_finds_large_angle_branch(self) -> None:
        u1 = 1.0 / 4.0
        u2 = 1.0 / 6.0
        found = audit.find_turning_c2(u1, u2, 4.0)

        self.assertIsNotNone(found)
        assert found is not None
        c2, phi_value, u_turn = found
        self.assertLess(c2, 1.0 / 27.0)
        self.assertGreater(u_turn, u1)
        self.assertLess(u_turn, 1.0 / 3.0)
        self.assertAlmostEqual(phi_value, 4.0, places=5)


if __name__ == "__main__":
    unittest.main()
