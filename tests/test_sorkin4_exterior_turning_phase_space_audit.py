from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "explore" / "sorkin4_schwarzschild_benchmark"))

from explore.sorkin4_schwarzschild_benchmark import audit_exterior_turning_phase_space as audit


AUDIT_JSON = (
    ROOT
    / "explore"
    / "sorkin4_schwarzschild_benchmark"
    / "schwarzschild_exterior_turning_phase_space_audit.json"
)


class SchwarzschildExteriorTurningPhaseSpaceAuditTests(unittest.TestCase):
    def test_default_artifact_records_disjoint_ranges(self) -> None:
        self.assertTrue(AUDIT_JSON.exists())
        payload = json.loads(AUDIT_JSON.read_text(encoding="utf-8"))
        summary = payload["summary"]

        self.assertEqual(summary["grid_points"], 11205)
        self.assertEqual(summary["status_counts"], {"disjoint_ranges": 11205})
        self.assertTrue(summary["all_ranges_disjoint"])
        self.assertGreater(summary["min_turning_minus_direct_phi_gap"], 0.0)

    def test_small_grid_has_positive_gap(self) -> None:
        rows = audit.audit_grid(r1_min=3.1, r1_max=3.5, r2_max=4.0, step=0.1)

        self.assertGreater(len(rows), 0)
        self.assertEqual({row["status"] for row in rows}, {"disjoint_ranges"})
        self.assertGreater(
            min(row["turning_minus_direct_phi_gap"] for row in rows),
            0.0,
        )


if __name__ == "__main__":
    unittest.main()
