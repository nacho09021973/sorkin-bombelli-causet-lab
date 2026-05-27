from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "explore" / "sorkin4_schwarzschild_benchmark"))

from explore.sorkin4_schwarzschild_benchmark import audit_exterior_turning_asymptotic as audit


AUDIT_JSON = (
    ROOT
    / "explore"
    / "sorkin4_schwarzschild_benchmark"
    / "schwarzschild_exterior_turning_asymptotic_audit.json"
)


class SchwarzschildExteriorTurningAsymptoticAuditTests(unittest.TestCase):
    def test_default_artifact_records_positive_asymptotic_gap(self) -> None:
        self.assertTrue(AUDIT_JSON.exists())
        payload = json.loads(AUDIT_JSON.read_text(encoding="utf-8"))
        summary = payload["summary"]

        self.assertEqual(summary["grid_points"], 70)
        self.assertEqual(summary["status_counts"], {"disjoint_ranges": 70})
        self.assertTrue(summary["all_ranges_disjoint"])
        self.assertGreater(summary["min_turning_minus_direct_phi_gap"], 0.0)

    def test_small_asymptotic_grid_has_positive_gap(self) -> None:
        rows = audit.audit_asymptotic_grid(
            r_values=[20.0, 50.0],
            eps_values=[0.001, 0.01],
        )

        self.assertEqual(len(rows), 4)
        self.assertEqual({row["status"] for row in rows}, {"disjoint_ranges"})
        self.assertGreater(
            min(row["turning_minus_direct_phi_gap"] for row in rows),
            0.0,
        )


if __name__ == "__main__":
    unittest.main()
