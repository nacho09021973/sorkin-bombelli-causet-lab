from __future__ import annotations

import json
import unittest
from pathlib import Path

from explore.sorkin4_kerr_benchmark import audit_kerr_l0_scaffold as audit


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT
    / "explore"
    / "sorkin4_kerr_benchmark"
    / "kerr_l0_scaffold_control_a0_a0p5_n12_seed1959.json"
)


class KerrL0ScaffoldAuditTests(unittest.TestCase):
    def test_kerr_l0_scaffold_audit_computation(self) -> None:
        payload = audit.run_audit()
        summary = payload["summary"]
        a0_checks = payload["a0_checks"]
        spin_checks = payload["spin_checks"]

        self.assertIs(summary["all_checks_pass"], True)
        self.assertIs(summary["a0_checks_pass"], True)
        self.assertIs(summary["spin_checks_pass"], True)
        self.assertTrue(all(a0_checks.values()))
        self.assertTrue(all(spin_checks.values()))

        self.assertEqual(summary["a0_true_relations"], 1)
        self.assertEqual(summary["a0_false_relations"], 64)
        self.assertEqual(summary["a0_undecided_pairs"], 1)
        self.assertEqual(summary["a0_decided_pairs"], 65)

        self.assertEqual(summary["spin_true_relations"], 0)
        self.assertEqual(summary["spin_false_relations"], 0)
        self.assertEqual(summary["spin_undecided_pairs"], 66)
        self.assertEqual(summary["spin_decided_pairs"], 0)

    def test_kerr_l0_scaffold_artifact_is_current(self) -> None:
        self.assertTrue(
            ARTIFACT.exists(),
            "missing Kerr L0 scaffold audit artifact; run "
            "`python explore/sorkin4_kerr_benchmark/audit_kerr_l0_scaffold.py`",
        )
        payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        summary = payload["summary"]

        self.assertEqual(summary["benchmark"], "S4-K0 Kerr L0 scaffold control audit")
        self.assertEqual(summary["N"], audit.DEFAULT_N)
        self.assertEqual(summary["seed"], audit.DEFAULT_SEED)
        self.assertIs(summary["all_checks_pass"], True)
        self.assertIs(summary["a0_events_match_schwarzschild"], True)
        self.assertIs(summary["a0_relation_states_match_schwarzschild"], True)
        self.assertEqual(summary["spin_status"], "kerr_scaffold_pairs_undecided")
        self.assertEqual(summary["spin_undecided_pairs"], summary["possible_pairs"])


if __name__ == "__main__":
    unittest.main()
