from __future__ import annotations

import json
import math
import unittest
from pathlib import Path

from explore.sorkin4_kerr_benchmark import audit_kerr_k1_control as audit


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT
    / "explore"
    / "sorkin4_kerr_benchmark"
    / "kerr_k1_control_spin_sweep_n12_seed1959.json"
)


class KerrK1ControlAuditTests(unittest.TestCase):
    def test_kerr_k1_control_sweep_computation(self) -> None:
        payload = audit.run_audit()
        aggregate = payload["aggregate"]
        rows = payload["rows"]
        cases = payload["cases"]

        self.assertIs(aggregate["all_checks_pass"], True)
        self.assertEqual(aggregate["M"], 1.0)
        self.assertEqual(aggregate["spins"], [0.0, 0.25, 0.5, 0.75])
        self.assertEqual(len(rows), 4)
        self.assertEqual(len(cases), 4)

        a0 = rows[0]
        self.assertEqual(a0["a"], 0.0)
        self.assertEqual(a0["true_relations"], 1)
        self.assertEqual(a0["false_relations"], 64)
        self.assertEqual(a0["undecided_pairs"], 1)
        self.assertIs(a0["case_checks_pass"], True)

        possible_pairs = aggregate["possible_pairs"]
        for row in rows[1:]:
            spin = row["a"]
            self.assertLess(abs(spin), aggregate["M"])
            self.assertAlmostEqual(row["r_plus"], 1.0 + math.sqrt(1.0 - spin * spin))
            self.assertEqual(row["true_relations"], 0)
            self.assertEqual(row["false_relations"], 0)
            self.assertEqual(row["undecided_pairs"], possible_pairs)
            self.assertEqual(row["decided_pairs"], 0)
            self.assertIs(row["case_checks_pass"], True)

    def test_kerr_k1_control_artifact_is_current(self) -> None:
        self.assertTrue(
            ARTIFACT.exists(),
            "missing Kerr K1 control artifact; run "
            "`python explore/sorkin4_kerr_benchmark/audit_kerr_k1_control.py`",
        )
        payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        aggregate = payload["aggregate"]
        rows = payload["rows"]

        self.assertEqual(aggregate["benchmark"], "S4-K1 Kerr scaffold invariant control audit")
        self.assertEqual(aggregate["N"], audit.DEFAULT_N)
        self.assertEqual(aggregate["seed"], audit.DEFAULT_SEED)
        self.assertEqual(aggregate["M"], audit.DEFAULT_MASS)
        self.assertEqual(aggregate["spins"], list(audit.DEFAULT_SPINS))
        self.assertIs(aggregate["all_checks_pass"], True)
        self.assertIs(aggregate["a0_exact_schwarzschild_control"], True)
        self.assertIs(aggregate["positive_spin_cases_all_undecided"], True)

        self.assertEqual([row["a"] for row in rows], list(audit.DEFAULT_SPINS))
        self.assertTrue(all(row["all_checks_pass"] for row in rows))


if __name__ == "__main__":
    unittest.main()
