"""Tests for S4-KERR-K6-ZAMO-OMEGA-HORIZON-001.

Test requirements:
 1.  Artifact CSV/JSON/MD/PNG files exist.
 2.  CSV is parseable and has exactly 25 rows (5 spins × 5 deltas).
 3.  JSON is parseable; has "aggregate", "spin_summaries", "rows".
 4.  MD exists.
 5.  PNG exists and has nonzero size.
 6.  JSON aggregate all_checks_pass=True.
 7.  Spin list is exactly [0.0, 0.25, 0.5, 0.75, 0.9].
 8.  Delta list is exactly [1e-1, 3e-2, 1e-2, 3e-3, 1e-3].
 9.  All CSV rows have exterior_pass=True.
10.  For a=0: convergence_monotone_pass=True, omega_zamo_below_Omega_H=True
     (both trivially True by convention).
11.  For a>0: convergence_monotone_pass=True (residuals decrease with delta).
12.  For a>0: omega_zamo_below_Omega_H=True.
13.  For a>0: global causal accounting invariant
     (global_true_relations=0, global_false_relations=0,
      global_undecided_pairs=N*(N-1)/2).
14.  Existing K1–K5 and Schwarzschild stability tests still pass
     (covered by the full validation command; not duplicated here).
"""
from __future__ import annotations

import csv
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "explore" / "sorkin4_kerr_benchmark"))

import audit_kerr_k6_zamo_omega_horizon_001 as audit

ARTIFACT_DIR  = ROOT / "explore" / "sorkin4_kerr_benchmark"
ARTIFACT_CSV  = ARTIFACT_DIR / "kerr_k6_zamo_omega_horizon_001_n12_seed1959.csv"
ARTIFACT_JSON = ARTIFACT_DIR / "kerr_k6_zamo_omega_horizon_001_n12_seed1959.json"
ARTIFACT_MD   = ARTIFACT_DIR / "kerr_k6_zamo_omega_horizon_001_n12_seed1959.md"
ARTIFACT_PNG  = ARTIFACT_DIR / "kerr_k6_zamo_omega_horizon_001_n12_seed1959.png"

EXPECTED_SPINS  = [0.0, 0.25, 0.5, 0.75, 0.9]
EXPECTED_DELTAS = [1e-1, 3e-2, 1e-2, 3e-3, 1e-3]
EXPECTED_ROWS   = len(EXPECTED_SPINS) * len(EXPECTED_DELTAS)  # 25

_MISSING_MSG = (
    "Missing K6 artifact. Run: python3 explore/sorkin4_kerr_benchmark/"
    "audit_kerr_k6_zamo_omega_horizon_001.py"
)


class KerrK6ArtifactFileTests(unittest.TestCase):
    """Tests 1–8: artifact files exist, are parseable, and have correct metadata."""

    def test_csv_exists_and_parseable(self) -> None:
        self.assertTrue(ARTIFACT_CSV.exists(), msg=_MISSING_MSG)
        with ARTIFACT_CSV.open(encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        self.assertEqual(len(rows), EXPECTED_ROWS)

    def test_json_exists_and_parseable(self) -> None:
        self.assertTrue(ARTIFACT_JSON.exists(), msg=_MISSING_MSG)
        with ARTIFACT_JSON.open(encoding="utf-8") as fh:
            payload = json.load(fh)
        self.assertIn("aggregate", payload)
        self.assertIn("spin_summaries", payload)
        self.assertIn("rows", payload)

    def test_md_exists(self) -> None:
        self.assertTrue(ARTIFACT_MD.exists(), msg=_MISSING_MSG)

    def test_png_exists_and_nonempty(self) -> None:
        self.assertTrue(ARTIFACT_PNG.exists(), msg=_MISSING_MSG)
        self.assertGreater(ARTIFACT_PNG.stat().st_size, 0)

    def test_json_all_checks_pass(self) -> None:
        self.assertTrue(ARTIFACT_JSON.exists(), msg=_MISSING_MSG)
        with ARTIFACT_JSON.open(encoding="utf-8") as fh:
            payload = json.load(fh)
        self.assertTrue(
            payload["aggregate"]["all_checks_pass"],
            msg="aggregate all_checks_pass is False in frozen artifact",
        )

    def test_json_spin_list(self) -> None:
        self.assertTrue(ARTIFACT_JSON.exists(), msg=_MISSING_MSG)
        with ARTIFACT_JSON.open(encoding="utf-8") as fh:
            payload = json.load(fh)
        self.assertEqual(
            payload["aggregate"]["spins"],
            EXPECTED_SPINS,
            msg="Spin list in frozen artifact does not match expected sweep",
        )

    def test_json_delta_list(self) -> None:
        self.assertTrue(ARTIFACT_JSON.exists(), msg=_MISSING_MSG)
        with ARTIFACT_JSON.open(encoding="utf-8") as fh:
            payload = json.load(fh)
        self.assertEqual(
            payload["aggregate"]["deltas"],
            EXPECTED_DELTAS,
            msg="Delta list in frozen artifact does not match expected sweep",
        )

    def test_json_row_count(self) -> None:
        self.assertTrue(ARTIFACT_JSON.exists(), msg=_MISSING_MSG)
        with ARTIFACT_JSON.open(encoding="utf-8") as fh:
            payload = json.load(fh)
        self.assertEqual(len(payload["rows"]), EXPECTED_ROWS)


class KerrK6ArtifactRowTests(unittest.TestCase):
    """Tests 9–13: per-row checks on the frozen JSON artifact."""

    def setUp(self) -> None:
        self.assertTrue(ARTIFACT_JSON.exists(), msg=_MISSING_MSG)
        with ARTIFACT_JSON.open(encoding="utf-8") as fh:
            self.payload = json.load(fh)
        self.rows          = self.payload["rows"]
        self.spin_summaries = self.payload["spin_summaries"]
        self.n             = self.payload["aggregate"]["N"]

    # Test 9
    def test_all_rows_exterior_pass(self) -> None:
        for row in self.rows:
            self.assertTrue(
                row["exterior_pass"],
                msg=f"exterior_pass=False for a={row['spin_a']}, delta={row['delta']}",
            )

    # Test 10
    def test_a0_convergence_mono_pass_by_convention(self) -> None:
        a0 = next(s for s in self.spin_summaries if s["spin_a"] == 0.0)
        self.assertTrue(a0["convergence_monotone_pass"],
                        msg="convergence_monotone_pass must be True at a=0 (by convention)")
        self.assertTrue(a0["omega_zamo_below_Omega_H"],
                        msg="omega_zamo_below_Omega_H must be True at a=0 (by convention)")
        self.assertTrue(a0["all_checks_pass"],
                        msg="all_checks_pass must be True at a=0")

    # Test 11
    def test_positive_spins_convergence_monotone(self) -> None:
        for s in self.spin_summaries:
            if s["spin_a"] > 0.0:
                self.assertTrue(
                    s["convergence_monotone_pass"],
                    msg=f"convergence_monotone_pass=False for a={s['spin_a']}",
                )

    # Test 12
    def test_positive_spins_omega_zamo_below_Omega_H(self) -> None:
        for s in self.spin_summaries:
            if s["spin_a"] > 0.0:
                self.assertTrue(
                    s["omega_zamo_below_Omega_H"],
                    msg=f"omega_zamo_below_Omega_H=False for a={s['spin_a']}",
                )

    # Test 13
    def test_global_causal_accounting_positive_spins(self) -> None:
        possible_pairs = self.n * (self.n - 1) // 2
        for s in self.spin_summaries:
            if s["spin_a"] > 0.0:
                self.assertEqual(
                    s["global_true_relations"], 0,
                    msg=f"global_true_relations != 0 for a={s['spin_a']}",
                )
                self.assertEqual(
                    s["global_false_relations"], 0,
                    msg=f"global_false_relations != 0 for a={s['spin_a']}",
                )
                self.assertEqual(
                    s["global_undecided_pairs"], possible_pairs,
                    msg=f"global_undecided_pairs != {possible_pairs} for a={s['spin_a']}",
                )

    def test_positive_spin_cases_all_undecided_aggregate(self) -> None:
        self.assertTrue(self.payload["aggregate"]["positive_spin_cases_all_undecided"])


class KerrK6ComputationTests(unittest.TestCase):
    """In-process tests: run the computation without relying on frozen artifacts."""

    def setUp(self) -> None:
        self.payload        = audit.run_audit()
        self.rows           = self.payload["rows"]
        self.spin_summaries = self.payload["spin_summaries"]
        self.aggregate      = self.payload["aggregate"]

    def test_all_checks_pass(self) -> None:
        self.assertTrue(self.aggregate["all_checks_pass"])

    def test_spin_list(self) -> None:
        self.assertEqual(self.aggregate["spins"], EXPECTED_SPINS)

    def test_delta_list(self) -> None:
        self.assertEqual(self.aggregate["deltas"], EXPECTED_DELTAS)

    def test_row_count(self) -> None:
        self.assertEqual(len(self.rows), EXPECTED_ROWS)

    def test_spin_summary_count(self) -> None:
        self.assertEqual(len(self.spin_summaries), len(EXPECTED_SPINS))

    def test_n_and_mass(self) -> None:
        self.assertEqual(self.aggregate["N"], 12)
        self.assertEqual(self.aggregate["M"], 1.0)

    def test_all_rows_exterior_pass(self) -> None:
        for row in self.rows:
            self.assertTrue(
                row["exterior_pass"],
                msg=f"a={row['spin_a']}, delta={row['delta']}",
            )

    def test_a0_checks_by_convention(self) -> None:
        a0 = next(s for s in self.spin_summaries if s["spin_a"] == 0.0)
        self.assertTrue(a0["convergence_monotone_pass"])
        self.assertTrue(a0["omega_zamo_below_Omega_H"])
        self.assertTrue(a0["all_checks_pass"])

    def test_positive_spins_convergence_monotone(self) -> None:
        for s in self.spin_summaries:
            if s["spin_a"] > 0.0:
                self.assertTrue(
                    s["convergence_monotone_pass"],
                    msg=f"a={s['spin_a']}",
                )

    def test_positive_spins_omega_zamo_below_Omega_H(self) -> None:
        for s in self.spin_summaries:
            if s["spin_a"] > 0.0:
                self.assertTrue(
                    s["omega_zamo_below_Omega_H"],
                    msg=f"a={s['spin_a']}",
                )

    def test_global_causal_accounting_positive_spins(self) -> None:
        n = self.aggregate["N"]
        possible = n * (n - 1) // 2
        for s in self.spin_summaries:
            if s["spin_a"] > 0.0:
                self.assertEqual(s["global_true_relations"],  0,       msg=f"a={s['spin_a']}")
                self.assertEqual(s["global_false_relations"], 0,       msg=f"a={s['spin_a']}")
                self.assertEqual(s["global_undecided_pairs"], possible, msg=f"a={s['spin_a']}")

    def test_positive_spin_cases_all_undecided(self) -> None:
        self.assertTrue(self.aggregate["positive_spin_cases_all_undecided"])

    # --- Unit tests on helpers ---

    def test_omega_H_analytic_a0(self) -> None:
        """Omega_H = 0 for a=0 (non-rotating)."""
        r_plus = 2.0  # Schwarzschild: r_+ = 2M
        self.assertEqual(audit.omega_H_analytic(r_plus, 0.0), 0.0)

    def test_omega_H_analytic_known_values(self) -> None:
        """Check a few analytic values."""
        M = 1.0
        for a in [0.25, 0.5, 0.75, 0.9]:
            r_plus = M + math.sqrt(M * M - a * a)
            expected = a / (r_plus ** 2 + a ** 2)
            result = audit.omega_H_analytic(r_plus, a)
            self.assertAlmostEqual(result, expected, places=13,
                                   msg=f"Omega_H mismatch for a={a}")

    def test_omega_zamo_a0_is_zero(self) -> None:
        """For a=0, omega_ZAMO = 0 everywhere outside the horizon."""
        M = 1.0
        for r in [2.5, 3.0, 4.0, 6.0, 10.0]:
            self.assertAlmostEqual(
                audit.omega_zamo_at_r(r, M, 0.0), 0.0, places=13,
                msg=f"omega_ZAMO != 0 at r={r} for a=0",
            )

    def test_omega_zamo_below_Omega_H_near_horizon(self) -> None:
        """For a>0 and r > r_+, omega_ZAMO < Omega_H."""
        M = 1.0
        for a in [0.25, 0.5, 0.75, 0.9]:
            r_plus = M + math.sqrt(M * M - a * a)
            Omega_H = audit.omega_H_analytic(r_plus, a)
            for delta in [1e-1, 3e-2, 1e-2, 3e-3, 1e-3]:
                r_eval = r_plus + delta
                omega_v = audit.omega_zamo_at_r(r_eval, M, a)
                self.assertLess(
                    omega_v, Omega_H,
                    msg=f"omega_ZAMO >= Omega_H for a={a}, delta={delta}",
                )

    def test_omega_zamo_converges_to_Omega_H(self) -> None:
        """Residuals decrease as delta -> 0."""
        M = 1.0
        deltas = [1e-1, 3e-2, 1e-2, 3e-3, 1e-3]
        for a in [0.25, 0.5, 0.75, 0.9]:
            r_plus = M + math.sqrt(M * M - a * a)
            Omega_H = audit.omega_H_analytic(r_plus, a)
            residuals = [
                abs(audit.omega_zamo_at_r(r_plus + d, M, a) - Omega_H)
                for d in deltas
            ]
            for i in range(len(residuals) - 1):
                self.assertGreater(
                    residuals[i], residuals[i + 1],
                    msg=f"residuals not monotone for a={a} at step {i}",
                )

    def test_omega_zamo_linear_convergence_rate(self) -> None:
        """residual/delta should be roughly constant (linear O(delta) convergence)."""
        M = 1.0
        deltas = [3e-2, 1e-2, 3e-3, 1e-3]   # skip largest delta for better linearity
        for a in [0.25, 0.5, 0.75]:
            r_plus = M + math.sqrt(M * M - a * a)
            Omega_H = audit.omega_H_analytic(r_plus, a)
            rates = [
                abs(audit.omega_zamo_at_r(r_plus + d, M, a) - Omega_H) / d
                for d in deltas
            ]
            # Rate should not vary by more than a factor of 2 across the range
            self.assertLess(max(rates) / min(rates), 2.0,
                            msg=f"linear rate not stable for a={a}: {rates}")

    def test_equatorial_metric_raises_inside_horizon(self) -> None:
        """equatorial_metric_at_r raises ValueError when Delta <= 0."""
        M = 1.0
        a = 0.5
        r_plus = M + math.sqrt(M * M - a * a)
        r_inside = r_plus - 0.1
        with self.assertRaises(ValueError):
            audit.equatorial_metric_at_r(r_inside, M, a)

    def test_equatorial_metric_schwarzschild_limit(self) -> None:
        """At a=0 the metric reduces to the Schwarzschild form."""
        M = 1.0
        TOL = 1.0e-12
        for r in [2.5, 3.0, 4.0, 6.0, 10.0]:
            m = audit.equatorial_metric_at_r(r, M, 0.0)
            self.assertAlmostEqual(m["g_tt"],   -(1.0 - 2.0 * M / r), delta=TOL, msg=f"g_tt at r={r}")
            self.assertAlmostEqual(m["g_tphi"],  0.0,                  delta=TOL, msg=f"g_tphi at r={r}")
            self.assertAlmostEqual(m["g_rr"],    1.0 / (1.0 - 2.0 * M / r), delta=TOL, msg=f"g_rr at r={r}")
            self.assertAlmostEqual(m["g_phiphi"], r * r,               delta=TOL, msg=f"g_phiphi at r={r}")

    def test_residuals_stored_in_spin_summaries(self) -> None:
        """Spin summaries carry residual lists of correct length."""
        for s in self.spin_summaries:
            self.assertEqual(len(s["residuals"]),        len(EXPECTED_DELTAS))
            self.assertEqual(len(s["omega_zamo_values"]), len(EXPECTED_DELTAS))
            self.assertEqual(len(s["deltas"]),            len(EXPECTED_DELTAS))


if __name__ == "__main__":
    unittest.main()
