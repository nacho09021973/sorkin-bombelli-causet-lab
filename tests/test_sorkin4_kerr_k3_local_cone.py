"""Tests for S4-KERR-K3-LOCAL-CONE-001.

Covers:
1.  Artifact files exist and are parseable.
2.  all_checks_pass is True in the JSON aggregate.
3.  a=0 schwarzschild_reduction_pass is True.
4.  For a=0, max_abs_g_tphi == 0 (within tolerance).
5.  For a>0, frame_dragging_sign_pass is True.
6.  For all spins, all_points_exterior is True.
7.  For all spins, min_Delta > 0, min_g_rr > 0, min_g_phiphi > 0.
8.  For a>0, global_true_relations=0, global_false_relations=0,
    global_undecided_pairs = N*(N-1)/2.
9.  local_timelike_count + local_nullish_count + local_spacelike_count
    == local_evaluated_pair_count (NOT necessarily N*(N-1)/2).
10. Existing K1/K2 tests: covered by the full validation command; not
    duplicated here to avoid redundancy.
"""
from __future__ import annotations

import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "explore" / "sorkin4_kerr_benchmark"))

import audit_kerr_k3_local_cone_001 as audit

ARTIFACT_DIR = ROOT / "explore" / "sorkin4_kerr_benchmark"
ARTIFACT_JSON = ARTIFACT_DIR / "kerr_k3_local_cone_001_n12_seed1959.json"
ARTIFACT_CSV  = ARTIFACT_DIR / "kerr_k3_local_cone_001_n12_seed1959.csv"
ARTIFACT_MD   = ARTIFACT_DIR / "kerr_k3_local_cone_001_n12_seed1959.md"

METRIC_TOL = 1.0e-12


class KerrK3ArtifactTests(unittest.TestCase):
    """Guard the committed frozen artifacts (Tests 1-9 via artifact)."""

    def setUp(self) -> None:
        for path in (ARTIFACT_JSON, ARTIFACT_CSV, ARTIFACT_MD):
            self.assertTrue(
                path.exists(),
                msg=(
                    f"Missing K3 artifact {path}. Run: "
                    "python explore/sorkin4_kerr_benchmark/"
                    "audit_kerr_k3_local_cone_001.py"
                ),
            )
        with ARTIFACT_JSON.open(encoding="utf-8") as fh:
            self.payload = json.load(fh)

    # Test 1: files exist & parseable — covered by setUp above

    # Test 2: aggregate all_checks_pass
    def test_all_checks_pass(self) -> None:
        self.assertTrue(self.payload["aggregate"]["all_checks_pass"])

    # Tests 3-9 on frozen rows
    def test_per_row_checks(self) -> None:
        rows = self.payload["rows"]
        n = self.payload["aggregate"]["N"]
        possible_pairs = n * (n - 1) // 2

        for row in rows:
            spin = row["spin_a"]

            # Test 6: all_points_exterior
            self.assertTrue(
                row["all_points_exterior"],
                msg=f"all_points_exterior failed for a={spin}",
            )
            # Test 7: metric positivity
            self.assertGreater(row["min_Delta"],    0.0, msg=f"min_Delta ≤ 0 for a={spin}")
            self.assertGreater(row["min_g_rr"],     0.0, msg=f"min_g_rr ≤ 0 for a={spin}")
            self.assertGreater(row["min_g_phiphi"], 0.0, msg=f"min_g_phiphi ≤ 0 for a={spin}")

            # Test 9: local count integrity
            local_sum = (
                row["local_timelike_count"]
                + row["local_nullish_count"]
                + row["local_spacelike_count"]
            )
            self.assertEqual(
                local_sum,
                row["local_evaluated_pair_count"],
                msg=(
                    f"local count mismatch for a={spin}: "
                    f"{local_sum} != {row['local_evaluated_pair_count']}"
                ),
            )

            if abs(spin) <= 0.0:
                # Test 3: Schwarzschild reduction
                self.assertTrue(
                    row["schwarzschild_reduction_pass"],
                    msg="schwarzschild_reduction_pass failed for a=0",
                )
                # Test 4: g_tphi = 0 for a=0
                self.assertLessEqual(
                    row["max_abs_g_tphi"],
                    METRIC_TOL,
                    msg=f"max_abs_g_tphi={row['max_abs_g_tphi']} > tol for a=0",
                )
            else:
                # Test 5: frame dragging sign for a>0
                self.assertTrue(
                    row["frame_dragging_sign_pass"],
                    msg=f"frame_dragging_sign_pass failed for a={spin}",
                )
                # Test 8: global causal accounting for a>0
                self.assertEqual(row["global_true_relations"],  0, msg=f"a={spin}")
                self.assertEqual(row["global_false_relations"], 0, msg=f"a={spin}")
                self.assertEqual(
                    row["global_undecided_pairs"],
                    possible_pairs,
                    msg=f"a={spin}",
                )


class KerrK3ComputationTests(unittest.TestCase):
    """In-process tests: run computation without relying on artifact files."""

    def setUp(self) -> None:
        self.payload   = audit.run_audit()
        self.rows      = self.payload["rows"]
        self.aggregate = self.payload["aggregate"]

    def test_aggregate_all_checks_pass(self) -> None:
        self.assertTrue(self.aggregate["all_checks_pass"])

    def test_spins_and_n(self) -> None:
        self.assertEqual(self.aggregate["spins"], [0.0, 0.25, 0.5, 0.75])
        self.assertEqual(self.aggregate["M"], 1.0)
        self.assertEqual(self.aggregate["N"], 12)
        self.assertEqual(len(self.rows), 4)

    def test_per_row_local_count_integrity(self) -> None:
        for row in self.rows:
            local_sum = (
                row["local_timelike_count"]
                + row["local_nullish_count"]
                + row["local_spacelike_count"]
            )
            self.assertEqual(
                local_sum,
                row["local_evaluated_pair_count"],
                msg=f"local count mismatch for a={row['spin_a']}",
            )

    def test_a0_schwarzschild_reduction(self) -> None:
        a0_row = next(r for r in self.rows if r["spin_a"] == 0.0)
        self.assertTrue(a0_row["schwarzschild_reduction_pass"])
        self.assertLessEqual(a0_row["max_abs_g_tphi"], METRIC_TOL)

    def test_positive_spin_frame_dragging(self) -> None:
        for row in self.rows:
            if row["spin_a"] > 0.0:
                self.assertTrue(
                    row["frame_dragging_sign_pass"],
                    msg=f"frame_dragging_sign_pass failed for a={row['spin_a']}",
                )

    def test_positive_spin_global_undecided(self) -> None:
        n = self.aggregate["N"]
        possible = n * (n - 1) // 2
        for row in self.rows:
            if row["spin_a"] > 0.0:
                self.assertEqual(row["global_true_relations"],  0, msg=f"a={row['spin_a']}")
                self.assertEqual(row["global_false_relations"], 0, msg=f"a={row['spin_a']}")
                self.assertEqual(row["global_undecided_pairs"], possible, msg=f"a={row['spin_a']}")

    def test_equatorial_metric_at_r_a0_schwarzschild(self) -> None:
        """a=0 must give exact Schwarzschild equatorial BL coefficients."""
        r, M = 4.0, 1.0
        m = audit.equatorial_metric_at_r(r, M, 0.0)
        self.assertAlmostEqual(m["g_tt"],     -(1.0 - 2 * M / r),  places=14)
        self.assertAlmostEqual(m["g_tphi"],   0.0,                  places=14)
        self.assertAlmostEqual(m["g_rr"],     1.0 / (1.0 - 2*M/r), places=14)
        self.assertAlmostEqual(m["g_phiphi"], r * r,                places=14)

    def test_equatorial_metric_at_r_frame_dragging_sign(self) -> None:
        """g_tphi < 0 for a>0, M>0 at any exterior r."""
        for a in (0.25, 0.5, 0.75):
            m = audit.equatorial_metric_at_r(4.0, 1.0, a)
            self.assertLess(m["g_tphi"], 0.0, msg=f"g_tphi >= 0 for a={a}")

    def test_local_pair_filter(self) -> None:
        """Only pairs with |dr|<=DR and |dphi|<=DPHI are evaluated."""
        from run_kerr_minimal_benchmark import Event
        p = Event(index=0, t=0.0, r=3.0, theta=math.pi / 2, phi=0.0)
        # Within threshold
        q_near = Event(index=1, t=1.0, r=3.2, theta=math.pi / 2, phi=0.1)
        self.assertTrue(audit.is_local_pair(p, q_near))
        # |dr| = 2.5 > DR_LOCAL_THRESHOLD = 1.0
        q_far_r = Event(index=2, t=1.0, r=5.5, theta=math.pi / 2, phi=0.1)
        self.assertFalse(audit.is_local_pair(p, q_far_r))
        # |dphi| = 0.6 > DPHI_LOCAL_THRESHOLD = 0.5
        q_far_phi = Event(index=3, t=1.0, r=3.1, theta=math.pi / 2, phi=0.6)
        self.assertFalse(audit.is_local_pair(p, q_far_phi))


if __name__ == "__main__":
    unittest.main()
