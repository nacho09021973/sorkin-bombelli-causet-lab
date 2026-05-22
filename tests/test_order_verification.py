"""Unit tests for the SORKIN-2 order verification infrastructure.

Tests ``induced_order_from_coords``, ``compare_causal_orders``, and
``verify_recovery`` from ``validation_suite``.  All tests are purely
order-theoretic or geometric: no Bombelli energy is evaluated, no
dimension estimator is invoked, and no benchmark data are required.

The ``TestVerifyRecovery`` class runs the simulator once on
``tesis_like_6.in`` with a short schedule.  It checks only the return
type and field names — not whether the annealer converged — so it
makes no claim about recoverability.
"""

from __future__ import annotations

import contextlib
import io
import math
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import cones
import validation_suite as vs


def make_chain(n: int) -> list[list[bool]]:
    return [[(i < j) for j in range(n)] for i in range(n)]


def make_antichain(n: int) -> list[list[bool]]:
    return [[False] * n for _ in range(n)]


class TestInducedOrderFromCoords(unittest.TestCase):
    """Geometric causal order computed from explicit coordinates."""

    def test_chain_4_timelike_column(self) -> None:
        # Four events on the r-axis: r increases with label index.
        # Cones convention: for z[i][j], zero energy requires r[i] <= r[j].
        # With r = [1, 1.5, 2, 2.5] and all x=0, every i<j pair is causal.
        coords = [(1.0,), (1.5,), (2.0,), (2.5,)]
        z = vs.induced_order_from_coords(coords)
        expected = make_chain(4)
        self.assertEqual(z, expected)

    def test_chain_4_exact_run_coordinates(self) -> None:
        coords = [
            (1.0, 0.0, 0.0),
            (1.5, 0.0, 0.0),
            (2.0, 0.0, 0.0),
            (2.5, 0.0, 0.0),
        ]
        z = vs.induced_order_from_coords(coords)
        expected = make_chain(4)
        self.assertEqual(z, expected)
        expected_pairs = {(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)}
        actual_pairs = {
            (i, j)
            for i in range(4)
            for j in range(i + 1, 4)
            if z[i][j]
        }
        self.assertEqual(actual_pairs, expected_pairs)

    def test_antichain_4_spacelike_column(self) -> None:
        # Equal-r pairs are null/timelike under the cones convention.
        # Use spacelike separation so no pair is causally related.
        coords = [(1.0, 0.0), (1.0, 2.0), (1.0, 4.0), (1.0, 6.0)]
        z = vs.induced_order_from_coords(coords)
        expected = make_antichain(4)
        self.assertEqual(z, expected)

    def test_equal_r_null_column_is_chain(self) -> None:
        coords = [(1.0,), (1.0,), (1.0,), (1.0,)]
        z = vs.induced_order_from_coords(coords)
        expected = make_chain(4)
        self.assertEqual(z, expected)

    def test_timelike_pair_is_related(self) -> None:
        # (r=1, x=0) precedes (r=2, x=0): rij=-1, dx_sq=0, 1 >= 0.
        coords = [(1.0, 0.0), (2.0, 0.0)]
        z = vs.induced_order_from_coords(coords)
        self.assertTrue(z[0][1])

    def test_spacelike_pair_not_related(self) -> None:
        # rij=-1, |dx|=2: 1^2 < 2^2, so not causal.
        coords = [(1.0, 0.0), (2.0, 2.0)]
        z = vs.induced_order_from_coords(coords)
        self.assertFalse(z[0][1])

    def test_null_pair_is_related(self) -> None:
        # rij=-1, |dx|=1: 1^2 >= 1^2, on the light cone -> causal.
        coords = [(1.0, 0.0), (2.0, 1.0)]
        z = vs.induced_order_from_coords(coords)
        self.assertTrue(z[0][1])

    def test_reversed_r_not_related(self) -> None:
        # r[0] > r[1]: event 0 cannot precede event 1 in the i<j matrix.
        coords = [(2.0, 0.0), (1.0, 0.0)]
        z = vs.induced_order_from_coords(coords)
        self.assertFalse(z[0][1])

    def test_decreasing_r_column_induces_no_upper_triangular_relations(self) -> None:
        coords = [(4.0, 0.0), (3.0, 0.0), (2.0, 0.0), (1.0, 0.0)]
        z = vs.induced_order_from_coords(coords)
        expected = make_antichain(4)
        self.assertEqual(z, expected)

    def test_returns_upper_triangular(self) -> None:
        # Lower triangle must be all False (only i<j pairs are set).
        coords = [(1.0, 0.0), (2.0, 0.0), (3.0, 0.0)]
        z = vs.induced_order_from_coords(coords)
        for i in range(3):
            for j in range(i):
                self.assertFalse(z[i][j], msg=f"lower triangle z[{i}][{j}] should be False")


class TestCompareCausalOrders(unittest.TestCase):
    """Element-by-element comparison of two order matrices."""

    def test_identical_chain_is_exact_match(self) -> None:
        z = make_chain(4)
        result = vs.compare_causal_orders(z, z)
        self.assertTrue(result.exact_match)
        self.assertEqual(result.missing_relations, ())
        self.assertEqual(result.extra_relations, ())

    def test_identical_antichain_is_exact_match(self) -> None:
        z = make_antichain(4)
        result = vs.compare_causal_orders(z, z)
        self.assertTrue(result.exact_match)
        self.assertEqual(result.missing_relations, ())
        self.assertEqual(result.extra_relations, ())

    def test_missing_one_relation(self) -> None:
        target = make_chain(4)
        # Induced order is missing the (0,1) relation.
        induced = [row[:] for row in target]
        induced[0][1] = False
        result = vs.compare_causal_orders(target, induced)
        self.assertFalse(result.exact_match)
        self.assertIn((0, 1), result.missing_relations)
        self.assertEqual(result.extra_relations, ())

    def test_extra_one_relation(self) -> None:
        target = make_antichain(4)
        # Induced order has an extra (1, 2) relation not in target.
        induced = [row[:] for row in target]
        induced[1][2] = True
        result = vs.compare_causal_orders(target, induced)
        self.assertFalse(result.exact_match)
        self.assertIn((1, 2), result.extra_relations)
        self.assertEqual(result.missing_relations, ())

    def test_n_field(self) -> None:
        z = make_chain(5)
        result = vs.compare_causal_orders(z, z)
        self.assertEqual(result.n, 5)

    def test_relation_counts(self) -> None:
        z = make_chain(4)  # 4*3/2 = 6 pairs all True
        result = vs.compare_causal_orders(z, z)
        self.assertEqual(result.total_relations_target, 6)
        self.assertEqual(result.total_relations_induced, 6)

    def test_size_mismatch_raises(self) -> None:
        with self.assertRaises(ValueError):
            vs.compare_causal_orders(make_chain(3), make_chain(4))


class TestVerifyRecovery(unittest.TestCase):
    """verify_recovery returns an OrderComparison from a live simulator.

    This test makes no assertion about exact_match: the historical
    annealer may or may not converge on a short schedule.  The test
    verifies only that the function runs, returns the correct type, and
    populates all fields consistently.
    """

    @classmethod
    def setUpClass(cls) -> None:
        z = cones.parse_cones_input(ROOT / "benchmarks" / "tesis_like_6.in")
        sim = cones.ConesSimulator(
            z=z,
            dim=2,
            seed=1959,
            interactive=False,
            max_data=35,
            plot_path=None,
            warmup_limit=10,
            anneal_limit=10,
            initial_temp=100.0,
            cooling_factor=0.9,
            acceptance_scale=4.0,
            backend="cpu",
        )
        with tempfile.TemporaryDirectory() as tmpdir, contextlib.redirect_stdout(io.StringIO()):
            sim.run(Path(tmpdir) / "out.txt")
        cls.sim = sim
        cls.result = vs.verify_recovery(sim)

    def test_returns_order_comparison(self) -> None:
        self.assertIsInstance(self.result, vs.OrderComparison)

    def test_n_matches_simulator(self) -> None:
        self.assertEqual(self.result.n, self.sim.n)

    def test_missing_and_extra_are_tuples(self) -> None:
        self.assertIsInstance(self.result.missing_relations, tuple)
        self.assertIsInstance(self.result.extra_relations, tuple)

    def test_exact_match_is_bool(self) -> None:
        self.assertIsInstance(self.result.exact_match, bool)

    def test_exact_match_consistent_with_empty_sets(self) -> None:
        expected = (
            len(self.result.missing_relations) == 0
            and len(self.result.extra_relations) == 0
        )
        self.assertEqual(self.result.exact_match, expected)

    def test_relation_counts_non_negative(self) -> None:
        self.assertGreaterEqual(self.result.total_relations_target, 0)
        self.assertGreaterEqual(self.result.total_relations_induced, 0)


if __name__ == "__main__":
    unittest.main()
