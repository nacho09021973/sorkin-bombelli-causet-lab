"""Regression tests for ``causet_invariants``.

The order-theoretic invariants are deliberately deterministic. Each
test fixes a concrete causet and checks the invariants against values
that can be derived by hand or from the closed-form Myrheim-Meyer
formula.

For the Poisson-sprinkled cases we use moderately large ``n`` so that
the empirical ordering fraction sits within a tolerance of the
expected continuum value. The tolerances are deliberately loose:
the point of the test is to detect implementation bugs, not to
characterize statistical fluctuations.
"""

from __future__ import annotations

import math
import unittest
from pathlib import Path

import causet_invariants as ci
import cones


ROOT = Path(__file__).resolve().parents[1]


def make_chain(n: int) -> list[list[bool]]:
    """Total order on n elements: i prec j for every i < j."""

    return [[(i < j) for j in range(n)] for i in range(n)]


def make_antichain(n: int) -> list[list[bool]]:
    """No relations at all: a width-n antichain."""

    return [[False] * n for _ in range(n)]


def make_two_chains(k: int) -> list[list[bool]]:
    """Two parallel chains of length k, no cross relations.

    Labeling interleaves the chains: events 0, 2, 4, ... form one chain
    and 1, 3, 5, ... form the other. Each pair of same-parity indices
    is related; different-parity pairs are not.
    """

    n = 2 * k
    z = [[False] * n for _ in range(n)]
    for i in range(n - 1):
        for j in range(i + 1, n):
            if (i % 2) == (j % 2):
                z[i][j] = True
    return z


class OrderingFractionTests(unittest.TestCase):
    def test_total_chain(self) -> None:
        z = make_chain(5)
        self.assertEqual(ci.relation_count(z), 10)
        self.assertEqual(ci.ordering_fraction(z), 1.0)

    def test_antichain(self) -> None:
        z = make_antichain(7)
        self.assertEqual(ci.relation_count(z), 0)
        self.assertEqual(ci.ordering_fraction(z), 0.0)

    def test_two_disjoint_chains(self) -> None:
        # Two chains of length 4 means 6 + 6 = 12 related pairs;
        # total pairs = C(8, 2) = 28; fraction = 12/28.
        z = make_two_chains(4)
        self.assertEqual(ci.relation_count(z), 12)
        self.assertAlmostEqual(ci.ordering_fraction(z), 12.0 / 28.0)


class MyrheimMeyerTests(unittest.TestCase):
    def test_known_values(self) -> None:
        # Closed-form values for the unordered-pair convention.
        # f(1) = 1: every pair related in 1D Minkowski (just time).
        self.assertAlmostEqual(ci._myrheim_meyer_f(1.0), 1.0, places=12)
        # f(2) = 1/2: 1+1 unit diamond, verified empirically.
        self.assertAlmostEqual(ci._myrheim_meyer_f(2.0), 0.5, places=12)
        # f(3) = 24/105 ~= 0.2286.
        self.assertAlmostEqual(ci._myrheim_meyer_f(3.0), 24.0 / 105.0, places=12)
        # f(4) = 1/10.
        self.assertAlmostEqual(ci._myrheim_meyer_f(4.0), 0.1, places=12)

    def test_chain_returns_d_min(self) -> None:
        # A total chain has ordering fraction 1.0 = f(1). The estimator
        # clamps to d_min (1.0) at this saturating value.
        z = make_chain(6)
        self.assertEqual(ci.myrheim_meyer_dimension(z), 1.0)

    def test_antichain_returns_infinity(self) -> None:
        z = make_antichain(6)
        self.assertEqual(ci.myrheim_meyer_dimension(z), float("inf"))

    def test_recovers_dimension_on_1plus1_sprinkling(self) -> None:
        # 1+1 Minkowski has expected ordering fraction f(2) = 1/2.
        # cones.generate_sprinkled_causet with spacetime_dim=1 IS a
        # canonical uniform sprinkling of the unit Minkowski diamond
        # (the (u, v) light-cone square maps to the diamond with
        # Jacobian 2, giving uniform measure in (t, x)). Sprinkle a
        # single causet of 80 events and check the estimator gets
        # close to d = 2.
        z, _ = cones.generate_sprinkled_causet(80, seed=1987, spacetime_dim=1)
        d_est = ci.myrheim_meyer_dimension(z)
        self.assertAlmostEqual(d_est, 2.0, delta=0.4)

    def test_recovers_dimension_on_canonical_diamond_ensemble(self) -> None:
        # Average over an ensemble of canonical sprinklings in 2+1 and
        # 3+1 and check that the recovered dimension is within 0.4 of
        # the truth. The bound is loose: at n = 80 the per-seed
        # standard deviation in d_est is ~0.2; with 8 seeds the
        # standard error of the mean is ~0.07, so a tolerance of 0.4
        # is comfortable.
        #
        # Uses validation_suite.sprinkle_minkowski_diamond, which is
        # the canonical (uniform-in-diamond) sprinkler. The legacy
        # cones.generate_sprinkled_causet does not sample uniformly
        # for d_spatial >= 2 (its marginal in t is uniform, not
        # proportional to the cross-section volume), so it cannot be
        # used here.
        import validation_suite as vs
        for d_spacetime, expected in [(3, 3.0), (4, 4.0)]:
            est = []
            for seed in range(2001, 2009):
                z, _ = vs.sprinkle_minkowski_diamond(
                    n=80, seed=seed, d_spacetime=d_spacetime
                )
                est.append(ci.myrheim_meyer_dimension(z))
            mean_est = sum(est) / len(est)
            self.assertAlmostEqual(
                mean_est, expected, delta=0.4,
                msg=f"d_spacetime={d_spacetime}: recovered {mean_est:.3f}"
            )


class HeightAndProfileTests(unittest.TestCase):
    def test_chain_has_height_n(self) -> None:
        z = make_chain(5)
        self.assertEqual(ci.height(z), 5)
        self.assertEqual(ci.antichain_profile(z), [1, 1, 1, 1, 1])

    def test_antichain_has_height_one(self) -> None:
        z = make_antichain(5)
        self.assertEqual(ci.height(z), 1)
        self.assertEqual(ci.antichain_profile(z), [5])

    def test_two_chains_have_height_k(self) -> None:
        z = make_two_chains(4)
        self.assertEqual(ci.height(z), 4)
        # Each level should contain exactly 2 elements (one from each chain).
        self.assertEqual(ci.antichain_profile(z), [2, 2, 2, 2])


class ChainCountsTests(unittest.TestCase):
    def test_chain_binomial(self) -> None:
        # For a total chain of length n, the number of k-chains is C(n, k).
        n = 6
        z = make_chain(n)
        counts = ci.chain_counts(z, k_max=4)
        self.assertEqual(counts[2], math.comb(n, 2))
        self.assertEqual(counts[3], math.comb(n, 3))
        self.assertEqual(counts[4], math.comb(n, 4))

    def test_antichain_has_no_chains(self) -> None:
        z = make_antichain(6)
        counts = ci.chain_counts(z, k_max=4)
        self.assertEqual(counts[2], 0)
        self.assertEqual(counts[3], 0)
        self.assertEqual(counts[4], 0)

    def test_three_chain_abundance_normalization(self) -> None:
        z = make_chain(6)
        self.assertAlmostEqual(ci.three_chain_abundance(z), 1.0)
        self.assertEqual(ci.three_chain_abundance(make_antichain(6)), 0.0)

    def test_three_chain_abundance_finite_on_phase1d_families(self) -> None:
        import validation_suite as vs

        matrices = [
            vs.sprinkle_minkowski_diamond(n=32, seed=1959, d_spacetime=3)[0],
            vs.generate_kleitman_rothschild(n=32, seed=1959),
            vs.generate_corona_poset(n=32, seed=1959),
        ]
        for z in matrices:
            abundance = ci.three_chain_abundance(z)
            self.assertTrue(math.isfinite(abundance))
            self.assertGreaterEqual(abundance, 0.0)
            self.assertLessEqual(abundance, 1.0)


class LinkCountTests(unittest.TestCase):
    def test_chain_has_n_minus_one_links(self) -> None:
        # In a total chain only consecutive elements are linked.
        z = make_chain(5)
        self.assertEqual(ci.link_count(z), 4)

    def test_antichain_has_no_links(self) -> None:
        z = make_antichain(5)
        self.assertEqual(ci.link_count(z), 0)

    def test_against_tesis_like_6(self) -> None:
        z = cones.parse_cones_input(ROOT / "benchmarks" / "tesis_like_6.in")
        # The benchmark records 4 relations; for this small causet the
        # transitive reduction equals the full relation graph (no
        # triangle could shorten anything beyond what is recorded).
        self.assertEqual(ci.relation_count(z), 4)
        self.assertEqual(ci.link_count(z), 4)


class MidpointScalingTests(unittest.TestCase):
    """Sanity checks for Meyer's midpoint interval scaler."""

    def test_chain_dimension_is_one(self) -> None:
        # A total chain has a unique midpoint that bisects the chain.
        # The ratio of full to half cardinality is approximately 2,
        # giving d ~ log_2(2) = 1.
        z = make_chain(8)
        d = ci.midpoint_scaling_dimension(z)
        self.assertAlmostEqual(d, 1.0, delta=0.6)

    def test_antichain_returns_nan(self) -> None:
        z = make_antichain(8)
        d = ci.midpoint_scaling_dimension(z)
        self.assertTrue(math.isnan(d))

    def test_tiny_causet_returns_nan(self) -> None:
        # A 2-element chain has only the trivial endpoint interval.
        z = make_chain(2)
        d = ci.midpoint_scaling_dimension(z)
        self.assertTrue(math.isnan(d))

    def test_finite_on_sprinkling(self) -> None:
        # On a moderate Minkowski sprinkling the estimator must
        # produce a finite, real number. The numerical agreement
        # with the true dimension is covered separately by the
        # Phase 1 atlas; this test only guards against blow-ups.
        z, _ = cones.generate_sprinkled_causet(
            n=64, seed=1987, spacetime_dim=1
        )
        d = ci.midpoint_scaling_dimension(z)
        self.assertTrue(math.isfinite(d))
        self.assertGreater(d, 0.0)

    def test_finite_on_kleitman_rothschild(self) -> None:
        # The KR generator lives in validation_suite; importing it
        # here exercises the cross-module diagnostic path the Phase 1
        # atlas relies on.
        import validation_suite as vs
        z = vs.generate_kleitman_rothschild(n=64, seed=2026)
        d = ci.midpoint_scaling_dimension(z)
        self.assertTrue(math.isfinite(d))
        self.assertGreater(d, 0.0)


class FingerprintTests(unittest.TestCase):
    def test_deterministic_on_fixed_causet(self) -> None:
        z = cones.parse_cones_input(ROOT / "benchmarks" / "tesis_like_12.in")
        fp1 = ci.invariants_fingerprint(z)
        fp2 = ci.invariants_fingerprint(z)
        self.assertEqual(fp1, fp2)
        # Sanity: the recorded benchmark has 19 relations.
        self.assertEqual(fp1["relation_count"], 19)
        self.assertEqual(fp1["n"], 12)


if __name__ == "__main__":
    unittest.main()
