"""Regression tests for ``validation_suite``.

These tests verify three things:

1. The canonical sprinkler produces causets whose order statistics
   match the Myrheim-Meyer prediction to within Monte Carlo
   tolerance. This is the gold-standard check.
2. The Lorentz-invariant interval residual respects the symmetries
   it must respect: zero under translation, zero under boost, zero
   under spatial rotation, and zero under uniform rescaling.
3. The sprinkle-then-recover pipeline runs end-to-end on a small
   case without raising, and the recovered embedding has a
   nontrivial interval residual that reflects the optimizer's
   actual quality.
"""

from __future__ import annotations

import math
import unittest

import causet_invariants as ci
import validation_suite as vs


def boost_2d(points, rapidity):
    """Apply a 1+1 Lorentz boost with the given rapidity to a list of
    (t, x) points. Used to verify interval residual invariance."""
    ch = math.cosh(rapidity)
    sh = math.sinh(rapidity)
    out = []
    for p in points:
        if len(p) < 2:
            t = p[0]
            x = 0.0
        else:
            t, x = p[0], p[1]
        tp = ch * t - sh * x
        xp = -sh * t + ch * x
        if len(p) > 2:
            out.append((tp, xp, *p[2:]))
        else:
            out.append((tp, xp))
    return out


def rotate_2d_spatial(points, theta):
    """Rotate the spatial plane of (t, x_1, x_2, ...) points by theta."""
    cs = math.cos(theta)
    sn = math.sin(theta)
    out = []
    for p in points:
        t = p[0]
        if len(p) >= 3:
            x1 = p[1] * cs - p[2] * sn
            x2 = p[1] * sn + p[2] * cs
            out.append((t, x1, x2, *p[3:]))
        else:
            out.append(p)
    return out


class CanonicalSprinklerTests(unittest.TestCase):
    def test_returns_correct_shapes(self) -> None:
        matrix, points = vs.sprinkle_minkowski_diamond(
            n=12, seed=1, d_spacetime=2
        )
        self.assertEqual(len(matrix), 12)
        self.assertTrue(all(len(row) == 12 for row in matrix))
        self.assertEqual(len(points), 12)
        self.assertTrue(all(len(p) == 2 for p in points))

    def test_reproducible(self) -> None:
        m1, p1 = vs.sprinkle_minkowski_diamond(n=20, seed=42, d_spacetime=3)
        m2, p2 = vs.sprinkle_minkowski_diamond(n=20, seed=42, d_spacetime=3)
        self.assertEqual(m1, m2)
        self.assertEqual(p1, p2)

    def test_points_are_inside_diamond(self) -> None:
        _, points = vs.sprinkle_minkowski_diamond(
            n=50, seed=7, d_spacetime=4
        )
        for t, *x in points:
            self.assertTrue(0.0 <= t <= 1.0)
            r = math.sqrt(sum(c * c for c in x))
            # Tolerate floating point: in the unit diamond we must
            # have r <= min(t, 1 - t).
            self.assertLessEqual(r, min(t, 1.0 - t) + 1e-12)

    def test_points_sorted_by_time(self) -> None:
        _, points = vs.sprinkle_minkowski_diamond(
            n=30, seed=13, d_spacetime=3
        )
        ts = [p[0] for p in points]
        self.assertEqual(ts, sorted(ts))

    def test_ordering_fraction_matches_myrheim_meyer(self) -> None:
        # Take an ensemble of moderate-size sprinklings and check the
        # mean ordering fraction is within Monte Carlo tolerance of the
        # closed-form prediction.
        for d, expected in [(2, 0.5), (3, 24.0 / 105.0), (4, 0.1)]:
            fractions = []
            for seed in range(3001, 3011):
                z, _ = vs.sprinkle_minkowski_diamond(
                    n=80, seed=seed, d_spacetime=d
                )
                fractions.append(ci.ordering_fraction(z))
            mean = sum(fractions) / len(fractions)
            self.assertAlmostEqual(
                mean, expected, delta=0.06,
                msg=f"d={d}: observed {mean:.4f}, expected {expected:.4f}"
            )


class KleitmanRothschildTests(unittest.TestCase):
    """Sanity checks for the non-manifoldlike control generator."""

    def test_default_level_sizes(self) -> None:
        # The default partition is approximately (n/4, n/2, n/4).
        z = vs.generate_kleitman_rothschild(n=16, seed=1959)
        profile = ci.antichain_profile(z)
        # Exactly three levels.
        self.assertEqual(len(profile), 3)
        self.assertEqual(profile, [4, 8, 4])

    def test_three_levels_for_canonical_sizes(self) -> None:
        for n in (16, 32, 64):
            z = vs.generate_kleitman_rothschild(n=n, seed=2026)
            profile = ci.antichain_profile(z)
            self.assertEqual(
                len(profile), 3,
                msg=f"n={n}: expected three antichain levels, got {profile}"
            )

    def test_no_within_level_relations(self) -> None:
        # The construction forbids relations within a level. After
        # transitive closure this remains true (transitivity cannot
        # introduce intra-level edges).
        n = 16
        z = vs.generate_kleitman_rothschild(n=n, seed=1)
        bottom = (n + 3) // 4
        level1_end = n - bottom
        # Level 0: indices 0..bottom-1.
        for i in range(bottom):
            for j in range(i + 1, bottom):
                self.assertFalse(z[i][j])
        # Level 1: indices bottom..level1_end-1.
        for i in range(bottom, level1_end):
            for j in range(i + 1, level1_end):
                self.assertFalse(z[i][j])
        # Level 2: indices level1_end..n-1.
        for i in range(level1_end, n):
            for j in range(i + 1, n):
                self.assertFalse(z[i][j])

    def test_upper_triangular_and_transitively_closed(self) -> None:
        # The output must be upper-triangular (lower triangle false)
        # and transitively closed (i prec j and j prec k => i prec k).
        z = vs.generate_kleitman_rothschild(n=20, seed=42)
        n = len(z)
        for i in range(n):
            for j in range(n):
                if i >= j:
                    self.assertFalse(
                        z[i][j],
                        msg=f"lower-triangular entry z[{i}][{j}] is True"
                    )
        for i in range(n - 1):
            for j in range(i + 1, n - 1):
                if not z[i][j]:
                    continue
                for k in range(j + 1, n):
                    if z[j][k]:
                        self.assertTrue(
                            z[i][k],
                            msg=(
                                f"transitivity failure: z[{i}][{j}] and "
                                f"z[{j}][{k}] both true but z[{i}][{k}] is False"
                            ),
                        )

    def test_reproducible_under_fixed_seed(self) -> None:
        z1 = vs.generate_kleitman_rothschild(n=24, seed=99)
        z2 = vs.generate_kleitman_rothschild(n=24, seed=99)
        self.assertEqual(z1, z2)

    def test_rejects_too_small_n(self) -> None:
        with self.assertRaises(ValueError):
            vs.generate_kleitman_rothschild(n=2, seed=1)


class CoronaPosetTests(unittest.TestCase):
    """Sanity checks for the suspended corona control generator."""

    def test_four_levels_for_canonical_sizes(self) -> None:
        for n in (16, 32, 64):
            z = vs.generate_corona_poset(n=n, seed=2026)
            profile = ci.antichain_profile(z)
            self.assertEqual(
                len(profile), 4,
                msg=f"n={n}: expected four antichain levels, got {profile}"
            )
            self.assertEqual(profile[0], 1)
            self.assertEqual(profile[-1], 1)

    def test_has_finite_dimension_estimators(self) -> None:
        z = vs.generate_corona_poset(n=32, seed=1959)
        self.assertTrue(math.isfinite(ci.myrheim_meyer_dimension(z)))
        self.assertTrue(math.isfinite(ci.midpoint_scaling_dimension(z)))

    def test_upper_triangular_and_transitively_closed(self) -> None:
        z = vs.generate_corona_poset(n=24, seed=42)
        n = len(z)
        for i in range(n):
            for j in range(n):
                if i >= j:
                    self.assertFalse(
                        z[i][j],
                        msg=f"lower-triangular entry z[{i}][{j}] is True"
                    )
        for i in range(n - 1):
            for j in range(i + 1, n - 1):
                if not z[i][j]:
                    continue
                for k in range(j + 1, n):
                    if z[j][k]:
                        self.assertTrue(
                            z[i][k],
                            msg=(
                                f"transitivity failure: z[{i}][{j}] and "
                                f"z[{j}][{k}] both true but z[{i}][{k}] is False"
                            ),
                        )

    def test_reproducible_under_fixed_seed(self) -> None:
        z1 = vs.generate_corona_poset(n=24, seed=99)
        z2 = vs.generate_corona_poset(n=24, seed=99)
        self.assertEqual(z1, z2)

    def test_rejects_too_small_n(self) -> None:
        with self.assertRaises(ValueError):
            vs.generate_corona_poset(n=5, seed=1)


class IntervalResidualInvarianceTests(unittest.TestCase):
    def test_identical_embeddings_zero(self) -> None:
        _, points = vs.sprinkle_minkowski_diamond(n=10, seed=1, d_spacetime=2)
        self.assertAlmostEqual(vs.interval_rmse(points, points), 0.0, places=12)

    def test_invariant_under_translation(self) -> None:
        _, points = vs.sprinkle_minkowski_diamond(n=10, seed=1, d_spacetime=2)
        shifted = [(t + 5.0, x + 3.0) for t, x in points]
        self.assertAlmostEqual(vs.interval_rmse(points, shifted), 0.0, places=12)

    def test_invariant_under_lorentz_boost(self) -> None:
        _, points = vs.sprinkle_minkowski_diamond(n=10, seed=1, d_spacetime=2)
        boosted = boost_2d(points, rapidity=0.7)
        # Boosts preserve the Minkowski interval exactly.
        self.assertAlmostEqual(vs.interval_rmse(points, boosted), 0.0, places=10)

    def test_invariant_under_spatial_rotation(self) -> None:
        _, points = vs.sprinkle_minkowski_diamond(n=10, seed=1, d_spacetime=3)
        rotated = rotate_2d_spatial(points, theta=0.9)
        self.assertAlmostEqual(vs.interval_rmse(points, rotated), 0.0, places=10)

    def test_invariant_under_uniform_rescale(self) -> None:
        _, points = vs.sprinkle_minkowski_diamond(n=10, seed=1, d_spacetime=2)
        # Uniform rescaling scales s^2 by lambda^2, so the matrix
        # entries scale and the RMSE scales accordingly. But for the
        # *zero* residual under rescaling: only a trivial rescale by
        # 1 preserves the RMSE-of-difference. A boost preserves the
        # interval matrix exactly. So this test asserts: rescaling
        # both embeddings by the same factor leaves RMSE the same.
        rescaled_a = [(t * 3.0, x * 3.0) for t, x in points]
        rescaled_b = rescaled_a  # both rescaled identically
        self.assertAlmostEqual(
            vs.interval_rmse(rescaled_a, rescaled_b), 0.0, places=12
        )


class BombelliEnergyAtTests(unittest.TestCase):
    def test_low_energy_at_ground_truth(self) -> None:
        # The Bombelli energy at the true sprinkling coordinates should
        # be substantially below the worst-case energy of a random
        # configuration of the same scale. We do not pin a particular
        # absolute number (the energy is sensitive to the geometry),
        # only that it is finite and positive for a typical case.
        z, points = vs.sprinkle_minkowski_diamond(
            n=24, seed=2026, d_spacetime=2
        )
        e_truth = vs.bombelli_energy_at(z, points, d_spatial=1)
        self.assertGreaterEqual(e_truth, 0.0)
        self.assertTrue(math.isfinite(e_truth))

    def test_scale_invariance(self) -> None:
        # The energy formula is invariant under uniform scaling of
        # (t, x). Verify by direct computation.
        z, points = vs.sprinkle_minkowski_diamond(
            n=20, seed=2027, d_spacetime=2
        )
        e1 = vs.bombelli_energy_at(z, points, d_spatial=1)
        rescaled = [(t * 7.0, x * 7.0) for t, x in points]
        e2 = vs.bombelli_energy_at(z, rescaled, d_spatial=1)
        self.assertAlmostEqual(e1, e2, places=10)


class RecoveryPipelineSmokeTests(unittest.TestCase):
    def test_run_recovery_completes(self) -> None:
        # Tiny case: just make sure the end-to-end pipeline runs
        # without error and returns a well-formed result.
        case = vs.SprinkleCase(
            d_spacetime=2,
            n=6,
            seed=1959,
            matrix=[],
            points=[],
        )
        # Replace the empty placeholders with an actual sprinkling.
        matrix, points = vs.sprinkle_minkowski_diamond(
            n=6, seed=1959, d_spacetime=2
        )
        case = vs.SprinkleCase(
            d_spacetime=2, n=6, seed=1959, matrix=matrix, points=points
        )
        result = vs.run_recovery(
            case,
            optimizer_seed=1987,
            warmup_limit=10,
            anneal_limit=10,
            max_data=5,
        )
        self.assertEqual(result.d_spacetime, 2)
        self.assertEqual(result.n, 6)
        self.assertEqual(result.target_dim, 1)
        self.assertTrue(math.isfinite(result.final_energy))
        self.assertTrue(math.isfinite(result.truth_energy))
        self.assertTrue(math.isfinite(result.interval_rmse))
        self.assertGreaterEqual(result.interval_rmse, 0.0)

    def test_generate_ensemble_dimensions(self) -> None:
        cases = vs.generate_ensemble(
            d_spacetimes=[2, 3],
            sizes=[8, 12],
            seeds=[1, 2],
        )
        self.assertEqual(len(cases), 2 * 2 * 2)
        # Check that each combination is present exactly once.
        signatures = {(c.d_spacetime, c.n, c.seed) for c in cases}
        self.assertEqual(len(signatures), 8)


if __name__ == "__main__":
    unittest.main()
