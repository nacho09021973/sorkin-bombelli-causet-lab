"""Invariance tests for the diagnostic foundation.

These check the structural symmetries that any honest diagnostic
must respect:

- **Label permutation**: re-labeling the events of a causet (while
  keeping the abstract poset the same) must not change any
  order-theoretic invariant. This catches bugs where an algorithm
  accidentally depends on the input ordering instead of on the
  partial order itself.

- **Lorentz invariance** of the interval RMSE: already exercised
  in ``test_validation_suite``; this file documents the property
  in one canonical place by exercising a composite transformation
  (translation + boost + rotation + uniform rescale).

The existing optimizer-seed reproducibility is already covered by
``tests/test_regression.py::test_thesis_like_*_regression``; we do
not duplicate those checks here.
"""

from __future__ import annotations

import math
import random
import unittest
from typing import List

import causet_invariants
import validation_suite as vs


def _topologically_sort_matrix(
    matrix: List[List[bool]],
    permutation: List[int],
) -> List[List[bool]]:
    """Apply a permutation to the labels of a causet, then re-canonicalize.

    The project's convention is to label events in a (any)
    topologically consistent order: ``i prec j`` implies ``i < j``.
    After applying a permutation, the labels may no longer respect
    the partial order. We compute a fresh topological order and
    return the matrix in the canonical convention so the result is
    a valid input to the rest of the toolchain.
    """

    n = len(matrix)
    # Permuted relation: new_rel[a][b] is True iff in the permuted
    # labeling, a prec b. Equivalently, original event permutation[a]
    # precedes permutation[b].
    permuted: List[List[bool]] = [[False] * n for _ in range(n)]
    for a in range(n):
        for b in range(n):
            if a == b:
                continue
            orig_a = permutation[a]
            orig_b = permutation[b]
            # The input matrix is upper-triangular under a topological
            # labeling. Event ``orig_a`` precedes ``orig_b`` iff
            # ``orig_a < orig_b`` and ``matrix[orig_a][orig_b]`` is True.
            if orig_a < orig_b and matrix[orig_a][orig_b]:
                permuted[a][b] = True

    # Now permuted[a][b] = True iff in the new labeling, the event
    # now labeled a precedes the event now labeled b. Some "a prec b"
    # relations may have a > b in the new labels. We need to relabel
    # so that a < b for all "a prec b" pairs, i.e., topologically
    # sort.
    in_degree = [0] * n
    for a in range(n):
        for b in range(n):
            if permuted[a][b]:
                in_degree[b] += 1
    order: List[int] = []
    available = sorted(i for i in range(n) if in_degree[i] == 0)
    while available:
        u = available.pop(0)
        order.append(u)
        for b in range(n):
            if permuted[u][b]:
                in_degree[b] -= 1
                if in_degree[b] == 0:
                    # Maintain sort order so the output is deterministic.
                    bisect_insort(available, b)
    if len(order) != n:
        raise RuntimeError("permuted relation is cyclic; this should not happen")

    # Build the canonical upper-triangular matrix in the new label order.
    new_index_of = {old: new for new, old in enumerate(order)}
    out: List[List[bool]] = [[False] * n for _ in range(n)]
    for a in range(n):
        for b in range(n):
            if permuted[a][b]:
                i_new = new_index_of[a]
                j_new = new_index_of[b]
                if i_new < j_new:
                    out[i_new][j_new] = True
    return out


def bisect_insort(seq: List[int], value: int) -> None:
    """In-place sorted insertion for the topological order helper.

    Avoids importing ``bisect`` just for one call so the helper is
    self-contained and easy to read inline.
    """

    lo, hi = 0, len(seq)
    while lo < hi:
        mid = (lo + hi) // 2
        if seq[mid] < value:
            lo = mid + 1
        else:
            hi = mid
    seq.insert(lo, value)


class LabelPermutationInvarianceTests(unittest.TestCase):
    def _check_invariance(self, matrix: List[List[bool]], n_perms: int) -> None:
        rng = random.Random(7)
        baseline = causet_invariants.invariants_fingerprint(matrix)
        n = len(matrix)
        for _ in range(n_perms):
            perm = list(range(n))
            rng.shuffle(perm)
            permuted = _topologically_sort_matrix(matrix, perm)
            other = causet_invariants.invariants_fingerprint(permuted)
            self.assertEqual(other, baseline)

    def test_chain_invariants_under_permutation(self) -> None:
        n = 6
        z = [[(i < j) for j in range(n)] for i in range(n)]
        self._check_invariance(z, n_perms=5)

    def test_sprinkled_invariants_under_permutation(self) -> None:
        matrix, _ = vs.sprinkle_minkowski_diamond(
            n=12, seed=1987, d_spacetime=2
        )
        self._check_invariance(matrix, n_perms=5)

    def test_higher_dim_sprinkled_invariants_under_permutation(self) -> None:
        matrix, _ = vs.sprinkle_minkowski_diamond(
            n=16, seed=2026, d_spacetime=3
        )
        self._check_invariance(matrix, n_perms=5)


class CompositePoincareInvarianceTests(unittest.TestCase):
    def test_interval_rmse_zero_under_composite_transform(self) -> None:
        # Translation + boost + spatial rotation + uniform rescale.
        # Each operation individually preserves the Lorentzian
        # interval (up to a global multiplicative factor for the
        # rescale, which is shared between both embeddings).
        _, points = vs.sprinkle_minkowski_diamond(
            n=12, seed=1959, d_spacetime=3
        )

        # Translation.
        translated = [(t + 1.0, x + 0.3, y - 0.2) for t, x, y in points]

        # 1+1 boost on (t, x).
        rapidity = 0.5
        ch = math.cosh(rapidity)
        sh = math.sinh(rapidity)
        boosted = []
        for t, x, y in translated:
            tp = ch * t - sh * x
            xp = -sh * t + ch * x
            boosted.append((tp, xp, y))

        # Spatial rotation on (x, y).
        theta = 0.7
        cs = math.cos(theta)
        sn = math.sin(theta)
        rotated = []
        for t, x, y in boosted:
            xp = x * cs - y * sn
            yp = x * sn + y * cs
            rotated.append((t, xp, yp))

        # Uniform rescale (factor of 2 on every coordinate).
        rescaled_a = [(2.0 * t, 2.0 * x, 2.0 * y) for t, x, y in points]
        rescaled_b = [(2.0 * t, 2.0 * x, 2.0 * y) for t, x, y in rotated]

        # Both embeddings now differ by translation+boost+rotation,
        # at the same global scale. Interval RMSE must be zero.
        self.assertAlmostEqual(
            vs.interval_rmse(rescaled_a, rescaled_b), 0.0, places=8
        )


if __name__ == "__main__":
    unittest.main()
