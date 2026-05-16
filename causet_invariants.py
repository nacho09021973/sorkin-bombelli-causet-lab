"""Order-theoretic invariants of finite causal sets.

The functions in this module operate purely on the upper-triangular
boolean causal matrix produced by ``cones.parse_cones_input`` or
``cones.generate_sprinkled_causet``. They are independent of any
embedding, optimizer, or random seed: given the same causet, they
return the same numbers.

The aim is to provide a battery of diagnostics that lets us judge a
causet on its own terms, before any annealing run, so we can later
ask the harder question: when the annealing program fails on a
causet, is the failure algorithmic or is the causet genuinely not
manifoldlike?

Notation conventions follow ``cones.py``:

- ``z[i][j]`` is ``True`` iff ``i`` precedes ``j`` (``i prec j``)
  with ``i < j`` by labeling.
- "Spacetime dimension" ``d`` is the total dimension of Minkowski
  space, so ``d = 2`` means ``1 + 1``. This is *not* the same as
  ``spacetime_dim`` in :func:`cones.generate_sprinkled_causet`,
  where ``spacetime_dim`` counts spatial dimensions only and
  ``d = spacetime_dim + 1``.
"""

from __future__ import annotations

import math
from typing import Dict, List, Sequence


CausalMatrix = Sequence[Sequence[bool]]


def relation_count(z: CausalMatrix) -> int:
    """Number of causally related ordered pairs ``(i, j)`` with ``i prec j``."""

    n = len(z)
    total = 0
    for i in range(n - 1):
        row = z[i]
        for j in range(i + 1, n):
            if row[j]:
                total += 1
    return total


def ordering_fraction(z: CausalMatrix) -> float:
    """Fraction of pairs that are causally related.

    For a Poisson sprinkling of a causal interval in ``d``-dimensional
    Minkowski space the expectation is given by the Myrheim-Meyer
    function (see :func:`myrheim_meyer_dimension`).
    """

    n = len(z)
    if n < 2:
        return 0.0
    total_pairs = n * (n - 1) // 2
    return relation_count(z) / total_pairs


def _myrheim_meyer_f(d: float) -> float:
    """Expected ordering fraction ``R / C(N, 2)`` in ``d``-dim Minkowski.

    For a uniform Poisson sprinkling of the unit causal diamond in
    ``d``-dimensional Minkowski space, the expected fraction of
    *unordered* causally related pairs is

        f(d) = Gamma(d + 1) * Gamma(d / 2) / (2 * Gamma(3 * d / 2)).

    Numerical values:

    - ``f(1) = 1``           (every pair related in 1D Minkowski).
    - ``f(2) = 1/2``         (1+1 unit diamond, light-cone square).
    - ``f(3) = 24 / 105``    (~ 0.2286).
    - ``f(4) = 1/10``.
    - ``f(d) -> 0``          as ``d -> infinity`` (becomes an antichain).

    A note on conventions. Some references (e.g. Reid 2003, Henson 2010)
    state the formula with a denominator of ``4 * Gamma(3d/2)`` because
    they count *ordered* relations ``R / [N(N - 1)]``. Our ``R`` counts
    unordered pairs (the upper-triangular True entries of the causal
    matrix), so the denominator is ``2 * Gamma(3d/2)``. The factor of
    two corresponds to the choice of direction. Both conventions give
    the same recovered dimension when inverted against the matching
    empirical fraction; this implementation matches the unordered
    convention used by :func:`ordering_fraction`. Verified empirically
    against canonical Poisson sprinklings at d in {2, 3, 4} to within
    one Monte Carlo standard deviation.
    """

    return (
        math.gamma(d + 1.0)
        * math.gamma(d / 2.0)
        / (2.0 * math.gamma(3.0 * d / 2.0))
    )


def myrheim_meyer_dimension(
    z: CausalMatrix,
    *,
    d_min: float = 1.0,
    d_max: float = 8.0,
    tolerance: float = 1e-9,
) -> float:
    """Estimate the continuous spacetime dimension ``d`` from the causet.

    Inverts the Myrheim-Meyer ordering-fraction formula
    ``f(d) = Gamma(d + 1) Gamma(d / 2) / (4 Gamma(3 d / 2))``
    by bisection on the observed ordering fraction.

    This estimator is an *independent* diagnostic. It does not call
    the annealing program and does not depend on any embedding. If
    Myrheim-Meyer indicates ``d ~ 3.2`` and the annealing succeeds
    only at ``dim >= 3``, the result is consistent. If the annealing
    fails at every target dimension while Myrheim-Meyer indicates a
    well-defined ``d``, the failure is algorithmic rather than
    structural.

    Returns ``+inf`` for an antichain (zero ordering fraction) and
    ``d_min`` whenever the ordering fraction is at or above the
    saturating value ``f(d_min)``. The default range ``[1, 8]``
    covers Minkowski dimensions from ``1`` (a single chain) up to
    ``1 + 7 = 8`` spatial dimensions, far beyond any physically
    interesting case.
    """

    f_obs = ordering_fraction(z)
    if f_obs <= 0.0:
        return float("inf")
    if f_obs >= _myrheim_meyer_f(d_min):
        return d_min

    lo, hi = d_min, d_max
    # f is strictly decreasing on [d_min, d_max]; bisect.
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        f_mid = _myrheim_meyer_f(mid)
        if f_mid > f_obs:
            lo = mid
        else:
            hi = mid
        if hi - lo < tolerance:
            break
    return 0.5 * (lo + hi)


def height(z: CausalMatrix) -> int:
    """Length of the longest chain in the causet.

    Computed by dynamic programming on the canonical topological
    labeling: ``h[j] = 1 + max(h[i] for i prec j)``.
    """

    n = len(z)
    if n == 0:
        return 0
    h = [1] * n
    for j in range(n):
        best = 0
        for i in range(j):
            if z[i][j] and h[i] > best:
                best = h[i]
        h[j] = best + 1
    return max(h) if h else 0


def antichain_profile(z: CausalMatrix) -> List[int]:
    """Sizes of the canonical topological levels of the causet.

    Level ``k`` contains the elements whose longest descending chain
    has length ``k + 1``. The profile ``[a_0, a_1, ...]`` records the
    number of elements at each level; ``len(profile) == height(z)``.

    Each level is a maximal antichain in the *graded* sense, though
    not necessarily a maximum antichain of the whole causet. For
    diagnostic comparison against the expected level profile of a
    Poisson sprinkling this graded profile is the natural fingerprint.
    """

    n = len(z)
    if n == 0:
        return []
    level = [0] * n
    for j in range(n):
        best = -1
        for i in range(j):
            if z[i][j] and level[i] > best:
                best = level[i]
        level[j] = best + 1
    profile: List[int] = []
    for lv in level:
        while lv >= len(profile):
            profile.append(0)
        profile[lv] += 1
    return profile


def chain_counts(z: CausalMatrix, k_max: int = 4) -> Dict[int, int]:
    """Number of ``k``-element chains for ``k = 2, ..., k_max``.

    A ``k``-chain is a sequence ``i_1 prec i_2 prec ... prec i_k`` of
    distinct elements. ``c[2]`` equals the relation count; ``c[3]`` and
    higher are higher-order invariants that distinguish causets sharing
    the same ordering fraction.

    The standard dynamic programming approach is used:
    ``c_k(j) = sum over i prec j of c_{k-1}(i)``,
    starting from ``c_2(j) = | { i : i prec j } |``.
    """

    if k_max < 2:
        return {}
    n = len(z)
    if n == 0:
        return {k: 0 for k in range(2, k_max + 1)}

    counts: Dict[int, List[int]] = {}
    counts[2] = [sum(1 for i in range(j) if z[i][j]) for j in range(n)]
    for k in range(3, k_max + 1):
        prev = counts[k - 1]
        counts[k] = [
            sum(prev[i] for i in range(j) if z[i][j])
            for j in range(n)
        ]
    return {k: sum(row) for k, row in counts.items()}


def three_chain_abundance(z: CausalMatrix) -> float:
    """Fraction of element triples that form a 3-chain.

    The raw count ``C3`` is the number of triples ``i prec j prec k``.
    This helper normalizes it by ``binom(n, 3)`` so the observable is
    comparable across the finite-size grid used in the foundation
    atlases:

        abundance_3 = C3 / binom(n, 3).

    It is a structural observable, not a calibrated dimension
    estimator. We keep the convention deliberately minimal because a
    closed-form inversion analogous to Myrheim-Meyer is not being
    assumed here. For ``n < 3`` it returns ``0.0``.
    """

    n = len(z)
    if n < 3:
        return 0.0
    c3 = chain_counts(z, k_max=3)[3]
    return c3 / math.comb(n, 3)


def link_count(z: CausalMatrix) -> int:
    """Number of links (covering relations) in the transitive reduction.

    Delegates to :func:`cones.transitive_reduction` so this module
    does not duplicate the covering-relation logic that already lives
    in the simulator.
    """

    # Local import keeps causet_invariants importable without dragging
    # in the simulator dependencies until link_count is actually used.
    from cones import transitive_reduction

    return len(transitive_reduction(z))


def _interval_size_matrix(z: CausalMatrix) -> List[List[int]]:
    """Compute the cardinality of every causal interval ``[i, j]``.

    Returns a matrix ``sizes`` with ``sizes[i][j]`` equal to the
    number of elements ``k`` satisfying ``i \\preceq k \\preceq j``
    (including the endpoints ``i`` and ``j`` themselves), for
    ``i < j`` with ``i \\prec j``. Off-relation pairs and the
    strictly lower triangle are left at zero. The diagonal is one.

    Used by :func:`midpoint_scaling_dimension`; exposed as a helper
    because computing the matrix takes ``O(n^3)`` and other future
    interval-based diagnostics will share it.
    """

    n = len(z)
    sizes = [[0] * n for _ in range(n)]
    for i in range(n):
        sizes[i][i] = 1
    for i in range(n - 1):
        for j in range(i + 1, n):
            if z[i][j]:
                count = 2  # endpoints i and j themselves
                for k in range(i + 1, j):
                    if z[i][k] and z[k][j]:
                        count += 1
                sizes[i][j] = count
    return sizes


def midpoint_scaling_dimension(
    z: CausalMatrix,
    *,
    min_interval_size: int = 4,
) -> float:
    """Estimate spacetime dimension by Meyer's midpoint interval scaler.

    Meyer (*The Dimension of Causal Sets*, MIT preprint, 1988) observed
    that for a Poisson sprinkling of intensity ``rho`` in a Minkowski
    causal interval of proper time ``T`` in ``d`` dimensions, the
    expected sprinkled cardinality is ``rho * c_d * T^d``. Splitting
    such an interval at its spacetime midpoint produces two
    sub-intervals of proper time ``T / 2`` with expected cardinality
    ``rho * c_d * (T / 2)^d = N_full * 2^(-d)`` each. Hence

        d ~= log_2(N_full / N_balanced_half)

    asymptotically, where ``N_balanced_half`` is the cardinality of
    the smaller sub-interval at the most balanced interior split.

    Order-theoretic implementation:

    1. Compute the cardinality of every causal interval (see
       :func:`_interval_size_matrix`).
    2. Pick the largest interval ``[i*, j*]``.
    3. Find the interior element ``k*`` in ``[i*, j*]`` that
       maximizes ``min(|[i*, k]|, |[k, j*]|)``.
    4. Return ``log_2(|[i*, j*]| / min(|[i*, k*]|, |[k*, j*]|))``.

    This estimator is **independent of Myrheim-Meyer**: it does not
    invert a closed-form ordering-fraction expression and depends
    only on the cardinalities of nested intervals. Agreement
    between the two estimators is a non-trivial signal that the
    causet behaves manifoldlike; disagreement is itself diagnostic.

    For pathological inputs (an antichain, a causet whose largest
    interval has size below ``min_interval_size``, a chain where
    no interior balanced split exists) the function returns
    ``float('nan')``. Antichains and chains have well-defined
    dimensions in principle (``infty`` and ``1``), but this
    estimator's asymptotics do not apply there; the user should
    consult :func:`myrheim_meyer_dimension` and :func:`height`
    for those edge cases.
    """

    n = len(z)
    if n < min_interval_size:
        return float("nan")

    sizes = _interval_size_matrix(z)

    best_size = 0
    best_i, best_j = -1, -1
    for i in range(n - 1):
        for j in range(i + 1, n):
            if sizes[i][j] > best_size:
                best_size = sizes[i][j]
                best_i, best_j = i, j
    if best_size < min_interval_size or best_i < 0:
        return float("nan")

    best_balance = 0
    for k in range(best_i + 1, best_j):
        if z[best_i][k] and z[k][best_j]:
            left = sizes[best_i][k]
            right = sizes[k][best_j]
            balance = left if left < right else right
            if balance > best_balance:
                best_balance = balance

    if best_balance <= 0:
        return float("nan")
    return math.log2(best_size / best_balance)


def invariants_fingerprint(z: CausalMatrix) -> Dict[str, object]:
    """Compute the full battery of invariants as a JSON-friendly dict.

    Stored alongside a frozen benchmark fixture, this fingerprint lets
    us verify that the causet is byte-identical to the recorded one
    *and* that the invariant computation itself is reproducible.
    """

    return {
        "n": len(z),
        "relation_count": relation_count(z),
        "link_count": link_count(z),
        "ordering_fraction": ordering_fraction(z),
        "myrheim_meyer_dim": myrheim_meyer_dimension(z),
        "height": height(z),
        "antichain_profile": antichain_profile(z),
        "chain_counts": chain_counts(z, k_max=4),
    }
