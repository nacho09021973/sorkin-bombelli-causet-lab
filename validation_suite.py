"""Sprinkle-then-recover validation for causet embedding methods.

This module supplies the ground-truth against which any embedding
algorithm in the project (the historical Bombelli annealing or any
v2 variant) is judged. The protocol is:

1. Generate a causet by uniform Poisson sprinkling of a known
   Minkowski causal diamond. The true coordinates are kept.
2. Run an embedding method on the causet (only the causal matrix
   is exposed to the method; the coordinates are hidden).
3. Compare the recovered embedding to the truth using a
   Lorentz-invariant residual (the squared-interval matrix RMSE).
4. Independently, compute the Myrheim-Meyer dimension on the
   causal matrix itself and verify it is consistent with the
   sprinkling dimension.

Two embeddings of the same causet that differ by a proper Poincare
transformation produce identical interval matrices. Comparing the
matrices avoids the need to solve the (open-form) Lorentzian
Procrustes problem explicitly. This is also the more Sorkinian
metric: it scores embeddings on Lorentz-invariant data.

A canonical sprinkler is provided here because the legacy
:func:`cones.generate_sprinkled_causet` does not sample the
diamond uniformly for ``d_spacetime >= 3``. Its marginal on the
time coordinate is uniform on ``[0, 1]``, but the canonical
uniform measure has marginal ``p(t) ~ min(t, 1 - t)^(d_spatial)``
(proportional to the spatial cross-section volume).
The two distributions agree only for ``d_spacetime = 2``, where
the (u, v) light-cone square maps to the diamond with constant
Jacobian.
"""

from __future__ import annotations

import contextlib
import io
import math
import random
import struct
import tempfile
import zlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Sequence, Tuple

import causet_invariants
import cones


Coord = Tuple[float, ...]
CausalMatrix = List[List[bool]]


# ---------------------------------------------------------------------
# Canonical Minkowski-diamond sprinkler
# ---------------------------------------------------------------------


def _sample_spatial_in_ball(rng: random.Random, d_spatial: int) -> Tuple[Tuple[float, ...], float]:
    """Sample a point uniformly in the open unit ball in R^d_spatial.

    Returns ``(x, r)`` where ``x`` is the point and ``r = |x|``.
    Uses rejection sampling against the unit cube.
    """

    while True:
        x = tuple(2.0 * rng.random() - 1.0 for _ in range(d_spatial))
        r2 = sum(c * c for c in x)
        if r2 <= 1.0:
            return x, math.sqrt(r2)


def sprinkle_minkowski_diamond(
    n: int,
    seed: int,
    d_spacetime: int,
) -> Tuple[CausalMatrix, List[Coord]]:
    """Uniform Poisson sprinkling of the unit Minkowski causal diamond.

    The diamond is the Alexandrov interval from ``(0, 0, ..., 0)`` to
    ``(1, 0, ..., 0)`` in ``d_spacetime``-dimensional Minkowski space.
    Sampling is by rejection against the bounding box
    ``[0, 1] x B^{d_spacetime - 1}`` (time interval and unit spatial
    ball), which guarantees the resulting distribution is exactly
    uniform over the diamond's Lebesgue measure.

    Parameters
    ----------
    n: number of sprinkled events.
    seed: integer seed for the random number generator.
    d_spacetime: total spacetime dimension. ``d_spacetime = 2`` means
        1+1 Minkowski (one time, one space).

    Returns
    -------
    A pair ``(matrix, points)``. ``matrix`` is the upper-triangular
    causal matrix in the convention used everywhere else in the
    project (``matrix[i][j] = True`` iff ``i prec j`` and ``i < j``),
    after sorting events by ``t``. ``points[i]`` is ``(t, x_1, ...,
    x_{d_spacetime - 1})``.

    Notes
    -----
    This function deliberately uses Python's :class:`random.Random`,
    *not* the simulator's PascalRNG, because the sprinkling is part
    of the experimental design rather than the embedding algorithm.
    Reproducibility within this sprinkler is per-seed, independent of
    any other random stream.
    """

    if n <= 0:
        raise ValueError("n must be positive")
    if d_spacetime < 2:
        raise ValueError("d_spacetime must be at least 2 (use a chain for d=1)")

    rng = random.Random(seed)
    d_spatial = d_spacetime - 1

    points: List[Coord] = []
    while len(points) < n:
        t = rng.random()
        x, r = _sample_spatial_in_ball(rng, d_spatial)
        if r <= min(t, 1.0 - t):
            points.append((t, *x))

    points.sort(key=lambda p: (p[0],) + p[1:])

    matrix: CausalMatrix = [[False] * n for _ in range(n)]
    for i in range(n - 1):
        ti = points[i][0]
        xi = points[i][1:]
        for j in range(i + 1, n):
            tj = points[j][0]
            xj = points[j][1:]
            dt = tj - ti
            if dt < 0.0:
                continue  # impossible after the sort, but guard anyway
            dx_sq = 0.0
            for a, b in zip(xi, xj):
                dx_sq += (b - a) * (b - a)
            if dt * dt >= dx_sq:
                matrix[i][j] = True
    return matrix, points


# ---------------------------------------------------------------------
# Non-manifoldlike control: Kleitman-Rothschild posets
# ---------------------------------------------------------------------


def generate_kleitman_rothschild(
    n: int,
    seed: int,
    *,
    edge_prob: float = 0.5,
    level_sizes: Tuple[int, int, int] | None = None,
) -> CausalMatrix:
    """Generate a three-level Kleitman-Rothschild-style poset.

    Kleitman and Rothschild (1975, *Asymptotic enumeration of partial
    orders on a finite set*, Trans. AMS 205, 205-220) proved that the
    fraction of labeled posets on ``n`` elements with exactly three
    antichain levels approaches one as ``n -> infty``, and that the
    typical level sizes are concentrated near ``(n/4, n/2, n/4)``.
    Such posets are emphatically *non-manifoldlike*: they cannot be
    embedded faithfully in a low-dimensional Minkowski causal
    interval. They are included here as a **diagnostic control**,
    not as a model of spacetime.

    The sampling procedure is the standard one used in numerical
    studies of non-manifoldlike causets:

    1. Partition ``n`` elements into three levels of sizes
       approximately ``(n/4, n/2, n/4)``, with the canonical
       topological labeling (level 0 first, then level 1, then
       level 2). Override via ``level_sizes`` if needed.
    2. For each pair ``(a, b)`` with ``a`` in level 0 and ``b`` in
       level 1, include the relation ``a prec b`` with probability
       ``edge_prob``.
    3. Likewise between level 1 and level 2.
    4. Take transitive closure.

    Parameters
    ----------
    n: number of elements, ``n >= 3``.
    seed: integer seed for the random number generator.
    edge_prob: probability of an inter-level edge, default 0.5.
    level_sizes: optional explicit ``(bottom, middle, top)``
        triple; must sum to ``n``. Default rounds to
        ``(ceil(n/4), n - 2 ceil(n/4), ceil(n/4))``.

    Returns
    -------
    The upper-triangular causal matrix in the project's convention.
    The labeling is already topologically sorted by construction:
    level-0 indices precede level-1, which precede level-2, and
    there are no relations within a level.

    Notes
    -----
    This generator deliberately fixes the level structure rather
    than sampling uniformly over all labeled posets. The uniform
    measure has no closed-form sampler; the fixed three-level
    construction with random inter-level edges captures the same
    asymptotic behavior and is the convention used throughout the
    poset-enumeration literature (see also Brightwell and
    Winkler 1991).
    """

    if n < 3:
        raise ValueError("n must be at least 3 for a three-level structure")
    if not (0.0 <= edge_prob <= 1.0):
        raise ValueError("edge_prob must be in [0, 1]")

    if level_sizes is None:
        bottom = (n + 3) // 4  # round-half-up of n / 4
        top = bottom
        middle = n - bottom - top
        if middle < 1:
            third = n // 3
            bottom = third
            top = third
            middle = n - bottom - top
        level_sizes = (bottom, middle, top)
    bottom, middle, top = level_sizes
    if bottom < 1 or middle < 1 or top < 1:
        raise ValueError(
            f"level_sizes {level_sizes} must each be at least 1"
        )
    if bottom + middle + top != n:
        raise ValueError(
            f"level_sizes {level_sizes} must sum to n={n}"
        )

    rng = random.Random(seed)

    # Indices 0..bottom-1 are level 0; bottom..bottom+middle-1 are
    # level 1; the rest are level 2. This labeling is topologically
    # sorted by construction.
    level0_end = bottom
    level1_end = bottom + middle

    matrix: CausalMatrix = [[False] * n for _ in range(n)]

    # Direct edges between consecutive levels.
    for a in range(level0_end):
        for b in range(level0_end, level1_end):
            if rng.random() < edge_prob:
                matrix[a][b] = True
    for b in range(level0_end, level1_end):
        for c in range(level1_end, n):
            if rng.random() < edge_prob:
                matrix[b][c] = True

    # Transitive closure: a prec c if there exists b in level 1 with
    # a prec b and b prec c.
    for a in range(level0_end):
        for c in range(level1_end, n):
            for b in range(level0_end, level1_end):
                if matrix[a][b] and matrix[b][c]:
                    matrix[a][c] = True
                    break
    return matrix


def generate_corona_poset(n: int, seed: int) -> CausalMatrix:
    """Generate a suspended corona/crown non-manifoldlike control.

    The construction is a two-level crown poset with a global
    minimum and maximum adjoined:

    - element ``0`` is below every other element;
    - element ``n - 1`` is above every other element;
    - the remaining elements are split into lower and upper corona
      layers;
    - each lower-layer element precedes every upper-layer element
      except one seed-dependent forbidden partner.

    This gives a cheap negative control with large causal intervals
    (so the midpoint scaler is finite) but no manifoldlike nested
    interval structure. The relation density is asymptotically
    stable, while the largest interval's best balanced split remains
    shallow, making it a useful check that the KR finite-size
    signature is not an artifact of one three-level generator.
    """

    if n < 6:
        raise ValueError("n must be at least 6 for a suspended corona")

    interior = n - 2
    lower = interior // 2
    upper = interior - lower
    if lower < 2 or upper < 2:
        raise ValueError("corona layers must each contain at least 2 elements")

    rng = random.Random(seed)
    forbidden = list(range(upper))
    rng.shuffle(forbidden)

    matrix: CausalMatrix = [[False] * n for _ in range(n)]
    top = n - 1
    lower_start = 1
    upper_start = lower_start + lower

    for j in range(1, n):
        matrix[0][j] = True
    for i in range(0, top):
        matrix[i][top] = True

    for a_offset in range(lower):
        a = lower_start + a_offset
        skipped = forbidden[a_offset % upper]
        for b_offset in range(upper):
            if b_offset == skipped:
                continue
            b = upper_start + b_offset
            matrix[a][b] = True

    return matrix


# ---------------------------------------------------------------------
# Lorentz-invariant embedding residual
# ---------------------------------------------------------------------


def minkowski_interval_matrix(points: Sequence[Coord]) -> List[List[float]]:
    """Matrix of squared Minkowski intervals between all event pairs.

    With signature ``(-, +, +, ..., +)`` the squared interval between
    events ``i`` and ``j`` is
    ``s^2 = -(t_j - t_i)^2 + sum_k (x_j[k] - x_i[k])^2``.

    The full symmetric matrix is returned. Diagonal entries are zero.

    Two embeddings of the same causet related by a proper Poincare
    transformation share identical squared-interval matrices. The
    pairwise residual of this matrix is therefore the canonical
    Lorentz-invariant scalar with which to compare an embedding to
    a ground truth without solving for the transformation.
    """

    n = len(points)
    s2: List[List[float]] = [[0.0] * n for _ in range(n)]
    if n < 2:
        return s2
    for i in range(n):
        ti = points[i][0]
        xi = points[i][1:]
        for j in range(i + 1, n):
            tj = points[j][0]
            xj = points[j][1:]
            dt = tj - ti
            dx_sq = 0.0
            for a, b in zip(xi, xj):
                dx_sq += (b - a) * (b - a)
            value = -dt * dt + dx_sq
            s2[i][j] = value
            s2[j][i] = value
    return s2


def interval_rmse(
    points_a: Sequence[Coord],
    points_b: Sequence[Coord],
) -> float:
    """Lorentz-invariant residual between two embeddings of one causet.

    Defined as the root-mean-square of the elementwise difference of
    the squared-interval matrices, averaged over unordered pairs.

    Returns 0.0 iff the two embeddings are related by an isometry of
    Minkowski space (proper or improper Poincare). It is *not*
    sensitive to which embedding "looks more like a sprinkling"; for
    that, use :func:`causet_invariants.myrheim_meyer_dimension` on
    the causal matrix and compare against the truth dimension.
    """

    if len(points_a) != len(points_b):
        raise ValueError(
            f"point lists differ in length: {len(points_a)} vs {len(points_b)}"
        )
    n = len(points_a)
    if n < 2:
        return 0.0
    s2_a = minkowski_interval_matrix(points_a)
    s2_b = minkowski_interval_matrix(points_b)
    sq_sum = 0.0
    pair_count = 0
    for i in range(n):
        for j in range(i + 1, n):
            diff = s2_a[i][j] - s2_b[i][j]
            sq_sum += diff * diff
            pair_count += 1
    return math.sqrt(sq_sum / pair_count)


# ---------------------------------------------------------------------
# Bombelli energy at a fixed embedding (no annealing)
# ---------------------------------------------------------------------


def bombelli_energy_at(
    z: CausalMatrix,
    points: Sequence[Coord],
    *,
    d_spatial: int | None = None,
) -> float:
    """Bombelli energy of a fixed configuration of coordinates.

    Replicates :meth:`cones.ConesSimulator.energy` exactly without
    starting an annealing run. Useful for asking the diagnostic
    question "is the ground-truth embedding actually a low-energy
    configuration?". If the ground-truth energy is comparable to a
    failed annealing's final energy, the optimizer reached the
    bottom and the energy function itself is the issue. If the
    ground-truth energy is much lower, the optimizer is stuck.

    The function is scale-invariant: rescaling all spatial and
    time coordinates by the same factor leaves the energy unchanged
    (the formula's terms each scale linearly and the ``rave``
    denominator absorbs the factor). The configuration must
    satisfy ``points[i][0] > 0`` for every ``i``, matching the
    simulator's positivity convention on the time-like radius.
    """

    n = len(z)
    if n == 0:
        return 0.0
    if d_spatial is None:
        d_spatial = max(len(p) - 1 for p in points)

    rave = sum(p[0] for p in points) / n
    if rave <= 0.0:
        raise ValueError("average time coordinate must be positive")

    roottwo = math.sqrt(2.0)
    total = 0.0
    for i in range(n - 1):
        ri = points[i][0]
        xi = points[i][1:]
        for j in range(i + 1, n):
            rj = points[j][0]
            xj = points[j][1:]
            # Match cones.energy exactly: rij = rnew[i] - rnew[j].
            rij = ri - rj
            xij_sq = 0.0
            for k in range(d_spatial):
                a = xi[k] if k < len(xi) else 0.0
                b = xj[k] if k < len(xj) else 0.0
                xij_sq += (a - b) * (a - b)
            s2 = -(rij ** 2) + xij_sq
            xij = math.sqrt(max(xij_sq, 0.0))
            if z[i][j]:
                if s2 > 0.0:
                    total += (xij + rij) / (roottwo * rave)
                elif rij > 0.0:
                    total += math.sqrt(s2 + 2.0 * (rij ** 2)) / rave
                else:
                    pass  # correctly timelike, no penalty
            else:
                if s2 > 0.0:
                    pass  # correctly spacelike, no penalty
                else:
                    total += (abs(rij) - xij) / (roottwo * rave)
    return total


# ---------------------------------------------------------------------
# Sprinkle-then-recover protocol
# ---------------------------------------------------------------------


@dataclass(frozen=True)
class SprinkleCase:
    """A single frozen sprinkling: causal structure and ground truth.

    The causal matrix is what the optimizer sees; the points are the
    ground truth against which the recovered embedding is compared.
    Both are reproducible from ``(d_spacetime, n, seed)``.
    """

    d_spacetime: int
    n: int
    seed: int
    matrix: CausalMatrix
    points: List[Coord]


@dataclass
class RecoveryResult:
    """One attempt to recover the embedding of one SprinkleCase."""

    d_spacetime: int
    n: int
    case_seed: int
    optimizer_seed: int
    target_dim: int
    initial_energy: float
    warmup_energy: float
    final_energy: float
    truth_energy: float
    interval_rmse: float
    mm_dim_truth: float
    mm_dim_recovered: float


def generate_ensemble(
    d_spacetimes: Sequence[int],
    sizes: Sequence[int],
    seeds: Sequence[int],
) -> List[SprinkleCase]:
    """Build a grid of canonical sprinklings.

    The Cartesian product of ``d_spacetimes x sizes x seeds`` is
    produced. Each cell is one :class:`SprinkleCase`.
    """

    cases: List[SprinkleCase] = []
    for d in d_spacetimes:
        for n in sizes:
            for s in seeds:
                matrix, points = sprinkle_minkowski_diamond(
                    n=n, seed=s, d_spacetime=d
                )
                cases.append(
                    SprinkleCase(
                        d_spacetime=d,
                        n=n,
                        seed=s,
                        matrix=matrix,
                        points=points,
                    )
                )
    return cases


def _recovered_embedding(sim: cones.ConesSimulator) -> List[Coord]:
    """Read the recovered Minkowski coordinates out of a simulator.

    The simulator stores the time-like radius in ``rnew`` and the
    spatial vector in ``xnew``. We return them as ``(t, x_1, ..., x_d)``
    so the result can be compared directly to a :class:`SprinkleCase`
    points list.
    """

    return [(sim.rnew[i], *sim.xnew[i]) for i in range(sim.n)]


def run_recovery(
    case: SprinkleCase,
    *,
    optimizer_seed: int,
    target_dim: int | None = None,
    warmup_limit: int = 100,
    anneal_limit: int = 100,
    initial_temp: float = 100.0,
    cooling_factor: float = 0.9,
    max_data: int = 35,
    backend: str = "cpu",
) -> RecoveryResult:
    """Run one annealing recovery on one sprinkled case.

    Records the Bombelli energy at the ground-truth coordinates as
    well as the optimizer's final energy, the Lorentz-invariant
    interval RMSE between recovered and truth, and the Myrheim-Meyer
    dimension of the causal matrix (an embedding-free diagnostic).
    """

    d_spatial = case.d_spacetime - 1
    if target_dim is None:
        target_dim = d_spatial

    with tempfile.TemporaryDirectory() as tmpdir, contextlib.redirect_stdout(io.StringIO()):
        sim = cones.ConesSimulator(
            z=case.matrix,
            dim=target_dim,
            seed=optimizer_seed,
            interactive=False,
            max_data=max_data,
            plot_path=None,
            warmup_limit=warmup_limit,
            anneal_limit=anneal_limit,
            initial_temp=initial_temp,
            cooling_factor=cooling_factor,
            backend=backend,
        )
        sim.run(Path(tmpdir) / "out.txt")

    recovered = _recovered_embedding(sim)

    # The simulator's recovered embedding is in (r, x) where r is the
    # rescaled time-like radius. The ground-truth points are in (t, x)
    # with t in [0, 1]. The two are at different scales but the
    # interval matrix RMSE is invariant under uniform rescaling, so we
    # rescale the ground truth to the same scale as the simulator
    # before comparison.
    rave_truth = sum(p[0] for p in case.points) / case.n
    rave_recovered = sum(p[0] for p in recovered) / sim.n
    if rave_truth > 0.0 and rave_recovered > 0.0:
        scale = rave_recovered / rave_truth
        truth_scaled = [
            tuple(c * scale for c in p) for p in case.points
        ]
    else:
        truth_scaled = list(case.points)

    interval_residual = interval_rmse(recovered, truth_scaled)
    truth_energy = bombelli_energy_at(
        case.matrix, truth_scaled, d_spatial=target_dim
    )

    return RecoveryResult(
        d_spacetime=case.d_spacetime,
        n=case.n,
        case_seed=case.seed,
        optimizer_seed=optimizer_seed,
        target_dim=target_dim,
        initial_energy=sim.initial_energy,
        warmup_energy=sim.warmup_energy,
        final_energy=sim.data[-1][1] if sim.data else float("nan"),
        truth_energy=truth_energy,
        interval_rmse=interval_residual,
        mm_dim_truth=float(case.d_spacetime),
        mm_dim_recovered=causet_invariants.myrheim_meyer_dimension(case.matrix),
    )


# ---------------------------------------------------------------------
# SORKIN-2 order verification
# ---------------------------------------------------------------------


@dataclass(frozen=True)
class OrderComparison:
    """Binary comparison of a target causal order against an induced order.

    Produced by :func:`compare_causal_orders`. All fields are exact
    relation counts or lists of offending pairs; no energy, no
    continuous observable, and no dimension estimator is involved.

    Fields
    ------
    n:
        Number of events.
    total_relations_target:
        Number of True entries in the upper-triangular target matrix.
    total_relations_induced:
        Number of True entries in the upper-triangular induced matrix.
    missing_relations:
        Pairs ``(i, j)`` with ``i < j`` that are present in the target
        order but absent in the induced order.
    extra_relations:
        Pairs ``(i, j)`` with ``i < j`` that are absent in the target
        order but present in the induced order.
    exact_match:
        ``True`` iff ``missing_relations`` and ``extra_relations`` are
        both empty, i.e. the induced order reproduces the target exactly.
    """

    n: int
    total_relations_target: int
    total_relations_induced: int
    missing_relations: Tuple[Tuple[int, int], ...]
    extra_relations: Tuple[Tuple[int, int], ...]
    exact_match: bool


@dataclass(frozen=True)
class OrderMatrixPlotPaths:
    """PNG outputs for an order-matrix comparison."""

    target_order_matrix: Path
    induced_order_matrix: Path
    order_difference_matrix: Path


def induced_order_from_coords(coords: Sequence[Coord]) -> CausalMatrix:
    """Causal order induced by a set of annealer coordinates.

    Accepts coordinates in the format used by :class:`~cones.ConesSimulator`:
    each entry is ``(r_i, x_i_1, ..., x_i_d)`` where ``r_i`` is the
    time-like radius (``sim.rold[i]``) and ``x_i_k`` are the spatial
    components (``sim.xold[i][k]``).

    In the annealer convention, event ``i`` causally precedes event ``j``
    iff both conditions hold:

    .. code-block:: text

        r_i <= r_j                                   (earlier label has no larger r)
        (r_i - r_j)^2 >= sum_k (x_ik - x_jk)^2     (timelike or null interval)

    This matches the zero-energy condition in
    :meth:`~cones.ConesSimulator.energy` exactly: for a pair with
    ``z[i][j] = True``, the Bombelli energy contribution is zero iff
    ``r_i <= r_j`` and the Minkowski interval is non-positive.

    The function is independent of the Bombelli energy: it uses only the
    geometric Minkowski causal relation applied to the coordinate values.
    It does not call the simulator, does not use the energy matrix, and
    does not invoke any invariant estimator.

    Returns the upper-triangular causal matrix: ``z[i][j]`` is ``True``
    for ``i < j`` (by label) iff event ``i`` causally precedes event ``j``
    in the given geometry.
    """
    n = len(coords)
    matrix: CausalMatrix = [[False] * n for _ in range(n)]
    for i in range(n - 1):
        ri = coords[i][0]
        xi = coords[i][1:]
        for j in range(i + 1, n):
            rj = coords[j][0]
            xj = coords[j][1:]
            rij = ri - rj
            if rij > 0.0:
                continue
            dx_sq = sum((a - b) * (a - b) for a, b in zip(xi, xj))
            if rij * rij >= dx_sq:
                matrix[i][j] = True
    return matrix


def compare_causal_orders(
    z_target: CausalMatrix,
    z_induced: CausalMatrix,
) -> OrderComparison:
    """Compare two upper-triangular causal order matrices element by element.

    Returns an :class:`OrderComparison` with the sets of missing and extra
    relations. This is a purely order-theoretic comparison: no energy,
    no coordinates, no dimension estimator is involved.

    Parameters
    ----------
    z_target:
        The target causal order matrix (e.g., from a ``.in`` file or a
        known sprinkling). ``z_target[i][j]`` is ``True`` iff event ``i``
        should causally precede event ``j``.
    z_induced:
        The causal order matrix computed from a recovered embedding
        (e.g., from :func:`induced_order_from_coords`).

    Returns
    -------
    An :class:`OrderComparison` with ``exact_match=True`` iff the two
    matrices agree on every pair ``(i, j)`` with ``i < j``.

    Raises
    ------
    ValueError:
        If the two matrices have different sizes.
    """
    n = len(z_target)
    if len(z_induced) != n:
        raise ValueError(
            f"matrix sizes differ: target n={n}, induced n={len(z_induced)}"
        )
    total_target = 0
    total_induced = 0
    missing: list[Tuple[int, int]] = []
    extra: list[Tuple[int, int]] = []
    for i in range(n - 1):
        for j in range(i + 1, n):
            t = bool(z_target[i][j])
            ind = bool(z_induced[i][j])
            if t:
                total_target += 1
            if ind:
                total_induced += 1
            if t and not ind:
                missing.append((i, j))
            elif ind and not t:
                extra.append((i, j))
    return OrderComparison(
        n=n,
        total_relations_target=total_target,
        total_relations_induced=total_induced,
        missing_relations=tuple(missing),
        extra_relations=tuple(extra),
        exact_match=len(missing) == 0 and len(extra) == 0,
    )


def verify_recovery(sim: cones.ConesSimulator) -> OrderComparison:
    """Check whether a simulator's final accepted state reproduces the target order.

    Extracts the last-accepted coordinates (``rold``, ``xold``) from the
    simulator, computes the induced causal order via
    :func:`induced_order_from_coords`, and compares it against the target
    order ``sim.z`` via :func:`compare_causal_orders`.

    Uses ``rold``/``xold`` (last accepted state) rather than
    ``rnew``/``xnew`` (last proposed state), consistent with what
    :meth:`~cones.ConesSimulator.writeout` reports as the final
    configuration.

    This is the canonical SORKIN-2 recoverability check. A result with
    ``exact_match=True`` means the recovered configuration induces
    exactly the target causal order, independently of the final Bombelli
    energy value.

    Parameters
    ----------
    sim:
        A :class:`~cones.ConesSimulator` after
        :meth:`~cones.ConesSimulator.run` has been called.

    Returns
    -------
    :class:`OrderComparison` with the exact-match verdict and the full
    lists of missing and extra relations.
    """
    coords: List[Coord] = [(sim.rold[i], *sim.xold[i]) for i in range(sim.n)]
    z_induced = induced_order_from_coords(coords)
    return compare_causal_orders(sim.z, z_induced)


# ---------------------------------------------------------------------
# SORKIN-2 order-matrix diagnostic figures
# ---------------------------------------------------------------------


Color = Tuple[int, int, int]


def _coerce_square_bool_matrix(matrix: Sequence[Sequence[object]], name: str) -> CausalMatrix:
    """Return a square matrix of bools from a list-like matrix."""

    n = len(matrix)
    rows: CausalMatrix = []
    for i, row in enumerate(matrix):
        if len(row) != n:
            raise ValueError(
                f"{name} must be square: row {i} has length {len(row)}, expected {n}"
            )
        rows.append([bool(cell) for cell in row])
    return rows


def _write_png_rgb(path: Path, pixels: Sequence[Sequence[Color]]) -> None:
    """Write an RGB PNG using only the standard library."""

    height = len(pixels)
    width = len(pixels[0]) if height else 0
    if width <= 0 or height <= 0:
        raise ValueError("PNG image must have positive width and height")
    for row in pixels:
        if len(row) != width:
            raise ValueError("PNG rows must all have the same width")

    def chunk(kind: bytes, data: bytes) -> bytes:
        return (
            struct.pack(">I", len(data))
            + kind
            + data
            + struct.pack(">I", zlib.crc32(kind + data) & 0xFFFFFFFF)
        )

    raw = bytearray()
    for row in pixels:
        raw.append(0)  # filter type 0
        for red, green, blue in row:
            raw.extend((red, green, blue))

    png = b"\x89PNG\r\n\x1a\n"
    png += chunk("IHDR".encode("ascii"), struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0))
    png += chunk("IDAT".encode("ascii"), zlib.compress(bytes(raw)))
    png += chunk("IEND".encode("ascii"), b"")
    path.write_bytes(png)


def _matrix_to_pixels(codes: Sequence[Sequence[int]], palette: dict[int, Color]) -> List[List[Color]]:
    """Expand a small matrix into a readable pixel grid."""

    n = len(codes)
    cell = max(18, min(48, 360 // max(n, 1)))
    grid = 2
    background = (245, 245, 242)
    line = (70, 70, 70)
    size = n * cell + (n + 1) * grid
    pixels: List[List[Color]] = [[background for _ in range(size)] for _ in range(size)]

    for r in range(size):
        for c in range(size):
            if r % (cell + grid) < grid or c % (cell + grid) < grid:
                pixels[r][c] = line

    for i, row in enumerate(codes):
        if len(row) != n:
            raise ValueError("plot matrix must be square")
        for j, code in enumerate(row):
            color = palette[code]
            y0 = grid + i * (cell + grid)
            x0 = grid + j * (cell + grid)
            for y in range(y0, y0 + cell):
                for x in range(x0, x0 + cell):
                    pixels[y][x] = color
    return pixels


def _write_matrix_png(
    matrix: Sequence[Sequence[int]],
    path: Path,
    *,
    title: str,
    palette: dict[int, Color],
    matplotlib_cmap: object | None = None,
    vmin: int = 0,
    vmax: int = 1,
) -> Path:
    """Write one matrix plot, preferring matplotlib when available."""

    try:
        import matplotlib

        matplotlib.use("Agg", force=True)
        from matplotlib.colors import ListedColormap
        import matplotlib.pyplot as plt
    except Exception:
        _write_png_rgb(path, _matrix_to_pixels(matrix, palette))
        return path

    if matplotlib_cmap is None:
        ordered = [palette[i] for i in range(vmin, vmax + 1)]
        matplotlib_cmap = ListedColormap(
            [(red / 255.0, green / 255.0, blue / 255.0) for red, green, blue in ordered]
        )

    fig, ax = plt.subplots(figsize=(4.0, 4.0), constrained_layout=True)
    ax.imshow(matrix, cmap=matplotlib_cmap, vmin=vmin, vmax=vmax, interpolation="nearest")
    ax.set_title(title)
    ax.set_xlabel("j")
    ax.set_ylabel("i")
    n = len(matrix)
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.grid(which="minor", color="black", linewidth=0.5)
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


def write_order_matrix_plots(
    z_target: Sequence[Sequence[object]],
    z_induced: Sequence[Sequence[object]],
    output_dir: Path,
) -> OrderMatrixPlotPaths:
    """Write the three permanent order-matrix diagnostic PNGs.

    The inputs are already-computed causal order matrices. This function does
    not run the annealer, does not evaluate Bombelli energy, and does not call
    any dimension estimator.

    The difference plot uses four discrete classes:

    - no relation in either matrix;
    - matching relation in both matrices;
    - missing relation, present in target but absent in induced;
    - extra relation, absent in target but present in induced.
    """

    target = _coerce_square_bool_matrix(z_target, "z_target")
    induced = _coerce_square_bool_matrix(z_induced, "z_induced")
    if len(target) != len(induced):
        raise ValueError(
            f"matrix sizes differ: target n={len(target)}, induced n={len(induced)}"
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    target_path = output_dir / "target_order_matrix.png"
    induced_path = output_dir / "induced_order_matrix.png"
    difference_path = output_dir / "order_difference_matrix.png"

    target_codes = [[1 if cell else 0 for cell in row] for row in target]
    induced_codes = [[1 if cell else 0 for cell in row] for row in induced]
    difference_codes: List[List[int]] = []
    for target_row, induced_row in zip(target, induced):
        diff_row: List[int] = []
        for target_cell, induced_cell in zip(target_row, induced_row):
            if target_cell and induced_cell:
                diff_row.append(1)
            elif target_cell and not induced_cell:
                diff_row.append(2)
            elif induced_cell and not target_cell:
                diff_row.append(3)
            else:
                diff_row.append(0)
        difference_codes.append(diff_row)

    binary_palette = {
        0: (245, 245, 242),
        1: (34, 94, 168),
    }
    difference_palette = {
        0: (245, 245, 242),  # no relation in either matrix
        1: (35, 118, 82),    # matching relation
        2: (217, 95, 2),     # missing relation
        3: (117, 112, 179),  # extra relation
    }

    _write_matrix_png(
        target_codes,
        target_path,
        title="Target causal order",
        palette=binary_palette,
        matplotlib_cmap="Blues",
        vmin=0,
        vmax=1,
    )
    _write_matrix_png(
        induced_codes,
        induced_path,
        title="Induced causal order",
        palette=binary_palette,
        matplotlib_cmap="Blues",
        vmin=0,
        vmax=1,
    )
    _write_matrix_png(
        difference_codes,
        difference_path,
        title="Order difference",
        palette=difference_palette,
        matplotlib_cmap=None,
        vmin=0,
        vmax=3,
    )

    return OrderMatrixPlotPaths(
        target_order_matrix=target_path,
        induced_order_matrix=induced_path,
        order_difference_matrix=difference_path,
    )


# ---------------------------------------------------------------------
# CSV / Markdown output
# ---------------------------------------------------------------------


_RESULT_CSV_HEADER = [
    "d_spacetime",
    "n",
    "case_seed",
    "optimizer_seed",
    "target_dim",
    "initial_energy",
    "warmup_energy",
    "final_energy",
    "truth_energy",
    "interval_rmse",
    "mm_dim_truth",
    "mm_dim_recovered",
]


def write_results_csv(results: Sequence[RecoveryResult], path: Path) -> None:
    """Write a CSV table of recovery results."""

    lines = [",".join(_RESULT_CSV_HEADER)]
    for r in results:
        lines.append(
            ",".join(
                [
                    str(r.d_spacetime),
                    str(r.n),
                    str(r.case_seed),
                    str(r.optimizer_seed),
                    str(r.target_dim),
                    f"{r.initial_energy:.6f}",
                    f"{r.warmup_energy:.6f}",
                    f"{r.final_energy:.6f}",
                    f"{r.truth_energy:.6f}",
                    f"{r.interval_rmse:.6f}",
                    f"{r.mm_dim_truth:.6f}",
                    f"{r.mm_dim_recovered:.6f}",
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def summarize(results: Sequence[RecoveryResult]) -> List[dict]:
    """Group results by ``(d_spacetime, n)`` and aggregate statistics."""

    buckets: dict = {}
    for r in results:
        key = (r.d_spacetime, r.n)
        buckets.setdefault(key, []).append(r)
    summary: List[dict] = []
    for (d, n), rows in sorted(buckets.items()):
        final_energies = [r.final_energy for r in rows]
        truth_energies = [r.truth_energy for r in rows]
        residuals = [r.interval_rmse for r in rows]
        mm_estimates = [r.mm_dim_recovered for r in rows]
        summary.append({
            "d_spacetime": d,
            "n": n,
            "runs": len(rows),
            "mean_final_energy": sum(final_energies) / len(final_energies),
            "mean_truth_energy": sum(truth_energies) / len(truth_energies),
            "mean_interval_rmse": sum(residuals) / len(residuals),
            "mean_mm_dim": sum(mm_estimates) / len(mm_estimates),
        })
    return summary


def write_summary_csv(summary: Sequence[dict], path: Path) -> None:
    headers = [
        "d_spacetime",
        "n",
        "runs",
        "mean_final_energy",
        "mean_truth_energy",
        "mean_interval_rmse",
        "mean_mm_dim",
    ]
    lines = [",".join(headers)]
    for row in summary:
        lines.append(
            ",".join(
                [
                    str(row["d_spacetime"]),
                    str(row["n"]),
                    str(row["runs"]),
                    f"{row['mean_final_energy']:.6f}",
                    f"{row['mean_truth_energy']:.6f}",
                    f"{row['mean_interval_rmse']:.6f}",
                    f"{row['mean_mm_dim']:.6f}",
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_summary_report(summary: Sequence[dict], path: Path) -> None:
    """Write a markdown report of recovery summary statistics."""

    if not summary:
        path.write_text("# Validation Suite\n\nNo runs.\n", encoding="utf-8")
        return
    lines = [
        "# Validation Suite",
        "",
        "Sprinkle-then-recover summary across ensembles. ",
        "`mean_final_energy` is the optimizer's final energy averaged",
        "across runs; `mean_truth_energy` is the Bombelli energy at",
        "the ground-truth coordinates; `mean_interval_rmse` is the",
        "Lorentz-invariant residual; `mean_mm_dim` is the",
        "Myrheim-Meyer dimension recovered from the causal matrix",
        "alone (no embedding).",
        "",
        "| d_spacetime | n | runs | final_E | truth_E | rmse | mm_dim |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in summary:
        lines.append(
            "| {d_spacetime} | {n} | {runs} | "
            "{mean_final_energy:.4f} | {mean_truth_energy:.4f} | "
            "{mean_interval_rmse:.4f} | {mean_mm_dim:.3f} |".format(**row)
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
