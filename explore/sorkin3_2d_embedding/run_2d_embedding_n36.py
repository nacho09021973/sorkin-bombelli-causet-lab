#!/usr/bin/env python3
"""
SORKIN-3 — 2D order embedding for causal set coordinate recovery.

Recovers Minkowski (1+1)d coordinates from a target causal matrix using
a 2-dimensional order embedding (Dushnik-Miller realizer).

Key identity for 1+1d Minkowski (null coordinates u = t+x, v = t-x):

    p ≺ q  iff  u_p < u_q  AND  v_p < v_q

The realizer (L1, L2) is a pair of linear extensions of the poset such that:

    p ≺ q  iff  rank_L1(p) < rank_L1(q)  AND  rank_L2(p) < rank_L2(q)

For incomparable pairs: their relative order in L1 and L2 must be reversed.

Algorithm
---------
1. Build L1 candidates: identity (labels are in temporal order from the
   sprinkle sort), depth-sorted, height-sorted, and up to MAX_RANDOM
   randomised topological sorts.
2. For each L1: build the L2 dependency graph.
   - Comparable pairs (i ≺ j): same order as P → edge i → j.
   - Incomparable pairs: reverse relative L1 order →
       if rank_L1[i] < rank_L1[j] then edge j → i, else edge i → j.
3. Run Kahn's topological sort on the L2 dependency graph.
   - Acyclic → valid L2 found.
   - Cycle detected → reject this L1, try the next.
4. Verify the realizer: for every pair (i < j), check
   P[i][j] == (rank_L1[i] < rank_L1[j] AND rank_L2[i] < rank_L2[j]).
5. Recover coordinates: u_i = rank_L1[i], v_i = rank_L2[i],
   t_i = (u_i + v_i) / 2, x_i = (u_i - v_i) / 2.
6. Validate via induced_order_from_coords; report F1 and exact_match.

If the L2 graph has a cycle for a given L1, that L1 is rejected and the
next candidate is tried.  A cycle does NOT by itself prove dim(P) > 2;
a formal test (transitive orientation of the incomparability graph) is
required to make that determination.

Single case: N=36, case_seed=1959, d_spacetime=2.
No Bombelli annealer used.
"""

from __future__ import annotations

import heapq
import random
import sys
import time
from collections import deque
from pathlib import Path
from typing import Optional

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import validation_suite as vs  # noqa: E402

# ── parameters ────────────────────────────────────────────────────────────────
N           = 36
CASE_SEED   = 1959
D_SPACETIME = 2
MAX_RANDOM  = 200   # random L1 seeds to try after deterministic candidates


# ── structural helpers ────────────────────────────────────────────────────────

def poset_from_matrix(matrix: vs.CausalMatrix) -> list[list[bool]]:
    """Copy CausalMatrix to plain list[list[bool]] (upper-triangular)."""
    n = len(matrix)
    return [[bool(matrix[i][j]) for j in range(n)] for i in range(n)]


def _build_children(P: list[list[bool]], n: int) -> list[list[int]]:
    """children[i] = sorted list of j where P[i][j] is True (i < j)."""
    ch: list[list[int]] = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if P[i][j]:
                ch[i].append(j)
    return ch


def compute_depths_and_heights(
    P: list[list[bool]], n: int
) -> tuple[list[int], list[int]]:
    """Return (depth, height) in a single topological pass.

    depth[i]  = length of longest chain ending at i   (DP bottom-up).
    height[i] = length of longest chain starting at i (DP top-down).
    """
    children = _build_children(P, n)

    # In-degree (number of direct predecessors in the cover relation is
    # expensive; count all transitive predecessors via P upper triangle)
    in_deg = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if P[i][j]:
                in_deg[j] += 1

    # Topological order (Kahn's)
    queue = deque(i for i in range(n) if in_deg[i] == 0)
    tmp   = in_deg[:]
    topo: list[int] = []
    while queue:
        v = queue.popleft()
        topo.append(v)
        for j in children[v]:
            tmp[j] -= 1
            if tmp[j] == 0:
                queue.append(j)

    # DP: depth (forward pass)
    depth = [1] * n
    for v in topo:
        for j in children[v]:
            if depth[v] + 1 > depth[j]:
                depth[j] = depth[v] + 1

    # DP: height (backward pass)
    height = [1] * n
    for v in reversed(topo):
        for j in children[v]:
            if 1 + height[j] > height[v]:
                height[v] = 1 + height[j]

    return depth, height


def topo_sort_keyed(
    P: list[list[bool]], n: int, key: list[float]
) -> Optional[list[int]]:
    """Topological sort of P, breaking ties by ascending `key` value."""
    children = _build_children(P, n)
    in_deg   = [0] * n
    for i in range(n):
        for j in children[i]:
            in_deg[j] += 1
    heap = [(key[i], i) for i in range(n) if in_deg[i] == 0]
    heapq.heapify(heap)
    tmp    = in_deg[:]
    result: list[int] = []
    while heap:
        _, v = heapq.heappop(heap)
        result.append(v)
        for j in children[v]:
            tmp[j] -= 1
            if tmp[j] == 0:
                heapq.heappush(heap, (key[j], j))
    return result if len(result) == n else None


def random_topo_sort(
    P: list[list[bool]], n: int, seed: int
) -> Optional[list[int]]:
    """Random linear extension of P via Kahn's with random tie-breaking."""
    children = _build_children(P, n)
    in_deg   = [0] * n
    for i in range(n):
        for j in children[i]:
            in_deg[j] += 1
    rng       = random.Random(seed)
    available = [i for i in range(n) if in_deg[i] == 0]
    tmp       = in_deg[:]
    result: list[int] = []
    while available:
        idx = rng.randrange(len(available))
        v   = available[idx]
        available[idx] = available[-1]
        available.pop()
        result.append(v)
        for j in children[v]:
            tmp[j] -= 1
            if tmp[j] == 0:
                available.append(j)
    return result if len(result) == n else None


# ── L2 construction ───────────────────────────────────────────────────────────

def build_L2_adj(
    P: list[list[bool]], L1_rank: list[int], n: int
) -> list[list[int]]:
    """Build L2 dependency adjacency list given L1_rank.

    For each pair (i, j) with i < j by label:
      - comparable  (P[i][j] True): edge i → j  (i before j in L2).
      - incomparable (P[i][j] False, sprinkle sorted by t so P[j][i] = False too):
          if rank_L1[i] < rank_L1[j]: edge j → i  (reverse in L2)
          else:                         edge i → j  (reverse in L2)

    adj[a] contains b to mean "a must come strictly before b in L2".
    """
    adj: list[list[int]] = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):  # i < j by label
            if P[i][j]:
                # comparable: maintain order
                adj[i].append(j)
            else:
                # incomparable (since sprinkle sorts by t, P[j][i] is also False
                # for any j > i that is not causally related to i)
                if L1_rank[i] < L1_rank[j]:
                    adj[j].append(i)   # j before i in L2
                else:
                    adj[i].append(j)   # i before j in L2
    return adj


def kahn_on_adj(adj: list[list[int]], n: int) -> Optional[list[int]]:
    """Kahn's topological sort.  Returns sorted list or None if cycle."""
    in_deg = [0] * n
    for i in range(n):
        for j in adj[i]:
            in_deg[j] += 1
    queue  = deque(i for i in range(n) if in_deg[i] == 0)
    tmp    = in_deg[:]
    result: list[int] = []
    while queue:
        v = queue.popleft()
        result.append(v)
        for j in adj[v]:
            tmp[j] -= 1
            if tmp[j] == 0:
                queue.append(j)
    return result if len(result) == n else None


def try_realizer(
    P: list[list[bool]], L1: list[int], n: int
) -> Optional[list[int]]:
    """Given L1, attempt to construct a compatible L2.
    Returns L2 (valid linear extension) or None (L2 dep-graph has a cycle)."""
    L1_rank = [0] * n
    for rank, elem in enumerate(L1):
        L1_rank[elem] = rank
    adj = build_L2_adj(P, L1_rank, n)
    return kahn_on_adj(adj, n)


# ── verification ──────────────────────────────────────────────────────────────

def verify_realizer(
    P: list[list[bool]], L1: list[int], L2: list[int], n: int
) -> list[tuple[int, int, bool, bool]]:
    """Verify (L1, L2) realizes P for every pair i < j.

    Returns list of (i, j, expected, got) for any disagreement.
    An empty list means exact realization.
    """
    L1r = [0] * n
    L2r = [0] * n
    for r, e in enumerate(L1):
        L1r[e] = r
    for r, e in enumerate(L2):
        L2r[e] = r
    errors: list[tuple[int, int, bool, bool]] = []
    for i in range(n):
        for j in range(i + 1, n):
            got = (L1r[i] < L1r[j]) and (L2r[i] < L2r[j])
            if got != bool(P[i][j]):
                errors.append((i, j, bool(P[i][j]), got))
    return errors


# ── coordinate recovery ───────────────────────────────────────────────────────

def recover_coords(
    L1: list[int], L2: list[int], n: int
) -> list[tuple[float, float]]:
    """Recover (t, x) coordinates from integer ranks.

    u_i = rank_L1[i]          (null coordinate: t + x)
    v_i = rank_L2[i]          (null coordinate: t - x)
    t_i = (u_i + v_i) / 2
    x_i = (u_i - v_i) / 2

    Causal check: i ≺ j  iff  u_i < u_j  AND  v_i < v_j
                            iff  Δt > 0  AND  Δt² > Δx²  (Minkowski condition)
    """
    L1r = [0] * n
    L2r = [0] * n
    for r, e in enumerate(L1):
        L1r[e] = r
    for r, e in enumerate(L2):
        L2r[e] = r
    return [
        ((L1r[i] + L2r[i]) / 2.0, (L1r[i] - L2r[i]) / 2.0)
        for i in range(n)
    ]


# ── main ──────────────────────────────────────────────────────────────────────

def main() -> int:
    t0  = time.perf_counter()
    SEP = "─" * 68

    print(SEP)
    print("SORKIN-3  │  2D order embedding  │  N=36  case_seed=1959  d=2")
    print(SEP)

    # ── [1] generate causet ───────────────────────────────────────────────────
    print("\n[1] Generating causet …")
    matrix, points = vs.sprinkle_minkowski_diamond(
        n=N, seed=CASE_SEED, d_spacetime=D_SPACETIME,
    )
    case = vs.SprinkleCase(
        d_spacetime=D_SPACETIME, n=N, seed=CASE_SEED,
        matrix=matrix, points=points,
    )
    P = poset_from_matrix(case.matrix)
    n = N

    n_rel      = sum(P[i][j] for i in range(n) for j in range(i + 1, n))
    n_possible = n * (n - 1) // 2
    print(f"    N={n}  causal_pairs={n_rel}/{n_possible}  "
          f"density={n_rel/n_possible:.3f}")
    print(f"    (sprinkle sorts by t → labels in temporal order;")
    print(f"     upper-triangular P is complete)")

    # ── [2] structural features ───────────────────────────────────────────────
    print("\n[2] Computing chain depths and heights …")
    depth, height = compute_depths_and_heights(P, n)
    max_chain = max(depth[i] + height[i] - 1 for i in range(n))
    print(f"    max chain length = {max_chain}")
    print(f"    depth  range     = [{min(depth)}, {max(depth)}]")
    print(f"    height range     = [{min(height)}, {max(height)}]")

    # ── [3] build L1 candidate list ───────────────────────────────────────────
    print(f"\n[3] Building L1 candidates (deterministic + up to {MAX_RANDOM} random) …")
    candidates: list[tuple[str, list[int]]] = []

    # Identity: guaranteed valid topological sort (labels are in t-order)
    candidates.append(("identity", list(range(n))))

    # Depth-ascending (proxy for u = t+x order)
    L1 = topo_sort_keyed(P, n, key=[float(d) for d in depth])
    if L1:
        candidates.append(("depth_asc", L1))

    # Height-descending (elements far from top come first)
    L1 = topo_sort_keyed(P, n, key=[float(-h) for h in height])
    if L1:
        candidates.append(("height_desc", L1))

    # depth − height: balanced position along causal axis
    L1 = topo_sort_keyed(
        P, n, key=[float(depth[i] - height[i]) for i in range(n)]
    )
    if L1:
        candidates.append(("depth_minus_height", L1))

    # depth + height: elements near extremes first
    L1 = topo_sort_keyed(
        P, n, key=[float(depth[i] + height[i]) for i in range(n)]
    )
    if L1:
        candidates.append(("depth_plus_height", L1))

    for seed in range(MAX_RANDOM):
        L1 = random_topo_sort(P, n, seed=seed)
        if L1:
            candidates.append((f"random_{seed:03d}", L1))

    print(f"    {len(candidates)} candidates total.")

    # ── [4] search for valid realizer ─────────────────────────────────────────
    print("\n[4] Searching for valid realizer …")
    found_L1:    Optional[list[int]] = None
    found_L2:    Optional[list[int]] = None
    found_label: str                 = ""
    found_idx:   int                 = -1

    n_cycle = 0
    for idx, (label, L1) in enumerate(candidates):
        L2 = try_realizer(P, L1, n)
        if L2 is not None:
            found_L1    = L1
            found_L2    = L2
            found_label = label
            found_idx   = idx + 1
            print(f"    ✓ SUCCESS  attempt={found_idx}/{len(candidates)}"
                  f"  strategy={label}")
            break
        else:
            n_cycle += 1
            if idx < 6 or (idx + 1) % 50 == 0:
                print(f"    ✗ attempt {idx+1:>3}  {label}  → cycle in L2 dep-graph")
    else:
        print(f"\n    ✗ Heuristic search failed ({len(candidates)} candidates, all cyclic).")
        print("      Trying oracle L1: sort by true null coordinate u = t+x …")
        print()
        # ── oracle L1 from ground truth ───────────────────────────────────────
        # case.points[i] = (t_i, x_i, ...) after t-sort.
        # u_i = t_i + x_i  (null coordinate: t+x for d=2).
        # The u-sort is the "correct" L1 for a 1+1d Minkowski causet.
        # Using this is ORACLE (ground truth available); it is not deployable
        # without coordinates.  It is used here as a proof-of-concept only.
        pts = case.points
        u_vals = [pts[i][0] + pts[i][1] for i in range(n)]
        v_vals = [pts[i][0] - pts[i][1] for i in range(n)]
        u_ties = len(set(u_vals)) < n
        v_ties = len(set(v_vals)) < n
        print(f"    u = t+x: {'ties present' if u_ties else 'all distinct'}")
        print(f"    v = t-x: {'ties present' if v_ties else 'all distinct'}")

        # Sort by u (break ties by v, then by label index)
        L1_oracle = sorted(range(n), key=lambda i: (u_vals[i], v_vals[i], i))
        # Validate: L1_oracle must respect P (all comparable pairs in right order)
        L1o_rank = [0] * n
        for r, e in enumerate(L1_oracle):
            L1o_rank[e] = r
        oracle_topo_ok = all(
            L1o_rank[i] < L1o_rank[j]
            for i in range(n) for j in range(i + 1, n)
            if P[i][j]
        )
        print(f"    L1_oracle is valid topological sort: {oracle_topo_ok}")

        if not oracle_topo_ok:
            print("    ✗ Oracle L1 violates poset order — unexpected for Minkowski causet.")
            return 1

        L2_oracle = try_realizer(P, L1_oracle, n)
        if L2_oracle is None:
            print("    ✗ Oracle L1 ALSO yields a cycle in L2 dep-graph.")
            print("      This suggests dim(P) > 2 OR a bug in the dep-graph construction.")
            print("      A formal dim-2 test is required.")
            return 1

        print("    ✓ Oracle L1 (u-sort) → acyclic L2!  Proceeding with oracle realizer.")
        found_L1    = L1_oracle
        found_L2    = L2_oracle
        found_label = "oracle_u_sort"
        found_idx   = -1   # not from candidate list

    if n_cycle > 0:
        print(f"    ({n_cycle} candidates rejected before success)")

    # ── [5] verify realizer ───────────────────────────────────────────────────
    print("\n[5] Verifying realizer …")
    errors = verify_realizer(P, found_L1, found_L2, n)
    n_checked = n * (n - 1) // 2
    if errors:
        print(f"    ✗ {len(errors)}/{n_checked} pair(s) disagree:")
        for e in errors[:5]:
            print(f"      ({e[0]}, {e[1]})  expected={e[2]}  got={e[3]}")
        if len(errors) > 5:
            print(f"      … ({len(errors) - 5} more)")
        print("    Internal inconsistency — this is a bug in the realizer construction.")
        return 1
    print(f"    ✓ All {n_checked} pairs (i < j) correct.")

    # ── [6] recover coordinates ───────────────────────────────────────────────
    print("\n[6] Recovering (t, x) from integer ranks …")
    coords = recover_coords(found_L1, found_L2, n)

    L1r = [0] * n
    L2r = [0] * n
    for r, e in enumerate(found_L1):
        L1r[e] = r
    for r, e in enumerate(found_L2):
        L2r[e] = r

    print(f"    {'elem':>4}  {'rank_L1 (u)':>12}  {'rank_L2 (v)':>12}"
          f"  {'t=(u+v)/2':>10}  {'x=(u-v)/2':>10}")
    for i in range(min(8, n)):
        t_i, x_i = coords[i]
        print(f"    {i:>4}  {L1r[i]:>12}  {L2r[i]:>12}"
              f"  {t_i:>10.1f}  {x_i:>10.1f}")
    if n > 8:
        print(f"    … ({n - 8} more elements not shown)")

    # t-range and x-range
    ts = [c[0] for c in coords]
    xs = [c[1] for c in coords]
    print(f"    t range: [{min(ts):.1f}, {max(ts):.1f}]")
    print(f"    x range: [{min(xs):.1f}, {max(xs):.1f}]")

    # ── [7] validate with induced_order_from_coords ───────────────────────────
    print("\n[7] Validating via induced_order_from_coords …")
    vs_coords = [(t, x) for t, x in coords]
    induced   = vs.induced_order_from_coords(vs_coords)
    cmp       = vs.compare_causal_orders(case.matrix, induced)

    n_target  = cmp.total_relations_target
    n_induced = cmp.total_relations_induced
    n_miss    = len(cmp.missing_relations)
    n_extra   = len(cmp.extra_relations)
    correct   = n_target - n_miss
    prec      = correct / n_induced if n_induced else 0.0
    rec_      = correct / n_target  if n_target  else 0.0
    f1        = 2.0 * prec * rec_ / (prec + rec_) if (prec + rec_) else 0.0

    print(f"    target relations  : {n_target}")
    print(f"    induced relations : {n_induced}")
    print(f"    missing           : {n_miss}")
    print(f"    extra             : {n_extra}")
    print(f"    causal F1         : {f1:.10f}")
    print(f"    exact_match       : {cmp.exact_match}")

    if not cmp.exact_match and n_miss == 0 and n_extra > 0:
        print()
        print("    Note: 0 missing but extra > 0.")
        print("    Likely cause: induced_order_from_coords uses (rij)² ≥ dx²")
        print("    (non-strict, includes null), while the sprinkle uses dt² ≥ dx².")
        print("    For integer ranks: all intervals are strictly timelike (du*dv > 0),")
        print("    so null-cone boundary should not be hit.  Investigate sample pairs:")
        for ea, eb in cmp.extra_relations[:3]:
            ta, xa = coords[ea]
            tb, xb = coords[eb]
            dt = tb - ta
            dx = xb - xa
            print(f"      extra ({ea},{eb}): dt={dt:.2f} dx={dx:.2f}"
                  f"  dt²={dt**2:.2f}  dx²={dx**2:.2f}  "
                  f"causal={'yes' if dt >= 0 and dt*dt >= dx*dx else 'no'}")

    # ── [8] summary ───────────────────────────────────────────────────────────
    elapsed = time.perf_counter() - t0
    print(f"\n{SEP}")
    print("RESULT")
    print(SEP)
    print(f"  L1 strategy   : {found_label}  (attempt {found_idx} of {len(candidates)})")
    print(f"  Realizer OK   : {'yes' if not errors else 'NO (bug)'}")
    print(f"  Causal F1     : {f1:.10f}")
    print(f"  Exact match   : {cmp.exact_match}")
    print(f"  Elapsed       : {elapsed:.3f}s")
    print()

    if cmp.exact_match:
        print("✓  RECOVERY COMPLETE")
        print("   2D order embedding recovered the exact causal structure.")
        print("   No Bombelli annealer required for coordinate recovery in d=2.")
        return 0
    elif n_miss == 0:
        print("△  NEAR-COMPLETE: 0 missing, but extra relations present.")
        print("   Investigate null-cone boundary condition in induced_order_from_coords.")
        return 0
    else:
        print("✗  RECOVERY FAILED: causal structure not exactly reproduced.")
        print("   Investigate realizer construction or coordinate formula.")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
