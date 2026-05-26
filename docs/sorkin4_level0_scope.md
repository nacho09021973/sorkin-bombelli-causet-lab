# SORKIN-4 Level 0 — Scope decision

## Decision

Level 0 (L0) is the Dou-Sorkin horizon-molecule program:
count horizon-crossing links, verify N_links ∝ A = 16πM².

L0 is **not** the exterior-only ordering-fraction program.

---

## Why not exterior-only

Schwarzschild exterior is Ricci-flat (R = 0, vacuum solution). This
rules out two of the three candidate observables immediately:

- **Myrheim-Meyer dimension**: inverts the ordering-fraction formula
  against the Minkowski expectation. It estimates d, not M. In the
  Schwarzschild exterior it returns d ≈ 4 for all M.

- **Benincasa-Dowker / Eichhorn curvature estimators**: these recover
  the Ricci scalar R. In vacuum R = 0 everywhere, independent of M.
  Schwarzschild curvature is purely Weyl (tidal), C² = 48M²/r⁶, and
  no standard causal-set estimator projects onto Weyl curvature.

The surviving observable is the **ordering fraction** (causal pairs /
possible pairs), which carries a weak M-signal through the tortoise
coordinate r* = r + 2M ln|r/2M − 1|. The signal exists because the
lightcone shape depends on M. But Schwarzschild exterior is invariant
under the rescaling r → λr, t → λt, M → λM: the causal structure
in rescaled coordinates is identical. A fixed coordinate box [r₁, r₂]
introduces a scale, but the signal is an adimensional ratio M/r_box,
not M itself.

Consequence: with a fixed coordinate box, ordering fraction recovers
the position of the box relative to the horizon (how deep in the
gravitational potential), not M. Discriminating M₁ from M₂ in the
same box requires N ~ several hundred (Poisson noise of order
√N_pairs drowns the signal at N = 12–50) and restricts M < r₁/2.

This is **not** L0 as originally framed. It is a legitimate
exterior-only preliminary (measuring Schwarzschild-vs-Minkowski
intensity), but it must be labeled separately if it is ever done.

---

## Why the horizon breaks the degeneracy

The Schwarzschild horizon at r = r_s = 2M is a physical, intrinsic
feature of the spacetime. It is not coordinate-dependent and cannot
be scaled away: the area A = 4πr_s² = 16πM² is an absolute
quantity (in Planck units). Counting covering relations that cross
the horizon reads this area directly. This is the Dou-Sorkin result:

    N_horizon-links ≈ A / l_Planck²

This observable is:
- purely ordinal (computed from Z alone, no metric input),
- proportional to M² (not to an adimensional ratio),
- derivable from the causal set without reference to embedding
  coordinates.

It requires that the sprinkling include events on both sides of the
horizon and a causal solver that handles interior-exterior pairs.

---

## What L0 requires

1. **Interior event generator**: sample events at r < 2M in
   ingoing Eddington-Finkelstein (IEF) coordinates (v, r, θ, φ).

2. **Causal criterion for interior-interior pairs**: in the interior
   r is the timelike coordinate. For radial pairs, the criterion is
   algebraic in IEF: pair (p, q) with r_p > r_q (p earlier) is
   causal iff v_q ≥ v_p (ingoing null) and v_q ≤ v_outgoing(r_q | p).

3. **Causal criterion for exterior→interior pairs**: radial case is
   algebraic in IEF (v_q ≥ v_p, same condition). Interior→exterior
   is always False (nothing escapes the future horizon).

4. **Horizon-link counter**: covering relations in the transitive
   reduction that have one end in r > 2M and one in r < 2M.

---

## What stays second-priority until the interior solver exists

- Shooting for the 5 undecided exterior pairs: corrects solver
  completeness but does not enable horizon-molecules.

- N and volume study for ordering-fraction discrimination: premature
  until the M signal is established.

- Myrheim-Meyer / chain-abundance analysis of the Schwarzschild
  exterior: valid diagnostics of manifoldlikeness but not M probes.

---

## Denominator convention

`ordering_fraction` = n_causal / (n*(n-1)//2), denominator is always
total possible pairs. Null/undecided entries count as unknown, not as
non-causal. The quantity `solver_coverage` = n_decided / n_possible is
a solver bookkeeping metric and must not be confused with the physical
ordering fraction.

---

## Provenance

- Decision reached in conversation 2026-05-26 after analysis of
  `explore/sorkin4_schwarzschild_benchmark/schwarzschild_minimal_benchmark.json`
  and `explore/sorkin4_kerr_benchmark/kerr_equatorial_scaffold_a0_n12_seed1959.json`.
- Scale-degeneracy argument: Malament theorem + Weyl-only exterior.
- Dou & Sorkin (2003): "Area and the Causal Set".
- He & Rideout (2009): "A Causal Set Black Hole".
