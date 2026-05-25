# S4-1 Schwarzschild Minimal Benchmark

## Objective

This directory is the first SORKIN-4 bridge from flat known-truth causal
orders toward black-hole causal-set diagnostics:

Schwarzschild benchmark before Kerr ordinal diagnostics.

The goal is not a new physical result.  The goal is a traceable benchmark
shape for exterior Schwarzschild sprinkling, causal-relation output files,
basic order checks, and later link diagnostics.

## Bibliographic Motivation

The benchmark is motivated primarily by He & Rideout (2009), "A Causal Set
Black Hole", as summarized in the root README section "Help from
bibliography".  That paper gives the relevant Schwarzschild causal-relation
strategy: bounded Schwarzschild sprinklings, cheap sufficient causal/spacelike
checks, and numerical null-geodesic integration for generic pairs.

Dou & Sorkin (2003) and Homšak & Veroni (2024) motivate later link and
horizon-molecule diagnostics once a real Schwarzschild causal matrix exists.

## Implemented

- A small reproducible known-coordinate event set with `N=12` and seed `1959`.
- Exterior-only regime: `r > r_s + margin`, with `r_s = 2.0` and margin `0.35`.
- Bounded coordinate ranges for ingoing Eddington-Finkelstein-like exterior
  coordinates.
- Radial/angular sampling using the `r^2 sin(theta)` volume-factor shape in a
  bounded exterior region.
- CSV and JSON outputs with `N`, seed, true relation count, false relation
  count, ordering fraction over decided pairs, link count, antisymmetry,
  transitivity checks, and undecided-pair count.
- A generic transitive-reduction link counter for the asserted order matrix.
- S4-2 partial causal relation extraction from He & Rideout: exact radial
  criteria plus sufficient exterior spacelike/timelike bounds.

## Not Implemented

- The full He & Rideout Schwarzschild pairwise causal-relation algorithm.
- Numerical null-geodesic shooting/integration for generic pairs.
- Horizon-crossing molecule counts.
- Any order-only horizon finder.
- Kerr geometry.
- Machine learning, symbolic regression, or spin estimation.

The current matrix is partial.  In `relation_states`, `true` means a relation
was decided by the implemented He-Rideout radial/bound tests, `false` means a
non-relation was decided by those tests, and `null` means the pair is generic
and still requires the numerical null-geodesic procedure.  In `causal_matrix`,
`false` includes both decided false pairs and undecided pairs, so
`relation_states` is the authoritative diagnostic field.

## S4-2 He-Rideout causal relation extraction

He & Rideout use Eddington-Finkelstein coordinates
`E=(t,r,theta,phi)` for Schwarzschild.  For a pair with `t1 <= t2`, they rotate
the angular coordinates so both events lie in an equatorial plane with
`phi1=0` and
`phi2 = arccos(cos(theta1) cos(theta2) + sin(theta1) sin(theta2) cos(phi1-phi2))`,
where `phi2` lies in `[0, pi]`.  The causal question is whether the fastest
future-directed null geodesic from `E1` reaches the stationary worldline
through `E2` no later than `t2`.

For radial pairs, He & Rideout give exact criteria.  If `r1 >= r2`, the pair
is related iff `dt >= r1-r2`.  If `r2 >= r1 > 2M`, the pair is related iff
`dt >= r2-r1 + 4M log((r2-2M)/(r1-2M))`.  If a point behind the horizon tries
to reach larger `r`, it is unrelated.  This benchmark implements these radial
criteria.

For non-radial pairs, He & Rideout give cheap sufficient tests before the
generic numerical step.  This benchmark implements:

- radial spacelike lower bounds: if `dt` is below the radial null travel time,
  the pair is decided unrelated;
- angular spacelike lower bounds using `f(r)=r/sqrt(1-2M/r)` outside the
  horizon, with the paper's minimizing `r0`;
- composed radial-plus-angular null-curve timelike bounds, which decide a pair
  related when the constructed null curve reaches the target worldline by
  `t2`.

The inputs needed by the implemented subset are only `M`, `t`, `r`, `theta`,
and `phi` for two events.  The current benchmark is exterior-only:
`r > r_s + margin`, with `r_s = 2M`.  It does not attempt interior or
horizon-crossing generic geodesics.

What remains pending is the generic-pair step from He & Rideout Appendix B:
solve for `c^2=(E/L)^2` so the numerical integral of `dphi/du` reaches the
rotated angular separation, then integrate `dt/du` along that null geodesic
and compare the arrival time with `dt`.  The paper uses composite Simpson
integration, root/domain checks for the cubic `2Mu^3-u^2+c^2`, and iterative
updates of `c^2`.  That procedure is too long for this minimal S4-2 patch, so
generic pairs remain `undecided` rather than being approximated.

Minimum checks for this partial stage:

- antisymmetry of the asserted true-relation matrix;
- no decided-false contradiction of a decided true chain;
- ordering fraction computed only as `true_relations / decided_pairs`;
- explicit `undecided_pairs` count and warning when any generic pair remains.

## Why Schwarzschild Before Kerr

The bibliography section identifies Schwarzschild as the strongest concrete
algorithmic baseline: He & Rideout provide the core pairwise causal-relation
machinery, and Homšak & Veroni scale that line into modern numerical work.
Kerr breaks spherical symmetry and requires additional care, so this project
needs a verified Schwarzschild causal-matrix benchmark before attempting any
rotating-black-hole ordinal diagnostic.

## Condition For Moving Past This Step

Before moving to Kerr, this benchmark needs a real implementation of
`causal_relation_schwarzschild(p, q)` following He & Rideout's Schwarzschild
logic, plus validation on small known-coordinate cases where causal matrices,
links, and basic horizon-local observables can be audited.
