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
- CSV and JSON outputs with `N`, seed, relation count, ordering fraction, link
  count, antisymmetry, transitivity, and undecided-pair count.
- A generic transitive-reduction link counter for the asserted order matrix.

## Not Implemented

- The He & Rideout Schwarzschild pairwise causal-relation algorithm.
- Numerical null-geodesic integration.
- Horizon-crossing molecule counts.
- Any order-only horizon finder.
- Kerr geometry.
- Machine learning, symbolic regression, or spin estimation.

The current matrix is a scaffold: no physical Schwarzschild pairwise relation
is asserted.  `False` entries in the JSON matrix mean "not asserted by this
S4-1 scaffold", not "physically spacelike".

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
