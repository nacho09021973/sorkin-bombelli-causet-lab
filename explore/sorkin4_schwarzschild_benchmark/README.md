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

## S4-2b partial-model robustness sweep

The runner accepts simple CLI parameters:

```bash
python3 explore/sorkin4_schwarzschild_benchmark/run_schwarzschild_minimal_benchmark.py \
  --N 12 \
  --seed 1959 \
  --num-seeds 10 \
  --out-prefix schwarzschild_sweep_n12_10seeds
```

`--num-seeds 1` preserves the single-case behavior and writes detailed CSV/JSON
outputs with events, relation states, the asserted causal matrix, and links.
When `--num-seeds > 1`, the script runs consecutive seeds starting at `--seed`
and writes one CSV row per seed plus an aggregate JSON summary.

Per-seed rows contain:

- `N`, `seed`;
- `true_relations`, `false_relations`, `undecided_pairs`, `decided_pairs`;
- `ordering_fraction_decided`;
- `links`;
- `antisymmetric`, `transitive_true_matrix`;
- `decided_transitivity_no_false_contradictions`.

The aggregate JSON reports min/mean/max values for `undecided_pairs`,
`true_relations`, and `ordering_fraction_decided`, plus a count of failed
internal checks.

This sweep does not validate the full Schwarzschild causal relation.  It only
measures how stable and decisive the current partial He-Rideout radial/bound
model is on small exterior samples.  If `undecided_pairs` remains nonzero, the
generic null-geodesic shooting step is still required before treating the
benchmark as a full Schwarzschild causal-matrix generator.

## S4-3a Generic Schwarzschild shooting plan

This section extracts the generic-pair procedure from He & Rideout (2009)
before implementing it.  It is a mathematical implementation plan, not an
executable integrator.

### Inputs

- Two events in Eddington-Finkelstein Schwarzschild coordinates:
  `E1=(t1,r1,theta1,phi1)` and `E2=(t2,r2,theta2,phi2)`.
- Schwarzschild mass `M`, with horizon radius `2M`.
- The pair is considered in time order, `t1 <= t2`.
- S4-3b should call this only after S4-2 radial exact tests and sufficient
  spacelike/timelike bounds have failed to decide the pair.

### Outputs

- `True` if the fastest null geodesic from `E1` reaches the stationary
  worldline through `E2` by time `t2`.
- `False` if that fastest null geodesic reaches after `t2`, or if the paper's
  no-curve condition is met.
- A numerical diagnostic record for audit: final `c^2`, angular integral,
  time integral, Simpson resolution, root/turning-point status, and convergence
  reason.

### Equations Needed

Angular reduction to the equatorial plane follows He & Rideout Section 2.1.
Rotate coordinates so both events have `vartheta=pi/2`, with

```text
phi2 = arccos(
  cos(theta1) cos(theta2)
  + sin(theta1) sin(theta2) cos(phi1 - phi2_original)
)
```

and use the reduced pair `E1=(t1,r1,pi/2,0)`,
`E2=(t2,r2,pi/2,phi2)` with `phi2 in [0,pi]`.  The worldline `gamma`
through `E2` holds `(r, vartheta, varphi)=(r2, pi/2, phi2)` fixed.  The
fastest relevant null geodesic is the one with `|Delta phi|=phi2`, as proved
in their Appendix A.

Constants of motion and null-geodesic equations follow He & Rideout Section
2.4.  For nonzero angular momentum, define

```text
E = (1 - 2M/r) dt_s/dtau
L = r^2 dvarphi/dtau
c = E/L
u = 1/r
```

The radial/angular equation is

```text
(du/dvarphi)^2 = 2M u^3 - u^2 + c^2
dvarphi/du = +/- (2M u^3 - u^2 + c^2)^(-1/2)
```

The EF-time equation can be integrated either as `dt/dvarphi` or `dt/du`.
For implementation, Appendix B uses the `u` form:

```text
dt/du = ( +/- c (2M u^3 - u^2 + c^2)^(-1/2) - 2M u )
        / (u^2 - 2M u^3)
```

where the sign convention is the same as for `dvarphi/du`: plus when
`dvarphi/du > 0`, minus when `dvarphi/du < 0`.

The shooting parameter is `c^2=(E/L)^2`.  For a candidate `c^2`, integrate
`dvarphi/du` from `u1=1/r1` to `u2=1/r2` along the chosen null-geodesic branch.
The target is

```text
Integral[dvarphi/du] = phi2
```

After finding `c^2` that hits `gamma`, integrate `dt/du` over the same path.
The final causal criterion from Section 2.4 is:

```text
related iff elapsed_time <= t2 - t1
unrelated iff elapsed_time > t2 - t1
```

### Turning Points And Roots

Appendix B centers the numerical risk on the cubic

```text
f(u) = 2M u^3 - u^2 + c^2
```

which appears under the square root.  Its root behavior controls allowed
domains and turning points:

- if `c^2 > 1/(27 M^2)`, there are no nonnegative real roots;
- if `c^2 = 1/(27 M^2)`, there is a double root at `u=1/(3M)` and an
  irrelevant negative root at `u=-1/(6M)`;
- as `c^2` decreases below `1/(27 M^2)`, the double root splits: one positive
  root moves toward the horizon `u=1/(2M)`, and the other moves toward
  asymptotic infinity `u=0`;
- the integrand becomes divergent or imaginary if the integration domain hits
  or crosses the forbidden root interval.

S4-3b must therefore compute the nonnegative real roots of `f(u)` for each
candidate `c^2` and reject or adjust candidates where `u1` or `u2` lies
between the nonnegative roots.  He & Rideout adjust `c^2` upward using a
linear extrapolation of `f(u)` at the roots until both endpoints escape the
root interval.

### S4-3a1 Branch and turning-point convention

He & Rideout's generic decision procedure uses a single direct `u` integral
from `u1=1/r1` to `u2=1/r2`.  The branch is not chosen by enumerating all
possible null geodesics.  Instead, their Appendix A selects the fastest
relevant geodesic: among future-directed null geodesics from `E1` to the
stationary worldline `gamma` through `E2`, the fastest one is the geodesic
with `|Delta phi|=phi2`, i.e. the smallest angular travel without extra
windings.  They then search over `c^2=(E/L)^2` for that direct geodesic.

Physical summary:

- `phi` is taken to increase from `0` to `phi2`, so `dphi >= 0`.
- If `u2 > u1` (`r2 < r1`), the geodesic is ingoing in the paper's sign
  convention: `du > 0`, `dphi/du > 0`, and the plus sign is used in
  `dvarphi/du` and `dt/du`.
- If `u2 < u1` (`r2 > r1`), the geodesic is outgoing: `du < 0`,
  `dphi/du < 0`, and the minus sign is used in `dvarphi/du` and `dt/du`.
- If `u2 == u1`, the generic shooting branch is not the first tool; S4-2's
  angular spacelike/timelike bounds should decide many constant-radius cases,
  and any remaining equal-radius generic case needs an explicit S4-3b policy
  because the paper's Appendix B is phrased as an integral from `u1` to `u2`.

Turning-point convention:

- A radial turning point corresponds to a nonnegative root of
  `f(u)=2M u^3-u^2+c^2`, where `(du/dvarphi)^2=0`.
- Appendix B does not split the integral at a turning point and return along a
  second radial branch.
- Instead, for each candidate `c^2`, the implementation must check the
  nonnegative real roots.  If the direct interval between `u1` and `u2` would
  make the integrand divergent or imaginary, the candidate is invalid for the
  direct fastest-branch search.
- He & Rideout raise `c^2` until the roots do not obstruct the direct
  integration domain.  Once `dvarphi/du` is real and finite on the whole direct
  interval, they integrate `Uphi` with Simpson's rule.
- Therefore S4-3b should not invent an ida-y-vuelta radial path through a
  turning point.  It should implement the paper's direct monotone-`u` branch
  first and treat obstructed candidates as `c^2` search failures/adjustments,
  not as instructions to add another integral segment.

Pseudocode for branch choice, still non-executable:

```text
choose_branch(u1, u2):
  if u2 > u1:
    direction = ingoing
    sign = +1       # dphi/du > 0, du > 0
    interval = [u1, u2]
  else if u2 < u1:
    direction = outgoing
    sign = -1       # dphi/du < 0, du < 0, dphi remains nonnegative
    interval = [u2, u1] for root-obstruction tests
    integrate from u1 to u2 with the negative-sign integrand
  else:
    defer to explicit equal-radius handling

candidate_is_valid(c2, interval):
  roots = nonnegative_real_roots(2M u^3 - u^2 + c2)
  if any root lies in the open interval between u1 and u2:
    return false
  if either endpoint is inside the forbidden interval between positive roots:
    return false
  return true

shoot:
  choose branch from sign(u2-u1)
  adjust c2 upward until candidate_is_valid(c2, interval)
  integrate dphi/du on the direct branch
  update c2 until angular integral reaches phi_target
  integrate dt/du on the same branch and same accepted c2
```

Checks for the future implementation:

- For `r2 < r1`, the selected branch must report `ingoing` and use the plus
  sign in both `dvarphi/du` and `dt/du`.
- For `r2 > r1`, the selected branch must report `outgoing` and use the minus
  sign in both equations.
- `dphi` accumulated along either branch must be nonnegative and should be
  compared against `phi2`.
- A candidate with a nonnegative root inside the direct interval must not be
  integrated as one smooth Simpson interval.
- No implementation should add a reflected turning-point segment unless a
  later note traces that extension to a specific source beyond the current
  He & Rideout Appendix B procedure.

What remains ambiguous:

- The paper gives the sign convention and root-obstruction rule, but not a
  full software-level specification for equal-radius generic pairs where
  `u1 == u2` and S4-2 bounds fail.
- The phrase about checking whether `u1` or `u2` lies between roots is paired
  with the requirement that the integrand be real and finite on the whole
  domain `[u1,u2]`; S4-3b should implement the stronger domain-obstruction
  test and document it.
- This resolves branch choice for the direct fastest-branch He-Rideout
  procedure.  It does not authorize adding multi-turn or lensing-like branches.

### Pseudocode

```text
generic_relation(E1, E2, M):
  require t1 <= t2
  rotate angular coordinates:
    phi_target = arccos(cos(theta1)cos(theta2)
                        + sin(theta1)sin(theta2)cos(phi1-phi2))
  set u1 = 1/r1, u2 = 1/r2

  if S4-2 radial exact tests or sufficient bounds decide the pair:
    return that decision

  initialize c2:
    if either event is behind horizon: c2 = 0
    else: c2 = 1/(27 M^2)

  n = 512 Simpson subintervals
  c2_min = unset
  previous guesses = empty

  loop until convergence or explicit numerical failure:
    roots = nonnegative_real_roots(2M u^3 - u^2 + c2)
    if u1 or u2 lies between nonnegative roots:
      raise c2 using the paper's upward adjustment / root extrapolation
      continue

    phi_integral = composite_simpson(dvarphi_du, u1, u2, c2, n, branch)

    if E2 is behind the horizon and c2 == 0 and phi_integral < phi_target:
      return False  # paper's no-curve-to-gamma condition

    if abs(phi_integral - phi_target) <= 5*epsilon:
      tentatively converged
    else if the two latest c2 guesses differ by <= 5*epsilon:
      tentatively converged
    else:
      choose next c2:
        second guess = c2 + 0.03 if phi_integral > phi_target
        second guess = c2 - 0.005 if phi_integral < phi_target
        later guesses use linear interpolation/extrapolation from last two
      if a proposed c2 is below c2_min:
        c2 = mean(c2_min, previous_c2)
      continue

    verify angular integral at 4n:
      phi_4n = composite_simpson(dvarphi_du, u1, u2, c2, 4n, branch)
      eta = max(phi_4n - phi_target, 5*epsilon)
      if abs(phi_4n - phi_integral) > 8*eta:
        double n and repeat the c2 iteration
      else:
        accept c2

  elapsed_time = composite_simpson(dt_du, u1, u2, c2, accepted_n, branch)
  return elapsed_time <= (t2 - t1)
```

The pseudocode intentionally preserves the paper's update rules and tolerance
structure.  It does not replace them with a generic root solver until S4-3b
decides, with documentation, whether a standard bracketing method is an
equivalent implementation of the same shooting target.

### Cases Already Covered By S4-2

- Angular separation effectively zero: exact radial null criteria.
- Exterior ingoing/outgoing radial spacelike lower bounds.
- Exterior angular spacelike lower bounds using `f(r)=r/sqrt(1-2M/r)`.
- Exterior composed radial-plus-angular null-curve timelike bounds.
- Behind-horizon escape attempts are rejected, though the current benchmark is
  exterior-only.

### Cases Requiring Shooting

- Non-radial pairs that pass both spacelike lower-bound tests but fail the
  composed-null-curve timelike sufficient test.
- Exterior generic pairs where the fastest null geodesic is not represented
  by the S4-2 radial/angular composed bound.
- Later, if the benchmark is extended beyond exterior-only, horizon-crossing
  and interior generic pairs require the same root-aware treatment plus the
  paper's no-curve-to-`gamma` condition.

### Numerical Risks

- The sign/branch of `dvarphi/du` and `dt/du` must match the actual direction
  in `u`; turning points may require piecewise integration rather than a single
  monotone interval.
- Near `c^2=1/(27M^2)`, the double root at `u=1/(3M)` makes the angular
  integrand singular.
- Candidate `c^2` values below the valid minimum can put roots inside the
  integration domain, making the square root imaginary.
- The Appendix B update constants (`+0.03`, `-0.005`, `5 epsilon`,
  `8 eta`, starting `n=512`) are paper-specific and should be reproduced first
  before replacing them.
- Floating-point root classification of a cubic with nearly repeated roots
  needs explicit tolerances and audit logging.
- Simpson integration assumes smooth intervals; if a path has a turning point,
  S4-3b must split the integral at the turning point or document why the
  selected branch avoids it for the exterior cases tested.

### Minimum Checks

- Reproduce the paper's Table 1 angular separations and Table 2 generic-pair
  classifications for at least the exterior generic example `1 6`.
- For S4-2b sweep cases, `undecided_pairs` should decrease when shooting is
  enabled; any remaining undecided pair must have a logged numerical failure
  reason.
- Preserve antisymmetry and the decided-transitivity contradiction check.
- Compare Simpson `n` vs `4n` angular integrals and record the accepted
  resolution.
- Log `c^2`, nonnegative roots of `f(u)`, angular integral, elapsed time, and
  final relation decision for every generic pair.

### Implementation blockers

- Resolved for S4-3b direct He-Rideout implementation: use the monotone direct
  `u1 -> u2` branch, plus sign for ingoing `du>0`, minus sign for outgoing
  `du<0`, and treat roots obstructing the direct interval as invalid
  candidates requiring upward `c^2` adjustment.  Do not split into a
  turning-point return branch in S4-3b.
- Need a robust cubic root helper using only Python standard-library tools or
  existing project dependencies.
- Need validation fixtures from the paper's Table 1/Table 2 values before
  trusting new generic decisions.

## S4-3b Diagnostic: `target_above_direct_max`

After implementing the direct He-Rideout shooting (S4-3a pseudocode → code),
a sweep of 10 seeds (`N=12`, seeds 1959–1968, `--enable-shooting`) gives:

```
undecided_pairs=min/mean/max 0/0.5/1
```

The 5 residual undecided pairs all share the same failure code:
`bracket_failure_type = target_above_direct_max`.

### What it means

The direct angular integral `Integral[dvarphi/du, u1 -> u2]` has a maximum
value over all valid `c^2 >= critical`.  That maximum is attained near
`c^2 = 1/(27M^2)`, the critical value, and decays as `c^2` grows.  When the
required angular separation `phi_target` exceeds this maximum, no valid `c^2`
on the direct branch can produce a null geodesic that rotates by the required
angle on a single monotone-`u` pass from `u1` to `u2`.

The five residual pairs from the sweep have:

| seed | pair | phi_target (rad) | phi_max_sampled (rad) |
|------|------|------------------|-----------------------|
| 1960 | 2–11 | 0.468 | 0.110 |
| 1963 | 0–11 | 0.420 | 0.348 |
| 1964 | 1–7  | 0.640 | 0.401 |
| 1966 | 1–7  | 0.341 | 0.273 |
| 1968 | 5–11 | 0.845 | 0.488 |

In every case `c2_at_phi_max = 1/(27M^2) ≈ 0.0370`, confirming the maximum
angular reach of the direct branch occurs at the critical impact parameter.

### Physical interpretation

These are pairs where the fastest null geodesic — the direct monotone-`u`
branch of He & Rideout Appendix B — cannot subtend the required angle.
Physically, a null geodesic from `r1` to `r2` on a single ingoing/outgoing
pass can only rotate up to some maximum angle set by the curvature geometry
and the radial range.  Pairs with larger angular separation require geodesics
that either:

- wrap more than once around the photon sphere (`|Delta phi| > pi` or
  multi-winding geodesics), or
- are connected by a causal curve that is not null but still within the light
  cone (i.e., by a timelike path, which is always possible for causally
  related points once the relation is established).

He & Rideout's Appendix B searches only the direct fastest-branch null
geodesic.  If that geodesic cannot reach `phi_target`, their procedure returns
no relation decision for the generic pair.  The paper's design handles this by
relying on the composed radial-plus-angular sufficient bound to decide most
pairs first; the generic shooting is a residual step for pairs that the bounds
do not decide.

### Status for this benchmark

For the current exterior-only `N=12` sprinklings, `target_above_direct_max`
failures leave between 0 and 1 pairs undecided per seed.  These pairs are not
decided as false; they are genuinely undecided by the implemented procedure.
The current implementation does not invent additional geodesic branches.  The
correct next step, if tighter coverage is required, is to check whether the
composed null-curve bound can be tightened for large-angular-separation pairs,
or whether He & Rideout's paper provides any additional decision rule for these
cases (outside of Appendix B).

Do not interpret residual `target_above_direct_max` undecided pairs as
non-causal relations.  A `None` return is the correct diagnostic output for
pairs that the algorithm cannot decide.

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
