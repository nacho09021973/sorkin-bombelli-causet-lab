# Schwarzschild Exterior Turning-Branch Asymptotic Note

This note records the weak-field, nearby-radius expansion suggested by the
turning-branch phase-space audits.  It is a local asymptotic argument for the
outgoing exterior regime

```text
r1 = R,  r2 = R + eps,  R >> M,  eps/R << 1.
```

The notation follows the implementation:

```text
u1 = a = 1/R
u2 = b = 1/(R + eps)
du = a - b = eps / (R (R + eps))
f(u; c2) = 2 M u^3 - u^2 + c2
c_*^2 = 1/(27 M^2)
u_ph = 1/(3M)
```

## Direct Branch

The maximum direct angular reach occurs as the no-root direct branch approaches
the photon-sphere threshold from above:

```text
phi_direct_max = int_b^a du / sqrt(f(u; c_*^2)).
```

For `a = 1/R << 1/M`, the critical cubic is nearly constant on the small
interval:

```text
f(u; c_*^2) = c_*^2 [1 + O((M/R)^2)].
```

Since `c_* = 1/(3 sqrt(3) M)`,

```text
phi_direct_max
  = du / c_* [1 + O(eps/R) + O((M/R)^2)]
  = 3 sqrt(3) M eps / R^2
    [1 + O(eps/R) + O((M/R)^2)].
```

Thus the direct angular maximum is `O(M eps / R^2)`.

## One-Turn Branch

For the outgoing one-turn competitor, let the lower turning root be `u_t`.  The
smallest one-turn angular reach is obtained as `u_t -> a+`.  In that limit the
first leg `a -> u_t` contributes zero, while the return leg tends to

```text
phi_turning_min = int_b^a du / sqrt(f_a(u)),
f_a(u) = 2 M u^3 - u^2 + a^2 - 2 M a^3.
```

The factorization at the endpoint is

```text
f_a(u) = (a - u) [(u + a) - 2M(u^2 + ua + a^2)].
```

Near `u = a`, write `s = a - u`.  Then

```text
f_a(a - s)
  = s [2a(1 - 3Ma) + O(s)].
```

Therefore

```text
phi_turning_min
  = int_0^du ds / sqrt{s [2a(1 - 3Ma) + O(s)]}
  = sqrt(2 du / (a (1 - 3Ma)))
    [1 + O(du/a)].
```

Using `a = 1/R` and `du/a = eps/(R + eps)`,

```text
phi_turning_min
  = sqrt( 2 eps / (R (1 - 3M/R)) )
    [1 + O(eps/R)].
```

Thus the one-turn angular minimum is `O(sqrt(eps/R))`.

## Angular Gap

Combining the two estimates gives

```text
gap = phi_turning_min - phi_direct_max

    = sqrt( 2 eps / (R (1 - 3M/R)) )
      - 3 sqrt(3) M eps / R^2
      + lower-order terms.
```

The ratio of the direct term to the turning term is

```text
O(M sqrt(eps) / R^(3/2)).
```

In the regime `R >> M` and `eps/R << 1`, this ratio is small.  Hence the
leading gap is positive:

```text
gap ~ sqrt(2 eps / R) > 0.
```

For fixed small `eps`, the gap falls from above as `R^(-1/2)`.  Equivalently,

```text
gap * sqrt(R / eps) -> sqrt(2)
```

up to the Schwarzschild correction `(1 - 3M/R)^(-1/2)` and the higher-order
terms above.  This matches the asymptotic audit columns, where
`gap * sqrt(R/eps)` stays close to `sqrt(2)` in the large-`R`, small-`eps`
corner.

## Consequence For The Diagnostic

This asymptotic argument explains why the one-turn outgoing competitor does not
overlap the direct branch in angular target near the numerically dangerous
corner.  The direct branch can sweep only `O(eps/R^2)` radians, while the
one-turn branch has a minimum angular sweep of order `O(sqrt(eps/R))`.

It remains a local asymptotic argument, not a global theorem for all exterior
parameters.
