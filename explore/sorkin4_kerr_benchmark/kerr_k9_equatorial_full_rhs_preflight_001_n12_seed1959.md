# S4-KERR-K9 equatorial full RHS preflight

K9 integrates the full equatorial RHS `t(lambda)`, `r(lambda)`, and `phi(lambda)`.
K9 is not point-to-point shooting.
K9 does not decide causal reachability.
K9 does not classify sprinkled event pairs.
K9 is a preflight for a future geodesic integrator/shooter.
The `b=0` cases are safe radial-flow control cases, not generic Kerr null-geodesic families.
Circular photon orbit cases are drift/hold diagnostics, not stability claims.
Advisory circular cases, when present, are not used as evidence of orbital stability.

- Total cases: 16
- Passed cases: 16
- Failed cases: 0
- Advisory cases: 2
- Global undecided pairs (a>0 control accounting): 66

## Artifact Set

- `kerr_k9_equatorial_full_rhs_preflight_001_n12_seed1959.csv`
- `kerr_k9_equatorial_full_rhs_preflight_001_n12_seed1959.json`
- `kerr_k9_equatorial_full_rhs_preflight_001_n12_seed1959.md`
- `kerr_k9_equatorial_full_rhs_preflight_001_n12_seed1959.png`

## Advisory Cases

- `circular_spin_0.75_pro`
- `circular_spin_0.75_retro`
