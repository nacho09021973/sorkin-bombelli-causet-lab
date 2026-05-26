# Schwarzschild Horizon-Shooting Branch Audit

This is a numerical audit of the direct plunging null branch used by `--horizon-shooting`.
It is evidence for the implementation branch, not a proof of the fastest-geodesic lemma.

- Seed range: 1..40
- Audited horizon-crossing links: 16
- Seeds with crossing links: [1, 4, 8, 16, 19, 21, 23, 28, 31, 32]
- All checks pass: True
- Max angular shooting error: 9.009e-10
- Max Simpson relative drift, 1024 vs 2048: 2.155e-12
- Max exterior raw-vs-regular time-integral error: 6.661e-16

Checks per audited link:

- `c2 > 1/(27M^2)`, so the cubic has no positive root obstruction.
- The rationalized EF-time integrand is finite at the horizon.
- The direct shooting solution reproduces the target angular separation.
- Simpson refinement from 1024 to 2048 intervals is stable.
- Local `phi(c2)` is monotone around the selected shooting parameter.
- The selected null arrives no later than the target event time.

Residual caveat: this does not prove that every possible non-direct or turning branch is slower.
