# Schwarzschild Exterior Turning-Branch Audit

This is a numerical audit of the outgoing one-turn competitor to the direct exterior shooting branch.
It is a falsifiable diagnostic, not a proof of the fastest-geodesic lemma.

- Seed range: 1959..1968
- N: 12
- Audited outgoing direct pairs: 64
- Checked turning branches: 0
- Outcome counts: {'no_turning_solution': 64}
- All checked turning branches later: None
- Min turning minus direct EF-time: None
- Min turning angular range gap: 0.6971976012521455

Scope notes:

- Ingoing pairs are skipped: the turning point is the maximum u, so a turned branch cannot later reach a larger u.
- The genuine competitor is outgoing and requires u2 < u1 < u_turn < 1/(3M), i.e. both endpoints outside the photon sphere.
- The code uses c2=(E/L)^2. Direct no-root branches are on c2 > 1/(27M^2); one-turn branches are on c2 < 1/(27M^2).
