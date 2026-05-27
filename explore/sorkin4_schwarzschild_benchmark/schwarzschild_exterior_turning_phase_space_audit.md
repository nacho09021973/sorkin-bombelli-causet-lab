# Schwarzschild Exterior Turning-Branch Phase-Space Audit

This is the C' numerical diagnostic for the outgoing one-turn competitor.
It checks angular-range overlap before any arrival-time comparison.

- r1 range: 3.1..12.0
- r2 max: 20.0
- step: 0.1
- Grid points: 11205
- Status counts: {'disjoint_ranges': 11205}
- All ranges disjoint: True
- Min angular gap: 0.1446478935306829
- Min-gap row: {'r1': 12.0, 'r2': 12.1, 'u1': 0.08333333333333333, 'u2': 0.08264462809917356, 'phi_direct_max': 0.003893230609996174, 'phi_turning_min': 0.14854112414067905, 'turning_minus_direct_phi_gap': 0.1446478935306829, 'status': 'disjoint_ranges'}

Definitions:

- `phi_direct_max`: direct branch angular reach as `c2 -> 1/(27M^2)+`.
- `phi_turning_min`: one-turn branch angular reach as `u_turn -> u1+`.
- Positive gap means the two angular ranges do not overlap at that grid point.
