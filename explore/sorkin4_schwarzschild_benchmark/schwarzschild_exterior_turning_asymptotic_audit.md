# Schwarzschild Exterior Turning-Branch Asymptotic Audit

This probes the weak-field, nearby-radius corner suggested by the phase-space audit.

- R values: [6.0, 8.0, 10.0, 12.0, 16.0, 20.0, 30.0, 50.0, 80.0, 120.0]
- eps values: [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0]
- Grid points: 70
- Status counts: {'disjoint_ranges': 70}
- All ranges disjoint: True
- Min angular gap: 0.004287200542867878
- Min-gap row: {'R': 120.0, 'eps': 0.001, 'eps_over_R': 8.333333333333334e-06, 'r1': 120.0, 'r2': 120.001, 'u1': 0.008333333333333333, 'u2': 0.008333263889467587, 'phi_direct_max': 3.611740192187456e-07, 'phi_turning_min': 0.004287561716887097, 'turning_minus_direct_phi_gap': 0.004287200542867878, 'gap_times_sqrt_R': 0.04696392891754261, 'gap_times_R': 0.5144640651441453, 'status': 'disjoint_ranges'}
- Largest-R / smallest-eps row: {'R': 120.0, 'eps': 0.001, 'eps_over_R': 8.333333333333334e-06, 'r1': 120.0, 'r2': 120.001, 'u1': 0.008333333333333333, 'u2': 0.008333263889467587, 'phi_direct_max': 3.611740192187456e-07, 'phi_turning_min': 0.004287561716887097, 'turning_minus_direct_phi_gap': 0.004287200542867878, 'gap_times_sqrt_R': 0.04696392891754261, 'gap_times_R': 0.5144640651441453, 'status': 'disjoint_ranges'}

Scaling columns:

- `gap_times_sqrt_R` checks whether the gap is roughly `O(R^-1/2)`.
- `gap_times_R` checks whether the gap is roughly `O(R^-1)`.
