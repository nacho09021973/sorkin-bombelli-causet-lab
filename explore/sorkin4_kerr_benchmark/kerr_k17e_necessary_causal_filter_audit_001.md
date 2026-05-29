# S4-KERR-K17E necessary causal filter audit

Necessary / pre-shooting filters applied to K17d selected pairs.
not_excluded is not reachability. rejected_by_filter is not proof of spacelike separation.

K17d CSV rows total: 2159
K17d selected pairs read: 2159

## Per-filter survivor counts

- time_order: 2159
- schwarzschild_radial_time_bound: 2159
- angular_sector: 2159
- radial_potential: 2159
- combined_not_excluded: 2159

## Independent rejection counts (a pair may be counted in multiple reasons)

- rejected_time_order: 0
- rejected_schwarzschild_radial_time: 0
- rejected_angular_sector: 0
- rejected_radial_potential: 0

## radial_time_filter_mode histogram

- not_applied_kerr: 1431
- schwarzschild_hard: 728

## Residual distribution before / after filtering

- before: n=2159, min=0.12546941635537312, median=0.9331903958462267, max=1.9785416634980981
- after:  n=2159, min=0.12546941635537312, median=0.9331903958462267, max=1.9785416634980981
- best_residual_before = 0.12546941635537312
- best_residual_after  = 0.12546941635537312
- hits_le_W_TOL before / after = 0 / 0
- near_hits before / after = 0 / 0

## Verdict: non-selective on K17d-selected pairs

K17E is a negative diagnostic. On the 2159 K17d-selected pairs,
every pair survives the implemented necessary/pre-shooting filters:

- time_order_pass: 2159 / 2159
- Schwarzschild radial-time hard filter for a=0: 728 / 728
- radial_time_filter_mode="not_applied_kerr" for a>0: 1431 / 1431
- angular_sector_admissibility_pass: 2159 / 2159
- radial_potential_admissibility_pass: 2159 / 2159
- combined_not_excluded: 2159 / 2159

The residual distribution is unchanged:

- best_residual_before = 0.12546941635537312
- best_residual_after = 0.12546941635537312
- median_residual_before = 0.9331903958462267
- median_residual_after = 0.9331903958462267
- near_hits_before / after = 0 / 0

Interpretation:
The implemented necessary filters are too weak to concentrate useful
cloud-cloud candidates beyond the existing K17d selection. In particular,
radial_potential_admissibility only checks that radial motion is not
obviously forbidden for some b in the grid; it does not enforce compatibility
with the observed Δt, Δr, and Δφ endpoint data. Therefore K17E does not
justify K18 naive cloud-cloud search.

Operational conclusion:
Do not pursue K17F as another threshold filter over the same K17d residuals.
K18 should not be naive i.i.d. cloud-cloud. It should be designed as a
top-K / structured candidate sandbox, with no causal labels and with explicit
caveats.

## Per (N, seed, spin) survival

| N | seed | spin | n_in | n_out | survival_rate |
|---|---|---|---|---|---|
| 12 | 1959 | 0.00 | 10 | 10 | 1.000 |
| 12 | 1959 | 0.25 | 10 | 10 | 1.000 |
| 12 | 1959 | 0.50 | 10 | 10 | 1.000 |
| 12 | 1960 | 0.00 | 9 | 9 | 1.000 |
| 12 | 1960 | 0.25 | 9 | 9 | 1.000 |
| 12 | 1960 | 0.50 | 9 | 9 | 1.000 |
| 12 | 1961 | 0.00 | 15 | 15 | 1.000 |
| 12 | 1961 | 0.25 | 13 | 13 | 1.000 |
| 12 | 1961 | 0.50 | 13 | 13 | 1.000 |
| 24 | 1959 | 0.00 | 39 | 39 | 1.000 |
| 24 | 1959 | 0.25 | 39 | 39 | 1.000 |
| 24 | 1959 | 0.50 | 38 | 38 | 1.000 |
| 24 | 1960 | 0.00 | 44 | 44 | 1.000 |
| 24 | 1960 | 0.25 | 44 | 44 | 1.000 |
| 24 | 1960 | 0.50 | 43 | 43 | 1.000 |
| 24 | 1961 | 0.00 | 46 | 46 | 1.000 |
| 24 | 1961 | 0.25 | 44 | 44 | 1.000 |
| 24 | 1961 | 0.50 | 44 | 44 | 1.000 |
| 48 | 1959 | 0.00 | 154 | 154 | 1.000 |
| 48 | 1959 | 0.25 | 152 | 152 | 1.000 |
| 48 | 1959 | 0.50 | 150 | 150 | 1.000 |
| 48 | 1960 | 0.00 | 195 | 195 | 1.000 |
| 48 | 1960 | 0.25 | 194 | 194 | 1.000 |
| 48 | 1960 | 0.50 | 191 | 191 | 1.000 |
| 48 | 1961 | 0.00 | 216 | 216 | 1.000 |
| 48 | 1961 | 0.25 | 214 | 214 | 1.000 |
| 48 | 1961 | 0.50 | 214 | 214 | 1.000 |

## Constants used

- MASS = 1.0
- ENERGY = 1.0
- W_TOL = 0.001
- R_TOL = 1e-10
- TIME_TOLERANCE_BAND = 0.05
- PHI_SECTOR_MAX = 1.6
- B_GRID = [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0]
- SECTORS = [-2, -1, 0, 1, 2]
- N_R_SAMPLES = 32

## Defaults documented

- B_GRID is the 7-value K17.B_GRID, not the 5-value K17d.PROBE_B_GRID — chosen so the necessary radial-potential filter does not over-reject.
- TIME_TOLERANCE_BAND reuses K17c (0.05) — same band that entered the K17d selection.
- R_TOL reuses K9.R_MIN_TOL (1e-10) — the same tolerance K17 uses internally for R >= 0.
- N_R_SAMPLES = 32 is a sampling resolution decision, not physics.
- For spin_a > 0, radial_time_filter_mode = not_applied_kerr. No Kerr lower bound invented.

## Caveats

- not_excluded is not reachability.
- rejected_by_filter is not proof of spacelike separation.
- no causal_true/false relations decided.
- no production classifier introduced.
- no global Kerr causal claim introduced.
- no Level-B Hawking/Bekenstein claim introduced.
- Schwarzschild radial-time bound applied as hard filter only for spin_a == 0.
- For spin_a > 0 the radial-time filter is recorded as not_applied_kerr; no Kerr lower bound invented.
- The radial_potential filter is necessary, not sufficient.
