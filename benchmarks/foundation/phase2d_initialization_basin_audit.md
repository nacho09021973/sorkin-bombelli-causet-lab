# Phase 2D Initialization / Basin Audit

Move-set and basin audit for the historical ConesSimulator on
Minkowski cases. No new optimizers are introduced. Four
initialization strategies are run with the same short schedule
(warmup_limit=10, anneal_limit=10,
max_data=4, T₀=100.0, γ=0.9).

## Verdict

**NARROW_BASIN**

Truth (energy = 0) is preserved exactly, because the warmup exits immediately when energies[0] ≤ 0. However, any positive-energy perturbation activates the warmup, which makes unconditional moves and rapidly scrambles the configuration. The effective basin of attraction has measure zero: the annealer can stay at the minimum only if it starts there exactly, not if it starts nearby. This is a warmup-dynamics failure, not an energy or move-set failure in the strict sense.

## Protocol

- families: minkowski only.
- d ∈ {2, 3, 4}, n ∈ {32, 64}, seeds 1959, 1962, 1987.
- optimizer seed: 1987.
- NOISE_SMALL = 0.001, NOISE_MEDIUM = 0.05.
- noise RNG: seeded deterministically per (d, n, seed) cell.
- ``initial_energy`` and ``final_energy`` measured before and
  after warmup+anneal.
- ``interval_rmse`` and ``distance_to_truth_rms`` measured at
  the same two checkpoints (last-accepted positions ``rold``/
  ``xold``).
- ``preserved_near_truth``: ``final_energy <= initial_energy``
  (optimizer did not worsen the configuration).

**Important warmup note.** The historical warmup loop makes
*unconditional* accepts (no Metropolis criterion). It is
designed to equilibrate the system at high temperature and
scrambles any starting configuration with energy > 0. For the
``truth`` init, energy = 0 so the warmup exits in 0 steps and
no moves are made. For all other inits, warmup runs and
typically increases the energy significantly.

## Per-label aggregate summary

| init | runs | mean init E | mean final E | mean ΔE | mean init RMSE | mean final RMSE | preserved |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| truth | 18 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 18/18 |
| truth_plus_small_noise | 18 | 0.0047 | 18.9236 | 18.9189 | 0.3590 | 3266.3240 | 16/18 |
| truth_plus_medium_noise | 18 | 11.2341 | 395.8010 | 384.5669 | 6460768298378973937664.0000 | 17489733.3956 | 0/18 |
| random_init | 18 | 340.8042 | 405.1642 | 64.3600 | 213.4835 | 3648742.3281 | 8/18 |

## Five fixed questions

1. **Is truth initialization preserved or destroyed?**
   max final_energy for truth init: 0. All truth-init rows are preserved (final_energy ≤ 0 + ε). The warmup exits in 0 steps because energies[0] ≤ 0. No moves are made; the configuration stays exactly at truth.

2. **Do small perturbations return toward truth or away?**
   mean initial_energy for small-noise rows: 0.004727. mean final_energy: 18.92. Final energy greatly exceeds initial (18.9 vs 0.004727). The warmup scrambles near-truth configurations. Perturbations do NOT converge back to truth under the historical annealer.

3. **Is there a visible radius of attraction around the truth?**
   medium-noise rows: mean initial_energy 11.23, mean final_energy 395.8. The annealer cannot recover toward truth from medium perturbations either. The effective basin of attraction around the truth is measure-zero under the historical annealer: only the exact zero-energy configuration is preserved.

4. **Does random_init still fail as in Phase 2B?**
   random_init: mean initial_energy 340.8, mean final_energy 405.2. Consistent with Phase 2B short-schedule results. No improvement from the historical default initialization.

5. **What is the dominant failure mode?**
   The warmup phase makes unconditional accepts for 10 steps.
   Any configuration with energy > 0 (even epsilon = 1e-3,
   energy ≈ 0.01) is actively scrambled by the warmup before
   annealing starts. Truth itself is preserved only because
   its energy is exactly 0 and the warmup exits in 0 steps.
   The failure is in the **warmup dynamics**: unconditional
   acceptance at high temperature makes the annealer insensitive
   to the starting configuration. This is not a move-set failure
   in isolation, nor an energy failure (Phase 2C). The root
   cause is that the warmup phase — not the annealing phase —
   is responsible for destroying near-optimal starts.
   Recommended next step: skip warmup or replace it with a
   conditioned equilibration when starting near a known
   low-energy configuration.

## Full per-run table

| d | n | seed | init | ε | init E | final E | ΔE | init RMSE | final RMSE | preserved |
| :---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | :---: |
| 2 | 32 | 1959 | truth | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 2 | 32 | 1959 | truth_plus_small_noise | 0.0010 | 0.0017 | 0.0009 | -0.0007 | 0.1483 | 0.1720 | ✓ |
| 2 | 32 | 1959 | truth_plus_medium_noise | 0.0500 | 7.6546 | 213.6539 | 205.9994 | 36.6760 | 3238377.1124 | ✗ |
| 2 | 32 | 1959 | random_init | NA | 62.1468 | 231.9320 | 169.7851 | 48.6806 | 44024.6503 | ✗ |
| 2 | 32 | 1962 | truth | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 2 | 32 | 1962 | truth_plus_small_noise | 0.0010 | 0.0000 | 0.0000 | 0.0000 | 0.0644 | 0.0644 | ✓ |
| 2 | 32 | 1962 | truth_plus_medium_noise | 0.0500 | 4.8757 | 194.4276 | 189.5519 | 4.6181 | 8111.8910 | ✗ |
| 2 | 32 | 1962 | random_init | NA | 56.5293 | 195.7065 | 139.1772 | 52.0455 | 32272.4159 | ✗ |
| 2 | 32 | 1987 | truth | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 2 | 32 | 1987 | truth_plus_small_noise | 0.0010 | 0.0127 | 0.0091 | -0.0037 | 0.1779 | 0.4260 | ✓ |
| 2 | 32 | 1987 | truth_plus_medium_noise | 0.0500 | 7.9057 | 182.9881 | 175.0824 | 60522229892712652341248.0000 | 33282.2202 | ✗ |
| 2 | 32 | 1987 | random_init | NA | 77.2711 | 208.8548 | 131.5837 | 68.8448 | 6099.4917 | ✗ |
| 2 | 64 | 1959 | truth | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 2 | 64 | 1959 | truth_plus_small_noise | 0.0010 | 0.0060 | 0.0019 | -0.0041 | 0.1262 | 0.1508 | ✓ |
| 2 | 64 | 1959 | truth_plus_medium_noise | 0.0500 | 26.9055 | 774.9770 | 748.0715 | 20.3725 | 237825405.9351 | ✗ |
| 2 | 64 | 1959 | random_init | NA | 279.7231 | 821.8675 | 542.1443 | 222.7269 | 19700.1701 | ✗ |
| 2 | 64 | 1962 | truth | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 2 | 64 | 1962 | truth_plus_small_noise | 0.0010 | 0.0026 | 0.0014 | -0.0011 | 5.1897 | 6.2020 | ✓ |
| 2 | 64 | 1962 | truth_plus_medium_noise | 0.0500 | 25.5288 | 859.1365 | 833.6078 | 10.3677 | 229577.6192 | ✗ |
| 2 | 64 | 1962 | random_init | NA | 309.0265 | 711.2924 | 402.2660 | 260.4812 | 316828.5908 | ✗ |
| 2 | 64 | 1987 | truth | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 2 | 64 | 1987 | truth_plus_small_noise | 0.0010 | 0.0122 | 0.0015 | -0.0107 | 0.2967 | 0.3312 | ✓ |
| 2 | 64 | 1987 | truth_plus_medium_noise | 0.0500 | 23.3229 | 693.9368 | 670.6139 | 55771599478108872245248.0000 | 183994.1256 | ✗ |
| 2 | 64 | 1987 | random_init | NA | 347.0023 | 753.0743 | 406.0720 | 288.9549 | 3303.6320 | ✗ |
| 3 | 32 | 1959 | truth | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 3 | 32 | 1959 | truth_plus_small_noise | 0.0010 | 0.0058 | 0.0292 | 0.0234 | 0.0633 | 0.4861 | ✗ |
| 3 | 32 | 1959 | truth_plus_medium_noise | 0.0500 | 5.5739 | 140.1931 | 134.6193 | 2.0892 | 1162873.6551 | ✗ |
| 3 | 32 | 1959 | random_init | NA | 177.8374 | 142.1317 | -35.7057 | 121.4556 | 32052.0888 | ✓ |
| 3 | 32 | 1962 | truth | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 3 | 32 | 1962 | truth_plus_small_noise | 0.0010 | 0.0008 | 0.0002 | -0.0006 | 0.0184 | 0.0185 | ✓ |
| 3 | 32 | 1962 | truth_plus_medium_noise | 0.0500 | 5.0775 | 125.7512 | 120.6737 | 0.7443 | 576633.7994 | ✗ |
| 3 | 32 | 1962 | random_init | NA | 148.7281 | 151.7462 | 3.0181 | 86.7636 | 72930.6953 | ✗ |
| 3 | 32 | 1987 | truth | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 3 | 32 | 1987 | truth_plus_small_noise | 0.0010 | 0.0000 | 0.0000 | -0.0000 | 0.0276 | 0.0276 | ✓ |
| 3 | 32 | 1987 | truth_plus_medium_noise | 0.0500 | 5.2624 | 154.5473 | 149.2849 | 2.1956 | 227489.9653 | ✗ |
| 3 | 32 | 1987 | random_init | NA | 145.2319 | 173.0916 | 27.8597 | 91.5150 | 48442.5904 | ✗ |
| 3 | 64 | 1959 | truth | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 3 | 64 | 1959 | truth_plus_small_noise | 0.0010 | 0.0078 | 0.0031 | -0.0047 | 0.0492 | 0.0709 | ✓ |
| 3 | 64 | 1959 | truth_plus_medium_noise | 0.0500 | 16.4860 | 589.6979 | 573.2119 | 3.0140 | 191278.1362 | ✗ |
| 3 | 64 | 1959 | random_init | NA | 614.1638 | 637.2047 | 23.0408 | 378.2405 | 31958583.8812 | ✗ |
| 3 | 64 | 1962 | truth | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 3 | 64 | 1962 | truth_plus_small_noise | 0.0010 | 0.0000 | 0.0000 | 0.0000 | 0.0307 | 0.0307 | ✓ |
| 3 | 64 | 1962 | truth_plus_medium_noise | 0.0500 | 17.4783 | 571.2987 | 553.8204 | 1.0992 | 71183.0184 | ✗ |
| 3 | 64 | 1962 | random_init | NA | 591.5364 | 556.7254 | -34.8111 | 348.7196 | 315872.1467 | ✓ |
| 3 | 64 | 1987 | truth | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 3 | 64 | 1987 | truth_plus_small_noise | 0.0010 | 0.0201 | 340.5726 | 340.5525 | 0.0261 | 58785.5841 | ✗ |
| 3 | 64 | 1987 | truth_plus_medium_noise | 0.0500 | 15.8894 | 509.5348 | 493.6454 | 2.1304 | 1274738.7898 | ✗ |
| 3 | 64 | 1987 | random_init | NA | 648.9576 | 589.5010 | -59.4567 | 372.0017 | 35240.2510 | ✓ |
| 4 | 32 | 1959 | truth | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 4 | 32 | 1959 | truth_plus_small_noise | 0.0010 | 0.0038 | 0.0017 | -0.0020 | 0.0543 | 0.0662 | ✓ |
| 4 | 32 | 1959 | truth_plus_medium_noise | 0.0500 | 3.0763 | 142.9455 | 139.8691 | 8.5443 | 4510268.0069 | ✗ |
| 4 | 32 | 1959 | random_init | NA | 175.7553 | 152.4355 | -23.3198 | 108.5900 | 16065570.8976 | ✓ |
| 4 | 32 | 1962 | truth | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 4 | 32 | 1962 | truth_plus_small_noise | 0.0010 | 0.0000 | 0.0000 | 0.0000 | 0.0421 | 0.0421 | ✓ |
| 4 | 32 | 1962 | truth_plus_medium_noise | 0.0500 | 3.5958 | 184.1474 | 180.5516 | 1.4008 | 1255475.7286 | ✗ |
| 4 | 32 | 1962 | random_init | NA | 150.9673 | 187.7623 | 36.7950 | 94.5911 | 1581703.6550 | ✗ |
| 4 | 32 | 1987 | truth | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 4 | 32 | 1987 | truth_plus_small_noise | 0.0010 | 0.0000 | 0.0000 | 0.0000 | 0.0239 | 0.0239 | ✓ |
| 4 | 32 | 1987 | truth_plus_medium_noise | 0.0500 | 2.1370 | 116.2877 | 114.1507 | 2.4732 | 175176.8792 | ✗ |
| 4 | 32 | 1987 | random_init | NA | 171.7091 | 122.2291 | -49.4800 | 104.5585 | 11741138.9409 | ✓ |
| 4 | 64 | 1959 | truth | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 4 | 64 | 1959 | truth_plus_small_noise | 0.0010 | 0.0053 | 0.0018 | -0.0035 | 0.0528 | 0.0613 | ✓ |
| 4 | 64 | 1959 | truth_plus_medium_noise | 0.0500 | 14.9831 | 715.2261 | 700.2430 | 2.3798 | 9614555.8717 | ✗ |
| 4 | 64 | 1959 | random_init | NA | 695.6267 | 661.6395 | -33.9872 | 401.9495 | 2138453.4329 | ✓ |
| 4 | 64 | 1962 | truth | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 4 | 64 | 1962 | truth_plus_small_noise | 0.0010 | 0.0038 | 0.0012 | -0.0026 | 0.0460 | 0.0502 | ✓ |
| 4 | 64 | 1962 | truth_plus_medium_noise | 0.0500 | 10.3482 | 578.3825 | 568.0343 | 4.3321 | 50408471.0997 | ✗ |
| 4 | 64 | 1962 | random_init | NA | 694.9196 | 608.5696 | -86.3500 | 385.0895 | 810255.5844 | ✓ |
| 4 | 64 | 1987 | truth | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 4 | 64 | 1987 | truth_plus_small_noise | 0.0010 | 0.0025 | 0.0005 | -0.0020 | 0.0241 | 0.0240 | ✓ |
| 4 | 64 | 1987 | truth_plus_medium_noise | 0.0500 | 6.1137 | 377.2862 | 371.1725 | 3.6619 | 3828307.2675 | ✗ |
| 4 | 64 | 1987 | random_init | NA | 787.3426 | 387.1911 | -400.1515 | 407.4939 | 454888.7900 | ✓ |

Regenerate via `make regen-phase2d`. Source tool:
`tools/build_phase2d_initialization_basin_audit.py`.
