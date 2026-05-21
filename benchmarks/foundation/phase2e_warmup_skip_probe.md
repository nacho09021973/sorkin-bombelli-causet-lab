# Phase 2E Warmup-Skip Probe

Paired comparison of ``with_warmup`` (Phase 2D baseline) vs
``skip_warmup`` (anneal-only) on the same Phase 2D grid.
No changes to ``cones.py``. The diagnostic question is:
is the warmup the primary cause of near-truth destruction,
or does the annealing phase also fail independently?

## Verdict

**WARMUP_IS_PRIMARY_FAILURE_ON_TESTED_GRID**

On this grid, skipping the warmup improves small-noise preservation: 17/18 vs 16/18 with warmup. The warmup phase is the primary contributor to near-truth destruction. Some residual instability in the annealing phase remains, but skipping warmup is a clear improvement.

## Protocol

- families: minkowski only.
- d ∈ {2, 3, 4}, n ∈ {32, 64}, seeds 1959, 1962, 1987.
- optimizer seed: 1987.
- NOISE_SMALL = 0.001, NOISE_MEDIUM = 0.05.
- anneal_limit=10, max_data=4, T₀=100.0, γ=0.9.
- ``with_warmup``: warmup_limit=10 unconditional steps, then anneal.
- ``skip_warmup``: anneal only — ``sim.warmup()`` not called.
- Initialization is identical for the two modes at each ``paired_key``.
- ``paired_key`` = ``target_dim|n|seed|init_label|noise_epsilon``.
- Metrics from last-accepted positions (``rold``/``xold``).
- ``preserved_near_truth``: ``delta_energy ≤ 1e-6``.

## Per-label aggregate by warmup mode

| init | warmup_mode | runs | mean init E | mean final E | mean ΔE | preserved |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| truth | with_warmup | 18 | 0.0000 | 0.0000 | 0.0000 | 18/18 |
| truth | skip_warmup | 18 | 0.0000 | 0.0000 | 0.0000 | 18/18 |
| truth_plus_small_noise | with_warmup | 18 | 0.0047 | 18.9236 | 18.9189 | 16/18 |
| truth_plus_small_noise | skip_warmup | 18 | 0.0047 | 12.1187 | 12.1140 | 17/18 |
| truth_plus_medium_noise | with_warmup | 18 | 11.2341 | 395.8010 | 384.5669 | 0/18 |
| truth_plus_medium_noise | skip_warmup | 18 | 11.2341 | 286.0342 | 274.8001 | 0/18 |
| random_init | with_warmup | 18 | 340.8042 | 405.1642 | 64.3600 | 8/18 |
| random_init | skip_warmup | 18 | 340.8042 | 307.8948 | -32.9094 | 11/18 |

## Paired deltas (skip_warmup − with_warmup)

Negative delta_final_energy means skip_warmup ends lower (better).

| init | mean Δ final E | mean Δ final RMSE | mean Δ final dist | skip pres | with pres |
| --- | ---: | ---: | ---: | ---: | ---: |
| truth | 0.0000 | 0.0000 | 0.0000 | 18/18 | 18/18 |
| truth_plus_small_noise | -6.8049 | 295767.4590 | 57.1209 | 17/18 | 16/18 |
| truth_plus_medium_noise | -109.7668 | -8037875.0926 | -290.3476 | 0/18 | 0/18 |
| random_init | -97.2694 | 59648818.5970 | 669.8424 | 11/18 | 8/18 |

## Five fixed questions

1. **Does skipping warmup preserve truth init?**
   with_warmup: 18/18 preserved, mean final_energy = 0.0000.
   skip_warmup: 18/18 preserved, mean final_energy = 0.0000.
   Truth init has energy = 0; the warmup exits in 0 steps in both modes so both results should be identical.

2. **Does skipping warmup preserve small-noise starts?**
   with_warmup: 16/18 preserved, mean final_energy = 18.9236.
   skip_warmup: 17/18 preserved, mean final_energy = 12.1187.

3. **Does skipping warmup help medium-noise starts?**
   with_warmup: 0/18 preserved, mean final_energy = 395.8010.
   skip_warmup: 0/18 preserved, mean final_energy = 286.0342.

4. **Does random_init benefit or suffer from skipping warmup?**
   with_warmup: 8/18 preserved, mean final_energy = 405.1642.
   skip_warmup: 11/18 preserved, mean final_energy = 307.8948.
   The warmup was designed to help random starts explore before cooling. If skip_warmup hurts random_init, that is the expected trade-off.

5. **What is the dominant failure mode without warmup?**
   17/18 small-noise starts preserved without warmup. Partial improvement: the warmup is a major contributor, but the annealing phase also has some instability near the truth minimum.

## Full per-run table

| d | n | seed | init | warmup | ε | init E | final E | ΔE | preserved |
| :---: | ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | :---: |
| 2 | 32 | 1959 | truth | with_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 2 | 32 | 1959 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 2 | 32 | 1959 | truth_plus_small_noise | with_warmup | 0.0010 | 0.0017 | 0.0009 | -0.0007 | ✓ |
| 2 | 32 | 1959 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0017 | 0.0001 | -0.0016 | ✓ |
| 2 | 32 | 1959 | truth_plus_medium_noise | with_warmup | 0.0500 | 7.6546 | 213.6539 | 205.9994 | ✗ |
| 2 | 32 | 1959 | truth_plus_medium_noise | skip_warmup | 0.0500 | 7.6546 | 190.6723 | 183.0177 | ✗ |
| 2 | 32 | 1959 | random_init | with_warmup | NA | 62.1468 | 231.9320 | 169.7851 | ✗ |
| 2 | 32 | 1959 | random_init | skip_warmup | NA | 62.1468 | 161.4165 | 99.2697 | ✗ |
| 2 | 32 | 1962 | truth | with_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 2 | 32 | 1962 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 2 | 32 | 1962 | truth_plus_small_noise | with_warmup | 0.0010 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 2 | 32 | 1962 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 2 | 32 | 1962 | truth_plus_medium_noise | with_warmup | 0.0500 | 4.8757 | 194.4276 | 189.5519 | ✗ |
| 2 | 32 | 1962 | truth_plus_medium_noise | skip_warmup | 0.0500 | 4.8757 | 147.8539 | 142.9781 | ✗ |
| 2 | 32 | 1962 | random_init | with_warmup | NA | 56.5293 | 195.7065 | 139.1772 | ✗ |
| 2 | 32 | 1962 | random_init | skip_warmup | NA | 56.5293 | 176.7375 | 120.2082 | ✗ |
| 2 | 32 | 1987 | truth | with_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 2 | 32 | 1987 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 2 | 32 | 1987 | truth_plus_small_noise | with_warmup | 0.0010 | 0.0127 | 0.0091 | -0.0037 | ✓ |
| 2 | 32 | 1987 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0127 | 0.0036 | -0.0092 | ✓ |
| 2 | 32 | 1987 | truth_plus_medium_noise | with_warmup | 0.0500 | 7.9057 | 182.9881 | 175.0824 | ✗ |
| 2 | 32 | 1987 | truth_plus_medium_noise | skip_warmup | 0.0500 | 7.9057 | 167.8099 | 159.9042 | ✗ |
| 2 | 32 | 1987 | random_init | with_warmup | NA | 77.2711 | 208.8548 | 131.5837 | ✗ |
| 2 | 32 | 1987 | random_init | skip_warmup | NA | 77.2711 | 150.4007 | 73.1297 | ✗ |
| 2 | 64 | 1959 | truth | with_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 2 | 64 | 1959 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 2 | 64 | 1959 | truth_plus_small_noise | with_warmup | 0.0010 | 0.0060 | 0.0019 | -0.0041 | ✓ |
| 2 | 64 | 1959 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0060 | 0.0011 | -0.0049 | ✓ |
| 2 | 64 | 1959 | truth_plus_medium_noise | with_warmup | 0.0500 | 26.9055 | 774.9770 | 748.0715 | ✗ |
| 2 | 64 | 1959 | truth_plus_medium_noise | skip_warmup | 0.0500 | 26.9055 | 562.6024 | 535.6969 | ✗ |
| 2 | 64 | 1959 | random_init | with_warmup | NA | 279.7231 | 821.8675 | 542.1443 | ✗ |
| 2 | 64 | 1959 | random_init | skip_warmup | NA | 279.7231 | 668.0234 | 388.3003 | ✗ |
| 2 | 64 | 1962 | truth | with_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 2 | 64 | 1962 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 2 | 64 | 1962 | truth_plus_small_noise | with_warmup | 0.0010 | 0.0026 | 0.0014 | -0.0011 | ✓ |
| 2 | 64 | 1962 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0026 | 0.0008 | -0.0018 | ✓ |
| 2 | 64 | 1962 | truth_plus_medium_noise | with_warmup | 0.0500 | 25.5288 | 859.1365 | 833.6078 | ✗ |
| 2 | 64 | 1962 | truth_plus_medium_noise | skip_warmup | 0.0500 | 25.5288 | 584.0654 | 558.5367 | ✗ |
| 2 | 64 | 1962 | random_init | with_warmup | NA | 309.0265 | 711.2924 | 402.2660 | ✗ |
| 2 | 64 | 1962 | random_init | skip_warmup | NA | 309.0265 | 563.1479 | 254.1214 | ✗ |
| 2 | 64 | 1987 | truth | with_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 2 | 64 | 1987 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 2 | 64 | 1987 | truth_plus_small_noise | with_warmup | 0.0010 | 0.0122 | 0.0015 | -0.0107 | ✓ |
| 2 | 64 | 1987 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0122 | 0.0010 | -0.0112 | ✓ |
| 2 | 64 | 1987 | truth_plus_medium_noise | with_warmup | 0.0500 | 23.3229 | 693.9368 | 670.6139 | ✗ |
| 2 | 64 | 1987 | truth_plus_medium_noise | skip_warmup | 0.0500 | 23.3229 | 527.2845 | 503.9616 | ✗ |
| 2 | 64 | 1987 | random_init | with_warmup | NA | 347.0023 | 753.0743 | 406.0720 | ✗ |
| 2 | 64 | 1987 | random_init | skip_warmup | NA | 347.0023 | 533.5615 | 186.5592 | ✗ |
| 3 | 32 | 1959 | truth | with_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 3 | 32 | 1959 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 3 | 32 | 1959 | truth_plus_small_noise | with_warmup | 0.0010 | 0.0058 | 0.0292 | 0.0234 | ✗ |
| 3 | 32 | 1959 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0058 | 0.0048 | -0.0010 | ✓ |
| 3 | 32 | 1959 | truth_plus_medium_noise | with_warmup | 0.0500 | 5.5739 | 140.1931 | 134.6193 | ✗ |
| 3 | 32 | 1959 | truth_plus_medium_noise | skip_warmup | 0.0500 | 5.5739 | 109.0441 | 103.4703 | ✗ |
| 3 | 32 | 1959 | random_init | with_warmup | NA | 177.8374 | 142.1317 | -35.7057 | ✓ |
| 3 | 32 | 1959 | random_init | skip_warmup | NA | 177.8374 | 106.3990 | -71.4384 | ✓ |
| 3 | 32 | 1962 | truth | with_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 3 | 32 | 1962 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 3 | 32 | 1962 | truth_plus_small_noise | with_warmup | 0.0010 | 0.0008 | 0.0002 | -0.0006 | ✓ |
| 3 | 32 | 1962 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0008 | 0.0002 | -0.0006 | ✓ |
| 3 | 32 | 1962 | truth_plus_medium_noise | with_warmup | 0.0500 | 5.0775 | 125.7512 | 120.6737 | ✗ |
| 3 | 32 | 1962 | truth_plus_medium_noise | skip_warmup | 0.0500 | 5.0775 | 105.3374 | 100.2599 | ✗ |
| 3 | 32 | 1962 | random_init | with_warmup | NA | 148.7281 | 151.7462 | 3.0181 | ✗ |
| 3 | 32 | 1962 | random_init | skip_warmup | NA | 148.7281 | 115.1311 | -33.5970 | ✓ |
| 3 | 32 | 1987 | truth | with_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 3 | 32 | 1987 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 3 | 32 | 1987 | truth_plus_small_noise | with_warmup | 0.0010 | 0.0000 | 0.0000 | -0.0000 | ✓ |
| 3 | 32 | 1987 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0000 | 0.0000 | -0.0000 | ✓ |
| 3 | 32 | 1987 | truth_plus_medium_noise | with_warmup | 0.0500 | 5.2624 | 154.5473 | 149.2849 | ✗ |
| 3 | 32 | 1987 | truth_plus_medium_noise | skip_warmup | 0.0500 | 5.2624 | 108.4301 | 103.1677 | ✗ |
| 3 | 32 | 1987 | random_init | with_warmup | NA | 145.2319 | 173.0916 | 27.8597 | ✗ |
| 3 | 32 | 1987 | random_init | skip_warmup | NA | 145.2319 | 126.5424 | -18.6895 | ✓ |
| 3 | 64 | 1959 | truth | with_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 3 | 64 | 1959 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 3 | 64 | 1959 | truth_plus_small_noise | with_warmup | 0.0010 | 0.0078 | 0.0031 | -0.0047 | ✓ |
| 3 | 64 | 1959 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0078 | 0.0020 | -0.0057 | ✓ |
| 3 | 64 | 1959 | truth_plus_medium_noise | with_warmup | 0.0500 | 16.4860 | 589.6979 | 573.2119 | ✗ |
| 3 | 64 | 1959 | truth_plus_medium_noise | skip_warmup | 0.0500 | 16.4860 | 400.7696 | 384.2835 | ✗ |
| 3 | 64 | 1959 | random_init | with_warmup | NA | 614.1638 | 637.2047 | 23.0408 | ✗ |
| 3 | 64 | 1959 | random_init | skip_warmup | NA | 614.1638 | 531.5656 | -82.5982 | ✓ |
| 3 | 64 | 1962 | truth | with_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 3 | 64 | 1962 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 3 | 64 | 1962 | truth_plus_small_noise | with_warmup | 0.0010 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 3 | 64 | 1962 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 3 | 64 | 1962 | truth_plus_medium_noise | with_warmup | 0.0500 | 17.4783 | 571.2987 | 553.8204 | ✗ |
| 3 | 64 | 1962 | truth_plus_medium_noise | skip_warmup | 0.0500 | 17.4783 | 412.0962 | 394.6179 | ✗ |
| 3 | 64 | 1962 | random_init | with_warmup | NA | 591.5364 | 556.7254 | -34.8111 | ✓ |
| 3 | 64 | 1962 | random_init | skip_warmup | NA | 591.5364 | 427.7928 | -163.7436 | ✓ |
| 3 | 64 | 1987 | truth | with_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 3 | 64 | 1987 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 3 | 64 | 1987 | truth_plus_small_noise | with_warmup | 0.0010 | 0.0201 | 340.5726 | 340.5525 | ✗ |
| 3 | 64 | 1987 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0201 | 218.1184 | 218.0983 | ✗ |
| 3 | 64 | 1987 | truth_plus_medium_noise | with_warmup | 0.0500 | 15.8894 | 509.5348 | 493.6454 | ✗ |
| 3 | 64 | 1987 | truth_plus_medium_noise | skip_warmup | 0.0500 | 15.8894 | 451.9343 | 436.0449 | ✗ |
| 3 | 64 | 1987 | random_init | with_warmup | NA | 648.9576 | 589.5010 | -59.4567 | ✓ |
| 3 | 64 | 1987 | random_init | skip_warmup | NA | 648.9576 | 413.0051 | -235.9525 | ✓ |
| 4 | 32 | 1959 | truth | with_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 4 | 32 | 1959 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 4 | 32 | 1959 | truth_plus_small_noise | with_warmup | 0.0010 | 0.0038 | 0.0017 | -0.0020 | ✓ |
| 4 | 32 | 1959 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0038 | 0.0007 | -0.0031 | ✓ |
| 4 | 32 | 1959 | truth_plus_medium_noise | with_warmup | 0.0500 | 3.0763 | 142.9455 | 139.8691 | ✗ |
| 4 | 32 | 1959 | truth_plus_medium_noise | skip_warmup | 0.0500 | 3.0763 | 120.9576 | 117.8812 | ✗ |
| 4 | 32 | 1959 | random_init | with_warmup | NA | 175.7553 | 152.4355 | -23.3198 | ✓ |
| 4 | 32 | 1959 | random_init | skip_warmup | NA | 175.7553 | 125.4376 | -50.3177 | ✓ |
| 4 | 32 | 1962 | truth | with_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 4 | 32 | 1962 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 4 | 32 | 1962 | truth_plus_small_noise | with_warmup | 0.0010 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 4 | 32 | 1962 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 4 | 32 | 1962 | truth_plus_medium_noise | with_warmup | 0.0500 | 3.5958 | 184.1474 | 180.5516 | ✗ |
| 4 | 32 | 1962 | truth_plus_medium_noise | skip_warmup | 0.0500 | 3.5958 | 135.0566 | 131.4608 | ✗ |
| 4 | 32 | 1962 | random_init | with_warmup | NA | 150.9673 | 187.7623 | 36.7950 | ✗ |
| 4 | 32 | 1962 | random_init | skip_warmup | NA | 150.9673 | 154.3684 | 3.4011 | ✗ |
| 4 | 32 | 1987 | truth | with_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 4 | 32 | 1987 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 4 | 32 | 1987 | truth_plus_small_noise | with_warmup | 0.0010 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 4 | 32 | 1987 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 4 | 32 | 1987 | truth_plus_medium_noise | with_warmup | 0.0500 | 2.1370 | 116.2877 | 114.1507 | ✗ |
| 4 | 32 | 1987 | truth_plus_medium_noise | skip_warmup | 0.0500 | 2.1370 | 85.9437 | 83.8067 | ✗ |
| 4 | 32 | 1987 | random_init | with_warmup | NA | 171.7091 | 122.2291 | -49.4800 | ✓ |
| 4 | 32 | 1987 | random_init | skip_warmup | NA | 171.7091 | 126.5390 | -45.1701 | ✓ |
| 4 | 64 | 1959 | truth | with_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 4 | 64 | 1959 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 4 | 64 | 1959 | truth_plus_small_noise | with_warmup | 0.0010 | 0.0053 | 0.0018 | -0.0035 | ✓ |
| 4 | 64 | 1959 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0053 | 0.0033 | -0.0020 | ✓ |
| 4 | 64 | 1959 | truth_plus_medium_noise | with_warmup | 0.0500 | 14.9831 | 715.2261 | 700.2430 | ✗ |
| 4 | 64 | 1959 | truth_plus_medium_noise | skip_warmup | 0.0500 | 14.9831 | 447.0284 | 432.0453 | ✗ |
| 4 | 64 | 1959 | random_init | with_warmup | NA | 695.6267 | 661.6395 | -33.9872 | ✓ |
| 4 | 64 | 1959 | random_init | skip_warmup | NA | 695.6267 | 439.3516 | -256.2751 | ✓ |
| 4 | 64 | 1962 | truth | with_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 4 | 64 | 1962 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 4 | 64 | 1962 | truth_plus_small_noise | with_warmup | 0.0010 | 0.0038 | 0.0012 | -0.0026 | ✓ |
| 4 | 64 | 1962 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0038 | 0.0011 | -0.0027 | ✓ |
| 4 | 64 | 1962 | truth_plus_medium_noise | with_warmup | 0.0500 | 10.3482 | 578.3825 | 568.0343 | ✗ |
| 4 | 64 | 1962 | truth_plus_medium_noise | skip_warmup | 0.0500 | 10.3482 | 362.7359 | 352.3877 | ✗ |
| 4 | 64 | 1962 | random_init | with_warmup | NA | 694.9196 | 608.5696 | -86.3500 | ✓ |
| 4 | 64 | 1962 | random_init | skip_warmup | NA | 694.9196 | 421.0576 | -273.8620 | ✓ |
| 4 | 64 | 1987 | truth | with_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 4 | 64 | 1987 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | ✓ |
| 4 | 64 | 1987 | truth_plus_small_noise | with_warmup | 0.0010 | 0.0025 | 0.0005 | -0.0020 | ✓ |
| 4 | 64 | 1987 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0025 | 0.0003 | -0.0022 | ✓ |
| 4 | 64 | 1987 | truth_plus_medium_noise | with_warmup | 0.0500 | 6.1137 | 377.2862 | 371.1725 | ✗ |
| 4 | 64 | 1987 | truth_plus_medium_noise | skip_warmup | 0.0500 | 6.1137 | 228.9941 | 222.8804 | ✗ |
| 4 | 64 | 1987 | random_init | with_warmup | NA | 787.3426 | 387.1911 | -400.1515 | ✓ |
| 4 | 64 | 1987 | random_init | skip_warmup | NA | 787.3426 | 301.6288 | -485.7138 | ✓ |

Regenerate via `make regen-phase2e`. Source tool:
`tools/build_phase2e_warmup_skip_probe.py`.
