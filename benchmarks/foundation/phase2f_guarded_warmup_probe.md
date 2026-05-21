# Phase 2F Guarded-Warmup Probe

Three-way paired comparison of ``legacy_warmup`` (Phase 2E baseline),
``skip_warmup`` (Phase 2E baseline), and ``guarded_warmup`` (new).
No changes to ``cones.py``. The diagnostic question:
can a non-destructive warmup preserve exploratory benefit (helping
random_init) without destroying near-truth starts?

## Verdict

**GUARDED_WARMUP_FIXES_PRIMARY_FAILURE_ON_TESTED_GRID**

guarded_warmup matches or beats skip_warmup on small-noise (18/18 vs 17/18) and improves random_init (12/18 vs 11/18). The guarded warmup preserves exploratory benefit for random starts while not destroying near-truth configurations. The unconditional warmup was the primary failure; the energy guard removes the observed small-noise failure on this grid without sacrificing warmup's utility.

## Protocol

- families: minkowski only.
- d ∈ {2, 3, 4}, n ∈ {32, 64}, seeds 1959, 1962, 1987.
- optimizer seed: 1987.
- NOISE_SMALL = 0.001, NOISE_MEDIUM = 0.05.
- anneal_limit=10, max_data=4, T₀=100.0, γ=0.9.
- warmup_limit=10 for legacy and guarded modes.
- ``paired_key`` = ``target_dim|n|seed|init_label|noise_epsilon``.
- Metrics from last-accepted positions (``rold``/``xold``).
- ``preserved_near_truth``: ``delta_energy ≤ 1e-6``.

### Guarded-warmup details

- GUARD_THRESHOLD = 0.0 (strictly non-worsening,
  applied to pre-normalization ``sim.deltae``).
- On each proposed move: accept iff ``sim.deltae <= 0``.
- On rejection: change flags cleared, ``sim.rave`` restored to ``sim.r``.
- ``rnew``/``xnew`` are NOT explicitly reset because ``reconfigure()``
  always overwrites them from ``rold``/``xold`` at the next call.
- Energy guard is pre-normalization; ``warmup_energy_after`` may differ
  from ``warmup_energy_before`` by a normalization factor even on
  all-accept runs. This is documented, not a defect.

## Per-label aggregate by warmup mode

| init | warmup_mode | runs | mean init E | mean final E | mean ΔE | mean W ΔE | preserved |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| truth | legacy_warmup | 18 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 18/18 |
| truth | skip_warmup | 18 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 18/18 |
| truth | guarded_warmup | 18 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 18/18 |
| truth_plus_small_noise | legacy_warmup | 18 | 0.0047 | 18.9236 | 18.9189 | 0.0564 | 16/18 |
| truth_plus_small_noise | skip_warmup | 18 | 0.0047 | 12.1187 | 12.1140 | 0.0000 | 17/18 |
| truth_plus_small_noise | guarded_warmup | 18 | 0.0047 | 0.0013 | -0.0034 | -0.0023 | 18/18 |
| truth_plus_medium_noise | legacy_warmup | 18 | 11.2341 | 395.8010 | 384.5669 | 1218.2232 | 0/18 |
| truth_plus_medium_noise | skip_warmup | 18 | 11.2341 | 286.0342 | 274.8001 | 0.0000 | 0/18 |
| truth_plus_medium_noise | guarded_warmup | 18 | 11.2341 | 255.3192 | 244.0850 | -0.0921 | 0/18 |
| random_init | legacy_warmup | 18 | 340.8042 | 405.1642 | 64.3600 | 811.7395 | 8/18 |
| random_init | skip_warmup | 18 | 340.8042 | 307.8948 | -32.9094 | 0.0000 | 11/18 |
| random_init | guarded_warmup | 18 | 340.8042 | 271.4181 | -69.3861 | -61.8779 | 12/18 |

## Six fixed questions

1. **Does guarded_warmup preserve truth init?**
   legacy_warmup: 18/18 preserved, mean final E = 0.0000.
   skip_warmup: 18/18 preserved, mean final E = 0.0000.
   guarded_warmup: 18/18 preserved, mean final E = 0.0000.

2. **Does guarded_warmup improve small-noise over skip_warmup?**
   legacy_warmup: 16/18 preserved, mean final E = 18.9236.
   skip_warmup: 17/18 preserved, mean final E = 12.1187.
   guarded_warmup: 18/18 preserved, mean final E = 0.0013.

3. **Does guarded_warmup improve medium-noise over skip_warmup?**
   legacy_warmup: 0/18 preserved, mean final E = 395.8010.
   skip_warmup: 0/18 preserved, mean final E = 286.0342.
   guarded_warmup: 0/18 preserved, mean final E = 255.3192.

4. **Does guarded_warmup help random_init without destroying near-truth?**
   legacy_warmup: 8/18 preserved, mean final E = 405.1642.
   skip_warmup: 11/18 preserved, mean final E = 307.8948.
   guarded_warmup: 12/18 preserved, mean final E = 271.4181.

5. **Is Phase 2E's improvement fully explained by eliminating unconditional accepts?**
   Yes: guarded_warmup (18/18) matches skip_warmup (17/18) on small-noise. The energy-gated warmup is sufficient to avoid the Phase 2D failure.

6. **Is there residual failure attributable to anneal/move-set/cooling?**
   Yes: 1/18 small-noise cases fail even with skip_warmup. The annealing phase has residual instability near the truth minimum independent of the warmup.

## Full per-run table

| d | n | seed | init | warmup | ε | init E | final E | ΔE | W-att | W-acc | W-rej | preserved |
| :---: | ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :---: |
| 2 | 32 | 1959 | truth | legacy_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 2 | 32 | 1959 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 2 | 32 | 1959 | truth | guarded_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 2 | 32 | 1959 | truth_plus_small_noise | legacy_warmup | 0.0010 | 0.0017 | 0.0009 | -0.0007 | 10 | 10 | 0 | ✓ |
| 2 | 32 | 1959 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0017 | 0.0001 | -0.0016 | 0 | 0 | 0 | ✓ |
| 2 | 32 | 1959 | truth_plus_small_noise | guarded_warmup | 0.0010 | 0.0017 | 0.0002 | -0.0015 | 10 | 3 | 7 | ✓ |
| 2 | 32 | 1959 | truth_plus_medium_noise | legacy_warmup | 0.0500 | 7.6546 | 213.6539 | 205.9994 | 10 | 10 | 0 | ✗ |
| 2 | 32 | 1959 | truth_plus_medium_noise | skip_warmup | 0.0500 | 7.6546 | 190.6723 | 183.0177 | 0 | 0 | 0 | ✗ |
| 2 | 32 | 1959 | truth_plus_medium_noise | guarded_warmup | 0.0500 | 7.6546 | 135.6582 | 128.0037 | 10 | 1 | 9 | ✗ |
| 2 | 32 | 1959 | random_init | legacy_warmup | NA | 62.1468 | 231.9320 | 169.7851 | 10 | 10 | 0 | ✗ |
| 2 | 32 | 1959 | random_init | skip_warmup | NA | 62.1468 | 161.4165 | 99.2697 | 0 | 0 | 0 | ✗ |
| 2 | 32 | 1959 | random_init | guarded_warmup | NA | 62.1468 | 138.1541 | 76.0073 | 10 | 2 | 8 | ✗ |
| 2 | 32 | 1962 | truth | legacy_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 2 | 32 | 1962 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 2 | 32 | 1962 | truth | guarded_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 2 | 32 | 1962 | truth_plus_small_noise | legacy_warmup | 0.0010 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 2 | 32 | 1962 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 2 | 32 | 1962 | truth_plus_small_noise | guarded_warmup | 0.0010 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 2 | 32 | 1962 | truth_plus_medium_noise | legacy_warmup | 0.0500 | 4.8757 | 194.4276 | 189.5519 | 10 | 10 | 0 | ✗ |
| 2 | 32 | 1962 | truth_plus_medium_noise | skip_warmup | 0.0500 | 4.8757 | 147.8539 | 142.9781 | 0 | 0 | 0 | ✗ |
| 2 | 32 | 1962 | truth_plus_medium_noise | guarded_warmup | 0.0500 | 4.8757 | 126.6257 | 121.7500 | 10 | 3 | 7 | ✗ |
| 2 | 32 | 1962 | random_init | legacy_warmup | NA | 56.5293 | 195.7065 | 139.1772 | 10 | 10 | 0 | ✗ |
| 2 | 32 | 1962 | random_init | skip_warmup | NA | 56.5293 | 176.7375 | 120.2082 | 0 | 0 | 0 | ✗ |
| 2 | 32 | 1962 | random_init | guarded_warmup | NA | 56.5293 | 135.6960 | 79.1668 | 10 | 1 | 9 | ✗ |
| 2 | 32 | 1987 | truth | legacy_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 2 | 32 | 1987 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 2 | 32 | 1987 | truth | guarded_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 2 | 32 | 1987 | truth_plus_small_noise | legacy_warmup | 0.0010 | 0.0127 | 0.0091 | -0.0037 | 10 | 10 | 0 | ✓ |
| 2 | 32 | 1987 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0127 | 0.0036 | -0.0092 | 0 | 0 | 0 | ✓ |
| 2 | 32 | 1987 | truth_plus_small_noise | guarded_warmup | 0.0010 | 0.0127 | 0.0013 | -0.0115 | 10 | 6 | 4 | ✓ |
| 2 | 32 | 1987 | truth_plus_medium_noise | legacy_warmup | 0.0500 | 7.9057 | 182.9881 | 175.0824 | 10 | 10 | 0 | ✗ |
| 2 | 32 | 1987 | truth_plus_medium_noise | skip_warmup | 0.0500 | 7.9057 | 167.8099 | 159.9042 | 0 | 0 | 0 | ✗ |
| 2 | 32 | 1987 | truth_plus_medium_noise | guarded_warmup | 0.0500 | 7.9057 | 130.7728 | 122.8671 | 10 | 2 | 8 | ✗ |
| 2 | 32 | 1987 | random_init | legacy_warmup | NA | 77.2711 | 208.8548 | 131.5837 | 10 | 10 | 0 | ✗ |
| 2 | 32 | 1987 | random_init | skip_warmup | NA | 77.2711 | 150.4007 | 73.1297 | 0 | 0 | 0 | ✗ |
| 2 | 32 | 1987 | random_init | guarded_warmup | NA | 77.2711 | 137.7986 | 60.5276 | 10 | 2 | 8 | ✗ |
| 2 | 64 | 1959 | truth | legacy_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 2 | 64 | 1959 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 2 | 64 | 1959 | truth | guarded_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 2 | 64 | 1959 | truth_plus_small_noise | legacy_warmup | 0.0010 | 0.0060 | 0.0019 | -0.0041 | 10 | 10 | 0 | ✓ |
| 2 | 64 | 1959 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0060 | 0.0011 | -0.0049 | 0 | 0 | 0 | ✓ |
| 2 | 64 | 1959 | truth_plus_small_noise | guarded_warmup | 0.0010 | 0.0060 | 0.0004 | -0.0056 | 10 | 5 | 5 | ✓ |
| 2 | 64 | 1959 | truth_plus_medium_noise | legacy_warmup | 0.0500 | 26.9055 | 774.9770 | 748.0715 | 10 | 10 | 0 | ✗ |
| 2 | 64 | 1959 | truth_plus_medium_noise | skip_warmup | 0.0500 | 26.9055 | 562.6024 | 535.6969 | 0 | 0 | 0 | ✗ |
| 2 | 64 | 1959 | truth_plus_medium_noise | guarded_warmup | 0.0500 | 26.9055 | 577.9384 | 551.0330 | 10 | 2 | 8 | ✗ |
| 2 | 64 | 1959 | random_init | legacy_warmup | NA | 279.7231 | 821.8675 | 542.1443 | 10 | 10 | 0 | ✗ |
| 2 | 64 | 1959 | random_init | skip_warmup | NA | 279.7231 | 668.0234 | 388.3003 | 0 | 0 | 0 | ✗ |
| 2 | 64 | 1959 | random_init | guarded_warmup | NA | 279.7231 | 488.5965 | 208.8734 | 10 | 1 | 9 | ✗ |
| 2 | 64 | 1962 | truth | legacy_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 2 | 64 | 1962 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 2 | 64 | 1962 | truth | guarded_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 2 | 64 | 1962 | truth_plus_small_noise | legacy_warmup | 0.0010 | 0.0026 | 0.0014 | -0.0011 | 10 | 10 | 0 | ✓ |
| 2 | 64 | 1962 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0026 | 0.0008 | -0.0018 | 0 | 0 | 0 | ✓ |
| 2 | 64 | 1962 | truth_plus_small_noise | guarded_warmup | 0.0010 | 0.0026 | 0.0002 | -0.0023 | 10 | 5 | 5 | ✓ |
| 2 | 64 | 1962 | truth_plus_medium_noise | legacy_warmup | 0.0500 | 25.5288 | 859.1365 | 833.6078 | 10 | 10 | 0 | ✗ |
| 2 | 64 | 1962 | truth_plus_medium_noise | skip_warmup | 0.0500 | 25.5288 | 584.0654 | 558.5367 | 0 | 0 | 0 | ✗ |
| 2 | 64 | 1962 | truth_plus_medium_noise | guarded_warmup | 0.0500 | 25.5288 | 500.4166 | 474.8879 | 10 | 1 | 9 | ✗ |
| 2 | 64 | 1962 | random_init | legacy_warmup | NA | 309.0265 | 711.2924 | 402.2660 | 10 | 10 | 0 | ✗ |
| 2 | 64 | 1962 | random_init | skip_warmup | NA | 309.0265 | 563.1479 | 254.1214 | 0 | 0 | 0 | ✗ |
| 2 | 64 | 1962 | random_init | guarded_warmup | NA | 309.0265 | 472.1448 | 163.1184 | 10 | 2 | 8 | ✗ |
| 2 | 64 | 1987 | truth | legacy_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 2 | 64 | 1987 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 2 | 64 | 1987 | truth | guarded_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 2 | 64 | 1987 | truth_plus_small_noise | legacy_warmup | 0.0010 | 0.0122 | 0.0015 | -0.0107 | 10 | 10 | 0 | ✓ |
| 2 | 64 | 1987 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0122 | 0.0010 | -0.0112 | 0 | 0 | 0 | ✓ |
| 2 | 64 | 1987 | truth_plus_small_noise | guarded_warmup | 0.0010 | 0.0122 | 0.0011 | -0.0111 | 10 | 6 | 4 | ✓ |
| 2 | 64 | 1987 | truth_plus_medium_noise | legacy_warmup | 0.0500 | 23.3229 | 693.9368 | 670.6139 | 10 | 10 | 0 | ✗ |
| 2 | 64 | 1987 | truth_plus_medium_noise | skip_warmup | 0.0500 | 23.3229 | 527.2845 | 503.9616 | 0 | 0 | 0 | ✗ |
| 2 | 64 | 1987 | truth_plus_medium_noise | guarded_warmup | 0.0500 | 23.3229 | 531.3812 | 508.0583 | 10 | 1 | 9 | ✗ |
| 2 | 64 | 1987 | random_init | legacy_warmup | NA | 347.0023 | 753.0743 | 406.0720 | 10 | 10 | 0 | ✗ |
| 2 | 64 | 1987 | random_init | skip_warmup | NA | 347.0023 | 533.5615 | 186.5592 | 0 | 0 | 0 | ✗ |
| 2 | 64 | 1987 | random_init | guarded_warmup | NA | 347.0023 | 461.7232 | 114.7210 | 10 | 2 | 8 | ✗ |
| 3 | 32 | 1959 | truth | legacy_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 3 | 32 | 1959 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 3 | 32 | 1959 | truth | guarded_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 3 | 32 | 1959 | truth_plus_small_noise | legacy_warmup | 0.0010 | 0.0058 | 0.0292 | 0.0234 | 10 | 10 | 0 | ✗ |
| 3 | 32 | 1959 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0058 | 0.0048 | -0.0010 | 0 | 0 | 0 | ✓ |
| 3 | 32 | 1959 | truth_plus_small_noise | guarded_warmup | 0.0010 | 0.0058 | 0.0010 | -0.0048 | 10 | 6 | 4 | ✓ |
| 3 | 32 | 1959 | truth_plus_medium_noise | legacy_warmup | 0.0500 | 5.5739 | 140.1931 | 134.6193 | 10 | 10 | 0 | ✗ |
| 3 | 32 | 1959 | truth_plus_medium_noise | skip_warmup | 0.0500 | 5.5739 | 109.0441 | 103.4703 | 0 | 0 | 0 | ✗ |
| 3 | 32 | 1959 | truth_plus_medium_noise | guarded_warmup | 0.0500 | 5.5739 | 75.0721 | 69.4982 | 10 | 3 | 7 | ✗ |
| 3 | 32 | 1959 | random_init | legacy_warmup | NA | 177.8374 | 142.1317 | -35.7057 | 10 | 10 | 0 | ✓ |
| 3 | 32 | 1959 | random_init | skip_warmup | NA | 177.8374 | 106.3990 | -71.4384 | 0 | 0 | 0 | ✓ |
| 3 | 32 | 1959 | random_init | guarded_warmup | NA | 177.8374 | 105.2643 | -72.5731 | 10 | 9 | 1 | ✓ |
| 3 | 32 | 1962 | truth | legacy_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 3 | 32 | 1962 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 3 | 32 | 1962 | truth | guarded_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 3 | 32 | 1962 | truth_plus_small_noise | legacy_warmup | 0.0010 | 0.0008 | 0.0002 | -0.0006 | 10 | 10 | 0 | ✓ |
| 3 | 32 | 1962 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0008 | 0.0002 | -0.0006 | 0 | 0 | 0 | ✓ |
| 3 | 32 | 1962 | truth_plus_small_noise | guarded_warmup | 0.0010 | 0.0008 | 0.0000 | -0.0008 | 10 | 4 | 6 | ✓ |
| 3 | 32 | 1962 | truth_plus_medium_noise | legacy_warmup | 0.0500 | 5.0775 | 125.7512 | 120.6737 | 10 | 10 | 0 | ✗ |
| 3 | 32 | 1962 | truth_plus_medium_noise | skip_warmup | 0.0500 | 5.0775 | 105.3374 | 100.2599 | 0 | 0 | 0 | ✗ |
| 3 | 32 | 1962 | truth_plus_medium_noise | guarded_warmup | 0.0500 | 5.0775 | 84.0160 | 78.9385 | 10 | 2 | 8 | ✗ |
| 3 | 32 | 1962 | random_init | legacy_warmup | NA | 148.7281 | 151.7462 | 3.0181 | 10 | 10 | 0 | ✗ |
| 3 | 32 | 1962 | random_init | skip_warmup | NA | 148.7281 | 115.1311 | -33.5970 | 0 | 0 | 0 | ✓ |
| 3 | 32 | 1962 | random_init | guarded_warmup | NA | 148.7281 | 108.2280 | -40.5001 | 10 | 6 | 4 | ✓ |
| 3 | 32 | 1987 | truth | legacy_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 3 | 32 | 1987 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 3 | 32 | 1987 | truth | guarded_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 3 | 32 | 1987 | truth_plus_small_noise | legacy_warmup | 0.0010 | 0.0000 | 0.0000 | -0.0000 | 10 | 10 | 0 | ✓ |
| 3 | 32 | 1987 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0000 | 0.0000 | -0.0000 | 0 | 0 | 0 | ✓ |
| 3 | 32 | 1987 | truth_plus_small_noise | guarded_warmup | 0.0010 | 0.0000 | 0.0000 | -0.0000 | 10 | 6 | 4 | ✓ |
| 3 | 32 | 1987 | truth_plus_medium_noise | legacy_warmup | 0.0500 | 5.2624 | 154.5473 | 149.2849 | 10 | 10 | 0 | ✗ |
| 3 | 32 | 1987 | truth_plus_medium_noise | skip_warmup | 0.0500 | 5.2624 | 108.4301 | 103.1677 | 0 | 0 | 0 | ✗ |
| 3 | 32 | 1987 | truth_plus_medium_noise | guarded_warmup | 0.0500 | 5.2624 | 83.5373 | 78.2748 | 10 | 0 | 10 | ✗ |
| 3 | 32 | 1987 | random_init | legacy_warmup | NA | 145.2319 | 173.0916 | 27.8597 | 10 | 10 | 0 | ✗ |
| 3 | 32 | 1987 | random_init | skip_warmup | NA | 145.2319 | 126.5424 | -18.6895 | 0 | 0 | 0 | ✓ |
| 3 | 32 | 1987 | random_init | guarded_warmup | NA | 145.2319 | 116.6746 | -28.5572 | 10 | 8 | 2 | ✓ |
| 3 | 64 | 1959 | truth | legacy_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 3 | 64 | 1959 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 3 | 64 | 1959 | truth | guarded_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 3 | 64 | 1959 | truth_plus_small_noise | legacy_warmup | 0.0010 | 0.0078 | 0.0031 | -0.0047 | 10 | 10 | 0 | ✓ |
| 3 | 64 | 1959 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0078 | 0.0020 | -0.0057 | 0 | 0 | 0 | ✓ |
| 3 | 64 | 1959 | truth_plus_small_noise | guarded_warmup | 0.0010 | 0.0078 | 0.0024 | -0.0053 | 10 | 6 | 4 | ✓ |
| 3 | 64 | 1959 | truth_plus_medium_noise | legacy_warmup | 0.0500 | 16.4860 | 589.6979 | 573.2119 | 10 | 10 | 0 | ✗ |
| 3 | 64 | 1959 | truth_plus_medium_noise | skip_warmup | 0.0500 | 16.4860 | 400.7696 | 384.2835 | 0 | 0 | 0 | ✗ |
| 3 | 64 | 1959 | truth_plus_medium_noise | guarded_warmup | 0.0500 | 16.4860 | 394.1235 | 377.6375 | 10 | 1 | 9 | ✗ |
| 3 | 64 | 1959 | random_init | legacy_warmup | NA | 614.1638 | 637.2047 | 23.0408 | 10 | 10 | 0 | ✗ |
| 3 | 64 | 1959 | random_init | skip_warmup | NA | 614.1638 | 531.5656 | -82.5982 | 0 | 0 | 0 | ✓ |
| 3 | 64 | 1959 | random_init | guarded_warmup | NA | 614.1638 | 449.0967 | -165.0671 | 10 | 7 | 3 | ✓ |
| 3 | 64 | 1962 | truth | legacy_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 3 | 64 | 1962 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 3 | 64 | 1962 | truth | guarded_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 3 | 64 | 1962 | truth_plus_small_noise | legacy_warmup | 0.0010 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 3 | 64 | 1962 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 3 | 64 | 1962 | truth_plus_small_noise | guarded_warmup | 0.0010 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 3 | 64 | 1962 | truth_plus_medium_noise | legacy_warmup | 0.0500 | 17.4783 | 571.2987 | 553.8204 | 10 | 10 | 0 | ✗ |
| 3 | 64 | 1962 | truth_plus_medium_noise | skip_warmup | 0.0500 | 17.4783 | 412.0962 | 394.6179 | 0 | 0 | 0 | ✗ |
| 3 | 64 | 1962 | truth_plus_medium_noise | guarded_warmup | 0.0500 | 17.4783 | 340.7128 | 323.2345 | 10 | 0 | 10 | ✗ |
| 3 | 64 | 1962 | random_init | legacy_warmup | NA | 591.5364 | 556.7254 | -34.8111 | 10 | 10 | 0 | ✓ |
| 3 | 64 | 1962 | random_init | skip_warmup | NA | 591.5364 | 427.7928 | -163.7436 | 0 | 0 | 0 | ✓ |
| 3 | 64 | 1962 | random_init | guarded_warmup | NA | 591.5364 | 390.9566 | -200.5798 | 10 | 6 | 4 | ✓ |
| 3 | 64 | 1987 | truth | legacy_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 3 | 64 | 1987 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 3 | 64 | 1987 | truth | guarded_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 3 | 64 | 1987 | truth_plus_small_noise | legacy_warmup | 0.0010 | 0.0201 | 340.5726 | 340.5525 | 10 | 10 | 0 | ✗ |
| 3 | 64 | 1987 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0201 | 218.1184 | 218.0983 | 0 | 0 | 0 | ✗ |
| 3 | 64 | 1987 | truth_plus_small_noise | guarded_warmup | 0.0010 | 0.0201 | 0.0158 | -0.0043 | 10 | 4 | 6 | ✓ |
| 3 | 64 | 1987 | truth_plus_medium_noise | legacy_warmup | 0.0500 | 15.8894 | 509.5348 | 493.6454 | 10 | 10 | 0 | ✗ |
| 3 | 64 | 1987 | truth_plus_medium_noise | skip_warmup | 0.0500 | 15.8894 | 451.9343 | 436.0449 | 0 | 0 | 0 | ✗ |
| 3 | 64 | 1987 | truth_plus_medium_noise | guarded_warmup | 0.0500 | 15.8894 | 305.1790 | 289.2896 | 10 | 2 | 8 | ✗ |
| 3 | 64 | 1987 | random_init | legacy_warmup | NA | 648.9576 | 589.5010 | -59.4567 | 10 | 10 | 0 | ✓ |
| 3 | 64 | 1987 | random_init | skip_warmup | NA | 648.9576 | 413.0051 | -235.9525 | 0 | 0 | 0 | ✓ |
| 3 | 64 | 1987 | random_init | guarded_warmup | NA | 648.9576 | 372.0618 | -276.8958 | 10 | 8 | 2 | ✓ |
| 4 | 32 | 1959 | truth | legacy_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 4 | 32 | 1959 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 4 | 32 | 1959 | truth | guarded_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 4 | 32 | 1959 | truth_plus_small_noise | legacy_warmup | 0.0010 | 0.0038 | 0.0017 | -0.0020 | 10 | 10 | 0 | ✓ |
| 4 | 32 | 1959 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0038 | 0.0007 | -0.0031 | 0 | 0 | 0 | ✓ |
| 4 | 32 | 1959 | truth_plus_small_noise | guarded_warmup | 0.0010 | 0.0038 | 0.0003 | -0.0034 | 10 | 4 | 6 | ✓ |
| 4 | 32 | 1959 | truth_plus_medium_noise | legacy_warmup | 0.0500 | 3.0763 | 142.9455 | 139.8691 | 10 | 10 | 0 | ✗ |
| 4 | 32 | 1959 | truth_plus_medium_noise | skip_warmup | 0.0500 | 3.0763 | 120.9576 | 117.8812 | 0 | 0 | 0 | ✗ |
| 4 | 32 | 1959 | truth_plus_medium_noise | guarded_warmup | 0.0500 | 3.0763 | 117.6867 | 114.6103 | 10 | 3 | 7 | ✗ |
| 4 | 32 | 1959 | random_init | legacy_warmup | NA | 175.7553 | 152.4355 | -23.3198 | 10 | 10 | 0 | ✓ |
| 4 | 32 | 1959 | random_init | skip_warmup | NA | 175.7553 | 125.4376 | -50.3177 | 0 | 0 | 0 | ✓ |
| 4 | 32 | 1959 | random_init | guarded_warmup | NA | 175.7553 | 131.2435 | -44.5118 | 10 | 9 | 1 | ✓ |
| 4 | 32 | 1962 | truth | legacy_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 4 | 32 | 1962 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 4 | 32 | 1962 | truth | guarded_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 4 | 32 | 1962 | truth_plus_small_noise | legacy_warmup | 0.0010 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 4 | 32 | 1962 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 4 | 32 | 1962 | truth_plus_small_noise | guarded_warmup | 0.0010 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 4 | 32 | 1962 | truth_plus_medium_noise | legacy_warmup | 0.0500 | 3.5958 | 184.1474 | 180.5516 | 10 | 10 | 0 | ✗ |
| 4 | 32 | 1962 | truth_plus_medium_noise | skip_warmup | 0.0500 | 3.5958 | 135.0566 | 131.4608 | 0 | 0 | 0 | ✗ |
| 4 | 32 | 1962 | truth_plus_medium_noise | guarded_warmup | 0.0500 | 3.5958 | 119.7721 | 116.1763 | 10 | 2 | 8 | ✗ |
| 4 | 32 | 1962 | random_init | legacy_warmup | NA | 150.9673 | 187.7623 | 36.7950 | 10 | 10 | 0 | ✗ |
| 4 | 32 | 1962 | random_init | skip_warmup | NA | 150.9673 | 154.3684 | 3.4011 | 0 | 0 | 0 | ✗ |
| 4 | 32 | 1962 | random_init | guarded_warmup | NA | 150.9673 | 127.8952 | -23.0721 | 10 | 8 | 2 | ✓ |
| 4 | 32 | 1987 | truth | legacy_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 4 | 32 | 1987 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 4 | 32 | 1987 | truth | guarded_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 4 | 32 | 1987 | truth_plus_small_noise | legacy_warmup | 0.0010 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 4 | 32 | 1987 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 4 | 32 | 1987 | truth_plus_small_noise | guarded_warmup | 0.0010 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 4 | 32 | 1987 | truth_plus_medium_noise | legacy_warmup | 0.0500 | 2.1370 | 116.2877 | 114.1507 | 10 | 10 | 0 | ✗ |
| 4 | 32 | 1987 | truth_plus_medium_noise | skip_warmup | 0.0500 | 2.1370 | 85.9437 | 83.8067 | 0 | 0 | 0 | ✗ |
| 4 | 32 | 1987 | truth_plus_medium_noise | guarded_warmup | 0.0500 | 2.1370 | 78.3608 | 76.2238 | 10 | 3 | 7 | ✗ |
| 4 | 32 | 1987 | random_init | legacy_warmup | NA | 171.7091 | 122.2291 | -49.4800 | 10 | 10 | 0 | ✓ |
| 4 | 32 | 1987 | random_init | skip_warmup | NA | 171.7091 | 126.5390 | -45.1701 | 0 | 0 | 0 | ✓ |
| 4 | 32 | 1987 | random_init | guarded_warmup | NA | 171.7091 | 112.3089 | -59.4002 | 10 | 8 | 2 | ✓ |
| 4 | 64 | 1959 | truth | legacy_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 4 | 64 | 1959 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 4 | 64 | 1959 | truth | guarded_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 4 | 64 | 1959 | truth_plus_small_noise | legacy_warmup | 0.0010 | 0.0053 | 0.0018 | -0.0035 | 10 | 10 | 0 | ✓ |
| 4 | 64 | 1959 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0053 | 0.0033 | -0.0020 | 0 | 0 | 0 | ✓ |
| 4 | 64 | 1959 | truth_plus_small_noise | guarded_warmup | 0.0010 | 0.0053 | 0.0005 | -0.0048 | 10 | 3 | 7 | ✓ |
| 4 | 64 | 1959 | truth_plus_medium_noise | legacy_warmup | 0.0500 | 14.9831 | 715.2261 | 700.2430 | 10 | 10 | 0 | ✗ |
| 4 | 64 | 1959 | truth_plus_medium_noise | skip_warmup | 0.0500 | 14.9831 | 447.0284 | 432.0453 | 0 | 0 | 0 | ✗ |
| 4 | 64 | 1959 | truth_plus_medium_noise | guarded_warmup | 0.0500 | 14.9831 | 446.4729 | 431.4898 | 10 | 0 | 10 | ✗ |
| 4 | 64 | 1959 | random_init | legacy_warmup | NA | 695.6267 | 661.6395 | -33.9872 | 10 | 10 | 0 | ✓ |
| 4 | 64 | 1959 | random_init | skip_warmup | NA | 695.6267 | 439.3516 | -256.2751 | 0 | 0 | 0 | ✓ |
| 4 | 64 | 1959 | random_init | guarded_warmup | NA | 695.6267 | 430.4733 | -265.1534 | 10 | 9 | 1 | ✓ |
| 4 | 64 | 1962 | truth | legacy_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 4 | 64 | 1962 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 4 | 64 | 1962 | truth | guarded_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 4 | 64 | 1962 | truth_plus_small_noise | legacy_warmup | 0.0010 | 0.0038 | 0.0012 | -0.0026 | 10 | 10 | 0 | ✓ |
| 4 | 64 | 1962 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0038 | 0.0011 | -0.0027 | 0 | 0 | 0 | ✓ |
| 4 | 64 | 1962 | truth_plus_small_noise | guarded_warmup | 0.0010 | 0.0038 | 0.0004 | -0.0034 | 10 | 2 | 8 | ✓ |
| 4 | 64 | 1962 | truth_plus_medium_noise | legacy_warmup | 0.0500 | 10.3482 | 578.3825 | 568.0343 | 10 | 10 | 0 | ✗ |
| 4 | 64 | 1962 | truth_plus_medium_noise | skip_warmup | 0.0500 | 10.3482 | 362.7359 | 352.3877 | 0 | 0 | 0 | ✗ |
| 4 | 64 | 1962 | truth_plus_medium_noise | guarded_warmup | 0.0500 | 10.3482 | 312.0815 | 301.7333 | 10 | 3 | 7 | ✗ |
| 4 | 64 | 1962 | random_init | legacy_warmup | NA | 694.9196 | 608.5696 | -86.3500 | 10 | 10 | 0 | ✓ |
| 4 | 64 | 1962 | random_init | skip_warmup | NA | 694.9196 | 421.0576 | -273.8620 | 0 | 0 | 0 | ✓ |
| 4 | 64 | 1962 | random_init | guarded_warmup | NA | 694.9196 | 362.0641 | -332.8555 | 10 | 8 | 2 | ✓ |
| 4 | 64 | 1987 | truth | legacy_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 4 | 64 | 1987 | truth | skip_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 4 | 64 | 1987 | truth | guarded_warmup | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0 | 0 | 0 | ✓ |
| 4 | 64 | 1987 | truth_plus_small_noise | legacy_warmup | 0.0010 | 0.0025 | 0.0005 | -0.0020 | 10 | 10 | 0 | ✓ |
| 4 | 64 | 1987 | truth_plus_small_noise | skip_warmup | 0.0010 | 0.0025 | 0.0003 | -0.0022 | 0 | 0 | 0 | ✓ |
| 4 | 64 | 1987 | truth_plus_small_noise | guarded_warmup | 0.0010 | 0.0025 | 0.0002 | -0.0023 | 10 | 5 | 5 | ✓ |
| 4 | 64 | 1987 | truth_plus_medium_noise | legacy_warmup | 0.0500 | 6.1137 | 377.2862 | 371.1725 | 10 | 10 | 0 | ✗ |
| 4 | 64 | 1987 | truth_plus_medium_noise | skip_warmup | 0.0500 | 6.1137 | 228.9941 | 222.8804 | 0 | 0 | 0 | ✗ |
| 4 | 64 | 1987 | truth_plus_medium_noise | guarded_warmup | 0.0500 | 6.1137 | 235.9372 | 229.8235 | 10 | 2 | 8 | ✗ |
| 4 | 64 | 1987 | random_init | legacy_warmup | NA | 787.3426 | 387.1911 | -400.1515 | 10 | 10 | 0 | ✓ |
| 4 | 64 | 1987 | random_init | skip_warmup | NA | 787.3426 | 301.6288 | -485.7138 | 0 | 0 | 0 | ✓ |
| 4 | 64 | 1987 | random_init | guarded_warmup | NA | 787.3426 | 345.1451 | -442.1975 | 10 | 10 | 0 | ✓ |

Regenerate via `make regen-phase2f`. Source tool:
`tools/build_phase2f_guarded_warmup_probe.py`.
