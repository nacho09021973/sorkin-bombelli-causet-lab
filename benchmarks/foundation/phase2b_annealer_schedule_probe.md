# Phase 2B Annealer Schedule Probe

Minkowski-only probe of the historical Bombelli-Sorkin annealer
across a small grid of iteration budgets. Phase 2 left the
annealer running with a deliberately tiny schedule and observed
large energy gaps and large interval residuals. Phase 2B asks:
how much of that is the schedule and how much is structural?

Scope and what this probe does *not* do:

- No new optimizer is introduced. Only ``warmup_limit``,
  ``anneal_limit`` and ``max_data`` are varied; the temperature
  schedule is held fixed.
- Only Minkowski sprinklings are run. KR and corona controls
  have no ground-truth coordinates and are excluded by
  construction; this probe is *not* a manifoldness classifier.
- A low final energy alone is not treated as a successful
  embedding; the diagnostic quantities are ``energy_gap`` and
  ``interval_rmse`` against the known ground truth.

Protocol:

- families: minkowski only.
- target spacetime dimensions: 2, 3, 4.
- sizes: n = 32, 64.
- case seeds: 1959, 1962, 1987 (first three Phase 1B atlas seeds).
- optimizer seed: 1987.
- temperature: initial_temp=100.0, cooling_factor=0.9 (fixed across schedules).
- success criterion (conservative): ``energy_gap <= 1.0`` with finite ``interval_rmse``.

Schedules:

| label | warmup_limit | anneal_limit | max_data | reconfigure budget |
| --- | ---: | ---: | ---: | ---: |
| short | 10 | 10 | 4 | 50 |
| medium | 20 | 20 | 6 | 140 |
| long | 30 | 30 | 10 | 330 |

The ``short`` schedule is exactly the Phase 2 configuration, so
the corresponding rows reproduce the Phase 2 Minkowski numbers
and are the baseline against which the larger budgets are read.

Per-schedule aggregates over all (d, n, seed) cells:

| schedule | runs | mean gap | min gap | max gap | mean RMSE | successes |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| short | 18 | 405.1642 | 122.2291 | 821.8675 | 2526264426.9669 | 0 |
| medium | 18 | 727.9536 | 253.3710 | 1476.0588 | 114348777.9043 | 0 |
| long | 18 | 634.8886 | 227.3569 | 1307.8658 | 899751371.7415 | 0 |

Per-dimension aggregates (mean across seeds and sizes):

| d | schedule | runs | mean gap | mean RMSE |
| :---: | --- | ---: | ---: | ---: |
| 2 | short | 6 | 487.1212 | 21335691.3565 |
| 2 | medium | 6 | 912.3535 | 220711533.8809 |
| 2 | long | 6 | 803.5933 | 508858844.3713 |
| 3 | short | 6 | 375.0668 | 15742201.9043 |
| 3 | medium | 6 | 647.5235 | 1222843.1825 |
| 3 | long | 6 | 579.2974 | 2188320371.4972 |
| 4 | short | 6 | 353.3045 | 7541715387.6400 |
| 4 | medium | 6 | 623.9839 | 121111956.6495 |
| 4 | long | 6 | 521.7750 | 2074899.3561 |

Full per-run table:

| d | n | seed | schedule | final E | truth E | gap | RMSE | success | t(s) |
| :---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | :---: | ---: |
| 2 | 32 | 1959 | short | 231.9320 | 0.0000 | 231.9320 | 296115.4204 | no | 0.0200 |
| 2 | 32 | 1959 | medium | 469.4480 | 0.0000 | 469.4480 | 100597897.9249 | no | 0.0556 |
| 2 | 32 | 1959 | long | 374.1263 | 0.0000 | 374.1263 | 283364.8252 | no | 0.1427 |
| 2 | 32 | 1962 | short | 195.7065 | 0.0000 | 195.7065 | 329162.0181 | no | 0.0181 |
| 2 | 32 | 1962 | medium | 389.3593 | 0.0000 | 389.3593 | 9988.8477 | no | 0.0551 |
| 2 | 32 | 1962 | long | 351.2579 | 0.0000 | 351.2579 | 2973706879.7954 | no | 0.1381 |
| 2 | 32 | 1987 | short | 208.8548 | 0.0000 | 208.8548 | 41453.6818 | no | 0.0190 |
| 2 | 32 | 1987 | medium | 416.9110 | 0.0000 | 416.9110 | 560523910.0672 | no | 0.0519 |
| 2 | 32 | 1987 | long | 382.8861 | 0.0000 | 382.8861 | 5498669.6226 | no | 0.1429 |
| 2 | 64 | 1959 | short | 821.8675 | 0.0000 | 821.8675 | 77567276.6828 | no | 0.1124 |
| 2 | 64 | 1959 | medium | 1361.7732 | 0.0000 | 1361.7732 | 93002880.4968 | no | 0.6135 |
| 2 | 64 | 1959 | long | 1307.8658 | 0.0000 | 1307.8658 | 42228090.9231 | no | 1.8209 |
| 2 | 64 | 1962 | short | 711.2924 | 0.0000 | 711.2924 | 43237555.4647 | no | 0.1395 |
| 2 | 64 | 1962 | medium | 1476.0588 | 0.0000 | 1476.0588 | 263309296.5863 | no | 0.4454 |
| 2 | 64 | 1962 | long | 1222.9199 | 0.0000 | 1222.9199 | 8654279.9765 | no | 2.3115 |
| 2 | 64 | 1987 | short | 753.0743 | 0.0000 | 753.0743 | 6542584.8711 | no | 0.0882 |
| 2 | 64 | 1987 | medium | 1360.5707 | 0.0000 | 1360.5707 | 306825229.3623 | no | 0.4603 |
| 2 | 64 | 1987 | long | 1182.5038 | 0.0000 | 1182.5038 | 22781781.0848 | no | 1.6417 |
| 3 | 32 | 1959 | short | 142.1317 | 0.0000 | 142.1317 | 1541467.7505 | no | 0.0230 |
| 3 | 32 | 1959 | medium | 253.3710 | 0.0000 | 253.3710 | 390429.1266 | no | 0.0612 |
| 3 | 32 | 1959 | long | 227.3569 | 0.0000 | 227.3569 | 2326173.5163 | no | 0.1430 |
| 3 | 32 | 1962 | short | 151.7462 | 0.0000 | 151.7462 | 15365289.6835 | no | 0.0213 |
| 3 | 32 | 1962 | medium | 280.3966 | 0.0000 | 280.3966 | 127088.8928 | no | 0.0598 |
| 3 | 32 | 1962 | long | 255.6875 | 0.0000 | 255.6875 | 2189037.5008 | no | 0.1488 |
| 3 | 32 | 1987 | short | 173.0916 | 0.0000 | 173.0916 | 23971832.3418 | no | 0.0223 |
| 3 | 32 | 1987 | medium | 310.9035 | 0.0000 | 310.9035 | 958507.1969 | no | 0.0602 |
| 3 | 32 | 1987 | long | 276.2730 | 0.0000 | 276.2730 | 786944.6560 | no | 0.1494 |
| 3 | 64 | 1959 | short | 637.2047 | 0.0000 | 637.2047 | 40941762.2346 | no | 0.1244 |
| 3 | 64 | 1959 | medium | 1063.9280 | 0.0000 | 1063.9280 | 4084574.2952 | no | 0.6528 |
| 3 | 64 | 1959 | long | 937.1529 | 0.0000 | 937.1529 | 13115748275.9812 | no | 1.6073 |
| 3 | 64 | 1962 | short | 556.7254 | 0.0000 | 556.7254 | 2716567.3915 | no | 0.1050 |
| 3 | 64 | 1962 | medium | 998.6971 | 0.0000 | 998.6971 | 200244.1987 | no | 0.3756 |
| 3 | 64 | 1962 | long | 873.1917 | 0.0000 | 873.1917 | 4803898.5920 | no | 1.3290 |
| 3 | 64 | 1987 | short | 589.5010 | 0.0000 | 589.5010 | 9916292.0238 | no | 0.1630 |
| 3 | 64 | 1987 | medium | 977.8445 | 0.0000 | 977.8445 | 1576215.3846 | no | 0.3135 |
| 3 | 64 | 1987 | long | 906.1221 | 0.0000 | 906.1221 | 4067898.7370 | no | 1.1459 |
| 4 | 32 | 1959 | short | 152.4355 | 0.0000 | 152.4355 | 44620513074.3680 | no | 0.0315 |
| 4 | 32 | 1959 | medium | 314.0936 | 0.0000 | 314.0936 | 3299673.0685 | no | 0.0840 |
| 4 | 32 | 1959 | long | 292.6916 | 0.0000 | 292.6916 | 3138277.3589 | no | 0.2167 |
| 4 | 32 | 1962 | short | 187.7623 | 0.0000 | 187.7623 | 5869180.0806 | no | 0.0309 |
| 4 | 32 | 1962 | medium | 339.1882 | 0.0000 | 339.1882 | 642042589.8706 | no | 0.0848 |
| 4 | 32 | 1962 | long | 319.3694 | 0.0000 | 319.3694 | 846978.7197 | no | 0.2344 |
| 4 | 32 | 1987 | short | 122.2291 | 0.0000 | 122.2291 | 11741138.9409 | no | 0.0289 |
| 4 | 32 | 1987 | medium | 311.8596 | 0.0000 | 311.8596 | 4779951.4089 | no | 0.0780 |
| 4 | 32 | 1987 | long | 242.6113 | 0.0000 | 242.6113 | 1201753.3772 | no | 0.1952 |
| 4 | 64 | 1959 | short | 661.6395 | 0.0000 | 661.6395 | 17146618.9655 | no | 0.2411 |
| 4 | 64 | 1959 | medium | 1088.9789 | 0.0000 | 1088.9789 | 57864184.1834 | no | 0.9338 |
| 4 | 64 | 1959 | long | 894.0180 | 0.0000 | 894.0180 | 3398333.7707 | no | 4.5924 |
| 4 | 64 | 1962 | short | 608.5696 | 0.0000 | 608.5696 | 592942713.9338 | no | 0.1761 |
| 4 | 64 | 1962 | medium | 966.6654 | 0.0000 | 966.6654 | 10959365.4267 | no | 0.7297 |
| 4 | 64 | 1962 | long | 813.9820 | 0.0000 | 813.9820 | 1277147.8290 | no | 2.6348 |
| 4 | 64 | 1987 | short | 387.1911 | 0.0000 | 387.1911 | 2079599.5512 | no | 0.1943 |
| 4 | 64 | 1987 | medium | 723.1177 | 0.0000 | 723.1177 | 7725975.9388 | no | 0.4975 |
| 4 | 64 | 1987 | long | 567.9778 | 0.0000 | 567.9778 | 2586905.0814 | no | 1.7594 |

Interpretation (conservative, framed as five fixed questions):

1. **Does ``energy_gap`` fall as the budget grows?**
   Mean gap rises (short = 405.164, long = 634.889). The minimum gap across all cells is 122.229 for the short schedule and 227.357 for the long schedule. Adding budget makes things worse on average. The annealer's warmup is pushing the state away from the neighborhood of the truth and the cooling stage is not recovering. The bottleneck is in the energy or move implementation, not in iteration budget. This is *not* a claim that Minkowski sprinklings are non-manifoldlike.

2. **Does ``interval_rmse`` drop coherently with budget?**
   Mean RMSE moves from 2.526e+09 (short) to 8.998e+08 (long), but the per-run table shows that RMSE varies by several orders of magnitude across neighboring cells under any single schedule, so a single ensemble mean is not a stable summary. The conservative reading is that the recovered coordinates are not approaching the ground truth under any of the three budgets tested.

3. **Is there a clear difference between d = 2, 3, 4?**
   - d = 2: gap rises with budget (487.1 -> 803.6).
   - d = 3: gap rises with budget (375.1 -> 579.3).
   - d = 4: gap rises with budget (353.3 -> 521.8).
   The same qualitative behaviour appears at every dimension:
   the gap does not collapse with budget. No dimension is
   distinguished as 'fixed by more iterations'.

4. **Does the failure look like budget/schedule or like
   something more structural in the annealer?**
   Across this small grid the gap is not budget-limited.
   That rules out 'short Phase 2 schedule' as the dominant
   cause of the Phase 2 Minkowski residual. The remaining
   candidates are the energy definition, its parametrization,
   or the historical move-set implementation. This probe
   does *not* discriminate between those candidates; it only
   removes the simplest budget explanation.

5. **Is it still invalid to read the annealer as a
   manifoldness classifier?**
   Yes. Across 54 runs the conservative ``success_flag`` is True in 0 cases.
   No Minkowski case has been recovered to within numerical
   tolerance of the truth at any budget tested, so a low
   ``final_energy`` cannot be cited as a successful embedding.
   The annealer is therefore not a validated manifoldness
   classifier and is never applied to KR or corona causets
   in this probe.

Side remarks:

- Phase 2B does **not** introduce a new optimizer, basin
  hopping, parallel tempering, ML, or PySR-driven search.
  Only the historical annealer's iteration knobs are varied.
- KR and corona causets are excluded by construction; the
  probe must not be cited as a manifoldness classifier.

Regenerate via `make regen-phase2b`. Source tool:
`tools/build_phase2b_annealer_schedule_probe.py`.
