# Phase 4C — Optimizer-seed multi-start probe

**Status:** exploratory probe of optimizer-seed variance.  No new optimizer is introduced.  `cones.py` is not modified.

## Objective

Phase 4A/4B/5 use a single hardcoded `OPTIMIZER_SEED = 1987` inside `ConesSimulator`.  Phase 4C asks whether the Phase 4B `MIXED` and Phase 5 `INSUFFICIENT` outcomes are limited by that single optimizer seed, or robust under K = 3 optimizer-seed perturbation.

Phase 4C does not aim to lower the loss.  It only quantifies how much of the existing pipeline's behaviour is reproducible across different optimizer seeds while holding the causet, the epsilon, the initialization noise, the schedule, and the energy formula unchanged.

## Grid

- Sizes: 32, 48, 64
- Spacetime dims: 2, 3, 4
- Epsilons: 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2
- Causet seeds (10): 1900, 1916, 1923, 1939, 1953, 1973, 1981, 1995, 2003, 2020
- Optimizer seeds (K = 3): 1987, 1990, 1993
- Wall-clock: 249s

## Sample bookkeeping

- Total runs: 2160
- Valid runs: 2160
- Cells: 9

## Global verdict

**OPTIMIZER_SEED_LIMITED**

Decision rule:

- `OPTIMIZER_SEED_LIMITED` if any of: `mean_label_stability < 0.5`, `mean(IQR_K)/mean(loss_K) > 0.1`, `floor_fraction_K < 0.5 * floor_fraction_phase4a`.
- `OPTIMIZER_SEED_ROBUST` if all of: `mean_label_stability > 0.9`, `mean(IQR_K)/mean(loss_K) < 0.01`, `|floor_fraction_K/floor_fraction_phase4a - 1| < 0.1` (or both zero).
- `INCONCLUSIVE` otherwise.

Verdict inputs:

- `mean_label_stability` = 0.444
- `mean_iqr_K` = 0.05682975582
- `mean_loss_K` = 0.06186082723
- `iqr_ratio` (IQR/mean) = 0.9186711263
- `floor_fraction_K` = 0.4743055556
- `floor_fraction_phase4a` = 0.4726273148

## Per-cell label stability

`label_stability_cell` = 1.0 iff all K optimizer-seed curves for a given (n, target_dim) cell receive the same `curve_shape` label.

| n | target_dim | label_stability | curve_shape_per_optimizer_seed |
| ---: | :---: | ---: | --- |
| 32 | 2 | 1.00 | `monotone_decay|monotone_decay|monotone_decay` |
| 32 | 3 | 0.00 | `v_shape|monotone_decay|v_shape` |
| 32 | 4 | 1.00 | `v_shape|v_shape|v_shape` |
| 48 | 2 | 0.00 | `monotone_decay|noisy|noisy` |
| 48 | 3 | 0.00 | `v_shape|v_shape|monotone_decay` |
| 48 | 4 | 0.00 | `noisy|monotone_decay|noisy` |
| 64 | 2 | 1.00 | `noisy|noisy|noisy` |
| 64 | 3 | 0.00 | `monotone_decay|v_shape|noisy` |
| 64 | 4 | 1.00 | `monotone_decay|monotone_decay|monotone_decay` |

## Per-cell-epsilon statistics (K runs combined with causet seeds)

`mean_loss_K`, `std_loss_K`, `IQR_loss_K`, `floor_K` are computed over K * len(causet_seeds) valid rows per (n, target_dim, epsilon).  Phase 4A baseline columns are NaN for n = 48 (no baseline).

| n | target_dim | epsilon | n_valid | mean_loss_K | std_loss_K | IQR_loss_K | floor_K | phase4a_floor | Δ_floor |
| ---: | :---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 32 | 2 | 0.01 | 30 | 0.3140915083 | 0.1549660608 | 0.2323679013 | 0 | 0 | 0 |
| 32 | 2 | 0.02 | 30 | 0.154320764 | 0.1054125009 | 0.1121841487 | 0 | 0 | 0 |
| 32 | 2 | 0.04 | 30 | 0.0512017406 | 0.06486774866 | 0.06355727049 | 0.3 | 0.2 | 0.1 |
| 32 | 2 | 0.06 | 30 | 0.006271576753 | 0.01984316278 | 0.0007547495971 | 0.7 | 0.7 | 0 |
| 32 | 2 | 0.08 | 30 | 0.01664754534 | 0.04177478267 | 0 | 0.7333333333 | 0.5 | 0.2333333333 |
| 32 | 2 | 0.1 | 30 | 0.01042466442 | 0.03411256787 | 0 | 0.8 | 0.6 | 0.2 |
| 32 | 2 | 0.15 | 30 | 0.002687394879 | 0.007719034715 | 0 | 0.8666666667 | 0.9 | -0.03333333333 |
| 32 | 2 | 0.2 | 30 | 0.01180834488 | 0.02768363659 | 0 | 0.7666666667 | 0.7 | 0.06666666667 |
| 32 | 3 | 0.01 | 30 | 0.2942485311 | 0.1804293811 | 0.1760394862 | 0 | 0 | 0 |
| 32 | 3 | 0.02 | 30 | 0.1824275354 | 0.1097039351 | 0.1244319416 | 0.03333333333 | 0.1 | -0.06666666667 |
| 32 | 3 | 0.04 | 30 | 0.09293109877 | 0.08858005135 | 0.1184725529 | 0.1666666667 | 0.2 | -0.03333333333 |
| 32 | 3 | 0.06 | 30 | 0.06158432669 | 0.08669339802 | 0.06905512698 | 0.2333333333 | 0.4 | -0.1666666667 |
| 32 | 3 | 0.08 | 30 | 0.05615566826 | 0.06403002628 | 0.08591411069 | 0.3666666667 | 0.5 | -0.1333333333 |
| 32 | 3 | 0.1 | 30 | 0.03424035415 | 0.04931442872 | 0.05077661879 | 0.4 | 0.3 | 0.1 |
| 32 | 3 | 0.15 | 30 | 0.02617531489 | 0.03232206173 | 0.03715042959 | 0.3666666667 | 0.3 | 0.06666666667 |
| 32 | 3 | 0.2 | 30 | 0.02788850035 | 0.04236342237 | 0.04609427154 | 0.5333333333 | 0.3 | 0.2333333333 |
| 32 | 4 | 0.01 | 30 | 0.3922232458 | 0.1852650701 | 0.2907930647 | 0 | 0 | 0 |
| 32 | 4 | 0.02 | 30 | 0.2725559821 | 0.1327633618 | 0.1663741695 | 0 | 0 | 0 |
| 32 | 4 | 0.04 | 30 | 0.1938173592 | 0.1212179252 | 0.2040368457 | 0 | 0 | 0 |
| 32 | 4 | 0.06 | 30 | 0.1134046765 | 0.09774461938 | 0.1138474312 | 0.06666666667 | 0.125 | -0.05833333333 |
| 32 | 4 | 0.08 | 30 | 0.09467885857 | 0.0948705417 | 0.1552930796 | 0.2666666667 | 0.375 | -0.1083333333 |
| 32 | 4 | 0.1 | 30 | 0.06088045465 | 0.05471257332 | 0.08105573065 | 0.1333333333 | 0 | 0.1333333333 |
| 32 | 4 | 0.15 | 30 | 0.09367870815 | 0.09849976039 | 0.08962370154 | 0.1 | 0.125 | -0.025 |
| 32 | 4 | 0.2 | 30 | 0.1305842256 | 0.1491811118 | 0.1630001055 | 0.2 | 0.25 | -0.05 |
| 48 | 2 | 0.01 | 30 | 0.1551851538 | 0.1031007184 | 0.1904161585 | 0.03333333333 | NA | NA |
| 48 | 2 | 0.02 | 30 | 0.0863638401 | 0.1080636272 | 0.1018233746 | 0.3 | NA | NA |
| 48 | 2 | 0.04 | 30 | 0.005250362676 | 0.01340770383 | 0 | 0.7666666667 | NA | NA |
| 48 | 2 | 0.06 | 30 | 0.001814220606 | 0.00976987696 | 0 | 0.9666666667 | NA | NA |
| 48 | 2 | 0.08 | 30 | 0.0001985595726 | 0.001069276022 | 0 | 0.9666666667 | NA | NA |
| 48 | 2 | 0.1 | 30 | 0 | 0 | 0 | 1 | NA | NA |
| 48 | 2 | 0.15 | 30 | 0.0008958464809 | 0.004028322108 | 0 | 0.9333333333 | NA | NA |
| 48 | 2 | 0.2 | 30 | 0.00190608472 | 0.0074739288 | 0 | 0.9333333333 | NA | NA |
| 48 | 3 | 0.01 | 30 | 0.1727266432 | 0.08439302331 | 0.09515986124 | 0 | NA | NA |
| 48 | 3 | 0.02 | 30 | 0.08655145602 | 0.07382927253 | 0.1217738847 | 0.1 | NA | NA |
| 48 | 3 | 0.04 | 30 | 0.0165283147 | 0.02380242988 | 0.02201971876 | 0.4666666667 | NA | NA |
| 48 | 3 | 0.06 | 30 | 0.01283765733 | 0.02179025259 | 0.01441163677 | 0.5666666667 | NA | NA |
| 48 | 3 | 0.08 | 30 | 0.006488544226 | 0.0210126959 | 0.0008058647895 | 0.7 | NA | NA |
| 48 | 3 | 0.1 | 30 | 0.003155551601 | 0.01057278261 | 0 | 0.8 | NA | NA |
| 48 | 3 | 0.15 | 30 | 0.001963226772 | 0.008042472758 | 0 | 0.9333333333 | NA | NA |
| 48 | 3 | 0.2 | 30 | 0.009471471183 | 0.01812230356 | 0.008921920221 | 0.7 | NA | NA |
| 48 | 4 | 0.01 | 30 | 0.2570969107 | 0.153647105 | 0.2459160688 | 0.03333333333 | NA | NA |
| 48 | 4 | 0.02 | 30 | 0.1363251487 | 0.1125692008 | 0.1689135303 | 0 | NA | NA |
| 48 | 4 | 0.04 | 30 | 0.06603536613 | 0.05914062947 | 0.0756569838 | 0.1333333333 | NA | NA |
| 48 | 4 | 0.06 | 30 | 0.01890510048 | 0.02550354144 | 0.02424073798 | 0.3666666667 | NA | NA |
| 48 | 4 | 0.08 | 30 | 0.02821638747 | 0.039294958 | 0.03198619801 | 0.3333333333 | NA | NA |
| 48 | 4 | 0.1 | 30 | 0.02201134772 | 0.0457942415 | 0.008482238064 | 0.6 | NA | NA |
| 48 | 4 | 0.15 | 30 | 0.03244079685 | 0.05331857904 | 0.03363866076 | 0.5 | NA | NA |
| 48 | 4 | 0.2 | 30 | 0.0246225851 | 0.03801636987 | 0.01915065663 | 0.3333333333 | NA | NA |
| 64 | 2 | 0.01 | 30 | 0.07711519203 | 0.07260000994 | 0.0819046911 | 0.1 | 0.2 | -0.1 |
| 64 | 2 | 0.02 | 30 | 0.02506546659 | 0.03544670502 | 0.02768120706 | 0.4333333333 | 0.6 | -0.1666666667 |
| 64 | 2 | 0.04 | 30 | 0.0005500040386 | 0.002216683009 | 0 | 0.9333333333 | 1 | -0.06666666667 |
| 64 | 2 | 0.06 | 30 | 0 | 0 | 0 | 1 | 1 | 0 |
| 64 | 2 | 0.08 | 30 | 0 | 0 | 0 | 1 | 1 | 0 |
| 64 | 2 | 0.1 | 30 | 0 | 0 | 0 | 1 | 1 | 0 |
| 64 | 2 | 0.15 | 30 | 0 | 0 | 0 | 1 | 1 | 0 |
| 64 | 2 | 0.2 | 30 | 0 | 0 | 0 | 1 | 1 | 0 |
| 64 | 3 | 0.01 | 30 | 0.1314978409 | 0.09926449803 | 0.110850865 | 0.03333333333 | 0 | 0.03333333333 |
| 64 | 3 | 0.02 | 30 | 0.03356873749 | 0.03765933186 | 0.04626149763 | 0.2666666667 | 0.2 | 0.06666666667 |
| 64 | 3 | 0.04 | 30 | 0.004625532466 | 0.01060975648 | 0.001224945401 | 0.6666666667 | 0.7 | -0.03333333333 |
| 64 | 3 | 0.06 | 30 | 0.001045315585 | 0.004772801763 | 0 | 0.9333333333 | 0.9 | 0.03333333333 |
| 64 | 3 | 0.08 | 30 | 0 | 0 | 0 | 1 | 1 | 0 |
| 64 | 3 | 0.1 | 30 | 0.001813759223 | 0.007657772035 | 0 | 0.9333333333 | 0.9 | 0.03333333333 |
| 64 | 3 | 0.15 | 30 | 0.001847393225 | 0.007202989368 | 0 | 0.9333333333 | 0.9 | 0.03333333333 |
| 64 | 3 | 0.2 | 30 | 0.007781938816 | 0.01948631641 | 0 | 0.7666666667 | 0.6 | 0.1666666667 |
| 64 | 4 | 0.01 | 30 | 0.1502328128 | 0.1078179152 | 0.1286368151 | 0.03333333333 | 0 | 0.03333333333 |
| 64 | 4 | 0.02 | 30 | 0.08689383263 | 0.06296867855 | 0.09550653406 | 0.1 | 0.2222222222 | -0.1222222222 |
| 64 | 4 | 0.04 | 30 | 0.01856650321 | 0.01982679401 | 0.02696483574 | 0.3 | 0.2222222222 | 0.07777777778 |
| 64 | 4 | 0.06 | 30 | 0.004938483724 | 0.009283668714 | 0.004924207363 | 0.6666666667 | 0.6666666667 | 0 |
| 64 | 4 | 0.08 | 30 | 0.005912468089 | 0.01358810725 | 0.0007740781457 | 0.7 | 0.6666666667 | 0.03333333333 |
| 64 | 4 | 0.1 | 30 | 0.02621345673 | 0.06383139427 | 0.009273436137 | 0.6 | 0.6666666667 | -0.06666666667 |
| 64 | 4 | 0.15 | 30 | 0.02076948288 | 0.04506319529 | 0.004390114818 | 0.7 | 0.8888888889 | -0.1888888889 |
| 64 | 4 | 0.2 | 30 | 0.01362238465 | 0.02459059458 | 0.009208961186 | 0.6333333333 | 0.7777777778 | -0.1444444444 |

## Interpretation

The Phase 4B/5 outcomes are sensitive to the choice of `OPTIMIZER_SEED`.  At least one of (a) curve-shape labels flip across optimizer seeds, (b) loss IQR is comparable to the mean loss, or (c) the floor-saturation rate drops materially below the Phase 4A baseline.  The Phase 3F `INTERMEDIATE` signal and the Phase 5 `INSUFFICIENT` censoring should be re-read as potentially contaminated by single-seed variance.  This justifies a follow-up moderate intervention (e.g., reheating, larger K, or a multi-start aggregate target) under pre-registered criteria; it does not by itself establish a physical claim.

## Scope

- `loss` retains the Phase 4A/4B definition `|warmup_delta_energy / initial_energy|`.  It is an optimizer/embedding-response diagnostic, not a physical observable.
- No Phase 4A, 4B, or 5 CSV, label, threshold, or outcome is modified.  Phase 3F is not rerun.  PySR is not invoked.
- K = 1 with `optimizer_seed = 1987` reproduces Phase 4A bit-for-bit (see `--verify-against-phase4a`).

## Reproducibility

Regenerate via `make regen-phase4c`.
Source: `tools/build_phase4c_optimizer_seed_probe.py`.
