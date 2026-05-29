# S4-KERR-K18C residual boundary-sensitivity diagnostic 001

## Status

Completed (single-candidate boundary-sensitivity diagnostic). Overall classification: **INTERIOR** (b-cut: INTERIOR, lambda-cut: INTERIOR).

## Fixed candidate

- N = 24
- seed = 1961
- spin_a = 0.5
- event_A = 15
- event_B = 4
- direction = ingoing (probe direction held fixed at K18B best_direction)
- K17d best_b = 1.0 (upper edge of old PROBE_B_GRID)
- K17d best_lambda = 0.5 (lower edge of old PROBE_LAMBDA_GRID)
- sector m: optimized internally by k17._eval_trial exactly as in K17d (recorded as best_sector_m).

## Input artifacts / code actually used

- `explore/sorkin4_kerr_benchmark/audit_kerr_k17_controlled_candidate_pair_sandbox_001.py` (`k17._eval_trial` — the residual evaluation, reused unchanged)
- `explore/sorkin4_kerr_benchmark/audit_kerr_k17d_cloud_size_seed_scan_001.py` (`_probe_best` call pattern and PROBE grids reproduced)
- `explore/sorkin4_kerr_benchmark/run_kerr_minimal_benchmark.py` (`kerr_horizon_radius`, `generate_exterior_events`, `Event`)
- `explore/sorkin4_schwarzschild_benchmark/run_schwarzschild_minimal_benchmark.py` (`EXTERIOR_MARGIN`)
- `explore/sorkin4_kerr_benchmark/kerr_k18b_single_candidate_geometric_audit_001.md` (candidate provenance / reference values)

cones.py and production geometry were NOT modified or imported beyond the existing K17 audit stack.

## Exact commands used

```bash
python3 -m explore.sorkin4_kerr_benchmark.audit_kerr_k18c_residual_boundary_sensitivity_001
```

## Provenance / reproduction guard

- regenerated cloud: N=24, seed=1961, spin_a=0.5, n_events=24, r_plus=1.86602540378
- event_A coords (t,r,phi) = (2.91035674173, 4.91768037017, 4.95722245691) — matches K18B: True
- event_B coords (t,r,phi) = (3.89832985076, 4.55058878914, 5.10652257915) — matches K18B: True
- recomputed residual at (b=1.0, lambda=0.5, ingoing) = 0.125469416355 (K18B reference 0.125469416355); reproduced = True

## Probe grid used for b and lambda

- B_CUT_GRID = (-0.5, 0.0, 0.5, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0) (lambda fixed at 0.5)
- LAMBDA_CUT_GRID = (0.03125, 0.0625, 0.125, 0.25, 0.5, 1.0, 2.0) (b fixed at 1.0)
- CORNER_B_GRID = (0.5, 1.0, 1.5, 2.0)
- CORNER_LAMBDA_GRID = (0.125, 0.25, 0.5, 1.0)
- (reference) K17d PROBE_B_GRID = (-1.0, -0.5, 0.0, 0.5, 1.0)
- (reference) K17d PROBE_LAMBDA_GRID = (0.5, 1.0, 2.0)

## 1D b-cut table

Fixed: lambda = 0.5, direction = ingoing. Extension direction: increasing b beyond old upper edge 1.0.

| b | weighted_residual | best_sector_m | t_residual | r_residual | phi_residual_sector_adj | note |
|---:|---:|---:|---:|---:|---:|---|
| -0.5 | 0.152064657187 | 0 | -0.110625405502 | -0.137836879073 | -0.152064657187 | in K17d grid |
| 0 | 0.140777070426 | 0 | -0.114942394894 | -0.137008915454 | -0.140777070426 | in K17d grid |
| 0.5 | 0.132908418971 | 0 | -0.119485382039 | -0.132908418971 | -0.129512926228 | in K17d grid |
| 1 | 0.125469416355 | 0 | -0.124239835393 | -0.125469416355 | -0.11829642045 | old best_b (upper edge); in K17d grid |
| 1.25 | 0.126692588759 | 0 | -0.126692588759 | -0.120460405047 | -0.112713845497 | extension |
| 1.5 | 0.129194058286 | 0 | -0.129194058286 | -0.11456280398 | -0.107152740234 | extension |
| 2 | 0.134339543208 | 0 | -0.134339543208 | -0.0999828745228 | -0.0961088055523 | extension |
| 2.5 | 0.139671755012 | 0 | -0.139671755012 | -0.0814232878119 | -0.0851943030434 | extension |
| 3 | 0.145191612238 | 0 | -0.145191612238 | -0.0584351183212 | -0.0744432660713 | extension |
| 4 | 0.156844649011 | 0 | -0.156844649011 | 0.00386133469726 | -0.053607511498 | extension |

b-cut classification: **INTERIOR**

## 1D lambda-cut table

Fixed: b = 1.0, direction = ingoing. Extension direction: decreasing lambda below old lower edge 0.5.

| lambda | weighted_residual | best_sector_m | t_residual | r_residual | phi_residual_sector_adj | note |
|---:|---:|---:|---:|---:|---:|---|
| 0.03125 | 0.935780055206 | 0 | -0.935780055206 | 0.336265777143 | -0.147581908817 | extension |
| 0.0625 | 0.883371269397 | 0 | -0.883371269397 | 0.305444965046 | -0.145837612179 | extension |
| 0.125 | 0.777888649986 | 0 | -0.777888649986 | 0.243818667618 | -0.142268303763 | extension |
| 0.25 | 0.564130267405 | 0 | -0.564130267405 | 0.120629955801 | -0.134788415312 | extension |
| 0.5 | 0.125469416355 | 0 | -0.124239835393 | -0.125469416355 | -0.11829642045 | old best_lambda (lower edge); in K17d grid |
| 1 | 0.817590889728 | 0 | 0.817590889728 | -0.616338892994 | -0.077421219285 | in K17d grid |
| 2 | 3.14968093396 | 0 | 3.14968093396 | -1.59012979095 | 0.0642904939157 | in K17d grid |

lambda-cut classification: **INTERIOR**

## 2D corner table

Direction = ingoing. Small 2D grid around the corner (b=1.0, lambda=0.5). Cell = weighted_residual.

| b \ lambda | 0.125 | 0.25 | 0.5 | 1 |
|---:|---:|---:|---:|---:|
| 0.5 | 0.776951190596 | 0.562106450398 | 0.132908418971 | 0.831244273734 |
| 1 | 0.777888649986 | 0.564130267405 | 0.125469416355 | 0.817590889728 |
| 1.5 | 0.778836510477 | 0.566198231774 | 0.129194058286 | 0.802894538159 |
| 2 | 0.779794587693 | 0.568309230433 | 0.134339543208 | 0.787250131769 |

## K19 exploratory null: pair-label shuffle

Minimal K19 exploratory add-on: same-cloud pair-label shuffle null.

- seed = 19001
- requested_n = 64
- evaluated_n = 64
- fixed b = 1.0
- fixed lambda = 0.5
- fixed direction = ingoing
- baseline weighted_residual = 0.125469416355
- finite null n = 64
- null min / median / max = 0.209695069865 / 2.55648242291 / 5.53503314056
- null samples with residual <= baseline = 0
- rank_fraction_lower_is_better = 0
- baseline_below_all_nulls = True

| A_index | B_index | weighted_residual | best_sector_m | t_residual | r_residual | phi_residual_sector_adj |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 13 | 1.92870570851 | 1 | 0.367490753893 | -1.86048958544 | 1.92870570851 |
| 2 | 16 | 2.94469499158 | 0 | 0.181132083155 | -0.290200002029 | -2.94469499158 |
| 13 | 23 | 2.14684957072 | -1 | -0.567083291331 | -0.258458952911 | -2.14684957072 |
| 20 | 7 | 2.43396258166 | 0 | 2.43396258166 | 0.029546866458 | 0.824510910778 |
| 8 | 18 | 4.27338642128 | 1 | 4.27338642128 | 0.322685764153 | 3.03186815169 |
| 6 | 18 | 3.28349501495 | 0 | 3.28349501495 | -0.00572305578326 | 0.887047430464 |
| 11 | 14 | 3.54146321854 | 0 | 2.68368717509 | -3.54146321854 | 2.17458540719 |
| 6 | 3 | 1.9328910562 | 0 | 1.65816135379 | -0.386966786563 | 1.9328910562 |
| 10 | 2 | 1.77057118562 | -1 | 1.37396221576 | -0.242978440129 | -1.77057118562 |
| 12 | 16 | 3.28708007464 | 0 | 3.28708007464 | 0.56478762173 | -2.67478562254 |
| 3 | 12 | 3.14108689474 | -1 | -0.895829416889 | 0.755452303507 | -3.14108689474 |
| 7 | 9 | 2.72900031509 | 1 | 1.002904314 | -0.331443352001 | 2.72900031509 |
| 19 | 20 | 3.08919882075 | -1 | 0.348592919924 | -2.72269629 | -3.08919882075 |
| 8 | 16 | 3.2644610552 | 0 | 3.2644610552 | 2.24591919163 | -1.80665785471 |
| 7 | 12 | 2.50709992044 | 0 | -2.50709992044 | 0.38340350768 | -0.208230096649 |
| 3 | 4 | 1.58173701619 | 0 | -1.01661791929 | 0.468547954903 | -1.58173701619 |
| 3 | 17 | 1.51666547147 | 0 | 1.50042157642 | 1.24780265242 | 1.51666547147 |
| 5 | 18 | 3.78888459768 | 1 | 3.78888459768 | 0.0499100846994 | 2.96397139004 |
| 10 | 5 | 2.90757560397 | -1 | -1.73483603962 | -2.51580413574 | -2.90757560397 |
| 22 | 19 | 2.46123285852 | 0 | -0.306123839137 | 2.02103619106 | -2.46123285852 |
| 23 | 12 | 0.209695069865 | 0 | -0.128474913747 | -0.209695069865 | 0.0648738567859 |
| 14 | 18 | 3.65166355642 | 0 | 3.65166355642 | 0.233676947821 | -0.714187848136 |
| 1 | 16 | 2.9817453412 | 0 | 2.9817453412 | 0.832401602847 | 0.763912542525 |
| 0 | 20 | 3.72173542733 | 0 | 3.72173542733 | -3.27862701095 | 0.200156652201 |
| 17 | 23 | 1.66279460513 | 0 | -0.230353323087 | -1.26156481521 | 1.66279460513 |
| 16 | 2 | 3.14151970748 | -1 | 2.55419107558 | -0.677814904077 | -3.14151970748 |
| 0 | 17 | 4.42127989943 | 0 | 4.42127989943 | -1.3857961355 | -0.858502325961 |
| 6 | 16 | 2.33170673124 | 0 | 2.27456964887 | 1.9175103717 | 2.33170673124 |
| 8 | 9 | 4.46262653071 | 0 | 4.46262653071 | 0.476095525477 | -2.4326148153 |
| 22 | 17 | 1.6762593583 | 0 | -0.539542206862 | 1.6762593583 | -0.412496250162 |
| 3 | 8 | 2.24811144785 | 0 | -1.08317087951 | -0.93036269154 | 2.24811144785 |
| 17 | 20 | 2.38040895795 | 0 | 0.418033464989 | -2.38040895795 | 1.12283988721 |
| 7 | 18 | 1.9102979749 | 1 | 0.813664204567 | -0.484853113325 | 1.9102979749 |
| 10 | 3 | 2.35275985999 | 0 | -0.357074729218 | -2.35275985999 | 1.08055612074 |
| 4 | 23 | 2.08803926717 | -1 | 2.08803926717 | -0.486082934493 | -1.54823110673 |
| 11 | 4 | 2.23171539364 | 0 | 2.12215538417 | -2.23171539364 | 0.881149895882 |
| 17 | 3 | 2.2294287832 | 0 | 0.416757407628 | -2.2294287832 | -1.42919975988 |
| 13 | 19 | 1.01426633657 | 0 | 1.01426633657 | 0.860304612561 | 0.488985432009 |
| 10 | 13 | 1.61857228515 | 0 | 0.445692236594 | -1.61857228515 | 0.0695871886956 |
| 16 | 4 | 1.89541421985 | 0 | -0.993064835651 | -1.82500032075 | -1.89541421985 |
| 1 | 15 | 1.5370911889 | 0 | 1.5370911889 | -0.876571133258 | -1.090624829 |
| 1 | 8 | 2.58992351279 | 0 | 0.482565119647 | -1.90839019867 | 2.58992351279 |
| 23 | 21 | 1.95787751325 | 0 | 1.95787751325 | -1.70205405343 | -0.899497419489 |
| 9 | 11 | 2.86665350438 | 0 | -2.86665350438 | 1.688715317 | -2.03870789749 |
| 5 | 20 | 2.16482699388 | 0 | 2.16482699388 | -0.482313820831 | 0.278669355688 |
| 20 | 22 | 2.79162739163 | 0 | 2.79162739163 | -0.772643630249 | -0.604958522532 |
| 2 | 7 | 2.22175600798 | 0 | 1.21319442739 | -2.22175600798 | 0.011243971438 |
| 18 | 0 | 2.68737209392 | -1 | -2.62244239107 | 1.7741051848 | -2.68737209392 |
| 10 | 23 | 2.11063482142 | -1 | -1.00418545993 | -1.384895892 | -2.11063482142 |
| 16 | 23 | 2.8016019639 | 0 | 0.176043399892 | -1.81973235595 | 2.8016019639 |
| 4 | 13 | 3.5379169637 | 0 | 3.5379169637 | -0.71975932764 | 0.63199090338 |
| 3 | 19 | 1.73383994414 | 0 | 1.73383994414 | 1.59257948518 | -0.532071136888 |
| 0 | 9 | 5.53503314056 | 0 | 5.53503314056 | -2.5929933441 | -2.57902428043 |
| 19 | 3 | 2.57171611525 | 0 | 0.347316862563 | -2.57171611525 | 0.641946839345 |
| 3 | 20 | 2.57532444963 | 0 | 0.800877104325 | -0.645028223029 | 2.57532444963 |
| 14 | 8 | 2.55648242291 | 0 | 0.143557968789 | -0.58388142622 | 2.55648242291 |
| 20 | 23 | 0.624489682139 | 0 | 0.138381211632 | 0.624489682139 | 0.561750637738 |
| 23 | 18 | 3.19228921126 | 1 | 3.19228921126 | -1.07795169087 | 2.18340192833 |
| 19 | 0 | 3.11909301495 | -1 | -1.48752839923 | 0.0767868938532 | -3.11909301495 |
| 3 | 2 | 2.82784250382 | -1 | 2.53063799194 | 1.61573337158 | -2.82784250382 |
| 5 | 22 | 4.1696863858 | 0 | 4.1696863858 | -0.760602990479 | -0.348085086654 |
| 7 | 14 | 2.64455529321 | 1 | -2.06635663193 | -1.21324866583 | 2.64455529321 |
| 12 | 15 | 1.84242592234 | 1 | 1.84242592234 | -1.14418511438 | 1.75386231311 |
| 17 | 5 | 2.39247305895 | 0 | -0.961003902774 | -2.39247305895 | 0.865853822584 |

Reading: this null only asks whether the reproduced baseline endpoint pair is exceptional under the same residual objective when A/B labels are shuffled inside the same generated cloud. It does not preserve boundary/sampling structure beyond using the same cloud, and it is not a Kerr confirmation test.

K19-001 verdict: if baseline_below_all_nulls is true and zero null samples match or improve the baseline, the reproduced K18B endpoint pair is exceptional against this same-cloud pair-label shuffle null. This supports continuing to a stronger boundary/sampling-preserving null, but does not establish causal reachability or Kerr structure.

## Classification: INTERIOR / RUNAWAY / INCONCLUSIVE

- b-cut: **INTERIOR**
- lambda-cut: **INTERIOR**
- overall: **INTERIOR**

Pre-registered rule:
- RUNAWAY if the minimum residual on a cut sits at the extreme grid point in the extension direction (residual keeps improving off the new edge).
- INTERIOR if the minimum is strictly interior / on the non-extension side (extending past the old K17d edge does not keep lowering the residual).
- INCONCLUSIVE if a cut has fewer than three finite residuals, or if the provenance/reproduction guard fails.
- overall = INCONCLUSIVE if any cut is INCONCLUSIVE; else RUNAWAY if any cut is RUNAWAY; else INTERIOR.

## Interpretation

The residual attains an interior minimum within the extended b/lambda range (it stops improving once the old K17d edge is passed). Under the pre-registered reading this indicates the old K17d probe grid was too narrow, so K17d/K18A residual rankings are NOT final and would need re-evaluation on a corrected grid. This is an objective-grid statement only; no causal/reachability claim follows (see guardrails).

## Guardrails

- residual profile != causal_true
- interior minimum != reachability
- runaway residual != spacelike separation
- K18C is a numerical objective diagnostic only
- no global causal claim is made
- boundary-sensitivity of a residual is a statement about the objective's parametrization, not about Kerr causal structure
- this audit is local to ONE K18A/K18B candidate family

## Next operational recommendation

Single next step: treat K17d/K18A residual rankings as non-final and re-evaluate the top candidates on a corrected (interior-containing) b/lambda grid before any further interpretation. Still no causal/reachability claim.

