# N=20 known-truth validation

This is a manufactured known-truth SORKIN-2 validation case. It does not use a community golden table.

## Configuration

- Command: `python3 explore/n20_known_truth_validation/run_n20_known_truth_validation.py`
- Generated at UTC: `2026-05-25T10:47:42+00:00`
- Runtime seconds: `7.820`
- Known-truth directory: `explore/known_truth_n20`
- family: `minkowski`
- N: `20`
- d_spacetime: `2`
- case_seed: `1959`
- optimizer_seeds: `1959, 1962, 1987, 2001`
- schedule_label: `gamma_0p5`
- cooling_factor: `0.5`
- T0: `100.0`
- budget_label: `medium_25_25_8`
- warmup_limit: `25`
- anneal_limit: `25`
- max_data: `8`
- total_relations_target: `111`
- Instrumentation: read-only `block_callback` over `sim.rold`/`sim.xold`.

## Per-seed checkpoint analysis

| optimizer_seed | final block | final F1 | final recall | final missing | final extra | final exact | best block | best temp | best F1 | best recall | best missing | best extra | any exact | delta best-final | min energy block | best=min energy |
| ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | --- |
| 1959 | 8 | 0.714286 | 0.675676 | 36 | 24 | false | 8 | 0.78125 | 0.714286 | 0.675676 | 36 | 24 | false | 0 | 8 | true |
| 1962 | 8 | 0.742268 | 0.648649 | 39 | 11 | false | 8 | 0.78125 | 0.742268 | 0.648649 | 39 | 11 | false | 0 | 8 | true |
| 1987 | 8 | 0.72 | 0.648649 | 39 | 17 | false | 8 | 0.78125 | 0.72 | 0.648649 | 39 | 17 | false | 0 | 8 | true |
| 2001 | 8 | 0.666667 | 0.585586 | 46 | 19 | false | 7 | 1.5625 | 0.71134 | 0.621622 | 42 | 14 | false | 0.0446735 | 8 | false |

## Global summary

| metric | value |
| --- | ---: |
| n_runs | 4 |
| avg_final_causal_f1 | 0.710805 |
| avg_best_checkpoint_causal_f1 | 0.721973 |
| count_exact_match_final | 0 |
| count_exact_match_any_checkpoint | 0 |
| count_best_gt_final | 1 |
| count_best_matches_min_energy | 3 |
| avg_delta_best_minus_final | 0.0111684 |
| max_delta_best_minus_final | 0.0446735 |

## Readout

1. N=20 recoverability: This budget does not show high-fidelity recovery for this N=20 case.
2. Exact match: endpoint exact matches occur in `0/4` runs; any-checkpoint exact matches occur in `0/4` runs.
3. Endpoint selection: The endpoint does not always preserve the best causal checkpoint.
4. Energy selection: Minimum energy does not consistently select the best causal F1 checkpoint.
5. N=20 vs N=36: Compared with the N=36 probes, this N=20 validation is smaller and uses a single gamma, so qualitative comparison is limited to endpoint/checkpoint behavior rather than a schedule matrix.
6. Implication for N=36: if N=20 preserves best checkpoints more often than N=36, the N=36 selection/parada hypothesis is strengthened; if not, checkpoint selection remains a broader annealer diagnostic rather than an N=36-only effect.

## Guardrails

This is one manufactured known-truth causal set, one family, one case_seed, one schedule, four optimizer seeds, and the unchanged historical Bombelli objective.
It is not a community benchmark, not an embeddability claim, not a physical gamma claim, not an N-transition claim, and not proof of general annealer failure.
