# N=24 Gamma Probe

This folder contains an exploratory gamma scan for SORKIN-2/Bombelli
annealer accessibility. It is a diagnostic of the historical annealer's
behavior on one known-truth case, not a confirmatory result.

Protocol:

- `N = 24`
- `T0 = 100`
- gamma grid: logarithmic grid from `0.5` to `0.9`
- budget: short Phase 2B schedule
- `warmup_limit = 10`
- `anneal_limit = 10`
- `max_data = 4`

The gamma values are fixed grid nodes. In particular,
`gamma = 0.562373` is the predefined logarithmic-grid value

```text
gamma = 0.5 * (0.9 / 0.5)^(2/10)
```

It must not be interpreted as special because it is close to
`10^(-1/4)`.

Current diagnostic readout:

- `truth_energy = 0` for every row, so `final_energy` and `energy_gap`
  are equivalent in this run.
- Every `success_flag` value is `false`.
- Ranking by `final_energy` selects `gamma = 0.562373`.
- Ranking by `interval_rmse` selects `gamma = 0.848624`.
- The rankings disagree.
- The apparent gamma optimum is not robust; it is metric-dependent.

Conservative conclusion: this run does not establish a robust gamma
optimum. It is an accessibility diagnostic only. It makes no physical
resonance claim, no embeddability claim, and no recovery claim.

Files:

- `run_gamma_n24_probe.py`: permanent runner used to generate the N=24
  exploratory scan.
- `gamma_n24_probe.csv`: original numerical result table.
- `gamma_n24_probe.md`: original run summary.
- `gamma_n24_probe.svg`: original final-energy plot.
- `build_gamma_n24_ranking_diagnostic.py`: post-analysis script that
  reads only `gamma_n24_probe.csv`.
- `gamma_n24_ranking_diagnostic.csv`: derived ranking table.
- `gamma_n24_ranking_diagnostic.md`: ranking diagnostic note.
- `gamma_n24_energy_vs_rmse.svg`: final-energy versus RMSE comparison
  figure.
