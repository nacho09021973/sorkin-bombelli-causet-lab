# Bombelli Thesis vs Revived Port

This note records the current comparison between the 1987 thesis program and the revived Python port in `cones.py`.

## What Is Faithfully Ported

- The input format is the original Pascal upper-triangular incidence matrix.
- The random generator structure follows the thesis `ran2` and `gasdev` routines.
- The annealing loop keeps the thesis shape:
  - warmup block
  - anneal blocks
  - acceptance test proportional to `4 * exp(-deltaE / T)`
  - multiplicative cooling

## Thesis Defaults

The historical defaults preserved in `cones.py` are:

- initial temperature: `100.0`
- cooling factor: `0.9`
- acceptance scale: `4.0`
- warmup limit: `100`
- anneal limit: `100`

## Benchmarks Used

- `benchmarks/tesis_like_6.in`
- `benchmarks/tesis_like_12.in`

Both were generated reproducibly from sprinklings with seed `1987`.

## Control Run On `tesis_like_12.in`

Using the thesis defaults on the 12-event benchmark gave:

- mean final energy: `20.021376`
- best tested seed: `1962`
- zero rate: `0.00`

This is a faithful reproduction of the algorithmic behavior, but not a particularly good optimization outcome for this benchmark.

## Empirical Schedule Search

A finer schedule sweep around the best coarse candidate found:

- initial temperature: `180.0`
- cooling factor: `0.8`
- mean final energy: `0.166158`
- zero rate: `0.00`

That schedule is far better on this benchmark than the thesis defaults.
The mean final energy is lower by about `99.17%` relative to the thesis-default control run on the same benchmark.

## Interpretation

- The thesis schedule is historically correct.
- The port is behaving as expected.
- The optimization landscape is highly schedule-sensitive.
- For this benchmark, the historical defaults are not near-optimal.

## Practical Conclusion

The revived code now serves two roles:

1. A faithful historical port of Bombelli's 1987 annealing experiment.
2. A reproducible laboratory for exploring better schedules, dimensions, and benchmark causets.

The forward research program is recorded in [`research_agenda_2026.md`](/home/adnac/sorkin/research_agenda_2026.md).
The first empirical findings are recorded in [`results_note_2026.md`](/home/adnac/sorkin/results_note_2026.md).

## Relevant Files

- `cones.py`
- `schedule_sweep.py`
- `dimension_sweep.py`
- `analyze_sweep.py`
- `research_agenda_2026.md`
- `results_note_2026.md`
- `benchmarks/README.md`
