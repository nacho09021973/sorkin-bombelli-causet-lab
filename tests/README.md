# Tests

This directory contains regression tests and frozen fixtures for the revived `cones.py` pipeline.

## Coverage

- `test_regression.py`
  - parses the benchmark inputs
  - checks fixed `cones.py` runs
  - checks the `--sprinkle` CLI path
  - checks `schedule_sweep.py`
  - checks `dimension_sweep.py`
  - checks `analyze_sweep.py`
  - checks `ensemble_scan.py` in `--gpu-first` mode
  - checks `phase_diagram.py`
  - compares generated CSV and SVG outputs against frozen fixtures

## Frozen Fixtures

- `fixtures/schedule_sweep_1959.csv`
- `fixtures/dimension_sweep_1959.csv`
- `fixtures/analyze_sweep_sample.csv`
- `fixtures/ensemble_scan_runs.csv`
- `fixtures/ensemble_scan_summary.csv`
- `fixtures/ensemble_scan.md`
- `fixtures/cones_cli_smoke.svg`
- `fixtures/cones_cli_sprinkle_edge.svg`
- `fixtures/cones_cli_sprinkle_causet.svg`
- `fixtures/schedule_sweep.svg`
- `fixtures/dimension_sweep.svg`
- `fixtures/analyze_sweep.svg`
- `fixtures/ensemble_scan.svg`

## Regeneration

If the rendering or output format changes intentionally, regenerate the fixtures from the current code and update the matching tests together.

Recommended workflow:

1. Run the relevant command in a temporary directory.
2. Copy the generated CSV or SVG into `tests/fixtures/`.
3. Update the corresponding assertion in `tests/test_regression.py`.
4. Run `python3 -m unittest discover -s tests -v`.

## Notes

- The fixtures are intentionally small and deterministic.
- Paths inside CSVs are normalized in the tests so temporary output directories do not matter.
- The ensemble markdown fixture also normalizes the temporary directory prefix.
- The SVG comparisons are exact to catch accidental visual regressions.
