.RECIPEPREFIX := >
.PHONY: test regen-fixtures regen-foundation regen-phase1 regen-phase1b regen-phase1c regen-phase1d smoke smoke-cuda ensemble phase refine

test:
>python3 -m unittest discover -s tests -v

regen-fixtures:
>python3 tools/update_test_fixtures.py
>python3 -m unittest discover -s tests -v

regen-foundation:
>python3 tools/build_foundation_benchmarks.py
>python3 -m unittest tests.test_foundation_benchmarks -v

regen-phase1:
>python3 tools/build_phase1_atlas.py
>python3 -m unittest tests.test_phase1_atlas -v

regen-phase1b:
>python3 tools/build_phase1b_scaling_atlas.py
>python3 -m unittest tests.test_phase1b_scaling_atlas -v

regen-phase1c:
>python3 tools/build_phase1c_scaling_atlas.py
>python3 -m unittest tests.test_phase1c_scaling_atlas -v

regen-phase1d:
>python3 tools/build_phase1d_structural_atlas.py
>python3 -m unittest tests.test_phase1d_structural_atlas -v

smoke:
>python3 cones.py benchmarks/tesis_like_6.in --dim 2 --output /tmp/cone.out --plot /tmp/cone.svg --max-data 5

smoke-cuda:
>python3 cones.py benchmarks/tesis_like_6.in --dim 2 --backend cuda --output /tmp/cone_cuda.out --plot /tmp/cone_cuda.svg --max-data 5

ensemble:
>python3 ensemble_scan.py benchmarks/tesis_like_6.in --dim 2 --gpu-first --backend auto --seed-start 1959 --seed-count 2 --seed-step 1 --temp-min 100 --temp-max 180 --temp-count 2 --cooling-min 0.8 --cooling-max 0.9 --cooling-count 2 --workers 1 --run-csv /tmp/ensemble_runs.csv --summary-csv /tmp/ensemble_summary.csv --report-md /tmp/ensemble.md --heatmap-svg /tmp/ensemble.svg

phase:
>python3 phase_diagram.py --n-values 6,12 --dim-min 1 --dim-max 4 --seed-start 1959 --seed-count 4 --gpu-first --backend auto --fast-frontier --runs-csv /tmp/phase_runs.csv --summary-csv /tmp/phase_summary.csv --report-md /tmp/phase.md --heatmap-svg /tmp/phase.svg

refine:
>python3 phase_refine.py --cells 6:3,6:4,12:2,16:2 --seed-start 1959 --seed-count 16 --gpu-first --backend auto --partial cell --runs-csv /tmp/phase_refine_runs.csv --summary-csv /tmp/phase_refine_summary.csv --report-md /tmp/phase_refine.md --heatmap-svg /tmp/phase_refine.svg
