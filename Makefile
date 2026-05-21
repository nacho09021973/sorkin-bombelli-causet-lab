.RECIPEPREFIX := >
.PHONY: test regen-fixtures regen-foundation regen-phase1 regen-phase1b regen-phase1c regen-phase1d regen-phase1e regen-phase2 regen-phase2b regen-phase2c regen-phase2d regen-phase2e regen-phase2f regen-phase2g regen-phase3a regen-phase3b regen-phase3c regen-phase3d regen-phase3e regen-phase3f regen-phase4a regen-phase4b test-phase4b regen-phase4c test-phase4c regen-phase4d test-phase4d regen-phase5 test-phase5 smoke smoke-cuda ensemble phase refine

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

regen-phase2:
>python3 tools/build_phase2_embedding_bridge.py
>python3 -m unittest tests.test_phase2_embedding_bridge -v

regen-phase2b:
>python3 tools/build_phase2b_annealer_schedule_probe.py
>python3 -m unittest tests.test_phase2b_annealer_schedule_probe -v

regen-phase2c:
>python3 tools/build_phase2c_oracle_embedding_audit.py
>python3 -m unittest tests.test_phase2c_oracle_embedding_audit -v

regen-phase2d:
>python3 tools/build_phase2d_initialization_basin_audit.py
>python3 -m unittest tests.test_phase2d_initialization_basin_audit -v

regen-phase2e:
>python3 tools/build_phase2e_warmup_skip_probe.py
>python3 -m unittest tests.test_phase2e_warmup_skip_probe -v

regen-phase2f:
>python3 tools/build_phase2f_guarded_warmup_probe.py
>python3 -m unittest tests.test_phase2f_guarded_warmup_probe -v

# Phase3, Phase4A/4B, and Phase5 targets are exploratory and are not provenance-freeze benchmark artifacts.
regen-phase3a:
>python3 tools/build_phase3a_pysr_warmup_rule.py

regen-phase3b:
>python3 tools/build_phase3b_pysr_order_only.py

regen-phase3c:
>python3 tools/build_phase3c_ablation.py

regen-phase3d:
>python3 tools/build_phase3d_pysr_residual_target.py

regen-phase3e:
>python3 tools/build_phase3e_pysr_residual_by_warmup_mode.py

regen-phase1e:
>python3 tools/build_phase1e_extended_structural_atlas.py

regen-phase2g:
>python3 tools/build_phase2g_extended_guarded_warmup_probe.py

regen-phase3f:
>python3 tools/build_phase3f_pysr_final_ablation.py

regen-phase4a:
>python3 tools/build_phase4a_epsilon_sweep.py

regen-phase4b:
>python3 tools/build_phase4b_survival_probe.py --grid pilot
>python3 -m unittest tests.test_phase4b_survival_probe -v

test-phase4b:
>python3 -m unittest tests.test_phase4b_survival_probe -v

regen-phase4c:
>python3 tools/build_phase4c_optimizer_seed_probe.py
>python3 -m unittest tests.test_phase4c_optimizer_seed_probe -v

test-phase4c:
>python3 -m unittest tests.test_phase4c_optimizer_seed_probe -v

regen-phase4d:
>python3 tools/build_phase4d_robustness_audit.py
>python3 -m unittest tests.test_phase4d_robustness_audit -v

test-phase4d:
>python3 -m unittest tests.test_phase4d_robustness_audit -v

regen-phase5:
>python3 tools/build_phase5_seed_curve_morphology.py
>python3 -m unittest tests.test_phase5_seed_curve_morphology -v

test-phase5:
>python3 -m unittest tests.test_phase5_seed_curve_morphology -v

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
