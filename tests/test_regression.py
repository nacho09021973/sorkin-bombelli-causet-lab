from __future__ import annotations

import contextlib
import csv
import io
import re
import tempfile
import unittest
from pathlib import Path

import analyze_sweep
import cones
import ensemble_scan
import dimension_sweep
import phase_diagram
import phase_refine
import schedule_sweep

try:
    from cuda_backend import cuda_available
except Exception:
    def cuda_available() -> bool:  # type: ignore
        return False


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = Path(__file__).resolve().parent / "fixtures"


def load_benchmark(name: str) -> list[list[bool]]:
    return cones.parse_cones_input(ROOT / "benchmarks" / name)


def load_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def normalize_generated_csv_rows(rows: list[dict[str, str]], fields: tuple[str, ...]) -> list[dict[str, str]]:
    normalized: list[dict[str, str]] = []
    for row in rows:
        item = dict(row)
        for field in fields:
            if field in item and item[field]:
                item[field] = "<tmp>"
        normalized.append(item)
    return normalized


def assert_file_matches(testcase: unittest.TestCase, generated: Path, fixture: Path) -> None:
    testcase.assertEqual(generated.read_text(encoding="utf-8"), fixture.read_text(encoding="utf-8"))


def normalize_text_paths(text: str, tmpdir: str, root: str = "") -> str:
    if root:
        text = text.replace(root, "")
    text = text.replace(tmpdir, "<tmp>")
    return re.sub(r"/tmp/ensemble_scan_[^/]+", "<tmp>", text)


def run_sim(
    benchmark: str,
    *,
    dim: int,
    seed: int,
    initial_temp: float = 100.0,
    cooling_factor: float = 0.9,
    acceptance_scale: float = 4.0,
    warmup_limit: int = 100,
    anneal_limit: int = 100,
    max_data: int = 35,
    backend: str = "cpu",
) -> cones.ConesSimulator:
    z = load_benchmark(benchmark)
    sim = cones.ConesSimulator(
        z=z,
        dim=dim,
        seed=seed,
        interactive=False,
        max_data=max_data,
        plot_path=None,
        warmup_limit=warmup_limit,
        anneal_limit=anneal_limit,
        initial_temp=initial_temp,
        cooling_factor=cooling_factor,
        acceptance_scale=acceptance_scale,
        backend=backend,
    )
    with tempfile.TemporaryDirectory() as tmpdir, contextlib.redirect_stdout(io.StringIO()):
        sim.run(Path(tmpdir) / "out.txt")
    return sim


class ConesRegressionTests(unittest.TestCase):
    def test_parse_benchmarks(self) -> None:
        z6 = load_benchmark("tesis_like_6.in")
        self.assertEqual(len(z6), 6)
        self.assertEqual(sum(1 for i in range(len(z6)) for j in range(i + 1, len(z6)) if z6[i][j]), 4)

        z12 = load_benchmark("tesis_like_12.in")
        self.assertEqual(len(z12), 12)
        self.assertEqual(sum(1 for i in range(len(z12)) for j in range(i + 1, len(z12)) if z12[i][j]), 19)

    def test_thesis_like_6_seed_1959_regression(self) -> None:
        sim = run_sim("tesis_like_6.in", dim=2, seed=1959)
        self.assertAlmostEqual(sim.initial_energy, 2.9698484809835, places=12)
        self.assertAlmostEqual(sim.warmup_energy, 2.0439711583916464, places=12)
        self.assertAlmostEqual(sim.data[-1][1], 1.645239674297646e-11, places=12)
        self.assertEqual(len(sim.data), 5)

    def test_thesis_like_12_seed_1962_regression(self) -> None:
        sim = run_sim("tesis_like_12.in", dim=2, seed=1962)
        self.assertAlmostEqual(sim.initial_energy, 14.407300666675905, places=12)
        self.assertAlmostEqual(sim.warmup_energy, 30.500528120981947, places=12)
        self.assertAlmostEqual(sim.data[-1][1], 16.80890143460476, places=12)
        self.assertEqual(len(sim.data), 35)

    def test_empirical_schedule_beats_thesis_defaults_on_benchmark_12(self) -> None:
        default_sim = run_sim("tesis_like_12.in", dim=2, seed=1962)
        tuned_sim = run_sim(
            "tesis_like_12.in",
            dim=2,
            seed=1962,
            initial_temp=180.0,
            cooling_factor=0.8,
        )

        self.assertGreater(default_sim.data[-1][1], 10.0)
        self.assertLess(tuned_sim.data[-1][1], 0.01)
        self.assertLess(tuned_sim.data[-1][1], default_sim.data[-1][1] / 1000.0)

    def test_schedule_sweep_integration(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir, contextlib.redirect_stdout(io.StringIO()):
            tmp = Path(tmpdir)
            run_csv = tmp / "runs.csv"
            summary_csv = tmp / "summary.csv"
            report_md = tmp / "report.md"
            heatmap_svg = tmp / "heatmap.svg"
            rc = schedule_sweep.main(
                [
                    str(ROOT / "benchmarks" / "tesis_like_6.in"),
                    "--dim",
                    "2",
                    "--seed-count",
                    "1",
                    "--seed-start",
                    "1959",
                    "--seed-step",
                    "1",
                    "--temp-min",
                    "100",
                    "--temp-max",
                    "100",
                    "--temp-count",
                    "1",
                    "--cooling-min",
                    "0.9",
                    "--cooling-max",
                    "0.9",
                    "--cooling-count",
                    "1",
                    "--run-csv",
                    str(run_csv),
                    "--summary-csv",
                    str(summary_csv),
                    "--report-md",
                    str(report_md),
                    "--heatmap-svg",
                    str(heatmap_svg),
                ]
            )

            self.assertEqual(rc, 0)
            self.assertTrue(run_csv.exists())
            self.assertTrue(summary_csv.exists())
            self.assertTrue(report_md.exists())
            self.assertTrue(heatmap_svg.exists())

            run_rows = normalize_generated_csv_rows(
                load_csv_rows(run_csv),
                ("input_file", "output_file", "plot_file", "causet_plot_file"),
            )
            self.assertEqual(run_rows, load_csv_rows(FIXTURES / "schedule_sweep_1959.csv"))

            summary_rows = load_csv_rows(summary_csv)
            self.assertEqual(len(summary_rows), 1)
            self.assertEqual(summary_rows[0]["mean_final"], "0.000000")
            assert_file_matches(self, heatmap_svg, FIXTURES / "schedule_sweep.svg")
            report = report_md.read_text(encoding="utf-8")
            self.assertIn("# Schedule Sweep", report)
            self.assertIn("Best schedule", report)
            self.assertIn("mean final energy", report)

    def test_cones_cli_smoke(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir, contextlib.redirect_stdout(io.StringIO()):
            tmp = Path(tmpdir)
            output_path = tmp / "cone.out"
            plot_path = tmp / "cone.svg"
            rc = cones.main(
                [
                    str(ROOT / "benchmarks" / "tesis_like_6.in"),
                    "--dim",
                    "2",
                    "--seed",
                    "1959",
                    "--output",
                    str(output_path),
                    "--plot",
                    str(plot_path),
                    "--max-data",
                    "5",
                ]
            )

            self.assertEqual(rc, 0)
            self.assertTrue(output_path.exists())
            self.assertTrue(plot_path.exists())
            self.assertIn("Total energy", output_path.read_text(encoding="utf-8"))
            assert_file_matches(self, plot_path, FIXTURES / "cones_cli_smoke.svg")

    def test_cones_cli_sprinkle_smoke(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir, contextlib.redirect_stdout(io.StringIO()):
            tmp = Path(tmpdir)
            output_path = tmp / "sprinkle.out"
            plot_path = tmp / "sprinkle.svg"
            causet_plot = tmp / "causet.svg"
            generated_input = tmp / "sprinkle.in"
            rc = cones.main(
                [
                    "--sprinkle",
                    "4",
                    "--sprinkle-spacetime-dim",
                    "1",
                    "--dim",
                    "2",
                    "--seed",
                    "1987",
                    "--output",
                    str(output_path),
                    "--plot",
                    str(plot_path),
                    "--causet-plot",
                    str(causet_plot),
                    "--save-generated-input",
                    str(generated_input),
                    "--max-data",
                    "5",
                ]
            )

            self.assertEqual(rc, 0)
            self.assertTrue(output_path.exists())
            self.assertTrue(plot_path.exists())
            self.assertTrue(causet_plot.exists())
            self.assertTrue(generated_input.exists())
            self.assertTrue(generated_input.read_text(encoding="utf-8").startswith("4\n"))
            assert_file_matches(self, causet_plot, FIXTURES / "cones_cli_sprinkle_edge.svg")

    def test_cones_cli_sprinkle_antichain_smoke(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir, contextlib.redirect_stdout(io.StringIO()):
            tmp = Path(tmpdir)
            output_path = tmp / "sprinkle_no_edges.out"
            plot_path = tmp / "sprinkle_no_edges.svg"
            causet_plot = tmp / "causet_no_edges.svg"
            rc = cones.main(
                [
                    "--sprinkle",
                    "4",
                    "--sprinkle-spacetime-dim",
                    "3",
                    "--dim",
                    "2",
                    "--seed",
                    "1986",
                    "--output",
                    str(output_path),
                    "--plot",
                    str(plot_path),
                    "--causet-plot",
                    str(causet_plot),
                    "--max-data",
                    "5",
                ]
            )

            self.assertEqual(rc, 0)
            self.assertTrue(output_path.exists())
            self.assertTrue(plot_path.exists())
            self.assertTrue(causet_plot.exists())
            assert_file_matches(self, causet_plot, FIXTURES / "cones_cli_sprinkle_causet.svg")

    @unittest.skipUnless(cuda_available(), "CUDA device unavailable")
    def test_cuda_backend_matches_cpu_on_small_run(self) -> None:
        cpu = run_sim("tesis_like_6.in", dim=2, seed=1959)
        gpu = run_sim("tesis_like_6.in", dim=2, seed=1959, backend="cuda")
        self.assertAlmostEqual(gpu.initial_energy, cpu.initial_energy, places=12)
        self.assertAlmostEqual(gpu.warmup_energy, cpu.warmup_energy, places=12)
        self.assertAlmostEqual(gpu.data[-1][1], cpu.data[-1][1], places=12)

    def test_dimension_sweep_integration(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir, contextlib.redirect_stdout(io.StringIO()):
            tmp = Path(tmpdir)
            output_csv = tmp / "dimension.csv"
            report_md = tmp / "dimension.md"
            heatmap_svg = tmp / "dimension.svg"
            rc = dimension_sweep.main(
                [
                    str(ROOT / "benchmarks" / "tesis_like_6.in"),
                    "--dim-min",
                    "1",
                    "--dim-max",
                    "2",
                    "--seed-start",
                    "1959",
                    "--seed-count",
                    "1",
                    "--seed-step",
                    "1",
                    "--output-csv",
                    str(output_csv),
                    "--report-md",
                    str(report_md),
                    "--heatmap-svg",
                    str(heatmap_svg),
                ]
            )

            self.assertEqual(rc, 0)
            self.assertTrue(output_csv.exists())
            self.assertTrue(report_md.exists())
            self.assertTrue(heatmap_svg.exists())

            rows = normalize_generated_csv_rows(load_csv_rows(output_csv), ("output_file",))
            self.assertEqual(rows, load_csv_rows(FIXTURES / "dimension_sweep_1959.csv"))
            assert_file_matches(self, heatmap_svg, FIXTURES / "dimension_sweep.svg")

            report = report_md.read_text(encoding="utf-8")
            self.assertIn("# Dimension Sweep", report)
            self.assertIn("Per-dimension summary", report)
            self.assertIn("All runs", report)

    def test_analyze_sweep_report_generation(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            csv_path = tmp / "sweep.csv"
            report_md = tmp / "analysis.md"
            plot_svg = tmp / "analysis.svg"
            csv_path.write_text((FIXTURES / "analyze_sweep_sample.csv").read_text(encoding="utf-8"), encoding="utf-8")

            rc = analyze_sweep.main([str(csv_path), "--report-md", str(report_md), "--plot-svg", str(plot_svg)])

            self.assertEqual(rc, 0)
            self.assertTrue(report_md.exists())
            self.assertTrue(plot_svg.exists())

            report = report_md.read_text(encoding="utf-8")
            self.assertIn("Mean final energy", report)
            self.assertIn("seed 1959", report)
            self.assertIn("0.000000", report)
            assert_file_matches(self, plot_svg, FIXTURES / "analyze_sweep.svg")

    def test_ensemble_scan_integration(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir, contextlib.redirect_stdout(io.StringIO()):
            tmp = Path(tmpdir)
            run_csv = tmp / "ensemble_runs.csv"
            summary_csv = tmp / "ensemble_summary.csv"
            report_md = tmp / "ensemble.md"
            heatmap_svg = tmp / "ensemble.svg"

            rc = ensemble_scan.main(
                [
                    str(ROOT / "benchmarks" / "tesis_like_6.in"),
                    "--dim",
                    "2",
                    "--gpu-first",
                    "--backend",
                    "auto",
                    "--seed-start",
                    "1959",
                    "--seed-count",
                    "2",
                    "--seed-step",
                    "1",
                    "--temp-min",
                    "100",
                    "--temp-max",
                    "180",
                    "--temp-count",
                    "2",
                    "--cooling-min",
                    "0.8",
                    "--cooling-max",
                    "0.9",
                    "--cooling-count",
                    "2",
                    "--workers",
                    "1",
                    "--run-csv",
                    str(run_csv),
                    "--summary-csv",
                    str(summary_csv),
                    "--report-md",
                    str(report_md),
                    "--heatmap-svg",
                    str(heatmap_svg),
                ]
            )

            self.assertEqual(rc, 0)
            self.assertTrue(run_csv.exists())
            self.assertTrue(summary_csv.exists())
            self.assertTrue(report_md.exists())
            self.assertTrue(heatmap_svg.exists())

            run_rows = load_csv_rows(run_csv)
            for row in run_rows:
                row["output_file"] = f"<tmp>/{Path(row['output_file']).name}"
            self.assertEqual(run_rows, load_csv_rows(FIXTURES / "ensemble_scan_runs.csv"))

            summary_rows = load_csv_rows(summary_csv)
            self.assertEqual(summary_rows, load_csv_rows(FIXTURES / "ensemble_scan_summary.csv"))

            report_text = normalize_text_paths(report_md.read_text(encoding="utf-8"), tmpdir, f"{ROOT}/")
            self.assertEqual(report_text, (FIXTURES / "ensemble_scan.md").read_text(encoding="utf-8"))
            assert_file_matches(self, heatmap_svg, FIXTURES / "ensemble_scan.svg")

            self.assertIn("Ensemble Scan", report_text)
            self.assertIn("Best schedule by success rate", report_text)
            self.assertIn("0.50", report_text)

    def test_phase_diagram_integration(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir, contextlib.redirect_stdout(io.StringIO()):
            tmp = Path(tmpdir)
            runs_csv = tmp / "phase_runs.csv"
            summary_csv = tmp / "phase_summary.csv"
            report_md = tmp / "phase.md"
            heatmap_svg = tmp / "phase.svg"

            rc = phase_diagram.main(
                [
                    "--n-values",
                    "4,6",
                    "--dim-min",
                    "1",
                    "--dim-max",
                    "2",
                    "--seed-start",
                    "1959",
                    "--seed-count",
                    "1",
                    "--warmup-limit",
                    "10",
                    "--anneal-limit",
                    "10",
                    "--max-data",
                    "3",
                    "--runs-csv",
                    str(runs_csv),
                    "--summary-csv",
                    str(summary_csv),
                    "--report-md",
                    str(report_md),
                    "--heatmap-svg",
                    str(heatmap_svg),
                ]
            )

            self.assertEqual(rc, 0)
            self.assertTrue(runs_csv.exists())
            self.assertTrue(summary_csv.exists())
            self.assertTrue(report_md.exists())
            self.assertTrue(heatmap_svg.exists())

            run_rows = load_csv_rows(runs_csv)
            summary_rows = load_csv_rows(summary_csv)
            self.assertEqual(len(run_rows), 4)
            self.assertEqual(len(summary_rows), 4)
            self.assertEqual({row["n"] for row in summary_rows}, {"4", "6"})
            self.assertEqual({row["dim"] for row in summary_rows}, {"1", "2"})
            self.assertTrue(all("success_rate" in row for row in summary_rows))

            report = report_md.read_text(encoding="utf-8")
            self.assertIn("# Phase Diagram", report)
            self.assertIn("Best cell", report)
            self.assertIn("success_rate", report)
            self.assertIn("Causet embeddability phase diagram", heatmap_svg.read_text(encoding="utf-8"))

    def test_phase_refine_explicit_cells_integration(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir, contextlib.redirect_stdout(io.StringIO()):
            tmp = Path(tmpdir)
            runs_csv = tmp / "refine_runs.csv"
            summary_csv = tmp / "refine_summary.csv"
            report_md = tmp / "refine.md"
            heatmap_svg = tmp / "refine.svg"

            rc = phase_refine.main(
                [
                    "--cells",
                    "4:1,4:2",
                    "--seed-start",
                    "1959",
                    "--seed-count",
                    "1",
                    "--warmup-limit",
                    "10",
                    "--anneal-limit",
                    "10",
                    "--max-data",
                    "3",
                    "--max-run-seconds",
                    "0",
                    "--runs-csv",
                    str(runs_csv),
                    "--summary-csv",
                    str(summary_csv),
                    "--report-md",
                    str(report_md),
                    "--heatmap-svg",
                    str(heatmap_svg),
                ]
            )

            self.assertEqual(rc, 0)
            self.assertTrue(runs_csv.exists())
            self.assertTrue(summary_csv.exists())
            self.assertTrue(report_md.exists())
            self.assertTrue(heatmap_svg.exists())

            run_rows = load_csv_rows(runs_csv)
            summary_rows = load_csv_rows(summary_csv)
            self.assertEqual(len(run_rows), 2)
            self.assertEqual(len(summary_rows), 2)
            self.assertEqual({row["dim"] for row in summary_rows}, {"1", "2"})
            self.assertTrue(all(row["n"] == "4" for row in summary_rows))
            self.assertTrue(all(row["status"] == "ok" for row in run_rows))
            self.assertIn("# Phase Diagram", report_md.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
