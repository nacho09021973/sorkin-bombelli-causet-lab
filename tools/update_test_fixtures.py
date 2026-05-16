#!/usr/bin/env python3
"""Regenerate frozen regression fixtures for the cones pipeline."""

from __future__ import annotations

import contextlib
import io
import re
import tempfile
from pathlib import Path

import analyze_sweep
import cones
import ensemble_scan
import dimension_sweep
import schedule_sweep


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "fixtures"


def run_silently(func, args):
    with contextlib.redirect_stdout(io.StringIO()):
        return func(args)


def copy_text(src: Path, dst: Path) -> None:
    dst.write_text(src.read_text(encoding="utf-8"), encoding="utf-8")


def copy_text_replacing(src: Path, dst: Path, old: str, new: str) -> None:
    dst.write_text(src.read_text(encoding="utf-8").replace(old, new), encoding="utf-8")


def copy_ensemble_report(src: Path, dst: Path, root: str, tmpdir: str) -> None:
    text = src.read_text(encoding="utf-8")
    text = text.replace(root, "")
    text = text.replace(tmpdir, "<tmp>")
    text = re.sub(r"/tmp/ensemble_scan_[^/]+", "<tmp>", text)
    dst.write_text(text, encoding="utf-8")


def main() -> int:
    FIXTURES.mkdir(parents=True, exist_ok=True)

    # cones.py CLI smoke fixture
    with tempfile.TemporaryDirectory(prefix="cones_fixture_") as tmpdir:
        tmp = Path(tmpdir)
        out = tmp / "cone.out"
        svg = tmp / "cone.svg"
        run_silently(
            cones.main,
            [
                str(ROOT / "benchmarks" / "tesis_like_6.in"),
                "--dim",
                "2",
                "--seed",
                "1959",
                "--output",
                str(out),
                "--plot",
                str(svg),
                "--max-data",
                "5",
            ],
        )
        copy_text(svg, FIXTURES / "cones_cli_smoke.svg")

    # sprinkle edge fixture
    with tempfile.TemporaryDirectory(prefix="cones_fixture_") as tmpdir:
        tmp = Path(tmpdir)
        out = tmp / "sprinkle.out"
        svg = tmp / "sprinkle.svg"
        causet = tmp / "causet.svg"
        generated = tmp / "sprinkle.in"
        run_silently(
            cones.main,
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
                str(out),
                "--plot",
                str(svg),
                "--causet-plot",
                str(causet),
                "--save-generated-input",
                str(generated),
                "--max-data",
                "5",
            ],
        )
        copy_text(causet, FIXTURES / "cones_cli_sprinkle_edge.svg")

    # sprinkle antichain fixture
    with tempfile.TemporaryDirectory(prefix="cones_fixture_") as tmpdir:
        tmp = Path(tmpdir)
        out = tmp / "sprinkle_no_edges.out"
        svg = tmp / "sprinkle_no_edges.svg"
        causet = tmp / "causet_no_edges.svg"
        run_silently(
            cones.main,
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
                str(out),
                "--plot",
                str(svg),
                "--causet-plot",
                str(causet),
                "--max-data",
                "5",
            ],
        )
        copy_text(causet, FIXTURES / "cones_cli_sprinkle_causet.svg")

    # schedule sweep fixture
    with tempfile.TemporaryDirectory(prefix="cones_fixture_") as tmpdir:
        tmp = Path(tmpdir)
        run_csv = tmp / "runs.csv"
        summary_csv = tmp / "summary.csv"
        report_md = tmp / "report.md"
        heatmap_svg = tmp / "heatmap.svg"
        run_silently(
            schedule_sweep.main,
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
            ],
        )
        copy_text(run_csv, FIXTURES / "schedule_sweep_1959.csv")
        copy_text(heatmap_svg, FIXTURES / "schedule_sweep.svg")

    # dimension sweep fixture
    with tempfile.TemporaryDirectory(prefix="cones_fixture_") as tmpdir:
        tmp = Path(tmpdir)
        output_csv = tmp / "dimension.csv"
        report_md = tmp / "dimension.md"
        heatmap_svg = tmp / "dimension.svg"
        run_silently(
            dimension_sweep.main,
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
            ],
        )
        copy_text(output_csv, FIXTURES / "dimension_sweep_1959.csv")
        copy_text(heatmap_svg, FIXTURES / "dimension_sweep.svg")

    # analyze sweep fixture
    with tempfile.TemporaryDirectory(prefix="cones_fixture_") as tmpdir:
        tmp = Path(tmpdir)
        csv_path = tmp / "sweep.csv"
        report_md = tmp / "analysis.md"
        plot_svg = tmp / "analysis.svg"
        copy_text(FIXTURES / "analyze_sweep_sample.csv", csv_path)
        run_silently(
            analyze_sweep.main,
            [str(csv_path), "--report-md", str(report_md), "--plot-svg", str(plot_svg)],
        )
        copy_text(plot_svg, FIXTURES / "analyze_sweep.svg")

    # ensemble scan fixture
    with tempfile.TemporaryDirectory(prefix="cones_fixture_") as tmpdir:
        tmp = Path(tmpdir)
        run_csv = tmp / "ensemble_runs.csv"
        summary_csv = tmp / "ensemble_summary.csv"
        report_md = tmp / "ensemble.md"
        heatmap_svg = tmp / "ensemble.svg"
        run_silently(
            ensemble_scan.main,
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
            ],
        )
        copy_text_replacing(run_csv, FIXTURES / "ensemble_scan_runs.csv", tmpdir, "<tmp>")
        copy_text(summary_csv, FIXTURES / "ensemble_scan_summary.csv")
        copy_ensemble_report(report_md, FIXTURES / "ensemble_scan.md", str(ROOT) + "/", tmpdir)
        copy_text(heatmap_svg, FIXTURES / "ensemble_scan.svg")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
