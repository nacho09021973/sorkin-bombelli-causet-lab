"""Tests for SORKIN-2 order-matrix diagnostic figures.

These tests use synthetic causal matrices and a temporary directory only.
They do not run the annealer, do not evaluate energy, and do not write
benchmark result data.
"""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import validation_suite as vs


def make_chain(n: int) -> list[list[bool]]:
    return [[i < j for j in range(n)] for i in range(n)]


def make_antichain(n: int) -> list[list[bool]]:
    return [[False] * n for _ in range(n)]


class TestOrderMatrixPlots(unittest.TestCase):
    def assert_png_exists_and_nonempty(self, path: Path) -> None:
        self.assertTrue(path.exists(), msg=f"{path} should exist")
        self.assertGreater(path.stat().st_size, 0, msg=f"{path} should be nonempty")
        self.assertEqual(path.read_bytes()[:8], b"\x89PNG\r\n\x1a\n")

    def test_writes_chain_vs_chain_figures(self) -> None:
        z = make_chain(4)
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = vs.write_order_matrix_plots(z, z, Path(tmpdir))

            self.assertEqual(paths.target_order_matrix.name, "target_order_matrix.png")
            self.assertEqual(paths.induced_order_matrix.name, "induced_order_matrix.png")
            self.assertEqual(paths.order_difference_matrix.name, "order_difference_matrix.png")
            self.assert_png_exists_and_nonempty(paths.target_order_matrix)
            self.assert_png_exists_and_nonempty(paths.induced_order_matrix)
            self.assert_png_exists_and_nonempty(paths.order_difference_matrix)

    def test_writes_chain_vs_antichain_difference_figure(self) -> None:
        target = make_chain(4)
        induced = make_antichain(4)
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = vs.write_order_matrix_plots(target, induced, Path(tmpdir))

            self.assert_png_exists_and_nonempty(paths.target_order_matrix)
            self.assert_png_exists_and_nonempty(paths.induced_order_matrix)
            self.assert_png_exists_and_nonempty(paths.order_difference_matrix)


if __name__ == "__main__":
    unittest.main()
