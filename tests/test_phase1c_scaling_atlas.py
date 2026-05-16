"""Integrity tests for the Phase 1C two-control scaling atlas."""

from __future__ import annotations

import csv
import math
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ATLAS_CSV = ROOT / "benchmarks" / "foundation" / "phase1c_scaling_atlas.csv"
ATLAS_MD = ROOT / "benchmarks" / "foundation" / "phase1c_scaling_atlas.md"

EXPECTED_HEADERS = (
    "family",
    "d_spacetime",
    "n",
    "seed",
    "mm_dim",
    "midpoint_dim",
    "abs_discrepancy",
)

EXPECTED_FAMILIES = {
    "minkowski",
    "kleitman_rothschild",
    "corona_poset",
}
EXPECTED_SIZES = {32, 64, 128, 256}
EXPECTED_MINKOWSKI_DIMS = {"2", "3", "4"}
EXPECTED_SEEDS = {1959, 1962, 1987, 2009, 2026}


class Phase1CScalingAtlasIntegrityTests(unittest.TestCase):
    def setUp(self) -> None:
        self.assertTrue(
            ATLAS_CSV.exists(),
            msg=(
                f"missing atlas CSV at {ATLAS_CSV}; "
                "run `make regen-phase1c`"
            ),
        )
        with ATLAS_CSV.open(newline="", encoding="utf-8") as fh:
            self.rows = list(csv.DictReader(fh))

    def test_markdown_exists(self) -> None:
        self.assertTrue(
            ATLAS_MD.exists(),
            msg=(
                f"missing atlas markdown at {ATLAS_MD}; "
                "run `make regen-phase1c`"
            ),
        )

    def test_columns_match_expected(self) -> None:
        with ATLAS_CSV.open(newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh)
            header = tuple(next(reader))
        self.assertEqual(header, EXPECTED_HEADERS)

    def test_row_count_and_families(self) -> None:
        # 4 sizes x 5 seeds x (3 Minkowski d's + 2 controls) = 100 rows.
        self.assertEqual(len(self.rows), 100)
        self.assertEqual(
            {r["family"] for r in self.rows}, EXPECTED_FAMILIES
        )

    def test_covers_full_grid(self) -> None:
        minkowski_keys = {
            (r["d_spacetime"], int(r["n"]), int(r["seed"]))
            for r in self.rows if r["family"] == "minkowski"
        }
        expected_minkowski = {
            (d, n, s)
            for d in EXPECTED_MINKOWSKI_DIMS
            for n in EXPECTED_SIZES
            for s in EXPECTED_SEEDS
        }
        self.assertEqual(minkowski_keys, expected_minkowski)

        for family in ("kleitman_rothschild", "corona_poset"):
            control_keys = {
                (int(r["n"]), int(r["seed"]))
                for r in self.rows if r["family"] == family
            }
            expected = {(n, s) for n in EXPECTED_SIZES for s in EXPECTED_SEEDS}
            self.assertEqual(control_keys, expected)

    def test_no_nan_or_inf_in_dimension_columns(self) -> None:
        for row in self.rows:
            for col in ("mm_dim", "midpoint_dim"):
                raw = row[col]
                self.assertNotIn(
                    raw.lower(), {"nan", "inf", "-inf"},
                    msg=(
                        f"{col}={raw} for cell "
                        f"family={row['family']}, d={row['d_spacetime']}, "
                        f"n={row['n']}, seed={row['seed']}"
                    ),
                )
                value = float(raw)
                self.assertTrue(
                    math.isfinite(value),
                    msg=(
                        f"non-finite {col}={value} for cell "
                        f"family={row['family']}, n={row['n']}, "
                        f"seed={row['seed']}"
                    ),
                )


if __name__ == "__main__":
    unittest.main()
