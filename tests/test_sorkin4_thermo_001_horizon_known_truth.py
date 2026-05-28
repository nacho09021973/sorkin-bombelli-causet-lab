"""Regression tests for S4-THERMO-001-HORIZON-KNOWN-TRUTH-AUDIT."""

from __future__ import annotations

import csv
import json
import math
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
THERMO_DIR = ROOT / "explore" / "sorkin4_thermo_benchmark"
PREFIX = "s4_thermo_001_horizon_known_truth"

CSV_PATH = THERMO_DIR / f"{PREFIX}.csv"
JSON_PATH = THERMO_DIR / f"{PREFIX}.json"
MD_PATH = THERMO_DIR / f"{PREFIX}.md"
PNG_PATH = THERMO_DIR / f"{PREFIX}.png"
README_PATH = THERMO_DIR / "README.md"
EXPECTED_CHI = [0.0, 0.25, 0.5, 0.75, 0.9, 0.99]


def _load_csv() -> list[dict[str, str]]:
    with CSV_PATH.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _load_json() -> dict:
    return json.loads(JSON_PATH.read_text(encoding="utf-8"))


class Thermo001KnownTruthTests(unittest.TestCase):
    def test_artifacts_exist(self):
        assert CSV_PATH.exists(), f"Missing {CSV_PATH}"
        assert JSON_PATH.exists(), f"Missing {JSON_PATH}"
        assert MD_PATH.exists(), f"Missing {MD_PATH}"
        assert PNG_PATH.exists(), f"Missing {PNG_PATH}"

    def test_csv_json_parseable(self):
        rows = _load_csv()
        payload = _load_json()
        assert len(rows) > 0
        assert "rows" in payload

    def test_png_nonzero_size(self):
        assert PNG_PATH.stat().st_size > 0

    def test_summary_all_checks_pass(self):
        assert _load_json()["global_summary"]["all_checks_pass"] is True

    def test_chi_list_exact(self):
        chi = [r["chi"] for r in _load_json()["rows"]]
        assert chi == EXPECTED_CHI, chi

    def test_chi0_schwarzschild(self):
        row0 = _load_json()["rows"][0]
        m = row0["M"]
        tol = 1.0e-12
        assert abs(row0["r_plus"] - 2.0 * m) <= tol
        assert abs(row0["area"] - 16.0 * math.pi * m * m) <= tol
        assert abs(row0["kappa"] - 1.0 / (4.0 * m)) <= tol
        assert abs(row0["T_H"] - 1.0 / (8.0 * math.pi * m)) <= tol
        assert abs(row0["entropy_BH"] - 4.0 * math.pi * m * m) <= tol
        assert abs(row0["Omega_H"]) <= tol

    def test_area_monotonic_decrease(self):
        area = [r["area"] for r in _load_json()["rows"]]
        assert all(area[i + 1] < area[i] for i in range(len(area) - 1))

    def test_entropy_monotonic_decrease(self):
        ent = [r["entropy_BH"] for r in _load_json()["rows"]]
        assert all(ent[i + 1] < ent[i] for i in range(len(ent) - 1))

    def test_kappa_monotonic_decrease(self):
        kap = [r["kappa"] for r in _load_json()["rows"]]
        assert all(kap[i + 1] < kap[i] for i in range(len(kap) - 1))

    def test_temperature_monotonic_decrease(self):
        temp = [r["T_H"] for r in _load_json()["rows"]]
        assert all(temp[i + 1] < temp[i] for i in range(len(temp) - 1))

    def test_omega_monotonic_increase(self):
        omg = [r["Omega_H"] for r in _load_json()["rows"]]
        assert all(omg[i + 1] > omg[i] for i in range(len(omg) - 1))

    def test_chi099_extremal_trend(self):
        rows = _load_json()["rows"]
        r09 = rows[4]
        r099 = rows[5]
        area_ext = 8.0 * math.pi
        assert r099["kappa"] < r09["kappa"]
        assert r099["T_H"] < r09["T_H"]
        assert abs(r099["area"] - area_ext) < abs(r09["area"] - area_ext)

    def test_fixed_m_caveat_in_md_and_readme(self):
        text = (
            MD_PATH.read_text(encoding="utf-8").lower()
            + "\n"
            + README_PATH.read_text(encoding="utf-8").lower()
        )
        assert "fixed-m kerr spin sweeps compare stationary solutions" in text
        assert "not physical dynamical evolution" in text

    def test_level_a_flags(self):
        summary = _load_json()["global_summary"]
        assert summary["level_A_only"] is True
        assert summary["level_B_discrete_rediscovery_claimed"] is False
