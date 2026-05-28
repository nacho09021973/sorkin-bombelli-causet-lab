"""Regression tests for S4-KERR-K9-EQUATORIAL-FULL-RHS-PREFLIGHT-001."""

from __future__ import annotations

import csv
import json
import math
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
KERR_DIR = ROOT / "explore" / "sorkin4_kerr_benchmark"
PREFIX = "kerr_k9_equatorial_full_rhs_preflight_001_n12_seed1959"

CSV_PATH = KERR_DIR / f"{PREFIX}.csv"
JSON_PATH = KERR_DIR / f"{PREFIX}.json"
MD_PATH = KERR_DIR / f"{PREFIX}.md"
PNG_PATH = KERR_DIR / f"{PREFIX}.png"
README_PATH = KERR_DIR / "README.md"

N_EVENTS = 12
EXPECTED_PAIRS = N_EVENTS * (N_EVENTS - 1) // 2


def _to_bool(value: str | bool) -> bool:
    if isinstance(value, bool):
        return value
    return value.strip().lower() == "true"


def _load_csv() -> list[dict[str, str]]:
    with CSV_PATH.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _load_json() -> dict:
    return json.loads(JSON_PATH.read_text(encoding="utf-8"))


class KerrK9ArtifactTests(unittest.TestCase):
    def test_artifacts_exist(self):
        assert CSV_PATH.exists(), f"Missing {CSV_PATH}"
        assert JSON_PATH.exists(), f"Missing {JSON_PATH}"
        assert MD_PATH.exists(), f"Missing {MD_PATH}"
        assert PNG_PATH.exists(), f"Missing {PNG_PATH}"

    def test_csv_json_parseable(self):
        rows = _load_csv()
        payload = _load_json()
        assert len(rows) > 0
        assert "cases" in payload

    def test_png_nonzero_size(self):
        assert PNG_PATH.stat().st_size > 0

    def test_summary_counts_and_pass(self):
        summary = _load_json()["global_summary"]
        assert summary["all_checks_pass"] is True
        assert summary["total_cases"] == 16
        assert summary["passed_cases"] == 16
        assert summary["failed_cases"] == 0
        assert summary["advisory_cases"] == 2

    def test_non_advisory_cases_all_checks_pass(self):
        for case in _load_json()["cases"]:
            if not case["advisory_only"]:
                assert case["all_checks_pass"] is True, case["case_id"]

    def test_all_cases_core_checks(self):
        for case in _load_json()["cases"]:
            assert case["all_points_exterior"] is True, case["case_id"]
            assert case["finite_rhs_all_steps"] is True, case["case_id"]
            assert case["finite_solution_all_steps"] is True, case["case_id"]
            assert case["null_condition_pass"] is True, case["case_id"]
            assert case["constants_consistency_pass"] is True, case["case_id"]

    def test_accepted_cases_monotonic_and_radial_rhs(self):
        for case in _load_json()["cases"]:
            if case["all_checks_pass"]:
                assert case["t_monotonic_future_pass"] is True, case["case_id"]
                assert case["radial_rhs_consistency_pass"] is True, case["case_id"]

    def test_a0_b0_schwarzschild_limit_pass(self):
        for case in _load_json()["cases"]:
            if math.isclose(case["spin_a"], 0.0) and math.isclose(case["b"], 0.0):
                assert case["schwarzschild_radial_limit_pass"] is True, case["case_id"]

    def test_a0_b0_residuals_match_expected_limits(self):
        for case in _load_json()["cases"]:
            if math.isclose(case["spin_a"], 0.0) and math.isclose(case["b"], 0.0):
                assert abs(case["max_abs_schwarzschild_radial_residual"]) <= 1.0e-6, case["case_id"]
                assert abs(case["max_abs_E_residual"]) <= 1.0e-7, case["case_id"]
                assert abs(case["max_abs_L_residual"]) <= 1.0e-7, case["case_id"]

    def test_circular_cases_drift_or_advisory(self):
        for case in _load_json()["cases"]:
            if "circular" in case["case_id"]:
                assert (
                    case["circular_orbit_radial_drift_pass"] is True
                    or case["advisory_only"] is True
                ), case["case_id"]

    def test_positive_spin_global_causal_accounting_controls(self):
        payload = _load_json()
        summary = payload["global_summary"]
        assert summary["global_true_relations"] == 0
        assert summary["global_false_relations"] == 0
        assert summary["global_undecided_pairs"] == EXPECTED_PAIRS
        for case in payload["cases"]:
            if case["spin_a"] > 0.0:
                assert case["all_checks_pass"] is True, case["case_id"]

    def test_md_and_readme_caveats(self):
        md = MD_PATH.read_text(encoding="utf-8").lower()
        readme = README_PATH.read_text(encoding="utf-8").lower()
        text = md + "\n" + readme
        required = [
            "not point-to-point shooting",
            "does not decide causal reachability",
            "does not classify sprinkled event pairs",
            "b=0",
            "safe radial-flow control",
            "not generic kerr null-geodesic families",
            "drift diagnostics, not stability claims",
            "not used as evidence of orbital stability",
        ]
        for needle in required:
            assert needle in text, needle

    def test_csv_and_json_case_count_match(self):
        rows = _load_csv()
        payload = _load_json()
        assert len(rows) == len(payload["cases"]) == 16

    def test_csv_flags_true_for_required_checks(self):
        for row in _load_csv():
            assert _to_bool(row["all_points_exterior"]), row["case_id"]
            assert _to_bool(row["finite_rhs_all_steps"]), row["case_id"]
            assert _to_bool(row["finite_solution_all_steps"]), row["case_id"]
            assert _to_bool(row["null_condition_pass"]), row["case_id"]
            assert _to_bool(row["constants_consistency_pass"]), row["case_id"]

