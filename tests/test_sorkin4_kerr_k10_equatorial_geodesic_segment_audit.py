"""Regression tests for S4-KERR-K10-EQUATORIAL-GEODESIC-SEGMENT-AUDIT-001."""

from __future__ import annotations

import csv
import json
import math
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
KERR_DIR = ROOT / "explore" / "sorkin4_kerr_benchmark"
PREFIX = "kerr_k10_equatorial_geodesic_segment_audit_001_n12_seed1959"

CSV_PATH = KERR_DIR / f"{PREFIX}.csv"
JSON_PATH = KERR_DIR / f"{PREFIX}.json"
MD_PATH = KERR_DIR / f"{PREFIX}.md"
PNG_PATH = KERR_DIR / f"{PREFIX}.png"
README_PATH = KERR_DIR / "README.md"


def _to_bool(value: str | bool) -> bool:
    if isinstance(value, bool):
        return value
    return value.strip().lower() == "true"


def _load_csv() -> list[dict[str, str]]:
    with CSV_PATH.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _load_json() -> dict:
    return json.loads(JSON_PATH.read_text(encoding="utf-8"))


class KerrK10ArtifactTests(unittest.TestCase):
    def test_artifacts_exist(self):
        assert CSV_PATH.exists(), f"Missing {CSV_PATH}"
        assert JSON_PATH.exists(), f"Missing {JSON_PATH}"
        assert MD_PATH.exists(), f"Missing {MD_PATH}"
        assert PNG_PATH.exists(), f"Missing {PNG_PATH}"

    def test_csv_json_parseable(self):
        rows = _load_csv()
        payload = _load_json()
        assert len(rows) > 0
        assert "segments" in payload

    def test_png_nonzero_size(self):
        assert PNG_PATH.stat().st_size > 0

    def test_summary_all_checks_pass(self):
        summary = _load_json()["global_summary"]
        assert summary["all_checks_pass"] is True

    def test_non_advisory_segments_pass(self):
        for seg in _load_json()["segments"]:
            if not seg.get("advisory_only", False):
                assert seg["all_checks_pass"] is True, seg["case_id"]

    def test_all_segments_core_checks(self):
        for seg in _load_json()["segments"]:
            assert seg["all_points_exterior"] is True, seg["case_id"]
            assert seg["finite_rhs_all_steps"] is True, seg["case_id"]
            assert seg["finite_solution_all_steps"] is True, seg["case_id"]
            assert seg["null_condition_pass"] is True, seg["case_id"]
            assert seg["constants_consistency_pass"] is True, seg["case_id"]

    def test_accepted_segments_monotonic_and_rhs(self):
        for seg in _load_json()["segments"]:
            if seg["all_checks_pass"]:
                assert seg["t_monotonic_future_pass"] is True, seg["case_id"]
                assert seg["radial_rhs_consistency_pass"] is True, seg["case_id"]

    def test_a0_b0_schwarzschild_limit_pass(self):
        for seg in _load_json()["segments"]:
            if math.isclose(seg["spin_a"], 0.0) and math.isclose(seg["b"], 0.0):
                assert seg["schwarzschild_radial_limit_pass"] is True, seg["case_id"]

    def test_delta_phi_delta_t_finite_for_accepted(self):
        for seg in _load_json()["segments"]:
            if seg["all_checks_pass"]:
                assert math.isfinite(seg["delta_phi"]), seg["case_id"]
                assert math.isfinite(seg["delta_t"]), seg["case_id"]

    def test_no_endpoint_targeting_pass(self):
        for seg in _load_json()["segments"]:
            assert seg["no_endpoint_targeting_pass"] is True, seg["case_id"]

    def test_positive_spin_global_causal_accounting(self):
        summary = _load_json()["global_summary"]
        assert summary["global_true_relations"] == 0
        assert summary["global_false_relations"] == 0
        assert summary["global_undecided_pairs"] == 66

    def test_md_and_readme_caveats(self):
        md = MD_PATH.read_text(encoding="utf-8").lower()
        readme = README_PATH.read_text(encoding="utf-8").lower()
        text = md + "\n" + readme
        required = [
            "does not do point-to-point shooting",
            "does not solve boundary-value problems",
            "does not decide causal reachability",
            "does not classify sprinkled event pairs",
        ]
        for needle in required:
            assert needle in text, needle

    def test_csv_required_flags_true(self):
        for row in _load_csv():
            assert _to_bool(row["all_points_exterior"]), row["case_id"]
            assert _to_bool(row["finite_rhs_all_steps"]), row["case_id"]
            assert _to_bool(row["finite_solution_all_steps"]), row["case_id"]
            assert _to_bool(row["null_condition_pass"]), row["case_id"]
            assert _to_bool(row["constants_consistency_pass"]), row["case_id"]
            assert _to_bool(row["no_endpoint_targeting_pass"]), row["case_id"]
