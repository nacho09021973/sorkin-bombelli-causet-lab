"""Regression tests for S4-KERR-K13-NEAR-PHOTON-WINDING-PROBE-001."""

from __future__ import annotations

import csv
import json
import math
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
KERR_DIR = ROOT / "explore" / "sorkin4_kerr_benchmark"
PREFIX = "kerr_k13_near_photon_winding_probe_001_n12_seed1959"

CSV_PATH = KERR_DIR / f"{PREFIX}.csv"
JSON_PATH = KERR_DIR / f"{PREFIX}.json"
MD_PATH = KERR_DIR / f"{PREFIX}.md"
PNG_PATH = KERR_DIR / f"{PREFIX}.png"
README_PATH = KERR_DIR / "README.md"


def _load_csv() -> list[dict[str, str]]:
    with CSV_PATH.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _load_json() -> dict:
    return json.loads(JSON_PATH.read_text(encoding="utf-8"))


def _to_bool(v: str | bool) -> bool:
    if isinstance(v, bool):
        return v
    return v.strip().lower() == "true"


class KerrK13NearPhotonProbeTests(unittest.TestCase):
    def test_artifacts_exist(self):
        assert CSV_PATH.exists(), f"Missing {CSV_PATH}"
        assert JSON_PATH.exists(), f"Missing {JSON_PATH}"
        assert MD_PATH.exists(), f"Missing {MD_PATH}"
        assert PNG_PATH.exists(), f"Missing {PNG_PATH}"

    def test_csv_json_parse(self):
        assert len(_load_csv()) > 0
        assert "cases" in _load_json()

    def test_png_nonzero(self):
        assert PNG_PATH.stat().st_size > 0

    def test_summary_all_checks_pass(self):
        assert _load_json()["global_summary"]["all_checks_pass"] is True

    def test_scanned_and_recorded_positive(self):
        s = _load_json()["global_summary"]
        assert s["total_candidates_scanned"] > 0
        assert s["total_cases_recorded"] > 0

    def test_max_abs_delta_phi_finite_positive(self):
        m = _load_json()["global_summary"]["max_abs_delta_phi"]
        assert math.isfinite(m)
        assert m > 0.0

    def test_non_advisory_non_unresolved_pass(self):
        for c in _load_json()["cases"]:
            if not c["advisory_only"] and not c["unresolved"]:
                assert c["all_checks_pass"] is True, c["case_id"]

    def test_accepted_invariants(self):
        for c in _load_json()["cases"]:
            if c["all_checks_pass"]:
                assert c["all_points_exterior"] is True, c["case_id"]
                assert c["null_condition_pass"] is True, c["case_id"]
                assert c["constants_consistency_pass"] is True, c["case_id"]

    def test_no_sprinkling_no_causal_classifier_flags(self):
        for c in _load_json()["cases"]:
            assert c["no_sprinkling_pair_used"] is True, c["case_id"]
            assert c["no_global_causal_relations_decided"] is True, c["case_id"]
            assert c["no_causal_classifier_introduced"] is True, c["case_id"]

    def test_positive_spin_global_accounting(self):
        s = _load_json()["global_summary"]
        assert s["global_true_relations"] == 0
        assert s["global_false_relations"] == 0
        assert s["global_undecided_pairs"] == 66

    def test_caveats_in_md_and_readme(self):
        text = (
            MD_PATH.read_text(encoding="utf-8").lower()
            + "\n"
            + README_PATH.read_text(encoding="utf-8").lower()
        )
        required = [
            "does not classify causal reachability",
            "does not use sprinkling event pairs",
            "does not implement a production kerr causal classifier",
            "not endpoint reachability",
        ]
        for needle in required:
            assert needle in text, needle

    def test_csv_flags(self):
        for r in _load_csv():
            assert _to_bool(r["no_sprinkling_pair_used"]), r["case_id"]
            assert _to_bool(r["no_global_causal_relations_decided"]), r["case_id"]
            assert _to_bool(r["no_causal_classifier_introduced"]), r["case_id"]
