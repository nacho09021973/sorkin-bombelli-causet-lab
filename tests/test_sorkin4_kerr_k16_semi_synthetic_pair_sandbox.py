"""Regression tests for S4-KERR-K16-SEMI-SYNTHETIC-PAIR-SANDBOX-001."""

from __future__ import annotations

import csv
import json
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
KERR_DIR = ROOT / "explore" / "sorkin4_kerr_benchmark"
SCRIPT_PATH = KERR_DIR / "audit_kerr_k16_semi_synthetic_pair_sandbox_001.py"
PREFIX = "kerr_k16_semi_synthetic_pair_sandbox_001_n12_seed1959"

CSV_PATH = KERR_DIR / f"{PREFIX}.csv"
JSON_PATH = KERR_DIR / f"{PREFIX}.json"
MD_PATH = KERR_DIR / f"{PREFIX}.md"
PNG_PATH = KERR_DIR / f"{PREFIX}.png"


def _load_csv() -> list[dict[str, str]]:
    with CSV_PATH.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _load_json() -> dict:
    return json.loads(JSON_PATH.read_text(encoding="utf-8"))


def _to_bool(v: str | bool) -> bool:
    if isinstance(v, bool):
        return v
    return v.strip().lower() == "true"


class KerrK16SemiSyntheticPairSandboxTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run(["python3", str(SCRIPT_PATH)], check=True, cwd=str(ROOT))

    def test_artifacts_exist(self):
        assert CSV_PATH.exists()
        assert JSON_PATH.exists()
        assert MD_PATH.exists()
        assert PNG_PATH.exists()

    def test_csv_json_parse(self):
        assert len(_load_csv()) > 0
        assert "cases" in _load_json()

    def test_png_nonzero(self):
        assert PNG_PATH.stat().st_size > 0

    def test_summary_all_checks_pass(self):
        assert _load_json()["global_summary"]["all_checks_pass"] is True

    def test_summary_counts(self):
        s = _load_json()["global_summary"]
        assert s["total_pairs"] > 0
        assert s["cloud_events_available"] > 0
        assert s["event_A_from_cloud_count"] > 0
        assert s["semi_synthetic_null_connected_pairs"] > 0
        assert s["semi_synthetic_no_match_pairs"] > 0
        assert s["recovered_known_answer_pairs"] > 0
        assert s["negative_controls_passed"] > 0

    def test_null_connected_rows(self):
        for c in _load_json()["cases"]:
            if (not c["advisory_only"]) and (not c["unresolved"]) and c["semi_synthetic_pair_classification"] == "semi_synthetic_null_connected":
                assert c["event_A_from_cloud"] is True, c["case_id"]
                assert c["event_B_forward_generated"] is True, c["case_id"]
                assert c["semi_synthetic_pair_recovered"] is True, c["case_id"]
                assert c["known_answer_null_connection_recovered"] is True, c["case_id"]
                assert c["correct_sector_recovered"] is True, c["case_id"]
                assert c["all_checks_pass"] is True, c["case_id"]

    def test_no_match_rows(self):
        for c in _load_json()["cases"]:
            if c["semi_synthetic_pair_classification"] == "semi_synthetic_no_match":
                assert c["semi_synthetic_pair_recovered"] is False, c["case_id"]
                assert c["known_answer_null_connection_recovered"] is False, c["case_id"]

    def test_guardrail_flags_all_rows(self):
        for c in _load_json()["cases"]:
            assert c["no_arbitrary_pair_used"] is True, c["case_id"]
            assert c["no_sprinkling_pair_reachability_claimed"] is True, c["case_id"]
            assert c["no_global_causal_relations_decided"] is True, c["case_id"]
            assert c["no_production_classifier_introduced"] is True, c["case_id"]

    def test_accepted_invariants(self):
        for c in _load_json()["cases"]:
            if c["semi_synthetic_pair_recovered"]:
                assert c["all_points_exterior"] is True, c["case_id"]
                assert c["null_condition_pass"] is True, c["case_id"]
                assert c["constants_consistency_pass"] is True, c["case_id"]

    def test_md_caveats(self):
        text = MD_PATH.read_text(encoding="utf-8").lower()
        required = [
            "semi-synthetic",
            "a comes from a deterministic event cloud, b is forward-generated",
            "does not classify arbitrary event pairs",
            "does not implement a production causal classifier",
            "not a physical/global kerr causal claim",
        ]
        for needle in required:
            assert needle in text, needle

    def test_global_causal_accounting(self):
        s = _load_json()["global_summary"]
        assert s["global_true_relations"] == 0
        assert s["global_false_relations"] == 0
        assert s["global_undecided_pairs"] == 66

    def test_csv_guardrails(self):
        for r in _load_csv():
            assert _to_bool(r["no_arbitrary_pair_used"]), r["case_id"]
            assert _to_bool(r["no_sprinkling_pair_reachability_claimed"]), r["case_id"]
            assert _to_bool(r["no_global_causal_relations_decided"]), r["case_id"]
            assert _to_bool(r["no_production_classifier_introduced"]), r["case_id"]


if __name__ == "__main__":
    unittest.main()
