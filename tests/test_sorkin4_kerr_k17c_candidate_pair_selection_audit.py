"""Regression tests for S4-KERR-K17C-CANDIDATE-PAIR-SELECTION-AUDIT-001."""

from __future__ import annotations

import csv
import json
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
KERR_DIR = ROOT / "explore" / "sorkin4_kerr_benchmark"
SCRIPT_PATH = KERR_DIR / "audit_kerr_k17c_candidate_pair_selection_audit_001.py"
PREFIX = "kerr_k17c_candidate_pair_selection_audit_001_n12_seed1959"

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


class KerrK17cCandidatePairSelectionAuditTests(unittest.TestCase):
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

    def test_pair_counts_positive(self):
        s = _load_json()["global_summary"]
        assert s["total_pairs_enumerated"] > 0
        assert s["forward_time_pairs"] > 0
        assert s["exterior_forward_pairs"] > 0

    def test_named_selection_outcomes_exist(self):
        rows = _load_json()["cases"]
        has_named = any(
            r["selection_label"].startswith("selected_")
            or r["selection_label"] in {"rejected_time_short", "rejected_angular_large", "rejected_radial_proxy"}
            for r in rows
        )
        assert has_named

    def test_no_causal_claim_flags(self):
        for r in _load_csv():
            assert _to_bool(r["no_causal_claim_introduced"]), r["case_id"]
            assert _to_bool(r["no_production_classifier_introduced"]), r["case_id"]
            assert _to_bool(r["no_global_causal_relations_decided"]), r["case_id"]

    def test_md_caveats(self):
        text = MD_PATH.read_text(encoding="utf-8").lower()
        assert "selection_candidate is not reachability" in text
        assert "rejected_by_selection is not proof of spacelike separation" in text


if __name__ == "__main__":
    unittest.main()
