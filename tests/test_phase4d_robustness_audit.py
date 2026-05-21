"""Regression tests for Phase 4D robustness-vs-invariants audit."""

from __future__ import annotations

import csv
import math
import unittest
from pathlib import Path

import validation_suite as vs

from tools import build_phase4a_epsilon_sweep as p4a
from tools import build_phase4d_robustness_audit as p4d


ROOT = Path(__file__).resolve().parents[1]
PHASE4C_PER_RUN_CSV = (
    ROOT / "benchmarks" / "foundation" / "phase4c_optimizer_seed_probe_per_run.csv"
)
PHASE4C_PER_CELL_EPS_CSV = (
    ROOT / "benchmarks" / "foundation" /
    "phase4c_optimizer_seed_probe_per_cell_epsilon.csv"
)
PHASE4D_PER_SEED_CSV = ROOT / "benchmarks" / "foundation" / "phase4d_robustness_per_seed.csv"
PHASE4D_PER_CELL_CSV = ROOT / "benchmarks" / "foundation" / "phase4d_robustness_per_cell.csv"
PHASE4D_MD           = ROOT / "benchmarks" / "foundation" / "phase4d_robustness_audit.md"


# ---------------------------------------------------------------------------
# Spearman / Pearson sanity
# ---------------------------------------------------------------------------

class CorrelationPrimitiveTests(unittest.TestCase):
    def test_pearson_perfect_positive(self) -> None:
        rho, n = p4d.pearson_r([1.0, 2.0, 3.0, 4.0, 5.0],
                               [10.0, 20.0, 30.0, 40.0, 50.0])
        self.assertEqual(n, 5)
        self.assertAlmostEqual(rho, 1.0, places=10)

    def test_pearson_perfect_negative(self) -> None:
        rho, _ = p4d.pearson_r([1.0, 2.0, 3.0], [3.0, 2.0, 1.0])
        self.assertAlmostEqual(rho, -1.0, places=10)

    def test_spearman_perfect_monotonic_not_linear(self) -> None:
        # y = exp(x) is monotonic-increasing → Spearman = 1, Pearson < 1
        xs = [1.0, 2.0, 3.0, 4.0, 5.0]
        ys = [math.exp(x) for x in xs]
        rho, _ = p4d.spearman_rho(xs, ys)
        self.assertAlmostEqual(rho, 1.0, places=10)
        pr, _ = p4d.pearson_r(xs, ys)
        self.assertLess(pr, 1.0)

    def test_spearman_with_ties(self) -> None:
        rho, n = p4d.spearman_rho([1.0, 2.0, 2.0, 3.0], [1.0, 2.0, 2.0, 3.0])
        self.assertEqual(n, 4)
        self.assertAlmostEqual(rho, 1.0, places=10)

    def test_pairwise_complete_drops_nan(self) -> None:
        rho, n = p4d.spearman_rho(
            [1.0, float("nan"), 3.0, 4.0],
            [1.0, 2.0, 3.0, float("inf")],
        )
        self.assertEqual(n, 2)

    def test_constant_input_gives_zero(self) -> None:
        rho, _ = p4d.pearson_r([1.0, 1.0, 1.0], [1.0, 2.0, 3.0])
        self.assertEqual(rho, 0.0)
        rho_s, _ = p4d.spearman_rho([1.0, 1.0, 1.0], [1.0, 2.0, 3.0])
        self.assertEqual(rho_s, 0.0)


# ---------------------------------------------------------------------------
# Verdict logic on synthetic correlation matrices.
# ---------------------------------------------------------------------------

class VerdictLogicTests(unittest.TestCase):
    def _corr(self, pairs: dict[tuple[str, str], float]) -> dict:
        out: dict[tuple[str, str], dict] = {}
        for (inv, tgt), rho in pairs.items():
            out[(inv, tgt)] = {"spearman": rho, "pearson": rho, "n": 90}
        return out

    def test_synthetic_zero_correlation_yields_no_robust(self) -> None:
        invs    = ("inv_a", "inv_b", "inv_c")
        targets = ("tgt_x", "tgt_y")
        pairs   = {(i, t): 0.05 * (1 if "_a" in i else -1) for i in invs for t in targets}
        verdict, info = p4d.compute_verdict(
            self._corr(pairs), invs, targets,
        )
        self.assertEqual(verdict, "NO_ROBUST_ORDER_THEORETIC_CORRELATE")
        self.assertLess(info["top_abs_spearman"], p4d.THRESHOLD_WEAK)

    def test_synthetic_single_invariant_above_weak_only_yields_weak(self) -> None:
        invs    = ("inv_a", "inv_b", "inv_c")
        targets = ("tgt_x", "tgt_y")
        pairs   = {(i, t): 0.10 for i in invs for t in targets}
        pairs[("inv_a", "tgt_x")] = 0.45
        verdict, info = p4d.compute_verdict(
            self._corr(pairs), invs, targets,
        )
        self.assertEqual(verdict, "WEAK_CORRELATE")
        self.assertAlmostEqual(info["top_abs_spearman"], 0.45, places=6)
        self.assertEqual(info["top_invariant"], "inv_a")
        self.assertEqual(info["top_target"], "tgt_x")

    def test_synthetic_two_invariants_above_detected_yields_detected(self) -> None:
        invs    = ("inv_a", "inv_b", "inv_c")
        targets = ("tgt_x", "tgt_y")
        pairs   = {(i, t): 0.10 for i in invs for t in targets}
        pairs[("inv_a", "tgt_x")] = 0.80
        pairs[("inv_b", "tgt_x")] = 0.70
        verdict, info = p4d.compute_verdict(
            self._corr(pairs), invs, targets,
        )
        self.assertEqual(verdict, "ORDER_THEORETIC_CORRELATE_DETECTED")
        self.assertEqual(info["detected_pair"], ("inv_a", "inv_b"))

    def test_synthetic_two_above_detected_opposite_signs_only_weak(self) -> None:
        invs    = ("inv_a", "inv_b", "inv_c")
        targets = ("tgt_x", "tgt_y")
        pairs   = {(i, t): 0.10 for i in invs for t in targets}
        pairs[("inv_a", "tgt_x")] = +0.80
        pairs[("inv_b", "tgt_x")] = -0.75
        verdict, _ = p4d.compute_verdict(
            self._corr(pairs), invs, targets,
        )
        # Same target, opposite sign → does not confirm a coherent coexistence.
        self.assertEqual(verdict, "WEAK_CORRELATE")

    def test_synthetic_two_high_but_different_targets_yields_weak(self) -> None:
        invs    = ("inv_a", "inv_b", "inv_c")
        targets = ("tgt_x", "tgt_y")
        pairs   = {(i, t): 0.10 for i in invs for t in targets}
        pairs[("inv_a", "tgt_x")] = +0.75
        pairs[("inv_b", "tgt_y")] = +0.70
        verdict, _ = p4d.compute_verdict(
            self._corr(pairs), invs, targets,
        )
        self.assertEqual(verdict, "WEAK_CORRELATE")


# ---------------------------------------------------------------------------
# Fixture-based regression tests (skip if Phase 4D not generated yet).
# ---------------------------------------------------------------------------

class Phase4DFixtureTests(unittest.TestCase):
    def test_per_seed_csv_schema_if_generated(self) -> None:
        if not PHASE4D_PER_SEED_CSV.exists():
            self.skipTest("Phase 4D per-seed CSV not generated yet; run make regen-phase4d")
        with PHASE4D_PER_SEED_CSV.open(newline="", encoding="utf-8") as fh:
            header = tuple(next(csv.reader(fh)))
        self.assertEqual(header, p4d.PER_SEED_HEADERS)

    def test_per_cell_csv_schema_if_generated(self) -> None:
        if not PHASE4D_PER_CELL_CSV.exists():
            self.skipTest("Phase 4D per-cell CSV not generated yet; run make regen-phase4d")
        with PHASE4D_PER_CELL_CSV.open(newline="", encoding="utf-8") as fh:
            header = tuple(next(csv.reader(fh)))
        self.assertEqual(header, p4d.PER_CELL_HEADERS)

    def test_per_seed_row_count_if_generated(self) -> None:
        if not PHASE4D_PER_SEED_CSV.exists():
            self.skipTest("Phase 4D per-seed CSV not generated yet; run make regen-phase4d")
        with PHASE4D_PER_SEED_CSV.open(newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        # Expect 9 (n,d) cells × 10 causet seeds = 90
        self.assertEqual(len(rows), 90)

    def test_per_cell_row_count_if_generated(self) -> None:
        if not PHASE4D_PER_CELL_CSV.exists():
            self.skipTest("Phase 4D per-cell CSV not generated yet; run make regen-phase4d")
        with PHASE4D_PER_CELL_CSV.open(newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        self.assertEqual(len(rows), 9)

    def test_label_stability_matches_phase4c_if_generated(self) -> None:
        if not (PHASE4D_PER_CELL_CSV.exists() and PHASE4C_PER_CELL_EPS_CSV.exists()):
            self.skipTest("Phase 4C/4D CSVs not generated yet")
        with PHASE4D_PER_CELL_CSV.open(newline="", encoding="utf-8") as fh:
            p4d_rows = list(csv.DictReader(fh))
        with PHASE4C_PER_CELL_EPS_CSV.open(newline="", encoding="utf-8") as fh:
            p4c_rows = list(csv.DictReader(fh))
        p4c_stab = {
            (int(r["n"]), int(r["target_dim"])): float(r["label_stability_cell"])
            for r in p4c_rows
        }
        for r in p4d_rows:
            key = (int(r["n"]), int(r["target_dim"]))
            self.assertIn(key, p4c_stab)
            self.assertAlmostEqual(
                float(r["label_stability_cell"]), p4c_stab[key], places=12,
                msg=f"label_stability mismatch for {key}",
            )

    def test_invariants_match_phase4a_compute_invariants_if_generated(self) -> None:
        """Regression check: invariants reported in Phase 4D per-seed agree
        with `p4a.compute_invariants` for the same (n, d, causet_seed)."""
        if not PHASE4D_PER_SEED_CSV.exists():
            self.skipTest("Phase 4D per-seed CSV not generated yet; run make regen-phase4d")
        with PHASE4D_PER_SEED_CSV.open(newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        # Spot-check a handful of rows across the grid.
        sample_keys = [(32, 2, 1900), (48, 3, 1916), (64, 4, 2020)]
        for (n, d, cseed) in sample_keys:
            row = next(
                (r for r in rows
                 if int(r["n"]) == n
                 and int(r["target_dim"]) == d
                 and int(r["causet_seed"]) == cseed),
                None,
            )
            self.assertIsNotNone(row, msg=f"missing row for {(n, d, cseed)}")
            matrix, _ = vs.sprinkle_minkowski_diamond(n=n, seed=cseed, d_spacetime=d)
            expected = p4a.compute_invariants(matrix, n)
            for inv in ("ordering_fraction", "chain3_abundance", "link_count", "height"):
                csv_val = float(row[inv])
                exp_val = float(expected[inv])
                self.assertAlmostEqual(
                    csv_val, exp_val, places=8,
                    msg=f"invariant {inv} mismatch for {(n,d,cseed)}",
                )

    def test_markdown_present_if_generated(self) -> None:
        if not PHASE4D_MD.exists():
            self.skipTest("Phase 4D markdown not generated yet; run make regen-phase4d")
        text = PHASE4D_MD.read_text(encoding="utf-8")
        for section in (
            "## Semantic caveats",
            "## Objective",
            "## Method",
            "## Targets",
            "## Verdict",
            "## Correlation matrix",
            "## Interpretation",
            "## Scope",
        ):
            self.assertIn(section, text, msg=f"missing section {section}")

    def test_markdown_has_no_obsolete_physical_language_if_generated(self) -> None:
        if not PHASE4D_MD.exists():
            self.skipTest("Phase 4D markdown not generated yet; run make regen-phase4d")
        text = PHASE4D_MD.read_text(encoding="utf-8")
        # The MD intentionally cites the forbidden words inside the caveats
        # section to negate them; we instead check that the *verdict text*
        # never affirms any of these.  Sentinel: the conservative conclusion
        # phrase, if NO_ROBUST, is the exact quote from the user's criterion.
        forbidden_assertions = (
            "physical transition",
            "manifoldlike evidence",
            "embeddability evidence",
            "new result in CST",
        )
        for phrase in forbidden_assertions:
            self.assertNotIn(phrase, text, msg=f"obsolete claim: {phrase!r}")


# ---------------------------------------------------------------------------
# n-control helper unit tests (Step 4 of phase4d_ncontrol_prompt.md)
# ---------------------------------------------------------------------------

class NControlHelperTests(unittest.TestCase):

    def _make_rows(self, n_values, x_vals, y_vals):
        return [{"n": n, "x": x, "y": y}
                for n, x, y in zip(n_values, x_vals, y_vals)]

    def test_partial_spearman_collapses_when_driven_by_n(self) -> None:
        """Pure-n: x = n + tiny noise, y = n + tiny noise → partial ρ(x,y|n) ≈ 0."""
        import random
        rng = random.Random(42)
        n_vals = [32] * 30 + [48] * 30 + [64] * 30
        x_vals = [n + rng.uniform(-0.1, 0.1) for n in n_vals]
        y_vals = [n + rng.uniform(-0.1, 0.1) for n in n_vals]
        rows = self._make_rows(n_vals, x_vals, y_vals)

        rho_raw, _ = p4d.spearman_rho(x_vals, y_vals)
        self.assertGreater(abs(rho_raw), 0.9, "raw Spearman must be large when both track n")

        rho_partial, method = p4d.partial_spearman_rho_n(rows, "x", "y", "n")
        # Must collapse substantially: raw |ρ| > 0.9 but partial should be < 0.3
        self.assertLess(abs(rho_partial), 0.3,
                        msg="partial Spearman must collapse when x,y are driven by n")
        self.assertEqual(method, p4d.N_CONTROL_METHOD)

    def test_partial_spearman_survives_genuine_signal(self) -> None:
        """Genuine signal: y = x + tiny noise, x independent of n → partial ρ stays large."""
        import random
        rng = random.Random(99)
        n_vals = [32] * 30 + [48] * 30 + [64] * 30
        x_vals = [rng.uniform(0.0, 1.0) for _ in n_vals]
        y_vals = [x + rng.uniform(-0.05, 0.05) for x in x_vals]
        rows = self._make_rows(n_vals, x_vals, y_vals)

        rho_raw, _ = p4d.spearman_rho(x_vals, y_vals)
        self.assertGreater(abs(rho_raw), 0.9)

        rho_partial, _ = p4d.partial_spearman_rho_n(rows, "x", "y", "n")
        self.assertGreater(abs(rho_partial), 0.8,
                           msg="partial Spearman must survive when signal is not n-driven")

    def test_stratified_collapses_when_driven_by_n(self) -> None:
        """Pure-n: within each stratum x,y are conditionally independent → min_abs ≈ 0."""
        import random
        rng = random.Random(42)
        n_vals = [32] * 30 + [48] * 30 + [64] * 30
        x_vals = [n + rng.uniform(-0.1, 0.1) for n in n_vals]
        y_vals = [n + rng.uniform(-0.1, 0.1) for n in n_vals]
        rows = self._make_rows(n_vals, x_vals, y_vals)

        strat_dict, min_abs = p4d.stratified_spearman_by_n(
            rows, "x", "y", (32, 48, 64), "n"
        )
        self.assertLess(min_abs, 0.3,
                        msg="min_abs must collapse when correlation is purely n-driven")
        for n_val in (32, 48, 64):
            self.assertIn(n_val, strat_dict)

    def test_stratified_survives_genuine_signal(self) -> None:
        """Genuine signal: within each stratum y tracks x → min_abs stays large."""
        import random
        rng = random.Random(99)
        n_vals = [32] * 30 + [48] * 30 + [64] * 30
        x_vals = [rng.uniform(0.0, 1.0) for _ in n_vals]
        y_vals = [x + rng.uniform(-0.05, 0.05) for x in x_vals]
        rows = self._make_rows(n_vals, x_vals, y_vals)

        strat_dict, min_abs = p4d.stratified_spearman_by_n(
            rows, "x", "y", (32, 48, 64), "n"
        )
        self.assertGreater(min_abs, 0.8,
                           msg="min_abs must stay large when signal is genuine")


class Phase4DNControlFixtureTests(unittest.TestCase):
    """Fixture tests for n-control outputs (skip if not yet generated)."""

    def test_per_seed_csv_has_ncontrol_columns_if_generated(self) -> None:
        if not PHASE4D_PER_SEED_CSV.exists():
            self.skipTest("Phase 4D per-seed CSV not generated yet; run make regen-phase4d")
        with PHASE4D_PER_SEED_CSV.open(newline="", encoding="utf-8") as fh:
            header = tuple(next(csv.reader(fh)))
        self.assertIn("nctrl_method", header, "nctrl_method column missing from per-seed CSV")
        partial_cols = [c for c in header if c.startswith("nctrl_partial__")]
        self.assertGreater(len(partial_cols), 0, "no nctrl_partial__ columns found")
        minabs_cols = [c for c in header if c.startswith("nctrl_minabs__")]
        self.assertGreater(len(minabs_cols), 0, "no nctrl_minabs__ columns found")

    def test_per_seed_csv_schema_matches_declared_headers_if_generated(self) -> None:
        if not PHASE4D_PER_SEED_CSV.exists():
            self.skipTest("Phase 4D per-seed CSV not generated yet; run make regen-phase4d")
        with PHASE4D_PER_SEED_CSV.open(newline="", encoding="utf-8") as fh:
            header = tuple(next(csv.reader(fh)))
        self.assertEqual(header, p4d.PER_SEED_HEADERS,
                         "per-seed CSV header does not match PER_SEED_HEADERS")

    def test_per_cell_csv_schema_unchanged_if_generated(self) -> None:
        if not PHASE4D_PER_CELL_CSV.exists():
            self.skipTest("Phase 4D per-cell CSV not generated yet; run make regen-phase4d")
        with PHASE4D_PER_CELL_CSV.open(newline="", encoding="utf-8") as fh:
            header = tuple(next(csv.reader(fh)))
        self.assertEqual(header, p4d.PER_CELL_HEADERS,
                         "per-cell CSV must be byte-equivalent: no nctrl columns expected")
        nctrl_cols = [c for c in header if c.startswith("nctrl_")]
        self.assertEqual(nctrl_cols, [], "per-cell CSV must not have nctrl columns")

    def test_markdown_has_ncontrol_section_if_generated(self) -> None:
        if not PHASE4D_MD.exists():
            self.skipTest("Phase 4D markdown not generated yet; run make regen-phase4d")
        text = PHASE4D_MD.read_text(encoding="utf-8")
        self.assertIn("n-control (interpretive layer)", text)
        self.assertIn("interpretive auditing, not mechanical", text)
        self.assertIn("extensive size-like invariant", text)
        self.assertIn("No n-control is computed at per-cell level", text)

    def test_verdict_unchanged_after_ncontrol_if_generated(self) -> None:
        """The n-control must not alter the Phase 4D verdict."""
        if not PHASE4D_MD.exists():
            self.skipTest("Phase 4D markdown not generated yet; run make regen-phase4d")
        text = PHASE4D_MD.read_text(encoding="utf-8")
        # Verdict section must still carry the original label
        self.assertIn("ORDER_THEORETIC_CORRELATE_DETECTED", text)
        # No new DETECTED/NO_DETECTED label introduced
        self.assertNotIn("NCONTROL_DETECTED", text)
        self.assertNotIn("N_CONTROLLED_DETECTED", text)


if __name__ == "__main__":
    unittest.main()
