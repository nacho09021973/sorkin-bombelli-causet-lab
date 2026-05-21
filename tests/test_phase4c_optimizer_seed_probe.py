"""Regression tests for Phase 4C optimizer-seed multi-start probe."""

from __future__ import annotations

import csv
import math
import unittest
from pathlib import Path

import validation_suite as vs

from tools import build_phase4a_epsilon_sweep as p4a
from tools import build_phase4c_optimizer_seed_probe as p4c


ROOT = Path(__file__).resolve().parents[1]
PHASE4A_CSV = ROOT / "benchmarks" / "foundation" / "phase4a_epsilon_sweep.csv"
PHASE4C_PER_RUN_CSV = (
    ROOT / "benchmarks" / "foundation" / "phase4c_optimizer_seed_probe_per_run.csv"
)
PHASE4C_PER_CELL_CSV = (
    ROOT / "benchmarks" / "foundation" /
    "phase4c_optimizer_seed_probe_per_cell_epsilon.csv"
)


# ---------------------------------------------------------------------------
# Reproducibility: K=1 with optimizer_seed=1987 must match Phase 4A exactly.
# ---------------------------------------------------------------------------

class Phase4CReproducibilityTests(unittest.TestCase):
    def test_K1_reproduces_phase4a_exactly(self) -> None:
        if not PHASE4A_CSV.exists():
            self.skipTest("Phase 4A CSV not generated yet; run regen-phase4a first")
        rows = p4a.load_rows_from_csv(PHASE4A_CSV)
        valid = [r for r in rows if r.get("row_valid", False)][:5]
        self.assertGreater(len(valid), 0, "Phase 4A CSV must contain valid rows")
        for row in valid:
            n = row["n"]; d = row["target_dim"]
            seed = row["seed"]; eps = row["epsilon"]
            matrix, points = vs.sprinkle_minkowski_diamond(
                n=n, seed=seed, d_spacetime=d,
            )
            sim_result = p4a._run_one(
                d, n, seed, eps, matrix, points,
                optimizer_seed=p4a.OPTIMIZER_SEED,
            )
            E0  = sim_result["initial_energy"]
            wdE = sim_result["warmup_delta_energy"]
            if not math.isfinite(E0) or abs(E0) < p4a.INITIAL_ENERGY_FLOOR:
                computed = float("nan")
            else:
                computed = abs(wdE / E0)
            expected = row["abs_relative_drift"]
            if math.isnan(expected) and math.isnan(computed):
                continue
            # The Phase 4A CSV stores floats with `.4f` precision (half-unit = 5e-5).
            # The simulator output itself is deterministic.
            self.assertTrue(
                math.isclose(expected, computed, rel_tol=0.0, abs_tol=5e-5),
                msg=(
                    f"K=1 default seed must reproduce Phase 4A to CSV precision: "
                    f"row (n={n}, d={d}, seed={seed}, eps={eps}) "
                    f"expected={expected!r}, computed={computed!r}"
                ),
            )

    def test_verify_against_phase4a_helper_passes(self) -> None:
        if not PHASE4A_CSV.exists():
            self.skipTest("Phase 4A CSV not generated yet; run regen-phase4a first")
        ok, failures = p4c.verify_against_phase4a(
            optimizer_seed=p4a.OPTIMIZER_SEED,
            sample_size=6,
        )
        self.assertTrue(ok, msg="verify_against_phase4a failures:\n" + "\n".join(failures))

    def test_run_one_signature_keeps_default(self) -> None:
        """The optimizer_seed argument must default to OPTIMIZER_SEED so existing
        Phase 4A/4B callers preserve bit-exact behaviour."""
        import inspect
        sig = inspect.signature(p4a._run_one)
        self.assertIn("optimizer_seed", sig.parameters)
        self.assertEqual(
            sig.parameters["optimizer_seed"].default,
            p4a.OPTIMIZER_SEED,
        )


# ---------------------------------------------------------------------------
# Verdict logic on synthetic per_cell_eps rows.
# ---------------------------------------------------------------------------

def _make_cell_row(
    n: int, d: int, eps: float, mean_loss: float, iqr: float,
    floor_K: float, floor_p4a: float, label_stability: float,
    shapes: str = "v_shape|v_shape|v_shape",
) -> dict:
    return {
        "phase": "phase4c_optimizer_seed_probe",
        "n": n, "target_dim": d, "epsilon": eps, "K": 3,
        "n_valid_runs": 30,
        "mean_loss_K": mean_loss, "std_loss_K": iqr / 1.349,
        "min_loss_K": max(0.0, mean_loss - iqr),
        "max_loss_K": mean_loss + iqr,
        "iqr_loss_K": iqr,
        "floor_saturated_fraction_K": floor_K,
        "phase4a_mean_loss": mean_loss,
        "phase4a_floor_saturated_fraction": floor_p4a,
        "delta_min_loss_vs_phase4a": 0.0,
        "delta_floor_saturated_vs_phase4a": floor_K - floor_p4a,
        "curve_shape_per_optimizer_seed": shapes,
        "label_stability_cell": label_stability,
    }


class Phase4CVerdictTests(unittest.TestCase):
    def test_synthetic_high_variance_flagged(self) -> None:
        """High IQR/loss ratio with label flips → OPTIMIZER_SEED_LIMITED."""
        rows = [
            _make_cell_row(
                n=32, d=3, eps=eps,
                mean_loss=0.20, iqr=0.10,
                floor_K=0.0, floor_p4a=0.0,
                label_stability=0.0,
                shapes="v_shape|monotone_decay|noisy",
            )
            for eps in (0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.15, 0.20)
        ]
        verdict, info = p4c.compute_verdict(rows)
        self.assertEqual(verdict, "OPTIMIZER_SEED_LIMITED")
        self.assertEqual(info["mean_label_stability"], 0.0)
        self.assertGreater(info["iqr_ratio"], p4c.IQR_LIMIT_RATIO)

    def test_synthetic_floor_drop_flags_limited(self) -> None:
        """floor_K << floor_phase4a (but >0 baseline) → OPTIMIZER_SEED_LIMITED."""
        rows = [
            _make_cell_row(
                n=64, d=2, eps=eps,
                mean_loss=0.10, iqr=0.0001,
                floor_K=0.05, floor_p4a=0.80,
                label_stability=1.0,
                shapes="monotone_decay|monotone_decay|monotone_decay",
            )
            for eps in (0.04, 0.06, 0.08)
        ]
        verdict, info = p4c.compute_verdict(rows)
        self.assertEqual(verdict, "OPTIMIZER_SEED_LIMITED")
        self.assertTrue(info["floor_limited_triggered"])

    def test_synthetic_low_variance_flagged(self) -> None:
        """Tight IQR, stable labels, matching floor → OPTIMIZER_SEED_ROBUST."""
        rows = [
            _make_cell_row(
                n=32, d=3, eps=eps,
                mean_loss=0.10, iqr=0.0001,
                floor_K=0.20, floor_p4a=0.20,
                label_stability=1.0,
                shapes="v_shape|v_shape|v_shape",
            )
            for eps in (0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.15, 0.20)
        ]
        verdict, info = p4c.compute_verdict(rows)
        self.assertEqual(verdict, "OPTIMIZER_SEED_ROBUST")
        self.assertGreater(info["mean_label_stability"], p4c.LABEL_ROBUST)
        self.assertLess(info["iqr_ratio"], p4c.IQR_ROBUST_RATIO)
        self.assertTrue(info["floor_robust_satisfied"])

    def test_synthetic_inconclusive_when_mixed(self) -> None:
        """Stable labels but moderate IQR ratio → INCONCLUSIVE."""
        rows = [
            _make_cell_row(
                n=32, d=3, eps=eps,
                mean_loss=0.10, iqr=0.005,  # ratio = 0.05, between bands
                floor_K=0.20, floor_p4a=0.20,
                label_stability=1.0,
                shapes="v_shape|v_shape|v_shape",
            )
            for eps in (0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.15, 0.20)
        ]
        verdict, _ = p4c.compute_verdict(rows)
        self.assertEqual(verdict, "INCONCLUSIVE")


# ---------------------------------------------------------------------------
# Negative-control: target_dim=2 must never receive v_shape labels in Phase 4C
# under any of the K optimizer seeds when the underlying losses are the
# real Phase 4A d=2 curve (which is monotone_decay).
# ---------------------------------------------------------------------------

class Phase4CTargetDim2ControlTests(unittest.TestCase):
    def test_target_dim_2_controls_remain_negative(self) -> None:
        """A synthetic d=2 cell with the real Phase 4A (32,2) loss curve
        (replicated identically across K optimizer seeds) must classify as
        a non-v_shape label for every optimizer seed."""
        losses_32_2 = [0.361, 0.143, 0.0444, 0.00164, 0.0461, 0.0107, 0.00205, 0.0247]
        epsilons = p4a.EPSILONS
        per_run_rows: list[dict] = []
        for opt_seed in p4c.PHASE4C_OPTIMIZER_SEEDS:
            for causet_seed in p4a.PHASE4A_SEEDS[:3]:
                for eps, loss in zip(epsilons, losses_32_2):
                    per_run_rows.append({
                        "phase": "phase4c_optimizer_seed_probe",
                        "n": 32, "target_dim": 2,
                        "causet_seed": causet_seed,
                        "epsilon": eps,
                        "optimizer_seed": opt_seed,
                        "valid": True, "failure_mode": "",
                        "loss": loss,
                        "initial_energy": 1.0,
                        "final_energy": 0.0,
                        "warmup_delta_energy": -loss,
                        "warmup_attempted_moves": 10,
                        "warmup_accepted_moves": 5,
                        "warmup_rejected_moves": 5,
                    })
        per_cell_eps = p4c.summarize_per_cell_epsilon(
            per_run_rows,
            p4c.PHASE4C_OPTIMIZER_SEEDS,
            phase4a_baseline={},
        )
        d2_rows = [r for r in per_cell_eps if r["target_dim"] == 2]
        self.assertGreater(len(d2_rows), 0)
        for row in d2_rows:
            shapes = row["curve_shape_per_optimizer_seed"].split("|")
            for shape in shapes:
                self.assertNotEqual(
                    shape, "v_shape",
                    msg=f"d=2 control became v_shape: {row}",
                )


# ---------------------------------------------------------------------------
# CSV schema (fixtures present only after a probe run).
# ---------------------------------------------------------------------------

class Phase4CFixtureTests(unittest.TestCase):
    def test_per_run_csv_schema_if_generated(self) -> None:
        if not PHASE4C_PER_RUN_CSV.exists():
            self.skipTest("Phase 4C per-run CSV not generated yet; run regen-phase4c")
        with PHASE4C_PER_RUN_CSV.open(newline="", encoding="utf-8") as fh:
            header = tuple(next(csv.reader(fh)))
        self.assertEqual(header, p4c.PER_RUN_HEADERS)

    def test_per_cell_epsilon_csv_schema_if_generated(self) -> None:
        if not PHASE4C_PER_CELL_CSV.exists():
            self.skipTest(
                "Phase 4C per-cell-epsilon CSV not generated yet; run regen-phase4c"
            )
        with PHASE4C_PER_CELL_CSV.open(newline="", encoding="utf-8") as fh:
            header = tuple(next(csv.reader(fh)))
        self.assertEqual(header, p4c.PER_CELL_EPSILON_HEADERS)

    def test_per_run_csv_covers_pilot_grid_if_generated(self) -> None:
        if not PHASE4C_PER_RUN_CSV.exists():
            self.skipTest("Phase 4C per-run CSV not generated yet; run regen-phase4c")
        with PHASE4C_PER_RUN_CSV.open(newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        sizes = {int(r["n"]) for r in rows}
        dims  = {int(r["target_dim"]) for r in rows}
        seeds = {int(r["optimizer_seed"]) for r in rows}
        self.assertTrue(set(p4c.PILOT_SIZES).issubset(sizes))
        self.assertTrue(set(p4c.PILOT_DIMS).issubset(dims))
        self.assertIn(p4a.OPTIMIZER_SEED, seeds)


if __name__ == "__main__":
    unittest.main()
