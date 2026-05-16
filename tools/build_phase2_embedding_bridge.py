#!/usr/bin/env python3
"""Build the Phase 2 embedding bridge.

Phase 2 connects pre-embedding order-theoretic diagnostics to a small,
fixed run of the historical Bombelli-Sorkin annealing code. It is not
an optimizer search. The grid is deliberately tiny: one size, one case
seed, and five families from Phase 1D. The question is whether the
structural atlas gives useful warning signs before invoking the
embedding algorithm.
"""

from __future__ import annotations

import contextlib
import io
import math
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import cones  # noqa: E402
from tools.build_phase1_atlas import _format_field  # noqa: E402
from tools.build_phase1b_scaling_atlas import FOUNDATION  # noqa: E402
from tools.build_phase1d_structural_atlas import _row_for_matrix  # noqa: E402
import validation_suite as vs  # noqa: E402


CASE_N = 64
CASE_SEED = 1959
OPTIMIZER_SEED = 1987
WARMUP_LIMIT = 10
ANNEAL_LIMIT = 10
MAX_DATA = 4
INITIAL_TEMP = 100.0
COOLING_FACTOR = 0.9
CONTROL_EMBEDDING_DIM = 2


CSV_HEADERS = (
    "family",
    "target_dim",
    "n",
    "seed",
    "optimizer_seed",
    "embedding_dim",
    "mm_dim",
    "midpoint_dim",
    "abs_discrepancy_mm_midpoint",
    "chain2_count",
    "chain3_count",
    "chain3_abundance",
    "initial_energy",
    "warmup_energy",
    "final_energy",
    "truth_energy",
    "energy_gap",
    "interval_rmse",
    "optimizer_status",
    "failure_mode",
    "runtime_seconds",
    "optimizer_steps",
)


def _format_csv_value(value) -> str:
    if value is None:
        return "NA"
    if isinstance(value, str):
        return value
    return _format_field(value)


def _minkowski_case(d_spacetime: int) -> dict:
    matrix, points = vs.sprinkle_minkowski_diamond(
        n=CASE_N, seed=CASE_SEED, d_spacetime=d_spacetime
    )
    case = vs.SprinkleCase(
        d_spacetime=d_spacetime,
        n=CASE_N,
        seed=CASE_SEED,
        matrix=matrix,
        points=points,
    )
    pre = _row_for_matrix(
        family="minkowski",
        target_dim=d_spacetime,
        n=CASE_N,
        seed=CASE_SEED,
        matrix=matrix,
    )

    start = time.perf_counter()
    try:
        result = vs.run_recovery(
            case,
            optimizer_seed=OPTIMIZER_SEED,
            target_dim=d_spacetime - 1,
            warmup_limit=WARMUP_LIMIT,
            anneal_limit=ANNEAL_LIMIT,
            max_data=MAX_DATA,
            initial_temp=INITIAL_TEMP,
            cooling_factor=COOLING_FACTOR,
            backend="cpu",
        )
        runtime = time.perf_counter() - start
        final_energy = result.final_energy
        truth_energy = result.truth_energy
        energy_gap = final_energy - truth_energy
        row = {
            **pre,
            "optimizer_seed": OPTIMIZER_SEED,
            "embedding_dim": d_spacetime - 1,
            "initial_energy": result.initial_energy,
            "warmup_energy": result.warmup_energy,
            "final_energy": final_energy,
            "truth_energy": truth_energy,
            "energy_gap": energy_gap,
            "interval_rmse": result.interval_rmse,
            "optimizer_status": "completed",
            "failure_mode": "",
            "runtime_seconds": runtime,
            "optimizer_steps": MAX_DATA,
        }
    except Exception as exc:  # pragma: no cover - exercised only on failure
        runtime = time.perf_counter() - start
        row = {
            **pre,
            "optimizer_seed": OPTIMIZER_SEED,
            "embedding_dim": d_spacetime - 1,
            "initial_energy": None,
            "warmup_energy": None,
            "final_energy": None,
            "truth_energy": None,
            "energy_gap": None,
            "interval_rmse": None,
            "optimizer_status": "failed",
            "failure_mode": type(exc).__name__,
            "runtime_seconds": runtime,
            "optimizer_steps": 0,
        }
    return row


def _control_case(family: str, matrix) -> dict:
    pre = _row_for_matrix(
        family=family,
        target_dim="",
        n=CASE_N,
        seed=CASE_SEED,
        matrix=matrix,
    )
    start = time.perf_counter()
    try:
        with tempfile.TemporaryDirectory() as tmpdir, contextlib.redirect_stdout(io.StringIO()):
            sim = cones.ConesSimulator(
                z=matrix,
                dim=CONTROL_EMBEDDING_DIM,
                seed=OPTIMIZER_SEED,
                interactive=False,
                max_data=MAX_DATA,
                plot_path=None,
                warmup_limit=WARMUP_LIMIT,
                anneal_limit=ANNEAL_LIMIT,
                initial_temp=INITIAL_TEMP,
                cooling_factor=COOLING_FACTOR,
                backend="cpu",
            )
            sim.run(Path(tmpdir) / "out.txt")
        runtime = time.perf_counter() - start
        row = {
            **pre,
            "optimizer_seed": OPTIMIZER_SEED,
            "embedding_dim": CONTROL_EMBEDDING_DIM,
            "initial_energy": sim.initial_energy,
            "warmup_energy": sim.warmup_energy,
            "final_energy": sim.data[-1][1] if sim.data else sim.eave,
            "truth_energy": None,
            "energy_gap": None,
            "interval_rmse": None,
            "optimizer_status": "completed",
            "failure_mode": "",
            "runtime_seconds": runtime,
            "optimizer_steps": len(sim.data),
        }
    except Exception as exc:  # pragma: no cover - exercised only on failure
        runtime = time.perf_counter() - start
        row = {
            **pre,
            "optimizer_seed": OPTIMIZER_SEED,
            "embedding_dim": CONTROL_EMBEDDING_DIM,
            "initial_energy": None,
            "warmup_energy": None,
            "final_energy": None,
            "truth_energy": None,
            "energy_gap": None,
            "interval_rmse": None,
            "optimizer_status": "failed",
            "failure_mode": type(exc).__name__,
            "runtime_seconds": runtime,
            "optimizer_steps": 0,
        }
    return row


def build_rows() -> list[dict]:
    rows = [_minkowski_case(d) for d in (2, 3, 4)]
    rows.append(
        _control_case(
            "kleitman_rothschild",
            vs.generate_kleitman_rothschild(n=CASE_N, seed=CASE_SEED),
        )
    )
    rows.append(
        _control_case(
            "corona_poset",
            vs.generate_corona_poset(n=CASE_N, seed=CASE_SEED),
        )
    )
    return rows


def write_csv(rows: list[dict], path: Path) -> None:
    lines = [",".join(CSV_HEADERS)]
    for row in rows:
        lines.append(",".join(_format_csv_value(row[h]) for h in CSV_HEADERS))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _display(value) -> str:
    if value is None:
        return "NA"
    if isinstance(value, str):
        return value or "-"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float) and math.isfinite(value):
        return _format_field(value)
    return "NA"


def write_markdown(rows: list[dict], path: Path) -> None:
    lines = [
        "# Phase 2 Embedding Bridge",
        "",
        "A minimal bridge between Phase 1D order-theoretic diagnostics",
        "and the historical Bombelli-Sorkin annealing code. This is a",
        "fixed, small probe rather than an optimizer search.",
        "",
        "Protocol:",
        "",
        f"- cases: n={CASE_N}, case seed={CASE_SEED}; Minkowski d=2,3,4,",
        "  Kleitman-Rothschild, and suspended corona.",
        f"- optimizer seed: {OPTIMIZER_SEED}.",
        f"- schedule: warmup_limit={WARMUP_LIMIT}, anneal_limit={ANNEAL_LIMIT},",
        f"  max_data={MAX_DATA}, initial_temp={INITIAL_TEMP},",
        f"  cooling_factor={COOLING_FACTOR}.",
        "- controls have no ground-truth coordinates, so `truth_energy`,",
        "  `energy_gap`, and `interval_rmse` are recorded as NA.",
        f"- controls are embedded at spatial dim={CONTROL_EMBEDDING_DIM};",
        "  Minkowski cases use spatial dim=d_spacetime-1.",
        "",
        "| family | d | embed dim | MM | midpoint | \\|disc\\| | C3 abundance | final E | truth E | gap | RMSE | status |",
        "| --- | :---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in rows:
        lines.append(
            "| {family} | {d} | {edim} | {mm} | {mid} | {disc} | {c3} | "
            "{final} | {truth} | {gap} | {rmse} | {status} |".format(
                family=row["family"],
                d=_display(row["target_dim"]),
                edim=_display(row["embedding_dim"]),
                mm=_display(row["mm_dim"]),
                mid=_display(row["midpoint_dim"]),
                disc=_display(row["abs_discrepancy_mm_midpoint"]),
                c3=_display(row["chain3_abundance"]),
                final=_display(row["final_energy"]),
                truth=_display(row["truth_energy"]),
                gap=_display(row["energy_gap"]),
                rmse=_display(row["interval_rmse"]),
                status=row["optimizer_status"],
            )
        )

    lines += [
        "",
        "Interpretation:",
        "",
        "- This probe is designed to expose qualitative alignment or",
        "  tension between structural diagnostics and annealing outcomes.",
        "- A low final energy alone is not treated as proof of a faithful",
        "  embedding; for Minkowski cases the useful checks are also",
        "  truth energy, energy gap, and interval RMSE.",
        "- If a Minkowski case has good structural diagnostics but a large",
        "  gap or RMSE, that points to optimizer/schedule failure rather",
        "  than immediate non-manifoldlikeness.",
        "- If a control reaches a deceptively low final energy, that is a",
        "  warning that the energy can reward non-geometric artifacts.",
        "",
        "Regenerate via `make regen-phase2`. Source tool:",
        "`tools/build_phase2_embedding_bridge.py`.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    rows = build_rows()
    FOUNDATION.mkdir(parents=True, exist_ok=True)
    write_csv(rows, FOUNDATION / "phase2_embedding_bridge.csv")
    write_markdown(rows, FOUNDATION / "phase2_embedding_bridge.md")
    print(f"Wrote {len(rows)} Phase 2 bridge rows to {FOUNDATION}")


if __name__ == "__main__":
    main()
