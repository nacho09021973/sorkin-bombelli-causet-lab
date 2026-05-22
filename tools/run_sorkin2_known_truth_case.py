#!/usr/bin/env python3
"""Run one enabled SORKIN-2 known-truth case with verifier artifacts.

Only the two trivial sanity checks are enabled here. This harness is
intentionally narrow: it runs one permanent input, verifies the recovered order with
``validation_suite.verify_recovery``, and writes the required order-matrix
figures into a traceable run directory.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import io
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

import cones
import validation_suite as vs


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_ROOT = ROOT / "results" / "sorkin2_known_truth"

PRIMARY_CRITERION_WARNING = (
    "exact_match from validation_suite.compare_causal_orders is the primary "
    "criterion for exact recovery; low final energy alone is not sufficient."
)


@dataclass(frozen=True)
class KnownTruthCase:
    case_id: str
    input_file: Path
    dim: int
    annealer_mode: str
    warmup_limit: int = 100
    anneal_limit: int = 100
    max_data: int = 35
    initial_temp: float = 100.0
    cooling_factor: float = 0.9
    acceptance_scale: float = 4.0
    backend: str = "cpu"


ENABLED_CASES = {
    "chain_4_d2": KnownTruthCase(
        case_id="chain_4_d2",
        input_file=ROOT / "benchmarks" / "known_truth" / "chain_4_d2.in",
        dim=2,
        annealer_mode="historical/default",
    ),
    "antichain_4_d2": KnownTruthCase(
        case_id="antichain_4_d2",
        input_file=ROOT / "benchmarks" / "known_truth" / "antichain_4_d2.in",
        dim=2,
        annealer_mode="historical/default",
    ),
    "chain_12_d2": KnownTruthCase(
        case_id="chain_12_d2",
        input_file=ROOT / "benchmarks" / "known_truth" / "n12_topology_panel" / "chain_12_d2.in",
        dim=2,
        annealer_mode="historical/default",
    ),
    "antichain_12_d2": KnownTruthCase(
        case_id="antichain_12_d2",
        input_file=ROOT / "benchmarks" / "known_truth" / "n12_topology_panel" / "antichain_12_d2.in",
        dim=2,
        annealer_mode="historical/default",
    ),
    "layered_4_4_4_d2": KnownTruthCase(
        case_id="layered_4_4_4_d2",
        input_file=ROOT / "benchmarks" / "known_truth" / "n12_topology_panel" / "layered_4_4_4_d2.in",
        dim=2,
        annealer_mode="historical/default",
    ),
    "layered_4_4_4_d2_T100_g08": KnownTruthCase(
        case_id="layered_4_4_4_d2_T100_g08",
        input_file=ROOT / "benchmarks" / "known_truth" / "n12_topology_panel" / "layered_4_4_4_d2.in",
        dim=2,
        annealer_mode="mechanism/T100_g08",
        initial_temp=100.0,
        cooling_factor=0.8,
    ),
    "minkowski_6_s1959_d2": KnownTruthCase(
        case_id="minkowski_6_s1959_d2",
        input_file=ROOT / "benchmarks" / "tesis_like_6.in",
        dim=2,
        annealer_mode="historical/default",
    ),
    "minkowski_12_s1962_d2_hist": KnownTruthCase(
        case_id="minkowski_12_s1962_d2_hist",
        input_file=ROOT / "benchmarks" / "tesis_like_12.in",
        dim=2,
        annealer_mode="historical/default",
    ),
    "minkowski_12_s1962_d2_tuned": KnownTruthCase(
        case_id="minkowski_12_s1962_d2_tuned",
        input_file=ROOT / "benchmarks" / "tesis_like_12.in",
        dim=2,
        annealer_mode="tuned/non-historical",
        initial_temp=180.0,
        cooling_factor=0.8,
    ),
    "minkowski_12_s1962_d2_T180_g09": KnownTruthCase(
        case_id="minkowski_12_s1962_d2_T180_g09",
        input_file=ROOT / "benchmarks" / "tesis_like_12.in",
        dim=2,
        annealer_mode="mechanism/T180_g09",
        initial_temp=180.0,
        cooling_factor=0.9,
    ),
    "minkowski_12_s1962_d2_T100_g08": KnownTruthCase(
        case_id="minkowski_12_s1962_d2_T100_g08",
        input_file=ROOT / "benchmarks" / "tesis_like_12.in",
        dim=2,
        annealer_mode="mechanism/T100_g08",
        initial_temp=100.0,
        cooling_factor=0.8,
    ),
}


def make_run_id(case_id: str, seed: int, now: datetime | None = None) -> str:
    """Build a stable, sortable run id."""

    if now is None:
        now = datetime.now(timezone.utc)
    timestamp = now.astimezone(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return f"{timestamp}_{case_id}_seed{seed}"


def build_run_dir(output_root: Path, case_id: str, run_id: str) -> Path:
    """Return the required SORKIN-2 known-truth run directory."""

    return output_root / case_id / run_id


def get_enabled_case(case_id: str) -> KnownTruthCase:
    """Return an enabled case or fail clearly for not-yet-enabled cases."""

    try:
        return ENABLED_CASES[case_id]
    except KeyError as exc:
        enabled = ", ".join(sorted(ENABLED_CASES))
        raise ValueError(
            f"unsupported SORKIN-2 known-truth case_id {case_id!r}; "
            f"enabled case_ids: {enabled}"
        ) from exc


def run_annealer(
    case: KnownTruthCase,
    *,
    seed: int,
    output_path: Path,
    block_callback: Any | None = None,
) -> cones.ConesSimulator:
    """Run the historical Bombelli annealer for one enabled known-truth case."""

    z = cones.parse_cones_input(case.input_file)
    sim = cones.ConesSimulator(
        z=z,
        dim=case.dim,
        seed=seed,
        interactive=False,
        max_data=case.max_data,
        plot_path=None,
        backend=case.backend,
        warmup_limit=case.warmup_limit,
        anneal_limit=case.anneal_limit,
        initial_temp=case.initial_temp,
        cooling_factor=case.cooling_factor,
        acceptance_scale=case.acceptance_scale,
        block_callback=block_callback,
    )
    with contextlib.redirect_stdout(io.StringIO()):
        sim.run(output_path)
    return sim


def _json_pair_list(pairs: Sequence[tuple[int, int]]) -> list[list[int]]:
    return [[i, j] for i, j in pairs]


def _target_pairs(z: Sequence[Sequence[bool]]) -> list[tuple[int, int]]:
    return [
        (i, j)
        for i in range(len(z))
        for j in range(i + 1, len(z))
        if bool(z[i][j])
    ]


def _induced_pairs(z_induced: Sequence[Sequence[bool]]) -> list[tuple[int, int]]:
    return [
        (i, j)
        for i in range(len(z_induced))
        for j in range(i + 1, len(z_induced))
        if bool(z_induced[i][j])
    ]


def _coordinates_record(sim: Any) -> list[dict[str, Any]]:
    return [
        {
            "i": i,
            "R": sim.rold[i],
            "X": list(sim.xold[i]),
        }
        for i in range(sim.n)
    ]


def _final_energy(sim: Any) -> float | None:
    if getattr(sim, "data", None):
        return float(sim.data[-1][1])
    if hasattr(sim, "eave"):
        return float(sim.eave)
    return None


_TRACE_FIELDS = [
    "block",
    "temperature",
    "energy",
    "induced_relations",
    "missing_relations_count",
    "extra_relations_count",
    "exact_match",
]

_TRACE_NOTE = (
    "All columns are populated at every block when run via run_case() (block_callback "
    "captures rold/xold per block). Legacy/synthetic runs without a callback have "
    "causal-order columns populated at the final block only."
)


def write_trace_artifacts(
    *,
    sim: Any,
    comparison: vs.OrderComparison,
    run_dir: Path,
    block_records: list[dict[str, Any]] | None = None,
) -> dict[str, str]:
    """Write trace.csv, trace_energy.png, and trace_relations.png for one run.

    When ``block_records`` is provided (full per-block data from a block callback),
    all causal-order columns are populated at every row. When omitted, causal-order
    columns are populated at the final block only, using the post-run comparison.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    target = comparison.total_relations_target

    # --- Build CSV rows ---
    if block_records is not None:
        rows: list[dict[str, Any]] = [
            {
                "block": rec["block"],
                "temperature": rec["temperature"],
                "energy": rec["energy"],
                "induced_relations": rec["induced_relations"],
                "missing_relations_count": rec["missing_relations_count"],
                "extra_relations_count": rec["extra_relations_count"],
                "exact_match": "true" if rec["exact_match"] else "false",
            }
            for rec in block_records
        ]
        final = block_records[-1] if block_records else {}
        final_induced = int(final.get("induced_relations", 0))
        final_missing = int(final.get("missing_relations_count", 0))
        final_extra = int(final.get("extra_relations_count", 0))
        final_exact = bool(final.get("exact_match", False))
    else:
        block_data = list(getattr(sim, "data", []))
        n_fallback = len(block_data)
        rows = []
        for idx, (temp, eave) in enumerate(block_data):
            is_final = idx == n_fallback - 1
            rows.append({
                "block": idx + 1,
                "temperature": temp,
                "energy": eave,
                "induced_relations": comparison.total_relations_induced if is_final else "",
                "missing_relations_count": len(comparison.missing_relations) if is_final else "",
                "extra_relations_count": len(comparison.extra_relations) if is_final else "",
                "exact_match": ("true" if comparison.exact_match else "false") if is_final else "",
            })
        final_induced = comparison.total_relations_induced
        final_missing = len(comparison.missing_relations)
        final_extra = len(comparison.extra_relations)
        final_exact = comparison.exact_match

    n_blocks = len(rows)

    # --- trace.csv ---
    trace_csv = run_dir / "trace.csv"
    with trace_csv.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=_TRACE_FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    # --- trace_pairs.jsonl and trace_coordinates.jsonl ---
    trace_pairs_jsonl = run_dir / "trace_pairs.jsonl"
    trace_coordinates_jsonl = run_dir / "trace_coordinates.jsonl"
    with trace_pairs_jsonl.open("w", encoding="utf-8") as pairs_fh, trace_coordinates_jsonl.open(
        "w", encoding="utf-8"
    ) as coords_fh:
        if block_records is not None:
            for rec in block_records:
                pairs_item = {
                    "block": rec["block"],
                    "temperature": rec["temperature"],
                    "energy": rec["energy"],
                    "induced_relations": rec["induced_relations"],
                    "correct_relations": rec.get(
                        "correct_relations",
                        rec["induced_relations"] - rec["extra_relations_count"],
                    ),
                    "missing_relations_count": rec["missing_relations_count"],
                    "extra_relations_count": rec["extra_relations_count"],
                    "exact_match": rec["exact_match"],
                    "correct_pairs": rec.get("correct_pairs", []),
                    "missing_pairs": rec.get("missing_pairs", []),
                    "extra_pairs": rec.get("extra_pairs", []),
                }
                coords_item = {
                    "block": rec["block"],
                    "coordinates": rec.get("coordinates", []),
                }
                pairs_fh.write(json.dumps(pairs_item, sort_keys=True) + "\n")
                coords_fh.write(json.dumps(coords_item, sort_keys=True) + "\n")
        else:
            # Fallback for synthetic/legacy callers that have no block snapshots:
            # keep one JSONL line per CSV row, with pair lists populated only
            # where the existing CSV summary has final-state order data.
            correct_pairs = [
                pair
                for pair in _target_pairs(sim.z)
                if pair not in set(comparison.missing_relations)
            ]
            for idx, row in enumerate(rows):
                has_order_data = row["induced_relations"] != ""
                pairs_item = {
                    "block": row["block"],
                    "temperature": row["temperature"],
                    "energy": row["energy"],
                    "induced_relations": row["induced_relations"],
                    "correct_relations": len(correct_pairs) if has_order_data else "",
                    "missing_relations_count": row["missing_relations_count"],
                    "extra_relations_count": row["extra_relations_count"],
                    "exact_match": row["exact_match"],
                    "correct_pairs": _json_pair_list(correct_pairs) if has_order_data else [],
                    "missing_pairs": _json_pair_list(comparison.missing_relations) if has_order_data else [],
                    "extra_pairs": _json_pair_list(comparison.extra_relations) if has_order_data else [],
                }
                coords_item = {
                    "block": row["block"],
                    "coordinates": _coordinates_record(sim) if idx == len(rows) - 1 else [],
                }
                pairs_fh.write(json.dumps(pairs_item, sort_keys=True) + "\n")
                coords_fh.write(json.dumps(coords_item, sort_keys=True) + "\n")

    blocks = [r["block"] for r in rows]
    energies = [r["energy"] for r in rows]

    # --- trace_energy.png ---
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(blocks, energies, "o-", color="#1f77b4", lw=1.5, ms=4)
    ax.set_xlabel("Cooling block")
    ax.set_ylabel("E_ave")
    ax.set_title("Energy per cooling block")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    energy_png = run_dir / "trace_energy.png"
    fig.savefig(energy_png, dpi=120)
    plt.close(fig)

    # --- trace_relations.png ---
    # Top panel: energy trace. Bottom panel: final-block relation bar.
    fig, (ax_e, ax_r) = plt.subplots(
        2, 1, figsize=(7, 6), gridspec_kw={"height_ratios": [3, 1]}
    )

    ax_e.plot(blocks, energies, "o-", color="#1f77b4", lw=1.5, ms=4)
    ax_e.set_ylabel("E_ave")
    ax_e.set_title("Energy trace — final-block relation counts below")
    ax_e.grid(True, alpha=0.3)

    ax_r.bar(
        ["induced", "missing", "extra"],
        [final_induced, final_missing, final_extra],
        color=["#2166ac", "#d73027", "#ff7f00"],
        width=0.5,
    )
    ax_r.axhline(target, color="#333", ls="--", lw=1, label=f"target={target}")
    ax_r.set_ylabel("count")
    ax_r.set_title(
        f"Block {n_blocks}: induced={final_induced}/{target}  "
        f"exact={'yes' if final_exact else 'no'}"
    )
    ax_r.legend(fontsize=8)
    ax_r.set_ylim(0, max(target + 2, final_induced + 2, final_missing + 2, final_extra + 2, 3))
    fig.text(0.5, 0.005, _TRACE_NOTE, ha="center", fontsize=6.5, color="#888",
             wrap=True)
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    relations_png = run_dir / "trace_relations.png"
    fig.savefig(relations_png, dpi=120)
    plt.close(fig)

    return {
        "trace_csv": trace_csv.name,
        "trace_pairs": trace_pairs_jsonl.name,
        "trace_coordinates": trace_coordinates_jsonl.name,
        "trace_energy": energy_png.name,
        "trace_relations": relations_png.name,
        "note": _TRACE_NOTE,
    }


def write_result_artifacts(
    *,
    case: KnownTruthCase,
    seed: int,
    run_id: str,
    run_dir: Path,
    sim: Any,
    block_records: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Write result.json and order-matrix figures for a completed run."""

    comparison = vs.verify_recovery(sim)
    recovered_coords = [(sim.rold[i], *sim.xold[i]) for i in range(sim.n)]
    induced = vs.induced_order_from_coords(recovered_coords)
    figure_paths = vs.write_order_matrix_plots(sim.z, induced, run_dir)

    generated_figures = {
        "target_order_matrix": figure_paths.target_order_matrix.name,
        "induced_order_matrix": figure_paths.induced_order_matrix.name,
        "order_difference_matrix": figure_paths.order_difference_matrix.name,
    }
    result: dict[str, Any] = {
        "case_id": case.case_id,
        "run_id": run_id,
        "input_file": str(case.input_file.relative_to(ROOT)),
        "annealer_mode": case.annealer_mode,
        "seed": seed,
        "n": comparison.n,
        "dim": case.dim,
        "schedule": {
            "warmup_limit": case.warmup_limit,
            "anneal_limit": case.anneal_limit,
            "max_data": case.max_data,
            "initial_temp": case.initial_temp,
            "cooling_factor": case.cooling_factor,
            "acceptance_scale": case.acceptance_scale,
            "backend": case.backend,
        },
        "final_energy": _final_energy(sim),
        "exact_match": comparison.exact_match,
        "total_relations_target": comparison.total_relations_target,
        "total_relations_induced": comparison.total_relations_induced,
        "missing_relations": _json_pair_list(comparison.missing_relations),
        "extra_relations": _json_pair_list(comparison.extra_relations),
        "generated_figures": generated_figures,
        "code_paths": [
            "tools/run_sorkin2_known_truth_case.py",
            "cones.py",
            "validation_suite.py",
        ],
        "verifier": "validation_suite.verify_recovery",
        "primary_recovery_criterion_warning": PRIMARY_CRITERION_WARNING,
    }

    trace_artifacts = write_trace_artifacts(
        sim=sim, comparison=comparison, run_dir=run_dir, block_records=block_records
    )
    result["trace_artifacts"] = trace_artifacts

    (run_dir / "result.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return result


def run_case(
    *,
    case_id: str,
    seed: int,
    output_root: Path = DEFAULT_OUTPUT_ROOT,
    run_id: str | None = None,
) -> Path:
    """Run one enabled known-truth case and return its run directory."""

    case = get_enabled_case(case_id)
    if run_id is None:
        run_id = make_run_id(case_id, seed)
    run_dir = build_run_dir(output_root, case_id, run_id)
    if run_dir.exists() and any(run_dir.iterdir()):
        raise FileExistsError(f"run directory already exists and is not empty: {run_dir}")
    run_dir.mkdir(parents=True, exist_ok=True)

    block_records: list[dict[str, Any]] = []

    def _block_callback(sim: Any, block_idx: int, temp: float, eave: float) -> None:
        coords = [(sim.rold[i], *sim.xold[i]) for i in range(sim.n)]
        induced = vs.induced_order_from_coords(coords)
        comp = vs.compare_causal_orders(sim.z, induced)
        target_pairs = set(_target_pairs(sim.z))
        induced_pairs = set(_induced_pairs(induced))
        correct_pairs = sorted(target_pairs & induced_pairs)
        block_records.append({
            "block": block_idx,
            "temperature": temp,
            "energy": eave,
            "induced_relations": comp.total_relations_induced,
            "correct_relations": len(correct_pairs),
            "missing_relations_count": len(comp.missing_relations),
            "extra_relations_count": len(comp.extra_relations),
            "exact_match": comp.exact_match,
            "correct_pairs": _json_pair_list(correct_pairs),
            "missing_pairs": _json_pair_list(comp.missing_relations),
            "extra_pairs": _json_pair_list(comp.extra_relations),
            "coordinates": _coordinates_record(sim),
        })

    annealer_output = run_dir / "annealer_output.txt"
    sim = run_annealer(case, seed=seed, output_path=annealer_output, block_callback=_block_callback)
    write_result_artifacts(
        case=case,
        seed=seed,
        run_id=run_id,
        run_dir=run_dir,
        sim=sim,
        block_records=block_records,
    )
    return run_dir


def build_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run one enabled SORKIN-2 known-truth case."
    )
    parser.add_argument("case_id", help="enabled case id")
    parser.add_argument("--seed", type=int, default=1959, help="optimizer seed")
    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help="root for results/sorkin2_known_truth style outputs",
    )
    parser.add_argument("--run-id", help="optional explicit run id")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_argparser().parse_args(argv)
    run_dir = run_case(
        case_id=args.case_id,
        seed=args.seed,
        output_root=args.output_root,
        run_id=args.run_id,
    )
    print(run_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
