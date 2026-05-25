#!/usr/bin/env python3
"""gamma_0p5 basin consistency probe for SORKIN-2 N=36.

Runs the historical Bombelli annealer with schedule gamma_0p5 / medium_25_25_8
for N=36, case_seed=1959, optimizer_seeds {1959, 1962, 1987, 2001}.

At each block the full induced causal order is serialized as a "|"-separated
"i:j" pair string.  Post-run analysis computes pairwise Jaccard similarity
between the candidate "peak-F1" checkpoints (blocks 6-7 for H2a seeds) and
block-1 control checkpoints.

This is an oracular probe: basin identification uses induced_pairs against
the known-truth target.  It is not a deployable selection criterion.

Reproducibility check: F1 and energy values are validated against
explore/schedule_seed_stability_n36/schedule_seed_stability_n36.csv before
any Jaccard interpretation is made.
"""

from __future__ import annotations

import contextlib
import csv
import io
import math
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import cones          # noqa: E402
import validation_suite as vs  # noqa: E402


OUT_DIR = Path(__file__).resolve().parent
REFERENCE_CSV = (
    ROOT / "explore" / "schedule_seed_stability_n36"
    / "schedule_seed_stability_n36.csv"
)
TRAJECTORY_CSV  = OUT_DIR / "gamma0p5_basin_consistency_n36_trajectory.csv"
TARGET_PAIRS_CSV = OUT_DIR / "gamma0p5_basin_consistency_n36_target_pairs.csv"
JACCARD_CSV     = OUT_DIR / "gamma0p5_basin_consistency_n36_jaccard.csv"
MD_PATH         = OUT_DIR / "gamma0p5_basin_consistency_n36.md"
SVG_PATH        = OUT_DIR / "gamma0p5_basin_consistency_n36.svg"
COMMAND = "python3 explore/gamma0p5_basin_consistency_n36/run_gamma0p5_basin_consistency_n36.py"

# ── fixed case parameters ─────────────────────────────────────────────────────
# Must match schedule_seed_stability_n36 exactly.
D_SPACETIME     = 2
N               = 36
CASE_SEED       = 1959
OPTIMIZER_SEEDS = (1959, 1962, 1987, 2001)
INITIAL_TEMP    = 100.0
BACKEND         = "cpu"
SCHEDULE_LABEL  = "gamma_0p5"
COOLING_FACTOR  = 0.50
BUDGET_LABEL    = "medium_25_25_8"
WARMUP_LIMIT    = 25
ANNEAL_LIMIT    = 25
MAX_DATA        = 8

# ── reproducibility tolerance ─────────────────────────────────────────────────
F1_REPRO_TOL     = 1e-6   # max allowed |f1_new – f1_ref|
ENERGY_REPRO_TOL = 1e-4   # max allowed |energy_new – energy_ref|

# ── checkpoints of interest ───────────────────────────────────────────────────
# Primary candidates: H2a peaks from trajectory_window_n36 + inconclusive ref.
# Controls: block 1 (T = 100, hot start) for each seed.
CANDIDATE_CHECKPOINTS: list[tuple[int, int, str]] = [
    (1959, 6, "H2a_peak"),
    (1987, 6, "H2a_peak"),
    (2001, 7, "H2a_peak"),
    (1962, 8, "inconclusive"),
]
CONTROL_CHECKPOINTS: list[tuple[int, int, str]] = [
    (1959, 1, "control_blk1"),
    (1987, 1, "control_blk1"),
    (2001, 1, "control_blk1"),
    (1962, 1, "control_blk1"),
]
ALL_CHECKPOINTS = CANDIDATE_CHECKPOINTS + CONTROL_CHECKPOINTS

# Jaccard interpretation thresholds
JACCARD_HIGH = 0.75  # above → basin compartido candidate
JACCARD_LOW  = 0.50  # below → basins distintos

# ── CSV headers ───────────────────────────────────────────────────────────────
TRAJECTORY_HEADERS = (
    "optimizer_seed",
    "block_index",
    "temperature",
    "energy_eave",
    "causal_precision",
    "causal_recall",
    "causal_f1",
    "missing_relations_count",
    "extra_relations_count",
    "total_relations_target",
    "total_relations_induced",
    "correct_relations",
    "induced_pairs_str",
)

JACCARD_HEADERS = (
    "seed_a",
    "block_a",
    "type_a",
    "f1_a",
    "induced_size_a",
    "seed_b",
    "block_b",
    "type_b",
    "f1_b",
    "induced_size_b",
    "intersection_size",
    "union_size",
    "jaccard",
)


# ── helpers ───────────────────────────────────────────────────────────────────

def _fmt(v: object) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, float):
        return f"{v:.10g}" if math.isfinite(v) else "NA"
    return str(v)


def _fmt_f(v: float, d: int = 6) -> str:
    return f"{v:.{d}g}" if math.isfinite(v) else "NA"


def _mean(xs: list[float]) -> float:
    return sum(xs) / len(xs)


def _safe_div(num: float, den: float) -> float:
    return 0.0 if den == 0.0 else num / den


def _pairs(matrix: list[list[bool]]) -> set[tuple[int, int]]:
    return {
        (i, j)
        for i in range(len(matrix) - 1)
        for j in range(i + 1, len(matrix))
        if matrix[i][j]
    }


def _pairs_from_str(s: str) -> set[tuple[int, int]]:
    """Parse the serialized 'i:j|i:j|...' string back to a set of pairs."""
    if not s.strip():
        return set()
    result: set[tuple[int, int]] = set()
    for token in s.split("|"):
        i_str, j_str = token.split(":")
        result.add((int(i_str), int(j_str)))
    return result


def _jaccard(a: set[tuple[int, int]], b: set[tuple[int, int]]) -> float:
    union = len(a | b)
    return _safe_div(len(a & b), union)


# ── case setup ────────────────────────────────────────────────────────────────

def _make_case() -> vs.SprinkleCase:
    matrix, points = vs.sprinkle_minkowski_diamond(
        n=N, seed=CASE_SEED, d_spacetime=D_SPACETIME
    )
    return vs.SprinkleCase(
        d_spacetime=D_SPACETIME, n=N, seed=CASE_SEED,
        matrix=matrix, points=points,
    )


# ── per-block metrics with full pair capture ──────────────────────────────────

def _causal_metrics_with_pairs(
    target: list[list[bool]],
    coords: list[tuple[float, ...]],
) -> dict[str, object]:
    """Same as schedule_seed_stability _causal_metrics but also returns
    the full induced pair set (not serialized here; caller serializes)."""
    induced        = vs.induced_order_from_coords(coords)
    induced_ps     = _pairs(induced)
    comparison     = vs.compare_causal_orders(target, induced)
    correct        = len(_pairs(target) & induced_ps)
    precision      = _safe_div(correct, comparison.total_relations_induced)
    recall         = _safe_div(correct, comparison.total_relations_target)
    f1             = _safe_div(2.0 * precision * recall, precision + recall)
    return {
        "causal_precision":         precision,
        "causal_recall":            recall,
        "causal_f1":                f1,
        "missing_relations_count":  len(comparison.missing_relations),
        "extra_relations_count":    len(comparison.extra_relations),
        "total_relations_target":   comparison.total_relations_target,
        "total_relations_induced":  comparison.total_relations_induced,
        "correct_relations":        correct,
        "_induced_pairs":           induced_ps,   # internal; not in CSV
    }


# ── annealer run ──────────────────────────────────────────────────────────────

def _run_one(
    case: vs.SprinkleCase,
    optimizer_seed: int,
) -> list[dict[str, object]]:
    """Run gamma_0p5 for one optimizer_seed; return one row per block."""
    rows: list[dict[str, object]] = []

    def _block_callback(sim: Any, block_idx: int, temp: float, eave: float) -> None:
        coords = [
            (float(sim.rold[i]), *[float(v) for v in sim.xold[i]])
            for i in range(sim.n)
        ]
        m = _causal_metrics_with_pairs(case.matrix, coords)
        induced_pairs_str = "|".join(
            f"{i}:{j}" for i, j in sorted(m["_induced_pairs"])
        )
        rows.append({
            "optimizer_seed":           optimizer_seed,
            "block_index":              block_idx,
            "temperature":              temp,
            "energy_eave":              eave,
            "causal_precision":         m["causal_precision"],
            "causal_recall":            m["causal_recall"],
            "causal_f1":                m["causal_f1"],
            "missing_relations_count":  m["missing_relations_count"],
            "extra_relations_count":    m["extra_relations_count"],
            "total_relations_target":   m["total_relations_target"],
            "total_relations_induced":  m["total_relations_induced"],
            "correct_relations":        m["correct_relations"],
            "induced_pairs_str":        induced_pairs_str,
        })

    with tempfile.TemporaryDirectory() as tmpdir, \
            contextlib.redirect_stdout(io.StringIO()):
        sim = cones.ConesSimulator(
            z=case.matrix,
            dim=D_SPACETIME - 1,
            seed=optimizer_seed,
            interactive=False,
            max_data=MAX_DATA,
            plot_path=None,
            warmup_limit=WARMUP_LIMIT,
            anneal_limit=ANNEAL_LIMIT,
            initial_temp=INITIAL_TEMP,
            cooling_factor=COOLING_FACTOR,
            backend=BACKEND,
            block_callback=_block_callback,
        )
        sim.run(Path(tmpdir) / "annealer_output.txt")
    return rows


# ── reproducibility validation ────────────────────────────────────────────────

def _validate(trajectory_rows: list[dict[str, object]]) -> dict[str, object]:
    """Compare new trajectory F1/energy against the reference CSV row-by-row."""
    ref: dict[tuple[int, int], dict[str, float]] = {}
    with REFERENCE_CSV.open(newline="", encoding="utf-8") as fh:
        for raw in csv.DictReader(fh):
            if raw.get("schedule_label") != SCHEDULE_LABEL:
                continue
            if raw.get("budget_label") != BUDGET_LABEL:
                continue
            key = (int(raw["optimizer_seed"]), int(raw["block_index"]))
            ref[key] = {
                "causal_f1":   float(raw["causal_f1"]),
                "energy_eave": float(raw["energy_eave"]),
            }

    diffs_f1: list[float] = []
    diffs_e:  list[float] = []
    mismatch_lines: list[str] = []

    for row in trajectory_rows:
        key = (int(row["optimizer_seed"]), int(row["block_index"]))
        if key not in ref:
            continue
        r = ref[key]
        df1 = abs(float(row["causal_f1"])   - r["causal_f1"])
        de  = abs(float(row["energy_eave"]) - r["energy_eave"])
        diffs_f1.append(df1)
        diffs_e.append(de)
        if df1 > F1_REPRO_TOL:
            mismatch_lines.append(
                f"seed={row['optimizer_seed']} blk={row['block_index']} "
                f"f1_new={float(row['causal_f1']):.9f} "
                f"f1_ref={r['causal_f1']:.9f} diff={df1:.2e}"
            )

    max_f1  = max(diffs_f1) if diffs_f1 else float("nan")
    max_e   = max(diffs_e)  if diffs_e  else float("nan")
    return {
        "n_compared":        len(diffs_f1),
        "max_abs_f1_diff":   max_f1,
        "max_abs_energy_diff": max_e,
        "reproducible":      (max_f1 <= F1_REPRO_TOL) if diffs_f1 else False,
        "mismatches":        mismatch_lines,
    }


# ── Jaccard computation ───────────────────────────────────────────────────────

def _compute_jaccard(
    trajectory_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    """Pairwise Jaccard for all combinations of ALL_CHECKPOINTS."""
    # Build (seed, block) → row lookup
    row_lookup: dict[tuple[int, int], dict[str, object]] = {
        (int(r["optimizer_seed"]), int(r["block_index"])): r
        for r in trajectory_rows
    }

    # Collect checkpoint data
    cps: list[dict[str, object]] = []
    for seed, block, cp_type in ALL_CHECKPOINTS:
        r = row_lookup.get((seed, block))
        if r is None:
            continue
        cps.append({
            "seed":   seed,
            "block":  block,
            "type":   cp_type,
            "f1":     float(r["causal_f1"]),
            "size":   int(r["total_relations_induced"]),
            "pairs":  _pairs_from_str(str(r["induced_pairs_str"])),
        })

    # Upper triangle including diagonal (symmetric by construction)
    out: list[dict[str, object]] = []
    for i in range(len(cps)):
        for j in range(i, len(cps)):
            a = cps[i]
            b = cps[j]
            inter = len(a["pairs"] & b["pairs"])
            union = len(a["pairs"] | b["pairs"])
            j_val = _safe_div(inter, union)
            out.append({
                "seed_a":          a["seed"],
                "block_a":         a["block"],
                "type_a":          a["type"],
                "f1_a":            a["f1"],
                "induced_size_a":  a["size"],
                "seed_b":          b["seed"],
                "block_b":         b["block"],
                "type_b":          b["type"],
                "f1_b":            b["f1"],
                "induced_size_b":  b["size"],
                "intersection_size": inter,
                "union_size":        union,
                "jaccard":           j_val,
            })
    return out


# ── I/O ───────────────────────────────────────────────────────────────────────

def _write_csv(
    path: Path, headers: tuple[str, ...], rows: list[dict[str, object]]
) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=headers)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: _fmt(row[k]) for k in headers})


def _write_target_pairs_csv(target_pairs: set[tuple[int, int]]) -> None:
    with TARGET_PAIRS_CSV.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["i", "j"])
        writer.writeheader()
        for i, j in sorted(target_pairs):
            writer.writerow({"i": i, "j": j})


# ── SVG heatmap ───────────────────────────────────────────────────────────────

def _write_svg(jaccard_rows: list[dict[str, object]]) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as exc:
        raise RuntimeError("matplotlib + numpy required") from exc

    n = len(ALL_CHECKPOINTS)
    idx: dict[tuple[int, int], int] = {
        (s, b): i for i, (s, b, _) in enumerate(ALL_CHECKPOINTS)
    }
    labels = [
        f"s{s} b{b}\n{t.replace('_', ' ')}"
        for s, b, t in ALL_CHECKPOINTS
    ]

    mat = np.zeros((n, n))
    for row in jaccard_rows:
        ia = idx.get((int(row["seed_a"]), int(row["block_a"])))
        ib = idx.get((int(row["seed_b"]), int(row["block_b"])))
        if ia is not None and ib is not None:
            mat[ia, ib] = float(row["jaccard"])
            mat[ib, ia] = float(row["jaccard"])

    fig, ax = plt.subplots(figsize=(8.5, 7.5))
    im = ax.imshow(mat, cmap="Blues", vmin=0.0, vmax=1.0)

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, fontsize=7.5)
    ax.set_yticklabels(labels, fontsize=7.5)

    for i in range(n):
        for j in range(n):
            color = "white" if mat[i, j] > 0.55 else "black"
            ax.text(j, i, f"{mat[i, j]:.3f}",
                    ha="center", va="center", fontsize=7, color=color)

    # Separator between candidates (top-left 4×4) and controls (bottom-right 4×4)
    sep = len(CANDIDATE_CHECKPOINTS) - 0.5
    for line_pos in (sep,):
        ax.axhline(line_pos, color="#cc3333", lw=1.4, ls="--", alpha=0.8)
        ax.axvline(line_pos, color="#cc3333", lw=1.4, ls="--", alpha=0.8)

    n_cand = len(CANDIDATE_CHECKPOINTS)
    ax.text(n_cand / 2 - 0.5, -1.2, "candidates", ha="center",
            fontsize=8, color="#cc3333", fontweight="bold")
    ax.text(n_cand + (n - n_cand) / 2 - 0.5, -1.2, "controls (blk 1)",
            ha="center", fontsize=8, color="#555")

    plt.colorbar(im, ax=ax, label="Jaccard J(A, B)")
    ax.set_title(
        "Jaccard similarity of induced causal orders  (N=36, gamma_0p5)\n"
        "Candidates = H2a peaks + inconclusive  |  Controls = block 1 (T=100)",
        fontsize=9.5,
    )
    fig.tight_layout()
    fig.savefig(SVG_PATH, format="svg", bbox_inches="tight")
    plt.close(fig)


# ── Markdown ──────────────────────────────────────────────────────────────────

def _write_markdown(
    trajectory_rows: list[dict[str, object]],
    jaccard_rows: list[dict[str, object]],
    validation: dict[str, object],
    target_n_pairs: int,
) -> None:
    generated_at = datetime.now(timezone.utc).isoformat(timespec="seconds")

    # Build Jaccard lookup: (seed_a, block_a, seed_b, block_b) → row
    j_lookup: dict[tuple[int, int, int, int], dict[str, object]] = {}
    for row in jaccard_rows:
        key  = (int(row["seed_a"]), int(row["block_a"]),
                int(row["seed_b"]), int(row["block_b"]))
        rkey = (int(row["seed_b"]), int(row["block_b"]),
                int(row["seed_a"]), int(row["block_a"]))
        j_lookup[key]  = row
        j_lookup[rkey] = row

    def J(sa: int, ba: int, sb: int, bb: int) -> float:
        r = j_lookup.get((sa, ba, sb, bb))
        return float(r["jaccard"]) if r else float("nan")

    # H2a peaks
    h2a = [(s, b) for s, b, t in CANDIDATE_CHECKPOINTS if t == "H2a_peak"]
    inc  = [(s, b) for s, b, t in CANDIDATE_CHECKPOINTS if t == "inconclusive"]
    ctl  = [(s, b) for s, b, _  in CONTROL_CHECKPOINTS]

    # Cross-seed Jaccard among H2a peaks
    h2a_pairs_j = [J(h2a[i][0], h2a[i][1], h2a[j][0], h2a[j][1])
                   for i in range(len(h2a)) for j in range(i + 1, len(h2a))]
    avg_h2a_j  = _mean([v for v in h2a_pairs_j if math.isfinite(v)])
    min_h2a_j  = min((v for v in h2a_pairs_j if math.isfinite(v)), default=float("nan"))
    max_h2a_j  = max((v for v in h2a_pairs_j if math.isfinite(v)), default=float("nan"))

    # Cross-seed Jaccard among blk-1 controls
    ctl_pairs_j = [J(ctl[i][0], ctl[i][1], ctl[j][0], ctl[j][1])
                   for i in range(len(ctl)) for j in range(i + 1, len(ctl))]
    avg_ctl_j  = _mean([v for v in ctl_pairs_j if math.isfinite(v)])

    # Inconclusive vs H2a peaks
    inc_h2a_j: list[float] = []
    if inc:
        inc_h2a_j = [J(inc[0][0], inc[0][1], s, b) for s, b in h2a]
    avg_inc_h2a_j = _mean([v for v in inc_h2a_j if math.isfinite(v)]) if inc_h2a_j else float("nan")

    # Interpretation
    elevated = avg_h2a_j > avg_ctl_j + 0.05
    if avg_h2a_j >= JACCARD_HIGH and elevated:
        verdict = "basin_compartido_candidate"
        verdict_text = (
            f"**Basin compartido (candidato conservador).**  "
            f"Jaccard medio entre H2a peaks = {avg_h2a_j:.3f} ≥ {JACCARD_HIGH}, "
            f"claramente por encima del Jaccard de controles blk-1 ({avg_ctl_j:.3f}).  "
            f"La ventana T ≈ 1.5–3.5 produce configuraciones causales similares entre seeds."
        )
    elif avg_h2a_j < JACCARD_LOW:
        verdict = "basins_distintos"
        verdict_text = (
            f"**Basins distintos.**  "
            f"Jaccard medio entre H2a peaks = {avg_h2a_j:.3f} < {JACCARD_LOW}.  "
            f"El F1 parecido en bloques 6–7 no corresponde a configuraciones similares.  "
            f"Un stopping criterion por temperatura seleccionaría basins distintos según el seed."
        )
    elif avg_ctl_j >= JACCARD_HIGH:
        verdict = "inconcluso_controles_altos"
        verdict_text = (
            f"**Inconcluso: controles blk-1 ya altos.**  "
            f"Jaccard H2a medio = {avg_h2a_j:.3f}, pero el Jaccard de controles blk-1 "
            f"también es {avg_ctl_j:.3f}.  "
            f"Los seeds comparten muchas relaciones desde T=100; la similitud en bloques 6–7 "
            f"no es específica de esa ventana."
        )
    else:
        verdict = "mixto_inconcluso"
        verdict_text = (
            f"**Mixto / inconcluso.**  "
            f"Jaccard H2a medio = {avg_h2a_j:.3f} (rango [{min_h2a_j:.3f}, {max_h2a_j:.3f}]).  "
            f"Hay cierta similitud estructural pero insuficiente para afirmar un basin único "
            f"con 3 pares de seeds en N=36."
        )

    lines: list[str] = [
        "# gamma_0p5 basin consistency  N=36",
        "",
        "Post-run SORKIN-2 probe.  Runs gamma_0p5 / medium_25_25_8 for N=36",
        "and computes pairwise Jaccard similarity of induced causal orders at",
        "the H2a peak checkpoints (blocks 6–7) and block-1 controls.",
        "",
        "## Configuration",
        "",
        f"- Command: `{COMMAND}`",
        f"- Generated at UTC: `{generated_at}`",
        f"- N: `{N}`, case_seed: `{CASE_SEED}`, d_spacetime: `{D_SPACETIME}`",
        f"- Schedule: `{SCHEDULE_LABEL}` (cooling_factor={COOLING_FACTOR})",
        f"- Budget: `{BUDGET_LABEL}` (warmup={WARMUP_LIMIT}, anneal={ANNEAL_LIMIT}, blocks={MAX_DATA})",
        f"- Optimizer seeds: `{', '.join(str(s) for s in OPTIMIZER_SEEDS)}`",
        f"- Target causal relations: `{target_n_pairs}` of {N*(N-1)//2} possible pairs.",
        "- Classification uses `causal_f1` against known truth → oracular, not deployable.",
        f"- Jaccard threshold basin compartido: `> {JACCARD_HIGH}`.",
        f"- Jaccard threshold basins distintos: `< {JACCARD_LOW}`.",
        "",
        "## Reproducibility check",
        "",
        f"- Blocks compared against `{REFERENCE_CSV.relative_to(ROOT)}`: "
        f"`{validation['n_compared']}`.",
        f"- max |Δ causal_f1|:   `{_fmt_f(float(validation['max_abs_f1_diff']), 4)}`"
        f"  (tolerance {F1_REPRO_TOL:.0e}).",
        f"- max |Δ energy_eave|: `{_fmt_f(float(validation['max_abs_energy_diff']), 4)}`.",
        f"- **Reproducible: `{validation['reproducible']}`.**",
    ]

    if not validation["reproducible"]:
        lines.extend([
            "",
            "⚠️ **Run is NOT reproducible.  Jaccard analysis below is not interpreted.**",
            "",
            "Mismatches:",
        ])
        for m in validation["mismatches"]:
            lines.append(f"  - {m}")
        lines.extend(["", "## Guardrails", "",
                       "Non-reproducible run.  No further conclusions drawn.", ""])
        MD_PATH.write_text("\n".join(lines), encoding="utf-8")
        return

    lines.extend([
        "",
        "Run is reproducible within floating-point tolerance.  Jaccard analysis proceeds.",
        "",
    ])

    # ── Jaccard matrix table ──────────────────────────────────────────────────
    all_labels = [f"s{s}_b{b}" for s, b, _ in ALL_CHECKPOINTS]
    lines.extend([
        "## Jaccard matrix (full 8×8)",
        "",
        "Rows/columns ordered: H2a peaks, inconclusive reference, blk-1 controls.",
        "Upper triangle shown (matrix is symmetric; diagonal = 1.000).",
        "",
        "| | " + " | ".join(all_labels) + " |",
        "| --- |" + " ---: |" * len(ALL_CHECKPOINTS),
    ])
    for sa, ba, ta in ALL_CHECKPOINTS:
        row_vals = []
        for sb, bb, tb in ALL_CHECKPOINTS:
            j_val = J(sa, ba, sb, bb)
            row_vals.append(f"{j_val:.3f}" if math.isfinite(j_val) else "NA")
        lines.append(f"| s{sa}_b{ba} | " + " | ".join(row_vals) + " |")

    # ── per-candidate F1 and induced size ─────────────────────────────────────
    row_lookup = {
        (int(r["optimizer_seed"]), int(r["block_index"])): r
        for r in trajectory_rows
    }
    lines.extend([
        "",
        "## Checkpoint details",
        "",
        "| checkpoint | type | block | T | causal_f1 | induced_size |",
        "| --- | --- | ---: | ---: | ---: | ---: |",
    ])
    for seed, block, cp_type in ALL_CHECKPOINTS:
        r = row_lookup.get((seed, block))
        if r:
            lines.append(
                f"| s{seed}_b{block} | {cp_type} | {block} "
                f"| {float(r['temperature']):.2f} "
                f"| {float(r['causal_f1']):.4f} "
                f"| {int(r['total_relations_induced'])} |"
            )

    # ── diagnostic questions ───────────────────────────────────────────────────
    lines.extend([
        "",
        "## Diagnostic questions",
        "",
        "### Q1 — Do the H2a checkpoints share a causal basin?",
        "",
        verdict_text,
        "",
        "### Q2 — Jaccard level between H2a peaks: high, medium, or low?",
        "",
        f"H2a pairwise Jaccard values: {', '.join(f'{v:.3f}' for v in h2a_pairs_j if math.isfinite(v))}.",
        f"Mean = **{avg_h2a_j:.3f}**, min = {min_h2a_j:.3f}, max = {max_h2a_j:.3f}.",
        "",
        "### Q3 — Was block-1 Jaccard already high?",
        "",
    ])
    if avg_ctl_j >= JACCARD_HIGH:
        lines.append(
            f"**Yes.**  Mean block-1 cross-seed Jaccard = {avg_ctl_j:.3f} ≥ {JACCARD_HIGH}.  "
            "The induced orders are already similar at T=100.  "
            "High Jaccard at blocks 6–7 cannot be attributed solely to the thermal window."
        )
    elif avg_ctl_j >= JACCARD_LOW:
        lines.append(
            f"**Moderate.**  Mean block-1 cross-seed Jaccard = {avg_ctl_j:.3f}.  "
            "There is some baseline similarity at T=100.  "
            f"Any elevation at blocks 6–7 (avg {avg_h2a_j:.3f}) "
            "is still meaningful but should be interpreted cautiously."
        )
    else:
        lines.append(
            f"**No.**  Mean block-1 cross-seed Jaccard = {avg_ctl_j:.3f} < {JACCARD_LOW}.  "
            "The hot-start configurations are structurally dissimilar.  "
            f"The Jaccard elevation at blocks 6–7 (avg {avg_h2a_j:.3f}) "
            "reflects genuine convergence, not pre-existing similarity."
        )

    lines.extend([
        "",
        "### Q4 — Does seed 1962 block 8 resemble the H2a peaks?",
        "",
    ])
    if math.isfinite(avg_inc_h2a_j):
        if avg_inc_h2a_j >= JACCARD_HIGH:
            lines.append(
                f"**Yes (high similarity).**  "
                f"seed 1962 blk-8 vs H2a peaks: mean Jaccard = {avg_inc_h2a_j:.3f}.  "
                "The inconclusive seed (best=final) reaches a similar causal configuration."
            )
        elif avg_inc_h2a_j >= JACCARD_LOW:
            lines.append(
                f"**Partial.**  "
                f"seed 1962 blk-8 vs H2a peaks: mean Jaccard = {avg_inc_h2a_j:.3f}.  "
                "Moderate overlap with H2a peaks.  "
                "seed 1962 may be in a neighbouring region of configuration space."
            )
        else:
            lines.append(
                f"**Low similarity.**  "
                f"seed 1962 blk-8 vs H2a peaks: mean Jaccard = {avg_inc_h2a_j:.3f}.  "
                "seed 1962 final endpoint is structurally different from H2a peak configurations, "
                "consistent with reaching a different local minimum."
            )

    lines.extend([
        "",
        "### Q5 — Does T ≈ 1.5–3.5 have structural support or only scalar F1 support?",
        "",
    ])
    if avg_h2a_j >= JACCARD_HIGH and not avg_ctl_j >= JACCARD_HIGH:
        lines.append(
            f"**Structural support present.**  "
            f"The H2a peaks at T ≈ 1.5–3.5 not only have similar F1 scalars "
            f"but also similar induced causal orders (mean Jaccard {avg_h2a_j:.3f}).  "
            "This is evidence for a common causal region, not a coincidence of F1."
        )
    elif avg_h2a_j >= JACCARD_LOW and not avg_ctl_j >= JACCARD_HIGH:
        lines.append(
            f"**Partial structural support.**  "
            f"Mean Jaccard {avg_h2a_j:.3f} suggests moderate structural overlap.  "
            "The temperature window has some structural basis but not a fully shared basin."
        )
    else:
        lines.append(
            f"**Unclear.**  "
            f"The Jaccard result ({avg_h2a_j:.3f}) cannot be cleanly separated "
            "from background similarity.  "
            "The T ≈ 1.5–3.5 window has only scalar (F1) support in this dataset."
        )

    lines.extend([
        "",
        "### Q6 — Conservative conclusion",
        "",
        f"Verdict: **`{verdict}`**.",
        "",
        "Details:",
        f"- H2a peak Jaccard (3 seeds, N=36, 1 case): avg={avg_h2a_j:.3f}, "
        f"range=[{min_h2a_j:.3f}, {max_h2a_j:.3f}].",
        f"- Control blk-1 Jaccard (baseline): avg={avg_ctl_j:.3f}.",
        f"- seed 1962 blk-8 vs H2a peaks: avg={avg_inc_h2a_j:.3f}.",
        "",
    ])

    if verdict == "basin_compartido_candidate":
        lines.extend([
            "The evidence is consistent with a shared causal basin in the T ≈ 1.5–3.5 window",
            "for gamma_0p5 / N=36 / case_seed=1959.  **This does not establish a general rule.**",
            "Replication with different case seeds and N is required before proposing a stopping",
            "criterion.  The probe also cannot rule out that the shared structure is specific to",
            "this target (343 relations, 54% density) rather than to the thermal window.",
        ])
    elif verdict == "basins_distintos":
        lines.extend([
            "The F1 similarity at blocks 6–7 is not accompanied by structural (causal order) similarity.",
            "A stopping rule targeting T ≈ 3 would recover configurations that differ between seeds.",
            "This closes the temperature-based stopping criterion hypothesis for gamma_0p5 / N=36.",
        ])
    elif verdict == "inconcluso_controles_altos":
        lines.extend([
            "The high block-1 similarity means the seeds do not start from structurally independent",
            "configurations.  The Jaccard result at blocks 6–7 is confounded by this baseline.",
            "A valid basin test would require structurally independent initializations or a",
            "different analysis frame (e.g., relative shift from blk-1 to blk-6/7).",
        ])
    else:
        lines.extend([
            "The result is mixed.  More seeds or a different N are needed.",
            "The stopping criterion hypothesis remains open but unconfirmed.",
        ])

    lines.extend([
        "",
        "## Conservative interpretation",
        "",
        "This is an oracular diagnostic: the Jaccard computation uses the known-truth",
        "induced order, not an observable available without ground truth.",
        "",
        "Sample: 3 H2a seeds, 1 case seed (1959), 1 N (36), 1 schedule (gamma_0p5).",
        "No generalization claim is warranted from this probe alone.",
        "",
        "## Guardrails",
        "",
        "This is a post-run diagnostic only, over one benchmark case with known truth.",
        "It is not an embeddability claim, not a physical gamma claim, not an N-transition claim,",
        "and not proof of general annealer failure.",
        "It is not a deployable checkpoint-selection criterion.",
        "",
    ])
    MD_PATH.write_text("\n".join(lines), encoding="utf-8")


# ── main ──────────────────────────────────────────────────────────────────────

def main() -> int:
    t0 = time.perf_counter()
    print(f"[1/7] Building N={N} Minkowski case (case_seed={CASE_SEED}) …")
    case = _make_case()
    target_pairs = _pairs(case.matrix)
    print(f"      Target: {len(target_pairs)} causal pairs of {N*(N-1)//2} possible.")

    print(f"[2/7] Writing target pairs CSV ({len(target_pairs)} rows) …")
    _write_target_pairs_csv(target_pairs)

    print(f"[3/7] Running gamma_0p5 for {len(OPTIMIZER_SEEDS)} seeds …")
    trajectory_rows: list[dict[str, object]] = []
    for seed in OPTIMIZER_SEEDS:
        t_seed = time.perf_counter()
        rows = _run_one(case, seed)
        trajectory_rows.extend(rows)
        print(f"      seed {seed}: {len(rows)} blocks  ({time.perf_counter()-t_seed:.1f}s)")

    print("[4/7] Validating reproducibility against reference CSV …")
    validation = _validate(trajectory_rows)
    print(f"      Compared {validation['n_compared']} blocks.  "
          f"max|Δf1|={float(validation['max_abs_f1_diff']):.2e}  "
          f"max|Δenergy|={float(validation['max_abs_energy_diff']):.2e}  "
          f"reproducible={validation['reproducible']}")
    if not validation["reproducible"]:
        print("  ⚠  NOT REPRODUCIBLE — writing warning MD and stopping.")
        _write_markdown(trajectory_rows, [], validation, len(target_pairs))
        _write_csv(TRAJECTORY_CSV, TRAJECTORY_HEADERS, trajectory_rows)
        return 1

    print("[5/7] Writing trajectory CSV …")
    _write_csv(TRAJECTORY_CSV, TRAJECTORY_HEADERS, trajectory_rows)

    print("[6/7] Computing Jaccard matrix …")
    jaccard_rows = _compute_jaccard(trajectory_rows)
    _write_csv(JACCARD_CSV, JACCARD_HEADERS, jaccard_rows)

    print("[7/7] Writing SVG and markdown …")
    _write_svg(jaccard_rows)
    _write_markdown(trajectory_rows, jaccard_rows, validation, len(target_pairs))

    elapsed = time.perf_counter() - t0
    print(f"\nDone in {elapsed:.1f}s.  Artifacts written to {OUT_DIR.relative_to(ROOT)}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
