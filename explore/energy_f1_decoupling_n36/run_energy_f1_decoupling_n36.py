#!/usr/bin/env python3
"""Energy–F1 decoupling diagnostic for SORKIN-2 N=36.

Reads the already-generated oracle checkpoint ceiling CSV and classifies
each (schedule, seed) group into one of two candidate failure modes:

  H2a_over_annealing_candidate
      The best-F1 checkpoint has *higher* energy than the final endpoint.
      The annealer cooled through the good causal region and froze in a
      lower-energy but causally worse state.

  H2b_escape_or_nonconvergence_candidate
      The best-F1 checkpoint has *lower* energy than the final endpoint.
      The system visited a good causal region but later escaped to states
      of higher energy and worse causal quality.

  inconclusive
      Neither condition holds (e.g. best == final, or delta_F1 == 0).

This script does NOT run the annealer.  It is a pure post-run diagnostic
that uses causal_f1 against the known-truth target, so classification is
oracular and not deployable without ground truth.
"""

from __future__ import annotations

import csv
import math
import statistics
import sys
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

OUT_DIR = Path(__file__).resolve().parent
SOURCE_CSV = (
    ROOT
    / "explore"
    / "oracle_checkpoint_ceiling_n36"
    / "oracle_checkpoint_ceiling_n36.csv"
)
CSV_PATH = OUT_DIR / "energy_f1_decoupling_n36.csv"
BY_SCHEDULE_CSV_PATH = OUT_DIR / "energy_f1_decoupling_n36_by_schedule.csv"
MD_PATH = OUT_DIR / "energy_f1_decoupling_n36.md"
SVG_PATH = OUT_DIR / "energy_f1_decoupling_n36.svg"
COMMAND = "python3 explore/energy_f1_decoupling_n36/run_energy_f1_decoupling_n36.py"

GROUP_HEADERS = (
    "optimizer_seed",
    "schedule_label",
    "cooling_factor",
    "budget_label",
    "final_causal_f1",
    "best_checkpoint_causal_f1",
    "delta_best_minus_final",
    "final_energy_eave",
    "best_checkpoint_energy_eave",
    "min_energy_eave",
    "delta_energy_best_minus_final",
    "block_index_final",
    "block_index_best_causal_f1",
    "block_index_min_energy",
    "best_matches_min_energy",
    "inferred_failure_mode",
)

BY_SCHEDULE_HEADERS = (
    "schedule_label",
    "cooling_factor",
    "budget_label",
    "n_groups",
    "avg_delta_energy_best_minus_final",
    "median_delta_energy_best_minus_final",
    "avg_delta_best_minus_final",
    "avg_best_checkpoint_causal_f1",
    "avg_final_causal_f1",
    "count_H2a_over_annealing_candidate",
    "count_H2b_escape_or_nonconvergence_candidate",
    "count_inconclusive",
    "count_best_matches_min_energy",
    "avg_block_index_best_causal_f1",
    "avg_block_index_min_energy",
)


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------

def _fmt(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return f"{value:.10g}" if math.isfinite(value) else "NA"
    return str(value)


def _fmt_f(value: float, digits: int = 6) -> str:
    return f"{value:.{digits}g}" if math.isfinite(value) else "NA"


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def _read_source() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with SOURCE_CSV.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for raw in reader:
            rows.append({
                "optimizer_seed": int(raw["optimizer_seed"]),
                "schedule_label": str(raw["schedule_label"]),
                "cooling_factor": float(raw["cooling_factor"]),
                "budget_label": str(raw["budget_label"]),
                "final_causal_f1": float(raw["final_causal_f1"]),
                "best_checkpoint_causal_f1": float(raw["best_checkpoint_causal_f1"]),
                "delta_best_minus_final": float(raw["delta_best_minus_final"]),
                "final_energy_eave": float(raw["final_energy_eave"]),
                "best_checkpoint_energy_eave": float(raw["best_checkpoint_energy_eave"]),
                "min_energy_eave": float(raw["min_energy_eave"]),
                "block_index_final": int(raw["block_index_final"]),
                "block_index_best_causal_f1": int(raw["block_index_best_causal_f1"]),
                "block_index_min_energy": int(raw["block_index_min_energy"]),
                "best_matches_min_energy": raw["best_matches_min_energy"].strip().lower() == "true",
            })
    return rows


def _write_csv(
    path: Path,
    headers: tuple[str, ...],
    rows: list[dict[str, object]],
) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=headers)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: _fmt(row[k]) for k in headers})


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------

def _infer_failure_mode(
    delta_energy: float,
    delta_f1: float,
) -> str:
    """Conservative two-regime classification.

    H2a: best-F1 checkpoint at higher energy than final, F1 is better.
         Annealer over-annealed past the good causal region.

    H2b: best-F1 checkpoint at lower energy than final, F1 is better.
         System escaped from the good causal region to higher-energy states.

    inconclusive: best == final (delta_f1 == 0) or the energy signal is
                  ambiguous (delta_energy == 0 with delta_f1 > 0 would be
                  unusual; treated as inconclusive).
    """
    best_is_better = delta_f1 > 0.0
    energy_positive = delta_energy > 0.0
    energy_negative = delta_energy < 0.0

    if best_is_better and energy_positive:
        return "H2a_over_annealing_candidate"
    if best_is_better and energy_negative:
        return "H2b_escape_or_nonconvergence_candidate"
    return "inconclusive"


# ---------------------------------------------------------------------------
# Group-level computation
# ---------------------------------------------------------------------------

def _build_group_rows(source: list[dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for r in source:
        delta_energy = r["best_checkpoint_energy_eave"] - r["final_energy_eave"]  # type: ignore[operator]
        rows.append({
            "optimizer_seed": r["optimizer_seed"],
            "schedule_label": r["schedule_label"],
            "cooling_factor": r["cooling_factor"],
            "budget_label": r["budget_label"],
            "final_causal_f1": r["final_causal_f1"],
            "best_checkpoint_causal_f1": r["best_checkpoint_causal_f1"],
            "delta_best_minus_final": r["delta_best_minus_final"],
            "final_energy_eave": r["final_energy_eave"],
            "best_checkpoint_energy_eave": r["best_checkpoint_energy_eave"],
            "min_energy_eave": r["min_energy_eave"],
            "delta_energy_best_minus_final": delta_energy,
            "block_index_final": r["block_index_final"],
            "block_index_best_causal_f1": r["block_index_best_causal_f1"],
            "block_index_min_energy": r["block_index_min_energy"],
            "best_matches_min_energy": r["best_matches_min_energy"],
            "inferred_failure_mode": _infer_failure_mode(
                delta_energy,
                float(r["delta_best_minus_final"]),  # type: ignore[arg-type]
            ),
        })
    # Sort by cooling_factor asc, then seed asc
    return sorted(rows, key=lambda x: (float(x["cooling_factor"]), int(x["optimizer_seed"])))


# ---------------------------------------------------------------------------
# Schedule-level aggregation
# ---------------------------------------------------------------------------

def _mean(values: list[float]) -> float:
    return sum(values) / len(values)


def _build_by_schedule(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, float, str], list[dict[str, object]]] = {}
    for row in rows:
        key = (str(row["schedule_label"]), float(row["cooling_factor"]), str(row["budget_label"]))
        grouped.setdefault(key, []).append(row)

    out: list[dict[str, object]] = []
    for (schedule_label, cooling_factor, budget_label), group in sorted(
        grouped.items(), key=lambda kv: kv[0][1]
    ):
        delta_energies = [float(r["delta_energy_best_minus_final"]) for r in group]
        out.append({
            "schedule_label": schedule_label,
            "cooling_factor": cooling_factor,
            "budget_label": budget_label,
            "n_groups": len(group),
            "avg_delta_energy_best_minus_final": _mean(delta_energies),
            "median_delta_energy_best_minus_final": statistics.median(delta_energies),
            "avg_delta_best_minus_final": _mean([float(r["delta_best_minus_final"]) for r in group]),
            "avg_best_checkpoint_causal_f1": _mean([float(r["best_checkpoint_causal_f1"]) for r in group]),
            "avg_final_causal_f1": _mean([float(r["final_causal_f1"]) for r in group]),
            "count_H2a_over_annealing_candidate": sum(
                1 for r in group if r["inferred_failure_mode"] == "H2a_over_annealing_candidate"
            ),
            "count_H2b_escape_or_nonconvergence_candidate": sum(
                1 for r in group if r["inferred_failure_mode"] == "H2b_escape_or_nonconvergence_candidate"
            ),
            "count_inconclusive": sum(
                1 for r in group if r["inferred_failure_mode"] == "inconclusive"
            ),
            "count_best_matches_min_energy": sum(
                1 for r in group if bool(r["best_matches_min_energy"])
            ),
            "avg_block_index_best_causal_f1": _mean([float(r["block_index_best_causal_f1"]) for r in group]),
            "avg_block_index_min_energy": _mean([float(r["block_index_min_energy"]) for r in group]),
        })
    return out


# ---------------------------------------------------------------------------
# SVG
# ---------------------------------------------------------------------------

def _write_svg(rows: list[dict[str, object]]) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("matplotlib is required for SVG output") from exc

    # Colour and marker per schedule
    schedule_order = ["gamma_0p5", "gamma_0p8", "gamma_0p9", "gamma_0p95"]
    palette = {
        "gamma_0p5":  ("#1f77b4", "o"),
        "gamma_0p8":  ("#2ca02c", "s"),
        "gamma_0p9":  ("#ff7f0e", "^"),
        "gamma_0p95": ("#d62728", "D"),
    }

    fig, ax = plt.subplots(figsize=(7.2, 5.4))

    plotted_labels: set[str] = set()
    for row in rows:
        label = str(row["schedule_label"])
        color, marker = palette.get(label, ("#888888", "x"))
        legend_label = label if label not in plotted_labels else None
        plotted_labels.add(label)
        ax.scatter(
            float(row["delta_energy_best_minus_final"]),
            float(row["delta_best_minus_final"]),
            color=color,
            marker=marker,
            s=52,
            label=legend_label,
            zorder=3,
        )

    ax.axvline(0.0, color="#555555", linewidth=1.2, linestyle="--", zorder=2)
    ax.axhline(0.0, color="#555555", linewidth=1.2, linestyle=":",  zorder=2)

    ax.set_xlabel("Δ energy  (best_checkpoint − final)", fontsize=11)
    ax.set_ylabel("Δ causal F1  (best_checkpoint − final)", fontsize=11)
    ax.set_title("Energy–F1 decoupling: failure-mode quadrants  (N=36)", fontsize=11)

    # Annotate quadrants
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.text(
        xlim[1] * 0.97, ylim[1] * 0.97,
        "H2a\n(over-annealing)",
        ha="right", va="top", fontsize=8, color="#1f77b4",
        style="italic",
    )
    ax.text(
        xlim[0] * 0.97, ylim[1] * 0.97,
        "H2b\n(escape/non-conv.)",
        ha="left", va="top", fontsize=8, color="#d62728",
        style="italic",
    )

    ax.legend(fontsize=9, framealpha=0.85)
    ax.grid(True, alpha=0.28)
    fig.tight_layout()
    fig.savefig(SVG_PATH, format="svg")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Markdown
# ---------------------------------------------------------------------------

def _write_markdown(
    rows: list[dict[str, object]],
    by_schedule: list[dict[str, object]],
) -> None:
    generated_at = datetime.now(timezone.utc).isoformat(timespec="seconds")

    n_total = len(rows)
    n_h2a = sum(1 for r in rows if r["inferred_failure_mode"] == "H2a_over_annealing_candidate")
    n_h2b = sum(1 for r in rows if r["inferred_failure_mode"] == "H2b_escape_or_nonconvergence_candidate")
    n_inc = sum(1 for r in rows if r["inferred_failure_mode"] == "inconclusive")

    # Fast vs slow schedule separation
    fast_schedules = {"gamma_0p5", "gamma_0p8"}
    slow_schedules = {"gamma_0p9", "gamma_0p95"}
    fast_rows = [r for r in rows if str(r["schedule_label"]) in fast_schedules]
    slow_rows = [r for r in rows if str(r["schedule_label"]) in slow_schedules]

    fast_h2a = sum(1 for r in fast_rows if r["inferred_failure_mode"] == "H2a_over_annealing_candidate")
    fast_h2b = sum(1 for r in fast_rows if r["inferred_failure_mode"] == "H2b_escape_or_nonconvergence_candidate")
    slow_h2a = sum(1 for r in slow_rows if r["inferred_failure_mode"] == "H2a_over_annealing_candidate")
    slow_h2b = sum(1 for r in slow_rows if r["inferred_failure_mode"] == "H2b_escape_or_nonconvergence_candidate")

    # avg delta_energy per regime
    fast_avg_de = _mean([float(r["delta_energy_best_minus_final"]) for r in fast_rows]) if fast_rows else float("nan")
    slow_avg_de = _mean([float(r["delta_energy_best_minus_final"]) for r in slow_rows]) if slow_rows else float("nan")

    all_delta_energies = [float(r["delta_energy_best_minus_final"]) for r in rows]
    sign_consistency = sum(1 for d in all_delta_energies if d > 0) / n_total if n_total else 0.0

    lines: list[str] = [
        "# Energy–F1 decoupling diagnostic  N=36",
        "",
        "Post-run SORKIN-2 diagnostic.  Reads the oracle checkpoint ceiling CSV",
        "and classifies each (schedule, seed) group by the sign of",
        "`delta_energy_best_minus_final = best_checkpoint_energy_eave − final_energy_eave`",
        "relative to `delta_best_minus_final` (causal F1 gain of the oracle checkpoint).",
        "",
        "## Configuration",
        "",
        f"- Command: `{COMMAND}`",
        f"- Generated at UTC: `{generated_at}`",
        f"- Source CSV: `{SOURCE_CSV.relative_to(ROOT)}`",
        f"- Main output CSV: `{CSV_PATH.relative_to(ROOT)}`",
        f"- By-schedule CSV: `{BY_SCHEDULE_CSV_PATH.relative_to(ROOT)}`",
        f"- SVG: `{SVG_PATH.relative_to(ROOT)}`",
        "- This script does not run the annealer.",
        "- Classification uses `causal_f1` against known-truth target → oracular, not deployable.",
        "",
        "## Failure-mode classification (per group)",
        "",
        "| seed | schedule | γ | budget | final F1 | best F1 | ΔF1 | final E | best E | ΔE | blk best | blk minE | best=minE | mode |",
        "| ---: | --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |",
    ]

    mode_short = {
        "H2a_over_annealing_candidate": "**H2a**",
        "H2b_escape_or_nonconvergence_candidate": "**H2b**",
        "inconclusive": "inconclusive",
    }
    for row in rows:
        lines.append(
            "| {seed} | {sched} | {gamma} | {budget} | {ff1:.4f} | {bf1:.4f} |"
            " {df1:+.4f} | {fe:.1f} | {be:.1f} | {de:+.1f} | {bb} | {mb} | {bm} | {mode} |".format(
                seed=row["optimizer_seed"],
                sched=row["schedule_label"],
                gamma=_fmt_f(float(row["cooling_factor"])),
                budget=row["budget_label"],
                ff1=float(row["final_causal_f1"]),
                bf1=float(row["best_checkpoint_causal_f1"]),
                df1=float(row["delta_best_minus_final"]),
                fe=float(row["final_energy_eave"]),
                be=float(row["best_checkpoint_energy_eave"]),
                de=float(row["delta_energy_best_minus_final"]),
                bb=row["block_index_best_causal_f1"],
                mb=row["block_index_min_energy"],
                bm="yes" if bool(row["best_matches_min_energy"]) else "no",
                mode=mode_short.get(str(row["inferred_failure_mode"]), str(row["inferred_failure_mode"])),
            )
        )

    lines.extend([
        "",
        "## By-schedule aggregates",
        "",
        "| schedule | γ | budget | N | avg ΔE | med ΔE | avg ΔF1 | avg best F1 | avg final F1"
        " | H2a | H2b | inconc | best=minE | avg blk best | avg blk minE |",
        "| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ])
    for s in by_schedule:
        lines.append(
            "| {sched} | {gamma} | {budget} | {n} | {ade:+.1f} | {mde:+.1f} |"
            " {adf:.4f} | {abf:.4f} | {aff:.4f} | {h2a} | {h2b} | {inc} | {bme} | {abb:.2f} | {amb:.2f} |".format(
                sched=s["schedule_label"],
                gamma=_fmt_f(float(s["cooling_factor"])),
                budget=s["budget_label"],
                n=s["n_groups"],
                ade=float(s["avg_delta_energy_best_minus_final"]),
                mde=float(s["median_delta_energy_best_minus_final"]),
                adf=float(s["avg_delta_best_minus_final"]),
                abf=float(s["avg_best_checkpoint_causal_f1"]),
                aff=float(s["avg_final_causal_f1"]),
                h2a=s["count_H2a_over_annealing_candidate"],
                h2b=s["count_H2b_escape_or_nonconvergence_candidate"],
                inc=s["count_inconclusive"],
                bme=s["count_best_matches_min_energy"],
                abb=float(s["avg_block_index_best_causal_f1"]),
                amb=float(s["avg_block_index_min_energy"]),
            )
        )

    lines.extend([
        "",
        "## Diagnostic questions",
        "",
        "### Q1 — Is there a clear separation between fast and slow schedules?",
        "",
    ])

    if fast_h2a >= 3 and slow_h2b >= 3:
        separation = (
            f"**Yes, with conservative support.**  "
            f"Fast schedules (gamma_0p5, gamma_0p8): {fast_h2a}/{len(fast_rows)} groups → H2a, "
            f"avg ΔE = {fast_avg_de:+.1f}.  "
            f"Slow schedules (gamma_0p9, gamma_0p95): {slow_h2b}/{len(slow_rows)} groups → H2b, "
            f"avg ΔE = {slow_avg_de:+.1f}.  "
            f"The sign of ΔE is predominantly positive for fast schedules and predominantly "
            f"negative for slow schedules.  With only 4 groups per schedule the separation is "
            f"suggestive, not conclusive."
        )
    elif fast_h2a >= 2 or slow_h2b >= 2:
        separation = (
            f"**Partial separation.**  "
            f"Fast schedules: {fast_h2a}/{len(fast_rows)} H2a, {fast_h2b}/{len(fast_rows)} H2b.  "
            f"Slow schedules: {slow_h2a}/{len(slow_rows)} H2a, {slow_h2b}/{len(slow_rows)} H2b.  "
            f"There is a directional tendency but not a clean split at this sample size."
        )
    else:
        separation = (
            f"**No clear separation.**  "
            f"Fast schedules: {fast_h2a}/{len(fast_rows)} H2a.  "
            f"Slow schedules: {slow_h2b}/{len(slow_rows)} H2b.  "
            f"The hypothesised two-regime structure is not supported by these {n_total} groups."
        )

    lines.append(separation)
    lines.append("")

    lines.extend([
        "### Q2 — How many groups fall into H2a?",
        "",
        f"**{n_h2a} of {n_total} groups** are classified as `H2a_over_annealing_candidate`.",
        "",
        "### Q3 — How many groups fall into H2b?",
        "",
        f"**{n_h2b} of {n_total} groups** are classified as `H2b_escape_or_nonconvergence_candidate`.",
        f"**{n_inc} of {n_total} groups** are `inconclusive` (best == final or ΔE = 0).",
        "",
        "### Q4 — Is the classification by sign of ΔE robust or marginal?",
        "",
    ])

    # Robustness assessment
    small_delta_threshold = 5.0
    near_zero = [r for r in rows if abs(float(r["delta_energy_best_minus_final"])) < small_delta_threshold]
    n_near = len(near_zero)
    if n_near == 0:
        robustness = (
            "All classified groups have |ΔE| ≥ {t:.0f} energy units.  "
            "No cases are close to the ΔE = 0 boundary; the sign-based classification "
            "appears robust to small measurement variation within this dataset.".format(t=small_delta_threshold)
        )
    else:
        near_modes = [str(r["inferred_failure_mode"]) for r in near_zero]
        robustness = (
            f"{n_near} group(s) have |ΔE| < {small_delta_threshold:.0f} energy units "
            f"(modes: {near_modes}).  "
            "These are near the boundary and their classification is fragile."
        )

    lines.append(robustness)
    lines.append("")

    lines.extend([
        "### Q5 — Is it still reasonable to treat the 16 groups as a single population?",
        "",
    ])

    if fast_h2a >= 3 and slow_h2b >= 3:
        pooling = (
            "**No.**  The failure-mode split correlates strongly with cooling speed.  "
            "Pooling all 16 groups flattens two qualitatively distinct dynamics: "
            "one in which the annealer over-cools through the good causal region (H2a) "
            "and one in which it visits then escapes that region (H2b).  "
            "Treating the 16 groups as a homogeneous sample would produce averages that "
            "describe neither regime accurately."
        )
    else:
        pooling = (
            "**Uncertain.**  The failure-mode separation is partial or absent at this sample size.  "
            "Pooling is not clearly wrong but may obscure regime structure.  "
            "More seeds or schedules would clarify."
        )

    lines.append(pooling)
    lines.append("")

    lines.extend([
        "### Q6 — What does this imply for the next probe?",
        "",
        "If the H2a/H2b separation is confirmed:",
        "",
        "- **H2a probes** should focus on the temperature window during which causal F1 peaks "
          "and whether it corresponds to a stable or transient basin.  "
          "A trajectory-dump probe for fast schedules (gamma_0p5 / gamma_0p8) recording "
          "causal F1 at each block across multiple seeds would establish whether the good-causal "
          "temperature range is reproducible.",
        "",
        "- **H2b probes** should focus on why the annealer escapes from the initially-good state.  "
          "This is a different failure: the energy landscape is not monotone in causal quality "
          "even at the start.  A probe recording which relations flip between the best block "
          "and later blocks would characterise the escape mechanism.",
        "",
        "- **Do not design a single stopping criterion** that targets both regimes simultaneously.  "
          "The two regimes require different interventions.",
        "",
        "## Conservative interpretation",
        "",
        "This diagnostic uses `causal_f1` against the known-truth partial order to identify "
        "the oracle-best checkpoint.  It is therefore not a deployable selection criterion "
        "for truth-free cases.",
        "",
        "The H2a / H2b classification is a hypothesis with `N_groups = 16` and "
        "only 4 seeds per schedule.  It may reflect true regime structure or schedule-specific "
        "noise.  The classification should be treated as a diagnostic candidate, not a "
        "confirmed physical finding.",
        "",
        "The diagnostic does not distinguish between:",
        "  - a fundamentally misaligned energy function (same function for all regimes),",
        "  - schedule-specific convergence behaviour (different dynamics per cooling rate).",
        "Both could produce the observed pattern.  Next probes should hold the energy function "
        "fixed (it is historical) and vary the observable: trajectory timing, relation-flip "
        "rate, or seed-to-seed configuration similarity.",
        "",
        "## Guardrails",
        "",
        "This is a post-run diagnostic only, over benchmark cases with known truth.",
        "It is not an embeddability claim, not a physical gamma claim, not an N-transition claim,",
        "and not proof of general annealer failure.",
        "It is not a deployable criterion for truth-free cases.",
        "",
    ])

    MD_PATH.write_text("\n".join(lines), encoding="utf-8")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    source = _read_source()
    rows = _build_group_rows(source)
    by_schedule = _build_by_schedule(rows)

    _write_csv(CSV_PATH, GROUP_HEADERS, rows)
    _write_csv(BY_SCHEDULE_CSV_PATH, BY_SCHEDULE_HEADERS, by_schedule)
    _write_svg(rows)
    _write_markdown(rows, by_schedule)

    n_h2a = sum(1 for r in rows if r["inferred_failure_mode"] == "H2a_over_annealing_candidate")
    n_h2b = sum(1 for r in rows if r["inferred_failure_mode"] == "H2b_escape_or_nonconvergence_candidate")
    n_inc = sum(1 for r in rows if r["inferred_failure_mode"] == "inconclusive")
    print(f"energy_f1_decoupling_n36: {len(rows)} groups  →  H2a={n_h2a}  H2b={n_h2b}  inconclusive={n_inc}")
    print(f"Artifacts written to {OUT_DIR.relative_to(ROOT)}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
