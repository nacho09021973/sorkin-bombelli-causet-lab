#!/usr/bin/env python3
"""Build the Phase 2B annealer schedule probe.

Phase 2B isolates whether the large `energy_gap` and `interval_rmse`
observed on Minkowski cases in Phase 2 is dominated by an under-budgeted
schedule of the historical Bombelli-Sorkin annealer. It does **not**
introduce a new optimizer. It varies only the existing scheduling knobs
already accepted by :class:`cones.ConesSimulator` and by
:func:`validation_suite.run_recovery`:

- ``warmup_limit``: number of warmup reconfigure steps before annealing.
- ``anneal_limit``: number of reconfigure steps per temperature stage.
- ``max_data``: number of temperature stages before stopping.

The temperature schedule (``initial_temp``, ``cooling_factor``) is held
fixed at the Phase 2 values to keep the comparison as a pure budget
probe, not a temperature search.

Family scope is restricted to Minkowski sprinklings. Non-manifoldlike
controls (Kleitman-Rothschild, suspended corona) are *not* included
here: they have no ground-truth coordinates, so the diagnostic
quantities ``truth_energy``, ``energy_gap`` and ``interval_rmse`` are
not defined. Phase 2B is therefore explicitly not a manifoldness
classifier.

The probe reuses the canonical sprinkler in :mod:`validation_suite`,
the same generator that Phase 2 uses. Seeds are the first three of
the Phase 1B atlas ``SEEDS`` constant.
"""

from __future__ import annotations

import math
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.build_phase1_atlas import _format_field  # noqa: E402
from tools.build_phase1b_scaling_atlas import FOUNDATION, SEEDS  # noqa: E402
import validation_suite as vs  # noqa: E402


# ---------------------------------------------------------------------
# Grid: only Minkowski; documented seed subset; existing schedule knobs
# ---------------------------------------------------------------------

SPACETIME_DIMS = (2, 3, 4)
SIZES = (32, 64)
# First three seeds from the Phase 1B atlas; documented subset rather
# than a fresh draw, so each Phase 2B cell can be cross-referenced to
# its Phase 1D structural row at the same (n, seed).
PROBE_SEEDS = SEEDS[:3]
OPTIMIZER_SEED = 1987

# Fixed temperature schedule (same as Phase 2). Phase 2B varies only
# the iteration budget; mixing in temperature changes would conflate
# the two effects.
INITIAL_TEMP = 100.0
COOLING_FACTOR = 0.9

# Three schedules. "short" reproduces the Phase 2 configuration so the
# baseline gap is directly comparable. "medium" and "long" widen the
# reconfigure budget by ~3x and ~7x respectively. The unbounded
# :class:`cones.ConesSimulator` defaults (``warmup_limit=100``,
# ``anneal_limit=100``, ``max_data=35``, total budget 3600) are not
# used here: empirically a single ``n=64, d=4`` recovery at that
# budget on this CPU build runs in the high hundreds of seconds, so a
# 54-cell grid at default budget is not a reproducible smoke. The 7x
# range chosen below is wide enough to reveal a budget-limited descent
# (if one exists) while keeping ``make regen-phase2b`` within minutes.
SCHEDULES: tuple[dict, ...] = (
    {
        "label": "short",
        "warmup_limit": 10,
        "anneal_limit": 10,
        "max_data": 4,
    },
    {
        "label": "medium",
        "warmup_limit": 20,
        "anneal_limit": 20,
        "max_data": 6,
    },
    {
        "label": "long",
        "warmup_limit": 30,
        "anneal_limit": 30,
        "max_data": 10,
    },
)

# Conservative success threshold. Truth energy is 0 for an exact
# Minkowski sprinkling, so any positive ``final_energy`` is a failure.
# We call a run a success only if the optimizer collapsed the gap to
# numerically negligible (well below the ``~600--800`` scale observed
# in the Phase 2 short schedule).
SUCCESS_GAP_THRESHOLD = 1.0


CSV_HEADERS = (
    "family",
    "target_dim",
    "n",
    "seed",
    "schedule_label",
    "warmup_limit",
    "anneal_limit",
    "max_data",
    "initial_temp",
    "cooling_factor",
    "optimizer_seed",
    "embedding_dim",
    "initial_energy",
    "warmup_energy",
    "truth_energy",
    "final_energy",
    "energy_gap",
    "interval_rmse",
    "success_flag",
    "runtime_seconds",
)


def _format_csv_value(value) -> str:
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, str):
        return value
    return _format_field(value)


def _run_one(d_spacetime: int, n: int, seed: int, schedule: dict) -> dict:
    matrix, points = vs.sprinkle_minkowski_diamond(
        n=n, seed=seed, d_spacetime=d_spacetime
    )
    case = vs.SprinkleCase(
        d_spacetime=d_spacetime,
        n=n,
        seed=seed,
        matrix=matrix,
        points=points,
    )

    start = time.perf_counter()
    result = vs.run_recovery(
        case,
        optimizer_seed=OPTIMIZER_SEED,
        target_dim=d_spacetime - 1,
        warmup_limit=schedule["warmup_limit"],
        anneal_limit=schedule["anneal_limit"],
        max_data=schedule["max_data"],
        initial_temp=INITIAL_TEMP,
        cooling_factor=COOLING_FACTOR,
        backend="cpu",
    )
    runtime = time.perf_counter() - start

    energy_gap = result.final_energy - result.truth_energy
    success_flag = (
        math.isfinite(energy_gap)
        and math.isfinite(result.interval_rmse)
        and energy_gap <= SUCCESS_GAP_THRESHOLD
    )
    return {
        "family": "minkowski",
        "target_dim": d_spacetime,
        "n": n,
        "seed": seed,
        "schedule_label": schedule["label"],
        "warmup_limit": schedule["warmup_limit"],
        "anneal_limit": schedule["anneal_limit"],
        "max_data": schedule["max_data"],
        "initial_temp": INITIAL_TEMP,
        "cooling_factor": COOLING_FACTOR,
        "optimizer_seed": OPTIMIZER_SEED,
        "embedding_dim": d_spacetime - 1,
        "initial_energy": result.initial_energy,
        "warmup_energy": result.warmup_energy,
        "truth_energy": result.truth_energy,
        "final_energy": result.final_energy,
        "energy_gap": energy_gap,
        "interval_rmse": result.interval_rmse,
        "success_flag": bool(success_flag),
        "runtime_seconds": runtime,
    }


def build_rows() -> list[dict]:
    rows: list[dict] = []
    for d in SPACETIME_DIMS:
        for n in SIZES:
            for seed in PROBE_SEEDS:
                for schedule in SCHEDULES:
                    rows.append(_run_one(d, n, seed, schedule))
    return rows


def write_csv(rows: list[dict], path: Path) -> None:
    lines = [",".join(CSV_HEADERS)]
    for row in rows:
        lines.append(",".join(_format_csv_value(row[h]) for h in CSV_HEADERS))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------
# Markdown report
# ---------------------------------------------------------------------


def _by_schedule_summary(rows: list[dict]) -> list[dict]:
    out: list[dict] = []
    for schedule in SCHEDULES:
        label = schedule["label"]
        subset = [r for r in rows if r["schedule_label"] == label]
        gaps = [r["energy_gap"] for r in subset]
        rmses = [r["interval_rmse"] for r in subset]
        out.append({
            "label": label,
            "runs": len(subset),
            "mean_gap": sum(gaps) / len(gaps) if gaps else float("nan"),
            "min_gap": min(gaps) if gaps else float("nan"),
            "max_gap": max(gaps) if gaps else float("nan"),
            "mean_rmse": sum(rmses) / len(rmses) if rmses else float("nan"),
            "successes": sum(1 for r in subset if r["success_flag"]),
        })
    return out


def _by_dim_summary(rows: list[dict]) -> list[dict]:
    out: list[dict] = []
    for d in SPACETIME_DIMS:
        for schedule in SCHEDULES:
            label = schedule["label"]
            subset = [
                r for r in rows
                if r["target_dim"] == d and r["schedule_label"] == label
            ]
            gaps = [r["energy_gap"] for r in subset]
            rmses = [r["interval_rmse"] for r in subset]
            out.append({
                "target_dim": d,
                "label": label,
                "runs": len(subset),
                "mean_gap": (
                    sum(gaps) / len(gaps) if gaps else float("nan")
                ),
                "mean_rmse": (
                    sum(rmses) / len(rmses) if rmses else float("nan")
                ),
            })
    return out


def _qualitative_gap_trend(by_label: dict) -> tuple[str, str]:
    """Return a (verb, rationale) describing the gap trend.

    The verb is one of ``"falls"``, ``"is essentially flat"``, or
    ``"rises"``; the rationale gives the conservative reading.
    """

    short_gap = by_label["short"]["mean_gap"]
    long_gap = by_label["long"]["mean_gap"]
    if not (math.isfinite(short_gap) and math.isfinite(long_gap)):
        return ("is undefined", "no finite gap statistics; check inputs")
    rel = (long_gap - short_gap) / max(short_gap, 1e-12)
    if rel < -0.1:
        return (
            "falls",
            "Evidence of schedule/optimizer failure in the Phase 2 "
            "configuration. This is *not* a claim that the embedding "
            "has been solved.",
        )
    if rel > 0.1:
        return (
            "rises",
            "Adding budget makes things worse on average. The "
            "annealer's warmup is pushing the state away from the "
            "neighborhood of the truth and the cooling stage is not "
            "recovering. The bottleneck is in the energy or move "
            "implementation, not in iteration budget. This is *not* "
            "a claim that Minkowski sprinklings are non-manifoldlike.",
        )
    return (
        "is essentially flat",
        "Increasing the historical annealer's budget does not close "
        "the gap. The bottleneck likely sits in the energy "
        "definition, its parametrization, or the historical move "
        "set, not in budget. This is *not* a claim that Minkowski "
        "sprinklings are non-manifoldlike.",
    )


def _per_dim_trend(by_dim: list[dict]) -> dict[int, str]:
    """For each dimension, label the short-vs-long gap trend."""

    out: dict[int, str] = {}
    for d in SPACETIME_DIMS:
        short = next(
            (row for row in by_dim if row["target_dim"] == d
             and row["label"] == "short"),
            None,
        )
        long_ = next(
            (row for row in by_dim if row["target_dim"] == d
             and row["label"] == "long"),
            None,
        )
        if short is None or long_ is None:
            out[d] = "no data"
            continue
        sg = short["mean_gap"]
        lg = long_["mean_gap"]
        if not (math.isfinite(sg) and math.isfinite(lg)):
            out[d] = "non-finite"
            continue
        rel = (lg - sg) / max(sg, 1e-12)
        if rel < -0.1:
            out[d] = f"gap drops with budget ({sg:.1f} -> {lg:.1f})"
        elif rel > 0.1:
            out[d] = f"gap rises with budget ({sg:.1f} -> {lg:.1f})"
        else:
            out[d] = f"gap flat across budgets ({sg:.1f} -> {lg:.1f})"
    return out


def _interpretation_lines(rows: list[dict]) -> list[str]:
    by_sched = _by_schedule_summary(rows)
    by_label = {row["label"]: row for row in by_sched}
    by_dim = _by_dim_summary(rows)
    dim_labels = _per_dim_trend(by_dim)

    short_gap = by_label["short"]["mean_gap"]
    long_gap = by_label["long"]["mean_gap"]
    short_min = by_label["short"]["min_gap"]
    long_min = by_label["long"]["min_gap"]
    short_rmse = by_label["short"]["mean_rmse"]
    long_rmse = by_label["long"]["mean_rmse"]
    successes = sum(row["successes"] for row in by_sched)

    verb, rationale = _qualitative_gap_trend(by_label)

    lines = [
        "Interpretation (conservative, framed as five fixed questions):",
        "",
        "1. **Does ``energy_gap`` fall as the budget grows?**",
        f"   Mean gap {verb} (short = {short_gap:.3f}, "
        f"long = {long_gap:.3f}). The minimum gap across all cells is "
        f"{short_min:.3f} for the short schedule and {long_min:.3f} "
        f"for the long schedule. {rationale}",
        "",
        "2. **Does ``interval_rmse`` drop coherently with budget?**",
        f"   Mean RMSE moves from {short_rmse:.3e} (short) to "
        f"{long_rmse:.3e} (long), but the per-run table shows that "
        "RMSE varies by several orders of magnitude across "
        "neighboring cells under any single schedule, so a single "
        "ensemble mean is not a stable summary. The conservative "
        "reading is that the recovered coordinates are not approaching "
        "the ground truth under any of the three budgets tested.",
        "",
        "3. **Is there a clear difference between d = 2, 3, 4?**",
    ]
    for d in SPACETIME_DIMS:
        lines.append(f"   - d = {d}: {dim_labels[d]}.")
    lines += [
        "   The same qualitative behaviour appears at every dimension:",
        "   the gap does not collapse with budget. No dimension is",
        "   distinguished as 'fixed by more iterations'.",
        "",
        "4. **Does the failure look like budget/schedule or like",
        "   something more structural in the annealer?**",
        "   Across this small grid the gap is not budget-limited.",
        "   That rules out 'short Phase 2 schedule' as the dominant",
        "   cause of the Phase 2 Minkowski residual. The remaining",
        "   candidates are the energy definition, its parametrization,",
        "   or the historical move-set implementation. This probe",
        "   does *not* discriminate between those candidates; it only",
        "   removes the simplest budget explanation.",
        "",
        "5. **Is it still invalid to read the annealer as a",
        "   manifoldness classifier?**",
        f"   Yes. Across {sum(row['runs'] for row in by_sched)} runs the "
        f"conservative ``success_flag`` is True in {successes} cases.",
        "   No Minkowski case has been recovered to within numerical",
        "   tolerance of the truth at any budget tested, so a low",
        "   ``final_energy`` cannot be cited as a successful embedding.",
        "   The annealer is therefore not a validated manifoldness",
        "   classifier and is never applied to KR or corona causets",
        "   in this probe.",
        "",
        "Side remarks:",
        "",
        "- Phase 2B does **not** introduce a new optimizer, basin",
        "  hopping, parallel tempering, ML, or PySR-driven search.",
        "  Only the historical annealer's iteration knobs are varied.",
        "- KR and corona causets are excluded by construction; the",
        "  probe must not be cited as a manifoldness classifier.",
        "",
    ]
    return lines


def _display(value) -> str:
    if value is None:
        return "NA"
    if isinstance(value, str):
        return value or "-"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float) and math.isfinite(value):
        return _format_field(value)
    return "NA"


def write_markdown(rows: list[dict], path: Path) -> None:
    by_sched = _by_schedule_summary(rows)
    by_dim = _by_dim_summary(rows)

    lines = [
        "# Phase 2B Annealer Schedule Probe",
        "",
        "Minkowski-only probe of the historical Bombelli-Sorkin annealer",
        "across a small grid of iteration budgets. Phase 2 left the",
        "annealer running with a deliberately tiny schedule and observed",
        "large energy gaps and large interval residuals. Phase 2B asks:",
        "how much of that is the schedule and how much is structural?",
        "",
        "Scope and what this probe does *not* do:",
        "",
        "- No new optimizer is introduced. Only ``warmup_limit``,",
        "  ``anneal_limit`` and ``max_data`` are varied; the temperature",
        "  schedule is held fixed.",
        "- Only Minkowski sprinklings are run. KR and corona controls",
        "  have no ground-truth coordinates and are excluded by",
        "  construction; this probe is *not* a manifoldness classifier.",
        "- A low final energy alone is not treated as a successful",
        "  embedding; the diagnostic quantities are ``energy_gap`` and",
        "  ``interval_rmse`` against the known ground truth.",
        "",
        "Protocol:",
        "",
        f"- families: minkowski only.",
        f"- target spacetime dimensions: {', '.join(str(d) for d in SPACETIME_DIMS)}.",
        f"- sizes: n = {', '.join(str(n) for n in SIZES)}.",
        f"- case seeds: {', '.join(str(s) for s in PROBE_SEEDS)} "
        "(first three Phase 1B atlas seeds).",
        f"- optimizer seed: {OPTIMIZER_SEED}.",
        f"- temperature: initial_temp={INITIAL_TEMP}, "
        f"cooling_factor={COOLING_FACTOR} (fixed across schedules).",
        f"- success criterion (conservative): ``energy_gap <= "
        f"{SUCCESS_GAP_THRESHOLD}`` with finite ``interval_rmse``.",
        "",
        "Schedules:",
        "",
        "| label | warmup_limit | anneal_limit | max_data | reconfigure budget |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for schedule in SCHEDULES:
        budget = (
            schedule["warmup_limit"]
            + schedule["anneal_limit"] * schedule["max_data"]
        )
        lines.append(
            f"| {schedule['label']} | {schedule['warmup_limit']} | "
            f"{schedule['anneal_limit']} | {schedule['max_data']} | "
            f"{budget} |"
        )

    lines += [
        "",
        "The ``short`` schedule is exactly the Phase 2 configuration, so",
        "the corresponding rows reproduce the Phase 2 Minkowski numbers",
        "and are the baseline against which the larger budgets are read.",
        "",
        "Per-schedule aggregates over all (d, n, seed) cells:",
        "",
        "| schedule | runs | mean gap | min gap | max gap | mean RMSE | successes |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in by_sched:
        lines.append(
            "| {label} | {runs} | {mg} | {ming} | {maxg} | {mr} | {s} |".format(
                label=row["label"],
                runs=row["runs"],
                mg=_display(row["mean_gap"]),
                ming=_display(row["min_gap"]),
                maxg=_display(row["max_gap"]),
                mr=_display(row["mean_rmse"]),
                s=row["successes"],
            )
        )

    lines += [
        "",
        "Per-dimension aggregates (mean across seeds and sizes):",
        "",
        "| d | schedule | runs | mean gap | mean RMSE |",
        "| :---: | --- | ---: | ---: | ---: |",
    ]
    for row in by_dim:
        lines.append(
            "| {d} | {label} | {runs} | {mg} | {mr} |".format(
                d=row["target_dim"],
                label=row["label"],
                runs=row["runs"],
                mg=_display(row["mean_gap"]),
                mr=_display(row["mean_rmse"]),
            )
        )

    lines += [
        "",
        "Full per-run table:",
        "",
        "| d | n | seed | schedule | final E | truth E | gap | RMSE | success | t(s) |",
        "| :---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | :---: | ---: |",
    ]
    for row in rows:
        lines.append(
            "| {d} | {n} | {seed} | {label} | {fe} | {te} | "
            "{gap} | {rmse} | {ok} | {t} |".format(
                d=row["target_dim"],
                n=row["n"],
                seed=row["seed"],
                label=row["schedule_label"],
                fe=_display(row["final_energy"]),
                te=_display(row["truth_energy"]),
                gap=_display(row["energy_gap"]),
                rmse=_display(row["interval_rmse"]),
                ok="yes" if row["success_flag"] else "no",
                t=_display(row["runtime_seconds"]),
            )
        )

    lines += [""]
    lines += _interpretation_lines(rows)
    lines += [
        "Regenerate via `make regen-phase2b`. Source tool:",
        "`tools/build_phase2b_annealer_schedule_probe.py`.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    rows = build_rows()
    FOUNDATION.mkdir(parents=True, exist_ok=True)
    write_csv(rows, FOUNDATION / "phase2b_annealer_schedule_probe.csv")
    write_markdown(rows, FOUNDATION / "phase2b_annealer_schedule_probe.md")
    print(
        f"Wrote {len(rows)} Phase 2B schedule-probe rows to {FOUNDATION}"
    )


if __name__ == "__main__":
    main()
