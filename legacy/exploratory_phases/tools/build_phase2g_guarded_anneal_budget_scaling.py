#!/usr/bin/env python3
"""Build the Phase 2G guarded-anneal budget-scaling probe.

Phase 2F showed GUARDED_WARMUP_FIXES_PRIMARY_FAILURE: replacing the
unconditional legacy warmup with an energy-gated warmup (GUARD_THRESHOLD=0)
completely fixes small-noise destruction and improves random_init. Medium-noise
and random_init remain unresolved at the short schedule (anneal_limit=10).

Phase 2G asks the follow-up question: once the warmup failure is corrected,
does increasing the anneal budget recover the residual failures?

The physical question: is the residual failure budget-limited (more annealing
closes the gap) or move-set/cooling-limited (more budget gives no traction)?

Protocol:
  - warmup_mode: guarded_warmup only (GUARD_THRESHOLD = 0.0, 10 warmup steps)
  - anneal schedule varies across four labels: short/medium/long/xlong
  - max_data = 4 (fixed), T0 = 100, gamma = 0.9 (fixed)
  - grid: d in {2,3,4}, n in {32,64}, seeds {1959,1962,1987}
  - init_labels: truth, small_noise, medium_noise, random_init

OverflowError handling: at high budgets from random_init, coordinates can
grow to extreme values causing floating-point overflow. These rows are
recorded with final_energy = NaN and notes documenting the overflow.
This is itself a diagnostic finding about coordinate-space stability.
"""

from __future__ import annotations

import contextlib
import io
import math
import random
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import cones  # noqa: E402
from tools.build_phase1_atlas import _format_field  # noqa: E402
from tools.build_phase1b_scaling_atlas import FOUNDATION, SEEDS  # noqa: E402
import validation_suite as vs  # noqa: E402


SPACETIME_DIMS = (2, 3, 4)
SIZES = (32, 64)
PROBE_SEEDS = SEEDS[:3]  # (1959, 1962, 1987)
OPTIMIZER_SEED = 1987

GUARDED_WARMUP_LIMIT = 10
GUARD_THRESHOLD = 0.0
MAX_DATA = 4
INITIAL_TEMP = 100.0
COOLING_FACTOR = 0.9

NOISE_SMALL = 1e-3
NOISE_MEDIUM = 5e-2

INIT_LABELS: tuple[tuple[str, float | None], ...] = (
    ("truth", 0.0),
    ("truth_plus_small_noise", NOISE_SMALL),
    ("truth_plus_medium_noise", NOISE_MEDIUM),
    ("random_init", None),
)

# Anneal budget sweep. warmup_limit is fixed; only anneal_limit varies.
SCHEDULES: tuple[tuple[str, int], ...] = (
    ("short", 10),
    ("medium", 30),
    ("long", 100),
    ("xlong", 300),
)

# Recovery thresholds. These are conservative and documented, not physics.
SMALL_NOISE_ENERGY_THRESHOLD = 0.10     # final_energy < 0.1 for small_noise
RECOVERY_ENERGY_THRESHOLD = 1.0         # final_energy < 1.0 for medium/random
RECOVERY_RELATIVE_THRESHOLD = 0.90     # or (initial-final)/initial > 90%

CSV_HEADERS = (
    "family",
    "target_dim",
    "n",
    "seed",
    "init_label",
    "schedule_label",
    "anneal_limit",
    "noise_epsilon",
    "paired_key",
    "initial_energy",
    "post_warmup_energy",
    "final_energy",
    "delta_energy",
    "initial_interval_rmse",
    "final_interval_rmse",
    "initial_distance_to_truth_rms",
    "final_distance_to_truth_rms",
    "improved_energy",
    "improved_interval_rmse",
    "preserved_near_truth",
    "recovered_flag",
    "warmup_attempted_moves",
    "warmup_accepted_moves",
    "warmup_rejected_moves",
    "runtime_seconds",
    "notes",
)


# ------------------------------------------------------------------
# Coordinate helpers
# ------------------------------------------------------------------


def _coord_distance_rms(
    coords: list[vs.Coord],
    truth_scaled: list[vs.Coord],
    n: int,
    d_spatial: int,
) -> float:
    sq_sum = 0.0
    for i in range(n):
        ci, ti = coords[i], truth_scaled[i]
        dr = ci[0] - ti[0]
        sq_sum += dr * dr
        for k in range(1, d_spatial + 1):
            dx = (ci[k] if k < len(ci) else 0.0) - (ti[k] if k < len(ti) else 0.0)
            sq_sum += dx * dx
    return math.sqrt(sq_sum / n)


def _coords_from_sim(sim: cones.ConesSimulator) -> list[vs.Coord]:
    return [(sim.rold[i], *sim.xold[i]) for i in range(sim.n)]


def _scale_truth(points: list[vs.Coord], rave_target: float, n: int) -> list[vs.Coord]:
    rave_truth = sum(p[0] for p in points) / n
    if rave_truth <= 0.0:
        return list(points)
    scale = rave_target / rave_truth
    return [tuple(c * scale for c in p) for p in points]


def _paired_key(d: int, n: int, seed: int, label: str, eps) -> str:
    return f"{d}|{n}|{seed}|{label}|{'NA' if eps is None else eps}"


# ------------------------------------------------------------------
# Custom startup
# ------------------------------------------------------------------


def _custom_startup(
    sim: cones.ConesSimulator,
    points: list[vs.Coord],
    noise_epsilon: float,
    noise_rng: random.Random,
) -> None:
    n = sim.n
    for i in range(n):
        sim.change[i] = True
        t_i = points[i][0]
        if noise_epsilon > 0.0:
            t_i = max(1e-12, t_i + noise_epsilon * noise_rng.gauss(0, 1))
        sim.rnew[i] = t_i
        for k in range(sim.dim):
            x_ik = points[i][k + 1] if k + 1 < len(points[i]) else 0.0
            if noise_epsilon > 0.0:
                x_ik += noise_epsilon * noise_rng.gauss(0, 1)
            sim.xnew[i][k] = x_ik
    sim.rave = sum(sim.rnew) / n
    sim.energy()
    sim.update()
    sim.initial_energy = sim.energies[0]


# ------------------------------------------------------------------
# Guarded warmup
# ------------------------------------------------------------------


def _guarded_warmup(sim: cones.ConesSimulator) -> dict:
    """GUARD_THRESHOLD=0: accept iff deltae <= 0 (pre-normalization)."""
    energy_before = sim.energies[0]
    attempted = accepted = rejected = 0
    for _ in range(GUARDED_WARMUP_LIMIT):
        if sim.energies[0] <= 0.0:
            break
        for i in range(sim.n):
            sim.change[i] = False
        rave_saved = sim.r
        attempted += 1
        sim.reconfigure()
        sim.energy()
        if sim.deltae <= GUARD_THRESHOLD:
            sim.update()
            accepted += 1
        else:
            for i in range(sim.n):
                sim.change[i] = False
            sim.rave = rave_saved
            sim.deltae = 0.0
            rejected += 1
    sim.statistics()
    sim.warmup_energy = sim.energies[0]
    return {
        "warmup_attempted_moves": attempted,
        "warmup_accepted_moves": accepted,
        "warmup_rejected_moves": rejected,
        "post_warmup_energy": sim.energies[0],
    }


# ------------------------------------------------------------------
# Recovery flag
# ------------------------------------------------------------------


def _recovered(init_label: str, initial_energy: float, final_energy: float) -> bool:
    if not math.isfinite(final_energy):
        return False
    if init_label == "truth":
        return final_energy < 1e-6
    if init_label == "truth_plus_small_noise":
        return final_energy < SMALL_NOISE_ENERGY_THRESHOLD
    rel_reduction = (
        (initial_energy - final_energy) / initial_energy
        if initial_energy > 1e-9 else 0.0
    )
    return (
        final_energy < RECOVERY_ENERGY_THRESHOLD
        or rel_reduction > RECOVERY_RELATIVE_THRESHOLD
    )


# ------------------------------------------------------------------
# Per-row computation
# ------------------------------------------------------------------


def _probe_row(
    d_spacetime: int,
    n: int,
    seed: int,
    init_label: str,
    noise_epsilon: float | None,
    schedule_label: str,
    anneal_limit: int,
) -> dict:
    d_spatial = d_spacetime - 1
    matrix, points = vs.sprinkle_minkowski_diamond(
        n=n, seed=seed, d_spacetime=d_spacetime
    )
    noise_seed = seed * 10007 + d_spacetime * 1009 + n * 97
    noise_rng = random.Random(noise_seed)

    overflow = False
    final_energy = float("nan")
    final_interval_rmse = float("nan")
    final_dist = float("nan")
    wmeta: dict = {}

    t0 = time.perf_counter()
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            sim = cones.ConesSimulator(
                z=matrix,
                dim=d_spatial,
                seed=OPTIMIZER_SEED,
                interactive=False,
                max_data=MAX_DATA,
                plot_path=None,
                warmup_limit=GUARDED_WARMUP_LIMIT,
                anneal_limit=anneal_limit,
                initial_temp=INITIAL_TEMP,
                cooling_factor=COOLING_FACTOR,
                backend="cpu",
            )
            buf = io.StringIO()

            if noise_epsilon is None:
                sim.startup(buf)
            else:
                _custom_startup(sim, points, noise_epsilon, noise_rng)

            initial_energy = sim.initial_energy
            initial_coords = _coords_from_sim(sim)
            rave_initial = sum(c[0] for c in initial_coords) / n
            truth_scaled_initial = _scale_truth(points, rave_initial, n)
            initial_interval_rmse = vs.interval_rmse(initial_coords, truth_scaled_initial)
            initial_dist = _coord_distance_rms(
                initial_coords, truth_scaled_initial, n, d_spatial
            )

            wmeta = _guarded_warmup(sim)
            sim.anneal(buf)

            final_energy = sim.data[-1][1] if sim.data else sim.eave
            final_coords = _coords_from_sim(sim)
            rave_final = sum(c[0] for c in final_coords) / n
            truth_scaled_final = _scale_truth(points, rave_final, n)
            final_interval_rmse = vs.interval_rmse(final_coords, truth_scaled_final)
            final_dist = _coord_distance_rms(
                final_coords, truth_scaled_final, n, d_spatial
            )

    except OverflowError:
        overflow = True
        if not wmeta:
            initial_energy = float("nan")
            initial_interval_rmse = float("nan")
            initial_dist = float("nan")

    runtime = time.perf_counter() - t0

    if overflow:
        delta_energy = float("nan")
        improved_energy = False
        improved_interval_rmse = False
        preserved_near_truth = False
        recovered = False
        note = (
            f"OverflowError: coordinate runaway during anneal "
            f"(d={d_spacetime}, n={n}, seed={seed}, "
            f"{schedule_label} budget, {init_label})"
        )
    else:
        delta_energy = final_energy - initial_energy
        improved_energy = delta_energy < -1e-9
        improved_interval_rmse = (
            math.isfinite(final_interval_rmse)
            and math.isfinite(initial_interval_rmse)
            and final_interval_rmse < initial_interval_rmse - 1e-9
        )
        preserved_near_truth = math.isfinite(delta_energy) and delta_energy <= 1e-6
        recovered = _recovered(init_label, initial_energy, final_energy)

        wword = "guarded+anneal"
        if init_label == "truth":
            note = (
                f"truth preserved: energy remains zero ({schedule_label})"
                if abs(final_energy) < 1e-9
                else f"truth not preserved ({schedule_label}): final E={final_energy:.4g}"
            )
        elif init_label in ("truth_plus_small_noise", "truth_plus_medium_noise"):
            if improved_energy:
                note = f"{wword} improved {initial_energy:.4g} -> {final_energy:.4g} ({schedule_label})"
            elif preserved_near_truth:
                note = f"{wword} preserved near-truth ({schedule_label}): {initial_energy:.4g} -> {final_energy:.4g}"
            else:
                note = f"{wword} lost near-truth ({schedule_label}): {initial_energy:.4g} -> {final_energy:.4g}"
        else:
            if improved_energy:
                note = f"{wword} improved random_init ({schedule_label}): {initial_energy:.4g} -> {final_energy:.4g}"
            else:
                note = f"{wword} did not improve random_init ({schedule_label}): {initial_energy:.4g} -> {final_energy:.4g}"

    return {
        "family": "minkowski",
        "target_dim": d_spacetime,
        "n": n,
        "seed": seed,
        "init_label": init_label,
        "schedule_label": schedule_label,
        "anneal_limit": anneal_limit,
        "noise_epsilon": noise_epsilon,
        "paired_key": _paired_key(d_spacetime, n, seed, init_label, noise_epsilon),
        "initial_energy": initial_energy,
        "post_warmup_energy": wmeta.get("post_warmup_energy", float("nan")),
        "final_energy": final_energy,
        "delta_energy": delta_energy,
        "initial_interval_rmse": initial_interval_rmse,
        "final_interval_rmse": final_interval_rmse,
        "initial_distance_to_truth_rms": initial_dist,
        "final_distance_to_truth_rms": final_dist,
        "improved_energy": improved_energy,
        "improved_interval_rmse": improved_interval_rmse,
        "preserved_near_truth": preserved_near_truth,
        "recovered_flag": recovered,
        "warmup_attempted_moves": wmeta.get("warmup_attempted_moves", 0),
        "warmup_accepted_moves": wmeta.get("warmup_accepted_moves", 0),
        "warmup_rejected_moves": wmeta.get("warmup_rejected_moves", 0),
        "runtime_seconds": runtime,
        "notes": note,
    }


def build_rows() -> list[dict]:
    rows: list[dict] = []
    for d in SPACETIME_DIMS:
        for n in SIZES:
            for seed in PROBE_SEEDS:
                for label, eps in INIT_LABELS:
                    for sched_label, anneal_limit in SCHEDULES:
                        rows.append(
                            _probe_row(d, n, seed, label, eps, sched_label, anneal_limit)
                        )
    return rows


# ------------------------------------------------------------------
# CSV output
# ------------------------------------------------------------------


def _fmt(value) -> str:
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float) and math.isnan(value):
        return "NaN"
    if isinstance(value, str):
        return value
    if isinstance(value, int):
        return str(value)
    return _format_field(value)


def write_csv(rows: list[dict], path: Path) -> None:
    lines = [",".join(CSV_HEADERS)]
    for row in rows:
        lines.append(",".join(_fmt(row[h]) for h in CSV_HEADERS))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ------------------------------------------------------------------
# Aggregate helpers
# ------------------------------------------------------------------


def _agg(rows: list[dict], label: str, sched: str) -> dict | None:
    subset = [
        r for r in rows
        if r["init_label"] == label and r["schedule_label"] == sched
        and math.isfinite(float(r["final_energy"])) if r["final_energy"] != "NaN" else False
    ]
    # Re-filter cleanly
    subset = [
        r for r in rows
        if r["init_label"] == label
        and r["schedule_label"] == sched
        and isinstance(r["final_energy"], float)
        and math.isfinite(r["final_energy"])
    ]
    if not subset:
        return None
    fe = [r["final_energy"] for r in subset]
    pnt = sum(1 for r in subset if r["preserved_near_truth"])
    rec = sum(1 for r in subset if r["recovered_flag"])
    total = len([
        r for r in rows
        if r["init_label"] == label and r["schedule_label"] == sched
    ])
    return {
        "valid": len(subset),
        "total": total,
        "mean_final_energy": sum(fe) / len(fe),
        "preserved_count": pnt,
        "recovered_count": rec,
    }


# ------------------------------------------------------------------
# Verdict
# ------------------------------------------------------------------


def _verdict(rows: list[dict]) -> tuple[str, str]:
    # Check if medium_noise improves monotonically with budget
    sched_order = [s for s, _ in SCHEDULES]

    def mean_fe(label: str, sched: str) -> float:
        subset = [
            r for r in rows
            if r["init_label"] == label
            and r["schedule_label"] == sched
            and isinstance(r["final_energy"], float)
            and math.isfinite(r["final_energy"])
        ]
        return sum(r["final_energy"] for r in subset) / len(subset) if subset else float("inf")

    def mean_rec(label: str, sched: str) -> float:
        subset = [
            r for r in rows
            if r["init_label"] == label and r["schedule_label"] == sched
        ]
        total = len(subset)
        rec = sum(1 for r in subset if isinstance(r["recovered_flag"], bool) and r["recovered_flag"])
        return rec / total if total else 0.0

    # Is there a clear improvement from short to long?
    med_short = mean_fe("truth_plus_medium_noise", "short")
    med_long = mean_fe("truth_plus_medium_noise", "long")
    ri_short = mean_fe("random_init", "short")
    ri_long = mean_fe("random_init", "long")

    med_improves = med_long < med_short * 0.7  # >30% reduction
    ri_improves = ri_long < ri_short * 0.7

    # Check monotonicity: does energy decrease with each budget step?
    med_energies = [mean_fe("truth_plus_medium_noise", s) for s in sched_order]
    ri_energies = [mean_fe("random_init", s) for s in sched_order]
    med_monotone = all(
        med_energies[i] >= med_energies[i + 1]
        for i in range(len(med_energies) - 1)
        if math.isfinite(med_energies[i]) and math.isfinite(med_energies[i + 1])
    )
    ri_monotone = all(
        ri_energies[i] >= ri_energies[i + 1]
        for i in range(len(ri_energies) - 1)
        if math.isfinite(ri_energies[i]) and math.isfinite(ri_energies[i + 1])
    )

    # Check d/n breakdown: does d=2/n=32 do better than d=4/n=64?
    def mean_fe_dn(label: str, d: int, n: int, sched: str) -> float:
        subset = [
            r for r in rows
            if r["init_label"] == label
            and r["target_dim"] == d and r["n"] == n
            and r["schedule_label"] == sched
            and isinstance(r["final_energy"], float)
            and math.isfinite(r["final_energy"])
        ]
        return sum(r["final_energy"] for r in subset) / len(subset) if subset else float("inf")

    d2n32_med_long = mean_fe_dn("truth_plus_medium_noise", 2, 32, "long")
    d4n64_med_long = mean_fe_dn("truth_plus_medium_noise", 4, 64, "long")
    dimension_barrier = (
        math.isfinite(d2n32_med_long) and math.isfinite(d4n64_med_long)
        and d4n64_med_long > 3 * d2n32_med_long
    )

    if med_improves and ri_improves and med_monotone and ri_monotone:
        return (
            "BUDGET_FIXES_RESIDUAL",
            (
                f"Medium-noise and random-init final energies decrease clearly "
                f"and monotonically with anneal budget after guarded_warmup "
                f"(medium_noise short={med_short:.1f} → long={med_long:.1f}; "
                f"random_init short={ri_short:.1f} → long={ri_long:.1f}). "
                "Increasing the anneal budget opens meaningful traction. "
                "The residual failure from Phase 2F was primarily budget-limited, "
                "not move-set or cooling limited."
            ),
        )

    if (med_improves or ri_improves) and dimension_barrier:
        return (
            "DIMENSION_OR_SIZE_BARRIER",
            (
                f"Some improvement with budget (medium_noise short={med_short:.1f} → "
                f"long={med_long:.1f}; random_init short={ri_short:.1f} → long={ri_long:.1f}), "
                f"but d=2/n=32 cells recover far better than d=4/n=64 "
                f"({d2n32_med_long:.1f} vs {d4n64_med_long:.1f} at long). "
                "Budget helps in lower dimensions, but higher dimensions or "
                "larger sizes hit a structural barrier."
            ),
        )

    if med_improves or ri_improves:
        return (
            "BUDGET_PARTIALLY_HELPS",
            (
                f"Partial improvement with budget: "
                f"medium_noise short={med_short:.1f} → long={med_long:.1f} "
                f"({'monotone' if med_monotone else 'non-monotone'}); "
                f"random_init short={ri_short:.1f} → long={ri_long:.1f} "
                f"({'monotone' if ri_monotone else 'non-monotone'}). "
                "Budget provides some traction after guarded_warmup but not "
                "robust recovery. Move-set or cooling may be a secondary limit."
            ),
        )

    return (
        "MOVE_SET_OR_COOLING_LIMIT",
        (
            f"Increasing anneal budget after guarded_warmup does not clearly "
            f"reduce final energy for medium-noise or random-init. "
            f"medium_noise short={med_short:.1f}, long={med_long:.1f}; "
            f"random_init short={ri_short:.1f}, long={ri_long:.1f}. "
            "The residual failure is not budget-limited. The move set or "
            "cooling schedule is insufficient to escape the current basins, "
            "independent of the warmup correction."
        ),
    )


# ------------------------------------------------------------------
# Markdown report
# ------------------------------------------------------------------


def write_markdown(rows: list[dict], path: Path) -> None:
    verdict_label, verdict_text = _verdict(rows)
    overflow_count = sum(1 for r in rows if not isinstance(r["final_energy"], float) or not math.isfinite(r["final_energy"]))

    lines = [
        "# Phase 2G Guarded-Anneal Budget-Scaling Probe",
        "",
        "Budget-scaling study with ``guarded_warmup`` fixed: does increasing",
        "the anneal budget after the corrected warmup recover medium-noise",
        "and random-init configurations?",
        "",
        "## Verdict",
        "",
        f"**{verdict_label}**",
        "",
        verdict_text,
        "",
        "## Protocol",
        "",
        "- families: minkowski only.",
        f"- d ∈ {{2, 3, 4}}, n ∈ {{32, 64}}, "
        f"seeds {', '.join(str(s) for s in PROBE_SEEDS)}.",
        f"- optimizer seed: {OPTIMIZER_SEED}.",
        f"- warmup: guarded_warmup, GUARD_THRESHOLD={GUARD_THRESHOLD}, "
        f"warmup_limit={GUARDED_WARMUP_LIMIT}.",
        f"- NOISE_SMALL = {NOISE_SMALL}, NOISE_MEDIUM = {NOISE_MEDIUM}.",
        f"- max_data={MAX_DATA}, T₀={INITIAL_TEMP}, γ={COOLING_FACTOR} (fixed).",
        "- ``paired_key`` = ``target_dim|n|seed|init_label|noise_epsilon``.",
        "- ``preserved_near_truth``: ``delta_energy ≤ 1e-6``.",
        f"- ``recovered_flag`` (documented thresholds, not physics):",
        f"  truth: final_energy < 1e-6;",
        f"  small_noise: final_energy < {SMALL_NOISE_ENERGY_THRESHOLD};",
        f"  medium_noise/random_init: final_energy < {RECOVERY_ENERGY_THRESHOLD}",
        f"  OR relative_reduction > {int(RECOVERY_RELATIVE_THRESHOLD * 100)}%.",
        f"- OverflowError: {overflow_count} rows hit coordinate overflow at high budget.",
        "  Recorded as NaN; this is a diagnostic, not a crash.",
        "",
        "## Anneal schedules",
        "",
        "| label | anneal_limit | budget (warmup+anneal×max_data) |",
        "| --- | ---: | ---: |",
    ]
    for slabel, alimit in SCHEDULES:
        budget = GUARDED_WARMUP_LIMIT + alimit * MAX_DATA
        lines.append(f"| {slabel} | {alimit} | {budget} |")

    lines += [
        "",
        "## Per-label aggregate by schedule",
        "",
        "| init | schedule | valid/total | mean final E | preserved | recovered |",
        "| --- | --- | ---: | ---: | ---: | ---: |",
    ]
    for label, _ in INIT_LABELS:
        for slabel, _ in SCHEDULES:
            b = _agg(rows, label, slabel)
            if b is None:
                lines.append(f"| {label} | {slabel} | 0/? | NaN | — | — |")
                continue
            lines.append(
                "| {label} | {sched} | {v}/{t} | {fe} | {pnt}/{v} | {rec}/{v} |".format(
                    label=label,
                    sched=slabel,
                    v=b["valid"],
                    t=b["total"],
                    fe=_format_field(b["mean_final_energy"]),
                    pnt=b["preserved_count"],
                    rec=b["recovered_count"],
                )
            )

    # Per-label budget curves
    lines += ["", "## Energy vs budget curves (mean over valid rows)"]
    for label, _ in INIT_LABELS:
        lines += ["", f"### {label}"]
        lines += [
            "| schedule | anneal_limit | mean final E | preserved | recovered |",
            "| --- | ---: | ---: | ---: | ---: |",
        ]
        for slabel, alimit in SCHEDULES:
            b = _agg(rows, label, slabel)
            total = len([r for r in rows if r["init_label"] == label and r["schedule_label"] == slabel])
            if b is None:
                lines.append(f"| {slabel} | {alimit} | NaN | —/{total} | —/{total} |")
            else:
                lines.append(
                    "| {s} | {a} | {fe} | {pnt}/{v} | {rec}/{v} |".format(
                        s=slabel,
                        a=alimit,
                        fe=_format_field(b["mean_final_energy"]),
                        pnt=b["preserved_count"],
                        v=b["valid"],
                        rec=b["recovered_count"],
                    )
                )

    lines += ["", "## Seven fixed questions"]
    sched_order = [s for s, _ in SCHEDULES]

    def mean_fe_label(label, sched):
        subset = [
            r for r in rows
            if r["init_label"] == label and r["schedule_label"] == sched
            and isinstance(r["final_energy"], float) and math.isfinite(r["final_energy"])
        ]
        return sum(r["final_energy"] for r in subset) / len(subset) if subset else float("nan")

    lines += ["", "1. **Does small_noise stay resolved as anneal budget increases?**"]
    for s in sched_order:
        mfe = mean_fe_label("truth_plus_small_noise", s)
        lines.append(f"   {s}: mean final E = {_format_field(mfe)}.")

    lines += ["", "2. **Does medium_noise start recovering with more anneal?**"]
    for s in sched_order:
        mfe = mean_fe_label("truth_plus_medium_noise", s)
        lines.append(f"   {s}: mean final E = {_format_field(mfe)}.")

    lines += ["", "3. **Does random_init improve monotonically with budget?**"]
    ri_vals = []
    for s in sched_order:
        mfe = mean_fe_label("random_init", s)
        ri_vals.append(mfe)
        lines.append(f"   {s}: mean final E = {_format_field(mfe)}.")
    finite_vals = [v for v in ri_vals if math.isfinite(v)]
    if len(finite_vals) >= 2:
        monotone = all(finite_vals[i] >= finite_vals[i + 1] for i in range(len(finite_vals) - 1))
        lines.append(f"   Monotone decreasing: {'yes' if monotone else 'no'}.")

    lines += ["", "4. **Is the energy trajectory coherent or erratic?**"]
    for label, _ in INIT_LABELS:
        vals = [mean_fe_label(label, s) for s in sched_order]
        finite = [v for v in vals if math.isfinite(v)]
        if len(finite) >= 2:
            diffs = [finite[i + 1] - finite[i] for i in range(len(finite) - 1)]
            n_dec = sum(1 for d in diffs if d < 0)
            n_inc = sum(1 for d in diffs if d > 0)
            lines.append(
                f"   {label}: {n_dec}/{len(diffs)} steps decreasing, "
                f"{n_inc}/{len(diffs)} increasing."
            )

    lines += ["", "5. **Does the residual failure point to budget, cooling, move-set, or global init?**"]
    lines.append(f"   Verdict: {verdict_label}. {verdict_text[:200]}...")

    lines += ["", "6. **Are there clear differences between d=2, d=3, d=4?**"]
    for label in ("truth_plus_medium_noise", "random_init"):
        for d in (2, 3, 4):
            vals = []
            for s in sched_order:
                subset = [
                    r for r in rows
                    if r["init_label"] == label and r["target_dim"] == d
                    and r["schedule_label"] == s
                    and isinstance(r["final_energy"], float) and math.isfinite(r["final_energy"])
                ]
                vals.append(
                    _format_field(sum(r["final_energy"] for r in subset) / len(subset))
                    if subset else "NaN"
                )
            lines.append(
                f"   {label} d={d}: " + " / ".join(f"{s}={v}" for s, v in zip(sched_order, vals))
            )

    lines += ["", "7. **Is n=64 systematically harder than n=32?**"]
    for label in ("truth_plus_medium_noise", "random_init"):
        for n in (32, 64):
            long_fe = []
            for d in (2, 3, 4):
                subset = [
                    r for r in rows
                    if r["init_label"] == label and r["n"] == n and r["target_dim"] == d
                    and r["schedule_label"] == "long"
                    and isinstance(r["final_energy"], float) and math.isfinite(r["final_energy"])
                ]
                long_fe.extend(r["final_energy"] for r in subset)
            if long_fe:
                lines.append(
                    f"   {label} n={n} (long schedule): "
                    f"mean final E = {_format_field(sum(long_fe) / len(long_fe))} "
                    f"over {len(long_fe)} valid cells."
                )

    lines += [
        "",
        "## Full per-run table",
        "",
        "| d | n | seed | init | sched | init E | post-W E | final E | ΔE | W-att | W-acc | W-rej | recov | t(s) |",
        "| :---: | ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :---: | ---: |",
    ]
    for r in rows:
        lines.append(
            "| {d} | {n} | {seed} | {label} | {sched} | {ie} | {pwe} | {fe} | {de}"
            " | {att} | {acc} | {rej} | {rec} | {rt} |".format(
                d=r["target_dim"],
                n=r["n"],
                seed=r["seed"],
                label=r["init_label"],
                sched=r["schedule_label"],
                ie=_format_field(r["initial_energy"]) if isinstance(r["initial_energy"], float) else "NaN",
                pwe=_format_field(r["post_warmup_energy"]) if isinstance(r["post_warmup_energy"], float) else "NaN",
                fe=_format_field(r["final_energy"]) if isinstance(r["final_energy"], float) and math.isfinite(r["final_energy"]) else "NaN",
                de=_format_field(r["delta_energy"]) if isinstance(r["delta_energy"], float) and math.isfinite(r["delta_energy"]) else "NaN",
                att=r["warmup_attempted_moves"],
                acc=r["warmup_accepted_moves"],
                rej=r["warmup_rejected_moves"],
                rec="✓" if r["recovered_flag"] else "✗",
                rt=_format_field(r["runtime_seconds"]),
            )
        )

    lines += [
        "",
        "Regenerate via `make regen-phase2g`. Source tool:",
        "`tools/build_phase2g_guarded_anneal_budget_scaling.py`.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    t0 = time.perf_counter()
    rows = build_rows()
    total_time = time.perf_counter() - t0

    FOUNDATION.mkdir(parents=True, exist_ok=True)
    write_csv(rows, FOUNDATION / "phase2g_guarded_anneal_budget_scaling.csv")
    write_markdown(rows, FOUNDATION / "phase2g_guarded_anneal_budget_scaling.md")

    verdict_label, _ = _verdict(rows)
    overflow_count = sum(
        1 for r in rows
        if not isinstance(r["final_energy"], float) or not math.isfinite(r["final_energy"])
    )
    recovered_total = sum(1 for r in rows if isinstance(r["recovered_flag"], bool) and r["recovered_flag"])

    print(
        f"Wrote {len(rows)} Phase 2G rows to {FOUNDATION} in {total_time:.1f}s. "
        f"recovered: {recovered_total}/{len(rows)}. "
        f"overflow_rows: {overflow_count}. "
        f"Verdict: {verdict_label}."
    )


if __name__ == "__main__":
    main()
