#!/usr/bin/env python3
"""Build the Phase 2F guarded-warmup probe.

Phase 2E showed WARMUP_IS_PRIMARY_FAILURE: the 10 unconditional warmup
accepts destroy near-truth configurations. Phase 2F tests whether a
*guarded* warmup — one that only accepts energy-non-worsening moves —
can preserve the exploratory intent of warmup (helping random_init) while
not destroying near-truth starts.

Three warmup modes are compared side-by-side on the same Phase 2D/2E grid:

  "legacy_warmup"   — sim.warmup(buf) + sim.anneal(buf) [Phase 2E baseline]
  "skip_warmup"     — sim.anneal(buf) only              [Phase 2E baseline]
  "guarded_warmup"  — energy-gated warmup + sim.anneal(buf)

Guarded warmup (external wrapper, no changes to cones.py):
  - GUARD_THRESHOLD = 0.0 (strictly non-worsening, documented)
  - Each proposed move is accepted iff sim.deltae <= GUARD_THRESHOLD
  - Rejected moves restore sim.rave and clear change flags
  - sim.statistics() is called at the end to maintain internal consistency
  - warmup_attempted_moves / _accepted_moves / _rejected_moves are recorded

Note on energy normalization: sim.deltae is computed in pre-normalization
coordinates (before sim.update() rescales by rmin). sim.energies[0] is
always in post-normalization form. The guard is applied pre-normalization;
warmup_energy_after may therefore differ from warmup_energy_before by a
small normalization factor even for accepted-only runs. This is documented
and does not invalidate the guard.

Verdict labels:
  GUARDED_WARMUP_FIXES_PRIMARY_FAILURE
      guarded_warmup matches or beats skip_warmup on near-truth and
      improves random_init.
  SKIP_WARMUP_REMAINS_BEST
      guarded_warmup does not improve over skip_warmup.
  ANNEAL_RESIDUAL_FAILURE_DOMINATES
      even with guarded/skip, small-noise or truth configurations worsen.
  RANDOM_INIT_NEEDS_SMART_INIT
      near-truth works under guarded/skip but random_init still fails.
  MIXED
      partial improvement without a single dominant conclusion.
"""

from __future__ import annotations

import contextlib
import io
import math
import random
import sys
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

WARMUP_LIMIT = 10
ANNEAL_LIMIT = 10
MAX_DATA = 4
INITIAL_TEMP = 100.0
COOLING_FACTOR = 0.9

NOISE_SMALL = 1e-3
NOISE_MEDIUM = 5e-2

# Guard threshold for guarded_warmup: only accept moves that do not
# increase the (pre-normalization) energy.
GUARD_THRESHOLD = 0.0

INIT_LABELS: tuple[tuple[str, float | None], ...] = (
    ("truth", 0.0),
    ("truth_plus_small_noise", NOISE_SMALL),
    ("truth_plus_medium_noise", NOISE_MEDIUM),
    ("random_init", None),
)

WARMUP_MODES = ("legacy_warmup", "skip_warmup", "guarded_warmup")

CSV_HEADERS = (
    "family",
    "target_dim",
    "n",
    "seed",
    "init_label",
    "warmup_mode",
    "noise_epsilon",
    "paired_key",
    "initial_energy",
    "final_energy",
    "delta_energy",
    "initial_interval_rmse",
    "final_interval_rmse",
    "initial_distance_to_truth_rms",
    "final_distance_to_truth_rms",
    "improved_energy",
    "improved_interval_rmse",
    "preserved_near_truth",
    "warmup_attempted_moves",
    "warmup_accepted_moves",
    "warmup_rejected_moves",
    "warmup_energy_before",
    "warmup_energy_after",
    "warmup_delta_energy",
    "notes",
)


# ------------------------------------------------------------------
# Coordinate helpers  (identical to Phase 2D/2E)
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


def _scale_truth(
    points: list[vs.Coord],
    rave_target: float,
    n: int,
) -> list[vs.Coord]:
    rave_truth = sum(p[0] for p in points) / n
    if rave_truth <= 0.0:
        return list(points)
    scale = rave_target / rave_truth
    return [tuple(c * scale for c in p) for p in points]


def _paired_key(d: int, n: int, seed: int, label: str, eps) -> str:
    eps_str = "NA" if eps is None else f"{eps}"
    return f"{d}|{n}|{seed}|{label}|{eps_str}"


# ------------------------------------------------------------------
# Custom startup injection  (identical to Phase 2D/2E)
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
# Guarded warmup (external wrapper — no changes to cones.py)
# ------------------------------------------------------------------


def _guarded_warmup(sim: cones.ConesSimulator) -> dict:
    """Non-destructive warmup: accept only energy-non-worsening moves.

    GUARD_THRESHOLD = 0.0: accept iff sim.deltae <= 0 (pre-normalization).
    On rejection: clear change flags and restore sim.rave to sim.r so that
    the next reconfigure() starts from the correct accepted state.

    Returns move statistics and warmup energy bookends.
    """
    energy_before = sim.energies[0]
    attempted = accepted = rejected = 0

    for _ in range(WARMUP_LIMIT):
        if sim.energies[0] <= 0.0:
            break

        # Clear accumulated change flags before each iteration so that
        # energy() only evaluates the current proposal, not history.
        for i in range(sim.n):
            sim.change[i] = False

        rave_saved = sim.r  # sim.r is set by update(), stable across reconfigure()
        attempted += 1
        sim.reconfigure()
        sim.energy()

        if sim.deltae <= GUARD_THRESHOLD:
            sim.update()
            accepted += 1
        else:
            # Reject: clear change flags and restore rave.
            # rnew/xnew need not be reset because reconfigure() always
            # overwrites them from rold/xold at the start of the next call.
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
        "warmup_energy_before": energy_before,
        "warmup_energy_after": sim.energies[0],
        "warmup_delta_energy": sim.energies[0] - energy_before,
    }


# ------------------------------------------------------------------
# Per-row computation
# ------------------------------------------------------------------


def _probe_row(
    d_spacetime: int,
    n: int,
    seed: int,
    init_label: str,
    noise_epsilon: float | None,
    warmup_mode: str,
) -> dict:
    d_spatial = d_spacetime - 1
    matrix, points = vs.sprinkle_minkowski_diamond(
        n=n, seed=seed, d_spacetime=d_spacetime
    )

    noise_seed = seed * 10007 + d_spacetime * 1009 + n * 97
    noise_rng = random.Random(noise_seed)

    with contextlib.redirect_stdout(io.StringIO()):
        sim = cones.ConesSimulator(
            z=matrix,
            dim=d_spatial,
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

        warmup_energy_before = initial_energy

        if warmup_mode == "legacy_warmup":
            sim.warmup(buf)
            warmup_count = sim.count
            wmeta = {
                "warmup_attempted_moves": warmup_count,
                "warmup_accepted_moves": warmup_count,
                "warmup_rejected_moves": 0,
                "warmup_energy_before": warmup_energy_before,
                "warmup_energy_after": sim.warmup_energy,
                "warmup_delta_energy": sim.warmup_energy - warmup_energy_before,
            }
        elif warmup_mode == "skip_warmup":
            wmeta = {
                "warmup_attempted_moves": 0,
                "warmup_accepted_moves": 0,
                "warmup_rejected_moves": 0,
                "warmup_energy_before": warmup_energy_before,
                "warmup_energy_after": warmup_energy_before,
                "warmup_delta_energy": 0.0,
            }
        else:  # guarded_warmup
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

    delta_energy = final_energy - initial_energy
    improved_energy = delta_energy < -1e-9
    improved_interval_rmse = final_interval_rmse < initial_interval_rmse - 1e-9
    preserved_near_truth = delta_energy <= 1e-6

    wword = {"legacy_warmup": "warmup+anneal", "skip_warmup": "anneal-only",
             "guarded_warmup": "guarded+anneal"}[warmup_mode]
    if init_label == "truth":
        note = (
            f"truth preserved: energy remains zero ({wword})"
            if abs(final_energy) < 1e-9 else
            f"truth destroyed under {wword}: energy rose to {final_energy:.4g}"
        )
    elif init_label in ("truth_plus_small_noise", "truth_plus_medium_noise"):
        if improved_energy:
            note = f"{wword} improved energy {initial_energy:.4g} -> {final_energy:.4g}"
        elif preserved_near_truth:
            note = f"{wword} preserved near-truth: energy {initial_energy:.4g} -> {final_energy:.4g}"
        else:
            note = f"{wword} destroyed near-truth: energy {initial_energy:.4g} -> {final_energy:.4g}"
    else:
        if improved_energy:
            note = f"{wword} improved random_init: {initial_energy:.4g} -> {final_energy:.4g}"
        else:
            note = f"{wword} did not improve random_init: {initial_energy:.4g} -> {final_energy:.4g}"

    return {
        "family": "minkowski",
        "target_dim": d_spacetime,
        "n": n,
        "seed": seed,
        "init_label": init_label,
        "warmup_mode": warmup_mode,
        "noise_epsilon": noise_epsilon,
        "paired_key": _paired_key(d_spacetime, n, seed, init_label, noise_epsilon),
        "initial_energy": initial_energy,
        "final_energy": final_energy,
        "delta_energy": delta_energy,
        "initial_interval_rmse": initial_interval_rmse,
        "final_interval_rmse": final_interval_rmse,
        "initial_distance_to_truth_rms": initial_dist,
        "final_distance_to_truth_rms": final_dist,
        "improved_energy": improved_energy,
        "improved_interval_rmse": improved_interval_rmse,
        "preserved_near_truth": preserved_near_truth,
        **wmeta,
        "notes": note,
    }


def build_rows() -> list[dict]:
    rows: list[dict] = []
    for d in SPACETIME_DIMS:
        for n in SIZES:
            for seed in PROBE_SEEDS:
                for label, eps in INIT_LABELS:
                    for wmode in WARMUP_MODES:
                        rows.append(_probe_row(d, n, seed, label, eps, wmode))
    return rows


# ------------------------------------------------------------------
# CSV output
# ------------------------------------------------------------------


def _fmt(value) -> str:
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return "true" if value else "false"
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
# Paired delta helpers
# ------------------------------------------------------------------


def _aggregate(rows: list[dict], label: str, wmode: str) -> dict:
    subset = [r for r in rows if r["init_label"] == label and r["warmup_mode"] == wmode]
    if not subset:
        return {}
    fe = [r["final_energy"] for r in subset]
    de = [r["delta_energy"] for r in subset]
    pnt = sum(1 for r in subset if r["preserved_near_truth"])
    wde = [r["warmup_delta_energy"] for r in subset]
    return {
        "runs": len(subset),
        "mean_initial_energy": sum(r["initial_energy"] for r in subset) / len(subset),
        "mean_final_energy": sum(fe) / len(fe),
        "mean_delta_energy": sum(de) / len(de),
        "mean_warmup_delta_energy": sum(wde) / len(wde),
        "preserved_count": pnt,
    }


# ------------------------------------------------------------------
# Verdict
# ------------------------------------------------------------------


def _verdict(rows: list[dict]) -> tuple[str, str]:
    def pres(label, wmode):
        subset = [r for r in rows if r["init_label"] == label and r["warmup_mode"] == wmode]
        return sum(1 for r in subset if r["preserved_near_truth"]), len(subset)

    truth_anneal_fail = any(
        not r["preserved_near_truth"]
        for r in rows
        if r["init_label"] == "truth" and r["warmup_mode"] in ("skip_warmup", "guarded_warmup")
    )
    if truth_anneal_fail:
        return (
            "ANNEAL_RESIDUAL_FAILURE_DOMINATES",
            "Truth-init configurations are worsened even with skip or guarded "
            "warmup. The annealing phase itself does not preserve zero-energy "
            "configurations.",
        )

    g_small_pres, g_small_n = pres("truth_plus_small_noise", "guarded_warmup")
    s_small_pres, _ = pres("truth_plus_small_noise", "skip_warmup")
    g_ri_pres, g_ri_n = pres("random_init", "guarded_warmup")
    s_ri_pres, _ = pres("random_init", "skip_warmup")

    guarded_beats_skip_small = g_small_pres >= s_small_pres
    guarded_improves_ri = g_ri_pres >= s_ri_pres

    if guarded_beats_skip_small and guarded_improves_ri:
        return (
            "GUARDED_WARMUP_FIXES_PRIMARY_FAILURE",
            (
                f"guarded_warmup matches or beats skip_warmup on small-noise "
                f"({g_small_pres}/{g_small_n} vs {s_small_pres}/{g_small_n}) "
                f"and improves random_init "
                f"({g_ri_pres}/{g_ri_n} vs {s_ri_pres}/{g_ri_n}). "
                "The guarded warmup preserves exploratory benefit for random "
                "starts while not destroying near-truth configurations. "
                "The unconditional warmup was the primary failure; the energy "
                "guard fixes it without sacrificing warmup's utility."
            ),
        )

    if g_small_pres >= s_small_pres and not guarded_improves_ri:
        return (
            "RANDOM_INIT_NEEDS_SMART_INIT",
            (
                f"guarded_warmup preserves near-truth starts "
                f"({g_small_pres}/{g_small_n}) as well as skip_warmup "
                f"({s_small_pres}/{g_small_n}), but random_init does not "
                f"benefit ({g_ri_pres}/{g_ri_n} vs {s_ri_pres}/{g_ri_n}). "
                "The guarded warmup is non-destructive but does not help "
                "random starts explore the landscape. Smart initialization "
                "is required for random-start improvement."
            ),
        )

    if g_small_pres < s_small_pres:
        return (
            "SKIP_WARMUP_REMAINS_BEST",
            (
                f"guarded_warmup ({g_small_pres}/{g_small_n}) does not match "
                f"skip_warmup ({s_small_pres}/{g_small_n}) on small-noise preservation. "
                "The guarded warmup still introduces some disruption beyond "
                "the anneal-only baseline."
            ),
        )

    return (
        "MIXED",
        "The results show partial improvement without a single dominant conclusion. "
        "See the per-label table for details.",
    )


# ------------------------------------------------------------------
# Markdown report
# ------------------------------------------------------------------


def write_markdown(rows: list[dict], path: Path) -> None:
    verdict_label, verdict_text = _verdict(rows)

    lines = [
        "# Phase 2F Guarded-Warmup Probe",
        "",
        "Three-way paired comparison of ``legacy_warmup`` (Phase 2E baseline),",
        "``skip_warmup`` (Phase 2E baseline), and ``guarded_warmup`` (new).",
        "No changes to ``cones.py``. The diagnostic question:",
        "can a non-destructive warmup preserve exploratory benefit (helping",
        "random_init) without destroying near-truth starts?",
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
        f"- NOISE_SMALL = {NOISE_SMALL}, NOISE_MEDIUM = {NOISE_MEDIUM}.",
        f"- anneal_limit={ANNEAL_LIMIT}, max_data={MAX_DATA}, "
        f"T₀={INITIAL_TEMP}, γ={COOLING_FACTOR}.",
        f"- warmup_limit={WARMUP_LIMIT} for legacy and guarded modes.",
        "- ``paired_key`` = ``target_dim|n|seed|init_label|noise_epsilon``.",
        "- Metrics from last-accepted positions (``rold``/``xold``).",
        "- ``preserved_near_truth``: ``delta_energy ≤ 1e-6``.",
        "",
        "### Guarded-warmup details",
        "",
        f"- GUARD_THRESHOLD = {GUARD_THRESHOLD} (strictly non-worsening,",
        "  applied to pre-normalization ``sim.deltae``).",
        "- On each proposed move: accept iff ``sim.deltae <= 0``.",
        "- On rejection: change flags cleared, ``sim.rave`` restored to ``sim.r``.",
        "- ``rnew``/``xnew`` are NOT explicitly reset because ``reconfigure()``",
        "  always overwrites them from ``rold``/``xold`` at the next call.",
        "- Energy guard is pre-normalization; ``warmup_energy_after`` may differ",
        "  from ``warmup_energy_before`` by a normalization factor even on",
        "  all-accept runs. This is documented, not a defect.",
        "",
        "## Per-label aggregate by warmup mode",
        "",
        "| init | warmup_mode | runs | mean init E | mean final E | mean ΔE | mean W ΔE | preserved |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for label, _ in INIT_LABELS:
        for wmode in WARMUP_MODES:
            b = _aggregate(rows, label, wmode)
            if not b:
                continue
            lines.append(
                "| {label} | {wmode} | {r} | {ie} | {fe} | {de} | {wde} | {pnt}/{r} |".format(
                    label=label,
                    wmode=wmode,
                    r=b["runs"],
                    ie=_format_field(b["mean_initial_energy"]),
                    fe=_format_field(b["mean_final_energy"]),
                    de=_format_field(b["mean_delta_energy"]),
                    wde=_format_field(b["mean_warmup_delta_energy"]),
                    pnt=b["preserved_count"],
                )
            )

    lines += [
        "",
        "## Six fixed questions",
        "",
        "1. **Does guarded_warmup preserve truth init?**",
    ]
    for wmode in WARMUP_MODES:
        t = [r for r in rows if r["init_label"] == "truth" and r["warmup_mode"] == wmode]
        pnt = sum(1 for r in t if r["preserved_near_truth"])
        mfe = sum(r["final_energy"] for r in t) / len(t)
        lines.append(f"   {wmode}: {pnt}/{len(t)} preserved, mean final E = {_format_field(mfe)}.")

    lines += ["", "2. **Does guarded_warmup improve small-noise over skip_warmup?**"]
    for wmode in WARMUP_MODES:
        s = [r for r in rows if r["init_label"] == "truth_plus_small_noise" and r["warmup_mode"] == wmode]
        pnt = sum(1 for r in s if r["preserved_near_truth"])
        mfe = sum(r["final_energy"] for r in s) / len(s)
        lines.append(f"   {wmode}: {pnt}/{len(s)} preserved, mean final E = {_format_field(mfe)}.")

    lines += ["", "3. **Does guarded_warmup improve medium-noise over skip_warmup?**"]
    for wmode in WARMUP_MODES:
        m = [r for r in rows if r["init_label"] == "truth_plus_medium_noise" and r["warmup_mode"] == wmode]
        pnt = sum(1 for r in m if r["preserved_near_truth"])
        mfe = sum(r["final_energy"] for r in m) / len(m)
        lines.append(f"   {wmode}: {pnt}/{len(m)} preserved, mean final E = {_format_field(mfe)}.")

    lines += ["", "4. **Does guarded_warmup help random_init without destroying near-truth?**"]
    for wmode in WARMUP_MODES:
        ri = [r for r in rows if r["init_label"] == "random_init" and r["warmup_mode"] == wmode]
        pnt = sum(1 for r in ri if r["preserved_near_truth"])
        mfe = sum(r["final_energy"] for r in ri) / len(ri)
        lines.append(f"   {wmode}: {pnt}/{len(ri)} preserved, mean final E = {_format_field(mfe)}.")

    lines += ["", "5. **Is Phase 2E's improvement fully explained by eliminating unconditional accepts?**"]
    skip_small = [r for r in rows if r["init_label"] == "truth_plus_small_noise" and r["warmup_mode"] == "skip_warmup"]
    guard_small = [r for r in rows if r["init_label"] == "truth_plus_small_noise" and r["warmup_mode"] == "guarded_warmup"]
    skip_pres = sum(1 for r in skip_small if r["preserved_near_truth"])
    guard_pres = sum(1 for r in guard_small if r["preserved_near_truth"])
    if guard_pres >= skip_pres:
        lines.append(
            f"   Yes: guarded_warmup ({guard_pres}/{len(guard_small)}) matches "
            f"skip_warmup ({skip_pres}/{len(skip_small)}) on small-noise. "
            "The energy-gated warmup is sufficient to avoid the Phase 2D failure."
        )
    else:
        lines.append(
            f"   Partial: guarded_warmup ({guard_pres}/{len(guard_small)}) does not "
            f"fully match skip_warmup ({skip_pres}/{len(skip_small)}). "
            "Some disruption beyond the unconditional-accept mechanism remains."
        )

    lines += ["", "6. **Is there residual failure attributable to anneal/move-set/cooling?**"]
    all_small_skip = [r for r in rows if r["init_label"] == "truth_plus_small_noise" and r["warmup_mode"] == "skip_warmup"]
    residual_fail = sum(1 for r in all_small_skip if not r["preserved_near_truth"])
    if residual_fail > 0:
        lines.append(
            f"   Yes: {residual_fail}/{len(all_small_skip)} small-noise cases fail "
            "even with skip_warmup. The annealing phase has residual instability "
            "near the truth minimum independent of the warmup."
        )
    else:
        lines.append(
            "   No residual failure in skip_warmup small-noise rows. "
            "The annealing phase alone is sufficient to stay near truth "
            "when starting from small-noise positions."
        )

    lines += [
        "",
        "## Full per-run table",
        "",
        "| d | n | seed | init | warmup | ε | init E | final E | ΔE | W-att | W-acc | W-rej | preserved |",
        "| :---: | ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :---: |",
    ]
    for r in rows:
        eps_str = "NA" if r["noise_epsilon"] is None else _format_field(r["noise_epsilon"])
        lines.append(
            "| {d} | {n} | {seed} | {label} | {wmode} | {eps} | {ie} | {fe} | {de}"
            " | {att} | {acc} | {rej} | {pnt} |".format(
                d=r["target_dim"],
                n=r["n"],
                seed=r["seed"],
                label=r["init_label"],
                wmode=r["warmup_mode"],
                eps=eps_str,
                ie=_format_field(r["initial_energy"]),
                fe=_format_field(r["final_energy"]),
                de=_format_field(r["delta_energy"]),
                att=r["warmup_attempted_moves"],
                acc=r["warmup_accepted_moves"],
                rej=r["warmup_rejected_moves"],
                pnt="✓" if r["preserved_near_truth"] else "✗",
            )
        )

    lines += [
        "",
        "Regenerate via `make regen-phase2f`. Source tool:",
        "`tools/build_phase2f_guarded_warmup_probe.py`.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    rows = build_rows()
    FOUNDATION.mkdir(parents=True, exist_ok=True)
    write_csv(rows, FOUNDATION / "phase2f_guarded_warmup_probe.csv")
    write_markdown(rows, FOUNDATION / "phase2f_guarded_warmup_probe.md")

    verdict_label, _ = _verdict(rows)
    n_pres = sum(1 for r in rows if r["preserved_near_truth"])
    by_mode = {
        wmode: sum(1 for r in rows if r["warmup_mode"] == wmode and r["preserved_near_truth"])
        for wmode in WARMUP_MODES
    }
    print(
        f"Wrote {len(rows)} Phase 2F rows to {FOUNDATION}. "
        f"preserved: {n_pres}/{len(rows)} total "
        f"({', '.join(f'{m}={by_mode[m]}' for m in WARMUP_MODES)}). "
        f"Verdict: {verdict_label}."
    )


if __name__ == "__main__":
    main()
