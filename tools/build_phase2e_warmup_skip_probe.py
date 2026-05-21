#!/usr/bin/env python3
"""Build the Phase 2E warmup-skip probe.

Phase 2D identified the dominant failure mode: the historical warmup loop
makes unconditional accepts, destroying any near-truth configuration with
energy > 0 before annealing even starts. Phase 2E tests the direct fix:
skip the warmup entirely and run only the annealing phase.

The probe repeats the Phase 2D grid (same d, n, seed, init_label) with
two warmup modes side-by-side:

  "with_warmup"  — sim.warmup(buf) then sim.anneal(buf)  [Phase 2D baseline]
  "skip_warmup"  — sim.anneal(buf) only, no warmup call

Each physical case (target_dim, n, seed, init_label, noise_epsilon) is
identified by a ``paired_key`` so the two modes can be compared row-for-row
without ambiguity. The markdown aggregated table shows the paired deltas:
  delta_final_energy_skip_minus_with
  delta_final_interval_rmse_skip_minus_with
  delta_final_distance_skip_minus_with

No changes are made to cones.py.

Verdict labels:
  WARMUP_IS_PRIMARY_FAILURE   — skip_warmup clearly preserves small-noise
                                better than with_warmup.
  WARMUP_NOT_PRIMARY          — skip_warmup does not appreciably improve.
  ANNEAL_DESTROYS_NEAR_TRUTH  — even without warmup, truth or small-noise
                                configurations worsen.
  NARROW_BASIN_WITHOUT_WARMUP — truth preserved but small/medium noise
                                still not recovered without warmup.
  INIT_PROBLEM_CONFIRMED      — near-truth works, random_init fails.
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

INIT_LABELS: tuple[tuple[str, float | None], ...] = (
    ("truth", 0.0),
    ("truth_plus_small_noise", NOISE_SMALL),
    ("truth_plus_medium_noise", NOISE_MEDIUM),
    ("random_init", None),
)

WARMUP_MODES = ("with_warmup", "skip_warmup")

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
    "notes",
)


# ------------------------------------------------------------------
# Coordinate helpers  (identical to Phase 2D)
# ------------------------------------------------------------------


def _coord_distance_rms(
    coords: list[vs.Coord],
    truth_scaled: list[vs.Coord],
    n: int,
    d_spatial: int,
) -> float:
    sq_sum = 0.0
    for i in range(n):
        ci = coords[i]
        ti = truth_scaled[i]
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
# Custom startup injection  (identical to Phase 2D)
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

        if warmup_mode == "with_warmup":
            sim.warmup(buf)
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

    warmup_word = "warmup+anneal" if warmup_mode == "with_warmup" else "anneal-only"
    if init_label == "truth":
        if abs(final_energy) < 1e-9:
            note = f"truth preserved: energy remains zero ({warmup_word})"
        else:
            note = (
                f"truth destroyed under {warmup_word}: "
                f"energy rose to {final_energy:.4g}"
            )
    elif init_label in ("truth_plus_small_noise", "truth_plus_medium_noise"):
        if improved_energy:
            note = (
                f"{warmup_word} improved energy "
                f"{initial_energy:.4g} -> {final_energy:.4g}"
            )
        elif preserved_near_truth:
            note = (
                f"{warmup_word} preserved near-truth: "
                f"energy {initial_energy:.4g} -> {final_energy:.4g}"
            )
        else:
            note = (
                f"{warmup_word} destroyed near-truth: "
                f"energy {initial_energy:.4g} -> {final_energy:.4g}"
            )
    else:
        if improved_energy:
            note = (
                f"{warmup_word} improved random_init: "
                f"{initial_energy:.4g} -> {final_energy:.4g}"
            )
        else:
            note = (
                f"{warmup_word} did not improve random_init: "
                f"{initial_energy:.4g} -> {final_energy:.4g}"
            )

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


def _paired_deltas(rows: list[dict]) -> list[dict]:
    """For each paired_key return the skip-minus-with delta row."""
    by_key: dict[str, dict] = {}
    for r in rows:
        key = r["paired_key"]
        wmode = r["warmup_mode"]
        by_key.setdefault(key, {})[wmode] = r

    deltas = []
    for key, pair in by_key.items():
        if "with_warmup" not in pair or "skip_warmup" not in pair:
            continue
        w = pair["with_warmup"]
        s = pair["skip_warmup"]
        deltas.append(
            {
                "paired_key": key,
                "init_label": w["init_label"],
                "target_dim": w["target_dim"],
                "n": w["n"],
                "seed": w["seed"],
                "delta_final_energy_skip_minus_with": s["final_energy"] - w["final_energy"],
                "delta_final_interval_rmse_skip_minus_with": s["final_interval_rmse"] - w["final_interval_rmse"],
                "delta_final_distance_skip_minus_with": s["final_distance_to_truth_rms"] - w["final_distance_to_truth_rms"],
                "preserved_with": w["preserved_near_truth"],
                "preserved_skip": s["preserved_near_truth"],
            }
        )
    return deltas


# ------------------------------------------------------------------
# Verdict
# ------------------------------------------------------------------


def _verdict(rows: list[dict]) -> tuple[str, str]:
    skip_truth = [
        r for r in rows
        if r["init_label"] == "truth" and r["warmup_mode"] == "skip_warmup"
    ]
    skip_small = [
        r for r in rows
        if r["init_label"] == "truth_plus_small_noise"
        and r["warmup_mode"] == "skip_warmup"
    ]
    with_small = [
        r for r in rows
        if r["init_label"] == "truth_plus_small_noise"
        and r["warmup_mode"] == "with_warmup"
    ]
    skip_ri = [
        r for r in rows
        if r["init_label"] == "random_init" and r["warmup_mode"] == "skip_warmup"
    ]
    with_ri = [
        r for r in rows
        if r["init_label"] == "random_init" and r["warmup_mode"] == "with_warmup"
    ]

    truth_destroyed = any(not r["preserved_near_truth"] for r in skip_truth)
    pres_skip_small = sum(1 for r in skip_small if r["preserved_near_truth"])
    pres_with_small = sum(1 for r in with_small if r["preserved_near_truth"])
    pres_skip_ri = sum(1 for r in skip_ri if r["preserved_near_truth"])
    pres_with_ri = sum(1 for r in with_ri if r["preserved_near_truth"])

    if truth_destroyed:
        return (
            "ANNEAL_DESTROYS_NEAR_TRUTH",
            (
                "Even without warmup, truth-init configurations are worsened "
                "by the annealing phase. The move set or acceptance rule does "
                "not preserve zero-energy configurations, independent of the "
                "warmup. This is a more severe failure than Phase 2D suggested."
            ),
        )

    if pres_skip_small == len(skip_small) and pres_skip_small > pres_with_small:
        return (
            "WARMUP_IS_PRIMARY_FAILURE",
            (
                f"Skipping the warmup fully recovers small-noise starts: "
                f"{pres_skip_small}/{len(skip_small)} preserved without warmup vs "
                f"{pres_with_small}/{len(with_small)} with warmup. "
                "The annealing phase alone, starting from a near-truth position, "
                "stays near the minimum. The warmup phase — with its unconditional "
                "accepts — is the primary cause of near-truth destruction. "
                "Skipping the warmup when initializing near a known minimum is "
                "the recommended practical fix. This fully confirms Phase 2D's "
                "NARROW_BASIN diagnosis."
            ),
        )

    if pres_skip_small > pres_with_small:
        return (
            "WARMUP_IS_PRIMARY_FAILURE",
            (
                f"Skipping the warmup improves small-noise preservation: "
                f"{pres_skip_small}/{len(skip_small)} vs "
                f"{pres_with_small}/{len(with_small)} with warmup. "
                "The warmup phase is the primary contributor to near-truth "
                "destruction. Some residual instability in the annealing phase "
                "remains, but skipping warmup is a clear improvement."
            ),
        )

    skip_medium = [
        r for r in rows
        if r["init_label"] == "truth_plus_medium_noise"
        and r["warmup_mode"] == "skip_warmup"
    ]
    pres_skip_med = sum(1 for r in skip_medium if r["preserved_near_truth"])
    if pres_skip_small == 0 and pres_skip_med == 0:
        if pres_skip_small > 0 or pres_skip_ri > pres_with_ri:
            return (
                "INIT_PROBLEM_CONFIRMED",
                (
                    "Near-truth starts are preserved by the anneal-only path, "
                    "but random_init does not benefit from skipping warmup. "
                    "The initialization strategy is a key factor: the optimizer "
                    "needs a good starting point to converge."
                ),
            )
        return (
            "NARROW_BASIN_WITHOUT_WARMUP",
            (
                f"Truth is preserved but small-noise starts are not recovered "
                f"even without warmup ({pres_skip_small}/{len(skip_small)}). "
                "The annealing phase itself has an extremely narrow basin of "
                "attraction. The warmup is not the sole cause: the move set or "
                "cooling schedule is also insufficient to stay near the truth."
            ),
        )

    return (
        "WARMUP_NOT_PRIMARY",
        (
            f"Skipping the warmup does not appreciably improve near-truth "
            f"preservation: {pres_skip_small}/{len(skip_small)} (skip) vs "
            f"{pres_with_small}/{len(with_small)} (with). "
            "The annealing phase itself is also failing near the truth minimum. "
            "The warmup is not the primary cause of the Phase 2D failure."
        ),
    )


# ------------------------------------------------------------------
# Markdown report
# ------------------------------------------------------------------


def _aggregate(rows: list[dict], label: str, wmode: str) -> dict:
    subset = [
        r for r in rows
        if r["init_label"] == label and r["warmup_mode"] == wmode
    ]
    if not subset:
        return {}
    ie = [r["initial_energy"] for r in subset]
    fe = [r["final_energy"] for r in subset]
    de = [r["delta_energy"] for r in subset]
    pnt = sum(1 for r in subset if r["preserved_near_truth"])
    return {
        "runs": len(subset),
        "mean_initial_energy": sum(ie) / len(ie),
        "mean_final_energy": sum(fe) / len(fe),
        "mean_delta_energy": sum(de) / len(de),
        "preserved_count": pnt,
    }


def write_markdown(rows: list[dict], path: Path) -> None:
    verdict_label, verdict_text = _verdict(rows)
    deltas = _paired_deltas(rows)

    lines = [
        "# Phase 2E Warmup-Skip Probe",
        "",
        "Paired comparison of ``with_warmup`` (Phase 2D baseline) vs",
        "``skip_warmup`` (anneal-only) on the same Phase 2D grid.",
        "No changes to ``cones.py``. The diagnostic question is:",
        "is the warmup the primary cause of near-truth destruction,",
        "or does the annealing phase also fail independently?",
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
        "- ``with_warmup``: warmup_limit=10 unconditional steps, then anneal.",
        "- ``skip_warmup``: anneal only — ``sim.warmup()`` not called.",
        "- Initialization is identical for the two modes at each ``paired_key``.",
        "- ``paired_key`` = ``target_dim|n|seed|init_label|noise_epsilon``.",
        "- Metrics from last-accepted positions (``rold``/``xold``).",
        "- ``preserved_near_truth``: ``delta_energy ≤ 1e-6``.",
        "",
        "## Per-label aggregate by warmup mode",
        "",
        "| init | warmup_mode | runs | mean init E | mean final E | mean ΔE | preserved |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for label, _ in INIT_LABELS:
        for wmode in WARMUP_MODES:
            b = _aggregate(rows, label, wmode)
            if not b:
                continue
            lines.append(
                "| {label} | {wmode} | {r} | {ie} | {fe} | {de} | {pnt}/{r} |".format(
                    label=label,
                    wmode=wmode,
                    r=b["runs"],
                    ie=_format_field(b["mean_initial_energy"]),
                    fe=_format_field(b["mean_final_energy"]),
                    de=_format_field(b["mean_delta_energy"]),
                    pnt=b["preserved_count"],
                )
            )

    # Paired delta table (per init_label, mean over all cells)
    lines += [
        "",
        "## Paired deltas (skip_warmup − with_warmup)",
        "",
        "Negative delta_final_energy means skip_warmup ends lower (better).",
        "",
        "| init | mean Δ final E | mean Δ final RMSE | mean Δ final dist | skip pres | with pres |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for label, _ in INIT_LABELS:
        subset = [d for d in deltas if d["init_label"] == label]
        if not subset:
            continue
        mean_de = sum(d["delta_final_energy_skip_minus_with"] for d in subset) / len(subset)
        mean_dr = sum(d["delta_final_interval_rmse_skip_minus_with"] for d in subset) / len(subset)
        mean_dd = sum(d["delta_final_distance_skip_minus_with"] for d in subset) / len(subset)
        skip_pres = sum(1 for d in subset if d["preserved_skip"])
        with_pres = sum(1 for d in subset if d["preserved_with"])
        lines.append(
            "| {label} | {de} | {dr} | {dd} | {sp}/{tot} | {wp}/{tot} |".format(
                label=label,
                de=_format_field(mean_de),
                dr=_format_field(mean_dr),
                dd=_format_field(mean_dd),
                sp=skip_pres,
                wp=with_pres,
                tot=len(subset),
            )
        )

    lines += [
        "",
        "## Five fixed questions",
        "",
        "1. **Does skipping warmup preserve truth init?**",
    ]
    for wmode in WARMUP_MODES:
        t = [r for r in rows if r["init_label"] == "truth" and r["warmup_mode"] == wmode]
        pnt = sum(1 for r in t if r["preserved_near_truth"])
        mean_fe = sum(r["final_energy"] for r in t) / len(t)
        lines.append(
            f"   {wmode}: {pnt}/{len(t)} preserved, "
            f"mean final_energy = {_format_field(mean_fe)}."
        )
    lines.append(
        "   Truth init has energy = 0; the warmup exits in 0 steps in both modes "
        "so both results should be identical."
    )

    lines += ["", "2. **Does skipping warmup preserve small-noise starts?**"]
    for wmode in WARMUP_MODES:
        s = [r for r in rows if r["init_label"] == "truth_plus_small_noise" and r["warmup_mode"] == wmode]
        pnt = sum(1 for r in s if r["preserved_near_truth"])
        mean_fe = sum(r["final_energy"] for r in s) / len(s)
        lines.append(
            f"   {wmode}: {pnt}/{len(s)} preserved, "
            f"mean final_energy = {_format_field(mean_fe)}."
        )

    lines += ["", "3. **Does skipping warmup help medium-noise starts?**"]
    for wmode in WARMUP_MODES:
        m = [r for r in rows if r["init_label"] == "truth_plus_medium_noise" and r["warmup_mode"] == wmode]
        pnt = sum(1 for r in m if r["preserved_near_truth"])
        mean_fe = sum(r["final_energy"] for r in m) / len(m)
        lines.append(
            f"   {wmode}: {pnt}/{len(m)} preserved, "
            f"mean final_energy = {_format_field(mean_fe)}."
        )

    lines += ["", "4. **Does random_init benefit or suffer from skipping warmup?**"]
    for wmode in WARMUP_MODES:
        ri = [r for r in rows if r["init_label"] == "random_init" and r["warmup_mode"] == wmode]
        pnt = sum(1 for r in ri if r["preserved_near_truth"])
        mean_fe = sum(r["final_energy"] for r in ri) / len(ri)
        lines.append(
            f"   {wmode}: {pnt}/{len(ri)} preserved, "
            f"mean final_energy = {_format_field(mean_fe)}."
        )
    lines.append(
        "   The warmup was designed to help random starts explore before cooling. "
        "If skip_warmup hurts random_init, that is the expected trade-off."
    )

    lines += ["", "5. **What is the dominant failure mode without warmup?**"]
    skip_small = [
        r for r in rows
        if r["init_label"] == "truth_plus_small_noise"
        and r["warmup_mode"] == "skip_warmup"
    ]
    pres = sum(1 for r in skip_small if r["preserved_near_truth"])
    if pres == len(skip_small):
        lines.append(
            "   All small-noise starts are preserved when warmup is skipped. "
            "The annealing phase alone is sufficient to stay near the minimum. "
            "The warmup is the sole identified cause of near-truth destruction."
        )
    elif pres > 0:
        lines.append(
            f"   {pres}/{len(skip_small)} small-noise starts preserved without warmup. "
            "Partial improvement: the warmup is a major contributor, but the "
            "annealing phase also has some instability near the truth minimum."
        )
    else:
        lines.append(
            "   No small-noise starts preserved even without warmup. "
            "The annealing phase itself is failing near the minimum. "
            "Move set or cooling schedule require further investigation."
        )

    lines += [
        "",
        "## Full per-run table",
        "",
        "| d | n | seed | init | warmup | ε | init E | final E | ΔE | preserved |",
        "| :---: | ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | :---: |",
    ]
    for r in rows:
        eps_str = "NA" if r["noise_epsilon"] is None else _format_field(r["noise_epsilon"])
        lines.append(
            "| {d} | {n} | {seed} | {label} | {wmode} | {eps} | {ie} | {fe} | {de} | {pnt} |".format(
                d=r["target_dim"],
                n=r["n"],
                seed=r["seed"],
                label=r["init_label"],
                wmode=r["warmup_mode"],
                eps=eps_str,
                ie=_format_field(r["initial_energy"]),
                fe=_format_field(r["final_energy"]),
                de=_format_field(r["delta_energy"]),
                pnt="✓" if r["preserved_near_truth"] else "✗",
            )
        )

    lines += [
        "",
        "Regenerate via `make regen-phase2e`. Source tool:",
        "`tools/build_phase2e_warmup_skip_probe.py`.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    rows = build_rows()
    FOUNDATION.mkdir(parents=True, exist_ok=True)
    write_csv(rows, FOUNDATION / "phase2e_warmup_skip_probe.csv")
    write_markdown(rows, FOUNDATION / "phase2e_warmup_skip_probe.md")

    verdict_label, _ = _verdict(rows)
    n_pres = sum(1 for r in rows if r["preserved_near_truth"])
    by_mode = {
        wmode: sum(1 for r in rows if r["warmup_mode"] == wmode and r["preserved_near_truth"])
        for wmode in WARMUP_MODES
    }
    print(
        f"Wrote {len(rows)} Phase 2E rows to {FOUNDATION}. "
        f"preserved: {n_pres}/{len(rows)} total "
        f"(with_warmup={by_mode['with_warmup']}, "
        f"skip_warmup={by_mode['skip_warmup']}). "
        f"Verdict: {verdict_label}."
    )


if __name__ == "__main__":
    main()
