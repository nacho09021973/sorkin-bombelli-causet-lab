#!/usr/bin/env python3
"""Post-analysis for the completed N=24 gamma probe.

Reads the existing probe CSV only. Does not run the annealer.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path


OUT_DIR = Path(__file__).resolve().parent
INPUT_CSV = OUT_DIR / "gamma_n24_probe.csv"
DERIVED_CSV = OUT_DIR / "gamma_n24_ranking_diagnostic.csv"
MD_PATH = OUT_DIR / "gamma_n24_ranking_diagnostic.md"
SVG_PATH = OUT_DIR / "gamma_n24_energy_vs_rmse.svg"


def _read_rows() -> list[dict[str, object]]:
    with INPUT_CSV.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))

    parsed: list[dict[str, object]] = []
    for row in rows:
        parsed.append(
            {
                "gamma": float(row["gamma"]),
                "final_energy": float(row["final_energy"]),
                "energy_gap": float(row["energy_gap"]),
                "interval_rmse": float(row["interval_rmse"]),
                "success_flag": row["success_flag"].strip().lower() == "true",
            }
        )
    return parsed


def _dense_ranks(rows: list[dict[str, object]], key: str) -> dict[float, int]:
    values = sorted({float(row[key]) for row in rows})
    return {value: index + 1 for index, value in enumerate(values)}


def _annotate(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    final_ranks = _dense_ranks(rows, "final_energy")
    gap_ranks = _dense_ranks(rows, "energy_gap")
    rmse_ranks = _dense_ranks(rows, "interval_rmse")
    best_final = min(float(row["final_energy"]) for row in rows)

    annotated = []
    for row in sorted(rows, key=lambda item: float(item["gamma"])):
        final_energy = float(row["final_energy"])
        energy_gap = float(row["energy_gap"])
        interval_rmse = float(row["interval_rmse"])
        annotated.append(
            {
                **row,
                "final_energy_rank": final_ranks[final_energy],
                "energy_gap_rank": gap_ranks[energy_gap],
                "interval_rmse_rank": rmse_ranks[interval_rmse],
                "success_rank": 1 if bool(row["success_flag"]) else 2,
                "within_10pct_final": final_energy <= best_final * 1.10 + 1e-9,
            }
        )
    return annotated


def _fmt_float(value: float) -> str:
    return f"{value:.6f}"


def write_csv(rows: list[dict[str, object]]) -> None:
    headers = (
        "gamma",
        "final_energy",
        "final_energy_rank",
        "energy_gap",
        "energy_gap_rank",
        "interval_rmse",
        "interval_rmse_rank",
        "success_flag",
        "success_rank",
        "within_10pct_final",
    )
    with DERIVED_CSV.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=headers)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "gamma": _fmt_float(float(row["gamma"])),
                    "final_energy": _fmt_float(float(row["final_energy"])),
                    "final_energy_rank": row["final_energy_rank"],
                    "energy_gap": _fmt_float(float(row["energy_gap"])),
                    "energy_gap_rank": row["energy_gap_rank"],
                    "interval_rmse": _fmt_float(float(row["interval_rmse"])),
                    "interval_rmse_rank": row["interval_rmse_rank"],
                    "success_flag": "true" if bool(row["success_flag"]) else "false",
                    "success_rank": row["success_rank"],
                    "within_10pct_final": (
                        "true" if bool(row["within_10pct_final"]) else "false"
                    ),
                }
            )


def write_svg(rows: list[dict[str, object]]) -> None:
    width = 980
    height = 560
    margin = 72
    split = width / 2
    panel_gap = 58
    panel_width = (width - 2 * margin - panel_gap) / 2
    panel_height = height - 2 * margin
    gammas = [float(row["gamma"]) for row in rows]
    energies = [float(row["final_energy"]) for row in rows]
    rmses = [float(row["interval_rmse"]) for row in rows]
    log_rmses = [math.log10(value) for value in rmses]

    def sx(gamma: float, x0: float) -> float:
        return x0 + (gamma - min(gammas)) * panel_width / (max(gammas) - min(gammas))

    def sy(value: float, values: list[float]) -> float:
        lo = min(values)
        hi = max(values)
        if lo == hi:
            lo -= 1.0
            hi += 1.0
        return height - margin - (value - lo) * panel_height / (hi - lo)

    left_x = margin
    right_x = split + panel_gap / 2
    best_energy = min(rows, key=lambda row: float(row["final_energy"]))
    best_rmse = min(rows, key=lambda row: float(row["interval_rmse"]))

    energy_points = " ".join(
        f"{sx(float(row['gamma']), left_x):.2f},{sy(float(row['final_energy']), energies):.2f}"
        for row in rows
    )
    rmse_points = " ".join(
        f"{sx(float(row['gamma']), right_x):.2f},{sy(math.log10(float(row['interval_rmse'])), log_rmses):.2f}"
        for row in rows
    )

    energy_circles = []
    rmse_circles = []
    labels = []
    for row in rows:
        gamma = float(row["gamma"])
        energy_fill = "#d1495b" if row is best_energy else "#2a9d8f"
        rmse_fill = "#d1495b" if row is best_rmse else "#2a9d8f"
        energy_circles.append(
            f"<circle cx='{sx(gamma, left_x):.2f}' cy='{sy(float(row['final_energy']), energies):.2f}' r='4.5' fill='{energy_fill}' />"
        )
        rmse_circles.append(
            f"<circle cx='{sx(gamma, right_x):.2f}' cy='{sy(math.log10(float(row['interval_rmse'])), log_rmses):.2f}' r='4.5' fill='{rmse_fill}' />"
        )
        labels.append(
            f"<text x='{sx(gamma, left_x):.2f}' y='{height - margin + 22}' text-anchor='middle' font-family='monospace' font-size='9'>{gamma:.3f}</text>"
        )
        labels.append(
            f"<text x='{sx(gamma, right_x):.2f}' y='{height - margin + 22}' text-anchor='middle' font-family='monospace' font-size='9'>{gamma:.3f}</text>"
        )

    svg = f"""<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>
  <rect width='100%' height='100%' fill='#f7f4ed'/>
  <text x='{width / 2:.0f}' y='34' text-anchor='middle' font-family='serif' font-size='22' fill='#222'>N=24 gamma ranking diagnostic</text>
  <text x='{left_x + panel_width / 2:.0f}' y='62' text-anchor='middle' font-family='monospace' font-size='13' fill='#222'>final energy</text>
  <text x='{right_x + panel_width / 2:.0f}' y='62' text-anchor='middle' font-family='monospace' font-size='13' fill='#222'>log10 interval RMSE</text>
  <line x1='{left_x}' y1='{height - margin}' x2='{left_x + panel_width}' y2='{height - margin}' stroke='#333' stroke-width='2'/>
  <line x1='{left_x}' y1='{margin}' x2='{left_x}' y2='{height - margin}' stroke='#333' stroke-width='2'/>
  <line x1='{right_x}' y1='{height - margin}' x2='{right_x + panel_width}' y2='{height - margin}' stroke='#333' stroke-width='2'/>
  <line x1='{right_x}' y1='{margin}' x2='{right_x}' y2='{height - margin}' stroke='#333' stroke-width='2'/>
  <polyline fill='none' stroke='#1f77b4' stroke-width='2.5' points='{energy_points}' />
  <polyline fill='none' stroke='#1f77b4' stroke-width='2.5' points='{rmse_points}' />
  {''.join(energy_circles)}
  {''.join(rmse_circles)}
  {''.join(labels)}
  <text x='{left_x}' y='{height - 18}' font-family='monospace' font-size='12' fill='#444'>red marks best by that panel's metric; lower is better</text>
</svg>
"""
    SVG_PATH.write_text(svg, encoding="utf-8")


def write_markdown(rows: list[dict[str, object]]) -> None:
    best_final = min(rows, key=lambda row: float(row["final_energy"]))
    best_rmse = min(rows, key=lambda row: float(row["interval_rmse"]))
    all_truth_zero = all(abs(float(row["energy_gap"]) - float(row["final_energy"])) < 1e-9 for row in rows)
    all_fail = all(not bool(row["success_flag"]) for row in rows)
    agree = float(best_final["gamma"]) == float(best_rmse["gamma"])

    lines = [
        "# N=24 Gamma Ranking Diagnostic",
        "",
        "**Status:** post-analysis of the completed exploratory N=24 CSV only. No annealer run is performed here.",
        "",
        "This note is an accessibility diagnostic for the historical annealer. It makes no embeddability claim, no physical gamma claim, and no success claim.",
        "",
        "## Provenance",
        "",
        f"- Input CSV: `{INPUT_CSV.relative_to(OUT_DIR.parents[1])}`",
        f"- Analysis script: `{Path(__file__).relative_to(OUT_DIR.parents[1])}`",
        f"- Derived CSV: `{DERIVED_CSV.relative_to(OUT_DIR.parents[1])}`",
        f"- Figure: `{SVG_PATH.relative_to(OUT_DIR.parents[1])}`",
        "- Exact command: `python3 explore/gamma_n24_probe/build_gamma_n24_ranking_diagnostic.py`",
        "",
        "## Metric Equivalence and Limits",
        "",
    ]
    if all_truth_zero:
        lines.append(
            "- `truth_energy` is zero for every row, so `final_energy` and `energy_gap` are equivalent in this run."
        )
    if all_fail:
        lines.append(
            "- Every `success_flag` value is false, so the success ranking cannot distinguish gamma values and cannot establish recovery."
        )
    lines += [
        "",
        "## Answers",
        "",
        f"- Best gamma by `final_energy`: `{float(best_final['gamma']):.6f}` with final energy `{float(best_final['final_energy']):.6f}`.",
        f"- Best gamma by `interval_rmse`: `{float(best_rmse['gamma']):.6f}` with interval RMSE `{float(best_rmse['interval_rmse']):.6f}`.",
        f"- Do the rankings agree? {'yes' if agree else 'no'}.",
        "- Apparent optimum: metric-dependent. The final-energy optimum sits in a broad 10% final-energy band, while interval RMSE favors a different gamma.",
        "- Next exploratory priority: stabilize the diagnostic metric first, before refining gamma or changing T0.",
        "",
        "## Ranking Table",
        "",
        "| gamma | final_energy | final rank | energy_gap | gap rank | interval_rmse | RMSE rank | success | within 10% final |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | :---: | :---: |",
    ]
    for row in rows:
        lines.append(
            "| {gamma:.6f} | {final:.6f} | {final_rank} | {gap:.6f} | {gap_rank} | {rmse:.6f} | {rmse_rank} | {success} | {window} |".format(
                gamma=float(row["gamma"]),
                final=float(row["final_energy"]),
                final_rank=int(row["final_energy_rank"]),
                gap=float(row["energy_gap"]),
                gap_rank=int(row["energy_gap_rank"]),
                rmse=float(row["interval_rmse"]),
                rmse_rank=int(row["interval_rmse_rank"]),
                success="yes" if bool(row["success_flag"]) else "no",
                window="yes" if bool(row["within_10pct_final"]) else "no",
            )
        )
    lines += [
        "",
        "## Conservative Readout",
        "",
        "- `final_energy` and `energy_gap` select the same gamma because the ground-truth energy is zero in every row.",
        "- `interval_rmse` selects a different gamma, so the apparent gamma optimum is not robust across diagnostics.",
        "- Since all success flags are false, this post-analysis does not establish recovery.",
    ]
    MD_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    rows = _annotate(_read_rows())
    write_csv(rows)
    write_svg(rows)
    write_markdown(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
