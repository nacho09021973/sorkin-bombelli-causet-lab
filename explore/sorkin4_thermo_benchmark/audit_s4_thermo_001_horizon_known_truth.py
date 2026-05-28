#!/usr/bin/env python3
"""S4-THERMO-001-HORIZON-KNOWN-TRUTH-AUDIT.

Level-A closed-form Schwarzschild/Kerr horizon geometry and thermodynamic
scalar guardrail. This is not Level-B discrete rediscovery.
"""

from __future__ import annotations

import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt


ARTIFACT_DIR = Path(__file__).resolve().parent
OUT_PREFIX = "s4_thermo_001_horizon_known_truth"
CHIS = [0.0, 0.25, 0.5, 0.75, 0.9, 0.99]
MASS = 1.0


def kerr_row(mass: float, chi: float) -> dict[str, float | bool]:
    a = chi * mass
    disc = math.sqrt(mass * mass - a * a)
    r_plus = mass + disc
    r_minus = mass - disc
    area = 4.0 * math.pi * (r_plus * r_plus + a * a)
    omega_h = a / (r_plus * r_plus + a * a)
    kappa = (r_plus - r_minus) / (2.0 * (r_plus * r_plus + a * a))
    t_h = kappa / (2.0 * math.pi)
    entropy = area / 4.0

    area_s = 16.0 * math.pi * mass * mass
    entropy_s = 4.0 * math.pi * mass * mass
    kappa_s = 1.0 / (4.0 * mass)
    t_s = 1.0 / (8.0 * math.pi * mass)
    return {
        "M": mass,
        "chi": chi,
        "a": a,
        "r_plus": r_plus,
        "r_minus": r_minus,
        "area": area,
        "entropy_BH": entropy,
        "kappa": kappa,
        "T_H": t_h,
        "Omega_H": omega_h,
        "area_over_16piM2": area / area_s,
        "entropy_over_4piM2": entropy / entropy_s,
        "kappa_over_schwarzschild": kappa / kappa_s,
        "T_over_schwarzschild": t_h / t_s,
        "omegaH_times_M": omega_h * mass,
    }


def compute_rows() -> list[dict[str, float | bool]]:
    rows = [kerr_row(MASS, chi) for chi in CHIS]
    sch = rows[0]
    tol = 1.0e-12
    sch_pass = (
        abs(sch["r_plus"] - 2.0 * MASS) <= tol
        and abs(sch["r_minus"] - 0.0) <= tol
        and abs(sch["area"] - 16.0 * math.pi * MASS * MASS) <= tol
        and abs(sch["kappa"] - 1.0 / (4.0 * MASS)) <= tol
        and abs(sch["T_H"] - 1.0 / (8.0 * math.pi * MASS)) <= tol
        and abs(sch["entropy_BH"] - 4.0 * math.pi * MASS * MASS) <= tol
        and abs(sch["Omega_H"] - 0.0) <= tol
    )

    areas = [float(r["area"]) for r in rows]
    ent = [float(r["entropy_BH"]) for r in rows]
    kap = [float(r["kappa"]) for r in rows]
    temp = [float(r["T_H"]) for r in rows]
    omg = [float(r["Omega_H"]) for r in rows]
    mono_area = all(areas[i + 1] < areas[i] for i in range(len(areas) - 1))
    mono_ent = all(ent[i + 1] < ent[i] for i in range(len(ent) - 1))
    mono_kap = all(kap[i + 1] < kap[i] for i in range(len(kap) - 1))
    mono_temp = all(temp[i + 1] < temp[i] for i in range(len(temp) - 1))
    mono_omg = all(omg[i + 1] > omg[i] for i in range(len(omg) - 1))

    area_ext = 8.0 * math.pi * MASS * MASS
    entropy_ext = 2.0 * math.pi * MASS * MASS
    omega_ext = 1.0 / (2.0 * MASS)
    row_09 = rows[4]
    row_099 = rows[5]
    extremal_trend_pass = (
        row_099["kappa"] < row_09["kappa"]
        and row_099["T_H"] < row_09["T_H"]
        and abs(row_099["area"] - area_ext) < abs(row_09["area"] - area_ext)
        and abs(row_099["entropy_BH"] - entropy_ext) < abs(row_09["entropy_BH"] - entropy_ext)
        and abs(row_099["Omega_H"] - omega_ext) < abs(row_09["Omega_H"] - omega_ext)
    )

    fixed_m_caveat = True
    mono_family = mono_area and mono_ent and mono_kap and mono_temp and mono_omg
    for row in rows:
        row["schwarzschild_limit_pass"] = sch_pass if row["chi"] == 0.0 else True
        row["extremal_trend_pass"] = extremal_trend_pass
        row["monotonic_family_pass"] = mono_family
        row["fixed_M_family_caveat_present"] = fixed_m_caveat
        row["all_checks_pass"] = bool(
            row["schwarzschild_limit_pass"]
            and row["extremal_trend_pass"]
            and row["monotonic_family_pass"]
            and row["fixed_M_family_caveat_present"]
        )
    return rows


def write_csv(rows: list[dict[str, float | bool]], path: Path) -> None:
    fields = [
        "M", "chi", "a", "r_plus", "r_minus", "area", "entropy_BH", "kappa", "T_H", "Omega_H",
        "area_over_16piM2", "entropy_over_4piM2", "kappa_over_schwarzschild",
        "T_over_schwarzschild", "omegaH_times_M", "schwarzschild_limit_pass",
        "extremal_trend_pass", "monotonic_family_pass", "fixed_M_family_caveat_present",
        "all_checks_pass",
    ]
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)


def write_json(rows: list[dict[str, float | bool]], path: Path) -> None:
    areas = [float(r["area"]) for r in rows]
    ent = [float(r["entropy_BH"]) for r in rows]
    kap = [float(r["kappa"]) for r in rows]
    temp = [float(r["T_H"]) for r in rows]
    omg = [float(r["Omega_H"]) for r in rows]
    payload = {
        "benchmark": "S4-THERMO-001 horizon known-truth guardrail",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "rows": rows,
        "global_summary": {
            "total_rows": len(rows),
            "all_checks_pass": all(bool(r["all_checks_pass"]) for r in rows),
            "schwarzschild_limit_pass": bool(rows[0]["schwarzschild_limit_pass"]),
            "extremal_trend_pass": bool(rows[0]["extremal_trend_pass"]),
            "monotonic_area_decrease_pass": all(areas[i + 1] < areas[i] for i in range(len(areas) - 1)),
            "monotonic_entropy_decrease_pass": all(ent[i + 1] < ent[i] for i in range(len(ent) - 1)),
            "monotonic_temperature_decrease_pass": all(temp[i + 1] < temp[i] for i in range(len(temp) - 1)),
            "monotonic_omegaH_increase_pass": all(omg[i + 1] > omg[i] for i in range(len(omg) - 1)),
            "level_A_only": True,
            "level_B_discrete_rediscovery_claimed": False,
        },
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def write_md(rows: list[dict[str, float | bool]], path: Path) -> None:
    lines = [
        "# S4-THERMO-001 horizon known-truth guardrail",
        "",
        "This is a Level-A closed-form identity audit.",
        "It checks standard Schwarzschild/Kerr horizon geometry and thermodynamic scalars.",
        "It does not derive Hawking radiation.",
        "It does not test discrete causal-set emergence.",
        "It does not constitute Level-B rediscovery.",
        "Fixed-M Kerr spin sweeps compare stationary solutions and are not physical dynamical evolution.",
        "",
        "Spin grid (chi=a/M):",
        ", ".join(str(r["chi"]) for r in rows),
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_png(rows: list[dict[str, float | bool]], path: Path) -> None:
    chi = [float(r["chi"]) for r in rows]
    area_norm = [float(r["area_over_16piM2"]) for r in rows]
    s_norm = [float(r["entropy_over_4piM2"]) for r in rows]
    k_norm = [float(r["kappa_over_schwarzschild"]) for r in rows]
    t_norm = [float(r["T_over_schwarzschild"]) for r in rows]
    om = [float(r["omegaH_times_M"]) for r in rows]

    fig, axs = plt.subplots(2, 2, figsize=(10, 7))
    fig.suptitle("S4-THERMO-001 horizon known-truth guardrail")

    axs[0, 0].plot(chi, area_norm, "o-")
    axs[0, 0].set_title("area/(16*pi*M^2)")
    axs[0, 0].set_xlabel("chi")

    axs[0, 1].plot(chi, s_norm, "o-")
    axs[0, 1].set_title("S_BH/(4*pi*M^2)")
    axs[0, 1].set_xlabel("chi")

    axs[1, 0].plot(chi, k_norm, "o-", label="kappa/kappa_schw")
    axs[1, 0].plot(chi, t_norm, "s-", label="T_H/T_schw")
    axs[1, 0].set_title("kappa and T_H normalized")
    axs[1, 0].set_xlabel("chi")
    axs[1, 0].legend(fontsize=8)

    axs[1, 1].plot(chi, om, "o-")
    axs[1, 1].set_title("Omega_H*M")
    axs[1, 1].set_xlabel("chi")

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main() -> None:
    rows = compute_rows()
    csv_path = ARTIFACT_DIR / f"{OUT_PREFIX}.csv"
    json_path = ARTIFACT_DIR / f"{OUT_PREFIX}.json"
    md_path = ARTIFACT_DIR / f"{OUT_PREFIX}.md"
    png_path = ARTIFACT_DIR / f"{OUT_PREFIX}.png"
    write_csv(rows, csv_path)
    write_json(rows, json_path)
    write_md(rows, md_path)
    write_png(rows, png_path)
    print(f"wrote {csv_path.name},{json_path.name},{md_path.name},{png_path.name}")


if __name__ == "__main__":
    main()
