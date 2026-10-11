"""산화 모듈의 그래프와 결과 표를 모두 다시 만든다.

실행: python scripts/run_oxidation.py
출력: figures/oxidation/*.png, results/oxidation/*.csv
"""

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from isrc_sim.constants import K_B_EV, celsius_to_kelvin  # noqa: E402
from isrc_sim.oxidation import (  # noqa: E402
    deal_grove_rates,
    deal_grove_thickness,
    deal_grove_time,
    massoud_thickness,
    massoud_time,
    regime_crossover_thickness,
    si_consumed,
)

FIG_DIR = ROOT / "figures" / "oxidation"
RES_DIR = ROOT / "results" / "oxidation"
REF_CSV = ROOT / "data" / "reference" / "nanofab_sithox_2026-10-11.csv"

# 색상: 범주형은 고정 순서, 온도처럼 순서가 있는 값은 한 색상의 명도 단계로 표현한다.
TEXT = "#0b0b0b"
TEXT_2 = "#52514e"
GRID = "#e4e3df"
SURFACE = "#fcfcfb"
CAT = ["#2a78d6", "#eb6834", "#1baf7a"]          # blue, orange, aqua
BLUE_RAMP = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": TEXT_2, "axes.labelcolor": TEXT, "xtick.color": TEXT_2, "ytick.color": TEXT_2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "lines.linewidth": 2.0, "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold",
    "legend.frameon": False, "figure.dpi": 150,
})

TEMPS_C = [900, 1000, 1100, 1200]


def save(fig, name):
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_DIR / name, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved figures/oxidation/{name}")


def fig1_thickness_vs_time():
    t = np.linspace(0, 4, 400)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, ambient, title in ((axes[0], "dry", "Dry O$_2$, (100) Si"), (axes[1], "wet", "Wet (H$_2$O), (100) Si")):
        for T, c in zip(TEMPS_C, BLUE_RAMP):
            x_nm = deal_grove_thickness(t, T, ambient) * 1e3
            ax.plot(t, x_nm, color=c, label=f"{T} °C")
            ax.annotate(f"{T} °C", (t[-1], x_nm[-1]), xytext=(4, 0), textcoords="offset points",
                        va="center", fontsize=8.5, color=TEXT_2)
        ax.set(title=title, xlabel="Oxidation time (h)", ylabel="Oxide thickness (nm)", xlim=(0, 4.6))
        ax.set_ylim(bottom=0)
    axes[0].legend(loc="upper left")
    fig.suptitle("Deal-Grove: oxide thickness vs. time", color=TEXT, fontweight="bold")
    fig.tight_layout()
    save(fig, "fig1_thickness_vs_time.png")


def fig2_arrhenius():
    T = np.linspace(800, 1200, 100)
    inv_T = 1000.0 / celsius_to_kelvin(T)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ambient, c in (("dry", CAT[0]), ("wet", CAT[1])):
        B, _, _ = deal_grove_rates(T, ambient, "100")
        axes[0].semilogy(inv_T, B, color=c, label=f"{ambient}")
        for orient, ls in (("100", "-"), ("111", "--")):
            _, BA, _ = deal_grove_rates(T, ambient, orient)
            axes[1].semilogy(inv_T, BA, color=c, ls=ls, label=f"{ambient}, ({orient})")
    axes[0].set(title="Parabolic constant B  (E$_a$: dry 1.23 eV, wet 0.78 eV)",
                xlabel="1000 / T (1/K)", ylabel="B (µm²/h)")
    axes[1].set(title="Linear constant B/A  (E$_a$: dry 2.00 eV, wet 2.05 eV)",
                xlabel="1000 / T (1/K)", ylabel="B/A (µm/h)")
    for ax in axes:
        ax.legend()
        top = ax.secondary_xaxis("top", functions=(lambda v: 1000.0 / np.asarray(v) - 273.15,
                                                   lambda v: 1000.0 / (np.asarray(v) + 273.15)))
        top.set_xlabel("T (°C)", color=TEXT_2)
    fig.tight_layout()
    save(fig, "fig2_arrhenius.png")


def fig3_universal_curve():
    """x^2 + A x = B t 를 u = x/A, s = B t / A^2 로 무차원화하면 u^2 + u = s 하나로 모인다."""
    s = np.logspace(-3, 3, 300)
    u = (-1 + np.sqrt(1 + 4 * s)) / 2
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    ax.loglog(s, s, color=TEXT_2, ls=":", lw=1.5)
    ax.loglog(s, np.sqrt(s), color=TEXT_2, ls="--", lw=1.5)
    ax.loglog(s, u, color=CAT[0], lw=2.5, label="Deal-Grove: u² + u = s")
    ax.text(1.3e-2, 1.5e-3, "linear limit  u = s\n(reaction-limited)", color=TEXT_2, fontsize=8.5)
    ax.text(15, 1.2, "parabolic limit  u = √s\n(diffusion-limited)", color=TEXT_2, fontsize=8.5)
    markers = {"dry": "o", "wet": "s"}
    for ambient, c in (("dry", CAT[0]), ("wet", CAT[1])):
        for T in TEMPS_C:
            B, _, A = deal_grove_rates(T, ambient)
            t = np.logspace(-2, 1.3, 7)
            x = deal_grove_thickness(t, T, ambient)
            ax.loglog(B * t / A**2, x / A, ls="none", marker=markers[ambient], ms=5,
                      mfc=SURFACE, mec=c, mew=1.4)
    ax.plot([], [], ls="none", marker="o", mfc=SURFACE, mec=CAT[0], label="dry, 900~1200 °C, 0.01~20 h")
    ax.plot([], [], ls="none", marker="s", mfc=SURFACE, mec=CAT[1], label="wet, 900~1200 °C, 0.01~20 h")
    ax.axvline(0.75, color=GRID, lw=1.2)
    ax.text(0.8, 2e-3, "x = A/2\n(crossover)", color=TEXT_2, fontsize=8.5)
    ax.set(title="Every condition collapses onto one curve",
           xlabel="s = B·t / A²  (dimensionless time)", ylabel="u = x / A  (dimensionless thickness)",
           xlim=(1e-3, 1e3), ylim=(1e-3, 50))
    ax.legend(loc="upper left")
    fig.tight_layout()
    save(fig, "fig3_universal_curve.png")


def fig4_recipe_map():
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    T_dry = np.linspace(800, 1000, 41)
    for x_nm, c in zip((10, 20, 50), BLUE_RAMP[1:]):
        t_min = [massoud_time(x_nm * 1e-3, T) * 60 for T in T_dry]
        axes[0].semilogy(T_dry, t_min, color=c, label=f"{x_nm} nm")
    axes[0].set(title="Gate oxide, dry O$_2$ (Massoud, 800~1000 °C)",
                xlabel="Temperature (°C)", ylabel="Required time (min)")
    T_wet = np.linspace(900, 1200, 61)
    for x_nm, c in zip((300, 500, 1000), BLUE_RAMP[1:]):
        t_min = deal_grove_time(x_nm * 1e-3, T_wet, "wet") * 60
        axes[1].semilogy(T_wet, t_min, color=c, label=f"{x_nm} nm")
    axes[1].set(title="Field oxide, wet H$_2$O (Deal-Grove)",
                xlabel="Temperature (°C)", ylabel="Required time (min)")
    for ax in axes:
        ax.legend(title="target", loc="upper right")
    fig.suptitle("Recipe map: time needed to reach a target thickness, (100) Si", color=TEXT, fontweight="bold")
    fig.tight_layout()
    save(fig, "fig4_recipe_map.png")


def fig5_thin_oxide(T=900):
    t = np.linspace(0, 4, 400)
    x_dg = deal_grove_thickness(t, T, "dry") * 1e3
    x_dg25 = deal_grove_thickness(t, T, "dry", x_i_um=0.025) * 1e3
    x_ma = massoud_thickness(t, T) * 1e3
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    ax = axes[0]
    ax.plot(t, x_ma, color=CAT[2], label="Massoud (thin-oxide model)")
    ax.plot(t, x_dg25, color=CAT[1], ls="--", label="Deal-Grove, x$_i$ = 25 nm offset")
    ax.plot(t, x_dg, color=CAT[0], label="Deal-Grove, x$_i$ = 0")
    ax.set(title=f"Dry O$_2$, (100) Si, {T} °C", xlabel="Oxidation time (h)", ylabel="Oxide thickness (nm)")
    ax.set_ylim(bottom=0)
    ax.legend(loc="upper left")
    # 오른쪽: 시간을 로그 축으로 넓혀 두 Deal-Grove 근사가 Massoud에 수렴하는 과정을 본다.
    ax = axes[1]
    t_long = np.logspace(-1, np.log10(20), 300)
    m = massoud_thickness(t_long, T)
    err0 = (deal_grove_thickness(t_long, T, "dry") - m) / m * 100
    err25 = (deal_grove_thickness(t_long, T, "dry", x_i_um=0.025) - m) / m * 100
    ax.semilogx(t_long, err25, color=CAT[1], ls="--", label="Deal-Grove, x$_i$ = 25 nm offset")
    ax.semilogx(t_long, err0, color=CAT[0], label="Deal-Grove, x$_i$ = 0")
    ax.axhline(0, color=TEXT_2, lw=1)
    ax.set_ylim(-100, 150)
    for t_mark in (1, 4, 20):
        ax.annotate(f"{massoud_thickness(t_mark, T) * 1e3:.0f} nm", (t_mark, -100), xytext=(0, 4),
                    textcoords="offset points", ha="center", fontsize=8, color=TEXT_2)
    ax.set(title=f"Deviation from Massoud, {T} °C",
           xlabel="Oxidation time (h, log scale)   [labels: Massoud thickness]", ylabel="Error vs. Massoud (%)")
    ax.legend(loc="upper right")
    fig.tight_layout()
    save(fig, "fig5_thin_oxide.png")


def validation():
    rows = []
    with open(REF_CSV, newline="") as f:
        for r in csv.DictReader(f):
            t_h = float(r["time_min"]) / 60
            T = float(r["T_C"])
            x_i = float(r["xi_A"]) * 1e-4
            if r["model"] == "Massoud":
                x = massoud_thickness(t_h, T, x_i_um=x_i)
            else:
                x = deal_grove_thickness(t_h, T, r["ambient"], r["orientation"], x_i_um=x_i,
                                         pressure_atm=float(r["pp_atm"]))
            ref = float(r["thickness_A"])
            rows.append({**r, "model_A": round(x * 1e4, 1), "error_pct": round((x * 1e4 - ref) / ref * 100, 3)})
    RES_DIR.mkdir(parents=True, exist_ok=True)
    with open(RES_DIR / "validation_vs_nanofab.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("  saved results/oxidation/validation_vs_nanofab.csv")

    fig, ax = plt.subplots(figsize=(5.2, 5))
    style = {"Deal-Grove": ("o", CAT[0]), "Massoud": ("D", CAT[2])}
    for model, (m, c) in style.items():
        sub = [r for r in rows if r["model"] == model]
        ax.loglog([float(r["thickness_A"]) / 10 for r in sub], [r["model_A"] / 10 for r in sub],
                  ls="none", marker=m, ms=7, mfc=c, mec=SURFACE, mew=1.5, label=f"{model} ({len(sub)} cases)")
    lim = (3, 1500)
    ax.plot(lim, lim, color=TEXT_2, lw=1, ls=":")
    worst = max(abs(r["error_pct"]) for r in rows)
    ax.set(xlim=lim, ylim=lim, xlabel="nanoFAB calculator (nm)", ylabel="This project (nm)",
           title=f"Validation: max |error| = {worst:.2f} %")
    ax.legend(loc="upper left")
    fig.tight_layout()
    save(fig, "fig6_validation_parity.png")
    return rows


def recipe_examples():
    """대표 공정 목표에 대한 레시피 계산 예시."""
    cases = [
        ("Gate oxide", "dry", 0.010, 900, "Massoud"),
        ("Gate oxide", "dry", 0.020, 950, "Massoud"),
        ("Pad / screen oxide", "dry", 0.050, 1000, "Massoud"),
        ("Field oxide", "wet", 0.300, 1000, "Deal-Grove"),
        ("Field oxide", "wet", 0.500, 1050, "Deal-Grove"),
        ("Field oxide", "wet", 1.000, 1100, "Deal-Grove"),
    ]
    out = []
    for name, ambient, x, T, model in cases:
        t_h = massoud_time(x, T) if model == "Massoud" else float(deal_grove_time(x, T, ambient))
        out.append({
            "purpose": name, "ambient": ambient, "target_nm": round(x * 1e3), "T_C": T, "model": model,
            "time_min": round(t_h * 60, 1),
            "deal_grove_time_min": round(float(deal_grove_time(x, T, ambient)) * 60, 1),
            "si_consumed_nm": round(float(si_consumed(x)) * 1e3, 1),
            "crossover_A_half_nm": round(regime_crossover_thickness(T, ambient) * 1e3, 1),
        })
    with open(RES_DIR / "recipe_examples.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)
    print("  saved results/oxidation/recipe_examples.csv")
    return out


def rate_table():
    out = []
    for ambient in ("dry", "wet"):
        for T in TEMPS_C:
            B, BA, A = deal_grove_rates(T, ambient, "100")
            out.append({"ambient": ambient, "T_C": T, "B_um2_h": f"{B:.4g}", "B_over_A_um_h": f"{BA:.4g}",
                        "A_um": f"{A:.4g}", "crossover_A_half_nm": round(A / 2 * 1e3, 1)})
    with open(RES_DIR / "rate_constants_100.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)
    print("  saved results/oxidation/rate_constants_100.csv")


def main():
    print("Oxidation module")
    fig1_thickness_vs_time()
    fig2_arrhenius()
    fig3_universal_curve()
    fig4_recipe_map()
    fig5_thin_oxide()
    validation()
    recipe_examples()
    rate_table()


if __name__ == "__main__":
    main()
