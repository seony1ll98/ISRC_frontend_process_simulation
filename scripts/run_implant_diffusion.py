"""이온주입 및 확산 모듈의 그래프와 결과 표를 모두 다시 만든다.

실행: python scripts/run_implant_diffusion.py
출력: figures/implant_diffusion/*.png, results/implant_diffusion/*.csv
"""

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from plot_style import BLUE_RAMP_6, CAT, SURFACE, TEXT, TEXT_2  # noqa: E402
from plot_style import save as _save  # noqa: E402

from isrc_sim.constants import IMPLANT_ENERGY_KEV, IMPLANT_RANGE_UM, MASETTI, celsius_to_kelvin  # noqa: E402
from isrc_sim.diffusion import (  # noqa: E402
    annealed_gaussian,
    diffuse_fd,
    dose_of,
    gaussian_drive_in_Dt_for_xj,
    intrinsic_diffusivity,
    junction_depth,
    predeposition_dose,
    predeposition_erfc,
    surface_fraction_lost,
    thermal_budget,
)
from isrc_sim.electrical import mobility, sheet_resistance  # noqa: E402
from isrc_sim.implant import implant_profile, peak_concentration, range_parameters  # noqa: E402

FIG_DIR = ROOT / "figures" / "implant_diffusion"
RES_DIR = ROOT / "results" / "implant_diffusion"
BYU_CSV = ROOT / "data" / "reference" / "byu_resistivity_2026-10-11.csv"

DOPANTS = ["B", "P", "As", "Sb"]
DOPANT_COLOR = dict(zip(DOPANTS, CAT))

# 시나리오: p형 기판(B 1e15)에 P를 주입해 n-well을 만든다.
NB = 1e15
WELL = dict(dopant="P", energy_keV=100, dose=1e13, T_C=1100)
WELL_TIMES_H = [0, 0.5, 1, 2, 4, 8]


def save(fig, name):
    _save(fig, FIG_DIR, name)


def write_csv(name, rows):
    RES_DIR.mkdir(parents=True, exist_ok=True)
    with open(RES_DIR / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"  saved results/implant_diffusion/{name}")


def fig1_range_table():
    E = np.linspace(10, 200, 200)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for d in DOPANTS:
        table = np.array(IMPLANT_RANGE_UM[d])
        Rp, dRp = range_parameters(d, E)
        for ax, y, yt in ((axes[0], Rp, table[:, 0]), (axes[1], dRp, table[:, 1])):
            ax.plot(E, y * 1e3, color=DOPANT_COLOR[d], label=d)
            ax.plot(IMPLANT_ENERGY_KEV, yt * 1e3, ls="none", marker="o", ms=4.5,
                    mfc=SURFACE, mec=DOPANT_COLOR[d], mew=1.3)
    axes[0].set(title="Projected range R$_p$", xlabel="Implant energy (keV)", ylabel="R$_p$ (nm)")
    axes[1].set(title="Straggle ΔR$_p$", xlabel="Implant energy (keV)", ylabel="ΔR$_p$ (nm)")
    for ax in axes:
        ax.set_xlim(0, 210)
        ax.set_ylim(bottom=0)
        ax.legend(loc="upper left")
    fig.suptitle("Range statistics in Si (markers: Plummer table, lines: PCHIP interpolation)",
                 color=TEXT, fontweight="bold")
    fig.tight_layout()
    save(fig, "fig1_range_table.png")


def fig2_implant_profiles():
    x = np.linspace(0, 0.8, 1600)
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    rows = []
    for d in ("B", "P", "As"):
        C = implant_profile(x, d, 100, 1e15)
        Rp, dRp = range_parameters(d, 100)
        ax.semilogy(x * 1e3, C, color=DOPANT_COLOR[d], label=f"{d}: R$_p$ = {Rp * 1e3:.0f} nm")
        rows.append({"dopant": d, "energy_keV": 100, "dose_cm2": "1e15", "Rp_nm": round(Rp * 1e3, 1),
                     "dRp_nm": round(dRp * 1e3, 1), "peak_cm3": f"{peak_concentration(1e15, dRp):.3g}",
                     "fraction_above_surface": f"{surface_fraction_lost(Rp, dRp):.2e}"})
    ax.set(title="As-implanted profiles, 100 keV, 10$^{15}$ cm$^{-2}$", xlabel="Depth (nm)",
           ylabel="Concentration (cm$^{-3}$)", ylim=(1e16, 3e21), xlim=(0, 800))
    ax.legend(loc="upper right")
    fig.tight_layout()
    save(fig, "fig2_implant_profiles.png")
    write_csv("implant_100keV.csv", rows)


def well_profiles(x):
    Rp, dRp = range_parameters(WELL["dopant"], WELL["energy_keV"])
    out = []
    for t in WELL_TIMES_H:
        Dt = thermal_budget(WELL["dopant"], [(WELL["T_C"], t)])
        out.append((t, Dt, annealed_gaussian(x, WELL["dose"], Rp, dRp, Dt)))
    return Rp, dRp, out


def fig3_well_anneal():
    x = np.linspace(0, 5, 2501)
    Rp, dRp, profiles = well_profiles(x)
    fig, ax = plt.subplots(figsize=(7, 4.6))
    for (t, Dt, C), c in zip(profiles, BLUE_RAMP_6):
        xj = junction_depth(x, C, NB)
        ax.semilogy(x, C, color=c, label=f"{t:g} h  (x$_j$ = {xj:.2f} µm)")
    # 수치해를 같은 그래프에 겹쳐 해석해와 일치하는지 보인다.
    C0 = annealed_gaussian(x, WELL["dose"], Rp, dRp, 0.0)
    for t in (1, 8):
        Dt = thermal_budget(WELL["dopant"], [(WELL["T_C"], t)])
        C_fd = diffuse_fd(x, C0, Dt, n_steps=400)
        ax.semilogy(x[::60], C_fd[::60], ls="none", marker="o", ms=5, mfc=SURFACE, mec=TEXT_2, mew=1.2)
    ax.plot([], [], ls="none", marker="o", ms=5, mfc=SURFACE, mec=TEXT_2, label="numerical (FD), 1 h and 8 h")
    ax.axhline(NB, color=TEXT_2, ls="--", lw=1.2)
    ax.text(4.95, NB * 1.3, "p-substrate N$_B$ = 10$^{15}$", ha="right", color=TEXT_2, fontsize=8.5)
    ax.set(title=f"N-well drive-in: P 100 keV, 10$^{{13}}$ cm$^{{-2}}$, {WELL['T_C']} °C",
           xlabel="Depth (µm)", ylabel="Concentration (cm$^{-3}$)", ylim=(1e14, 2e18), xlim=(0, 5))
    ax.legend(title="anneal time", loc="upper right", fontsize=8.5)
    fig.tight_layout()
    save(fig, "fig3_well_anneal.png")


def fig4_well_design():
    x = np.linspace(0, 8, 8001)
    Rp, dRp = range_parameters(WELL["dopant"], WELL["energy_keV"])
    t = np.linspace(0, 10, 81)
    xj, Rs, peak = [], [], []
    for ti in t:
        Dt = thermal_budget(WELL["dopant"], [(WELL["T_C"], ti)])
        C = annealed_gaussian(x, WELL["dose"], Rp, dRp, Dt)
        xj.append(junction_depth(x, C, NB))
        Rs.append(sheet_resistance(x, C, WELL["dopant"], NB))
        peak.append(C.max())
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].plot(t, xj, color=CAT[0])
    axes[0].set(title="Junction depth", xlabel="Drive-in time at 1100 °C (h)", ylabel="x$_j$ (µm)")
    axes[0].set_ylim(bottom=0)
    axes[1].plot(t, Rs, color=CAT[0])
    axes[1].set(title="Sheet resistance of the well", xlabel="Drive-in time at 1100 °C (h)",
                ylabel="R$_s$ (Ω/sq)")
    axes[1].set_ylim(bottom=0)
    for ax, y in ((axes[0], xj), (axes[1], Rs)):
        for tm in (2, 8):
            i = int(np.argmin(np.abs(t - tm)))
            ax.plot(t[i], y[i], marker="o", ms=6, color=CAT[0], mec=SURFACE, mew=1.5)
            ax.annotate(f"{y[i]:.2f}" if ax is axes[0] else f"{y[i]:.0f}", (t[i], y[i]),
                        xytext=(6, -12), textcoords="offset points", fontsize=8.5, color=TEXT_2)
    fig.suptitle("N-well design curves (same dose, longer anneal)", color=TEXT, fontweight="bold")
    fig.tight_layout()
    save(fig, "fig4_well_design.png")

    rows = []
    for ti in WELL_TIMES_H:
        Dt = thermal_budget(WELL["dopant"], [(WELL["T_C"], ti)])
        C = annealed_gaussian(x, WELL["dose"], Rp, dRp, Dt)
        rows.append({"time_h": ti, "Dt_cm2": f"{Dt:.3g}", "diffusion_length_2sqrtDt_um": round(2 * np.sqrt(Dt) * 1e4, 3),
                     "peak_cm3": f"{C.max():.3g}", "xj_um": round(junction_depth(x, C, NB), 3),
                     "Rs_ohm_sq": round(sheet_resistance(x, C, WELL["dopant"], NB))})
    write_csv("nwell_drive_in.csv", rows)


def fig5_dopant_comparison():
    T = np.linspace(900, 1200, 100)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for d in DOPANTS:
        axes[0].semilogy(1000 / celsius_to_kelvin(T), intrinsic_diffusivity(d, T), color=DOPANT_COLOR[d], label=d)
    axes[0].set_title("Intrinsic diffusivity", pad=28)
    axes[0].set(xlabel="1000 / T (1/K)", ylabel="D (cm²/s)")
    axes[0].legend(loc="upper right")
    top = axes[0].secondary_xaxis("top", functions=(lambda v: 1000.0 / np.asarray(v) - 273.15,
                                                    lambda v: 1000.0 / (np.asarray(v) + 273.15)))
    top.set_xlabel("T (°C)", color=TEXT_2)

    x = np.linspace(0, 6, 6001)
    rows = []
    xj0, xj1 = [], []
    for d in DOPANTS:
        Rp, dRp = range_parameters(d, 100)
        Dt = thermal_budget(d, [(1100, 2)])
        a = junction_depth(x, annealed_gaussian(x, 1e13, Rp, dRp, 0.0), NB)
        b = junction_depth(x, annealed_gaussian(x, 1e13, Rp, dRp, Dt), NB)
        xj0.append(a)
        xj1.append(b)
        rows.append({"dopant": d, "D_1100C_cm2_s": f"{intrinsic_diffusivity(d, 1100):.3g}", "Rp_um": round(Rp, 4),
                     "xj_as_implanted_um": round(a, 3), "xj_after_1100C_2h_um": round(b, 3)})
    pos = np.arange(len(DOPANTS))
    w = 0.38
    axes[1].bar(pos - w / 2 - 0.01, xj0, w, color="#b7d3f6", label="as-implanted")
    axes[1].bar(pos + w / 2 + 0.01, xj1, w, color=CAT[0], label="after 1100 °C, 2 h")
    for i, (a, b) in enumerate(zip(xj0, xj1)):
        axes[1].text(i + w / 2 + 0.01, b + 0.03, f"{b:.2f}", ha="center", fontsize=8.5, color=TEXT_2)
    axes[1].set_xticks(pos, DOPANTS)
    axes[1].set(title="Junction depth, 100 keV, 10$^{13}$ cm$^{-2}$, N$_B$ = 10$^{15}$",
                ylabel="x$_j$ (µm)")
    axes[1].grid(axis="x", visible=False)
    axes[1].legend(loc="upper right")
    fig.tight_layout()
    save(fig, "fig5_dopant_comparison.png")
    write_csv("dopant_comparison.csv", rows)


def validation():
    rows = []
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))

    # (a) 주입 Gaussian 열처리: 수치해 vs 해석해
    x = np.linspace(0, 4, 2001)
    Rp, dRp = range_parameters("B", 100)
    C0 = annealed_gaussian(x, 1e13, Rp, dRp, 0.0)
    Dt = thermal_budget("B", [(1100, 2.0)])
    C_an = annealed_gaussian(x, 1e13, Rp, dRp, Dt)
    C_fd = diffuse_fd(x, C0, Dt, n_steps=400)
    m = C_an > 1e-3 * C_an.max()
    err_a = np.max(np.abs(C_fd[m] - C_an[m]) / C_an[m]) * 100
    ax = axes[0]
    ax.semilogy(x, C_an, color=CAT[0], label="analytic (Gaussian + image)")
    ax.semilogy(x[::50], C_fd[::50], ls="none", marker="o", ms=5, mfc=SURFACE, mec=CAT[1], mew=1.3,
                label="numerical (Crank-Nicolson)")
    ax.set(title=f"B implant + 1100 °C 2 h: max error {err_a:.3f} %", xlabel="Depth (µm)",
           ylabel="Concentration (cm$^{-3}$)", ylim=(1e14, 1e18), xlim=(0, 2.5))
    ax.legend(loc="upper right")
    rows.append({"case": "B 100 keV 1e13 + 1100C 2h (reflecting surface)", "max_rel_error_pct": round(err_a, 3),
                 "dose_error_pct": round((dose_of(x, C_fd) / 1e13 - 1) * 100, 4)})

    # (b) predeposition: 수치해 vs erfc
    x = np.linspace(0, 2, 2001)
    Dt = 1e-9
    C_an = predeposition_erfc(x, 1e20, Dt)
    C_fd = diffuse_fd(x, np.zeros_like(x), Dt, n_steps=2000, surface="fixed", Cs=1e20)
    m = C_an > 1e-3 * 1e20
    err_b = np.max(np.abs(C_fd[m] - C_an[m]) / C_an[m]) * 100
    ax = axes[1]
    ax.semilogy(x, C_an, color=CAT[0], label="analytic (erfc)")
    ax.semilogy(x[::40], C_fd[::40], ls="none", marker="o", ms=5, mfc=SURFACE, mec=CAT[1], mew=1.3,
                label="numerical (Crank-Nicolson)")
    ax.set(title=f"Constant-source predep: max error {err_b:.3f} %", xlabel="Depth (µm)",
           ylabel="Concentration (cm$^{-3}$)", ylim=(1e16, 3e20), xlim=(0, 1.5))
    ax.legend(loc="upper right")
    rows.append({"case": "predeposition Cs=1e20, Dt=1e-9 (fixed surface)", "max_rel_error_pct": round(err_b, 3),
                 "dose_error_pct": round((dose_of(x, C_fd) / predeposition_dose(1e20, Dt) - 1) * 100, 4)})

    # (c) 이동도: Masetti 구현 vs BYU 계산기
    ax = axes[2]
    N = np.logspace(14, 20, 200)
    byu = list(csv.DictReader(open(BYU_CSV, newline="")))
    mob_rows = []
    for d, c in (("P", CAT[0]), ("As", CAT[2]), ("B", CAT[1])):
        ax.semilogx(N, mobility(N, d), color=c, label=f"{d} (this project)")
        pts = [r for r in byu if r["dopant"] == d]
        ax.semilogx([float(r["N_cm3"]) for r in pts], [float(r["mobility_cm2_Vs"]) for r in pts],
                    ls="none", marker="D", ms=5.5, mfc=SURFACE, mec=c, mew=1.4)
        for r in pts:
            ours = float(mobility(float(r["N_cm3"]), d))
            mob_rows.append({"dopant": d, "N_cm3": r["N_cm3"], "byu_cm2_Vs": r["mobility_cm2_Vs"],
                             "this_project_cm2_Vs": round(ours, 3),
                             "diff_pct": round((ours / float(r["mobility_cm2_Vs"]) - 1) * 100, 3)})
    ax.plot([], [], ls="none", marker="D", mfc=SURFACE, mec=TEXT_2, label="BYU calculator")
    ax.annotate("B: BYU uses μ$_{min2}$ = μ$_{min1}$\n(original Masetti: μ$_{min2}$ = 0)", (1e15, 417),
                xytext=(2e14, 620), fontsize=8, color=TEXT_2,
                arrowprops=dict(arrowstyle="-", color=TEXT_2, lw=0.8))
    ax.set(title="Masetti mobility vs. BYU calculator", xlabel="Doping (cm$^{-3}$)",
           ylabel="Mobility (cm²/V·s)", ylim=(0, 1500))
    ax.legend(loc="upper right", fontsize=8.5)
    fig.tight_layout()
    save(fig, "fig6_validation.png")
    write_csv("validation_fd_vs_analytic.csv", rows)
    write_csv("validation_mobility_vs_byu.csv", mob_rows)


def plummer_example():
    """Plummer Ch. 7 예제(B p-well, xj = 3 um, Rs = 900 ohm/sq, NB = 1e15)를 다시 계산한다."""
    xj, NB_ex, Rs_target, T = 3.0, 1e15, 900.0, 1100
    x = np.linspace(0, xj, 30001)

    def Rs_for(Cs):
        Dt = gaussian_drive_in_Dt_for_xj(xj, Cs, NB_ex)
        return sheet_resistance(x, Cs * np.exp(-((x * 1e-4) ** 2) / (4 * Dt)), "B", NB_ex)

    D = intrinsic_diffusivity("B", T)
    Cs_fit = brentq(lambda c: Rs_for(c) - Rs_target, 5e16, 4e17)
    rows = []
    for label, Cs, D_used in (("Plummer (Irvin curve, D from text)", 4e17, 1.5e-13),
                              ("This project, Cs from text", 4e17, D),
                              ("This project, Cs from Masetti integration", Cs_fit, D)):
        Dt = gaussian_drive_in_Dt_for_xj(xj, Cs, NB_ex)
        rows.append({"case": label, "Cs_cm3": f"{Cs:.3g}", "D_cm2_s": f"{D_used:.3g}", "Dt_cm2": f"{Dt:.3g}",
                     "drive_in_h": round(Dt / D_used / 3600, 2), "Rs_ohm_sq_masetti": round(Rs_for(Cs))})
    write_csv("plummer_example.csv", rows)
    return rows


def main():
    print("Implant & diffusion module")
    fig1_range_table()
    fig2_implant_profiles()
    fig3_well_anneal()
    fig4_well_design()
    fig5_dopant_comparison()
    validation()
    plummer_example()


if __name__ == "__main__":
    main()
