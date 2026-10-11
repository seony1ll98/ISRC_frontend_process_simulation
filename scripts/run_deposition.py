"""증착 및 금속배선 모듈의 그래프와 결과 표를 모두 다시 만든다.

실행: python scripts/run_deposition.py
출력: figures/deposition/*.png, results/deposition/*.csv
"""

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from plot_style import BLUE_RAMP, CAT, GRID, SURFACE, TEXT, TEXT_2  # noqa: E402
from plot_style import save as _save  # noqa: E402

from isrc_sim.constants import celsius_to_kelvin  # noqa: E402
from isrc_sim.deposition import (  # noqa: E402
    direct_flux,
    grove_growth_rate,
    lpcvd_hG_gain,
    solve_trench_flux,
    step_coverage_metrics,
    surface_rate_constant,
    temperature_sensitivity,
    transition_temperature_shift,
    trench_elements,
    view_factor_matrix,
    view_factor_to_opening,
)
from isrc_sim.metallization import (  # noqa: E402
    INFINITE_SHEET_CF,
    RESISTIVITY_UOHM_CM,
    current_density_MA_cm2,
    four_point_probe_cf_disk,
    four_point_probe_fd,
    rc_delay,
    resistivity_from_sheet,
    sheet_from_resistivity,
    smits_center_cf,
)

FIG_DIR = ROOT / "figures" / "deposition"
RES_DIR = ROOT / "results" / "deposition"

EA = 1.6                     # Plummer: Si 에피 증착의 표면 반응 활성화 에너지 [eV]
T_TRANSITION_AP = 1100.0     # 상압 CVD에서 ks = hG가 되는 온도 (시나리오 가정) [C]
HG_GAIN = lpcvd_hG_gain(1.0)

# Plummer Ch. 9 표: 증착 방법별 입사각 지수 n과 부착 확률 S_c (범위의 대표값)
METHODS = [
    ("Sputter (n = 1, S$_c$ = 1)", "sputter", 1, 1.0),
    ("Collimated PVD (n = 20, S$_c$ = 1)", "collimated", 20, 1.0),
    ("LPCVD SiO$_2$, silane (S$_c$ = 0.3)", "silane", 1, 0.3),
    ("LPCVD SiO$_2$, TEOS (S$_c$ = 0.075)", "TEOS", 1, 0.075),
    ("LPCVD W (S$_c$ = 0.01)", "W", 1, 0.01),
]


def save(fig, name):
    _save(fig, FIG_DIR, name)


def write_csv(name, rows):
    RES_DIR.mkdir(parents=True, exist_ok=True)
    with open(RES_DIR / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"  saved results/deposition/{name}")


def fig1_cvd():
    T = np.linspace(500, 1300, 300)
    inv = 1000 / celsius_to_kelvin(T)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    ax = axes[0]
    ks = surface_rate_constant(T, EA, T_TRANSITION_AP, hG=1.0)
    ax.semilogy(inv, ks, color=TEXT_2, ls=":", lw=1.3, label="k$_s$ (surface reaction)")
    for hG, c, lab in ((1.0, CAT[1], "APCVD, 760 torr"), (HG_GAIN, CAT[0], "LPCVD, 1 torr")):
        ax.semilogy(inv, grove_growth_rate(ks, hG), color=c, label=f"{lab}: growth rate")
        ax.axhline(hG, color=c, ls="--", lw=1)
    T_lp = transition_temperature_shift(T_TRANSITION_AP, EA, HG_GAIN)
    ax.text(0.70, 1.4, "h$_G$ (760 torr)", color=CAT[1], fontsize=8)
    ax.text(0.70, HG_GAIN * 1.4, f"h$_G$ (1 torr) ≈ {HG_GAIN:.0f}×", color=CAT[0], fontsize=8)
    ax.set(title="Grove model: reaction- vs. transport-limited", xlabel="1000 / T (1/K)",
           ylabel="Growth rate (normalized to h$_G$ at 760 torr)", ylim=(1e-4, 1e3))
    top = ax.secondary_xaxis("top", functions=(lambda v: 1000.0 / np.asarray(v) - 273.15,
                                               lambda v: 1000.0 / (np.asarray(v) + 273.15)))
    top.set_xlabel("T (°C)", color=TEXT_2)
    ax.legend(loc="lower left", fontsize=8)

    # (b) 온도와 가스 흐름(h_G)에 대한 민감도
    ax = axes[1]
    rows = []
    for hG, c, lab in ((1.0, CAT[1], "APCVD"), (HG_GAIN, CAT[0], "LPCVD")):
        k = surface_rate_constant(T, EA, T_TRANSITION_AP, hG=1.0)
        frac_rxn = hG / (k + hG)                    # 반응 율속 비중
        s_T = frac_rxn * temperature_sensitivity(T, EA) * 100   # %/K
        s_h = k / (k + hG) * 10                     # h_G 10 % 변화에 대한 % 변화
        ax.plot(T, s_T, color=c, label=f"{lab}: per +1 °C")
        ax.plot(T, s_h, color=c, ls="--", label=f"{lab}: per +10 % h$_G$ (gas flow)")
        for Ti in (700, 900, 1100):
            i = int(np.argmin(np.abs(T - Ti)))
            rows.append({"system": lab, "T_C": Ti, "rate_change_pct_per_K": round(float(s_T[i]), 2),
                         "rate_change_pct_per_10pct_hG": round(float(s_h[i]), 2)})
    ax.set(title="What the wafer-to-wafer uniformity depends on", xlabel="Temperature (°C)",
           ylabel="Growth-rate change (%)", ylim=(0, 10.5))
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    save(fig, "fig1_cvd_regimes.png")
    rows.append({"system": "transition", "T_C": round(T_lp, 1), "rate_change_pct_per_K": "LPCVD transition T",
                 "rate_change_pct_per_10pct_hG": f"hG gain {HG_GAIN:.1f}"})
    write_csv("cvd_sensitivity.csv", rows)
    return T_lp


def draw_film(ax, el, rate, width, depth, scale, color):
    """트렌치 단면 위에 증착 두께(상대 속도 x scale)를 법선 방향으로 그린다."""
    w2 = width / 2
    ax.fill([-1.2 * width, -w2, -w2, w2, w2, 1.2 * width, 1.2 * width, -1.2 * width],
            [0, 0, -depth, -depth, 0, 0, -depth - 0.3, -depth - 0.3], color="#d9d7d0", lw=0)
    outer = el["mid"] + el["normal"] * (rate[:, None] * scale)
    kind = el["kind"]
    poly = np.vstack([[-1.2 * width, scale], [-w2 - 0 * scale, scale], [-w2 + rate[kind == "left"][0] * scale, 0],
                      outer, [w2 - rate[kind == "right"][-1] * scale, 0], [w2, scale], [1.2 * width, scale],
                      [1.2 * width, 0], [w2, 0], [w2, -depth], [-w2, -depth], [-w2, 0], [-1.2 * width, 0]])
    ax.add_patch(Polygon(poly, closed=True, color=color, lw=0, alpha=0.85))
    ax.set_xlim(-1.2 * width, 1.2 * width)
    ax.set_ylim(-depth - 0.3, scale + 0.15)
    ax.set_aspect("equal")
    ax.axis("off")


def fig2_step_coverage(width=1.0, depth=2.0):
    el = trench_elements(width, depth, n_per_side=60)
    V = view_factor_matrix(el)
    fig = plt.figure(figsize=(15, 6.2))
    gs = fig.add_gridspec(2, 5, height_ratios=[1.25, 1])
    rows = []
    ax_line = fig.add_subplot(gs[1, :])
    # 둘레를 따라 펼친 좌표: 왼쪽 벽(위→아래), 바닥, 오른쪽 벽(아래→위)
    s = np.concatenate([[0], np.cumsum(el["length"])])[:-1] + el["length"] / 2
    for k, ((label, key, n, Sc), c) in enumerate(zip(METHODS, [CAT[1], CAT[3], BLUE_RAMP[1], BLUE_RAMP[2], BLUE_RAMP[3]])):
        F, _, _ = solve_trench_flux(el, width, n=n, Sc=Sc, V=V)
        m = step_coverage_metrics(el, F)
        ax = fig.add_subplot(gs[0, k])
        draw_film(ax, el, F, width, depth, 0.3, c)
        ax.set_title(f"{label}\nbottom {m['bottom_center']:.2f}, wall {m['sidewall_mid']:.2f}", fontsize=9)
        ax_line.plot(s, F, color=c, label=label)
        rows.append({"method": key, "n": n, "Sc": Sc, "aspect_ratio": depth / width,
                     **{k2: round(v, 3) for k2, v in m.items()}})
    L = el["length"]
    for x in (depth, depth + width):
        ax_line.axvline(x, color=GRID, lw=1.2)
    ax_line.text(depth / 2, 1.03, "left wall (top → bottom)", ha="center", color=TEXT_2, fontsize=8.5)
    ax_line.text(depth + width / 2, 1.03, "bottom", ha="center", color=TEXT_2, fontsize=8.5)
    ax_line.text(depth + width + depth / 2, 1.03, "right wall (bottom → top)", ha="center", color=TEXT_2, fontsize=8.5)
    ax_line.set(xlabel="Position along the trench surface (µm)", ylabel="Deposition rate / flat field",
                xlim=(0, L.sum()), ylim=(0, 1.1))
    ax_line.legend(loc="lower center", ncol=5, fontsize=8, bbox_to_anchor=(0.5, -0.42))
    fig.suptitle(f"Step coverage in a trench, aspect ratio {depth / width:g} (film drawn 0.3 × field thickness)",
                 color=TEXT, fontweight="bold")
    fig.tight_layout()
    save(fig, "fig2_step_coverage.png")
    write_csv("step_coverage_by_method.csv", rows)


def fig3_coverage_maps():
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    Scs = np.logspace(-3, 0, 25)
    rows = []
    ax = axes[0]
    for ar, c in zip((1, 2, 4), BLUE_RAMP[1:]):
        el = trench_elements(1.0, float(ar), n_per_side=40 * ar)
        V = view_factor_matrix(el)
        mins = [solve_trench_flux(el, 1.0, n=1, Sc=s, V=V)[0].min() for s in Scs]
        ax.semilogx(Scs, mins, color=c, label=f"aspect ratio {ar}")
        for s_ in (1.0, 0.3, 0.075, 0.01, 0.001):
            rows.append({"aspect_ratio": ar, "Sc": s_, "n": 1,
                         "min_coverage": round(float(solve_trench_flux(el, 1.0, n=1, Sc=s_, V=V)[0].min()), 3)})
    for lo, hi, lab in ((0.2, 0.4, "silane"), (0.05, 0.1, "TEOS"), (0.0005, 0.01, "W, poly")):
        ax.axvspan(lo, hi, color=GRID, lw=0)
        ax.text(np.sqrt(lo * hi), 0.03, lab, ha="center", color=TEXT_2, fontsize=8)
    ax.set(title="Isotropic arrival (n = 1): minimum coverage vs. S$_c$", xlabel="Sticking coefficient S$_c$",
           ylabel="Thinnest point / flat field", ylim=(0, 1.05))
    ax.legend(loc="center left", fontsize=8.5)

    ax = axes[1]
    ns = [1, 2, 4, 8, 20, 40, 80]
    el = trench_elements(1.0, 2.0, n_per_side=80)
    V = view_factor_matrix(el)
    bot, wall = [], []
    for n in ns:
        m = step_coverage_metrics(el, solve_trench_flux(el, 1.0, n=n, Sc=1.0, V=V)[0])
        bot.append(m["bottom_center"])
        wall.append(m["sidewall_mid"])
        rows.append({"aspect_ratio": 2, "Sc": 1.0, "n": n, "min_coverage": round(m["min_over_trench"], 3)})
    ax.semilogx(ns, bot, marker="o", color=CAT[0], label="bottom center")
    ax.semilogx(ns, wall, marker="s", color=CAT[1], label="sidewall middle")
    ax.axvspan(1, 4, color=GRID, lw=0)
    ax.text(2, 0.55, "sputter", ha="center", color=TEXT_2, fontsize=8)
    ax.axvspan(8, 80, color=GRID, lw=0, alpha=0.6)
    ax.text(25, 0.55, "collimated /\nionized PVD", ha="center", color=TEXT_2, fontsize=8)
    ax.set(title="PVD (S$_c$ = 1), aspect ratio 2: effect of arrival angle", xlabel="Arrival exponent n (cos$^n$θ)",
           ylabel="Deposition rate / flat field", ylim=(0, 1.05))
    ax.set_xticks(ns, [str(v) for v in ns])
    ax.minorticks_off()
    ax.legend(loc="upper left", fontsize=8.5)
    fig.tight_layout()
    save(fig, "fig3_coverage_maps.png")
    write_csv("coverage_maps.csv", rows)


def fig4_four_point_probe():
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    ax = axes[0]
    r = np.logspace(np.log10(4), 3, 200)
    ax.semilogx(r, [four_point_probe_cf_disk(v, 1.0) for v in r], color=CAT[0], label="image method (this project)")
    ax.semilogx(r[::12], [smits_center_cf(v, 1.0) for v in r[::12]], ls="none", marker="o", ms=5, mfc=SURFACE,
                mec=CAT[1], mew=1.3, label="Smits (1958) closed form")
    rr = [5, 10, 20, 40]
    ax.semilogx(rr, [four_point_probe_fd(v, 1.0, 0, 20) for v in rr], ls="none", marker="s", ms=6, mfc=SURFACE,
                mec=CAT[2], mew=1.4, label="finite difference")
    ax.axhline(INFINITE_SHEET_CF, color=TEXT_2, ls=":", lw=1.2)
    ax.text(900, INFINITE_SHEET_CF + 0.04, "π / ln 2 = 4.532", ha="right", color=TEXT_2, fontsize=8)
    ax.set(title="Correction factor at the wafer center", xlabel="Sample diameter / probe spacing (D/s)",
           ylabel="R$_s$ / (V/I)", ylim=(3.2, 4.7))
    ax.legend(loc="lower right", fontsize=8.5)

    ax = axes[1]
    rows = []
    for D, c, lab in ((100.0, CAT[0], "100 mm wafer"), (20.0, CAT[1], "20 mm coupon")):
        x0 = np.linspace(0, D / 2 - 1.6, 200)
        err = [(INFINITE_SHEET_CF / four_point_probe_cf_disk(D, 1.0, v) - 1) * 100 for v in x0]
        ax.plot(x0 / (D / 2) * 100, err, color=c, label=f"{lab}, s = 1 mm")
        for frac in (0.0, 0.5, 0.8, 0.9, 0.95):
            xv = frac * (D / 2 - 1.6)
            rows.append({"sample": lab, "probe_center_from_middle_mm": round(xv, 2),
                         "CF": round(four_point_probe_cf_disk(D, 1.0, xv), 4),
                         "error_if_4.532_used_pct": round((INFINITE_SHEET_CF / four_point_probe_cf_disk(D, 1.0, xv) - 1) * 100, 2)})
    ax.axhspan(-1, 1, color=GRID, lw=0)
    ax.set(title="Overestimate of R$_s$ if 4.532 is used everywhere", xlabel="Probe center position (% of radius)",
           ylabel="Error in R$_s$ (%)", xlim=(0, 100))
    ax.legend(loc="upper left", fontsize=8.5)
    fig.tight_layout()
    save(fig, "fig4_four_point_probe.png")
    write_csv("four_point_probe_correction.csv", rows)


def fig5_interconnect():
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    ax = axes[0]
    Lmm = np.logspace(-1, 1.3, 100)
    geom = dict(W_um=1.0, H_um=1.0, x_ox_um=1.0, L_s_um=1.0)
    colors = {"poly-Si": "#86b6ef", "WSi2": "#5598e7", "W": "#2a78d6", "Al": CAT[1], "Cu": CAT[2]}
    for mat, rho in RESISTIVITY_UOHM_CM.items():
        ax.loglog(Lmm, rc_delay(rho, Lmm, K_ox=3.9, **geom), color=colors[mat], label=f"{mat} ({rho:g} µΩ·cm), SiO$_2$")
    ax.loglog(Lmm, rc_delay(2.2, Lmm, K_ox=2.7, **geom), color=CAT[2], ls="--", label="Cu, low-k (K = 2.7)")
    ax.axhspan(1e-11, 1e-10, color=GRID, lw=0)
    ax.text(0.12, 1.4e-11, "10~100 ps", color=TEXT_2, fontsize=8)
    ax.set(title="Line delay τ = 0.89 RC (1 µm lines, K$_I$ = 1)", xlabel="Line length (mm)", ylabel="Delay (s)")
    ax.legend(loc="upper left", fontsize=7.5)

    ax = axes[1]
    cases = [("Al / SiO$_2$", 3.0, 3.9), ("Cu / SiO$_2$", 2.2, 3.9), ("Cu / low-k 2.7", 2.2, 2.7)]
    taus = [rc_delay(r, 10.0, K_ox=k, **geom) * 1e12 for _, r, k in cases]
    ax.bar(range(3), taus, color=[CAT[1], CAT[2], "#0e7d54"], width=0.55)
    for i, t in enumerate(taus):
        ax.text(i, t + 3, f"{t:.0f} ps ({t / taus[0] * 100:.0f} %)", ha="center", fontsize=9, color=TEXT)
    ax.set_xticks(range(3), [c[0] for c in cases])
    ax.set(title="10 mm global line", ylabel="Delay (ps)", ylim=(0, max(taus) * 1.2))
    ax.grid(axis="x", visible=False)
    fig.tight_layout()
    save(fig, "fig5_interconnect_rc.png")
    rows = [{"stack": c[0].replace("$_2$", "2"), "rho_uohm_cm": c[1], "K": c[2], "delay_ps_10mm": round(t, 1),
             "relative_to_Al_SiO2": round(t / taus[0], 3)} for c, t in zip(cases, taus)]
    for J in (0.1, 0.5):
        rows.append({"stack": f"max current at {J} MA/cm2 (1x1 um line)", "rho_uohm_cm": "", "K": "",
                     "delay_ps_10mm": f"{J * 1e6 * 1e-8 * 1e3:.1f} mA", "relative_to_Al_SiO2": ""})
    write_csv("interconnect_rc.csv", rows)


def validation():
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    rows = []
    w, h = 1.0, 2.0
    el = trench_elements(w, h, n_per_side=60)
    V = view_factor_matrix(el)
    vo = view_factor_to_opening(el, w)
    s = np.concatenate([[0], np.cumsum(el["length"])])[:-1] + el["length"] / 2

    ax = axes[0]
    fd = direct_flux(el, w, n=1)
    ax.plot(s, vo, color=CAT[0], label="view factor to opening (crossed strings)")
    ax.plot(s[::5], fd[::5], ls="none", marker="o", ms=5, mfc=SURFACE, mec=CAT[1], mew=1.3,
            label="direct flux, angular integration")
    err = np.max(np.abs(fd - vo))
    ax.set(title=f"S$_c$ = 1, n = 1: max difference {err:.1e}", xlabel="Position along trench surface (µm)",
           ylabel="Flux / flat field", ylim=(0, 0.6))
    ax.legend(loc="upper center", fontsize=8)
    rows.append({"case": "direct flux (n=1) vs view factor, AR 2", "value": f"{err:.1e}"})
    rows.append({"case": "enclosure sum rule max error", "value": f"{np.max(np.abs(V.sum(1) + vo - 1)):.1e}"})

    ax = axes[1]
    Scs = np.logspace(-4, 0, 17)
    bal, mins = [], []
    for Sc in Scs:
        F, _, _ = solve_trench_flux(el, w, n=1, Sc=Sc, V=V)
        stuck = np.sum(Sc * F * el["length"])
        esc = np.sum((1 - Sc) * F * el["length"] * vo)
        bal.append(abs(stuck + esc - w))
        mins.append(F.min())
    ax.semilogx(Scs, mins, color=CAT[0], marker="o", ms=4, label="minimum coverage → 1 as S$_c$ → 0")
    ax.set(title=f"Limits and particle balance (max balance error {max(bal):.0e})", xlabel="Sticking coefficient S$_c$",
           ylabel="Thinnest point / flat field", ylim=(0, 1.05))
    ax.legend(loc="lower left", fontsize=8.5)
    rows.append({"case": "particle balance max error (n=1, Sc 1e-4..1)", "value": f"{max(bal):.1e}"})
    rows.append({"case": "min coverage at Sc = 1e-4", "value": f"{mins[0]:.4f}"})

    ax = axes[2]
    ks = [4, 10, 20, 40]
    for D, c in ((10.0, CAT[0]), (20.0, CAT[1])):
        ref = four_point_probe_cf_disk(D, 1.0)
        errs = [abs(four_point_probe_fd(D, 1.0, 0, k) / ref - 1) * 100 for k in ks]
        ax.loglog(ks, errs, marker="o", color=c, label=f"D/s = {D:g}")
        rows.append({"case": f"4PP FD vs image, D/s={D:g}, nodes per s {ks}", "value": " / ".join(f"{e:.3f}%" for e in errs)})
    ax.set_xticks(ks, [str(k) for k in ks])
    ax.minorticks_off()
    ax.set(title="4-point probe: finite difference vs. image method", xlabel="Grid points per probe spacing",
           ylabel="Difference (%)")
    ax.legend(loc="upper right", fontsize=8.5)
    fig.tight_layout()
    save(fig, "fig6_validation.png")
    write_csv("validation.csv", rows)


def main():
    print("Deposition & metallization module")
    T_lp = fig1_cvd()
    print(f"  LPCVD transition temperature: {T_lp:.0f} C (APCVD {T_TRANSITION_AP:.0f} C)")
    fig2_step_coverage()
    fig3_coverage_maps()
    fig4_four_point_probe()
    fig5_interconnect()
    validation()
    # ISRC 실습과 연결: 스퍼터 Al 박막의 면저항 예시
    for t in (300, 500, 1000):
        print(f"  Al {t} nm: Rs = {sheet_from_resistivity(3.0, t) * 1e3:.1f} mOhm/sq "
              f"(check rho = {resistivity_from_sheet(sheet_from_resistivity(3.0, t), t):.1f} uOhm cm)")
    print(f"  0.5 MA/cm2 in 1x1 um line: {0.5 * 1e6 * 1e-8 * 1e3:.1f} mA, check {current_density_MA_cm2(5, 1, 1):.2f}")


if __name__ == "__main__":
    main()
