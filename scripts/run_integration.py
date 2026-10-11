"""공정 통합 모듈의 그래프와 결과 표를 모두 다시 만든다.

Plummer Ch. 2의 CMOS 공정 흐름을 따라 01~05 모듈을 차례로 적용하고, 문턱전압과 유효 채널 길이를 계산한다.
실행: python scripts/run_integration.py
출력: figures/integration/*.png, results/integration/*.csv
"""

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon, Rectangle
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from plot_style import BLUE_RAMP, BLUE_RAMP_6, CAT, GRID, SURFACE, TEXT, TEXT_2  # noqa: E402
from plot_style import save as _save  # noqa: E402

from isrc_sim.constants import Q_E, SI_CONSUMPTION_RATIO  # noqa: E402
from isrc_sim.deposition import solve_trench_flux, step_coverage_metrics, trench_elements  # noqa: E402
from isrc_sim.diffusion import intrinsic_diffusivity  # noqa: E402
from isrc_sim.electrical import resistivity  # noqa: E402
from isrc_sim.etching import required_selectivity  # noqa: E402
from isrc_sim.implant import implant_profile, peak_concentration, range_parameters  # noqa: E402
from isrc_sim.integration import (  # noqa: E402
    apply_thermal_steps,
    diffuse_2d_fd,
    implanted_and_annealed,
    junction_depth_between,
    ldd_2d_analytic,
    lateral_junction_extent,
    oxide_capacitance,
    vt_nonuniform,
    vt_uniform,
)
from isrc_sim.lithography import depth_of_focus, k1_factor, rayleigh_resolution  # noqa: E402
from isrc_sim.metallization import rc_delay, sheet_from_resistivity  # noqa: E402
from isrc_sim.oxidation import deal_grove_thickness, massoud_thickness, massoud_time  # noqa: E402

FIG_DIR = ROOT / "figures" / "integration"
RES_DIR = ROOT / "results" / "integration"

# ---------------------------------------------------------------------------
# 공정 조건 (Plummer Ch. 2의 대표 범위 안에서 선택)
# ---------------------------------------------------------------------------
SUB_RHO = 10.0                       # p형 기판 5~50 ohm cm
PAD_OX_UM = 0.04                     # 패드 산화막 약 40 nm
LOCOS = (1000, 1.5)                  # 90 min @ 1000 C, H2O
PWELL = ("B", 180, 1e13)             # 1e13 @ 150~200 keV
NWELL = ("P", 200, 1e13)             # 교재 300+ keV, 이 프로젝트의 비정 표는 200 keV까지
DRIVE_IN = [(1050, 5.0)]             # 4~6 h @ 1000~1100 C
GOX_NM, GOX_T = 4.0, 800             # 3~5 nm @ 800 C, O2
POLY_UM, GATE_L_UM = 0.5, 0.5
LDD_N = ("P", 50, 5e13)
LDD_P = ("B", 50, 5e13)
SPACER_DEP_UM = 0.5
SD_N = ("As", 75, 3e15)
SD_P = ("B", 50, 2e15)
ANNEALS = {"furnace 900 °C 30 min": [(900, 0.5)], "RTA 1000 °C 1 min": [(1000, 1 / 60)]}
SILICIDE = [(650, 1 / 60)]
VT_TARGET = 0.30                     # |Vt| 목표 [V] (시나리오)
VT_ADJ_N = ("B", 10)                 # 얕은 주입이 필요해 저에너지 선택 (3.3 참고)
VT_ADJ_P = ("As", 30)

X = np.linspace(0, 6, 3001)
I_BULK = int(np.argmin(np.abs(X - 1.0)))


def save(fig, name):
    _save(fig, FIG_DIR, name)


def write_csv(name, rows):
    RES_DIR.mkdir(parents=True, exist_ok=True)
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(RES_DIR / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"  saved results/integration/{name}")


def substrate_doping():
    return brentq(lambda N: resistivity(N, "B") - SUB_RHO, 1e13, 1e17)


def gate_ox_time(T=GOX_T, tox=GOX_NM):
    return massoud_time(tox * 1e-3, T)


def channel_profiles(anneal, vt_dose_n=0.0, vt_dose_p=0.0, drive=DRIVE_IN, gox_T=GOX_T):
    """NMOS(p-well)와 PMOS(n-well) 채널의 순 도핑 분포. 게이트 산화 시간은 4 nm 기준으로 고정한다."""
    NA = substrate_doping()
    post = [(gox_T, gate_ox_time())] + anneal + SILICIDE
    Bw = implanted_and_annealed(X, PWELL[0], PWELL[1], PWELL[2], drive + post)
    Pw = implanted_and_annealed(X, NWELL[0], NWELL[1], NWELL[2], drive + post)
    Bvt = implanted_and_annealed(X, VT_ADJ_N[0], VT_ADJ_N[1], vt_dose_n, post) if vt_dose_n else 0 * X
    Asvt = implanted_and_annealed(X, VT_ADJ_P[0], VT_ADJ_P[1], vt_dose_p, post) if vt_dose_p else 0 * X
    return {"NA": NA, "Bw": Bw, "Pw": Pw, "Bvt": Bvt, "Asvt": Asvt,
            "nmos": Bw + NA + Bvt, "pmos": Pw - NA + Asvt,
            "nmos_bulk": (Bw + NA)[I_BULK], "pmos_bulk": (Pw - NA)[I_BULK]}


def vts(prof, tox=GOX_NM):
    vn, Wn = vt_nonuniform(X, prof["nmos"], tox, "n", N_bulk=prof["nmos_bulk"], gate="n+")
    vp, Wp = vt_nonuniform(X, prof["pmos"], tox, "p", N_bulk=prof["pmos_bulk"], gate="p+")
    return vn, vp, Wn, Wp


FURNACE = ANNEALS["furnace 900 °C 30 min"]


def design_vt_doses():
    def f_n(d):
        return vts(channel_profiles(FURNACE, vt_dose_n=d))[0] - VT_TARGET

    def f_p(d):
        return vts(channel_profiles(FURNACE, vt_dose_p=d))[1] + VT_TARGET

    return brentq(f_n, 1e11, 5e13), brentq(f_p, 1e10, 5e13)


# ---------------------------------------------------------------------------
# 그래프
# ---------------------------------------------------------------------------
def fig1_thermal_budget():
    stages = [("as implanted", []), ("after well drive-in", DRIVE_IN),
              ("after gate oxidation", DRIVE_IN + [(GOX_T, gate_ox_time())]),
              ("after final anneal (furnace)", DRIVE_IN + [(GOX_T, gate_ox_time())] + FURNACE + SILICIDE)]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    ax = axes[0]
    for (lab, steps), c in zip(stages, BLUE_RAMP):
        C = implanted_and_annealed(X, PWELL[0], PWELL[1], PWELL[2], steps)
        ax.semilogy(X, C, color=c, label=lab)
    ax.axhline(substrate_doping(), color=TEXT_2, ls="--", lw=1)
    ax.text(2.9, substrate_doping() * 1.3, "p-substrate (10 Ω·cm)", ha="right", color=TEXT_2, fontsize=8)
    ax.set(title="P-well (B 10$^{13}$, 180 keV) through the flow", xlabel="Depth (µm)",
           ylabel="Concentration (cm$^{-3}$)", xlim=(0, 3), ylim=(1e14, 1e18))
    ax.legend(loc="upper right", fontsize=8)
    ax.text(1.6, 3e16, "curves after drive-in overlap:\nthe drive-in dominates Dt", color=TEXT_2, fontsize=8)

    # (b) 단계별 Dt 기여: 같은 도펀트라도 어느 열처리가 분포를 정하는가
    ax = axes[1]
    steps = [("well drive-in\n1050 °C\n5 h", DRIVE_IN), ("gate ox\n800 °C\n37 min", [(GOX_T, gate_ox_time())]),
             ("furnace\n900 °C\n30 min", FURNACE), ("RTA\n1000 °C\n1 min", ANNEALS["RTA 1000 °C 1 min"]),
             ("silicide\n650 °C\n1 min", SILICIDE)]
    rows = []
    w = 0.38
    for k, (dop, c) in enumerate((("B", CAT[0]), ("As", CAT[1]))):
        vals = [sum(intrinsic_diffusivity(dop, T) * t * 3600 for T, t in st) for _, st in steps]
        ax.bar(np.arange(len(steps)) + (k - 0.5) * w, vals, w, color=c, label=dop)
        for (lab, st), v in zip(steps, vals):
            rows.append({"step": lab.replace("\n", " "), "dopant": dop, "Dt_cm2": f"{v:.3g}",
                         "diffusion_length_nm": round(2 * np.sqrt(v) * 1e7, 1)})
    ax.set_yscale("log")
    ax.set_xticks(range(len(steps)), [s[0] for s in steps], fontsize=8)
    ax.set(title="Thermal budget Dt of each step", ylabel="Dt (cm²)", ylim=(1e-20, 1e-8))
    ax.grid(axis="x", visible=False)
    ax.legend(loc="upper right", fontsize=8.5)
    fig.tight_layout()
    save(fig, "fig1_thermal_budget.png")
    write_csv("thermal_budget.csv", rows)


def fig2_vt_design():
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    doses = np.linspace(0, 6e12, 13)
    rows = []
    for ax, ch, dop, energies in ((axes[0], "n", "B", (10, 20, 50, 75)), (axes[1], "p", "As", (10, 30, 75, 100))):
        for E, c in zip(energies, BLUE_RAMP):
            v = []
            for d in doses:
                if ch == "n":
                    Badj = implanted_and_annealed(X, dop, E, d, [(GOX_T, gate_ox_time())] + FURNACE + SILICIDE) if d else 0 * X
                    base = channel_profiles(FURNACE)
                    v.append(vt_nonuniform(X, base["nmos"] + Badj, GOX_NM, "n", N_bulk=base["nmos_bulk"], gate="n+")[0])
                else:
                    Aadj = implanted_and_annealed(X, dop, E, d, [(GOX_T, gate_ox_time())] + FURNACE + SILICIDE) if d else 0 * X
                    base = channel_profiles(FURNACE)
                    v.append(vt_nonuniform(X, base["pmos"] + Aadj, GOX_NM, "p", N_bulk=base["pmos_bulk"], gate="p+")[0])
            Rp, _ = range_parameters(dop, E)
            ax.plot(doses / 1e12, v, color=c, marker="o", ms=3.5, label=f"{dop} {E} keV (R$_p$ {Rp * 1e3:.0f} nm)")
            for d, vv in zip(doses[::3], v[::3]):
                rows.append({"device": "NMOS" if ch == "n" else "PMOS", "dopant": dop, "energy_keV": E,
                             "dose_cm2": f"{d:.1e}", "Vt_V": round(vv, 3)})
        ax.axvspan(1, 5, color=GRID, lw=0)
        ax.axhline(VT_TARGET if ch == "n" else -VT_TARGET, color=TEXT_2, ls="--", lw=1)
        ax.text(5.9, (VT_TARGET if ch == "n" else -VT_TARGET) + (0.02 if ch == "n" else -0.05),
                f"target {'+' if ch == 'n' else '−'}{VT_TARGET:.2f} V", ha="right", color=TEXT_2, fontsize=8)
        ax.text(3, (0.008 if ch == "n" else -0.03), "textbook dose range", ha="center", color=TEXT_2, fontsize=8)
    axes[0].set(title=f"NMOS (n+ poly, t$_{{ox}}$ {GOX_NM:g} nm): B V$_t$ adjust", xlabel="Dose (10$^{12}$ cm$^{-2}$)",
                ylabel="V$_t$ (V)", ylim=(0, 0.5))
    axes[1].set(title=f"PMOS (p+ poly, t$_{{ox}}$ {GOX_NM:g} nm): As V$_t$ adjust", xlabel="Dose (10$^{12}$ cm$^{-2}$)",
                ylabel="V$_t$ (V)", ylim=(-1.0, 0))
    axes[0].legend(loc="upper left", fontsize=8)
    axes[1].legend(loc="lower left", fontsize=8)
    fig.tight_layout()
    save(fig, "fig2_vt_design.png")
    write_csv("vt_design_map.csv", rows)


def fig3_profiles_and_sensitivity(dn, dp):
    prof = channel_profiles(FURNACE, dn, dp)
    vn, vp, Wn, Wp = vts(prof)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    ax = axes[0]
    ax.semilogy(X, prof["nmos"], color=CAT[0], label=f"NMOS body: p-well + sub + B adjust (V$_t$ = {vn:+.2f} V)")
    ax.semilogy(X, np.clip(prof["pmos"], 1e13, None), color=CAT[1],
                label=f"PMOS body: n-well − sub + As adjust (V$_t$ = {vp:+.2f} V)")
    ax.axvline(Wn, color=CAT[0], ls=":", lw=1.2)
    ax.axvline(Wp, color=CAT[1], ls=":", lw=1.2)
    ax.text(Wn + 0.01, 3e18, f"W$_n$ = {Wn * 1e3:.0f} nm", color=CAT[0], fontsize=8)
    ax.text(Wp + 0.01, 1.2e18, f"W$_p$ = {Wp * 1e3:.0f} nm", color=CAT[1], fontsize=8)
    ax.set(title="Final channel doping (depletion width at threshold dotted)", xlabel="Depth (µm)",
           ylabel="Net doping (cm$^{-3}$)", xlim=(0, 0.6), ylim=(1e15, 1e19))
    ax.legend(loc="lower left", fontsize=7.5)

    # (b) 민감도
    ax = axes[1]
    base = (vn, vp)
    cases = []
    d_tox = (massoud_thickness(gate_ox_time(), GOX_T + 2) * 1e3 - GOX_NM)
    for lab, kw in (
        (f"gate-ox furnace +2 °C (t$_{{ox}}$ +{d_tox:.2f} nm)", dict(tox=GOX_NM + d_tox)),
        ("V$_t$ adjust dose +5 %", dict(dn=dn * 1.05, dp=dp * 1.05)),
        ("well drive-in +30 min", dict(drive=[(1050, 5.5)])),
        ("final anneal: RTA instead of furnace", dict(anneal=ANNEALS["RTA 1000 °C 1 min"])),
    ):
        p = channel_profiles(kw.get("anneal", FURNACE), kw.get("dn", dn), kw.get("dp", dp), drive=kw.get("drive", DRIVE_IN))
        v = vts(p, tox=kw.get("tox", GOX_NM))
        cases.append((lab, (v[0] - base[0]) * 1e3, (v[1] - base[1]) * 1e3))
    labels = [c[0] for c in cases]
    yy = np.arange(len(cases))
    ax.barh(yy - 0.18, [c[1] for c in cases], 0.34, color=CAT[0], label="ΔV$_{tn}$")
    ax.barh(yy + 0.18, [c[2] for c in cases], 0.34, color=CAT[1], label="ΔV$_{tp}$")
    for y_, c in zip(yy, cases):
        ax.text(c[1] + (1 if c[1] >= 0 else -1), y_ - 0.18, f"{c[1]:+.1f}", va="center",
                ha="left" if c[1] >= 0 else "right", fontsize=8, color=TEXT)
        ax.text(c[2] + (1 if c[2] >= 0 else -1), y_ + 0.18, f"{c[2]:+.1f}", va="center",
                ha="left" if c[2] >= 0 else "right", fontsize=8, color=TEXT)
    ax.axvline(0, color=TEXT_2, lw=1)
    ax.set_yticks(yy, labels, fontsize=8)
    ax.invert_yaxis()
    ax.set(title="Process variation → threshold shift", xlabel="ΔV$_t$ (mV)", xlim=(-35, 35))
    ax.grid(axis="y", visible=False)
    ax.legend(loc="lower right", fontsize=8.5)
    fig.tight_layout()
    save(fig, "fig3_channel_and_sensitivity.png")
    rows = [{"variation": c[0].replace("$", "").replace("_{ox}", "ox").replace("_t", "t"), "dVtn_mV": round(c[1], 1),
             "dVtp_mV": round(c[2], 1)} for c in cases]
    write_csv("vt_sensitivity.csv", rows)
    return vn, vp, Wn, Wp, cases


def ldd_case(anneal_steps, x, y):
    Dt = sum(intrinsic_diffusivity(LDD_N[0], T) * t * 3600 for T, t in anneal_steps + SILICIDE)
    C = ldd_2d_analytic(x, y, LDD_N[0], LDD_N[1], LDD_N[2], Dt)
    return C, Dt


def fig4_ldd_leff(dn, dp):
    x = np.linspace(-0.3, 0.3, 601)
    y = np.linspace(0, 0.4, 401)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.3), sharey=True)
    rows = []
    for ax, (name, steps) in zip(axes, ANNEALS.items()):
        Nbody = np.interp(y, X, channel_profiles(steps, dn, dp)["nmos"])   # 같은 열처리를 거친 채널 도핑
        C, Dt = ldd_case(steps, x, y)
        net = np.log10(np.clip(C, 1e14, None))
        cs = ax.contourf(x, y, net, levels=np.arange(15, 19.6, 0.5), cmap="Blues")
        ax.contour(x, y, C - Nbody[:, None], levels=[0], colors=[CAT[1]], linewidths=2)
        reach = lateral_junction_extent(x, y, C, Nbody)
        xj = junction_depth_between(y, C[:, 0], Nbody)
        ax.add_patch(Rectangle((0, -0.06), 0.3, 0.06, color="#c9c7c0", clip_on=False))
        ax.text(0.15, -0.035, "n+ poly gate", ha="center", va="center", fontsize=8, color=TEXT)
        ax.axvline(0, color=TEXT_2, lw=1, ls=":")
        ax.set(title=f"{name}\nlateral reach {reach * 1e3:.0f} nm, x$_j$ {xj * 1e3:.0f} nm",
               xlabel="x (µm, gate edge at 0)", xlim=(-0.3, 0.3), ylim=(0.25, -0.06))
        leff = GATE_L_UM - 2 * reach
        rows.append({"anneal": name, "Dt_P_cm2": f"{Dt:.3g}", "lateral_reach_nm": round(reach * 1e3, 1),
                     "LDD_xj_nm": round(xj * 1e3, 1), "Lgate_um": GATE_L_UM, "Leff_um": round(leff, 3)})
    axes[0].set_ylabel("Depth (µm)")
    fig.colorbar(cs, ax=axes, label="log$_{10}$ LDD P (cm$^{-3}$)", shrink=0.85)
    fig.suptitle("NMOS LDD (P 5×10$^{13}$, 50 keV) under the gate edge (orange: junction)", color=TEXT,
                 fontweight="bold", y=1.03)
    save(fig, "fig4_ldd_leff.png")
    write_csv("ldd_leff.csv", rows)
    return rows


def spacer_width():
    """고립된 poly 선의 측벽에 LPCVD TEOS 산화막이 쌓이는 두께 (넓은 트렌치로 근사)."""
    el = trench_elements(20 * POLY_UM, POLY_UM, n_per_side=40, n_bottom=200)
    F, _, _ = solve_trench_flux(el, 20 * POLY_UM, n=1, Sc=0.075)
    m = step_coverage_metrics(el, F)
    return m["sidewall_mid"] * SPACER_DEP_UM, m["sidewall_mid"]


def fig5_cross_section(vn, vp, leff, sw):
    """계산한 치수로 그린 CMOS 단면 (깊이 방향 치수는 계산값, 게이트 산화막은 보이도록 과장)."""
    NA = substrate_doping()
    t_field = deal_grove_thickness(LOCOS[1], LOCOS[0], "wet", x_i_um=PAD_OX_UM)
    prof = channel_profiles(FURNACE)
    xj_nwell = junction_depth_between(X, prof["Pw"], prof["NA"] + 0 * X)
    # p-well은 기판과 같은 p형이므로 접합이 없다. 웰 B 농도가 기판 농도와 같아지는 깊이를 웰 깊이로 둔다.
    xj_pwell = junction_depth_between(X, prof["Bw"], prof["NA"] + 0 * X)
    fig, ax = plt.subplots(figsize=(12, 5.2))
    W = 8.0
    ax.add_patch(Rectangle((0, -3.0), W, 3.0, color="#efe9df", lw=0))
    ax.add_patch(Rectangle((0, -xj_pwell), W / 2, xj_pwell, color="#f6d3c2", lw=0))
    ax.add_patch(Rectangle((W / 2, -xj_nwell), W / 2, xj_nwell, color="#d6e6fb", lw=0))
    ax.text(0.15, -2.8, "p-substrate (10 Ω·cm)", fontsize=8, color=TEXT_2)
    ax.text(0.15, -xj_pwell + 0.1, f"P-well (depth {xj_pwell:.2f} µm)", fontsize=8, color=TEXT_2)
    ax.text(W / 2 + 0.15, -xj_nwell + 0.1, f"N-well (x$_j$ = {xj_nwell:.2f} µm)", fontsize=8, color=TEXT_2)
    # LOCOS 필드 산화막: Si 소모 비율만큼 아래로 파고든다
    for x0 in (0.0, W / 2 - 0.6, W - 0.6):
        wdt = 1.2 if 0 < x0 < W - 1 else 0.6
        ax.add_patch(Rectangle((x0, -SI_CONSUMPTION_RATIO * t_field), wdt, t_field, color="#bcd9c8", lw=0))
    ax.text(W / 2, t_field * 0.55, f"LOCOS {t_field:.2f} µm", ha="center", fontsize=8, color=TEXT)
    active = {"NMOS": (0.6, W / 2 - 0.6), "PMOS": (W / 2 + 0.6, W - 0.6)}
    for xc, ch, sd_col in ((2.0, "NMOS", "#2a78d6"), (6.0, "PMOS", "#eb6834")):
        Lg = GATE_L_UM
        ax.add_patch(Rectangle((xc - Lg / 2, 0.0), Lg, 0.06, color="#7d7a73", lw=0))           # gate oxide (exaggerated)
        ax.add_patch(Rectangle((xc - Lg / 2, 0.06), Lg, POLY_UM, color="#9b968c", lw=0))     # poly
        for s in (-1, 1):
            xs = xc + s * Lg / 2
            ax.add_patch(Polygon([[xs, 0.06], [xs, 0.06 + POLY_UM], [xs + s * sw, 0.06]], color="#bcd9c8", lw=0))
            xin = xs + s * sw
            edge = active[ch][0] if s < 0 else active[ch][1]
            ax.add_patch(Rectangle((min(xin, edge), -0.18), abs(edge - xin), 0.18, color=sd_col, alpha=0.75, lw=0))
            ax.add_patch(Rectangle((min(xs, xin), -0.07), abs(sw), 0.07, color=sd_col, alpha=0.4, lw=0))
        ax.text(xc, 0.06 + POLY_UM + 0.08, f"{ch}\nV$_t$ = {vn if ch == 'NMOS' else vp:+.2f} V", ha="center",
                fontsize=8.5, color=TEXT)
    ax.annotate(f"spacer {sw * 1e3:.0f} nm", (2.0 + GATE_L_UM / 2 + sw / 2, 0.15), xytext=(2.9, 0.9),
                fontsize=8, color=TEXT_2, arrowprops=dict(arrowstyle="-", color=TEXT_2, lw=0.8))
    ax.annotate(f"L$_g$ {GATE_L_UM:g} µm → L$_{{eff}}$ {leff:.2f} µm", (2.0, -0.02), xytext=(0.4, 1.0),
                fontsize=8, color=TEXT_2, arrowprops=dict(arrowstyle="-", color=TEXT_2, lw=0.8))
    ax.set(xlim=(0, W), ylim=(-3.0, 1.4), xlabel="Lateral position (µm, schematic)", ylabel="Depth (µm)")
    ax.set_title("CMOS cross-section assembled from module results (depths to scale, gate oxide exaggerated)",
                 fontsize=10)
    ax.grid(False)
    fig.tight_layout()
    save(fig, "fig5_cmos_cross_section.png")
    return t_field, xj_nwell


def fig6_validation():
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    rows = []
    Xv = np.linspace(0, 3, 30001)
    ax = axes[0]
    Ns = np.logspace(15, 18, 13)
    for ch, gate, c, lab in (("n", "n+", CAT[0], "NMOS, n+ poly"), ("p", "p+", CAT[1], "PMOS, p+ poly")):
        an = [vt_uniform(N, GOX_NM, ch, gate=gate) for N in Ns]
        nu = [vt_nonuniform(Xv, np.full_like(Xv, N), GOX_NM, ch, gate=gate)[0] for N in Ns]
        ax.semilogx(Ns, an, color=c, label=f"{lab}: analytic")
        ax.semilogx(Ns, nu, ls="none", marker="o", ms=5, mfc=SURFACE, mec=c, mew=1.3, label=f"{lab}: depletion solver")
        rows.append({"case": f"uniform doping {lab}", "metric": "max |dVt| V", "value": f"{np.max(np.abs(np.array(an) - nu)):.1e}"})
    ax.set(title="Uniform doping: solver vs. analytic V$_t$", xlabel="Body doping (cm$^{-3}$)", ylabel="V$_t$ (V)")
    ax.legend(loc="center left", fontsize=7.5)

    ax = axes[1]
    N, Q = 1e16, 1e10
    v0, W = vt_nonuniform(Xv, np.full_like(Xv, N), GOX_NM, "n", N_bulk=N)
    depths = np.linspace(0.005, W * 0.95, 15)
    sim = []
    for dpt in depths:
        prof = np.full_like(Xv, N)
        prof[int(round(dpt / (Xv[1] - Xv[0])))] += Q / ((Xv[1] - Xv[0]) * 1e-4)
        sim.append((vt_nonuniform(Xv, prof, GOX_NM, "n", N_bulk=N)[0] - v0) * 1e3)
    qq = Q_E * Q / oxide_capacitance(GOX_NM) * 1e3
    ax.plot(depths * 1e3, qq * (1 - depths / W), color=CAT[0], label="(qQ/C$_{ox}$)(1 − x$_i$/W)")
    ax.plot(depths * 1e3, sim, ls="none", marker="o", ms=5, mfc=SURFACE, mec=CAT[1], mew=1.3, label="depletion solver")
    ax.set(title=f"Thin implanted layer at depth x$_i$ (W = {W * 1e3:.0f} nm)", xlabel="Depth of layer x$_i$ (nm)",
           ylabel="ΔV$_t$ (mV)")
    ax.legend(loc="upper right", fontsize=8.5)
    err = np.max(np.abs(np.array(sim) - qq * (1 - depths / W)) / qq) * 100
    rows.append({"case": "buried sheet efficiency (1 - x/W)", "metric": "max error % of qQ/Cox", "value": round(err, 2)})

    ax = axes[2]
    x = np.linspace(-0.4, 0.4, 321)
    y = np.linspace(0, 0.5, 201)
    Dt = intrinsic_diffusivity("P", 900) * 1800
    C0 = ldd_2d_analytic(x, y, "P", 50, 5e13, 0.0)
    C_fd = diffuse_2d_fd(x, y, C0, Dt, n_steps=200)
    C_an = ldd_2d_analytic(x, y, "P", 50, 5e13, Dt)
    lv = [17, 18]
    ax.contour(x, y, np.log10(np.clip(C_an, 1e10, None)), levels=lv, colors=[CAT[0]], linewidths=2.2)
    ax.contour(x, y, np.log10(np.clip(C_fd, 1e10, None)), levels=lv, colors=[CAT[1]], linewidths=1.3,
               linestyles="--")
    ax.plot([], [], color=CAT[0], lw=2.2, label="analytic (separable)")
    ax.plot([], [], color=CAT[1], lw=1.3, ls="--", label="2D numerical")
    rel = np.max(np.abs(C_fd - C_an)) / C_an.max() * 100
    ax.set(title=f"2D LDD diffusion: max diff {rel:.2f} % of peak", xlabel="x (µm)", ylabel="Depth (µm)",
           ylim=(0.3, 0))
    ax.legend(loc="lower right", fontsize=8.5)
    rows.append({"case": "2D LDD analytic vs numerical (900 C 30 min)", "metric": "max diff % of peak", "value": round(rel, 3)})
    fig.tight_layout()
    save(fig, "fig6_validation.png")
    write_csv("validation.csv", rows)


def recipe_sheet(dn, dp, vn, vp, ldd_rows, sw, sw_frac, t_field, xj_nwell):
    """Plummer Ch. 2 공정 흐름의 각 단계를 이 프로젝트 모듈로 계산한 결과."""
    NA = substrate_doping()
    t_g = gate_ox_time()
    S_poly = required_selectivity(POLY_UM, GOX_NM / 2 * 1e-3, 0.05, 0.05, 0.2)
    Rp_n, dRp_n = range_parameters(*SD_N[:2])
    sheet = [
        ("Substrate", "p-type (100), 5~50 Ω·cm", f"10 Ω·cm → N_A = {NA:.2e} cm^-3", "02 (Masetti)"),
        ("LOCOS", "90 min @ 1000 °C H2O → ≈ 0.5 µm", f"{t_field:.2f} µm (Si 소모 {SI_CONSUMPTION_RATIO * t_field:.2f} µm)", "01"),
        ("Wells", "B, P 1e13; drive-in 4~6 h @ 1000~1100 °C", f"5 h @ 1050 °C → N-well x_j {xj_nwell:.2f} µm", "02"),
        ("Gate oxide", "3~5 nm, 0.5~1 h @ 800 °C O2", f"4 nm → {t_g * 60:.0f} min @ 800 °C (Massoud)", "01"),
        ("Vt adjust", "B / As 1~5e12 @ 50~100 keV", f"B {VT_ADJ_N[1]} keV {dn:.2e}, As {VT_ADJ_P[1]} keV {dp:.2e} → Vt {vn:+.2f} / {vp:+.2f} V", "02 + 06"),
        ("Gate litho", "-", f"L_g {GATE_L_UM} µm, i-line NA 0.5: k1 {k1_factor(GATE_L_UM, 0.365, 0.5):.2f}, DOF ±{depth_of_focus(0.365, 0.5):.2f} µm", "03"),
        ("Poly etch", "anisotropic plasma etch, poly 0.5 µm", f"gate oxide loss ≤ 2 nm needs selectivity ≥ {S_poly:.0f} (OE 20 %)", "04"),
        ("Spacer", "conformal SiO2 0.5 µm + anisotropic etch", f"TEOS sidewall coverage {sw_frac:.2f} → spacer {sw * 1e3:.0f} nm", "05 + 04"),
        ("LDD + final anneal", "P/B 5e13 @ 50 keV; 30 min @ 900 °C or 1 min RTA @ 1000 °C",
         "; ".join(f"{r['anneal']}: L_eff {r['Leff_um']:.2f} µm" for r in ldd_rows), "02 + 06"),
        ("S/D", "As 2~4e15 @ 75 keV, B 1~3e15 @ 50 keV",
         f"As 3e15: R_p {Rp_n * 1e3:.0f} nm, peak {peak_concentration(SD_N[2], dRp_n):.1e} cm^-3 (진성 확산 범위 밖)", "02"),
        ("Contact W fill", "TiN barrier + W CVD", "W (S_c 0.01) min coverage 0.95 at AR 2", "05"),
        ("Metal 1", "Al sputter + plasma etch", f"Al 0.5 µm: R_s {sheet_from_resistivity(3.0, 500) * 1e3:.0f} mΩ/sq; 10 mm line {rc_delay(3.0, 10, 1, 1, 1, 1, 3.9) * 1e12:.0f} ps", "05"),
    ]
    rows = [{"step": a, "textbook_condition": b, "this_project": c, "module": d} for a, b, c, d in sheet]
    write_csv("recipe_sheet.csv", rows)
    return rows


def main():
    print("Process integration module")
    fig1_thermal_budget()
    dn, dp = design_vt_doses()
    print(f"  Vt adjust doses for |Vt| = {VT_TARGET} V: B {VT_ADJ_N[1]} keV {dn:.2e}, As {VT_ADJ_P[1]} keV {dp:.2e}")
    fig2_vt_design()
    vn, vp, Wn, Wp, cases = fig3_profiles_and_sensitivity(dn, dp)
    print(f"  Vtn {vn:+.3f} V (W {Wn * 1e3:.0f} nm), Vtp {vp:+.3f} V (W {Wp * 1e3:.0f} nm)")
    ldd_rows = fig4_ldd_leff(dn, dp)
    sw, sw_frac = spacer_width()
    t_field, xj_nwell = fig5_cross_section(vn, vp, ldd_rows[0]["Leff_um"], sw)
    fig6_validation()
    for r in recipe_sheet(dn, dp, vn, vp, ldd_rows, sw, sw_frac, t_field, xj_nwell):
        print(f"  {r['step']:<18} {r['this_project']}")


if __name__ == "__main__":
    main()
