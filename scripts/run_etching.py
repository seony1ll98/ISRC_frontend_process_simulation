"""건식 식각 모듈의 그래프와 결과 표를 모두 다시 만든다.

실행: python scripts/run_etching.py
출력: figures/etching/*.png, results/etching/*.csv
"""

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from plot_style import BLUE_RAMP, CAT, GRID, SURFACE, TEXT, TEXT_2  # noqa: E402
from plot_style import save as _save  # noqa: E402

from isrc_sim.etching import (  # noqa: E402
    MASK,
    VOID,
    anisotropy_from_bias,
    directional_flux_factor,
    make_structure,
    profile_metrics,
    rate_linear,
    rate_saturation,
    required_mask_thickness,
    required_selectivity,
    simulate_profile,
    underlayer_loss,
)

FIG_DIR = ROOT / "figures" / "etching"
RES_DIR = ROOT / "results" / "etching"

OPENING, FILM_T, MASK_T, UNDER_T = 0.5, 0.5, 0.3, 0.25
DX = 0.0025
# 재료 색: 빈 공간, 마스크(PR), 식각할 막, 하부막
MAT_CMAP = ListedColormap([SURFACE, "#f6c9b0", "#cde2fb", "#d9d7d0"])

# 과식각 시나리오: 01 모듈에서 설계한 10 nm 게이트 산화막 위의 poly-Si 게이트 식각
POLY_T, THK_NU, RATE_NU = 0.30, 0.05, 0.05


def save(fig, name):
    _save(fig, FIG_DIR, name)


def write_csv(name, rows):
    RES_DIR.mkdir(parents=True, exist_ok=True)
    with open(RES_DIR / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"  saved results/etching/{name}")


def fig1_rate_models():
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    ax = axes[0]
    ion = np.linspace(0, 4, 200)
    for chem, c in zip((1.0, 0.5, 0.2, 0.0), BLUE_RAMP[::-1]):
        ax.plot(ion, rate_saturation(chem, ion), color=c, label=f"S$_c$F$_c$ = {chem:g}")
    ax.set(title="Saturation-adsorption model (ion-enhanced)", xlabel="Ion flux × K$_i$",
           ylabel="Etch rate × N", xlim=(0, 4), ylim=(0, 1.05))
    ax.legend(loc="upper left", fontsize=8.5)

    ax = axes[1]
    ratio = np.logspace(-1.3, 1.3, 200)
    Rc = 1.0
    rows = []
    for label, n, model, c, ls in (("linear, n = 40", 40, "linear", CAT[1], "--"),
                                   ("saturation, n = 10", 10, "saturation", BLUE_RAMP[1], "-"),
                                   ("saturation, n = 40", 40, "saturation", BLUE_RAMP[2], "-"),
                                   ("saturation, n = 80", 80, "saturation", BLUE_RAMP[3], "-")):
        G = directional_flux_factor(np.pi / 2, n)
        f = rate_saturation if model == "saturation" else rate_linear
        af = 1 - f(Rc, ratio * Rc * G) / f(Rc, ratio * Rc)
        ax.semilogx(ratio, af, color=c, ls=ls, label=label)
        for r in (0.25, 1, 4):
            rows.append({"model": model, "ion_angular_n": n, "sidewall_ion_fraction": round(G, 4), "R_ion_over_R_chem": r,
                         "anisotropy_Af": round(float(1 - f(Rc, r * Rc * G) / f(Rc, r * Rc)), 3)})
    ax.set(title="Initial anisotropy A$_f$ = 1 − r$_{lat}$/r$_{ver}$", xlabel="Ion / chemical component (R$_i$ / R$_c$)",
           ylabel="A$_f$", ylim=(0, 1.02))
    ax.legend(loc="lower right", fontsize=8.5)
    fig.tight_layout()
    save(fig, "fig1_rate_models.png")
    write_csv("anisotropy_design.csv", rows)


def run_case(Rc, Ri, model, ion_n=None, sel_under=np.inf, t_extra=0.0, n_snapshots=5):
    x, y, mat = make_structure(OPENING, FILM_T, mask_um=MASK_T, under_um=UNDER_T, width_um=1.8, dx_um=DX)
    if ion_n is None:
        v_floor = (rate_saturation if model == "saturation" else rate_linear)(Rc, Ri)
    else:
        v_floor = (rate_saturation if model == "saturation" else rate_linear)(Rc, Ri)
    t_clear = FILM_T / float(v_floor)
    times, snaps = simulate_profile(x, y, mat, Rc, Ri, t_clear * (1 + t_extra), model=model, ion_n=ion_n,
                                    sel_under=sel_under, n_snapshots=n_snapshots)
    return x, y, mat, times, snaps, t_clear


def draw_profile(ax, x, y, mat, snaps, title):
    final = np.where(snaps[-1] & (mat != MASK), VOID, mat)
    ax.imshow(final, origin="lower", extent=(x[0], x[-1], y[0], y[-1]), cmap=MAT_CMAP, vmin=0, vmax=3,
              interpolation="nearest", aspect="equal")
    for s, c in zip(snaps[1:-1], BLUE_RAMP):
        ax.contour(x, y, s.astype(float), levels=[0.5], colors=[c], linewidths=1.0)
    ax.contour(x, y, snaps[-1].astype(float), levels=[0.5], colors=[TEXT], linewidths=1.4)
    ax.set_title(title, fontsize=10)
    ax.set_xlim(-0.75, 0.75)
    ax.set_ylim(-FILM_T - UNDER_T, MASK_T + 0.05)
    ax.set_xlabel("x (µm)")
    ax.grid(False)


def fig2_profiles():
    cases = [
        ("(a) chemical only", dict(Rc=0.1, Ri=0.0, model="linear")),
        ("(b) ions only (vertical)", dict(Rc=0.0, Ri=0.1, model="linear")),
        ("(c) linear, half and half", dict(Rc=0.05, Ri=0.05, model="linear")),
        ("(d) ion-enhanced (saturation)", dict(Rc=0.1, Ri=0.1, model="saturation")),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(15, 4.4), sharey=True)
    rows = []
    for ax, (title, kw) in zip(axes, cases):
        x, y, mat, times, snaps, t_clear = run_case(kw["Rc"], kw["Ri"], kw["model"], kw.get("ion_n"))
        m = profile_metrics(x, y, snaps[-1], at_depth_um=0.1)
        bias = max(m["half_width_um"] - OPENING / 2, 0.0)
        af = anisotropy_from_bias(bias, m["depth_um"])
        draw_profile(ax, x, y, mat, snaps, f"{title}\nbias {bias * 1e3:.0f} nm, A$_f$ = {af:.2f}, {t_clear:.0f} min")
        rows.append({"case": title, "R_chem": kw["Rc"], "R_ion": kw["Ri"], "model": kw["model"],
                     "time_to_clear_min": round(t_clear, 2), "depth_um": round(m["depth_um"], 3),
                     "bias_nm_at_0.1um": round(bias * 1e3), "Af": round(af, 2)})
    axes[0].set_ylabel("y (µm)")
    axes[1].text(-0.72, MASK_T / 2, "PR mask", fontsize=8, color=TEXT_2, va="center")
    axes[1].text(-0.72, -0.08, "film", fontsize=8, color=TEXT_2)
    axes[1].text(-0.72, -FILM_T - 0.1, "underlayer", fontsize=8, color=TEXT_2)
    fig.suptitle("Etch profiles until the film clears at the center, vertical ions (thin lines: earlier times)", color=TEXT,
                 fontweight="bold")
    fig.tight_layout()
    save(fig, "fig2_profiles.png")
    write_csv("profiles.csv", rows)


def fig3_ion_angle():
    fig, ax = plt.subplots(figsize=(6.6, 4.3))
    ns = [10, 20, 40, 80]
    a = np.linspace(0, np.pi / 2, 91)
    rows = []
    for n, c in zip(ns, BLUE_RAMP):
        g = directional_flux_factor(a, n)
        ax.plot(np.degrees(a), g, color=c, label=f"n = {n}: sidewall {g[-1]:.3f}")
        rows.append({"ion_n": n, "flux_at_45deg": round(float(np.interp(np.pi / 4, a, g)), 3),
                     "flux_at_80deg": round(float(np.interp(np.radians(80), a, g)), 3),
                     "sidewall_flux_90deg": round(float(g[-1]), 4)})
    ax.plot(np.degrees(a), np.cos(a), color=TEXT_2, ls=":", lw=1.2, label="cos α (perfectly vertical)")
    ax.set(title="Ion flux on a tilted surface, cos$^n$θ arrival (flat = 1)", xlabel="Surface tilt α (deg)",
           ylabel="Relative ion flux", xlim=(0, 90), ylim=(0, 1.02))
    ax.legend(loc="lower left", fontsize=8.5)
    fig.tight_layout()
    save(fig, "fig3_ion_angle.png")
    write_csv("ion_angular_flux.csv", rows)


def fig4_overetch():
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    S = np.logspace(0.7, 2.3, 200)
    rows = []
    ax = axes[0]
    for oe, c in zip((0.1, 0.2, 0.5), BLUE_RAMP[1:]):
        loss = [underlayer_loss(POLY_T, 1.0, s, THK_NU, RATE_NU, oe) * 1e3 for s in S]
        ax.loglog(S, loss, color=c, label=f"over-etch {oe * 100:.0f} %")
        rows.append({"overetch_pct": oe * 100,
                     "S_for_loss_le_5nm": round(required_selectivity(POLY_T, 0.005, THK_NU, RATE_NU, oe), 1),
                     "S_for_loss_le_2nm": round(required_selectivity(POLY_T, 0.002, THK_NU, RATE_NU, oe), 1),
                     "loss_at_S25_nm": round(underlayer_loss(POLY_T, 1.0, 25, THK_NU, RATE_NU, oe) * 1e3, 2),
                     "loss_at_S50_nm": round(underlayer_loss(POLY_T, 1.0, 50, THK_NU, RATE_NU, oe) * 1e3, 2)})
    ax.axhline(10, color=TEXT_2, lw=1.2)
    ax.text(5.2, 11, "gate oxide 10 nm (module 01)", color=TEXT_2, fontsize=8)
    ax.axhline(5, color=TEXT_2, ls="--", lw=1)
    ax.text(5.2, 5.5, "half of the oxide", color=TEXT_2, fontsize=8)
    ax.axvspan(25, 50, color=GRID, lw=0)
    ax.text(26, 0.6, "textbook target\nS = 25~50", color=TEXT_2, fontsize=8)
    ax.set(title="Gate oxide loss in poly-Si gate etch (300 nm, ±5 %)", xlabel="Selectivity poly : oxide",
           ylabel="Worst-case oxide loss (nm)", xlim=(S[0], S[-1]), ylim=(0.5, 100))
    ax.legend(loc="upper right", fontsize=8.5)

    ax = axes[1]
    Sm = np.linspace(1.5, 10, 100)
    for oe, c in zip((0.1, 0.2, 0.5), BLUE_RAMP[1:]):
        ax.plot(Sm, [required_mask_thickness(POLY_T, s, THK_NU, RATE_NU, oe) * 1e3 for s in Sm], color=c,
                label=f"over-etch {oe * 100:.0f} %")
    ax.set(title="Minimum PR thickness (fully consumed at the end)", xlabel="Selectivity poly : PR",
           ylabel="PR thickness (nm)", ylim=(0, 300))
    ax.legend(loc="upper right", fontsize=8.5)
    fig.tight_layout()
    save(fig, "fig4_overetch.png")
    write_csv("overetch_selectivity.csv", rows)


def validation():
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))
    rows = []
    # (a) 등방성 식각: level set vs 정확해(개구에서 거리 r 이내)
    ax = axes[0]
    x, y, mat = make_structure(OPENING, FILM_T, mask_um=MASK_T, under_um=UNDER_T, width_um=1.8, dx_um=DX)
    r = 0.4
    _, snaps = simulate_profile(x, y, mat, 0.1, 0.0, r / 0.1, n_snapshots=2)
    X, Y = np.meshgrid(x, y)
    exact = (Y < 0) & (np.hypot(np.clip(np.abs(X) - OPENING / 2, 0, None), Y) <= r)
    sim = snaps[-1] & (Y < 0)
    shift = np.sum(sim ^ exact) * DX**2 / (OPENING + np.pi * r)
    ax.contour(x, y, exact.astype(float), levels=[0.5], colors=[CAT[0]], linewidths=2.2)
    ax.contour(x, y, sim.astype(float), levels=[0.5], colors=[CAT[1]], linewidths=1.4, linestyles="--")
    ax.plot([], [], color=CAT[0], lw=2.2, label="exact (distance ≤ R·t)")
    ax.plot([], [], color=CAT[1], lw=1.4, ls="--", label="level set")
    ax.fill_between([OPENING / 2, 0.9], 0, MASK_T, color="#f6c9b0", lw=0)
    ax.fill_between([-0.9, -OPENING / 2], 0, MASK_T, color="#f6c9b0", lw=0)
    ax.set(title=f"Isotropic etch: mean boundary error {shift * 1e3:.1f} nm", xlabel="x (µm)", ylabel="y (µm)",
           xlim=(-0.8, 0.8), ylim=(-0.5, 0.35), aspect="equal")
    ax.legend(loc="upper center", fontsize=8)
    rows.append({"case": f"isotropic etch vs exact (r = 0.4 um, dx = {DX * 1e3:g} nm)", "metric": "mean boundary error nm",
                 "value": round(shift * 1e3, 2)})

    # (b) 격자 수렴: 등방성 경계 오차, 선형 모델의 깊이와 bias 오차
    ax = axes[1]
    dxs = [0.02, 0.01, 0.005, 0.0025]
    e_iso, e_depth, e_bias = [], [], []
    for d in dxs:
        x, y, mat = make_structure(OPENING, FILM_T, mask_um=MASK_T, under_um=UNDER_T, width_um=1.8, dx_um=d)
        X, Y = np.meshgrid(x, y)
        _, sn = simulate_profile(x, y, mat, 0.1, 0.0, 4.0, n_snapshots=2)
        ex = (Y < 0) & (np.hypot(np.clip(np.abs(X) - OPENING / 2, 0, None), Y) <= 0.4)
        e_iso.append(np.sum((sn[-1] & (Y < 0)) ^ ex) * d**2 / (OPENING + np.pi * 0.4) * 1e3)
        _, sn = simulate_profile(x, y, mat, 0.05, 0.05, 4.0, n_snapshots=2)
        m = profile_metrics(x, y, sn[-1], at_depth_um=0.1)
        e_depth.append(abs(m["depth_um"] - 0.4) * 1e3)
        e_bias.append(abs(m["half_width_um"] - OPENING / 2 - 0.2) * 1e3)
    g = np.array(dxs) * 1e3
    ax.loglog(g, e_iso, marker="o", color=CAT[0], label="isotropic: mean boundary error")
    ax.loglog(g, e_bias, marker="s", color=CAT[1], label="linear model: bias error")
    ax.loglog(g, g, color=TEXT_2, ls=":", lw=1.2, label="= grid size (1st order)")
    ax.set(title="Grid convergence against exact solutions", xlabel="Grid size (nm)", ylabel="Absolute error (nm)")
    ax.set_xticks(g, [f"{v:g}" for v in g])
    ax.minorticks_off()
    ax.legend(loc="upper left", fontsize=8)
    for d, a_, b_, c_ in zip(dxs, e_iso, e_depth, e_bias):
        rows.append({"case": f"dx = {d * 1e3:g} nm", "metric": "isotropic boundary / linear depth / linear bias error nm",
                     "value": f"{a_:.1f} / {b_:.2f} / {c_:.1f}"})

    # (c) 이온 플럭스 계수: n = 1 해석해, 큰 n은 cos
    ax = axes[2]
    a = np.linspace(0, np.pi, 61)
    ax.plot(np.degrees(a), (1 + np.cos(a)) / 2, color=CAT[0], label="analytic n = 1: (1 + cos α)/2")
    g1 = directional_flux_factor(a, 1)
    ax.plot(np.degrees(a[::3]), g1[::3], ls="none", marker="o", ms=5, mfc=SURFACE, mec=CAT[1], mew=1.3,
            label="numerical n = 1")
    ax.plot(np.degrees(a), np.clip(np.cos(a), 0, None), color=TEXT_2, ls=":", label="limit n → ∞: cos α")
    g400 = directional_flux_factor(a, 400)
    ax.plot(np.degrees(a[::3]), g400[::3], ls="none", marker="s", ms=4.5, mfc=SURFACE, mec=CAT[2], mew=1.3,
            label="numerical n = 400")
    ax.set(title=f"Flux factor: max error n = 1 is {np.max(np.abs(g1 - (1 + np.cos(a)) / 2)):.1e}",
           xlabel="Surface tilt α (deg)", ylabel="Relative flux", xlim=(0, 180))
    ax.legend(loc="upper right", fontsize=8)
    rows.append({"case": "flux factor n=1 vs (1+cos)/2", "metric": "max abs error",
                 "value": f"{np.max(np.abs(g1 - (1 + np.cos(a)) / 2)):.1e}"})
    fig.tight_layout()
    save(fig, "fig5_validation.png")
    write_csv("validation.csv", rows)


def main():
    print("Etching module")
    fig1_rate_models()
    fig2_profiles()
    fig3_ion_angle()
    fig4_overetch()
    validation()


if __name__ == "__main__":
    main()
