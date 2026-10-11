"""포토리소그래피 모듈의 그래프와 결과 표를 모두 다시 만든다.

실행: python scripts/run_lithography.py
출력: figures/lithography/*.png, results/lithography/*.csv
"""

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from plot_style import BLUE_RAMP, CAT, GRID, SURFACE, TEXT, TEXT_2  # noqa: E402
from plot_style import save as _save  # noqa: E402

from isrc_sim.lithography import (  # noqa: E402
    aerial_image,
    angular_spectrum_propagate,
    coherent_defocus_contrast_analytic,
    cutoff_pitch,
    dense_ls_contrast_analytic,
    depth_of_focus,
    dose_to_size,
    image_contrast,
    k1_factor,
    proximity_min_feature,
    rayleigh_resolution,
    resist_cmtf,
    resist_line_cd,
    slit_field,
    slit_fresnel_intensity,
    threshold_width,
)

FIG_DIR = ROOT / "figures" / "lithography"
RES_DIR = ROOT / "results" / "lithography"

I_LINE = 0.365
# Plummer Ch. 5에 나오는 노광 장비 예시 (g-line, i-line은 ATHENA 시뮬레이션 예, KrF는 Rayleigh 예제)
TOOLS = [
    {"name": "g-line", "wavelength": 0.436, "NA": 0.43},
    {"name": "i-line", "wavelength": 0.365, "NA": 0.50},
    {"name": "KrF", "wavelength": 0.248, "NA": 0.60},
]
TOOL_COLOR = {"g-line": CAT[0], "i-line": CAT[1], "KrF": CAT[2]}
SIGMA = 0.6                      # Plummer: 부분 결맞음 S = 0.5~0.7
CMTF_ILINE = resist_cmtf(2.5)    # g/i-line 레지스트 gamma 2~3
CMTF_DUV = resist_cmtf(7.0)      # DUV 레지스트 gamma 5~10


def save(fig, name):
    _save(fig, FIG_DIR, name)


def write_csv(name, rows):
    RES_DIR.mkdir(parents=True, exist_ok=True)
    with open(RES_DIR / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"  saved results/lithography/{name}")


def fig1_rayleigh():
    NA = np.linspace(0.3, 0.95, 200)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for name, lam, c in (("g-line 436 nm", 0.436, CAT[0]), ("i-line 365 nm", 0.365, CAT[1]),
                         ("KrF 248 nm", 0.248, CAT[2]), ("ArF 193 nm", 0.193, CAT[3])):
        axes[0].plot(NA, rayleigh_resolution(lam, NA) * 1e3, color=c, label=name)
        axes[1].plot(NA, depth_of_focus(lam, NA) * 1e3, color=c, label=name)
    rows = []
    for t in TOOLS:
        R = rayleigh_resolution(t["wavelength"], t["NA"])
        D = depth_of_focus(t["wavelength"], t["NA"])
        for ax, y in ((axes[0], R), (axes[1], D)):
            ax.plot(t["NA"], y * 1e3, marker="o", ms=7, color=TOOL_COLOR[t["name"]], mec=SURFACE, mew=1.5)
        axes[0].annotate(f"{t['name']} NA {t['NA']}", (t["NA"], R * 1e3), xytext=(6, 4),
                         textcoords="offset points", fontsize=8, color=TEXT_2)
        rows.append({"tool": t["name"], "wavelength_nm": round(t["wavelength"] * 1e3), "NA": t["NA"],
                     "R_k1_0.61_nm": round(R * 1e3), "DOF_k2_0.5_um": round(D, 3)})
    axes[0].set(title="Resolution R = 0.61 λ / NA", xlabel="Numerical aperture NA", ylabel="R (nm)")
    axes[1].set(title="Depth of focus ± 0.5 λ / NA²", xlabel="Numerical aperture NA", ylabel="DOF, one side (nm)")
    for ax in axes:
        ax.set_ylim(bottom=0)
        ax.legend(loc="upper right")
    fig.tight_layout()
    save(fig, "fig1_rayleigh.png")
    write_csv("rayleigh_tools.csv", rows)


def fig2_proximity():
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    # (a) 2 um 슬릿의 간격별 세기 (정확한 각스펙트럼 전파)
    x = np.linspace(-150, 150, 2**17, endpoint=False)
    W = 2.0
    ax = axes[0]
    ax.fill_between([-W / 2, W / 2], 0, 1, color=GRID, lw=0, label="mask opening (2 µm)")
    for g, c in zip((2, 10, 40), (BLUE_RAMP[0], BLUE_RAMP[2], BLUE_RAMP[3])):
        I = np.abs(angular_spectrum_propagate(x, slit_field(x, W), g, I_LINE)) ** 2
        ax.plot(x, I, color=c, label=f"gap {g} µm  (√(λg) = {proximity_min_feature(I_LINE, g):.1f} µm)")
    ax.axhline(0.25, color=TEXT_2, ls=":", lw=1)
    ax.text(5.8, 0.27, "print threshold 0.25", ha="right", color=TEXT_2, fontsize=8)
    ax.set(title="i-line proximity printing of a 2 µm slit", xlabel="Position (µm)",
           ylabel="Intensity (clear field = 1)", xlim=(-6, 6), ylim=(0, 1.8))
    ax.legend(loc="upper left", fontsize=8)

    # (b) CD/W를 W / sqrt(lambda g)로 정리하면 간격과 무관하게 한 곡선이 된다.
    ax = axes[1]
    r = np.linspace(0.6, 5.0, 89)
    rows = []
    for g, c in zip((5, 10, 20, 40), BLUE_RAMP):
        cd = []
        for ri in r:
            Wi = ri * np.sqrt(I_LINE * g)
            I = np.abs(angular_spectrum_propagate(x, slit_field(x, Wi), g, I_LINE)) ** 2
            cd.append(threshold_width(x, I, 0.25) / Wi)
        ax.plot(r, cd, color=c, label=f"gap {g} µm")
        for ri in (1.0, 1.5, 2.0, 3.0):
            k = int(np.argmin(np.abs(r - ri)))
            rows.append({"gap_um": g, "W_over_sqrt_lambda_g": ri, "W_um": round(ri * np.sqrt(I_LINE * g), 3),
                         "printed_width_over_W": round(cd[k], 3)})
    ax.axhspan(0.9, 1.1, color=GRID, lw=0)
    ax.axvline(1.0, color=TEXT_2, ls="--", lw=1)
    ax.text(1.05, 1.62, "W = √(λg)", color=TEXT_2, fontsize=8.5)
    ax.text(4.95, 1.12, "±10 %", ha="right", color=TEXT_2, fontsize=8.5)
    ax.set(title="Printed width collapses onto one curve", xlabel="W / √(λg)",
           ylabel="Printed width / mask width", ylim=(0, 1.8))
    ax.legend(loc="lower right", fontsize=8.5)
    fig.tight_layout()
    save(fig, "fig2_proximity.png")
    write_csv("proximity_printed_width.csv", rows)


def fig3_projection_images():
    cases = [("g-line", 0.436, 0.43, 1.0), ("g-line", 0.436, 0.43, 0.5), ("i-line", 0.365, 0.50, 0.5)]
    fig, axes = plt.subplots(1, 3, figsize=(15, 3.9), sharey=True)
    rows = []
    for ax, (name, lam, NA, hp) in zip(axes, cases):
        p = 2 * hp
        x = np.linspace(-p, p, 1601)
        I = aerial_image(x, p, hp, lam, NA, SIGMA, n_across=61)
        mask = (np.abs(((x + p / 2) % p) - p / 2) > hp / 2).astype(float)
        ax.fill_between(x, 0, mask, color=GRID, lw=0, step="mid", label="mask (clear = 1)")
        ax.plot(x, I, color=TOOL_COLOR[name], label="aerial image")
        con = image_contrast(I)
        ax.set(title=f"{name} NA {NA}, {hp:g} µm L/S: contrast {con:.2f}", xlabel="Position (µm)",
               xlim=(-p, p), ylim=(0, 1.3))
        rows.append({"tool": name, "wavelength_nm": round(lam * 1e3), "NA": NA, "sigma": SIGMA, "half_pitch_um": hp,
                     "k1": round(k1_factor(hp, lam, NA), 3), "contrast": round(con, 3),
                     "resolved_by_iline_resist": con > CMTF_ILINE})
    axes[0].set_ylabel("Intensity (clear field = 1)")
    axes[0].legend(loc="upper right", fontsize=8)
    fig.suptitle("Projection printing, partially coherent (σ = 0.6): cases from Plummer's simulation example",
                 color=TEXT, fontweight="bold")
    fig.tight_layout()
    save(fig, "fig3_projection_images.png")
    write_csv("projection_cases.csv", rows)


def fig4_contrast_vs_pitch():
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    hps = np.linspace(0.1, 1.0, 91)
    rows = []
    for t in TOOLS:
        con = []
        for hp in hps:
            p = 2 * hp
            x = np.linspace(-p / 2, p / 2, 401)
            con.append(image_contrast(aerial_image(x, p, hp, t["wavelength"], t["NA"], SIGMA)))
        con = np.array(con)
        c = TOOL_COLOR[t["name"]]
        ax.plot(hps, con, color=c, label=f"{t['name']}, NA {t['NA']}")
        hp_min = hps[con >= CMTF_ILINE].min()
        R = rayleigh_resolution(t["wavelength"], t["NA"])
        ax.plot(R, np.interp(R, hps, con), marker="o", ms=6, color=c, mec=SURFACE, mew=1.5)
        rows.append({"tool": t["name"], "rayleigh_R_um": round(R, 3),
                     "cutoff_half_pitch_um": round(cutoff_pitch(t["wavelength"], t["NA"], SIGMA) / 2, 3),
                     "min_half_pitch_iline_resist_um": round(hp_min, 3),
                     "k1_at_min": round(k1_factor(hp_min, t["wavelength"], t["NA"]), 3),
                     "min_half_pitch_duv_resist_um": round(hps[con >= CMTF_DUV].min(), 3)})
    for v, lab in ((CMTF_ILINE, f"g/i-line resist CMTF {CMTF_ILINE:.2f} (γ = 2.5)"),
                   (CMTF_DUV, f"DUV resist CMTF {CMTF_DUV:.2f} (γ = 7)")):
        ax.axhline(v, color=TEXT_2, ls="--", lw=1)
        ax.text(0.99, v + 0.02, lab, ha="right", color=TEXT_2, fontsize=8)
    ax.plot([], [], ls="none", marker="o", color=TEXT_2, label="Rayleigh R = 0.61 λ/NA")
    ax.set(title="Image contrast vs. half pitch (dense L/S, σ = 0.6)", xlabel="Half pitch (µm)",
           ylabel="Image contrast (MTF)", xlim=(0.1, 1.0), ylim=(0, 1.05))
    ax.legend(loc="center right", bbox_to_anchor=(1.0, 0.62), fontsize=8.5)
    fig.tight_layout()
    save(fig, "fig4_contrast_vs_pitch.png")
    write_csv("resolution_limits.csv", rows)


def fig5_process_window():
    lam, NA, hp = I_LINE, 0.5, 0.5
    p = 2 * hp
    x = np.linspace(-p / 2, p / 2, 801)
    zs = np.linspace(-1.6, 1.6, 65)
    Is = [aerial_image(x, p, hp, lam, NA, SIGMA, z) for z in zs]
    E0 = dose_to_size(x, Is[len(zs) // 2], hp)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    ax = axes[0]
    doses = [0.90, 0.95, 1.00, 1.05, 1.10]
    for d, c in zip(doses, ["#86b6ef", "#5598e7", "#2a78d6", "#184f95", "#0d366b"]):
        cd = [resist_line_cd(x, I, d * E0) for I in Is]
        ax.plot(zs, np.array(cd) * 1e3, color=c, label=f"dose {d:.2f} E₀")
    ax.axhspan(450, 550, color=GRID, lw=0)
    ax.text(-1.55, 553, "CD spec 500 nm ± 10 %", color=TEXT_2, fontsize=8)
    ax.set(title="Bossung curves: i-line NA 0.5, 0.5 µm L/S", xlabel="Defocus (µm)", ylabel="Line CD (nm)",
           ylim=(300, 600))
    ax.legend(loc="lower center", fontsize=8, ncol=2)

    # 초점별로 CD가 규격 안에 드는 노광량 범위
    ds = np.linspace(0.7, 1.4, 701)
    lo, hi = [], []
    for I in Is:
        ok = np.array([abs(resist_line_cd(x, I, d * E0) - hp) <= 0.1 * hp for d in ds])
        lo.append(ds[ok].min() if ok.any() else np.nan)
        hi.append(ds[ok].max() if ok.any() else np.nan)
    lo, hi = np.array(lo), np.array(hi)
    ax = axes[1]
    ax.fill_between(zs, lo, hi, color="#b7d3f6", lw=0, label="CD within ±10 %")
    ax.plot(zs, lo, color=CAT[0], lw=1.5)
    ax.plot(zs, hi, color=CAT[0], lw=1.5)
    # 노광 여유도(EL)별 DOF: 초점 범위 [-z, z]에서 공통으로 허용되는 노광량 폭
    rows = []
    mid = len(zs) // 2
    for k in range(0, mid + 1):
        sl = slice(mid - k, mid + k + 1)
        a, b = np.nanmax(lo[sl]), np.nanmin(hi[sl])
        el = (b - a) / ((a + b) / 2) * 100 if b > a else 0.0
        rows.append({"focus_range_um": round(2 * zs[mid + k], 3), "dose_window_low": round(a, 3),
                     "dose_window_high": round(b, 3), "exposure_latitude_pct": round(el, 1)})
    dof10 = max(r["focus_range_um"] for r in rows if r["exposure_latitude_pct"] >= 10)
    ray = depth_of_focus(lam, NA)
    ax.axvspan(-ray, ray, color=GRID, alpha=0.6, lw=0)
    ax.axvline(-dof10 / 2, color=CAT[1], ls="--", lw=1.5)
    ax.axvline(dof10 / 2, color=CAT[1], ls="--", lw=1.5, label=f"DOF at 10 % EL: ±{dof10 / 2:.2f} µm")
    ax.plot([], [], lw=8, color=GRID, label=f"Rayleigh ±0.5 λ/NA²: ±{ray:.2f} µm")
    ax.set(title="Process window (ideal optics, thin resist)", xlabel="Defocus (µm)", ylabel="Dose / E₀",
           ylim=(0.7, 1.4), xlim=(zs[0], zs[-1]))
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    save(fig, "fig5_process_window.png")
    write_csv("process_window_el_dof.csv", rows)
    return dof10


def validation():
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    rows = []

    # (a) Abbe 수치 결상 vs 원 겹침 면적 해석해 (0차 + 1차 영역)
    ax = axes[0]
    lam, NA, s = I_LINE, 0.5, SIGMA
    ps = np.linspace(cutoff_pitch(lam, NA, s) * 1.01, lam / NA * 0.99, 30)
    an = [dense_ls_contrast_analytic(p, p / 2, lam, NA, s) for p in ps]
    nu = []
    for p in ps:
        x = np.linspace(-p / 2, p / 2, 801)
        nu.append(image_contrast(aerial_image(x, p, p / 2, lam, NA, s, n_across=81)))
    ax.plot(ps, an, color=CAT[0], label="analytic (circle overlap)")
    ax.plot(ps[::2], nu[::2], ls="none", marker="o", ms=5, mfc=SURFACE, mec=CAT[1], mew=1.3, label="numerical (Abbe)")
    # 차단 피치 근처는 contrast가 0에 가까워 상대 오차가 커지므로 절대 오차(contrast 단위)로 본다.
    err_a = np.max(np.abs(np.array(nu) - an))
    ax.set(title=f"Partial coherence: max |Δcontrast| = {err_a:.4f}", xlabel="Pitch (µm)", ylabel="Image contrast")
    ax.legend(loc="upper left", fontsize=8.5)
    rows.append({"case": "Abbe vs circle-overlap analytic (i-line NA0.5 sigma0.6)", "metric": "max abs contrast diff",
                 "value": round(err_a, 4)})

    # (b) 결맞음 조명 초점 이탈 vs 해석해
    ax = axes[1]
    zs = np.linspace(-2, 2, 81)
    an = [coherent_defocus_contrast_analytic(1.0, 0.5, lam, NA, z) for z in zs]
    x = np.linspace(-0.5, 0.5, 2001)
    nu = [image_contrast(aerial_image(x, 1.0, 0.5, lam, NA, 0.0, z)) for z in zs]
    ax.plot(zs, an, color=CAT[0], label="analytic (3-beam phase)")
    ax.plot(zs[::3], nu[::3], ls="none", marker="o", ms=5, mfc=SURFACE, mec=CAT[1], mew=1.3, label="numerical")
    err_b = np.max(np.abs(np.array(nu) - an) / np.array(an)) * 100
    ax.set(title=f"Coherent defocus: max error {err_b:.4f} %", xlabel="Defocus (µm)", ylabel="Image contrast")
    ax.legend(loc="lower center", fontsize=8.5)
    rows.append({"case": "coherent defocus vs analytic (1 um pitch)", "metric": "max rel error pct", "value": round(err_b, 5)})

    # (c) Fresnel 적분 해석해 vs 각스펙트럼 수치 전파
    ax = axes[2]
    W, g = 2.0, 2.0
    L = 400.0
    N = int(L**2 / (I_LINE * g))
    xs = -L / 2 + np.arange(N) * L / N
    Ip = np.abs(angular_spectrum_propagate(xs, slit_field(xs, W), g, I_LINE, paraxial=True)) ** 2
    xe = np.linspace(-150, 150, 2**17, endpoint=False)
    Ie = np.abs(angular_spectrum_propagate(xe, slit_field(xe, W), g, I_LINE)) ** 2
    xa = np.linspace(-3, 3, 1201)
    ax.plot(xa, slit_fresnel_intensity(xa, W, g, I_LINE), color=CAT[0], label="Fresnel integral (paraxial)")
    m = np.abs(xs) < 3
    ax.plot(xs[m][::40], Ip[m][::40], ls="none", marker="o", ms=4.5, mfc=SURFACE, mec=CAT[1], mew=1.2,
            label="FFT, paraxial")
    me = np.abs(xe) < 3
    ax.plot(xe[me], Ie[me], color=CAT[2], ls="--", lw=1.6, label="FFT, exact (non-paraxial)")
    mm = np.abs(xs) < W + 4
    err_p = np.max(np.abs(Ip[mm] - slit_fresnel_intensity(xs[mm], W, g, I_LINE)))
    mme = np.abs(xe) < W + 4
    err_e = np.max(np.abs(Ie[mme] - slit_fresnel_intensity(xe[mme], W, g, I_LINE)))
    ax.set(title=f"2 µm slit at 2 µm gap: paraxial diff {err_p:.0e}", xlabel="Position (µm)",
           ylabel="Intensity", xlim=(-3, 3), ylim=(0, 2.0))
    ax.legend(loc="upper left", fontsize=8)
    rows.append({"case": "Fresnel integral vs paraxial FFT (W 2um, g 2um)", "metric": "max abs intensity diff",
                 "value": f"{err_p:.1e}"})
    rows.append({"case": "Fresnel approximation vs exact FFT (W 2um, g 2um)", "metric": "max abs intensity diff",
                 "value": round(err_e, 3)})
    for gg in (10.0, 40.0):
        Ie2 = np.abs(angular_spectrum_propagate(xe, slit_field(xe, W), gg, I_LINE)) ** 2
        err = np.max(np.abs(Ie2[mme] - slit_fresnel_intensity(xe[mme], W, gg, I_LINE)))
        rows.append({"case": f"Fresnel approximation vs exact FFT (W 2um, g {gg:g}um)", "metric": "max abs intensity diff",
                     "value": round(err, 3)})
    fig.tight_layout()
    save(fig, "fig6_validation.png")
    write_csv("validation.csv", rows)


def main():
    print("Lithography module")
    fig1_rayleigh()
    fig2_proximity()
    fig3_projection_images()
    fig4_contrast_vs_pitch()
    dof10 = fig5_process_window()
    print(f"  DOF at 10 % EL: +/-{dof10 / 2:.2f} um")
    validation()


if __name__ == "__main__":
    main()
