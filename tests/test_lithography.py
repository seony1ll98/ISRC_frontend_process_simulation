import numpy as np
import pytest

from isrc_sim.lithography import (
    aerial_image,
    angular_spectrum_propagate,
    coherent_defocus_contrast_analytic,
    cutoff_pitch,
    dense_ls_contrast_analytic,
    depth_of_focus,
    dose_to_size,
    fresnel_number,
    image_contrast,
    proximity_min_feature,
    rayleigh_resolution,
    resist_cmtf,
    resist_gamma,
    resist_line_cd,
    slit_field,
    slit_fresnel_intensity,
    threshold_width,
)

I_LINE = 0.365


# ---------------------------------------------------------------------------
# 교재 예제
# ---------------------------------------------------------------------------
def test_plummer_krf_example():
    # KrF, NA 0.6 -> R ~ 0.3 um (k1 = 0.75), DOF ~ +/-0.35 um (k2 = 0.5)
    assert rayleigh_resolution(0.248, 0.6, k1=0.75) == pytest.approx(0.3, rel=0.05)
    assert depth_of_focus(0.248, 0.6, k2=0.5) == pytest.approx(0.35, rel=0.05)


def test_plummer_proximity_example():
    # i-line, g = 10 um -> W_min ~ 2 um
    assert proximity_min_feature(I_LINE, 10.0) == pytest.approx(2.0, rel=0.1)


def test_resist_cmtf_matches_textbook_ranges():
    # g/i-line gamma 2~3 -> CMTF ~ 0.4,  DUV gamma 5~10 -> CMTF 0.1~0.2
    assert resist_cmtf(2.5) == pytest.approx(0.4, abs=0.05)
    assert 0.1 <= resist_cmtf(10) < resist_cmtf(5) <= 0.23
    assert resist_gamma(10, 100) == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# 근접 노광 (Fresnel)
# ---------------------------------------------------------------------------
def test_wide_slit_center_and_edge_intensity():
    # 아주 넓은 슬릿: 중앙은 1, 기하학적 가장자리는 반평면 회절 값 0.25
    W, g = 400.0, 10.0
    assert slit_fresnel_intensity(0.0, W, g, I_LINE) == pytest.approx(1.0, abs=0.02)
    assert slit_fresnel_intensity(W / 2, W, g, I_LINE) == pytest.approx(0.25, abs=0.01)


def test_fresnel_conserves_power():
    W, g = 2.0, 10.0
    x = np.linspace(-300, 300, 600001)
    I = slit_fresnel_intensity(x, W, g, I_LINE)
    assert np.trapezoid(I, x) == pytest.approx(W, rel=0.01)


@pytest.mark.parametrize("W,g", [(2.0, 10.0), (2.0, 2.0), (5.0, 20.0)])
def test_paraxial_angular_spectrum_matches_fresnel_integral(W, g):
    # 전달함수 방법은 임계 샘플링 dx = lambda g / L에서 aliasing 없이 정확하다.
    L = 400.0
    N = int(L**2 / (I_LINE * g))
    x = -L / 2 + np.arange(N) * L / N
    I_num = np.abs(angular_spectrum_propagate(x, slit_field(x, W), g, I_LINE, paraxial=True)) ** 2
    I_an = slit_fresnel_intensity(x, W, g, I_LINE)
    m = np.abs(x) < W + 4
    assert np.max(np.abs(I_num[m] - I_an[m])) < 1e-3


def test_proximity_cd_depends_only_on_fresnel_number():
    # 같은 W / sqrt(lambda g)이면 간격이 달라도 CD/W가 같다.
    x = np.linspace(-150, 150, 2**17, endpoint=False)
    ratios = []
    for g in (5.0, 40.0):
        W = np.sqrt(I_LINE * g)
        I = np.abs(angular_spectrum_propagate(x, slit_field(x, W), g, I_LINE)) ** 2
        ratios.append(threshold_width(x, I, 0.25) / W)
    assert ratios[0] == pytest.approx(ratios[1], rel=0.02)
    assert fresnel_number(np.sqrt(I_LINE * 10), I_LINE, 10) == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# 축소 투영 노광 (Abbe)
# ---------------------------------------------------------------------------
def test_open_mask_gives_unit_intensity():
    x = np.linspace(-0.5, 0.5, 101)
    I = aerial_image(x, 1.0, 0.0, I_LINE, 0.5, sigma=0.6)
    assert I == pytest.approx(np.ones_like(x))


@pytest.mark.parametrize("pitch", [0.50, 0.55, 0.60, 0.70])
def test_abbe_matches_circle_overlap_analytic(pitch):
    x = np.linspace(-pitch / 2, pitch / 2, 801)
    I = aerial_image(x, pitch, pitch / 2, I_LINE, 0.5, sigma=0.6, n_across=81)
    assert image_contrast(I) == pytest.approx(dense_ls_contrast_analytic(pitch, pitch / 2, I_LINE, 0.5, 0.6), rel=0.01)


def test_no_image_beyond_cutoff_pitch():
    p = 0.98 * cutoff_pitch(I_LINE, 0.5, 0.6)
    x = np.linspace(-p / 2, p / 2, 401)
    assert image_contrast(aerial_image(x, p, p / 2, I_LINE, 0.5, sigma=0.6)) < 1e-9


@pytest.mark.parametrize("z", [0.0, 0.3, 0.8, 1.2])
def test_coherent_defocus_matches_analytic(z):
    p = 1.0
    x = np.linspace(-p / 2, p / 2, 2001)
    I = aerial_image(x, p, p / 2, I_LINE, 0.5, sigma=0.0, defocus_um=z)
    assert image_contrast(I) == pytest.approx(coherent_defocus_contrast_analytic(p, p / 2, I_LINE, 0.5, z), rel=1e-5)


def test_defocus_is_symmetric_without_aberration():
    x = np.linspace(-0.5, 0.5, 401)
    a = aerial_image(x, 1.0, 0.5, I_LINE, 0.5, sigma=0.6, defocus_um=0.7)
    b = aerial_image(x, 1.0, 0.5, I_LINE, 0.5, sigma=0.6, defocus_um=-0.7)
    assert a == pytest.approx(b)


def test_contrast_drops_with_defocus_and_smaller_pitch():
    def c(p, z=0.0):
        x = np.linspace(-p / 2, p / 2, 401)
        return image_contrast(aerial_image(x, p, p / 2, I_LINE, 0.5, sigma=0.6, defocus_um=z))
    assert c(1.0) > c(1.0, 1.0)
    assert c(1.0) > c(0.7) > c(0.55)


# ---------------------------------------------------------------------------
# 레지스트 문턱값 모델
# ---------------------------------------------------------------------------
def test_dose_to_size_and_cd_trend():
    p = 1.0
    x = np.linspace(-p / 2, p / 2, 801)
    I = aerial_image(x, p, 0.5, I_LINE, 0.5, sigma=0.6)
    E0 = dose_to_size(x, I, 0.5)
    assert resist_line_cd(x, I, E0) == pytest.approx(0.5, abs=1e-4)
    # 양성 레지스트: 노광량이 늘면 남는 선이 가늘어진다.
    assert resist_line_cd(x, I, 1.1 * E0) < 0.5 < resist_line_cd(x, I, 0.9 * E0)
