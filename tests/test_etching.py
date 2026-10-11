import numpy as np
import pytest

from isrc_sim.etching import (
    MASK,
    UNDER,
    anisotropy,
    anisotropy_from_bias,
    clear_time_window,
    directional_flux_factor,
    make_structure,
    profile_metrics,
    rate_linear,
    rate_saturation,
    required_mask_thickness,
    required_selectivity,
    selectivity,
    simulate_profile,
    underlayer_loss,
)

DX = 0.005


# ---------------------------------------------------------------------------
# 지표와 설계식
# ---------------------------------------------------------------------------
def test_anisotropy_definitions():
    assert anisotropy(1.0, 1.0) == 0.0                   # 등방성
    assert anisotropy(0.0, 1.0) == 1.0                   # 완전 이방성
    assert anisotropy_from_bias(0.1, 0.5) == pytest.approx(0.8)
    assert selectivity(0.5, 0.02) == pytest.approx(25)


def test_clear_time_window():
    t1, t2 = clear_time_window(0.3, 0.1, 0.05, 0.05)
    assert t1 == pytest.approx(0.3 * 0.95 / 0.105)
    assert t2 == pytest.approx(0.3 * 1.05 / 0.095)


def test_underlayer_loss_scaling_and_required_selectivity():
    args = (0.3, 0.1, 0.05, 0.05, 0.2)
    assert underlayer_loss(0.3, 0.1, 10, 0.05, 0.05, 0.2) == pytest.approx(2 * underlayer_loss(0.3, 0.1, 20, 0.05, 0.05, 0.2))
    S = required_selectivity(0.3, 0.005, 0.05, 0.05, 0.2)
    assert underlayer_loss(0.3, 0.1, S, 0.05, 0.05, 0.2) == pytest.approx(0.005)
    # 식각 속도와 무관하다
    assert underlayer_loss(0.3, 0.5, S, 0.05, 0.05, 0.2) == pytest.approx(0.005)
    assert underlayer_loss(0.3, 0.1, np.inf, *args[2:]) == 0.0


def test_required_mask_thickness():
    # 마스크가 정확히 0이 되는 두께 = 최대 식각 속도 x 공정 시간 / 선택비
    _, t_last = clear_time_window(0.3, 1.0, 0.05, 0.05)
    assert required_mask_thickness(0.3, 4, 0.05, 0.05, 0.2) == pytest.approx(t_last * 1.2 * 1.05 / 4)


# ---------------------------------------------------------------------------
# 식각 속도 모델
# ---------------------------------------------------------------------------
def test_saturation_model_needs_both_species():
    assert rate_saturation(1.0, 0.0) == 0.0
    assert rate_saturation(0.0, 1.0) == 0.0
    assert rate_saturation(1.0, 1.0) == pytest.approx(0.5)
    ion = np.logspace(-2, 3, 50)
    r = rate_saturation(1.0, ion)
    assert np.all(r <= np.minimum(1.0, ion) + 1e-12)
    assert r[-1] == pytest.approx(1.0, rel=1e-2)          # 이온이 충분하면 화학 성분으로 포화


def test_linear_model_is_sum():
    assert rate_linear(0.3, 0.2) == pytest.approx(0.5)


# ---------------------------------------------------------------------------
# 이온 입사각 분포
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("tilt_deg", [0, 30, 60, 90, 135, 180])
def test_flux_factor_cos1_matches_view_factor(tilt_deg):
    a = np.radians(tilt_deg)
    assert directional_flux_factor(a, 1) == pytest.approx((1 + np.cos(a)) / 2, abs=2e-3)


def test_flux_factor_narrow_distribution_approaches_cosine():
    for deg in (0, 20, 45, 70):
        a = np.radians(deg)
        assert directional_flux_factor(a, 400) == pytest.approx(np.cos(a), abs=0.01)
    # 측벽(90도)이 받는 이온은 분포가 좁을수록 줄어든다
    s = [directional_flux_factor(np.pi / 2, n) for n in (10, 20, 40, 80)]
    assert np.all(np.diff(s) < 0)


# ---------------------------------------------------------------------------
# 단면 형상 (level set)
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def structure():
    return make_structure(0.5, 0.5, under_um=0.4, width_um=2.0, dx_um=DX)


def test_isotropic_profile_matches_analytic(structure):
    x, y, mat = structure
    r = 0.4
    _, snaps = simulate_profile(x, y, mat, R_chem=0.1, R_ion=0.0, t_end=r / 0.1)
    void = snaps[-1] & (mat != MASK)
    X, Y = np.meshgrid(x, y)
    # 정확해: 개구(|x| <= w/2, y = 0)에서 거리 r 이내인 막 영역
    dist = np.hypot(np.clip(np.abs(X) - 0.25, 0, None), Y)
    exact = (Y < 0) & (dist <= r)
    sim = void & (Y < 0)
    mean_shift = np.sum(sim ^ exact) * DX**2 / (0.5 + np.pi * r)
    assert mean_shift < 2 * DX


def _isotropic_error(dx):
    x, y, mat = make_structure(0.5, 0.5, under_um=0.4, width_um=2.0, dx_um=dx)
    _, snaps = simulate_profile(x, y, mat, R_chem=0.1, R_ion=0.0, t_end=4.0, n_snapshots=2)
    X, Y = np.meshgrid(x, y)
    exact = (Y < 0) & (np.hypot(np.clip(np.abs(X) - 0.25, 0, None), Y) <= 0.4)
    return np.sum((snaps[-1] & (Y < 0)) ^ exact) * dx**2 / (0.5 + np.pi * 0.4)


def test_level_set_converges_with_grid():
    # 1차 정확도: 격자를 절반으로 줄이면 오차도 대략 절반이 된다.
    assert _isotropic_error(0.02) / _isotropic_error(0.01) > 1.5


def test_ion_only_profile_is_vertical(structure):
    x, y, mat = structure
    _, snaps = simulate_profile(x, y, mat, R_chem=0.0, R_ion=0.1, t_end=4.0)
    m = profile_metrics(x, y, snaps[-1])
    assert m["depth_um"] == pytest.approx(0.4, abs=2 * DX)
    assert m["half_width_um"] == pytest.approx(0.25, abs=2 * DX)


@pytest.mark.parametrize("Rc,Ri", [(0.05, 0.05), (0.02, 0.08)])
def test_linear_model_depth_and_bias(structure, Rc, Ri):
    x, y, mat = structure
    t = 4.0
    _, snaps = simulate_profile(x, y, mat, R_chem=Rc, R_ion=Ri, t_end=t)
    m = profile_metrics(x, y, snaps[-1], at_depth_um=0.1)
    # 정확해: (개구 x [-R_i t, 0] 직사각형) + 반지름 R_c t 원의 Minkowski 합 -> 측벽 위치 = w/2 + R_c t
    assert m["depth_um"] == pytest.approx((Rc + Ri) * t, abs=2 * DX)
    assert m["half_width_um"] - 0.25 == pytest.approx(Rc * t, abs=3 * DX)


def test_saturation_model_with_vertical_ions_has_no_undercut(structure):
    x, y, mat = structure
    _, snaps = simulate_profile(x, y, mat, R_chem=0.1, R_ion=0.1, t_end=4.0, model="saturation")
    m = profile_metrics(x, y, snaps[-1], at_depth_um=0.1)
    assert m["depth_um"] == pytest.approx(0.05 * 4.0, abs=2 * DX)
    assert m["half_width_um"] == pytest.approx(0.25, abs=2 * DX)


def test_mask_is_not_etched_and_underlayer_selectivity(structure):
    x, y, mat = structure
    # 막 0.5 um를 0.1 um/min으로 뚫은 뒤 1 min 더 식각, 하부막 선택비 10 -> 0.01 um 손실
    _, snaps = simulate_profile(x, y, mat, R_chem=0.0, R_ion=0.1, t_end=6.0, sel_under=10)
    void = snaps[-1]
    assert not np.any(void & (mat == MASK))
    m = profile_metrics(x, y, void)
    assert m["depth_um"] - 0.5 == pytest.approx(0.01, abs=2 * DX)
    assert np.any(void & (mat == UNDER))
