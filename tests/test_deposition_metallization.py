import numpy as np
import pytest

from isrc_sim.deposition import (
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
    view_factor_segments,
    view_factor_to_opening,
)
from isrc_sim.metallization import (
    INFINITE_SHEET_CF,
    current_density_MA_cm2,
    four_point_probe_cf_disk,
    four_point_probe_fd,
    line_resistance_capacitance,
    rc_delay,
    resistivity_from_sheet,
    sheet_from_resistivity,
    smits_center_cf,
)


# ---------------------------------------------------------------------------
# CVD (Grove)
# ---------------------------------------------------------------------------
def test_grove_limits():
    hG = 1.0
    assert grove_growth_rate(1e-4, hG) == pytest.approx(1e-4, rel=1e-3)   # 반응 율속: v ~ ks
    assert grove_growth_rate(1e4, hG) == pytest.approx(1.0, rel=1e-3)     # 물질 전달 율속: v ~ hG
    assert grove_growth_rate(1.0, hG) == pytest.approx(0.5)


def test_transition_temperature_and_shift():
    ks = surface_rate_constant(1000, 1.6, 1000, hG=1.0)
    assert ks == pytest.approx(1.0)
    T2 = transition_temperature_shift(1000, 1.6, 100)
    assert surface_rate_constant(T2, 1.6, 1000, hG=1.0) == pytest.approx(100.0, rel=1e-9)
    assert T2 > 1000


def test_lpcvd_gain_matches_textbook():
    # Plummer: 1 torr에서 D_G는 760배, 경계층은 약 7배 -> h_G 약 100배
    assert lpcvd_hG_gain(1.0) == pytest.approx(108.6, rel=1e-3)


def test_temperature_sensitivity():
    assert temperature_sensitivity(1000, 1.6) == pytest.approx(1.6 / (8.617333262e-5 * 1273.15**2))


# ---------------------------------------------------------------------------
# 단차 피복
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def trench():
    el = trench_elements(1.0, 2.0, n_per_side=60)
    return el, view_factor_matrix(el), view_factor_to_opening(el, 1.0)


def test_view_factor_parallel_plates():
    # 마주 보는 길이 1, 간격 1인 두 판: F = sqrt(2) - 1
    assert view_factor_segments((0, 0), (1, 0), (1, 1), (0, 1)) == pytest.approx(np.sqrt(2) - 1)


def test_enclosure_sum_rule(trench):
    el, V, vo = trench
    assert V.sum(axis=1) + vo == pytest.approx(np.ones(len(vo)), abs=1e-12)


def test_isotropic_direct_flux_matches_view_factor(trench):
    el, V, vo = trench
    assert direct_flux(el, 1.0, n=1) == pytest.approx(vo, abs=2e-3)


def test_bottom_center_analytic():
    # 바닥 중앙의 미소 요소가 입구를 보는 view factor = sin(atan(w / 2h))
    w, h = 1.0, 2.0
    el = trench_elements(w, h, n_per_side=80, n_bottom=81)
    F, _, _ = solve_trench_flux(el, w, n=1, Sc=1.0)
    assert step_coverage_metrics(el, F, 1.0)["bottom_center"] == pytest.approx(np.sin(np.arctan(w / (2 * h))), rel=2e-3)


@pytest.mark.parametrize("Sc", [1.0, 0.3, 0.05, 0.005])
@pytest.mark.parametrize("n", [1, 4, 40])
def test_particle_balance(trench, Sc, n):
    # 입구로 들어온 입자 = 붙은 입자 + 다시 빠져나간 입자
    el, V, vo = trench
    F, _, _ = solve_trench_flux(el, 1.0, n=n, Sc=Sc, V=V)
    stuck = np.sum(Sc * F * el["length"])
    escaped = np.sum((1 - Sc) * F * el["length"] * vo)
    assert stuck + escaped == pytest.approx(1.0, abs=2e-3)


def test_low_sticking_gives_conformal_coverage(trench):
    el, V, _ = trench
    F, _, _ = solve_trench_flux(el, 1.0, n=1, Sc=1e-4, V=V)
    assert F.min() > 0.99


def test_lower_sticking_is_more_conformal(trench):
    el, V, _ = trench
    mins = [solve_trench_flux(el, 1.0, n=1, Sc=s, V=V)[0].min() for s in (1.0, 0.3, 0.05, 0.005)]
    assert np.all(np.diff(mins) > 0)


def test_directional_arrival_fills_bottom_but_starves_sidewall(trench):
    el, V, _ = trench
    m1 = step_coverage_metrics(el, solve_trench_flux(el, 1.0, n=1, Sc=1.0, V=V)[0], 1.0)
    m80 = step_coverage_metrics(el, solve_trench_flux(el, 1.0, n=80, Sc=1.0, V=V)[0], 1.0)
    assert m80["bottom_center"] > 0.95 > m1["bottom_center"]
    assert m80["sidewall_mid"] < m1["sidewall_mid"]


# ---------------------------------------------------------------------------
# 4-point probe
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("D_over_s", [5, 10, 20, 40, 100])
def test_image_method_equals_smits_at_center(D_over_s):
    assert four_point_probe_cf_disk(D_over_s, 1.0) == pytest.approx(smits_center_cf(D_over_s, 1.0), rel=1e-10)


def test_large_wafer_approaches_infinite_sheet():
    assert INFINITE_SHEET_CF == pytest.approx(4.532, abs=1e-3)
    assert four_point_probe_cf_disk(1e4, 1.0) == pytest.approx(INFINITE_SHEET_CF, rel=1e-6)


@pytest.mark.parametrize("D_over_s,x0", [(10, 0.0), (20, 0.0), (20, 6.0), (40, 17.0)])
def test_finite_difference_matches_image_method(D_over_s, x0):
    assert four_point_probe_fd(D_over_s, 1.0, x0, nodes_per_spacing=20) == pytest.approx(
        four_point_probe_cf_disk(D_over_s, 1.0, x0), rel=3e-3)


def test_cf_drops_toward_edge():
    cf = [four_point_probe_cf_disk(40, 1.0, x0) for x0 in (0, 5, 10, 15, 18)]
    assert np.all(np.diff(cf) < 0)


def test_resistivity_conversion_roundtrip():
    assert sheet_from_resistivity(3.0, 500) == pytest.approx(0.06)
    assert resistivity_from_sheet(0.06, 500) == pytest.approx(3.0)


# ---------------------------------------------------------------------------
# RC 지연
# ---------------------------------------------------------------------------
def test_rc_delay_equals_0p89_RC():
    R, C = line_resistance_capacitance(3.0, 10, 1, 1, 1, 1, 3.9)
    assert rc_delay(3.0, 10, 1, 1, 1, 1, 3.9) == pytest.approx(0.89 * R * C)


def test_rc_delay_scales_with_length_squared():
    assert rc_delay(3.0, 20, 1, 1, 1, 1, 3.9) / rc_delay(3.0, 10, 1, 1, 1, 1, 3.9) == pytest.approx(4.0)


def test_current_density():
    assert current_density_MA_cm2(5.0, 1.0, 1.0) == pytest.approx(0.5)
