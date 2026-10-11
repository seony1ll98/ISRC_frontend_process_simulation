import numpy as np
import pytest

from isrc_sim.constants import Q_E
from isrc_sim.integration import (
    diffuse_2d_fd,
    junction_depth_between,
    ldd_2d_analytic,
    oxide_capacitance,
    vt_nonuniform,
    vt_uniform,
    work_function_difference,
)
from isrc_sim.diffusion import intrinsic_diffusivity
from isrc_sim.oxidation import deal_grove_thickness, massoud_thickness

X = np.linspace(0, 3, 30001)


# ---------------------------------------------------------------------------
# 교재 공정 흐름(Plummer Ch. 2)의 대표 조건 재현
# ---------------------------------------------------------------------------
def test_locos_field_oxide_matches_textbook():
    # 90 min @ 1000 C, H2O -> 약 0.5 um (패드 산화막 40 nm 위에서 성장)
    assert deal_grove_thickness(1.5, 1000, "wet", x_i_um=0.04) == pytest.approx(0.5, rel=0.1)


def test_gate_oxide_matches_textbook():
    # 0.5~1 h @ 800 C, O2 -> 3~5 nm
    assert 2.5 < massoud_thickness(0.5, 800) * 1e3 < 4.0
    assert 4.5 < massoud_thickness(1.0, 800) * 1e3 < 6.0


# ---------------------------------------------------------------------------
# 문턱전압
# ---------------------------------------------------------------------------
def test_work_function_difference_signs():
    phiF = 0.4
    assert work_function_difference(phiF, "n", "n+") == pytest.approx(-(0.56 + 0.4))
    assert work_function_difference(phiF, "p", "p+") == pytest.approx(0.56 + 0.4)
    assert work_function_difference(phiF, "p", "n+") == pytest.approx(-(0.56 - 0.4))


@pytest.mark.parametrize("N", [1e15, 1e16, 1e17, 5e17])
@pytest.mark.parametrize("channel,gate", [("n", "n+"), ("p", "p+"), ("p", "n+")])
def test_nonuniform_solver_matches_uniform_analytic(N, channel, gate):
    v, _ = vt_nonuniform(X, np.full_like(X, N), 4.0, channel, gate=gate)
    assert v == pytest.approx(vt_uniform(N, 4.0, channel, gate=gate), abs=1e-6)


def test_dual_poly_is_symmetric():
    assert vt_uniform(1e17, 4.0, "n", gate="n+") == pytest.approx(-vt_uniform(1e17, 4.0, "p", gate="p+"))


def _with_sheet(N, Q, depth_um):
    """균일 도핑 N에 깊이 depth에 면밀도 Q인 얇은 층을 더한다(격자 한 칸에 담는다)."""
    out = np.full_like(X, N)
    dx = (X[1] - X[0]) * 1e-4
    out[int(round(depth_um / (X[1] - X[0])))] += Q / dx
    return out


def test_surface_implant_shift_equals_qQ_over_Cox():
    N, Q = 1e16, 1e11
    v0, _ = vt_nonuniform(X, np.full_like(X, N), 4.0, "n", N_bulk=N)
    v1, _ = vt_nonuniform(X, _with_sheet(N, Q, 0.0) + 0 * X, 4.0, "n", N_bulk=N)
    # 표면 격자점은 사다리꼴 적분에서 절반만 반영되므로 2배를 넣은 효과와 같다
    assert (v1 - v0) == pytest.approx(0.5 * Q_E * Q / oxide_capacitance(4.0), rel=1e-3)


@pytest.mark.parametrize("depth", [0.02, 0.05, 0.1])
def test_buried_implant_efficiency_factor(depth):
    # 공핍 영역 안 깊이 x_i의 얇은 층: dVt = (qQ/Cox)(1 - x_i/W)  (작은 Q의 선형 근사)
    N, Q = 1e16, 1e10
    v0, W = vt_nonuniform(X, np.full_like(X, N), 4.0, "n", N_bulk=N)
    v1, _ = vt_nonuniform(X, _with_sheet(N, Q, depth), 4.0, "n", N_bulk=N)
    expected = Q_E * Q / oxide_capacitance(4.0) * (1 - depth / W)
    assert (v1 - v0) == pytest.approx(expected, rel=2e-2)


def test_vt_increases_with_oxide_thickness_for_nmos():
    assert vt_uniform(1e17, 5.0, "n") > vt_uniform(1e17, 4.0, "n")


# ---------------------------------------------------------------------------
# 접합과 2D 측면 확산
# ---------------------------------------------------------------------------
def test_junction_depth_between_simple():
    x = np.linspace(0, 1, 1001)
    assert junction_depth_between(x, 1e17 * np.exp(-x / 0.1), np.full_like(x, 1e15)) == pytest.approx(
        0.1 * np.log(100), rel=1e-3)


def test_ldd_2d_analytic_matches_2d_numerical():
    x = np.linspace(-0.4, 0.4, 321)
    y = np.linspace(0, 0.5, 201)
    Dt = intrinsic_diffusivity("P", 900) * 1800
    C0 = ldd_2d_analytic(x, y, "P", 50, 5e13, 0.0)
    C_fd = diffuse_2d_fd(x, y, C0, Dt, n_steps=200)
    C_an = ldd_2d_analytic(x, y, "P", 50, 5e13, Dt)
    m = C_an > 1e-3 * C_an.max()
    assert np.max(np.abs(C_fd[m] - C_an[m]) / C_an.max()) < 5e-3


def test_ldd_lateral_profile_is_half_at_mask_edge():
    x = np.array([-1.0, 0.0, 1.0])
    y = np.linspace(0, 0.3, 301)
    C = ldd_2d_analytic(x, y, "P", 50, 5e13, 1e-12)
    assert C[:, 1] == pytest.approx(0.5 * C[:, 0], rel=1e-6)
