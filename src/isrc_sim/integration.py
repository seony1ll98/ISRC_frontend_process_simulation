"""공정 통합: 01~05 모듈의 결과를 CMOS 공정 순서대로 연결하고 소자 특성(문턱전압, 유효 채널 길이)을 계산한다.

- 누적 열처리(thermal budget)에 따른 1D 도핑 분포 변화 (02 모듈의 수치 확산 사용)
- MOS 문턱전압: 공핍 근사, 불균일 도핑
- LDD 접합의 측면 확산(2D)과 유효 채널 길이

단위: 깊이와 길이 [um], 농도 [cm^-3], 전압 [V]
"""

import numpy as np
from scipy.optimize import brentq
from scipy.special import erfc

from .constants import Q_E
from .diffusion import diffuse_fd, intrinsic_diffusivity
from .implant import implant_profile, range_parameters

# 300 K 실리콘 상수 (S. M. Sze and K. K. Ng, Physics of Semiconductor Devices, 3rd ed., 2007)
KT_Q_300 = 0.025852          # V
NI_300 = 9.65e9              # cm^-3
EG_300 = 1.12                # eV
EPS0_F_CM = 8.8541878128e-14
EPS_SI = 11.7 * EPS0_F_CM
EPS_OX = 3.9 * EPS0_F_CM
UM = 1e-4


# ---------------------------------------------------------------------------
# 열처리 순서에 따른 1D 도핑 분포
# ---------------------------------------------------------------------------
def apply_thermal_steps(x_um, C, dopant, steps, n_steps_per=200):
    """열처리 단계 [(T_C, t_h), ...]를 차례로 적용한다. 단계마다 D가 다르므로 Dt를 따로 계산한다."""
    for T, t_h in steps:
        Dt = intrinsic_diffusivity(dopant, T) * t_h * 3600.0
        if Dt > 0:
            C = diffuse_fd(x_um, C, Dt, n_steps=n_steps_per)
    return C


def implanted_and_annealed(x_um, dopant, energy_keV, dose, steps):
    """주입 직후 Gaussian 분포에 열처리 단계를 순서대로 적용한 분포."""
    C0 = implant_profile(x_um, dopant, energy_keV, dose)
    return apply_thermal_steps(x_um, C0, dopant, steps)


def junction_depth_between(x_um, C_n, C_p):
    """순 도핑(C_n - C_p)의 부호가 바뀌는 가장 얕은 깊이(표면 쪽 도펀트가 우세한 영역의 끝) [um]."""
    net = np.asarray(C_n) - np.asarray(C_p)
    s = np.sign(net)
    idx = np.nonzero(s[:-1] * s[1:] < 0)[0]
    if idx.size == 0:
        return float("nan")
    i = idx[0]
    return float(x_um[i] + net[i] / (net[i] - net[i + 1]) * (x_um[i + 1] - x_um[i]))


# ---------------------------------------------------------------------------
# MOS 문턱전압 (공핍 근사)
# ---------------------------------------------------------------------------
def oxide_capacitance(tox_nm):
    """C_ox [F/cm^2]."""
    return EPS_OX / (tox_nm * 1e-7)


def fermi_potential(N):
    """phi_F = kT/q ln(N/ni) [V]. N은 순 도핑 농도의 크기."""
    return KT_Q_300 * np.log(N / NI_300)


def work_function_difference(phiF, channel, gate):
    """게이트(n+ 또는 p+ poly)와 채널 기판의 일함수 차 phi_ms [V].

    n+ poly의 페르미 준위는 전도대 끝, p+ poly는 가전자대 끝에 있다고 근사한다.
    기판 페르미 준위는 진성 준위에서 phi_F만큼 떨어져 있다(p형은 아래, n형은 위).
    """
    gate_offset = +EG_300 / 2 if gate == "n+" else -EG_300 / 2      # 게이트 E_F - E_i
    body_offset = -phiF if channel == "n" else +phiF                 # 기판 E_F - E_i
    return -(gate_offset - body_offset)


def vt_uniform(N_body, tox_nm, channel="n", Qf_cm2=0.0, gate="n+"):
    """균일 도핑 기판의 장채널 문턱전압 (해석해, 검증용).

    gate는 "n+" 또는 "p+" poly다.
      NMOS (p형 기판): Vt = V_FB + 2 phi_F + sqrt(2 eps q N 2 phi_F) / C_ox
      PMOS (n형 기판): Vt = V_FB - 2 phi_F - sqrt(2 eps q N 2 phi_F) / C_ox
      V_FB = phi_ms - q Q_f / C_ox
    """
    Cox = oxide_capacitance(tox_nm)
    phiF = fermi_potential(N_body)
    Qd = np.sqrt(2 * EPS_SI * Q_E * N_body * 2 * phiF)
    Vfb = work_function_difference(phiF, channel, gate) - Q_E * Qf_cm2 / Cox
    if channel == "n":
        return Vfb + 2 * phiF + Qd / Cox
    return Vfb - 2 * phiF - Qd / Cox


def vt_nonuniform(x_um, N_net, tox_nm, channel="n", N_bulk=None, Qf_cm2=0.0, gate="n+"):
    """불균일 도핑 기판의 문턱전압 (공핍 근사).

    N_net: 채널 반도체의 순 도핑(다수 캐리어 쪽이 양수) [cm^-3]
    표면 전위가 psi_s = 2 phi_F(N_bulk)가 될 때를 문턱으로 보고, 공핍 폭 W는
      psi_s = (q / eps) * int_0^W x N(x) dx
    를 만족하도록 구한다. 공핍 전하 Q_d = q int_0^W N dx 이고
      Vt = V_FB +/- (psi_s + Q_d / C_ox) 이다.
    N_bulk가 None이면 분포의 가장 깊은 값을 기판 농도로 쓴다.
    """
    x = np.asarray(x_um, dtype=float) * UM
    N = np.clip(np.asarray(N_net, dtype=float), 0, None)
    if N_bulk is None:
        N_bulk = N[-1]
    phiF = fermi_potential(N_bulk)
    psi_t = 2 * phiF
    xN = np.concatenate([[0.0], np.cumsum(0.5 * (x[1:] * N[1:] + x[:-1] * N[:-1]) * np.diff(x))])
    Nint = np.concatenate([[0.0], np.cumsum(0.5 * (N[1:] + N[:-1]) * np.diff(x))])
    psi = Q_E / EPS_SI * xN
    if psi[-1] < psi_t:
        raise ValueError("계산 영역이 공핍 폭보다 얕습니다. 깊이를 늘리세요.")
    W = float(np.interp(psi_t, psi, x))
    Qd = Q_E * float(np.interp(W, x, Nint))
    Cox = oxide_capacitance(tox_nm)
    Vfb = work_function_difference(phiF, channel, gate) - Q_E * Qf_cm2 / Cox
    if channel == "n":
        return Vfb + psi_t + Qd / Cox, W / UM
    return Vfb - psi_t - Qd / Cox, W / UM


# ---------------------------------------------------------------------------
# LDD 접합의 측면 확산 (2D)
# ---------------------------------------------------------------------------
def ldd_2d_analytic(x_um, y_um, dopant, energy_keV, dose, Dt_cm2, lateral_straggle_um=None):
    """게이트 모서리(x = 0, 게이트는 x > 0)로 가려진 주입의 2D 분포와 열처리 후 분포.

    초기 분포를 깊이 방향 Gaussian x 측면 계단(측면 straggle로 번짐)의 곱으로 두면, D가 일정할 때
    확산 후에도 곱의 형태가 유지된다.
      깊이: 분산 dRp^2 + 2Dt (표면 반사 경계는 영상 항으로)
      측면: (1/2) erfc(x / sqrt(2 (dR_lat^2 + 2Dt)))
    측면 straggle을 주지 않으면 dRp와 같다고 둔다(Plummer 표에 측면 값이 없어 둔 가정).
    """
    Rp, dRp = range_parameters(dopant, energy_keV)
    lat = dRp if lateral_straggle_um is None else lateral_straggle_um
    X, Y = np.meshgrid(np.asarray(x_um, dtype=float), np.asarray(y_um, dtype=float), indexing="xy")
    var_y = (dRp * UM) ** 2 + 2 * Dt_cm2
    var_x = (lat * UM) ** 2 + 2 * Dt_cm2
    yy = Y * UM
    vert = dose / np.sqrt(2 * np.pi * var_y) * (np.exp(-(yy - Rp * UM) ** 2 / (2 * var_y))
                                               + np.exp(-(yy + Rp * UM) ** 2 / (2 * var_y)))
    side = 0.5 * erfc((X * UM) / np.sqrt(2 * var_x)) if var_x > 0 else (X <= 0).astype(float)
    return vert * side


def diffuse_2d_fd(x_um, y_um, C0, Dt_cm2, n_steps=200):
    """2D 확산을 깊이 방향과 측면 방향 1D Crank-Nicolson으로 번갈아 푼다 (검증용).

    D가 일정하면 두 방향 연산자가 교환 가능하므로 방향 분할에 따른 시간 오차가 없다.
    경계는 모두 반사 경계다. C0의 행은 깊이(y), 열은 측면(x)이다.
    """
    C = np.asarray(C0, dtype=float).copy()
    for j in range(C.shape[1]):
        C[:, j] = diffuse_fd(y_um, C[:, j], Dt_cm2, n_steps=n_steps)
    for i in range(C.shape[0]):
        C[i, :] = diffuse_fd(x_um - x_um[0], C[i, :], Dt_cm2, n_steps=n_steps)
    return C


def lateral_junction_extent(x_um, y_um, C_ldd, N_body_of_y):
    """게이트 아래(x > 0)로 접합이 가장 멀리 들어간 거리 [um]."""
    Nb = np.asarray(N_body_of_y)[:, None]
    n_type = C_ldd > Nb
    reach = 0.0
    for i in range(C_ldd.shape[0]):
        cols = np.nonzero(n_type[i])[0]
        if cols.size == 0:
            continue
        j = cols.max()
        if j + 1 < C_ldd.shape[1]:
            c0, c1 = C_ldd[i, j] - Nb[i, 0], C_ldd[i, j + 1] - Nb[i, 0]
            xj = x_um[j] + c0 / (c0 - c1) * (x_um[j + 1] - x_um[j])
        else:
            xj = x_um[j]
        reach = max(reach, float(xj))
    return reach
