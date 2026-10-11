"""도펀트 확산(열처리) 모델.

- 진성 확산계수 (Plummer)
- 해석해: 주입 Gaussian의 열처리, 표면 delta 소스 drive-in, 일정 표면 농도 predeposition
- 수치해: Crank-Nicolson 유한차분 (해석해 검증과 임의 초기 분포용)

공통 단위: 깊이 [um], 시간 [h], 온도 [C], 농도 [cm^-3], Dt [cm^2]
"""

import numpy as np
from scipy.linalg import solve_banded
from scipy.special import erf, erfc

from .constants import INTRINSIC_DIFFUSIVITY, K_B_EV, celsius_to_kelvin

UM_TO_CM = 1e-4
H_TO_S = 3600.0


def intrinsic_diffusivity(dopant, T_C):
    """진성 확산계수 D [cm^2/s]."""
    if dopant not in INTRINSIC_DIFFUSIVITY:
        raise ValueError(f"지원하는 도펀트는 {sorted(INTRINSIC_DIFFUSIVITY)}입니다: {dopant!r}")
    D0, Ea = INTRINSIC_DIFFUSIVITY[dopant]
    return D0 * np.exp(-Ea / (K_B_EV * celsius_to_kelvin(T_C)))


def thermal_budget(dopant, steps):
    """여러 열처리 단계의 유효 Dt [cm^2] = sum(D_i t_i).

    steps: [(T_C, t_h), ...]
    """
    return float(sum(intrinsic_diffusivity(dopant, T) * t * H_TO_S for T, t in steps))


# ---------------------------------------------------------------------------
# 해석해
# ---------------------------------------------------------------------------
def annealed_gaussian(x_um, dose_cm2, Rp_um, dRp_um, Dt_cm2):
    """주입 Gaussian을 열처리한 분포 [cm^-3].

    분산이 dRp^2 + 2Dt로 넓어진다(Plummer Ch. 8, Eq. 7).
    표면을 반사 경계(외부 확산 없음)로 두기 위해 x = -Rp에 영상(image) 항을 더한다.
    이 해는 주입 분포 중 표면 위로 넘어가는 꼬리가 무시할 만큼 작을 때 정확하다.
    """
    x = np.asarray(x_um, dtype=float) * UM_TO_CM
    Rp = Rp_um * UM_TO_CM
    var = (dRp_um * UM_TO_CM) ** 2 + 2.0 * Dt_cm2
    pref = dose_cm2 / np.sqrt(2.0 * np.pi * var)
    return pref * (np.exp(-((x - Rp) ** 2) / (2 * var)) + np.exp(-((x + Rp) ** 2) / (2 * var)))


def surface_source_drive_in(x_um, dose_cm2, Dt_cm2):
    """표면의 dose Q가 drive-in으로 퍼진 반쪽 Gaussian [cm^-3] (Plummer Ch. 8, Eq. 8)."""
    x = np.asarray(x_um, dtype=float) * UM_TO_CM
    return dose_cm2 / np.sqrt(np.pi * Dt_cm2) * np.exp(-(x**2) / (4.0 * Dt_cm2))


def predeposition_erfc(x_um, Cs, Dt_cm2):
    """표면 농도 Cs가 일정한 predeposition 분포 [cm^-3]."""
    x = np.asarray(x_um, dtype=float) * UM_TO_CM
    return Cs * erfc(x / (2.0 * np.sqrt(Dt_cm2)))


def predeposition_dose(Cs, Dt_cm2):
    """predeposition으로 들어간 dose Q = 2 Cs sqrt(Dt/pi) [cm^-2]."""
    return 2.0 * Cs * np.sqrt(Dt_cm2 / np.pi)


def surface_fraction_lost(Rp_um, dRp_um):
    """주입 Gaussian 중 표면 위(x < 0)에 해당하는 비율. 영상 근사의 오차 지표다."""
    return 0.5 * (1.0 - erf(Rp_um / (np.sqrt(2.0) * dRp_um)))


# ---------------------------------------------------------------------------
# 수치해 (Crank-Nicolson)
# ---------------------------------------------------------------------------
def diffuse_fd(x_um, C0, Dt_cm2, n_steps=400, surface="reflecting", Cs=None):
    """일정한 D에 대한 1D 확산을 Crank-Nicolson으로 푼다.

    D와 t를 따로 줄 필요 없이 Dt만 있으면 된다(무차원 시간 tau = Dt).

    Parameters
    ----------
    x_um : 균일 격자 [um], x[0] = 0이 표면
    C0 : 초기 농도 [cm^-3]
    surface : "reflecting"(외부 확산 없음) 또는 "fixed"(표면 농도 Cs 고정)
    아래쪽 경계는 반사 경계다. 도메인을 확산 길이보다 충분히 깊게 잡아야 한다.
    """
    x = np.asarray(x_um, dtype=float) * UM_TO_CM
    dx = x[1] - x[0]
    if not np.allclose(np.diff(x), dx):
        raise ValueError("x는 균일 격자여야 합니다.")
    n = x.size
    r = (Dt_cm2 / n_steps) / dx**2

    # 반사 경계는 ghost node 대칭(C[-1] = C[1])으로 처리한다.
    L = np.zeros((3, n))                      # (I - r/2 * Lap), banded
    L[0, 1:] = -r / 2
    L[1, :] = 1 + r
    L[2, :-1] = -r / 2
    L[0, 1] = -r                              # 표면 행: C[1] 계수가 두 배
    L[2, n - 2] = -r                          # 바닥 행
    if surface == "fixed":
        if Cs is None:
            raise ValueError("surface='fixed'에는 Cs가 필요합니다.")
        L[1, 0] = 1.0
        L[0, 1] = 0.0
    elif surface != "reflecting":
        raise ValueError("surface는 'reflecting' 또는 'fixed'입니다.")

    C = np.asarray(C0, dtype=float).copy()
    if surface == "fixed":
        C[0] = Cs
    for _ in range(n_steps):
        lap = np.empty_like(C)
        lap[1:-1] = C[:-2] - 2 * C[1:-1] + C[2:]
        lap[0] = 2 * (C[1] - C[0])
        lap[-1] = 2 * (C[-2] - C[-1])
        rhs = C + (r / 2) * lap
        if surface == "fixed":
            rhs[0] = Cs
        C = solve_banded((1, 1), L, rhs)
    return C


# ---------------------------------------------------------------------------
# 분포 해석
# ---------------------------------------------------------------------------
def dose_of(x_um, C):
    """분포를 적분한 dose [cm^-2]."""
    return float(np.trapezoid(C, np.asarray(x_um) * UM_TO_CM))


def junction_depth(x_um, C, N_background):
    """도핑 분포가 기판 농도와 같아지는 가장 깊은 지점 [um]. 없으면 nan.

    교차 구간은 로그 농도로 선형 보간한다.
    """
    x = np.asarray(x_um, dtype=float)
    C = np.asarray(C, dtype=float)
    above = C > N_background
    if not above.any() or above[-1]:
        return float("nan")
    i = np.nonzero(above)[0][-1]              # 마지막으로 기판 농도보다 높은 점
    y0, y1 = np.log(C[i]), np.log(max(C[i + 1], 1e-300))
    f = (np.log(N_background) - y0) / (y1 - y0)
    return float(x[i] + f * (x[i + 1] - x[i]))


def gaussian_drive_in_Dt_for_xj(xj_um, Cs, N_background):
    """표면 Gaussian 분포에서 접합 깊이 xj를 만드는 Dt [cm^2] (Plummer Ch. 7 예제의 식)."""
    xj = xj_um * UM_TO_CM
    return xj**2 / (4.0 * np.log(Cs / N_background))
