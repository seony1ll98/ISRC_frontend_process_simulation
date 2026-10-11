"""도핑 분포의 전기적 특성: 이동도, 비저항, 면저항."""

import numpy as np

from .constants import MASETTI, Q_E

UM_TO_CM = 1e-4


def mobility(N, dopant):
    """Masetti 모델 이동도 [cm^2/Vs] (300 K). N은 전체 이온화 불순물 농도 [cm^-3]."""
    if dopant not in MASETTI:
        raise ValueError(f"지원하는 도펀트는 {sorted(MASETTI)}입니다: {dopant!r}")
    p = MASETTI[dopant]
    N = np.asarray(N, dtype=float)
    return (p["mu_min1"] * np.exp(-p["Pc"] / N)
            + (p["mu_max"] - p["mu_min2"]) / (1.0 + (N / p["Cr"]) ** p["alpha"])
            - p["mu_1"] / (1.0 + (p["Cs"] / N) ** p["beta"]))


def resistivity(N, dopant):
    """균일 도핑 Si의 비저항 [ohm cm]."""
    return 1.0 / (Q_E * np.asarray(N, dtype=float) * mobility(N, dopant))


def sheet_resistance(x_um, C, dopant, N_background):
    """접합 위 도핑층의 면저항 [ohm/sq].

    Rs = 1 / (q * integral_0^xj (C - NB) mu(C + NB) dx)  (Plummer Ch. 7, Eq. 3)
    순 캐리어 농도는 C - NB이고, 이동도는 전체 불순물 농도 C + NB로 계산한다.
    """
    x = np.asarray(x_um, dtype=float) * UM_TO_CM
    C = np.asarray(C, dtype=float)
    net = C - N_background
    mask = net > 0
    integrand = np.where(mask, net * mobility(C + N_background, dopant), 0.0)
    return 1.0 / (Q_E * np.trapezoid(integrand, x))
