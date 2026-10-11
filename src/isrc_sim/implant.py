"""이온주입 프로파일 모델.

공통 단위: 깊이 [um], 에너지 [keV], dose [cm^-2], 농도 [cm^-3]
"""

import numpy as np
from scipy.interpolate import PchipInterpolator

from .constants import IMPLANT_ENERGY_KEV, IMPLANT_RANGE_UM

UM_TO_CM = 1e-4


def _check_dopant(dopant):
    if dopant not in IMPLANT_RANGE_UM:
        raise ValueError(f"지원하는 도펀트는 {sorted(IMPLANT_RANGE_UM)}입니다: {dopant!r}")


def range_parameters(dopant, energy_keV):
    """투영 비정 Rp와 표준편차 dRp [um]를 반환한다.

    Plummer 표(10~200 keV)를 단조 보간(PCHIP)한다. 표 범위 밖은 외삽하지 않는다.
    """
    _check_dopant(dopant)
    E = np.asarray(energy_keV, dtype=float)
    lo, hi = IMPLANT_ENERGY_KEV[0], IMPLANT_ENERGY_KEV[-1]
    if np.any(E < lo) or np.any(E > hi):
        raise ValueError(f"에너지는 {lo}~{hi} keV 범위여야 합니다.")
    table = np.array(IMPLANT_RANGE_UM[dopant])
    Rp = PchipInterpolator(IMPLANT_ENERGY_KEV, table[:, 0])(E)
    dRp = PchipInterpolator(IMPLANT_ENERGY_KEV, table[:, 1])(E)
    if E.ndim == 0:
        return float(Rp), float(dRp)
    return Rp, dRp


def peak_concentration(dose_cm2, dRp_um):
    """Gaussian 프로파일의 최대 농도 Cp = Q / (sqrt(2 pi) dRp) [cm^-3]."""
    return dose_cm2 / (np.sqrt(2.0 * np.pi) * dRp_um * UM_TO_CM)


def gaussian_profile(x_um, dose_cm2, Rp_um, dRp_um):
    """주입 직후 Gaussian 농도 분포 [cm^-3]."""
    x = np.asarray(x_um, dtype=float)
    return peak_concentration(dose_cm2, dRp_um) * np.exp(-((x - Rp_um) ** 2) / (2.0 * dRp_um**2))


def implant_profile(x_um, dopant, energy_keV, dose_cm2):
    """도펀트, 에너지, dose로 주입 직후 농도 분포를 계산한다."""
    Rp, dRp = range_parameters(dopant, energy_keV)
    return gaussian_profile(x_um, dose_cm2, Rp, dRp)
