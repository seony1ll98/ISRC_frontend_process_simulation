"""열산화(Thermal Oxidation) 모델.

- Deal-Grove 모델: dry O2와 wet(H2O), (100)과 (111)
- Massoud 모델: 얇은 dry 산화막, (100), 800~1000 C

공통 단위: 두께 [um], 시간 [h], 온도 [C]
"""

import warnings

import numpy as np
from scipy.optimize import brentq

from .constants import (
    DEAL_GROVE,
    K_B_EV,
    MASSOUD_DRY_100,
    SI_CONSUMPTION_RATIO,
    celsius_to_kelvin,
)


def _arrhenius(prefactor, E_eV, T_C, sign=-1.0):
    return prefactor * np.exp(sign * E_eV / (K_B_EV * celsius_to_kelvin(T_C)))


def _check_inputs(ambient, orientation):
    if ambient not in DEAL_GROVE:
        raise ValueError(f"ambient는 'dry' 또는 'wet'이어야 합니다: {ambient!r}")
    if orientation not in ("100", "111"):
        raise ValueError(f"orientation은 '100' 또는 '111'이어야 합니다: {orientation!r}")


# ---------------------------------------------------------------------------
# Deal-Grove
# ---------------------------------------------------------------------------
def deal_grove_rates(T_C, ambient="dry", orientation="100", pressure_atm=None):
    """Deal-Grove 속도 상수를 반환한다.

    Returns
    -------
    B : 포물선 상수 [um^2/h]
    B_over_A : 선형 상수 [um/h]
    A : B / (B/A) [um]

    pressure_atm을 주면 B와 B/A를 기준 분압 대비 선형으로 보정한다(Henry 법칙).
    """
    _check_inputs(ambient, orientation)
    p = DEAL_GROVE[ambient]
    B = _arrhenius(p["B0_um2_h"], p["E_B_eV"], T_C)
    B_over_A = _arrhenius(p["BA0_um_h"][orientation], p["E_BA_eV"], T_C)
    if pressure_atm is not None:
        scale = pressure_atm / p["p_ref_atm"]
        B = B * scale
        B_over_A = B_over_A * scale
    return B, B_over_A, B / B_over_A


def deal_grove_thickness(t_h, T_C, ambient="dry", orientation="100", x_i_um=0.0, pressure_atm=None):
    """산화 시간 t_h [h] 후의 전체 산화막 두께 [um] (초기 산화막 포함)."""
    B, _, A = deal_grove_rates(T_C, ambient, orientation, pressure_atm)
    tau = (x_i_um**2 + A * x_i_um) / B
    t = np.asarray(t_h, dtype=float)
    return (-A + np.sqrt(A**2 + 4.0 * B * (t + tau))) / 2.0


def deal_grove_time(x_um, T_C, ambient="dry", orientation="100", x_i_um=0.0, pressure_atm=None):
    """목표 두께 x_um [um]까지 필요한 산화 시간 [h]. 레시피 설계용 역계산이다."""
    B, _, A = deal_grove_rates(T_C, ambient, orientation, pressure_atm)
    x = np.asarray(x_um, dtype=float)
    if np.any(x < x_i_um):
        raise ValueError("목표 두께는 초기 두께보다 커야 합니다.")
    tau = (x_i_um**2 + A * x_i_um) / B
    return (x**2 + A * x) / B - tau


def deal_grove_growth_rate(x_um, T_C, ambient="dry", orientation="100", pressure_atm=None):
    """두께 x에서의 성장 속도 dx/dt = B / (A + 2x) [um/h]."""
    B, _, A = deal_grove_rates(T_C, ambient, orientation, pressure_atm)
    return B / (A + 2.0 * np.asarray(x_um, dtype=float))


def regime_crossover_thickness(T_C, ambient="dry", orientation="100"):
    """계면 반응 저항과 산화막 확산 저항이 같아지는 두께 x = A/2 [um].

    이보다 얇으면 반응 율속(선형), 두꺼우면 확산 율속(포물선) 영역이다.
    """
    _, _, A = deal_grove_rates(T_C, ambient, orientation)
    return A / 2.0


# ---------------------------------------------------------------------------
# Massoud (얇은 dry 산화막)
# ---------------------------------------------------------------------------
def massoud_parameters(T_C):
    """Massoud 모델 파라미터를 nm, min 단위로 반환한다."""
    lo, hi = MASSOUD_DRY_100["valid_T_C"]
    if np.any(np.asarray(T_C) < lo) or np.any(np.asarray(T_C) > hi):
        warnings.warn(f"Massoud 계수는 {lo:.0f}~{hi:.0f} C 범위에서 구한 값입니다.", stacklevel=3)
    m = MASSOUD_DRY_100
    B = _arrhenius(*m["B"], T_C)
    B_over_A = _arrhenius(*m["BA"], T_C)
    K1 = _arrhenius(*m["K1"], T_C)
    K2 = _arrhenius(*m["K2"], T_C)
    tau1 = _arrhenius(*m["tau1"], T_C, sign=+1.0)
    tau2 = _arrhenius(*m["tau2"], T_C, sign=+1.0)
    return {"B": B, "B_over_A": B_over_A, "A": B / B_over_A,
            "K1": K1, "K2": K2, "tau1": tau1, "tau2": tau2}


def massoud_thickness(t_h, T_C, x_i_um=0.0):
    """Massoud 모델로 계산한 dry (100) 산화막 두께 [um]."""
    p = massoud_parameters(T_C)
    t_min = np.asarray(t_h, dtype=float) * 60.0
    x_i_nm = x_i_um * 1e3
    A = p["A"]
    rhs = (p["B"] * t_min
           + p["K1"] * p["tau1"] * (1.0 - np.exp(-t_min / p["tau1"]))
           + p["K2"] * p["tau2"] * (1.0 - np.exp(-t_min / p["tau2"]))
           + x_i_nm**2 + A * x_i_nm)
    x_nm = (-A + np.sqrt(A**2 + 4.0 * rhs)) / 2.0
    return x_nm * 1e-3


def massoud_time(x_um, T_C, x_i_um=0.0):
    """Massoud 모델에서 목표 두께까지 필요한 시간 [h]. 수치 해법(brentq)을 쓴다."""
    if x_um <= x_i_um:
        raise ValueError("목표 두께는 초기 두께보다 커야 합니다.")
    f = lambda t: massoud_thickness(t, T_C, x_i_um) - x_um
    t_hi = 1.0
    while f(t_hi) < 0:
        t_hi *= 2.0
    return brentq(f, 0.0, t_hi, xtol=1e-10)


# ---------------------------------------------------------------------------
# 공통
# ---------------------------------------------------------------------------
def si_consumed(x_ox_um):
    """산화막 두께 x_ox_um [um]를 만들 때 소모되는 Si 두께 [um]."""
    return SI_CONSUMPTION_RATIO * np.asarray(x_ox_um, dtype=float)
