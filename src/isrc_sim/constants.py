"""물리 상수와 교과서 공정 파라미터.

모든 값에는 출처를 함께 적는다. 단위는 변수 이름 끝에 붙인다.
"""

# Boltzmann 상수 [eV/K]
K_B_EV = 8.617333262e-5


def celsius_to_kelvin(T_C):
    return T_C + 273.15


# ---------------------------------------------------------------------------
# Si 소모 비율
# 산화막 두께 D를 키울 때 소모되는 Si 두께 = ratio * D
# 밀도와 몰질량: CRC Handbook of Chemistry and Physics, 77th ed. (1996)
# ---------------------------------------------------------------------------
RHO_SI_G_CM3 = 2.329      # 결정질 Si
RHO_SIO2_G_CM3 = 2.196    # 비정질 SiO2
M_SI_G_MOL = 28.0855
M_O_G_MOL = 15.9994
M_SIO2_G_MOL = M_SI_G_MOL + 2 * M_O_G_MOL

SI_CONSUMPTION_RATIO = (M_SI_G_MOL / RHO_SI_G_CM3) / (M_SIO2_G_MOL / RHO_SIO2_G_CM3)  # 약 0.4407


# ---------------------------------------------------------------------------
# Deal-Grove 계수 (1 atm 기준)
#   B   = B0  * exp(-E_B  / kT)   [um^2/h]  포물선 상수, 결정 방향과 무관
#   B/A = BA0 * exp(-E_BA / kT)   [um/h]    선형 상수, 결정 방향에 따라 다름
# 출처: B. E. Deal and A. S. Grove, J. Appl. Phys. 36, 3770 (1965)
#       R. C. Jaeger, Introduction to Microelectronic Fabrication, 2nd ed. (2002)
#       J. D. Plummer et al., Silicon VLSI Technology (2000)
# 교재 표는 (111) 기준이며, (100)은 BA0를 1.68로 나눈다.
# wet은 H2O(steam) 계수다. 이 계수의 기준 분압은 0.92 atm(pyrogenic)으로 둔다.
# nanoFAB 계산기 출력과 비교해 이 가정이 맞음을 확인했다(docs/01_oxidation.md).
# ---------------------------------------------------------------------------
ORIENTATION_BA_FACTOR_111_OVER_100 = 1.68

DEAL_GROVE = {
    "dry": {
        "B0_um2_h": 7.72e2,
        "E_B_eV": 1.23,
        "BA0_um_h": {"111": 6.23e6, "100": 6.23e6 / ORIENTATION_BA_FACTOR_111_OVER_100},
        "E_BA_eV": 2.00,
        "p_ref_atm": 1.00,
    },
    "wet": {
        "B0_um2_h": 3.86e2,
        "E_B_eV": 0.78,
        "BA0_um_h": {"111": 1.63e8, "100": 1.63e8 / ORIENTATION_BA_FACTOR_111_OVER_100},
        "E_BA_eV": 2.05,
        "p_ref_atm": 0.92,
    },
}

# 두꺼운 dry 산화막을 Deal-Grove로 맞출 때 쓰는 관례적 초기 두께 [um]
# 실제 초기 성장은 더 빠르므로, 처음부터 약 25 nm가 있었다고 가정하면 잘 맞는다.
DRY_INITIAL_OFFSET_UM = 0.025


# ---------------------------------------------------------------------------
# Massoud 얇은 산화막 모델, dry O2, (100), 800~1000 C, 1 atm
#   x^2 + A x = B t + M1 (1 - exp(-t/tau1)) + M2 (1 - exp(-t/tau2)) + M0
#   M1 = K1 tau1, M2 = K2 tau2, M0 = x_i^2 + A x_i
#   B, B/A, K1, K2 = C0 * exp(-E/kT)
#   tau1, tau2     = C0 * exp(+E/kT)   (온도가 높을수록 짧아진다)
# 출처: H. Z. Massoud, J. D. Plummer, and E. A. Irene,
#       J. Electrochem. Soc. 132, 1745 (1985)
#       값은 이 논문을 인용한 L. Filipovic, PhD thesis, TU Wien (2012), Table 2.3, 2.4에서 옮겼다.
# 단위: nm, min
# tau의 지수 부호를 exp(-E/kT)로 읽으면 tau가 1e-12 min 수준이 되어 비물리적이다.
# exp(+E/kT)로 계산한 값이 nanoFAB 계산기와 일치함을 확인했다.
# ---------------------------------------------------------------------------
MASSOUD_DRY_100 = {
    "B": (1.70e11, 2.22),       # nm^2/min, eV
    "BA": (7.35e6, 1.76),       # nm/min, eV
    "K1": (2.49e11, 2.18),      # nm^2/min, eV
    "K2": (3.72e11, 2.28),      # nm^2/min, eV
    "tau1": (4.14e-6, 1.38),    # min, eV
    "tau2": (2.71e-7, 1.88),    # min, eV
    "valid_T_C": (800.0, 1000.0),
}
