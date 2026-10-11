import csv
from pathlib import Path

import numpy as np
import pytest

from isrc_sim.constants import IMPLANT_ENERGY_KEV, IMPLANT_RANGE_UM
from isrc_sim.diffusion import (
    annealed_gaussian,
    diffuse_fd,
    dose_of,
    gaussian_drive_in_Dt_for_xj,
    intrinsic_diffusivity,
    junction_depth,
    predeposition_dose,
    predeposition_erfc,
    surface_fraction_lost,
    surface_source_drive_in,
    thermal_budget,
)
from isrc_sim.electrical import mobility, resistivity, sheet_resistance
from isrc_sim.implant import implant_profile, peak_concentration, range_parameters

BYU_CSV = Path(__file__).resolve().parents[1] / "data" / "reference" / "byu_resistivity_2026-10-11.csv"


# ---------------------------------------------------------------------------
# 이온주입
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("dopant", ["B", "P", "As", "Sb"])
def test_range_interpolation_hits_table(dopant):
    Rp, dRp = range_parameters(dopant, np.array(IMPLANT_ENERGY_KEV, dtype=float))
    table = np.array(IMPLANT_RANGE_UM[dopant])
    assert Rp == pytest.approx(table[:, 0])
    assert dRp == pytest.approx(table[:, 1])


def test_range_increases_with_energy_and_lighter_ions_go_deeper():
    E = np.linspace(10, 200, 50)
    for d in ("B", "P", "As"):
        assert np.all(np.diff(range_parameters(d, E)[0]) > 0)
    assert range_parameters("B", 100)[0] > range_parameters("P", 100)[0] > range_parameters("As", 100)[0]


def test_range_outside_table_raises():
    with pytest.raises(ValueError):
        range_parameters("B", 5)


def test_implanted_dose_inside_silicon():
    # Gaussian 꼬리 중 표면 위(x < 0) 부분은 Si 안에 없으므로 그만큼 dose가 줄어든다.
    x = np.linspace(0, 2, 20001)
    C = implant_profile(x, "P", 100, 1e15)
    Rp, dRp = range_parameters("P", 100)
    assert dose_of(x, C) == pytest.approx(1e15 * (1 - surface_fraction_lost(Rp, dRp)), rel=1e-4)


def test_peak_concentration_formula():
    # Cp = 0.4 Q / dRp
    assert peak_concentration(1e15, 0.01) == pytest.approx(0.3989 * 1e15 / 1e-6, rel=1e-3)


# ---------------------------------------------------------------------------
# 확산계수
# ---------------------------------------------------------------------------
def test_boron_diffusivity_matches_plummer_example():
    # Plummer Ch. 7 예제: 1100 C에서 B의 D는 약 1.5e-13 cm^2/s
    assert intrinsic_diffusivity("B", 1100) == pytest.approx(1.5e-13, rel=0.1)


def test_diffusivity_ordering_at_1000C():
    # P와 B는 빠른 확산자, As와 Sb는 느린 확산자
    D = {d: intrinsic_diffusivity(d, 1000) for d in ("B", "P", "As", "Sb")}
    assert min(D["B"], D["P"]) > max(D["As"], D["Sb"])


def test_thermal_budget_adds_linearly():
    a = thermal_budget("B", [(1000, 1.0)])
    b = thermal_budget("B", [(1100, 0.5)])
    assert thermal_budget("B", [(1000, 1.0), (1100, 0.5)]) == pytest.approx(a + b)


# ---------------------------------------------------------------------------
# 해석해와 수치해 비교
# ---------------------------------------------------------------------------
def test_annealed_gaussian_conserves_dose():
    x = np.linspace(0, 6, 30001)
    for Dt in (0.0, 1e-10, 1e-9, 1e-8):
        C = annealed_gaussian(x, 1e13, 0.3, 0.086, Dt)
        assert dose_of(x, C) == pytest.approx(1e13, rel=1e-3)


def test_fd_matches_annealed_gaussian():
    x = np.linspace(0, 4, 2001)
    Rp, dRp = range_parameters("B", 100)
    C0 = annealed_gaussian(x, 1e13, Rp, dRp, 0.0)
    Dt = thermal_budget("B", [(1100, 2.0)])
    C_fd = diffuse_fd(x, C0, Dt, n_steps=400)
    C_an = annealed_gaussian(x, 1e13, Rp, dRp, Dt)
    mask = C_an > 1e-3 * C_an.max()
    assert np.max(np.abs(C_fd[mask] - C_an[mask]) / C_an[mask]) < 1e-2
    assert dose_of(x, C_fd) == pytest.approx(1e13, rel=1e-3)


def test_fd_matches_erfc_predeposition():
    x = np.linspace(0, 2, 2001)
    Dt = 1e-9
    C_fd = diffuse_fd(x, np.zeros_like(x), Dt, n_steps=2000, surface="fixed", Cs=1e20)
    C_an = predeposition_erfc(x, 1e20, Dt)
    mask = C_an > 1e-3 * 1e20
    assert np.max(np.abs(C_fd[mask] - C_an[mask]) / C_an[mask]) < 2e-2
    assert dose_of(x, C_fd) == pytest.approx(predeposition_dose(1e20, Dt), rel=1e-2)


def test_surface_source_drive_in_conserves_dose():
    x = np.linspace(0, 10, 50001)
    C = surface_source_drive_in(x, 1e13, 1e-8)
    assert dose_of(x, C) == pytest.approx(1e13, rel=1e-3)


# ---------------------------------------------------------------------------
# 접합 깊이
# ---------------------------------------------------------------------------
def test_junction_depth_matches_closed_form():
    Cs, NB, xj = 4e17, 1e15, 3.0
    Dt = gaussian_drive_in_Dt_for_xj(xj, Cs, NB)
    x = np.linspace(0, 6, 60001)
    C = Cs * np.exp(-((x * 1e-4) ** 2) / (4 * Dt))
    assert junction_depth(x, C, NB) == pytest.approx(xj, rel=1e-4)


def test_plummer_example_dt():
    # Plummer Ch. 7 예제: Cs = 4e17, NB = 1e15, xj = 3 um -> Dt = 3.7e-9 cm^2
    assert gaussian_drive_in_Dt_for_xj(3.0, 4e17, 1e15) == pytest.approx(3.7e-9, rel=0.02)


def test_junction_depth_nan_when_no_junction():
    x = np.linspace(0, 1, 101)
    assert np.isnan(junction_depth(x, np.full_like(x, 1e14), 1e15))


# ---------------------------------------------------------------------------
# 이동도와 면저항
# ---------------------------------------------------------------------------
def _load_byu():
    with open(BYU_CSV, newline="") as f:
        return list(csv.DictReader(f))


@pytest.mark.parametrize("row", [r for r in _load_byu() if r["dopant"] in ("As", "P")],
                         ids=lambda r: f"{r['dopant']}-{r['N_cm3']}")
def test_mobility_matches_byu_for_electrons(row):
    N = float(row["N_cm3"])
    assert mobility(N, row["dopant"]) == pytest.approx(float(row["mobility_cm2_Vs"]), rel=1e-4)
    assert resistivity(N, row["dopant"]) == pytest.approx(float(row["resistivity_ohm_cm"]), rel=1e-4)


def test_hole_mobility_approaches_mu_max_at_low_doping():
    # 원 논문(mu_min2 = 0)에서는 저농도 정공 이동도가 mu_max = 470.5에 접근해야 한다.
    assert mobility(1e12, "B") == pytest.approx(470.5, rel=1e-3)


def test_mobility_decreases_with_doping():
    N = np.logspace(14, 20, 50)
    for d in ("B", "P", "As"):
        assert np.all(np.diff(mobility(N, d)) < 0)


def test_sheet_resistance_of_uniform_layer():
    # 균일 도핑층: Rs = rho / t
    x = np.linspace(0, 1.0, 1001)
    C = np.full_like(x, 1e17)
    NB = 1e10
    expected = resistivity(1e17, "P") / 1e-4
    assert sheet_resistance(x, C, "P", NB) == pytest.approx(expected, rel=1e-3)
