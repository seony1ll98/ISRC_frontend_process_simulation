import csv
from pathlib import Path

import numpy as np
import pytest

from isrc_sim.constants import SI_CONSUMPTION_RATIO
from isrc_sim.oxidation import (
    deal_grove_rates,
    deal_grove_thickness,
    deal_grove_time,
    massoud_thickness,
    massoud_time,
    regime_crossover_thickness,
    si_consumed,
)

REF_CSV = Path(__file__).resolve().parents[1] / "data" / "reference" / "nanofab_sithox_2026-10-11.csv"


# ---------------------------------------------------------------------------
# Deal-Grove: 물리적 극한
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("ambient", ["dry", "wet"])
def test_short_time_is_linear(ambient):
    B, B_over_A, _ = deal_grove_rates(1000, ambient)
    t = 1e-4  # h
    x = deal_grove_thickness(t, 1000, ambient)
    assert x == pytest.approx(B_over_A * t, rel=1e-2)


@pytest.mark.parametrize("ambient", ["dry", "wet"])
def test_long_time_is_parabolic(ambient):
    B, _, _ = deal_grove_rates(1000, ambient)
    t = 1e5  # h
    x = deal_grove_thickness(t, 1000, ambient)
    assert x == pytest.approx(np.sqrt(B * t), rel=1e-2)


def test_thickness_increases_with_time_and_temperature():
    t = np.linspace(0, 10, 200)
    for T in (900, 1000, 1100):
        assert np.all(np.diff(deal_grove_thickness(t, T, "dry")) > 0)
    xs = [deal_grove_thickness(1.0, T, "dry") for T in (900, 1000, 1100, 1200)]
    assert np.all(np.diff(xs) > 0)


def test_wet_faster_than_dry_and_111_faster_than_100():
    assert deal_grove_thickness(1.0, 1000, "wet") > deal_grove_thickness(1.0, 1000, "dry")
    assert deal_grove_thickness(1.0, 1000, "dry", "111") > deal_grove_thickness(1.0, 1000, "dry", "100")


def test_orientation_changes_only_linear_constant():
    B100, BA100, _ = deal_grove_rates(1000, "dry", "100")
    B111, BA111, _ = deal_grove_rates(1000, "dry", "111")
    assert B100 == pytest.approx(B111)
    assert BA111 / BA100 == pytest.approx(1.68)


def test_inverse_roundtrip():
    for ambient, x_i in (("dry", 0.0), ("dry", 0.025), ("wet", 0.0)):
        t = np.array([0.1, 1.0, 5.0])
        x = deal_grove_thickness(t, 1050, ambient, x_i_um=x_i)
        assert deal_grove_time(x, 1050, ambient, x_i_um=x_i) == pytest.approx(t, rel=1e-9)


def test_initial_oxide_is_returned_at_time_zero():
    assert deal_grove_thickness(0.0, 1000, "dry", x_i_um=0.025) == pytest.approx(0.025)


def test_crossover_rates_are_equal():
    # x = A/2에서 확산 저항(2x/B)과 반응 저항(A/B)이 같다.
    B, _, A = deal_grove_rates(1000, "dry")
    x = regime_crossover_thickness(1000, "dry")
    assert 2 * x / B == pytest.approx(A / B)


# ---------------------------------------------------------------------------
# Massoud
# ---------------------------------------------------------------------------
def test_massoud_thicker_than_deal_grove_for_thin_dry_oxide():
    for T in (800, 900, 1000):
        assert massoud_thickness(0.5, T) > deal_grove_thickness(0.5, T, "dry")


def test_massoud_inverse_roundtrip():
    x = massoud_thickness(1.3, 900)
    assert massoud_time(x, 900) == pytest.approx(1.3, rel=1e-6)


def test_massoud_warns_outside_valid_range():
    with pytest.warns(UserWarning):
        massoud_thickness(1.0, 1100)


# ---------------------------------------------------------------------------
# Si 소모
# ---------------------------------------------------------------------------
def test_si_consumption_ratio():
    assert SI_CONSUMPTION_RATIO == pytest.approx(0.4407, abs=1e-4)
    assert si_consumed(1.0) == pytest.approx(0.4407, abs=1e-4)


# ---------------------------------------------------------------------------
# 독립 구현(nanoFAB 계산기)과의 비교
# ---------------------------------------------------------------------------
def _load_reference():
    with open(REF_CSV, newline="") as f:
        return list(csv.DictReader(f))


@pytest.mark.parametrize("row", _load_reference(), ids=lambda r: f"{r['model']}-{r['ambient']}-{r['T_C']}C-{r['time_min']}min")
def test_matches_nanofab_calculator(row):
    t_h = float(row["time_min"]) / 60.0
    T = float(row["T_C"])
    x_i = float(row["xi_A"]) * 1e-4
    if row["model"] == "Massoud":
        x = massoud_thickness(t_h, T, x_i_um=x_i)
    else:
        x = deal_grove_thickness(t_h, T, row["ambient"], row["orientation"],
                                 x_i_um=x_i, pressure_atm=float(row["pp_atm"]))
    assert x * 1e4 == pytest.approx(float(row["thickness_A"]), rel=5e-3)
