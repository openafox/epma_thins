import math

import pytest

from epma_thins.data.lines import XrayLine
from epma_thins.data.mac import XraylibMac
from epma_thins.data.mac30 import Mac30, mac30_energy


def test_matches_stratagem_mac30_value():
    # STRATAGem (MAC30) uses 2201 cm2/g for Al Ka in Ti [Moy21a]p13
    assert Mac30()(XrayLine("Al", "Ka1"), "Ti") == pytest.approx(2201, rel=1e-3)


@pytest.mark.parametrize(
    "element, line, absorber",
    [
        ("Al", "Ka1", "Al"), ("Ti", "Ka1", "Al"), ("Ti", "Ka1", "Ti"), ("Si", "Ka1", "Si"),
        ("Si", "Ka1", "Fe"), ("Fe", "Ka1", "Fe"), ("Fe", "Ka1", "Si"), ("Cu", "Ka1", "Cu"),
        ("Cu", "Ka1", "Au"), ("Au", "Ma", "Cu"), ("Au", "La1", "Ag"), ("O", "Ka", "Al"),
        ("O", "Ka", "Si"), ("Ni", "Ka1", "Cr"), ("Cr", "Ka1", "Ni"), ("Gd", "La", "Fe"),
        ("Pt", "Ma", "Gd"), ("Zr", "La1", "B"),
    ],
)
def test_agrees_with_tabulated_macs_away_from_edges(element, line, absorber):
    # MAC30 is a fit to tabulated MACs; typical agreement is a few percent
    xl = XrayLine(element, line)
    assert Mac30()(xl, absorber) == pytest.approx(XraylibMac()(xl, absorber), rel=0.08)


@pytest.mark.parametrize("z", range(3, 96))
def test_finite_and_positive_across_periodic_table(z):
    for energy_ev in (277, 525, 1487, 4510, 8048, 14958):
        mu = mac30_energy(z, energy_ev)
        assert math.isfinite(mu) and mu > 0, (z, energy_ev)


def test_known_failure_raises():
    # MAC30 goes negative for B Ka in Zr; PAP replaced it with 4400 [PAP91]p63
    assert mac30_energy(40, 183.3) < 0
    with pytest.raises(ValueError, match="OverriddenMac"):
        Mac30()(XrayLine("B", "Ka"), "Zr")


def test_dtsa_variant_only_differs_for_heavy_absorbers_above_k_edge():
    cu_ka = 8048.0
    assert mac30_energy(29, cu_ka, dtsa=True) == mac30_energy(29, cu_ka)
    # Above the K edge of a Z >= 50 absorber the DTSA conditions apply
    assert mac30_energy(50, 40000.0, dtsa=True) != mac30_energy(50, 40000.0)
