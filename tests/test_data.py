import pytest

from epma_thins.data.atomic import atomic_number, atomic_weight, edge_energy, symbol
from epma_thins.data.lines import XrayLine
from epma_thins.data.mac import OverriddenMac, XraylibMac, mixture_mac


def test_atomic_lookup():
    assert atomic_number("cu") == 29
    assert symbol(29) == "Cu"
    assert atomic_weight("Cu") == pytest.approx(63.55, abs=0.02)
    assert edge_energy("Cu", "K") == pytest.approx(8.979, abs=0.002)
    with pytest.raises(ValueError):
        atomic_number("Xx")


@pytest.mark.parametrize(
    "element, line, canonical, shell, energy",
    [
        ("Cu", "Ka1", "Ka1", "K", 8.048),
        ("Cu", "ka", "Ka", "K", 8.041),
        ("Si", "K-L3", "K-L3", "K", 1.740),
        ("Au", "La1", "La1", "L3", 9.713),
        ("Au", "Ma", "Ma", "M5", 2.123),
        ("Gd", "La", "La", "L3", 6.05),
    ],
)
def test_lines(element, line, canonical, shell, energy):
    xl = XrayLine(element, line)
    assert xl.line == canonical
    assert xl.shell == shell
    assert xl.energy == pytest.approx(energy, abs=0.01)
    assert xl.edge > xl.energy


def test_line_equality_and_hash():
    assert XrayLine("cu", "ka1") == XrayLine("Cu", "Ka1")
    assert len({XrayLine("Cu", "Ka1"), XrayLine("Cu", "Ka1")}) == 1


def test_invalid_lines():
    with pytest.raises(ValueError):
        XrayLine("Cu", "Zz")
    with pytest.raises(ValueError):
        XrayLine("C", "La1")  # carbon has no L lines


def test_xraylib_macs_reasonable():
    mac = XraylibMac()
    # Ti Ka in Ti: ~110 cm2/g (legacy mac.py returned ~4e-7)
    assert mac(XrayLine("Ti", "Ka1"), "Ti") == pytest.approx(110, rel=0.05)
    # Al Ka in Ti: PENELOPE-2018 gives 2132 cm2/g [Moy21a]p13
    assert mac(XrayLine("Al", "Ka1"), "Ti") == pytest.approx(2132, rel=0.02)
    assert XraylibMac("photo")(XrayLine("Al", "Ka1"), "Ti") < mac(XrayLine("Al", "Ka1"), "Ti")


def test_mac_override_and_mixture():
    al = XrayLine("Al", "Ka1")
    mac = OverriddenMac(XraylibMac(), {(al, "ti"): 2201.0})
    assert mac(al, "Ti") == 2201.0
    assert mac(al, "Al") == XraylibMac()(al, "Al")
    mix = mixture_mac(mac, al, {"Al": 0.398, "Ti": 0.602})
    assert mix == pytest.approx(0.398 * mac(al, "Al") + 0.602 * 2201.0)


@pytest.mark.parametrize(
    "element, line, family",
    [("Cu", "Ka1", "Ka"), ("Cu", "Ka", "Ka"), ("Cu", "K-L3", "Ka"), ("Fe", "La1", "La"),
     ("Au", "Ma", "Ma"), ("Au", "M4-N6", "Mb"), ("Cu", "Lb1", "Lb")],
)
def test_line_family(element, line, family):
    assert XrayLine(element, line).family == family


def test_override_by_family_and_exact_precedence():
    mac = OverriddenMac(XraylibMac(), {(("Fe", "La"), "Fe"): 3350.0,
                                       (XrayLine("Fe", "La2"), "Fe"): 1.0})
    assert mac(XrayLine("Fe", "La1"), "Fe") == 3350.0
    assert mac(XrayLine("Fe", "La2"), "Fe") == 1.0
    assert mac(XrayLine("Fe", "Ka1"), "Fe") == XraylibMac()(XrayLine("Fe", "Ka1"), "Fe")


def test_pap91_mac_set():
    from epma_thins.data.mac import pap91_mac
    from epma_thins.data.mac30 import Mac30

    mac = pap91_mac()
    # [PAP91]p63 Appendix 5 values
    assert mac(XrayLine("B", "Ka"), "Zr") == 4400  # MAC30 itself fails here
    assert mac(XrayLine("C", "Ka"), "C") == 2170
    assert mac(XrayLine("N", "Ka"), "Ti") == 4270
    assert mac(XrayLine("Fe", "La1"), "Fe") == 3350
    assert mac(XrayLine("Si", "Ka1"), "Ta") == 1490
    assert mac(XrayLine("Au", "Mb"), "Pt") == 2550
    # everything else is MAC30
    al_ka = XrayLine("Al", "Ka1")
    assert mac(al_ka, "Ti") == Mac30()(al_ka, "Ti")
