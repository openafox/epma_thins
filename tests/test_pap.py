import math

import pytest
from scipy.integrate import quad

from epma_thins.data.lines import XrayLine
from epma_thins.models import pap

AL_KA = XrayLine("Al", "Ka1")

# [Moy21b]p4 Table 1: fictitious compositions (wt%) for 100 nm Al2O3 on Cu, Al Ka1, 15 kV
INITIAL = {"Al": 26.46, "O": 23.54, "Cu": 50.00}
ITERATION_1 = {"Al": 3.44, "O": 3.06, "Cu": 93.50}
ITERATION_2 = {"Al": 3.37, "O": 3.00, "Cu": 93.64}
BACKSCATTER = {"Al": 7.20, "O": 6.40, "Cu": 86.40}
DECELERATION = {"Al": 6.35, "O": 5.64, "Cu": 88.01}
U0 = 15.0 / AL_KA.edge


class TestMoy21bTable1:
    def test_ionization_range(self):
        assert pap.ionization_range(INITIAL, AL_KA, 15.0) * 1e6 == pytest.approx(686.7, abs=0.1)
        for comp in (ITERATION_1, ITERATION_2):
            assert pap.ionization_range(comp, AL_KA, 15.0) * 1e6 == pytest.approx(700.3, abs=0.1)

    def test_backscatter(self):
        assert pap.backscatter_z(BACKSCATTER) == pytest.approx(25.94, abs=0.01)
        assert pap.backscatter_factor(BACKSCATTER, U0) == pytest.approx(0.8471, abs=2e-4)
        assert pap.surface_ionization(BACKSCATTER, U0) == pytest.approx(1.6933, abs=2e-4)

    def test_deceleration_and_area(self):
        assert pap.mean_z(DECELERATION) == pytest.approx(26.8, abs=0.05)
        inv_s = pap.deceleration_factor(DECELERATION, AL_KA, 15.0)
        assert inv_s == pytest.approx(1.39e-4, rel=0.01)  # printed truncated; F needs 1.399e-4
        area = (pap.backscatter_factor(BACKSCATTER, U0) * inv_s
                / pap.ionization_cross_section(AL_KA, 15.0))
        assert area == pytest.approx(8.92e-4, abs=0.005e-4)

    def test_depth_of_maximum(self):
        z_av = (pap.backscatter_z(BACKSCATTER) + pap.mean_z(DECELERATION)) / 2
        assert z_av == pytest.approx(26.37, abs=0.01)
        assert pap.depth_of_maximum(z_av, U0, 700.3e-6) * 1e6 == pytest.approx(106.7, abs=0.1)


def test_mean_ionization_potential_pure_element():
    # Zeller's J_i [PAP91]p35e7 for Cu: 1e-3 * 29 * (10.04 + 8.25 exp(-29/11.22))
    expected = 1e-3 * 29 * (10.04 + 8.25 * math.exp(-29 / 11.22))
    assert pap.mean_ionization_potential({"Cu": 1.0}) == pytest.approx(expected)


def test_hutchins_exponent():
    assert pap.hutchins_m(XrayLine("Cu", "Ka1")) == pytest.approx(0.86, abs=1e-6)
    assert pap.hutchins_m(XrayLine("C", "Ka")) == pytest.approx(0.86 + 0.12 * math.exp(-1.44))
    assert pap.hutchins_m(XrayLine("Au", "La1")) == 0.82
    assert pap.hutchins_m(XrayLine("Au", "Ma")) == 0.78


@pytest.fixture(params=[("Cu", "Ka1", 15.0), ("Al", "Ka1", 5.0), ("Au", "Ma", 20.0),
                        ("Si", "Ka1", 30.0), ("Fe", "Ka1", 10.0)])
def bulk(request):
    el, line, e0 = request.param
    return pap.bulk_distribution({el: 1.0}, XrayLine(el, line), e0)


class TestDistributionShape:
    def test_area_equals_f(self, bulk):
        assert bulk.integral(0, bulk.rx) == pytest.approx(bulk.area, rel=1e-9)

    def test_boundary_conditions(self, bulk):
        assert bulk.phi(0.0) == pytest.approx(bulk.phi0, rel=1e-9)
        assert 0 < bulk.rm < bulk.rc < bulk.rx
        eps = 1e-6 * bulk.rx
        assert bulk.phi(bulk.rc - eps) == pytest.approx(bulk.phi(bulk.rc + eps), rel=1e-4)
        assert bulk.phi(bulk.rx - eps) == pytest.approx(0.0, abs=1e-9)
        assert bulk.phi(bulk.rm) >= bulk.phi(bulk.rm * 0.9)
        assert bulk.phi(bulk.rm) >= bulk.phi(bulk.rm * 1.1)

    def test_integral_is_additive(self, bulk):
        mid = 0.37 * bulk.rx
        chi = 1500.0
        whole = bulk.integral(0, bulk.rx, chi)
        parts = bulk.integral(0, mid, chi) + bulk.integral(mid, bulk.rx, chi)
        assert parts == pytest.approx(whole)

    @pytest.mark.parametrize("chi", [1e-3, 10.0, 300.0, 3000.0, 30000.0])
    def test_emitted_matches_quadrature(self, bulk, chi):
        expected = quad(lambda x: bulk.phi(x) * math.exp(-chi * x), 0, bulk.rx,
                        points=[bulk.rc], epsabs=0, epsrel=1e-12, limit=200)[0]
        assert bulk.emitted(chi) == pytest.approx(expected, rel=1e-9)
        assert bulk.emitted(chi) <= bulk.area

    @pytest.mark.parametrize("chi", [10.0, 300.0, 3000.0, 30000.0])
    def test_emitted_matches_corrected_closed_form(self, bulk, chi):
        # [PAP91]p38e21 as corrected in NIST EPQ PAP1991.java (sign of 2/chi^2 in F2)
        a1, a2, rm, rc, rx = bulk.a1, bulk.a2, bulk.rm, bulk.rc, bulk.rx
        f1 = (a1 / chi) * (((rc - rm) * (rx - rc - 2 / chi) - 2 / chi**2) * math.exp(-chi * rc)
                           - (rc - rm) * rx + rm * (rc - 2 / chi) + 2 / chi**2)
        f2 = (a2 / chi) * (((rx - rc) * (rx - rc - 2 / chi) + 2 / chi**2) * math.exp(-chi * rc)
                           - (2 / chi**2) * math.exp(-chi * rx))
        if chi * rx > 0.5:  # the closed form cancels badly at small chi
            assert bulk.emitted(chi) == pytest.approx(f1 + f2, rel=1e-8)

    def test_series_and_closed_form_agree_at_switchover(self, bulk):
        chi = 0.5 / bulk.rx
        below, above = bulk.emitted(chi * (1 - 1e-9)), bulk.emitted(chi * (1 + 1e-9))
        assert below == pytest.approx(above, rel=1e-7)


def test_below_edge_rejected():
    with pytest.raises(ValueError):
        pap.bulk_distribution({"Cu": 1.0}, XrayLine("Cu", "Ka1"), 8.0)


def test_low_overvoltage_fallback_keeps_valid_shape():
    # U0 = 1.25: the Rm(Rx) relation gives no root in (Rm, Rx); Appendix 3 applies
    d = pap.bulk_distribution({"Fe": 1.0}, XrayLine("Fe", "Ka1"), 8.9)
    assert 0 < d.rm < d.rc < d.rx
    assert d.integral(0, d.rx) == pytest.approx(d.area, rel=1e-9)


def test_very_low_overvoltage_not_implemented():
    # U0 = 1.12: φ(ρz) has no interior maximum; needs the exponential forms of Appendix 3
    with pytest.raises(NotImplementedError):
        pap.bulk_distribution({"Fe": 1.0}, XrayLine("Fe", "Ka1"), 8.0)
