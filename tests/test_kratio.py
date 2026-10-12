import pytest

from epma_thins.data.lines import XrayLine
from epma_thins.kratio import bulk_kratio

AL_TI = {"Al": 0.398, "Ti": 0.602}

# [Moy21a]p14 Tables 1-2: Al0.398Ti0.602 (wt), TOA 40 deg, MAC30.
# (kV, experimental k [PAP91] Appendix 6, CalcZAF PAP k)
TI_KA = [
    (7.0, None, 0.561), (10.5, None, 0.564), (15.7, 0.579, 0.562), (20.9, 0.561, 0.557),
    (26.1, 0.558, 0.550), (31.3, 0.544, 0.542), (36.5, 0.544, 0.533), (40.0, 0.538, 0.526),
]
AL_KA = [
    (7.0, 0.385, 0.396), (10.5, 0.354, 0.362), (15.7, 0.301, 0.310), (20.8, 0.259, 0.263),
    (26.0, 0.216, 0.223), (31.3, 0.187, 0.190), (36.4, 0.158, 0.166), (40.0, 0.143, 0.152),
]


@pytest.mark.parametrize("line, kv, k_exp, k_calczaf",
                         [("Ti", *row) for row in TI_KA] + [("Al", *row) for row in AL_KA])
def test_al_ti_matches_calczaf_pap(line, kv, k_exp, k_calczaf):
    k = bulk_kratio(AL_TI, XrayLine(line, "Ka1"), kv)
    # Published to 3 decimals; PAP + MAC30 in CalcZAF and GMRFilm agree with this
    assert k == pytest.approx(k_calczaf, abs=1e-3)
    if k_exp is not None:
        # PAP itself is within a few percent of experiment (Al high-kV is ~6% high,
        # as for every program in [Moy21a])
        assert k == pytest.approx(k_exp, rel=0.07)


def test_pure_element_is_one():
    assert bulk_kratio({"Cu": 1.0}, XrayLine("Cu", "Ka1"), 20.0) == pytest.approx(1.0)


def test_compound_standard():
    # k against a compound standard equals the ratio of k's against the pure element
    line = XrayLine("Al", "Ka1")
    std = {"Al": 0.5293, "O": 0.4707}
    k_vs_pure = bulk_kratio(AL_TI, line, 15.0)
    k_std_vs_pure = bulk_kratio(std, line, 15.0)
    assert bulk_kratio(AL_TI, line, 15.0, standard=std) == pytest.approx(k_vs_pure / k_std_vs_pure)
