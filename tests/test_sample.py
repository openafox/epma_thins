import pytest

from epma_thins import Layer, Measurement, Sample, XrayLine
from epma_thins.sample import mass_thickness_to_nm, nm_to_mass_thickness


def test_formula_weight_fractions():
    # Al2O3 = 52.93 wt% Al, 47.07 wt% O [Moy21b]p4
    layer = Layer.from_formula("Al2O3")
    assert layer.composition["Al"] == pytest.approx(0.5293, abs=2e-4)
    assert layer.composition["O"] == pytest.approx(0.4707, abs=2e-4)
    assert layer.total == pytest.approx(1.0)


def test_thickness_conversion():
    # 100 nm Al2O3 at 3.95 g/cm3 = 39.5 ug/cm2
    assert nm_to_mass_thickness(100, 3.95) == pytest.approx(39.5e-6)
    assert mass_thickness_to_nm(39.5e-6, 3.95) == pytest.approx(100)
    layer = Layer.from_nm("Al2O3", 100, 3.95)
    assert layer.mass_thickness == pytest.approx(39.5e-6)
    assert layer.thickness_nm == pytest.approx(100)


def test_composition_symbols_normalized():
    layer = Layer({"cu": 0.5, "AL": 0.25, "Cu": 0.25})
    assert layer.composition == {"Cu": 0.75, "Al": 0.25}
    with pytest.raises(ValueError):
        Layer({"Cu": -0.1})


def test_sample_structure():
    film = Layer.from_nm("Al2O3", 100, 3.95)
    sub = Layer({"Cu": 1.0})
    sample = Sample([film, sub])
    assert sample.substrate is sub
    assert sample.films == [film]
    assert not sample.is_bulk
    assert sample.elements == {"Al", "O", "Cu"}
    assert sample.boundaries() == pytest.approx([0.0, 39.5e-6])
    with pytest.raises(ValueError):
        Sample([Layer({"Al": 1.0}), sub])  # film without thickness
    assert Sample.bulk("SiO2").is_bulk


def test_measurement_defaults_and_checks():
    m = Measurement(XrayLine("Al", "Ka1"), beam_energy=15.0, kratio=0.07)
    assert m.standard.layers[0].composition == {"Al": 1.0}
    with pytest.raises(ValueError):
        Measurement(XrayLine("Cu", "Ka1"), beam_energy=5.0, kratio=0.1)  # below Cu K edge
    with pytest.raises(ValueError):
        Measurement(XrayLine("Al", "Ka1"), 15.0, 0.1, standard=Sample.bulk({"Si": 1.0}))
