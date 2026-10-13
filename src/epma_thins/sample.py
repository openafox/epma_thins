"""Sample description: layers, samples and k-ratio measurements.

Layers are ordered from the surface down; the last layer is the substrate
(semi-infinite, no thickness). A bulk sample is a Sample with one layer.
"""

from collections.abc import Mapping
from dataclasses import dataclass, field

import xraylib

from epma_thins.data.atomic import symbol
from epma_thins.data.lines import XrayLine


def _normalized(composition: Mapping[str, float]) -> dict[str, float]:
    comp = {}
    for el, w in composition.items():
        if w < 0:
            raise ValueError(f"negative weight fraction for {el}: {w}")
        key = symbol(el)
        comp[key] = comp.get(key, 0.0) + float(w)
    if not comp:
        raise ValueError("composition is empty")
    return comp


def nm_to_mass_thickness(thickness_nm: float, density: float) -> float:
    """Linear thickness [nm] and density [g/cm^3] -> mass thickness [g/cm^2]."""
    return thickness_nm * 1e-7 * density


def mass_thickness_to_nm(mass_thickness: float, density: float) -> float:
    """Mass thickness [g/cm^2] and density [g/cm^3] -> linear thickness [nm]."""
    return mass_thickness / density * 1e7


@dataclass
class Layer:
    """A homogeneous layer.

    composition     weight fractions by element symbol (need not sum to 1 —
                    unnormalized totals are meaningful in EPMA)
    mass_thickness  g/cm^2; None for the substrate
    density         g/cm^3; only needed to convert to/from linear thickness
    """

    composition: dict[str, float]
    mass_thickness: float | None = None
    density: float | None = None
    name: str | None = None

    def __post_init__(self):
        self.composition = _normalized(self.composition)
        if self.mass_thickness is not None and self.mass_thickness < 0:
            raise ValueError("mass_thickness must be >= 0")

    @classmethod
    def from_formula(cls, formula: str, **kwargs) -> "Layer":
        """Layer from a chemical formula, e.g. Layer.from_formula('Al2O3', density=3.95)."""
        parsed = xraylib.CompoundParser(formula)
        comp = {xraylib.AtomicNumberToSymbol(z): w
                for z, w in zip(parsed["Elements"], parsed["massFractions"], strict=True)}
        kwargs.setdefault("name", formula)
        return cls(comp, **kwargs)

    @classmethod
    def from_nm(cls, composition: Mapping[str, float] | str, thickness_nm: float,
                density: float, **kwargs) -> "Layer":
        """Film given its linear thickness [nm] and density [g/cm^3]."""
        mt = nm_to_mass_thickness(thickness_nm, density)
        if isinstance(composition, str):
            return cls.from_formula(composition, mass_thickness=mt, density=density, **kwargs)
        return cls(dict(composition), mass_thickness=mt, density=density, **kwargs)

    @property
    def total(self) -> float:
        return sum(self.composition.values())

    @property
    def thickness_nm(self) -> float | None:
        if self.mass_thickness is None or self.density is None:
            return None
        return mass_thickness_to_nm(self.mass_thickness, self.density)


@dataclass
class Sample:
    """A stack of layers, surface first; the last layer is the substrate."""

    layers: list[Layer]

    def __post_init__(self):
        if not self.layers:
            raise ValueError("a sample needs at least one layer")
        for i, layer in enumerate(self.layers[:-1]):
            if layer.mass_thickness is None:
                raise ValueError(f"layer {i} is not the substrate but has no mass_thickness")

    @classmethod
    def bulk(cls, composition: Mapping[str, float] | str, **kwargs) -> "Sample":
        if isinstance(composition, str):
            return cls([Layer.from_formula(composition, **kwargs)])
        return cls([Layer(dict(composition), **kwargs)])

    @property
    def substrate(self) -> Layer:
        return self.layers[-1]

    @property
    def films(self) -> list[Layer]:
        return self.layers[:-1]

    @property
    def is_bulk(self) -> bool:
        return len(self.layers) == 1

    @property
    def elements(self) -> set[str]:
        return {el for layer in self.layers for el in layer.composition}

    def boundaries(self) -> list[float]:
        """Mass depth [g/cm^2] of the top of each layer (first is 0)."""
        depths, depth = [], 0.0
        for layer in self.layers:
            depths.append(depth)
            depth += layer.mass_thickness or 0.0
        return depths


@dataclass
class Measurement:
    """An experimental k-ratio.

    standard defaults to the pure element of `line` (bulk).
    uncertainty is the 1-sigma absolute uncertainty on k.
    """

    line: XrayLine
    beam_energy: float
    kratio: float
    uncertainty: float | None = None
    standard: Sample | None = None
    takeoff_angle: float = 40.0
    meta: dict = field(default_factory=dict)

    def __post_init__(self):
        if self.standard is None:
            self.standard = Sample.bulk({self.line.element: 1.0})
        if self.line.element not in self.standard.elements:
            raise ValueError(f"standard does not contain {self.line.element}")
        if self.beam_energy <= self.line.edge:
            raise ValueError(f"beam energy {self.beam_energy} keV is below the "
                             f"{self.line} edge ({self.line.edge:.3f} keV)")
