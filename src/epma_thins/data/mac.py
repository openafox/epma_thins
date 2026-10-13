"""Mass absorption coefficients (MACs).

A MAC source maps (X-ray line, absorbing element) -> cm^2/g. Sources are
pluggable so results can be compared against programs using MAC30 (Heinrich
1987), PENELOPE/EPDL, or user-measured values.
"""

from collections.abc import Mapping
from typing import Protocol

import xraylib

from epma_thins.data.atomic import atomic_number, symbol
from epma_thins.data.lines import XrayLine
from epma_thins.data.mac30 import Mac30


class MacSource(Protocol):
    def __call__(self, line: XrayLine, absorber: str | int) -> float:
        """MAC of `line` radiation in pure `absorber` [cm^2/g]."""
        ...


class XraylibMac:
    """MACs from xraylib's cross-section tables.

    kind='total' gives the total mass attenuation coefficient (photoelectric +
    coherent + incoherent); kind='photo' gives the photoelectric part only,
    which is what PENELOPE-2018 / BadgerFilm use by default.
    """

    def __init__(self, kind: str = "total"):
        funcs = {"total": xraylib.CS_Total, "photo": xraylib.CS_Photo}
        if kind not in funcs:
            raise ValueError(f"kind must be one of {sorted(funcs)}")
        self.kind = kind
        self._cs = funcs[kind]

    def __call__(self, line: XrayLine, absorber: str | int) -> float:
        return self._cs(atomic_number(absorber), line.energy)


class OverriddenMac:
    """Wrap a MAC source with user-specified values for particular pairs.

    overrides: {(line, absorber): MAC} where line is an XrayLine or an
    (element, line name) pair. A family name such as ('Fe', 'La') applies to
    every line of that family (La1, La2, ...); an exact line takes precedence.
    """

    def __init__(self, base: MacSource,
                 overrides: Mapping[tuple[XrayLine | tuple[str, str], str], float]):
        self.base = base
        self.overrides = {}
        for (line, absorber), value in overrides.items():
            element, name = (line.element, line.line) if isinstance(line, XrayLine) else line
            self.overrides[(symbol(element), name, symbol(absorber))] = value

    def __call__(self, line: XrayLine, absorber: str | int) -> float:
        absorber = symbol(absorber)
        for name in (line.line, line.family):
            value = self.overrides.get((line.element, name, absorber))
            if value is not None:
                return value
        return self.base(line, absorber)


# [PAP91]p63 Appendix 5, first table: Ka of very light emitters. Values in
# parentheses in the original (B Ka in Co, Ni) are included as printed.
_PAP91_LIGHT = {
    "B": {"B": 3500, "C": 6750, "N": 11000, "O": 16500, "Al": 64000, "Si": 80000,
          "Ti": 15000, "V": 18000, "Cr": 20700, "Fe": 27800, "Co": 32000, "Ni": 37000,
          "Zr": 4400, "Nb": 4500, "Mo": 4600, "La": 2500, "Ta": 23000, "W": 21000, "U": 7400},
    "C": {"B": 39000, "C": 2170, "Si": 35000, "Ti": 8100, "V": 8850, "Cr": 10700,
          "Fe": 13500, "Zr": 25000, "Nb": 24000, "Mo": 20500, "Hf": 18000, "Ta": 17000,
          "W": 18000},
    "N": {"B": 15800, "N": 1640, "Al": 13800, "Si": 17000, "Ti": 4270, "V": 4950,
          "Cr": 5650, "Fe": 7190, "Zr": 24000, "Nb": 25000, "Mo": 25800, "Hf": 14000,
          "Ta": 15500},
}
# [PAP91]p63 Appendix 5, second table: lines near edges or with resonances.
# Au Mb in Pt is printed as 2550; NIST EPQ's Pouchou1991 table has 2250.
_PAP91_SPECIAL = {
    ("Si", "Ka", "Ta"): 1490, ("S", "Ka", "Au"): 2200, ("Cu", "Lb", "Cu"): 6750,
    ("As", "La", "Ga"): 7000, ("Mo", "La", "Au"): 2200, ("Gd", "Mb", "Gd"): 4700,
    ("Hf", "Mb", "Hf"): 3000, ("Ta", "Mb", "Ta"): 2500, ("W", "Mb", "W"): 2080,
    ("Au", "Mb", "Pt"): 2550, ("Hg", "Mb", "Au"): 2170, ("Sc", "La", "Sc"): 4750,
    ("Ti", "La", "Ti"): 4550, ("V", "La", "V"): 4370, ("Cr", "La", "Cr"): 3850,
    ("Mn", "La", "Mn"): 3340, ("Fe", "La", "Fe"): 3350, ("Co", "La", "Co"): 3260,
    ("Ni", "La", "Ni"): 3560, ("Cu", "La", "Cu"): 1755,
}

PAP91_MAC_OVERRIDES = {
    **{((emitter, "Ka"), absorber): value
       for emitter, table in _PAP91_LIGHT.items() for absorber, value in table.items()},
    **{((emitter, line), absorber): value
       for (emitter, line, absorber), value in _PAP91_SPECIAL.items()},
}


def pap91_mac() -> OverriddenMac:
    """The MAC set used with the PAP model: MAC30 plus the [PAP91] Appendix 5 values.

    PAP used Henke's values for other very light emitters (Z < 10); this falls
    back to MAC30 for those, as NIST EPQ does.
    """
    return OverriddenMac(Mac30(), PAP91_MAC_OVERRIDES)


def mixture_mac(source: MacSource, line: XrayLine, composition: Mapping[str, float]) -> float:
    """MAC of a material: weight-fraction average of the elemental MACs [cm^2/g]."""
    return sum(w * source(line, el) for el, w in composition.items())
