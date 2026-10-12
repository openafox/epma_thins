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

    overrides: {(XrayLine, absorber symbol): MAC}
    """

    def __init__(self, base: MacSource, overrides: Mapping[tuple[XrayLine, str], float]):
        self.base = base
        self.overrides = {(line, symbol(absorber)): value
                          for (line, absorber), value in overrides.items()}

    def __call__(self, line: XrayLine, absorber: str | int) -> float:
        key = (line, symbol(absorber))
        if key in self.overrides:
            return self.overrides[key]
        return self.base(line, absorber)


def mixture_mac(source: MacSource, line: XrayLine, composition: Mapping[str, float]) -> float:
    """MAC of a material: weight-fraction average of the elemental MACs [cm^2/g]."""
    return sum(w * source(line, el) for el, w in composition.items())
