"""Characteristic X-ray lines.

Lines may be given in Siegbahn notation ('Ka1', 'La1', 'Ma') or IUPAC notation
('K-L3', 'L3-M5'). 'Ka', 'Kb', 'La' and 'Lb' are xraylib's intensity-weighted
composites; 'Ma' is treated as Ma1 (M5-N7).
"""

from dataclasses import dataclass
from functools import cached_property

import xraylib

from epma_thins.data.atomic import atomic_number, edge_energy, symbol

# Siegbahn name -> (IUPAC transition, xraylib line macro)
_SIEGBAHN = {
    "Ka1": ("K-L3", "KL3_LINE"),
    "Ka2": ("K-L2", "KL2_LINE"),
    "Kb1": ("K-M3", "KM3_LINE"),
    "Kb3": ("K-M2", "KM2_LINE"),
    "Ka": ("K-L2,3", "KA_LINE"),
    "Kb": ("K-M2,3", "KB_LINE"),
    "La1": ("L3-M5", "L3M5_LINE"),
    "La2": ("L3-M4", "L3M4_LINE"),
    "Lb1": ("L2-M4", "L2M4_LINE"),
    "Lb2": ("L3-N5", "L3N5_LINE"),
    "Lb3": ("L1-M3", "L1M3_LINE"),
    "Lb4": ("L1-M2", "L1M2_LINE"),
    "Lg1": ("L2-N4", "L2N4_LINE"),
    "Lg2": ("L1-N2", "L1N2_LINE"),
    "Lg3": ("L1-N3", "L1N3_LINE"),
    "Ll": ("L3-M1", "L3M1_LINE"),
    "Ln": ("L2-M1", "L2M1_LINE"),
    "La": ("L3-M4,5", "LA_LINE"),
    "Lb": ("L-M/N", "LB_LINE"),
    "Ma1": ("M5-N7", "M5N7_LINE"),
    "Ma2": ("M5-N6", "M5N6_LINE"),
    "Ma": ("M5-N7", "M5N7_LINE"),
    "Mb": ("M4-N6", "M4N6_LINE"),
    "Mg": ("M3-N5", "M3N5_LINE"),
}
# Composite lines whose initial vacancy is not in the first characters of the IUPAC name
_COMPOSITE_SHELL = {"Ka": "K", "Kb": "K", "La": "L3", "Lb": "L2"}


def _normalize(line: str) -> tuple[str, str, str]:
    """Return (canonical name, shell, xraylib macro) for a line name."""
    text = line.strip()
    for name, (iupac, macro) in _SIEGBAHN.items():
        if text.lower() == name.lower():
            shell = _COMPOSITE_SHELL.get(name, iupac.split("-")[0])
            return name, shell, macro
    if "-" in text:
        initial, final = (s.strip().upper() for s in text.split("-", 1))
        macro = f"{initial}{final}_LINE"
        if hasattr(xraylib, macro):
            return f"{initial}-{final}", initial, macro
    raise ValueError(f"unknown X-ray line {line!r}")


@dataclass(frozen=True)
class XrayLine:
    """A characteristic X-ray line of an element, e.g. XrayLine('Cu', 'Ka1')."""

    element: str
    line: str

    def __post_init__(self):
        object.__setattr__(self, "element", symbol(self.element))
        name, _, _ = _normalize(self.line)
        object.__setattr__(self, "line", name)
        self.energy  # noqa: B018 -- validate that the line exists for this element

    @property
    def z(self) -> int:
        return atomic_number(self.element)

    @property
    def shell(self) -> str:
        """Shell whose ionization produces the line (e.g. 'K', 'L3')."""
        return _normalize(self.line)[1]

    @cached_property
    def energy(self) -> float:
        """Line energy [keV]."""
        macro = _normalize(self.line)[2]
        energy = xraylib.LineEnergy(self.z, getattr(xraylib, macro))
        if energy <= 0:
            raise ValueError(f"line {self.line} does not exist for {self.element}")
        return energy

    @property
    def edge(self) -> float:
        """Critical ionization energy Ec of the line's shell [keV]."""
        return edge_energy(self.z, self.shell)

    def __str__(self) -> str:
        return f"{self.element} {self.line}"
