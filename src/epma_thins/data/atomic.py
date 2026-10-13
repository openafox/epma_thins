"""Element-level atomic data (backed by xraylib)."""

from functools import cache

import xraylib

SHELLS = ("K", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5")


def atomic_number(element: str | int) -> int:
    """Atomic number from a symbol ('Cu') or pass through an int."""
    if isinstance(element, int):
        if not 1 <= element <= 103:
            raise ValueError(f"atomic number out of range: {element}")
        return element
    z = xraylib.SymbolToAtomicNumber(element.strip().capitalize())
    if z <= 0:
        raise ValueError(f"unknown element symbol: {element!r}")
    return z


def symbol(element: str | int) -> str:
    return xraylib.AtomicNumberToSymbol(atomic_number(element))


@cache
def atomic_weight(element: str | int) -> float:
    """Standard atomic weight [g/mol]."""
    return xraylib.AtomicWeight(atomic_number(element))


@cache
def edge_energy(element: str | int, shell: str) -> float:
    """Ionization (absorption edge) energy of a shell [keV]."""
    if shell not in SHELLS:
        raise ValueError(f"unknown shell {shell!r}; expected one of {SHELLS}")
    energy = xraylib.EdgeEnergy(atomic_number(element), getattr(xraylib, f"{shell}_SHELL"))
    if energy <= 0:
        raise ValueError(f"no {shell} edge for {symbol(element)}")
    return energy
