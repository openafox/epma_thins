"""Emitted intensities and k-ratios (primary intensity only, no fluorescence yet)."""

import math
from collections.abc import Mapping

from epma_thins.data.lines import XrayLine
from epma_thins.data.mac import MacSource, mixture_mac, pap91_mac
from epma_thins.models import pap


def absorption_chi(mac: MacSource, line: XrayLine, comp: Mapping[str, float],
                   takeoff_angle: float) -> float:
    """χ = (μ/ρ) csc(θ) [cm^2/g]; composition normalized  [PAP91]p38."""
    total = sum(comp.values())
    normalized = {el: w / total for el, w in comp.items()}
    return mixture_mac(mac, line, normalized) / math.sin(math.radians(takeoff_angle))


def bulk_intensity(comp: Mapping[str, float], line: XrayLine, e0: float,
                   takeoff_angle: float, mac: MacSource) -> float:
    """Emitted primary intensity of `line` from a bulk sample, up to factors that
    cancel in a k-ratio: C_A · F(χ)  [PAP91]p38e19."""
    distribution = pap.bulk_distribution(comp, line, e0)
    chi = absorption_chi(mac, line, comp, takeoff_angle)
    return comp.get(line.element, 0.0) * distribution.emitted(chi)


def bulk_kratio(comp: Mapping[str, float], line: XrayLine, e0: float, *,
                takeoff_angle: float = 40.0, standard: Mapping[str, float] | None = None,
                mac: MacSource | None = None) -> float:
    """k-ratio of `line` for a bulk sample against a bulk standard (pure element by default)."""
    mac = mac or pap91_mac()
    standard = standard or {line.element: 1.0}
    return (bulk_intensity(comp, line, e0, takeoff_angle, mac)
            / bulk_intensity(standard, line, e0, takeoff_angle, mac))
