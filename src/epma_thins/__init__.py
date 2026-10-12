"""Thin-film and multilayer quantification for electron probe microanalysis.

Units used throughout the package:
    energies        keV
    mass thickness  g/cm^2  (helpers convert to/from nm given a density)
    density         g/cm^3
    composition     weight fractions
    angles          degrees
    MACs            cm^2/g
"""

from epma_thins.data.lines import XrayLine
from epma_thins.sample import Layer, Measurement, Sample

__all__ = ["Layer", "Measurement", "Sample", "XrayLine"]
__version__ = "0.1.0.dev0"
