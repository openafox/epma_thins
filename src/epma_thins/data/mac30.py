"""Heinrich's MAC30 parameterization of mass absorption coefficients [Heinrich87].

Ported from the NIST EPQ library (`MassAbsorptionCoefficient.HeinrichBase`,
github.com/usnistgov/EPQ), which notes it is "Heinrich's formula as implemented
by Myklebust". EPQ notice: "This software was developed by employees of the
National Institute of Standards and Technology (NIST), an agency of the Federal
Government. Pursuant to title 17 United States Code Section 105, works of NIST
employees are not subject to copyright protection in the United States and are
considered to be in the public domain." Provided "as is" without warranty.

Differences from EPQ: edge energies come from xraylib rather than EPQ's own
table, which changes results only for lines within a few eV of an edge.
TODO: verify coefficients against [Heinrich87] once a copy is obtained.
"""

import math

import xraylib

from epma_thins.data.atomic import atomic_number, atomic_weight
from epma_thins.data.lines import XrayLine


def _edge_ev(z: int, shell: str) -> float:
    """Edge energy [eV], or NaN if the shell does not exist (EPQ semantics)."""
    try:
        energy = xraylib.EdgeEnergy(z, getattr(xraylib, f"{shell}_SHELL"))
    except ValueError:
        return math.nan
    return energy * 1000.0 if energy > 0 else math.nan


def mac30_energy(z: int, energy_ev: float, dtsa: bool = False) -> float:
    """MAC30 for photons of `energy_ev` [eV] in element `z` [cm^2/g].

    dtsa=True reproduces the extra high-Z conditions and cutoff of the DTSA
    implementation, which are not in the ICXOM-11 article.
    """
    if energy_ev <= 10.0:
        return 1e6
    if z < 3 or z > 95:
        return 0.001

    nm = cc = az = bias = 0.0
    ee_k = _edge_ev(z, "K")
    ee_n1 = _edge_ev(z, "N1")

    if energy_ev > ee_k:  # above the K edge
        if z < 6:
            cc = 1.808599e-3 * z - 2.87536e-4
            az = (-14.15422 * z + 155.6055) * z + 24.4545
            bias = 18.2 * z - 103
            nm = (-0.01273815 * z + 0.02652873) * z + 3.34745
        else:
            cc = 5.253e-3 + z * (1.33257e-3 + z * (-7.5937e-5 + z * (1.69357e-6 + -1.3975e-8 * z)))
            az = ((-0.152624 * z + 6.52) * z + 47) * z
            nm = 3.112 - 0.0121 * z
            if dtsa:
                if z >= 50:
                    az = ((-0.015 * z + 3.52) * z + 47) * z
                if z >= 57:
                    cc = 2.0e-4 + (1.0e-4 - z) * z
    else:
        ee_l3 = _edge_ev(z, "L3")
        if energy_ev > ee_l3:  # between the L3 and K edges
            cc = -0.0924e-3 + z * (0.141478e-3 + z * (-0.00524999e-3 + z * (
                9.85296e-8 + z * (-9.07306e-10 + z * 3.19245e-12))))
            az = (((-1.16286e-4 * z + 0.01253775) * z + 0.067429) * z + 17.8096) * z
            nm = (-4.982e-5 * z + 1.889e-3) * z + 2.7575
            ee_l2 = _edge_ev(z, "L2")
            if energy_ev < _edge_ev(z, "L1") and energy_ev >= ee_l2:
                cc *= 0.858
            if energy_ev < ee_l2:
                cc *= 0.8933 + z * (-8.29e-3 + 6.38e-5 * z)
        else:
            ee_m1 = _edge_ev(z, "M1")
            if ee_m1 < energy_ev <= ee_l3:  # between the M1 and L3 edges
                nm = (((4.4509e-6 * z - 1.08246e-3) * z + 0.084597) * z) + 0.5385
                if z < 30:
                    cc = ((((7.2773258e-9 * z - 1.1641145e-6) * z + 6.9602789e-5) * z
                           - 1.8517159e-3) * z) + 1.889757e-2
                else:
                    cc = ((((1.497763e-10 * z - 4.0585911e-8) * z + 4.0424792e-6) * z
                           - 1.73663566e-4) * z) + 3.0039e-3
                az = (((-1.8641019e-4 * z + 2.63199611e-2) * z - 0.822863477) * z
                      + 10.2575657) * z
                if z < 61:
                    bias = (((-1.683474e-4 * z + 0.018972278) * z - 0.536839169) * z + 5.654) * z
                else:
                    bias = (((3.1779619e-3 * z - 0.699473097) * z + 51.114164) * z
                            - 1232.4022) * z
            else:
                ee_m5 = _edge_ev(z, "M5")
                if energy_ev >= ee_m5:  # between the M5 and M1 edges
                    ee_m4 = _edge_ev(z, "M4")
                    az = (4.62 - 0.04 * z) * z
                    cc = (((-1.29086e-9 * z + 2.209365e-7) * z - 7.83544e-6) * z) + 7.7708e-5
                    cc *= (((4.865e-6 * z - 0.0006561) * z + 0.0162) * z) + 1.406
                    bias = ((3.78e-4 * z - 0.052) * z + 2.51) * ee_m4
                    nm = 3 - 0.004 * z
                    if energy_ev >= _edge_ev(z, "M2"):
                        cc *= (-0.0001285 * z + 0.01955) * z + 0.584
                    elif energy_ev >= _edge_ev(z, "M3"):
                        cc *= 0.001366 * z + 1.082
                    elif energy_ev >= ee_m4:
                        cc *= 0.95
                    else:
                        cc *= (4.0664e-4 * z - 4.8e-2) * z + 1.6442
                else:  # below the M5 edge
                    cc = 1.08 * ((((-6.69827e-9 * z + 1.707073e-6) * z - 1.4653e-4) * z)
                                 + 4.3156e-3)
                    az = ((5.39309e-3 * z - 0.61239) * z + 19.64) * z
                    bias = 4.5 * z - 113.0
                    nm = 0.3736 + 0.02401 * z

    base = cc * (12397.0 / energy_ev) ** nm * z**4 / atomic_weight(z)
    if energy_ev > ee_n1 or math.isnan(ee_n1):
        return base * (1 - math.exp((bias - energy_ev) / az))
    mu = base * (1 - math.exp((bias - ee_n1) / az))
    cutoff = 10.0 if dtsa else (0.252 * z - 31.1812) * z + 1042.0
    if energy_ev > cutoff:
        mu = 1.02 * mu * (energy_ev - cutoff) / (ee_n1 - cutoff)
    return mu


class Mac30:
    """MAC source using Heinrich's MAC30 parameterization."""

    def __init__(self, dtsa: bool = False):
        self.dtsa = dtsa

    def __call__(self, line: XrayLine, absorber: str | int) -> float:
        mu = mac30_energy(atomic_number(absorber), line.energy * 1000.0, self.dtsa)
        if not mu > 0:
            # e.g. B Ka in Zr: the parameterization fails just above some M5 edges
            raise ValueError(f"MAC30 is invalid for {line} in {absorber} ({mu:.4g} cm2/g); "
                             "supply a value with OverriddenMac")
        return mu
