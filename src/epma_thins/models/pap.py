"""The PAP φ(ρz) model of Pouchou & Pichoir [PAP91].

φ(ρz) is described by two parabolic branches joined at Rc:

    φ1(ρz) = A1 (ρz - Rm)^2 + B1     0  <= ρz <= Rc
    φ2(ρz) = A2 (ρz - Rx)^2          Rc <= ρz <= Rx

with its area F and the form parameters φ(0), Rm and Rx given by empirical
expressions of composition and overvoltage. Units: keV, g/cm^2.

Every quantity is computed from an explicit composition (weight fractions,
normalized internally) so that thin-film calculations can supply a different
weighted ("fictitious") composition for each parameter [PAP91]p53-55, [Moy21b]p3.
"""

import math
from collections.abc import Mapping
from dataclasses import dataclass

from epma_thins.data.atomic import atomic_number, atomic_weight
from epma_thins.data.lines import XrayLine

Composition = Mapping[str, float]

# Deceleration law f(V) = sum D_k V^P_k  [PAP91]p35e8 (J in keV)
_P1, _P2 = 0.78, 0.1


def _pk(j: float) -> tuple[float, float, float]:
    return (_P1, _P2, -(0.5 - 0.25 * j))


def _dk(j: float) -> tuple[float, float, float]:
    return (6.6e-6, 1.12e-5 * (1.35 - 0.45 * j**2), 2.2e-6 / j)


def _normalized(comp: Composition) -> dict[str, float]:
    total = sum(comp.values())
    if total <= 0:
        raise ValueError("composition must have a positive total")
    return {el: w / total for el, w in comp.items()}


# --- averaging rules ---------------------------------------------------------

def mean_z(comp: Composition) -> float:
    """Z̄ = Σ C_i Z_i  [PAP91]p60."""
    return sum(c * atomic_number(el) for el, c in _normalized(comp).items())


def backscatter_z(comp: Composition) -> float:
    """Z_b = (Σ C_i Z_i^0.5)^2  [PAP91]p59."""
    return sum(c * math.sqrt(atomic_number(el)) for el, c in _normalized(comp).items()) ** 2


def log_mean_z(comp: Composition) -> float:
    """Z_n with ln Z_n = Σ C_i ln Z_i  [PAP91]p60."""
    return math.exp(sum(c * math.log(atomic_number(el)) for el, c in _normalized(comp).items()))


def stopping_m(comp: Composition) -> float:
    """M = Σ C_i Z_i / A_i  [PAP91]p35e6."""
    return sum(c * atomic_number(el) / atomic_weight(el) for el, c in _normalized(comp).items())


def mean_ionization_potential(comp: Composition) -> float:
    """J [keV]: ln J = Σ (C_i Z_i/A_i) ln J_i / M  [PAP91]p35e6, J_i from Zeller e7."""
    comp = _normalized(comp)
    total = 0.0
    for el, c in comp.items():
        z = atomic_number(el)
        j_i = 1e-3 * z * (10.04 + 8.25 * math.exp(-z / 11.22))
        total += c * z / atomic_weight(el) * math.log(j_i)
    return math.exp(total / stopping_m(comp))


# --- area F -----------------------------------------------------------------

def hutchins_m(line: XrayLine) -> float:
    """Exponent m of the Hutchins cross section  [PAP91]p36."""
    shell = line.shell[0]
    if shell == "K":
        return 0.86 + 0.12 * math.exp(-((line.z / 5.0) ** 2))
    if shell == "L":
        return 0.82
    if shell == "M":
        return 0.78
    raise ValueError(f"no Hutchins exponent for {line} ({line.shell} shell)")


def ionization_cross_section(line: XrayLine, e0: float) -> float:
    """Hutchins Q(E0) ∝ ln(U0) / (U0^m Ec^2), constant omitted  [PAP91]p36e10."""
    u0 = e0 / line.edge
    return math.log(u0) / (u0 ** hutchins_m(line) * line.edge**2)


def deceleration_factor(comp: Composition, line: XrayLine, e0: float) -> float:
    """1/S = ∫ Q(E)/(dE/dρs) dE from Ec to E0  [PAP91]p36e11."""
    j = mean_ionization_potential(comp)
    m = hutchins_m(line)
    u0 = e0 / line.edge
    v0 = e0 / j
    total = 0.0
    for dk, pk in zip(_dk(j), _pk(j), strict=True):
        tk = 1.0 + pk - m
        total += dk * (v0 / u0) ** pk * (tk * u0**tk * math.log(u0) - u0**tk + 1.0) / tk**2
    return u0 / (v0 * stopping_m(comp)) * total


def backscatter_coefficient(comp: Composition) -> float:
    """Mean backscatter coefficient η̄(Z_b)  [PAP91]p59."""
    zb = backscatter_z(comp)
    return 1.75e-3 * zb + 0.37 * (1.0 - math.exp(-0.015 * zb**1.3))


def backscatter_factor(comp: Composition, u0: float) -> float:
    """R = 1 - η̄ W̄ [1 - G(U0)]  [PAP91]p37e12, p59."""
    eta = backscatter_coefficient(comp)
    w = 0.595 + eta / 3.7 + eta**4.55
    q = (2.0 * w - 1.0) / (1.0 - w)
    j_u = 1.0 + u0 * (math.log(u0) - 1.0)
    g = (u0 - 1.0 - (1.0 - u0 ** -(1.0 + q)) / (1.0 + q)) / ((2.0 + q) * j_u)
    return 1.0 - eta * w * (1.0 - g)


# --- form parameters ----------------------------------------------------------

def surface_ionization(comp: Composition, u0: float) -> float:
    """φ(0) = 1 + 3.3 (1 - U0^-r) η̄^1.2, r = 2 - 2.3 η̄  [PAP91]p60."""
    eta = backscatter_coefficient(comp)
    r = 2.0 - 2.3 * eta
    return 1.0 + 3.3 * (1.0 - u0**-r) * eta**1.2


def total_trajectory(comp: Composition, e0: float, ec: float) -> float:
    """R0: path length while slowing from E0 to Ec [g/cm^2]  [PAP91]p35e9, p61; [Moy21a]e29."""
    j = mean_ionization_potential(comp)
    total = sum(j ** (1.0 - pk) * dk * (e0 ** (1.0 + pk) - ec ** (1.0 + pk)) / (1.0 + pk)
                for dk, pk in zip(_dk(j), _pk(j), strict=True))
    return total / stopping_m(comp)


def ionization_range(comp: Composition, line: XrayLine, e0: float) -> float:
    """Rx = Q D R0 [g/cm^2]  [PAP91]p60-61."""
    u0 = e0 / line.edge
    zbar = mean_z(comp)
    zn = log_mean_z(comp)
    q0 = 1.0 - 0.535 * math.exp(-((21.0 / zn) ** 1.2)) - 2.5e-4 * (zn / 20.0) ** 3.5
    q = q0 + (1.0 - q0) * math.exp(-(u0 - 1.0) / (40.0 / zbar))
    d = 1.0 + 1.0 / u0 ** (zbar**0.45)
    return q * d * total_trajectory(comp, e0, line.edge)


def depth_of_maximum(zbar: float, u0: float, rx: float) -> float:
    """Rm = G1(Z̄) G2(U0) G3(U0, Z̄) Rx  [PAP91]p61."""
    g1 = 0.11 + 0.41 * math.exp(-((zbar / 12.75) ** 0.75))
    g2 = 1.0 - math.exp(-((u0 - 1.0) ** 0.35) / 1.19)
    g3 = 1.0 - math.exp(-(u0 - 0.5) * zbar**0.4 / 4.0)
    return g1 * g2 * g3 * rx


# --- the distribution --------------------------------------------------------

def _poly_exp_integral(c2: float, c1: float, c0: float, chi: float, a: float, b: float) -> float:
    """∫_a^b (c2 x^2 + c1 x + c0) e^(-χx) dx, stable for small χ."""
    if b <= a:
        return 0.0
    if chi * b <= 0.5:  # Taylor series of the exponential; converges fast here
        total, term = 0.0, 1.0
        for n in range(25):
            moments = sum(c * (b ** (k + n + 1) - a ** (k + n + 1)) / (k + n + 1)
                          for k, c in ((2, c2), (1, c1), (0, c0)))
            total += term * moments
            term *= -chi / (n + 1)
        return total

    def antiderivative(x: float) -> float:
        p = c2 * x * x + c1 * x + c0
        dp = 2.0 * c2 * x + c1
        return -math.exp(-chi * x) * (p / chi + dp / chi**2 + 2.0 * c2 / chi**3)

    return antiderivative(b) - antiderivative(a)


@dataclass(frozen=True)
class PAPDistribution:
    """A parabolic PAP φ(ρz) distribution (depths in g/cm^2)."""

    area: float  # F
    phi0: float  # φ(0)
    rm: float
    rc: float
    rx: float
    a1: float
    b1: float
    a2: float

    @classmethod
    def from_form_parameters(cls, area: float, phi0: float, rm: float, rx: float
                             ) -> "PAPDistribution":
        """Solve for Rc and the branch coefficients  [PAP91]p38e16-18, p61 Appendix 3."""
        d = (rx - rm) * (area - phi0 * rx / 3.0) * ((rx - rm) * area - phi0 * rx * (rm + rx / 3.0))
        rc = math.nan
        if d >= 0.0:
            rc = 1.5 * ((area - phi0 * rx / 3.0) / phi0 - math.sqrt(d) / (phi0 * (rx - rm)))
        if not rm < rc < rx:
            # No satisfactory root (d < 0, or Rc outside (Rm, Rx)) at low overvoltage:
            # drop the Rm(Rx) relation and impose d = 0 [PAP91]p61 Appendix 3
            rm = rx * (area - phi0 * rx / 3.0) / (area + phi0 * rx)
            rc = 1.5 * rm * (area + phi0 * rx) / (phi0 * rx)
        if not 0.0 < rm < rx:
            # [PAP91]p61-62 describes exponential forms for this case; not yet implemented
            raise NotImplementedError(
                f"φ(ρz) has no interior maximum (Rm={rm:.3g}, Rx={rx:.3g}); "
                "the very-low-overvoltage forms of PAP Appendix 3 are not implemented")
        a1 = phi0 / (rm * (rc - rx * (rc / rm - 1.0)))
        b1 = phi0 - a1 * rm**2
        a2 = a1 * (rc - rm) / (rc - rx)
        return cls(area, phi0, rm, rc, rx, a1, b1, a2)

    def phi(self, rhoz: float) -> float:
        if rhoz < 0.0 or rhoz >= self.rx:
            return 0.0
        if rhoz <= self.rc:
            return self.a1 * (rhoz - self.rm) ** 2 + self.b1
        return self.a2 * (rhoz - self.rx) ** 2

    def integral(self, start: float, stop: float, chi: float = 0.0) -> float:
        """∫ φ(ρz) e^(-χ ρz) dρz from `start` to `stop` (clipped to [0, Rx])."""
        start, stop = max(start, 0.0), min(stop, self.rx)
        if stop <= start:
            return 0.0
        total = 0.0
        lo, hi = start, min(stop, self.rc)
        if hi > lo:  # A1 (x - Rm)^2 + B1
            total += _poly_exp_integral(self.a1, -2 * self.a1 * self.rm,
                                        self.a1 * self.rm**2 + self.b1, chi, lo, hi)
        lo, hi = max(start, self.rc), stop
        if hi > lo:  # A2 (x - Rx)^2
            total += _poly_exp_integral(self.a2, -2 * self.a2 * self.rx,
                                        self.a2 * self.rx**2, chi, lo, hi)
        return total

    def emitted(self, chi: float) -> float:
        """F(χ) = ∫_0^∞ φ(ρz) e^(-χ ρz) dρz  [PAP91]p38e19."""
        return self.integral(0.0, self.rx, chi)


def distribution(line: XrayLine, e0: float, *, range_comp: Composition,
                 backscatter_comp: Composition, deceleration_comp: Composition,
                 rm_z: float) -> PAPDistribution:
    """φ(ρz) from the compositions used for each parameter.

    range_comp          -> Rx                [Moy21b] step i
    backscatter_comp    -> R and φ(0)        step ii
    deceleration_comp   -> 1/S (and so F)    step iii
    rm_z                -> Z used for Rm     step iv
    For a bulk sample all compositions are the same and rm_z = Z̄.
    """
    if e0 <= line.edge:
        raise ValueError(f"E0 = {e0} keV is not above the {line} edge ({line.edge} keV)")
    u0 = e0 / line.edge
    area = (backscatter_factor(backscatter_comp, u0)
            * deceleration_factor(deceleration_comp, line, e0)
            / ionization_cross_section(line, e0))
    phi0 = surface_ionization(backscatter_comp, u0)
    rx = ionization_range(range_comp, line, e0)
    rm = depth_of_maximum(rm_z, u0, rx)
    return PAPDistribution.from_form_parameters(area, phi0, rm, rx)


def bulk_distribution(comp: Composition, line: XrayLine, e0: float) -> PAPDistribution:
    """φ(ρz) for a homogeneous bulk sample."""
    return distribution(line, e0, range_comp=comp, backscatter_comp=comp,
                        deceleration_comp=comp, rm_z=mean_z(comp))
