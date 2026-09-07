from __future__ import annotations

from fractions import Fraction
from math import gcd
from typing import TYPE_CHECKING

from ...arithmetic.algebra import NeronDgon
from ...arithmetic.common import divisors
from ...arithmetic.function import phi
from ..data import FiberRecord
from ..modular_curve import CurveFiber

if TYPE_CHECKING:
    from .curve import ModularCurveFq


CuspDatum = tuple[int, NeronDgon]


class CuspFqFiber(CurveFiber):
    """The cusp contribution for a chosen level-N moduli problem."""

    def __init__(
        self,
        gamma: "LevelStructure",
        t: int,
        dgon: NeronDgon,
    ):
        super().__init__(gamma)
        self.t = t
        self.dgon = dgon

    def count(self) -> Fraction:
        d = self.dgon.d
        if self.gamma.type == 0:
            g = gcd(d, self.gamma.N // d)
            return Fraction(phi(1)(g), 2)
        d = self.dgon.d
        Nd = self.gamma.N // d

        val = phi(1)(d) * phi(1)(Nd)
        return Fraction(val, 2)

    def snapshot(self) -> FiberRecord:
        total = self.count()
        return FiberRecord(
            kind="cusp",
            t=self.t,
            d=self.dgon.d,
            total_count=total,
            normalized_count=total,
            total_mass=total,
        )

def enum_d_gons_0(curve):  # Gamma0(N)
    """Enumerate split and nonsplit cusp data for `Gamma_0(N)`."""

    dgons = []
    for d in divisors(curve.N):
        g = gcd(d, curve.N // d)
        if (curve.q - 1) % g == 0:
            dgons.append((1, NeronDgon(d)))
            dgons.append((-1, NeronDgon(d)))
    return dgons

def enum_d_gons(
    curve: "ModularCurveFq",
    level_type: int,
) -> list[CuspDatum]:
    """Enumerate cusp strata compatible with the chosen level type."""

    dgons: list[CuspDatum] = []
    if level_type == 0:
        return enum_d_gons_0(curve)
    if level_type == 2:
        # TODO: implement enumeration for Gamma2(N)
        return []
    dgons: list[CuspDatum] = []
    for d in divisors(curve.N):
        Nd = curve.N // d
        if (curve.q - 1) % Nd == 0:
            dgons.append((1, NeronDgon(d)))
    for d in (1, 2):
        if curve.N % d != 0:
            continue
        Nd = curve.N // d
        if (curve.q + 1) % Nd == 0:
            dgons.append((-1, NeronDgon(d)))
    return dgons
