from __future__ import annotations

from fractions import Fraction
from math import gcd, lcm
from typing import TYPE_CHECKING

from ...arithmetic.algebra import NeronDgon
from ...arithmetic.common import divisors, factorize, valuation as vl
from ...arithmetic.function import phi
from ..data import FiberRecord
from ..modular_curve import CurveFiber
from ...utils.logging import Colors, Logger

if TYPE_CHECKING:
    from .curve import ModularCurveFq


class NeronDgonFq(NeronDgon):
    def __init__(self, d: int, q: int, t:int):
        super().__init__(d)
        self.q = q
        self.t = t

    def stable_subgroup(self, n: int):
        if n % self.d != 0:
            return (0, 0)
        e = n // self.d
        g = gcd(self.d, e)
        if (self.q - 1) % g != 0:
            return (0, 0)
        # Frob stable subgroup, diamond
        return (1,g)

    def fixed_subgroup(self, n: int) -> "NeronDgon":
        """Return the n-torsion subgroup of the NeronDgon."""
        if n % self.d != 0:
            return (0,0)
        e = n // self.d
        # If split type 1 we use defualt, splitt ype -1 fixed by Frob twisted
        if self.t == 1:
            return (gcd(e, self.q-1), self.d)
        else:
            return (gcd(e, self.q+1), gcd(self.d, 2))

    def invariants(self):
        if self.t == 1:
            a, b = self.q - 1, self.d
        else:
            a, b = self.q + 1, gcd(self.d, 2)
        n1, n2 = gcd(a, b), lcm(a, b)  # invariant factors, n1 | n2, n1*n2 = a*b
        return n1, n2


CuspDatum = tuple[int, NeronDgonFq]

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
        n1, n2 = self.dgon.fixed_subgroup(self.gamma.N) if self.gamma.type > 0 else self.dgon.stable_subgroup(self.gamma.N)
        if self.gamma.type > 0 and n1 * n2 < self.gamma.N:
            return Fraction(0, 1)
        val = phi(1)(n1) * phi(1)(n2)
        if self.gamma.type == 2: # if we are here then d = N and e | q-1 hence full incl
            val *= phi(-1)(self.gamma.N)
        clr = Colors.GREEN if val > 0 else Colors.RED
        Logger.cprint(
            f"d={self.dgon.d}, s={self.dgon.t} | C^sm_d = {(n1, n2)}, #Gamma{self.gamma.type}={val}",
            clr,
        )
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

def enum_d_gons(
    curve: "ModularCurveFq",
    level_type: int,
) -> list[CuspDatum]:
    """Enumerate cusp strata compatible with the chosen level type."""
    dgons: list[CuspDatum] = []
    divs = divisors(curve.N) if level_type != 2 else [curve.N]
    for d in divs:
        for s in (-1, 1):
            dgons.append((s, NeronDgonFq(d, curve.q, s)))
    return dgons

"""
def enum_d_gons_0(curve):  # Gamma0(N)
    dgons = []
    _sum = 0
    _sum2 = 0
    for d in divisors(curve.N):
        e = curve.N // d
        g = gcd(d, e)
        if (curve.q - 1) % g == 0:
            # diamond operator
            dgons.append((1, NeronDgon(d)))
            dgons.append((-1, NeronDgon(d)))
            print(f"d={d}, e={e}, g={g}, d1={gcd(g, curve.q - 1)}")
            val = phi(1)(g)
            _sum += val
        for s in (-1, 1):
            dgon = NeronDgonFq(d, curve.q, s)
            n1, n2 = dgon.fixed_subgroup(curve.N)
            if n1 * n2 < curve.N:
                continue
            _sum2 += phi(1)(n1) * phi(1)(n2)
    _sum2 *= Fraction(1, phi(1)(curve.N)*2)
    print(_sum2)
    return _sum"""
