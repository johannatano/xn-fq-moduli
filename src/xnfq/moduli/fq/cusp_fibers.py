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


CuspDatum = tuple[int, NeronDgonFq, list[int]]

class CuspFqFiber(CurveFiber):
    """The cusp contribution for a chosen level-N moduli problem."""

    def __init__(
        self,
        gamma: "LevelStructure",
        t: int,
        dgon: NeronDgonFq,
        eigen_vals: list[int],
    ):
        super().__init__(gamma)
        self.t = t
        self.dgon = dgon
        self.eigen_vals = eigen_vals

    def count(self) -> Fraction:
        val = 0

        

        if self.gamma.type == 2:
            val = phi(-1)(
                self.gamma.N
            )  # we have full torsion, if we are here d=N, e = 1 etc, no mixed heights
        else:
            # count varying stable lines counts from eigen values, depends on d
            e = self.gamma.N // self.dgon.d
            g = gcd(e, self.dgon.d)
            val = (
                len(self.eigen_vals) * Fraction(self.dgon.d, g) * phi(1)(g)
            )  # num stable lines count

            print(
                f"d={self.dgon.d}, s={self.dgon.t} lines={val} scaled={Fraction(val, self.dgon.d*2) * self.gamma.weight} g={g}, inv={self.dgon.fixed_subgroup(self.gamma.N)}| #Gamma{self.gamma.type}={val}, eigen_vals={self.eigen_vals}"
            )
        return Fraction(val, self.dgon.d*2) * self.gamma.weight

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
    strata: list[CuspDatum] = []

    def get_eigenvalues_over_d(d:int, t:int, n: int) -> list[int]:
        e = n // d
        result = []
        lambdas = range(n) if level_type == 0 else [1]
        for lam in lambdas:
            if (lam - t * curve.q) % e == 0 and (lam - t) % d == 0:
                result.append(lam)
        return result

    divs = divisors(curve.N) if level_type != 2 else [curve.N]
    for d in divs:
        for s in (-1, 1):
            eigen_vals = get_eigenvalues_over_d(d, s, curve.N)
            if len(eigen_vals) == 0:
                continue
            strata.append((s, NeronDgonFq(d, curve.q, s), eigen_vals))
    return strata

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
