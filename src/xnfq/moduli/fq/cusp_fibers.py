from __future__ import annotations

from fractions import Fraction
from math import gcd
from typing import TYPE_CHECKING

from ...arithmetic.algebra import NeronDgon
from ...arithmetic.common import divisors
from ...arithmetic.function import phi
from ..data import FiberRecord
from ..modular_curve import CurveFiber
from ...utils.logging import Colors, Logger

if TYPE_CHECKING:
    from .curve import ModularCurveFq

CuspDatum = tuple[int, NeronDgonFq]


class NeronDgonFq(NeronDgon):
    def __init__(self, d: int, q: int, t:int):
        super().__init__(d)
        self.q = q
        self.type = t

    def torsion_subgroup(self, n: int) -> "NeronDgon":
        """Return the n-torsion subgroup of the NeronDgon."""
        if n % self.d != 0:
            return (0,0)
        if self.t == 1:
            return (gcd(n // self.d, self.q-1), self.d)
        else:
            return (gcd(n // self.d, self.q+1), gcd(self.d, 2))

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
        # d = self.dgon.d

        n1, n2 = self.dgon.torsion_subgroup(self.gamma.N)
        val = phi(1)(n1) * phi(1)(n2)
        #if self.gamma.type == 0:
        #    g = gcd(d, e)
        #    val = phi(1)(g)  # unconditionally rational — no q-test, ever
        # elif self.gamma.type == 2:
        #    val = phi(1)(self.gamma.N) * phi(-1)(self.gamma.N)
        # else:
        #    val = phi(1)(d) * phi(1)(e)
        return Fraction(val, 2)

    """def count(self) -> Fraction:
        d = self.dgon.d
        e = self.gamma.N // d
        # this is really d*phi(1)(e) * phi(1)(d) but d cancels in the aut
        val = phi(1)(e) * phi(1)(d)
        if self.gamma.type == 2:
            # note, the extra N factor cancels from the aut = d*2
            val *= phi(1)(self.gamma.N)
        print(
            f"count from d={d}, e={e}, val={val}, phi(1)(d)={phi(1)(d)}, phi(1)(e)={phi(1)(e)}"
        )
        return Fraction(val, 2)"""

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


def enum_d_gons_2(curve):  # Gamma0(N)
    """Enumerate split and nonsplit cusp data for `Gamma_0(N)`."""
    dgons = []
    d = curve.N

    print(f"curve.N={curve.N}, curve.q={curve.q}, d={d}, ((curve.q - 1) % d)={(curve.q - 1) % d}")
    if (curve.q - 1) % d == 0:
        dgons.append((1, NeronDgon(d)))
        print("Appending split cusp for d =", d)
        if (curve.q + 1) % d == 0:
            # N == 2 case
            dgons.append((-1, NeronDgon(d)))
    return dgons

def enum_d_gons_0(curve):  # Gamma0(N)
    """Enumerate split and nonsplit cusp data for `Gamma_0(N)`."""
    dgons = []
    for d in divisors(curve.N):
        g = gcd(d, curve.N // d)
        if (curve.q - 1) % g == 0:
            # diamond operator
            dgons.append((1, NeronDgon(d)))
            dgons.append((-1, NeronDgon(d)))
            print(
                f"Appending cusp for d={d}, g={g}, curve.q={curve.q}, (curve.q - 1)={(curve.q - 1)% g}"
            )
        else:
            print(f"Not appending cusp for d={d}, g={g}, curve.q={curve.q}")
    return dgons

def enum_d_gons(
    curve: "ModularCurveFq",
    level_type: int,
) -> list[CuspDatum]:
    """Enumerate cusp strata compatible with the chosen level type."""
    # dgons: list[CuspDatum] = []
    # if level_type == 0:
    #    return enum_d_gons_0(curve)
    # if level_type == 2:
    #    # TODO: implement enumeration for Gamma2(N)
    #    return enum_d_gons_2(curve)
    #    return []
    dgons: list[CuspDatum] = []

    # splits = []
    # non_splits=[]
    # TODO: later we optimize, now just enum all
    #divs = divisors(curve.N) # if level_type != 2 else [curve.N]
    for d in divisors(curve.N):
        for s in (-1, 1):
            dgon = NeronDgonFq(d, curve.q, s)
        # the N torsion part is C_d[N] = mu_(N/d) x Z/dZ
        # n1, n2 = dgon.torsion_subgroup(curve.N)
        # we want to know, does anything survive in the Gm component, in particular
        # for the Frob -> Id: gcd(N/d, q-1)
        e = gcd(n1, q-1)
        # for the Frob -> -1: gcd(N/d, q+1) and d in Z/2Z
        # if self.gamma.type == 0: # we have looser crit for cyclic subgrps
        #    e = gcd(d, curve.N // d)
        # the frob split type 1 gives then: mu_(N/d) in Fq
        if (curve.q - 1) % e == 0:
            dgons.append((1, NeronDgonFq(d, curve.q, 1)))
        # the frob non split gives mu(N/d) in Uq+1: that is: if (curve.q + 1) % e == 0: and d in Z/2Z
        if (curve.q + 1) % e == 0 and d in (1, 2):
            dgons.append((-1, NeronDgonFq(d, curve.q, -1)))
        """if (curve.q - 1) % e == 0:
            dgons.append((1, NeronDgonFq(d, curve.q, 1)))
            splits.append(d)
        if d in (1, 2):
            print(
                f"Checking non-split cusp for d={d}, e={e} (curve.q + 1)% e={(curve.q + 1)% e}"
            )
            if (curve.q + 1) % e == 0:
                dgons.append((-1, NeronDgonFq(d, curve.q, -1        )))
                non_splits.append(d)"""
        # for s in (-1,1):
        # if (curve.q-s) % Nd == 0:
        # print(f"-------------TEST Appending cusp for d={d}, s={s}, Nd={Nd}, curve.q={curve.q}")
    # print(f"Enumerated dgons: {len(dgons)}")
    for s, neron in dgons:
        print(f"Enumerated cusp: {neron.d} with sign={s}")
    return dgons
