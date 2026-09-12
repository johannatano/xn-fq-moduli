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

    def stable_lines_count(self, l: int, a:int, e1: int, e2: int, scalar_only=False) -> int:
        """Count stable lines in the local quotient defined by `alpha`."""
        if e1 + e2 < a:
            return 0
        if e1 * e2 == 0:
            return l**e2
        else:
            return phi(1)(l**e2) 

    def count(self) -> Fraction:
        N = self.gamma.N
        d = self.dgon.d
        t = self.t
        q = self.dgon.q
        e = self.gamma.N // self.dgon.d
        g = gcd(e, self.dgon.d)

        
        print(f"d={d}, t={t}, e={e}, g={g}, #eigen={len(self.eigen_vals)}")
        
        num_eigen_lines = len(self.eigen_vals) * Fraction(self.dgon.d, g)*phi(1)(g)
        val = num_eigen_lines if self.gamma.type < 2 else (phi(-1)(self.gamma.N) if num_eigen_lines > 0 else 0)
        
        _val = 0
        _val_new = 1
        for l,a in factorize(N):
            _ev_val = 0
            lambdas = range(l**a) if self.gamma.type == 0 else [1]
            for ev in lambdas:
                e1 = min(vl(e, l), vl(t * q - ev, l))
                e2 = min(vl(d, l), vl(t - ev, l))
                n1 = min(e1, e2)
                n2 = max(e1, e2)
                stable_lines = self.stable_lines_count(l, a, e1=e1, e2=e2)
                _ev_val += stable_lines
                print(
                    f"d={self.dgon.d}, e={e}, lam={ev}, inv={e1,e2}, snf_inv={n1,n2}, stable_lines={stable_lines}, _ev_val={_ev_val}"
                )
            _val_new *= _ev_val

        """lambdas = range(N) if self.gamma.type == 0 else [1]
        for ev in lambdas:
            # we reduce (tq-lam)x = -ejk mod N ==> since e = N/d, we solve mod e and then lift to mod N
            d1 = gcd(e, t * q - ev)
            d2 = gcd(d, t - ev)
            #if d2 != d or d1 != e:
            #    continue # ample gate
            _lines = Fraction(self.dgon.d, g) * phi(1)(g)
            # Num lifted solutions (tq-lam)x = -ejk mod N with x in (Z/(e)Z)^x
            num_solutions = d * phi(1)(e)
            num_lines = Fraction(num_solutions*phi(1)(d), phi(1)(N))
            _val += num_lines
            print(
                f"d={self.dgon.d}, e={e}, lam={ev}, d2={d2}, num_solutions={num_solutions}, num_lines={num_lines}, _lines={_lines}"
            )"""
            

        _val *= Fraction(1, self.dgon.d * 2) * self.gamma.weight(smooth=False)
        _val_new *= Fraction(1, self.dgon.d * 2) * self.gamma.weight(smooth=False)
        clr = Colors.RED if _val_new != val_aut else Colors.GREEN
        Logger.cprint(f"Final: correc_val={val_aut}, _val_new={_val_new}", clr)
        return Fraction(val, self.dgon.d * 2) * self.gamma.weight(smooth=False)

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
            #if len(eigen_vals) == 0:
            #    continue
            strata.append((s, NeronDgonFq(d, curve.q, s), eigen_vals))
    return strata
