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

    def stable_lines_count_ell(self, l: int, a:int, e1: int, e2: int, scalar_only=False) -> int:
        """Count stable lines in the local quotient defined by `alpha`."""
        if e1 + e2 < a:
            return 0
        if e1 * e2 == 0:
            return l**e2
        else:
            return phi(1)(l**e2)

    def stable_lines_count(
        self, d1: int, d2: int, scalar_only=False
    ) -> int:
        """Count stable lines in the local quotient defined by `alpha`."""
        if d1*d2 < self.gamma.N:
            return 0
        if d1 == 1 or d2 == 1:
            return d2
        else:
            g = gcd(d1, d2)
            return Fraction(d2, g) * phi(1)(g)

    def count(self) -> Fraction:
        N = self.gamma.N
        d = self.dgon.d
        t = self.t
        q = self.dgon.q
        e = self.gamma.N // self.dgon.d
        g = gcd(e, self.dgon.d)
        #m = g if self.gamma.type == 0 else 1
        # num_e = m * Fraction(self.dgon.d, g) * phi(1)(g)
        val = len(self.eigen_vals) * Fraction(self.dgon.d, g) * phi(1)(g)
        # TODO: for gamma type == 0, we do not have to search the eigenvals, CRT gives #eigenvals g always.
        # an optimized version, computationally faster but hides the structure
        # num_eigen_lines = len(self.eigen_vals) * Fraction(self.dgon.d, g)*phi(1)(g)
        # val = num_eigen_lines if self.gamma.type < 2 else (phi(-1)(self.gamma.N) if num_eigen_lines > 0 else 0)
        # prime power version
        """
        for l,a in factorize(N):
            _ev_val = 0
            lambdas = range(l**a) if self.gamma.type == 0 else [1]
            for ev in lambdas:
                e1 = min(vl(e, l), vl(t * q - ev, l))
                e2 = min(vl(d, l), vl(t - ev, l))
                stable_lines = self.stable_lines_count(l, a, e1=e1, e2=e2)
                _ev_val += stable_lines
                print(
                    f"d={self.dgon.d}, e={e}, lam={ev}, inv={e1,e2}, stable_lines={stable_lines}, _ev_val={_ev_val}"
                )
            _val *= _ev_val"""
        _val = 0
        lambdas = range(self.gamma.N) if self.gamma.type == 0 else [1]
        for ev in lambdas:
            d1 = gcd(e, t * q - ev)
            d2 = gcd(d, t - ev)
            stable_lines = self.stable_lines_count(d1,d2)
            _ell_val = 1
            for l,a in factorize(self.gamma.N):
                _ell_val *= self.stable_lines_count_ell(l, a, vl(d1, l), vl(d2, l))
            true_val = Fraction(self.dgon.d, g) * phi(1)(g) if d1 == e and d2 == d else 0
            clr = Colors.DIM if stable_lines == 0 else Colors.BOLD
            _val += stable_lines
            if stable_lines != 0:
                Logger.cprint(
                    f"d={self.dgon.d}, e={e}, g={g}, lam={ev}, inv={d1,d2}, #lines={stable_lines}, d/g*phi(1)({g})={true_val}, _ell_val={_ell_val}, num_eigens={len(self.eigen_vals)}",
                    clr,
                )
        val_aut = Fraction(_val, self.dgon.d * 2) * self.gamma.weight(smooth=False)
        clr = Colors.RED if _val != val else Colors.GREEN if _val > 0 else Colors.DIM
        compare_str = f"CORRECT VALUE: {val}, COMPUTED VALUE: {_val}" if _val != val else ""
        Logger.cprint(
            f"Final count for cusp at d={self.dgon.d}. (unweighted): {_val},  (weighted) = {val_aut} {compare_str}",
            clr,
        )

        return val_aut

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
            """if t == 1:
                # note: for t = -1, lam = 1 gives 1-(-1)=2, we need d | 2 hence d in (1,2)
                print(
                    f"split cusp d={d} checking={lam}, (e | q-1)={(curve.q - 1) % e == 0}, (lam - t)={(lam - t)}"
                )
            else:
                print(
                    f"non split cusp d={d} checking={lam}, (e | q-1)={(curve.q - 1) % e == 0}, (lam - t)={(lam - t)}"
                )"""
            if (lam - t * curve.q) % e == 0 and (lam - t) % d == 0:
                result.append(lam)
                """print(f"accepted lambda={lam} at d={d}, t={t}")"""
        return result
    divs = divisors(curve.N) if level_type != 2 else [curve.N]
    lambdas = range(curve.N) if level_type == 0 else [1]
    _count = 0
    for t in (1, -1):
        for lam in lambdas:
            d1 = gcd(lam - t * curve.q, curve.N)
            d2 = gcd(lam - t, curve.N)
            if d1*d2 % curve.N != 0:
                continue
            m = curve.N // d1
            for k in divisors(d2 // m):
                d = m * k
                e = curve.N // d
                g = gcd(e, d)
                _count += Fraction(1, g) * phi(1)(g)
                #print(
                #    f"d={d}, e={e}, lambda={lam}, t={t}, d1={d1}, d2={d2}, d={d}, _count={_count}, Fraction(1, g) * phi(1)(g)={Fraction(1, g) * phi(1)(g)}"
                #)
    _count *= Fraction(1, 2) * curve.level_structure.weight(smooth=False)
    print(f"_count={_count}")
    for d in divs:
        for s in (-1, 1):
            eigen_vals = get_eigenvalues_over_d(d, s, curve.N)
            # if len(eigen_vals) == 0:
            #    continue
            strata.append((s, NeronDgonFq(d, curve.q, s), eigen_vals))
    return strata
