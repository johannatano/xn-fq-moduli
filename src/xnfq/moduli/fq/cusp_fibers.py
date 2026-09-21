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
        g = gcd(self.gamma.N // self.dgon.d, self.dgon.d)
        val = len(self.eigen_vals) * Fraction(1, g) * phi(1)(g) # simplified, removed d/d
        # TODO: for gamma type == 0, we do not have to search the eigenvals, CRT gives #eigenvals g always.
        # val = num_eigen_lines if self.gamma.type < 2 else (phi(-1)(self.gamma.N) if num_eigen_lines > 0 else 0)
        val_aut = Fraction(val, 2) * self.gamma.weight(smooth=False)
        return val_aut

    def snapshot(self) -> FiberRecord:
        total = self.count()
        """for ev in self.eigen_vals:
            d1 = gcd(e, t * q - ev)
            d2 = gcd(d, t - ev)
            stable_lines = self.stable_lines_count(d1,d2)
            print(
                f"d={self.dgon.d}, e={e}, lam={ev}, inv={d1,d2}, stable_lines={stable_lines}, self.eigen_vals={self.eigen_vals}"
            )
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
                )"""
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
    lambdas = range(curve.N) if level_type == 0 else [1]
    eigenvalues_by_stratum: dict[tuple[int, int], set[int]] = {}
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
                eigenvalues_by_stratum.setdefault((t, d), set()).add(lam)
    for (t, d), eigen_vals in sorted(eigenvalues_by_stratum.items()):
        strata.append((t, NeronDgonFq(d, curve.q, t), sorted(eigen_vals)))
    return strata
