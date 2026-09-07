from __future__ import annotations

from fractions import Fraction
from math import comb, gcd
from typing import TYPE_CHECKING

from ...arithmetic.forms import BinaryQuadraticForm
from ...arithmetic.common import legendre
from ...config import get_pari
from ..level_structures import LevelStructure
from ..modular_curve import ModularCurve

if TYPE_CHECKING:
    from ...arithmetic.algebra import NeronDgon
    from ...arithmetic.quadratic import LatticeTower
    from .smooth_fibers import EigenTrace


def _hK(DK: int) -> int:
    """Class number of the quadratic order with discriminant `DK`."""

    pari = get_pari()
    if pari:
        return int(pari.qfbclassno(DK))
    return BinaryQuadraticForm.h(DK)

def _mK(DK: int, p: int) -> Fraction:
    return _hK(DK) if DK != 0 else (p - 1)

def _uK(DK: int) -> Fraction:
    return 2 * (2 if DK == -4 else 3 if DK == -3 else 12 if DK == 0 else 1)

def _d0(DK: int, p: int, t: int) -> int:
    return 1 - legendre(DK, p) if t % p == 0 else 1

WeilDatum = tuple[int, "LatticeTower", "EigenTrace"]
CuspDatum = tuple[int, "NeronDgon"]

class ModularCurveFq(ModularCurve):
    """The `F_q`-rational point count of a level-N modular curve."""

    def __init__(self, level_structure: LevelStructure, p: int, n: int):
        super().__init__(level_structure)
        self.p, self.n, self.q = p, n, p**n

    def change_base(self, p: int, n: int) -> "ModularCurveFq":
        """Rebuild the same level problem over a new finite field."""
        return type(self)(self.level_structure, p, n)

    def smooth_fibers(self):
        """Yield the non-cuspidal Frobenius strata over `F_q`."""
        from .smooth_fibers import WeilqFiber, enum_weil_q

        for t, tower, eigen_trace in enum_weil_q(self):
            mass = Fraction(
                _mK(tower.DK, self.p) * _d0(tower.DK, self.p, t), _uK(tower.DK)
            )
            if mass == 0:
                continue
            yield WeilqFiber(self.level_structure, t, tower, eigen_trace, mass)

    def cusps(self):
        """Yield cusp fibers allowed by the level structure over `F_q`."""
        from .cusp_fibers import CuspFqFiber, enum_d_gons

        for t, dgon in enum_d_gons(self, self.level_structure.type):
            yield CuspFqFiber(self.level_structure, t, dgon)

    def hk(self, t: int, q: int, k: int) -> int:
        """Trace polynomial for the `k`th symmetric power at Frobenius trace `t`."""
        return sum(
            comb(k - j, j) * (-q) ** j * t ** (k - 2 * j) for j in range(k // 2 + 1)
        )

    def tr_frob_symk(self, k: int) -> int:
        """Compute the trace of Frobenius on the `k`th symmetric power."""
        val = 0
        # TODO: not sure if this always holds?
        if k == 0:
            val = self.q + (1 if gcd(self.q, self.N) == 1 else 0)
        curves_term = sum(self.hk(c.t, self.q, k) * c.count() for c in self.smooth_fibers())
        cusp_term = sum((c.t ** (k + 2)) * c.count() for c in self.cusps())
        return val - curves_term - cusp_term
