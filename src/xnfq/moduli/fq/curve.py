from __future__ import annotations

from fractions import Fraction
from math import comb, gcd, isqrt

from ...arithmetic.forms import BinaryQuadraticForm
from ...arithmetic.common import legendre
from ..level_structures import LevelStructure
from ..modular_curve import ModularCurve
from .data import TrFqTraceRecord
from ...config import get_config, get_pari

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


class ModularCurveFq(ModularCurve):
    def __init__(self, level_structure: LevelStructure, p: int, n: int):
        super().__init__(level_structure)
        self.p, self.n, self.q = p, n, p**n

    def change_base(self, p: int, n: int) -> "ModularCurveFq":
        return type(self)(self.level_structure, p, n)

    def smooth_fibers(self):
        from .smooth_fibers import SmoothFiberFq, enum_weil_q
        config = get_config()
        HB = 2 * isqrt(self.q)
        fast_trace_enum = config.fast_trace and self.level_structure.type > 0
        total_enum = HB // self.N

        strata = enum_weil_q(self, fast_trace_enum)
        with self._progress_ctx(total=len(strata), desc="smooth fibers") as _p:
            for t, tower, frob_data in strata:
                mass = Fraction(
                    _mK(tower.DK, self.p) * _d0(tower.DK, self.p, t), _uK(tower.DK)
                )
                if mass == 0:
                    try:
                        _p.update(1)
                    except Exception:
                        pass
                    continue
                yield SmoothFiberFq(self.level_structure, t, tower, frob_data, mass)
                try:
                    _p.update(1)
                except Exception:
                    pass

    def cusps(self):
        from .cusp_fibers import CuspFiberFq, enum_d_gons
        strata = enum_d_gons(self, self.level_structure.type)
        with self._progress_ctx(total=len(strata), desc="cusp fibers") as _p:
            for t, frob_data in strata:
                yield CuspFiberFq(self.level_structure, t, frob_data)
                try:
                    _p.update(1)
                except Exception:
                    pass

    def hk(self, t: int, q: int, k: int) -> int:
        """Complete homogeneous polynomial for the `k`th symmetric power at Frobenius trace `t`."""
        return sum(
            comb(k - j, j) * (-q) ** j * t ** (k - 2 * j) for j in range(k // 2 + 1)
        )
    def tr_fq(self, k: int) -> TrFqTraceRecord:
        assert k >= 2, "tr_fq is only defined for k >= 2"
        eps0 = 0
        # TODO: double check eps0 for n > 1
        if k == 2: # we recover internal k = 0
            eps0 = self.q + (1 if gcd(self.q, self.N) == 1 else 0)
        curves_term = sum(self.hk(c.t, self.q, k-2) * c.count() for c in self.smooth_fibers())
        cusp_term = sum((c.t ** (k)) * c.count() for c in self.cusps())
        
        return TrFqTraceRecord(
            val=eps0 - curves_term - cusp_term,
            smooth=curves_term,
            cusp=cusp_term,
            eps0=eps0,
        )

    def info(self) -> str:
        base = super().info()
        return f"{base}(F_{self.p}^{self.n})"
