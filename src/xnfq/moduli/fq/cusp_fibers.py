from __future__ import annotations

from fractions import Fraction
from math import gcd

from ...arithmetic.common import divisors
from ...arithmetic.function import phi
from .data import (
    EigenForm,
    EigenFormRecord,
    CuspFiberRecordFq,
    LevelStructureRecord,
    CuspForm,
)
from ..modular_curve import CurveFiber

class CuspFiberFq(CurveFiber):
    def __init__(
        self,
        gamma: "LevelStructure",
        t: int,
        eigen_forms: list[EigenForm],
    ):
        super().__init__(gamma)
        self.t = t
        self.eigen_forms = eigen_forms
    def count(self) -> Fraction:
        val = 0
        for eigen_form in self.eigen_forms:
            form = eigen_form.form
            m = self.gamma.N // form.d1
            for k in divisors(form.d2 // m):
                d = m * k
                
                if self.gamma.scalar_only and d != self.gamma.N:
                    continue
                
                e = self.gamma.N // d
                g = gcd(self.gamma.N // d, d)
                # note, we cancel the d in the numerator, since the true aut size  is 2d
                val += Fraction(1, g) * phi(1)(g)

        return val * Fraction(1, 2) * self.gamma.weight(smooth=False)

    def snapshot(self) -> CuspFiberRecordFq:
        eigen_records = []
        for eigen_form in self.eigen_forms:
            form = eigen_form.form
            m = self.gamma.N // form.d1
            levels = []
            for k in divisors(form.d2 // m):
                d = m * k
                e = self.gamma.N // d
                g = gcd(self.gamma.N // d, d)
                num_lines = Fraction(d, g) * phi(1)(g)
                gamma_count = num_lines * self.gamma.weight(smooth=False) * Fraction(1, d)
                if d != self.gamma.N and self.gamma.scalar_only:
                    gamma_count = 0
                levels.append(LevelStructureRecord(inv=(e, d), mass=1, num_lines=num_lines, gamma_count=gamma_count))
            eigen_records.append(
                EigenFormRecord(
                    eigenform=eigen_form,
                    levels=levels
                )
            )
        return CuspFiberRecordFq(
            t=self.t,
            m0 = Fraction(1, 2),
            count=self.count(),
            eigen_records=tuple(eigen_records),
        )

def enum_d_gons(
    curve: "ModularCurveFq",
    level_type: int,
) -> list[tuple[int, list[EigenForm]]]:
    strata: list[tuple[int, list[EigenForm]]] = []
    lambdas = range(curve.N) if level_type == 0 else [1]

    def is_valid(d1:int, d2:int) -> bool:
        if curve.level_structure.type < 2:
            return d1 * d2 % curve.N == 0
        return d2 == curve.N

    for t in (1, -1):
        eigen_forms = list[EigenForm]()
        for lam in lambdas:
            d1 = gcd(lam - t * curve.q, curve.N)
            d2 = gcd(lam - t, curve.N)
            if d1 * d2 % curve.N != 0:
                continue
            eigen_forms.append(EigenForm(value=lam, form=CuspForm(d1, d2)))

        if eigen_forms:
            strata.append((t, eigen_forms))
    return strata
