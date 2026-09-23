from __future__ import annotations

from fractions import Fraction
from math import gcd

from ...arithmetic.algebra import NeronDgon
from ...arithmetic.common import divisors
from ...arithmetic.function import phi
from ..data import EigenForm, EigenFormRecord, FiberRecord
from ..modular_curve import CurveFiber
from .data import FrobData


CuspDatum = tuple[int, FrobData]

class CuspFiberFq(CurveFiber):
    """One cusp trace stratum containing one form per ``(lambda, d)``."""

    def __init__(
        self,
        gamma: "LevelStructure",
        t: int,
        frob_data: FrobData,
    ):
        super().__init__(gamma)
        self.t = t
        self.frob_data = frob_data

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
        return sum(self.level_counts().values(), Fraction(0))

    def level_counts(self) -> dict[int, Fraction]:
        """Return the contribution of each polygon level in this stratum."""

        counts: dict[int, Fraction] = {}
        for eigen_form in self.frob_data.eigen_forms:
            d = eigen_form.form.d
            g = gcd(self.gamma.N // d, d)
            contribution = Fraction(1, g) * phi(1)(g)
            counts[d] = counts.get(d, Fraction(0)) + contribution
        weight = Fraction(self.gamma.weight(smooth=False), 2)
        return {d: value * weight for d, value in counts.items()}

    def snapshot(self) -> FiberRecord:
        total = self.count()
        level_counts = self.level_counts()
        eigen_records = tuple(
            EigenFormRecord(eigenform=eigen_form)
            for eigen_form in self.frob_data.eigen_forms
        )
        return FiberRecord(
            kind="cusp",
            t=self.t,
            d=None,
            total_count=total,
            normalized_count=total,
            total_mass=total,
            eigenvalues=tuple(
                eigen_form.eigenvalue
                for eigen_form in self.frob_data.eigen_forms
            ),
            lattice_levels=tuple(sorted(level_counts)),
            lattice_counts=tuple(sorted(level_counts.items())),
            eigen_records=eigen_records,
        )

def enum_d_gons(
    curve: "ModularCurveFq",
    level_type: int,
) -> list[CuspDatum]:
    """Enumerate cusp strata compatible with the chosen level type."""
    strata: list[CuspDatum] = []
    lambdas = range(curve.N) if level_type == 0 else [1]

    for t in (1, -1):
        frob_data = FrobData(trace=t)
        for lam in lambdas:
            d1 = gcd(lam - t * curve.q, curve.N)
            d2 = gcd(lam - t, curve.N)
            if d1*d2 % curve.N != 0:
                continue
            # reverse iteration for the e | q-lambda
            """for d in divisors(d2):
                e = curve.N // d
                if d1 % e != 0:
                    continue
                print(f"--t={t}, lam={lam}, inv={d1, d2}, d={d}, e={e}")"""
            # we instead simply enumerate the divisible by d1
            m = curve.N // d1
            for k in divisors(d2 // m):
                d = m * k
                e = curve.N // d
                frob_data.eigen_forms.append(EigenForm(value=lam, form=NeronDgon(d)))

        if frob_data.eigen_forms:
            strata.append((t, frob_data))
    return strata
