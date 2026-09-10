from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from math import gcd, isqrt
from typing import Protocol
import math

from .common import factorize

Exponent = tuple[int, ...]


class SupportsQuadDisc(Protocol):
    def quaddisc(self, D: int) -> int: ...

class Form:
    def __init__(self, rank: int, degree: int, coeffs) -> None:
        self.rank = rank
        self.degree = degree
        self.coeffs = self._normalize_coeffs(coeffs)

    def _normalize_coeffs(self, coeffs) -> dict[Exponent, int | Fraction]:
        normalized: dict[Exponent, int | Fraction] = {}
        for exp, coeff in dict(coeffs).items():
            key = tuple(exp)
            if coeff:
                normalized[key] = coeff
        return normalized

    @property
    def content(self) -> int:
        values = [coeff for coeff in self.coeffs.values() if coeff]
        if not values:
            return 0
        numerators = [Fraction(coeff).numerator for coeff in values]
        return math.gcd(*numerators)

    def is_primitive(self) -> bool:
        return self.content == 1


class BinaryQuadraticForm(Form):
    """Helper functions, naive iteration of finding representatives of the SL2(Z) orbits, that is the reduced forms, as well as the primitive representations"""

    @staticmethod
    def _D0(D: int, pari: SupportsQuadDisc | None = None) -> tuple[int, int]:
        """Return (DK, f) where DK is the fundamental discriminant and D = DK * f^2."""
        if D == 0:
            return 0, 0
        if pari is not None:
            DK = pari.quaddisc(D)
            f = isqrt(abs(D) // abs(DK))
            return int(DK), int(f)
        sign = -1 if D < 0 else 1
        sf = 1

        for p, e in factorize(abs(D)):
            if e % 2 == 1:
                sf *= p
        sf *= sign
        DK = sf if sf % 4 == 1 else 4 * sf
        f = math.isqrt(abs(D) // abs(DK))
        return DK, f

    def reduced_triples(D: int) -> list[tuple[int, int, int]]:
        """B~(D) of Schoof (2.5)."""
        out = []
        amax = isqrt(abs(D) // 3)
        for a in range(1, amax + 1):
            for b in range(-a, a + 1):                 # |b| <= a
                if (b * b - D) % (4 * a):
                    continue
                c = (b * b - D) // (4 * a)
                if c < a:                              # a <= c
                    continue
                if b < 0 and (a == -b or a == c):      # b >= 0 on the boundary
                    continue
                out.append((a, b, c))
        return out

    def H(D: int) -> int:
        return len(BinaryQuadraticForm.reduced_triples(D))

    def h(D: int) -> int:
        return sum(
            1
            for a, b, c in BinaryQuadraticForm.reduced_triples(D)
            if gcd(gcd(a, b), c) == 1
        )

    def __init__(self, A: int, B: int, C: int) -> None:
        self._A = A
        self._B = B
        self._C = C
        super().__init__(2, 2, {(2, 0): A, (1, 1): B, (0, 2): C})

    @property
    def A(self) -> int:
        return self._A

    @property
    def B(self) -> int:
        return self._B

    @property
    def C(self) -> int:
        return self._C

    @property
    def discriminant(self) -> int:
        return self.B * self.B - 4 * self.A * self.C

    @property
    def coords(self) -> tuple[int, int, int]:
        return self.A, self.B, self.C

    def shift(self, a: int) -> "BinaryQuadraticForm":
        """Translate the form by `x -> x + a y`.
        For a monic form `(1, t, n)`, this is the norm form of the shifted
        element `pi + a`.
        """
        return type(self)(self.A, self.B + 2 * self.A * a, self.A * a * a + self.B * a + self.C)

    def to_D0_basis(self) -> "BinaryQuadraticForm":
        if self.C == 0:
            return (0, 0)
        h = self.B // 2
        trace_zero_form = self.shift(-h) # we move to (1,{0,1}, N')
        even_parity = trace_zero_form.B % 2 == 0
        _D0 = 1
        _c = 1
        # reduced discr, remove 4 factor, so either odd or we add abck 4 later
        disc = trace_zero_form.C if even_parity else (1 - 4 * trace_zero_form.C)
        for l, a in factorize(abs(disc)):
            c_exp = a // 2
            d_exp = a % 2
            _c *= l ** c_exp
            _D0 *= l ** d_exp
        if even_parity:
            _D0 *= 4
        return (-_D0, _c)
