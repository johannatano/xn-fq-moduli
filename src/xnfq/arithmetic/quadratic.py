from __future__ import annotations
from dataclasses import dataclass, field
from fractions import Fraction
from math import gcd, isqrt, prod

from .algebra import AbGrp, Lattice, NumberField, PSylow
from .forms import BinaryQuadraticForm
from .common import factorize, legendre, valuation as vl
from .function import KroneckerCharacter, phi
from ..config import get_pari


class QuadraticNumberField(NumberField):
    def __init__(self, d: int):
        self.d = d
        self.is_1_mod_4 = d % 4 == 1
        self.d_K = d if self.is_1_mod_4 else 4 * d
        super().__init__(degree=2, basis=[(1, 0), (0, 1)])

    @classmethod
    def from_discriminant(cls, DK: int) -> "QuadraticNumberField":
        d = DK if DK % 4 == 1 else DK // 4
        return cls(d)


@dataclass(frozen=True, init=False)
class BicyclicPSylow(PSylow):
    """A bicyclic p-primary group recorded by exponent data."""

    exponents: tuple[int, int]

    def __init__(self, prime: int, e1: int, e2: int, generator: object | None = None):
        if e1 < 0 or e2 < 0:
            raise ValueError("BicyclicPSylow exponents must be non-negative")
        if e1 > e2:
            raise ValueError("BicyclicPSylow exponents must satisfy e1 <= e2")
        invariants = tuple(prime**e for e in (e1, e2))
        super().__init__(prime, invariants, generator=generator)
        object.__setattr__(self, "exponents", (e1, e2))

    @property
    def e1(self) -> int:
        return self.exponents[0]

    @property
    def e2(self) -> int:
        return self.exponents[1]

    @property
    def d1(self) -> int:
        return self.prime**self.e1

    @property
    def d2(self) -> int:
        return self.prime**self.e2

    def C(self, a: int) -> int:
        if self.e2 < a:
            return 0
        if self.e1 < a:
            return self.prime**self.e1
        return phi(-1)(self.prime**a)

    def C_full(self, a: int) -> tuple[int, bool]:
        if self.e2 < a:
            return (0, False)
        if self.e1 < a:
            return (self.prime**self.e1, False)
        return (phi(-1)(self.prime**a), True)

    def P(self, a: int) -> int:
        s1 = min(a, self.e1)
        s2 = min(a, self.e2)
        return self.prime ** (s1 + s2) - self.prime ** (min(a - 1, s1) + min(a - 1, s2))


@dataclass(frozen=True, init=False)
class BicyclicGroup(AbGrp):
    """A finite abelian group presented by invariant factors."""

    def __init__(self, d1: int, d2: int, generator: object | None = None):
        object.__setattr__(self, "invariants", (d1, d2))
        object.__setattr__(self, "generator", generator)

    @property
    def d1(self) -> int:
        return self.invariants[0]

    @property
    def d2(self) -> int:
        return self.invariants[1]

    def p_primary_part(self, p: int) -> BicyclicPSylow:
        return BicyclicPSylow(p, vl(self.d1, p), vl(self.d2, p), generator=self.generator)

    def sylow(self, p: int) -> BicyclicPSylow:
        return self.p_primary_part(p)

    def ell_invariants(self, l: int) -> tuple[int, int]:
        return self.sylow(l).exponents

    def C_ell(self, l: int, a: int) -> int:
        return self.sylow(l).C(a)

    def P_ell(self, ell: int, a: int) -> int:
        return self.sylow(ell).P(a)

    def P_N(self, N: int) -> int:
        if self.d2 < N:
            return 0
        return prod(self.sylow(ell).P(a) for ell, a in factorize(N))

    def C_N(self, N: int) -> int:
        if self.d2 < N:
            return 0
        return prod(self.sylow(ell).C(a) for ell, a in factorize(N))
    
    def tor_subgroup(self, N: int) -> int:
        return (gcd(self.d1,N), gcd(self.d2,N))


class QuadraticLattice(Lattice):
    order: "QuadraticOrder"

    def __init__(self, order: "QuadraticOrder"):
        self.order = order
        super().__init__(ambient_space=order.ambient_space, basis=[(1, 0), (0, order.f)])

    def quotient_by(self, element: "QuadraticOrderElement") -> BicyclicGroup:
        if element.order is not self.order:
            raise ValueError("element does not belong to this lattice order")
        d1 = gcd(element.u, element.v)
        d2 = element.norm // d1
        return BicyclicGroup(d1, d2, generator=element)

    def quotient_by_sylow(
        self, element: "QuadraticOrderElement", l: int, e: int = 0
    ) -> BicyclicPSylow:
        if element.order is not self.order:
            raise ValueError("element does not belong to this lattice order")
        if e < 0:
            raise ValueError("conductor exponent must be non-negative")
        v = vl(element.v, l)
        if v < e:
            raise ValueError("element does not belong to the requested l-power suborder")
        e1 = min(vl(element.u, l), v - e)
        e2 = vl(element.norm, l) - e1
        return BicyclicPSylow(l, e1, e2, generator=element)

class QuadraticOrder:
    """Quadratic order of conductor f in the basis [1, f*omega].

    Carries the norm-form data (1, B, C) with B^2 - 4C = f^2 * d_K.
    """

    def __init__(self, ambient_field: "QuadraticNumberField", conductor_f: int = 1):
        self.ambient_space = ambient_field
        self.f = conductor_f
        DK, f2 = ambient_field.d_K, conductor_f * conductor_f
        if DK % 4 == 0:
            self.B, self.C = 0, -(f2 * DK) // 4
        else:
            self.B, self.C = conductor_f, (f2 * (1 - DK)) // 4
        self.discriminant = f2 * DK
        self.lattice = QuadraticLattice(self)

    def __repr__(self) -> str:
        return f"QuadraticOrder(d_K={self.ambient_space.d_K}, f={self.f})"

    @property
    def form_coefficients(self) -> tuple[int, int, int]:
        return (1, self.B, self.C)

    def omega(self) -> "QuadraticOrderElement":
        return self.element(0, 1)

    def element(self, u: int, v: int) -> "QuadraticOrderElement":
        return QuadraticOrderElement(self, (u, v))

    def suborder(self, conductor_f: int) -> "QuadraticOrder":
        return QuadraticOrder(self.ambient_space, conductor_f)

    def embed_suborder(
        self, max_element: "QuadraticOrderElement | tuple[int, int]"
    ) -> "QuadraticOrderElement":
        if isinstance(max_element, QuadraticOrderElement):
            x, y = max_element.coords
        else:
            x, y = max_element
        if y % self.f:
            raise ValueError(
                f"Element does not belong to the order of conductor {self.f}"
            )
        return QuadraticOrderElement(self, (x, y // self.f))

    def embed_into_suborder(
        self, max_element: "QuadraticOrderElement | tuple[int, int]"
    ) -> "QuadraticOrderElement":
        return self.embed_suborder(max_element)

@dataclass(frozen=True)
class QuadraticOrderElement:
    """An element of a quadratic order written in the basis `[1, f * omega]`."""

    order: "QuadraticOrder"
    coords: tuple[int, int]
    u: int = field(init=False, repr=False, compare=False)
    v: int = field(init=False, repr=False, compare=False)
    norm: int = field(init=False, repr=False, compare=False)
    trace: int = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if len(self.coords) != 2:
            raise ValueError("quadratic order elements require two coordinates")
        u, v = self.coords
        set_ = object.__setattr__
        set_(self, "u", u)
        set_(self, "v", v)
        set_(self, "norm", u * u + self.order.B * u * v + self.order.C * v * v)
        set_(self, "trace", 2 * u + self.order.B * v)

    @classmethod
    def from_norm_form(
        cls,
        order: "QuadraticOrder",
        bqf: "BinaryQuadraticForm",
        shift: int = 0,
    ) -> "QuadraticOrderElement":
        if bqf.A != 1:
            raise ValueError("from_norm_form requires a monic form (1, t, n)")
        if order.discriminant == 0:
            return cls(order, (bqf.B // 2 + shift, 0))
        v = isqrt(bqf.discriminant // order.discriminant)
        return cls(order, ((bqf.B - order.B * v) // 2 + shift, v))

    def content(self) -> int:
        return gcd(self.u, self.v)

    def ell_content(self, l: int) -> int:
        return min(vl(self.u, l), vl(self.v, l))

    def ell_invariants(self, l: int) -> tuple[int, int]:
        d1 = self.ell_content(l)
        d2 = vl(self.norm, l) - d1
        return (d1, d2)

    def ell_kernel(self, l: int, a:int) -> tuple[int, int]:
        d1 = min(a, self.ell_content(l))
        d2 = min(a, vl(self.norm, l) - d1)
        return (d1, d2)

    @property
    def norm_form_coeffs(self) -> tuple[int, int, int]:
        return (1, self.trace, self.norm)

    def to_norm_form(self) -> BinaryQuadraticForm:
        return BinaryQuadraticForm(1, self.trace, self.norm)

    @property
    def maximal_order_coords(self) -> tuple[int, int]:
        """Coordinates of the same element in the maximal-order basis `[1, omega]`."""
        return (self.u, self.v * self.order.f)

    def coords_in_basis(
        self,
        generator: "QuadraticOrderElement",
    ) -> tuple[Fraction, Fraction]:
        """Write `self` in the basis `[1, generator]` inside the ambient field.

        Returns `(a, b)` with `self = a + b * generator`. The coefficients are
        generally rational unless `generator` itself gives an integral basis.
        """

        if self.order.ambient_space.d_K != generator.order.ambient_space.d_K:
            raise ValueError("basis change requires elements from the same ambient field")

        gen_u, gen_v = generator.maximal_order_coords
        if gen_v == 0:
            return self.coords
            raise ValueError("generator must be non-rational to define a quadratic basis")

        self_u, self_v = self.coords
        b = Fraction(self_v, gen_v)
        a = Fraction(self_u) - b * gen_u
        return (a, b)

    def basis_index(self) -> int:
        """Index of the lattice `Z + Z*self` inside the current order."""

        return abs(self.v)

    def shift(self, a: int) -> "QuadraticOrderElement":
        return type(self)(self.order, (self.u + a, self.v))


@dataclass(init=False)
class LatticeTower:
    DK: int
    f_max: int
    _min_form: "BinaryQuadraticForm"
    exclude: int | None = None
    include: int | None = None
    conductors: dict[int, int] = field(init=False)
    ambient_field: QuadraticNumberField = field(init=False, repr=False)
    _orders: dict[int, QuadraticOrder] = field(init=False, repr=False)
    OK: QuadraticOrder = field(init=False)
    chi_K: KroneckerCharacter = field(init=False, repr=False)

    def __init__(
        self,
        DK: int,
        f_max: int,
        min_form: "BinaryQuadraticForm" = None,
        exclude: int | None = None,
        include: int | None = None,
    ):
        self.DK = DK
        self.f_max = f_max
        self._min_form = min_form
        self.exclude = exclude
        self.include = include
        self.__post_init__()

    @classmethod
    def from_form(
        cls,
        min_form: "BinaryQuadraticForm",
        exclude: int | None = None,
        include: int | None = None,
        pari=None,
    ) -> "LatticeTower":
        from .forms import BinaryQuadraticForm
        DK, f = BinaryQuadraticForm._D0(min_form.discriminant, pari=pari)
        return cls(DK, f, min_form, exclude=exclude, include=include)

    def __post_init__(self) -> None:
        self.ambient_field = QuadraticNumberField.from_discriminant(self.DK)
        self._orders = {}
        self.OK = self.order(1)
        self.chi_K = KroneckerCharacter(self.DK)

        if self.f_max == 0:
            self.conductors = {1: 0}
            return

        def include_prime(p: int) -> bool:
            if self.exclude is not None and p == self.exclude:
                return False
            if self.include is not None:
                return p == self.include
            return True

        self.conductors = {p: a for p, a in factorize(self.f_max) if include_prime(p)}
        self.conductors[1] = 0

    @property
    def O_K(self) -> QuadraticOrder:
        return self.OK

    def order(self, conductor_f: int = 1) -> QuadraticOrder:
        order = self._orders.get(conductor_f)
        if order is None:
            order = QuadraticOrder(self.ambient_field, conductor_f)
            self._orders[conductor_f] = order
        return order

    def supported_divisors(self, support: int | None = None) -> list[int]:
        if support is None:
            return self.conductor_divisors()

        supported, _ = self.split(support)
        return self.divisors(supported)

    def split(self, N: int) -> tuple[dict[int, int], dict[int, int]]:
        coprime = {}
        supported = {}
        for l, k in self.conductors.items():
            if N % l == 0:
                supported[l] = k
            else:
                coprime[l] = k
        return supported, coprime

    @staticmethod
    def divisors(conductors: dict[int, int]) -> list[int]:
        divs = [1]
        for l, k in conductors.items():
            if l == 1:
                continue
            divs = [d * (l**i) for d in divs for i in range(k + 1)]
        return sorted(set(divs))

    def conductor_divisors(self) -> list[int]:
        return self.divisors(self.conductors)

    def max_coprime_conductor(self, N: int) -> int:
        _, coprime = self.split(N)
        return prod(l**a for l, a in coprime.items())

    def class_size(self, f: int) -> int:
        return phi(self.chi_K)(f)

    def min_form(self, sign: int) -> "BinaryQuadraticForm":
        from .forms import BinaryQuadraticForm
        if sign > 0:
            return self._min_form
        return BinaryQuadraticForm(1, -self._min_form.B, self._min_form.C)

    def quotient_group(self, element: QuadraticOrderElement, conductor_f: int) -> BicyclicGroup:
        order = self.order(conductor_f)
        return order.lattice.quotient_by(order.embed_suborder(element))

    def quotient_sylow(
        self, element: QuadraticOrderElement, l: int, conductor_exp: int
    ) -> BicyclicPSylow:
        return self.OK.lattice.quotient_by_sylow(element, l, e=conductor_exp)
