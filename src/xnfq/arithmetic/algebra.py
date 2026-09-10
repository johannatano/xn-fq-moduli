from __future__ import annotations

from dataclasses import dataclass
from math import prod

from .common import factorize, valuation as vl


class Set:
    elements: type


class AdditiveStructure(Set):
    def add(self, x, y):
        pass

    def negate(self, x):
        pass


class Ring(AdditiveStructure):
    def multiply(self, x, y):
        pass


class IntegerRing(Ring):
    pass


class Module:
    ring: Ring
    group: AdditiveStructure

    def act(self, r, m):
        pass


class ZModule(Module):
    ring = IntegerRing()

    def act(self, z, m):
        if isinstance(m, tuple):
            return tuple(z * coord for coord in m)
        return z * m


@dataclass(frozen=True)
class AbGrp(Set):
    """A finite abelian group presented by invariant factors."""

    invariants: tuple[int, ...]
    generator: object | None = None

    @property
    def size(self) -> int:
        return prod(self.invariants)

    def p_primary_part(self, p: int) -> "PSylow":
        factors = factorize(p)
        if len(factors) != 1 or factors[0][0] != p or factors[0][1] != 1:
            raise ValueError("p must be prime")
        invariants = tuple(p**vp for d in self.invariants for vp in [vl(d, p)] if vp > 0)
        return PSylow(p, invariants, generator=self.generator)

    def sylow(self, p: int) -> "PSylow":
        return self.p_primary_part(p)


@dataclass(frozen=True, init=False)
class PSylow(AbGrp):
    """The maximal p-primary part of a finite abelian group."""

    prime: int

    def __init__(self, prime: int, invariants: tuple[int, ...], generator: object | None = None):
        factors = factorize(prime)
        #if len(factors) != 1 or factors[0][0] != prime or factors[0][1] != 1:
        #    raise ValueError("prime must be prime")
        # any(d <= 0 or vl(d, prime) == 0 or d != prime ** vl(d, prime) for d in invariants):
        #    raise ValueError("PSylow invariants must be positive powers of the given prime")
        object.__setattr__(self, "prime", prime)
        object.__setattr__(self, "invariants", invariants)
        object.__setattr__(self, "generator", generator)

    @classmethod
    def from_group(cls, group: AbGrp, prime: int) -> "PSylow":
        return group.p_primary_part(prime)


class Lattice(ZModule):
    ambient_space: Module | None = None
    basis: list

    def __init__(self, ambient_space: Module | None = None, basis: list | None = None):
        self.ambient_space = ambient_space
        self.basis = [] if basis is None else list(basis)

    def quotient(self, matrix):
        raise NotImplementedError


class Field(Ring):
    def divide(self, x, y):
        pass


class RationalField(Field):
    pass


class FieldExtension(Field, Module):
    def __init__(self, base_field: Field, degree: int, basis: list):
        self.ring = base_field
        self.degree = degree
        self.basis = basis

    def act(self, scalar, element):
        pass


class NumberField(FieldExtension):
    def __init__(self, degree: int, basis: list):
        super().__init__(base_field=RationalField(), degree=degree, basis=basis)

class NeronDgon(AbGrp):
    def __init__(self, d: int):
        self.d = d

    def torsion_subgroup(self, n: int) -> "NeronDgon":
        """Return the n-torsion subgroup of the NeronDgon."""
        if n % self.d != 0:
            return (0,0)
        return (n // self.d, self.d)
