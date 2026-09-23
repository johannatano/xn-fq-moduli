from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from ..arithmetic.quadratic import QuadraticOrder

if TYPE_CHECKING:
    from ..arithmetic.forms import BinaryQuadraticForm
    from ..arithmetic.algebra import NeronDgon


@dataclass
class EigenForm:
    """One global eigenvalue together with its model form."""
    value: int
    form: "BinaryQuadraticForm | NeronDgon"

@dataclass(frozen=True)
class LevelStructureRecord:
    """One local group-structure contribution inside a smooth fiber."""
    index: int
    mass: int
    coords: tuple[int, int]
    coords_pi_basis: tuple[int, int] = field(default_factory=tuple)
    im: tuple[int, int] = field(default_factory=tuple)
    ker: tuple[int, int] = field(default_factory=tuple)
    num_lines: int = 0
    scalar: bool = False
    order: QuadraticOrder = None

@dataclass(frozen=True)
class EigenFormRecord:
    """A global eigenvalue together with its local structure across all primes."""
    eigenform: EigenForm
    level_records: tuple[LevelStructureRecord, ...] = ()

@dataclass(frozen=True)
class FiberRecord:
    """Serializable summary of either a Weil fiber or a cusp fiber."""

    kind: str
    t: int
    d: int | None
    total_count: Fraction = Fraction(0)
    normalized_count: Fraction = Fraction(0)
    total_mass: Fraction = Fraction(0)
    conductors: tuple[int, ...] = ()
    coprime: tuple[int, ...] = ()
    eigen_records: tuple[EigenFormRecord, ...] = ()
    eigenvalues: tuple[int, ...] = ()
    lattice_levels: tuple[int, ...] = ()
    lattice_counts: tuple[tuple[int, Fraction], ...] = ()

    @property
    def trace(self) -> int:
        return self.t

    @property
    def discriminant(self) -> int | None:
        return self.d

@dataclass(frozen=True)
class StructureReport:
    """Top-level structure summary for a fixed level problem."""

    problem: object
    fibers: tuple[FiberRecord, ...]

    @property
    def total(self) -> Fraction:
        return sum((fiber.total_count for fiber in self.fibers), start=Fraction(0))
