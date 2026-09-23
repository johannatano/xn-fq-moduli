from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction

from ...arithmetic.quadratic import QuadraticOrder
from ...arithmetic.forms import BinaryQuadraticForm

from ..data import FiberRecord

class CuspForm:
    def __init__(self, d1: int, d2:int):
        self.d1 = d1
        self.d2 = d2

@dataclass(frozen=True)
class LevelStructureRecord:
    index: int = 0
    mass: int = 0
    coords: tuple[int, int] = field(default_factory=tuple)
    inv: tuple[int, int] = field(default_factory=tuple)
    num_lines: Fraction = Fraction(0)
    scalar: bool = False
    order: QuadraticOrder = None

@dataclass
class EigenForm:
    value: int
    form: "BinaryQuadraticForm | CuspForm"

@dataclass(frozen=True)
class EigenFormRecord:
    eigenform: EigenForm
    levels: tuple[LevelStructureRecord, ...] = ()

@dataclass(frozen=True)
class FiberRecordFq(FiberRecord):
    t: int = 0
    m0: Fraction = Fraction(0)
    eigen_records: tuple[EigenFormRecord, ...] = ()


@dataclass(frozen=True)
class SmoothFiberRecordFq(FiberRecordFq):
    """Serializable summary of a smooth fiber over Fq."""
    kind: str = "smooth"
    D_K: int = 0

@dataclass(frozen=True)
class CuspFiberRecordFq(FiberRecordFq):
    """Serializable summary of a cusp fiber over Fq."""
    kind: str = "cusp"
