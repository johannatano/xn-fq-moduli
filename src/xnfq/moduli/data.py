from dataclasses import dataclass, field
from fractions import Fraction
@dataclass(frozen=True)
class FiberRecord:
    kind: str
    count: Fraction = Fraction(0)
