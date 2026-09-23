from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from ..data import EigenForm

if TYPE_CHECKING:
    from ...arithmetic.forms import BinaryQuadraticForm


@dataclass
class FrobData:
    """One trace stratum containing its global eigenforms."""

    trace: int
    base_form: "BinaryQuadraticForm | None" = None
    eigen_forms: list[EigenForm] = field(default_factory=list)

    def survives(self) -> bool:
        return bool(self.eigen_forms)
