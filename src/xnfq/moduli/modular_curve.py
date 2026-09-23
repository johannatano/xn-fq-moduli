from __future__ import annotations

from fractions import Fraction
from .data import FiberRecord, StructureReport
from .level_structures import Gamma, Gamma0, Gamma1, LevelStructure

class ModularCurve:
    """Base moduli problem at level `N`, independent of a finite field."""

    def __init__(self, level_structure: LevelStructure):
        self.level_structure = level_structure
        self.N = level_structure.N

    def over(self, p: int, n: int) -> "ModularCurveFq":
        from .fq.curve import ModularCurveFq

        return ModularCurveFq(self.level_structure, p, n)

    def change_base(self, p: int, n: int) -> "ModularCurveFq":
        return self.over(p, n)

    def get_structure(self) -> StructureReport:
        """Materialize the current fiber decomposition into serializable records."""
        fibers = tuple(self.fibers())
        return StructureReport(
            problem=self.level_structure,
            fibers=tuple(fiber.snapshot() for fiber in fibers),
        )

    def fibers(self):
        yield from self.smooth_fibers()
        yield from self.cusps()

    def smooth_fibers(self):
        pass

    def cusps(self):
        pass

    def count(self):
        """Return the smooth and cusp counts separately."""
        return (sum(f.count() for f in self.smooth_fibers()), sum(f.count() for f in self.cusps()))


class CurveFiber:
    """Common interface for one fiber contribution to the point count."""

    def __init__(
        self,
        gamma: "LevelStructure",
    ):
        self.gamma = gamma
        self.N = gamma.N

    def snapshot(self) -> FiberRecord:
        raise NotImplementedError
    def count(self) -> Fraction:
        raise NotImplementedError

class X0(ModularCurve):
    """Convenience wrapper for the `Gamma_0(N)` moduli problem."""

    def __init__(self, N: int):
        super().__init__(Gamma0(N))

class X1(ModularCurve):
    """Convenience wrapper for the `Gamma_1(N)` moduli problem."""

    def __init__(self, N: int):
        super().__init__(Gamma1(N))

class X(ModularCurve):
    """Convenience wrapper for the full level-`N` moduli problem."""

    def __init__(self, N: int):
        super().__init__(Gamma(N))
