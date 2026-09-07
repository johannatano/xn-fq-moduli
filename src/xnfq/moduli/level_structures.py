class LevelStructure:
    """Abstract level-structure descriptor used by the counting code."""

    def __init__(self, N: int):
        self.N = N
        self._class = None
        self._scalar_only = False

    @property
    def type(self) -> int:
        return self._class

    def eigenvalues(self, m: int) -> list[int]:
        """Return admissible eigenvalue shifts modulo `m`."""
        return list(range(0, m))

    @property
    def scalar_only(self) -> bool:
        return self._scalar_only


class Gamma0(LevelStructure):
    """Level structure for cyclic subgroups of order `N`."""

    def __init__(self, N: int):
        super().__init__(N)
        self._class = 0


class Gamma1(LevelStructure):
    """Level structure for a chosen point of exact order `N`."""

    def __init__(self, N: int):
        super().__init__(N)
        self._class = 1

    def eigenvalues(self, m: int) -> list[int]:
        return [1]

class Gamma(Gamma1):
    """Full level structure, restricted to scalar-compatible eigenvalues."""

    def __init__(self, N: int):
        super().__init__(N)
        self._class = 2
        self._scalar_only = True
