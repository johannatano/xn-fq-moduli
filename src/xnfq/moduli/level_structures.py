from ..arithmetic.function import phi

class LevelStructure:
    """Abstract level-structure descriptor used by the counting code."""

    def __init__(self, N: int):
        from ..arithmetic.function import phi
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
    
    def weight(self, smooth: bool = True) -> int:
        return 1


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
    
    def weight(self, smooth: bool = True) -> int:
        return phi(1)(self.N)


class Gamma(Gamma1):
    """Full level structure, restricted to scalar-compatible eigenvalues."""
    def __init__(self, N: int):
        super().__init__(N)
        self._class = 2
        self._scalar_only = True
        
    def weight(self, smooth:bool = True) -> int:
        if not smooth:
            return phi(-1)(self.N) * phi(1)(self.N)
        return self.N*phi(1)(self.N)
