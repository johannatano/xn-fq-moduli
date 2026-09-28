from ..arithmetic.function import phi

class LevelStructure:
    """Abstract level-structure descriptor used by the counting code."""
    def __init__(self, N: int):
        self.N = N
        self._class = None
        self._scalar_only = False
    @property
    def type(self) -> int:
        return self._class
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
    
    def weight(self, smooth: bool = True) -> int:
        return phi(1)(self.N)

class Gamma2(Gamma1):
    """Full level structure (internal name: Gamma2)."""
    def __init__(self, N: int):
        super().__init__(N)
        self._class = 2
        self._scalar_only = True

    def weight(self, smooth: bool = True) -> int:
        if not smooth:
            return phi(-1)(self.N) * phi(1)(self.N)
        return self.N * phi(1)(self.N)


def Gamma(N: int, type: int = 2) -> LevelStructure:
    """Factory returning a LevelStructure for given `N` and `type`.
    - `type=0` -> `Gamma0(N)`
    - `type=1` -> `Gamma1(N)`
    - `type=2` -> full level `Gamma2(N)` (default)
    """
    if type == 0:
        return Gamma0(N)
    if type == 1:
        return Gamma1(N)
    return Gamma2(N)
