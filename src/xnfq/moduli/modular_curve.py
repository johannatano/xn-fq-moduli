from __future__ import annotations

from fractions import Fraction
from typing import Callable
from .data import FiberRecord
from .level_structures import Gamma, Gamma0, Gamma1, LevelStructure


class ModularCurve:
    """Base moduli problem at level `N`, independent of a finite field."""
    def __init__(self, level_structure: LevelStructure):
        self.level_structure = level_structure
        self.N = level_structure.N
        # By default do not emit progress logging. Callers may enable
        # console progress with `enable_console_progress(True)` or register
        # a custom handler via `set_progress_handler`.
        self._progress_handler = None

    class _NoopProgress:
        def __init__(self, total=None, desc=None):
            pass

        def __enter__(self):
            return self

        def update(self, n=1):
            return None

        def close(self):
            return None

        def __exit__(self, exc_type, exc, tb):
            return False

        @staticmethod
        def write(msg):
            return None

    class _ConsoleProgress:
        def __init__(self, total=None, desc=None):
            self.count = 0
            self.desc = desc or ""
            # prefer explicit total, otherwise use owner's q if present
            self.total = total
            self._closed = False
            self._last_len = 0

        def __enter__(self):
            if self.total:
                msg = f"{self.desc} 0.0% (0/{self.total})"
            else:
                msg = f"{self.desc} 0 (0/?)"
            print(msg, end="\r", flush=True)
            self._last_len = len(msg)
            return self

        def update(self, n=1):
            self.count += n
            if self.total:
                pct = 100.0 * self.count / float(self.total)
                msg = f"{self.desc} {pct:5.1f}% ({self.count}/{self.total})"
            else:
                msg = f"{self.desc} {self.count} items"
            pad = max(0, self._last_len - len(msg))
            print(msg + (" " * pad), end="\r", flush=True)
            self._last_len = len(msg)

        def close(self):
            if not self._closed:
                if self.total:
                    pct = 100.0 * self.count / float(self.total)
                    final = f"{self.desc} {pct:5.1f}% ({self.count}/{self.total})"
                else:
                    final = f"{self.desc} {self.count} items"
                print(final + (" " * max(0, self._last_len - len(final))))
                self._closed = True

        def __exit__(self, exc_type, exc, tb):
            self.close()
            return False

        @staticmethod
        def write(msg):
            print(msg)

    def set_progress_handler(self, handler: Callable[[int | None, str | None], object] | None) -> None:
        """Register a progress handler callable. Handler will be called as
        `handler(total=..., desc=...)` and should return an object with
        `update(n)`, `close()` and context-manager support. Pass `None` to
        disable progress reporting.
        """
        self._progress_handler = handler

    def _progress_ctx(self, total: int | None = None, desc: str | None = None):
        """Return a context manager progress object from the registered handler
        or a no-op progress object when none is registered."""
        if self._progress_handler is None:
            # Default: no-op progress (suppress console output). Callers may
            # enable console progress explicitly via `enable_console_progress`.
            return ModularCurve._NoopProgress()
        try:
            return self._progress_handler(total=total, desc=desc)
        except TypeError:
            # fallback: call with positional args
            return self._progress_handler(total, desc)

    def enable_console_progress(self, enable: bool = True) -> None:
        """Enable or disable the built-in console progress printer.

        When enabled, the curve will print progress messages to stdout for
        operations that use `_progress_ctx`. When disabled, progress is
        suppressed. Clients can still register a custom handler via
        `set_progress_handler`.
        """
        if enable:
            self._progress_handler = lambda total=None, desc=None: ModularCurve._ConsoleProgress(total=total, desc=desc)
        else:
            self._progress_handler = None

    def over(self, p: int, n: int) -> "ModularCurveFq":
        from .fq.curve import ModularCurveFq

        return ModularCurveFq(self.level_structure, p, n)

    def _decompose_q(self, q: int) -> tuple[int, int]:
        """Decompose q into p**n. If q is not a pure power, return (q,1).
        Try integer roots for n from high to low to find exact decomposition.
        """
        if q <= 1:
            return (q, 1)
        import math
        # Try exponents from high to low (max 64 to be safe)
        max_exp = int(math.log(q, 2)) + 1
        for n in range(max_exp, 1, -1):
            p = int(round(q ** (1.0 / n)))
            if p > 1 and pow(p, n) == q:
                return (p, n)
        return (q, 1)

    def F(self, q: int) -> "ModularCurveFq":
        """Convenience wrapper: set base field by full size `q = p**n`.

        Examples:
            `curve.F(5**4)` is equivalent to `curve.over(5,4)`.
        """
        p, n = self._decompose_q(q)
        return self.over(p, n)

    def change_base(self, p: int, n: int) -> "ModularCurveFq":
        return self.over(p, n)

    def get_structure(self) -> list[FiberRecord]:
        return [fiber.snapshot() for fiber in self.fibers()]

    def info(self) -> str:
        return f"X{self.level_structure.type}({self.N})"

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

""" Convenience wrappers"""
class X0(ModularCurve):
    def __init__(self, N: int):
        super().__init__(Gamma0(N))

class X1(ModularCurve):
    def __init__(self, N: int):
        super().__init__(Gamma1(N))

class X(ModularCurve):
    def __init__(self, N: int):
        super().__init__(Gamma(N))


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