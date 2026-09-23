from __future__ import annotations

from fractions import Fraction
from math import isqrt
import time

from sympy import nextprime, prevprime, primerange

from xnfq.config import apply_config_from_args
from xnfq.moduli.modular_curve import X, X0, X1
from utils.args import parse_example_args
from utils.ui.dashboard import Dashboard, Param

MAX_PRIME_START = 100_000
MAX_WINDOW_SIZE = 50
MAX_TRACE_WORK = 1_000


def XN_over(level_type: int, N: int, p: int):
    if level_type == 0:
        return X0(N).over(p, 1)
    if level_type == 1:
        return X1(N).over(p, 1)
    if level_type == 2:
        return X(N).over(p, 1)
    raise ValueError(f"Unsupported type: {level_type}")


def curve_label(level_type: int) -> str:
    return {0: "X_0", 1: "X_1", 2: "X"}[level_type]


class HeckeTracePlotView(Dashboard):
    """Plot the symmetric-power Frobenius trace as the prime varies."""

    title = "Hecke traces over prime fields"
    figsize = (11.0, 6.5)

    def params(self):
        return [
            Param("N", 1, 100, 11, step=1, label="level N"),
            Param("k", 0, 20, 2, step=1, label="symmetric power k"),
            Param("type", 0, 2, 1, step=1, label="curve type"),
            Param(
                "prime_start",
                2,
                MAX_PRIME_START,
                5,
                step=1,
                label="starting prime",
            ),
            Param(
                "window_size",
                1,
                MAX_WINDOW_SIZE,
                25,
                step=1,
                label="primes per window",
            ),
        ]

    def panels(self):
        return ["trace"]

    def actions(self):
        return {
            "Previous prime window": self.previous_window,
            "Next prime window": self.next_window,
        }

    @staticmethod
    def prime_window(start: int, size: int) -> list[int]:
        """Return ``size`` consecutive primes starting at or above ``start``."""
        primes: list[int] = []
        candidate = max(2, start)
        upper = max(candidate + 2, candidate * 2)
        while len(primes) < size:
            primes.extend(primerange(candidate, upper))
            candidate = upper
            upper *= 2
        return primes[:size]

    def previous_window(self) -> None:
        primes = self.prime_window(self["prime_start"], self["window_size"])
        start = primes[0]
        for _ in range(self["window_size"]):
            if start <= 2:
                break
            start = prevprime(start)
        self.set("prime_start", start, redraw=False)

    def next_window(self) -> None:
        primes = self.prime_window(self["prime_start"], self["window_size"])
        self.set("prime_start", nextprime(primes[-1]), redraw=False)

    def _snapshot(self):
        key = (
            self["N"],
            self["k"],
            self["type"],
            self["prime_start"],
            self["window_size"],
        )
        if getattr(self, "_cache_key", None) == key:
            return self._cache_data

        N, k, level_type, prime_start, window_size = key
        estimated_work = window_size * max(1, isqrt(prime_start)) * max(1, N // 10)
        if estimated_work > MAX_TRACE_WORK:
            message = (
                "Window skipped: reduce starting prime, window size, or level "
                f"(estimated work {estimated_work:,}, limit {MAX_TRACE_WORK:,})."
            )
            self._cache_key = key
            self._cache_data = {
                "p": prime_start,
                "N": N,
                "k": k,
                "type": level_type,
                "prime_start": prime_start,
                "window_size": window_size,
                "traces": [],
                "elapsed": 0.0,
                "error": message,
            }
            return self._cache_data

        primes = self.prime_window(prime_start, window_size)
        started = time.perf_counter()
        traces: list[tuple[int, int | Fraction]] = []
        
        for p in primes:
            curve = XN_over(level_type, N, p)
            # we use k+2, k-1+2=k+1
            norm = curve.tr_frob_symk(k) / p**((k+1)//2)
            traces.append((p, norm))

        self._cache_key = key
        self._cache_data = {
            "N": N,
            "k": k,
            "type": level_type,
            "prime_start": prime_start,
            "window_size": window_size,
            "traces": traces,
            "elapsed": time.perf_counter() - started,
            "error": None,
        }
        return self._cache_data

    def draw(self, name: str, ax) -> None:
        snapshot = self._snapshot()
        N = snapshot["N"]
        k = snapshot["k"]
        level_type = snapshot["type"]
        prime_start = snapshot["prime_start"]
        window_size = snapshot["window_size"]
        traces = snapshot["traces"]
        label = curve_label(level_type)
        error = snapshot.get("error")

        ax.set_title(f"Tr Sym^{k}(Frob_p) for {label}({N})")
        if error:
            ax.text(
                0.5,
                0.5,
                error,
                ha="center",
                va="center",
                wrap=True,
                transform=ax.transAxes,
            )
            ax.set_axis_off()
            return

        if traces:
            primes = [p for p, _ in traces]
            values = [float(trace) for _, trace in traces]
            ax.plot(primes, values, marker="o", markersize=4, linewidth=1.2, color="steelblue")
            ax.set_xlim(primes[0], primes[-1])
            y_min = min(values)
            y_max = max(values)
            if y_min == y_max:
                padding = max(1.0, abs(y_min) * 0.08)
            else:
                padding = (y_max - y_min) * 0.08
            ax.set_ylim(y_min - padding, y_max + padding)
        else:
            ax.text(0.5, 0.5, "no primes in selected range", ha="center", va="center")

        ax.set_xlabel("q = p")
        ax.set_ylabel("final trace value")
        ax.grid(color="0.92", linewidth=0.8)
        ax.text(
            0.98,
            0.98,
            (
                f"primes={len(traces)}\n"
                f"start={prime_start}\n"
                f"window={window_size} primes\n"
                f"range=[{traces[0][0]}, {traces[-1][0]}]\n"
                f"time={snapshot['elapsed']:.3f}s"
            ),
            transform=ax.transAxes,
            ha="right",
            va="top",
            bbox={"facecolor": "white", "edgecolor": "0.8", "boxstyle": "round,pad=0.4"},
        )


def run() -> None:
    args = parse_args()
    apply_config_from_args(args)
    HeckeTracePlotView(
        N=args.N,
        k=args.k,
        type=args.type,
        prime_start=args.p,
    ).show()


def parse_args():
    return parse_example_args(
        "Plot Tr Sym^k(Frob_p) over a prime range",
        include_sym_power=True,
    )


if __name__ == "__main__":
    run()
