from __future__ import annotations

import sys
import time

from sympy import primerange

from xnfq.config import apply_config_from_args
from xnfq.moduli import ModularCurve, Gamma
from utils.args import parse_example_args
from utils.logging import Logger, Colors

from utils.ui.dashboard import Dashboard, Param


def XN_over(level_type: int, N: int, p: int):
    return ModularCurve(Gamma(N, level_type)).over(p, 1)


class HeckeTraceInteractive(Dashboard):
    """Interactive dashboard for Tr(F_p) (symmetric powers) as p varies.

    This view mirrors the previous `tr_tp_range` dashboard but uses the
    `tr_fq` computation for the curves over F_p (n=1)."""

    title = "Hecke Trace (interactive)"
    figsize = (11.0, 6.5)
    use_text_inputs = True
    prime_step = 1000
    prime_window_size = 1000

    def __init__(self, **overrides):
        self.pmin = max(2, int(overrides.pop("pmin", 2)))
        self.pmax = max(
            self.pmin,
            int(overrides.pop("pmax", self.pmin + self.prime_window_size - 1)),
        )
        super().__init__(**overrides)

    def params(self):
        return [
            Param("N", 1, 10_000, 100, step=1, label="N"),
            Param("k", 0, 10_000, 2, step=1, label="k"),
            Param("type", 0, 2, 1, step=1, label="Gamma"),
        ]

    def panels(self):
        return ["trace"]

    def actions(self):
        return {"<< prev": self.previous_window, "next >>": self.next_window}

    def primes(self) -> list[int]:
        return list(primerange(self.pmin, self.pmax + 1))

    def previous_window(self) -> None:
        self.shift_window(-1)

    def next_window(self) -> None:
        self.shift_window(1)

    def shift_window(self, direction: int) -> None:
        new_pmin = max(2, self.pmin + direction * self.prime_step)
        actual_shift = new_pmin - self.pmin
        self.pmin = new_pmin
        self.pmax += actual_shift
        self._cache_key = None

    def _snapshot(self):
        key = (self["N"], self["k"], self["type"], self.pmin, self.pmax)
        if getattr(self, "_cache_key", None) == key:
            return self._cache_data

        N, k, level_type, pmin, pmax = key
        primes = list(primerange(pmin, pmax + 1))
        started = time.perf_counter()
        traces: list[tuple[int, float]] = []

        Logger.cprint(
            f"Tr(T_p, S[Gamma{level_type}({N}), k={k})) (p range {pmin}-{pmax})",
            Colors.NEON_PURPLE,
        )

        for p in primes:
            curve = XN_over(level_type, N, p)
            norm = curve.tr_fq(k) / p ** ((k - 1) / 2)
            traces.append((p, norm))

        self._cache_key = key
        self._cache_data = {
            "N": N,
            "k": k,
            "type": level_type,
            "pmin": pmin,
            "pmax": pmax,
            "traces": traces,
            "elapsed": time.perf_counter() - started,
        }
        return self._cache_data

    def draw(self, name: str, ax) -> None:
        snapshot = self._snapshot()
        N = snapshot["N"]
        k = snapshot["k"]
        pmin = snapshot["pmin"]
        pmax = snapshot["pmax"]
        traces = snapshot["traces"]

        ax.set_title(f"Tr(T_p | S_{k}({N}))")

        if traces:
            primes = [p for p, _ in traces]
            values = [float(trace) for _, trace in traces]
            ax.plot(primes, values, markersize=1, linewidth=1.0, color="black")
            ax.set_xlim(pmin, pmax)
            y_min = min(values)
            y_max = max(values)
            if y_min == y_max:
                padding = max(1.0, abs(y_min) * 0.08)
            else:
                padding = (y_max - y_min) * 0.08
            ax.set_ylim(y_min - padding, y_max + padding)
        else:
            ax.text(0.5, 0.5, "no primes in selected range", ha="center", va="center")

        ax.set_xlabel("p")
        ax.tick_params(axis="y", labelleft=False)
        ax.set_ylabel("a_p normalized")
        ax.grid(color="grey", linewidth=.5)


def run() -> None:
    args = parse_example_args("Hecke Trace (interactive)", include_weight=True)
    apply_config_from_args(args)
    overrides = {"pmin": args.p, "pmax": args.p + HeckeTraceInteractive.prime_window_size - 1}
    for name in ("N", "k", "type"):
        value = getattr(args, name)
        if value is not None:
            overrides[name] = value
    HeckeTraceInteractive(**overrides).show()


if __name__ == "__main__":
    run()
