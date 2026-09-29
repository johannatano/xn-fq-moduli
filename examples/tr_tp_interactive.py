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
    title = ""
    figsize = (8.0, 4.75)
    use_text_inputs = True
    prime_step = 5000
    prime_window_size = 5000

    def __init__(self, *, N: int = 100, k: int = 2, type: int = 1, pmin: int = 2, pmax: int | None = None):
        self.initial_N = int(N)
        self.initial_k = int(k)
        self.initial_type = int(type)
        self.pmin = max(2, int(pmin))
        self.pmax = max(
            self.pmin,
            int(pmax if pmax is not None else self.pmin + self.prime_window_size - 1),
        )
        super().__init__()

    def params(self):
        return [
            Param("N", 1, 10_000, self.initial_N, step=1, label="level"),
            Param("k", 2, 10_000, self.initial_k, step=1, label="weight"),
            Param("type", 0, 2, self.initial_type, step=1, label=r"$\Gamma$"),
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
            f"Tr(T_p, S[Gamma{level_type}({N}), k={k}]) (p range {pmin}-{pmax})",
            Colors.NEON_PURPLE,
        )

        for p in primes:
            curve = XN_over(level_type, N, p)
            val = curve.tr_fq(k).val
            norm = val / p ** ((k - 1) / 2)
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

        ax.set_title(
            rf"$\mathrm{{Tr}}(T_p \mid S_{{{k}}}(\Gamma_{{{snapshot['type']}}}({N})))$"
        )

        if traces:
            primes = [p for p, _ in traces]
            values = [float(trace) for _, trace in traces]
            ax.scatter(primes, values, s=4, color="black", marker="o")  # scatter
            # ax.plot(primes, values, markersize=1, linewidth=1.0, color="black")  # line
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
        ax.set_ylabel("$a_p$ normalized")
        ax.grid(color="grey", linewidth=.5)


def run() -> None:
    args = parse_example_args(
        "",
        include_weight=True,
        overrides={
            "p": 50000,
            "k": 4,
            "N": 1000,
        },
    )
    apply_config_from_args(args)
    HeckeTraceInteractive(
        N=args.N,
        k=args.k,
        type=args.type,
        pmin=args.p,
        pmax=args.p + HeckeTraceInteractive.prime_window_size - 1,
    ).show()


if __name__ == "__main__":
    run()
