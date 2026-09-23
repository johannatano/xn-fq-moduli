from __future__ import annotations

import sys
import time

from sympy import primerange

from xnfq.config import apply_config_from_args
from xnfq.moduli.modular_curve import X, X0, X1
from utils.args import parse_example_args
from utils.ui.dashboard import Dashboard, Param

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

    title = "Hecke Trace"
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
        return {
            "<< prev": self.previous_window,
            "next >>": self.next_window,
        }

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
        key = (
            self["N"],
            self["k"],
            self["type"],
            self.pmin,
            self.pmax,
        )
        if getattr(self, "_cache_key", None) == key:
            return self._cache_data

        N, k, level_type, pmin, pmax = key
        primes = list(primerange(pmin, pmax + 1))
        started = time.perf_counter()
        traces: list[tuple[int, float]] = []

        for p in primes:
            curve = XN_over(level_type, N, p)
            norm = curve.tr_frob_symk(k) / p ** ((k + 1) / 2)
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
        level_type = snapshot["type"]
        pmin = snapshot["pmin"]
        pmax = snapshot["pmax"]
        traces = snapshot["traces"]
        label = curve_label(level_type)

        ax.set_title(f"Tr(T_p | S_{k+2}({N}))")

        if traces:
            primes = [p for p, _ in traces]
            values = [float(trace) for _, trace in traces]
            ax.plot(
                primes,
                values,
                markersize=1,
                linewidth=1.0,
                color="black",
            )
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
    args = parse_args()
    apply_config_from_args(args)
    overrides = {
        "pmin": args.p,
        "pmax": args.p + HeckeTracePlotView.prime_window_size - 1,
    }
    for name in ("N", "k", "type"):
        value = getattr(args, name)
        if value is not None:
            overrides[name] = value
    HeckeTracePlotView(**overrides).show()


def parse_args():
    args = parse_example_args(
        "Hecke Trace",
        include_sym_power=True,
    )
    options = set(sys.argv[1:])
    if not {"-N", "--N"} & options:
        args.N = None
    if not {"-type", "--type"} & options:
        args.type = None
    if "-k" not in options:
        args.k = None
    return args


if __name__ == "__main__":
    run()
