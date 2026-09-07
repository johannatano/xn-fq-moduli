from __future__ import annotations

from fractions import Fraction
from math import sqrt
import time

from xnfq.config import apply_config_from_args
from xnfq.moduli.modular_curve import X1
from xnfq.arithmetic.quadratic import LatticeTower
from utils.args import parse_example_args
from utils.fmt import fmt_magnitude
from utils.ui.dashboard import Dashboard, Param


def _prime_choices(stop: int) -> list[int]:
    primes: list[int] = []
    for candidate in range(2, stop + 1):
        is_prime = True
        for prime in primes:
            if prime * prime > candidate:
                break
            if candidate % prime == 0:
                is_prime = False
                break
        if is_prime:
            primes.append(candidate)
    return primes


PRIME_CHOICES = _prime_choices(97)


class X1FqPlotView(Dashboard):
    title = "X1(N)(F_q) strata"
    figsize = (13.0, 6.5)

    def params(self):
        return [
            Param("p", 2, 97, 5, step=1, label="prime p", choices=PRIME_CHOICES),
            Param("n", 1, 8, 1, step=1, label="degree n"),
            Param("N", 1, 5000, 11, step=1, label="level N"),
        ]

    def panels(self):
        return ["strata", "contrib"]

    def _snapshot(self):
        key = (self["p"], self["n"], self["N"])
        if getattr(self, "_cache_key", None) == key:
            return self._cache_data

        p, n, N = key
        q = p**n

        start_t = time.time()
        x1n = X1(N).over(p, n)
        grouped_map = {}
        for fiber in x1n.yfq():
            entry = grouped_map.get(fiber.lattice.DK)
            if entry is None:
                grouped_map[fiber.lattice.DK] = (fiber.lattice, [fiber])
            else:
                entry[1].append(fiber)
        grouped = [
            (DK, lattice, fibers) for DK, (lattice, fibers) in grouped_map.items()
        ]
        y1_count = sum(
            (sum(stratum.count() for stratum in strata) for _, _, strata in grouped),
            start=Fraction(0),
        )
        cusp_count = sum((cusp.count() for cusp in x1n.cusps()), start=Fraction(0))
        run_t = time.time() - start_t

        self._cache_key = key
        self._cache_data = {
            "p": p,
            "n": n,
            "N": N,
            "q": q,
            "grouped": grouped,
            "y1_count": y1_count,
            "cusp_count": cusp_count,
            "run_t": run_t,
        }
        return self._cache_data

    def draw(self, name: str, ax) -> None:
        snapshot = self._snapshot()
        p = snapshot["p"]
        n = snapshot["n"]
        N = snapshot["N"]
        q = snapshot["q"]
        grouped = snapshot["grouped"]
        y1_count = snapshot["y1_count"]
        cusp_count = snapshot["cusp_count"]
        run_t = snapshot["run_t"]

        if name == "strata":
            points: list[tuple[int, float]] = []
            for DK, lattice, _strata in grouped:
                for conductor in sorted(LatticeTower.divisors(lattice.conductors)):
                    print(f"Conductor: {conductor}, lattice.f_max={lattice.f_max}")
                    normalized = 1.0 if DK == 0 else conductor / lattice.f_max
                    points.append((DK, normalized))

            if points:
                x_values = [DK for DK, _ in points]
                y_values = [normalized for _, normalized in points]
                ax.scatter(x_values, y_values, s=24, color="steelblue", alpha=0.8)
                ax.set_xscale("symlog", linthresh=1, base=10)
                ax.set_xlim(min(x_values) - 1, max(x_values) + 1)
                ax.set_ylim(-0.02, 1.02)
            else:
                ax.text(0.5, 0.5, "no Y1 strata", ha="center", va="center")

            ax.set_title(f"Y1 strata over F_{p}^{n} for X1({N})")
            ax.set_xlabel("D_K")
            ax.set_ylabel(r"$f/f_{\max}$")
            ax.grid(color="0.92", linewidth=0.8)

            summary = (
                f"q={q} ({fmt_magnitude(q)})\n"
                f"Y1(N)={y1_count}\n"
                f"Cusps={cusp_count}\n"
                f"fields={len(grouped)}\n"
                f"time={run_t:.6f}s"
            )
            ax.text(
                0.98,
                0.98,
                summary,
                transform=ax.transAxes,
                ha="right",
                va="top",
                bbox={"facecolor": "white", "edgecolor": "0.8", "boxstyle": "round,pad=0.4"},
            )
            return

        contrib_points = [
            (DK, sum((stratum.count() for stratum in strata), start=Fraction(0)))
            for DK, _lattice, strata in grouped
        ]
        if contrib_points:
            x_values = [DK for DK, _ in contrib_points]
            y_values = [float(total) for _, total in contrib_points]
            ax.scatter(x_values, y_values, s=32, color="darkorange", alpha=0.85)
            ax.vlines(x_values, 0.0, y_values, color="darkorange", alpha=0.45, linewidth=1.2)
            ax.set_xscale("symlog", linthresh=1, base=10)
            ax.set_xlim(min(x_values) - 1, max(x_values) + 1)
            ax.set_ylim(0.0, max(y_values) * 1.08 if y_values else 1.0)
        else:
            ax.text(0.5, 0.5, "no Y1 contributions", ha="center", va="center")

        ax.set_title(f"Total Y1 contribution by D_K over F_{p}^{n}")
        ax.set_xlabel("D_K")
        ax.set_ylabel("total contribution")
        ax.grid(color="0.92", linewidth=0.8)


def run() -> None:
    args = parse_args()
    apply_config_from_args(args)
    X1FqPlotView(p=args.p, n=args.n, N=args.N).show()


def parse_args():
    return parse_example_args("Interactive X1(N)(F_q) plot")


if __name__ == "__main__":
    run()
