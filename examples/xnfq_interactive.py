from __future__ import annotations

import sys
from fractions import Fraction
import time

from xnfq.config import apply_config_from_args
from xnfq.moduli.modular_curve import X0, X1, X
from utils.args import parse_example_args
from utils.fmt import fmt_magnitude
from utils.ui.dashboard import Dashboard, Param

def XN_over(_type: int, N: int, p: int, n: int):
    if _type == 0:
        return X0(N).over(p, n)
    elif _type == 1:
        return X1(N).over(p, n)
    elif _type == 2:
        return X(N).over(p, n)
    else:
        raise ValueError(f"Unsupported type: {_type}")


def curve_label(_type: int) -> str:
    return {0: "X_0", 1: "X_1", 2: "X"}[_type]


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


class XNFqPlotView(Dashboard):
    figsize = (14.0, 6.5)

    def params(self):
        return [
            Param("p", 2, 97, 43, step=1, label="p", choices=PRIME_CHOICES),
            Param("n", 1, 8, 1, step=1, label="n"),
            Param("N", 1, 100, 3, step=1, label="N"),
            Param("type", 0, 2, 1, step=1, label="Gamma"),
        ]

    def panels(self):
        return ["smooth", "cusps"]

    def header(self) -> str:
        q = self["p"] ** self["n"]
        return f"{curve_label(self['type'])}({self['N']}) over F_{q} strata"

    def _snapshot(self):
        key = (self["p"], self["n"], self["N"], self["type"])
        if getattr(self, "_cache_key", None) == key:
            return self._cache_data

        p, n, N, _type = key
        q = p**n

        if q > 1_000_000:
            raise ValueError(f"q = {q} is too large")

        XN = XN_over(_type, N, p, n)
        report = XN.get_structure()
        smooth_weight = XN.level_structure.weight(smooth=True)
        cusp_weight = XN.level_structure.weight(smooth=False)
        smooth_by_dk: dict[int, Fraction] = {}
        cusp_by_d: dict[int, Fraction] = {}

        for fiber in report:
            if fiber.kind == "smooth":
                smooth_by_dk[fiber.D_K] = (
                    smooth_by_dk.get(fiber.D_K, Fraction(0))
                    + fiber.count
                )
            elif fiber.kind == "cusp":
                for eigen_record in fiber.eigen_records:
                    for level in eigen_record.levels:
                        d = level.inv[1]
                        cusp_by_d[d] = cusp_by_d.get(d, Fraction(0)) + (
                            Fraction(level.num_lines, d)
                            * Fraction(cusp_weight, 2)
                        )

        yn_count = sum(smooth_by_dk.values(), start=Fraction(0))
        cusp_count = sum(cusp_by_d.values(), start=Fraction(0))
        smooth_points = sorted(
            (d_k, total / smooth_weight)
            for d_k, total in smooth_by_dk.items()
            if total != 0
        )
        cusp_points = sorted(
            (d, total / cusp_weight)
            for d, total in cusp_by_d.items()
            if total != 0
        )

        self._cache_key = key
        self._cache_data = {
            "p": p,
            "n": n,
            "N": N,
            "q": q,
            "smooth_by_dk": smooth_points,
            "cusp_by_d": cusp_points,
            "yn_count": yn_count,
            "cusp_count": cusp_count,
        }
        return self._cache_data

    def draw(self, name: str, ax) -> None:
        snapshot = self._snapshot()
        p = snapshot["p"]
        n = snapshot["n"]
        N = snapshot["N"]
        q = snapshot["q"]
        smooth_by_dk = snapshot["smooth_by_dk"]
        cusp_by_d = snapshot["cusp_by_d"]
        yn_count = snapshot["yn_count"]
        cusp_count = snapshot["cusp_count"]
        label = curve_label(self["type"])
        if name == "smooth":
            if smooth_by_dk:
                x_values = [d_k for d_k, _ in smooth_by_dk]
                y_values = [float(total) for _, total in smooth_by_dk]
                # ax.scatter(x_values, y_values, s=24, color="black", alpha=1)
                ax.vlines(
                    x_values, 0.0, y_values, color="black", alpha=.5, linewidth=1.5
                )
                ax.set_xscale("symlog", linthresh=1, base=10)
                ax.set_xlim(min(x_values) - 1, max(x_values) + 1)
                y_max = max(y_values, default=0.0)
                ax.set_ylim(0.0, max(1.0, y_max * 1.08))
            else:
                ax.text(0.5, 0.5, "no smooth fibers", ha="center", va="center")

            ax.set_title(f"Smooth Fibers")
            ax.set_xlabel("D_K")
            ax.tick_params(axis="y", left=False, labelleft=False)
            ax.grid(color="0.92", linewidth=0.8)

            summary = (
                f"total count={yn_count}"
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

        if cusp_by_d:
            x_values = [d for d, _ in cusp_by_d]
            y_values = [float(total) for _, total in cusp_by_d]
            # ax.scatter(x_values, y_values, s=32, color="black", alpha=1)
            ax.vlines(x_values, 0.0, y_values, color="black", alpha=.5, linewidth=1.5)
            ax.set_xticks(x_values)
            ax.set_xlim(min(x_values) - 1, max(x_values) + 1)
            y_max = max(y_values, default=0.0)
            ax.set_ylim(0.0, max(1.0, y_max * 1.08))
        else:
            ax.text(0.5, 0.5, "no cusp fibers", ha="center", va="center")

        ax.set_title(f"Cusp Fibers")
        ax.set_xlabel("d")
        ax.tick_params(axis="y", left=False, labelleft=False)
        ax.grid(color="0.92", linewidth=0.8)
        summary = (
            f"total count={cusp_count}"
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


def run() -> None:
    args = parse_args()
    apply_config_from_args(args)
    overrides = {}
    for name in ("p", "n", "N", "type"):
        value = getattr(args, name)
        if value is not None:
            overrides[name] = value
    XNFqPlotView(**overrides).show()

# ===== args parsing ============================================
def parse_args():
    args = parse_example_args(
        "Compute Fq rational Points on X1(N)", include_prime_range=True
    )
    options = set(sys.argv[1:])
    if not {"-p", "--p"} & options:
        args.p = None
    if "-n" not in options and "--n" not in options:
        args.n = None
    if "-N" not in options and "--N" not in options:
        args.N = None
    if "-type" not in options and "--type" not in options:
        args.type = None
    return args
# ================================================================

if __name__ == "__main__":
    run()
