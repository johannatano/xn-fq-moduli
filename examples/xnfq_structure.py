# Generic imports
from __future__ import annotations
from math import sqrt, pi, gcd, prod, isqrt, comb
from fractions import Fraction
import time

# X1 Library Specific Imports
from xnfq.config import apply_config_from_args

from xnfq.moduli import EigenFormRecord, CuspFiberRecordFq, SmoothFiberRecordFq, X, X0, X1

from xnfq.arithmetic.forms import BinaryQuadraticForm
from xnfq.arithmetic.common import factorize

from utils.args import parse_example_args
from utils.fmt import fmt_magnitude, fmt_factored, fmt_invariants
from utils.logging import Logger, Colors
from utils.data import ResultData

from sympy import primerange

def format_eigenform_data(
    fiber: CuspFiberRecordFq | SmoothFiberRecordFq, eigenform_rec: EigenFormRecord, q: int, N: int
) -> list[list[ResultData]]:
    rows: list[list[ResultData]] = []
    form = eigenform_rec.eigenform

    for level in eigenform_rec.levels:

        if fiber.kind == "smooth":
            rows.append(
                [
                    ResultData("f", level.index),
                    ResultData("mass", level.mass),
                    ResultData("num_lines", level.num_lines),
                    #ResultData(f"coords", level.coords),
                    ResultData(f"inv", level.inv),
                ]
            )
        else:
            rows.append(
                [
                    ResultData("d", level.inv[1]),
                    ResultData("mass", level.mass),
                    ResultData("num_lines", level.num_lines),
                    ResultData(f"inv", level.inv),
                ]
            )
    return rows


def XN_over(_type:int, N:int, p:int, n:int):
    if _type == 0:
        return X0(N).over(p, n)
    elif _type == 1:
        return X1(N).over(p, n)
    elif _type == 2:
        return X(N).over(p, n)
    else:
        raise ValueError(f"Unsupported type: {_type}")
def run():
    # ===== init ====================================================
    args = parse_args()
    fast_trace = args.fast_trace
    using_pari = args.use_pari
    apply_config_from_args(args)
    N, p, n = args.N, args.p, args.n
    q = p**n
    start_t = time.time()
    # ================================================================

    # ===== Main Computation =========================================
    Logger.cprint(
        f"ARGS: Optimized trace enum: {fast_trace and args.type == 1}, Using PARI: {using_pari} | q mag={fmt_magnitude(q)}",
        Colors.HEADER,
    )
    XN = XN_over(args.type, N, p, n)
    # Count the number of rational points on Y1(N) including cusps
    cusp_fibers = []

    for fiber in XN.get_structure():
        # set this to the minimum number of levels to probe
        probe_len = -1
        if fiber.count == 0 or not any(
            len(record.levels) >= probe_len for record in fiber.eigen_records
        ):
            continue
        fiber_data = f", D_K: {fiber.D_K}" if fiber.kind == "smooth" else ""
        Logger.cprint(
            f"{fiber.kind.capitalize()} fiber over t: {fiber.t}{fiber_data}, size: {fiber.count}",
            Colors.YELLOW,
        )
        for record in fiber.eigen_records:
            if args.type == 0:
                Logger.cprint(
                    f"Probe level count: {len(record.levels)}",
                    Colors.CYAN,
                )
            Logger.print_results(
                format_eigenform_data(fiber, record, q, N)
            )
            print()

    # ================================================================
    # ===== Results ==================================================
    # ================================================================


# ===== args parsing ============================================
def parse_args():
    return parse_example_args(
        "Compute Fq rational Points on X1(N)", include_prime_range=True
    )
# ================================================================


if __name__ == "__main__":
    run()
