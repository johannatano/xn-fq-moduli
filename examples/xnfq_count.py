# Generic imports
from __future__ import annotations
from math import sqrt, pi, gcd, prod, isqrt, comb
from fractions import Fraction
import time

# X1 Library Specific Imports
from xnfq.config import apply_config_from_args
from xnfq.moduli.modular_curve import X, X0, X1
from utils.args import parse_example_args
from utils.fmt import fmt_magnitude
from utils.logging import Logger, Colors


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
    YN, CN = XN.count()
    # ================================================================

    # ===== Results ==================================================
    run_t = time.time() - start_t
    print(
        f"X_{args.type}({N}) fiber count over F_{p}^{n} : #Y_{args.type}({N}) = {YN}, #Cusps_{args.type}({N}) = {CN} | (t={(run_t):.2f}s)"
    )
    # ================================================================


# ===== args parsing ============================================
def parse_args():
    return parse_example_args(
        "Compute Fq rational Points on X1(N)", include_prime_range=True
    )
# ================================================================


if __name__ == "__main__":
    run()
