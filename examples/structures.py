from __future__ import annotations
import time
from math import comb, gcd

from xnfq.config import apply_config_from_args
from xnfq.moduli.level_structures import Gamma, Gamma0, Gamma1
from xnfq.moduli.modular_curve import X1
from xnfq.arithmetic.function import phi
from xnfq.arithmetic.common import kronecker

from utils.args import parse_example_args
from utils.fmt import fmt_magnitude
from utils.lmfdb_api import fetch_traces

from utils.logging import Logger, Colors

def hk(t: int, q: int, k: int) -> int:
    return sum(comb(k - j, j) * (-q) ** j * t ** (k - 2 * j) for j in range(k // 2 + 1))

def tr_frob_symk(x1n, k: int) -> int:
    """Compute the trace of Frobenius on the `k`th symmetric power."""
    val = 0
    #TODO: not sure if this always holds?
    if k == 0:
        val = x1n.q + (1 if gcd(x1n.q, x1n.N) == 1 else 0)
    curves_term = sum(hk(c.t, x1n.q, k) * c.count() for c in x1n.yfq())
    cusp_term = sum((c.t ** (k + 2)) * c.count() for c in x1n.cusps())
    return val - curves_term - cusp_term


def run():
    args = parse_args()
    apply_config_from_args(args)
    p = args.p
    n = args.n
    N = args.N
    q = p ** n

    fast_trace = args.fast_trace
    using_pari = args.use_pari if hasattr(args, "use_pari") else False
    Logger.cprint(
        f"Computing S_{args.k+2}(Frob_{q} | Sym^{args.k}) | Optimized trace enum: {fast_trace} | Using PARI: {using_pari} | q mag={fmt_magnitude(q)}",
        Colors.HEADER,
    )

    start_t = time.time()
    x1n = X1(N).over(p, n)
    trace = tr_frob_symk(x1n, args.k)
    stop_t = time.time()

    #Y0N = x1n.y_count(Gamma0)
    #print("Y0N:", Y0N)
    #Y1N = x1n.y_count(Gamma1)
    #print("Y1N:", Y1N)
    # full structure
    full_structure = x1n.get_structure()
    max_val = 0
    max_str = None
    max_c = None
    for s_rec in full_structure.fibers:
        if s_rec.kind != "weil":
            continue
        normalized_contrib = s_rec.total_count / (phi(1)(N) * phi(-1)(N))
        if normalized_contrib > max_val:
            max_val = normalized_contrib
            max_str = s_rec
        # if normalized_contrib < 1:
        #    continue
        # split_type=kronecker(s_rec.discriminant, N)
        # clr = Colors.CYAN if split_type == 1 else Colors.MAGENTA
        '''Logger.cprint(
            f"Stratum DK: {s_rec.discriminant}, Conductors: {s_rec.conductors}, Total Mass: {s_rec.total_mass}, Total Count: {s_rec.total_count}, normalized contribution: {normalized_contrib:.1f}, split_type={split_type}",
            clr
        )'''

    if max_str is not None:
        Logger.cprint(
            f"Maximum Structure: {max_val}, Discriminant: {max_str.discriminant} trace={max_str.trace}, split_type={kronecker(max_str.discriminant, N)}",
            Colors.MAGENTA,
        )
        for level_rec in max_str.level_records:
            clr = (
                Colors.GREEN
                if level_rec.full
                else Colors.RED if level_rec.cyclic == 0 else Colors.YELLOW
            )
            Logger.cprint(f"Level record:{level_rec}", clr)
            normalized_contrib += level_rec.mass * level_rec.cyclic
    # y = x1n.open_count(Gamma)
    # print("Y:", y)

    ref = None
    if args.sage:
        from sage.all import CuspForms, Gamma1 as SageGamma1
        start_t_sage = time.time()
        form = CuspForms(SageGamma1(N), args.k + 2)
        ref = int(form.hecke_operator(q).trace())
        stop_t_sage = time.time()
        error = abs(trace - ref)
        clr = Colors.YELLOW if error == 0 else Colors.RED
        Logger.cprint(f"Result: trace={trace}, sage ref={ref}, error={error}, time={(stop_t - start_t):.9f}, sage time={(stop_t_sage - start_t_sage):.9f}", clr)
    else:
        error = None
        clr = Colors.YELLOW
        Logger.cprint(
            f"Result: trace={trace}, time={(stop_t - start_t):.9f}", clr
        )


def parse_args():
    return parse_example_args(
        "Compute Tr Sym^k on X1(N)", include_sym_power=True, include_sage=True
    )


if __name__ == "__main__":
    run()
