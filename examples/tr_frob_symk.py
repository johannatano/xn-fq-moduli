from __future__ import annotations
import time
from math import comb, gcd

from xnfq.config import apply_config_from_args
from xnfq.moduli.level_structures import Gamma, Gamma0, Gamma1
from xnfq.moduli.modular_curve import ModularCurve, X, X0, X1
from xnfq.arithmetic.function import phi
from xnfq.arithmetic.common import kronecker

from utils.args import parse_example_args
from utils.fmt import fmt_magnitude
from utils.lmfdb_api import fetch_traces

from utils.logging import Logger, Colors
from sympy import primerange


def run():
    args = parse_args()
    apply_config_from_args(args)
    p = args.p
    n = args.n
    N = args.N
    q = p ** n

    #primes = list(primerange(2, 500))
    #for _p in primes:
    #    if (_p-1) % N == 0:
    #        Logger.cprint(f"Prime {_p} satisfies (p-1) % N={(_p-1) % N}, (q-1) % N={(q-1) % N}", Colors.GREEN)

    fast_trace = args.fast_trace
    using_pari = args.use_pari if hasattr(args, "use_pari") else False
    
    Logger.cprint(
        f"Computing S_{args.k+2}(Frob_{q} | Sym^{args.k}) | Optimized trace enum: {fast_trace} | Using PARI: {using_pari} | q mag={fmt_magnitude(q)}",
        Colors.HEADER,
    )

    start_t = time.time()
    q = p**n
    k = args.k
    mod_curve = None
    if args.type == 1:
        mod_curve = X1(N).over(p, n)
    elif args.type == 0:
        mod_curve = X0(N).over(p, n)
    else:
        mod_curve = X(N).over(p, n)

    trace = mod_curve.tr_frob_symk(args.k)
    stop_t = time.time()

    ref = None
    if args.sage:
        from sage.all import (
            CuspForms,
            Gamma1 as SageGamma1,
            Gamma0 as SageGamma0,
            ModularSymbols,
            GammaH,
        )

        start_t_sage = time.time()
        congruence_subgroup = None
        if args.type == 1:
            congruence_subgroup = SageGamma1(N)
            form = CuspForms(congruence_subgroup, args.k + 2)
            ref = int(form.hecke_operator(q).trace())
        elif args.type == 0:
            congruence_subgroup = SageGamma0(N)
            form = CuspForms(congruence_subgroup, args.k + 2)
            ref = int(form.hecke_operator(q).trace())
        else:
            H = GammaH(N**2, [1 + N])
            form = CuspForms(H, args.k + 2)
            ref = int(form.hecke_operator(q).trace())

        stop_t_sage = time.time()
        error = abs(trace - ref)
        clr = Colors.GREEN if error == 0 else Colors.RED
        Logger.cprint(f"Result: trace={trace}, sage ref={ref}, error={error}, time={(stop_t - start_t):.9f}, sage time={(stop_t_sage - start_t_sage):.9f}", clr)
    else:
        error = None
        clr = Colors.GREEN
        Logger.cprint(
            f"Result: trace={trace}, time={(stop_t - start_t):.9f}", clr
        )

def parse_args():
    return parse_example_args(
        "Compute Tr Sym^k on X1(N)", include_sym_power=True, include_sage=True
    )

if __name__ == "__main__":
    run()
