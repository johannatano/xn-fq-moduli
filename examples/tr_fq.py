from __future__ import annotations
import sys
import time
from math import comb, gcd

from xnfq.config import apply_config_from_args
from xnfq.moduli import ModularCurve, Gamma

from utils.args import parse_example_args
from utils.fmt import fmt_magnitude
from utils.logging import Logger, Colors
from utils.data import ResultData
from utils.sage import compare_trace_with_sage, get_sage_ref

def run():
    # ===== init ====================================================
    args = parse_args()
    # no plotting in this example; it's a simple console trace table
    using_pari = args.use_pari if hasattr(args, "use_pari") else False
    apply_config_from_args(args)
    p, n, N = args.p, args.n, args.N
    q = p ** n
    start_t = time.time()
    # ================================================================

    # ===== Main Computation =========================================
    pari_info = " | Using PARI" if using_pari else ""
    Logger.cprint(
        f"Computing Tr(F_{p}^{n}, S[Gamma{args.type}({N}), k={args.k})) | q mag={fmt_magnitude(q)}{pari_info}",
        Colors.NEON_PURPLE,
    )

    # Determine whether to run across a prime range or just the baseline `-p`.
    # Behavior: `-p` is always the baseline. Use `--range R` to include primes
    # in the interval [p, p+R). If `--range` is not provided, run only for `p`.
    primes: list[int]
    if getattr(args, "range", None) is not None:
        # Interpret `--range` as a count: include `range` primes starting at `p`.
        from utils.primes import next_primes

        cnt = max(1, int(args.range))
        primes = next_primes(int(args.p), cnt)
    else:
        primes = [int(args.p)]

    results = []
    sage_msgs: list[str] = []
    for pp in primes:
        XN = ModularCurve(Gamma(N, args.type)).over(pp, n)
        # trace = XN.tr_fq(args.k) - we rebuild the trace manually using smooth fibers and cusps for transparency
        curves_term = sum(XN.hk(c.t, XN.q, args.k - 2) * c.count() for c in XN.smooth_fibers())
        cusp_term = sum((c.t ** (args.k)) * c.count() for c in XN.cusps())
        # correction term if k-2 = 0
        eps0 = 0
        if args.k == 2:
            eps0 = XN.q + (1 if gcd(XN.q, XN.N) == 1 else 0)

        # final value
        trace = eps0 - curves_term - cusp_term

        # optionally obtain a Sage reference and error when requested and on prime fields
        sage_ref = None
        sage_err = None
        sage_time = 0.0
        if getattr(args, "sage", False):
            if args.n != 1:
                sage_msgs.append("Sage comparison only supports prime fields (n=1); skipping.")
            else:
                ref, msg, stime = get_sage_ref(args, pp, N, args.k)
                sage_time = stime
                if msg is not None:
                    sage_msgs.append(msg)
                else:
                    sage_ref = ref
                    sage_err = abs(trace - ref)

        results.append((pp, trace, curves_term, cusp_term, eps0, sage_ref, sage_err, sage_time))

    # build ResultData rows matching the structure report style
    rows: list[list[ResultData]] = []
    label_name = "p" if args.n == 1 else "q"
    show_eps0 = args.k == 2
    for pp, trace, curves_term, cusp_term, eps0, sage_ref, sage_err, sage_time in results:
        display = pp if args.n == 1 else (pp ** args.n)
        sage_ref_col = None
        sage_err_col = None
        if getattr(args, "sage", False) and args.n == 1:
            sage_ref_col = sage_ref if sage_ref is not None else "N/A"
            sage_err_col = sage_err if sage_err is not None else "N/A"
        rows.append(
            [
                ResultData(label_name, display),
                ResultData("trace", trace),
                ResultData("smooth", curves_term),
                ResultData("cusp", cusp_term),
                ResultData("eps0", eps0, show=show_eps0),
                ResultData("sage ref", sage_ref_col, show=(sage_ref_col is not None)),
                ResultData("error", sage_err_col, show=(sage_err_col is not None)),
            ]
        )
    Logger.header(f"Traces across prime range (N={N}, k={args.k}, type={args.type})", Colors.NEON_LIME)
    Logger.print_results(rows, color_fn=lambda row: Colors.BOLD)

    # Results timing
    stop_t = time.time()

    # Print any Sage import / congruence messages collected while building rows
    for m in sage_msgs:
        Logger.cprint(m, Colors.WARNING)

    if not getattr(args, "sage", False):
        Logger.cprint(f"Computed {len(results)} traces in {(stop_t - start_t):.6f}s", Colors.NEON_LIME)
    else:
        # If Sage was used, summarize timing
        total_sage = sum(r[7] for r in results)
        Logger.cprint(f"Computed {len(results)} traces (sage time {total_sage:.6f}s) in {(stop_t - start_t):.6f}s", Colors.NEON_LIME)


# ===== args parsing ============================================
def parse_args():
    return parse_example_args(
        "",
        include_weight=True,
        include_prime_range=True,
        include_sage=True,
        overrides={
            "p": 7,
            "n": 1,
            "N": 3,
        },
    )


if __name__ == "__main__":
    run()
