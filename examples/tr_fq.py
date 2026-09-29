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
from utils.sage_api import tr_fq


def format_trace_result(
    p: int,
    n:int,
    rec: TrFqTraceRecord,
    rec_ref: TrFqTraceRecord
) -> list[ResultData]:
    cols = [
        ResultData("p" if n == 1 else "q", p**n),
        ResultData("trace", rec.val),
        ResultData("smooth", rec.smooth),
        ResultData("cusp", rec.cusp),
        ResultData("eps0", rec.eps0, show=rec.eps0 > 0),
    ]
    if rec_ref is not None:
        cols.append(ResultData("sage ref", rec_ref.val))
        cols.append(ResultData("error", abs(rec.val - rec_ref.val)))
    return cols


def run():
    # ===== init ====================================================
    args = parse_args()
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

    result_rows: list[list[ResultData]] = []
    logs: list[str] = []
    for pp in args.primes:
        XN = ModularCurve(Gamma(N, args.type)).over(pp, n)
        trace_rec = XN.tr_fq(args.k)
        trace_rec_ref = None
        # automatic verification using sage, TODO: use LMFDB API instead
        if args.sage:
            trace_rec_ref, msg = tr_fq(args, pp, n, N, args.k)
            logs.append(msg)
            
        result_rows.append(format_trace_result(pp, n, trace_rec, trace_rec_ref))

    Logger.header(f"Traces across prime range (N={N}, k={args.k}, type={args.type})", Colors.NEON_LIME)
    Logger.print_results(result_rows, color_fn=lambda row: Colors.BOLD)

    stop_t = time.time()

    for m in logs:
        if m is not None and m != "":
            Logger.cprint(m, Colors.WARNING)

    Logger.cprint(f"Computed {len(result_rows)} traces in {(stop_t - start_t):.6f}s", Colors.NEON_LIME)

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
