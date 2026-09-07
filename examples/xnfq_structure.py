# Generic imports
from __future__ import annotations
from math import sqrt, pi, gcd, prod, isqrt, comb
from fractions import Fraction
import time

# X1 Library Specific Imports
from xnfq.config import apply_config_from_args
from xnfq.moduli.data import EigenFormRecord, FiberRecord
from xnfq.moduli.modular_curve import X, X0, X1
from utils.args import parse_example_args
from utils.fmt import fmt_magnitude, fmt_factored, fmt_invariants
from utils.logging import Logger, Colors
from utils.data import ResultData

"""
def format_weil_fiber_data(data: FiberRecord, N: int) -> list[list[ResultData]]:
    rows: list[list[ResultData]] = []
    for record in data.level_records:
        rows.append(
            [
                ResultData("f", record.index),
                ResultData("l", record.l),
                ResultData("Im(pi-lam)", record.inv),
                ResultData("ker(pi-lam)", record.tor_subgrp),
                ResultData("mass", record.mass),
                ResultData("cyclic", record.cyclic),
                ResultData("level count", record.count),
                ResultData("full", record.full),
            ]
        )
    return rows"""


def format_eigenform_fiber_data(data: EigenFormRecord, N: int) -> list[list[ResultData]]:
    rows: list[list[ResultData]] = []
    eigenform = data.eigenform

    for level in data.level_records:
        rows.append(
            [
                ResultData("f", level.index),
                # ResultData(f"Im(pi-{eigenform.eigenvalue})", level.im),
                ResultData(f"ker(pi-{eigenform.eigenvalue})", level.ker),
                ResultData("mass", level.mass),
                ResultData("local-basis", level.coords),
                #ResultData("local-basis-pi", level.coords_pi_basis),
                ResultData("---form---", ((eigenform.form.B, (eigenform.form.C)))),
                ResultData("---order (B,C)---", ((level.order.B, (level.order.C)))),
                ResultData(f"C_{N}", level.num_lines),
            ]
        )
    """rows.append(
        [
            #ResultData("f", data.index),
            ResultData("l", eigenform.l),
            ResultData(f"Im(pi-{eigenform.eigenvalue})", eigenform.im),
            ResultData(f"ker(pi-{eigenform.eigenvalue})", eigenform.ker),
            ResultData("mass", data.mass),
            ResultData("cyclic", eigenform.num_lines),
        ]
    )"""
    return rows


def format_cusp_fiber_data(data: FiberRecord) -> list[ResultData]:
    return [
        ResultData("row", "cusp"),
        ResultData("d", data.d),
        ResultData("t", data.t),
        ResultData("count", data.total_count),
    ]


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
    full_structure = XN.get_structure()
    cusp_fibers = []
    for s_rec in full_structure.fibers:
        if s_rec.kind == "weil" and s_rec.total_count >= 0:
            
            probe_len = 2
            if not any(len(record.level_records) >= probe_len for record in s_rec.eigen_records):
                continue
                
            Logger.cprint(
                f"Fiber over DK={s_rec.d}, t={s_rec.t}, coprime={s_rec.coprime},s_rec.total_count={s_rec.total_count}",
                Colors.YELLOW,
            )
            for record in s_rec.eigen_records:
                
                Logger.print_results(format_eigenform_fiber_data(record, N))
                print()
        elif s_rec.kind == "cusp":
            cusp_fibers.append(format_cusp_fiber_data(s_rec))
        else:
            continue
    # ================================================================

    # ===== Results ==================================================
    run_t = time.time() - start_t
    Logger.cprint(
        f"Cusp Fibers:",
        Colors.YELLOW,
    )
    Logger.print_results(cusp_fibers)
    # print(
    #    f"X_{args.type}({N}) fiber count over F_{p}^{n} : #Y_{args.type}({N}) = {YN}, #Cusps_{args.type}({N}) = {CN} | (t={(run_t):.2f}s)"
    # )
    # ================================================================


# ===== args parsing ============================================
def parse_args():
    return parse_example_args(
        "Compute Fq rational Points on X1(N)", include_prime_range=True
    )
# ================================================================


if __name__ == "__main__":
    run()
