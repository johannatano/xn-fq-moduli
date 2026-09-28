from __future__ import annotations
from fractions import Fraction
import time

from xnfq.config import apply_config_from_args

from xnfq.moduli import EigenFormRecord, CuspFiberRecordFq, SmoothFiberRecordFq, ModularCurve, Gamma

from xnfq.arithmetic.forms import BinaryQuadraticForm
from xnfq.arithmetic.common import factorize

from utils.args import parse_example_args
from utils.fmt import fmt_magnitude, fmt_factored, fmt_invariants
from utils.logging import Logger, Colors
from utils.data import ResultData


def format_eigenform_data(
    fiber: CuspFiberRecordFq | SmoothFiberRecordFq, eigenform_rec: EigenFormRecord, q: int, gamma: Gamma
) -> list[list[ResultData]]:
    table: list[list[ResultData]] = []
    form = eigenform_rec.eigenform
    for level in eigenform_rec.levels:
        row = []
        if fiber.kind == "smooth":
            row.append(ResultData("f", level.index))
        else:
            row.append(ResultData("d", level.inv[1]))
        row.append(ResultData("mass", level.mass))
        row.append(ResultData("inv", level.inv))
        row.append(ResultData("stable lines", level.num_lines))
        clr = Colors.DIM if level.gamma_count == 0 else Colors.BOLD
        row.append(ResultData("structure count", level.gamma_count))
        row.append(
            ResultData(
                "level contrib", level.gamma_count * level.mass * fiber.m0, fmt=clr
            )
        )
        table.append(row)
    return table

def run():
    # ===== init ====================================================
    args = parse_args()
    apply_config_from_args(args)
    pari_info = " | Using PARI" if args.use_pari else ""
    N, p, n = args.N, args.p, args.n
    q = p**n
    start_t = time.time()
    # ================================================================

    # ===== Main Computation =========================================
    gamma = Gamma(N, args.type)
    XN = ModularCurve(gamma).over(p, n)
    Logger.cprint(
        f"Begin Structure Report of {XN.info()} | q mag={fmt_magnitude(q)}{pari_info}",
        Colors.NEON_PURPLE,
    )
    # Count the number of rational points on Y1(N) including cusps
    cusp_fibers = []

    for fiber in XN.get_structure():
        # set this to the minimum number of levels to probe, eg only display if num levels more than 3
        probe_len = -1
        if fiber.count == 0 or not any(
            len(record.levels) >= probe_len for record in fiber.eigen_records
        ):
            continue
        fiber_data = f", D_K: {fiber.D_K}, m0: {fiber.m0}" if fiber.kind == "smooth" else ""
        clr = Colors.NEON_LIME if fiber.kind == "smooth" else Colors.NEON_PINK
        Logger.cprint(
            f"{fiber.kind.capitalize()} fiber of X{args.type}({N})(F_{p}^{n}) | t: {fiber.t}{fiber_data}, size: {fiber.count}",
            clr,
        )
        for record in fiber.eigen_records:
            if args.type == 0:
                Logger.cprint(
                    f"λ={record.eigenform.value}",
                    Colors.NEON_PURPLE,
                )
            Logger.print_results(
                format_eigenform_data(fiber, record, q, gamma)
            )
            print()

    # ================================================================

# ===== args parsing ============================================
def parse_args():
    return parse_example_args(
        "",
        overrides={
            "p": 13,
            "n": 4,
            "N": 12,
        }
    )
# ================================================================


if __name__ == "__main__":
    run()
