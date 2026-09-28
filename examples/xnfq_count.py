from __future__ import annotations
import time

from xnfq.config import apply_config_from_args
from xnfq.moduli import ModularCurve, Gamma
from utils.args import parse_example_args
from utils.fmt import fmt_magnitude
from utils.logging import Logger, Colors

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
    # construct the modular curve over finite field
    XN = ModularCurve(Gamma(N, args.type)).over(p, n)
    Logger.cprint(
        f"Counting points on {XN.info()} | q mag={fmt_magnitude(q)}{pari_info}",
        Colors.NEON_PURPLE,
    )
    # count the points on the modular curve
    YN, CN = XN.count()
    # ================================================================

    # ===== Results ==================================================
    run_t = time.time() - start_t
    Logger.cprint(
        f"Result: #X{args.type}({N})(F_{p}^{n})={YN + CN}, #Smooth={YN}, #Cusps={CN}, time={(run_t):.9f}",
        Colors.NEON_LIME,
    )
    # ================================================================

# ===== args parsing ============================================
def parse_args():
    return parse_example_args(
        "",
        overrides={
            "p": 13,
            "n": 4,
            "N": 5,
        },
    )
# ================================================================

if __name__ == "__main__":
    run()
