from __future__ import annotations

import time
from sympy import primerange

from utils.logging import Logger, Colors


def compare_trace_with_sage(args, p: int, N: int, trace: int, start_t: float, stop_t: float) -> bool:
    """Run Sage comparison and print results. Returns True when handled.

    This helper prints the result (or error) using `Logger` and returns True
    to indicate the caller should exit after the comparison.
    """
    try:
        from sage.all import (
            CuspForms,
            Gamma1 as SageGamma1,
            Gamma0 as SageGamma0,
            GammaH,
        )
    except Exception as e:
        Logger.cprint(f"Unable to import Sage: {e}", Colors.FAIL)
        Logger.cprint(f"Result: trace={trace}, time={(stop_t - start_t):.9f}", Colors.NEON_LIME)
        return True

    start_t_sage = time.time()
    ref = None
    if args.type == 1:
        congruence_subgroup = SageGamma1(N)
        form = CuspForms(congruence_subgroup, args.k)
        ref = int(form.hecke_operator(p).trace())
    elif args.type == 0:
        congruence_subgroup = SageGamma0(N)
        form = CuspForms(congruence_subgroup, args.k)
        ref = int(form.hecke_operator(p).trace())
    else:
        if args.type == 2 and (p - 1) % N != 0:
            valid_primes = [pp for pp in primerange(2, 10**2) if (pp - 1) % N == 0]
            Logger.cprint(
                "Sage ref only valid for N | (p-1), try primes: " + ", ".join(map(str, valid_primes)),
                Colors.FAIL,
            )
            return True
        H = GammaH(N ** 2, [1 + N])
        form = CuspForms(H, args.k)
        ref = int(form.hecke_operator(p).trace())

    stop_t_sage = time.time()
    error = abs(trace - ref)
    clr = Colors.NEON_LIME if error == 0 else Colors.FAIL
    Logger.cprint(
        f"Result: trace={trace}, sage ref={ref}, error={error}, time={(stop_t - start_t):.9f}, sage time={(stop_t_sage - start_t_sage):.9f}",
        clr,
    )
    return True


def get_sage_ref(args, p: int, N: int, k: int):
    """Return (ref, msg, sage_time).

    - `ref` is the integer Sage trace if available, else None.
    - `msg` is a short explanatory message when ref is None (e.g. import error or congruence issue).
    - `sage_time` is the time spent in Sage (0 when not run).
    """
    try:
        from sage.all import (
            CuspForms,
            Gamma1 as SageGamma1,
            Gamma0 as SageGamma0,
            GammaH,
        )
    except Exception as e:
        return None, f"Unable to import Sage: {e}", 0.0

    t0 = time.time()
    ref = None
    # type 1 and 0 are straightforward
    if args.type == 1:
        congruence_subgroup = SageGamma1(N)
        form = CuspForms(congruence_subgroup, k)
        ref = int(form.hecke_operator(p).trace())
    elif args.type == 0:
        congruence_subgroup = SageGamma0(N)
        form = CuspForms(congruence_subgroup, k)
        ref = int(form.hecke_operator(p).trace())
    else:
        # For GammaH/Gamma2 we require N | (p-1) for a valid reference
        if args.type == 2 and (p - 1) % N != 0:
            valid_primes = [pp for pp in primerange(2, 10**2) if (pp - 1) % N == 0]
            msg = (
                "Sage ref only valid for N | (p-1); try primes: "
                + ", ".join(map(str, valid_primes))
            )
            return None, msg, 0.0
        H = GammaH(N ** 2, [1 + N])
        form = CuspForms(H, k)
        ref = int(form.hecke_operator(p).trace())

    t1 = time.time()
    return ref, None, (t1 - t0)
