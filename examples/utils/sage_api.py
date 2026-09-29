from __future__ import annotations

from sympy import primerange

from utils.logging import Logger, Colors
from xnfq.moduli.fq.data import TrFqTraceRecord
def tr_fq(args, p: int, n: int, N: int, k: int) -> tuple[TrFqTraceRecord, str]:
    rec = TrFqTraceRecord(0, 0, 0, 0)
    if n != 1:
        return rec, "Sage comparison only supports prime fields (n=1); skipping."
    try:
        from sage.all import (
            CuspForms,
            Gamma1 as SageGamma1,
            Gamma0 as SageGamma0,
            GammaH,
        )
    except Exception as e:
        return rec, f"Unable to import Sage: {e}"
    # type 1 and 0 are straightforward
    if args.type == 1:
        congruence_subgroup = SageGamma1(N)
        form = CuspForms(congruence_subgroup, k)
        rec = TrFqTraceRecord(val=int(form.hecke_operator(p).trace()))
    elif args.type == 0:
        congruence_subgroup = SageGamma0(N)
        form = CuspForms(congruence_subgroup, k)
        rec = TrFqTraceRecord(val=int(form.hecke_operator(p).trace()))
    else:
        # For GammaH/Gamma2 we require N | (p-1) for a valid reference
        if args.type == 2 and (p - 1) % N != 0:
            valid_primes = [pp for pp in primerange(2, 10**2) if (pp - 1) % N == 0]
            msg = (
                "Sage ref only valid for N | (p-1); try primes: "
                + ", ".join(map(str, valid_primes))
            )
            return rec, msg
        H = GammaH(N ** 2, [1 + N])
        form = CuspForms(H, k)
        rec = TrFqTraceRecord(val=int(form.hecke_operator(p).trace()))
    return rec, None
