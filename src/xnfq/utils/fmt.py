from ..arithmetic.common import factorize

def fmt_magnitude(n: int | float) -> str:
    """Fast order-of-magnitude label, e.g. '~10^6' or '~10^-3'."""
    import math
    if n == 0:
        return "0"
    exp = math.floor(math.log10(abs(n)))
    return f"~10^{exp}"

def fmt_invariants(inv: tuple, color: bool = True) -> str:
    a, b = fmt_factored(inv[0]), fmt_factored(inv[1])
    x = "\033[33mx\033[0m" if color else "x"
    return f"({a} {x} {b})"

def fmt_factored(n):
    if n == 0:
        return "0"
    if n > 10**20:
        return str("")
    factors = factorize(abs(n))
    if not factors:
        return str(n)
    s = "·".join(f"{p}^{e}" if e > 1 else str(p) for p, e in factors)
    return f"-{s}" if n < 0 else s
