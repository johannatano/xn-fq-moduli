from __future__ import annotations
from math import gcd
from ..config import get_pari

""" Number theory basics """
def kronecker(D, n):
    """Kronecker symbol (D/n) for n >= 1."""
    if n == 0:
        return 1 if D in (1, -1) else 0
    if gcd(D, n) != 1:
        return 0
    r = 1
    while n % 2 == 0:  # (D/2) rule
        n //= 2
        if D % 8 in (3, 5):
            r = -r
    D %= n  # Jacobi by reciprocity
    while D:
        while D % 2 == 0:
            D //= 2
            if n % 8 in (3, 5):
                r = -r
        D, n = n, D
        if D % 4 == 3 and n % 4 == 3:
            r = -r
        D %= n
    return r if n == 1 else 0


def legendre(a, p):
    return kronecker(a, p)

""" Factorizations """
def factorize(n):
    """Return list of (prime, exponent) pairs for n > 1."""
    pari = get_pari()
    if pari:
        return [
            (int(prime), int(exponent))
            for prime, exponent in pari.factorint(n).python()
        ]
        #from sympy import factorint
        #return list(factorint(n).items())
    factors = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            e = 0
            while n % d == 0:
                n //= d
                e += 1
            factors.append((d, e))
        d += 1
    if n > 1:
        factors.append((n, 1))
    return factors

def divisors(n):
    divs = []
    i = 1
    while i * i <= n:
        if n % i == 0:
            divs.append(i)
            if i != n // i:
                divs.append(n // i)
        i += 1
    return sorted(divs)


""" p-adic valuation, supports n = 0 returns inf """
def valuation(n, p):
    if n == 0:
        return float("inf")
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v
