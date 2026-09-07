from __future__ import annotations
from math import gcd
from ..config import get_config, get_pari

pari = get_pari()
""" Number theory basics """
def euler_phi(n):
    result = n
    for p, _ in factorize(n):
        result -= result // p
    return result

def jordan_totient(n):
    """Number of elements of exact order n in (Z/n x Z/n) — counts pairs of generators."""
    result = n**2
    for p, _ in factorize(n):
        result *= p**2 - 1
        result //= p**2
    return result

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
    if p == 2:
        return 0 if a % 2 == 0 else (1 if a % 8 in (1, 7) else -1)
    a = a % p
    if a == 0:
        return 0
    r = pow(a, (p - 1) // 2, p)
    return -1 if r == p - 1 else 1


""" Factorizations """
def coprime_part(n: int, m: int) -> int:
    """Largest divisor of n coprime to m — strips all primes of m from n."""
    import math
    g = math.gcd(n, m)
    while g > 1:
        n //= g
        g = math.gcd(n, m)
    return n

def factorize(n):
    """Return list of (prime, exponent) pairs for n > 1."""
    #if pari:
    #return list(pari.factorint(n).items())
    #    from sympy import factorint
    #    return list(factorint(n).items())
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


""" integer sign and equivalence functions """
def sgn(x):
    return 1 if x > 0 else -1 if x < 0 else 0

def equiv(a, b):
    val = a % b
    return -1 if val == b - 1 else val

""" p-adic valuation, supports n = 0 returns inf """
def valuation(n, p):
    if n == 0:
        return float("inf")
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v

""" other arithmetic functions """
def convolve(a, b):
    def c(n):
        return sum(a(d) * b(n // d) for d in divisors(n))
    return c
