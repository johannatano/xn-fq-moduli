from .common import factorize, kronecker, divisors

""" Arithmetic functions and Dirichlet characters """
class ArithmeticFunction:
    def __call__(self, n): ...

class Multiplicative(ArithmeticFunction):
    pass

class CompletelyMultiplicative(Multiplicative):
    pass

class DirichletCharacter(CompletelyMultiplicative):
    def __init__(self, modulus):
        self.modulus = modulus

class phi(Multiplicative):
    """f * prod_{p|f} (1 - chi(p)/p) -- the chi-twisted totient."""

    def __init__(self, chi: DirichletCharacter | int = 1):
        if chi in (1, -1):
            self.chi = lambda _n, value=chi: value
        else:
            self.chi = chi

    '''def __call__(self, f):
        r = Fraction(f)
        for p, _ in factorize(f):
            r *= Fraction(p - self.chi(p), p)
        return r'''

    def __call__(self, f):
        if f == 0:
            return 0
        value = 1
        for p, a in factorize(f):
            value *= p ** (a - 1) * (p - self.chi(p))
        return value

class phi_ell(Multiplicative):
    def __init__(self, chi: DirichletCharacter | int = 1):
        if chi in (1, -1):
            self.chi = lambda _n, value=chi: value
        else:
            self.chi = chi

    def __call__(self, ell: int, a: int):
        if a == 0:
            return 1
        return ell ** (a - 1) * (ell - self.chi(ell))

class Mobius(Multiplicative):
    def __call__(self, n):
        if n == 1:
            return 1
        r, p = 1, 2
        while p * p <= n:
            if n % p == 0:
                n //= p
                if n % p == 0:
                    return 0
                r = -r
            p += 1
        return -r if n > 1 else r

class KroneckerCharacter(DirichletCharacter):
    def __init__(self, D):
        super().__init__(abs(D))
        self.D = D

    def __call__(self, n):
        return kronecker(self.D, n)


class IdealCount(Multiplicative):
    def __init__(self, chi: DirichletCharacter):
        self.chi = chi

    def __call__(self, n):
        return sum(self.chi(d) for d in range(1, n + 1) if n % d == 0)


def dirichlet_series(a: ArithmeticFunction, s, N=10000):
    res = 0
    for n in range(1, N):
        term = a(n) * n ** (-s)
        # print(f"term={Fraction(a(n), n ** (s))}")
        res += term
    return res

def Phi_ell(chi: DirichletCharacter | int, ell: int, a: int):
    if a < 0:
        raise ValueError("a must be non-negative")
    chi_ell = chi if chi in (1, -1) else chi(ell)
    geometric_sum = (ell**a - 1) // (ell - 1) if a > 0 else 0
    return 1 + (ell - chi_ell) * geometric_sum


def Phi(chi: DirichletCharacter, f: int):
    return sum(phi(chi)(e) for e in divisors(f))
