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

class KroneckerCharacter(DirichletCharacter):
    def __init__(self, D):
        super().__init__(abs(D))
        self.D = D

    def __call__(self, n):
        return kronecker(self.D, n)


def Phi(chi: DirichletCharacter, f: int):
    return sum(phi(chi)(e) for e in divisors(f))
