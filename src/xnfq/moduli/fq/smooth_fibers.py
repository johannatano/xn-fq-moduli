from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from math import prod, sqrt, gcd

from ...arithmetic.common import (
    factorize,
    kronecker,
    legendre,
    valuation as vl,
    divisors,
)
from ...arithmetic.forms import BinaryQuadraticForm
from ...arithmetic.function import Phi, phi
from ...arithmetic.quadratic import LatticeTower, QuadraticOrderElement
from ...config import get_config, get_pari
from .data import EigenForm, EigenFormRecord, LevelStructureRecord, SmoothFiberRecordFq
from ..modular_curve import CurveFiber

from ..level_structures import LevelStructure

class SmoothFiberFq(CurveFiber):
    """A minimal wrapper around a Frobenius lattice tower and one trace sign."""

    def __init__(
        self,
        gamma: LevelStructure,
        t: int,
        frob_tower: LatticeTower,
        eigen_forms: list[EigenForm],
        mass: Fraction,
    ):
        super().__init__(gamma)
        self.t = t
        self.frob_tower = frob_tower
        self.m0 = mass  # * (1 if self.base.level_structure.type > 0 else 2)
        self.eigen_forms = eigen_forms
        self.chi_K = frob_tower.chi_K

    def stable_lines_count(
        self, l: int, a: int, alpha: QuadraticOrderElement, scalar_only=False
    ) -> int:
        """Count stable lines in the local quotient defined by `alpha`."""
        e1, e2 = alpha.ell_invariants(l)
        if e2 < a:
            return 0
        if e1 < a:
            return l**e1 if not scalar_only else 0
        else:
            return phi(-1)(l**a)

    def local_level_record(
        self,
        prime_powers: tuple[tuple[int, int], ...],
        f: int,
        alpha: QuadraticOrderElement,
    ) -> LevelStructureRecord:
        """Build the conductor-`f` record for one global eigenvalue modulo `N`."""
        alpha_local = self.frob_tower.order(f).embed_suborder(alpha)
        """Build the conductor-`f` record for one global eigenform modulo `N`."""
        ker = [1, 1]
        num_lines = 1
        scalar = True
        for l, a in prime_powers:
            ker_local = alpha_local.ell_kernel(l, a)
            for idx, exp in enumerate(ker_local):
                ker[idx] *= l**exp
            num_lines *= self.stable_lines_count(l, a, alpha_local)
            scalar = scalar and (vl(alpha_local.v, l) >= a)

        return LevelStructureRecord(
            index=f,
            order=self.frob_tower.order(f),
            mass=self.local_mass(f),
            coords=alpha_local.coords,
            inv=tuple(ker),
            num_lines=num_lines,
            scalar=scalar,
        )

    def local_mass(self, f: int) -> int:
        return self.frob_tower.class_size(f)

    def count(self) -> Fraction:
        count = 0

        c_N_support, c_coprime = self.frob_tower.split(self.N)

        coprime = Phi(
            self.chi_K,
            prod(l**a for l, a in c_coprime.items()),
        )

        if self.frob_tower.DK == 0:
            c_N = 1
        else:
            c_N = prod(l**a for l, a in c_N_support.items())

        for eigen_form in self.eigen_forms:
            alpha = QuadraticOrderElement.from_norm_form(self.frob_tower.OK, eigen_form.form)
            for d in divisors(c_N):
                e = self.frob_tower.f_max // d  # TODO: check correct for (p,N) neq 1, ie might use c_N here
                n1 = gcd(alpha.u, e, self.N)  
                n2 = gcd(alpha.norm // n1, self.N)
                if n2 < self.N or (n1 < self.N and self.gamma.type == 2):
                    # NOTE: we still need this, since the exp might drop below in the conductor tower even if N | norm(alpha)
                    continue
                count += (
                    Fraction(phi(-1)(self.N), phi(-1)(self.N // n1))
                    * self.local_mass(d)
                    * coprime
                )
        return count * self.m0 * self.gamma.weight()

    def get_structure(self) -> tuple[EigenFormRecord, ...]:
        eigen_records: list[EigenFormRecord] = []
        prime_powers = tuple(factorize(self.N))
        conductors = tuple(self.frob_tower.conductor_divisors())

        for eigen_form in self.eigen_forms:
            alpha = QuadraticOrderElement.from_norm_form(self.frob_tower.OK, eigen_form.form)
            level_records = [
                self.local_level_record(prime_powers, f, alpha) for f in conductors
            ]
            eigen_records.append(
                EigenFormRecord(
                    eigenform=eigen_form,
                    levels=tuple(level_records),
                )
            )
        return tuple(eigen_records)

    def snapshot(self) -> SmoothFiberRecordFq:
        eigen_records = self.get_structure()
        return SmoothFiberRecordFq(
            t=self.t,
            D_K=self.frob_tower.DK,
            count=self.count(),
            eigen_records=eigen_records,
        )


def enum_weil_q(
    curve: "ModularCurveFq",
) -> list[tuple[int, LatticeTower, list[EigenForm]]]:
    """Enumerate the trace strata `(t, tower)` over `F_{p^n}`."""
    strata: list[tuple[int, LatticeTower, list[EigenForm]]] = []
    p, q = curve.p, curve.q
    HB = int(2 * sqrt(q))
    config = get_config()
    pari = get_pari()

    def trace_class_residues() -> tuple[tuple[int, int], ...]:
        return (((q + 1) % curve.N, 1), ((-(q + 1)) % curve.N, -1))

    def has_level(t: int) -> bool:
        return (q + 1 - t) % curve.N == 0

    def trace_classes(HB_local: int | None = None) -> dict[int, list[int]]:
        if HB_local is None:
            HB_local = HB
        classes: dict[int, list[int]] = {}
        for residue, sign in trace_class_residues():
            for t in range(residue, HB_local + 1, curve.N):
                if t == 0 and sign == -1:
                    continue
                trace = sign * t
                if has_level(trace):
                    classes.setdefault(t, []).append(sign)
        return classes

    def t_signs(t: int) -> tuple[int]:
        return (1, -1) if t != 0 else (0,)

    def is_valid(D: int, t: int) -> bool:
        if t % p != 0:
            return True
        return legendre(D, p) != -1

    def get_frob_data(
        t: int,
        signs: list[int],
    ) -> list[tuple[int, list[EigenForm]]]:
        data: list[tuple[int, list[EigenForm]]] = []
        D = t * t - 4 * q
        for s in signs:
            trace = s * t
            if not is_valid(D, trace):
                continue
            lambdas = range(curve.N) if curve.level_structure.type == 0 else [1]
            base_form = BinaryQuadraticForm(1, trace, q)
            eigen_forms: list[EigenForm] = []
            for lam in lambdas:
                norm = base_form.A * lam * lam - base_form.B * lam + base_form.C
                if norm % curve.N != 0:
                    continue
                eigen_forms.append(EigenForm(lam, base_form.shift(-lam)))

            if not eigen_forms:
                continue
            data.append((trace, eigen_forms))
        return data

    def get_fibers_over_t(t: int, signs: list[int] = None) -> list[WeilDatum]:

        if signs is None:
            signs = t_signs(t)
        data = get_frob_data(t, signs)
        
        if not data:
            return []
        # NOTE: We only construct ONE tower and reuse, since L(pi)=L(-pi) in tower structure, and it is expensive to find D0.
        base_form = BinaryQuadraticForm(1, t, q)
        DK, f = BinaryQuadraticForm._D0(base_form.discriminant, pari=pari)
        tower = LatticeTower(DK, f, base_form, exclude=p)
        return [(trace, tower, eigen_forms) for trace, eigen_forms in data]

    # in the case of Gamma1, we might speed up t^2 <= 4q enumeration by solving explicit residues for pm(q+1) % N
    if config.fast_trace and curve.level_structure.type == 1:
        for t, signs in trace_classes(HB).items():
            strata.extend(get_fibers_over_t(t, signs))
    else:
        for t in range(0, HB + 1):
            strata.extend(get_fibers_over_t(t))
    return strata
