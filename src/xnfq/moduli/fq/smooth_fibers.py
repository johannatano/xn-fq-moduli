from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from math import prod, sqrt, gcd
from typing import TYPE_CHECKING

from ...arithmetic.common import (
    factorize,
    kronecker,
    legendre,
    valuation as vl,
    divisors,
)
from ...arithmetic.forms import BinaryQuadraticForm
from ...arithmetic.function import Phi, phi
from ...arithmetic.quadratic import BicyclicGroup, LatticeTower, QuadraticOrderElement
from ...config import get_config, get_pari
from ...utils.logging import Colors, Logger
from ...utils.fmt import fmt_factored
from ..data import EigenForm, EigenFormRecord, FiberRecord, LevelStructureRecord
from ..modular_curve import CurveFiber

from ..level_structures import LevelStructure


@dataclass
class EigenValues:
    """Local eigenvalue data for one prime-power divisor of `N`."""
    l: int
    a: int
    split_type: int
    values: list[int] = field(default_factory=list)
    @property
    def modulus(self) -> int:
        return self.l**self.a
    def is_empty(self) -> bool:
        return len(self.values) == 0

@dataclass
class StableLinesData:
    """One Frobenius trace together with its admissible local eigenvalues."""
    trace: int
    base_form: BinaryQuadraticForm
    eigenvalues: dict[int, EigenValues] = field(default_factory=dict)
    global_eigenvalues: list[int] = field(default_factory=list)

    def compute_global_eigenvalues(self, N: int) -> None:
        """Combine the local eigenvalue sets into residues modulo ``N``."""
        combined = [0]
        modulus = 1
        for local_data in self.eigenvalues.values():
            next_combined = []
            local_modulus = local_data.modulus
            for current in combined:
                for local_value in local_data.values:
                    offset = (
                        (local_value - current)
                        * pow(modulus, -1, local_modulus)
                    ) % local_modulus
                    next_combined.append(current + modulus * offset)
            combined = next_combined
            modulus *= local_modulus
        self.global_eigenvalues = sorted({value % N for value in combined})

    def survives(self) -> bool:
        return all(not data.is_empty() for data in self.eigenvalues.values())


class WeilqFiber(CurveFiber):
    """A minimal wrapper around a Frobenius lattice tower and one trace sign."""

    def __init__(
        self,
        gamma: LevelStructure,
        t: int,
        frob_tower: LatticeTower,
        frob_lines: StableLinesData,
        mass: Fraction,
    ):
        super().__init__(gamma)
        self.t = t
        self.frob_tower = frob_tower
        pi_form = frob_lines.base_form
        self.frob = QuadraticOrderElement.from_norm_form(frob_tower.OK, pi_form)
        self.m0 = mass  # * (1 if self.base.level_structure.type > 0 else 2)
        self.frob_lines = frob_lines
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
        lam: int,
    ) -> LevelStructureRecord:
        """Build the conductor-`f` record for one global eigenvalue modulo `N`."""
        pi_local = self.frob_tower.order(f).embed_suborder(self.frob)
        pi_local_shifted = pi_local.shift(-lam)
        im = [1, 1]
        ker = [1, 1]
        num_lines = 1
        scalar = True
        for l, a in prime_powers:
            im_local = pi_local_shifted.ell_invariants(l)
            ker_local = pi_local_shifted.ell_kernel(l, a)
            for idx, exp in enumerate(im_local):
                im[idx] *= l**exp
            for idx, exp in enumerate(ker_local):
                ker[idx] *= l**exp
            num_lines *= self.stable_lines_count(l, a, pi_local_shifted)
            scalar = scalar and (vl(pi_local_shifted.v, l) >= a)

        return LevelStructureRecord(
            index=f,
            order=self.frob_tower.order(f),
            mass=self.local_mass(f),
            coords=pi_local_shifted.coords,
            im=tuple(im),
            ker=tuple(ker),
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
        if self.frob.v == 0:
            c_N = 1
        else:
            q = self.frob.norm
            #c_N = prod(p ** vl(self.frob.v, p) for p, _ in factorize(self.N) if q % p != 0) # N supported, DOES NOT CHANGE over lambdas
            c_N = prod(l**a for l, a in c_N_support.items())
            #print(f"c_N_2={c_N_2}, c_N={c_N}")

        for lam in self.frob_lines.global_eigenvalues:
            alpha = self.frob.shift(-lam)
            for d in divisors(c_N):
                n1 = gcd(alpha.u, self.frob.v // d, self.N) #TODO: check correct for (p,N) neq 1
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
        for lam in self.frob_lines.global_eigenvalues:
            alpha = self.frob.shift(-lam)
            level_records = [
                self.local_level_record(prime_powers, f, lam) for f in conductors
            ]
            eigen_records.append(
                EigenFormRecord(
                    base_form=self.frob.to_norm_form(),
                    eigenform=EigenForm(eigenvalue=lam, form=alpha.to_norm_form()),
                    level_records=tuple(level_records),
                )
            )
        return tuple(eigen_records)

    def snapshot(self) -> FiberRecord:
        eigen_records = self.get_structure()
        return FiberRecord(
            kind="weil",
            t=self.t,
            d=self.frob_tower.DK,
            total_count=self.count(),
            normalized_count=sum(
                record.num_lines
                for eigen_record in eigen_records
                for record in eigen_record.level_records
            ),
            total_mass=self.m0 * sum(self.local_mass(f) for f in self.frob_tower.conductor_divisors()),
            conductors=tuple(self.frob_tower.conductor_divisors()),
            coprime=tuple(self.frob_tower.split(self.N)[1]),
            eigen_records=eigen_records,
        )


def enum_weil_q(curve: "ModularCurveFq") -> list[WeilDatum]:
    """Enumerate the trace strata `(t, tower)` over `F_{p^n}`."""
    strata: list[WeilDatum] = []
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
                # only want one sign for t == 0
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

    def get_eigenvals(
        base_form: BinaryQuadraticForm,
        l: int,
        a: int,
    ) -> list[int]:
        modulus = l**a
        shifts: list[int] = []
        for lam in curve.level_structure.eigenvalues(modulus):
            # For a shifted monic form, divisibility by l^a is detected on the
            # constant term alone: C' = lam^2 - t*lam + q.
            shifted_c = base_form.A * lam * lam - base_form.B * lam + base_form.C
            if shifted_c % modulus != 0:
                continue
            shifts.append(lam)
        return shifts

    def get_eigenforms(
        t: int,
        signs: list[int],
    ) -> dict[int, StableLinesData]:
        eigenforms: dict[int, StableLinesData] = {}
        allowed_split_types = [0, 1] if curve.level_structure.scalar_only else [0, 1]
        D = t * t - 4 * q
        for s in signs:
            trace = s * t
            if not is_valid(D, trace):
                continue
            base_form = BinaryQuadraticForm(1, trace, q)
            eigenforms[trace] = StableLinesData(trace, base_form)
        for l, a in factorize(curve.N):
            split_type = kronecker(D, l)
            for trace_data in eigenforms.values():
                values = []
                if split_type in allowed_split_types:
                    values = get_eigenvals(trace_data.base_form, l, a)
                trace_data.eigenvalues[l] = EigenValues(l, a, split_type, values)
        for trace_data in eigenforms.values():
            trace_data.compute_global_eigenvalues(curve.N)
        return eigenforms

    def get_fibers_over_t(t: int, signs: list[int] = None) -> list[WeilDatum]:
        if signs is None:
            signs = t_signs(t)
        eigenforms = get_eigenforms(t, signs)
        base_form = BinaryQuadraticForm(1, t, q)
        surviving = [
            trace_data for trace_data in eigenforms.values() if trace_data.survives()
        ]
        if not surviving:
            return []
        # DK, f = base_form.to_D0_basis()
        # NOTE: We only construct ONE tower and reuse, since L(pi)=L(-pi) in tower structure, and it is expensive to find D0.
        DK, f = BinaryQuadraticForm._D0(base_form.discriminant, pari=pari)
        tower = LatticeTower(DK, f, base_form, exclude=p)
        
        
        return [(trace_data.trace, tower, trace_data) for trace_data in surviving]

    # in the case of Gamma1, we might speed up t^2 <= 4q enumeration by solving explicit residues for pm(q+1) % N
    if config.fast_trace and curve.level_structure.type == 1:
        for t, signs in trace_classes(HB).items():
            strata.extend(get_fibers_over_t(t, signs))
    else:
        for t in range(0, HB + 1):
            strata.extend(get_fibers_over_t(t))
    return strata
