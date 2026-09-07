from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from math import prod, sqrt
from typing import TYPE_CHECKING

from ...arithmetic.common import factorize, kronecker, legendre, valuation as vl
from ...arithmetic.forms import BinaryQuadraticForm
from ...arithmetic.function import Phi, phi
from ...arithmetic.quadratic import BicyclicGroup, LatticeTower, QuadraticOrderElement
from ...config import get_config, get_pari
from ...utils.logging import Colors
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
class EigenTrace:
    """One Frobenius trace together with its admissible local eigenvalues."""
    trace: int
    base_form: BinaryQuadraticForm
    eigenvalues: dict[int, EigenValues] = field(default_factory=dict)

    def survives(self) -> bool:
        return all(not data.is_empty() for data in self.eigenvalues.values())


class WeilqFiber(CurveFiber):
    """A minimal wrapper around a Frobenius lattice tower and one trace sign."""

    def __init__(
        self,
        gamma: LevelStructure,
        t: int,
        frob_tower: LatticeTower,
        eigen_trace: EigenTrace,
        mass: Fraction,
    ):
        super().__init__(gamma)
        self.t = t
        self.frob_tower = frob_tower
        # sign = -1 if t < 0 else 1
        pi_form = eigen_trace.base_form
        self.frob = QuadraticOrderElement.from_norm_form(frob_tower.OK, pi_form)


        self.m0 = mass  # * (1 if self.base.level_structure.type > 0 else 2)
        self.eigen_trace = eigen_trace
        self.chi_K = frob_tower.chi_K

    def stable_lines_count(
        self, l: int, a: int, alpha: QuadraticOrderElement, scalar_only=False
    ) -> int:
        """Count stable lines in the local quotient defined by `alpha`."""
        # e1 = min(vl(alpha.u, l), vl(alpha.v, l))
        # e2 = vl(alpha.norm, l) - e1
        e1, e2 = alpha.ell_invariants(l)
        if e2 < a:
            return 0
        if e1 < a:
            return l**e1
        else:
            return phi(-1)(l**a)

    def local_structure_count(self, l: int, a: int, f: int = 0) -> int:
        """Count local level structures above the conductor slice `f`."""
        eigen_data = self.eigen_trace.eigenvalues[l]
        num_eigen = len(eigen_data.values)
        if num_eigen == 0:
            return 0
        pi_local = self.frob_tower.order(f).embed_suborder(self.frob)
        # some early exits
        if self.gamma.type != 1 and vl(pi_local.v, l) >= a:
            total_val = phi(-1)(l**a)  # we reach maximum valuation
        elif self.gamma.type == 2:
            total_val = 0  # we only want max
        else:  # we have partial inclusion, we have to check each eigenvalue
            total_val = 0
            for lam in eigen_data.values:
                pi_local_shifted = pi_local.shift(-lam)
                total_val += self.stable_lines_count(l, a, pi_local_shifted)
        return total_val

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
            #coords_pi_basis=pi_local_shifted.coords_in_basis(self.frob),
            im=tuple(im),
            ker=tuple(ker),
            num_lines=num_lines,
            scalar=scalar,
        )

    def global_structure_count(self) -> int:
        """Global multiplicity contributed by the chosen level problem."""
        if self.gamma.type == 1:
            return phi(1)(self.N)
        if self.gamma.type == 2:
            return self.N * phi(1)(self.N)
        return 1

    def local_mass(self, f: int) -> int:
        return self.frob_tower.class_size(f)

    def tower_sum(self, l: int, a: int) -> int:
        """Aggregate the local counts across the conductor tower at `l^a`."""
        k = self.frob_tower.conductors.get(l, 0)
        return sum(
            self.local_mass(l**i) * self.local_structure_count(l, a, l**i)
            for i in range(0, k + 1)
        )

    """def local_structure(self, f: int) -> BicyclicGroup:
        return self.frob_tower.quotient_group(self.frob_1, f)"""
    """def tower_quotients(
        self, element: QuadraticOrderElement
    ) -> tuple[tuple[int, BicyclicGroup, int], ...]:
        return tuple(
            (f, self.local_structure(f), self.local_mass(f))
            for f in self.frob_tower.conductor_divisors()
        )"""

    def count(self) -> Fraction:
        """Compute the weighted point-count contribution of this smooth fiber."""

        lattice_sum = (
            prod(self.tower_sum(l, a) for l, a in factorize(self.N))
            * self.global_structure_count()
            * Phi(
                self.chi_K,
                self.frob_tower.max_coprime_conductor(self.N),
            )
        )

        return lattice_sum * self.m0

    def get_eigen_structure(self) -> tuple[EigenFormRecord, ...]:
        eigen_records: list[EigenFormRecord] = []

        prime_powers = tuple(factorize(self.N))
        conductors = tuple(self.frob_tower.supported_divisors(self.N))
        conductors = tuple(self.frob_tower.conductor_divisors())
        # coprime = tuple(self.frob_tower.split(self.N)[1])

        for lam in self.gamma.eigenvalues(self.N):
            is_valid = True
            for l, a in prime_powers:
                modulus = l**a
                residue = lam % modulus
                eigendata = self.eigen_trace.eigenvalues[l]
                if residue not in eigendata.values:
                    is_valid = False
                    break
            if not is_valid:
                continue

            alpha = self.frob.shift(-lam)
            level_records = [
                self.local_level_record(prime_powers, f, lam) for f in conductors
            ]

            eigen_records.append(
                EigenFormRecord(
                    eigenform=EigenForm(eigenvalue=lam, form=alpha.to_norm_form()),
                    level_records=tuple(level_records),
                )
            )

        return tuple(eigen_records)

    def snapshot(self) -> FiberRecord:
        eigen_records = self.get_eigen_structure()
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
    ) -> dict[int, EigenTrace]:
        eigenforms: dict[int, EigenTrace] = {}
        allowed_split_types = [0] if curve.level_structure.scalar_only else [0, 1]
        D = t * t - 4 * q
        for s in signs:
            trace = s * t
            if not is_valid(D, trace):
                continue
            base_form = BinaryQuadraticForm(1, trace, q)
            eigenforms[trace] = EigenTrace(trace, base_form)
        for l, a in factorize(curve.N):
            split_type = kronecker(D, l)
            for trace_data in eigenforms.values():
                values = []
                if split_type in allowed_split_types:
                    values = get_eigenvals(trace_data.base_form, l, a)
                trace_data.eigenvalues[l] = EigenValues(l, a, split_type, values)
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
