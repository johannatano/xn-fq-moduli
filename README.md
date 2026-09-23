# xn_fq_moduli

`xn_fq_moduli` is a computational framework for studying level-
`N` modular curves and their point-counting data over finite fields. The
central design goal is to keep the moduli problem, the arithmetic of the
finite field, and the internal structure of each geometric fiber separate.
This makes the package useful both as a counting program and as an
experimental companion to a theoretical description of the strata of
`X_N`.

## Run locally

The project can be run directly from this repository. Requirements are Git and
Python 3.11 or newer.

Clone the repository, create an isolated environment, and install this local
checkout together with the example dependencies:

```bash
git clone YOUR_REPOSITORY_URL
cd xn-fq-moduli
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[examples]"
```

Run the examples directly from the repository:

```bash
python examples/xnfq_count.py --type 1 --N 11 --p 5 --n 1
python examples/xnfq_interactive.py
```

The interactive plotting example requires a graphical environment. The
editable install points Python at the checked-out source, so changes to the
repository are immediately available.

The optional PARI-backed path can be added with:

```bash
python -m pip install -e ".[pari]"
```

Use the PARI path by enabling the session option before computing:

```python
from xnfq.config import configure_session

configure_session(use_pari=True)
```

## Mathematical viewpoint

The package treats a modular curve as a level structure together with a
base field. The level structure is specified independently of the field and
is represented by an object such as

```python
from xnfq.moduli.modular_curve import X, X0, X1

X0(N)       # Gamma_0(N)
X1(N)       # Gamma_1(N)
X(N)        # full level Gamma(N)
```

The common abstraction is `LevelStructure`. It records the level `N`, the
type of level problem, the admissible eigenvalue shifts, and the weighting
factor that converts the underlying geometric count into the count for the
chosen moduli problem. The concrete level types currently used by the
package are:

- `Gamma0(N)`: a cyclic subgroup of order `N`; eigenvalue shifts are
  generally a family of residues modulo `N`.
- `Gamma1(N)`: a chosen point of exact order `N`; the distinguished shift is
  `lambda = 1`.
- `Gamma(N)`: full level structure, with the scalar-compatibility condition
  built into the level object and its weighting.

The point of this interface is that the geometric and arithmetic machinery
does not need a separate implementation for every Gamma type. A level
problem supplies the admissible local or global shifts and the relevant
weights, while the fiber code performs the common calculations.

For a finite field, the same level problem is base-changed by

```python
curve = X0(N).over(p, n)
```

which studies the curve over `F_q`, where `q = p^n`. The finite-field curve
then exposes the decomposition

```python
curve.smooth_fibers()
curve.cusps()
curve.fibers()
```

and can aggregate the corresponding contributions to the point count or to
Frobenius-symmetric-power traces.

## The first stratification: smooth and cuspidal fibers

The compactified modular curve naturally separates into its non-cuspidal
and cuspidal parts. The package mirrors this geometry directly.

### Smooth fibers: Frobenius or Weil strata

The smooth part is enumerated by Frobenius trace data. For a fixed finite
field `F_q`, the relevant trace values satisfy the usual Hasse bound

```text
|t| <= 2 sqrt(q),
```

together with the congruence and level-survival conditions imposed by the
chosen level problem. A surviving trace class determines a Frobenius
polynomial, represented here by a binary quadratic form such as

```text
x^2 - t*x + q.
```

The discriminant of this polynomial determines the relevant imaginary
quadratic field or order. The associated `LatticeTower` organizes the
orders containing the Frobenius endomorphism as the conductor varies. Thus
the basic smooth datum is a Weil/Frobenius fiber

```text
(t, lattice tower, eigenvalue data).
```

The tower carries the arithmetic information needed to evaluate class
masses, embed Frobenius into orders of different conductors, and compute the
local structure of the quotient by a shifted Frobenius element.

### Cuspidal fibers: Neron polygons

The cuspidal part is represented by degenerate fibers, modeled by Neron
`d`-gons. A cusp datum has the form

```text
(t, NeronDgonFq(d, q, t), admissible eigenvalues).
```

Here `t` records the split or nonsplit Frobenius sign and `d` is the
geometric polygon parameter. The complementary integer

```text
e = N / d
```

controls the other cyclic direction in the level structure. For each
candidate eigenvalue `lambda`, the compatibility conditions determine the
possible `d` values through the divisibility data

```text
d1 = gcd(lambda - t*q, N)
d2 = gcd(lambda - t, N).
```

The cusp enumerator uses these conditions to discover the compatible
`(t, d)` strata and accumulates all admissible eigenvalues on the same
stratum. In particular, the output contains one cusp fiber per geometric
stratum rather than one fiber for every `(lambda, d)` incidence.

## The second stratification: arithmetic structure inside a fiber

The first decomposition says which geometric fiber is being considered.
The next decomposition records how the level structure interacts with the
arithmetic object defining that fiber. This is where the package is intended
to expose structure rather than only produce a final integer.

### Global eigenvalues and local eigenvalues

Let

```text
N = product over l of l^a.
```

For a fixed smooth Frobenius class, the package first computes admissible
eigenvalue shifts at each prime power `l^a`. These are stored as local
`EigenValues` records. The global eigenvalues are the residues

```text
lambda mod N
```

whose reductions are admissible at every prime power dividing `N`. Since the
prime powers are pairwise coprime, the Chinese remainder theorem identifies
the global set with the product of the local sets:

```text
Lambda_N ~= product over l of Lambda_(l^a).
```

The `EigenTrace` object therefore keeps both levels of information:

- local eigenvalue sets, together with the prime, exponent, and splitting
  type;
- the CRT-combined global eigenvalue set for the whole level `N`.

This distinction matters conceptually. Local eigenvalues explain the
prime-by-prime constraints, while a global `lambda` is the actual shift used
to form

```text
alpha_lambda = pi - lambda
```

globally. The global set is computed once for a fixed Frobenius/isogeny
class and is then reused by all subsequent conductor calculations.

### Smooth fibers: conductor and order strata

Fix a smooth Frobenius class and one global eigenvalue `lambda`. The shifted
endomorphism `alpha_lambda` is then examined through the conductor tower.
For each relevant conductor `f`, the code embeds the Frobenius element into
the order of conductor `f` and computes the local quotient structure at all
prime powers dividing `N`.

The natural dependency is therefore

```text
Frobenius class -> global lambda -> conductor f -> local level structure.
```

The conductor is not a second way of discovering the global eigenvalues. It
is a refinement of a fixed shifted Frobenius object. This is why the
serializable smooth-fiber records are organized as an eigenvalue record
containing a collection of conductor-level records. Each conductor record
can retain data such as:

- the order containing the shifted Frobenius;
- its class-group mass;
- coordinates of the shifted element;
- image and kernel invariant factors at level `N`;
- the number of stable lines;
- whether the action is scalar at the relevant prime powers.

This organization supports both computation and exposition. The total
smooth contribution is obtained by summing over global eigenvalues and over
the conductor tower, but the intermediate records preserve the arithmetic
reason for each contribution.

### Cuspidal fibers: `d`-gon strata and eigenvalue incidences

For cusps, `d` is not merely an auxiliary order parameter. It identifies the
geometric boundary stratum itself. The natural organization is therefore

```text
cuspidal sign t -> d-gon parameter d -> admissible global lambdas.
```

The enumerator discovers a compatible `(t, d)` pair while scanning the
level-compatibility conditions and appends each compatible eigenvalue to the
single record for that pair. Repeated discoveries of the same `(t, d)` do
not create duplicate fibers. The resulting `CuspFqFiber` can count the
stable level structures using the size of its admissible eigenvalue list and
the Neron polygon invariants.

Thus the smooth and cuspidal parts have deliberately different inner
stratifications:

```text
smooth:  Frobenius class -> global lambda -> conductor/order data
cusps:   (t, d)-gon stratum -> global lambda incidences
```

This asymmetry reflects the geometry. A smooth fiber is controlled by a
Frobenius endomorphism and its orders; a cusp fiber is controlled first by
the boundary polygon type.

## Fiber records and structural reports

Every fiber implements the same basic interface: it can be counted and it
can be converted into a serializable `FiberRecord`. A complete structure
report is obtained from

```python
report = curve.get_structure()
```

The report retains the level problem and the list of smooth and cusp fiber
records. For smooth fibers, a record includes the trace, discriminant,
conductors, coprime conductor contribution, total mass, and the nested
eigenvalue/conductor records. For cusp fibers, it includes the sign, the
`d`-gon parameter, and the cusp contribution.

This is useful for more than debugging. It gives a basis for tables and
experiments in which one wants to ask questions such as:

- Which trace classes survive a given level problem?
- How many global eigenvalues occur in a fixed smooth class?
- Which conductors contribute for a particular shifted Frobenius element?
- How do stable-line counts change between different orders?
- Which `d`-gon strata carry the cuspidal mass for `Gamma_0(N)`,
  `Gamma_1(N)`, or full level?

The record-oriented interface also keeps the expensive arithmetic objects
separate from the final presentation layer. A computation can construct the
tower once, reuse it for all signs or eigenvalues in the same class, and
then materialize a compact report.

## Frobenius symmetric-power traces

The trace parameter `t` is retained on every fiber because it is precisely
the variable needed for later Frobenius symmetric-power calculations. If
`h_k(t, q)` denotes the trace polynomial for the `k`th symmetric power,
then the smooth contribution has the schematic form

```text
sum over smooth fibers  h_k(t, q) * fiber_count,
```

while the cusp contribution is weighted by the corresponding cusp power of
the sign parameter. The finite-field curve exposes this through

```python
curve.tr_frob_symk(k)
```

The important architectural choice is that the fiber enumeration does not
need to be repeated for every `k`. The expensive work produces a stratified
collection of fibers carrying their `t` values and arithmetic records. Each
later symmetric-power calculation only changes the trace polynomial applied
to those records.

In this sense, the package separates two stages:

1. **Geometric and arithmetic enumeration:** determine the smooth Frobenius
   fibers, cusp `d`-gon fibers, admissible global eigenvalues, and conductor
   or local structure data.
2. **Cohomological or representation-theoretic aggregation:** apply the
   appropriate polynomial in `t` to the stored fiber counts to obtain
   Frobenius trace expressions.

## Design summary

The package is organized around the following chain:

```text
level problem Gamma at N
	|
	v
base change to F_q
	|
	+----------------------+----------------------+
	|                                             |
	v                                             v
smooth Frobenius fibers                         cusp d-gon fibers
	|                                             |
	v                                             v
global lambda set                               (t, d) stratum
	|                                             |
	v                                             v
conductor/order refinements                    lambda incidences
	|                                             |
	+----------------------+----------------------+
			       v
		      fiber records and counts
			       |
			       v
		    Frobenius symmetric powers
```

The main principle is modularity of the moduli problem and explicitness of
the strata. Changing from `Gamma_0` to `Gamma_1` or full level changes the
admissibility and weighting rules, but not the overall fiber architecture.
Changing the base field changes the available trace classes and arithmetic
masses, but not the representation of a fiber. Finally, refining a fiber by
global eigenvalue, conductor, or `d`-gon type reveals the arithmetic content
behind the aggregate point count and makes the same data reusable for later
trace computations.