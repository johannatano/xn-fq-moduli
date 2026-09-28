# xn_fq_moduli

Python software for counting points and studying the structure of modular curves over finite fields.

## Contents
- [01. Overview](#01-overview)
- [02. Installation](#02-installation)
- [03. Documentation](#03-documentation)
- [04. Examples](#04-examples)
- [05. TODO](#05-todo)

## 01. Overview
We study $\mathbb F_q$-rational points $(E,\gamma)$ on modular curves, where $E$ is a generalized elliptic curve and $\gamma$ is a level-$N$ structure of unspecified type. By sending each point
$$
(E,\gamma)\;\mapsto\;E\;\mapsto\;t, \qquad t = a_p(E)
$$
we naturally stratify the point count by $t$ allowing efficient computation of traces of modular forms.

The fibers $\gamma$ over a fixed $E$ fall into two primary types — smooth fibers and cusp fibers. Over $\overline{\mathbb F}_q$ smooth fibers are modelled by level structures on $\mathbb Z/N\mathbb Z\times\mathbb Z/N\mathbb Z$, while cusp fibers are modelled via Néron $d$-gons $\mathbb F_q^{\times}\times\mathbb Z/d\mathbb Z$.

The $\mathbb F_q$-rational points in a fiber correspond to Frobenius-stable eigenforms, determined by admissible eigenvalues $\lambda$ depending on the chosen level structure type. For a smooth fiber the eigenform is modelled by the binary quadratic form $(1,t,q)$ shifted by $-\lambda$, giving the algebraic integer $\pi-\lambda$ which satisfies
$N\mid\mathrm{Norm}(\pi-\lambda)$. (When $\pi\in\mathbb Z$ this degenerates to an integral lattice in $\mathbb Z^2$). Cuspidal fibers are described by congruence conditions on Néron polygons.


## 02. Installation
From the repository root:
```bash
python -m pip install -e ".[examples]"
```
The optional PARI support can be installed with:

```bash
python -m pip install -e ".[pari]"
```
PARI **significantly** speeds up computations for large $q$. The project does not
yet use a cached database for class numbers (TODO). After installing the optional
PARI support, pass `--use-pari` to an example:

```bash
python examples/xnfq_count.py --use-pari
```

## 03. Documentation
### 3.1 Level Structure
```python
from xnfq.moduli import Gamma0, Gamma1, Gamma2, Gamma
N = 11
# Gamma0: subgroup of upper-triangular/cyclic N-torsion (not scalar-only)
Gamma0(N)
# Gamma1: choice of a point of exact order N (stabilizer subgroup, not scalar-only)
Gamma1(N)
# Gamma2: full N-torsion basis; marked as scalar-only in the code
Gamma2(N)
# convenience factory
gamma = Gamma(N, type=1)  # returns Gamma1(N)
gamma = Gamma(N)          # defaults to full level (Gamma2)
```
### 3.2 Modular Curve
```python
from xnfq.moduli.level_structures import Gamma
from xnfq.moduli.modular_curve import X, X0, X1
N = 13
# General constructor
mod_curve = ModularCurve(Gamma(N, type=1))
# Set base field
mod_curve.over(p=5, n=2)

# Convenience: set base field using full field size `q = p**n`
# (equivalent to `over(p,n)`)
mod_curve.F(5**2)  # same as mod_curve.over(5, 2)

# Convenience wrappers
X0(N) # ModularCurve(Gamma(N, type=0)) 
X1(N) # ModularCurve(Gamma(N, type=1)) 
X(N) # ModularCurve(Gamma(N, type=2))

```
### 3.3 Methods
#### 3.3.1 Count
The count method returns the total number of $\mathbb{F}_q$-rational points over each fiber type. 
```python
Y1, Cusp1 = X1(11).F(5**4).count()
```
Both have the following universal structure on the internal strata over $t$

$$
\mathrm{eigen\_count}(t)=\sum_{\lambda}\;\sum_{L\in\mathrm{Levels}(\lambda)}\;\mathrm{num\_lines}(\lambda,L)\cdot\mathrm{mass}(L).
$$

local levels are either sublattices indexed by conductor, with torsion invariant $\mathcal L_f/(\pi-\lambda)\mathcal L_f$, or Neron $d$-gon with invariants $\mu_{N/d}\times\mathbb Z/d\mathbb Z$.

The final fiber count of $(E, \gamma)$ over $t$ is then obtained by multiplying the total eigen-count by the global weights

$$
\#\mathrm{fiber}(t)=m_0\cdot\mathrm{eigen\_count}(t)\cdot w_{\Gamma_i}.
$$

where $m_0$ is global mass (class-number / automorphism factor) and $w_{\Gamma_i}$ gives the multiplicity of level $N$ structures per Frobenius stabale line as

$$
w_{\Gamma_0}=1,\qquad w_{\Gamma_1}=\varphi(N),\qquad w_{\Gamma(N)}=\varphi_{-1}(N).
$$

#### 3.3.2 Generate Structure Report
The fiber strata reveal deeper structure within each eigenform. The full report — including torsion invariants at each local level within the $t$-strata — can be inspected without additional expensive computations.

```python
mod_curve = X1(12).F(5**8)
for fiber in mod_curve.get_structure():
    print(fiber)
```
Each fiber record has the full structure
```text
FiberRecordFq
├── kind: smooth | cusp
├── t: Frobenius trace
├── count: final fiber contribution
└── eigen_records
    └── EigenFormRecord
        ├── eigenform
        │   ├── value: lambda
        │   └── form: BinaryQuadraticForm | CuspForm
        └── levels
            └── LevelStructureRecord
                ├── mass: local multiplicity
                ├── inv: rational N-torsion subgroup
                ├── num_lines: number of stable lines
                └── smooth-only
                    ├── scalar: whether pi acts as a scalar
                    ├── coords: coordinates of (pi - lambda) in an
                    │          OK basis
                    └── order: imaginary quadratic order
```
#### 3.3.2.1 Example : Smooth fiber of $X_0(12)(\mathbb{F}_{13^4})$ | $t=322$, $D_K=-660$, size: 72
```bash
python examples/xnfq_structure.py -p 13 -n 4 -N 12 --type 0
```
outputs

Smooth fiber of $X_0(12)(\mathbb{F}_{13^4})$ | $t=322$, $D_K=-660$, $m_0=4$, size: 72

λ = 5

| f | mass | inv       | stable lines | structure count | level contrib |
|---:|:----:|:----------:|------------:|----------------:|--------------:|
| 1 | 1    | (4, 12)   | 6           | 6               | 24            |
| 2 | 2    | (2, 12)   | 2           | 2               | 16            |
| 4 | 4    | (1, 12)   | 1           | 1               | 16            |

λ = 11

| f | mass | inv       | stable lines | structure count | level contrib |
|---:|:----:|:----------:|------------:|----------------:|--------------:|
| 1 | 1    | (2, 6)    | 0           | 0               | 0             |
| 2 | 2    | (2, 6)    | 0           | 0               | 0             |
| 4 | 4    | (1, 12)   | 1           | 1               | 16            |

#### 3.3.3 Trace
Compute Frobenius trace on the space of cusp forms of weight $k$ and level $N$ via
```python
X1(12).F(13**3).tr_fq(k)
```
For an instantiated $X_i(N)(\mathbb{F}_q)$ we compute trace via the formula

$$
\mathrm{tr\_fq}(k)=
-\sum_{t^2 \leq 4q} h_k(t,q,k-2)\cdot\#\mathrm{fiber}(t)
-\sum_{t^2 = 1} t^{k}\cdot\#\mathrm{fiber}(t),
$$

Here $h_k(t,q,m)$ denotes the complete homogeneous polynomial used above (implemented as `hk(t,q,k)` in `xnfq.moduli.fq.curve`).

#### 3.3.3.1 Example : Basic usage
```bash
python examples/tr_fq.py -p 43 -N 11 -k 2 --type 1
```
#### 3.3.3.2 Example : Prime-range
```bash
python examples/tr_fq.py -p 10000 --range 5 -N 1 -k 12 --type 1 --use-pari
```

Recovering coefficients for the Ramanujan tau function (k=12).

| p     | trace                       | smooth                      | cusp |
|:-----:|:---------------------------:|:---------------------------:|:----:|
| 10007 | -5758585642481476962744    | 5758585642481476962743     | 1    |
| 10009 |  2193555277314013067210    | -2193555277314013067211    | 1    |
| 10037 | -16325223441888896359314   | 16325223441888896359313    | 1    |
| 10039 |  16290359320023605948840   | -16290359320023605948841   | 1    |
| 10061 |  -6292973131198206403338   | 6292973131198206403337     | 1    |

Computed 5 traces in 0.286272s

#### 3.3.3.2 Example : Sage verification (primes only)
```bash
python examples/tr_fq.py -p 13 -N 5 -k 2 --sage
```
If `--sage` is given and Sage is available, `tr_fq.py` will attempt a reference trace per-prime (only for prime fields `n==1`). When `--sage` is used the table includes two extra columns `sage_ref` and `sage_err` showing the Sage reference trace and the absolute error. WARNING: The sage computation is very slow and is not feasable for $p > 100$ or $N > 20$

## 04. Examples

This repository includes a set of example scripts under the `examples/` folder. Below is a concise, organized guide to the common examples, their arguments, and recommended usage.

Common CLI arguments (supported by most examples):

```bash
-N <int>      # level N
-p <int>      # base prime p
-n <int>      # extension degree n (so q = p**n)
-k <int>      # weight k (integer >= 2)
--type <0|1|2> # Gamma type (0,1,2)
--use-pari     # enable optional PARI acceleration
--sage         # run a Sage comparison when supported
--help         # show script-specific options
```

### 04.1 Base
- `examples/xnfq_count.py` — point counts for a given level/field:

- `examples/xnfq_structure.py` — detailed structure report for a single curve/fiber:

- `examples/tr_fq.py` — compute canonical `tr_fq(k)` for a single F_q specialization:

### 04.2 Interactive / plotting
- `examples/xnfq_interactive.py` — Probe the distribution of level structures varying $N$ and extention field $p^n$, normalized fiber count plots by lattice $D_K$ (smooth) and $d$ (cusps)

![Smooth and cusp fiber counts](Figure_1.png)

- `examples/tr_tp_interactive.py` — Interactive prime-range trace plot, vary weight $k$, level $N$ and step through primes, recovering Hecke traces

![Hecke trace](Figure_2.png)


## 05. TODO

- Add precomputed class numbers database.
- Probe individual lattice structures $(N, \pi-\lambda)$ to recover isogeny graphs.
- Embed lattice orders $\mathbb{Z} \oplus \mathbb{Z}\tau$ into the full modular curve $\mathbb{H} \setminus \Gamma$.
