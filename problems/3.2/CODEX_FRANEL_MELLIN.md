# Franel--Mellin bounded object: normalization, monodromy, and the residual obstruction

This ledger separates three levels which must not be conflated:

1. an identity of functions with values in `F_p`;
2. the reduction modulo a prime above `p` of an algebraic-integer Frobenius trace;
3. the complex value obtained from an `ell`-adic trace after an embedding into `C`.

The finite-field decomposition in the specification is exact.  The main new
structural conclusion is negative: its two displayed terms are not independent
geometric objects.  At the Picard--Fuchs level the pushforward splits into the
untwisted Apéry constituent and its quadratic companion, so the companion in
the decomposition cancels one constituent and leaves the original Apéry
Mellin transform.

## 0. Notation and status summary

Put

\[
 q(t)=t^2-34t+1,\qquad
 \phi(x)=\frac{x(1-8x)}{1+x},
\]

\[
 F_n=\sum_k\binom nk^3,\quad
 h(x)=\sum_{n\geq0}F_nx^n,\quad
 H_p(x)=\sum_{0\leq n<p}F_nx^n,
\]

and let `A_p(t)` be the degree-`p-1` Apéry truncation.  For a multiplicative
character `rho` write

\[
 \operatorname{Mell}_p(g;\rho)=\sum_{t\in\mathbf F_p^\times}g(t)\rho(t).
\]

The ledger conclusions are:

- **[VERIFIED-0.1]** The inherited pointwise, fiber, and Mellin identities
  reproduce in the new scripts, including all `1 <= r <= p-2` for nine test
  primes.
- **[VERIFIED-0.2]** The Franel period is the period of an explicit toric
  elliptic family.  Its Hasse--Witt scalar is `H_p`; on smooth fibers its
  integral Frobenius trace is congruent to `H_p (mod p)`.
- **[VERIFIED-0.3]** The characteristic-zero rank-two equations and their
  pullback/twist relation are exact.  Their local exponents and the resulting
  tame conductor bookkeeping are explicit.
- **[NEGATIVE-PAIR]** The fixed pair is geometrically redundant:
  `phi_* Sym^2(F)` is the sum of the Apéry object and its `q`-quadratic twist.
  The virtual trace in the specification therefore collapses to minus the
  untwisted Apéry trace.
- **[NEGATIVE-KATZ-p]** Katz's Mellin equidistribution concerns pure
  `ell`-adic traces and complex embeddings.  It does not turn reduction
  modulo the defining characteristic into complex cancellation and gives no
  zero-density theorem for `b_r = 0 in F_p`.
- **[GAP-1]** A source-level arithmetic comparison at every middle-extension
  stalk (especially `q(t)=0`) is still needed before the displayed
  characteristic-zero local-system splitting can be quoted literally as an
  isomorphism of the chosen integral compatible systems with the asserted
  mod-`p` trace at every `t`.
- **[GAP-2]** The full arithmetic Tannakian monodromy of the Mellin object has
  not been computed.  Local monodromy and absence of Kummer self-twists do not
  by themselves determine Katz's group.
- **[GAP-3]** Residual anti-concentration for the reductions of the Mellin
  traces is the actual missing zero-density input.  It is not a theorem in the
  cited Katz machinery.
- **[GAP-4]** The first Hasse--Witt congruence has not been upgraded here to a
  literal unit-root formula with a chosen Frobenius lift and higher Dwork
  limit.
- **[GAP-5]** Katz's extension-field equidistribution has not been replaced by
  a uniform horizontal theorem for the varying-prime family in this problem.

## 1. Sheaf-theoretic normalization

### 1.1 The Franel rank-two object

Use the fixed Laurent polynomial

\[
 \Lambda_F(u,v)=(1+u)(1+v)(1+(uv)^{-1}).
\]

Since `F_n = CT Lambda_F^n`, the holomorphic toric period of

\[
 E_x:\quad 1-x\Lambda_F(u,v)=0
\]

is `h(x)`.  Compactifying in the toric surface of the reflexive hexagon
gives a genus-one family over

\[
 U_x=\mathbf P^1\setminus\{0,-1,1/8,\infty\}.
\]

Let

\[
 \mathcal F=R^1\pi_*\overline{\mathbf Q}_\ell
\]

on `U_x` (with the usual middle extension understood).  It has rank two and
weight one.  If `a_{p,x}` denotes its integral Frobenius trace, then

\[
 \operatorname{Tr}(\operatorname{Frob}_x\mid\operatorname{Sym}^2\mathcal F)
 =a_{p,x}^2-p
 \equiv H_p(x)^2\pmod p.                                      \tag{1.1}
\]

The last congruence uses the Hasse--Witt congruence
`a_{p,x} = H_p(x) (mod p)`.  The toric truncation statement is an instance of
Huang--Lian--Yau--Yu, Theorem 1.2; the script also checks the corresponding
point-count congruence at 74 smooth fibers.  The same family gives the formal
Lucas--Dwork congruence

\[
 h(x)=H_p(x)h(x)^p\quad\text{in }\mathbf F_p[[x]].              \tag{1.2}
\]

**[GAP-4]** Equation (1.2) is the first Hasse--Witt congruence.  A literal
unit-root formula requires a chosen Frobenius lift and the higher Dwork limit;
it is not correct to call `H_p(x)` itself the `p`-adic unit root.  On the
ordinary locus `H_p(x) != 0`, it is the first approximation to that unit root.

Combining the two Lucas identities with
`f_alpha(phi(x))=(1+x)h(x)^2` gives the stronger rational-function identity

\[
 A_p(\phi(x))=(1+x)^{1-p}H_p(x)^2\quad\text{in }\mathbf F_p(x). \tag{1.2a}
\]

The exponent `1-p` is forced by direct substitution and is checked after
clearing denominators for nine primes.  (The display `(1+x)^{p-1}H^2` in the
checked CFVZ source has the exponent sign reversed.)  For
`x in F_p \ {-1}`, Fermat reduces (1.2a) to the inherited pointwise identity.
It does not do so over arbitrary finite extensions, so that base-field
identity is not being used here as evidence for an isomorphism of sheaves at
all closed points; the geometric comparison remains **[GAP-1]**.

Define the rank-six pushforward

\[
 \mathcal P=\phi_*\operatorname{Sym}^2\mathcal F.
\]

For an unramified rational `t`, the trace of a finite pushforward is the sum
over rational points of the fiber (a Frobenius-swapped non-rational pair has
trace zero).  Therefore (1.1) gives

\[
 \operatorname{Tr}(\operatorname{Frob}_t\mid\mathcal P)
 \equiv\sum_{\substack{x\in\mathbf F_p\\\phi(x)=t}}H_p(x)^2
 =(1+\chi_2(q(t)))A_p(t)\pmod p.                              \tag{1.3}
\]

The last equality is the inherited pointwise identity plus the exact fiber
count.

### 1.2 The companion

A fixed toric model for the Apéry period is obtained from

\[
 \Lambda_A(u,v,w)
 =\frac{(u+v)(w+1)(u+v+w)(v+w+1)}{uvw},\qquad
 A_n=\operatorname{CT}\Lambda_A^n.
\]

This constant-term realization is among the sporadic Apéry-like models proved
by Gorodetsky; the script checks the displayed Laurent polynomial directly for
`0 <= n <= 6` as a guard against convention changes.

After compactification/resolution, `1-t Lambda_A=0` is the usual Apéry K3
variation.  Denote its rank-three transcendental local system by `mathcal K`.
The Hasse--Witt truncation is `A_p(t)`.  The rank-one quadratic sheaf

\[
 \mathcal L_q=\mathcal L_{\chi_2(q(t))}
\]

has trace `chi_2(q(t))`, extended by zero at the two roots of `q`.  Thus the
companion is

\[
 \mathcal Q=\mathcal K\otimes\mathcal L_q,
 \qquad
 \operatorname{Tr}(\operatorname{Frob}_t\mid\mathcal Q)
 \equiv\chi_2(q(t))A_p(t)\pmod p.                             \tag{1.4}
\]

The CFVZ convention is exactly

\[
 A_p(t)=q(t)^{\epsilon_p}B_p(t)^2,
\quad
 \epsilon_p=\frac{1-(\frac{-6}{p})}{2},                      \tag{1.5}
\]

where `epsilon_p=0` for `p mod 24` in `{1,5,7,11}` and `epsilon_p=1`
for `{13,17,19,23}`.  CFVZ prove (1.5) as a factorization of the reduced
generating series; they do **not** state an `ell`-adic trace-sheaf theorem or a
unit-root theorem.  The two square roots are the truncations of the solutions
`S_+` and `S_-` below, with `S_-=S_+/sqrt(q)`.

### 1.3 The exact Mellin statement, with its comparison gap exposed

The negative sign makes the natural object virtual.  In the Grothendieck group
of constructible complexes set

\[
 [\mathcal G]=[\mathcal Q]-[\mathcal P],
\]

or represent its trace function by the bounded constructible complex
`mathcal P[1] direct-sum mathcal Q`.  This distinction matters: it is not a
single pure lisse sheaf, so Katz must be applied to the two pure constituents
and their joint Tannakian group, with the difference taken afterward.

Choose the prime `mathfrak p` of `Q(mu_{p-1})` for which the Teichmuller
character satisfies `omega(t)=t (mod mathfrak p)`.  Conditional on **[GAP-1]**
at the middle-extension stalks, (1.3)--(1.4) give the precise statement

\[
 \boxed{
 b_r\equiv
 \operatorname{Mell}_p
 \left(t\mapsto\operatorname{Tr}(\operatorname{Frob}_t\mid\mathcal G);
       \omega^{-r}\right)
 \pmod{\mathfrak p}}
 \qquad(1\le r\le p-2).                                      \tag{1.6}
\]

For the conductor convention

\[
 \operatorname{cond}(\mathcal V)
 =\operatorname{rank}(\mathcal V)
  +\sum_s(\operatorname{drop}_s(\mathcal V)+\operatorname{Swan}_s(\mathcal V)),
\]

the local table below gives

\[
 \operatorname{cond}(\mathcal P)=20,\qquad
 \operatorname{cond}(\mathcal Q)=11,
 \qquad \operatorname{cond}(\mathcal G)\le31.                \tag{1.7}
\]

Hence the requested constant can be taken to be `C=31` for the displayed
direct-sum model.  After the cancellation in Section 2.3, the surviving
rank-three object has conductor `9`.

## 2. Local monodromy, twists, and the collapse of the pair

### 2.1 Rank-two equations and exponents

The Franel period satisfies

\[
 x(x+1)(8x-1)h''+(24x^2+14x-1)h'+(8x+2)h=0.                 \tag{2.1}
\]

The two normalized square-root periods satisfy

\[
 4tqS_+''+4(2t^2-51t+1)S_+'+(t-10)S_+=0,                   \tag{2.2}
\]

\[
 4tqS_-''+4(4t^2-85t+1)S_-'+3(3t-26)S_-=0.                \tag{2.3}
\]

The symbolic checks prove

\[
 S_+(\phi(x))=\sqrt{1+x}\,h(x),\qquad S_-=S_+/\sqrt q,     \tag{2.4}
\]

at the level of differential modules.  Their Riemann exponents are:

| rank-two object | first finite point | second finite point | third finite point | `infinity` |
|---|---:|---:|---:|---:|
| `mathcal F` on the `x`-line | `x=0: (0,0)` | `x=-1: (0,0)` | `x=1/8: (0,0)` | `(1,1)` |
| `mathcal S_+` on the `t`-line | `t=0: (0,0)` | first `q`-root: `(0,1/2)` | second `q`-root: `(0,1/2)` | `(1/2,1/2)` |
| `mathcal S_-` on the `t`-line | `t=0: (0,0)` | first `q`-root: `(0,-1/2)` | second `q`-root: `(0,-1/2)` | `(3/2,3/2)` |

The repeated exponents at `0` and `infinity` give nontrivial unipotent
blocks; the half-integral differences at the roots of `q` give eigenvalues
`{1,-1}`.

### 2.2 The two sheaves in (1.6)

For `mathcal P=phi_* Sym^2 mathcal F`:

| point | local monodromy on rank 6 | drop |
|---|---|---:|
| `0` | `Unip(3) direct-sum Unip(3)` | 4 |
| each root of `q` | three `+1` and three `-1` eigenvalues | 3 |
| `infinity` | `Unip(3) direct-sum Unip(3)` | 4 |

For `mathcal Q=mathcal K tensor mathcal L_q`, with
`mathcal K = Sym^2 mathcal S_+` at the differential-module level:

| point | local monodromy on rank 3 | drop |
|---|---|---:|
| `0` | `Unip(3)` | 2 |
| each root of `q` | eigenvalues `{-1,1,-1}` | 2 |
| `infinity` | `Unip(3)` | 2 |

All entries are tame for primes of good reduction, yielding (1.7).

### 2.3 Descent and the negative result

Equation (2.4) implies at the differential-module level

\[
 \phi^*\mathcal K\simeq\operatorname{Sym}^2\mathcal F.       \tag{2.5}
\]

The degree-two cover has discriminant `q`, so

\[
 \phi_*\mathbf 1\simeq\mathbf 1\oplus\mathcal L_q.
\]

Projection formula then gives

\[
 \boxed{
 \mathcal P\simeq\mathcal K\oplus(\mathcal K\otimes\mathcal L_q)
 =\mathcal K\oplus\mathcal Q.}                              \tag{2.6}
\]

This identity already appears at the finite-field trace level as

\[
 (1+\chi_2(q(t)))A_p(t)=A_p(t)+\chi_2(q(t))A_p(t).
\]

Consequently the virtual trace in (1.6) is

\[
 -\mathcal P+\mathcal Q=-\mathcal K                           \tag{2.7}
\]

in the Grothendieck group.  This is **[NEGATIVE-PAIR]**: the new decomposition
does not create two independent sheaves whose correlations can be separated.
It reconstructs the original Apéry Hasse--Witt Mellin transform.

The upgrade of (2.5)--(2.6) from the verified differential modules to the
specific arithmetic middle extensions, including their bad stalks, is the
precise content of **[GAP-1]**.  It is not supplied by CFVZ, whose theorem is
about Kummer extensions of reduced generating series.

### 2.4 Kummer self-twists and mutual twists

A Kummer sheaf `mathcal L_rho(t)` on `G_m` is ramified only at `0` and
`infinity`.  Let its local scalar at `0` be `lambda`.

- At `0`, every eigenvalue of `mathcal P` and `mathcal Q` is `1`.  Twisting
  changes every one to `lambda`.  Isomorphism therefore forces `lambda=1`,
  hence `rho` is trivial.
- Thus neither displayed sheaf has a nontrivial Kummer self-twist.  In
  particular, no twist of unbounded order fixes either one.
- `mathcal P` and `mathcal Q` have ranks 6 and 3, so no Kummer twist can make
  them isomorphic.
- The two rank-three constituents `mathcal K` and `mathcal Q` differ by the
  quadratic sheaf `mathcal L_q`, but this is **not** a Kummer sheaf on `G_m`:
  it is ramified at the two finite roots of `q`.  Moreover their eigenvalue
  multiplicities at a root of `q` are `{1,1,-1}` and `{-1,-1,1}`, so no
  Kummer twist (unramified there) identifies them.

This proves the requested bounded-order statement in the strongest possible
form at the verified local-system level: the Kummer self-twist group is
trivial.

## 3. What Katz-style Mellin equidistribution does and does not give

### 3.1 The theorem that is actually available

Katz starts with an `ell`-adic middle-extension/perverse sheaf on `G_m/k`,
pure of weight zero after normalization.  His Mellin sums are algebraic
numbers

\[
 S(\mathcal N,E,\rho)
 =\sum_{t\in E^\times}\rho(t)
   \operatorname{Tr}(\operatorname{Frob}_{E,t}\mid\mathcal N).
\]

Theorem 1.1 / Theorem 7.2 of *Convolution and Equidistribution* gives
equidistribution of the Frobenius conjugacy classes as `E/k` runs through
larger finite extensions, under purity, arithmetic semisimplicity, property
`P`, and the stated relation between arithmetic and geometric Tannakian
groups.  It uses an embedding of the `ell`-adic coefficient field into `C`.

Local monodromy and the self-twist calculation are useful inputs, but they do
not prove the required Tannakian group.  That remaining computation is
**[GAP-2]**.

Even granting **[GAP-2]**, if the Haar pushforward under the trace map has no
atom at zero, Katz implies zero density for the event

\[
 S(\mathcal N,E,\rho)=0\quad\text{as an algebraic/complex number}. \tag{3.1}
\]

It does not imply zero density for reduction of a nonzero algebraic integer
modulo a prime.

There is also a quantifier mismatch: the main theorem is an extension-field
limit over a fixed base field.  The set requested here is the one field
`E=F_p`, all `p-1` characters, followed implicitly by `p -> infinity`.
A compatible cross-prime family plus a uniform horizontal theorem would be
needed even for the complex version of that limit.  No such theorem is
supplied by the checked sources; this is **[GAP-5]**.

### 3.2 The defining-characteristic caveat

The desired event is

\[
 b_r=0\text{ in }\mathbf F_p
 \quad\Longleftrightarrow\quad
 S(\mathcal G,\mathbf F_p,\omega^{-r})\equiv0\pmod{\mathfrak p}. \tag{3.2}
\]

Katz does not run the Weil-II purity/equidistribution argument with
`F_p`-coefficient sheaves in characteristic `p`.  One must first form the
`ell`-adic/algebraic-integer trace using Teichmuller lifts, and only then reduce
it at `mathfrak p | p`.  The reduction map is not detected by its chosen
complex embedding.

The exact Gaussian-integer check in `codex_fm_residual_caveat.py` is the
smallest counterexample to the inference:

\[
 2-i\ne0\text{ in }\mathbf C,\qquad N(2-i)=5,
 \qquad 2-i\equiv0\pmod{(5,i-2)}.                              \tag{3.3}
\]

It is realized by a four-term Teichmuller Mellin sum over `F_5^*`.  Thus
complex non-cancellation is compatible with residual vanishing.

### 3.3 The exact additional hypothesis that would give zero density

Let

\[
 M_{p,r}=S(\mathcal G,\mathbf F_p,\omega^{-r}).
\]

Any one of the following genuinely residual hypotheses would suffice:

\[
 \#\{1\le r\le p-2:M_{p,r}\equiv0\pmod{\mathfrak p}\}=o(p),   \tag{RAC}
\]

which is the desired conclusion itself in trace language, or the stronger
uniform local limit

\[
 \#\{r:M_{p,r}\equiv a\pmod{\mathfrak p}\}
 =\frac{p-2}{p}+O(p^{1/2})\quad(a\in\mathbf F_p),             \tag{RLL}
\]

which would give `#zeros = O(sqrt(p))`, hence density zero.  A weaker maximum
atom bound `max_a #(...) = o(p)` also suffices.

No theorem cited here supplies `(RAC)` or `(RLL)`.  This is **[GAP-3]**, and
the cancellation (2.7) shows that it is exactly residual anti-concentration
for the original Apéry/K3 Mellin coordinate, not for a newly independent
Franel pair.

## 4. Verification ledger

All commands were run from the repository root.

```text
$ python3 problems/3.2/research/scripts/codex_fm_geometry.py
VERIFIED characteristic-zero pullback through O(x^24)
VERIFIED CT Lambda_A^n=A_n for 0<=n<=6
VERIFIED cover discriminant q(t)=t^2-34t+1 and its square pullback
VERIFIED A_p(phi)=(1+x)^(1-p)H_p^2 in F_p(x) for p=5,7,11,13,17,19,23,29,37
VERIFIED h(x)=H_p(x)h(x)^p mod p through degree 4p-1 for p=5,7,11,13,17,19,29,37
VERIFIED Franel toric Hasse/point-count congruence at 74 smooth fibers
VERIFIED pushforward=(1+chi(q))A_p, virtual cancellation=-A_p, and Mellin=b_r for 143 (p,r) pairs

$ python3 problems/3.2/research/scripts/codex_fm_factorization.py
VERIFIED CFVZ A_p=q^epsilon B_p^2 and the corresponding S_+/S_- rank-two recurrence for 44 primes 5<=p<=199
VERIFIED epsilon=0 on p mod 24 in {1,5,7,11} and epsilon=1 on {13,17,19,23}

$ python3 problems/3.2/research/scripts/codex_fm_local_monodromy.py
VERIFIED phi-pullback of S_+ is sqrt(1+x) times the Franel rank-two equation
VERIFIED S_-=S_+/sqrt(q) transforms the plus equation into the minus equation
VERIFIED Sym^2(S_+) is annihilated by the third-order Apery operator
VERIFIED Franel exponents: 0:(0,0), -1:(0,0), 1/8:(0,0), infinity:(1,1)
VERIFIED S_+ exponents: 0:(0,0), q-roots:(0,1/2), infinity:(1/2,1/2)
VERIFIED S_- exponents: 0:(0,0), q-roots:(0,-1/2), infinity:(3/2,3/2)
VERIFIED tame conductor bookkeeping: cond(phi_* Sym^2 F)=20, cond(K tensor L_chi(q))=11
VERIFIED local eigenvalue test at t=0 forces every Kummer self-twist scalar to be 1
VERIFIED ranks 6 and 3 rule out a mutual Kummer twist of the two displayed sheaves

$ python3 problems/3.2/research/scripts/codex_fm_residual_caveat.py
VERIFIED residual Mellin sum is 0 in F_5
VERIFIED its Teichmueller lift is 2-i != 0 in C and has norm 5
VERIFIED 2-i reduces to 0 at the prime (5, i-2)
```

The scripts are finite verification, not substitutes for **[GAP-1]** or a
residual equidistribution theorem.

## 5. Sources actually checked

1. X. Caruso, F. Fuernsinn, D. Vargas-Montoya, W. Zudilin,
   [*Galois Groups of Apery-like Series Modulo Primes*](https://arxiv.org/abs/2510.23298),
   especially Theorem 2 and the proof using
   `f_alpha(t(x))=(1+x)h(x)^2`.  This source proves the polynomial
   factorization, not a sheaf or unit-root statement.
2. N. Katz,
   [*Convolution and Equidistribution: Sato-Tate Theorems for Finite-Field Mellin Transforms*](https://web.math.princeton.edu/~nmk/mellin331.pdf),
   especially the setup in Chapter 1 and Theorems 1.1 and 7.2.  The coefficient
   objects are `ell`-adic and the distribution is read through a complex
   embedding as extension fields grow.
3. A. Huang, B. Lian, S.-T. Yau, C. Yu,
   [*Hasse-Witt matrices, unit roots and period integrals*](https://arxiv.org/abs/1801.01189),
   Theorem 1.2, for the toric Hasse--Witt/truncated-period comparison.
4. J. Stienstra, F. Beukers,
   [*On the Picard-Fuchs Equation and the Formal Brauer Group of Certain Elliptic K3-Surfaces*](https://doi.org/10.1007/BF01455990),
   for the Apéry Picard--Fuchs/K3 and formal-Brauer setting.  It is not used
   here as a citation for the unverified all-stalk isomorphism in **[GAP-1]**.
5. O. Gorodetsky,
   [*New representations for all sporadic Apéry-like sequences, with applications to congruences*](https://doi.org/10.1080/10586458.2021.1982080),
   for the exact constant-term representations of the sporadic Apéry-like
   sequences, including the Laurent-polynomial model used for `Lambda_A`.

## Least-confident step

**[GAP-1] is the least-confident step:** the differential-module descent and
all smooth-fiber congruences are clear, but the exact arithmetic
middle-extension/projector normalization at the two `q(t)=0` stalks has not
been located in the checked sources and should not be silently inferred from
CFVZ's polynomial factorization.
