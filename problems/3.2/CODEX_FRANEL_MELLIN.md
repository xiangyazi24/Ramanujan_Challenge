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
- **[VERIFIED-0.2]** The Franel period is the period of the explicit elliptic
  family `E_u: y^2+(1-2u)xy+u^2y=x^3`.  Its Hasse--Witt scalar is `H_p`; on
  smooth fibers its integral Frobenius trace is congruent to `H_p (mod p)`.
  On the ordinary locus the HLYY/Dwork quotient gives the actual unit root;
  `H_p` itself is only its first mod-`p` approximation.
- **[VERIFIED-0.3]** The characteristic-zero rank-two equations and their
  pullback/twist relation are exact.  Their local exponents and the resulting
  tame conductor bookkeeping are explicit.
- **[NEGATIVE-PAIR]** The fixed pair is geometrically redundant:
  `phi_* Sym^2(F)` is the sum of the Apéry object and its `q`-quadratic twist.
  The virtual trace in the specification therefore collapses to minus the
  untwisted Apéry trace.
- **[NEGATIVE-CFVZ-R2]** CFVZ's factorization `A_p=B_p^2` or `qB_p^2` is a
  statement in `F_p[t]`.  It does not construct a rank-two `ell`-adic sheaf,
  and twisting a rank-two sheaf by a quadratic character cannot create the
  companion after symmetric square because `Sym^2(S tensor L)=Sym^2(S)`.
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
- **[VERIFIED-KATZ-H]** Katz, Theorem 27.1, already allows a sequence of
  finite fields of different characteristics.  Its hypotheses require a
  common arithmetic/geometric Tannakian group and uniform generic-rank and
  bad-character bounds; verifying those hypotheses here remains part of
  **[GAP-2]**.

## 1. Sheaf-theoretic normalization

### 1.1 The Franel rank-two object

The toric Franel cubic is

\[
 C_x:\quad XYZ=x(X+Y)(Y+Z)(Z+X).
\]

Its invariant period at `x=0` is
`h(x)=sum_n CT(((X+Y)(Y+Z)(Z+X)/(XYZ))^n)x^n`.  The projective
linear substitution

\[
 [X:Y:Z]\longmapsto[U:V:W]
 =[x^2(X+Y+Z):x^2(xX-Y+xZ):x(X+Y)-Z]
\]

puts it into the explicit Weierstrass form

\[
 E_x:\quad V^2+(1-2x)UV+x^2V=U^3.                            \tag{1.1}
\]

The inverse substitution and both cleared cubic identities are checked in
`codex_fm_geometry.py`.  The invariants of (1.1) are

\[
 \Delta=x^6(1+x)^2(1-8x),\qquad
 c_4=(1-2x)^4-24(1-2x)x^2.
\]

Thus, in every good characteristic `p>3`, the fibers at
`0,-1,1/8,infinity` have Kodaira types `I_6,I_2,I_1,I_3`.  Put

\[
 U_x=\mathbf P^1\setminus\{0,-1,1/8,\infty\},\qquad
 \mathcal F=R^1\pi_*\overline{\mathbf Q}_\ell.
\]

This is a rank-two, weight-one lisse sheaf.  Completing the square in (1.1)
shows that its Hasse invariant is the coefficient of `U^(p-1)` in

\[
 \left(4U^3+((1-2x)U+x^2)^2\right)^{(p-1)/2}.
\]

Direct expansion, equivalently HLYY Theorem 1.2 applied to `C_x`, gives
exactly `H_p(x)`: if `P=(X+Y)(Y+Z)(Z+X)`, then

\[
 [(XYZ)^{p-1}](XYZ-xP)^{p-1}
 \equiv\sum_{n=0}^{p-1}[(XYZ)^n]P^n\,x^n
 =\sum_{n=0}^{p-1}F_nx^n.
\]

Here `binom(p-1,n)(-1)^n=1 mod p`.  Hence, for the integral Frobenius
trace `a_{p,x}`,

\[
 a_{p,x}\equiv H_p(x)\pmod p.                                \tag{1.2}
\]

The script checks (1.2) at every smooth parameter for every prime
`5<=p<=101`, including `x=1/2`, which is hidden by the usual rational
hypergeometric coordinate.

There are two useful meanings of `FranelSquare`:

\[
 \mathcal F^{\square}=\mathcal F\otimes\mathcal F,
 \qquad
 \mathcal F^{(2)}=\operatorname{Sym}^2\mathcal F.
\]

The first has exact trace `a_{p,x}^2`.  The second is the minimal rank-three
object, because

\[
 \mathcal F^{\square}
 =\mathcal F^{(2)}\oplus\det(\mathcal F),\qquad
 \operatorname{Tr}(\operatorname{Frob}_x\mid\mathcal F^{(2)})
 =a_{p,x}^2-p\equiv H_p(x)^2\pmod p.                         \tag{1.3}
\]

HLYY Theorem 1.5 also supplies the literal unit-root statement on the smooth
ordinary locus.  Its `p`-adic continuation

\[
 g_p(x)=h(x)/h(x^p),\qquad g_p(x)\equiv H_p(x)\pmod p,
\]

satisfies, for `x_bar in F_{p^a}` and its Teichmuller lift `x_hat`,

\[
 \alpha_{\rm unit}(\bar x)=
 \prod_{i=0}^{a-1}g_p(\widehat{x}^{p^i}).                    \tag{1.4}
\]

Thus `H_p` is the first approximation, not the unit root itself; there is no
unit root at a supersingular parameter.  The elementary Lucas identity
`h=H_p h^p` is only the reduction of (1.4).

Combining the two Lucas identities with
`f_alpha(phi(x))=(1+x)h(x)^2` gives

\[
 A_p(\phi(x))=(1+x)^{1-p}H_p(x)^2\quad\text{in }\mathbf F_p(x). \tag{1.5}
\]

The exponent `1-p` is forced by direct substitution and is checked after
clearing denominators for nine primes.  The opposite sign printed in the
checked CFVZ display is a typo.  Over `F_p`, Fermat reduces (1.5) to the
inherited pointwise identity; it does not do so over arbitrary extensions.

Since `phi^{-1}(G_m)=U_x`, define

\[
 \mathcal P^\square=\phi_*\mathcal F^\square\quad(\text{rank }8),
 \qquad
 \mathcal P=\phi_*\mathcal F^{(2)}\quad(\text{rank }6).
\]

The trace of finite pushforward is the rational-fiber sum; a
Frobenius-swapped non-rational pair contributes zero.  The inherited identity
and the exact fiber count therefore give, for every `t in F_p^*`,

\[
 \operatorname{Tr}(\operatorname{Frob}_t\mid\mathcal P^\square)
 \equiv
 \operatorname{Tr}(\operatorname{Frob}_t\mid\mathcal P)
 \equiv(1+\chi_2(q(t)))A_p(t)\pmod p.                        \tag{1.6}
\]

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
 \equiv\chi_2(q(t))A_p(t)\pmod p.                             \tag{1.7}
\]

The congruence on the smooth K3 locus follows from its Hasse--Witt scalar and
crystalline--`ell`-adic compatibility; extension by zero makes both sides zero
at the roots of `q`.  A source proving the chosen integral compatible-system
normalization and every middle-extension stalk was not located: that exact
arithmetic upgrade is **[GAP-1]**.

Thus, for the companion sum called `T(r)` in the specification, the precise
bounded-object statement is

\[
 \boxed{
 T(r)\equiv\operatorname{Mell}_p
 \left(t\mapsto\operatorname{Tr}(\operatorname{Frob}_t\mid\mathcal G_T);
       \omega^{-r}\right)\pmod{\mathfrak p},
 \quad \mathcal G_T:=\mathcal Q,
 \quad\operatorname{cond}(\mathcal G_T)=11.}                 \tag{1.8}
\]

Here and below equality with the original `F_p`-valued sum means reduction at
the chosen prime `mathfrak p` above `p`; it is not equality of the complex
Mellin value with an element of `F_p`.

The CFVZ convention is exactly

\[
 A_p(t)=q(t)^{\epsilon_p}B_p(t)^2,
\quad
 \epsilon_p=\frac{1-(\frac{-6}{p})}{2}.                      \tag{1.9}
\]

where `epsilon_p=0` for `p mod 24` in `{1,5,7,11}` and `epsilon_p=1`
for `{13,17,19,23}`.  CFVZ prove (1.9) as a factorization of the reduced
generating series; they do **not** state an `ell`-adic trace-sheaf theorem or a
unit-root theorem.  The two square roots are the truncations of the solutions
`S_+` and `S_-` below, with `S_-=S_+/sqrt(q)`.

In particular, (1.9) does not by itself make `B_p` a Frobenius trace.  Nor can
the companion be obtained by first quadratically twisting a rank-two object,
because

\[
 \operatorname{Sym}^2(\mathcal S\otimes\mathcal L_q)
 \simeq\operatorname{Sym}^2\mathcal S.
\]

The safe companion is the rank-three K3 object `mathcal Q` above.  Peters's
rank-two modular differential module `mathcal S_+`, whose symmetric square is
the Apéry module, explains the characteristic-zero factorization but does not
close **[GAP-1]** arithmetically.  This is **[NEGATIVE-CFVZ-R2]**.

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

The literal tensor-square version is

\[
 [\mathcal G^\square]=[\mathcal Q]-[\mathcal P^\square].
\]

Its trace has the same reduction modulo `p`, because the difference is the
pushforward of the determinant/Tate summand `det(F)`, and every Frobenius
trace of that summand is divisible by `p`.

Choose the prime `mathfrak p` of `Q(mu_{p-1})` for which the Teichmuller
character satisfies `omega(t)=t (mod mathfrak p)`.  Conditional on **[GAP-1]**
for the chosen K3 compatible system, (1.6)--(1.7) give

\[
 \boxed{
 b_r\equiv
 \operatorname{Mell}_p
 \left(t\mapsto\operatorname{Tr}(\operatorname{Frob}_t\mid\mathcal G);
       \omega^{-r}\right)
 \pmod{\mathfrak p}}
 \qquad(1\le r\le p-2).                                     \tag{1.10}
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
 \operatorname{cond}(\mathcal P^\square)=24,\qquad
 \operatorname{cond}(\mathcal Q)=11,
 \qquad \operatorname{cond}(\mathcal G)\le31,
 \qquad \operatorname{cond}(\mathcal G^\square)\le35.       \tag{1.11}
\]

Thus the canonical exact tensor-square choice has the absolute bound `C=35`;
after discarding the Tate determinant, which is zero modulo `p`, the reduced
choice has `C=31`.  The companion `T(r)` alone has `C=11`.  After the
cancellation in Section 2.3, the surviving rank-three object has conductor
`9`.

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

The Kodaira fibers in Section 1.1 and the local Frobenius recurrences checked
by the script show that the repeated exponents really have logarithmic second
solutions; they are nontrivial unipotent blocks, not apparent singularities.
The half-integral differences at the roots of `q` give eigenvalues `{1,-1}`.

### 2.2 The sheaves in (1.10)

For the literal square `mathcal P^square=phi_*(mathcal F tensor mathcal F)`:

| point | local monodromy on rank 8 | drop |
|---|---|---:|
| `0` | two copies of `Unip(3) direct-sum 1` | 4 |
| each root of `q` | four `+1` and four `-1` eigenvalues | 4 |
| `infinity` | two copies of `Unip(3) direct-sum 1` | 4 |

For `mathcal P=phi_* Sym^2 mathcal F`:

| point | local monodromy on rank 6 | drop |
|---|---|---:|
| `0` | `Unip(3) direct-sum Unip(3)` | 4 |
| each root of `q` | three `+1` and three `-1` eigenvalues | 3 |
| `infinity` | `Unip(3) direct-sum Unip(3)` | 4 |

The untwisted Apéry object `mathcal K=Sym^2 mathcal S_+` has `Unip(3)`
at `0` and `infinity`, and eigenvalues `{1,-1,1}` at each root of `q`.
Its drops are therefore `2,1,1,2` and its conductor is `9`.

For `mathcal Q=mathcal K tensor mathcal L_q`, with
`mathcal K = Sym^2 mathcal S_+` at the differential-module level:

| point | local monodromy on rank 3 | drop |
|---|---|---:|
| `0` | `Unip(3)` | 2 |
| each root of `q` | eigenvalues `{-1,1,-1}` | 2 |
| `infinity` | `Unip(3)` | 2 |

All entries are tame for primes of good reduction, yielding (1.11).

#### Global geometric monodromy versus Mellin monodromy

The Franel elliptic pencil is the Beauville `Gamma_1(6)` modular family, so
the rank-two monodromy is a finite-index subgroup of `SL_2(Z)` and is Zariski
dense in `SL_2`.  Its symmetric square has connected Zariski closure `SO_3`.
Independently, Peters, Theorem 7.2.1, identifies the Apéry rank-three variation
with the modular variation on `Z(6)` and its lattice monodromy with
`Gamma_0(6)^*`; its connected Zariski closure is again `SO_3`.  Because of
(2.5), the connected monodromy of the two rank-three constituents is the same
diagonal `SO_3`, not a product of two independent groups.

This ordinary geometric monodromy is not automatically Katz's Tannakian group
for multiplicative convolution.  Computing that convolution group, including
its arithmetic components, is exactly **[GAP-2]**.

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

Consequently the reduced virtual trace in (1.10) is

\[
 -\mathcal P+\mathcal Q=-\mathcal K                           \tag{2.7}
\]

in the Grothendieck group.  This is **[NEGATIVE-PAIR]**: the new decomposition
does not create two independent sheaves whose correlations can be separated.
It reconstructs the original Apéry Hasse--Witt Mellin transform.

For the literal tensor square one instead has

\[
 -\mathcal P^\square+\mathcal Q
 =-\mathcal K-\phi_*\det(\mathcal F).                         \tag{2.8}
\]

The last term is Tate of weight two and has trace divisible by `p`, so (2.8)
has the same residual Mellin value as `-mathcal K`.

The upgrade of (2.5)--(2.6) from the verified differential modules to the
specific arithmetic middle extensions, including their bad stalks, is the
precise content of **[GAP-1]**.  It is not supplied by CFVZ, whose theorem is
about Kummer extensions of reduced generating series.

### 2.4 Kummer self-twists and mutual twists

A geometric Kummer sheaf `mathcal L_rho(t)` on `G_m` is ramified only at `0`
and `infinity`.  Let its local scalar at `0` be `lambda`.

- At `0`, every eigenvalue of `mathcal P^square`, `mathcal P`, and `mathcal Q`
  is `1`.  Twisting changes every one to `lambda`.  Isomorphism therefore
  forces `lambda=1`, hence `rho` is trivial.
- Thus none of the displayed sheaves has a nontrivial Kummer self-twist.  In
  particular, no twist of unbounded order fixes any of them.
- Independently, if a geometrically irreducible rank-`n` sheaf `V` satisfies
  `V tensor L_rho = V`, determinants force `L_rho^n=1`.  For a semisimple
  rank-`n` sheaf, twisting permutes irreducible constituents; following an
  orbit of length `m` and constituent rank `d` gives `L_rho^(md)=1` with
  `md<=n`.  Fixed rank therefore rules out unbounded-order self-twists even
  before using the sharper local calculation.
- `mathcal P^square`, `mathcal P`, and `mathcal Q` have ranks `8,6,3`, so no
  Kummer twist can identify either pushforward with the companion.
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

Theorem 7.2 of *Convolution and Equidistribution* treats extensions `E/k` of a
fixed finite field.  More importantly here, Theorem 27.1 is already a
horizontal theorem: it allows any sequence `k_i` of finite fields, explicitly
including prime fields of different characteristics, with `#k_i -> infinity`.
It fixes `ell`, omits the one characteristic `ell`, and assumes:

1. arithmetically semisimple perverse objects `N_i`, pure of weight zero;
2. a single reductive group `G` with
   `G_geom,N_i=G_arith,N_i=G` in a fixed faithful representation;
3. uniform bounds for generic rank and the number of bad characters.

Under these hypotheses the good-character Frobenius conjugacy classes become
Haar equidistributed in a compact form of `G`.  Thus there is no missing
cross-prime Katz theorem.  What is missing here is the common arithmetic and
geometric Tannakian group (and the compatible integral realization needed in
Section 1).  Local monodromy and the absence of Kummer self-twists constrain
that group but do not determine it.  This is **[GAP-2]**.

Suppose **[GAP-2]** were closed and the trace (or virtual trace) were not
identically zero on any Haar-supported component.  Its zero set is then a
proper real-analytic subset and has Haar measure zero.  Weak equidistribution
and the closed-set inequality would give

\[
 \frac{\#\{\rho:S(\mathcal N,k_i,\rho)=0
                    \text{ as an algebraic number}\}}
      {\#\operatorname{Good}(k_i,\mathcal N)}\longrightarrow0. \tag{3.1}
\]

This is the exact complex/algebraic zero-density conclusion available under
standard Katz hypotheses.  It says nothing about reduction of a nonzero
algebraic integer modulo a moving prime.

### 3.2 The defining-characteristic caveat

The maps `t -> t^r` are genuine `F_p`-valued multiplicative characters, and
as `r` varies they are all characters of `F_p^*`.  The issue is not their
existence.  The issue is the coefficient category: Weil II and Katz use
`ell`-adic sheaves with `ell != p`.

Choose the Teichmuller lift `omega:F_p^* -> mu_{p-1}` and the prime
`mathfrak p | p` for which `omega(t)` reduces to `t`.  The desired event is

\[
 b_r=0\text{ in }\mathbf F_p
 \quad\Longleftrightarrow\quad
 S(\mathcal G,\mathbf F_p,\omega^{-r})\equiv0\pmod{\mathfrak p}. \tag{3.2}
\]

One must first form the algebraic-integer trace using `ell`-adic Kummer sheaves
and Teichmuller lifts, and only then reduce it at `mathfrak p`.  Katz studies
the same algebraic integer through an `ell`-adic-to-complex embedding.  That
archimedean image does not detect its valuation at the moving prime
`mathfrak p`.

The exact Gaussian-integer check in `codex_fm_residual_caveat.py` is the
smallest counterexample to the inference:

\[
 2-i\ne0\text{ in }\mathbf C,\qquad N(2-i)=5,
 \qquad 2-i\equiv0\pmod{(5,i-2)}.                              \tag{3.3}
\]

It is realized by a four-term Teichmuller Mellin sum over `F_5^*`.  More
decisively, for integer values of natural size `p^(3/2)`, changing each by at
most `p` changes the normalized complex value by `O(p^(-1/2))` while allowing
all values to be made `0 mod p`; adding `1` instead makes none of them zero.
The two modified families have the same limiting complex distribution and
residual zero densities `1` and `0`.  The script checks this construction.

Sawin's actual Witt-vector twist theorem (arXiv:1805.04330, Theorems 1.2--1.3)
also uses `ell`-adic Galois representations/middle-extension sheaves and
unitary Frobenius classes.  It varies highly ramified Artin--Schreier--Witt
characters over function fields; it is not a theorem about the reduction of
Mellin values modulo the defining characteristic.  No Sawin result with the
needed residual conclusion was located.

### 3.3 The exact additional hypothesis that would give zero density

Let the residual value be

\[
 m_p(r)=S(\mathcal G,\mathbf F_p,\omega^{-r})\bmod\mathfrak p
 \in\mathbf F_p.
\]

Additive-character orthogonality gives the exact formula

\[
 Z_p:=\#\{1\le r\le p-2:m_p(r)=0\}
 =\frac1p\left((p-2)+
   \sum_{s\in\mathbf F_p^*}\sum_{r=1}^{p-2}e_p(s m_p(r))\right). \tag{3.4}
\]

Consequently the defining-characteristic estimate

\[
 \max_{s\in\mathbf F_p^*}
 \left|\sum_{r=1}^{p-2}e_p(s m_p(r))\right|\le C\sqrt p     \tag{AS}
\]

would imply `Z_p=O(sqrt(p))`, hence the desired zero density.  A uniform local
limit or merely `max_a #\{r:m_p(r)=a\}=o(p)` would also suffice.

A Deligne proof of `(AS)` would require a bounded-conductor geometric object
in the `r`-variable whose trace is `m_p(r)`, and geometric nontriviality of all
its Artin--Schreier pullbacks by `s != 0`.  Katz's Mellin Tannakian object
controls the lifted Frobenius conjugacy classes, not such a residual
`r`-parameter object.  No bounded-complexity construction of this kind is
known here.

No theorem cited here supplies `(AS)` or an equivalent atom bound.  This is
**[GAP-3]**, and
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
VERIFIED Delta(E_u)=u^6(1+u)^2(1-8u) and fiber types I6,I2,I1,I3
VERIFIED toric cubic Hasse coefficient equals H_p by constant terms
VERIFIED projective linear isomorphism from the toric Franel cubic to E_u
VERIFIED a_p(E_u)=H_p(u) mod p at every smooth u over all primes 5<=p<=101 (1084 fibers, including u=1/2)
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
VERIFIED repeated exponents have logarithmic second solutions and nontrivial unipotent blocks
VERIFIED tame conductor bookkeeping: cond(phi_* Sym^2 F)=20, cond(K tensor L_chi(q))=11
VERIFIED exact tensor-square bookkeeping: cond(phi_*(F tensor F))=24 and total C=35
VERIFIED reduced virtual object has C=31 and surviving Apery object has conductor 9
VERIFIED local eigenvalue test at t=0 forces every Kummer self-twist scalar to be 1
VERIFIED ranks 6 and 3 rule out a mutual Kummer twist of the two displayed sheaves

$ python3 problems/3.2/research/scripts/codex_fm_residual_caveat.py
VERIFIED residual Mellin sum is 0 in F_5
VERIFIED its Teichmueller lift is 2-i != 0 in C and has norm 5
VERIFIED 2-i reduces to 0 at the prime (5, i-2)
VERIFIED O(p^-1/2) normalized perturbations can force residual zero density 1 or 0
VERIFIED additive-character orthogonality formula for residual zero counts
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
   especially Chapter 2 for `ell != p`, Theorem 7.2 for the fixed-base-field
   theorem, and Theorem 27.1 for sequences of finite fields of possibly
   different characteristics.  All distributions are read through a complex
   embedding.
3. A. Huang, B. Lian, S.-T. Yau, C. Yu,
   [*Hasse-Witt matrices, unit roots and period integrals*](https://arxiv.org/abs/1801.01189),
   Theorem 1.2 for the toric Hasse--Witt/truncated-period comparison and
   Theorem 1.5 for the ordinary unit-root formula.
4. J. Stienstra, F. Beukers,
   [*On the Picard-Fuchs Equation and the Formal Brauer Group of Certain Elliptic K3-Surfaces*](https://doi.org/10.1007/BF01455990),
   for the Apéry Picard--Fuchs/K3 and formal-Brauer setting.  It is not used
   here as a citation for the unverified all-stalk isomorphism in **[GAP-1]**.
5. O. Gorodetsky,
   [*New representations for all sporadic Apéry-like sequences, with applications to congruences*](https://doi.org/10.1080/10586458.2021.1982080),
   for the exact constant-term representations of the sporadic Apéry-like
   sequences, including the Laurent-polynomial model used for `Lambda_A`.
6. C. Peters,
   [*Monodromy and Picard-Fuchs equations for families of K3-surfaces and elliptic curves*](https://www.numdam.org/item/ASENS_1986_4_19_4_583_0/),
   especially Sections 6.5--6.6 for symmetric squares and the rank-two
   modular equation underlying the Apéry rank-three variation.
7. D. Zagier,
   [*Integral solutions of Apéry-like recurrence equations*](https://people.mpim-bonn.mpg.de/zagier/files/tex/AperylikeRecEqs/fulltext.pdf),
   Section 7 for the Beauville `Gamma_1(6)` cubic family.  The exact linear
   change to (1.1) is verified in the script rather than attributed to Zagier.
8. W. Sawin,
   [*The equidistribution of L-functions of twists by Witt vector Dirichlet characters over function fields*](https://arxiv.org/abs/1805.04330),
   Theorems 1.2--1.3.  These are explicitly `ell`-adic/unitary
   equidistribution results, not defining-characteristic residual estimates.

## Least-confident step

**[GAP-3] is the least-confident prospective step:** no checked Katz or Sawin
theorem controls the moving-prime reductions, and no bounded-conductor object
in the character-index variable `r` is presently available from which `(AS)`
could follow.  **[GAP-1]** is narrower but still real: the exact arithmetic
middle-extension/projector normalization should not be inferred from CFVZ's
polynomial factorization.
