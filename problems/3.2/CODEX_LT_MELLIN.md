# The defining-characteristic Mellin frontier: orbit norms and the T-adic audit

## Executive verdict

No power saving for

\[
 Z_p=\{1\le r\le p-2:\mathfrak p\mid \widetilde M(r)-\widetilde T(r)\}
\]

follows from either proposed route.

- **[VERIFIED-GALOIS]** The character convention is covariant:
  `sigma_a(chi_r(t))=chi_{ar}(t)`.  The elliptic traces and the pointwise
  traces defining `T` are rational integers, hence are fixed.  Exact
  cyclotomic computations at `p=13,17` verify the action on every `r` and
  every automorphism.
- **[VERIFIED-NORM]** For one Galois orbit, its product is an integer.  If
  `k` members vanish modulo the same chosen `mathfrak p` and the product is
  nonzero, then `p^k` divides that integer.
- **[NEGATIVE-NORM]** The sharp available Archimedean scale is
  `p^(3/2)` per conjugate.  It yields `k <= (3/2+o(1))|O|`, which is weaker
  than `k<=|O|`.  Orbitwise AM--GM becomes nontrivial only when the orbit
  root-mean-square is `<p`, a factor `sqrt(p)` below the natural scale.
- **[CORRECTION-PARSEVAL]** The Parseval formula proposed in the task omits
  cross terms inside the degree-two fibers of `phi`.  The exact right side
  uses the fiber sums `B_t=sum_{phi(x)=t}a_{p,x}^2`.  This correction makes
  the second moment larger, not smaller.
- **[GAP-CHAR0-ZERO]** If an orbit product is zero, the whole orbit consists
  of characteristic-zero trace zeros.  Both Mellin objects have weight
  three; purity does not exclude this case or bound its orbit size.
- **[NEGATIVE-TADIC]** The `T` in the checked T-adic exponential-sum theory
  interpolates wild additive characters of `p`-power conductor.  A tame
  multiplicative character is fixed in the “twisted T-adic” papers.  These
  theorems do not interpolate the `p-1` twists `omega^{-r}` and give no
  `o(p)` bad-twist bound.
- **[VERIFIED-ISOGENY]** The missing Beauville-IV correspondence can be made
  explicit by quotienting the displayed Weierstrass model by its rational
  cyclic subgroup of order six.  The quotient is a constant `-3` quadratic
  twist of the deck-conjugate fiber; the twist disappears after `Sym^2`.
- **[FALSE-3F2]** The literal
  `3F2(1/3,2/3,1;1,1)` cancels to the rank-two
  `2F1(1/3,2/3;1)`.  CFVZ use this rank-two Franel object; the Apéry object
  is its rank-three symmetric square after a rational pullback.

## 1. The Galois action and the integer orbit product

Put `N=p-1`, choose a primitive root `g` modulo `p`, and realize the lifted
Teichmüller character by

\[
 \omega(g)=\zeta_N,\qquad
 \chi_r(g^j)=\zeta_N^{-rj}.
\]

For `a in (Z/NZ)^times`, let `sigma_a(zeta_N)=zeta_N^a`.  Then

\[
 \sigma_a(\chi_r(t))=\chi_{ar}(t).                         \tag{1.1}
\]

The automorphism acts on the character value, not on `t`; the notation
`chi_r circ sigma_a` would be type-incorrect.  Since every `a_{p,x}` is a
rational integer,

\[
 \sigma_a(\widetilde M(r))=\widetilde M(ar).               \tag{1.2}
\]

The rank-three Apéry companion is defined over `Z` and its quadratic twist
has values in `{0,+1,-1}`, so its pointwise Frobenius traces are also rational
integers.  Therefore the same argument gives

\[
 \sigma_a(\widetilde T(r))=\widetilde T(ar),\qquad
 \sigma_a(\Delta_r)=\Delta_{ar},\quad
 \Delta_r:=\widetilde M(r)-\widetilde T(r).                \tag{1.3}
\]

Let `d=gcd(r,N)` and `q=N/d`.  The orbit is

\[
 O_r=\{s\bmod N:\gcd(s,N)=d\},\qquad |O_r|=\varphi(q).    \tag{1.4}
\]

It is the set of characters of exact order `q`; consequently there is one
orbit for every divisor `q>1` of `N`.  The product

\[
 P_r=\prod_{s\in O_r}\Delta_s                             \tag{1.5}
\]

is fixed by `Gal(Q(zeta_N)/Q)`.  Each factor is an algebraic integer, so

\[
 \boxed{P_r\in\mathbf Z}.                                  \tag{1.6}
\]

This remains true if `Delta_r` happens to lie in a proper subfield: (1.5)
then repeats the smaller-field norm by a fixed multiplicity.

### 1.1 Exact transfer from `mathfrak p` to `p`

Fix the split prime `mathfrak p` selected by `zeta_N -> g mod p`.  Suppose
exactly `k` factors in (1.5) lie in `mathfrak p`.  Then

\[
 P_r\in\mathfrak p^k.                                     \tag{1.7}
\]

Because `p=1 mod N`, it splits completely in `Q(zeta_N)`: the ramification
index and residue degree of `mathfrak p` are both one.  For the rational
integer `P_r`, therefore,

\[
 v_p(P_r)=v_{\mathfrak p}(P_r)\ge k.
\]

Thus, provided `P_r != 0`,

\[
 \boxed{p^k\mid P_r}.                                     \tag{1.8}
\]

Equivalently, after restricting to `Q(zeta_q)`, the `k` hits of conjugates at
the fixed prime correspond to `k` distinct conjugate primes dividing
`Delta_r`, and its norm contains `p^k`.

## 2. The Archimedean bound is intrinsically vacuous

For the generic nonexceptional characters, the inherited cohomology ledger
has Mellin dimensions six for `M` and four for `T`.  Purity gives

\[
 |\widetilde M(r)|\le6p^{3/2},\qquad
 |\widetilde T(r)|\le4p^{3/2},
\]

hence one may take `c=10` in the direct presentation.  The conductor ledger
gives the coarser uniform constant `31` when all middle-extension corrections
are retained.  After the characteristic-zero descent cancellation, the
generic surviving `A_+` Mellin cohomology has dimension two, giving the
strongest constant

\[
 |\Delta_r|\le2p^{3/2}.                                   \tag{2.1}
\]

Any of these constants has the same fatal exponent.  Writing `m=|O_r|`,
(1.8) and the bound `|Delta_s|<=c p^(3/2)` give

\[
 k\le \log_p|P_r|
 \le m\left(\frac32+\log_p c\right).                      \tag{2.2}
\]

Even with `c=2`, the right side is larger than `m`, so (2.2) is weaker than
the trivial `k<=m`.

### 2.1 Orbitwise AM--GM: the exact threshold

Put

\[
 S_O=\sum_{s\in O}|\Delta_s|^2.
\]

AM--GM gives the sharp moment version

\[
 |P_O|\le\left(\frac{S_O}{m}\right)^{m/2},\qquad
 k\le\frac m2\log_p\left(\frac{S_O}{m}\right).            \tag{2.3}
\]

This improves on `k<=m` if and only if

\[
 \frac{S_O}{m}<p^2,                                       \tag{2.4}
\]

i.e. the orbit root-mean-square is `<p`.  More generally,

\[
 S_O/m\le p^{2-2\eta}\quad\Longrightarrow\quad
 k\le(1-\eta)m.                                           \tag{2.5}
\]

The natural and observed root-mean-square is instead of order `p^(3/2)`.
No orbit size changes the exponent in (2.3): small orbits have too little
averaging, while an orbit of size comparable to `p` still has mean square of
order `p^3`.  The available global second moment does not force (2.4) for a
single orbit.  Hence **no orbit structure `q|p-1` is presently
nonvacuous**.

## 3. Parseval: the required fiber correction

The sum defining `M` is not a Fourier transform of independently indexed
`x`-values, because `phi` has degree two.  Define

\[
 B_t=\sum_{\substack{x:\phi(x)=t}}a_{p,x}^2,
 \qquad t\in\mathbf F_p^\times.
\]

Then, including all `N` multiplicative characters,

\[
 \boxed{
 \sum_{r=0}^{N-1}|\widetilde M(r)|^2
   =N\sum_{t\in\mathbf F_p^\times}|B_t|^2.}               \tag{3.1}
\]

In general this is not

\[
 N\sum_x a_{p,x}^4.                                       \tag{3.2}
\]

The difference consists exactly of the cross terms
`2 a_{p,x}^2 a_{p,y}^2` with `x!=y` and `phi(x)=phi(y)`.  The explicit
degree-six correspondence below shows that paired fibers are isogenous up to
a quadratic twist, hence their squared traces agree.  Away from branch
fibers, grouping therefore doubles the contribution relative to (3.2).
Thus the heuristic `2p^4` for the ungrouped expression becomes approximately
`4p^4` for the actual Parseval right side.  The exact identity (3.1), not the
asymptotic constant, is what is used here.

For the surviving rank-three amplitude `d_t` one likewise has

\[
 \sum_{r=0}^{N-1}|\Delta_r|^2=N\sum_t|d_t|^2.              \tag{3.3}
\]

The pointwise Weil bound `|d_t|<=3p` gives only a global mean-square of order
`p^3` per character.  Substituting that scale into (2.3) again gives the
vacuous exponent `3m/2`.  Parseval confirms the obstruction rather than
removing it.

## 4. The zero orbit product cannot be removed by purity

If `P_r=0`, then one `Delta_s` is zero as an algebraic number.  Galois
covariance implies

\[
 \Delta_s=0\quad\text{for every }s\in O_r.                 \tag{4.1}
\]

There is no weight mismatch to exclude (4.1).  Both `M` and `T` come from
weight-two input trace functions followed by `H_c^1`, so their generic
Mellin eigenvalues and their difference have weight three.  After descent,
`Delta_s` is the trace on a two-dimensional weight-three Mellin cohomology
space; a trace of two nonzero Weil numbers can be exactly zero.

Katz's complex equidistribution may predict that exact trace zero has measure
zero in an appropriate extension-field limit, but it supplies no uniform
bound for these characteristic-zero zero orbits as `p` and the tame character
group vary diagonally.  This is **[GAP-CHAR0-ZERO]**.  It is logically prior
to applying a nonzero integer norm bound and could contain an orbit of size
`phi(q)` comparable to `p`.

## 5. Vertical family / T-adic audit

### 5.1 What the variable `T` actually interpolates

Liu--Wan define a locally constant additive character
`psi:Z_p -> C_p^times` whose image has order `p^m`.  Their specialization
`T=pi_psi=psi(1)-1` interpolates wild additive characters of increasing
`p`-power conductor.  Their T-adic entire-function, Hodge-bound, transfer,
and ordinarity theorems all vary this additive character.

In Liu--Niu's “twisted T-adic” theory, the multiplicative character is
written `chi=omega^{-u}` and is fixed.  The variable `T` still changes only
the additive `p`-power character.  The Hasse polynomial `H_u` depends on the
fixed exponent `u`; the theorems do not assemble the `q-1` possible values of
`u` into one family or count the exceptional `u`.

Davis--Wan--Xiao similarly study characters of the Galois group `Z_p` of an
Artin--Schreier--Witt tower.  Their slope stability as the conductor grows is
a wild pro-`p` statement.  Our characters

\[
 \chi_r:\mathbf F_p^\times\longrightarrow\mu_{p-1}
\]

are tame and prime-to-`p`.  They are not specializations of that weight
space.  Consequently, the checked T-adic/eigencurve results give **no bound
better than `O(p)`**, and in particular neither `o(p)` nor `O(polylog p)`.

### 5.2 A direct obstruction to a `Z_p` interpolation of `r`

For a fixed nontrivial tame root `zeta`, suppose a continuous function on
`Z_p` agreed with `n -> zeta^n` at every integer.  Since `p^j -> 0` in
`Z_p`, continuity would require `zeta^(p^j) -> 1`.  But

\[
 p^j\equiv1\pmod{p-1},\qquad \zeta^{p^j}=\zeta\ne1.        \tag{5.1}
\]

Thus the tame character indices do not lie on one natural rigid disk.  They
form a zero-dimensional finite etale character scheme of length `p-1`.
A function on that scheme can have `Theta(p)` zeros, and Strassmann's theorem
cannot help without a separately constructed analytic function of uniformly
bounded Weierstrass degree.

### 5.3 “Non-ordinary divisor” needs four missing inputs

Even the identification “trace zero modulo `p` = non-ordinary fiber” is not
automatic.  It requires:

- **[GAP-CRYS]** an integral overconvergent `F`-crystal companion and lattice
  for the rank-two Mellin object, with its trace reducing to the displayed
  `b_r`;
- **[GAP-HODGE]** a determinant/Hodge-polygon calculation proving that trace
  cancellation modulo `p` is equivalent to a Newton-polygon jump;
- **[GAP-FAMILY]** one rigid/formal parameter family containing all tame
  twists `omega^{-r}`;
- **[GAP-DEGREE]** an `o(p)` bound for the degree of its Hasse or Weierstrass
  function.

The common heuristic that a Hasse locus on a curve has conductor-bounded
degree is false in defining characteristic without further normalization.
For an elliptic family the Hasse invariant is a section of
`omega_Hodge^(p-1)`, so its divisor has degree
`(p-1)deg(omega_Hodge)`.  The Legendre Hasse polynomial has degree
`(p-1)/2`, an explicit `Theta(p)` example.  Even a genuine curve family would
therefore not by itself yield `o(p)`.

## 6. Explicit Beauville-IV correspondence

Write

\[
 E_u:\quad y^2+(1-2u)xy+u^2y=x^3,\qquad
 \iota(u)=\frac{1-8u}{8(1+u)}.                             \tag{6.1}
\]

On `E_u`,

\[
 P=(0,0),\quad T=(-u^2,-u^3),\quad Q=P+T=(u,-u),           \tag{6.2}
\]

have orders `3`, `2`, and `6`, respectively.  The cyclic kernel generated by
`Q` has nonzero points with `x`-coordinates

\[
 0,\quad -u^2,\quad u.                                    \tag{6.3}
\]

Applying Velu's formulas gives a rational degree-six map

\[
 E_u\longrightarrow C_u=E_u/\langle Q\rangle.             \tag{6.4}
\]

Explicitly,

\[
 C_u:\ y^2+(1-2u)xy+u^2y=x^3+A_4x+A_6,
\]

\[
 A_4=-5u^4-5u^3-20u^2-5u,
\quad
 A_6=3u^6+7u^5-41u^4-24u^3-14u^2-u.                     \tag{6.5}
\]

The Velu map is `(x,y) -> (X,Y)` with

\[
 X=\frac{N_X}{x^2(x-u)^2(x+u^2)},\qquad
 Y=\frac{N_Y}{x^3(x-u)^3(x+u^2)^2},                      \tag{6.6}
\]

where

\[
\begin{aligned}
N_X={}&x^6+(u^2-2u)x^5+(u^4-u^3+5u^2+u)x^4\\
 &+(-2u^5+7u^4-2u^3)x^3+(5u^6-6u^5+u^4)x^2\\
 &+(-4u^7+2u^6)x+u^8,
\end{aligned}                                             \tag{6.7}
\]

and

\[
\begin{aligned}
N_Y={}&x^8y +(2u^2-3u)x^7y
 +(2u^5+u^4+7u^3-2u^2-u)x^7\\
&+(-7u^3-u^2-u)x^6y
 +(u^7-6u^6+24u^5-18u^4+2u^3)x^6\\
&+(-7u^4-u^3-u^2)x^5y
 +(-3u^8+38u^7-54u^6+22u^5-4u^4)x^5\\
&+(-20u^6+10u^5-5u^4)x^4y
 +(15u^9-67u^8+55u^7-17u^6+u^5)x^4\\
&+(-8u^8+26u^7-14u^6+u^5)x^3y
 +(-25u^{10}+56u^9-30u^8+4u^7)x^3\\
&+(12u^9-19u^8+4u^7)x^2y
 +(19u^{11}-24u^{10}+6u^9)x^2\\
&+(-8u^{10}+5u^9)xy+(-7u^{12}+4u^{11})x
 +2u^{11}y+u^{13}.
\end{aligned}                                             \tag{6.8}
\]

Let `w^2=-3`, set `v=iota(u)`, and put

\[
\begin{aligned}
 q&=\frac{4(u+1)w}{3},
 &r&=-\frac{13u^2+5u+1}{3},\\
 s&=\frac{4uw+2u+w-1}{2},
 &h&=-\frac{64u^3w+78u^3+48u^2w-15uw-9u+w-3}{18}.
\end{aligned}                                             \tag{6.9}
\]

The coordinate change

\[
 X_C=q^2X_v+r,\qquad
 Y_C=q^3Y_v+q^2sX_v+h                                    \tag{6.10}
\]

identifies `C_u` with `E_v` over `Q(w)(u)`.  All formulas (6.2)--(6.10),
including the target equation after substitution, are checked symbolically
in `codex_lt_isogeny.py`.
The exact comparison is

\[
 C_u\simeq E_{\iota(u)}^{(-3)}\quad\text{over }\mathbf Q(u),
 \qquad
 C_u\simeq E_{\iota(u)}\quad\text{over }\mathbf Q(\sqrt{-3})(u). \tag{6.11}
\]

The constant quadratic twist is load-bearing at rank two and corrects the
untwisted wording in Q6394.  It disappears after symmetric square:

\[
 \operatorname{Sym}^2(V\otimes\chi_{-3})
   \simeq\operatorname{Sym}^2(V).                          \tag{6.12}
\]

Thus (6.11) supplies precisely the Frobenius-compatible correspondence needed
for the rank-three descent.  At five good split primes, the script evaluates
the rational map and the coordinate change at every eligible affine point.

## 7. CFVZ hypergeometric convention

CFVZ define the Franel series `h` and state that it satisfies a second-order
equation.  They identify it, after rational pullback and gauge, with

\[
 {}_2F_1\left(\frac13,\frac23;1;y\right),
\]

and express the Apéry period by

\[
 f_A(\phi(x))=(1+x)h(x)^2.                                \tag{7.1}
\]

They do not name the Apéry period as the literal
`3F2(1/3,2/3,1;1,1)`.  Indeed, coefficientwise,

\[
 {}_3F_2\left(\frac13,\frac23,1;1,1;z\right)
 ={}_2F_1\left(\frac13,\frac23;1;z\right),                \tag{7.2}
\]

because one upper `1` cancels one lower `1`.  It remains rank two.  A genuine
Clausen symmetric-square presentation is

\[
 {}_2F_1\left(\frac13,\frac23;1;y\right)^2
 ={}_3F_2\left(\frac13,\frac23,\frac12;1,1;4y(1-y)\right). \tag{7.3}
\]

Hence the adjudication is:

- rank two: the Franel elliptic / arithmetic `2F1` object;
- rank three: the Apéry object, obtained from its symmetric square and a
  rational pullback;
- the suspected literal `3F2(...,1;1,1)` identification is false as a
  rank-three identification.

## 8. Verification ledger

All commands were run from the repository root.

```text
$ python3 problems/3.2/research/scripts/codex_lt_galois_action.py
VERIFIED p=13: exact Q(zeta_12) Galois action sigma_a(M(r))=M(ar) for 44 (a,r) pairs; generator=2, smooth Mellin fibers=10, character checks=528
VERIFIED p=17: exact Q(zeta_16) Galois action sigma_a(M(r))=M(ar) for 120 (a,r) pairs; generator=3, smooth Mellin fibers=14, character checks=1920
VERIFIED exact cyclotomic Galois action at p=13,17: 164 Mellin identities and 2448 point-character identities

$ python3 problems/3.2/research/scripts/codex_lt_orbit_product.py
VERIFIED p=13: exact cyclotomic Galois action, integral orbit products, and p-adic hit transfer
VERIFIED p=13: grouped Parseval identity; ungrouped RHS differs by 13440
VERIFIED p=17: exact cyclotomic Galois action, integral orbit products, and p-adic hit transfer
VERIFIED p=17: grouped Parseval identity; ungrouped RHS differs by 207360

$ python3 problems/3.2/research/scripts/codex_lt_tadic_obstruction.py
VERIFIED tame exponent obstruction: p^j=1 mod (p-1), so zeta^(p^j)=zeta != 1 although p^j tends p-adically to 0
VERIFIED Legendre Hasse polynomial degree=(p-1)/2 for p=5,7,11,13,17,19,23,29,31

$ python3 problems/3.2/research/scripts/codex_lt_cfvz_rank_check.py
VERIFIED literal 3F2(1/3,2/3,1;1,1) coefficient cancellation through degree 29: it is the rank-2 2F1 series
VERIFIED Clausen square through degree 23: 2F1(1/3,2/3;1;y)^2 = 3F2(1/3,2/3,1/2;1,1;4y(1-y))
VERIFIED CFVZ Apéry/Franel pullback f_A(phi(x))=(1+x)h(x)^2 through degree 13

$ python3 problems/3.2/research/scripts/codex_lt_isogeny.py
VERIFIED exact Franel model and deck involution phi(iota(u))=phi(u)
VERIFIED Q=(u,-u) has exact order 6 and kernel x-coordinates 0,u,-u^2
VERIFIED exact formulas for the cyclic 6-isogeny E_u -> C_u
VERIFIED C_u becomes E_iota(u) over Q(sqrt(-3))(u); twist class is -3
VERIFIED p=7: rational map lands on E_iota(u) at 12 affine points
VERIFIED p=13: rational map lands on E_iota(u) at 84 affine points
VERIFIED p=19: rational map lands on E_iota(u) at 228 affine points
VERIFIED p=31: rational map lands on E_iota(u) at 732 affine points
VERIFIED p=37: rational map lands on E_iota(u) at 1092 affine points
VERIFIED explicit isogeny at five primes (2148 point evaluations total)
```

## 9. Sources checked

1. C. Liu and D. Wan,
   [*T-adic exponential sums over finite fields*](https://doi.org/10.2140/ant.2009.3.489),
   Definitions 1.1--1.2 and Theorems 1.4--1.5, 2.3, 2.9.
2. C. Liu and C. Niu,
   [*Generic twisted T-adic exponential sums of polynomials*](https://arxiv.org/abs/0911.4213),
   Definition 1.1 and Theorems 1.9--1.12.
3. C. Niu, W. Liu, and C. Liu,
   [*On L-functions of twisted T-adic exponential sums*](https://doi.org/10.1016/j.jnt.2017.10.016).
4. C. Davis, D. Wan, and L. Xiao,
   [*Newton slopes for Artin--Schreier--Witt towers*](https://arxiv.org/abs/1310.5311),
   especially Theorem 1.2.
5. X. Caruso, F. Fuernsinn, D. Vargas-Montoya, and W. Zudilin,
   [*Galois Groups of Apéry-like Series Modulo Primes*](https://arxiv.org/abs/2510.23298),
   especially pp. 2--4 and p. 9 of v1.
6. D. Zagier,
   [*Integral solutions of Apéry-like recurrence equations*](https://people.mpim-bonn.mpg.de/zagier/files/tex/AperylikeRecEqs/fulltext.pdf),
   Section 7, for Beauville family IV and its `Gamma_1(6)` model.  The explicit
   cyclic degree-six map in Section 6 above is a new Velu calculation checked
   by the accompanying script; it is not attributed to CFVZ or Zagier.

## Least-confident step

**[GAP-CRYS] is the least-confident step.**  The integral `ell`-adic descent,
the explicit degree-six isogeny (including its constant quadratic twist), and
all finite-field trace congruences are concrete.  What has not been
constructed is one integral overconvergent rank-two Mellin `F`-crystal with a
lattice whose reduction realizes every `b_r` and whose Hodge data turn trace
cancellation into a Newton-polygon jump.  Without that object, “non-ordinary
twist” is a heuristic renaming of the original divisibility event, not a
theorem to which a Hasse-locus degree bound can be applied.
