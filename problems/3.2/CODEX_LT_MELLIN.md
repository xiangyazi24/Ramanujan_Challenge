# The defining-characteristic Mellin frontier: orbit norms and the T-adic audit

## Executive verdict

No power saving for

\[
 Z_p=\{1\le r\le p-2:\mathfrak p\mid \widetilde M(r)-\widetilde T(r)\}
\]

follows from either proposed route.

- **[VERIFIED-GALOIS-M]** The character convention is covariant:
  `sigma_a(chi_r(t))=chi_{ar}(t)`.  Exact cyclotomic computations at
  `p=13,17` verify the action on `M` for every `r` and every automorphism.
  The identical statement for `T` follows if its pointwise characteristic-zero
  traces are rational integers.  The intended Apéry/K3 construction has this
  normalization, but the inherited phrase “integral compatible system” alone
  only implies algebraic integrality; this documentary issue is
  **[GAP-LT-RATIONALITY]**.
- **[CONDITIONAL-NORM]** Under that rationality input, the product over one
  Galois orbit is an integer.  If `k` members vanish modulo the same chosen
  `mathfrak p` and the product is nonzero, then `p^k` divides that integer.
  Even granting this strongest coefficient-field case, the resulting bound
  is vacuous.
- **[NEGATIVE-NORM]** The sharp available Archimedean scale is
  `p^(3/2)` per conjugate.  It yields `k <= (3/2+o(1))|O|`, which is weaker
  than `k<=|O|`.  Orbitwise AM--GM becomes nontrivial only when the orbit
  root-mean-square is `<p`, a factor `sqrt(p)` below the natural scale.
- **[CORRECTION-PARSEVAL]** The Parseval formula proposed in the task omits
  cross terms inside the degree-two fibers of `phi`.  The exact right side
  uses the fiber sums `B_t=sum_{phi(x)=t}a_{p,x}^2`.  This correction makes
  the second moment larger, not smaller.
- **[GAP-LT-EXACT-ZERO]** If an orbit product is zero, the whole orbit consists
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

## 1. The Galois action and the conditional integer orbit product

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

For `T`, the same conclusion requires the coefficient-field assertion

\[
 \tau_{p,t}:=\operatorname{Tr}(\operatorname{Frob}_t\mid\mathcal A_-)
 \in\mathbf Z.                                             \tag{1.3a}
\]

It is compatible with the intended construction: the Apéry rank-three object
is defined over `Q`, and the quadratic twist has values in `{0,+1,-1}`.
However, Q6394 states only “integral compatible system”; algebraic integrality
without `Q`-rationality would give `tau_{p,t} in O_E` for a coefficient field
`E`, not necessarily an integer.  Until the rational coefficient field and
all middle-extension stalks are written explicitly, this is
**[GAP-LT-RATIONALITY]**.  Granting (1.3a),

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

This remains true if `Delta_r` happens to lie in a proper cyclotomic
subfield: (1.5) then repeats the smaller-field norm by a fixed multiplicity.
If instead the coefficient field is `E != Q`, the same cyclotomic product is
only in `O_E`; applying `N_{E/Q}` recovers an integer but multiplies the
Archimedean exponent by `[E:Q]`.  Thus the coefficient-field caveat can only
weaken the method and does not affect the negative conclusion below.

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

For the generic nonexceptional characters, the literal `a_{p,x}^2` lift in
the task is the rank-eight tensor-square pushforward, while `T` has Mellin
dimension four.  The inherited local-monodromy ledger therefore gives

\[
 |\widetilde M(r)|\le8p^{3/2},\qquad
 |\widetilde T(r)|\le4p^{3/2},
\]

hence one may take `c=12` in the direct presentation.  The conductor ledger
gives the coarser uniform constant `35` when all middle-extension corrections
are retained.  The rank-two `A_+` Mellin trace is only the *residual* survivor:
the genuine characteristic-zero difference also contains the Tate summand
`phi_* det(F)`, whose trace vanishes modulo `p` but not in characteristic
zero.  Generically its nontrivial Mellin part has dimension two, so the exact
descent can at best improve the constant to

\[
 |\Delta_r|\le4p^{3/2}.                                   \tag{2.1}
\]

Any of these constants has the same fatal exponent.  Writing `m=|O_r|`,
(1.8) and the bound `|Delta_s|<=c p^(3/2)` give

\[
 k\le \log_p|P_r|
 \le m\left(\frac32+\log_p c\right).                      \tag{2.2}
\]

Even with `c=4`, the right side is larger than `m`, so (2.2) is weaker than
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

The natural and observed root-mean-square is of order `p^(3/2)`.
No orbit size changes the exponent in (2.3): small orbits have too little
averaging, while an orbit of size comparable to `p` still has mean square of
order `p^3`.  The available global second moment does not force (2.4) for a
single orbit.  The precise worst-case calculation appears in Section 3.
Hence **no orbit structure `q|p-1` is presently nonvacuous**.

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

For the literal difference in the target, write

\[
 \tau_t=\operatorname{Tr}(\operatorname{Frob}_t\mid\mathcal A_-),
 \qquad D_t=B_t-\tau_t,
 \qquad \Delta_r=\sum_tD_t\chi_r(t).                       \tag{3.3}
\]

Under the rationality normalization in Section 1, the degree-two fibers and
Hasse bound give `|B_t|<=8p`; rank three and purity give `|tau_t|<=3p`.
Hence `|D_t|<=11p` and exact full-character Parseval gives

\[
 \sum_{r=0}^{N-1}|\Delta_r|^2
   =N\sum_t|D_t|^2 \le121N^2p^2<121p^4.                   \tag{3.4}
\]

For `r=1,...,N-1`, the exact left side is the right side of (3.4) minus
`|sum_t D_t|^2`.  If an orbit has size `m`, the global bound alone yields

\[
 \sqrt{S_O/m}\le {11Np\over\sqrt m},\qquad
 {k\over m}\le\log_p\left({11Np\over\sqrt m}\right).       \tag{3.5}
\]

If `m=p^theta`, the latter is `2-theta/2+o(1)`.  Even for the largest
possible `m` comparable to `p`, it is `3/2+o(1)`.  Equivalently, forcing the
nontrivial threshold `S_O/m<p^2` from (3.4) would require `m>121N^2`, whereas
`m<=N`.
Thus Parseval makes **no orbit size or divisor structure nonvacuous**.  A
useful orbit argument would require the genuinely new orbitwise input

\[
 S_O<mp^2,                                                  \tag{3.6}
\]

marked **[GAP-LT-ORBIT-ENERGY]**; the global `O(p^4)` moment does not imply
it.

## 4. The zero orbit product cannot be removed by purity

If `P_r=0`, then one `Delta_s` is zero as an algebraic number.  Galois
covariance implies

\[
 \Delta_s=0\quad\text{for every }s\in O_r.                 \tag{4.1}
\]

There is no weight mismatch to exclude (4.1).  Both `M` and `T` come from
weight-two input trace functions followed by `H_c^1`, so their generic
Mellin eigenvalues and their difference have weight three.  After reduction,
the surviving `A_+` block has two-dimensional weight-three Mellin cohomology;
the genuine characteristic-zero difference also has the Tate-pushforward
block described above.  Already a trace of two nonzero Weil numbers can be
exactly zero, so neither presentation excludes exact cancellation.

Katz's complex equidistribution may predict that exact trace zero has measure
zero in an appropriate extension-field limit, but it supplies no uniform
bound for these characteristic-zero zero orbits as `p` and the tame character
group vary diagonally.  This is **[GAP-LT-EXACT-ZERO]**.  It is logically prior
to applying a nonzero integer norm bound and could contain an orbit of size
`phi(q)` comparable to `p`.

## 5. Vertical family / T-adic audit

### 5.1 Exact literature verdict: `T` is the wrong variable

The named papers were checked at their definitions and main theorems.

- Liu--Wan, Definitions 1.1--1.2, take a locally constant **additive**
  character `psi:Z_p -> C_p^times` of order `p^m` and specialize
  `T=psi(1)-1`.  Theorems 1.4--1.5 and 2.9 concern analytic continuation,
  the Hodge bound, and ordinarity transfer across those additive
  specializations.
- Davis--Wan--Xiao take finite characters of the additive Galois group
  `Z_p` of an Artin--Schreier--Witt tower.  Theorem 1.2 varies their
  `p`-power conductor and proves eventual slope progressions.
- Ren--Wan--Xiao--Yu replace `Z_p` by the additive group
  `Z_{p^ell}`, topologically `Z_p^ell`.  Their universal character has
  variables `T_1,...,T_ell`; Theorem 1.4 is an
  `I=(T_1,...,T_ell)`-adic Hodge/leading
  term result and Theorem 1.6 is a spectral-halo statement on an admissible
  locus.  These are several wild variables, not a tame one.
- In Liu--Niu's twisted theory, `chi=omega^{-u}` is fixed and `T` still
  varies the additive character.  The Hasse polynomial `H_u` is constructed
  separately for each fixed `u`.
- Ren makes the separation explicit: Definition 1.2 defines
  `L_f(omega^u,chi,s)` with fixed `u` and additive
  `chi:Z_p -> C_p^times`; Theorem 1.6 varies only the conductor of `chi`.
  Definition 2.7 specializes `T=chi(1)-1`.  Lemma 2.8 packages all `u` only
  through the discrete product

  \[
   \prod_{u=0}^{q-2}L_f(\omega^u,T,s)
       =L_{f(x^{q-1})}(1,T,s),                             \tag{5.1}
  \]

  not through analytic interpolation in `u`.  For `q=p`, the degree of
  `f(x^{p-1})` and the associated Hodge/conductor complexity are already of
  order `p`.

Our characters `chi_r:F_p^times -> mu_{p-1}` are purely tame.  None of these
theorems compares them, so the exact verdict is **[NEGATIVE-TADIC]**: the
checked T-adic/Newton--Hodge results give neither `o(p)` nor `O(polylog p)`.

### 5.2 The tame base is `p-1` disconnected points

The correct integral tame-character base at a fixed prime is

\[
 X_p=\operatorname{Spec}\mathbf Z_p[U]/(U^{p-1}-1)
     \simeq\coprod_{r\bmod p-1}\operatorname{Spec}\mathbf Z_p, \tag{5.2}
\]

because the Teichmüller roots are distinct modulo `p`.  Its ring of functions
is the product of `p-1` copies of `Z_p`; a section may vanish on an arbitrary
subset of the components.  On a reduced zero-dimensional base a nonempty
zero locus is a union of connected components, not a bounded-degree Cartier
divisor.

The same obstruction is visible in the full multiplicative weight space:
`Hom_cont(Z_p^times,C_p^times)` is a disjoint union of `p-1` open disks,
indexed by the restriction to `mu_{p-1}`.  The pure tame characters
`omega^{-r}` are the wild-trivial center points of different disks.  Halo
theorems work inside one fixed disk and do not compare these centers.

There is also a direct continuity obstruction.  For a nontrivial tame root
`zeta`, if `n -> zeta^n` extended continuously to `Z_p`, then `p^j -> 0`
would imply `zeta^{p^j} -> 1`; instead

\[
 p^j\equiv1\pmod{p-1},\qquad \zeta^{p^j}=\zeta\ne1.        \tag{5.3}
\]

Thus no connected rigid reparametrization joins all `r`.  The disconnected
package (5.2) exists, but a general Hasse tuple on it has `p-1` independent
coordinates and only the tautological `O(p)` zero bound.

### 5.3 Trace divisibility is not automatically non-ordinarity

For an integral rank-two Frobenius module with characteristic polynomial

\[
 X^2-BX+D,\qquad v_p(D)=h>0,                               \tag{5.4}
\]

and Hodge slopes `{0,h}`, ordinarity means a unique unit root.  In this
narrow normalization, `B` is a unit exactly when that unit root exists, so
`p|B` is equivalent to non-ordinarity.  Every hypothesis is load-bearing:

- `diag(1,-1)` has two unit roots and trace zero but is ordinary for Hodge
  slopes `{0,0}`;
- `diag(p,p^2)` is ordinary for Hodge slopes `{1,2}` and its trace is
  automatically divisible by `p`;
- in higher rank, two unit eigenvalues may cancel in the full trace.

Q6394 supplies an `ell`-adic rank-two generic Mellin cohomology trace after
residual cancellation.  It does not yet print an integral overconvergent
`F`-crystal, stable lattice, determinant valuation, or Hodge polygon.  Hence
the asserted equivalence between `p|b_r` and non-ordinarity of that rank-two
Mellin block is **[GAP-LT-TRACE-ORD]**, not an inherited consequence of
complex weight three.

Even for a genuine family over a curve, “the Hasse locus is a divisor” gives
no uniform bounded degree in defining characteristic.  The Legendre Hasse
polynomial has degree `(p-1)/2`; equivalently the elliptic Hasse invariant is
a section of `omega_Hodge^(p-1)`.  Thus a curve by itself would still permit
linearly many nonordinary points.

### 5.4 Explicit computation plan for the rank-three objects

There is a useful computation, but it tests a new conjecture rather than
invoking an existing T-adic theorem.

1. Construct integral overconvergent `F`-crystals for `A_-` (conductor at
   most 11) and for the residual `A_+` block.  Print their connections,
   Frobenius matrices, stable lattices, local residues, and determinants.
   Rigid/log-de Rham cohomology must reproduce generic Mellin ranks four for
   `A_-` and two for `A_+`; the finite exceptional low-order characters and
   their `H_c^0/H_c^2` corrections must be enumerated separately.  For the
   generic `A_+` block, the desired calculation is

   \[
    \det(1-X\Phi_r)=1-B_rX+p^3\epsilon_rX^2,
    \quad \epsilon_r\in\mathbf Z_p^\times,
    \quad B_r\bmod p=b_r,                                 \tag{5.5}
   \]

   together with Hodge slopes `{0,3}`.  Equation (5.5) and those slopes are
   **[GAP-LT-CRYS]**; they are the precise input needed to close
   **[GAP-LT-TRACE-ORD]**.
2. Pull the crystal back along the Lang cover
   `L_p:G_m -> G_m`, `y -> y^(p-1)`, whose deck group is `F_p^times`.
   Projection and the integral idempotents

   \[
    e_r={1\over p-1}\sum_{a\in\mathbf F_p^\times}
          \omega^r(a)[a]                                  \tag{5.6}
   \]

   split the global rigid cohomology into the individual Kummer-twist
   blocks.  Compute Frobenius on the pullback to precision at least `p^4`,
   apply a Teichmüller DFT to the deck action, and compare every block trace
   modulo `p` with the existing `b_r` table.
3. Compute each Newton polygon.  If (5.5) holds, put
   `delta_r=NP_r(1)-HP_r(1)`.  A hit has `delta_r>=1`; hence

   \[
    \#Z_p\le\sum_r\delta_r.                                \tag{5.7}
   \]

   This makes the total Newton--Hodge excess a concrete experimental target.

The pullback has degree `p-1` and total cohomological rank linear in `p`.
Standard Dwork, Grothendieck--Ogg--Shafarevich, and Hodge estimates therefore
give only `sum_r delta_r=O(p)`, exactly as Ren's product (5.1) acquires degree
linear in `p`.

### 5.5 The theorem that would actually finish the route

What is missing is not another fixed-component halo theorem but a
cross-component residual theorem.  A sufficient statement is:

> **[GAP-LT-TAME-COUNT]** For a fixed integral bounded-conductor
> `F`-isocrystal `A`, there is `delta>0` such that, uniformly at good primes,
> the number of tame `chi:F_p^times -> mu_{p-1}` for which
> `Tr(Frob | H_c^1(G_m,A tensor L_chi))` is zero modulo `p` is
> `O_A(p^(1-delta))`.

Equivalently for the computation above, one could prove
`sum_r delta_r=o(p)`.  None of Liu--Wan, Davis--Wan--Xiao,
Ren--Wan--Xiao--Yu, Liu--Niu, or Ren contains such a conclusion.  Constructing
a relative Kummer--Mellin crystal over (5.2) would only place the same `p-1`
blocks side by side; without **[GAP-LT-TAME-COUNT]** it is tautological.

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

Moreover the coordinate scale in (6.9) satisfies

\[
 -{q^2\over6}={8\over9}(1+u)^2=:g(u),\qquad
 g(u)g(\iota(u))=1.                                      \tag{6.13}
\]

Thus the degree-six normalized symmetric square realizes exactly the
connection cocycle from Q6394.  This supplies the generic-fiber
Frobenius-compatible correspondence needed for the rank-three descent.  Its
integral/crystalline normalization at all bad middle-extension stalks is not
checked by these rational maps and remains **[GAP-LT-BAD-STALKS]**.  At five
good split primes, the script evaluates the rational map and coordinate
change at every eligible affine point.
At both split and nonsplit primes it separately verifies the exact trace law

\[
 a_p(E_{\iota(u)})=\left(\frac{-3}{p}\right)a_p(E_u),      \tag{6.14}
\]

so the traces themselves need not be equal; their squares always are.  This
also corrects any untwisted “same Frobenius trace” formulation of the Fricke
pair elsewhere in the working notes.

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
VERIFIED tame character special fiber splits into p-1 independent points at p=5,7,13,17,29
VERIFIED trace mod p is not a general ordinarity test; the unique-unit-root Hodge/determinant hypotheses are load-bearing
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
VERIFIED normalized Sym^2 scalar g(u)=8(1+u)^2/9 and g(u)g(iota(u))=1
VERIFIED p=7: rational map lands on E_iota(u) at 12 affine points
VERIFIED p=13: rational map lands on E_iota(u) at 84 affine points
VERIFIED p=19: rational map lands on E_iota(u) at 228 affine points
VERIFIED p=31: rational map lands on E_iota(u) at 732 affine points
VERIFIED p=37: rational map lands on E_iota(u) at 1092 affine points
VERIFIED explicit isogeny at five primes (2148 point evaluations total)
VERIFIED a_p(E_iota(u))=(-3|p)*a_p(E_u) at p=5,7,11,13,17 (38 smooth fibers; split and nonsplit twists)
```

## 9. Sources checked

1. C. Liu and D. Wan,
   [*T-adic exponential sums over finite fields*](https://doi.org/10.2140/ant.2009.3.489),
   Definitions 1.1--1.2 and Theorems 1.4--1.5, 2.3, 2.9.
2. C. Liu and C. Niu,
   [*Generic twisted T-adic exponential sums of polynomials*](https://arxiv.org/abs/0911.4213),
   Definition 1.1 and Theorems 1.4 and 1.10.
3. C. Davis, D. Wan, and L. Xiao,
   [*Newton slopes for Artin--Schreier--Witt towers*](https://arxiv.org/abs/1310.5311),
   Theorem 1.2 and Section 4.
4. R. Ren, D. Wan, L. Xiao, and M. Yu,
   [*Slopes for higher rank Artin--Schreier--Witt towers*](https://arxiv.org/abs/1605.02254),
   Theorems 1.1, 1.4, 1.6 and Section 2.3.
5. R. Ren,
   [*Newton slopes for twisted Artin--Schreier--Witt towers*](https://arxiv.org/abs/1704.07017),
   Definition 1.2, Theorem 1.6, Definition 2.7, and Lemma 2.8.
6. R. Liu, D. Wan, and L. Xiao,
   [*The eigencurve over the boundary of weight space*](https://arxiv.org/abs/1412.2584),
   Introduction pp. 2--4 and Theorems 1.3 and 1.5; the tame component is fixed.
7. X. Caruso, F. Fuernsinn, D. Vargas-Montoya, and W. Zudilin,
   [*Galois Groups of Apéry-like Series Modulo Primes*](https://arxiv.org/abs/2510.23298),
   especially pp. 2--4 and p. 9 of v1.
8. D. Zagier,
   [*Integral solutions of Apéry-like recurrence equations*](https://people.mpim-bonn.mpg.de/zagier/files/tex/AperylikeRecEqs/fulltext.pdf),
   Section 7, for Beauville family IV and its `Gamma_1(6)` model.  The explicit
   cyclic degree-six map in Section 6 above is a new Velu calculation checked
   by the accompanying script; it is not attributed to CFVZ or Zagier.

## Least-confident step

**[GAP-LT-CRYS] is the least-confident step.**  The generic-fiber degree-six
isogeny, its constant quadratic twist, and the displayed finite-field trace
checks are concrete.  What has not been constructed is one integral
overconvergent rank-two Mellin `F`-crystal with a stable lattice whose
reduction realizes every `b_r` and whose Hodge data turn trace cancellation
into a Newton-polygon jump.  The coefficient-field issue
**[GAP-LT-RATIONALITY]** and bad-stalk normalization
**[GAP-LT-BAD-STALKS]** are part of that arithmetic upgrade.  Without it,
“non-ordinary twist” is a heuristic renaming of the original divisibility
event, not a theorem to which a Hasse-locus argument can be applied.
