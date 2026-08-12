# Transverse Apéry common primes: a level-six Fricke/Sturm audit

## Executive verdict

**[THEOREM — modular reconstruction]** With

\[
t(\tau)=\left(\frac{\eta(\tau)\eta(6\tau)}{\eta(2\tau)\eta(3\tau)}\right)^{12},
\qquad
E(\tau)=\frac{\eta(2\tau)^7\eta(3\tau)^7}
{\eta(\tau)^5\eta(6\tau)^5}=F(t(\tau)),
\]

and

\[
H=\frac{q\,dt/dq}{t},\qquad \Psi=EH,
\]

one has

\[
H^2=(1-34t+t^2)E^2,
\qquad
\Psi=E^2\sqrt{1-34t+t^2},
\qquad
g(t)=\frac1{\Psi(\tau)}.
\]

The square root is the branch with constant term `1`.  The four cusps of
`Gamma_0(6)` are `infinity, 0, 1/2, 1/3`; their widths are `1,6,3,2`.
The cusp divisors are

\[
\operatorname{ord}(t)=(1,1,-1,-1),\qquad
\operatorname{ord}(E)=(0,0,1,1),\qquad
\operatorname{ord}(H)=(0,0,0,0).
\]

Under the normalized Fricke involution `W_6`,

\[
t:+1,\qquad E:-1,\qquad H:+1,\qquad \Psi:-1.
\]

Thus

\[
\Omega_m:=\Psi t^{-m}
\]

is a weight-four Fricke-odd weakly holomorphic form with cusp orders

\[
(-m,-m,m+1,m+1).
\]

**[THEOREM — coefficient interpretation]** The Apéry coordinate is literally
one Fourier coefficient of the Laurent grid:

\[
[q^0]\Omega_m=b_m.
\]

Moreover, if

\[
\mathcal X_r:=\Omega_0+5\sum_{m=1}^r g_m\Omega_{m-1},
\]

then

\[
[q^0]\mathcal X_r=-\Xi_r.
\]

**[THEOREM — exact transverse Cartier form]** For every prime `p>=5` and
`0<=j<p`, the level-six Hecke/Laurent-grid argument gives

\[
U_p\Omega_j\equiv b_j\Omega_0\pmod {p^3}
\]

coefficientwise in `Z((q))`.  Consequently, for `p>r`,

\[
U_p\mathcal X_r\equiv-\Xi_r\Omega_0\pmod {p^3}.
\]

Hence the transverse common condition is exactly two first-order Cartier
kernel conditions modulo `p`:

\[
p\mid b_r,\ p\mid\Xi_r
\quad\Longrightarrow\quad
U_p\Omega_r\equiv U_p\mathcal X_r\equiv0\pmod p.
\]

This is stronger and more geometric than saying that two unrelated integers
vanish, but it still does **not** produce a consecutive Fourier gap.

**[THEOREM — Sturm scale]** Since

\[
[\mathrm{SL}_2(\mathbf Z):\Gamma_0(6)]=12,
\]

the Sturm bound for integral weight `k` on `Gamma_0(6)` is

\[
B_6(k)=\left\lfloor\frac{k}{12}\,12\right\rfloor=k.
\]

The standard weight-four cusp form

\[
\Delta_6=\eta(\tau)^2\eta(2\tau)^2\eta(3\tau)^2\eta(6\tau)^2
\]

has order one at all four cusps and is Fricke-even.  Thus

\[
\Delta_6^r\Omega_r\in M_{4r+4}(\Gamma_0(6)),
\]

and a **new** theorem forcing a nonzero holomorphic form of weight `4r+4`
to have order at least `p-r` at infinity would imply

\[
p-r\le4r+4,
\qquad\boxed{p\le5r+4}.
\]

This is the precise place where the attractive constant `5` can arise from a
level-six Sturm calculation.

**[REMAINING GAP — load bearing]** The common condition does not currently
force that `p-r` gap, or any comparable consecutive gap.  It forces vanishing
of coefficients in a `p`-dissection (equivalently Cartier kernel), whereas
Sturm controls consecutive coefficients from a cusp.  For `p>r`, the
principal parts of `Omega_r` and `mathcal X_r` contain no negative exponent
divisible by `p`, so `U_p` simply discards those principal parts.  There is no
known theorem converting the two Cartier-kernel relations into `p-r`
consecutive zero Fourier coefficients of a nonzero bounded-weight form.

There is an even more elementary obstruction to the naive pole-clearing
argument.  The weight-two form

\[
C:=Et
\]

has cusp orders `(1,1,0,0)`, so it is the minimal holomorphic clearer of the
two poles.  But

\[
C^r\Omega_r=E^{r+1}H
\]

has constant term `1`.  Likewise `Delta_6^r Omega_r` has constant term `1`.
Thus `p|b_r` does not make either naturally cleared form begin with zero; the
coordinate `b_r` sits at `q^0` **before** clearing, not at the leading
coefficient after clearing.

**[FINITE EVIDENCE — not a theorem]** The exact characteristic-zero
factorization supplied in the project through `r=10000` has only

\[
(p,r)=(17,13),\qquad (2237,492),
\]

among common pairs with `p>r`.  Both satisfy `p<=5r`, and there is no repeated
high common prime in that finite box.  This is compatible with the conjectural
support bound but does not prove it.

**Verdict.** I do not obtain an unconditional absolute linear-support theorem,
and I do not obtain an exact counterexample.  What survives hostile audit is a
sharp obstruction: **`C=5` is a valid Sturm consequence of one missing
Cartier-to-consecutive-gap (or equivalent pole-drop) theorem, but the actual
Fricke/Laurent identities presently give only sparse `p`-dissection
vanishing.**  The truth of an absolute `C` therefore remains open.  A natural
two-coordinate random model actually predicts that every fixed `C` should
eventually fail, so the linear-support statement, if true, requires a genuine
Apéry-specific deterministic exclusion beyond the current modular/Sturm
package.

---

## 1. Arithmetic normalization

The sequence is

\[
b_0=1,\qquad b_1=5,
\]

\[
(n+1)^3b_{n+1}
=(34n^3+51n^2+27n+5)b_n-n^3b_{n-1}.
\]

Equivalently its generating series satisfies

\[
\left[
\theta^3-t(2\theta+1)(17\theta^2+17\theta+5)
+t^2(\theta+1)^3
\right]F=0,
\qquad \theta=t\frac d{dt}.
\]

The first values are

\[
b_0=1,\quad b_1=5,\quad b_2=73,\quad b_3=1445,\quad b_4=33001.
\]

The transverse series is

\[
g(t)=\frac1{F(t)^2\sqrt{D(t)}}=\sum_{m\ge0}g_mt^m,
\qquad D(t)=1-34t+t^2,
\]

with the square root chosen to have constant term one.  The first coefficients
are

\[
g_0=1,\qquad g_1=7,\qquad g_2=192.
\]

Thus

\[
\Xi_0=-1,\qquad \Xi_1=-36,\qquad \Xi_2=-4836.
\]

### 1.1 Integrality of `g`

**[THEOREM]** `g(t)` lies in `Z[[t]]`.

**Proof.** The eta product for `t` has the form

\[
t(q)=q+q^2\mathbf Z[[q]],
\]

and `Psi(q)=1+q Z[[q]]`.  Therefore its compositional inverse satisfies

\[
q(t)=t+t^2\mathbf Z[[t]],
\]

recursively with integral coefficients.  Inverting `Psi(q)` in `Z[[q]]` and
composing with `q(t)` gives

\[
\frac1{\Psi(q(t))}\in\mathbf Z[[t]].
\]

Section 3 proves that this last series is exactly
`1/(F^2 sqrt(D))`.  Hence every `g_m` is integral. `square`

A useful endpoint consequence is immediate:

**[THEOREM — endpoint `p=5`]** For every `r>=0`,

\[
\Xi_r\equiv-1\pmod5.
\]

So `p=5` is never a transverse common prime.  For `p=2,3` and `p>r`, the only
possible rows are `r<2` or `r<3`; `b_0=1,b_1=5,b_2=73` exclude them.  Also
`r=0` is excluded by `b_0=1`.  Therefore every actual common pair with `p>r`
has `p>=7` and `r>=1`.

---

## 2. The exact level-six modular objects

Set `q=exp(2 pi i tau)`.  Define

\[
t=\left(\frac{\eta(\tau)\eta(6\tau)}
{\eta(2\tau)\eta(3\tau)}\right)^{12},
\]

\[
E=\frac{\eta(2\tau)^7\eta(3\tau)^7}
{\eta(\tau)^5\eta(6\tau)^5},
\]

\[
H=\frac{q\,dt/dq}{t},\qquad \Psi=EH.
\]

Beukers' modular parametrization is

\[
E(\tau)=F(t(\tau)).
\]

The eta quotient congruence conditions show that `t` is a weight-zero modular
function and `E` a weight-two modular form with trivial character on
`Gamma_0(6)`.

### 2.1 Cusps, widths, and eta orders

**[THEOREM]** `Gamma_0(6)` has index `12`, four cusps, no elliptic points, and
genus zero.  Representatives and widths are

| cusp | infinity | 0 | 1/2 | 1/3 |
|---|---:|---:|---:|---:|
| width | 1 | 6 | 3 | 2 |

For an eta quotient

\[
f=\prod_{\delta\mid6}\eta(\delta\tau)^{r_\delta},
\]

the Ligozat cusp-order formula at the cusp represented by `1/c` is

\[
\operatorname{ord}_{1/c}(f)=
\frac{6}{24c\gcd(c,6/c)}
\sum_{\delta\mid6}
\frac{r_\delta\gcd(c,\delta)^2}{\delta}.
\]

Applying it to the exponent vectors

\[
t:(12,-12,-12,12),\qquad
E:(-5,7,7,-5)
\]

gives

| object | infinity | 0 | 1/2 | 1/3 |
|---|---:|---:|---:|---:|
| `t` | 1 | 1 | -1 | -1 |
| `E` | 0 | 0 | 1 | 1 |

The logarithmic derivative has the explicit Eisenstein expression

\[
H=\frac12\bigl(E_2(\tau)+6E_2(6\tau)
-2E_2(2\tau)-3E_2(3\tau)\bigr).
\]

At a cusp where `t=c z^e(1+O(z))`, `e=+/-1`, one has
`dt/t=e dz/z+O(dz)`.  Thus the local weight-two coefficient of `H` has
nonzero constant term.  Hence

\[
\operatorname{ord}_{\mathfrak a}H=0
\]

at all four cusps.  Consequently

\[
\operatorname{ord}(\Psi)=(0,0,1,1)
\]

and

\[
\boxed{
\operatorname{ord}(\Omega_m)=(-m,-m,m+1,m+1).}
\]

### 2.2 Interior zeros and the quadratic differential identity

The polar divisor of `t` has degree two, so

\[
t:X_0(6)\longrightarrow\mathbf P^1
\]

has degree two.  The Fricke involution preserves `t`, and therefore `t` is a
coordinate on the quotient by `W_6`.  Its two ramification points are the two
interior Fricke fixed points.  `H` has a simple zero at each and no other
interior zero.  `E`, being an eta quotient, has no zeros in the upper
half-plane.  Hence `Psi` has exactly those two simple interior zeros in
addition to its cusp zeros at `1/2,1/3`.

**[THEOREM]**

\[
\boxed{H^2=D(t)E^2,\qquad D(t)=1-34t+t^2.}
\]

**Proof.** Both sides are holomorphic weight-four forms on `Gamma_0(6)` with
trivial character.  The eta product gives `t=q+O(q^2)`, `E=1+O(q)`, and
`H=1+O(q)`.  Expanding through the level-six weight-four Sturm bound `4`
checks equality.  Equivalently, divisor comparison shows that `H/E` is a
weight-zero square root of the degree-two branch polynomial, with the constant
term selecting the positive branch. `square`

Therefore

\[
H=E\sqrt D,
\qquad
\boxed{\Psi=E^2\sqrt D},
\]

and, using `E=F(t)`,

\[
\boxed{g(t)=\frac1{\Psi(\tau)}}.
\]

---

## 3. Fricke signs and the odd Laurent grid

Use the normalized Fricke slash

\[
(f|_kW_6)(\tau)=6^{-k/2}\tau^{-k}
 f\left(-\frac1{6\tau}\right).
\]

The eta transformation formula gives

\[
t|W_6=t,
\qquad
E|_2W_6=-E.
\]

Differentiating `t(-1/(6tau))=t(tau)` gives

\[
H|_2W_6=+H.
\]

Hence

\[
\boxed{\Psi|_4W_6=-\Psi},
\qquad
\boxed{\Omega_m|_4W_6=-\Omega_m}.
\]

### 3.1 Exact grid

For integers `M>=L>=0`, let `G(M,L)` be the space of weight-four
Fricke-odd weak forms whose cusp orders satisfy

\[
\operatorname{ord}_{\infty},\operatorname{ord}_0\ge-M,
\qquad
\operatorname{ord}_{1/2},\operatorname{ord}_{1/3}\ge L+1.
\]

**[THEOREM]**

\[
\boxed{
G(M,L)=\operatorname{span}\{\Omega_L,\Omega_{L+1},\ldots,\Omega_M\}.}
\]

**Proof.** Any Fricke-odd weight-four form vanishes at each Fricke fixed point:
in a local coordinate in which `W_6:z -> -z`, a weight-four form is a
quadratic differential `phi(z)(dz)^2`, and oddness forces
`phi(-z)=-phi(z)`.  Thus division by the simple fixed-point zeros of `Psi`
introduces no interior poles.  For `f in G(M,L)`, the quotient `f/Psi` is a
Fricke-invariant weight-zero function with possible poles only above
`t=0,infinity`; because the invariant function field is `C(t)`, the quotient
is a Laurent polynomial

\[
\sum_{m=L}^M c_m t^{-m}.
\]

Multiplication by `Psi` gives the claimed basis.  Distinct leading powers
`q^{-m}` give linear independence. `square`

The relevant cusp form is

\[
\Delta_6=\eta(\tau)^2\eta(2\tau)^2\eta(3\tau)^2\eta(6\tau)^2.
\]

Its cusp order is one at all four cusps, and

\[
\Delta_6|_4W_6=+\Delta_6.
\]

Thus it is not a hidden odd grid direction.

### 3.2 The holomorphic weight-four companion audit

**[THEOREM]** `dim M_4(Gamma_0(6))=5` and `dim S_4=1`.  A standard Eisenstein
basis is `E_4(d tau)`, `d|6`, and the cusp direction is `Delta_6`.
The Fricke action on the four old Eisenstein forms is

\[
E_4(d\tau)|W_6=\frac{36}{d^4}E_4((6/d)\tau).
\]

Hence a convenient basis of the Fricke-odd holomorphic Eisenstein subspace is

\[
B_1=E_4(\tau)-36E_4(6\tau),
\]

\[
B_2=4E_4(2\tau)-9E_4(3\tau).
\]

Their first two Fourier coordinates are

\[
([q^0]B_1,[q^1]B_1)=(-35,240),
\]

\[
([q^0]B_2,[q^1]B_2)=(-5,0).
\]

For `p>=7` these coordinates separate the two odd Eisenstein directions.  This
is also why a union-of-coordinate supports is not legitimate: the canonical
Eisenstein directions admit integral linear combinations with the small unit
coefficient `240`; divisibility can cancel in coordinates without forcing one
named companion individually to vanish.

---

## 4. The two transverse coordinates as Fourier residues

The identity `g=1/Psi` makes both coordinates intrinsic to the same Laurent
grid.

### 4.1 `b_m` is the constant term of `Omega_m`

**[THEOREM]**

\[
\boxed{[q^0]\Omega_m=b_m.}
\]

**Proof.** Since

\[
H=\frac{q}{t}\frac{dt}{dq},
\]

one has

\[
\frac{dq}{q}=\frac{dt}{tH}.
\]

Therefore

\[
\Omega_m\frac{dq}{q}
=EHt^{-m}\frac{dt}{tH}
=E t^{-m-1}dt
=F(t)t^{-m-1}dt.
\]

The residue at `q=0` of the left side is `[q^0]Omega_m`; the residue at `t=0`
of the right side is `b_m`. `square`

### 4.2 The exact Xi companion

Define, for `r>=0`,

\[
\boxed{
\mathcal X_r
:=\Omega_0+5\sum_{m=1}^r g_m\Omega_{m-1}.}
\]

The empty sum at `r=0` is allowed.

**[THEOREM]**

\[
\boxed{[q^0]\mathcal X_r=-\Xi_r.}
\]

Indeed

\[
[q^0]\mathcal X_r
=1+5\sum_{m=1}^r g_m b_{m-1}
=-\Xi_r.
\]

Moreover `mathcal X_r` is weight four, Fricke odd, integral as a Laurent
series, has poles only at `infinity,0`, and those poles have order at most
`r-1`.

This formulation is useful because it does not identify `r` with a modular
fiber.  Both `b_r` and `Xi_r` are coefficient/Mellin coordinates of actual
modular objects, exactly as required by the guard in the question.

---

## 5. Hecke descent and the exact Cartier interpretation

Let

\[
U_p\left(\sum_n a_nq^n\right)=\sum_n a_{pn}q^n,
\qquad
V_p\left(\sum_n a_nq^n\right)=\sum_n a_nq^{pn}.
\]

At weight four and `p` prime to `6`,

\[
T_p=U_p+p^3V_p.
\]

### 5.1 Laurent theorem in the range `j<p`

**[THEOREM]** For every prime `p>=5` and every `0<=j<p`,

\[
\boxed{U_p\Omega_j\equiv b_j\Omega_0\pmod {p^3}}
\]

coefficientwise in `Z((q))`.

**Proof.** `T_p` commutes with the Atkin-Lehner involutions.  At `infinity`,
`p^3V_p Omega_j` contributes the leading pole `p^3q^{-jp}`; the `U_p` term has
no negative exponent because `j<p`.  At `0` the same holds by `W_6`.  At
`1/2,1/3`, the positive order `j+1` is carried by `U_p` to order at least
`ceil((j+1)/p)=1`.  Thus

\[
T_p\Omega_j\in G(jp,0)
=\operatorname{span}\{\Omega_0,\ldots,\Omega_{jp}\}.
\]

Write

\[
T_p\Omega_j=\sum_{m=0}^{jp}C_m\Omega_m.
\]

The grid is integral and triangular because `Omega_m=q^{-m}+...`.  Reducing
modulo `p^3` replaces `T_p` by `U_p`.  Since `U_p Omega_j` has no negative
powers, descending from the coefficient of `q^{-jp}` shows

\[
C_m\equiv0\pmod{p^3}\qquad(m>0).
\]

At `q^0`, Section 4 gives

\[
[q^0]U_p\Omega_j=[q^0]\Omega_j=b_j,
\]

while `[q^0]Omega_0=1`; hence `C_0=b_j mod p^3`.  The exact grid identity then
upgrades the principal-part calculation to the full Laurent series. `square`

### 5.2 The Xi Cartier theorem

By linearity and `m-1<r<p`,

\[
U_p\mathcal X_r
\equiv
\left(1+5\sum_{m=1}^r g_m b_{m-1}\right)\Omega_0
\pmod {p^3}.
\]

Thus:

**[THEOREM]** For every prime `p>=5` and `p>r`,

\[
\boxed{U_p\mathcal X_r\equiv-\Xi_r\Omega_0\pmod{p^3}.}
\]

Consequently:

**[THEOREM — exact transverse reformulation]** If `p>r` and
`p|b_r, p|Xi_r`, then

\[
\boxed{U_p\Omega_r\equiv0\pmod p,\qquad
U_p\mathcal X_r\equiv0\pmod p.}
\]

No Hasse-polynomial evaluation is used here.  The row `r` is a Laurent-grid
index and the congruence is a Cartier/Hecke coefficient statement.

---

## 6. Where a linear support theorem would have to come from

### 6.1 Sturm bound

For level six,

\[
[SL_2(Z):Gamma_0(6)]=12.
\]

**[THEOREM]** If

\[
f=\sum_{n\ge0}a_nq^n\in M_k(\Gamma_0(6);\mathbf Z_{(p)})
\]

has `a_0=...=a_{B}=0 mod p` through

\[
B=\left\lfloor\frac{k}{12}[SL_2(Z):Gamma_0(6)]\right\rfloor=k,
\]

then `f=0 mod p`.  Equivalently a nonzero reduction has

\[
\operatorname{ord}_{\infty}(f\bmod p)\le k.
\]

This is the only place in a bare Sturm argument where a numerical support
constant can appear.

### 6.2 Why the number five is natural — conditionally

Multiplication by `Delta_6^r` clears the order-`r` poles of `Omega_r`, giving
weight

\[
4+4r.
\]

Therefore:

**[CONDITIONAL THEOREM — exact `C=5` gateway]** Suppose that for every common
pair `(p,r)` with `p>r` one can construct from the actual Apéry/Fricke/Eichler
data a form

\[
K_{p,r}\in M_{4r+4}(\Gamma_0(6);\mathbf Z_{(p)})
\]

such that

1. `K_{p,r} mod p` is nonzero; and
2. `ord_infinity(K_{p,r} mod p) >= p-r`.

Then

\[
p-r\le4r+4,
\]

so

\[
\boxed{p\le5r+4.}
\]

If the construction gave `ord >= p-r+1`, the endpoint improves to
`p<=5r+3`.  The `+4` is not cosmetic: it is exactly the level-six Sturm bound
for weight `4(r+1)`.

This conditional theorem identifies precisely how `p` and `r` must enter:

- `r` enters the **weight** through the cost of clearing an order-`r` pole;
- `p` must enter the **vanishing order**, through a new Cartier/Eichler
  pole-drop theorem.

The existing results provide the first bullet, not the second.

### 6.3 A smaller clearing weight does not help the actual condition

The weight-two eta quotient

\[
C=Et
\]

has cusp orders `(1,1,0,0)`.  It is therefore a more economical pole clearer:

\[
C^r\Omega_r=E^{r+1}H\in M_{2r+4}(\Gamma_0(6)).
\]

Its constant term is exactly `1`.

**[THEOREM — obstruction to naive Sturm]** Neither `p|b_r` nor the transverse
common condition can imply any positive order of vanishing at infinity for
this naturally cleared form, because

\[
[q^0](C^r\Omega_r)=1
\]

in every characteristic.

Likewise

\[
[q^0](\Delta_6^r\Omega_r)=1.
\]

Thus the missing `K_{p,r}` in the conditional theorem cannot simply be the
pole-cleared `Omega_r`; it must be a genuinely new **p-dependent defect or
Eichler/Cartier descendant** whose nonzero reduction and long gap are both
proved.

---

## 7. Why the two Cartier zeros do not presently imply the missing gap

This is the central negative audit.

Write

\[
\Omega_r=\sum_{n\ge-r}a_nq^n,
\qquad
\mathcal X_r=\sum_{n\ge-(r-1)}x_nq^n.
\]

For `p>r`, there is no negative integer in either principal part divisible by
`p`.  Hence

\[
U_p\Omega_r=\sum_{n\ge0}a_{pn}q^n,
\qquad
U_p\mathcal X_r=\sum_{n\ge0}x_{pn}q^n.
\]

The common condition, through Section 5, forces all coefficients in these two
`p`-spaced subsequences to vanish modulo `p`.  It says nothing directly about

\[
a_1,a_2,\ldots,a_{p-r-1}
\]

or the corresponding consecutive block for `mathcal X_r`.

**[THEOREM — distinction of predicates]** `U_p f=0 mod p` is the condition

\[
a_{pn}=0\pmod p\quad\text{for every }n,
\]

whereas a Sturm gap of length `L` is

\[
a_0=a_1=\cdots=a_{L-1}=0\pmod p.
\]

For `L<p` the former tests only `a_0` inside the interval
`0,...,L-1`.  Therefore no formal implication from Cartier kernel to a
length-`p-r` Sturm gap exists at the level of `q`-series alone.

The modular structure does strengthen the Cartier statement — the whole
`U_p` image is zero, not only its constant term — but it still samples one
residue class modulo `p`.  Fricke oddness repeats the same phenomenon at the
paired cusp `0`; it does not turn the arithmetic progression into consecutive
orders at `infinity`.

### 7.1 Why dimension counting does not insert `p`

The spaces `G(M,L)` have bases indexed by the Laurent pole interval.  For
`p>r`, the operator `U_p` erases the principal part of every `Omega_j` with
`j<=r`.  The exact Laurent theorem then collapses the entire span
`<Omega_0,...,Omega_r>` modulo `p` to the single line generated by `Omega_0`,
with coefficient `b_j`.  Thus its Cartier kernel has large dimension; two
independent kernel vectors are not contradictory.

A raw Riemann-Roch or Fricke dimension count therefore sees `r`, but not the
size of the prime `p`.  The only occurrence of `p` is the spacing in `U_p`.
To obtain a support theorem one needs a new result which converts that spacing
into a geometric order of vanishing or a drop in a pole divisor.

### 7.2 The exact missing statement

The strongest useful formulation I can isolate is:

**[REMAINING GAP — Cartier/Eichler pole drop]** Construct a canonical,
nonzero mod-`p` holomorphic form `K_{p,r}` from the pair

\[
(\Omega_r,\mathcal X_r)
\]

and the integral level-six Gauss-Manin/Eichler normalization such that common
Cartier vanishing implies

\[
\operatorname{ord}_{\infty}K_{p,r}\ge p-r
\]

with weight at most `4r+4` (or any comparable linear weight).

Equivalently, prove a theorem that the simultaneous Cartier kernel has a
`p`-dependent **pole drop** in the actual two-coordinate Apéry submodule.

This is not an issue of choosing a better Eisenstein basis.  The four
Eisenstein companions are linearly related by small integral coefficients,
and the `240` coefficient above makes coordinate-union support arguments
invalid.  The needed conclusion is geometric: a long gap or pole drop for one
nonzero section.

---

## 8. Is an absolute linear support theorem likely true?

No exact counterexample is known in the stated finite range, so the following
is deliberately not called a theorem.

### 8.1 Actual-Apéry counterexample mechanism

For each row define the characteristic-zero integer

\[
G_r:=\gcd(b_r,\Xi_r).
\]

**[THEOREM]** A prime `p>r` is a transverse common prime at row `r` if and only
if `p` is a prime divisor of `G_r` with `p>r`.

Therefore a counterexample to the constant `C` is not a toy model: it is an
exact prime factor

\[
p\mid G_r,\qquad p>Cr+A.
\]

There is no a priori reason from the recurrence, the eta quotient, or Sturm
alone for every prime factor of the integer `G_r` to be `O(r)`.  Such a result
would be a strong smoothness theorem for this specific Apéry/Eichler gcd
sequence.

The exact script

```text
problems/3.2/research/scripts/q7693_transverse_support.sage
```

computes `b_r`, `g_r`, `Xi_r`, factors `G_r` with Sage proof mode enabled, and
reports all common pairs, violations of `p<=Cr+A`, and repeated prime labels.
The command

```bash
sage problems/3.2/research/scripts/q7693_transverse_support.sage \
  --max-r 10000 --expect-known-10000
```

is a regression against the supplied finite ledger.  Larger `--max-r` values
are a direct exact search for a real Apéry counterexample.  A returned
violation is self-certifying: the script proves primality of `p` and checks
`b_r mod p = Xi_r mod p = 0`.

### 8.2 A heuristic reason not to assume any fixed `C`

**[HEURISTIC — not a theorem]** If the two residues `b_r mod p` and
`Xi_r mod p` behaved like independent nondegenerate coordinates as `r` ranges
through a positive proportion of `0,...,p-1`, then the chance of a common zero
at one row would be about `p^{-2}`.  Among the approximately `p/C` rows
`r<p/C`, the expected number of violations of `p<=Cr` at a fixed prime would
be about

\[
\frac1{Cp}.
\]

Since

\[
\sum_p\frac1p
\]

diverges, this naive model predicts infinitely many violations of every fixed
`C`, albeit extremely sparsely.  Conditional on a rare common hit, a positive
fraction would lie in `r<p/C`.

This is not evidence that the Apéry pair is independent; the entire point of
the modular problem is that it may have hidden rigidity.  It does show that an
absolute support theorem is **not** an automatic consequence of generic
Deligne/Katz equidistribution.  If `C=5` is true, it should come from a
specific deterministic exclusion such as the missing Cartier/Eichler pole-drop
theorem, not from probability or a finite scan.

### 8.3 A second exact search frontier

If full factorization of `G_r` becomes expensive, one can run a rectangular
finite-field search instead:

1. choose a prime bound `P`;
2. for every prime `p<=P`, compute the exact reductions of the already-defined
   integer sequences `b_r,Xi_r` for `0<=r<p/C`;
3. report rows where both vanish.

This certifies the absence of `p>Cr` only for `p<=P`, whereas factoring `G_r`
for `r<=R` certifies **all** prime sizes for those rows.  The two scans are
therefore complementary.

---

## 9. Consequences of any weaker absolute constant `C`

Let `A` be a fixed additive constant and assume the support statement

\[
p\mid b_r,\quad p\mid\Xi_r,\quad p>r
\quad\Longrightarrow\quad
p\le Cr+A.
\tag{LS(C,A)}
\]

For a dyadic row block `R<r<=2R`, define the label multiplicity

\[
\nu_p(R)=\#\{r:\ R<r\le2R,\ p>r,\ p\mid b_r,\ p\mid\Xi_r\}.
\]

It is essential not to replace `nu_p` by the indicator `nu_p>0` when the
application counts repeated labels.

### 9.1 Union/radical mass

Define the deduplicated prime mass

\[
M_{\rm union}(R)=\sum_{p:\nu_p(R)>0}\log p.
\]

**[THEOREM — conditional on `LS(C,A)`]** Every such prime satisfies

\[
R<p\le2CR+A,
\]

so Chebyshev's estimate gives

\[
M_{\rm union}(R)
\le\vartheta(2CR+A)
=O_C(R).
\]

Thus any fixed linear support constant, not specifically `5`, collapses the
unbounded prime tail to a bounded comparable-prime interval.

### 9.2 Repeated-label/content mass

Define instead

\[
M_{\rm lab}(R)=\sum_p\nu_p(R)\log p
=\sum_{\substack{R<r\le2R\\p>r\\p\mid b_r,\Xi_r}}\log p.
\]

**[THEOREM — exact warning]** `LS(C,A)` alone does **not** imply
`M_lab(R)=O(R)` or `o(R)`, because one prime may label many rows.  It gives only

\[
M_{\rm lab}(R)
\le \log(2CR+A)\sum_{p\le2CR+A}\nu_p(R).
\]

So repeated labels require a genuine incidence/multiplicity bound.  The finite
observation that no high prime repeats through `r=10000` cannot be promoted to
such a theorem.

There is a weak universal separation: two consecutive Apéry zeros modulo a
prime `p` with indices below `p` are impossible, because the second-order
recurrence has invertible transfer determinant there.  But this only gives a
multiplicity of order `p`, far too large to replace a real incidence estimate.

### 9.3 What a fixed `C` buys in a dyadic incidence argument

**[CONDITIONAL]** Suppose a separate theorem available for the common-content
problem gives, uniformly for each prime dyadic block `Y<p<=2Y` with `Y` a
constant multiple of `R`, an incidence-weight estimate

\[
\sum_{Y<p\le2Y}\nu_p(R)\log p
\ll R^{5/3+o(1)}.
\]

Then `LS(C,A)` partitions the entire high-prime contribution into only

\[
O(1+\log C)
\]

such prime blocks.  Therefore

\[
M_{\rm lab}(R)\ll_C R^{5/3+o(1)}.
\]

This implication **does count repeated labels**, because `nu_p(R)` remains in
every dyadic sum.  Hence even a poor absolute constant, say `C=20` or `100`,
would be structurally valuable: it would turn the open unbounded tail into
finitely many comparable-scale bands.  What it would not do by itself is
control multiplicity inside those bands.

---

## 10. Endpoint and normalization audit

1. **[THEOREM] `r=0`.** No common prime exists because `b_0=1`.
2. **[THEOREM] `p=2,3`.** Under `p>r`, the possible rows lie among
   `r=0,1,2`; `b_0=1,b_1=5,b_2=73` exclude them.
3. **[THEOREM] `p=5`.** `Xi_r=-1 mod 5` for every `r`, so there is no common
   pair.
4. **[THEOREM] square-root sign.** `sqrt(D)=1-17t-144t^2+...`; the constant
   term `+1` is forced by `H/E=1+O(q)` at infinity.
5. **[THEOREM] Fricke sign.** The factor `(-i tau)^2=-tau^2` in the eta
   transformation is essential: `E` is Fricke **odd**, not even.  `H` is even,
   hence `Psi` and every `Omega_m` are odd.
6. **[THEOREM] cusp form parity.** `Delta_6` is Fricke even; it cannot be
   inserted as an extra odd Laurent-grid basis vector.
7. **[THEOREM] coefficient versus value.** `b_r` and `Xi_r` are constant
   Fourier coefficients of `Omega_r` and `mathcal X_r`.  No step treats `r`
   as a geometric fiber parameter or evaluates a Hasse polynomial at `r`.
8. **[THEOREM] Sturm endpoint.** Nonzero weight `4r+4` reduction can have
   order exactly `4r+4`; therefore the conditional support conclusion is
   `p<=5r+4`, not automatically `p<5r+4`.
9. **[FINITE EVIDENCE]** `(17,13)` and `(2237,492)` are the supplied exact
   characteristic-zero hits through `10000`; their existence is not used in
   any proof above except as a sanity check against overstrong auxiliary
   claims.

---

## 11. The theorem DAG

```text
eta quotient t, E
      |
      +--> cusp divisors, widths, genus zero
      +--> W6 signs: t even, E odd
      |
H = q t'/t
      +--> H even, cusp order 0
      +--> H^2 = D(t) E^2               [Sturm weight 4]
      |
Psi = E H = E^2 sqrt(D)
      |
      +--> g(t)=1/Psi
      +--> Omega_m = Psi t^{-m}
      |       |
      |       +--> CT(Omega_m)=b_m      [residue theorem]
      |       +--> exact odd Laurent grid G(M,L)
      |       +--> U_p Omega_j = b_j Omega_0 mod p^3, j<p
      |
      +--> X_r = Omega_0 + 5 sum g_m Omega_{m-1}
              |
              +--> CT(X_r)=-Xi_r
              +--> U_p X_r = -Xi_r Omega_0 mod p^3, p>r

p | b_r, Xi_r
      |
      +--> U_p Omega_r = U_p X_r = 0 mod p
      |
      +--> [MISSING] Cartier/Eichler pole drop or consecutive q-gap
                    |
                    +--> nonzero K_{p,r} in M_{4r+4}
                         with ord_infty >= p-r
                    |
                    +--> Sturm: p-r <= 4r+4
                    |
                    +--> p <= 5r+4
```

The dashed/missing arrow is the only genuinely new theorem needed by this
Sturm architecture.  Everything before it is exact level-six structure; the
last implication is an elementary Sturm calculation.

---

## 12. Final assessment

**[THEOREM]** The level-six eta/Fricke/Laurent package puts the two common
coordinates into one exact modular framework and upgrades them to two Cartier
kernel relations.

**[THEOREM]** A nonzero weight-`4r+4` form with a common-condition-induced
`p-r` gap would prove `p<=5r+4`.

**[THEOREM]** The obvious pole-cleared forms do not have such a gap; their
leading constant is `1`.  `U_p` vanishing is sparse `p`-dissection vanishing,
not consecutive Sturm vanishing.

**[REMAINING GAP]** No currently established level-six Fricke/Sturm identity
converts the simultaneous Cartier kernel of `(Omega_r, X_r)` into the required
`p`-dependent pole drop.  This is where an actual Eichler/Gauss-Manin theorem,
not a finite scan or coordinate union, must enter.

**[FINITE EVIDENCE]** No violation of `p<=5r` is known through `r=10000`, and
the two known hits are nonrepeating in that box.

**[HEURISTIC]** Generic two-coordinate behavior would predict extremely sparse
but eventually unbounded ratios `p/r`; therefore an absolute `C`, if true,
is a genuinely rigid Apéry phenomenon.

Accordingly, the correct publication status is: **absolute linear support is
open; `C=5` is not presently proved; the exact missing lemma is a
Cartier/Eichler-to-Sturm pole-drop theorem.**

## References for standard ingredients

- F. Beukers, *Irrationality proofs using modular forms*, Asterisque 147-148
  (1987), 271-283.
- J. Stienstra and F. Beukers, *On the Picard-Fuchs equation and the formal
  Brauer group of certain elliptic K3-surfaces*, Math. Ann. 271 (1985),
  269-304.
- J. Sturm, *On the congruence of modular forms*, in Number Theory,
  Lecture Notes in Math. 1240 (1987), 275-280.
- Standard eta-quotient cusp-order criterion of Ligozat/Newman, used above
  only in the explicitly displayed level-six calculation.
