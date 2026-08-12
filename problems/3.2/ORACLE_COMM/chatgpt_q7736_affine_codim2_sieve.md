# Q7736 — defining-characteristic affine codimension-two sieve

## Verdict

**[THEOREM — exact residual affine transfer]** For every prime `p >= 7` and every `1 <= r < p`, the canonical Apéry period and the canonically normalized Eichler coordinate are two literal matrix coefficients of one affine transfer element

\[
\mathcal A_{p,r}=\begin{pmatrix}T_r&v_r\\0&1\end{pmatrix}\in \operatorname{AGL}_2(\mathbb F_p),
\qquad \det T_r=r^{-3}.
\]

With

\[
\widehat B=(5,1,0)^t,\qquad
\widehat K=(-36,0,1)^t,\qquad
\widehat e_1=(1,0,0)^t,
\]

one has exactly

\[
\boxed{
 b_r=\widehat e_1^t\mathcal A_{p,r}\widehat B,
 \qquad
 \kappa_r=\widehat e_1^t\mathcal A_{p,r}\widehat K
 \quad\text{in }\mathbb F_p.}
\tag{A}
\]

This is an unconditional residual **affine transfer** statement.  I do not call it a geometric affine-monodromy theorem: the repository has not yet constructed an integral lisse/overconvergent extension sheaf whose Frobenius matrix coefficient is the `kappa` coordinate.

**[THEOREM — exact ambient codimension two]** On every determinant coset

\[
\mathscr A_d(\mathbb F_p)
=\left\{
\begin{pmatrix}T&v\\0&1\end{pmatrix}:
T\in\operatorname{GL}_2(\mathbb F_p),\ \det T=d
\right\},\qquad d\ne0,
\]

the two functions appearing in (A) cut out an everywhere-smooth codimension-two complete intersection.  More precisely,

\[
\boxed{
\#\{A\in\mathscr A_d(\mathbb F_p):F(A)=H(A)=0\}=p^2(p-1),}
\tag{B}
\]

where

\[
F(A)=\widehat e_1^tA\widehat B,
\qquad H(A)=\widehat e_1^tA\widehat K.
\]

Since

\[
\#\mathscr A_d(\mathbb F_p)=p^3(p^2-1),
\]

the exact ambient density is

\[
\boxed{\frac1{p(p+1)}.}
\tag{C}
\]

Thus the two residual conditions really are transverse: the obstruction is not a hidden algebraic dependence between `b` and `kappa`.  What is missing is a theorem saying that the **actual nonautonomous samples** `A_{p,r}` see this codimension-two locus with anything like the ambient frequency.

**[THEOREM — bounded character order sector]** If

\[
d_{p,r}=\operatorname{ord}(\omega_p^r)
=\frac{p-1}{\gcd(p-1,r)}\le D,
\]

then for fixed `r` there are at most `D tau(r)` possible primes `p`: indeed

\[
p=d g+1,\qquad d\le D,\quad g\mid r.
\]

Consequently, for `R < r <= 2R`, the logarithmically weighted transverse mass in this sector satisfies

\[
\boxed{
W_D(R)\ll D R\log R\,\log(2DR+1).}
\tag{D}
\]

In particular every fixed-order sector is `O_D(R log^2 R)=o(R^2)`, and more generally the same holds whenever `D log R log(2DR)=o(R)`.  This is unconditional and requires no monodromy theorem.

**[AUXILIARY-ELL ONLY]** The level-6/Symmetric-square geometry makes large *linear* residual image at auxiliary primes `ell != p` plausible and standard after the exact modular-cover monodromy map is packaged.  A full affine image would additionally require proving that the Eichler extension class remains nonzero modulo `ell`; irreducibility would then force the translation subgroup to be all of the residual Sym^2 module.  Neither this affine extension-sheaf realization nor its mod-`ell` nonvanishing is currently banked in the repository.  More importantly, even a perfect auxiliary-`ell` theorem would not imply the target event at `ell=p`.

**[NO-GO — diagonal step]** Existing compatible-system Frobenius sieves and current bilinear trace-function estimates work with an auxiliary `ell` invertible in the base field.  Here the residual condition is taken at the **defining characteristic itself**.  There is no legitimate specialization `ell=p` of those theorems.  I therefore obtain no unconditional power-saving bilinear dispersion estimate in `(p,r)` from the standard trace-function literature.

**[CONDITIONAL — smallest new estimate]** Let

\[
N_R(P)=
\sum_{\substack{R<r\le2R\\P<p\le2P\\p>2R\\p\mid b_r,\ p\mid K_r}}
\log p,
\tag{E}
\]

where `K_r` is the primitive numerator of `kappa_r`.  The ambient density (C) predicts the main term `~ R/P`.  A diagonal affine-dispersion theorem of the shape

\[
\boxed{
N_R(P)\ll \frac RP+R^{1-\delta}P^A
}
\tag{F}
\]

would therefore be the correct codimension-two target **locally in a defining-prime dyadic block**.

There is an important summation qualification: if `P` in (F) denotes the defining-prime scale and `A>0`, (F) alone does **not** imply `o(R^2)` over arbitrarily large prime factors, because the error grows with `P`.  A closing theorem must either have `A=0` (or a summable negative `P`-dependence), or prove (F) only up to `P<=R^B` with `AB<1+delta` together with a separate `o(R^2)` tail above `R^B`.  Under that qualified form, summing dyadic blocks gives `o(R^2)`.  Calling the desired diagonal estimate a “large sieve” would hide precisely the new step.

The exact verifier is

`problems/3.2/research/scripts/q7736_affine_transversality.sage`.

A Sage 10.6 GitHub Actions run (`31583508887`) passed all matrix-coefficient identities, the locked common pairs `(17,13)` and `(2237,492)`, the nontransverse control `(11,5)`, and exhaustive ambient count/Jacobian checks for `p=7,11,13`.

---

## 1. Canonical normalization: why `kappa` is a legitimate second coordinate

Write

\[
P(n)=34n^3+51n^2+27n+5.
\]

The repository has already fixed the homogeneous companion `u` by

\[
u_0=0,\qquad u_1=1,
\]

with the Apéry recurrence, and proved the unit Casoratian

\[
\boxed{r^3(b_{r-1}u_r-b_ru_{r-1})=1.}
\tag{1.1}
\]

It also fixes the inhomogeneous coordinate by

\[
\kappa_0=0,\qquad \kappa_1=-36,
\]

and, for `r>=2`,

\[
\boxed{
r^3\kappa_r-P(r-1)\kappa_{r-1}+(r-1)^3\kappa_{r-2}=-5g_r.}
\tag{1.2}
\]

The first row has the real boundary defect

\[
A_1\kappa=-36=-5g_1-1,
\]

so it must be incorporated in the initial vector, not erased.

The scale-sensitive cross-row identity is

\[
\boxed{
\Xi_r=r^3(b_{r-1}\kappa_r-b_r\kappa_{r-1}).}
\tag{1.3}
\]

If `kappa_r=K_r/d_r` is in lowest terms, the repository proves that every prime divisor of `d_r` is at most `r`.  Therefore, for every prime `p>r` with `p|b_r`, (1.1) makes `r^3b_{r-1}` a `p`-unit and (1.3) gives

\[
\boxed{
p\mid\Xi_r\iff p\mid K_r\iff \kappa_r=0\text{ in }\mathbb F_p.}
\tag{1.4}
\]

Thus the requested pair `(b_r,kappa_r)` is not a rescaled surrogate.  On the high-prime target it is exactly the canonical transverse pair.

---

## 2. The exact residual affine transfer object

For `p>r`, all denominators below are invertible in `F_p`.  Put

\[
M_j=
\begin{pmatrix}
P(j)/(j+1)^3&-j^3/(j+1)^3\\
1&0
\end{pmatrix}
\quad(1\le j<p-1).
\tag{2.1}
\]

For a homogeneous solution `y`, the state vector obeys

\[
\binom{y_{j+1}}{y_j}=M_j\binom{y_j}{y_{j-1}}.
\]

Define

\[
T_1=I_2,
\qquad
T_r=M_{r-1}\cdots M_1\quad(r\ge2).
\tag{2.2}
\]

Since

\[
\det M_j=\frac{j^3}{(j+1)^3},
\]

the determinant telescopes:

\[
\boxed{\det T_r=r^{-3}.}
\tag{2.3}
\]

Let

\[
B=(5,1)^t,
\qquad K=(-36,0)^t.
\]

Then

\[
\binom{b_r}{b_{r-1}}=T_rB.
\tag{2.4}
\]

For the inhomogeneous coordinate define a translation vector by

\[
v_1=0,
\]

and

\[
\boxed{
v_{j+1}=M_jv_j-rac{5g_{j+1}}{(j+1)^3}e_1.}
\tag{2.5}
\]

Induction using (1.2) gives

\[
\boxed{
\binom{\kappa_r}{\kappa_{r-1}}=T_rK+v_r.}
\tag{2.6}
\]

Now form the affine matrix

\[
\mathcal A_{p,r}
=\begin{pmatrix}
T_r&v_r\\
0&1
\end{pmatrix}.
\tag{2.7}
\]

With the augmented vectors from the verdict,

\[
\widehat B=(5,1,0)^t,
\qquad
\widehat K=(-36,0,1)^t,
\]

(2.4) and (2.6) become the two exact matrix coefficients (A).

### What “monodromy” means here

The object (2.7) is the strongest unconditional object I can construct directly from the banked arithmetic.  It is a nonautonomous affine cocycle in the row variable `r`.  A true geometric affine-monodromy representation would require an extension object over the level-6 parameter curve, schematically

\[
0\longrightarrow \mathcal V_\ell
\longrightarrow \mathcal E_\ell
\longrightarrow \mathbf Q_\ell
\longrightarrow0,
\tag{2.8}
\]

where `V_ell` is the primitive rank-three K3/Sym^2 local system and the extension class realizes the Eichler coordinate.  Its monodromy would land in

\[
V_\ell\rtimes G_\ell,
\]

with `G_ell` of orthogonal/Sym^2 type.

The repository currently proves the K3/Sym^2 geometry and the scalar Eichler recurrence, but it does not yet prove that (2.8) exists integrally with `kappa` as the desired Frobenius/Mellin coefficient.  I therefore do not use (2.8) as an input.

---

## 3. Unconditional codimension-two transversality

Fix `d in F_p^*` and write

\[
T=\begin{pmatrix}a&b\\c&e\end{pmatrix},
\qquad v=(x,y)^t.
\]

On `A_d` the two functions are

\[
F=5a+b,
\qquad
H=-36a+x.
\tag{3.1}
\]

The common-zero equations are therefore

\[
b=-5a,
\qquad x=36a,
\qquad ae-bc=d.
\tag{3.2}
\]

Substituting `b=-5a` into the determinant equation gives

\[
a(e+5c)=d.
\tag{3.3}
\]

Since `d != 0`, every common point has `a != 0`.  Conversely:

- choose `a in F_p^*`: `p-1` choices;
- choose `c in F_p`: `p` choices;
- equation (3.3) determines `e` uniquely;
- `b=-5a` and `x=36a` are forced;
- choose `y in F_p`: `p` choices.

This proves (B).

For smoothness regard the locus in the six affine coordinates `(a,b,c,e,x,y)`.  The three equations are

\[
ae-bc-d=0,
\qquad 5a+b=0,
\qquad -36a+x=0.
\]

Their Jacobian is

\[
\begin{pmatrix}
e&-c&-b&a&0&0\\
5&1&0&0&0&0\\
-36&0&0&0&1&0
\end{pmatrix}.
\tag{3.4}
\]

The third row is independent of the first two because it is the only one with an `x` component.  The determinant row cannot be a multiple of the second row at a common point because its `(c,e)`-coordinate pair `(-b,a)` cannot vanish when `a != 0`.  Thus the rank is exactly three everywhere.  Since the determinant coset itself has dimension five, the two additional equations cut codimension two smoothly.

Finally

\[
\#\{T:\det T=d\}=\#SL_2(\mathbb F_p)=p(p^2-1),
\]

and there are `p^2` translations, proving (C).

A useful refinement is

\[
\Pr_{A\in\mathscr A_d}[F(A)=0]=\frac1{p+1},
\]

and, conditional on `F(A)=0`,

\[
\Pr[H(A)=0\mid F(A)=0]=\frac1p.
\]

Thus the second equation supplies a literal residual factor `1/p` after the first.  The missing theorem is to transfer any part of that ambient conditional density to the actual affine cocycle samples.

---

## 4. The dyadic codimension-two target and the summation audit

For the block `(R,2R]` and defining-prime dyadic scale `(P,2P]`, take `P>2R` and define `N_R(P)` by (E).  Because `p>2R>=r`, (1.4) makes this exactly the original transverse common-prime event.

If the actual affine samples behaved with the ambient density (C), then for a fixed `r` the logarithmically weighted expected contribution of primes `p~P` would be

\[
\sum_{P<p\le2P}\frac{\log p}{p^2}\asymp\frac1P,
\]

and over `R` rows the predicted size is `R/P`.  This is the origin of the main term in (F); it is not guessed from an unrelated large-sieve formula.

### Conditional local estimate

A concrete new theorem to aim for is:

> **Diagonal affine dispersion, local form.** There exist fixed `delta>0` and `A>=0` such that, uniformly for all dyadic `R>=R_0` and all defining-prime blocks in the stated range,
> \[
> N_R(P)\ll R/P+R^{1-\delta}P^A.
> \]

This is a statement about the two-dimensional arithmetic array `(p,r)` and the actual cocycle (2.7).  It is not a consequence of the ambient count and it is not being renamed “large sieve”.

### The literal shape needs a tail qualification

If `A>0`, the local statement by itself does not close the all-prime height sum: the error grows as the defining prime grows.  A correct closing formulation is, for example:

> There exist `delta>0`, `A>=0`, `B>0` with `AB<1+delta` such that (F) holds for every `2R<P<=R^B`, and independently
> \[
> \sum_{\substack{R<r\le2R\\p>R^B\\p\mid b_r,K_r}}\log p=o(R^2).
> \tag{4.1}
> \]

Then the dyadic main terms sum geometrically, while the error terms over `P<=R^B` are

\[
\ll R^{1-\delta+AB}\log R=o(R^2),
\]

and (4.1) finishes the tail.

If one could obtain `A=0`, then even the exponentially many possible defining-prime scales allowed by the elementary row heights can be summed with a power saving: `O(R)` dyadic scales times `R^{1-delta}` is `O(R^{2-delta})`.

This distinction matters.  If `P` were instead an independently chosen *auxiliary sieve level*, one could optimize `R/P+R^{1-delta}P^A`; but that is not the same variable as the defining-prime dyadic scale in (E).

---

## 5. What can actually be proved now

### 5.1 Residual primitivity/transversality: yes, in the exact ambient sense

Theorem (B) is stronger than merely saying that the two linear forms are not proportional.  The two equations remain independent after the determinant constraint, the common zero is everywhere smooth, and it has the exact expected `p^{-2}` density.  Moreover at an actual Apéry zero the first row of `T_r` cannot degenerate: if

\[
5a+b=0,
\]

then `a=0` would force `b=0`, contradicting invertibility.  At an actual common zero one therefore has

\[
\boxed{a\ne0,\qquad b=-5a,\qquad x=36a.}
\tag{5.1}
\]

This is the unconditional residual primitivity statement I am willing to bank.

What it does not prove is that the finite set

\[
\{\mathcal A_{p,r}:R<r\le2R\}
\]

has expansion, equidistribution, or even positive-dimensional Zariski closure inside the ambient coset for one `p`, let alone uniformly as `p` moves.

### 5.2 Bounded-character-order sector: yes, unconditionally

For a Kummer character `omega_p^r`, put

\[
d_{p,r}=\frac{p-1}{\gcd(p-1,r)}.
\]

If `d_{p,r}<=D`, let `g=gcd(p-1,r)`.  Then

\[
g\mid r,
\qquad p=d_{p,r}g+1.
\]

For fixed `r`, each pair `(d,g)` with `d<=D` and `g|r` determines at most one integer `p`, hence at most one prime.  Therefore

\[
\#\{p>r:d_{p,r}\le D\}\le D\tau(r).
\tag{5.2}
\]

If `R<r<=2R`, every such candidate satisfies

\[
p\le Dr+1\le2DR+1.
\]

Hence, even before imposing `p|b_r,K_r`,

\[
\begin{aligned}
W_D(R)
&\le \log(2DR+1)\sum_{R<r\le2R}D\tau(r)\\
&\ll DR\log R\,\log(2DR+1),
\end{aligned}
\]

using the standard divisor-sum bound `sum_{n<=x} tau(n) = O(x log x)`.  This proves (D).

This sector is therefore not where the main difficulty lives.  Any putative bad mass of quadratic size must come from characters whose order grows with `R`.

### 5.3 Full affine image for an auxiliary family: not yet, but reduced to one precise extension-class problem

The repository identifies the primitive K3 system, after the level-6 quadratic cover, with a symmetric-square/Asai construction from the elliptic family

\[
E_x:\quad Y^2+(1-2x)XY+x^2Y=X^3,
\]

and exhibits an exact rational point of order six.  Thus the linear rank-three system is of the expected `Sym^2` / orthogonal type.

At an auxiliary prime `ell` prime to the base characteristic and to the level, the standard modular-cover route is:

1. prove the exact level-6 monodromy subgroup on `H^1(E_x)`;
2. reduce it modulo `ell`, obtaining the expected `SL_2(F_ell)` image on the two-dimensional system;
3. apply `Sym^2`, whose projective image is `PSL_2(F_ell)` and whose three-dimensional module is irreducible for `ell>=5`.

The **affine** upgrade needs more.  Suppose an extension sheaf (2.8) has been constructed and let

\[
H_\ell\subset V_\ell\rtimes G_\ell
\]

be its residual image.  If

- the projection of `H_ell` onto `G_ell` is surjective;
- `V_ell` is an irreducible `G_ell`-module; and
- `H_ell cap V_ell != 0`, equivalently the residual extension/cocycle contributes a nonzero translation,

then `H_ell cap V_ell` is a nonzero invariant subspace, hence equals `V_ell`, and therefore

\[
H_\ell=V_\ell\rtimes G_\ell.
\tag{5.3}
\]

So the smallest missing **auxiliary** lemma is not “prove a large sieve”.  It is:

> construct the integral Eichler extension local system and prove that its extension class is nonzero modulo all but finitely many auxiliary `ell`.

I do not have this theorem from the current repository state, and I do not use fixed-`S` results with an `S` that grows with `R` to manufacture it.

### 5.4 Bilinear dispersion in `(p,r)`: no nontrivial theorem obtained

Current bilinear trace-function theorems are the wrong quantifier for the diagonal event.  They start with a sheaf over one finite field `F_q` and an auxiliary `ell` invertible in `F_q`, then estimate complex trace-function bilinear forms inside that fixed field.  They can be relevant to a future fixed-`p` row-dispersion input once a suitable trace-function realization of the affine cocycle is known.  They do not estimate

\[
\mathbf 1_{\{b_r=0\bmod p\}}
\mathbf 1_{\{\kappa_r=0\bmod p\}}
\]

while the same rational prime `p` simultaneously changes the field, the residual coefficient prime, and the row recurrence.

I therefore obtain no unconditional `R^{1-delta}` bilinear saving in the two moving variables.

---

## 6. Auxiliary `ell` versus diagonal `ell=p`

This is the main logical firewall.

### Auxiliary-`ell` statement

A compatible-system/Frobenius sieve fixes the geometric family in characteristic `p_0` and reduces its `ell`-adic representations at coefficient primes `ell` with

\[
\ell\ne p_0.
\]

Large residual image can then make a codimension-two target have local density about `ell^{-2}` across many auxiliary reductions.  Hall's big orthogonal/symplectic monodromy criteria and Frobenius large-sieve frameworks of Kowalski/Zywina/Perret-Gentil belong to this category.

### Our diagonal statement

Here, for every pair `(p,r)`, the event itself is

\[
b_r\equiv\kappa_r\equiv0\pmod p.
\]

The coefficient prime is the characteristic of the finite field in which the residual recurrence is being evaluated.  Setting the auxiliary coefficient prime equal to the base characteristic is not an allowed specialization of an etale compatible-system theorem.

The same issue persists if one lets the base characteristic vary: the forbidden coefficient prime moves diagonally with it.  A crystalline companion at each `p` provides an object *at that p*; it does not supply independence between the diagonal objects for different rational primes.

This is why an auxiliary full-image theorem, even if completed, would be evidence and infrastructure but not the missing all-index theorem.

---

## 7. Exact finite-field verifier

The committed Sage script is

```text
problems/3.2/research/scripts/q7736_affine_transversality.sage
```

It independently performs the following checks.

1. Reconstructs `b_0,...,b_r` modulo `p` from the Apéry recurrence.
2. Reconstructs `Q=(1-34t+t^2)^(-1/2)`, then `g=Q/F^2` as a truncated exact power series over `F_p`.
3. Reconstructs `kappa` from the inhomogeneous recurrence with the genuine initial value `kappa_1=-36`.
4. Builds `T_r` and `v_r` from (2.2) and (2.5).
5. Checks, row by row,
   \[
   \det T_r=r^{-3},\quad
   b_r=e_1^tT_rB,\quad
   \kappa_r=e_1^t(T_rK+v_r).
   \]
6. At every tested `b_r=0`, checks `a!=0`; at common zeros checks `v_{r,1}=36a`.
7. Checks the locked rows `(11,5)`, `(17,13)`, `(2237,492)`.
8. Exhaustively enumerates one nonzero determinant coset for `p=7,11,13`, counts the first-zero and common-zero loci, and checks the Jacobian rank at every common point.
9. Checks the exact divisor parametrization behind the bounded-character-order sector on a finite sample.

The executed Sage 10.6 run printed

```text
LOCKED 11 5 b0 1 k0 0 det 3 a 3 v0 2
LOCKED 17 13 b0 1 k0 1 det 13 a 16 v0 15
LOCKED 2237 492 b0 1 k0 1 det 1703 a 722 v0 1385
AMBIENT 7 total 16464 first_zero 2058 common 294 ratio_den 56 smooth 294
AMBIENT 11 total 159720 first_zero 13310 common 1210 ratio_den 132 smooth 1210
AMBIENT 13 total 369096 first_zero 26364 common 2028 ratio_den 182 smooth 2028
BOUNDED_ORDER R 80 D 8 checked 290
Q7736_AFFINE_TRANSVERSALITY PASS
```

The ambient counts agree exactly with

\[
\#\mathscr A_d=p^3(p^2-1),
\quad
\#\{F=0\}=p^3(p-1),
\quad
\#\{F=H=0\}=p^2(p-1).
\]

The temporary workflow used only to execute this audit was removed after the successful run.

---

## 8. Literature contacts and exact scope

The following contacts are useful only with their real quantifiers.

- **Hain, _The Hodge-de Rham Theory of Modular Groups_, arXiv:1403.6443.** Relative completion of modular groups packages extensions of symmetric-power variations by modular/Eisenstein classes and normal functions.  This is the right conceptual home for an Eichler affine extension, but it is not a residual defining-characteristic zero-density theorem.
- **Hall, _Big symplectic or orthogonal monodromy modulo l_, arXiv:math/0608718.** Gives criteria for uniformly big residual orthogonal/symplectic monodromy in compatible systems.  This is auxiliary-`ell` infrastructure, not `ell=p`.
- **Perret-Gentil, _Exponential sums over finite fields and the large sieve_, arXiv:1703.06965.** Gives zero-density results for algebraic conditions on trace functions through compatible systems and auxiliary residual primes.  The construction keeps the coefficient primes in the etale range; it does not authorize reduction at the defining characteristic.
- **Zywina, _The Large Sieve and Galois Representations_, arXiv:0812.2222.** Nonabelian Frobenius sieve for systems of Galois representations, again auxiliary-prime technology.
- **Fouvry--Kowalski--Michel--Sawin, _Bilinear forms with trace functions_, arXiv:2511.09459.** Gives strong bilinear bounds for trace functions over a fixed finite field under monodromy hypotheses.  It does not supply the moving-field, same-prime residual indicator required here.
- **Fouvry--Kowalski--Michel, _Algebraic trace functions over the primes_, arXiv:1211.6043.** Studies prime arguments fed to a trace function modulo a fixed finite-field characteristic.  Here the modulus/field itself is the prime being tested, so this is not the same horizontal diagonal.

No fixed-`S` theorem is used with a set `S` depending on `R`.

---

## 9. Smallest new estimate that would move the theorem

The exact algebraic work has reduced the problem to one quantitative statement.

> **Diagonal affine anti-concentration.**  Prove that the actual nonautonomous affine cocycle samples
> \[
> \mathcal A_{p,r},\qquad R<r\le2R,
> \]
> do not concentrate on the smooth codimension-two locus
> \[
> F=H=0
> \]
> as the defining characteristic `p` varies, with a uniform power saving whose error is summable over the defining-prime scales.

A local model is (F), with the summation qualification in Section 4.  An even cleaner closing formulation would be a directly summable estimate such as

\[
N_R(P)
\ll \frac RP
+R^{1-\delta}\left(\frac RP\right)^\eta
\qquad(P>2R)
\tag{9.1}
\]

for some fixed `delta,eta>0`; its dyadic errors sum geometrically.  I do not prove (9.1).

The ambient theorem shows that such an estimate would exploit a real codimension-two structure rather than an artificial choice of coordinates.  The bounded-order theorem removes one natural exceptional sector.  What remains is precisely the growing-order, defining-characteristic dispersion problem.
