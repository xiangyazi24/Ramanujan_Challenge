# Q7694: scale-sensitive cross-row identity for the transverse Apéry content

## Verdict

There is an exact Apéry-specific cross-row identity that is stronger than the
projective pair minors and removes the arbitrary independent row scaling:

\[
\boxed{\Xi_r=r^3\bigl(b_{r-1}\kappa_r-b_r\kappa_{r-1}\bigr).}
\tag{1}
\]

It comes from the **unit-Casoratian normalized Apéry frame**, not from generic
linear algebra.  It has two useful consequences.

1. For every prime \(p>r\), on the actual transverse target \(p\mid b_r\), the
   factor \(r^3b_{r-1}\) is a \(p\)-adic unit.  Hence, after reducing
   \(\kappa_r=K_r/d_r\) to lowest terms,
   \[
   v_p\gcd(b_r,\Xi_r)=v_p\gcd(b_r,K_r).
   \tag{2}
   \]
   Thus the high-prime content is recovered **exactly** from a canonical
   scale-fixed one-row integer coordinate \(K_r\).
2. This does **not** solve the height problem.  The primitive numerator carrier
   \(K_r\) has quadratic raw block height, and exact computation shows that
   removing every prime factor \(q\le 2R\) barely changes that height through
   \(R=256\).  So (1) is a genuine scale-sensitive identity, but its natural
   rowwise integer carrier is not a hidden subquadratic certificate.

The rigorous no-go proved below is deliberately scoped: **raw / primitive
numerator rowwise height is \(\Omega(R^2)\)**.  After maximal small-prime
saturation, I do not have an asymptotic lower bound or upper bound; the exact
finite audit is strongly negative for this candidate, but is not promoted to a
theorem.  Therefore the surviving class is the adjacent unit-Casoratian / banded
source identities, to be used inside a genuinely multirow determinant/Fitting
construction (or a modular/Hecke identity) before one multiplies one exponential
integer per row.

No reflected-depth law is used.  The fixed six-slope passport is not used.
BFH is not used to infer anything about all primes dividing \(b_r\): all target
statements below remain conditioned on the simultaneous event
\(p\mid b_r,\Xi_r\), exactly as required.

An independent standard-library verifier is

`problems/3.2/research/scripts/q7694_scale_sensitive_crossrow_verify.py`.

A remote exact run at `--N 520` passed all identities, denominator checks, the
locked regressions `(17,13)`, `(2237,492)`, `(11,5)`, and all block target
comparisons.

---

## 1. Canonical homogeneous companion and the unit Casoratian

Write

\[
P(n)=34n^3+51n^2+27n+5.
\]

Besides the Apéry solution \(b\), define the canonical homogeneous companion
\(u\) by

\[
u_0=0,\qquad u_1=1,
\]
\[
(n+1)^3u_{n+1}=P(n)u_n-n^3u_{n-1}.
\tag{3}
\]

Set

\[
W_r=b_{r-1}u_r-b_ru_{r-1}.
\]

Applying the two homogeneous recurrences at \(n=r-1\) gives

\[
r^3W_r=(r-1)^3W_{r-1}.
\]

Since \(W_1=b_0u_1-b_1u_0=1\), induction yields the exact unit normalization

\[
\boxed{r^3\bigl(b_{r-1}u_r-b_ru_{r-1}\bigr)=1.}
\tag{4}
\]

This is the scale-sensitive input unavailable to the generic row-scaling
argument: \(u_0=0,u_1=1\) fixes the second solution, and (4) fixes the
Casoratian to the literal unit \(1\).

---

## 2. Variation-of-parameters frame and a synchronized adjacent row

Define

\[
\Phi_0=0,\qquad
\Phi_r=5\sum_{m=1}^r g_m u_{m-1},
\tag{5}
\]

and use the exact frame

\[
\kappa_r=\Xi_r u_r+\Phi_r b_r.
\tag{6}
\]

The prefix definitions give

\[
\Xi_r-\Xi_{r-1}=-5g_rb_{r-1},\qquad
\Phi_r-\Phi_{r-1}=+5g_ru_{r-1}.
\tag{7}
\]

The two increments cancel when evaluated on the previous homogeneous row:

\[
\begin{aligned}
\Xi_ru_{r-1}+\Phi_rb_{r-1}
&=(\Xi_{r-1}-5g_rb_{r-1})u_{r-1}
 +(\Phi_{r-1}+5g_ru_{r-1})b_{r-1}\\
&=\kappa_{r-1}.
\end{aligned}
\tag{8}
\]

Hence the frame is not merely a one-row decomposition; it gives the exact
adjacent system

\[
\binom{\kappa_r}{\kappa_{r-1}}
=
\begin{pmatrix}
u_r&b_r\\u_{r-1}&b_{r-1}\end{pmatrix}
\binom{\Xi_r}{\Phi_r}.
\tag{9}
\]

By (4), the determinant of the matrix in (9) is exactly \(1/r^3\).  Inverting
it gives **both** cross-row identities

\[
\boxed{\Xi_r=r^3(b_{r-1}\kappa_r-b_r\kappa_{r-1}),}
\tag{10}
\]

\[
\boxed{\Phi_r=r^3(u_r\kappa_{r-1}-u_{r-1}\kappa_r).}
\tag{11}
\]

Equation (10) is the requested new scale-sensitive transverse identity.
Unlike \(b_s\Xi_r-b_r\Xi_s\), it is an adjacent determinant in a **canonically
normalized unit-Casoratian frame**; no independent rescaling of row \(r\) is
available.

### Small exact check

The reconstructed sequences begin

\[
(b_0,\ldots,b_5)=(1,5,73,1445,33001,819005),
\]

\[
(g_0,\ldots,g_4)=(1,7,192,5520,165168),
\]

\[
(\Xi_0,\Xi_1,\Xi_2,\Xi_3)=(-1,-36,-4836,-2019636),
\]

and

\[
\kappa_0=0,\quad \kappa_1=-36,\quad
\kappa_2=-\frac{1293}{2},\quad
\kappa_3=-\frac{82931}{6}.
\]

Then

\[
8\left(5\left(-\frac{1293}{2}\right)-73(-36)\right)=-4836=\Xi_2,
\]

and

\[
27\left(73\left(-\frac{82931}{6}\right)
 -1445\left(-\frac{1293}{2}\right)\right)
=-2019636=\Xi_3.
\]

The verifier checks (4), (6), (8), (10), and (11) exactly for every row through
80, with \(\kappa\) reconstructed independently from the source recurrence.

---

## 3. The source recurrence, including the boundary defect

Subtract (10) at consecutive rows.  For \(r\ge2\),

\[
\begin{aligned}
\Xi_r-\Xi_{r-1}
={}&r^3(b_{r-1}\kappa_r-b_r\kappa_{r-1})\\
&-(r-1)^3(b_{r-2}\kappa_{r-1}-b_{r-1}\kappa_{r-2}).
\end{aligned}
\]

Using

\[
r^3b_r+(r-1)^3b_{r-2}=P(r-1)b_{r-1},
\]

this factors as

\[
\Xi_r-\Xi_{r-1}
=b_{r-1}\Bigl(
 r^3\kappa_r-P(r-1)\kappa_{r-1}+(r-1)^3\kappa_{r-2}
\Bigr).
\]

But the left side is \(-5g_rb_{r-1}\), and \(b_{r-1}>0\).  Therefore

\[
\boxed{
A_r\kappa:={r^3\kappa_r-P(r-1)\kappa_{r-1}+(r-1)^3\kappa_{r-2}}
=-5g_r,\qquad r\ge2.}
\tag{12}
\]

There is a real normalization defect at the first row:

\[
A_1\kappa=\kappa_1-P(0)\kappa_0=-36,
\]

whereas \(-5g_1=-35\).  Thus a uniform statement is

\[
A_r\kappa=-5g_r-\mathbf 1_{r=1}.
\tag{13}
\]

The verifier explicitly rejects the incorrect boundary equation
\(A_1\kappa=-5g_1\).

This matters because a formal three-term source transport that starts at row 1
can otherwise hide a one-unit source term.

---

## 4. Exact summation-by-parts identity — useful, but not a second height gain

For \(2\le a\le s\), multiplying (12) by \(b_{m-1}\) and summing gives

\[
\boxed{
\sum_{m=a}^s b_{m-1}A_m\kappa
=-5\sum_{m=a}^s g_mb_{m-1}
=\Xi_s-\Xi_{a-1}.}
\tag{14}
\]

Substituting (10) at the two endpoints gives the boundary form

\[
\begin{aligned}
\Xi_s-\Xi_{a-1}
={}&s^3(b_{s-1}\kappa_s-b_s\kappa_{s-1})\\
&-(a-1)^3(b_{a-2}\kappa_{a-1}-b_{a-1}\kappa_{a-2}).
\end{aligned}
\tag{15}
\]

This is the exact discrete Green / summation-by-parts identity suggested by the
prefix definition.  It is checked by the verifier for many intervals.

However, (14)-(15) are algebraically equivalent to the adjacent identity (10)
plus the Apéry recurrence.  They **must not be counted as a separate height
improvement**.  Their value is structural: all interior source terms telescope
to two normalized boundary determinants.

---

## 5. Denominator control and exact high-prime target preservation

Write \(\kappa_r=K_r/d_r\) in lowest terms with \(d_r>0\).

From (12), starting with \(\kappa_0=0,\kappa_1=-36\), induction shows:

> Every prime divisor of \(d_r\) is at most \(r\).

Indeed the only new denominator introduced at row \(r\) is division by
\(r^3\); previous denominators have prime support at most \(r-1\).  The verifier
checks this assertion exactly through row 520 after reducing every fraction to
lowest terms.

Now let \(p>r\) and assume the **actual target condition** \(p\mid b_r\).
Reduce the unit Casoratian (4) in \(\mathbf Z_{(p)}\):

\[
r^3b_{r-1}u_r\equiv1\pmod p.
\]

Thus \(p\nmid r b_{r-1}\), and all relevant denominators are \(p\)-units.
Equation (10) becomes

\[
\Xi_r=(r^3b_{r-1})\kappa_r-b_r(r^3\kappa_{r-1}),
\tag{16}
\]

with the first coefficient a \(p\)-adic unit and the second term a multiple of
\(b_r\).  Therefore, including prime powers,

\[
\boxed{
\min(v_p(b_r),v_p(\Xi_r))
=
\min(v_p(b_r),v_p(\kappa_r))
=
\min(v_p(b_r),v_p(K_r)).}
\tag{17}
\]

For \(r\in(R,2R]\) and \(p>2R\), (17) applies automatically.  Hence

\[
\boxed{
[C_R]_{p>2R}
=
\prod_{R<r\le2R}
[\gcd(b_r,K_r)]_{p>2R}.}
\tag{18}
\]

This is an **exact target-preserving integer carrier**.  It is not a heuristic
or a slope passport.

Crucially, (17) was derived only after imposing \(p\mid b_r\).  Nothing here
claims control of all primes dividing \(b_r\).  BFH, if used later, may only be
applied to the simultaneous transverse event \(p\mid b_r,\Xi_r\) and in the
row-aggregate form actually available; no pointwise BFH consequence is used or
needed here.

### Locked regressions

The independent run verifies:

- \((p,r)=(17,13)\): \(p\mid b_r,\Xi_r,K_r\);
- \((p,r)=(2237,492)\): \(p\mid b_r,\Xi_r,K_r\);
- \((p,r)=(11,5)\): \(11\mid b_5\) but
  \(11\nmid\Xi_5\) and \(11\nmid K_5\).

Thus the conversion preserves both known transverse pairs and the locked
nontransverse control.

---

## 6. Saturation before height: the natural carrier fails

A raw numerator height comparison is not sufficient.  Since only primes
\(p>2R\) matter, define the maximally legitimate small-prime saturation

\[
\operatorname{sat}_{>2R}(n)
=
\frac{|n|}{\prod_{q\le2R}q^{v_q(n)}}.
\tag{19}
\]

Because \(K_r/d_r\) is already reduced and \(p>2R\ge r\) never divides
\(d_r\), replacing \(K_r\) by \(\operatorname{sat}_{>2R}(K_r)\) does not change
any valuation in (18).  This strips **all** small-prime content, not merely an
ad hoc list of visible factors.

The verifier reports

\[
H_K(R)=\sum_{R<r\le2R}\log\max(1,|K_r|),
\]

\[
H_K^{\rm sat}(R)=
\sum_{R<r\le2R}\log\max(1,\operatorname{sat}_{>2R}(K_r)),
\]

and the exact target height from (18).  The remote exact run produced:

| R | raw \(H_K/R^2\) | saturated \(H_K^{sat}/R^2\) | target/R² | raw \(\Xi\)/R² | target rows |
|---:|---:|---:|---:|---:|---:|
| 8 | 8.630968350 | 8.515956789 | 0.044268959 | 9.802723869 | 1 |
| 16 | 9.174146750 | 9.147579238 | 0 | 10.103221465 | 0 |
| 32 | 9.369866830 | 9.361902890 | 0 | 10.296596390 | 0 |
| 64 | 9.532169880 | 9.528850052 | 0 | 10.414895411 | 0 |
| 128 | 9.649731045 | 9.646820987 | 0 | 10.484862951 | 0 |
| 256 | 9.703487841 | 9.702098228 | 0.000117689 | 10.525258824 | 1 |

Every row of the target table was independently computed twice,
from \(\gcd(b_r,\Xi_r)\) and from \(\gcd(b_r,K_r)\), saturated at \(2R\), and
exact equality was asserted.

The two nonzero target rows are exactly what one expects from the locked finite
regressions: \((17,13)\) appears in the \(R=8\) block and \((2237,492)\) in the
\(R=256\) block.

**Finite conclusion only.**  The saturation does not reveal a hidden height
saving for this carrier; at \(R=256\) it changes the normalized height from
9.70349 to only 9.70210.  This is strong negative evidence, not an asymptotic
lower bound for the saturated carrier.

---

## 7. A rigorous no-go for the unsaturated / primitive-numerator rowwise carrier

The raw failure is not merely numerical.

Let

\[
\alpha=17+12\sqrt2,\qquad \rho=\alpha^{-1}=17-12\sqrt2.
\]

The standard Poincare asymptotic for the Apéry recurrence is

\[
b_r=c_b\alpha^r r^{-3/2}(1+O(r^{-1})),\qquad c_b>0.
\tag{20}
\]

Consequently \(F(\rho)=\sum b_r\rho^r\) converges to a positive finite value.
Since

\[
D(t)=(1-\alpha t)(1-\alpha^{-1}t),
\]

\(D\) has a simple zero at \(t=\rho\).  Thus

\[
g(t)=\frac1{F(t)^2\sqrt{D(t)}}
=c_g(1-t/\rho)^{-1/2}+O(1)
\]

in the leading singular scale, with \(c_g>0\), and singularity transfer gives

\[
g_r=c'_g\alpha^r r^{-1/2}(1+o(1)),\qquad c'_g>0.
\tag{21}
\]

Now (12) implies, for large \(r\),

\[
5g_r
\le
\bigl(r^3+P(r-1)+(r-1)^3\bigr)
\max(|\kappa_r|,|\kappa_{r-1}|,|\kappa_{r-2}|).
\]

The coefficient in parentheses is \(O(r^3)\), so by (21)

\[
\max(|\kappa_r|,|\kappa_{r-1}|,|\kappa_{r-2}|)
\ge c\alpha^r r^{-7/2}
\tag{22}
\]

for some \(c>0\).

Choose \(r=R+3,R+6,\ldots\le2R\).  The triples
\(\{r-2,r-1,r\}\) are disjoint and all lie in \((R,2R]\).  From each triple,
(22) supplies one index \(j\) with

\[
\log\max(1,|\kappa_j|)
\ge r\log\alpha-O(\log R).
\]

Since \(K_j=d_j\kappa_j\) and \(d_j\ge1\),
\(|K_j|\ge|\kappa_j|\).  Summing over \(\asymp R\) disjoint triples gives

\[
\boxed{
\sum_{R<j\le2R}\log\max(1,|K_j|)=\Omega(R^2).}
\tag{23}
\]

So **the primitive numerator carrier furnished by the new identity cannot have
subquadratic raw block height**.  This is a theorem, unlike the finite
saturation table.

Equation (23) does *not* survive arbitrary removal of all primes \(\le2R\): an
integer of large Archimedean size can in principle be very smooth.  Therefore I
do not claim a saturated \(\Omega(R^2)\) no-go from (23).  The exact table above
is the evidence for saturation, and the asymptotic saturated-height question
remains open for this carrier.

---

## 8. Precisely scoped failure and the smallest surviving identity class

### What has been achieved

The generic row-scale obstruction has been removed.  The pair
\((b_r,\kappa_r)\) sits in a canonical unit-Casoratian frame and (10) converts
the transverse coordinate \(\Xi_r\) into an adjacent determinant with literal
determinant \(1/r^3\).  High target content is exactly the high content of
\(\gcd(b_r,K_r)\).

This is materially stronger than rephrasing \(b_r+Y\Xi_r\) or the pair minors.
It fixes the scale and preserves prime-power target valuations.

### What has failed

If one now takes one integer \(K_r\) per row and multiplies their heights, the
quadratic barrier returns.  Raw quadratic height is proved by (23), and maximal
small-prime saturation shows no finite saving in the tested blocks.

The Green identity (14) is not another saving: it is the telescoped form of the
same unit-Casoratian identity.

### Smallest surviving algebraic class

The smallest exact data that survive this audit are therefore

1. the adjacent determinant
   \[
   \Xi_r=r^3(b_{r-1}\kappa_r-b_r\kappa_{r-1});
   \]
2. the banded inhomogeneous source law
   \[
   A_r\kappa=-5g_r\quad(r\ge2);
   \]
3. their interval Green identity (14)-(15).

Any further **algebraic** height gain has to use these relations *before*
forming one scalar per row: for example, a multirow resultant/Fitting minor of
the banded system whose determinant couples many target labels without
multiplying an exponential primitive numerator for every row.  I have not
found such a determinant here, and no theorem below claims it exists.

A modular/eta/Hecke route also remains genuinely different: it could constrain
the normalized \(K_r\) residues globally across \(r\), rather than paying their
Archimedean size row by row.  The present audit neither proves nor refutes that
possibility.

What is ruled out is narrower and precise: **the new normalized identity does
not by itself turn the target into a subquadratic product of rowwise primitive
integer carriers.**

---

## 9. BFH scope discipline

Nothing in the derivation uses an assertion of the form “BFH controls all
\(p\mid b_r\).”  It does not.

For a prime \(p>2R\), the chain is

\[
p\mid b_r,\Xi_r
\Longleftrightarrow
p\mid b_r,K_r,
\]

with equality of the corresponding gcd valuations by (17).  This is an exact
local algebraic equivalence on a row already in the transverse target.

If BFH is subsequently inserted, it may only bound the aggregate set of such
simultaneous rows/primes in the form BFH actually supplies.  No step here turns
that row-aggregate information into a pointwise statement, and no step assigns
BFH control to a prime that merely divides \(b_r\).

---

## 10. Verifier contract

Run

```bash
python problems/3.2/research/scripts/q7694_scale_sensitive_crossrow_verify.py --N 520
```

The script independently reconstructs \(b\), \(\sqrt D\), \(F^2\sqrt D\),
\(g\), \(\Xi\), and \(\kappa\).  It then reconstructs \(u,\Phi\) separately and
checks:

- the unit Casoratian (4);
- the frame (6) and synchronized row (8);
- both inverse identities (10)-(11);
- the source law (12) and the nonzero boundary defect (13);
- the Green identity (14);
- reduced-denominator prime support \(\le r\) through row 520;
- exact high-target equality after complete \(q\le2R\) saturation for every row
  in the reported blocks;
- common locked pairs `(17,13)`, `(2237,492)`;
- nontransverse control `(11,5)`;
- explicit non-use of the failed `R=1024,p=4013` six-slope passport and of all
  reflected-depth laws.

The remote run completed successfully with final line

```text
Q7694 VERIFIER: PASS
```

## Bottom line

The actual new identity is (10).  It is Apéry-specific, cross-row,
scale-sensitive, denominator-safe for \(p>2R\), and exactly target-preserving.
It therefore closes the *projective scaling* objection to using a transverse
coordinate.

It does **not** close the high-prime content theorem: after proper denominator
reduction and maximal small-prime saturation, its natural rowwise numerator
carrier still exhibits quadratic height, and raw quadratic height is provably
unavoidable.  The next viable algebraic move must compress many rows at once
inside a banded Fitting/resultant object, or else use genuinely modular/Hecke
information; multiplying the normalized \(K_r\) one row at a time is not an
improvement.