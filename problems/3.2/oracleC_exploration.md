# Oracle C exploration: a fixed toric marked coordinate

This file records exact finite computations for the algebraic marked-coordinate
prerequisite.  It does **not** test, assume, or prove the two-prime dispersion
statement.  The computations were produced by
`python3 problems/3.2/oracleC_explore.py --limit 2000 --workers 1` from the
repository root.  The full exploration scan self-checks the transfer and
support inequalities while producing all displayed aggregate statistics.
`oracleC_verify.py` independently organizes the symbolic identities, every
transfer marker through $p=2000$, and representative finite-statistic checks
as PASS/FAIL groups; it deliberately does not duplicate the aggregate
BM/interpolation scan.

## 1. The successful coordinate

Define the fixed Laurent polynomial

$$
\Lambda(x,y,z)=
\frac{(1+x)(1+y)(1+z)((1+y)(1+z)+xyz)}{xyz}.
$$

Expanding according to the number $k$ of `xyz` choices in the final factor
gives the symbolic identity

$$
\begin{aligned}
\operatorname{CT}\Lambda^n
 &=\sum_{k=0}^n \binom nk^2\binom{2n-k}{n-k}^2\\
 &=\sum_{r=0}^n \binom nr^2\binom{n+r}{r}^2=b_n.
\end{aligned}
$$

The verifier also constructs $\Lambda$ as a Laurent-polynomial dictionary,
raises it by independent exact multiplication, and checks the constant term
against the Apéry sum.  This is separate from the displayed binomial
derivation.

Every exponent of $\Lambda$ in every coordinate lies in $[-1,1]$.  Hence,
for $0\le j\le p-2$, $\Lambda^j$ has no nonzero exponent divisible by
$p-1$.  Three applications of torus orthogonality give

$$
c^{\mathrm{tor}}_{p,j}
:=-\sum_{(x,y,z)\in(\mathbf F_p^\times)^3}\Lambda(x,y,z)^j
=\operatorname{CT}\Lambda^j=b_j\pmod p.
$$

The direct finite-field implementation checks every $j$ in this range for
$p=5,7,11,13,17,19,23,29,31$.  These checks are diagnostics; the exponent
box and orthogonality are the proof.

The endpoint $j=p-1$ must be defined separately as $c_{p,p-1}=1$.  Indeed,
inclusion--exclusion gives

$$
\#V_T(\Lambda)=4p^2-14p+13,
\qquad
\#U(\mathbf F_p)=p^3-7p^2+17p-14,
$$

where $T=(\mathbf G_m)^3$ and $U=T\setminus V(\Lambda)$.  Consequently the
raw power coordinate at $j=p-1$ is $-\#U\equiv14\pmod p$, not $1$; it is
zero at $p=7$.  At $j=0$ we instead interpret $\Lambda^0$ as the constant
Laurent polynomial one on all of $T$, which gives $-(p-1)^3=1$.

### Geometric complexity

For $1\le j\le p-2$, the cyclotomic/Teichmüller lift is the negative
Frobenius supertrace on

$$
R\Gamma_c\left(
U_{\overline{\mathbf F}_p},
\Lambda^*\mathcal K_{\chi_{p,j}}
\right).
$$

The coefficient sheaf has rank one.  In $(\mathbf P^1)^3$, the divisor of
$\Lambda$ has four simple zero components

$$
x=-1,\quad y=-1,\quad z=-1,\quad
X_0(Y_0+Y_1)(Z_0+Z_1)+X_1Y_1Z_1=0
$$

and the six coordinate-boundary pole components.  All Kummer monodromy is
tame and all Swan conductors are zero.  The original ten-component divisor
is not simple normal crossings, so `10` is input component count, not a
fabricated exact higher-dimensional conductor.  One fixed log resolution
and a fixed stratification give an absolute bound for the total compactly
supported Betti dimension, uniformly in $(p,j)$.
This is geometric complexity (rank, divisor, ramification, and Betti/L-degree);
it does not count the degree of the varying cyclotomic coefficient field.

The coordinate change

$$
u=\frac{1+y}{y},\qquad v=\frac{1+z}{z},\qquad
w=-\frac{x}{uv}
$$

identifies $U$ with the fixed toric-arrangement complement

$$
(\mathbf G_m)^3\setminus\{u=1,\ v=1,\ w=1,\ uvw=1\}.
$$

This also explains the point-count polynomial above and
$\chi_c(U)=-3$.  We do not infer concentration in middle cohomology from
the Euler characteristic; resonant boundary strata prevent that shortcut.

The vanishing statement is only after reduction at a prime above $p$.
It is not characteristic-zero trace vanishing.  At $(p,j)=(11,5)$, for
example, the exact quadratic-character sum is $33$, so the negative lifted
trace is $-33\ne0$ although $b_5\equiv0\pmod {11}$.

This completes the literal algebraic G1.  It does not produce a sheaf on a
`j`-line, and it proves no two-prime dispersion estimate.

## 2. Transfer coordinate

Write

\[
 P(n)=34n^3+51n^2+27n+5,
 \qquad
 M(n)=\begin{pmatrix}P(n)&-n^6\\1&0\end{pmatrix},
 \qquad B_n=(n!)^3b_n.
\]

For $T_0=I$ and $T_j=M(j-1)\cdots M(0)$, the scan built every
prefix product and checked

\[
 T_j=\begin{pmatrix}B_j&0\\B_{j-1}&0\end{pmatrix}
 \quad(j\geq1).
\]

Thus the marked entry is `(1,1)`, not `(1,2)`.  For all 301 primes
$5\leq p\leq1999$, all 277045 pairs $0\leq j<p$ passed both the
entry identity and

\[
 (T_j)_{11}=0\quad\Longleftrightarrow\quad b_j=0\pmod p.
\]

The product containing $M(0)$ has rank one.  The full fundamental matrix is

\[
 M(j-1)\cdots M(1)
 \begin{pmatrix}5&1\\1&0\end{pmatrix}
 =\begin{pmatrix}B_j&D_j\\B_{j-1}&D_{j-1}\end{pmatrix},
 \qquad \det=-((j-1)!)^6.
\]

The fixed matrix size is only a fixed local state dimension.  It does not
make the length-$j$ word a bounded-degree family.

## 3. Recurrence complexity in the index

The sequence has a fixed order-two recurrence with polynomial coefficients:

\[
 (n+2)^3b_{n+2}-P(n+1)b_{n+1}+(n+1)^3b_n=0.
\]

The primitive coefficient vector, ordered by $(b_n,b_{n+1},b_{n+2})$ and
in increasing powers of $n$, is

\[
 (1,3,3,1;\ -117,-231,-153,-34;\ 8,12,6,1).
\]

This degree $3$ is minimal in characteristic zero among order-two
recurrences having three polynomial coefficients of a common degree bound.
The degree-$2$ ansatz gives the exact determinant

\[
 -2^{26}3^{12}5^6 7^2\cdot26309\cdot50077\cdot171131\ne0.
\]

Over every $\mathbf F_p$, $11\leq p\leq1999$, its coefficient matrix has
rank $9$.  For every $13\leq p\leq1999$, the degree-$3$ matrix has rank
$11$ and its one-dimensional nullspace is generated by the displayed Apéry
recurrence.

There is no constant-coefficient order-two recurrence even on the full first
characteristic window.  The first two equations would force
$A=45/2$, $C=-79/2$; the next two integral residuals are

\[
 6744,\qquad267120,\qquad\gcd(6744,267120)=24.
\]

Hence no prime $p\geq5$ makes both residuals vanish (the first already
excludes $p=5$).

For the finite prefix $(b_0,\ldots,b_{p-1})$, Berlekamp--Massey gave

| sequence | generic value | exceptions |
|---|---:|---|
| $b_j$ | $L_p=(p+1)/2$ for 297 primes | $p=23,47,67,827$, where $L_p=(p-1)/2$ |
| $B_j$ | $L_p=(p+1)/2$ for 300 primes | $p=19$, where $L_p=(p-1)/2$ |

Selected values are

| $p$ | 5 | 7 | 11 | 23 | 181 | 827 | 1999 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| $L_p(b)$ | 3 | 4 | 6 | 11 | 91 | 413 | 1000 |

This is a finite-prefix diagnostic.  High Berlekamp--Massey complexity is not
a lower bound for an étale-sheaf conductor.

The unique ordinary interpolation polynomial
$I_p\in\mathbf F_p[J]$, $\deg I_p<p$, with $I_p(j)=b_j$ had only the
following two degrees:

| degree | number of primes |
|---:|---:|
| $p-1$ | 144 |
| $p-3$ | 157 |

The second case agreed exactly with $H_p(1)=\sum_jb_j=0$.  Again this
rules out only a direct bounded-degree regular-function representation; it
does not diagnose sheaf conductor.

## 4. Shifted products and the factorial gauge

For

\[
 R_\ell(x)=M(x+\ell-1)\cdots M(x),
\]

an exact leading-coefficient induction gives, for $\ell\geq2$, the entrywise degree
matrix

\[
 \deg R_\ell=
 \begin{pmatrix}3\ell&3\ell+3\\3\ell-3&3\ell\end{pmatrix}.
\]

Writing $A_\ell=(R_\ell)_{11}$, its leading coefficients obey
$\lambda_{\ell+1}=34\lambda_\ell-\lambda_{\ell-1}$ with
$(\lambda_0,\lambda_1)=(1,34)$, and hence remain positive and strictly
increasing.  The other three entries are shifted copies of $A_{\ell-1}$
or $A_{\ell-2}$ multiplied by $-x^6$.  Thus these are exact degrees over
$\mathbf Z$, not upper bounds inferred from finite multiplication.  The
verifier checks the initial symbolic products independently.

The standard recurrence gauge does not remove the growth.  Its marked
coefficient is

\[
 C_h(m)=\frac{N_h(m)}{\prod_{s=2}^{h}(m+s)^3}.
\]

The endpoint identity

\[
 N_h(-s)=(-1)^{s-1}B_{s-1}B_{h-s}\ne0
 \qquad(1\leq s\leq h)
\]

shows that no displayed denominator factor cancels.  Thus this direct gauge
has exactly $h-1$ distinct finite poles, each of order $3$.  This is a
no-go for the direct transfer word and this gauge, not for every possible
geometric realization.

## 5. Raw Mellin support

For every scanned prime, the program evaluated $H_p(t)$ at every
$t\in\mathbf F_p^\times$.  If $g$ is a primitive root, the corrected
group-algebra marker is

\[
 Q_p(X)=-\sum_{r=0}^{p-2}H_p(g^r)X^{-r}.
\]

The omitted $r=0$ term in the older interpolation formula is essential.
For example, at $p=11$, $H_{11}(1)=3$ and $b_5=0$; deleting $t=1$
makes the purported marker take the nonzero value $3$ at that marked zero.

All primes passed

\[
 |\operatorname{supp}Q_p|
 =\#\{t:H_p(t)\ne0\}
 \geq \frac{p-1}{2}-\varepsilon_p
 \geq\frac{p-3}{2}.
\]

Selected exact support counts are

| $p$ | 5 | 7 | 11 | 23 | 31 | 181 | 827 | 1999 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| support | 2 | 5 | 5 | 14 | 27 | 173 | 793 | 1979 |
| proved lower bound | 2 | 3 | 5 | 10 | 15 | 89 | 413 | 999 |

The proof of the lower bound uses the factorization
$H_p=\Delta^{\varepsilon_p}B_p^2$; it does not use squarefreeness or
coprimality.  The same support is the exact order of the minimal cyclic
constant-coefficient recurrence for the **raw** periodic Mellin marker.

## 6. What the original-line Kummer computation does and does not compute

For the Kummer factor alone, with geometric conductor normalized as

\[
 \operatorname{rank}+\#\{\text{singular points}\}
   +\sum\operatorname{Swan},
\]

the trivial character has conductor $1$, while every nontrivial tame
Kummer character on $\mathbf G_m$ has conductor $3$.  This is uniform in
$p$ and $j$.  It is a conductor on the **$t$-line**.

No sheaf on a $j$-line whose trace is $j\mapsto b_j$ is specified.
Consequently the requested number of $j$-singularities and its $L$-degree
are not defined and were not fabricated.  The script reports

```
j_trace_sheaf_conductor=NOT_DEFINED:
  no j-space sheaf with trace j -> b_j is specified
```

The toric construction in Section 1 does not change this answer.  It gives
one uniformly bounded cohomological space for each character $j$; it does
not make $j\mapsto b_j$ the trace function of one sheaf on an additive
$j$-line.  This is exactly why the algebraic marked coordinate can now be
declared constructed while the moving-character dispersion problem remains
untouched.

Conditionally, if a rank-$r$ object $\mathcal V$ on
$U=\mathbf P^1\setminus S$, $|S|=s$, has total Swan conductor
$\Sigma$, then tame Kummer twisting preserves $\Sigma$, and
Grothendieck--Ogg--Shafarevich gives

\[
 \dim H_c^1=r(s-2)+\Sigma+h_\chi,
 \qquad 0\leq h_\chi=\dim H_c^2\leq r.
\]

Thus each individual twist has bounded cohomological dimension.  For the
standard four-puncture Apéry operator and $\Sigma=0$, the generic degree is
$6$ for the rank-$3$ transcendental system and $44$ for the full
rank-$22$ K3 system.  The relevant system is not rank $2$: the Apéry
operator is the rank-$3$ symmetric square of an underlying rank-$2$
operator.

## 7. Smooth-locus boundary failure

The K3 smooth-locus twist is an alternative to, not the source of, the toric
marked coordinate.  It fails to preserve the marked zero set.  At $p=31$,
exact reduction gives

\[
 \Delta(t)=(t-14)(t-20),\qquad
 H_{31}(14)=H_{31}(20)=7,
 \qquad b_8=b_{22}=0.
\]

For $j=8$, the full Mellin sum is zero, but the omitted bad-fiber term is

\[
 7(14^{-8}+20^{-8})=7(19+18)=11\pmod{31}.
\]

Therefore

\[
 \sum_{\substack{t\in\mathbf F_{31}^\times\\\Delta(t)\ne0}}
 H_{31}(t)t^{-8}=20\ne0.
\]

The same values, with the two inverse powers exchanged, occur at $j=22$.
Hence a compactly supported twist on the smooth locus does not merely differ
by an innocuous constant: without geometrically controlled bad-fiber terms it
fails to preserve the marked zero locus.

## Exact program output

```text
ORACLE C FINITE EXPLORATION
limit=2000 workers=1
primes=301 first=5 last=1999 transfer_cases=277045
transfer_marker=(M(j-1)...M(0))[1,1]=B_j (one-based matrix entry)
constant_order2_residuals=6744,267120 gcd=24
BM_b_default=297 exceptions=23,47,67,827
BM_B_default=300 exceptions=19
interpolation_degree_b_drop_counts=0:144,2:157
H_p(1)_zero_count=157
minimum_mellin_support_ratio=2/4 at p=5
raw_support_factorization_bound=PASS
toric_marker_samples=5,7,11,13,17,19,23,29,31 PASS
toric_complexity=dimension=3 half_width=1 rank=1 input_divisor_components=10 Swan=0
toric_U_point_count=q^3-7q^2+17q-14
toric_raw_endpoint=-#U=14(mod p), separate endpoint=1
kummer_factor_geometric_conductor=1(trivial),3(nontrivial)
j_trace_sheaf_conductor=NOT_DEFINED: no j-space sheaf with trace j -> b_j is specified
p31_smooth_boundary=roots=(14, 20) H=(7, 7) sums={8: (0, 11, 20), 22: (0, 11, 20)}
samples=p,Z,L_b,L_B,deg_I_b,support,lower_bound
5,2,3,3,2,2,2
7,0,4,4,4,5,3
11,1,6,6,10,5,5
23,0,11,12,20,14,10
31,2,16,16,28,27,15
181,4,91,91,178,173,89
827,0,413,414,826,793,413
1999,0,1000,1000,1996,1979,999
```
