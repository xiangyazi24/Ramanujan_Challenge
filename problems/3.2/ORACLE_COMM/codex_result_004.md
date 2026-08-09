# Codex Result 004: two-flip identity, line sections, and (AT)

## Main finding

For all 115 primes \(p\le 2000\) with \(Z(p)\ge2\), every
nonzero-slope line through the origin contains at most one off-diagonal
point of the incidence cloud. Equivalently, throughout this range,

\[
M_p(k,k')\le1
\qquad
(k,k'\in\mathbb F_p^\times,\ k\ne k').
\]

Thus every tested Apéry zero set is multiplicatively Sidon. Its
line-section energy attains the combinatorial minimum

\[
\boxed{E_{\mathrm{line}}(p)=2Z(p)^2-Z(p).}
\]

The zero sets used below were computed with the division-free recurrence
for \(\beta_n=(n!)^3b_n\), and Task B independently checks every one of
them against the original divided Apéry recurrence.

## Task A: two-flip reciprocity

Let

\[
u_p\equiv(q\ell)^{-1}\pmod p,
\qquad
u_\ell\equiv(pq)^{-1}\pmod\ell,
\]

and define the integer

\[
m^*
=s+(r-s)q\ell u_p+(t-s)pq u_\ell.
\]

Reduction modulo the three primes gives

\[
m^*\equiv r\pmod p,\qquad
m^*\equiv s\pmod q,\qquad
m^*\equiv t\pmod\ell.
\]

Therefore \(m^*\equiv m(r,s,t)\pmod{pq\ell}\), and

\[
\begin{aligned}
e(km/pq\ell)
&=e(km^*/pq\ell)\\
&=e_p\!\left(k(r-s)(q\ell)^{-1}\right)
  e_\ell\!\left(k(t-s)(pq)^{-1}\right)
  e(ks/pq\ell).
\end{aligned}
\]

The script
[task004_two_flip_verify.py](task004_two_flip_verify.py) exhausts all
\(5\cdot7\cdot11=385\) residue triples and all
\(-25\le k\le25\):

~~~text
moduli=(5,7,11), product=385
exact checks=19635
max floating phase error=3.000e-14
~~~

The primary check is the exact integer congruence; the floating-point
phase comparison is only a secondary diagnostic.

## Task B: line-section statistics

### Exact dictionary

For an ordered zero pair \((r,r')\in\mathcal Z_p^2\), set

\[
h\equiv r'-r\pmod p,
\]

representing the zero residue by \(h=p\). The no-consecutive-zeros
lemma excludes \(h=1\), while \(h=p\) is precisely the diagonal
\(r'=r\). Hence

\[
|\mathcal C_p|=Z(p)^2.
\]

Also \(0\notin\mathcal Z_p\), so each cloud point lies on a unique
finite-slope line \(L_c:h=cr\); the vertical line is empty. Since
\(r'=r+h\),

\[
kr=k'r'
\quad\Longleftrightarrow\quad
k'h=(k-k')r.
\]

For \(k'\ne0\), \(M_p(k,k')\) is therefore the section size on the line
of slope

\[
c=\frac{k-k'}{k'}.
\]

When \(k,k'\ne0\) and \(k\ne k'\), the admissible slopes exclude both
\(0\) and \(-1\).

### Observed distribution

The script
[task004_line_sections.py](task004_line_sections.py) computed every
prime \(p\le2000\). There are 115 primes with \(Z(p)\ge2\):

| \(Z(p)\) | number of primes | energy per prime | largest nonzero section |
|---:|---:|---:|---:|
| 2 | 94 | 6 | 1 |
| 4 | 17 | 28 | 1 |
| 6 | 4 | 66 | 1 |

For each active prime, the finite-slope pencil consists exactly of:

- one slope-\(0\) line containing the \(Z(p)\) diagonal points;
- \(Z(p)(Z(p)-1)\) distinct nonzero lines containing one point each;
- all remaining finite-slope lines containing no points.

Aggregating all 115 primes gives:

| section size | finite-slope lines | full projective pencil |
|---:|---:|---:|
| 0 | 102,356 | 102,471 |
| 1 | 512 | 512 |
| 2 | 94 | 94 |
| 4 | 17 | 17 |
| 6 | 4 | 4 |

The extra 115 zero sections in the last column are the empty vertical
lines. The section sizes \(2,4,6\) occur only on the forced
slope-\(0\) lines. In particular,

\[
\sum_{c\ne0}|\mathcal C_p\cap L_c|^2
=Z(p)(Z(p)-1),
\]

and

\[
\sum_p E_{\mathrm{line}}(p)=1304,\qquad
\sum_p\sum_{c\ne0}|\mathcal C_p\cap L_c|^2=512.
\]

The collision-free finding can also be written as

\[
\frac{r_1'}{r_1}=\frac{r_2'}{r_2},
\quad
r_i,r_i'\in\mathcal Z_p,\quad r_i\ne r_i'
\quad\Longrightarrow\quad
(r_1,r_1')=(r_2,r_2')
\]

for every tested prime. This is computational evidence, not a proof
beyond \(p=2000\).

### Random comparison

Conditioning only on \(r\ne0\), a uniform cloud of
\(M=Z(p)^2\) points has \(p\) finite-slope bins and expected energy

\[
M+\frac{M(M-1)}p.
\]

Its aggregate prediction is \(812.1403\), versus the observed \(1304\).
This null model is inappropriate because it ignores the
\(Z(p)\) diagonal points forced onto slope \(0\).

After fixing the diagonal, the remaining
\(M_*=Z(p)(Z(p)-1)\) ordered pairs have \(p-2\) admissible slopes:
slope \(0\) is diagonal, while slope \(-1\) would force \(r'=0\).
The conditioned random prediction is therefore

\[
E_{\mathrm{cond}}(p)
=Z(p)^2+M_*+\frac{M_*(M_*-1)}{p-2}.
\tag{1}
\]

| statistic | observed | conditioned prediction | observed/predicted |
|---|---:|---:|---:|
| full line energy | 1304 | 1314.4735 | 0.9920 |
| nonzero-slope energy | 512 | 522.4735 | 0.9800 |

The aggregate ratio is close to one because the predicted collision
excess is only \(10.4735\). The sharper fact is that this excess is
entirely absent: no two off-diagonal ordered pairs have the same ratio
within any tested prime.

Per-prime data and aggregate output are in
[task004_line_sections.csv](task004_line_sections.csv) and
[task004_line_summary.json](task004_line_summary.json).

## Task C: precise (AT) and numerical test

Write

\[
K_X(m)=\#\{p\in(X,2X]:m\bmod p\in\mathcal Z_p\},
\qquad
\lambda_X=\sum_{X<p\le2X}\frac{Z(p)}p,
\]

and

\[
M_X=\max_{0\le m<X^2}K_X(m).
\]

A quantifier-safe form of the proposed anti-tail hypothesis is:

> **(AT).** For every \(\varepsilon>0\), there are constants
> \(C_\varepsilon,X_\varepsilon>0\) such that every dyadic
> \(X\ge X_\varepsilon\) satisfies
> \[
> M_X\le C_\varepsilon
> \left(1+X^{2/3+\varepsilon}\lambda_X\right).
> \]

The additive \(1\) handles empty or exceptionally sparse blocks. The
available second-moment argument gives only

\[
M_X\ll1+X\lambda_X.
\]

The exhaustive scan of every integer \(0\le m<X^2\) gives:

| \(X\) | \(\lambda_X\) | \(M_X\) | \(X^{2/3}\lambda_X\) | \(M_X/(X^{2/3}\lambda_X)\) | \(X\lambda_X\) | \(M_X/(X\lambda_X)\) |
|---:|---:|---:|---:|---:|---:|---:|
| 16 | 0.287426 | 2 | 1.825 | 1.096 | 4.599 | 0.435 |
| 32 | 0.169520 | 2 | 1.709 | 1.171 | 5.425 | 0.369 |
| 64 | 0.048016 | 2 | 0.768 | 2.603 | 3.073 | 0.651 |
| 128 | 0.157665 | 3 | 4.004 | 0.749 | 20.181 | 0.149 |
| 256 | 0.126286 | 3 | 5.092 | 0.589 | 32.329 | 0.0928 |
| 512 | 0.092772 | 3 | 5.937 | 0.505 | 47.499 | 0.0632 |
| 1024 | 0.082244 | 4 | 8.355 | 0.479 | 84.218 | 0.0475 |
| 2048 | 0.077521 | 4 | 12.502 | 0.320 | 158.763 | 0.0252 |
| 4096 | 0.076828 | 4 | 19.668 | 0.203 | 314.689 | 0.0127 |

The observed maximum is always \(2\), \(3\), or \(4\). At \(X=4096\),
it is \(78.7\) times smaller than \(X\lambda_X\) and \(4.92\) times
smaller than \(X^{2/3}\lambda_X\). The ratio \(2.603\) at \(X=64\)
comes from an unusually sparse block and is absorbed by the constant
and additive \(1\) in the precise formulation.

The data strongly support (AT), and the observed maxima look closer to
a sparse occupancy extreme-value scale than to either power scale.
They do not by themselves supply the missing proof of the
\(X^{1/3}\) improvement.

The exhaustive data and generator are
[task004_at_measurements.csv](task004_at_measurements.csv),
[task004_at_summary.json](task004_at_summary.json), and
[task004_at_measure.py](task004_at_measure.py).
