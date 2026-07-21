# Oracle B exploration: coefficient zeros, fixed anchors, and the Hasse square

This file records finite reconnaissance only.  No numerical observation below
is used as a proof of cross-prime decorrelation.  The rigorous output and the
quantifier audit are in `oracleB_result.tex`.

## Reproduction

The full run was

```text
python3 problems/3.2/oracleB_explore.py \
  --limit 20000 --permutations 999 --hasse-mode all \
  --json-output /tmp/oracleB_full.json
```

The ordinary Python phase recomputes the coefficient zero-set

\[
\mathcal Z_p=\{0\le j<p:b_j\equiv0\pmod p\}
\]

using the division-free recurrence for \((n!)^3b_n\).  Every recomputed pair
was compared with `data_zp_pairs.bin`; the file has SHA-256
`8746d0b400c1b669b001eae955c602908a10c9ee4cb3cac62c6676ea2ddd874d`.
The Sage worker independently constructs the divided coefficients of
\(H_p(t)=\sum_{j<p}b_jt^j\), checks the polynomial square identity over
\(\mathbf F_p\), and factors only a sparse, predetermined sample of square
roots.  All 2260 primes \(5\le p\le20000\) were checked.

## 1. Coefficient zero-sets

The exact histogram is:

| \(Z(p)=|\mathcal Z_p|\) | 0 | 1 | 2 | 4 | 6 | 8 |
|---:|---:|---:|---:|---:|---:|---:|
| number of primes | 1356 | 2 | 695 | 176 | 27 | 4 |

There are 904 nonempty zero-sets.  The odd cases are exactly \(p=11,3137\)
in this range.  All sets passed reflection and no-consecutive-zero checks.

The phrase “first zero, equivalently the whole zero-set” is false outside the
doublet subfamily.  For example,

\[
\mathcal Z_{181}=\{19,47,133,161\}.
\]

Thus \(\rho_p=\min\mathcal Z_p\) is undefined for 1356 primes and fails to
determine \(\mathcal Z_p\) for 209 of the 904 active primes.  In the tables
below, \(h_p=p-1-2\rho_p\) always refers only to the outermost reflected pair.

For all active primes, the empirical quantiles are:

| quantile | 0 | 0.10 | 0.25 | 0.50 | 0.75 | 0.90 | 1 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| \(\rho_p/p\) | 0.00034 | 0.03548 | 0.09205 | 0.21430 | 0.34001 | 0.43382 | 0.49984 |
| \(h_p/p\) | 0 | 0.13096 | 0.31973 | 0.57048 | 0.81580 | 0.92896 | 0.99926 |

To account for multi-pair sets, the null experiment keeps each prime's actual
number of reflected pairs and its center status, but chooses the noncentral
pairs uniformly.  The KS distance for \(2\rho_p/(p-1)\) is 0.03436; the
distance for \(h_p/(p-1)\) is the same.  Restricting to the 695 genuine
doublets gives KS distance 0.03894 from the continuous uniform limit, with the
usual asymptotic diagnostic p-value 0.2379.  The p-value ignores the very fine,
prime-dependent discrete grids and is reported only as reconnaissance.

## 2. Fixed anchors

For each displayed fixed nonnegative integer \(c\), the run counted primes
with \(c\bmod p\in\mathcal Z_p\) and separately checked that every hit divides
the fixed integer \(b_c\):

| \(c\) | hit primes up to 20000 | decimal digits of \(b_c\) |
|---:|:---|---:|
| 0 | none | 1 |
| 1 | 5 | 1 |
| 2 | 73 | 2 |
| 3 | 5, 17 | 4 |
| 5 | 11, 14891 | 6 |
| 10 | 19, 41 | 14 |
| 20 | 17, 163 | 29 |
| 50 | none | 74 |
| 100 | 211 | 150 |
| 1000 | 31, 2411 | 1526 |

The maximum cumulative count among these anchors, compared with \(\pi(X)\),
is:

| \(X\) | 100 | 500 | 1000 | 5000 | 10000 | 20000 |
|---:|---:|---:|---:|---:|---:|---:|
| maximum count | 2 | 2 | 2 | 2 | 2 | 2 |
| maximum count / \(\pi(X)\) | 0.0800 | 0.0211 | 0.0119 | 0.0030 | 0.0016 | 0.00088 |

The plateau is explained exactly by Lucas, not statistically: a fixed anchor
hit implies \(p\mid b_c\).  This finite computation merely checks examples of
the theorem in `oracleB_result.tex`.

## 3. Small-modulus correlation of outer gaps

Pearson's statistic was calibrated by 999 deterministic-seed permutations of
the gap labels.  Odd moduli avoid the automatic evenness of \(h_p\).  The
literal all-active test, using the gap of the outermost reflected pair, is:

| \(q\) | sample size | Cramér \(V\) | permutation p-value |
|---:|---:|---:|---:|
| 3 | 904 | 0.0469 | 0.398 |
| 5 | 903 | 0.0628 | 0.559 |
| 7 | 903 | 0.0878 | 0.257 |
| 11 | 902 | 0.1015 | 0.675 |
| 13 | 902 | 0.1128 | 0.627 |

For comparison, restricting to \(Z(p)=2\) is conceptually cleaner because
only there does the first zero determine the entire noncentral reflected
zero-set:

| \(q\) | sample size | Cramér \(V\) | permutation p-value |
|---:|---:|---:|---:|
| 3 | 695 | 0.0289 | 0.716 |
| 5 | 694 | 0.0760 | 0.440 |
| 7 | 694 | 0.0992 | 0.278 |
| 11 | 694 | 0.1222 | 0.401 |
| 13 | 694 | 0.1180 | 0.957 |

No tested table rejects independence.  This is much weaker than AMTD: fixed
small moduli do not resolve moving collision hyperplanes or minor arcs whose
denominators grow with the dyadic scale.

## 4. The Hasse square and its complexity

For every tested prime, the exact identity

\[
H_p(t)=\Delta(t)^{\varepsilon_p}B_p(t)^2,
\qquad \Delta(t)=t^2-34t+1,
\]

held, with

\[
\deg B_p=\frac{p-1-2\varepsilon_p}{2}.
\]

The computed \(B_p\) was squarefree and coprime to \(\Delta\) at all 2260
primes.  These latter two finite checks are not substituted for a published
proof.  The normalized constant-one square root is not uniformly
anti-palindromic: it is reciprocal for
\(p\bmod24\in\{1,11,17,19\}\) and anti-reciprocal for
\(p\bmod24\in\{5,7,13,23\}\) throughout the tested range.

For every tested prime, the worker also verified

\[
B_p(t)\mid t^{p^2}-t
\]

and used \(\deg\gcd(B_p,t^p-t)\) to extract the exact numbers of linear and
quadratic irreducible factors.  Thus all irreducible factors have degree one
or two at every one of the 2260 tested primes.  Complete Sage/FLINT
factorizations, including the factor polynomials rather than just their degree
types, were performed only at predetermined cross-check samples.  Selected
rows are:

| \(p\) | \(\deg B_p\) | linear factors | quadratic factors | number of factors |
|---:|---:|---:|---:|---:|
| 101 | 50 | 10 | 20 | 30 |
| 503 | 250 | 46 | 102 | 148 |
| 1009 | 504 | 20 | 242 | 262 |
| 5003 | 2501 | 67 | 1217 | 1284 |
| 10007 | 5002 | 188 | 2407 | 2595 |
| 19997 | 9998 | 118 | 4940 | 5058 |

The degree-one/two statement is an exact finite-range result, not an
extrapolation to larger primes.  The conclusion relevant to B2 does not
depend on this factor pattern: the degree of \(B_p\) itself grows linearly.
The locus \(B_p(t)=0\) in \(\mathbf A^1\) is zero-dimensional, so it has no
genus.  If one instead manufactures the curve \(y^2=B_p(t)\), then its genus is
\(\lfloor(\deg B_p-1)/2\rfloor=\Theta(p)\), equal to 4998 for \(p=19997\).
Neither construction has uniformly bounded geometric complexity.

There is also a source-level distinction.  Theorem 1.2 of
Caruso--Fürnsinn--Vargas-Montoya--Zudilin proves the displayed
square/quadratic-times-square factorization.  Its proof uses

\[
f_\alpha(t(x))=(1+x)h(x)^2,
\qquad t=\frac{x(1-8x)}{1+x},
\]

and the quadratic extension with discriminant \(t^2-34t+1\).  The paper does
not formulate the result in terms of Hasse invariants, Picard groups, or an
elliptic pencil, and Theorem 1.2 does not state that \(B_p\) is squarefree.
Those stronger interpretations need a separate reference or proof.

## 5. The decisive object mismatch

The factorization controls evaluation roots of \(H_p\).  Oracle B needs zero
coefficients.  They are not the same statistic.  At \(p=7\), exact arithmetic
gives

\[
\begin{aligned}
H_7(t)
 &=1+5t+3t^2+3t^3+3t^4+5t^5+t^6\\
 &=(t-1)^2(t^2+1)^2\pmod7.
\end{aligned}
\]

Thus \(\mathcal Z_7=\varnothing\), while \(H_7\) and its square root have the
evaluation root \(t=1\).  In general,

\[
b_j=-\sum_{t\in\mathbf F_p^\times}H_p(t)t^{-j}
\qquad(1\le j\le p-2),
\]

so \(\mathcal Z_p\) is a zero-set in the finite Mellin dual, not the root
divisor of \(B_p\).

The explicit equation proposed in the oracle specification,

\[
E_t:y^2=x(x-1)(x-t(1-t)),
\]

is not identified with the Apéry truncation in the cited material.  In the
printed \(t\)-coordinate its raw discriminant is supported on
\(t=0,1,\infty\) and \(t^2-t+1=0\), not on
\(t^2-34t+1=0\).  At \(p=7\), substituting
\(\lambda=t(1-t)\) into the Deuring polynomial
\(1+2\lambda+2\lambda^2+\lambda^3\) gives the root \(t=4\), whereas the
Apéry \(H_7\) has root \(t=1\).  This does not rule out an undisplayed base
change and gauge; it proves that the equation as printed cannot simply be
declared to produce \(B_p\), much less the coefficient-first-zero \(\rho_p\).

## 6. What the computation does and does not support

The data are consistent with a reflected random model and show no detectable
small-modulus dependence.  They do not prove minor-arc cancellation,
two-prime dispersion, or AMTD.  In particular:

- fixed \(h\) Chebotarev for a gap polynomial does not cover
  \(h=h_p=\Theta(p)\);
- fixed \(q\) marginal equidistribution does not control two-prime collision
  hyperplanes;
- a fixed geometric family may have fixed rank and conductor while its
  characteristic-\(p\) Hasse section has degree \(\Theta(p)\);
- the literal fixed-anchor theorem has quantifiers
  \(\forall c\,\exists C_c\), whereas the no-go star chooses
  \(c=c_N\) after the scale \(N\).

The exact missing input remains a horizontal, two-characteristic crystalline
Mellin anti-concentration estimate, such as the marked dispersion statement
isolated in `oracleB_result.tex`.
