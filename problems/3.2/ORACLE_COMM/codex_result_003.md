# Codex Result 003: multiplicative-correlation numerics

## Normalization audit

The two random predictions stated in the task do not match the
statistics being summed.

For \(1\le k,k'\le K\), define

\[
T_p^{\ne}(K)
=\sum_{\substack{1\le k,k'\le K\\k\ne k'}}M_p(k,k').
\]

If

\[
H_{p,K}(x)
=\#\{(k,r):1\le k\le K,\ r\in\mathcal Z_p,\ kr=x\},
\]

then exactly

\[
T_p^{\ne}(K)=\sum_xH_{p,K}(x)^2-KZ(p).
\tag{1}
\]

The proposed expression \(K^2Z(p)^2/p+KZ(p)\) is a
diagonal-containing total-energy benchmark, while \(T_p^{\ne}\)
explicitly excludes \(k=k'\).

For a uniformly random fixed-cardinality subset
\(\mathcal Z_p\subset\mathbb F_p^\times\), each fixed pair \(k\ne k'\)
has

\[
\mathbb E M_p(k,k')=\frac{Z(p)(Z(p)-1)}{p-2}.
\]

Consequently the exact expectation of the requested off-diagonal
statistic is

\[
\boxed{
\mathbb E T_p^{\ne}(K)
=K(K-1)\frac{Z(p)(Z(p)-1)}{p-2}.
}
\tag{2}
\]

The multiplicative-energy benchmark \(pZ(p)^2\) is impossible under
the stated definition. Indeed,

\[
\sum_{k=1}^{p-1}M_p(k,1)=Z(p)^2,
\qquad
M_p(k,1)\le Z(p),
\]

because every ordered pair \((r,r')\in\mathcal Z_p^2\) determines the
unique ratio \(k=r'/r\). Therefore

\[
\boxed{
E_{\mathrm{mult}}(p)
=\sum_{k=1}^{p-1}M_p(k,1)^2
\le Z(p)^3.
}
\tag{3}
\]

Since \(Z(p)<p\), the task benchmark \(pZ(p)^2\) exceeds even this
deterministic upper bound.

For comparison with a genuine random model, let \(N=p-1\) and let
\((x)_t=x(x-1)\cdots(x-t+1)\). Among all ratio-energy quadruples in
\((\mathbb F_p^\times)^4\), the numbers having exactly \(t=1,2,3,4\)
distinct entries are respectively

\[
Q_1=N,\quad
Q_2=2N^2-N,\quad
Q_3=2N(N-2),\quad
Q_4=N(N-2)^2.
\]

Thus the exact mean over uniformly random \(Z\)-element subsets is

\[
\boxed{
\mathbb E_{\mathrm{fixed}\ Z}E_{\mathrm{mult}}
=\sum_{t=1}^4Q_t\frac{(Z)_t}{(N)_t}.
}
\tag{4}
\]

The simpler occupancy expression

\[
Z^2+Z(Z-1)
+\frac{Z(Z-1)(Z(Z-1)-1)}{p-2}
\tag{5}
\]

is a useful independent-balls heuristic, but unlike (4) it is not the
exact fixed-subset expectation when \(Z\ge3\).

## Computation

[mp_numerics.py](mp_numerics.py) recomputed every Apéry zero set for
the 666 primes \(7\le p\le5000\). It uses the recurrence for
\(\beta_n=(n!)^3b_n\), and checks the resulting zero set against the
original divided Apéry recurrence for every prime. There are 246 primes
with \(Z(p)\ge2\). All collision counts are exact integers; no floating
Fourier evaluation is used.

The complete per-prime ratios requested in the task are in
[mp_numerics.csv](mp_numerics.csv), with aggregate values in
[mp_numerics_summary.json](mp_numerics_summary.json).

### Small positive frequencies

| scale | observed \(T_p^{\ne}\) | task benchmark | observed/task | exact prediction (2) | observed/predicted | total energy | corrected total prediction | total ratio |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| \(K=\lfloor\sqrt p\rfloor\) | 974 | 28,217.680 | 0.0345 | 1,064.206 | **0.9152** | 27,512 | 27,602.206 | **0.9967** |
| \(K=\lfloor p^{1/3}\rfloor\) | 60 | 7,244.986 | 0.00828 | 81.113 | **0.7397** | 7,164 | 7,185.113 | **0.9971** |

The very small observed/task ratios come from comparing an
off-diagonal count with a diagonal-dominated benchmark. With the
normalization corrected, the square-root-scale aggregate is within
\(8.5\%\) of the exact fixed-margin prediction. At cube-root scale the
expected aggregate is only \(81.1\), so the observed \(60\) remains a
sparse statistic; 225 of the 246 active primes have no off-diagonal
collision.

The near-one total-energy ratios are less informative because the
forced diagonal \(KZ(p)\) dominates, especially at cube-root scale.

There is substantial primewise variation. For example,

\[
p=3797,\qquad \mathcal Z_p=\{1035,2761\},
\]

has corrected ratios \(10.37\) at square-root scale and \(36.14\) at
cube-root scale, caused by a small rational relation between its two
zero positions. Hence the data support an averaged law, not a uniform
per-prime random estimate.

### Full multiplicative energy

The aggregate comparison is:

| benchmark | observed | prediction | observed/predicted |
|---|---:|---:|---:|
| task \(pZ^2\) | 2,856 | 3,975,216 | 0.000718 |
| independent occupancy heuristic (5) | 2,856 | 2,864.898 | 0.9969 |
| exact fixed-subset expectation (4) | 2,856 | 2,861.212 | **0.9982** |

For 245 of the 246 active primes, the energy attains its absolute
minimum

\[
E_{\mathrm{mult}}(p)=2Z(p)^2-Z(p),
\]

meaning that all ordered ratios \(r'/r\) with \(r\ne r'\) are
distinct. The only exception is

\[
p=3727,\qquad
\mathcal Z_{3727}
=\{99,868,1011,1294,2432,2715,2858,3627\}.
\]

Here \(Z=8\), the minimum energy is \(120\), and the observed energy is
\(124\). Exactly two nonidentity ratio buckets have multiplicity \(2\);
they are inverse to one another:

\[
\frac{1294}{1011}
=\frac{2858}{1294}
=3614\pmod{3727},
\]

\[
\frac{1011}{1294}
=\frac{1294}{2858}
=1847\pmod{3727}.
\]

No nonidentity ratio has multiplicity above \(2\) for \(p\le5000\).

## Interpretation

The full-ratio statistic is not merely random-scale: almost every
tested zero set is multiplicatively Sidon, with one explicit failure at
\(p=3727\). This is strong evidence for low line-section energy. It
does not imply uniform small-window bounds, because a single
nonidentity ratio can have many representations as \(k/k'\) with small
\(k,k'\); the exceptional small-frequency behavior at \(p=3797\)
illustrates this distinction.
