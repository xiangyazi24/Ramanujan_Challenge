# Codex Result 002: trilinear CRT, Fourier non-concentration, and \(M_p\)

## Verdict

- **Task A:** The trilinear Fourier identity and the unit-group
  orthogonality identity both verify. The latter requires
  \(v\in(\mathbb Z/pp'\mathbb Z)^\times\) and \(k,k'\ne0\).
- **Task B:** Lemma 1 holds with the explicit absolute constant
  \(C=4\pi^2\). A direct finite-group Fejér argument needs a two-half
  partition of \(\mathcal Z_p\) to eliminate cyclic wrap-around gaps.
  A standalone proof is in [lemma1_proof.tex](lemma1_proof.tex).
- **Task C:** After matching the benchmark to the requested
  off-diagonal signed-frequency statistic, the aggregate ratios are
  \(0.9924\) for \(K=\lfloor\sqrt p\rfloor\) and \(0.9661\) for
  \(K=\lfloor p^{1/3}\rfloor\).

## Task B: Fourier non-concentration

Let

\[
F_p(k)=\sum_{r\in\mathcal Z_p}e_p(kr),
\qquad Z=Z(p),
\]

and define the ordinary, non-wrapping gap count

\[
A_p(h)
=\#\{r:0\le r<r+h<p,\ r,r+h\in\mathcal Z_p\}.
\]

The no-consecutive-zeros theorem gives \(A_p(1)=0\). The gap
polynomial \(N_h\), its nonvanishing modulo \(p\), and
\(\deg N_h=3(h-1)\) give

\[
A_p(h)\le3(h-1)\qquad(2\le h<p).
\tag{1}
\]

### Why the direct Fejér argument needs correction

A finite-group Fejér expansion applied to all of \(\mathcal Z_p\)
counts cyclic differences. Its coefficient at a small positive \(h\)
therefore sees both ordinary gaps \(h\) and \(p-h\). Bound (1) controls
the former but gives no useful estimate for the latter.

Partition

\[
\mathcal Z_p=Z_0\sqcup Z_1,\qquad
Z_0\subseteq[0,(p-1)/2],\qquad
Z_1\subseteq[(p+1)/2,p-1],
\]

and put \(F_i(k)=\sum_{r\in Z_i}e_p(kr)\). Then

\[
|F_p(k)|^2\le2|F_0(k)|^2+2|F_1(k)|^2.
\tag{2}
\]

Within either half, a cyclic difference of absolute size below \(p/4\)
is the corresponding ordinary integer difference.

### Fejér majorant

First assume \(K\le p/8\), and set

\[
H=\left\lfloor\frac p{4K}\right\rfloor,
\qquad
\frac p{8K}\le H\le\frac p{4K}.
\]

For the finite Fejér kernel

\[
\mathcal F_H(k)
=\frac1H\left|\sum_{u=0}^{H-1}e_p(ku)\right|^2
=\sum_{|h|<H}\left(1-\frac{|h|}{H}\right)e_p(kh),
\]

the sine formula gives, whenever \(1\le|k|\le K\),

\[
\mathcal F_H(k)\ge\frac{4H}{\pi^2}.
\tag{3}
\]

Finite Fourier orthogonality, applied separately to \(Z_0\) and
\(Z_1\), now yields from (2)--(3)

\[
\begin{aligned}
\sum_{1\le|k|\le K}|F_p(k)|^2
&\le \frac{\pi^2p}{2H}
\left(Z+2\sum_{h=1}^{H-1}A_p(h)\right)\\
&\le \frac{\pi^2p}{2H}(Z+3H^2)\\
&\le4\pi^2KZ+\frac{3\pi^2}{4}\frac{p^2}{K}.
\end{aligned}
\]

For \(p/8<K\le(p-1)/2\), Parseval gives

\[
\sum_{k\bmod p}|F_p(k)|^2=pZ,
\]

so the requested sum is at most \(pZ\le8KZ\). Combining the two
ranges proves

\[
\boxed{
\sum_{1\le|k|\le K}|F_p(k)|^2
\le4\pi^2\left(KZ(p)+\frac{p^2}{K}\right)
}
\]

for every prime \(p\ge7\) and \(1\le K\le(p-1)/2\).

The standalone file [lemma1_proof.tex](lemma1_proof.tex) contains the
fully labelled proof. It compiles independently. A numerical check at
four values of \(K\) for every prime \(p\le2000\) passed; the largest
observed value of

\[
\frac{\sum_{1\le|k|\le K}|F_p(k)|^2}
{KZ(p)+p^2/K}
\]

was \(1.19349\), at \((p,K,Z)=(1069,534,6)\).

## Task A: exact trilinear expansion

Put \(L=X^2\), \(n=pq\ell\), and

\[
D_L(\alpha)=\sum_{0\le m<L}e(m\alpha).
\]

Fourier inversion in the three prime fields gives the exact identity

\[
\begin{aligned}
C_{p,q,\ell}(L)
&:=\#\{0\le m<L:
  m\bmod p\in\mathcal Z_p,\
  m\bmod q\in\mathcal Z_q,\
  m\bmod\ell\in\mathcal Z_\ell\}\\
&=\frac1{pq\ell}
\sum_{a\bmod p}\sum_{b\bmod q}\sum_{c\bmod\ell}
F_p(a)F_q(b)F_\ell(c)
D_L\!\left(-\frac ap-\frac bq-\frac c\ell\right).
\end{aligned}
\tag{4}
\]

The zero-frequency term is

\[
\frac{L}{pq\ell}Z(p)Z(q)Z(\ell),
\]

so the exact Fourier error is

\[
\boxed{
E_{p,q,\ell}
=\frac1{pq\ell}
\sum_{(a,b,c)\ne(0,0,0)}
F_p(a)F_q(b)F_\ell(c)
D_L\!\left(-\frac ap-\frac bq-\frac c\ell\right).
}
\tag{5}
\]

The script
[task002_trilinear_verify.py](task002_trilinear_verify.py) computes the
left side by integer CRT, stores the main term and error as exact
rationals, and independently reconstructs (4) using a length-\(pq\ell\)
DFT.

For \(X=50\), the primes are

\[
53,59,61,67,71,73,79,83,89,97.
\]

Only four have nonempty zero sets:

\[
\mathcal Z_{59}=\{9,49\},\quad
\mathcal Z_{61}=\{4,56\},\quad
\mathcal Z_{73}=\{2,70\},\quad
\mathcal Z_{97}=\{25,71\}.
\]

All \(\binom{10}{3}=120\) triples were checked. The four nontrivial
rows are:

| \((p,q,\ell)\) | exact count | exact main term | exact \(E\) |
|---|---:|---:|---:|
| \((59,61,73)\) | 0 | \(20000/262727\) | \(-20000/262727\) |
| \((59,61,97)\) | 0 | \(20000/349103\) | \(-20000/349103\) |
| \((59,73,97)\) | 0 | \(20000/417779\) | \(-20000/417779\) |
| \((61,73,97)\) | 0 | \(20000/431941\) | \(-20000/431941\) |

The other 116 triples have zero main term, count, and error. Numerical
DFT diagnostics were:

~~~text
max |Fourier count - exact count|: 1.510e-16
max reconstructed-indicator error: 5.193e-17
max imaginary reconstruction error: 2.723e-17
~~~

The exact-rational output is
[task002_trilinear_X50.csv](task002_trilinear_X50.csv).

### Unit-group orthogonality

The correctly quantified identity is, for distinct primes \(p,p'\) and
nonzero \(k\bmod p\), \(k'\bmod p'\),

\[
\boxed{
\sum_{v\in(\mathbb Z/pp'\mathbb Z)^\times}
F_p(kv^{-1})\overline{F_{p'}(k'v^{-1})}
=Z(p)Z(p').
}
\tag{6}
\]

CRT factors the left side into

\[
\left(\sum_{u\in\mathbb F_p^\times}F_p(ku)\right)
\overline{
\left(\sum_{u'\in\mathbb F_{p'}^\times}F_{p'}(k'u')\right)}.
\]

Since \(0\notin\mathcal Z_p\),

\[
\sum_{u\ne0}F_p(ku)
=\sum_{r\in\mathcal Z_p}\sum_{u\ne0}e_p(kur)
=-Z(p).
\]

The product is therefore \(Z(p)Z(p')\). The script checked this exact
integer calculation for all ordered \(p\ne p'\) in the \(X=50\) block
and every nonzero \(k,k'\): 467,408 checks. Direct numerical
enumeration over \(v\) for representative frequencies had maximum
error \(6.274\times10^{-13}\).

## Task C: multiplicative collision measurement

Let

\[
\mathcal K_K=\{-K,\ldots,-1,1,\ldots,K\}\pmod p,
\qquad L_K=|\mathcal K_K|=2K,
\]

and

\[
H_p(x)
=\#\{(k,r)\in\mathcal K_K\times\mathcal Z_p:kr=x\}.
\]

Then

\[
\sum_{\substack{k,k'\in\mathcal K_K\\k\ne k'}}M_p(k,k')
=\sum_{x\in\mathbb F_p}H_p(x)^2-L_KZ(p).
\tag{7}
\]

The subtracted diagonal is exact because \(M_p(k,k)=Z(p)\). This
reveals two inconsistencies in the benchmark printed in the task:

1. \(K^2Z^2/p+KZ\) contains a diagonal term, while (7) excludes
   \(k=k'\).
2. The signed window has \(2K\), rather than \(K\), frequencies.

For a uniformly random fixed-cardinality subset of
\(\mathbb F_p^\times\), the exact expectation of the literal ordered
off-diagonal statistic is

\[
\mathbb E T_p^{\ne}
=L_K(L_K-1)\frac{Z(p)(Z(p)-1)}{p-2}.
\tag{8}
\]

The exhaustive computation for \(p\le10000\), restricted to the 468
primes with \(Z(p)\ge2\), gives:

| scale | observed off-diagonal | task benchmark | observed/task | exact prediction (8) | observed/predicted |
|---|---:|---:|---:|---:|---:|
| \(K=\lfloor\sqrt p\rfloor\) | 8,740 | 78,308.612 | 0.1116 | 8,806.635 | **0.9924** |
| \(K=\lfloor p^{1/3}\rfloor\) | 544 | 18,066.666 | 0.0301 | 563.117 | **0.9661** |

Restoring the exact diagonal gives total-energy
observed/predicted ratios \(0.99958\) and \(0.99947\), respectively.
Thus the aggregate data agree closely with the fixed-cardinality random
scale.

There is substantial primewise variation. At cube-root scale, 386 of
the 468 active primes have no off-diagonal collision. Conversely,

\[
p=6151,\qquad\mathcal Z_{6151}=\{2460,3690\},
\]

has 24 cube-root-scale and 104 square-root-scale collisions. The ratio
between the two zero positions is \(3/2\), so it is represented by many
small frequency pairs. Its individual corrected cube-root ratio is
\(58.6\). The experiment therefore supports an averaged
multiplicative-correlation law, not a uniform per-prime bound.

The generator is [task002_mp_measure.py](task002_mp_measure.py).
Per-prime data and aggregate summaries are in
[task002_mp_measurements.csv](task002_mp_measurements.csv) and
[task002_mp_summary.json](task002_mp_summary.json).
