# Arbitrary-cutoff Lucas blocks and the loss-count quotient filter

## 1. Definitions

Let

\[
L(n,k)=\binom nk\binom{n+k}{k},\qquad
Q_n(t)=\sum_{k=0}^n L(n,k)t^k,
\]

and let

\[
F_m=\sum_{i=0}^m\binom mi^3.
\]

For a center \(c\), put

\[
K_{n,m}(c)=[y^m]Q_n(c+y)
\]

and

\[
g_m(c)=\sum_{i=0}^m\binom mi(-c)^{m-i}F_i.
\]

The cutoff Legendre--Euler transform is

\[
T_{n,J}(c)=\sum_{m=0}^J K_{n,m}(c)g_m(c),
\qquad
\Gamma_{n,J}=\operatorname{content}_c T_{n,J}(c).
\]

The full Strehl--Franel identity and the proved Franel-tail theorem give

\[
T_{n,n}(c)=A_n,\qquad \Gamma_{n,J}\mid A_n.
\]

## 2. The arbitrary-cutoff block identity

Let \(p\) be an odd prime and write

\[
n=qp+r,\qquad J=Ap+B,\qquad 0\le r,B<p.
\]

### Theorem 2.1

Coefficientwise modulo \(p\),

\[
\boxed{
T_{qp+r,Ap+B}(c)
\equiv
A_rT_{q,A-1}(c^p)
+K_{q,A}(c^p)g_A(c^p)T_{r,B}(c).
}
\tag{2.1}
\]

Here \(T_{q,-1}=0\), and a cutoff beyond the polynomial degree is
understood to be saturated.

### Proof

Lucas's theorem applied to both binomial factors in \(L(n,k)\) gives

\[
Q_{qp+r}(t)\equiv Q_q(t^p)Q_r(t)\pmod p.
\]

After substituting \(t=c+y\) and using
\((c+y)^p=c^p+y^p\), the low factor has degree less than \(p\).
Therefore

\[
K_{qp+r,ap+b}(c)
\equiv K_{q,a}(c^p)K_{r,b}(c)\pmod p.
\tag{2.2}
\]

Franel Lucas gives \(F_{ap+b}=F_aF_b\pmod p\). Splitting the index in the
binomial translate into its two base-\(p\) digits gives

\[
g_{ap+b}(c)\equiv g_a(c^p)g_b(c)\pmod p.
\tag{2.3}
\]

Split the defining sum for \(T_{n,J}\) into the complete blocks
\(a=0,\ldots,A-1\) and the final partial block \(a=A\). In every complete
block,

\[
\sum_{b=0}^{p-1}K_{r,b}(c)g_b(c)
=T_{r,r}(c)=A_r.
\]

Equations (2.2) and (2.3) now give (2.1). \(\square\)

The script `q32_allq_cutoff_block.py` checks 2,984 complete polynomial
identities for \(p=5,7,11,13\), every quotient \(q\le3\), every residue,
and every cutoff.

## 3. Folding the residue

For \(0\le r<p\), put

\[
j=\min(r,p-1-r).
\]

If \(r=p-1-j\), then

\[
Q_r(t)\equiv Q_j(t)\pmod p.
\]

Indeed, the coefficients through degree \(j\) are obtained by the two
standard reflected binomial congruences, and all later coefficients have a
carry. Hence

\[
T_{r,B}(c)\equiv T_{j,B}(c)\pmod p.
\tag{3.1}
\]

For \(B\ge j\), this transform is saturated and equals \(A_j\). For
\(B<j\), its degree-\(j\) coefficient is

\[
[c^j]T_{j,B}(c)
=(-1)^B\binom{2j}{j}\binom{j-1}{B},
\tag{3.2}
\]

which is nonzero modulo \(p\). Consequently, under \(p\mid A_j\),

\[
T_{r,B}(c)\equiv0\pmod p
\quad\Longleftrightarrow\quad B\ge j.
\tag{3.3}
\]

For \(q=1\), Theorem 2.1 immediately gives the exact two-plateau profile

\[
p\mid\Gamma_{n,J}
\quad\Longleftrightarrow\quad
j\le J\le p-1
\ \text{ or }\
p+j\le J\le n.
\tag{3.4}
\]

## 4. Effective high digit and the number of losses

For \(0\le q<p\), define

\[
e_p(q)=\min(q,p-1-q).
\]

The effective degree of \(Q_q(t)\) modulo \(p\) is exactly \(e_p(q)\).
Indeed,

\[
L(q,k)\not\equiv0\pmod p
\quad\Longleftrightarrow\quad
k\le q\ \text{ and }\ q+k<p.
\]

It follows that

\[
K_{q,A}(c)\not\equiv0\pmod p
\quad\Longleftrightarrow\quad A\le e_p(q).
\tag{4.1}
\]

For the forward implication, the coefficient contributed by the effective
top degree is a unit. For \(A>e_p(q)\), the required derivative coefficient
of \(Q_q\) is identically zero.

Assume \(p\mid A_j\). In (2.1), the completed-block term vanishes because
\(A_r=A_j=0\pmod p\). Since \(g_A\) is primitive, (3.3) and (4.1) show:

- for \(0\le A\le e_p(q)\), block \(A\) is nonzero for \(B<j\) and zero
  for \(B\ge j\);
- for \(A>e_p(q)\), the transform is zero for the entire block.

Thus the content loses the prime precisely at

\[
J=p,2p,\ldots,e_p(q)p.
\tag{4.2}
\]

Define

\[
R_{n,J}=\operatorname{rad}\Gamma_{n,J}
\]

and the squarefree loss quotient

\[
\Lambda_{n,J}
=\frac{R_{n,J-1}}{\gcd(R_{n,J-1},R_{n,J})}.
\tag{4.3}
\]

Finally put

\[
P_n=\prod_{J=1}^n\Lambda_{n,J}.
\]

### Theorem 4.2

For every bad prime in the two-digit range,

\[
\boxed{v_p(P_n)=\min(q,p-1-q).}
\tag{4.4}
\]

This valuation counts zero-to-nonzero cutoff transitions. It is independent
of the direct/reflected choice of the low digit.

## 5. Extracting the q=1 support

Let

\[
\operatorname{Once}(N)
=\frac{\operatorname{rad}N}
{\operatorname{rad}\gcd(N,N/\operatorname{rad}N)}.
\tag{5.1}
\]

This is the product of the primes occurring in \(N\) to exponent exactly
one. By (4.4), a bad prime occurs in \(\operatorname{Once}(P_n)\) precisely
when

\[
q=1\qquad\text{or}\qquad q=p-2.
\tag{5.2}
\]

The second case is a pointwise \(O(\log n)\) nuisance. It implies

\[
p^2-2p\le n\le p^2-p-1.
\]

These intervals are disjoint for consecutive \(p\), so a fixed \(n\)
belongs to at most one of them. Therefore, after the already controlled
small-prime part, \(\operatorname{Once}(P_n)\) has exactly the q=1 bad-prime
support, up to one near-square-root prime.

If an exact magnitude filter is desired, intersect with

\[
C_n^{\mathrm{top}}
=\frac{\operatorname{lcm}(1,\ldots,n)}
{\operatorname{lcm}(1,\ldots,\lfloor n/2\rfloor)}.
\]

Above \(\sqrt n\), this carrier has exactly the prime support
\((n/2,n]\), so it removes the \(q=p-2\) case and leaves precisely q=1.

This is a support theorem, not a height theorem. The missing estimate is

\[
\log\operatorname{Once}(P_n)=o(n).
\tag{5.3}
\]

No known bound on individual cutoff contents or on the raw product \(P_n\)
implies (5.3).

## 6. Why matching two jumps aliases higher quotients

The earlier same-index construction used

\[
H=\left\lfloor\frac{n-1}{3}\right\rfloor,\qquad
\sigma(h)=\left\lceil\frac{n+1+3h}{2}\right\rceil.
\]

For a reflected zero \(n=qp+p-1-j\), its low-tail jump, when present, is

\[
h_{\mathrm{low}}
=p\left\lfloor\frac Hp\right\rfloor+j,
\]

while its high-tail jump is

\[
h_{\mathrm{high}}
=j+\left\lfloor\frac{(q-1)p+1}{3}\right\rfloor.
\]

If \(q=3s+1\), both values are \(sp+j\). Thus the same-index construction
aliases \(q=1,4,7,\ldots\). The exact counterexamples from \(17\mid A_3\)
are

\[
(q,n,h)=(4,81,20),\qquad(7,132,37).
\]

The script `q32_reflected_jump_alias.py` verifies every cutoff in these
examples and also checks the loss-count formula (4.4).
