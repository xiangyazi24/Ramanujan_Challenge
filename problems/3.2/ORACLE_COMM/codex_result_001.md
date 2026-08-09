# Codex Result 001: the exact bilinear CRT target

## 1. Exact dispersion expansion

Fix \(I_N=(N,2N]\) and a dyadic prime block

\[
\mathcal P(P)=\{p\text{ prime}:P<p\le2P\}.
\]

For \(p\in\mathcal P(P)\), put

\[
X_p(m)=\mathbf1_{\{m\bmod p\in\mathcal Z_p\}},
\qquad
A_p=\sum_{m\in I_N}X_p(m),
\]

\[
J_{p,q}
=\sum_{m\in I_N}X_p(m)X_q(m),
\qquad
B_P(m)=\sum_{p\in\mathcal P(P)}X_p(m),
\]

and

\[
S=S(P,N)=\sum_pA_p.
\]

All sums over \(p\ne q\) below are ordered. Expanding around the
empirical mean gives

\[
\boxed{
V^\circ(P,N)=D_{P,N}+E^\circ_{P,N},
}
\]

where

\[
D_{P,N}
=\sum_pA_p\left(1-\frac{A_p}{N}\right)
\le S
\]

and the exact cross-prime bilinear sum is

\[
\boxed{
E^\circ_{P,N}
=\sum_{p\ne q}
\left(J_{p,q}-\frac{A_pA_q}{N}\right)
=\sum_{p\ne q}\sum_{m\in I_N}
\left(X_p(m)-\frac{A_p}{N}\right)
\left(X_q(m)-\frac{A_q}{N}\right).
}
\tag{1}
\]

Thus the required one-sided estimate is

\[
E^\circ_{P,N}\ll N^{o(1)}S.
\tag{2}
\]

No pairwise absolute bound and no universal negative-covariance
statement is needed.

Equivalently, define the ordered collision count

\[
F^{\mathrm{ord}}_{P,N}
=\sum_{p\ne q}J_{p,q}
=\sum_{m\in I_N}B_P(m)(B_P(m)-1).
\]

Then

\[
\boxed{
V^\circ(P,N)
=S+F^{\mathrm{ord}}_{P,N}-\frac{S^2}{N}.
}
\tag{3}
\]

Consequently \(V^\circ\le CS\) is exactly equivalent to

\[
\boxed{
F^{\mathrm{ord}}_{P,N}
\le\frac{S^2}{N}+(C-1)S.
}
\tag{4}
\]

This is the block second-factorial-moment, or projected-fiber-square,
form of AMTD.

For logarithmic weights \(w_p=\log p\), the exact off-diagonal term is

\[
E^{\circ,w}_{P,N}
=\sum_{p\ne q}w_pw_q
\left(J_{p,q}-\frac{A_pA_q}{N}\right),
\tag{5}
\]

and the diagonal is

\[
D^w_{P,N}
=\sum_pw_p^2A_p(1-A_p/N)
\le\sum_pw_p^2A_p.
\]

Inside one dyadic block, \(w_p\asymp\log P\), so this is the same
collision problem up to dyadic weight constants.

### Honest shell Fourier form

Define

\[
g_p(m)=X_p(m)-A_p/N,
\qquad
\Phi_p(\theta)=\sum_{m\in I_N}g_p(m)e(m\theta).
\]

Continuous Parseval gives

\[
V^\circ(P,N)
=\int_{\mathbb T}
\left|\sum_p\Phi_p(\theta)\right|^2\,d\theta
\]

and hence

\[
\boxed{
E^\circ_{P,N}
=\sum_{p\ne q}
\int_{\mathbb T}
\Phi_p(\theta)\overline{\Phi_q(\theta)}\,d\theta.
}
\tag{6}
\]

The uncentered part of \(\Phi_p\) is exactly

\[
\sum_{r\in\mathcal Z_p}
\sum_{\substack{j\in\mathbb Z\\N<r+jp\le2N}}
e\!\left(\theta(r+jp)\right).
\]

This includes every admissible lift and automatically handles the
top-block clipping omitted by the unrestricted two-lift formula in the
palindromic-Fourier remark of proof.tex.

## 2. Exact short-arc/Linnik bilinear form

Let

\[
F_p(a)=\sum_{r\in\mathcal Z_p}e_p(ar),
\qquad
c_{p,a}=\frac{(\log p)F_p(a)}p,
\]

\[
Q_P
=\sum_{P<p\le2P}\sum_{a=1}^{p-1}|c_{p,a}|^2
=\sum_{P<p\le2P}(\log p)^2
\left(\frac{Z(p)}p-\frac{Z(p)^2}{p^2}\right),
\]

and

\[
K_M(t)=\left(M^{-1}-\|t\|\right)_+.
\]

After the same-prime terms are separated, the exact unequal-prime
short-arc form from oracleA_result.tex is

\[
\boxed{
R_{\ne}(P,M)
=\sum_{\substack{P<p,q\le2P\\p\ne q}}
\sum_{a=1}^{p-1}\sum_{b=1}^{q-1}
K_M\!\left(\frac ap-\frac bq\right)
\frac{(\log p)(\log q)}{pq}
F_p(a)\overline{F_q(b)}.
}
\tag{7}
\]

The precise missing scalar estimate is one-sided:

\[
\boxed{
\operatorname{Re}R_{\ne}(P,M)
\ll M^{-1+o(1)}Q_P.
}
\tag{SDC}
\]

Let \(\langle d\rangle_{pq}\) denote the least-absolute-value
representative modulo \(pq\). Setting

\[
k=\langle aq-bp\rangle_{pq}
\]

gives the exact reciprocal form

\[
R_{\ne}(P,M)=M^{-1}\mathcal B_{P,M},
\]

\[
\boxed{
\mathcal B_{P,M}
=\sum_{\substack{P<p,q\le2P\\p\ne q}}
\frac{(\log p)(\log q)}{pq}
\sum_{\substack{0<|k|<pq/M\\(k,pq)=1}}
\left(1-\frac{M|k|}{pq}\right)
F_p(k\bar q)F_q(k\bar p).
}
\tag{8}
\]

Here \(\bar q\) and \(\bar p\) are inverses modulo \(p\) and \(q\).
The target is

\[
\boxed{
\operatorname{Re}\mathcal B_{P,M}
\ll M^{o(1)}Q_P.
}
\tag{9}
\]

This is the moving reciprocal-fraction bilinear sum. The arguments
\(k\bar q\bmod p\), \(k\bar p\bmod q\), and both zero sets vary jointly
with \(p,q\). Separate bounds for individual primes do not address this
coupling.

There is an important distinction between (1) and (7). Formula (1)
uses empirical shell means \(A_p/N\) and is literally equivalent to
AMTD. The short-arc calculation uses the vertical means \(Z(p)/p\).
Gallagher's lemma with the necessary dilation \(M=4N\) gives

\[
\sum_{N<m\le2N}|D_P(m)|^2
\ll N^{1+o(1)}Q_P,
\]

where

\[
D_P(m)=\sum_p(\log p)
\left(X_p(m)-\frac{Z(p)}p\right).
\tag{10}
\]

Empirical centering can only decrease this energy. The weighted
majorant \(NQ_P\) is sufficient for the all-index conclusion, but is
not literally \(S(P,N)\), particularly when \(p>N\) and shell clipping
removes some residues. Thus (SDC) is the exact missing short-arc lemma,
while (1) or (4) is the exact literal AMTD target.

## 3. What the gap-polynomial projection gives

For an ordinary non-wrapping gap \(h\),

\[
b_m\equiv b_{m+h}\equiv0\pmod p
\quad\Longrightarrow\quad
N_h(m)\equiv0\pmod p,
\tag{11}
\]

where \(\deg_mN_h=3(h-1)\). This supplies the vertical estimate

\[
\#\{m:b_m\equiv b_{m+h}\equiv0\pmod p\}
\le3(h-1)
\]

whenever the reduction of \(N_h\) is nonzero.

Three limitations prevent this from implying (2) or (SDC).

1. Equation (11) is only a necessary condition unless the condition
   \(b_m=0\) is retained. A root of \(N_h\) is a candidate return of
   the projective recurrence, not automatically an Apéry zero pair.
2. The degree is bounded only for fixed \(h\). In the moving range
   \(h\asymp p\), both degree and conductor grow with \(p\).
3. For a zero-dimensional root set, the bound by \(\deg N_h\) is
   already elementary. A square-root Weil error becomes relevant only
   after constructing a positive-dimensional trace-function family;
   its conductor here grows with \(h\). Summing such pointwise errors
   absolutely over \(h,p,q\) is too costly.

The projection becomes useful for dispersion only through its fiber
square. Define the actual hit incidence

\[
\mathscr X
=\{(m,p):m\in I_N,\ p\in\mathcal P(P),\
m\bmod p\in\mathcal Z_p\},
\qquad
\pi(m,p)=m.
\]

Then

\[
B_P(m)=|\pi^{-1}(m)|
\]

and

\[
\boxed{
F^{\mathrm{ord}}_{P,N}
=\#\left(
(\mathscr X\times_{I_N}\mathscr X)_{p\ne q}
\right).
}
\tag{12}
\]

Thus projection helps only if one can show that this off-diagonal fiber
square has no excess component and that its aggregate point-count error
is \(O(N^{o(1)}S)\). Fixed-\(h\) degree bounds give vertical information
inside one characteristic; (12) asks for horizontal information across
two varying characteristics.

A sheaf-theoretic input strong enough to address (8) would need all of
the following:

- the centered indicator
  \(\mathbf1_{\mathcal Z_p}-Z(p)/p\) arises from a compatible
  bounded-conductor pushforward family;
- after the reciprocal substitutions in (8), unequal-prime tensor
  products have no geometrically constant or exceptional rank-one
  constituent;
- the joint horizontal family satisfies a bilinear large-sieve bound in
  \(p,q\), rather than merely a separate Weil bound for each pair.

The last item is indispensable. On the doublet sector,
\(a\mapsto F_p(a)\) is already a rank-two bounded-conductor trace
function, yet the reflected anchored model described below still
violates (SDC). Local bounded conductor alone is therefore insufficient.

## 4. Concrete sufficient incidence lemmas

The exact fiber-square criterion is:

> **No-excess fiber-square lemma.** Suppose that, uniformly in \(N,P\),
> \[
> F^{\mathrm{ord}}_{P,N}
> \le
> \frac{S^2-\sum_pA_p^2}{N}
> +C_0N^{o(1)}S.
> \tag{PF2}
> \]
> Then
> \[
> V^\circ(P,N)
> \le(1+C_0N^{o(1)})S(P,N).
> \]

Indeed,

\[
E^\circ_{P,N}
=F^{\mathrm{ord}}_{P,N}
-\frac{S^2-\sum_pA_p^2}{N}
\le C_0N^{o(1)}S,
\]

and \(D_{P,N}\le S\).

A stronger, simpler algebraic-incidence condition is uniformly bounded
degree of the horizontal hit projection:

\[
\max_{m\in I_N}|\pi^{-1}(m)|
=\max_mB_P(m)\le d.
\tag{13}
\]

It immediately gives

\[
V^\circ
\le\sum_mB_P(m)^2
\le d\sum_mB_P(m)
=dS.
\]

Condition (13) is generally too strong; (PF2) is the sharp useful
version. It permits rare large fibers but requires their second
factorial moment to remain at the independent main-term scale.

## 5. Palindromic Fourier structure

With

\[
F_p(a)=\sum_{r\in\mathcal Z_p}e_p(ar),
\]

the reflection \(r\leftrightarrow p-1-r\) gives

\[
F_p(a)
=e\!\left(\frac{a(p-1)}{2p}\right)R_p(a),
\qquad R_p(a)\in\mathbb R.
\tag{14}
\]

Therefore the ordinary real rotation is

\[
e\!\left(-\frac{a(p-1)}{2p}\right)F_p(a)\in\mathbb R.
\tag{15}
\]

Equivalently,

\[
e_p(a/2)F_p(a)\in\mathbb R,
\]

where \(1/2=(p+1)/2\) in \(\mathbb F_p\). The statement using the
positive ordinary phase
\(e(a(p-1)/(2p))F_p(a)\) is generally false: that positive factor is
the phase already present in \(F_p(a)\).

For a doublet

\[
\mathcal Z_p=\{r_p,p-1-r_p\},
\qquad
h_p=p-1-2r_p,
\]

formula (14) becomes

\[
F_p(a)
=2e\!\left(\frac{a(p-1)}{2p}\right)
\cos\!\left(\frac{\pi ah_p}{p}\right).
\tag{16}
\]

This factorization does not force cancellation in (7). If
\(\delta=a/p-b/q\), then

\[
F_p(a)\overline{F_q(b)}
=(-1)^{a-b}e(-\delta/2)R_p(a)R_q(b).
\tag{17}
\]

On a short arc, the extra phase is nearly constant after choosing the
appropriate cyclic representative. Palindromy converts the amplitudes
to real cosines, but supplies no sign cancellation as \(p,q,k\) vary.
After conjugate pairing it mainly makes the bilinear form real.

The rigorous obstruction is the reflected anchored family

\[
\mathcal Z_p^*
=\{\,2N\bmod p,\,-2N-1\bmod p\,\},
\qquad
N<p\le2N,
\]

with the finitely many primes dividing
\((4N+1)(2N+1)\) removed. Every set is a nonadjacent reflected doublet
with bounded local Fourier conductor, but every column contains the
same shell point \(2N\). Hence

\[
B_P(2N)\asymp\#\mathcal P(P),
\]

while \(S\asymp\#\mathcal P(P)\), so

\[
V^\circ\gg(\#\mathcal P(P))^2.
\]

The same example violates (SDC) by \(N^{1-o(1)}\). Palindromy,
doublet cardinality, vertical degree bounds, and bounded local conductor
therefore do not imply the required bilinear cancellation.

## Conclusion

The exact missing object for AMTD is the signed ordered covariance sum
(1), equivalently the collision estimate (4). In the short-arc route it
is the moving reciprocal sum (8) with the one-sided bound (9). Gap
polynomials provide fixed-prime, fixed-gap necessary conditions.
Projecting their incidence helps only after controlling the
off-diagonal fiber square of the actual cross-prime hit projection.
The missing algebraic input must therefore be horizontal and
family-level—strong enough to imply (PF2) or (SDC)—rather than another
pointwise Weil bound for one fixed \(N_h\).
