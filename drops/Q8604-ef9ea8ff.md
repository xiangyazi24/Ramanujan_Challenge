ANSWER Q8604 ef9ea8ff

# Q8604 — Short-CRT gateway audit

## Verdict

The short-CRT formulation is a useful reparameterization of the full-defect FDTG3 constraint, not a strictly smaller theorem by itself. It removes the explicit shifted variable \(N,t\), but the missing cancellation remains the same three-prime Apéry correlation.

## 1. Exact equivalence and extra conditions

Fix ordered primes \(p_1<p_2<p_3\), signs \(s_i\), defects \(d_1=0,d_2,d_3\), and

\[
u_i(t)=s_i(t+d_i+\tau_i),\qquad \tau_i=(1-s_i)/2.
\]

The product variables satisfy

\[
n_1=c_1p_1,\qquad n_i=c_ip_i,
\]

and

\[
n_1-n_i=d_i.
\]

Hence

\[
c_1p_1-c_2p_2=d_2,\qquad c_1p_1-c_3p_3=d_3.
\]

Reducing modulo \(p_2,p_3\) gives the unique residue

\[
\kappa\equiv d_2p_1^{-1}\pmod {p_2},\qquad
\kappa\equiv d_3p_1^{-1}\pmod {p_3},
\]

with \(0\leq\kappa<p_2p_3\).

For an actual physical triple, \(c_1\in\mathbb Z\) and the target range gives

\[
1\le c_1\le X-1.
\]

Therefore

\[
c_1=\kappa.
\]

This implication is exact only after imposing the missing physical conditions:

1. \(p_1,p_2,p_3\) are distinct shell primes and are ordered.
2. The row is the common physical row
\[
 m=N+t,
\]
with \(N=c_1p_1\).
3. The remaining coefficients satisfy
\[
 c_i=(c_1p_1-d_i)/p_i
\]
and lie in the required positive/negative sign sectors.
4. The three folded roots are the actual Apéry folded zeros modulo their primes.
5. Quotient regularity and target cutoff \(p_i\le m\) hold.

Thus

\[
T_{\rm target}(E,s,d)\le SCRT_X(s,d)
\]

is valid, but equality is false because SCRT keeps triples satisfying only the first coefficient congruence test.

## 2. Additive-character expansion

The selector is

\[
1_{1\le \kappa<X}.
\]

Modulo \(M=p_2p_3\), write

\[
1_{1\le\kappa<X}=\frac1M\sum_{a\bmod M}
\left(\sum_{c=1}^{X-1}e(ac/M)\right)e(-a\kappa/M).
\]

Therefore

\[
SCRT=\sum_t\sum_{p_1<p_2<p_3}
\prod_i1_{p_i\in S_i(t)}
\frac1{p_2p_3}
\sum_{a\bmod p_2p_3}
\left(\sum_{c=1}^{X-1}e(ac/(p_2p_3))\right)e(-a\kappa/(p_2p_3)).
\]

The zero frequency is

\[
a=0:\qquad \frac{X-1}{p_2p_3}.
\]

Summing the zero mode gives

\[
\sum_t\sum_{p_1<p_2<p_3}
\frac{X}{p_2p_3}
1_{p_i\in S_i(t)}.
\]

This is not automatically \(X^2\lambda_X^3\); the missing factor is the size of the \(t\)-sum. The required estimate would need the Apéry-root average

\[
\sum_t\prod_i1_{p_i\in S_i(t)}
\ll X^2\prod_i\frac{|Z_{p_i}|}{p_i},
\]

which is precisely the three-character incidence statement.

## 3. Nonzero frequencies

The nonzero term is

\[
\sum_t\sum_{p_1<p_2<p_3}
\prod_i1_{p_i\in S_i(t)}
\sum_{a\ne0}
W_X(a;p_2p_3)e(-a\kappa/(p_2p_3)),
\]

where

\[
W_X(a;M)=\frac1M\sum_{c=1}^{X-1}e(ac/M).
\]

The phase contains

\[
\kappa\equiv d_2p_1^{-1}\pmod{p_2},
\quad
\kappa\equiv d_3p_1^{-1}\pmod{p_3}.
\]

CRT splitting gives two reciprocal phases, not a one-variable Kloosterman sum with independent coefficients. The remaining sums are still over simultaneous Apéry zero conditions:

\[
p_i\mid B_{u_i(t)}.
\]

Hence the exponential sum is a three-character Apéry correlation, not a separated product.

## 4. X=1024 packet check

For the packet

\[
m=15468,
\]

with

\[
p=1069,1381,1847,2011,
\]

one has

\[
m/p=(14,11,8,7)\text{ quotients}
\]

with the corresponding remainders giving the four shell hits. The short-CRT selector only chooses triples among these hits whose coefficient congruence has residue below \(X\). It does not remove the common-row coincidence.

The packet therefore demonstrates that the CRT selector is a filter, not an independent arithmetic constraint.

## 5. Conclusion

The exact statement is:

\[
T_{\rm target}\subseteq SCRT.
\]

But SCRT is not proven smaller than FDTG3. Its zero-frequency term already requires the same missing triple Apéry incidence estimate, and its nonzero frequencies remain the same correlated root problem.

A genuine saving would require a theorem of the form

\[
\sum_t\prod_{i=1}^3 1_{p_i\mid B_{u_i(t)}}
\ll
X^{2+o(1)}\prod_i\frac{|Z_{p_i}|}{p_i},
\]

or an equivalent cancellation statement for the nonzero-frequency sums. Without such an input, short-CRT is a reparameterization, not a new gate.
