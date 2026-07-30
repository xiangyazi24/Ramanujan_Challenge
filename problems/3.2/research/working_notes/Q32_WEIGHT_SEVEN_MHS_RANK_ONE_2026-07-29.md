# The effective weight-seven endpoint MHS quotient is one-dimensional

Date: 2026-07-29. Owner: Codex.

## 0. Result

For every prime \(p\ge 11\), use strict finite multiple harmonic sums
\[
 H(a_1,\ldots,a_d)
 =\sum_{1\le k_1<\cdots<k_d\le p-1}
   \frac1{k_1^{a_1}\cdots k_d^{a_d}},
\]
and put
\[
 \xi=\frac{H(6)}p,\qquad
 \eta=\frac{H(2,4)}p,\qquad
 A=H(2,2,3),\qquad B=H(2,5)
 \quad\pmod p.
\]
The divisions by \(p\) are integral.  The four coordinates satisfy
\[
 \boxed{
  3\eta=2\xi,\qquad 3A=14\xi,\qquad 2B=-7\xi
  \pmod p.}
\tag{0.1}
\]
Thus the conservative four-generator presentation
\[
 \langle \xi,\eta,A,B\rangle
\]
which occurs in the effective weight-seven endpoint calculation has
actual image of dimension at most one:
\[
 (\xi,\eta,A,B)
 =\xi\left(1,\frac23,\frac{14}3,-\frac72\right).
\tag{0.2}
\]

The point missed by ordinary reversal is the first lifted correction
\[
 H(2,4)-H(4,2)
 \equiv p\{2H(4,3)+4H(5,2)\}\pmod {p^2}.
\tag{0.3}
\]
Keeping (0.3), rather than dividing a congruence known only modulo
\(p\), proves the first relation in (0.1).

## 1. The depth-two weight-seven values

Let
\[
 \beta=B_{p-7}\pmod p
\]
with the convention \(B_1=-1/2\).  For \(a+b=w<p\), the standard
depth-two congruence is
\[
 H(a,b)
 \equiv
 \frac{(-1)^b}{w}\binom wa B_{p-w}\pmod p.
\tag{1.1}
\]
For completeness, write
\[
 H(a,b)
 =\sum_{j=1}^{p-1}j^{-b}
   \sum_{i=1}^{j-1}i^{-a}.
\]
Replace \(i^{-a}\) and \(j^{-b}\) modulo \(p\) by
\(i^{p-1-a}\) and \(j^{p-1-b}\), and apply Faulhaber's formula to
the inner sum.  In the resulting Bernoulli-polynomial expansion, the
only power whose sum over \(\mathbf F_p^\times\) is nonzero has index
\(k=p-w\).  Its coefficient is
\[
 -\frac1{p-a}\binom{p-a}{p-w}
 \equiv
 \frac1a\binom{-a}{b}
 =\frac{(-1)^b}{w}\binom wa\pmod p,
\]
which proves (1.1).

At \(w=7\), (1.1) gives
\[
 H(2,5)=-3\beta,\qquad
 H(4,3)=-5\beta,\qquad
 H(5,2)=3\beta
 \pmod p.
\tag{1.2}
\]
The standard even full-sum congruence is
\[
 H(6)\equiv \frac67pB_{p-7}\pmod {p^2}.
\tag{1.3}
\]
One derivation of (1.3) replaces \(k^{-6}\) modulo \(p^2\) by
\(k^{p(p-1)-6}\), applies Faulhaber's formula modulo \(p^2\), and
then uses Kummer's congruence between the even indices
\(p(p-1)-6\) and \(p-7\).  Consequently
\[
 \xi=\frac67\beta.
\tag{1.4}
\]
Equations (1.2) and (1.4) already prove
\[
 2B=-7\xi.
\tag{1.5}
\]

## 2. Lifted reversal determines \(\eta\)

Changing variables
\[
 (i,j)\longmapsto(p-j,p-i)
\]
and expanding both factors to first order gives, when \(a+b\) is even,
\[
 H(a,b)
 \equiv H(b,a)
 +p\{aH(b,a+1)+bH(b+1,a)\}\pmod {p^2}.
\tag{2.1}
\]
The case \((a,b)=(2,4)\), together with (1.2), is
\[
 H(2,4)-H(4,2)
 \equiv 2p\beta\pmod {p^2}.
\tag{2.2}
\]

Stuffle gives the exact identity
\[
 H(2)H(4)=H(2,4)+H(4,2)+H(6).
\tag{2.3}
\]
Both full sums on the left are divisible by \(p\), while \(H(6)\) is
also divisible by \(p\).  Ordinary reversal modulo \(p\), followed by
(2.3), first shows that both \(H(2,4)\) and \(H(4,2)\) are divisible
by \(p\).  Dividing (2.2)--(2.3) by \(p\) is now legitimate and gives
\[
 \begin{aligned}
  \frac{H(2,4)}p-\frac{H(4,2)}p&=2\beta,\\
  \frac{H(2,4)}p+\frac{H(4,2)}p&=-\xi
  \pmod p.
 \end{aligned}
\]
Hence
\[
 \eta=\beta-\frac12\xi=\frac47\beta=\frac23\xi,
\tag{2.4}
\]
which is the first relation of (0.1).

## 3. The depth-three coordinate

Odd-weight reversal gives
\[
 H(2,3,2)=0\pmod p.
\tag{3.1}
\]
The stuffle product with a repeated first exponent is
\[
 H(2)H(2,3)
 =2H(2,2,3)+H(2,3,2)+H(4,3)+H(2,5).
\tag{3.2}
\]
Since \(H(2)=0\pmod p\), equations (1.2), (3.1), and (3.2) yield
\[
 2A-5\beta-3\beta=0,
\qquad A=4\beta=\frac{14}3\xi.
\tag{3.3}
\]
This completes the proof of (0.1).

## 4. Scope

This removes a possible new finite-MHS dimension at effective weight
seven.  It does **not** by itself prove the all-\(m\) precision-eight
endpoint law.  That law is written after lower endpoint coordinates
and their lifted digits have been subtracted.  A proof still has to
perform that full coordinate change and telescope the resulting
direct/reflected block sums.  Naively projecting only the primitive
four-vector loses those lower-coordinate lift terms.

The exact audit is

```text
../scripts/q32_weight_seven_mhs_rank_one_audit.py
```

Through all \(164\) primes \(11\le p\le1000\), it checks the lifted
reversal, stuffle identities, Bernoulli normalizations, and all three
relations in (0.1), with no failure.
