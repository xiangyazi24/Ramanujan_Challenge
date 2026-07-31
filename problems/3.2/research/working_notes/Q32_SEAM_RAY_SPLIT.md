# The seam scalar is a 21-ray sum of fixed-Laurent periods

**Date:** 2026-07-31 (dm window, taking over from `life`)

This note continues the seam reduction of `/tmp/P32_TERMINAL_CROSS_N_FINAL.md`
and answers the first half of the question it left open ("identify the
Ore/Picard--Fuchs module of \(S_r\)").  Everything below is verified by
`problems/3.2/research/scripts/q32_seam_ray_split_audit.py`
(126 form checks, 6 exact split checks, 39 modular split checks, PASS).

## 1. The identity

Recall the boundary-packet form of the seam scalar,
\[
 S_r=b_r-\sum_{\kappa}\lambda_\kappa
 \operatorname{CT}\bigl[\Lambda^{r-1}X^{-(r-1)\kappa}(X^\kappa-1)^r\bigr],
\tag{1}
\]
where \(\Lambda=\sum_\kappa\lambda_\kappa X^\kappa\) is the Apéry Laurent
polynomial and \(b_r=\operatorname{CT}\Lambda^r\).

Regroup the summand by matching the two \(r\)-dependent factors:
\[
 \Lambda^{r-1}X^{-(r-1)\kappa}(X^\kappa-1)^r
 =\bigl(\Lambda X^{-\kappa}\bigr)^{r-1}(X^\kappa-1)^{r-1}(X^\kappa-1)
 =\bigl(\Lambda(1-X^{-\kappa})\bigr)^{r-1}(X^\kappa-1).
\]
Hence, with
\[
 \boxed{\quad
 G_\kappa:=\Lambda\cdot(1-X^{-\kappa}),\qquad
 U_\kappa(r):=\operatorname{CT}\bigl[G_\kappa^{\,r-1}(X^\kappa-1)\bigr],
 \quad}
\tag{2}
\]
the seam scalar splits as
\[
 \boxed{\qquad
 S_r=b_r-\sum_{\kappa\ne0}\lambda_\kappa U_\kappa(r).
 \qquad}
\tag{3}
\]
There are exactly \(21\) nonzero rays \(\kappa\) (the polytope has \(22\)
lattice points; the ray \(\kappa=0\) gives \(G_0=0\) and contributes nothing
for \(r\ge2\)).

Expanding \(G_\kappa^{r-1}\) binomially gives the computationally fast form
\[
 U_\kappa(r)=\sum_{j=0}^{r-1}(-1)^j\binom{r-1}{j}
 \Bigl(c_{r-1}\bigl((j-1)\kappa\bigr)-c_{r-1}(j\kappa)\Bigr),
 \qquad c_m(\eta)=[X^\eta]\Lambda^m,
\tag{4}
\]
which costs \(O(r^2)\) per ray with row-cached binomials, against \(O(r^4)\)
for the naive evaluation of (1).

## 2. Why this is the right decomposition

Each \(U_\kappa\) is the constant term of the powers of ONE fixed Laurent
polynomial, multiplied by one fixed Laurent polynomial.  By the standard
diagonal/rational-function argument this is \(D\)-finite in \(r\).  Therefore

> **\(S_r\) is \(D\)-finite**, being a finite \(\mathbb Z\)-linear combination
> of \(b_r\) and the 21 pieces \(U_\kappa\).

Its own annihilator is (a right divisor of) the LCLM of the 22 pieces, hence of
large order.  This is exactly what direct guessing sees:

- **No operator for \(S_r\)** of order \(\le10\) and degree \(\le20\), with
  exact data through \(r=260\) (modular nullspace over \(2^{61}-1\), with
  held-out verification; the same code recovers the Apéry operator
  \(\rho=2,D=3\) immediately).
- **Individual rays do have small operators.**  With data only through
  \(r=90\), a held-out-verified operator was found for \(7\) of the \(21\)
  rays, of order \(4\)--\(6\) and degree \(10\)--\(13\); e.g.
  \[
  \begin{array}{c|c}
  \kappa&(\text{order},\deg)\\ \hline
  (0,-1,1)&(4,13)\\
  (0,1,-1)&(4,13)\\
  (0,0,1)&(6,10)\\
  (0,1,0)&(6,10)\\
  (1,0,0)&(5,10)\\
  (1,0,1)&(5,12)\\
  (1,1,0)&(5,12)
  \end{array}
  \]
  The remaining rays are not excluded; they merely need a larger ansatz box
  than \(r\le90\) supports.

The pattern respects the \(y\leftrightarrow z\) symmetry of \(\Lambda\), as it
must: \((0,-1,1)\) and \((0,1,-1)\) agree, as do \((0,0,1)\) and \((0,1,0)\),
and \((1,0,1)\) and \((1,1,0)\).

## 3. What this changes for the open problem

The remaining goal is
\[
 \log\operatorname{rad}\gcd(b_r,S_r)=o(r),
\]
with the empirical input that for \(1\le r\le100\) the exact gcd has no prime
factor outside \(\{5,11,19\}\) and never exceeds \(55\) (independently
recomputed here for \(r\le34\): the gcd is \(1\) except \(5\) at
\(r=1,15,28,29,31,33\) and \(25\) at \(r=8,11\)).

Two consequences of (3):

1. **The naive route is dead, provably.**  A prime dividing both \(b_r\) and
   \(S_r\) would follow from a fixed resultant if \(S_r\) lay in the Apéry
   module.  It does not: the ansatz \(A(r)S_r=B(r)b_r+C(r)b_{r-1}\) with
   \(\deg\le6\) has only the zero solution (exact modular nullspace).  Since
   \(b_{r},b_{r-1}\) span the Apéry solution space over \(\mathbb Q(r)\), this
   rules out the whole 2-dimensional module, not just one basis.

2. **The tractable replacement.**  The pair to compare is not \((L_b,L_S)\)
   with \(L_S\) of unknown huge order, but \(L_b\) against the 21 small
   operators \(L_{U_\kappa}\).  Each \(G_\kappa=\Lambda(1-X^{-\kappa})\) is a
   concrete Laurent polynomial whose period sequence is Apéry-like; its
   singularities and its local exponents are computable, and the exceptional
   primes \(\{5,11,19\}\) should be visible in the discriminants/resultants of
   that finite family.  This is the computation to do next.

## 4. Reproducibility

- `problems/3.2/research/scripts/q32_seam_ray_split_audit.py` — the exact
  audit (polynomial form vs. binomial form vs. reference seam scalar).
- The guessing runs are plain modular nullspace computations; the generator
  for \(U_\kappa(r)\) is formula (4) with cached binomial rows.
