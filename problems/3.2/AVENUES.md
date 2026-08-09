# Autonomous avenues: saturated corank and Cartier atoms

## Main goal

Prove the full Apéry GCD conjecture (G_n=e^{o(n)}), beginning by either
closing or rigorously refuting the live saturated quadruple-corank mechanism
and then attacking the pointwise atom tail with genuinely transverse data.

## Ranked avenues

### (a) Skipped-triple structural saturation of quadruple certificates

**Terminal outcome for the certificate architecture (2026-08-09): repaired.**
The former adjacent-pair mass

\[
 \sum_{a+b+c\le H}\log\gcd(S^*_{a,b},S^*_{b,c})
\]

contains a forced self-gcd whenever \(a=c\), because
\(|S_{a,b}|=|S_{b,a}|\).  The earlier decreasing-slope computation had
silently excluded this diagonal; an exact recomputation shows the total
mass at \(H=32\) is \(8.22348H^3\), with the diagonal alone equal to
\(0.252843H^4\).  Thus this is not the transverse mass required by the
corank argument.

Four zeros at gaps \((a,b,c)\) instead force the two resultants
\(S_{a,b}\) and \(S_{a,b+c}\).  They share the same starting roots, so
the Sylvester-corank valuation certificate is unchanged, while \(c\ge2\)
prevents a self-gcd even after swapping the resultant parameters.  The new
live hypothesis is

\[
 \mathrm{GM}^{\dagger}:\qquad
 \sum_{a+b+c\le H}\log\gcd(S^*_{a,b},S^*_{a,b+c})
 \ll H^{3+o(1)}.
\]

The exact ratio of this mass to \(H^3\) is \(0.123287\) at \(H=20\) and
\(0.168730\) at \(H=32\); the mean logarithmic gcd is respectively
\(1.450436\) and \(1.513120\).  This is evidence, not a proof.  The full
conditional implication has been repaired around
\(\mathrm{GM}^{\dagger}\); proving this new mass bound remains the live
arithmetic problem.

Replace the prime-support saturation

\[
S^*_{h,k}=S_{h,k}/\gcd(S_{h,k},\mathfrak U_H^\infty)
\]

by an exact valuation-level factorization.  Determine whether four actual
Apéry zeros force a (p)-adic valuation strictly above the predictable
factorial/Apéry/Lucas content in both adjacent resultants.  This would retain
the detecting prime after structural factors are divided out.

- Success: prove a primitive quadruple certificate whose total logarithmic
  gcd mass is (H^{3+o(1)}), including primes dividing the structural carrier,
  and derive the full (X^{3/2+o(1)}) dyadic zero-mass bound with every step
  written in `atom_tail_section.tex`.
- Proof-of-failure: exhibit an infinite or parametrically growing family in
  which an actual four-zero witness contributes no excess valuation beyond
  structural content, or prove that the primitive valuation statement is
  equivalent to an already open atom/zero-mass bound.

### (b) Weighted exceptional-class absorption

**Terminal outcome (2026-08-09): success by level-adaptive absorption.**
The fixed-scale estimate displayed below is not known.  It is also not needed:
on the dyadic level $T<Z(p)\le 2T$, take $H=\lceil16X/T\rceil$.
Then the exceptional-prime count gives level mass

\[
 2T\cdot O(H^2/\log X)=O(X^2/(T\log X)),
\]

whose sum over $T=2^k\sqrt X$ is $O(X^{3/2}/\log X)$.  The
generic part is $O(H^{3+o(1)}/\log X)$ under
$\mathrm{GM}^{\dagger}$, and its level sum is $X^{3/2+o(1)}$.
This closes the exceptional-mass gap in the conditional implication; the
complete proof is now in `atom_tail_section.tex`.

Attack directly the missing estimate

\[
\sum_{\substack{X<p\le2X\\p\mid\mathfrak U_H}} Z(p)
 \ll X^{3/2+o(1)}
\qquad(H\asymp X^{1/2}),
\]

splitting factorial, Apéry-value, and Lucas rank-of-apparition factors.  A
mere bound on the number of exceptional primes is insufficient because the
summand is $Z(p)$.

- Success: prove the displayed weighted estimate without assuming the desired
  average-zero theorem, and combine it with reduced gcd mass.
- Proof-of-failure: give a quantifier-preserving reduction showing that the
  weighted estimate already contains the original average-zero or atom-tail
  problem.

### (c) Unbounded Cartier--Frobenius compression

Search the full block jet, not another fixed (p)-adic digit, for a determinant
or Casoratian eliminating the one new local coordinate introduced at each
grade.  The output must be a nonzero global carrier whose logarithmic height
per recovered grade tends to zero.

- Success: construct an all-order recurrence/determinant yielding either the
  atom bound or an apparition law for pairs ((p,m\bmod p)).
- Proof-of-failure: prove a rank/dimension conservation theorem for the entire
  Cartier jet, not just the already checked orders through (p^8).

### (d) Singular-digit and diagonal atom decomposition

Separate (q=0), regular (p\nmid b_q), singular (p\mid b_q, q\ne r),
and repeated-digit (q=r) channels in (m=qp+r).  Exploit the exact divisor
condition (m=q(p+1)) on the diagonal and test whether the remaining singular
branch has a low-height cross-prime carrier beyond pair certificates.

- Success: obtain a pointwise (K_X(m)\le\lambda_X X^{o(1)}) bound or a strict
  power improvement on the unresolved regular/singular mass.
- Proof-of-failure: reduce every branch explicitly to CRT-compatible local
  residues or an already documented pair-level certificate.

## Fallback after all four terminal verdicts

Keep the conditional ($\mathrm{GM}^{\dagger}$) statement only with a complete proof of
its implication, correct any unconditional wording, and record the exact
weighted-exceptional or unbounded-grade lemma that remains.  Do not replace the
full conjecture by the weaker (p^{1/2}) average-zero milestone.

## Existing terminal facts that constrain every avenue

- Fixed Cartier/cofactor depth is already computed through precision (p^8)
  and remains CRT reconstruction.
- The fourth- and fifth-digit companion laws are rank one and CRT-compatible.
- Pair certificates, reflection, Mellin twists, fixed-depth height arguments,
  and unsaturated resultant mass do not break certificate parity.
- The current paper's exceptional-class lemma counts primes but does not by
  itself bound their (Z(p))-weighted mass.
