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

Four zeros at gaps \((a,b,c)\) also force the two resultants
\(S_{a,b}\) and \(S_{a,b+c}\).  This removes the literal self-gcd, but
does not yet remove algebraic dependence.  On the slice \(c=a\), every
even \(b\) contributes the propagated center factor
\(T_a^{(b)}\) to both resultants.  Its exact mass is \(H^4\)-scale:
at \(H=32\), the \(c=a\) slice alone has mass
\(0.00483348H^4\).  Thus the intermediate
\(\mathrm{GM}^{\dagger}\) formulation is also withdrawn.

For even \(b\), divide the central linear factor from the polynomial,
\[
 N_b^\circ(x)=N_b(x)/(2x+b+1),\qquad
 D_{a,b}=\operatorname{Res}(N_a(x),N_b^\circ(x+a));
\]
for odd \(b\), put \(N_b^\circ=N_b\).  A palindromic four-zero
configuration not detected by \(D_{a,b}\) has endpoints summing to
\(p-1\); among globally consecutive zero quadruples there is at most one
such centered configuration per prime.  The remaining slice
\(a=b=c\) has only one parameter and is controlled unconditionally by
\(\sum_{3a\le H}\log|S_{a,a}|\ll H^3\log H\).

The scalar repair still charges a prime that divides its two resultants at
unrelated roots.  The exact same-root invariant is instead

\[
 \mathcal C_{a,b,c}=\operatorname{cont}_T
 \operatorname{Res}_x(F_{a,b,c},G_{a,b,c}+T J_{a,b,c}),
\]

where off the palindromic slice
\((F,G,J)=(N_a,N_b(\cdot+a),N_{b+c}(\cdot+a))\), while on
\(c=a\ne b\) the middle polynomial is
\(N_b^\circ(\cdot+a)\).  A Gauss-DVR Sylvester argument proves that
\(t\) distinct common roots force \(v_p(\mathcal C)\ge t\).  Moreover
\(\mathcal C^*\) divides the corresponding saturated scalar gcd, because
the constant and top coefficients of the pencil resultant are the two old
resultants (up to a leading-coefficient power supported on
\(\mathfrak U_H\)).

The new live hypothesis is therefore the weaker, correctly aligned mass

\[
 \mathrm{AC}^{\mathrm{tr}}:\quad
 \sum_{\substack{a,b,c\ge2,\ a+b+c\le H\\\neg(a=b=c)}}
 \log \mathcal C^*_{a,b,c}\ll H^{3+o(1)}.
\]

Exact \(\mathbb Z[T]\) computation gives zero reduced mass through
\(H=14\).  At \(H=20\), only 2 of 675 terms are nontrivial and the ratio
of the mass to \(H^3\) is \(0.001620394\); at \(H=24\), only 18 of 1323
terms are nontrivial and the ratio is \(0.007248253\).  This is evidence, not a
proof.  The full conditional implication is now written around
\(\mathrm{AC}^{\mathrm{tr}}\); proving this content bound remains the live
arithmetic problem.

- Success: prove the aligned-content mass is \(H^{3+o(1)}\), or dominate it
  by a global carrier of that height.  The level-adaptive proof then gives
  the full \(X^{3/2+o(1)}\) dyadic zero-mass bound.
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
$\mathrm{AC}^{\mathrm{tr}}$, and its level sum is $X^{3/2+o(1)}$.
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

**Fixed-depth subavenue: terminal failure (Q7104, 2026-08-09).**  Lucas
gives the exact Boolean-OR law

\[
 b_{\sum_i n_ip^i}\equiv0\pmod p
 \quad\Longleftrightarrow\quad
 \text{some }n_i\in Z(p),
\]

and hence zero density \(1-(1-Z(p)/p)^d\) at depth~\(d\).  Singular
blocks are entirely zero and regular blocks merely copy the base atom.
Thus every fixed-depth Cartier tower describes descendants of the atom but
places no restriction on the base set \(Z(p)\).  Only an unbounded-order
compression could still be new.

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

Keep the conditional ($\mathrm{AC}^{\mathrm{tr}}$) statement only with a complete proof of
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
