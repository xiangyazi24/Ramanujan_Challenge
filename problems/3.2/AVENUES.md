# Autonomous avenues: saturated corank and Cartier atoms

## Main goal

Prove the full Apéry GCD conjecture (G_n=e^{o(n)}), beginning by either
closing or rigorously refuting the live saturated quadruple-corank mechanism
and then attacking the pointwise atom tail with genuinely transverse data.

## Breakthrough run begun 2026-08-09

The present run keeps the live target \(\mathrm{AC}^{\mathrm{tr}}\) fixed and ranks
four new attack avenues.  Claims inherited from earlier runs are inputs only
after exact algebraic or computational re-verification.

### (a0) Local-algebra formula for aligned content

**Exact formula proved; compression outcome negative (2026-08-09).**  Over a
complete DVR with perfect residue field and unit leading coefficient, factor
\(F=\prod_i f_i^{m_i}\), let \(L_i=K[x]/(f_i)\), and let \(\mathcal O_i\)
be its valuation ring.  Then
\[
 v\!\left(\operatorname{cont}_T\operatorname{Res}(F,G+TJ)\right)
 =\sum_i m_i\,\ell_R\bigl(\mathcal O_i/(G(\alpha_i),J(\alpha_i))\bigr).
\]
This is a length on the normalization, not generally on the order
\(R[x]/(F)\).  It identifies the local invariant exactly but supplies no
cross-triple cancellation: primitive off-center phantom roots already occur
at \((a,b,c,p,r)=(17,14,4,157,38)\), while \(Z(157)=\varnothing\).

Regard
\[
 \operatorname{cont}_T\operatorname{Res}_x(F,G+TJ)
\]
as a Fitting invariant of the finite \(\mathbf Z_p\)-algebra
\(\mathbf Z_p[x]/(F)\).  Derive an exact local formula, including inseparable
and non-monic cases, in terms of the simultaneous ideal \((G,J)\), and then
use the continuant identities to seek cancellation or compression when the
formula is summed over \((a,b,c)\).

- Success: prove the \(H^{3+o(1)}\) total aligned-content bound, or reduce it
  to a strictly lower-dimensional height sum of that size.
- Proof-of-failure: produce a verified family showing that local intersection
  lengths can have \(H^4\)-scale total mass even after the existing structural
  saturation, or prove that the Fitting reformulation is tautologically
  equivalent to the uncompressed content sum.

### (b0) Residual-prime rank amplification

**Terminal negative outcome (2026-08-09).**  The paper's \(\rho(p)\) is the
Lucas leading-coefficient rank, not an Ap\'ery zero rank.  Reflection only
adds the reversed collision chain and has order two; the exact residual
examples show no \(H\)-scale proliferation.  The attempted rank argument was
therefore based on a category error and cannot compress the mass.

For every residual prime found by the exact aligned pencils, compute its
rank of apparition and its entire collision graph up to several multiples of
that rank.  Test and then prove whether one aligned four-zero collision forces
many pair collisions, a large rank, or a long orbit under the Ap\'ery transfer
matrix.  Any amplification by a factor comparable to \(H\) would convert the
current \(H^4\)-parameter census into an \(H^3\)-scale carrier.

- Success: prove a uniform amplification/rank lemma strong enough to charge
  every residual valuation to \(O(H^{3+o(1)})\) orbit data.
- Proof-of-failure: give an infinite or parametrically growing family of
  isolated aligned collisions with no such orbit proliferation.

### (c0) Arithmetic-intersection height compression

**Generic route exhausted (2026-08-09).**  Prime-first Fitting/Smith
invariants can remove excess pencil valuation, but direct sums over the
\(\asymp H^3\) labelled triples still add their heights.  Product-formula,
Ekedahl-sieve, large-sieve, and block-diagonal determinant arguments do not
save a gap parameter without an Ap\'ery-specific cross-gap identity.  The
fully deflated adjacent pencil is retained as a cleaner alternative mass
hypothesis, not as an unconditional content-one theorem.

Combine the complex root-strip separation already used for the First Lemma
with the non-archimedean intersection lengths from (a0).  Look for an
adelic/product-formula estimate, arithmetic B\'ezout inequality, or determinant
normalization that bounds the total finite-place overlap without paying for
the full triangular product of cut-edge values.

- Success: obtain a global height carrier with logarithmic height
  \(H^{3+o(1)}\) that dominates all aligned local lengths.
- Proof-of-failure: prove that every such direct product-formula carrier must
  retain \(\Omega(H^3\log H)\) height and cannot be absorbed at the critical
  dyadic level.

### (d0) Transfer-matrix codegree theorem

**Generic codegree route exhausted; Ap\'ery-specific gcd tail remains open
(2026-08-09).**  The affine-plane construction in
`gcdtail_result.tex` shows that the available row-degree and column bounds,
even strengthened to codegree one, still permit \(\asymp H^2\) deep pairs.
The scalar resultant-gcd formulation also loses the common-root label, while
the aligned second moment misses singleton fibers.  Thus a uniform codegree
bound by itself cannot supply the required saving.  What remains meaningful
is the strictly stronger, arithmetic split affine gcd-tail estimate displayed
in that file, together with a separate estimate for the complementary
annulus.

Translate repeated common roots and left-pair codegrees into identities among
the \(\mathrm{PGL}_2\) continuant transfer matrices.  Use the recurrence and
reflection involution to seek a uniform bound on separated collision
multiplicity or codegree, matching the exhaustive census through \(p=20000\).

- Success: prove the Ap\'ery-specific weighted split affine gcd-tail bound
  (and the complementary-annulus estimate) strong enough to give
  \(E_p(H)\ll H^{3/2+o(1)}\).
- Proof-of-failure: construct an Ap\'ery realization of the abstract
  affine-plane obstruction, or derive a quantifier-preserving reduction to
  the original zero-set problem.

### (e0) Distinguished-orbit coupling after phantom removal

**Current active avenue.**  Fully deflating all three adjacent gaps gives a
rigorous alternative hypothesis \(\mathrm{FDAC}^{\mathrm{tr}}\) and loses at
most three centered sliding quadruples per prime.  The height-36 endpoint-gcd
classification finds a primitive reflected phantom pair
\[
 (a,b,c,p,r)=(5,20,10,1297,360),\quad(10,20,5,1297,901).
\]
Here the span is \(35\) and \(1297>35^2\), but the actual zero set is
\(\{459,530,766,837\}\).  Thus even the primitive formal-projective analogue
of the short-window statement is false; the distinguished-orbit condition is
essential.  Among all nine height-36 endpoint candidates in the short range,
the remaining cases are three endpoint false positives and four nonprimitive
roots, and none is an actual Ap\'ery-zero start.

The sharper active subtarget is the **short-window reflection principle**:
there is a fixed \(\eta>0\) such that four consecutive actual Ap\'ery zeros
of span \(s\le\eta\sqrt p\) contain a centered adjacent pair.  The extended
\(350{,}104\)-pair census verifies the stronger choice \(\eta=1\) through
\(p\le5\cdot10^6\); among its \(1{,}418\) off-center windows, the minimum
observed \(s^2/p\) is still \(1428025/3727>383\).  For any such fixed
\(\eta\), partition into blocks of
cardinality \(\lfloor\eta\sqrt p\rfloor+1\) to obtain
\[
 Z(p)\le
 3\left\lceil
   \frac{p}{\lfloor\eta\sqrt p\rfloor+1}
  \right\rceil+3
 <3\eta^{-1}\sqrt p+6.
\]
Indeed, at most three selected sliding quadruples can contain the unique
centered adjacent pair.  Without assuming the principle, the same span count
gives the exact reduction
\(Z(p)<3\eta^{-1}\sqrt p+6+E_{\mathrm{sw}}(p;\eta)\), where
\(E_{\mathrm{sw}}(p;\eta)\) counts off-center windows of span at most
\(\eta\sqrt p\).  Thus the
orbit-coupling problem can be
replaced by this precise actual-zero statement; it cannot be replaced by the
unit-scale formal-projective statement even after retaining primitivity.

There is also an alternative purely algebraic target, **quadratic primitive-return
support (QPRS)**.  It asks only for an absolute \(C\) such that every
primitive off-center return chain with a non-all-equal gap vector outside
\(\mathfrak U_s\) satisfies \(p\le Cs^2\).  The constant need not be~1:
on zero-count levels
\(T\ge K\sqrt X\), take \(K\gg\sqrt C\) and
\(H=\lceil16X/T\rceil\).  Then \(CH^2<X<p\), so QPRS rules out every
generic selected quadruple; the remaining \(\mathfrak U_H\)-supported primes
have level mass \(O(X^2/(T\log X))\), whose dyadic sum is at the target
scale.  Hence QPRS alone implies the square-root prime average.

The formal evidence now makes this a secondary conditional target rather
than the main avenue.  An exact prime-first census of all projective fibers
for every \(p\le500000\) finds twenty primitive off-center chains with
\(p>s^2\), four with \(p>2s^2\), and maximum
\[
 \frac{p}{s^2}=\frac{128047}{164^2}=4.760819\ldots
\]
at gaps \((41,86,37)\); all twenty are phantoms.  Thus \(C=2\) is false.
Moreover, the random-fiber occupancy heuristic gives an expected
\(s^3/p^2\) short four-collisions per prime.  At
\(s=\sqrt{p/C}\), its sum through \(P\) has order
\(\sqrt P/(C^{3/2}\log P)\), suggesting that every fixed \(C\) may
eventually fail.  This is not a proof against QPRS, but it removes the
formal-projective condition from the favored route.

- Success: prove the short-window reflection principle, or construct a
  bounded-height orbit-coupled module whose local corank dominates actual
  primitive starts and whose aggregate height is \(H^{3+o(1)}\).  A proof of
  QPRS would still suffice, but the formal census makes it less credible.
- Proof-of-failure: either construct primitive formal chains with unbounded
  \(p/s^2\), thereby retiring QPRS, or give a quantifier-preserving reduction
  showing that every distinguished-orbit coupling is equivalent to the
  parity-barrier atom problem.

If all live avenues terminate negatively, retain the fully verified
conditional implications under the four stated hypotheses, record the
aggregate content and distinguished-orbit inputs that remain open, and
continue with genuinely new invariants;
finite-depth Cartier lifting and the two withdrawn scalar gcd formulations
remain closed and are not to be revived.

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

Keep only conditional statements whose implications are proved completely,
correct any unconditional wording, and record the exact aggregate-content or
orbit-coupling hypothesis that remains.  Do not replace the full conjecture by
the weaker ($p^{1/2}$) average-zero milestone.

## Existing terminal facts that constrain every avenue

- Fixed Cartier/cofactor depth is already computed through precision (p^8)
  and remains CRT reconstruction.
- The fourth- and fifth-digit companion laws are rank one and CRT-compatible.
- Pair certificates, reflection, Mellin twists, fixed-depth height arguments,
  and unsaturated resultant mass do not break certificate parity.
- The current paper's exceptional-class lemma counts primes but does not by
  itself bound their (Z(p))-weighted mass.
