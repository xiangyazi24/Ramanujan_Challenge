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

The companion-label refinement is now terminal in its basis-covariant form.
Let \(B=(1,5)^t\), \(C=(0,1)^t\), let \(M_r\) transport the initial state to
index \(r\), and put \(R_r=e_1^tM_r\) and
\(L(x,y)=y-5x\).  At a distinguished zero, \(R_rB=b_r=0\), while the
Wronskian gives \(R_rC=c_r\ne0\).  Since \((B,C)\) is a basis and
\(L(B)=0,L(C)=1\), one has the exact marked row identity
\[
 R_r=c_rL.
\]
For two returns \(r,r+h\), composing the intervening transfer block recovers
only
\[
 \frac{c_{r+h}}{c_r}
 =-\frac{\prod_{j=1}^{h-1}(r+j)^3}
 {(r+h)^3N_{h-1}(r)}.
\]
For four returns, every direct-versus-segmented comparison is the cocycle
law for these ratios.  It is independent of the finite projective fiber
\(b-qc\).  Differentiation in the marked coordinate adds only
\(\partial_q(b_r-qc_r)|_{q=0}=-c_r\ne0\): the distinguished fiber is a
simple root for each fixed return equation, but no further relation among
the return indices follows.  More precisely, for four evaluation pairs
\((B_i,C_i)\), localization at \(C_0\) gives
\[
 (B_i-qC_i:0\le i\le3)
 =
 (C_0q-B_0,\ B_iC_0-B_0C_i:1\le i\le3).
\]
Thus marking contributes exactly the pin \(B_0=0\), followed by the three
ordinary projective-return minors; the nonzero \(q\)-derivative supplies no
fifth equation.  The carrier-clean certificate
\(p=709,q=67,r=282,300\) is especially sharp: the complete normalized pair
\((b_r,c_r)\) repeats across a nonreflection gap \(18<\sqrt{709}\), with the
reflected repeat at \(408,426\).  Thus neither companion ratios nor equality
of normalized labels can distinguish the Ap\'ery fiber; any successful
marked argument must import arithmetic information beyond fixed-return
transfer algebra.

An independent index-first computation gives a finite, prime-unbounded
check of the actual rather than projective statement.  Computing the exact
integers \(b_0,\ldots,b_{5000}\), factoring every
\(\gcd(b_i,b_j)\) with \(2\le j-i\le500\), and reconstructing every relevant
zero prefix produces \(197\) candidate primes and \(22\) consecutive
four-zero windows.  Every window contains a centered adjacent pair; there
is no off-center window in this box.  The smaller \((3000,300)\) run also
records \(46\) four-subsets inside short blocks but only \(16\) consecutive
windows, quantifying the overcount in an all-four-subsets statistic.  This is
a useful complement to the prime-first census, but it is a finite certificate
rather than an asymptotic input.

The deterministic pair caps and reflection symmetry cannot by themselves
prove any subquadratic bound for the short off-center windows.  The explicit
construction in `pair_cap_extremal.tex` takes \(H=12m\), repeats the gap word
\(3m,3m+1,\ldots,4m-1\) for \(m/100\) periods, places its reflection beyond
a central gap, and chooses a prime \(p>4(L+1)+1>H^2\).  It omits both
endpoints, has no adjacent points,
obeys \(A(h)\le3(h-1)\) for every gap, but has exactly
\[
 2\left(\frac{m^2}{100}-2\right)
 =\frac{H^2}{7200}-4
\]
consecutive short off-center windows.  For \(h<m^2\), every admissible block
has length below one period and each length represents \(h\) at most \(2R\)
times, giving \(A(h)\le h/75\).  For \(h\ge m^2\), the trivial one-endpoint
per start bound is already below \(3(h-1)\); the cross-half differences have
the same property.  Therefore a proof of the short-window principle must use
an Ap\'ery-specific four-point constraint, not just reflection, endpoint
exclusion, no adjacency, and the individual continuant degree bounds.

The same construction gives a sharper negative result for any argument that
uses only these set-theoretic constraints and treats the state labels as
otherwise unconstrained.
For every sufficiently large prescribed prime \(p\), one may take a multiple
\(m\asymp p^{1/3}\) of \(100\) small enough that the two reflected packets
fit inside \([0,p-1]\).  Declaring the entire set to be one abstract state,
the \(2(N-2)\) short windows form \(k=N-2\asymp p^{2/3}\) reflection
orbits.  Two orbit supports overlap exactly when their left-window indices
differ by at most three, and hence
\[
 E_p^{\mathrm{sep}}=k(k-1)-(6k-12)=(k-3)(k-4)\asymp p^{4/3}.
\]
Doing this independently at every \(p\in(X,2X]\) permits abstract dyadic
separated energy \(\gg X^{7/3}/\log X\).  This is not an actual Ap\'ery
orbit and therefore gives no lower bound for the genuine weak-DPLS variance:
the common label was assigned artificially, rather than produced by the
Ap\'ery cocycle.  Its precise conclusion is that pair caps, reflection, no
adjacency, off-centeredness, and the non-AP filter cannot suffice without a
state-realizability axiom.  The exact construction and proof are in
`pair_cap_extremal.tex`; an Ap\'ery-specific carrier or transfer constraint
is indispensable.  ChatGPT audit Q7195 identified this distinction and
confirmed the packet geometry and pair-cap calculation.

The leading prime-averaged formulation is now a projective variance problem
that retains exact primitivity.  For each state
\(q\in\mathbf P^1(\mathbf F_p)\), let \(C_p(q)\) count the exact
carrier-clean, non-all-equal, off-center primitive chains in that state with
span at most \(\lfloor\sqrt p\rfloor\), and put
\[
 M_p=\sum_q C_p(q),\qquad
 V_p=\sum_q\left(C_p(q)-\frac{M_p}{p+1}\right)^2.
\]
Each start supports at most one such chain, so \(M_p\le p\).  The actual
Ap\'ery contribution is \(A_p=C_p(0)\), and hence
\[
 \sum_{X<p\le2X}A_p
 \le \sum_{X<p\le2X}\frac{M_p}{p+1}
 +\left(\frac{X}{\log X}\sum_{X<p\le2X}V_p\right)^{1/2}.
\]
Consequently the absolute weak-DPLS estimate
\[
 \sum_{X<p\le2X}V_p\ll\frac{X^{2+o(1)}}{\log X}
\]
already suffices; the relative diagonal estimate
\(\sum V_p\ll X^{o(1)}\sum M_p\) is stronger than necessary.

It is cleaner to remove the automatic reflection diagonal before attacking
this estimate.  The involution
\[
 (x;a,b,c)\longmapsto(p-1-x-a-b-c;c,b,a)
\]
preserves the projective state by the exact orbit reflection
\(\pi(n)=\pi(p-1-n)\).  It has no fixed point on the off-center chains:
a fixed chain would have \(a=c\) and \(2x+a+b+c=p-1\), making its middle
adjacent pair centered.  If \(\bar C_p(q)\) counts reflection orbits and
\(\bar M_p,\bar V_p,\bar A_p\) denote the corresponding statistics, then
\[
 C_p(q)=2\bar C_p(q),\qquad
 M_p=2\bar M_p,\qquad V_p=4\bar V_p,\qquad A_p=2\bar A_p.
\]
Thus weak-DPLS is equivalent, up to an absolute factor, to the same estimate
for \(\bar V_p\).

The exact overlap contribution is already within budget.  If
\(r_1<\cdots<r_t\) are the occurrences of a fixed state, every primitive
chain is a filtered sliding window
\(\gamma_i=(r_i,r_{i+1},r_{i+2},r_{i+3})\).  A fixed \(\gamma_i\) can have
overlapping closed return interval only with
\(\gamma_{i+j}\), \(0<|j|\le3\).  Hence the ordered distinct overlap energy
is at most \(6M_p\), and the diagonal plus overlap contribution over
\(X<p\le2X\) is
\[
 \le7\sum_{X<p\le2X}M_p\ll\frac{X^2}{\log X}.
\]
The same constant survives the reflection quotient: after fixing one
oriented representative of an orbit, any overlapping second orbit has a
representative among those same six neighboring windows.  Thus the quotient
square energy is at most
\(7\bar M_p+E_p^{\mathrm{sep}}\), where the last term counts ordered pairs
of distinct reflection orbits for which no two representatives have
intersecting closed intervals.  The complete statement, including the
Cauchy--Schwarz reduction of the distinguished fiber, is proved in
`projective_variance_reduction.tex`.
Unlike the distinguished value, the variance is invariant under a change
of solution basis, so formal projective tools are legitimate here.  After
the reflection quotient and overlap estimate, the exact obstruction is only
the energy of separated distinct chains.  Such a pair adds the long return
\(N_L(x)=0\), with \(L\) as large as \(p\); this defeats all currently
bounded-gap continuant certificates.

The near part of this obstruction is nevertheless unconditional.  For a
separated ordered orbit pair, choose canonically a representative of the
first orbit lying before a representative of the second and write \(s\) for
the first span and \(G\) for the intervening bridge.  Consecutiveness makes
\((p,x,s,G)\) injective: the second chain, including its span, is already
determined by its start \(x+s+G\).  Hence it is enough to count common roots
of
\[
 N_s(x),\qquad N_G(x+s).
\]
Their integer resultant is nonzero by the disjoint complex root strips, and
its logarithmic height is \(O(sG\log H)\).  Summing Smith nullities proves
\[
 \sum_{X<p\le2X}E_p^{\mathrm{near}}(K)\ll H^2K^2,
 \qquad H=\lfloor\sqrt{2X}\rfloor.
\]
Thus the full weak-DPLS budget is already met for
\(G\le H/\sqrt{\log X}\).  The proof is in `near_bridge_energy.tex`; the
remaining theorem concerns only genuinely long bridges.  An independent
tmux-11 audit confirmed the canonical orientation, the Smith-nullity
inequality under degree drop, and the resultant-height calculation.

The long part admits a sharper exact reduction before any resultant is
introduced.  Let \(B_p(K)\) count a selected oriented first chain together
with an arbitrary later occurrence of its endpoint state at bridge greater
than~\(K\), and let \(A_p^{\rm nw}(K)\) count nonwrapping triples
\[
 x<y<z,\qquad 2\le y-x\le\lfloor\sqrt p\rfloor,\qquad z-y>K,
 \qquad \pi_p(x)=\pi_p(y)=\pi_p(z).
\]
Canonical orientation and consecutiveness give injective maps, with no
reflection factor,
\[
 E_p^{\rm far}(K)\le B_p^{[4]}(K)\le B_p^\circ(K)
 \le B_p(K)\le A_p^{\rm nw}(K).
\]
Here \(B_p^\circ\) removes every later point in the first chain's own
reflection-invariant support, and \(B_p^{[4]}\) further requires the external
point to begin a four-consecutive-occurrence window of span at most
\(\lfloor\sqrt p\rfloor\).  Both restrictions are necessary for the second
member of a separated quotient pair.  Moreover,
\[
 0\le B_p(K)-B_p^\circ(K)\le4M_p,
\]
so the removed automatic reflection baseline is already within the diagonal
mass budget.
The final majorant is exactly
\[
 \sum_{s=2}^{\lfloor\sqrt p\rfloor}
 \sum_{G=K+1}^{p-1-s}
 \#\{0\le x\le p-1-s-G:
       N_s(x)=N_G(x+s)=0\pmod p\}.
\]
Only after dropping the displayed interval restriction may it be bounded by
the corresponding full-\(\mathbf F_p\) common-root mass.  In particular,
\(\deg\gcd_{\mathbf F_p[x]}(N_s(x),N_G(x+s))\) is not the exact observable:
it can also count nonsplit irreducible factors.  The complete proof is in
`far_bridge_incidence.tex`.

There is a stronger unconditional bound than the termwise degree estimate.
For one projective fiber, let $e_h$ count consecutive occurrence gaps of
length $h$.  Then $e_h\le3(h-1)$ and
$\sum_h h e_h\le p-1$.  Splitting at $h=p^{1/3}$ gives
$\max_q|\pi_p^{-1}(q)|\le4p^{2/3}$.  Since the total number of
same-fiber pairs of gap at most $\sqrt p$ is at most
$3\sum_{h\le\sqrt p}(h-1)\le3p/2$, this proves
\[
 A_p^{\rm nw}(K)\le6p^{5/3},
 \qquad
 \sum_{X<p\le2X}A_p^{\rm nw}(K)
 \ll \frac{X^{8/3}}{\log X}.
\]
This improves the raw $X^3/\log X$ degree bound, but remains a factor
$X^{2/3}$ above the separated-energy budget.

The full orbit-energy theorem and a reciprocal-gap refinement improve this
three-point majorant further.  If the consecutive gaps in the $q$-fiber are
$g_{q,i}$, put
\[
 T_q=\sum_i g_{q,i}^{-2}.
\]
The global fixed-gap caps give
$\sum_qT_q\le3\log p$, while H\"older gives
$m_q-1\le p^{2/3}T_q^{1/3}$.  Partitioning the index interval into blocks
of length $H=\lfloor\sqrt p\rfloor$ and applying the same inequality inside
each block controls the number $P_q(H)$ of short pairs by
\[
 P_q(H)\le3m_q
   +3H^{4/3}\lceil p/H\rceil^{1/3}T_q^{2/3}.
\]
Together with the full fiber-square estimate
$\sum_qm_q^2\ll p^{3/2}$, this proves, uniformly in~$K$,
\[
 A_p^{\rm nw}(K)\ll p^{3/2}\log(2p),
 \qquad
 \sum_{X<p\le2X}A_p^{\rm nw}(K)\ll X^{5/2}.
\]
The proof is in `far_bridge_incidence.tex`.  This is still a factor
$X^{1/2}\log X$ above the target and is one logarithm weaker than the
direct short-chain energy estimate below.

Keeping the two short four-windows gives a better unconditional energy
bound.  If \(e_q(h)\) counts consecutive fiber edges of length \(h\le H\)
and \(E_q=\sum_{h\le H}e_q(h)\), then every selected sliding four-window
charges three edges and every edge is charged at most three times.  Hence
\(C_p(q)\le E_q\).  The global gap caps and the total length of one fiber
give
\[
 \sum_qe_q(h)\le3(h-1),\qquad
 \sum_{h\le H}he_q(h)\le p-1.
\]
For \(t\ge4p/H\), splitting at \(L=\lfloor2p/t\rfloor\) shows that
\[
 \#\{q:E_q\ge t\}\le\frac{12p^2}{t^3}.
\]
A dyadic level decomposition then proves
\[
 \sum_qC_p(q)^2\le30pH,\qquad
 E_p^{\rm sep}\le\frac{15}{2}pH.
\]
At \(H=\lfloor\sqrt p\rfloor\), this gives the dyadic bound
\[
 \sum_{X<p\le2X}E_p^{\rm sep}
 \ll\frac{X^{5/2}}{\log X}.
\]
This is stronger than the three-point \(X^{8/3}/\log X\) majorant, but
still exceeds the required scale by \(X^{1/2}\).  The proof is in
projective_variance_reduction.tex.  A direct regression through all
92 primes \(7\le p\le499\) has fail-closed digest
0796e0f34777aba4a631767076dbedb71503eef4422e19553ad09e702be6566c.

There is a mass-sensitive interpolation, but it identifies a substantially
stronger missing estimate rather than closing the gap.  If
\(M_p=\sum_qC_p(q)\), the same cubic level-set bound gives
\[
 \sum_qC_p(q)^2\le20p\sqrt{M_p},\qquad
 E_p^{\rm sep}\le8p\sqrt{\bar M_p}.
\]
Consequently, with
\(S_X=\sum_{X<p\le2X}\bar M_p\),
\[
 \sum_{X<p\le2X}E_p^{\rm sep}
 \le16X\sqrt{\#\mathcal P_X\,S_X}.
\]
This reaches the target if
\(S_X\ll X^{1+o(1)}/\log X\), namely if the average number of selected
short chains is \(X^{o(1)}\).  By contrast, even a natural
\(H^{3+o(1)}\)-scale prime-averaged certificate bound at
\(H\asymp\sqrt X\) would give only
\(S_X\ll X^{3/2+o(1)}/\log X\) and hence separated energy
\(X^{9/4+o(1)}/\log X\).  Thus ordinary three-gap content control is still
a factor \(X^{1/2}\) too weak at the mass level; the required statement is
genuine short-chain sparsity.

This exponent cannot be improved from the preceding set-theoretic inputs.
For infinitely many primes, `pair_cap_extremal.tex` constructs a complete
reflection-invariant coloring satisfying every global fixed-gap cap
\[
 \sum_q\#\{x:\kappa_p(x)=\kappa_p(x+h)=q\}\le3(h-1)
\]
and selected four-consecutive-occurrence windows with gaps at least~\(3\),
nonconstant gap triples, and no centered adjacent pair, but with
\[
 E_p^{\rm sep}\ge \frac{pH}{1920}.
\]
Thus any saving over \(pH\) must use Ap\'ery-specific arithmetic, not only
the global pair caps and local window axioms.  The row-shift construction
and all-pair verification have fail-closed digest
707f006bbafe1da79d69fb859f76b78b885d238c9bacda97e47043981600e39e;
an independent tmux-11 audit checked the partition, every gap range,
reflection quotient, and constant.

The same complete coloring also closes the purely combinatorial route
through the strongest conditioned far majorant.  Its full fiber-square
energy is $O(p^{3/2})$.  With
$H_p=\lfloor\sqrt p\rfloor$ and
$K=\lfloor H_p/\sqrt{\log p}\rfloor$, take a left-packet selected window
and a later left-packet occurrence, at least five start rows away, which
begins another selected window.  These raw incidences are external to the
first reflection support and give
\[
 B_{\kappa_p}^{[4]}(K)
 \ge \frac{pH}{3840}\gg p^{3/2}.
\]
Thus retaining a short four-window at both ends, the global pair caps, and
the full $p^{3/2}$ fiber energy still cannot yield a sub-$p^{3/2}$ bound
without Ap\'ery-specific arithmetic.

Direct full-cycle multiplication supplies no further bridge equation.
For the factorially renormalized solution row
\(\mathbf v_n=(n!^3u_n,n!^3w_n)\), a return \(N_h(x)=0\) has the exact
affine multiplier
\[
 \mathbf v_{x+h}=\mu_h(x)\mathbf v_x,\qquad
 \mu_h(x)=-(x+1)^6N_{h-1}(x+1).
\]
Writing \(x^\vee=p-1-x-h\), reflection and Wilson's theorem give
\[
 \mu_h(x)\mu_h(x^\vee)
 =(-1)^h\left(\prod_{j=1}^h(x+j)\right)^6.
\]
Thus the bridge and its reflected bridge multiply to the determinant unit
already prescribed by the factorial gauge, rather than to a second
vanishing condition.  The transfer step across \(n=0\) is rank one, so it
cannot be inverted to turn this identity into a complementary-gap return.
This rules out the direct reflection/full-cycle rank-amplification attempt;
it does not rule out an Apéry-specific relation carrying genuinely new
higher-order data.  The proof is in full_cycle_bridge_reciprocity.tex.
The independent verifier checks every return for the 24 primes
\(5\le p\le101\) and all four known return offsets at
\((p,x)=(1297,360)\); its fail-closed digest is
4faf19b270b89c14f9e6e584fd2631bd79a242fc65704f902e9dfbcd78f98b07.

The occurrence-list identity makes these majorants directly computable.
For a fiber \(r_0<\cdots<r_t\), put
\[
 L_i=\#\{j<i:2\le r_i-r_j\le\lfloor\sqrt p\rfloor\},
 \qquad
 R_i(K)=\#\{k>i:r_k-r_i>K\}.
\]
Then its contribution to \(A_p^{\rm nw}(K)\) is
\(\sum_iL_iR_i(K)\), while its contribution to \(B_p(K)\) is
\(\sum_i\mathbf1_{\{(r_{i-3},\ldots,r_i)\ \mathrm{selected}\}}R_i(K)\).
Deleting the at most four positions in the same reflection orbit gives
\(B_p^\circ(K)\); retaining only external positions \(r_j\) with
\(j+3\le t\) and \(r_{j+3}-r_j\le\lfloor\sqrt p\rfloor\) gives
\(B_p^{[4]}(K)\).
The exact C++ implementation `long_bridge_incidence_scan.cpp` gives the
following complete dyadic data, with
\(H=\lfloor\sqrt{2X}\rfloor\) and
\(K=\lfloor H/\sqrt{\log X}\rfloor\):

| \(X\) | primes | \(K\) | \(A^{\rm nw}(K)\) | \(B(K)\) | \(B^\circ(K)\) | \(B^{[4]}(K)\) | raw chains | \(A^{\rm nw}/(\#p\sqrt X)\) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1,000 | 135 | 16 | 8,356 | 4 | 0 | 0 | 2 | 1.9573 |
| 5,000 | 560 | 34 | 81,018 | 6 | 2 | 0 | 2 | 2.0460 |
| 10,000 | 1,033 | 46 | 214,495 | 12 | 4 | 0 | 4 | 2.0764 |
| 20,000 | 1,941 | 63 | 575,531 | 0 | 0 | 0 | 0 | 2.0967 |
| 50,000 | 4,459 | 96 | 2,099,977 | 8 | 4 | 0 | 2 | 2.1062 |
| 100,000 | 8,392 | 131 | 5,614,535 | 16 | 8 | 0 | 4 | 2.1157 |

The source digest is
`8a9b8614e27d0a841a70a692daa77bd25b331c930666d646cc9b2bae493191d6`.
It passes strict `clang++` compilation and ASan/UBSan.  The independent
quadratic verifier `long_bridge_incidence_verify.py` does not use the
occurrence-list counting algorithm: it enumerates every pair and triple and
also checks pointwise, for every nonwrapping \((x,h)\) in its default dyad,
that \(N_h(x)=0\) if and only if \(\pi_p(x)=\pi_p(x+h)\).  Its \(X=70\)
regression passes with output digest
`9125a62abc9121a31e4190f0a10c1650fd7742b0360b79036345ace34a410fe4`.
Finally, `long_bridge_selected_extension_verify.py` recomputes the fibers of
all 20 selected raw records in the complete \(p\le500000\) census.  It finds
22 external later occurrences after reflection-support deletion, but none
begins another short four-occurrence window; its fail-closed digest is
`13f2f869de772d2f30de9d78a1818ad149e9cef4f92adec9ec2eef7c4e25ea1b`.
The stable normalization in the table is consistent with
\(A^{\rm nw}(K)\asymp(\#p)\sqrt X\), already far below the required budget,
and the conditioned count is much smaller still.  This remains finite
evidence: neither scanner supplies the missing uniform arithmetic estimate.

Two possible height/compression shortcuts have now been audited exactly.
First, `separated_resultant_deflation_probe.py` distinguishes the full
resultant, the resultant after deleting every forced even center factor, and
the remaining integer after deleting the maximal part supported on the
structural carrier.  At \(H=12\) and \(2\le G\le36\), the aggregate bit
lengths are respectively

\[
 2699303,\qquad 2604380,\qquad 2113651.
\]

Thus the tested residual retains \(0.7830\) of the full bit length, uniformly
across the three bridge buckets, and no tested resultant is removed
completely.  On the diagonal \(2\le h\le12\), the formal root-product term
\(\operatorname{lc}(N_h)^{2\deg N_h}\) never divides the integer resultant;
it is an archimedean decomposition, not a removable carrier factor.  The
center-and-carrier residual has
\(\log|R^*_{h,h}|/(h^2\log h)=10.26,10.47,10.69\) at \(h=8,10,12\).
These are finite exact data, not an asymptotic lower bound, and bit length
alone gives no information about support on primes \(p\asymp X\).

Second, a hidden low-order polynomial bispectral operator is absent in a
large explicit ansatz.  If

\[
 L=\sum_{j=-r}^{r}A_j(x)T^j,\qquad \deg A_j\le d,\qquad
 LN_h=\lambda_hN_h,
\]

then the exact coefficient system for \(1\le h\le20\), reduced modulo
\(1000003\), has rank \(546\) in \(547\) columns for \(r=8,d=30\).
Its kernel is therefore exactly the scalar identity over \(\mathbb Q\), which
also excludes every sub-ansatz \(r\le8,d\le30\).  The reproducible certificate
is `bispectral_operator_scan.py`.  The same script also handles a common
rational denominator without solving a bilinear system.  After clearing the
denominator, any rational eigenoperator must satisfy
\[
 N_h(x)\mid\sum_{j=-r}^{r}A_j(x)N_h(x+j).
\]
The corresponding modular remainder matrix has shape \(570\)-by-\(527\),
rank \(496\), and kernel dimension \(31=d+1\), exactly the subspace of
multiplication numerators \(A_0(x)\).  Hence every common-denominator rational
operator whose cleared numerators have \(r\le8,d\le30\) is scalar.  This does
not exclude a gauge whose conjugated cleared numerators exceed the degree
bound, or operators of larger order.

There is nevertheless an exact fixed-order differential compression of one
moving gap, but it is an inhomogeneous Green kernel rather than a bispectral
eigenfunction.  Put
\[
 V_m(a)=\frac{U_m(a)}{(a+1)_m^3},\qquad
 F_a(z)=\sum_{m\ge0}V_m(a)z^m,
 \qquad \theta=z\frac{d}{dz}.
\]
Then
\[
 (a+m+1)^3V_{m+1}=P(a+m)V_m-(a+m)^3V_{m-1}
\]
and exact summation, including the initial term, gives
\[
 \bigl((\theta+a)^3-zP(\theta+a)
       +z^2(\theta+a+1)^3\bigr)F_a=a^3.
\]
Thus, for the usual Ap\'ery operator
\(\mathscr L_0=\theta^3-zP(\theta)+z^2(\theta+1)^3\),
\[
 \mathscr L_0(z^aF_a)=a^3z^a.
\]
The homogeneous eigenfunction equation proposed in ChatGPT Q7231 is false:
it omitted the boundary term at \(m=0\).  After adjoining the constant
solution, the corrected equation is a rank-four system with fixed singular
support
\(z(1-34z+z^2)=0\), independent of \(a\).  Moreover
\[
 V_m(a)=a^3(b_{a-1}c_{a+m}-c_{a-1}b_{a+m}),
\]
so it is exactly the two-point Green--Casoratian kernel.  For a long bridge,
\[
 N_G(u)=((u+2)_{G-1})^3V_{G-1}(u+1).
\]
This identifies a fixed differential object, but does not yet control the
moving coefficient \([z^{G-1}]F_{u+1}\): a naive coefficient-jet
algebraization still has order growing with \(G\).  The corrected theorem and
proof are in `gap_kernel_green.tex`; `gap_kernel_green_verify.py` checks the
symbolic recurrence, every differential coefficient through degree seven,
the rejected eigenfunction equation, and 110 exact Casoratian values.

The starting index and gap index can be compressed simultaneously at the
formal level.  If
\[
 \mathcal G(s,z)=\sum_{a\ge1,m\ge0}V_m(a)s^az^m
\]
and
\[
 \mathcal K(X,Z)=X\bigl((\theta_X+1)^3B(X)C(Z)
                  -(\theta_X+1)^3C(X)B(Z)\bigr),
\]
then coefficientwise
\[
 \mathcal G(s,z)=\mathcal P_z^{\ge0}\mathcal K(s/z,z)
 =\operatorname{CT}_u\frac{\mathcal K(su/z,z/u)}{1-u}.
\]
Writing \(\Theta=\theta_s+\theta_z\), it also satisfies the exact PDE
\[
 \bigl(\Theta^3-zP(\Theta)+z^2(\Theta+1)^3\bigr)\mathcal G
 =\frac{s(1+4s+s^2)}{(1-s)^4}.
\]
This is a genuine fixed formal two-index encoding, but not a rational-period
shortcut.  For every prime \(p\ge5\), reflection and the recurrence at
\(p-1\) give
\[
 v_p(c_p)=-3,\qquad p^3c_p\equiv1\pmod p.
\]
Since \([s z^{p-1}]\mathcal G=c_p\), infinitely many new denominator primes
occur.  Hence neither \(C\) nor \(\mathcal G\) is globally bounded, and
\(\mathcal G\) cannot be the ordinary diagonal or constant term of a fixed
rational function regular at the expansion point.  This does not exclude a
mixed/regulator realization or a characteristic-dependent finite-field
model.  The proof is in `gap_kernel_bivariate.tex`; the exact verifier and an
independent tmux-11 audit both pass.

The finite nonwrapping truncation has an exact Frobenius-boundary law.  For
\[
 F_{p,a}^{\triangle}(z)=
 \sum_{m=0}^{p-a-1}V_m(a)z^m,
\]
reduced on the strict triangle, one has
\[
 \mathscr D_aF_{p,a}^{\triangle}
 =a^3-a^3b_{a-1}z^{p-a},\qquad
 \mathscr L_0(z^aF_{p,a}^{\triangle})
 =a^3z^a-a^3b_{a-1}z^p.
\]
Thus one Cartier extraction reads back \(-a^3b_{a-1}z\).  Summing over
\(a\) places the entire differentiated Hasse polynomial on total degree
\(p\).  On a marked fiber \(d_{a-1}=0\), the same Green column is exactly
\(d_{a+m}/d_a\), so disappearance of the boundary is the original future
zero problem, not a new condition.  Complete Mellin and Cartier extraction
are exact but use a moving character/frame, while
\(\deg F_{p,a}^{\triangle}\ge p-a-2\).  This proves that the natural
rank-four Green/Dwork compression is circular at its first Frobenius
boundary; it does not exclude a genuinely new characteristic-dependent
mixed crystal.  The theorem and its independent verifier are
`gap_kernel_cartier_defect.tex` and
`gap_kernel_cartier_defect_verify.py`.

Distributing the necessary three factorial Borel orders between the gap and
total-degree directions does not remove this denominator.  The four natural
gauges are
\[
 W_m^{(r)}(a)=
 \frac{U_m(a)}{(m!)^{3-r}(a+1)_m^r}
 =\binom{a+m}{m}^{3-r}V_m(a),\qquad0\le r\le3.
\]
They are mutually unit-equivalent on every strict triangle \(a+m<p\), but
\[
 W_p^{(r)}(1)=(p+1)^{3-r}c_{p+1},
 \qquad v_p(c_{p+1})=-3
\]
for every prime \(p\ge7\).  Hence none is globally bounded, algebraic, or an
ordinary rational diagonal/constant term regular at the origin.  This no-go
is not universal.  With \(\lambda_n=\operatorname{lcm}(1,\ldots,n)\), the
arithmetic gauge
\[
 \Xi_m(a)=6\lambda_{a+m}^3V_m(a)
\]
is integral, has only exponential growth, and is still unit-equivalent to
\(U_m(a)\) whenever \(a+m<p\).  Its lcm multiplier is not a finite
Gamma/Borel gauge, since \(\lambda_{n+1}/\lambda_n\) equals one off prime
powers and has a prime spike at every prime power.  These spikes give a
strict obstruction: a CRT--Dirichlet isolation argument proves that
\((6\lambda_n^3)\) is not P-recursive, so the bivariate lcm-gauged series is
not D-finite.  Moreover, any polynomial over \(\mathbb F_p\) representing
its zeroth gap slice has degree at least
\[
 p-\#\{q^e:2\le q^e<p\}=p-O(p/\log p).
\]
Thus an integral positive-radius encoding of the exact nonwrapping zero
relation exists, but it has asymptotically full direct finite-field
complexity.  Since the gauge is a unit on the strict triangle, it leaves
every far-bridge incidence and selected short-chain window exactly
unchanged.  The theorem and verifier are
`gap_kernel_borel_gauges.tex` and `gap_kernel_borel_gauges_verify.py`; an
independent audit passes.

The toric fiber admits a second exact collapse, now in the coefficient
variable.  For
\[
 A_j(t)=\sum_{k=0}^j\binom jk^2t^k
\]
and \(1\le j\le p-2\), direct elimination of two torus coordinates gives
\[
 b_j=-\sum_{t\in\mathbb F_p^\times}
       t^jA_j(t)A_{p-1-j}(t).
\]
The one-dimensional character sum behind this formula is exactly
\[
 \sum_{r\ne0}r^j
 \chi\bigl((t-1)^2r^2-2(t+1)r+1\bigr)
 =-A_{p-1-j}(t).
\]
For \(1\le j\le(p-1)/2\), Legendre duality has no extra sign and yields
\[
 A_{p-1-j}(t)=(1-t)^{p-1-2j}A_j(t),
 \qquad
 b_j=-\sum_{t\ne0}t^j(1-t)^{p-1-2j}A_j(t)^2.
\]
This is an indefinite finite-field square pairing, not a positive norm.
Expanding the \(t\)-sum recovers the original Ap\'ery binomial formula
coefficient by coefficient, so the collapse alone supplies no new zero-set
bound.  The exact derivation is in `toric_mellin_square.tex`;
`toric_mellin_square_verify.py` checks 66 fibers, 7440 inner character sums,
252 pairings, 132 dualities, and 252 tautology expansions.  An independent
tmux-11 audit checked the coordinate exclusions, both minus signs, the
middle endpoint, and the absence of a spurious \((-1)^j\) factor.

The same toric coordinates also improve the pointwise fiber estimate.  The
closure of a nonzero fiber in \((\mathbb P^1)^3\) is the explicit
anticanonical \((2,2,2)\) surface
\[
 U_1V_1(W_1-W_0)(U_0V_0W_0-U_1V_1W_1)
 -aU_0V_0W_0W_1(U_1-U_0)(V_1-V_0)=0.
\]
Its boundary is a union of sixteen rational curves with exactly \(16p-11\)
rational points.  It has four fixed \(A_1\)'s, two fixed \(A_2\)'s, a central
\(A_1\) which becomes \(A_2\) at \(a=1\), and one additional \(A_1\) exactly
when \(a^2-34a+1=0\).  Resolving these split ADE points gives a K3 surface
\(\widetilde X_a\).  If \(\tau_p(a)\) is its \(H^2\) Frobenius trace and
\(\varepsilon_a,\delta_a\) are the two displayed exceptional indicators,
then the character sum in the quadratic fiber formula is exactly
\[
 \tau_p(a)-(21+\varepsilon_a+\delta_a)p+8.
\]
The strict transforms of the sixteen boundary curves together with the
split ADE exceptional curves have an explicit intersection graph.  A
nineteen-class principal minor has determinant \(48\); adjoining one central
\(A_2\) component at \(a=1\), or the open exceptional curve when
\(a^2-34a+1=0\), gives determinants \(-32\) and \(-96\), respectively.
Thus they contain a nondegenerate Frobenius-fixed subspace of rank
\(19+\varepsilon_a+\delta_a\).  Removing it leaves a weight-two trace
\(\theta_p(a)\) of dimension \(3-\varepsilon_a-\delta_a\), and
\[
 \mu_p(a)=p^2-6p+12+\theta_p(a),\qquad
 \mu_p(a)-(p-2)^2=\theta_p(a)-2p+8.
\]
Deligne therefore gives the uniform bound \(5p+8\), improving the termwise
elliptic bound \(O(p^{3/2})\) to \(O(p)\).  This still does not control the vanishing
defining-characteristic Mellin coefficients: pointwise trace bounds allow
exact cancellation modulo \(p\).  The full compactification and correction
ledger are in `toric_fiber_k3.tex`; the standard-library verifier checks all
fibers and singularities for the 88 pairs \((p,a)\) with \(5\le p\le23\),
as well as the three exact lattice minors.  The full-paper compile and the
independent tmux-11 seven-point audit pass.

The optimized exact census in
\`primitive_projective_prime_scan.cpp\` groups both raw and
reflection-quotient statistics in \(O(p)\) operations per prime; its exact
pair decomposition takes an additional \(O(E_p)\) operations.  An exact
16-shard census through \(p\le500000\) checks all \(41{,}535\) primes in
the range and finds only \(20\) raw short chains, forming ten reflected
phantom pairs.  Every quotient state has count at most one, there is no
nonreflection same-state pair, and \(\bar A_p=0\) throughout.  The local
single-process implementation independently reproduced the same 20 records
and aggregate, and `primitive_projective_prime_verify.py` recomputed all 20
from fixed inputs using a separate recurrence and carrier check.  The local
scanner independently enforces the first raw pair
\[
 p=1297,\quad q=454,\quad
 (360;5,20,10)\longleftrightarrow(901;10,20,5).
\]
It also decomposes the quotient energy exactly into the diagonal, overlapping,
and separated parts and asserts that the three contributions sum to the full
energy; both distinct-pair contributions vanish throughout this census.
The random-fiber parameter count predicts one quotient chain with mean
\(p^{-1/2+o(1)}\).  Two overlapping windows can share three returns and
have heuristic mass \(p^{-1+o(1)}\), but their total energy is already
covered by the six-partner lemma.  Two separated chain orbits in one state
have heuristic mass \(p^{-2+o(1)}\), summable over primes.  This makes a
summability theorem for the separated energy a sharper target than
pointwise uniqueness of all chain orbits.  The census is unusually strong
evidence for projective dispersion, but remains a finite certificate;
weak-DPLS is not a consequence of the existing local algebra.

A separate toric calculation closes the bounded-degree interpolation
shortcut.  Let \(d_p\) be the degree of the reduced polynomial representing
\(r\mapsto b_r\) on \(\mathbb F_p\).  The fixed Laurent polynomial
\(\Lambda\) in \`oracleC_result.tex\` has nonnegative integral coefficients,
\(\operatorname{CT}\Lambda^n=b_n\), constant coefficient \(5\), and
\(\Lambda(1,1,1)=40\).  Hence \(\Lambda-1\) is still coefficientwise
nonnegative, with constant coefficient \(4\).  If \(d_p<p-1\), put
\(m=d_p+1\).  Finite differences give
\[
 p\mid\Delta^m b_0=\operatorname{CT}(\Lambda-1)^m,
 \qquad
 0<4^m\le\operatorname{CT}(\Lambda-1)^m\le39^m,
\]
and therefore
\[
 d_p\ge\frac{\log p}{\log39}-1.
\]
This rules out bounded interpolation degree but is only logarithmic and
gives no useful upper bound for the zero set.  The square/twisted-square
Hasse factorization is also terminal as a purely formal input.  If
\(B(t)=1+x_1t+\cdots+x_{d-1}t^{d-1}+t^d\) is reciprocal and
\(H=\Delta^\varepsilon B^2\), then for every \(m<d/2\)
\[
 H_j=2x_j+Q_j(x_1,\ldots,x_{j-1})\qquad(1\le j\le m).
\]
In odd characteristic this triangular map is a polynomial automorphism, so
one can prescribe an arbitrary clustered zero pattern in the first \(m\)
coefficients; reciprocity supplies its reflected copy.  The proof is in
`hasse_square_no_go.tex`.  Thus squareness and reciprocity alone give no
short-window exclusion.  Any use of the actual Hasse square root needs new
arithmetic information distinguishing its coefficient vector from the full
reciprocal-square parameter space.

There is an equivalent marked arithmetic observable.  The endpoint identity
\[
 N_{x+1}(-1)=(x!)^3b_x
\]
shows that, for \(p>x+a+b+c\), an actual four-zero tuple at gaps
\((a,b,c)\) is exactly detected by the large-prime part of
\[
 \gcd\!\left(
 N_{x+1}(-1),\,N_a(x),\,N_b(x+a),\,N_c(x+a+b)
 \right).
\]
The full addition identity is
\[
 N_{x+1+h}(-1)
 =N_{x+2}(-1)N_h(x)
 -(x+1)^6N_{x+1}(-1)N_{h-1}(x+1).
\]
After imposing \(p\mid N_{x+1}(-1)\), the first factor on the right is
a unit because consecutive Ap\'ery zeros are impossible.  Hence this
addition law proves the exact support equivalence but supplies no additional
congruence or height saving.
This removes every phantom before elimination, but the first entry has depth
\(x+1\asymp p\).  A direct height estimate therefore does not save a gap
parameter.  The missing arithmetic input can be stated as a uniform
\(1/X\)-thinning of the fully deflated short-chain content after imposing
this long endpoint divisor; it is the anchored moving-gcd counterpart of
weak-DPLS.

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
