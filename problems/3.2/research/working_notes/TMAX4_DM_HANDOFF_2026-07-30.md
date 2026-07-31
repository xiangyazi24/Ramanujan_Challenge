# P3.2 authority handoff to TMax4 `dm`

Date: 2026-07-30  
From: Codex `zinan:3` (`life`)  
To: Claude `zinan:4` (`dm`)

## Objective and authoritative source

Prove, unconditionally and for every \(n\),
\[
\log G_n=o(n),\qquad
G_n=\gcd(d_n^3a_n,d_n^3b_n)
\]
for the Apéry \(\zeta(3)\) sequences.

The authoritative accumulated analysis is:

`problems/3.2/research/working_notes/Q32_SEPARATION_ANALYSIS.md`

In particular, read the executive ledger and Sections 26, 30, and
68.19--68.23.  Do not rely on the older `dmm` prose without rechecking it.

## Exact remaining problem

The small-prime and outer-prime parts are proved.  After the Casoratian,
Lucas, and companion reductions, it is enough to prove
\[
\sum_{\substack{n/2<p\le n\\p\mid b_n}}\log p=o(n).
\]
Equivalently,
\[
\log\operatorname{rad}\gcd(b_n,\operatorname{lcm}(1,\ldots,n))=o(n).
\]
Write
\[
{\cal T}_n=\{p:n/2<p\le n,\ p\mid b_n\}.
\]
For \(p=n-r\), Apéry--Lucas gives
\[
p\mid b_n\iff p\mid b_r.
\]
No published standard theorem found so far proves this moving
large-prime radical estimate.

## New terminal carrier: exact theorem, not yet a height bound

Put \(M=n-1\),
\[
d_0=\lfloor M/2\rfloor+1,\quad L_0=M-d_0,\quad
f_L=G_{M-L,L}(C_M),\quad F_j=f_{L_0-j}.
\]
For \(1\le j<L_0\), with \(L=L_0-j\), define
\[
{\cal E}_j=
\frac{F_{j-1}F_{j+1}-F_j^2}
{\gcd\{\binom nL,\binom n{L+1}\}}.
\]
This is integral.  Every target
\[
p>d_0+K+1
\]
divides \(\gcd({\cal E}_1,\ldots,{\cal E}_K)\).  Thus \(K=2\) misses
only three integers immediately above the half boundary.

Exact hostile rows:
\[
\begin{array}{c|c}
n&\gcd({\cal E}_1,{\cal E}_2)\\ \hline
200&2^2\cdot5\cdot139\cdot181\\
272&2\cdot191\cdot233\\
300&11^2\cdot191\cdot227\\
321&179\cdot193\cdot211.
\end{array}
\]
The individual \({\cal E}_j\)'s have thousands of bits.  Therefore a
proof of
\[
\log\gcd({\cal E}_1,{\cal E}_2)=o(n)
\]
would close the top-half channel up to \(O(\log n)\).

Reproducible files:

- `problems/3.2/research/working_notes/Q32_TERMINAL_TURAN_HANKEL.md`
- `problems/3.2/research/scripts/q32_terminal_turan_hankel_audit.py`

The audit now also checks the rank-one formula for every safe candidate,
not only targets.  A dense run \(12\le n\le60\) checked 829 cases
(33 target, 796 non-target), all exactly.

## Sharp scalar no-go

Let
\[
a_j=F_{j-1}-F_j,\quad
A_j=F_{j-1}F_{j+1}-F_j^2,\quad A_j=q_j{\cal E}_j.
\]
Then
\[
A_j=F_j(a_j-a_{j+1})-a_ja_{j+1}
\]
and
\[
(a_{j+1}-a_{j+2})A_j-(a_j-a_{j+1})A_{j+1}
=a_{j+1}(a_ja_{j+2}-a_{j+1}^2).
\]
At every common candidate prime, write
\[
a_j=p\alpha_j,\qquad q_j=pu_j.
\]
All \(F_j\) have the same marked residue, and
\[
{\cal E}_j\equiv
F_0u_j^{-1}(\alpha_j-\alpha_{j+1})\pmod p.
\]
Hence all scalar Turán eliminations remain rank one: eliminating \(F_0\)
reintroduces the entire candidate primorial of linear logarithmic height.
Larger Hankel determinants give no independent equation by
Desnanot--Jacobi.

At the moving onset \(L=r\), the once-divided quotient is exactly
\[
\frac{\widetilde{\cal E}_r}{p}
\equiv(-1)^{r+1}(r+1)B_r\,\frac{f_{r+1}}p\pmod p.
\]
Thus the surviving coordinate is the genuine Newton--Cartier/Fermat
ghost digit \(f_{r+1}/p\), not a fixed polynomial in \(n,p,r\).

Cross-\(n\) residue checks at \(n=200,272,300,321\) found no simple
relation between the two Turán residues at \(n,n-1,\ldots,n-4\).
The prefix gcd itself is nevertheless striking:

- \(n=200\): \(\gcd(F_0,F_1)=5\cdot139\cdot181\);
- \(n=272\): \(191\cdot233\);
- \(n=300\): \(11\cdot191\cdot227\);
- \(n=321\): \(179\cdot193\cdot211\).

Adding more \(F_j\)'s usually does not change that small gcd.  This is
evidence only; the missing result is its characteristic-zero height
bound.

### New post-selector seam reduction

The final cross-\(n\) audit produced one further exact positive result.
Its complete report and verifier are:

- `/tmp/P32_TERMINAL_CROSS_N_FINAL.md`
- `problems/3.2/research/scripts/q32_terminal_cross_n_audit.py`

Define
\[
J_m(i)=\sum_{\kappa\in P\cap\mathbb Z^3}
\lambda_\kappa c_m(i\kappa)
=\operatorname{CT}\{\Lambda(X)^m\Lambda(X^{-i})\}.
\]
For \(M=p+r-1\) and \(0\le i\le r-1\), Freshman's dream plus the
polytope support gives the exact congruence
\[
C_{p+r-1}(p+i)\equiv J_{r-1}(i)\pmod p.
\]
If \(s=p-1-d_0\), then \(F_s^{(n)}\) is the last terminal window
containing the marked node \(p-1\), and the first post-seam value is
\[
F_{s+1}^{(n)}\equiv S_r\pmod p,
\qquad
S_r=\sum_{i=0}^{r-1}(-1)^i\binom r{i+1}J_{r-1}(i).
\]
Equivalently,
\[
S_r=b_r-\sum_\kappa\lambda_\kappa\operatorname{CT}
\left[
\Lambda^{r-1}X^{-(r-1)\kappa}(X^\kappa-1)^r
\right].
\]
Thus the uncontrolled \(n\)-scale seam has been reduced to the
lower-index pair
\[
\boxed{\gcd(b_r,S_r)}.
\]
The proof ingredients were checked in 1,595 exact coefficient cases
for every \(5\le p\le31\), \(1\le r<p\), \(0\le i<r\), and on all nine
hostile targets.  For \(1\le r\le100\), the exact gcd has no prime
factor outside \(\{5,11,19\}\), and its maximum is \(55\).

This is not yet closure.  A theorem
\[
\log\operatorname{rad}\gcd(b_r,S_r)=o(r)
\]
(in particular, a fixed resultant/content bound) would supply the
post-seam boundary observability datum needed by the finite-rank
terminal propagation route.  The immediate task is to compute and
compare the Ore/Picard--Fuchs modules of \(b_r\) and \(S_r\), and test
whether the observed \(\{5,11,19\}\) support comes from a fixed
resultant.

The naive cross-\(n\) \(2\times2\) Casoratians do not work: every
candidate prime divides them, and every tested target has valuation
exactly one.  Dividing the universal factor removes the target too.

## Actual Fermat ghost: newest exact local result

See:

- `/tmp/P32_ACTUAL_FERMAT_GHOST_RANK_AUDIT.md`
- `/tmp/p32_actual_fermat_hessian_audit.sage`

In the smooth critical algebra
\[
\mathbb Z[1/2,z]/(2z^2-1),\qquad
\lambda=17+24z=17+12\sqrt2,
\]
the scalar Fermat ghost value can vanish: exact quadratic-Wieferich
examples occur at \(p=13,31\).  Therefore any proof using only the
critical value is false.

The full local two-jet is stronger:
\[
\Phi_p\equiv q_{p,\chi}(\lambda)
+\lambda^{p-1}(\Lambda-\lambda)\pmod{(p,\mathfrak m^3)}.
\]
Moreover
\[
\det\operatorname{Hess}_\theta\Lambda=4u\lambda^2
\]
has fixed norm \(-8\).  Hence the quadratic ghost jet is uniformly
nondegenerate for every odd \(p\), including \(13,31\).  An exact
central \(p=31\) computation has actual rank two.  This identifies the
quadratic/conormal Hessian layer as the only currently credible source
of a second direction.

However, targetness has not yet been shown to annihilate an independent
Hessian observation.  The precise positive task is:

> Prove a target-preserving transversality identity between the
> nondegenerate quadratic marked jet and an actual boundary/terminal
> packet, with only fixed-prime denominators and subexponential
> characteristic-zero height.

The scalar Turán quotient does not do this: every hostile target has
\(v_p({\cal E}_j)=1\), so \({\cal E}_j/p\) is a unit.

## Global ghost/norm routes already ruled out

See `/tmp/P32_ONSET_GLOBAL_HEIGHT_AUDIT.md`.

1. Every fixed \(A_n(P,R)\in\mathbb Z[P,R]\), under
   \(P=p,R=n-p\), reduces modulo \(p\) to \(A_n(0,n)\).  Moving points
   collapse, so fixed-degree root counting gives nothing.
2. The actual raw Fermat ghosts have a triangular boundary coefficient
   of size \(\exp(\Omega(p))\).  Nontrivial raw linear combinations and
   products have logarithmic height \(\Omega(n)\).
3. The literal cyclotomic norm vanishes for every odd candidate because
   \(\Lambda(-1,y,z)=0\); it is target-blind.
4. Clearing denominators in a product of onset quotients reintroduces
   exactly the unknown target radical.

Do not spend time on raw ghost products, fixed polynomials in \((p,r)\),
or literal root-of-unity norms unless a genuinely new Apéry-specific
functional cancels the exponential boundary while preserving the onset
digit.

## Other proved reductions and no-go results

- Casoratian: no primes \(p>n\); primes \(p\le\sqrt n\) contribute
  \(O(\sqrt n)\).
- Density one is unconditional from \(|Z_p|\ll p^{2/3}\), but sparse
  exceptional columns remain.
- Pointwise \(|Z_p|=O(1)\) is not supported and is modularity-hard.
- Plain binomial windows and general factorial-ratio carriers have a
  rigid linear Chebyshev-height tradeoff.
- \(p^2/p^3\) two-coordinate block lifts are rank one on targets; the
  first transverse quotient is a unit, not another target zero.
- Boundary-blind Ore recurrences retain a cumulative integration
  constant; increasing recurrence order alone cannot remove it.
- Irrationality-measure arguments cannot yield \(o(n)\).

A sufficient averaged alternative remains
\[
H_2(N):=\sum_{N<n\le2N}|{\cal T}_n|(|{\cal T}_n|-1)
=o(N^2/\log^2N).
\]
Numerically \(H_2(N)\asymp N/\log^2N\), but no proof is known.

## Immediate recommended continuation

1. Start from the exact Hessian two-jet report, not from the scalar ghost.
2. Express the quadratic/conormal jet in the actual terminal
   constant-term packet and test whether targetness kills a second
   divided observation.
3. If that fails, preserve the failure certificate and return to the
   weak global second-moment criterion; do not reopen saturated
   scalar/Pascal routes.
4. A Mathematica 15 installation is available at
   `/opt/homebrew/bin/wolframscript`.  It is useful for recurrence and
   elimination discovery, but every claim must be converted to an exact
   replayable identity.

The completed Mathematica diagnostic is
`/tmp/P32_MATHEMATICA_TERMINAL_RECURRENCE.md`.  On the exact parity
subsequences through \(n=100\), constant-recurrence Hankel matrices are
full rank through order \(10\).  Modular polynomial-coefficient
recurrence ansätze through the reported order/degree boxes (including a
common \(F_0,F_1\) operator through all 71 decidable pairs with
order/degree at most \(8\)) also have full column rank modulo
\(10^9+7\).  The exhibited nonzero modular minors rigorously rule out
only those finite ansätze; they do not rule out a larger finite
\(\Lambda^2\) packet transfer.  A scalar low-order parity recurrence is
therefore not the next route.

## Repository state

The worktree contains many unrelated user changes.  Do not stage them.
At handoff time `main` is one commit ahead and four behind `origin/main`.
Stage only P3.2 files after reconciling upstream.  No full proof has been
obtained, so no Telegram breakthrough message and no formalization have
been started.
