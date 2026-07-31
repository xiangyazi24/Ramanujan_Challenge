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
with the empirical input that the exact gcd is uniformly tiny.  **Caution:** the
recorded claim that its prime support is \(\{5,11,19\}\) holds only for
\(r\le100\); see the correction in \S3b, where the support is shown to grow.

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
   singularities and its local exponents are computable.  (The original motivation —
   that \(\{5,11,19\}\) would show up as a discriminant/resultant of that finite
   family — has since been refuted: see \S3b for the growing support and the absence
   of any corresponding rank drop.)

## 3b. The raywise Smith content is trivial: it is exactly the \(y\leftrightarrow z\) symmetry

Form the integer matrix whose rows are consecutive \(r\) and whose \(22\) columns are the
\(21\) pieces \(U_\kappa(r)\) together with \(b_r\).  For a \(26\times22\) window starting at
\(r=2\):
\[
 \operatorname{rank}_{\mathbb F_q}=15\quad\text{for every prime }5\le q\le113,
\]
with the only smaller values at \(q=2\) (rank \(3\)) and \(q=3\) (rank \(14\)).  The
deficiency \(22-15=7\) is therefore a **characteristic-zero** deficiency, not a
prime-specific one, and it is completely explained:

\(\Lambda\) is invariant under \(\sigma:y\leftrightarrow z\) (checked: \(\lambda_{\sigma\kappa}
=\lambda_\kappa\) for every lattice point), hence
\[
 \boxed{\;U_{\sigma\kappa}(r)=U_\kappa(r)\quad\text{for all }r.\;}
\]
Of the \(21\) nonzero rays, exactly \(7\) are \(\sigma\)-fixed and \(7\) form swapped pairs;
the pair identities were verified at \(r=5,9,13\) with no violation.  So the rank is
\(22-7=15\) on the nose.

Two consequences.

1. **No hidden content.**  The exceptional primes \(\{5,11,19\}\) of \(\gcd(b_r,S_r)\) are
   *not* elementary divisors of the raywise value matrix — there is no rank drop at
   \(5,11,19\) beyond the universal one.  The "fixed resultant/Smith content" guess has to
   be tested on the operators, not on the values.
2. **The family halves.**  Only \(14\) distinct pieces occur (\(7\) fixed rays and \(7\)
   pair representatives with doubled weight), which halves every subsequent raywise
   computation.

Separately, two mechanisms that would have been convenient are ruled out by exact
computation:

- the Casoratian-style determinant \(W_r=b_rS_{r+1}-b_{r+1}S_r\) has **no** small
  factorisation — for \(r\le16\) it carries big prime factors
  (\(9.5\cdot10^7\), \(6.6\cdot10^{12}\), \(1.06\cdot10^{19}\), \(3.6\cdot10^{45}\), …), so
  \((b,S)\) is not a pair with an explicit hypergeometric Casoratian;
- the Apéry-operator residual
  \(R_r=(r+1)^3S_{r+1}-(34r^3+51r^2+27r+5)S_r+r^3S_{r-1}\) is likewise not simple:
  \(|R_r|\approx|b_r|\cdot e^{0.19r}\) with large prime factors, so applying \(L_b\) to
  \(S\) does not produce a hypergeometric right-hand side.

### Correction to the recorded \(\gcd(b_r,S_r)\) data

The statement "for \(1\le r\le100\), \(\gcd(b_r,S_r)\) has no prime factor outside
\(\{5,11,19\}\), maximum \(55\)" is **an artifact of the range**.  An independent exact
recomputation through \(r=300\) (78 nontrivial values) gives

\[
\begin{array}{c|l}
 r\text{-window}&\text{primes occurring in }\gcd(b_r,S_r)\\ \hline
 [1,100)&5,\ 11,\ 19\\
 [100,200)&5,\ 19\\
 [200,300)&5,\ 11,\ 17,\ 19,\ 31,\ 37,\ 61
\end{array}
\]

with maximum \(\gcd=305=5\cdot61\) at \(r=257\) (and \(125\) at \(r=171,286\)).  So the
support is **not** a fixed finite set, and the "fixed resultant / fixed discriminant"
reading of the \(\{5,11,19\}\) data is not supported.  Together with §3b (no rank drop at
\(5,11,19\)), the fixed-content hypothesis should be considered dead.

What does survive, and is the statement the application actually needs (there the relevant
prime satisfies \(p>r\)):

> Over \(1\le r\le300\), **every** prime factor of \(\gcd(b_r,S_r)\) is \(<r\), the only
> exception being the trivial \(r=1\).  The gcd itself never exceeds \(305\).

So the empirical picture is a slowly growing support with uniformly tiny total size —
consistent with \(\log\gcd(b_r,S_r)=O(\log r)\), but with no fixed-resultant mechanism
behind it.

## 3c. The marked scalar is a moment of point counts

The one mod-\(p\) scalar that the whole terminal programme keeps aliasing has an exact
exponential-sum form.  For every prime \(p\ge5\) and every moment \(M\),
\[
 \boxed{\quad
 C_M(p-1)\equiv-\sum_{x,y,z\in\mathbb F_p^\times}\Lambda(x,y,z)^M\pmod p.
 \quad}
\]
Equivalently, with \(N_p(t)=\#\{(x,y,z)\in(\mathbb F_p^\times)^3:\Lambda=t\}\) and
\(r=M\bmod(p-1)\),
\[
 C_M(p-1)\equiv-\sum_{t\in\mathbb F_p^\times}t^{\,r}N_p(t)\pmod p,
\]
so that, with \(n=p+r\) (whence \(M=n-1\equiv r\bmod p-1\)), **the target condition
\(p\mid F_0\) is precisely the vanishing modulo \(p\) of the \(r\)-th moment of the
point-count function of the Apéry family**, in agreement with Apéry--Lucas
(\(p\mid b_n\iff p\mid b_r\)).

Audited by `problems/3.2/research/scripts/q32_marked_scalar_character_sum.py`:
95 shell-vs-exponential-sum checks, 95 shell-vs-moment checks and 44
shell-vs-Apéry checks, all exact.

Two things this buys.

1. **It explains the rank-one wall structurally.**  Modulo \(p\) the only surviving datum
   of the family is the one-dimensional unit-root/Hasse--Witt invariant; every mod-\(p\)
   observable built from the coefficient array is therefore a multiple of that single
   scalar.  This is why no linear elimination inside any terminal family produced a second
   condition, and it predicts that a genuinely independent condition needs the mod-\(p^2\)
   datum, where the full Frobenius matrix and not just its unit root appears.
2. **It puts the target condition in standard language.**  The column problem becomes:
   for fixed \(n\), bound
   \(\#\{p\in(n/2,n]:\sum_{t}t^{\,n-p}N_p(t)\equiv0\ (p)\}\).
   The family is the modular one (weight-4 level-8 newform 8.4.a.a), so this is where
   Hasse--Witt/unit-root and Lang--Trotter-type technology, rather than Pascal elimination,
   would have to act.

## 3d. The M-direction Casoratian route fails on height

With the Pascal-normalised carrier
\(\widetilde{\cal E}_j(M)=(k+1)\bigl(F_{j-1}F_{j+1}-F_j^2\bigr)/\binom nk\),
\(k=L_0-j\) (integrality confirmed at every tested level), the \(2\times2\) step-two
Casoratian
\[
 W(M)=\widetilde{\cal E}_1(M)\widetilde{\cal E}_2(M+2)
      -\widetilde{\cal E}_1(M+2)\widetilde{\cal E}_2(M)
\]
is a genuine carrier: every base-level target divides it (verified at the only levels with
targets in the scanned window, \(n=54,56,68\)).  But
\[
 \log_2|W(M)|\approx21\,n,
\]
i.e. the height is \(\Theta(n)\), and \(W\) admits **no** first-order rational ratio
\(Q(M)W(M+2)=P(M)W(M)\) with \(\deg\le12\) on either parity class.  So \(W\) is not
hypergeometric and no "balanced ratio" (equal degrees and equal leading absolute values,
which would have given \(\log|W|=O(\log M)\)) can exist.  The route is closed at the
\(2\times2\) level, and larger \(K\times K\) Casoratians only increase the height.

## 3e. The target count itself, over 20 000 indices

The reduction ends at the top-half sum
\[
 T(n)=\sum_{\substack{n/2<p\le n\\ p\mid b_n}}\log p,
 \qquad
 K(n)=\#\{p\in(n/2,n]:p\mid b_n\},
\]
(using Apéry--Lucas, \(p\mid b_n\iff p\mid b_{n-p}\)), and the open statement is
\(T(n)=o(n)\).  Computing \(Z_p\) for every prime \(p\le N\) by iterating the recurrence
modulo \(p\) and assembling gives, for \(N=60\,000\):

\[
\boxed{\;K(n)\le3\ \text{ for every }n\le60\,000.\;}
\]
The running maximum increases only three times over the whole range: \(K=1\) first at
\(n=6\), \(K=2\) first at \(n=200\), \(K=3\) first at \(n=321\), and never \(4\).
Exactly \(144\) indices have \(K\ge2\) and exactly \(7\) have \(K=3\).  The mean is
\[
 \overline{K}=0.0733,
\]
against the Mertens/Poisson prediction
\(\sum_{n/2<p\le n}p^{-1}\cdot\mathbb E|Z_p|\approx\log2/\log n=0.075\) at this scale —
a quantitative confirmation of the Poisson model for the actual target count, not merely
for \(|Z_p|\).

So the empirical picture is that \(K(n)=O(1)\) on the tested range, hence
\(T(n)=O(\log n)\), vastly stronger than the required \(o(n)\); and the difficulty is
entirely that the trivial bound \(K(n)\le\log b_n/\log(n/2)=O(n/\log n)\) is the best one
provable unconditionally.  Script:
`problems/3.2/research/scripts/q32_top_window_target_counts.py`.

## 3f. Empirical audit of the reduction chain: where \(G_n\) actually lives

Computing the headline quantity directly (rational Apéry recursion, \(d_n=\mathrm{lcm}(1..n)^3\),
`q32_actual_Gn_audit.py`) for \(n\le330\), and splitting
\[
 G_n=\gcd(d_na_n,d_nb_n)=\frac{d_n}{D_n}\cdot\gcd(A_n,b_n),
 \qquad a_n=\frac{A_n}{D_n}\ \text{in lowest terms},
\]
gives two facts worth recording.

1. **The intrinsic numerator gcd is trivial**:
   \[
    \gcd(A_n,b_n)=1\quad\text{for every }n\le330 .
   \]
   So all of \(G_n\) is the over-clearance factor \(d_n/D_n\); the classical denominator
   envelope \(d_n\) is *not* wasteful in the direction that would break the reduction.
2. **The top-window targets are exactly the large primes of \(G_n\)**.  At the two indices
   in range with more than one target,
   \[
   \begin{aligned}
    G_{200}&=2\cdot3^3\cdot5\cdot17\cdot19\cdot\mathbf{139}\cdot\mathbf{181}
      &&(32\text{ bits}),\\
    G_{321}&=2\cdot3\cdot5\cdot7^2\cdot\mathbf{179}\cdot\mathbf{193}\cdot\mathbf{211}
      &&(34\text{ bits}),
   \end{aligned}
   \]
   and the bold primes are precisely the targets \(p\in(n/2,n]\) with \(p\mid b_n\).
   Every one of them divides \(G_n\).

This is an independent confirmation that the chain
"\(G_n\rightarrow\) top-window radical of \(b_n\)" is attacking the right object, and it
localises the difficulty sharply: over the whole tested range \(G_n\) never exceeds
\(34\) bits, its large part is exactly the \(K(n)\le3\) targets, and its small part is
a handful of tiny primes.

## 3g. The reflection law and the parity of \(|Z_p|\), derived

The moment identity of \S3c explains, rather than observes, the reflection phenomenon the
notes had recorded empirically.

**Step 1 (geometry).**  The point count is inversion-symmetric:
\[
 N_p(t)=N_p(t^{-1})\qquad\text{for every }t\in\mathbb F_p^\times
\]
(verified for all \(p\le23\): 420 checks, no failures).

**Step 2 (moments).**  Hence
\(\sum_t t^{\,r}N_p(t)=\sum_t t^{-r}N_p(t)\), and since \(t^{-r}=t^{\,p-1-r}\),
\S3c gives the **palindromy congruence**
\[
 \boxed{\;b_{\,p-1-r}\equiv b_r\pmod p\quad\text{for all }0\le r<p-1.\;}
\]
Verified exactly (ratio \(b_{p-1-r}/b_r\equiv1\), never merely \(\pm1\)) for every
\(p\le23\) — 420 checks — and spot-checked at \(p=101,211\).

**Step 3 (zero sets).**  Therefore \(Z_p\) is stable under the involution
\(r\mapsto p-1-r\).  Checked on **all 548 primes \(5\le p<4000\)**: no failure.
Consequently
\[
 |Z_p|\ \text{is even, unless the fixed point }\tfrac{p-1}2\in Z_p .
\]
Over \(p<4000\) exactly two primes have odd \(|Z_p|\), namely \(p=11\) with
\(Z_{11}=\{5\}\) and \(p=3137\) with \(Z_{3137}=\{1568\}\) — in both cases the zero **is**
the fixed point \((p-1)/2\).  This identifies the "central exceptions \(p=11,3137\)"
recorded in the executive ledger as exactly the fixed-point case, and it explains why the
pair count \(|Z_p|/2\) — not \(|Z_p|\) — is the quantity that follows a Poisson law.

Audit: `q32_marked_scalar_character_sum.py` (420 inversion checks, 420 palindromy checks,
zero-set involution checks, and the odd-count classification).

## 4. Reproducibility

- `problems/3.2/research/scripts/q32_seam_ray_split_audit.py` — the exact
  audit (polynomial form vs. binomial form vs. reference seam scalar).
- The guessing runs are plain modular nullspace computations; the generator
  for \(U_\kappa(r)\) is formula (4) with cached binomial rows.
