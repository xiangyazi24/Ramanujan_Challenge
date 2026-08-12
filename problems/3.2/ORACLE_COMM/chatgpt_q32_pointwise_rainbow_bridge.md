# Q7697 — pointwise rainbow bridge for the all-index Apéry theorem

## Executive verdict

There is a clean pointwise bridge, but it is **not** the existing row-aggregate
BFH statement and it is **not** the independence-scale factorial moment
`(HM)_3`.

The right intermediate statement is a **short parabolic centered-variance / early-CRT
covariance bound**.  Put

\[
 Z_p=\{0\le r<p:p\mid b_r\},\qquad
 K_P(m)=\#\{P<p\le2P:m\bmod p\in Z_p\},
\]

\[
 \lambda_P=\sum_{P<p\le2P}\frac{|Z_p|}{p},\qquad
 M=2LP,\qquad
 \mu_{P,L}=\frac1M\sum_{0\le m<M}K_P(m).
\]

The banked Apéry bound is

\[
 \lambda_P\ll \frac{P^{2/3}}{\log P}.
\]

The precise sufficient theorem is the following.

> **Short parabolic rainbow variance theorem `PV_eta`.**  There are fixed
> `eta>0` and `C` such that, for every sufficiently large `P` and every integer
> `1<=L<=(log P)^20`, one has
> \[
> V(P,L):=\sum_{0\le m<M}(K_P(m)-\mu_{P,L})^2
> \ll M\lambda_P P^{1/3-\eta}(\log P)^C.                 \tag{PV_eta}
> \]

This theorem is a genuinely cross-prime statement and it includes both the
transverse and nontransverse coefficient-zero channels.  It implies
`W_n=o(n)` **for every n**.

I do **not** find a proof of `(PV_eta)` from the currently banked Apéry
identities.  The exact obstruction is now sharply localized: one must prove
anti-alignment of the actual coefficient-zero sets in *different
characteristics* against the initial interval of CRT representatives.  Lucas,
reflection, the one-sided carrier, the level-six Hasse square, and the current
transverse BFH content each stop before that mixed-characteristic step.

The smallest unconditional theorem that survives is a hard-core reduction:
after removing small primes, small reflected rows, and a polylogarithmic
quotient-exception carrier, every remaining target is a fresh prime divisor of
the single integer `b_n`.  This reduces the pointwise problem to a moving-anchor
rainbow divisor set and cleanly separates its transverse and nontransverse
parts, but gives no sublinear estimate by itself.

---

## 1. Why `(PV_eta)` implies the full pointwise theorem

For each `p` in `(P,2P]`, let

\[
 A_p=\#\{0\le m<M:m\bmod p\in Z_p\}.
\]

Periodicity gives

\[
 A_p=M\frac{|Z_p|}{p}+O(|Z_p|).
\]

Since `p<=2P`,

\[
 \sum_{P<p\le2P}|Z_p|\le 2P\lambda_P,
\]

and therefore, because `M=2LP`,

\[
 |\mu_{P,L}-\lambda_P|
 \le \frac{2P\lambda_P}{M}
 =\frac{\lambda_P}{L}.                                  \tag{1}
\]

Thus `mu_{P,L}<=2 lambda_P` for `L>=1`.

Now `(PV_eta)` bounds **every individual** point in `[0,M)`, not merely a
density-one subset: for each fixed `m<M`,

\[
 |K_P(m)-\mu_{P,L}|^2\le V(P,L).
\]

Hence

\[
 K_P(m)
 \ll \lambda_P+
 \sqrt{M\lambda_P P^{1/3-\eta}(\log P)^C}.
\]

Using `lambda_P << P^{2/3}/log P` and `M=2LP`,

\[
 K_P(m)
 \ll \frac{P^{2/3}}{\log P}
 +P^{1-\eta/2}\,L^{1/2}(\log P)^{(C-1)/2}.               \tag{2}
\]

So the logarithmic prime weight contributed by one dyadic block is

\[
 \sum_{\substack{P<p\le2P\\p\mid b_{m\bmod p}}}\log p
 \ll P^{2/3}
 +P^{1-\eta/2}L^{1/2}(\log P)^{(C+1)/2}.                 \tag{3}
\]

Fix `n` and take

\[
 L_n=(\log n)^{10}.
\]

The primes `p<=n/L_n` have total logarithmic weight

\[
 \sum_{p\le n/L_n}\log p=O(n/L_n)=o(n)
\]

by Chebyshev.  Decompose the remaining interval `(n/L_n,n]` into dyadic
blocks `(P,2P]`.  Each such block has `P>n/(2L_n)`, and therefore

\[
 n<2L_nP=M,
\]

so `m=n` is an admissible point in `(PV_eta)`.  Also
`L_n<=(log P)^20` for all sufficiently large `n`.  Summing (3) over the
`O(log L_n)` blocks is a geometric sum dominated by its largest block:

\[
 \sum_{\substack{n/L_n<p\le n\\p\mid b_{n\bmod p}}}\log p
 \ll n^{2/3+o(1)}
 +n^{1-\eta/2+o(1)}L_n^{1/2}(\log n)^{(C+1)/2}
 =o(n).
\]

Adding the given `O(n^{2/3})` remainder proves `W_n=o(n)` pointwise.

### Exact exponent threshold

More generally, if one could prove

\[
 V(P,L)\ll M\lambda_P P^{\beta+o(1)},
\]

then the pointwise deviation has exponent

\[
 \frac{1+2/3+\beta}{2}=\frac56+\frac\beta2.
\]

Thus the sharp threshold for this route is

\[
 \boxed{\beta<1/3.}
\]

This is why the `P^{1/3}` boundary in `(PV_eta)` is the right quantitative
target.

---

## 2. The equivalent early-CRT correlation statement

For distinct primes `p,q` in `(P,2P]`, define

\[
 J_{p,q}(M)=\#\{0\le m<M:m\bmod p\in Z_p,\ m\bmod q\in Z_q\}.
\]

Expanding the variance gives the exact identity

\[
 V(P,L)=
 \sum_p\left(A_p-\frac{A_p^2}{M}\right)
 +2\sum_{p<q}\left(J_{p,q}(M)-\frac{A_pA_q}{M}\right).   \tag{4}
\]

The diagonal term is `O(M lambda_P)`.  Therefore it is enough to prove the
one-sided signed covariance estimate

\[
 \left(
 2\sum_{p<q}\left(J_{p,q}(M)-\frac{A_pA_q}{M}\right)
 \right)_+
 \ll M\lambda_PP^{1/3-\eta}(\log P)^C.                  \tag{EC_eta}
\]

In the intended range `L=(log P)^{O(1)}`, one has `M<<P^2`.  Consequently
`pq>M` for every distinct `p,q` in the block, so every residue pair
`(r,s) in Z_p x Z_q` has at most one representative below `M`.  Hence
`J_{p,q}(M)` is **exactly** the number of actual Apéry zero pairs whose unique
CRT representative lies unusually early.

This is the pointwise rainbow object.  The experiment below measures (4) and
`(EC_eta)` directly.

### Why this is not another wrong-quantifier row average

A row aggregate counts incidences essentially linearly.  It can be small while
one integer `m` carries many different primes.  By contrast, a high-load
integer contributes the square of its overload to `V(P,L)`.  Thus the global
variance estimate has a deterministic maximal consequence
`|K_P(m)-mu|<=sqrt(V)` for **every** `m`.  No exceptional-set/density-one
argument is used.

---

## 3. Why independence-scale `(HM)_3` does not close the pointwise theorem

The existing global high-moment framework uses `m<P^2` and

\[
 \sum_{m<P^2}(K_P(m))_k\ll P^{2+o(1)}\lambda_P^k.
\]

It gives

\[
 \max_{m<P^2}K_P(m)
 \ll P^{2/k+o(1)}\lambda_P.
\]

With `lambda_P<<P^{2/3}/log P`, `k=3` gives

\[
 \max K_P(m)\ll P^{4/3+o(1)}/\log P,
\]

which is worse than the trivial prime count.  The banked theorem correctly
needs `k>6` for the all-index conclusion.

Even shortening the interval to `M=2LP` does not make an **uncentered
factorial** third moment sufficient at the naive independence scale:

\[
 \sum_{m<M}(K_P(m))_3\ll M\lambda_P^3
 \quad\Longrightarrow\quad
 \max K_P(m)\ll M^{1/3}\lambda_P
 \asymp L^{1/3}\frac{P}{\log P},
\]

whose logarithmic weight is still `asymp L^{1/3}P`.

What helps is **centering**.  The mean `lambda_P` is already harmless; only a
large positive deviation must be ruled out.  A centered third-moment theorem
would also work, e.g.

\[
 \sum_{m<M}|K_P(m)-\mu|^3
 \ll M\lambda_P P^{4/3-\eta+o(1)}
\]

would yield a pointwise saving.  But `(PV_eta)` is the lower-moment and cleaner
target.

---

## 4. Attempted proof from actual Apéry identities

### 4.1 Lucas: exact reduction to fresh prime divisors of one `b_m`

For `m<M=2LP` and `p in (P,2P]`, write

\[
 m=ap+r,\qquad 0\le r<p.
\]

Then `0<=a<2L`.  The exact Apéry Lucas congruence gives

\[
 b_m\equiv b_a b_r\pmod p.                              \tag{5}
\]

Let

\[
 B_L=\prod_{0\le a<2L}b_a.
\]

For a prime `p` with `p not| B_L`, (5) gives the equivalence

\[
 p\mid b_{m\bmod p}\quad\Longleftrightarrow\quad p\mid b_m. \tag{6}
\]

This is a real structural gain: outside a small quotient-exception set, the
rainbow load is exactly the number of dyadic prime divisors of a **single**
Apéry integer `b_m`.

The exceptional primes are harmless in logarithmic weight.  From the standard
positive Apéry sum

\[
 b_a=\sum_{k=0}^a\binom ak^2\binom{a+k}{k}^2
\]

one has `b_a<=(a+1)64^a`, hence `log b_a=O(a)` and

\[
 \sum_{p\mid B_L}\log p\le\log B_L=O(L^2).              \tag{7}
\]

For polylogarithmic `L`, this is negligible.

#### Exact obstruction

Equation (6) does **not** bound concentration.  It only converts it to

\[
 \prod_{\substack{P<p\le2P\\p\mid b_{m\bmod p}\\p\nmid B_L}}p\mid b_m.
\]

The height estimate `log b_m=O(m)` permits `O(P)` logarithmic prime content at
a moving point `m`; on the short interval it permits even `O(LP)`.  Thus the
Lucas carrier is target-exact but height-borderline.  It reproduces the
one-sided global-carrier obstruction rather than beating it.

This is also why a fixed-width linear-form carrier cannot repair the argument:
its universal prime content has exponential height, so one simply repackages
the target without obtaining the needed `o(P)` block weight.

### 4.2 Reflection: useful endpoint stripping, no cross-prime coupling

The exact reflection law

\[
 b_r\equiv b_{p-1-r}\pmod p                              \tag{8}
\]

implies that a target with

\[
 j_p:=\min(r,p-1-r)\le R
\]

has `p|b_{j_p}`.  Therefore all such distinct primes have total weight

\[
 \ll\sum_{j\le R}\log b_j=O(R^2).                       \tag{9}
\]

This is valuable for the hard-core reduction below.  But (8) is a
same-characteristic symmetry; it imposes no relation between the early CRT
positions for `p` and `q`.  The anchored-rainbow obstruction in the existing
`hm3_result.tex` already shows that reflection, no-consecutive marks, local
interval bounds, and bounded row codegrees can coexist with severe
cross-prime bunching.

### 4.3 Generic CRT / parabolic large sieve: the `P^2` barrier

For `M<<P^2`, generic CRT says each residue pair contributes either zero or one
point below `M`.  That is exactly the problem, not an estimate: the desired
random proportion is `M/(pq)~L/P`, and generic CRT provides no reason the
actual Apéry pairs should realize that proportion.

Likewise, a standard large sieve on moduli of size `P` contains the geometric
term `M+P^2`; here `M=2LP<<P^2`, so the `P^2` term dominates.  Without an
Apéry-specific cancellation identity, this cannot prove `(EC_eta)`.

### 4.4 Level-six Hasse square / exact ODE: wrong object unless Mellin is mixed horizontally

The banked exact finite-field identity has the form

\[
 H_p(t)=\Delta(t)^{\varepsilon_p}B_p(t)^2,
 \qquad \Delta(t)=t^2-34t+1.
\]

But `Z_p` is a **coefficient** zero set, not the evaluation-root divisor of
`H_p`.  The exact finite Mellin relation is

\[
 b_r=-\sum_{t\in\mathbf F_p^\times}H_p(t)t^{-r}
 \qquad(1\le r\le p-2).                                 \tag{10}
\]

Thus the condition in `J_{p,q}(M)` asks for simultaneous vanishing of two
Mellin coefficients in two different characteristics, with their indices tied
by one integer `m` through `r=m mod p`, `s=m mod q`.

A Weil bound, the Hasse square, or the exact ODE at one prime controls only one
characteristic.  Taking a product over `p` and `q` tensorizes by CRT and leaves
the early-representative indicator untouched.  The missing theorem would have
to be a **mixed-characteristic horizontal Mellin correlation** whose diagonal
is strictly smaller than the occupancy-one atom term.  No such identity is
currently banked.

### 4.5 One-sided carrier / rank-one connection

The quotient carrier and Lucas (6) faithfully carry every fresh target, but
again their logarithmic height is linear in the moving index.  The previously
banked rank-one/full-period transfer phenomenon and SNF `diag(1,p^3)` explain
why the same connection does not automatically supply a second independent
mod-`p` direction.  A second target condition must come from genuinely new
mixed-characteristic information, not another presentation of the same
rank-one carrier.

---

## 5. Transverse versus nontransverse, and row versus pointwise

Use the terminology in the current Q7697 prompt:

\[
 T(p,r):\quad p\mid b_r\ \text{and}\ p\mid\Xi_r,
\]

\[
 NT(p,r):\quad p\mid b_r\ \text{and}\ p\nmid\Xi_r.
\]

Define `K_P^T(m)` and `K_P^NT(m)` by restricting the incidences in `K_P(m)`
to these two channels.  Then

\[
 K_P=K_P^T+K_P^{NT}.
\]

A channel-wise version of `(PV_eta)` for both centered loads would imply the
full theorem (by `(x+y)^2<=2x^2+2y^2`), but the current information is much
weaker:

* The BFH/transverse content concerns `T(p,r)` and is aggregated over row
  labels `r`.  Even a subquadratic row sum yields at most density-one-style
  information; it does not bound the maximal early-CRT load at a specified
  integer `m=n`.
* The nontransverse channel is genuinely nonempty; `(p,r)=(11,5)` is the
  explicit example in the prompt.  A proof that only treats `T` can therefore
  never be a proof of the full coefficient-zero statement.
* The covariance in (4) must include `T-T`, `T-NT`, and `NT-NT` pairs.  The
  mixed terms are not seen by a theorem whose hypothesis already requires
  `p|Xi_r` at both ends.

This is the precise quantifier/channel mismatch: BFH is a linear row aggregate
on a subrelation; `(EC_eta)` is a signed quadratic correlation on the full
relation, localized to the moving CRT graph generated by one point `m`.

---

## 6. Smallest unconditional theorem that survives

The following reduction uses actual Apéry identities and is pointwise.

> **Hard-core moving-anchor theorem.**  Let `L>=2` and `1<=R<sqrt(n)`.  For
> each target prime in the given lower-digit sum write
> `n=q_pp+r_p`, `0<=r_p<p`, and put
> `j_p=min(r_p,p-1-r_p)`.  Let
> \[
> B_L=\prod_{0\le a<2L}b_a
> \]
> and let `H_n(L,R)` be the set of primes satisfying
> \[
> \sqrt n<p\le n,\quad p>n/L,\quad p\mid b_{r_p},\quad
> j_p>R,\quad p\nmid B_L.
> \]
> Then
> \[
> W_n=
> \sum_{p\in H_n(L,R)}\log p
> +O\!\left(\frac nL+R^2+L^2+n^{2/3}\right).            \tag{11}
> \]
> Moreover, for every `p in H_n(L,R)`, one has `q_p<L` and
> \[
> p\mid b_{r_p}\quad\Longleftrightarrow\quad p\mid b_n. \tag{12}
> \]

### Proof

1. Primes `p<=n/L` have total weight `O(n/L)` by Chebyshev.
2. For `p>n/L`, `q_p=floor(n/p)<L`.
3. If `j_p<=R`, reflection (8) gives `p|b_{j_p}`.  The product of all
   distinct primes assigned to a fixed `j` divides `b_j`; summing over
   `j<=R` gives the `O(R^2)` bound (9).
4. The primes dividing `B_L` have total weight `O(L^2)` by (7).
5. On the remaining set, `q_p<L<2L` and `p not|B_L`, so Lucas gives (12).
6. Insert these removals into the exact lower-digit reduction from the prompt;
   its existing remainder contributes `O(n^{2/3})`.

Taking, for example,

\[
 L=(\log n)^{10},\qquad R=\frac{\sqrt n}{(\log n)^{10}},
\]

makes every error in (11) `o(n)`.  Therefore the entire all-index problem is
reduced to showing

\[
 \sum_{p\in H_n(L,R)}\log p=o(n),                       \tag{13}
\]

where every remaining prime is simultaneously

* a large fresh divisor of the moving integer `b_n`;
* represented by a row far from both reflected endpoints;
* attached to a polylogarithmically small quotient `q_p`;
* and either transverse or nontransverse in the `Xi` sense.

This is, in my view, the smallest honest surviving theorem.  It removes every
part that Lucas/reflection/height can actually remove, but does not rename the
remaining pointwise rainbow obstruction as an estimate.

---

## 7. Exact Sage experiment: measure the missing correlation, not zero counts

Files added:

* `problems/3.2/research/scripts/q7697_pointwise_rainbow_correlation.sage`
* `problems/3.2/research/scripts/q7697_pointwise_rainbow_verify.sage`
* `problems/3.2/research/scripts/q7697_pointwise_rainbow_reference.json`

The experiment recomputes `Z_p` exactly from the division-free recurrence

\[
 A_{r+1}=(34r^3+51r^2+27r+5)A_r-r^6A_{r-1},
 \qquad A_r=(r!)^3b_r,
\]

but it does **not** report a zero-set histogram as its target statistic.  It
constructs the full load vector `K_P(m)` for `m<2LP` and records:

1. exact `lambda_P` and exact finite mean `mu`;
2. exact centered variance `V(P,L)`;
3. exact centered absolute cubic moment;
4. exact `F_2=sum(K)_2` and `F_3=sum(K)_3`;
5. exact aggregate early-CRT pair discrepancy
   \[
   \sum_{p<q}\left(J_{p,q}(M)-A_pA_q/M\right);
   \]
6. the largest positive individual pair discrepancies;
7. every prime/quotient/residue incidence at points attaining `max K`;
8. the diagnostic ratios `V/(M lambda_P)` and
   `V/(M lambda_P P^{1/3})`.

The core correlation calculation is exact Sage arithmetic:

```sage
state = build_state(P, L)
M = state["M"]
A = state["A"]
K = state["K"]
mu = QQ(sum(K), M)
V = sum((QQ(k)-mu)^2 for k in K)

pair_excess = QQ(0)
for p, q in combinations(state["primes"], 2):
    J = len(set(state["positions"][p]).intersection(state["positions"][q]))
    pair_excess += QQ(J) - QQ(A[p]*A[q], M)

assert V == sum(QQ(A[p])-QQ(A[p]^2, M) for p in state["primes"]) \
            + 2*pair_excess
```

Run a scale by

```text
sage problems/3.2/research/scripts/q7697_pointwise_rainbow_correlation.sage \
  --P 1000 --L 8 --output /tmp/q7697-P1000-L8.json
```

and run the exact verifier by

```text
sage problems/3.2/research/scripts/q7697_pointwise_rainbow_verify.sage
```

The verifier independently enumerates CRT representatives, rather than merely
reusing the intersection calculation.  It checks:

* the banked exact zero sets `Z_7=[]`, `Z_11=[5]`, and
  `Z_181=[19,47,133,161]`;
* every pair `J_{p,q}(M)` against an independent CRT enumeration on the
  reference blocks;
* the exact variance decomposition (4);
* `F_2=2*(unordered pair count)`;
* and, on the smallest block, `F_3=6*(unordered triple CRT count)`.

No fitted random model is used in these assertions.  The JSON reference file
contains only exact verifier cases/known zero sets, not fabricated correlation
results.

### What would count as positive experimental evidence

The key quantity is not `|Z_p|`.  Across increasing `P` with polylogarithmic
`L`, inspect

\[
 \frac{V(P,L)}{M\lambda_P P^{1/3}}
\]

and the largest positive pair discrepancies.  Decay by a power of `P` would
match `(PV_eta)`.  Bounded `V/(M lambda_P)` would be much stronger (Poisson
scale).  Conversely, growth near `P^{1/3}` or a persistent family of large
positive early-CRT pair excesses identifies the exact place where the proposed
bridge fails.

---

## 8. Final theorem status

**Proved here:**

* `(PV_eta) => W_n=o(n)` pointwise, with the sharp variance-loss threshold
  `beta<1/3` checked explicitly;
* the exact early-CRT covariance identity (4);
* the Lucas fresh-divisor reduction (6) outside an `O(L^2)` logarithmic
  exception set;
* the reflection endpoint strip bound `O(R^2)`;
* the hard-core moving-anchor reduction (11)-(12).

**Not proved:** `(PV_eta)` / `(EC_eta)` itself.

**Exact obstruction:** the available Apéry identities are single-prime or
row-aggregate.  None supplies cancellation for

\[
 J_{p,q}(M)-A_pA_q/M
\]

when `p` and `q` vary and the residues are tied by the same moving integer
`m`.  Lucas collapses the incidences to prime divisors of `b_m`, but the
carrier height is only linear and therefore borderline; reflection is
same-prime; BFH omits the nontransverse channel and has the wrong row/pointwise
quantifier; the Hasse/ODE structure is a primewise evaluation statement whose
coefficient zeros are Mellin-dual and whose cross-characteristic product
still tensorizes by CRT.

So the missing theorem can be stated without ambiguity:

> **Prove a power-saving one-sided early-CRT covariance for the full Apéry
> coefficient-zero relation in short parabolic intervals `m<2LP`, with loss
> exponent strictly below `1/3`.**

That is enough for every `n`, includes `(11,5)`-type nontransverse targets, and
does not confuse a row-density theorem with a pointwise maximal theorem.

## Provenance note

The GitHub connector's visible `main` branch does not yet contain the Q76xx
files named in the current prompt.  I therefore used the prompt's current
`Xi`/transverse definitions as authoritative and did not invent a formula for
`Xi_r`.  The repo-backed inputs used above are the banked high-moment theorem,
`hm3_result.tex`, the exact Lucas/reflection machinery, and the existing Hasse
coefficient/evaluation audit.  To avoid overwriting a remote `main` that is
behind the stated working tree, these Q7697 files were placed on the dedicated
branch `chatgpt-q7697-pointwise-rainbow`.
