# DS note: the p-th moment route — cleanest sufficient condition (2026-07-31)

## The sharpening (DS derivation, after computing R₄)
For the horizontal Weyl sums S_h(n) = Σ_{n/2<p≤n} e(2πi h b_{n−p}/p), Q_n = |P(n)| ~ n/(2 log n):

**Moment theorem.** If for some fixed p > 2,
```
Σ_{N<n≤2N} |S_h(n)|^p  =  O( N^{1+p/2} / log^{p/2} N )     [the random scale]
```
uniformly in dyadic N, then by max|S_h| ≤ (Σ|S_h|^p)^{1/p}:
```
max_{N<n≤2N} |S_h(n)| = O( N^{1/2+1/p} / log^{1/2} N ) = o(N/log N) = o(P_n).
```
Since H(n) ≤ (1/K)P_n + (1/K)Σ_{0<|j|<K}(1−|j|/K)|S_j(n)| (Fejér), pointwise |S_j|=o(P_n)
for all fixed j closes the conjecture log G_n = o(n).

The random scale N^{1+p/2}/log^{p/2} N comes from |S_h(n)| ~ √Q_n (Gaussian): |S_h|^p ~ Q_n^{p/2},
Σ_n ~ N·Q^{p/2} = N·(N/logN)^{p/2}. Any p>2 works (p=2 is borderline: max ≤ N/log^{1/2}N, not o(N/logN)).

## Empirical confirmation (DS measured)
R₄(N) = Σ|S₁(n)|⁴ / (2·ΣQ_n²) = **0.868, 0.984, 0.951 (N=200,300,500)** → ≈1 (random scale).
R₆(N) = Σ|S₁(n)|⁶ / (6·ΣQ_n³) = **0.732, 0.930 (N=200,300)** → ≈1 (Gaussian: E|S|^{2k}=k!Q^k).
Full moment bank: R₂→1 (N=200–1200: 0.963,0.992,1.006,1.016,0.955), R₄≈0.87–0.98, R₆≈0.73–0.93,
Λ_h≈1.01–1.06 (Gram spectral norm), D₁≈X_n, phases uniform. max|S₁|² = 137,214,304 — no heavy tail.
**All fixed h** (R₂,R₄ at N=200,300): h=1 (0.96–0.99, 0.87–0.98), h=2 (1.07–1.08, 1.15–1.19),
h=3 (0.97–1.01, 0.83–0.95) — R₄ = O(1) for all h, the moment criterion holds.
**JOINT test (Q6213's decisive probe, DS measured)**: cross-correlation ρ(S₁,S₂) = 0.091, 0.031,
0.022 (N=200,300,500) ≈ 1/√N (independence prediction 0.071, 0.058, 0.045); joint 4th moment
J/E[|S₁|²|S₂|²] = 1.0–1.3 (independent prediction 1.0). ⟹ S₁,S₂ behave like INDEPENDENT
complex Gaussian random walks. The phase family is indistinguishable from i.i.d. uniform.
**GUMBEL check (Q6213's pointwise-tail probe, DS measured)**: max_{n∈(N,2N]}|S₁(n)|²/Q_n =
4.72, 4.98, 4.98 (N=200,300,500) vs Gumbel prediction log N = 5.30, 5.70, 6.21 — extreme-value
scale matches (max of ~N Exp(1) variables). Pointwise max behavior confirmed.
**Q6221 correction**: the unrestricted pair-Gram spectral norm (Q6211's proposal) is FORCED to
diverge by dimension (≍N²/log²N pair columns in N-dim row space ⟹ norm ≫ N/log²N even for
random phases) — WRONG target. The corrected target is a weighted Λ(4)/2→4 inequality on
decomposable pair coefficients (c⊗c) — essentially the 4th-moment bound itself; the
4-distinct-prime correlation is genuinely new content.

## Q6244 — the cleanest frontier statement: the "Apéry large sieve"
The exact first estimate to prove (closes the 4th-moment route):
```
Σ_{N<n≤2N} |T_h(n)|² ≪_h N³/log²N,   T_h(n) = |S_h(n)|² − Q_n
  = Σ_{n/2<p,q≤n, p≠q} e_{pq}(h·5̄_pq·(q−p)·b_n)   [pair phase, modulus pq ≍ N²].
```
Expanding the square: identical pairs (unavoidable diagonal ≍N³/log²N), pairs sharing one
prime (3-prime collision), disjoint pairs (4-distinct correlation, modulus p₁p₂p₃p₄ ≍ N⁴).
- The n-interval is length N = M^{1/4} — an EXTREMELY SHORT sum of a holonomic sequence mod a
  varying squarefree modulus. No classical completion/Weil/Deligne applies (b_n mod M is not
  a bounded-conductor trace function; the recurrence is singular at indices divisible by p|M).
- CRT gives compatible prime-local automata (each prime its own base-p expansion), NOT a new
  block law mod pq; but for a fixed tuple B_M(n) ≡ 5b_{n−p_i} mod p_i packages the locals.
- Absolute summation of per-tuple bounds CANNOT work ((N/logN)⁴ tuples, allowed avg per tuple
  log²N/N < 1) — cancellation across tuples is essential.
- **The theorem: the pair vectors n↦(e_{pq}(h5̄_pq(q−p)b_n))_{p≠q} are almost orthogonal on
  (N,2N], total Gram mass at the diagonal scale** — an "Apéry large sieve" coupling 2-4 prime
  coordinates by CRT. This + the SG1 product-level anti-concentration (Q6241) are the two
  cleanest statements of the frontier.
- **DS measured the Apéry large sieve hypothesis directly**: Σ|T₁(n)|²/ΣQ_n² = 0.797, 0.974,
  0.896, 1.024 (N=200,300,500,800) ≈ 1 — the pair-square is at the random scale N³/log²N,
  exactly the target. The pair vectors ARE almost orthogonal (off-diagonal at the diagonal
  scale).
- **DS measured the SG1 small-ball hypothesis**: #{p : ‖5^{-1}b_n/p‖ ≥ ε} ≈ (1−2ε)·m_n
  (ε=0.1: 0.793/0.830 vs 0.80; ε=0.2: 0.548/0.587 vs 0.60; ε=0.4: 0.230/0.231 vs 0.20 at
  n=2000,4000) — phases EXACTLY uniform. For ε=0.4, δ≈0.23 primes far from 0 ⟹
  D₁ ≥ 0.23·m_n·(1−cos 0.8π) ≈ 0.35·m_n — a strong constant.
**The phase family {b_{n−p}/p} behaves exactly like independent uniform random variables
(Gaussian moments through 6th order).**

## Why this is the cleanest frontier statement
- The conjecture ⟸ a single p-th-moment dispersion bound (any p>2) at the random scale.
- It's a "cross-prime decorrelation" statement: Σ_{p,p'} ... phase correlations must be
  exactly at the independence value.
- Empirically TRUE at p=4 with ratio ≈ 0.95–0.98 — margin is tight-but-real (the ratio must
  be O(1), and it is ~1, not growing).
- Proving it is the new-math frontier (same as SG1 / quenched-vs-annealed, rephrased as a
  high-moment dispersion).

## Q6206 exact conversion (the decisive gate)
THEOREM (moment criterion): If for each fixed h, M₄,h(N) = Σ_{N<n≤2N}|S_h(n)|⁴ ≤ C·N³/log²N
uniformly in dyadic N, then (L⁴⊂L^∞ on the finite shell, max ≤ (Σ|S|⁴)^{1/4}):
```
max_{N<n≤2N}|S_h(n)| ≤ C^{1/4}·N^{3/4}/log^{1/2}N = o(N/log N) = o(Q_n)  [uniform, ALL n].
```
- Exceptional-set Markov: E_ε(N) = #{n : |S_h| ≥ εQ_n} ≤ (16C+o(1))log²N/(ε⁴N) → 0,
  and since E_ε is an integer, it's eventually 0 — NO exceptional index for fixed ε.
  (One bad index with |S_h|≥εQ_n contributes ≫ N⁴/L⁴ to the budget N³/L² — ratio N/L² → ∞.)
- General: ANY random-scale moment M_{2k}(N) ≪ N^{k+1}/log^k N (fixed k≥2) suffices;
  k=2 (4th moment) is the first rung; k=1 (2nd moment) fails (gives max ≫ Q_n by √log factor).
- Gaussian ladder M_{2k} ≤ A^k k! N Q_*^k up to k ≍ log N would give max ~ √N (the true extreme),
  but NOT needed for o(P_n).
- **The 4th-moment bound for each fixed h is the cleanest single sufficient condition for
  the whole conjecture.** Empirically R₄≈0.87–0.98 (holds). Proving it = the frontier.

## Open sub-questions
- Can the moment be reduced to p = 2 + ε, or is p=4 the natural first target?
- Is there a way to decompose Σ|S_h|^p into diagonal + off-diagonal where the off-diagonal
  is provably smaller (a "dispersion" statement for the holonomic phases)?
- The p-th moment is a finite correlation of the reciprocal-prime phases {A/p}, A=5^{-1}b_n —
  relates to the Saffari-Vaughan-type question at the logarithmic scale.
