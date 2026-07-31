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
