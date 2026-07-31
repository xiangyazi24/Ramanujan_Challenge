# DS audit: gaps in Claude4.6's SUBMIT/3.2 (cross-check 2026-07-31, after push)

Purpose: items from the DS synthesis (Q6123–Q6221, 33 ChatGPT consultations + DS
verification) that would sharpen Claude's submission. All are VERIFIED. Claude's
submission is honest and correct as far as it goes; these are additions, not corrections
to errors (one terminology fix excepted).

## 1. THE missing headline: the 4th-moment moment criterion (Q6206, exact conversion)
Claude's solution.tex uses fourth-moment only for the covariance/collision problems, NOT
for the horizontal Weyl sums. The cleanest sufficient condition for the whole conjecture:
```
THEOREM (moment criterion). If for each fixed h, Σ_{N<n≤2N}|S_h(n)|⁴ ≤ C·N³/log²N
uniformly in dyadic N, then (max ≤ (Σ|S|⁴)^{1/4}, trivial L⁴⊂L^∞ on the shell):
  max_{N<n≤2N}|S_h(n)| ≤ C^{1/4}·N^{3/4}/log^{1/2}N = o(n/log n)  [uniform, ALL n].
Fejér (H(n) ≤ P_n/K + (1/K)Σ_{0<|j|<K}(1−|j|/K)|S_j|) then closes the conjecture.
```
- Markov gives E_ε(N) = #{n : |S_h|≥εQ_n} ≤ (16C+o(1))log²N/(ε⁴N) → 0, i.e. eventually
  NO exceptional index for fixed ε (one bad index would contribute ≫N⁴/L⁴ to the budget N³/L²).
- Any fixed random-scale moment M_{2k}(N) ≪ N^{k+1}/log^k N (k≥2) suffices; k=2 (4th) is
  the first rung; k=1 (2nd) fails by a √log factor.
- This is the sharpest statement of the frontier: the conjecture ⟸ the 4th moment of S_h
  at the random scale. (Not in the submission; the submission states the obstruction but
  not this decisive gate.)

## 2. SG1 — the minimal first-constant theorem (Q6177)
```
SG1: ∃η>0, for all large n: D₁(n) = Σ_{n/2<p≤n}(1−cos 2πθ_p) ≥ η·n/log n
⟹ H(n) ≤ (1/2−η/2+o(1))·n/log n   (K=2 Fejér exact: H(n) ≤ m_n − D₁/2).
```
- Empirically D₁ ≈ n/log n (η≈1, phases exactly uniform). SG1 is the minimal assertion a
  compensated star violates; all local statistics are compatible with it.
- 2D₁(n) = Σ_p |e(θ_p)−1|² (Q6207). Inverse theorems only conclude phases cluster near 1
  (the star), no contradiction.
- (Not in submission; would give the submission a clean "first constant" target.)

## 3. The p⁸-carrier threshold (Q6129 calibration) — Claude's own grade route's target
```
Λ = log(17+12√2) = 3.52549.  A target-selective p^k carrier of height Λn gives
H(n) ≤ (Λ/k+o(1))·n/log n.  Λ/7 = 0.5036 > 1/2 (useless); Λ/8 = 0.44069 < 1/2 (FIRST useful).
```
- The grade route needs p⁸ at Λ-height (not p⁷, not deeper-but-superlinear scaled-index
  congruences, which don't furnish the carrier). This is the quantitative target for
  Claude's own active grade-g supercongruence direction. (Not in submission.)

## 4. Why counting methods fail — two citable facts (Q6125/Q6181)
- **Perret-Gentil blindness**: auxiliary-ℓ-adic Frobenius sieves reduce the trace at
  ℓ ≠ p and are LITERALLY BLIND to p-divisibility — the defining characteristic is the only
  valuation that sees it. (Not in submission.)
- **Saffari–Vaughan regime**: the closest theorem for the reciprocal-prime sum
  S_h(n)=Σ_p e(h·5^{-1}b_n/p) needs prime cutoff Y > X^{6/11+ε} (polynomial range); for
  Apéry X=b_n, Y=n=X^{o(1)} — misses exponentially. (Not in submission.)

## 5. Terminology fix (Q6125) — the one correction
- Eigenvalues of size p^{3/2} have **Weil weight 3**, not "weight 4" ("weight 4" is the
  classical modular-form weight). The session files mix these; the submission should use
  Weil weight 3 for the Frobenius size and modular weight 4 for the newform 8.4.a.a.

## 6. Strengthened empirical bank for the "why true" section (all DS-verified)
- max H(n) = 3 for n ≤ 8·10⁵ (extends the submission's 5·10⁵).
- R₂(N)=Σ|S_h|²/ΣQ_n → 1 (N=200–1200); R₄≈0.87–0.98 (h=1), 1.15–1.19 (h=2), 0.83–0.95
  (h=3); R₆≈0.73–0.93; Λ_h (Gram spectral norm) ≈ 1.01–1.06.
- **Joint S₁–S₂ independence**: ρ(S₁,S₂) = 0.091/0.031/0.022 ≈ 1/√N (independence pred
  0.071/0.058/0.045); joint 4th moment J/Jpred ≈ 1.0–1.3.
- **Gumbel**: max_n|S₁|²/Q_n = 4.72/4.98/4.98 vs log N = 5.30/5.70/6.21 (extreme-value scale).
- Off-diagonal 4th-moment excess slightly NEGATIVE (−1.6% to −13%).
⟹ the phase family {b_{n−p}/p} is empirically indistinguishable from i.i.d. uniform; the
  random-phase model is the exact content of "why the conjecture is true."

## 7. Other recorded closures (in DS notes; submission can cite)
- p-curvature route refuted (det=1, not nilpotent; C(p−1) not a conjugacy invariant).
- CD (Christoffel–Darboux) coboundary absent — the mod-24 square factorization is
  analytically exhausted for S_h (sympy, 174–524 overdetermined eqs).
- Pair-Gram spectral norm (unrestricted) FORCED to diverge by dimension (≍N²/log²N pair
  columns in N-dim space) — wrong target; corrected = weighted Λ(4)/2→4 inequality.
- Randomization typicality (Q6173): the conjecture = quenched-vs-annealed comparison;
  the random no-star is an elementary Chernoff/occupancy theorem.

## Files with the full details
DS_NOTES_SYNTHESIS.md, DS_NOTES_moment_route.md (4th-moment + joint + Gumbel),
DS_NOTES_Fejer_constant.md (SG1), DS_NOTES_reciprocal_prime.md (Saffari-Vaughan),
DS_NOTES_MellinTraces.md (Perret-Gentil, Weil weight), DS_CLAUDE_COLLAB.md (all Q#s).
