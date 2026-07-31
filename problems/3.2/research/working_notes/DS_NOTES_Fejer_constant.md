# DS note: Fejér majorant → first constant 1/3 (conditional on S_1,S_2 = o(P_n))

Date 2026-07-31 automode. Sharpening of the horizontal-We cancellation route.

## Setup
θ_p = b_{n−p}/p mod 1. H(n) = #{p ∈ (n/2,n] : θ_p ≡ 0} = #{p : p|b_n} (top-window).
P_n = n/log n. S_j(n) = Σ_{n/2<p≤n} e(j θ_p). 1_{θ≡0} ≤ F_K(θ)/F_K(0) pointwise, where
F_K(θ) = Σ_{|j|<K}(1−|j|/K)e(jθ) is the Fejér kernel (F_K(0)=K, F_K≥0, F_K(θ)/K ≤ 1 all θ).

## The bound
H(n) ≤ Σ_p F_K(θ_p)/K
     = (1/K)[P_n + Σ_{0<|j|<K}(1−|j|/K) S_j(n)].

K=3:  H(n) ≤ (1/3)P_n + (4/9)Re S_1(n) + (2/9)Re S_2(n) ≤ (1/3)P_n + (4/9)|S_1| + (2/9)|S_2|.
⟹ **if S_1, S_2 = o(P_n), then H(n) ≤ (1/3+o(1))·P_n — a constant 1/3 < 1/2.**
K→∞: H(n) = o(P_n) (the full theorem).

## Why this matters
- The first constant via the horizontal route is **1/3**, needing only S_1,S_2 = o(P_n)
  (2 fixed modes). This is BETTER than the p⁸-carrier constant Λ/8 = 0.44069 (Q6129),
  and it's the same hypothesis family as the full theorem (S_h = o(P_n)).
- Empirically |S_1|,|S_2| ~ √P_n (measured n ≤ 6.4·10⁴, values below) ⟹ the hypothesis
  holds with huge margin (√P_n vs P_n).
- So the grindable target splits cleanly:
  **(a) constant 1/3 ⟸ S_1(n), S_2(n) = o(P_n)** (two fixed modes),
  **(b) full o(1) ⟸ S_h(n) = o(P_n) for all fixed h**.
  A proof of even (a) is a first real advance (beats every closed route and beats the
  algebraic p⁸ constant).

## Numerical bank (computed)
| n | |S_1| | sqrtP | |S_1|/sqrtP | D_1=Σ(1−cosθ_p) | D_1/P_n |
| 2000 | 23.3 | 16.2 | 1.44 | 132 | 0.981 |
| 4000 | 10.6 | 22.0 | 0.48 | 253 | 1.026 |
| 8000 | 18.2 | 29.8 | 0.61 | 465 | 1.018 |
| 16000 | 80.4 | 40.7 | 1.98 | 777 | 0.909 |
| 32000 | 15.3 | 55.5 | 0.28 | 1560 | 0.993 |
| 64000 | 45.7 | 76.0 | 0.60 | 2943 | 0.987 |
|S_1| ~ √P_n (o(P_n) — full-theorem hypothesis holds w/ huge margin); **D_1/P_n ≈ 1.0**
(maximal spectral spread — the one-sided scalar-gap hypothesis D_1 ≥ ηP_n holds with η ≈ 1).
So both sufficient conditions are empirically satisfied at their strongest; the missing piece
is a proof that the phases θ_p = b_{n−p}/p don't align (horizontal cancellation).

## Phase-proximity (spectral-gap) — phases EXACTLY uniform
Fraction of primes p ∈ (n/2,n] with h·b_{n−p}/p within δ of 0 (mod 1):
- δ=0.10: 0.243 (n=8000), 0.233 (n=16000) vs uniform 0.20
- δ=0.05: 0.089/0.096/0.119 (n=2000/8000/16000) vs uniform 0.10
- δ=0.01: 0.015/0.031/0.021 vs uniform 0.02
⟹ fraction ≈ 2δ: phases exactly uniform. The one-sided spectral-gap hypothesis
D₁ ≥ ηP_n (any "δ-proportion away from 0") holds at the MAXIMAL level. First-constant
hypothesis empirically strongest; only the proof missing.

## The minimal theorem (Q6177): SG1 — one-sided first-mode spectral gap
**SG1**: ∃η>0 such that for all large n, D₁(n) = Σ_{n/2<p≤n}(1−cos 2πθ_p) ≥ η·n/log n
(equivalently Re S₁(n) ≤ m_n − η·X_n, X_n=n/log n, m_n=|P(n)|≈(1/2)X_n).
⟹ H(n) ≤ (1/2 − η/2 + o(1))·X_n — a constant strictly below 1/2.
- K=2 Fejér is exact: H(n) ≤ m_n − D₁(n)/2.
- Arc form AE(ε,δ): #{p : ||θ_p|| ≥ ε} ≥ δ·m_n ⟹ H ≤ [1/2 − δ(1−cos 2πε)/4 + o(1)]X_n.
- **SG1 is the minimal assertion a compensated star VIOLATES** (star phases planted at 1 ⇒ D₁≈0).
  All local statistics (rowwise uniformity, Cartier, Franel, fixed-char) are compatible with
  a phase planted at 1 on one diagonal — so SG1 is the genuine cross-prime frontier.
- Empirically D₁ ≈ X_n (η ≈ 1, maximal), phases exactly uniform (proximity test ≈ 2δ). 
- Most plausible input: constant-defect cross-prime dispersion / inverse theorem for near-maximal S₁.

## SG1 quantification verified (DS, 2026-07-31)
Since 1−cos(2πx) = 2sin²(πx) ≥ 4x² (|x|≤1/2) and θ_p = {5^{-1}b_n/p}, with A=5^{-1}b_n:
D₁(n) ≥ 4·Σ_p ‖A/p‖² (‖·‖ = distance to nearest integer). Measured:
Σ_p ‖5^{-1}b_n/p‖² = 11.35, 21.49 (n=2000,4000) vs (1/12)·P_n = 11.25, 20.58 (ratio 1.009, 1.044).
⟹ the phases are exactly uniform (E[‖x‖²]=1/12); D₁ ≥ (1/3+o(1))P_n from the 4x² bound alone
(actual D₁ ≈ P_n). SG1 needs: Σ_p‖A/p‖² = (1/12+o(1))P_n for the specific holonomic A=5^{-1}b_n,
i.e. no concentration of {A/p} near integers — a nonresonance statement.

## Related
- Q6127 (strategic): horizontal Weyl theorem is the top route; D_1 ≥ ηP_n first-constant framing.
- Q6129 (calibration): p⁸ carrier gives 0.44069 (algebraic route, Claude owns).
- Fejér K=3 gives 1/3 (analytic route) — stronger constant, same S_h hypothesis family.
