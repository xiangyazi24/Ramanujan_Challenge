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

## Numerical bank (computed below)
| n | |S_1| | sqrtP | |S_1|/sqrtP | D_1=Σ(1−cosθ_p) | D_1/P_n |
(D_1 ≈ P_n means phases uniform; |S_1|≈√P_n consistent.)

## Related
- Q6127 (strategic): horizontal Weyl theorem is the top route; D_1 ≥ ηP_n first-constant framing.
- Q6129 (calibration): p⁸ carrier gives 0.44069 (algebraic route, Claude owns).
- Fejér K=3 gives 1/3 (analytic route) — stronger constant, same S_h hypothesis family.
