# DS AVENUES — autonomous run 2026-07-31 (parallel track to Claude4.6's p⁸ route)

**Goal (one sentence):** prove `H(n) = #{n/2<p≤n : p|b_n} = o(n/log n)` (⟺ P3.2 log G_n=o(n)),
via the horizontal-We-cancellation / F₂-dispersion route — DS drives ChatGPT ds1–ds5 + its own
numerical/Ore-algebra grind; Claude4.6 (dm window) owns the p⁸-carrier algebraic route in parallel.

## State (2026-07-31, from DS synthesis — see DS_NOTES_SYNTHESIS.md)
Five independent ChatGPT consultations (Q6123–Q6130) + DS numerical verification:
- Every route from fibre hypotheses is PROVABLY dead (compensated-star countermodel, CRT √S
  threshold, Perret-Gentil blindness, p-curvature det=1). No existing theorem does pointwise.
- **Target theorem (verified true in data): horizontal Weyl cancellation**
  `S_h(n) = Σ_{n/2<p≤n} e(2πi h·b_{n−p}/p) = o(P_n)`, P_n=n/log n. Fejér ⇒ H(n)=o(P_n).
  DS measured |S_1|,|S_2| ~ O(√P_n) for n=10³..1.6·10⁴ (pushing to 6.4·10⁴).
- **Fejér sharpening (DS, banked in DS_NOTES_Fejer_constant.md)**: H(n) ≤ (1/K)P_n +
  (1/K)Σ_{0<|j|<K}(1−|j|/K)S_j(n). **K=3 ⇒ if S_1,S_2=o(P_n) then H(n) ≤ (1/3+o(1))P_n**
  — first constant 1/3, BETTER than the p⁸ route's Λ/8=0.44069. Full o(1) via K→∞.
- **Verified structure:** truncated Apéry poly A_p(T)=Σ_{k<p}b_kT^k is square iff p≡1,5,7,11
  (mod 24), fixed quadratic×square otherwise (CFVZ arXiv:2510.23298; DS verified 7≤p≤400).
  A_p = Hasse–Witt scalar CT(1−tΛ)^{p−1}. Value: rank-lowering to Franel rank-2 orbit.
- Claude's route (parallel): target-selective p⁸ carrier at Λ-height ⇒ first constant Λ/8=0.44069.

## DS avenues (ranked; each complete attack plan)
- **(a) Christoffel–Darboux coboundary gate test** (Q6130's falsifiable test): build the two
  rank-two Franel-root sequences (branches (2+)/(2−)) whose reductions give B_p; set up the
  convolution Σ_i g_ε(i)g_ε(r−i); decide whether a rational ρ_ε(r) + matrix M_ε give the
  order-zero coboundary collapse (creative telescoping / Ore algebra, Sage ore_algebra or
  Mathematica). Terminal: coboundary exists (new analytic variable for S_h) OR no rational
  solution on both branches (square structure analytically exhausted for S_h — record + move).
- **(b) F₂ cross-characteristic dispersion** (Q6129): measure Σ_{p<q} C_{p,q}(N) with
  C_{p,q}(N)=#{s : p|b_{s+q−p}, q|b_s, N<q+s≤2N}; confirm = o(N²/log²N) with the predicted
  N/log²N scale; look for structure in the pair correlation. Terminal: prediction confirmed
  (banked + next attack) OR anomaly (new direction).
- **(c) S_h class-split by mod 24**: measure S_{h,a}(n)=Σ_{p≡a (24)} e(h·b_{n−p}/p) per class a,
  check the square classes behave differently (spectral gap). Terminal: class-dependent signal
  (target for (a)) OR uniform (no gain from mod-24).
- **(d) Standing ChatGPT saturation (duty, not avenue): keep ds1–ds5 full** with the backlog:
  coboundary construction details, F₂ numerical probes, DGK Test B Sage recipe, two-polynomial
  inverse theorem candidates, horizontal-Weyl proof strategies. Dispatch-before-process.

## Fallbacks
(a) no rational solution ⇒ the horizontal route needs a genuinely new cross-prime input; hand the
exact negative certificate to Claude + bank it. (b)/(c) clean ⇒ the analytic evidence is complete;
theorem awaits the new equidistribution input.

## Terminal conditions (any avenue)
Success = a banked, verified theorem/certificate advancing H(n)=o(P_n) or the first constant.
Proof-of-failure = a written, reproducible certificate (e.g. Ore algebra says "no rational
solution", or a counterexample), not a feeling.
