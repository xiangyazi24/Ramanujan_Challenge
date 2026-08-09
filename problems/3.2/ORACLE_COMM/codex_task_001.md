# Codex Task 001: Bilinear CRT Decorrelation — Precise Formulation

## Context
We are attacking the full unconditional Apéry GCD conjecture: G_n = e^{o(n)} for ALL n.

## What's landed
- Theorem A: b_r ≡ [x^r y^r z^r w^r] F^{p-1} (mod p), F = (1-x-y)(1-z-w) - xyzw
- Z(p) = O(p^{2/3}), gap polynomials N_h have degree 3(h-1), SL₂-tiling
- (HM)_2 proved unconditionally
- Separated-block resultant: Res(N_d, N_r(·+d)) ≠ 0 (verified for d+r ≤ 41)
- Mellin: b_r = Mellin transform of K3 fiber point-counts
- c(k) trade-off: need c + 2/k < 1 for full conjecture

## The precise CRT obstacle
(HM)_k for k ≥ 3 has CRT error (Σ Z(p))^k vs target X^2 λ_X^k — factor X^{k-2} gap.

## YOUR TASK

1. Read problems/3.2/proof.tex, specifically:
   - Section "AMTD" (§10.8, around line 3940)  
   - Section "The dispersion target" (around line 3253)
   - Section "Palm decorrelation" (around line 3634)
   - Remark "palindromic-fourier" (around line 4066)

2. Extract the EXACT bilinear sum that needs to be bounded for the dispersion inequality V° ≤ C·S.

3. Identify: what algebraic input about Z_p would make the bilinear sum cancellation work?
   Specifically: the gap polynomials N_h have bounded degree → their roots satisfy Weil-type bounds.
   Z_p is a PROJECTION of the (m,h)-incidence variety. How does this projection help?

4. State a CONCRETE lemma: "If [specific algebraic property of Z_p], then V° ≤ C·S."

5. Check: does the palindromic Fourier structure (the fact that e(a(p-1)/(2p)) S_p(a) ∈ ℝ) 
   give any cancellation in the bilinear expansion?

## Output format
Write your analysis to problems/3.2/ORACLE_COMM/codex_result_001.md.
Be specific — name exact sums, exact bounds, exact lemmas.
