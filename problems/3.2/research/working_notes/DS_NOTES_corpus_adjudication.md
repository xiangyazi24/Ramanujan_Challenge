# DS adjudication: the E(p) → pointwise corpus contradiction (2026-07-31, Fable-flagged)

## The contradiction
- Q3.2_density_theorem.md (~line 1069): "Since E(p)=O(p) (indeed ≪ p^{1+o(1)}) suffices for
  mass X^{3/2}=URE=pointwise Problem 3.2".
- HANDOFF singleton-aligned countermodel: S_p = {N−p, p−1−(N−p)} for p∈(N/2,N] has all
  vertical properties (|Z_p|=2, reflection, no consecutive) and small vertical energy, yet
  H(N)=P_N (level N aligned).

## VERDICT: the density-doc claim is WRONG; the countermodel is RIGHT.
E(p) = Σ_c #{m∈[0,p−1]: b_m≡c}² is a VERTICAL quantity (per-prime collision energy). It
enters the gcd chain ONLY through |Z_p| ≤ √E(p). The countermodel has Σ_p|Z_p|² = O(N/log N)
(small vertical energy) yet H(N)=P_N. **Vertical data (small |Z_p|, small E(p)) cannot give
the pointwise bound** — the singleton alignment is invisible to any per-prime/vertical
quantity. VERIFIED: |Z_p|=2 all p, Σ|Z_p|²=540=O(N/log N) at N=2000, H(N)=135=P_N.

## Correct statement
E(p) = O(p) suffices for the **DENSITY-ONE / exceptional-set** version (the mass X^{3/2}
argument via the first moment / Markov — the URE chain), NOT the pointwise. The pointwise
needs a **NON-VERTICAL (horizontal / cross-prime) input** — exactly the missing cross-prime
decorrelation theorem (Apéry large sieve / SG1 small-ball / local-limit law). Any claim that
vertical E(p) alone closes pointwise Problem 3.2 is invalid.

## Corpus fix
- Q3.2_density_theorem.md ~line 1069: annotate "suffices for mass X^{3/2} = density-one /
  exceptional-set, NOT pointwise; the pointwise requires a horizontal input (the countermodel
  is a valid obstruction)."
- STATUS.md / any doc repeating "E(p)=O(p) ⟹ pointwise": same annotation.

## Implication for the campaign
E(p) = 3p + O(√p) is empirically true (verified p=101..307) and is a genuine Weil-signature
fact, but it is a VERTICAL statement. It strengthens the exceptional-set results; it does NOT
close the pointwise conjecture. The pointwise frontier remains the horizontal/cross-prime
decorrelation (the same wall every route hits).
