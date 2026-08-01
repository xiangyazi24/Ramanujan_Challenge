# DS note: p=11 Dwork/Cartier computation — setup verified (2026-07-31)

Tool: **Sage** (correct for finite-field/p-adic Dwork work; Mathematica is unsuitable for mod-p).

## Verified setup (mod 11)
- Apéry b_m mod 11 correct (exact-integer computation; naive mod-p recurrence fails at m=p−1
  because (m+1)³ vanishes mod p — must fall back to exact values, as the session file noted).
- **PF operator annihilates F(t)=Σb_m t^m mod 11** (coefficient-wise m=2..22) — the rank-3
  Beukers-Peters connection matrix A(t) (a0,a1,a2 as rational functions mod 11) is correct.
- **Block law b_{p+r} ≡ 5·b_r mod 11** holds for r=1..10 (valid range).
- **Cartier operator acts as IDENTITY on the Apéry coefficient array**: b_{pm} ≡ b_m mod p
  (verified p=7,11,13, m=1..19) — from the vector Lucas b_{qp+r} ≡ b_q b_r (verified q=2,3 all r).
  This is the coefficient-level Frobenius/Cartier structure: the Apéry coefficient vector is
  a Frobenius eigenvector with eigenvalue 1.

## The decisive next step (Q6276 §5 — heavy, research-level)
Tr(Frob|H_Dwork) = Σ_{m<M} e_11(b_m) requires the FULL Griffiths-Dwork Cartier matrix on the
relative Dwork cohomology of F(t)=CT(1−tΛ)^{-1}. This is the Beukers-Vlasenko rational-form
Cartier operator (Prop 3.3) + Griffiths-Dwork reduction, which Q6174 honestly noted is NOT yet
executed. The recipe is in Q6153 (Cartier formula, precision N=2 then N=4) and Q6276 §5.
- The overconvergent F-isocrystal structure + additive-character twist on the coefficient
  coordinate is the object; its Frobenius slopes/dimensions bounded as p varies would confirm
  the p-adic framework.
- This is a substantial implementation (multi-hour), not a quick grind.

## Katz-sheaf moment bridge verified (p=11, DS)
- CT Λ^m = b_m mod 11 (m=0..8) — the constant-term realization.
- **sum_x Λ(x)^m ≡ (p−1)³·b_m ≡ −b_m (mod 11)** for m=1..7 — Claude's Katz-sheaf bridge
  (c_m = moment = (p−1)^d·b_m, d=3) CONFIRMED. (My first comparison failed only by missing
  the mod-11 reduction.)
- Cartier operator on the coefficient array: C(b)_m = b_{pm} ≡ b_m (m=1..4) — the coefficient
  Dwork/Frobenius eigenvector structure.
- Exponential sums |Σ_{m<M} e_11(b_m)| = 2.73, 6.06, 7.79, 13.99 (M=11,22,33,55) ~ √M — the
  random-walk scale (what the Dwork framework must explain, not sheaf-cancellation).

## Griffiths-Dwork Cartier matrix (p=11) — DS computed
- **CT(C(e₀)) = b_m EXACTLY** (m=0..25, no mismatches) — the Cartier operator on the e₀ period
  reproduces the Apéry numbers.
- **Cartier matrix on the period vectors**: M = [[1,0,0],[0,1,0],[0,1,0]] (mod 11), verified
  over m=1..30 with ZERO mismatches. Meaning: C(e₀ period)=e₀, C(e₁ period)=e₁,
  C(e₂ period)=e₁ (the e₂ correction's Cartier CT cancels into e₁).
- Interpretation: the Cartier/Frobenius action on the Apéry periods is the IDENTITY (+ a
  collapse) — Frobenius-simple, slopes 0, bounded. This is the "bounded slopes/dimensions"
  structure Q6276 asked about, at the PERIOD level. (The full cohomology Cartier on the forms
  mod exact-forms is beyond this CT-level; but the period level is what controls the Apéry
  sequence / the exponential sums' target values.)
- The decisive identity Tr(Frob|H_Dwork)=Σe_p(b_m) still needs the p-adic overconvergent
  machinery (splitting function + Banach space + finite-dim extraction) — the next heavy step.
