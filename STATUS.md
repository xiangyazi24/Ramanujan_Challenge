# Ramanujan Challenge — Status (2026-07-14)

## Score: 10/10 addressed (9 unconditional + 1 conditional on explicit hypothesis)

**Deadline:** August 1, 2026

## SOLVED (9 problems, unconditional proofs)

| Problem | Topic | Method | Status |
|---------|-------|--------|--------|
| 2.1 | PCF → π | Q(√5) gauge, Poincaré roots | ✅ proof.tex |
| 2.2 | γ Apéry | Aptekarev recurrence (index shift m=n+3) | ✅ proof.tex |
| 2.3 | π+e | Series identification | ✅ proof.tex |
| 2.4 | harmonic+polylog | Weight-4 HPLs symbolic summation | ✅ proof.tex |
| 2.5 | Catalan CMF | Rivoal-Zudilin connection | ✅ proof.tex |
| 2.6 | ζ(2)+ζ(3) | GF ODE connection formula | ✅ proof.tex |
| 2.7 | 4-term ζ(2)+ζ(3) | Adjoint certificate (Lagrange bracket) | ✅ proof.tex |
| 2.8 | √10005/π | Chudnovsky formula in CMF disguise | ✅ proof.tex |
| 3.1 | knot π² | A-polynomial / Mahler measure | ✅ proof.tex |

## P3.2 — gcd(d_n a_n, d_n b_n) = e^{o(n)} for Apéry sequences

**Status:** ADDRESSED (Section 3 = open problem). Layered proof: unconditional reduction + conditional conclusion.

**What the proof establishes:**
1. **Unconditional:** log G_n = O(√n) + O(log n) · B(n) where B(n) counts bad primes
2. **Unconditional:** Small primes (p ≤ √n) contribute O(√n) via Wronskian
3. **Unconditional:** Denominator connection lemma: for p ∈ (n/2, n], v_p(G_n)=0 iff p∤b_{n-p}
4. **Conditional on Hypothesis Z (Z(p)=o(p) avg):** B(n) = O(1), hence log G_n = O(√n) = o(n)
5. **Equivalence:** conjecture ⟺ average Hypothesis Z

**Hypothesis Z evidence (p ≤ 10^4, 1227 primes):**
- Z(p) ∈ {0,1,2,4,6,8}, max=8, mean=0.957
- P(Z=0) = 61.6% ≈ e^{-1/2} = 60.7% (Poisson model)
- Power law: Z(p) ~ 2.0 p^{0.02} ≈ O(1)
- Symmetry b_j ≡ b_{p-1-j} verified 100%

**Why Hypothesis Z is open:** Contains the non-ordinary prime density problem for
weight-4 non-CM modular forms (Jin-Ma-Ono 2016) as a special case.

**Key scripts:**
- `scripts/p32_gcd_analysis.py` — GCD evidence to n=200
- `scripts/p32_extended_analysis.py` — GCD evidence to n=500
- `scripts/p32_lucas_verify.py` — Lucas congruence verification
- `scripts/p32_zero_count.py` — Z(p) for p ≤ 997
- `scripts/p32_zp_extended.py` — Z(p) for p ≤ 10^4
- `scripts/p32_denom_verify.py` — Denominator connection verification
- `problems/3.2/proof.tex` — 5-page proof (compiles clean)
