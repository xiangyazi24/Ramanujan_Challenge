# Ramanujan Challenge — Status (2026-07-14)

## Score: 10/10 addressed (9 unconditional + 1 partial: density-1 unconditional, all-n conditional)

**Deadline:** August 1, 2026
**Total pages:** 51 (all compile clean)

## SOLVED (9 problems, unconditional proofs)

| Problem | Topic | Method | Pages | Status |
|---------|-------|--------|-------|--------|
| 2.1 | PCF → π | Q(√5) gauge, Poincaré roots | 4 | ✅ |
| 2.2 | γ Apéry | Aptekarev recurrence (index shift m=n+3) | 3 | ✅ |
| 2.3 | π+e | Lambert × derangement series | 3 | ✅ |
| 2.4 | harmonic+polylog | Weight-4 HPL symbolic summation | 3 | ✅ |
| 2.5 | Catalan CMF | Delannoy decomposition + k-recurrence | 10 | ✅ |
| 2.6 | ζ(2)+ζ(3) | GF ODE connection formula | 3 | ✅ |
| 2.7 | 4-term ζ(2)+ζ(3) | Adjoint certificate (Lagrange bracket) | 7 | ✅ |
| 2.8 | √10005/π | Chudnovsky formula in CMF disguise | 4 | ✅ |
| 3.1 | knot π² | A-polynomial / Mahler measure | 7 | ✅ |

## P3.2 — gcd(d_n a_n, d_n b_n) = e^{o(n)} for Apéry sequences

**Status:** ADDRESSED (7 pages). Two-tier proof.

### Unconditional results
1. **Z(p) = O(p^{2/3}) for ALL primes** — effective bound, no excluded primes
   - No consecutive zeros of b_j mod p (backward induction from b_0=1)
   - Gap-h polynomial C_h(m) has degree 3(h-1), leading coeff U_{h-1}(17) (Chebyshev)
   - **Nonvanishing lemma** (two-point proof): N_h(-1) = b_{h-1}·((h-1)!)³,
     N_h(-2) = -5·b_{h-2}·((h-2)!)³ → if N_h ≡ 0 mod p, consecutive zeros → ⊥
   - Effective: Z(p) ≤ p/H + 3(H-1)(H-2)/2 + 1, optimize H = (p/3)^{1/3}
   - Reflection law: N_h(-m-h-1) = (-1)^{h-1} N_h(m)
2. **Content restriction:** p | cont(N_h) ⟹ h ≥ 2p (for p ≥ 7)
3. **Resultant coprimality:** Res(N_2, N_3) = -5^6 (coprime mod p for ALL p ≥ 7)
4. **Density-1 unconditional:** G_n = e^{o(n)} for density-1 of n
   - First-moment argument over dyadic intervals using Z(p) = o(p)

### Conditional result (under Hypothesis Z)
5. Hypothesis Z (average Z(p) = O(1)): B(n) = O(1) → log G_n = O(√n) for ALL n

### Computational evidence (p ≤ 10^5, 9590 primes)
- Z(p) ∈ {0, 1, 2, 4, 6, 8, 10}, max = 10, mean ≈ 0.990
- P(Z=0) ≈ 61% ≈ e^{-1/2}; pair count Z/2 ∼ Poisson(1/2)
- Z(p) = 1 at exactly 2 non-ordinary primes (p | a_p(f), f = 8.4.a.a)
- Non-CM modular form → no systematic forced zeros (contrast: ζ(2) Apéry is CM)

### Key scripts
- `scripts/p32_zp_extended.py` — Z(p) for p ≤ 10^5
- `problems/3.2/proof.tex` — 7-page proof (compiles clean, 12 references)
