# Ramanujan Challenge — Status (2026-07-13)

## Strategy (Fable oracle, confirmed)

**Breadth-first, CAS-first.** Submit each problem as solved. Target: 5-7 problems.
Section 2 = reverse-lookup (find which known result is disguised). 2.2 confirmed this.
Build ore_algebra pipeline on uisai2 (days 1-2), amortize across 6 problems.
Lean only as garnish after CAS proofs.

## Progress

### SOLVED
- **2.2 (γ Apéry)** ✅ — Aptekarev recurrence (index shift m=n+3). Proof + PDF written.

### In Progress  
- **2.3 (π+e)**: Poincaré roots 0,0,1±√2. Fable predicts LCLM of two order-2 ops (π+e).
  Need ore_algebra factorization on uisai2. HIGHEST LEVERAGE TEST.
- **2.1 (PCF→π)**: Gauge (n!)³, roots -110±50√5 ∈ Q(√5). ChatGPT: no named match yet.
- **2.4 (harmonic+polylog)**: Two-stage symbolic summation → weight-4 HPLs at 1/2. CAS job.
- **2.6 (ζ(2)+ζ(3))**: Poincaré roots 1, 1/4. GF ODE known. Need connection formula.
- **2.5 (Catalan CMF)**: Order-3 = summation lift of order-2. Rivoal-Zudilin connection.
- **2.7 (4-term ζ(2)+ζ(3))**: Related to 2.6. Needs numerical verification first.

### Bounded Budget / Dropped
- **2.8 (√10005/π)**: Possibly Chudnovsky-related (discriminant). Bounded exploration only.
- **3.1 (knot π²)**: One-shot literature search only.
- **3.2 (Apéry measure)**: DROPPED — genuinely open, zero EV in 19 days.

## Attack Order (Fable-recommended)
1. ~~2.2~~ ✅ DONE
2. 2.3 factorization test (day 1-2, cheap, huge information value)
3. 2.4 (days 3-7, most mechanical)
4. 2.6 (days 4-8, numerics-then-prove)
5. 2.1 (days 6-10, with pipeline)
6. 2.7 (days 9-13, template from 2.6)
7. 2.5 (days 10-15)
8. 2.3 full solve + 2.8 attempt (days 12-17)

## Infrastructure TODO
- [ ] ore_algebra + Mathematica pipeline on uisai2
- [ ] 500-digit numerical verification for all problems
- [ ] OEIS lookups for all sequences
- [ ] Poincaré root fingerprinting against known families
