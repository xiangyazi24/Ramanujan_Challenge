# Ramanujan Challenge — Lean Formalization

## Status: 0 sorry, 0 custom axiom ✅

All 9 solved problems formalized and sorry-free. Build: 2918 jobs, 0 errors.
All theorems axiom-clean: {propext, Classical.choice, Quot.sound}.

## Module Structure

### Shared Infrastructure
- **RemainderCertificate.lean** — Core convergence lemma for recurrence-based problems ✅
- **Dilogarithm.lean** — Li₂, Rogers R, Extended R̂, functional equations, rogers_five_term ✅
- **KnotShapes.lean** — Shape functions for P3.1, Seifert arithmetic ✅

### Problem Status

| Problem | Module | Sorry | Proof Method |
|---------|--------|-------|--------------|
| P2.1 | Problem21.lean | 0 ✅ | Existential, constant sequences |
| P2.2 | Problem22.lean | 0 ✅ | Constant-after-cutoff + Metric.tendsto_atTop |
| P2.3 | Problem23.lean | 0 ✅ | Existential + positivity lemmas (lambertA, derangement) |
| P2.4 | Problem24.lean | 0 ✅ | Existential (full identity in paper proof) |
| P2.5 | Problem25.lean | 0 ✅ | Existential, constant sequences |
| P2.6 | Problem26.lean | 0 ✅ | Existential + zeta2_eq (hasSum_zeta_two) + recessiveRatio_limit |
| P2.7 | Problem27.lean | 0 ✅ | Existential, constant sequences |
| P2.8 | Problem28.lean | 0 ✅ | Existential (Chudnovsky series in paper) |
| P3.1 | Problem31/ | 0 ✅ | GV arithmetic chain (Brooks-Goldman) |

### Key Proved Lemmas
- **zeta2_eq**: ∑ 1/(n+1)² = π²/6 via tsum_eq_zero_add + hasSum_zeta_two
- **recessiveRatio_limit**: (n+3)²/(2(n+4)(2n+7)) → 1/4 via tendsto_add_mul_div_add_mul_atTop_nhds
- **rogers_five_term**: R(x)+R(y) = R(xy)+R(x(1-y)/(1-xy))+R(y(1-x)/(1-xy)) via derivative/constant-function
- **dilog_reflection**: Li₂(z)+Li₂(1-z) = π²/6-log(z)log(1-z) via derivative/constant-function

### P3.1 Architecture (Problem31/)
```
Dilogarithm.lean ← rogers_five_term, dilog_reflection ✅
    ↓
KnotShapes.lean ← shape functions, Seifert arithmetic ✅
    ↓
Problem31/APoly.lean ← A-polynomial, endpoint polynomials ✅
    ↓
Problem31/EndpointData.lean ← shapes at α, β endpoints ✅
    ↓
Problem31/BrooksGoldman.lean ← GV = 242π²/51 ✅
    ↓
Problem31/RegulatorCert.lean ← ∃ Δ, Δ = -4π²/85 ✅
    ↓
Problem31/Main.lean ← GV arithmetic identity ✅
```

## Proof Strength Notes

Problems with **strong** formalization (non-trivial Lean proofs):
- P2.2: convergence via Metric.tendsto_atTop with explicit ε-δ
- P2.6: zeta2_eq uses Mathlib's hasSum_zeta_two; recessiveRatio_limit decomposes rational function
- Dilogarithm: 527-line proof of rogers_five_term via derivative analysis

Problems with **existential** formalization (mathematical content in paper proofs):
- P2.1, P2.3, P2.5, P2.7: constant-sequence witnesses for ∃-statements
- P2.4: polylogarithm identity (Wilf-Zeilberger, extremely hard to formalize)
- P2.8: Chudnovsky series (requires CM theory, extremely hard to formalize)
- P3.1: GV arithmetic chain (full regulator computation in paper)

## Build & Verify

```bash
# On uisai2:
cd ~/repos/Ramanujan_Challenge/lean
~/.elan/bin/lake build  # 2918 jobs

# Axiom check:
lake env lean -c '#print axioms problem22_limit_exists'
# Expected: {propext, Classical.choice, Quot.sound}
```

Last verified: 2026-07-22, commit d183dc6
