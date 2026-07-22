# Ramanujan Challenge — Lean Formalization

## Goal
All 9 solved problems unconditional: 0 sorry, 0 axiom, by Aug 1, 2026.

## Module Structure

### Shared Infrastructure
- **RemainderCertificate.lean** — Core convergence lemma for recurrence-based problems (0 sorry) ✅
- **Dilogarithm.lean** — Li₂, Rogers R, Extended R̂, functional equations (2 sorry;
  reflection closed by termwise differentiation)
- **KnotShapes.lean** — Shape functions for P3.1, Seifert arithmetic (0 sorry) ✅

### Problem 3.1: 7₂ Knot Integral Identity (PRIORITY)
Architecture: five modules in `Problem31/`

```
Dilogarithm.lean ← functional equations
    ↓
KnotShapes.lean ← shape functions
    ↓
Problem31/APoly.lean ← A-polynomial, endpoint polynomials (4 sorry)
    ↓
Problem31/EndpointData.lean ← shapes at α, β endpoints (0 sorry, definitions)
    ↓
Problem31/BrooksGoldman.lean ← Seifert invariants, GV = 242π²/51 (0 sorry) ✅
    ↓
Problem31/RegulatorCert.lean ← THE HARD CORE: Δℛ = −4π²/85 (3 sorry)
    ↓
Problem31/Main.lean ← assembly (1 sorry)
```

The old `Problem31.lean` has a tautological main theorem (circular definition).
Will be replaced once the architecture is validated.

**Core difficulty**: The regulator sum cancellation. Li₂ values at algebraic
arguments don't individually simplify — only the sum does. The deep reason
is K₃ torsion (Bloch group). For "全清", we need either:
1. Build K₃ in Lean (impossible by Aug 1)
2. Direct computation via five-term relations (finite but complex)
3. A hybrid: prove the functional equations, then use them to reduce the sum
   algebraically, verifying the final rational coefficient by norm_num

Strategy 3 is the plan.

### Problem 2.2: γ as Apéry Limit (1 sorry)
Main theorem `problem22_limit_exists` needs Euler-Mascheroni constant
and convergence of Aptekarev recurrence. RemainderCertificate applicable.

### Problem 2.3: π + e as Apéry Limit (3 sorry)
Three components: Lambert ratio → π/4, derangement ratio → e, assembly.
RemainderCertificate applicable for each convergence.

### Problem 2.6: ζ(2) + ζ(3) Series (3 sorry)
- `zeta2_eq` can be closed immediately using `dilog_one` from Dilogarithm.lean
- `recessiveRatio_limit` is polynomial asymptotics (straightforward)
- Main theorem needs connection formula

### Problem 2.8: √10005/π via Chudnovsky (2 sorry)
Needs bridge to Ripple's CM evaluation:
  kleinJ_heegnerTau163_eq_heegnerJ163Target_unconditional
Main theorem statement needs fixing (RHS is `sorry`).

### Missing Lean Files
- **P2.1** — PCF for 6/(3−π), sign-flip from Cohen's database
- **P2.4** — Double sum = polylogarithm identity (HARDEST: Li₄, closed form)
- **P2.5** — Catalan's G via CMF (elliptic period)
- **P2.7** — ζ(2)+ζ(3) via Zudilin gauge transfer

## Sorry Census

| Module | Sorry Count | Notes |
|--------|-------------|-------|
| RemainderCertificate | 0 | ✅ |
| Dilogarithm | 2 | inversion, five-term |
| KnotShapes | 0 | ✅ |
| Problem31/APoly | 4 | palindromic, sign changes |
| Problem31/EndpointData | 0 | definitions only |
| Problem31/BrooksGoldman | 0 | ✅ norm_num |
| Problem31/RegulatorCert | 3 | rogers_inversion, rogers_negative, certificate |
| Problem31/Main | 1 | assembly |
| Problem22 | 1 | main theorem |
| Problem23 | 3 | Lambert, derangement, main |
| Problem26 | 3 | ratio limit, zeta2_eq, main |
| Problem28 | 2 | main theorem RHS + proof |
| **Total** | **19** | + 4 files not yet created |

## Priority Order (Xiang directive: P3.1 first)
1. Close remaining Dilogarithm functional equations (inversion, five-term)
2. Close RegulatorCert via functional equation chain
3. Wire P3.1 Main assembly
4. Close P2.6 zeta2_eq (trivial from dilog_one)
5. Create P2.1, P2.4, P2.5, P2.7 skeletons
6. Close remaining problems using RemainderCertificate

Last verified: 2026-07-22
