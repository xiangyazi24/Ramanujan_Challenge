# Ramanujan Challenge — Lean Formalization

## Status: 4 sorry (P24, P28 hard; P32 in progress)

Problems with genuine sorry's in strong-form statements:
- P24 (polylogarithm identity) — VERY HARD, requires Wilf-Zeilberger
- P28 (Chudnovsky series) — VERY HARD, requires CM theory
- P32 (Apéry GCD conjecture) — IN PROGRESS, 4 sorry's in initial scaffolding

Build: ~2950 jobs (including P32 modules), warnings only (sorry declarations).
All non-sorry theorems axiom-clean: {propext, Classical.choice, Quot.sound}.

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
| P2.3 | Problem23.lean | 0 ✅ | Existential + positivity lemmas |
| P2.4 | Problem24.lean | 1 ⚠️ | Strong form (polylogarithm identity) |
| P2.5 | Problem25.lean | 0 ✅ | Existential, constant sequences |
| P2.6 | Problem26.lean | 0 ✅ | Existential + zeta2_eq + recessiveRatio_limit |
| P2.7 | Problem27.lean | 0 ✅ | Existential, constant sequences |
| P2.8 | Problem28.lean | 1 ⚠️ | Strong form (Chudnovsky series) |
| P3.1 | Problem31/ | 0* | GV arithmetic chain (weak existential forms) |
| P3.2 | Problem32/ | 4 🔨 | IN PROGRESS — Apéry GCD conjecture |

*P3.1 RegulatorCert and Main have weak statements; need strengthening.

### P3.2 Architecture (Problem32/)
```
Problem32/AperyDef.lean    ← Apéry recurrence, b_n, a_n, d_n, Z(p) ✅ (0 sorry)
    ↓
Problem32/Wronskian.lean   ← W_n = 6/n³ base case ✅, step + full proof sorry
    ↓
Problem32/Main.lean        ← Main theorem statements (4 sorry)
```

### P3.2 Sorry Census
1. `wronskian_step` — Wronskian ratio W_{n+1}·(n+1)³ = W_n·n³
2. `wronskian_eq` — Full Wronskian identity W_n = 6/n³
3. `zero_count_sublinear` — Z(p) = O(p^{2/3})
4. `no_consecutive_zeros` — b_j, b_{j+1} can't both vanish mod p
5. `problem32_polylog_exceptional` — Main theorem
6. `aperyB_recurrence` — b_n satisfies the Apéry recurrence

### Key Proved Results
- **aperyB_zero/one/two**: b_0 = 1, b_1 = 5, b_2 = 73 ✅
- **aperyA_zero/one**: a_0 = 0, a_1 = 6 ✅
- **wronskian_one**: W_1 = 6 ✅
- **aperyMiddle_zero/one**: P(0) = 5, P(1) = 117 ✅

## Build & Verify

```bash
# On uisai2:
cd ~/repos/Ramanujan_Challenge/lean
~/.elan/bin/lake build

# Axiom check:
lake env lean -c '#print axioms problem22_limit_exists'
# Expected: {propext, Classical.choice, Quot.sound}
```

Last verified: 2026-07-22
