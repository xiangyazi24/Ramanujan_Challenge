# Ramanujan Challenge — Lean Formalization

## ⚠️ Read this first: "0 sorry" was NOT the same as "proved" (corrected 2026-07-30)

Until 2026-07-30 this table reported P2.1/2.2/2.3/2.5/2.6/2.7 as "0 sorry ✅".
That was misleading. Their main theorems were **vacuous**: statements of the form

```lean
theorem problem25_identity : ∃ (p q : ℕ → ℝ), Tendsto (fun n => p n / q n) atTop (𝓝 L) :=
  ⟨fun _ => L, fun _ => 1, by simp⟩
```

witnessed by constant sequences, saying nothing about the challenge recurrence.
P2.2's version was worse — it asserted the initial values and then converged to
`179/306`, not γ. `Problem21.lean` also carried `sign_flip_P … : True := trivial`
with unused hypotheses.

**All of these have been deleted** (not weakened, not `sorry`-ed) and each file
now carries an explicit note saying the limit is not formalized there. A `sorry`
is honest; a vacuous existential is not, and shipping one in a submission would
read as an attempt to fake a formalization.

P2.3 has been **rebuilt with real content** — see below.

## Status

- **Genuine `sorry`s** (declared-hard statements): P24 (1), P28 (1), P32 (3).
- **Genuinely proved**: P2.3 (full chain modulo one cited classical input), plus
  real auxiliary content in P2.1/2.2/2.6 and the P3.1 and P3.2 layers.
- **Not formalized**: the main limits of P2.1, P2.2, P2.5, P2.6, P2.7.

Build: 3427 jobs, 0 errors; warnings are exactly the five `sorry` declarations.
All non-`sorry` theorems axiom-clean: {propext, Classical.choice, Quot.sound}.

## Module Structure

### Shared Infrastructure
- **RemainderCertificate.lean** — Core convergence lemma for recurrence-based problems ✅
- **Dilogarithm.lean** — Li₂, Rogers R, Extended R̂, functional equations, rogers_five_term ✅
- **KnotShapes.lean** — Shape functions for P3.1, Seifert arithmetic ✅

### Problem Status

| Problem | Module | Sorry | What is actually proved |
|---------|--------|-------|--------------------------|
| P2.1 | Problem21.lean | 0 | coefficient definitions + the index-shift identities `a_n = −α(n+1)`, `b_n = β(n)` (`ring`). **Limit not formalized.** |
| P2.2 | Problem22.lean | 0 | the four shift identities `c̃_j(m) = c_j(m−3)` (`ring`). **Limit not formalized.** |
| P2.3 | Problem23.lean | 0 ✅ | **full chain**: tensor-annihilation theorem for arbitrary solutions, closed forms solve the recurrence, all 8 initial values, `c₀ ≠ 0` + uniqueness, exact ratio splitting, `m!/D_m → e` from Mathlib, main theorem — conditional only on Lambert's CF value, which is an explicit hypothesis |
| P2.4 | Problem24.lean | 1 ⚠️ | strong form (polylogarithm identity) left as `sorry` |
| P2.5 | Problem25.lean | 0 | only the definition of Catalan's constant. **Nothing substantive.** |
| P2.6 | Problem26.lean | 0 | `recessiveRatio_limit` (ratio → 1/4) and `zeta2_eq` (ζ(2) = π²/6). **Identity not formalized.** |
| P2.7 | Problem27.lean | 0 | nothing beyond documentation. **Nothing substantive.** |
| P2.8 | Problem28.lean | 1 ⚠️ | strong form (Chudnovsky series) left as `sorry`; the *submitted* 2.8 package uses the separate Ripple extraction, which is far stronger |
| P3.1 | Problem31/ | 0* | GV arithmetic chain |
| P3.2 | Problem32/ | 3 🔨 | IN PROGRESS — Apéry GCD conjecture |

*P3.1 RegulatorCert and Main have weak statements; need strengthening.

**Rule going forward:** a formal statement counts only if it mentions the
challenge's own recurrence or sequences. An `∃ p q, …` with free `p, q` is not a
formalization of anything.

### P3.2 Architecture (Problem32/)
```
Problem32/AperyDef.lean    ← Apéry recurrence, b_n, a_n, d_n, Z(p) ✅ (0 sorry)
    ↓
Problem32/Wronskian.lean   ← W_n = 6/n³ ✅ (0 sorry, all axiom-clean)
    ↓
Problem32/Main.lean        ← Main theorem (3 sorry)
```

### P3.2 Sorry Census (3 sorry)
1. `aperyB_recurrence_int` — b_n (closed form) satisfies the Apéry recurrence (WZ identity)
2. `zero_count_sublinear` — Z(p) = O(p^{2/3}) (gap polynomial argument)
3. `problem32_polylog_exceptional` — Main theorem

### P3.2 Proved Results
- **aperyB_zero/one/two**: b_0 = 1, b_1 = 5, b_2 = 73 ✅
- **aperyA_zero/one**: a_0 = 0, a_1 = 6 ✅
- **aperyA/BQ_recurrence**: definitional from recurrence construction ✅
- **wronskian_one**: W_1 = 6 ✅
- **wronskian_step**: (n+1)³W_{n+1} = n³W_n via `linear_combination` ✅
- **wronskian_mul**: W_n·n³ = 6 by induction ✅
- **wronskian_eq**: W_n = 6/n³ ✅
- **no_consecutive_zeros**: b_j, b_{j+1} can't both vanish mod p (modulo aperyB_recurrence_int) ✅
- **aperyP_zero/one**: P(0) = 5, P(1) = 117 ✅

## Build & Verify

```bash
# On uisai2:
cd ~/repos/Ramanujan_Challenge/lean
~/.elan/bin/lake build

# Axiom check:
lake env lean -c '#print axioms problem22_limit_exists'
# Expected: {propext, Classical.choice, Quot.sound}
```

Last verified: 2026-07-22 (Wronskian layer complete, no_consecutive_zeros proved)
