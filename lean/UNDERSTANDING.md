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
- **Genuinely proved**: P2.3 (full chain modulo one cited classical input),
  P2.6 (unconditional full chain), plus real auxiliary content in P2.1/2.2
  and the P3.1 and P3.2 layers.
- **Not formalized**: the main limits of P2.1, P2.2, P2.5, P2.7.
  P2.5 now has a substantial exact algebraic layer, but the Catalan connection
  coefficient is still missing.

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
| P2.1 | Problem21.lean | 0 ✅ | **full chain**: Cohen Entry 5.3.22 sign-flip, even contraction, moment integral J_n→0, unconditional `problem21_pcf_value`. |
| P2.2 | Problem22.lean | 0 | the four shift identities `c̃_j(m) = c_j(m−3)` (`ring`). **Limit not formalized.** |
| P2.3 | Problem23.lean | 0 ✅ | **full chain**: tensor-annihilation theorem for arbitrary solutions, closed forms solve the recurrence, all 8 initial values, `c₀ ≠ 0` + uniqueness, exact ratio splitting, `m!/D_m → e` from Mathlib, main theorem — conditional only on Lambert's CF value, which is an explicit hypothesis |
| P2.4 | Problem24.lean | 1 ⚠️ | strong form (polylogarithm identity) left as `sorry` |
| P2.5 | Problem25.lean + Problem25Connection.lean + Problem25Moment.lean | 0 | ~2400 lines across 3 files. CMF transcription, row recurrences, det factorization, Pochhammer gauge, Poincaré cubic, sign conjugation, positive cocycle, nonvanishing denominators, convex-hull nesting, Catalan error recurrence + initial signs. **Problem25Connection.lean** (757 lines): unconditional projective contraction → all 3 ratios converge to `commonLimit` (2/3 geometric rate). `Problem25Claim ↔ commonLimit = catalanConstant`. **Problem25Moment.lean** (WIP): moment formulas ∫wq·R = Q, ∫wp·R = P proved, remainder polynomial defined. **Gap: Catalan integral representation G = -∫log(t)/(1+t²)dt + remainder integral identity → commonLimit = G.** |
| P2.6 | Problem26.lean + 5 aux files | 0 ✅ | **Complete unconditional formalization** (6246 lines). Exact recurrence and initial values, uniqueness, Ore factorization and reduction of order, summability, generating-function/integral bridge, weight-2 and nested weight-3 inverse-binomial evaluations, and the cyclotomic logarithmic integral are all proved. `problem26_hasSum_of_spec` states convergence and the value for every rational sequence satisfying the printed recurrence and initial values. Lean v4.29; axiom audit: `{propext, Classical.choice, Quot.sound}`. |
| P2.7 | Problem27.lean + Problem27Barnes.lean | 0 | ~2800 lines. Full 3×3 rational gauge (column identities all proved), analytic transfer complete, Barnes midpoint ≤ 2·(1/4)^n, Euler cosh product, integrable envelope, **Barnes error integral → 0 proved**. **Gap: source normalization (Zudilin recurrence error = Barnes contour integral).** |
| P2.8 | **Problem28_SUBMIT/** | 0 ✅ | **Complete formalization synced from SUBMIT** (30080 lines, 27 files, Lean v4.29). Full Ripple CM extraction: j(τ₁₆₃) = −640320³, Chudnovsky hypergeometric, period bridge. `Problem28.lean` (107 lines, 1 sorry) is the old stub. |
| P3.1 | **Problem31_SUBMIT/** | 0 ✅ | **Complete formalization synced from SUBMIT** (4724 lines, 28 files, Lean v4.29). Bloch-Wigner four-shape vanishing, rational reconstruction. `Problem31/` is the old stub. See `SUBMIT/3.1/lean/UNDERSTANDING.md` for atom checklist. |
| P3.2 | Problem32/ | 3 🔨 | IN PROGRESS — Apéry GCD conjecture |

**Synced SUBMIT packages (2026-08-01):** `Problem28_SUBMIT/` and `Problem31_SUBMIT/` are verbatim copies of `SUBMIT/2.8/lean/Ripple/` and `SUBMIT/3.1/lean/Ramanujan31/` respectively. Both are complete, 0-sorry formalizations on Lean v4.29. The old `Problem28.lean` and `Problem31/` stubs remain for reference but are superseded.

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

## P2.5 AZ-certificate audit (2026-07-30)

The polynomial certificate produced by `scripts/p25_az_correct.py` does verify
the cleared differential identity

```text
Phi P_(n+1) - P_n M_H = D_u A + D_v B + D_t C.
```

It does **not** yield the Delannoy-square period claimed in the exploratory
notes.  The covariant derivatives are divergence terms for the ordinary form
`du dv dt`.  After the factor `u v (1+t^2)` in `P` cancels the carrier
denominator, small `u`- and `v`-contours extract the coefficient of
`u^-1 v^-1`, not the constant term.  Replacing these by logarithmic measures
`du/u` and `dv/v` changes the functional and invalidates the integration-by-
parts step.  At `eps = 0` the remaining `t`-integrand is constant, so its
ordinary closed-contour integral is zero.

There is also a direct exact recurrence countercheck: the proposed carrier row
`r_0 = [1, 1, 1/4]` gives

```text
r_0 M_H(0) = [15/8, -11/8, 1/4],
```

not `[81/4, 27/2, 9/4]`.  Therefore this certificate must not be cited as an
exact connection coefficient or a proof of the Catalan limit.  The same warning
is recorded at the top of `scripts/p25_az_verify.py`.
