# Spec: Close `dilog_reflection` in Dilogarithm.lean

## Target
Prove the following in `RamanujanChallenge/Dilogarithm.lean`:
```lean
theorem dilog_reflection {z : ℝ} (hz0 : 0 < z) (hz1 : z < 1) :
    dilog z + dilog (1 - z) = Real.pi ^ 2 / 6 - Real.log z * Real.log (1 - z)
```

where `dilog z = ∑' n : ℕ, z ^ (n + 1) / (↑(n + 1) : ℝ) ^ 2`.

## Proof Strategy

Use the derivative approach:
1. Define `f(z) = dilog z + dilog (1-z) + log(z) * log(1-z)` for 0 < z < 1.
2. Show `f'(z) = 0` on (0,1):
   - `d/dz [dilog z] = -log(1-z)/z` (from term-by-term differentiation of the power series)
   - `d/dz [dilog(1-z)] = log(z)/(1-z)` (chain rule)
   - `d/dz [log(z)*log(1-z)] = log(1-z)/z - log(z)/(1-z)` (product rule)
   - Sum = 0
3. `f` is continuous on (0,1) and constant (by step 2).
4. Evaluate: as z → 1⁻, f(z) → dilog(1) + dilog(0) + 0 = π²/6.
   (Use `dilog_one` and `dilog_zero` already proved in the file.)
5. Rearrange: `dilog z + dilog (1-z) = π²/6 - log(z)*log(1-z)`.

## Key Mathlib lemmas needed

- `hasDerivAt_tsum` or term-by-term differentiation of power series
- `Real.hasDerivAt_log` : derivative of log
- `Real.log_mul_continuous_on` or product rule for log*log
- `IsConst.eq` or `isConstOn_of_derivWithin_eq_zero` : f' = 0 → f constant
- `dilog_one` and `dilog_zero` (already proved in the file)
- `-log(1-z) = ∑ z^n/n` for |z| < 1 (Mathlib: `hasSum_pow_div_log_of_abs_lt_one` or similar)

## Alternative approach (if derivative approach is too hard)

Use the integral representation:
- `dilog z = -∫_0^z (log(1-t)/t) dt`
- `dilog(1-z) = -∫_0^{1-z} (log(1-t)/t) dt`
- Change variables in the second integral: t → 1-u
- Combine and use `∫_0^1 (log(1-t)/t) dt = -π²/6`

## Hard constraints
- No sorry, no axiom, no native_decide
- Must build with `lake env lean RamanujanChallenge/Dilogarithm.lean` on Lean 4 v4.30.0 + Mathlib v4.30.0
- Work only in `~/repos/Ramanujan_Challenge/lean/`
- The existing proved theorems (`dilog_zero`, `dilog_summable`, `dilog_one`, `rogers_zero`, `rogers_one`) must remain proved

## What does NOT count as done
- A proof that uses `sorry` or `axiom`
- A proof that changes the theorem statement
- A proof that breaks other theorems in the file
