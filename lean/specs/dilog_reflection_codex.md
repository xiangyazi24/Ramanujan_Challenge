# Codex Spec: Close `dilog_reflection` in Dilogarithm.lean

## File
`~/Ramanujan_Challenge/lean/RamanujanChallenge/Dilogarithm.lean`

## Target
```lean
theorem dilog_reflection {z : ℝ} (hz0 : 0 < z) (hz1 : z < 1) :
    dilog z + dilog (1 - z) = Real.pi ^ 2 / 6 - Real.log z * Real.log (1 - z)
```

where `dilog z = ∑' n : ℕ, z ^ (n + 1) / (↑(n + 1) : ℝ) ^ 2`.

## Already proved in the file (DO NOT modify these)
- `dilog_zero : dilog 0 = 0`
- `dilog_summable {z : ℝ} (hz : |z| ≤ 1) : Summable ...`
- `dilog_one : dilog 1 = Real.pi ^ 2 / 6`
- `summable_one_div_succ_sq : Summable (fun n => 1 / (↑(n+1))^2)`

## Proof strategy

**Approach A (derivative + constant function):**
1. Define f(z) = dilog(z) + dilog(1-z) + log(z)·log(1-z) on (0,1)
2. Show f has derivative 0 on (0,1):
   - d/dz[dilog(z)] = -log(1-z)/z  (term-by-term differentiation of power series)
   - d/dz[dilog(1-z)] = log(z)/(1-z)  (chain rule)
   - d/dz[log(z)·log(1-z)] = log(1-z)/z - log(z)/(1-z)  (product rule)
   - Sum = 0
3. f is continuous on (0,1), f'=0 implies f is constant
4. Show f(z) → π²/6 as z → 1⁻ (using dilog_one, dilog_zero, and log(1)·log(0⁺) → 0)
5. Therefore f(z) = π²/6 for all z ∈ (0,1)
6. Rearrange: dilog(z) + dilog(1-z) = π²/6 - log(z)·log(1-z)

Key Mathlib tools:
- `hasDerivAt_tsum` or `HasFPowerSeriesAt` for term-by-term differentiation
- `Real.hasDerivAt_log` for log derivative
- `is_const_of_derivWithin_eq_zero` or `isConstOfDerivWithinEqZero`
- Power series for -log(1-z) = ∑ z^n/n: `Real.hasSum_log_one_sub_of_abs_lt` or similar

**Approach B (if approach A is too hard):**
Use the integral representation directly:
- dilog(z) = -∫₀ᶻ log(1-t)/t dt
- This is provable from the power series by term-by-term integration
- Then dilog(z) + dilog(1-z) can be computed by substitution in the integral

## Constraints
- No sorry, no axiom, no native_decide
- Must not modify existing proved theorems
- Lean 4 v4.30.0, Mathlib v4.30.0
- Build command: `source ~/.profile && cd ~/Ramanujan_Challenge/lean && lake env lean RamanujanChallenge/Dilogarithm.lean`
- If `lake env lean` gives stale errors, use `lake build RamanujanChallenge.Dilogarithm`
- Work only in `~/Ramanujan_Challenge/lean/`

## What does NOT count as done
- A proof that uses sorry or axiom
- Changing the theorem statement
- Breaking existing theorems
- A proof that only works for specific values of z
