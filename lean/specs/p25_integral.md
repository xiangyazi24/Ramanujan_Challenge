# Task: Prove Catalan integral identity in Lean 4

## Context
Working in this repo at `lean/`. Lean v4.29.0 + Mathlib v4.29.0.

`catalanConstant` is defined in `RamanujanChallenge/Problem25.lean`:
```lean
def catalanConstant : ℝ := ∑' n : ℕ, (-1 : ℝ) ^ n / (2 * (n : ℝ) + 1) ^ 2
```

Also proved there: `catalanSeries_summable` and `catalanPartialSum_tendsto`.

## Goal
Create `RamanujanChallenge/Problem25Integral.lean` proving:

```lean
theorem catalanConstant_eq_integral :
    catalanConstant = ∫ t in (0:ℝ)..1, (-Real.log t) / (1 + t ^ 2)
```

## Strategy

### Step 1: ∫₀¹ t^k · (-log t) dt = 1/(k+1)²
Use integration by parts: u = -log t, dv = t^k dt.
Mathlib has `integral_log_from_zero : ∫ s in 0..b, log s = b * log b - b`.
Also `hasDerivAt_mul_log` and FTC tools.

### Step 2: dominated convergence swap
Partial sums S_N(t) = Σ_{n<N} (-1)^n t^{2n} converge to 1/(1+t²) pointwise.
|S_N(t)| ≤ 1 for t ∈ [0,1]. And (-log t) is integrable on [0,1].
So ∫(-log t)·S_N → ∫(-log t)/(1+t²) by dominated convergence.
And ∫(-log t)·S_N = Σ_{n<N} (-1)^n/(2n+1)² → catalanConstant.

## Build
```bash
~/.elan/toolchains/leanprover--lean4---v4.29.0/bin/lake build RamanujanChallenge.Problem25Integral
```

## Rules
- Create only `RamanujanChallenge/Problem25Integral.lean`
- 0 sorry, 0 errors
- Use `set_option maxHeartbeats 0` if needed
- Import `RamanujanChallenge.Problem25` and relevant Mathlib modules
