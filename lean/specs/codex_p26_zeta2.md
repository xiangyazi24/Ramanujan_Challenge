# Codex Spec: Close `zeta2_eq` in Problem26.lean

## Target file
`RamanujanChallenge/Problem26.lean` — line 88-89

## Current state
```lean
theorem zeta2_eq : ∑' n : ℕ, (1 : ℝ) / (↑n + 1) ^ 2 = Real.pi ^ 2 / 6 := by
  sorry
```

## Available imports (already in the file)
```lean
import Mathlib.NumberTheory.ZetaValues  -- has hasSum_zeta_two
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
```

## Proof strategy

`hasSum_zeta_two` from Mathlib gives:
```
HasSum (fun n : ℕ => 1 / (↑n : ℝ) ^ 2) (π²/6)
```

The target sum is shifted: `∑' n, 1/(n+1)²` vs `∑' n, 1/n²`.
They differ by the n=0 term, which is `1/0² = 0` (since `zero_pow two_ne_zero` gives `0^2 = 0`, so `1/0 = 0` by `div_zero`).

Approach: use `hasSum_zeta_two.tsum_eq` to get `∑' n, 1/n² = π²/6`, then
use the index shift to relate the two sums.

Specifically: `∑' n, 1/n² = 1/0² + ∑' n, 1/(n+1)²`. Since 1/0² = 0,
we get `∑' n, 1/(n+1)² = ∑' n, 1/n² = π²/6`.

Key Mathlib API:
- `hasSum_zeta_two : HasSum (fun n : ℕ => 1 / (↑n : ℝ) ^ 2) (π²/6)`
- `HasSum.tsum_eq` to get the tsum value
- `Summable.tsum_eq_zero_add` or `tsum_eq_zero_add` to split off n=0 term
- `Nat.cast_zero`, `zero_pow`, `div_zero` for simplification

Alternative: use `Equiv.sumCompl` or just show the two functions are equal
up to a shift by `tsum_eq_tsum_of_hasSum_iff_hasSum`.

## Build and verify
```bash
cd ~/Ramanujan_Challenge/lean
source ~/.profile 2>/dev/null
lake build RamanujanChallenge.Problem26 2>&1 | tail -20
```

## Constraints
- No sorry, no axiom, no native_decide
- Only modify the `zeta2_eq` proof (lines 88-89)
- Do not change any definitions or other theorems
- Lean 4 v4.30.0, Mathlib v4.30.0
