# Codex Spec: Close `recessiveRatio_limit` in Problem26.lean

## Target file
`RamanujanChallenge/Problem26.lean` — lines 67-69

## Current state
```lean
def recessiveRatio (n : ℕ) : ℚ :=
  (↑n + 3) ^ 2 / (2 * (↑n + 4) * (2 * ↑n + 7))

theorem recessiveRatio_limit :
    Filter.Tendsto recessiveRatio Filter.atTop (nhds (1 / 4 : ℚ)) := by
  sorry
```

## Proof strategy

This is a rational function limit: (n+3)²/(2(n+4)(2n+7)) → 1/4 as n → ∞.

Expand: numerator = n² + 6n + 9, denominator = 4n² + 22n + 56.

Approach 1 (algebraic decomposition): Show recessiveRatio n = 1/4 + error(n)
where error(n) → 0. Concretely:
  recessiveRatio n - 1/4 = [(n+3)² - (n+4)(2n+7)/2] / [2(n+4)(2n+7)]
  = [n² + 6n + 9 - (2n² + 15n + 28)/2] / [2(n+4)(2n+7)]
  = [2n² + 12n + 18 - 2n² - 15n - 28] / [4(n+4)(2n+7)]
  = [-3n - 10] / [4(n+4)(2n+7)]

So |recessiveRatio n - 1/4| = (3n + 10) / (4(n+4)(2n+7)).

For large n, this is O(1/n) → 0.

Approach 2 (Tendsto API): Use `Filter.Tendsto.div` after showing
numerator/n² → 1 and denominator/n² → 4. This requires working with
`Filter.Tendsto` and `nhds` in ℚ.

Approach 3 (squeeze): Show |recessiveRatio n - 1/4| ≤ C/n for some C,
then use `tendsto_const_div_atTop_nhds_0_nat`.

NOTE: This is over ℚ (not ℝ). The `nhds` topology on ℚ comes from the
absolute value. Use `Metric.tendsto_atTop` or `NormedAddCommGroup` structure.

Mathlib lemmas that may help:
- `Filter.Tendsto.div`
- `Rat.cast_injective` to transfer from ℝ if needed
- `tendsto_const_div_atTop_nhds_0_nat` (might need ℝ version)
- Consider casting to ℝ, proving the limit there, then pulling back

## Build and verify
```bash
cd ~/Ramanujan_Challenge/lean
source ~/.profile 2>/dev/null
lake build RamanujanChallenge.Problem26 2>&1 | tail -20
```

## Constraints
- No sorry, no axiom, no native_decide
- Only modify the `recessiveRatio_limit` proof (line 69)
- Do not change any definitions or other theorems
- Lean 4 v4.30.0, Mathlib v4.30.0
