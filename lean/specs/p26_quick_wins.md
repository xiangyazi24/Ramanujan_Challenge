# Spec: Close zeta2_eq and recessiveRatio_limit in Problem26.lean

## Working directory
`~/Ramanujan_Challenge/lean`

## Target 1: zeta2_eq (line 88-89)

```lean
theorem zeta2_eq : ∑' n : ℕ, (1 : ℝ) / (↑n + 1) ^ 2 = Real.pi ^ 2 / 6 := by
  sorry
```

### Proof strategy
Use `hasSum_zeta_two : HasSum (fun n => 1 / (n : ℝ) ^ 2) (π²/6)` from Mathlib.

Key insight: the n=0 term of the zeta-two series is `1/0² = 0`, so:
∑' n, 1/n² = 0 + ∑' n, 1/(n+1)²

Use `Summable.tsum_eq_zero_add` to split off the n=0 term:
```lean
have hs := hasSum_zeta_two.summable
have h := hs.tsum_eq_zero_add  -- ∑' n, f n = f 0 + ∑' n, f(n+1)
```
Then simplify f(0) = 1/0² = 0, and show the shifted sum equals our target.

Note: `(n+1 : ℕ) : ℝ` needs `Nat.cast_succ` or `push_cast` to match `↑n + 1`.

## Target 2: recessiveRatio_limit (line 67-69)

```lean
def recessiveRatio (n : ℕ) : ℚ :=
  (↑n + 3) ^ 2 / (2 * (↑n + 4) * (2 * ↑n + 7))

theorem recessiveRatio_limit :
    Filter.Tendsto recessiveRatio Filter.atTop (nhds (1 / 4 : ℚ)) := by
  sorry
```

### Proof strategy
This is a rational function of n with leading coefficient ratio = 1/4:
- Numerator: n² + 6n + 9 (degree 2, leading coeff 1)
- Denominator: 4n² + 22n + 56 (degree 2, leading coeff 4)

Standard approach: divide numerator and denominator by n², take limits.

In Lean/Mathlib, use `Filter.Tendsto.div` and show the numerator/n² → 1
and denominator/n² → 4. Or use `Tendsto.ratCast_atTop` type lemmas.

Alternative: show `recessiveRatio n = 1/4 + (negative terms)/denominator` and
the error term → 0.

## BUILD COMMAND
```
rsync -az --exclude=.lake --exclude=lake-packages . uisai2:~/Ramanujan_Challenge/lean/
ssh uisai2 'cd ~/Ramanujan_Challenge/lean && source ~/.profile 2>/dev/null && lake build RamanujanChallenge.Problem26 2>&1 | tail -20'
```

## CONSTRAINTS
- No sorry, no axiom
- Lean 4 v4.30.0, Mathlib v4.30.0
- Do NOT modify any other file
- Only edit Problem26.lean
