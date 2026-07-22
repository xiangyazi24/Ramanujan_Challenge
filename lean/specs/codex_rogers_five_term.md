# Codex Spec: Close `rogers_five_term` in Dilogarithm.lean

## Target file
`RamanujanChallenge/Dilogarithm.lean` — lines 261-267

## Current state
```lean
theorem rogers_five_term {x y : ℝ}
    (hx0 : 0 < x) (hx1 : x < 1) (hy0 : 0 < y) (hy1 : y < 1) :
    rogersDialogarithm x + rogersDialogarithm y =
      rogersDialogarithm (x * y) +
      rogersDialogarithm (x * (1 - y) / (1 - x * y)) +
      rogersDialogarithm (y * (1 - x) / (1 - x * y)) := by
  sorry
```

where
```lean
def rogersDialogarithm (z : ℝ) : ℝ :=
  dilog z + (1 / 2) * Real.log z * Real.log (1 - z)
```

## Already proved infrastructure (in the same file)

```lean
-- dilog and its derivative on (0,1)
private theorem dilog_hasDerivAt {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasDerivAt dilog (-(Real.log (1 - x)) / x) x

-- dilog reflection (also proved via derivative/constant approach)
theorem dilog_reflection {z : ℝ} (hz0 : 0 < z) (hz1 : z < 1) :
    dilog z + dilog (1 - z) = Real.pi ^ 2 / 6 - Real.log z * Real.log (1 - z)
```

## Proof strategy

Use the SAME approach that proved `dilog_reflection`: the derivative/constant-function method.

1. Define an auxiliary function of ONE variable (fix one of x, y):
   F(x) = LHS(x,y) - RHS(x,y) for fixed y ∈ (0,1).

2. Show F'(x) = 0 for all x ∈ (0,1):
   - Compute derivatives of all five Rogers terms w.r.t. x
   - The Rogers derivative is: R'(z) = dilog'(z) + ½[log'(z)·log(1-z) + log(z)·log'(1-z)]
     = -log(1-z)/z + ½[log(1-z)/z - log(z)/(1-z)]
     = -½·log(1-z)/z - ½·log(z)/(1-z)
   - Use chain rule for the composed terms: R(xy), R(x(1-y)/(1-xy)), R(y(1-x)/(1-xy))
   - After cancellation, F'(x) = 0

3. Show F is continuous on [a, 1] for some 0 < a < 1 (use ContinuousOn).

4. By `constant_of_has_deriv_right_zero`, F is constant on [a, 1].

5. Evaluate at x = 0⁺ (or take limit): as x → 0, Rogers(x) → 0, Rogers(xy) → 0,
   Rogers(x(1-y)/(1-xy)) → 0, Rogers(y(1-x)/(1-xy)) → Rogers(y), so F → 0.
   Or evaluate at some special point where the identity is easy to verify.

Alternative: Evaluate at x → 0⁺. As x → 0:
- R(x) → R(0) = 0
- R(xy) → R(0) = 0  
- R(x(1-y)/(1-xy)) → R(0) = 0
- R(y(1-x)/(1-xy)) → R(y)
So LHS → 0 + R(y) = R(y) and RHS → 0 + 0 + R(y) = R(y). F → 0.

Already proved: `rogers_zero : rogersDialogarithm 0 = 0`

## Key technical challenges
- The derivative computation has many terms and requires careful field_simp/ring
- Need ContinuousOn for Rogers on appropriate intervals
- The arguments x(1-y)/(1-xy) and y(1-x)/(1-xy) need to be shown in (0,1)
  when x,y ∈ (0,1): since 0 < 1-xy (because xy < 1), the arguments are positive
  and less than 1.

## Imports available
The file already imports everything needed (calculus, log, trig, zeta values).

## Build and verify
```bash
cd ~/Ramanujan_Challenge/lean
source ~/.profile 2>/dev/null
lake build RamanujanChallenge.Dilogarithm 2>&1 | tail -20
```

## Constraints
- No sorry, no axiom, no native_decide
- Only modify the `rogers_five_term` proof (line 267)
- May add private helper lemmas above the theorem (between line 258 and 261)
- Do not change any existing definitions or proved theorems
- Lean 4 v4.30.0, Mathlib v4.30.0
