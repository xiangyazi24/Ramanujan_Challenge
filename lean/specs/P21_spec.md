# Codex Spec: Problem 2.1 (6/(3-π)) — Full Lean Formalization

## Target
Create `RamanujanChallenge/Problem21.lean` with an unconditional proof that the PCF
with a_n = -220n³ - 484n² - 301n - 42 and b_n = 4n²(2n+1)²(5n-4)(5n+6)
converges to 6/(3-π).

## Working directory
`~/repos/Ramanujan_Challenge/lean/`

## Mathematical proof (from proof.tex)

The proof has three steps:
1. **Index shift:** a_n = -α(n+1) and b_n = β(n), where α(n) = 220n³-176n²-7n+5
   and β(n) = 4n²(2n+1)²(5n-4)(5n+6) are Cohen's CF coefficients for π.
2. **Sign-flip lemma:** If a CF with partial quotients c_k and partial numerators d_k
   converges to S, then the CF with -c_k and d_k converges to -S.
   Proof: by induction, P̃_n = (-1)^{n+1} P_n and Q̃_n = (-1)^n Q_n.
3. **Cohen's CF:** π = 3 + 6/(α(1) + β(1)/(α(2) + ...)).
   So the first tail T = 6/(π-3), and our CF = -T = 6/(3-π).

## Lean structure

```lean
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Data.Int.Basic

-- Definitions
def challenge_a (n : ℤ) : ℤ := -220 * n ^ 3 - 484 * n ^ 2 - 301 * n - 42
def challenge_b (n : ℤ) : ℤ := 4 * n ^ 2 * (2 * n + 1) ^ 2 * (5 * n - 4) * (5 * n + 6)
def cohen_alpha (n : ℤ) : ℤ := 220 * n ^ 3 - 176 * n ^ 2 - 7 * n + 5
def cohen_beta (n : ℤ) : ℤ := 4 * n ^ 2 * (2 * n + 1) ^ 2 * (5 * n - 4) * (5 * n + 6)

-- Step 1: polynomial identity (by ring)
theorem shift_identity (n : ℤ) : challenge_a n = -cohen_alpha (n + 1) := by ring
theorem beta_identity (n : ℤ) : challenge_b n = cohen_beta n := by ring

-- Step 2: General sign-flip lemma for continued fractions
-- Define CF convergents P_n, Q_n by recurrence:
--   P_n = c_n * P_{n-1} + d_n * P_{n-2}, P_{-1}=1, P_0=c_0
--   Q_n = c_n * Q_{n-1} + d_n * Q_{n-2}, Q_{-1}=0, Q_0=1
-- Prove: if P̃, Q̃ use -c_n instead of c_n, then
--   P̃_n = (-1)^{n+1} P_n and Q̃_n = (-1)^n Q_n
-- Hence P̃_n/Q̃_n = -P_n/Q_n.

-- Step 3: If Cohen's CF converges to π, then:
-- First tail T = α(1) + β(1)/(α(2) + ...) satisfies π = 3 + 6/T
-- So T = 6/(π-3)
-- Our CF = sign-flipped(T) = -T = -6/(π-3) = 6/(3-π)

-- The convergence of Cohen's CF itself needs:
-- Either cite it as a classical result, or prove via the Bauer-Muir chain.
-- For formalization, the cleanest route is the remainder-certificate pattern:
-- Define Cohen's convergents, show the remainder q_n·π - p_n decays geometrically.
```

## Key Mathlib APIs
- `Real.pi_gt_three : 3 < π` — ensures 3-π ≠ 0
- `Real.arctan_one : arctan 1 = π / 4`
- `Filter.Tendsto` for limits
- `Nat.rec` / structural induction for convergent sequences

## What needs to be built
1. General CF convergent machinery (could be shared with P2.3)
2. Sign-flip lemma (pure algebra, induction + ring)
3. Cohen's CF convergence to π (the main analytical content)
4. Final algebra: -T = 6/(3-π)

## Hard rules
- No sorry, no axiom, no native_decide
- Use `ring` for polynomial identities
- Use `norm_num` for numerical checks
- If Cohen's CF convergence is too hard to prove from scratch in one dispatch,
  accept a formulation where Cohen's convergence is the target theorem
  (prove everything else, state Cohen's convergence, and close the main theorem from it).

## Verification
```bash
export PATH="$HOME/.elan/bin:$PATH"
lake env lean RamanujanChallenge/Problem21.lean 2>&1
```

## Stall protocol
If stuck on proving Cohen's CF converges to π, deliver the file with
the sign-flip infrastructure complete + Cohen convergence as the one remaining sorry.
Report what approach you tried and the exact goal state.
