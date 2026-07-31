/-
  Ramanujan Challenge Problem 2.2: Euler's constant γ as an Apéry limit.

  This file gives a faithful formal model of the two sequences in the
  challenge.  It proves the index shift, defines the order-three recurrence,
  verifies both triples of initial values, proves the recurrence and its
  uniqueness, and records γ = 1 - Γ'(2).  It also identifies the challenge
  sequences exactly as a first-order Ore transform of the factorial-scaled
  numerator and denominator recurrence in Rivoal's construction.

  An earlier draft identified these sequences with Aptekarev's published
  approximants.  That identification is false: Aptekarev--Tulyakov's
  four-term recurrence starts with (0,2,31) and (1,3,50) and has different
  coefficients.  Positivity of Rivoal's denominator is proved below directly
  from its recurrence.  The remaining analytic input from Rivoal (convergence
  of his ratio) is isolated explicitly rather than being silently treated as
  a Lean theorem.
-/
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Data.Nat.Factorial.Basic
import Mathlib.NumberTheory.Harmonic.GammaDeriv

noncomputable section

open Filter Topology Real
open scoped BigOperators

namespace RamanujanChallenge.P22

/-! ## The shifted challenge recurrence

The challenge recurrence (after shift m = n + 3) has coefficients:
  c̃₀(m) = -8m³ + 21m² - 15m + 4
  c̃₁(m) = 24m⁵ - 23m⁴ - 51m³ + 39m² - m - 12
  c̃₂(m) = -m(m-1)(24m⁵ - 87m⁴ + 34m³ + 66m² - 3m - 20)
  c̃₃(m) = m(m-2)(m-1)⁴(8m³ + 3m² - 3m - 2)

Leading asymptotics: -8m³, 24m⁵, -24m⁷, 8m⁹ (factorial gauge k = 2,
characteristic polynomial -8(r-1)³ — maximally resonant).
-/

def c0_shifted (m : ℤ) : ℤ := -8 * m ^ 3 + 21 * m ^ 2 - 15 * m + 4

def c1_shifted (m : ℤ) : ℤ := 24 * m ^ 5 - 23 * m ^ 4 - 51 * m ^ 3 + 39 * m ^ 2 - m - 12

def c2_shifted (m : ℤ) : ℤ :=
  -m * (m - 1) * (24 * m ^ 5 - 87 * m ^ 4 + 34 * m ^ 3 + 66 * m ^ 2 - 3 * m - 20)

def c3_shifted (m : ℤ) : ℤ :=
  m * (m - 2) * (m - 1) ^ 4 * (8 * m ^ 3 + 3 * m ^ 2 - 3 * m - 2)

/-! ## The original challenge coefficients

The original coefficients c₀(n), c₁(n), c₂(n), c₃(n) are related by
the shift n = m - 3, i.e., c_j(n) = c̃_j(n + 3).
-/

def c0_orig (n : ℤ) : ℤ := -8 * n ^ 3 - 51 * n ^ 2 - 105 * n - 68

def c1_orig (n : ℤ) : ℤ := 24 * n ^ 5 + 337 * n ^ 4 + 1833 * n ^ 3 +
  4818 * n ^ 2 + 6092 * n + 2928

def c2_orig (n : ℤ) : ℤ :=
  -(n + 2) * (n + 3) * (24 * n ^ 5 + 273 * n ^ 4 + 1150 * n ^ 3 +
    2154 * n ^ 2 + 1635 * n + 268)

def c3_orig (n : ℤ) : ℤ :=
  (n + 1) * (n + 2) ^ 4 * (n + 3) * (8 * n ^ 3 + 75 * n ^ 2 + 231 * n + 232)

/-! ## Verification: shifted coefficients = original at n = m - 3

This is a finite polynomial identity, verified by comparing coefficients.
-/

theorem shift_c0 (m : ℤ) : c0_shifted m = c0_orig (m - 3) := by
  simp only [c0_shifted, c0_orig]; ring

theorem shift_c1 (m : ℤ) : c1_shifted m = c1_orig (m - 3) := by
  simp only [c1_shifted, c1_orig]; ring

theorem shift_c2 (m : ℤ) : c2_shifted m = c2_orig (m - 3) := by
  simp only [c2_shifted, c2_orig]; ring

theorem shift_c3 (m : ℤ) : c3_shifted m = c3_orig (m - 3) := by
  simp only [c3_shifted, c3_orig]; ring

/-! ## Initial values

(P₀, P₁, P₂) = (0, 7, 179)
(Q₀, Q₁, Q₂) = (1, 12, 306)

These are the challenge's `(p₋₃,p₋₂,p₋₁)` and `(q₋₃,q₋₂,q₋₁)` after shifting
the index by three.
-/

def challengeP22Init : Fin 3 → ℤ := ![0, 7, 179]
def challengeQ22Init : Fin 3 → ℤ := ![1, 12, 306]

/-! ## The gauge structure

The degree pattern (3, 5, 7, 9) with factorial gauge k = 2 gives
characteristic polynomial -8(r-1)³ (triple root at r = 1).

This algebraic identity alone does not establish asymptotics in the resonant
case, so no asymptotic conclusion is drawn from it here.
-/

theorem gauge_degree : (3 : ℕ) = 5 - 2 ∧ (5 : ℕ) = 7 - 2 ∧ (7 : ℕ) = 9 - 2 := by
  omega

/-! ## The actual challenge sequences -/

def c0Q22 (n : ℕ) : ℚ :=
  -8 * (n : ℚ) ^ 3 - 51 * (n : ℚ) ^ 2 - 105 * (n : ℚ) - 68

def c1Q22 (n : ℕ) : ℚ :=
  24 * (n : ℚ) ^ 5 + 337 * (n : ℚ) ^ 4 + 1833 * (n : ℚ) ^ 3 +
    4818 * (n : ℚ) ^ 2 + 6092 * (n : ℚ) + 2928

def c2Q22 (n : ℕ) : ℚ :=
  -((n : ℚ) + 2) * ((n : ℚ) + 3) *
    (24 * (n : ℚ) ^ 5 + 273 * (n : ℚ) ^ 4 + 1150 * (n : ℚ) ^ 3 +
      2154 * (n : ℚ) ^ 2 + 1635 * (n : ℚ) + 268)

def c3Q22 (n : ℕ) : ℚ :=
  ((n : ℚ) + 1) * ((n : ℚ) + 2) ^ 4 * ((n : ℚ) + 3) *
    (8 * (n : ℚ) ^ 3 + 75 * (n : ℚ) ^ 2 + 231 * (n : ℚ) + 232)

/-- The challenge recurrence, with the three negative initial indices shifted
to indices `0,1,2`. -/
def SatisfiesRecurrence22 (u : ℕ → ℚ) : Prop :=
  ∀ n : ℕ,
    c0Q22 n * u (n + 3) + c1Q22 n * u (n + 2) +
      c2Q22 n * u (n + 1) + c3Q22 n * u n = 0

theorem c0Q22_neg (n : ℕ) : c0Q22 n < 0 := by
  simp only [c0Q22]
  have hn : (0 : ℚ) ≤ n := by positivity
  have hn₂ : (0 : ℚ) ≤ (n : ℚ) ^ 2 := sq_nonneg _
  have hn₃ : (0 : ℚ) ≤ (n : ℚ) ^ 3 := by positivity
  nlinarith

theorem c0Q22_ne_zero (n : ℕ) : c0Q22 n ≠ 0 := (c0Q22_neg n).ne

/-- The unique solution of the recurrence with a prescribed initial triple. -/
def recurrenceSolution22 (u₀ u₁ u₂ : ℚ) : ℕ → ℚ
  | 0 => u₀
  | 1 => u₁
  | 2 => u₂
  | n + 3 =>
      -(c1Q22 n * recurrenceSolution22 u₀ u₁ u₂ (n + 2) +
          c2Q22 n * recurrenceSolution22 u₀ u₁ u₂ (n + 1) +
          c3Q22 n * recurrenceSolution22 u₀ u₁ u₂ n) / c0Q22 n

def challengeP22 : ℕ → ℚ := recurrenceSolution22 0 7 179

def challengeQ22 : ℕ → ℚ := recurrenceSolution22 1 12 306

@[simp] theorem challengeP22_zero : challengeP22 0 = 0 := rfl
@[simp] theorem challengeP22_one : challengeP22 1 = 7 := rfl
@[simp] theorem challengeP22_two : challengeP22 2 = 179 := rfl
@[simp] theorem challengeQ22_zero : challengeQ22 0 = 1 := rfl
@[simp] theorem challengeQ22_one : challengeQ22 1 = 12 := rfl
@[simp] theorem challengeQ22_two : challengeQ22 2 = 306 := rfl

theorem recurrenceSolution22_satisfies (u₀ u₁ u₂ : ℚ) :
    SatisfiesRecurrence22 (recurrenceSolution22 u₀ u₁ u₂) := by
  intro n
  simp only [recurrenceSolution22]
  field_simp [c0Q22_ne_zero]
  ring

theorem challengeP22_recurrence : SatisfiesRecurrence22 challengeP22 :=
  recurrenceSolution22_satisfies 0 7 179

theorem challengeQ22_recurrence : SatisfiesRecurrence22 challengeQ22 :=
  recurrenceSolution22_satisfies 1 12 306

/-- A transcription check beyond the three initial entries. -/
theorem challengeP22_three : challengeP22 3 = 7542 := by
  norm_num [challengeP22, recurrenceSolution22, c0Q22, c1Q22, c2Q22, c3Q22]

/-- A transcription check beyond the three initial entries. -/
theorem challengeQ22_three : challengeQ22 3 = 13056 := by
  norm_num [challengeQ22, recurrenceSolution22, c0Q22, c1Q22, c2Q22, c3Q22]

/-- Three initial values determine a solution of the challenge recurrence. -/
theorem recurrence22_unique
    (u v : ℕ → ℚ)
    (h₀ : u 0 = v 0) (h₁ : u 1 = v 1) (h₂ : u 2 = v 2)
    (hu : SatisfiesRecurrence22 u) (hv : SatisfiesRecurrence22 v) :
    ∀ n : ℕ, u n = v n := by
  intro n
  induction n using Nat.strong_induction_on with
  | h n ih =>
      by_cases hn₀ : n = 0
      · simpa [hn₀] using h₀
      by_cases hn₁ : n = 1
      · simpa [hn₁] using h₁
      by_cases hn₂ : n = 2
      · simpa [hn₂] using h₂
      have hn : 3 ≤ n := by omega
      let k := n - 3
      have hkn : k + 3 = n := by
        dsimp [k]
        omega
      have hk₀ : k < n := by omega
      have hk₁ : k + 1 < n := by omega
      have hk₂ : k + 2 < n := by omega
      have hru := hu k
      have hrv := hv k
      rw [ih k hk₀, ih (k + 1) hk₁, ih (k + 2) hk₂] at hru
      have hmul :
          c0Q22 k * (u (k + 3) - v (k + 3)) = 0 := by
        rw [mul_sub]
        linarith
      have heq : u (k + 3) = v (k + 3) :=
        sub_eq_zero.mp ((mul_eq_zero.mp hmul).resolve_left (c0Q22_ne_zero k))
      simpa [hkn] using heq

/-! ## Rivoal's companion recurrence and the Ore intertwiner

Rivoal's rational approximants to Euler's constant satisfy

```
(n+3)²(8n+11)(8n+19) R_{n+3}
  - (n+3)(8n+11)(24n²+145n+215) R_{n+2}
  + (8n+27)(24n³+105n²+124n+25) R_{n+1}
  - (n+2)²(8n+19)(8n+27) R_n = 0.
```

This is the `x = 1` specialization in T. Rivoal, *Rational
approximations for values of derivatives of the Gamma function*,
Trans. Amer. Math. Soc. 361 (2009), Theorem 1 and Corollary 1.  The
initial triples in that normalization are
`(-1, 4, 77/4)` and `(1, 7, 65/2)`, and the cited limit has direction
`numerator / denominator → γ`.

It is more convenient here to use the factorial gauge
`T_n = (n!)² R_n`.  After cancelling a common factor `(n+3)²`, the
recurrence has the coefficients below.  This section is purely algebraic:
it proves that the challenge recurrence is a first-order Ore transform of
this scaled Rivoal recurrence.  The analytic convergence theorem for
Rivoal's approximants is kept separate below.
-/

def rivoalC0Q22 (n : ℕ) : ℚ :=
  (8 * (n : ℚ) + 11) * (8 * (n : ℚ) + 19)

def rivoalC1Q22 (n : ℕ) : ℚ :=
  -((n : ℚ) + 3) * (8 * (n : ℚ) + 11) *
    (24 * (n : ℚ) ^ 2 + 145 * (n : ℚ) + 215)

def rivoalC2Q22 (n : ℕ) : ℚ :=
  ((n : ℚ) + 2) ^ 2 * (8 * (n : ℚ) + 27) *
    (24 * (n : ℚ) ^ 3 + 105 * (n : ℚ) ^ 2 + 124 * (n : ℚ) + 25)

def rivoalC3Q22 (n : ℕ) : ℚ :=
  -((n : ℚ) + 1) ^ 2 * ((n : ℚ) + 2) ^ 4 *
    (8 * (n : ℚ) + 19) * (8 * (n : ℚ) + 27)

def SatisfiesRivoalScaled22 (u : ℕ → ℚ) : Prop :=
  ∀ n : ℕ,
    rivoalC0Q22 n * u (n + 3) + rivoalC1Q22 n * u (n + 2) +
      rivoalC2Q22 n * u (n + 1) + rivoalC3Q22 n * u n = 0

theorem rivoalC0Q22_pos (n : ℕ) : 0 < rivoalC0Q22 n := by
  simp only [rivoalC0Q22]
  positivity

theorem rivoalC0Q22_ne_zero (n : ℕ) : rivoalC0Q22 n ≠ 0 :=
  (rivoalC0Q22_pos n).ne'

/-- The unique solution of the factorial-scaled Rivoal recurrence with a
prescribed initial triple. -/
def rivoalScaledSolution22 (u₀ u₁ u₂ : ℚ) : ℕ → ℚ
  | 0 => u₀
  | 1 => u₁
  | 2 => u₂
  | n + 3 =>
      -(rivoalC1Q22 n * rivoalScaledSolution22 u₀ u₁ u₂ (n + 2) +
          rivoalC2Q22 n * rivoalScaledSolution22 u₀ u₁ u₂ (n + 1) +
          rivoalC3Q22 n * rivoalScaledSolution22 u₀ u₁ u₂ n) / rivoalC0Q22 n

/-- `(n!)²` times Rivoal's numerator sequence
`(-1, 4, 77/4, ...)`. -/
def rivoalScaledP22 : ℕ → ℚ := rivoalScaledSolution22 (-1) 4 77

/-- `(n!)²` times Rivoal's denominator sequence
`(1, 7, 65/2, ...)`. -/
def rivoalScaledQ22 : ℕ → ℚ := rivoalScaledSolution22 1 7 130

theorem rivoalScaledSolution22_satisfies (u₀ u₁ u₂ : ℚ) :
    SatisfiesRivoalScaled22 (rivoalScaledSolution22 u₀ u₁ u₂) := by
  intro n
  simp only [rivoalScaledSolution22]
  field_simp [rivoalC0Q22_ne_zero]
  ring

theorem rivoalScaledP22_recurrence :
    SatisfiesRivoalScaled22 rivoalScaledP22 :=
  rivoalScaledSolution22_satisfies (-1) 4 77

theorem rivoalScaledQ22_recurrence :
    SatisfiesRivoalScaled22 rivoalScaledQ22 :=
  rivoalScaledSolution22_satisfies 1 7 130

@[simp] theorem rivoalScaledP22_zero : rivoalScaledP22 0 = -1 := rfl
@[simp] theorem rivoalScaledP22_one : rivoalScaledP22 1 = 4 := rfl
@[simp] theorem rivoalScaledP22_two : rivoalScaledP22 2 = 77 := rfl
@[simp] theorem rivoalScaledQ22_zero : rivoalScaledQ22 0 = 1 := rfl
@[simp] theorem rivoalScaledQ22_one : rivoalScaledQ22 1 = 7 := rfl
@[simp] theorem rivoalScaledQ22_two : rivoalScaledQ22 2 = 130 := rfl

theorem rivoalScaledP22_three : rivoalScaledP22 3 = 2523 := by
  norm_num [rivoalScaledP22, rivoalScaledSolution22, rivoalC0Q22,
    rivoalC1Q22, rivoalC2Q22, rivoalC3Q22]

theorem rivoalScaledQ22_three : rivoalScaledQ22 3 = 4362 := by
  norm_num [rivoalScaledQ22, rivoalScaledSolution22, rivoalC0Q22,
    rivoalC1Q22, rivoalC2Q22, rivoalC3Q22]

/-! ## Positivity of Rivoal's denominator

Undoing the factorial gauge exposes a simple invariant cone.  If `Q_n` denotes
the unscaled denominator, then positivity of

```
Q_n,   Q_{n+1} - Q_n,   Q_{n+2} - 2 Q_{n+1} + Q_n
```

is preserved by the recurrence.  This proves positivity directly from the
fixed recurrence and its initial values; no asymptotic theorem is used.
-/

/-- Rivoal's denominator before multiplication by `(n!)²`. -/
def rivoalQ22 (n : ℕ) : ℚ :=
  rivoalScaledQ22 n / ((n.factorial : ℚ) ^ 2)

/-- Rivoal's numerator before multiplication by `(n!)²`. -/
def rivoalP22 (n : ℕ) : ℚ :=
  rivoalScaledP22 n / ((n.factorial : ℚ) ^ 2)

def rivoalUnscaledC0Q22 (n : ℕ) : ℚ :=
  ((n : ℚ) + 3) ^ 2 * (8 * (n : ℚ) + 11) * (8 * (n : ℚ) + 19)

def rivoalUnscaledC1Q22 (n : ℕ) : ℚ :=
  -((n : ℚ) + 3) * (8 * (n : ℚ) + 11) *
    (24 * (n : ℚ) ^ 2 + 145 * (n : ℚ) + 215)

def rivoalUnscaledC2Q22 (n : ℕ) : ℚ :=
  (8 * (n : ℚ) + 27) *
    (24 * (n : ℚ) ^ 3 + 105 * (n : ℚ) ^ 2 + 124 * (n : ℚ) + 25)

def rivoalUnscaledC3Q22 (n : ℕ) : ℚ :=
  -((n : ℚ) + 2) ^ 2 * (8 * (n : ℚ) + 19) * (8 * (n : ℚ) + 27)

def SatisfiesRivoalUnscaled22 (u : ℕ → ℚ) : Prop :=
  ∀ n : ℕ,
    rivoalUnscaledC0Q22 n * u (n + 3) +
      rivoalUnscaledC1Q22 n * u (n + 2) +
      rivoalUnscaledC2Q22 n * u (n + 1) +
      rivoalUnscaledC3Q22 n * u n = 0

@[simp] theorem rivoalQ22_zero : rivoalQ22 0 = 1 := by
  norm_num [rivoalQ22]

@[simp] theorem rivoalQ22_one : rivoalQ22 1 = 7 := by
  norm_num [rivoalQ22]

@[simp] theorem rivoalQ22_two : rivoalQ22 2 = 65 / 2 := by
  norm_num [rivoalQ22]

@[simp] theorem rivoalP22_zero : rivoalP22 0 = -1 := by
  norm_num [rivoalP22]

@[simp] theorem rivoalP22_one : rivoalP22 1 = 4 := by
  norm_num [rivoalP22]

@[simp] theorem rivoalP22_two : rivoalP22 2 = 77 / 4 := by
  norm_num [rivoalP22]

theorem rivoalScaledQ22_eq_factorial_sq_mul (n : ℕ) :
    rivoalScaledQ22 n = ((n.factorial : ℚ) ^ 2) * rivoalQ22 n := by
  rw [rivoalQ22]
  field_simp

theorem rivoalScaledP22_eq_factorial_sq_mul (n : ℕ) :
    rivoalScaledP22 n = ((n.factorial : ℚ) ^ 2) * rivoalP22 n := by
  rw [rivoalP22]
  field_simp

theorem rivoalQ22_recurrence : SatisfiesRivoalUnscaled22 rivoalQ22 := by
  intro n
  have h := rivoalScaledQ22_recurrence n
  rw [rivoalScaledQ22_eq_factorial_sq_mul (n + 3),
    rivoalScaledQ22_eq_factorial_sq_mul (n + 2),
    rivoalScaledQ22_eq_factorial_sq_mul (n + 1),
    rivoalScaledQ22_eq_factorial_sq_mul n] at h
  simp only [rivoalUnscaledC0Q22, rivoalUnscaledC1Q22,
    rivoalUnscaledC2Q22, rivoalUnscaledC3Q22]
  simp only [rivoalC0Q22, rivoalC1Q22, rivoalC2Q22, rivoalC3Q22] at h
  norm_num [Nat.factorial_succ, Nat.cast_add, Nat.cast_one,
    Nat.add_assoc] at h
  let K : ℚ :=
    (n.factorial : ℚ) ^ 2 * ((n : ℚ) + 1) ^ 2 * ((n : ℚ) + 2) ^ 2
  have hK : K ≠ 0 := by
    dsimp [K]
    positivity
  apply mul_left_cancel₀ hK
  dsimp [K]
  linear_combination h

theorem rivoalP22_recurrence : SatisfiesRivoalUnscaled22 rivoalP22 := by
  intro n
  have h := rivoalScaledP22_recurrence n
  rw [rivoalScaledP22_eq_factorial_sq_mul (n + 3),
    rivoalScaledP22_eq_factorial_sq_mul (n + 2),
    rivoalScaledP22_eq_factorial_sq_mul (n + 1),
    rivoalScaledP22_eq_factorial_sq_mul n] at h
  simp only [rivoalUnscaledC0Q22, rivoalUnscaledC1Q22,
    rivoalUnscaledC2Q22, rivoalUnscaledC3Q22]
  simp only [rivoalC0Q22, rivoalC1Q22, rivoalC2Q22, rivoalC3Q22] at h
  norm_num [Nat.factorial_succ, Nat.cast_add, Nat.cast_one,
    Nat.add_assoc] at h
  let K : ℚ :=
    (n.factorial : ℚ) ^ 2 * ((n : ℚ) + 1) ^ 2 * ((n : ℚ) + 2) ^ 2
  have hK : K ≠ 0 := by
    dsimp [K]
    positivity
  apply mul_left_cancel₀ hK
  dsimp [K]
  linear_combination h

def rivoalConeA22 (n : ℕ) : ℚ :=
  512 * (n : ℚ) ^ 3 + 3776 * (n : ℚ) ^ 2 + 8872 * (n : ℚ) + 6591

def rivoalConeB22 (n : ℕ) : ℚ :=
  2 * (320 * (n : ℚ) ^ 3 + 2336 * (n : ℚ) ^ 2 +
    5415 * (n : ℚ) + 3936)

def rivoalConeC22 (n : ℕ) : ℚ :=
  ((n : ℚ) + 3) * (8 * (n : ℚ) + 11) *
    (8 * (n : ℚ) ^ 2 + 59 * (n : ℚ) + 101)

theorem rivoalUnscaledC0Q22_pos (n : ℕ) :
    0 < rivoalUnscaledC0Q22 n := by
  simp only [rivoalUnscaledC0Q22]
  positivity

theorem rivoalConeA22_pos (n : ℕ) : 0 < rivoalConeA22 n := by
  simp only [rivoalConeA22]
  positivity

theorem rivoalConeB22_pos (n : ℕ) : 0 < rivoalConeB22 n := by
  simp only [rivoalConeB22]
  positivity

theorem rivoalConeC22_pos (n : ℕ) : 0 < rivoalConeC22 n := by
  simp only [rivoalConeC22]
  positivity

/-- The recurrence preserves positivity of the second forward difference:
after clearing its positive leading coefficient, the new difference is a
positive linear combination of the preceding value and two differences. -/
theorem rivoalSecondDifference22
    (u : ℕ → ℚ) (hu : SatisfiesRivoalUnscaled22 u) (n : ℕ) :
    rivoalUnscaledC0Q22 n *
        (u (n + 3) - 2 * u (n + 2) + u (n + 1)) =
      rivoalConeA22 n * u n +
        rivoalConeB22 n * (u (n + 1) - u n) +
        rivoalConeC22 n * (u (n + 2) - 2 * u (n + 1) + u n) := by
  have h := hu n
  simp only [rivoalUnscaledC0Q22, rivoalUnscaledC1Q22,
    rivoalUnscaledC2Q22, rivoalUnscaledC3Q22, rivoalConeA22,
    rivoalConeB22, rivoalConeC22] at h ⊢
  linear_combination h

def RivoalPositiveCone22 (n : ℕ) : Prop :=
  0 < rivoalQ22 n ∧
    0 < rivoalQ22 (n + 1) - rivoalQ22 n ∧
    0 < rivoalQ22 (n + 2) - 2 * rivoalQ22 (n + 1) + rivoalQ22 n

theorem rivoalPositiveCone22_zero : RivoalPositiveCone22 0 := by
  norm_num [RivoalPositiveCone22]

theorem rivoalPositiveCone22_step (n : ℕ)
    (h : RivoalPositiveCone22 n) :
    RivoalPositiveCone22 (n + 1) := by
  rcases h with ⟨hq, hd, he⟩
  have hrec := rivoalSecondDifference22
    rivoalQ22 rivoalQ22_recurrence n
  have hRhs :
      0 <
        rivoalConeA22 n * rivoalQ22 n +
          rivoalConeB22 n * (rivoalQ22 (n + 1) - rivoalQ22 n) +
          rivoalConeC22 n *
            (rivoalQ22 (n + 2) - 2 * rivoalQ22 (n + 1) +
              rivoalQ22 n) := by
    exact add_pos
      (add_pos (mul_pos (rivoalConeA22_pos n) hq)
        (mul_pos (rivoalConeB22_pos n) hd))
      (mul_pos (rivoalConeC22_pos n) he)
  have hProd :
      0 < rivoalUnscaledC0Q22 n *
        (rivoalQ22 (n + 3) - 2 * rivoalQ22 (n + 2) +
          rivoalQ22 (n + 1)) := by
    rw [hrec]
    exact hRhs
  have hnew :
      0 < rivoalQ22 (n + 3) - 2 * rivoalQ22 (n + 2) +
        rivoalQ22 (n + 1) :=
    pos_of_mul_pos_right hProd (rivoalUnscaledC0Q22_pos n).le
  constructor
  · linarith
  constructor
  · linarith
  · simpa [Nat.add_assoc] using hnew

theorem rivoalPositiveCone22 (n : ℕ) : RivoalPositiveCone22 n := by
  induction n with
  | zero => exact rivoalPositiveCone22_zero
  | succ n ih =>
      simpa [Nat.succ_eq_add_one] using rivoalPositiveCone22_step n ih

theorem rivoalQ22_pos (n : ℕ) : 0 < rivoalQ22 n :=
  (rivoalPositiveCone22 n).1

theorem rivoalScaledQ22_pos (n : ℕ) : 0 < rivoalScaledQ22 n := by
  rw [rivoalScaledQ22_eq_factorial_sq_mul]
  have hq := rivoalQ22_pos n
  positivity

/-! ## Explicit hypergeometric formulas

The unscaled Rivoal solutions have finite-sum formulas.  We use the
hypergeometric weight recursively in `k`; below it is identified with
`choose n k ^ 2 / k!`.
-/

def rivoalWeight22 (n : ℕ) : ℕ → ℚ
  | 0 => 1
  | k + 1 =>
      rivoalWeight22 n k * ((n : ℚ) - (k : ℚ)) ^ 2 /
        ((k : ℚ) + 1) ^ 3

def rivoalHarmonicKernel22 (n k : ℕ) : ℚ :=
  3 * harmonic k - 2 * harmonic (n - k)

def rivoalExplicitQTerm22 (n k : ℕ) : ℚ :=
  (2 * (n : ℚ) + (k : ℚ) + 1) * rivoalWeight22 n k

def rivoalExplicitPTerm22 (n k : ℕ) : ℚ :=
  rivoalWeight22 n k *
    ((2 * (n : ℚ) + (k : ℚ) + 1) * rivoalHarmonicKernel22 n k - 1)

def rivoalExplicitQ22 (n : ℕ) : ℚ :=
  ∑ k ∈ Finset.range (n + 1), rivoalExplicitQTerm22 n k

def rivoalExplicitP22 (n : ℕ) : ℚ :=
  ∑ k ∈ Finset.range (n + 1), rivoalExplicitPTerm22 n k

def rivoalApproxValue22 (n k : ℕ) : ℚ :=
  rivoalHarmonicKernel22 n k -
    1 / (2 * (n : ℚ) + (k : ℚ) + 1)

theorem rivoalWeight22_succ (n k : ℕ) :
    rivoalWeight22 n (k + 1) =
      rivoalWeight22 n k * ((n : ℚ) - (k : ℚ)) ^ 2 /
        ((k : ℚ) + 1) ^ 3 :=
  rfl

@[simp] theorem rivoalExplicitQ22_zero : rivoalExplicitQ22 0 = 1 := by
  norm_num [rivoalExplicitQ22, rivoalExplicitQTerm22, rivoalWeight22]

@[simp] theorem rivoalExplicitQ22_one : rivoalExplicitQ22 1 = 7 := by
  norm_num [rivoalExplicitQ22, rivoalExplicitQTerm22, rivoalWeight22,
    Finset.sum_range_succ]

@[simp] theorem rivoalExplicitQ22_two : rivoalExplicitQ22 2 = 65 / 2 := by
  norm_num [rivoalExplicitQ22, rivoalExplicitQTerm22, rivoalWeight22,
    Finset.sum_range_succ]

@[simp] theorem rivoalExplicitP22_zero : rivoalExplicitP22 0 = -1 := by
  norm_num [rivoalExplicitP22, rivoalExplicitPTerm22, rivoalWeight22,
    rivoalHarmonicKernel22, harmonic]

@[simp] theorem rivoalExplicitP22_one : rivoalExplicitP22 1 = 4 := by
  norm_num [rivoalExplicitP22, rivoalExplicitPTerm22, rivoalWeight22,
    rivoalHarmonicKernel22, harmonic, Finset.sum_range_succ]

@[simp] theorem rivoalExplicitP22_two : rivoalExplicitP22 2 = 77 / 4 := by
  norm_num [rivoalExplicitP22, rivoalExplicitPTerm22, rivoalWeight22,
    rivoalHarmonicKernel22, harmonic, Finset.sum_range_succ]

theorem rivoalWeight22_eq_zero_of_lt (n k : ℕ) (h : n < k) :
    rivoalWeight22 n k = 0 := by
  induction k with
  | zero => omega
  | succ k ih =>
      simp only [rivoalWeight22]
      by_cases hkn : k = n
      · subst k
        norm_num
      · rw [ih (by omega)]
        simp

theorem rivoalExplicitQTerm22_eq_zero_of_lt (n k : ℕ) (h : n < k) :
    rivoalExplicitQTerm22 n k = 0 := by
  simp [rivoalExplicitQTerm22, rivoalWeight22_eq_zero_of_lt n k h]

theorem rivoalExplicitPTerm22_eq_zero_of_lt (n k : ℕ) (h : n < k) :
    rivoalExplicitPTerm22 n k = 0 := by
  simp [rivoalExplicitPTerm22, rivoalWeight22_eq_zero_of_lt n k h]

theorem rivoalWeight22_shift (n k : ℕ) :
    ((n : ℚ) + 1 - (k : ℚ)) ^ 2 * rivoalWeight22 (n + 1) k =
      ((n : ℚ) + 1) ^ 2 * rivoalWeight22 n k := by
  induction k with
  | zero =>
      simp [rivoalWeight22]
  | succ k ih =>
      simp only [rivoalWeight22]
      norm_num [Nat.cast_add, Nat.cast_one] at ih ⊢
      field_simp
      linear_combination ((n : ℚ) - (k : ℚ)) ^ 2 * ih

theorem rivoalWeight22_eq_ratio_succ (n k : ℕ) :
    rivoalWeight22 n k =
      (((n : ℚ) + 1 - (k : ℚ)) ^ 2 / ((n : ℚ) + 1) ^ 2) *
        rivoalWeight22 (n + 1) k := by
  have h := rivoalWeight22_shift n k
  field_simp
  nlinarith

theorem rivoalWeight22_eq_choose (n k : ℕ) :
    rivoalWeight22 n k =
      ((n.choose k : ℕ) : ℚ) ^ 2 / (k.factorial : ℚ) := by
  induction k with
  | zero =>
      simp [rivoalWeight22]
  | succ k ih =>
      by_cases hk : k ≤ n
      · have hcNat := Nat.choose_succ_right_eq n k
        have hc :
            ((n.choose (k + 1) : ℕ) : ℚ) =
              ((n.choose k : ℕ) : ℚ) *
                ((n : ℚ) - (k : ℚ)) / ((k : ℚ) + 1) := by
          apply (eq_div_iff (by positivity)).2
          exact_mod_cast hcNat
        rw [rivoalWeight22, ih, hc, Nat.factorial_succ]
        norm_num [Nat.cast_add, Nat.cast_one]
        field_simp
      · have hnk : n < k := by omega
        rw [rivoalWeight22, rivoalWeight22_eq_zero_of_lt n k hnk,
          Nat.choose_eq_zero_of_lt (by omega)]
        simp

def rivoalCertificateDen22 (n : ℕ) : ℚ :=
  ((n : ℚ) + 1) ^ 2 * ((n : ℚ) + 2) ^ 2 * ((n : ℚ) + 3) ^ 2

def rivoalQCertificatePoly22 (n k : ℕ) : ℚ :=
  let x : ℚ := n
  let y : ℚ := k
  (64 * x ^ 4 + 624 * x ^ 3 + 2241 * x ^ 2 + 3524 * x + 2052) * y ^ 4 +
    (-256 * x ^ 5 - 3136 * x ^ 4 - 15204 * x ^ 3 - 36506 * x ^ 2 -
      43448 * x - 20520) * y ^ 3 +
    (1536 * x ^ 5 + 16880 * x ^ 4 + 72882 * x ^ 3 + 155328 * x ^ 2 +
      164166 * x + 69093) * y ^ 2 +
    (1024 * x ^ 7 + 12864 * x ^ 6 + 64352 * x ^ 5 + 159865 * x ^ 4 +
      191178 * x ^ 3 + 58269 * x ^ 2 - 82610 * x - 59562) * y +
    (-1024 * x ^ 8 - 16512 * x ^ 7 - 113104 * x ^ 6 - 428290 * x ^ 5 -
      975748 * x ^ 4 - 1359516 * x ^ 3 - 1117711 * x ^ 2 -
      484508 * x - 80403)

def rivoalQCertificate22 (n k : ℕ) : ℚ :=
  rivoalWeight22 (n + 3) k * (k : ℚ) ^ 3 *
    rivoalQCertificatePoly22 n k / rivoalCertificateDen22 n

def rivoalPCertificatePoly22 (n k : ℕ) : ℚ :=
  let x : ℚ := n
  let y : ℚ := k
  (-7 * (x + 2) ^ 2 * (8 * x + 19) * (8 * x + 27)) * y ^ 4 +
    (12 * (x + 2) ^ 2 * (2 * x + 5) * (8 * x + 19) *
      (8 * x + 27)) * y ^ 3 +
    (-5 * (8 * x + 27) *
      (192 * x ^ 4 + 1462 * x ^ 3 + 4176 * x ^ 2 + 5322 * x + 2559)) *
      y ^ 2 +
    (-4 * (8 * x + 27) *
      (128 * x ^ 6 + 1176 * x ^ 5 + 4075 * x ^ 4 + 6230 * x ^ 3 +
        2871 * x ^ 2 - 2406 * x - 2206)) * y +
    3 * (x + 3) *
      (1024 * x ^ 7 + 13440 * x ^ 6 + 72784 * x ^ 5 + 209938 * x ^ 4 +
        345934 * x ^ 3 + 321714 * x ^ 2 + 152569 * x + 26801)

def rivoalPCorrectionCertificate22 (n k : ℕ) : ℚ :=
  rivoalWeight22 (n + 3) k * (k : ℚ) ^ 2 *
    rivoalPCertificatePoly22 n k / rivoalCertificateDen22 n

def rivoalPCertificate22 (n k : ℕ) : ℚ :=
  rivoalQCertificate22 n k * rivoalHarmonicKernel22 (n + 3) k +
    rivoalPCorrectionCertificate22 n k

theorem rivoalCertificateDen22_pos (n : ℕ) :
    0 < rivoalCertificateDen22 n := by
  simp only [rivoalCertificateDen22]
  positivity

theorem rivoalCertificateDen22_ne_zero (n : ℕ) :
    rivoalCertificateDen22 n ≠ 0 :=
  (rivoalCertificateDen22_pos n).ne'

theorem rivoalHarmonicKernel22_eq_succ
    (n k : ℕ) (h : k ≤ n) :
    rivoalHarmonicKernel22 n k =
      rivoalHarmonicKernel22 (n + 1) k +
        2 / ((n : ℚ) + 1 - (k : ℚ)) := by
  have hsub : n + 1 - k = (n - k) + 1 := by omega
  simp only [rivoalHarmonicKernel22, hsub, harmonic_succ]
  norm_num [Nat.cast_add, Nat.cast_one, Nat.cast_sub h]
  field_simp
  ring

theorem rivoalHarmonicKernel22_succ_right
    (n k : ℕ) (h : k < n) :
    rivoalHarmonicKernel22 n (k + 1) =
      rivoalHarmonicKernel22 n k +
        3 / ((k : ℚ) + 1) + 2 / ((n : ℚ) - (k : ℚ)) := by
  have hsub : n - k = (n - (k + 1)) + 1 := by omega
  simp only [rivoalHarmonicKernel22, harmonic_succ, hsub]
  norm_num [Nat.cast_add, Nat.cast_one, Nat.cast_sub (Nat.le_of_lt h),
    Nat.cast_sub (by omega : k + 1 ≤ n)]
  field_simp
  ring

set_option maxRecDepth 10000 in
theorem rivoalExplicitQTerm22_telescope (n k : ℕ) :
    rivoalUnscaledC0Q22 n * rivoalExplicitQTerm22 (n + 3) k +
        rivoalUnscaledC1Q22 n * rivoalExplicitQTerm22 (n + 2) k +
        rivoalUnscaledC2Q22 n * rivoalExplicitQTerm22 (n + 1) k +
        rivoalUnscaledC3Q22 n * rivoalExplicitQTerm22 n k =
      rivoalQCertificate22 n (k + 1) - rivoalQCertificate22 n k := by
  simp only [rivoalExplicitQTerm22]
  rw [rivoalWeight22_eq_ratio_succ n k,
    rivoalWeight22_eq_ratio_succ (n + 1) k,
    rivoalWeight22_eq_ratio_succ (n + 2) k]
  simp only [rivoalQCertificate22, rivoalQCertificatePoly22,
    rivoalUnscaledC0Q22, rivoalUnscaledC1Q22,
    rivoalUnscaledC2Q22, rivoalUnscaledC3Q22]
  simp only [rivoalWeight22]
  norm_num [Nat.cast_add, Nat.cast_one, Nat.add_assoc]
  simp only [rivoalCertificateDen22]
  field_simp [rivoalCertificateDen22_ne_zero]
  ring

theorem rivoalExplicitPTerm22_decompose (n k : ℕ) :
    rivoalExplicitPTerm22 n k =
      rivoalExplicitQTerm22 n k * rivoalHarmonicKernel22 n k -
        rivoalWeight22 n k := by
  simp only [rivoalExplicitPTerm22, rivoalExplicitQTerm22]
  ring

theorem rivoalExplicitPTerm22_eq_q_mul_approx (n k : ℕ) :
    rivoalExplicitPTerm22 n k =
      rivoalExplicitQTerm22 n k * rivoalApproxValue22 n k := by
  have hden : 2 * (n : ℚ) + (k : ℚ) + 1 ≠ 0 := by
    positivity
  simp only [rivoalExplicitPTerm22, rivoalExplicitQTerm22,
    rivoalApproxValue22]
  field_simp [hden]

theorem rivoalExplicitP22_eq_weighted (n : ℕ) :
    rivoalExplicitP22 n =
      ∑ k ∈ Finset.range (n + 1),
        rivoalExplicitQTerm22 n k * rivoalApproxValue22 n k := by
  simp only [rivoalExplicitP22]
  apply Finset.sum_congr rfl
  intro k _
  exact rivoalExplicitPTerm22_eq_q_mul_approx n k

def rivoalResidualA22 (n k : ℕ) : ℚ :=
  (n : ℚ) + 1 - (k : ℚ)

def rivoalResidualB22 (n k : ℕ) : ℚ :=
  (n : ℚ) + 2 - (k : ℚ)

def rivoalResidualC22 (n k : ℕ) : ℚ :=
  (n : ℚ) + 3 - (k : ℚ)

def rivoalResidualD22 (k : ℕ) : ℚ :=
  (k : ℚ) + 1

def rivoalPResidualLeftCoefficient22 (n k : ℕ) : ℚ :=
  let a := rivoalResidualA22 n k
  let b := rivoalResidualB22 n k
  let c := rivoalResidualC22 n k
  rivoalUnscaledC0Q22 n * (-1) +
    rivoalUnscaledC1Q22 n *
      (c ^ 2 / ((n : ℚ) + 3) ^ 2) *
        (2 * (2 * ((n : ℚ) + 2) + (k : ℚ) + 1) / c - 1) +
    rivoalUnscaledC2Q22 n *
      ((b ^ 2 / ((n : ℚ) + 2) ^ 2) *
        (c ^ 2 / ((n : ℚ) + 3) ^ 2)) *
        ((2 * ((n : ℚ) + 1) + (k : ℚ) + 1) *
            (2 / b + 2 / c) - 1) +
    rivoalUnscaledC3Q22 n *
      ((a ^ 2 / ((n : ℚ) + 1) ^ 2) *
        ((b ^ 2 / ((n : ℚ) + 2) ^ 2) *
          (c ^ 2 / ((n : ℚ) + 3) ^ 2))) *
        ((2 * (n : ℚ) + (k : ℚ) + 1) *
            (2 / a + 2 / b + 2 / c) - 1)

def rivoalPResidualRightCoefficient22 (n k : ℕ) : ℚ :=
  let c := rivoalResidualC22 n k
  let d := rivoalResidualD22 k
  (c ^ 2 / d ^ 3) * d ^ 3 *
        rivoalQCertificatePoly22 n (k + 1) /
        rivoalCertificateDen22 n *
      (3 / d + 2 / c) +
    (c ^ 2 / d ^ 3) * d ^ 2 *
        rivoalPCertificatePoly22 n (k + 1) /
        rivoalCertificateDen22 n -
    (k : ℚ) ^ 2 * rivoalPCertificatePoly22 n k /
      rivoalCertificateDen22 n

set_option maxRecDepth 10000 in
theorem rivoalPResidualCoefficient22
    (n k : ℕ) (h : k ≤ n) :
    rivoalPResidualLeftCoefficient22 n k =
      rivoalPResidualRightCoefficient22 n k := by
  have ha : rivoalResidualA22 n k ≠ 0 := by
    simp only [rivoalResidualA22]
    have hkQ : (k : ℚ) ≤ (n : ℚ) := by exact_mod_cast h
    linarith
  have hb : rivoalResidualB22 n k ≠ 0 := by
    simp only [rivoalResidualB22]
    have hkQ : (k : ℚ) ≤ (n : ℚ) := by exact_mod_cast h
    linarith
  have hc : rivoalResidualC22 n k ≠ 0 := by
    simp only [rivoalResidualC22]
    have hkQ : (k : ℚ) ≤ (n : ℚ) := by exact_mod_cast h
    linarith
  have hd : rivoalResidualD22 k ≠ 0 := by
    simp only [rivoalResidualD22]
    positivity
  simp only [rivoalPResidualLeftCoefficient22,
    rivoalPResidualRightCoefficient22]
  field_simp [rivoalCertificateDen22_ne_zero, ha, hb, hc, hd]
  simp only [rivoalResidualA22, rivoalResidualB22,
    rivoalResidualC22, rivoalResidualD22,
    rivoalQCertificatePoly22, rivoalPCertificatePoly22,
    rivoalUnscaledC0Q22, rivoalUnscaledC1Q22,
    rivoalUnscaledC2Q22, rivoalUnscaledC3Q22,
    rivoalCertificateDen22]
  norm_num [Nat.cast_add, Nat.cast_one, Nat.add_assoc]
  ring

def rivoalPResidualLeft22 (n k : ℕ) : ℚ :=
  rivoalUnscaledC0Q22 n *
        (rivoalExplicitQTerm22 (n + 3) k *
            (rivoalHarmonicKernel22 (n + 3) k -
              rivoalHarmonicKernel22 (n + 3) k) -
          rivoalWeight22 (n + 3) k) +
    rivoalUnscaledC1Q22 n *
        (rivoalExplicitQTerm22 (n + 2) k *
            (rivoalHarmonicKernel22 (n + 2) k -
              rivoalHarmonicKernel22 (n + 3) k) -
          rivoalWeight22 (n + 2) k) +
    rivoalUnscaledC2Q22 n *
        (rivoalExplicitQTerm22 (n + 1) k *
            (rivoalHarmonicKernel22 (n + 1) k -
              rivoalHarmonicKernel22 (n + 3) k) -
          rivoalWeight22 (n + 1) k) +
    rivoalUnscaledC3Q22 n *
        (rivoalExplicitQTerm22 n k *
            (rivoalHarmonicKernel22 n k -
              rivoalHarmonicKernel22 (n + 3) k) -
          rivoalWeight22 n k)

def rivoalPResidualRight22 (n k : ℕ) : ℚ :=
  rivoalQCertificate22 n (k + 1) *
      (rivoalHarmonicKernel22 (n + 3) (k + 1) -
        rivoalHarmonicKernel22 (n + 3) k) +
    rivoalPCorrectionCertificate22 n (k + 1) -
    rivoalPCorrectionCertificate22 n k

set_option maxRecDepth 10000 in
theorem rivoalExplicitPResidual22_telescope_of_le
    (n k : ℕ) (h : k ≤ n) :
    rivoalPResidualLeft22 n k = rivoalPResidualRight22 n k := by
  have hc := rivoalPResidualCoefficient22 n k h
  simp only [rivoalPResidualLeft22, rivoalPResidualRight22]
  simp only [rivoalExplicitQTerm22, rivoalQCertificate22,
    rivoalPCorrectionCertificate22]
  rw [rivoalWeight22_eq_ratio_succ n k,
    rivoalWeight22_eq_ratio_succ (n + 1) k,
    rivoalWeight22_eq_ratio_succ (n + 2) k]
  rw [rivoalHarmonicKernel22_eq_succ n k h,
    rivoalHarmonicKernel22_eq_succ (n + 1) k (by omega),
    rivoalHarmonicKernel22_eq_succ (n + 2) k (by omega),
    rivoalHarmonicKernel22_succ_right (n + 3) k (by omega)]
  simp only [rivoalWeight22]
  norm_num [Nat.cast_add, Nat.cast_one, Nat.add_assoc]
  simp only [rivoalPResidualLeftCoefficient22,
    rivoalPResidualRightCoefficient22,
    rivoalResidualA22, rivoalResidualB22,
    rivoalResidualC22, rivoalResidualD22] at hc
  convert congrArg
      (fun z : ℚ => rivoalWeight22 (n + 3) k * z) hc using 1 <;>
    ring_nf (config := { red := .instances })

theorem rivoalExplicitPTerm22_telescope_of_residual
    (n k : ℕ)
    (hr : rivoalPResidualLeft22 n k =
      rivoalPResidualRight22 n k) :
    rivoalUnscaledC0Q22 n * rivoalExplicitPTerm22 (n + 3) k +
        rivoalUnscaledC1Q22 n * rivoalExplicitPTerm22 (n + 2) k +
        rivoalUnscaledC2Q22 n * rivoalExplicitPTerm22 (n + 1) k +
        rivoalUnscaledC3Q22 n * rivoalExplicitPTerm22 n k =
      rivoalPCertificate22 n (k + 1) - rivoalPCertificate22 n k := by
  rw [rivoalExplicitPTerm22_decompose,
    rivoalExplicitPTerm22_decompose,
    rivoalExplicitPTerm22_decompose,
    rivoalExplicitPTerm22_decompose]
  have hq := rivoalExplicitQTerm22_telescope n k
  simp only [rivoalPResidualLeft22, rivoalPResidualRight22] at hr
  simp only [rivoalPCertificate22]
  linear_combination
    rivoalHarmonicKernel22 (n + 3) k * hq + hr

theorem rivoalExplicitPTerm22_telescope_of_le
    (n k : ℕ) (h : k ≤ n) :
    rivoalUnscaledC0Q22 n * rivoalExplicitPTerm22 (n + 3) k +
        rivoalUnscaledC1Q22 n * rivoalExplicitPTerm22 (n + 2) k +
        rivoalUnscaledC2Q22 n * rivoalExplicitPTerm22 (n + 1) k +
        rivoalUnscaledC3Q22 n * rivoalExplicitPTerm22 n k =
      rivoalPCertificate22 n (k + 1) - rivoalPCertificate22 n k :=
  rivoalExplicitPTerm22_telescope_of_residual n k
    (rivoalExplicitPResidual22_telescope_of_le n k h)

set_option maxRecDepth 10000 in
theorem rivoalExplicitPResidual22_telescope_add_three (n : ℕ) :
    rivoalPResidualLeft22 n (n + 3) =
      rivoalPResidualRight22 n (n + 3) := by
  simp only [rivoalPResidualLeft22, rivoalPResidualRight22]
  rw [rivoalExplicitQTerm22_eq_zero_of_lt (n + 2) (n + 3) (by omega),
    rivoalExplicitQTerm22_eq_zero_of_lt (n + 1) (n + 3) (by omega),
    rivoalExplicitQTerm22_eq_zero_of_lt n (n + 3) (by omega),
    rivoalWeight22_eq_zero_of_lt (n + 2) (n + 3) (by omega),
    rivoalWeight22_eq_zero_of_lt (n + 1) (n + 3) (by omega),
    rivoalWeight22_eq_zero_of_lt n (n + 3) (by omega)]
  simp only [sub_self, mul_zero, zero_sub, rivoalQCertificate22,
    rivoalPCorrectionCertificate22]
  rw [rivoalWeight22_eq_zero_of_lt (n + 3) (n + 3 + 1) (by omega)]
  simp only [zero_mul, zero_div, zero_add]
  field_simp [rivoalCertificateDen22_ne_zero]
  simp only [rivoalPCertificatePoly22, rivoalUnscaledC0Q22,
    rivoalCertificateDen22]
  norm_num [Nat.cast_add, Nat.cast_one, Nat.add_assoc]
  ring

set_option maxRecDepth 10000 in
theorem rivoalExplicitPResidual22_telescope_add_two (n : ℕ) :
    rivoalPResidualLeft22 n (n + 2) =
      rivoalPResidualRight22 n (n + 2) := by
  simp only [rivoalPResidualLeft22, rivoalPResidualRight22]
  rw [rivoalExplicitQTerm22_eq_zero_of_lt (n + 1) (n + 2) (by omega),
    rivoalExplicitQTerm22_eq_zero_of_lt n (n + 2) (by omega),
    rivoalWeight22_eq_zero_of_lt (n + 1) (n + 2) (by omega),
    rivoalWeight22_eq_zero_of_lt n (n + 2) (by omega)]
  rw [rivoalHarmonicKernel22_eq_succ (n + 2) (n + 2) le_rfl,
    rivoalHarmonicKernel22_succ_right (n + 3) (n + 2) (by omega)]
  simp only [sub_self, mul_zero, zero_sub, rivoalExplicitQTerm22,
    rivoalQCertificate22, rivoalPCorrectionCertificate22]
  rw [rivoalWeight22_eq_ratio_succ (n + 2) (n + 2),
    rivoalWeight22_succ (n + 3) (n + 2)]
  norm_num [Nat.cast_add, Nat.cast_one, Nat.add_assoc]
  field_simp [rivoalCertificateDen22_ne_zero]
  simp only [rivoalQCertificatePoly22, rivoalPCertificatePoly22,
    rivoalUnscaledC0Q22, rivoalUnscaledC1Q22,
    rivoalCertificateDen22]
  norm_num [Nat.cast_add, Nat.cast_one, Nat.add_assoc]
  ring

set_option maxRecDepth 10000 in
theorem rivoalExplicitPResidual22_telescope_add_one (n : ℕ) :
    rivoalPResidualLeft22 n (n + 1) =
      rivoalPResidualRight22 n (n + 1) := by
  simp only [rivoalPResidualLeft22, rivoalPResidualRight22]
  rw [rivoalExplicitQTerm22_eq_zero_of_lt n (n + 1) (by omega),
    rivoalWeight22_eq_zero_of_lt n (n + 1) (by omega)]
  rw [rivoalHarmonicKernel22_eq_succ (n + 1) (n + 1) le_rfl,
    rivoalHarmonicKernel22_eq_succ (n + 2) (n + 1) (by omega),
    rivoalHarmonicKernel22_succ_right (n + 3) (n + 1) (by omega)]
  simp only [sub_self, mul_zero, zero_sub, rivoalExplicitQTerm22,
    rivoalQCertificate22, rivoalPCorrectionCertificate22]
  rw [rivoalWeight22_eq_ratio_succ (n + 1) (n + 1),
    rivoalWeight22_eq_ratio_succ (n + 2) (n + 1),
    rivoalWeight22_succ (n + 3) (n + 1)]
  norm_num [Nat.cast_add, Nat.cast_one, Nat.add_assoc]
  field_simp [rivoalCertificateDen22_ne_zero]
  simp only [rivoalQCertificatePoly22, rivoalPCertificatePoly22,
    rivoalUnscaledC0Q22, rivoalUnscaledC1Q22,
    rivoalUnscaledC2Q22, rivoalCertificateDen22]
  norm_num [Nat.cast_add, Nat.cast_one, Nat.add_assoc]
  ring

theorem rivoalExplicitPTerm22_telescope
    (n k : ℕ) (h : k < n + 4) :
    rivoalUnscaledC0Q22 n * rivoalExplicitPTerm22 (n + 3) k +
        rivoalUnscaledC1Q22 n * rivoalExplicitPTerm22 (n + 2) k +
        rivoalUnscaledC2Q22 n * rivoalExplicitPTerm22 (n + 1) k +
        rivoalUnscaledC3Q22 n * rivoalExplicitPTerm22 n k =
      rivoalPCertificate22 n (k + 1) - rivoalPCertificate22 n k := by
  by_cases hk : k ≤ n
  · exact rivoalExplicitPTerm22_telescope_of_le n k hk
  · have hk' : k = n + 1 ∨ k = n + 2 ∨ k = n + 3 := by
      omega
    rcases hk' with rfl | rfl | rfl
    · exact rivoalExplicitPTerm22_telescope_of_residual n (n + 1)
        (rivoalExplicitPResidual22_telescope_add_one n)
    · exact rivoalExplicitPTerm22_telescope_of_residual n (n + 2)
        (rivoalExplicitPResidual22_telescope_add_two n)
    · exact rivoalExplicitPTerm22_telescope_of_residual n (n + 3)
        (rivoalExplicitPResidual22_telescope_add_three n)

theorem rivoalExplicitQ22_sum_pad (n M : ℕ) (h : n + 1 ≤ M) :
    rivoalExplicitQ22 n =
      ∑ k ∈ Finset.range M, rivoalExplicitQTerm22 n k := by
  rw [rivoalExplicitQ22]
  apply Finset.sum_subset (Finset.range_mono h)
  intro k hkM hkn
  apply rivoalExplicitQTerm22_eq_zero_of_lt
  simp only [Finset.mem_range] at hkM hkn
  omega

theorem rivoalExplicitP22_sum_pad (n M : ℕ) (h : n + 1 ≤ M) :
    rivoalExplicitP22 n =
      ∑ k ∈ Finset.range M, rivoalExplicitPTerm22 n k := by
  rw [rivoalExplicitP22]
  apply Finset.sum_subset (Finset.range_mono h)
  intro k hkM hkn
  apply rivoalExplicitPTerm22_eq_zero_of_lt
  simp only [Finset.mem_range] at hkM hkn
  omega

theorem rivoalExplicitQ22_recurrence :
    SatisfiesRivoalUnscaled22 rivoalExplicitQ22 := by
  intro n
  rw [rivoalExplicitQ22_sum_pad (n + 3) (n + 4) (by omega),
    rivoalExplicitQ22_sum_pad (n + 2) (n + 4) (by omega),
    rivoalExplicitQ22_sum_pad (n + 1) (n + 4) (by omega),
    rivoalExplicitQ22_sum_pad n (n + 4) (by omega)]
  simp only [Finset.mul_sum]
  rw [← Finset.sum_add_distrib, ← Finset.sum_add_distrib,
    ← Finset.sum_add_distrib]
  calc
    (∑ k ∈ Finset.range (n + 4),
        (rivoalUnscaledC0Q22 n * rivoalExplicitQTerm22 (n + 3) k +
          rivoalUnscaledC1Q22 n * rivoalExplicitQTerm22 (n + 2) k +
          rivoalUnscaledC2Q22 n * rivoalExplicitQTerm22 (n + 1) k +
          rivoalUnscaledC3Q22 n * rivoalExplicitQTerm22 n k)) =
        ∑ k ∈ Finset.range (n + 4),
          (rivoalQCertificate22 n (k + 1) -
            rivoalQCertificate22 n k) := by
      apply Finset.sum_congr rfl
      intro k hk
      rw [← rivoalExplicitQTerm22_telescope]
    _ = rivoalQCertificate22 n (n + 4) -
        rivoalQCertificate22 n 0 := by
      rw [Finset.sum_range_sub]
    _ = 0 := by
      simp [rivoalQCertificate22,
        rivoalWeight22_eq_zero_of_lt (n + 3) (n + 4) (by omega)]

theorem rivoalExplicitP22_recurrence :
    SatisfiesRivoalUnscaled22 rivoalExplicitP22 := by
  intro n
  rw [rivoalExplicitP22_sum_pad (n + 3) (n + 4) (by omega),
    rivoalExplicitP22_sum_pad (n + 2) (n + 4) (by omega),
    rivoalExplicitP22_sum_pad (n + 1) (n + 4) (by omega),
    rivoalExplicitP22_sum_pad n (n + 4) (by omega)]
  simp only [Finset.mul_sum]
  rw [← Finset.sum_add_distrib, ← Finset.sum_add_distrib,
    ← Finset.sum_add_distrib]
  calc
    (∑ k ∈ Finset.range (n + 4),
        (rivoalUnscaledC0Q22 n * rivoalExplicitPTerm22 (n + 3) k +
          rivoalUnscaledC1Q22 n * rivoalExplicitPTerm22 (n + 2) k +
          rivoalUnscaledC2Q22 n * rivoalExplicitPTerm22 (n + 1) k +
          rivoalUnscaledC3Q22 n * rivoalExplicitPTerm22 n k)) =
        ∑ k ∈ Finset.range (n + 4),
          (rivoalPCertificate22 n (k + 1) -
            rivoalPCertificate22 n k) := by
      apply Finset.sum_congr rfl
      intro k hk
      rw [← rivoalExplicitPTerm22_telescope n k
        (Finset.mem_range.mp hk)]
    _ = rivoalPCertificate22 n (n + 4) -
        rivoalPCertificate22 n 0 := by
      rw [Finset.sum_range_sub]
    _ = 0 := by
      simp [rivoalPCertificate22, rivoalQCertificate22,
        rivoalPCorrectionCertificate22,
        rivoalWeight22_eq_zero_of_lt (n + 3) (n + 4) (by omega)]

/-- Three initial values determine a solution of the unscaled Rivoal
recurrence. -/
theorem rivoalUnscaled22_unique
    (u v : ℕ → ℚ)
    (h₀ : u 0 = v 0) (h₁ : u 1 = v 1) (h₂ : u 2 = v 2)
    (hu : SatisfiesRivoalUnscaled22 u)
    (hv : SatisfiesRivoalUnscaled22 v) :
    ∀ n : ℕ, u n = v n := by
  intro n
  induction n using Nat.strong_induction_on with
  | h n ih =>
      by_cases hn₀ : n = 0
      · simpa [hn₀] using h₀
      by_cases hn₁ : n = 1
      · simpa [hn₁] using h₁
      by_cases hn₂ : n = 2
      · simpa [hn₂] using h₂
      have hn : 3 ≤ n := by omega
      let k := n - 3
      have hkn : k + 3 = n := by
        dsimp [k]
        omega
      have hk₀ : k < n := by omega
      have hk₁ : k + 1 < n := by omega
      have hk₂ : k + 2 < n := by omega
      have hru := hu k
      have hrv := hv k
      rw [ih k hk₀, ih (k + 1) hk₁, ih (k + 2) hk₂] at hru
      have hmul :
          rivoalUnscaledC0Q22 k * (u (k + 3) - v (k + 3)) = 0 := by
        rw [mul_sub]
        linarith
      have heq : u (k + 3) = v (k + 3) :=
        sub_eq_zero.mp ((mul_eq_zero.mp hmul).resolve_left
          (rivoalUnscaledC0Q22_pos k).ne')
      simpa [hkn] using heq

theorem rivoalQ22_eq_explicit (n : ℕ) :
    rivoalQ22 n = rivoalExplicitQ22 n := by
  apply rivoalUnscaled22_unique
    rivoalQ22 rivoalExplicitQ22
    (by simp) (by simp) (by simp)
    rivoalQ22_recurrence rivoalExplicitQ22_recurrence n

theorem rivoalP22_eq_explicit (n : ℕ) :
    rivoalP22 n = rivoalExplicitP22 n := by
  apply rivoalUnscaled22_unique
    rivoalP22 rivoalExplicitP22
    (by simp) (by simp) (by simp)
    rivoalP22_recurrence rivoalExplicitP22_recurrence n

theorem rivoalQ22_sum_choose (n : ℕ) :
    rivoalQ22 n =
      ∑ k ∈ Finset.range (n + 1),
        (2 * (n : ℚ) + (k : ℚ) + 1) *
          ((n.choose k : ℕ) : ℚ) ^ 2 / (k.factorial : ℚ) := by
  rw [rivoalQ22_eq_explicit]
  simp only [rivoalExplicitQ22, rivoalExplicitQTerm22,
    rivoalWeight22_eq_choose]
  apply Finset.sum_congr rfl
  intro k _
  ring

theorem rivoalP22_sum_choose (n : ℕ) :
    rivoalP22 n =
      ∑ k ∈ Finset.range (n + 1),
        (((n.choose k : ℕ) : ℚ) ^ 2 / (k.factorial : ℚ)) *
          ((2 * (n : ℚ) + (k : ℚ) + 1) *
            (3 * harmonic k - 2 * harmonic (n - k)) - 1) := by
  rw [rivoalP22_eq_explicit]
  simp only [rivoalExplicitP22, rivoalExplicitPTerm22,
    rivoalHarmonicKernel22, rivoalWeight22_eq_choose]

theorem rivoalExplicitQ22_pos (n : ℕ) :
    0 < rivoalExplicitQ22 n := by
  rw [← rivoalQ22_eq_explicit]
  exact rivoalQ22_pos n

theorem rivoalExplicitQ22_ne_zero (n : ℕ) :
    rivoalExplicitQ22 n ≠ 0 :=
  (rivoalExplicitQ22_pos n).ne'

theorem rivoalScaledRatio22_eq_explicit (n : ℕ) :
    rivoalScaledP22 n / rivoalScaledQ22 n =
      rivoalExplicitP22 n / rivoalExplicitQ22 n := by
  rw [rivoalScaledP22_eq_factorial_sq_mul,
    rivoalScaledQ22_eq_factorial_sq_mul,
    rivoalP22_eq_explicit, rivoalQ22_eq_explicit]
  field_simp [rivoalExplicitQ22_ne_zero]

/-- The first-order Ore transform from the factorial-scaled Rivoal module to
the challenge module:

`C_n = (T_{n+1} + (n+1)(3n+4)T_n) / (8n+11)`.
-/
def rivoalOre22 (u : ℕ → ℚ) (n : ℕ) : ℚ :=
  (u (n + 1) + ((n : ℚ) + 1) * (3 * (n : ℚ) + 4) * u n) /
    (8 * (n : ℚ) + 11)

/- The displayed Ore operator intertwines Rivoal's scaled recurrence with
the challenge recurrence. -/
set_option maxRecDepth 10000 in
theorem rivoalOre22_satisfies
    (u : ℕ → ℚ) (hu : SatisfiesRivoalScaled22 u) :
    SatisfiesRecurrence22 (rivoalOre22 u) := by
  intro n
  have h₀ := hu n
  have h₁ := hu (n + 1)
  simp only [rivoalC0Q22, rivoalC1Q22, rivoalC2Q22, rivoalC3Q22] at h₀ h₁
  simp only [c0Q22, c1Q22, c2Q22, c3Q22, rivoalOre22]
  norm_num [Nat.cast_add, Nat.cast_one, Nat.add_assoc] at h₁ ⊢
  field_simp
  linear_combination
    -((n : ℚ) + 3) * (3 * (n : ℚ) + 4) * (8 * (n : ℚ) + 35) *
        (8 * (n : ℚ) ^ 3 + 75 * (n : ℚ) ^ 2 + 231 * (n : ℚ) + 232) * h₀ +
      -(8 * (n : ℚ) + 11) *
        (8 * (n : ℚ) ^ 3 + 51 * (n : ℚ) ^ 2 + 105 * (n : ℚ) + 68) * h₁

/-- The transformed Rivoal numerator has exactly the challenge initial
triple. -/
theorem rivoalOre22_P_initial :
    rivoalOre22 rivoalScaledP22 0 = 0 ∧
      rivoalOre22 rivoalScaledP22 1 = 7 ∧
      rivoalOre22 rivoalScaledP22 2 = 179 := by
  rw [rivoalOre22]
  norm_num
  rw [rivoalOre22]
  norm_num
  rw [rivoalOre22, rivoalScaledP22_three]
  norm_num

/-- The transformed Rivoal denominator has exactly the challenge initial
triple. -/
theorem rivoalOre22_Q_initial :
    rivoalOre22 rivoalScaledQ22 0 = 1 ∧
      rivoalOre22 rivoalScaledQ22 1 = 12 ∧
      rivoalOre22 rivoalScaledQ22 2 = 306 := by
  rw [rivoalOre22]
  norm_num
  rw [rivoalOre22]
  norm_num
  rw [rivoalOre22, rivoalScaledQ22_three]
  norm_num

/-- Exact identification of the challenge numerator with the transformed
Rivoal numerator. -/
theorem challengeP22_eq_rivoalOre (n : ℕ) :
    challengeP22 n = rivoalOre22 rivoalScaledP22 n := by
  apply recurrence22_unique challengeP22 (rivoalOre22 rivoalScaledP22)
    (challengeP22_zero.trans rivoalOre22_P_initial.1.symm)
    (challengeP22_one.trans rivoalOre22_P_initial.2.1.symm)
    (challengeP22_two.trans rivoalOre22_P_initial.2.2.symm)
    challengeP22_recurrence
    (rivoalOre22_satisfies rivoalScaledP22 rivoalScaledP22_recurrence)
    n

/-- Exact identification of the challenge denominator with the transformed
Rivoal denominator. -/
theorem challengeQ22_eq_rivoalOre (n : ℕ) :
    challengeQ22 n = rivoalOre22 rivoalScaledQ22 n := by
  apply recurrence22_unique challengeQ22 (rivoalOre22 rivoalScaledQ22)
    (challengeQ22_zero.trans rivoalOre22_Q_initial.1.symm)
    (challengeQ22_one.trans rivoalOre22_Q_initial.2.1.symm)
    (challengeQ22_two.trans rivoalOre22_Q_initial.2.2.symm)
    challengeQ22_recurrence
    (rivoalOre22_satisfies rivoalScaledQ22 rivoalScaledQ22_recurrence)
    n

theorem challengeQ22_pos (n : ℕ) : 0 < challengeQ22 n := by
  rw [challengeQ22_eq_rivoalOre, rivoalOre22]
  have h₀ := rivoalScaledQ22_pos n
  have h₁ := rivoalScaledQ22_pos (n + 1)
  have ha :
      0 ≤ ((n : ℚ) + 1) * (3 * (n : ℚ) + 4) := by
    positivity
  exact div_pos
    (add_pos_of_pos_of_nonneg h₁ (mul_nonneg ha h₀.le))
    (by positivity)

theorem challengeQ22_ne_zero (n : ℕ) : challengeQ22 n ≠ 0 :=
  (challengeQ22_pos n).ne'

/-! ## The target constant and the still-open analytic step -/

/-- The Mathlib Euler--Mascheroni constant is `1 - Γ'(2)`.  This is the
identity from `Ripple.Number.EulerGamma` specialized to the part that depends
only on Mathlib. -/
theorem euler_gamma_eq_one_sub_deriv_Gamma_two :
    eulerMascheroniConstant = 1 - deriv Gamma 2 := by
  have h := deriv_Gamma_nat 1
  simp only [Nat.factorial, harmonic_succ, harmonic_zero, Nat.cast_one] at h
  norm_num at h
  linarith

/-- The elementary estimate that turns an absolute first moment into
convergence of a finite weighted average. -/
theorem abs_weightedAverage_sub_le
    {ι : Type*} (s : Finset ι) (w v : ι → ℝ) (L : ℝ)
    (hw : ∀ i ∈ s, 0 ≤ w i)
    (hW : 0 < ∑ i ∈ s, w i) :
    |(∑ i ∈ s, w i * v i) / (∑ i ∈ s, w i) - L| ≤
      (∑ i ∈ s, w i * |v i - L|) / (∑ i ∈ s, w i) := by
  have hdiff :
      (∑ i ∈ s, w i * (v i - L)) =
        (∑ i ∈ s, w i * v i) - L * (∑ i ∈ s, w i) := by
    calc
      (∑ i ∈ s, w i * (v i - L)) =
          ∑ i ∈ s, (w i * v i - w i * L) := by
        apply Finset.sum_congr rfl
        intro i _
        ring
      _ = (∑ i ∈ s, w i * v i) -
          (∑ i ∈ s, w i * L) := by
        rw [Finset.sum_sub_distrib]
      _ = (∑ i ∈ s, w i * v i) -
          L * (∑ i ∈ s, w i) := by
        rw [← Finset.sum_mul]
        ring
  have hratio :
      (∑ i ∈ s, w i * v i) / (∑ i ∈ s, w i) - L =
        (∑ i ∈ s, w i * (v i - L)) / (∑ i ∈ s, w i) := by
    rw [hdiff]
    field_simp [hW.ne']
  rw [hratio, abs_div, abs_of_pos hW]
  apply div_le_div_of_nonneg_right _ hW.le
  calc
    |∑ i ∈ s, w i * (v i - L)| ≤
        ∑ i ∈ s, |w i * (v i - L)| :=
      Finset.abs_sum_le_sum_abs _ _
    _ = ∑ i ∈ s, w i * |v i - L| := by
      apply Finset.sum_congr rfl
      intro i hi
      rw [abs_mul, abs_of_nonneg (hw i hi)]

def rivoalRealWeight22 (n k : ℕ) : ℝ :=
  ((rivoalExplicitQTerm22 n k : ℚ) : ℝ)

def rivoalRealApproxValue22 (n k : ℕ) : ℝ :=
  ((rivoalApproxValue22 n k : ℚ) : ℝ)

def rivoalRealHarmonicValue22 (n k : ℕ) : ℝ :=
  ((rivoalHarmonicKernel22 n k : ℚ) : ℝ)

def rivoalCorrectionBound22 (n : ℕ) : ℝ :=
  1 / (2 * (n : ℝ) + 1)

theorem rivoalRealApproxValue22_eq (n k : ℕ) :
    rivoalRealApproxValue22 n k =
      rivoalRealHarmonicValue22 n k -
        1 / (2 * (n : ℝ) + (k : ℝ) + 1) := by
  simp only [rivoalRealApproxValue22, rivoalApproxValue22,
    rivoalRealHarmonicValue22, Rat.cast_sub, Rat.cast_div,
    Rat.cast_one, Rat.cast_add, Rat.cast_mul, Rat.cast_natCast]
  norm_num

theorem rivoalRealWeight22_nonneg (n k : ℕ) :
    0 ≤ rivoalRealWeight22 n k := by
  simp only [rivoalRealWeight22, rivoalExplicitQTerm22,
    rivoalWeight22_eq_choose, Rat.cast_mul, Rat.cast_add,
    Rat.cast_natCast, Rat.cast_one, Rat.cast_div, Rat.cast_pow]
  positivity

theorem rivoalRealWeight22_sum (n : ℕ) :
    (∑ k ∈ Finset.range (n + 1), rivoalRealWeight22 n k) =
      ((rivoalExplicitQ22 n : ℚ) : ℝ) := by
  simp only [rivoalRealWeight22, rivoalExplicitQ22, Rat.cast_sum]

theorem rivoalRealWeightedValue22_sum (n : ℕ) :
    (∑ k ∈ Finset.range (n + 1),
      rivoalRealWeight22 n k * rivoalRealApproxValue22 n k) =
      ((rivoalExplicitP22 n : ℚ) : ℝ) := by
  rw [rivoalExplicitP22_eq_weighted]
  simp only [rivoalRealWeight22, rivoalRealApproxValue22,
    Rat.cast_sum, Rat.cast_mul]

/-- The one analytic quantity left by the exact finite-sum reduction.  Its
vanishing says that the positive hypergeometric weights concentrate where
the harmonic sample value is close to `γ`. -/
def rivoalWeightedError22 (n : ℕ) : ℝ :=
  (∑ k ∈ Finset.range (n + 1),
      rivoalRealWeight22 n k *
        |rivoalRealApproxValue22 n k - eulerMascheroniConstant|) /
    ((rivoalExplicitQ22 n : ℚ) : ℝ)

def rivoalWeightedHarmonicError22 (n : ℕ) : ℝ :=
  (∑ k ∈ Finset.range (n + 1),
      rivoalRealWeight22 n k *
        |rivoalRealHarmonicValue22 n k - eulerMascheroniConstant|) /
    ((rivoalExplicitQ22 n : ℚ) : ℝ)

theorem rivoalWeightedError22_nonneg (n : ℕ) :
    0 ≤ rivoalWeightedError22 n := by
  apply div_nonneg
  · apply Finset.sum_nonneg
    intro k _
    exact mul_nonneg (rivoalRealWeight22_nonneg n k) (abs_nonneg _)
  · exact (by exact_mod_cast (rivoalExplicitQ22_pos n).le)

theorem rivoalRealApproxValue22_error_le (n k : ℕ) :
    |rivoalRealApproxValue22 n k - eulerMascheroniConstant| ≤
      |rivoalRealHarmonicValue22 n k - eulerMascheroniConstant| +
        rivoalCorrectionBound22 n := by
  rw [rivoalRealApproxValue22_eq]
  calc
    |(rivoalRealHarmonicValue22 n k -
          1 / (2 * (n : ℝ) + (k : ℝ) + 1)) -
        eulerMascheroniConstant| =
        |(rivoalRealHarmonicValue22 n k -
            eulerMascheroniConstant) +
          (-1 / (2 * (n : ℝ) + (k : ℝ) + 1))| := by
      congr 1
      ring
    _ ≤ |rivoalRealHarmonicValue22 n k -
          eulerMascheroniConstant| +
        |-1 / (2 * (n : ℝ) + (k : ℝ) + 1)| :=
      abs_add_le _ _
    _ = |rivoalRealHarmonicValue22 n k -
          eulerMascheroniConstant| +
        1 / (2 * (n : ℝ) + (k : ℝ) + 1) := by
      rw [show -1 / (2 * (n : ℝ) + (k : ℝ) + 1) =
          -(1 / (2 * (n : ℝ) + (k : ℝ) + 1)) by ring,
        abs_neg,
        abs_of_pos (by
          positivity :
          0 < 1 / (2 * (n : ℝ) + (k : ℝ) + 1))]
    _ ≤ |rivoalRealHarmonicValue22 n k -
          eulerMascheroniConstant| +
        rivoalCorrectionBound22 n := by
      simp only [rivoalCorrectionBound22]
      have hrec :
          1 / (2 * (n : ℝ) + (k : ℝ) + 1) ≤
            1 / (2 * (n : ℝ) + 1) := by
        apply one_div_le_one_div_of_le
        · positivity
        · have hk : (0 : ℝ) ≤ k := by positivity
          linarith
      exact add_le_add_right hrec _

theorem rivoalWeightedError22_le_harmonic (n : ℕ) :
    rivoalWeightedError22 n ≤
      rivoalWeightedHarmonicError22 n + rivoalCorrectionBound22 n := by
  have hQ :
      0 < ((rivoalExplicitQ22 n : ℚ) : ℝ) := by
    exact_mod_cast rivoalExplicitQ22_pos n
  have hsum :
      (∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k *
            |rivoalRealApproxValue22 n k - eulerMascheroniConstant|) ≤
        ∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k *
            (|rivoalRealHarmonicValue22 n k -
                eulerMascheroniConstant| +
              rivoalCorrectionBound22 n) := by
    apply Finset.sum_le_sum
    intro k _
    exact mul_le_mul_of_nonneg_left
      (rivoalRealApproxValue22_error_le n k)
      (rivoalRealWeight22_nonneg n k)
  rw [rivoalWeightedError22, rivoalWeightedHarmonicError22]
  apply (div_le_div_of_nonneg_right hsum hQ.le).trans_eq
  have hsplit :
      (∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k *
            (|rivoalRealHarmonicValue22 n k -
                eulerMascheroniConstant| +
              rivoalCorrectionBound22 n)) =
        (∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k *
            |rivoalRealHarmonicValue22 n k -
              eulerMascheroniConstant|) +
        rivoalCorrectionBound22 n *
          ((rivoalExplicitQ22 n : ℚ) : ℝ) := by
    calc
      _ = (∑ k ∈ Finset.range (n + 1),
            (rivoalRealWeight22 n k *
                |rivoalRealHarmonicValue22 n k -
                  eulerMascheroniConstant| +
              rivoalRealWeight22 n k * rivoalCorrectionBound22 n)) := by
          apply Finset.sum_congr rfl
          intro k _
          ring
      _ = (∑ k ∈ Finset.range (n + 1),
            rivoalRealWeight22 n k *
              |rivoalRealHarmonicValue22 n k -
                eulerMascheroniConstant|) +
          ∑ k ∈ Finset.range (n + 1),
            rivoalRealWeight22 n k * rivoalCorrectionBound22 n := by
          rw [Finset.sum_add_distrib]
      _ = _ := by
          rw [← Finset.sum_mul, rivoalRealWeight22_sum]
          ring
  rw [hsplit]
  field_simp [hQ.ne']

theorem rivoalCorrectionBound22_tendsto :
    Tendsto rivoalCorrectionBound22 atTop (𝓝 0) := by
  have hden :
      Tendsto (fun n : ℕ => 2 * (n : ℝ) + 1) atTop atTop := by
    exact
      (tendsto_natCast_atTop_atTop.const_mul_atTop
        (by norm_num : (0 : ℝ) < 2)).atTop_add tendsto_const_nhds
  unfold rivoalCorrectionBound22
  simpa only [one_div, Function.comp_apply] using
    tendsto_inv_atTop_zero.comp hden

theorem rivoalExplicitRatio22_error_le (n : ℕ) :
    |((rivoalExplicitP22 n : ℚ) : ℝ) /
          ((rivoalExplicitQ22 n : ℚ) : ℝ) -
        eulerMascheroniConstant| ≤
      rivoalWeightedError22 n := by
  have hW :
      0 < ∑ k ∈ Finset.range (n + 1), rivoalRealWeight22 n k := by
    rw [rivoalRealWeight22_sum]
    exact_mod_cast rivoalExplicitQ22_pos n
  have h := abs_weightedAverage_sub_le
    (Finset.range (n + 1))
    (rivoalRealWeight22 n) (rivoalRealApproxValue22 n)
    eulerMascheroniConstant
    (by
      intro k _
      exact rivoalRealWeight22_nonneg n k)
    hW
  rw [rivoalRealWeight22_sum, rivoalRealWeightedValue22_sum] at h
  exact h

/-- Minimal concentration/harmonic interface for the analytic part of
Problem 2.2. -/
def RivoalConcentrationHarmonicClaim22 : Prop :=
  Tendsto rivoalWeightedError22 atTop (𝓝 0)

/-- After the elementary rational correction is removed, only the weighted
harmonic concentration remains. -/
def RivoalHarmonicConcentrationClaim22 : Prop :=
  Tendsto rivoalWeightedHarmonicError22 atTop (𝓝 0)

theorem rivoalConcentrationHarmonicClaim22_of_harmonic
    (h : RivoalHarmonicConcentrationClaim22) :
    RivoalConcentrationHarmonicClaim22 := by
  rw [RivoalHarmonicConcentrationClaim22] at h
  rw [RivoalConcentrationHarmonicClaim22]
  apply squeeze_zero
  · exact rivoalWeightedError22_nonneg
  · exact rivoalWeightedError22_le_harmonic
  · simpa using h.add rivoalCorrectionBound22_tendsto

theorem rivoalExplicitRatio22_tendsto_of_concentration
    (h : RivoalConcentrationHarmonicClaim22) :
    Tendsto
      (fun n : ℕ =>
        ((rivoalExplicitP22 n : ℚ) : ℝ) /
          ((rivoalExplicitQ22 n : ℚ) : ℝ))
      atTop (𝓝 eulerMascheroniConstant) := by
  rw [RivoalConcentrationHarmonicClaim22] at h
  rw [Metric.tendsto_atTop] at h ⊢
  intro ε hε
  obtain ⟨N, hN⟩ := h ε hε
  refine ⟨N, ?_⟩
  intro n hn
  have herr := hN n hn
  rw [Real.dist_eq] at herr ⊢
  have herrlt : rivoalWeightedError22 n < ε := by
    simpa [sub_zero, abs_of_nonneg (rivoalWeightedError22_nonneg n)] using
      herr
  exact lt_of_le_of_lt (rivoalExplicitRatio22_error_le n) herrlt

/-- The exact limit assertion posed by Problem 2.2.  It is stated here so that
future work cannot accidentally prove a limit about unrelated witnesses. -/
def Problem22Claim : Prop :=
  Tendsto
    (fun n : ℕ =>
      ((challengeP22 (n + 3) : ℚ) : ℝ) /
        ((challengeQ22 (n + 3) : ℚ) : ℝ))
    atTop (𝓝 eulerMascheroniConstant)

/-- The analytic limit assertion for the exact factorial-scaled Rivoal
sequences used above.  Multiplying numerator and denominator by `(n!)²`
does not change their ratio. -/
def RivoalLimitClaim22 : Prop :=
  Tendsto
    (fun n : ℕ =>
      ((rivoalScaledP22 n : ℚ) : ℝ) /
        ((rivoalScaledQ22 n : ℚ) : ℝ))
    atTop (𝓝 eulerMascheroniConstant)

theorem rivoalLimitClaim22_of_concentration
    (h : RivoalConcentrationHarmonicClaim22) :
    RivoalLimitClaim22 := by
  have hexplicit :=
    rivoalExplicitRatio22_tendsto_of_concentration h
  rw [RivoalLimitClaim22]
  apply hexplicit.congr'
  filter_upwards [] with n
  norm_cast
  exact (rivoalScaledRatio22_eq_explicit n).symm

theorem rivoalLimitClaim22_of_harmonic_concentration
    (h : RivoalHarmonicConcentrationClaim22) :
    RivoalLimitClaim22 :=
  rivoalLimitClaim22_of_concentration
    (rivoalConcentrationHarmonicClaim22_of_harmonic h)

/-- Positivity of Rivoal's exact denominator sequence.  This is separated
from `RivoalLimitClaim22` because it is the algebraic hypothesis needed to
show that the Ore transform preserves the ratio limit. -/
def RivoalDenominatorPositivity22 : Prop :=
  ∀ n : ℕ, 0 < ((rivoalScaledQ22 n : ℚ) : ℝ)

/-- Rivoal's exact denominator is positive at every index. -/
theorem rivoalDenominatorPositivity22 :
    RivoalDenominatorPositivity22 := by
  intro n
  exact_mod_cast rivoalScaledQ22_pos n

/-- A positive weighted average of two adjacent ratios has the same limit as
the original ratios.  This is the analytic mechanism behind preservation of
the Rivoal limit by the first-order Ore transform. -/
theorem tendsto_adjacent_weighted_ratio22
    (p q a : ℕ → ℝ) (L : ℝ)
    (hq : ∀ n, 0 < q n) (ha : ∀ n, 0 ≤ a n)
    (hlim : Tendsto (fun n => p n / q n) atTop (𝓝 L)) :
    Tendsto
      (fun n => (p (n + 1) + a n * p n) /
        (q (n + 1) + a n * q n))
      atTop (𝓝 L) := by
  rw [Metric.tendsto_atTop] at hlim ⊢
  intro ε hε
  obtain ⟨N, hN⟩ := hlim ε hε
  refine ⟨N, ?_⟩
  intro n hn
  have hn₁ : N ≤ n + 1 := by omega
  have h₀ := hN n hn
  have h₁ := hN (n + 1) hn₁
  rw [Real.dist_eq] at h₀ h₁ ⊢
  rw [abs_lt] at h₀ h₁ ⊢
  have hq₀ : 0 < q n := hq n
  have hq₁ : 0 < q (n + 1) := hq (n + 1)
  have ha₀ : 0 ≤ a n := ha n
  have hden : 0 < q (n + 1) + a n * q n :=
    add_pos_of_pos_of_nonneg hq₁ (mul_nonneg ha₀ hq₀.le)
  constructor
  · have hl₀ : L - ε < p n / q n := by linarith [h₀.1]
    have hl₁ : L - ε < p (n + 1) / q (n + 1) := by linarith [h₁.1]
    have hl₀' : (L - ε) * q n < p n := (lt_div_iff₀ hq₀).mp hl₀
    have hl₁' : (L - ε) * q (n + 1) < p (n + 1) :=
      (lt_div_iff₀ hq₁).mp hl₁
    have hmul₀ : a n * ((L - ε) * q n) ≤ a n * p n :=
      mul_le_mul_of_nonneg_left hl₀'.le ha₀
    rw [lt_sub_iff_add_lt]
    apply (lt_div_iff₀ hden).2
    nlinarith
  · have hu₀ : p n / q n < L + ε := by linarith [h₀.2]
    have hu₁ : p (n + 1) / q (n + 1) < L + ε := by linarith [h₁.2]
    have hu₀' : p n < (L + ε) * q n := (div_lt_iff₀ hq₀).mp hu₀
    have hu₁' : p (n + 1) < (L + ε) * q (n + 1) :=
      (div_lt_iff₀ hq₁).mp hu₁
    have hmul₀ : a n * p n ≤ a n * ((L + ε) * q n) :=
      mul_le_mul_of_nonneg_left hu₀'.le ha₀
    rw [sub_lt_iff_lt_add]
    apply (div_lt_iff₀ hden).2
    nlinarith

/-- Conditional closure of Problem 2.2 from the exact limit statement for
Rivoal's approximants.  Everything in the transfer from those fixed sequences
to the challenge sequences, including denominator positivity, is proved in
Lean; no citation is masquerading as a formal proof. -/
theorem problem22_of_rivoal
    (hlim : RivoalLimitClaim22) :
    Problem22Claim := by
  let p : ℕ → ℝ := fun n => ((rivoalScaledP22 n : ℚ) : ℝ)
  let q : ℕ → ℝ := fun n => ((rivoalScaledQ22 n : ℚ) : ℝ)
  let a : ℕ → ℝ := fun n => ((n : ℝ) + 1) * (3 * (n : ℝ) + 4)
  have hq : ∀ n, 0 < q n := by
    intro n
    exact rivoalDenominatorPositivity22 n
  have ha : ∀ n, 0 ≤ a n := by
    intro n
    dsimp [a]
    positivity
  have hbase : Tendsto (fun n => p n / q n) atTop
      (𝓝 eulerMascheroniConstant) := by
    exact hlim
  have hw := tendsto_adjacent_weighted_ratio22
    p q a eulerMascheroniConstant hq ha hbase
  have hOre : ∀ n : ℕ,
      (((rivoalOre22 rivoalScaledP22 n : ℚ) : ℝ) /
          ((rivoalOre22 rivoalScaledQ22 n : ℚ) : ℝ)) =
        ((p (n + 1) + a n * p n) /
          (q (n + 1) + a n * q n)) := by
    intro n
    dsimp [p, q, a]
    simp only [rivoalOre22, Rat.cast_div, Rat.cast_add, Rat.cast_mul,
      Rat.cast_natCast]
    rw [div_div_div_cancel_right₀]
    · norm_num
    · positivity
  have hOreLimit : Tendsto
      (fun n : ℕ =>
        (((rivoalOre22 rivoalScaledP22 n : ℚ) : ℝ) /
          ((rivoalOre22 rivoalScaledQ22 n : ℚ) : ℝ)))
      atTop (𝓝 eulerMascheroniConstant) := by
    apply hw.congr'
    filter_upwards [] with n
    exact (hOre n).symm
  have hChallengeLimit : Tendsto
      (fun n : ℕ =>
        (((challengeP22 n : ℚ) : ℝ) /
          ((challengeQ22 n : ℚ) : ℝ)))
      atTop (𝓝 eulerMascheroniConstant) := by
    apply hOreLimit.congr'
    filter_upwards [] with n
    rw [challengeP22_eq_rivoalOre, challengeQ22_eq_rivoalOre]
  rw [Problem22Claim]
  exact (tendsto_add_atTop_iff_nat 3).2 hChallengeLimit

theorem problem22_of_harmonic_concentration
    (h : RivoalHarmonicConcentrationClaim22) :
    Problem22Claim :=
  problem22_of_rivoal
    (rivoalLimitClaim22_of_harmonic_concentration h)

/-! ## Status of the main statement

The exact limit is `Problem22Claim`.  The explicit hypergeometric formulas,
their WZ certificates (including all boundary terms), the Ore identification,
the complete limit-transfer argument, and denominator positivity are proved
above.  The remaining formalization boundary is
`RivoalHarmonicConcentrationClaim22`: the positive weights

```
(2n+k+1) * choose(n,k)^2 / k!
```

must concentrate where `3 H_k - 2 H_(n-k)` is close to `γ`.  The harmless
`1 / (2n+k+1)` correction has already been bounded and removed in Lean.

An earlier version of this file carried a theorem of the shape

```
theorem ..._identity : ∃ (p q : ℕ → ℝ), Tendsto (fun n => p n / q n) atTop (𝓝 L) :=
  ⟨fun _ => L, fun _ => 1, by simp⟩
```

which is vacuous: it is witnessed by constant sequences and says nothing about
the challenge recurrence.  It has been removed rather than shipped.  The
present file instead proves the exact algebraic reduction to Rivoal's
approximants and exposes the one remaining analytic obligation as a named
proposition.
-/

end RamanujanChallenge.P22

end
