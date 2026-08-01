/-
  Problem 2.5 — exact algebra behind the proposed Delannoy error split.

  This file deliberately separates the identities that follow from triangular
  inversion from the two analytic inputs still requiring proof: an all-index
  recurrence certificate and a quantitative Poincare theorem for that
  recurrence.  In particular, a numerically guessed recurrence is not used as
  a theorem here.
-/
import RamanujanChallenge.Problem25Delannoy
import RamanujanChallenge.Problem25TightBounds

noncomputable section

namespace RamanujanChallenge.P25

open Nat Finset

/-- The coefficientwise ratio error from the proposed Delannoy argument. -/
def delannoyEpsilon (k : ℕ) : ℝ :=
  (delannoyG k : ℝ) / (delannoyF k : ℝ) - catalanConstant

/-- The division-free error coefficient. -/
def delannoyErrorCoefficient (k : ℕ) : ℝ :=
  (delannoyG k : ℝ) - catalanConstant * (delannoyF k : ℝ)

theorem delannoyErrorCoefficient_eq_mul_epsilon (k : ℕ)
    (hf : delannoyF k ≠ 0) :
    delannoyErrorCoefficient k = (delannoyF k : ℝ) * delannoyEpsilon k := by
  have hfr : (delannoyF k : ℝ) ≠ 0 := by exact_mod_cast hf
  simp only [delannoyErrorCoefficient, delannoyEpsilon]
  field_simp

/-- Exact Delannoy splitting of the normalized first-column error. -/
theorem normalized_error_decomposition (N : ℕ) :
    (normalizedNumerator N : ℝ) -
        catalanConstant * (normalizedDenominator N : ℝ) =
      ∑ k ∈ range (N + 1),
        delannoyErrorCoefficient k * (delannoyB N k : ℝ) := by
  have hp := congrArg (fun x : ℚ => (x : ℝ))
    (normalizedNumerator_decomposition N)
  have hq := congrArg (fun x : ℚ => (x : ℝ))
    (normalizedDenominator_decomposition N)
  push_cast at hp hq
  rw [hp, hq, mul_sum, ← sum_sub_distrib]
  apply sum_congr rfl
  intro k hk
  simp only [delannoyErrorCoefficient]
  ring

/-- Bounds on `G` turn into bounds on a coefficient ratio without any
asymptotic argument. -/
theorem delannoyEpsilon_bounds (k : ℕ) (lower upper : ℝ)
    (hlower : lower < catalanConstant)
    (hupper : catalanConstant < upper) :
    (delannoyG k : ℝ) / (delannoyF k : ℝ) - upper < delannoyEpsilon k ∧
      delannoyEpsilon k <
        (delannoyG k : ℝ) / (delannoyF k : ℝ) - lower := by
  simp only [delannoyEpsilon]
  constructor <;> linarith

/-- The first coefficient receives a fully rational rigorous error interval. -/
theorem delannoyEpsilon_zero_bounds :
    (30921 : ℝ) / 33750 - 21390206625 / 23352603750 < delannoyEpsilon 0 ∧
      delannoyEpsilon 0 <
        (30921 : ℝ) / 33750 - 1590511050 / 1736437500 := by
  simpa using delannoyEpsilon_bounds 0
    ((1590511050 : ℝ) / 1736437500)
    ((21390206625 : ℝ) / 23352603750)
    catalan_tight_lower catalan_tight_upper

/-! ## The guessed operator's limiting polynomial -/

/-- Coefficient form of the limiting polynomial printed for the guessed
order-nine LCLM.  This definition records only the checked polynomial
arithmetic; it does not assert that the guessed finite-data operator
annihilates `delannoyF` or `delannoyG` for every index. -/
def kPoincarePolynomial (xi : ℚ) : ℚ :=
  xi ^ 9 - (9 / 8) * xi ^ 8 - (27 / 64) * xi ^ 7 +
    (147 / 512) * xi ^ 6 + (819 / 4096) * xi ^ 5 +
    (1701 / 32768) * xi ^ 4 + (1911 / 262144) * xi ^ 3 +
    (1233 / 2097152) * xi ^ 2 + (27 / 1048576) * xi +
    1 / 2097152

theorem kPoincarePolynomial_factor (xi : ℚ) :
    kPoincarePolynomial xi = (xi - 1) ^ 2 * (xi + 1 / 8) ^ 7 := by
  unfold kPoincarePolynomial
  ring

@[simp] theorem kPoincarePolynomial_one : kPoincarePolynomial 1 = 0 := by
  rw [kPoincarePolynomial_factor]
  norm_num

@[simp] theorem kPoincarePolynomial_neg_one_eighth :
    kPoincarePolynomial (-1 / 8) = 0 := by
  rw [kPoincarePolynomial_factor]
  norm_num

theorem neg_one_eighth_abs_lt_one : |(-1 / 8 : ℝ)| < 1 := by norm_num

/-! ## Exact consequence of the proposed dominant-mode cancellation -/

/-- If a nonzero rational dominant projection `A` and a rational projection
`B` satisfied the cancellation equation claimed in the route specification,
then Catalan's constant would equal the rational number `B/A`.

Thus finite rational Casorati data cannot by itself establish cancellation:
the missing premise is precisely an exact evaluation of Catalan's constant as
that rational number. -/
theorem catalan_eq_rational_of_dominant_cancellation (A B : ℚ) (hA : A ≠ 0)
    (hcancel : (B : ℝ) - catalanConstant * (A : ℝ) = 0) :
    catalanConstant = ((B / A : ℚ) : ℝ) := by
  have hAr : (A : ℝ) ≠ 0 := by exact_mod_cast hA
  rw [Rat.cast_div]
  apply (eq_div_iff hAr).2
  linarith

theorem catalan_is_rational_of_dominant_cancellation (A B : ℚ) (hA : A ≠ 0)
    (hcancel : (B : ℝ) - catalanConstant * (A : ℝ) = 0) :
    ∃ q : ℚ, catalanConstant = (q : ℝ) := by
  exact ⟨B / A, catalan_eq_rational_of_dominant_cancellation A B hA hcancel⟩

end RamanujanChallenge.P25

end
