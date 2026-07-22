/-
  Ramanujan Challenge Problem 2.4: Polylogarithm Identity

  The double sum identity:
    Σ_{m≥0} Σ_{k=0}^{m} C(m,k)² H_k² / ((m+1)² C(2m,m))
    = 20·Li₄(1/2) + (5/6)·log⁴2 + 10·ζ(2) − (65/9)·ζ(2)²
      − log²2·(12 + 5ζ(2)) + ½·ζ(3) + log2·(35/2·ζ(3) − 16)

  where H_k = Σ_{j=1}^{k} 1/j is the k-th harmonic number.

  This is a CLOSED-FORM identity (not a convergence result).
  The proof uses the Wilf-Zeilberger method for the inner sum,
  polylogarithm ladder identities, and known special values.

  Reference: Xiang Huang, "Solution to Problem 2.4", July 2026.
-/
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Data.Nat.Choose.Basic

noncomputable section

open Real

/-! ## Harmonic numbers -/

def harmonicNumber (k : ℕ) : ℝ := ∑ j ∈ Finset.range k, (1 : ℝ) / (↑j + 1)

theorem harmonicNumber_zero : harmonicNumber 0 = 0 := by
  simp [harmonicNumber]

theorem harmonicNumber_one : harmonicNumber 1 = 1 := by
  simp [harmonicNumber, Finset.sum_range_one]

/-! ## Polylogarithm Li₄ -/

def polylog4 (z : ℝ) : ℝ := ∑' n : ℕ, z ^ (n + 1) / (↑(n + 1) : ℝ) ^ 4

/-! ## The double sum (LHS) -/

def lhs_24 : ℝ := ∑' m : ℕ, ∑ k ∈ Finset.range (m + 1),
  (Nat.choose m k : ℝ) ^ 2 * (harmonicNumber k) ^ 2 /
  ((↑m + 1) ^ 2 * (Nat.choose (2 * m) m : ℝ))

/-! ## The RHS -/

def zeta3_24 : ℝ := ∑' n : ℕ, (1 : ℝ) / (↑n + 1) ^ 3

def rhs_24 : ℝ :=
  20 * polylog4 (1/2) + (5/6) * (Real.log 2) ^ 4 +
  10 * (Real.pi ^ 2 / 6) - (65/9) * (Real.pi ^ 2 / 6) ^ 2 -
  (Real.log 2) ^ 2 * (12 + 5 * (Real.pi ^ 2 / 6)) +
  (1/2) * zeta3_24 +
  Real.log 2 * ((35/2) * zeta3_24 - 16)

/-! ## Main theorem -/

theorem problem24_identity : lhs_24 = rhs_24 := by
  sorry

end
