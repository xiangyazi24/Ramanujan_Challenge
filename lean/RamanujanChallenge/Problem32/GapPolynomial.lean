/-
  Problem 3.2 — Apéry gap polynomials.

  The polynomial `gapPolynomial h` is the numerator obtained by propagating
  the Apéry recurrence across a gap of length `h`.  We include the harmless
  seed `N₀ = 0`; then the recurrence is uniform from `h = 1` onward and gives
  `N₁ = 1`, `N₂(X) = P(X+1)`.
-/
import RamanujanChallenge.Problem32.WZCertificate
import Mathlib.Algebra.Polynomial.Degree.Operations
import Mathlib.Algebra.Polynomial.Eval.Defs

noncomputable section

open Polynomial

/-! ## The shifted Apéry coefficient polynomial -/

/-- The integer middle coefficient in the Apéry recurrence. -/
def aperyPInt (n : ℤ) : ℤ :=
  34 * n ^ 3 + 51 * n ^ 2 + 27 * n + 5

/-- `P(X+s)`, regarded as a polynomial over the integers. -/
def shiftedAperyP (s : ℤ) : ℤ[X] :=
  34 * (X + C s) ^ 3 + 51 * (X + C s) ^ 2 + 27 * (X + C s) + 5

@[simp] theorem eval_shiftedAperyP (s x : ℤ) :
    (shiftedAperyP s).eval x = aperyPInt (x + s) := by
  simp [shiftedAperyP, aperyPInt]

@[simp] theorem natDegree_shiftedAperyP (s : ℤ) :
    (shiftedAperyP s).natDegree = 3 := by
  unfold shiftedAperyP
  compute_degree <;> norm_num

/-! ## Gap-polynomial recurrence -/

/-- The cleared numerator for an Apéry recurrence gap of length `h`.

The `N₀ = 0` seed is only an indexing convenience.  The mathematical family
used in the paper begins at `N₁ = 1`.
-/
def gapPolynomial : ℕ → ℤ[X]
  | 0 => 0
  | 1 => 1
  | h + 2 =>
      shiftedAperyP (h + 1) * gapPolynomial (h + 1) -
        (X + C (h + 1 : ℤ)) ^ 6 * gapPolynomial h

@[simp] theorem gapPolynomial_zero : gapPolynomial 0 = 0 := rfl

@[simp] theorem gapPolynomial_one : gapPolynomial 1 = 1 := rfl

@[simp] theorem gapPolynomial_two : gapPolynomial 2 = shiftedAperyP 1 := by
  simp [gapPolynomial]

/-- The paper's uniform `h`-direction recurrence. -/
theorem gapPolynomial_succ (h : ℕ) (hh : 1 ≤ h) :
    gapPolynomial (h + 1) =
      shiftedAperyP h * gapPolynomial h -
        (X + C (h : ℤ)) ^ 6 * gapPolynomial (h - 1) := by
  cases h with
  | zero => omega
  | succ k => simp [gapPolynomial]

/-- Pointwise form of the gap-polynomial recurrence. -/
theorem eval_gapPolynomial_succ (h : ℕ) (hh : 1 ≤ h) (x : ℤ) :
    (gapPolynomial (h + 1)).eval x =
      aperyPInt (x + h) * (gapPolynomial h).eval x -
        (x + h) ^ 6 * (gapPolynomial (h - 1)).eval x := by
  rw [gapPolynomial_succ h hh]
  simp

/-! ## Degree bound -/

/-- The root-counting argument only needs this upper bound.  It deliberately
avoids using the manuscript's informal leading-coefficient argument. -/
theorem gapPolynomial_natDegree_le (h : ℕ) :
    (gapPolynomial h).natDegree ≤ 3 * (h - 1) := by
  induction h using Nat.strong_induction_on with
  | h h ih =>
      cases h with
      | zero => simp
      | succ h =>
          cases h with
          | zero => simp
          | succ h =>
              cases h with
              | zero => simp [gapPolynomial]
              | succ k =>
                  rw [gapPolynomial_succ (k + 2) (by omega)]
                  apply (Polynomial.natDegree_sub_le _ _).trans
                  apply max_le
                  · exact Polynomial.natDegree_mul_le.trans (by
                      rw [natDegree_shiftedAperyP]
                      have hprev := ih (k + 2) (by omega)
                      omega)
                  · exact Polynomial.natDegree_mul_le.trans (by
                      rw [Polynomial.natDegree_pow_X_add_C]
                      have hprev := ih (k + 1) (by omega)
                      have hprev' :
                          (gapPolynomial (k + 2 - 1)).natDegree ≤
                            3 * ((k + 1) - 1) := by
                        simpa only [Nat.add_sub_cancel] using hprev
                      omega)

end
