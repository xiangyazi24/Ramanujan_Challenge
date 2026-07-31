import RamanujanChallenge.Problem25
import Mathlib.Analysis.SpecificLimits.Basic
import Mathlib.Topology.Order.MonotoneConvergence

noncomputable section

namespace RamanujanChallenge.P25

/-!
# Problem 2.5: unconditional projective contraction

This file proves that the three columns in Problem 2.5 converge
unconditionally to one common real number.  The proof uses an explicit
invariant cone for the positive denominator row.  Inside this cone, each
normalized column weight gives mass at least `1 / 6` to each of the first two
old columns.  Removing these two common masses leaves total mass `2 / 3`, so
the common enclosing interval contracts by that factor at every step.

Thus the remaining Catalan connection is a single scalar identity: the common
limit constructed below equals `catalanConstant`.
-/

private theorem affine_rectangle_nonneg
    (A B C x y : ℝ)
    (hx₀ : 1 ≤ x) (hx₁ : x ≤ 3 / 2)
    (hy₀ : 0 ≤ y) (hy₁ : y ≤ 2)
    (h₀₀ : 0 ≤ A + B * 1 + C * 0)
    (h₀₁ : 0 ≤ A + B * 1 + C * 2)
    (h₁₀ : 0 ≤ A + B * (3 / 2) + C * 0)
    (h₁₁ : 0 ≤ A + B * (3 / 2) + C * 2) :
    0 ≤ A + B * x + C * y := by
  rcases le_total 0 B with hB | hB <;>
    rcases le_total 0 C with hC | hC <;> nlinarith

/-- The first projective coordinate of the positive denominator row, with the
balancing power forced by the degrees of the matrix entries. -/
def denominatorX (n : ℕ) : ℝ :=
  ((n : ℝ) + 1) * (positiveDenominator n 1 : ℝ) /
    (positiveDenominator n 0 : ℝ)

/-- The second projective coordinate of the positive denominator row. -/
def denominatorY (n : ℕ) : ℝ :=
  ((n : ℝ) + 1) ^ 3 * (positiveDenominator n 2 : ℝ) /
    (positiveDenominator n 0 : ℝ)

private def stepCoefficient (n : ℕ) (x y : ℝ) (j : Fin 3) : ℝ :=
  (positiveMatrix (n : ℤ) 0 j : ℝ) +
    x / ((n : ℝ) + 1) * (positiveMatrix (n : ℤ) 1 j : ℝ) +
    y / ((n : ℝ) + 1) ^ 3 * (positiveMatrix (n : ℤ) 2 j : ℝ)

private theorem positiveDenominator_succ_factor (n : ℕ) (j : Fin 3) :
    (positiveDenominator (n + 1) j : ℝ) =
      (positiveDenominator n 0 : ℝ) *
        stepCoefficient n (denominatorX n) (denominatorY n) j := by
  have hq : (positiveDenominator n 0 : ℝ) ≠ 0 := by
    exact_mod_cast (positiveDenominator_pos n 0).ne'
  rw [positiveDenominator_succ]
  push_cast
  rw [Fin.sum_univ_three]
  simp only [stepCoefficient, denominatorX, denominatorY]
  field_simp

private theorem stepCoefficient_pos (n : ℕ) (x y : ℝ)
    (hx : 0 ≤ x) (hy : 0 ≤ y) (j : Fin 3) :
    0 < stepCoefficient n x y j := by
  have hn : (0 : ℝ) < (n : ℝ) + 1 := by positivity
  have h0 : (0 : ℝ) < positiveMatrix (n : ℤ) 0 j := by
    exact_mod_cast positiveMatrix_pos n 0 j
  have h1 : (0 : ℝ) < positiveMatrix (n : ℤ) 1 j := by
    exact_mod_cast positiveMatrix_pos n 1 j
  have h2 : (0 : ℝ) < positiveMatrix (n : ℤ) 2 j := by
    exact_mod_cast positiveMatrix_pos n 2 j
  have hx' :
      0 ≤ x / ((n : ℝ) + 1) * (positiveMatrix (n : ℤ) 1 j : ℝ) := by
    exact mul_nonneg (div_nonneg hx hn.le) h1.le
  have hy' :
      0 ≤ y / ((n : ℝ) + 1) ^ 3 * (positiveMatrix (n : ℤ) 2 j : ℝ) := by
    exact mul_nonneg (div_nonneg hy (pow_nonneg hn.le 3)) h2.le
  simp only [stepCoefficient]
  linarith

private theorem step_x_lower (n : ℕ) (x y : ℝ)
    (hx₀ : 1 ≤ x) (hx₁ : x ≤ 3 / 2)
    (hy₀ : 0 ≤ y) (hy₁ : y ≤ 2) :
    stepCoefficient n x y 0 ≤
      ((n : ℝ) + 2) * stepCoefficient n x y 1 := by
  let A : ℝ :=
    ((n : ℝ) + 2) * (positiveMatrix (n : ℤ) 0 1 : ℝ) -
      (positiveMatrix (n : ℤ) 0 0 : ℝ)
  let B : ℝ :=
    (((n : ℝ) + 2) * (positiveMatrix (n : ℤ) 1 1 : ℝ) -
      (positiveMatrix (n : ℤ) 1 0 : ℝ)) / ((n : ℝ) + 1)
  let C : ℝ :=
    (((n : ℝ) + 2) * (positiveMatrix (n : ℤ) 2 1 : ℝ) -
      (positiveMatrix (n : ℤ) 2 0 : ℝ)) / ((n : ℝ) + 1) ^ 3
  have hrect : 0 ≤ A + B * x + C * y := by
    apply affine_rectangle_nonneg A B C x y hx₀ hx₁ hy₀ hy₁
    all_goals
      dsimp [A, B, C]
      norm_num [positiveMatrix, Matrix.cons_val_two]
      ring_nf
      positivity
  suffices 0 ≤
      ((n : ℝ) + 2) * stepCoefficient n x y 1 -
        stepCoefficient n x y 0 by linarith
  convert hrect using 1 <;>
    simp only [stepCoefficient] <;>
    dsimp [A, B, C] <;>
    ring

private theorem step_x_upper (n : ℕ) (x y : ℝ)
    (hx₀ : 1 ≤ x) (hx₁ : x ≤ 3 / 2)
    (hy₀ : 0 ≤ y) (hy₁ : y ≤ 2) :
    2 * (((n : ℝ) + 2) * stepCoefficient n x y 1) ≤
      3 * stepCoefficient n x y 0 := by
  let A : ℝ :=
    3 * (positiveMatrix (n : ℤ) 0 0 : ℝ) -
      2 * ((n : ℝ) + 2) * (positiveMatrix (n : ℤ) 0 1 : ℝ)
  let B : ℝ :=
    (3 * (positiveMatrix (n : ℤ) 1 0 : ℝ) -
      2 * ((n : ℝ) + 2) * (positiveMatrix (n : ℤ) 1 1 : ℝ)) /
        ((n : ℝ) + 1)
  let C : ℝ :=
    (3 * (positiveMatrix (n : ℤ) 2 0 : ℝ) -
      2 * ((n : ℝ) + 2) * (positiveMatrix (n : ℤ) 2 1 : ℝ)) /
        ((n : ℝ) + 1) ^ 3
  have hrect : 0 ≤ A + B * x + C * y := by
    apply affine_rectangle_nonneg A B C x y hx₀ hx₁ hy₀ hy₁
    all_goals
      dsimp [A, B, C]
      norm_num [positiveMatrix, Matrix.cons_val_two]
      ring_nf
      positivity
  suffices 0 ≤
      3 * stepCoefficient n x y 0 -
        2 * (((n : ℝ) + 2) * stepCoefficient n x y 1) by linarith
  convert hrect using 1 <;>
    simp only [stepCoefficient] <;>
    dsimp [A, B, C] <;>
    ring

private theorem step_y_upper (n : ℕ) (x y : ℝ)
    (hx₀ : 1 ≤ x) (hx₁ : x ≤ 3 / 2)
    (hy₀ : 0 ≤ y) (hy₁ : y ≤ 2) :
    (((n : ℝ) + 2) ^ 3) * stepCoefficient n x y 2 ≤
      2 * stepCoefficient n x y 0 := by
  let A : ℝ :=
    2 * (positiveMatrix (n : ℤ) 0 0 : ℝ) -
      ((n : ℝ) + 2) ^ 3 * (positiveMatrix (n : ℤ) 0 2 : ℝ)
  let B : ℝ :=
    (2 * (positiveMatrix (n : ℤ) 1 0 : ℝ) -
      ((n : ℝ) + 2) ^ 3 * (positiveMatrix (n : ℤ) 1 2 : ℝ)) /
        ((n : ℝ) + 1)
  let C : ℝ :=
    (2 * (positiveMatrix (n : ℤ) 2 0 : ℝ) -
      ((n : ℝ) + 2) ^ 3 * (positiveMatrix (n : ℤ) 2 2 : ℝ)) /
        ((n : ℝ) + 1) ^ 3
  have hrect : 0 ≤ A + B * x + C * y := by
    apply affine_rectangle_nonneg A B C x y hx₀ hx₁ hy₀ hy₁
    ·
      dsimp [A, B, C]
      norm_num [positiveMatrix, Matrix.cons_val_two]
      ring_nf
      positivity
    · rw [show A + B * 1 + C * 2 =
          2 * (48 * (n : ℝ) ^ 10 + 2152 * (n : ℝ) ^ 9 +
              31224 * (n : ℝ) ^ 8 + 234534 * (n : ℝ) ^ 7 +
              1070705 * (n : ℝ) ^ 6 + 3185887 * (n : ℝ) ^ 5 +
              6347745 * (n : ℝ) ^ 4 + 8441595 * (n : ℝ) ^ 3 +
              7222041 * (n : ℝ) ^ 2 + 3610850 * (n : ℝ) + 805305) /
            ((n : ℝ) + 1) ^ 3 by
          dsimp [A, B, C]
          norm_num [positiveMatrix, Matrix.cons_val_two]
          field_simp
          ring]
      positivity
    ·
      dsimp [A, B, C]
      norm_num [positiveMatrix, Matrix.cons_val_two]
      ring_nf
      positivity
    · rw [show A + B * (3 / 2) + C * 2 =
          (128 * (n : ℝ) ^ 10 + 5442 * (n : ℝ) ^ 9 +
              77495 * (n : ℝ) ^ 8 + 575780 * (n : ℝ) ^ 7 +
              2606623 * (n : ℝ) ^ 6 + 7695660 * (n : ℝ) ^ 5 +
              15206998 * (n : ℝ) ^ 4 + 20033148 * (n : ℝ) ^ 3 +
              16948912 * (n : ℝ) ^ 2 + 8362836 * (n : ℝ) + 1836810) /
            ((n : ℝ) + 1) ^ 3 by
          dsimp [A, B, C]
          norm_num [positiveMatrix, Matrix.cons_val_two]
          field_simp
          ring]
      positivity
  suffices 0 ≤
      2 * stepCoefficient n x y 0 -
        ((n : ℝ) + 2) ^ 3 * stepCoefficient n x y 2 by linarith
  convert hrect using 1 <;>
    simp only [stepCoefficient] <;>
    dsimp [A, B, C] <;>
    ring

/-- The exact invariant rectangle for the balanced positive denominator row. -/
theorem denominator_projective_cone (n : ℕ) :
    1 ≤ denominatorX n ∧ denominatorX n ≤ 3 / 2 ∧
      0 ≤ denominatorY n ∧ denominatorY n ≤ 2 := by
  induction n with
  | zero =>
      norm_num [denominatorX, denominatorY, positiveDenominator, denominator,
        approximants, initialMatrix, coordinateSign, Matrix.cons_val_two]
  | succ n ih =>
      rcases ih with ⟨hx₀, hx₁, hy₀, hy₁⟩
      have hA0 :
          0 < stepCoefficient n (denominatorX n) (denominatorY n) 0 :=
        stepCoefficient_pos n _ _ (le_trans (by norm_num) hx₀) hy₀ 0
      have hA2 :
          0 < stepCoefficient n (denominatorX n) (denominatorY n) 2 :=
        stepCoefficient_pos n _ _ (le_trans (by norm_num) hx₀) hy₀ 2
      have hq : (0 : ℝ) < positiveDenominator n 0 := by
        exact_mod_cast positiveDenominator_pos n 0
      have hq' : (0 : ℝ) < positiveDenominator (n + 1) 0 := by
        exact_mod_cast positiveDenominator_pos (n + 1) 0
      have hxlow := step_x_lower n (denominatorX n) (denominatorY n)
        hx₀ hx₁ hy₀ hy₁
      have hxhigh := step_x_upper n (denominatorX n) (denominatorY n)
        hx₀ hx₁ hy₀ hy₁
      have hyhigh := step_y_upper n (denominatorX n) (denominatorY n)
        hx₀ hx₁ hy₀ hy₁
      have hX :
          denominatorX (n + 1) =
            ((n : ℝ) + 2) *
                stepCoefficient n (denominatorX n) (denominatorY n) 1 /
              stepCoefficient n (denominatorX n) (denominatorY n) 0 := by
        rw [denominatorX, positiveDenominator_succ_factor,
          positiveDenominator_succ_factor]
        field_simp
        norm_num [Nat.cast_add, Nat.cast_one]
        ring
      have hY :
          denominatorY (n + 1) =
            ((n : ℝ) + 2) ^ 3 *
                stepCoefficient n (denominatorX n) (denominatorY n) 2 /
              stepCoefficient n (denominatorX n) (denominatorY n) 0 := by
        rw [denominatorY, positiveDenominator_succ_factor,
          positiveDenominator_succ_factor]
        field_simp
        norm_num [Nat.cast_add, Nat.cast_one]
        ring
      rw [hX, hY]
      constructor
      · apply (le_div_iff₀ hA0).2
        nlinarith
      constructor
      · apply (div_le_iff₀ hA0).2
        nlinarith
      constructor
      · positivity
      · apply (div_le_iff₀ hA0).2
        nlinarith

private theorem rowZeroCorner_nonneg (n : ℕ) (j : Fin 3) :
    0 ≤ 5 * (positiveMatrix (n : ℤ) 0 j : ℝ) -
      (3 / 2 : ℝ) / ((n : ℝ) + 1) *
        (positiveMatrix (n : ℤ) 1 j : ℝ) -
      2 / ((n : ℝ) + 1) ^ 3 *
        (positiveMatrix (n : ℤ) 2 j : ℝ) := by
  fin_cases j
  · change 0 ≤
      5 * (positiveMatrix (n : ℤ) 0 0 : ℝ) -
        (3 / 2 : ℝ) / ((n : ℝ) + 1) *
          (positiveMatrix (n : ℤ) 1 0 : ℝ) -
        2 / ((n : ℝ) + 1) ^ 3 *
          (positiveMatrix (n : ℤ) 2 0 : ℝ)
    rw [show
        5 * (positiveMatrix (n : ℤ) 0 0 : ℝ) -
            (3 / 2 : ℝ) / ((n : ℝ) + 1) *
              (positiveMatrix (n : ℤ) 1 0 : ℝ) -
            2 / ((n : ℝ) + 1) ^ 3 *
              (positiveMatrix (n : ℤ) 2 0 : ℝ) =
          ((n : ℝ) + 3) ^ 2 * (2 * (n : ℝ) + 5) *
            (536 * (n : ℝ) ^ 7 + 7010 * (n : ℝ) ^ 6 +
              37549 * (n : ℝ) ^ 5 + 105798 * (n : ℝ) ^ 4 +
              166809 * (n : ℝ) ^ 3 + 143271 * (n : ℝ) ^ 2 +
              58534 * (n : ℝ) + 7149) / ((n : ℝ) + 1) ^ 3 by
        norm_num [positiveMatrix, Matrix.cons_val_two]
        field_simp
        ring]
    positivity
  · change 0 ≤
      5 * (positiveMatrix (n : ℤ) 0 1 : ℝ) -
        (3 / 2 : ℝ) / ((n : ℝ) + 1) *
          (positiveMatrix (n : ℤ) 1 1 : ℝ) -
        2 / ((n : ℝ) + 1) ^ 3 *
          (positiveMatrix (n : ℤ) 2 1 : ℝ)
    rw [show
        5 * (positiveMatrix (n : ℤ) 0 1 : ℝ) -
            (3 / 2 : ℝ) / ((n : ℝ) + 1) *
              (positiveMatrix (n : ℤ) 1 1 : ℝ) -
            2 / ((n : ℝ) + 1) ^ 3 *
              (positiveMatrix (n : ℤ) 2 1 : ℝ) =
          (3024 * (n : ℝ) ^ 9 + 58152 * (n : ℝ) ^ 8 +
              484644 * (n : ℝ) ^ 7 + 2285508 * (n : ℝ) ^ 6 +
              6674247 * (n : ℝ) ^ 5 + 12393322 * (n : ℝ) ^ 4 +
              14413315 * (n : ℝ) ^ 3 + 9859510 * (n : ℝ) ^ 2 +
              3399360 * (n : ℝ) + 372798) /
            (2 * ((n : ℝ) + 1) ^ 3) by
        norm_num [positiveMatrix, Matrix.cons_val_two]
        field_simp
        ring]
    positivity
  · change 0 ≤
      5 * (positiveMatrix (n : ℤ) 0 2 : ℝ) -
        (3 / 2 : ℝ) / ((n : ℝ) + 1) *
          (positiveMatrix (n : ℤ) 1 2 : ℝ) -
        2 / ((n : ℝ) + 1) ^ 3 *
          (positiveMatrix (n : ℤ) 2 2 : ℝ)
    rw [show
        5 * (positiveMatrix (n : ℤ) 0 2 : ℝ) -
            (3 / 2 : ℝ) / ((n : ℝ) + 1) *
              (positiveMatrix (n : ℤ) 1 2 : ℝ) -
            2 / ((n : ℝ) + 1) ^ 3 *
              (positiveMatrix (n : ℤ) 2 2 : ℝ) =
          (1888 * (n : ℝ) ^ 7 + 24466 * (n : ℝ) ^ 6 +
              129719 * (n : ℝ) ^ 5 + 361132 * (n : ℝ) ^ 4 +
              560879 * (n : ℝ) ^ 3 + 471566 * (n : ℝ) ^ 2 +
              185268 * (n : ℝ) + 19890) / ((n : ℝ) + 1) ^ 3 by
        norm_num [positiveMatrix, Matrix.cons_val_two]
        field_simp
        ring]
    positivity

private theorem rowOneCorner_nonneg (n : ℕ) (j : Fin 3) :
    0 ≤ 5 / ((n : ℝ) + 1) *
        (positiveMatrix (n : ℤ) 1 j : ℝ) -
      (positiveMatrix (n : ℤ) 0 j : ℝ) -
      2 / ((n : ℝ) + 1) ^ 3 *
        (positiveMatrix (n : ℤ) 2 j : ℝ) := by
  fin_cases j
  · change 0 ≤
      5 / ((n : ℝ) + 1) * (positiveMatrix (n : ℤ) 1 0 : ℝ) -
        (positiveMatrix (n : ℤ) 0 0 : ℝ) -
        2 / ((n : ℝ) + 1) ^ 3 *
          (positiveMatrix (n : ℤ) 2 0 : ℝ)
    rw [show
        5 / ((n : ℝ) + 1) * (positiveMatrix (n : ℤ) 1 0 : ℝ) -
            (positiveMatrix (n : ℤ) 0 0 : ℝ) -
            2 / ((n : ℝ) + 1) ^ 3 *
              (positiveMatrix (n : ℤ) 2 0 : ℝ) =
          ((n : ℝ) + 3) ^ 2 * (2 * (n : ℝ) + 5) *
            (344 * (n : ℝ) ^ 7 + 4780 * (n : ℝ) ^ 6 +
              27622 * (n : ℝ) ^ 5 + 85655 * (n : ℝ) ^ 4 +
              152784 * (n : ℝ) ^ 3 + 154802 * (n : ℝ) ^ 2 +
              80698 * (n : ℝ) + 16011) / ((n : ℝ) + 1) ^ 3 by
        norm_num [positiveMatrix, Matrix.cons_val_two]
        field_simp
        ring]
    positivity
  · change 0 ≤
      5 / ((n : ℝ) + 1) * (positiveMatrix (n : ℤ) 1 1 : ℝ) -
        (positiveMatrix (n : ℤ) 0 1 : ℝ) -
        2 / ((n : ℝ) + 1) ^ 3 *
          (positiveMatrix (n : ℤ) 2 1 : ℝ)
    rw [show
        5 / ((n : ℝ) + 1) * (positiveMatrix (n : ℤ) 1 1 : ℝ) -
            (positiveMatrix (n : ℤ) 0 1 : ℝ) -
            2 / ((n : ℝ) + 1) ^ 3 *
              (positiveMatrix (n : ℤ) 2 1 : ℝ) =
          (976 * (n : ℝ) ^ 9 + 19480 * (n : ℝ) ^ 8 +
              169804 * (n : ℝ) ^ 7 + 845942 * (n : ℝ) ^ 6 +
              2644162 * (n : ℝ) ^ 5 + 5349328 * (n : ℝ) ^ 4 +
              6951524 * (n : ℝ) ^ 3 + 5530419 * (n : ℝ) ^ 2 +
              2396694 * (n : ℝ) + 415341) / ((n : ℝ) + 1) ^ 3 by
        norm_num [positiveMatrix, Matrix.cons_val_two]
        field_simp
        ring]
    positivity
  · change 0 ≤
      5 / ((n : ℝ) + 1) * (positiveMatrix (n : ℤ) 1 2 : ℝ) -
        (positiveMatrix (n : ℤ) 0 2 : ℝ) -
        2 / ((n : ℝ) + 1) ^ 3 *
          (positiveMatrix (n : ℤ) 2 2 : ℝ)
    rw [show
        5 / ((n : ℝ) + 1) * (positiveMatrix (n : ℤ) 1 2 : ℝ) -
            (positiveMatrix (n : ℤ) 0 2 : ℝ) -
            2 / ((n : ℝ) + 1) ^ 3 *
              (positiveMatrix (n : ℤ) 2 2 : ℝ) =
          2 * (544 * (n : ℝ) ^ 7 + 7468 * (n : ℝ) ^ 6 +
              42622 * (n : ℝ) ^ 5 + 130366 * (n : ℝ) ^ 4 +
              228682 * (n : ℝ) ^ 3 + 226473 * (n : ℝ) ^ 2 +
              113884 * (n : ℝ) + 21075) / ((n : ℝ) + 1) ^ 3 by
        norm_num [positiveMatrix, Matrix.cons_val_two]
        field_simp
        ring]
    positivity

private theorem step_rowZero_six_le (n : ℕ) (x y : ℝ)
    (hx : x ≤ 3 / 2) (hy : y ≤ 2) (j : Fin 3) :
    stepCoefficient n x y j ≤
      6 * (positiveMatrix (n : ℤ) 0 j : ℝ) := by
  have hn : (0 : ℝ) < (n : ℝ) + 1 := by positivity
  have h1 : (0 : ℝ) < positiveMatrix (n : ℤ) 1 j := by
    exact_mod_cast positiveMatrix_pos n 1 j
  have h2 : (0 : ℝ) < positiveMatrix (n : ℤ) 2 j := by
    exact_mod_cast positiveMatrix_pos n 2 j
  have hx' :
      x / ((n : ℝ) + 1) * (positiveMatrix (n : ℤ) 1 j : ℝ) ≤
        (3 / 2 : ℝ) / ((n : ℝ) + 1) *
          (positiveMatrix (n : ℤ) 1 j : ℝ) := by
    exact mul_le_mul_of_nonneg_right
      (div_le_div_of_nonneg_right hx hn.le) h1.le
  have hy' :
      y / ((n : ℝ) + 1) ^ 3 * (positiveMatrix (n : ℤ) 2 j : ℝ) ≤
        2 / ((n : ℝ) + 1) ^ 3 *
          (positiveMatrix (n : ℤ) 2 j : ℝ) := by
    exact mul_le_mul_of_nonneg_right
      (div_le_div_of_nonneg_right hy (pow_nonneg hn.le 3)) h2.le
  have hc := rowZeroCorner_nonneg n j
  simp only [stepCoefficient]
  linarith

private theorem step_rowOne_six_le (n : ℕ) (x y : ℝ)
    (hx : 1 ≤ x) (hy : y ≤ 2) (j : Fin 3) :
    stepCoefficient n x y j ≤
      6 * (x / ((n : ℝ) + 1) *
        (positiveMatrix (n : ℤ) 1 j : ℝ)) := by
  have hn : (0 : ℝ) < (n : ℝ) + 1 := by positivity
  have h1 : (0 : ℝ) < positiveMatrix (n : ℤ) 1 j := by
    exact_mod_cast positiveMatrix_pos n 1 j
  have h2 : (0 : ℝ) < positiveMatrix (n : ℤ) 2 j := by
    exact_mod_cast positiveMatrix_pos n 2 j
  have hx' :
      1 / ((n : ℝ) + 1) * (positiveMatrix (n : ℤ) 1 j : ℝ) ≤
        x / ((n : ℝ) + 1) * (positiveMatrix (n : ℤ) 1 j : ℝ) := by
    exact mul_le_mul_of_nonneg_right
      (div_le_div_of_nonneg_right hx hn.le) h1.le
  have hy' :
      y / ((n : ℝ) + 1) ^ 3 * (positiveMatrix (n : ℤ) 2 j : ℝ) ≤
        2 / ((n : ℝ) + 1) ^ 3 *
          (positiveMatrix (n : ℤ) 2 j : ℝ) := by
    exact mul_le_mul_of_nonneg_right
      (div_le_div_of_nonneg_right hy (pow_nonneg hn.le 3)) h2.le
  have hx5 :
      5 / ((n : ℝ) + 1) * (positiveMatrix (n : ℤ) 1 j : ℝ) ≤
        5 * (x / ((n : ℝ) + 1) *
          (positiveMatrix (n : ℤ) 1 j : ℝ)) := by
    calc
      5 / ((n : ℝ) + 1) * (positiveMatrix (n : ℤ) 1 j : ℝ) =
          5 * (1 / ((n : ℝ) + 1) *
            (positiveMatrix (n : ℤ) 1 j : ℝ)) := by ring
      _ ≤ _ := mul_le_mul_of_nonneg_left hx' (by norm_num)
  have hc := rowOneCorner_nonneg n j
  simp only [stepCoefficient]
  linarith

private def columnWeightZero (n : ℕ) (j : Fin 3) : ℝ :=
  (positiveMatrix (n : ℤ) 0 j : ℝ) /
    stepCoefficient n (denominatorX n) (denominatorY n) j

private def columnWeightOne (n : ℕ) (j : Fin 3) : ℝ :=
  (denominatorX n / ((n : ℝ) + 1) *
      (positiveMatrix (n : ℤ) 1 j : ℝ)) /
    stepCoefficient n (denominatorX n) (denominatorY n) j

private def columnWeightTwo (n : ℕ) (j : Fin 3) : ℝ :=
  (denominatorY n / ((n : ℝ) + 1) ^ 3 *
      (positiveMatrix (n : ℤ) 2 j : ℝ)) /
    stepCoefficient n (denominatorX n) (denominatorY n) j

private theorem columnWeightZero_nonneg (n : ℕ) (j : Fin 3) :
    0 ≤ columnWeightZero n j := by
  have hcone := denominator_projective_cone n
  have hA := stepCoefficient_pos n (denominatorX n) (denominatorY n)
    (le_trans (by norm_num) hcone.1) hcone.2.2.1 j
  have hP : (0 : ℝ) < positiveMatrix (n : ℤ) 0 j := by
    exact_mod_cast positiveMatrix_pos n 0 j
  exact div_nonneg hP.le hA.le

private theorem columnWeightOne_nonneg (n : ℕ) (j : Fin 3) :
    0 ≤ columnWeightOne n j := by
  have hcone := denominator_projective_cone n
  have hn : (0 : ℝ) < (n : ℝ) + 1 := by positivity
  have hA := stepCoefficient_pos n (denominatorX n) (denominatorY n)
    (le_trans (by norm_num) hcone.1) hcone.2.2.1 j
  have hP : (0 : ℝ) < positiveMatrix (n : ℤ) 1 j := by
    exact_mod_cast positiveMatrix_pos n 1 j
  exact div_nonneg (mul_nonneg (div_nonneg (le_trans (by norm_num) hcone.1)
    hn.le) hP.le) hA.le

private theorem columnWeightTwo_nonneg (n : ℕ) (j : Fin 3) :
    0 ≤ columnWeightTwo n j := by
  have hcone := denominator_projective_cone n
  have hn : (0 : ℝ) < (n : ℝ) + 1 := by positivity
  have hA := stepCoefficient_pos n (denominatorX n) (denominatorY n)
    (le_trans (by norm_num) hcone.1) hcone.2.2.1 j
  have hP : (0 : ℝ) < positiveMatrix (n : ℤ) 2 j := by
    exact_mod_cast positiveMatrix_pos n 2 j
  exact div_nonneg (mul_nonneg (div_nonneg hcone.2.2.1
    (pow_nonneg hn.le 3)) hP.le) hA.le

private theorem columnWeights_sum (n : ℕ) (j : Fin 3) :
    columnWeightZero n j + columnWeightOne n j + columnWeightTwo n j = 1 := by
  have hcone := denominator_projective_cone n
  have hA : stepCoefficient n (denominatorX n) (denominatorY n) j ≠ 0 :=
    (stepCoefficient_pos n _ _
      (le_trans (by norm_num) hcone.1) hcone.2.2.1 j).ne'
  simp only [columnWeightZero, columnWeightOne, columnWeightTwo]
  rw [← add_div, ← add_div]
  exact (div_eq_one_iff_eq hA).2 rfl

private theorem one_sixth_le_columnWeightZero (n : ℕ) (j : Fin 3) :
    (1 / 6 : ℝ) ≤ columnWeightZero n j := by
  rcases denominator_projective_cone n with ⟨hx₀, hx₁, hy₀, hy₁⟩
  have hA := stepCoefficient_pos n (denominatorX n) (denominatorY n)
    (le_trans (by norm_num) hx₀) hy₀ j
  rw [columnWeightZero]
  apply (le_div_iff₀ hA).2
  have h := step_rowZero_six_le n (denominatorX n) (denominatorY n)
    hx₁ hy₁ j
  nlinarith

private theorem one_sixth_le_columnWeightOne (n : ℕ) (j : Fin 3) :
    (1 / 6 : ℝ) ≤ columnWeightOne n j := by
  rcases denominator_projective_cone n with ⟨hx₀, hx₁, hy₀, hy₁⟩
  have hA := stepCoefficient_pos n (denominatorX n) (denominatorY n)
    (le_trans (by norm_num) hx₀) hy₀ j
  rw [columnWeightOne]
  apply (le_div_iff₀ hA).2
  have h := step_rowOne_six_le n (denominatorX n) (denominatorY n)
    hx₀ hy₁ j
  nlinarith

private theorem positiveRatio_succ_weighted (n : ℕ) (j : Fin 3) :
    positiveRatio (n + 1) j =
      columnWeightZero n j * positiveRatio n 0 +
        columnWeightOne n j * positiveRatio n 1 +
        columnWeightTwo n j * positiveRatio n 2 := by
  have hq : (positiveDenominator n 0 : ℝ) ≠ 0 := by
    exact_mod_cast (positiveDenominator_pos n 0).ne'
  have hn : (0 : ℝ) < (n : ℝ) + 1 := by positivity
  have hqOne :
      (positiveDenominator n 1 : ℝ) =
        (positiveDenominator n 0 : ℝ) *
          (denominatorX n / ((n : ℝ) + 1)) := by
    simp only [denominatorX]
    field_simp
  have hqTwo :
      (positiveDenominator n 2 : ℝ) =
        (positiveDenominator n 0 : ℝ) *
          (denominatorY n / ((n : ℝ) + 1) ^ 3) := by
    simp only [denominatorY]
    field_simp
  rw [positiveRatio_succ_convex, Fin.sum_univ_three, Fin.sum_univ_three]
  rw [hqOne, hqTwo]
  simp only [columnWeightZero, columnWeightOne, columnWeightTwo,
    stepCoefficient]
  field_simp

private theorem refined_weighted_bounds
    (w₀ w₁ w₂ r₀ r₁ r₂ l u : ℝ)
    (hsum : w₀ + w₁ + w₂ = 1)
    (hw₀ : 1 / 6 ≤ w₀) (hw₁ : 1 / 6 ≤ w₁) (hw₂ : 0 ≤ w₂)
    (hl₀ : l ≤ r₀) (hl₁ : l ≤ r₁) (hl₂ : l ≤ r₂)
    (hu₀ : r₀ ≤ u) (hu₁ : r₁ ≤ u) (hu₂ : r₂ ≤ u) :
    (r₀ + r₁) / 6 + (2 / 3) * l ≤
        w₀ * r₀ + w₁ * r₁ + w₂ * r₂ ∧
      w₀ * r₀ + w₁ * r₁ + w₂ * r₂ ≤
        (r₀ + r₁) / 6 + (2 / 3) * u := by
  have hres₀ : 0 ≤ w₀ - 1 / 6 := sub_nonneg.mpr hw₀
  have hres₁ : 0 ≤ w₁ - 1 / 6 := sub_nonneg.mpr hw₁
  have hlo₀ := mul_nonneg hres₀ (sub_nonneg.mpr hl₀)
  have hlo₁ := mul_nonneg hres₁ (sub_nonneg.mpr hl₁)
  have hlo₂ := mul_nonneg hw₂ (sub_nonneg.mpr hl₂)
  have hup₀ := mul_nonneg hres₀ (sub_nonneg.mpr hu₀)
  have hup₁ := mul_nonneg hres₁ (sub_nonneg.mpr hu₁)
  have hup₂ := mul_nonneg hw₂ (sub_nonneg.mpr hu₂)
  have hw₂eq : w₂ = 1 - w₀ - w₁ := by linarith
  have hloSum :
      0 ≤ (w₀ - 1 / 6) * (r₀ - l) +
        (w₁ - 1 / 6) * (r₁ - l) + w₂ * (r₂ - l) :=
    add_nonneg (add_nonneg hlo₀ hlo₁) hlo₂
  have hupSum :
      0 ≤ (w₀ - 1 / 6) * (u - r₀) +
        (w₁ - 1 / 6) * (u - r₁) + w₂ * (u - r₂) :=
    add_nonneg (add_nonneg hup₀ hup₁) hup₂
  have hloEq :
      w₀ * r₀ + w₁ * r₁ + w₂ * r₂ -
          ((r₀ + r₁) / 6 + (2 / 3) * l) =
        (w₀ - 1 / 6) * (r₀ - l) +
          (w₁ - 1 / 6) * (r₁ - l) + w₂ * (r₂ - l) := by
    rw [hw₂eq]
    ring
  have hupEq :
      (r₀ + r₁) / 6 + (2 / 3) * u -
          (w₀ * r₀ + w₁ * r₁ + w₂ * r₂) =
        (w₀ - 1 / 6) * (u - r₀) +
          (w₁ - 1 / 6) * (u - r₁) + w₂ * (u - r₂) := by
    rw [hw₂eq]
    ring
  constructor <;> linarith

private theorem positiveRatio_succ_refined
    (n : ℕ) (j : Fin 3) (l u : ℝ)
    (h : ∀ i : Fin 3, l ≤ positiveRatio n i ∧ positiveRatio n i ≤ u) :
    (positiveRatio n 0 + positiveRatio n 1) / 6 + (2 / 3) * l ≤
        positiveRatio (n + 1) j ∧
      positiveRatio (n + 1) j ≤
        (positiveRatio n 0 + positiveRatio n 1) / 6 + (2 / 3) * u := by
  rw [positiveRatio_succ_weighted]
  apply refined_weighted_bounds
  · exact columnWeights_sum n j
  · exact one_sixth_le_columnWeightZero n j
  · exact one_sixth_le_columnWeightOne n j
  · exact columnWeightTwo_nonneg n j
  · exact (h 0).1
  · exact (h 1).1
  · exact (h 2).1
  · exact (h 0).2
  · exact (h 1).2
  · exact (h 2).2

/-- A nested lower envelope for all three positive ratios. -/
def lowerEnvelope : ℕ → ℝ
  | 0 => (8240 : ℝ) / 9000
  | n + 1 =>
      (positiveRatio n 0 + positiveRatio n 1) / 6 +
        (2 / 3) * lowerEnvelope n

/-- A nested upper envelope for all three positive ratios. -/
def upperEnvelope : ℕ → ℝ
  | 0 => (30921 : ℝ) / 33750
  | n + 1 =>
      (positiveRatio n 0 + positiveRatio n 1) / 6 +
        (2 / 3) * upperEnvelope n

@[simp] theorem lowerEnvelope_zero :
    lowerEnvelope 0 = (8240 : ℝ) / 9000 := rfl

@[simp] theorem upperEnvelope_zero :
    upperEnvelope 0 = (30921 : ℝ) / 33750 := rfl

@[simp] theorem lowerEnvelope_succ (n : ℕ) :
    lowerEnvelope (n + 1) =
      (positiveRatio n 0 + positiveRatio n 1) / 6 +
        (2 / 3) * lowerEnvelope n := rfl

@[simp] theorem upperEnvelope_succ (n : ℕ) :
    upperEnvelope (n + 1) =
      (positiveRatio n 0 + positiveRatio n 1) / 6 +
        (2 / 3) * upperEnvelope n := rfl

/-- Every ratio lies in the recursively contracting envelope. -/
theorem positiveRatio_envelope (n : ℕ) (j : Fin 3) :
    lowerEnvelope n ≤ positiveRatio n j ∧
      positiveRatio n j ≤ upperEnvelope n := by
  induction n generalizing j with
  | zero =>
      have h := challengeRatio_bounds 0 j
      rw [challengeRatio_eq_positiveRatio] at h
      simpa using h
  | succ n ih =>
      exact positiveRatio_succ_refined n j (lowerEnvelope n) (upperEnvelope n)
        (fun i => ih i)

theorem lowerEnvelope_monotone : Monotone lowerEnvelope := by
  apply monotone_nat_of_le_succ
  intro n
  have h₀ := (positiveRatio_envelope n 0).1
  have h₁ := (positiveRatio_envelope n 1).1
  rw [lowerEnvelope_succ]
  nlinarith

theorem upperEnvelope_antitone : Antitone upperEnvelope := by
  apply antitone_nat_of_succ_le
  intro n
  have h₀ := (positiveRatio_envelope n 0).2
  have h₁ := (positiveRatio_envelope n 1).2
  rw [upperEnvelope_succ]
  nlinarith

/-- Width of the common enclosing interval. -/
def envelopeGap (n : ℕ) : ℝ := upperEnvelope n - lowerEnvelope n

@[simp] theorem envelopeGap_succ (n : ℕ) :
    envelopeGap (n + 1) = (2 / 3 : ℝ) * envelopeGap n := by
  simp only [envelopeGap, upperEnvelope_succ, lowerEnvelope_succ]
  ring

theorem envelopeGap_eq (n : ℕ) :
    envelopeGap n = (2 / 3 : ℝ) ^ n * envelopeGap 0 := by
  induction n with
  | zero => simp
  | succ n ih =>
      rw [envelopeGap_succ, ih, pow_succ]
      ring

theorem envelopeGap_tendsto_zero :
    Filter.Tendsto envelopeGap Filter.atTop (nhds 0) := by
  have hp :
      Filter.Tendsto (fun n : ℕ => (2 / 3 : ℝ) ^ n)
        Filter.atTop (nhds 0) :=
    tendsto_pow_atTop_nhds_zero_of_lt_one (by norm_num) (by norm_num)
  have hmul := hp.mul_const (envelopeGap 0)
  simpa only [zero_mul, ← envelopeGap_eq] using hmul

/-- The common projective limit selected by the positive cocycle. -/
def commonLimit : ℝ := ⨆ n : ℕ, lowerEnvelope n

private theorem lowerEnvelope_bddAbove :
    BddAbove (Set.range lowerEnvelope) := by
  refine ⟨(30921 : ℝ) / 33750, ?_⟩
  rintro _ ⟨n, rfl⟩
  have h₁ := (positiveRatio_envelope n 0).1
  have h₂ := (challengeRatio_bounds n 0).2
  rw [challengeRatio_eq_positiveRatio] at h₂
  exact h₁.trans h₂

theorem lowerEnvelope_tendsto_common :
    Filter.Tendsto lowerEnvelope Filter.atTop (nhds commonLimit) := by
  exact tendsto_atTop_ciSup lowerEnvelope_monotone lowerEnvelope_bddAbove

theorem upperEnvelope_tendsto_common :
    Filter.Tendsto upperEnvelope Filter.atTop (nhds commonLimit) := by
  have hsum := lowerEnvelope_tendsto_common.add envelopeGap_tendsto_zero
  have heq : (fun n => lowerEnvelope n + envelopeGap n) = upperEnvelope := by
    funext n
    simp only [envelopeGap]
    ring
  rw [heq] at hsum
  simpa using hsum

/-- All three positive ratios converge unconditionally to the same limit. -/
theorem positiveRatio_tendsto_common (j : Fin 3) :
    Filter.Tendsto (fun n : ℕ => positiveRatio n j)
      Filter.atTop (nhds commonLimit) := by
  exact tendsto_of_tendsto_of_tendsto_of_le_of_le
    lowerEnvelope_tendsto_common upperEnvelope_tendsto_common
    (fun n => (positiveRatio_envelope n j).1)
    (fun n => (positiveRatio_envelope n j).2)

/-- Unconditional common convergence of the three literal challenge ratios. -/
theorem challengeRatio_tendsto_common (j : Fin 3) :
    Filter.Tendsto (fun n : ℕ => challengeRatio n j)
      Filter.atTop (nhds commonLimit) := by
  simpa only [challengeRatio_eq_positiveRatio] using
    positiveRatio_tendsto_common j

/-- The original Problem 2.5 is now equivalent to one scalar connection
identity, rather than two vector-valued dominant-mode assumptions. -/
theorem problem25Claim_iff_commonLimit_eq_catalan :
    Problem25Claim ↔ commonLimit = catalanConstant := by
  constructor
  · intro h
    exact tendsto_nhds_unique (challengeRatio_tendsto_common 0) (h 0)
  · intro h
    intro j
    simpa only [h] using challengeRatio_tendsto_common j

end RamanujanChallenge.P25
