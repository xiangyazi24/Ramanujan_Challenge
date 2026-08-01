/-
  Problem 2.5 — finite positivity and Catalan-error checks for the
  Delannoy-basis coefficients.
-/
import RamanujanChallenge.Problem25Delannoy
import RamanujanChallenge.Problem25TightBounds

noncomputable section

namespace RamanujanChallenge.P25

open Nat Finset

set_option maxHeartbeats 0
set_option maxRecDepth 10000

/-- One step of each triangular inversion, phrased using the named coefficient
sequences.  This lets the exact evaluations below reuse earlier rows. -/
private theorem delannoyF_step (k : ℕ) :
    delannoyF k =
      (normalizedDenominator k - ∑ j : Fin k, delannoyF j * delannoyB k j) /
        delannoyB k k := by
  rw [delannoyF, delannoyCoeff]
  rfl

private theorem delannoyG_step (k : ℕ) :
    delannoyG k =
      (normalizedNumerator k - ∑ j : Fin k, delannoyG j * delannoyB k j) /
        delannoyB k k := by
  rw [delannoyG, delannoyCoeff]
  rfl

macro "delannoy_norm" : tactic =>
  `(tactic|
    (norm_num [normalizedDenominator, normalizedNumerator, denominator, numerator,
      approximants, matrixProduct, initialMatrix, challengeMatrix, m11, m12, m13,
      m21, m22, m23, m31, m32, m33, pochhammerProduct, gaugeDelta, delannoyB,
      Fin.sum_univ_succ] <;>
      norm_num [Nat.choose]))

@[simp] theorem delannoyF_one : delannoyF 1 = 5160375 / 32 := by
  rw [delannoyF_step]
  delannoy_norm

@[simp] theorem delannoyF_two : delannoyF 2 = 59423875 / 192 := by
  rw [delannoyF_step]
  delannoy_norm

@[simp] theorem delannoyF_three : delannoyF 3 = 3822850785 / 8192 := by
  rw [delannoyF_step]
  delannoy_norm

@[simp] theorem delannoyF_four : delannoyF 4 = 15766041135 / 25088 := by
  rw [delannoyF_step]
  delannoy_norm

@[simp] theorem delannoyF_five : delannoyF 5 = 39287063252875 / 49545216 := by
  rw [delannoyF_step]
  delannoy_norm

@[simp] theorem delannoyF_six : delannoyF 6 = 8945487057855125 / 9325510656 := by
  rw [delannoyF_step]
  delannoy_norm

@[simp] theorem delannoyF_seven :
    delannoyF 7 = 9277524675251255375 / 8233854959616 := by
  rw [delannoyF_step]
  delannoy_norm

@[simp] theorem delannoyF_eight :
    delannoyF 8 = 5998406532546392455 / 4631543414784 := by
  rw [delannoyF_step]
  delannoy_norm

@[simp] theorem delannoyF_nine :
    delannoyF 9 = 844607807995210534995 / 576868868685824 := by
  rw [delannoyF_step]
  delannoy_norm

@[simp] theorem delannoyF_ten :
    delannoyF 10 = 7484382923587075329187875 / 4581492555102814208 := by
  rw [delannoyF_step]
  delannoy_norm

@[simp] theorem delannoyG_one : delannoyG 1 = 463216893 / 3136 := by
  rw [delannoyG_step]
  delannoy_norm

@[simp] theorem delannoyG_two : delannoyG 2 = 432067238033 / 1524096 := by
  rw [delannoyG_step]
  delannoy_norm

@[simp] theorem delannoyG_three : delannoyG 3 = 5605465888674661 / 13113999360 := by
  rw [delannoyG_step]
  delannoy_norm

@[simp] theorem delannoyG_four :
    delannoyG 4 = 19534564802667839663 / 33936571468800 := by
  rw [delannoyG_step]
  delannoy_norm

@[simp] theorem delannoyG_five :
    delannoyG 5 = 4543258282031622755015 / 6255188853129216 := by
  rw [delannoyG_step]
  delannoy_norm

@[simp] theorem delannoyG_six :
    delannoyG 6 =
      384382993659231804297403895 / 437475398010151108608 := by
  rw [delannoyG_step]
  delannoy_norm

@[simp] theorem delannoyG_seven :
    delannoyG 7 =
      143912823236440808231730214671005 / 139441301196806838779117568 := by
  rw [delannoyG_step]
  delannoy_norm

@[simp] theorem delannoyG_eight :
    delannoyG 8 =
      22154093314780495296120468510901 / 18675174267429487336488960 := by
  rw [delannoyG_step]
  delannoy_norm

@[simp] theorem delannoyG_nine :
    delannoyG 9 =
      635315685819580344498952919076642678953 /
        473731679294867251542006064742400 := by
  rw [delannoyG_step]
  delannoy_norm

@[simp] theorem delannoyG_ten :
    delannoyG 10 =
      511797075843849345312422854111349547518443 /
        342034272450894155613328378744012800 := by
  rw [delannoyG_step]
  delannoy_norm

/-- The first eleven denominator coefficients in the Delannoy basis are positive. -/
theorem delannoyF_pos_upto_ten (k : ℕ) (hk : k ≤ 10) : 0 < delannoyF k := by
  interval_cases k <;> norm_num

/-- On the first eleven Delannoy coefficients the numerator/denominator ratio
is within `10⁻³` of Catalan's constant at `k = 0`, and within `10⁻⁵`
thereafter. -/
theorem delannoy_error_ratio_upto_ten (k : ℕ) (hk : k ≤ 10) :
    |(delannoyG k : ℝ) / (delannoyF k : ℝ) - catalanConstant| <
      if k = 0 then (1 : ℝ) / 1000 else 1 / 100000 := by
  have hlo := catalan_tight_lower
  have hup := catalan_tight_upper
  interval_cases k <;>
    rw [abs_lt] <;>
    constructor <;>
    norm_num <;>
    linarith

end RamanujanChallenge.P25

end
