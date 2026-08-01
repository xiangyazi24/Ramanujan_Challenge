import RamanujanChallenge.Problem22Concentration

open Filter Topology Real
open scoped BigOperators

noncomputable section

namespace RamanujanChallenge.P22.Q6502

/-- Uniform harmonic-minus-log remainder using `log (m+1)`, valid also at `m=0`. -/
theorem abs_harmonic_sub_log_succ_sub_gamma_le (m : ℕ) :
    |((harmonic m : ℚ) : ℝ) - Real.log ((m : ℝ) + 1) -
        Real.eulerMascheroniConstant| ≤
      1 / ((m : ℝ) + 1) := by
  have hlo := Real.eulerMascheroniSeq_lt_eulerMascheroniConstant m
  have hhi := Real.eulerMascheroniConstant_lt_eulerMascheroniSeq' (m + 1)
  simp only [Real.eulerMascheroniSeq, Real.eulerMascheroniSeq',
    Nat.succ_ne_zero, if_false, harmonic_succ, Rat.cast_add,
    Rat.cast_inv, Rat.cast_natCast, Rat.cast_one, Nat.cast_add,
    Nat.cast_one] at hlo hhi
  have hneg :
      ((harmonic m : ℚ) : ℝ) - Real.log ((m : ℝ) + 1) -
          Real.eulerMascheroniConstant ≤ 0 := by
    linarith
  rw [abs_of_nonpos hneg]
  have hbound :
      -(((harmonic m : ℚ) : ℝ) - Real.log ((m : ℝ) + 1) -
          Real.eulerMascheroniConstant) <
        (((m : ℝ) + 1)⁻¹) := by
    linarith
  simpa only [one_div] using hbound.le

/-- Lipschitz bound for the logarithm on a positive interval. -/
theorem abs_log_sub_log_le_abs_sub_div_min
    {a b : ℝ} (ha : 0 < a) (hb : 0 < b) :
    |Real.log a - Real.log b| ≤ |a - b| / min a b := by
  rcases le_total a b with hab | hba
  · have hlog : Real.log a ≤ Real.log b := Real.log_le_log ha hab
    rw [abs_of_nonpos (sub_nonpos.mpr hlog), abs_of_nonpos (sub_nonpos.mpr hab),
      min_eq_left hab]
    have hratio : 0 < b / a := div_pos hb ha
    have hle := Real.log_le_sub_one_of_pos hratio
    rw [Real.log_div hb.ne' ha.ne'] at hle
    have ha0 : a ≠ 0 := ha.ne'
    field_simp [ha0] at hle ⊢
    nlinarith
  · have hlog : Real.log b ≤ Real.log a := Real.log_le_log hb hba
    rw [abs_of_nonneg (sub_nonneg.mpr hlog), abs_of_nonneg (sub_nonneg.mpr hba),
      min_eq_right hba]
    have hratio : 0 < a / b := div_pos ha hb
    have hle := Real.log_le_sub_one_of_pos hratio
    rw [Real.log_div ha.ne' hb.ne'] at hle
    have hb0 : b ≠ 0 := hb.ne'
    field_simp [hb0] at hle ⊢
    nlinarith

#print axioms abs_harmonic_sub_log_succ_sub_gamma_le
#print axioms abs_log_sub_log_le_abs_sub_div_min

end RamanujanChallenge.P22.Q6502
