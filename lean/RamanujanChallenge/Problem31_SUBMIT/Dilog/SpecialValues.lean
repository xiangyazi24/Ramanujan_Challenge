import Ramanujan31.Dilog.FiveTerm

/-!
# Exact special values of the Rogers dilogarithm

The results in this file are consequences of the two functional equations
proved in `Rogers.lean` and `FiveTerm.lean`; no numerical approximation enters.
-/

namespace Real

set_option maxHeartbeats 200000

/-- If `ρ` is the root in `(0,1)` of `ρ² + ρ = 1`, then
`R(ρ) = π² / 10`. -/
theorem rogers_eq_pi_sq_div_ten_of_sq_eq_one_sub
    {ρ : ℝ} (hρ0 : 0 < ρ) (hρ1 : ρ < 1) (hρsq : ρ ^ 2 = 1 - ρ) :
    rogers ρ = Real.pi ^ 2 / 10 := by
  have hρne : ρ ≠ 0 := ne_of_gt hρ0
  have hden : 1 - ρ * ρ = ρ := by
    calc
      1 - ρ * ρ = 1 - ρ ^ 2 := by ring
      _ = ρ := by rw [hρsq]; ring
  have harg :
      ρ * (1 - ρ) / (1 - ρ * ρ) = ρ ^ 2 := by
    rw [hden]
    calc
      ρ * (1 - ρ) / ρ = 1 - ρ := by field_simp [hρne]
      _ = ρ ^ 2 := hρsq.symm
  have hfive := rogers_five_term hρ0 hρ1 hρ0 hρ1
  rw [harg] at hfive
  rw [show ρ * ρ = ρ ^ 2 by ring] at hfive
  have heuler := rogers_add_rogers_one_sub hρ0 hρ1
  rw [← hρsq] at heuler
  linarith

/-- The explicit positive golden-ratio conjugate. -/
noncomputable def goldenRatioConjugate : ℝ := (Real.sqrt 5 - 1) / 2

theorem goldenRatioConjugate_pos : 0 < goldenRatioConjugate := by
  have hsqrt : 1 < Real.sqrt 5 := by
    nlinarith [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 5),
      Real.sqrt_nonneg (5 : ℝ)]
  unfold goldenRatioConjugate
  linarith

theorem goldenRatioConjugate_lt_one : goldenRatioConjugate < 1 := by
  have hsqrt : Real.sqrt 5 < 3 := by
    nlinarith [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 5),
      Real.sqrt_nonneg (5 : ℝ)]
  unfold goldenRatioConjugate
  linarith

theorem goldenRatioConjugate_sq :
    goldenRatioConjugate ^ 2 = 1 - goldenRatioConjugate := by
  have hsqrt := Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 5)
  unfold goldenRatioConjugate
  nlinarith

/-- The Rogers value at the explicit golden-ratio conjugate. -/
theorem rogers_goldenRatioConjugate :
    rogers goldenRatioConjugate = Real.pi ^ 2 / 10 :=
  rogers_eq_pi_sq_div_ten_of_sq_eq_one_sub
    goldenRatioConjugate_pos goldenRatioConjugate_lt_one
    goldenRatioConjugate_sq

end Real
