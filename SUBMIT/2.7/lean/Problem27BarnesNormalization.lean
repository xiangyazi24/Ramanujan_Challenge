import RamanujanChallenge.Problem27Barnes
import Mathlib.Analysis.SpecialFunctions.Gamma.Beta
import Mathlib.Analysis.SpecialFunctions.Sigmoid
import Mathlib.MeasureTheory.Function.JacobianOneDim

/-!
# Problem 2.7: normalization of the Barnes integral

This file supplies the analytic normalization deliberately kept separate from
the midpoint estimate in `Problem27Barnes`.  The first ingredient is the
Fourier transform of the logistic density, proved by a real change of
variables to Euler's beta integral.
-/

open Filter Set MeasureTheory

noncomputable section

namespace RamanujanChallenge.P27

private def scaledSigmoid27 (y : ℝ) : ℝ :=
  Real.sigmoid (2 * Real.pi * y)

private def betaFourierIntegrand27 (ξ x : ℝ) : ℂ :=
  (2 / Real.pi : ℝ) *
    (x : ℂ) ^ (-(Complex.I * (ξ : ℂ))) *
      (1 - (x : ℂ)) ^ (Complex.I * (ξ : ℂ))

private theorem range_scaledSigmoid27 :
    Set.range scaledSigmoid27 = Ioo 0 1 := by
  rw [← Real.range_sigmoid]
  apply Set.Subset.antisymm
  · rintro _ ⟨y, rfl⟩
    exact ⟨2 * Real.pi * y, rfl⟩
  · rintro _ ⟨x, rfl⟩
    refine ⟨x / (2 * Real.pi), ?_⟩
    unfold scaledSigmoid27
    congr 1
    field_simp [Real.pi_ne_zero]

private theorem scaledSigmoid27_injective :
    Function.Injective scaledSigmoid27 := by
  intro x y h
  unfold scaledSigmoid27 at h
  have hxy := Real.sigmoid_injective h
  exact mul_left_cancel₀ (mul_ne_zero (by norm_num) Real.pi_ne_zero) hxy

private theorem hasDerivAt_scaledSigmoid27 (y : ℝ) :
    HasDerivAt scaledSigmoid27
      ((2 * Real.pi) * scaledSigmoid27 y * (1 - scaledSigmoid27 y)) y := by
  unfold scaledSigmoid27
  convert (Real.hasDerivAt_sigmoid (2 * Real.pi * y)).comp y
    ((hasDerivAt_id y).const_mul (2 * Real.pi)) using 1 <;> ring

private theorem scaledSigmoid27_mul_one_sub (y : ℝ) :
    4 * scaledSigmoid27 y * (1 - scaledSigmoid27 y) =
      1 / Real.cosh (Real.pi * y) ^ 2 := by
  have hcosh : Real.cosh (Real.pi * y) =
      Real.exp (Real.pi * y) *
        (1 + Real.exp (-((2 * Real.pi) * y))) / 2 := by
    rw [Real.cosh_eq]
    have hprod : Real.exp (Real.pi * y) *
        Real.exp (-((2 * Real.pi) * y)) = Real.exp (-(Real.pi * y)) := by
      rw [← Real.exp_add]
      congr 1
      ring
    rw [mul_add, mul_one, hprod]
  simp only [scaledSigmoid27, Real.sigmoid_def]
  rw [hcosh]
  field_simp [Real.exp_ne_zero]
  simp only [add_sub_cancel_left]
  rw [pow_two, ← Real.exp_add, mul_assoc, ← Real.exp_add]
  rw [show -(2 * Real.pi * y) + (Real.pi * y + Real.pi * y) = 0 by ring,
    Real.exp_zero]
  norm_num

private theorem scaledSigmoid27_logit (y : ℝ) :
    Real.log (scaledSigmoid27 y) -
        Real.log (1 - scaledSigmoid27 y) = 2 * Real.pi * y := by
  unfold scaledSigmoid27
  rw [← Real.log_div
    (ne_of_gt (Real.sigmoid_pos (2 * Real.pi * y)))
    (ne_of_gt (sub_pos.mpr (Real.sigmoid_lt_one (2 * Real.pi * y))))]
  simp only [scaledSigmoid27, Real.sigmoid_def]
  have he : 0 < Real.exp (-((2 * Real.pi) * y)) := Real.exp_pos _
  have hsum : 0 < 1 + Real.exp (-((2 * Real.pi) * y)) := by positivity
  rw [show (1 + Real.exp (-((2 * Real.pi) * y)))⁻¹ /
      (1 - (1 + Real.exp (-((2 * Real.pi) * y)))⁻¹) =
      Real.exp ((2 * Real.pi) * y) by
        field_simp
        simp only [add_sub_cancel_left]
        rw [← Real.exp_add]
        ring_nf
        rw [Real.exp_zero]]
  exact Real.log_exp _

private theorem betaFourierIntegrand_scaledSigmoid27 (ξ y : ℝ) :
    betaFourierIntegrand27 ξ (scaledSigmoid27 y) =
      (2 / Real.pi : ℝ) *
        Complex.exp (-(2 * Real.pi * y * ξ) * Complex.I) := by
  have hx0 : 0 < scaledSigmoid27 y := Real.sigmoid_pos _
  have hx1 : 0 < 1 - scaledSigmoid27 y := sub_pos.mpr (Real.sigmoid_lt_one _)
  unfold betaFourierIntegrand27
  rw [Complex.cpow_def_of_ne_zero (Complex.ofReal_ne_zero.mpr hx0.ne'),
    show 1 - (scaledSigmoid27 y : ℂ) =
      ((1 - scaledSigmoid27 y : ℝ) : ℂ) by exact (Complex.ofReal_sub _ _).symm,
    Complex.cpow_def_of_ne_zero (Complex.ofReal_ne_zero.mpr hx1.ne')]
  rw [← Complex.ofReal_log hx0.le, ← Complex.ofReal_log hx1.le,
    mul_assoc, ← Complex.exp_add]
  congr 2
  have hlog :
      ((Real.log (scaledSigmoid27 y) : ℝ) : ℂ) -
          (Real.log (1 - scaledSigmoid27 y) : ℂ) =
        (2 * Real.pi * y : ℝ) := by
    exact_mod_cast scaledSigmoid27_logit y
  calc
    ((Real.log (scaledSigmoid27 y) : ℂ) * -(Complex.I * ξ) +
        (Real.log (1 - scaledSigmoid27 y) : ℂ) * (Complex.I * ξ)) =
        -(Complex.I * ξ) *
          ((Real.log (scaledSigmoid27 y) : ℂ) -
            (Real.log (1 - scaledSigmoid27 y) : ℂ)) := by ring
    _ = -(Complex.I * ξ) * (2 * Real.pi * y : ℝ) := by rw [hlog]
    _ = -(2 * (Real.pi : ℂ) * (y : ℂ) * (ξ : ℂ)) * Complex.I := by
      push_cast
      ring

/-- Fourier transform of `sech² (π y)`, in beta-integral form. -/
theorem integral_sechSq_cexp_eq_beta27 (ξ : ℝ) :
    (∫ y : ℝ,
      Complex.exp (-(2 * Real.pi * y * ξ) * Complex.I) /
        (Real.cosh (Real.pi * y) : ℂ) ^ 2) =
      (2 / Real.pi : ℝ) *
        Complex.betaIntegral (1 - Complex.I * ξ) (1 + Complex.I * ξ) := by
  let f : ℝ → ℝ := scaledSigmoid27
  let f' : ℝ → ℝ := fun y ↦
    (2 * Real.pi) * scaledSigmoid27 y * (1 - scaledSigmoid27 y)
  let g : ℝ → ℂ := betaFourierIntegrand27 ξ
  have hchange := integral_image_eq_integral_abs_deriv_smul
    (s := Set.univ) (f := f) (f' := f') MeasurableSet.univ
    (fun y _ ↦ (hasDerivAt_scaledSigmoid27 y).hasDerivWithinAt)
    scaledSigmoid27_injective.injOn g
  have hfimage : f '' Set.univ = Ioo (0 : ℝ) 1 := by
    simpa only [Set.image_univ] using range_scaledSigmoid27
  rw [hfimage] at hchange
  have hleft :
      (∫ x in Ioo (0 : ℝ) 1, g x) =
        (2 / Real.pi : ℝ) *
          Complex.betaIntegral (1 - Complex.I * ξ) (1 + Complex.I * ξ) := by
    rw [← MeasureTheory.integral_Ioc_eq_integral_Ioo,
      ← intervalIntegral.integral_of_le (by norm_num : (0 : ℝ) ≤ 1)]
    rw [Complex.betaIntegral]
    calc
      (∫ x : ℝ in 0..1, g x) =
          ∫ x : ℝ in 0..1,
            ((2 / Real.pi : ℝ) : ℂ) *
              ((x : ℂ) ^ (1 - Complex.I * ξ - 1) *
                (1 - (x : ℂ)) ^ (1 + Complex.I * ξ - 1)) := by
            apply intervalIntegral.integral_congr
            intro x hx
            unfold g betaFourierIntegrand27
            ring_nf
      _ = ((2 / Real.pi : ℝ) : ℂ) *
          ∫ x : ℝ in 0..1,
            (x : ℂ) ^ (1 - Complex.I * ξ - 1) *
              (1 - (x : ℂ)) ^ (1 + Complex.I * ξ - 1) :=
        intervalIntegral.integral_const_mul _ _
  rw [hleft] at hchange
  simp only [Measure.restrict_univ] at hchange
  rw [hchange]
  apply integral_congr_ae
  filter_upwards with y
  have hfpos : 0 < f' y := by
    dsimp only [f']
    exact mul_pos
      (mul_pos (mul_pos (by norm_num) Real.pi_pos) (Real.sigmoid_pos _))
      (sub_pos.mpr (Real.sigmoid_lt_one _))
  rw [abs_of_pos hfpos]
  dsimp only [f, f', g]
  rw [betaFourierIntegrand_scaledSigmoid27]
  have hw :
      1 / (Real.cosh (Real.pi * y) : ℂ) ^ 2 =
        4 * (scaledSigmoid27 y : ℂ) * (1 - scaledSigmoid27 y) := by
    exact_mod_cast (scaledSigmoid27_mul_one_sub y).symm
  rw [div_eq_mul_inv, ← one_div, hw]
  change
    Complex.exp (-(2 * (Real.pi : ℂ) * y * ξ) * Complex.I) *
        (4 * (scaledSigmoid27 y : ℂ) * (1 - scaledSigmoid27 y)) =
      (((2 * Real.pi) * scaledSigmoid27 y *
          (1 - scaledSigmoid27 y) : ℝ) : ℂ) *
        (((2 / Real.pi : ℝ) : ℂ) *
          Complex.exp (-(2 * (Real.pi : ℂ) * y * ξ) * Complex.I))
  push_cast
  field_simp [Real.pi_ne_zero]
  ring

private theorem betaIntegral_one_sub_I_mul_one_add_I_mul27
    (ξ : ℝ) :
    Complex.betaIntegral (1 - Complex.I * ξ) (1 + Complex.I * ξ) =
      if ξ = 0 then 1 else
        ((Real.pi * ξ / Real.sinh (Real.pi * ξ) : ℝ) : ℂ) := by
  by_cases hξ : ξ = 0
  · subst ξ
    simp [Complex.betaIntegral_eval_one_right]
  · rw [if_neg hξ]
    have hIξ : Complex.I * (ξ : ℂ) ≠ 0 :=
      mul_ne_zero Complex.I_ne_zero (Complex.ofReal_ne_zero.mpr hξ)
    rw [Complex.betaIntegral_eq_Gamma_mul_div]
    · rw [show (1 - Complex.I * ξ) + (1 + Complex.I * ξ) = 2 by ring]
      norm_num [Complex.Gamma_nat_eq_factorial]
      rw [show 1 + Complex.I * ξ = Complex.I * ξ + 1 by ring,
        Complex.Gamma_add_one _ hIξ]
      have hreflect := Complex.Gamma_mul_Gamma_one_sub (Complex.I * (ξ : ℂ))
      rw [show 1 - Complex.I * (ξ : ℂ) = 1 - Complex.I * ξ by rfl] at hreflect
      calc
        Complex.Gamma (1 - Complex.I * ξ) *
              (Complex.I * ξ * Complex.Gamma (Complex.I * ξ)) =
            (Complex.I * ξ) *
              (Complex.Gamma (Complex.I * ξ) *
                Complex.Gamma (1 - Complex.I * ξ)) := by ring
        _ = (Complex.I * ξ) *
              ((Real.pi : ℂ) /
                Complex.sin ((Real.pi : ℂ) * (Complex.I * ξ))) := by
              rw [hreflect]
        _ = (Real.pi : ℂ) * (ξ : ℂ) /
              Complex.sinh ((Real.pi : ℂ) * (ξ : ℂ)) := by
              rw [show (Real.pi : ℂ) * (Complex.I * ξ) =
                ((Real.pi * ξ : ℝ) : ℂ) * Complex.I by push_cast; ring,
                Complex.sin_mul_I]
              field_simp [Complex.I_ne_zero]
              rw [show (ξ : ℂ) * (Real.pi : ℂ) =
                ((Real.pi * ξ : ℝ) : ℂ) by push_cast; ring]
    · norm_num
    · norm_num

/-- Closed form of the Fourier transform of `sech² (π y)`. -/
theorem integral_sechSq_cexp27 (ξ : ℝ) :
    (∫ y : ℝ,
      Complex.exp (-(2 * Real.pi * y * ξ) * Complex.I) /
        (Real.cosh (Real.pi * y) : ℂ) ^ 2) =
      if ξ = 0 then (2 / Real.pi : ℝ) else
        (2 * ξ / Real.sinh (Real.pi * ξ) : ℝ) := by
  rw [integral_sechSq_cexp_eq_beta27,
    betaIntegral_one_sub_I_mul_one_add_I_mul27]
  split_ifs with hξ
  · push_cast
    ring
  · push_cast
    field_simp [Real.pi_ne_zero]

end RamanujanChallenge.P27
