import RamanujanChallenge.Problem27PoleBlockCore6383
import Mathlib.MeasureTheory.Integral.Prod

open Filter Set MeasureTheory Topology
open scoped Topology

noncomputable section

namespace RamanujanChallenge.P27.Q6383

def poleAbscissa (m : ℕ) : ℝ := (m : ℝ) + 1 / 2

def polePoint (m : ℕ) (y : ℝ) : ℂ :=
  (poleAbscissa m : ℂ) + (y : ℂ) * Complex.I

@[simp] theorem polePoint_re (m : ℕ) (y : ℝ) :
    (polePoint m y).re = poleAbscissa m := by
  simp [polePoint]

def poleBlock (m : ℕ) (y : ℝ) : ℂ :=
  (polePoint m y)⁻¹ + ((polePoint m y) ^ 2)⁻¹ / 2

def laplaceWeight (m : ℕ) (t : ℝ) : ℝ :=
  (1 + t / 2) * Real.exp (-(poleAbscissa m * t))

def sechSq (y : ℝ) : ℝ := 1 / Real.cosh (Real.pi * y) ^ 2

def fourierKernel (t y : ℝ) : ℂ :=
  Complex.exp (-(y * t : ℝ) * Complex.I) /
    (Real.cosh (Real.pi * y) : ℂ) ^ 2

theorem poleAbscissa_pos (m : ℕ) : 0 < poleAbscissa m := by
  unfold poleAbscissa
  positivity

theorem integrableOn_laplaceWeight (m : ℕ) :
    IntegrableOn (laplaceWeight m) (Ioi 0) := by
  have h0 := integrableOn_pow_mul_exp_neg_mul 0 (poleAbscissa_pos m)
  have h1 := integrableOn_pow_mul_exp_neg_mul 1 (poleAbscissa_pos m)
  have hsum := h0.add (h1.const_mul (1 / 2 : ℝ))
  refine hsum.congr ?_
  filter_upwards [ae_restrict_mem measurableSet_Ioi] with t ht
  simp [laplaceWeight]
  ring

theorem integrable_sechSq : Integrable sechSq := by
  apply integrable_zudilinBarnesEnvelope27.mono'
  · have hc : Continuous fun y : ℝ => Real.cosh (Real.pi * y) := by
      fun_prop
    exact (continuous_const.div (hc.pow 2)
      (fun y => pow_ne_zero 2 (ne_of_gt (Real.cosh_pos _)))).aestronglyMeasurable
  · filter_upwards with y
    let C : ℝ := Real.cosh (Real.pi * y)
    have hCpos : 0 < C := Real.cosh_pos _
    have hC1 : 1 ≤ C := Real.one_le_cosh _
    have hsqrt0 : 0 ≤ Real.sqrt C := Real.sqrt_nonneg _
    have hsqrt_sq : (Real.sqrt C) ^ 2 = C := Real.sq_sqrt hCpos.le
    have hsqrt_le : Real.sqrt C ≤ C ^ 2 := by
      nlinarith [sq_nonneg (Real.sqrt C - 1), sq_nonneg (C - 1)]
    rw [Real.norm_eq_abs, abs_of_nonneg (by unfold sechSq; positivity)]
    unfold sechSq zudilinBarnesEnvelope27
    exact one_div_le_one_div_of_le (Real.sqrt_pos.2 hCpos) hsqrt_le

theorem integral_fourierKernel {t : ℝ} (ht : 0 < t) :
    (∫ y : ℝ, fourierKernel t y) =
      ((t / (Real.pi * Real.sinh (t / 2)) : ℝ) : ℂ) := by
  let ξ : ℝ := t / (2 * Real.pi)
  have hξ : ξ ≠ 0 := by
    dsimp only [ξ]
    exact div_ne_zero ht.ne' (mul_ne_zero (by norm_num) Real.pi_ne_zero)
  have h := integral_sechSq_cexp27 ξ
  rw [if_neg hξ] at h
  calc
    (∫ y : ℝ, fourierKernel t y) =
        ∫ y : ℝ,
          Complex.exp (-(2 * Real.pi * y * ξ) * Complex.I) /
            (Real.cosh (Real.pi * y) : ℂ) ^ 2 := by
              apply integral_congr_ae
              filter_upwards with y
              unfold fourierKernel
              congr 2
              norm_cast
              dsimp only [ξ]
              field_simp [Real.pi_ne_zero]
    _ = ((2 * ξ / Real.sinh (Real.pi * ξ) : ℝ) : ℂ) := h
    _ = ((t / (Real.pi * Real.sinh (t / 2)) : ℝ) : ℂ) := by
      norm_cast
      dsimp only [ξ]
      rw [show Real.pi * (t / (2 * Real.pi)) = t / 2 by
        field_simp [Real.pi_ne_zero]]
      field_simp [Real.pi_ne_zero]

#print axioms integral_fourierKernel

end RamanujanChallenge.P27.Q6383
