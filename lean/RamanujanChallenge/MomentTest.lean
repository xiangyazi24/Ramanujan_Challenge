import Mathlib.Analysis.SpecialFunctions.Integrals.Basic

open MeasureTheory

noncomputable section

private theorem poly_integrable_01 {f : ℝ → ℝ} (hf : Continuous f) :
    IntervalIntegrable f MeasureTheory.volume 0 1 :=
  hf.intervalIntegrable 0 1

-- ∫₀¹ (10 - 18x) dx = 1
theorem moment_wq_one :
    ∫ x in (0:ℝ)..1, (10 - 18 * x) = (1 : ℝ) := by
  have h : (fun x : ℝ => (10 : ℝ) - 18 * x) = (fun x => 10 * (1 : ℝ) + (-18) * x) := by
    ext; ring
  rw [h, intervalIntegral.integral_add (poly_integrable_01 (by fun_prop))
    (poly_integrable_01 (by fun_prop)),
    intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul,
    integral_one, integral_id]
  norm_num

-- ∫₀¹ x(10-18x) dx = -1
theorem moment_wq_x :
    ∫ x in (0:ℝ)..1, x * (10 - 18 * x) = (-1 : ℝ) := by
  have h : (fun x : ℝ => x * (10 - 18 * x)) = (fun x => 10 * x + (-18) * x ^ 2) := by
    ext x; ring
  rw [h, intervalIntegral.integral_add (poly_integrable_01 (by fun_prop))
    (poly_integrable_01 (by fun_prop)),
    intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul,
    integral_id, integral_pow]
  norm_num

-- ∫₀¹ (6 - 12x) dx = 0
theorem moment_wp_zero :
    ∫ x in (0:ℝ)..1, (6 - 12 * x) = (0 : ℝ) := by
  have h : (fun x : ℝ => (6 : ℝ) - 12 * x) = (fun x => 6 * (1 : ℝ) + (-12) * x) := by
    ext; ring
  rw [h, intervalIntegral.integral_add (poly_integrable_01 (by fun_prop))
    (poly_integrable_01 (by fun_prop)),
    intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul,
    integral_one, integral_id]
  norm_num

-- ∫₀¹ x(6-12x) dx = -1
theorem moment_wp_x :
    ∫ x in (0:ℝ)..1, x * (6 - 12 * x) = (-1 : ℝ) := by
  have h : (fun x : ℝ => x * (6 - 12 * x)) = (fun x => 6 * x + (-12) * x ^ 2) := by
    ext x; ring
  rw [h, intervalIntegral.integral_add (poly_integrable_01 (by fun_prop))
    (poly_integrable_01 (by fun_prop)),
    intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul,
    integral_id, integral_pow]
  norm_num
