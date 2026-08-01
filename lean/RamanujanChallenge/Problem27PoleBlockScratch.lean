import RamanujanChallenge.Problem27BarnesNormalization
import Mathlib.Analysis.Calculus.ParametricIntegral
import Mathlib.Analysis.PSeries
import Mathlib.Analysis.SpecialFunctions.Gamma.Basic
import Mathlib.Analysis.SpecialFunctions.PolynomialExp
import Mathlib.MeasureTheory.Integral.DominatedConvergence
import Mathlib.MeasureTheory.Integral.ExpDecay
import Mathlib.MeasureTheory.Integral.Prod
import Mathlib.NumberTheory.ZetaValues
import Mathlib.Topology.Algebra.InfiniteSum.NatInt

open Filter Set MeasureTheory Topology
open scoped Topology

noncomputable section

namespace RamanujanChallenge.P27

private theorem integrableOn_pow_mul_exp_neg_mul27
    (n : ℕ) {a : ℝ} (ha : 0 < a) :
    IntegrableOn (fun t : ℝ => t ^ n * Real.exp (-(a * t))) (Ioi 0) := by
  have hbase :
      IntegrableOn (fun x : ℝ => Real.exp (-x) * x ^ n) (Ioi 0) := by
    convert Real.GammaIntegral_convergent
      (s := (n : ℝ) + 1) (by positivity) using 1
    ext x
    rw [show (n : ℝ) + 1 - 1 = n by ring, Real.rpow_natCast]
  have hbase0 :
      IntegrableOn (fun x : ℝ => Real.exp (-x) * x ^ n) (Ioi (a * 0)) := by
    simpa only [mul_zero] using hbase
  have hscaled :
      IntegrableOn
        (fun t : ℝ => Real.exp (-(a * t)) * (a * t) ^ n) (Ioi 0) := by
    exact (integrableOn_Ioi_comp_mul_left_iff
      (fun x : ℝ => Real.exp (-x) * x ^ n) 0 ha).2 hbase0
  have hc := hscaled.const_mul ((a ^ n)⁻¹)
  refine hc.congr ?_
  filter_upwards [ae_restrict_mem measurableSet_Ioi] with t ht
  rw [mul_pow]
  field_simp [ha.ne']

private theorem integrableOn_cpow_mul_cexp_neg_mul27
    (n : ℕ) {z : ℂ} (hz : 0 < z.re) :
    IntegrableOn
      (fun t : ℝ => (t : ℂ) ^ n * Complex.exp (-z * (t : ℂ)))
      (Ioi 0) := by
  have hreal := integrableOn_pow_mul_exp_neg_mul27 n hz
  apply hreal.mono'
  · apply Continuous.aestronglyMeasurable
    fun_prop
  · filter_upwards [ae_restrict_mem measurableSet_Ioi] with t ht
    rw [norm_mul, norm_pow, Complex.norm_real, Real.norm_eq_abs,
      abs_of_pos ht, Complex.norm_exp]
    simp

private theorem integral_cexp_neg_mul_Ioi27
    {z : ℂ} (hz : 0 < z.re) :
    (∫ t : ℝ in Ioi 0, Complex.exp (-z * (t : ℂ))) = z⁻¹ := by
  have h := integral_exp_mul_complex_Ioi (a := -z) (by simpa using hz) 0
  simpa [div_eq_mul_inv] using h

private theorem tendsto_self_mul_exp_neg_mul27
    {a : ℝ} (ha : 0 < a) :
    Tendsto (fun t : ℝ => t * Real.exp (-(a * t))) atTop (𝓝 0) := by
  have h := (Real.tendsto_pow_mul_exp_neg_atTop_nhds_zero 1).comp
    (tendsto_id.const_mul_atTop' ha)
  have hc := h.const_mul a⁻¹
  convert hc using 1
  · funext t
    simp only [Function.comp_apply, id_eq, pow_one]
    field_simp [ha.ne']
  · simp

private theorem integral_self_mul_cexp_neg_mul_Ioi27
    {z : ℂ} (hz : 0 < z.re) :
    (∫ t : ℝ in Ioi 0, (t : ℂ) * Complex.exp (-z * (t : ℂ))) =
      (z ^ 2)⁻¹ := by
  have hz0 : z ≠ 0 := by
    intro h
    rw [h] at hz
    norm_num at hz
  let u : ℝ → ℂ := fun t => t
  let v : ℝ → ℂ := fun t => -Complex.exp (-z * (t : ℂ)) / z
  let u' : ℝ → ℂ := fun _ => 1
  let v' : ℝ → ℂ := fun t => Complex.exp (-z * (t : ℂ))
  have hu (t : ℝ) : HasDerivAt u 1 t := by
    exact (hasDerivAt_id t).ofReal_comp
  have hv (t : ℝ) : HasDerivAt v (Complex.exp (-z * (t : ℂ))) t := by
    have hlin :
        HasDerivAt (fun x : ℝ => -z * (x : ℂ)) (-z) t := by
      convert (hasDerivAt_id t).ofReal_comp.const_mul (-z) using 1 <;> simp
    convert hlin.cexp.neg.div_const z using 1
    field_simp [hz0]
  have huv' : IntegrableOn (u * v') (Ioi 0) := by
    simpa [u, v'] using integrableOn_cpow_mul_cexp_neg_mul27 1 hz
  have hu'v : IntegrableOn (u' * v) (Ioi 0) := by
    have h := (integrableOn_cpow_mul_cexp_neg_mul27 0 hz).neg.div_const z
    refine h.congr ?_
    filter_upwards [ae_restrict_mem measurableSet_Ioi] with t ht
    simp [u', v]
  have hzero : Tendsto (u * v) (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    have hc : ContinuousAt (u * v) 0 := by
      dsimp [u, v]
      fun_prop
    simpa [u, v] using hc.continuousWithinAt.tendsto
  have hinfty : Tendsto (u * v) atTop (𝓝 0) := by
    rw [tendsto_zero_iff_norm_tendsto_zero]
    have hreal :
        Tendsto (fun t : ℝ =>
          t * Real.exp (-(z.re * t)) / ‖z‖) atTop (𝓝 0) := by
      simpa using (tendsto_self_mul_exp_neg_mul27 hz).div_const ‖z‖
    refine hreal.congr' ?_
    filter_upwards [eventually_ge_atTop (0 : ℝ)] with t ht
    dsimp [u, v]
    rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg ht,
      norm_div, norm_neg, Complex.norm_exp]
    simp
    ring
  have hibp := MeasureTheory.integral_Ioi_mul_deriv_eq_deriv_mul
    (a := (0 : ℝ)) (u := u) (v := v) (u' := u') (v' := v')
    (fun t _ => hu t) (fun t _ => hv t) huv' hu'v hzero hinfty
  have hvint :
      (∫ t : ℝ in Ioi 0, v t) = -(z⁻¹ * z⁻¹) := by
    calc
      (∫ t : ℝ in Ioi 0, v t) =
          ∫ t : ℝ in Ioi 0,
            (-z⁻¹) * Complex.exp (-z * (t : ℂ)) := by
              apply integral_congr_ae
              filter_upwards with t
              simp [v, div_eq_mul_inv]
              ring
      _ = (-z⁻¹) * ∫ t : ℝ in Ioi 0,
            Complex.exp (-z * (t : ℂ)) :=
        MeasureTheory.integral_const_mul _ _
      _ = -(z⁻¹ * z⁻¹) := by rw [integral_cexp_neg_mul_Ioi27 hz]; ring
  have hibp' :
      (∫ t : ℝ in Ioi 0,
        (t : ℂ) * Complex.exp (-z * (t : ℂ))) =
        -(∫ t : ℝ in Ioi 0, v t) := by
    simpa [u, v', u'] using hibp
  calc
    (∫ t : ℝ in Ioi 0,
        (t : ℂ) * Complex.exp (-z * (t : ℂ))) =
        -(∫ t : ℝ in Ioi 0, v t) := hibp'
    _ = -(-(z⁻¹ * z⁻¹)) := by rw [hvint]
    _ = (z ^ 2)⁻¹ := by field_simp [hz0]

private theorem integral_laplace_poleBlock27
    {z : ℂ} (hz : 0 < z.re) :
    (∫ t : ℝ in Ioi 0,
      ((1 + t / 2 : ℝ) : ℂ) * Complex.exp (-z * (t : ℂ))) =
      z⁻¹ + (z ^ 2)⁻¹ / 2 := by
  have h0 : IntegrableOn
      (fun t : ℝ => Complex.exp (-z * (t : ℂ))) (Ioi 0) := by
    simpa using integrableOn_cpow_mul_cexp_neg_mul27 0 hz
  have h1 : IntegrableOn
      (fun t : ℝ => (t : ℂ) * Complex.exp (-z * (t : ℂ))) (Ioi 0) := by
    simpa using integrableOn_cpow_mul_cexp_neg_mul27 1 hz
  calc
    (∫ t : ℝ in Ioi 0,
        ((1 + t / 2 : ℝ) : ℂ) * Complex.exp (-z * (t : ℂ))) =
      ∫ t : ℝ in Ioi 0,
        Complex.exp (-z * (t : ℂ)) +
          (1 / 2 : ℂ) * ((t : ℂ) * Complex.exp (-z * (t : ℂ))) := by
            apply integral_congr_ae
            filter_upwards with t
            push_cast
            ring
    _ = (∫ t : ℝ in Ioi 0, Complex.exp (-z * (t : ℂ))) +
        ∫ t : ℝ in Ioi 0,
          (1 / 2 : ℂ) * ((t : ℂ) * Complex.exp (-z * (t : ℂ))) := by
            rw [integral_add h0 (h1.const_mul (1 / 2 : ℂ))]
    _ = (∫ t : ℝ in Ioi 0, Complex.exp (-z * (t : ℂ))) +
        (1 / 2 : ℂ) *
          ∫ t : ℝ in Ioi 0, (t : ℂ) * Complex.exp (-z * (t : ℂ)) := by
            congr 1
            exact MeasureTheory.integral_const_mul _ _
    _ = z⁻¹ + (z ^ 2)⁻¹ / 2 := by
      rw [integral_cexp_neg_mul_Ioi27 hz,
        integral_self_mul_cexp_neg_mul_Ioi27 hz]
      ring

private def poleAbscissa27 (m : ℕ) : ℝ :=
  (m : ℝ) + 1 / 2

private def polePoint27 (m : ℕ) (y : ℝ) : ℂ :=
  (poleAbscissa27 m : ℂ) + (y : ℂ) * Complex.I

@[simp] private theorem polePoint_re27 (m : ℕ) (y : ℝ) :
    (polePoint27 m y).re = poleAbscissa27 m := by
  simp [polePoint27]

private def poleBlock27 (m : ℕ) (y : ℝ) : ℂ :=
  (polePoint27 m y)⁻¹ + ((polePoint27 m y) ^ 2)⁻¹ / 2

private def laplaceWeight27 (m : ℕ) (t : ℝ) : ℝ :=
  (1 + t / 2) * Real.exp (-(poleAbscissa27 m * t))

private def sechSq27 (y : ℝ) : ℝ :=
  1 / Real.cosh (Real.pi * y) ^ 2

private def fourierKernel27 (t y : ℝ) : ℂ :=
  Complex.exp (-(y * t : ℝ) * Complex.I) /
    (Real.cosh (Real.pi * y) : ℂ) ^ 2

private theorem poleAbscissa_pos27 (m : ℕ) : 0 < poleAbscissa27 m := by
  unfold poleAbscissa27
  positivity

private theorem integrableOn_laplaceWeight27 (m : ℕ) :
    IntegrableOn (laplaceWeight27 m) (Ioi 0) := by
  have h0 := integrableOn_pow_mul_exp_neg_mul27 0 (poleAbscissa_pos27 m)
  have h1 := integrableOn_pow_mul_exp_neg_mul27 1 (poleAbscissa_pos27 m)
  have hsum := h0.add (h1.const_mul (1 / 2 : ℝ))
  refine hsum.congr ?_
  filter_upwards [ae_restrict_mem measurableSet_Ioi] with t ht
  simp [laplaceWeight27]
  ring

private theorem integrable_sechSq27 : Integrable sechSq27 := by
  apply integrable_zudilinBarnesEnvelope27.mono'
  · apply Continuous.aestronglyMeasurable
    fun_prop
  · filter_upwards with y
    let C : ℝ := Real.cosh (Real.pi * y)
    have hCpos : 0 < C := Real.cosh_pos _
    have hC1 : 1 ≤ C := Real.one_le_cosh _
    have hsqrt0 : 0 ≤ Real.sqrt C := Real.sqrt_nonneg _
    have hsqrt_sq : (Real.sqrt C) ^ 2 = C := Real.sq_sqrt hCpos.le
    have hsqrt_le : Real.sqrt C ≤ C ^ 2 := by
      nlinarith [sq_nonneg (Real.sqrt C - 1), sq_nonneg (C - 1)]
    rw [Real.norm_eq_abs, abs_of_nonneg (by unfold sechSq27; positivity)]
    unfold sechSq27 zudilinBarnesEnvelope27
    exact one_div_le_one_div_of_le (Real.sqrt_pos.2 hCpos) hsqrt_le

private theorem integral_fourierKernel27 {t : ℝ} (ht : 0 < t) :
    (∫ y : ℝ, fourierKernel27 t y) =
      ((t / (Real.pi * Real.sinh (t / 2)) : ℝ) : ℂ) := by
  let ξ : ℝ := t / (2 * Real.pi)
  have hξ : ξ ≠ 0 := by
    dsimp only [ξ]
    exact div_ne_zero ht.ne' (mul_ne_zero (by norm_num) Real.pi_ne_zero)
  have h := integral_sechSq_cexp27 ξ
  rw [if_neg hξ] at h
  calc
    (∫ y : ℝ, fourierKernel27 t y) =
        ∫ y : ℝ,
          Complex.exp (-(2 * Real.pi * y * ξ) * Complex.I) /
            (Real.cosh (Real.pi * y) : ℂ) ^ 2 := by
              apply integral_congr_ae
              filter_upwards with y
              unfold fourierKernel27
              congr 2
              push_cast
              dsimp only [ξ]
              field_simp [Real.pi_ne_zero]
              ring
    _ = ((2 * ξ / Real.sinh (Real.pi * ξ) : ℝ) : ℂ) := h
    _ = ((t / (Real.pi * Real.sinh (t / 2)) : ℝ) : ℂ) := by
      norm_cast
      dsimp only [ξ]
      rw [show Real.pi * (t / (2 * Real.pi)) = t / 2 by
        field_simp [Real.pi_ne_zero]]
      field_simp [Real.pi_ne_zero]

#print axioms integral_laplace_poleBlock27
#print axioms integral_fourierKernel27

end RamanujanChallenge.P27
