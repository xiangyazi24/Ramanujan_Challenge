import RamanujanChallenge.Problem27BarnesNormalization
import Mathlib.Analysis.SpecialFunctions.Gamma.Basic
import Mathlib.Analysis.SpecialFunctions.PolynomialExp

open Filter Set MeasureTheory Topology
open scoped Topology

noncomputable section

namespace RamanujanChallenge.P27.Q6383

theorem integrableOn_pow_mul_exp_neg_mul
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

theorem integrableOn_cpow_mul_cexp_neg_mul
    (n : ℕ) {z : ℂ} (hz : 0 < z.re) :
    IntegrableOn
      (fun t : ℝ => (t : ℂ) ^ n * Complex.exp (-z * (t : ℂ)))
      (Ioi 0) := by
  have hreal := integrableOn_pow_mul_exp_neg_mul n hz
  apply hreal.mono'
  · apply Continuous.aestronglyMeasurable
    fun_prop
  · filter_upwards [ae_restrict_mem measurableSet_Ioi] with t ht
    rw [norm_mul, norm_pow, Complex.norm_real, Real.norm_eq_abs,
      abs_of_pos ht, Complex.norm_exp]
    simp

theorem integral_cexp_neg_mul_Ioi
    {z : ℂ} (hz : 0 < z.re) :
    (∫ t : ℝ in Ioi 0, Complex.exp (-z * (t : ℂ))) = z⁻¹ := by
  have h := integral_exp_mul_complex_Ioi (a := -z) (by simpa using hz) 0
  simpa [div_eq_mul_inv] using h

theorem tendsto_self_mul_exp_neg_mul
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

theorem integral_self_mul_cexp_neg_mul_Ioi
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
    simpa [u, v'] using integrableOn_cpow_mul_cexp_neg_mul 1 hz
  have hu'v : IntegrableOn (u' * v) (Ioi 0) := by
    have h := (integrableOn_cpow_mul_cexp_neg_mul 0 hz).neg.div_const z
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
      simpa using (tendsto_self_mul_exp_neg_mul hz).div_const ‖z‖
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
      _ = -(z⁻¹ * z⁻¹) := by rw [integral_cexp_neg_mul_Ioi hz]; ring
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

theorem integral_laplace_poleBlock
    {z : ℂ} (hz : 0 < z.re) :
    (∫ t : ℝ in Ioi 0,
      ((1 + t / 2 : ℝ) : ℂ) * Complex.exp (-z * (t : ℂ))) =
      z⁻¹ + (z ^ 2)⁻¹ / 2 := by
  have h0 : IntegrableOn
      (fun t : ℝ => Complex.exp (-z * (t : ℂ))) (Ioi 0) := by
    simpa using integrableOn_cpow_mul_cexp_neg_mul 0 hz
  have h1 : IntegrableOn
      (fun t : ℝ => (t : ℂ) * Complex.exp (-z * (t : ℂ))) (Ioi 0) := by
    simpa using integrableOn_cpow_mul_cexp_neg_mul 1 hz
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
      rw [integral_cexp_neg_mul_Ioi hz,
        integral_self_mul_cexp_neg_mul_Ioi hz]
      ring

#print axioms integral_laplace_poleBlock

end RamanujanChallenge.P27.Q6383
