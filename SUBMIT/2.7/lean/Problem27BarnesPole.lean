import RamanujanChallenge.Problem27BarnesMoments
import Mathlib.Analysis.SpecialFunctions.ImproperIntegrals
import Mathlib.MeasureTheory.Integral.Prod
import Mathlib.MeasureTheory.Integral.DominatedConvergence
import Mathlib.NumberTheory.ZetaValues
import Mathlib.Analysis.PSeries

open Filter Set MeasureTheory Asymptotics Topology

noncomputable section

namespace RamanujanChallenge.P27

private theorem integrableOn_cexp_neg_mul27 {z : ℂ} (hz : 0 < z.re) :
    IntegrableOn (fun t : ℝ => Complex.exp (-z * t)) (Ioi 0) := by
  simpa only [neg_mul] using
    integrableOn_exp_mul_complex_Ioi (a := -z) (by simpa using hz) 0

private theorem integrableOn_mul_cexp_neg27 {z : ℂ} (hz : 0 < z.re) :
    IntegrableOn (fun t : ℝ => (t : ℂ) * Complex.exp (-z * t)) (Ioi 0) := by
  have hr := integrableOn_rpow_mul_exp_neg_mul_rpow
    (s := (1 : ℝ)) (p := (1 : ℝ)) (b := z.re)
    (by norm_num) (by norm_num) hz
  have hr' : IntegrableOn (fun t : ℝ => t * Real.exp (-z.re * t)) (Ioi 0) := by
    apply hr.congr_fun _ measurableSet_Ioi
    intro t ht
    simp only [Real.rpow_one]
  apply hr'.mono'
  · fun_prop
  · filter_upwards [ae_restrict_mem measurableSet_Ioi] with t ht
    have ht0 : 0 < t := ht
    rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_pos ht0,
      Complex.norm_exp]
    simp only [Complex.mul_re, Complex.neg_re, Complex.neg_im, Complex.ofReal_re,
      Complex.ofReal_im, mul_zero, sub_zero]
    exact le_rfl

private theorem integral_cexp_neg_mul_Ioi27 {z : ℂ} (hz : 0 < z.re) :
    (∫ t : ℝ in Ioi 0, Complex.exp (-z * t)) = 1 / z := by
  rw [show (fun t : ℝ => Complex.exp (-z * t)) =
      fun t : ℝ => Complex.exp ((-z) * t) by rfl,
    integral_exp_mul_complex_Ioi (by simpa using hz)]
  simp only [mul_zero, Complex.ofReal_zero, Complex.exp_zero]
  have hz0 : z ≠ 0 := by
    intro h
    rw [h] at hz
    norm_num at hz
  field_simp [hz0]

private theorem tendsto_mul_cexp_neg_atTop27 {z : ℂ} (hz : 0 < z.re) :
    Tendsto (fun t : ℝ => (t : ℂ) * Complex.exp (-z * t)) atTop (𝓝 0) := by
  rw [Metric.tendsto_atTop]
  intro ε hε
  have hreal : Tendsto (fun t : ℝ => t * Real.exp (-z.re * t)) atTop (𝓝 0) := by
    have hbase := Real.tendsto_pow_mul_exp_neg_atTop_nhds_zero 1
    have hscale : Tendsto (fun t : ℝ => z.re * t) atTop atTop :=
      tendsto_id.const_mul_atTop' hz
    have h := hbase.comp hscale
    have hc : Tendsto (fun _ : ℝ => 1 / z.re) atTop (𝓝 (1 / z.re)) :=
      tendsto_const_nhds
    convert (hc.mul h) using 1
    · funext t
      simp only [Function.comp_apply, pow_one]
      field_simp [ne_of_gt hz]
    · simp
  rw [Metric.tendsto_atTop] at hreal
  obtain ⟨a, ha⟩ := hreal ε hε
  refine ⟨max a 0, fun b hb => ?_⟩
  have hba : a ≤ b := le_trans (le_max_left _ _) hb
  have hb0 : 0 ≤ b := le_trans (le_max_right _ _) hb
  specialize ha b hba
  simpa only [dist_zero_right, norm_mul, Complex.norm_real, Real.norm_eq_abs,
    abs_of_nonneg hb0, Complex.norm_exp, Complex.mul_re, Complex.neg_re,
    Complex.neg_im, Complex.ofReal_re, Complex.ofReal_im, mul_zero, sub_zero,
    Real.norm_eq_abs, abs_of_nonneg (Real.exp_pos _).le] using ha

private theorem integral_mul_cexp_neg_mul_Ioi27 {z : ℂ} (hz : 0 < z.re) :
    (∫ t : ℝ in Ioi 0, (t : ℂ) * Complex.exp (-z * t)) = 1 / z ^ 2 := by
  have hz0 : z ≠ 0 := by
    intro h
    rw [h] at hz
    norm_num at hz
  have hvint : IntegrableOn
      (fun t : ℝ => -Complex.exp (-z * t) / z) (Ioi 0) :=
    (integrableOn_cexp_neg_mul27 hz).neg.div_const z
  have hparts := MeasureTheory.integral_Ioi_mul_deriv_eq_deriv_mul
    (a := (0 : ℝ))
    (u := fun t : ℝ => (t : ℂ))
    (u' := fun _ : ℝ => (1 : ℂ))
    (v := fun t : ℝ => -Complex.exp (-z * t) / z)
    (v' := fun t : ℝ => Complex.exp (-z * t))
    (a' := (0 : ℂ)) (b' := (0 : ℂ))
    (fun t _ => (hasDerivAt_id t).ofReal_comp)
    (fun t _ => by
      have he := ((hasDerivAt_id t).ofReal_comp.const_mul (-z)).cexp
      convert he.neg.div_const z using 1 <;>
        simp only [id_eq, one_mul] <;> field_simp [hz0] <;> ring <;> norm_num)
    (integrableOn_mul_cexp_neg27 hz)
    (by
      exact IntegrableOn.congr_fun hvint
        (fun t _ => by simp only [Pi.mul_apply, one_mul]) measurableSet_Ioi)
    (by
      have hc : ContinuousAt
          (fun t : ℝ => (t : ℂ) * (-Complex.exp (-z * t) / z)) 0 := by
        fun_prop
      simpa using hc.tendsto.mono_left inf_le_left)
    (by
      have ht := (tendsto_mul_cexp_neg_atTop27 hz).neg.div_const z
      convert ht using 1
      · funext t
        simp only [Pi.mul_apply, div_eq_mul_inv]
        ring
      · simp)
  simp only [Pi.mul_apply, one_mul, zero_sub] at hparts
  calc
    (∫ t : ℝ in Ioi 0, (t : ℂ) * Complex.exp (-z * t)) =
        -(∫ x : ℝ in Ioi 0, -Complex.exp (-z * x) / z) := by
      simpa using hparts
    _ = -(-(∫ x : ℝ in Ioi 0, Complex.exp (-z * x)) / z) := by
      have hv :
          (∫ x : ℝ in Ioi 0, -Complex.exp (-z * x) * z⁻¹) =
            (-(∫ x : ℝ in Ioi 0, Complex.exp (-z * x))) * z⁻¹ := by
        calc
          (∫ x : ℝ in Ioi 0, -Complex.exp (-z * x) * z⁻¹) =
              (∫ x : ℝ in Ioi 0, -Complex.exp (-z * x)) * z⁻¹ :=
            MeasureTheory.integral_mul_const z⁻¹ _
          _ = (-(∫ x : ℝ in Ioi 0, Complex.exp (-z * x))) * z⁻¹ := by
            rw [MeasureTheory.integral_neg]
      simpa only [div_eq_mul_inv] using congrArg Neg.neg hv
    _ = 1 / z ^ 2 := by
      rw [integral_cexp_neg_mul_Ioi27 hz]
      field_simp [hz0]

private theorem integrableOn_exp_mul_t_pow27 (c : ℝ) (hc : 0 < c) (n : ℕ) :
    IntegrableOn (fun t : ℝ => Real.exp (-c * t) * t ^ n) (Ioi 0) := by
  have h := integrableOn_rpow_mul_exp_neg_mul_rpow
    (s := (n : ℝ)) (p := (1 : ℝ)) (b := c)
    (by have hn : (0 : ℝ) ≤ n := Nat.cast_nonneg n; linarith)
    (by norm_num) hc
  apply h.congr_fun _ measurableSet_Ioi
  intro t ht
  change t ^ (n : ℝ) * Real.exp (-c * t ^ (1 : ℝ)) =
    Real.exp (-c * t) * t ^ n
  rw [Real.rpow_natCast, Real.rpow_one]
  ring

private theorem integral_exp_mul_t_add_sq27 (c : ℝ) (hc : 0 < c) :
    (∫ t : ℝ in Ioi 0, Real.exp (-c * t) * (t + t ^ 2 / 2)) =
      1 / c ^ 2 + 1 / c ^ 3 := by
  have h1 := Real.integral_rpow_mul_exp_neg_mul_Ioi
    (a := (2 : ℝ)) (r := c) (by norm_num) hc
  have h2 := Real.integral_rpow_mul_exp_neg_mul_Ioi
    (a := (3 : ℝ)) (r := c) (by norm_num) hc
  have hi1 := integrableOn_exp_mul_t_pow27 c hc 1
  have hi2 := integrableOn_exp_mul_t_pow27 c hc 2
  have hi1' : IntegrableOn
      (fun t : ℝ => Real.exp (-c * t) * t) (Ioi 0) := by
    simpa only [pow_one] using hi1
  have hi2' : IntegrableOn
      (fun t : ℝ => (Real.exp (-c * t) * t ^ 2) / 2) (Ioi 0) :=
    hi2.div_const 2
  have hA :
      (∫ t : ℝ in Ioi 0, Real.exp (-c * t) * t) =
        ∫ t : ℝ in Ioi 0, t ^ ((2 : ℝ) - 1) * Real.exp (-(c * t)) := by
    apply setIntegral_congr_fun measurableSet_Ioi
    intro t ht
    change Real.exp (-c * t) * t =
      t ^ ((2 : ℝ) - 1) * Real.exp (-(c * t))
    rw [show t ^ ((2 : ℝ) - 1) = t by norm_num]
    ring
  have hB :
      (∫ t : ℝ in Ioi 0, Real.exp (-c * t) * t ^ 2) =
        ∫ t : ℝ in Ioi 0, t ^ ((3 : ℝ) - 1) * Real.exp (-(c * t)) := by
    apply setIntegral_congr_fun measurableSet_Ioi
    intro t ht
    change Real.exp (-c * t) * t ^ 2 =
      t ^ ((3 : ℝ) - 1) * Real.exp (-(c * t))
    rw [show t ^ ((3 : ℝ) - 1) = t ^ 2 by norm_num [Real.rpow_two]]
    ring
  calc
    (∫ t : ℝ in Ioi 0, Real.exp (-c * t) * (t + t ^ 2 / 2)) =
        ∫ t : ℝ in Ioi 0,
          Real.exp (-c * t) * t + (Real.exp (-c * t) * t ^ 2) / 2 := by
      apply setIntegral_congr_fun measurableSet_Ioi
      intro t ht
      ring
    _ =
        (∫ t : ℝ in Ioi 0, Real.exp (-c * t) * t) +
          (∫ t : ℝ in Ioi 0, (Real.exp (-c * t) * t ^ 2) / 2) := by
      exact MeasureTheory.integral_add hi1' hi2'
    _ = (∫ t : ℝ in Ioi 0, t ^ ((2 : ℝ) - 1) * Real.exp (-(c * t))) +
          (∫ t : ℝ in Ioi 0, t ^ ((3 : ℝ) - 1) * Real.exp (-(c * t))) / 2 := by
      rw [MeasureTheory.integral_div]
      rw [hA, hB]
    _ = 1 / c ^ 2 + 1 / c ^ 3 := by
      rw [h1, h2]
      norm_num [Real.Gamma_nat_eq_factorial, Real.rpow_two,
        Real.rpow_natCast]

private def barnesTailTerm27 (m j : ℕ) (t : ℝ) : ℝ :=
  Real.exp (-(((m + 1 + j : ℕ) : ℝ) * t)) * (t + t ^ 2 / 2)

private theorem integrableOn_barnesTailTerm27 (m j : ℕ) :
    IntegrableOn (barnesTailTerm27 m j) (Ioi 0) := by
  have hc : 0 < ((m + 1 + j : ℕ) : ℝ) := by positivity
  have h1 := integrableOn_exp_mul_t_pow27 ((m + 1 + j : ℕ) : ℝ) hc 1
  have h2 := integrableOn_exp_mul_t_pow27 ((m + 1 + j : ℕ) : ℝ) hc 2
  apply (h1.add (h2.div_const 2)).congr_fun _ measurableSet_Ioi
  intro t ht
  simp only [barnesTailTerm27, pow_one, Pi.add_apply, Pi.div_apply]
  ring

private theorem integral_barnesTailTerm27 (m j : ℕ) :
    (∫ t : ℝ in Ioi 0, barnesTailTerm27 m j t) =
      1 / (((m + 1 + j : ℕ) : ℝ) ^ 2) +
        1 / (((m + 1 + j : ℕ) : ℝ) ^ 3) := by
  unfold barnesTailTerm27
  simpa only [neg_mul] using
    integral_exp_mul_t_add_sq27 ((m + 1 + j : ℕ) : ℝ) (by positivity)

private theorem hasSum_barnesTailTerm27 (m : ℕ) {t : ℝ} (ht : 0 < t) :
    HasSum (fun j : ℕ => barnesTailTerm27 m j t)
      (Real.exp (-(((m + 1 : ℕ) : ℝ) * t)) * (t + t ^ 2 / 2) /
        (1 - Real.exp (-t))) := by
  have hq : |Real.exp (-t)| < 1 := by
    rw [abs_of_pos (Real.exp_pos _), Real.exp_lt_one_iff]
    linarith
  have hgeom := (hasSum_geometric_of_abs_lt_one hq).mul_left
    (Real.exp (-(((m + 1 : ℕ) : ℝ) * t)) * (t + t ^ 2 / 2))
  exact hgeom.congr_fun fun j => by
    simp only [barnesTailTerm27]
    rw [← Real.exp_nat_mul]
    have hexp :
        Real.exp (-(((m + 1 + j : ℕ) : ℝ) * t)) =
          Real.exp (-(((m + 1 : ℕ) : ℝ) * t)) *
            Real.exp ((j : ℝ) * -t) := by
      rw [← Real.exp_add]
      congr 1
      push_cast
      ring
    rw [hexp]
    ring

private theorem summable_barnesTailValues27 (m : ℕ) :
    Summable (fun j : ℕ =>
      1 / (((m + 1 + j : ℕ) : ℝ) ^ 2) +
        1 / (((m + 1 + j : ℕ) : ℝ) ^ 3)) := by
  have h2base : Summable (fun n : ℕ => 1 / (n : ℝ) ^ 2) :=
    hasSum_zeta_two.summable
  have h3base : Summable (fun n : ℕ => 1 / (n : ℝ) ^ 3) :=
    Real.summable_one_div_nat_pow.mpr (by norm_num)
  have h2 := (summable_nat_add_iff (m + 1)).2 h2base
  have h3 := (summable_nat_add_iff (m + 1)).2 h3base
  apply (h2.add h3).congr
  intro j
  change
    1 / (((j + (m + 1) : ℕ) : ℝ) ^ 2) +
        1 / (((j + (m + 1) : ℕ) : ℝ) ^ 3) =
      1 / (((m + 1 + j : ℕ) : ℝ) ^ 2) +
        1 / (((m + 1 + j : ℕ) : ℝ) ^ 3)
  push_cast
  ring

private theorem integral_barnesTail27 (m : ℕ) :
    (∫ t : ℝ in Ioi 0,
      Real.exp (-(((m + 1 : ℕ) : ℝ) * t)) * (t + t ^ 2 / 2) /
        (1 - Real.exp (-t))) =
      ∑' j : ℕ, (1 / (((m + 1 + j : ℕ) : ℝ) ^ 2) +
        1 / (((m + 1 + j : ℕ) : ℝ) ^ 3)) := by
  have hInt (j : ℕ) : Integrable (barnesTailTerm27 m j)
      (volume.restrict (Ioi 0)) := integrableOn_barnesTailTerm27 m j
  have hNorm (j : ℕ) :
      (∫ t : ℝ, ‖barnesTailTerm27 m j t‖ ∂volume.restrict (Ioi 0)) =
        1 / (((m + 1 + j : ℕ) : ℝ) ^ 2) +
          1 / (((m + 1 + j : ℕ) : ℝ) ^ 3) := by
    rw [← integral_barnesTailTerm27 m j]
    apply integral_congr_ae
    filter_upwards [ae_restrict_mem measurableSet_Ioi] with t ht
    rw [Real.norm_eq_abs, abs_of_nonneg]
    unfold barnesTailTerm27
    exact mul_nonneg (Real.exp_pos _).le
      (add_nonneg ht.le (div_nonneg (sq_nonneg t) (by norm_num)))
  have hSwap := MeasureTheory.integral_tsum_of_summable_integral_norm hInt
    ((summable_barnesTailValues27 m).congr (fun j => hNorm j |>.symm))
  calc
    (∫ t : ℝ in Ioi 0,
        Real.exp (-(((m + 1 : ℕ) : ℝ) * t)) * (t + t ^ 2 / 2) /
          (1 - Real.exp (-t))) =
        ∫ t : ℝ, (∑' j : ℕ, barnesTailTerm27 m j t)
          ∂volume.restrict (Ioi 0) := by
      apply integral_congr_ae
      filter_upwards [ae_restrict_mem measurableSet_Ioi] with t ht
      exact (hasSum_barnesTailTerm27 m ht).tsum_eq.symm
    _ = ∑' j : ℕ,
        ∫ t : ℝ, barnesTailTerm27 m j t ∂volume.restrict (Ioi 0) := hSwap.symm
    _ = ∑' j : ℕ, (1 / (((m + 1 + j : ℕ) : ℝ) ^ 2) +
        1 / (((m + 1 + j : ℕ) : ℝ) ^ 3)) := by
      apply tsum_congr
      intro j
      exact integral_barnesTailTerm27 m j

def barnesPoleBlock27 (z : ℂ) : ℂ :=
  1 / z + 1 / (2 * z ^ 2)

private theorem integral_laplacePoleBlock27 {z : ℂ} (hz : 0 < z.re) :
    (∫ t : ℝ in Ioi 0,
      Complex.exp (-z * t) * (1 + (t : ℂ) / 2)) = barnesPoleBlock27 z := by
  have h0 := integrableOn_cexp_neg_mul27 hz
  have h1 := integrableOn_mul_cexp_neg27 hz
  calc
    (∫ t : ℝ in Ioi 0,
        Complex.exp (-z * t) * (1 + (t : ℂ) / 2)) =
        (∫ t : ℝ in Ioi 0, Complex.exp (-z * t)) +
          ∫ t : ℝ in Ioi 0,
            ((t : ℂ) * Complex.exp (-z * t)) / 2 := by
      rw [← MeasureTheory.integral_add h0 (h1.div_const 2)]
      apply setIntegral_congr_fun measurableSet_Ioi
      intro t ht
      ring
    _ = 1 / z + (1 / z ^ 2) / 2 := by
      rw [integral_cexp_neg_mul_Ioi27 hz]
      congr 1
      calc
        (∫ t : ℝ in Ioi 0,
            ((t : ℂ) * Complex.exp (-z * t)) / 2) =
            (∫ t : ℝ in Ioi 0,
              (t : ℂ) * Complex.exp (-z * t)) / 2 :=
          MeasureTheory.integral_div 2 _
        _ = (1 / z ^ 2) / 2 := by
          rw [integral_mul_cexp_neg_mul_Ioi27 hz]
    _ = barnesPoleBlock27 z := by
      unfold barnesPoleBlock27
      ring

def barnesPolePoint27 (m : ℕ) (y : ℝ) : ℂ :=
  (((m : ℝ) + 1 / 2 : ℝ) : ℂ) + (y : ℂ) * Complex.I

private def barnesPoleFubini27 (m : ℕ) (y t : ℝ) : ℂ :=
  sechSq27 y * Complex.exp (-barnesPolePoint27 m y * t) *
    (1 + (t : ℂ) / 2)

private theorem integrable_barnesPoleFubini27 (m : ℕ) :
    Integrable (Function.uncurry (barnesPoleFubini27 m))
      (volume.prod (volume.restrict (Ioi 0))) := by
  let a : ℝ := (m : ℝ) + 1 / 2
  have ha : 0 < a := by dsimp [a]; positivity
  have hK : Integrable sechSq27 := by
    apply (integrable_pow_sechSq27 0).congr
    filter_upwards with y
    simp only [pow_zero, Complex.real_smul, Complex.ofReal_one, one_mul]
  have h0 : IntegrableOn (fun t : ℝ => Complex.exp (-(a : ℂ) * t)) (Ioi 0) :=
    integrableOn_cexp_neg_mul27 (z := (a : ℂ)) (by simpa using ha)
  have h1 : IntegrableOn
      (fun t : ℝ => (t : ℂ) * Complex.exp (-(a : ℂ) * t)) (Ioi 0) :=
    integrableOn_mul_cexp_neg27 (z := (a : ℂ)) (by simpa using ha)
  have hG : Integrable
      (fun t : ℝ => Complex.exp (-(a : ℂ) * t) * (1 + (t : ℂ) / 2))
      (volume.restrict (Ioi 0)) := by
    apply (h0.add (h1.div_const 2)).congr_fun _ measurableSet_Ioi
    intro t ht
    change Complex.exp (-(a : ℂ) * t) +
        ((t : ℂ) * Complex.exp (-(a : ℂ) * t)) / 2 =
      Complex.exp (-(a : ℂ) * t) * (1 + (t : ℂ) / 2)
    ring
  have hProd := hK.mul_prod hG
  apply hProd.congr'
  · have hsech : Continuous (fun p : ℝ × ℝ => sechSq27 p.1) := by
      unfold sechSq27
      apply Continuous.div continuous_const
        ((Complex.continuous_ofReal.comp <|
          Real.continuous_cosh.comp
            (continuous_const.mul continuous_fst)).pow 2)
      intro p
      exact pow_ne_zero _
        (Complex.ofReal_ne_zero.mpr (Real.cosh_pos _).ne')
    apply Continuous.aestronglyMeasurable
    unfold Function.uncurry barnesPoleFubini27 barnesPolePoint27
    apply (hsech.mul ?_).mul
    · fun_prop
    · fun_prop
  · filter_upwards with p
    rcases p with ⟨y, t⟩
    dsimp only [Function.uncurry_apply_pair, barnesPoleFubini27,
      barnesPolePoint27]
    simp only [norm_mul, Complex.norm_exp, Complex.mul_re, Complex.neg_re,
      Complex.neg_im, Complex.add_re, Complex.add_im, Complex.ofReal_re,
      Complex.ofReal_im, Complex.I_re, Complex.I_im]
    dsimp only [a]
    ring_nf

private theorem integral_sechSq_phase27 {t : ℝ} (ht : 0 < t) :
    (∫ y : ℝ, Complex.exp (-(y * t) * Complex.I) * sechSq27 y) =
      ((t / (Real.pi * Real.sinh (t / 2)) : ℝ) : ℂ) := by
  have hξ : t / (2 * Real.pi) ≠ 0 := by
    exact div_ne_zero ht.ne' (mul_ne_zero (by norm_num) Real.pi_ne_zero)
  have h := integral_sechSq_cexp27 (t / (2 * Real.pi))
  rw [if_neg hξ] at h
  calc
    (∫ y : ℝ, Complex.exp (-(y * t) * Complex.I) * sechSq27 y) =
        ∫ y : ℝ,
          Complex.exp (-(2 * Real.pi * y * (t / (2 * Real.pi))) * Complex.I) /
            (Real.cosh (Real.pi * y) : ℂ) ^ 2 := by
      apply integral_congr_ae
      filter_upwards with y
      unfold sechSq27
      rw [div_eq_mul_inv]
      congr 1
      · congr 1
        field_simp [Real.pi_ne_zero]
      · simp
    _ = ((2 * (t / (2 * Real.pi)) /
          Real.sinh (Real.pi * (t / (2 * Real.pi))) : ℝ) : ℂ) := by
      convert h using 1 <;> push_cast <;> ring
    _ = ((t / (Real.pi * Real.sinh (t / 2)) : ℝ) : ℂ) := by
      push_cast
      field_simp [Real.pi_ne_zero]

private theorem integral_barnesPoleFubini_inner_left27 (m : ℕ) (y : ℝ) :
    (∫ t : ℝ in Ioi 0, barnesPoleFubini27 m y t) =
      sechSq27 y * barnesPoleBlock27 (barnesPolePoint27 m y) := by
  have hz : 0 < (barnesPolePoint27 m y).re := by
    unfold barnesPolePoint27
    simp only [Complex.add_re, Complex.ofReal_re, Complex.mul_re,
      Complex.ofReal_im, Complex.I_re, Complex.I_im, mul_zero, zero_mul, sub_zero]
    positivity
  calc
    (∫ t : ℝ in Ioi 0, barnesPoleFubini27 m y t) =
        ∫ t : ℝ in Ioi 0, sechSq27 y *
          (Complex.exp (-barnesPolePoint27 m y * t) * (1 + (t : ℂ) / 2)) := by
      unfold barnesPoleFubini27
      apply setIntegral_congr_fun measurableSet_Ioi
      intro t ht
      ring
    _ =
        sechSq27 y * ∫ t : ℝ in Ioi 0,
          Complex.exp (-barnesPolePoint27 m y * t) * (1 + (t : ℂ) / 2) := by
      exact MeasureTheory.integral_const_mul (sechSq27 y) _
    _ = sechSq27 y * barnesPoleBlock27 (barnesPolePoint27 m y) := by
      rw [integral_laplacePoleBlock27 hz]

/-- Absolute integrability of each pole block on a Barnes line. -/
theorem integrable_sechSq_mul_barnesPoleBlock27 (m : ℕ) :
    Integrable (fun y : ℝ =>
      sechSq27 y * barnesPoleBlock27 (barnesPolePoint27 m y)) := by
  have h := (integrable_barnesPoleFubini27 m).integral_prod_left
  apply h.congr
  filter_upwards with y
  exact integral_barnesPoleFubini_inner_left27 m y

private theorem integral_barnesPoleFubini_inner_right27
    (m : ℕ) {t : ℝ} (ht : 0 < t) :
    (∫ y : ℝ, barnesPoleFubini27 m y t) =
      (((Real.exp (-(((m : ℝ) + 1 / 2) * t)) * (1 + t / 2) *
        (t / (Real.pi * Real.sinh (t / 2)))) : ℝ) : ℂ) := by
  let c : ℂ := ((Real.exp (-(((m : ℝ) + 1 / 2) * t)) *
    (1 + t / 2) : ℝ) : ℂ)
  calc
    (∫ y : ℝ, barnesPoleFubini27 m y t) =
        ∫ y : ℝ, c *
          (Complex.exp (-(y * t) * Complex.I) * sechSq27 y) := by
      change (∫ y : ℝ,
        sechSq27 y * Complex.exp (-barnesPolePoint27 m y * t) *
          (1 + (t : ℂ) / 2)) = _
      apply integral_congr_ae
      filter_upwards with y
      unfold barnesPolePoint27 c
      have hexp :
          Complex.exp (-(((((m : ℝ) + 1 / 2 : ℝ) : ℂ) +
            (y : ℂ) * Complex.I)) * (t : ℂ)) =
            Complex.exp (((-(((m : ℝ) + 1 / 2) * t) : ℝ) : ℂ)) *
              Complex.exp (-(y * t) * Complex.I) := by
        rw [← Complex.exp_add]
        congr 1
        push_cast
        ring
      rw [hexp]
      push_cast
      ring
    _ = c * ∫ y : ℝ, Complex.exp (-(y * t) * Complex.I) * sechSq27 y :=
      MeasureTheory.integral_const_mul c _
    _ = c * ((t / (Real.pi * Real.sinh (t / 2)) : ℝ) : ℂ) := by
      rw [integral_sechSq_phase27 ht]
    _ = (((Real.exp (-(((m : ℝ) + 1 / 2) * t)) * (1 + t / 2) *
        (t / (Real.pi * Real.sinh (t / 2)))) : ℝ) : ℂ) := by
      unfold c
      push_cast
      ring

private theorem barnesFourier_to_tail27 (m : ℕ) {t : ℝ} (ht : 0 < t) :
    (Real.pi / 2) *
        (Real.exp (-(((m : ℝ) + 1 / 2) * t)) * (1 + t / 2) *
          (t / (Real.pi * Real.sinh (t / 2)))) =
      Real.exp (-(((m + 1 : ℕ) : ℝ) * t)) * (t + t ^ 2 / 2) /
        (1 - Real.exp (-t)) := by
  have hsinh : Real.sinh (t / 2) ≠ 0 :=
    Real.sinh_ne_zero.mpr (div_ne_zero ht.ne' (by norm_num))
  have hden : 1 - Real.exp (-t) ≠ 0 := by
    have he : Real.exp (-t) < 1 := (Real.exp_lt_one_iff.mpr (by linarith))
    linarith
  have hprod1 :
      Real.exp (-(((m : ℝ) + 1 / 2) * t)) * Real.exp (-t) =
        Real.exp (-(((m : ℝ) + 3 / 2) * t)) := by
    rw [← Real.exp_add]
    congr 1
    ring
  have hprod2 :
      Real.exp (-(((m + 1 : ℕ) : ℝ) * t)) * Real.exp (t / 2) =
        Real.exp (-(((m : ℝ) + 1 / 2) * t)) := by
    rw [← Real.exp_add]
    congr 1
    push_cast
    ring
  have hprod3 :
      Real.exp (-(((m + 1 : ℕ) : ℝ) * t)) * Real.exp (-(t / 2)) =
        Real.exp (-(((m : ℝ) + 3 / 2) * t)) := by
    rw [← Real.exp_add]
    congr 1
    push_cast
    ring
  have hcore :
      Real.exp (-(((m : ℝ) + 1 / 2) * t)) *
          (1 - Real.exp (-t)) =
        2 * Real.exp (-(((m + 1 : ℕ) : ℝ) * t)) *
          Real.sinh (t / 2) := by
    calc
      Real.exp (-(((m : ℝ) + 1 / 2) * t)) *
          (1 - Real.exp (-t)) =
        Real.exp (-(((m : ℝ) + 1 / 2) * t)) -
          Real.exp (-(((m : ℝ) + 1 / 2) * t)) * Real.exp (-t) := by ring
      _ = Real.exp (-(((m : ℝ) + 1 / 2) * t)) -
          Real.exp (-(((m : ℝ) + 3 / 2) * t)) := by rw [hprod1]
      _ = Real.exp (-(((m + 1 : ℕ) : ℝ) * t)) * Real.exp (t / 2) -
          Real.exp (-(((m + 1 : ℕ) : ℝ) * t)) * Real.exp (-(t / 2)) := by
        rw [hprod2, hprod3]
      _ = 2 * Real.exp (-(((m + 1 : ℕ) : ℝ) * t)) *
          Real.sinh (t / 2) := by
        rw [Real.sinh_eq]
        ring
  field_simp [Real.pi_ne_zero, hsinh, hden]
  convert hcore using 1 <;> push_cast <;> ring

/-- Universal pole-block evaluation used in the exact Barnes initial values. -/
theorem integral_barnesPoleBlock27 (m : ℕ) :
    ((Real.pi / 2 : ℝ) : ℂ) *
        ∫ y : ℝ, sechSq27 y * barnesPoleBlock27 (barnesPolePoint27 m y) =
      ((∑' j : ℕ, (1 / (((m + 1 + j : ℕ) : ℝ) ^ 2) +
        1 / (((m + 1 + j : ℕ) : ℝ) ^ 3)) : ℝ) : ℂ) := by
  have hSwap := MeasureTheory.integral_integral_swap
    (integrable_barnesPoleFubini27 m)
  calc
    ((Real.pi / 2 : ℝ) : ℂ) *
        ∫ y : ℝ, sechSq27 y * barnesPoleBlock27 (barnesPolePoint27 m y) =
      ((Real.pi / 2 : ℝ) : ℂ) *
        ∫ y : ℝ, ∫ t : ℝ in Ioi 0, barnesPoleFubini27 m y t := by
      congr 1
      apply integral_congr_ae
      filter_upwards with y
      exact (integral_barnesPoleFubini_inner_left27 m y).symm
    _ = ((Real.pi / 2 : ℝ) : ℂ) *
        ∫ t : ℝ in Ioi 0, ∫ y : ℝ, barnesPoleFubini27 m y t := by
      rw [hSwap]
    _ = ((Real.pi / 2 : ℝ) : ℂ) *
        ∫ t : ℝ in Ioi 0,
          (((Real.exp (-(((m : ℝ) + 1 / 2) * t)) * (1 + t / 2) *
            (t / (Real.pi * Real.sinh (t / 2)))) : ℝ) : ℂ) := by
      congr 1
      apply integral_congr_ae
      filter_upwards [ae_restrict_mem measurableSet_Ioi] with t ht
      exact integral_barnesPoleFubini_inner_right27 m ht
    _ = ∫ t : ℝ in Ioi 0,
        ((Real.pi / 2 : ℝ) : ℂ) *
          (((Real.exp (-(((m : ℝ) + 1 / 2) * t)) * (1 + t / 2) *
            (t / (Real.pi * Real.sinh (t / 2)))) : ℝ) : ℂ) := by
      exact (MeasureTheory.integral_const_mul _ _).symm
    _ = ∫ t : ℝ in Ioi 0,
        (((Real.pi / 2) *
          (Real.exp (-(((m : ℝ) + 1 / 2) * t)) * (1 + t / 2) *
            (t / (Real.pi * Real.sinh (t / 2)))) : ℝ) : ℂ) := by
      apply setIntegral_congr_fun measurableSet_Ioi
      intro t ht
      push_cast
      ring
    _ = ∫ t : ℝ in Ioi 0,
        ((Real.exp (-(((m + 1 : ℕ) : ℝ) * t)) * (t + t ^ 2 / 2) /
          (1 - Real.exp (-t)) : ℝ) : ℂ) := by
      apply setIntegral_congr_fun measurableSet_Ioi
      intro t ht
      exact congrArg (fun x : ℝ => (x : ℂ)) (barnesFourier_to_tail27 m ht)
    _ = (((∫ t : ℝ in Ioi 0,
        Real.exp (-(((m + 1 : ℕ) : ℝ) * t)) * (t + t ^ 2 / 2) /
          (1 - Real.exp (-t))) : ℝ) : ℂ) := by
      exact integral_ofReal
    _ = ((∑' j : ℕ, (1 / (((m + 1 + j : ℕ) : ℝ) ^ 2) +
        1 / (((m + 1 + j : ℕ) : ℝ) ^ 3)) : ℝ) : ℂ) := by
      rw [integral_barnesTail27]

private theorem barnesPoleTail_tsum_eq27 (m : ℕ) :
    (∑' j : ℕ, (1 / (((m + 1 + j : ℕ) : ℝ) ^ 2) +
      1 / (((m + 1 + j : ℕ) : ℝ) ^ 3)) : ℝ) =
      Real.pi ^ 2 / 6 + zeta3 -
        ∑ n ∈ Finset.range (m + 1),
          (1 / (n : ℝ) ^ 2 + 1 / (n : ℝ) ^ 3) := by
  let f2 : ℕ → ℝ := fun n => 1 / (n : ℝ) ^ 2
  let f3 : ℕ → ℝ := fun n => 1 / (n : ℝ) ^ 3
  let f : ℕ → ℝ := fun n => f2 n + f3 n
  have h2 : Summable f2 := by
    exact hasSum_zeta_two.summable
  have h3 : Summable f3 := by
    exact Real.summable_one_div_nat_pow.mpr (by norm_num)
  have hf : Summable f := h2.add h3
  have htotal2 : (∑' n : ℕ, f2 n) = Real.pi ^ 2 / 6 :=
    hasSum_zeta_two.tsum_eq
  have htotal3 : (∑' n : ℕ, f3 n) = zeta3 := by
    have hshift := h3.sum_add_tsum_nat_add 1
    dsimp only [f3] at hshift
    norm_num at hshift
    change (∑' n : ℕ, 1 / (n : ℝ) ^ 3) = zeta3
    simpa only [zeta3, Nat.cast_add, Nat.cast_one, one_div] using hshift.symm
  have htotal : (∑' n : ℕ, f n) = Real.pi ^ 2 / 6 + zeta3 := by
    rw [show (∑' n : ℕ, f n) = (∑' n : ℕ, f2 n) + ∑' n : ℕ, f3 n by
      exact (h2.hasSum.add h3.hasSum).tsum_eq]
    rw [htotal2, htotal3]
  have hsplit := hf.sum_add_tsum_nat_add (m + 1)
  rw [htotal] at hsplit
  have htail :
      (∑' j : ℕ, f (j + (m + 1))) =
        Real.pi ^ 2 / 6 + zeta3 - ∑ n ∈ Finset.range (m + 1), f n := by
    linarith
  rw [← htail]
  apply tsum_congr
  intro j
  unfold f f2 f3
  congr 2 <;> push_cast <;> ring

theorem barnesPoleTail_zero27 :
    (∑' j : ℕ, (1 / (((0 + 1 + j : ℕ) : ℝ) ^ 2) +
      1 / (((0 + 1 + j : ℕ) : ℝ) ^ 3)) : ℝ) =
      Real.pi ^ 2 / 6 + zeta3 := by
  rw [barnesPoleTail_tsum_eq27 0]
  norm_num

theorem barnesPoleTail_one27 :
    (∑' j : ℕ, (1 / (((1 + 1 + j : ℕ) : ℝ) ^ 2) +
      1 / (((1 + 1 + j : ℕ) : ℝ) ^ 3)) : ℝ) =
      Real.pi ^ 2 / 6 + zeta3 - 2 := by
  rw [barnesPoleTail_tsum_eq27 1]
  norm_num [Finset.sum_range_succ]

theorem barnesPoleTail_two27 :
    (∑' j : ℕ, (1 / (((2 + 1 + j : ℕ) : ℝ) ^ 2) +
      1 / (((2 + 1 + j : ℕ) : ℝ) ^ 3)) : ℝ) =
      Real.pi ^ 2 / 6 + zeta3 - 19 / 8 := by
  rw [barnesPoleTail_tsum_eq27 2]
  norm_num [Finset.sum_range_succ]

theorem barnesPoleTail_three27 :
    (∑' j : ℕ, (1 / (((3 + 1 + j : ℕ) : ℝ) ^ 2) +
      1 / (((3 + 1 + j : ℕ) : ℝ) ^ 3)) : ℝ) =
      Real.pi ^ 2 / 6 + zeta3 - 545 / 216 := by
  rw [barnesPoleTail_tsum_eq27 3]
  norm_num [Finset.sum_range_succ]

theorem barnesPoleTail_four27 :
    (∑' j : ℕ, (1 / (((4 + 1 + j : ℕ) : ℝ) ^ 2) +
      1 / (((4 + 1 + j : ℕ) : ℝ) ^ 3)) : ℝ) =
      Real.pi ^ 2 / 6 + zeta3 - 4495 / 1728 := by
  rw [barnesPoleTail_tsum_eq27 4]
  norm_num [Finset.sum_range_succ]

end RamanujanChallenge.P27
