import RamanujanChallenge.Problem27BarnesNormalization
import Mathlib.Analysis.Fourier.FourierTransformDeriv
import Mathlib.Analysis.SpecialFunctions.Gaussian.GaussianIntegral
import Mathlib.Analysis.Calculus.Taylor

open Filter Set MeasureTheory Asymptotics Topology
open scoped FourierTransform

noncomputable section

namespace RamanujanChallenge.P27

def sechSq27 (y : ℝ) : ℂ :=
  1 / (Real.cosh (Real.pi * y) : ℂ) ^ 2

private theorem sechSq27_norm_le_exp (y : ℝ) :
    ‖sechSq27 y‖ ≤ 16 * Real.exp (-2 * Real.pi * |y|) := by
  have hpos : 0 < Real.cosh (Real.pi * y) := Real.cosh_pos _
  have henv := zudilinBarnesEnvelope_le_exp27 y
  have henv0 : 0 ≤ zudilinBarnesEnvelope27 y := by
    unfold zudilinBarnesEnvelope27
    positivity
  calc
    ‖sechSq27 y‖ = zudilinBarnesEnvelope27 y ^ 4 := by
      unfold sechSq27 zudilinBarnesEnvelope27
      rw [norm_div, norm_one, norm_pow, Complex.norm_real,
        Real.norm_eq_abs, abs_of_pos hpos]
      field_simp [Real.sqrt_ne_zero'.mpr hpos]
      rw [show Real.sqrt (Real.cosh (Real.pi * y)) ^ 4 =
          (Real.sqrt (Real.cosh (Real.pi * y)) ^ 2) ^ 2 by ring,
        Real.sq_sqrt hpos.le]
    _ ≤ (2 * Real.exp (-(Real.pi / 2) * |y|)) ^ 4 := by
      exact pow_le_pow_left₀ henv0 henv 4
    _ = 16 * Real.exp (-2 * Real.pi * |y|) := by
      rw [mul_pow]
      norm_num
      rw [← Real.exp_nat_mul]
      congr 2
      ring

theorem integrable_pow_sechSq27 (n : ℕ) :
    Integrable (fun y : ℝ ↦ y ^ n • sechSq27 y) := by
  have hpos : IntegrableOn
      (fun y : ℝ ↦ y ^ (n : ℝ) * Real.exp (-(2 * Real.pi) * y ^ (1 : ℝ)))
      (Ioi 0) := by
    exact integrableOn_rpow_mul_exp_neg_mul_rpow
      (s := (n : ℝ)) (p := 1) (b := 2 * Real.pi)
      (by have hn : (0 : ℝ) ≤ (n : ℝ) := Nat.cast_nonneg n; linarith) (by norm_num)
      (mul_pos (by norm_num) Real.pi_pos)
  let g : ℝ → ℝ := fun y ↦
    16 * (|y| ^ n * Real.exp (-2 * Real.pi * |y|))
  have hmajor : Integrable g := by
    rw [← integrableOn_univ, ← Iic_union_Ioi (a := (0 : ℝ)), integrableOn_union]
    constructor
    · rw [integrableOn_Iic_iff_integrableOn_Iio]
      rw [← (Measure.measurePreserving_neg (volume : Measure ℝ)).integrableOn_comp_preimage
        (Homeomorph.neg ℝ).measurableEmbedding]
      simp only [Function.comp_def, neg_preimage, neg_Iio, neg_zero]
      apply (hpos.const_mul 16).congr
      filter_upwards [ae_restrict_mem measurableSet_Ioi] with y hy
      dsimp only [g]
      rw [abs_neg, abs_of_pos hy]
      simp only [Real.rpow_natCast, Real.rpow_one]
      congr 2
      ring
    · apply (hpos.const_mul 16).congr
      filter_upwards [ae_restrict_mem measurableSet_Ioi] with y hy
      dsimp only [g]
      rw [abs_of_pos hy]
      simp only [Real.rpow_natCast, Real.rpow_one]
      congr 2
      ring
  apply hmajor.mono'
  · simp only [Complex.real_smul]
    apply Continuous.aestronglyMeasurable
    apply Continuous.mul
      (Complex.ofRealCLM.continuous.comp (continuous_id.pow n))
    unfold sechSq27
    apply Continuous.div continuous_const
      ((Complex.ofRealCLM.continuous.comp <|
        Real.continuous_cosh.comp (continuous_const.mul continuous_id)).pow 2)
    intro y
    exact pow_ne_zero _ (Complex.ofReal_ne_zero.mpr (Real.cosh_pos _).ne')
  · filter_upwards with y
    rw [Complex.real_smul, norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_pow]
    calc
      |y| ^ n * ‖sechSq27 y‖ ≤
          |y| ^ n * (16 * Real.exp (-2 * Real.pi * |y|)) := by
        exact mul_le_mul_of_nonneg_left (sechSq27_norm_le_exp y)
          (pow_nonneg (abs_nonneg _) _)
      _ = g y := by
        have hg0 : 0 ≤ g y := by
          dsimp only [g]
          exact mul_nonneg (by norm_num) <|
            mul_nonneg (pow_nonneg (abs_nonneg _) _) (Real.exp_pos _).le
        dsimp only [g]
        ring

private theorem sinh_taylor_three27 :
    (fun x : ℝ ↦ Real.sinh x - (x + x ^ 3 / 6)) =o[nhds 0]
      fun x ↦ (x - 0) ^ 3 := by
  have h := taylor_isLittleO_univ Real.contDiff_sinh
    (x₀ := (0 : ℝ)) (n := 3)
  simp only [taylorWithinEval_succ, taylor_within_zero_eval,
    iteratedDerivWithin_univ] at h
  norm_num [Real.iteratedDeriv_odd_sinh] at h
  convert h using 1 <;> ring_nf

private theorem sinh_sub_id_mul_cosh_taylor_three27 :
    (fun x : ℝ ↦
      (Real.sinh x - x * Real.cosh x) - (-(x ^ 3) / 3)) =o[nhds 0]
      fun x ↦ (x - 0) ^ 3 := by
  have hcont : ContDiff ℝ ⊤ (fun x : ℝ ↦ Real.sinh x - x * Real.cosh x) := by
    fun_prop
  have h := taylor_isLittleO_univ (hcont.of_le (by simp))
    (x₀ := (0 : ℝ)) (n := 3)
  simp only [taylorWithinEval_succ, taylor_within_zero_eval,
    iteratedDerivWithin_univ] at h
  let g : ℝ → ℝ := fun x ↦ Real.sinh x - x * Real.cosh x
  have hg1 (x : ℝ) : HasDerivAt g (-x * Real.sinh x) x := by
    dsimp only [g]
    convert (Real.hasDerivAt_sinh x).sub
      ((hasDerivAt_id x).mul (Real.hasDerivAt_cosh x)) using 1 <;>
      simp only [id_eq] <;> ring
  have hd1 : deriv g = fun x ↦ -x * Real.sinh x := by
    funext x
    exact (hg1 x).deriv
  have hg2 (x : ℝ) : HasDerivAt (fun x ↦ -x * Real.sinh x)
      (-(Real.sinh x + x * Real.cosh x)) x := by
    convert (hasDerivAt_id x).neg.mul (Real.hasDerivAt_sinh x) using 1 <;>
      (try simp only [Pi.neg_apply, id_eq]) <;> ring
  have hd2 : deriv (deriv g) 0 = 0 := by
    rw [hd1]
    simpa using (hg2 0).deriv
  have hg3 (x : ℝ) : HasDerivAt
      (fun x ↦ -(Real.sinh x + x * Real.cosh x))
      (-(2 * Real.cosh x + x * Real.sinh x)) x := by
    convert ((Real.hasDerivAt_sinh x).add
      ((hasDerivAt_id x).mul (Real.hasDerivAt_cosh x))).neg using 1 <;>
      (try simp only [id_eq]) <;> ring
  have hd2fun : deriv (fun x ↦ -x * Real.sinh x) =
      fun x ↦ -(Real.sinh x + x * Real.cosh x) := by
    funext x
    exact (hg2 x).deriv
  have hd3 : deriv (deriv (deriv g)) 0 = -2 := by
    rw [hd1, hd2fun]
    simpa using (hg3 0).deriv
  norm_num [iteratedDeriv_succ] at h
  change deriv (deriv (fun x ↦ Real.sinh x - x * Real.cosh x)) 0 = 0 at hd2
  change deriv (deriv (deriv (fun x ↦ Real.sinh x - x * Real.cosh x))) 0 = -2 at hd3
  rw [hd2, hd3] at h
  norm_num at h
  convert h using 1 <;> ring_nf

private theorem id_sub_sinh_isLittleO_sq27 :
    (fun x : ℝ ↦ x - Real.sinh x) =o[nhds 0] fun x ↦ x ^ 2 := by
  have h := taylor_isLittleO_univ Real.contDiff_sinh
    (x₀ := (0 : ℝ)) (n := 2)
  simp only [taylorWithinEval_succ, taylor_within_zero_eval,
    iteratedDerivWithin_univ] at h
  norm_num [Real.iteratedDeriv_odd_sinh] at h
  have hn := h.neg_left
  convert hn using 1
  · funext x
    ring

private theorem tendsto_id_div_sinh27 :
    Tendsto (fun x : ℝ ↦ x / Real.sinh x)
      (nhdsWithin 0 ({0} : Set ℝ)ᶜ) (nhds 1) := by
  have hdiff := Real.isEquivalent_sinh.isLittleO.tendsto_div_nhds_zero
  have hsdiv : Tendsto (fun x : ℝ ↦ Real.sinh x / x)
      (nhdsWithin 0 ({0} : Set ℝ)ᶜ) (nhds 1) := by
    have hdiff' : Tendsto (fun x : ℝ ↦
        (Real.sinh x - x) / x)
        (nhdsWithin 0 ({0} : Set ℝ)ᶜ) (nhds 0) :=
      by simpa only [Pi.sub_apply, id_eq] using hdiff.mono_left inf_le_left
    have ht := hdiff'.add_const 1
    have heq : (fun x : ℝ ↦ (Real.sinh x - x) / x + 1) =ᶠ[
        nhdsWithin 0 ({0} : Set ℝ)ᶜ] fun x ↦ Real.sinh x / x := by
      filter_upwards [self_mem_nhdsWithin] with x hx
      simp only [mem_compl_iff, mem_singleton_iff] at hx
      field_simp [hx]
      ring
    simpa using ht.congr' heq
  have hinv := hsdiv.inv₀ (by norm_num : (1 : ℝ) ≠ 0)
  have heq : (fun x : ℝ ↦ (Real.sinh x / x)⁻¹) =ᶠ[
      nhdsWithin 0 ({0} : Set ℝ)ᶜ] fun x ↦ x / Real.sinh x := by
    filter_upwards [self_mem_nhdsWithin] with x hx
    simp only [mem_compl_iff, mem_singleton_iff] at hx
    field_simp [hx, Real.sinh_ne_zero.mpr hx]
  simpa using hinv.congr' heq

private def sinhInvRatio27 (x : ℝ) : ℝ :=
  if x = 0 then 1 else x / Real.sinh x

private def sinhInvRatioDeriv27 (x : ℝ) : ℝ :=
  if x = 0 then 0 else
    (Real.sinh x - x * Real.cosh x) / Real.sinh x ^ 2

private theorem hasDerivAt_sinhInvRatio_zero27 :
    HasDerivAt sinhInvRatio27 0 0 := by
  rw [hasDerivAt_iff_tendsto_slope_zero]
  have hrem := id_sub_sinh_isLittleO_sq27.tendsto_div_nhds_zero
  have hprod := (hrem.mono_left inf_le_left).mul tendsto_id_div_sinh27
  convert hprod using 1
  · funext x
    by_cases hx : x = 0
    · subst x
      simp [sinhInvRatio27]
    · simp only [zero_add, sinhInvRatio27, if_neg hx, if_pos, smul_eq_mul]
      field_simp [hx, Real.sinh_ne_zero.mpr hx]
  · norm_num

private theorem hasDerivAt_sinhInvRatio27 (x : ℝ) :
    HasDerivAt sinhInvRatio27 (sinhInvRatioDeriv27 x) x := by
  by_cases hx : x = 0
  · subst x
    simpa [sinhInvRatioDeriv27] using hasDerivAt_sinhInvRatio_zero27
  · rw [sinhInvRatioDeriv27, if_neg hx]
    have hbase := (hasDerivAt_id x).div (Real.hasDerivAt_sinh x)
      (Real.sinh_ne_zero.mpr hx)
    apply (hbase.congr_deriv (by simp only [id_eq]; ring)).congr_of_eventuallyEq
      (f₁ := sinhInvRatio27)
    · filter_upwards [eventually_ne_nhds hx] with y hy
      simp [sinhInvRatio27, hy]

private theorem hasDerivAt_sinhInvRatioDeriv_zero27 :
    HasDerivAt sinhInvRatioDeriv27 (-(1 : ℝ) / 3) 0 := by
  rw [hasDerivAt_iff_tendsto_slope_zero]
  have hrem :=
    sinh_sub_id_mul_cosh_taylor_three27.tendsto_div_nhds_zero
  have hlead : Tendsto
      (fun x : ℝ ↦ (Real.sinh x - x * Real.cosh x) / x ^ 3)
      (nhdsWithin 0 ({0} : Set ℝ)ᶜ) (nhds (-(1 : ℝ) / 3)) := by
    have hrem' : Tendsto
        (fun x : ℝ ↦
          ((Real.sinh x - x * Real.cosh x) - (-(x ^ 3) / 3)) /
            (x - 0) ^ 3)
        (nhdsWithin 0 ({0} : Set ℝ)ᶜ) (nhds 0) :=
      hrem.mono_left inf_le_left
    have ht := hrem'.add_const (-(1 : ℝ) / 3)
    have heq : (fun x : ℝ ↦
        ((Real.sinh x - x * Real.cosh x) - (-(x ^ 3) / 3)) /
          (x - 0) ^ 3 + (-(1 : ℝ) / 3)) =ᶠ[
          nhdsWithin 0 ({0} : Set ℝ)ᶜ]
        fun x ↦ (Real.sinh x - x * Real.cosh x) / x ^ 3 := by
      filter_upwards [self_mem_nhdsWithin] with x hx
      simp only [mem_compl_iff, mem_singleton_iff] at hx
      field_simp [hx]
      ring
    simpa using ht.congr' heq
  have hprod := hlead.mul (tendsto_id_div_sinh27.pow 2)
  convert hprod using 1
  · funext x
    by_cases hx : x = 0
    · subst x
      simp [sinhInvRatioDeriv27]
    · simp only [zero_add, sinhInvRatioDeriv27, if_neg hx, if_pos, smul_eq_mul]
      field_simp [hx, Real.sinh_ne_zero.mpr hx]
      ring
  · norm_num

private theorem iteratedDeriv_two_sinhInvRatio27 :
    iteratedDeriv 2 sinhInvRatio27 0 = -(1 : ℝ) / 3 := by
  have hderiv : deriv sinhInvRatio27 = sinhInvRatioDeriv27 := by
    funext x
    exact (hasDerivAt_sinhInvRatio27 x).deriv
  calc
    iteratedDeriv 2 sinhInvRatio27 0 = deriv (deriv sinhInvRatio27) 0 := by
      rw [show (2 : ℕ) = 1 + 1 by norm_num, iteratedDeriv_succ,
        iteratedDeriv_one]
    _ = deriv sinhInvRatioDeriv27 0 := by rw [hderiv]
    _ = -(1 : ℝ) / 3 := hasDerivAt_sinhInvRatioDeriv_zero27.deriv

private def sechSqFourierClosed27 (x : ℝ) : ℝ :=
  (2 / Real.pi) * sinhInvRatio27 (Real.pi * x)

private theorem iteratedDeriv_two_sechSqFourierClosed27 :
    iteratedDeriv 2 sechSqFourierClosed27 0 = -2 * Real.pi / 3 := by
  have hfirst (x : ℝ) : HasDerivAt sechSqFourierClosed27
      ((2 / Real.pi) *
        (sinhInvRatioDeriv27 (Real.pi * x) * Real.pi)) x := by
    unfold sechSqFourierClosed27
    simpa only [Function.comp_apply, mul_one] using
      (((hasDerivAt_sinhInvRatio27 (Real.pi * x)).comp x
        ((hasDerivAt_id x).const_mul Real.pi)).const_mul (2 / Real.pi))
  have hderiv : deriv sechSqFourierClosed27 = fun x ↦
      (2 / Real.pi) *
        (sinhInvRatioDeriv27 (Real.pi * x) * Real.pi) := by
    funext x
    exact (hfirst x).deriv
  have hsecond : HasDerivAt
      (fun x : ℝ ↦ (2 / Real.pi) *
        (sinhInvRatioDeriv27 (Real.pi * x) * Real.pi))
      (-2 * Real.pi / 3) 0 := by
    have hc : HasDerivAt
        (fun x : ℝ ↦ sinhInvRatioDeriv27 (Real.pi * x))
        (-(1 : ℝ) / 3 * Real.pi) 0 := by
      have hinner : HasDerivAt (fun x : ℝ ↦ Real.pi * x) Real.pi 0 := by
        simpa only [id_eq, mul_one] using
          (hasDerivAt_id (0 : ℝ)).const_mul Real.pi
      have hz : HasDerivAt sinhInvRatioDeriv27 (-(1 : ℝ) / 3)
          (Real.pi * 0) := by
        simpa only [mul_zero] using hasDerivAt_sinhInvRatioDeriv_zero27
      simpa only [Function.comp_apply, mul_zero] using
        hz.comp 0 hinner
    convert ((hc.mul_const Real.pi).const_mul (2 / Real.pi)) using 1
    field_simp [Real.pi_ne_zero]
  calc
    iteratedDeriv 2 sechSqFourierClosed27 0 =
        deriv (deriv sechSqFourierClosed27) 0 := by
      rw [show (2 : ℕ) = 1 + 1 by norm_num, iteratedDeriv_succ,
        iteratedDeriv_one]
    _ = deriv (fun x ↦ (2 / Real.pi) *
        (sinhInvRatioDeriv27 (Real.pi * x) * Real.pi)) 0 := by rw [hderiv]
    _ = -2 * Real.pi / 3 := hsecond.deriv

private theorem fourier_sechSq27_eq_closed :
    𝓕 sechSq27 = fun x ↦ (sechSqFourierClosed27 x : ℂ) := by
  funext x
  rw [Real.fourier_eq']
  calc
    (∫ v : ℝ, Complex.exp ((↑(-2 * Real.pi * inner ℝ v x) * Complex.I)) •
        sechSq27 v) =
        ∫ y : ℝ, Complex.exp (-(2 * Real.pi * y * x) * Complex.I) /
          (Real.cosh (Real.pi * y) : ℂ) ^ 2 := by
      apply integral_congr_ae
      filter_upwards with y
      rw [show inner ℝ y x = y * x by change x * y = y * x; ring]
      simp only [sechSq27, smul_eq_mul]
      change Complex.exp (((-2 * Real.pi * (y * x) : ℝ) : ℂ) * Complex.I) *
          (1 / (Real.cosh (Real.pi * y) : ℂ) ^ 2) =
        Complex.exp (-(2 * Real.pi * y * x) * Complex.I) /
          (Real.cosh (Real.pi * y) : ℂ) ^ 2
      rw [div_eq_mul_inv]
      congr 1
      · push_cast
        ring
      · simp
    _ = (if x = 0 then (2 / Real.pi : ℝ) else
          2 * x / Real.sinh (Real.pi * x)) := integral_sechSq_cexp27 x
    _ = (sechSqFourierClosed27 x : ℂ) := by
      by_cases hx : x = 0
      · subst x
        simp [sechSqFourierClosed27, sinhInvRatio27]
      · rw [if_neg hx]
        unfold sechSqFourierClosed27 sinhInvRatio27
        rw [if_neg (mul_ne_zero Real.pi_ne_zero hx)]
        push_cast
        field_simp [Real.pi_ne_zero]

private theorem iteratedDeriv_two_coe_sechSqFourierClosed27 :
    iteratedDeriv 2 (fun x ↦ (sechSqFourierClosed27 x : ℂ)) 0 =
      ((-2 * Real.pi / 3 : ℝ) : ℂ) := by
  have hfirst (x : ℝ) : HasDerivAt sechSqFourierClosed27
      ((2 / Real.pi) *
        (sinhInvRatioDeriv27 (Real.pi * x) * Real.pi)) x := by
    unfold sechSqFourierClosed27
    simpa only [Function.comp_apply, mul_one] using
      (((hasDerivAt_sinhInvRatio27 (Real.pi * x)).comp x
        ((hasDerivAt_id x).const_mul Real.pi)).const_mul (2 / Real.pi))
  have hderiv : deriv (fun x ↦ (sechSqFourierClosed27 x : ℂ)) = fun x ↦
      (((2 / Real.pi) *
        (sinhInvRatioDeriv27 (Real.pi * x) * Real.pi) : ℝ) : ℂ) := by
    funext x
    exact (hfirst x).ofReal_comp.deriv
  have hsecond : HasDerivAt
      (fun x : ℝ ↦ (2 / Real.pi) *
        (sinhInvRatioDeriv27 (Real.pi * x) * Real.pi))
      (-2 * Real.pi / 3) 0 := by
    have hc : HasDerivAt
        (fun x : ℝ ↦ sinhInvRatioDeriv27 (Real.pi * x))
        (-(1 : ℝ) / 3 * Real.pi) 0 := by
      have hinner : HasDerivAt (fun x : ℝ ↦ Real.pi * x) Real.pi 0 := by
        simpa only [id_eq, mul_one] using
          (hasDerivAt_id (0 : ℝ)).const_mul Real.pi
      have hz : HasDerivAt sinhInvRatioDeriv27 (-(1 : ℝ) / 3)
          (Real.pi * 0) := by
        simpa only [mul_zero] using hasDerivAt_sinhInvRatioDeriv_zero27
      simpa only [Function.comp_apply, mul_zero] using
        hz.comp 0 hinner
    convert ((hc.mul_const Real.pi).const_mul (2 / Real.pi)) using 1
    field_simp [Real.pi_ne_zero]
  calc
    iteratedDeriv 2 (fun x ↦ (sechSqFourierClosed27 x : ℂ)) 0 =
        deriv (deriv (fun x ↦ (sechSqFourierClosed27 x : ℂ))) 0 := by
      rw [show (2 : ℕ) = 1 + 1 by norm_num, iteratedDeriv_succ,
        iteratedDeriv_one]
    _ = deriv (fun x ↦
        (((2 / Real.pi) *
          (sinhInvRatioDeriv27 (Real.pi * x) * Real.pi) : ℝ) : ℂ)) 0 := by
      rw [hderiv]
    _ = ((-2 * Real.pi / 3 : ℝ) : ℂ) := hsecond.ofReal_comp.deriv

private theorem iteratedDeriv_two_fourier_sechSq27 :
    iteratedDeriv 2 (𝓕 sechSq27) 0 = ((-2 * Real.pi / 3 : ℝ) : ℂ) := by
  rw [fourier_sechSq27_eq_closed]
  exact iteratedDeriv_two_coe_sechSqFourierClosed27

private theorem fourier_iteratedDeriv_two_sechSq27 :
    iteratedDeriv 2 (𝓕 sechSq27) =
      𝓕 (fun x : ℝ ↦ (-2 * Real.pi * Complex.I * x) ^ 2 • sechSq27 x) := by
  apply Real.iteratedDeriv_fourier (N := (2 : ℕ∞))
  · intro n hn
    exact integrable_pow_sechSq27 n
  · norm_num

theorem integral_sq_sechSq27 :
    (∫ y : ℝ, (y : ℂ) ^ 2 * sechSq27 y) =
      ((1 / (6 * Real.pi) : ℝ) : ℂ) := by
  have h := congrFun fourier_iteratedDeriv_two_sechSq27 0
  rw [iteratedDeriv_two_fourier_sechSq27, Real.fourier_eq'] at h
  simp only [inner_zero_right, mul_zero, Complex.ofReal_zero, zero_mul,
    Complex.exp_zero, one_mul, one_smul] at h
  have hInt :
      (∫ v : ℝ, (-2 * Real.pi * Complex.I * v) ^ 2 • sechSq27 v) =
        (-4 * (Real.pi : ℂ) ^ 2) *
          ∫ v : ℝ, (v : ℂ) ^ 2 * sechSq27 v := by
    calc
      (∫ v : ℝ, (-2 * Real.pi * Complex.I * v) ^ 2 • sechSq27 v) =
          ∫ v : ℝ, (-4 * (Real.pi : ℂ) ^ 2) *
            ((v : ℂ) ^ 2 * sechSq27 v) := by
        apply integral_congr_ae
        filter_upwards with v
        simp only [smul_eq_mul]
        rw [show
          (-2 * (Real.pi : ℂ) * Complex.I * (v : ℂ)) ^ 2 * sechSq27 v =
            ((-2 * (Real.pi : ℂ) * (v : ℂ)) ^ 2 * Complex.I ^ 2) *
              sechSq27 v by ring]
        rw [Complex.I_sq]
        ring
      _ = (-4 * (Real.pi : ℂ) ^ 2) *
            ∫ v : ℝ, (v : ℂ) ^ 2 * sechSq27 v :=
        MeasureTheory.integral_const_mul _ _
  rw [hInt] at h
  calc
    (∫ y : ℝ, (y : ℂ) ^ 2 * sechSq27 y) =
        ((-2 * Real.pi / 3 : ℝ) : ℂ) /
          (-4 * (Real.pi : ℂ) ^ 2) := by
      apply (eq_div_iff (mul_ne_zero (by norm_num)
        (pow_ne_zero 2 (Complex.ofReal_ne_zero.mpr Real.pi_ne_zero)))).2
      simpa [mul_comm] using h.symm
    _ = ((1 / (6 * Real.pi) : ℝ) : ℂ) := by
      push_cast
      field_simp [Real.pi_ne_zero]
      ring

/-- Total mass of the squared hyperbolic secant kernel. -/
theorem integral_sechSq27 :
    (∫ y : ℝ, sechSq27 y) = ((2 / Real.pi : ℝ) : ℂ) := by
  have h := integral_sechSq_cexp27 0
  simpa [sechSq27] using h

end RamanujanChallenge.P27
