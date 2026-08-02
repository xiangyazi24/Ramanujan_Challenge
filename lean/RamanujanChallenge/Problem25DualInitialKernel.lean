import RamanujanChallenge.Problem25Integral
import Mathlib.MeasureTheory.Integral.IntervalIntegral.IntegrationByParts
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Chebyshev.Orthogonality
import Mathlib.Analysis.SpecialFunctions.Log.NegMulLog

noncomputable section

namespace RamanujanChallenge.P25

open MeasureTheory Set Filter Topology
open scoped Interval

private def logOneIntegrand (x : ℝ) : ℝ :=
  Real.log (1 + x) / (1 + x ^ 2)

private def cauchyIntegrand (x : ℝ) : ℝ :=
  1 / (1 + x ^ 2)

private theorem mobius_hasDerivAt (x : ℝ) (hx : x ≠ -1) :
    HasDerivAt (fun y : ℝ => (1 - y) / (1 + y))
      (-2 / (1 + x) ^ 2) x := by
  have hden : 1 + x ≠ 0 := by
    intro h
    apply hx
    linarith
  have h := ((hasDerivAt_const x 1).sub (hasDerivAt_id x)).div
    ((hasDerivAt_const x 1).add (hasDerivAt_id x)) hden
  convert h using 1
  simp only [Pi.add_apply, Pi.sub_apply, id_eq]
  field_simp [hden]
  ring

private theorem logOneIntegral_symmetry :
    2 * (∫ x in (0 : ℝ)..1, logOneIntegrand x) =
      Real.log 2 * (∫ x in (0 : ℝ)..1, cauchyIntegrand x) := by
  let f : ℝ → ℝ := fun x => (1 - x) / (1 + x)
  let f' : ℝ → ℝ := fun x => -2 / (1 + x) ^ 2
  have hf : ContinuousOn f [[(0 : ℝ), 1]] := by
    apply ContinuousOn.div
    · fun_prop
    · fun_prop
    · intro x hx
      simp only [uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1), mem_Icc] at hx
      linarith
  have hderiv : ∀ x ∈ Ioo (min (0 : ℝ) 1) (max (0 : ℝ) 1),
      HasDerivAt f (f' x) x := by
    intro x hx
    exact mobius_hasDerivAt x (by norm_num at hx; linarith)
  have hnonpos : ∀ x ∈ Ioo (min (0 : ℝ) 1) (max (0 : ℝ) 1),
      f' x ≤ 0 := by
    intro x hx
    dsimp [f']
    exact div_nonpos_of_nonpos_of_nonneg (by norm_num) (sq_nonneg (1 + x))
  have hsub := intervalIntegral.integral_comp_mul_deriv_of_deriv_nonpos
    (g := logOneIntegrand) hf hderiv hnonpos
  have hends : f 0 = 1 ∧ f 1 = 0 := by norm_num [f]
  rw [hends.1, hends.2] at hsub
  rw [intervalIntegral.integral_symm (f := logOneIntegrand) 0 1] at hsub
  have hleft :
      (∫ x in (0 : ℝ)..1, (logOneIntegrand ∘ f) x * f' x) =
        (∫ x in (0 : ℝ)..1,
          logOneIntegrand x - Real.log 2 * cauchyIntegrand x) := by
    apply intervalIntegral.integral_congr
    intro x hx
    simp only [uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1), mem_Icc] at hx
    have hx1 : 1 + x ≠ 0 := by linarith
    have htwo : (2 : ℝ) ≠ 0 := by norm_num
    dsimp [logOneIntegrand, cauchyIntegrand, f, f', Function.comp_apply]
    have hlogarg : 1 + (1 - x) / (1 + x) = 2 / (1 + x) := by
      field_simp [hx1]
      ring
    rw [hlogarg, Real.log_div (by norm_num) hx1]
    field_simp [hx1]
    ring
  rw [hleft] at hsub
  rw [intervalIntegral.integral_sub] at hsub
  · rw [intervalIntegral.integral_const_mul] at hsub
    linarith
  · apply ContinuousOn.intervalIntegrable
    unfold logOneIntegrand
    apply ContinuousOn.div
    · apply ContinuousOn.log (continuousOn_const.add continuousOn_id)
      intro x hx
      simp only [uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1), mem_Icc] at hx
      change 1 + x ≠ 0
      linarith
    · fun_prop
    · intro x hx
      positivity
  · exact hcauchyIntegrable.const_mul (Real.log 2)

where
  hcauchyIntegrable : IntervalIntegrable cauchyIntegrand volume 0 1 := by
    apply ContinuousOn.intervalIntegrable
    unfold cauchyIntegrand
    apply ContinuousOn.div
    · fun_prop
    · fun_prop
    · intro x hx
      positivity

private def catalanKernelIntegrand (t : ℝ) : ℝ :=
  Real.log ((1 + t) / t) / Real.sqrt (1 - t ^ 2)

def quadraticMap (x : ℝ) : ℝ :=
  2 * x / (1 + x ^ 2)

def quadraticMapDeriv (x : ℝ) : ℝ :=
  2 * (1 - x ^ 2) / (1 + x ^ 2) ^ 2

theorem quadraticMap_hasDerivAt (x : ℝ) :
    HasDerivAt quadraticMap (quadraticMapDeriv x) x := by
  have hden : 1 + x ^ 2 ≠ 0 := by positivity
  have h := ((hasDerivAt_const x 2).mul (hasDerivAt_id x)).div
    ((hasDerivAt_const x 1).add ((hasDerivAt_id x).pow 2)) hden
  convert h using 1
  simp only [Pi.add_apply, Pi.mul_apply, id_eq]
  dsimp [quadraticMap, quadraticMapDeriv]
  field_simp [hden]
  ring

private theorem catalanKernelIntegral :
    (∫ t in (0 : ℝ)..1, catalanKernelIntegrand t) = 2 * catalanConstant := by
  have hmap : ContinuousOn quadraticMap [[(0 : ℝ), 1]] := by
    unfold quadraticMap
    apply ContinuousOn.div
    · fun_prop
    · fun_prop
    · intro x hx
      positivity
  have hderiv : ∀ x ∈ Ioo (min (0 : ℝ) 1) (max (0 : ℝ) 1),
      HasDerivAt quadraticMap (quadraticMapDeriv x) x := by
    intro x hx
    exact quadraticMap_hasDerivAt x
  have hnonneg : ∀ x ∈ Ioo (min (0 : ℝ) 1) (max (0 : ℝ) 1),
      0 ≤ quadraticMapDeriv x := by
    intro x hx
    norm_num at hx
    dsimp [quadraticMapDeriv]
    exact div_nonneg (by nlinarith [sq_nonneg x]) (sq_nonneg _)
  have hsub := intervalIntegral.integral_comp_mul_deriv_of_deriv_nonneg
    (g := catalanKernelIntegrand) hmap hderiv hnonneg
  have hends : quadraticMap 0 = 0 ∧ quadraticMap 1 = 1 := by
    norm_num [quadraticMap]
  rw [hends.1, hends.2] at hsub
  have hpoint : ∀ᵐ x : ℝ ∂volume,
      x ∈ Ι (0 : ℝ) 1 →
        (catalanKernelIntegrand ∘ quadraticMap) x * quadraticMapDeriv x =
          2 * (2 * Real.log (1 + x) - Real.log 2 - Real.log x) /
            (1 + x ^ 2) := by
    filter_upwards [MeasureTheory.Measure.ae_ne volume (0 : ℝ),
      MeasureTheory.Measure.ae_ne volume (1 : ℝ)] with x hx0 hx1 hxmem
    simp only [uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1), mem_Ioc] at hxmem
    have hxpos : 0 < x := hxmem.1
    have hxlt : x < 1 := lt_of_le_of_ne hxmem.2 hx1
    have hx1p : 0 < 1 + x := by linarith
    have hden : 0 < 1 + x ^ 2 := by positivity
    have hphi : 0 < quadraticMap x := by
      dsimp [quadraticMap]
      positivity
    have hphi1 : quadraticMap x < 1 := by
      dsimp [quadraticMap]
      rw [div_lt_one hden]
      nlinarith [sq_pos_of_pos (sub_pos.mpr hxlt)]
    have hsquare :
        1 - quadraticMap x ^ 2 = ((1 - x ^ 2) / (1 + x ^ 2)) ^ 2 := by
      dsimp [quadraticMap]
      field_simp [hden.ne']
      ring
    have hsqrt :
        Real.sqrt (1 - quadraticMap x ^ 2) =
          (1 - x ^ 2) / (1 + x ^ 2) := by
      rw [hsquare, Real.sqrt_sq_eq_abs, abs_of_nonneg]
      exact div_nonneg (by nlinarith [sq_nonneg x]) hden.le
    have hlogarg :
        (1 + quadraticMap x) / quadraticMap x =
          (1 + x) ^ 2 / (2 * x) := by
      dsimp [quadraticMap]
      field_simp [hden.ne', hxpos.ne']
      ring
    have hlog :
        Real.log ((1 + quadraticMap x) / quadraticMap x) =
          2 * Real.log (1 + x) - Real.log 2 - Real.log x := by
      rw [hlogarg, Real.log_div (pow_ne_zero 2 hx1p.ne')
        (mul_ne_zero (by norm_num) hxpos.ne'), Real.log_pow]
      rw [Real.log_mul (by norm_num) hxpos.ne']
      ring
    dsimp [catalanKernelIntegrand, Function.comp_apply]
    rw [hsqrt, hlog]
    dsimp [quadraticMapDeriv]
    have hxsub : 1 - x ^ 2 ≠ 0 := by nlinarith
    field_simp [hden.ne', hxsub] <;> ring
  have hcomp :
      (∫ x in (0 : ℝ)..1,
        (catalanKernelIntegrand ∘ quadraticMap) x * quadraticMapDeriv x) =
      ∫ x in (0 : ℝ)..1,
        2 * (2 * Real.log (1 + x) - Real.log 2 - Real.log x) /
          (1 + x ^ 2) :=
    intervalIntegral.integral_congr_ae hpoint
  rw [hcomp] at hsub
  have hlogOneInt : IntervalIntegrable logOneIntegrand volume 0 1 := by
    apply ContinuousOn.intervalIntegrable
    unfold logOneIntegrand
    apply ContinuousOn.div
    · apply ContinuousOn.log (continuousOn_const.add continuousOn_id)
      intro x hx
      simp only [uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1), mem_Icc] at hx
      change 1 + x ≠ 0
      linarith
    · fun_prop
    · intro x hx
      positivity
  have hcauchyInt : IntervalIntegrable cauchyIntegrand volume 0 1 := by
    apply ContinuousOn.intervalIntegrable
    unfold cauchyIntegrand
    apply ContinuousOn.div
    · fun_prop
    · fun_prop
    · intro x hx
      positivity
  have hcatalanInt : IntervalIntegrable
      (fun x : ℝ => (-Real.log x) / (1 + x ^ 2)) volume 0 1 := by
    have hbase : IntervalIntegrable (fun x : ℝ => -Real.log x) volume 0 1 :=
      intervalIntegral.intervalIntegrable_log'.neg
    have hcont : ContinuousOn (fun x : ℝ => (1 + x ^ 2)⁻¹) [[(0 : ℝ), 1]] := by
      apply ContinuousOn.inv₀
      · fun_prop
      · intro x hx
        positivity
    simpa only [div_eq_mul_inv] using hbase.mul_continuousOn hcont
  have hexpand :
      (∫ x in (0 : ℝ)..1,
        2 * (2 * Real.log (1 + x) - Real.log 2 - Real.log x) /
          (1 + x ^ 2)) =
        4 * (∫ x in (0 : ℝ)..1, logOneIntegrand x) -
          2 * Real.log 2 * (∫ x in (0 : ℝ)..1, cauchyIntegrand x) +
          2 * (∫ x in (0 : ℝ)..1, (-Real.log x) / (1 + x ^ 2)) := by
    have hfun : (fun x : ℝ =>
        2 * (2 * Real.log (1 + x) - Real.log 2 - Real.log x) /
          (1 + x ^ 2)) =
        (fun x => 4 * logOneIntegrand x -
          (2 * Real.log 2) * cauchyIntegrand x +
          2 * ((-Real.log x) / (1 + x ^ 2))) := by
      funext x
      dsimp [logOneIntegrand, cauchyIntegrand]
      ring
    have hfirst : IntervalIntegrable
        (fun x => 4 * logOneIntegrand x) volume 0 1 := hlogOneInt.const_mul 4
    have hsecond : IntervalIntegrable
        (fun x => (2 * Real.log 2) * cauchyIntegrand x) volume 0 1 :=
      hcauchyInt.const_mul (2 * Real.log 2)
    have hthird : IntervalIntegrable
        (fun x => 2 * ((-Real.log x) / (1 + x ^ 2))) volume 0 1 :=
      hcatalanInt.const_mul 2
    rw [hfun, intervalIntegral.integral_add (hfirst.sub hsecond) hthird,
      intervalIntegral.integral_sub hfirst hsecond,
      intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul,
      intervalIntegral.integral_const_mul]
  rw [hexpand] at hsub
  rw [← catalanConstant_eq_integral] at hsub
  have hsym := logOneIntegral_symmetry
  linarith

private theorem catalanKernel_intervalIntegrable :
    IntervalIntegrable catalanKernelIntegrand volume 0 1 := by
  let w : ℝ → ℝ := fun x => (Real.sqrt (1 - x ^ 2))⁻¹
  have hw : IntervalIntegrable w volume 0 1 := by
    have h := Polynomial.Chebyshev.intervalIntegrable_sqrt_one_sub_sq_inv.mono_set
      (c := (0 : ℝ)) (d := 1) (by
        intro x hx
        simp only [uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1), mem_Icc] at hx
        simp only [uIcc_of_le (by norm_num : (-1 : ℝ) ≤ 1), mem_Icc]
        constructor <;> linarith)
    simpa [w] using h
  have hlogOne : ContinuousOn (fun x : ℝ => Real.log (1 + x)) [[(0 : ℝ), 1]] := by
    apply ContinuousOn.log (continuousOn_const.add continuousOn_id)
    intro x hx
    simp only [uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1), mem_Icc] at hx
    change 1 + x ≠ 0
    linarith
  have hfirst : IntervalIntegrable (fun x => Real.log (1 + x) * w x)
      volume 0 1 := hw.continuousOn_mul hlogOne
  have hwLeft : ContinuousOn w [[(0 : ℝ), (1 / 2 : ℝ)]] := by
    dsimp [w]
    apply ContinuousOn.inv₀
    · fun_prop
    · intro x hx
      simp only [uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2), mem_Icc] at hx
      have : 0 < 1 - x ^ 2 := by nlinarith [sq_nonneg x]
      exact (Real.sqrt_pos.2 this).ne'
  have hlogwLeft : IntervalIntegrable (fun x => Real.log x * w x)
      volume 0 (1 / 2) :=
    intervalIntegral.intervalIntegrable_log'.mul_continuousOn hwLeft
  have hwRight : IntervalIntegrable w volume (1 / 2) 1 := by
    apply hw.mono_set
    intro x hx
    simp only [uIcc_of_le (by norm_num : (1 / 2 : ℝ) ≤ 1), mem_Icc] at hx
    simp only [uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1), mem_Icc]
    constructor <;> linarith
  have hlogRight : ContinuousOn Real.log [[(1 / 2 : ℝ), 1]] := by
    exact Real.continuousOn_log.mono (by
      intro x hx
      simp only [uIcc_of_le (by norm_num : (1 / 2 : ℝ) ≤ 1), mem_Icc] at hx
      exact (by linarith : x ≠ 0))
  have hlogwRight : IntervalIntegrable (fun x => Real.log x * w x)
      volume (1 / 2) 1 := hwRight.continuousOn_mul hlogRight
  have hsecond : IntervalIntegrable (fun x => Real.log x * w x)
      volume 0 1 := hlogwLeft.trans hlogwRight
  have hdiff := hfirst.sub hsecond
  convert hdiff using 1
  funext x
  dsimp [catalanKernelIntegrand, w]
  by_cases hx : x = 0
  · subst x
    norm_num
  by_cases hxneg : x = -1
  · subst x
    norm_num
  rw [Real.log_div (by intro h; apply hxneg; linarith) hx]
  ring

def logRatio (t : ℝ) : ℝ :=
  Real.log (1 + t) - Real.log t

private def denA (t : ℝ) : ℝ :=
  -75 * (300 * t ^ 8 + 750 * t ^ 7 + 550 * t ^ 6 + 75 * t ^ 5 -
      15 * t ^ 4 + 5 * t ^ 3 + 7 * t ^ 2 + 6 * t + 2) /
    (32 * (1 + t) ^ 3)

private def denAPrime (t : ℝ) : ℝ :=
  -75 * t * (750 * t ^ 7 + 2700 * t ^ 6 + 3450 * t ^ 5 +
      1725 * t ^ 4 + 180 * t ^ 3 - 30 * t ^ 2 + 4 * t + 1) /
    (16 * (1 + t) ^ 4)

private def denB (t : ℝ) : ℝ := 5625 / 8 * t ^ 6

private def denBPrime (t : ℝ) : ℝ := 16875 / 4 * t ^ 5

def denR (t : ℝ) : ℝ :=
  (-375 : ℝ) / 32 * t ^ 3 *
      (-420 * t ^ 6 - 1050 * t ^ 5 - 410 * t ^ 4 + 795 * t ^ 3 +
        681 * t ^ 2 + 83 * t - 15) /
    (1 + t) ^ 3

def denS (t : ℝ) : ℝ :=
  (-375 : ℝ) / 32 * t ^ 3 *
      (420 * t ^ 7 + 1260 * t ^ 6 + 900 * t ^ 5 - 660 * t ^ 4 -
        1080 * t ^ 3 - 360 * t ^ 2) /
    (1 + t) ^ 3

def denReducedIntegrand (t : ℝ) : ℝ :=
  (denR t + denS t * logRatio t) / Real.sqrt (1 - t ^ 2)

private def denPrimitive (t : ℝ) : ℝ :=
  Real.sqrt (1 - t ^ 2) * (denA t + denB t * logRatio t)

private theorem denB_operator (t : ℝ) (ht : 1 + t ≠ 0) :
    (1 - t ^ 2) * denBPrime t - t * denB t = denS t := by
  unfold denBPrime denB denS
  field_simp [ht]
  ring

private theorem denA_operator (t : ℝ) (ht0 : t ≠ 0) (ht1 : 1 + t ≠ 0) :
    (1 - t ^ 2) * denAPrime t - t * denA t -
        (1 - t) * denB t / t = denR t := by
  unfold denAPrime denA denB denR
  field_simp [ht0, ht1]
  ring

private theorem denA_hasDerivAt (t : ℝ) (ht : 1 + t ≠ 0) :
    HasDerivAt denA (denAPrime t) t := by
  let P : ℝ → ℝ := fun x =>
    300 * x ^ 8 + 750 * x ^ 7 + 550 * x ^ 6 + 75 * x ^ 5 -
      15 * x ^ 4 + 5 * x ^ 3 + 7 * x ^ 2 + 6 * x + 2
  have hP : HasDerivAt P
      (2400 * t ^ 7 + 5250 * t ^ 6 + 3300 * t ^ 5 + 375 * t ^ 4 -
        60 * t ^ 3 + 15 * t ^ 2 + 14 * t + 6) t := by
    have h0 := (hasDerivAt_pow 8 t).const_mul 300
    have h1 := h0.add ((hasDerivAt_pow 7 t).const_mul 750)
    have h2 := h1.add ((hasDerivAt_pow 6 t).const_mul 550)
    have h3 := h2.add ((hasDerivAt_pow 5 t).const_mul 75)
    have h4 := h3.sub ((hasDerivAt_pow 4 t).const_mul 15)
    have h5 := h4.add ((hasDerivAt_pow 3 t).const_mul 5)
    have h6 := h5.add ((hasDerivAt_pow 2 t).const_mul 7)
    have h7 := h6.add ((hasDerivAt_id t).const_mul 6)
    have hraw := h7.add_const 2
    convert hraw using 1
    · norm_num
      ring
  have hden := (((hasDerivAt_const t 1).add (hasDerivAt_id t)).pow 3).const_mul 32
  have hquot := (hP.const_mul (-75)).div hden
    (mul_ne_zero (by norm_num) (pow_ne_zero 3 ht))
  unfold denA denAPrime
  dsimp [P] at hquot
  convert hquot using 1
  field_simp [ht]
  ring

private theorem denB_hasDerivAt (t : ℝ) :
    HasDerivAt denB (denBPrime t) t := by
  unfold denB denBPrime
  convert (hasDerivAt_pow 6 t).const_mul (5625 / 8) using 1 <;>
    norm_num
  ring

private theorem logRatio_hasDerivAt (t : ℝ) (ht0 : t ≠ 0)
    (ht1 : 1 + t ≠ 0) :
    HasDerivAt logRatio (-1 / (t * (1 + t))) t := by
  have harg := (hasDerivAt_const t 1).add (hasDerivAt_id t)
  have h := (Real.hasDerivAt_log ht1).comp t harg
  have h' := h.sub (Real.hasDerivAt_log ht0)
  unfold logRatio
  convert h' using 1
  field_simp [ht0, ht1]
  ring

private theorem sqrtOneSubSq_hasDerivAt (t : ℝ) (ht : 1 - t ^ 2 ≠ 0) :
    HasDerivAt (fun x : ℝ => Real.sqrt (1 - x ^ 2))
      (-t / Real.sqrt (1 - t ^ 2)) t := by
  have hinner : HasDerivAt (fun x : ℝ => 1 - x ^ 2) (-2 * t) t := by
    convert (hasDerivAt_const t 1).sub ((hasDerivAt_id t).pow 2) using 1 <;>
      norm_num [id_eq]
  have h := (Real.hasDerivAt_sqrt ht).comp t hinner
  convert h using 1
  ring

private theorem denPrimitive_hasDerivAt (t : ℝ) (ht : t ∈ Ioo (0 : ℝ) 1) :
    HasDerivAt denPrimitive (denReducedIntegrand t) t := by
  have ht0 : t ≠ 0 := ht.1.ne'
  have ht1 : 1 + t ≠ 0 := by linarith [ht.1]
  have hspos : 0 < 1 - t ^ 2 := by
    nlinarith [ht.1, ht.2, sq_nonneg t, sq_nonneg (1 - t)]
  have hsqrt := sqrtOneSubSq_hasDerivAt t hspos.ne'
  have hinside := (denA_hasDerivAt t ht1).add
    ((denB_hasDerivAt t).mul (logRatio_hasDerivAt t ht0 ht1))
  have h := hsqrt.mul hinside
  have hroot : 0 < Real.sqrt (1 - t ^ 2) := Real.sqrt_pos.2 hspos
  unfold denPrimitive denReducedIntegrand
  convert h using 1
  simp only [Pi.add_apply, Pi.mul_apply]
  unfold denA denAPrime denB denBPrime denR denS
  field_simp [ht0, ht1, hroot.ne']
  rw [Real.sq_sqrt hspos.le]
  ring

private theorem denB_mul_log_continuous :
    Continuous (fun t : ℝ => denB t * Real.log t) := by
  have h : Continuous (fun t : ℝ =>
      (5625 / 8 : ℝ) * t ^ 5 * (t * Real.log t)) :=
    (continuous_const.mul (continuous_id.pow 5)).mul Real.continuous_mul_log
  convert h using 1
  funext t
  unfold denB
  ring

private theorem denS_mul_log_continuousOn :
    ContinuousOn (fun t : ℝ => denS t * Real.log t) (Icc 0 1) := by
  let C : ℝ → ℝ := fun t =>
    (-375 : ℝ) / 32 * t ^ 4 *
        (420 * t ^ 5 + 1260 * t ^ 4 + 900 * t ^ 3 - 660 * t ^ 2 -
          1080 * t - 360) /
      (1 + t) ^ 3
  have hC : ContinuousOn C (Icc (0 : ℝ) 1) := by
    dsimp [C]
    apply ContinuousOn.div
    · fun_prop
    · fun_prop
    · intro t ht
      have ht1 : 1 + t ≠ 0 := by linarith [ht.1]
      exact pow_ne_zero 3 ht1
  have hprod := hC.mul Real.continuous_mul_log.continuousOn
  apply hprod.congr
  intro t ht
  have ht1 : 1 + t ≠ 0 := by linarith [ht.1]
  dsimp [C]
  unfold denS
  field_simp [ht1]

private theorem denPrimitive_continuousOn :
    ContinuousOn denPrimitive (Icc (0 : ℝ) 1) := by
  have hA : ContinuousOn denA (Icc (0 : ℝ) 1) := by
    unfold denA
    apply ContinuousOn.div
    · fun_prop
    · fun_prop
    · intro t ht
      have ht1 : 1 + t ≠ 0 := by linarith [ht.1]
      exact mul_ne_zero (by norm_num) (pow_ne_zero 3 ht1)
  have hB : ContinuousOn denB (Icc (0 : ℝ) 1) := by
    unfold denB
    fun_prop
  have hlogOne : ContinuousOn (fun t : ℝ => Real.log (1 + t)) (Icc 0 1) := by
    apply ContinuousOn.log (continuousOn_const.add continuousOn_id)
    intro t ht
    have ht1 : 1 + t ≠ 0 := by linarith [ht.1]
    exact ht1
  have hBlogOne : ContinuousOn (fun t : ℝ => denB t * Real.log (1 + t))
      (Icc 0 1) := hB.mul hlogOne
  have hBlog : ContinuousOn (fun t : ℝ => denB t * logRatio t) (Icc 0 1) := by
    have hdiff := hBlogOne.sub denB_mul_log_continuous.continuousOn
    convert hdiff using 1
    funext t
    unfold logRatio
    ring
  unfold denPrimitive
  exact (Real.continuous_sqrt.comp
      (continuous_const.sub (continuous_id.pow 2))).continuousOn.mul (hA.add hBlog)

theorem denReducedIntegrand_intervalIntegrable :
    IntervalIntegrable denReducedIntegrand volume 0 1 := by
  let w : ℝ → ℝ := fun t => (Real.sqrt (1 - t ^ 2))⁻¹
  have hw : IntervalIntegrable w volume 0 1 := by
    have h := Polynomial.Chebyshev.intervalIntegrable_sqrt_one_sub_sq_inv.mono_set
      (c := (0 : ℝ)) (d := 1) (by
        intro t ht
        simp only [uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1), mem_Icc] at ht
        simp only [uIcc_of_le (by norm_num : (-1 : ℝ) ≤ 1), mem_Icc]
        constructor <;> linarith)
    simpa [w] using h
  have hR : ContinuousOn denR (Icc (0 : ℝ) 1) := by
    unfold denR
    apply ContinuousOn.div
    · fun_prop
    · fun_prop
    · intro t ht
      exact pow_ne_zero 3 (by linarith [ht.1])
  have hS : ContinuousOn denS (Icc (0 : ℝ) 1) := by
    unfold denS
    apply ContinuousOn.div
    · fun_prop
    · fun_prop
    · intro t ht
      exact pow_ne_zero 3 (by linarith [ht.1])
  have hlogOne : ContinuousOn (fun t : ℝ => Real.log (1 + t)) (Icc 0 1) := by
    apply ContinuousOn.log (continuousOn_const.add continuousOn_id)
    intro t ht
    exact (by linarith [ht.1] : 1 + t ≠ 0)
  have hSlogOne : ContinuousOn (fun t : ℝ => denS t * Real.log (1 + t))
      (Icc 0 1) := hS.mul hlogOne
  have hSlog : ContinuousOn (fun t : ℝ => denS t * logRatio t) (Icc 0 1) := by
    have hdiff := hSlogOne.sub denS_mul_log_continuousOn
    convert hdiff using 1
    unfold logRatio
    ring
  have hnum := hR.add hSlog
  have hnum' : ContinuousOn (denR + fun t : ℝ => denS t * logRatio t)
      [[(0 : ℝ), 1]] := by
    simpa [uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hnum
  have hint := hw.continuousOn_mul hnum'
  simpa [denReducedIntegrand, w, div_eq_mul_inv] using hint

theorem denReducedIntegral :
    (∫ t in (0 : ℝ)..1, denReducedIntegrand t) = 75 / 16 := by
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le (by norm_num)
    denPrimitive_continuousOn denPrimitive_hasDerivAt
    denReducedIntegrand_intervalIntegrable]
  norm_num [denPrimitive, denA, denB, logRatio]

private def numA (t : ℝ) : ℝ :=
  -t * (123780 * t ^ 7 + 309210 * t ^ 6 + 226330 * t ^ 5 +
      30205 * t ^ 4 - 6999 * t ^ 3 + 1075 * t ^ 2 + 1125 * t + 450) /
    (192 * (1 + t) ^ 3)

private def numAPrime (t : ℝ) : ℝ :=
  -(309450 * t ^ 8 + 1113540 * t ^ 7 + 1421730 * t ^ 6 +
      709195 * t ^ 5 + 72013 * t ^ 4 - 13998 * t ^ 3 +
      1050 * t ^ 2 + 675 * t + 225) /
    (96 * (1 + t) ^ 4)

private def numB (t : ℝ) : ℝ :=
  5 / 32 * t * (4126 * t ^ 5 - 8 * t ^ 4 - 10 * t ^ 2 - 15)

private def numBPrime (t : ℝ) : ℝ :=
  5 / 32 * (24756 * t ^ 5 - 40 * t ^ 4 - 30 * t ^ 2 - 15)

def numR (t : ℝ) : ℝ :=
  (-1 : ℝ) / 192 * t ^ 3 *
      (-866460 * t ^ 6 - 2164710 * t ^ 5 - 842230 * t ^ 4 +
        1642725 * t ^ 3 + 1405263 * t ^ 2 + 171157 * t - 30921) /
    (1 + t) ^ 3

def numS (t : ℝ) : ℝ :=
  (-1 : ℝ) / 192 * t ^ 3 *
      (866460 * t ^ 7 + 2597940 * t ^ 6 + 1852380 * t ^ 5 -
        1365900 * t ^ 4 - 2229480 * t ^ 3 - 742680 * t ^ 2) /
    (1 + t) ^ 3

def numReducedIntegrand (t : ℝ) : ℝ :=
  (numR t + numS t * logRatio t) / Real.sqrt (1 - t ^ 2)

private def numPrimitive (t : ℝ) : ℝ :=
  Real.sqrt (1 - t ^ 2) * (numA t + numB t * logRatio t)

private theorem numB_operator (t : ℝ) (ht : 1 + t ≠ 0) :
    75 / 32 + (1 - t ^ 2) * numBPrime t - t * numB t = numS t := by
  unfold numBPrime numB numS
  field_simp [ht]
  ring

private theorem numA_operator (t : ℝ) (ht0 : t ≠ 0) (ht1 : 1 + t ≠ 0) :
    (1 - t ^ 2) * numAPrime t - t * numA t -
        (1 - t) * numB t / t = numR t := by
  unfold numAPrime numA numB numR
  field_simp [ht0, ht1]
  ring

private theorem numA_hasDerivAt (t : ℝ) (ht : 1 + t ≠ 0) :
    HasDerivAt numA (numAPrime t) t := by
  let P : ℝ → ℝ := fun x =>
    123780 * x ^ 7 + 309210 * x ^ 6 + 226330 * x ^ 5 +
      30205 * x ^ 4 - 6999 * x ^ 3 + 1075 * x ^ 2 + 1125 * x + 450
  have h0 := (hasDerivAt_pow 7 t).const_mul 123780
  have h1 := h0.add ((hasDerivAt_pow 6 t).const_mul 309210)
  have h2 := h1.add ((hasDerivAt_pow 5 t).const_mul 226330)
  have h3 := h2.add ((hasDerivAt_pow 4 t).const_mul 30205)
  have h4 := h3.sub ((hasDerivAt_pow 3 t).const_mul 6999)
  have h5 := h4.add ((hasDerivAt_pow 2 t).const_mul 1075)
  have h6 := h5.add ((hasDerivAt_id t).const_mul 1125)
  have hP := h6.add_const 450
  have hnum := (hasDerivAt_id t).neg.mul hP
  have hden := (((hasDerivAt_const t 1).add (hasDerivAt_id t)).pow 3).const_mul 192
  have hquot := hnum.div hden
    (mul_ne_zero (by norm_num) (pow_ne_zero 3 ht))
  unfold numA numAPrime
  dsimp [P] at hquot
  convert hquot using 1
  field_simp [ht]
  ring

private theorem numB_hasDerivAt (t : ℝ) :
    HasDerivAt numB (numBPrime t) t := by
  have hraw := (((((hasDerivAt_pow 5 t).const_mul 4126).sub
      ((hasDerivAt_pow 4 t).const_mul 8)).sub
      ((hasDerivAt_pow 2 t).const_mul 10)).sub_const 15)
  have hP : HasDerivAt
      (fun x : ℝ => 4126 * x ^ 5 - 8 * x ^ 4 - 10 * x ^ 2 - 15)
      (20630 * t ^ 4 - 32 * t ^ 3 - 20 * t) t := by
    convert hraw using 1
    norm_num
    ring
  have h := ((hasDerivAt_id t).mul hP).const_mul (5 / 32)
  unfold numB numBPrime
  convert h using 1
  · funext x
    simp only [Pi.mul_apply, id_eq]
    ring
  · norm_num [id_eq]
    ring

private theorem numPrimitive_hasDerivAt (t : ℝ) (ht : t ∈ Ioo (0 : ℝ) 1) :
    HasDerivAt numPrimitive
      (numReducedIntegrand t - 75 / 32 * catalanKernelIntegrand t) t := by
  have ht0 : t ≠ 0 := ht.1.ne'
  have ht1 : 1 + t ≠ 0 := by linarith [ht.1]
  have hspos : 0 < 1 - t ^ 2 := by
    nlinarith [ht.1, ht.2, sq_nonneg t, sq_nonneg (1 - t)]
  have hsqrt := sqrtOneSubSq_hasDerivAt t hspos.ne'
  have hinside := (numA_hasDerivAt t ht1).add
    ((numB_hasDerivAt t).mul (logRatio_hasDerivAt t ht0 ht1))
  have h := hsqrt.mul hinside
  have hroot : 0 < Real.sqrt (1 - t ^ 2) := Real.sqrt_pos.2 hspos
  have hkernel : catalanKernelIntegrand t =
      logRatio t / Real.sqrt (1 - t ^ 2) := by
    unfold catalanKernelIntegrand logRatio
    rw [Real.log_div ht1 ht0]
  rw [hkernel]
  unfold numPrimitive numReducedIntegrand
  convert h using 1
  simp only [Pi.add_apply, Pi.mul_apply]
  unfold numA numAPrime numB numBPrime numR numS
  field_simp [ht0, ht1, hroot.ne']
  rw [Real.sq_sqrt hspos.le]
  ring

private theorem numB_mul_log_continuous :
    Continuous (fun t : ℝ => numB t * Real.log t) := by
  have hpoly : Continuous (fun t : ℝ =>
      (5 / 32 : ℝ) * (4126 * t ^ 5 - 8 * t ^ 4 - 10 * t ^ 2 - 15)) := by
    fun_prop
  have h := hpoly.mul Real.continuous_mul_log
  convert h using 1
  funext t
  simp only [Pi.mul_apply]
  unfold numB
  ring

private theorem numS_mul_log_continuousOn :
    ContinuousOn (fun t : ℝ => numS t * Real.log t) (Icc 0 1) := by
  let C : ℝ → ℝ := fun t =>
    (-1 : ℝ) / 192 * t ^ 4 *
        (866460 * t ^ 5 + 2597940 * t ^ 4 + 1852380 * t ^ 3 -
          1365900 * t ^ 2 - 2229480 * t - 742680) /
      (1 + t) ^ 3
  have hC : ContinuousOn C (Icc (0 : ℝ) 1) := by
    dsimp [C]
    apply ContinuousOn.div
    · fun_prop
    · fun_prop
    · intro t ht
      exact pow_ne_zero 3 (by linarith [ht.1])
  have hprod := hC.mul Real.continuous_mul_log.continuousOn
  apply hprod.congr
  intro t ht
  have ht1 : 1 + t ≠ 0 := by linarith [ht.1]
  dsimp [C]
  unfold numS
  field_simp [ht1]

private theorem numPrimitive_continuousOn :
    ContinuousOn numPrimitive (Icc (0 : ℝ) 1) := by
  have hA : ContinuousOn numA (Icc (0 : ℝ) 1) := by
    unfold numA
    apply ContinuousOn.div
    · fun_prop
    · fun_prop
    · intro t ht
      exact mul_ne_zero (by norm_num) (pow_ne_zero 3 (by linarith [ht.1]))
  have hB : ContinuousOn numB (Icc (0 : ℝ) 1) := by
    unfold numB
    fun_prop
  have hlogOne : ContinuousOn (fun t : ℝ => Real.log (1 + t)) (Icc 0 1) := by
    apply ContinuousOn.log (continuousOn_const.add continuousOn_id)
    intro t ht
    exact (by linarith [ht.1] : 1 + t ≠ 0)
  have hBlogOne : ContinuousOn (fun t : ℝ => numB t * Real.log (1 + t))
      (Icc 0 1) := hB.mul hlogOne
  have hBlog : ContinuousOn (fun t : ℝ => numB t * logRatio t) (Icc 0 1) := by
    have hdiff := hBlogOne.sub numB_mul_log_continuous.continuousOn
    convert hdiff using 1
    funext t
    unfold logRatio
    ring
  unfold numPrimitive
  exact (Real.continuous_sqrt.comp
      (continuous_const.sub (continuous_id.pow 2))).continuousOn.mul (hA.add hBlog)

theorem numReducedIntegrand_intervalIntegrable :
    IntervalIntegrable numReducedIntegrand volume 0 1 := by
  let w : ℝ → ℝ := fun t => (Real.sqrt (1 - t ^ 2))⁻¹
  have hw : IntervalIntegrable w volume 0 1 := by
    have h := Polynomial.Chebyshev.intervalIntegrable_sqrt_one_sub_sq_inv.mono_set
      (c := (0 : ℝ)) (d := 1) (by
        intro t ht
        simp only [uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1), mem_Icc] at ht
        simp only [uIcc_of_le (by norm_num : (-1 : ℝ) ≤ 1), mem_Icc]
        constructor <;> linarith)
    simpa [w] using h
  have hR : ContinuousOn numR (Icc (0 : ℝ) 1) := by
    unfold numR
    apply ContinuousOn.div
    · fun_prop
    · fun_prop
    · intro t ht
      exact pow_ne_zero 3 (by linarith [ht.1])
  have hS : ContinuousOn numS (Icc (0 : ℝ) 1) := by
    unfold numS
    apply ContinuousOn.div
    · fun_prop
    · fun_prop
    · intro t ht
      exact pow_ne_zero 3 (by linarith [ht.1])
  have hlogOne : ContinuousOn (fun t : ℝ => Real.log (1 + t)) (Icc 0 1) := by
    apply ContinuousOn.log (continuousOn_const.add continuousOn_id)
    intro t ht
    exact (by linarith [ht.1] : 1 + t ≠ 0)
  have hSlogOne : ContinuousOn (fun t : ℝ => numS t * Real.log (1 + t))
      (Icc 0 1) := hS.mul hlogOne
  have hSlog : ContinuousOn (fun t : ℝ => numS t * logRatio t) (Icc 0 1) := by
    have hdiff := hSlogOne.sub numS_mul_log_continuousOn
    convert hdiff using 1
    unfold logRatio
    ring
  have hnum := hR.add hSlog
  have hnum' : ContinuousOn (numR + fun t : ℝ => numS t * logRatio t)
      [[(0 : ℝ), 1]] := by
    simpa [uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hnum
  have hint := hw.continuousOn_mul hnum'
  simpa [numReducedIntegrand, w, div_eq_mul_inv] using hint

theorem numReducedIntegral :
    (∫ t in (0 : ℝ)..1, numReducedIntegrand t) =
      75 / 16 * catalanConstant := by
  have hscaled : IntervalIntegrable
      (fun t : ℝ => 75 / 32 * catalanKernelIntegrand t) volume 0 1 :=
    catalanKernel_intervalIntegrable.const_mul (75 / 32)
  have hdiff := numReducedIntegrand_intervalIntegrable.sub hscaled
  have hzero :
      (∫ t in (0 : ℝ)..1,
        (numReducedIntegrand t - 75 / 32 * catalanKernelIntegrand t)) = 0 := by
    rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le (by norm_num)
      numPrimitive_continuousOn numPrimitive_hasDerivAt hdiff]
    norm_num [numPrimitive, numA, numB, logRatio]
  rw [intervalIntegral.integral_sub numReducedIntegrand_intervalIntegrable hscaled,
    intervalIntegral.integral_const_mul, catalanKernelIntegral] at hzero
  linarith
end RamanujanChallenge.P25

end

