import RamanujanChallenge.Problem25DualInitialKernel
import RamanujanChallenge.Problem25DualProductReduction
import Mathlib.Analysis.Calculus.Deriv.Polynomial

/-! The exact initial pairing for the positive adjoint moment solution. -/

noncomputable section

namespace RamanujanChallenge.P25

open MeasureTheory Set
open scoped Interval

private def polyOfCoeffs : List ℝ → Polynomial ℝ
  | [] => 0
  | c :: cs => Polynomial.C c + Polynomial.X * polyOfCoeffs cs

private def evalCoeffs : List ℝ → ℝ → ℝ
  | [], _ => 0
  | c :: cs, x => c + x * evalCoeffs cs x

private def derivEvalCoeffs : List ℝ → ℝ → ℝ
  | [], _ => 0
  | _ :: cs, x => evalCoeffs cs x + x * derivEvalCoeffs cs x

@[simp] private theorem polyOfCoeffs_eval (cs : List ℝ) (x : ℝ) :
    (polyOfCoeffs cs).eval x = evalCoeffs cs x := by
  induction cs with
  | nil => simp [polyOfCoeffs, evalCoeffs]
  | cons c cs ih => simp [polyOfCoeffs, evalCoeffs, ih]

@[simp] private theorem polyOfCoeffs_derivative_eval (cs : List ℝ) (x : ℝ) :
    (polyOfCoeffs cs).derivative.eval x = derivEvalCoeffs cs x := by
  induction cs with
  | nil => simp [polyOfCoeffs, derivEvalCoeffs]
  | cons c cs ih =>
      simp [polyOfCoeffs, derivEvalCoeffs, Polynomial.derivative_add,
        Polynomial.derivative_mul, ih]

private def denAPolyP (t : ℝ) : Polynomial ℝ :=
  polyOfCoeffs [
    214 * t ^ 9 + 960 * t ^ 8 + 198 * t ^ 7,
    650 * t ^ 8 + 4800 * t ^ 7 + 1350 * t ^ 6,
    250 * t ^ 7 + 9600 * t ^ 6 + 3600 * t ^ 5,
    -1150 * t ^ 6 + 9600 * t ^ 5 + 4800 * t ^ 4,
    -1625 * t ^ 5 + 4800 * t ^ 4 + 3300 * t ^ 3,
    -745 * t ^ 4 + 960 * t ^ 3 + 1020 * t ^ 2,
    -70 * t ^ 3 - 66 * t,
    10 * t ^ 2 + 216 * t - 90,
    -97 * t + 120,
    -45]

private def denAPolyQ (t : ℝ) : Polynomial ℝ :=
  polyOfCoeffs [
    420 * t ^ 9 - 360 * t ^ 7,
    2100 * t ^ 8 - 1800 * t ^ 6,
    4200 * t ^ 7 - 3600 * t ^ 5,
    4200 * t ^ 6 - 3600 * t ^ 4,
    2100 * t ^ 5 - 1800 * t ^ 3,
    420 * t ^ 4 - 360 * t ^ 2]

private def denAIntegrand (t a : ℝ) : ℝ :=
  1125 / 8 * a ^ 5 * t ^ 3 * (a - 1) ^ 2 *
      (15 * a ^ 2 + 58 * a * t + 63 * t ^ 2) /
    (a + t) ^ 6

private def denAPrimitive (t a : ℝ) : ℝ :=
  (-375 * t ^ 3 / 32) *
      ((denAPolyP t).eval a + (denAPolyQ t).eval a * Real.log (a + t)) /
    (a + t) ^ 5

private theorem log_add_hasDerivAt (t a : ℝ) (h : a + t ≠ 0) :
    HasDerivAt (fun x : ℝ => Real.log (x + t)) (a + t)⁻¹ a := by
  have harg := (hasDerivAt_id a).add_const t
  simpa [Function.comp_def, id_eq] using (Real.hasDerivAt_log h).comp a harg

private theorem denAPrimitive_hasDerivAt (t a : ℝ) (h : a + t ≠ 0) :
    HasDerivAt (denAPrimitive t) (denAIntegrand t a) a := by
  have hP := (denAPolyP t).hasDerivAt a
  have hQ := (denAPolyQ t).hasDerivAt a
  have hlog := log_add_hasDerivAt t a h
  have hnum := hP.add (hQ.mul hlog)
  have hden := ((hasDerivAt_id a).add_const t).pow 5
  have hquot := hnum.div hden (pow_ne_zero 5 h)
  have hfull := hquot.const_mul (-375 * t ^ 3 / 32)
  unfold denAPrimitive denAIntegrand
  convert hfull using 1
  · funext x
    simp only [Pi.add_apply, Pi.mul_apply, Pi.div_apply, Pi.pow_apply, id_eq]
    ring
  · simp [denAPolyP, denAPolyQ, evalCoeffs, derivEvalCoeffs, id_eq]
    field_simp [h]
    ring

private theorem denAPrimitive_boundary (t : ℝ) (ht : 0 < t) :
    denAPrimitive t 1 - denAPrimitive t 0 =
      denR t + denS t * logRatio t := by
  have ht0 : t ≠ 0 := ht.ne'
  have ht1 : 1 + t ≠ 0 := by positivity
  unfold denAPrimitive denR denS logRatio
  simp [denAPolyP, denAPolyQ, evalCoeffs]
  field_simp [ht0, ht1]
  ring

theorem denAIntegral (t : ℝ) (ht : 0 < t) :
    (∫ a in (0 : ℝ)..1, denAIntegrand t a) =
      denR t + denS t * logRatio t := by
  have hne : ∀ a ∈ Icc (0 : ℝ) 1, a + t ≠ 0 := by
    intro a ha
    nlinarith [ha.1, ht]
  have hcont : ContinuousOn (denAPrimitive t) (Icc (0 : ℝ) 1) := by
    intro a ha
    exact (denAPrimitive_hasDerivAt t a (hne a ha)).continuousAt.continuousWithinAt
  have hint : IntervalIntegrable (denAIntegrand t) volume 0 1 := by
    apply ContinuousOn.intervalIntegrable
    unfold denAIntegrand
    apply ContinuousOn.div
    · fun_prop
    · fun_prop
    · intro a ha
      exact pow_ne_zero 6 (hne a (by
        simpa [uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using ha))
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le (by norm_num)
    hcont (fun a ha => denAPrimitive_hasDerivAt t a
      (hne a ⟨ha.1.le, ha.2.le⟩)) hint]
  exact denAPrimitive_boundary t ht

private def numAPolyP (t : ℝ) : Polynomial ℝ :=
  polyOfCoeffs [
    444002 * t ^ 9 + 1979592 * t ^ 8 + 407754 * t ^ 7,
    1353550 * t ^ 8 + 9899400 * t ^ 7 + 2781450 * t ^ 6,
    540950 * t ^ 7 + 19802400 * t ^ 6 + 7419600 * t ^ 5,
    -2347250 * t ^ 6 + 19807200 * t ^ 5 + 9895200 * t ^ 4,
    -3339775 * t ^ 5 + 9907200 * t ^ 4 + 6804300 * t ^ 3,
    -1534415 * t ^ 4 + 1982880 * t ^ 3 + 2103540 * t ^ 2,
    -144410 * t ^ 3 + 240 * t ^ 2 - 135870 * t,
    20630 * t ^ 2 + 445080 * t - 185526,
    -199895 * t + 247368,
    -92763]

private def numAPolyQ (t : ℝ) : Polynomial ℝ :=
  polyOfCoeffs [
    866460 * t ^ 9 - 1440 * t ^ 8 - 742680 * t ^ 7,
    4332300 * t ^ 8 - 7200 * t ^ 7 - 3713400 * t ^ 6,
    8664600 * t ^ 7 - 14400 * t ^ 6 - 7426800 * t ^ 5,
    8664600 * t ^ 6 - 14400 * t ^ 5 - 7426800 * t ^ 4,
    4332300 * t ^ 5 - 7200 * t ^ 4 - 3713400 * t ^ 3,
    866460 * t ^ 4 - 1440 * t ^ 3 - 742680 * t ^ 2]

private def numAIntegrand (t a : ℝ) : ℝ :=
  1 / 16 * a ^ 5 * t ^ 3 * (a - 1) ^ 2 *
      (30921 * a ^ 2 + 119546 * a * t + 129825 * t ^ 2) /
    (a + t) ^ 6

private def numAPrimitive (t a : ℝ) : ℝ :=
  (-t ^ 3 / 192) *
      ((numAPolyP t).eval a + (numAPolyQ t).eval a * Real.log (a + t)) /
    (a + t) ^ 5

private theorem numAPrimitive_hasDerivAt (t a : ℝ) (h : a + t ≠ 0) :
    HasDerivAt (numAPrimitive t) (numAIntegrand t a) a := by
  have hP := (numAPolyP t).hasDerivAt a
  have hQ := (numAPolyQ t).hasDerivAt a
  have hlog := log_add_hasDerivAt t a h
  have hnum := hP.add (hQ.mul hlog)
  have hden := ((hasDerivAt_id a).add_const t).pow 5
  have hquot := hnum.div hden (pow_ne_zero 5 h)
  have hfull := hquot.const_mul (-t ^ 3 / 192)
  unfold numAPrimitive numAIntegrand
  convert hfull using 1
  · funext x
    simp only [Pi.add_apply, Pi.mul_apply, Pi.div_apply, Pi.pow_apply, id_eq]
    ring
  · simp [numAPolyP, numAPolyQ, evalCoeffs, derivEvalCoeffs, id_eq]
    field_simp [h]
    ring

private theorem numAPrimitive_boundary (t : ℝ) (ht : 0 < t) :
    numAPrimitive t 1 - numAPrimitive t 0 =
      numR t + numS t * logRatio t := by
  have ht0 : t ≠ 0 := ht.ne'
  have ht1 : 1 + t ≠ 0 := by positivity
  unfold numAPrimitive numR numS logRatio
  simp [numAPolyP, numAPolyQ, evalCoeffs]
  field_simp [ht0, ht1]
  ring

theorem numAIntegral (t : ℝ) (ht : 0 < t) :
    (∫ a in (0 : ℝ)..1, numAIntegrand t a) =
      numR t + numS t * logRatio t := by
  have hne : ∀ a ∈ Icc (0 : ℝ) 1, a + t ≠ 0 := by
    intro a ha
    nlinarith [ha.1, ht]
  have hcont : ContinuousOn (numAPrimitive t) (Icc (0 : ℝ) 1) := by
    intro a ha
    exact (numAPrimitive_hasDerivAt t a (hne a ha)).continuousAt.continuousWithinAt
  have hint : IntervalIntegrable (numAIntegrand t) volume 0 1 := by
    apply ContinuousOn.intervalIntegrable
    unfold numAIntegrand
    apply ContinuousOn.div
    · fun_prop
    · fun_prop
    · intro a ha
      exact pow_ne_zero 6 (hne a (by
        simpa [uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using ha))
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le (by norm_num)
    hcont (fun a ha => numAPrimitive_hasDerivAt t a
      (hne a ⟨ha.1.le, ha.2.le⟩)) hint]
  exact numAPrimitive_boundary t ht

private theorem unit_integral_eq_interval (f : ℝ → ℝ) :
    (∫ x, f x ∂unitMeasure) = ∫ x in (0 : ℝ)..1, f x := by
  rw [intervalIntegral.integral_of_le (by norm_num : (0 : ℝ) ≤ 1)]

private def denH (v a : ℝ) : ℝ :=
  16 * (33750 * v ^ 3 / (a * (1 + v ^ 2) + 2 * v) ^ 4 +
    126000 * v ^ 4 / (a * (1 + v ^ 2) + 2 * v) ^ 5 +
    180000 * v ^ 5 / (a * (1 + v ^ 2) + 2 * v) ^ 6)

private def numH (v a : ℝ) : ℝ :=
  16 * (30921 * v ^ 3 / (a * (1 + v ^ 2) + 2 * v) ^ 4 +
    115408 * v ^ 4 / (a * (1 + v ^ 2) + 2 * v) ^ 5 +
    164800 * v ^ 5 / (a * (1 + v ^ 2) + 2 * v) ^ 6)

private theorem quadraticMap_pos_lt_one {v : ℝ} (hv0 : 0 < v) (hv1 : v < 1) :
    0 < quadraticMap v ∧ quadraticMap v < 1 := by
  have hden : 0 < 1 + v ^ 2 := by positivity
  constructor
  · dsimp [quadraticMap]
    positivity
  · dsimp [quadraticMap]
    rw [div_lt_one hden]
    nlinarith [sq_pos_of_pos (sub_pos.mpr hv1)]

private theorem sqrt_one_sub_quadraticMap_sq {v : ℝ}
    (hv0 : 0 < v) (hv1 : v < 1) :
    Real.sqrt (1 - quadraticMap v ^ 2) =
      (1 - v ^ 2) / (1 + v ^ 2) := by
  have hden : 0 < 1 + v ^ 2 := by positivity
  have hsquare :
      1 - quadraticMap v ^ 2 = ((1 - v ^ 2) / (1 + v ^ 2)) ^ 2 := by
    dsimp [quadraticMap]
    field_simp [hden.ne']
    ring
  rw [hsquare, Real.sqrt_sq_eq_abs, abs_of_nonneg]
  exact div_nonneg (by nlinarith [sq_nonneg v]) hden.le

private theorem denH_transformed {v a : ℝ}
    (hv0 : 0 < v) (hv1 : v < 1) (ha0 : 0 ≤ a) :
    a ^ 5 * (1 - a) ^ 2 * denH v a =
      (16 * quadraticMapDeriv v / Real.sqrt (1 - quadraticMap v ^ 2)) *
        denAIntegrand (quadraticMap v) a := by
  have hden : 0 < 1 + v ^ 2 := by positivity
  have hvsub : 1 - v ^ 2 ≠ 0 := by nlinarith
  have hD : a * (1 + v ^ 2) + 2 * v ≠ 0 := by positivity
  have hat : a + 2 * v / (1 + v ^ 2) ≠ 0 := by positivity
  rw [sqrt_one_sub_quadraticMap_sq hv0 hv1]
  unfold denH denAIntegrand quadraticMap quadraticMapDeriv
  field_simp [hden.ne', hvsub, hD, hat]
  ring

private theorem numH_transformed {v a : ℝ}
    (hv0 : 0 < v) (hv1 : v < 1) (ha0 : 0 ≤ a) :
    a ^ 5 * (1 - a) ^ 2 * numH v a =
      (16 * quadraticMapDeriv v / Real.sqrt (1 - quadraticMap v ^ 2)) *
        numAIntegrand (quadraticMap v) a := by
  have hden : 0 < 1 + v ^ 2 := by positivity
  have hvsub : 1 - v ^ 2 ≠ 0 := by nlinarith
  have hD : a * (1 + v ^ 2) + 2 * v ≠ 0 := by positivity
  have hat : a + 2 * v / (1 + v ^ 2) ≠ 0 := by positivity
  rw [sqrt_one_sub_quadraticMap_sq hv0 hv1]
  unfold numH numAIntegrand quadraticMap quadraticMapDeriv
  field_simp [hden.ne', hvsub, hD, hat]
  ring

private theorem denSliceIntegral {v : ℝ} (hv0 : 0 < v) (hv1 : v < 1) :
    (∫ a, a ^ 5 * (1 - a) ^ 2 * denH v a ∂unitMeasure) =
      16 * denReducedIntegrand (quadraticMap v) * quadraticMapDeriv v := by
  rw [unit_integral_eq_interval]
  calc
    (∫ a in (0 : ℝ)..1, a ^ 5 * (1 - a) ^ 2 * denH v a) =
        ∫ a in (0 : ℝ)..1,
          (16 * quadraticMapDeriv v / Real.sqrt (1 - quadraticMap v ^ 2)) *
            denAIntegrand (quadraticMap v) a := by
      apply intervalIntegral.integral_congr
      intro a ha
      apply denH_transformed hv0 hv1
      simpa [uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using ha.1
    _ = (16 * quadraticMapDeriv v / Real.sqrt (1 - quadraticMap v ^ 2)) *
        ∫ a in (0 : ℝ)..1, denAIntegrand (quadraticMap v) a := by
      rw [intervalIntegral.integral_const_mul]
    _ = _ := by
      rw [denAIntegral (quadraticMap v) (quadraticMap_pos_lt_one hv0 hv1).1]
      unfold denReducedIntegrand
      ring

private theorem numSliceIntegral {v : ℝ} (hv0 : 0 < v) (hv1 : v < 1) :
    (∫ a, a ^ 5 * (1 - a) ^ 2 * numH v a ∂unitMeasure) =
      16 * numReducedIntegrand (quadraticMap v) * quadraticMapDeriv v := by
  rw [unit_integral_eq_interval]
  calc
    (∫ a in (0 : ℝ)..1, a ^ 5 * (1 - a) ^ 2 * numH v a) =
        ∫ a in (0 : ℝ)..1,
          (16 * quadraticMapDeriv v / Real.sqrt (1 - quadraticMap v ^ 2)) *
            numAIntegrand (quadraticMap v) a := by
      apply intervalIntegral.integral_congr
      intro a ha
      apply numH_transformed hv0 hv1
      simpa [uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using ha.1
    _ = (16 * quadraticMapDeriv v / Real.sqrt (1 - quadraticMap v ^ 2)) *
        ∫ a in (0 : ℝ)..1, numAIntegrand (quadraticMap v) a := by
      rw [intervalIntegral.integral_const_mul]
    _ = _ := by
      rw [numAIntegral (quadraticMap v) (quadraticMap_pos_lt_one hv0 hv1).1]
      unfold numReducedIntegrand
      ring

private theorem quadraticMap_integral (g : ℝ → ℝ) :
    (∫ v, (g ∘ quadraticMap) v * quadraticMapDeriv v ∂unitMeasure) =
      ∫ t in (0 : ℝ)..1, g t := by
  rw [unit_integral_eq_interval]
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
    (g := g) hmap hderiv hnonneg
  have hends : quadraticMap 0 = 0 ∧ quadraticMap 1 = 1 := by
    norm_num [quadraticMap]
  rwa [hends.1, hends.2] at hsub

private theorem denOuterIntegral :
    (∫ v, 16 * denReducedIntegrand (quadraticMap v) * quadraticMapDeriv v
        ∂unitMeasure) = 75 := by
  calc
    _ = 16 * ∫ v,
        (denReducedIntegrand ∘ quadraticMap) v * quadraticMapDeriv v
        ∂unitMeasure := by
      rw [← MeasureTheory.integral_const_mul]
      apply integral_congr_ae
      filter_upwards [] with v
      simp only [Function.comp_apply]
      ring
    _ = 16 * ∫ t in (0 : ℝ)..1, denReducedIntegrand t := by
      rw [quadraticMap_integral]
    _ = 75 := by rw [denReducedIntegral]; norm_num

private theorem numOuterIntegral :
    (∫ v, 16 * numReducedIntegrand (quadraticMap v) * quadraticMapDeriv v
        ∂unitMeasure) = 75 * catalanConstant := by
  calc
    _ = 16 * ∫ v,
        (numReducedIntegrand ∘ quadraticMap) v * quadraticMapDeriv v
        ∂unitMeasure := by
      rw [← MeasureTheory.integral_const_mul]
      apply integral_congr_ae
      filter_upwards [] with v
      simp only [Function.comp_apply]
      ring
    _ = 16 * ∫ t in (0 : ℝ)..1, numReducedIntegrand t := by
      rw [quadraticMap_integral]
    _ = 75 * catalanConstant := by rw [numReducedIntegral]; ring

private theorem denH_continuousOn {v : ℝ} (hv : 0 < v) :
    ContinuousOn (denH v) (Icc (0 : ℝ) 1) := by
  have hD : ∀ a ∈ Icc (0 : ℝ) 1,
      a * (1 + v ^ 2) + 2 * v ≠ 0 := by
    intro a ha
    have ha0 : 0 ≤ a := ha.1
    positivity
  have hterm (c : ℝ) (m k : ℕ) : ContinuousOn
      (fun a : ℝ => c * v ^ m / (a * (1 + v ^ 2) + 2 * v) ^ k)
      (Icc (0 : ℝ) 1) := by
    apply ContinuousOn.div
    · fun_prop
    · fun_prop
    · intro a ha
      exact pow_ne_zero k (hD a ha)
  unfold denH
  exact continuousOn_const.mul
    (((hterm 33750 3 4).add (hterm 126000 4 5)).add (hterm 180000 5 6))

private theorem numH_continuousOn {v : ℝ} (hv : 0 < v) :
    ContinuousOn (numH v) (Icc (0 : ℝ) 1) := by
  have hD : ∀ a ∈ Icc (0 : ℝ) 1,
      a * (1 + v ^ 2) + 2 * v ≠ 0 := by
    intro a ha
    have ha0 : 0 ≤ a := ha.1
    positivity
  have hterm (c : ℝ) (m k : ℕ) : ContinuousOn
      (fun a : ℝ => c * v ^ m / (a * (1 + v ^ 2) + 2 * v) ^ k)
      (Icc (0 : ℝ) 1) := by
    apply ContinuousOn.div
    · fun_prop
    · fun_prop
    · intro a ha
      exact pow_ne_zero k (hD a ha)
  unfold numH
  exact continuousOn_const.mul
    (((hterm 30921 3 4).add (hterm 115408 4 5)).add (hterm 164800 5 6))

private theorem denProductMomentReduction {v : ℝ} (hv : 0 < v) :
    (∫ p, ∫ q, p ^ 6 * q ^ 5 * (1 - q ^ 2) * denH v (p * q)
        ∂unitMeasure ∂unitMeasure) =
      ∫ a, a ^ 5 * (1 - a) ^ 2 * denH v a ∂unitMeasure := by
  apply product_moment_reduction
  obtain ⟨C, hC⟩ := isCompact_Icc.exists_bound_of_continuousOn
    (denH_continuousOn hv)
  apply triangleKernel_integrable_of_bounded (denH v)
    (by unfold denH; fun_prop) (max 0 C) (le_max_left _ _)
  intro u hu
  exact (hC u ⟨hu.1.le, hu.2⟩).trans (le_max_right _ _)

private theorem numProductMomentReduction {v : ℝ} (hv : 0 < v) :
    (∫ p, ∫ q, p ^ 6 * q ^ 5 * (1 - q ^ 2) * numH v (p * q)
        ∂unitMeasure ∂unitMeasure) =
      ∫ a, a ^ 5 * (1 - a) ^ 2 * numH v a ∂unitMeasure := by
  apply product_moment_reduction
  obtain ⟨C, hC⟩ := isCompact_Icc.exists_bound_of_continuousOn
    (numH_continuousOn hv)
  apply triangleKernel_integrable_of_bounded (numH v)
    (by unfold numH; fun_prop) (max 0 C) (le_max_left _ _)
  intro u hu
  exact (hC u ⟨hu.1.le, hu.2⟩).trans (le_max_right _ _)

private theorem cube_integrable_cycle (f : ℝ × (ℝ × ℝ) → ℝ)
    (hf : Integrable f cubeMeasure) :
    Integrable (fun x : ℝ × (ℝ × ℝ) => f (x.2.1, x.2.2, x.1))
      cubeMeasure := by
  let g : (ℝ × ℝ) × ℝ → ℝ := fun x => f (x.1.1, x.1.2, x.2)
  have hg : Integrable g ((unitMeasure.prod unitMeasure).prod unitMeasure) := by
    have h := (MeasureTheory.measurePreserving_prodAssoc
      unitMeasure unitMeasure unitMeasure).integrable_comp_of_integrable hf
    simpa [g, cubeMeasure, Function.comp_def] using h
  have hswap := (MeasureTheory.Measure.measurePreserving_swap
    (μ := unitMeasure) (ν := unitMeasure.prod unitMeasure)).integrable_comp_of_integrable hg
  simpa [g, cubeMeasure, Function.comp_def] using hswap

private def denRawCombination (x : ℝ × (ℝ × ℝ)) : ℝ :=
  33750 * rawMomentIntegrand 0 0 0 0 0 x +
    126000 * rawMomentIntegrand 0 0 0 1 1 x +
    180000 * rawMomentIntegrand 0 0 0 2 2 x

private theorem denRawCombination_point (p q v : ℝ) :
    denRawCombination (p, q, v) =
      p ^ 6 * q ^ 5 * (1 - q ^ 2) * denH v (p * q) := by
  simp [denRawCombination, rawMomentIntegrand, denH, dualD]
  ring

private theorem denRawCombination_integrable :
    Integrable denRawCombination cubeMeasure := by
  have h0 := rawMomentIntegrand_integrable 0 0 0 0 0 (by omega)
  have h1 := rawMomentIntegrand_integrable 0 0 0 1 1 (by omega)
  have h2 := rawMomentIntegrand_integrable 0 0 0 2 2 (by omega)
  simpa [denRawCombination] using
    ((h0.const_mul 33750).add (h1.const_mul 126000)).add (h2.const_mul 180000)

private def denCycled (x : ℝ × (ℝ × ℝ)) : ℝ :=
  x.2.1 ^ 6 * x.2.2 ^ 5 * (1 - x.2.2 ^ 2) *
    denH x.1 (x.2.1 * x.2.2)

private theorem denCycled_integrable : Integrable denCycled cubeMeasure := by
  have h := cube_integrable_cycle denRawCombination denRawCombination_integrable
  simpa only [denCycled, denRawCombination_point] using h

private theorem denRawCombination_outer :
    (∫ x, denRawCombination x ∂cubeMeasure) = 75 := by
  calc
    (∫ x, denRawCombination x ∂cubeMeasure) =
        ∫ x, denRawCombination (x.2.1, x.2.2, x.1) ∂cubeMeasure :=
      cube_integral_cycle denRawCombination
    _ = ∫ x, denCycled x ∂cubeMeasure := by
      apply integral_congr_ae
      filter_upwards [] with x
      exact denRawCombination_point x.2.1 x.2.2 x.1
    _ = ∫ v, ∫ y, denCycled (v, y) ∂(unitMeasure.prod unitMeasure)
          ∂unitMeasure := by
      rw [MeasureTheory.integral_prod denCycled denCycled_integrable]
    _ = ∫ v, ∫ p, ∫ q,
          p ^ 6 * q ^ 5 * (1 - q ^ 2) * denH v (p * q)
          ∂unitMeasure ∂unitMeasure ∂unitMeasure := by
      apply integral_congr_ae
      filter_upwards [denCycled_integrable.prod_right_ae] with v hv
      rw [MeasureTheory.integral_prod _ hv]
      rfl
    _ = ∫ v, ∫ a, a ^ 5 * (1 - a) ^ 2 * denH v a
          ∂unitMeasure ∂unitMeasure := by
      apply integral_congr_ae
      filter_upwards [unit_ae_bounds] with v hv
      exact denProductMomentReduction hv.1
    _ = ∫ v, 16 * denReducedIntegrand (quadraticMap v) * quadraticMapDeriv v
          ∂unitMeasure := by
      apply integral_congr_ae
      filter_upwards [unit_ae_bounds,
        MeasureTheory.Measure.ae_ne unitMeasure (1 : ℝ)] with v hv hvne
      exact denSliceIntegral hv.1 (lt_of_le_of_ne hv.2 hvne)
    _ = 75 := denOuterIntegral

theorem denInitialMomentIdentity :
    33750 * dualMoment 0 0 0 0 0 +
      126000 * dualMoment 0 0 0 1 1 +
      180000 * dualMoment 0 0 0 2 2 = 75 := by
  have h0 := rawMomentIntegrand_integrable 0 0 0 0 0 (by omega)
  have h1 := rawMomentIntegrand_integrable 0 0 0 1 1 (by omega)
  have h2 := rawMomentIntegrand_integrable 0 0 0 2 2 (by omega)
  have hexpand : (∫ x, denRawCombination x ∂cubeMeasure) =
      33750 * dualMoment 0 0 0 0 0 +
        126000 * dualMoment 0 0 0 1 1 +
        180000 * dualMoment 0 0 0 2 2 := by
    have hsum01 := MeasureTheory.integral_add
      (h0.const_mul 33750) (h1.const_mul 126000)
    have hsum012 := MeasureTheory.integral_add
      ((h0.const_mul 33750).add (h1.const_mul 126000)) (h2.const_mul 180000)
    calc
      (∫ x, denRawCombination x ∂cubeMeasure) =
          (∫ x, 33750 * rawMomentIntegrand 0 0 0 0 0 x +
            126000 * rawMomentIntegrand 0 0 0 1 1 x ∂cubeMeasure) +
          ∫ x, 180000 * rawMomentIntegrand 0 0 0 2 2 x ∂cubeMeasure := by
        simpa only [denRawCombination, Pi.add_apply] using hsum012
      _ = ((∫ x, 33750 * rawMomentIntegrand 0 0 0 0 0 x ∂cubeMeasure) +
          ∫ x, 126000 * rawMomentIntegrand 0 0 0 1 1 x ∂cubeMeasure) +
          ∫ x, 180000 * rawMomentIntegrand 0 0 0 2 2 x ∂cubeMeasure := by
        rw [hsum01]
      _ = _ := by
        rw [MeasureTheory.integral_const_mul, MeasureTheory.integral_const_mul,
          MeasureTheory.integral_const_mul]
        rfl
  rw [← hexpand]
  exact denRawCombination_outer

private def numRawCombination (x : ℝ × (ℝ × ℝ)) : ℝ :=
  30921 * rawMomentIntegrand 0 0 0 0 0 x +
    115408 * rawMomentIntegrand 0 0 0 1 1 x +
    164800 * rawMomentIntegrand 0 0 0 2 2 x

private theorem numRawCombination_point (p q v : ℝ) :
    numRawCombination (p, q, v) =
      p ^ 6 * q ^ 5 * (1 - q ^ 2) * numH v (p * q) := by
  simp [numRawCombination, rawMomentIntegrand, numH, dualD]
  ring

private theorem numRawCombination_integrable :
    Integrable numRawCombination cubeMeasure := by
  have h0 := rawMomentIntegrand_integrable 0 0 0 0 0 (by omega)
  have h1 := rawMomentIntegrand_integrable 0 0 0 1 1 (by omega)
  have h2 := rawMomentIntegrand_integrable 0 0 0 2 2 (by omega)
  simpa [numRawCombination] using
    ((h0.const_mul 30921).add (h1.const_mul 115408)).add (h2.const_mul 164800)

private def numCycled (x : ℝ × (ℝ × ℝ)) : ℝ :=
  x.2.1 ^ 6 * x.2.2 ^ 5 * (1 - x.2.2 ^ 2) *
    numH x.1 (x.2.1 * x.2.2)

private theorem numCycled_integrable : Integrable numCycled cubeMeasure := by
  have h := cube_integrable_cycle numRawCombination numRawCombination_integrable
  simpa only [numCycled, numRawCombination_point] using h

private theorem numRawCombination_outer :
    (∫ x, numRawCombination x ∂cubeMeasure) = 75 * catalanConstant := by
  calc
    (∫ x, numRawCombination x ∂cubeMeasure) =
        ∫ x, numRawCombination (x.2.1, x.2.2, x.1) ∂cubeMeasure :=
      cube_integral_cycle numRawCombination
    _ = ∫ x, numCycled x ∂cubeMeasure := by
      apply integral_congr_ae
      filter_upwards [] with x
      exact numRawCombination_point x.2.1 x.2.2 x.1
    _ = ∫ v, ∫ y, numCycled (v, y) ∂(unitMeasure.prod unitMeasure)
          ∂unitMeasure := by
      rw [MeasureTheory.integral_prod numCycled numCycled_integrable]
    _ = ∫ v, ∫ p, ∫ q,
          p ^ 6 * q ^ 5 * (1 - q ^ 2) * numH v (p * q)
          ∂unitMeasure ∂unitMeasure ∂unitMeasure := by
      apply integral_congr_ae
      filter_upwards [numCycled_integrable.prod_right_ae] with v hv
      rw [MeasureTheory.integral_prod _ hv]
      rfl
    _ = ∫ v, ∫ a, a ^ 5 * (1 - a) ^ 2 * numH v a
          ∂unitMeasure ∂unitMeasure := by
      apply integral_congr_ae
      filter_upwards [unit_ae_bounds] with v hv
      exact numProductMomentReduction hv.1
    _ = ∫ v, 16 * numReducedIntegrand (quadraticMap v) * quadraticMapDeriv v
          ∂unitMeasure := by
      apply integral_congr_ae
      filter_upwards [unit_ae_bounds,
        MeasureTheory.Measure.ae_ne unitMeasure (1 : ℝ)] with v hv hvne
      exact numSliceIntegral hv.1 (lt_of_le_of_ne hv.2 hvne)
    _ = 75 * catalanConstant := numOuterIntegral

theorem numInitialMomentIdentity :
    30921 * dualMoment 0 0 0 0 0 +
      115408 * dualMoment 0 0 0 1 1 +
      164800 * dualMoment 0 0 0 2 2 = 75 * catalanConstant := by
  have h0 := rawMomentIntegrand_integrable 0 0 0 0 0 (by omega)
  have h1 := rawMomentIntegrand_integrable 0 0 0 1 1 (by omega)
  have h2 := rawMomentIntegrand_integrable 0 0 0 2 2 (by omega)
  have hexpand : (∫ x, numRawCombination x ∂cubeMeasure) =
      30921 * dualMoment 0 0 0 0 0 +
        115408 * dualMoment 0 0 0 1 1 +
        164800 * dualMoment 0 0 0 2 2 := by
    have hsum01 := MeasureTheory.integral_add
      (h0.const_mul 30921) (h1.const_mul 115408)
    have hsum012 := MeasureTheory.integral_add
      ((h0.const_mul 30921).add (h1.const_mul 115408)) (h2.const_mul 164800)
    calc
      (∫ x, numRawCombination x ∂cubeMeasure) =
          (∫ x, 30921 * rawMomentIntegrand 0 0 0 0 0 x +
            115408 * rawMomentIntegrand 0 0 0 1 1 x ∂cubeMeasure) +
          ∫ x, 164800 * rawMomentIntegrand 0 0 0 2 2 x ∂cubeMeasure := by
        simpa only [numRawCombination, Pi.add_apply] using hsum012
      _ = ((∫ x, 30921 * rawMomentIntegrand 0 0 0 0 0 x ∂cubeMeasure) +
          ∫ x, 115408 * rawMomentIntegrand 0 0 0 1 1 x ∂cubeMeasure) +
          ∫ x, 164800 * rawMomentIntegrand 0 0 0 2 2 x ∂cubeMeasure := by
        rw [hsum01]
      _ = _ := by
        rw [MeasureTheory.integral_const_mul, MeasureTheory.integral_const_mul,
          MeasureTheory.integral_const_mul]
        rfl
  rw [← hexpand]
  exact numRawCombination_outer

theorem dualVector_initial_denominator_pair :
    ∑ j : Fin 3, (positiveDenominator 0 j : ℝ) * dualVector 0 j = 75 := by
  simp only [Fin.sum_univ_three]
  norm_num [positiveDenominator, denominator, approximants, initialMatrix,
    coordinateSign, dualVector, Matrix.cons_val_two]
  linarith [denInitialMomentIdentity]

theorem dualVector_initial_numerator_pair :
    ∑ j : Fin 3, (positiveNumerator 0 j : ℝ) * dualVector 0 j =
      75 * catalanConstant := by
  simp only [Fin.sum_univ_three]
  norm_num [positiveNumerator, numerator, approximants, initialMatrix,
    coordinateSign, dualVector, Matrix.cons_val_two]
  linarith [numInitialMomentIdentity]

theorem dualVector_initial_error_pair :
    ∑ j : Fin 3, positiveCatalanError 0 j * dualVector 0 j = 0 := by
  have hden := dualVector_initial_denominator_pair
  have hnum := dualVector_initial_numerator_pair
  simp only [Fin.sum_univ_three] at hden hnum ⊢
  rw [positiveCatalanError_eq, positiveCatalanError_eq,
    positiveCatalanError_eq]
  calc
    _ = catalanConstant *
          ((positiveDenominator 0 0 : ℝ) * dualVector 0 0 +
            (positiveDenominator 0 1 : ℝ) * dualVector 0 1 +
            (positiveDenominator 0 2 : ℝ) * dualVector 0 2) -
        ((positiveNumerator 0 0 : ℝ) * dualVector 0 0 +
          (positiveNumerator 0 1 : ℝ) * dualVector 0 1 +
          (positiveNumerator 0 2 : ℝ) * dualVector 0 2) := by ring
    _ = catalanConstant * 75 - 75 * catalanConstant := by rw [hden, hnum]
    _ = 0 := by ring

end RamanujanChallenge.P25

end
