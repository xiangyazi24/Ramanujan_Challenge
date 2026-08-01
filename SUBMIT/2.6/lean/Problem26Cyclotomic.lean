import RamanujanChallenge.Problem26WeightThree
import Mathlib.MeasureTheory.Integral.IntervalIntegral.IntegrationByParts

/-!
  The remaining cubic cyclotomic special value for Problem 2.6.

  The proof stays over the reals.  The substitution `u = x / (1 + x)`
  moves the cubic polynomial `1 + x + x²` to the reflection-symmetric
  polynomial `1 - u + u²` on `[0, 1/2]`.
-/

open Filter Set Topology
open scoped Interval

noncomputable section

namespace RamanujanChallenge.P26

private theorem logSquare_intervalIntegrable26 :
    IntervalIntegrable (fun x : ℝ => Real.log x ^ 2)
      MeasureTheory.volume 0 1 := by
  let f : ℝ → ℝ := fun x => x ^ 2
  let f' : ℝ → ℝ := fun x => 2 * x
  let g : ℝ → ℝ := fun x => Real.log x ^ 2
  have hf : ContinuousOn f [[(0 : ℝ), 1]] := by
    unfold f
    fun_prop
  have hff' :
      ∀ x ∈ Ioo (min (0 : ℝ) 1) (max (0 : ℝ) 1),
        HasDerivAt f (f' x) x := by
    intro x hx
    unfold f f'
    convert (hasDerivAt_id x).pow 2 using 1
    simp
  have hf' :
      ∀ x ∈ Ioo (min (0 : ℝ) 1) (max (0 : ℝ) 1),
        0 ≤ f' x := by
    intro x hx
    norm_num at hx
    unfold f'
    exact mul_nonneg (by norm_num) hx.1.le
  have hlog :
      IntervalIntegrable (fun x : ℝ => 8 * Real.log x)
        MeasureTheory.volume 0 1 :=
    intervalIntegral.intervalIntegrable_log'.const_mul 8
  have hpull :
      IntervalIntegrable (fun x : ℝ => (g ∘ f) x * f' x)
        MeasureTheory.volume 0 1 := by
    apply hlog.mono_fun
    · exact
        (((Real.measurable_log.comp (measurable_id.pow_const 2)).pow_const 2).mul
          (measurable_const.mul measurable_id)).aestronglyMeasurable
    · filter_upwards [
        MeasureTheory.ae_restrict_mem measurableSet_uIoc
      ] with x hx
      have hx' : x ∈ Ioc (0 : ℝ) 1 := by
        simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hx
      have hx0 : 0 < x := hx'.1
      have hx1 : x ≤ 1 := hx'.2
      have hbound : |Real.log x * x| ≤ 1 :=
        (Real.abs_log_mul_self_lt x hx0 hx1).le
      dsimp [f, f', g, Function.comp_apply]
      rw [Real.log_pow]
      norm_num only [Nat.cast_ofNat]
      have hxabs : |x| = x := abs_of_pos hx0
      rw [abs_mul, hxabs] at hbound
      simp only [abs_mul, abs_pow, hxabs,
        abs_of_nonneg (by norm_num : (0 : ℝ) ≤ 2),
        abs_of_nonneg (by norm_num : (0 : ℝ) ≤ 8)]
      nlinarith [abs_nonneg (Real.log x)]
  have hiff :=
    intervalIntegral.integrable_comp_mul_deriv_iff_of_deriv_nonneg
      (g := g) hf hff' hf'
  simpa [f, g] using hiff.mp hpull

private def logAtOneSlope26 (x : ℝ) : ℝ :=
  Function.update (fun y : ℝ => Real.log y / (1 - y)) 1 (-1) x

private theorem logAtOneSlope26_continuousOn :
    ContinuousOn logAtOneSlope26 (Icc (1 / 2 : ℝ) 1) := by
  intro x hx
  by_cases hxone : x = 1
  · subst x
    have hlog : HasDerivAt Real.log 1 1 := by
      simpa using Real.hasDerivAt_log (by norm_num : (1 : ℝ) ≠ 0)
    have hc := hlog.continuousAt_div
    have hc' : ContinuousAt logAtOneSlope26 1 := by
      convert hc.neg using 1
      funext y
      by_cases hy : y = 1
      · subst y
        simp [logAtOneSlope26]
      · have h1y : 1 - y ≠ 0 := sub_ne_zero.mpr (Ne.symm hy)
        have hy1 : y - 1 ≠ 0 := sub_ne_zero.mpr hy
        simp [logAtOneSlope26, hy]
        field_simp [h1y, hy1]
        ring
    exact hc'.continuousWithinAt
  · have hbase :
        ContinuousAt (fun y : ℝ => Real.log y / (1 - y)) x := by
      have hx0 : 0 < x := by linarith [hx.1]
      exact
          (Real.continuousAt_log (ne_of_gt hx0)).div
          (continuousAt_const.sub continuousAt_id)
          (sub_ne_zero.mpr (Ne.symm hxone))
    have heq :
        logAtOneSlope26 =ᶠ[𝓝 x]
          (fun y : ℝ => Real.log y / (1 - y)) := by
      filter_upwards [eventually_ne_nhds hxone] with y hy
      simp [logAtOneSlope26, hy]
    exact (hbase.congr_of_eventuallyEq heq).continuousWithinAt

private theorem logSquareOverOneAdd_intervalIntegrable26 :
    IntervalIntegrable (fun x : ℝ => Real.log x ^ 2 / (1 + x))
      MeasureTheory.volume 0 1 := by
  have hfactor :
      ContinuousOn (fun x : ℝ => 1 / (1 + x)) [[(0 : ℝ), 1]] := by
    apply ContinuousOn.div continuousOn_const
      (continuousOn_const.add continuousOn_id)
    intro x hx
    have hx' : x ∈ Icc (0 : ℝ) 1 := by
      simpa [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hx
    exact ne_of_gt (show 0 < 1 + x by linarith [hx'.1])
  have h :=
    logSquare_intervalIntegrable26.mul_continuousOn hfactor
  apply h.congr
  intro x _
  ring

private theorem logSquareOverOneSub_intervalIntegrable26 :
    IntervalIntegrable (fun x : ℝ => Real.log x ^ 2 / (1 - x))
      MeasureTheory.volume 0 1 := by
  apply IntervalIntegrable.trans (b := (1 / 2 : ℝ))
  · have hfactor :
        ContinuousOn (fun x : ℝ => 1 / (1 - x))
          [[(0 : ℝ), (1 / 2 : ℝ)]] := by
      apply ContinuousOn.div continuousOn_const
        (continuousOn_const.sub continuousOn_id)
      intro x hx
      have hx' : x ∈ Icc (0 : ℝ) (1 / 2) := by
        simpa [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2)] using hx
      exact ne_of_gt (show 0 < 1 - x by linarith [hx'.2])
    have hsquare :=
      logSquare_intervalIntegrable26.mono_set
        (show [[(0 : ℝ), (1 / 2 : ℝ)]] ⊆ [[(0 : ℝ), 1]] by
          rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2),
            Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
          intro x hx
          exact ⟨hx.1, hx.2.trans (by norm_num)⟩)
    have h := hsquare.mul_continuousOn hfactor
    apply h.congr
    intro x _
    ring
  · have hcont :
        ContinuousOn
          (fun x : ℝ => (1 - x) * logAtOneSlope26 x ^ 2)
          (Icc (1 / 2 : ℝ) 1) :=
      (continuousOn_const.sub continuousOn_id).mul
        (logAtOneSlope26_continuousOn.pow 2)
    have h :
        IntervalIntegrable
          (fun x : ℝ => (1 - x) * logAtOneSlope26 x ^ 2)
          MeasureTheory.volume (1 / 2) 1 := by
      apply ContinuousOn.intervalIntegrable
      rw [Set.uIcc_of_le (by norm_num : (1 / 2 : ℝ) ≤ 1)]
      exact hcont
    apply h.congr
    intro x hx
    have hx' : x ∈ Ioc (1 / 2 : ℝ) 1 := by
      rw [Set.uIoc_of_le (by norm_num : (1 / 2 : ℝ) ≤ 1)] at hx
      exact hx
    by_cases hxone : x = 1
    · subst x
      simp
    · simp [logAtOneSlope26, hxone]
      field_simp [sub_ne_zero.mpr hxone]

private def logSquareOverOneSubPrimitive26 (x : ℝ) : ℝ :=
  -(Real.log x ^ 2 * Real.log (1 - x)) -
    2 * Real.log x * dilog x + 2 * trilog26 x

private theorem logSquareOverOneSubPrimitive26_hasDerivAt
    {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasDerivAt logSquareOverOneSubPrimitive26
      (Real.log x ^ 2 / (1 - x)) x := by
  have hxne : x ≠ 0 := ne_of_gt hx0
  have h1xne : 1 - x ≠ 0 := ne_of_gt (sub_pos.mpr hx1)
  have hlogx : HasDerivAt Real.log (1 / x) x := by
    simpa [one_div] using Real.hasDerivAt_log hxne
  have hsub :
      HasDerivAt (fun y : ℝ => 1 - y) (-1) x := by
    convert (hasDerivAt_const x 1).sub (hasDerivAt_id x) using 1
    simp
  have hlog1 :
      HasDerivAt (fun y : ℝ => Real.log (1 - y))
        (-1 / (1 - x)) x := by
    convert hsub.log h1xne using 1
  have hd := dilog_hasDerivAt hx0 hx1
  have ht := trilog26_hasDerivAt_of_abs_lt_one
    (by rw [abs_of_pos hx0]; exact hx1) hxne
  unfold logSquareOverOneSubPrimitive26
  have htotal :=
    ((hlogx.pow 2).mul hlog1).neg.sub
      ((hlogx.mul hd).const_mul 2) |>.add (ht.const_mul 2)
  convert htotal using 1
  · funext y
    simp only [Pi.add_apply, Pi.sub_apply, Pi.mul_apply, Pi.pow_apply,
      Pi.neg_apply]
    ring
  · simp only [Pi.pow_apply]
    field_simp [hxne, h1xne]
    ring

private def logSquareOverOneAddPrimitive26 (x : ℝ) : ℝ :=
  Real.log x ^ 2 * Real.log (1 + x) +
    2 * Real.log x * dilog (-x) - 2 * trilog26 (-x)

private theorem logSquareOverOneAddPrimitive26_hasDerivAt
    {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasDerivAt logSquareOverOneAddPrimitive26
      (Real.log x ^ 2 / (1 + x)) x := by
  have hxne : x ≠ 0 := ne_of_gt hx0
  have h1xne : 1 + x ≠ 0 := by linarith
  have hlogx : HasDerivAt Real.log (1 / x) x := by
    simpa [one_div] using Real.hasDerivAt_log hxne
  have hadd :
      HasDerivAt (fun y : ℝ => 1 + y) 1 x := by
    convert (hasDerivAt_const x 1).add (hasDerivAt_id x) using 1
    simp
  have hlog1 :
      HasDerivAt (fun y : ℝ => Real.log (1 + y))
        (1 / (1 + x)) x := by
    convert hadd.log h1xne using 1
  have hneg :
      HasDerivAt (fun y : ℝ => -y) (-1) x := by
    simpa using (hasDerivAt_id x).neg
  have hd :
      HasDerivAt (fun y : ℝ => dilog (-y))
        (-(Real.log (1 + x) / x)) x := by
    convert
      (dilog_hasDerivAt_of_abs_lt_one
        (x := -x)
        (by rw [abs_neg, abs_of_pos hx0]; exact hx1)
        (neg_ne_zero.mpr hxne)).comp x hneg using 1
    ring
  have ht :
      HasDerivAt (fun y : ℝ => trilog26 (-y))
        (dilog (-x) / x) x := by
    convert
      (trilog26_hasDerivAt_of_abs_lt_one
        (x := -x)
        (by rw [abs_neg, abs_of_pos hx0]; exact hx1)
        (neg_ne_zero.mpr hxne)).comp x hneg using 1
    ring
  unfold logSquareOverOneAddPrimitive26
  have htotal :=
    ((hlogx.pow 2).mul hlog1).add
      ((hlogx.mul hd).const_mul 2) |>.sub (ht.const_mul 2)
  convert htotal using 1
  · funext y
    simp only [Pi.add_apply, Pi.sub_apply, Pi.mul_apply, Pi.pow_apply]
    ring
  · simp only [Pi.pow_apply]
    field_simp [hxne, h1xne]
    ring

private theorem logSquare_mul_self_tendsto_zero26 :
    Tendsto (fun x : ℝ => Real.log x ^ 2 * x)
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  have h :=
    tendsto_log_mul_rpow_nhdsGT_zero
      (show (0 : ℝ) < 1 / 2 by norm_num)
  have hsq := h.mul h
  have hsq0 :
      Tendsto
        (fun x : ℝ =>
          (Real.log x * x ^ (1 / 2 : ℝ)) *
            (Real.log x * x ^ (1 / 2 : ℝ)))
        (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    simpa using hsq
  apply Filter.Tendsto.congr' _ hsq0
  filter_upwards [self_mem_nhdsWithin] with x hx
  have hx0 : 0 < x := hx
  have hrpow :
      x ^ (1 / 2 : ℝ) * x ^ (1 / 2 : ℝ) = x := by
    rw [← Real.rpow_add hx0]
    norm_num
  calc
    (Real.log x * x ^ (1 / 2 : ℝ)) *
        (Real.log x * x ^ (1 / 2 : ℝ)) =
      Real.log x ^ 2 *
        (x ^ (1 / 2 : ℝ) * x ^ (1 / 2 : ℝ)) := by ring
    _ = Real.log x ^ 2 * x := by rw [hrpow]

private theorem log_mul_dilog_tendsto_zero26 :
    Tendsto (fun x : ℝ => Real.log x * dilog x)
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  have hslope :
      Tendsto (fun x : ℝ => x⁻¹ * dilog x)
        (𝓝[>] (0 : ℝ)) (𝓝 1) := by
    simpa [dilog_zero] using
      dilog_hasDerivAt_zero26.tendsto_slope_zero_right
  have hlogx :
      Tendsto (fun x : ℝ => Real.log x * x)
        (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    simpa [Real.rpow_one] using
      (tendsto_log_mul_rpow_nhdsGT_zero
        (show (0 : ℝ) < 1 by norm_num))
  have hraw :=
    hlogx.mul hslope
  have hraw0 :
      Tendsto
        (fun x : ℝ => (Real.log x * x) * (x⁻¹ * dilog x))
        (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    simpa using hraw
  apply Filter.Tendsto.congr' _ hraw0
  filter_upwards [self_mem_nhdsWithin] with x hx
  have hxne : x ≠ 0 := ne_of_gt hx
  field_simp [hxne]

private theorem log_mul_dilog_neg_tendsto_zero26 :
    Tendsto (fun x : ℝ => Real.log x * dilog (-x))
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  have hneg :
      HasDerivAt (fun x : ℝ => dilog (-x)) (-1) 0 := by
    have hinner :
        HasDerivAt (fun x : ℝ => -x) (-1) 0 := by
      simpa using (hasDerivAt_id (0 : ℝ)).neg
    have houter :
        HasDerivAt dilog 1 (-((0 : ℝ))) := by
      simpa using dilog_hasDerivAt_zero26
    convert houter.comp 0 hinner using 1 <;> norm_num
  have hslope :
      Tendsto (fun x : ℝ => x⁻¹ * dilog (-x))
        (𝓝[>] (0 : ℝ)) (𝓝 (-1)) := by
    simpa [dilog_zero] using hneg.tendsto_slope_zero_right
  have hlogx :
      Tendsto (fun x : ℝ => Real.log x * x)
        (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    simpa [Real.rpow_one] using
      (tendsto_log_mul_rpow_nhdsGT_zero
        (show (0 : ℝ) < 1 by norm_num))
  have hraw := hlogx.mul hslope
  have hraw0 :
      Tendsto
        (fun x : ℝ => (Real.log x * x) * (x⁻¹ * dilog (-x)))
        (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    simpa using hraw
  apply Filter.Tendsto.congr' _ hraw0
  filter_upwards [self_mem_nhdsWithin] with x hx
  have hxne : x ≠ 0 := ne_of_gt hx
  field_simp [hxne]

private theorem logSquare_mul_log_one_add_tendsto_zero26 :
    Tendsto
      (fun x : ℝ => Real.log x ^ 2 * Real.log (1 + x))
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  have hslope :
      Tendsto (fun x : ℝ => x⁻¹ * Real.log (1 + x))
        (𝓝[>] (0 : ℝ)) (𝓝 1) := by
    have hinner :
        HasDerivAt (fun x : ℝ => 1 + x) 1 0 := by
      convert (hasDerivAt_const (0 : ℝ) 1).add
        (hasDerivAt_id 0) using 1
      norm_num
    have hlog :
        HasDerivAt (fun x : ℝ => Real.log (1 + x)) 1 0 := by
      convert hinner.log (by norm_num) using 1 <;> norm_num
    simpa using hlog.tendsto_slope_zero_right
  have hraw := logSquare_mul_self_tendsto_zero26.mul hslope
  have hraw0 :
      Tendsto
        (fun x : ℝ =>
          (Real.log x ^ 2 * x) *
            (x⁻¹ * Real.log (1 + x)))
        (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    simpa using hraw
  apply Filter.Tendsto.congr' _ hraw0
  filter_upwards [self_mem_nhdsWithin] with x hx
  have hxne : x ≠ 0 := ne_of_gt hx
  field_simp [hxne]

private theorem logSquare_mul_log_one_sub_tendsto_zero26 :
    Tendsto
      (fun x : ℝ => Real.log x ^ 2 * Real.log (1 - x))
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  have hslope :
      Tendsto (fun x : ℝ => x⁻¹ * Real.log (1 - x))
        (𝓝[>] (0 : ℝ)) (𝓝 (-1)) := by
    have hinner :
        HasDerivAt (fun x : ℝ => 1 - x) (-1) 0 := by
      convert (hasDerivAt_const (0 : ℝ) 1).sub
        (hasDerivAt_id 0) using 1
      norm_num
    have hlog :
        HasDerivAt (fun x : ℝ => Real.log (1 - x)) (-1) 0 := by
      convert hinner.log (by norm_num) using 1 <;> norm_num
    simpa using hlog.tendsto_slope_zero_right
  have hraw := logSquare_mul_self_tendsto_zero26.mul hslope
  have hraw0 :
      Tendsto
        (fun x : ℝ =>
          (Real.log x ^ 2 * x) *
            (x⁻¹ * Real.log (1 - x)))
        (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    simpa using hraw
  apply Filter.Tendsto.congr' _ hraw0
  filter_upwards [self_mem_nhdsWithin] with x hx
  have hxne : x ≠ 0 := ne_of_gt hx
  field_simp [hxne]

private theorem oneSub_tendsto_nhdsLT_one_nhdsGT_zero26 :
    Tendsto (fun x : ℝ => 1 - x)
      (𝓝[<] (1 : ℝ)) (𝓝[>] (0 : ℝ)) := by
  rw [tendsto_nhdsWithin_iff]
  constructor
  · have hcont : ContinuousAt (fun x : ℝ => 1 - x) 1 := by
      fun_prop
    simpa using hcont.tendsto.mono_left
      (show (𝓝[<] (1 : ℝ)) ≤ 𝓝 1 from inf_le_left)
  · filter_upwards [self_mem_nhdsWithin] with x hx
    exact sub_pos.mpr (show x < 1 from hx)

private theorem logSquare_mul_log_one_sub_tendsto_one26 :
    Tendsto
      (fun x : ℝ => Real.log x ^ 2 * Real.log (1 - x))
      (𝓝[<] (1 : ℝ)) (𝓝 0) := by
  have hslope :
      Tendsto logAtOneSlope26
        (𝓝[<] (1 : ℝ)) (𝓝 (-1)) := by
    have hWithin :
        ContinuousWithinAt logAtOneSlope26
          (Iio (1 : ℝ)) 1 :=
      (logAtOneSlope26_continuousOn 1
        (by norm_num : (1 : ℝ) ∈ Icc (1 / 2) 1)).mono_of_mem_nhdsWithin
        (Icc_mem_nhdsLT (show (1 / 2 : ℝ) < 1 by norm_num))
    simpa [logAtOneSlope26] using hWithin.tendsto
  have hsub :
      Tendsto (fun x : ℝ => 1 - x)
        (𝓝[<] (1 : ℝ)) (𝓝 0) :=
    oneSub_tendsto_nhdsLT_one_nhdsGT_zero26.mono_right inf_le_left
  have hmulLog :
      Tendsto (fun x : ℝ => (1 - x) * Real.log (1 - x))
        (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    simpa using
      Real.continuous_mul_log.continuousAt.tendsto.comp hsub
  have hraw :=
    hmulLog.mul (hsub.mul (hslope.mul hslope))
  have hraw0 :
      Tendsto
        (fun x : ℝ =>
          ((1 - x) * Real.log (1 - x)) *
            ((1 - x) *
              (logAtOneSlope26 x * logAtOneSlope26 x)))
        (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    simpa using hraw
  apply Filter.Tendsto.congr' _ hraw0
  filter_upwards [
    Ioo_mem_nhdsLT (show (1 / 2 : ℝ) < 1 by norm_num)
  ] with x hx
  have hxone : x ≠ 1 := ne_of_lt hx.2
  simp [logAtOneSlope26, hxone]
  field_simp [sub_ne_zero.mpr (Ne.symm hxone)]

private theorem logSquareOverOneSubPrimitive26_tendsto_zero :
    Tendsto logSquareOverOneSubPrimitive26
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  have htri :
      Tendsto trilog26 (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    have hcont :
        ContinuousAt trilog26 0 :=
      trilog26_continuousOn_unit.continuousAt
        (Icc_mem_nhds (by norm_num : (-1 : ℝ) < 0)
          (by norm_num : (0 : ℝ) < 1))
    simpa using tendsto_nhdsWithin_of_tendsto_nhds hcont.tendsto
  unfold logSquareOverOneSubPrimitive26
  convert
    logSquare_mul_log_one_sub_tendsto_zero26.neg.sub
      (log_mul_dilog_tendsto_zero26.const_mul 2) |>.add
        (htri.const_mul 2) using 1
  · funext x
    ring
  · ring

private theorem logSquareOverOneAddPrimitive26_tendsto_zero :
    Tendsto logSquareOverOneAddPrimitive26
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  have htri :
      Tendsto (fun x : ℝ => trilog26 (-x))
        (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    have houter :
        ContinuousAt trilog26 (0 : ℝ) :=
      trilog26_continuousOn_unit.continuousAt
        (Icc_mem_nhds (by norm_num : (-1 : ℝ) < 0)
          (by norm_num : (0 : ℝ) < 1))
    have hneg :
        Tendsto (fun x : ℝ => -x)
          (𝓝[>] (0 : ℝ)) (𝓝 (0 : ℝ)) :=
      tendsto_nhdsWithin_of_tendsto_nhds
        (by
          simpa using
            (show ContinuousAt (fun x : ℝ => -x) 0 by
              simpa using
                (continuousAt_id.neg :
                  ContinuousAt (fun x : ℝ => -x) 0)).tendsto)
    simpa [Function.comp_def, trilog26_zero] using houter.tendsto.comp hneg
  unfold logSquareOverOneAddPrimitive26
  convert
    logSquare_mul_log_one_add_tendsto_zero26.add
      (log_mul_dilog_neg_tendsto_zero26.const_mul 2) |>.sub
        (htri.const_mul 2) using 1
  · funext x
    ring
  · ring

private theorem logSquareOverOneSubPrimitive26_tendsto_one :
    Tendsto logSquareOverOneSubPrimitive26
      (𝓝[<] (1 : ℝ)) (𝓝 (2 * zeta3)) := by
  have hlog :
      Tendsto Real.log (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    simpa using tendsto_nhdsWithin_of_tendsto_nhds
      (Real.continuousAt_log (by norm_num : (1 : ℝ) ≠ 0)).tendsto
  have hdWithin :
      ContinuousWithinAt dilog (Iio (1 : ℝ)) 1 :=
    (dilog_continuousOn_unit 1 (by norm_num)).mono_of_mem_nhdsWithin
      (Icc_mem_nhdsLT (show (-1 : ℝ) < 1 by norm_num))
  have hd :
      Tendsto dilog (𝓝[<] (1 : ℝ)) (𝓝 (dilog 1)) :=
    hdWithin.tendsto
  have hproduct :
      Tendsto (fun x : ℝ => Real.log x * dilog x)
        (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    simpa using hlog.mul hd
  have htWithin :
      ContinuousWithinAt trilog26 (Iio (1 : ℝ)) 1 :=
    (trilog26_continuousOn_unit 1 (by norm_num)).mono_of_mem_nhdsWithin
      (Icc_mem_nhdsLT (show (-1 : ℝ) < 1 by norm_num))
  have ht :
      Tendsto trilog26 (𝓝[<] (1 : ℝ)) (𝓝 zeta3) := by
    simpa [trilog26_one] using htWithin.tendsto
  unfold logSquareOverOneSubPrimitive26
  convert
    logSquare_mul_log_one_sub_tendsto_one26.neg.sub
      (hproduct.const_mul 2) |>.add (ht.const_mul 2) using 1
  · funext x
    ring
  · ring

private theorem logSquareOverOneAddPrimitive26_tendsto_one :
    Tendsto logSquareOverOneAddPrimitive26
      (𝓝[<] (1 : ℝ)) (𝓝 ((3 / 2 : ℝ) * zeta3)) := by
  have hlog :
      Tendsto Real.log (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    simpa using tendsto_nhdsWithin_of_tendsto_nhds
      (Real.continuousAt_log (by norm_num : (1 : ℝ) ≠ 0)).tendsto
  have hlogOneAdd :
      Tendsto (fun x : ℝ => Real.log (1 + x))
        (𝓝[<] (1 : ℝ)) (𝓝 (Real.log 2)) := by
    have hcont :
        ContinuousAt (fun x : ℝ => Real.log (1 + x)) 1 :=
      (continuousAt_const.add continuousAt_id).log
        (by norm_num : (1 + (1 : ℝ)) ≠ 0)
    have ht :
        Tendsto (fun x : ℝ => Real.log (1 + x))
          (𝓝[<] (1 : ℝ)) (𝓝 (Real.log (1 + 1))) :=
      hcont.tendsto.mono_left
        (show (𝓝[<] (1 : ℝ)) ≤ 𝓝 1 from inf_le_left)
    norm_num at ht ⊢
    exact ht
  have hfirst :
      Tendsto
        (fun x : ℝ => Real.log x ^ 2 * Real.log (1 + x))
        (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    simpa using (hlog.pow 2).mul hlogOneAdd
  have hDcont :
      ContinuousOn (fun x : ℝ => dilog (-x)) (Icc (0 : ℝ) 1) := by
    apply dilog_continuousOn_unit.comp
    · fun_prop
    · intro x hx
      constructor <;> linarith [hx.1, hx.2]
  have hd :
      Tendsto (fun x : ℝ => dilog (-x))
        (𝓝[<] (1 : ℝ)) (𝓝 (dilog (-1))) := by
    exact
      ((hDcont 1 (by norm_num)).mono_of_mem_nhdsWithin
        (Icc_mem_nhdsLT (show (0 : ℝ) < 1 by norm_num))).tendsto
  have hproduct :
      Tendsto (fun x : ℝ => Real.log x * dilog (-x))
        (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    simpa using hlog.mul hd
  have hTcont :
      ContinuousOn (fun x : ℝ => trilog26 (-x)) (Icc (0 : ℝ) 1) := by
    apply trilog26_continuousOn_unit.comp
    · fun_prop
    · intro x hx
      constructor <;> linarith [hx.1, hx.2]
  have ht :
      Tendsto (fun x : ℝ => trilog26 (-x))
        (𝓝[<] (1 : ℝ)) (𝓝 (-(3 : ℝ) / 4 * zeta3)) := by
    have hwithin :=
      (hTcont 1 (by norm_num)).mono_of_mem_nhdsWithin
        (Icc_mem_nhdsLT (show (0 : ℝ) < 1 by norm_num))
    simpa [trilog26_neg_one] using hwithin.tendsto
  unfold logSquareOverOneAddPrimitive26
  convert
    hfirst.add (hproduct.const_mul 2) |>.sub (ht.const_mul 2) using 1
  · funext x
    ring
  · ring

theorem logSquareOverOneSubIntegral26 :
    (∫ x : ℝ in 0..1, Real.log x ^ 2 / (1 - x)) =
      2 * zeta3 := by
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_tendsto
    (f := logSquareOverOneSubPrimitive26)
    (fa := (0 : ℝ)) (fb := 2 * zeta3)
    (by norm_num)
    (fun x hx =>
      logSquareOverOneSubPrimitive26_hasDerivAt hx.1 hx.2)
    logSquareOverOneSub_intervalIntegrable26
    logSquareOverOneSubPrimitive26_tendsto_zero
    logSquareOverOneSubPrimitive26_tendsto_one]
  ring

theorem logSquareOverOneAddIntegral26 :
    (∫ x : ℝ in 0..1, Real.log x ^ 2 / (1 + x)) =
      (3 / 2 : ℝ) * zeta3 := by
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_tendsto
    (f := logSquareOverOneAddPrimitive26)
    (fa := (0 : ℝ)) (fb := (3 / 2 : ℝ) * zeta3)
    (by norm_num)
    (fun x hx =>
      logSquareOverOneAddPrimitive26_hasDerivAt hx.1 hx.2)
    logSquareOverOneAdd_intervalIntegrable26
    logSquareOverOneAddPrimitive26_tendsto_zero
    logSquareOverOneAddPrimitive26_tendsto_one]
  ring

private theorem cubicLogSquareKernel_intervalIntegrable26 :
    IntervalIntegrable
      (fun x : ℝ => x ^ 2 * Real.log x ^ 2 / (1 - x ^ 3))
      MeasureTheory.volume 0 1 := by
  let f : ℝ → ℝ := fun x => x ^ 3
  let f' : ℝ → ℝ := fun x => 3 * x ^ 2
  let g : ℝ → ℝ := fun x => Real.log x ^ 2 / (1 - x)
  have hf : ContinuousOn f [[(0 : ℝ), 1]] := by
    unfold f
    fun_prop
  have hff' :
      ∀ x ∈ Ioo (min (0 : ℝ) 1) (max (0 : ℝ) 1),
        HasDerivAt f (f' x) x := by
    intro x hx
    unfold f f'
    convert (hasDerivAt_id x).pow 3 using 1
    norm_num
  have hf' :
      ∀ x ∈ Ioo (min (0 : ℝ) 1) (max (0 : ℝ) 1),
        0 ≤ f' x := by
    intro x hx
    unfold f'
    positivity
  have hpull :
      IntervalIntegrable (fun x : ℝ => (g ∘ f) x * f' x)
        MeasureTheory.volume 0 1 := by
    apply
      (intervalIntegral.integrable_comp_mul_deriv_iff_of_deriv_nonneg
        (g := g) hf hff' hf').mpr
    simpa [f, g] using logSquareOverOneSub_intervalIntegrable26
  have hscaled := hpull.const_mul (1 / 27 : ℝ)
  apply hscaled.congr
  intro x hx
  dsimp [f, f', g, Function.comp_apply]
  rw [Real.log_pow]
  norm_num only [Nat.cast_ofNat]
  ring

theorem cubicLogSquareIntegral26 :
    (∫ x : ℝ in 0..1,
      x ^ 2 * Real.log x ^ 2 / (1 - x ^ 3)) =
      (2 / 27 : ℝ) * zeta3 := by
  let f : ℝ → ℝ := fun x => x ^ 3
  let f' : ℝ → ℝ := fun x => 3 * x ^ 2
  let g : ℝ → ℝ := fun x => Real.log x ^ 2 / (1 - x)
  have hf : ContinuousOn f [[(0 : ℝ), 1]] := by
    unfold f
    fun_prop
  have hff' :
      ∀ x ∈ Ioo (min (0 : ℝ) 1) (max (0 : ℝ) 1),
        HasDerivAt f (f' x) x := by
    intro x hx
    unfold f f'
    convert (hasDerivAt_id x).pow 3 using 1
    norm_num
  have hf' :
      ∀ x ∈ Ioo (min (0 : ℝ) 1) (max (0 : ℝ) 1),
        0 ≤ f' x := by
    intro x hx
    unfold f'
    positivity
  have hsubst :=
    intervalIntegral.integral_comp_mul_deriv_of_deriv_nonneg
      (a := (0 : ℝ)) (b := 1) (g := g) hf hff' hf'
  have hleft :
      (∫ x : ℝ in 0..1, (g ∘ f) x * f' x) =
        27 * ∫ x : ℝ in 0..1,
          x ^ 2 * Real.log x ^ 2 / (1 - x ^ 3) := by
    rw [← intervalIntegral.integral_const_mul]
    apply intervalIntegral.integral_congr_ae
    exact Filter.Eventually.of_forall fun x hx => by
      have hx' : x ∈ Ioc (0 : ℝ) 1 := by
        simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hx
      have hx0 : 0 < x := hx'.1
      dsimp [f, f', g, Function.comp_apply]
      rw [Real.log_pow]
      norm_num only [Nat.cast_ofNat]
      ring
  rw [hleft] at hsubst
  have hright :
      (∫ x : ℝ in f 0..f 1, g x) = 2 * zeta3 := by
    simpa [f, g] using logSquareOverOneSubIntegral26
  rw [hright] at hsubst
  linarith

private theorem cyclotomicLogSquareKernel_intervalIntegrable26 :
    IntervalIntegrable
      (fun x : ℝ =>
        (2 * x + 1) / (1 + x + x ^ 2) * Real.log x ^ 2)
      MeasureTheory.volume 0 1 := by
  have hdiff :=
    logSquareOverOneSub_intervalIntegrable26.sub
      (cubicLogSquareKernel_intervalIntegrable26.const_mul 3)
  apply hdiff.congr
  intro x hx
  have hx' : x ∈ Ioc (0 : ℝ) 1 := by
    simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hx
  by_cases hx1 : x = 1
  · subst x
    norm_num
  have hone : 1 - x ≠ 0 := sub_ne_zero.mpr (Ne.symm hx1)
  have hq : 1 + x + x ^ 2 ≠ 0 :=
    (cyclotomicQuadratic_pos26 x).ne'
  have hcube : 1 - x ^ 3 ≠ 0 := by
    rw [show 1 - x ^ 3 = (1 - x) * (1 + x + x ^ 2) by ring]
    exact mul_ne_zero hone hq
  field_simp [hone, hq, hcube]
  ring

theorem cyclotomicLogSquareIntegral26 :
    (∫ x : ℝ in 0..1,
      (2 * x + 1) / (1 + x + x ^ 2) * Real.log x ^ 2) =
      (16 / 9 : ℝ) * zeta3 := by
  have hcongr :
      (∫ x : ℝ in 0..1,
        (2 * x + 1) / (1 + x + x ^ 2) * Real.log x ^ 2) =
        ∫ x : ℝ in 0..1,
          Real.log x ^ 2 / (1 - x) -
            3 * (x ^ 2 * Real.log x ^ 2 / (1 - x ^ 3)) := by
    apply intervalIntegral.integral_congr_ae
    exact Filter.Eventually.of_forall fun x hx => by
      have hx' : x ∈ Ioc (0 : ℝ) 1 := by
        simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hx
      by_cases hx1 : x = 1
      · subst x
        norm_num
      have hone : 1 - x ≠ 0 := sub_ne_zero.mpr (Ne.symm hx1)
      have hq : 1 + x + x ^ 2 ≠ 0 :=
        (cyclotomicQuadratic_pos26 x).ne'
      have hcube : 1 - x ^ 3 ≠ 0 := by
        rw [show 1 - x ^ 3 = (1 - x) * (1 + x + x ^ 2) by ring]
        exact mul_ne_zero hone hq
      field_simp [hone, hq, hcube]
      ring
  rw [hcongr,
    intervalIntegral.integral_sub
      logSquareOverOneSub_intervalIntegrable26
      (cubicLogSquareKernel_intervalIntegrable26.const_mul 3),
    intervalIntegral.integral_const_mul,
    logSquareOverOneSubIntegral26,
    cubicLogSquareIntegral26]
  ring

private def minusCyclotomicQuadratic26 (u : ℝ) : ℝ :=
  1 - u + u ^ 2

private def minusCyclotomicDerivative26 (u : ℝ) : ℝ :=
  (2 * u - 1) / minusCyclotomicQuadratic26 u

private theorem minusCyclotomicQuadratic_pos26 (u : ℝ) :
    0 < minusCyclotomicQuadratic26 u := by
  unfold minusCyclotomicQuadratic26
  nlinarith [sq_nonneg (u - 1 / 2)]

private theorem minusCyclotomicDerivative_continuous26 :
    Continuous minusCyclotomicDerivative26 := by
  apply Continuous.div
  · fun_prop
  · unfold minusCyclotomicQuadratic26
    fun_prop
  · intro u
    exact (minusCyclotomicQuadratic_pos26 u).ne'

private def halfCyclotomicKernel26 (u : ℝ) : ℝ :=
  minusCyclotomicDerivative26 u *
    (Real.log (1 - u) ^ 2 -
      Real.log u * Real.log (1 - u))

private def halfCyclotomicReflectedKernel26 (u : ℝ) : ℝ :=
  minusCyclotomicDerivative26 u *
    (Real.log u ^ 2 -
      Real.log u * Real.log (1 - u))

private theorem logOneSub_continuousOn_half26 :
    ContinuousOn (fun u : ℝ => Real.log (1 - u))
      (Icc (0 : ℝ) (1 / 2)) := by
  apply (continuousOn_const.sub continuousOn_id).log
  intro u hu
  exact ne_of_gt (show 0 < 1 - u by linarith [hu.2])

private theorem logSquare_intervalIntegrable_half26 :
    IntervalIntegrable (fun u : ℝ => Real.log u ^ 2)
      MeasureTheory.volume 0 (1 / 2) :=
  logSquare_intervalIntegrable26.mono_set (by
    rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2),
      Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
    intro u hu
    exact ⟨hu.1, hu.2.trans (by norm_num)⟩)

private theorem halfCyclotomicKernel_intervalIntegrable26 :
    IntervalIntegrable halfCyclotomicKernel26
      MeasureTheory.volume 0 (1 / 2) := by
  have hk :
      ContinuousOn minusCyclotomicDerivative26
        (Icc (0 : ℝ) (1 / 2)) :=
    minusCyclotomicDerivative_continuous26.continuousOn
  have hkU :
      ContinuousOn minusCyclotomicDerivative26
        [[(0 : ℝ), (1 / 2 : ℝ)]] := by
    simpa [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2)] using hk
  have hlogU :
      ContinuousOn (fun u : ℝ => Real.log (1 - u))
        [[(0 : ℝ), (1 / 2 : ℝ)]] := by
    simpa [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2)] using
      logOneSub_continuousOn_half26
  have hfirst :
      IntervalIntegrable
        (fun u : ℝ =>
          minusCyclotomicDerivative26 u *
            Real.log (1 - u) ^ 2)
        MeasureTheory.volume 0 (1 / 2) :=
    (hkU.mul (hlogU.pow 2)).intervalIntegrable
  have hcross :
      IntervalIntegrable
        (fun u : ℝ =>
          (minusCyclotomicDerivative26 u *
            Real.log (1 - u)) * Real.log u)
        MeasureTheory.volume 0 (1 / 2) :=
    intervalIntegral.intervalIntegrable_log'.continuousOn_mul
      (hkU.mul hlogU)
  unfold halfCyclotomicKernel26
  convert hfirst.sub hcross using 1
  funext u
  ring

private theorem halfCyclotomicReflectedKernel_intervalIntegrable26 :
    IntervalIntegrable halfCyclotomicReflectedKernel26
      MeasureTheory.volume 0 (1 / 2) := by
  have hk :
      ContinuousOn minusCyclotomicDerivative26
        (Icc (0 : ℝ) (1 / 2)) :=
    minusCyclotomicDerivative_continuous26.continuousOn
  have hkU :
      ContinuousOn minusCyclotomicDerivative26
        [[(0 : ℝ), (1 / 2 : ℝ)]] := by
    simpa [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2)] using hk
  have hlogU :
      ContinuousOn (fun u : ℝ => Real.log (1 - u))
        [[(0 : ℝ), (1 / 2 : ℝ)]] := by
    simpa [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2)] using
      logOneSub_continuousOn_half26
  have hfirst :=
    logSquare_intervalIntegrable_half26.mul_continuousOn hkU
  have hcross :
      IntervalIntegrable
        (fun u : ℝ =>
          (minusCyclotomicDerivative26 u *
            Real.log (1 - u)) * Real.log u)
        MeasureTheory.volume 0 (1 / 2) :=
    intervalIntegral.intervalIntegrable_log'.continuousOn_mul
      (hkU.mul hlogU)
  unfold halfCyclotomicReflectedKernel26
  convert hfirst.sub hcross using 1
  funext u
  ring

private def halfMap26 (x : ℝ) : ℝ :=
  x / (1 + x)

private def halfMapDeriv26 (x : ℝ) : ℝ :=
  1 / (1 + x) ^ 2

private theorem halfMap26_hasDerivAt
    {x : ℝ} (hx : -1 < x) :
    HasDerivAt halfMap26 (halfMapDeriv26 x) x := by
  have hden : 1 + x ≠ 0 := ne_of_gt (by linarith)
  unfold halfMap26 halfMapDeriv26
  convert
    (hasDerivAt_id x).div
      ((hasDerivAt_const x 1).add (hasDerivAt_id x)) hden using 1
  · simp

private theorem halfMap_change_of_variables26 (g : ℝ → ℝ) :
    (∫ x : ℝ in 0..1, (g ∘ halfMap26) x * halfMapDeriv26 x) =
      ∫ u : ℝ in 0..(1 / 2 : ℝ), g u := by
  have hf : ContinuousOn halfMap26 [[(0 : ℝ), 1]] := by
    intro x hx
    have hx' : x ∈ Icc (0 : ℝ) 1 := by
      simpa [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hx
    exact (halfMap26_hasDerivAt (by linarith [hx'.1])).continuousAt.continuousWithinAt
  have hderiv :
      ∀ x ∈ Ioo (min (0 : ℝ) 1) (max (0 : ℝ) 1),
        HasDerivAt halfMap26 (halfMapDeriv26 x) x := by
    intro x hx
    norm_num at hx
    exact halfMap26_hasDerivAt (by linarith [hx.1])
  have hnonneg :
      ∀ x ∈ Ioo (min (0 : ℝ) 1) (max (0 : ℝ) 1),
        0 ≤ halfMapDeriv26 x := by
    intro x hx
    unfold halfMapDeriv26
    positivity
  have hsubst :=
    intervalIntegral.integral_comp_mul_deriv_of_deriv_nonneg
      (a := (0 : ℝ)) (b := 1) (g := g) hf hderiv hnonneg
  convert hsubst using 1 <;> norm_num [halfMap26]

private theorem halfMap_log_difference26
    {x : ℝ} (hx : 0 < x) :
    Real.log (1 - halfMap26 x) - Real.log (halfMap26 x) =
      -Real.log x := by
  have hxne : x ≠ 0 := ne_of_gt hx
  have hden : 1 + x ≠ 0 := ne_of_gt (by linarith)
  have hsub : 1 - halfMap26 x = 1 / (1 + x) := by
    unfold halfMap26
    field_simp [hden]
    ring
  rw [hsub, Real.log_div (by norm_num : (1 : ℝ) ≠ 0) hden,
    show halfMap26 x = x / (1 + x) by rfl,
    Real.log_div hxne hden]
  simp
  ring

private theorem halfMap_cyclotomic_derivative26
    {x : ℝ} (hx : 0 ≤ x) :
    minusCyclotomicDerivative26 (halfMap26 x) *
        halfMapDeriv26 x =
      (2 * x + 1) / (1 + x + x ^ 2) - 2 / (1 + x) := by
  have hden : 1 + x ≠ 0 := ne_of_gt (by linarith)
  have hq : 1 + x + x ^ 2 ≠ 0 :=
    (cyclotomicQuadratic_pos26 x).ne'
  have hminus :
      minusCyclotomicQuadratic26 (halfMap26 x) ≠ 0 :=
    (minusCyclotomicQuadratic_pos26 (halfMap26 x)).ne'
  unfold minusCyclotomicDerivative26 minusCyclotomicQuadratic26
    halfMap26 halfMapDeriv26
  field_simp [hden, hq, hminus]
  ring

private theorem halfMap_sum_integrand26
    {x : ℝ} (hx : 0 < x) :
    ((fun u : ℝ =>
        halfCyclotomicKernel26 u +
          halfCyclotomicReflectedKernel26 u) ∘ halfMap26) x *
        halfMapDeriv26 x =
      (2 * x + 1) / (1 + x + x ^ 2) * Real.log x ^ 2 -
        2 * (Real.log x ^ 2 / (1 + x)) := by
  have hlog := halfMap_log_difference26 hx
  have hrat := halfMap_cyclotomic_derivative26 hx.le
  simp only [Function.comp_apply]
  calc
    (halfCyclotomicKernel26 (halfMap26 x) +
          halfCyclotomicReflectedKernel26 (halfMap26 x)) *
        halfMapDeriv26 x =
      (minusCyclotomicDerivative26 (halfMap26 x) *
          halfMapDeriv26 x) *
        (Real.log (1 - halfMap26 x) -
          Real.log (halfMap26 x)) ^ 2 := by
      unfold halfCyclotomicKernel26 halfCyclotomicReflectedKernel26
      ring
    _ =
      ((2 * x + 1) / (1 + x + x ^ 2) - 2 / (1 + x)) *
        (-Real.log x) ^ 2 := by rw [hrat, hlog]
    _ =
      (2 * x + 1) / (1 + x + x ^ 2) * Real.log x ^ 2 -
        2 * (Real.log x ^ 2 / (1 + x)) := by ring

private theorem halfCyclotomic_sum_integral26 :
    (∫ u : ℝ in 0..(1 / 2 : ℝ),
      halfCyclotomicKernel26 u +
        halfCyclotomicReflectedKernel26 u) =
      -(11 / 9 : ℝ) * zeta3 := by
  have hsubst :=
    halfMap_change_of_variables26
      (fun u : ℝ =>
        halfCyclotomicKernel26 u +
          halfCyclotomicReflectedKernel26 u)
  have hleft :
      (∫ x : ℝ in 0..1,
        (((fun u : ℝ =>
            halfCyclotomicKernel26 u +
              halfCyclotomicReflectedKernel26 u) ∘ halfMap26) x) *
          halfMapDeriv26 x) =
        ∫ x : ℝ in 0..1,
          (2 * x + 1) / (1 + x + x ^ 2) * Real.log x ^ 2 -
            2 * (Real.log x ^ 2 / (1 + x)) := by
    apply intervalIntegral.integral_congr_ae
    exact Filter.Eventually.of_forall fun x hx => by
      have hx' : x ∈ Ioc (0 : ℝ) 1 := by
        simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hx
      exact halfMap_sum_integrand26 hx'.1
  rw [hleft,
    intervalIntegral.integral_sub
      cyclotomicLogSquareKernel_intervalIntegrable26
      (logSquareOverOneAdd_intervalIntegrable26.const_mul 2),
    intervalIntegral.integral_const_mul,
    cyclotomicLogSquareIntegral26,
    logSquareOverOneAddIntegral26] at hsubst
  convert hsubst.symm using 1 <;> ring

private theorem minusCyclotomicDerivative_reflect26 (u : ℝ) :
    minusCyclotomicDerivative26 (1 - u) =
      -minusCyclotomicDerivative26 u := by
  have hq : minusCyclotomicQuadratic26 u ≠ 0 :=
    (minusCyclotomicQuadratic_pos26 u).ne'
  have hqr : minusCyclotomicQuadratic26 (1 - u) ≠ 0 :=
    (minusCyclotomicQuadratic_pos26 (1 - u)).ne'
  unfold minusCyclotomicDerivative26 minusCyclotomicQuadratic26
  field_simp [hq, hqr]
  ring

private def fullCyclotomicReflectionKernel26 (u : ℝ) : ℝ :=
  minusCyclotomicDerivative26 u * Real.log (1 - u) ^ 2

private theorem fullCyclotomicReflectionKernel_intervalIntegrable26 :
    IntervalIntegrable fullCyclotomicReflectionKernel26
      MeasureTheory.volume 0 1 := by
  have hkU :
      ContinuousOn minusCyclotomicDerivative26 [[(0 : ℝ), 1]] :=
    minusCyclotomicDerivative_continuous26.continuousOn
  have hbase :
      IntervalIntegrable
        (fun u : ℝ =>
          minusCyclotomicDerivative26 u * Real.log u ^ 2)
        MeasureTheory.volume 0 1 := by
    convert logSquare_intervalIntegrable26.mul_continuousOn hkU using 1
    funext u
    ring
  have href :=
    (hbase.comp_sub_left 1).symm.neg
  convert href using 1
  · funext u
    unfold fullCyclotomicReflectionKernel26
    simp only [Pi.neg_apply]
    rw [minusCyclotomicDerivative_reflect26]
    ring
  · norm_num
  · norm_num

private theorem halfCyclotomic_diff_eq_full26 :
    (∫ u : ℝ in 0..(1 / 2 : ℝ),
      halfCyclotomicKernel26 u -
        halfCyclotomicReflectedKernel26 u) =
      ∫ u : ℝ in 0..1, fullCyclotomicReflectionKernel26 u := by
  have hE0 :
      IntervalIntegrable fullCyclotomicReflectionKernel26
        MeasureTheory.volume 0 (1 / 2) :=
    fullCyclotomicReflectionKernel_intervalIntegrable26.mono_set (by
      rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2),
        Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
      intro u hu
      exact ⟨hu.1, hu.2.trans (by norm_num)⟩)
  have hE1 :
      IntervalIntegrable fullCyclotomicReflectionKernel26
        MeasureTheory.volume (1 / 2) 1 :=
    fullCyclotomicReflectionKernel_intervalIntegrable26.mono_set (by
      rw [Set.uIcc_of_le (by norm_num : (1 / 2 : ℝ) ≤ 1),
        Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
      intro u hu
      exact ⟨(by linarith [hu.1]), hu.2⟩)
  have hEr :
      IntervalIntegrable
        (fun u : ℝ => fullCyclotomicReflectionKernel26 (1 - u))
        MeasureTheory.volume 0 (1 / 2) := by
    have h := (hE1.comp_sub_left 1).symm
    convert h using 1 <;> norm_num
  have hcongr :
      (∫ u : ℝ in 0..(1 / 2 : ℝ),
        halfCyclotomicKernel26 u -
          halfCyclotomicReflectedKernel26 u) =
        ∫ u : ℝ in 0..(1 / 2 : ℝ),
          fullCyclotomicReflectionKernel26 u +
            fullCyclotomicReflectionKernel26 (1 - u) := by
    apply intervalIntegral.integral_congr
    intro u hu
    unfold halfCyclotomicKernel26 halfCyclotomicReflectedKernel26
      fullCyclotomicReflectionKernel26
    change
      minusCyclotomicDerivative26 u *
            (Real.log (1 - u) ^ 2 -
              Real.log u * Real.log (1 - u)) -
          minusCyclotomicDerivative26 u *
            (Real.log u ^ 2 -
              Real.log u * Real.log (1 - u)) =
        minusCyclotomicDerivative26 u * Real.log (1 - u) ^ 2 +
          minusCyclotomicDerivative26 (1 - u) *
            Real.log (1 - (1 - u)) ^ 2
    rw [show 1 - (1 - u) = u by ring]
    rw [minusCyclotomicDerivative_reflect26]
    ring
  have hreflect :
      (∫ u : ℝ in 0..(1 / 2 : ℝ),
        fullCyclotomicReflectionKernel26 (1 - u)) =
        ∫ u : ℝ in (1 / 2 : ℝ)..1,
          fullCyclotomicReflectionKernel26 u := by
    convert
      intervalIntegral.integral_comp_sub_left
        (a := (0 : ℝ)) (b := (1 / 2 : ℝ))
        fullCyclotomicReflectionKernel26 1 using 1 <;> norm_num
  rw [hcongr,
    intervalIntegral.integral_add hE0 hEr,
    hreflect,
    intervalIntegral.integral_add_adjacent_intervals hE0 hE1]

private theorem cube_change_of_variables26 (g : ℝ → ℝ) :
    (∫ x : ℝ in 0..1,
      (g ∘ fun y : ℝ => y ^ 3) x * (3 * x ^ 2)) =
      ∫ y : ℝ in 0..1, g y := by
  let f : ℝ → ℝ := fun x => x ^ 3
  let f' : ℝ → ℝ := fun x => 3 * x ^ 2
  have hf : ContinuousOn f [[(0 : ℝ), 1]] := by
    unfold f
    fun_prop
  have hderiv :
      ∀ x ∈ Ioo (min (0 : ℝ) 1) (max (0 : ℝ) 1),
        HasDerivAt f (f' x) x := by
    intro x hx
    unfold f f'
    convert (hasDerivAt_id x).pow 3 using 1
    norm_num
  have hnonneg :
      ∀ x ∈ Ioo (min (0 : ℝ) 1) (max (0 : ℝ) 1),
        0 ≤ f' x := by
    intro x hx
    unfold f'
    positivity
  simpa [f, f'] using
    intervalIntegral.integral_comp_mul_deriv_of_deriv_nonneg
      (a := (0 : ℝ)) (b := 1) (g := g) hf hderiv hnonneg

private def cubicRadialKernel26 (x : ℝ) : ℝ :=
  Real.log (1 + x ^ 3) * Real.log x / x

private theorem cubicRadialKernel_intervalIntegrable26 :
    IntervalIntegrable cubicRadialKernel26
      MeasureTheory.volume 0 1 := by
  let f : ℝ → ℝ := fun x => x ^ 3
  let f' : ℝ → ℝ := fun x => 3 * x ^ 2
  have hf : ContinuousOn f [[(0 : ℝ), 1]] := by
    unfold f
    fun_prop
  have hderiv :
      ∀ x ∈ Ioo (min (0 : ℝ) 1) (max (0 : ℝ) 1),
        HasDerivAt f (f' x) x := by
    intro x hx
    unfold f f'
    convert (hasDerivAt_id x).pow 3 using 1
    norm_num
  have hnonneg :
      ∀ x ∈ Ioo (min (0 : ℝ) 1) (max (0 : ℝ) 1),
        0 ≤ f' x := by
    intro x hx
    unfold f'
    positivity
  have hpull :
      IntervalIntegrable
        (fun x : ℝ => (radialWeightThreeKernel26 ∘ f) x * f' x)
        MeasureTheory.volume 0 1 := by
    apply
      (intervalIntegral.integrable_comp_mul_deriv_iff_of_deriv_nonneg
        (g := radialWeightThreeKernel26) hf hderiv hnonneg).mpr
    simpa [f] using radialWeightThreeKernel_intervalIntegrable26
  have hscaled := hpull.const_mul (1 / 9 : ℝ)
  apply hscaled.congr
  intro x hx
  have hx' : x ∈ Ioc (0 : ℝ) 1 := by
    simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hx
  have hxne : x ≠ 0 := ne_of_gt hx'.1
  dsimp [f, f', Function.comp_apply]
  unfold radialWeightThreeKernel26 cubicRadialKernel26
  rw [Real.log_pow]
  norm_num only [Nat.cast_ofNat]
  field_simp [hxne]
  ring

private theorem cubicRadialIntegral26 :
    (∫ x : ℝ in 0..1, cubicRadialKernel26 x) =
      -(1 / 12 : ℝ) * zeta3 := by
  have hsubst := cube_change_of_variables26 radialWeightThreeKernel26
  have hleft :
      (∫ x : ℝ in 0..1,
        (radialWeightThreeKernel26 ∘ fun y : ℝ => y ^ 3) x *
          (3 * x ^ 2)) =
        9 * ∫ x : ℝ in 0..1, cubicRadialKernel26 x := by
    rw [← intervalIntegral.integral_const_mul]
    apply intervalIntegral.integral_congr_ae
    exact Filter.Eventually.of_forall fun x hx => by
      have hx' : x ∈ Ioc (0 : ℝ) 1 := by
        simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hx
      have hxne : x ≠ 0 := ne_of_gt hx'.1
      simp only [Function.comp_apply]
      unfold radialWeightThreeKernel26 cubicRadialKernel26
      rw [Real.log_pow]
      norm_num only [Nat.cast_ofNat]
      field_simp [hxne]
      ring
  rw [hleft, radialWeightThreeIntegral26] at hsubst
  linarith

private def qLogRadialKernel26 (x : ℝ) : ℝ :=
  Real.log (minusCyclotomicQuadratic26 x) * Real.log x / x

private theorem qLogRadialKernel_intervalIntegrable26 :
    IntervalIntegrable qLogRadialKernel26
      MeasureTheory.volume 0 1 := by
  have hdiff :=
    cubicRadialKernel_intervalIntegrable26.sub
      radialWeightThreeKernel_intervalIntegrable26
  apply hdiff.congr
  intro x hx
  have hx' : x ∈ Ioc (0 : ℝ) 1 := by
    simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hx
  have hx0 : 0 < x := hx'.1
  have hxne : x ≠ 0 := ne_of_gt hx0
  have hplus : 1 + x ≠ 0 := ne_of_gt (by linarith)
  have hcubic : 1 + x ^ 3 ≠ 0 := ne_of_gt (by positivity)
  have hfactor :
      minusCyclotomicQuadratic26 x = (1 + x ^ 3) / (1 + x) := by
    unfold minusCyclotomicQuadratic26
    field_simp [hplus]
    ring
  unfold qLogRadialKernel26 cubicRadialKernel26
    radialWeightThreeKernel26
  rw [hfactor, Real.log_div hcubic hplus]
  ring

private theorem qLogRadialIntegral26 :
    (∫ x : ℝ in 0..1, qLogRadialKernel26 x) =
      (2 / 3 : ℝ) * zeta3 := by
  have hcongr :
      (∫ x : ℝ in 0..1, qLogRadialKernel26 x) =
        ∫ x : ℝ in 0..1,
          cubicRadialKernel26 x - radialWeightThreeKernel26 x := by
    apply intervalIntegral.integral_congr
    intro x hx
    have hx' : x ∈ Icc (0 : ℝ) 1 := by
      simpa [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hx
    by_cases hxzero : x = 0
    · subst x
      simp [qLogRadialKernel26, cubicRadialKernel26,
        radialWeightThreeKernel26, minusCyclotomicQuadratic26]
    have hx0 : 0 < x := lt_of_le_of_ne hx'.1 (Ne.symm hxzero)
    have hplus : 1 + x ≠ 0 := ne_of_gt (by linarith)
    have hcubic : 1 + x ^ 3 ≠ 0 := ne_of_gt (by positivity)
    have hfactor :
        minusCyclotomicQuadratic26 x = (1 + x ^ 3) / (1 + x) := by
      unfold minusCyclotomicQuadratic26
      field_simp [hplus]
      ring
    unfold qLogRadialKernel26 cubicRadialKernel26
      radialWeightThreeKernel26
    rw [hfactor, Real.log_div hcubic hplus]
    ring
  rw [hcongr,
    intervalIntegral.integral_sub
      cubicRadialKernel_intervalIntegrable26
      radialWeightThreeKernel_intervalIntegrable26,
    cubicRadialIntegral26, radialWeightThreeIntegral26]
  ring

private theorem logMinusCyclotomicQuadratic_hasDerivAt26 (x : ℝ) :
    HasDerivAt
      (fun u : ℝ => Real.log (minusCyclotomicQuadratic26 u))
      (minusCyclotomicDerivative26 x) x := by
  have hpoly :
      HasDerivAt minusCyclotomicQuadratic26 (2 * x - 1) x := by
    unfold minusCyclotomicQuadratic26
    convert
      ((hasDerivAt_const x 1).sub (hasDerivAt_id x)).add
        ((hasDerivAt_id x).pow 2) using 1 <;> simp <;> ring
  have hne : minusCyclotomicQuadratic26 x ≠ 0 :=
    (minusCyclotomicQuadratic_pos26 x).ne'
  unfold minusCyclotomicDerivative26
  exact hpoly.log hne

private def qLogSquareProduct26 (x : ℝ) : ℝ :=
  Real.log (minusCyclotomicQuadratic26 x) * Real.log x ^ 2

private def qLogSquareDerivativeKernel26 (x : ℝ) : ℝ :=
  minusCyclotomicDerivative26 x * Real.log x ^ 2 +
    2 * qLogRadialKernel26 x

private theorem qLogSquareProduct_hasDerivAt26
    {x : ℝ} (hx : 0 < x) :
    HasDerivAt qLogSquareProduct26
      (qLogSquareDerivativeKernel26 x) x := by
  have hxne : x ≠ 0 := ne_of_gt hx
  have hq := logMinusCyclotomicQuadratic_hasDerivAt26 x
  have hlog : HasDerivAt Real.log (1 / x) x := by
    simpa [one_div] using Real.hasDerivAt_log hxne
  unfold qLogSquareProduct26 qLogSquareDerivativeKernel26
    qLogRadialKernel26
  convert hq.mul (hlog.pow 2) using 1 <;>
    simp only [Pi.pow_apply] <;>
    field_simp [hxne] <;> ring

private theorem qLogSquareDerivativeKernel_intervalIntegrable26 :
    IntervalIntegrable qLogSquareDerivativeKernel26
      MeasureTheory.volume 0 1 := by
  have hkU :
      ContinuousOn minusCyclotomicDerivative26 [[(0 : ℝ), 1]] :=
    minusCyclotomicDerivative_continuous26.continuousOn
  have hfirst :
      IntervalIntegrable
        (fun x : ℝ =>
          minusCyclotomicDerivative26 x * Real.log x ^ 2)
        MeasureTheory.volume 0 1 := by
    convert logSquare_intervalIntegrable26.mul_continuousOn hkU using 1
    funext x
    ring
  unfold qLogSquareDerivativeKernel26
  exact hfirst.add (qLogRadialKernel_intervalIntegrable26.const_mul 2)

private theorem qLogSquareProduct_tendsto_zero26 :
    Tendsto qLogSquareProduct26
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  have hderiv :
      HasDerivAt
        (fun x : ℝ => Real.log (minusCyclotomicQuadratic26 x))
        (-1) 0 := by
    simpa [minusCyclotomicDerivative26,
      minusCyclotomicQuadratic26] using
        logMinusCyclotomicQuadratic_hasDerivAt26 0
  have hslope :
      Tendsto
        (fun x : ℝ =>
          x⁻¹ * Real.log (minusCyclotomicQuadratic26 x))
        (𝓝[>] (0 : ℝ)) (𝓝 (-1)) := by
    simpa [minusCyclotomicQuadratic26] using
      hderiv.tendsto_slope_zero_right
  have hraw := logSquare_mul_self_tendsto_zero26.mul hslope
  have hraw0 :
      Tendsto
        (fun x : ℝ =>
          (Real.log x ^ 2 * x) *
            (x⁻¹ * Real.log (minusCyclotomicQuadratic26 x)))
        (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    simpa using hraw
  unfold qLogSquareProduct26
  apply Filter.Tendsto.congr' _ hraw0
  filter_upwards [self_mem_nhdsWithin] with x hx
  have hxne : x ≠ 0 := ne_of_gt hx
  field_simp [hxne]

private theorem qLogSquareProduct_tendsto_one26 :
    Tendsto qLogSquareProduct26
      (𝓝[<] (1 : ℝ)) (𝓝 0) := by
  have hq :
      Tendsto
        (fun x : ℝ => Real.log (minusCyclotomicQuadratic26 x))
        (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    have hcont :=
      (logMinusCyclotomicQuadratic_hasDerivAt26 1).continuousAt
    have ht :=
      hcont.tendsto.mono_left
        (show (𝓝[<] (1 : ℝ)) ≤ 𝓝 1 from inf_le_left)
    simpa [minusCyclotomicQuadratic26] using ht
  have hlog :
      Tendsto Real.log (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    simpa using tendsto_nhdsWithin_of_tendsto_nhds
      (Real.continuousAt_log (by norm_num : (1 : ℝ) ≠ 0)).tendsto
  unfold qLogSquareProduct26
  simpa using hq.mul (hlog.pow 2)

private theorem qLogSquareDerivativeIntegral26 :
    (∫ x : ℝ in 0..1, qLogSquareDerivativeKernel26 x) = 0 := by
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_tendsto
    (f := qLogSquareProduct26)
    (fa := (0 : ℝ)) (fb := 0)
    (by norm_num)
    (fun x hx => qLogSquareProduct_hasDerivAt26 hx.1)
    qLogSquareDerivativeKernel_intervalIntegrable26
    qLogSquareProduct_tendsto_zero26
    qLogSquareProduct_tendsto_one26]
  ring

private theorem minusCyclotomicLogSquareIntegral26 :
    (∫ x : ℝ in 0..1,
      minusCyclotomicDerivative26 x * Real.log x ^ 2) =
      -(4 / 3 : ℝ) * zeta3 := by
  have hkU :
      ContinuousOn minusCyclotomicDerivative26 [[(0 : ℝ), 1]] :=
    minusCyclotomicDerivative_continuous26.continuousOn
  have hfirst :
      IntervalIntegrable
        (fun x : ℝ =>
          minusCyclotomicDerivative26 x * Real.log x ^ 2)
        MeasureTheory.volume 0 1 := by
    convert logSquare_intervalIntegrable26.mul_continuousOn hkU using 1
    funext x
    ring
  have h := qLogSquareDerivativeIntegral26
  unfold qLogSquareDerivativeKernel26 at h
  rw [intervalIntegral.integral_add
      hfirst (qLogRadialKernel_intervalIntegrable26.const_mul 2),
    intervalIntegral.integral_const_mul,
    qLogRadialIntegral26] at h
  linarith

private theorem fullCyclotomicReflectionIntegral26 :
    (∫ x : ℝ in 0..1, fullCyclotomicReflectionKernel26 x) =
      (4 / 3 : ℝ) * zeta3 := by
  let f : ℝ → ℝ :=
    fun x =>
      minusCyclotomicDerivative26 x * Real.log x ^ 2
  have hreflect :
      (∫ x : ℝ in 0..1, f (1 - x)) =
        ∫ x : ℝ in 0..1, f x := by
    simpa using
      (intervalIntegral.integral_comp_sub_left
        (a := (0 : ℝ)) (b := 1) f 1)
  have hcongr :
      (∫ x : ℝ in 0..1, fullCyclotomicReflectionKernel26 x) =
        -(∫ x : ℝ in 0..1, f (1 - x)) := by
    rw [← intervalIntegral.integral_neg]
    apply intervalIntegral.integral_congr
    intro x hx
    dsimp [f]
    unfold fullCyclotomicReflectionKernel26
    rw [minusCyclotomicDerivative_reflect26]
    ring
  rw [hcongr, hreflect]
  simpa [f] using congrArg Neg.neg minusCyclotomicLogSquareIntegral26

private theorem halfCyclotomicIntegral26 :
    (∫ u : ℝ in 0..(1 / 2 : ℝ), halfCyclotomicKernel26 u) =
      (1 / 18 : ℝ) * zeta3 := by
  have hsum := halfCyclotomic_sum_integral26
  rw [intervalIntegral.integral_add
    halfCyclotomicKernel_intervalIntegrable26
    halfCyclotomicReflectedKernel_intervalIntegrable26] at hsum
  have hdiff := halfCyclotomic_diff_eq_full26
  rw [intervalIntegral.integral_sub
      halfCyclotomicKernel_intervalIntegrable26
      halfCyclotomicReflectedKernel_intervalIntegrable26,
    fullCyclotomicReflectionIntegral26] at hdiff
  linarith

private theorem halfMap_log_one_sub26
    {x : ℝ} (hx : 0 ≤ x) :
    Real.log (1 - halfMap26 x) = -Real.log (1 + x) := by
  have hden : 1 + x ≠ 0 := ne_of_gt (by linarith)
  have hsub : 1 - halfMap26 x = 1 / (1 + x) := by
    unfold halfMap26
    field_simp [hden]
    ring
  rw [hsub, Real.log_div (by norm_num : (1 : ℝ) ≠ 0) hden]
  simp

private theorem halfMap_main_integrand26
    {x : ℝ} (hx : 0 < x) :
    (halfCyclotomicKernel26 ∘ halfMap26) x *
        halfMapDeriv26 x =
      cyclotomicWeightThreeKernel26 x -
        2 * alternatingWeightThreeKernel26 x := by
  have hdiff := halfMap_log_difference26 hx
  have hone := halfMap_log_one_sub26 hx.le
  have hrat := halfMap_cyclotomic_derivative26 hx.le
  simp only [Function.comp_apply]
  calc
    halfCyclotomicKernel26 (halfMap26 x) *
        halfMapDeriv26 x =
      (minusCyclotomicDerivative26 (halfMap26 x) *
          halfMapDeriv26 x) *
        (Real.log (1 - halfMap26 x) *
          (Real.log (1 - halfMap26 x) -
            Real.log (halfMap26 x))) := by
      unfold halfCyclotomicKernel26
      ring
    _ =
      ((2 * x + 1) / (1 + x + x ^ 2) - 2 / (1 + x)) *
        ((-Real.log (1 + x)) * (-Real.log x)) := by
      rw [hrat, hdiff, hone]
    _ =
      cyclotomicWeightThreeKernel26 x -
        2 * alternatingWeightThreeKernel26 x := by
      unfold cyclotomicWeightThreeKernel26
        alternatingWeightThreeKernel26
      ring

/-- The cubic cyclotomic logarithmic special value required by the
Problem 2.6 assembly. -/
theorem cyclotomicLogIntegralEvaluation26 :
    CyclotomicLogIntegralEvaluation26 := by
  unfold CyclotomicLogIntegralEvaluation26
  have hsubst := halfMap_change_of_variables26 halfCyclotomicKernel26
  have hleft :
      (∫ x : ℝ in 0..1,
        (halfCyclotomicKernel26 ∘ halfMap26) x *
          halfMapDeriv26 x) =
        ∫ x : ℝ in 0..1,
          cyclotomicWeightThreeKernel26 x -
            2 * alternatingWeightThreeKernel26 x := by
    apply intervalIntegral.integral_congr_ae
    exact Filter.Eventually.of_forall fun x hx => by
      have hx' : x ∈ Ioc (0 : ℝ) 1 := by
        simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hx
      exact halfMap_main_integrand26 hx'.1
  rw [hleft,
    intervalIntegral.integral_sub
      cyclotomicWeightThreeKernel_intervalIntegrable26
      (alternatingWeightThreeKernel_intervalIntegrable26.const_mul 2),
    intervalIntegral.integral_const_mul,
    halfCyclotomicIntegral26,
    alternatingWeightThreeIntegral26] at hsubst
  linarith

end RamanujanChallenge.P26
