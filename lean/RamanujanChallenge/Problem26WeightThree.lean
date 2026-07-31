import RamanujanChallenge.Problem26Nested

/-!
  Real polylogarithm evaluations for the weight-three part of Problem 2.6.

  This file proves the alternating logarithmic integral unconditionally.  It
  is downstream from `Problem26Nested`, whose remaining cyclotomic integral is
  treated separately.
-/

open Filter Set Topology
open scoped Interval

noncomputable section

namespace RamanujanChallenge.P26

private def dilogLandenAux26 (x : ℝ) : ℝ :=
  dilog (-x / (1 - x)) + dilog x +
    (1 / 2 : ℝ) * Real.log (1 - x) ^ 2

private theorem dilogLandenAux_hasDerivAt_zero26
    {x : ℝ} (hx0 : 0 < x) (hxhalf : x < 1 / 2) :
    HasDerivAt dilogLandenAux26 0 x := by
  have hxne : x ≠ 0 := ne_of_gt hx0
  have h1x0 : 0 < 1 - x := by linarith
  have h1xne : 1 - x ≠ 0 := ne_of_gt h1x0
  let y : ℝ := -x / (1 - x)
  have hyneg : y < 0 := div_neg_of_neg_of_pos (neg_lt_zero.mpr hx0) h1x0
  have hyne : y ≠ 0 := ne_of_lt hyneg
  have hygt : -1 < y := by
    dsimp [y]
    apply (lt_div_iff₀ h1x0).2
    linarith
  have hyabs : |y| < 1 := (abs_lt).2 ⟨hygt, hyneg.trans zero_lt_one⟩
  have hyderiv :
      HasDerivAt (fun z : ℝ => -z / (1 - z))
        (-1 / (1 - x) ^ 2) x := by
    convert
      (hasDerivAt_id x).neg.div
        ((hasDerivAt_const x 1).sub (hasDerivAt_id x)) h1xne using 1
    simp only [Pi.sub_apply, Pi.neg_apply, id_eq]
    field_simp [h1xne]
    ring
  have hdY :
      HasDerivAt (fun z : ℝ => dilog (-z / (1 - z)))
        ((-(Real.log (1 - y)) / y) * (-1 / (1 - x) ^ 2)) x :=
    by
      simpa [y] using
        (dilog_hasDerivAt_of_abs_lt_one hyabs hyne).comp x hyderiv
  have hdX := dilog_hasDerivAt hx0 (by linarith : x < 1)
  have hlog :
      HasDerivAt (fun z : ℝ => Real.log (1 - z))
        (-1 / (1 - x)) x :=
    by
      convert
        ((hasDerivAt_const x 1).sub (hasDerivAt_id x)).log h1xne using 1
      simp
  unfold dilogLandenAux26
  convert (hdY.add hdX).add ((hlog.pow 2).const_mul (1 / 2 : ℝ)) using 1
  have hone :
      1 - y = (1 - x)⁻¹ := by
    dsimp [y]
    field_simp [h1xne]
    ring
  rw [hone, Real.log_inv]
  dsimp [y]
  field_simp [hxne, h1xne]
  ring

private theorem dilogLandenAux_continuousOn26 :
    ContinuousOn dilogLandenAux26 (Icc (0 : ℝ) (1 / 2)) := by
  have hy : ContinuousOn (fun x : ℝ => -x / (1 - x))
      (Icc (0 : ℝ) (1 / 2)) := by
    apply (continuousOn_id.neg).div (continuousOn_const.sub continuousOn_id)
    intro x hx
    simpa only [Pi.sub_apply, Pi.one_apply, id_eq] using
      (ne_of_gt (show 0 < 1 - x by linarith [hx.2]))
  have hymem :
      MapsTo (fun x : ℝ => -x / (1 - x))
        (Icc (0 : ℝ) (1 / 2)) (Icc (-1 : ℝ) 1) := by
    intro x hx
    have h1x0 : 0 < 1 - x := by linarith [hx.2]
    constructor
    · apply (le_div_iff₀ h1x0).2
      linarith [hx.2]
    · have hnonpos : -x / (1 - x) ≤ 0 :=
        div_nonpos_of_nonpos_of_nonneg (neg_nonpos.mpr hx.1) h1x0.le
      linarith
  have hdY :
      ContinuousOn (fun x : ℝ => dilog (-x / (1 - x)))
        (Icc (0 : ℝ) (1 / 2)) :=
    dilog_continuousOn_unit.comp hy hymem
  have hdX :
      ContinuousOn dilog (Icc (0 : ℝ) (1 / 2)) :=
    dilog_continuousOn_unit.mono (by
      intro x hx
      constructor <;> linarith [hx.1, hx.2])
  have hlog :
      ContinuousOn (fun x : ℝ => Real.log (1 - x))
        (Icc (0 : ℝ) (1 / 2)) := by
    apply (continuousOn_const.sub continuousOn_id).log
    intro x hx
    simpa only [Pi.sub_apply, Pi.one_apply, id_eq] using
      (ne_of_gt (show 0 < 1 - x by linarith [hx.2]))
  unfold dilogLandenAux26
  exact (hdY.add hdX).add (continuousOn_const.mul (hlog.pow 2))

theorem dilog_landen_half26
    {x : ℝ} (hx0 : 0 < x) (hxhalf : x ≤ 1 / 2) :
    dilog (-x / (1 - x)) =
      -dilog x - (1 / 2 : ℝ) * Real.log (1 - x) ^ 2 := by
  have hconst := constant_of_has_deriv_right_zero
    dilogLandenAux_continuousOn26
    (fun y hy => by
      rcases eq_or_lt_of_le hy.1 with rfl | hy0
      · have hyderiv :
            HasDerivAt (fun z : ℝ => -z / (1 - z)) (-1) 0 := by
          convert
            (hasDerivAt_id (0 : ℝ)).neg.div
              ((hasDerivAt_const (0 : ℝ) 1).sub (hasDerivAt_id 0))
              (by norm_num) using 1
          norm_num
        have hdY :
            HasDerivAt (fun z : ℝ => dilog (-z / (1 - z))) (-1) 0 := by
          have houter :
              HasDerivAt dilog 1 (-((0 : ℝ)) / (1 - 0)) := by
            simpa using dilog_hasDerivAt_zero26
          convert houter.comp 0 hyderiv using 1 <;> norm_num
        have hdX : HasDerivAt dilog 1 0 := dilog_hasDerivAt_zero26
        have hlog :
            HasDerivAt (fun z : ℝ => Real.log (1 - z)) (-1) 0 := by
          convert
            ((hasDerivAt_const (0 : ℝ) 1).sub (hasDerivAt_id 0)).log
              (by norm_num) using 1
          norm_num
        have htotal : HasDerivAt dilogLandenAux26 0 0 := by
          unfold dilogLandenAux26
          convert
            (hdY.add hdX).add ((hlog.pow 2).const_mul (1 / 2 : ℝ)) using 1
          norm_num
        exact htotal.hasDerivWithinAt
      · exact
          (dilogLandenAux_hasDerivAt_zero26 hy0 hy.2).hasDerivWithinAt)
  have hxmem : x ∈ Icc (0 : ℝ) (1 / 2) := ⟨hx0.le, hxhalf⟩
  have heq := hconst x hxmem
  simp [dilogLandenAux26, dilog_zero] at heq
  linarith

private def logOneMinusSlope26 (x : ℝ) : ℝ :=
  Function.update (fun y : ℝ => Real.log (1 - y) / y) 0 (-1) x

private theorem logOneMinusSlope26_continuousOn :
    ContinuousOn logOneMinusSlope26 (Icc (0 : ℝ) (1 / 2)) := by
  intro x hx
  by_cases hxzero : x = 0
  · subst x
    have hlog :
        HasDerivAt (fun y : ℝ => Real.log (1 - y)) (-1) 0 := by
      have hone :
          HasDerivAt (fun y : ℝ => 1 - y) (-1) 0 := by
        convert (hasDerivAt_const (0 : ℝ) 1).sub (hasDerivAt_id 0) using 1
        norm_num
      simpa [Function.comp_def] using
        (HasDerivAt.comp (h := fun y : ℝ => 1 - y) 0
          (Real.hasDerivAt_log
            (by norm_num : (fun y : ℝ => 1 - y) 0 ≠ 0)) hone)
    have hc := hlog.continuousAt_div
    have hc' : ContinuousAt logOneMinusSlope26 0 := by
      convert hc using 1
      funext y
      simp [logOneMinusSlope26]
    exact hc'.continuousWithinAt
  · have hbase :
        ContinuousAt (fun y : ℝ => Real.log (1 - y) / y) x := by
      have hone : ContinuousAt (fun y : ℝ => 1 - y) x := by fun_prop
      have hlog :
          ContinuousAt (fun y : ℝ => Real.log (1 - y)) x :=
        (Real.continuousAt_log
          (by linarith [hx.2] : 1 - x ≠ 0)).comp hone
      exact hlog.div continuousAt_id hxzero
    have heq :
        logOneMinusSlope26 =ᶠ[𝓝 x]
          (fun y : ℝ => Real.log (1 - y) / y) := by
      filter_upwards [eventually_ne_nhds hxzero] with y hy
      simp [logOneMinusSlope26, hy]
    exact hbase.congr_of_eventuallyEq heq |>.continuousWithinAt

private theorem log_mul_sq_log_one_sub_continuousOn26 :
    ContinuousOn
      (fun x : ℝ => Real.log x * Real.log (1 - x) ^ 2)
      (Icc (0 : ℝ) (1 / 2)) := by
  have hlog :
      ContinuousOn (fun x : ℝ => Real.log (1 - x))
        (Icc (0 : ℝ) (1 / 2)) := by
    apply (continuousOn_const.sub continuousOn_id).log
    intro x hx
    simpa only [Pi.sub_apply, Pi.one_apply, id_eq] using
      (ne_of_gt (show 0 < 1 - x by linarith [hx.2]))
  have hright :
      ContinuousOn
        (fun x : ℝ =>
          (x * Real.log x) * logOneMinusSlope26 x * Real.log (1 - x))
        (Icc (0 : ℝ) (1 / 2)) :=
    (Real.continuous_mul_log.continuousOn.mul
      logOneMinusSlope26_continuousOn).mul hlog
  apply hright.congr
  intro x hx
  by_cases hxzero : x = 0
  · subst x
    simp
  · simp [logOneMinusSlope26, hxzero]
    field_simp [hxzero]

private def trilogLandenAux26 (x : ℝ) : ℝ :=
  trilog26 x + trilog26 (-x / (1 - x)) + trilog26 (1 - x) -
      zeta3 - Real.pi ^ 2 / 6 * Real.log (1 - x) +
    (1 / 2 : ℝ) * Real.log x * Real.log (1 - x) ^ 2 -
    (1 / 6 : ℝ) * Real.log (1 - x) ^ 3

private theorem trilogLandenAux_continuousOn26 :
    ContinuousOn trilogLandenAux26 (Icc (0 : ℝ) (1 / 2)) := by
  have hxcont : ContinuousOn trilog26 (Icc (0 : ℝ) (1 / 2)) :=
    trilog26_continuousOn_unit.mono (by
      intro x hx
      constructor <;> linarith [hx.1, hx.2])
  have hycont : ContinuousOn (fun x : ℝ => -x / (1 - x))
      (Icc (0 : ℝ) (1 / 2)) := by
    apply (continuousOn_id.neg).div (continuousOn_const.sub continuousOn_id)
    intro x hx
    simpa only [Pi.sub_apply, Pi.one_apply, id_eq] using
      (ne_of_gt (show 0 < 1 - x by linarith [hx.2]))
  have hymem :
      MapsTo (fun x : ℝ => -x / (1 - x))
        (Icc (0 : ℝ) (1 / 2)) (Icc (-1 : ℝ) 1) := by
    intro x hx
    have h1x0 : 0 < 1 - x := by linarith [hx.2]
    constructor
    · apply (le_div_iff₀ h1x0).2
      linarith [hx.2]
    · have hnonpos : -x / (1 - x) ≤ 0 :=
        div_nonpos_of_nonpos_of_nonneg (neg_nonpos.mpr hx.1) h1x0.le
      linarith
  have hytri :
      ContinuousOn (fun x : ℝ => trilog26 (-x / (1 - x)))
        (Icc (0 : ℝ) (1 / 2)) :=
    trilog26_continuousOn_unit.comp hycont hymem
  have hsubcont : ContinuousOn (fun x : ℝ => 1 - x)
      (Icc (0 : ℝ) (1 / 2)) :=
    continuousOn_const.sub continuousOn_id
  have hsubmem :
      MapsTo (fun x : ℝ => 1 - x)
        (Icc (0 : ℝ) (1 / 2)) (Icc (-1 : ℝ) 1) := by
    intro x hx
    constructor <;> linarith [hx.1, hx.2]
  have hsubtri :
      ContinuousOn (fun x : ℝ => trilog26 (1 - x))
        (Icc (0 : ℝ) (1 / 2)) :=
    trilog26_continuousOn_unit.comp hsubcont hsubmem
  have hlog :
      ContinuousOn (fun x : ℝ => Real.log (1 - x))
        (Icc (0 : ℝ) (1 / 2)) := by
    apply hsubcont.log
    intro x hx
    simpa only [Pi.sub_apply, Pi.one_apply, id_eq] using
      (ne_of_gt (show 0 < 1 - x by linarith [hx.2]))
  have hpi :
      ContinuousOn
        (fun x : ℝ => Real.pi ^ 2 / 6 * Real.log (1 - x))
        (Icc (0 : ℝ) (1 / 2)) :=
    continuousOn_const.mul hlog
  have hhalfprod :
      ContinuousOn
        (fun x : ℝ =>
          (1 / 2 : ℝ) * (Real.log x * Real.log (1 - x) ^ 2))
        (Icc (0 : ℝ) (1 / 2)) :=
    continuousOn_const.mul log_mul_sq_log_one_sub_continuousOn26
  have hcub :
      ContinuousOn
        (fun x : ℝ => (1 / 6 : ℝ) * Real.log (1 - x) ^ 3)
        (Icc (0 : ℝ) (1 / 2)) :=
    continuousOn_const.mul (hlog.pow 3)
  have hzeta :
      ContinuousOn (fun _ : ℝ => zeta3) (Icc (0 : ℝ) (1 / 2)) :=
    continuousOn_const
  unfold trilogLandenAux26
  have htotal :=
    (((((hxcont.add hytri).add hsubtri).sub hzeta).sub hpi).add
      hhalfprod).sub hcub
  simpa only [Pi.add_apply, Pi.sub_apply, mul_assoc] using htotal

private theorem trilogLandenAux_hasDerivAt_zero26
    {x : ℝ} (hx0 : 0 < x) (hxhalf : x < 1 / 2) :
    HasDerivAt trilogLandenAux26 0 x := by
  have hxne : x ≠ 0 := ne_of_gt hx0
  have hx1 : x < 1 := by linarith
  have h1x0 : 0 < 1 - x := by linarith
  have h1x1 : 1 - x < 1 := by linarith
  have h1xne : 1 - x ≠ 0 := ne_of_gt h1x0
  let y : ℝ := -x / (1 - x)
  have hyneg : y < 0 :=
    div_neg_of_neg_of_pos (neg_lt_zero.mpr hx0) h1x0
  have hyne : y ≠ 0 := ne_of_lt hyneg
  have hygt : -1 < y := by
    dsimp [y]
    apply (lt_div_iff₀ h1x0).2
    linarith
  have hyabs : |y| < 1 :=
    (abs_lt).2 ⟨hygt, hyneg.trans zero_lt_one⟩
  have hyderiv :
      HasDerivAt (fun z : ℝ => -z / (1 - z))
        (-1 / (1 - x) ^ 2) x := by
    convert
      (hasDerivAt_id x).neg.div
        ((hasDerivAt_const x 1).sub (hasDerivAt_id x)) h1xne using 1
    simp only [Pi.sub_apply, Pi.neg_apply, id_eq]
    field_simp [h1xne]
    ring
  have hTx := trilog26_hasDerivAt_of_abs_lt_one
    (by rw [abs_of_pos hx0]; exact hx1) hxne
  have hTy :
      HasDerivAt (fun z : ℝ => trilog26 (-z / (1 - z)))
        ((dilog y / y) * (-1 / (1 - x) ^ 2)) x := by
    simpa [y] using
      (trilog26_hasDerivAt_of_abs_lt_one hyabs hyne).comp x hyderiv
  have hsub :
      HasDerivAt (fun z : ℝ => 1 - z) (-1) x := by
    convert (hasDerivAt_const x 1).sub (hasDerivAt_id x) using 1
    simp
  have hTsub :
      HasDerivAt (fun z : ℝ => trilog26 (1 - z))
        (-(dilog (1 - x) / (1 - x))) x := by
    convert
      (trilog26_hasDerivAt_of_abs_lt_one
        (by rw [abs_of_pos h1x0]; exact h1x1) h1xne).comp x hsub using 1
    ring
  have hlogx : HasDerivAt (fun z : ℝ => Real.log z) (1 / x) x := by
    simpa [one_div] using Real.hasDerivAt_log hxne
  have hlog1 :
      HasDerivAt (fun z : ℝ => Real.log (1 - z))
        (-1 / (1 - x)) x := by
    convert hsub.log h1xne using 1
  have htotal :=
    (((((hTx.add hTy).add hTsub).sub_const zeta3).sub
      ((hlog1.const_mul (Real.pi ^ 2 / 6)))).add
      ((hlogx.mul (hlog1.pow 2)).const_mul (1 / 2 : ℝ))).sub
      ((hlog1.pow 3).const_mul (1 / 6 : ℝ))
  unfold trilogLandenAux26
  convert htotal using 1
  · funext z
    simp only [Pi.add_apply, Pi.sub_apply, Pi.mul_apply, Pi.pow_apply]
    ring
  · have hLanden := dilog_landen_half26 hx0 hxhalf.le
    have hReflect := dilog_reflection hx0 hx1
    dsimp [y] at hLanden ⊢
    field_simp [hxne, h1xne] at hLanden hReflect ⊢
    nlinarith [hLanden, hReflect]

theorem trilog26_half :
    trilog26 (1 / 2) =
      (7 : ℝ) / 8 * zeta3 -
        Real.pi ^ 2 * Real.log 2 / 12 +
        Real.log 2 ^ 3 / 6 := by
  let a : ℕ → ℝ := fun n => 1 / ((n : ℝ) + 2)
  have ha0 : ∀ n, 0 < a n := by
    intro n
    dsimp [a]
    positivity
  have hahalf : ∀ n, a n ≤ 1 / 2 := by
    intro n
    dsimp [a]
    apply (div_le_iff₀ (by positivity : (0 : ℝ) < (n : ℝ) + 2)).2
    have hn : (0 : ℝ) ≤ n := by positivity
    nlinarith
  have heq : ∀ n, trilogLandenAux26 (1 / 2) = trilogLandenAux26 (a n) := by
    intro n
    have hsubset :
        Icc (a n) (1 / 2) ⊆ Icc (0 : ℝ) (1 / 2) := by
      intro x hx
      exact ⟨le_trans (ha0 n).le hx.1, hx.2⟩
    have hcont :
        ContinuousOn trilogLandenAux26 (Icc (a n) (1 / 2)) :=
      trilogLandenAux_continuousOn26.mono hsubset
    have hconst := constant_of_has_deriv_right_zero hcont
      (fun x hx =>
        (trilogLandenAux_hasDerivAt_zero26
          (lt_of_lt_of_le (ha0 n) hx.1) hx.2).hasDerivWithinAt)
    have hhalfmem : (1 / 2 : ℝ) ∈ Icc (a n) (1 / 2) :=
      ⟨hahalf n, le_rfl⟩
    exact hconst (1 / 2) hhalfmem
  have habase :
      Tendsto (fun n : ℕ => 1 / ((n : ℝ) + 1))
        atTop (𝓝 0) :=
    tendsto_one_div_add_atTop_nhds_zero_nat
  have hatendsto : Tendsto a atTop (𝓝 0) := by
    have hshift := (tendsto_add_atTop_iff_nat 1).2 habase
    convert hshift using 1
    funext n
    dsimp [a]
    push_cast
    ring
  have haWithin : Tendsto a atTop (𝓝[Icc (0 : ℝ) (1 / 2)] 0) := by
    rw [tendsto_nhdsWithin_iff]
    refine ⟨hatendsto, ?_⟩
    exact Eventually.of_forall fun n => ⟨(ha0 n).le, hahalf n⟩
  have haux :
      Tendsto (fun n => trilogLandenAux26 (a n))
        atTop (𝓝 (trilogLandenAux26 0)) :=
    (trilogLandenAux_continuousOn26 0
      (by norm_num : (0 : ℝ) ∈ Icc 0 (1 / 2))).tendsto.comp haWithin
  have hconst :
      Tendsto (fun _ : ℕ => trilogLandenAux26 (1 / 2))
        atTop (𝓝 (trilogLandenAux26 (1 / 2))) :=
    tendsto_const_nhds
  have heqEventually :
      (fun _ : ℕ => trilogLandenAux26 (1 / 2)) =ᶠ[atTop]
        (fun n => trilogLandenAux26 (a n)) :=
    Eventually.of_forall heq
  have hzero :
      trilogLandenAux26 (1 / 2) = trilogLandenAux26 0 :=
    tendsto_nhds_unique hconst (haux.congr' heqEventually.symm)
  have hauxzero : trilogLandenAux26 0 = 0 := by
    have htri1 : trilog26 1 = zeta3 := by
      unfold trilog26 zeta3
      apply tsum_congr
      intro n
      simp [Nat.cast_add]
    simp [trilogLandenAux26, htri1]
  rw [hauxzero] at hzero
  have hloghalf : Real.log (1 / 2 : ℝ) = -Real.log 2 := by
    rw [one_div, Real.log_inv]
  simp only [trilogLandenAux26] at hzero
  rw [show -(1 / 2 : ℝ) / (1 - 1 / 2) = -1 by norm_num,
    show 1 - (1 / 2 : ℝ) = 1 / 2 by norm_num,
    trilog26_neg_one, hloghalf] at hzero
  linarith

theorem dilog26_half :
    dilog (1 / 2) =
      Real.pi ^ 2 / 12 - Real.log 2 ^ 2 / 2 := by
  have h := dilog_reflection
    (z := (1 / 2 : ℝ)) (by norm_num) (by norm_num)
  have hloghalf : Real.log (1 / 2 : ℝ) = -Real.log 2 := by
    rw [one_div, Real.log_inv]
  rw [show 1 - (1 / 2 : ℝ) = 1 / 2 by norm_num, hloghalf] at h
  linarith

private def logSquareTransform26 (u : ℝ) : ℝ :=
  (1 / 3 : ℝ) * Real.log u ^ 3 -
      Real.log u ^ 2 * Real.log (1 - u) -
    2 * Real.log u * dilog u + 2 * trilog26 u

private theorem logSquareTransform_hasDerivAt26
    {u : ℝ} (hu0 : 0 < u) (hu1 : u < 1) :
    HasDerivAt logSquareTransform26
      (Real.log u ^ 2 / (u * (1 - u))) u := by
  have hune : u ≠ 0 := ne_of_gt hu0
  have h1u0 : 0 < 1 - u := sub_pos.mpr hu1
  have h1une : 1 - u ≠ 0 := ne_of_gt h1u0
  have hlu : HasDerivAt (fun x : ℝ => Real.log x) (1 / u) u := by
    simpa [one_div] using Real.hasDerivAt_log hune
  have hsub : HasDerivAt (fun x : ℝ => 1 - x) (-1) u := by
    convert (hasDerivAt_const u 1).sub (hasDerivAt_id u) using 1
    simp
  have hl1 :
      HasDerivAt (fun x : ℝ => Real.log (1 - x))
        (-1 / (1 - u)) u := by
    convert hsub.log h1une using 1
  have hd := dilog_hasDerivAt hu0 hu1
  have ht := trilog26_hasDerivAt_of_abs_lt_one
    (by rw [abs_of_pos hu0]; exact hu1) hune
  unfold logSquareTransform26
  have htotal :=
    ((((hlu.pow 3).const_mul (1 / 3 : ℝ)).sub
      ((hlu.pow 2).mul hl1)).sub
      ((hlu.mul hd).const_mul 2)).add (ht.const_mul 2)
  convert htotal using 1
  · funext x
    simp only [Pi.add_apply, Pi.sub_apply, Pi.mul_apply, Pi.pow_apply]
    ring
  · simp only [Pi.pow_apply]
    field_simp [hune, h1une]
    ring

private theorem logSquareTransform_continuousOn26 :
    ContinuousOn logSquareTransform26 (Icc (1 / 2 : ℝ) 1) := by
  have hlu :
      ContinuousOn Real.log (Icc (1 / 2 : ℝ) 1) := by
    apply continuousOn_id.log
    intro u hu
    simpa only [id_eq] using (ne_of_gt (show 0 < u by linarith [hu.1]))
  have hsub :
      ContinuousOn (fun u : ℝ => 1 - u) (Icc (1 / 2 : ℝ) 1) :=
    continuousOn_const.sub continuousOn_id
  have hsubmem :
      MapsTo (fun u : ℝ => 1 - u)
        (Icc (1 / 2 : ℝ) 1) (Icc (0 : ℝ) (1 / 2)) := by
    intro u hu
    constructor <;> linarith [hu.1, hu.2]
  have hsingular :
      ContinuousOn
        (fun u : ℝ => Real.log u ^ 2 * Real.log (1 - u))
        (Icc (1 / 2 : ℝ) 1) := by
    have hcomp :=
      log_mul_sq_log_one_sub_continuousOn26.comp hsub hsubmem
    convert hcomp using 1
    funext u
    simp [Function.comp_apply, mul_comm]
  have hd :
      ContinuousOn dilog (Icc (1 / 2 : ℝ) 1) :=
    dilog_continuousOn_unit.mono (by
      intro u hu
      constructor <;> linarith [hu.1, hu.2])
  have ht :
      ContinuousOn trilog26 (Icc (1 / 2 : ℝ) 1) :=
    trilog26_continuousOn_unit.mono (by
      intro u hu
      constructor <;> linarith [hu.1, hu.2])
  have hfirst :
      ContinuousOn
        (fun u : ℝ => (1 / 3 : ℝ) * Real.log u ^ 3)
        (Icc (1 / 2 : ℝ) 1) :=
    continuousOn_const.mul (hlu.pow 3)
  have hthird :
      ContinuousOn
        (fun u : ℝ => 2 * Real.log u * dilog u)
        (Icc (1 / 2 : ℝ) 1) :=
    (continuousOn_const.mul hlu).mul hd
  have hfourth :
      ContinuousOn
        (fun u : ℝ => 2 * trilog26 u)
        (Icc (1 / 2 : ℝ) 1) :=
    continuousOn_const.mul ht
  unfold logSquareTransform26
  have htotal := ((hfirst.sub hsingular).sub hthird).add hfourth
  simpa only [Pi.add_apply, Pi.sub_apply, Pi.mul_apply, mul_assoc] using htotal

private theorem logSquareTransform_one26 :
    logSquareTransform26 1 = 2 * zeta3 := by
  have htri1 : trilog26 1 = zeta3 := by
    unfold trilog26 zeta3
    apply tsum_congr
    intro n
    simp [Nat.cast_add]
  simp [logSquareTransform26, htri1]

private theorem logSquareTransform_half26 :
    logSquareTransform26 (1 / 2) = (7 : ℝ) / 4 * zeta3 := by
  have hloghalf : Real.log (1 / 2 : ℝ) = -Real.log 2 := by
    rw [one_div, Real.log_inv]
  rw [logSquareTransform26,
    show 1 - (1 / 2 : ℝ) = 1 / 2 by norm_num,
    hloghalf, dilog26_half, trilog26_half]
  ring

private def logSquareOnePlusPrimitive26 (x : ℝ) : ℝ :=
  2 * zeta3 - logSquareTransform26 (1 / (1 + x))

private theorem logSquareOnePlusPrimitive_hasDerivAt26
    {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasDerivAt logSquareOnePlusPrimitive26
      (Real.log (1 + x) ^ 2 / x) x := by
  have h1x0 : 0 < 1 + x := by linarith
  have h1xne : 1 + x ≠ 0 := ne_of_gt h1x0
  let u : ℝ := 1 / (1 + x)
  have hu0 : 0 < u := by dsimp [u]; positivity
  have hu1 : u < 1 := by
    dsimp [u]
    exact (div_lt_one h1x0).2 (by linarith)
  have huderiv :
      HasDerivAt (fun y : ℝ => 1 / (1 + y))
        (-1 / (1 + x) ^ 2) x := by
    have hden :
        HasDerivAt (fun y : ℝ => 1 + y) 1 x := by
      convert (hasDerivAt_const x 1).add (hasDerivAt_id x) using 1
      simp
    convert (hasDerivAt_const x 1).div hden h1xne using 1
    field_simp [h1xne]
    ring
  have hcomp :
      HasDerivAt (fun y : ℝ => logSquareTransform26 (1 / (1 + y)))
        ((Real.log u ^ 2 / (u * (1 - u))) *
          (-1 / (1 + x) ^ 2)) x := by
    simpa [u] using
      (logSquareTransform_hasDerivAt26 hu0 hu1).comp x huderiv
  unfold logSquareOnePlusPrimitive26
  convert (hcomp.const_sub (2 * zeta3)) using 1
  have hlogu : Real.log u = -Real.log (1 + x) := by
    dsimp [u]
    rw [one_div, Real.log_inv]
  rw [hlogu]
  dsimp [u]
  field_simp [h1xne, ne_of_gt hx0]
  ring

private theorem logSquareOnePlusPrimitive_continuousOn26 :
    ContinuousOn logSquareOnePlusPrimitive26 (Icc (0 : ℝ) 1) := by
  have hu :
      ContinuousOn (fun x : ℝ => 1 / (1 + x)) (Icc (0 : ℝ) 1) := by
    apply continuousOn_const.div (continuousOn_const.add continuousOn_id)
    intro x hx
    simpa only [Pi.add_apply, Pi.one_apply, id_eq] using
      (ne_of_gt (show 0 < 1 + x by linarith [hx.1]))
  have humem :
      MapsTo (fun x : ℝ => 1 / (1 + x))
        (Icc (0 : ℝ) 1) (Icc (1 / 2 : ℝ) 1) := by
    intro x hx
    have h1x0 : 0 < 1 + x := by linarith [hx.1]
    constructor
    · apply (le_div_iff₀ h1x0).2
      nlinarith [hx.2]
    · exact (div_le_one h1x0).2 (by linarith [hx.1])
  unfold logSquareOnePlusPrimitive26
  exact continuousOn_const.sub
    (logSquareTransform_continuousOn26.comp hu humem)

@[simp] private theorem logSquareOnePlusPrimitive_zero26 :
    logSquareOnePlusPrimitive26 0 = 0 := by
  simp [logSquareOnePlusPrimitive26, logSquareTransform_one26]

private theorem logSquareOnePlusPrimitive_one26 :
    logSquareOnePlusPrimitive26 1 = (1 : ℝ) / 4 * zeta3 := by
  rw [logSquareOnePlusPrimitive26]
  norm_num
  rw [logSquareTransform_half26]
  ring

private theorem log_mul_sq_log_one_add_continuousOn26 :
    ContinuousOn
      (fun x : ℝ => Real.log x * Real.log (1 + x) ^ 2)
      (Icc (0 : ℝ) 1) := by
  have hlog :
      ContinuousOn (fun x : ℝ => Real.log (1 + x))
        (Icc (0 : ℝ) 1) := by
    apply (continuousOn_const.add continuousOn_id).log
    intro x hx
    simpa only [Pi.add_apply, Pi.one_apply, id_eq] using
      (ne_of_gt (show 0 < 1 + x by linarith [hx.1]))
  have hright :
      ContinuousOn
        (fun x : ℝ =>
          (x * Real.log x) * logOnePlusSlope26 x * Real.log (1 + x))
        (Icc (0 : ℝ) 1) :=
    (Real.continuous_mul_log.continuousOn.mul
      logOnePlusSlope26_continuousOn).mul hlog
  apply hright.congr
  intro x hx
  by_cases hxzero : x = 0
  · subst x
    simp
  · simp [logOnePlusSlope26, hxzero]
    field_simp [hxzero]

private def alternatingWeightThreePrimitive26 (x : ℝ) : ℝ :=
  (1 / 2 : ℝ) * Real.log x * Real.log (1 + x) ^ 2 -
    (1 / 2 : ℝ) * logSquareOnePlusPrimitive26 x

private theorem alternatingWeightThreePrimitive_hasDerivAt26
    {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasDerivAt alternatingWeightThreePrimitive26
      (alternatingWeightThreeKernel26 x) x := by
  have hxne : x ≠ 0 := ne_of_gt hx0
  have h1xne : 1 + x ≠ 0 := by linarith
  have hlogx : HasDerivAt (fun y : ℝ => Real.log y) (1 / x) x := by
    simpa [one_div] using Real.hasDerivAt_log hxne
  have hlog1 :
      HasDerivAt (fun y : ℝ => Real.log (1 + y))
        (1 / (1 + x)) x := by
    have hinner :
        HasDerivAt (fun y : ℝ => 1 + y) 1 x := by
      convert (hasDerivAt_const x 1).add (hasDerivAt_id x) using 1
      simp
    convert hinner.log h1xne using 1
  have hK := logSquareOnePlusPrimitive_hasDerivAt26 hx0 hx1
  unfold alternatingWeightThreePrimitive26
  have htotal :=
    ((hlogx.mul (hlog1.pow 2)).const_mul (1 / 2 : ℝ)).sub
      (hK.const_mul (1 / 2 : ℝ))
  convert htotal using 1
  · funext y
    simp only [Pi.sub_apply, Pi.mul_apply, Pi.pow_apply]
    ring
  · unfold alternatingWeightThreeKernel26
    simp only [Pi.pow_apply]
    field_simp [hxne, h1xne]
    ring

private theorem alternatingWeightThreePrimitive_continuousOn26 :
    ContinuousOn alternatingWeightThreePrimitive26 (Icc (0 : ℝ) 1) := by
  unfold alternatingWeightThreePrimitive26
  have hfirst :
      ContinuousOn
        (fun x : ℝ =>
          (1 / 2 : ℝ) * (Real.log x * Real.log (1 + x) ^ 2))
        (Icc (0 : ℝ) 1) :=
    continuousOn_const.mul log_mul_sq_log_one_add_continuousOn26
  have hsecond :
      ContinuousOn
        (fun x : ℝ => (1 / 2 : ℝ) * logSquareOnePlusPrimitive26 x)
        (Icc (0 : ℝ) 1) :=
    continuousOn_const.mul logSquareOnePlusPrimitive_continuousOn26
  simpa only [mul_assoc] using hfirst.sub hsecond

@[simp] private theorem alternatingWeightThreePrimitive_zero26 :
    alternatingWeightThreePrimitive26 0 = 0 := by
  simp [alternatingWeightThreePrimitive26]

private theorem alternatingWeightThreePrimitive_one26 :
    alternatingWeightThreePrimitive26 1 = -(1 : ℝ) / 8 * zeta3 := by
  rw [alternatingWeightThreePrimitive26,
    logSquareOnePlusPrimitive_one26]
  simp
  ring

theorem alternatingWeightThreeIntegral26 :
    (∫ x : ℝ in 0..1, alternatingWeightThreeKernel26 x) =
      -(1 : ℝ) / 8 * zeta3 := by
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le
    (f := alternatingWeightThreePrimitive26)
    (f' := alternatingWeightThreeKernel26)
    (by norm_num)
    alternatingWeightThreePrimitive_continuousOn26
    (fun x hx =>
      alternatingWeightThreePrimitive_hasDerivAt26 hx.1 hx.2)
    alternatingWeightThreeKernel_intervalIntegrable26]
  rw [alternatingWeightThreePrimitive_one26,
    alternatingWeightThreePrimitive_zero26]
  ring

/-- The companion logarithmic-square integral.  This is useful for the
level-two Euler-sum reductions in Problem 2.4. -/
theorem logSquareOnePlusIntegral26 :
    (∫ x : ℝ in 0..1, Real.log (1 + x) ^ 2 / x) =
      (1 : ℝ) / 4 * zeta3 := by
  have hcont :
      ContinuousOn
        (fun x : ℝ => Real.log (1 + x) ^ 2 / x)
        (Icc (0 : ℝ) 1) := by
    have haux :
        ContinuousOn
          (fun x : ℝ => x * logOnePlusSlope26 x ^ 2)
          (Icc (0 : ℝ) 1) :=
      continuousOn_id.mul (logOnePlusSlope26_continuousOn.pow 2)
    apply haux.congr
    intro x hx
    by_cases hxzero : x = 0
    · subst x
      simp
    · simp [logOnePlusSlope26, hxzero]
      field_simp [hxzero]
  have hint :
      IntervalIntegrable
        (fun x : ℝ => Real.log (1 + x) ^ 2 / x)
        MeasureTheory.volume 0 1 := by
    apply ContinuousOn.intervalIntegrable
    simpa [Set.uIcc_of_le (show (0 : ℝ) ≤ 1 by norm_num)] using hcont
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le
    (f := logSquareOnePlusPrimitive26)
    (f' := fun x : ℝ => Real.log (1 + x) ^ 2 / x)
    (by norm_num)
    logSquareOnePlusPrimitive_continuousOn26
    (fun x hx => logSquareOnePlusPrimitive_hasDerivAt26 hx.1 hx.2)
    hint]
  rw [logSquareOnePlusPrimitive_one26,
    logSquareOnePlusPrimitive_zero26]
  ring

/-- The only remaining special value after the real Landen reduction: the
cyclotomic logarithmic integral at the cubic polynomial `1 + x + x²`. -/
def CyclotomicLogIntegralEvaluation26 : Prop :=
  (∫ x : ℝ in 0..1, cyclotomicWeightThreeKernel26 x) =
    -(7 : ℝ) / 36 * zeta3

theorem cyclotomicWeightThreeEvaluation26_of_cyclotomicLog
    (hCyclotomic : CyclotomicLogIntegralEvaluation26) :
    CyclotomicWeightThreeEvaluation26 := by
  unfold CyclotomicLogIntegralEvaluation26 at hCyclotomic
  unfold CyclotomicWeightThreeEvaluation26
  rw [hCyclotomic, alternatingWeightThreeIntegral26]
  ring

theorem nestedCyclotomicIntegral26_of_cyclotomicLog
    (hCyclotomic : CyclotomicLogIntegralEvaluation26) :
    (∫ x : ℝ in 0..1, nestedCyclotomicKernel26 x) =
      Real.pi ^ 2 / 18 - zeta3 / 3 :=
  nestedCyclotomicIntegral26_of_weightThree
    (cyclotomicWeightThreeEvaluation26_of_cyclotomicLog hCyclotomic)

end RamanujanChallenge.P26
