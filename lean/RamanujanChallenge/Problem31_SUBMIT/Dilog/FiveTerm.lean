import Ramanujan31.Dilog.Rogers

/-!
# Abel's five-term relation for the real Rogers dilogarithm

The proof is elementary on the open unit square.  It differentiates the
five-term defect in its second variable, proves the derivative vanishes by
logarithm algebra, and determines the constant from the limit at zero.
-/

open Set Filter
open scoped Topology

namespace Real

private noncomputable def rogersDerivative (z : ℝ) : ℝ :=
  (-log z / (1 - z) - log (1 - z) / z) / 2

private noncomputable def abelB (x y : ℝ) : ℝ :=
  x * (1 - y) / (1 - x * y)

private noncomputable def abelC (x y : ℝ) : ℝ :=
  y * (1 - x) / (1 - x * y)

private theorem abel_arguments
    {x y : ℝ} (hx0 : 0 < x) (hx1 : x < 1) (hy0 : 0 < y) (hy1 : y < 1) :
    0 < x * y ∧ x * y < 1 ∧
      0 < abelB x y ∧ abelB x y < 1 ∧
      0 < abelC x y ∧ abelC x y < 1 := by
  have hxy0 : 0 < x * y := mul_pos hx0 hy0
  have hxy1 : x * y < 1 :=
    mul_lt_one_of_nonneg_of_lt_one_right hx1.le hy0.le hy1
  have hD : 0 < 1 - x * y := sub_pos.mpr hxy1
  have h1x : 0 < 1 - x := sub_pos.mpr hx1
  have h1y : 0 < 1 - y := sub_pos.mpr hy1
  have hB0 : 0 < abelB x y := by
    exact div_pos (mul_pos hx0 h1y) hD
  have hB1 : abelB x y < 1 := by
    rw [abelB, div_lt_one hD]
    nlinarith
  have hC0 : 0 < abelC x y := by
    exact div_pos (mul_pos hy0 h1x) hD
  have hC1 : abelC x y < 1 := by
    rw [abelC, div_lt_one hD]
    nlinarith
  exact ⟨hxy0, hxy1, hB0, hB1, hC0, hC1⟩

private theorem rogersDerivative_fiveTerm
    {x y : ℝ} (hx0 : 0 < x) (hx1 : x < 1) (hy0 : 0 < y) (hy1 : y < 1) :
    rogersDerivative y
      - x * rogersDerivative (x * y)
      - (-x * (1 - x) / (1 - x * y) ^ 2) * rogersDerivative (abelB x y)
      - ((1 - x) / (1 - x * y) ^ 2) * rogersDerivative (abelC x y) = 0 := by
  have hxne : x ≠ 0 := ne_of_gt hx0
  have hyne : y ≠ 0 := ne_of_gt hy0
  have h1x : 0 < 1 - x := sub_pos.mpr hx1
  have h1y : 0 < 1 - y := sub_pos.mpr hy1
  have h1xne : 1 - x ≠ 0 := ne_of_gt h1x
  have h1yne : 1 - y ≠ 0 := ne_of_gt h1y
  have hxy1 : x * y < 1 :=
    mul_lt_one_of_nonneg_of_lt_one_right hx1.le hy0.le hy1
  have hD : 0 < 1 - x * y := sub_pos.mpr hxy1
  have hDne : 1 - x * y ≠ 0 := ne_of_gt hD
  have hDne' : 1 - y * x ≠ 0 := by
    simpa [mul_comm] using hDne
  have hB0 := (abel_arguments hx0 hx1 hy0 hy1).2.2.1
  have hC0 := (abel_arguments hx0 hx1 hy0 hy1).2.2.2.2.1
  have hBne : abelB x y ≠ 0 := ne_of_gt hB0
  have hCne : abelC x y ≠ 0 := ne_of_gt hC0
  have h1B :
      1 - abelB x y = (1 - x) / (1 - x * y) := by
    rw [abelB]
    field_simp
    ring
  have h1C :
      1 - abelC x y = (1 - y) / (1 - x * y) := by
    rw [abelC]
    field_simp [hDne, hDne']
    ring
  have h1Bne : 1 - abelB x y ≠ 0 := by
    rw [h1B]
    exact div_ne_zero h1xne hDne
  have h1Cne : 1 - abelC x y ≠ 0 := by
    rw [h1C]
    exact div_ne_zero h1yne hDne
  unfold rogersDerivative
  rw [h1B, h1C]
  unfold abelB abelC
  rw [Real.log_mul hxne hyne]
  rw [Real.log_div (mul_ne_zero hxne h1yne) hDne,
    Real.log_mul hxne h1yne]
  rw [Real.log_div (mul_ne_zero hyne h1xne) hDne,
    Real.log_mul hyne h1xne]
  rw [Real.log_div h1xne hDne, Real.log_div h1yne hDne]
  field_simp [hDne, hDne']
  ring

private theorem hasDerivAt_abelB
    {x y : ℝ} (hxy : 1 - x * y ≠ 0) :
    HasDerivAt (abelB x)
      (-x * (1 - x) / (1 - x * y) ^ 2) y := by
  have hnum :=
    (hasDerivAt_const y x).mul
      ((hasDerivAt_const y 1).sub (hasDerivAt_id y))
  have hden :=
    (hasDerivAt_const y 1).sub
      ((hasDerivAt_const y x).mul (hasDerivAt_id y))
  unfold abelB
  convert hnum.div hden hxy using 1
  simp
  field_simp [hxy]
  ring

private theorem hasDerivAt_abelC
    {x y : ℝ} (hxy : 1 - x * y ≠ 0) :
    HasDerivAt (abelC x)
      ((1 - x) / (1 - x * y) ^ 2) y := by
  have hnum :=
    (hasDerivAt_id y).mul
      ((hasDerivAt_const y 1).sub (hasDerivAt_const y x))
  have hden :=
    (hasDerivAt_const y 1).sub
      ((hasDerivAt_const y x).mul (hasDerivAt_id y))
  unfold abelC
  convert hnum.div hden hxy using 1
  simp
  field_simp [hxy]
  ring

private noncomputable def abelDefect (x y : ℝ) : ℝ :=
  rogers x + rogers y - rogers (x * y) - rogers (abelB x y) - rogers (abelC x y)

private theorem hasDerivAt_abelDefect
    {x y : ℝ} (hx0 : 0 < x) (hx1 : x < 1) (hy0 : 0 < y) (hy1 : y < 1) :
    HasDerivAt (abelDefect x) 0 y := by
  obtain ⟨hxy0, hxy1, hB0, hB1, hC0, hC1⟩ :=
    abel_arguments hx0 hx1 hy0 hy1
  have hDne : 1 - x * y ≠ 0 := ne_of_gt (sub_pos.mpr hxy1)
  have hxy' : HasDerivAt (fun z : ℝ => x * z) x y := by
    convert (hasDerivAt_const y x).mul (hasDerivAt_id y) using 1
    ring
  have hRy :
      HasDerivAt rogers (rogersDerivative y) y := by
    simpa [rogersDerivative] using hasDerivAt_rogers hy0 hy1
  have hRxy :
      HasDerivAt (fun z : ℝ => rogers (x * z))
        (rogersDerivative (x * y) * x) y :=
    by
      have hbase :
          HasDerivAt rogers (rogersDerivative (x * y)) (x * y) := by
        simpa [rogersDerivative] using hasDerivAt_rogers hxy0 hxy1
      exact hbase.comp y hxy'
  have hRB :
      HasDerivAt (fun z : ℝ => rogers (abelB x z))
        (rogersDerivative (abelB x y) *
          (-x * (1 - x) / (1 - x * y) ^ 2)) y :=
    by
      have hbase :
          HasDerivAt rogers (rogersDerivative (abelB x y)) (abelB x y) := by
        simpa [rogersDerivative] using hasDerivAt_rogers hB0 hB1
      exact hbase.comp y (hasDerivAt_abelB hDne)
  have hRC :
      HasDerivAt (fun z : ℝ => rogers (abelC x z))
        (rogersDerivative (abelC x y) *
          ((1 - x) / (1 - x * y) ^ 2)) y :=
    by
      have hbase :
          HasDerivAt rogers (rogersDerivative (abelC x y)) (abelC x y) := by
        simpa [rogersDerivative] using hasDerivAt_rogers hC0 hC1
      exact hbase.comp y (hasDerivAt_abelC hDne)
  have h :=
    (((hasDerivAt_const y (rogers x)).add hRy).sub hRxy).sub hRB |>.sub hRC
  unfold abelDefect
  convert h using 1
  simpa [mul_comm] using
    (rogersDerivative_fiveTerm hx0 hx1 hy0 hy1).symm

private theorem tendsto_rogers_zero :
    Tendsto rogers (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  have hfilter0 :
      𝓝[>] (0 : ℝ) ≤ 𝓝[Icc 0 1] (0 : ℝ) :=
    nhdsWithin_le_iff.mpr (Icc_mem_nhdsGT zero_lt_one)
  have hdilog0 :
      Tendsto dilog (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    simpa using
      (continuousOn_dilog_Icc.continuousWithinAt
        (by norm_num : (0 : ℝ) ∈ Icc 0 1)).mono_left hfilter0
  have hcorrection :
      Tendsto (fun z : ℝ => (1 / 2 : ℝ) * (log z * log (1 - z)))
        (𝓝[>] 0) (𝓝 0) := by
    simpa using
      (tendsto_const_nhds.mul tendsto_log_mul_log_one_sub_zero)
  unfold rogers
  simpa [mul_assoc] using hdilog0.add hcorrection

private theorem tendsto_abelDefect_zero
    {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    Tendsto (abelDefect x) (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  have hto0 :
      Tendsto (fun y : ℝ => y) (𝓝[>] (0 : ℝ)) (𝓝 0) :=
    (tendsto_id :
      Tendsto (fun y : ℝ => y) (𝓝 0) (𝓝 0)).mono_left inf_le_left
  have hxy :
      Tendsto (fun y : ℝ => x * y) (𝓝[>] 0) (𝓝 0) := by
    simpa using tendsto_const_nhds.mul hto0
  have hxyPos :
      Tendsto (fun y : ℝ => x * y) (𝓝[>] 0) (𝓝[>] 0) := by
    refine tendsto_nhdsWithin_iff.mpr ⟨hxy, ?_⟩
    filter_upwards [self_mem_nhdsWithin] with y hy
    change 0 < y at hy
    exact mul_pos hx0 hy
  have hD :
      Tendsto (fun y : ℝ => 1 - x * y) (𝓝[>] 0) (𝓝 1) := by
    simpa using tendsto_const_nhds.sub hxy
  have hOneSub :
      Tendsto (fun y : ℝ => 1 - y) (𝓝[>] 0) (𝓝 1) := by
    simpa using tendsto_const_nhds.sub hto0
  have hB :
      Tendsto (abelB x) (𝓝[>] 0) (𝓝 x) := by
    unfold abelB
    simpa using (tendsto_const_nhds.mul hOneSub).div hD one_ne_zero
  have hC :
      Tendsto (abelC x) (𝓝[>] 0) (𝓝 0) := by
    unfold abelC
    simpa using (hto0.mul tendsto_const_nhds).div hD one_ne_zero
  have hCPos :
      Tendsto (abelC x) (𝓝[>] 0) (𝓝[>] 0) := by
    refine tendsto_nhdsWithin_iff.mpr ⟨hC, ?_⟩
    filter_upwards
      [self_mem_nhdsWithin,
        mem_inf_of_left (Iio_mem_nhds (show (0 : ℝ) < 1 by norm_num))]
      with y hy0 hy1
    change 0 < y at hy0
    change y < 1 at hy1
    exact (abel_arguments hx0 hx1 hy0 hy1).2.2.2.2.1
  have hRB :
      Tendsto (fun y : ℝ => rogers (abelB x y))
        (𝓝[>] 0) (𝓝 (rogers x)) := by
    have hcont : Tendsto rogers (𝓝 x) (𝓝 (rogers x)) :=
      (hasDerivAt_rogers hx0 hx1).continuousAt
    exact hcont.comp hB
  have hRxy :
      Tendsto (fun y : ℝ => rogers (x * y)) (𝓝[>] 0) (𝓝 0) :=
    tendsto_rogers_zero.comp hxyPos
  have hRC :
      Tendsto (fun y : ℝ => rogers (abelC x y)) (𝓝[>] 0) (𝓝 0) :=
    tendsto_rogers_zero.comp hCPos
  have hRx :
      Tendsto (fun _ : ℝ => rogers x) (𝓝[>] 0) (𝓝 (rogers x)) :=
    tendsto_const_nhds
  unfold abelDefect
  simpa using (((hRx.add tendsto_rogers_zero).sub hRxy).sub hRB).sub hRC

/-- Abel's five-term relation for the real Rogers dilogarithm. -/
theorem rogers_five_term
    {x y : ℝ} (hx0 : 0 < x) (hx1 : x < 1) (hy0 : 0 < y) (hy1 : y < 1) :
    rogers x + rogers y =
      rogers (x * y)
        + rogers (x * (1 - y) / (1 - x * y))
        + rogers (y * (1 - x) / (1 - x * y)) := by
  have hdiff : DifferentiableOn ℝ (abelDefect x) (Ioo 0 1) := by
    intro z hz
    exact (hasDerivAt_abelDefect hx0 hx1 hz.1 hz.2).differentiableAt
      |>.differentiableWithinAt
  have hzero : (Ioo (0 : ℝ) 1).EqOn (deriv (abelDefect x)) 0 := by
    intro z hz
    exact (hasDerivAt_abelDefect hx0 hx1 hz.1 hz.2).deriv
  have hconst :
      ∀ ⦃z : ℝ⦄, z ∈ Ioo 0 1 → abelDefect x y = abelDefect x z := by
    intro z hz
    exact isOpen_Ioo.is_const_of_deriv_eq_zero
      isPreconnected_Ioo hdiff hzero ⟨hy0, hy1⟩ hz
  have hevent :
      abelDefect x =ᶠ[𝓝[>] (0 : ℝ)] fun _ => abelDefect x y := by
    filter_upwards
      [self_mem_nhdsWithin,
        mem_inf_of_left (Iio_mem_nhds (show (0 : ℝ) < 1 by norm_num))]
      with z hz0 hz1
    change 0 < z at hz0
    change z < 1 at hz1
    exact (hconst ⟨hz0, hz1⟩).symm
  have hlim :
      Tendsto (abelDefect x) (𝓝[>] 0) (𝓝 (abelDefect x y)) :=
    tendsto_const_nhds.congr' hevent.symm
  have heq : abelDefect x y = 0 :=
    tendsto_nhds_unique hlim (tendsto_abelDefect_zero hx0 hx1)
  unfold abelDefect abelB abelC at heq
  linarith

end Real
