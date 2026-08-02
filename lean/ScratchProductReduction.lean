import RamanujanChallenge.Problem25DualIntegral

noncomputable section

namespace RamanujanChallenge.P25

open MeasureTheory Set
open scoped Interval

private theorem unit_integral_eq_interval (f : ℝ → ℝ) :
    (∫ x, f x ∂unitMeasure) = ∫ x in (0 : ℝ)..1, f x := by
  rw [intervalIntegral.integral_of_le (by norm_num : (0 : ℝ) ≤ 1)]

private theorem unit_integral_indicator_Ici_eq_interval
    (f : ℝ → ℝ) {u : ℝ} (hu0 : 0 < u) (hu1 : u ≤ 1) :
    (∫ p, Set.indicator {p : ℝ | u ≤ p} f p ∂unitMeasure) =
      ∫ p in u..1, f p := by
  change (∫ p, (Ici u).indicator f p ∂unitMeasure) = _
  rw [MeasureTheory.integral_indicator measurableSet_Ici]
  change (∫ p, f p ∂(volume.restrict (Ioc (0 : ℝ) 1)).restrict (Ici u)) = _
  rw [Measure.restrict_restrict measurableSet_Ici]
  have hset : Ici u ∩ Ioc (0 : ℝ) 1 = Icc u 1 := by
    ext p
    simp only [mem_inter_iff, mem_Ioc, mem_Ici, mem_Icc]
    constructor
    · intro hp
      exact ⟨hp.1, hp.2.2⟩
    · intro hp
      exact ⟨hp.1, ⟨hu0.trans_le hp.1, hp.2⟩⟩
  rw [hset, ← Measure.restrict_congr_set Ioc_ae_eq_Icc]
  rw [intervalIntegral.integral_of_le hu1]

private def oneSubSquarePrimitive (u p : ℝ) : ℝ :=
  p + u ^ 2 / p

private theorem oneSubSquarePrimitive_hasDerivAt
    (u p : ℝ) (hp : p ≠ 0) :
    HasDerivAt (oneSubSquarePrimitive u)
      (1 - u ^ 2 / p ^ 2) p := by
  have h := (hasDerivAt_id p).add
    ((hasDerivAt_const p (u ^ 2)).div (hasDerivAt_id p) hp)
  convert h using 1 <;> simp [oneSubSquarePrimitive, id_eq] <;>
    field_simp [hp] <;> ring

private theorem integral_one_sub_square_div
    {u : ℝ} (hu0 : 0 < u) (hu1 : u ≤ 1) :
    (∫ p in u..1, (1 - u ^ 2 / p ^ 2)) = (1 - u) ^ 2 := by
  have hcont : ContinuousOn (oneSubSquarePrimitive u) (Icc u 1) := by
    unfold oneSubSquarePrimitive
    apply ContinuousOn.add continuousOn_id
    apply ContinuousOn.div continuousOn_const continuousOn_id
    intro p hp
    change p ≠ 0
    linarith [hp.1]
  have hint : IntervalIntegrable (fun p : ℝ => 1 - u ^ 2 / p ^ 2)
      volume u 1 := by
    apply ContinuousOn.intervalIntegrable
    apply ContinuousOn.sub continuousOn_const
    apply ContinuousOn.div continuousOn_const
    · fun_prop
    · intro p hp
      have hp' : p ≠ 0 := by
        simp only [uIcc_of_le hu1, mem_Icc] at hp
        linarith [hp.1]
      exact pow_ne_zero 2 hp'
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le hu1 hcont
    (fun p hp => oneSubSquarePrimitive_hasDerivAt u p (by linarith [hp.1])) hint]
  unfold oneSubSquarePrimitive
  field_simp [hu0.ne']
  ring

private def triangleKernel (H : ℝ → ℝ) (p u : ℝ) : ℝ :=
  Set.indicator {u : ℝ | u ≤ p}
    (fun u => u ^ 5 * (1 - u ^ 2 / p ^ 2) * H u) u

theorem triangleKernel_integrable_of_bounded
    (H : ℝ → ℝ) (hH : Measurable H) (C : ℝ) (hC : 0 ≤ C)
    (hbound : ∀ u ∈ Ioc (0 : ℝ) 1, ‖H u‖ ≤ C) :
    Integrable (fun x : ℝ × ℝ => triangleKernel H x.1 x.2)
      (unitMeasure.prod unitMeasure) := by
  have hbase : Measurable (fun x : ℝ × ℝ =>
      x.2 ^ 5 * (1 - x.2 ^ 2 / x.1 ^ 2) * H x.2) := by
    fun_prop
  have hmeas : Measurable (fun x : ℝ × ℝ => triangleKernel H x.1 x.2) := by
    change Measurable (({x : ℝ × ℝ | x.2 ≤ x.1}).indicator
      (fun x => x.2 ^ 5 * (1 - x.2 ^ 2 / x.1 ^ 2) * H x.2))
    exact hbase.indicator (measurableSet_le measurable_snd measurable_fst)
  apply (integrable_const C).mono' hmeas.aestronglyMeasurable
  have hmem : ∀ᵐ x : ℝ × ℝ ∂unitMeasure.prod unitMeasure,
      x ∈ Ioc (0 : ℝ) 1 ×ˢ Ioc (0 : ℝ) 1 := by
    rw [Measure.ae_prod_mem_iff_ae_ae_mem
      (measurableSet_Ioc.prod measurableSet_Ioc)]
    filter_upwards [unit_ae_bounds] with p hp
    filter_upwards [unit_ae_bounds] with u hu
    exact ⟨hp, hu⟩
  filter_upwards [hmem] with x hx
  rcases hx with ⟨hp, hu⟩
  by_cases hle : x.2 ≤ x.1
  · rw [triangleKernel, Set.indicator_of_mem
      (show x.2 ∈ {u : ℝ | u ≤ x.1} by exact hle)]
    have hp2 : 0 < x.1 ^ 2 := sq_pos_of_pos hp.1
    have hu2 : 0 ≤ x.2 ^ 2 := sq_nonneg x.2
    have hsq : x.2 ^ 2 ≤ x.1 ^ 2 := by
      simpa only [pow_two] using mul_self_le_mul_self hu.1.le hle
    have hratio0 : 0 ≤ x.2 ^ 2 / x.1 ^ 2 := div_nonneg hu2 hp2.le
    have hratio1 : x.2 ^ 2 / x.1 ^ 2 ≤ 1 := (div_le_one hp2).2 hsq
    have hfactor0 : 0 ≤ 1 - x.2 ^ 2 / x.1 ^ 2 := by linarith
    have hfactor1 : 1 - x.2 ^ 2 / x.1 ^ 2 ≤ 1 := by linarith
    have hu50 : 0 ≤ x.2 ^ 5 := pow_nonneg hu.1.le _
    have hu5 : x.2 ^ 5 ≤ 1 := pow_le_one₀ hu.1.le hu.2
    have hprod : x.2 ^ 5 * (1 - x.2 ^ 2 / x.1 ^ 2) ≤ 1 :=
      mul_le_one₀ hu5 hfactor0 hfactor1
    have hHnorm := hbound x.2 hu
    rw [Real.norm_eq_abs] at hHnorm ⊢
    rw [abs_mul, abs_mul, abs_of_nonneg hu50,
      abs_of_nonneg hfactor0]
    simpa using mul_le_mul hprod hHnorm (abs_nonneg _) (by norm_num : (0 : ℝ) ≤ 1)
  · rw [triangleKernel, Set.indicator_of_notMem
      (show x.2 ∉ {u : ℝ | u ≤ x.1} by exact hle)]
    simp [hC]

private theorem triangleKernel_inner_integral (H : ℝ → ℝ)
    {u : ℝ} (hu0 : 0 < u) (hu1 : u ≤ 1) :
    (∫ p, triangleKernel H p u ∂unitMeasure) =
      u ^ 5 * (1 - u) ^ 2 * H u := by
  rw [show (fun p => triangleKernel H p u) =
      fun p => Set.indicator {p : ℝ | u ≤ p}
        (fun p => u ^ 5 * (1 - u ^ 2 / p ^ 2) * H u) p by
    funext p
    simp only [triangleKernel]
    rfl]
  rw [unit_integral_indicator_Ici_eq_interval _ hu0 hu1]
  have hfun : (fun p : ℝ => u ^ 5 * (1 - u ^ 2 / p ^ 2) * H u) =
      fun p => (u ^ 5 * H u) * (1 - u ^ 2 / p ^ 2) := by
    funext p
    ring
  rw [hfun, intervalIntegral.integral_const_mul,
    integral_one_sub_square_div hu0 hu1]
  ring

private theorem product_inner_eq_triangle (H : ℝ → ℝ)
    {p : ℝ} (hp0 : 0 < p) (hp1 : p ≤ 1) :
    (∫ q, p ^ 6 * q ^ 5 * (1 - q ^ 2) * H (p * q) ∂unitMeasure) =
      ∫ u, triangleKernel H p u ∂unitMeasure := by
  rw [unit_integral_eq_interval, unit_integral_eq_interval]
  let F : ℝ → ℝ := fun u => u ^ 5 * (1 - u ^ 2 / p ^ 2) * H u
  calc
    (∫ q in (0 : ℝ)..1, p ^ 6 * q ^ 5 * (1 - q ^ 2) * H (p * q)) =
        p * ∫ q in (0 : ℝ)..1, F (q * p) := by
          rw [← intervalIntegral.integral_const_mul]
          apply intervalIntegral.integral_congr
          intro q hq
          dsimp [F]
          rw [mul_comm q p]
          field_simp [hp0.ne']
    _ = ∫ u in (0 : ℝ) * p..1 * p, F u :=
      intervalIntegral.mul_integral_comp_mul_right p
    _ = ∫ u in (0 : ℝ)..p, F u := by norm_num
    _ = ∫ u in (0 : ℝ)..1, triangleKernel H p u := by
      change (∫ u in (0 : ℝ)..p, F u) =
        ∫ u in (0 : ℝ)..1, Set.indicator {u : ℝ | u ≤ p} F u
      exact (intervalIntegral.integral_indicator
        (f := F) (μ := volume) ⟨hp0.le, hp1⟩).symm

theorem product_moment_reduction (H : ℝ → ℝ)
    (htri : Integrable (fun x : ℝ × ℝ => triangleKernel H x.1 x.2)
      (unitMeasure.prod unitMeasure)) :
    (∫ p, ∫ q, p ^ 6 * q ^ 5 * (1 - q ^ 2) * H (p * q)
          ∂unitMeasure ∂unitMeasure) =
      ∫ u, u ^ 5 * (1 - u) ^ 2 * H u ∂unitMeasure := by
  calc
    _ = ∫ p, ∫ u, triangleKernel H p u ∂unitMeasure ∂unitMeasure := by
      apply integral_congr_ae
      filter_upwards [unit_ae_bounds] with p hp
      exact product_inner_eq_triangle H hp.1 hp.2
    _ = ∫ u, ∫ p, triangleKernel H p u ∂unitMeasure ∂unitMeasure := by
      exact MeasureTheory.integral_integral_swap htri
    _ = _ := by
      apply integral_congr_ae
      filter_upwards [unit_ae_bounds] with u hu
      exact triangleKernel_inner_integral H hu.1 hu.2

end RamanujanChallenge.P25

end
