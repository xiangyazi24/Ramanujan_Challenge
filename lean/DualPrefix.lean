import RamanujanChallenge.Problem25Integral
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic
import Mathlib.MeasureTheory.Integral.Prod

noncomputable section

namespace RamanujanChallenge.P25

open MeasureTheory Set Filter Topology
open scoped Interval

set_option maxRecDepth 10000

private abbrev unitMeasure : Measure ℝ := volume.restrict (Ioc 0 1)

private abbrev cubeMeasure : Measure (ℝ × (ℝ × ℝ)) :=
  unitMeasure.prod (unitMeasure.prod unitMeasure)

example : IsFiniteMeasure unitMeasure := inferInstance
example : IsFiniteMeasure cubeMeasure := inferInstance

private def dualD (p q v : ℝ) : ℝ := p * q * (1 + v ^ 2) + 2 * v

private def rawMomentIntegrand
    (n A B C k : ℕ) (x : ℝ × (ℝ × ℝ)) : ℝ :=
  let p := x.1
  let q := x.2.1
  let v := x.2.2
  16 * p ^ (2 * n + 6 + A) * q ^ (2 * n + 5 + B) *
      (1 - p ^ 2) ^ n * (1 - q ^ 2) ^ (n + 1) *
      v ^ (2 * n + 3 + C) /
    dualD p q v ^ (2 * n + 4 + k)

private def boundedMomentIntegrand
    (n A B C k : ℕ) (x : ℝ × (ℝ × ℝ)) : ℝ :=
  let p := x.1
  let q := x.2.1
  let v := x.2.2
  (16 * 4 ^ (4 - k) / 2 ^ (2 * n + 3) : ℝ) *
    p ^ (2 * n + 1 + A) * q ^ (2 * n + B) * v ^ C *
    (1 - p ^ 2) ^ n * (1 - q ^ 2) ^ (n + 1) *
    (p * q / dualD p q v) ^ 5 *
    (2 * v / dualD p q v) ^ (2 * n + 3) *
    (dualD p q v / 4) ^ (4 - k)

private theorem cube_ae_bounds :
    ∀ᵐ x : ℝ × (ℝ × ℝ) ∂cubeMeasure,
      0 < x.1 ∧ x.1 ≤ 1 ∧
      0 < x.2.1 ∧ x.2.1 ≤ 1 ∧
      0 < x.2.2 ∧ x.2.2 ≤ 1 := by
  have hu : ∀ᵐ u : ℝ ∂unitMeasure, 0 < u ∧ u ≤ 1 :=
    ae_restrict_mem (measurableSet_Ioc : MeasurableSet (Ioc (0 : ℝ) 1))
  have hmem : ∀ᵐ x : ℝ × (ℝ × ℝ) ∂cubeMeasure,
      x ∈ Ioc (0 : ℝ) 1 ×ˢ (Ioc (0 : ℝ) 1 ×ˢ Ioc (0 : ℝ) 1) := by
    rw [Measure.ae_prod_mem_iff_ae_ae_mem
      (measurableSet_Ioc.prod (measurableSet_Ioc.prod measurableSet_Ioc))]
    filter_upwards [hu] with p hp
    have hqv : ∀ᵐ y : ℝ × ℝ ∂unitMeasure.prod unitMeasure,
        y ∈ Ioc (0 : ℝ) 1 ×ˢ Ioc (0 : ℝ) 1 := by
      rw [Measure.ae_prod_mem_iff_ae_ae_mem
        (measurableSet_Ioc.prod measurableSet_Ioc)]
      filter_upwards [hu] with q hq
      filter_upwards [hu] with v hv
      exact ⟨hq, hv⟩
    filter_upwards [hqv] with y hy
    exact ⟨hp, hy⟩
  filter_upwards [hmem] with x hx
  exact ⟨hx.1.1, hx.1.2, hx.2.1.1, hx.2.1.2, hx.2.2.1, hx.2.2.2⟩

private theorem rawMomentIntegrand_ae_eq_bounded
    (n A B C k : ℕ) (hk : k ≤ 4) :
    rawMomentIntegrand n A B C k =ᵐ[cubeMeasure]
      boundedMomentIntegrand n A B C k := by
  filter_upwards [cube_ae_bounds] with x hx
  rcases hx with ⟨hp, hp1, hq, hq1, hv, hv1⟩
  have hD : dualD x.1 x.2.1 x.2.2 ≠ 0 := by
    have : 0 < dualD x.1 x.2.1 x.2.2 := by
      dsimp [dualD]
      positivity
    exact this.ne'
  dsimp [rawMomentIntegrand, boundedMomentIntegrand]
  have hExp : (2 * n + 4 + k) + (4 - k) = 5 + (2 * n + 3) := by omega
  have hDpow :
      dualD x.1 x.2.1 x.2.2 ^ (2 * n + 4 + k) *
          dualD x.1 x.2.1 x.2.2 ^ (4 - k) =
        dualD x.1 x.2.1 x.2.2 ^ 5 *
          dualD x.1 x.2.1 x.2.2 ^ (2 * n + 3) := by
    rw [← pow_add, ← pow_add, hExp]
  rw [show 2 * n + 6 + A = (2 * n + 1 + A) + 5 by omega,
    show 2 * n + 5 + B = (2 * n + B) + 5 by omega,
    show 2 * n + 3 + C = C + (2 * n + 3) by omega,
    pow_add, pow_add, pow_add]
  simp only [div_pow]
  field_simp [hD, pow_ne_zero]
  rw [mul_pow]
  ring_nf at hDpow ⊢
  rw [hDpow]

private theorem boundedMomentIntegrand_measurable (n A B C k : ℕ) :
    Measurable (boundedMomentIntegrand n A B C k) := by
  unfold boundedMomentIntegrand dualD
  fun_prop

private theorem boundedMomentIntegrand_norm_le
    (n A B C k : ℕ) (hk : k ≤ 4) :
    ∀ᵐ x ∂cubeMeasure,
      ‖boundedMomentIntegrand n A B C k x‖ ≤
        ‖(16 * 4 ^ (4 - k) / 2 ^ (2 * n + 3) : ℝ)‖ := by
  filter_upwards [cube_ae_bounds] with x hx
  rcases hx with ⟨hp, hp1, hq, hq1, hv, hv1⟩
  have hp0 : 0 ≤ x.1 := hp.le
  have hq0 : 0 ≤ x.2.1 := hq.le
  have hv0 : 0 ≤ x.2.2 := hv.le
  have hp2 : x.1 ^ 2 ≤ 1 := by nlinarith [sq_nonneg (1 - x.1)]
  have hq2 : x.2.1 ^ 2 ≤ 1 := by nlinarith [sq_nonneg (1 - x.2.1)]
  have hDpos : 0 < dualD x.1 x.2.1 x.2.2 := by
    dsimp [dualD]
    positivity
  have hpqD : 0 ≤ x.1 * x.2.1 / dualD x.1 x.2.1 x.2.2 :=
    div_nonneg (mul_nonneg hp0 hq0) hDpos.le
  have hpqD1 : x.1 * x.2.1 / dualD x.1 x.2.1 x.2.2 ≤ 1 := by
    rw [div_le_one hDpos]
    dsimp [dualD]
    nlinarith [mul_nonneg (mul_nonneg hp0 hq0) (sq_nonneg x.2.2)]
  have hvD : 0 ≤ 2 * x.2.2 / dualD x.1 x.2.1 x.2.2 :=
    div_nonneg (mul_nonneg (by norm_num) hv0) hDpos.le
  have hvD1 : 2 * x.2.2 / dualD x.1 x.2.1 x.2.2 ≤ 1 := by
    rw [div_le_one hDpos]
    dsimp [dualD]
    nlinarith [mul_nonneg (mul_nonneg hp0 hq0) (by positivity : 0 ≤ 1 + x.2.2 ^ 2)]
  have hD0 : 0 ≤ dualD x.1 x.2.1 x.2.2 / 4 :=
    div_nonneg hDpos.le (by norm_num)
  have hD1 : dualD x.1 x.2.1 x.2.2 / 4 ≤ 1 := by
    rw [div_le_one (by norm_num : (0 : ℝ) < 4)]
    have hpq : x.1 * x.2.1 ≤ 1 := by
      nlinarith [mul_le_mul hp1 hq1 hq0 (by norm_num : (0 : ℝ) ≤ 1)]
    have hv2 : x.2.2 ^ 2 ≤ 1 := by nlinarith [sq_nonneg (1 - x.2.2)]
    have hterm : x.1 * x.2.1 * (1 + x.2.2 ^ 2) ≤ 2 := by
      calc
        x.1 * x.2.1 * (1 + x.2.2 ^ 2) ≤ 1 * 2 := by
          apply mul_le_mul hpq (by linarith) (by positivity) (by norm_num)
        _ = 2 := by norm_num
    dsimp [dualD]
    nlinarith
  have hnonneg : 0 ≤ boundedMomentIntegrand n A B C k x := by
    dsimp [boundedMomentIntegrand]
    have hpm : 0 ≤ 1 - x.1 ^ 2 := by linarith
    have hqm : 0 ≤ 1 - x.2.1 ^ 2 := by linarith
    positivity
  rw [Real.norm_eq_abs, abs_of_nonneg hnonneg]
  dsimp [boundedMomentIntegrand]
  have hconst : 0 ≤ (16 * 4 ^ (4 - k) / 2 ^ (2 * n + 3) : ℝ) := by positivity
  change _ ≤ |(16 * 4 ^ (4 - k) / 2 ^ (2 * n + 3) : ℝ)|
  rw [abs_of_nonneg hconst]
  have hpPow : x.1 ^ (2 * n + 1 + A) ≤ 1 := pow_le_one₀ hp0 hp1
  have hqPow : x.2.1 ^ (2 * n + B) ≤ 1 := pow_le_one₀ hq0 hq1
  have hvPow : x.2.2 ^ C ≤ 1 := pow_le_one₀ hv0 hv1
  have hpm0 : 0 ≤ 1 - x.1 ^ 2 := by linarith
  have hqm0 : 0 ≤ 1 - x.2.1 ^ 2 := by linarith
  have hpm1 : (1 - x.1 ^ 2) ^ n ≤ 1 :=
    pow_le_one₀ hpm0 (by nlinarith [sq_nonneg x.1])
  have hqm1 : (1 - x.2.1 ^ 2) ^ (n + 1) ≤ 1 :=
    pow_le_one₀ hqm0 (by nlinarith [sq_nonneg x.2.1])
  have hr1 : (x.1 * x.2.1 / dualD x.1 x.2.1 x.2.2) ^ 5 ≤ 1 :=
    pow_le_one₀ hpqD hpqD1
  have hr2 : (2 * x.2.2 / dualD x.1 x.2.1 x.2.2) ^ (2 * n + 3) ≤ 1 :=
    pow_le_one₀ hvD hvD1
  have hr3 : (dualD x.1 x.2.1 x.2.2 / 4) ^ (4 - k) ≤ 1 :=
    pow_le_one₀ hD0 hD1
  have h1 := mul_le_one₀ hpPow (pow_nonneg hq0 _) hqPow
  have h2 := mul_le_one₀ h1 (pow_nonneg hv0 _) hvPow
  have h3 := mul_le_one₀ h2 (pow_nonneg hpm0 _) hpm1
  have h4 := mul_le_one₀ h3 (pow_nonneg hqm0 _) hqm1
  have h5 := mul_le_one₀ h4 (pow_nonneg hpqD _) hr1
  have h6 := mul_le_one₀ h5 (pow_nonneg hvD _) hr2
  have h7 := mul_le_one₀ h6 (pow_nonneg hD0 _) hr3
  calc
    16 * 4 ^ (4 - k) / 2 ^ (2 * n + 3) * x.1 ^ (2 * n + 1 + A) *
          x.2.1 ^ (2 * n + B) * x.2.2 ^ C * (1 - x.1 ^ 2) ^ n *
          (1 - x.2.1 ^ 2) ^ (n + 1) *
          (x.1 * x.2.1 / dualD x.1 x.2.1 x.2.2) ^ 5 *
          (2 * x.2.2 / dualD x.1 x.2.1 x.2.2) ^ (2 * n + 3) *
          (dualD x.1 x.2.1 x.2.2 / 4) ^ (4 - k) =
        (16 * 4 ^ (4 - k) / 2 ^ (2 * n + 3) : ℝ) *
          (x.1 ^ (2 * n + 1 + A) * x.2.1 ^ (2 * n + B) * x.2.2 ^ C *
          (1 - x.1 ^ 2) ^ n * (1 - x.2.1 ^ 2) ^ (n + 1) *
          (x.1 * x.2.1 / dualD x.1 x.2.1 x.2.2) ^ 5 *
          (2 * x.2.2 / dualD x.1 x.2.1 x.2.2) ^ (2 * n + 3) *
          (dualD x.1 x.2.1 x.2.2 / 4) ^ (4 - k)) := by ring
    _ ≤ (16 * 4 ^ (4 - k) / 2 ^ (2 * n + 3) : ℝ) * 1 :=
      mul_le_mul_of_nonneg_left h7 hconst
    _ = (16 * 4 ^ (4 - k) / 2 ^ (2 * n + 3) : ℝ) := by ring

private theorem rawMomentIntegrand_integrable
    (n A B C k : ℕ) (hk : k ≤ 4) :
    Integrable (rawMomentIntegrand n A B C k) cubeMeasure := by
  rw [integrable_congr (rawMomentIntegrand_ae_eq_bounded n A B C k hk)]
  apply (integrable_const
    (μ := cubeMeasure) (‖(16 * 4 ^ (4 - k) / 2 ^ (2 * n + 3) : ℝ)‖)).mono'
  · exact (boundedMomentIntegrand_measurable n A B C k).aestronglyMeasurable
  · exact boundedMomentIntegrand_norm_le n A B C k hk

private def dualMoment (n A B C k : ℕ) : ℝ :=
  ∫ x, rawMomentIntegrand n A B C k x ∂cubeMeasure

private theorem rawMomentIntegrand_nonneg_ae (n A B C k : ℕ) :
    0 ≤ᵐ[cubeMeasure] rawMomentIntegrand n A B C k := by
  filter_upwards [cube_ae_bounds] with x hx
  rcases hx with ⟨hp, hp1, hq, hq1, hv, hv1⟩
  have hp2 : x.1 ^ 2 ≤ 1 := by nlinarith [sq_nonneg (1 - x.1)]
  have hq2 : x.2.1 ^ 2 ≤ 1 := by nlinarith [sq_nonneg (1 - x.2.1)]
  have hpm : 0 ≤ 1 - x.1 ^ 2 := by linarith
  have hqm : 0 ≤ 1 - x.2.1 ^ 2 := by linarith
  have hD : 0 < dualD x.1 x.2.1 x.2.2 := by
    dsimp [dualD]
    positivity
  dsimp [rawMomentIntegrand]
  positivity

private theorem subcube_measure_pos :
    0 < cubeMeasure
      (Ioc (1 / 2 : ℝ) (3 / 4) ×ˢ
        (Ioc (1 / 2 : ℝ) (3 / 4) ×ˢ Ioc (1 / 2 : ℝ) (3 / 4))) := by
  rw [Measure.prod_prod, Measure.prod_prod]
  have hmeas : MeasurableSet (Ioc (1 / 2 : ℝ) (3 / 4)) := measurableSet_Ioc
  have hsub : Ioc (1 / 2 : ℝ) (3 / 4) ⊆ Ioc 0 1 := by
    intro x hx
    constructor <;> linarith [hx.1, hx.2]
  rw [Measure.restrict_apply hmeas]
  rw [inter_eq_self_of_subset_left hsub, Real.volume_Ioc]
  norm_num

private theorem dualMoment_pos (n A B C k : ℕ) (hk : k ≤ 4) :
    0 < dualMoment n A B C k := by
  rw [dualMoment, MeasureTheory.integral_pos_iff_support_of_nonneg_ae
    (rawMomentIntegrand_nonneg_ae n A B C k)
    (rawMomentIntegrand_integrable n A B C k hk)]
  apply lt_of_lt_of_le subcube_measure_pos
  apply measure_mono
  intro x hx
  rcases hx with ⟨hp, hq, hv⟩
  have hp0 : 0 < x.1 := by linarith [hp.1]
  have hp1 : x.1 < 1 := by linarith [hp.2]
  have hq0 : 0 < x.2.1 := by linarith [hq.1]
  have hq1 : x.2.1 < 1 := by linarith [hq.2]
  have hv0 : 0 < x.2.2 := by linarith [hv.1]
  have hpm : 0 < 1 - x.1 ^ 2 := by nlinarith
  have hqm : 0 < 1 - x.2.1 ^ 2 := by nlinarith
  have hD : 0 < dualD x.1 x.2.1 x.2.2 := by
    dsimp [dualD]
    positivity
  change rawMomentIntegrand n A B C k x ≠ 0
  apply ne_of_gt
  dsimp [rawMomentIntegrand]
  positivity

private theorem rawMomentIntegrand_shift_two (n C k : ℕ) (x : ℝ × (ℝ × ℝ)) :
    rawMomentIntegrand n 2 2 C k x =
      rawMomentIntegrand n 0 0 C k x * x.1 ^ 2 * x.2.1 ^ 2 := by
  dsimp [rawMomentIntegrand]
  rw [show 2 * n + 6 + 2 = (2 * n + 6) + 2 by omega,
    show 2 * n + 5 + 2 = (2 * n + 5) + 2 by omega,
    pow_add, pow_add]
  ring

private theorem dualMoment_three_sub_shift_pos (n : ℕ) :
    0 < dualMoment n 0 0 3 3 - dualMoment n 2 2 3 3 := by
  let f : ℝ × (ℝ × ℝ) → ℝ := fun x =>
    rawMomentIntegrand n 0 0 3 3 x - rawMomentIntegrand n 2 2 3 3 x
  have hfint : Integrable f cubeMeasure :=
    (rawMomentIntegrand_integrable n 0 0 3 3 (by omega)).sub
      (rawMomentIntegrand_integrable n 2 2 3 3 (by omega))
  have hfnonneg : 0 ≤ᵐ[cubeMeasure] f := by
    filter_upwards [cube_ae_bounds] with x hx
    rcases hx with ⟨hp, hp1, hq, hq1, hv, hv1⟩
    change 0 ≤ rawMomentIntegrand n 0 0 3 3 x -
      rawMomentIntegrand n 2 2 3 3 x
    rw [rawMomentIntegrand_shift_two]
    have hraw : 0 ≤ rawMomentIntegrand n 0 0 3 3 x := by
      have hp2 : x.1 ^ 2 ≤ 1 := by nlinarith [sq_nonneg (1 - x.1)]
      have hq2 : x.2.1 ^ 2 ≤ 1 := by nlinarith [sq_nonneg (1 - x.2.1)]
      have hpm : 0 ≤ 1 - x.1 ^ 2 := by linarith
      have hqm : 0 ≤ 1 - x.2.1 ^ 2 := by linarith
      have hD : 0 < dualD x.1 x.2.1 x.2.2 := by
        dsimp [dualD]
        positivity
      dsimp [rawMomentIntegrand]
      positivity
    have hpq : x.1 ^ 2 * x.2.1 ^ 2 ≤ 1 := by
      nlinarith [mul_le_mul (show x.1 ^ 2 ≤ 1 by nlinarith [sq_nonneg (1 - x.1)])
        (show x.2.1 ^ 2 ≤ 1 by nlinarith [sq_nonneg (1 - x.2.1)])
        (sq_nonneg x.2.1) (by norm_num : (0 : ℝ) ≤ 1)]
    nlinarith
  have hfpos : 0 < ∫ x, f x ∂cubeMeasure := by
    rw [MeasureTheory.integral_pos_iff_support_of_nonneg_ae hfnonneg hfint]
    apply lt_of_lt_of_le subcube_measure_pos
    apply measure_mono
    intro x hx
    rcases hx with ⟨hp, hq, hv⟩
    change f x ≠ 0
    apply ne_of_gt
    dsimp [f]
    rw [rawMomentIntegrand_shift_two]
    have hp0 : 0 < x.1 := by linarith [hp.1]
    have hp1 : x.1 < 1 := by linarith [hp.2]
    have hq0 : 0 < x.2.1 := by linarith [hq.1]
    have hq1 : x.2.1 < 1 := by linarith [hq.2]
    have hv0 : 0 < x.2.2 := by linarith [hv.1]
    have hpm : 0 < 1 - x.1 ^ 2 := by nlinarith
    have hqm : 0 < 1 - x.2.1 ^ 2 := by nlinarith
    have hD : 0 < dualD x.1 x.2.1 x.2.2 := by
      dsimp [dualD]
      positivity
    have hraw : 0 < rawMomentIntegrand n 0 0 3 3 x := by
      dsimp [rawMomentIntegrand]
      positivity
    have hpq : x.1 ^ 2 * x.2.1 ^ 2 < 1 := by
      have hp2 : x.1 ^ 2 < 1 := by nlinarith [sq_nonneg x.1]
      have hq2 : x.2.1 ^ 2 ≤ 1 := by nlinarith [sq_nonneg (1 - x.2.1)]
      nlinarith [mul_lt_mul_of_pos_right hp2 (sq_pos_of_pos hq0),
        mul_le_mul_of_nonneg_left hq2 (sq_nonneg x.1)]
    nlinarith
  rw [dualMoment, dualMoment, ← MeasureTheory.integral_sub
    (rawMomentIntegrand_integrable n 0 0 3 3 (by omega))
    (rawMomentIntegrand_integrable n 2 2 3 3 (by omega))]
  exact hfpos


end RamanujanChallenge.P25
end
