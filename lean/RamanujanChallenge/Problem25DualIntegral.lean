import RamanujanChallenge.Problem25Integral
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic
import Mathlib.MeasureTheory.Integral.Prod

noncomputable section

namespace RamanujanChallenge.P25

open MeasureTheory Set Filter Topology
open scoped Interval

set_option maxRecDepth 10000

abbrev unitMeasure : Measure ℝ := volume.restrict (Ioc 0 1)

abbrev cubeMeasure : Measure (ℝ × (ℝ × ℝ)) :=
  unitMeasure.prod (unitMeasure.prod unitMeasure)

example : IsFiniteMeasure unitMeasure := inferInstance
example : IsFiniteMeasure cubeMeasure := inferInstance

theorem unit_ae_bounds :
    ∀ᵐ u : ℝ ∂unitMeasure, 0 < u ∧ u ≤ 1 :=
  ae_restrict_mem (measurableSet_Ioc : MeasurableSet (Ioc (0 : ℝ) 1))

def dualD (p q v : ℝ) : ℝ := p * q * (1 + v ^ 2) + 2 * v

def rawMomentIntegrand
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
  have hScalar :
      (16 * (4 : ℝ) ^ (4 - k) / 2 ^ (2 * n + 3)) *
          2 ^ (2 * n + 3) * (1 / 4) ^ (4 - k) = 16 := by
    field_simp
    rw [← mul_pow]
    norm_num
  have hDen :
      (dualD x.1 x.2.1 x.2.2 ^ 5)⁻¹ *
            (dualD x.1 x.2.1 x.2.2 ^ (2 * n + 3))⁻¹ *
            dualD x.1 x.2.1 x.2.2 ^ (4 - k) =
          (dualD x.1 x.2.1 x.2.2 ^ (2 * n + 4 + k))⁻¹ := by
    field_simp [hD]
    rw [← pow_add, ← pow_add]
    congr 1
    omega
  rw [show 2 * n + 6 + A = (2 * n + 1 + A) + 5 by omega,
    show 2 * n + 5 + B = (2 * n + B) + 5 by omega,
    show 2 * n + 3 + C = C + (2 * n + 3) by omega,
    pow_add, pow_add, pow_add]
  simp only [div_eq_mul_inv, mul_pow]
  calc
    _ =
        ((16 * (4 : ℝ) ^ (4 - k) / 2 ^ (2 * n + 3)) *
            2 ^ (2 * n + 3) * (1 / 4) ^ (4 - k)) *
          x.1 ^ (2 * n + 1 + A) * x.2.1 ^ (2 * n + B) *
          x.2.2 ^ C * (1 - x.1 ^ 2) ^ n *
          (1 - x.2.1 ^ 2) ^ (n + 1) * x.1 ^ 5 * x.2.1 ^ 5 *
          x.2.2 ^ (2 * n + 3) *
          ((dualD x.1 x.2.1 x.2.2 ^ 5)⁻¹ *
            (dualD x.1 x.2.1 x.2.2 ^ (2 * n + 3))⁻¹ *
            dualD x.1 x.2.1 x.2.2 ^ (4 - k)) := by
              rw [hScalar, hDen]
              ring
    _ = _ := by ring

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

theorem rawMomentIntegrand_integrable
    (n A B C k : ℕ) (hk : k ≤ 4) :
    Integrable (rawMomentIntegrand n A B C k) cubeMeasure := by
  rw [integrable_congr (rawMomentIntegrand_ae_eq_bounded n A B C k hk)]
  apply (integrable_const
    (μ := cubeMeasure) (‖(16 * 4 ^ (4 - k) / 2 ^ (2 * n + 3) : ℝ)‖)).mono'
  · exact (boundedMomentIntegrand_measurable n A B C k).aestronglyMeasurable
  · exact boundedMomentIntegrand_norm_le n A B C k hk

def dualMoment (n A B C k : ℕ) : ℝ :=
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

private def dualKernelD (u v : ℝ) : ℝ := u * (1 + v ^ 2) + 2 * v

private def dualKernelBoundary (a : ℕ) (u v : ℝ) : ℝ :=
  2 * u * v ^ (a + 1) * (1 - v ^ 2) /
    dualKernelD u v ^ (a + 2)

private def dualKernelDerivative (a : ℕ) (u v : ℝ) : ℝ :=
  ((a : ℝ) + 1) * (2 * v ^ a / dualKernelD u v ^ (a + 1)) -
    (2 * (a : ℝ) + 3) *
      (4 * v ^ (a + 1) / dualKernelD u v ^ (a + 2)) +
    ((a : ℝ) + 2) * 8 * (1 - u ^ 2) *
      v ^ (a + 2) / dualKernelD u v ^ (a + 3)

private def dualBaseConstant (n : ℕ) (p q : ℝ) : ℝ :=
  16 * p ^ (2 * n + 6) * q ^ (2 * n + 5) *
    (1 - p ^ 2) ^ n * (1 - q ^ 2) ^ (n + 1)

private def dualVBoundary (n : ℕ) (p q v : ℝ) : ℝ :=
  dualBaseConstant n p q * dualKernelBoundary (2 * n + 4) (p * q) v

private def dualVCombination (n : ℕ) (p q v : ℝ) : ℝ :=
  (2 * (n : ℝ) + 5) *
      (2 * rawMomentIntegrand n 0 0 1 1 (p, q, v)) -
    (4 * (n : ℝ) + 11) *
      (4 * rawMomentIntegrand n 0 0 2 2 (p, q, v)) +
    (2 * (n : ℝ) + 6) * 8 *
      (rawMomentIntegrand n 0 0 3 3 (p, q, v) -
        rawMomentIntegrand n 2 2 3 3 (p, q, v))

private theorem dualKernelBoundary_hasDerivAt (a : ℕ) (u v : ℝ)
    (hD : dualKernelD u v ≠ 0) :
    HasDerivAt (dualKernelBoundary a u)
      (dualKernelDerivative a u v) v := by
  have hden : HasDerivAt (dualKernelD u) (2 * u * v + 2) v := by
    unfold dualKernelD
    convert (((hasDerivAt_const v u).mul
      ((hasDerivAt_const v 1).add ((hasDerivAt_id v).pow 2))).add
        ((hasDerivAt_const v 2).mul (hasDerivAt_id v))) using 1 <;>
      simp only [id_eq] <;> ring
  have hnum : HasDerivAt
      (fun z : ℝ => 2 * u * (z ^ (a + 1) * (1 - z ^ 2)))
      (2 * u * (((a + 1 : ℕ) : ℝ) * v ^ a * (1 - v ^ 2) +
        v ^ (a + 1) * (-2 * v))) v := by
    have hpow : HasDerivAt (fun z : ℝ => z ^ (a + 1))
        (((a + 1 : ℕ) : ℝ) * v ^ a) v := by
      convert (hasDerivAt_id v).pow (a + 1) using 1
      simp only [id_eq]
      rw [show a + 1 - 1 = a by omega]
      ring
    have hone : HasDerivAt (fun z : ℝ => 1 - z ^ 2) (-2 * v) v := by
      convert (hasDerivAt_const v 1).sub ((hasDerivAt_id v).pow 2) using 1 <;>
        norm_num [id_eq]
    convert (hasDerivAt_const v (2 * u)).mul (hpow.mul hone) using 1 <;> ring
  have hfun :
      (fun z : ℝ =>
        (2 * u * (z ^ (a + 1) * (1 - z ^ 2))) /
          dualKernelD u z ^ (a + 2)) = dualKernelBoundary a u := by
    funext z
    unfold dualKernelBoundary
    ring
  rw [← hfun]
  convert hnum.div (hden.pow (a + 2)) (pow_ne_zero _ hD) using 1
  unfold dualKernelDerivative
  simp only [Pi.pow_apply]
  simp only [Nat.cast_add, Nat.cast_one]
  rw [show a + 2 - 1 = a + 1 by omega]
  simp only [pow_succ]
  unfold dualKernelD at hD ⊢
  field_simp [hD]
  ring

private theorem rawMomentIntegrand_zero_factor (n C k : ℕ) (p q v : ℝ) :
    rawMomentIntegrand n 0 0 C k (p, q, v) =
      dualBaseConstant n p q *
        (v ^ (2 * n + 3 + C) / dualD p q v ^ (2 * n + 4 + k)) := by
  dsimp [rawMomentIntegrand, dualBaseConstant]
  ring

private theorem rawMomentIntegrand_shift_factor (n C k : ℕ) (p q v : ℝ) :
    rawMomentIntegrand n 2 2 C k (p, q, v) =
      dualBaseConstant n p q * p ^ 2 * q ^ 2 *
        (v ^ (2 * n + 3 + C) / dualD p q v ^ (2 * n + 4 + k)) := by
  dsimp [rawMomentIntegrand, dualBaseConstant]
  rw [show 2 * n + 6 + 2 = (2 * n + 6) + 2 by omega,
    show 2 * n + 5 + 2 = (2 * n + 5) + 2 by omega,
    pow_add, pow_add]
  ring

private theorem dualVCombination_factor (n : ℕ) (p q v : ℝ) :
    dualVCombination n p q v =
      dualBaseConstant n p q *
        dualKernelDerivative (2 * n + 4) (p * q) v := by
  rw [dualVCombination, rawMomentIntegrand_zero_factor,
    rawMomentIntegrand_zero_factor, rawMomentIntegrand_zero_factor,
    rawMomentIntegrand_shift_factor]
  rw [show 2 * n + 3 + 1 = 2 * n + 4 by omega,
    show 2 * n + 4 + 1 = (2 * n + 4) + 1 by omega,
    show 2 * n + 3 + 2 = (2 * n + 4) + 1 by omega,
    show 2 * n + 4 + 2 = (2 * n + 4) + 2 by omega,
    show 2 * n + 3 + 3 = (2 * n + 4) + 2 by omega,
    show 2 * n + 4 + 3 = (2 * n + 4) + 3 by omega]
  rw [show (2 * (n : ℝ) + 5) = ((2 * n + 4 : ℕ) : ℝ) + 1 by
      push_cast; ring,
    show (4 * (n : ℝ) + 11) = 2 * ((2 * n + 4 : ℕ) : ℝ) + 3 by
      push_cast; ring,
    show (2 * (n : ℝ) + 6) = ((2 * n + 4 : ℕ) : ℝ) + 2 by
      push_cast; ring]
  unfold dualKernelDerivative
  rw [show dualKernelD (p * q) v = dualD p q v by rfl]
  ring

private theorem dualVBoundary_hasDerivAt (n : ℕ) (p q v : ℝ)
    (hD : dualD p q v ≠ 0) :
    HasDerivAt (dualVBoundary n p q)
      (dualVCombination n p q v) v := by
  rw [dualVCombination_factor n p q v]
  unfold dualVBoundary
  apply HasDerivAt.const_mul
  apply dualKernelBoundary_hasDerivAt
  simpa only [dualKernelD, dualD] using hD

private theorem dualVCombination_interval_integral_zero (n : ℕ)
    (p q : ℝ) (hp : 0 < p) (hq : 0 < q) :
    ∫ v in (0 : ℝ)..1, dualVCombination n p q v = 0 := by
  have hDpos : ∀ v ∈ Icc (0 : ℝ) 1, 0 < dualD p q v := by
    intro v hv
    dsimp [dualD]
    have hpq : 0 < p * q := mul_pos hp hq
    have hv2 : 0 ≤ v ^ 2 := sq_nonneg v
    have hv0 : 0 ≤ v := hv.1
    positivity
  have hboundary : ContinuousOn (dualVBoundary n p q) (Icc (0 : ℝ) 1) := by
    unfold dualVBoundary dualKernelBoundary dualKernelD
    apply ContinuousOn.mul continuousOn_const
    apply ContinuousOn.div
    · fun_prop
    · fun_prop
    · intro v hv
      exact pow_ne_zero _ (hDpos v hv).ne'
  have hcombination : ContinuousOn (dualVCombination n p q) (Icc (0 : ℝ) 1) := by
    have hraw : ∀ A B C k : ℕ,
        ContinuousOn (fun v : ℝ => rawMomentIntegrand n A B C k (p, q, v))
          (Icc (0 : ℝ) 1) := by
      intro A B C k
      unfold rawMomentIntegrand dualD
      apply ContinuousOn.div
      · fun_prop
      · fun_prop
      · intro v hv
        exact pow_ne_zero _ (hDpos v hv).ne'
    unfold dualVCombination
    simpa [mul_assoc] using ((((hraw 0 0 1 1).const_mul 2).const_mul (2 * (n : ℝ) + 5)).sub
      (((hraw 0 0 2 2).const_mul 4).const_mul (4 * (n : ℝ) + 11))).add
        ((((hraw 0 0 3 3).sub (hraw 2 2 3 3)).const_mul 8).const_mul
          (2 * (n : ℝ) + 6))
  calc
    (∫ v in (0 : ℝ)..1, dualVCombination n p q v) =
        dualVBoundary n p q 1 - dualVBoundary n p q 0 := by
      apply intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le (by norm_num)
        hboundary
      · intro v hv
        exact dualVBoundary_hasDerivAt n p q v
          (hDpos v ⟨hv.1.le, hv.2.le⟩).ne'
      · have hc : ContinuousOn (dualVCombination n p q) (uIcc (0 : ℝ) 1) := by
          simpa [uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hcombination
        exact hc.intervalIntegrable
    _ = 0 := by
      simp [dualVBoundary, dualKernelBoundary]

private theorem dualVCombination_integrable (n : ℕ) :
    Integrable (fun x : ℝ × (ℝ × ℝ) =>
      dualVCombination n x.1 x.2.1 x.2.2) cubeMeasure := by
  have hraw (A B C k : ℕ) (hk : k ≤ 4) :
      Integrable (fun x : ℝ × (ℝ × ℝ) =>
        rawMomentIntegrand n A B C k x) cubeMeasure :=
    rawMomentIntegrand_integrable n A B C k hk
  simpa [dualVCombination, mul_assoc] using
    ((((hraw 0 0 1 1 (by omega)).const_mul 2).const_mul (2 * (n : ℝ) + 5)).sub
      (((hraw 0 0 2 2 (by omega)).const_mul 4).const_mul (4 * (n : ℝ) + 11))).add
        ((((hraw 0 0 3 3 (by omega)).sub (hraw 2 2 3 3 (by omega))).const_mul 8).const_mul
          (2 * (n : ℝ) + 6))

private theorem dualVCombination_cube_integral_zero (n : ℕ) :
    (∫ x : ℝ × (ℝ × ℝ), dualVCombination n x.1 x.2.1 x.2.2 ∂cubeMeasure) = 0 := by
  let f : ℝ × (ℝ × ℝ) → ℝ := fun x =>
    dualVCombination n x.1 x.2.1 x.2.2
  have hf : Integrable f cubeMeasure := dualVCombination_integrable n
  rw [MeasureTheory.integral_prod f hf]
  apply integral_eq_zero_of_ae
  filter_upwards [unit_ae_bounds, hf.prod_right_ae] with p hp hfp
  rw [MeasureTheory.integral_prod _ hfp]
  apply integral_eq_zero_of_ae
  filter_upwards [unit_ae_bounds, hfp.prod_right_ae] with q hq hfq
  have hz := dualVCombination_interval_integral_zero n p q hp.1 hq.1
  rw [intervalIntegral.integral_of_le (by norm_num : (0 : ℝ) ≤ 1)] at hz
  simpa [f, unitMeasure] using hz

def dualVector (n : ℕ) : Fin 3 → ℝ :=
  ![dualMoment n 0 0 0 0,
    2 * ((n : ℝ) + 2) * dualMoment n 0 0 1 1,
    -((n : ℝ) + 2) * dualMoment n 0 0 1 1 +
      2 * ((n : ℝ) + 2) * (2 * (n : ℝ) + 5) *
        dualMoment n 0 0 2 2]
theorem dualVector_zero_pos (n : ℕ) : 0 < dualVector n 0 := by
  simpa [dualVector] using dualMoment_pos n 0 0 0 0 (by omega)

theorem dualVector_one_pos (n : ℕ) : 0 < dualVector n 1 := by
  simp only [dualVector, Matrix.cons_val_one]
  exact mul_pos (by positivity) (dualMoment_pos n 0 0 1 1 (by omega))

theorem dualVector_two_pos (n : ℕ) : 0 < dualVector n 2 := by
  have hzero := dualVCombination_cube_integral_zero n
  have h11 := rawMomentIntegrand_integrable n 0 0 1 1 (by omega)
  have h22 := rawMomentIntegrand_integrable n 0 0 2 2 (by omega)
  have h33 := rawMomentIntegrand_integrable n 0 0 3 3 (by omega)
  have hshift := rawMomentIntegrand_integrable n 2 2 3 3 (by omega)
  have h11' : Integrable (fun x : ℝ × (ℝ × ℝ) =>
      rawMomentIntegrand n 0 0 1 1 (x.1, x.2.1, x.2.2)) cubeMeasure := by
    simpa using h11
  have h22' : Integrable (fun x : ℝ × (ℝ × ℝ) =>
      rawMomentIntegrand n 0 0 2 2 (x.1, x.2.1, x.2.2)) cubeMeasure := by
    simpa using h22
  have h33' : Integrable (fun x : ℝ × (ℝ × ℝ) =>
      rawMomentIntegrand n 0 0 3 3 (x.1, x.2.1, x.2.2)) cubeMeasure := by
    simpa using h33
  have hshift' : Integrable (fun x : ℝ × (ℝ × ℝ) =>
      rawMomentIntegrand n 2 2 3 3 (x.1, x.2.1, x.2.2)) cubeMeasure := by
    simpa using hshift
  have hterm1 := (h11'.const_mul 2).const_mul (2 * (n : ℝ) + 5)
  have hterm2 := (h22'.const_mul 4).const_mul (4 * (n : ℝ) + 11)
  have hleft := hterm1.sub hterm2
  have hright := ((h33'.sub hshift').const_mul 8).const_mul (2 * (n : ℝ) + 6)
  simp only [dualVCombination] at hzero
  simp only [mul_assoc] at hzero
  have hzero1 :
      (∫ x, (2 * (n : ℝ) + 5) *
          (2 * rawMomentIntegrand n 0 0 1 1 (x.1, x.2.1, x.2.2)) -
          (4 * (n : ℝ) + 11) *
            (4 * rawMomentIntegrand n 0 0 2 2 (x.1, x.2.1, x.2.2)) ∂cubeMeasure) +
        (∫ x, (2 * (n : ℝ) + 6) * (8 *
          (rawMomentIntegrand n 0 0 3 3 (x.1, x.2.1, x.2.2) -
            rawMomentIntegrand n 2 2 3 3 (x.1, x.2.1, x.2.2))) ∂cubeMeasure) = 0 := by
    calc
      _ = ∫ x, (((fun y : ℝ × (ℝ × ℝ) =>
            (2 * (n : ℝ) + 5) *
              (2 * rawMomentIntegrand n 0 0 1 1 (y.1, y.2.1, y.2.2))) -
          (fun y : ℝ × (ℝ × ℝ) =>
            (4 * (n : ℝ) + 11) *
              (4 * rawMomentIntegrand n 0 0 2 2 (y.1, y.2.1, y.2.2)))) +
          (fun y : ℝ × (ℝ × ℝ) => (2 * (n : ℝ) + 6) * (8 *
            (rawMomentIntegrand n 0 0 3 3 (y.1, y.2.1, y.2.2) -
              rawMomentIntegrand n 2 2 3 3 (y.1, y.2.1, y.2.2))))) x ∂cubeMeasure :=
        (MeasureTheory.integral_add hleft hright).symm
      _ = 0 := by simpa only [Pi.add_apply, Pi.sub_apply, mul_assoc] using hzero
  rw [MeasureTheory.integral_sub hterm1 hterm2] at hzero1
  simp only [MeasureTheory.integral_const_mul] at hzero1
  rw [MeasureTheory.integral_sub h33' hshift'] at hzero1
  change
    (2 * (n : ℝ) + 5) * (2 * dualMoment n 0 0 1 1) -
        (4 * (n : ℝ) + 11) * (4 * dualMoment n 0 0 2 2) +
      (2 * (n : ℝ) + 6) *
        (8 * (dualMoment n 0 0 3 3 - dualMoment n 2 2 3 3)) = 0 at hzero1
  have hmoment : 0 < dualMoment n 0 0 2 2 :=
    dualMoment_pos n 0 0 2 2 (by omega)
  have hdelta := dualMoment_three_sub_shift_pos n
  have ht : (0 : ℝ) < 2 * (n : ℝ) + 5 := by positivity
  have hn : (0 : ℝ) ≤ n := by positivity
  let d := -((n : ℝ) + 2) * dualMoment n 0 0 1 1 +
    2 * ((n : ℝ) + 2) * (2 * (n : ℝ) + 5) * dualMoment n 0 0 2 2
  have hid :
      (2 * (n : ℝ) + 5) * d = ((n : ℝ) + 2) *
        ((8 * (n : ℝ) ^ 2 + 32 * (n : ℝ) + 28) * dualMoment n 0 0 2 2 +
          (8 * (n : ℝ) + 24) *
            (dualMoment n 0 0 3 3 - dualMoment n 2 2 3 3)) := by
    dsimp [d]
    linear_combination -(((n : ℝ) + 2) / 2) * hzero1
  have hcoef : 0 < 8 * (n : ℝ) ^ 2 + 32 * (n : ℝ) + 28 := by
    nlinarith [sq_nonneg (n : ℝ)]
  have hrhs : 0 < ((n : ℝ) + 2) *
      ((8 * (n : ℝ) ^ 2 + 32 * (n : ℝ) + 28) * dualMoment n 0 0 2 2 +
        (8 * (n : ℝ) + 24) *
          (dualMoment n 0 0 3 3 - dualMoment n 2 2 3 3)) := by
    exact mul_pos (by positivity)
      (add_pos (mul_pos hcoef hmoment) (mul_pos (by positivity) hdelta))
  have hprod : 0 < (2 * (n : ℝ) + 5) * d := hid.symm ▸ hrhs
  have hthird : 0 < d := pos_of_mul_pos_right hprod ht.le
  dsimp [d] at hthird
  simpa [dualVector] using hthird


end RamanujanChallenge.P25

end
