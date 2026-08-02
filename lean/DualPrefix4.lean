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
  simp only [div_pow, div_eq_mul_inv, mul_pow]
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


end RamanujanChallenge.P25
end
