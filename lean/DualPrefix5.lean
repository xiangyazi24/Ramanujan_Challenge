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


end RamanujanChallenge.P25
end
