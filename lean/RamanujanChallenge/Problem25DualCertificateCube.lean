import RamanujanChallenge.Problem25DualCertificateIntegral
import RamanujanChallenge.Problem25DualCertificateIntegralQV

noncomputable section

namespace RamanujanChallenge.P25

open MeasureTheory Set Filter Topology
open scoped Interval

set_option maxRecDepth 10000

theorem dualCertOpP_cube_integral_zero (n : ℕ) (P : DualCertPoly) :
    ∫ x : ℝ × (ℝ × ℝ),
      dualCertWeightedIntegrand n (dualCertOpPPoly P) x ∂cubeMeasure = 0 := by
  let f : ℝ × (ℝ × ℝ) → ℝ := fun x =>
    dualCertWeightedIntegrand n (dualCertOpPPoly P) x
  have hf : Integrable f cubeMeasure :=
    dualCertWeightedIntegrand_integrable n (dualCertOpPPoly P)
  rw [MeasureTheory.integral_prod_symm f hf]
  apply integral_eq_zero_of_ae
  have hqv : ∀ᵐ y : ℝ × ℝ ∂unitMeasure.prod unitMeasure,
      0 < y.1 ∧ y.1 ≤ 1 ∧ 0 < y.2 ∧ y.2 ≤ 1 := by
    have hmem : ∀ᵐ y : ℝ × ℝ ∂unitMeasure.prod unitMeasure,
        y ∈ Ioc (0 : ℝ) 1 ×ˢ Ioc (0 : ℝ) 1 := by
      rw [Measure.ae_prod_mem_iff_ae_ae_mem
        (measurableSet_Ioc.prod measurableSet_Ioc)]
      filter_upwards [unit_ae_bounds] with q hq
      filter_upwards [unit_ae_bounds] with v hv
      exact ⟨hq, hv⟩
    filter_upwards [hmem] with y hy
    exact ⟨hy.1.1, hy.1.2, hy.2.1, hy.2.2⟩
  filter_upwards [hqv] with y hy
  have hz := dualCertOpP_interval_integral_zero n P y.1 y.2 hy.1 hy.2.2.1
  rw [intervalIntegral.integral_of_le (by norm_num : (0 : ℝ) ≤ 1)] at hz
  simpa [f, unitMeasure] using hz

theorem dualCertOpQ_cube_integral_zero (n : ℕ) (P : DualCertPoly) :
    ∫ x : ℝ × (ℝ × ℝ),
      dualCertWeightedIntegrand n (dualCertOpQPoly P) x ∂cubeMeasure = 0 := by
  let f : ℝ × (ℝ × ℝ) → ℝ := fun x =>
    dualCertWeightedIntegrand n (dualCertOpQPoly P) x
  have hf : Integrable f cubeMeasure :=
    dualCertWeightedIntegrand_integrable n (dualCertOpQPoly P)
  rw [MeasureTheory.integral_prod f hf]
  apply integral_eq_zero_of_ae
  filter_upwards [unit_ae_bounds, hf.prod_right_ae] with p hp hfp
  rw [MeasureTheory.integral_prod_symm _ hfp]
  apply integral_eq_zero_of_ae
  filter_upwards [unit_ae_bounds] with v hv
  have hz := dualCertOpQ_interval_integral_zero n P p v hp.1 hv.1
  rw [intervalIntegral.integral_of_le (by norm_num : (0 : ℝ) ≤ 1)] at hz
  simpa [f, unitMeasure] using hz

theorem dualCertOpV_cube_integral_zero (n : ℕ) (P : DualCertPoly) :
    ∫ x : ℝ × (ℝ × ℝ),
      dualCertWeightedIntegrand n (dualCertOpVPoly P) x ∂cubeMeasure = 0 := by
  let f : ℝ × (ℝ × ℝ) → ℝ := fun x =>
    dualCertWeightedIntegrand n (dualCertOpVPoly P) x
  have hf : Integrable f cubeMeasure :=
    dualCertWeightedIntegrand_integrable n (dualCertOpVPoly P)
  rw [MeasureTheory.integral_prod f hf]
  apply integral_eq_zero_of_ae
  filter_upwards [unit_ae_bounds, hf.prod_right_ae] with p hp hfp
  rw [MeasureTheory.integral_prod _ hfp]
  apply integral_eq_zero_of_ae
  filter_upwards [unit_ae_bounds, hfp.prod_right_ae] with q hq hfq
  have hz := dualCertOpV_interval_integral_zero n P p q hp.1 hq.1
  rw [intervalIntegral.integral_of_le (by norm_num : (0 : ℝ) ≤ 1)] at hz
  simpa [f, unitMeasure] using hz

end RamanujanChallenge.P25

end
