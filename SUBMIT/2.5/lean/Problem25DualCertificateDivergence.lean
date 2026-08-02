import RamanujanChallenge.Problem25DualCertificateIntegral
import RamanujanChallenge.Problem25DualCertificateFluxQV

noncomputable section

namespace RamanujanChallenge.P25

open MeasureTheory Set Filter Topology
open scoped Interval

set_option maxRecDepth 10000

theorem dualCertOpQ_interval_integral_zero (n : ℕ)
    (P : DualCertPoly) (p v : ℝ) (hp : 0 < p) (hv : 0 < v) :
    ∫ q in (0 : ℝ)..1,
      dualCertWeightedIntegrand n (dualCertOpQPoly P) (p, q, v) = 0 := by
  have hDpos : ∀ q ∈ Icc (0 : ℝ) 1, 0 < dualD p q v := by
    intro q hq
    unfold dualD
    exact add_pos_of_nonneg_of_pos
      (mul_nonneg (mul_nonneg hp.le hq.1) (by positivity)) (by positivity)
  have hboundary : ContinuousOn (fun q => dualCertFluxQ n P p q v)
      (Icc (0 : ℝ) 1) := by
    unfold dualCertFluxQ dualD
    apply ContinuousOn.div
    · have hP : Continuous (fun q : ℝ => dualCertEval P n p q v) :=
        continuous_iff_continuousAt.mpr fun q =>
          (dualCertEval_q_hasDerivAt P n p v q).continuousAt
      fun_prop
    · fun_prop
    · intro q hq
      exact pow_ne_zero _ (hDpos q hq).ne'
  have hintegrand : ContinuousOn
      (fun q : ℝ => dualCertWeightedIntegrand n (dualCertOpQPoly P) (p, q, v))
      (Icc (0 : ℝ) 1) := by
    intro q hq
    have hpoly : ContinuousAt
        (fun z : ℝ => dualCertEval (dualCertOpQPoly P) n p z v) q :=
      (dualCertEval_q_hasDerivAt (dualCertOpQPoly P) n p v q).continuousAt
    unfold dualCertWeightedIntegrand rawMomentIntegrand
    dsimp only
    apply ContinuousAt.continuousWithinAt
    apply ContinuousAt.mul
    · apply ContinuousAt.div
      · fun_prop
      · unfold dualD
        fun_prop
      · exact pow_ne_zero _ (hDpos q hq).ne'
    · exact hpoly
  calc
    _ = dualCertFluxQ n P p 1 v - dualCertFluxQ n P p 0 v := by
      apply intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le (by norm_num)
        hboundary
      · intro q hq
        exact dualCertFluxQ_hasDerivAt n P p q v
          (hDpos q ⟨hq.1.le, hq.2.le⟩).ne'
      · have hc : ContinuousOn
            (fun q : ℝ => dualCertWeightedIntegrand n (dualCertOpQPoly P) (p, q, v))
            (uIcc (0 : ℝ) 1) := by
          simpa [uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hintegrand
        exact hc.intervalIntegrable
    _ = 0 := by simp [dualCertFluxQ]

theorem dualCertOpV_interval_integral_zero (n : ℕ)
    (P : DualCertPoly) (p q : ℝ) (hp : 0 < p) (hq : 0 < q) :
    ∫ v in (0 : ℝ)..1,
      dualCertWeightedIntegrand n (dualCertOpVPoly P) (p, q, v) = 0 := by
  have hDpos : ∀ v ∈ Icc (0 : ℝ) 1, 0 < dualD p q v := by
    intro v hv
    unfold dualD
    have hpq : 0 < p * q := mul_pos hp hq
    exact add_pos_of_pos_of_nonneg
      (mul_pos hpq (by positivity)) (mul_nonneg (by norm_num) hv.1)
  have hboundary : ContinuousOn (fun v => dualCertFluxV n P p q v)
      (Icc (0 : ℝ) 1) := by
    unfold dualCertFluxV dualD
    apply ContinuousOn.div
    · have hP : Continuous (fun v : ℝ => dualCertEval P n p q v) :=
        continuous_iff_continuousAt.mpr fun v =>
          (dualCertEval_v_hasDerivAt P n p q v).continuousAt
      fun_prop
    · fun_prop
    · intro v hv
      exact pow_ne_zero _ (hDpos v hv).ne'
  have hintegrand : ContinuousOn
      (fun v : ℝ => dualCertWeightedIntegrand n (dualCertOpVPoly P) (p, q, v))
      (Icc (0 : ℝ) 1) := by
    intro v hv
    have hpoly : ContinuousAt
        (fun z : ℝ => dualCertEval (dualCertOpVPoly P) n p q z) v :=
      (dualCertEval_v_hasDerivAt (dualCertOpVPoly P) n p q v).continuousAt
    unfold dualCertWeightedIntegrand rawMomentIntegrand
    dsimp only
    apply ContinuousAt.continuousWithinAt
    apply ContinuousAt.mul
    · apply ContinuousAt.div
      · fun_prop
      · unfold dualD
        fun_prop
      · exact pow_ne_zero _ (hDpos v hv).ne'
    · exact hpoly
  calc
    _ = dualCertFluxV n P p q 1 - dualCertFluxV n P p q 0 := by
      apply intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le (by norm_num)
        hboundary
      · intro v hv
        exact dualCertFluxV_hasDerivAt n P p q v
          (hDpos v ⟨hv.1.le, hv.2.le⟩).ne'
      · have hc : ContinuousOn
            (fun v : ℝ => dualCertWeightedIntegrand n (dualCertOpVPoly P) (p, q, v))
            (uIcc (0 : ℝ) 1) := by
          simpa [uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hintegrand
        exact hc.intervalIntegrable
    _ = 0 := by simp [dualCertFluxV]

theorem dualCertOpP_cube_integral_zero (n : ℕ) (P : DualCertPoly) :
    (∫ x, dualCertWeightedIntegrand n (dualCertOpPPoly P) x ∂cubeMeasure) = 0 := by
  let f : ℝ × (ℝ × ℝ) → ℝ :=
    dualCertWeightedIntegrand n (dualCertOpPPoly P)
  have hf : Integrable f cubeMeasure := dualCertWeightedIntegrand_integrable _ _
  rw [MeasureTheory.integral_prod_symm f hf]
  apply integral_eq_zero_of_ae
  have hqv : ∀ᵐ y : ℝ × ℝ ∂unitMeasure.prod unitMeasure,
      0 < y.1 ∧ y.1 ≤ 1 ∧ 0 < y.2 ∧ y.2 ≤ 1 := by
    rw [Measure.ae_prod_iff_ae_ae (by measurability : MeasurableSet
      {y : ℝ × ℝ | 0 < y.1 ∧ y.1 ≤ 1 ∧ 0 < y.2 ∧ y.2 ≤ 1})]
    filter_upwards [unit_ae_bounds] with q hq
    filter_upwards [unit_ae_bounds] with v hv
    exact ⟨hq.1, hq.2, hv.1, hv.2⟩
  filter_upwards [hqv] with y hy
  have hz := dualCertOpP_interval_integral_zero n P y.1 y.2 hy.1 hy.2.2.1
  rw [intervalIntegral.integral_of_le (by norm_num : (0 : ℝ) ≤ 1)] at hz
  simpa [f, unitMeasure] using hz

theorem dualCertOpQ_cube_integral_zero (n : ℕ) (P : DualCertPoly) :
    (∫ x, dualCertWeightedIntegrand n (dualCertOpQPoly P) x ∂cubeMeasure) = 0 := by
  let f : ℝ × (ℝ × ℝ) → ℝ :=
    dualCertWeightedIntegrand n (dualCertOpQPoly P)
  have hf : Integrable f cubeMeasure := dualCertWeightedIntegrand_integrable _ _
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
    (∫ x, dualCertWeightedIntegrand n (dualCertOpVPoly P) x ∂cubeMeasure) = 0 := by
  let f : ℝ × (ℝ × ℝ) → ℝ :=
    dualCertWeightedIntegrand n (dualCertOpVPoly P)
  have hf : Integrable f cubeMeasure := dualCertWeightedIntegrand_integrable _ _
  rw [MeasureTheory.integral_prod f hf]
  apply integral_eq_zero_of_ae
  filter_upwards [unit_ae_bounds, hf.prod_right_ae] with p hp hfp
  rw [MeasureTheory.integral_prod _ hfp]
  apply integral_eq_zero_of_ae
  filter_upwards [unit_ae_bounds] with q hq
  have hz := dualCertOpV_interval_integral_zero n P p q hp.1 hq.1
  rw [intervalIntegral.integral_of_le (by norm_num : (0 : ℝ) ≤ 1)] at hz
  simpa [f, unitMeasure] using hz

end RamanujanChallenge.P25

end
