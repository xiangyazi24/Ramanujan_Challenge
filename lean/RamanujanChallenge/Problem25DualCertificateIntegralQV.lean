import RamanujanChallenge.Problem25DualCertificateAnalyticQV

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
    (∫ q in (0 : ℝ)..1,
        dualCertWeightedIntegrand n (dualCertOpQPoly P) (p, q, v)) =
        dualCertFluxQ n P p 1 v - dualCertFluxQ n P p 0 v := by
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
    exact add_pos_of_pos_of_nonneg
      (mul_pos (mul_pos hp hq) (by positivity)) (mul_nonneg (by norm_num) hv.1)
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
    (∫ v in (0 : ℝ)..1,
        dualCertWeightedIntegrand n (dualCertOpVPoly P) (p, q, v)) =
        dualCertFluxV n P p q 1 - dualCertFluxV n P p q 0 := by
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

end RamanujanChallenge.P25

end
