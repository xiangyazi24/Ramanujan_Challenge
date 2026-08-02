import RamanujanChallenge.Problem25DualCertificateAnalytic

noncomputable section

namespace RamanujanChallenge.P25

open MeasureTheory Set Filter Topology
open scoped Interval

set_option maxRecDepth 10000

theorem dualCertOpP_interval_integral_zero (n : ℕ)
    (P : DualCertPoly) (q v : ℝ) (hq : 0 < q) (hv : 0 < v) :
    ∫ p in (0 : ℝ)..1,
      dualCertWeightedIntegrand n (dualCertOpPPoly P) (p, q, v) = 0 := by
  have hDpos : ∀ p ∈ Icc (0 : ℝ) 1, 0 < dualD p q v := by
    intro p hp
    unfold dualD
    exact add_pos_of_nonneg_of_pos
      (mul_nonneg (mul_nonneg hp.1 hq.le) (by positivity)) (by positivity)
  have hboundary : ContinuousOn (fun p => dualCertFluxP n P p q v)
      (Icc (0 : ℝ) 1) := by
    unfold dualCertFluxP dualD
    apply ContinuousOn.div
    · have hP : Continuous (fun p : ℝ => dualCertEval P n p q v) :=
        continuous_iff_continuousAt.mpr fun p =>
          (dualCertEval_p_hasDerivAt P n q v p).continuousAt
      fun_prop
    · fun_prop
    · intro p hp
      exact pow_ne_zero _ (hDpos p hp).ne'
  have hintegrand : ContinuousOn
      (fun p : ℝ => dualCertWeightedIntegrand n (dualCertOpPPoly P) (p, q, v))
      (Icc (0 : ℝ) 1) := by
    intro p hp
    have hpoly : ContinuousAt
        (fun z : ℝ => dualCertEval (dualCertOpPPoly P) n z q v) p :=
      (dualCertEval_p_hasDerivAt (dualCertOpPPoly P) n q v p).continuousAt
    unfold dualCertWeightedIntegrand rawMomentIntegrand
    dsimp only
    apply ContinuousAt.continuousWithinAt
    apply ContinuousAt.mul
    · apply ContinuousAt.div
      · fun_prop
      · unfold dualD
        fun_prop
      · exact pow_ne_zero _ (hDpos p hp).ne'
    · exact hpoly
  calc
    (∫ p in (0 : ℝ)..1,
        dualCertWeightedIntegrand n (dualCertOpPPoly P) (p, q, v)) =
        dualCertFluxP n P 1 q v - dualCertFluxP n P 0 q v := by
      apply intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le (by norm_num)
        hboundary
      · intro p hp
        exact dualCertFluxP_hasDerivAt n P p q v
          (hDpos p ⟨hp.1.le, hp.2.le⟩).ne'
      · have hc : ContinuousOn
            (fun p : ℝ => dualCertWeightedIntegrand n (dualCertOpPPoly P) (p, q, v))
            (uIcc (0 : ℝ) 1) := by
          simpa [uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hintegrand
        exact hc.intervalIntegrable
    _ = 0 := by simp [dualCertFluxP]

end RamanujanChallenge.P25

end
