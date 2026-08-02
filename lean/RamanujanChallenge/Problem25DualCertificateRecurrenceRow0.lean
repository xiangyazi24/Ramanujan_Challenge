import RamanujanChallenge.Problem25DualCertificateMoments

noncomputable section

namespace RamanujanChallenge.P25

open MeasureTheory Set Filter Topology

set_option maxRecDepth 10000

theorem dualCert_identity_zero_weighted (n : ℕ)
    (x : ℝ × (ℝ × ℝ)) :
    (4 * (2 * (n : ℝ) + 3) * ((n : ℝ) + 2)) *
        ((positiveMatrix (n : ℤ) 0 0 : ℝ) * dualCertWeightedNextIntegrand n 0 x +
          (positiveMatrix (n : ℤ) 0 1 : ℝ) * dualCertWeightedNextIntegrand n 1 x +
          (positiveMatrix (n : ℤ) 0 2 : ℝ) * dualCertWeightedNextIntegrand n 2 x -
          dualCertLambda (n : ℝ) * dualCertWeightedCurIntegrand n 0 x) =
      dualCertWeightedIntegrand n (dualCertOpPPoly dualCertPp0Poly) x +
        dualCertWeightedIntegrand n (dualCertOpQPoly dualCertPq0Poly) x +
        dualCertWeightedIntegrand n (dualCertOpVPoly dualCertPv0Poly) x := by
  have h := congrArg
    (fun z : ℝ => rawMomentIntegrand n 0 0 0 4 x * z)
    (dualCert_identity_zero n x.1 x.2.1 x.2.2)
  unfold dualCertPp0 dualCertPq0 dualCertPv0 at h
  unfold dualCertDPp0 dualCertDPq0 dualCertDPv0 at h
  unfold dualCertWeightedNextIntegrand dualCertWeightedCurIntegrand
  unfold dualCertWeightedIntegrand
  rw [dualCertEval_OpPPoly, dualCertEval_OpQPoly, dualCertEval_OpVPoly]
  convert h using 1 <;> ring

theorem dualVector_recurrence_row_zero (n : ℕ) :
    (positiveMatrix (n : ℤ) 0 0 : ℝ) * dualVector (n + 1) 0 +
        (positiveMatrix (n : ℤ) 0 1 : ℝ) * dualVector (n + 1) 1 +
        (positiveMatrix (n : ℤ) 0 2 : ℝ) * dualVector (n + 1) 2 =
      dualCertLambda (n : ℝ) * dualVector n 0 := by
  let a : ℝ := (positiveMatrix (n : ℤ) 0 0 : ℝ)
  let b : ℝ := (positiveMatrix (n : ℤ) 0 1 : ℝ)
  let c : ℝ := (positiveMatrix (n : ℤ) 0 2 : ℝ)
  let lam : ℝ := dualCertLambda (n : ℝ)
  let delta : ℝ := 4 * (2 * (n : ℝ) + 3) * ((n : ℝ) + 2)
  let A : ℝ × (ℝ × ℝ) → ℝ := fun x =>
    a * dualCertWeightedNextIntegrand n 0 x
  let B : ℝ × (ℝ × ℝ) → ℝ := fun x =>
    b * dualCertWeightedNextIntegrand n 1 x
  let C : ℝ × (ℝ × ℝ) → ℝ := fun x =>
    c * dualCertWeightedNextIntegrand n 2 x
  let E : ℝ × (ℝ × ℝ) → ℝ := fun x =>
    lam * dualCertWeightedCurIntegrand n 0 x
  let P : ℝ × (ℝ × ℝ) → ℝ :=
    dualCertWeightedIntegrand n (dualCertOpPPoly dualCertPp0Poly)
  let Q : ℝ × (ℝ × ℝ) → ℝ :=
    dualCertWeightedIntegrand n (dualCertOpQPoly dualCertPq0Poly)
  let V : ℝ × (ℝ × ℝ) → ℝ :=
    dualCertWeightedIntegrand n (dualCertOpVPoly dualCertPv0Poly)
  have hA : Integrable A cubeMeasure :=
    (dualCertWeightedNextIntegrand_integrable n 0).const_mul a
  have hB : Integrable B cubeMeasure :=
    (dualCertWeightedNextIntegrand_integrable n 1).const_mul b
  have hC : Integrable C cubeMeasure :=
    (dualCertWeightedNextIntegrand_integrable n 2).const_mul c
  have hE : Integrable E cubeMeasure :=
    (dualCertWeightedCurIntegrand_integrable n 0).const_mul lam
  have hP : Integrable P cubeMeasure :=
    dualCertWeightedIntegrand_integrable n (dualCertOpPPoly dualCertPp0Poly)
  have hQ : Integrable Q cubeMeasure :=
    dualCertWeightedIntegrand_integrable n (dualCertOpQPoly dualCertPq0Poly)
  have hV : Integrable V cubeMeasure :=
    dualCertWeightedIntegrand_integrable n (dualCertOpVPoly dualCertPv0Poly)
  have hIntegralEq :
      (∫ x, delta * (A x + B x + C x - E x) ∂cubeMeasure) =
        ∫ x, P x + Q x + V x ∂cubeMeasure := by
    apply MeasureTheory.integral_congr_ae
    exact Eventually.of_forall fun x => by
      simpa [a, b, c, lam, delta, A, B, C, E, P, Q, V] using
        dualCert_identity_zero_weighted n x
  have hRight : (∫ x, P x + Q x + V x ∂cubeMeasure) = 0 := by
    change (∫ x, (P + Q) x + V x ∂cubeMeasure) = 0
    rw [MeasureTheory.integral_add (hP.add hQ) hV]
    change (∫ x, P x + Q x ∂cubeMeasure) + (∫ x, V x ∂cubeMeasure) = 0
    rw [MeasureTheory.integral_add hP hQ]
    rw [show (∫ x, P x ∂cubeMeasure) = 0 by
          simpa [P] using dualCertOpP_cube_integral_zero n dualCertPp0Poly,
      show (∫ x, Q x ∂cubeMeasure) = 0 by
          simpa [Q] using dualCertOpQ_cube_integral_zero n dualCertPq0Poly,
      show (∫ x, V x ∂cubeMeasure) = 0 by
          simpa [V] using dualCertOpV_cube_integral_zero n dualCertPv0Poly]
    ring
  have hLeft :
      (∫ x, delta * (A x + B x + C x - E x) ∂cubeMeasure) =
        delta * (a * dualVector (n + 1) 0 + b * dualVector (n + 1) 1 +
          c * dualVector (n + 1) 2 - lam * dualVector n 0) := by
    rw [MeasureTheory.integral_const_mul]
    change delta * (∫ x, ((A + B) + C) x - E x ∂cubeMeasure) = _
    rw [MeasureTheory.integral_sub ((hA.add hB).add hC) hE]
    change delta * ((∫ x, (A + B) x + C x ∂cubeMeasure) -
      ∫ x, E x ∂cubeMeasure) = _
    rw [MeasureTheory.integral_add (hA.add hB) hC]
    change delta * (((∫ x, A x + B x ∂cubeMeasure) +
      ∫ x, C x ∂cubeMeasure) - ∫ x, E x ∂cubeMeasure) = _
    rw [MeasureTheory.integral_add hA hB]
    simp only [A, B, C, E, MeasureTheory.integral_const_mul]
    rw [dualCertWeightedNextIntegrand_integral n 0,
      dualCertWeightedNextIntegrand_integral n 1,
      dualCertWeightedNextIntegrand_integral n 2,
      dualCertWeightedCurIntegrand_integral n 0]
  have hz : delta * (a * dualVector (n + 1) 0 + b * dualVector (n + 1) 1 +
      c * dualVector (n + 1) 2 - lam * dualVector n 0) = 0 := by
    rw [← hLeft, hIntegralEq, hRight]
  have hdelta : delta ≠ 0 := by
    unfold delta
    positivity
  have hres := (mul_eq_zero.mp hz).resolve_left hdelta
  unfold a b c lam at hres
  linarith

end RamanujanChallenge.P25

end
