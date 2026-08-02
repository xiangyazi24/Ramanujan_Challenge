import RamanujanChallenge.Problem25DualCertificateRow0

namespace RamanujanChallenge.P25

set_option maxHeartbeats 0 in
set_option maxRecDepth 100000 in
example (n : ℕ) (p q v : ℝ) :
    (4 * (2 * (n : ℝ) + 3) * ((n : ℝ) + 2)) *
        ((positiveMatrix (n : ℤ) 0 0 : ℝ) * dualCertNextNum (n : ℝ) 0 p q v +
          (positiveMatrix (n : ℤ) 0 1 : ℝ) * dualCertNextNum (n : ℝ) 1 p q v +
          (positiveMatrix (n : ℤ) 0 2 : ℝ) * dualCertNextNum (n : ℝ) 2 p q v -
          dualCertLambda (n : ℝ) * dualCertCurNum (n : ℝ) 0 p q v) =
      dualCertOpP (n : ℝ) (dualCertPp0 (n : ℝ) p q v)
          (dualCertDPp0 (n : ℝ) p q v) p q v +
        dualCertOpQ (n : ℝ) (dualCertPq0 (n : ℝ) p q v)
          (dualCertDPq0 (n : ℝ) p q v) p q v +
        dualCertOpV (n : ℝ) (dualCertPv0 (n : ℝ) p q v)
          (dualCertDPv0 (n : ℝ) p q v) p q v := by
  simp [positiveMatrix, Matrix.cons_val_two,
    dualCertD, dualCertSnum, dualCertNextNum, dualCertCurNum, dualCertLambda,
    dualCertOpP, dualCertOpQ, dualCertOpV,
    dualCertPp0, dualCertPq0, dualCertPv0,
    dualCertDPp0, dualCertDPq0, dualCertDPv0,
    dualCertPp0Poly, dualCertPq0Poly, dualCertPv0Poly,
    dualCertPp0Terms, dualCertPq0Terms, dualCertPv0Terms,
    dualCertPDeriv, dualCertQDeriv, dualCertVDeriv,
    dualCertPDerivMonomial, dualCertQDerivMonomial, dualCertVDerivMonomial,
    dualCertOfTerms, dualCertOfTerm, dualCertOfNCoeffs,
    dualCertEval, dualCertMonomialEval]
  <;> ring

end RamanujanChallenge.P25
