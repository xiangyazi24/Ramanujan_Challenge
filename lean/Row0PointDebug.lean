import RamanujanChallenge.Problem25DualCertificateSemantics
import RamanujanChallenge.Problem25DualCertificateRow0

noncomputable section
namespace RamanujanChallenge.P25

macro "row0_eval_num" name:ident value:num : command =>
  `(@[simp] private theorem $name (n p q v : ℝ) :
      dualCertEval ($value : DualCertPoly) n p q v = ($value : ℝ) := by
    change dualCertEval (dualCertConst ($value : ℤ)) n p q v = _
    simp)

row0_eval_num eval136 136
row0_eval_num eval1424 1424
row0_eval_num eval5548 5548
row0_eval_num eval9551 9551
row0_eval_num eval6141 6141
row0_eval_num eval384 384
row0_eval_num eval6384 6384
row0_eval_num eval44168 44168
row0_eval_num eval162698 162698
row0_eval_num eval336377 336377
row0_eval_num eval369933 369933
row0_eval_num eval169011 169011
row0_eval_num eval480 480
row0_eval_num eval4980 4980
row0_eval_num eval19210 19210
row0_eval_num eval32690 32690
row0_eval_num eval20730 20730

example (n : ℕ) (p q v : ℝ) :
    (4 * (2 * (n : ℝ) + 3) * ((n : ℝ) + 2)) *
        ((positiveMatrix (n : ℤ) 0 0 : ℝ) * dualCertNextNum (n : ℝ) 0 p q v +
          (positiveMatrix (n : ℤ) 0 1 : ℝ) * dualCertNextNum (n : ℝ) 1 p q v +
          (positiveMatrix (n : ℤ) 0 2 : ℝ) * dualCertNextNum (n : ℝ) 2 p q v -
          dualCertLambda (n : ℝ) * dualCertCurNum (n : ℝ) 0 p q v) =
      dualCertOpP (n : ℝ)
          (dualCertEval dualCertPp0Poly (n : ℝ) p q v)
          (dualCertEval (dualCertPDeriv dualCertPp0Poly) (n : ℝ) p q v) p q v +
        dualCertOpQ (n : ℝ)
          (dualCertEval dualCertPq0Poly (n : ℝ) p q v)
          (dualCertEval (dualCertQDeriv dualCertPq0Poly) (n : ℝ) p q v) p q v +
        dualCertOpV (n : ℝ)
          (dualCertEval dualCertPv0Poly (n : ℝ) p q v)
          (dualCertEval (dualCertVDeriv dualCertPv0Poly) (n : ℝ) p q v) p q v := by
  have h := congrArg (fun P => dualCertEval P (n : ℝ) p q v)
    dualCert_identity_zero_poly
  dsimp only at h
  rw [dualCertEval_normalize, dualCertEval_normalize] at h
  simp only [dualCertDeltaPoly, dualCertM00Poly, dualCertM01Poly, dualCertM02Poly,
    dualCertNextNumPoly, dualCertCurNumPoly, dualCertLambdaPoly,
    dualCertOpPPoly, dualCertOpQPoly, dualCertOpVPoly,
    dualCertDPoly, dualCertSnumPoly] at h
  simp only [dualCertEval_add, dualCertEval_sub, dualCertEval_neg,
    dualCertEval_mul, dualCertEval_pow, dualCertEval_const,
    dualCertEval_N, dualCertEval_P, dualCertEval_Q, dualCertEval_V] at h
  simp at h
  simpa [dualCertD, dualCertSnum, dualCertNextNum, dualCertCurNum,
    dualCertLambda, dualCertOpP, dualCertOpQ, dualCertOpV,
    positiveMatrix, Matrix.cons_val_two] using h

end RamanujanChallenge.P25
end
