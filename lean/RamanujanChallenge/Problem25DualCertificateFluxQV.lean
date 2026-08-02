import RamanujanChallenge.Problem25DualCertificateAnalytic

noncomputable section

namespace RamanujanChallenge.P25

open MeasureTheory Set Filter Topology
open scoped Interval

set_option maxRecDepth 10000

def dualCertFluxQ (n : ℕ) (P : DualCertPoly)
    (p q v : ℝ) : ℝ :=
  16 * p ^ (2 * n + 6) * (1 - p ^ 2) ^ n *
      v ^ (2 * n + 3) *
      (q ^ (2 * n + 6) * (1 - q ^ 2) ^ (n + 2) *
        dualCertEval P n p q v) /
    dualD p q v ^ (2 * n + 7)

theorem dualCertFluxQ_algebra (n : ℕ) (q D F dF Dq C : ℝ)
    (hD : D ≠ 0) :
    (C *
          (((2 * (n : ℝ) + 6) * q ^ (2 * n + 5) * (1 - q ^ 2) ^ (n + 2) +
                q ^ (2 * n + 6) * ((n : ℝ) + 2) *
                  (1 - q ^ 2) ^ (n + 1) * (-2 * q)) * F +
            q ^ (2 * n + 6) * (1 - q ^ 2) ^ (n + 2) * dF) *
          D ^ (2 * n + 7) -
        C * (q ^ (2 * n + 6) * (1 - q ^ 2) ^ (n + 2) * F) *
          ((2 * (n : ℝ) + 7) * D ^ (2 * n + 6) * Dq)) /
        (D ^ (2 * n + 7)) ^ 2 =
      (C * q ^ (2 * n + 5) * (1 - q ^ 2) ^ (n + 1) /
          (D ^ (2 * n + 7) * D)) *
        (D * (q * (1 - q ^ 2) * dF +
              ((2 * (n : ℝ) + 6) - (4 * (n : ℝ) + 10) * q ^ 2) * F) -
          (2 * (n : ℝ) + 7) * q * (1 - q ^ 2) * Dq * F) := by
  have hqSucc : q ^ (2 * n + 6) = q ^ (2 * n + 5) * q := by
    rw [show 2 * n + 6 = 2 * n + 5 + 1 by omega, pow_succ]
  have hOneSucc : (1 - q ^ 2) ^ (n + 2) =
      (1 - q ^ 2) ^ (n + 1) * (1 - q ^ 2) := by
    rw [show n + 2 = (n + 1) + 1 by omega, pow_succ]
  have hDSucc : D ^ (2 * n + 7) = D ^ (2 * n + 6) * D := by
    rw [show 2 * n + 7 = 2 * n + 6 + 1 by omega, pow_succ]
  rw [hqSucc, hOneSucc, hDSucc]
  field_simp [hD]
  ring

theorem dualCertFluxQ_hasDerivAt (n : ℕ) (P : DualCertPoly)
    (p q v : ℝ) (hD : dualD p q v ≠ 0) :
    HasDerivAt (fun z => dualCertFluxQ n P p z v)
      (dualCertWeightedIntegrand n (dualCertOpQPoly P) (p, q, v)) q := by
  have hpow : HasDerivAt (fun z : ℝ => z ^ (2 * n + 6))
      (((2 * n + 6 : ℕ) : ℝ) * q ^ (2 * n + 5)) q := by
    convert (hasDerivAt_id q).pow (2 * n + 6) using 1
    simp only [id_eq]
    rw [show 2 * n + 6 - 1 = 2 * n + 5 by omega]
    ring
  have hone : HasDerivAt (fun z : ℝ => 1 - z ^ 2) (-2 * q) q := by
    convert (hasDerivAt_const q 1).sub ((hasDerivAt_id q).pow 2) using 1 <;>
      norm_num [id_eq]
  have hD' : HasDerivAt (fun z : ℝ => dualD p z v) (p * (1 + v ^ 2)) q := by
    unfold dualD
    convert (((hasDerivAt_id q).const_mul p).mul_const
      (1 + v ^ 2)).add_const (2 * v) using 1 <;> ring
  have hP := dualCertEval_q_hasDerivAt P n p v q
  let C : ℝ := 16 * p ^ (2 * n + 6) * (1 - p ^ 2) ^ n *
    v ^ (2 * n + 3)
  have hnum := ((hpow.mul (hone.pow (n + 2))).mul hP).const_mul C
  have hquot := hnum.div (hD'.pow (2 * n + 7)) (pow_ne_zero _ hD)
  convert hquot using 1
  have hweighted :
      dualCertWeightedIntegrand n (dualCertOpQPoly P) (p, q, v) =
        (C * q ^ (2 * n + 5) * (1 - q ^ 2) ^ (n + 1) /
          (dualD p q v ^ (2 * n + 7) * dualD p q v)) *
          dualCertOpQ (n : ℝ) (dualCertEval P n p q v)
            (dualCertEval (dualCertQDeriv P) n p q v) p q v := by
    unfold dualCertWeightedIntegrand rawMomentIntegrand C
    dsimp only
    rw [dualCertEval_OpQPoly]
    have hDSucc : dualD p q v ^ (2 * n + 8) =
        dualD p q v ^ (2 * n + 7) * dualD p q v := by
      rw [show 2 * n + 8 = 2 * n + 7 + 1 by omega, pow_succ]
    rw [show 2 * n + 4 + 4 = 2 * n + 8 by omega, hDSucc]
    ring
  rw [hweighted]
  simp only [Pi.pow_apply, Pi.mul_apply, Nat.cast_add, Nat.cast_mul, Nat.cast_ofNat]
  rw [show n + 2 - 1 = n + 1 by omega,
    show 2 * n + 7 - 1 = 2 * n + 6 by omega]
  unfold dualCertOpQ dualCertD
  convert (dualCertFluxQ_algebra n q (dualD p q v)
    (dualCertEval P n p q v) (dualCertEval (dualCertQDeriv P) n p q v)
    (p * (1 + v ^ 2)) C hD).symm using 1 <;>
    unfold dualD <;> ring

def dualCertFluxV (n : ℕ) (P : DualCertPoly)
    (p q v : ℝ) : ℝ :=
  16 * p ^ (2 * n + 6) * q ^ (2 * n + 5) *
      (1 - p ^ 2) ^ n * (1 - q ^ 2) ^ (n + 1) *
      (v ^ (2 * n + 4) * (1 - v ^ 2) *
        dualCertEval P n p q v) /
    dualD p q v ^ (2 * n + 7)

theorem dualCertFluxV_algebra (n : ℕ) (v D F dF Dv C : ℝ)
    (hD : D ≠ 0) :
    (C *
          (((2 * (n : ℝ) + 4) * v ^ (2 * n + 3) * (1 - v ^ 2) +
                v ^ (2 * n + 4) * (-2 * v)) * F +
            v ^ (2 * n + 4) * (1 - v ^ 2) * dF) *
          D ^ (2 * n + 7) -
        C * (v ^ (2 * n + 4) * (1 - v ^ 2) * F) *
          ((2 * (n : ℝ) + 7) * D ^ (2 * n + 6) * Dv)) /
        (D ^ (2 * n + 7)) ^ 2 =
      (C * v ^ (2 * n + 3) /
          (D ^ (2 * n + 7) * D)) *
        (D * (v * (1 - v ^ 2) * dF +
              ((2 * (n : ℝ) + 4) - (2 * (n : ℝ) + 6) * v ^ 2) * F) -
          (2 * (n : ℝ) + 7) * v * (1 - v ^ 2) * Dv * F) := by
  have hvSucc : v ^ (2 * n + 4) = v ^ (2 * n + 3) * v := by
    rw [show 2 * n + 4 = 2 * n + 3 + 1 by omega, pow_succ]
  have hDSucc : D ^ (2 * n + 7) = D ^ (2 * n + 6) * D := by
    rw [show 2 * n + 7 = 2 * n + 6 + 1 by omega, pow_succ]
  rw [hvSucc, hDSucc]
  field_simp [hD]
  ring

theorem dualCertFluxV_hasDerivAt (n : ℕ) (P : DualCertPoly)
    (p q v : ℝ) (hD : dualD p q v ≠ 0) :
    HasDerivAt (fun z => dualCertFluxV n P p q z)
      (dualCertWeightedIntegrand n (dualCertOpVPoly P) (p, q, v)) v := by
  have hpow : HasDerivAt (fun z : ℝ => z ^ (2 * n + 4))
      (((2 * n + 4 : ℕ) : ℝ) * v ^ (2 * n + 3)) v := by
    convert (hasDerivAt_id v).pow (2 * n + 4) using 1
    simp only [id_eq]
    rw [show 2 * n + 4 - 1 = 2 * n + 3 by omega]
    ring
  have hone : HasDerivAt (fun z : ℝ => 1 - z ^ 2) (-2 * v) v := by
    convert (hasDerivAt_const v 1).sub ((hasDerivAt_id v).pow 2) using 1 <;>
      norm_num [id_eq]
  have hD' : HasDerivAt (fun z : ℝ => dualD p q z) (2 * p * q * v + 2) v := by
    unfold dualD
    convert ((hasDerivAt_const v (p * q)).mul
      ((hasDerivAt_const v 1).add ((hasDerivAt_id v).pow 2))).add
        ((hasDerivAt_const v 2).mul (hasDerivAt_id v)) using 1 <;>
      simp only [id_eq] <;> ring
  have hP := dualCertEval_v_hasDerivAt P n p q v
  let C : ℝ := 16 * p ^ (2 * n + 6) * q ^ (2 * n + 5) *
    (1 - p ^ 2) ^ n * (1 - q ^ 2) ^ (n + 1)
  have hnum := ((hpow.mul hone).mul hP).const_mul C
  have hquot := hnum.div (hD'.pow (2 * n + 7)) (pow_ne_zero _ hD)
  convert hquot using 1
  have hweighted :
      dualCertWeightedIntegrand n (dualCertOpVPoly P) (p, q, v) =
        (C * v ^ (2 * n + 3) /
          (dualD p q v ^ (2 * n + 7) * dualD p q v)) *
          dualCertOpV (n : ℝ) (dualCertEval P n p q v)
            (dualCertEval (dualCertVDeriv P) n p q v) p q v := by
    unfold dualCertWeightedIntegrand rawMomentIntegrand C
    dsimp only
    rw [dualCertEval_OpVPoly]
    have hDSucc : dualD p q v ^ (2 * n + 8) =
        dualD p q v ^ (2 * n + 7) * dualD p q v := by
      rw [show 2 * n + 8 = 2 * n + 7 + 1 by omega, pow_succ]
    rw [show 2 * n + 4 + 4 = 2 * n + 8 by omega, hDSucc]
  rw [hweighted]
  simp only [Pi.pow_apply, Pi.mul_apply, Nat.cast_add, Nat.cast_mul, Nat.cast_ofNat]
  rw [show 2 * n + 7 - 1 = 2 * n + 6 by omega]
  unfold dualCertOpV dualCertD
  convert (dualCertFluxV_algebra n v (dualD p q v)
    (dualCertEval P n p q v) (dualCertEval (dualCertVDeriv P) n p q v)
    (2 * p * q * v + 2) C hD).symm using 1 <;>
    unfold dualD <;> ring

end RamanujanChallenge.P25

end
