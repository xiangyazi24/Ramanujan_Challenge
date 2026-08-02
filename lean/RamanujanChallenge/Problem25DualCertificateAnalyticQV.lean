import RamanujanChallenge.Problem25DualCertificateFluxAlgebra

noncomputable section

namespace RamanujanChallenge.P25

open MeasureTheory Set Filter Topology
open scoped Interval

set_option maxRecDepth 10000
set_option maxHeartbeats 2000000

def dualCertFluxQ (n : ℕ) (P : DualCertPoly)
    (p q v : ℝ) : ℝ :=
  16 * p ^ (2 * n + 6) * (1 - p ^ 2) ^ n * v ^ (2 * n + 3) *
      (q ^ (2 * n + 6) * (1 - q ^ 2) ^ (n + 2) *
        dualCertEval P n p q v) /
    dualD p q v ^ (2 * n + 7)

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
    convert (((hasDerivAt_id q).const_mul p).mul_const (1 + v ^ 2)).add_const
      (2 * v) using 1 <;> ring
  have hP := dualCertEval_q_hasDerivAt P n p v q
  let C : ℝ := 16 * p ^ (2 * n + 6) * (1 - p ^ 2) ^ n * v ^ (2 * n + 3)
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
  unfold dualCertOpQ
  have ha : 2 * n + 5 + 1 = 2 * n + 6 := by omega
  have hb : n + 1 + 1 = n + 2 := by omega
  have hr : 2 * n + 6 + 1 = 2 * n + 7 := by omega
  have hA : (((2 * n + 5 + 1 : ℕ) : ℝ)) = 2 * (n : ℝ) + 6 := by
    push_cast
    ring
  have hB : (((n + 1 + 1 : ℕ) : ℝ)) = (n : ℝ) + 2 := by
    push_cast
    ring
  have hAB : (((2 * n + 5 + 1 : ℕ) : ℝ)) +
      2 * (((n + 1 + 1 : ℕ) : ℝ)) = 4 * (n : ℝ) + 10 := by
    push_cast
    ring
  have hAB' : 2 * (n : ℝ) + 6 + 2 * ((n : ℝ) + 2) =
      4 * (n : ℝ) + 10 := by ring
  have hR : (((2 * n + 6 + 1 : ℕ) : ℝ)) = 2 * (n : ℝ) + 7 := by
    push_cast
    ring
  have hAssoc : q ^ (2 * n + 6) * ((n : ℝ) + 2) *
      (1 - q ^ 2) ^ (n + 1) * (-2 * q) =
      q ^ (2 * n + 6) *
        (((n : ℝ) + 2) * (1 - q ^ 2) ^ (n + 1) * (-2 * q)) := by
    ring
  simpa only [ha, hb, hr, hA, hB, hAB, hAB', hR, dualCertD, dualD,
    hAssoc] using
    (dualCertFlux_algebra (2 * n + 5) (n + 1) (2 * n + 6)
    q (dualD p q v) (dualCertEval P n p q v)
    (dualCertEval (dualCertQDeriv P) n p q v)
    (p * (1 + v ^ 2)) C hD).symm

def dualCertFluxV (n : ℕ) (P : DualCertPoly)
    (p q v : ℝ) : ℝ :=
  16 * p ^ (2 * n + 6) * q ^ (2 * n + 5) *
      (1 - p ^ 2) ^ n * (1 - q ^ 2) ^ (n + 1) *
      (v ^ (2 * n + 4) * (1 - v ^ 2) ^ 1 * dualCertEval P n p q v) /
    dualD p q v ^ (2 * n + 7)

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
  have honePlus : HasDerivAt (fun z : ℝ => 1 + z ^ 2) (2 * v) v := by
    convert (hasDerivAt_const v 1).add ((hasDerivAt_id v).pow 2) using 1 <;>
      norm_num [id_eq]
  have hD' : HasDerivAt (fun z : ℝ => dualD p q z) (2 * p * q * v + 2) v := by
    unfold dualD
    convert (honePlus.const_mul (p * q)).add ((hasDerivAt_id v).const_mul 2) using 1 <;>
      ring
  have hP := dualCertEval_v_hasDerivAt P n p q v
  let C : ℝ := 16 * p ^ (2 * n + 6) * q ^ (2 * n + 5) *
    (1 - p ^ 2) ^ n * (1 - q ^ 2) ^ (n + 1)
  have hnum := ((hpow.mul (hone.pow 1)).mul hP).const_mul C
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
  rw [show 1 - 1 = 0 by omega,
    show 2 * n + 7 - 1 = 2 * n + 6 by omega]
  simp only [pow_zero, Nat.cast_one, one_mul]
  unfold dualCertOpV
  have ha : 2 * n + 3 + 1 = 2 * n + 4 := by omega
  have hb : 0 + 1 = 1 := by omega
  have hr : 2 * n + 6 + 1 = 2 * n + 7 := by omega
  have hA : (((2 * n + 3 + 1 : ℕ) : ℝ)) = 2 * (n : ℝ) + 4 := by
    push_cast
    ring
  have hAB : (((2 * n + 3 + 1 : ℕ) : ℝ)) +
      2 * (((0 + 1 : ℕ) : ℝ)) = 2 * (n : ℝ) + 6 := by
    push_cast
    ring
  have hAB' : 2 * (n : ℝ) + 4 + 2 = 2 * (n : ℝ) + 6 := by ring
  have hR : (((2 * n + 6 + 1 : ℕ) : ℝ)) = 2 * (n : ℝ) + 7 := by
    push_cast
    ring
  simpa only [ha, hb, hr, hA, hAB, hAB', hR, Nat.cast_one, pow_zero, pow_one,
    one_mul, mul_one, dualCertD, dualD] using
    (dualCertFlux_algebra (2 * n + 3) 0 (2 * n + 6)
    v (dualD p q v) (dualCertEval P n p q v)
    (dualCertEval (dualCertVDeriv P) n p q v)
    (2 * p * q * v + 2) C hD).symm

end RamanujanChallenge.P25

end
