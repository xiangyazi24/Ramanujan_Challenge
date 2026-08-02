import RamanujanChallenge.Problem25DualIntegral
import RamanujanChallenge.Problem25DualCertificateSemantics
import RamanujanChallenge.Problem25DualCertificateRow0

noncomputable section

namespace RamanujanChallenge.P25

open MeasureTheory Set Filter Topology
open scoped Interval

set_option maxRecDepth 10000

def dualCertWeightedIntegrand (n : ℕ) (P : DualCertPoly)
    (x : ℝ × (ℝ × ℝ)) : ℝ :=
  rawMomentIntegrand n 0 0 0 4 x *
    dualCertEval P n x.1 x.2.1 x.2.2

private theorem dualCertWeightedMonomial_integrable (n : ℕ)
    (t : DualCertMonomial) :
    Integrable (fun x : ℝ × (ℝ × ℝ) =>
      rawMomentIntegrand n 0 0 0 4 x *
        dualCertMonomialEval n x.1 x.2.1 x.2.2 t) cubeMeasure := by
  let c : ℝ := (t.coeff : ℝ) * (n : ℝ) ^ t.exp.nExp
  have h := (rawMomentIntegrand_integrable n t.exp.pExp t.exp.qExp
    t.exp.vExp 4 (by omega)).const_mul c
  convert h using 1
  funext x
  dsimp [c, rawMomentIntegrand, dualCertMonomialEval]
  rw [show 2 * n + 6 + t.exp.pExp = (2 * n + 6) + t.exp.pExp by omega,
    show 2 * n + 5 + t.exp.qExp = (2 * n + 5) + t.exp.qExp by omega,
    show 2 * n + 3 + t.exp.vExp = (2 * n + 3) + t.exp.vExp by omega,
    pow_add, pow_add, pow_add]
  ring

theorem dualCertWeightedIntegrand_integrable (n : ℕ) (P : DualCertPoly) :
    Integrable (dualCertWeightedIntegrand n P) cubeMeasure := by
  rcases P with ⟨ts⟩
  unfold dualCertWeightedIntegrand dualCertEval
  induction ts with
  | nil => simp
  | cons t ts ih =>
      simp only [List.map_cons, List.sum_cons, mul_add]
      exact (dualCertWeightedMonomial_integrable n t).add ih

def dualCertFluxP (n : ℕ) (P : DualCertPoly)
    (p q v : ℝ) : ℝ :=
  16 * q ^ (2 * n + 5) * (1 - q ^ 2) ^ (n + 1) *
      v ^ (2 * n + 3) *
      (p ^ (2 * n + 7) * (1 - p ^ 2) ^ (n + 1) *
        dualCertEval P n p q v) /
    dualD p q v ^ (2 * n + 7)

theorem dualCertFluxP_algebra (n : ℕ) (p D F dF Dp C : ℝ)
    (hD : D ≠ 0) :
    (C *
          (((2 * (n : ℝ) + 7) * p ^ (2 * n + 6) * (1 - p ^ 2) ^ (n + 1) +
                p ^ (2 * n + 7) * ((n : ℝ) + 1) * (1 - p ^ 2) ^ n * (-2 * p)) * F +
            p ^ (2 * n + 7) * (1 - p ^ 2) ^ (n + 1) * dF) *
          D ^ (2 * n + 7) -
        C * (p ^ (2 * n + 7) * (1 - p ^ 2) ^ (n + 1) * F) *
          ((2 * (n : ℝ) + 7) * D ^ (2 * n + 6) * Dp)) /
        (D ^ (2 * n + 7)) ^ 2 =
      (C * p ^ (2 * n + 6) * (1 - p ^ 2) ^ n /
          (D ^ (2 * n + 7) * D)) *
        (D * (p * (1 - p ^ 2) * dF +
              ((2 * (n : ℝ) + 7) - (4 * (n : ℝ) + 9) * p ^ 2) * F) -
          (2 * (n : ℝ) + 7) * p * (1 - p ^ 2) * Dp * F) := by
  have hpSucc : p ^ (2 * n + 7) = p ^ (2 * n + 6) * p := by
    rw [show 2 * n + 7 = 2 * n + 6 + 1 by omega, pow_succ]
  have hOneSucc : (1 - p ^ 2) ^ (n + 1) =
      (1 - p ^ 2) ^ n * (1 - p ^ 2) := by
    rw [pow_succ]
  have hDSucc : D ^ (2 * n + 7) = D ^ (2 * n + 6) * D := by
    rw [show 2 * n + 7 = 2 * n + 6 + 1 by omega, pow_succ]
  rw [hpSucc, hOneSucc, hDSucc]
  field_simp [hD]
  ring

theorem dualCertFluxP_hasDerivAt (n : ℕ) (P : DualCertPoly)
    (p q v : ℝ) (hD : dualD p q v ≠ 0) :
    HasDerivAt (fun z => dualCertFluxP n P z q v)
      (dualCertWeightedIntegrand n (dualCertOpPPoly P) (p, q, v)) p := by
  have hpow : HasDerivAt (fun z : ℝ => z ^ (2 * n + 7))
      (((2 * n + 7 : ℕ) : ℝ) * p ^ (2 * n + 6)) p := by
    convert (hasDerivAt_id p).pow (2 * n + 7) using 1
    simp only [id_eq]
    rw [show 2 * n + 7 - 1 = 2 * n + 6 by omega]
    ring
  have hone : HasDerivAt (fun z : ℝ => 1 - z ^ 2) (-2 * p) p := by
    convert (hasDerivAt_const p 1).sub ((hasDerivAt_id p).pow 2) using 1 <;>
      norm_num [id_eq]
  have hD' : HasDerivAt (fun z : ℝ => dualD z q v) (q * (1 + v ^ 2)) p := by
    unfold dualD
    convert (((hasDerivAt_id p).mul_const q).mul_const (1 + v ^ 2)).add_const
      (2 * v) using 1 <;> ring
  have hP := dualCertEval_p_hasDerivAt P n q v p
  let C : ℝ := 16 * q ^ (2 * n + 5) * (1 - q ^ 2) ^ (n + 1) *
    v ^ (2 * n + 3)
  have hnum := ((hpow.mul (hone.pow (n + 1))).mul hP).const_mul C
  have hquot := hnum.div (hD'.pow (2 * n + 7)) (pow_ne_zero _ hD)
  convert hquot using 1
  have hweighted :
      dualCertWeightedIntegrand n (dualCertOpPPoly P) (p, q, v) =
        (C * p ^ (2 * n + 6) * (1 - p ^ 2) ^ n /
          (dualD p q v ^ (2 * n + 7) * dualD p q v)) *
          dualCertOpP (n : ℝ) (dualCertEval P n p q v)
            (dualCertEval (dualCertPDeriv P) n p q v) p q v := by
    unfold dualCertWeightedIntegrand rawMomentIntegrand C
    dsimp only
    rw [dualCertEval_OpPPoly]
    have hDSucc : dualD p q v ^ (2 * n + 8) =
        dualD p q v ^ (2 * n + 7) * dualD p q v := by
      rw [show 2 * n + 8 = 2 * n + 7 + 1 by omega, pow_succ]
    rw [show 2 * n + 4 + 4 = 2 * n + 8 by omega, hDSucc]
    ring
  rw [hweighted]
  simp only [Pi.pow_apply, Pi.mul_apply, Nat.cast_add, Nat.cast_mul, Nat.cast_ofNat]
  rw [show n + 1 - 1 = n by omega,
    show 2 * n + 7 - 1 = 2 * n + 6 by omega]
  unfold dualCertOpP dualCertD
  convert (dualCertFluxP_algebra n p (dualD p q v)
    (dualCertEval P n p q v) (dualCertEval (dualCertPDeriv P) n p q v)
    (q * (1 + v ^ 2)) C hD).symm using 1 <;>
    unfold dualD <;> ring

end RamanujanChallenge.P25

end
