import RamanujanChallenge.Problem25DualCertificateCube
import RamanujanChallenge.Problem25DualCertificateRow0

noncomputable section

namespace RamanujanChallenge.P25

open MeasureTheory Set Filter Topology

set_option maxRecDepth 10000
set_option maxHeartbeats 1000000

def dualCertWeightedNextIntegrand (n : ℕ) (j : Fin 3)
    (x : ℝ × (ℝ × ℝ)) : ℝ :=
  rawMomentIntegrand n 0 0 0 4 x *
    dualCertNextNum (n : ℝ) j x.1 x.2.1 x.2.2

def dualCertWeightedCurIntegrand (n : ℕ) (i : Fin 3)
    (x : ℝ × (ℝ × ℝ)) : ℝ :=
  rawMomentIntegrand n 0 0 0 4 x *
    dualCertCurNum (n : ℝ) i x.1 x.2.1 x.2.2

private def dualCertBaseConstant (n : ℕ) (p q : ℝ) : ℝ :=
  16 * p ^ (2 * n + 6) * q ^ (2 * n + 5) *
    (1 - p ^ 2) ^ n * (1 - q ^ 2) ^ (n + 1)

private theorem rawMomentIntegrand_zero_factor (n C k : ℕ) (p q v : ℝ) :
    rawMomentIntegrand n 0 0 C k (p, q, v) =
      dualCertBaseConstant n p q *
        (v ^ (2 * n + 3 + C) / dualD p q v ^ (2 * n + 4 + k)) := by
  dsimp [rawMomentIntegrand, dualCertBaseConstant]
  ring

private theorem dualCertBaseConstant_succ (n : ℕ) (p q : ℝ) :
    dualCertBaseConstant n p q * p ^ 2 * q ^ 2 *
        (1 - p ^ 2) * (1 - q ^ 2) =
      dualCertBaseConstant (n + 1) p q := by
  unfold dualCertBaseConstant
  rw [show 2 * (n + 1) + 6 = (2 * n + 6) + 2 by omega,
    show 2 * (n + 1) + 5 = (2 * n + 5) + 2 by omega,
    show n + 1 + 1 = (n + 1) + 1 by omega,
    pow_add, pow_add, pow_succ, pow_succ]
  ring

private theorem div_pow_add_mul_pow_cancel (A Z D : ℝ) (m k : ℕ)
    (hm : 0 < m) :
    A / D ^ (m + k) * (Z * D ^ k) = A * Z / D ^ m := by
  by_cases hD : D = 0
  · simp [hD, ne_of_gt hm]
  · rw [pow_add]
    field_simp [hD]

private theorem mul_div_pow_add_mul_pow_cancel (B A Z D : ℝ) (m k : ℕ)
    (hm : 0 < m) :
    B * (A / D ^ (m + k)) * (Z * D ^ k) =
      B * (A * Z / D ^ m) := by
  rw [mul_assoc, div_pow_add_mul_pow_cancel A Z D m k hm]

private theorem dualCertD_eq_dualD (p q v : ℝ) :
    dualCertD p q v = dualD p q v := by rfl

theorem dualCertWeightedNextIntegrand_eq (n : ℕ) (j : Fin 3)
    (x : ℝ × (ℝ × ℝ)) :
    dualCertWeightedNextIntegrand n j x =
      ![rawMomentIntegrand (n + 1) 0 0 0 0 x,
        2 * ((n : ℝ) + 3) * rawMomentIntegrand (n + 1) 0 0 1 1 x,
        -((n : ℝ) + 3) * rawMomentIntegrand (n + 1) 0 0 1 1 x +
          2 * ((n : ℝ) + 3) * (2 * (n : ℝ) + 7) *
            rawMomentIntegrand (n + 1) 0 0 2 2 x] j := by
  rcases x with ⟨p, q, v⟩
  fin_cases j
  · simp [dualCertWeightedNextIntegrand, dualCertNextNum]
    rw [rawMomentIntegrand_zero_factor n 0 4,
      rawMomentIntegrand_zero_factor (n + 1) 0 0,
      dualCertD_eq_dualD]
    unfold dualCertSnum
    simp only [Nat.add_zero]
    rw [show 2 * n + 4 + 4 = (2 * n + 6) + 2 by omega,
      show 2 * (n + 1) + 3 = 2 * n + 5 by omega,
      show 2 * (n + 1) + 4 = 2 * n + 6 by omega]
    rw [mul_div_pow_add_mul_pow_cancel
      (dualCertBaseConstant n p q) (v ^ (2 * n + 3))
      (p ^ 2 * q ^ 2 * (1 - p ^ 2) * (1 - q ^ 2) * v ^ 2)
      (dualD p q v) (2 * n + 6) 2 (by omega)]
    rw [← dualCertBaseConstant_succ n p q,
      show v ^ (2 * n + 5) = v ^ (2 * n + 3) * v ^ 2 by
        rw [show 2 * n + 5 = (2 * n + 3) + 2 by omega, pow_add]]
    ring
  · simp [dualCertWeightedNextIntegrand, dualCertNextNum]
    rw [rawMomentIntegrand_zero_factor n 0 4,
      rawMomentIntegrand_zero_factor (n + 1) 1 1,
      dualCertD_eq_dualD]
    unfold dualCertSnum
    simp only [Nat.add_zero]
    rw [show 2 * n + 4 + 4 = (2 * n + 7) + 1 by omega,
      show 2 * (n + 1) + 3 + 1 = 2 * n + 6 by omega,
      show 2 * n + 6 + 1 = 2 * n + 7 by omega]
    have hc := mul_div_pow_add_mul_pow_cancel
      (dualCertBaseConstant n p q) (v ^ (2 * n + 3))
      (p ^ 2 * q ^ 2 * (1 - p ^ 2) * (1 - q ^ 2) * v ^ 2 *
        (2 * ((n : ℝ) + 3)) * v)
      (dualD p q v) (2 * n + 7) 1 (by omega)
    simp only [pow_one] at hc
    calc
      _ = dualCertBaseConstant n p q *
            (v ^ (2 * n + 3) / dualD p q v ^ (2 * n + 7 + 1)) *
          ((p ^ 2 * q ^ 2 * (1 - p ^ 2) * (1 - q ^ 2) * v ^ 2 *
            (2 * ((n : ℝ) + 3)) * v) * dualD p q v) := by ring
      _ = dualCertBaseConstant n p q *
          (v ^ (2 * n + 3) *
            (p ^ 2 * q ^ 2 * (1 - p ^ 2) * (1 - q ^ 2) * v ^ 2 *
              (2 * ((n : ℝ) + 3)) * v) / dualD p q v ^ (2 * n + 7)) := hc
      _ = _ := by
        rw [← dualCertBaseConstant_succ n p q,
          show v ^ (2 * n + 6) = v ^ (2 * n + 3) * v ^ 3 by
            rw [show 2 * n + 6 = (2 * n + 3) + 3 by omega, pow_add]]
        ring
  · simp [dualCertWeightedNextIntegrand, dualCertNextNum]
    rw [rawMomentIntegrand_zero_factor n 0 4,
      rawMomentIntegrand_zero_factor (n + 1) 1 1,
      rawMomentIntegrand_zero_factor (n + 1) 2 2,
      dualCertD_eq_dualD]
    unfold dualCertSnum
    simp only [Nat.add_zero]
    rw [show 2 * n + 4 + 4 = (2 * n + 7) + 1 by omega,
      show 2 * (n + 1) + 3 + 1 = 2 * n + 6 by omega,
      show 2 * n + 6 + 1 = 2 * n + 7 by omega,
      show 2 * n + 6 + 2 = 2 * n + 8 by omega]
    have h1 := mul_div_pow_add_mul_pow_cancel
      (dualCertBaseConstant n p q) (v ^ (2 * n + 3))
      (p ^ 2 * q ^ 2 * (1 - p ^ 2) * (1 - q ^ 2) * v ^ 2 *
        (-((n : ℝ) + 3) * v))
      (dualD p q v) (2 * n + 7) 1 (by omega)
    have h2 := mul_div_pow_add_mul_pow_cancel
      (dualCertBaseConstant n p q) (v ^ (2 * n + 3))
      (p ^ 2 * q ^ 2 * (1 - p ^ 2) * (1 - q ^ 2) * v ^ 2 *
        (2 * ((n : ℝ) + 3) * (2 * (n : ℝ) + 7) * v ^ 2))
      (dualD p q v) (2 * n + 8) 0 (by omega)
    simp only [pow_one, pow_zero, mul_one] at h1 h2
    rw [show 2 * n + 8 = (2 * n + 7) + 1 by omega]
    calc
      _ =
          dualCertBaseConstant n p q *
              (v ^ (2 * n + 3) / dualD p q v ^ (2 * n + 7 + 1)) *
                ((p ^ 2 * q ^ 2 * (1 - p ^ 2) * (1 - q ^ 2) * v ^ 2 *
                    (-((n : ℝ) + 3) * v)) * dualD p q v) +
          dualCertBaseConstant n p q *
              (v ^ (2 * n + 3) / dualD p q v ^ (2 * n + 8 + 0)) *
                (p ^ 2 * q ^ 2 * (1 - p ^ 2) * (1 - q ^ 2) * v ^ 2 *
                  (2 * ((n : ℝ) + 3) * (2 * (n : ℝ) + 7) * v ^ 2)) := by
            ring
      _ =
          dualCertBaseConstant n p q *
              (v ^ (2 * n + 3) *
                (p ^ 2 * q ^ 2 * (1 - p ^ 2) * (1 - q ^ 2) * v ^ 2 *
                  (-((n : ℝ) + 3) * v)) / dualD p q v ^ (2 * n + 7)) +
          dualCertBaseConstant n p q *
              (v ^ (2 * n + 3) *
                (p ^ 2 * q ^ 2 * (1 - p ^ 2) * (1 - q ^ 2) * v ^ 2 *
                  (2 * ((n : ℝ) + 3) * (2 * (n : ℝ) + 7) * v ^ 2)) /
                dualD p q v ^ (2 * n + 8)) := by rw [h1, h2]
      _ = _ := by
        rw [← dualCertBaseConstant_succ n p q,
          show 2 * (n + 1) + 3 + 2 = 2 * n + 7 by omega,
          show 2 * n + 7 + 1 = 2 * n + 8 by omega,
          show v ^ (2 * n + 6) = v ^ (2 * n + 3) * v ^ 3 by
            rw [show 2 * n + 6 = (2 * n + 3) + 3 by omega, pow_add],
          show v ^ (2 * n + 7) = v ^ (2 * n + 3) * v ^ 4 by
            rw [show 2 * n + 7 = (2 * n + 3) + 4 by omega, pow_add]]
        ring

theorem dualCertWeightedCurIntegrand_eq (n : ℕ) (i : Fin 3)
    (x : ℝ × (ℝ × ℝ)) :
    dualCertWeightedCurIntegrand n i x =
      ![rawMomentIntegrand n 0 0 0 0 x,
        2 * ((n : ℝ) + 2) * rawMomentIntegrand n 0 0 1 1 x,
        -((n : ℝ) + 2) * rawMomentIntegrand n 0 0 1 1 x +
          2 * ((n : ℝ) + 2) * (2 * (n : ℝ) + 5) *
            rawMomentIntegrand n 0 0 2 2 x] i := by
  rcases x with ⟨p, q, v⟩
  fin_cases i
  · simp [dualCertWeightedCurIntegrand, dualCertCurNum]
    rw [rawMomentIntegrand_zero_factor n 0 4,
      rawMomentIntegrand_zero_factor n 0 0,
      dualCertD_eq_dualD]
    simp only [Nat.add_zero]
    rw [show 2 * n + 4 + 4 = (2 * n + 4) + 4 by omega]
    simpa using mul_div_pow_add_mul_pow_cancel
      (dualCertBaseConstant n p q) (v ^ (2 * n + 3)) 1
      (dualD p q v) (2 * n + 4) 4 (by omega)
  · simp [dualCertWeightedCurIntegrand, dualCertCurNum]
    rw [rawMomentIntegrand_zero_factor n 0 4,
      rawMomentIntegrand_zero_factor n 1 1,
      dualCertD_eq_dualD]
    simp only [Nat.add_zero]
    rw [show 2 * n + 4 + 4 = (2 * n + 5) + 3 by omega]
    rw [mul_div_pow_add_mul_pow_cancel
      (dualCertBaseConstant n p q) (v ^ (2 * n + 3))
      (2 * ((n : ℝ) + 2) * v) (dualD p q v) (2 * n + 5) 3 (by omega)]
    rw [show v ^ (2 * n + 3 + 1) = v ^ (2 * n + 3) * v by
      rw [pow_succ]]
    ring
  · simp [dualCertWeightedCurIntegrand, dualCertCurNum]
    rw [rawMomentIntegrand_zero_factor n 0 4,
      rawMomentIntegrand_zero_factor n 1 1,
      rawMomentIntegrand_zero_factor n 2 2,
      dualCertD_eq_dualD]
    simp only [Nat.add_zero]
    have h1 := mul_div_pow_add_mul_pow_cancel
      (dualCertBaseConstant n p q) (v ^ (2 * n + 3))
      (-((n : ℝ) + 2) * v) (dualD p q v) (2 * n + 5) 3 (by omega)
    have h2 := mul_div_pow_add_mul_pow_cancel
      (dualCertBaseConstant n p q) (v ^ (2 * n + 3))
      (2 * ((n : ℝ) + 2) * (2 * (n : ℝ) + 5) * v ^ 2)
      (dualD p q v) (2 * n + 6) 2 (by omega)
    rw [show 2 * n + 4 + 4 = (2 * n + 5) + 3 by omega]
    calc
      _ =
          dualCertBaseConstant n p q *
              (v ^ (2 * n + 3) / dualD p q v ^ (2 * n + 5 + 3)) *
                ((-((n : ℝ) + 2) * v) * dualD p q v ^ 3) +
          dualCertBaseConstant n p q *
              (v ^ (2 * n + 3) / dualD p q v ^ (2 * n + 6 + 2)) *
                ((2 * ((n : ℝ) + 2) * (2 * (n : ℝ) + 5) * v ^ 2) *
                  dualD p q v ^ 2) := by ring
      _ =
          dualCertBaseConstant n p q *
              (v ^ (2 * n + 3) * (-((n : ℝ) + 2) * v) /
                dualD p q v ^ (2 * n + 5)) +
          dualCertBaseConstant n p q *
              (v ^ (2 * n + 3) *
                (2 * ((n : ℝ) + 2) * (2 * (n : ℝ) + 5) * v ^ 2) /
                dualD p q v ^ (2 * n + 6)) := by rw [h1, h2]
      _ = _ := by
        rw [show v ^ (2 * n + 3 + 1) = v ^ (2 * n + 3) * v by rw [pow_succ],
          show v ^ (2 * n + 3 + 2) = v ^ (2 * n + 3) * v ^ 2 by
            rw [pow_add]]
        ring

def dualVectorIntegrand (n : ℕ) (i : Fin 3)
    (x : ℝ × (ℝ × ℝ)) : ℝ :=
  ![rawMomentIntegrand n 0 0 0 0 x,
    2 * ((n : ℝ) + 2) * rawMomentIntegrand n 0 0 1 1 x,
    -((n : ℝ) + 2) * rawMomentIntegrand n 0 0 1 1 x +
      2 * ((n : ℝ) + 2) * (2 * (n : ℝ) + 5) *
        rawMomentIntegrand n 0 0 2 2 x] i

theorem dualCertWeightedNextIntegrand_eq_dualVectorIntegrand
    (n : ℕ) (j : Fin 3) (x : ℝ × (ℝ × ℝ)) :
    dualCertWeightedNextIntegrand n j x = dualVectorIntegrand (n + 1) j x := by
  rw [dualCertWeightedNextIntegrand_eq]
  fin_cases j
  · rfl
  · change 2 * ((n : ℝ) + 3) * rawMomentIntegrand (n + 1) 0 0 1 1 x =
      2 * (((n + 1 : ℕ) : ℝ) + 2) * rawMomentIntegrand (n + 1) 0 0 1 1 x
    have h : (((n + 1 : ℕ) : ℝ) + 2) = (n : ℝ) + 3 := by
      push_cast
      ring
    rw [h]
  · change
      -((n : ℝ) + 3) * rawMomentIntegrand (n + 1) 0 0 1 1 x +
          2 * ((n : ℝ) + 3) * (2 * (n : ℝ) + 7) *
            rawMomentIntegrand (n + 1) 0 0 2 2 x =
        -(((n + 1 : ℕ) : ℝ) + 2) * rawMomentIntegrand (n + 1) 0 0 1 1 x +
          2 * (((n + 1 : ℕ) : ℝ) + 2) *
            (2 * ((n + 1 : ℕ) : ℝ) + 5) *
              rawMomentIntegrand (n + 1) 0 0 2 2 x
    have hA : (((n + 1 : ℕ) : ℝ) + 2) = (n : ℝ) + 3 := by
      push_cast
      ring
    have hB : 2 * ((n + 1 : ℕ) : ℝ) + 5 = 2 * (n : ℝ) + 7 := by
      push_cast
      ring
    rw [hA, hB]

theorem dualCertWeightedCurIntegrand_eq_dualVectorIntegrand
    (n : ℕ) (i : Fin 3) (x : ℝ × (ℝ × ℝ)) :
    dualCertWeightedCurIntegrand n i x = dualVectorIntegrand n i x := by
  simpa [dualVectorIntegrand] using dualCertWeightedCurIntegrand_eq n i x

theorem dualVectorIntegrand_integrable (n : ℕ) (i : Fin 3) :
    Integrable (dualVectorIntegrand n i) cubeMeasure := by
  have h00 := rawMomentIntegrand_integrable n 0 0 0 0 (by omega)
  have h11 := rawMomentIntegrand_integrable n 0 0 1 1 (by omega)
  have h22 := rawMomentIntegrand_integrable n 0 0 2 2 (by omega)
  fin_cases i
  · simpa [dualVectorIntegrand] using h00
  · simpa [dualVectorIntegrand] using
      h11.const_mul (2 * ((n : ℝ) + 2))
  · have h := (h11.const_mul (-((n : ℝ) + 2))).add
        (h22.const_mul (2 * ((n : ℝ) + 2) * (2 * (n : ℝ) + 5)))
    refine h.congr (Eventually.of_forall fun x => ?_)
    simp [dualVectorIntegrand]

theorem dualVectorIntegrand_integral (n : ℕ) (i : Fin 3) :
    (∫ x, dualVectorIntegrand n i x ∂cubeMeasure) = dualVector n i := by
  fin_cases i
  · simp [dualVectorIntegrand, dualVector, dualMoment]
  · simp [dualVectorIntegrand, dualVector, dualMoment,
      MeasureTheory.integral_const_mul]
  · have h11 := rawMomentIntegrand_integrable n 0 0 1 1 (by omega)
    have h22 := rawMomentIntegrand_integrable n 0 0 2 2 (by omega)
    change (∫ x,
        -((n : ℝ) + 2) * rawMomentIntegrand n 0 0 1 1 x +
          2 * ((n : ℝ) + 2) * (2 * (n : ℝ) + 5) *
            rawMomentIntegrand n 0 0 2 2 x ∂cubeMeasure) =
      -((n : ℝ) + 2) * dualMoment n 0 0 1 1 +
        2 * ((n : ℝ) + 2) * (2 * (n : ℝ) + 5) * dualMoment n 0 0 2 2
    rw [MeasureTheory.integral_add
      (h11.const_mul (-((n : ℝ) + 2)))
      (h22.const_mul (2 * ((n : ℝ) + 2) * (2 * (n : ℝ) + 5))),
      MeasureTheory.integral_const_mul, MeasureTheory.integral_const_mul]
    rfl

theorem dualCertWeightedNextIntegrand_integrable (n : ℕ) (j : Fin 3) :
    Integrable (dualCertWeightedNextIntegrand n j) cubeMeasure :=
  (dualVectorIntegrand_integrable (n + 1) j).congr
    (Eventually.of_forall fun x =>
      (dualCertWeightedNextIntegrand_eq_dualVectorIntegrand n j x).symm)

theorem dualCertWeightedCurIntegrand_integrable (n : ℕ) (i : Fin 3) :
    Integrable (dualCertWeightedCurIntegrand n i) cubeMeasure :=
  (dualVectorIntegrand_integrable n i).congr
    (Eventually.of_forall fun x =>
      (dualCertWeightedCurIntegrand_eq_dualVectorIntegrand n i x).symm)

theorem dualCertWeightedNextIntegrand_integral (n : ℕ) (j : Fin 3) :
    (∫ x, dualCertWeightedNextIntegrand n j x ∂cubeMeasure) =
      dualVector (n + 1) j := by
  rw [MeasureTheory.integral_congr_ae (Eventually.of_forall fun x =>
    dualCertWeightedNextIntegrand_eq_dualVectorIntegrand n j x)]
  exact dualVectorIntegrand_integral (n + 1) j

theorem dualCertWeightedCurIntegrand_integral (n : ℕ) (i : Fin 3) :
    (∫ x, dualCertWeightedCurIntegrand n i x ∂cubeMeasure) =
      dualVector n i := by
  rw [MeasureTheory.integral_congr_ae (Eventually.of_forall fun x =>
    dualCertWeightedCurIntegrand_eq_dualVectorIntegrand n i x)]
  exact dualVectorIntegrand_integral n i

end RamanujanChallenge.P25

end
