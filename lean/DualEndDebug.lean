import RamanujanChallenge.Problem25Integral
import Mathlib.MeasureTheory.Integral.Prod

noncomputable section
namespace RamanujanChallenge.P25
open MeasureTheory

private abbrev debugMeasure : Measure (ℝ × (ℝ × ℝ)) := volume
private opaque debugRaw (n A B C k : ℕ) : ℝ × (ℝ × ℝ) → ℝ
private def debugMoment (n A B C k : ℕ) : ℝ :=
  ∫ x, debugRaw n A B C k x ∂debugMeasure
private def debugCombination (n : ℕ) (p q v : ℝ) : ℝ :=
  (2 * (n : ℝ) + 5) * (2 * debugRaw n 0 0 1 1 (p, q, v)) -
    (4 * (n : ℝ) + 11) * (4 * debugRaw n 0 0 2 2 (p, q, v)) +
    (2 * (n : ℝ) + 6) * 8 *
      (debugRaw n 0 0 3 3 (p, q, v) - debugRaw n 2 2 3 3 (p, q, v))
private def debugVector (n : ℕ) : Fin 3 → ℝ :=
  ![debugMoment n 0 0 0 0,
    2 * ((n : ℝ) + 2) * debugMoment n 0 0 1 1,
    -((n : ℝ) + 2) * debugMoment n 0 0 1 1 +
      2 * ((n : ℝ) + 2) * (2 * (n : ℝ) + 5) * debugMoment n 0 0 2 2]

private axiom debug_integrable (n A B C k : ℕ) :
  Integrable (debugRaw n A B C k) debugMeasure
private axiom debug_zero (n : ℕ) :
  (∫ x, debugCombination n x.1 x.2.1 x.2.2 ∂debugMeasure) = 0
private axiom debug_moment_pos (n A B C k : ℕ) : 0 < debugMoment n A B C k
private axiom debug_delta_pos (n : ℕ) :
  0 < debugMoment n 0 0 3 3 - debugMoment n 2 2 3 3

example (n : ℕ) : 0 < debugVector n 2 := by
  have hzero := debug_zero n
  have h11 := debug_integrable n 0 0 1 1
  have h22 := debug_integrable n 0 0 2 2
  have h33 := debug_integrable n 0 0 3 3
  have hshift := debug_integrable n 2 2 3 3
  have h11' : Integrable (fun x : ℝ × (ℝ × ℝ) =>
      debugRaw n 0 0 1 1 (x.1, x.2.1, x.2.2)) debugMeasure := by
    simpa using h11
  have h22' : Integrable (fun x : ℝ × (ℝ × ℝ) =>
      debugRaw n 0 0 2 2 (x.1, x.2.1, x.2.2)) debugMeasure := by
    simpa using h22
  have h33' : Integrable (fun x : ℝ × (ℝ × ℝ) =>
      debugRaw n 0 0 3 3 (x.1, x.2.1, x.2.2)) debugMeasure := by
    simpa using h33
  have hshift' : Integrable (fun x : ℝ × (ℝ × ℝ) =>
      debugRaw n 2 2 3 3 (x.1, x.2.1, x.2.2)) debugMeasure := by
    simpa using hshift
  have hterm1 := (h11'.const_mul 2).const_mul (2 * (n : ℝ) + 5)
  have hterm2 := (h22'.const_mul 4).const_mul (4 * (n : ℝ) + 11)
  have hleft := hterm1.sub hterm2
  have hright := ((h33'.sub hshift').const_mul 8).const_mul (2 * (n : ℝ) + 6)
  simp only [debugCombination] at hzero
  simp only [mul_assoc] at hzero
  have hzero1 :
      (∫ x, (2 * (n : ℝ) + 5) * (2 * debugRaw n 0 0 1 1 (x.1, x.2.1, x.2.2)) -
          (4 * (n : ℝ) + 11) * (4 * debugRaw n 0 0 2 2 (x.1, x.2.1, x.2.2)) ∂debugMeasure) +
        (∫ x, (2 * (n : ℝ) + 6) * (8 *
          (debugRaw n 0 0 3 3 (x.1, x.2.1, x.2.2) -
            debugRaw n 2 2 3 3 (x.1, x.2.1, x.2.2))) ∂debugMeasure) = 0 := by
    calc
      _ = ∫ x, (((fun y : ℝ × (ℝ × ℝ) =>
            (2 * (n : ℝ) + 5) * (2 * debugRaw n 0 0 1 1 (y.1, y.2.1, y.2.2))) -
          (fun y : ℝ × (ℝ × ℝ) =>
            (4 * (n : ℝ) + 11) * (4 * debugRaw n 0 0 2 2 (y.1, y.2.1, y.2.2)))) +
          (fun y : ℝ × (ℝ × ℝ) => (2 * (n : ℝ) + 6) * (8 *
            (debugRaw n 0 0 3 3 (y.1, y.2.1, y.2.2) -
              debugRaw n 2 2 3 3 (y.1, y.2.1, y.2.2))))) x ∂debugMeasure :=
        (MeasureTheory.integral_add hleft hright).symm
      _ = 0 := by simpa only [Pi.add_apply, Pi.sub_apply, mul_assoc] using hzero
  rw [MeasureTheory.integral_sub hterm1 hterm2] at hzero1
  simp only [MeasureTheory.integral_const_mul] at hzero1
  rw [MeasureTheory.integral_sub h33' hshift'] at hzero1
  change
    (2 * (n : ℝ) + 5) * (2 * debugMoment n 0 0 1 1) -
        (4 * (n : ℝ) + 11) * (4 * debugMoment n 0 0 2 2) +
      (2 * (n : ℝ) + 6) *
        (8 * (debugMoment n 0 0 3 3 - debugMoment n 2 2 3 3)) = 0 at hzero1
  have hmoment : 0 < debugMoment n 0 0 2 2 := debug_moment_pos n 0 0 2 2
  have hdelta := debug_delta_pos n
  have ht : (0 : ℝ) < 2 * (n : ℝ) + 5 := by positivity
  have hn : (0 : ℝ) ≤ n := by positivity
  let d := -((n : ℝ) + 2) * debugMoment n 0 0 1 1 +
    2 * ((n : ℝ) + 2) * (2 * (n : ℝ) + 5) * debugMoment n 0 0 2 2
  have hid :
      (2 * (n : ℝ) + 5) * d = ((n : ℝ) + 2) *
        ((8 * (n : ℝ) ^ 2 + 32 * (n : ℝ) + 28) * debugMoment n 0 0 2 2 +
          (8 * (n : ℝ) + 24) *
            (debugMoment n 0 0 3 3 - debugMoment n 2 2 3 3)) := by
    dsimp [d]
    linear_combination -(((n : ℝ) + 2) / 2) * hzero1
  have hcoef : 0 < 8 * (n : ℝ) ^ 2 + 32 * (n : ℝ) + 28 := by
    nlinarith [sq_nonneg (n : ℝ)]
  have hrhs : 0 < ((n : ℝ) + 2) *
      ((8 * (n : ℝ) ^ 2 + 32 * (n : ℝ) + 28) * debugMoment n 0 0 2 2 +
        (8 * (n : ℝ) + 24) *
          (debugMoment n 0 0 3 3 - debugMoment n 2 2 3 3)) := by
    exact mul_pos (by positivity)
      (add_pos (mul_pos hcoef hmoment) (mul_pos (by positivity) hdelta))
  have hprod : 0 < (2 * (n : ℝ) + 5) * d := hid.symm ▸ hrhs
  have hthird : 0 < d := pos_of_mul_pos_right hprod ht.le
  dsimp [d] at hthird
  simpa [debugVector] using hthird

end RamanujanChallenge.P25
end
