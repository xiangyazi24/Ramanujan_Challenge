/-
  Sign verification of positiveCatalanError at N=1.
  Uses the tight Catalan bounds to show the (-, +, +) sign pattern persists.
-/
import RamanujanChallenge.Problem25TightBounds

noncomputable section

namespace RamanujanChallenge.P25

open MeasureTheory Set Filter Topology

set_option maxHeartbeats 0 in
private theorem positiveDenominator_one_zero :
    positiveDenominator 1 0 = 23352603750 := by
  rw [positiveDenominator_succ]
  have h0 := congrFun positiveDenominator_zero (0 : Fin 3)
  have h1 := congrFun positiveDenominator_zero (1 : Fin 3)
  have h2 := congrFun positiveDenominator_zero (2 : Fin 3)
  norm_num at h0 h1 h2
  rw [Fin.sum_univ_three, h0, h1, h2]
  norm_num [positiveMatrix, challengeMatrix, m11, m12, m13, m21, m22, m23,
    m31, m32, m33, coordinateSign]

set_option maxHeartbeats 0 in
private theorem positiveNumerator_one_zero :
    positiveNumerator 1 0 = 21390206625 := by
  rw [positiveNumerator_succ]
  have h0 := congrFun positiveNumerator_zero (0 : Fin 3)
  have h1 := congrFun positiveNumerator_zero (1 : Fin 3)
  have h2 := congrFun positiveNumerator_zero (2 : Fin 3)
  norm_num at h0 h1 h2
  rw [Fin.sum_univ_three, h0, h1, h2]
  norm_num [positiveMatrix, challengeMatrix, m11, m12, m13, m21, m22, m23,
    m31, m32, m33, coordinateSign]

set_option maxHeartbeats 0 in
private theorem positiveDenominator_one_two :
    positiveDenominator 1 2 = 1736437500 := by
  rw [positiveDenominator_succ]
  have h0 := congrFun positiveDenominator_zero (0 : Fin 3)
  have h1 := congrFun positiveDenominator_zero (1 : Fin 3)
  have h2 := congrFun positiveDenominator_zero (2 : Fin 3)
  norm_num at h0 h1 h2
  rw [Fin.sum_univ_three, h0, h1, h2]
  norm_num [positiveMatrix, challengeMatrix, m11, m12, m13, m21, m22, m23,
    m31, m32, m33, coordinateSign]

set_option maxHeartbeats 0 in
private theorem positiveNumerator_one_two :
    positiveNumerator 1 2 = 1590511050 := by
  rw [positiveNumerator_succ]
  have h0 := congrFun positiveNumerator_zero (0 : Fin 3)
  have h1 := congrFun positiveNumerator_zero (1 : Fin 3)
  have h2 := congrFun positiveNumerator_zero (2 : Fin 3)
  norm_num at h0 h1 h2
  rw [Fin.sum_univ_three, h0, h1, h2]
  norm_num [positiveMatrix, challengeMatrix, m11, m12, m13, m21, m22, m23,
    m31, m32, m33, coordinateSign]

/-- At N=1, the j=0 component of the Catalan error is negative.
This means catalanConstant < positiveRatio 1 0. -/
theorem positiveCatalanError_one_zero_neg : positiveCatalanError 1 0 < 0 := by
  rw [positiveCatalanError_eq, positiveDenominator_one_zero, positiveNumerator_one_zero]
  have h := catalan_tight_upper
  push_cast
  nlinarith

/-- At N=1, the j=2 component of the Catalan error is positive.
This means catalanConstant > positiveRatio 1 2. -/
theorem positiveCatalanError_one_two_pos : 0 < positiveCatalanError 1 2 := by
  rw [positiveCatalanError_eq, positiveDenominator_one_two, positiveNumerator_one_two]
  have h := catalan_tight_lower
  push_cast
  nlinarith

/-- G is trapped in the N=1 envelope: lowerEnvelope 1 < G < upperEnvelope 1.
Combined with the N=0 signs (already proved) and the contraction, this gives
|G - commonLimit| < envelopeGap 1 ≈ 4.5e-6. -/
theorem catalan_in_envelope_one :
    positiveRatio 1 2 < catalanConstant ∧ catalanConstant < positiveRatio 1 0 := by
  have h0 := positiveCatalanError_one_zero_neg
  have h2 := positiveCatalanError_one_two_pos
  rw [positiveCatalanError_eq] at h0 h2
  have hq0 : (0 : ℝ) < (positiveDenominator 1 0 : ℝ) := by
    exact_mod_cast positiveDenominator_pos 1 0
  have hq2 : (0 : ℝ) < (positiveDenominator 1 2 : ℝ) := by
    exact_mod_cast positiveDenominator_pos 1 2
  simp only [positiveRatio]
  constructor
  · rw [lt_div_iff₀ hq2]; linarith
  · rw [div_lt_iff₀ hq0]; linarith

end RamanujanChallenge.P25

end
