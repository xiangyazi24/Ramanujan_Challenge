/-
  Sign verification of positiveCatalanError at N=1.
  Uses the tight Catalan bounds to show G is trapped in the N=1 envelope.
-/
import RamanujanChallenge.Problem25TightBounds

noncomputable section

namespace RamanujanChallenge.P25

set_option maxHeartbeats 0 in
set_option maxRecDepth 2000 in
private theorem positiveDenominator_one :
    (fun j => positiveDenominator 1 j) = ![23352603750, 14261609250, 1736437500] := by
  funext j
  fin_cases j <;>
    norm_num [positiveDenominator, denominator, approximants, matrixProduct,
      initialMatrix, challengeMatrix, m11, m12, m13, m21, m22, m23, m31, m32, m33,
      coordinateSign, Matrix.mul_apply, Fin.sum_univ_succ]

set_option maxHeartbeats 0 in
set_option maxRecDepth 2000 in
private theorem positiveNumerator_one :
    (fun j => positiveNumerator 1 j) = ![21390206625, 13063139595, 1590511050] := by
  funext j
  fin_cases j <;>
    norm_num [positiveNumerator, numerator, approximants, matrixProduct,
      initialMatrix, challengeMatrix, m11, m12, m13, m21, m22, m23, m31, m32, m33,
      coordinateSign, Matrix.mul_apply, Fin.sum_univ_succ]

/-- At N=1, the j=0 component of the Catalan error is negative. -/
theorem positiveCatalanError_one_zero_neg : positiveCatalanError 1 0 < 0 := by
  rw [positiveCatalanError_eq]
  have hq := congrFun positiveDenominator_one (0 : Fin 3)
  have hp := congrFun positiveNumerator_one (0 : Fin 3)
  norm_num at hq hp
  rw [hq, hp]
  have h := catalan_tight_upper
  push_cast
  nlinarith

/-- At N=1, the j=2 component of the Catalan error is positive. -/
theorem positiveCatalanError_one_two_pos : 0 < positiveCatalanError 1 2 := by
  rw [positiveCatalanError_eq]
  have hq := congrFun positiveDenominator_one (2 : Fin 3)
  have hp := congrFun positiveNumerator_one (2 : Fin 3)
  norm_num [Matrix.cons_val_two] at hq hp
  rw [hq, hp]
  have h := catalan_tight_lower
  push_cast
  nlinarith

end RamanujanChallenge.P25

end
