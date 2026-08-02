import Mathlib

example (x : ℝ) : HasDerivAt (fun y : ℝ => y ^ 2 + 3 * y) (2 * x + 3) x := by
  fun_prop

example (x : ℝ) (hx : x ≠ 0) :
    HasDerivAt (fun y : ℝ => y ^ 3 / y ^ 2 + Real.log y)
      ((3 * x ^ 2 * x ^ 2 - x ^ 3 * (2 * x)) / (x ^ 2) ^ 2 + x⁻¹) x := by
  fun_prop
