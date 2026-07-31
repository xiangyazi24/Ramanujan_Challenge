import RamanujanChallenge.Problem26Generating
import RamanujanChallenge.Problem26Cyclotomic

/-!
# Problem 2.6: analytic assembly

This module joins the generating-function bridge with the real
weight-three and cyclotomic reductions.
-/

namespace RamanujanChallenge.P26

/-- The exact Problem 2.6 assertion follows from the single remaining
cyclotomic logarithmic special value. -/
theorem problem26_of_cyclotomic_log_integral
    (hCyclotomic : CyclotomicLogIntegralEvaluation26) :
    Problem26Claim := by
  apply problem26_of_nested_inverse_binomial_evaluation
  convert inverseBinomialDTerm26_hasSum_cyclotomic_integral using 1
  rw [nestedCyclotomicIntegral26_of_cyclotomicLog hCyclotomic]
  ring_nf

/-- The exact Problem 2.6 challenge assertion. -/
theorem problem26 : Problem26Claim :=
  problem26_of_cyclotomic_log_integral
    cyclotomicLogIntegralEvaluation26

end RamanujanChallenge.P26
