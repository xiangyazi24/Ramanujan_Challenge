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

/-- The printed series converges to the claimed value for every rational
sequence with the two stated initial values and the stated three-term
recurrence.  `HasSum` records convergence as well as the value of the sum. -/
theorem problem26_hasSum_of_spec
    (u : ℕ → ℚ)
    (hu₁ : u 1 = -93 / 4480)
    (hu₂ : u 2 = -117 / 14000)
    (hu : SatisfiesRecurrence26 u) :
    HasSum (fun n : ℕ => ((u (n + 1) : ℚ) : ℝ))
      (Real.pi ^ 2 / 6 + zeta3 - (2077 : ℝ) / 720) := by
  have hterms (n : ℕ) :
      ((u (n + 1) : ℚ) : ℝ) =
        ((challengeU26 (n + 1) : ℚ) : ℝ) := by
    exact_mod_cast eq_challengeU26_of_spec u hu₁ hu₂ hu n
  have hvalue :
      (∑' n : ℕ, ((challengeU26 (n + 1) : ℚ) : ℝ)) =
        Real.pi ^ 2 / 6 + zeta3 - (2077 : ℝ) / 720 := by
    have h := problem26
    unfold Problem26Claim at h
    linarith
  have hchallenge :
      HasSum (fun n : ℕ => ((challengeU26 (n + 1) : ℚ) : ℝ))
        (Real.pi ^ 2 / 6 + zeta3 - (2077 : ℝ) / 720) := by
    rw [← hvalue]
    exact challengeU26_summable.hasSum
  exact hchallenge.congr_fun hterms

/-- The challenge identity in the displayed-equality formulation of the
printed problem. -/
theorem problem26_of_spec
    (u : ℕ → ℚ)
    (hu₁ : u 1 = -93 / 4480)
    (hu₂ : u 2 = -117 / 14000)
    (hu : SatisfiesRecurrence26 u) :
    (2077 : ℝ) / 720 +
        (∑' n : ℕ, ((u (n + 1) : ℚ) : ℝ)) =
      Real.pi ^ 2 / 6 + zeta3 := by
  rw [(problem26_hasSum_of_spec u hu₁ hu₂ hu).tsum_eq]
  ring

end RamanujanChallenge.P26
