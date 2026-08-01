/-
  Problem 2.5 — final, assumption-transparent closing step.

  The projective contraction already proves that all three challenge ratios
  have the same limit.  Consequently, a proof that the first-column Catalan
  error is negligible relative to its denominator closes the challenge.
  The Delannoy/Poincare route is intended to supply exactly that premise.
-/
import RamanujanChallenge.Problem25Connection
import RamanujanChallenge.Problem25EpsilonDecay

noncomputable section

namespace RamanujanChallenge.P25

/-- The normalized Delannoy error ratio is the negative of the original
first-column Catalan-error ratio; the Pochhammer gauge cancels exactly. -/
theorem first_column_error_ratio_eq_normalized (N : ℕ) :
    catalanError N 0 / (denominator N 0 : ℝ) =
      -(((normalizedNumerator N : ℝ) -
          catalanConstant * (normalizedDenominator N : ℝ)) /
        (normalizedDenominator N : ℝ)) := by
  have hp : (pochhammerProduct N : ℝ) ≠ 0 := by
    exact_mod_cast pochhammerProduct_ne_zero N
  have hq : (denominator N 0 : ℝ) ≠ 0 := by
    exact_mod_cast denominator_ne_zero N 0
  simp only [catalanError, normalizedNumerator, normalizedDenominator]
  push_cast
  field_simp [hp, hq]
  ring

/-- Once the first-column relative Catalan error tends to zero, the existing
projective-contraction theorem closes all three columns. -/
theorem problem25Claim_of_first_column_error_decay
    (herror : Filter.Tendsto
      (fun N : ℕ => catalanError N 0 / (denominator N 0 : ℝ))
      Filter.atTop (nhds 0)) :
    Problem25Claim := by
  have hratio : Filter.Tendsto (fun N : ℕ => challengeRatio N 0)
      Filter.atTop (nhds catalanConstant) := by
    have heq :
        (fun N : ℕ => challengeRatio N 0) =
          (fun N : ℕ =>
            catalanConstant - catalanError N 0 / (denominator N 0 : ℝ)) := by
      funext N
      have hq : (denominator N 0 : ℝ) ≠ 0 := by
        exact_mod_cast denominator_ne_zero N 0
      simp only [challengeRatio, catalanError]
      field_simp [hq]
      ring
    rw [heq]
    simpa using
      Filter.Tendsto.sub
        (tendsto_const_nhds :
          Filter.Tendsto (fun _ : ℕ => catalanConstant)
            Filter.atTop (nhds catalanConstant))
        herror
  have hcommon : commonLimit = catalanConstant :=
    tendsto_nhds_unique (challengeRatio_tendsto_common 0) hratio
  rw [problem25Claim_iff_commonLimit_eq_catalan]
  exact hcommon

/-- Final interface for the Delannoy route: it is enough to prove decay of
the normalized weighted error supplied by `normalized_error_decomposition`. -/
theorem problem25Claim_of_normalized_delannoy_error_decay
    (herror : Filter.Tendsto
      (fun N : ℕ =>
        ((normalizedNumerator N : ℝ) -
            catalanConstant * (normalizedDenominator N : ℝ)) /
          (normalizedDenominator N : ℝ))
      Filter.atTop (nhds 0)) :
    Problem25Claim := by
  apply problem25Claim_of_first_column_error_decay
  have hneg : Filter.Tendsto
      (fun N : ℕ =>
        -(((normalizedNumerator N : ℝ) -
            catalanConstant * (normalizedDenominator N : ℝ)) /
          (normalizedDenominator N : ℝ)))
      Filter.atTop (nhds 0) := by
    simpa using herror.neg
  apply hneg.congr'
  filter_upwards [] with N
  exact (first_column_error_ratio_eq_normalized N).symm

end RamanujanChallenge.P25

end
