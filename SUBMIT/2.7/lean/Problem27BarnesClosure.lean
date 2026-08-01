import RamanujanChallenge.Problem27BarnesRecurrence
import RamanujanChallenge.Problem27BarnesInitial

/-!
# Problem 2.7: unconditional Barnes normalization and challenge theorem

The analytic Barnes sequence and the recurrence-defined source error satisfy
the same third-order recurrence and the same first three values.  Recurrence
uniqueness therefore supplies the exact normalization omitted in the source
paper, after which the already proved Barnes decay closes the challenge.
-/

open Filter Topology

noncomputable section

namespace RamanujanChallenge.P27

private theorem ctAlpha_ne_zeroC27 (k : ℕ) :
    (ctAlpha27 k : ℂ) ≠ 0 := by
  have hq : ctAlpha27 k ≠ 0 := by
    simp only [ctAlpha27]
    exact mul_ne_zero (by norm_num) (ne_of_gt (zudilinA_pos k))
  exact_mod_cast hq

theorem ZudilinSatisfiesRecC27.eq_of_initial
    {u v : ℕ → ℂ}
    (hu : ZudilinSatisfiesRecC27 u)
    (hv : ZudilinSatisfiesRecC27 v)
    (h0 : u 0 = v 0)
    (h1 : u 1 = v 1)
    (h2 : u 2 = v 2) :
    u = v := by
  funext n
  induction n using Nat.strong_induction_on with
  | h n ih =>
      match n with
      | 0 => exact h0
      | 1 => exact h1
      | 2 => exact h2
      | k + 3 =>
          have huk := hu k
          have hvk := hv k
          rw [ih (k + 2) (by omega),
            ih (k + 1) (by omega), ih k (by omega)] at huk
          apply mul_left_cancel₀ (ctAlpha_ne_zeroC27 k)
          linear_combination huk - hvk

theorem ZudilinSatisfiesRec.castComplex27
    {u : ℕ → ℚ} (hu : ZudilinSatisfiesRec u) :
    ZudilinSatisfiesRecC27 (fun n => (u n : ℂ)) := by
  intro k
  have hq := (zudilinSatisfiesRec_iff_ct27 u).1 hu k
  have hc := congrArg (fun q : ℚ => (q : ℂ)) hq
  push_cast at hc
  exact hc

def zudilinL27 : ℂ :=
  ((Real.pi ^ 2 / 6 + zeta3 : ℝ) : ℂ)

def zudilinSourceErrorC27 (n : ℕ) : ℂ :=
  zudilinL27 * (zudilinB n : ℂ) - (zudilinM23 n : ℂ)

theorem zudilinB_recC27 :
    ZudilinSatisfiesRecC27 (fun n => (zudilinB n : ℂ)) :=
  ZudilinSatisfiesRec.castComplex27 zudilinB_rec

theorem zudilinM23_recC27 :
    ZudilinSatisfiesRecC27 (fun n => (zudilinM23 n : ℂ)) :=
  ZudilinSatisfiesRec.castComplex27 zudilinM23_rec

theorem zudilinSourceErrorC_rec27 :
    ZudilinSatisfiesRecC27 zudilinSourceErrorC27 := by
  intro k
  have hB := zudilinB_recC27 k
  have hM := zudilinM23_recC27 k
  simp only [zudilinSourceErrorC27]
  linear_combination zudilinL27 * hB - hM

theorem zudilinSourceErrorC_zero27 :
    zudilinSourceErrorC27 0 = zudilinL27 := by
  norm_num [zudilinSourceErrorC27, zudilinL27, zudilinB,
    zudilinM23, zudilinB₂, zudilinB₃, zudilinSolution]

theorem zudilinSourceErrorC_one27 :
    zudilinSourceErrorC27 1 = 7 * zudilinL27 - 20 := by
  norm_num [zudilinSourceErrorC27, zudilinL27, zudilinB,
    zudilinM23, zudilinB₂, zudilinB₃, zudilinSolution]
  ring

theorem zudilinSourceErrorC_two27 :
    zudilinSourceErrorC27 2 =
      163 * zudilinL27 - (7425 / 16 : ℂ) := by
  norm_num [zudilinSourceErrorC27, zudilinL27, zudilinB,
    zudilinM23, zudilinB₂, zudilinB₃, zudilinSolution]
  ring

theorem coe_zudilinCombinedError_eq_sourceErrorC27 (n : ℕ) :
    (zudilinCombinedError n : ℂ) = zudilinSourceErrorC27 n := by
  simp only [zudilinCombinedError, zudilinSourceErrorC27,
    zudilinL27, zudilinBReal, zudilinM23Real]
  push_cast
  ring

theorem zudilinBarnesErrorIntegral_zeroC27 :
    zudilinBarnesErrorIntegral27 0 = zudilinL27 := by
  simpa only [zudilinL27] using zudilinBarnesErrorIntegral_zero27

theorem zudilinBarnesErrorIntegral_oneC27 :
    zudilinBarnesErrorIntegral27 1 = 7 * zudilinL27 - 20 := by
  rw [zudilinBarnesErrorIntegral_one27]
  simp only [zudilinL27]
  push_cast
  ring

theorem zudilinBarnesErrorIntegral_twoC27 :
    zudilinBarnesErrorIntegral27 2 =
      163 * zudilinL27 - (7425 / 16 : ℂ) := by
  rw [zudilinBarnesErrorIntegral_two27]
  simp only [zudilinL27]
  push_cast
  ring

/-- Exact, assumption-free identification of the recurrence-defined source
error with the analytic Barnes integral. -/
theorem zudilinCombinedError_eq_barnes27 (n : ℕ) :
    (zudilinCombinedError n : ℂ) =
      zudilinBarnesErrorIntegral27 n := by
  have hfun : zudilinSourceErrorC27 = zudilinBarnesErrorIntegral27 :=
    ZudilinSatisfiesRecC27.eq_of_initial
      zudilinSourceErrorC_rec27 zudilinBarnesErrorIntegral_recC27
      (zudilinSourceErrorC_zero27.trans
        zudilinBarnesErrorIntegral_zeroC27.symm)
      (zudilinSourceErrorC_one27.trans
        zudilinBarnesErrorIntegral_oneC27.symm)
      (zudilinSourceErrorC_two27.trans
        zudilinBarnesErrorIntegral_twoC27.symm)
  rw [coe_zudilinCombinedError_eq_sourceErrorC27, hfun]

theorem zudilinCombinedError_tendsto_zero_direct27 :
    Tendsto zudilinCombinedError atTop (𝓝 0) :=
  zudilinCombinedError_tendsto_of_barnes_representation27
    zudilinCombinedError_eq_barnes27

/-- Direct theorem corresponding to the challenge statement, with no analytic
or normalization hypotheses. -/
theorem problem27_zeta2_add_zeta3_direct :
    Tendsto
      (fun n => (challengeP n : ℝ) / (challengeQ n : ℝ))
      atTop (𝓝 (Real.pi ^ 2 / 6 + zeta3)) :=
  problem27_zeta2_add_zeta3_of_barnes_representation
    zudilinCombinedError_eq_barnes27

#print axioms ctRKernelRaw_one_strip27
#print axioms ctSKernelRaw_one_strip27
#print axioms ctRPhiVerticalIntegral_telescoper27
#print axioms zudilinBarnesErrorIntegral_recC27
#print axioms zudilinCombinedError_eq_barnes27
#print axioms problem27_zeta2_add_zeta3_direct

end RamanujanChallenge.P27
