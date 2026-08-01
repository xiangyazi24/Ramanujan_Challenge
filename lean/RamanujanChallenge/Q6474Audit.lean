import RamanujanChallenge.Problem27BarnesTelescoper

open Filter Topology

noncomputable section

namespace RamanujanChallenge.P27.Q6474

/-- Homogeneous Zudilin source recurrence over `ℂ`.  Using the homogeneous
form avoids all division/cast issues in the uniqueness bridge. -/
def ZudilinSatisfiesRecC27 (u : ℕ → ℂ) : Prop :=
  ∀ k : ℕ,
    (ctAlpha27 k : ℂ) * u (k + 3)
      - (ctBeta27 k : ℂ) * u (k + 2)
      + (ctGamma27 k : ℂ) * u (k + 1)
      - (ctDelta27 k : ℂ) * u k = 0

private theorem ctAlpha_ne_zero27 (k : ℕ) : (ctAlpha27 k : ℂ) ≠ 0 := by
  have hA : zudilinA (k + 2) ≠ 0 := ne_of_gt (zudilinA_pos k)
  have hq : ctAlpha27 k ≠ 0 := by
    simp only [ctAlpha27]
    exact mul_ne_zero (by norm_num) hA
  exact_mod_cast hq

/-- Third-order recurrence uniqueness over `ℂ`. -/
theorem ZudilinSatisfiesRecC27.eq_of_initial
    {u v : ℕ → ℂ}
    (hu : ZudilinSatisfiesRecC27 u)
    (hv : ZudilinSatisfiesRecC27 v)
    (h0 : u 0 = v 0) (h1 : u 1 = v 1) (h2 : u 2 = v 2) :
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
          rw [ih (k + 2) (by omega), ih (k + 1) (by omega),
            ih k (by omega)] at huk
          apply mul_left_cancel₀ (ctAlpha_ne_zero27 k)
          linear_combination huk - hvk

/-- Rational source recurrence transported to `ℂ`. -/
theorem ZudilinSatisfiesRec.castComplex27
    {u : ℕ → ℚ} (hu : ZudilinSatisfiesRec u) :
    ZudilinSatisfiesRecC27 (fun n => (u n : ℂ)) := by
  intro k
  have hq := (zudilinSatisfiesRec_iff_ct27 u).1 hu k
  exact_mod_cast hq

/-- The exact complex constant used by both the source and Barnes sequences. -/
def L27 : ℂ := ((Real.pi ^ 2 / 6 + zeta3 : ℝ) : ℂ)

/-- Complexification of the recurrence-defined combined source error. -/
def sourceErrorC27 (n : ℕ) : ℂ :=
  L27 * (zudilinB n : ℂ) - (zudilinM23 n : ℂ)

theorem zudilinB_recC27 :
    ZudilinSatisfiesRecC27 (fun n => (zudilinB n : ℂ)) :=
  ZudilinSatisfiesRec.castComplex27 zudilinB_rec

theorem zudilinM23_recC27 :
    ZudilinSatisfiesRecC27 (fun n => (zudilinM23 n : ℂ)) :=
  ZudilinSatisfiesRec.castComplex27 zudilinM23_rec

theorem sourceErrorC_rec27 : ZudilinSatisfiesRecC27 sourceErrorC27 := by
  intro k
  have hB := zudilinB_recC27 k
  have hM := zudilinM23_recC27 k
  simp only [sourceErrorC27]
  linear_combination L27 * hB - hM

theorem sourceErrorC_zero27 : sourceErrorC27 0 = L27 := by
  norm_num [sourceErrorC27, L27, zudilinB, zudilinM23,
    zudilinB₂, zudilinB₃, zudilinSolution]

theorem sourceErrorC_one27 : sourceErrorC27 1 = 7 * L27 - 20 := by
  norm_num [sourceErrorC27, L27, zudilinB, zudilinM23,
    zudilinB₂, zudilinB₃, zudilinSolution]
  ring

/-- Audit point: `2145/8 + 3135/16 = 7425/16`, not `7425/8`. -/
theorem sourceErrorC_two27 :
    sourceErrorC27 2 = 163 * L27 - (7425 / 16 : ℂ) := by
  norm_num [sourceErrorC27, L27, zudilinB, zudilinM23,
    zudilinB₂, zudilinB₃, zudilinSolution]
  ring

/-- The repository real source error coerces to the complex source error. -/
theorem coe_zudilinCombinedError_eq_sourceErrorC27 (n : ℕ) :
    (zudilinCombinedError n : ℂ) = sourceErrorC27 n := by
  simp only [zudilinCombinedError, sourceErrorC27, L27,
    zudilinBReal, zudilinM23Real]
  push_cast
  ring

/-- The entire normalization bridge.  The four hypotheses are intended to be
filled by the proved Barnes recurrence and the already compiled exact values
at `0,1,2`. -/
theorem source_eq_barnes_of_rec_initial27
    (hBarnesRec : ZudilinSatisfiesRecC27 zudilinBarnesErrorIntegral27)
    (h0 : zudilinBarnesErrorIntegral27 0 = L27)
    (h1 : zudilinBarnesErrorIntegral27 1 = 7 * L27 - 20)
    (h2 : zudilinBarnesErrorIntegral27 2 =
      163 * L27 - (7425 / 16 : ℂ)) :
    ∀ n : ℕ,
      (zudilinCombinedError n : ℂ) = zudilinBarnesErrorIntegral27 n := by
  have hfun : sourceErrorC27 = zudilinBarnesErrorIntegral27 :=
    ZudilinSatisfiesRecC27.eq_of_initial sourceErrorC_rec27 hBarnesRec
      (sourceErrorC_zero27.trans h0.symm)
      (sourceErrorC_one27.trans h1.symm)
      (sourceErrorC_two27.trans h2.symm)
  intro n
  rw [coe_zudilinCombinedError_eq_sourceErrorC27 n, hfun]

/-- Once the four normalization theorems are installed, Barnes decay closes
`zudilinCombinedError` directly. -/
theorem source_tendsto_zero_of_rec_initial27
    (hBarnesRec : ZudilinSatisfiesRecC27 zudilinBarnesErrorIntegral27)
    (h0 : zudilinBarnesErrorIntegral27 0 = L27)
    (h1 : zudilinBarnesErrorIntegral27 1 = 7 * L27 - 20)
    (h2 : zudilinBarnesErrorIntegral27 2 =
      163 * L27 - (7425 / 16 : ℂ)) :
    Tendsto zudilinCombinedError atTop (𝓝 0) :=
  zudilinCombinedError_tendsto_of_barnes_representation27
    (source_eq_barnes_of_rec_initial27 hBarnesRec h0 h1 h2)

/-- Same bridge all the way to the challenge conclusion. -/
theorem challenge_of_rec_initial27
    (hBarnesRec : ZudilinSatisfiesRecC27 zudilinBarnesErrorIntegral27)
    (h0 : zudilinBarnesErrorIntegral27 0 = L27)
    (h1 : zudilinBarnesErrorIntegral27 1 = 7 * L27 - 20)
    (h2 : zudilinBarnesErrorIntegral27 2 =
      163 * L27 - (7425 / 16 : ℂ)) :
    Tendsto
      (fun n => (challengeP n : ℝ) / (challengeQ n : ℝ))
      atTop (𝓝 (Real.pi ^ 2 / 6 + zeta3)) :=
  problem27_zeta2_add_zeta3_of_barnes_representation
    (source_eq_barnes_of_rec_initial27 hBarnesRec h0 h1 h2)

#print axioms ZudilinSatisfiesRecC27.eq_of_initial
#print axioms sourceErrorC_rec27
#print axioms sourceErrorC_two27
#print axioms source_eq_barnes_of_rec_initial27
#print axioms challenge_of_rec_initial27

end RamanujanChallenge.P27.Q6474
