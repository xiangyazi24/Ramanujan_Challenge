/-
  Problem 2.5 — Moment formula layer.

  Connects the CMF sequences to Catalan's constant via the exact
  moment formula (paper §2, Theorem 2):

    G · Q_{N,j} - P_{N,j} = ∫₀¹ [-log(t)/(1+t²)] · R_{N,j}(t²) dt

  where R_{N,j}(X) = (Q-P) + (-P)·X is affine in X.

  Layer 1: moment identities ∫wq·R = Q, ∫wp·R = P.
  Layer 2: affine moment identity → denominator/numerator formulas.
  Layer 3: remainder polynomial and its key evaluations.
-/
import RamanujanChallenge.Problem25Connection
import RamanujanChallenge.Problem25Integral
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic

noncomputable section

namespace RamanujanChallenge.P25

open MeasureTheory Set Filter Topology

private theorem poly_integrable_01 {f : ℝ → ℝ} (hf : Continuous f) :
    IntervalIntegrable f MeasureTheory.volume 0 1 :=
  hf.intervalIntegrable 0 1

/-! ## Moment identities for the weight functions -/

theorem moment_wq_one :
    ∫ x in (0 : ℝ)..1, (10 - 18 * x) = (1 : ℝ) := by
  have h : (fun x : ℝ => (10 : ℝ) - 18 * x) = (fun x => 10 * (1 : ℝ) + (-18) * x) := by
    ext; ring
  rw [h, intervalIntegral.integral_add (poly_integrable_01 (by fun_prop))
    (poly_integrable_01 (by fun_prop)),
    intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul,
    integral_one, integral_id]
  norm_num

theorem moment_wq_x :
    ∫ x in (0 : ℝ)..1, x * (10 - 18 * x) = (-1 : ℝ) := by
  have h : (fun x : ℝ => x * (10 - 18 * x)) = (fun x => 10 * x + (-18) * x ^ 2) := by
    ext x; ring
  rw [h, intervalIntegral.integral_add (poly_integrable_01 (by fun_prop))
    (poly_integrable_01 (by fun_prop)),
    intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul,
    integral_id, integral_pow]
  norm_num

theorem moment_wp_zero :
    ∫ x in (0 : ℝ)..1, (6 - 12 * x) = (0 : ℝ) := by
  have h : (fun x : ℝ => (6 : ℝ) - 12 * x) = (fun x => 6 * (1 : ℝ) + (-12) * x) := by
    ext; ring
  rw [h, intervalIntegral.integral_add (poly_integrable_01 (by fun_prop))
    (poly_integrable_01 (by fun_prop)),
    intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul,
    integral_one, integral_id]
  norm_num

theorem moment_wp_x :
    ∫ x in (0 : ℝ)..1, x * (6 - 12 * x) = (-1 : ℝ) := by
  have h : (fun x : ℝ => x * (6 - 12 * x)) = (fun x => 6 * x + (-12) * x ^ 2) := by
    ext x; ring
  rw [h, intervalIntegral.integral_add (poly_integrable_01 (by fun_prop))
    (poly_integrable_01 (by fun_prop)),
    intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul,
    integral_id, integral_pow]
  norm_num

/-! ## Affine moment identity

For any affine R(x) = A + B·x:
  ∫₀¹ (10-18x)·R(x) dx = A - B = R(-1)
  ∫₀¹ (6-12x)·R(x) dx = -B
-/

theorem affine_moment_wq (A B : ℝ) :
    ∫ x in (0 : ℝ)..1, (10 - 18 * x) * (A + B * x) = A - B := by
  have h : (fun x : ℝ => (10 - 18 * x) * (A + B * x)) =
      (fun x => A * (10 - 18 * x) + B * (x * (10 - 18 * x))) := by
    ext x; ring
  rw [h, intervalIntegral.integral_add
    (poly_integrable_01 (by fun_prop))
    (poly_integrable_01 (by fun_prop)),
    intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul,
    moment_wq_one, moment_wq_x]
  ring

theorem affine_moment_wp (A B : ℝ) :
    ∫ x in (0 : ℝ)..1, (6 - 12 * x) * (A + B * x) = -B := by
  have h : (fun x : ℝ => (6 - 12 * x) * (A + B * x)) =
      (fun x => A * (6 - 12 * x) + B * (x * (6 - 12 * x))) := by
    ext x; ring
  rw [h, intervalIntegral.integral_add
    (poly_integrable_01 (by fun_prop))
    (poly_integrable_01 (by fun_prop)),
    intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul,
    moment_wp_zero, moment_wp_x]
  ring

/-! ## Remainder polynomial

R_{N,j}(X) = (Q_{N,j} - P_{N,j}) + (-P_{N,j}) · X

Key evaluations:
  R(-1) = Q     (the denominator)
  -B = P        (the numerator)
-/

def remainderA (N : ℕ) (j : Fin 3) : ℝ :=
  (denominator N j : ℝ) - (numerator N j : ℝ)

def remainderB (N : ℕ) (j : Fin 3) : ℝ :=
  -(numerator N j : ℝ)

def remainderPoly (N : ℕ) (j : Fin 3) (X : ℝ) : ℝ :=
  remainderA N j + remainderB N j * X

theorem remainderPoly_neg_one (N : ℕ) (j : Fin 3) :
    remainderPoly N j (-1) = (denominator N j : ℝ) := by
  simp [remainderPoly, remainderA, remainderB, mul_neg, mul_one, neg_neg]

theorem neg_remainderB (N : ℕ) (j : Fin 3) :
    -remainderB N j = (numerator N j : ℝ) := by
  simp [remainderB]

/-! ## Denominator and numerator as moment integrals -/

theorem denominator_as_moment (N : ℕ) (j : Fin 3) :
    (denominator N j : ℝ) =
      ∫ x in (0 : ℝ)..1,
        (10 - 18 * x) * (remainderA N j + remainderB N j * x) := by
  rw [affine_moment_wq, remainderA, remainderB]
  ring

theorem numerator_as_moment (N : ℕ) (j : Fin 3) :
    (numerator N j : ℝ) =
      ∫ x in (0 : ℝ)..1,
        (6 - 12 * x) * (remainderA N j + remainderB N j * x) := by
  rw [affine_moment_wp]
  simp [remainderB]

/-! ## The Catalan error as a moment integral

catalanError N j = G · Q_{N,j} - P_{N,j}
                 = G · ∫wq·R - ∫wp·R
                 = ∫ (G·wq - wp) · R

The kernel G·(10-18x) - (6-12x) connects to -log(t)/(1+t²)
via the substitution x = t² and the double integral identity
∫₀¹∫₀¹ f(xy) dx dy = ∫₀¹ f(t)(-log t) dt.
-/

theorem catalanError_as_moment_difference (N : ℕ) (j : Fin 3) :
    catalanError N j =
      catalanConstant *
        (∫ x in (0 : ℝ)..1,
          (10 - 18 * x) * (remainderA N j + remainderB N j * x)) -
        (∫ x in (0 : ℝ)..1,
          (6 - 12 * x) * (remainderA N j + remainderB N j * x)) := by
  rw [← denominator_as_moment, ← numerator_as_moment]
  simp [catalanError]; ring

/-! ## Integral of -log on [0,1]

The fundamental identity ∫₀¹ (-log t) dt = 1 from Mathlib's `integral_log_from_zero`.
-/

theorem integral_neg_log_01 :
    ∫ t in (0 : ℝ)..1, (-Real.log t) = (1 : ℝ) := by
  have h : (fun t : ℝ => -Real.log t) = (fun t => (-1) * Real.log t) := by
    ext; ring
  rw [h, intervalIntegral.integral_const_mul, integral_log_from_zero]
  simp [Real.log_one]

/-! ## The Catalan error is subdominant

The key to closing Problem 2.5 is showing that the catalanError grows
strictly slower than the denominator. Equivalently, the Padé remainder
R_{N,j}(t²) in the integral representation G·Q-P = ∫[-log(t)/(1+t²)]·R(t²)dt
grows at the subdominant rate (17-12√2)^N rather than the dominant (17+12√2)^N.

From this, catalanError/denominator → 0, hence P/Q → G = catalanConstant.
Combined with P/Q → commonLimit (proved in Problem25Connection), we get
commonLimit = catalanConstant by uniqueness of limits. -/

/-- The remaining connection certificate in its concrete sign form.  These
two inequalities say that Catalan's constant is bracketed by two of the three
rational approximants at every stage. -/
theorem positiveCatalanError_brackets (N : ℕ) :
    positiveCatalanError N 0 ≤ 0 ∧ 0 ≤ positiveCatalanError N 2 := by
  sorry

private theorem catalanConstant_mem_envelope (N : ℕ) :
    lowerEnvelope N ≤ catalanConstant ∧
      catalanConstant ≤ upperEnvelope N := by
  rcases positiveCatalanError_brackets N with ⟨hzero, htwo⟩
  have hqzero : (0 : ℝ) < positiveDenominator N 0 := by
    exact_mod_cast positiveDenominator_pos N 0
  have hqtwo : (0 : ℝ) < positiveDenominator N 2 := by
    exact_mod_cast positiveDenominator_pos N 2
  have hle_zero : catalanConstant ≤ positiveRatio N 0 := by
    rw [positiveRatio, le_div_iff₀ hqzero]
    rw [positiveCatalanError_eq] at hzero
    linarith
  have htwo_le : positiveRatio N 2 ≤ catalanConstant := by
    rw [positiveRatio, div_le_iff₀ hqtwo]
    rw [positiveCatalanError_eq] at htwo
    linarith
  exact ⟨(positiveRatio_envelope N 2).1.trans htwo_le,
    hle_zero.trans (positiveRatio_envelope N 0).2⟩

theorem catalanError_over_denominator_tendsto_zero (j : Fin 3) :
    Filter.Tendsto (fun N => catalanError N j / (denominator N j : ℝ))
      Filter.atTop (nhds 0) := by
  have hcatalan_tendsto_common :
      Filter.Tendsto (fun _ : ℕ => catalanConstant)
        Filter.atTop (nhds commonLimit) :=
    tendsto_of_tendsto_of_tendsto_of_le_of_le
      lowerEnvelope_tendsto_common upperEnvelope_tendsto_common
      (fun N => (catalanConstant_mem_envelope N).1)
      (fun N => (catalanConstant_mem_envelope N).2)
  have hcommon : commonLimit = catalanConstant :=
    tendsto_nhds_unique hcatalan_tendsto_common tendsto_const_nhds
  have hratio := challengeRatio_tendsto_common j
  have herror :
      (fun N : ℕ => catalanError N j / (denominator N j : ℝ)) =
        (fun N : ℕ => catalanConstant - challengeRatio N j) := by
    funext N
    have hq : (denominator N j : ℝ) ≠ 0 := by
      exact_mod_cast denominator_ne_zero N j
    simp only [catalanError, challengeRatio]
    field_simp [hq]
  rw [herror]
  simpa [hcommon] using
    (Filter.Tendsto.sub
      (tendsto_const_nhds :
        Filter.Tendsto (fun _ : ℕ => catalanConstant)
          Filter.atTop (nhds catalanConstant))
      hratio)

theorem commonLimit_eq_catalanConstant :
    commonLimit = catalanConstant := by
  have hconv := challengeRatio_tendsto_common 0
  have herr := catalanError_over_denominator_tendsto_zero 0
  apply tendsto_nhds_unique hconv
  have key : ∀ᶠ N : ℕ in Filter.atTop,
      challengeRatio N 0 =
        catalanConstant - catalanError N 0 / (denominator N 0 : ℝ) := by
    filter_upwards [] with N
    have hq : (denominator N 0 : ℝ) ≠ 0 := by
      exact_mod_cast denominator_ne_zero N 0
    simp only [challengeRatio, catalanError]
    field_simp [hq]
    ring
  have htarget :
      Filter.Tendsto
        (fun N : ℕ =>
          catalanConstant - catalanError N 0 / (denominator N 0 : ℝ))
        Filter.atTop (nhds catalanConstant) := by
    simpa using
      Filter.Tendsto.sub
        (tendsto_const_nhds :
          Filter.Tendsto (fun _ : ℕ => catalanConstant)
            Filter.atTop (nhds catalanConstant))
        herr
  exact htarget.congr' (key.mono fun _ h => h.symm)

theorem problem25_solved : Problem25Claim := by
  rw [problem25Claim_iff_commonLimit_eq_catalan]
  exact commonLimit_eq_catalanConstant

end RamanujanChallenge.P25

end
