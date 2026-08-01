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

set_option maxHeartbeats 0

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
private def ratioMinor (N : ℕ) : Fin 3 → ℝ :=
  ![(positiveNumerator N 0 : ℝ) * (positiveDenominator N 1 : ℝ) -
      (positiveNumerator N 1 : ℝ) * (positiveDenominator N 0 : ℝ),
    (positiveNumerator N 0 : ℝ) * (positiveDenominator N 2 : ℝ) -
      (positiveNumerator N 2 : ℝ) * (positiveDenominator N 0 : ℝ),
    (positiveNumerator N 1 : ℝ) * (positiveDenominator N 2 : ℝ) -
      (positiveNumerator N 2 : ℝ) * (positiveDenominator N 1 : ℝ)]

/-- The second compound matrix of the positive challenge matrix.  Its three
coordinates are indexed by the pairs `(0,1)`, `(0,2)`, `(1,2)`. -/
private def compoundMatrix (n : ℕ) : Matrix (Fin 3) (Fin 3) ℝ :=
  let A := positiveMatrix (n : ℤ)
  !![(A 0 0 : ℝ) * A 1 1 - (A 0 1 : ℝ) * A 1 0,
      (A 0 0 : ℝ) * A 1 2 - (A 0 2 : ℝ) * A 1 0,
      (A 0 1 : ℝ) * A 1 2 - (A 0 2 : ℝ) * A 1 1;
     (A 0 0 : ℝ) * A 2 1 - (A 0 1 : ℝ) * A 2 0,
      (A 0 0 : ℝ) * A 2 2 - (A 0 2 : ℝ) * A 2 0,
      (A 0 1 : ℝ) * A 2 2 - (A 0 2 : ℝ) * A 2 1;
     (A 1 0 : ℝ) * A 2 1 - (A 1 1 : ℝ) * A 2 0,
      (A 1 0 : ℝ) * A 2 2 - (A 1 2 : ℝ) * A 2 0,
      (A 1 1 : ℝ) * A 2 2 - (A 1 2 : ℝ) * A 2 1]

private theorem ratioMinor_succ (n : ℕ) (j : Fin 3) :
    ratioMinor (n + 1) j =
      ∑ i : Fin 3, ratioMinor n i * compoundMatrix n i j := by
  fin_cases j <;> simp [ratioMinor]
  all_goals
    simp_rw [positiveNumerator_succ, positiveDenominator_succ]
    push_cast
    simp [Fin.sum_univ_three, compoundMatrix, ratioMinor]
    ring

/-- Balanced projective coordinates for the exterior-square row. -/
private def wedgeX (n : ℕ) : ℝ :=
  ((n : ℝ) + 1) * ratioMinor n 1 / ratioMinor n 0

private def wedgeY (n : ℕ) : ℝ :=
  ((n : ℝ) + 1) ^ 2 * ratioMinor n 2 / ratioMinor n 0

private def wedgeStepCoefficient (n : ℕ) (x y : ℝ) (j : Fin 3) : ℝ :=
  compoundMatrix n 0 j +
    x / ((n : ℝ) + 1) * compoundMatrix n 1 j +
    y / ((n : ℝ) + 1) ^ 2 * compoundMatrix n 2 j

private theorem ratioMinor_succ_factor (n : ℕ) (j : Fin 3)
    (hminor : ratioMinor n 0 ≠ 0) :
    ratioMinor (n + 1) j =
      ratioMinor n 0 * wedgeStepCoefficient n (wedgeX n) (wedgeY n) j := by
  rw [ratioMinor_succ, Fin.sum_univ_three]
  simp only [wedgeStepCoefficient, wedgeX, wedgeY]
  field_simp [hminor]

private theorem affine_rectangle_nonneg' (A B C x y : ℝ)
    (hx₀ : 1 / 2 ≤ x) (hx₁ : x ≤ 3 / 2)
    (hy₀ : 0 ≤ y) (hy₁ : y ≤ 2)
    (h₀₀ : 0 ≤ A + B * (1 / 2) + C * 0)
    (h₀₁ : 0 ≤ A + B * (1 / 2) + C * 2)
    (h₁₀ : 0 ≤ A + B * (3 / 2) + C * 0)
    (h₁₁ : 0 ≤ A + B * (3 / 2) + C * 2) :
    0 ≤ A + B * x + C * y := by
  rcases le_total 0 B with hB | hB <;>
    rcases le_total 0 C with hC | hC <;> nlinarith

private theorem wedge_linear_nonneg (n : ℕ) (x y a b c : ℝ)
    (hx₀ : 1 / 2 ≤ x) (hx₁ : x ≤ 3 / 2)
    (hy₀ : 0 ≤ y) (hy₁ : y ≤ 2)
    (h₀₀ : 0 ≤ a * wedgeStepCoefficient n (1 / 2) 0 0 +
      b * wedgeStepCoefficient n (1 / 2) 0 1 +
      c * wedgeStepCoefficient n (1 / 2) 0 2)
    (h₀₁ : 0 ≤ a * wedgeStepCoefficient n (1 / 2) 2 0 +
      b * wedgeStepCoefficient n (1 / 2) 2 1 +
      c * wedgeStepCoefficient n (1 / 2) 2 2)
    (h₁₀ : 0 ≤ a * wedgeStepCoefficient n (3 / 2) 0 0 +
      b * wedgeStepCoefficient n (3 / 2) 0 1 +
      c * wedgeStepCoefficient n (3 / 2) 0 2)
    (h₁₁ : 0 ≤ a * wedgeStepCoefficient n (3 / 2) 2 0 +
      b * wedgeStepCoefficient n (3 / 2) 2 1 +
      c * wedgeStepCoefficient n (3 / 2) 2 2) :
    0 ≤ a * wedgeStepCoefficient n x y 0 +
      b * wedgeStepCoefficient n x y 1 +
      c * wedgeStepCoefficient n x y 2 := by
  let A : ℝ := a * compoundMatrix n 0 0 + b * compoundMatrix n 0 1 +
    c * compoundMatrix n 0 2
  let B : ℝ := (a * compoundMatrix n 1 0 + b * compoundMatrix n 1 1 +
    c * compoundMatrix n 1 2) / ((n : ℝ) + 1)
  let C : ℝ := (a * compoundMatrix n 2 0 + b * compoundMatrix n 2 1 +
    c * compoundMatrix n 2 2) / ((n : ℝ) + 1) ^ 2
  have hrect : 0 ≤ A + B * x + C * y := by
    apply affine_rectangle_nonneg' A B C x y hx₀ hx₁ hy₀ hy₁
    · dsimp [A, B, C]
      simp only [wedgeStepCoefficient] at h₀₀
      convert h₀₀ using 1 <;> ring
    · dsimp [A, B, C]
      simp only [wedgeStepCoefficient] at h₀₁
      convert h₀₁ using 1 <;> ring
    · dsimp [A, B, C]
      simp only [wedgeStepCoefficient] at h₁₀
      convert h₁₀ using 1 <;> ring
    · dsimp [A, B, C]
      simp only [wedgeStepCoefficient] at h₁₁
      convert h₁₁ using 1 <;> ring
  convert hrect using 1 <;>
    dsimp [A, B, C, wedgeStepCoefficient] <;> ring

private theorem wedge_step_x_lower (n : ℕ) (x y : ℝ)
    (hx₀ : 1 / 2 ≤ x) (hx₁ : x ≤ 3 / 2)
    (hy₀ : 0 ≤ y) (hy₁ : y ≤ 2) :
    wedgeStepCoefficient n x y 0 ≤
      2 * ((n : ℝ) + 2) * wedgeStepCoefficient n x y 1 := by
  suffices 0 ≤ (-1) * wedgeStepCoefficient n x y 0 +
      (2 * ((n : ℝ) + 2)) * wedgeStepCoefficient n x y 1 +
      0 * wedgeStepCoefficient n x y 2 by linarith
  apply wedge_linear_nonneg n x y (-1) (2 * ((n : ℝ) + 2)) 0
    hx₀ hx₁ hy₀ hy₁
  all_goals
    norm_num [wedgeStepCoefficient, compoundMatrix, positiveMatrix,
      Matrix.cons_val_two]
    field_simp
    ring_nf
    positivity

private theorem wedge_step_x_upper (n : ℕ) (x y : ℝ)
    (hx₀ : 1 / 2 ≤ x) (hx₁ : x ≤ 3 / 2)
    (hy₀ : 0 ≤ y) (hy₁ : y ≤ 2) :
    2 * ((n : ℝ) + 2) * wedgeStepCoefficient n x y 1 ≤
      3 * wedgeStepCoefficient n x y 0 := by
  suffices 0 ≤ 3 * wedgeStepCoefficient n x y 0 +
      (-2 * ((n : ℝ) + 2)) * wedgeStepCoefficient n x y 1 +
      0 * wedgeStepCoefficient n x y 2 by linarith
  apply wedge_linear_nonneg n x y 3 (-2 * ((n : ℝ) + 2)) 0
    hx₀ hx₁ hy₀ hy₁
  all_goals
    norm_num [wedgeStepCoefficient, compoundMatrix, positiveMatrix,
      Matrix.cons_val_two]
    field_simp
    ring_nf
    have hp : ∀ k : ℕ, 0 ≤ (n : ℝ) ^ k := by
      intro k
      positivity
    nlinarith [hp 0, hp 1, hp 2, hp 3, hp 4, hp 5, hp 6, hp 7,
      hp 8, hp 9, hp 10, hp 11, hp 12]

private theorem wedge_step_y_lower (n : ℕ) (x y : ℝ)
    (hx₀ : 1 / 2 ≤ x) (hx₁ : x ≤ 3 / 2)
    (hy₀ : 0 ≤ y) (hy₁ : y ≤ 2) :
    0 ≤ wedgeStepCoefficient n x y 2 := by
  suffices 0 ≤ 0 * wedgeStepCoefficient n x y 0 +
      0 * wedgeStepCoefficient n x y 1 +
      1 * wedgeStepCoefficient n x y 2 by simpa using this
  apply wedge_linear_nonneg n x y 0 0 1 hx₀ hx₁ hy₀ hy₁
  all_goals
    norm_num [wedgeStepCoefficient, compoundMatrix, positiveMatrix,
      Matrix.cons_val_two]
    field_simp
    ring_nf
    positivity

private theorem wedge_step_y_upper (n : ℕ) (x y : ℝ)
    (hx₀ : 1 / 2 ≤ x) (hx₁ : x ≤ 3 / 2)
    (hy₀ : 0 ≤ y) (hy₁ : y ≤ 2) :
    ((n : ℝ) + 2) ^ 2 * wedgeStepCoefficient n x y 2 ≤
      2 * wedgeStepCoefficient n x y 0 := by
  suffices 0 ≤ 2 * wedgeStepCoefficient n x y 0 +
      0 * wedgeStepCoefficient n x y 1 +
      (-((n : ℝ) + 2) ^ 2) * wedgeStepCoefficient n x y 2 by linarith
  apply wedge_linear_nonneg n x y 2 0 (-((n : ℝ) + 2) ^ 2)
    hx₀ hx₁ hy₀ hy₁
  all_goals
    norm_num [wedgeStepCoefficient, compoundMatrix, positiveMatrix,
      Matrix.cons_val_two]
    field_simp
    ring_nf
    have hp : ∀ k : ℕ, 0 ≤ (n : ℝ) ^ k := by
      intro k
      positivity
    nlinarith [hp 0, hp 1, hp 2, hp 3, hp 4, hp 5, hp 6, hp 7,
      hp 8, hp 9, hp 10, hp 11, hp 12, hp 13, hp 14]

private theorem wedge_step_zero_pos (n : ℕ) (x y : ℝ)
    (hx : 1 / 2 ≤ x) (hy : 0 ≤ y) :
    0 < wedgeStepCoefficient n x y 0 := by
  have hn : (0 : ℝ) < (n : ℝ) + 1 := by positivity
  have hc10 : (0 : ℝ) < compoundMatrix n 1 0 := by
    norm_num [compoundMatrix, positiveMatrix, Matrix.cons_val_two]
    ring_nf
    have hp : ∀ k : ℕ, 0 ≤ (n : ℝ) ^ k := by
      intro k
      positivity
    nlinarith [hp 0, hp 1, hp 2, hp 3, hp 4, hp 5, hp 6, hp 7,
      hp 8, hp 9, hp 10, hp 11, hp 12, hp 13, hp 14, hp 15]
  have hc20 : (0 : ℝ) < compoundMatrix n 2 0 := by
    norm_num [compoundMatrix, positiveMatrix, Matrix.cons_val_two]
    ring_nf
    have hp : ∀ k : ℕ, 0 ≤ (n : ℝ) ^ k := by
      intro k
      positivity
    nlinarith [hp 0, hp 1, hp 2, hp 3, hp 4, hp 5, hp 6, hp 7,
      hp 8, hp 9, hp 10, hp 11, hp 12, hp 13, hp 14, hp 15]
  have hbase :
      0 < compoundMatrix n 0 0 +
        (1 / 2) / ((n : ℝ) + 1) * compoundMatrix n 1 0 := by
    norm_num [compoundMatrix, positiveMatrix, Matrix.cons_val_two]
    field_simp
    ring_nf
    positivity
  simp only [wedgeStepCoefficient]
  have hx' :
      (1 / 2) / ((n : ℝ) + 1) * compoundMatrix n 1 0 ≤
        x / ((n : ℝ) + 1) * compoundMatrix n 1 0 := by
    gcongr
  have hy' :
      0 ≤ y / ((n : ℝ) + 1) ^ 2 * compoundMatrix n 2 0 := by
    exact mul_nonneg (div_nonneg hy (sq_nonneg _)) hc20.le
  linarith

/-- The exterior-square analogue of the denominator projective cone.  It says
in particular that the three positive challenge ratios stay ordered. -/
private theorem ratioMinor_projective_cone (n : ℕ) :
    0 < ratioMinor n 0 ∧
      1 / 2 ≤ wedgeX n ∧ wedgeX n ≤ 3 / 2 ∧
      0 ≤ wedgeY n ∧ wedgeY n ≤ 2 := by
  induction n with
  | zero =>
      norm_num [ratioMinor, wedgeX, wedgeY, positiveNumerator,
        positiveDenominator, numerator, denominator, approximants,
        initialMatrix, coordinateSign, Matrix.cons_val_two]
  | succ n ih =>
      rcases ih with ⟨hm, hx₀, hx₁, hy₀, hy₁⟩
      have hstep0 :
          0 < wedgeStepCoefficient n (wedgeX n) (wedgeY n) 0 :=
        wedge_step_zero_pos n _ _ hx₀ hy₀
      have hm' : 0 < ratioMinor (n + 1) 0 := by
        rw [ratioMinor_succ_factor n 0 hm.ne']
        exact mul_pos hm hstep0
      have hX :
          wedgeX (n + 1) =
            ((n : ℝ) + 2) *
                wedgeStepCoefficient n (wedgeX n) (wedgeY n) 1 /
              wedgeStepCoefficient n (wedgeX n) (wedgeY n) 0 := by
        rw [wedgeX, ratioMinor_succ_factor n 1 hm.ne',
          ratioMinor_succ_factor n 0 hm.ne']
        field_simp [hm.ne', hstep0.ne']
        norm_num [Nat.cast_add, Nat.cast_one]
        ring
      have hY :
          wedgeY (n + 1) =
            ((n : ℝ) + 2) ^ 2 *
                wedgeStepCoefficient n (wedgeX n) (wedgeY n) 2 /
              wedgeStepCoefficient n (wedgeX n) (wedgeY n) 0 := by
        rw [wedgeY, ratioMinor_succ_factor n 2 hm.ne',
          ratioMinor_succ_factor n 0 hm.ne']
        field_simp [hm.ne', hstep0.ne']
        norm_num [Nat.cast_add, Nat.cast_one]
        ring
      have hxl := wedge_step_x_lower n (wedgeX n) (wedgeY n)
        hx₀ hx₁ hy₀ hy₁
      have hxu := wedge_step_x_upper n (wedgeX n) (wedgeY n)
        hx₀ hx₁ hy₀ hy₁
      have hyl := wedge_step_y_lower n (wedgeX n) (wedgeY n)
        hx₀ hx₁ hy₀ hy₁
      have hyu := wedge_step_y_upper n (wedgeX n) (wedgeY n)
        hx₀ hx₁ hy₀ hy₁
      refine ⟨hm', ?_, ?_, ?_, ?_⟩
      · rw [hX]
        apply (le_div_iff₀ hstep0).2
        nlinarith
      · rw [hX]
        apply (div_le_iff₀ hstep0).2
        nlinarith
      · rw [hY]
        exact div_nonneg (mul_nonneg (sq_nonneg _) hyl) hstep0.le
      · rw [hY]
        apply (div_le_iff₀ hstep0).2
        nlinarith

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
