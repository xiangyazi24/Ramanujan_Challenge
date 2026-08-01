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
  simp [remainderPoly, remainderA, remainderB, mul_neg, mul_one, neg_neg, sub_neg_eq_add,
    add_sub_cancel]

theorem neg_remainderB (N : ℕ) (j : Fin 3) :
    -remainderB N j = (numerator N j : ℝ) := by
  simp [remainderB]

/-! ## Denominator and numerator as moment integrals -/

theorem denominator_as_moment (N : ℕ) (j : Fin 3) :
    (denominator N j : ℝ) =
      ∫ x in (0 : ℝ)..1,
        (10 - 18 * x) * (remainderA N j + remainderB N j * x) := by
  rw [affine_moment_wq]
  simp [remainderA, remainderB]; ring

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

end RamanujanChallenge.P25

end
