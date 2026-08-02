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
import RamanujanChallenge.Problem25DualAlgebra
import RamanujanChallenge.Problem25DualInitialMoment
import RamanujanChallenge.Problem25DualCertificateDivergence
import RamanujanChallenge.Problem25DualCertificateSparseRow0
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

/-! ## A fast majorant for the Miller solution

The companion proof below only needs the first coordinate of `dualVector`.
The coarse integrability majorant from `Problem25DualIntegral` gives a factor
`4⁻ⁿ`.  Keeping the two elementary beta factors

`p²(1-p²) ≤ 1/4` and `q²(1-q²) ≤ 1/4`

gives the stronger `64⁻ⁿ` bound needed to recognize the recessive solution. -/

private theorem miller_cube_ae_bounds :
    ∀ᵐ x : ℝ × (ℝ × ℝ) ∂cubeMeasure,
      0 < x.1 ∧ x.1 ≤ 1 ∧
      0 < x.2.1 ∧ x.2.1 ≤ 1 ∧
      0 < x.2.2 ∧ x.2.2 ≤ 1 := by
  have hmem : ∀ᵐ x : ℝ × (ℝ × ℝ) ∂cubeMeasure,
      x ∈ Ioc (0 : ℝ) 1 ×ˢ (Ioc (0 : ℝ) 1 ×ˢ Ioc (0 : ℝ) 1) := by
    rw [Measure.ae_prod_mem_iff_ae_ae_mem
      (measurableSet_Ioc.prod (measurableSet_Ioc.prod measurableSet_Ioc))]
    filter_upwards [unit_ae_bounds] with p hp
    have hqv : ∀ᵐ y : ℝ × ℝ ∂unitMeasure.prod unitMeasure,
        y ∈ Ioc (0 : ℝ) 1 ×ˢ Ioc (0 : ℝ) 1 := by
      rw [Measure.ae_prod_mem_iff_ae_ae_mem
        (measurableSet_Ioc.prod measurableSet_Ioc)]
      filter_upwards [unit_ae_bounds] with q hq
      filter_upwards [unit_ae_bounds] with v hv
      exact ⟨hq, hv⟩
    filter_upwards [hqv] with y hy
    exact ⟨hp, hy⟩
  filter_upwards [hmem] with x hx
  exact ⟨hx.1.1, hx.1.2, hx.2.1.1, hx.2.1.2, hx.2.2.1, hx.2.2.2⟩

private def millerBoundedIntegrand
    (n : ℕ) (x : ℝ × (ℝ × ℝ)) : ℝ :=
  let p := x.1
  let q := x.2.1
  let v := x.2.2
  (512 / 64 ^ n : ℝ) *
    p * (1 - q ^ 2) *
    (4 * p ^ 2 * (1 - p ^ 2)) ^ n *
    (4 * q ^ 2 * (1 - q ^ 2)) ^ n *
    (p * q / dualD p q v) ^ 5 *
    (2 * v / dualD p q v) ^ (2 * n + 3) *
    (dualD p q v / 4) ^ 4

private theorem rawMomentIntegrand_ae_eq_millerBounded (n : ℕ) :
    rawMomentIntegrand n 0 0 0 0 =ᵐ[cubeMeasure]
      millerBoundedIntegrand n := by
  filter_upwards [miller_cube_ae_bounds] with x hx
  rcases hx with ⟨hp, hp1, hq, hq1, hv, hv1⟩
  have hD : dualD x.1 x.2.1 x.2.2 ≠ 0 := by
    apply ne_of_gt
    dsimp [dualD]
    positivity
  have hScalar :
      (512 / 64 ^ n : ℝ) * 4 ^ n * 4 ^ n * 2 ^ (2 * n + 3) *
          (1 / 4) ^ 4 = 16 := by
    rw [show (64 : ℝ) ^ n = 4 ^ n * 4 ^ n * 4 ^ n by
      rw [← mul_pow, ← mul_pow]
      norm_num]
    rw [show (2 : ℝ) ^ (2 * n + 3) = 8 * 4 ^ n by
      rw [pow_add, pow_mul]
      norm_num
      ring]
    field_simp
    norm_num
  have hDen :
      (dualD x.1 x.2.1 x.2.2 ^ 5)⁻¹ *
            (dualD x.1 x.2.1 x.2.2 ^ (2 * n + 3))⁻¹ *
            dualD x.1 x.2.1 x.2.2 ^ 4 =
          (dualD x.1 x.2.1 x.2.2 ^ (2 * n + 4))⁻¹ := by
    field_simp [hD]
    rw [show 2 * n + 4 = 1 + (2 * n + 3) by omega, pow_add]
    simp
  have hpExp :
      x.1 ^ (2 * n + 6) = x.1 * x.1 ^ (2 * n) * x.1 ^ 5 := by
    rw [show 2 * n + 6 = 1 + 2 * n + 5 by omega, pow_add, pow_add]
    ring
  have hqExp :
      x.2.1 ^ (2 * n + 5) = x.2.1 ^ (2 * n) * x.2.1 ^ 5 := by
    rw [show 2 * n + 5 = 2 * n + 5 by rfl, pow_add]
  have hqRem :
      (1 - x.2.1 ^ 2) ^ (n + 1) =
        (1 - x.2.1 ^ 2) * (1 - x.2.1 ^ 2) ^ n := by
    rw [pow_succ]
    ring
  have hpBeta :
      (4 * x.1 ^ 2 * (1 - x.1 ^ 2)) ^ n =
        4 ^ n * x.1 ^ (2 * n) * (1 - x.1 ^ 2) ^ n := by
    rw [mul_pow, mul_pow, ← pow_mul]
  have hqBeta :
      (4 * x.2.1 ^ 2 * (1 - x.2.1 ^ 2)) ^ n =
        4 ^ n * x.2.1 ^ (2 * n) * (1 - x.2.1 ^ 2) ^ n := by
    rw [mul_pow, mul_pow, ← pow_mul]
  dsimp [rawMomentIntegrand, millerBoundedIntegrand]
  rw [hpExp, hqExp, hqRem, hpBeta, hqBeta]
  simp only [div_eq_mul_inv, mul_pow]
  have hScalar' :
      512 * (64 ^ n)⁻¹ * 4 ^ n * 4 ^ n * 2 ^ (2 * n + 3) *
          4⁻¹ ^ 4 = (16 : ℝ) := by
    simpa only [div_eq_mul_inv, one_div, one_mul] using hScalar
  have hDen' :
      (dualD x.1 x.2.1 x.2.2)⁻¹ ^ 5 *
            (dualD x.1 x.2.1 x.2.2)⁻¹ ^ (2 * n + 3) *
            dualD x.1 x.2.1 x.2.2 ^ 4 =
          (dualD x.1 x.2.1 x.2.2 ^ (2 * n + 4))⁻¹ := by
    simpa only [inv_pow] using hDen
  calc
    16 * (x.1 * x.1 ^ (2 * n) * x.1 ^ 5) *
            (x.2.1 ^ (2 * n) * x.2.1 ^ 5) *
            (1 - x.1 ^ 2) ^ n *
            ((1 - x.2.1 ^ 2) * (1 - x.2.1 ^ 2) ^ n) *
            x.2.2 ^ (2 * n + 3) *
            (dualD x.1 x.2.1 x.2.2 ^ (2 * n + 4))⁻¹ =
        (512 * (64 ^ n)⁻¹ * 4 ^ n * 4 ^ n *
              2 ^ (2 * n + 3) * 4⁻¹ ^ 4) *
            (x.1 * x.1 ^ (2 * n) * x.1 ^ 5) *
            (x.2.1 ^ (2 * n) * x.2.1 ^ 5) *
            (1 - x.1 ^ 2) ^ n *
            ((1 - x.2.1 ^ 2) * (1 - x.2.1 ^ 2) ^ n) *
            x.2.2 ^ (2 * n + 3) *
            ((dualD x.1 x.2.1 x.2.2)⁻¹ ^ 5 *
              (dualD x.1 x.2.1 x.2.2)⁻¹ ^ (2 * n + 3) *
              dualD x.1 x.2.1 x.2.2 ^ 4) := by
          rw [hScalar', hDen']
    _ = _ := by ring

private theorem millerBoundedIntegrand_le (n : ℕ) :
    ∀ᵐ x ∂cubeMeasure,
      millerBoundedIntegrand n x ≤ (512 / 64 ^ n : ℝ) := by
  filter_upwards [miller_cube_ae_bounds] with x hx
  rcases hx with ⟨hp, hp1, hq, hq1, hv, hv1⟩
  have hp0 : 0 ≤ x.1 := hp.le
  have hq0 : 0 ≤ x.2.1 := hq.le
  have hv0 : 0 ≤ x.2.2 := hv.le
  have hp2 : x.1 ^ 2 ≤ 1 := by nlinarith [sq_nonneg (1 - x.1)]
  have hq2 : x.2.1 ^ 2 ≤ 1 := by nlinarith [sq_nonneg (1 - x.2.1)]
  have hpm0 : 0 ≤ 1 - x.1 ^ 2 := by linarith
  have hqm0 : 0 ≤ 1 - x.2.1 ^ 2 := by linarith
  have hDpos : 0 < dualD x.1 x.2.1 x.2.2 := by
    dsimp [dualD]
    positivity
  have hpPair0 : 0 ≤ 4 * x.1 ^ 2 * (1 - x.1 ^ 2) := by positivity
  have hqPair0 : 0 ≤ 4 * x.2.1 ^ 2 * (1 - x.2.1 ^ 2) := by positivity
  have hpPair1 : 4 * x.1 ^ 2 * (1 - x.1 ^ 2) ≤ 1 := by
    nlinarith [sq_nonneg (2 * x.1 ^ 2 - 1)]
  have hqPair1 : 4 * x.2.1 ^ 2 * (1 - x.2.1 ^ 2) ≤ 1 := by
    nlinarith [sq_nonneg (2 * x.2.1 ^ 2 - 1)]
  have hpqD : 0 ≤ x.1 * x.2.1 / dualD x.1 x.2.1 x.2.2 :=
    div_nonneg (mul_nonneg hp0 hq0) hDpos.le
  have hpqD1 : x.1 * x.2.1 / dualD x.1 x.2.1 x.2.2 ≤ 1 := by
    rw [div_le_one hDpos]
    dsimp [dualD]
    nlinarith [mul_nonneg (mul_nonneg hp0 hq0) (sq_nonneg x.2.2)]
  have hvD : 0 ≤ 2 * x.2.2 / dualD x.1 x.2.1 x.2.2 :=
    div_nonneg (mul_nonneg (by norm_num) hv0) hDpos.le
  have hvD1 : 2 * x.2.2 / dualD x.1 x.2.1 x.2.2 ≤ 1 := by
    rw [div_le_one hDpos]
    dsimp [dualD]
    nlinarith [mul_nonneg (mul_nonneg hp0 hq0)
      (by positivity : 0 ≤ 1 + x.2.2 ^ 2)]
  have hD0 : 0 ≤ dualD x.1 x.2.1 x.2.2 / 4 :=
    div_nonneg hDpos.le (by norm_num)
  have hD1 : dualD x.1 x.2.1 x.2.2 / 4 ≤ 1 := by
    rw [div_le_one (by norm_num : (0 : ℝ) < 4)]
    have hpq : x.1 * x.2.1 ≤ 1 := by
      nlinarith [mul_le_mul hp1 hq1 hq0 (by norm_num : (0 : ℝ) ≤ 1)]
    have hv2 : x.2.2 ^ 2 ≤ 1 := by
      nlinarith [sq_nonneg (1 - x.2.2)]
    have hterm : x.1 * x.2.1 * (1 + x.2.2 ^ 2) ≤ 2 := by
      calc
        x.1 * x.2.1 * (1 + x.2.2 ^ 2) ≤ 1 * 2 := by
          apply mul_le_mul hpq (by linarith) (by positivity) (by norm_num)
        _ = 2 := by norm_num
    dsimp [dualD]
    nlinarith
  have hpPow : x.1 ≤ 1 := hp1
  have hqm : 1 - x.2.1 ^ 2 ≤ 1 := by nlinarith [sq_nonneg x.2.1]
  have hpBeta :
      (4 * x.1 ^ 2 * (1 - x.1 ^ 2)) ^ n ≤ 1 :=
    pow_le_one₀ hpPair0 hpPair1
  have hqBeta :
      (4 * x.2.1 ^ 2 * (1 - x.2.1 ^ 2)) ^ n ≤ 1 :=
    pow_le_one₀ hqPair0 hqPair1
  have hr1 : (x.1 * x.2.1 / dualD x.1 x.2.1 x.2.2) ^ 5 ≤ 1 :=
    pow_le_one₀ hpqD hpqD1
  have hr2 : (2 * x.2.2 / dualD x.1 x.2.1 x.2.2) ^ (2 * n + 3) ≤ 1 :=
    pow_le_one₀ hvD hvD1
  have hr3 : (dualD x.1 x.2.1 x.2.2 / 4) ^ 4 ≤ 1 :=
    pow_le_one₀ hD0 hD1
  have h1 := mul_le_one₀ hpPow hqm0 hqm
  have h2 := mul_le_one₀ h1 (pow_nonneg hpPair0 _) hpBeta
  have h3 := mul_le_one₀ h2 (pow_nonneg hqPair0 _) hqBeta
  have h4 := mul_le_one₀ h3 (pow_nonneg hpqD _) hr1
  have h5 := mul_le_one₀ h4 (pow_nonneg hvD _) hr2
  have h6 := mul_le_one₀ h5 (pow_nonneg hD0 _) hr3
  have hconst : 0 ≤ (512 / 64 ^ n : ℝ) := by positivity
  dsimp [millerBoundedIntegrand]
  calc
    (512 / 64 ^ n : ℝ) * x.1 * (1 - x.2.1 ^ 2) *
          (4 * x.1 ^ 2 * (1 - x.1 ^ 2)) ^ n *
          (4 * x.2.1 ^ 2 * (1 - x.2.1 ^ 2)) ^ n *
          (x.1 * x.2.1 / dualD x.1 x.2.1 x.2.2) ^ 5 *
          (2 * x.2.2 / dualD x.1 x.2.1 x.2.2) ^ (2 * n + 3) *
          (dualD x.1 x.2.1 x.2.2 / 4) ^ 4 =
        (512 / 64 ^ n : ℝ) *
          (x.1 * (1 - x.2.1 ^ 2) *
            (4 * x.1 ^ 2 * (1 - x.1 ^ 2)) ^ n *
            (4 * x.2.1 ^ 2 * (1 - x.2.1 ^ 2)) ^ n *
            (x.1 * x.2.1 / dualD x.1 x.2.1 x.2.2) ^ 5 *
            (2 * x.2.2 / dualD x.1 x.2.1 x.2.2) ^ (2 * n + 3) *
            (dualD x.1 x.2.1 x.2.2 / 4) ^ 4) := by ring
    _ ≤ (512 / 64 ^ n : ℝ) * 1 :=
      mul_le_mul_of_nonneg_left h6 hconst
    _ = (512 / 64 ^ n : ℝ) := by ring

private theorem dualMoment_zero_le_fast (n : ℕ) :
    dualMoment n 0 0 0 0 ≤ (512 / 64 ^ n : ℝ) := by
  rw [dualMoment, integral_congr_ae (rawMomentIntegrand_ae_eq_millerBounded n)]
  have hconstInt :
      Integrable (fun _ : ℝ × (ℝ × ℝ) => (512 / 64 ^ n : ℝ)) cubeMeasure :=
    integrable_const _
  have hmillerInt : Integrable (millerBoundedIntegrand n) cubeMeasure :=
    (rawMomentIntegrand_integrable n 0 0 0 0 (by omega)).congr
      (rawMomentIntegrand_ae_eq_millerBounded n)
  have hunit : unitMeasure Set.univ = 1 := by
    simp [unitMeasure]
  have hcube : cubeMeasure Set.univ = 1 := by
    change (unitMeasure.prod (unitMeasure.prod unitMeasure)) Set.univ = 1
    calc
      _ = unitMeasure Set.univ *
          (unitMeasure.prod unitMeasure) Set.univ := by
            rw [← Set.univ_prod_univ, MeasureTheory.Measure.prod_prod]
      _ = unitMeasure Set.univ *
          (unitMeasure Set.univ * unitMeasure Set.univ) := by
            rw [← Set.univ_prod_univ, MeasureTheory.Measure.prod_prod]
      _ = 1 := by rw [hunit]; norm_num
  calc
    (∫ x, millerBoundedIntegrand n x ∂cubeMeasure) ≤
        ∫ _x, (512 / 64 ^ n : ℝ) ∂cubeMeasure :=
      integral_mono_ae hmillerInt hconstInt (millerBoundedIntegrand_le n)
    _ = (512 / 64 ^ n : ℝ) := by
      rw [integral_const]
      change ENNReal.toReal (cubeMeasure Set.univ) *
          (512 / 64 ^ n : ℝ) = _
      rw [hcube]
      norm_num

/-! ## Kernel-checked sparse certificate semantics

The imported creative-telescoping recurrence originally used generated
decision certificates for two of its three polynomial rows.  The Miller proof
below uses a sparse polynomial checker instead.  This section connects that
checker to the ordinary polynomial semantics inside Lean's kernel. -/

private def millerSpCoeffEval (n : ℝ) : List ℤ → ℝ
  | [] => 0
  | c :: cs => c + n * millerSpCoeffEval n cs

private theorem millerSpCoeffEval_normalize (n : ℝ) (cs : List ℤ) :
    millerSpCoeffEval n (spCoeffNormalize cs) =
      millerSpCoeffEval n cs := by
  induction cs with
  | nil => rfl
  | cons c cs ih =>
      cases h : spCoeffNormalize cs with
      | nil =>
          rw [h] at ih
          simp only [millerSpCoeffEval] at ih
          by_cases hc : c = 0
          · simp [spCoeffNormalize, h, hc, millerSpCoeffEval, ← ih]
          · simp [spCoeffNormalize, h, hc, millerSpCoeffEval, ← ih]
      | cons d ds =>
          rw [h] at ih
          simp only [millerSpCoeffEval] at ih
          simp only [spCoeffNormalize, h, millerSpCoeffEval]
          rw [ih]

private theorem millerSpCoeffEval_addRaw
    (n : ℝ) (as bs : List ℤ) :
    millerSpCoeffEval n (spCoeffAddRaw as bs) =
      millerSpCoeffEval n as + millerSpCoeffEval n bs := by
  induction as generalizing bs with
  | nil => simp [spCoeffAddRaw, millerSpCoeffEval]
  | cons a as ih =>
      cases bs with
      | nil => simp [spCoeffAddRaw, millerSpCoeffEval]
      | cons b bs =>
          simp [spCoeffAddRaw, millerSpCoeffEval, ih]
          ring

@[simp] private theorem millerSpCoeffEval_add
    (n : ℝ) (as bs : List ℤ) :
    millerSpCoeffEval n (spCoeffAdd as bs) =
      millerSpCoeffEval n as + millerSpCoeffEval n bs := by
  rw [spCoeffAdd, millerSpCoeffEval_normalize,
    millerSpCoeffEval_addRaw]

private theorem millerSpCoeffEval_map_neg
    (n : ℝ) (as : List ℤ) :
    millerSpCoeffEval n (as.map (-·)) =
      -millerSpCoeffEval n as := by
  induction as with
  | nil => simp [millerSpCoeffEval]
  | cons a as ih =>
      simp [millerSpCoeffEval, ih]
      ring

@[simp] private theorem millerSpCoeffEval_neg
    (n : ℝ) (as : List ℤ) :
    millerSpCoeffEval n (spCoeffNeg as) =
      -millerSpCoeffEval n as := by
  rw [spCoeffNeg, millerSpCoeffEval_normalize,
    millerSpCoeffEval_map_neg]

@[simp] private theorem millerSpCoeffEval_scale
    (n : ℝ) (a : ℤ) (bs : List ℤ) :
    millerSpCoeffEval n (spCoeffScale a bs) =
      a * millerSpCoeffEval n bs := by
  induction bs with
  | nil => simp [spCoeffScale, millerSpCoeffEval]
  | cons b bs ih =>
      simp only [spCoeffScale, List.map_cons, millerSpCoeffEval]
      change
        millerSpCoeffEval n (List.map (fun x => a * x) bs) =
          (a : ℝ) * millerSpCoeffEval n bs at ih
      rw [ih]
      push_cast
      ring

private theorem millerSpCoeffEval_mulRaw
    (n : ℝ) (as bs : List ℤ) :
    millerSpCoeffEval n (spCoeffMulRaw as bs) =
      millerSpCoeffEval n as * millerSpCoeffEval n bs := by
  induction as with
  | nil => simp [spCoeffMulRaw, millerSpCoeffEval]
  | cons a as ih =>
      simp [spCoeffMulRaw, millerSpCoeffEval_addRaw,
        millerSpCoeffEval_scale, millerSpCoeffEval, ih]
      ring

@[simp] private theorem millerSpCoeffEval_mul
    (n : ℝ) (as bs : List ℤ) :
    millerSpCoeffEval n (spCoeffMul as bs) =
      millerSpCoeffEval n as * millerSpCoeffEval n bs := by
  rw [spCoeffMul, millerSpCoeffEval_normalize,
    millerSpCoeffEval_mulRaw]

private def millerSparseTermEval
    (n p q v : ℝ) (t : SparseTerm) : ℝ :=
  millerSpCoeffEval n t.nCoeffs *
    p ^ t.pExp * q ^ t.qExp * v ^ t.vExp

private def millerSparseEval
    (P : SparsePoly) (n p q v : ℝ) : ℝ :=
  (P.terms.map (millerSparseTermEval n p q v)).sum

private theorem millerSparseExpCompare_eq_iff
    (a b : SparseTerm) :
    dualCertExpCompare (sparseExp a) (sparseExp b) = .eq ↔
      sparseExp a = sparseExp b := by
  constructor
  · intro h
    by_contra hab
    simp only [dualCertExpCompare, if_neg hab] at h
    split at h <;> simp_all
  · intro h
    simp [dualCertExpCompare, h]

private theorem millerSparseTermEval_add_coeff
    (n p q v : ℝ) (a b : SparseTerm)
    (h : sparseExp a = sparseExp b) :
    millerSparseTermEval n p q v
        { a with nCoeffs := spCoeffAdd a.nCoeffs b.nCoeffs } =
      millerSparseTermEval n p q v a +
        millerSparseTermEval n p q v b := by
  rcases a with ⟨ac, ap, aq, av⟩
  rcases b with ⟨bc, bp, bq, bv⟩
  simp only [sparseExp] at h
  injection h with _ hp hq hv
  subst bp
  subst bq
  subst bv
  simp [millerSparseTermEval]
  ring

private theorem millerSparseEval_insert
    (n p q v : ℝ) (t : SparseTerm) :
    ∀ ts : List SparseTerm,
      ((sparseInsert t ts).map
          (millerSparseTermEval n p q v)).sum =
        millerSparseTermEval n p q v t +
          (ts.map (millerSparseTermEval n p q v)).sum := by
  intro ts
  induction ts with
  | nil =>
      simp only [sparseInsert]
      split_ifs with h
      · simp [millerSparseTermEval, h, millerSpCoeffEval]
      · simp
  | cons u us ih =>
      simp only [sparseInsert]
      cases hcmp :
          dualCertExpCompare (sparseExp t) (sparseExp u) with
      | lt =>
          simp only [hcmp]
          split_ifs with h
          · simp [millerSparseTermEval, h, millerSpCoeffEval]
          · simp
      | eq =>
          simp only [hcmp]
          have hexp : sparseExp t = sparseExp u :=
            (millerSparseExpCompare_eq_iff t u).mp hcmp
          split_ifs with hsum
          · have heval :
                millerSparseTermEval n p q v t +
                    millerSparseTermEval n p q v u = 0 := by
              rw [← millerSparseTermEval_add_coeff
                n p q v t u hexp]
              simp [millerSparseTermEval, hsum,
                millerSpCoeffEval]
            simp only [List.map_cons, List.sum_cons]
            linarith
          · simp only [List.map_cons, List.sum_cons]
            rw [millerSparseTermEval_add_coeff
              n p q v t u hexp]
            ring
      | gt =>
          simp only [hcmp, List.map_cons, List.sum_cons, ih]
          ring

private theorem millerSparseEval_normalize
    (P : SparsePoly) (n p q v : ℝ) :
    millerSparseEval (sparseNormalize P) n p q v =
      millerSparseEval P n p q v := by
  rcases P with ⟨ts⟩
  unfold sparseNormalize millerSparseEval
  induction ts with
  | nil => simp
  | cons t ts ih =>
      simp only [List.foldr_cons, List.map_cons, List.sum_cons]
      rw [millerSparseEval_insert, ih]

@[simp] private theorem millerSparseEval_const
    (z : ℤ) (n p q v : ℝ) :
    millerSparseEval (SparsePoly.const z) n p q v = z := by
  rw [SparsePoly.const, millerSparseEval_normalize]
  by_cases hz : z = 0
  · simp [millerSparseEval, millerSparseTermEval,
      millerSpCoeffEval, hz]
  · simp [millerSparseEval, millerSparseTermEval,
      millerSpCoeffEval, hz]

@[simp] private theorem millerSparseEval_n (n p q v : ℝ) :
    millerSparseEval SparsePoly.n n p q v = n := by
  simp [SparsePoly.n, millerSparseEval, millerSparseTermEval,
    millerSpCoeffEval]

@[simp] private theorem millerSparseEval_p (n p q v : ℝ) :
    millerSparseEval SparsePoly.p n p q v = p := by
  simp [SparsePoly.p, millerSparseEval, millerSparseTermEval,
    millerSpCoeffEval]

@[simp] private theorem millerSparseEval_q (n p q v : ℝ) :
    millerSparseEval SparsePoly.q n p q v = q := by
  simp [SparsePoly.q, millerSparseEval, millerSparseTermEval,
    millerSpCoeffEval]

@[simp] private theorem millerSparseEval_v (n p q v : ℝ) :
    millerSparseEval SparsePoly.v n p q v = v := by
  simp [SparsePoly.v, millerSparseEval, millerSparseTermEval,
    millerSpCoeffEval]

@[simp] private theorem millerSparseEval_ofNat
    (k : ℕ) (n p q v : ℝ) :
    millerSparseEval (OfNat.ofNat k : SparsePoly) n p q v = k := by
  rw [sparse_ofNat_eq, millerSparseEval_const]
  norm_num

@[simp] private theorem millerSparseEval_add
    (P Q : SparsePoly) (n p q v : ℝ) :
    millerSparseEval (P + Q) n p q v =
      millerSparseEval P n p q v +
        millerSparseEval Q n p q v := by
  rw [sparse_add_eq, millerSparseEval_normalize]
  simp [millerSparseEval]

@[simp] private theorem millerSparseEval_neg
    (P : SparsePoly) (n p q v : ℝ) :
    millerSparseEval (-P) n p q v =
      -millerSparseEval P n p q v := by
  rw [sparse_neg_eq, millerSparseEval_normalize]
  rcases P with ⟨ts⟩
  induction ts with
  | nil => simp [millerSparseEval]
  | cons t ts ih =>
      simp only [millerSparseEval, List.map_cons,
        List.sum_cons] at ih ⊢
      simp only [List.map_map, Function.comp_apply] at ih ⊢
      rw [show
          millerSparseTermEval n p q v
              { t with nCoeffs := spCoeffNeg t.nCoeffs } =
            -millerSparseTermEval n p q v t by
          simp [millerSparseTermEval]]
      rw [ih]
      ring

@[simp] private theorem millerSparseEval_sub
    (P Q : SparsePoly) (n p q v : ℝ) :
    millerSparseEval (P - Q) n p q v =
      millerSparseEval P n p q v -
        millerSparseEval Q n p q v := by
  rw [sparse_sub_eq, millerSparseEval_add,
    millerSparseEval_neg]
  ring

private theorem millerSparseTermEval_mul
    (n p q v : ℝ) (a b : SparseTerm) :
    millerSparseTermEval n p q v (sparseMulTerm a b) =
      millerSparseTermEval n p q v a *
        millerSparseTermEval n p q v b := by
  rcases a with ⟨ac, ap, aq, av⟩
  rcases b with ⟨bc, bp, bq, bv⟩
  simp [sparseMulTerm, millerSparseTermEval, pow_add]
  ring

private theorem millerSparseEval_map_mul_left
    (a : SparseTerm) (qs : List SparseTerm) (n p q v : ℝ) :
    ((qs.map (sparseMulTerm a)).map
        (millerSparseTermEval n p q v)).sum =
      millerSparseTermEval n p q v a *
        (qs.map (millerSparseTermEval n p q v)).sum := by
  induction qs with
  | nil => simp
  | cons b qs ih =>
      simp only [List.map_cons, List.sum_cons, ih,
        millerSparseTermEval_mul]
      ring

private theorem millerSparseEval_mul_raw
    (P Q : SparsePoly) (n p q v : ℝ) :
    ((P.terms.flatMap fun a =>
        Q.terms.map (sparseMulTerm a)).map
          (millerSparseTermEval n p q v)).sum =
      millerSparseEval P n p q v *
        millerSparseEval Q n p q v := by
  rcases P with ⟨ps⟩
  rcases Q with ⟨qs⟩
  simp only [millerSparseEval]
  induction ps with
  | nil => simp
  | cons a ps ih =>
      simp only [List.flatMap_cons, List.map_append, List.sum_append,
        List.map_cons, List.sum_cons, ih]
      rw [millerSparseEval_map_mul_left]
      ring

@[simp] private theorem millerSparseEval_mul
    (P Q : SparsePoly) (n p q v : ℝ) :
    millerSparseEval (P * Q) n p q v =
      millerSparseEval P n p q v *
        millerSparseEval Q n p q v := by
  rw [sparse_mul_eq, millerSparseEval_normalize]
  exact millerSparseEval_mul_raw P Q n p q v

@[simp] private theorem millerSparseEval_pow
    (P : SparsePoly) (k : ℕ) (n p q v : ℝ) :
    millerSparseEval (P ^ k) n p q v =
      millerSparseEval P n p q v ^ k := by
  induction k with
  | zero => simp
  | succ k ih =>
      change millerSparseEval ((P ^ k) * P) n p q v =
        millerSparseEval P n p q v ^ k *
          millerSparseEval P n p q v
      rw [millerSparseEval_mul, ih]

private theorem millerDualEval_ofNCoeffs
    (n p q v : ℝ) (pExp qExp vExp k : ℕ) (cs : List ℤ) :
    ((dualCertOfNCoeffs pExp qExp vExp k cs).map
        (dualCertMonomialEval n p q v)).sum =
      n ^ k * millerSpCoeffEval n cs *
        p ^ pExp * q ^ qExp * v ^ vExp := by
  induction cs generalizing k with
  | nil =>
      simp [dualCertOfNCoeffs, millerSpCoeffEval]
  | cons c cs ih =>
      simp only [dualCertOfNCoeffs, List.map_cons, List.sum_cons,
        millerSpCoeffEval]
      rw [ih (k + 1)]
      simp [dualCertMonomialEval, pow_succ]
      push_cast
      ring

private theorem millerSparseTerm_eval_ofDual
    (n p q v : ℝ) (t : DualCertTerm) :
    millerSparseTermEval n p q v (sparseTermOfDual t) =
      ((dualCertOfTerm t).map
        (dualCertMonomialEval n p q v)).sum := by
  rcases t with ⟨cs, pExp, qExp, vExp⟩
  unfold dualCertOfTerm
  rw [millerDualEval_ofNCoeffs]
  simp [millerSparseTermEval, sparseTermOfDual, dualCertOfTerm,
    millerSpCoeffEval_normalize]

private theorem millerSparseEval_ofDualTerms_raw
    (ts : List DualCertTerm) (n p q v : ℝ) :
    millerSparseEval ⟨ts.map sparseTermOfDual⟩ n p q v =
      dualCertEval (dualCertOfTerms ts) n p q v := by
  induction ts with
  | nil =>
      simp [millerSparseEval, dualCertEval, dualCertOfTerms]
  | cons t ts ih =>
      simp only [List.map_cons]
      unfold millerSparseEval dualCertEval dualCertOfTerms
      unfold millerSparseEval dualCertEval dualCertOfTerms at ih
      simp only [List.flatMap_cons, List.map_append, List.sum_append,
        List.map_cons, List.sum_cons]
      rw [millerSparseTerm_eval_ofDual, ih]

private theorem millerSparseEval_ofDualTerms
    (ts : List DualCertTerm) (n p q v : ℝ) :
    millerSparseEval (SparsePoly.ofDualTerms ts) n p q v =
      dualCertEval (dualCertOfTerms ts) n p q v := by
  rw [SparsePoly.ofDualTerms, millerSparseEval_normalize]
  exact millerSparseEval_ofDualTerms_raw ts n p q v

private theorem millerSparseTerm_hasDerivAt_p
    (t : SparseTerm) (n p q v : ℝ) :
    HasDerivAt (fun x => millerSparseTermEval n x q v t)
      (millerSparseTermEval n p q v
        (if t.pExp = 0 then { t with nCoeffs := [] }
        else { t with
          nCoeffs := spCoeffScale t.pExp t.nCoeffs
          pExp := t.pExp - 1 })) p := by
  rcases t with ⟨cs, pExp, qExp, vExp⟩
  cases pExp with
  | zero =>
      simpa [millerSparseTermEval, millerSpCoeffEval] using
        (hasDerivAt_const p
          (millerSpCoeffEval n cs * q ^ qExp * v ^ vExp))
  | succ pExp =>
      have hd := (hasDerivAt_pow (pExp + 1) p).const_mul
        (millerSpCoeffEval n cs * q ^ qExp * v ^ vExp)
      convert hd using 1
      · funext x
        simp [millerSparseTermEval]
        ring
      · simp [millerSparseTermEval, millerSpCoeffEval_scale,
          Nat.cast_add, Nat.cast_one, Nat.succ_eq_add_one]
        ring

private theorem millerSparseTerm_hasDerivAt_q
    (t : SparseTerm) (n p q v : ℝ) :
    HasDerivAt (fun x => millerSparseTermEval n p x v t)
      (millerSparseTermEval n p q v
        (if t.qExp = 0 then { t with nCoeffs := [] }
        else { t with
          nCoeffs := spCoeffScale t.qExp t.nCoeffs
          qExp := t.qExp - 1 })) q := by
  rcases t with ⟨cs, pExp, qExp, vExp⟩
  cases qExp with
  | zero =>
      simpa [millerSparseTermEval, millerSpCoeffEval] using
        (hasDerivAt_const q
          (millerSpCoeffEval n cs * p ^ pExp * v ^ vExp))
  | succ qExp =>
      have hd := (hasDerivAt_pow (qExp + 1) q).const_mul
        (millerSpCoeffEval n cs * p ^ pExp * v ^ vExp)
      convert hd using 1
      · funext x
        simp [millerSparseTermEval]
        ring
      · simp [millerSparseTermEval, millerSpCoeffEval_scale,
          Nat.cast_add, Nat.cast_one, Nat.succ_eq_add_one]
        ring

private theorem millerSparseTerm_hasDerivAt_v
    (t : SparseTerm) (n p q v : ℝ) :
    HasDerivAt (fun x => millerSparseTermEval n p q x t)
      (millerSparseTermEval n p q v
        (if t.vExp = 0 then { t with nCoeffs := [] }
        else { t with
          nCoeffs := spCoeffScale t.vExp t.nCoeffs
          vExp := t.vExp - 1 })) v := by
  rcases t with ⟨cs, pExp, qExp, vExp⟩
  cases vExp with
  | zero =>
      simpa [millerSparseTermEval, millerSpCoeffEval] using
        (hasDerivAt_const v
          (millerSpCoeffEval n cs * p ^ pExp * q ^ qExp))
  | succ vExp =>
      have hd := (hasDerivAt_pow (vExp + 1) v).const_mul
        (millerSpCoeffEval n cs * p ^ pExp * q ^ qExp)
      convert hd using 1 <;>
        simp [millerSparseTermEval, millerSpCoeffEval_scale,
          Nat.cast_add, Nat.cast_one, Nat.succ_eq_add_one] <;>
        ring

private theorem millerSparseEval_hasDerivAt_p
    (P : SparsePoly) (n p q v : ℝ) :
    HasDerivAt (fun x => millerSparseEval P n x q v)
      (millerSparseEval (spDerivP P) n p q v) p := by
  rcases P with ⟨ts⟩
  unfold spDerivP
  rw [millerSparseEval_normalize]
  unfold millerSparseEval
  induction ts with
  | nil => simpa using hasDerivAt_const p (0 : ℝ)
  | cons t ts ih =>
      simp only [List.map_cons, List.sum_cons]
      exact (millerSparseTerm_hasDerivAt_p t n p q v).add ih

private theorem millerSparseEval_hasDerivAt_q
    (P : SparsePoly) (n p q v : ℝ) :
    HasDerivAt (fun x => millerSparseEval P n p x v)
      (millerSparseEval (spDerivQ P) n p q v) q := by
  rcases P with ⟨ts⟩
  unfold spDerivQ
  rw [millerSparseEval_normalize]
  unfold millerSparseEval
  induction ts with
  | nil => simpa using hasDerivAt_const q (0 : ℝ)
  | cons t ts ih =>
      simp only [List.map_cons, List.sum_cons]
      exact (millerSparseTerm_hasDerivAt_q t n p q v).add ih

private theorem millerSparseEval_hasDerivAt_v
    (P : SparsePoly) (n p q v : ℝ) :
    HasDerivAt (fun x => millerSparseEval P n p q x)
      (millerSparseEval (spDerivV P) n p q v) v := by
  rcases P with ⟨ts⟩
  unfold spDerivV
  rw [millerSparseEval_normalize]
  unfold millerSparseEval
  induction ts with
  | nil => simpa using hasDerivAt_const v (0 : ℝ)
  | cons t ts ih =>
      simp only [List.map_cons, List.sum_cons]
      exact (millerSparseTerm_hasDerivAt_v t n p q v).add ih

private theorem millerSparseEval_pDeriv_ofDualTerms
    (ts : List DualCertTerm) (n p q v : ℝ) :
    millerSparseEval
        (spDerivP (SparsePoly.ofDualTerms ts)) n p q v =
      dualCertEval
        (dualCertPDeriv (dualCertOfTerms ts)) n p q v := by
  have hs := millerSparseEval_hasDerivAt_p
    (SparsePoly.ofDualTerms ts) n p q v
  have hs' :
      HasDerivAt
        (fun x => dualCertEval (dualCertOfTerms ts) n x q v)
        (millerSparseEval
          (spDerivP (SparsePoly.ofDualTerms ts)) n p q v) p := by
    simpa only [millerSparseEval_ofDualTerms] using hs
  exact hs'.unique
    (dualCertEval_hasDerivAt_p (dualCertOfTerms ts) n p q v)

private theorem millerSparseEval_qDeriv_ofDualTerms
    (ts : List DualCertTerm) (n p q v : ℝ) :
    millerSparseEval
        (spDerivQ (SparsePoly.ofDualTerms ts)) n p q v =
      dualCertEval
        (dualCertQDeriv (dualCertOfTerms ts)) n p q v := by
  have hs := millerSparseEval_hasDerivAt_q
    (SparsePoly.ofDualTerms ts) n p q v
  have hs' :
      HasDerivAt
        (fun x => dualCertEval (dualCertOfTerms ts) n p x v)
        (millerSparseEval
          (spDerivQ (SparsePoly.ofDualTerms ts)) n p q v) q := by
    simpa only [millerSparseEval_ofDualTerms] using hs
  exact hs'.unique
    (dualCertEval_hasDerivAt_q (dualCertOfTerms ts) n p q v)

private theorem millerSparseEval_vDeriv_ofDualTerms
    (ts : List DualCertTerm) (n p q v : ℝ) :
    millerSparseEval
        (spDerivV (SparsePoly.ofDualTerms ts)) n p q v =
      dualCertEval
        (dualCertVDeriv (dualCertOfTerms ts)) n p q v := by
  have hs := millerSparseEval_hasDerivAt_v
    (SparsePoly.ofDualTerms ts) n p q v
  have hs' :
      HasDerivAt
        (fun x => dualCertEval (dualCertOfTerms ts) n p q x)
        (millerSparseEval
          (spDerivV (SparsePoly.ofDualTerms ts)) n p q v) v := by
    simpa only [millerSparseEval_ofDualTerms] using hs
  exact hs'.unique
    (dualCertEval_hasDerivAt_v (dualCertOfTerms ts) n p q v)

private def millerDualTermEval
    (t : DualCertTerm) (n p q v : ℝ) : ℝ :=
  millerSpCoeffEval n t.nCoeffs *
    p ^ t.pExp * q ^ t.qExp * v ^ t.vExp

private def millerDualTermPDerivEval
    (t : DualCertTerm) (n p q v : ℝ) : ℝ :=
  (t.pExp : ℝ) * millerSpCoeffEval n t.nCoeffs *
    p ^ (t.pExp - 1) * q ^ t.qExp * v ^ t.vExp

private def millerDualTermQDerivEval
    (t : DualCertTerm) (n p q v : ℝ) : ℝ :=
  (t.qExp : ℝ) * millerSpCoeffEval n t.nCoeffs *
    p ^ t.pExp * q ^ (t.qExp - 1) * v ^ t.vExp

private def millerDualTermVDerivEval
    (t : DualCertTerm) (n p q v : ℝ) : ℝ :=
  (t.vExp : ℝ) * millerSpCoeffEval n t.nCoeffs *
    p ^ t.pExp * q ^ t.qExp * v ^ (t.vExp - 1)

private def millerDualTermsEval
    (ts : List DualCertTerm) (n p q v : ℝ) : ℝ :=
  (ts.map fun t => millerDualTermEval t n p q v).sum

private def millerDualTermsPDerivEval
    (ts : List DualCertTerm) (n p q v : ℝ) : ℝ :=
  (ts.map fun t => millerDualTermPDerivEval t n p q v).sum

private def millerDualTermsQDerivEval
    (ts : List DualCertTerm) (n p q v : ℝ) : ℝ :=
  (ts.map fun t => millerDualTermQDerivEval t n p q v).sum

private def millerDualTermsVDerivEval
    (ts : List DualCertTerm) (n p q v : ℝ) : ℝ :=
  (ts.map fun t => millerDualTermVDerivEval t n p q v).sum

private theorem millerSparseTermEval_ofDual_eq
    (t : DualCertTerm) (n p q v : ℝ) :
    millerSparseTermEval n p q v (sparseTermOfDual t) =
      millerDualTermEval t n p q v := by
  simp [millerSparseTermEval, millerDualTermEval,
    sparseTermOfDual, millerSpCoeffEval_normalize]

private theorem millerDualCertEval_ofTerms_eq
    (ts : List DualCertTerm) (n p q v : ℝ) :
    dualCertEval (dualCertOfTerms ts) n p q v =
      millerDualTermsEval ts n p q v := by
  rw [← millerSparseEval_ofDualTerms]
  rw [SparsePoly.ofDualTerms, millerSparseEval_normalize]
  unfold millerSparseEval millerDualTermsEval
  induction ts with
  | nil => rfl
  | cons t ts ih =>
      simp only [List.map_cons, List.sum_cons]
      rw [millerSparseTermEval_ofDual_eq, ih]

private theorem millerDualTerm_hasDerivAt_p
    (t : DualCertTerm) (n p q v : ℝ) :
    HasDerivAt (fun x => millerDualTermEval t n x q v)
      (millerDualTermPDerivEval t n p q v) p := by
  rcases t with ⟨cs, pe, qe, ve⟩
  cases pe with
  | zero =>
      simpa [millerDualTermEval, millerDualTermPDerivEval] using
        hasDerivAt_const p
          (millerSpCoeffEval n cs * q ^ qe * v ^ ve)
  | succ pe =>
      have h := (hasDerivAt_pow (pe + 1) p).const_mul
        (millerSpCoeffEval n cs * q ^ qe * v ^ ve)
      convert h using 1
      all_goals
        first
        | (funext x
           simp [millerDualTermEval]
           ring)
        | (simp [millerDualTermPDerivEval,
             Nat.cast_add, Nat.cast_one]
           ring)

private theorem millerDualTerm_hasDerivAt_q
    (t : DualCertTerm) (n p q v : ℝ) :
    HasDerivAt (fun x => millerDualTermEval t n p x v)
      (millerDualTermQDerivEval t n p q v) q := by
  rcases t with ⟨cs, pe, qe, ve⟩
  cases qe with
  | zero =>
      simpa [millerDualTermEval, millerDualTermQDerivEval] using
        hasDerivAt_const q
          (millerSpCoeffEval n cs * p ^ pe * v ^ ve)
  | succ qe =>
      have h := (hasDerivAt_pow (qe + 1) q).const_mul
        (millerSpCoeffEval n cs * p ^ pe * v ^ ve)
      convert h using 1
      all_goals
        first
        | (funext x
           simp [millerDualTermEval]
           ring)
        | (simp [millerDualTermQDerivEval,
             Nat.cast_add, Nat.cast_one]
           ring)

private theorem millerDualTerm_hasDerivAt_v
    (t : DualCertTerm) (n p q v : ℝ) :
    HasDerivAt (fun x => millerDualTermEval t n p q x)
      (millerDualTermVDerivEval t n p q v) v := by
  rcases t with ⟨cs, pe, qe, ve⟩
  cases ve with
  | zero =>
      simpa [millerDualTermEval, millerDualTermVDerivEval] using
        hasDerivAt_const v
          (millerSpCoeffEval n cs * p ^ pe * q ^ qe)
  | succ ve =>
      have h := (hasDerivAt_pow (ve + 1) v).const_mul
        (millerSpCoeffEval n cs * p ^ pe * q ^ qe)
      convert h using 1
      all_goals
        first
        | (funext x
           simp [millerDualTermEval]
           ring)
        | (simp [millerDualTermVDerivEval,
             Nat.cast_add, Nat.cast_one]
           ring)

private theorem millerDualTerms_hasDerivAt_p
    (ts : List DualCertTerm) (n p q v : ℝ) :
    HasDerivAt (fun x => millerDualTermsEval ts n x q v)
      (millerDualTermsPDerivEval ts n p q v) p := by
  induction ts with
  | nil =>
      simpa [millerDualTermsEval, millerDualTermsPDerivEval] using
        hasDerivAt_const p (0 : ℝ)
  | cons t ts ih =>
      simp only [millerDualTermsEval, millerDualTermsPDerivEval,
        List.map_cons, List.sum_cons]
      exact (millerDualTerm_hasDerivAt_p t n p q v).add ih

private theorem millerDualTerms_hasDerivAt_q
    (ts : List DualCertTerm) (n p q v : ℝ) :
    HasDerivAt (fun x => millerDualTermsEval ts n p x v)
      (millerDualTermsQDerivEval ts n p q v) q := by
  induction ts with
  | nil =>
      simpa [millerDualTermsEval, millerDualTermsQDerivEval] using
        hasDerivAt_const q (0 : ℝ)
  | cons t ts ih =>
      simp only [millerDualTermsEval, millerDualTermsQDerivEval,
        List.map_cons, List.sum_cons]
      exact (millerDualTerm_hasDerivAt_q t n p q v).add ih

private theorem millerDualTerms_hasDerivAt_v
    (ts : List DualCertTerm) (n p q v : ℝ) :
    HasDerivAt (fun x => millerDualTermsEval ts n p q x)
      (millerDualTermsVDerivEval ts n p q v) v := by
  induction ts with
  | nil =>
      simpa [millerDualTermsEval, millerDualTermsVDerivEval] using
        hasDerivAt_const v (0 : ℝ)
  | cons t ts ih =>
      simp only [millerDualTermsEval, millerDualTermsVDerivEval,
        List.map_cons, List.sum_cons]
      exact (millerDualTerm_hasDerivAt_v t n p q v).add ih

private theorem millerDualCertPDerivEval_ofTerms_eq
    (ts : List DualCertTerm) (n p q v : ℝ) :
    dualCertEval (dualCertPDeriv (dualCertOfTerms ts)) n p q v =
      millerDualTermsPDerivEval ts n p q v := by
  have h :
      HasDerivAt
        (fun x => dualCertEval (dualCertOfTerms ts) n x q v)
        (millerDualTermsPDerivEval ts n p q v) p := by
    simpa only [millerDualCertEval_ofTerms_eq] using
      millerDualTerms_hasDerivAt_p ts n p q v
  exact (dualCertEval_hasDerivAt_p
    (dualCertOfTerms ts) n p q v).unique h

private theorem millerDualCertQDerivEval_ofTerms_eq
    (ts : List DualCertTerm) (n p q v : ℝ) :
    dualCertEval (dualCertQDeriv (dualCertOfTerms ts)) n p q v =
      millerDualTermsQDerivEval ts n p q v := by
  have h :
      HasDerivAt
        (fun x => dualCertEval (dualCertOfTerms ts) n p x v)
        (millerDualTermsQDerivEval ts n p q v) q := by
    simpa only [millerDualCertEval_ofTerms_eq] using
      millerDualTerms_hasDerivAt_q ts n p q v
  exact (dualCertEval_hasDerivAt_q
    (dualCertOfTerms ts) n p q v).unique h

private theorem millerDualCertVDerivEval_ofTerms_eq
    (ts : List DualCertTerm) (n p q v : ℝ) :
    dualCertEval (dualCertVDeriv (dualCertOfTerms ts)) n p q v =
      millerDualTermsVDerivEval ts n p q v := by
  have h :
      HasDerivAt
        (fun x => dualCertEval (dualCertOfTerms ts) n p q x)
        (millerDualTermsVDerivEval ts n p q v) v := by
    simpa only [millerDualCertEval_ofTerms_eq] using
      millerDualTerms_hasDerivAt_v ts n p q v
  exact (dualCertEval_hasDerivAt_v
    (dualCertOfTerms ts) n p q v).unique h

@[simp] private theorem millerSparseEval_spD
    (n p q v : ℝ) :
    millerSparseEval spD n p q v = dualCertD p q v := by
  simp only [spD, millerSparseEval_add, millerSparseEval_mul,
    millerSparseEval_pow, millerSparseEval_ofNat,
    millerSparseEval_p, millerSparseEval_q, millerSparseEval_v]
  simp [dualCertD]

@[simp] private theorem millerSparseEval_spSnum
    (n p q v : ℝ) :
    millerSparseEval spSnum n p q v = dualCertSnum p q v := by
  simp only [spSnum, millerSparseEval_sub, millerSparseEval_mul,
    millerSparseEval_pow, millerSparseEval_ofNat,
    millerSparseEval_p, millerSparseEval_q, millerSparseEval_v]
  simp [dualCertSnum]

@[simp] private theorem millerSparseEval_spNext
    (j : Fin 3) (n p q v : ℝ) :
    millerSparseEval (spNext j) n p q v =
      dualCertNextNum n j p q v := by
  fin_cases j
  · change millerSparseEval (spSnum * spD ^ 2) n p q v =
      dualCertNextNum n 0 p q v
    rw [millerSparseEval_mul, millerSparseEval_spSnum,
      millerSparseEval_pow, millerSparseEval_spD]
    simp [dualCertNextNum]
  · change millerSparseEval
      (spSnum * 2 * (SparsePoly.n + 3) *
        SparsePoly.v * spD) n p q v =
      dualCertNextNum n 1 p q v
    simp only [millerSparseEval_mul, millerSparseEval_add,
      millerSparseEval_ofNat, millerSparseEval_n,
      millerSparseEval_v, millerSparseEval_spD,
      millerSparseEval_spSnum]
    simp [dualCertNextNum]
  · change millerSparseEval
      (spSnum *
        (-(SparsePoly.n + 3) * SparsePoly.v * spD +
          2 * (SparsePoly.n + 3) *
            (2 * SparsePoly.n + 7) * SparsePoly.v ^ 2))
        n p q v =
      dualCertNextNum n 2 p q v
    simp only [millerSparseEval_add, millerSparseEval_neg,
      millerSparseEval_mul, millerSparseEval_pow,
      millerSparseEval_ofNat, millerSparseEval_n,
      millerSparseEval_v, millerSparseEval_spD,
      millerSparseEval_spSnum]
    simp [dualCertNextNum]

@[simp] private theorem millerSparseEval_spCur
    (i : Fin 3) (n p q v : ℝ) :
    millerSparseEval (spCur i) n p q v =
      dualCertCurNum n i p q v := by
  fin_cases i
  · change millerSparseEval (spD ^ 4) n p q v =
      dualCertCurNum n 0 p q v
    rw [millerSparseEval_pow, millerSparseEval_spD]
    simp [dualCertCurNum]
  · change millerSparseEval
      (2 * (SparsePoly.n + 2) * SparsePoly.v * spD ^ 3)
        n p q v =
      dualCertCurNum n 1 p q v
    simp only [millerSparseEval_add, millerSparseEval_mul,
      millerSparseEval_pow, millerSparseEval_ofNat,
      millerSparseEval_n, millerSparseEval_v,
      millerSparseEval_spD]
    simp [dualCertCurNum]
  · change millerSparseEval
      ((-(SparsePoly.n + 2) * SparsePoly.v * spD +
          2 * (SparsePoly.n + 2) *
            (2 * SparsePoly.n + 5) * SparsePoly.v ^ 2) *
        spD ^ 2) n p q v =
      dualCertCurNum n 2 p q v
    simp only [millerSparseEval_add, millerSparseEval_neg,
      millerSparseEval_mul, millerSparseEval_pow,
      millerSparseEval_ofNat, millerSparseEval_n,
      millerSparseEval_v, millerSparseEval_spD]
    simp [dualCertCurNum]

@[simp] private theorem millerSparseEval_spLambda
    (n p q v : ℝ) :
    millerSparseEval spLambda n p q v =
      dualCertLambda n := by
  simp only [spLambda, millerSparseEval_add, millerSparseEval_mul,
    millerSparseEval_pow, millerSparseEval_ofNat,
    millerSparseEval_n]
  simp [dualCertLambda]

@[simp] private theorem millerSparseEval_spDelta
    (n p q v : ℝ) :
    millerSparseEval spDelta n p q v =
      4 * (2 * n + 3) * (n + 2) := by
  change millerSparseEval
    (SparsePoly.const 4 *
      (SparsePoly.const 2 * SparsePoly.n + SparsePoly.const 3) *
      (SparsePoly.n + SparsePoly.const 2)) n p q v = _
  simp only [millerSparseEval_mul, millerSparseEval_add,
    millerSparseEval_const, millerSparseEval_n]
  ring

private theorem millerSparseEval_spOpP_ofDualTerms
    (ts : List DualCertTerm) (n p q v : ℝ) :
    millerSparseEval
        (spOpP (SparsePoly.ofDualTerms ts)) n p q v =
      dualCertOpP n
        (dualCertEval (dualCertOfTerms ts) n p q v)
        (dualCertEval
          (dualCertPDeriv (dualCertOfTerms ts)) n p q v)
        p q v := by
  change millerSparseEval
      (spD *
          (SparsePoly.p * (1 - SparsePoly.p ^ 2) *
              spDerivP (SparsePoly.ofDualTerms ts) +
            ((2 * SparsePoly.n + 7) -
                (4 * SparsePoly.n + 9) * SparsePoly.p ^ 2) *
              SparsePoly.ofDualTerms ts) -
        (2 * SparsePoly.n + 7) * SparsePoly.p *
          (1 - SparsePoly.p ^ 2) *
          (SparsePoly.q * (1 + SparsePoly.v ^ 2)) *
          SparsePoly.ofDualTerms ts) n p q v = _
  simp only [millerSparseEval_sub, millerSparseEval_add,
    millerSparseEval_mul, millerSparseEval_pow,
    millerSparseEval_ofNat, millerSparseEval_n,
    millerSparseEval_p, millerSparseEval_q,
    millerSparseEval_v, millerSparseEval_spD,
    millerSparseEval_ofDualTerms,
    millerSparseEval_pDeriv_ofDualTerms]
  simp [dualCertOpP, dualCertD]

private theorem millerSparseEval_spOpQ_ofDualTerms
    (ts : List DualCertTerm) (n p q v : ℝ) :
    millerSparseEval
        (spOpQ (SparsePoly.ofDualTerms ts)) n p q v =
      dualCertOpQ n
        (dualCertEval (dualCertOfTerms ts) n p q v)
        (dualCertEval
          (dualCertQDeriv (dualCertOfTerms ts)) n p q v)
        p q v := by
  change millerSparseEval
      (spD *
          (SparsePoly.q * (1 - SparsePoly.q ^ 2) *
              spDerivQ (SparsePoly.ofDualTerms ts) +
            ((2 * SparsePoly.n + 6) -
                (4 * SparsePoly.n + 10) * SparsePoly.q ^ 2) *
              SparsePoly.ofDualTerms ts) -
        (2 * SparsePoly.n + 7) * SparsePoly.q *
          (1 - SparsePoly.q ^ 2) *
          (SparsePoly.p * (1 + SparsePoly.v ^ 2)) *
          SparsePoly.ofDualTerms ts) n p q v = _
  simp only [millerSparseEval_sub, millerSparseEval_add,
    millerSparseEval_mul, millerSparseEval_pow,
    millerSparseEval_ofNat, millerSparseEval_n,
    millerSparseEval_p, millerSparseEval_q,
    millerSparseEval_v, millerSparseEval_spD,
    millerSparseEval_ofDualTerms,
    millerSparseEval_qDeriv_ofDualTerms]
  simp [dualCertOpQ, dualCertD]

private theorem millerSparseEval_spOpV_ofDualTerms
    (ts : List DualCertTerm) (n p q v : ℝ) :
    millerSparseEval
        (spOpV (SparsePoly.ofDualTerms ts)) n p q v =
      dualCertOpV n
        (dualCertEval (dualCertOfTerms ts) n p q v)
        (dualCertEval
          (dualCertVDeriv (dualCertOfTerms ts)) n p q v)
        p q v := by
  change millerSparseEval
      (spD *
          (SparsePoly.v * (1 - SparsePoly.v ^ 2) *
              spDerivV (SparsePoly.ofDualTerms ts) +
            ((2 * SparsePoly.n + 4) -
                (2 * SparsePoly.n + 6) * SparsePoly.v ^ 2) *
              SparsePoly.ofDualTerms ts) -
        (2 * SparsePoly.n + 7) * SparsePoly.v *
          (1 - SparsePoly.v ^ 2) *
          (2 * SparsePoly.p * SparsePoly.q *
              SparsePoly.v + 2) *
          SparsePoly.ofDualTerms ts) n p q v = _
  simp only [millerSparseEval_sub, millerSparseEval_add,
    millerSparseEval_mul, millerSparseEval_pow,
    millerSparseEval_ofNat, millerSparseEval_n,
    millerSparseEval_p, millerSparseEval_q,
    millerSparseEval_v, millerSparseEval_spD,
    millerSparseEval_ofDualTerms,
    millerSparseEval_vDeriv_ofDualTerms]
  simp [dualCertOpV, dualCertD]

/-! Data-only copies of the two remaining certificate rows.  Keeping these
lists here avoids importing the legacy dense normalization theorems, whose
elaboration is prohibitively expensive on a memory-constrained machine. -/

private def dualCertPp1Terms : List DualCertTerm := [
    { nCoeffs := [-1518912, -6060960, -10667088, -10868952, -7066344, -3040072, -865512, -157248, -16544, -768], pExp := 3, qExp := 5, vExp := 4 },
    { nCoeffs := [-1518912, -6060960, -10667088, -10868952, -7066344, -3040072, -865512, -157248, -16544, -768], pExp := 3, qExp := 5, vExp := 2 },
    { nCoeffs := [1518912, 6060960, 10667088, 10868952, 7066344, 3040072, 865512, 157248, 16544, 768], pExp := 3, qExp := 3, vExp := 4 },
    { nCoeffs := [1518912, 6060960, 10667088, 10868952, 7066344, 3040072, 865512, 157248, 16544, 768], pExp := 3, qExp := 3, vExp := 2 },
    { nCoeffs := [-4959360, -19528704, -33896064, -34040672, -21799576, -9232360, -2585856, -461888, -47744, -2176], pExp := 2, qExp := 4, vExp := 3 },
    { nCoeffs := [-1486080, -7805952, -18162048, -24674816, -21695376, -12908416, -5266352, -1455424, -260864, -27392, -1280], pExp := 2, qExp := 4, vExp := 1 },
    { nCoeffs := [4959360, 19528704, 33896064, 34040672, 21799576, 9232360, 2585856, 461888, 47744, 2176], pExp := 2, qExp := 2, vExp := 3 },
    { nCoeffs := [297216, 1442304, 3055488, 3712768, 2853968, 1440096, 477232, 100192, 12096, 640], pExp := 2, qExp := 2, vExp := 1 },
    { nCoeffs := [1012608, 4378176, 8233248, 8868528, 6044512, 2707936, 798432, 149552, 16160, 768], pExp := 1, qExp := 5, vExp := 4 },
    { nCoeffs := [1012608, 4378176, 8233248, 8868528, 6044512, 2707936, 798432, 149552, 16160, 768], pExp := 1, qExp := 5, vExp := 2 },
    { nCoeffs := [-1012608, -4378176, -8233248, -8868528, -6044512, -2707936, -798432, -149552, -16160, -768], pExp := 1, qExp := 3, vExp := 4 },
    { nCoeffs := [374400, 1824192, 3813856, 4500304, 3311808, 1580256, 489984, 95376, 10592, 512], pExp := 1, qExp := 3, vExp := 2 },
    { nCoeffs := [7925760, 40310784, 90806272, 119493632, 101845760, 58802048, 23308288, 6267008, 1094400, 112128, 5120], pExp := 1, qExp := 3, vExp := 0 },
    { nCoeffs := [198144, 1622016, 4999936, 8251392, 8339424, 5486208, 2403584, 697088, 128672, 13696, 640], pExp := 1, qExp := 1, vExp := 2 },
    { nCoeffs := [792576, 4110336, 9452032, 12708864, 11075712, 6542208, 2654048, 730496, 130592, 13696, 640], pExp := 1, qExp := 1, vExp := 0 },
    { nCoeffs := [-211680, -805392, -1348008, -1302908, -801568, -325564, -87312, -14912, -1472, -64], pExp := 0, qExp := 4, vExp := 5 },
    { nCoeffs := [4558464, 19411776, 36052896, 38431536, 25961168, 11541104, 3379984, 629344, 67648, 3200], pExp := 0, qExp := 4, vExp := 3 },
    { nCoeffs := [-211680, -805392, -1348008, -1302908, -801568, -325564, -87312, -14912, -1472, -64], pExp := 0, qExp := 4, vExp := 1 },
    { nCoeffs := [127008, 474768, 779976, 739356, 445780, 177340, 46564, 7784, 752, 32], pExp := 0, qExp := 2, vExp := 5 },
    { nCoeffs := [-15088896, -71267904, -149119968, -182516592, -145001592, -78248592, -29081272, -7356688, -1213152, -117824, -5120], pExp := 0, qExp := 2, vExp := 3 },
    { nCoeffs := [39755808, 198065808, 436637384, 561891356, 468085844, 264059068, 102251236, 26857064, 4582128, 458784, 20480], pExp := 0, qExp := 2, vExp := 1 },
    { nCoeffs := [7133184, 34219008, 72597504, 89921024, 72139136, 39219008, 14647616, 3713856, 612160, 59264, 2560], pExp := 0, qExp := 0, vExp := 3 },
    { nCoeffs := [-24569856, -120683520, -261839872, -331076096, -270568576, -149509824, -56626048, -14526784, -2417408, -235776, -10240], pExp := 0, qExp := 0, vExp := 1 }
  ]

private def dualCertPq1Terms : List DualCertTerm := [
    { nCoeffs := [1188864, 6363648, 15106560, 20962048, 18841408, 11468320, 4789120, 1355232, 248768, 26752, 1280], pExp := 4, qExp := 2, vExp := 1 },
    { nCoeffs := [-5944320, -31223808, -72648192, -98699264, -86781504, -51633664, -21065408, -5821696, -1043456, -109568, -5120], pExp := 3, qExp := 1, vExp := 0 },
    { nCoeffs := [84672, 372960, 712176, 775496, 531592, 238216, 69864, 12944, 1376, 64], pExp := 2, qExp := 2, vExp := 5 },
    { nCoeffs := [-2674944, -11791872, -22654848, -24975616, -17450192, -8025888, -2432752, -469024, -52224, -2560], pExp := 2, qExp := 2, vExp := 3 },
    { nCoeffs := [-509760, -2511648, -5398800, -6650040, -5176344, -2641976, -884600, -187440, -22816, -1216], pExp := 2, qExp := 2, vExp := 1 },
    { nCoeffs := [10699776, 52517376, 114203136, 145212160, 119752000, 67003936, 25782784, 6741600, 1146944, 114688, 5120], pExp := 2, qExp := 0, vExp := 3 },
    { nCoeffs := [-35665920, -179417088, -399541248, -519563264, -437511552, -249544960, -97718912, -25958656, -4479488, -453632, -20480], pExp := 2, qExp := 0, vExp := 1 },
    { nCoeffs := [168768, 673440, 1147728, 1095512, 642248, 237240, 53992, 6928, 384], pExp := 1, qExp := 3, vExp := 4 },
    { nCoeffs := [168768, 673440, 1147728, 1095512, 642248, 237240, 53992, 6928, 384], pExp := 1, qExp := 3, vExp := 2 },
    { nCoeffs := [275760, 744816, 800988, 411904, 76372, -20856, -13768, -2704, -192], pExp := 1, qExp := 1, vExp := 4 },
    { nCoeffs := [-1292976, -5298432, -10471820, -13031700, -11067028, -6549436, -2688160, -747520, -134016, -13952, -640], pExp := 1, qExp := 1, vExp := 2 },
    { nCoeffs := [1740384, 7211664, 12761896, 12412956, 6978472, 2011340, 18184, -199536, -69696, -10624, -640], pExp := 1, qExp := 1, vExp := 0 },
    { nCoeffs := [-105840, -402696, -674004, -651454, -400784, -162782, -43656, -7456, -736, -32], pExp := 0, qExp := 2, vExp := 5 },
    { nCoeffs := [2659104, 11196912, 20533848, 21582596, 14355072, 6274148, 1803792, 329184, 34624, 1600], pExp := 0, qExp := 2, vExp := 3 },
    { nCoeffs := [-105840, -402696, -674004, -651454, -400784, -162782, -43656, -7456, -736, -32], pExp := 0, qExp := 2, vExp := 1 },
    { nCoeffs := [-8434944, -40041984, -83954304, -102631040, -81151440, -43423424, -15939376, -3966080, -640576, -60672, -2560], pExp := 0, qExp := 0, vExp := 3 },
    { nCoeffs := [27981216, 136898208, 295079816, 369875488, 299101832, 163259056, 60976224, 15399552, 2518272, 240896, 10240], pExp := 0, qExp := 0, vExp := 1 }
  ]

private def dualCertPv1Terms : List DualCertTerm := [
    { nCoeffs := [594432, 3181824, 7553280, 10481024, 9420704, 5734160, 2394560, 677616, 124384, 13376, 640], pExp := 3, qExp := 1, vExp := 0 },
    { nCoeffs := [3566592, 18298368, 41605632, 55291136, 47612864, 27784928, 11135168, 3027744, 534784, 55424, 2560], pExp := 2, qExp := 0, vExp := 1 },
    { nCoeffs := [-241056, -1875312, -5396184, -8381412, -8085784, -5157044, -2224696, -644848, -120640, -13184, -640], pExp := 1, qExp := 3, vExp := 0 },
    { nCoeffs := [-14112, -283248, -877720, -1255332, -1031664, -527380, -171120, -34368, -3904, -192], pExp := 1, qExp := 1, vExp := 2 },
    { nCoeffs := [-198144, -928512, -1898752, -2227840, -1655008, -808048, -259488, -52880, -6208, -320], pExp := 1, qExp := 1, vExp := 0 },
    { nCoeffs := [-105840, -402696, -674004, -651454, -400784, -162782, -43656, -7456, -736, -32], pExp := 0, qExp := 4, vExp := 3 },
    { nCoeffs := [-105840, -402696, -674004, -651454, -400784, -162782, -43656, -7456, -736, -32], pExp := 0, qExp := 4, vExp := 1 },
    { nCoeffs := [21168, 72072, 105972, 87902, 44996, 14558, 2908, 328, 16], pExp := 0, qExp := 2, vExp := 3 },
    { nCoeffs := [-1184112, -9183960, -25997556, -39559762, -37322916, -23258274, -9798244, -2772728, -506352, -54016, -2560], pExp := 0, qExp := 2, vExp := 1 },
    { nCoeffs := [-2050848, -8286336, -14985416, -15899928, -10885128, -4973960, -1513264, -294912, -33344, -1664], pExp := 0, qExp := 0, vExp := 1 }
  ]

private def dualCertPp1Poly : DualCertPoly :=
  dualCertOfTerms dualCertPp1Terms

private def dualCertPq1Poly : DualCertPoly :=
  dualCertOfTerms dualCertPq1Terms

private def dualCertPv1Poly : DualCertPoly :=
  dualCertOfTerms dualCertPv1Terms

private def dualCertPp1 (n p q v : ℝ) : ℝ :=
  dualCertEval dualCertPp1Poly n p q v

private def dualCertPq1 (n p q v : ℝ) : ℝ :=
  dualCertEval dualCertPq1Poly n p q v

private def dualCertPv1 (n p q v : ℝ) : ℝ :=
  dualCertEval dualCertPv1Poly n p q v

private def dualCertDPp1 (n p q v : ℝ) : ℝ :=
  dualCertEval (dualCertPDeriv dualCertPp1Poly) n p q v

private def dualCertDPq1 (n p q v : ℝ) : ℝ :=
  dualCertEval (dualCertQDeriv dualCertPq1Poly) n p q v

private def dualCertDPv1 (n p q v : ℝ) : ℝ :=
  dualCertEval (dualCertVDeriv dualCertPv1Poly) n p q v

private def dualCertPp2Terms : List DualCertTerm := [
    { nCoeffs := [-1404864, -6624288, -13926672, -17193624, -13806192, -7535248, -2831296, -723256, -120224, -11744, -512], pExp := 3, qExp := 5, vExp := 4 },
    { nCoeffs := [-1404864, -6624288, -13926672, -17193624, -13806192, -7535248, -2831296, -723256, -120224, -11744, -512], pExp := 3, qExp := 5, vExp := 2 },
    { nCoeffs := [1404864, 6624288, 13926672, 17193624, 13806192, 7535248, 2831296, 723256, 120224, 11744, 512], pExp := 3, qExp := 3, vExp := 4 },
    { nCoeffs := [1404864, 6624288, 13926672, 17193624, 13806192, 7535248, 2831296, 723256, 120224, 11744, 512], pExp := 3, qExp := 3, vExp := 2 },
    { nCoeffs := [-4838400, -22570560, -46914720, -57227536, -45373112, -24435080, -9052944, -2278624, -372928, -35840, -1536], pExp := 2, qExp := 4, vExp := 3 },
    { nCoeffs := [-1157760, -6963264, -18785184, -30017968, -31580688, -22975168, -11797328, -4276624, -1072832, -177408, -17408, -768], pExp := 2, qExp := 4, vExp := 1 },
    { nCoeffs := [4838400, 22570560, 46914720, 57227536, 45373112, 24435080, 9052944, 2278624, 372928, 35840, 1536], pExp := 2, qExp := 2, vExp := 3 },
    { nCoeffs := [231552, 1300032, 3237024, 4708784, 4432624, 2821984, 1230672, 363056, 69344, 7744, 384], pExp := 2, qExp := 2, vExp := 1 },
    { nCoeffs := [936576, 4728384, 10548384, 13714608, 11523472, 6545312, 2547488, 671376, 114736, 11488, 512], pExp := 1, qExp := 5, vExp := 4 },
    { nCoeffs := [936576, 4728384, 10548384, 13714608, 11523472, 6545312, 2547488, 671376, 114736, 11488, 512], pExp := 1, qExp := 5, vExp := 2 },
    { nCoeffs := [-936576, -4728384, -10548384, -13714608, -11523472, -6545312, -2547488, -671376, -114736, -11488, -512], pExp := 1, qExp := 3, vExp := 4 },
    { nCoeffs := [144000, 926784, 2520992, 3862832, 3722192, 2370336, 1014528, 289104, 52624, 5536, 256], pExp := 1, qExp := 3, vExp := 2 },
    { nCoeffs := [6174720, 36108288, 94684160, 147066880, 150433664, 106460288, 53212672, 18793088, 4597376, 742144, 71168, 3072], pExp := 1, qExp := 3, vExp := 0 },
    { nCoeffs := [154368, 1381248, 4858304, 9399072, 11403808, 9238304, 5141952, 1977184, 516704, 87648, 8704, 384], pExp := 1, qExp := 1, vExp := 2 },
    { nCoeffs := [617472, 3672576, 9804800, 15518976, 16189120, 11692064, 5967456, 2152960, 538208, 88800, 8704, 384], pExp := 1, qExp := 1, vExp := 0 },
    { nCoeffs := [105840, 402696, 674004, 651454, 400784, 162782, 43656, 7456, 736, 32], pExp := 0, qExp := 4, vExp := 5 },
    { nCoeffs := [3492288, 17964000, 40633584, 53379720, 45202328, 25826568, 10096824, 2669936, 457440, 45888, 2048], pExp := 0, qExp := 4, vExp := 3 },
    { nCoeffs := [105840, 402696, 674004, 651454, 400784, 162782, 43656, 7456, 736, 32], pExp := 0, qExp := 4, vExp := 1 },
    { nCoeffs := [-63504, -237384, -389988, -369678, -222890, -88670, -23282, -3892, -376, -16], pExp := 0, qExp := 2, vExp := 5 },
    { nCoeffs := [-11912832, -65558880, -161363664, -234956776, -225248644, -149494248, -70171108, -23318056, -5380368, -821536, -74752, -3072], pExp := 0, qExp := 2, vExp := 3 },
    { nCoeffs := [30810096, 177216696, 456520348, 695877874, 697955414, 483990946, 236920846, 81914828, 19613320, 3098608, 290816, 12288], pExp := 0, qExp := 2, vExp := 1 },
    { nCoeffs := [5557248, 30892032, 76881408, 113208832, 109703872, 73520256, 34795456, 11636800, 2696448, 412480, 37504, 1536], pExp := 0, qExp := 0, vExp := 3 },
    { nCoeffs := [-19141632, -108601344, -275664896, -413571328, -407771328, -277620928, -133255296, -45117248, -10564800, -1630208, -149248, -6144], pExp := 0, qExp := 0, vExp := 1 }
  ]

private def dualCertPq2Terms : List DualCertTerm := [
    { nCoeffs := [926208, 5663232, 15548160, 25309184, 27148064, 20153184, 10566656, 3913568, 1003488, 169664, 17024, 768], pExp := 4, qExp := 2, vExp := 1 },
    { nCoeffs := [-4631040, -27853056, -75140736, -120071872, -126322752, -91900672, -47189312, -17106496, -4291328, -709632, -69632, -3072], pExp := 3, qExp := 1, vExp := 0 },
    { nCoeffs := [-42336, -186480, -356088, -387748, -265796, -119108, -34932, -6472, -688, -32], pExp := 2, qExp := 2, vExp := 5 },
    { nCoeffs := [-2083968, -10774080, -24653472, -32926640, -28459312, -16650208, -6683280, -1818608, -321248, -33280, -1536], pExp := 2, qExp := 2, vExp := 3 },
    { nCoeffs := [-505440, -2786544, -6830136, -9805316, -9131044, -5763076, -2496276, -732584, -139376, -15520, -768], pExp := 2, qExp := 2, vExp := 1 },
    { nCoeffs := [8335872, 47264256, 120162048, 181013504, 179690528, 123519456, 60033536, 20640992, 4922208, 775616, 72704, 3072], pExp := 2, qExp := 0, vExp := 3 },
    { nCoeffs := [-27786240, -160943616, -417823488, -642257536, -649956480, -454960384, -224897408, -78545536, -19001600, -3033600, -287744, -12288], pExp := 2, qExp := 0, vExp := 1 },
    { nCoeffs := [156096, 736032, 1512720, 1781528, 1326736, 648640, 208368, 42440, 4976, 256], pExp := 1, qExp := 3, vExp := 4 },
    { nCoeffs := [156096, 736032, 1512720, 1781528, 1326736, 648640, 208368, 42440, 4976, 256], pExp := 1, qExp := 3, vExp := 2 },
    { nCoeffs := [-378360, -1445160, -2487078, -2535236, -1686046, -756832, -228480, -44552, -5072, -256], pExp := 1, qExp := 1, vExp := 4 },
    { nCoeffs := [97848, -1038000, -5487154, -11445514, -13835758, -10908094, -5875888, -2188240, -555328, -91744, -8896, -384], pExp := 1, qExp := 1, vExp := 2 },
    { nCoeffs := [1016208, 5374008, 12130124, 15297426, 11682552, 5304202, 1121264, -153832, -173040, -49760, -6848, -384], pExp := 1, qExp := 1, vExp := 0 },
    { nCoeffs := [52920, 201348, 337002, 325727, 200392, 81391, 21828, 3728, 368, 16], pExp := 0, qExp := 2, vExp := 5 },
    { nCoeffs := [2037168, 10381992, 23236260, 30165014, 25209640, 14196006, 5462240, 1419520, 238656, 23456, 1024], pExp := 0, qExp := 2, vExp := 3 },
    { nCoeffs := [52920, 201348, 337002, 325727, 200392, 81391, 21828, 3728, 368, 16], pExp := 0, qExp := 2, vExp := 1 },
    { nCoeffs := [-6709248, -36877056, -90650112, -131703168, -125783424, -82983600, -38616192, -12682320, -2882112, -431808, -38400, -1536], pExp := 0, qExp := 0, vExp := 3 },
    { nCoeffs := [22138992, 124365456, 312733420, 464749112, 453638508, 305480832, 144866320, 48395744, 11165120, 1694592, 152320, 6144], pExp := 0, qExp := 0, vExp := 1 }
  ]

private def dualCertPv2Terms : List DualCertTerm := [
    { nCoeffs := [463104, 2831616, 7774080, 12654592, 13574032, 10076592, 5283328, 1956784, 501744, 84832, 8512, 384], pExp := 3, qExp := 1, vExp := 0 },
    { nCoeffs := [2778624, 16372224, 43280640, 67804672, 69986656, 49999200, 25238080, 9004192, 2225760, 363136, 35200, 1536], pExp := 2, qExp := 0, vExp := 1 },
    { nCoeffs := [-527472, -2881224, -7413300, -11697582, -12428360, -9255414, -4902096, -1840424, -479088, -82272, -8384, -384], pExp := 1, qExp := 3, vExp := 0 },
    { nCoeffs := [-134064, -559944, -1101908, -1343822, -1113628, -646070, -262148, -72688, -13088, -1376, -64], pExp := 1, qExp := 1, vExp := 2 },
    { nCoeffs := [-154368, -840960, -2030720, -2864384, -2615088, -1615472, -684128, -196176, -36464, -3968, -192], pExp := 1, qExp := 1, vExp := 0 },
    { nCoeffs := [52920, 201348, 337002, 325727, 200392, 81391, 21828, 3728, 368, 16], pExp := 0, qExp := 4, vExp := 3 },
    { nCoeffs := [52920, 201348, 337002, 325727, 200392, 81391, 21828, 3728, 368, 16], pExp := 0, qExp := 4, vExp := 1 },
    { nCoeffs := [-10584, -36036, -52986, -43951, -22498, -7279, -1454, -164, -8], pExp := 0, qExp := 2, vExp := 3 },
    { nCoeffs := [-2647944, -14178420, -35810742, -55479583, -57841646, -42231495, -21910654, -8051876, -2050440, -344320, -34304, -1536], pExp := 0, qExp := 2, vExp := 1 },
    { nCoeffs := [-578736, -4010976, -10862764, -16183788, -15106340, -9363572, -3934672, -1112112, -202944, -21632, -1024], pExp := 0, qExp := 0, vExp := 1 }
  ]

private def dualCertPp2Poly : DualCertPoly :=
  dualCertOfTerms dualCertPp2Terms

private def dualCertPq2Poly : DualCertPoly :=
  dualCertOfTerms dualCertPq2Terms

private def dualCertPv2Poly : DualCertPoly :=
  dualCertOfTerms dualCertPv2Terms

private def dualCertPp2 (n p q v : ℝ) : ℝ :=
  dualCertEval dualCertPp2Poly n p q v

private def dualCertPq2 (n p q v : ℝ) : ℝ :=
  dualCertEval dualCertPq2Poly n p q v

private def dualCertPv2 (n p q v : ℝ) : ℝ :=
  dualCertEval dualCertPv2Poly n p q v

private def dualCertDPp2 (n p q v : ℝ) : ℝ :=
  dualCertEval (dualCertPDeriv dualCertPp2Poly) n p q v

private def dualCertDPq2 (n p q v : ℝ) : ℝ :=
  dualCertEval (dualCertQDeriv dualCertPq2Poly) n p q v

private def dualCertDPv2 (n p q v : ℝ) : ℝ :=
  dualCertEval (dualCertVDeriv dualCertPv2Poly) n p q v

/-! The remaining sparse certificate rows.  Row zero is supplied by
`dualCert_sparse_identity_zero`; row two is normalized here with the same
kernel-checked sparse operations. -/

private def millerSpM20 : SparsePoly :=
  (4 * SparsePoly.n + 10) * (SparsePoly.n + 2) ^ 2 *
    (SparsePoly.n + 3) ^ 2 *
    (32 * SparsePoly.n ^ 4 + 302 * SparsePoly.n ^ 3 +
      1037 * SparsePoly.n ^ 2 + 1530 * SparsePoly.n + 813)

private def millerSpM21 : SparsePoly :=
  (SparsePoly.n + 2) ^ 2 *
    (192 * SparsePoly.n ^ 6 + 2984 * SparsePoly.n ^ 5 +
      19116 * SparsePoly.n ^ 4 + 64452 * SparsePoly.n ^ 3 +
      120256 * SparsePoly.n ^ 2 + 117279 * SparsePoly.n + 46476)

private def millerSpM22 : SparsePoly :=
  (SparsePoly.n + 2) ^ 2 *
    (16 * SparsePoly.n ^ 5 + 408 * SparsePoly.n ^ 4 +
      2912 * SparsePoly.n ^ 3 + 8884 * SparsePoly.n ^ 2 +
      12254 * SparsePoly.n + 6240)

private def millerSpPp2 : SparsePoly :=
  SparsePoly.ofDualTerms dualCertPp2Terms

private def millerSpPq2 : SparsePoly :=
  SparsePoly.ofDualTerms dualCertPq2Terms

private def millerSpPv2 : SparsePoly :=
  SparsePoly.ofDualTerms dualCertPv2Terms

private theorem millerSparseEval_spM00
    (n : ℕ) (p q v : ℝ) :
    millerSparseEval spM00 (n : ℝ) p q v =
      (positiveMatrix (n : ℤ) 0 0 : ℝ) := by
  simp only [spM00, millerSparseEval_add, millerSparseEval_mul,
    millerSparseEval_pow, millerSparseEval_ofNat,
    millerSparseEval_n]
  simp [positiveMatrix, Matrix.cons_val_zero]

private theorem millerSparseEval_spM01
    (n : ℕ) (p q v : ℝ) :
    millerSparseEval spM01 (n : ℝ) p q v =
      (positiveMatrix (n : ℤ) 0 1 : ℝ) := by
  simp only [spM01, millerSparseEval_add, millerSparseEval_mul,
    millerSparseEval_pow, millerSparseEval_ofNat,
    millerSparseEval_n]
  simp [positiveMatrix, Matrix.cons_val_one]

private theorem millerSparseEval_spM02
    (n : ℕ) (p q v : ℝ) :
    millerSparseEval spM02 (n : ℝ) p q v =
      (positiveMatrix (n : ℤ) 0 2 : ℝ) := by
  simp only [spM02, millerSparseEval_add, millerSparseEval_mul,
    millerSparseEval_pow, millerSparseEval_ofNat,
    millerSparseEval_n]
  simp [positiveMatrix, Matrix.cons_val_two,
    Matrix.head_cons, Matrix.tail_cons]

private theorem millerDualCert_identity_zero_pointwise
    (n : ℕ) (p q v : ℝ) :
    (4 * (2 * (n : ℝ) + 3) * ((n : ℝ) + 2)) *
        ((positiveMatrix (n : ℤ) 0 0 : ℝ) *
              dualCertNextNum (n : ℝ) 0 p q v +
          (positiveMatrix (n : ℤ) 0 1 : ℝ) *
              dualCertNextNum (n : ℝ) 1 p q v +
          (positiveMatrix (n : ℤ) 0 2 : ℝ) *
              dualCertNextNum (n : ℝ) 2 p q v -
          dualCertLambda (n : ℝ) *
              dualCertCurNum (n : ℝ) 0 p q v) =
      dualCertOpP (n : ℝ) (dualCertPp0 (n : ℝ) p q v)
          (dualCertDPp0 (n : ℝ) p q v) p q v +
        dualCertOpQ (n : ℝ) (dualCertPq0 (n : ℝ) p q v)
          (dualCertDPq0 (n : ℝ) p q v) p q v +
        dualCertOpV (n : ℝ) (dualCertPv0 (n : ℝ) p q v)
          (dualCertDPv0 (n : ℝ) p q v) p q v := by
  have h := congrArg
    (fun P : SparsePoly =>
      millerSparseEval P (n : ℝ) p q v)
    dualCert_sparse_identity_zero
  simp only [millerSparseEval_mul, millerSparseEval_add,
    millerSparseEval_sub, millerSparseEval_spDelta,
    millerSparseEval_spM00, millerSparseEval_spM01,
    millerSparseEval_spM02, millerSparseEval_spNext,
    millerSparseEval_spLambda, millerSparseEval_spCur,
    spPp0, spPq0, spPv0,
    millerSparseEval_spOpP_ofDualTerms,
    millerSparseEval_spOpQ_ofDualTerms,
    millerSparseEval_spOpV_ofDualTerms] at h
  simpa [dualCertPp0, dualCertPq0, dualCertPv0,
    dualCertDPp0, dualCertDPq0, dualCertDPv0,
    dualCertPp0Poly, dualCertPq0Poly, dualCertPv0Poly] using h

set_option maxHeartbeats 0 in
set_option maxRecDepth 1000000 in
private theorem millerDualCert_identity_one_pointwise
    (n : ℕ) (p q v : ℝ) :
    (4 * (2 * (n : ℝ) + 3) * ((n : ℝ) + 2)) *
        ((((n : ℝ) + 2) ^ 2 * ((n : ℝ) + 3) ^ 2 *
              (4 * (n : ℝ) + 10) *
              (48 * (n : ℝ) ^ 3 + 386 * (n : ℝ) ^ 2 +
                1017 * (n : ℝ) + 879)) *
              dualCertNextNum (n : ℝ) 0 p q v +
          (((n : ℝ) + 2) ^ 2 *
              (272 * (n : ℝ) ^ 5 + 3848 * (n : ℝ) ^ 4 +
                21732 * (n : ℝ) ^ 3 + 61184 * (n : ℝ) ^ 2 +
                85761 * (n : ℝ) + 47808)) *
              dualCertNextNum (n : ℝ) 1 p q v +
          (((n : ℝ) + 2) ^ 2 *
              (320 * (n : ℝ) ^ 3 + 2540 * (n : ℝ) ^ 2 +
                6610 * (n : ℝ) + 5640)) *
              dualCertNextNum (n : ℝ) 2 p q v -
          dualCertLambda (n : ℝ) *
              dualCertCurNum (n : ℝ) 1 p q v) =
      dualCertOpP (n : ℝ) (dualCertPp1 (n : ℝ) p q v)
          (dualCertDPp1 (n : ℝ) p q v) p q v +
        dualCertOpQ (n : ℝ) (dualCertPq1 (n : ℝ) p q v)
          (dualCertDPq1 (n : ℝ) p q v) p q v +
        dualCertOpV (n : ℝ) (dualCertPv1 (n : ℝ) p q v)
          (dualCertDPv1 (n : ℝ) p q v) p q v := by
  simp only [dualCertPp1, dualCertPq1, dualCertPv1,
    dualCertDPp1, dualCertDPq1, dualCertDPv1,
    dualCertPp1Poly, dualCertPq1Poly, dualCertPv1Poly,
    millerDualCertEval_ofTerms_eq,
    millerDualCertPDerivEval_ofTerms_eq,
    millerDualCertQDerivEval_ofTerms_eq,
    millerDualCertVDerivEval_ofTerms_eq]
  simp only [millerDualTermsEval, millerDualTermsPDerivEval,
    millerDualTermsQDerivEval, millerDualTermsVDerivEval,
    dualCertPp1Terms, dualCertPq1Terms, dualCertPv1Terms,
    List.map_cons, List.map_nil, List.sum_cons, List.sum_nil,
    millerDualTermEval, millerDualTermPDerivEval,
    millerDualTermQDerivEval, millerDualTermVDerivEval,
    dualCertNextNum, dualCertCurNum, dualCertLambda,
    dualCertOpP, dualCertOpQ, dualCertOpV,
    dualCertD, dualCertSnum,
    Matrix.cons_val_zero, Matrix.cons_val_one,
    Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]
  simp only [millerSpCoeffEval]
  ring

set_option maxHeartbeats 0 in
set_option maxRecDepth 1000000 in
private theorem millerDualCert_identity_two_pointwise
    (n : ℕ) (p q v : ℝ) :
    (4 * (2 * (n : ℝ) + 3) * ((n : ℝ) + 2)) *
        (((4 * (n : ℝ) + 10) * ((n : ℝ) + 2) ^ 2 *
              ((n : ℝ) + 3) ^ 2 *
              (32 * (n : ℝ) ^ 4 + 302 * (n : ℝ) ^ 3 +
                1037 * (n : ℝ) ^ 2 + 1530 * (n : ℝ) + 813)) *
              dualCertNextNum (n : ℝ) 0 p q v +
          (((n : ℝ) + 2) ^ 2 *
              (192 * (n : ℝ) ^ 6 + 2984 * (n : ℝ) ^ 5 +
                19116 * (n : ℝ) ^ 4 + 64452 * (n : ℝ) ^ 3 +
                120256 * (n : ℝ) ^ 2 + 117279 * (n : ℝ) + 46476)) *
              dualCertNextNum (n : ℝ) 1 p q v +
          (((n : ℝ) + 2) ^ 2 *
              (16 * (n : ℝ) ^ 5 + 408 * (n : ℝ) ^ 4 +
                2912 * (n : ℝ) ^ 3 + 8884 * (n : ℝ) ^ 2 +
                12254 * (n : ℝ) + 6240)) *
              dualCertNextNum (n : ℝ) 2 p q v -
          dualCertLambda (n : ℝ) *
              dualCertCurNum (n : ℝ) 2 p q v) =
      dualCertOpP (n : ℝ) (dualCertPp2 (n : ℝ) p q v)
          (dualCertDPp2 (n : ℝ) p q v) p q v +
        dualCertOpQ (n : ℝ) (dualCertPq2 (n : ℝ) p q v)
          (dualCertDPq2 (n : ℝ) p q v) p q v +
        dualCertOpV (n : ℝ) (dualCertPv2 (n : ℝ) p q v)
          (dualCertDPv2 (n : ℝ) p q v) p q v := by
  simp only [dualCertPp2, dualCertPq2, dualCertPv2,
    dualCertDPp2, dualCertDPq2, dualCertDPv2,
    dualCertPp2Poly, dualCertPq2Poly, dualCertPv2Poly,
    millerDualCertEval_ofTerms_eq,
    millerDualCertPDerivEval_ofTerms_eq,
    millerDualCertQDerivEval_ofTerms_eq,
    millerDualCertVDerivEval_ofTerms_eq]
  simp only [millerDualTermsEval, millerDualTermsPDerivEval,
    millerDualTermsQDerivEval, millerDualTermsVDerivEval,
    dualCertPp2Terms, dualCertPq2Terms, dualCertPv2Terms,
    List.map_cons, List.map_nil, List.sum_cons, List.sum_nil,
    millerDualTermEval, millerDualTermPDerivEval,
    millerDualTermQDerivEval, millerDualTermVDerivEval,
    dualCertNextNum, dualCertCurNum, dualCertLambda,
    dualCertOpP, dualCertOpQ, dualCertOpV,
    dualCertD, dualCertSnum,
    Matrix.cons_val_zero, Matrix.cons_val_one,
    Matrix.cons_val_succ, Matrix.cons_val_two,
    Matrix.head_cons, Matrix.tail_cons]
  simp only [millerSpCoeffEval]
  ring

private theorem millerDualCert_identity_one
    (n : ℕ) (p q v : ℝ) :
    (4 * (2 * (n : ℝ) + 3) * ((n : ℝ) + 2)) *
        ((positiveMatrix (n : ℤ) 1 0 : ℝ) *
              dualCertNextNum (n : ℝ) 0 p q v +
          (positiveMatrix (n : ℤ) 1 1 : ℝ) *
              dualCertNextNum (n : ℝ) 1 p q v +
          (positiveMatrix (n : ℤ) 1 2 : ℝ) *
              dualCertNextNum (n : ℝ) 2 p q v -
          dualCertLambda (n : ℝ) *
              dualCertCurNum (n : ℝ) 1 p q v) =
      dualCertOpP (n : ℝ) (dualCertPp1 (n : ℝ) p q v)
          (dualCertDPp1 (n : ℝ) p q v) p q v +
        dualCertOpQ (n : ℝ) (dualCertPq1 (n : ℝ) p q v)
          (dualCertDPq1 (n : ℝ) p q v) p q v +
        dualCertOpV (n : ℝ) (dualCertPv1 (n : ℝ) p q v)
          (dualCertDPv1 (n : ℝ) p q v) p q v := by
  simpa [positiveMatrix, Matrix.cons_val_two,
    Matrix.head_cons, Matrix.tail_cons] using
    millerDualCert_identity_one_pointwise n p q v

private theorem millerDualCert_identity_two
    (n : ℕ) (p q v : ℝ) :
    (4 * (2 * (n : ℝ) + 3) * ((n : ℝ) + 2)) *
        ((positiveMatrix (n : ℤ) 2 0 : ℝ) *
              dualCertNextNum (n : ℝ) 0 p q v +
          (positiveMatrix (n : ℤ) 2 1 : ℝ) *
              dualCertNextNum (n : ℝ) 1 p q v +
          (positiveMatrix (n : ℤ) 2 2 : ℝ) *
              dualCertNextNum (n : ℝ) 2 p q v -
          dualCertLambda (n : ℝ) *
              dualCertCurNum (n : ℝ) 2 p q v) =
      dualCertOpP (n : ℝ) (dualCertPp2 (n : ℝ) p q v)
          (dualCertDPp2 (n : ℝ) p q v) p q v +
        dualCertOpQ (n : ℝ) (dualCertPq2 (n : ℝ) p q v)
          (dualCertDPq2 (n : ℝ) p q v) p q v +
        dualCertOpV (n : ℝ) (dualCertPv2 (n : ℝ) p q v)
          (dualCertDPv2 (n : ℝ) p q v) p q v := by
  simpa [positiveMatrix, Matrix.cons_val_two,
    Matrix.head_cons, Matrix.tail_cons] using
    millerDualCert_identity_two_pointwise n p q v

private theorem millerCube_ae_bounds :
    ∀ᵐ x : ℝ × (ℝ × ℝ) ∂cubeMeasure,
      0 < x.1 ∧ x.1 ≤ 1 ∧
      0 < x.2.1 ∧ x.2.1 ≤ 1 ∧
      0 < x.2.2 ∧ x.2.2 ≤ 1 := by
  have hmem : ∀ᵐ x : ℝ × (ℝ × ℝ) ∂cubeMeasure,
      x ∈ Ioc (0 : ℝ) 1 ×ˢ
        (Ioc (0 : ℝ) 1 ×ˢ Ioc (0 : ℝ) 1) := by
    rw [Measure.ae_prod_mem_iff_ae_ae_mem
      (measurableSet_Ioc.prod
        (measurableSet_Ioc.prod measurableSet_Ioc))]
    filter_upwards [unit_ae_bounds] with p hp
    have hqv :
        ∀ᵐ y : ℝ × ℝ ∂unitMeasure.prod unitMeasure,
          y ∈ Ioc (0 : ℝ) 1 ×ˢ Ioc (0 : ℝ) 1 := by
      rw [Measure.ae_prod_mem_iff_ae_ae_mem
        (measurableSet_Ioc.prod measurableSet_Ioc)]
      filter_upwards [unit_ae_bounds] with q hq
      filter_upwards [unit_ae_bounds] with v hv
      exact ⟨hq, hv⟩
    filter_upwards [hqv] with y hy
    exact ⟨hp, hy⟩
  filter_upwards [hmem] with x hx
  exact ⟨hx.1.1, hx.1.2, hx.2.1.1, hx.2.1.2,
    hx.2.2.1, hx.2.2.2⟩

@[simp] private theorem millerDualCertEval_nextNumPoly
    (n p q v : ℝ) (j : Fin 3) :
    dualCertEval (dualCertNextNumPoly j) n p q v =
      dualCertNextNum n j p q v := by
  fin_cases j <;>
    simp [dualCertNextNumPoly, dualCertNextNum,
      dualCertSnumPoly, dualCertDPoly, dualCertSnum,
      dualCertD, Matrix.cons_val_two]

@[simp] private theorem millerDualCertEval_curNumPoly
    (n p q v : ℝ) (i : Fin 3) :
    dualCertEval (dualCertCurNumPoly i) n p q v =
      dualCertCurNum n i p q v := by
  fin_cases i <;>
    simp [dualCertCurNumPoly, dualCertCurNum,
      dualCertDPoly, dualCertD, Matrix.cons_val_two]

private theorem millerDualCertD_eq_dualD (p q v : ℝ) :
    dualCertD p q v = dualD p q v := rfl

private theorem millerDualCertNext_weight_ae
    (n : ℕ) (j : Fin 3) :
    dualCertWeightedIntegrand n (dualCertNextNumPoly j)
        =ᵐ[cubeMeasure]
      ![rawMomentIntegrand (n + 1) 0 0 0 0,
        fun x => (2 * ((n : ℝ) + 3)) *
          rawMomentIntegrand (n + 1) 0 0 1 1 x,
        fun x => -((n : ℝ) + 3) *
            rawMomentIntegrand (n + 1) 0 0 1 1 x +
          (2 * ((n : ℝ) + 3) * (2 * (n : ℝ) + 7)) *
            rawMomentIntegrand (n + 1) 0 0 2 2 x] j := by
  fin_cases j
  all_goals
    filter_upwards [millerCube_ae_bounds] with x hx
    rcases hx with ⟨hp, hp1, hq, hq1, hv, hv1⟩
    have hD : dualD x.1 x.2.1 x.2.2 ≠ 0 := by
      apply ne_of_gt
      unfold dualD
      positivity
    unfold dualCertWeightedIntegrand rawMomentIntegrand
    simp only [millerDualCertEval_nextNumPoly]
    simp [dualCertNextNum, Matrix.cons_val_two]
    simp only [dualCertSnum]
    rw [millerDualCertD_eq_dualD]
  all_goals
    rw [show 2 * (n + 1) + 6 = (2 * n + 6) + 2 by omega,
      show 2 * (n + 1) + 5 = (2 * n + 5) + 2 by omega,
      show 2 * (n + 1) + 3 = (2 * n + 3) + 2 by omega,
      show 2 * n + 4 + 4 = (2 * n + 6) + 2 by omega,
      pow_add, pow_succ]
    field_simp [hD]
    ring

private theorem millerDualCertCur_weight_ae
    (n : ℕ) (i : Fin 3) :
    dualCertWeightedIntegrand n (dualCertCurNumPoly i)
        =ᵐ[cubeMeasure]
      ![rawMomentIntegrand n 0 0 0 0,
        fun x => (2 * ((n : ℝ) + 2)) *
          rawMomentIntegrand n 0 0 1 1 x,
        fun x => -((n : ℝ) + 2) *
            rawMomentIntegrand n 0 0 1 1 x +
          (2 * ((n : ℝ) + 2) * (2 * (n : ℝ) + 5)) *
            rawMomentIntegrand n 0 0 2 2 x] i := by
  fin_cases i
  all_goals
    filter_upwards [millerCube_ae_bounds] with x hx
    rcases hx with ⟨hp, hp1, hq, hq1, hv, hv1⟩
    have hD : dualD x.1 x.2.1 x.2.2 ≠ 0 := by
      apply ne_of_gt
      unfold dualD
      positivity
    unfold dualCertWeightedIntegrand rawMomentIntegrand
    simp only [millerDualCertEval_curNumPoly]
    simp [dualCertCurNum, Matrix.cons_val_two]
    rw [millerDualCertD_eq_dualD]
  all_goals
    rw [show 2 * n + 4 + 4 = (2 * n + 4) + 4 by
      omega, pow_add]
    field_simp [hD]
    ring

private theorem millerDualCertNext_integral
    (n : ℕ) (j : Fin 3) :
    (∫ x, dualCertWeightedIntegrand n
        (dualCertNextNumPoly j) x ∂cubeMeasure) =
      dualVector (n + 1) j := by
  rw [integral_congr_ae (millerDualCertNext_weight_ae n j)]
  fin_cases j
  · rfl
  · change (∫ x, (2 * ((n : ℝ) + 3)) *
        rawMomentIntegrand (n + 1) 0 0 1 1 x
        ∂cubeMeasure) =
      2 * (((n + 1 : ℕ) : ℝ) + 2) *
        dualMoment (n + 1) 0 0 1 1
    rw [MeasureTheory.integral_const_mul]
    change 2 * ((n : ℝ) + 3) *
      dualMoment (n + 1) 0 0 1 1 = _
    congr 2
    push_cast
    ring
  · change (∫ x, -((n : ℝ) + 3) *
          rawMomentIntegrand (n + 1) 0 0 1 1 x +
        (2 * ((n : ℝ) + 3) * (2 * (n : ℝ) + 7)) *
          rawMomentIntegrand (n + 1) 0 0 2 2 x
          ∂cubeMeasure) =
      -(((n + 1 : ℕ) : ℝ) + 2) *
          dualMoment (n + 1) 0 0 1 1 +
        2 * (((n + 1 : ℕ) : ℝ) + 2) *
          (2 * ((n + 1 : ℕ) : ℝ) + 5) *
            dualMoment (n + 1) 0 0 2 2
    have h1 := rawMomentIntegrand_integrable
      (n + 1) 0 0 1 1 (by omega)
    have h2 := rawMomentIntegrand_integrable
      (n + 1) 0 0 2 2 (by omega)
    rw [MeasureTheory.integral_add
          (h1.const_mul _) (h2.const_mul _),
      MeasureTheory.integral_const_mul,
      MeasureTheory.integral_const_mul]
    change -((n : ℝ) + 3) *
        dualMoment (n + 1) 0 0 1 1 +
      2 * ((n : ℝ) + 3) * (2 * (n : ℝ) + 7) *
        dualMoment (n + 1) 0 0 2 2 = _
    push_cast
    ring

private theorem millerDualCertCur_integral
    (n : ℕ) (i : Fin 3) :
    (∫ x, dualCertWeightedIntegrand n
        (dualCertCurNumPoly i) x ∂cubeMeasure) =
      dualVector n i := by
  rw [integral_congr_ae (millerDualCertCur_weight_ae n i)]
  fin_cases i
  · rfl
  · change (∫ x, (2 * ((n : ℝ) + 2)) *
        rawMomentIntegrand n 0 0 1 1 x ∂cubeMeasure) =
      2 * ((n : ℝ) + 2) * dualMoment n 0 0 1 1
    rw [MeasureTheory.integral_const_mul]
    rfl
  · change (∫ x, -((n : ℝ) + 2) *
          rawMomentIntegrand n 0 0 1 1 x +
        (2 * ((n : ℝ) + 2) * (2 * (n : ℝ) + 5)) *
          rawMomentIntegrand n 0 0 2 2 x
          ∂cubeMeasure) =
      -((n : ℝ) + 2) * dualMoment n 0 0 1 1 +
        2 * ((n : ℝ) + 2) * (2 * (n : ℝ) + 5) *
          dualMoment n 0 0 2 2
    have h1 := rawMomentIntegrand_integrable
      n 0 0 1 1 (by omega)
    have h2 := rawMomentIntegrand_integrable
      n 0 0 2 2 (by omega)
    rw [MeasureTheory.integral_add
          (h1.const_mul _) (h2.const_mul _),
      MeasureTheory.integral_const_mul,
      MeasureTheory.integral_const_mul]
    rfl

private theorem millerDualCert_row_adjoint
    (n : ℕ) (i : Fin 3) (Pp Pq Pv : DualCertPoly)
    (hidentity : ∀ p q v : ℝ,
      (4 * (2 * (n : ℝ) + 3) * ((n : ℝ) + 2)) *
          ((∑ j : Fin 3, (positiveMatrix (n : ℤ) i j : ℝ) *
              dualCertEval (dualCertNextNumPoly j)
                (n : ℝ) p q v) -
            dualCertLambda (n : ℝ) *
              dualCertEval (dualCertCurNumPoly i)
                (n : ℝ) p q v) =
        dualCertEval (dualCertOpPPoly Pp) (n : ℝ) p q v +
          dualCertEval (dualCertOpQPoly Pq) (n : ℝ) p q v +
          dualCertEval (dualCertOpVPoly Pv) (n : ℝ) p q v) :
    ∑ j : Fin 3,
        (positiveMatrix (n : ℤ) i j : ℝ) *
          dualVector (n + 1) j =
      dualCertLambda (n : ℝ) * dualVector n i := by
  let δ : ℝ := 4 * (2 * (n : ℝ) + 3) * ((n : ℝ) + 2)
  let a0 : ℝ := positiveMatrix (n : ℤ) i 0
  let a1 : ℝ := positiveMatrix (n : ℤ) i 1
  let a2 : ℝ := positiveMatrix (n : ℤ) i 2
  let L : ℝ := dualCertLambda (n : ℝ)
  let f0 := dualCertWeightedIntegrand n (dualCertNextNumPoly 0)
  let f1 := dualCertWeightedIntegrand n (dualCertNextNumPoly 1)
  let f2 := dualCertWeightedIntegrand n (dualCertNextNumPoly 2)
  let fc := dualCertWeightedIntegrand n (dualCertCurNumPoly i)
  let fp := dualCertWeightedIntegrand n (dualCertOpPPoly Pp)
  let fq := dualCertWeightedIntegrand n (dualCertOpQPoly Pq)
  let fv := dualCertWeightedIntegrand n (dualCertOpVPoly Pv)
  have hf0 : Integrable f0 cubeMeasure :=
    dualCertWeightedIntegrand_integrable _ _
  have hf1 : Integrable f1 cubeMeasure :=
    dualCertWeightedIntegrand_integrable _ _
  have hf2 : Integrable f2 cubeMeasure :=
    dualCertWeightedIntegrand_integrable _ _
  have hfc : Integrable fc cubeMeasure :=
    dualCertWeightedIntegrand_integrable _ _
  have hfp : Integrable fp cubeMeasure :=
    dualCertWeightedIntegrand_integrable _ _
  have hfq : Integrable fq cubeMeasure :=
    dualCertWeightedIntegrand_integrable _ _
  have hfv : Integrable fv cubeMeasure :=
    dualCertWeightedIntegrand_integrable _ _
  have hpoint : ∀ x : ℝ × (ℝ × ℝ),
      δ * (a0 * f0 x + a1 * f1 x + a2 * f2 x -
        L * fc x) =
        fp x + fq x + fv x := by
    intro x
    have h := hidentity x.1 x.2.1 x.2.2
    simp only [Fin.sum_univ_three] at h
    unfold δ a0 a1 a2 L f0 f1 f2 fc fp fq fv
    unfold dualCertWeightedIntegrand
    linear_combination
      (rawMomentIntegrand n 0 0 0 4 x) * h
  have hsum :
      Integrable
        (fun x => a0 * f0 x + a1 * f1 x + a2 * f2 x)
        cubeMeasure :=
    ((hf0.const_mul a0).add (hf1.const_mul a1)).add
      (hf2.const_mul a2)
  have hsum_integral :
      (∫ x, a0 * f0 x + a1 * f1 x + a2 * f2 x
          ∂cubeMeasure) =
        a0 * (∫ x, f0 x ∂cubeMeasure) +
          a1 * (∫ x, f1 x ∂cubeMeasure) +
          a2 * (∫ x, f2 x ∂cubeMeasure) := by
    calc
      _ = (∫ x, a0 * f0 x + a1 * f1 x ∂cubeMeasure) +
          ∫ x, a2 * f2 x ∂cubeMeasure := by
        exact MeasureTheory.integral_add
          ((hf0.const_mul a0).add (hf1.const_mul a1))
          (hf2.const_mul a2)
      _ = _ := by
        rw [MeasureTheory.integral_add
              (hf0.const_mul a0) (hf1.const_mul a1),
          MeasureTheory.integral_const_mul,
          MeasureTheory.integral_const_mul,
          MeasureTheory.integral_const_mul]
  have hleft :
      (∫ x, δ *
          (a0 * f0 x + a1 * f1 x + a2 * f2 x -
            L * fc x) ∂cubeMeasure) =
        δ * ((a0 * dualVector (n + 1) 0 +
          a1 * dualVector (n + 1) 1 +
          a2 * dualVector (n + 1) 2) -
            L * dualVector n i) := by
    rw [MeasureTheory.integral_const_mul,
      MeasureTheory.integral_sub hsum (hfc.const_mul L),
      hsum_integral, MeasureTheory.integral_const_mul,
      millerDualCertNext_integral,
      millerDualCertNext_integral,
      millerDualCertNext_integral,
      millerDualCertCur_integral]
  have hzero :
      (∫ x, δ *
          (a0 * f0 x + a1 * f1 x + a2 * f2 x -
            L * fc x) ∂cubeMeasure) = 0 := by
    calc
      _ = ∫ x, fp x + fq x + fv x ∂cubeMeasure :=
        integral_congr_ae
          (Filter.Eventually.of_forall hpoint)
      _ = (∫ x, fp x ∂cubeMeasure) +
          (∫ x, fq x ∂cubeMeasure) +
          (∫ x, fv x ∂cubeMeasure) := by
        calc
          _ = (∫ x, fp x + fq x ∂cubeMeasure) +
              ∫ x, fv x ∂cubeMeasure :=
            MeasureTheory.integral_add (hfp.add hfq) hfv
          _ = _ := by
            rw [MeasureTheory.integral_add hfp hfq]
      _ = 0 := by
        rw [dualCertOpP_cube_integral_zero,
          dualCertOpQ_cube_integral_zero,
          dualCertOpV_cube_integral_zero]
        ring
  rw [hleft] at hzero
  have hδ : δ ≠ 0 := by
    unfold δ
    positivity
  have hinner :
      a0 * dualVector (n + 1) 0 +
          a1 * dualVector (n + 1) 1 +
          a2 * dualVector (n + 1) 2 -
        L * dualVector n i = 0 :=
    (mul_eq_zero.mp hzero).resolve_left hδ
  simp only [Fin.sum_univ_three]
  change a0 * dualVector (n + 1) 0 +
      a1 * dualVector (n + 1) 1 +
      a2 * dualVector (n + 1) 2 =
    L * dualVector n i
  linarith

private theorem millerDualVector_adjoint
    (n : ℕ) (i : Fin 3) :
    ∑ j : Fin 3,
        (positiveMatrix (n : ℤ) i j : ℝ) *
          dualVector (n + 1) j =
      dualCertLambda (n : ℝ) * dualVector n i := by
  fin_cases i
  · apply millerDualCert_row_adjoint n 0
      dualCertPp0Poly dualCertPq0Poly dualCertPv0Poly
    intro p q v
    simpa [Fin.sum_univ_three, dualCertPp0,
      dualCertPq0, dualCertPv0, dualCertDPp0,
      dualCertDPq0, dualCertDPv0] using
      millerDualCert_identity_zero_pointwise n p q v
  · apply millerDualCert_row_adjoint n 1
      dualCertPp1Poly dualCertPq1Poly dualCertPv1Poly
    intro p q v
    simpa [Fin.sum_univ_three, dualCertPp1,
      dualCertPq1, dualCertPv1, dualCertDPp1,
      dualCertDPq1, dualCertDPv1] using
      millerDualCert_identity_one n p q v
  · apply millerDualCert_row_adjoint n 2
      dualCertPp2Poly dualCertPq2Poly dualCertPv2Poly
    intro p q v
    simpa [Fin.sum_univ_three, dualCertPp2,
      dualCertPq2, dualCertPv2, dualCertDPp2,
      dualCertDPq2, dualCertDPv2] using
      millerDualCert_identity_two n p q v

/-! ## A rational companion and its backward chart

The Miller argument is carried out in a Krylov companion basis for the
inverse-transpose trajectory.  Its three scalar coefficients admit uniform
rational bounds from the very first step. -/

private def millerD (n : ℕ) : ℝ :=
  3072 * (n : ℝ) ^ 7 + 61824 * (n : ℝ) ^ 6 +
    531184 * (n : ℝ) ^ 5 + 2525522 * (n : ℝ) ^ 4 +
    7175793 * (n : ℝ) ^ 3 + 12183785 * (n : ℝ) ^ 2 +
    11446123 * (n : ℝ) + 4589916

private def millerP (n : ℕ) : ℝ :=
  430080 * (n : ℝ) ^ 11 + 15954432 * (n : ℝ) ^ 10 +
    267573056 * (n : ℝ) ^ 9 + 2677696216 * (n : ℝ) ^ 8 +
    17764567836 * (n : ℝ) ^ 7 + 82030514302 * (n : ℝ) ^ 6 +
    269008975593 * (n : ℝ) ^ 5 + 626475469357 * (n : ℝ) ^ 4 +
    1015303508704 * (n : ℝ) ^ 3 + 1090544343906 * (n : ℝ) ^ 2 +
    698695833276 * (n : ℝ) + 202285840752

private def millerQ (n : ℕ) : ℝ :=
  1720320 * (n : ℝ) ^ 11 + 62195712 * (n : ℝ) ^ 10 +
    1017208064 * (n : ℝ) ^ 9 + 9933863520 * (n : ℝ) ^ 8 +
    64362354608 * (n : ℝ) ^ 7 + 290489793800 * (n : ℝ) ^ 6 +
    931925873076 * (n : ℝ) ^ 5 + 2125092069586 * (n : ℝ) ^ 4 +
    3375535586927 * (n : ℝ) ^ 3 + 3557018693358 * (n : ℝ) ^ 2 +
    2237964273360 * (n : ℝ) + 636907494879

private def millerC0 (n : ℕ) : ℝ :=
  ((n : ℝ) + 2) ^ 2 * (2 * (n : ℝ) + 7) * millerD (n + 1) /
    (((n : ℝ) + 4) * ((n : ℝ) + 5) * (2 * (n : ℝ) + 9) * millerD n)

private def millerC1 (n : ℕ) : ℝ :=
  -millerP n /
    (((n : ℝ) + 4) * ((n : ℝ) + 5) * (2 * (n : ℝ) + 5) *
      (2 * (n : ℝ) + 9) * millerD n)

private def millerC2 (n : ℕ) : ℝ :=
  millerQ n /
    (4 * ((n : ℝ) + 3) * ((n : ℝ) + 5) * (2 * (n : ℝ) + 7) *
      (2 * (n : ℝ) + 9) * millerD n)

private theorem millerD_pos (n : ℕ) : 0 < millerD n := by
  unfold millerD
  positivity

private theorem miller_coeff_bounds (n : ℕ) :
    1 ≤ millerC0 n ∧ millerC0 n ≤ 4 / 3 ∧
      -(49 : ℝ) ≤ millerC1 n ∧ millerC1 n ≤ -35 ∧
      35 ≤ millerC2 n ∧ millerC2 n ≤ 37 := by
  have hn : (0 : ℝ) ≤ n := by positivity
  have hp : ∀ k : ℕ, 0 ≤ (n : ℝ) ^ k := by
    intro k
    positivity
  have hd := millerD_pos n
  have hd0 :
      0 < ((n : ℝ) + 4) * ((n : ℝ) + 5) *
        (2 * (n : ℝ) + 9) * millerD n := by positivity
  have hd1 :
      0 < ((n : ℝ) + 4) * ((n : ℝ) + 5) *
        (2 * (n : ℝ) + 5) * (2 * (n : ℝ) + 9) * millerD n := by
    positivity
  have hd2 :
      0 < 4 * ((n : ℝ) + 3) * ((n : ℝ) + 5) *
        (2 * (n : ℝ) + 7) * (2 * (n : ℝ) + 9) * millerD n := by
    positivity
  constructor
  · rw [millerC0, le_div_iff₀ hd0]
    norm_num [millerD, Nat.cast_add, Nat.cast_one]
    ring_nf
    nlinarith [hp 0, hp 1, hp 2, hp 3, hp 4, hp 5, hp 6, hp 7,
      hp 8, hp 9, hp 10, hp 11]
  constructor
  · rw [millerC0, div_le_iff₀ hd0]
    norm_num [millerD, Nat.cast_add, Nat.cast_one]
    ring_nf
    nlinarith [hp 0, hp 1, hp 2, hp 3, hp 4, hp 5, hp 6, hp 7,
      hp 8, hp 9, hp 10, hp 11]
  constructor
  · rw [millerC1, le_div_iff₀ hd1]
    norm_num [millerP, millerD]
    ring_nf
    nlinarith [hp 0, hp 1, hp 2, hp 3, hp 4, hp 5, hp 6, hp 7,
      hp 8, hp 9, hp 10, hp 11]
  constructor
  · rw [millerC1, div_le_iff₀ hd1]
    norm_num [millerP, millerD]
    ring_nf
    nlinarith [hp 0, hp 1, hp 2, hp 3, hp 4, hp 5, hp 6, hp 7,
      hp 8, hp 9, hp 10, hp 11]
  constructor
  · rw [millerC2, le_div_iff₀ hd2]
    norm_num [millerQ, millerD]
    ring_nf
    nlinarith [hp 0, hp 1, hp 2, hp 3, hp 4, hp 5, hp 6, hp 7,
      hp 8, hp 9, hp 10, hp 11]
  · rw [millerC2, div_le_iff₀ hd2]
    norm_num [millerQ, millerD]
    ring_nf
    nlinarith [hp 0, hp 1, hp 2, hp 3, hp 4, hp 5, hp 6, hp 7,
      hp 8, hp 9, hp 10, hp 11]

private def millerA (n : ℕ) : ℝ := -millerC1 n / millerC0 n
private def millerB (n : ℕ) : ℝ := -millerC2 n / millerC0 n
private def millerR (n : ℕ) : ℝ := 1 / millerC0 n

private theorem miller_pull_coeff_bounds (n : ℕ) :
    105 / 4 ≤ millerA n ∧ millerA n ≤ 49 ∧
      -37 ≤ millerB n ∧ millerB n ≤ -(105 / 4) ∧
      3 / 4 ≤ millerR n ∧ millerR n ≤ 1 := by
  rcases miller_coeff_bounds n with ⟨hc0l, hc0u, hc1l, hc1u, hc2l, hc2u⟩
  have hc0 : 0 < millerC0 n := lt_of_lt_of_le (by norm_num) hc0l
  constructor
  · rw [millerA, le_div_iff₀ hc0]
    nlinarith
  constructor
  · rw [millerA, div_le_iff₀ hc0]
    nlinarith
  constructor
  · rw [millerB, le_div_iff₀ hc0]
    nlinarith
  constructor
  · rw [millerB, div_le_iff₀ hc0]
    nlinarith
  constructor
  · rw [millerR, le_div_iff₀ hc0]
    nlinarith
  · rw [millerR, div_le_iff₀ hc0]
    nlinarith

private abbrev MillerChart := Fin 2 → ℝ
private abbrev MillerVec := Fin 3 → ℝ

private def millerInBox (z : MillerChart) : Prop :=
  0 ≤ z 0 ∧ z 0 ≤ 1 / 16 ∧
    0 ≤ z 1 ∧ z 1 ≤ 1 / 256

private def millerZeroChart : MillerChart := fun _ => 0

private def millerLift (z : MillerChart) : MillerVec :=
  ![1, z 0, z 1]

private def millerDen (n : ℕ) (z : MillerChart) : ℝ :=
  millerA n + millerB n * z 0 + millerR n * z 1

private def millerPull (n : ℕ) (z : MillerChart) : MillerChart :=
  ![1 / millerDen n z, z 0 / millerDen n z]

private def millerDist (x y : MillerChart) : ℝ :=
  max |x 0 - y 0| |x 1 - y 1|

private theorem miller_zero_mem_box : millerInBox millerZeroChart := by
  norm_num [millerInBox, millerZeroChart]

private theorem miller_dist_nonneg (x y : MillerChart) :
    0 ≤ millerDist x y := by
  simp [millerDist]

private theorem miller_coord_le_dist (x y : MillerChart) (i : Fin 2) :
    |x i - y i| ≤ millerDist x y := by
  fin_cases i
  · exact le_max_left _ _
  · exact le_max_right _ _

private theorem miller_box_diameter {x y : MillerChart}
    (hx : millerInBox x) (hy : millerInBox y) :
    millerDist x y ≤ 1 := by
  rcases hx with ⟨hx0, hx1, hx2, hx3⟩
  rcases hy with ⟨hy0, hy1, hy2, hy3⟩
  rw [millerDist, max_le_iff]
  constructor <;> rw [abs_le] <;> constructor <;> nlinarith

private theorem miller_den_bounds (n : ℕ) {z : MillerChart}
    (hz : millerInBox z) :
    23 ≤ millerDen n z ∧ millerDen n z ≤ 50 := by
  rcases miller_pull_coeff_bounds n with
    ⟨ha0, ha1, hb0, hb1, hr0, hr1⟩
  rcases hz with ⟨hz0, hz1, hz2, hz3⟩
  have hbz_lower :
      -(37 : ℝ) / 16 ≤ millerB n * z 0 := by
    have hprod : 0 ≤ (millerB n + 37) * z 0 :=
      mul_nonneg (by linarith) hz0
    have hzgap : 0 ≤ 37 * (1 / 16 - z 0) :=
      mul_nonneg (by norm_num) (by linarith)
    nlinarith
  have hrnonneg : 0 ≤ millerR n := le_trans (by norm_num) hr0
  have hrz_lower : 0 ≤ millerR n * z 1 :=
    mul_nonneg hrnonneg hz2
  have hbz_upper : millerB n * z 0 ≤ 0 :=
    mul_nonpos_of_nonpos_of_nonneg (by linarith) hz0
  have hrz_upper : millerR n * z 1 ≤ 1 / 256 := by
    have hprod := mul_le_mul hr1 hz3 hz2 (by linarith : 0 ≤ (1 : ℝ))
    nlinarith
  unfold millerDen
  constructor <;> nlinarith

private theorem miller_pull_mem_box (n : ℕ) {z : MillerChart}
    (hz : millerInBox z) :
    millerInBox (millerPull n z) := by
  rcases miller_den_bounds n hz with ⟨hd0, hd1⟩
  have hd : 0 < millerDen n z := lt_of_lt_of_le (by norm_num) hd0
  rcases hz with ⟨hz0, hz1, hz2, hz3⟩
  simp only [millerInBox, millerPull, Matrix.cons_val_zero,
    Matrix.cons_val_one]
  constructor
  · positivity
  constructor
  · rw [div_le_iff₀ hd]
    norm_num
    linarith
  constructor
  · exact div_nonneg hz0 hd.le
  · rw [div_le_iff₀ hd]
    have hgap : 0 ≤ (1 / 256 : ℝ) * (millerDen n z - 23) :=
      mul_nonneg (by norm_num) (by linarith)
    nlinarith

private theorem miller_den_sub (n : ℕ) (x y : MillerChart) :
    millerDen n x - millerDen n y =
      millerB n * (x 0 - y 0) + millerR n * (x 1 - y 1) := by
  simp only [millerDen]
  ring

private theorem miller_pull_contract (n : ℕ) {x y : MillerChart}
    (hx : millerInBox x) (hy : millerInBox y) :
    millerDist (millerPull n x) (millerPull n y) ≤
      (1 / 8 : ℝ) * millerDist x y := by
  let δ := millerDist x y
  let dx := millerDen n x
  let dy := millerDen n y
  have hδ : 0 ≤ δ := miller_dist_nonneg x y
  have hx0 := miller_coord_le_dist x y (0 : Fin 2)
  have hx1 := miller_coord_le_dist x y (1 : Fin 2)
  change |x 0 - y 0| ≤ δ at hx0
  change |x 1 - y 1| ≤ δ at hx1
  rcases miller_den_bounds n hx with ⟨hdx0, hdx1⟩
  rcases miller_den_bounds n hy with ⟨hdy0, hdy1⟩
  change 23 ≤ dx at hdx0
  change dx ≤ 50 at hdx1
  change 23 ≤ dy at hdy0
  change dy ≤ 50 at hdy1
  have hdx : 0 < dx := lt_of_lt_of_le (by norm_num) hdx0
  have hdy : 0 < dy := lt_of_lt_of_le (by norm_num) hdy0
  have hprod : 529 ≤ dx * dy := by
    nlinarith [mul_nonneg (sub_nonneg.mpr hdx0) (sub_nonneg.mpr hdy0)]
  have hprodPos : 0 < dx * dy := mul_pos hdx hdy
  rcases miller_pull_coeff_bounds n with
    ⟨ha0, ha1, hb0, hb1, hr0, hr1⟩
  have hbnonpos : millerB n ≤ 0 := by linarith
  have hbabs : |millerB n| ≤ 37 := by
    rw [abs_of_nonpos hbnonpos]
    linarith
  have hrabs : |millerR n| ≤ 1 := by
    rw [abs_of_nonneg (le_trans (by norm_num) hr0)]
    exact hr1
  have hddiff : |dx - dy| ≤ 38 * δ := by
    change |millerDen n x - millerDen n y| ≤ _
    rw [miller_den_sub]
    calc
      |millerB n * (x 0 - y 0) + millerR n * (x 1 - y 1)| ≤
          |millerB n * (x 0 - y 0)| +
            |millerR n * (x 1 - y 1)| := abs_add_le _ _
      _ = |millerB n| * |x 0 - y 0| +
          |millerR n| * |x 1 - y 1| := by rw [abs_mul, abs_mul]
      _ ≤ 37 * δ + 1 * δ := by
        gcongr
      _ = 38 * δ := by ring
  have hratio : |dy - dx| / (dx * dy) ≤ δ := by
    have hab : |dy - dx| = |dx - dy| := abs_sub_comm _ _
    rw [hab, div_le_iff₀ hprodPos]
    have hgap : 0 ≤ δ * (dx * dy - 529) :=
      mul_nonneg hδ (sub_nonneg.mpr hprod)
    nlinarith
  have hfirst :
      |1 / dx - 1 / dy| ≤ (1 / 8 : ℝ) * δ := by
    rw [show 1 / dx - 1 / dy = (dy - dx) / (dx * dy) by
      field_simp [hdx.ne', hdy.ne']]
    have h38 : |dy - dx| ≤ 38 * δ := by
      simpa [abs_sub_comm] using hddiff
    rw [abs_div, abs_of_pos hprodPos]
    rw [div_le_iff₀ hprodPos]
    have hgap : 0 ≤ δ * (dx * dy - 529) :=
      mul_nonneg hδ (sub_nonneg.mpr hprod)
    nlinarith
  have hyabs : |y 0| ≤ 1 / 16 := by
    rw [abs_of_nonneg hy.1]
    exact hy.2.1
  have hsecond :
      |x 0 / dx - y 0 / dy| ≤ (1 / 8 : ℝ) * δ := by
    rw [show x 0 / dx - y 0 / dy =
        (x 0 - y 0) / dx + y 0 * (dy - dx) / (dx * dy) by
      field_simp [hdx.ne', hdy.ne']
      ring]
    calc
      |(x 0 - y 0) / dx + y 0 * (dy - dx) / (dx * dy)| ≤
          |(x 0 - y 0) / dx| +
            |y 0 * (dy - dx) / (dx * dy)| := abs_add_le _ _
      _ = |x 0 - y 0| / dx +
          |y 0| * (|dy - dx| / (dx * dy)) := by
        rw [abs_div, abs_of_pos hdx, abs_div, abs_mul, abs_of_pos hprodPos]
        ring
      _ ≤ δ / 16 + (1 / 16) * δ := by
        have hfirstTerm : |x 0 - y 0| / dx ≤ δ / 16 := by
          rw [div_le_div_iff₀ hdx (by norm_num : (0 : ℝ) < 16)]
          have hgap : 0 ≤ δ * (dx - 23) :=
            mul_nonneg hδ (sub_nonneg.mpr hdx0)
          nlinarith
        gcongr
      _ = (1 / 8 : ℝ) * δ := by ring
  simp only [millerDist, millerPull, Matrix.cons_val_zero,
    Matrix.cons_val_one]
  exact max_le hfirst hsecond

/-! ### The pullback graph

`millerPullN k h z` pulls a terminal chart at time `k + h` back to
time `k`.  The invariant box and the uniform `1 / 8` contraction produce
the unique infinite pullback graph. -/

private def millerPullN (k : ℕ) : ℕ → MillerChart → MillerChart
  | 0, z => z
  | h + 1, z => millerPull k (millerPullN (k + 1) h z)

private theorem millerPullN_mem_box (k h : ℕ) {z : MillerChart}
    (hz : millerInBox z) :
    millerInBox (millerPullN k h z) := by
  induction h generalizing k with
  | zero => simpa [millerPullN] using hz
  | succ h ih =>
      simpa [millerPullN] using
        miller_pull_mem_box k (ih (k := k + 1))

private theorem millerPullN_add (k h t : ℕ) (z : MillerChart) :
    millerPullN k (h + t) z =
      millerPullN k h (millerPullN (k + h) t z) := by
  induction h generalizing k with
  | zero => simp [millerPullN]
  | succ h ih =>
      simp only [Nat.succ_add, millerPullN]
      rw [ih (k := k + 1)]
      rw [show (k + 1) + h = k + (h + 1) by omega]

private theorem millerPullN_contract (k h : ℕ) {x y : MillerChart}
    (hx : millerInBox x) (hy : millerInBox y) :
    millerDist (millerPullN k h x) (millerPullN k h y) ≤
      (1 / 8 : ℝ) ^ h * millerDist x y := by
  induction h generalizing k with
  | zero =>
      simp [millerPullN]
  | succ h ih =>
      calc
        millerDist (millerPullN k (h + 1) x)
            (millerPullN k (h + 1) y) =
            millerDist
              (millerPull k (millerPullN (k + 1) h x))
              (millerPull k (millerPullN (k + 1) h y)) := by
                rfl
        _ ≤ (1 / 8 : ℝ) *
              millerDist (millerPullN (k + 1) h x)
                (millerPullN (k + 1) h y) :=
            miller_pull_contract k
              (millerPullN_mem_box (k + 1) h hx)
              (millerPullN_mem_box (k + 1) h hy)
        _ ≤ (1 / 8 : ℝ) *
              ((1 / 8 : ℝ) ^ h * millerDist x y) :=
            mul_le_mul_of_nonneg_left (ih (k := k + 1))
              (by norm_num)
        _ = (1 / 8 : ℝ) ^ (h + 1) * millerDist x y := by
            rw [pow_succ]
            ring

private def millerApprox (k h : ℕ) : MillerChart :=
  millerPullN k h millerZeroChart

private theorem millerApprox_mem_box (k h : ℕ) :
    millerInBox (millerApprox k h) :=
  millerPullN_mem_box k h miller_zero_mem_box

private theorem millerApprox_add_dist (k h t : ℕ) :
    millerDist (millerApprox k h) (millerApprox k (h + t)) ≤
      (1 / 8 : ℝ) ^ h := by
  change millerDist
    (millerPullN k h millerZeroChart)
    (millerPullN k (h + t) millerZeroChart) ≤ _
  rw [millerPullN_add]
  calc
    millerDist (millerPullN k h millerZeroChart)
        (millerPullN k h
          (millerPullN (k + h) t millerZeroChart)) ≤
        (1 / 8 : ℝ) ^ h *
          millerDist millerZeroChart
            (millerPullN (k + h) t millerZeroChart) :=
      millerPullN_contract k h miller_zero_mem_box
        (millerPullN_mem_box (k + h) t miller_zero_mem_box)
    _ ≤ (1 / 8 : ℝ) ^ h * 1 :=
      mul_le_mul_of_nonneg_left
        (miller_box_diameter miller_zero_mem_box
          (millerPullN_mem_box (k + h) t miller_zero_mem_box))
        (pow_nonneg (by norm_num) h)
    _ = (1 / 8 : ℝ) ^ h := by ring

private theorem miller_pow_tendsto :
    Tendsto (fun h : ℕ => (1 / 8 : ℝ) ^ h) atTop (nhds 0) :=
  tendsto_pow_atTop_nhds_zero_of_lt_one (by norm_num) (by norm_num)

private theorem millerApprox_coord_cauchy (k : ℕ) (i : Fin 2) :
    CauchySeq (fun h => millerApprox k h i) := by
  rw [Metric.cauchySeq_iff']
  intro ε hε
  have hp :
      ∀ᶠ h : ℕ in atTop, (1 / 8 : ℝ) ^ h < ε :=
    miller_pow_tendsto.eventually (gt_mem_nhds hε)
  rcases eventually_atTop.1 hp with ⟨N, hN⟩
  refine ⟨N, fun n hn => ?_⟩
  obtain ⟨t, rfl⟩ := Nat.exists_eq_add_of_le hn
  rw [Real.dist_eq]
  calc
    |millerApprox k (N + t) i - millerApprox k N i| =
        |millerApprox k N i - millerApprox k (N + t) i| :=
      abs_sub_comm _ _
    _ ≤ millerDist (millerApprox k N) (millerApprox k (N + t)) :=
      miller_coord_le_dist _ _ i
    _ ≤ (1 / 8 : ℝ) ^ N := millerApprox_add_dist k N t
    _ < ε := hN N le_rfl

private def millerGraph (k : ℕ) : MillerChart :=
  fun i => Classical.choose
    (cauchySeq_tendsto_of_complete (millerApprox_coord_cauchy k i))

private theorem millerApprox_coord_tendsto (k : ℕ) (i : Fin 2) :
    Tendsto (fun h => millerApprox k h i) atTop
      (nhds (millerGraph k i)) :=
  Classical.choose_spec
    (cauchySeq_tendsto_of_complete (millerApprox_coord_cauchy k i))

private theorem millerApprox_tendsto (k : ℕ) :
    Tendsto (millerApprox k) atTop (nhds (millerGraph k)) := by
  rw [tendsto_pi_nhds]
  exact millerApprox_coord_tendsto k

private theorem millerGraph_mem_box (k : ℕ) :
    millerInBox (millerGraph k) := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · exact ge_of_tendsto (millerApprox_coord_tendsto k 0)
      (Eventually.of_forall fun h => (millerApprox_mem_box k h).1)
  · exact le_of_tendsto (millerApprox_coord_tendsto k 0)
      (Eventually.of_forall fun h => (millerApprox_mem_box k h).2.1)
  · exact ge_of_tendsto (millerApprox_coord_tendsto k 1)
      (Eventually.of_forall fun h => (millerApprox_mem_box k h).2.2.1)
  · exact le_of_tendsto (millerApprox_coord_tendsto k 1)
      (Eventually.of_forall fun h => (millerApprox_mem_box k h).2.2.2)

private theorem millerDen_tendsto (n : ℕ) {f : ℕ → MillerChart}
    {z : MillerChart} (hf : Tendsto f atTop (nhds z)) :
    Tendsto (fun h => millerDen n (f h)) atTop (nhds (millerDen n z)) := by
  have h0 := tendsto_pi_nhds.mp hf 0
  have h1 := tendsto_pi_nhds.mp hf 1
  have ha :
      Tendsto (fun _ : ℕ => millerA n) atTop (nhds (millerA n)) :=
    tendsto_const_nhds
  have hb :
      Tendsto (fun _ : ℕ => millerB n) atTop (nhds (millerB n)) :=
    tendsto_const_nhds
  have hr :
      Tendsto (fun _ : ℕ => millerR n) atTop (nhds (millerR n)) :=
    tendsto_const_nhds
  simpa only [millerDen] using
    (ha.add (hb.mul h0)).add (hr.mul h1)

private theorem millerPull_tendsto (n : ℕ) {f : ℕ → MillerChart}
    {z : MillerChart} (hf : Tendsto f atTop (nhds z))
    (hz : millerInBox z) :
    Tendsto (fun h => millerPull n (f h)) atTop
      (nhds (millerPull n z)) := by
  have h0 := tendsto_pi_nhds.mp hf 0
  have hden := millerDen_tendsto n hf
  have hden_ne : millerDen n z ≠ 0 := by
    have := (miller_den_bounds n hz).1
    linarith
  rw [tendsto_pi_nhds]
  intro i
  fin_cases i
  · have hone :
        Tendsto (fun _ : ℕ => (1 : ℝ)) atTop (nhds 1) :=
        tendsto_const_nhds
    simpa only [millerPull, Matrix.cons_val_zero] using
      hone.div hden hden_ne
  · simpa only [millerPull, Matrix.cons_val_one] using
      h0.div hden hden_ne

private theorem millerGraph_eq_pull (k : ℕ) :
    millerGraph k = millerPull k (millerGraph (k + 1)) := by
  have hleft :
      Tendsto (fun h => millerApprox k (h + 1)) atTop
        (nhds (millerGraph k)) :=
    (millerApprox_tendsto k).comp (tendsto_add_atTop_nat 1)
  have hright :
      Tendsto (fun h => millerPull k (millerApprox (k + 1) h)) atTop
        (nhds (millerPull k (millerGraph (k + 1)))) :=
    millerPull_tendsto k (millerApprox_tendsto (k + 1))
      (millerGraph_mem_box (k + 1))
  exact tendsto_nhds_unique hleft (by
    simpa [millerApprox, millerPullN] using hright)

private theorem millerPullN_moving_tendsto (k : ℕ)
    (terminal : ℕ → MillerChart)
    (hterminal : ∀ h, millerInBox (terminal h)) :
    Tendsto (fun h => millerPullN k h (terminal h)) atTop
      (nhds (millerGraph k)) := by
  rw [tendsto_pi_nhds]
  intro i
  have herr (h : ℕ) :
      |millerPullN k h (terminal h) i - millerApprox k h i| ≤
        (1 / 8 : ℝ) ^ h := by
    calc
      |millerPullN k h (terminal h) i - millerApprox k h i| ≤
          millerDist (millerPullN k h (terminal h))
            (millerApprox k h) :=
        miller_coord_le_dist _ _ i
      _ ≤ (1 / 8 : ℝ) ^ h *
          millerDist (terminal h) millerZeroChart := by
        simpa only [millerApprox] using
          millerPullN_contract k h (hterminal h) miller_zero_mem_box
      _ ≤ (1 / 8 : ℝ) ^ h * 1 :=
        mul_le_mul_of_nonneg_left
          (miller_box_diameter (hterminal h) miller_zero_mem_box)
          (pow_nonneg (by norm_num) h)
      _ = (1 / 8 : ℝ) ^ h := by ring
  have hlower (h : ℕ) :
      -((1 / 8 : ℝ) ^ h) ≤
        millerPullN k h (terminal h) i - millerApprox k h i :=
    (abs_le.mp (herr h)).1
  have hupper (h : ℕ) :
      millerPullN k h (terminal h) i - millerApprox k h i ≤
        (1 / 8 : ℝ) ^ h :=
    (abs_le.mp (herr h)).2
  have hdiff :
      Tendsto
        (fun h =>
          millerPullN k h (terminal h) i - millerApprox k h i)
        atTop (nhds 0) :=
    tendsto_of_tendsto_of_tendsto_of_le_of_le
      (by simpa using miller_pow_tendsto.neg)
      miller_pow_tendsto hlower hupper
  convert hdiff.add (millerApprox_coord_tendsto k i) using 1
  · ext h
    ring
  · ring

/-! ### Transverse expansion

The forward companion step expands every direction transverse to the
pullback graph.  This is the elementary Miller uniqueness estimate. -/

private def millerStep (n : ℕ) (x : MillerVec) : MillerVec :=
  ![x 1, x 2,
    millerC0 n * x 0 + millerC1 n * x 1 + millerC2 n * x 2]

private def millerDefect (z : MillerChart) (x : MillerVec) : MillerChart :=
  ![x 1 - z 0 * x 0, x 2 - z 1 * x 0]

private def millerDefectNorm (z : MillerChart) (x : MillerVec) : ℝ :=
  max |millerDefect z x 0| |millerDefect z x 1|

private theorem millerDefectNorm_nonneg (z : MillerChart) (x : MillerVec) :
    0 ≤ millerDefectNorm z x := by
  simp [millerDefectNorm]

private theorem miller_den_cross (n : ℕ) (z : MillerChart) :
    millerC0 n * millerDen n z =
      -millerC1 n + z 1 - millerC2 n * z 0 := by
  have hc0 : millerC0 n ≠ 0 := by
    have := (miller_coeff_bounds n).1
    linarith
  simp only [millerDen, millerA, millerB, millerR]
  field_simp [hc0]
  ring

private theorem miller_pull_graph_relations (n : ℕ) {z : MillerChart}
    (hz : millerInBox z) :
    millerPull n z 1 = millerPull n z 0 * z 0 ∧
      millerC0 n + (millerC1 n - z 1) * millerPull n z 0 +
        millerC2 n * millerPull n z 1 = 0 := by
  have hd : 0 < millerDen n z :=
    lt_of_lt_of_le (by norm_num) (miller_den_bounds n hz).1
  have hcross := miller_den_cross n z
  constructor
  · simp only [millerPull, Matrix.cons_val_zero, Matrix.cons_val_one]
    field_simp [hd.ne']
  · simp only [millerPull, Matrix.cons_val_zero, Matrix.cons_val_one]
    field_simp [hd.ne']
    nlinarith

private theorem millerGraph_relations (n : ℕ) :
    millerGraph n 1 = millerGraph n 0 * millerGraph (n + 1) 0 ∧
      millerC0 n +
          (millerC1 n - millerGraph (n + 1) 1) * millerGraph n 0 +
        millerC2 n * millerGraph n 1 = 0 := by
  rw [millerGraph_eq_pull n]
  exact miller_pull_graph_relations n (millerGraph_mem_box (n + 1))

private theorem miller_defect_step (n : ℕ) (x : MillerVec) :
    millerDefect (millerGraph (n + 1)) (millerStep n x) 0 =
        -millerGraph (n + 1) 0 *
            millerDefect (millerGraph n) x 0 +
          millerDefect (millerGraph n) x 1 ∧
      millerDefect (millerGraph (n + 1)) (millerStep n x) 1 =
        (millerC1 n - millerGraph (n + 1) 1) *
            millerDefect (millerGraph n) x 0 +
          millerC2 n * millerDefect (millerGraph n) x 1 := by
  rcases millerGraph_relations n with ⟨hmul, hcancel⟩
  constructor
  · simp only [millerDefect, millerStep, Matrix.cons_val_zero,
      Matrix.cons_val_one, Matrix.cons_val_two]
    rw [hmul]
    ring
  · change
      millerC0 n * x 0 + millerC1 n * x 1 + millerC2 n * x 2 -
          millerGraph (n + 1) 1 * x 1 =
        (millerC1 n - millerGraph (n + 1) 1) *
            (x 1 - millerGraph n 0 * x 0) +
          millerC2 n * (x 2 - millerGraph n 1 * x 0)
    linear_combination x 0 * hcancel

private theorem miller_defect_sup_lower
    (c1 c2 un vn q0 q1 Q0 Q1 : ℝ)
    (hc1_lo : -49 ≤ c1) (hc1_hi : c1 ≤ -35)
    (hc2_lo : 35 ≤ c2) (hc2_hi : c2 ≤ 37)
    (hun_lo : 0 ≤ un) (hun_hi : un ≤ (1 : ℝ) / 16)
    (hvn_lo : 0 ≤ vn) (hvn_hi : vn ≤ (1 : ℝ) / 256)
    (hQ0 : Q0 = -un * q0 + q1)
    (hQ1 : Q1 = (c1 - vn) * q0 + c2 * q1) :
    (1 / 2 : ℝ) * max |q0| |q1| ≤ max |Q0| |Q1| := by
  let Δ : ℝ := -c1 + vn - c2 * un
  let K : ℝ := max |Q0| |Q1|
  have hc2_nonneg : 0 ≤ c2 := by linarith
  have hvc_nonneg : 0 ≤ vn - c1 := by linarith
  have hcu : c2 * un ≤ (37 : ℝ) / 16 := by
    calc
      c2 * un ≤ 37 * un :=
        mul_le_mul_of_nonneg_right hc2_hi hun_lo
      _ ≤ 37 * ((1 : ℝ) / 16) :=
        mul_le_mul_of_nonneg_left hun_hi (by norm_num)
      _ = (37 : ℝ) / 16 := by ring
  have hdet : (523 : ℝ) / 16 ≤ Δ := by
    dsimp [Δ]
    nlinarith
  have h32 : (32 : ℝ) ≤ Δ := by
    have : (32 : ℝ) ≤ 523 / 16 := by norm_num
    exact this.trans hdet
  have hΔ_nonneg : 0 ≤ Δ := by linarith
  have hQ0K : |Q0| ≤ K := le_max_left _ _
  have hQ1K : |Q1| ≤ K := le_max_right _ _
  have hK_nonneg : 0 ≤ K := (abs_nonneg Q0).trans hQ0K
  have hinv0 : Δ * q0 = c2 * Q0 - Q1 := by
    dsimp [Δ]
    rw [hQ0, hQ1]
    ring
  have hinv1 : Δ * q1 = (vn - c1) * Q0 - un * Q1 := by
    dsimp [Δ]
    rw [hQ0, hQ1]
    ring
  have habs0 : Δ * |q0| = |c2 * Q0 - Q1| := by
    have h := congrArg abs hinv0
    simpa [abs_mul, abs_of_nonneg hΔ_nonneg] using h
  have habs1 : Δ * |q1| = |(vn - c1) * Q0 - un * Q1| := by
    have h := congrArg abs hinv1
    simpa [abs_mul, abs_of_nonneg hΔ_nonneg] using h
  have hc2Q : c2 * |Q0| ≤ 37 * K :=
    mul_le_mul hc2_hi hQ0K (abs_nonneg Q0) (by norm_num)
  have hrow0 : Δ * |q0| ≤ 38 * K := by
    calc
      Δ * |q0| = |c2 * Q0 - Q1| := habs0
      _ ≤ |c2 * Q0| + |Q1| := abs_sub _ _
      _ = c2 * |Q0| + |Q1| := by
        rw [abs_mul, abs_of_nonneg hc2_nonneg]
      _ ≤ 37 * K + K := add_le_add hc2Q hQ1K
      _ = 38 * K := by ring
  have hvc_le : vn - c1 ≤ 50 := by nlinarith
  have hun_le_one : un ≤ 1 := by nlinarith
  have hvcQ : (vn - c1) * |Q0| ≤ 50 * K :=
    mul_le_mul hvc_le hQ0K (abs_nonneg Q0) (by norm_num)
  have hunQ : un * |Q1| ≤ 1 * K :=
    mul_le_mul hun_le_one hQ1K (abs_nonneg Q1) (by norm_num)
  have hrow1 : Δ * |q1| ≤ 51 * K := by
    calc
      Δ * |q1| = |(vn - c1) * Q0 - un * Q1| := habs1
      _ ≤ |(vn - c1) * Q0| + |un * Q1| := abs_sub _ _
      _ = (vn - c1) * |Q0| + un * |Q1| := by
        rw [abs_mul, abs_of_nonneg hvc_nonneg,
          abs_mul, abs_of_nonneg hun_lo]
      _ ≤ 50 * K + 1 * K := add_le_add hvcQ hunQ
      _ = 51 * K := by ring
  have hq0_scaled : 32 * |q0| ≤ Δ * |q0| :=
    mul_le_mul_of_nonneg_right h32 (abs_nonneg q0)
  have hq1_scaled : 32 * |q1| ≤ Δ * |q1| :=
    mul_le_mul_of_nonneg_right h32 (abs_nonneg q1)
  have h38 : 38 * K ≤ 64 * K :=
    mul_le_mul_of_nonneg_right (by norm_num) hK_nonneg
  have h51 : 51 * K ≤ 64 * K :=
    mul_le_mul_of_nonneg_right (by norm_num) hK_nonneg
  have hq0_two : |q0| ≤ 2 * K := by
    have h : 32 * |q0| ≤ 64 * K :=
      hq0_scaled.trans (hrow0.trans h38)
    nlinarith
  have hq1_two : |q1| ≤ 2 * K := by
    have h : 32 * |q1| ≤ 64 * K :=
      hq1_scaled.trans (hrow1.trans h51)
    nlinarith
  have hmax : max |q0| |q1| ≤ 2 * K :=
    max_le hq0_two hq1_two
  calc
    (1 / 2 : ℝ) * max |q0| |q1| ≤
        (1 / 2 : ℝ) * (2 * K) :=
      mul_le_mul_of_nonneg_left hmax (by norm_num)
    _ = K := by ring
    _ = max |Q0| |Q1| := by rfl

private theorem millerDefectNorm_step_lower (n : ℕ) (x : MillerVec) :
    (1 / 2 : ℝ) * millerDefectNorm (millerGraph n) x ≤
      millerDefectNorm (millerGraph (n + 1)) (millerStep n x) := by
  rcases miller_coeff_bounds n with
    ⟨hc0l, hc0u, hc1l, hc1u, hc2l, hc2u⟩
  rcases millerGraph_mem_box (n + 1) with
    ⟨hu0, hu1, hv0, hv1⟩
  rcases miller_defect_step n x with ⟨hQ0, hQ1⟩
  exact miller_defect_sup_lower
    (millerC1 n) (millerC2 n)
    (millerGraph (n + 1) 0) (millerGraph (n + 1) 1)
    (millerDefect (millerGraph n) x 0)
    (millerDefect (millerGraph n) x 1)
    (millerDefect (millerGraph (n + 1)) (millerStep n x) 0)
    (millerDefect (millerGraph (n + 1)) (millerStep n x) 1)
    hc1l hc1u hc2l hc2u hu0 hu1 hv0 hv1 hQ0 hQ1

private theorem millerDefectNorm_initial_le (x : ℕ → MillerVec)
    (hstep : ∀ n, x (n + 1) = millerStep n (x n)) (n : ℕ) :
    millerDefectNorm (millerGraph 0) (x 0) ≤
      (2 : ℝ) ^ n * millerDefectNorm (millerGraph n) (x n) := by
  induction n with
  | zero => simp
  | succ n ih =>
      have hs := millerDefectNorm_step_lower n (x n)
      rw [← hstep n] at hs
      have hdouble :
          millerDefectNorm (millerGraph n) (x n) ≤
            2 * millerDefectNorm (millerGraph (n + 1)) (x (n + 1)) := by
        nlinarith
      calc
        millerDefectNorm (millerGraph 0) (x 0) ≤
            (2 : ℝ) ^ n * millerDefectNorm (millerGraph n) (x n) := ih
        _ ≤ (2 : ℝ) ^ n *
              (2 * millerDefectNorm (millerGraph (n + 1)) (x (n + 1))) :=
          mul_le_mul_of_nonneg_left hdouble (pow_nonneg (by norm_num) n)
        _ = (2 : ℝ) ^ (n + 1) *
              millerDefectNorm (millerGraph (n + 1)) (x (n + 1)) := by
          rw [pow_succ]
          ring

private theorem miller_fast_defect_zero (x : ℕ → MillerVec)
    (hstep : ∀ n, x (n + 1) = millerStep n (x n))
    (hfast :
      Tendsto
        (fun n => (2 : ℝ) ^ n *
          millerDefectNorm (millerGraph n) (x n))
        atTop (nhds 0)) :
    millerDefectNorm (millerGraph 0) (x 0) = 0 := by
  apply le_antisymm
  · exact ge_of_tendsto hfast
      (Eventually.of_forall fun n =>
        millerDefectNorm_initial_le x hstep n)
  · exact millerDefectNorm_nonneg _ _

private theorem miller_fast_eq_graph (x : ℕ → MillerVec)
    (hstep : ∀ n, x (n + 1) = millerStep n (x n))
    (hfast :
      Tendsto
        (fun n => (2 : ℝ) ^ n *
          millerDefectNorm (millerGraph n) (x n))
        atTop (nhds 0)) :
    x 0 = x 0 0 • millerLift (millerGraph 0) := by
  have hnorm := miller_fast_defect_zero x hstep hfast
  have habs0 :
      |millerDefect (millerGraph 0) (x 0) 0| = 0 := by
    apply le_antisymm
    · exact (le_max_left _ _).trans_eq hnorm
    · exact abs_nonneg _
  have habs1 :
      |millerDefect (millerGraph 0) (x 0) 1| = 0 := by
    apply le_antisymm
    · exact (le_max_right _ _).trans_eq hnorm
    · exact abs_nonneg _
  have hq0 := abs_eq_zero.mp habs0
  have hq1 := abs_eq_zero.mp habs1
  funext i
  fin_cases i
  · simp [millerLift]
  · change x 0 1 = x 0 0 * millerGraph 0 0
    change x 0 1 - millerGraph 0 0 * x 0 0 = 0 at hq0
    linarith
  · change x 0 2 = x 0 0 * millerGraph 0 1
    change x 0 2 - millerGraph 0 1 * x 0 0 = 0 at hq1
    linarith

/-! ### The rational trajectory and its Krylov basis -/

private abbrev MillerMat := Matrix (Fin 3) (Fin 3) ℝ

private def millerChallengeR (n : ℕ) : MillerMat :=
  fun i j => (challengeMatrix (n : ℤ) i j : ℝ)

private def millerTrajectory (n : ℕ) : MillerMat :=
  let q : ℝ := n
  !![
    4 * (q + 2) * (17 * q ^ 3 + 111 * q ^ 2 + 240 * q + 171) /
      ((q + 1) * (q + 3) * (2 * q + 3) * (2 * q + 5)),
    (q + 2) * (24 * q ^ 2 + 101 * q + 102) /
      ((q + 1) * (2 * q + 3)),
    (q + 2) * (2 * q + 5) * (16 * q ^ 2 + 81 * q + 90) /
      (2 * (q + 1) * (2 * q + 3));
    (96 * q ^ 4 + 780 * q ^ 3 + 2384 * q ^ 2 + 3273 * q + 1723) /
      ((q + 1) * (q + 2) * (q + 3) * (2 * q + 3) * (2 * q + 5)),
    (68 * q ^ 3 + 398 * q ^ 2 + 778 * q + 523) /
      (2 * (q + 1) * (q + 2) * (2 * q + 3)),
    (96 * q ^ 4 + 884 * q ^ 3 + 2970 * q ^ 2 + 4360 * q + 2403) /
      (4 * (q + 1) * (q + 2) * (2 * q + 3));
    -5 * (24 * q ^ 2 + 117 * q + 143) /
      ((q + 1) * (q + 2) * (q + 3) * (2 * q + 3) * (2 * q + 5)),
    -5 * (16 * q + 41) /
      (2 * (q + 1) * (q + 2) * (2 * q + 3)),
    (8 * q ^ 3 - 44 * q ^ 2 - 478 * q - 801) /
      (4 * (q + 1) * (q + 2) * (2 * q + 3))
  ]

private def millerSigma (n : ℕ) : ℝ :=
  -2 * ((n : ℝ) + 2) ^ 2 * ((n : ℝ) + 3) ^ 2 *
    (2 * (n : ℝ) + 5) * (2 * (n : ℝ) + 7) ^ 2

private def millerCompanion (n : ℕ) : MillerMat :=
  !![0, 0, millerC0 n;
     1, 0, millerC1 n;
     0, 1, millerC2 n]

private def millerPullMatrix (n : ℕ) : MillerMat :=
  !![millerA n, millerB n, millerR n;
     1, 0, 0;
     0, 1, 0]

private def millerE0 : MillerVec := ![1, 0, 0]
private def millerE2 : MillerVec := ![0, 0, 1]

private def millerPVec (n : ℕ) : MillerVec :=
  Matrix.mulVec (millerTrajectory n) millerE0

private def millerQVec (n : ℕ) : MillerVec :=
  Matrix.mulVec (millerTrajectory n) (millerPVec (n + 1))

private def millerRVec (n : ℕ) : MillerVec :=
  Matrix.mulVec (millerTrajectory n) (millerQVec (n + 1))

private def millerKrylov (n : ℕ) : MillerMat :=
  fun i j =>
    ![millerE0, millerPVec n, millerQVec n] j i

private theorem millerPullMatrix_mul_companion_transpose (n : ℕ) :
    millerPullMatrix n * (millerCompanion n).transpose =
      (1 : MillerMat) := by
  have hc0 : millerC0 n ≠ 0 := by
    have := (miller_coeff_bounds n).1
    linarith
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [millerPullMatrix, millerCompanion, millerA, millerB, millerR,
      Matrix.mul_apply, Fin.sum_univ_three, Matrix.one_apply] <;>
    field_simp [hc0] <;> ring

private theorem millerTrajectory_adjoint (n : ℕ) :
    (millerTrajectory n).transpose * millerChallengeR n =
      millerSigma n • (1 : MillerMat) := by
  have h1 : (n : ℝ) + 1 ≠ 0 := by positivity
  have h2 : (n : ℝ) + 2 ≠ 0 := by positivity
  have h3 : (n : ℝ) + 3 ≠ 0 := by positivity
  have h23 : 2 * (n : ℝ) + 3 ≠ 0 := by positivity
  have h25 : 2 * (n : ℝ) + 5 ≠ 0 := by positivity
  ext i j
  fin_cases i <;> fin_cases j
  all_goals
    simp only [Matrix.mul_apply, Matrix.transpose_apply,
      Matrix.smul_apply, Fin.sum_univ_three, Matrix.one_apply]
    norm_num [millerTrajectory, millerChallengeR, millerSigma,
      challengeMatrix, m11, m12, m13, m21, m22, m23, m31, m32, m33,
      Matrix.cons_val_two]
    push_cast
    field_simp [h1, h2, h3, h23, h25]
    ring

set_option maxRecDepth 100000 in
private theorem millerRVec_reconstruct (n : ℕ) (i : Fin 3) :
    millerRVec n i =
      millerC0 n * millerE0 i +
        millerC1 n * millerPVec n i +
        millerC2 n * millerQVec n i := by
  have h1 : (n : ℝ) + 1 ≠ 0 := by positivity
  have h2 : (n : ℝ) + 2 ≠ 0 := by positivity
  have h3 : (n : ℝ) + 3 ≠ 0 := by positivity
  have h4 : (n : ℝ) + 4 ≠ 0 := by positivity
  have h5 : (n : ℝ) + 5 ≠ 0 := by positivity
  have h23 : 2 * (n : ℝ) + 3 ≠ 0 := by positivity
  have h25 : 2 * (n : ℝ) + 5 ≠ 0 := by positivity
  have h27 : 2 * (n : ℝ) + 7 ≠ 0 := by positivity
  have h29 : 2 * (n : ℝ) + 9 ≠ 0 := by positivity
  have hd0 : millerD n ≠ 0 := (millerD_pos n).ne'
  have hd1 : millerD (n + 1) ≠ 0 := (millerD_pos (n + 1)).ne'
  fin_cases i
  all_goals
    simp [millerRVec, millerQVec, millerPVec, millerE0,
      millerTrajectory, millerC0, millerC1, millerC2,
      millerD, millerP, millerQ, Matrix.mulVec, dotProduct,
      Fin.sum_univ_three, Matrix.cons_val_zero, Matrix.cons_val_one,
      Matrix.cons_val_two, Nat.cast_add,
      Nat.cast_one]
    field_simp [h1, h2, h3, h4, h5, h23, h25, h27, h29, hd0, hd1]
    ring

private theorem millerTrajectory_mul_krylov_succ (n : ℕ) :
    millerTrajectory n * millerKrylov (n + 1) =
      millerKrylov n * millerCompanion n := by
  ext i j
  fin_cases j
  · simp [millerKrylov, millerCompanion, millerPVec,
      Matrix.mul_apply, Matrix.mulVec, dotProduct, Fin.sum_univ_three,
      millerE0, Matrix.cons_val_zero, Matrix.cons_val_one,
      Matrix.cons_val_two]
  · simp [millerKrylov, millerCompanion, millerQVec,
      Matrix.mul_apply, Matrix.mulVec, dotProduct, Fin.sum_univ_three,
      millerE0, Matrix.cons_val_zero, Matrix.cons_val_one,
      Matrix.cons_val_two]
  · have h := millerRVec_reconstruct n i
    simpa [millerKrylov, millerCompanion, millerRVec,
      Matrix.mul_apply, Matrix.mulVec, dotProduct, Fin.sum_univ_three,
      Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
      add_assoc, mul_comm, mul_left_comm, mul_assoc] using h

/-! ### The recessive row solution -/

private def millerDualRow (n : ℕ) : MillerVec :=
  fun j => (coordinateSign j : ℝ) * dualVector n j

private def millerRowMul (x : MillerVec) (A : MillerMat) : MillerVec :=
  fun j => ∑ i : Fin 3, x i * A i j

private theorem millerChallenge_mul_dualRow_succ (n : ℕ) (i : Fin 3) :
    ∑ j : Fin 3, millerChallengeR n i j * millerDualRow (n + 1) j =
      -dualCertLambda (n : ℝ) * millerDualRow n i := by
  have h := millerDualVector_adjoint n i
  have hsign (r c : Fin 3) :
      (positiveMatrix (n : ℤ) r c : ℝ) =
        -(coordinateSign r : ℝ) *
          millerChallengeR n r c * (coordinateSign c : ℝ) := by
    unfold millerChallengeR
    exact_mod_cast positiveMatrix_eq_sign_conjugate (n : ℤ) r c
  rw [Fin.sum_univ_three] at h ⊢
  rw [hsign i 0, hsign i 1, hsign i 2] at h
  fin_cases i <;>
    norm_num [millerDualRow, coordinateSign] at h ⊢ <;>
    linarith

private theorem millerDualRow_mul_challengeTranspose
    (n : ℕ) (i : Fin 3) :
    ∑ j : Fin 3, millerDualRow (n + 1) j * millerChallengeR n i j =
      -dualCertLambda (n : ℝ) * millerDualRow n i := by
  simpa [mul_comm] using millerChallenge_mul_dualRow_succ n i

private theorem millerChallengeTranspose_mul_trajectory
    (n : ℕ) (j k : Fin 3) :
    ∑ i : Fin 3, millerChallengeR n i j * millerTrajectory n i k =
      millerSigma n * (1 : MillerMat) j k := by
  have h := congrArg
    (fun A : MillerMat => A k j) (millerTrajectory_adjoint n)
  simp only [Matrix.mul_apply, Matrix.transpose_apply,
    Matrix.smul_apply] at h
  simpa [mul_comm, Matrix.one_apply, eq_comm] using h

private theorem millerSigma_mul_dualRow_succ (n : ℕ) (k : Fin 3) :
    millerSigma n * millerDualRow (n + 1) k =
      -dualCertLambda (n : ℝ) *
        millerRowMul (millerDualRow n) (millerTrajectory n) k := by
  have h0 := millerDualRow_mul_challengeTranspose n 0
  have h1 := millerDualRow_mul_challengeTranspose n 1
  have h2 := millerDualRow_mul_challengeTranspose n 2
  have t0 := millerChallengeTranspose_mul_trajectory n 0 k
  have t1 := millerChallengeTranspose_mul_trajectory n 1 k
  have t2 := millerChallengeTranspose_mul_trajectory n 2 k
  simp only [Fin.sum_univ_three] at h0 h1 h2 t0 t1 t2
  unfold millerRowMul
  rw [Fin.sum_univ_three]
  fin_cases k
  · simp [Matrix.one_apply] at t0 t1 t2 ⊢
    linear_combination
      (millerTrajectory n 0 0 * h0 +
        millerTrajectory n 1 0 * h1 +
        millerTrajectory n 2 0 * h2) -
      (millerDualRow (n + 1) 0 * t0 +
        millerDualRow (n + 1) 1 * t1 +
        millerDualRow (n + 1) 2 * t2)
  · simp [Matrix.one_apply] at t0 t1 t2 ⊢
    linear_combination
      (millerTrajectory n 0 1 * h0 +
        millerTrajectory n 1 1 * h1 +
        millerTrajectory n 2 1 * h2) -
      (millerDualRow (n + 1) 0 * t0 +
        millerDualRow (n + 1) 1 * t1 +
        millerDualRow (n + 1) 2 * t2)
  · simp [Matrix.one_apply] at t0 t1 t2 ⊢
    linear_combination
      (millerTrajectory n 0 2 * h0 +
        millerTrajectory n 1 2 * h1 +
        millerTrajectory n 2 2 * h2) -
      (millerDualRow (n + 1) 0 * t0 +
        millerDualRow (n + 1) 1 * t1 +
        millerDualRow (n + 1) 2 * t2)

private theorem miller_neg_sigma_div_lambda (n : ℕ) :
    -millerSigma n / dualCertLambda (n : ℝ) =
      2 * (2 * (n : ℝ) + 5) / ((n : ℝ) + 1) := by
  have h1 : (n : ℝ) + 1 ≠ 0 := by positivity
  unfold millerSigma dualCertLambda
  field_simp [h1]

private theorem millerDualRow_mul_trajectory (n : ℕ) (k : Fin 3) :
    millerRowMul (millerDualRow n) (millerTrajectory n) k =
      (2 * (2 * (n : ℝ) + 5) / ((n : ℝ) + 1)) *
        millerDualRow (n + 1) k := by
  have h := millerSigma_mul_dualRow_succ n k
  have hlam : dualCertLambda (n : ℝ) ≠ 0 := by
    unfold dualCertLambda
    positivity
  rw [← miller_neg_sigma_div_lambda n]
  field_simp [hlam]
  linarith

private def millerScaleFactor (n : ℕ) : ℝ :=
  2 * (2 * (n : ℝ) + 5) / ((n : ℝ) + 1)

private def millerScale : ℕ → ℝ
  | 0 => 1
  | n + 1 => millerScale n * millerScaleFactor n

private theorem millerScaleFactor_pos (n : ℕ) :
    0 < millerScaleFactor n := by
  unfold millerScaleFactor
  positivity

private theorem millerScaleFactor_le_ten (n : ℕ) :
    millerScaleFactor n ≤ 10 := by
  have hn : 0 < (n : ℝ) + 1 := by positivity
  rw [millerScaleFactor, div_le_iff₀ hn]
  have hnn : (0 : ℝ) ≤ n := by positivity
  nlinarith

private theorem millerScale_pos : ∀ n : ℕ, 0 < millerScale n
  | 0 => by norm_num [millerScale]
  | n + 1 => by
      rw [millerScale]
      exact mul_pos (millerScale_pos n) (millerScaleFactor_pos n)

private theorem millerScale_le_pow_ten : ∀ n : ℕ,
    millerScale n ≤ (10 : ℝ) ^ n
  | 0 => by norm_num [millerScale]
  | n + 1 => by
      calc
        millerScale (n + 1) =
            millerScale n * millerScaleFactor n := rfl
        _ ≤ (10 : ℝ) ^ n * millerScaleFactor n :=
          mul_le_mul_of_nonneg_right (millerScale_le_pow_ten n)
            (millerScaleFactor_pos n).le
        _ ≤ (10 : ℝ) ^ n * 10 :=
          mul_le_mul_of_nonneg_left (millerScaleFactor_le_ten n)
            (pow_nonneg (by norm_num) n)
        _ = (10 : ℝ) ^ (n + 1) := by rw [pow_succ]

private def millerFastRow (n : ℕ) : MillerVec :=
  fun j => millerScale n * millerDualRow n j

private theorem millerFastRow_succ (n : ℕ) :
    millerFastRow (n + 1) =
      millerRowMul (millerFastRow n) (millerTrajectory n) := by
  funext k
  calc
    millerFastRow (n + 1) k =
        millerScale n * millerScaleFactor n *
          millerDualRow (n + 1) k := by rfl
    _ = millerScale n *
          millerRowMul (millerDualRow n) (millerTrajectory n) k := by
      rw [millerDualRow_mul_trajectory]
      unfold millerScaleFactor
      ring
    _ = millerRowMul (millerFastRow n) (millerTrajectory n) k := by
      simp only [millerRowMul, millerFastRow, Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro i hi
      ring

private def millerFastScalar (n : ℕ) : ℝ :=
  millerScale n * dualVector n 0

private def millerFastState (n : ℕ) : MillerVec :=
  ![millerFastScalar n, millerFastScalar (n + 1),
    millerFastScalar (n + 2)]

private theorem millerFastRow_zero (n : ℕ) :
    millerFastRow n 0 = millerFastScalar n := by
  simp [millerFastRow, millerFastScalar, millerDualRow, coordinateSign]

private theorem millerRowMul_mul (x : MillerVec) (A B : MillerMat) :
    millerRowMul (millerRowMul x A) B =
      millerRowMul x (A * B) := by
  funext j
  fin_cases j <;>
    simp [millerRowMul, Matrix.mul_apply, Fin.sum_univ_three] <;>
    ring

private theorem millerFastState_eq_rowMul_krylov (n : ℕ) :
    millerFastState n =
      millerRowMul (millerFastRow n) (millerKrylov n) := by
  funext j
  fin_cases j
  · simp [millerFastState, millerRowMul, millerKrylov,
      millerFastRow_zero, millerFastRow, millerFastScalar, millerDualRow,
      coordinateSign, millerE0, Fin.sum_univ_three]
  · have h := congrFun (millerFastRow_succ n) 0
    simpa [millerFastState, millerRowMul, millerKrylov,
      millerFastRow_zero, millerPVec, Matrix.mulVec, dotProduct,
      Fin.sum_univ_three, millerE0, mul_comm, mul_left_comm, mul_assoc]
      using h
  · have h := congrFun (millerFastRow_succ (n + 1)) 0
    rw [millerFastRow_succ n] at h
    simp [millerFastState, millerRowMul, millerKrylov,
      millerFastRow_zero, millerQVec, millerPVec, Matrix.mulVec,
      Matrix.mul_apply, dotProduct, Fin.sum_univ_three, millerE0,
      Nat.add_assoc,
      mul_comm, mul_left_comm, mul_assoc] at h ⊢
    ring_nf at h ⊢
    exact h

private theorem millerFastState_mul_companion (n : ℕ) :
    millerFastState (n + 1) =
      millerRowMul (millerFastState n) (millerCompanion n) := by
  calc
    millerFastState (n + 1) =
        millerRowMul (millerFastRow (n + 1))
          (millerKrylov (n + 1)) :=
      millerFastState_eq_rowMul_krylov (n + 1)
    _ = millerRowMul
          (millerRowMul (millerFastRow n) (millerTrajectory n))
          (millerKrylov (n + 1)) := by rw [millerFastRow_succ]
    _ = millerRowMul (millerFastRow n)
          (millerTrajectory n * millerKrylov (n + 1)) :=
      millerRowMul_mul _ _ _
    _ = millerRowMul (millerFastRow n)
          (millerKrylov n * millerCompanion n) := by
      rw [millerTrajectory_mul_krylov_succ]
    _ = millerRowMul
          (millerRowMul (millerFastRow n) (millerKrylov n))
          (millerCompanion n) :=
      (millerRowMul_mul _ _ _).symm
    _ = millerRowMul (millerFastState n) (millerCompanion n) := by
      rw [← millerFastState_eq_rowMul_krylov]

private theorem millerFastState_succ (n : ℕ) :
    millerFastState (n + 1) = millerStep n (millerFastState n) := by
  rw [millerFastState_mul_companion]
  funext i
  fin_cases i <;>
    simp [millerRowMul, millerCompanion, millerStep,
      Fin.sum_univ_three] <;>
    ring

private def millerFastRho : ℝ := 5 / 32

private theorem millerFastScalar_eq_moment (n : ℕ) :
    millerFastScalar n =
      millerScale n * dualMoment n 0 0 0 0 := by
  simp [millerFastScalar, dualVector]

private theorem millerFastScalar_pos (n : ℕ) :
    0 < millerFastScalar n := by
  exact mul_pos (millerScale_pos n) (dualVector_zero_pos n)

private theorem millerFastScalar_abs_le (n : ℕ) :
    |millerFastScalar n| ≤ 512 * millerFastRho ^ n := by
  rw [abs_of_pos (millerFastScalar_pos n), millerFastScalar_eq_moment]
  calc
    millerScale n * dualMoment n 0 0 0 0 ≤
        millerScale n * (512 / (64 : ℝ) ^ n) :=
      mul_le_mul_of_nonneg_left (dualMoment_zero_le_fast n)
        (millerScale_pos n).le
    _ ≤ (10 : ℝ) ^ n * (512 / (64 : ℝ) ^ n) :=
      mul_le_mul_of_nonneg_right (millerScale_le_pow_ten n) (by positivity)
    _ = 512 * ((10 : ℝ) / 64) ^ n := by
      rw [div_pow]
      ring
    _ = 512 * millerFastRho ^ n := by
      norm_num [millerFastRho]

private theorem millerFastRho_nonneg : 0 ≤ millerFastRho := by
  norm_num [millerFastRho]

private theorem millerFastRho_le_one : millerFastRho ≤ 1 := by
  norm_num [millerFastRho]

private theorem millerFastRho_pow_succ_le (n : ℕ) :
    millerFastRho ^ (n + 1) ≤ millerFastRho ^ n := by
  rw [pow_succ]
  have hp : 0 ≤ millerFastRho ^ n :=
    pow_nonneg millerFastRho_nonneg n
  nlinarith [millerFastRho_le_one]

private theorem millerFastState_coord_abs_le (n : ℕ) (i : Fin 3) :
    |millerFastState n i| ≤ 512 * millerFastRho ^ n := by
  fin_cases i
  · simpa [millerFastState] using millerFastScalar_abs_le n
  · calc
      |millerFastState n 1| = |millerFastScalar (n + 1)| := by
        simp [millerFastState]
      _ ≤ 512 * millerFastRho ^ (n + 1) :=
        millerFastScalar_abs_le (n + 1)
      _ ≤ 512 * millerFastRho ^ n :=
        mul_le_mul_of_nonneg_left (millerFastRho_pow_succ_le n)
          (by norm_num)
  · have h21 :
        millerFastRho ^ (n + 2) ≤ millerFastRho ^ (n + 1) := by
      simpa [Nat.add_assoc] using millerFastRho_pow_succ_le (n + 1)
    have h20 : millerFastRho ^ (n + 2) ≤ millerFastRho ^ n :=
      h21.trans (millerFastRho_pow_succ_le n)
    calc
      |millerFastState n 2| = |millerFastScalar (n + 2)| := by
        simp [millerFastState]
      _ ≤ 512 * millerFastRho ^ (n + 2) :=
        millerFastScalar_abs_le (n + 2)
      _ ≤ 512 * millerFastRho ^ n :=
        mul_le_mul_of_nonneg_left h20 (by norm_num)

private theorem millerDefectNorm_le_of_coord
    (z : MillerChart) (x : MillerVec) (K : ℝ)
    (hz : millerInBox z) (hK : 0 ≤ K)
    (hx0 : |x 0| ≤ K) (hx1 : |x 1| ≤ K) (hx2 : |x 2| ≤ K) :
    millerDefectNorm z x ≤ (17 / 16 : ℝ) * K := by
  rcases hz with ⟨hz0, hz1, hz2, hz3⟩
  have hzx0 : z 0 * |x 0| ≤ (1 / 16 : ℝ) * K := by
    calc
      z 0 * |x 0| ≤ (1 / 16 : ℝ) * |x 0| :=
        mul_le_mul_of_nonneg_right hz1 (abs_nonneg _)
      _ ≤ (1 / 16 : ℝ) * K :=
        mul_le_mul_of_nonneg_left hx0 (by norm_num)
  have hzx1 : z 1 * |x 0| ≤ (1 / 256 : ℝ) * K := by
    calc
      z 1 * |x 0| ≤ (1 / 256 : ℝ) * |x 0| :=
        mul_le_mul_of_nonneg_right hz3 (abs_nonneg _)
      _ ≤ (1 / 256 : ℝ) * K :=
        mul_le_mul_of_nonneg_left hx0 (by norm_num)
  have hq0 :
      |millerDefect z x 0| ≤ (17 / 16 : ℝ) * K := by
    calc
      |millerDefect z x 0| = |x 1 - z 0 * x 0| := rfl
      _ ≤ |x 1| + |z 0 * x 0| := abs_sub _ _
      _ = |x 1| + z 0 * |x 0| := by
        rw [abs_mul, abs_of_nonneg hz0]
      _ ≤ K + (1 / 16 : ℝ) * K := add_le_add hx1 hzx0
      _ = (17 / 16 : ℝ) * K := by ring
  have hq1 :
      |millerDefect z x 1| ≤ (17 / 16 : ℝ) * K := by
    calc
      |millerDefect z x 1| = |x 2 - z 1 * x 0| := rfl
      _ ≤ |x 2| + |z 1 * x 0| := abs_sub _ _
      _ = |x 2| + z 1 * |x 0| := by
        rw [abs_mul, abs_of_nonneg hz2]
      _ ≤ K + (1 / 256 : ℝ) * K := add_le_add hx2 hzx1
      _ ≤ (17 / 16 : ℝ) * K := by nlinarith
  exact max_le hq0 hq1

private theorem millerFastState_defect_le (n : ℕ) :
    millerDefectNorm (millerGraph n) (millerFastState n) ≤
      (17 / 16 : ℝ) * (512 * millerFastRho ^ n) := by
  exact millerDefectNorm_le_of_coord
    (millerGraph n) (millerFastState n) (512 * millerFastRho ^ n)
    (millerGraph_mem_box n)
    (mul_nonneg (by norm_num) (pow_nonneg millerFastRho_nonneg n))
    (millerFastState_coord_abs_le n 0)
    (millerFastState_coord_abs_le n 1)
    (millerFastState_coord_abs_le n 2)

private theorem millerFastState_scaled_defect_le (n : ℕ) :
    (2 : ℝ) ^ n *
        millerDefectNorm (millerGraph n) (millerFastState n) ≤
      544 * ((5 : ℝ) / 16) ^ n := by
  calc
    (2 : ℝ) ^ n *
          millerDefectNorm (millerGraph n) (millerFastState n) ≤
        (2 : ℝ) ^ n *
          ((17 / 16 : ℝ) * (512 * millerFastRho ^ n)) :=
      mul_le_mul_of_nonneg_left (millerFastState_defect_le n)
        (pow_nonneg (by norm_num) n)
    _ = 544 * ((5 : ℝ) / 16) ^ n := by
      rw [show (5 : ℝ) / 16 = 2 * millerFastRho by
        norm_num [millerFastRho], mul_pow]
      ring

private theorem millerFastState_scaled_defect_tendsto_zero :
    Tendsto
      (fun n : ℕ =>
        (2 : ℝ) ^ n *
          millerDefectNorm (millerGraph n) (millerFastState n))
      atTop (nhds 0) := by
  have hgeom :
      Tendsto (fun n : ℕ => 544 * ((5 : ℝ) / 16) ^ n)
        atTop (nhds 0) := by
    simpa using
      (tendsto_pow_atTop_nhds_zero_of_lt_one
        (by norm_num : (0 : ℝ) ≤ 5 / 16)
        (by norm_num : (5 : ℝ) / 16 < 1)).const_mul 544
  exact tendsto_of_tendsto_of_tendsto_of_le_of_le
    tendsto_const_nhds hgeom
    (fun n => mul_nonneg (pow_nonneg (by norm_num) n)
      (millerDefectNorm_nonneg _ _))
    millerFastState_scaled_defect_le

private theorem millerFastState_zero_eq_graph :
    millerFastState 0 =
      millerFastState 0 0 • millerLift (millerGraph 0) :=
  miller_fast_eq_graph millerFastState millerFastState_succ
    millerFastState_scaled_defect_tendsto_zero

/-! ### Fixed covectors at the initial Krylov basis -/

private def millerEndDot (x y : MillerVec) : ℝ :=
  ∑ i : Fin 3, x i * y i

private def millerInitialNum : MillerVec :=
  ![30921, -32972, 8240]

private def millerInitialDen : MillerVec :=
  ![33750, -36000, 9000]

private def millerEllNum : MillerVec :=
  ![30102216645, -19896711516, 525137760]

private def millerEllDen : MillerVec :=
  ![32863650750, -21712357800, 573048000]

private theorem millerEndDot_smul_right
    (c : ℝ) (x y : MillerVec) :
    millerEndDot x (c • y) = c * millerEndDot x y := by
  simp [millerEndDot, Fin.sum_univ_three]
  ring

private theorem millerEndDot_transpose_mulVec
    (A : MillerMat) (x y : MillerVec) :
    millerEndDot x (Matrix.mulVec A.transpose y) =
      millerEndDot (Matrix.mulVec A x) y := by
  simp [millerEndDot, Matrix.mulVec, dotProduct, Fin.sum_univ_three]
  ring

private theorem millerKrylov_zero :
    millerKrylov 0 =
      !![1, 152 / 5, 195477 / 175;
         0, 1723 / 90, 1963751 / 2800;
         0, -143 / 18, -165201 / 560] := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    norm_num [millerKrylov, millerPVec, millerQVec, millerE0,
      millerTrajectory, Matrix.mulVec, dotProduct, Fin.sum_univ_three,
      Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two]

private theorem millerKrylov_zero_mul_ellNum :
    Matrix.mulVec (millerKrylov 0) millerEllNum =
      (382493 : ℝ) • millerInitialNum := by
  rw [millerKrylov_zero]
  funext i
  fin_cases i <;>
    norm_num [millerEllNum, millerInitialNum, Matrix.mulVec, dotProduct,
      Fin.sum_univ_three, Matrix.cons_val_zero, Matrix.cons_val_one,
      Matrix.cons_val_two]

private theorem millerKrylov_zero_mul_ellDen :
    Matrix.mulVec (millerKrylov 0) millerEllDen =
      (382493 : ℝ) • millerInitialDen := by
  rw [millerKrylov_zero]
  funext i
  fin_cases i <;>
    norm_num [millerEllDen, millerInitialDen, Matrix.mulVec, dotProduct,
      Fin.sum_univ_three, Matrix.cons_val_zero, Matrix.cons_val_one,
      Matrix.cons_val_two]

private theorem millerEllNum_covector (y : MillerVec) :
    millerEndDot millerEllNum
        (Matrix.mulVec (millerKrylov 0).transpose y) =
      382493 * millerEndDot millerInitialNum y := by
  rw [millerEndDot_transpose_mulVec, millerKrylov_zero_mul_ellNum]
  simp [millerEndDot, Fin.sum_univ_three]
  ring

private theorem millerEllDen_covector (y : MillerVec) :
    millerEndDot millerEllDen
        (Matrix.mulVec (millerKrylov 0).transpose y) =
      382493 * millerEndDot millerInitialDen y := by
  rw [millerEndDot_transpose_mulVec, millerKrylov_zero_mul_ellDen]
  simp [millerEndDot, Fin.sum_univ_three]
  ring

private theorem millerInitialNum_dualRow :
    millerEndDot millerInitialNum (millerDualRow 0) =
      75 * catalanConstant := by
  simpa [millerEndDot, millerInitialNum, millerDualRow,
    positiveNumerator, numerator, approximants, initialMatrix, matrixProduct,
    coordinateSign, Matrix.mul_apply, Fin.sum_univ_three,
    Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two]
    using dualVector_initial_numerator_pair

private theorem millerInitialDen_dualRow :
    millerEndDot millerInitialDen (millerDualRow 0) = 75 := by
  simpa [millerEndDot, millerInitialDen, millerDualRow,
    positiveDenominator, denominator, approximants, initialMatrix,
    matrixProduct, coordinateSign, Matrix.mul_apply, Fin.sum_univ_three,
    Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two]
    using dualVector_initial_denominator_pair

private theorem millerFastState_zero_as_krylov :
    millerFastState 0 =
      Matrix.mulVec (millerKrylov 0).transpose (millerDualRow 0) := by
  rw [millerFastState_eq_rowMul_krylov]
  funext i
  simp [millerRowMul, millerFastRow, millerScale, Matrix.mulVec, dotProduct,
    Fin.sum_univ_three]
  ring

private theorem millerEllNum_fastState_zero :
    millerEndDot millerEllNum (millerFastState 0) =
      (382493 : ℝ) * 75 * catalanConstant := by
  rw [millerFastState_zero_as_krylov, millerEllNum_covector,
    millerInitialNum_dualRow]
  ring

private theorem millerEllDen_fastState_zero :
    millerEndDot millerEllDen (millerFastState 0) =
      (382493 : ℝ) * 75 := by
  rw [millerFastState_zero_as_krylov, millerEllDen_covector,
    millerInitialDen_dualRow]

private def millerFixedCovectorRatio (z : MillerChart) : ℝ :=
  millerEndDot millerEllNum (millerLift z) /
    millerEndDot millerEllDen (millerLift z)

private theorem millerFixedCovectorDen_pos
    {z : MillerChart} (hz : millerInBox z) :
    0 < millerEndDot millerEllDen (millerLift z) := by
  rcases hz with ⟨hz0, hz1, hz2, hz3⟩
  simp [millerEndDot, millerEllDen, millerLift, Fin.sum_univ_three]
  nlinarith

private theorem millerFixedCovectorRatio_graph_eq_catalan :
    millerFixedCovectorRatio (millerGraph 0) = catalanConstant := by
  have hnum := millerEllNum_fastState_zero
  have hden := millerEllDen_fastState_zero
  rw [millerFastState_zero_eq_graph, millerEndDot_smul_right] at hnum hden
  have hk : millerFastState 0 0 ≠ 0 := by
    simpa [millerFastState] using (millerFastScalar_pos 0).ne'
  have htarget :
      millerEndDot millerEllNum (millerLift (millerGraph 0)) =
        catalanConstant *
          millerEndDot millerEllDen (millerLift (millerGraph 0)) := by
    have hprod :
        millerFastState 0 0 *
          (millerEndDot millerEllNum (millerLift (millerGraph 0)) -
            catalanConstant *
              millerEndDot millerEllDen (millerLift (millerGraph 0))) =
            0 := by
      calc
        _ = millerFastState 0 0 *
                millerEndDot millerEllNum (millerLift (millerGraph 0)) -
              catalanConstant *
                (millerFastState 0 0 *
                  millerEndDot millerEllDen (millerLift (millerGraph 0))) := by
              ring
        _ = (382493 : ℝ) * 75 * catalanConstant -
              catalanConstant * ((382493 : ℝ) * 75) := by
              rw [hnum, hden]
        _ = 0 := by ring
    exact sub_eq_zero.mp ((mul_eq_zero.mp hprod).resolve_left hk)
  unfold millerFixedCovectorRatio
  exact (div_eq_iff
    (millerFixedCovectorDen_pos (millerGraph_mem_box 0)).ne').2 htarget

private theorem millerEndDot_lift_tendsto
    (ell : MillerVec) {zseq : ℕ → MillerChart} {z : MillerChart}
    (hz : Tendsto zseq atTop (nhds z)) :
    Tendsto
      (fun h => millerEndDot ell (millerLift (zseq h)))
      atTop (nhds (millerEndDot ell (millerLift z))) := by
  have h0 := tendsto_pi_nhds.mp hz 0
  have h1 := tendsto_pi_nhds.mp hz 1
  have hc :
      Tendsto (fun _ : ℕ => ell 0) atTop (nhds (ell 0)) :=
    tendsto_const_nhds
  have ht0 :
      Tendsto (fun h : ℕ => ell 1 * zseq h 0)
        atTop (nhds (ell 1 * z 0)) :=
    tendsto_const_nhds.mul h0
  have ht1 :
      Tendsto (fun h : ℕ => ell 2 * zseq h 1)
        atTop (nhds (ell 2 * z 1)) :=
    tendsto_const_nhds.mul h1
  simpa [millerEndDot, millerLift, Fin.sum_univ_three] using
    (hc.add ht0).add ht1

private theorem millerFixedCovectorRatio_tendsto
    {zseq : ℕ → MillerChart} {z : MillerChart}
    (hz : Tendsto zseq atTop (nhds z)) (hbox : millerInBox z) :
    Tendsto (fun h => millerFixedCovectorRatio (zseq h))
      atTop (nhds (millerFixedCovectorRatio z)) := by
  have hnum := millerEndDot_lift_tendsto millerEllNum hz
  have hden := millerEndDot_lift_tendsto millerEllDen hz
  simpa [millerFixedCovectorRatio] using
    hnum.div hden (millerFixedCovectorDen_pos hbox).ne'

private theorem millerFixedCovectorRatio_pullN_moving_tendsto_catalan
    (terminal : ℕ → MillerChart)
    (hterminal : ∀ h, millerInBox (terminal h)) :
    Tendsto
      (fun h =>
        millerFixedCovectorRatio (millerPullN 0 h (terminal h)))
      atTop (nhds catalanConstant) := by
  have hpull :=
    millerPullN_moving_tendsto 0 terminal hterminal
  have hratio :=
    millerFixedCovectorRatio_tendsto hpull (millerGraph_mem_box 0)
  simpa only [millerFixedCovectorRatio_graph_eq_catalan] using hratio

/-! ### Exact matrix telescope -/

private def millerMProd : ℕ → MillerMat
  | 0 => 1
  | n + 1 => millerMProd n * millerChallengeR n

private def millerAProd : ℕ → MillerMat
  | 0 => 1
  | n + 1 => millerAProd n * millerPullMatrix n

private def millerSigmaProd : ℕ → ℝ
  | 0 => 1
  | n + 1 => millerSigmaProd n * millerSigma n

private theorem millerSigma_ne (n : ℕ) : millerSigma n ≠ 0 := by
  unfold millerSigma
  positivity

private theorem millerSigmaProd_ne : ∀ n : ℕ, millerSigmaProd n ≠ 0
  | 0 => by norm_num [millerSigmaProd]
  | n + 1 => by
      rw [millerSigmaProd]
      exact mul_ne_zero (millerSigmaProd_ne n) (millerSigma_ne n)

private theorem millerCompanionT_mul_krylovT (n : ℕ) :
    (millerCompanion n).transpose * (millerKrylov n).transpose =
      (millerKrylov (n + 1)).transpose *
        (millerTrajectory n).transpose := by
  have h := congrArg Matrix.transpose
    (millerTrajectory_mul_krylov_succ n)
  simpa using h.symm

private theorem millerKrylovT_mul_challenge (n : ℕ) :
    (millerKrylov n).transpose * millerChallengeR n =
      millerSigma n •
        (millerPullMatrix n * (millerKrylov (n + 1)).transpose) := by
  calc
    (millerKrylov n).transpose * millerChallengeR n =
        (millerPullMatrix n * (millerCompanion n).transpose) *
          ((millerKrylov n).transpose * millerChallengeR n) := by
      rw [millerPullMatrix_mul_companion_transpose]
      simp
    _ = millerPullMatrix n *
          (((millerCompanion n).transpose *
              (millerKrylov n).transpose) * millerChallengeR n) := by
      simp only [Matrix.mul_assoc]
    _ = millerPullMatrix n *
          (((millerKrylov (n + 1)).transpose *
              (millerTrajectory n).transpose) * millerChallengeR n) := by
      rw [millerCompanionT_mul_krylovT]
    _ = millerPullMatrix n *
          ((millerKrylov (n + 1)).transpose *
            ((millerTrajectory n).transpose * millerChallengeR n)) := by
      simp only [Matrix.mul_assoc]
    _ = millerPullMatrix n *
          ((millerKrylov (n + 1)).transpose *
            (millerSigma n • (1 : MillerMat))) := by
      rw [millerTrajectory_adjoint]
    _ = millerSigma n •
          (millerPullMatrix n * (millerKrylov (n + 1)).transpose) := by
      ext i j
      simp [Matrix.mul_apply, Fin.sum_univ_three]
      ring

private theorem miller_matrix_telescope : ∀ N : ℕ,
    (millerKrylov 0).transpose * millerMProd N =
      millerSigmaProd N •
        (millerAProd N * (millerKrylov N).transpose)
  | 0 => by
      simp [millerMProd, millerAProd, millerSigmaProd]
  | n + 1 => by
      calc
        (millerKrylov 0).transpose * millerMProd (n + 1) =
            ((millerKrylov 0).transpose * millerMProd n) *
              millerChallengeR n := by
          rw [millerMProd, ← Matrix.mul_assoc]
        _ = (millerSigmaProd n •
              (millerAProd n * (millerKrylov n).transpose)) *
                millerChallengeR n := by
          rw [miller_matrix_telescope n]
        _ = millerSigmaProd n •
              ((millerAProd n * (millerKrylov n).transpose) *
                millerChallengeR n) := by
          rw [Matrix.smul_mul]
        _ = millerSigmaProd n •
              (millerAProd n *
                ((millerKrylov n).transpose * millerChallengeR n)) := by
          rw [Matrix.mul_assoc]
        _ = millerSigmaProd n •
              (millerAProd n *
                (millerSigma n •
                  (millerPullMatrix n *
                    (millerKrylov (n + 1)).transpose))) := by
          rw [millerKrylovT_mul_challenge]
        _ = millerSigmaProd (n + 1) •
              (millerAProd (n + 1) *
                (millerKrylov (n + 1)).transpose) := by
          ext i j
          simp [millerSigmaProd, millerAProd, Matrix.mul_apply,
            Fin.sum_univ_three]
          ring

private def millerMatrixProductR (N : ℕ) : MillerMat :=
  fun i j => (matrixProduct N i j : ℝ)

private theorem millerMatrixProductR_succ (n : ℕ) :
    millerMatrixProductR (n + 1) =
      millerMatrixProductR n * millerChallengeR n := by
  ext i j
  simp [millerMatrixProductR, matrixProduct_succ, millerChallengeR,
    Matrix.mul_apply]

private theorem millerMatrixProductR_eq_millerMProd : ∀ N : ℕ,
    millerMatrixProductR N = millerMProd N
  | 0 => by
      ext i j
      fin_cases i <;> fin_cases j <;>
        norm_num [millerMatrixProductR, millerMProd, matrixProduct,
          Matrix.one_apply]
  | n + 1 => by
      rw [millerMatrixProductR_succ, millerMProd,
        millerMatrixProductR_eq_millerMProd n]

private theorem miller_matrix_telescope_cast (N : ℕ) :
    (millerKrylov 0).transpose * millerMatrixProductR N =
      millerSigmaProd N •
        (millerAProd N * (millerKrylov N).transpose) := by
  rw [millerMatrixProductR_eq_millerMProd]
  exact miller_matrix_telescope N

/-! ### Projective action of the pull matrices -/

private theorem millerDen_pos (n : ℕ) {z : MillerChart}
    (hz : millerInBox z) :
    0 < millerDen n z :=
  lt_of_lt_of_le (by norm_num) (miller_den_bounds n hz).1

private theorem millerPullMatrix_mulVec_lift
    (n : ℕ) {z : MillerChart} (hz : millerInBox z) :
    Matrix.mulVec (millerPullMatrix n) (millerLift z) =
      millerDen n z • millerLift (millerPull n z) := by
  have hden := (millerDen_pos n hz).ne'
  funext i
  fin_cases i
  · simp [millerPullMatrix, millerLift, millerDen, Matrix.mulVec,
      dotProduct, Fin.sum_univ_three]
  · simp [millerPullMatrix, millerLift, millerPull, Matrix.mulVec,
      dotProduct, Fin.sum_univ_three]
    field_simp [hden]
  · simp [millerPullMatrix, millerLift, millerPull, Matrix.mulVec,
      dotProduct, Fin.sum_univ_three]
    field_simp [hden]

private def millerAIntervalProd : ℕ → ℕ → MillerMat
  | _, 0 => 1
  | k, h + 1 => millerPullMatrix k * millerAIntervalProd (k + 1) h

private def millerPullScale : ℕ → ℕ → MillerChart → ℝ
  | _, 0, _ => 1
  | k, h + 1, z =>
      millerDen k (millerPullN (k + 1) h z) *
        millerPullScale (k + 1) h z

private theorem millerPullScale_pos
    (k h : ℕ) {z : MillerChart} (hz : millerInBox z) :
    0 < millerPullScale k h z := by
  induction h generalizing k with
  | zero =>
      norm_num [millerPullScale]
  | succ h ih =>
      rw [millerPullScale]
      exact mul_pos
        (millerDen_pos k (millerPullN_mem_box (k + 1) h hz))
        (ih (k := k + 1))

private theorem millerAIntervalProd_mulVec_lift
    (k h : ℕ) {z : MillerChart} (hz : millerInBox z) :
    Matrix.mulVec (millerAIntervalProd k h) (millerLift z) =
      millerPullScale k h z • millerLift (millerPullN k h z) := by
  induction h generalizing k with
  | zero =>
      simp [millerAIntervalProd, millerPullScale, millerPullN]
  | succ h ih =>
      rw [millerAIntervalProd]
      rw [← Matrix.mulVec_mulVec]
      rw [ih (k := k + 1)]
      rw [Matrix.mulVec_smul]
      rw [millerPullMatrix_mulVec_lift k
        (millerPullN_mem_box (k + 1) h hz)]
      simp [millerPullScale, millerPullN, smul_smul, mul_comm]

private theorem millerAIntervalProd_succ_right
    (k h : ℕ) :
    millerAIntervalProd k (h + 1) =
      millerAIntervalProd k h * millerPullMatrix (k + h) := by
  induction h generalizing k with
  | zero =>
      simp [millerAIntervalProd]
  | succ h ih =>
      calc
        millerAIntervalProd k (h + 2) =
            millerPullMatrix k *
              millerAIntervalProd (k + 1) (h + 1) := rfl
        _ = millerPullMatrix k *
              (millerAIntervalProd (k + 1) h *
                millerPullMatrix ((k + 1) + h)) := by
              rw [ih (k := k + 1)]
        _ = (millerPullMatrix k * millerAIntervalProd (k + 1) h) *
              millerPullMatrix (k + (h + 1)) := by
              rw [show (k + 1) + h = k + (h + 1) by omega]
              rw [← Matrix.mul_assoc]
        _ = millerAIntervalProd k (h + 1) *
              millerPullMatrix (k + (h + 1)) := rfl

private theorem millerAProd_eq_interval_zero : ∀ h : ℕ,
    millerAProd h = millerAIntervalProd 0 h
  | 0 => by simp [millerAProd, millerAIntervalProd]
  | h + 1 => by
      rw [millerAProd, millerAIntervalProd_succ_right,
        millerAProd_eq_interval_zero h]
      simp

private theorem millerAProd_mulVec_lift
    (h : ℕ) {z : MillerChart} (hz : millerInBox z) :
    Matrix.mulVec (millerAProd h) (millerLift z) =
      millerPullScale 0 h z • millerLift (millerPullN 0 h z) := by
  rw [millerAProd_eq_interval_zero]
  exact millerAIntervalProd_mulVec_lift 0 h hz

private theorem millerAProd_add_two (h : ℕ) :
    millerAProd (h + 2) =
      (millerAProd h * millerPullMatrix h) *
        millerPullMatrix (h + 1) := by
  rw [show h + 2 = (h + 1) + 1 by omega, millerAProd]
  rw [millerAProd]

/-! ### The terminal third Krylov column -/

private def millerTerminalR2 (N : ℕ) : ℝ :=
  6528 * (N : ℝ) ^ 8 + 54480 * (N : ℝ) ^ 7 +
    187144 * (N : ℝ) ^ 6 + 344616 * (N : ℝ) ^ 5 +
    371338 * (N : ℝ) ^ 4 + 239883 * (N : ℝ) ^ 3 +
    91396 * (N : ℝ) ^ 2 + 19200 * (N : ℝ) + 1800

private def millerTerminalNum2 (N : ℕ) : ℝ :=
  (N : ℝ) ^ 2 * (2 * (N : ℝ) + 1) * (2 * (N : ℝ) + 3) *
    (48 * (N : ℝ) ^ 4 + 306 * (N : ℝ) ^ 3 +
      715 * (N : ℝ) ^ 2 + 729 * (N : ℝ) + 275)

private def millerTerminalW1 (N : ℕ) : ℝ :=
  -5 * (24 * (N : ℝ) ^ 2 + 117 * (N : ℝ) + 143) /
    (((N : ℝ) + 1) * ((N : ℝ) + 2) * ((N : ℝ) + 3) *
      (2 * (N : ℝ) + 3) * (2 * (N : ℝ) + 5))

private def millerTerminalW2 (N : ℕ) : ℝ :=
  -5 * (13056 * (N : ℝ) ^ 6 + 220704 * (N : ℝ) ^ 5 +
      1545344 * (N : ℝ) ^ 4 + 5735344 * (N : ℝ) ^ 3 +
      11897116 * (N : ℝ) ^ 2 + 13075529 * (N : ℝ) + 5947236) /
    (4 * ((N : ℝ) + 1) * ((N : ℝ) + 2) ^ 2 *
      ((N : ℝ) + 3) * ((N : ℝ) + 4) *
      (2 * (N : ℝ) + 3) * (2 * (N : ℝ) + 5) ^ 2 *
      (2 * (N : ℝ) + 7))

private def millerTerminalHead (N : ℕ) : ℝ :=
  5 * (48 * (N : ℝ) ^ 4 + 306 * (N : ℝ) ^ 3 +
      715 * (N : ℝ) ^ 2 + 729 * (N : ℝ) + 275) /
    (((N : ℝ) + 1) ^ 2 * ((N : ℝ) + 2) ^ 2 *
      (2 * (N : ℝ) + 3) * (2 * (N : ℝ) + 5) ^ 2)

private def millerTerminalScalar (N : ℕ) : ℝ :=
  5 * millerTerminalR2 N /
    ((N : ℝ) ^ 2 * ((N : ℝ) + 1) ^ 2 *
      ((N : ℝ) + 2) ^ 2 * (2 * (N : ℝ) + 1) *
      (2 * (N : ℝ) + 3) ^ 2 * (2 * (N : ℝ) + 5) ^ 2)

private def millerTerminalChart (N : ℕ) : MillerChart :=
  ![millerTerminalNum2 N / millerTerminalR2 N, 0]

private def millerTerminalVec0 (N : ℕ) : MillerVec :=
  ![0, millerTerminalW1 N, millerTerminalW2 N]

private def millerTerminalVec1 (N : ℕ) : MillerVec :=
  ![millerTerminalHead N, 0, millerTerminalW1 N]

private theorem millerTerminalR2_pos (N : ℕ) :
    0 < millerTerminalR2 N := by
  unfold millerTerminalR2
  positivity

private theorem millerTerminalNum2_nonneg (N : ℕ) :
    0 ≤ millerTerminalNum2 N := by
  unfold millerTerminalNum2
  positivity

private theorem millerTerminal_sub_nonneg (N : ℕ) :
    0 ≤ millerTerminalR2 N - 16 * millerTerminalNum2 N := by
  unfold millerTerminalR2 millerTerminalNum2
  ring_nf
  positivity

private theorem millerTerminalChart_mem_box (N : ℕ) :
    millerInBox (millerTerminalChart N) := by
  have hR := millerTerminalR2_pos N
  have hnum := millerTerminalNum2_nonneg N
  have hupper :
      millerTerminalNum2 N / millerTerminalR2 N ≤ (1 : ℝ) / 16 := by
    rw [div_le_iff₀ hR]
    nlinarith [millerTerminal_sub_nonneg N]
  exact ⟨div_nonneg hnum hR.le, hupper,
    by norm_num [millerTerminalChart],
    by norm_num [millerTerminalChart]⟩

private theorem millerTerminalScalar_pos (h : ℕ) :
    0 < millerTerminalScalar (h + 2) := by
  unfold millerTerminalScalar
  have hR := millerTerminalR2_pos (h + 2)
  positivity

set_option maxRecDepth 100000 in
private theorem millerTerminal_e2_formula (N : ℕ) :
    Matrix.mulVec (millerKrylov N).transpose millerE2 =
      millerTerminalVec0 N := by
  have h1 : (N : ℝ) + 1 ≠ 0 := by positivity
  have h2 : (N : ℝ) + 2 ≠ 0 := by positivity
  have h3 : (N : ℝ) + 3 ≠ 0 := by positivity
  have h4 : (N : ℝ) + 4 ≠ 0 := by positivity
  have h23 : 2 * (N : ℝ) + 3 ≠ 0 := by positivity
  have h25 : 2 * (N : ℝ) + 5 ≠ 0 := by positivity
  have h27 : 2 * (N : ℝ) + 7 ≠ 0 := by positivity
  funext i
  fin_cases i
  · simp [millerKrylov, millerTerminalVec0, millerE0, millerE2,
      Matrix.mulVec, dotProduct, Fin.sum_univ_three]
  · simp [millerKrylov, millerTerminalVec0, millerTerminalW1,
      millerPVec, millerE0, millerE2, millerTrajectory, Matrix.mulVec,
      dotProduct, Fin.sum_univ_three, Matrix.cons_val_zero,
      Matrix.cons_val_one, Matrix.cons_val_two]
  · simp [millerKrylov, millerTerminalVec0, millerTerminalW2,
      millerQVec, millerPVec, millerE0, millerE2, millerTrajectory,
      Matrix.mulVec, dotProduct, Fin.sum_univ_three,
      Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
      Nat.cast_add, Nat.cast_one]
    field_simp [h1, h2, h3, h4, h23, h25, h27]
    ring

set_option maxRecDepth 100000 in
private theorem millerTerminal_first_pull (h : ℕ) :
    Matrix.mulVec (millerPullMatrix (h + 1))
        (millerTerminalVec0 (h + 2)) =
      millerTerminalVec1 (h + 2) := by
  funext i
  fin_cases i
  · norm_num [millerPullMatrix, millerA, millerB, millerR,
      millerC0, millerC1, millerC2, millerP, millerQ, millerD,
      millerTerminalVec0, millerTerminalVec1, millerTerminalW1,
      millerTerminalW2, millerTerminalHead, Matrix.mulVec, dotProduct,
      Fin.sum_univ_three, Matrix.cons_val_two, Nat.cast_add, Nat.cast_one]
    field_simp
    ring_nf
  · simp [millerPullMatrix, millerTerminalVec0, millerTerminalVec1,
      Matrix.mulVec, dotProduct, Fin.sum_univ_three]
  · simp [millerPullMatrix, millerTerminalVec0, millerTerminalVec1,
      Matrix.mulVec, dotProduct, Fin.sum_univ_three]

set_option maxRecDepth 100000 in
private theorem millerTerminal_second_pull (h : ℕ) :
    Matrix.mulVec (millerPullMatrix h)
        (millerTerminalVec1 (h + 2)) =
      ![millerTerminalScalar (h + 2),
        millerTerminalHead (h + 2), 0] := by
  funext i
  fin_cases i
  · norm_num [millerPullMatrix, millerA, millerB, millerR,
      millerC0, millerC1, millerC2, millerP, millerQ, millerD,
      millerTerminalVec1, millerTerminalW1, millerTerminalHead,
      millerTerminalScalar, millerTerminalR2, Matrix.mulVec, dotProduct,
      Fin.sum_univ_three, Matrix.cons_val_two, Nat.cast_add, Nat.cast_one]
    field_simp
    ring_nf
  · simp [millerPullMatrix, millerTerminalVec1, Matrix.mulVec, dotProduct,
      Fin.sum_univ_three]
  · simp [millerPullMatrix, millerTerminalVec1, Matrix.mulVec, dotProduct,
      Fin.sum_univ_three]

set_option maxRecDepth 100000 in
private theorem millerTerminal_vector_lift (h : ℕ) :
    ![millerTerminalScalar (h + 2),
        millerTerminalHead (h + 2), 0] =
      millerTerminalScalar (h + 2) •
        millerLift (millerTerminalChart (h + 2)) := by
  funext i
  fin_cases i
  · simp [millerLift]
  · simp [millerLift, millerTerminalChart, millerTerminalHead,
      millerTerminalScalar, millerTerminalNum2, millerTerminalR2,
      Nat.cast_add, Nat.cast_one]
    field_simp
  · simp [millerLift, millerTerminalChart]

private theorem miller_terminal_e2 (h : ℕ) :
    Matrix.mulVec (millerPullMatrix h)
        (Matrix.mulVec (millerPullMatrix (h + 1))
          (Matrix.mulVec (millerKrylov (h + 2)).transpose millerE2)) =
      millerTerminalScalar (h + 2) •
        millerLift (millerTerminalChart (h + 2)) := by
  rw [millerTerminal_e2_formula, millerTerminal_first_pull,
    millerTerminal_second_pull, millerTerminal_vector_lift]

/-! ### The literal third column and its projective normalization -/

private def millerChallengeCol2 (N : ℕ) : MillerVec :=
  Matrix.mulVec (millerMatrixProductR N) millerE2

private def millerCompanionCol2 (N : ℕ) : MillerVec :=
  Matrix.mulVec (millerKrylov 0).transpose (millerChallengeCol2 N)

private def millerFinalChart (h : ℕ) : MillerChart :=
  millerPullN 0 h (millerTerminalChart (h + 2))

private def millerCommonScale (h : ℕ) : ℝ :=
  millerSigmaProd (h + 2) * millerTerminalScalar (h + 2) *
    millerPullScale 0 h (millerTerminalChart (h + 2))

private theorem millerChallengeCol2_apply
    (N : ℕ) (i : Fin 3) :
    millerChallengeCol2 N i = (matrixProduct N i 2 : ℝ) := by
  simp [millerChallengeCol2, millerMatrixProductR, millerE2,
    Matrix.mulVec, dotProduct, Fin.sum_univ_three]

private theorem millerInitialNum_challengeCol2 (N : ℕ) :
    millerEndDot millerInitialNum (millerChallengeCol2 N) =
      (numerator N 2 : ℝ) := by
  simp [millerEndDot, millerInitialNum, millerChallengeCol2_apply,
    numerator, approximants, initialMatrix, Matrix.mul_apply,
    Fin.sum_univ_three] <;>
    push_cast <;> ring

private theorem millerInitialDen_challengeCol2 (N : ℕ) :
    millerEndDot millerInitialDen (millerChallengeCol2 N) =
      (denominator N 2 : ℝ) := by
  simp [millerEndDot, millerInitialDen, millerChallengeCol2_apply,
    denominator, approximants, initialMatrix, Matrix.mul_apply,
    Fin.sum_univ_three] <;>
    push_cast <;> ring

private theorem millerCompanionCol2_normalized (h : ℕ) :
    millerCompanionCol2 (h + 2) =
      millerCommonScale h • millerLift (millerFinalChart h) := by
  have htelescope :=
    congrArg (fun X : MillerMat => Matrix.mulVec X millerE2)
      (miller_matrix_telescope_cast (h + 2))
  have htelescopeVec :
      millerCompanionCol2 (h + 2) =
        millerSigmaProd (h + 2) •
          Matrix.mulVec (millerAProd (h + 2))
            (Matrix.mulVec (millerKrylov (h + 2)).transpose millerE2) := by
    calc
      millerCompanionCol2 (h + 2) =
          Matrix.mulVec
            ((millerKrylov 0).transpose *
              millerMatrixProductR (h + 2)) millerE2 := by
            simp [millerCompanionCol2, millerChallengeCol2,
              Matrix.mulVec_mulVec]
      _ = Matrix.mulVec
            (millerSigmaProd (h + 2) •
              (millerAProd (h + 2) *
                (millerKrylov (h + 2)).transpose)) millerE2 := htelescope
      _ = millerSigmaProd (h + 2) •
            Matrix.mulVec
              (millerAProd (h + 2) *
                (millerKrylov (h + 2)).transpose) millerE2 := by
            rw [Matrix.smul_mulVec]
      _ = millerSigmaProd (h + 2) •
            Matrix.mulVec (millerAProd (h + 2))
              (Matrix.mulVec
                (millerKrylov (h + 2)).transpose millerE2) := by
            rw [Matrix.mulVec_mulVec]
  rw [htelescopeVec, millerAProd_add_two]
  rw [← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec]
  rw [miller_terminal_e2 h]
  rw [Matrix.mulVec_smul]
  rw [millerAProd_mulVec_lift h
    (millerTerminalChart_mem_box (h + 2))]
  simp [millerCommonScale, millerFinalChart, smul_smul]
  ring

private theorem millerNumerator_common_pairing (h : ℕ) :
    millerCommonScale h *
        millerEndDot millerEllNum
          (millerLift (millerFinalChart h)) =
      (382493 : ℝ) * (numerator (h + 2) 2 : ℝ) := by
  have hpair := millerEllNum_covector
    (millerChallengeCol2 (h + 2))
  change
    millerEndDot millerEllNum
        (millerCompanionCol2 (h + 2)) =
      382493 *
        millerEndDot millerInitialNum
          (millerChallengeCol2 (h + 2)) at hpair
  rw [millerCompanionCol2_normalized,
    millerEndDot_smul_right,
    millerInitialNum_challengeCol2] at hpair
  exact hpair

private theorem millerDenominator_common_pairing (h : ℕ) :
    millerCommonScale h *
        millerEndDot millerEllDen
          (millerLift (millerFinalChart h)) =
      (382493 : ℝ) * (denominator (h + 2) 2 : ℝ) := by
  have hpair := millerEllDen_covector
    (millerChallengeCol2 (h + 2))
  change
    millerEndDot millerEllDen
        (millerCompanionCol2 (h + 2)) =
      382493 *
        millerEndDot millerInitialDen
          (millerChallengeCol2 (h + 2)) at hpair
  rw [millerCompanionCol2_normalized,
    millerEndDot_smul_right,
    millerInitialDen_challengeCol2] at hpair
  exact hpair

private theorem challengeRatio_two_eq_fixedCovectorRatio (h : ℕ) :
    challengeRatio (h + 2) 2 =
      millerFixedCovectorRatio
        (millerPullN 0 h (millerTerminalChart (h + 2))) := by
  have hphysical :
      (denominator (h + 2) 2 : ℝ) ≠ 0 := by
    exact_mod_cast denominator_ne_zero (h + 2) 2
  have hterminalBox :
      millerInBox (millerTerminalChart (h + 2)) :=
    millerTerminalChart_mem_box (h + 2)
  have hfinalBox :
      millerInBox (millerFinalChart h) := by
    exact millerPullN_mem_box 0 h hterminalBox
  have hfixed :
      millerEndDot millerEllDen
        (millerLift (millerFinalChart h)) ≠ 0 :=
    (millerFixedCovectorDen_pos hfinalBox).ne'
  have hnum := millerNumerator_common_pairing h
  have hden := millerDenominator_common_pairing h
  unfold challengeRatio millerFixedCovectorRatio
  change
    (numerator (h + 2) 2 : ℝ) /
        (denominator (h + 2) 2 : ℝ) =
      millerEndDot millerEllNum (millerLift (millerFinalChart h)) /
        millerEndDot millerEllDen (millerLift (millerFinalChart h))
  apply (div_eq_div_iff hphysical hfixed).2
  have hscaled :
      (382493 : ℝ) *
          ((numerator (h + 2) 2 : ℝ) *
              millerEndDot millerEllDen
                (millerLift (millerFinalChart h)) -
            millerEndDot millerEllNum
                (millerLift (millerFinalChart h)) *
              (denominator (h + 2) 2 : ℝ)) = 0 := by
    calc
      _ = ((382493 : ℝ) * (numerator (h + 2) 2 : ℝ)) *
              millerEndDot millerEllDen
                (millerLift (millerFinalChart h)) -
            millerEndDot millerEllNum
                (millerLift (millerFinalChart h)) *
              ((382493 : ℝ) *
                (denominator (h + 2) 2 : ℝ)) := by
            ring
      _ = (millerCommonScale h *
              millerEndDot millerEllNum
                (millerLift (millerFinalChart h))) *
              millerEndDot millerEllDen
                (millerLift (millerFinalChart h)) -
            millerEndDot millerEllNum
                (millerLift (millerFinalChart h)) *
              (millerCommonScale h *
                millerEndDot millerEllDen
                  (millerLift (millerFinalChart h))) := by
            rw [← hnum, ← hden]
      _ = 0 := by ring
  have hzero :
      (numerator (h + 2) 2 : ℝ) *
          millerEndDot millerEllDen
            (millerLift (millerFinalChart h)) -
        millerEndDot millerEllNum
            (millerLift (millerFinalChart h)) *
          (denominator (h + 2) 2 : ℝ) = 0 :=
    (mul_eq_zero.mp hscaled).resolve_left (by norm_num)
  exact sub_eq_zero.mp hzero

private theorem challengeRatio_two_shift_tendsto_catalan :
    Filter.Tendsto (fun h : ℕ => challengeRatio (h + 2) 2)
      Filter.atTop (nhds catalanConstant) := by
  have hratio :=
    millerFixedCovectorRatio_pullN_moving_tendsto_catalan
      (fun h : ℕ => millerTerminalChart (h + 2))
      (fun h => millerTerminalChart_mem_box (h + 2))
  simpa only [challengeRatio_two_eq_fixedCovectorRatio] using hratio

private theorem challengeRatio_two_tendsto_catalan :
    Filter.Tendsto (fun N : ℕ => challengeRatio N 2)
      Filter.atTop (nhds catalanConstant) := by
  exact (tendsto_add_atTop_iff_nat 2).mp
    challengeRatio_two_shift_tendsto_catalan

private theorem commonLimit_eq_catalan_miller :
    commonLimit = catalanConstant :=
  tendsto_nhds_unique
    (challengeRatio_tendsto_common 2)
    challengeRatio_two_tendsto_catalan

theorem catalanError_over_denominator_tendsto_zero_miller (j : Fin 3) :
    Filter.Tendsto (fun N => catalanError N j / (denominator N j : ℝ))
      Filter.atTop (nhds 0) := by
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
  simpa [commonLimit_eq_catalan_miller] using
    (Filter.Tendsto.sub
      (tendsto_const_nhds :
        Filter.Tendsto (fun _ : ℕ => catalanConstant)
          Filter.atTop (nhds catalanConstant))
      hratio)

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
  have hpairSum := positiveCatalanError_pairing_zero_of_adjoint
    dualVector (fun n => dualCertLambda (n : ℝ))
      millerDualVector_adjoint
      dualVector_initial_error_pair N
  have hpair :
      positiveCatalanError N 0 * dualVector N 0 +
          positiveCatalanError N 1 * dualVector N 1 +
        positiveCatalanError N 2 * dualVector N 2 = 0 := by
    simpa only [Fin.sum_univ_three] using hpairSum
  rcases ratioMinor_projective_cone N with
    ⟨hm0, hxLower, hxUpper, hyLower, hyUpper⟩
  have hm1 : 0 < ratioMinor N 1 := by
    by_contra h
    have hm1le : ratioMinor N 1 ≤ 0 := le_of_not_gt h
    have hnum : ((N : ℝ) + 1) * ratioMinor N 1 ≤ 0 :=
      mul_nonpos_of_nonneg_of_nonpos (by positivity) hm1le
    have hquot :
        ((N : ℝ) + 1) * ratioMinor N 1 / ratioMinor N 0 ≤ 0 :=
      div_nonpos_of_nonpos_of_nonneg hnum hm0.le
    unfold wedgeX at hxLower
    linarith
  have hm2 : 0 ≤ ratioMinor N 2 := by
    by_contra h
    have hm2neg : ratioMinor N 2 < 0 := lt_of_not_ge h
    have hnum : ((N : ℝ) + 1) ^ 2 * ratioMinor N 2 < 0 :=
      mul_neg_of_pos_of_neg (sq_pos_of_pos (by positivity)) hm2neg
    have hquot :
        ((N : ℝ) + 1) ^ 2 * ratioMinor N 2 / ratioMinor N 0 < 0 :=
      div_neg_of_neg_of_pos hnum hm0
    unfold wedgeY at hyLower
    linarith
  have h01 :
      0 < positiveCatalanError N 1 * (positiveDenominator N 0 : ℝ) -
        positiveCatalanError N 0 * (positiveDenominator N 1 : ℝ) := by
    calc
      _ = ratioMinor N 0 := by
        simp [positiveCatalanError_eq, ratioMinor]
        ring
      _ > 0 := hm0
  have h02 :
      0 < positiveCatalanError N 2 * (positiveDenominator N 0 : ℝ) -
        positiveCatalanError N 0 * (positiveDenominator N 2 : ℝ) := by
    calc
      _ = ratioMinor N 1 := by
        simp [positiveCatalanError_eq, ratioMinor]
        ring
      _ > 0 := hm1
  have h12 :
      0 ≤ positiveCatalanError N 2 * (positiveDenominator N 1 : ℝ) -
        positiveCatalanError N 1 * (positiveDenominator N 2 : ℝ) := by
    calc
      _ = ratioMinor N 2 := by
        simp [positiveCatalanError_eq, ratioMinor, Matrix.cons_val_two]
        ring
      _ ≥ 0 := hm2
  have hq0 : 0 < (positiveDenominator N 0 : ℝ) := by
    exact_mod_cast positiveDenominator_pos N 0
  have hq1 : 0 < (positiveDenominator N 1 : ℝ) := by
    exact_mod_cast positiveDenominator_pos N 1
  have hq2 : 0 < (positiveDenominator N 2 : ℝ) := by
    exact_mod_cast positiveDenominator_pos N 2
  exact endpoint_bracket_of_positive_pair
    (positiveCatalanError N) (fun j => (positiveDenominator N j : ℝ))
      (dualVector N) hq0 hq1 hq2
      (dualVector_zero_pos N) (dualVector_one_pos N) (dualVector_two_pos N)
      hpair h01 h02 h12

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

theorem catalanError_over_denominator_tendsto_zero_exterior (j : Fin 3) :
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

theorem catalanError_over_denominator_tendsto_zero (j : Fin 3) :
    Filter.Tendsto (fun N => catalanError N j / (denominator N j : ℝ))
      Filter.atTop (nhds 0) :=
  catalanError_over_denominator_tendsto_zero_miller j

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
