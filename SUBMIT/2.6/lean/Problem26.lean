/-
  Ramanujan Challenge Problem 2.6: A Series for ζ(2) + ζ(3)

  The three-term recurrence has Poincaré roots 1 and 1/4. Its recessive
  solution is hypergeometric, with exact ratio v_{n+1}/v_n = (n+3)²/(2(n+4)(2n+7)).
  The Casorati determinant and variation of constants yield a closed form
  for the dominant solution u_n, and the connection formula gives
  2077/720 + Σ u_j = ζ(2) + ζ(3).

  Reference: Xiang Huang, "Solution to Ramanujan Challenge Problem 2.6", July 2026.
-/
import Mathlib.NumberTheory.ZetaValues
import Mathlib.Data.Nat.Factorial.Basic
import Mathlib.Data.Nat.Choose.Basic
import Mathlib.Data.Nat.Choose.Central
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Analysis.SpecialFunctions.Gamma.Beta
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic
import Mathlib.MeasureTheory.Integral.DominatedConvergence
import Mathlib.Analysis.SpecialFunctions.Log.Deriv
import RamanujanChallenge.Dilogarithm

noncomputable section

namespace RamanujanChallenge.P26

open Set TopologicalSpace

/-! ## The recurrence coefficients

A(n) u_n = B(n) u_{n-1} - C(n) u_{n-2}  for n ≥ 3, where:
  A(n) = 2(n+3)³(2n+5)(3n+5)
  B(n) = (n+2)²(15n³ + 85n² + 155n + 93)
  C(n) = (n+1)³(n+2)(3n+8)
-/

def coeff_A (n : ℤ) : ℤ := 2 * (n + 3) ^ 3 * (2 * n + 5) * (3 * n + 5)

def coeff_B (n : ℤ) : ℤ := (n + 2) ^ 2 * (15 * n ^ 3 + 85 * n ^ 2 + 155 * n + 93)

def coeff_C (n : ℤ) : ℤ := (n + 1) ^ 3 * (n + 2) * (3 * n + 8)

/-- The first coefficient of the challenge recurrence, over `ℚ`. -/
def coeffAQ (n : ℕ) : ℚ :=
  2 * ((n : ℚ) + 3) ^ 3 * (2 * (n : ℚ) + 5) * (3 * (n : ℚ) + 5)

/-- The middle coefficient of the challenge recurrence, over `ℚ`. -/
def coeffBQ (n : ℕ) : ℚ :=
  ((n : ℚ) + 2) ^ 2 *
    (15 * (n : ℚ) ^ 3 + 85 * (n : ℚ) ^ 2 + 155 * (n : ℚ) + 93)

/-- The last coefficient of the challenge recurrence, over `ℚ`. -/
def coeffCQ (n : ℕ) : ℚ :=
  ((n : ℚ) + 1) ^ 3 * ((n : ℚ) + 2) * (3 * (n : ℚ) + 8)

/-- A sequence satisfies the recurrence printed in Problem 2.6.

The natural-number parameter is shifted by three, so no truncated subtraction
appears in the statement. -/
def SatisfiesRecurrence26 (u : ℕ → ℚ) : Prop :=
  ∀ n : ℕ,
    coeffAQ (n + 3) * u (n + 3) =
      coeffBQ (n + 3) * u (n + 2) - coeffCQ (n + 3) * u (n + 1)

/-! ## Poincaré analysis

The leading coefficients give:
  12n⁵ r² - 15n⁵ r + 3n⁵ = 0
  ⟹ 4r² - 5r + 1 = 0
  ⟹ (4r - 1)(r - 1) = 0
  ⟹ r = 1 or r = 1/4

Root r = 1 gives the dominant solution: u_n ~ c₁ n⁻⁴
Root r = 1/4 gives the recessive solution: u_n^{(1/4)} ~ c₂ · 4⁻ⁿ n⁻³/²
-/

theorem poincare_poly_roots :
    ∀ r : ℚ, 4 * r ^ 2 - 5 * r + 1 = 0 ↔ r = 1 ∨ r = 1 / 4 := by
  intro r
  constructor
  · intro h
    have : (4 * r - 1) * (r - 1) = 0 := by ring_nf; linarith
    rcases mul_eq_zero.mp this with h1 | h2
    · right; linarith
    · left; linarith
  · rintro (rfl | rfl) <;> ring

/-! ## Hypergeometric recessive solution

The recurrence has a hypergeometric solution v_n with ratio:
  v_{n+1} / v_n = (n+3)² / (2(n+4)(2n+7)) → 1/4

With normalization:
  v_n = 1120 · ((n+2)!)² / (2n+6)!
-/

def recessiveRatio (n : ℕ) : ℚ :=
  (↑n + 3) ^ 2 / (2 * (↑n + 4) * (2 * ↑n + 7))

theorem recessiveRatio_limit :
    Filter.Tendsto recessiveRatio Filter.atTop (nhds (1 / 4 : ℚ)) := by
  have h₁ :
      Filter.Tendsto (fun n : ℕ ↦ ((3 : ℚ) + 1 * n) / (8 + 2 * n))
        Filter.atTop (nhds ((1 : ℚ) / 2)) :=
    tendsto_add_mul_div_add_mul_atTop_nhds (𝕜 := ℚ) 3 8 1
      (by norm_num : (2 : ℚ) ≠ 0)
  have h₂ :
      Filter.Tendsto (fun n : ℕ ↦ ((3 : ℚ) + 1 * n) / (7 + 2 * n))
        Filter.atTop (nhds ((1 : ℚ) / 2)) :=
    tendsto_add_mul_div_add_mul_atTop_nhds (𝕜 := ℚ) 3 7 1
      (by norm_num : (2 : ℚ) ≠ 0)
  convert h₁.mul h₂ using 1
  · funext n
    simp only [recessiveRatio, one_mul]
    field_simp (discharger := positivity)
    ring
  · norm_num

/-! ## Initial values

u₁ = -93/4480, u₂ = -117/14000
-/

/-! ## Exact first-order reduction

The second-order operator has a right first-order factor.  For the particular
initial values in the challenge, the factor is not zero but the following
explicit rational forcing:

`2(2n+7)(n+4)u_{n+1} - (n+3)²u_n
    = -3(3n+11)/((n+3)(n+4)²)`.

The compatibility identity for the forcing is what makes every solution of
this first-order equation satisfy the printed second-order recurrence.
-/

/-- The right first-order factor of the recurrence operator. -/
def firstOrderResidual26 (u : ℕ → ℚ) (n : ℕ) : ℚ :=
  2 * (2 * (n : ℚ) + 7) * ((n : ℚ) + 4) * u (n + 1) -
    ((n : ℚ) + 3) ^ 2 * u n

/-- The rational forcing selected by the challenge's initial values. -/
def firstOrderForcing26 (n : ℕ) : ℚ :=
  -3 * (3 * (n : ℚ) + 11) / (((n : ℚ) + 3) * ((n : ℚ) + 4) ^ 2)

/-- The normalized recessive hypergeometric solution.  The normalization makes
`recessiveTerm26 1 = 1`. -/
def recessiveTerm26 (n : ℕ) : ℚ :=
  1120 * ((n + 2).factorial : ℚ) ^ 2 / ((2 * n + 6).factorial : ℚ)

/-- The hypergeometric term is killed by the right first-order factor. -/
theorem recessiveTerm26_firstOrder (n : ℕ) :
    2 * (2 * (n : ℚ) + 7) * ((n : ℚ) + 4) * recessiveTerm26 (n + 1) =
      ((n : ℚ) + 3) ^ 2 * recessiveTerm26 n := by
  have hnum :
      (((n + 3).factorial : ℕ) : ℚ) =
        (n + 3 : ℚ) * ((n + 2).factorial : ℚ) := by
    norm_cast
  have hden :
      (((2 * (n + 1) + 6).factorial : ℕ) : ℚ) =
        (2 * n + 8 : ℚ) * (2 * n + 7 : ℚ) *
          ((2 * n + 6).factorial : ℚ) := by
    norm_cast
    rw [show 2 * (n + 1) + 6 = (2 * n + 7) + 1 by omega,
      Nat.factorial_succ,
      show 2 * n + 7 = (2 * n + 6) + 1 by omega,
      Nat.factorial_succ]
    ring
  simp only [recessiveTerm26, hnum, hden]
  field_simp
  ring

theorem recessiveTerm26_one : recessiveTerm26 1 = 1 := by
  norm_num [recessiveTerm26]

/-- The central-binomial summand in the reduction-of-order formula. -/
def reductionSummand26 (j : ℕ) : ℚ :=
  (3 * (j : ℚ) + 8) / (((j : ℚ) + 2) * ((j : ℚ) + 3) ^ 2) *
    (Nat.choose (2 * j + 4) (j + 2) : ℚ)

/-- The finite coefficient multiplying the recessive solution.  For `n ≥ 1`
the range is exactly `j = 2, …, n`. -/
def reductionCoefficient26 (n : ℕ) : ℚ :=
  -93 / 4480 -
    3 / 1120 * ∑ k ∈ Finset.range (n - 1), reductionSummand26 (k + 2)

/-- Central-binomial term in the normalized inverse-binomial double sum. -/
def inverseBinomialA26 (m : ℕ) : ℚ :=
  (Nat.choose (2 * m) m : ℚ) * (3 * (m : ℚ) + 2) /
    ((m : ℚ) * ((m : ℚ) + 1) ^ 2)

/-- The inverse central-binomial factor. -/
def inverseBinomialC26 (k : ℕ) : ℚ :=
  1 / ((k : ℚ) ^ 2 * (Nat.choose (2 * k) k : ℚ))

/-- The beta moment underlying the inverse central-binomial summand. -/
theorem inverseBinomialBetaMoment26 (k : ℕ) (hk : 0 < k) :
    (∫ x : ℝ in 0..1, x ^ (k - 1) * (1 - x) ^ k) =
      1 / ((k : ℝ) * (Nat.choose (2 * k) k : ℝ)) := by
  have hbeta :=
    Complex.betaIntegral_eq_Gamma_mul_div (k : ℂ) (k + 1 : ℂ)
      (by simp; positivity) (by simp; positivity)
  have hkcast : (k : ℂ) = ((k - 1 : ℕ) : ℂ) + 1 := by
    norm_cast
    omega
  have htwok :
      (k : ℂ) + (k + 1 : ℂ) = ((2 * k : ℕ) : ℂ) + 1 := by
    push_cast
    ring
  have hGk :
      Complex.Gamma (k : ℂ) = ((k - 1).factorial : ℂ) := by
    rw [hkcast, Complex.Gamma_nat_eq_factorial (k - 1)]
  have hGk1 :
      Complex.Gamma ((k : ℂ) + 1) = (k.factorial : ℂ) := by
    exact Complex.Gamma_nat_eq_factorial k
  have hG2k1 :
      Complex.Gamma ((k : ℂ) + ((k : ℂ) + 1)) =
        ((2 * k).factorial : ℂ) := by
    rw [htwok, Complex.Gamma_nat_eq_factorial (2 * k)]
  rw [hGk, hGk1, hG2k1] at hbeta
  have hexp₁ : (k : ℂ) - 1 = ((k - 1 : ℕ) : ℂ) := by
    rw [hkcast]
    ring
  have hexp₂ : (k : ℂ) + 1 - 1 = (k : ℂ) := by ring
  rw [Complex.betaIntegral, hexp₁, hexp₂] at hbeta
  simp_rw [Complex.cpow_natCast] at hbeta
  have hcast :
      (∫ x : ℝ in 0..1,
          ((x : ℂ) ^ (k - 1) * (1 - (x : ℂ)) ^ k)) =
        ((∫ x : ℝ in 0..1, x ^ (k - 1) * (1 - x) ^ k : ℝ) : ℂ) := by
    rw [← intervalIntegral.integral_ofReal]
    congr 1
    funext x
    push_cast
    rfl
  rw [hcast] at hbeta
  have hfact :
      (Nat.choose (2 * k) k : ℕ) * k.factorial * k.factorial =
        (2 * k).factorial := by
    convert Nat.add_choose_mul_factorial_mul_factorial k k using 1 <;> ring
  have hkfact :
      k.factorial = k * (k - 1).factorial :=
    (Nat.mul_factorial_pred hk.ne').symm
  have hnat :
      (k - 1).factorial * k.factorial *
          (k * Nat.choose (2 * k) k) =
        (2 * k).factorial := by
    calc
      (k - 1).factorial * k.factorial *
            (k * Nat.choose (2 * k) k) =
          Nat.choose (2 * k) k * k.factorial * k.factorial := by
            rw [hkfact]
            ring
      _ = (2 * k).factorial := hfact
  have hreal :
      (∫ x : ℝ in 0..1, x ^ (k - 1) * (1 - x) ^ k) =
        ((k - 1).factorial : ℝ) * (k.factorial : ℝ) /
          ((2 * k).factorial : ℝ) := by
    rw [← Complex.ofReal_inj]
    push_cast
    exact hbeta
  rw [hreal]
  have hfacpos : (0 : ℝ) < ((2 * k).factorial : ℕ) := by positivity
  have hchoosepos : (0 : ℝ) < Nat.choose (2 * k) k := by
    exact_mod_cast Nat.choose_pos (by omega : k ≤ 2 * k)
  field_simp
  exact_mod_cast (by simpa [mul_assoc] using hnat)

/-- Integral representation of one positive-index inverse-binomial term. -/
theorem inverseBinomialC26_eq_integral (k : ℕ) (hk : 0 < k) :
    ((inverseBinomialC26 k : ℚ) : ℝ) =
      ∫ x : ℝ in 0..1,
        (x ^ (k - 1) * (1 - x) ^ k) / (k : ℝ) := by
  rw [intervalIntegral.integral_div, inverseBinomialBetaMoment26 k hk]
  simp only [inverseBinomialC26]
  push_cast
  have hkR : (k : ℝ) ≠ 0 := by positivity
  have hchooseR : (Nat.choose (2 * k) k : ℝ) ≠ 0 := by
    exact_mod_cast (Nat.choose_pos (by omega : k ≤ 2 * k)).ne'
  field_simp

/-- A uniformly summable continuous beta kernel, indexed from zero. -/
def inverseBinomialKernel26 (n : ℕ) : C(ℝ, ℝ) where
  toFun x := x ^ n * (1 - x) ^ (n + 1) / (n + 1 : ℝ)
  continuous_toFun := by fun_prop

theorem inverseBinomialKernel26_integral (n : ℕ) :
    (∫ x : ℝ in 0..1, inverseBinomialKernel26 n x) =
      ((inverseBinomialC26 (n + 1) : ℚ) : ℝ) := by
  simpa [inverseBinomialKernel26, Nat.cast_add, Nat.cast_one] using
    (inverseBinomialC26_eq_integral (n + 1) (by omega)).symm

theorem inverseBinomialKernel26_norm_le (n : ℕ) :
    ‖(inverseBinomialKernel26 n).restrict
        (⟨Set.uIcc (0 : ℝ) 1, isCompact_uIcc⟩ : Compacts ℝ)‖
      ≤ ((1 : ℝ) / 4) ^ n := by
  rw [ContinuousMap.norm_le _ (by positivity)]
  intro x
  rcases x.property with ⟨hx0, hx1⟩
  simp only [min_eq_left (by norm_num : (0 : ℝ) ≤ 1),
    max_eq_right (by norm_num : (0 : ℝ) ≤ 1)] at hx0 hx1
  simp only [inverseBinomialKernel26, ContinuousMap.coe_mk, Compacts.coe_mk,
    ContinuousMap.restrict_apply, Real.norm_eq_abs]
  have hx : 0 ≤ x.1 := hx0
  have h1x : 0 ≤ 1 - x.1 := sub_nonneg.mpr hx1
  have hkernel :
      0 ≤ x.1 ^ n * (1 - x.1) ^ (n + 1) / (n + 1 : ℝ) := by
    positivity
  rw [abs_of_nonneg hkernel]
  have hprodnonneg : 0 ≤ x.1 * (1 - x.1) := mul_nonneg hx h1x
  have hprod : x.1 * (1 - x.1) ≤ (1 : ℝ) / 4 := by
    nlinarith [sq_nonneg (x.1 - 1 / 2)]
  have hratio : (1 - x.1) / (n + 1 : ℝ) ≤ 1 := by
    apply (div_le_one (by positivity)).2
    have hn : (1 : ℝ) ≤ n + 1 := by
      exact_mod_cast Nat.succ_le_succ (Nat.zero_le n)
    linarith
  calc
    x.1 ^ n * (1 - x.1) ^ (n + 1) / (n + 1 : ℝ) =
        (x.1 * (1 - x.1)) ^ n *
          ((1 - x.1) / (n + 1 : ℝ)) := by
            rw [pow_succ, mul_pow]
            ring
    _ ≤ ((1 : ℝ) / 4) ^ n *
          ((1 - x.1) / (n + 1 : ℝ)) := by
            gcongr
    _ ≤ ((1 : ℝ) / 4) ^ n * 1 :=
      mul_le_mul_of_nonneg_left hratio (by positivity)
    _ = ((1 : ℝ) / 4) ^ n := mul_one _

theorem inverseBinomialKernel26_norm_summable :
    Summable (fun n : ℕ =>
      ‖(inverseBinomialKernel26 n).restrict
        (⟨Set.uIcc (0 : ℝ) 1, isCompact_uIcc⟩ : Compacts ℝ)‖) := by
  apply Summable.of_nonneg_of_le
  · intro n
    positivity
  · exact inverseBinomialKernel26_norm_le
  · exact summable_geometric_of_norm_lt_one (by norm_num)

/-- Summation and beta integration commute for the positive-index
inverse-binomial series. -/
theorem inverseBinomialC26_succ_hasSum_integral :
    HasSum (fun n : ℕ => ((inverseBinomialC26 (n + 1) : ℚ) : ℝ))
      (∫ x : ℝ in 0..1, ∑' n : ℕ, inverseBinomialKernel26 n x) := by
  have h := intervalIntegral.hasSum_intervalIntegral_of_summable_norm
    (a := (0 : ℝ)) (b := 1) inverseBinomialKernel26_norm_summable
  apply h.congr_fun
  intro n
  exact (inverseBinomialKernel26_integral n).symm

/-- Integral representation of the whole inverse-binomial series. -/
theorem inverseBinomialC26_hasSum_integral :
    HasSum (fun k : ℕ => ((inverseBinomialC26 k : ℚ) : ℝ))
      (∫ x : ℝ in 0..1, ∑' n : ℕ, inverseBinomialKernel26 n x) := by
  have h :
      HasSum (fun k : ℕ => ((inverseBinomialC26 k : ℚ) : ℝ))
        ((∫ x : ℝ in 0..1, ∑' n : ℕ, inverseBinomialKernel26 n x) +
          ∑ i ∈ Finset.range 1, ((inverseBinomialC26 i : ℚ) : ℝ)) :=
    (hasSum_nat_add_iff (f :=
      fun k : ℕ => ((inverseBinomialC26 k : ℚ) : ℝ)) 1).mp
        inverseBinomialC26_succ_hasSum_integral
  simpa [inverseBinomialC26] using h

/-- On the open right half of the unit interval, the beta kernels sum to an
elementary logarithmic function. -/
theorem inverseBinomialKernel26_hasSum_log
    (x : ℝ) (hx0 : 0 < x) (hx1 : x ≤ 1) :
    HasSum (fun n : ℕ => inverseBinomialKernel26 n x)
      (-Real.log (1 - x + x ^ 2) / x) := by
  have hprod0 : 0 ≤ x * (1 - x) :=
    mul_nonneg hx0.le (sub_nonneg.mpr hx1)
  have hprod4 : x * (1 - x) ≤ (1 : ℝ) / 4 := by
    nlinarith [sq_nonneg (x - 1 / 2)]
  have habs : |x * (1 - x)| < 1 := by
    rw [abs_of_nonneg hprod0]
    linarith
  have h :=
    (Real.hasSum_pow_div_log_of_abs_lt_one habs).div_const x
  convert h using 1
  · funext n
    simp only [inverseBinomialKernel26, ContinuousMap.coe_mk]
    have hxne : x ≠ 0 := ne_of_gt hx0
    rw [mul_pow, pow_succ]
    field_simp
    ring
  · congr 2
    ring

theorem inverseBinomialKernel26_tsum_log
    (x : ℝ) (hx0 : 0 < x) (hx1 : x ≤ 1) :
    (∑' n : ℕ, inverseBinomialKernel26 n x) =
      -Real.log (1 - x + x ^ 2) / x :=
  (inverseBinomialKernel26_hasSum_log x hx0 hx1).tsum_eq

/-- The classical inverse-binomial sum is reduced to one elementary
logarithmic integral. -/
theorem inverseBinomialC26_hasSum_logIntegral :
    HasSum (fun k : ℕ => ((inverseBinomialC26 k : ℚ) : ℝ))
      (∫ x : ℝ in 0..1, -Real.log (1 - x + x ^ 2) / x) := by
  convert inverseBinomialC26_hasSum_integral using 1
  apply intervalIntegral.integral_congr_ae
  exact Filter.Eventually.of_forall fun x hx => by
    have hx' : x ∈ Set.Ioc (0 : ℝ) 1 := by
      simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hx
    exact (inverseBinomialKernel26_tsum_log x hx'.1 hx'.2).symm

/-- The normalized sequence whose first three terms absorb `2077/720` and
whose tail is the challenge sequence. -/
def normalizedQ26 (k : ℕ) : ℚ :=
  3 * inverseBinomialC26 k *
    (2 - ∑ m ∈ Finset.range k, inverseBinomialA26 m)

/-- The nested inverse-binomial term that contains the genuinely new
weight-three evaluation. -/
def inverseBinomialDTerm26 (k : ℕ) : ℚ :=
  inverseBinomialC26 k *
    ∑ m ∈ Finset.range k, inverseBinomialA26 m

/-- First-order hypergeometric ratio for the inverse central-binomial
coefficient. -/
theorem inverseBinomialC26_succ (k : ℕ) (hk : 0 < k) :
    inverseBinomialC26 (k + 1) =
      (k : ℚ) ^ 2 / (2 * (k + 1) * (2 * k + 1)) *
        inverseBinomialC26 k := by
  have hchooseNat := Nat.succ_mul_centralBinom_succ k
  rw [Nat.centralBinom_eq_two_mul_choose,
    Nat.centralBinom_eq_two_mul_choose] at hchooseNat
  have hchoose :
      ((k + 1 : ℕ) : ℚ) *
          (Nat.choose (2 * (k + 1)) (k + 1) : ℚ) =
        2 * (2 * (k : ℚ) + 1) *
          (Nat.choose (2 * k) k : ℚ) := by
    exact_mod_cast hchooseNat
  have hkQ : (k : ℚ) ≠ 0 := by positivity
  have hk1Q : ((k : ℚ) + 1) ≠ 0 := by positivity
  have hchooseK : (Nat.choose (2 * k) k : ℚ) ≠ 0 := by
    exact_mod_cast (Nat.choose_pos (by omega : k ≤ 2 * k)).ne'
  have hchooseK1 :
      (Nat.choose (2 * (k + 1)) (k + 1) : ℚ) ≠ 0 := by
    exact_mod_cast
      (Nat.choose_pos (by omega : k + 1 ≤ 2 * (k + 1))).ne'
  unfold inverseBinomialC26
  push_cast at hchoose ⊢
  field_simp [hkQ, hk1Q, hchooseK, hchooseK1]
  nlinarith [hchoose]

/-- The nested term itself satisfies an inhomogeneous first-order recurrence
with rational forcing.  This is the local algebra behind the remaining
weight-three summation problem. -/
theorem inverseBinomialDTerm26_succ (k : ℕ) (hk : 0 < k) :
    inverseBinomialDTerm26 (k + 1) =
      (k : ℚ) ^ 2 / (2 * (k + 1) * (2 * k + 1)) *
          inverseBinomialDTerm26 k +
        (3 * (k : ℚ) + 2) /
          (2 * k * (k + 1) ^ 3 * (2 * k + 1)) := by
  rw [inverseBinomialDTerm26, inverseBinomialDTerm26,
    Finset.sum_range_succ, inverseBinomialC26_succ k hk]
  have hkQ : (k : ℚ) ≠ 0 := by positivity
  have hk1Q : ((k : ℚ) + 1) ≠ 0 := by positivity
  have hoddQ : (2 * (k : ℚ) + 1) ≠ 0 := by positivity
  have hchooseK : (Nat.choose (2 * k) k : ℚ) ≠ 0 := by
    exact_mod_cast (Nat.choose_pos (by omega : k ≤ 2 * k)).ne'
  simp only [inverseBinomialA26, inverseBinomialC26]
  field_simp [hkQ, hk1Q, hoddQ, hchooseK]

theorem normalizedQ26_decomposition (k : ℕ) :
    normalizedQ26 k =
      6 * inverseBinomialC26 k - 3 * inverseBinomialDTerm26 k := by
  unfold normalizedQ26 inverseBinomialDTerm26
  ring

theorem inverseBinomialA26_zero : inverseBinomialA26 0 = 0 := by
  norm_num [inverseBinomialA26]

theorem inverseBinomialA26_one : inverseBinomialA26 1 = 5 / 2 := by
  norm_num [inverseBinomialA26]

theorem inverseBinomialA26_two : inverseBinomialA26 2 = 8 / 3 := by
  norm_num [inverseBinomialA26, Nat.choose]

theorem inverseBinomialA26_three : inverseBinomialA26 3 = 55 / 12 := by
  norm_num [inverseBinomialA26, Nat.choose]

theorem inverseBinomialA26_initial_sum :
    (∑ m ∈ Finset.range 4, inverseBinomialA26 m) = 39 / 4 := by
  norm_num [Finset.sum_range_succ, inverseBinomialA26, Nat.choose]

theorem reductionSummand26_eq_inverseBinomialA26 (j : ℕ) :
    reductionSummand26 j = inverseBinomialA26 (j + 2) := by
  simp only [reductionSummand26, inverseBinomialA26]
  rw [show 2 * (j + 2) = 2 * j + 4 by omega]
  push_cast
  ring

/-- The factorial form of the recessive solution is exactly the inverse
central-binomial factor after shifting the index by three. -/
theorem recessiveTerm26_eq_inverseBinomialC26 (n : ℕ) :
    recessiveTerm26 n = 1120 * inverseBinomialC26 (n + 3) := by
  have hchooseNat :
      Nat.choose (2 * n + 6) (n + 3) * (n + 3).factorial *
          (n + 3).factorial =
        (2 * n + 6).factorial := by
    simpa [show (n + 3) + (n + 3) = 2 * n + 6 by omega] using
      Nat.add_choose_mul_factorial_mul_factorial (n + 3) (n + 3)
  have hchooseQ :
      (Nat.choose (2 * n + 6) (n + 3) : ℚ) *
          ((n + 3).factorial : ℚ) * ((n + 3).factorial : ℚ) =
        ((2 * n + 6).factorial : ℚ) := by
    exact_mod_cast hchooseNat
  have hfac_ne :
      ((n + 3).factorial : ℚ) * ((n + 3).factorial : ℚ) ≠ 0 := by
    positivity
  have hchooseDiv :
      (Nat.choose (2 * n + 6) (n + 3) : ℚ) =
        ((2 * n + 6).factorial : ℚ) /
          (((n + 3).factorial : ℚ) * ((n + 3).factorial : ℚ)) := by
    apply (eq_div_iff hfac_ne).2
    simpa [mul_comm, mul_left_comm, mul_assoc] using hchooseQ
  have hnum :
      (((n + 3).factorial : ℕ) : ℚ) =
        (n + 3 : ℚ) * ((n + 2).factorial : ℚ) := by
    norm_cast
  simp only [recessiveTerm26, inverseBinomialC26]
  rw [show 2 * (n + 3) = 2 * n + 6 by omega, hchooseDiv, hnum]
  push_cast
  field_simp

theorem normalizedQ26_one : normalizedQ26 1 = 3 := by
  norm_num [normalizedQ26, inverseBinomialC26, inverseBinomialA26,
    Finset.sum_range_succ]

theorem normalizedQ26_zero : normalizedQ26 0 = 0 := by
  norm_num [normalizedQ26, inverseBinomialC26]

theorem normalizedQ26_two : normalizedQ26 2 = -1 / 16 := by
  norm_num [normalizedQ26, inverseBinomialC26, inverseBinomialA26,
    Finset.sum_range_succ, Nat.choose]

theorem normalizedQ26_three : normalizedQ26 3 = -19 / 360 := by
  norm_num [normalizedQ26, inverseBinomialC26, inverseBinomialA26,
    Finset.sum_range_succ, Nat.choose]

theorem normalizedQ26_initial_sum :
    normalizedQ26 1 + normalizedQ26 2 + normalizedQ26 3 = 2077 / 720 := by
  rw [normalizedQ26_one, normalizedQ26_two, normalizedQ26_three]
  norm_num

theorem normalizedQ26_range_four_real :
    (∑ k ∈ Finset.range 4, ((normalizedQ26 k : ℚ) : ℝ)) =
      (2077 : ℝ) / 720 := by
  norm_num [Finset.sum_range_succ, normalizedQ26_zero, normalizedQ26_one,
    normalizedQ26_two, normalizedQ26_three]

theorem inverseBinomialA26_sum_split (n : ℕ) (hn : 1 ≤ n) :
    (∑ m ∈ Finset.range (n + 3), inverseBinomialA26 m) =
      39 / 4 +
        ∑ k ∈ Finset.range (n - 1), reductionSummand26 (k + 2) := by
  rw [show n + 3 = 4 + (n - 1) by omega,
    Finset.sum_range_add, inverseBinomialA26_initial_sum]
  congr 1
  apply Finset.sum_congr rfl
  intro k _
  rw [reductionSummand26_eq_inverseBinomialA26]
  congr 1
  omega

theorem reductionCoefficient26_one : reductionCoefficient26 1 = -93 / 4480 := by
  norm_num [reductionCoefficient26]

theorem reductionCoefficient26_succ (n : ℕ) (hn : 1 ≤ n) :
    reductionCoefficient26 (n + 1) - reductionCoefficient26 n =
      -3 / 1120 * reductionSummand26 (n + 1) := by
  simp only [reductionCoefficient26]
  rw [Nat.add_sub_cancel]
  have hsum :
      (∑ k ∈ Finset.range n, reductionSummand26 (k + 2)) =
        (∑ k ∈ Finset.range (n - 1), reductionSummand26 (k + 2)) +
          reductionSummand26 (n + 1) := by
    rw [show n = (n - 1) + 1 by omega, Finset.sum_range_succ]
    congr 1
  rw [hsum]
  ring

/-- The central-binomial summand converts the homogeneous solution into the
explicit rational forcing. -/
theorem reductionSummand26_forcing (n : ℕ) :
    ((n : ℚ) + 3) ^ 2 * recessiveTerm26 n *
        (-3 / 1120 * reductionSummand26 (n + 1)) =
      firstOrderForcing26 n := by
  have hchooseNat :
      Nat.choose (2 * n + 6) (n + 3) * (n + 3).factorial *
          (n + 3).factorial =
        (2 * n + 6).factorial := by
    simpa [show (n + 3) + (n + 3) = 2 * n + 6 by omega] using
      Nat.add_choose_mul_factorial_mul_factorial (n + 3) (n + 3)
  have hchooseQ :
      (Nat.choose (2 * n + 6) (n + 3) : ℚ) *
          ((n + 3).factorial : ℚ) * ((n + 3).factorial : ℚ) =
        ((2 * n + 6).factorial : ℚ) := by
    exact_mod_cast hchooseNat
  have hfac_ne :
      ((n + 3).factorial : ℚ) * ((n + 3).factorial : ℚ) ≠ 0 := by
    positivity
  have hchooseDiv :
      (Nat.choose (2 * n + 6) (n + 3) : ℚ) =
        ((2 * n + 6).factorial : ℚ) /
          (((n + 3).factorial : ℚ) * ((n + 3).factorial : ℚ)) := by
    apply (eq_div_iff hfac_ne).2
    simpa [mul_comm, mul_left_comm, mul_assoc] using hchooseQ
  have hnum :
      (((n + 3).factorial : ℕ) : ℚ) =
        (n + 3 : ℚ) * ((n + 2).factorial : ℚ) := by
    norm_cast
  simp only [recessiveTerm26, reductionSummand26, firstOrderForcing26]
  push_cast
  rw [show 2 * (n + 1) + 4 = 2 * n + 6 by omega,
    show n + 1 + 2 = n + 3 by omega]
  rw [hchooseDiv, hnum]
  field_simp
  ring

/-- Exact Ore-factorization identity, written without Ore-polynomial syntax. -/
theorem recurrence26_factorization (u : ℕ → ℚ) (n : ℕ) :
    coeffAQ (n + 2) * u (n + 2) -
        coeffBQ (n + 2) * u (n + 1) +
        coeffCQ (n + 2) * u n =
      ((n : ℚ) + 5) ^ 2 * (3 * (n : ℚ) + 11) *
          firstOrderResidual26 u (n + 1) -
        ((n : ℚ) + 3) * ((n : ℚ) + 4) * (3 * (n : ℚ) + 14) *
          firstOrderResidual26 u n := by
  simp only [coeffAQ, coeffBQ, coeffCQ, firstOrderResidual26]
  push_cast
  ring

/-- The explicit forcing is annihilated by the left first-order factor. -/
theorem firstOrderForcing26_compatibility (n : ℕ) :
    ((n : ℚ) + 5) ^ 2 * (3 * (n : ℚ) + 11) *
        firstOrderForcing26 (n + 1) -
      ((n : ℚ) + 3) * ((n : ℚ) + 4) * (3 * (n : ℚ) + 14) *
        firstOrderForcing26 n = 0 := by
  simp only [firstOrderForcing26]
  push_cast
  field_simp
  ring

/-- Any solution of the factored first-order equation satisfies the challenge
second-order recurrence. -/
theorem satisfiesRecurrence26_of_firstOrder
    (u : ℕ → ℚ)
    (hu : ∀ n : ℕ, 1 ≤ n →
      firstOrderResidual26 u n = firstOrderForcing26 n) :
    SatisfiesRecurrence26 u := by
  intro n
  have h₀ := hu (n + 1) (by omega)
  have h₁ := hu (n + 2) (by omega)
  have hfac := recurrence26_factorization u (n + 1)
  have hcompat := firstOrderForcing26_compatibility (n + 1)
  rw [h₀, h₁] at hfac
  rw [hcompat] at hfac
  linarith

/-- The challenge sequence, defined by its exact first-order reduction.

The value at index zero is outside the challenge and is set to zero. -/
def challengeU26 : ℕ → ℚ
  | 0 => 0
  | 1 => -93 / 4480
  | n + 2 =>
      (((n : ℚ) + 4) ^ 2 * challengeU26 (n + 1) +
          firstOrderForcing26 (n + 1)) /
        (2 * (2 * ((n : ℚ) + 1) + 7) * ((n : ℚ) + 5))

theorem challengeU26_zero : challengeU26 0 = 0 := rfl

theorem challengeU26_one : challengeU26 1 = -93 / 4480 := rfl

theorem challengeU26_firstOrder (n : ℕ) (hn : 1 ≤ n) :
    firstOrderResidual26 challengeU26 n = firstOrderForcing26 n := by
  cases n with
  | zero => omega
  | succ k =>
      simp only [firstOrderResidual26]
      rw [show k + 1 + 1 = k + 2 by omega, challengeU26]
      simp only [firstOrderForcing26]
      push_cast
      field_simp
      ring

theorem firstOrderForcing26_neg (n : ℕ) : firstOrderForcing26 n < 0 := by
  simp only [firstOrderForcing26]
  apply div_neg_of_neg_of_pos
  · have h : (0 : ℚ) < 3 * (3 * (n : ℚ) + 11) := by positivity
    nlinarith
  · positivity

/-- Every challenge summand is negative. -/
theorem challengeU26_neg (n : ℕ) (hn : 1 ≤ n) : challengeU26 n < 0 := by
  obtain ⟨k, rfl⟩ := Nat.exists_eq_add_of_le hn
  induction k with
  | zero => norm_num [challengeU26]
  | succ k ih =>
      have hrec := challengeU26_firstOrder (1 + k) (by omega)
      have ha :
          0 < 2 * (2 * ((1 + k : ℕ) : ℚ) + 7) *
            (((1 + k : ℕ) : ℚ) + 4) := by positivity
      have hb : 0 < (((1 + k : ℕ) : ℚ) + 3) ^ 2 := by positivity
      have hg := firstOrderForcing26_neg (1 + k)
      have hrhs :
          (((1 + k : ℕ) : ℚ) + 3) ^ 2 * challengeU26 (1 + k) +
              firstOrderForcing26 (1 + k) < 0 :=
        add_neg (mul_neg_of_pos_of_neg hb (ih (by omega))) hg
      simp only [firstOrderResidual26] at hrec
      have heq :
          2 * (2 * ((1 + k : ℕ) : ℚ) + 7) *
                (((1 + k : ℕ) : ℚ) + 4) * challengeU26 (1 + k + 1) =
            (((1 + k : ℕ) : ℚ) + 3) ^ 2 * challengeU26 (1 + k) +
              firstOrderForcing26 (1 + k) := by
        linarith
      by_contra hnonneg
      have hu : 0 ≤ challengeU26 (1 + (k + 1)) := le_of_not_gt hnonneg
      have hprod :
          0 ≤ 2 * (2 * ((1 + k : ℕ) : ℚ) + 7) *
              (((1 + k : ℕ) : ℚ) + 4) * challengeU26 (1 + k + 1) :=
        mul_nonneg ha.le (by simpa [add_assoc] using hu)
      linarith

/-- The elementary induction step for the `O(n⁻⁴)` bound. -/
theorem challengeU26_decay_step (n : ℕ) (hn : 2 ≤ n) :
    ((n : ℚ) + 3) ^ 2 * (-10 / (n : ℚ) ^ 4) +
        firstOrderForcing26 n ≥
      2 * (2 * (n : ℚ) + 7) * ((n : ℚ) + 4) *
        (-10 / ((n : ℚ) + 1) ^ 4) := by
  obtain ⟨k, rfl⟩ := Nat.exists_eq_add_of_le hn
  simp only [firstOrderForcing26]
  push_cast
  field_simp (discharger := positivity)
  ring_nf
  have hk₁ : (0 : ℚ) ≤ (k : ℚ) := by positivity
  have hk₂ : (0 : ℚ) ≤ (k : ℚ) ^ 2 := by positivity
  have hk₃ : (0 : ℚ) ≤ (k : ℚ) ^ 3 := by positivity
  have hk₄ : (0 : ℚ) ≤ (k : ℚ) ^ 4 := by positivity
  have hk₅ : (0 : ℚ) ≤ (k : ℚ) ^ 5 := by positivity
  have hk₆ : (0 : ℚ) ≤ (k : ℚ) ^ 6 := by positivity
  have hk₇ : (0 : ℚ) ≤ (k : ℚ) ^ 7 := by positivity
  have hk₈ : (0 : ℚ) ≤ (k : ℚ) ^ 8 := by positivity
  have hk₉ : (0 : ℚ) ≤ (k : ℚ) ^ 9 := by positivity
  nlinarith

theorem challengeU26_two : challengeU26 2 = -117 / 14000 := by
  have h := challengeU26_firstOrder 1 (by omega)
  norm_num [firstOrderResidual26, firstOrderForcing26, challengeU26_one] at h ⊢
  linarith

/-- A simple explicit `O(n⁻⁴)` bound, sufficient for absolute convergence. -/
theorem challengeU26_lower_bound (n : ℕ) (hn : 2 ≤ n) :
    -10 / (n : ℚ) ^ 4 ≤ challengeU26 n := by
  obtain ⟨k, rfl⟩ := Nat.exists_eq_add_of_le hn
  induction k with
  | zero => norm_num [challengeU26_two]
  | succ k ih =>
      have ih' :
          -10 / ((2 + k : ℕ) : ℚ) ^ 4 ≤ challengeU26 (2 + k) :=
        ih (by omega)
      have hrec := challengeU26_firstOrder (2 + k) (by omega)
      have hstep := challengeU26_decay_step (2 + k) (by omega)
      have ha :
          0 < 2 * (2 * ((2 + k : ℕ) : ℚ) + 7) *
            (((2 + k : ℕ) : ℚ) + 4) := by positivity
      have hb : 0 ≤ (((2 + k : ℕ) : ℚ) + 3) ^ 2 := by positivity
      have hbu :
          (((2 + k : ℕ) : ℚ) + 3) ^ 2 *
              (-10 / ((2 + k : ℕ) : ℚ) ^ 4) ≤
            (((2 + k : ℕ) : ℚ) + 3) ^ 2 * challengeU26 (2 + k) :=
        mul_le_mul_of_nonneg_left ih' hb
      simp only [firstOrderResidual26] at hrec
      have heq :
          2 * (2 * ((2 + k : ℕ) : ℚ) + 7) *
                (((2 + k : ℕ) : ℚ) + 4) * challengeU26 (2 + k + 1) =
            (((2 + k : ℕ) : ℚ) + 3) ^ 2 * challengeU26 (2 + k) +
              firstOrderForcing26 (2 + k) := by
        linarith
      have hmul :
          2 * (2 * ((2 + k : ℕ) : ℚ) + 7) *
                (((2 + k : ℕ) : ℚ) + 4) *
                (-10 / (((2 + k : ℕ) : ℚ) + 1) ^ 4) ≤
            2 * (2 * ((2 + k : ℕ) : ℚ) + 7) *
              (((2 + k : ℕ) : ℚ) + 4) * challengeU26 (2 + k + 1) := by
        rw [heq]
        exact hstep.trans (by
          simpa [add_comm] using
            (add_le_add_right hbu (firstOrderForcing26 (2 + k))))
      have hbound :
          -10 / (((2 + k : ℕ) : ℚ) + 1) ^ 4 ≤ challengeU26 (2 + k + 1) := by
        exact le_of_mul_le_mul_left (by simpa [mul_assoc] using hmul) ha
      simpa [Nat.add_assoc, add_assoc] using hbound

/-- Absolute form of the `O(n⁻⁴)` bound. -/
theorem challengeU26_abs_bound (n : ℕ) (hn : 2 ≤ n) :
    |challengeU26 n| ≤ 10 / (n : ℚ) ^ 4 := by
  rw [abs_of_neg (challengeU26_neg n (by omega))]
  have h := challengeU26_lower_bound n hn
  calc
    -challengeU26 n ≤ -(-10 / (n : ℚ) ^ 4) := neg_le_neg h
    _ = 10 / (n : ℚ) ^ 4 := by ring

/-- The challenge series is absolutely summable. -/
theorem challengeU26_summable :
    Summable (fun n : ℕ => ((challengeU26 (n + 1) : ℚ) : ℝ)) := by
  have hmajorant :
      Summable (fun n : ℕ => 10 * ((((n + 1 : ℕ) : ℝ) ^ 4)⁻¹)) := by
    apply Summable.mul_left
    rw [summable_nat_add_iff (f := fun n : ℕ => (((n : ℝ) ^ 4)⁻¹)) 1]
    exact Real.summable_nat_pow_inv.mpr (by norm_num)
  have htail :
      Summable (fun n : ℕ => ((challengeU26 (n + 2) : ℚ) : ℝ)) := by
    apply Summable.of_norm_bounded hmajorant
    intro n
    have hq := challengeU26_abs_bound (n + 2) (by omega)
    have hr :
        |((challengeU26 (n + 2) : ℚ) : ℝ)| ≤
          10 / (((n + 2 : ℕ) : ℝ) ^ 4) := by
      have hcast := Rat.cast_mono (K := ℝ) hq
      simpa using hcast
    rw [Real.norm_eq_abs]
    calc
      |((challengeU26 (n + 2) : ℚ) : ℝ)| ≤
          10 / (((n + 2 : ℕ) : ℝ) ^ 4) := hr
      _ ≤ 10 * ((((n + 1 : ℕ) : ℝ) ^ 4)⁻¹) := by
        rw [div_eq_mul_inv]
        gcongr
        norm_num
  have hshift :
      Summable
        (fun n : ℕ => ((challengeU26 ((n + 1) + 1) : ℚ) : ℝ)) := by
    simpa [Nat.add_assoc, Nat.add_comm, Nat.add_left_comm] using htail
  exact (summable_nat_add_iff
    (f := fun n : ℕ => ((challengeU26 (n + 1) : ℚ) : ℝ)) 1).mp hshift

/-- The recursively defined sequence satisfies the exact challenge recurrence. -/
theorem challengeU26_recurrence : SatisfiesRecurrence26 challengeU26 := by
  exact satisfiesRecurrence26_of_firstOrder challengeU26 challengeU26_firstOrder

/-- The leading recurrence coefficient never vanishes on natural indices. -/
theorem coeffAQ_pos (n : ℕ) : 0 < coeffAQ n := by
  simp only [coeffAQ]
  positivity

/-- The two printed initial values uniquely determine a sequence satisfying the
challenge recurrence (from index one onward). -/
theorem recurrence26_unique
    (u v : ℕ → ℚ)
    (hu₁ : u 1 = -93 / 4480) (hu₂ : u 2 = -117 / 14000)
    (hv₁ : v 1 = -93 / 4480) (hv₂ : v 2 = -117 / 14000)
    (hu : SatisfiesRecurrence26 u) (hv : SatisfiesRecurrence26 v) :
    ∀ n : ℕ, u (n + 1) = v (n + 1) := by
  intro n
  induction n using Nat.twoStepInduction with
  | zero => linarith
  | one => linarith
  | more n h₁ h₂ =>
      have hru := hu n
      have hrv := hv n
      rw [h₁, h₂] at hru
      have hmul :
          coeffAQ (n + 3) * (u (n + 3) - v (n + 3)) = 0 := by
        rw [mul_sub]
        linarith
      have hA : coeffAQ (n + 3) ≠ 0 := (coeffAQ_pos (n + 3)).ne'
      exact sub_eq_zero.mp ((mul_eq_zero.mp hmul).resolve_left hA)

/-- Any sequence with the challenge initial data and recurrence is the explicit
first-order sequence above. -/
theorem eq_challengeU26_of_spec
    (u : ℕ → ℚ)
    (hu₁ : u 1 = -93 / 4480)
    (hu₂ : u 2 = -117 / 14000)
    (hu : SatisfiesRecurrence26 u) :
    ∀ n : ℕ, u (n + 1) = challengeU26 (n + 1) :=
  recurrence26_unique u challengeU26 hu₁ hu₂
    challengeU26_one challengeU26_two hu challengeU26_recurrence

/-- The exact initial values and recurrence from the challenge. -/
theorem challengeU26_spec :
    challengeU26 1 = -93 / 4480 ∧
      challengeU26 2 = -117 / 14000 ∧
      SatisfiesRecurrence26 challengeU26 :=
  ⟨challengeU26_one, challengeU26_two, challengeU26_recurrence⟩

/-- The reduction-of-order expression, with a harmless value at index zero. -/
def reductionU26 : ℕ → ℚ
  | 0 => 0
  | n + 1 => recessiveTerm26 (n + 1) * reductionCoefficient26 (n + 1)

theorem reductionU26_of_pos (n : ℕ) (hn : 1 ≤ n) :
    reductionU26 n = recessiveTerm26 n * reductionCoefficient26 n := by
  cases n with
  | zero => omega
  | succ n => rfl

theorem reductionU26_one : reductionU26 1 = -93 / 4480 := by
  rw [reductionU26_of_pos 1 (by omega), recessiveTerm26_one,
    reductionCoefficient26_one, one_mul]

/-- Direct verification of the first-order equation by variation of constants. -/
theorem reductionU26_firstOrder (n : ℕ) (hn : 1 ≤ n) :
    firstOrderResidual26 reductionU26 n = firstOrderForcing26 n := by
  rw [firstOrderResidual26, reductionU26_of_pos n hn,
    reductionU26_of_pos (n + 1) (by omega)]
  have hv := recessiveTerm26_firstOrder n
  have hc := reductionCoefficient26_succ n hn
  have hf := reductionSummand26_forcing n
  calc
    2 * (2 * (n : ℚ) + 7) * ((n : ℚ) + 4) *
          (recessiveTerm26 (n + 1) * reductionCoefficient26 (n + 1)) -
        ((n : ℚ) + 3) ^ 2 *
          (recessiveTerm26 n * reductionCoefficient26 n) =
      ((n : ℚ) + 3) ^ 2 * recessiveTerm26 n *
        (reductionCoefficient26 (n + 1) - reductionCoefficient26 n) := by
          rw [show
            2 * (2 * (n : ℚ) + 7) * ((n : ℚ) + 4) *
                (recessiveTerm26 (n + 1) * reductionCoefficient26 (n + 1)) =
              (2 * (2 * (n : ℚ) + 7) * ((n : ℚ) + 4) *
                recessiveTerm26 (n + 1)) * reductionCoefficient26 (n + 1) by ring]
          rw [hv]
          ring
    _ = ((n : ℚ) + 3) ^ 2 * recessiveTerm26 n *
        (-3 / 1120 * reductionSummand26 (n + 1)) := by rw [hc]
    _ = firstOrderForcing26 n := hf

theorem reductionU26_two : reductionU26 2 = -117 / 14000 := by
  have h := reductionU26_firstOrder 1 (by omega)
  norm_num [firstOrderResidual26, firstOrderForcing26, reductionU26_one] at h ⊢
  linarith

theorem reductionU26_recurrence : SatisfiesRecurrence26 reductionU26 :=
  satisfiesRecurrence26_of_firstOrder reductionU26 reductionU26_firstOrder

/-- The exact finite reduction-of-order formula for the challenge sequence. -/
theorem challengeU26_reduction_formula (n : ℕ) (hn : 1 ≤ n) :
    challengeU26 n = recessiveTerm26 n * reductionCoefficient26 n := by
  have heq := recurrence26_unique reductionU26 challengeU26
    reductionU26_one reductionU26_two challengeU26_one challengeU26_two
    reductionU26_recurrence challengeU26_recurrence (n - 1)
  have hidx : n - 1 + 1 = n := Nat.sub_add_cancel hn
  rw [hidx] at heq
  rw [← heq]
  exact reductionU26_of_pos n hn

/-- After three initial correction terms, the normalized inverse-binomial
sequence is literally the challenge sequence. -/
theorem normalizedQ26_eq_challengeU26 (n : ℕ) (hn : 1 ≤ n) :
    normalizedQ26 (n + 3) = challengeU26 n := by
  rw [challengeU26_reduction_formula n hn,
    recessiveTerm26_eq_inverseBinomialC26]
  unfold normalizedQ26 reductionCoefficient26
  rw [inverseBinomialA26_sum_split n hn]
  ring

/-! ## Main Theorem

2077/720 + Σ_{j=1}^∞ u_j = ζ(2) + ζ(3)

where ζ(2) = π²/6 and ζ(3) is Apéry's constant.
-/

-- ζ(2) = π²/6
theorem zeta2_eq : ∑' n : ℕ, (1 : ℝ) / (↑n + 1) ^ 2 = Real.pi ^ 2 / 6 := by
  calc
    ∑' n : ℕ, (1 : ℝ) / (↑n + 1) ^ 2 = ∑' n : ℕ, (1 : ℝ) / (n : ℝ) ^ 2 := by
      rw [hasSum_zeta_two.summable.tsum_eq_zero_add]
      simp [Nat.cast_add, Nat.cast_one]
    _ = Real.pi ^ 2 / 6 := hasSum_zeta_two.tsum_eq

def zeta2Term26 (n : ℕ) : ℝ :=
  1 / ((n : ℝ) + 1) ^ 2

theorem zeta2Term26_summable : Summable zeta2Term26 := by
  have h :
      Summable (fun n : ℕ => (((n : ℝ) ^ 2)⁻¹)) :=
    Real.summable_nat_pow_inv.mpr (by norm_num)
  have hshift :
      Summable (fun n : ℕ => ((((n + 1 : ℕ) : ℝ) ^ 2)⁻¹)) :=
    (summable_nat_add_iff
      (f := fun n : ℕ => (((n : ℝ) ^ 2)⁻¹)) 1).mpr h
  convert hshift using 1
  funext n
  simp only [zeta2Term26, Nat.cast_add, Nat.cast_one, one_div]

theorem zeta2Term26_hasSum :
    HasSum zeta2Term26 (Real.pi ^ 2 / 6) := by
  apply zeta2Term26_summable.hasSum_iff.mpr
  exact zeta2_eq

def alternatingZeta2Term26 (n : ℕ) : ℝ :=
  (-1 : ℝ) ^ (n + 1) * zeta2Term26 n

/-- The alternating square series is one half of the ordinary zeta-two
series, with the sign convention beginning at `-1`. -/
theorem alternatingZeta2Term26_hasSum :
    HasSum alternatingZeta2Term26
      (-(1 : ℝ) / 2 * (Real.pi ^ 2 / 6)) := by
  have hOdd :
      HasSum (fun k : ℕ => zeta2Term26 (2 * k + 1))
        ((1 : ℝ) / 4 * (Real.pi ^ 2 / 6)) := by
    convert zeta2Term26_hasSum.mul_left ((1 : ℝ) / 4) using 1
    funext k
    unfold zeta2Term26
    push_cast
    have hk : (k : ℝ) + 1 ≠ 0 := by positivity
    field_simp
    ring
  have hEvenSummable :
      Summable (fun k : ℕ => zeta2Term26 (2 * k)) :=
    zeta2Term26_summable.comp_injective
      (mul_right_injective₀ (by omega : (2 : ℕ) ≠ 0))
  have hEvenTsum :
      (∑' k : ℕ, zeta2Term26 (2 * k)) =
        Real.pi ^ 2 / 6 - (1 : ℝ) / 4 * (Real.pi ^ 2 / 6) := by
    have hsplit :=
      tsum_even_add_odd hEvenSummable hOdd.summable
    rw [hOdd.tsum_eq, zeta2Term26_hasSum.tsum_eq] at hsplit
    linarith
  have hEven :
      HasSum (fun k : ℕ => zeta2Term26 (2 * k))
        (Real.pi ^ 2 / 6 - (1 : ℝ) / 4 * (Real.pi ^ 2 / 6)) :=
    hEvenSummable.hasSum_iff.mpr hEvenTsum
  have hAltEven :
      HasSum (fun k : ℕ => alternatingZeta2Term26 (2 * k))
        (-(Real.pi ^ 2 / 6 - (1 : ℝ) / 4 * (Real.pi ^ 2 / 6))) := by
    convert hEven.neg using 1
    funext k
    simp [alternatingZeta2Term26, pow_add, pow_mul]
  have hAltOdd :
      HasSum (fun k : ℕ => alternatingZeta2Term26 (2 * k + 1))
        ((1 : ℝ) / 4 * (Real.pi ^ 2 / 6)) := by
    convert hOdd using 1
    funext k
    simp [alternatingZeta2Term26, pow_add, pow_mul]
  convert HasSum.even_add_odd hAltEven hAltOdd using 1
  ring

theorem dilog_neg_one26 :
    dilog (-1) = -(Real.pi ^ 2) / 12 := by
  have h := alternatingZeta2Term26_hasSum.tsum_eq
  unfold dilog
  calc
    (∑' n : ℕ, (-1 : ℝ) ^ (n + 1) / (↑(n + 1) : ℝ) ^ 2) =
        ∑' n : ℕ, alternatingZeta2Term26 n := by
          apply tsum_congr
          intro n
          simp [alternatingZeta2Term26, zeta2Term26, div_eq_mul_inv]
    _ = -(1 : ℝ) / 2 * (Real.pi ^ 2 / 6) := h
    _ = -(Real.pi ^ 2) / 12 := by ring

/-- An antiderivative for the logarithmic kernel obtained from the beta
representation. -/
def inverseBinomialPrimitive26 (x : ℝ) : ℝ :=
  -dilog (-x) + (1 / 3 : ℝ) * dilog (-(x ^ 3))

theorem inverseBinomialPrimitive26_hasDerivAt
    (x : ℝ) (hx0 : 0 < x) (hx1 : x < 1) :
    HasDerivAt inverseBinomialPrimitive26
      (-Real.log (1 - x + x ^ 2) / x) x := by
  have hxne : x ≠ 0 := ne_of_gt hx0
  have hnegxabs : |-x| < 1 := by
    rw [abs_neg, abs_of_pos hx0]
    exact hx1
  have hnegxne : -x ≠ 0 := neg_ne_zero.mpr hxne
  have hd1base :=
    dilog_hasDerivAt_of_abs_lt_one hnegxabs hnegxne
  have hd1 :
      HasDerivAt (fun y : ℝ => -dilog (-y))
        (Real.log (1 + x) / x) x := by
    have hcomp :=
      hd1base.comp x ((hasDerivAt_id x).neg)
    convert hcomp.neg using 1
    ring
  have hx3pos : 0 < x ^ 3 := pow_pos hx0 3
  have hx3lt : x ^ 3 < 1 :=
    pow_lt_one₀ hx0.le hx1 (by norm_num)
  have hnegx3abs : |-(x ^ 3)| < 1 := by
    rw [abs_neg, abs_of_pos hx3pos]
    exact hx3lt
  have hnegx3ne : -(x ^ 3) ≠ 0 :=
    neg_ne_zero.mpr (ne_of_gt hx3pos)
  have hd3base :=
    dilog_hasDerivAt_of_abs_lt_one hnegx3abs hnegx3ne
  have hd3 :
      HasDerivAt (fun y : ℝ => (1 / 3 : ℝ) * dilog (-(y ^ 3)))
        (-Real.log (1 + x ^ 3) / x) x := by
    have hinner :
        HasDerivAt (fun y : ℝ => -(y ^ 3)) (-(3 * x ^ 2)) x :=
      (hasDerivAt_pow 3 x).neg
    have hcomp := hd3base.comp x hinner
    have hscaled := hcomp.const_mul (1 / 3 : ℝ)
    convert hscaled using 1
    · field_simp [hxne]
      ring
  have hlog :
      Real.log (1 - x + x ^ 2) =
        Real.log (1 + x ^ 3) - Real.log (1 + x) := by
    rw [show 1 - x + x ^ 2 = (1 + x ^ 3) / (1 + x) by
      field_simp
      ring]
    rw [Real.log_div (by positivity) (by positivity)]
  unfold inverseBinomialPrimitive26
  convert hd1.add hd3 using 1
  rw [hlog]
  field_simp
  ring

theorem inverseBinomialPrimitive26_continuousOn :
    ContinuousOn inverseBinomialPrimitive26 (Set.Icc (0 : ℝ) 1) := by
  have hneg :
      ContinuousOn (fun x : ℝ => dilog (-x)) (Set.Icc (0 : ℝ) 1) := by
    apply dilog_continuousOn_unit.comp
    · fun_prop
    · intro x hx
      constructor <;> linarith [hx.1, hx.2]
  have hcube :
      ContinuousOn (fun x : ℝ => dilog (-(x ^ 3))) (Set.Icc (0 : ℝ) 1) := by
    apply dilog_continuousOn_unit.comp
    · fun_prop
    · intro x hx
      have hpow0 : 0 ≤ x ^ 3 := pow_nonneg hx.1 3
      have hpow1 : x ^ 3 ≤ 1 := pow_le_one₀ hx.1 hx.2
      constructor
      · linarith
      · linarith
  unfold inverseBinomialPrimitive26
  fun_prop

theorem inverseBinomialKernel26_tsum_continuousOn :
    ContinuousOn (fun x : ℝ => ∑' n : ℕ, inverseBinomialKernel26 n x)
      (Set.Icc (0 : ℝ) 1) := by
  apply continuousOn_tsum
    (u := fun n : ℕ => ((1 : ℝ) / 4) ^ n)
  · intro n
    exact (inverseBinomialKernel26 n).continuous.continuousOn
  · exact summable_geometric_of_norm_lt_one (by norm_num)
  · intro n x hx
    have hx' : x ∈ Set.uIcc (0 : ℝ) 1 := by
      simpa [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hx
    have hpoint :=
      ContinuousMap.norm_coe_le_norm
        ((inverseBinomialKernel26 n).restrict
          (⟨Set.uIcc (0 : ℝ) 1, isCompact_uIcc⟩ : Compacts ℝ))
        ⟨x, hx'⟩
    exact hpoint.trans (inverseBinomialKernel26_norm_le n)

theorem inverseBinomialLog26_intervalIntegrable :
    IntervalIntegrable
      (fun x : ℝ => -Real.log (1 - x + x ^ 2) / x)
      MeasureTheory.volume 0 1 := by
  have hcont :
      ContinuousOn (fun x : ℝ =>
        ∑' n : ℕ, inverseBinomialKernel26 n x)
        (Set.uIcc (0 : ℝ) 1) := by
    simpa [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using
      inverseBinomialKernel26_tsum_continuousOn
  have hint : IntervalIntegrable
      (fun x : ℝ => ∑' n : ℕ, inverseBinomialKernel26 n x)
      MeasureTheory.volume 0 1 :=
    hcont.intervalIntegrable
  apply IntervalIntegrable.congr
    (f := fun x : ℝ => ∑' n : ℕ, inverseBinomialKernel26 n x) ?_ hint
  intro x hx
  have hx' : x ∈ Set.Ioc (0 : ℝ) 1 := by
    simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hx
  exact inverseBinomialKernel26_tsum_log x hx'.1 hx'.2

theorem inverseBinomialLogIntegral26 :
    (∫ x : ℝ in 0..1, -Real.log (1 - x + x ^ 2) / x) =
      Real.pi ^ 2 / 18 := by
  have hftc :=
    intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le
      (a := (0 : ℝ)) (b := 1)
      (by norm_num)
      inverseBinomialPrimitive26_continuousOn
      (fun x hx =>
        inverseBinomialPrimitive26_hasDerivAt x hx.1 hx.2)
      inverseBinomialLog26_intervalIntegrable
  rw [hftc]
  simp [inverseBinomialPrimitive26, dilog_zero, dilog_neg_one26]
  ring

/-- The classical weight-two inverse central-binomial evaluation, now proved
without an analytic hypothesis. -/
theorem inverseBinomialC26_hasSum :
    HasSum (fun k : ℕ => ((inverseBinomialC26 k : ℚ) : ℝ))
      ((Real.pi ^ 2 / 6) / 3) := by
  convert inverseBinomialC26_hasSum_logIntegral using 1
  rw [inverseBinomialLogIntegral26]
  ring

def zeta3 : ℝ := ∑' n : ℕ, (1 : ℝ) / (↑n + 1) ^ 3

/-- The positive-index summand used in the definition of `zeta3`. -/
def zeta3Term26 (n : ℕ) : ℝ :=
  1 / ((n : ℝ) + 1) ^ 3

theorem zeta3Term26_summable : Summable zeta3Term26 := by
  have h :
      Summable (fun n : ℕ => (((n : ℝ) ^ 3)⁻¹)) :=
    Real.summable_nat_pow_inv.mpr (by norm_num)
  have hshift :
      Summable (fun n : ℕ => ((((n + 1 : ℕ) : ℝ) ^ 3)⁻¹)) :=
    (summable_nat_add_iff
      (f := fun n : ℕ => (((n : ℝ) ^ 3)⁻¹)) 1).mpr h
  convert hshift using 1
  funext n
  simp only [zeta3Term26, Nat.cast_add, Nat.cast_one, one_div]

theorem zeta3Term26_hasSum : HasSum zeta3Term26 zeta3 := by
  unfold zeta3
  exact zeta3Term26_summable.hasSum

/-- The alternating cubic series, indexed so that its first sign is negative:
`-1 + 1/2^3 - 1/3^3 + ⋯`. -/
def alternatingZeta3Term26 (n : ℕ) : ℝ :=
  (-1 : ℝ) ^ (n + 1) * zeta3Term26 n

/-- The alternating cubic series is `-3/4 * ζ(3)`.  This proof only splits
the absolutely convergent zeta series into its even and odd indices. -/
theorem alternatingZeta3Term26_hasSum :
    HasSum alternatingZeta3Term26 (-(3 : ℝ) / 4 * zeta3) := by
  have hOdd :
      HasSum (fun k : ℕ => zeta3Term26 (2 * k + 1)) ((1 : ℝ) / 8 * zeta3) := by
    convert zeta3Term26_hasSum.mul_left ((1 : ℝ) / 8) using 1
    · funext k
      unfold zeta3Term26
      push_cast
      have hk : (k : ℝ) + 1 ≠ 0 := by positivity
      field_simp
      ring
  have hEvenSummable :
      Summable (fun k : ℕ => zeta3Term26 (2 * k)) :=
    zeta3Term26_summable.comp_injective
      (mul_right_injective₀ (by omega : (2 : ℕ) ≠ 0))
  have hEvenTsum :
      (∑' k : ℕ, zeta3Term26 (2 * k)) =
        zeta3 - (1 : ℝ) / 8 * zeta3 := by
    have hsplit :=
      tsum_even_add_odd hEvenSummable hOdd.summable
    rw [hOdd.tsum_eq, zeta3Term26_hasSum.tsum_eq] at hsplit
    linarith
  have hEven :
      HasSum (fun k : ℕ => zeta3Term26 (2 * k))
        (zeta3 - (1 : ℝ) / 8 * zeta3) :=
    hEvenSummable.hasSum_iff.mpr hEvenTsum
  have hAltEven :
      HasSum (fun k : ℕ => alternatingZeta3Term26 (2 * k))
        (-(zeta3 - (1 : ℝ) / 8 * zeta3)) := by
    convert hEven.neg using 1
    funext k
    simp [alternatingZeta3Term26, pow_add, pow_mul]
  have hAltOdd :
      HasSum (fun k : ℕ => alternatingZeta3Term26 (2 * k + 1))
        ((1 : ℝ) / 8 * zeta3) := by
    convert hOdd using 1
    funext k
    simp [alternatingZeta3Term26, pow_add, pow_mul]
  convert HasSum.even_add_odd hAltEven hAltOdd using 1
  ring

/-- The terms of `ζ(3)` whose positive index is a multiple of three. -/
def everyThirdZeta3Term26 (n : ℕ) : ℝ :=
  if n % 3 = 2 then zeta3Term26 n else 0

theorem everyThirdZeta3Term26_hasSum :
    HasSum everyThirdZeta3Term26 ((1 : ℝ) / 27 * zeta3) := by
  let g : ℕ → ℕ := fun k => 3 * k + 2
  have hg : Function.Injective g := by
    intro k l hkl
    dsimp [g] at hkl
    omega
  have hcomp :
      HasSum (everyThirdZeta3Term26 ∘ g) ((1 : ℝ) / 27 * zeta3) := by
    convert zeta3Term26_hasSum.mul_left ((1 : ℝ) / 27) using 1
    funext k
    simp only [Function.comp_apply, g, everyThirdZeta3Term26,
      show (3 * k + 2) % 3 = 2 by omega, ↓reduceIte, zeta3Term26]
    push_cast
    have hk : (k : ℝ) + 1 ≠ 0 := by positivity
    field_simp
    ring
  apply (hg.hasSum_iff ?_).mp hcomp
  intro n hn
  unfold everyThirdZeta3Term26
  split_ifs with h
  · exfalso
    apply hn
    refine ⟨n / 3, ?_⟩
    dsimp [g]
    omega
  · rfl

/-- The alternating cubic terms whose positive index is a multiple of three. -/
def everyThirdAlternatingZeta3Term26 (n : ℕ) : ℝ :=
  if n % 3 = 2 then alternatingZeta3Term26 n else 0

theorem everyThirdAlternatingZeta3Term26_hasSum :
    HasSum everyThirdAlternatingZeta3Term26
      ((1 : ℝ) / 27 * (-(3 : ℝ) / 4 * zeta3)) := by
  let g : ℕ → ℕ := fun k => 3 * k + 2
  have hg : Function.Injective g := by
    intro k l hkl
    dsimp [g] at hkl
    omega
  have hcomp :
      HasSum (everyThirdAlternatingZeta3Term26 ∘ g)
        ((1 : ℝ) / 27 * (-(3 : ℝ) / 4 * zeta3)) := by
    convert alternatingZeta3Term26_hasSum.mul_left ((1 : ℝ) / 27) using 1
    funext k
    simp only [Function.comp_apply, g, everyThirdAlternatingZeta3Term26,
      show (3 * k + 2) % 3 = 2 by omega, ↓reduceIte,
      alternatingZeta3Term26, zeta3Term26]
    push_cast
    have hk : (k : ℝ) + 1 ≠ 0 := by positivity
    rw [show 3 * k + 2 + 1 = 3 * (k + 1) by omega, pow_mul]
    norm_num
    field_simp
    ring
  apply (hg.hasSum_iff ?_).mp hcomp
  intro n hn
  unfold everyThirdAlternatingZeta3Term26
  split_ifs with h
  · exfalso
    apply hn
    refine ⟨n / 3, ?_⟩
    dsimp [g]
    omega
  · rfl

/-- The real cubic polylogarithm weight at a primitive cube root:
`1` on multiples of three and `-1/2` otherwise. -/
def omegaCubicTerm26 (n : ℕ) : ℝ :=
  (3 : ℝ) / 2 * everyThirdZeta3Term26 n -
    (1 : ℝ) / 2 * zeta3Term26 n

theorem omegaCubicTerm26_hasSum :
    HasSum omegaCubicTerm26 (-(4 : ℝ) / 9 * zeta3) := by
  have h :=
    (everyThirdZeta3Term26_hasSum.mul_left ((3 : ℝ) / 2)).sub
      (zeta3Term26_hasSum.mul_left ((1 : ℝ) / 2))
  have h' :
      HasSum omegaCubicTerm26
        ((3 : ℝ) / 2 * ((1 : ℝ) / 27 * zeta3) -
          (1 : ℝ) / 2 * zeta3) := by
    simpa only [omegaCubicTerm26] using h
  convert h' using 1
  ring

/-- The real cubic polylogarithm weight at a primitive sixth root. -/
def rhoCubicTerm26 (n : ℕ) : ℝ :=
  (3 : ℝ) / 2 * everyThirdAlternatingZeta3Term26 n -
    (1 : ℝ) / 2 * alternatingZeta3Term26 n

theorem rhoCubicTerm26_hasSum :
    HasSum rhoCubicTerm26 ((1 : ℝ) / 3 * zeta3) := by
  have h :=
    (everyThirdAlternatingZeta3Term26_hasSum.mul_left ((3 : ℝ) / 2)).sub
      (alternatingZeta3Term26_hasSum.mul_left ((1 : ℝ) / 2))
  have h' :
      HasSum rhoCubicTerm26
        ((3 : ℝ) / 2 * ((1 : ℝ) / 27 * (-(3 : ℝ) / 4 * zeta3)) -
          (1 : ℝ) / 2 * (-(3 : ℝ) / 4 * zeta3)) := by
    simpa only [rhoCubicTerm26] using h
  convert h' using 1
  ring

def reductionSeriesTerm26 (n : ℕ) : ℝ :=
  ((recessiveTerm26 (n + 1) * reductionCoefficient26 (n + 1) : ℚ) : ℝ)

/-- The normalized inverse-binomial series is summable.  Its tail after four
terms is exactly the absolutely summable challenge series. -/
theorem normalizedQ26_summable :
    Summable (fun k : ℕ => ((normalizedQ26 k : ℚ) : ℝ)) := by
  have htail :
      Summable (fun n : ℕ => ((normalizedQ26 (n + 4) : ℚ) : ℝ)) := by
    convert challengeU26_summable using 1
    funext n
    exact_mod_cast normalizedQ26_eq_challengeU26 (n + 1) (by omega)
  exact (summable_nat_add_iff
    (f := fun k : ℕ => ((normalizedQ26 k : ℚ) : ℝ)) 4).mp htail

/-- The single inverse-binomial series is absolutely summable, by comparison
with the ordinary inverse-square series. -/
theorem inverseBinomialC26_summable :
    Summable (fun k : ℕ => ((inverseBinomialC26 k : ℚ) : ℝ)) := by
  have hmajorant :
      Summable (fun k : ℕ => (((k : ℝ) ^ 2)⁻¹)) :=
    Real.summable_nat_pow_inv.mpr (by norm_num)
  apply Summable.of_norm_bounded hmajorant
  intro k
  rw [Real.norm_eq_abs]
  by_cases hk : k = 0
  · simp [hk, inverseBinomialC26]
  have hkpos : 0 < k := Nat.pos_of_ne_zero hk
  have hchooseNat : 1 ≤ Nat.choose (2 * k) k :=
    Nat.choose_pos (by omega)
  have hchoose : (1 : ℝ) ≤ (Nat.choose (2 * k) k : ℝ) := by
    exact_mod_cast hchooseNat
  have hkR : (0 : ℝ) < k := by exact_mod_cast hkpos
  have hnonneg : (0 : ℝ) ≤ ((inverseBinomialC26 k : ℚ) : ℝ) := by
    simp only [inverseBinomialC26]
    positivity
  rw [abs_of_nonneg hnonneg]
  simp only [inverseBinomialC26]
  push_cast
  rw [one_div]
  apply inv_anti₀ (sq_pos_of_pos hkR)
  exact le_mul_of_one_le_right (sq_nonneg (k : ℝ)) hchoose

/-- The nested inverse-binomial series is also absolutely summable.  This is
deduced algebraically from the already summable normalized challenge series
and the single inverse-binomial series. -/
theorem inverseBinomialDTerm26_summable :
    Summable (fun k : ℕ => ((inverseBinomialDTerm26 k : ℚ) : ℝ)) := by
  have h :=
    (inverseBinomialC26_summable.mul_left (2 : ℝ)).sub
      (normalizedQ26_summable.mul_left ((1 : ℝ) / 3))
  convert h using 1
  funext k
  norm_cast
  rw [normalizedQ26_decomposition]
  push_cast
  ring

/-- The challenge series is termwise the explicit nested central-binomial
series obtained by reduction of order. -/
theorem challengeU26_series_term (n : ℕ) :
    ((challengeU26 (n + 1) : ℚ) : ℝ) = reductionSeriesTerm26 n := by
  unfold reductionSeriesTerm26
  exact_mod_cast challengeU26_reduction_formula (n + 1) (by omega)

theorem challengeU26_tsum_eq_reduction :
    (∑' n : ℕ, ((challengeU26 (n + 1) : ℚ) : ℝ)) =
      ∑' n : ℕ, reductionSeriesTerm26 n :=
  tsum_congr challengeU26_series_term

/-- The exact assertion posed by Problem 2.6. -/
def Problem26Claim : Prop :=
  (2077 : ℝ) / 720 +
      (∑' n : ℕ, ((challengeU26 (n + 1) : ℚ) : ℝ)) =
    Real.pi ^ 2 / 6 + zeta3

/-- The challenge claim is exactly the evaluation of the normalized
inverse-binomial series.  This theorem performs all finite-prefix and index
bookkeeping; it does not assume the challenge statement itself. -/
theorem problem26_of_normalizedQ_hasSum
    (hEval :
      HasSum (fun k : ℕ => ((normalizedQ26 k : ℚ) : ℝ))
        (Real.pi ^ 2 / 6 + zeta3)) :
    Problem26Claim := by
  have htail (n : ℕ) :
      ((challengeU26 (n + 1) : ℚ) : ℝ) =
        ((normalizedQ26 (n + 4) : ℚ) : ℝ) := by
    exact_mod_cast (normalizedQ26_eq_challengeU26 (n + 1) (by omega)).symm
  unfold Problem26Claim
  calc
    (2077 : ℝ) / 720 +
          (∑' n : ℕ, ((challengeU26 (n + 1) : ℚ) : ℝ)) =
        (∑ k ∈ Finset.range 4, ((normalizedQ26 k : ℚ) : ℝ)) +
          ∑' n : ℕ, ((normalizedQ26 (n + 4) : ℚ) : ℝ) := by
            rw [normalizedQ26_range_four_real]
            congr 1
            exact tsum_congr htail
    _ = ∑' k : ℕ, ((normalizedQ26 k : ℚ) : ℝ) :=
      hEval.summable.sum_add_tsum_nat_add 4
    _ = Real.pi ^ 2 / 6 + zeta3 := hEval.tsum_eq

/-- It remains enough to evaluate two classical inverse-binomial series:
the weight-two single sum and the weight-three nested sum. -/
theorem problem26_of_inverse_binomial_evaluations
    (hC :
      HasSum (fun k : ℕ => ((inverseBinomialC26 k : ℚ) : ℝ))
        ((Real.pi ^ 2 / 6) / 3))
    (hD :
      HasSum (fun k : ℕ => ((inverseBinomialDTerm26 k : ℚ) : ℝ))
        ((Real.pi ^ 2 / 6 - zeta3) / 3)) :
    Problem26Claim := by
  have hEval :
      HasSum (fun k : ℕ => ((normalizedQ26 k : ℚ) : ℝ))
        (Real.pi ^ 2 / 6 + zeta3) := by
    have h :=
      (hC.mul_left (6 : ℝ)).sub (hD.mul_left (3 : ℝ))
    have h' :
        HasSum (fun k : ℕ => ((normalizedQ26 k : ℚ) : ℝ))
          (6 * ((Real.pi ^ 2 / 6) / 3) -
            3 * ((Real.pi ^ 2 / 6 - zeta3) / 3)) := by
      apply h.congr_fun
      intro k
      exact_mod_cast normalizedQ26_decomposition k
    convert h' using 1
    ring
  exact problem26_of_normalizedQ_hasSum hEval

/-- After the unconditional weight-two evaluation above, only the nested
weight-three inverse-binomial sum remains as an analytic input. -/
theorem problem26_of_nested_inverse_binomial_evaluation
    (hD :
      HasSum (fun k : ℕ => ((inverseBinomialDTerm26 k : ℚ) : ℝ))
        ((Real.pi ^ 2 / 6 - zeta3) / 3)) :
    Problem26Claim :=
  problem26_of_inverse_binomial_evaluations
    inverseBinomialC26_hasSum hD

/-- A precise description of the remaining analytic input: evaluation of the
explicit reduction-of-order series proves the challenge statement. -/
theorem problem26_of_reduction_series_evaluation
    (hEval :
      (2077 : ℝ) / 720 + (∑' n : ℕ, reductionSeriesTerm26 n) =
        Real.pi ^ 2 / 6 + zeta3) :
    Problem26Claim := by
  unfold Problem26Claim
  rw [challengeU26_tsum_eq_reduction]
  exact hEval

/-! ## Status of the main statement

The exact challenge assertion is `Problem26Claim`.  Its unconditional proof
is assembled in `Problem26Assembly` after the nested weight-three evaluation.

An earlier version carried a theorem of the shape

```
theorem problem26_identity : ∃ (p q : ℕ → ℝ), Tendsto (fun n => p n / q n) atTop (𝓝 L) :=
  ⟨fun _ => L, fun _ => 1, by simp⟩
```

which is vacuous: it is witnessed by constant sequences and says nothing about
the challenge recurrence.  It has been removed rather than shipped.  Here the
actual challenge sequence is defined and proved unique, while the exact Ore
factorization and finite reduction-of-order formula reduce the problem to the
single nested weight-three evaluation used by
`problem26_of_nested_inverse_binomial_evaluation`.
-/

end RamanujanChallenge.P26

end
