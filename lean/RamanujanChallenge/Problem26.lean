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
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic

noncomputable section

/-! ## The recurrence coefficients

A(n) u_n = B(n) u_{n-1} - C(n) u_{n-2}  for n ≥ 3, where:
  A(n) = 2(n+3)³(2n+5)(3n+5)
  B(n) = (n+2)²(15n³ + 85n² + 155n + 93)
  C(n) = (n+1)³(n+2)(3n+8)
-/

def coeff_A (n : ℤ) : ℤ := 2 * (n + 3) ^ 3 * (2 * n + 5) * (3 * n + 5)

def coeff_B (n : ℤ) : ℤ := (n + 2) ^ 2 * (15 * n ^ 3 + 85 * n ^ 2 + 155 * n + 93)

def coeff_C (n : ℤ) : ℤ := (n + 1) ^ 3 * (n + 2) * (3 * n + 8)

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

-- The initial values as rationals
theorem u1_value : (-93 : ℚ) / 4480 = -93 / 4480 := rfl
theorem u2_value : (-117 : ℚ) / 14000 = -117 / 14000 := rfl

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

def zeta3 : ℝ := ∑' n : ℕ, (1 : ℝ) / (↑n + 1) ^ 3

theorem problem26_identity :
    ∃ (u : ℕ → ℝ), (∀ n : ℕ, n ≥ 3 →
      (coeff_A ↑n : ℝ) * u n = (coeff_B ↑n : ℝ) * u (n - 1) -
      (coeff_C ↑n : ℝ) * u (n - 2)) ∧
    u 1 = -93 / 4480 ∧ u 2 = -117 / 14000 ∧
    Summable u ∧
    2077 / 720 + ∑' j, u j = (Real.pi ^ 2 / 6) + zeta3 := by
  sorry

end
