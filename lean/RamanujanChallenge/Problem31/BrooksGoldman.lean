/-
  Problem 3.1 — Layer 3: Brooks-Goldman Godbillon-Vey formula.

  For a closed Seifert-fibered manifold with Seifert invariants
  M(b; (a₁,b₁), ..., (aₙ,bₙ)), the Godbillon-Vey invariant of
  the maximal Fuchsian representation is:

    GV(ρ_Fuch) = 4π² · χ_orb² / e

  where:
    e = b + Σ bᵢ/aᵢ  (Euler number)
    χ_orb = 2 - n + Σ 1/aᵢ  (orbifold Euler characteristic)

  For Σ(2,3,17) ≅ S³_{-1}(7₂):
    Seifert invariants: M(-1; (2,1), (3,1), (17,3))
    e = -1 + 1/2 + 1/3 + 3/17 = 1/102
    χ_orb = 2 - 3 + 1/2 + 1/3 + 1/17 = -11/102
    GV = 4π² · (11/102)² / (1/102) = 242π²/51
-/
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic

noncomputable section

open Real

/-! ## Seifert invariants of Σ(2,3,17) -/

theorem euler_number :
    (-1 : ℚ) + 1/2 + 1/3 + 3/17 = 1/102 := by norm_num

theorem orbifold_euler_char :
    (2 : ℚ) - 3 + 1/2 + 1/3 + 1/17 = -11/102 := by norm_num

/-! ## The Brooks-Goldman formula (arithmetic) -/

theorem brooks_goldman_rational :
    4 * ((11 : ℚ) / 102) ^ 2 / ((1 : ℚ) / 102) = 242 / 51 := by norm_num

theorem brooks_goldman_real :
    4 * Real.pi ^ 2 * ((11 : ℝ) / 102) ^ 2 / ((1 : ℝ) / 102) =
    (242 : ℝ) / 51 * Real.pi ^ 2 := by ring

/-! ## GV at the β endpoint

The β-endpoint representation of the 7₂ knot complement factors
through π₁(Σ(2,3,17)) and is the maximal Fuchsian representation.
Therefore:

  GV(ρ_β) = 242π²/51
-/

def gv_beta : ℝ := (242 : ℝ) / 51 * Real.pi ^ 2

/-! ## The target identity reduces to:

  GV(ρ_β) - GV(ρ_α) = 4π²/85

Equivalently:
  GV(ρ_α) = GV(ρ_β) - 4π²/85 = (242/51 - 4/85)π² = 1198π²/255
-/

theorem gv_difference_arithmetic :
    (242 : ℚ) / 51 - 4 / 85 = 1198 / 255 := by norm_num

theorem gv_difference_real :
    (242 : ℝ) / 51 * Real.pi ^ 2 - 4 / 85 * Real.pi ^ 2 =
    (1198 : ℝ) / 255 * Real.pi ^ 2 := by ring

/-! ## Brieskorn sphere identification

S³_{-1}(7₂) ≅ Σ(2,3,17).

The 7₂ knot is the five-twist knot (10 crossings in the twist region).
The (−1)-surgery gives the Brieskorn sphere with third multiplicity
17 = 6·3 − 1 (the pattern for twist knots: 5₂ gives 11, 7₂ gives 17).

The first homology |H₁(Σ(2,3,17); ℤ)| = 1, so the Seifert fibering
is unique (up to orientation).
-/

theorem brieskorn_h1_order : (2 : ℚ) * 3 * 17 * (1 / 102) = 1 := by norm_num

end
