/-
  Problem 2.5 — Delannoy basis layer.

  The Clausen–Delannoy summand B(N,k) = 2^k · C(2k,k) · C(N,k) · C(N+k,k)
  forms a lower-triangular basis for sequences indexed by N.
  The CMF sequences Q̂_N and P̂_N decompose uniquely in this basis:
    Q̂_N = Σ_{k=0}^N f(k) · B(N,k)
    P̂_N = Σ_{k=0}^N g(k) · B(N,k)

  The positivity f(k) > 0 and the convergence g(k)/f(k) → G at geometric
  rate -1/8 are the two inputs to the splitting argument that closes P2.5.
-/
import RamanujanChallenge.Problem25
import Mathlib.Data.Nat.Choose.Factorization

noncomputable section

namespace RamanujanChallenge.P25

open Nat Finset

/-- The Clausen–Delannoy summand. -/
def delannoyB (N k : ℕ) : ℚ :=
  2 ^ k * (Nat.choose (2 * k) k : ℚ) * (Nat.choose N k : ℚ) * (Nat.choose (N + k) k : ℚ)

theorem delannoyB_nonneg (N k : ℕ) : 0 ≤ delannoyB N k := by
  unfold delannoyB; positivity

theorem delannoyB_pos (N k : ℕ) (hk : k ≤ N) : 0 < delannoyB N k := by
  unfold delannoyB
  apply mul_pos (mul_pos (mul_pos _ _) _) _
  · exact pow_pos (by norm_num : (0:ℚ) < 2) k
  · exact Nat.cast_pos.mpr (Nat.choose_pos (by omega))
  · exact Nat.cast_pos.mpr (Nat.choose_pos hk)
  · exact Nat.cast_pos.mpr (Nat.choose_pos (by omega))

theorem delannoyB_eq_zero_of_lt (N k : ℕ) (hk : N < k) : delannoyB N k = 0 := by
  unfold delannoyB
  have : Nat.choose N k = 0 := Nat.choose_eq_zero_of_lt hk
  simp [this]

@[simp] theorem delannoyB_zero_zero : delannoyB 0 0 = 1 := by
  simp [delannoyB]

theorem delannoyB_diag (k : ℕ) :
    delannoyB k k = 2 ^ k * (Nat.choose (2 * k) k : ℚ) ^ 2 := by
  unfold delannoyB
  rw [Nat.choose_self, show k + k = 2 * k by omega]
  push_cast
  ring

theorem delannoyB_diag_pos (k : ℕ) : 0 < delannoyB k k := by
  exact delannoyB_pos k k le_rfl

theorem delannoyB_diag_ne_zero (k : ℕ) : delannoyB k k ≠ 0 :=
  ne_of_gt (delannoyB_diag_pos k)

/-- Triangular inversion: given a sequence `s`, extract the coefficient of
`delannoyB k k`.  Indexing the preceding sum by `Fin k` exposes the strict
inequality required by the termination checker. -/
def delannoyCoeff (s : ℕ → ℚ) (k : ℕ) : ℚ :=
  (s k - ∑ j : Fin k, delannoyCoeff s j * delannoyB k j) / delannoyB k k
  termination_by k
  decreasing_by exact j.isLt

/-- Every rational sequence has a Delannoy decomposition. -/
theorem delannoy_decomposition (s : ℕ → ℚ) (N : ℕ) :
    s N = ∑ k ∈ range (N + 1), delannoyCoeff s k * delannoyB N k := by
  rw [sum_range_succ, delannoyCoeff, ← Fin.sum_univ_eq_sum_range]
  rw [div_mul_cancel₀ _ (delannoyB_diag_ne_zero N)]
  ring

/-- The Delannoy coefficients are uniquely determined by the decomposition. -/
theorem delannoyCoeff_unique (s c : ℕ → ℚ)
    (hc : ∀ N, s N = ∑ k ∈ range (N + 1), c k * delannoyB N k) :
    c = delannoyCoeff s := by
  funext k
  induction k using Nat.strong_induction_on with
  | h k ih =>
      have hc' := hc k
      have hd := delannoy_decomposition s k
      rw [sum_range_succ] at hc' hd
      have hprev :
          ∑ x ∈ range k, c x * delannoyB k x =
            ∑ x ∈ range k, delannoyCoeff s x * delannoyB k x := by
        apply sum_congr rfl
        intro x hx
        rw [ih x (mem_range.mp hx)]
      rw [hprev] at hc'
      have hmul : c k * delannoyB k k = delannoyCoeff s k * delannoyB k k := by
        linarith
      exact mul_right_cancel₀ (delannoyB_diag_ne_zero k) hmul

/-! ## The normalized CMF sequences -/

/-- First denominator column after division by the cumulative Pochhammer
gauge used by the numerical Delannoy decomposition. -/
def normalizedDenominator (N : ℕ) : ℚ :=
  (denominator N 0 : ℚ) / (pochhammerProduct N : ℚ)

/-- First numerator column after division by the cumulative Pochhammer
gauge used by the numerical Delannoy decomposition. -/
def normalizedNumerator (N : ℕ) : ℚ :=
  (numerator N 0 : ℚ) / (pochhammerProduct N : ℚ)

/-- Delannoy coefficient sequence for the normalized denominator. -/
def delannoyF (k : ℕ) : ℚ :=
  delannoyCoeff normalizedDenominator k

/-- Delannoy coefficient sequence for the normalized numerator. -/
def delannoyG (k : ℕ) : ℚ :=
  delannoyCoeff normalizedNumerator k

theorem normalizedDenominator_decomposition (N : ℕ) :
    normalizedDenominator N =
      ∑ k ∈ range (N + 1), delannoyF k * delannoyB N k := by
  exact delannoy_decomposition normalizedDenominator N

theorem normalizedNumerator_decomposition (N : ℕ) :
    normalizedNumerator N =
      ∑ k ∈ range (N + 1), delannoyG k * delannoyB N k := by
  exact delannoy_decomposition normalizedNumerator N

theorem delannoyF_unique (c : ℕ → ℚ)
    (hc : ∀ N, normalizedDenominator N =
      ∑ k ∈ range (N + 1), c k * delannoyB N k) :
    c = delannoyF := by
  exact delannoyCoeff_unique normalizedDenominator c hc

theorem delannoyG_unique (c : ℕ → ℚ)
    (hc : ∀ N, normalizedNumerator N =
      ∑ k ∈ range (N + 1), c k * delannoyB N k) :
    c = delannoyG := by
  exact delannoyCoeff_unique normalizedNumerator c hc

@[simp] theorem delannoyF_zero : delannoyF 0 = 33750 := by
  rw [delannoyF, delannoyCoeff]
  norm_num [normalizedDenominator, denominator, approximants, initialMatrix,
    delannoyB]

@[simp] theorem delannoyG_zero : delannoyG 0 = 30921 := by
  rw [delannoyG, delannoyCoeff]
  norm_num [normalizedNumerator, numerator, approximants, initialMatrix,
    delannoyB]

theorem delannoyF_zero_pos : 0 < delannoyF 0 := by norm_num

end RamanujanChallenge.P25

end
