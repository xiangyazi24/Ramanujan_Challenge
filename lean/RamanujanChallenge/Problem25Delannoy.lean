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

theorem delannoyB_diag_pos (k : ℕ) : 0 < delannoyB k k := by
  exact delannoyB_pos k k le_rfl

theorem delannoyB_diag_ne_zero (k : ℕ) : delannoyB k k ≠ 0 :=
  ne_of_gt (delannoyB_diag_pos k)

/-- Triangular inversion: given a sequence s, extract the unique coefficients
    c(k) such that s(N) = Σ_{k≤N} c(k) · B(N,k).
    We use an auxiliary function to avoid well-foundedness issues with Finset.sum. -/
/-- Accumulator: computes (delannoyCoeff s 0, ..., delannoyCoeff s (n-1)). -/
private noncomputable def delannoyCoeffAux (s : ℕ → ℚ) : ℕ → List ℚ
  | 0 => []
  | n + 1 =>
      let prev := delannoyCoeffAux s n
      let cn := (s n - (prev.enum.map fun ⟨j, cj⟩ => cj * delannoyB n j).sum) / delannoyB n n
      prev ++ [cn]

noncomputable def delannoyCoeff (s : ℕ → ℚ) (k : ℕ) : ℚ :=
  (delannoyCoeffAux s (k + 1)).getLast (by simp [delannoyCoeffAux])

end RamanujanChallenge.P25

end
