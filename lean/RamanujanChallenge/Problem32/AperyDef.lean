/-
  Problem 3.2 — Layer 1: Apéry Sequence Definitions.

  The Apéry sequences (a_n), (b_n) satisfy the three-term recurrence:
    (n+1)³ u_{n+1} = P(n) u_n - n³ u_{n-1},  n ≥ 1
  where P(n) = 34n³ + 51n² + 27n + 5.

  Initial values: a_0 = 0, a_1 = 6, b_0 = 1, b_1 = 5.

  We define BOTH sequences by recurrence over ℚ. The Wronskian identity
  W_n = 6/n³ then follows directly from the shared recurrence structure.

  The closed form b_n = Σ C(n,k)² C(n+k,k)² ∈ ℤ and the integrality
  d_n · a_n ∈ ℤ are stated as separate theorems.

  Reference: Xiang Huang, "On Ramanujan Challenge Problem 3.2:
  The Apéry GCD Conjecture", July 2026.
-/
import Mathlib.Data.Nat.Choose.Basic
import Mathlib.Data.Nat.Factorial.Basic
import Mathlib.Data.Rat.Defs
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Tactic

noncomputable section

open Finset

/-! ## The middle coefficient polynomial -/

def aperyP (n : ℚ) : ℚ := 34 * n ^ 3 + 51 * n ^ 2 + 27 * n + 5

theorem aperyP_zero : aperyP 0 = 5 := by unfold aperyP; ring

theorem aperyP_one : aperyP 1 = 117 := by unfold aperyP; ring

/-! ## Apéry recurrence step

Given u_n, u_{n-1}, compute u_{n+1} via:
  (n+1)³ u_{n+1} = P(n) u_n - n³ u_{n-1}
-/

def aperyStep (n : ℕ) (u_n u_nm1 : ℚ) : ℚ :=
  (aperyP (n : ℚ) * u_n - (n : ℚ) ^ 3 * u_nm1) / ((n : ℚ) + 1) ^ 3

/-! ## Both sequences via recurrence -/

def aperyA : ℕ → ℚ
  | 0 => 0
  | 1 => 6
  | (n + 2) => aperyStep (n + 1) (aperyA (n + 1)) (aperyA n)

def aperyBQ : ℕ → ℚ
  | 0 => 1
  | 1 => 5
  | (n + 2) => aperyStep (n + 1) (aperyBQ (n + 1)) (aperyBQ n)

theorem aperyA_zero : aperyA 0 = 0 := rfl
theorem aperyA_one : aperyA 1 = 6 := rfl

theorem aperyBQ_zero : aperyBQ 0 = 1 := rfl
theorem aperyBQ_one : aperyBQ 1 = 5 := rfl

/-! ## Recurrence identity (by definition)

Since both sequences are DEFINED by the recurrence, the recurrence
identity holds definitionally.
-/

theorem aperyA_recurrence (n : ℕ) (hn : n ≥ 1) :
    ((n : ℚ) + 1) ^ 3 * aperyA (n + 1) =
      aperyP (n : ℚ) * aperyA n - (n : ℚ) ^ 3 * aperyA (n - 1) := by
  cases n with
  | zero => omega
  | succ m =>
    simp only [aperyA, aperyStep, Nat.succ_sub_one]
    have hne : ((m : ℚ) + 1 + 1) ^ 3 ≠ 0 := by positivity
    field_simp

theorem aperyBQ_recurrence (n : ℕ) (hn : n ≥ 1) :
    ((n : ℚ) + 1) ^ 3 * aperyBQ (n + 1) =
      aperyP (n : ℚ) * aperyBQ n - (n : ℚ) ^ 3 * aperyBQ (n - 1) := by
  cases n with
  | zero => omega
  | succ m =>
    simp only [aperyBQ, aperyStep, Nat.succ_sub_one]
    have hne : ((m : ℚ) + 1 + 1) ^ 3 ≠ 0 := by positivity
    field_simp

/-! ## The Wronskian -/

def wronskian (n : ℕ) : ℚ :=
  aperyA n * aperyBQ (n - 1) - aperyA (n - 1) * aperyBQ n

/-! ## Apéry numbers b_n (closed form, integer) -/

def aperyB (n : ℕ) : ℤ :=
  ∑ k ∈ Finset.range (n + 1),
    (Nat.choose n k : ℤ) ^ 2 * (Nat.choose (n + k) k : ℤ) ^ 2

theorem aperyB_zero : aperyB 0 = 1 := by simp [aperyB]

theorem aperyB_one : aperyB 1 = 5 := by
  unfold aperyB
  simp [Finset.sum_range_succ]

theorem aperyB_two : aperyB 2 = 73 := by
  unfold aperyB
  simp [Finset.sum_range_succ, Nat.choose]

/-! ## The zero-count function Z(p) -/

def zeroCountApery (p : ℕ) : ℕ :=
  ((Finset.range p).filter fun j => (aperyB j) % (p : ℤ) = 0).card

/-! ## d_n = lcm(1, ..., n)³ -/

def lcmUpTo : ℕ → ℕ
  | 0 => 1
  | (n + 1) => Nat.lcm (lcmUpTo n) (n + 1)

def cubed_lcm (n : ℕ) : ℕ := (lcmUpTo n) ^ 3

theorem lcmUpTo_zero : lcmUpTo 0 = 1 := rfl

theorem lcmUpTo_succ (n : ℕ) : lcmUpTo (n + 1) = Nat.lcm (lcmUpTo n) (n + 1) := rfl

end
