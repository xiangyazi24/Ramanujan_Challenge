/-
  Problem 3.2 — Layer 1: Apéry Sequence Definitions.

  The Apéry sequences (a_n), (b_n) satisfy the three-term recurrence:
    (n+1)³ u_{n+1} = P(n) u_n - n³ u_{n-1},  n ≥ 1
  where P(n) = 34n³ + 51n² + 27n + 5.

  Initial values: a_0 = 0, a_1 = 6, b_0 = 1, b_1 = 5.

  b_n = Σ_{k=0}^n C(n,k)² C(n+k,k)² ∈ ℤ  (the Apéry numbers).
  a_n ∈ ℚ in general; d_n · a_n ∈ ℤ where d_n = lcm(1,...,n)³.

  G_n = gcd(d_n · a_n, d_n · b_n).
  The conjecture: G_n = e^{o(n)}, i.e., log G_n = o(n).

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

def aperyMiddle (n : ℤ) : ℤ := 34 * n ^ 3 + 51 * n ^ 2 + 27 * n + 5

theorem aperyMiddle_zero : aperyMiddle 0 = 5 := by unfold aperyMiddle; ring

theorem aperyMiddle_one : aperyMiddle 1 = 117 := by unfold aperyMiddle; ring

/-! ## Apéry numbers b_n (closed form) -/

def aperyB (n : ℕ) : ℤ :=
  ∑ k ∈ Finset.range (n + 1),
    (Nat.choose n k : ℤ) ^ 2 * (Nat.choose (n + k) k : ℤ) ^ 2

theorem aperyB_zero : aperyB 0 = 1 := by
  simp [aperyB]

theorem aperyB_one : aperyB 1 = 5 := by
  unfold aperyB
  simp [Finset.sum_range_succ]

theorem aperyB_two : aperyB 2 = 73 := by
  unfold aperyB
  simp [Finset.sum_range_succ, Nat.choose]

/-! ## Apéry companion sequence a_n (over ℚ, via recurrence) -/

def aperyA : ℕ → ℚ
  | 0 => 0
  | 1 => 6
  | (n + 2) =>
    let m : ℚ := ↑(n + 1)
    ((34 * m ^ 3 + 51 * m ^ 2 + 27 * m + 5) * aperyA (n + 1) -
     m ^ 3 * aperyA n) / (m + 1) ^ 3

theorem aperyA_zero : aperyA 0 = 0 := rfl

theorem aperyA_one : aperyA 1 = 6 := rfl

/-! ## The Wronskian W_n = a_n · b_{n-1} - a_{n-1} · b_n

The Casorati determinant satisfies W_n = 6/n³ for all n ≥ 1.
-/

def wronskian (n : ℕ) : ℚ :=
  aperyA n * ↑(aperyB (n - 1)) - aperyA (n - 1) * ↑(aperyB n)

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
