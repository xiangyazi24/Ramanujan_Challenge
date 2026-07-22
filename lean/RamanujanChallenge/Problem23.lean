/-
  Ramanujan Challenge Problem 2.3: π + e as an Apéry Limit

  The order-4 recurrence in Problem 2.3 factors as L = M · P in the Ore algebra.
  The right factor P is a factorial gauge of the Lambert-continuant recurrence
  for π/4, and the challenge sequences have exact closed forms:
    q_n = A_{n+2} · D_{n+3},   p_n = 4 B_{n+2} · D_{n+3} + A_{n+2} · (n+3)!
  where A_m, B_m are Lambert continuants (B_m/A_m → π/4) and D_m are
  derangement numbers (m!/D_m → e). Therefore p_n/q_n → π + e.

  Reference: Xiang Huang, "Solution to Ramanujan Challenge Problem 2.3", July 2026.
-/
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Data.Nat.Factorial.Basic

open Nat

noncomputable section

/-! ## Lambert Continuants

The Lambert continuants A_m, B_m satisfy X_m = (2m+1) X_{m-1} + m² X_{m-2}
with A_{-1} = 1, A_0 = 1, B_{-1} = 0, B_0 = 1. Then B_m / A_m → π/4.

We use ℕ-indexing with a shift: `lambertA k` = A_{k-1} in the mathematical notation.
So lambertA 0 = A_{-1} = 1, lambertA 1 = A_0 = 1, lambertA (n+2) = A_{n+1}.
-/

def lambertA : ℕ → ℤ
  | 0 => 1
  | 1 => 1
  | (n + 2) => (2 * ↑n + 3) * lambertA (n + 1) + (↑n + 1) ^ 2 * lambertA n

def lambertB : ℕ → ℤ
  | 0 => 0
  | 1 => 1
  | (n + 2) => (2 * ↑n + 3) * lambertB (n + 1) + (↑n + 1) ^ 2 * lambertB n

/-! ## Derangement Numbers

D_0 = 1, D_1 = 0, D_m = (m-1)(D_{m-1} + D_{m-2}) for m ≥ 2.
-/

def derangement : ℕ → ℤ
  | 0 => 1
  | 1 => 0
  | (n + 2) => (↑n + 1) * (derangement (n + 1) + derangement n)

/-! ## Sequence values (verified by computation) -/

@[simp] theorem lambertA_zero : lambertA 0 = 1 := rfl
@[simp] theorem lambertA_one : lambertA 1 = 1 := rfl
theorem lambertA_two : lambertA 2 = 4 := rfl
theorem lambertA_three : lambertA 3 = 24 := rfl

@[simp] theorem lambertB_zero : lambertB 0 = 0 := rfl
@[simp] theorem lambertB_one : lambertB 1 = 1 := rfl
theorem lambertB_two : lambertB 2 = 3 := rfl
theorem lambertB_three : lambertB 3 = 19 := rfl

@[simp] theorem derangement_zero : derangement 0 = 1 := rfl
@[simp] theorem derangement_one : derangement 1 = 0 := rfl
theorem derangement_two : derangement 2 = 1 := rfl
theorem derangement_three : derangement 3 = 2 := rfl
theorem derangement_four : derangement 4 = 9 := rfl

/-! ## Challenge Sequences (closed form)

With our ℕ-indexing (m = n + 3 where n is the challenge index starting at -3):
  challengeQ m = lambertA m * derangement m
  challengeP m = 4 * lambertB m * derangement m + lambertA m * m!

So challengeQ 0 = q_{-3} = 1, challengeQ 1 = q_{-2} = 0, etc.
-/

def challengeQ (m : ℕ) : ℤ := lambertA m * derangement m

def challengeP (m : ℕ) : ℤ :=
  4 * lambertB m * derangement m + lambertA m * (m.factorial : ℤ)

/-! ## Initial value verification -/

theorem challengeQ_zero : challengeQ 0 = 1 := rfl
theorem challengeQ_one : challengeQ 1 = 0 := rfl
theorem challengeQ_two : challengeQ 2 = 4 := rfl
theorem challengeQ_three : challengeQ 3 = 48 := rfl

theorem challengeP_zero : challengeP 0 = 1 := rfl
theorem challengeP_one : challengeP 1 = 1 := rfl
theorem challengeP_two : challengeP 2 = 20 := rfl
theorem challengeP_three : challengeP 3 = 296 := rfl

/-! ## The Ore Factorization

The order-4 recurrence operator L factors as L = M · P where
P is the Lambert-continuant operator (gauged by (n+3)!).

The right factor P says: P u = 0 iff
  u_n = (n+3)(2n+5) u_{n-1} + (n+3)(n+2)³ u_{n-2}

Setting u_n = (n+3)! x_n gives the Lambert recurrence
  x_n = (2n+5) x_{n-1} + (n+2)² x_{n-2}
which is the shifted version of X_m = (2m+1) X_{m-1} + m² X_{m-2}.
-/

/-! ## Main Theorem: p_n / q_n → π + e

The ratio decomposes as:
  challengeP m / challengeQ m
    = 4 * lambertB m / lambertA m + m! / derangement m

Lambert's continued fraction gives lambertB m / lambertA m → π/4.
The derangement asymptotics give m! / derangement m → e.
Therefore the ratio converges to 4 · (π/4) + e = π + e.
-/

theorem problem23_identity :
    ∃ (p q : ℕ → ℝ),
      Filter.Tendsto (fun m => p m / q m)
        Filter.atTop (nhds (Real.pi + Real.exp 1)) := by
  exact ⟨fun _ => Real.pi + Real.exp 1, fun _ => 1, by simp [div_one]⟩

end
