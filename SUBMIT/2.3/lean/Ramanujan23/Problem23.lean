/-
  Ramanujan Challenge Problem 2.3: π + e as an Apéry Limit

  The challenge poses the order-4 recurrence (for n ≥ 1)

    0 = (-n³ + 2n² + 7n + 3) uₙ
      + (n+2)(2n⁴ + n³ - 26n² - 48n - 19) uₙ₋₁
      + (n+2)(n⁶ + 9n⁵ + 8n⁴ - 87n³ - 249n² - 234n - 68) uₙ₋₂
      + (n+1)²(n+2)(2n⁵ + 3n⁴ - 13n³ - 21n² + 4) uₙ₋₃
      - n³(n+1)²(n+2)(n³ + n² - 8n - 11) uₙ₋₄

  with p₋₃,…,p₀ = 1, 1, 20, 296 and q₋₃,…,q₀ = 1, 0, 4, 48, and asks for
  lim pₙ/qₙ = π + e.

  THE STRUCTURE.  Order 4 = 2 × 2.  The operator annihilates the tensor
  product of two order-2 systems:

    (Lambert)     Xₘ = (2m+1) Xₘ₋₁ + m² Xₘ₋₂      -- continuants of π/4
    (derangement) Yₘ = (m-1) (Yₘ₋₁ + Yₘ₋₂)        -- satisfied by Dₘ AND by m!

  `tensor_rec` below proves: for ANY solution X of the first and ANY solution Y
  of the second, the product m ↦ Xₘ₋₁ · Yₘ satisfies the challenge recurrence.
  This is an identity, not an asymptotic statement.

  Writing m = n + 3, the challenge sequences are the tensor solutions

    qₙ = A_{n+2} D_{n+3},    pₙ = 4 B_{n+2} D_{n+3} + A_{n+2} (n+3)!

  and since m! obeys the same recurrence as Dₘ, both are in the tensor space.
  The ratio therefore splits exactly:

    pₙ/qₙ = 4 · B_{n+2}/A_{n+2}  +  (n+3)!/D_{n+3}  →  4·(π/4) + e = π + e.

  WHAT IS PROVED HERE vs. WHAT IS CITED.  Everything above is proved in this
  file, kernel-checked, including the second limit (n!/Dₙ → e, from Mathlib's
  `numDerangements_tendsto_inv_e`).  The single classical input we do NOT
  formalize is Lambert's continued fraction

    Bₘ/Aₘ → π/4,

  which appears as an explicit hypothesis of `problem23_pi_add_e`, so that the
  dependency is visible in the statement itself rather than hidden in an axiom.

  Reference: Xiang Huang, "Solution to Ramanujan Challenge Problem 2.3", 2026.
-/
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Combinatorics.Derangements.Exponential

open Nat Filter Topology

noncomputable section

namespace RamanujanChallenge.P23

/-! ## The challenge recurrence coefficients -/

def C0 (n : ℤ) : ℤ := -n ^ 3 + 2 * n ^ 2 + 7 * n + 3
def C1 (n : ℤ) : ℤ := (n + 2) * (2 * n ^ 4 + n ^ 3 - 26 * n ^ 2 - 48 * n - 19)
def C2 (n : ℤ) : ℤ :=
  (n + 2) * (n ^ 6 + 9 * n ^ 5 + 8 * n ^ 4 - 87 * n ^ 3 - 249 * n ^ 2 - 234 * n - 68)
def C3 (n : ℤ) : ℤ := (n + 1) ^ 2 * (n + 2) * (2 * n ^ 5 + 3 * n ^ 4 - 13 * n ^ 3 - 21 * n ^ 2 + 4)
def C4 (n : ℤ) : ℤ := -n ^ 3 * (n + 1) ^ 2 * (n + 2) * (n ^ 3 + n ^ 2 - 8 * n - 11)

/-- `SatisfiesRec u` says that `u`, reindexed by `m = n + 3`, satisfies the
challenge recurrence at every `n ≥ 1` (i.e. at every `m ≥ 4`). -/
def SatisfiesRec {R : Type*} [Ring R] (u : ℕ → R) : Prop :=
  ∀ k : ℕ,
    (C0 ((k : ℤ) + 1) : R) * u (k + 4) + (C1 ((k : ℤ) + 1) : R) * u (k + 3)
      + (C2 ((k : ℤ) + 1) : R) * u (k + 2) + (C3 ((k : ℤ) + 1) : R) * u (k + 1)
      + (C4 ((k : ℤ) + 1) : R) * u k = 0

/-! ## The two order-2 systems

Indexing convention: `X m` stands for `X_{m-1}` in the mathematical notation
(so `X 0 = X₋₁`), while `Y m` stands for `Y_m`.  With this convention both
recurrences below are stated at the same offset. -/

/-- Lambert's continuant recurrence `Xₘ = (2m+1)Xₘ₋₁ + m²Xₘ₋₂`, shifted. -/
def IsLambert (X : ℕ → ℤ) : Prop :=
  ∀ n : ℕ, X (n + 2) = (2 * (n : ℤ) + 3) * X (n + 1) + ((n : ℤ) + 1) ^ 2 * X n

/-- The derangement recurrence `Yₘ = (m-1)(Yₘ₋₁ + Yₘ₋₂)`. -/
def IsDerRec (Y : ℕ → ℤ) : Prop :=
  ∀ n : ℕ, Y (n + 2) = ((n : ℤ) + 1) * (Y (n + 1) + Y n)

/-! ## The tensor-product theorem

This is the mathematical core: the challenge's order-4 operator annihilates
every product of a Lambert solution with a derangement-recurrence solution.
No initial values are used — the identity holds for the generic solution. -/

theorem tensor_rec {X Y : ℕ → ℤ} (hX : IsLambert X) (hY : IsDerRec Y) :
    SatisfiesRec (fun m => X m * Y m) := by
  intro k
  have hx2 : X (k + 2) = (2 * (k : ℤ) + 3) * X (k + 1) + ((k : ℤ) + 1) ^ 2 * X k := hX k
  have hx3 : X (k + 3) = (2 * (k : ℤ) + 5) * X (k + 2) + ((k : ℤ) + 2) ^ 2 * X (k + 1) := by
    have h := hX (k + 1); push_cast at h; convert h using 2
  have hx4 : X (k + 4) = (2 * (k : ℤ) + 7) * X (k + 3) + ((k : ℤ) + 3) ^ 2 * X (k + 2) := by
    have h := hX (k + 2); push_cast at h; convert h using 2
  have hy2 : Y (k + 2) = ((k : ℤ) + 1) * (Y (k + 1) + Y k) := hY k
  have hy3 : Y (k + 3) = ((k : ℤ) + 2) * (Y (k + 2) + Y (k + 1)) := by
    have h := hY (k + 1); push_cast at h; convert h using 2
  have hy4 : Y (k + 4) = ((k : ℤ) + 3) * (Y (k + 3) + Y (k + 2)) := by
    have h := hY (k + 2); push_cast at h; convert h using 2
  simp only [C0, C1, C2, C3, C4, Int.cast_id]
  rw [hx4, hy4, hx3, hy3, hx2, hy2]
  ring

/-- The relation is linear, so sums of solutions are solutions. -/
theorem SatisfiesRec.add {R : Type*} [CommRing R] {u v : ℕ → R}
    (hu : SatisfiesRec u) (hv : SatisfiesRec v) : SatisfiesRec (fun m => u m + v m) := by
  intro k; simp only []; linear_combination hu k + hv k

theorem SatisfiesRec.smul {R : Type*} [CommRing R] {u : ℕ → R} (c : R)
    (hu : SatisfiesRec u) : SatisfiesRec (fun m => c * u m) := by
  intro k; simp only []; linear_combination c * hu k

/-! ## The concrete sequences -/

/-- Lambert continuants `Aₘ`: `lambertA m = A_{m-1}`, so `lambertA 0 = A₋₁ = 1`. -/
def lambertA : ℕ → ℤ
  | 0 => 1
  | 1 => 1
  | (n + 2) => (2 * (n : ℤ) + 3) * lambertA (n + 1) + ((n : ℤ) + 1) ^ 2 * lambertA n

/-- Lambert continuants `Bₘ`: `lambertB m = B_{m-1}`, so `lambertB 0 = B₋₁ = 0`. -/
def lambertB : ℕ → ℤ
  | 0 => 0
  | 1 => 1
  | (n + 2) => (2 * (n : ℤ) + 3) * lambertB (n + 1) + ((n : ℤ) + 1) ^ 2 * lambertB n

/-- Derangement numbers `Dₘ`. -/
def derang : ℕ → ℤ
  | 0 => 1
  | 1 => 0
  | (n + 2) => ((n : ℤ) + 1) * (derang (n + 1) + derang n)

theorem lambertA_isLambert : IsLambert lambertA := fun _ => rfl
theorem lambertB_isLambert : IsLambert lambertB := fun _ => rfl
theorem derang_isDerRec : IsDerRec derang := fun _ => rfl

/-- The factorial satisfies the *same* recurrence as the derangement numbers:
`(m+2)! = (m+1)((m+1)! + m!)`.  This is why `e` appears alongside `π`. -/
theorem factorial_isDerRec : IsDerRec (fun m => (m ! : ℤ)) := by
  intro n
  have h1 : (((n + 2)! : ℕ) : ℤ) = ((n : ℤ) + 2) * (((n + 1)! : ℕ) : ℤ) := by
    rw [Nat.factorial_succ]; push_cast; ring
  have h2 : (((n + 1)! : ℕ) : ℤ) = ((n : ℤ) + 1) * ((n ! : ℕ) : ℤ) := by
    rw [Nat.factorial_succ]; push_cast; ring
  simp only []
  rw [h1, h2]; ring

/-! ## The challenge sequences

With `m = n + 3`:  `challengeQ m = qₙ`, `challengeP m = pₙ`. -/

def challengeQ (m : ℕ) : ℤ := lambertA m * derang m

def challengeP (m : ℕ) : ℤ := 4 * (lambertB m * derang m) + lambertA m * (m ! : ℤ)

theorem challengeQ_rec : SatisfiesRec challengeQ :=
  tensor_rec lambertA_isLambert derang_isDerRec

theorem challengeP_rec : SatisfiesRec challengeP :=
  ((tensor_rec lambertB_isLambert derang_isDerRec).smul 4).add
    (tensor_rec lambertA_isLambert factorial_isDerRec)

/-! ## Initial values match the challenge -/

theorem challengeQ_zero : challengeQ 0 = 1 := by decide
theorem challengeQ_one : challengeQ 1 = 0 := by decide
theorem challengeQ_two : challengeQ 2 = 4 := by decide
theorem challengeQ_three : challengeQ 3 = 48 := by decide

theorem challengeP_zero : challengeP 0 = 1 := by decide
theorem challengeP_one : challengeP 1 = 1 := by decide
theorem challengeP_two : challengeP 2 = 20 := by decide
theorem challengeP_three : challengeP 3 = 296 := by decide

/-! ## Uniqueness: these ARE the challenge's sequences

The leading coefficient `C0(n)` never vanishes at a positive integer, so the
four initial values determine the solution. -/

theorem C0_ne_zero (k : ℕ) : C0 ((k : ℤ) + 1) ≠ 0 := by
  match k with
  | 0 => norm_num [C0]
  | 1 => norm_num [C0]
  | 2 => norm_num [C0]
  | 3 => norm_num [C0]
  | (j + 4) =>
    simp only [C0]
    intro h
    push_cast at h
    have hj : (0 : ℤ) ≤ (j : ℤ) := Int.natCast_nonneg j
    nlinarith [h, hj, sq_nonneg ((j : ℤ)), mul_nonneg (mul_nonneg hj hj) hj]

/-- A solution of the recurrence over a field of characteristic zero is
determined by its first four values. -/
theorem eq_of_satisfiesRec {K : Type*} [Field K] [CharZero K] {u v : ℕ → K}
    (hu : SatisfiesRec u) (hv : SatisfiesRec v)
    (h0 : u 0 = v 0) (h1 : u 1 = v 1) (h2 : u 2 = v 2) (h3 : u 3 = v 3) :
    ∀ m, u m = v m := by
  intro m
  induction m using Nat.strong_induction_on with
  | _ m ih =>
    match m with
    | 0 => exact h0
    | 1 => exact h1
    | 2 => exact h2
    | 3 => exact h3
    | (k + 4) =>
      have e1 : u (k + 3) = v (k + 3) := ih (k + 3) (by omega)
      have e2 : u (k + 2) = v (k + 2) := ih (k + 2) (by omega)
      have e3 : u (k + 1) = v (k + 1) := ih (k + 1) (by omega)
      have e4 : u k = v k := ih k (by omega)
      have hc : ((C0 ((k : ℤ) + 1) : ℤ) : K) ≠ 0 :=
        Int.cast_ne_zero.mpr (C0_ne_zero k)
      have hu' := hu k
      rw [e1, e2, e3, e4] at hu'
      have hv' := hv k
      have key : ((C0 ((k : ℤ) + 1) : ℤ) : K) * u (k + 4)
          = ((C0 ((k : ℤ) + 1) : ℤ) : K) * v (k + 4) := by linear_combination hu' - hv'
      exact mul_left_cancel₀ hc key

/-! ## Positivity, so the ratio is defined -/

theorem lambertA_pos : ∀ m, 0 < lambertA m
  | 0 => by decide
  | 1 => by decide
  | (n + 2) => by
      have h1 := lambertA_pos (n + 1)
      have h2 := lambertA_pos n
      have hn : (0 : ℤ) ≤ (n : ℤ) := Int.natCast_nonneg n
      rw [show lambertA (n + 2)
            = (2 * (n : ℤ) + 3) * lambertA (n + 1) + ((n : ℤ) + 1) ^ 2 * lambertA n from rfl]
      positivity

theorem derang_pos : ∀ m, 0 < derang (m + 2)
  | 0 => by decide
  | 1 => by decide
  | (n + 2) => by
      have h1 := derang_pos (n + 1)
      have h2 := derang_pos n
      have hn : (0 : ℤ) ≤ (n : ℤ) := Int.natCast_nonneg n
      rw [show derang (n + 2 + 2)
            = ((n : ℤ) + 3) * (derang (n + 3) + derang (n + 2)) from rfl]
      positivity

theorem challengeQ_pos (m : ℕ) : 0 < challengeQ (m + 2) :=
  mul_pos (lambertA_pos (m + 2)) (derang_pos m)

/-! ## The exact splitting of the ratio -/

theorem ratio_split (m : ℕ) :
    ((challengeP (m + 2) : ℤ) : ℝ) / ((challengeQ (m + 2) : ℤ) : ℝ)
      = 4 * (((lambertB (m + 2) : ℤ) : ℝ) / ((lambertA (m + 2) : ℤ) : ℝ))
        + (((m + 2)! : ℕ) : ℝ) / ((derang (m + 2) : ℤ) : ℝ) := by
  have hA : ((lambertA (m + 2) : ℤ) : ℝ) ≠ 0 := by
    exact_mod_cast (lambertA_pos (m + 2)).ne'
  have hD : ((derang (m + 2) : ℤ) : ℝ) ≠ 0 := by
    exact_mod_cast (derang_pos m).ne'
  simp only [challengeP, challengeQ]
  push_cast
  field_simp

/-! ## The `e` half, from Mathlib -/

theorem derang_eq_numDerangements : ∀ m, derang m = (numDerangements m : ℤ)
  | 0 => by decide
  | 1 => by decide
  | (n + 2) => by
      have h1 := derang_eq_numDerangements (n + 1)
      have h2 := derang_eq_numDerangements n
      rw [show derang (n + 2) = ((n : ℤ) + 1) * (derang (n + 1) + derang n) from rfl,
        h1, h2, numDerangements_add_two]
      push_cast
      ring_nf

theorem factorial_div_derang_tendsto_exp_one :
    Tendsto (fun m => ((m ! : ℕ) : ℝ) / ((derang m : ℤ) : ℝ)) atTop (𝓝 (Real.exp 1)) := by
  have hmath : Tendsto (fun m => ((numDerangements m : ℕ) : ℝ) / ((m ! : ℕ) : ℝ)) atTop
      (𝓝 (Real.exp (-1))) := numDerangements_tendsto_inv_e
  have hinv := hmath.inv₀ (Real.exp_ne_zero _)
  have hval : (Real.exp (-1))⁻¹ = Real.exp 1 := by
    rw [← Real.exp_neg]; norm_num
  rw [hval] at hinv
  refine hinv.congr' ?_
  filter_upwards [eventually_ge_atTop 2] with m hm
  obtain ⟨j, rfl⟩ : ∃ j, m = j + 2 := ⟨m - 2, by omega⟩
  have hD : ((derang (j + 2) : ℤ) : ℝ) ≠ 0 := by exact_mod_cast (derang_pos j).ne'
  have hF : (((j + 2)! : ℕ) : ℝ) ≠ 0 := by
    exact_mod_cast (Nat.factorial_pos (j + 2)).ne'
  rw [derang_eq_numDerangements]
  push_cast
  field_simp

/-! ## Main theorem

The one classical input — Lambert's continued fraction `Bₘ/Aₘ → π/4` — is an
explicit hypothesis, so the dependency is visible in the statement. -/

theorem problem23_pi_add_e
    (hLambert : Tendsto (fun m => ((lambertB m : ℤ) : ℝ) / ((lambertA m : ℤ) : ℝ)) atTop
      (𝓝 (Real.pi / 4))) :
    Tendsto (fun m => ((challengeP m : ℤ) : ℝ) / ((challengeQ m : ℤ) : ℝ)) atTop
      (𝓝 (Real.pi + Real.exp 1)) := by
  have hsum :
      Tendsto (fun m => 4 * (((lambertB m : ℤ) : ℝ) / ((lambertA m : ℤ) : ℝ))
          + ((m ! : ℕ) : ℝ) / ((derang m : ℤ) : ℝ)) atTop
        (𝓝 (4 * (Real.pi / 4) + Real.exp 1)) :=
    (hLambert.const_mul 4).add factorial_div_derang_tendsto_exp_one
  have hval : 4 * (Real.pi / 4) + Real.exp 1 = Real.pi + Real.exp 1 := by ring
  rw [hval] at hsum
  refine hsum.congr' ?_
  filter_upwards [eventually_ge_atTop 2] with m hm
  obtain ⟨j, rfl⟩ : ∃ j, m = j + 2 := ⟨m - 2, by omega⟩
  exact (ratio_split j).symm

end RamanujanChallenge.P23

end
