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

  EVERYTHING IS PROVED HERE.  The Lambert limit is obtained from a positive
  moment representation for its alternating remainder; the moment recurrence,
  its initial values, and a geometric bound are all proved below.  The second
  limit n!/Dₙ → e follows from Mathlib's
  `numDerangements_tendsto_inv_e`.  Thus `problem23_pi_add_e` has no mathematical
  hypotheses.

  Reference: Xiang Huang, "Solution to Ramanujan Challenge Problem 2.3", 2026.
-/
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic
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

/-! ## The Lambert half

The classical Lambert evaluation can be proved directly from a positive moment
representation.  Put `s = √2`, `c = 2 - s`, `P(x) = x(1-x)`, and
`D(x) = 1 - cP(x)`.  The moments

`Kₙ = ∫₀¹ P(x)ⁿ / D(x)ⁿ⁺¹ dx`

satisfy

`2n Kₙ = (n-1)Kₙ₋₂ - (2n-1)√2 Kₙ₋₁`.

Consequently `n! (√2)ⁿ Kₙ / √2` is exactly the alternating Lambert
remainder.  Positivity of the kernel and
`P(x)/D(x) ≤ 1/(2+√2)` give a geometric error bound. -/

private def lambertP (x : ℝ) : ℝ := x * (1 - x)

private def lambertD (x : ℝ) : ℝ :=
  1 - (2 - Real.sqrt 2) * lambertP x

private def lambertMoment (n : ℕ) : ℝ :=
  ∫ x in (0 : ℝ)..1, lambertP x ^ n / lambertD x ^ (n + 1)

private lemma sqrt_two_pos : 0 < Real.sqrt 2 := Real.sqrt_pos.2 (by norm_num)

private lemma sqrt_two_sq : (Real.sqrt 2) ^ 2 = 2 := by norm_num

private lemma sqrt_two_lt_two : Real.sqrt 2 < 2 := by nlinarith [sqrt_two_sq, sqrt_two_pos]

private lemma one_lt_sqrt_two : 1 < Real.sqrt 2 := by nlinarith [sqrt_two_sq, sqrt_two_pos]

private lemma lambertP_nonneg {x : ℝ} (hx : x ∈ Set.Icc (0 : ℝ) 1) :
    0 ≤ lambertP x := by
  exact mul_nonneg hx.1 (sub_nonneg.mpr hx.2)

private lemma four_mul_lambertP_le_one (x : ℝ) : 4 * lambertP x ≤ 1 := by
  dsimp [lambertP]
  nlinarith [sq_nonneg (2 * x - 1)]

private lemma lambertD_pos (x : ℝ) : 0 < lambertD x := by
  have hp : lambertP x ≤ 1 / 4 := by
    nlinarith [four_mul_lambertP_le_one x]
  have hc : 0 ≤ 2 - Real.sqrt 2 := sub_nonneg.mpr sqrt_two_lt_two.le
  have hmul :
      (2 - Real.sqrt 2) * lambertP x ≤ (2 - Real.sqrt 2) * (1 / 4) :=
    mul_le_mul_of_nonneg_left hp hc
  dsimp [lambertD]
  nlinarith [sqrt_two_pos]

private lemma lambertMoment_integrand_continuousOn (n : ℕ) :
    ContinuousOn (fun x : ℝ => lambertP x ^ n / lambertD x ^ (n + 1))
      (Set.Icc (0 : ℝ) 1) := by
  apply ContinuousOn.div
  · unfold lambertP
    fun_prop
  · unfold lambertD lambertP
    fun_prop
  · intro x _
    exact pow_ne_zero _ (lambertD_pos x).ne'

private lemma continuous_lambertP : Continuous lambertP := by
  unfold lambertP
  fun_prop

private lemma continuous_lambertD : Continuous lambertD := by
  simpa only [lambertD] using
    continuous_const.sub (continuous_lambertP.const_mul (2 - Real.sqrt 2))

private lemma arctan_sqrt_two_sub_one :
    Real.arctan (Real.sqrt 2 - 1) = Real.pi / 8 := by
  have h₁ : -1 < Real.sqrt 2 - 1 := by nlinarith [sqrt_two_pos]
  have h₂ : Real.sqrt 2 - 1 < 1 := by nlinarith [sqrt_two_lt_two]
  have h := Real.two_mul_arctan h₁ h₂
  have hden : 1 - (Real.sqrt 2 - 1) ^ 2 ≠ 0 := by
    nlinarith [sqrt_two_sq, sqrt_two_pos]
  have hfrac :
      2 * (Real.sqrt 2 - 1) / (1 - (Real.sqrt 2 - 1) ^ 2) = 1 := by
    field_simp
    nlinarith [sqrt_two_sq]
  rw [hfrac, Real.arctan_one] at h
  linarith

private lemma lambert_antiderivative_hasDerivAt (x : ℝ) :
    HasDerivAt
      (fun y : ℝ =>
        (2 / Real.sqrt 2) * Real.arctan ((Real.sqrt 2 - 1) * (2 * y - 1)))
      (1 / lambertD x) x := by
  have hs0 : Real.sqrt 2 ≠ 0 := sqrt_two_pos.ne'
  have hD : lambertD x ≠ 0 := (lambertD_pos x).ne'
  have harg :
      HasDerivAt (fun y : ℝ => (Real.sqrt 2 - 1) * (2 * y - 1))
        (2 * (Real.sqrt 2 - 1)) x := by
    convert ((((hasDerivAt_id x).const_mul 2).sub_const 1).const_mul
      (Real.sqrt 2 - 1)) using 1
    all_goals ring
  have hcalc :
      (2 / Real.sqrt 2) *
          (1 / (1 + ((Real.sqrt 2 - 1) * (2 * x - 1)) ^ 2) *
            (2 * (Real.sqrt 2 - 1))) =
        1 / lambertD x := by
    let Q : ℝ := 1 + ((Real.sqrt 2 - 1) * (2 * x - 1)) ^ 2
    have hQ : Q ≠ 0 := by
      dsimp [Q]
      positivity
    change (2 / Real.sqrt 2) * (1 / Q * (2 * (Real.sqrt 2 - 1))) =
      1 / lambertD x
    rw [show (2 / Real.sqrt 2) * (1 / Q * (2 * (Real.sqrt 2 - 1))) =
        (4 * (Real.sqrt 2 - 1) / Real.sqrt 2) / Q by ring]
    apply (eq_div_iff hD).2
    rw [div_mul_eq_mul_div]
    apply (div_eq_iff hQ).2
    rw [div_mul_eq_mul_div]
    apply (div_eq_iff hs0).2
    have hs3 : (Real.sqrt 2) ^ 3 = 2 * Real.sqrt 2 := by
      rw [pow_succ, sqrt_two_sq]
    dsimp [lambertD, lambertP, Q]
    ring_nf
    rw [hs3, sqrt_two_sq]
    ring
  convert harg.arctan.const_mul (2 / Real.sqrt 2) using 1
  exact hcalc.symm

private lemma lambertMoment_zero :
    lambertMoment 0 = Real.sqrt 2 * (Real.pi / 4) := by
  let F : ℝ → ℝ := fun x =>
    (2 / Real.sqrt 2) * Real.arctan ((Real.sqrt 2 - 1) * (2 * x - 1))
  have hFTC :
      (∫ x in (0 : ℝ)..1, 1 / lambertD x) = F 1 - F 0 := by
    apply intervalIntegral.integral_eq_sub_of_hasDerivAt
    · intro x _
      exact lambert_antiderivative_hasDerivAt x
    · apply ContinuousOn.intervalIntegrable
      apply ContinuousOn.div
      · fun_prop
      · unfold lambertD lambertP
        fun_prop
      · intro x _
        exact (lambertD_pos x).ne'
  rw [show lambertMoment 0 = ∫ x in (0 : ℝ)..1, 1 / lambertD x by
    simp [lambertMoment]]
  rw [hFTC]
  dsimp [F]
  have hpos :
      (Real.sqrt 2 - 1) * (2 * (1 : ℝ) - 1) = Real.sqrt 2 - 1 := by ring
  have hneg :
      (Real.sqrt 2 - 1) * (2 * (0 : ℝ) - 1) = -(Real.sqrt 2 - 1) := by ring
  rw [hpos, hneg, Real.arctan_neg, arctan_sqrt_two_sub_one]
  have hs0 : Real.sqrt 2 ≠ 0 := sqrt_two_pos.ne'
  field_simp
  nlinarith [sqrt_two_sq]

private lemma lambert_first_aux_hasDerivAt (x : ℝ) :
    HasDerivAt (fun y : ℝ => (1 - 2 * y) / lambertD y)
      (-Real.sqrt 2 * (1 / lambertD x) -
        2 * (lambertP x / lambertD x ^ 2)) x := by
  have hP : HasDerivAt lambertP (1 - 2 * x) x := by
    unfold lambertP
    have h :=
      (hasDerivAt_id x).mul ((hasDerivAt_const x 1).sub (hasDerivAt_id x))
    convert h using 1
    simp only [Pi.sub_apply, id_eq]
    ring
  have hD :
      HasDerivAt lambertD (-(2 - Real.sqrt 2) * (1 - 2 * x)) x := by
    unfold lambertD
    have h :=
      (hasDerivAt_const x 1).sub (hP.const_mul (2 - Real.sqrt 2))
    convert h using 1
    ring
  have hN :
      HasDerivAt (fun y : ℝ => 1 - 2 * y) (-2) x := by
    have h :=
      ((hasDerivAt_const x 1).sub (hasDerivAt_id x)).sub (hasDerivAt_id x)
    convert h using 1
    · funext y
      change 1 - 2 * y = (1 - y) - y
      ring
    · norm_num
  have hD0 : lambertD x ≠ 0 := (lambertD_pos x).ne'
  have hcalc :
      ((-2) * lambertD x -
          (1 - 2 * x) * (-(2 - Real.sqrt 2) * (1 - 2 * x))) /
          lambertD x ^ 2 =
        -Real.sqrt 2 * (1 / lambertD x) -
          2 * (lambertP x / lambertD x ^ 2) := by
    field_simp [hD0]
    dsimp [lambertD, lambertP]
    ring_nf
    rw [sqrt_two_sq]
    ring
  convert hN.div hD hD0 using 1
  exact hcalc.symm

private lemma lambertMoment_one :
    lambertMoment 1 = 1 - Real.pi / 4 := by
  let f : ℝ → ℝ := fun x => (1 - 2 * x) / lambertD x
  let g : ℝ → ℝ := fun x =>
    -Real.sqrt 2 * (1 / lambertD x) -
      2 * (lambertP x / lambertD x ^ 2)
  have hgcont : Continuous g := by
    dsimp [g]
    have hc1 : Continuous (fun x : ℝ => 1 / lambertD x) :=
      continuous_const.div continuous_lambertD fun x => (lambertD_pos x).ne'
    have hc2 : Continuous (fun x : ℝ => lambertP x / lambertD x ^ 2) :=
      continuous_lambertP.div (continuous_lambertD.pow 2) fun x =>
        pow_ne_zero _ (lambertD_pos x).ne'
    fun_prop
  have hFTC : (∫ x in (0 : ℝ)..1, g x) = -2 := by
    calc
      (∫ x in (0 : ℝ)..1, g x) = f 1 - f 0 := by
        apply intervalIntegral.integral_eq_sub_of_hasDerivAt
        · intro x _
          exact lambert_first_aux_hasDerivAt x
        · exact hgcont.continuousOn.intervalIntegrable
      _ = -2 := by
        dsimp [f, lambertD, lambertP]
        norm_num
  have hK :
      (∫ x in (0 : ℝ)..1, g x) =
        -Real.sqrt 2 * lambertMoment 0 - 2 * lambertMoment 1 := by
    have h0 :
        IntervalIntegrable (fun x : ℝ => 1 / lambertD x) MeasureTheory.volume 0 1 :=
      (continuous_const.div continuous_lambertD fun x =>
        (lambertD_pos x).ne').continuousOn.intervalIntegrable
    have h1 :
        IntervalIntegrable (fun x : ℝ => lambertP x / lambertD x ^ 2)
          MeasureTheory.volume 0 1 :=
      (continuous_lambertP.div (continuous_lambertD.pow 2) fun x =>
        pow_ne_zero _ (lambertD_pos x).ne').continuousOn.intervalIntegrable
    dsimp [g]
    rw [intervalIntegral.integral_sub]
    · rw [intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul]
      simp only [lambertMoment]
      norm_num
    · exact h0.const_mul _
    · exact h1.const_mul _
  rw [hFTC] at hK
  rw [lambertMoment_zero] at hK
  have hsprod :
      Real.sqrt 2 * (Real.sqrt 2 * (Real.pi / 4)) =
        2 * (Real.pi / 4) := by
    calc
      Real.sqrt 2 * (Real.sqrt 2 * (Real.pi / 4)) =
          (Real.sqrt 2) ^ 2 * (Real.pi / 4) := by ring
      _ = 2 * (Real.pi / 4) := by rw [sqrt_two_sq]
  have hsprod_neg :
      -Real.sqrt 2 * (Real.sqrt 2 * (Real.pi / 4)) =
        -2 * (Real.pi / 4) := by
    calc
      -Real.sqrt 2 * (Real.sqrt 2 * (Real.pi / 4)) =
          -(Real.sqrt 2 * (Real.sqrt 2 * (Real.pi / 4))) := by ring
      _ = -(2 * (Real.pi / 4)) := by rw [hsprod]
      _ = -2 * (Real.pi / 4) := by ring
  rw [hsprod_neg] at hK
  linarith

private lemma continuous_lambertMoment_integrand (n : ℕ) :
    Continuous (fun x : ℝ => lambertP x ^ n / lambertD x ^ (n + 1)) :=
  (continuous_lambertP.pow n).div (continuous_lambertD.pow (n + 1)) fun x =>
    pow_ne_zero _ (lambertD_pos x).ne'

private lemma lambert_recurrence_aux_hasDerivAt (n : ℕ) (x : ℝ) :
    HasDerivAt
      (fun y : ℝ =>
        (1 - 2 * y) * lambertP y ^ (n + 1) / lambertD y ^ (n + 2))
      (((n + 1 : ℕ) : ℝ) *
          (lambertP x ^ n / lambertD x ^ (n + 1)) -
        (2 * (n : ℝ) + 3) * Real.sqrt 2 *
          (lambertP x ^ (n + 1) / lambertD x ^ (n + 2)) -
        2 * ((n + 2 : ℕ) : ℝ) *
          (lambertP x ^ (n + 2) / lambertD x ^ (n + 3))) x := by
  have hP : HasDerivAt lambertP (1 - 2 * x) x := by
    unfold lambertP
    have h :=
      (hasDerivAt_id x).mul ((hasDerivAt_const x 1).sub (hasDerivAt_id x))
    convert h using 1
    simp only [Pi.sub_apply, id_eq]
    ring
  have hD :
      HasDerivAt lambertD (-(2 - Real.sqrt 2) * (1 - 2 * x)) x := by
    unfold lambertD
    have h :=
      (hasDerivAt_const x 1).sub (hP.const_mul (2 - Real.sqrt 2))
    convert h using 1
    ring
  have hN :
      HasDerivAt (fun y : ℝ => 1 - 2 * y) (-2) x := by
    have h :=
      ((hasDerivAt_const x 1).sub (hasDerivAt_id x)).sub (hasDerivAt_id x)
    convert h using 1
    · funext y
      change 1 - 2 * y = (1 - y) - y
      ring
    · norm_num
  have hU :
      HasDerivAt
        (fun y : ℝ => (1 - 2 * y) * lambertP y ^ (n + 1))
        ((-2) * lambertP x ^ (n + 1) +
          (1 - 2 * x) * (((n + 1 : ℕ) : ℝ) * lambertP x ^ n *
            (1 - 2 * x))) x := by
    convert hN.mul (hP.pow (n + 1)) using 1
  have hV :
      HasDerivAt (fun y : ℝ => lambertD y ^ (n + 2))
        (((n + 2 : ℕ) : ℝ) * lambertD x ^ (n + 1) *
          (-(2 - Real.sqrt 2) * (1 - 2 * x))) x := by
    convert hD.pow (n + 2) using 1
  have hD0 : lambertD x ≠ 0 := (lambertD_pos x).ne'
  have hV0 : lambertD x ^ (n + 2) ≠ 0 := pow_ne_zero _ hD0
  have hcalc :
      (((-2) * lambertP x ^ (n + 1) +
            (1 - 2 * x) * (((n + 1 : ℕ) : ℝ) * lambertP x ^ n *
              (1 - 2 * x))) *
          lambertD x ^ (n + 2) -
        ((1 - 2 * x) * lambertP x ^ (n + 1)) *
          (((n + 2 : ℕ) : ℝ) * lambertD x ^ (n + 1) *
            (-(2 - Real.sqrt 2) * (1 - 2 * x)))) /
          (lambertD x ^ (n + 2)) ^ 2 =
        ((n + 1 : ℕ) : ℝ) *
            (lambertP x ^ n / lambertD x ^ (n + 1)) -
          (2 * (n : ℝ) + 3) * Real.sqrt 2 *
            (lambertP x ^ (n + 1) / lambertD x ^ (n + 2)) -
          2 * ((n + 2 : ℕ) : ℝ) *
            (lambertP x ^ (n + 2) / lambertD x ^ (n + 3)) := by
    have hp1 : lambertP x ^ (n + 1) = lambertP x ^ n * lambertP x := by
      rw [pow_succ]
    have hp2 : lambertP x ^ (n + 2) = lambertP x ^ n * lambertP x ^ 2 := by
      rw [show n + 2 = (n + 1) + 1 by omega, pow_succ, hp1]
      ring
    have hd1 : lambertD x ^ (n + 1) = lambertD x ^ n * lambertD x := by
      rw [pow_succ]
    have hd2 : lambertD x ^ (n + 2) = lambertD x ^ n * lambertD x ^ 2 := by
      rw [show n + 2 = (n + 1) + 1 by omega, pow_succ, hd1]
      ring
    have hd3 : lambertD x ^ (n + 3) = lambertD x ^ n * lambertD x ^ 3 := by
      rw [show n + 3 = (n + 2) + 1 by omega, pow_succ, hd2]
      ring
    have hn2 :
        (1 - 2 * x) ^ 2 = 1 - 4 * lambertP x := by
      dsimp [lambertP]
      ring
    have hbase :
        ((n + 1 : ℕ) : ℝ) * lambertD x * (1 - 2 * x) ^ 2 -
            2 * lambertD x * lambertP x +
            ((n + 2 : ℕ) : ℝ) * (2 - Real.sqrt 2) *
              (1 - 2 * x) ^ 2 * lambertP x =
          ((n + 1 : ℕ) : ℝ) * lambertD x ^ 2 -
            (2 * (n : ℝ) + 3) * Real.sqrt 2 * lambertP x * lambertD x -
            2 * ((n + 2 : ℕ) : ℝ) * lambertP x ^ 2 := by
      rw [hn2]
      dsimp [lambertD]
      push_cast
      ring_nf
      rw [sqrt_two_sq]
      ring
    have hleft :
        (((-2) * lambertP x ^ (n + 1) +
              (1 - 2 * x) * (((n + 1 : ℕ) : ℝ) * lambertP x ^ n *
                (1 - 2 * x))) *
            lambertD x ^ (n + 2) -
          ((1 - 2 * x) * lambertP x ^ (n + 1)) *
            (((n + 2 : ℕ) : ℝ) * lambertD x ^ (n + 1) *
              (-(2 - Real.sqrt 2) * (1 - 2 * x)))) /
            (lambertD x ^ (n + 2)) ^ 2 =
          (lambertP x ^ n / lambertD x ^ (n + 3)) *
            (((n + 1 : ℕ) : ℝ) * lambertD x * (1 - 2 * x) ^ 2 -
              2 * lambertD x * lambertP x +
              ((n + 2 : ℕ) : ℝ) * (2 - Real.sqrt 2) *
                (1 - 2 * x) ^ 2 * lambertP x) := by
      field_simp [hD0]
      rw [hp1, hd1, hd2, hd3]
      ring
    have hright :
        ((n + 1 : ℕ) : ℝ) *
              (lambertP x ^ n / lambertD x ^ (n + 1)) -
            (2 * (n : ℝ) + 3) * Real.sqrt 2 *
              (lambertP x ^ (n + 1) / lambertD x ^ (n + 2)) -
            2 * ((n + 2 : ℕ) : ℝ) *
              (lambertP x ^ (n + 2) / lambertD x ^ (n + 3)) =
          (lambertP x ^ n / lambertD x ^ (n + 3)) *
            (((n + 1 : ℕ) : ℝ) * lambertD x ^ 2 -
              (2 * (n : ℝ) + 3) * Real.sqrt 2 *
                lambertP x * lambertD x -
              2 * ((n + 2 : ℕ) : ℝ) * lambertP x ^ 2) := by
      field_simp [hD0]
      rw [hp1, hp2, hd1, hd2, hd3]
      ring
    rw [hleft, hright, hbase]
  convert hU.div hV hV0 using 1
  exact hcalc.symm

private lemma lambertMoment_recurrence (n : ℕ) :
    2 * ((n + 2 : ℕ) : ℝ) * lambertMoment (n + 2) =
      ((n + 1 : ℕ) : ℝ) * lambertMoment n -
        (2 * (n : ℝ) + 3) * Real.sqrt 2 * lambertMoment (n + 1) := by
  let f0 : ℝ → ℝ := fun x =>
    lambertP x ^ n / lambertD x ^ (n + 1)
  let f1 : ℝ → ℝ := fun x =>
    lambertP x ^ (n + 1) / lambertD x ^ (n + 2)
  let f2 : ℝ → ℝ := fun x =>
    lambertP x ^ (n + 2) / lambertD x ^ (n + 3)
  let g : ℝ → ℝ := fun x =>
    ((n + 1 : ℕ) : ℝ) * f0 x -
      (2 * (n : ℝ) + 3) * Real.sqrt 2 * f1 x -
      2 * ((n + 2 : ℕ) : ℝ) * f2 x
  have hf0 : Continuous f0 := continuous_lambertMoment_integrand n
  have hf1 : Continuous f1 := continuous_lambertMoment_integrand (n + 1)
  have hf2 : Continuous f2 := continuous_lambertMoment_integrand (n + 2)
  have hg : Continuous g := by
    dsimp [g]
    exact
      ((hf0.const_mul (((n + 1 : ℕ) : ℝ))).sub
        (hf1.const_mul ((2 * (n : ℝ) + 3) * Real.sqrt 2))).sub
          (hf2.const_mul (2 * ((n + 2 : ℕ) : ℝ)))
  have hFTC : (∫ x in (0 : ℝ)..1, g x) = 0 := by
    calc
      (∫ x in (0 : ℝ)..1, g x) =
          ((fun x : ℝ =>
            (1 - 2 * x) * lambertP x ^ (n + 1) /
              lambertD x ^ (n + 2)) 1 -
            (fun x : ℝ =>
              (1 - 2 * x) * lambertP x ^ (n + 1) /
                lambertD x ^ (n + 2)) 0) := by
        refine intervalIntegral.integral_eq_sub_of_hasDerivAt
          (f := fun x : ℝ =>
            (1 - 2 * x) * lambertP x ^ (n + 1) /
              lambertD x ^ (n + 2))
          (f' := g) ?_ ?_
        · intro x _
          simpa [g, f0, f1, f2] using lambert_recurrence_aux_hasDerivAt n x
        · exact hg.continuousOn.intervalIntegrable
      _ = 0 := by
        simp [lambertP, lambertD]
  have hi0 : IntervalIntegrable f0 MeasureTheory.volume 0 1 :=
    hf0.continuousOn.intervalIntegrable
  have hi1 : IntervalIntegrable f1 MeasureTheory.volume 0 1 :=
    hf1.continuousOn.intervalIntegrable
  have hi2 : IntervalIntegrable f2 MeasureTheory.volume 0 1 :=
    hf2.continuousOn.intervalIntegrable
  have hInt :
      (∫ x in (0 : ℝ)..1, g x) =
        ((n + 1 : ℕ) : ℝ) * lambertMoment n -
          (2 * (n : ℝ) + 3) * Real.sqrt 2 * lambertMoment (n + 1) -
          2 * ((n + 2 : ℕ) : ℝ) * lambertMoment (n + 2) := by
    dsimp [g]
    rw [intervalIntegral.integral_sub
      ((hi0.const_mul _).sub (hi1.const_mul _)) (hi2.const_mul _)]
    rw [intervalIntegral.integral_sub (hi0.const_mul _) (hi1.const_mul _)]
    rw [intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul,
      intervalIntegral.integral_const_mul]
    rfl
  rw [hFTC] at hInt
  linarith

private def lambertErrorMagnitude (n : ℕ) : ℝ :=
  ((n ! : ℕ) : ℝ) * (Real.sqrt 2) ^ n / Real.sqrt 2 * lambertMoment n

private lemma lambertErrorMagnitude_zero :
    lambertErrorMagnitude 0 = Real.pi / 4 := by
  rw [lambertErrorMagnitude, lambertMoment_zero]
  norm_num

private lemma lambertErrorMagnitude_one :
    lambertErrorMagnitude 1 = 1 - Real.pi / 4 := by
  rw [lambertErrorMagnitude, lambertMoment_one]
  norm_num [Nat.factorial]

private lemma lambertErrorMagnitude_recurrence (n : ℕ) :
    lambertErrorMagnitude (n + 2) =
      ((n + 1 : ℕ) : ℝ) ^ 2 * lambertErrorMagnitude n -
        (2 * (n : ℝ) + 3) * lambertErrorMagnitude (n + 1) := by
  have hrec := lambertMoment_recurrence n
  have hf1 :
      ((((n + 1)! : ℕ) : ℝ)) =
        ((n + 1 : ℕ) : ℝ) * ((n ! : ℕ) : ℝ) := by
    exact_mod_cast Nat.factorial_succ n
  have hf2 :
      ((((n + 2)! : ℕ) : ℝ)) =
        ((n + 2 : ℕ) : ℝ) * ((n + 1 : ℕ) : ℝ) * ((n ! : ℕ) : ℝ) := by
    have hnat :
        (n + 2)! = (n + 2) * (n + 1) * n ! := by
      rw [show n + 2 = (n + 1) + 1 by omega, Nat.factorial_succ,
        Nat.factorial_succ]
      ring
    exact_mod_cast hnat
  have hp1 :
      (Real.sqrt 2) ^ (n + 1) = (Real.sqrt 2) ^ n * Real.sqrt 2 := by
    rw [pow_succ]
  have hp2 :
      (Real.sqrt 2) ^ (n + 2) = 2 * (Real.sqrt 2) ^ n := by
    calc
      (Real.sqrt 2) ^ (n + 2) =
          (Real.sqrt 2) ^ n * Real.sqrt 2 * Real.sqrt 2 := by
            rw [show n + 2 = (n + 1) + 1 by omega, pow_succ, hp1]
      _ = (Real.sqrt 2) ^ n * (Real.sqrt 2) ^ 2 := by ring
      _ = (Real.sqrt 2) ^ n * 2 := by rw [sqrt_two_sq]
      _ = 2 * (Real.sqrt 2) ^ n := by ring
  unfold lambertErrorMagnitude
  rw [hf1, hf2, hp1, hp2]
  linear_combination
    (((n + 1 : ℕ) : ℝ) * ((n ! : ℕ) : ℝ) *
      (Real.sqrt 2) ^ n / Real.sqrt 2) * hrec

private lemma lambert_remainder_eq_error :
    ∀ n : ℕ,
      ((lambertA n : ℤ) : ℝ) * (Real.pi / 4) - ((lambertB n : ℤ) : ℝ) =
        (-1 : ℝ) ^ n * lambertErrorMagnitude n
  | 0 => by
      rw [lambertErrorMagnitude_zero]
      norm_num [lambertA, lambertB]
  | 1 => by
      rw [lambertErrorMagnitude_one]
      norm_num [lambertA, lambertB]
  | (n + 2) => by
      rw [show lambertA (n + 2) =
          (2 * (n : ℤ) + 3) * lambertA (n + 1) +
            ((n : ℤ) + 1) ^ 2 * lambertA n from rfl]
      rw [show lambertB (n + 2) =
          (2 * (n : ℤ) + 3) * lambertB (n + 1) +
            ((n : ℤ) + 1) ^ 2 * lambertB n from rfl]
      push_cast
      rw [show
          ((2 * (n : ℝ) + 3) * (lambertA (n + 1) : ℝ) +
                ((n : ℝ) + 1) ^ 2 * (lambertA n : ℝ)) * (Real.pi / 4) -
              ((2 * (n : ℝ) + 3) * (lambertB (n + 1) : ℝ) +
                ((n : ℝ) + 1) ^ 2 * (lambertB n : ℝ)) =
            (2 * (n : ℝ) + 3) *
                ((lambertA (n + 1) : ℝ) * (Real.pi / 4) -
                  (lambertB (n + 1) : ℝ)) +
              ((n : ℝ) + 1) ^ 2 *
                ((lambertA n : ℝ) * (Real.pi / 4) -
                  (lambertB n : ℝ)) by ring]
      rw [lambert_remainder_eq_error (n + 1), lambert_remainder_eq_error n,
        lambertErrorMagnitude_recurrence n]
      rw [show n + 2 = (n + 1) + 1 by omega, pow_succ, pow_succ]
      push_cast
      ring

private lemma lambertMoment_nonneg (n : ℕ) : 0 ≤ lambertMoment n := by
  apply intervalIntegral.integral_nonneg (by norm_num)
  intro x hx
  exact div_nonneg (pow_nonneg (lambertP_nonneg hx) _)
    (pow_nonneg (lambertD_pos x).le _)

private lemma lambertP_div_lambertD_le {x : ℝ}
    (_hx : x ∈ Set.Icc (0 : ℝ) 1) :
    lambertP x / lambertD x ≤ 1 / (2 + Real.sqrt 2) := by
  have hD := lambertD_pos x
  have hs : 0 < 2 + Real.sqrt 2 := by positivity
  apply (div_le_div_iff₀ hD hs).2
  dsimp [lambertD]
  nlinarith [four_mul_lambertP_le_one x]

private lemma lambertMoment_le_geometric (n : ℕ) :
    lambertMoment n ≤
      (1 / (2 + Real.sqrt 2)) ^ n * lambertMoment 0 := by
  let r : ℝ := 1 / (2 + Real.sqrt 2)
  have hr0 : 0 ≤ r := by
    dsimp [r]
    positivity
  have hleft :
      IntervalIntegrable
        (fun x : ℝ => lambertP x ^ n / lambertD x ^ (n + 1))
        MeasureTheory.volume 0 1 :=
    (continuous_lambertMoment_integrand n).continuousOn.intervalIntegrable
  have hone : Continuous (fun x : ℝ => 1 / lambertD x) :=
    continuous_const.div continuous_lambertD fun x => (lambertD_pos x).ne'
  have hright :
      IntervalIntegrable (fun x : ℝ => r ^ n * (1 / lambertD x))
        MeasureTheory.volume 0 1 :=
    (hone.const_mul (r ^ n)).continuousOn.intervalIntegrable
  have hpoint :
      ∀ x ∈ Set.Icc (0 : ℝ) 1,
        lambertP x ^ n / lambertD x ^ (n + 1) ≤
          r ^ n * (1 / lambertD x) := by
    intro x hx
    have hD := lambertD_pos x
    have hp := lambertP_nonneg hx
    have hratio0 : 0 ≤ lambertP x / lambertD x :=
      div_nonneg hp hD.le
    have hratiole : lambertP x / lambertD x ≤ r := by
      exact lambertP_div_lambertD_le hx
    have hpow : (lambertP x / lambertD x) ^ n ≤ r ^ n :=
      pow_le_pow_left₀ hratio0 hratiole n
    have heq :
        lambertP x ^ n / lambertD x ^ (n + 1) =
          (lambertP x / lambertD x) ^ n * (1 / lambertD x) := by
      rw [div_pow]
      field_simp [(lambertD_pos x).ne']
      ring
    rw [heq]
    exact mul_le_mul_of_nonneg_right hpow (one_div_nonneg.mpr hD.le)
  calc
    lambertMoment n =
        ∫ x in (0 : ℝ)..1, lambertP x ^ n / lambertD x ^ (n + 1) := rfl
    _ ≤ ∫ x in (0 : ℝ)..1, r ^ n * (1 / lambertD x) :=
      intervalIntegral.integral_mono_on (by norm_num) hleft hright hpoint
    _ = r ^ n * (∫ x in (0 : ℝ)..1, 1 / lambertD x) := by
      rw [intervalIntegral.integral_const_mul]
    _ = r ^ n * lambertMoment 0 := by
      simp [lambertMoment]

private lemma sqrt_two_mul_inv_add :
    Real.sqrt 2 * (1 / (2 + Real.sqrt 2)) = Real.sqrt 2 - 1 := by
  have hden : 2 + Real.sqrt 2 ≠ 0 := by positivity
  rw [mul_one_div]
  apply (div_eq_iff hden).2
  nlinarith [sqrt_two_sq]

private lemma lambertErrorMagnitude_nonneg (n : ℕ) :
    0 ≤ lambertErrorMagnitude n := by
  unfold lambertErrorMagnitude
  positivity [lambertMoment_nonneg n, sqrt_two_pos]

private lemma lambertErrorMagnitude_le_geometric (n : ℕ) :
    lambertErrorMagnitude n ≤
      ((n ! : ℕ) : ℝ) * (Real.sqrt 2 - 1) ^ n * (Real.pi / 4) := by
  have hm := lambertMoment_le_geometric n
  have hc :
      0 ≤ ((n ! : ℕ) : ℝ) * (Real.sqrt 2) ^ n / Real.sqrt 2 := by
    positivity
  unfold lambertErrorMagnitude
  calc
    ((n ! : ℕ) : ℝ) * (Real.sqrt 2) ^ n / Real.sqrt 2 * lambertMoment n ≤
        ((n ! : ℕ) : ℝ) * (Real.sqrt 2) ^ n / Real.sqrt 2 *
          ((1 / (2 + Real.sqrt 2)) ^ n * lambertMoment 0) :=
      mul_le_mul_of_nonneg_left hm hc
    _ = ((n ! : ℕ) : ℝ) * (Real.sqrt 2 - 1) ^ n * (Real.pi / 4) := by
      rw [lambertMoment_zero]
      calc
        ((n ! : ℕ) : ℝ) * (Real.sqrt 2) ^ n / Real.sqrt 2 *
            ((1 / (2 + Real.sqrt 2)) ^ n *
              (Real.sqrt 2 * (Real.pi / 4))) =
            ((n ! : ℕ) : ℝ) *
              ((Real.sqrt 2) ^ n * (1 / (2 + Real.sqrt 2)) ^ n) *
              (Real.pi / 4) := by
                field_simp [sqrt_two_pos.ne']
        _ = ((n ! : ℕ) : ℝ) *
              (Real.sqrt 2 * (1 / (2 + Real.sqrt 2))) ^ n *
              (Real.pi / 4) := by rw [mul_pow]
        _ = ((n ! : ℕ) : ℝ) * (Real.sqrt 2 - 1) ^ n *
              (Real.pi / 4) := by rw [sqrt_two_mul_inv_add]

private lemma factorial_le_lambertA_real :
    ∀ n : ℕ, ((n ! : ℕ) : ℝ) ≤ ((lambertA n : ℤ) : ℝ)
  | 0 => by norm_num [lambertA]
  | 1 => by norm_num [lambertA]
  | (n + 2) => by
      have ih := factorial_le_lambertA_real (n + 1)
      have hcoef : 0 ≤ 2 * (n : ℝ) + 3 := by positivity
      have hmul :
          (2 * (n : ℝ) + 3) * (((n + 1)! : ℕ) : ℝ) ≤
            (2 * (n : ℝ) + 3) * ((lambertA (n + 1) : ℤ) : ℝ) :=
        mul_le_mul_of_nonneg_left ih hcoef
      have hcoefle : ((n + 2 : ℕ) : ℝ) ≤ 2 * (n : ℝ) + 3 := by
        push_cast
        linarith
      have hfacnon : 0 ≤ ((((n + 1)! : ℕ) : ℝ)) := by positivity
      have hfac :
          ((((n + 2)! : ℕ) : ℝ)) =
            ((n + 2 : ℕ) : ℝ) * (((n + 1)! : ℕ) : ℝ) := by
        exact_mod_cast Nat.factorial_succ (n + 1)
      have hextra :
          0 ≤ ((n : ℝ) + 1) ^ 2 * ((lambertA n : ℤ) : ℝ) := by
        have hA : 0 ≤ ((lambertA n : ℤ) : ℝ) := by
          exact_mod_cast (lambertA_pos n).le
        positivity
      rw [show lambertA (n + 2) =
          (2 * (n : ℤ) + 3) * lambertA (n + 1) +
            ((n : ℤ) + 1) ^ 2 * lambertA n from rfl]
      push_cast
      calc
        (((n + 2)! : ℕ) : ℝ) =
            ((n + 2 : ℕ) : ℝ) * (((n + 1)! : ℕ) : ℝ) := hfac
        _ ≤ (2 * (n : ℝ) + 3) * (((n + 1)! : ℕ) : ℝ) :=
          mul_le_mul_of_nonneg_right hcoefle hfacnon
        _ ≤ (2 * (n : ℝ) + 3) * ((lambertA (n + 1) : ℤ) : ℝ) := hmul
        _ ≤ (2 * (n : ℝ) + 3) * ((lambertA (n + 1) : ℤ) : ℝ) +
            ((n : ℝ) + 1) ^ 2 * ((lambertA n : ℤ) : ℝ) :=
          le_add_of_nonneg_right hextra

private lemma lambert_remainder_div_abs_le (n : ℕ) :
    |(((lambertA n : ℤ) : ℝ) * (Real.pi / 4) - ((lambertB n : ℤ) : ℝ)) /
        ((lambertA n : ℤ) : ℝ)| ≤
      (Real.pi / 4) * (Real.sqrt 2 - 1) ^ n := by
  have hApos : 0 < ((lambertA n : ℤ) : ℝ) := by
    exact_mod_cast lambertA_pos n
  have hfac := factorial_le_lambertA_real n
  have hT := lambertErrorMagnitude_le_geometric n
  have hT0 := lambertErrorMagnitude_nonneg n
  have hr0 : 0 ≤ (Real.sqrt 2 - 1) ^ n := by
    exact pow_nonneg (sub_nonneg.mpr one_lt_sqrt_two.le) n
  have hpi0 : 0 ≤ Real.pi / 4 := by positivity
  have hdiv :
      lambertErrorMagnitude n / ((lambertA n : ℤ) : ℝ) ≤
        (Real.pi / 4) * (Real.sqrt 2 - 1) ^ n := by
    apply (div_le_iff₀ hApos).2
    calc
      lambertErrorMagnitude n ≤
          ((n ! : ℕ) : ℝ) * (Real.sqrt 2 - 1) ^ n *
            (Real.pi / 4) := hT
      _ = ((n ! : ℕ) : ℝ) *
          ((Real.sqrt 2 - 1) ^ n * (Real.pi / 4)) := by ring
      _ ≤ ((lambertA n : ℤ) : ℝ) *
          ((Real.sqrt 2 - 1) ^ n * (Real.pi / 4)) :=
        mul_le_mul_of_nonneg_right hfac (mul_nonneg hr0 hpi0)
      _ = (Real.pi / 4) * (Real.sqrt 2 - 1) ^ n *
          ((lambertA n : ℤ) : ℝ) := by ring
  calc
    |(((lambertA n : ℤ) : ℝ) * (Real.pi / 4) - ((lambertB n : ℤ) : ℝ)) /
        ((lambertA n : ℤ) : ℝ)| =
        lambertErrorMagnitude n / ((lambertA n : ℤ) : ℝ) := by
      rw [lambert_remainder_eq_error]
      rw [abs_div, abs_mul, abs_pow, abs_neg, abs_one, one_pow,
        abs_of_nonneg hT0, abs_of_pos hApos]
      ring
    _ ≤ (Real.pi / 4) * (Real.sqrt 2 - 1) ^ n := hdiv

private lemma lambert_remainder_div_tendsto_zero :
    Tendsto
      (fun n =>
        (((lambertA n : ℤ) : ℝ) * (Real.pi / 4) - ((lambertB n : ℤ) : ℝ)) /
          ((lambertA n : ℤ) : ℝ))
      atTop (𝓝 0) := by
  rw [tendsto_zero_iff_norm_tendsto_zero]
  refine squeeze_zero
    (g := fun n => (Real.pi / 4) * (Real.sqrt 2 - 1) ^ n)
    (fun _ => norm_nonneg _) (fun n => ?_) ?_
  · rw [Real.norm_eq_abs]
    exact lambert_remainder_div_abs_le n
  · simpa using
      (tendsto_pow_atTop_nhds_zero_of_lt_one
        (sub_nonneg.mpr one_lt_sqrt_two.le) (by nlinarith [sqrt_two_lt_two])).const_mul
          (Real.pi / 4)

theorem lambertB_div_lambertA_tendsto_pi_div_four :
    Tendsto
      (fun n => ((lambertB n : ℤ) : ℝ) / ((lambertA n : ℤ) : ℝ))
      atTop (𝓝 (Real.pi / 4)) := by
  have hconst :
      Tendsto (fun _ : ℕ => Real.pi / 4) atTop (𝓝 (Real.pi / 4)) :=
    tendsto_const_nhds
  have h := hconst.sub lambert_remainder_div_tendsto_zero
  have h' :
      Tendsto
        (fun n =>
          Real.pi / 4 -
            ((((lambertA n : ℤ) : ℝ) * (Real.pi / 4) -
              ((lambertB n : ℤ) : ℝ)) / ((lambertA n : ℤ) : ℝ)))
        atTop (𝓝 (Real.pi / 4)) := by
    simpa using h
  refine h'.congr' (Eventually.of_forall fun n => ?_)
  have hA : ((lambertA n : ℤ) : ℝ) ≠ 0 := by
    exact_mod_cast (lambertA_pos n).ne'
  field_simp
  ring

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

/-! ## Main theorem -/

theorem problem23_pi_add_e :
    Tendsto (fun m => ((challengeP m : ℤ) : ℝ) / ((challengeQ m : ℤ) : ℝ)) atTop
      (𝓝 (Real.pi + Real.exp 1)) := by
  have hsum :
      Tendsto (fun m => 4 * (((lambertB m : ℤ) : ℝ) / ((lambertA m : ℤ) : ℝ))
          + ((m ! : ℕ) : ℝ) / ((derang m : ℤ) : ℝ)) atTop
        (𝓝 (4 * (Real.pi / 4) + Real.exp 1)) :=
    (lambertB_div_lambertA_tendsto_pi_div_four.const_mul 4).add
      factorial_div_derang_tendsto_exp_one
  have hval : 4 * (Real.pi / 4) + Real.exp 1 = Real.pi + Real.exp 1 := by ring
  rw [hval] at hsum
  refine hsum.congr' ?_
  filter_upwards [eventually_ge_atTop 2] with m hm
  obtain ⟨j, rfl⟩ : ∃ j, m = j + 2 := ⟨m - 2, by omega⟩
  exact (ratio_split j).symm

end RamanujanChallenge.P23

end
