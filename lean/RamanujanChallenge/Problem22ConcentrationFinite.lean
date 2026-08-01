import RamanujanChallenge.Problem22Concentration
import Mathlib.Algebra.Order.BigOperators.Ring.Finset

/-!
# Problem 2.2: finite-moment concentration closure

This file completes the finite Stein/moment route for the positive Rivoal
weights.  The constants are deliberately loose: `3`, `4`, and `100` replace
the sharper paper constants `5/2`, `7/2`, and `81`.
-/

noncomputable section

open Filter Topology Real
open scoped BigOperators

namespace RamanujanChallenge.P22

private theorem rivoalRealQ22_pos (n : ℕ) :
    0 < ((rivoalExplicitQ22 n : ℚ) : ℝ) := by
  exact_mod_cast rivoalExplicitQ22_pos n

private theorem rivoalSupportBounds22 {n k : ℕ}
    (hk : k ∈ Finset.range (n + 1)) :
    0 ≤ (k : ℝ) ∧ (k : ℝ) ≤ (n : ℝ) := by
  constructor
  · positivity
  · exact_mod_cast (Nat.lt_succ_iff.mp (Finset.mem_range.mp hk))

private theorem rivoalBirth22_nonneg_on_support {n k : ℕ}
    (_hk : k ∈ Finset.range (n + 1)) :
    0 ≤ rivoalBirth22 n k := by
  unfold rivoalBirth22
  exact mul_nonneg (by positivity) (sq_nonneg _)

private theorem rivoalBirth22_le_five_cube {n k : ℕ}
    (hn : 1 ≤ n) (hk : k ∈ Finset.range (n + 1)) :
    rivoalBirth22 n k ≤ 5 * (n : ℝ) ^ 3 := by
  obtain ⟨hk0, hkn⟩ := rivoalSupportBounds22 hk
  have hnR : (1 : ℝ) ≤ (n : ℝ) := by exact_mod_cast hn
  have hd0 : 0 ≤ (n : ℝ) - (k : ℝ) := sub_nonneg.mpr hkn
  have hdle : (n : ℝ) - (k : ℝ) ≤ (n : ℝ) := by linarith
  have hsquare :
      ((n : ℝ) - (k : ℝ)) ^ 2 ≤ (n : ℝ) ^ 2 := by
    have hprod :
        0 ≤ ((n : ℝ) - ((n : ℝ) - (k : ℝ))) *
          ((n : ℝ) + ((n : ℝ) - (k : ℝ))) :=
      mul_nonneg (sub_nonneg.mpr hdle) (add_nonneg (by positivity) hd0)
    nlinarith
  have hlinear :
      2 * (n : ℝ) + (k : ℝ) + 2 ≤ 5 * (n : ℝ) := by
    linarith
  unfold rivoalBirth22
  calc
    (2 * (n : ℝ) + (k : ℝ) + 2) *
        ((n : ℝ) - (k : ℝ)) ^ 2 ≤
      (5 * (n : ℝ)) * (n : ℝ) ^ 2 :=
        mul_le_mul hlinear hsquare (sq_nonneg _) (by positivity)
    _ = 5 * (n : ℝ) ^ 3 := by ring

/-- The constant-function Stein identity is exact birth/death balance. -/
theorem rivoalBirthDeathBalance22 (n : ℕ) :
    (∑ k ∈ Finset.range (n + 1),
        rivoalRealWeight22 n k * rivoalDeath22 n k) =
      ∑ k ∈ Finset.range (n + 1),
        rivoalRealWeight22 n k * rivoalBirth22 n k := by
  have h := rivoalWeightSteinShift22 n (fun _ => (1 : ℝ))
  simpa using h.symm

/-- Weighted cubic moment.  The loose constant `3` avoids division by `2n`. -/
theorem rivoalWeightedCube22_le {n : ℕ} (hn : 1 ≤ n) :
    (∑ k ∈ Finset.range (n + 1),
      rivoalRealWeight22 n k * (k : ℝ) ^ 3) ≤
      3 * (n : ℝ) ^ 2 * ((rivoalExplicitQ22 n : ℚ) : ℝ) := by
  have hnR : (0 : ℝ) < (n : ℝ) := by
    exact_mod_cast (lt_of_lt_of_le Nat.zero_lt_one hn)
  have hQ0 : 0 ≤ ((rivoalExplicitQ22 n : ℚ) : ℝ) :=
    (rivoalRealQ22_pos n).le
  have hlower :
      (2 * (n : ℝ)) *
          (∑ k ∈ Finset.range (n + 1),
            rivoalRealWeight22 n k * (k : ℝ) ^ 3) ≤
        ∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k * rivoalDeath22 n k := by
    calc
      (2 * (n : ℝ)) *
          (∑ k ∈ Finset.range (n + 1),
            rivoalRealWeight22 n k * (k : ℝ) ^ 3) =
        ∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k *
            ((2 * (n : ℝ)) * (k : ℝ) ^ 3) := by
              rw [Finset.mul_sum]
              apply Finset.sum_congr rfl
              intro k hk
              ring
      _ ≤ ∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k * rivoalDeath22 n k := by
            apply Finset.sum_le_sum
            intro k hk
            apply mul_le_mul_of_nonneg_left
            · unfold rivoalDeath22
              have hk0 : 0 ≤ (k : ℝ) := by positivity
              exact mul_le_mul_of_nonneg_right
                (by linarith : 2 * (n : ℝ) ≤ 2 * (n : ℝ) + (k : ℝ))
                (pow_nonneg hk0 3)
            · exact rivoalRealWeight22_nonneg n k
  have hupper :
      (∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k * rivoalBirth22 n k) ≤
        5 * (n : ℝ) ^ 3 * ((rivoalExplicitQ22 n : ℚ) : ℝ) := by
    calc
      (∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k * rivoalBirth22 n k) ≤
        ∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k * (5 * (n : ℝ) ^ 3) := by
            apply Finset.sum_le_sum
            intro k hk
            exact mul_le_mul_of_nonneg_left
              (rivoalBirth22_le_five_cube hn hk)
              (rivoalRealWeight22_nonneg n k)
      _ = 5 * (n : ℝ) ^ 3 * ((rivoalExplicitQ22 n : ℚ) : ℝ) := by
            rw [← Finset.sum_mul, rivoalRealWeight22_sum]
            ring
  have hmain :
      (2 * (n : ℝ)) *
          (∑ k ∈ Finset.range (n + 1),
            rivoalRealWeight22 n k * (k : ℝ) ^ 3) ≤
        5 * (n : ℝ) ^ 3 * ((rivoalExplicitQ22 n : ℚ) : ℝ) :=
    hlower.trans ((rivoalBirthDeathBalance22 n).trans_le hupper)
  have hscale :
      5 * (n : ℝ) ^ 3 * ((rivoalExplicitQ22 n : ℚ) : ℝ) ≤
        (2 * (n : ℝ)) *
          (3 * (n : ℝ) ^ 2 * ((rivoalExplicitQ22 n : ℚ) : ℝ)) := by
    have hnonneg :
        0 ≤ (n : ℝ) ^ 3 * ((rivoalExplicitQ22 n : ℚ) : ℝ) :=
      mul_nonneg (pow_nonneg (by positivity) 3) hQ0
    calc
      5 * (n : ℝ) ^ 3 * ((rivoalExplicitQ22 n : ℚ) : ℝ) ≤
          6 * (n : ℝ) ^ 3 * ((rivoalExplicitQ22 n : ℚ) : ℝ) := by
            nlinarith
      _ = (2 * (n : ℝ)) *
          (3 * (n : ℝ) ^ 2 * ((rivoalExplicitQ22 n : ℚ) : ℝ)) := by
            ring
  exact le_of_mul_le_mul_of_pos_left (hmain.trans hscale)
    (show 0 < 2 * (n : ℝ) by positivity)

private theorem sq_le_sq_add_cube_div22 {x a : ℝ}
    (hx : 0 ≤ x) (ha : 0 < a) :
    x ^ 2 ≤ a ^ 2 + x ^ 3 / a := by
  rcases le_total x a with hxa | hax
  · have hprod : 0 ≤ (a - x) * (a + x) :=
      mul_nonneg (sub_nonneg.mpr hxa) (add_nonneg ha.le hx)
    have hsq : x ^ 2 ≤ a ^ 2 := by nlinarith
    exact hsq.trans (le_add_of_nonneg_right
      (div_nonneg (pow_nonneg hx 3) ha.le))
  · have hmul : x ^ 2 * a ≤ x ^ 3 := by
      calc
        x ^ 2 * a ≤ x ^ 2 * x :=
          mul_le_mul_of_nonneg_left hax (sq_nonneg x)
        _ = x ^ 3 := by ring
    have hdiv : x ^ 2 ≤ x ^ 3 / a := (le_div_iff₀ ha).2 hmul
    exact hdiv.trans (le_add_of_nonneg_left (sq_nonneg a))

/-- Weighted quadratic moment obtained by interpolation at `sqrt n`. -/
theorem rivoalWeightedSquare22_le {n : ℕ} (hn : 1 ≤ n) :
    (∑ k ∈ Finset.range (n + 1),
      rivoalRealWeight22 n k * (k : ℝ) ^ 2) ≤
      4 * (n : ℝ) * Real.sqrt (n : ℝ) *
        ((rivoalExplicitQ22 n : ℚ) : ℝ) := by
  let s : ℝ := Real.sqrt (n : ℝ)
  have hnR : (1 : ℝ) ≤ (n : ℝ) := by exact_mod_cast hn
  have hs : 0 < s := Real.sqrt_pos.2 (by linarith)
  have hs2 : s ^ 2 = (n : ℝ) := by
    exact Real.sq_sqrt (by linarith)
  have hs1 : 1 ≤ s := by
    nlinarith [hs2]
  have hQ0 : 0 ≤ ((rivoalExplicitQ22 n : ℚ) : ℝ) :=
    (rivoalRealQ22_pos n).le
  have hsum :
      (∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k * (k : ℝ) ^ 2) ≤
        ∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k *
            (s ^ 2 + (k : ℝ) ^ 3 / s) := by
    apply Finset.sum_le_sum
    intro k hk
    exact mul_le_mul_of_nonneg_left
      (sq_le_sq_add_cube_div22 (by positivity) hs)
      (rivoalRealWeight22_nonneg n k)
  have hdecomp :
      (∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k *
            (s ^ 2 + (k : ℝ) ^ 3 / s)) =
        s ^ 2 * ((rivoalExplicitQ22 n : ℚ) : ℝ) +
          (∑ k ∈ Finset.range (n + 1),
            rivoalRealWeight22 n k * (k : ℝ) ^ 3) / s := by
    calc
      _ = (∑ k ∈ Finset.range (n + 1),
          (rivoalRealWeight22 n k * s ^ 2 +
            (rivoalRealWeight22 n k * (k : ℝ) ^ 3) / s)) := by
            apply Finset.sum_congr rfl
            intro k hk
            ring
      _ = (∑ k ∈ Finset.range (n + 1),
            rivoalRealWeight22 n k * s ^ 2) +
          ∑ k ∈ Finset.range (n + 1),
            (rivoalRealWeight22 n k * (k : ℝ) ^ 3) / s := by
              rw [Finset.sum_add_distrib]
      _ = _ := by
            rw [← Finset.sum_mul, ← Finset.sum_div,
              rivoalRealWeight22_sum]
            ring
  have hcube := rivoalWeightedCube22_le hn
  have hdiv :
      (∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k * (k : ℝ) ^ 3) / s ≤
        (3 * (n : ℝ) ^ 2 * ((rivoalExplicitQ22 n : ℚ) : ℝ)) / s :=
    div_le_div_of_nonneg_right hcube hs.le
  have hid :
      (3 * (n : ℝ) ^ 2 * ((rivoalExplicitQ22 n : ℚ) : ℝ)) / s =
        3 * (n : ℝ) * s * ((rivoalExplicitQ22 n : ℚ) : ℝ) := by
    rw [← hs2]
    field_simp [hs.ne']
  have hbase :
      (n : ℝ) * ((rivoalExplicitQ22 n : ℚ) : ℝ) ≤
        (n : ℝ) * s * ((rivoalExplicitQ22 n : ℚ) : ℝ) := by
    have hfac :
        0 ≤ (n : ℝ) * ((rivoalExplicitQ22 n : ℚ) : ℝ) :=
      mul_nonneg (by positivity) hQ0
    nlinarith
  calc
    _ ≤ _ := hsum
    _ = (n : ℝ) * ((rivoalExplicitQ22 n : ℚ) : ℝ) +
          (∑ k ∈ Finset.range (n + 1),
            rivoalRealWeight22 n k * (k : ℝ) ^ 3) / s := by
          rw [hdecomp, hs2]
    _ ≤ (n : ℝ) * ((rivoalExplicitQ22 n : ℚ) : ℝ) +
        (3 * (n : ℝ) ^ 2 * ((rivoalExplicitQ22 n : ℚ) : ℝ)) / s :=
      add_le_add (le_refl _) hdiv
    _ = (n : ℝ) * ((rivoalExplicitQ22 n : ℚ) : ℝ) +
        3 * (n : ℝ) * s * ((rivoalExplicitQ22 n : ℚ) : ℝ) := by rw [hid]
    _ ≤ 4 * (n : ℝ) * s * ((rivoalExplicitQ22 n : ℚ) : ℝ) := by
      nlinarith
    _ = 4 * (n : ℝ) * Real.sqrt (n : ℝ) *
        ((rivoalExplicitQ22 n : ℚ) : ℝ) := rfl

/-- Polynomial increment of the saddle observable. -/
theorem rivoalSaddleError22_succ (n k : ℕ) :
    rivoalSaddleError22 n (k + 1) - rivoalSaddleError22 n k =
      3 * (k : ℝ) ^ 2 + (k : ℝ) + 2 * (n : ℝ) := by
  simp only [rivoalSaddleError22, Nat.cast_add, Nat.cast_one]
  ring

/-- Drift factorization needed by the exact second-moment identity. -/
theorem rivoalBirthSubDeath22 (n k : ℕ) :
    rivoalBirth22 n k - rivoalDeath22 n k =
      2 * ((n : ℝ) - (k : ℝ)) ^ 2 -
        (2 * (n : ℝ) + (k : ℝ)) * rivoalSaddleError22 n k := by
  simp only [rivoalBirth22, rivoalDeath22, rivoalSaddleError22]
  ring

/-- Exact finite saddle second-moment certificate. -/
theorem rivoalSaddleSecondMoment22_exact (n : ℕ) :
    (∑ k ∈ Finset.range (n + 1),
      rivoalRealWeight22 n k *
        ((2 * (n : ℝ) + (k : ℝ)) * rivoalSaddleError22 n k ^ 2)) =
      ∑ k ∈ Finset.range (n + 1),
        rivoalRealWeight22 n k *
          (2 * ((n : ℝ) - (k : ℝ)) ^ 2 * rivoalSaddleError22 n k +
            rivoalBirth22 n k *
              (3 * (k : ℝ) ^ 2 + (k : ℝ) + 2 * (n : ℝ))) := by
  have h := rivoalWeightStein22 n (rivoalSaddleError22 n)
  have h' :
      (∑ k ∈ Finset.range (n + 1),
        rivoalRealWeight22 n k *
          ((2 * ((n : ℝ) - (k : ℝ)) ^ 2 * rivoalSaddleError22 n k +
              rivoalBirth22 n k *
                (3 * (k : ℝ) ^ 2 + (k : ℝ) + 2 * (n : ℝ))) -
            ((2 * (n : ℝ) + (k : ℝ)) *
              rivoalSaddleError22 n k ^ 2))) = 0 := by
    calc
      _ = (∑ k ∈ Finset.range (n + 1),
        rivoalRealWeight22 n k *
          (rivoalBirth22 n k *
              (rivoalSaddleError22 n (k + 1) -
                rivoalSaddleError22 n k) +
            (rivoalBirth22 n k - rivoalDeath22 n k) *
              rivoalSaddleError22 n k)) := by
              apply Finset.sum_congr rfl
              intro k hk
              rw [rivoalSaddleError22_succ, rivoalBirthSubDeath22]
              ring
      _ = 0 := h
  rw [← sub_eq_zero, ← Finset.sum_sub_distrib]
  calc
    (∑ k ∈ Finset.range (n + 1),
      (rivoalRealWeight22 n k *
          ((2 * (n : ℝ) + (k : ℝ)) * rivoalSaddleError22 n k ^ 2) -
        rivoalRealWeight22 n k *
          (2 * ((n : ℝ) - (k : ℝ)) ^ 2 * rivoalSaddleError22 n k +
            rivoalBirth22 n k *
              (3 * (k : ℝ) ^ 2 + (k : ℝ) + 2 * (n : ℝ))))) =
      - ∑ k ∈ Finset.range (n + 1),
        rivoalRealWeight22 n k *
          ((2 * ((n : ℝ) - (k : ℝ)) ^ 2 * rivoalSaddleError22 n k +
              rivoalBirth22 n k *
                (3 * (k : ℝ) ^ 2 + (k : ℝ) + 2 * (n : ℝ))) -
            ((2 * (n : ℝ) + (k : ℝ)) *
              rivoalSaddleError22 n k ^ 2)) := by
                rw [← Finset.sum_neg_distrib]
                apply Finset.sum_congr rfl
                intro k hk
                ring
    _ = 0 := by rw [h', neg_zero]

private theorem two_mul_le_mul_sq_add_sq_div22
    {a b c : ℝ} (hc : 0 < c) :
    2 * a * b ≤ c * a ^ 2 + b ^ 2 / c := by
  have hmul :
      c * (2 * a * b) ≤ c * (c * a ^ 2 + b ^ 2 / c) := by
    rw [mul_add, mul_div_cancel₀ _ hc.ne']
    nlinarith [sq_nonneg (c * a - b)]
  exact le_of_mul_le_mul_of_pos_left hmul hc

/-- Explicit weighted saddle second moment. -/
theorem rivoalSaddleSecondMoment22_le {n : ℕ} (hn : 1 ≤ n) :
    (∑ k ∈ Finset.range (n + 1),
      rivoalRealWeight22 n k * rivoalSaddleError22 n k ^ 2) ≤
      100 * (n : ℝ) ^ 3 * Real.sqrt (n : ℝ) *
        ((rivoalExplicitQ22 n : ℚ) : ℝ) := by
  let Q : ℝ := ((rivoalExplicitQ22 n : ℚ) : ℝ)
  let S2 : ℝ := ∑ k ∈ Finset.range (n + 1),
    rivoalRealWeight22 n k * rivoalSaddleError22 n k ^ 2
  let K2 : ℝ := ∑ k ∈ Finset.range (n + 1),
    rivoalRealWeight22 n k * (k : ℝ) ^ 2
  let s : ℝ := Real.sqrt (n : ℝ)
  have hnR : (0 : ℝ) < (n : ℝ) := by
    exact_mod_cast (lt_of_lt_of_le Nat.zero_lt_one hn)
  have hQ : 0 < Q := by exact rivoalRealQ22_pos n
  have hs : 0 < s := Real.sqrt_pos.2 hnR
  have hs2 : s ^ 2 = (n : ℝ) := Real.sq_sqrt hnR.le
  have hs1 : 1 ≤ s := by nlinarith [hs2]
  have hleft :
      (2 * (n : ℝ)) * S2 ≤
        ∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k *
            ((2 * (n : ℝ) + (k : ℝ)) * rivoalSaddleError22 n k ^ 2) := by
    dsimp only [S2]
    calc
      (2 * (n : ℝ)) *
          (∑ k ∈ Finset.range (n + 1),
            rivoalRealWeight22 n k * rivoalSaddleError22 n k ^ 2) =
        ∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k *
            ((2 * (n : ℝ)) * rivoalSaddleError22 n k ^ 2) := by
              rw [Finset.mul_sum]
              apply Finset.sum_congr rfl
              intro k hk
              ring
      _ ≤ _ := by
        apply Finset.sum_le_sum
        intro k hk
        apply mul_le_mul_of_nonneg_left
        · have hk0 : 0 ≤ (k : ℝ) := by positivity
          exact mul_le_mul_of_nonneg_right
            (by linarith : 2 * (n : ℝ) ≤ 2 * (n : ℝ) + (k : ℝ))
            (sq_nonneg _)
        · exact rivoalRealWeight22_nonneg n k
  have hterm1 :
      (∑ k ∈ Finset.range (n + 1),
        rivoalRealWeight22 n k *
          (2 * ((n : ℝ) - (k : ℝ)) ^ 2 * rivoalSaddleError22 n k)) ≤
        (n : ℝ) * S2 + (n : ℝ) ^ 3 * Q := by
    have hp :
        (∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k *
            (2 * ((n : ℝ) - (k : ℝ)) ^ 2 *
              rivoalSaddleError22 n k)) ≤
          ∑ k ∈ Finset.range (n + 1),
            rivoalRealWeight22 n k *
              ((n : ℝ) * rivoalSaddleError22 n k ^ 2 +
                (n : ℝ) ^ 3) := by
      apply Finset.sum_le_sum
      intro k hk
      obtain ⟨hk0, hkn⟩ := rivoalSupportBounds22 hk
      have hd0 : 0 ≤ (n : ℝ) - (k : ℝ) := sub_nonneg.mpr hkn
      have hdle : (n : ℝ) - (k : ℝ) ≤ (n : ℝ) := by linarith
      have hd4 : ((n : ℝ) - (k : ℝ)) ^ 4 ≤ (n : ℝ) ^ 4 :=
        pow_le_pow_left₀ hd0 hdle 4
      have hdiv :
          (((n : ℝ) - (k : ℝ)) ^ 2) ^ 2 / (n : ℝ) ≤
            (n : ℝ) ^ 3 := by
        apply (div_le_iff₀ hnR).2
        calc
          (((n : ℝ) - (k : ℝ)) ^ 2) ^ 2 =
              ((n : ℝ) - (k : ℝ)) ^ 4 := by ring
          _ ≤ (n : ℝ) ^ 4 := hd4
          _ = (n : ℝ) ^ 3 * (n : ℝ) := by ring
      apply mul_le_mul_of_nonneg_left
      · exact (two_mul_le_mul_sq_add_sq_div22
          (a := rivoalSaddleError22 n k)
          (b := ((n : ℝ) - (k : ℝ)) ^ 2) hnR).trans
            (add_le_add (le_refl _) hdiv)
      · exact rivoalRealWeight22_nonneg n k
    calc
      _ ≤ _ := hp
      _ = (n : ℝ) * S2 + (n : ℝ) ^ 3 * Q := by
        dsimp only [S2, Q]
        calc
          _ = (∑ k ∈ Finset.range (n + 1),
              (rivoalRealWeight22 n k *
                  ((n : ℝ) * rivoalSaddleError22 n k ^ 2) +
                rivoalRealWeight22 n k * (n : ℝ) ^ 3)) := by
                  apply Finset.sum_congr rfl
                  intro k hk
                  ring
          _ = (∑ k ∈ Finset.range (n + 1),
                rivoalRealWeight22 n k *
                  ((n : ℝ) * rivoalSaddleError22 n k ^ 2)) +
              ∑ k ∈ Finset.range (n + 1),
                rivoalRealWeight22 n k * (n : ℝ) ^ 3 := by
                  rw [Finset.sum_add_distrib]
          _ = _ := by
            rw [← Finset.sum_mul, rivoalRealWeight22_sum]
            rw [show (∑ k ∈ Finset.range (n + 1),
                rivoalRealWeight22 n k *
                  ((n : ℝ) * rivoalSaddleError22 n k ^ 2)) =
                (n : ℝ) *
                  ∑ k ∈ Finset.range (n + 1),
                    rivoalRealWeight22 n k *
                      rivoalSaddleError22 n k ^ 2 by
                rw [Finset.mul_sum]
                apply Finset.sum_congr rfl
                intro k hk
                ring]
            ring
  have hterm2 :
      (∑ k ∈ Finset.range (n + 1),
        rivoalRealWeight22 n k *
          (rivoalBirth22 n k *
            (3 * (k : ℝ) ^ 2 + (k : ℝ) + 2 * (n : ℝ)))) ≤
        15 * (n : ℝ) ^ 3 * K2 + 15 * (n : ℝ) ^ 4 * Q := by
    have hp :
        (∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k *
            (rivoalBirth22 n k *
              (3 * (k : ℝ) ^ 2 + (k : ℝ) + 2 * (n : ℝ)))) ≤
          ∑ k ∈ Finset.range (n + 1),
            rivoalRealWeight22 n k *
              (15 * (n : ℝ) ^ 3 * (k : ℝ) ^ 2 +
                15 * (n : ℝ) ^ 4) := by
      apply Finset.sum_le_sum
      intro k hk
      obtain ⟨hk0, hkn⟩ := rivoalSupportBounds22 hk
      have hb := rivoalBirth22_le_five_cube hn hk
      have hdelta :
          3 * (k : ℝ) ^ 2 + (k : ℝ) + 2 * (n : ℝ) ≤
            3 * (k : ℝ) ^ 2 + 3 * (n : ℝ) := by linarith
      have hprod :
          rivoalBirth22 n k *
              (3 * (k : ℝ) ^ 2 + (k : ℝ) + 2 * (n : ℝ)) ≤
            (5 * (n : ℝ) ^ 3) *
              (3 * (k : ℝ) ^ 2 + 3 * (n : ℝ)) :=
        mul_le_mul hb hdelta (by positivity) (by positivity)
      apply mul_le_mul_of_nonneg_left
      · calc
          _ ≤ _ := hprod
          _ = 15 * (n : ℝ) ^ 3 * (k : ℝ) ^ 2 +
              15 * (n : ℝ) ^ 4 := by ring
      · exact rivoalRealWeight22_nonneg n k
    calc
      _ ≤ _ := hp
      _ = 15 * (n : ℝ) ^ 3 * K2 + 15 * (n : ℝ) ^ 4 * Q := by
        dsimp only [K2, Q]
        calc
          _ = (∑ k ∈ Finset.range (n + 1),
              (15 * (n : ℝ) ^ 3 *
                  (rivoalRealWeight22 n k * (k : ℝ) ^ 2) +
                rivoalRealWeight22 n k * (15 * (n : ℝ) ^ 4))) := by
                  apply Finset.sum_congr rfl
                  intro k hk
                  ring
          _ = 15 * (n : ℝ) ^ 3 *
                (∑ k ∈ Finset.range (n + 1),
                  rivoalRealWeight22 n k * (k : ℝ) ^ 2) +
              ∑ k ∈ Finset.range (n + 1),
                rivoalRealWeight22 n k * (15 * (n : ℝ) ^ 4) := by
                  rw [Finset.sum_add_distrib, Finset.mul_sum]
          _ = _ := by
            rw [← Finset.sum_mul, rivoalRealWeight22_sum]
            ring
  have hexact := rivoalSaddleSecondMoment22_exact n
  have hright :
      (∑ k ∈ Finset.range (n + 1),
        rivoalRealWeight22 n k *
          (2 * ((n : ℝ) - (k : ℝ)) ^ 2 * rivoalSaddleError22 n k +
            rivoalBirth22 n k *
              (3 * (k : ℝ) ^ 2 + (k : ℝ) + 2 * (n : ℝ)))) ≤
        ((n : ℝ) * S2 + (n : ℝ) ^ 3 * Q) +
          (15 * (n : ℝ) ^ 3 * K2 + 15 * (n : ℝ) ^ 4 * Q) := by
    calc
      _ = (∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k *
            (2 * ((n : ℝ) - (k : ℝ)) ^ 2 * rivoalSaddleError22 n k)) +
        ∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k *
            (rivoalBirth22 n k *
              (3 * (k : ℝ) ^ 2 + (k : ℝ) + 2 * (n : ℝ))) := by
                rw [← Finset.sum_add_distrib]
                apply Finset.sum_congr rfl
                intro k hk
                ring
      _ ≤ _ := add_le_add hterm1 hterm2
  have hcore :
      (2 * (n : ℝ)) * S2 ≤
        ((n : ℝ) * S2 + (n : ℝ) ^ 3 * Q) +
          (15 * (n : ℝ) ^ 3 * K2 + 15 * (n : ℝ) ^ 4 * Q) := by
    exact hleft.trans (hexact.trans_le hright)
  have hK2 : K2 ≤ 4 * (n : ℝ) * s * Q := by
    dsimp only [K2, s, Q]
    exact rivoalWeightedSquare22_le hn
  have hKterm :
      15 * (n : ℝ) ^ 3 * K2 ≤
        60 * (n : ℝ) ^ 4 * s * Q := by
    calc
      _ ≤ 15 * (n : ℝ) ^ 3 * (4 * (n : ℝ) * s * Q) :=
        mul_le_mul_of_nonneg_left hK2 (by positivity)
      _ = _ := by ring
  have hA : (n : ℝ) ^ 3 * Q ≤ (n : ℝ) ^ 4 * s * Q := by
    have hfac : 1 ≤ (n : ℝ) * s := by nlinarith
    have hnonneg : 0 ≤ (n : ℝ) ^ 3 * Q :=
      mul_nonneg (pow_nonneg (by positivity) 3) hQ.le
    calc
      (n : ℝ) ^ 3 * Q = ((n : ℝ) ^ 3 * Q) * 1 := by ring
      _ ≤ ((n : ℝ) ^ 3 * Q) * ((n : ℝ) * s) :=
        mul_le_mul_of_nonneg_left hfac hnonneg
      _ = (n : ℝ) ^ 4 * s * Q := by ring
  have hC : 15 * (n : ℝ) ^ 4 * Q ≤ 15 * (n : ℝ) ^ 4 * s * Q := by
    have hnonneg : 0 ≤ 15 * (n : ℝ) ^ 4 * Q := by positivity
    nlinarith
  have hsub :
      (n : ℝ) * S2 ≤
        (n : ℝ) ^ 3 * Q +
          (15 * (n : ℝ) ^ 3 * K2 + 15 * (n : ℝ) ^ 4 * Q) := by
    nlinarith [hcore]
  have h76 :
      (n : ℝ) ^ 3 * Q +
          (15 * (n : ℝ) ^ 3 * K2 + 15 * (n : ℝ) ^ 4 * Q) ≤
        76 * (n : ℝ) ^ 4 * s * Q := by
    calc
      _ ≤ (n : ℝ) ^ 4 * s * Q +
          (60 * (n : ℝ) ^ 4 * s * Q +
            15 * (n : ℝ) ^ 4 * s * Q) :=
        add_le_add hA (add_le_add hKterm hC)
      _ = 76 * (n : ℝ) ^ 4 * s * Q := by ring
  have h100 :
      76 * (n : ℝ) ^ 4 * s * Q ≤
        100 * (n : ℝ) ^ 4 * s * Q := by
    have hnonneg : 0 ≤ (n : ℝ) ^ 4 * s * Q := by positivity
    nlinarith
  have hmul :
      (n : ℝ) * S2 ≤ 100 * (n : ℝ) ^ 4 * s * Q :=
    hsub.trans (h76.trans h100)
  apply le_of_mul_le_mul_of_pos_left _ hnR
  dsimp only [S2, s, Q] at hmul ⊢
  convert hmul using 1 <;> ring

/-- Weighted Cauchy--Schwarz for the saddle observable. -/
theorem rivoalWeightedCauchySchwarz22 (n : ℕ) :
    (∑ k ∈ Finset.range (n + 1),
        rivoalRealWeight22 n k * |rivoalSaddleError22 n k|) ^ 2 ≤
      ((rivoalExplicitQ22 n : ℚ) : ℝ) *
        (∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k * rivoalSaddleError22 n k ^ 2) := by
  have hcs := Finset.sq_sum_div_le_sum_sq_div
    (Finset.range (n + 1))
    (fun k => rivoalRealWeight22 n k * |rivoalSaddleError22 n k|)
    (g := fun k => rivoalRealWeight22 n k)
    (fun k hk => rivoalRealWeight22_pos n k
      (Nat.lt_succ_iff.mp (Finset.mem_range.mp hk)))
  have hright :
      (∑ k ∈ Finset.range (n + 1),
        (rivoalRealWeight22 n k * |rivoalSaddleError22 n k|) ^ 2 /
          rivoalRealWeight22 n k) =
        ∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k * rivoalSaddleError22 n k ^ 2 := by
    apply Finset.sum_congr rfl
    intro k hk
    have hw := rivoalRealWeight22_pos n k
      (Nat.lt_succ_iff.mp (Finset.mem_range.mp hk))
    field_simp [hw.ne']
  rw [rivoalRealWeight22_sum, hright] at hcs
  have hQ := rivoalRealQ22_pos n
  have hm := (div_le_iff₀ hQ).1 hcs
  simpa [mul_comm] using hm

/-- Squared normalized saddle mean. -/
theorem rivoalSaddleMean22_sq_le {n : ℕ} (hn : 1 ≤ n) :
    (((∑ k ∈ Finset.range (n + 1),
        rivoalRealWeight22 n k * |rivoalSaddleError22 n k|) /
      (((rivoalExplicitQ22 n : ℚ) : ℝ) * (n : ℝ) ^ 2)) ^ 2) ≤
      100 / Real.sqrt (n : ℝ) := by
  let Q : ℝ := ((rivoalExplicitQ22 n : ℚ) : ℝ)
  let M : ℝ := ∑ k ∈ Finset.range (n + 1),
    rivoalRealWeight22 n k * |rivoalSaddleError22 n k|
  let V : ℝ := ∑ k ∈ Finset.range (n + 1),
    rivoalRealWeight22 n k * rivoalSaddleError22 n k ^ 2
  let s : ℝ := Real.sqrt (n : ℝ)
  have hnR : (0 : ℝ) < (n : ℝ) := by
    exact_mod_cast (lt_of_lt_of_le Nat.zero_lt_one hn)
  have hQ : 0 < Q := rivoalRealQ22_pos n
  have hs : 0 < s := Real.sqrt_pos.2 hnR
  have hs2 : s ^ 2 = (n : ℝ) := Real.sq_sqrt hnR.le
  have hcs : M ^ 2 ≤ Q * V := by
    dsimp only [M, Q, V]
    exact rivoalWeightedCauchySchwarz22 n
  have hV : V ≤ 100 * (n : ℝ) ^ 3 * s * Q := by
    dsimp only [V, s, Q]
    exact rivoalSaddleSecondMoment22_le hn
  have hnum : M ^ 2 ≤ 100 * (n : ℝ) ^ 3 * s * Q ^ 2 := by
    calc
      M ^ 2 ≤ Q * V := hcs
      _ ≤ Q * (100 * (n : ℝ) ^ 3 * s * Q) :=
        mul_le_mul_of_nonneg_left hV hQ.le
      _ = 100 * (n : ℝ) ^ 3 * s * Q ^ 2 := by ring
  rw [div_pow]
  apply (div_le_iff₀ (sq_pos_of_pos (mul_pos hQ (sq_pos_of_pos hnR)))).2
  calc
    M ^ 2 ≤ 100 * (n : ℝ) ^ 3 * s * Q ^ 2 := hnum
    _ = (100 / s) * (Q * (n : ℝ) ^ 2) ^ 2 := by
      rw [← hs2]
      field_simp [hs.ne']

private theorem rivoalSaddleMean22_nonneg (n : ℕ) :
    0 ≤
      (∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k * |rivoalSaddleError22 n k|) /
        (((rivoalExplicitQ22 n : ℚ) : ℝ) * (n : ℝ) ^ 2) := by
  apply div_nonneg
  · apply Finset.sum_nonneg
    intro k hk
    exact mul_nonneg (rivoalRealWeight22_nonneg n k) (abs_nonneg _)
  · exact mul_nonneg (rivoalRealQ22_pos n).le (sq_nonneg _)

/-- The normalized weighted absolute saddle error tends to zero. -/
theorem tendsto_saddle_mean22 :
    Tendsto
      (fun n : ℕ =>
        (∑ k ∈ Finset.range (n + 1),
            rivoalRealWeight22 n k * |rivoalSaddleError22 n k|) /
          (((rivoalExplicitQ22 n : ℚ) : ℝ) * (n : ℝ) ^ 2))
      atTop (𝓝 0) := by
  have hsqrtTop :
      Tendsto (fun n : ℕ => Real.sqrt (n : ℝ)) atTop atTop :=
    Real.tendsto_sqrt_atTop.comp tendsto_natCast_atTop_atTop
  have hinv :
      Tendsto (fun n : ℕ => (Real.sqrt (n : ℝ))⁻¹) atTop (𝓝 0) :=
    tendsto_inv_atTop_zero.comp hsqrtTop
  have hinner :
      Tendsto (fun n : ℕ => 100 / Real.sqrt (n : ℝ)) atTop (𝓝 0) := by
    simpa only [div_eq_mul_inv, mul_zero] using hinv.const_mul 100
  have hroot :
      Tendsto (fun n : ℕ => Real.sqrt (100 / Real.sqrt (n : ℝ)))
        atTop (𝓝 0) := by
    simpa using (Real.continuous_sqrt.tendsto 0).comp hinner
  apply squeeze_zero'
  · exact Eventually.of_forall rivoalSaddleMean22_nonneg
  · filter_upwards [eventually_atTop.2 ⟨1, fun n hn => hn⟩] with n hn
    have hsquare := rivoalSaddleMean22_sq_le hn
    have hmean := rivoalSaddleMean22_nonneg n
    have hright : 0 ≤ 100 / Real.sqrt (n : ℝ) := by positivity
    have hsqrt0 := Real.sqrt_nonneg (100 / Real.sqrt (n : ℝ))
    have hsqrt2 := Real.sq_sqrt hright
    change
      (∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k * |rivoalSaddleError22 n k|) /
        (((rivoalExplicitQ22 n : ℚ) : ℝ) * (n : ℝ) ^ 2) ≤
      Real.sqrt (100 / Real.sqrt (n : ℝ))
    nlinarith
  · exact hroot

#print axioms rivoalWeightedCube22_le
#print axioms rivoalWeightedSquare22_le
#print axioms rivoalSaddleSecondMoment22_exact
#print axioms rivoalSaddleSecondMoment22_le
#print axioms rivoalWeightedCauchySchwarz22
#print axioms tendsto_saddle_mean22

end RamanujanChallenge.P22

end
