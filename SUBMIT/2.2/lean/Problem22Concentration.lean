import RamanujanChallenge.Problem22
import Mathlib.NumberTheory.Harmonic.EulerMascheroni

/-!
# Problem 2.2: explicit concentration infrastructure

This module starts the remaining analytic layer of Problem 2.2.  The exact
positive weights are those already identified in `Problem22`; the first
results below expose their adjacent ratio in a form suitable for elementary
geometric tail estimates.
-/

noncomputable section

open Filter Topology Real
open scoped BigOperators

namespace RamanujanChallenge.P22

/-- The exact adjacent multiplier for the positive Rivoal weights. -/
def rivoalWeightRatio22 (n k : ℕ) : ℝ :=
  (2 * (n : ℝ) + (k : ℝ) + 2) /
      (2 * (n : ℝ) + (k : ℝ) + 1) *
    (((n : ℝ) - (k : ℝ)) ^ 2 / ((k : ℝ) + 1) ^ 3)

/-- Exact adjacent ratio for the rational hypergeometric weight, including
the linear factor occurring in Rivoal's denominator sum. -/
theorem rivoalExplicitQTerm22_succ_eq (n k : ℕ) :
    rivoalExplicitQTerm22 n (k + 1) =
      rivoalExplicitQTerm22 n k *
        (2 * (n : ℚ) + (k : ℚ) + 2) /
          (2 * (n : ℚ) + (k : ℚ) + 1) *
        (((n : ℚ) - (k : ℚ)) ^ 2 / ((k : ℚ) + 1) ^ 3) := by
  rw [rivoalExplicitQTerm22, rivoalExplicitQTerm22,
    rivoalWeight22_succ]
  norm_num [Nat.cast_add, Nat.cast_one]
  have hk : (k : ℚ) + 1 ≠ 0 := by positivity
  field_simp [hk]
  ring

/-- Real form of the exact adjacent ratio. -/
theorem rivoalRealWeight22_succ_eq (n k : ℕ) :
    rivoalRealWeight22 n (k + 1) =
      rivoalRealWeight22 n k * rivoalWeightRatio22 n k := by
  have h := congrArg (fun q : ℚ => (q : ℝ))
    (rivoalExplicitQTerm22_succ_eq n k)
  norm_num [rivoalRealWeight22, rivoalWeightRatio22,
    Nat.cast_add, Nat.cast_one] at h ⊢
  rw [h]
  ring

/-- Cross-multiplied adjacent-ratio identity.  This is often more convenient
than the quotient form for `nlinarith` and positivity arguments. -/
theorem rivoalRealWeight22_succ_cross (n k : ℕ) :
    rivoalRealWeight22 n (k + 1) *
          (2 * (n : ℝ) + (k : ℝ) + 1) *
          ((k : ℝ) + 1) ^ 3 =
      rivoalRealWeight22 n k *
          (2 * (n : ℝ) + (k : ℝ) + 2) *
          ((n : ℝ) - (k : ℝ)) ^ 2 := by
  rw [rivoalRealWeight22_succ_eq]
  have hlin : 2 * (n : ℝ) + (k : ℝ) + 1 ≠ 0 := by positivity
  have hk : (k : ℝ) + 1 ≠ 0 := by positivity
  rw [rivoalWeightRatio22]
  field_simp [hlin, hk]

/-- On the support of the finite sum the real weights are strictly positive. -/
theorem rivoalRealWeight22_pos (n k : ℕ) (hk : k ≤ n) :
    0 < rivoalRealWeight22 n k := by
  simp only [rivoalRealWeight22, rivoalExplicitQTerm22,
    rivoalWeight22_eq_choose, Rat.cast_mul, Rat.cast_add,
    Rat.cast_natCast, Rat.cast_one, Rat.cast_div, Rat.cast_pow]
  have hchoose : 0 < n.choose k := Nat.choose_pos hk
  positivity

/-- The harmless linear factor makes the exact ratio slightly larger than
the pure hypergeometric ratio. -/
theorem one_lt_rivoalLinearRatio22 (n k : ℕ) :
    1 <
      (2 * (n : ℝ) + (k : ℝ) + 2) /
        (2 * (n : ℝ) + (k : ℝ) + 1) := by
  apply (lt_div_iff₀ (by positivity :
    0 < 2 * (n : ℝ) + (k : ℝ) + 1)).2
  linarith

/-- A uniform upper bound for the same linear correction. -/
theorem rivoalLinearRatio22_le_two (n k : ℕ) :
    (2 * (n : ℝ) + (k : ℝ) + 2) /
        (2 * (n : ℝ) + (k : ℝ) + 1) ≤ 2 := by
  apply (div_le_iff₀ (by positivity :
    0 < 2 * (n : ℝ) + (k : ℝ) + 1)).2
  have hn : (0 : ℝ) ≤ n := by positivity
  have hk : (0 : ℝ) ≤ k := by positivity
  linarith

/-- A cross-multiplied estimate implies a lower bound for the adjacent
weight multiplier. -/
theorem le_rivoalWeightRatio22_of_cross
    (n k : ℕ) (q : ℝ)
    (hcross :
      q * (2 * (n : ℝ) + (k : ℝ) + 1) * ((k : ℝ) + 1) ^ 3 ≤
        (2 * (n : ℝ) + (k : ℝ) + 2) *
          ((n : ℝ) - (k : ℝ)) ^ 2) :
    q ≤ rivoalWeightRatio22 n k := by
  have hlin : 0 < 2 * (n : ℝ) + (k : ℝ) + 1 := by positivity
  have hk : 0 < ((k : ℝ) + 1) ^ 3 := by positivity
  have hform :
      rivoalWeightRatio22 n k =
        ((2 * (n : ℝ) + (k : ℝ) + 2) *
          ((n : ℝ) - (k : ℝ)) ^ 2) /
          ((2 * (n : ℝ) + (k : ℝ) + 1) * ((k : ℝ) + 1) ^ 3) := by
    rw [rivoalWeightRatio22]
    field_simp [hlin.ne', hk.ne']
  rw [hform]
  apply (le_div_iff₀ (mul_pos hlin hk)).2
  nlinarith

/-- A cross-multiplied estimate implies an upper bound for the adjacent
weight multiplier. -/
theorem rivoalWeightRatio22_le_of_cross
    (n k : ℕ) (q : ℝ)
    (hcross :
      (2 * (n : ℝ) + (k : ℝ) + 2) *
          ((n : ℝ) - (k : ℝ)) ^ 2 ≤
        q * (2 * (n : ℝ) + (k : ℝ) + 1) * ((k : ℝ) + 1) ^ 3) :
    rivoalWeightRatio22 n k ≤ q := by
  have hlin : 0 < 2 * (n : ℝ) + (k : ℝ) + 1 := by positivity
  have hk : 0 < ((k : ℝ) + 1) ^ 3 := by positivity
  have hform :
      rivoalWeightRatio22 n k =
        ((2 * (n : ℝ) + (k : ℝ) + 2) *
          ((n : ℝ) - (k : ℝ)) ^ 2) /
          ((2 * (n : ℝ) + (k : ℝ) + 1) * ((k : ℝ) + 1) ^ 3) := by
    rw [rivoalWeightRatio22]
    field_simp [hlin.ne', hk.ne']
  rw [hform]
  apply (div_le_iff₀ (mul_pos hlin hk)).2
  nlinarith

/-- Geometric control of a finite left tail from adjacent-step bounds. -/
theorem sum_range_le_geometric_last22
    (w : ℕ → ℝ) (q : ℝ)
    (hq0 : 0 ≤ q) (hq1 : q < 1)
    (hw : ∀ k, 0 ≤ w k)
    (hstep : ∀ k, w k ≤ q * w (k + 1)) :
    ∀ K : ℕ,
      (∑ k ∈ Finset.range K, w k) ≤ q / (1 - q) * w K := by
  intro K
  induction K with
  | zero =>
      simp
      exact mul_nonneg (div_nonneg hq0 (by linarith)) (hw 0)
  | succ K ih =>
      rw [Finset.sum_range_succ]
      have hden : 0 < 1 - q := by linarith
      have hscale : 0 ≤ 1 / (1 - q) := by positivity
      have hstep' := mul_le_mul_of_nonneg_left (hstep K) hscale
      calc
        (∑ k ∈ Finset.range K, w k) + w K ≤
            q / (1 - q) * w K + w K := by
              simpa [add_comm] using add_le_add_right ih (w K)
        _ = (1 / (1 - q)) * w K := by field_simp; ring
        _ ≤ (1 / (1 - q)) * (q * w (K + 1)) := hstep'
        _ = q / (1 - q) * w (K + 1) := by ring

/-! ## A finite Stein identity for the same weights -/

/-- The birth polynomial in the cross-multiplied adjacent-weight identity. -/
def rivoalBirth22 (n k : ℕ) : ℝ :=
  (2 * (n : ℝ) + (k : ℝ) + 2) * ((n : ℝ) - (k : ℝ)) ^ 2

/-- The death polynomial in the cross-multiplied adjacent-weight identity. -/
def rivoalDeath22 (n k : ℕ) : ℝ :=
  (2 * (n : ℝ) + (k : ℝ)) * (k : ℝ) ^ 3

/-- The polynomial error whose zero set is the cubic saddle of the weights. -/
def rivoalSaddleError22 (n k : ℕ) : ℝ :=
  (k : ℝ) ^ 3 - ((n : ℝ) - (k : ℝ)) ^ 2

/-- Reindexing form of the finite birth--death identity.  The two missing
boundary terms vanish because `rivoalBirth22 n n = 0` and
`rivoalDeath22 n 0 = 0`. -/
theorem rivoalWeightSteinShift22 (n : ℕ) (f : ℕ → ℝ) :
    (∑ k ∈ Finset.range (n + 1),
        rivoalRealWeight22 n k * rivoalBirth22 n k * f (k + 1)) =
      ∑ k ∈ Finset.range (n + 1),
        rivoalRealWeight22 n k * rivoalDeath22 n k * f k := by
  rw [Finset.sum_range_succ, Finset.sum_range_succ']
  have hbirth :
      rivoalRealWeight22 n n * rivoalBirth22 n n * f (n + 1) = 0 := by
    simp [rivoalBirth22]
  have hdeath :
      rivoalRealWeight22 n 0 * rivoalDeath22 n 0 * f 0 = 0 := by
    simp [rivoalDeath22]
  simp only [hbirth, hdeath, add_zero]
  apply Finset.sum_congr rfl
  intro k hk
  have hcross := rivoalRealWeight22_succ_cross n k
  calc
    rivoalRealWeight22 n k * rivoalBirth22 n k * f (k + 1) =
        (rivoalRealWeight22 n k *
          (2 * (n : ℝ) + (k : ℝ) + 2) *
          ((n : ℝ) - (k : ℝ)) ^ 2) * f (k + 1) := by
            rw [rivoalBirth22]
            ring
    _ =
        (rivoalRealWeight22 n (k + 1) *
          (2 * (n : ℝ) + (k : ℝ) + 1) *
          ((k : ℝ) + 1) ^ 3) * f (k + 1) := by
            rw [hcross]
    _ = rivoalRealWeight22 n (k + 1) *
          rivoalDeath22 n (k + 1) * f (k + 1) := by
            have hdeathpoly :
                (2 * (n : ℝ) + (k : ℝ) + 1) * ((k : ℝ) + 1) ^ 3 =
                  rivoalDeath22 n (k + 1) := by
              simp only [rivoalDeath22, Nat.cast_add, Nat.cast_one]
              ring
            calc
              rivoalRealWeight22 n (k + 1) *
                    (2 * (n : ℝ) + (k : ℝ) + 1) *
                    ((k : ℝ) + 1) ^ 3 * f (k + 1) =
                  rivoalRealWeight22 n (k + 1) *
                    ((2 * (n : ℝ) + (k : ℝ) + 1) *
                      ((k : ℝ) + 1) ^ 3) * f (k + 1) := by ring
              _ = rivoalRealWeight22 n (k + 1) *
                    rivoalDeath22 n (k + 1) * f (k + 1) := by
                      rw [hdeathpoly]

/-- Finite Stein identity for the positive Rivoal weights. -/
theorem rivoalWeightStein22 (n : ℕ) (f : ℕ → ℝ) :
    (∑ k ∈ Finset.range (n + 1),
      rivoalRealWeight22 n k *
        (rivoalBirth22 n k * (f (k + 1) - f k) +
          (rivoalBirth22 n k - rivoalDeath22 n k) * f k)) = 0 := by
  calc
    (∑ k ∈ Finset.range (n + 1),
      rivoalRealWeight22 n k *
        (rivoalBirth22 n k * (f (k + 1) - f k) +
          (rivoalBirth22 n k - rivoalDeath22 n k) * f k)) =
        (∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k * rivoalBirth22 n k * f (k + 1)) -
        (∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k * rivoalDeath22 n k * f k) := by
            rw [← Finset.sum_sub_distrib]
            apply Finset.sum_congr rfl
            intro k hk
            ring
    _ = 0 := sub_eq_zero.mpr (rivoalWeightSteinShift22 n f)

end RamanujanChallenge.P22

end
