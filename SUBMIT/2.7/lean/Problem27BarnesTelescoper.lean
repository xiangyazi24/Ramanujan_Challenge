import RamanujanChallenge.Problem27BarnesNormalization
import Mathlib.Analysis.Calculus.Deriv.Shift

open Filter Set Topology
open scoped BigOperators

noncomputable section

namespace RamanujanChallenge.P27

/-! ## Homogeneous source-recurrence coefficients -/

def ctAlpha27 (k : ℕ) : ℚ :=
  2 * zudilinA (k + 2)

def ctBeta27 (k : ℕ) : ℚ :=
  2 * zudilinM (k + 2)

def ctGamma27 (k : ℕ) : ℚ :=
  2 * (k + 2) * zudilinN (k + 2)

def ctDelta27 (k : ℕ) : ℚ :=
  zudilinR (k + 2) * (k + 2) * (k + 1) ^ 3

/-- Exact match with the existing `ZudilinSatisfiesRec` convention. -/
theorem zudilinSatisfiesRec_iff_ct27 (v : ℕ → ℚ) :
    ZudilinSatisfiesRec v ↔
      ∀ k : ℕ,
        ctAlpha27 k * v (k + 3)
          - ctBeta27 k * v (k + 2)
          + ctGamma27 k * v (k + 1)
          - ctDelta27 k * v k = 0 := by
  constructor
  · intro hv k
    rw [hv k]
    simp only [ctAlpha27, ctBeta27, ctGamma27, ctDelta27]
    field_simp [ne_of_gt (zudilinA_pos k)]
    ring
  · intro h k
    have hk := h k
    simp only [ctAlpha27, ctBeta27, ctGamma27, ctDelta27] at hk
    have hA : zudilinA (k + 2) ≠ 0 := ne_of_gt (zudilinA_pos k)
    field_simp [hA]
    linear_combination hk

/-! ## Universal complex coefficient polynomials -/

def ctAlphaC27 (N : ℂ) : ℂ :=
  2 * (946 * N ^ 2 - 731 * N + 153) *
    (2 * N + 1) * (N + 1) ^ 3

def ctAlphaBarC27 (N : ℂ) : ℂ :=
  2 * (946 * N ^ 2 - 731 * N + 153) *
    (2 * N + 1) * (N + 1)

def ctBetaC27 (N : ℂ) : ℂ :=
  2 * (104060 * N ^ 6 + 127710 * N ^ 5 + 12788 * N ^ 4
    - 34525 * N ^ 3 - 8482 * N ^ 2 + 3298 * N + 1071)

def ctGammaC27 (N : ℂ) : ℂ :=
  2 * N * (3784 * N ^ 5 - 1032 * N ^ 4 - 1925 * N ^ 3
    + 853 * N ^ 2 + 328 * N - 184)

def ctDeltaC27 (N : ℂ) : ℂ :=
  (946 * N ^ 2 + 1161 * N + 368) * N * (N - 1) ^ 3

def ctQ0C27 (N : ℂ) : ℂ :=
  -4 * (N - 1) ^ 2 * N ^ 4 * (2 * N - 1) *
    (946 * N ^ 2 + 1161 * N + 368)

def ctQ1C27 (N : ℂ) : ℂ :=
  2 * (N - 1) * N ^ 3 * (10 * N ^ 2 - 2 * N - 1) *
    (946 * N ^ 2 + 1161 * N + 368)

def ctQ2C27 (N : ℂ) : ℂ :=
  -N ^ 3 * (38 * N ^ 2 - N - 5) *
    (946 * N ^ 2 + 1161 * N + 368)

def ctQ3C27 (N : ℂ) : ℂ :=
  2448 + 5440 * N - 26056 * N ^ 2 - 43516 * N ^ 3
    + 69621 * N ^ 4 + 156993 * N ^ 5 + 74734 * N ^ 6

def ctQ4C27 (N : ℂ) : ℂ :=
  -2 * (612 + 2971 * N - 12132 * N ^ 2 - 14629 * N ^ 3
    + 41538 * N ^ 4 + 45408 * N ^ 5)

def ctQ5C27 (N : ℂ) : ℂ :=
  -3 * (2 * N ^ 2 + 6 * N + 3) *
    (946 * N ^ 2 - 731 * N + 153)

def ctQ6C27 (N : ℂ) : ℂ :=
  (2 * N + 1) * (946 * N ^ 2 - 731 * N + 153)

def ctQhatC27 (N x : ℂ) : ℂ :=
  ctQ0C27 N + ctQ1C27 N * x + ctQ2C27 N * x ^ 2
    + ctQ3C27 N * x ^ 3 + ctQ4C27 N * x ^ 4
    + ctQ5C27 N * x ^ 5 + ctQ6C27 N * x ^ 6

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem ctQhat_poly_identity27 (N x : ℂ) :
    (N + x) ^ 4 * ctQhatC27 N (x + 1)
      - (x + 1) ^ 3 * (x + 2 * N + 1) * ctQhatC27 N x =
    ctAlphaBarC27 N * (x - 1) ^ 3 * x ^ 3 * (x + 1) ^ 3
      - ctBetaC27 N * (x + 1) ^ 3 * x ^ 3 * (x + 2 * N + 1)
      + ctGammaC27 N * N ^ 2 * (x + 2 * N) *
          (x + 1) ^ 3 * (x + 2 * N + 1)
      - ctDeltaC27 N * N ^ 2 * (N - 1) ^ 2 *
          (x + 2 * N) * (x + 2 * N - 1) *
          (x + 2 * N + 1) := by
  simp only [ctQhatC27, ctQ0C27, ctQ1C27, ctQ2C27, ctQ3C27,
    ctQ4C27, ctQ5C27, ctQ6C27, ctAlphaBarC27, ctBetaC27,
    ctGammaC27, ctDeltaC27]
  ring

theorem ctAlphaC_eq27 (N : ℂ) :
    ctAlphaC27 N = (N + 1) ^ 2 * ctAlphaBarC27 N := by
  simp only [ctAlphaC27, ctAlphaBarC27]
  ring

@[simp] theorem ctAlphaC_nat27 (k : ℕ) :
    ctAlphaC27 (((k + 2 : ℕ) : ℂ)) = (ctAlpha27 k : ℂ) := by
  simp only [ctAlphaC27, ctAlpha27, zudilinA, zudilinQ]
  push_cast
  ring

@[simp] theorem ctBetaC_nat27 (k : ℕ) :
    ctBetaC27 (((k + 2 : ℕ) : ℂ)) = (ctBeta27 k : ℂ) := by
  simp only [ctBetaC27, ctBeta27, zudilinM]
  push_cast
  ring

@[simp] theorem ctGammaC_nat27 (k : ℕ) :
    ctGammaC27 (((k + 2 : ℕ) : ℂ)) = (ctGamma27 k : ℂ) := by
  simp only [ctGammaC27, ctGamma27, zudilinN]
  push_cast
  ring

@[simp] theorem ctDeltaC_nat27 (k : ℕ) :
    ctDeltaC27 (((k + 2 : ℕ) : ℂ)) = (ctDelta27 k : ℂ) := by
  simp only [ctDeltaC27, ctDelta27, zudilinR]
  push_cast
  ring

/-! ## Direct `t`-coordinate rational functions -/

def ctNumerator27 (n : ℕ) (t : ℂ) : ℂ :=
  ∏ r ∈ Finset.range n, (t - (((r + 1 : ℕ) : ℂ))) ^ 3

def ctPoleProduct27 (m : ℕ) (t : ℂ) : ℂ :=
  ∏ j ∈ Finset.range m, (t + (j : ℂ))

def ctR27 (n : ℕ) (t : ℂ) : ℂ :=
  ctNumerator27 n t /
    (((n.factorial : ℂ) ^ 2) * ctPoleProduct27 (n + 1) t)

/-- The `(t-n)^3` factor is cancelled in the definition itself. -/
def ctS27 (n : ℕ) (t : ℂ) : ℂ :=
  ctQhatC27 (n : ℂ) (t - (n : ℂ)) * ctNumerator27 (n - 1) t /
    (((n.factorial : ℂ) ^ 2) * ctPoleProduct27 (n + 1) t)

@[simp] theorem ctNumerator_succ27 (n : ℕ) (t : ℂ) :
    ctNumerator27 (n + 1) t =
      ctNumerator27 n t * (t - (((n + 1 : ℕ) : ℂ))) ^ 3 := by
  simp only [ctNumerator27, Finset.prod_range_succ]

@[simp] theorem ctPoleProduct_succ27 (m : ℕ) (t : ℂ) :
    ctPoleProduct27 (m + 1) t =
      ctPoleProduct27 m t * (t + (m : ℂ)) := by
  simp only [ctPoleProduct27, Finset.prod_range_succ]

theorem ctNumerator_shift_succ27 (m : ℕ) (t : ℂ) :
    ctNumerator27 (m + 1) (t + 1) =
      t ^ 3 * ctNumerator27 m t := by
  induction m with
  | zero =>
      simp [ctNumerator27]
  | succ m ih =>
      rw [ctNumerator_succ27 (m + 1) (t + 1), ih,
        ctNumerator_succ27 m t]
      push_cast
      ring

theorem ctPoleProduct_shift_succ27 (m : ℕ) (t : ℂ) :
    ctPoleProduct27 (m + 1) t =
      t * ctPoleProduct27 m (t + 1) := by
  induction m with
  | zero =>
      simp [ctPoleProduct27]
  | succ m ih =>
      rw [ctPoleProduct_succ27 (m + 1) t,
        ctPoleProduct_succ27 m (t + 1), ih]
      push_cast
      ring

theorem ctPoleProduct_ne_zero_of_le27
    {m M : ℕ} {t : ℂ} (hm : m ≤ M)
    (hM : ctPoleProduct27 M t ≠ 0) :
    ctPoleProduct27 m t ≠ 0 := by
  unfold ctPoleProduct27 at hM ⊢
  rw [Finset.prod_ne_zero_iff] at hM ⊢
  intro j hj
  exact hM j (Finset.mem_range.mpr
    ((Finset.mem_range.mp hj).trans_le hm))

@[simp] theorem ct_factorial_succ_cast27 (n : ℕ) :
    (((n + 1).factorial : ℂ)) =
      (((n + 1 : ℕ) : ℂ)) * (n.factorial : ℂ) := by
  exact_mod_cast Nat.factorial_succ n

/-! ## Six common-denominator normal forms -/

def ctCommonDen27 (k : ℕ) (t : ℂ) : ℂ :=
  (((k + 2).factorial : ℂ) ^ 2) * ctPoleProduct27 (k + 4) t

private theorem ctDen_ne_zero27
    {n m : ℕ} {t : ℂ} (hP : ctPoleProduct27 m t ≠ 0) :
    ((n.factorial : ℂ) ^ 2) * ctPoleProduct27 m t ≠ 0 :=
  mul_ne_zero
    (pow_ne_zero _ (Nat.cast_ne_zero.mpr (Nat.factorial_ne_zero n))) hP

private theorem ct_common_mul_alpha_R_next27
    (k : ℕ) (t : ℂ)
    (hP : ctPoleProduct27 (k + 4) t ≠ 0) :
    ctCommonDen27 k t *
        (ctAlphaC27 (((k + 2 : ℕ) : ℂ)) * ctR27 (k + 3) t) =
      ctAlphaBarC27 (((k + 2 : ℕ) : ℂ)) * ctNumerator27 k t *
        (t - ((k + 2 : ℕ) : ℂ) - 1) ^ 3 *
        (t - ((k + 2 : ℕ) : ℂ)) ^ 3 *
        (t - ((k + 2 : ℕ) : ℂ) + 1) ^ 3 := by
  have hden :
      (((k + 3).factorial : ℂ) ^ 2) *
          ctPoleProduct27 (k + 4) t ≠ 0 :=
    ctDen_ne_zero27 (n := k + 3) hP
  have hfac : ((k + 3).factorial : ℂ) ≠ 0 :=
    Nat.cast_ne_zero.mpr (Nat.factorial_ne_zero (k + 3))
  unfold ctCommonDen27 ctR27
  rw [ctNumerator_succ27 (k + 2) t,
    ctNumerator_succ27 (k + 1) t,
    ctNumerator_succ27 k t, ctAlphaC_eq27]
  field_simp [hden, hfac]
  rw [ct_factorial_succ_cast27 (k + 2)]
  push_cast
  ring

private theorem ct_common_mul_beta_R27
    (k : ℕ) (t : ℂ)
    (hP : ctPoleProduct27 (k + 4) t ≠ 0) :
    ctCommonDen27 k t *
        (ctBetaC27 (((k + 2 : ℕ) : ℂ)) * ctR27 (k + 2) t) =
      ctBetaC27 (((k + 2 : ℕ) : ℂ)) * ctNumerator27 k t *
        (t - ((k + 2 : ℕ) : ℂ) + 1) ^ 3 *
        (t - ((k + 2 : ℕ) : ℂ)) ^ 3 *
        (t + ((k + 2 : ℕ) : ℂ) + 1) := by
  have hsub : ctPoleProduct27 (k + 3) t ≠ 0 :=
    ctPoleProduct_ne_zero_of_le27
      (m := k + 3) (M := k + 4) (t := t) (by omega) hP
  have hden :
      (((k + 2).factorial : ℂ) ^ 2) *
          ctPoleProduct27 (k + 3) t ≠ 0 :=
    ctDen_ne_zero27 (n := k + 2) hsub
  have hfac : ((k + 2).factorial : ℂ) ≠ 0 :=
    Nat.cast_ne_zero.mpr (Nat.factorial_ne_zero (k + 2))
  unfold ctCommonDen27 ctR27
  rw [ctNumerator_succ27 (k + 1) t,
    ctNumerator_succ27 k t, ctPoleProduct_succ27 (k + 3) t]
  field_simp [hden, hfac]
  push_cast
  ring

private theorem ct_common_mul_gamma_R_prev27
    (k : ℕ) (t : ℂ)
    (hP : ctPoleProduct27 (k + 4) t ≠ 0) :
    ctCommonDen27 k t *
        (ctGammaC27 (((k + 2 : ℕ) : ℂ)) * ctR27 (k + 1) t) =
      ctGammaC27 (((k + 2 : ℕ) : ℂ)) *
        (((k + 2 : ℕ) : ℂ) ^ 2) * ctNumerator27 k t *
        (t + ((k + 2 : ℕ) : ℂ)) *
        (t - ((k + 2 : ℕ) : ℂ) + 1) ^ 3 *
        (t + ((k + 2 : ℕ) : ℂ) + 1) := by
  have hsub : ctPoleProduct27 (k + 2) t ≠ 0 :=
    ctPoleProduct_ne_zero_of_le27
      (m := k + 2) (M := k + 4) (t := t) (by omega) hP
  have hden :
      (((k + 1).factorial : ℂ) ^ 2) *
          ctPoleProduct27 (k + 2) t ≠ 0 :=
    ctDen_ne_zero27 (n := k + 1) hsub
  have hfac : ((k + 1).factorial : ℂ) ≠ 0 :=
    Nat.cast_ne_zero.mpr (Nat.factorial_ne_zero (k + 1))
  unfold ctCommonDen27 ctR27
  rw [ctNumerator_succ27 k t,
    ctPoleProduct_succ27 (k + 3) t,
    ctPoleProduct_succ27 (k + 2) t]
  field_simp [hden, hfac]
  rw [ct_factorial_succ_cast27 (k + 1)]
  push_cast
  ring

private theorem ct_common_mul_delta_R_prev2_27
    (k : ℕ) (t : ℂ)
    (hP : ctPoleProduct27 (k + 4) t ≠ 0) :
    ctCommonDen27 k t *
        (ctDeltaC27 (((k + 2 : ℕ) : ℂ)) * ctR27 k t) =
      ctDeltaC27 (((k + 2 : ℕ) : ℂ)) *
        (((k + 2 : ℕ) : ℂ) ^ 2) *
        ((((k + 2 : ℕ) : ℂ) - 1) ^ 2) * ctNumerator27 k t *
        (t + ((k + 2 : ℕ) : ℂ)) *
        (t + ((k + 2 : ℕ) : ℂ) - 1) *
        (t + ((k + 2 : ℕ) : ℂ) + 1) := by
  have hsub : ctPoleProduct27 (k + 1) t ≠ 0 :=
    ctPoleProduct_ne_zero_of_le27
      (m := k + 1) (M := k + 4) (t := t) (by omega) hP
  have hden :
      ((k.factorial : ℂ) ^ 2) * ctPoleProduct27 (k + 1) t ≠ 0 :=
    ctDen_ne_zero27 (n := k) hsub
  have hfac : (k.factorial : ℂ) ≠ 0 :=
    Nat.cast_ne_zero.mpr (Nat.factorial_ne_zero k)
  unfold ctCommonDen27 ctR27
  rw [ctPoleProduct_succ27 (k + 3) t,
    ctPoleProduct_succ27 (k + 2) t,
    ctPoleProduct_succ27 (k + 1) t]
  field_simp [hden, hfac]
  rw [ct_factorial_succ_cast27 (k + 1),
    ct_factorial_succ_cast27 k]
  push_cast
  ring

private theorem ct_common_mul_S27
    (k : ℕ) (t : ℂ)
    (hP : ctPoleProduct27 (k + 4) t ≠ 0) :
    ctCommonDen27 k t * ctS27 (k + 2) t =
      ctNumerator27 k t *
        (t - ((k + 2 : ℕ) : ℂ) + 1) ^ 3 *
        (t + ((k + 2 : ℕ) : ℂ) + 1) *
        ctQhatC27 (((k + 2 : ℕ) : ℂ))
          (t - ((k + 2 : ℕ) : ℂ)) := by
  have hsub : ctPoleProduct27 (k + 3) t ≠ 0 :=
    ctPoleProduct_ne_zero_of_le27
      (m := k + 3) (M := k + 4) (t := t) (by omega) hP
  have hden :
      (((k + 2).factorial : ℂ) ^ 2) *
          ctPoleProduct27 (k + 3) t ≠ 0 :=
    ctDen_ne_zero27 (n := k + 2) hsub
  have hfac : ((k + 2).factorial : ℂ) ≠ 0 :=
    Nat.cast_ne_zero.mpr (Nat.factorial_ne_zero (k + 2))
  unfold ctCommonDen27 ctS27
  rw [show k + 2 - 1 = k + 1 by omega,
    ctNumerator_succ27 k t, ctPoleProduct_succ27 (k + 3) t]
  field_simp [hden, hfac]
  push_cast
  ring

private theorem ct_common_mul_S_shift27
    (k : ℕ) (t : ℂ)
    (hP : ctPoleProduct27 (k + 4) t ≠ 0) :
    ctCommonDen27 k t * ctS27 (k + 2) (t + 1) =
      ctNumerator27 k t * t ^ 4 *
        ctQhatC27 (((k + 2 : ℕ) : ℂ))
          (t - ((k + 2 : ℕ) : ℂ) + 1) := by
  have hshift : ctPoleProduct27 (k + 3) (t + 1) ≠ 0 := by
    intro hz
    apply hP
    rw [ctPoleProduct_shift_succ27 (k + 3) t, hz, mul_zero]
  have hden :
      (((k + 2).factorial : ℂ) ^ 2) *
          ctPoleProduct27 (k + 3) (t + 1) ≠ 0 :=
    ctDen_ne_zero27 (n := k + 2) hshift
  have hfac : ((k + 2).factorial : ℂ) ≠ 0 :=
    Nat.cast_ne_zero.mpr (Nat.factorial_ne_zero (k + 2))
  unfold ctCommonDen27 ctS27
  rw [show k + 2 - 1 = k + 1 by omega,
    ctNumerator_shift_succ27 k t,
    ctPoleProduct_shift_succ27 (k + 3) t]
  field_simp [hden, hfac]
  push_cast
  ring

/-! ## Final rational telescoper -/

set_option maxRecDepth 100000 in
theorem ctR_telescoper_stepC27
    (k : ℕ) (t : ℂ)
    (hP : ctPoleProduct27 (k + 4) t ≠ 0) :
    ctAlphaC27 (((k + 2 : ℕ) : ℂ)) * ctR27 (k + 3) t
      - ctBetaC27 (((k + 2 : ℕ) : ℂ)) * ctR27 (k + 2) t
      + ctGammaC27 (((k + 2 : ℕ) : ℂ)) * ctR27 (k + 1) t
      - ctDeltaC27 (((k + 2 : ℕ) : ℂ)) * ctR27 k t =
    ctS27 (k + 2) (t + 1) - ctS27 (k + 2) t := by
  have hD : ctCommonDen27 k t ≠ 0 := by
    simpa only [ctCommonDen27] using
      (ctDen_ne_zero27 (n := k + 2) hP)
  apply mul_left_cancel₀ hD
  simp only [mul_sub, mul_add]
  rw [ct_common_mul_alpha_R_next27 k t hP,
    ct_common_mul_beta_R27 k t hP,
    ct_common_mul_gamma_R_prev27 k t hP,
    ct_common_mul_delta_R_prev2_27 k t hP,
    ct_common_mul_S_shift27 k t hP,
    ct_common_mul_S27 k t hP]
  have hpoly := ctQhat_poly_identity27
    (((k + 2 : ℕ) : ℂ))
    (t - (((k + 2 : ℕ) : ℂ)))
  have hmul := congrArg
    (fun z : ℂ => ctNumerator27 k t * z) hpoly.symm
  convert hmul using 1 <;> ring

theorem ctR_telescoper_step27
    (k : ℕ) (t : ℂ)
    (hP : ctPoleProduct27 (k + 4) t ≠ 0) :
    (ctAlpha27 k : ℂ) * ctR27 (k + 3) t
      - (ctBeta27 k : ℂ) * ctR27 (k + 2) t
      + (ctGamma27 k : ℂ) * ctR27 (k + 1) t
      - (ctDelta27 k : ℂ) * ctR27 k t =
    ctS27 (k + 2) (t + 1) - ctS27 (k + 2) t := by
  simpa only [ctAlphaC_nat27, ctBetaC_nat27,
    ctGammaC_nat27, ctDeltaC_nat27] using
      ctR_telescoper_stepC27 k t hP

/-- Paper-indexed wrapper, `n ≥ 2`. -/
theorem ctR_telescoper27
    {n : ℕ} (hn : 2 ≤ n) (t : ℂ)
    (hP : ctPoleProduct27 (n + 2) t ≠ 0) :
    (((2 * zudilinA n : ℚ) : ℂ)) * ctR27 (n + 1) t
      - (((2 * zudilinM n : ℚ) : ℂ)) * ctR27 n t
      + (((2 * (n : ℚ) * zudilinN n : ℚ) : ℂ)) *
          ctR27 (n - 1) t
      - (((zudilinR n * (n : ℚ) * (((n - 1 : ℕ) : ℚ) ^ 3) : ℚ) : ℂ)) *
          ctR27 (n - 2) t =
    ctS27 n (t + 1) - ctS27 n t := by
  obtain ⟨k, rfl⟩ : ∃ k, n = k + 2 := ⟨n - 2, by omega⟩
  simpa [ctAlpha27, ctBeta27, ctGamma27, ctDelta27] using
    ctR_telescoper_step27 k t hP

/-! ## Differentiability on the pole-free domain -/

theorem ctNumerator_differentiableAt27 (n : ℕ) (t : ℂ) :
    DifferentiableAt ℂ (ctNumerator27 n) t := by
  unfold ctNumerator27
  exact DifferentiableAt.fun_finset_prod
    (u := Finset.range n) fun r _ =>
      (differentiableAt_id.sub_const _).fun_pow 3

theorem ctPoleProduct_differentiableAt27 (m : ℕ) (t : ℂ) :
    DifferentiableAt ℂ (ctPoleProduct27 m) t := by
  unfold ctPoleProduct27
  exact DifferentiableAt.fun_finset_prod
    (u := Finset.range m) fun j _ =>
      differentiableAt_id.add_const _

theorem ctR_differentiableAt27
    {n : ℕ} {t : ℂ}
    (hP : ctPoleProduct27 (n + 1) t ≠ 0) :
    DifferentiableAt ℂ (ctR27 n) t := by
  unfold ctR27
  apply DifferentiableAt.div
  · exact ctNumerator_differentiableAt27 n t
  · exact (differentiableAt_const _).mul
      (ctPoleProduct_differentiableAt27 (n + 1) t)
  · exact ctDen_ne_zero27 (n := n) hP

theorem ctS_differentiableAt27
    {n : ℕ} {t : ℂ}
    (hP : ctPoleProduct27 (n + 1) t ≠ 0) :
    DifferentiableAt ℂ (ctS27 n) t := by
  unfold ctS27
  apply DifferentiableAt.div
  · apply DifferentiableAt.mul
    · unfold ctQhatC27
      fun_prop
    · exact ctNumerator_differentiableAt27 (n - 1) t
  · exact (differentiableAt_const _).mul
      (ctPoleProduct_differentiableAt27 (n + 1) t)
  · exact ctDen_ne_zero27 (n := n) hP

/-! ## Local derivative congruence -/

theorem ctR_telescoper_deriv_step27
    (k : ℕ) (t : ℂ)
    (hP : ctPoleProduct27 (k + 4) t ≠ 0) :
    (ctAlpha27 k : ℂ) * deriv (ctR27 (k + 3)) t
      - (ctBeta27 k : ℂ) * deriv (ctR27 (k + 2)) t
      + (ctGamma27 k : ℂ) * deriv (ctR27 (k + 1)) t
      - (ctDelta27 k : ℂ) * deriv (ctR27 k) t =
    deriv (ctS27 (k + 2)) (t + 1)
      - deriv (ctS27 (k + 2)) t := by
  let L : ℂ → ℂ := fun z =>
    (ctAlpha27 k : ℂ) * ctR27 (k + 3) z
      - (ctBeta27 k : ℂ) * ctR27 (k + 2) z
      + (ctGamma27 k : ℂ) * ctR27 (k + 1) z
      - (ctDelta27 k : ℂ) * ctR27 k z
  let T : ℂ → ℂ := fun z =>
    ctS27 (k + 2) (z + 1) - ctS27 (k + 2) z

  have hlocalP :
      ∀ᶠ z in 𝓝 t, ctPoleProduct27 (k + 4) z ≠ 0 :=
    (ctPoleProduct_differentiableAt27 (k + 4) t).continuousAt.eventually_ne hP
  have hEq : L =ᶠ[𝓝 t] T := by
    filter_upwards [hlocalP] with z hz
    simpa only [L, T] using ctR_telescoper_step27 k z hz

  have hP3 : ctPoleProduct27 (k + 3) t ≠ 0 :=
    ctPoleProduct_ne_zero_of_le27
      (m := k + 3) (M := k + 4) (t := t) (by omega) hP
  have hP2 : ctPoleProduct27 (k + 2) t ≠ 0 :=
    ctPoleProduct_ne_zero_of_le27
      (m := k + 2) (M := k + 4) (t := t) (by omega) hP
  have hP1 : ctPoleProduct27 (k + 1) t ≠ 0 :=
    ctPoleProduct_ne_zero_of_le27
      (m := k + 1) (M := k + 4) (t := t) (by omega) hP
  have hshiftP : ctPoleProduct27 (k + 3) (t + 1) ≠ 0 := by
    intro hz
    apply hP
    rw [ctPoleProduct_shift_succ27 (k + 3) t, hz, mul_zero]

  have hR3 : DifferentiableAt ℂ (ctR27 (k + 3)) t :=
    ctR_differentiableAt27 (n := k + 3) (t := t) (by simpa using hP)
  have hR2 : DifferentiableAt ℂ (ctR27 (k + 2)) t :=
    ctR_differentiableAt27 (n := k + 2) (t := t) (by simpa using hP3)
  have hR1 : DifferentiableAt ℂ (ctR27 (k + 1)) t :=
    ctR_differentiableAt27 (n := k + 1) (t := t) (by simpa using hP2)
  have hR0 : DifferentiableAt ℂ (ctR27 k) t :=
    ctR_differentiableAt27 (n := k) (t := t) (by simpa using hP1)
  have hS0 : DifferentiableAt ℂ (ctS27 (k + 2)) t :=
    ctS_differentiableAt27 (n := k + 2) (t := t) (by simpa using hP3)
  have hS1 : DifferentiableAt ℂ (ctS27 (k + 2)) (t + 1) :=
    ctS_differentiableAt27 (n := k + 2) (t := t + 1)
      (by simpa using hshiftP)

  have hL : HasDerivAt L
      ((ctAlpha27 k : ℂ) * deriv (ctR27 (k + 3)) t
        - (ctBeta27 k : ℂ) * deriv (ctR27 (k + 2)) t
        + (ctGamma27 k : ℂ) * deriv (ctR27 (k + 1)) t
        - (ctDelta27 k : ℂ) * deriv (ctR27 k) t) t := by
    dsimp only [L]
    exact (((hR3.hasDerivAt.const_mul (ctAlpha27 k : ℂ)).sub
      (hR2.hasDerivAt.const_mul (ctBeta27 k : ℂ))).add
      (hR1.hasDerivAt.const_mul (ctGamma27 k : ℂ))).sub
      (hR0.hasDerivAt.const_mul (ctDelta27 k : ℂ))

  have hScomp : HasDerivAt
      (fun z : ℂ => ctS27 (k + 2) (z + 1))
      (deriv (ctS27 (k + 2)) (t + 1)) t :=
    HasDerivAt.comp_add_const t 1 hS1.hasDerivAt

  have hT : HasDerivAt T
      (deriv (ctS27 (k + 2)) (t + 1)
        - deriv (ctS27 (k + 2)) t) t := by
    dsimp only [T]
    exact hScomp.sub hS0.hasDerivAt

  calc
    _ = deriv L t := hL.deriv.symm
    _ = deriv T t := hEq.deriv_eq
    _ = _ := hT.deriv

/-! ## Phi identity -/

def ctRPhi27 (n : ℕ) (t : ℂ) : ℂ :=
  ctR27 n t - deriv (ctR27 n) t / 2

def ctSPhi27 (n : ℕ) (t : ℂ) : ℂ :=
  ctS27 n t - deriv (ctS27 n) t / 2

theorem ctPhi_telescoper_step27
    (k : ℕ) (t : ℂ)
    (hP : ctPoleProduct27 (k + 4) t ≠ 0) :
    (ctAlpha27 k : ℂ) * ctRPhi27 (k + 3) t
      - (ctBeta27 k : ℂ) * ctRPhi27 (k + 2) t
      + (ctGamma27 k : ℂ) * ctRPhi27 (k + 1) t
      - (ctDelta27 k : ℂ) * ctRPhi27 k t =
    ctSPhi27 (k + 2) (t + 1) - ctSPhi27 (k + 2) t := by
  have h0 := ctR_telescoper_step27 k t hP
  have h1 := ctR_telescoper_deriv_step27 k t hP
  simp only [ctRPhi27, ctSPhi27]
  linear_combination h0 - h1 / 2

/-! ## Translation to the repository's shifted Barnes coordinate -/

theorem ctNumerator_translate27 (n : ℕ) (s : ℂ) :
    ctNumerator27 n (s + (((n + 1 : ℕ) : ℂ))) =
      (∏ k ∈ Finset.range n,
        zudilinBarnesNumeratorFactor27 k s) ^ 3 := by
  induction n with
  | zero =>
      simp [ctNumerator27]
  | succ n ih =>
      rw [show s + ((((n + 1) + 1 : ℕ) : ℂ)) =
          (s + (((n + 1 : ℕ) : ℂ))) + 1 by
            push_cast
            ring]
      rw [ctNumerator_shift_succ27 n
        (s + (((n + 1 : ℕ) : ℂ))), ih,
        Finset.prod_range_succ]
      simp only [zudilinBarnesNumeratorFactor27]
      push_cast
      ring

theorem ctPoleProduct_translate27 (n : ℕ) (s : ℂ) :
    ctPoleProduct27 (n + 1) (s + (((n + 1 : ℕ) : ℂ))) =
      ∏ k ∈ Finset.range (n + 1),
        zudilinBarnesDenominatorFactor27 n k s := by
  unfold ctPoleProduct27 zudilinBarnesDenominatorFactor27
  apply Finset.prod_congr rfl
  intro k hk
  push_cast
  ring

theorem ctR_translate_eq_zudilinBarnesF27 (n : ℕ) (s : ℂ) :
    ctR27 n (s + (((n + 1 : ℕ) : ℂ))) =
      zudilinBarnesF27 n s := by
  unfold ctR27 zudilinBarnesF27
  rw [ctNumerator_translate27, ctPoleProduct_translate27]

theorem zudilinBarnesF_eq_ctR_translate27 (n : ℕ) (s : ℂ) :
    zudilinBarnesF27 n s =
      ctR27 n (s + (((n + 1 : ℕ) : ℂ))) :=
  (ctR_translate_eq_zudilinBarnesF27 n s).symm

theorem zudilinBarnesFPrime_eq_ctR_deriv_translate27
    (n : ℕ) (s : ℂ) :
    zudilinBarnesFPrime27 n s =
      deriv (ctR27 n) (s + (((n + 1 : ℕ) : ℂ))) := by
  unfold zudilinBarnesFPrime27
  have hfun : zudilinBarnesF27 n =
      fun z : ℂ => ctR27 n (z + (((n + 1 : ℕ) : ℂ))) := by
    funext z
    exact zudilinBarnesF_eq_ctR_translate27 n z
  rw [hfun]
  simpa using
    (deriv_comp_add_const (ctR27 n) (((n + 1 : ℕ) : ℂ)) s)

theorem zudilinBarnesPhi_eq_ctRPhi_translate27
    (n : ℕ) (s : ℂ) :
    zudilinBarnesPhi27 n s =
      ctRPhi27 n (s + (((n + 1 : ℕ) : ℂ))) := by
  unfold zudilinBarnesPhi27 ctRPhi27
  rw [zudilinBarnesF_eq_ctR_translate27,
    zudilinBarnesFPrime_eq_ctR_deriv_translate27]

theorem zudilinBarnesPhi_shift_back27 (n : ℕ) (t : ℂ) :
    zudilinBarnesPhi27 n
        (t - (((n + 1 : ℕ) : ℂ))) =
      ctRPhi27 n t := by
  rw [zudilinBarnesPhi_eq_ctRPhi_translate27]
  congr 1
  ring

/-- The common-pole premise is automatic on the Barnes half-plane. -/
theorem ctPoleProduct_translate_ne_zero27
    {n : ℕ} {s : ℂ} (hs : -(1 / 2 : ℝ) ≤ s.re) :
    ctPoleProduct27 (n + 1)
        (s + (((n + 1 : ℕ) : ℂ))) ≠ 0 := by
  rw [ctPoleProduct_translate27]
  exact Finset.prod_ne_zero_iff.mpr fun k hk =>
    zudilinBarnesDenominatorFactor_ne_zero27
      (n := n) (k := k) (s := s) hs

/-! ## Pointwise recurrence in the shifted Barnes coordinate -/

theorem zudilinBarnesPhi_telescoper_commonT27
    (k : ℕ) (t : ℂ)
    (hP : ctPoleProduct27 (k + 4) t ≠ 0) :
    (ctAlpha27 k : ℂ) *
        zudilinBarnesPhi27 (k + 3) (t - (((k + 4 : ℕ) : ℂ)))
      - (ctBeta27 k : ℂ) *
        zudilinBarnesPhi27 (k + 2) (t - (((k + 3 : ℕ) : ℂ)))
      + (ctGamma27 k : ℂ) *
        zudilinBarnesPhi27 (k + 1) (t - (((k + 2 : ℕ) : ℂ)))
      - (ctDelta27 k : ℂ) *
        zudilinBarnesPhi27 k (t - (((k + 1 : ℕ) : ℂ))) =
    ctSPhi27 (k + 2) (t + 1) - ctSPhi27 (k + 2) t := by
  simpa only [zudilinBarnesPhi_shift_back27] using
    ctPhi_telescoper_step27 k t hP

theorem zudilinBarnesPhi_telescoper_fixedS27
    (k : ℕ) (s : ℂ)
    (hP : ctPoleProduct27 (k + 4)
      (s + (((k + 4 : ℕ) : ℂ))) ≠ 0) :
    (ctAlpha27 k : ℂ) * zudilinBarnesPhi27 (k + 3) s
      - (ctBeta27 k : ℂ) * zudilinBarnesPhi27 (k + 2) (s + 1)
      + (ctGamma27 k : ℂ) * zudilinBarnesPhi27 (k + 1) (s + 2)
      - (ctDelta27 k : ℂ) * zudilinBarnesPhi27 k (s + 3) =
    ctSPhi27 (k + 2) (s + (((k + 5 : ℕ) : ℂ)))
      - ctSPhi27 (k + 2) (s + (((k + 4 : ℕ) : ℂ))) := by
  have h := zudilinBarnesPhi_telescoper_commonT27 k
    (s + (((k + 4 : ℕ) : ℂ))) hP
  convert h using 1 <;> push_cast <;> ring

theorem zudilinBarnesPhi_telescoper_halfplane27
    (k : ℕ) (s : ℂ) (hs : -(1 / 2 : ℝ) ≤ s.re) :
    (ctAlpha27 k : ℂ) * zudilinBarnesPhi27 (k + 3) s
      - (ctBeta27 k : ℂ) * zudilinBarnesPhi27 (k + 2) (s + 1)
      + (ctGamma27 k : ℂ) * zudilinBarnesPhi27 (k + 1) (s + 2)
      - (ctDelta27 k : ℂ) * zudilinBarnesPhi27 k (s + 3) =
    ctSPhi27 (k + 2) (s + (((k + 5 : ℕ) : ℂ)))
      - ctSPhi27 (k + 2) (s + (((k + 4 : ℕ) : ℂ))) := by
  apply zudilinBarnesPhi_telescoper_fixedS27 k s
  have h := ctPoleProduct_translate_ne_zero27 (n := k + 3) (s := s) hs
  convert h using 1

#print axioms ctQhat_poly_identity27
#print axioms ctR_telescoper_stepC27
#print axioms ctR_telescoper_step27
#print axioms ctR_telescoper_deriv_step27
#print axioms ctPhi_telescoper_step27
#print axioms zudilinBarnesPhi_telescoper_halfplane27

end RamanujanChallenge.P27
