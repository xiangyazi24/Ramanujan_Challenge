import RamanujanChallenge.Problem27
import Mathlib.Analysis.Calculus.LogDeriv
import Mathlib.Analysis.SpecialFunctions.ImproperIntegrals
import Mathlib.Analysis.SpecialFunctions.Trigonometric.DerivHyp
import Mathlib.Analysis.SpecialFunctions.Trigonometric.EulerSineProd
import Mathlib.MeasureTheory.Integral.IntervalIntegral.IntegrationByParts

/-!
# Problem 2.7: midpoint estimates for the Zudilin Barnes kernel

This file isolates the analytic certificate used to prove decay of
`zudilinCombinedError`.  The first layer is elementary: the exact midpoint
normalization has a one-step ratio bounded by `1 / 4`.
-/

open Filter Set Topology

noncomputable section

namespace RamanujanChallenge.P27

/-- The exact one-step ratio of the rational function at the midpoint
`n + 1/2` of its Barnes contour. -/
def zudilinMidpointRatio27 (n : ℕ) : ℝ :=
  ((n : ℝ) + 1 / 2) ^ 4 /
    ((2 * n + 3 / 2 : ℝ) * (2 * n + 5 / 2 : ℝ) *
      ((n : ℝ) + 1) ^ 2)

/-- The positive midpoint normalization.  Its product form is established
later from the finite-product definition of the Barnes rational function. -/
def zudilinMidpointValue27 : ℕ → ℝ
  | 0 => 2
  | n + 1 => zudilinMidpointValue27 n * zudilinMidpointRatio27 n

@[simp] theorem zudilinMidpointValue_zero27 :
    zudilinMidpointValue27 0 = 2 := rfl

@[simp] theorem zudilinMidpointValue_succ27 (n : ℕ) :
    zudilinMidpointValue27 (n + 1) =
      zudilinMidpointValue27 n * zudilinMidpointRatio27 n := rfl

theorem zudilinMidpointRatio_pos27 (n : ℕ) :
    0 < zudilinMidpointRatio27 n := by
  unfold zudilinMidpointRatio27
  positivity

theorem zudilinMidpointRatio_le_quarter27 (n : ℕ) :
    zudilinMidpointRatio27 n ≤ 1 / 4 := by
  have hn : 0 ≤ (n : ℝ) := Nat.cast_nonneg n
  have hden :
      0 < (2 * n + 3 / 2 : ℝ) * (2 * n + 5 / 2 : ℝ) *
        ((n : ℝ) + 1) ^ 2 := by
    positivity
  unfold zudilinMidpointRatio27
  rw [div_le_iff₀ hden]
  nlinarith [sq_nonneg ((n : ℝ) ^ 2), sq_nonneg (n : ℝ)]

theorem zudilinMidpointValue_pos27 (n : ℕ) :
    0 < zudilinMidpointValue27 n := by
  induction n with
  | zero => norm_num
  | succ n ih =>
      rw [zudilinMidpointValue_succ27]
      exact mul_pos ih (zudilinMidpointRatio_pos27 n)

theorem zudilinMidpointValue_le27 (n : ℕ) :
    zudilinMidpointValue27 n ≤ 2 * (1 / 4 : ℝ) ^ n := by
  induction n with
  | zero => norm_num
  | succ n ih =>
      rw [zudilinMidpointValue_succ27, pow_succ]
      calc
        zudilinMidpointValue27 n * zudilinMidpointRatio27 n ≤
            zudilinMidpointValue27 n * (1 / 4 : ℝ) :=
          mul_le_mul_of_nonneg_left
            (zudilinMidpointRatio_le_quarter27 n)
            (zudilinMidpointValue_pos27 n).le
        _ ≤ (2 * (1 / 4 : ℝ) ^ n) * (1 / 4 : ℝ) :=
          mul_le_mul_of_nonneg_right ih (by norm_num)
        _ = 2 * ((1 / 4 : ℝ) ^ n * (1 / 4 : ℝ)) := by ring

theorem zudilinMidpointValue_tendsto_zero27 :
    Tendsto zudilinMidpointValue27 atTop (nhds 0) := by
  have hnonneg : ∀ n, 0 ≤ zudilinMidpointValue27 n :=
    fun n => (zudilinMidpointValue_pos27 n).le
  have hgeom :
      Tendsto (fun n : ℕ => 2 * (1 / 4 : ℝ) ^ n) atTop (nhds 0) := by
    simpa using
      (tendsto_const_nhds.mul
        (tendsto_pow_atTop_nhds_zero_of_norm_lt_one
          (by norm_num : ‖(1 / 4 : ℝ)‖ < 1)))
  exact squeeze_zero' (Eventually.of_forall hnonneg)
    (Eventually.of_forall zudilinMidpointValue_le27) hgeom

/-! ## The fixed Barnes rational function -/

/-- A numerator factor of the fixed-line translate of Zudilin's rational
function. -/
def zudilinBarnesNumeratorFactor27 (k : ℕ) (s : ℂ) : ℂ :=
  s + (k + 1 : ℕ)

/-- A denominator factor of the fixed-line translate of Zudilin's rational
function. -/
def zudilinBarnesDenominatorFactor27 (n k : ℕ) (s : ℂ) : ℂ :=
  s + (n + 1 + k : ℕ)

/-- The fixed-line translate `R_n(s+n+1)` of Zudilin's Barnes rational
function.  Its zeros and poles lie strictly to the left of `re s = -1/2`. -/
def zudilinBarnesF27 (n : ℕ) (s : ℂ) : ℂ :=
  ((∏ k ∈ Finset.range n, zudilinBarnesNumeratorFactor27 k s) ^ 3) /
    (((n.factorial : ℂ) ^ 2) *
      ∏ k ∈ Finset.range (n + 1), zudilinBarnesDenominatorFactor27 n k s)

/-- The actual complex derivative of `zudilinBarnesF27`. -/
def zudilinBarnesFPrime27 (n : ℕ) (s : ℂ) : ℂ :=
  deriv (zudilinBarnesF27 n) s

/-- The logarithmic-derivative expression, used only where every displayed
factor is nonzero. -/
def zudilinBarnesLogDerivative27 (n : ℕ) (s : ℂ) : ℂ :=
  3 * ∑ k ∈ Finset.range n, (zudilinBarnesNumeratorFactor27 k s)⁻¹ -
    ∑ k ∈ Finset.range (n + 1),
      (zudilinBarnesDenominatorFactor27 n k s)⁻¹

/-- The combined Barnes factor for the sum of the zeta-two and zeta-three
linear forms. -/
def zudilinBarnesPhi27 (n : ℕ) (s : ℂ) : ℂ :=
  zudilinBarnesF27 n s - zudilinBarnesFPrime27 n s / 2

private theorem complex_ne_zero_of_re_pos27 {z : ℂ} (hz : 0 < z.re) : z ≠ 0 := by
  intro h
  rw [h] at hz
  norm_num at hz

theorem zudilinBarnesNumeratorFactor_re27 (k : ℕ) (s : ℂ) :
    (zudilinBarnesNumeratorFactor27 k s).re = s.re + (k + 1 : ℕ) := by
  simp [zudilinBarnesNumeratorFactor27]

theorem zudilinBarnesDenominatorFactor_re27 (n k : ℕ) (s : ℂ) :
    (zudilinBarnesDenominatorFactor27 n k s).re =
      s.re + (n + 1 + k : ℕ) := by
  simp [zudilinBarnesDenominatorFactor27]

theorem zudilinBarnesNumeratorFactor_ne_zero27
    {k : ℕ} {s : ℂ} (hs : -(1 / 2 : ℝ) ≤ s.re) :
    zudilinBarnesNumeratorFactor27 k s ≠ 0 := by
  apply complex_ne_zero_of_re_pos27
  rw [zudilinBarnesNumeratorFactor_re27]
  have hk : (0 : ℝ) ≤ k := Nat.cast_nonneg k
  norm_num at hs ⊢
  linarith

theorem zudilinBarnesDenominatorFactor_ne_zero27
    {n k : ℕ} {s : ℂ} (hs : -(1 / 2 : ℝ) ≤ s.re) :
    zudilinBarnesDenominatorFactor27 n k s ≠ 0 := by
  apply complex_ne_zero_of_re_pos27
  rw [zudilinBarnesDenominatorFactor_re27]
  have hn : (0 : ℝ) ≤ n := Nat.cast_nonneg n
  have hk : (0 : ℝ) ≤ k := Nat.cast_nonneg k
  norm_num at hs ⊢
  linarith

theorem zudilinBarnesF_ne_zero27
    {n : ℕ} {s : ℂ} (hs : -(1 / 2 : ℝ) ≤ s.re) :
    zudilinBarnesF27 n s ≠ 0 := by
  unfold zudilinBarnesF27
  apply div_ne_zero
  · exact pow_ne_zero _ (Finset.prod_ne_zero_iff.mpr fun k hk ↦
      zudilinBarnesNumeratorFactor_ne_zero27 hs)
  · apply mul_ne_zero
    · exact pow_ne_zero _ (Nat.cast_ne_zero.mpr (Nat.factorial_ne_zero n))
    · exact Finset.prod_ne_zero_iff.mpr fun k hk ↦
        zudilinBarnesDenominatorFactor_ne_zero27 hs

theorem zudilinBarnesNumeratorFactor_differentiableAt27 (k : ℕ) (s : ℂ) :
    DifferentiableAt ℂ (zudilinBarnesNumeratorFactor27 k) s := by
  unfold zudilinBarnesNumeratorFactor27
  exact differentiableAt_id.add_const _

theorem zudilinBarnesDenominatorFactor_differentiableAt27
    (n k : ℕ) (s : ℂ) :
    DifferentiableAt ℂ (zudilinBarnesDenominatorFactor27 n k) s := by
  unfold zudilinBarnesDenominatorFactor27
  exact differentiableAt_id.add_const _

theorem zudilinBarnesF_differentiableAt27
    {n : ℕ} {s : ℂ} (hs : -(1 / 2 : ℝ) ≤ s.re) :
    DifferentiableAt ℂ (zudilinBarnesF27 n) s := by
  unfold zudilinBarnesF27
  apply DifferentiableAt.div
  · exact (DifferentiableAt.fun_finset_prod (u := Finset.range n) fun k hk ↦
      zudilinBarnesNumeratorFactor_differentiableAt27 k s).fun_pow 3
  · exact DifferentiableAt.mul (differentiableAt_const _)
      (DifferentiableAt.fun_finset_prod (u := Finset.range (n + 1)) fun k hk ↦
        zudilinBarnesDenominatorFactor_differentiableAt27 n k s)
  · apply mul_ne_zero
    · exact pow_ne_zero _ (Nat.cast_ne_zero.mpr (Nat.factorial_ne_zero n))
    · exact Finset.prod_ne_zero_iff.mpr fun k hk ↦
        zudilinBarnesDenominatorFactor_ne_zero27 hs

private theorem zudilinBarnesNumeratorFactor_logDeriv27 (k : ℕ) (s : ℂ) :
    logDeriv (zudilinBarnesNumeratorFactor27 k) s =
      (zudilinBarnesNumeratorFactor27 k s)⁻¹ := by
  rw [logDeriv_apply]
  have hderiv : deriv (zudilinBarnesNumeratorFactor27 k) s = 1 := by
    unfold zudilinBarnesNumeratorFactor27
    simp
  rw [hderiv]
  simp [one_div]

private theorem zudilinBarnesDenominatorFactor_logDeriv27
    (n k : ℕ) (s : ℂ) :
    logDeriv (zudilinBarnesDenominatorFactor27 n k) s =
      (zudilinBarnesDenominatorFactor27 n k s)⁻¹ := by
  rw [logDeriv_apply]
  have hderiv : deriv (zudilinBarnesDenominatorFactor27 n k) s = 1 := by
    unfold zudilinBarnesDenominatorFactor27
    simp
  rw [hderiv]
  simp [one_div]

theorem zudilinBarnesF_logDeriv27
    {n : ℕ} {s : ℂ} (hs : -(1 / 2 : ℝ) ≤ s.re) :
    logDeriv (zudilinBarnesF27 n) s =
      zudilinBarnesLogDerivative27 n s := by
  let N : ℂ → ℂ := fun z ↦
    ∏ k ∈ Finset.range n, zudilinBarnesNumeratorFactor27 k z
  let D : ℂ → ℂ := fun z ↦
    ∏ k ∈ Finset.range (n + 1), zudilinBarnesDenominatorFactor27 n k z
  have hN0 : N s ≠ 0 := by
    exact Finset.prod_ne_zero_iff.mpr fun k hk ↦
      zudilinBarnesNumeratorFactor_ne_zero27 hs
  have hD0 : D s ≠ 0 := by
    exact Finset.prod_ne_zero_iff.mpr fun k hk ↦
      zudilinBarnesDenominatorFactor_ne_zero27 hs
  have hNd : DifferentiableAt ℂ N s := by
    simpa only [N] using
      (DifferentiableAt.fun_finset_prod (u := Finset.range n) fun k hk ↦
        zudilinBarnesNumeratorFactor_differentiableAt27 k s)
  have hDd : DifferentiableAt ℂ D s := by
    simpa only [D] using
      (DifferentiableAt.fun_finset_prod (u := Finset.range (n + 1)) fun k hk ↦
        zudilinBarnesDenominatorFactor_differentiableAt27 n k s)
  have hNlog :
      logDeriv N s =
        ∑ k ∈ Finset.range n,
          (zudilinBarnesNumeratorFactor27 k s)⁻¹ := by
    rw [show N = fun z ↦
      ∏ k ∈ Finset.range n, zudilinBarnesNumeratorFactor27 k z by rfl]
    rw [logDeriv_prod]
    · exact Finset.sum_congr rfl fun k hk ↦
        zudilinBarnesNumeratorFactor_logDeriv27 k s
    · exact fun k hk ↦ zudilinBarnesNumeratorFactor_ne_zero27 hs
    · exact fun k hk ↦
        zudilinBarnesNumeratorFactor_differentiableAt27 k s
  have hDlog :
      logDeriv D s =
        ∑ k ∈ Finset.range (n + 1),
          (zudilinBarnesDenominatorFactor27 n k s)⁻¹ := by
    rw [show D = fun z ↦
      ∏ k ∈ Finset.range (n + 1),
        zudilinBarnesDenominatorFactor27 n k z by rfl]
    rw [logDeriv_prod]
    · exact Finset.sum_congr rfl fun k hk ↦
        zudilinBarnesDenominatorFactor_logDeriv27 n k s
    · exact fun k hk ↦ zudilinBarnesDenominatorFactor_ne_zero27 hs
    · exact fun k hk ↦
        zudilinBarnesDenominatorFactor_differentiableAt27 n k s
  have hfac0 : ((n.factorial : ℂ) ^ 2) ≠ 0 :=
    pow_ne_zero _ (Nat.cast_ne_zero.mpr (Nat.factorial_ne_zero n))
  have hden0 : ((n.factorial : ℂ) ^ 2) * D s ≠ 0 :=
    mul_ne_zero hfac0 hD0
  have hdenDiff :
      DifferentiableAt ℂ (fun z ↦ ((n.factorial : ℂ) ^ 2) * D z) s :=
    DifferentiableAt.mul (differentiableAt_const _) hDd
  have hNpow0 : N s ^ 3 ≠ 0 := pow_ne_zero 3 hN0
  have hPowLog :
      logDeriv (N ^ 3) s = 3 * logDeriv N s := by
    simpa only [Pi.pow_apply] using (logDeriv_fun_pow hNd 3)
  have hDenLog :
      logDeriv (fun z ↦ ((n.factorial : ℂ) ^ 2) * D z) s = logDeriv D s :=
    logDeriv_const_mul s ((n.factorial : ℂ) ^ 2) hfac0
  have hQuot :=
    logDeriv_div s hNpow0 hden0 (hNd.pow 3) hdenDiff
  rw [hPowLog, hDenLog, hNlog, hDlog] at hQuot
  simpa only [zudilinBarnesF27, N, D, zudilinBarnesLogDerivative27] using hQuot

theorem zudilinBarnesFPrime_eq_mul_logDerivative27
    {n : ℕ} {s : ℂ} (hs : -(1 / 2 : ℝ) ≤ s.re) :
    zudilinBarnesFPrime27 n s =
      zudilinBarnesF27 n s * zudilinBarnesLogDerivative27 n s := by
  have hlog := zudilinBarnesF_logDeriv27 (n := n) hs
  have hF0 := zudilinBarnesF_ne_zero27 (n := n) hs
  simp only [logDeriv_apply, zudilinBarnesFPrime27] at hlog ⊢
  apply (div_eq_iff hF0).mp at hlog
  rw [hlog, mul_comm]

/-! ## Bounds on the fixed vertical line -/

/-- The fixed Barnes line `re s = -1/2`. -/
def zudilinBarnesLine27 (y : ℝ) : ℂ :=
  (-1 / 2 : ℝ) + y * Complex.I

@[simp] theorem zudilinBarnesLine_re27 (y : ℝ) :
    (zudilinBarnesLine27 y).re = -(1 / 2 : ℝ) := by
  norm_num [zudilinBarnesLine27]

@[simp] theorem zudilinBarnesLine_im27 (y : ℝ) :
    (zudilinBarnesLine27 y).im = y := by
  simp [zudilinBarnesLine27]

private theorem norm_inv_le_two_of_half_le_re27 {z : ℂ}
    (hz : (1 / 2 : ℝ) ≤ z.re) : ‖z⁻¹‖ ≤ 2 := by
  have hnorm : (1 / 2 : ℝ) ≤ ‖z‖ := by
    calc
      (1 / 2 : ℝ) ≤ z.re := hz
      _ ≤ |z.re| := le_abs_self _
      _ ≤ ‖z‖ := Complex.abs_re_le_norm z
  have hnormpos : 0 < ‖z‖ := lt_of_lt_of_le (by norm_num) hnorm
  rw [norm_inv]
  apply (inv_le_iff_one_le_mul₀' hnormpos).2
  nlinarith

theorem norm_inv_zudilinBarnesNumeratorFactor_le_two27
    (k : ℕ) (y : ℝ) :
    ‖(zudilinBarnesNumeratorFactor27 k (zudilinBarnesLine27 y))⁻¹‖ ≤ 2 := by
  apply norm_inv_le_two_of_half_le_re27
  rw [zudilinBarnesNumeratorFactor_re27, zudilinBarnesLine_re27]
  have hk : (0 : ℝ) ≤ k := Nat.cast_nonneg k
  norm_num at *

theorem norm_inv_zudilinBarnesDenominatorFactor_le_two27
    (n k : ℕ) (y : ℝ) :
    ‖(zudilinBarnesDenominatorFactor27 n k (zudilinBarnesLine27 y))⁻¹‖ ≤ 2 := by
  apply norm_inv_le_two_of_half_le_re27
  rw [zudilinBarnesDenominatorFactor_re27, zudilinBarnesLine_re27]
  have hn : (0 : ℝ) ≤ n := Nat.cast_nonneg n
  have hk : (0 : ℝ) ≤ k := Nat.cast_nonneg k
  norm_num at *
  linarith

private theorem norm_sum_zudilinBarnesNumerator_inv_le27
    (n : ℕ) (y : ℝ) :
    ‖∑ k ∈ Finset.range n,
        (zudilinBarnesNumeratorFactor27 k (zudilinBarnesLine27 y))⁻¹‖ ≤
      2 * n := by
  calc
    ‖∑ k ∈ Finset.range n,
        (zudilinBarnesNumeratorFactor27 k (zudilinBarnesLine27 y))⁻¹‖ ≤
        ∑ k ∈ Finset.range n,
          ‖(zudilinBarnesNumeratorFactor27 k (zudilinBarnesLine27 y))⁻¹‖ :=
      norm_sum_le _ _
    _ ≤ ∑ _k ∈ Finset.range n, (2 : ℝ) := by
      exact Finset.sum_le_sum fun k hk ↦
        norm_inv_zudilinBarnesNumeratorFactor_le_two27 k y
    _ = 2 * n := by simp [mul_comm]

private theorem norm_sum_zudilinBarnesDenominator_inv_le27
    (n : ℕ) (y : ℝ) :
    ‖∑ k ∈ Finset.range (n + 1),
        (zudilinBarnesDenominatorFactor27 n k (zudilinBarnesLine27 y))⁻¹‖ ≤
      2 * (n + 1) := by
  calc
    ‖∑ k ∈ Finset.range (n + 1),
        (zudilinBarnesDenominatorFactor27 n k (zudilinBarnesLine27 y))⁻¹‖ ≤
        ∑ k ∈ Finset.range (n + 1),
          ‖(zudilinBarnesDenominatorFactor27 n k (zudilinBarnesLine27 y))⁻¹‖ :=
      norm_sum_le _ _
    _ ≤ ∑ _k ∈ Finset.range (n + 1), (2 : ℝ) := by
      exact Finset.sum_le_sum fun k hk ↦
        norm_inv_zudilinBarnesDenominatorFactor_le_two27 n k y
    _ = 2 * (n + 1) := by simp [mul_comm]

theorem norm_zudilinBarnesLogDerivative_le27 (n : ℕ) (y : ℝ) :
    ‖zudilinBarnesLogDerivative27 n (zudilinBarnesLine27 y)‖ ≤
      8 * n + 2 := by
  unfold zudilinBarnesLogDerivative27
  calc
    ‖3 * ∑ k ∈ Finset.range n,
          (zudilinBarnesNumeratorFactor27 k (zudilinBarnesLine27 y))⁻¹ -
        ∑ k ∈ Finset.range (n + 1),
          (zudilinBarnesDenominatorFactor27 n k (zudilinBarnesLine27 y))⁻¹‖ ≤
        ‖3 * ∑ k ∈ Finset.range n,
          (zudilinBarnesNumeratorFactor27 k (zudilinBarnesLine27 y))⁻¹‖ +
        ‖∑ k ∈ Finset.range (n + 1),
          (zudilinBarnesDenominatorFactor27 n k (zudilinBarnesLine27 y))⁻¹‖ :=
      norm_sub_le _ _
    _ ≤ 3 * (2 * n) + 2 * (n + 1) := by
      rw [norm_mul]
      norm_num
      exact add_le_add
        (mul_le_mul_of_nonneg_left
          (norm_sum_zudilinBarnesNumerator_inv_le27 n y) (by norm_num))
        (norm_sum_zudilinBarnesDenominator_inv_le27 n y)
    _ = 8 * n + 2 := by ring

theorem norm_zudilinBarnesFPrime_le27 (n : ℕ) (y : ℝ) :
    ‖zudilinBarnesFPrime27 n (zudilinBarnesLine27 y)‖ ≤
      (8 * n + 2) * ‖zudilinBarnesF27 n (zudilinBarnesLine27 y)‖ := by
  rw [zudilinBarnesFPrime_eq_mul_logDerivative27 (by simp)]
  rw [norm_mul]
  rw [mul_comm ‖zudilinBarnesF27 n (zudilinBarnesLine27 y)‖]
  exact mul_le_mul_of_nonneg_right
    (norm_zudilinBarnesLogDerivative_le27 n y) (norm_nonneg _)

theorem norm_zudilinBarnesPhi_le27 (n : ℕ) (y : ℝ) :
    ‖zudilinBarnesPhi27 n (zudilinBarnesLine27 y)‖ ≤
      (4 * n + 2) * ‖zudilinBarnesF27 n (zudilinBarnesLine27 y)‖ := by
  rw [zudilinBarnesPhi27]
  calc
    ‖zudilinBarnesF27 n (zudilinBarnesLine27 y) -
        zudilinBarnesFPrime27 n (zudilinBarnesLine27 y) / 2‖ ≤
      ‖zudilinBarnesF27 n (zudilinBarnesLine27 y)‖ +
        ‖zudilinBarnesFPrime27 n (zudilinBarnesLine27 y) / 2‖ :=
      norm_sub_le _ _
    _ ≤ ‖zudilinBarnesF27 n (zudilinBarnesLine27 y)‖ +
        ((8 * n + 2) * ‖zudilinBarnesF27 n (zudilinBarnesLine27 y)‖) / 2 := by
      rw [norm_div]
      norm_num
      exact div_le_div_of_nonneg_right (norm_zudilinBarnesFPrime_le27 n y)
        (by norm_num : (0 : ℝ) ≤ 2)
    _ = (4 * n + 2) * ‖zudilinBarnesF27 n (zudilinBarnesLine27 y)‖ := by
      ring

/-! ## Identification of the midpoint normalization -/

def zudilinHalfFactor27 (k : ℕ) : ℝ :=
  (k : ℝ) + 1 / 2

/-- The midpoint value in a product form whose one-step quotient is
transparent. -/
def zudilinMidpointProduct27 (n : ℕ) : ℝ :=
  ((∏ k ∈ Finset.range n, zudilinHalfFactor27 k) ^ 4) /
    (((n.factorial : ℝ) ^ 2) *
      ∏ k ∈ Finset.range (2 * n + 1), zudilinHalfFactor27 k)

/-- The product form obtained directly by substituting the midpoint into the
fixed Barnes rational function. -/
def zudilinMidpointDirect27 (n : ℕ) : ℝ :=
  ((∏ k ∈ Finset.range n, zudilinHalfFactor27 k) ^ 3) /
    (((n.factorial : ℝ) ^ 2) *
      ∏ k ∈ Finset.range (n + 1), zudilinHalfFactor27 (n + k))

theorem zudilinHalfFactor_pos27 (k : ℕ) :
    0 < zudilinHalfFactor27 k := by
  unfold zudilinHalfFactor27
  positivity

theorem zudilinHalfFactor_ne_zero27 (k : ℕ) :
    zudilinHalfFactor27 k ≠ 0 :=
  ne_of_gt (zudilinHalfFactor_pos27 k)

theorem zudilinMidpointProduct_pos27 (n : ℕ) :
    0 < zudilinMidpointProduct27 n := by
  unfold zudilinMidpointProduct27
  apply div_pos
  · exact pow_pos (Finset.prod_pos fun k hk ↦ zudilinHalfFactor_pos27 k) _
  · exact mul_pos (pow_pos (Nat.cast_pos.mpr (Nat.factorial_pos n)) _)
      (Finset.prod_pos fun k hk ↦ zudilinHalfFactor_pos27 k)

@[simp] theorem zudilinMidpointProduct_zero27 :
    zudilinMidpointProduct27 0 = 2 := by
  norm_num [zudilinMidpointProduct27, zudilinHalfFactor27]

theorem zudilinMidpointProduct_succ27 (n : ℕ) :
    zudilinMidpointProduct27 (n + 1) =
      zudilinMidpointProduct27 n * zudilinMidpointRatio27 n := by
  have hindex : 2 * (n + 1) + 1 = (2 * n + 1) + 2 := by omega
  have hfac : ((n + 1).factorial : ℝ) = (n + 1) * (n.factorial : ℝ) := by
    exact_mod_cast Nat.factorial_succ n
  rw [zudilinMidpointProduct27, zudilinMidpointProduct27,
    zudilinMidpointRatio27, hindex]
  rw [show (2 * n + 1) + 2 = ((2 * n + 1) + 1) + 1 by omega]
  simp only [Finset.prod_range_succ, hfac, zudilinHalfFactor27]
  field_simp [Nat.cast_ne_zero.mpr (Nat.factorial_ne_zero n)]
  push_cast
  ring

theorem zudilinMidpointProduct_eq_value27 (n : ℕ) :
    zudilinMidpointProduct27 n = zudilinMidpointValue27 n := by
  induction n with
  | zero => simp
  | succ n ih =>
      rw [zudilinMidpointProduct_succ27, zudilinMidpointValue_succ27, ih]

theorem zudilinMidpointDirect_eq_product27 (n : ℕ) :
    zudilinMidpointDirect27 n = zudilinMidpointProduct27 n := by
  let A : ℝ := ∏ k ∈ Finset.range n, zudilinHalfFactor27 k
  let D : ℝ :=
    ∏ k ∈ Finset.range (n + 1), zudilinHalfFactor27 (n + k)
  have hA0 : A ≠ 0 := by
    exact Finset.prod_ne_zero_iff.mpr fun k hk ↦ zudilinHalfFactor_ne_zero27 k
  have hsplit :
      (∏ k ∈ Finset.range (2 * n + 1), zudilinHalfFactor27 k) = A * D := by
    have h := Finset.prod_range_add zudilinHalfFactor27 n (n + 1)
    simpa only [A, D, show n + (n + 1) = 2 * n + 1 by omega] using h
  rw [zudilinMidpointDirect27, zudilinMidpointProduct27]
  change A ^ 3 / (((n.factorial : ℝ) ^ 2) * D) =
    A ^ 4 /
      (((n.factorial : ℝ) ^ 2) *
        ∏ k ∈ Finset.range (2 * n + 1), zudilinHalfFactor27 k)
  rw [hsplit]
  field_simp [hA0, Nat.cast_ne_zero.mpr (Nat.factorial_ne_zero n)]

theorem zudilinBarnesF_line_zero_eq_midpoint27 (n : ℕ) :
    zudilinBarnesF27 n (zudilinBarnesLine27 0) =
      (zudilinMidpointValue27 n : ℂ) := by
  have hdirect :
      zudilinBarnesF27 n (zudilinBarnesLine27 0) =
        (zudilinMidpointDirect27 n : ℂ) := by
    simp only [zudilinBarnesF27, zudilinBarnesLine27,
      zudilinBarnesNumeratorFactor27, zudilinBarnesDenominatorFactor27,
      zudilinMidpointDirect27, zudilinHalfFactor27]
    push_cast
    ring
  rw [hdirect, zudilinMidpointDirect_eq_product27,
    zudilinMidpointProduct_eq_value27]

theorem norm_zudilinBarnesF_line_zero27 (n : ℕ) :
    ‖zudilinBarnesF27 n (zudilinBarnesLine27 0)‖ =
      zudilinMidpointValue27 n := by
  rw [zudilinBarnesF_line_zero_eq_midpoint27, Complex.norm_real,
    Real.norm_eq_abs, abs_of_pos (zudilinMidpointValue_pos27 n)]

/-! ## Euler's product on the half-integers

The fixed-line estimate uses the half-integer Euler product

`prod (1 + y² / (k + 1/2)²) = cosh (π y)`.

We derive it from Mathlib's Euler sine product rather than taking it as an
additional analytic input. -/

private def zudilinOddEulerProduct27 (n : ℕ) (z : ℂ) : ℂ :=
  ∏ k ∈ Finset.range n,
    (1 - 4 * z ^ 2 / ((2 * k + 1 : ℕ) : ℂ) ^ 2)

def zudilinHalfEulerProduct27 (n : ℕ) (y : ℝ) : ℝ :=
  ∏ k ∈ Finset.range n,
    (1 + y ^ 2 / ((k : ℝ) + 1 / 2) ^ 2)

private theorem zudilinEulerProduct_even_odd27 (n : ℕ) (z : ℂ) :
    (∏ j ∈ Finset.range (2 * n),
        (1 - (2 * z) ^ 2 / ((j : ℂ) + 1) ^ 2)) =
      (∏ k ∈ Finset.range n,
        (1 - z ^ 2 / ((k : ℂ) + 1) ^ 2)) *
          zudilinOddEulerProduct27 n z := by
  induction n with
  | zero => simp [zudilinOddEulerProduct27]
  | succ n ih =>
      rw [show 2 * (n + 1) = (2 * n + 1) + 1 by omega]
      simp only [Finset.prod_range_succ]
      have hop :
          zudilinOddEulerProduct27 (n + 1) z =
            zudilinOddEulerProduct27 n z *
              (1 - 4 * z ^ 2 / (((2 * n + 1 : ℕ) : ℂ) ^ 2)) := by
        unfold zudilinOddEulerProduct27
        rw [Finset.prod_range_succ]
      have hodd :
          1 - (2 * z) ^ 2 / (((2 * n : ℕ) : ℂ) + 1) ^ 2 =
            1 - 4 * z ^ 2 / (((2 * n + 1 : ℕ) : ℂ) ^ 2) := by
        push_cast
        ring
      have heven :
          1 - (2 * z) ^ 2 / (((2 * n + 1 : ℕ) : ℂ) + 1) ^ 2 =
            1 - z ^ 2 / ((n : ℂ) + 1) ^ 2 := by
        have hcast : (((2 * n + 1 : ℕ) : ℂ) + 1) =
            2 * ((n : ℂ) + 1) := by push_cast; ring
        rw [hcast]
        field_simp [Nat.cast_add_one_ne_zero (R := ℂ) n]
      rw [ih, hodd, heven, hop]
      ring

private theorem zudilinOddEulerProduct_eq_half27 (n : ℕ) (y : ℝ) :
    zudilinOddEulerProduct27 n ((y : ℂ) * Complex.I) =
      (zudilinHalfEulerProduct27 n y : ℂ) := by
  unfold zudilinOddEulerProduct27 zudilinHalfEulerProduct27
  rw [Complex.ofReal_prod]
  apply Finset.prod_congr rfl
  intro k hk
  push_cast
  field_simp [show (2 * (k : ℂ) + 1) ≠ 0 by
    simpa only [Nat.cast_mul, Nat.cast_ofNat] using
      (Nat.cast_add_one_ne_zero (R := ℂ) (2 * k))]
  rw [Complex.I_sq]
  ring

theorem zudilinHalfEulerProduct_tendsto_cosh27 (y : ℝ) :
    Tendsto (fun n : ℕ ↦ zudilinHalfEulerProduct27 n y) atTop
      (𝓝 (Real.cosh (Real.pi * y))) := by
  by_cases hy : y = 0
  · subst y
    simp [zudilinHalfEulerProduct27]
  let z : ℂ := (y : ℂ) * Complex.I
  have hz : z ≠ 0 := by
    exact mul_ne_zero (Complex.ofReal_ne_zero.mpr hy) Complex.I_ne_zero
  have hsin : Complex.sin ((Real.pi : ℂ) * z) ≠ 0 := by
    rw [show (Real.pi : ℂ) * z = ((Real.pi * y : ℝ) : ℂ) * Complex.I by
      simp [z]; ring]
    rw [Complex.sin_mul_I]
    exact mul_ne_zero (by
      rw [← Complex.ofReal_sinh]
      simpa only [Complex.ofReal_ne_zero] using
        (Real.sinh_ne_zero.mpr (mul_ne_zero Real.pi_ne_zero hy)))
      Complex.I_ne_zero
  have hnum := Complex.tendsto_euler_sin_prod (2 * z)
  have hnum' := hnum.comp
    (tendsto_id.const_mul_atTop' (by norm_num : 0 < (2 : ℕ)))
  have hden := Complex.tendsto_euler_sin_prod z
  have hquot :
      Tendsto
        (fun n : ℕ ↦
          ((Real.pi : ℂ) * (2 * z) *
              ∏ j ∈ Finset.range (2 * n),
                (1 - (2 * z) ^ 2 / ((j : ℂ) + 1) ^ 2)) /
            (2 * ((Real.pi : ℂ) * z *
              ∏ k ∈ Finset.range n,
                (1 - z ^ 2 / ((k : ℂ) + 1) ^ 2))))
        atTop
        (𝓝 (Complex.sin ((Real.pi : ℂ) * (2 * z)) /
          (2 * Complex.sin ((Real.pi : ℂ) * z)))) := by
    apply Filter.Tendsto.div hnum'
      (tendsto_const_nhds.mul hden)
    exact mul_ne_zero (by norm_num) hsin
  have hodd :
      Tendsto (fun n : ℕ ↦ zudilinOddEulerProduct27 n z) atTop
        (𝓝 (Complex.sin ((Real.pi : ℂ) * (2 * z)) /
          (2 * Complex.sin ((Real.pi : ℂ) * z)))) := by
    apply hquot.congr'
    filter_upwards [] with n
    rw [zudilinEulerProduct_even_odd27]
    have hbase :
        (Real.pi : ℂ) * z *
            ∏ k ∈ Finset.range n,
              (1 - z ^ 2 / ((k : ℂ) + 1) ^ 2) ≠ 0 := by
      apply mul_ne_zero
        (mul_ne_zero (Complex.ofReal_ne_zero.mpr Real.pi_ne_zero) hz)
      apply Finset.prod_ne_zero_iff.mpr
      intro k hk
      have hfactor :
          1 - z ^ 2 / ((k : ℂ) + 1) ^ 2 =
            ((1 + y ^ 2 / ((k : ℝ) + 1) ^ 2 : ℝ) : ℂ) := by
        dsimp only [z]
        push_cast
        field_simp [Nat.cast_add_one_ne_zero (R := ℂ) k]
        rw [Complex.I_sq]
        ring
      rw [hfactor]
      exact Complex.ofReal_ne_zero.mpr (by positivity)
    have hrewrite :
        (Real.pi : ℂ) * (2 * z) *
            ((∏ k ∈ Finset.range n,
              (1 - z ^ 2 / ((k : ℂ) + 1) ^ 2)) *
                zudilinOddEulerProduct27 n z) =
          2 * ((Real.pi : ℂ) * z *
            ∏ k ∈ Finset.range n,
              (1 - z ^ 2 / ((k : ℂ) + 1) ^ 2)) *
                zudilinOddEulerProduct27 n z := by
      ring
    rw [hrewrite]
    simp [hbase]
  have hlimit :
      Complex.sin ((Real.pi : ℂ) * (2 * z)) /
          (2 * Complex.sin ((Real.pi : ℂ) * z)) =
        (Real.cosh (Real.pi * y) : ℂ) := by
    rw [show (Real.pi : ℂ) * (2 * z) =
      2 * ((Real.pi : ℂ) * z) by ring, Complex.sin_two_mul]
    field_simp [hsin]
    rw [show (Real.pi : ℂ) * z =
      ((Real.pi * y : ℝ) : ℂ) * Complex.I by simp [z]; ring]
    rw [Complex.cos_mul_I]
    simp
  rw [hlimit] at hodd
  have hcast :
      Tendsto (fun n : ℕ ↦ (zudilinHalfEulerProduct27 n y : ℂ)) atTop
        (𝓝 (Real.cosh (Real.pi * y) : ℂ)) := by
    simpa only [z, zudilinOddEulerProduct_eq_half27] using hodd
  convert (Complex.continuous_re.tendsto _).comp hcast using 1 <;> simp

theorem zudilinHalfEulerProduct_monotone27 (y : ℝ) :
    Monotone (fun n : ℕ ↦ zudilinHalfEulerProduct27 n y) := by
  intro n m hnm
  unfold zudilinHalfEulerProduct27
  apply Finset.prod_le_prod_of_subset_of_one_le (Finset.range_mono hnm)
  · intro k hk
    have hq : 0 ≤ y ^ 2 / ((k : ℝ) + 1 / 2) ^ 2 := by positivity
    linarith
  · intro k hk hkn
    have hq : 0 ≤ y ^ 2 / ((k : ℝ) + 1 / 2) ^ 2 := by positivity
    linarith

theorem zudilinHalfEulerProduct_le_cosh27 (n : ℕ) (y : ℝ) :
    zudilinHalfEulerProduct27 n y ≤ Real.cosh (Real.pi * y) :=
  (zudilinHalfEulerProduct_monotone27 y).ge_of_tendsto
    (zudilinHalfEulerProduct_tendsto_cosh27 y) n

/-! ## A uniform majorant for the fixed-line rational function -/

private theorem norm_zudilinBarnesNumeratorFactor_eq27 (k : ℕ) (y : ℝ) :
    ‖zudilinBarnesNumeratorFactor27 k (zudilinBarnesLine27 y)‖ =
      zudilinHalfFactor27 k *
        Real.sqrt (1 + y ^ 2 / zudilinHalfFactor27 k ^ 2) := by
  rw [Complex.norm_def]
  simp only [zudilinBarnesNumeratorFactor27, zudilinBarnesLine27,
    zudilinHalfFactor27, Complex.normSq_apply, Complex.add_re, Complex.add_im,
    Complex.mul_re, Complex.mul_im, Complex.ofReal_re, Complex.ofReal_im,
    Complex.I_re, Complex.I_im]
  norm_num
  have hk : 0 < (k : ℝ) + 1 / 2 := by positivity
  rw [show ((-(1 / 2 : ℝ) + ((k : ℝ) + 1)) *
      (-(1 / 2 : ℝ) + ((k : ℝ) + 1)) + y * y) =
    (((k : ℝ) + 1 / 2) ^ 2) *
      (1 + y ^ 2 / ((k : ℝ) + 1 / 2) ^ 2) by
      field_simp [ne_of_gt hk]; ring]
  rw [Real.sqrt_mul (sq_nonneg _), Real.sqrt_sq hk.le]

private theorem zudilinHalfSqrtProduct_eq27 (n : ℕ) (y : ℝ) :
    (∏ k ∈ Finset.range n,
      Real.sqrt (1 + y ^ 2 / zudilinHalfFactor27 k ^ 2)) =
      Real.sqrt (zudilinHalfEulerProduct27 n y) := by
  rw [zudilinHalfEulerProduct27, ← Real.sqrt_prod]
  · unfold zudilinHalfFactor27
    rfl
  · intro k hk
    positivity

private theorem zudilinHalfFactor_le_norm_denominator27
    (n k : ℕ) (y : ℝ) :
    zudilinHalfFactor27 (n + k) ≤
      ‖zudilinBarnesDenominatorFactor27 n k (zudilinBarnesLine27 y)‖ := by
  rw [Complex.norm_def]
  simp only [zudilinBarnesDenominatorFactor27, zudilinBarnesLine27,
    zudilinHalfFactor27, Complex.normSq_apply, Complex.add_re, Complex.add_im,
    Complex.mul_re, Complex.mul_im, Complex.ofReal_re, Complex.ofReal_im,
    Complex.I_re, Complex.I_im]
  norm_num
  have h : 0 ≤ (n : ℝ) + k + 1 / 2 := by positivity
  apply (Real.le_sqrt h (by
    nlinarith [sq_nonneg (-(1 / 2 : ℝ) + ((n : ℝ) + k + 1)),
      sq_nonneg y])).2
  nlinarith [sq_nonneg y]

theorem norm_zudilinBarnesF_le_midpoint_cosh27 (n : ℕ) (y : ℝ) :
    ‖zudilinBarnesF27 n (zudilinBarnesLine27 y)‖ ≤
      zudilinMidpointValue27 n *
        (Real.cosh (Real.pi * y) *
          Real.sqrt (Real.cosh (Real.pi * y))) := by
  let A : ℝ := ∏ k ∈ Finset.range n, zudilinHalfFactor27 k
  let D : ℝ :=
    ∏ k ∈ Finset.range (n + 1), zudilinHalfFactor27 (n + k)
  let S : ℝ := ∏ k ∈ Finset.range n,
    Real.sqrt (1 + y ^ 2 / zudilinHalfFactor27 k ^ 2)
  have hA : 0 ≤ A :=
    Finset.prod_nonneg fun k hk ↦ (zudilinHalfFactor_pos27 k).le
  have hD : 0 < D :=
    Finset.prod_pos fun k hk ↦ zudilinHalfFactor_pos27 (n + k)
  have hS : 0 ≤ S :=
    Finset.prod_nonneg fun k hk ↦ Real.sqrt_nonneg _
  have hnum :
      ‖∏ k ∈ Finset.range n,
          zudilinBarnesNumeratorFactor27 k (zudilinBarnesLine27 y)‖ =
        A * S := by
    rw [norm_prod]
    simp_rw [norm_zudilinBarnesNumeratorFactor_eq27]
    rw [Finset.prod_mul_distrib]
  have hden :
      ((n.factorial : ℝ) ^ 2) * D ≤
        ‖((n.factorial : ℂ) ^ 2) *
          ∏ k ∈ Finset.range (n + 1),
            zudilinBarnesDenominatorFactor27 n k
              (zudilinBarnesLine27 y)‖ := by
    rw [norm_mul, norm_pow, Complex.norm_natCast, norm_prod]
    apply mul_le_mul_of_nonneg_left _ (sq_nonneg _)
    exact Finset.prod_le_prod
      (fun k hk ↦ (zudilinHalfFactor_pos27 (n + k)).le)
      fun k hk ↦ zudilinHalfFactor_le_norm_denominator27 n k y
  have hdenpos : 0 < ((n.factorial : ℝ) ^ 2) * D :=
    mul_pos (sq_pos_of_pos (Nat.cast_pos.mpr (Nat.factorial_pos n))) hD
  rw [zudilinBarnesF27, norm_div, norm_pow, hnum]
  calc
    (A * S) ^ 3 /
        ‖((n.factorial : ℂ) ^ 2) *
          ∏ k ∈ Finset.range (n + 1),
            zudilinBarnesDenominatorFactor27 n k
              (zudilinBarnesLine27 y)‖ ≤
      (A * S) ^ 3 / (((n.factorial : ℝ) ^ 2) * D) := by
        exact div_le_div_of_nonneg_left
          (pow_nonneg (mul_nonneg hA hS) _) hdenpos hden
    _ = zudilinMidpointDirect27 n * S ^ 3 := by
      rw [zudilinMidpointDirect27]
      change (A * S) ^ 3 / (((n.factorial : ℝ) ^ 2) * D) =
        (A ^ 3 / (((n.factorial : ℝ) ^ 2) * D)) * S ^ 3
      ring
    _ = zudilinMidpointValue27 n * S ^ 3 := by
      rw [zudilinMidpointDirect_eq_product27,
        zudilinMidpointProduct_eq_value27]
    _ ≤ zudilinMidpointValue27 n *
        (Real.cosh (Real.pi * y) *
          Real.sqrt (Real.cosh (Real.pi * y))) := by
      apply mul_le_mul_of_nonneg_left _ (zudilinMidpointValue_pos27 n).le
      rw [show S = Real.sqrt (zudilinHalfEulerProduct27 n y) by
        exact zudilinHalfSqrtProduct_eq27 n y]
      have hP : 0 ≤ zudilinHalfEulerProduct27 n y := by
        unfold zudilinHalfEulerProduct27
        positivity
      have hC : 0 ≤ Real.cosh (Real.pi * y) := (Real.cosh_pos _).le
      have hPC := zudilinHalfEulerProduct_le_cosh27 n y
      change (Real.sqrt (zudilinHalfEulerProduct27 n y)) ^ 2 *
          Real.sqrt (zudilinHalfEulerProduct27 n y) ≤ _
      rw [Real.sq_sqrt hP]
      exact mul_le_mul hPC (Real.sqrt_le_sqrt hPC)
        (Real.sqrt_nonneg _) hC

/-! ## Integrability and decay of the Barnes integral -/

def zudilinBarnesEnvelope27 (y : ℝ) : ℝ :=
  1 / Real.sqrt (Real.cosh (Real.pi * y))

private theorem cosh_lower_exp_abs27 (x : ℝ) :
    Real.exp |x| / 2 ≤ Real.cosh x := by
  rw [Real.cosh_eq]
  rcases le_total 0 x with hx | hx
  · rw [abs_of_nonneg hx]
    linarith [Real.exp_pos (-x)]
  · rw [abs_of_nonpos hx]
    linarith [Real.exp_pos x]

theorem zudilinBarnesEnvelope_le_exp27 (y : ℝ) :
    zudilinBarnesEnvelope27 y ≤
      2 * Real.exp (-(Real.pi / 2) * |y|) := by
  let a : ℝ := |Real.pi * y|
  have ha : 0 ≤ a := abs_nonneg _
  have hC : Real.exp a / 2 ≤ Real.cosh (Real.pi * y) :=
    cosh_lower_exp_abs27 (Real.pi * y)
  have hCpos : 0 < Real.cosh (Real.pi * y) := Real.cosh_pos _
  have hsq :
      (Real.exp (a / 2) / 2) ^ 2 ≤ Real.cosh (Real.pi * y) := by
    calc
      (Real.exp (a / 2) / 2) ^ 2 = Real.exp a / 4 := by
        rw [div_pow, pow_two (Real.exp (a / 2)), ← Real.exp_add]
        norm_num
      _ ≤ Real.exp a / 2 := by nlinarith [Real.exp_pos a]
      _ ≤ Real.cosh (Real.pi * y) := hC
  have hsqrt :
      Real.exp (a / 2) / 2 ≤ Real.sqrt (Real.cosh (Real.pi * y)) :=
    (Real.le_sqrt (by positivity) hCpos.le).2 hsq
  have hinv := one_div_le_one_div_of_le
    (by positivity : 0 < Real.exp (a / 2) / 2) hsqrt
  unfold zudilinBarnesEnvelope27
  calc
    1 / Real.sqrt (Real.cosh (Real.pi * y)) ≤
        1 / (Real.exp (a / 2) / 2) := hinv
    _ = 2 * Real.exp (-a / 2) := by
      rw [show -a / 2 = -(a / 2) by ring, Real.exp_neg]
      field_simp [ne_of_gt (Real.exp_pos (a / 2))]
    _ = 2 * Real.exp (-(Real.pi / 2) * |y|) := by
      congr 2
      simp only [a, abs_mul, abs_of_pos Real.pi_pos]
      ring

private theorem integrable_zudilinBarnesExpEnvelope27 :
    MeasureTheory.Integrable
      (fun y : ℝ ↦ 2 * Real.exp (-(Real.pi / 2) * |y|)) := by
  rw [← MeasureTheory.integrableOn_univ,
    ← Set.Iic_union_Ioi (a := (0 : ℝ)),
    MeasureTheory.integrableOn_union]
  constructor
  · have h := integrableOn_exp_mul_Iic
      (a := Real.pi / 2) (by positivity) 0
    have h' := h.const_mul 2
    apply h'.congr
    filter_upwards [MeasureTheory.ae_restrict_mem measurableSet_Iic] with y hy
    rw [abs_of_nonpos hy]
    congr 2
    ring
  · have h := integrableOn_exp_mul_Ioi
      (a := -(Real.pi / 2)) (by linarith [Real.pi_pos]) 0
    have h' := h.const_mul 2
    apply h'.congr
    filter_upwards [MeasureTheory.ae_restrict_mem measurableSet_Ioi] with y hy
    rw [abs_of_pos hy]

theorem integrable_zudilinBarnesEnvelope27 :
    MeasureTheory.Integrable zudilinBarnesEnvelope27 := by
  apply MeasureTheory.Integrable.mono' integrable_zudilinBarnesExpEnvelope27
  · apply Continuous.aestronglyMeasurable
    unfold zudilinBarnesEnvelope27
    apply Continuous.div continuous_const
      (Real.continuous_sqrt.comp
        (Real.continuous_cosh.comp (continuous_const.mul continuous_id)))
    intro y
    exact ne_of_gt (Real.sqrt_pos.2 (Real.cosh_pos _))
  · filter_upwards with y
    rw [Real.norm_eq_abs,
      abs_of_nonneg (by unfold zudilinBarnesEnvelope27; positivity)]
    exact zudilinBarnesEnvelope_le_exp27 y

def zudilinBarnesCoefficient27 (n : ℕ) : ℝ :=
  (4 * n + 2) * zudilinMidpointValue27 n

theorem zudilinBarnesCoefficient_tendsto_zero27 :
    Tendsto zudilinBarnesCoefficient27 atTop (𝓝 0) := by
  have hnonneg : ∀ n, 0 ≤ zudilinBarnesCoefficient27 n := by
    intro n
    unfold zudilinBarnesCoefficient27
    exact mul_nonneg (by positivity) (zudilinMidpointValue_pos27 n).le
  have hbound : ∀ n, zudilinBarnesCoefficient27 n ≤
      (8 * (n : ℝ) + 4) * (1 / 4 : ℝ) ^ n := by
    intro n
    unfold zudilinBarnesCoefficient27
    calc
      (4 * n + 2) * zudilinMidpointValue27 n ≤
          (4 * n + 2) * (2 * (1 / 4 : ℝ) ^ n) := by
        exact mul_le_mul_of_nonneg_left (zudilinMidpointValue_le27 n)
          (by positivity)
      _ = (8 * (n : ℝ) + 4) * (1 / 4 : ℝ) ^ n := by
        ring
  have hmajor :
      Tendsto (fun n : ℕ ↦ (8 * (n : ℝ) + 4) * (1 / 4 : ℝ) ^ n)
        atTop (𝓝 0) := by
    have hlin := tendsto_self_mul_const_pow_of_lt_one
      (r := (1 / 4 : ℝ)) (by norm_num) (by norm_num)
    have hpow := tendsto_pow_atTop_nhds_zero_of_lt_one
      (by norm_num : (0 : ℝ) ≤ 1 / 4)
      (by norm_num : (1 / 4 : ℝ) < 1)
    convert (hlin.const_mul 8).add (hpow.const_mul 4) using 1 <;> ring
  exact squeeze_zero' (Eventually.of_forall hnonneg)
    (Eventually.of_forall hbound) hmajor

/-- The real-line parameterization of the combined Barnes integrand. -/
def zudilinBarnesRealIntegrand27 (n : ℕ) (y : ℝ) : ℂ :=
  zudilinBarnesPhi27 n (zudilinBarnesLine27 y) /
    (Real.cosh (Real.pi * y) : ℂ) ^ 2

def zudilinBarnesSquaredSineKernel27 (s : ℂ) : ℂ :=
  ((Real.pi : ℂ) / Complex.sin ((Real.pi : ℂ) * s)) ^ 2

theorem zudilinBarnesSin_line27 (y : ℝ) :
    Complex.sin ((Real.pi : ℂ) * zudilinBarnesLine27 y) =
      -(Real.cosh (Real.pi * y) : ℂ) := by
  rw [show (Real.pi : ℂ) * zudilinBarnesLine27 y =
      -((Real.pi : ℂ) / 2) +
        ((Real.pi * y : ℝ) : ℂ) * Complex.I by
    simp [zudilinBarnesLine27]
    ring]
  rw [Complex.sin_add, Complex.cos_mul_I, Complex.sin_mul_I]
  simp

theorem zudilinBarnesSquaredSineKernel_line27 (y : ℝ) :
    zudilinBarnesSquaredSineKernel27 (zudilinBarnesLine27 y) =
      (Real.pi : ℂ) ^ 2 / (Real.cosh (Real.pi * y) : ℂ) ^ 2 := by
  rw [zudilinBarnesSquaredSineKernel27, zudilinBarnesSin_line27]
  ring

theorem norm_zudilinBarnesRealIntegrand_le27 (n : ℕ) (y : ℝ) :
    ‖zudilinBarnesRealIntegrand27 n y‖ ≤
      zudilinBarnesCoefficient27 n * zudilinBarnesEnvelope27 y := by
  let C : ℝ := Real.cosh (Real.pi * y)
  have hC : 0 < C := Real.cosh_pos _
  have hphi :
      ‖zudilinBarnesPhi27 n (zudilinBarnesLine27 y)‖ ≤
        zudilinBarnesCoefficient27 n * (C * Real.sqrt C) := by
    calc
      ‖zudilinBarnesPhi27 n (zudilinBarnesLine27 y)‖ ≤
          (4 * n + 2) *
            ‖zudilinBarnesF27 n (zudilinBarnesLine27 y)‖ :=
        norm_zudilinBarnesPhi_le27 n y
      _ ≤ (4 * n + 2) *
          (zudilinMidpointValue27 n * (C * Real.sqrt C)) := by
        apply mul_le_mul_of_nonneg_left
        · exact norm_zudilinBarnesF_le_midpoint_cosh27 n y
        · positivity
      _ = zudilinBarnesCoefficient27 n * (C * Real.sqrt C) := by
        unfold zudilinBarnesCoefficient27 C
        ring
  unfold zudilinBarnesRealIntegrand27 zudilinBarnesEnvelope27
  rw [norm_div, norm_pow, Complex.norm_real, Real.norm_eq_abs, abs_of_pos hC]
  calc
    ‖zudilinBarnesPhi27 n (zudilinBarnesLine27 y)‖ / C ^ 2 ≤
        (zudilinBarnesCoefficient27 n * (C * Real.sqrt C)) / C ^ 2 := by
      exact div_le_div_of_nonneg_right hphi (sq_nonneg C)
    _ = zudilinBarnesCoefficient27 n * (1 / Real.sqrt C) := by
      have hs : Real.sqrt C ≠ 0 := ne_of_gt (Real.sqrt_pos.2 hC)
      field_simp [ne_of_gt hC, hs]
      rw [Real.sq_sqrt hC.le]

/-- The combined, endpoint Barnes integral.  The factor `π/2` is the
Jacobian and normalization obtained from the vertical contour. -/
def zudilinBarnesErrorIntegral27 (n : ℕ) : ℂ :=
  ((Real.pi / 2 : ℝ) : ℂ) *
    ∫ y : ℝ, zudilinBarnesRealIntegrand27 n y

theorem zudilinBarnesErrorIntegral_eq_fixedLine27 (n : ℕ) :
    zudilinBarnesErrorIntegral27 n =
      (1 / (2 * (Real.pi : ℂ))) *
        ∫ y : ℝ,
          zudilinBarnesPhi27 n (zudilinBarnesLine27 y) *
            zudilinBarnesSquaredSineKernel27 (zudilinBarnesLine27 y) := by
  have hintegrand :
      (fun y : ℝ ↦
        zudilinBarnesPhi27 n (zudilinBarnesLine27 y) *
          zudilinBarnesSquaredSineKernel27 (zudilinBarnesLine27 y)) =
      fun y : ℝ ↦ (Real.pi : ℂ) ^ 2 •
        zudilinBarnesRealIntegrand27 n y := by
    funext y
    rw [zudilinBarnesSquaredSineKernel_line27]
    unfold zudilinBarnesRealIntegrand27
    simp only [smul_eq_mul]
    ring
  rw [zudilinBarnesErrorIntegral27, hintegrand,
    MeasureTheory.integral_smul]
  simp only [smul_eq_mul]
  push_cast
  field_simp [Complex.ofReal_ne_zero.mpr Real.pi_ne_zero]

theorem norm_zudilinBarnesErrorIntegral_le27 (n : ℕ) :
    ‖zudilinBarnesErrorIntegral27 n‖ ≤
      zudilinBarnesCoefficient27 n *
        ((Real.pi / 2) * ∫ y : ℝ, zudilinBarnesEnvelope27 y) := by
  have hcoeff : 0 ≤ zudilinBarnesCoefficient27 n := by
    unfold zudilinBarnesCoefficient27
    exact mul_nonneg (by positivity) (zudilinMidpointValue_pos27 n).le
  have hg : MeasureTheory.Integrable
      (fun y : ℝ ↦ zudilinBarnesCoefficient27 n *
        zudilinBarnesEnvelope27 y) :=
    integrable_zudilinBarnesEnvelope27.const_mul _
  have hint :
      ‖∫ y : ℝ, zudilinBarnesRealIntegrand27 n y‖ ≤
        ∫ y : ℝ,
          zudilinBarnesCoefficient27 n * zudilinBarnesEnvelope27 y := by
    exact MeasureTheory.norm_integral_le_of_norm_le hg
      (Filter.Eventually.of_forall
        (norm_zudilinBarnesRealIntegrand_le27 n))
  unfold zudilinBarnesErrorIntegral27
  rw [norm_mul, Complex.norm_real, Real.norm_eq_abs,
    abs_of_pos (by positivity : 0 < Real.pi / 2)]
  calc
    (Real.pi / 2) *
        ‖∫ y : ℝ, zudilinBarnesRealIntegrand27 n y‖ ≤
      (Real.pi / 2) *
        ∫ y : ℝ,
          zudilinBarnesCoefficient27 n * zudilinBarnesEnvelope27 y := by
      exact mul_le_mul_of_nonneg_left hint (by positivity)
    _ = zudilinBarnesCoefficient27 n *
        ((Real.pi / 2) * ∫ y : ℝ, zudilinBarnesEnvelope27 y) := by
      rw [MeasureTheory.integral_const_mul]
      ring

theorem zudilinBarnesErrorIntegral_tendsto_zero27 :
    Tendsto zudilinBarnesErrorIntegral27 atTop (𝓝 0) := by
  rw [tendsto_zero_iff_norm_tendsto_zero]
  let M : ℝ :=
    (Real.pi / 2) * ∫ y : ℝ, zudilinBarnesEnvelope27 y
  have hM : 0 ≤ M := by
    apply mul_nonneg (by positivity)
    exact MeasureTheory.integral_nonneg fun y ↦ by
      unfold zudilinBarnesEnvelope27
      positivity
  have hmajor :
      Tendsto (fun n ↦ zudilinBarnesCoefficient27 n * M) atTop (𝓝 0) := by
    simpa using zudilinBarnesCoefficient_tendsto_zero27.mul_const M
  apply squeeze_zero'
  · exact Eventually.of_forall fun n ↦ norm_nonneg _
  · exact Eventually.of_forall fun n ↦ by
      simpa only [M] using norm_zudilinBarnesErrorIntegral_le27 n
  · exact hmajor

/-! ## The remaining source-normalization interface

All asymptotic analysis is now internal.  The only remaining classical input
is the exact equality between Zudilin's recurrence-defined source error and
the specialized Barnes integral.  This is strictly stronger and more
auditable than assuming the error tends to zero. -/

theorem zudilinCombinedError_tendsto_of_barnes_representation27
    (hBarnes : ∀ n : ℕ,
      (zudilinCombinedError n : ℂ) = zudilinBarnesErrorIntegral27 n) :
    Tendsto zudilinCombinedError atTop (𝓝 0) := by
  have hcomplex :
      Tendsto (fun n ↦ (zudilinCombinedError n : ℂ)) atTop (𝓝 0) := by
    apply zudilinBarnesErrorIntegral_tendsto_zero27.congr'
    exact Eventually.of_forall fun n ↦ (hBarnes n).symm
  convert (Complex.continuous_re.tendsto _).comp hcomplex using 1 <;> simp

theorem problem27_zeta2_add_zeta3_of_barnes_representation
    (hBarnes : ∀ n : ℕ,
      (zudilinCombinedError n : ℂ) = zudilinBarnesErrorIntegral27 n) :
    Tendsto
      (fun n ↦ (challengeP n : ℝ) / (challengeQ n : ℝ))
      atTop (𝓝 (Real.pi ^ 2 / 6 + zeta3)) := by
  apply problem27_zeta2_add_zeta3
  exact zudilinCombinedError_tendsto_of_barnes_representation27 hBarnes

end RamanujanChallenge.P27
