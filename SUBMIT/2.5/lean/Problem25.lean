/-
  Ramanujan Challenge Problem 2.5: the 3×3 CMF for Catalan's constant.

  This file transcribes the challenge's nine matrix coefficients, its initial
  2×3 matrix, and the ordered matrix product.  It proves:

  * the numerator and denominator rows obey the printed CMF recurrence;
  * an exact first-step transcription check;
  * the complete factorization of det M(n);
  * the determinant-compatible Pochhammer gauge identity;
  * the exact factorization and three roots of the associated Poincaré cubic;
  * sign conjugation to a strictly positive cocycle, hence nonzero finite
    denominators and an exact convex-hull nesting law for all three ratios;
  * rigorous rational bounds on Catalan's constant and the signs of all three
    initial transformed Catalan errors;
  * the Catalan-error row obeys the same positive CMF recurrence.

  These are unconditional algebraic facts about the challenge data.  The final
  analytic step identifying the dominant connection coefficient with Catalan's
  constant is not claimed here: the current paper draft supports that step only
  by high-precision numerical evaluation of an adjoint scalar product.
-/
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Analysis.PSeries
import Mathlib.LinearAlgebra.Matrix.Determinant.Basic

noncomputable section

namespace RamanujanChallenge.P25

/-! ## Catalan's constant -/

/-- Catalan's constant in the series form used by the challenge. -/
def catalanConstant : ℝ := ∑' n : ℕ, (-1 : ℝ) ^ n / (2 * (n : ℝ) + 1) ^ 2

/-- The defining Catalan series is absolutely summable. -/
theorem catalanSeries_summable :
    Summable (fun n : ℕ => (-1 : ℝ) ^ n / (2 * (n : ℝ) + 1) ^ 2) := by
  apply Summable.of_norm_bounded (g := fun n : ℕ => (((n + 1 : ℕ) : ℝ) ^ 2)⁻¹)
  · rw [summable_nat_add_iff (f := fun n : ℕ => (((n : ℝ) ^ 2)⁻¹))]
    exact Real.summable_nat_pow_inv.mpr (by norm_num)
  · intro n
    rw [norm_div, norm_pow]
    norm_num [abs_of_nonneg]
    gcongr
    have hn : (0 : ℝ) ≤ n := by positivity
    linarith

/-- The ordinary partial sums converge to `catalanConstant`. -/
theorem catalanPartialSum_tendsto :
    Filter.Tendsto
      (fun N : ℕ => ∑ n ∈ Finset.range N, (-1 : ℝ) ^ n / (2 * (n : ℝ) + 1) ^ 2)
      Filter.atTop (nhds catalanConstant) :=
  catalanSeries_summable.hasSum.tendsto_sum_nat

/-- The positive magnitudes in the alternating Catalan series. -/
def catalanMagnitude (n : ℕ) : ℝ := 1 / (2 * (n : ℝ) + 1) ^ 2

theorem catalanMagnitude_antitone : Antitone catalanMagnitude := by
  apply antitone_nat_of_succ_le
  intro n
  apply one_div_le_one_div_of_le (by positivity)
  gcongr
  omega

theorem catalanMagnitude_partialSum_tendsto :
    Filter.Tendsto
      (fun N : ℕ => ∑ n ∈ Finset.range N, (-1 : ℝ) ^ n * catalanMagnitude n)
      Filter.atTop (nhds catalanConstant) := by
  simpa [catalanMagnitude, div_eq_mul_inv] using catalanPartialSum_tendsto

/-- A rational lower bound, proved from the 42-term alternating partial sum. -/
theorem catalan_lower_bound : (32972 : ℝ) / 36000 < catalanConstant := by
  calc
    (32972 : ℝ) / 36000 <
        ∑ i ∈ Finset.range (2 * 21), (-1 : ℝ) ^ i * catalanMagnitude i := by
      norm_num [catalanMagnitude, Finset.sum_range_succ]
    _ ≤ catalanConstant :=
      catalanMagnitude_antitone.alternating_series_le_tendsto
        catalanMagnitude_partialSum_tendsto 21

/-- A rational upper bound, proved from the 25-term alternating partial sum. -/
theorem catalan_upper_bound : catalanConstant < (30921 : ℝ) / 33750 := by
  calc
    catalanConstant ≤
        ∑ i ∈ Finset.range (2 * 12 + 1), (-1 : ℝ) ^ i * catalanMagnitude i :=
      catalanMagnitude_antitone.tendsto_le_alternating_series
        catalanMagnitude_partialSum_tendsto 12
    _ < (30921 : ℝ) / 33750 := by
      norm_num [catalanMagnitude, Finset.sum_range_succ]

/-! ## The nine printed matrix coefficients -/

def m11 (n : ℤ) : ℤ :=
  (-2 * n - 5) * (n + 3) ^ 2 *
    (136 * n ^ 4 + 1424 * n ^ 3 + 5548 * n ^ 2 + 9551 * n + 6141)

def m12 (n : ℤ) : ℤ :=
  384 * n ^ 6 + 6384 * n ^ 5 + 44168 * n ^ 4 + 162698 * n ^ 3 +
    336377 * n ^ 2 + 369933 * n + 169011

def m13 (n : ℤ) : ℤ :=
  -480 * n ^ 4 - 4980 * n ^ 3 - 19210 * n ^ 2 - 32690 * n - 20730

def m21 (n : ℤ) : ℤ :=
  (n + 2) ^ 2 * (n + 3) ^ 2 * (4 * n + 10) *
    (48 * n ^ 3 + 386 * n ^ 2 + 1017 * n + 879)

def m22 (n : ℤ) : ℤ :=
  (n + 2) ^ 2 *
    (-272 * n ^ 5 - 3848 * n ^ 4 - 21732 * n ^ 3 - 61184 * n ^ 2 -
      85761 * n - 47808)

def m23 (n : ℤ) : ℤ :=
  (n + 2) ^ 2 * (320 * n ^ 3 + 2540 * n ^ 2 + 6610 * n + 5640)

def m31 (n : ℤ) : ℤ :=
  (-4 * n - 10) * (n + 2) ^ 2 * (n + 3) ^ 2 *
    (32 * n ^ 4 + 302 * n ^ 3 + 1037 * n ^ 2 + 1530 * n + 813)

def m32 (n : ℤ) : ℤ :=
  (n + 2) ^ 2 *
    (192 * n ^ 6 + 2984 * n ^ 5 + 19116 * n ^ 4 + 64452 * n ^ 3 +
      120256 * n ^ 2 + 117279 * n + 46476)

def m33 (n : ℤ) : ℤ :=
  (n + 2) ^ 2 *
    (-16 * n ^ 5 - 408 * n ^ 4 - 2912 * n ^ 3 - 8884 * n ^ 2 -
      12254 * n - 6240)

/-- The coefficient matrix `M(n)` printed in Problem 2.5. -/
def challengeMatrix (n : ℤ) : Matrix (Fin 3) (Fin 3) ℤ :=
  !![m11 n, m12 n, m13 n;
     m21 n, m22 n, m23 n;
     m31 n, m32 n, m33 n]

/-- The challenge's initial 2×3 matrix. -/
def initialMatrix : Matrix (Fin 2) (Fin 3) ℤ :=
  !![30921, -32972, 8240;
     33750, -36000, 9000]

/-! ## The ordered CMF product and challenge sequences -/

/-- `matrixProduct N = M(0) M(1) ... M(N-1)`. -/
def matrixProduct : ℕ → Matrix (Fin 3) (Fin 3) ℤ
  | 0 => 1
  | n + 1 => matrixProduct n * challengeMatrix (n : ℤ)

/-- The 2×3 matrix `A M(0) ... M(N-1)` from the challenge. -/
def approximants (N : ℕ) : Matrix (Fin 2) (Fin 3) ℤ :=
  initialMatrix * matrixProduct N

/-- The numerator in column `j` after `N` steps. -/
def numerator (N : ℕ) (j : Fin 3) : ℤ := approximants N 0 j

/-- The denominator in column `j` after `N` steps. -/
def denominator (N : ℕ) (j : Fin 3) : ℤ := approximants N 1 j

@[simp] theorem matrixProduct_zero : matrixProduct 0 = 1 := rfl

@[simp] theorem matrixProduct_succ (n : ℕ) :
    matrixProduct (n + 1) = matrixProduct n * challengeMatrix (n : ℤ) := rfl

@[simp] theorem approximants_zero : approximants 0 = initialMatrix := by
  simp [approximants]

/-- The literal right-multiplication recurrence from the challenge. -/
theorem approximants_succ (n : ℕ) :
    approximants (n + 1) = approximants n * challengeMatrix (n : ℤ) := by
  simp [approximants, Matrix.mul_assoc]

/-- The numerator row satisfies the challenge CMF recurrence. -/
theorem numerator_succ (n : ℕ) (j : Fin 3) :
    numerator (n + 1) j =
      ∑ i : Fin 3, numerator n i * challengeMatrix (n : ℤ) i j := by
  change approximants (n + 1) 0 j = _
  rw [approximants_succ, Matrix.mul_apply]
  rfl

/-- The denominator row satisfies the challenge CMF recurrence. -/
theorem denominator_succ (n : ℕ) (j : Fin 3) :
    denominator (n + 1) j =
      ∑ i : Fin 3, denominator n i * challengeMatrix (n : ℤ) i j := by
  change approximants (n + 1) 1 j = _
  rw [approximants_succ, Matrix.mul_apply]
  rfl

theorem numerator_succ_zero (n : ℕ) :
    numerator (n + 1) 0 =
      numerator n 0 * m11 n + numerator n 1 * m21 n + numerator n 2 * m31 n := by
  rw [numerator_succ]
  simp [challengeMatrix, Fin.sum_univ_succ]
  ring

theorem numerator_succ_one (n : ℕ) :
    numerator (n + 1) 1 =
      numerator n 0 * m12 n + numerator n 1 * m22 n + numerator n 2 * m32 n := by
  rw [numerator_succ]
  simp [challengeMatrix, Fin.sum_univ_succ]
  ring

theorem numerator_succ_two (n : ℕ) :
    numerator (n + 1) 2 =
      numerator n 0 * m13 n + numerator n 1 * m23 n + numerator n 2 * m33 n := by
  rw [numerator_succ]
  simp [challengeMatrix, Fin.sum_univ_succ]
  ring

theorem denominator_succ_zero (n : ℕ) :
    denominator (n + 1) 0 =
      denominator n 0 * m11 n + denominator n 1 * m21 n + denominator n 2 * m31 n := by
  rw [denominator_succ]
  simp [challengeMatrix, Fin.sum_univ_succ]
  ring

theorem denominator_succ_one (n : ℕ) :
    denominator (n + 1) 1 =
      denominator n 0 * m12 n + denominator n 1 * m22 n + denominator n 2 * m32 n := by
  rw [denominator_succ]
  simp [challengeMatrix, Fin.sum_univ_succ]
  ring

theorem denominator_succ_two (n : ℕ) :
    denominator (n + 1) 2 =
      denominator n 0 * m13 n + denominator n 1 * m23 n + denominator n 2 * m33 n := by
  rw [denominator_succ]
  simp [challengeMatrix, Fin.sum_univ_succ]
  ring

/-! ## Sign conjugation to a strictly positive cocycle -/

/-- The alternating coordinate signs occurring in every row and column of the
printed matrix. -/
def coordinateSign : Fin 3 → ℤ := ![1, -1, 1]

/-- After conjugating by `diag(1,-1,1)` and changing the overall sign, every
entry of the challenge matrix is positive at a natural-number index. -/
def positiveMatrix (n : ℤ) : Matrix (Fin 3) (Fin 3) ℤ :=
  !![(2 * n + 5) * (n + 3) ^ 2 *
        (136 * n ^ 4 + 1424 * n ^ 3 + 5548 * n ^ 2 + 9551 * n + 6141),
      384 * n ^ 6 + 6384 * n ^ 5 + 44168 * n ^ 4 + 162698 * n ^ 3 +
        336377 * n ^ 2 + 369933 * n + 169011,
      480 * n ^ 4 + 4980 * n ^ 3 + 19210 * n ^ 2 + 32690 * n + 20730;
     (n + 2) ^ 2 * (n + 3) ^ 2 * (4 * n + 10) *
        (48 * n ^ 3 + 386 * n ^ 2 + 1017 * n + 879),
      (n + 2) ^ 2 *
        (272 * n ^ 5 + 3848 * n ^ 4 + 21732 * n ^ 3 + 61184 * n ^ 2 +
          85761 * n + 47808),
      (n + 2) ^ 2 * (320 * n ^ 3 + 2540 * n ^ 2 + 6610 * n + 5640);
     (4 * n + 10) * (n + 2) ^ 2 * (n + 3) ^ 2 *
        (32 * n ^ 4 + 302 * n ^ 3 + 1037 * n ^ 2 + 1530 * n + 813),
      (n + 2) ^ 2 *
        (192 * n ^ 6 + 2984 * n ^ 5 + 19116 * n ^ 4 + 64452 * n ^ 3 +
          120256 * n ^ 2 + 117279 * n + 46476),
      (n + 2) ^ 2 *
        (16 * n ^ 5 + 408 * n ^ 4 + 2912 * n ^ 3 + 8884 * n ^ 2 +
          12254 * n + 6240)]

/-- Entrywise form of
`positiveMatrix n = -diag(1,-1,1) * challengeMatrix n * diag(1,-1,1)`.
Stating the identity entrywise avoids introducing a second matrix product into
later row recurrences. -/
theorem positiveMatrix_eq_sign_conjugate (n : ℤ) (i j : Fin 3) :
    positiveMatrix n i j =
      -coordinateSign i * challengeMatrix n i j * coordinateSign j := by
  fin_cases i <;> fin_cases j <;>
    norm_num [positiveMatrix, coordinateSign, challengeMatrix,
      m11, m12, m13, m21, m22, m23, m31, m32, m33] <;> ring

/-- All nine entries of the sign-conjugated transition are strictly positive. -/
theorem positiveMatrix_pos (n : ℕ) (i j : Fin 3) :
    0 < positiveMatrix (n : ℤ) i j := by
  fin_cases i <;> fin_cases j <;>
    simp [positiveMatrix]
  all_goals positivity

/-- The numerator after the simultaneous parity and coordinate sign change. -/
def positiveNumerator (N : ℕ) (j : Fin 3) : ℤ :=
  (-1 : ℤ) ^ N * numerator N j * coordinateSign j

/-- The denominator after the simultaneous parity and coordinate sign change. -/
def positiveDenominator (N : ℕ) (j : Fin 3) : ℤ :=
  (-1 : ℤ) ^ N * denominator N j * coordinateSign j

/-- The transformed numerator evolves by the strictly positive cocycle. -/
theorem positiveNumerator_succ (n : ℕ) (j : Fin 3) :
    positiveNumerator (n + 1) j =
      ∑ i : Fin 3, positiveNumerator n i * positiveMatrix (n : ℤ) i j := by
  rw [positiveNumerator, numerator_succ]
  simp only [positiveNumerator, positiveMatrix_eq_sign_conjugate]
  simp [coordinateSign, Fin.sum_univ_succ, pow_succ]
  ring

/-- The transformed denominator evolves by the strictly positive cocycle. -/
theorem positiveDenominator_succ (n : ℕ) (j : Fin 3) :
    positiveDenominator (n + 1) j =
      ∑ i : Fin 3, positiveDenominator n i * positiveMatrix (n : ℤ) i j := by
  rw [positiveDenominator, denominator_succ]
  simp only [positiveDenominator, positiveMatrix_eq_sign_conjugate]
  simp [coordinateSign, Fin.sum_univ_succ, pow_succ]
  ring

theorem positiveNumerator_zero :
    (fun j => positiveNumerator 0 j) = ![30921, 32972, 8240] := by
  funext j
  fin_cases j <;>
    norm_num [positiveNumerator, numerator, approximants, initialMatrix, coordinateSign]

theorem positiveDenominator_zero :
    (fun j => positiveDenominator 0 j) = ![33750, 36000, 9000] := by
  funext j
  fin_cases j <;>
    norm_num [positiveDenominator, denominator, approximants, initialMatrix, coordinateSign]

/-- Every transformed numerator entry stays strictly positive. -/
theorem positiveNumerator_pos : ∀ N : ℕ, ∀ j : Fin 3, 0 < positiveNumerator N j
  | 0, j => by
      fin_cases j <;>
        norm_num [positiveNumerator, numerator, approximants, initialMatrix, coordinateSign]
  | n + 1, j => by
      rw [positiveNumerator_succ]
      exact Finset.sum_pos
        (fun i _ => mul_pos (positiveNumerator_pos n i) (positiveMatrix_pos n i j))
        Finset.univ_nonempty

/-- Every transformed denominator entry stays strictly positive. -/
theorem positiveDenominator_pos : ∀ N : ℕ, ∀ j : Fin 3, 0 < positiveDenominator N j
  | 0, j => by
      fin_cases j <;>
        norm_num [positiveDenominator, denominator, approximants, initialMatrix, coordinateSign]
  | n + 1, j => by
      rw [positiveDenominator_succ]
      exact Finset.sum_pos
        (fun i _ => mul_pos (positiveDenominator_pos n i) (positiveMatrix_pos n i j))
        Finset.univ_nonempty

/-- In particular every denominator occurring in the challenge ratio is
nonzero, at every finite stage and in every column. -/
theorem denominator_ne_zero (N : ℕ) (j : Fin 3) : denominator N j ≠ 0 := by
  intro h
  have hp := positiveDenominator_pos N j
  simp [positiveDenominator, h] at hp

/-! ## Exact convex-hull structure of the rational approximants -/

/-- The challenge ratio, now well-defined at every finite stage by
`denominator_ne_zero`. -/
def challengeRatio (N : ℕ) (j : Fin 3) : ℝ :=
  (numerator N j : ℝ) / (denominator N j : ℝ)

/-- The same ratio expressed using the positive numerator and denominator. -/
def positiveRatio (N : ℕ) (j : Fin 3) : ℝ :=
  (positiveNumerator N j : ℝ) / (positiveDenominator N j : ℝ)

theorem challengeRatio_eq_positiveRatio (N : ℕ) (j : Fin 3) :
    challengeRatio N j = positiveRatio N j := by
  have hq : (denominator N j : ℝ) ≠ 0 := by
    exact_mod_cast denominator_ne_zero N j
  have hs : ((-1 : ℝ) ^ N * (coordinateSign j : ℝ)) ≠ 0 := by
    fin_cases j <;> simp [coordinateSign]
  have hsign : (coordinateSign j : ℝ) ≠ 0 := by
    fin_cases j <;> simp [coordinateSign]
  simp only [challengeRatio, positiveRatio, positiveNumerator, positiveDenominator]
  push_cast
  field_simp [hsign]

/-- Each new column ratio is a weighted average of the three preceding
ratios, with strictly positive weights. -/
theorem positiveRatio_succ_convex (n : ℕ) (j : Fin 3) :
    positiveRatio (n + 1) j =
      (∑ i : Fin 3,
          ((positiveDenominator n i : ℤ) : ℝ) *
            (positiveMatrix (n : ℤ) i j : ℝ) * positiveRatio n i) /
        (∑ i : Fin 3,
          ((positiveDenominator n i : ℤ) : ℝ) *
            (positiveMatrix (n : ℤ) i j : ℝ)) := by
  simp only [positiveRatio]
  rw [positiveNumerator_succ, positiveDenominator_succ]
  push_cast
  congr 1
  apply Finset.sum_congr rfl
  intro i _
  have hq : (positiveDenominator n i : ℝ) ≠ 0 := by
    exact_mod_cast (positiveDenominator_pos n i).ne'
  field_simp

private theorem positiveWeight_pos (n : ℕ) (i j : Fin 3) :
    0 < ((positiveDenominator n i : ℤ) : ℝ) *
      (positiveMatrix (n : ℤ) i j : ℝ) := by
  have hq : (0 : ℝ) < positiveDenominator n i := by
    exact_mod_cast positiveDenominator_pos n i
  have hm : (0 : ℝ) < positiveMatrix (n : ℤ) i j := by
    exact_mod_cast positiveMatrix_pos n i j
  positivity

private theorem positiveWeight_sum_pos (n : ℕ) (j : Fin 3) :
    0 < ∑ i : Fin 3,
      ((positiveDenominator n i : ℤ) : ℝ) *
        (positiveMatrix (n : ℤ) i j : ℝ) :=
  Finset.sum_pos (fun i _ => positiveWeight_pos n i j) Finset.univ_nonempty

/-- A common upper bound on the three ratios is preserved by one challenge
step. -/
theorem positiveRatio_succ_le (n : ℕ) (j : Fin 3) (u : ℝ)
    (h : ∀ i : Fin 3, positiveRatio n i ≤ u) :
    positiveRatio (n + 1) j ≤ u := by
  rw [positiveRatio_succ_convex]
  apply (div_le_iff₀ (positiveWeight_sum_pos n j)).2
  calc
    (∑ i : Fin 3,
        ((positiveDenominator n i : ℤ) : ℝ) *
          (positiveMatrix (n : ℤ) i j : ℝ) * positiveRatio n i) ≤
        ∑ i : Fin 3,
          (((positiveDenominator n i : ℤ) : ℝ) *
            (positiveMatrix (n : ℤ) i j : ℝ)) * u := by
      apply Finset.sum_le_sum
      intro i _
      exact mul_le_mul_of_nonneg_left (h i) (positiveWeight_pos n i j).le
    _ = u * ∑ i : Fin 3,
        ((positiveDenominator n i : ℤ) : ℝ) *
          (positiveMatrix (n : ℤ) i j : ℝ) := by
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro i _
      ring

/-- A common lower bound on the three ratios is preserved by one challenge
step. -/
theorem positiveRatio_succ_ge (n : ℕ) (j : Fin 3) (l : ℝ)
    (h : ∀ i : Fin 3, l ≤ positiveRatio n i) :
    l ≤ positiveRatio (n + 1) j := by
  rw [positiveRatio_succ_convex]
  apply (le_div_iff₀ (positiveWeight_sum_pos n j)).2
  calc
    l * ∑ i : Fin 3,
        ((positiveDenominator n i : ℤ) : ℝ) *
          (positiveMatrix (n : ℤ) i j : ℝ) =
        ∑ i : Fin 3,
          (((positiveDenominator n i : ℤ) : ℝ) *
            (positiveMatrix (n : ℤ) i j : ℝ)) * l := by
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro i _
      ring
    _ ≤ ∑ i : Fin 3,
        ((positiveDenominator n i : ℤ) : ℝ) *
          (positiveMatrix (n : ℤ) i j : ℝ) * positiveRatio n i := by
      apply Finset.sum_le_sum
      intro i _
      exact mul_le_mul_of_nonneg_left (h i) (positiveWeight_pos n i j).le

/-- Every finite challenge ratio remains in the exact convex hull of the three
initial ratios.  This is a genuine nesting result, but by itself does not
identify the limiting point inside that interval. -/
theorem challengeRatio_bounds (N : ℕ) (j : Fin 3) :
    (8240 : ℝ) / 9000 ≤ challengeRatio N j ∧
      challengeRatio N j ≤ (30921 : ℝ) / 33750 := by
  rw [challengeRatio_eq_positiveRatio]
  constructor
  · induction N generalizing j with
    | zero =>
        fin_cases j <;>
          norm_num [positiveRatio, positiveNumerator, positiveDenominator,
            numerator, denominator, approximants, initialMatrix, coordinateSign]
    | succ n ih =>
        apply positiveRatio_succ_ge n j
        intro i
        exact ih i
  · induction N generalizing j with
    | zero =>
        fin_cases j <;>
          norm_num [positiveRatio, positiveNumerator, positiveDenominator,
            numerator, denominator, approximants, initialMatrix, coordinateSign]
    | succ n ih =>
        apply positiveRatio_succ_le n j
        intro i
        exact ih i

/-! ## Transcription checks -/

theorem numerator_zero :
    (fun j => numerator 0 j) = ![30921, -32972, 8240] := by
  funext j
  fin_cases j <;> norm_num [numerator, approximants, initialMatrix]

theorem denominator_zero :
    (fun j => denominator 0 j) = ![33750, -36000, 9000] := by
  funext j
  fin_cases j <;> norm_num [denominator, approximants, initialMatrix]

/-- Exact forward multiplication through `M(0)`, checking all nine coefficient
definitions and the product orientation against the challenge statement. -/
theorem approximants_one :
    approximants 1 =
      !![-21390206625, 13063139595, -1590511050;
         -23352603750, 14261609250, -1736437500] := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    norm_num [approximants, matrixProduct, initialMatrix, challengeMatrix,
      m11, m12, m13, m21, m22, m23, m31, m32, m33,
      Matrix.mul_apply, Fin.sum_univ_succ]

/-! ## Exact determinant and normalization identities -/

/-- The factorization of the determinant of the printed matrix. -/
theorem challengeMatrix_det (n : ℤ) :
    (challengeMatrix n).det =
      -8 * (n + 1) * (n + 2) ^ 6 * (n + 3) ^ 5 *
        (2 * n + 3) ^ 2 * (2 * n + 5) ^ 3 * (2 * n + 7) ^ 4 := by
  rw [Matrix.det_fin_three]
  change
    m11 n * m22 n * m33 n - m11 n * m23 n * m32 n -
        m12 n * m21 n * m33 n + m12 n * m23 n * m31 n +
        m13 n * m21 n * m32 n - m13 n * m22 n * m31 n =
      -8 * (n + 1) * (n + 2) ^ 6 * (n + 3) ^ 5 *
        (2 * n + 3) ^ 2 * (2 * n + 5) ^ 3 * (2 * n + 7) ^ 4
  simp only [m11, m12, m13, m21, m22, m23, m31, m32, m33]
  ring

/-- Every challenge step has nonzero determinant at a natural-number index. -/
theorem challengeMatrix_det_ne_zero (n : ℕ) :
    (challengeMatrix (n : ℤ)).det ≠ 0 := by
  rw [challengeMatrix_det]
  have hp : 0 < ((n : ℤ) + 1) * ((n : ℤ) + 2) ^ 6 * ((n : ℤ) + 3) ^ 5 *
      (2 * (n : ℤ) + 3) ^ 2 * (2 * (n : ℤ) + 5) ^ 3 *
      (2 * (n : ℤ) + 7) ^ 4 := by
    positivity
  have heq :
      -8 * ((n : ℤ) + 1) * ((n : ℤ) + 2) ^ 6 * ((n : ℤ) + 3) ^ 5 *
          (2 * (n : ℤ) + 3) ^ 2 * (2 * (n : ℤ) + 5) ^ 3 *
          (2 * (n : ℤ) + 7) ^ 4 =
        -8 * (((n : ℤ) + 1) * ((n : ℤ) + 2) ^ 6 * ((n : ℤ) + 3) ^ 5 *
          (2 * (n : ℤ) + 3) ^ 2 * (2 * (n : ℤ) + 5) ^ 3 *
          (2 * (n : ℤ) + 7) ^ 4) := by ring
  rw [heq]
  exact mul_ne_zero (by norm_num) hp.ne'

/-- The complete ordered CMF product is nonsingular at every finite stage. -/
theorem matrixProduct_det_ne_zero : ∀ N : ℕ, (matrixProduct N).det ≠ 0
  | 0 => by simp
  | n + 1 => by
      rw [matrixProduct_succ, Matrix.det_mul]
      exact mul_ne_zero (matrixProduct_det_ne_zero n) (challengeMatrix_det_ne_zero n)

/-- The minimal determinant-compatible one-step normalization from the
Pochhammer analysis. -/
def gaugeDelta (n : ℤ) : ℤ :=
  -2 * (n + 2) ^ 2 * (n + 3) ^ 2 * (2 * n + 5) * (2 * n + 7) ^ 2

/-- The determinant correction used to balance the normalized system. -/
def determinantGauge (n : ℤ) : ℤ :=
  (n + 1) * (n + 2) * (2 * n + 3) ^ 2 * (2 * n + 5) ^ 2

/-- Cross-multiplied form of
`det M(n) * h(n+1) / h(n) = gaugeDelta(n)^3`; it has no denominator
or nonvanishing side condition. -/
theorem determinant_gauge_identity (n : ℤ) :
    (challengeMatrix n).det * determinantGauge (n + 1) =
      gaugeDelta n ^ 3 * determinantGauge n := by
  rw [challengeMatrix_det]
  simp only [gaugeDelta, determinantGauge]
  ring

/-! ## The Poincaré cubic -/

/-- The cubic obtained from the leading coefficients of the scalar recurrence
recorded in the challenge investigation. -/
def poincarePolynomial (c : ℝ) : ℝ :=
  c ^ 3 + 560 * c ^ 2 + 8960 * c + 4096

theorem poincare_factor (c : ℝ) :
    poincarePolynomial c = (c + 16) * (c ^ 2 + 544 * c + 256) := by
  simp only [poincarePolynomial]
  ring

theorem poincare_root_neutral : poincarePolynomial (-16) = 0 := by
  norm_num [poincarePolynomial]

theorem poincare_root_plus :
    poincarePolynomial (-16 * (17 + 12 * Real.sqrt 2)) = 0 := by
  have hs : (Real.sqrt 2) ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  simp only [poincarePolynomial]
  nlinarith

theorem poincare_root_minus :
    poincarePolynomial (-16 * (17 - 12 * Real.sqrt 2)) = 0 := by
  have hs : (Real.sqrt 2) ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  simp only [poincarePolynomial]
  nlinarith

/-- The two silver-ratio factors are reciprocal. -/
theorem silver_ratio_reciprocal :
    (17 + 12 * Real.sqrt 2) * (17 - 12 * Real.sqrt 2) = 1 := by
  have hs : (Real.sqrt 2) ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  nlinarith

/-! ## The exact Catalan-error recurrence -/

/-- The error row whose decay would imply the challenge limit. -/
def catalanError (N : ℕ) (j : Fin 3) : ℝ :=
  (denominator N j : ℝ) * catalanConstant - (numerator N j : ℝ)

theorem catalanError_succ_zero (n : ℕ) :
    catalanError (n + 1) 0 =
      catalanError n 0 * m11 n + catalanError n 1 * m21 n +
        catalanError n 2 * m31 n := by
  simp only [catalanError]
  rw [denominator_succ_zero, numerator_succ_zero]
  push_cast
  ring

theorem catalanError_succ_one (n : ℕ) :
    catalanError (n + 1) 1 =
      catalanError n 0 * m12 n + catalanError n 1 * m22 n +
        catalanError n 2 * m32 n := by
  simp only [catalanError]
  rw [denominator_succ_one, numerator_succ_one]
  push_cast
  ring

theorem catalanError_succ_two (n : ℕ) :
    catalanError (n + 1) 2 =
      catalanError n 0 * m13 n + catalanError n 1 * m23 n +
        catalanError n 2 * m33 n := by
  simp only [catalanError]
  rw [denominator_succ_two, numerator_succ_two]
  push_cast
  ring

/-- The Catalan-error row after the same parity and coordinate sign change used
for the numerator and denominator. -/
def positiveCatalanError (N : ℕ) (j : Fin 3) : ℝ :=
  (-1 : ℝ) ^ N * catalanError N j * (coordinateSign j : ℝ)

/-- The transformed error is the same exact linear form in the transformed
denominator and numerator. -/
theorem positiveCatalanError_eq (N : ℕ) (j : Fin 3) :
    positiveCatalanError N j =
      catalanConstant * (positiveDenominator N j : ℝ) -
        (positiveNumerator N j : ℝ) := by
  simp only [positiveCatalanError, catalanError, positiveDenominator,
    positiveNumerator]
  push_cast
  ring

/-- The transformed Catalan error evolves by the positive matrix as well.  This
is the exact starting point for a cone/minimal-solution proof of the missing
connection identity. -/
theorem positiveCatalanError_succ (n : ℕ) (j : Fin 3) :
    positiveCatalanError (n + 1) j =
      ∑ i : Fin 3,
        positiveCatalanError n i * (positiveMatrix (n : ℤ) i j : ℝ) := by
  fin_cases j
  · change positiveCatalanError (n + 1) 0 =
      ∑ i : Fin 3, positiveCatalanError n i * (positiveMatrix (n : ℤ) i 0 : ℝ)
    rw [positiveCatalanError, catalanError_succ_zero]
    simp [positiveCatalanError, positiveMatrix, coordinateSign,
      Fin.sum_univ_succ, pow_succ, m11, m21, m31]
    ring
  · change positiveCatalanError (n + 1) 1 =
      ∑ i : Fin 3, positiveCatalanError n i * (positiveMatrix (n : ℤ) i 1 : ℝ)
    rw [positiveCatalanError, catalanError_succ_one]
    simp [positiveCatalanError, positiveMatrix, coordinateSign,
      Fin.sum_univ_succ, pow_succ, m12, m22, m32]
    ring
  · change positiveCatalanError (n + 1) 2 =
      ∑ i : Fin 3, positiveCatalanError n i * (positiveMatrix (n : ℤ) i 2 : ℝ)
    rw [positiveCatalanError, catalanError_succ_two]
    simp [positiveCatalanError, positiveMatrix, coordinateSign,
      Fin.sum_univ_succ, pow_succ, m13, m23, m33]
    ring

/-- The first transformed error coordinate is rigorously negative. -/
theorem positiveCatalanError_zero_zero_neg : positiveCatalanError 0 0 < 0 := by
  rw [positiveCatalanError_eq]
  have hq := congrFun positiveDenominator_zero (0 : Fin 3)
  have hp := congrFun positiveNumerator_zero (0 : Fin 3)
  norm_num at hq hp
  rw [hq, hp]
  have h := catalan_upper_bound
  norm_num at h
  norm_num at ⊢
  nlinarith

/-- The second transformed error coordinate is rigorously positive. -/
theorem positiveCatalanError_zero_one_pos : 0 < positiveCatalanError 0 1 := by
  rw [positiveCatalanError_eq]
  have hq := congrFun positiveDenominator_zero (1 : Fin 3)
  have hp := congrFun positiveNumerator_zero (1 : Fin 3)
  norm_num at hq hp
  rw [hq, hp]
  have h := catalan_lower_bound
  norm_num at h
  norm_num at ⊢
  nlinarith

/-- The third transformed error coordinate is rigorously positive. -/
theorem positiveCatalanError_zero_two_pos : 0 < positiveCatalanError 0 2 := by
  rw [positiveCatalanError_eq]
  have hq := congrFun positiveDenominator_zero (2 : Fin 3)
  have hp := congrFun positiveNumerator_zero (2 : Fin 3)
  norm_num [Matrix.cons_val_two] at hq hp
  rw [hq, hp]
  have h := catalan_lower_bound
  norm_num at h
  norm_num at ⊢
  nlinarith

/-! ## Faithful main claim and the remaining connection certificate -/

/-- The cumulative Pochhammer gauge forced by the determinant identity. -/
def pochhammerProduct : ℕ → ℤ
  | 0 => 1
  | n + 1 => pochhammerProduct n * gaugeDelta n

@[simp] theorem pochhammerProduct_zero : pochhammerProduct 0 = 1 := rfl

@[simp] theorem pochhammerProduct_succ (n : ℕ) :
    pochhammerProduct (n + 1) = pochhammerProduct n * gaugeDelta n := rfl

/-- None of the factors in the cumulative gauge vanishes at a
natural-number index. -/
theorem pochhammerProduct_ne_zero : ∀ N : ℕ, pochhammerProduct N ≠ 0
  | 0 => by norm_num
  | n + 1 => by
      rw [pochhammerProduct_succ]
      apply mul_ne_zero (pochhammerProduct_ne_zero n)
      have hneg : gaugeDelta (n : ℤ) < 0 := by
        have hpos :
            0 < ((n : ℤ) + 2) ^ 2 * ((n : ℤ) + 3) ^ 2 *
              (2 * (n : ℤ) + 5) * (2 * (n : ℤ) + 7) ^ 2 := by
          positivity
        rw [show gaugeDelta (n : ℤ) =
            -2 * (((n : ℤ) + 2) ^ 2 * ((n : ℤ) + 3) ^ 2 *
              (2 * (n : ℤ) + 5) * (2 * (n : ℤ) + 7) ^ 2) by
          simp only [gaugeDelta]
          ring]
        exact mul_neg_of_neg_of_pos (by norm_num) hpos
      exact ne_of_lt hneg

/-- The explicit dominant scale singled out by the Pochhammer normalization
and the largest root `17 + 12√2` of the normalized Poincaré cubic. -/
def dominantScale (N : ℕ) : ℝ :=
  (pochhammerProduct N : ℝ) * (17 + 12 * Real.sqrt 2) ^ N

theorem dominantScale_ne_zero (N : ℕ) : dominantScale N ≠ 0 := by
  apply mul_ne_zero
  · exact_mod_cast pochhammerProduct_ne_zero N
  · apply pow_ne_zero
    have hs : 0 ≤ Real.sqrt 2 := Real.sqrt_nonneg 2
    nlinarith

/-- The literal assertion of Problem 2.5: all three ratios built from the
printed matrix and printed initial values converge to Catalan's constant. -/
def Problem25Claim : Prop :=
  ∀ j : Fin 3,
    Filter.Tendsto (fun N : ℕ => challengeRatio N j)
      Filter.atTop (nhds catalanConstant)

/-- A precise reduction of the challenge to its missing connection
coefficient.  The first asymptotic hypothesis says that each denominator has
a nonzero component in the dominant Poincaré mode.  The second says that the
same dominant component vanishes for the exact Catalan-error row.  Both are
stated against the fixed, challenge-derived scale `dominantScale`; no
arbitrary witness sequence is introduced.

An integral, hypergeometric gauge, or adjoint-pairing certificate proving
these two limits would therefore close `Problem25Claim`. -/
theorem problem25Claim_of_dominant_connection
    (connection : Fin 3 → ℝ)
    (hconnection : ∀ j : Fin 3, connection j ≠ 0)
    (hdenominator : ∀ j : Fin 3,
      Filter.Tendsto
        (fun N : ℕ => (denominator N j : ℝ) / dominantScale N)
        Filter.atTop (nhds (connection j)))
    (herror : ∀ j : Fin 3,
      Filter.Tendsto
        (fun N : ℕ => catalanError N j / dominantScale N)
        Filter.atTop (nhds 0)) :
    Problem25Claim := by
  intro j
  have hscaled :
      Filter.Tendsto
        (fun N : ℕ =>
          (catalanError N j / dominantScale N) /
            ((denominator N j : ℝ) / dominantScale N))
        Filter.atTop (nhds 0) :=
    by simpa using (herror j).div (hdenominator j) (hconnection j)
  have hquotient :
      Filter.Tendsto
        (fun N : ℕ => catalanError N j / (denominator N j : ℝ))
        Filter.atTop (nhds 0) := by
    apply hscaled.congr'
    filter_upwards [] with N
    have hs : dominantScale N ≠ 0 := dominantScale_ne_zero N
    have hq : (denominator N j : ℝ) ≠ 0 := by
      exact_mod_cast denominator_ne_zero N j
    field_simp [hs, hq]
  have key : ∀ᶠ N : ℕ in Filter.atTop,
      challengeRatio N j =
        catalanConstant - catalanError N j / (denominator N j : ℝ) := by
    filter_upwards [] with N
    have hq : (denominator N j : ℝ) ≠ 0 := by
      exact_mod_cast denominator_ne_zero N j
    simp only [challengeRatio, catalanError]
    field_simp [hq]
    ring
  have htarget :
      Filter.Tendsto
        (fun N : ℕ =>
          catalanConstant - catalanError N j / (denominator N j : ℝ))
        Filter.atTop (nhds catalanConstant) := by
    simpa using
      (Filter.Tendsto.sub
        (tendsto_const_nhds :
          Filter.Tendsto (fun _ : ℕ => catalanConstant)
            Filter.atTop (nhds catalanConstant))
        hquotient)
  exact htarget.congr' (key.mono fun _ h => h.symm)

end RamanujanChallenge.P25

end
