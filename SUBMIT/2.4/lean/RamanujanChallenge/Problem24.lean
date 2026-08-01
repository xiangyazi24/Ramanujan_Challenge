/-
  Ramanujan Challenge Problem 2.4: Polylogarithm Identity

  The double sum identity:
    Σ_{m≥0} Σ_{k=0}^{m} C(m,k)² H_k² / ((m+1)² C(2m,m))
    = 20·Li₄(1/2) + (5/6)·log⁴2 + 10·ζ(2) − (65/9)·ζ(2)²
      − log²2·(12 + 5ζ(2)) + ½·ζ(3) + log2·(35/2·ζ(3) − 16)

  where H_k = Σ_{j=1}^{k} 1/j is the k-th harmonic number.

  The challenge series itself is defined below.  This file proves the exact
  finite inner-sum formula, nonnegativity, and summability, and reduces the
  double series to the corresponding weight-four scalar series.  It does not
  claim the displayed special-value evaluation: the draft HPL argument does
  not contain the certificates needed for that final step.
-/
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecialFunctions.Gamma.Beta
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Analysis.PSeries
import Mathlib.Data.Nat.Choose.Vandermonde
import Mathlib.NumberTheory.Harmonic.Bounds
import Mathlib.NumberTheory.Harmonic.EulerMascheroni
import Mathlib.NumberTheory.ZetaValues
import Mathlib.Tactic
import RamanujanChallenge.Problem26

noncomputable section

open Filter Real Topology

/-! ## Harmonic numbers -/

def harmonicNumber (k : ℕ) : ℝ := ∑ j ∈ Finset.range k, (1 : ℝ) / (↑j + 1)

theorem harmonicNumber_zero : harmonicNumber 0 = 0 := by
  simp [harmonicNumber]

theorem harmonicNumber_one : harmonicNumber 1 = 1 := by
  simp [harmonicNumber]

theorem harmonicNumber_succ (k : ℕ) :
    harmonicNumber (k + 1) = harmonicNumber k + 1 / (k + 1 : ℝ) := by
  simp [harmonicNumber, Finset.sum_range_succ]

theorem harmonicNumber_eq_cast_harmonic (k : ℕ) :
    harmonicNumber k = (harmonic k : ℝ) := by
  simp [harmonicNumber, harmonic, one_div]

theorem harmonicNumber_nonneg (k : ℕ) : 0 ≤ harmonicNumber k := by
  unfold harmonicNumber
  positivity

theorem harmonicNumber_mono {k m : ℕ} (hkm : k ≤ m) :
    harmonicNumber k ≤ harmonicNumber m := by
  unfold harmonicNumber
  exact Finset.sum_le_sum_of_subset_of_nonneg (Finset.range_mono hkm) (by
    intro i _ _
    positivity)

theorem harmonicNumber_two_mul_le (k : ℕ) :
    harmonicNumber (2 * k) ≤ 2 * harmonicNumber k := by
  unfold harmonicNumber
  rw [show 2 * k = k + k by omega, Finset.sum_range_add]
  push_cast
  calc
    (∑ x ∈ Finset.range k, (1 : ℝ) / (x + 1 : ℝ)) +
          ∑ x ∈ Finset.range k, (1 : ℝ) / (k + x + 1 : ℝ) ≤
        (∑ x ∈ Finset.range k, (1 : ℝ) / (x + 1 : ℝ)) +
          ∑ x ∈ Finset.range k, (1 : ℝ) / (x + 1 : ℝ) := by
      gcongr with x hx
      norm_cast
      omega
    _ = 2 * ∑ x ∈ Finset.range k, (1 : ℝ) / (x + 1 : ℝ) := by ring

/-! ## Polylogarithm Li₄ -/

def polylog4 (z : ℝ) : ℝ := ∑' n : ℕ, z ^ (n + 1) / (↑(n + 1) : ℝ) ^ 4

/-! ## The challenge's finite inner sums and infinite outer sum -/

def innerSum24 (m : ℕ) : ℝ := ∑ k ∈ Finset.range (m + 1),
  (Nat.choose m k : ℝ) ^ 2 * (harmonicNumber k) ^ 2

def firstMoment24 (m : ℕ) : ℝ := ∑ k ∈ Finset.range (m + 1),
  (Nat.choose m k : ℝ) ^ 2 * harmonicNumber k

def harmonicRemainder24 (m : ℕ) : ℝ :=
  2 * harmonicNumber m - harmonicNumber (2 * m)

def harmonicRemainderKernel24 (m : ℕ) (t : ℝ) : ℝ :=
  ∑ j ∈ Finset.range m, (t ^ j - t ^ (m + j))

/-- The signed finite harmonic sum `∑_{k=1}^n (-1)^k / k`. -/
def signedHarmonic24 (n : ℕ) : ℝ :=
  ∑ k ∈ Finset.range n, (-1 : ℝ) ^ (k + 1) / (k + 1 : ℝ)

/-- The level-two harmonic combination `H_n + 2∑_{k≤n} (-1)^k/k`.
At an even index `2m` this is exactly `2H_m - H_{2m}`. -/
def parityRemainder24 (n : ℕ) : ℝ :=
  harmonicNumber n + 2 * signedHarmonic24 n

def harmonicSquare24 (n : ℕ) : ℝ :=
  ∑ j ∈ Finset.range n, 1 / (j + 1 : ℝ) ^ 2

def inverseCentralSquareSum24 (m : ℕ) : ℝ :=
  ∑ j ∈ Finset.range m,
    1 / ((j + 1 : ℝ) ^ 2 *
      (Nat.choose (2 * (j + 1)) (j + 1) : ℝ))

def inverseCentralCoefficient24 (j : ℕ) : ℝ :=
  1 / ((j + 1 : ℝ) ^ 2 *
    (Nat.choose (2 * (j + 1)) (j + 1) : ℝ))

def inverseCentralFourthCoefficient24 (j : ℕ) : ℝ :=
  inverseCentralCoefficient24 j / (j + 1 : ℝ) ^ 2

/-- The coefficient of `a²` in the Bailey--Borwein--Bradley generating
identity for the even zeta values, with the index shifted to start at zero. -/
def bbbWeightFourTerm24 (j : ℕ) : ℝ :=
  3 * inverseCentralFourthCoefficient24 j -
    9 * harmonicSquare24 j * inverseCentralCoefficient24 j

/-- The derivative at `u = 0` of the Leshchiner WZ summand.  With
`k = j + 1`, this is
`2/(k⁴ C(2k,k)) - (3/2) H_{k-1}^{(2)}/(k² C(2k,k))`. -/
def leshchinerWeightFourTerm24 (j : ℕ) : ℝ :=
  2 * inverseCentralFourthCoefficient24 j -
    (3 / 2) * harmonicSquare24 j * inverseCentralCoefficient24 j

def inverseCentralDoubleTerm24 (m j : ℕ) : ℝ :=
  if j < m then inverseCentralCoefficient24 j / (m + 1 : ℝ) ^ 2 else 0

def closedInnerSum24 (m : ℕ) : ℝ :=
  (Nat.choose (2 * m) m : ℝ) *
    (harmonicRemainder24 m ^ 2 - harmonicSquare24 (2 * m) +
      3 * inverseCentralSquareSum24 m)

def elementaryOuterTerm24 (m : ℕ) : ℝ :=
  (harmonicRemainder24 m ^ 2 - harmonicSquare24 (2 * m)) /
    (m + 1 : ℝ) ^ 2

/-- The nonalternating quadratic Euler summand, indexed from `n = 1`. -/
def quadraticEulerTerm24 (n : ℕ) : ℝ :=
  (parityRemainder24 (n + 1) ^ 2 - harmonicSquare24 (n + 1)) /
    (n + 1 : ℝ) ^ 2

/-- The outer-alternating companion of `quadraticEulerTerm24`. -/
def alternatingQuadraticEulerTerm24 (n : ℕ) : ℝ :=
  (-1 : ℝ) ^ (n + 1) * quadraticEulerTerm24 n

/-- The nonalternating linear Euler summand of total weight four. -/
def cubicLinearEulerTerm24 (n : ℕ) : ℝ :=
  parityRemainder24 (n + 1) / (n + 1 : ℝ) ^ 3

/-- The outer-alternating companion of `cubicLinearEulerTerm24`. -/
def alternatingCubicLinearEulerTerm24 (n : ℕ) : ℝ :=
  (-1 : ℝ) ^ (n + 1) * cubicLinearEulerTerm24 n

/-- The linear level-two Euler summand of total weight two. -/
def linearEulerTerm24 (n : ℕ) : ℝ :=
  parityRemainder24 (n + 1) / (n + 1 : ℝ)

/-- The outer-alternating companion of `linearEulerTerm24`. -/
def alternatingLinearEulerTerm24 (n : ℕ) : ℝ :=
  (-1 : ℝ) ^ (n + 1) * linearEulerTerm24 n

/-- Consecutive-pair version of the conditionally convergent alternating
weight-two Euler series.  Pairing is part of the definition, so its
`HasSum` statement is an ordinary (indeed absolutely convergent) series
rather than an invalid unconditional rearrangement of the raw terms. -/
def pairedAlternatingLinearEulerTerm24 (m : ℕ) : ℝ :=
  alternatingLinearEulerTerm24 (2 * m) +
    alternatingLinearEulerTerm24 (2 * m + 1)

/-- The ordinary-harmonic part of the preceding paired series. -/
def pairedAlternatingHarmonicEulerTerm24 (m : ℕ) : ℝ :=
  -harmonicNumber (2 * m + 1) / (2 * m + 1 : ℝ) +
    harmonicNumber (2 * m + 2) / (2 * m + 2 : ℝ)

/-- The signed-harmonic part of the preceding paired series. -/
def pairedAlternatingSignedHarmonicEulerTerm24 (m : ℕ) : ℝ :=
  -signedHarmonic24 (2 * m + 1) / (2 * m + 1 : ℝ) +
    signedHarmonic24 (2 * m + 2) / (2 * m + 2 : ℝ)

/-- The linear level-two Euler summand of total weight three. -/
def quadraticLinearEulerTerm24 (n : ℕ) : ℝ :=
  parityRemainder24 (n + 1) / (n + 1 : ℝ) ^ 2

/-- The outer-alternating companion of `quadraticLinearEulerTerm24`. -/
def alternatingQuadraticLinearEulerTerm24 (n : ℕ) : ℝ :=
  (-1 : ℝ) ^ (n + 1) * quadraticLinearEulerTerm24 n

/-- The shifted linear term left by replacing the even index `2m` by
`2(m+1)` in the elementary series. -/
def shiftedLinearEulerTerm24 (m : ℕ) : ℝ :=
  parityRemainder24 (2 * (m + 1)) /
    ((2 * (m + 1) : ℝ) ^ 2 * (2 * m + 1 : ℝ))

/-- The purely rational correction in the same even-index shift. -/
def rationalCorrectionTerm24 (m : ℕ) : ℝ :=
  let n : ℝ := 2 * (m + 1)
  8 / (n ^ 2 * (n - 1) ^ 2) -
    24 / (n ^ 3 * (n - 1)) +
    40 / n ^ 4

def inverseCentralOuterTerm24 (m : ℕ) : ℝ :=
  inverseCentralSquareSum24 m / (m + 1 : ℝ) ^ 2

def outerTerm24 (m : ℕ) : ℝ :=
  innerSum24 m / ((↑m + 1) ^ 2 * (Nat.choose (2 * m) m : ℝ))

def lhs_24 : ℝ := ∑' m : ℕ, outerTerm24 m

/-- A finite-sum form of the harmonic remainder.  It is the termwise
integral of `harmonicRemainderKernel24`; unlike the beta-integral formula in
the draft, it has no exceptional case at `m = 0`. -/
theorem harmonicRemainder24_eq_sum (m : ℕ) :
    harmonicRemainder24 m =
      ∑ j ∈ Finset.range m,
        (1 / (j + 1 : ℝ) - 1 / (m + j + 1 : ℝ)) := by
  unfold harmonicRemainder24 harmonicNumber
  rw [show 2 * m = m + m by omega, Finset.sum_range_add]
  push_cast
  calc
    2 * (∑ j ∈ Finset.range m, 1 / (j + 1 : ℝ)) -
          ((∑ j ∈ Finset.range m, 1 / (j + 1 : ℝ)) +
            ∑ j ∈ Finset.range m, 1 / (m + j + 1 : ℝ)) =
        (∑ j ∈ Finset.range m, 1 / (j + 1 : ℝ)) -
          ∑ j ∈ Finset.range m, 1 / (m + j + 1 : ℝ) := by ring
    _ = ∑ j ∈ Finset.range m,
          (1 / (j + 1 : ℝ) - 1 / (m + j + 1 : ℝ)) := by
      rw [← Finset.sum_sub_distrib]

/-- Away from its removable endpoint, the polynomial kernel is exactly the
quotient used in the informal integral representation. -/
theorem harmonicRemainderKernel24_eq (m : ℕ) {t : ℝ} (ht : t ≠ 1) :
    harmonicRemainderKernel24 m t =
      (1 - t ^ m) ^ 2 / (1 - t) := by
  have hkernel :
      harmonicRemainderKernel24 m t =
        (1 - t ^ m) * ∑ j ∈ Finset.range m, t ^ j := by
    unfold harmonicRemainderKernel24
    rw [Finset.sum_sub_distrib]
    have hshift :
        (∑ j ∈ Finset.range m, t ^ (m + j)) =
          t ^ m * ∑ j ∈ Finset.range m, t ^ j := by
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro j _
      rw [pow_add]
    rw [hshift]
    ring
  rw [hkernel]
  apply (eq_div_iff (sub_ne_zero.mpr ht.symm)).2
  rw [mul_assoc, geom_sum_mul_neg]
  ring

/-- Rigorous polynomial-kernel version of
`r_m = ∫₀¹ (1 - t^m)² / (1 - t) dt`. -/
theorem integral_harmonicRemainderKernel24 (m : ℕ) :
    (∫ t : ℝ in 0..1, harmonicRemainderKernel24 m t) =
      harmonicRemainder24 m := by
  rw [harmonicRemainder24_eq_sum]
  unfold harmonicRemainderKernel24
  have hsum : ∀ n : ℕ,
      (∫ t : ℝ in 0..1, ∑ j ∈ Finset.range n,
          (t ^ j - t ^ (m + j))) =
        ∑ j ∈ Finset.range n,
          (1 / (j + 1 : ℝ) - 1 / (m + j + 1 : ℝ)) ∧
      IntervalIntegrable
        (fun t : ℝ ↦ ∑ j ∈ Finset.range n,
          (t ^ j - t ^ (m + j)))
        MeasureTheory.volume 0 1 := by
    intro n
    induction n with
    | zero =>
        simp
    | succ n ih =>
        have hleft :
            IntervalIntegrable (fun t : ℝ ↦ t ^ n)
              MeasureTheory.volume 0 1 :=
          (continuous_pow n).intervalIntegrable 0 1
        have hright :
            IntervalIntegrable (fun t : ℝ ↦ t ^ (m + n))
              MeasureTheory.volume 0 1 :=
          (continuous_pow (m + n)).intervalIntegrable 0 1
        have hterm :
            IntervalIntegrable (fun t : ℝ ↦ t ^ n - t ^ (m + n))
              MeasureTheory.volume 0 1 :=
          hleft.sub hright
        have hfun :
            (fun t : ℝ ↦ ∑ j ∈ Finset.range (n + 1),
              (t ^ j - t ^ (m + j))) =
              fun t : ℝ ↦
                (∑ j ∈ Finset.range n, (t ^ j - t ^ (m + j))) +
                  (t ^ n - t ^ (m + n)) := by
          funext t
          rw [Finset.sum_range_succ]
        constructor
        · rw [hfun, intervalIntegral.integral_add ih.2 hterm, ih.1,
            Finset.sum_range_succ,
            intervalIntegral.integral_sub hleft hright]
          simp only [integral_pow]
          norm_num
        · rw [hfun]
          exact ih.2.add hterm
  exact (hsum m).1

/-- The corrected beta-integral identity.  The draft used an extra factor
`1 / 2`; the exact denominator is `(2m + 1) * C(2m,m)`. -/
theorem betaIntegral_centralChoose24 (m : ℕ) :
    Complex.betaIntegral (m + 1) (m + 1) =
      1 / (((2 * m + 1 : ℕ) : ℂ) * (Nat.choose (2 * m) m : ℂ)) := by
  rw [Complex.betaIntegral_eq_Gamma_mul_div]
  · simp only [Nat.cast_add, Nat.cast_one]
    rw [Complex.Gamma_nat_eq_factorial m]
    rw [show (m : ℂ) + 1 + ((m : ℂ) + 1) =
        ((2 * m + 1 : ℕ) : ℂ) + 1 by
      push_cast
      ring]
    rw [Complex.Gamma_nat_eq_factorial (2 * m + 1)]
    have hchoose :=
      Nat.choose_mul_factorial_mul_factorial (n := 2 * m) (k := m)
        (by omega : m ≤ 2 * m)
    rw [show 2 * m - m = m by omega] at hchoose
    have hfact : Nat.factorial (2 * m + 1) =
        (2 * m + 1) * Nat.choose (2 * m) m *
          Nat.factorial m * Nat.factorial m := by
      rw [Nat.factorial_succ, ← hchoose]
      ring
    rw [hfact]
    push_cast
    field_simp [Nat.factorial_ne_zero]
  · norm_num [Complex.add_re]
    positivity
  · norm_num [Complex.add_re]
    positivity

theorem innerSum24_zero : innerSum24 0 = 0 := by
  norm_num [innerSum24, harmonicNumber]

theorem innerSum24_one : innerSum24 1 = 1 := by
  norm_num [innerSum24, harmonicNumber, Finset.sum_range_succ]

theorem innerSum24_two : innerSum24 2 = 25 / 4 := by
  norm_num [innerSum24, harmonicNumber, Finset.sum_range_succ]

theorem innerSum24_three : innerSum24 3 = 587 / 18 := by
  norm_num [innerSum24, harmonicNumber, Finset.sum_range_succ]

theorem sum_choose_sq_real (m : ℕ) :
    (∑ k ∈ Finset.range (m + 1), (Nat.choose m k : ℝ) ^ 2) =
      (Nat.choose (2 * m) m : ℝ) := by
  exact_mod_cast Nat.sum_range_choose_sq m

theorem sum_choose_mul_choose_succ (m : ℕ) :
    (∑ k ∈ Finset.range (m + 1),
      Nat.choose m k * Nat.choose (m + 1) (k + 1)) =
        Nat.choose (2 * m + 1) m := by
  rw [show 2 * m + 1 = m + (m + 1) by omega, Nat.add_choose_eq,
    Finset.Nat.sum_antidiagonal_eq_sum_range_succ_mk]
  apply Finset.sum_congr rfl
  intro k hk
  rw [Nat.choose_symm_of_eq_add (by
      have hkm : k ≤ m := Finset.mem_range_succ_iff.mp hk
      omega : m + 1 = (m - k) + (k + 1))]

theorem choose_div_succ (m k : ℕ) :
    (Nat.choose m k : ℝ) / (k + 1 : ℝ) =
      (Nat.choose (m + 1) (k + 1) : ℝ) / (m + 1 : ℝ) := by
  field_simp
  exact_mod_cast (by
    simpa [mul_comm] using Nat.add_one_mul_choose_eq m k)

theorem sum_choose_sq_div_succ (m : ℕ) :
    (∑ k ∈ Finset.range (m + 1),
      (Nat.choose m k : ℝ) ^ 2 / (k + 1 : ℝ)) =
        (Nat.choose (2 * m + 1) m : ℝ) / (m + 1 : ℝ) := by
  calc
    (∑ k ∈ Finset.range (m + 1),
        (Nat.choose m k : ℝ) ^ 2 / (k + 1 : ℝ))
        = (1 / (m + 1 : ℝ)) *
            ∑ k ∈ Finset.range (m + 1),
              (Nat.choose m k : ℝ) *
                (Nat.choose (m + 1) (k + 1) : ℝ) := by
          rw [Finset.mul_sum]
          apply Finset.sum_congr rfl
          intro k _
          rw [show (Nat.choose m k : ℝ) ^ 2 / (k + 1 : ℝ) =
              (Nat.choose m k : ℝ) *
                ((Nat.choose m k : ℝ) / (k + 1 : ℝ)) by ring,
            choose_div_succ]
          ring
    _ = (Nat.choose (2 * m + 1) m : ℝ) / (m + 1 : ℝ) := by
      rw [show (∑ k ∈ Finset.range (m + 1),
          (Nat.choose m k : ℝ) *
            (Nat.choose (m + 1) (k + 1) : ℝ)) =
          (Nat.choose (2 * m + 1) m : ℝ) by
        exact_mod_cast sum_choose_mul_choose_succ m]
      ring

theorem sum_choose_sq_shift (m : ℕ) :
    (∑ k ∈ Finset.range (m + 1),
      (Nat.choose m (k + 1) : ℝ) ^ 2) =
        (Nat.choose (2 * m) m : ℝ) - 1 := by
  calc
    (∑ k ∈ Finset.range (m + 1),
        (Nat.choose m (k + 1) : ℝ) ^ 2) =
        ∑ k ∈ Finset.range m,
          (Nat.choose m (k + 1) : ℝ) ^ 2 := by
            rw [Finset.sum_range_succ]
            simp
    _ = (Nat.choose (2 * m) m : ℝ) - 1 := by
      have h := sum_choose_sq_real m
      rw [Finset.sum_range_succ'] at h
      norm_num at h ⊢
      linarith

private def momentTelescoper24 (m : ℕ) : ℕ → ℝ
  | 0 => 0
  | k + 1 =>
      (Nat.choose m k : ℝ) ^ 2 *
        (2 * (k : ℝ) - 3 * (m : ℝ) - 1)

private def weightedMomentTelescoper24 (m k : ℕ) : ℝ :=
  (k : ℝ) * (Nat.choose m k : ℝ) ^ 2

private theorem choose_sq_difference_telescopes (m k : ℕ) (hk : k ≤ m) :
    (m + 1 : ℝ) * (Nat.choose (m + 1) k : ℝ) ^ 2 -
        2 * (2 * (m : ℝ) + 1) * (Nat.choose m k : ℝ) ^ 2 =
      momentTelescoper24 m (k + 1) - momentTelescoper24 m k := by
  cases k with
  | zero =>
      norm_num [momentTelescoper24]
      ring
  | succ k =>
      have hkm : k ≤ m := Nat.le_trans (Nat.le_succ k) hk
      have hchoose :
          (Nat.choose m (k + 1) : ℝ) * (k + 1 : ℝ) =
            (Nat.choose m k : ℝ) * ((m : ℝ) - k) := by
        rw [← Nat.cast_sub hkm]
        exact_mod_cast Nat.choose_succ_right_eq m k
      rw [Nat.choose_succ_succ]
      simp only [momentTelescoper24]
      push_cast
      linear_combination
        (-2 * ((Nat.choose m (k + 1) : ℝ) - (Nat.choose m k : ℝ))) * hchoose

private theorem weighted_moment_difference_telescopes
    (m k : ℕ) (hk : k ≤ m) :
    (((m + 1 : ℝ) ^ 2 / (k + 1 : ℝ) -
        (2 * (m : ℝ) + 1)) * (Nat.choose m k : ℝ) ^ 2) =
      weightedMomentTelescoper24 m (k + 1) -
        weightedMomentTelescoper24 m k := by
  have hchoose :
      (Nat.choose m (k + 1) : ℝ) * (k + 1 : ℝ) =
        (Nat.choose m k : ℝ) * ((m : ℝ) - k) := by
    rw [← Nat.cast_sub hk]
    exact_mod_cast Nat.choose_succ_right_eq m k
  have hsq := congrArg (fun x : ℝ ↦ x ^ 2) hchoose
  unfold weightedMomentTelescoper24
  field_simp
  push_cast
  linear_combination -hsq

private theorem moment_telescoper_sum (m : ℕ) :
    (∑ k ∈ Finset.range (m + 1),
        (momentTelescoper24 m (k + 1) - momentTelescoper24 m k) *
          harmonicNumber k) +
        (m + 1 : ℝ) * harmonicNumber (m + 1) =
      3 * (Nat.choose (2 * m + 1) m : ℝ) -
        2 * (Nat.choose (2 * m) m : ℝ) := by
  have hqend :
      momentTelescoper24 m (m + 1) = -(m + 1 : ℝ) := by
    simp [momentTelescoper24]
    ring
  calc
    (∑ k ∈ Finset.range (m + 1),
        (momentTelescoper24 m (k + 1) - momentTelescoper24 m k) *
          harmonicNumber k) +
        (m + 1 : ℝ) * harmonicNumber (m + 1)
        = (∑ k ∈ Finset.range (m + 1),
            ((momentTelescoper24 m (k + 1) * harmonicNumber (k + 1) -
              momentTelescoper24 m k * harmonicNumber k) -
            momentTelescoper24 m (k + 1) / (k + 1 : ℝ))) +
          (m + 1 : ℝ) * harmonicNumber (m + 1) := by
            congr 1
            apply Finset.sum_congr rfl
            intro k _
            rw [harmonicNumber_succ]
            ring
    _ = -∑ k ∈ Finset.range (m + 1),
          momentTelescoper24 m (k + 1) / (k + 1 : ℝ) := by
      rw [Finset.sum_sub_distrib,
        Finset.sum_range_sub
          (fun k ↦ momentTelescoper24 m k * harmonicNumber k) (m + 1),
        hqend]
      simp [momentTelescoper24, harmonicNumber_zero]
      ring
    _ = ∑ k ∈ Finset.range (m + 1),
          (3 * (m + 1 : ℝ) *
              ((Nat.choose m k : ℝ) ^ 2 / (k + 1 : ℝ)) -
            2 * (Nat.choose m k : ℝ) ^ 2) := by
      rw [← Finset.sum_neg_distrib]
      apply Finset.sum_congr rfl
      intro k _
      simp only [momentTelescoper24]
      field_simp
      ring
    _ = 3 * (m + 1 : ℝ) *
          (∑ k ∈ Finset.range (m + 1),
            (Nat.choose m k : ℝ) ^ 2 / (k + 1 : ℝ)) -
        2 * (∑ k ∈ Finset.range (m + 1),
          (Nat.choose m k : ℝ) ^ 2) := by
      rw [Finset.sum_sub_distrib, Finset.mul_sum, Finset.mul_sum]
    _ = 3 * (Nat.choose (2 * m + 1) m : ℝ) -
        2 * (Nat.choose (2 * m) m : ℝ) := by
      rw [sum_choose_sq_div_succ, sum_choose_sq_real]
      field_simp

theorem weightedFirstMoment24_identity (m : ℕ) :
    (m + 1 : ℝ) ^ 2 *
        (∑ k ∈ Finset.range (m + 1),
          (Nat.choose m k : ℝ) ^ 2 * harmonicNumber k / (k + 1 : ℝ)) -
      (2 * (m : ℝ) + 1) * firstMoment24 m =
        1 - (Nat.choose (2 * m) m : ℝ) := by
  calc
    (m + 1 : ℝ) ^ 2 *
          (∑ k ∈ Finset.range (m + 1),
            (Nat.choose m k : ℝ) ^ 2 * harmonicNumber k / (k + 1 : ℝ)) -
        (2 * (m : ℝ) + 1) * firstMoment24 m
        = ∑ k ∈ Finset.range (m + 1),
            ((((m + 1 : ℝ) ^ 2 / (k + 1 : ℝ) -
                (2 * (m : ℝ) + 1)) *
              (Nat.choose m k : ℝ) ^ 2) * harmonicNumber k) := by
          rw [firstMoment24, Finset.mul_sum, Finset.mul_sum,
            ← Finset.sum_sub_distrib]
          apply Finset.sum_congr rfl
          intro k _
          field_simp
    _ = ∑ k ∈ Finset.range (m + 1),
          (weightedMomentTelescoper24 m (k + 1) -
            weightedMomentTelescoper24 m k) * harmonicNumber k := by
      apply Finset.sum_congr rfl
      intro k hk
      rw [weighted_moment_difference_telescopes m k
        (Finset.mem_range_succ_iff.mp hk)]
    _ = -∑ k ∈ Finset.range (m + 1),
          weightedMomentTelescoper24 m (k + 1) / (k + 1 : ℝ) := by
      calc
        (∑ k ∈ Finset.range (m + 1),
            (weightedMomentTelescoper24 m (k + 1) -
              weightedMomentTelescoper24 m k) * harmonicNumber k)
            = (∑ k ∈ Finset.range (m + 1),
                ((weightedMomentTelescoper24 m (k + 1) *
                    harmonicNumber (k + 1) -
                  weightedMomentTelescoper24 m k * harmonicNumber k) -
                weightedMomentTelescoper24 m (k + 1) / (k + 1 : ℝ))) := by
                  apply Finset.sum_congr rfl
                  intro k _
                  rw [harmonicNumber_succ]
                  ring
        _ = -∑ k ∈ Finset.range (m + 1),
              weightedMomentTelescoper24 m (k + 1) / (k + 1 : ℝ) := by
          rw [Finset.sum_sub_distrib,
            Finset.sum_range_sub
              (fun k ↦ weightedMomentTelescoper24 m k * harmonicNumber k)
              (m + 1)]
          simp [weightedMomentTelescoper24]
    _ = -∑ k ∈ Finset.range (m + 1),
          (Nat.choose m (k + 1) : ℝ) ^ 2 := by
      congr 1
      apply Finset.sum_congr rfl
      intro k _
      simp [weightedMomentTelescoper24]
      field_simp
    _ = 1 - (Nat.choose (2 * m) m : ℝ) := by
      rw [sum_choose_sq_shift]
      ring

theorem weightedFirstMoment24_closed (m : ℕ) :
    (∑ k ∈ Finset.range (m + 1),
      (Nat.choose m k : ℝ) ^ 2 * harmonicNumber k / (k + 1 : ℝ)) =
        ((2 * (m : ℝ) + 1) * firstMoment24 m -
          (Nat.choose (2 * m) m : ℝ) + 1) / (m + 1 : ℝ) ^ 2 := by
  have h := weightedFirstMoment24_identity m
  rw [eq_div_iff (by positivity : (m + 1 : ℝ) ^ 2 ≠ 0)]
  linarith

theorem sum_choose_sq_div_succ_sq (m : ℕ) :
    (∑ k ∈ Finset.range (m + 1),
      (Nat.choose m k : ℝ) ^ 2 / (k + 1 : ℝ) ^ 2) =
        ((Nat.choose (2 * (m + 1)) (m + 1) : ℝ) - 1) /
          (m + 1 : ℝ) ^ 2 := by
  have hshift := sum_choose_sq_shift (m + 1)
  rw [Finset.sum_range_succ] at hshift
  simp at hshift
  calc
    (∑ k ∈ Finset.range (m + 1),
        (Nat.choose m k : ℝ) ^ 2 / (k + 1 : ℝ) ^ 2)
        = (1 / (m + 1 : ℝ) ^ 2) *
            ∑ k ∈ Finset.range (m + 1),
              (Nat.choose (m + 1) (k + 1) : ℝ) ^ 2 := by
          rw [Finset.mul_sum]
          apply Finset.sum_congr rfl
          intro k _
          have h := choose_div_succ m k
          calc
            (Nat.choose m k : ℝ) ^ 2 / (k + 1 : ℝ) ^ 2 =
                ((Nat.choose m k : ℝ) / (k + 1 : ℝ)) ^ 2 := by
                  rw [div_pow]
            _ = ((Nat.choose (m + 1) (k + 1) : ℝ) /
                (m + 1 : ℝ)) ^ 2 := by rw [h]
            _ = (1 / (m + 1 : ℝ) ^ 2) *
                (Nat.choose (m + 1) (k + 1) : ℝ) ^ 2 := by
                  field_simp
    _ = ((Nat.choose (2 * (m + 1)) (m + 1) : ℝ) - 1) /
          (m + 1 : ℝ) ^ 2 := by
      rw [hshift]
      ring

theorem harmonicNumber_sq_succ (k : ℕ) :
    harmonicNumber (k + 1) ^ 2 =
      harmonicNumber k ^ 2 +
        2 * harmonicNumber k / (k + 1 : ℝ) +
        1 / (k + 1 : ℝ) ^ 2 := by
  rw [harmonicNumber_succ]
  field_simp
  ring

private theorem second_moment_telescoper_sum (m : ℕ) :
    (∑ k ∈ Finset.range (m + 1),
        (momentTelescoper24 m (k + 1) - momentTelescoper24 m k) *
          harmonicNumber k ^ 2) +
        (m + 1 : ℝ) * harmonicNumber (m + 1) ^ 2 =
      6 * (m + 1 : ℝ) *
          (∑ k ∈ Finset.range (m + 1),
            (Nat.choose m k : ℝ) ^ 2 * harmonicNumber k / (k + 1 : ℝ)) -
        4 * firstMoment24 m +
        3 * (m + 1 : ℝ) *
          (∑ k ∈ Finset.range (m + 1),
            (Nat.choose m k : ℝ) ^ 2 / (k + 1 : ℝ) ^ 2) -
        2 * (∑ k ∈ Finset.range (m + 1),
          (Nat.choose m k : ℝ) ^ 2 / (k + 1 : ℝ)) := by
  have hqend :
      momentTelescoper24 m (m + 1) = -(m + 1 : ℝ) := by
    simp [momentTelescoper24]
    ring
  calc
    (∑ k ∈ Finset.range (m + 1),
        (momentTelescoper24 m (k + 1) - momentTelescoper24 m k) *
          harmonicNumber k ^ 2) +
        (m + 1 : ℝ) * harmonicNumber (m + 1) ^ 2 =
      (∑ k ∈ Finset.range (m + 1),
          ((momentTelescoper24 m (k + 1) * harmonicNumber (k + 1) ^ 2 -
              momentTelescoper24 m k * harmonicNumber k ^ 2) -
            momentTelescoper24 m (k + 1) *
              (2 * harmonicNumber k / (k + 1 : ℝ) +
                1 / (k + 1 : ℝ) ^ 2))) +
        (m + 1 : ℝ) * harmonicNumber (m + 1) ^ 2 := by
          congr 1
          apply Finset.sum_congr rfl
          intro k _
          rw [harmonicNumber_sq_succ]
          ring
    _ = -∑ k ∈ Finset.range (m + 1),
          momentTelescoper24 m (k + 1) *
            (2 * harmonicNumber k / (k + 1 : ℝ) +
              1 / (k + 1 : ℝ) ^ 2) := by
      rw [Finset.sum_sub_distrib,
        Finset.sum_range_sub
          (fun k ↦ momentTelescoper24 m k * harmonicNumber k ^ 2) (m + 1),
        hqend]
      simp [momentTelescoper24, harmonicNumber_zero]
      ring
    _ = ∑ k ∈ Finset.range (m + 1),
          (6 * (m + 1 : ℝ) *
              ((Nat.choose m k : ℝ) ^ 2 * harmonicNumber k /
                (k + 1 : ℝ)) -
            4 * ((Nat.choose m k : ℝ) ^ 2 * harmonicNumber k) +
            3 * (m + 1 : ℝ) *
              ((Nat.choose m k : ℝ) ^ 2 / (k + 1 : ℝ) ^ 2) -
            2 * ((Nat.choose m k : ℝ) ^ 2 / (k + 1 : ℝ))) := by
      rw [← Finset.sum_neg_distrib]
      apply Finset.sum_congr rfl
      intro k _
      simp only [momentTelescoper24]
      field_simp
      ring
    _ = 6 * (m + 1 : ℝ) *
          (∑ k ∈ Finset.range (m + 1),
            (Nat.choose m k : ℝ) ^ 2 * harmonicNumber k / (k + 1 : ℝ)) -
        4 * firstMoment24 m +
        3 * (m + 1 : ℝ) *
          (∑ k ∈ Finset.range (m + 1),
            (Nat.choose m k : ℝ) ^ 2 / (k + 1 : ℝ) ^ 2) -
        2 * (∑ k ∈ Finset.range (m + 1),
          (Nat.choose m k : ℝ) ^ 2 / (k + 1 : ℝ)) := by
      rw [firstMoment24, Finset.sum_sub_distrib, Finset.sum_add_distrib,
        Finset.sum_sub_distrib, Finset.mul_sum, Finset.mul_sum,
        Finset.mul_sum, Finset.mul_sum]

theorem innerSum24_succ (m : ℕ) :
    innerSum24 (m + 1) =
      (∑ k ∈ Finset.range (m + 1),
        (Nat.choose (m + 1) k : ℝ) ^ 2 * harmonicNumber k ^ 2) +
        harmonicNumber (m + 1) ^ 2 := by
  simp [innerSum24, Finset.sum_range_succ]

theorem innerSum24_recurrence_raw (m : ℕ) :
    (m + 1 : ℝ) * innerSum24 (m + 1) -
        2 * (2 * (m : ℝ) + 1) * innerSum24 m =
      6 * (m + 1 : ℝ) *
          (∑ k ∈ Finset.range (m + 1),
            (Nat.choose m k : ℝ) ^ 2 * harmonicNumber k / (k + 1 : ℝ)) -
        4 * firstMoment24 m +
        3 * (m + 1 : ℝ) *
          (∑ k ∈ Finset.range (m + 1),
            (Nat.choose m k : ℝ) ^ 2 / (k + 1 : ℝ) ^ 2) -
        2 * (∑ k ∈ Finset.range (m + 1),
          (Nat.choose m k : ℝ) ^ 2 / (k + 1 : ℝ)) := by
  rw [innerSum24_succ, innerSum24]
  calc
    (m + 1 : ℝ) *
          ((∑ k ∈ Finset.range (m + 1),
            (Nat.choose (m + 1) k : ℝ) ^ 2 * harmonicNumber k ^ 2) +
            harmonicNumber (m + 1) ^ 2) -
        2 * (2 * (m : ℝ) + 1) *
          (∑ k ∈ Finset.range (m + 1),
            (Nat.choose m k : ℝ) ^ 2 * harmonicNumber k ^ 2) =
      (∑ k ∈ Finset.range (m + 1),
          (((m + 1 : ℝ) * (Nat.choose (m + 1) k : ℝ) ^ 2 -
            2 * (2 * (m : ℝ) + 1) * (Nat.choose m k : ℝ) ^ 2) *
            harmonicNumber k ^ 2)) +
        (m + 1 : ℝ) * harmonicNumber (m + 1) ^ 2 := by
          rw [mul_add, Finset.mul_sum, Finset.mul_sum,
            ← sub_add_eq_add_sub, ← Finset.sum_sub_distrib]
          congr 1
          apply Finset.sum_congr rfl
          intro k _
          ring
    _ = (∑ k ∈ Finset.range (m + 1),
          (momentTelescoper24 m (k + 1) - momentTelescoper24 m k) *
            harmonicNumber k ^ 2) +
        (m + 1 : ℝ) * harmonicNumber (m + 1) ^ 2 := by
      congr 1
      apply Finset.sum_congr rfl
      intro k hk
      rw [choose_sq_difference_telescopes m k
        (Finset.mem_range_succ_iff.mp hk)]
    _ = _ := second_moment_telescoper_sum m

theorem firstMoment24_succ (m : ℕ) :
    firstMoment24 (m + 1) =
      (∑ k ∈ Finset.range (m + 1),
        (Nat.choose (m + 1) k : ℝ) ^ 2 * harmonicNumber k) +
        harmonicNumber (m + 1) := by
  simp [firstMoment24, Finset.sum_range_succ]

theorem firstMoment24_recurrence (m : ℕ) :
    (m + 1 : ℝ) * firstMoment24 (m + 1) -
        2 * (2 * (m : ℝ) + 1) * firstMoment24 m =
      3 * (Nat.choose (2 * m + 1) m : ℝ) -
        2 * (Nat.choose (2 * m) m : ℝ) := by
  rw [firstMoment24_succ, firstMoment24]
  calc
    (m + 1 : ℝ) *
          ((∑ k ∈ Finset.range (m + 1),
            (Nat.choose (m + 1) k : ℝ) ^ 2 * harmonicNumber k) +
            harmonicNumber (m + 1)) -
        2 * (2 * (m : ℝ) + 1) *
          (∑ k ∈ Finset.range (m + 1),
            (Nat.choose m k : ℝ) ^ 2 * harmonicNumber k)
        = (∑ k ∈ Finset.range (m + 1),
            (((m + 1 : ℝ) * (Nat.choose (m + 1) k : ℝ) ^ 2 -
              2 * (2 * (m : ℝ) + 1) * (Nat.choose m k : ℝ) ^ 2) *
              harmonicNumber k)) +
            (m + 1 : ℝ) * harmonicNumber (m + 1) := by
          rw [mul_add, Finset.mul_sum, Finset.mul_sum,
            ← sub_add_eq_add_sub,
            ← Finset.sum_sub_distrib]
          congr 1
          apply Finset.sum_congr rfl
          intro k _
          ring
    _ = (∑ k ∈ Finset.range (m + 1),
          (momentTelescoper24 m (k + 1) - momentTelescoper24 m k) *
            harmonicNumber k) +
          (m + 1 : ℝ) * harmonicNumber (m + 1) := by
      congr 1
      apply Finset.sum_congr rfl
      intro k hk
      rw [choose_sq_difference_telescopes m k
        (Finset.mem_range_succ_iff.mp hk)]
    _ = 3 * (Nat.choose (2 * m + 1) m : ℝ) -
        2 * (Nat.choose (2 * m) m : ℝ) :=
      moment_telescoper_sum m

theorem signedHarmonic24_succ (n : ℕ) :
    signedHarmonic24 (n + 1) =
      signedHarmonic24 n + (-1 : ℝ) ^ (n + 1) / (n + 1 : ℝ) := by
  simp [signedHarmonic24, Finset.sum_range_succ]

theorem parityRemainder24_succ (n : ℕ) :
    parityRemainder24 (n + 1) =
      parityRemainder24 n + 1 / (n + 1 : ℝ) +
        2 * ((-1 : ℝ) ^ (n + 1) / (n + 1 : ℝ)) := by
  rw [parityRemainder24, parityRemainder24, harmonicNumber_succ,
    signedHarmonic24_succ]
  ring

theorem abs_signedHarmonic24_le_harmonicNumber (n : ℕ) :
    |signedHarmonic24 n| ≤ harmonicNumber n := by
  unfold signedHarmonic24 harmonicNumber
  calc
    |∑ k ∈ Finset.range n,
        (-1 : ℝ) ^ (k + 1) / (k + 1 : ℝ)| ≤
        ∑ k ∈ Finset.range n,
          |(-1 : ℝ) ^ (k + 1) / (k + 1 : ℝ)| :=
      Finset.abs_sum_le_sum_abs _ _
    _ = ∑ k ∈ Finset.range n, 1 / (k + 1 : ℝ) := by
      apply Finset.sum_congr rfl
      intro k _
      rw [abs_div, abs_pow, abs_neg, abs_one, one_pow,
        abs_of_pos (by positivity : (0 : ℝ) < k + 1)]

theorem abs_parityRemainder24_le (n : ℕ) :
    |parityRemainder24 n| ≤ 3 * harmonicNumber n := by
  unfold parityRemainder24
  calc
    |harmonicNumber n + 2 * signedHarmonic24 n| ≤
        |harmonicNumber n| + |2 * signedHarmonic24 n| :=
      abs_add_le _ _
    _ = harmonicNumber n + 2 * |signedHarmonic24 n| := by
      rw [abs_of_nonneg (harmonicNumber_nonneg n), abs_mul]
      norm_num
    _ ≤ harmonicNumber n + 2 * harmonicNumber n := by
      gcongr
      exact abs_signedHarmonic24_le_harmonicNumber n
    _ = 3 * harmonicNumber n := by ring

theorem harmonicRemainder24_succ (m : ℕ) :
    harmonicRemainder24 (m + 1) =
      harmonicRemainder24 m +
        2 / (m + 1 : ℝ) -
        1 / (2 * (m : ℝ) + 1) -
        1 / (2 * (m : ℝ) + 2) := by
  have heven :
      harmonicNumber (2 * (m + 1)) =
        harmonicNumber (2 * m) +
          1 / (2 * (m : ℝ) + 1) +
          1 / (2 * (m : ℝ) + 2) := by
    rw [show 2 * (m + 1) = (2 * m + 1) + 1 by omega,
      harmonicNumber_succ, harmonicNumber_succ]
    push_cast
    ring
  rw [harmonicRemainder24, harmonicRemainder24, harmonicNumber_succ, heven]
  ring

theorem parityRemainder24_even (m : ℕ) :
    parityRemainder24 (2 * m) = harmonicRemainder24 m := by
  induction m with
  | zero =>
      simp [parityRemainder24, signedHarmonic24, harmonicRemainder24,
        harmonicNumber]
  | succ m ih =>
      rw [show 2 * (m + 1) = (2 * m + 1) + 1 by omega,
        parityRemainder24_succ, parityRemainder24_succ, ih,
        harmonicRemainder24_succ]
      have hodd : (-1 : ℝ) ^ (2 * m + 1) = -1 := by
        rw [pow_add, pow_mul]
        norm_num
      have heven : (-1 : ℝ) ^ (2 * m + 1 + 1) = 1 := by
        rw [show 2 * m + 1 + 1 = 2 * (m + 1) by omega, pow_mul]
        norm_num
      rw [hodd, heven]
      push_cast
      field_simp
      ring

theorem centralChoose24_succ (m : ℕ) :
    (Nat.choose (2 * (m + 1)) (m + 1) : ℝ) =
      2 * (Nat.choose (2 * m + 1) m : ℝ) := by
  exact_mod_cast (by
    rw [show 2 * (m + 1) = (2 * m + 1) + 1 by omega,
      Nat.choose_succ_succ', Nat.choose_symm_half]
    omega)

theorem centralChoose24_odd_relation (m : ℕ) :
    (m + 1 : ℝ) * (Nat.choose (2 * m + 1) m : ℝ) =
      (2 * (m : ℝ) + 1) * (Nat.choose (2 * m) m : ℝ) := by
  exact_mod_cast (by
    have hsub : 2 * m + 1 - m = m + 1 := by omega
    calc
      (m + 1) * Nat.choose (2 * m + 1) m =
          Nat.choose (2 * m + 1) m * (m + 1) := by
            exact Nat.mul_comm _ _
      _ = Nat.choose (2 * m) m * (2 * m + 1) := by
        simpa [hsub] using (Nat.choose_mul_succ_eq (2 * m) m).symm
      _ = (2 * m + 1) * Nat.choose (2 * m) m := by
        exact Nat.mul_comm _ _)

theorem centralChoose24_pos (m : ℕ) :
    0 < (Nat.choose (2 * m) m : ℝ) := by
  exact_mod_cast Nat.choose_pos (by omega : m ≤ 2 * m)

theorem harmonicSquare24_succ (n : ℕ) :
    harmonicSquare24 (n + 1) =
      harmonicSquare24 n + 1 / (n + 1 : ℝ) ^ 2 := by
  simp [harmonicSquare24, Finset.sum_range_succ]

theorem harmonicSquare24_le_harmonicNumber (n : ℕ) :
    harmonicSquare24 n ≤ harmonicNumber n := by
  unfold harmonicSquare24 harmonicNumber
  apply Finset.sum_le_sum
  intro j hj
  have hx : (1 : ℝ) ≤ j + 1 := by
    exact_mod_cast Nat.succ_le_succ (Nat.zero_le j)
  apply one_div_le_one_div_of_le (by positivity)
  nlinarith [mul_nonneg (sub_nonneg.mpr hx)
    (show (0 : ℝ) ≤ j + 1 by positivity)]

theorem harmonicSquare24_double_succ (m : ℕ) :
    harmonicSquare24 (2 * (m + 1)) =
      harmonicSquare24 (2 * m) +
        1 / (2 * (m : ℝ) + 1) ^ 2 +
        1 / (2 * (m : ℝ) + 2) ^ 2 := by
  rw [show 2 * (m + 1) = (2 * m + 1) + 1 by omega,
    harmonicSquare24_succ, harmonicSquare24_succ]
  push_cast
  ring

/-- Extract the positive even denominators from a series and its
outer-alternating companion.  The zero-based position `2m + 1`
corresponds to the positive integer `2m + 2`. -/
theorem hasSum_even_position24 {f : ℕ → ℝ} {a b : ℝ}
    (hplus : HasSum f a)
    (hminus : HasSum (fun n : ℕ ↦ (-1 : ℝ) ^ (n + 1) * f n) b) :
    HasSum (fun m : ℕ ↦ f (2 * m + 1)) ((a + b) / 2) := by
  let g : ℕ → ℝ := fun n ↦
    (1 / 2 : ℝ) * (f n + (-1 : ℝ) ^ (n + 1) * f n)
  have hg : HasSum g ((a + b) / 2) := by
    change HasSum
      (fun n ↦ (1 / 2 : ℝ) *
        (f n + (-1 : ℝ) ^ (n + 1) * f n)) ((a + b) / 2)
    convert (hplus.add hminus).mul_left (1 / 2 : ℝ) using 1
    all_goals ring
  have hgs : Summable g := hg.summable
  have he : Summable (fun m : ℕ ↦ g (2 * m)) :=
    hgs.comp_injective
      (mul_right_injective₀ (by norm_num : (2 : ℕ) ≠ 0))
  have ho : Summable (fun m : ℕ ↦ g (2 * m + 1)) :=
    hgs.comp_injective ((add_left_injective 1).comp
      (mul_right_injective₀ (by norm_num : (2 : ℕ) ≠ 0)))
  have hsplit := tsum_even_add_odd he ho
  have hezero : (∑' m : ℕ, g (2 * m)) = 0 := by
    calc
      (∑' m : ℕ, g (2 * m)) = ∑' _m : ℕ, (0 : ℝ) := by
        apply tsum_congr
        intro m
        simp [g, pow_add, pow_mul]
      _ = 0 := tsum_zero
  have hodd : (∑' m : ℕ, g (2 * m + 1)) = (a + b) / 2 := by
    rw [hezero, zero_add, hg.tsum_eq] at hsplit
    exact hsplit
  have hhas := ho.hasSum
  rw [hodd] at hhas
  convert hhas using 1
  funext m
  simp [g, pow_add, pow_mul]
  ring

/-- Group a convergent real series into consecutive pairs. -/
theorem HasSum.pair_consecutive24 {f : ℕ → ℝ} {a : ℝ}
    (h : HasSum f a) :
    HasSum (fun m : ℕ ↦ f (2 * m) + f (2 * m + 1)) a := by
  have hs : Summable f := h.summable
  have he : Summable (fun m : ℕ ↦ f (2 * m)) :=
    hs.comp_injective
      (mul_right_injective₀ (by norm_num : (2 : ℕ) ≠ 0))
  have ho : Summable (fun m : ℕ ↦ f (2 * m + 1)) :=
    hs.comp_injective ((add_left_injective 1).comp
      (mul_right_injective₀ (by norm_num : (2 : ℕ) ≠ 0)))
  have hsplit := tsum_even_add_odd he ho
  have hpairs : Summable (fun m : ℕ ↦ f (2 * m) + f (2 * m + 1)) :=
    he.add ho
  have hpair_tsum :
      (∑' m : ℕ, (f (2 * m) + f (2 * m + 1))) = a := by
    calc
      (∑' m : ℕ, (f (2 * m) + f (2 * m + 1))) =
          (∑' m : ℕ, f (2 * m)) + ∑' m : ℕ, f (2 * m + 1) :=
        he.tsum_add ho
      _ = ∑' m : ℕ, f m := hsplit
      _ = a := h.tsum_eq
  have hh := hpairs.hasSum
  rw [hpair_tsum] at hh
  exact hh

/-- Pointwise even-index decomposition behind the Euler-sum certificate for
the elementary scalar series. -/
theorem elementaryOuterTerm24_euler_decomposition (m : ℕ) :
    elementaryOuterTerm24 m =
      4 * quadraticEulerTerm24 (2 * m + 1) +
        8 * shiftedLinearEulerTerm24 m -
        24 * cubicLinearEulerTerm24 (2 * m + 1) +
        rationalCorrectionTerm24 m := by
  unfold elementaryOuterTerm24 quadraticEulerTerm24
    shiftedLinearEulerTerm24 cubicLinearEulerTerm24
    rationalCorrectionTerm24
  dsimp only
  rw [show 2 * m + 1 + 1 = 2 * (m + 1) by omega,
    parityRemainder24_even, harmonicRemainder24_succ,
    harmonicSquare24_double_succ]
  push_cast
  field_simp
  ring_nf
  field_simp
  ring

/-- Partial-fraction decomposition of the shifted linear sum into one
alternating weight-two sum, the grouped harmonic series, and the even part
of two weight-three sums. -/
theorem shiftedLinearEulerTerm24_lower_weight_decomposition (m : ℕ) :
    shiftedLinearEulerTerm24 m =
      -(alternatingLinearEulerTerm24 (2 * m) +
          alternatingLinearEulerTerm24 (2 * m + 1)) +
        3 * (1 / (2 * m + 1 : ℝ) - 1 / (2 * m + 2 : ℝ)) -
        quadraticLinearEulerTerm24 (2 * m + 1) := by
  unfold shiftedLinearEulerTerm24 alternatingLinearEulerTerm24
    linearEulerTerm24 quadraticLinearEulerTerm24
  have hodd : (-1 : ℝ) ^ (2 * m + 1) = -1 := by
    rw [pow_add, pow_mul]
    norm_num
  have heven : (-1 : ℝ) ^ (2 * m + 1 + 1) = 1 := by
    rw [show 2 * m + 1 + 1 = 2 * (m + 1) by omega, pow_mul]
    norm_num
  have hrem :
      parityRemainder24 (2 * (m + 1)) =
        parityRemainder24 (2 * m + 1) +
          3 / (2 * m + 2 : ℝ) := by
    rw [show 2 * (m + 1) = (2 * m + 1) + 1 by omega,
      parityRemainder24_succ, heven]
    push_cast
    ring
  rw [hodd, heven,
    show 2 * m + 1 + 1 = 2 * (m + 1) by omega, hrem]
  push_cast
  field_simp
  ring

theorem pairedAlternatingLinearEulerTerm24_decomposition (m : ℕ) :
    pairedAlternatingLinearEulerTerm24 m =
      pairedAlternatingHarmonicEulerTerm24 m +
        2 * pairedAlternatingSignedHarmonicEulerTerm24 m := by
  unfold pairedAlternatingLinearEulerTerm24
    alternatingLinearEulerTerm24 linearEulerTerm24
    pairedAlternatingHarmonicEulerTerm24
    pairedAlternatingSignedHarmonicEulerTerm24 parityRemainder24
  have hodd : (-1 : ℝ) ^ (2 * m + 1) = -1 := by
    rw [pow_add, pow_mul]
    norm_num
  have heven : (-1 : ℝ) ^ (2 * m + 1 + 1) = 1 := by
    rw [show 2 * m + 1 + 1 = 2 * (m + 1) by omega, pow_mul]
    norm_num
  rw [hodd, heven]
  push_cast
  ring

theorem pairedAlternatingHarmonicEulerTerm24_formula (m : ℕ) :
    pairedAlternatingHarmonicEulerTerm24 m =
      -harmonicNumber (2 * m + 1) *
          (1 / (2 * m + 1 : ℝ) - 1 / (2 * m + 2 : ℝ)) +
        1 / (2 * m + 2 : ℝ) ^ 2 := by
  unfold pairedAlternatingHarmonicEulerTerm24
  have hsucc :
      harmonicNumber (2 * m + 2) =
        harmonicNumber (2 * m + 1) +
          1 / (2 * m + 2 : ℝ) := by
    rw [show 2 * m + 2 = (2 * m + 1) + 1 by omega,
      harmonicNumber_succ]
    push_cast
    ring
  rw [hsucc]
  field_simp
  ring

theorem shifted_zeta_two_hasSum :
    HasSum (fun n : ℕ ↦ 1 / (n + 1 : ℝ) ^ 2) (Real.pi ^ 2 / 6) := by
  have h := (hasSum_nat_add_iff' 1).2 hasSum_zeta_two
  simpa using h

theorem shifted_zeta_four_hasSum24 :
    HasSum (fun n : ℕ ↦ 1 / (n + 1 : ℝ) ^ 4) (Real.pi ^ 4 / 90) := by
  have h := (hasSum_nat_add_iff' 1).2 hasSum_zeta_four
  simpa using h

theorem odd_square_hasSum24 :
    HasSum (fun m : ℕ ↦ 1 / (2 * m + 1 : ℝ) ^ 2)
      ((3 / 4) * (Real.pi ^ 2 / 6)) := by
  have hpairs := shifted_zeta_two_hasSum.pair_consecutive24
  have heven :
      HasSum (fun m : ℕ ↦ 1 / (2 * m + 2 : ℝ) ^ 2)
        ((1 / 4) * (Real.pi ^ 2 / 6)) := by
    convert shifted_zeta_two_hasSum.mul_left (1 / 4) using 1
    · funext m
      field_simp
      ring
  convert hpairs.sub heven using 1
  · funext m
    push_cast
    ring
  · ring

/-- The alternating harmonic series, grouped into consecutive positive
pairs.  This version is convenient because every summand is nonnegative. -/
theorem groupedAlternatingHarmonic_hasSum24 :
    HasSum
      (fun m : ℕ ↦
        1 / (2 * m + 1 : ℝ) - 1 / (2 * m + 2 : ℝ))
      (Real.log 2) := by
  have hpartial (N : ℕ) :
      (∑ m ∈ Finset.range N,
          (1 / (2 * m + 1 : ℝ) - 1 / (2 * m + 2 : ℝ))) =
        harmonicNumber (2 * N) - harmonicNumber N := by
    induction N with
    | zero =>
        simp [harmonicNumber_zero]
    | succ N ih =>
        rw [Finset.sum_range_succ, ih,
          show 2 * (N + 1) = (2 * N + 1) + 1 by omega,
          harmonicNumber_succ, harmonicNumber_succ,
          harmonicNumber_succ]
        push_cast
        field_simp
        ring
  have htwo : Tendsto (fun n : ℕ ↦ 2 * n) atTop atTop := by
    exact (strictMono_mul_left_of_pos (by norm_num : 0 < (2 : ℕ))).tendsto_atTop
  have hzero :
      Tendsto
        (fun n : ℕ ↦
          ((harmonic (2 * n) : ℝ) - Real.log (2 * n)) -
            ((harmonic n : ℝ) - Real.log n))
        atTop (𝓝 0) := by
    simpa [Function.comp_def, mul_comm] using
      (Real.tendsto_harmonic_sub_log.comp htwo).sub
        Real.tendsto_harmonic_sub_log
  have hlimit :
      Tendsto
        (fun n : ℕ ↦ harmonicNumber (2 * n) - harmonicNumber n)
        atTop (𝓝 (Real.log 2)) := by
    have hbase :
        Tendsto
          (fun n : ℕ ↦
            (((harmonic (2 * n) : ℝ) - Real.log (2 * n)) -
                ((harmonic n : ℝ) - Real.log n)) +
              Real.log 2)
          atTop (𝓝 (Real.log 2)) := by
      simpa using hzero.add
        (tendsto_const_nhds :
          Tendsto (fun _ : ℕ ↦ Real.log 2) atTop (𝓝 (Real.log 2)))
    apply hbase.congr'
    filter_upwards [eventually_ne_atTop 0] with n hn
    rw [harmonicNumber_eq_cast_harmonic,
      harmonicNumber_eq_cast_harmonic]
    rw [Real.log_mul (by norm_num : (2 : ℝ) ≠ 0)
      (Nat.cast_ne_zero.mpr hn)]
    ring
  rw [hasSum_iff_tendsto_nat_of_nonneg]
  · simpa only [hpartial] using hlimit
  · intro m
    apply sub_nonneg.mpr
    apply one_div_le_one_div_of_le
    · positivity
    · norm_num

/-- Even partial sums of the raw alternating harmonic series are the
negatives of the positive grouped series above. -/
theorem signedHarmonic24_even_partial (N : ℕ) :
    signedHarmonic24 (2 * N) =
      -(∑ m ∈ Finset.range N,
        (1 / (2 * m + 1 : ℝ) - 1 / (2 * m + 2 : ℝ))) := by
  induction N with
  | zero =>
      simp [signedHarmonic24]
  | succ N ih =>
      rw [show 2 * (N + 1) = (2 * N + 1) + 1 by omega,
        signedHarmonic24_succ, signedHarmonic24_succ,
        Finset.sum_range_succ, ih]
      have hodd : (-1 : ℝ) ^ (2 * N + 1) = -1 := by
        rw [pow_add, pow_mul]
        norm_num
      have heven : (-1 : ℝ) ^ (2 * N + 1 + 1) = 1 := by
        rw [show 2 * N + 1 + 1 = 2 * (N + 1) by omega, pow_mul]
        norm_num
      rw [hodd, heven]
      push_cast
      ring

theorem tendsto_signedHarmonic24_even :
    Tendsto (fun N : ℕ ↦ signedHarmonic24 (2 * N)) atTop
      (𝓝 (-Real.log 2)) := by
  have h := groupedAlternatingHarmonic_hasSum24.tendsto_sum_nat.neg
  simpa only [signedHarmonic24_even_partial] using h

/-- Finite triangular-sum identity for the signed harmonic numbers. -/
theorem signedHarmonic24_triangular_partial (N : ℕ) :
    (∑ n ∈ Finset.range N,
        (-1 : ℝ) ^ (n + 1) * signedHarmonic24 (n + 1) /
          (n + 1 : ℝ)) =
      (1 / 2) *
        (signedHarmonic24 N ^ 2 + harmonicSquare24 N) := by
  induction N with
  | zero =>
      simp [signedHarmonic24, harmonicSquare24]
  | succ N ih =>
      rw [Finset.sum_range_succ, ih, signedHarmonic24_succ,
        harmonicSquare24_succ]
      have hsignsq : ((-1 : ℝ) ^ (N + 1)) ^ 2 = 1 := by
        rw [← pow_mul]
        norm_num
      have hden : (N + 1 : ℝ) ≠ 0 := by positivity
      field_simp [hden]
      nlinarith [hsignsq]

/-- Pairing the raw finite triangular series gives the summand used in the
absolutely convergent signed-harmonic component. -/
theorem pairedAlternatingSignedHarmonic_partial (N : ℕ) :
    (∑ m ∈ Finset.range N,
        pairedAlternatingSignedHarmonicEulerTerm24 m) =
      ∑ n ∈ Finset.range (2 * N),
        (-1 : ℝ) ^ (n + 1) * signedHarmonic24 (n + 1) /
          (n + 1 : ℝ) := by
  induction N with
  | zero =>
      simp
  | succ N ih =>
      rw [Finset.sum_range_succ, ih,
        show 2 * (N + 1) = (2 * N + 1) + 1 by omega,
        Finset.sum_range_succ, Finset.sum_range_succ]
      unfold pairedAlternatingSignedHarmonicEulerTerm24
      have hodd : (-1 : ℝ) ^ (2 * N + 1) = -1 := by
        rw [pow_add, pow_mul]
        norm_num
      have heven : (-1 : ℝ) ^ (2 * N + 1 + 1) = 1 := by
        rw [show 2 * N + 1 + 1 = 2 * (N + 1) by omega, pow_mul]
        norm_num
      rw [hodd, heven]
      push_cast
      ring

theorem pairedAlternatingSignedHarmonicEulerTerm24_nonneg (m : ℕ) :
    0 ≤ pairedAlternatingSignedHarmonicEulerTerm24 m := by
  let G : ℝ :=
    ∑ i ∈ Finset.range m,
      (1 / (2 * i + 1 : ℝ) - 1 / (2 * i + 2 : ℝ))
  have hgroup (i : ℕ) :
      0 ≤ 1 / (2 * i + 1 : ℝ) - 1 / (2 * i + 2 : ℝ) := by
    apply sub_nonneg.mpr
    apply one_div_le_one_div_of_le
    · positivity
    · norm_num
  have hG : 0 ≤ G := by
    unfold G
    exact Finset.sum_nonneg fun i _ ↦ hgroup i
  have hform :
      pairedAlternatingSignedHarmonicEulerTerm24 m =
        G * (1 / (2 * m + 1 : ℝ) - 1 / (2 * m + 2 : ℝ)) +
          (1 / (2 * m + 1 : ℝ)) *
            (1 / (2 * m + 1 : ℝ) - 1 / (2 * m + 2 : ℝ)) +
          1 / (2 * m + 2 : ℝ) ^ 2 := by
    unfold pairedAlternatingSignedHarmonicEulerTerm24
    rw [signedHarmonic24_succ (2 * m),
      signedHarmonic24_even_partial m,
      show 2 * m + 2 = 2 * (m + 1) by omega,
      signedHarmonic24_even_partial (m + 1),
      Finset.sum_range_succ]
    have hodd : (-1 : ℝ) ^ (2 * m + 1) = -1 := by
      rw [pow_add, pow_mul]
      norm_num
    rw [hodd]
    push_cast
    simp only [G]
    field_simp
    ring
  rw [hform]
  positivity [hgroup m]

/-- The signed-harmonic triangular component evaluates without any HPL
input; it is one half of `log² 2 + ζ(2)`. -/
theorem pairedAlternatingSignedHarmonic_hasSum24 :
    HasSum pairedAlternatingSignedHarmonicEulerTerm24
      ((1 / 2) *
        (Real.log 2 ^ 2 + (Real.pi ^ 2 / 6))) := by
  have hpartial (N : ℕ) :
      (∑ m ∈ Finset.range N,
          pairedAlternatingSignedHarmonicEulerTerm24 m) =
        (1 / 2) *
          (signedHarmonic24 (2 * N) ^ 2 +
            harmonicSquare24 (2 * N)) := by
    rw [pairedAlternatingSignedHarmonic_partial,
      signedHarmonic24_triangular_partial]
  have htwo : Tendsto (fun n : ℕ ↦ 2 * n) atTop atTop :=
    (strictMono_mul_left_of_pos (by norm_num : 0 < (2 : ℕ))).tendsto_atTop
  have hharmonicSquare :
      Tendsto (fun N : ℕ ↦ harmonicSquare24 (2 * N)) atTop
        (𝓝 (Real.pi ^ 2 / 6)) := by
    simpa only [harmonicSquare24] using
      shifted_zeta_two_hasSum.tendsto_sum_nat.comp htwo
  have hlimit :
      Tendsto
        (fun N : ℕ ↦
          (1 / 2) *
            (signedHarmonic24 (2 * N) ^ 2 +
              harmonicSquare24 (2 * N)))
        atTop
        (𝓝 ((1 / 2) *
          (Real.log 2 ^ 2 + (Real.pi ^ 2 / 6)))) := by
    convert tendsto_const_nhds.mul
      ((tendsto_signedHarmonic24_even.pow 2).add hharmonicSquare) using 1
    · ring
  rw [hasSum_iff_tendsto_nat_of_nonneg]
  · simpa only [hpartial] using hlimit
  · exact pairedAlternatingSignedHarmonicEulerTerm24_nonneg

theorem shifted_zeta_two_tail (j : ℕ) :
    (∑' m : ℕ, if j < m then 1 / (m + 1 : ℝ) ^ 2 else 0) =
      Real.pi ^ 2 / 6 - harmonicSquare24 (j + 1) := by
  let c : ℕ → ℝ := fun m ↦ 1 / (m + 1 : ℝ) ^ 2
  let p : ℕ → ℝ := fun m ↦ if m ≤ j then c m else 0
  have hc : Summable c := by
    simpa [c] using shifted_zeta_two_hasSum.summable
  have hp : Summable p := by
    apply summable_of_ne_finset_zero (s := Finset.range (j + 1))
    intro m hm
    rw [Finset.mem_range, not_lt] at hm
    have hjm : j < m := by omega
    simp [p, Nat.not_le.mpr hjm]
  have hp_tsum : (∑' m : ℕ, p m) = harmonicSquare24 (j + 1) := by
    calc
      (∑' m : ℕ, p m) =
          ∑ m ∈ Finset.range (j + 1), p m := by
            apply tsum_eq_sum
            intro m hm
            rw [Finset.mem_range, not_lt] at hm
            have hjm : j < m := by omega
            simp [p, Nat.not_le.mpr hjm]
      _ = ∑ m ∈ Finset.range (j + 1), c m := by
        apply Finset.sum_congr rfl
        intro m hm
        have hmj : m ≤ j := Finset.mem_range_succ_iff.mp hm
        simp [p, hmj]
      _ = harmonicSquare24 (j + 1) := by
        rfl
  calc
    (∑' m : ℕ, if j < m then 1 / (m + 1 : ℝ) ^ 2 else 0) =
        ∑' m : ℕ, (c m - p m) := by
          apply tsum_congr
          intro m
          by_cases hjm : j < m
          · simp [c, p, hjm, Nat.not_le.mpr hjm]
          · have hmj : m ≤ j := Nat.le_of_not_gt hjm
            simp [c, p, hjm, hmj]
    _ = (∑' m : ℕ, c m) - ∑' m : ℕ, p m := hc.tsum_sub hp
    _ = Real.pi ^ 2 / 6 - harmonicSquare24 (j + 1) := by
      rw [show (∑' m : ℕ, c m) = Real.pi ^ 2 / 6 by
        simpa [c] using shifted_zeta_two_hasSum.tsum_eq, hp_tsum]

theorem summable_shifted_zeta_two_tail (j : ℕ) :
    Summable (fun m : ℕ ↦
      if j < m then 1 / (m + 1 : ℝ) ^ 2 else 0) := by
  exact shifted_zeta_two_hasSum.summable.of_norm_bounded (fun m ↦ by
    rw [Real.norm_eq_abs]
    by_cases hjm : j < m
    · simp only [hjm, if_true]
      rw [abs_of_nonneg (by positivity)]
    · simp only [hjm, if_false, abs_zero]
      positivity)

theorem inverseCentralSquareSum24_succ (m : ℕ) :
    inverseCentralSquareSum24 (m + 1) =
      inverseCentralSquareSum24 m +
        1 / ((m + 1 : ℝ) ^ 2 *
          (Nat.choose (2 * (m + 1)) (m + 1) : ℝ)) := by
  simp [inverseCentralSquareSum24, Finset.sum_range_succ]

theorem firstMoment24_closed_recurrence (m : ℕ) :
    (m + 1 : ℝ) *
        ((Nat.choose (2 * (m + 1)) (m + 1) : ℝ) *
          harmonicRemainder24 (m + 1)) -
      2 * (2 * (m : ℝ) + 1) *
        ((Nat.choose (2 * m) m : ℝ) * harmonicRemainder24 m) =
      3 * (Nat.choose (2 * m + 1) m : ℝ) -
        2 * (Nat.choose (2 * m) m : ℝ) := by
  rw [centralChoose24_succ, harmonicRemainder24_succ]
  have hrel := centralChoose24_odd_relation m
  have hcentral :
      (Nat.choose (2 * m) m : ℝ) =
        (m + 1 : ℝ) * (Nat.choose (2 * m + 1) m : ℝ) /
          (2 * (m : ℝ) + 1) := by
    rw [eq_div_iff (by positivity)]
    simpa [mul_comm] using hrel.symm
  rw [hcentral]
  field_simp
  ring

/-- The first harmonic moment of the squared-binomial distribution.  This is
the load-bearing lower-moment identity used by the creative-telescoping
recurrence for the challenge's inner sum. -/
theorem firstMoment24_closed (m : ℕ) :
    firstMoment24 m =
      (Nat.choose (2 * m) m : ℝ) * harmonicRemainder24 m := by
  induction m with
  | zero =>
      norm_num [firstMoment24, harmonicRemainder24, harmonicNumber]
  | succ m ih =>
      have hrec := firstMoment24_recurrence m
      rw [ih] at hrec
      have hclosed := firstMoment24_closed_recurrence m
      have hmul :
          (m + 1 : ℝ) * firstMoment24 (m + 1) =
            (m + 1 : ℝ) *
              ((Nat.choose (2 * (m + 1)) (m + 1) : ℝ) *
                harmonicRemainder24 (m + 1)) := by
        linarith
      exact mul_left_cancel₀ (by positivity : (m + 1 : ℝ) ≠ 0) hmul

theorem innerSum24_recurrence (m : ℕ) :
    (m + 1 : ℝ) * innerSum24 (m + 1) -
        2 * (2 * (m : ℝ) + 1) * innerSum24 m =
      (2 * (4 * (m : ℝ) + 1) / (m + 1 : ℝ)) *
          (Nat.choose (2 * m) m : ℝ) * harmonicRemainder24 m +
        (2 * ((m : ℝ) - 1) / (m + 1 : ℝ) ^ 2) *
          (Nat.choose (2 * m) m : ℝ) +
        3 / (m + 1 : ℝ) := by
  rw [innerSum24_recurrence_raw, weightedFirstMoment24_closed,
    firstMoment24_closed, sum_choose_sq_div_succ_sq,
    sum_choose_sq_div_succ, centralChoose24_succ]
  have hrel := centralChoose24_odd_relation m
  have hodd :
      (Nat.choose (2 * m + 1) m : ℝ) =
        (2 * (m : ℝ) + 1) * (Nat.choose (2 * m) m : ℝ) /
          (m + 1 : ℝ) := by
    rw [eq_div_iff (by positivity)]
    simpa [mul_comm] using hrel
  rw [hodd]
  field_simp
  ring

theorem closedInnerSum24_recurrence (m : ℕ) :
    (m + 1 : ℝ) * closedInnerSum24 (m + 1) -
        2 * (2 * (m : ℝ) + 1) * closedInnerSum24 m =
      (2 * (4 * (m : ℝ) + 1) / (m + 1 : ℝ)) *
          (Nat.choose (2 * m) m : ℝ) * harmonicRemainder24 m +
        (2 * ((m : ℝ) - 1) / (m + 1 : ℝ) ^ 2) *
          (Nat.choose (2 * m) m : ℝ) +
        3 / (m + 1 : ℝ) := by
  rw [closedInnerSum24, closedInnerSum24, harmonicRemainder24_succ,
    harmonicSquare24_double_succ, inverseCentralSquareSum24_succ,
    centralChoose24_succ]
  have hrel := centralChoose24_odd_relation m
  have hodd :
      (Nat.choose (2 * m + 1) m : ℝ) =
        (2 * (m : ℝ) + 1) * (Nat.choose (2 * m) m : ℝ) /
          (m + 1 : ℝ) := by
    rw [eq_div_iff (by positivity)]
    simpa [mul_comm] using hrel
  rw [hodd]
  field_simp [(centralChoose24_pos m).ne',
    (centralChoose24_pos (m + 1)).ne']
  ring

/-- Closed form for the finite binomial-harmonic square moment. -/
theorem innerSum24_closed (m : ℕ) :
    innerSum24 m =
      (Nat.choose (2 * m) m : ℝ) *
        (harmonicRemainder24 m ^ 2 - harmonicSquare24 (2 * m) +
          3 * inverseCentralSquareSum24 m) := by
  change innerSum24 m = closedInnerSum24 m
  induction m with
  | zero =>
      norm_num [innerSum24, closedInnerSum24, harmonicRemainder24,
        harmonicNumber, harmonicSquare24, inverseCentralSquareSum24]
  | succ m ih =>
      have hrec := innerSum24_recurrence m
      have hclosed := closedInnerSum24_recurrence m
      rw [ih] at hrec
      have hmul :
          (m + 1 : ℝ) * innerSum24 (m + 1) =
            (m + 1 : ℝ) * closedInnerSum24 (m + 1) := by
        linarith
      exact mul_left_cancel₀ (by positivity : (m + 1 : ℝ) ≠ 0) hmul

theorem inverseCentralSquareSum24_nonneg (m : ℕ) :
    0 ≤ inverseCentralSquareSum24 m := by
  unfold inverseCentralSquareSum24
  positivity

theorem inverseCentralSquareSum24_le_harmonic (m : ℕ) :
    inverseCentralSquareSum24 m ≤ harmonicNumber m := by
  unfold inverseCentralSquareSum24 harmonicNumber
  apply Finset.sum_le_sum
  intro j hj
  have hj_one : (1 : ℝ) ≤ (j + 1 : ℝ) := by exact_mod_cast Nat.succ_le_succ (Nat.zero_le j)
  have hj_pos : 0 < (j + 1 : ℝ) := lt_of_lt_of_le zero_lt_one hj_one
  have hchoose_nat :
      1 ≤ Nat.choose (2 * (j + 1)) (j + 1) :=
    Nat.choose_pos (by omega)
  have hchoose :
      (1 : ℝ) ≤ (Nat.choose (2 * (j + 1)) (j + 1) : ℝ) := by
    exact_mod_cast hchoose_nat
  have hden :
      0 < (j + 1 : ℝ) ^ 2 *
        (Nat.choose (2 * (j + 1)) (j + 1) : ℝ) := by
    positivity
  apply (div_le_div_iff₀ hden hj_pos).2
  have hsq : (j + 1 : ℝ) ≤ (j + 1 : ℝ) ^ 2 := by
    nlinarith [sq_nonneg ((j + 1 : ℝ) - 1)]
  calc
    (1 : ℝ) * (j + 1 : ℝ) = (j + 1 : ℝ) := by ring
    _ ≤ (j + 1 : ℝ) ^ 2 := hsq
    _ ≤ (j + 1 : ℝ) ^ 2 *
        (Nat.choose (2 * (j + 1)) (j + 1) : ℝ) :=
      le_mul_of_one_le_right (sq_nonneg _) hchoose
    _ = (1 : ℝ) *
        ((j + 1 : ℝ) ^ 2 *
          (Nat.choose (2 * (j + 1)) (j + 1) : ℝ)) := by ring

theorem outerTerm24_closed (m : ℕ) :
    outerTerm24 m =
      (harmonicRemainder24 m ^ 2 - harmonicSquare24 (2 * m) +
        3 * inverseCentralSquareSum24 m) / (m + 1 : ℝ) ^ 2 := by
  rw [outerTerm24, innerSum24_closed]
  field_simp [(centralChoose24_pos m).ne']

theorem lhs_24_eq_reduced_series :
    lhs_24 =
      ∑' m : ℕ,
        (harmonicRemainder24 m ^ 2 - harmonicSquare24 (2 * m) +
          3 * inverseCentralSquareSum24 m) / (m + 1 : ℝ) ^ 2 := by
  unfold lhs_24
  apply tsum_congr
  exact outerTerm24_closed

theorem innerSum24_nonneg (m : ℕ) : 0 ≤ innerSum24 m := by
  unfold innerSum24
  positivity

theorem innerSum24_le (m : ℕ) :
    innerSum24 m ≤
      (Nat.choose (2 * m) m : ℝ) * (harmonicNumber m) ^ 2 := by
  unfold innerSum24
  calc
    (∑ k ∈ Finset.range (m + 1),
        (Nat.choose m k : ℝ) ^ 2 * (harmonicNumber k) ^ 2)
        ≤ ∑ k ∈ Finset.range (m + 1),
            (Nat.choose m k : ℝ) ^ 2 * (harmonicNumber m) ^ 2 := by
          gcongr with k hk
          · exact harmonicNumber_nonneg k
          · exact harmonicNumber_mono
              (Nat.le_of_lt_succ (Finset.mem_range.mp hk))
    _ = (∑ k ∈ Finset.range (m + 1), (Nat.choose m k : ℝ) ^ 2) *
          (harmonicNumber m) ^ 2 := by
          rw [Finset.sum_mul]
    _ = (Nat.choose (2 * m) m : ℝ) * (harmonicNumber m) ^ 2 := by
          rw [sum_choose_sq_real]

theorem outerTerm24_nonneg (m : ℕ) : 0 ≤ outerTerm24 m := by
  unfold outerTerm24
  exact div_nonneg (innerSum24_nonneg m) (by positivity [centralChoose24_pos m])

theorem outerTerm24_le_harmonic (m : ℕ) :
    outerTerm24 m ≤ (harmonicNumber m) ^ 2 / (m + 1 : ℝ) ^ 2 := by
  unfold outerTerm24
  have hden :
      0 ≤ (m + 1 : ℝ) ^ 2 * (Nat.choose (2 * m) m : ℝ) := by
    positivity [centralChoose24_pos m]
  calc
    innerSum24 m /
        ((↑m + 1) ^ 2 * (Nat.choose (2 * m) m : ℝ))
        ≤ ((Nat.choose (2 * m) m : ℝ) * (harmonicNumber m) ^ 2) /
            ((↑m + 1) ^ 2 * (Nat.choose (2 * m) m : ℝ)) :=
      div_le_div_of_nonneg_right (innerSum24_le m) hden
    _ = (harmonicNumber m) ^ 2 / (m + 1 : ℝ) ^ 2 := by
      field_simp [(centralChoose24_pos m).ne']

theorem harmonicNumber_le_five_rpow {m : ℕ} (hm : 1 ≤ m) :
    harmonicNumber m ≤ 5 * (m : ℝ) ^ (1 / 4 : ℝ) := by
  have hm_real : (1 : ℝ) ≤ m := by exact_mod_cast hm
  have hrpow_one : (1 : ℝ) ≤ (m : ℝ) ^ (1 / 4 : ℝ) :=
    Real.one_le_rpow hm_real (by norm_num)
  have hlog :=
    Real.log_natCast_le_rpow_div m (by norm_num : (0 : ℝ) < 1 / 4)
  rw [harmonicNumber_eq_cast_harmonic]
  calc
    (harmonic m : ℝ) ≤ 1 + Real.log (m : ℝ) :=
      harmonic_le_one_add_log m
    _ ≤ 1 + 4 * (m : ℝ) ^ (1 / 4 : ℝ) := by
      norm_num at hlog ⊢
      linarith
    _ ≤ 5 * (m : ℝ) ^ (1 / 4 : ℝ) := by
      linarith

theorem harmonicNumber_sq_le_twentyfive_rpow {m : ℕ} (hm : 1 ≤ m) :
    (harmonicNumber m) ^ 2 ≤ 25 * (m : ℝ) ^ (1 / 2 : ℝ) := by
  have hpow_nonneg : 0 ≤ (m : ℝ) ^ (1 / 4 : ℝ) :=
    Real.rpow_nonneg (Nat.cast_nonneg m) _
  have hsquare :
      (harmonicNumber m) ^ 2 ≤
        (5 * (m : ℝ) ^ (1 / 4 : ℝ)) ^ 2 :=
    (sq_le_sq₀ (harmonicNumber_nonneg m) (by positivity)).2
      (harmonicNumber_le_five_rpow hm)
  have hpow_sq :
      ((m : ℝ) ^ (1 / 4 : ℝ)) ^ 2 =
        (m : ℝ) ^ (1 / 2 : ℝ) := by
    rw [← Real.rpow_mul_natCast (Nat.cast_nonneg m)]
    norm_num
  calc
    (harmonicNumber m) ^ 2
        ≤ (5 * (m : ℝ) ^ (1 / 4 : ℝ)) ^ 2 := hsquare
    _ = 25 * (m : ℝ) ^ (1 / 2 : ℝ) := by
      rw [mul_pow, hpow_sq]
      norm_num

theorem harmonicNumber_sq_div_le_pseries (m : ℕ) :
    (harmonicNumber m) ^ 2 / (m + 1 : ℝ) ^ 2 ≤
      25 / |(m : ℝ) + 1| ^ (3 / 2 : ℝ) := by
  by_cases hm0 : m = 0
  · simp [hm0, harmonicNumber_zero]
  have hm : 1 ≤ m := Nat.one_le_iff_ne_zero.mpr hm0
  have hm_le : (m : ℝ) ≤ (m : ℝ) + 1 := by linarith
  have hrpow_le :
      (m : ℝ) ^ (1 / 2 : ℝ) ≤ ((m : ℝ) + 1) ^ (1 / 2 : ℝ) :=
    Real.rpow_le_rpow (Nat.cast_nonneg m) hm_le (by norm_num)
  have hx : 0 < (m : ℝ) + 1 := by positivity
  calc
    (harmonicNumber m) ^ 2 / (m + 1 : ℝ) ^ 2
        ≤ (25 * (m : ℝ) ^ (1 / 2 : ℝ)) / (m + 1 : ℝ) ^ 2 :=
      div_le_div_of_nonneg_right
        (harmonicNumber_sq_le_twentyfive_rpow hm) (by positivity)
    _ ≤ (25 * ((m : ℝ) + 1) ^ (1 / 2 : ℝ)) /
          ((m : ℝ) + 1) ^ 2 := by
      gcongr
    _ = 25 / |(m : ℝ) + 1| ^ (3 / 2 : ℝ) := by
      rw [abs_of_pos hx]
      field_simp [(Real.rpow_pos_of_pos hx (1 / 2 : ℝ)).ne',
        (Real.rpow_pos_of_pos hx (3 / 2 : ℝ)).ne']
      rw [← Real.rpow_add hx]
      norm_num [Real.rpow_natCast]

theorem summable_harmonicNumber_sq_div :
    Summable (fun m : ℕ ↦
      (harmonicNumber m) ^ 2 / (m + 1 : ℝ) ^ 2) := by
  have hp :
      Summable (fun m : ℕ ↦
        25 / |(m : ℝ) + 1| ^ (3 / 2 : ℝ)) := by
    have hbase :=
      (Real.summable_one_div_nat_add_rpow 1 (3 / 2 : ℝ)).2 (by norm_num)
    simpa [div_eq_mul_inv] using hbase.mul_left 25
  exact hp.of_norm_bounded (fun m ↦ by
    rw [Real.norm_eq_abs, abs_of_nonneg (by positivity)]
    exact harmonicNumber_sq_div_le_pseries m)

theorem harmonicNumber_succ_sq_div_le_pseries (n : ℕ) :
    harmonicNumber (n + 1) ^ 2 / (n + 1 : ℝ) ^ 2 ≤
      25 / |(n : ℝ) + 1| ^ (3 / 2 : ℝ) := by
  have hbound :
      harmonicNumber (n + 1) ^ 2 ≤
        25 * ((n : ℝ) + 1) ^ (1 / 2 : ℝ) := by
    simpa using
      (harmonicNumber_sq_le_twentyfive_rpow
        (m := n + 1) (by omega : 1 ≤ n + 1))
  have hx : 0 < (n : ℝ) + 1 := by positivity
  calc
    harmonicNumber (n + 1) ^ 2 / (n + 1 : ℝ) ^ 2 ≤
        (25 * ((n : ℝ) + 1) ^ (1 / 2 : ℝ)) /
          ((n : ℝ) + 1) ^ 2 :=
      div_le_div_of_nonneg_right hbound (by positivity)
    _ = 25 / |(n : ℝ) + 1| ^ (3 / 2 : ℝ) := by
      rw [abs_of_pos hx]
      field_simp [(Real.rpow_pos_of_pos hx (1 / 2 : ℝ)).ne',
        (Real.rpow_pos_of_pos hx (3 / 2 : ℝ)).ne']
      rw [← Real.rpow_add hx]
      norm_num [Real.rpow_natCast]

theorem summable_harmonicNumber_succ_sq_div :
    Summable (fun n : ℕ ↦
      harmonicNumber (n + 1) ^ 2 / (n + 1 : ℝ) ^ 2) := by
  have hp :
      Summable (fun n : ℕ ↦
        25 / |(n : ℝ) + 1| ^ (3 / 2 : ℝ)) := by
    have hbase :=
      (Real.summable_one_div_nat_add_rpow 1 (3 / 2 : ℝ)).2 (by norm_num)
    simpa [div_eq_mul_inv] using hbase.mul_left 25
  exact hp.of_norm_bounded (fun n ↦ by
    rw [Real.norm_eq_abs, abs_of_nonneg (by positivity)]
    exact harmonicNumber_succ_sq_div_le_pseries n)

theorem one_le_harmonicNumber_succ (n : ℕ) :
    1 ≤ harmonicNumber (n + 1) := by
  simpa [harmonicNumber_one] using
    (harmonicNumber_mono (k := 1) (m := n + 1) (by omega))

theorem harmonicNumber_succ_le_sq (n : ℕ) :
    harmonicNumber (n + 1) ≤ harmonicNumber (n + 1) ^ 2 := by
  have h := one_le_harmonicNumber_succ n
  nlinarith

theorem quadraticLinearEulerTerm24_norm_le (n : ℕ) :
    ‖quadraticLinearEulerTerm24 n‖ ≤
      3 * (harmonicNumber (n + 1) ^ 2 / (n + 1 : ℝ) ^ 2) := by
  unfold quadraticLinearEulerTerm24
  rw [Real.norm_eq_abs, abs_div, abs_pow,
    abs_of_pos (by positivity : (0 : ℝ) < n + 1)]
  calc
    |parityRemainder24 (n + 1)| / (n + 1 : ℝ) ^ 2 ≤
        (3 * harmonicNumber (n + 1)) / (n + 1 : ℝ) ^ 2 :=
      div_le_div_of_nonneg_right
        (abs_parityRemainder24_le (n + 1)) (by positivity)
    _ ≤ (3 * harmonicNumber (n + 1) ^ 2) /
          (n + 1 : ℝ) ^ 2 := by
      gcongr
      exact harmonicNumber_succ_le_sq n
    _ = 3 * (harmonicNumber (n + 1) ^ 2 /
          (n + 1 : ℝ) ^ 2) := by ring

theorem summable_quadraticLinearEulerTerm24 :
    Summable quadraticLinearEulerTerm24 := by
  exact (summable_harmonicNumber_succ_sq_div.mul_left 3).of_norm_bounded
    quadraticLinearEulerTerm24_norm_le

theorem cubicLinearEulerTerm24_norm_le (n : ℕ) :
    ‖cubicLinearEulerTerm24 n‖ ≤
      3 * (harmonicNumber (n + 1) ^ 2 / (n + 1 : ℝ) ^ 2) := by
  unfold cubicLinearEulerTerm24
  rw [Real.norm_eq_abs, abs_div, abs_pow,
    abs_of_pos (by positivity : (0 : ℝ) < n + 1)]
  calc
    |parityRemainder24 (n + 1)| / (n + 1 : ℝ) ^ 3 ≤
        (3 * harmonicNumber (n + 1)) / (n + 1 : ℝ) ^ 3 :=
      div_le_div_of_nonneg_right
        (abs_parityRemainder24_le (n + 1)) (by positivity)
    _ ≤ (3 * harmonicNumber (n + 1) ^ 2) /
          (n + 1 : ℝ) ^ 3 := by
      gcongr
      exact harmonicNumber_succ_le_sq n
    _ ≤ (3 * harmonicNumber (n + 1) ^ 2) /
          (n + 1 : ℝ) ^ 2 := by
      gcongr
      all_goals norm_num
    _ = 3 * (harmonicNumber (n + 1) ^ 2 /
          (n + 1 : ℝ) ^ 2) := by ring

theorem summable_cubicLinearEulerTerm24 :
    Summable cubicLinearEulerTerm24 := by
  exact (summable_harmonicNumber_succ_sq_div.mul_left 3).of_norm_bounded
    cubicLinearEulerTerm24_norm_le

theorem quadraticEulerTerm24_norm_le (n : ℕ) :
    ‖quadraticEulerTerm24 n‖ ≤
      10 * (harmonicNumber (n + 1) ^ 2 / (n + 1 : ℝ) ^ 2) := by
  have hb := abs_parityRemainder24_le (n + 1)
  have hb_sq :
      parityRemainder24 (n + 1) ^ 2 ≤
        9 * harmonicNumber (n + 1) ^ 2 := by
    have hsquare :=
      (sq_le_sq₀ (abs_nonneg (parityRemainder24 (n + 1)))
        (mul_nonneg (by norm_num) (harmonicNumber_nonneg (n + 1)))).2 hb
    calc
      parityRemainder24 (n + 1) ^ 2 =
          |parityRemainder24 (n + 1)| ^ 2 := by rw [sq_abs]
      _ ≤ (3 * harmonicNumber (n + 1)) ^ 2 := hsquare
      _ = 9 * harmonicNumber (n + 1) ^ 2 := by ring
  have hsquare_le :
      harmonicSquare24 (n + 1) ≤ harmonicNumber (n + 1) ^ 2 :=
    (harmonicSquare24_le_harmonicNumber (n + 1)).trans
      (harmonicNumber_succ_le_sq n)
  have hsquare_nonneg : 0 ≤ harmonicSquare24 (n + 1) := by
    unfold harmonicSquare24
    positivity
  unfold quadraticEulerTerm24
  rw [Real.norm_eq_abs, abs_div, abs_pow,
    abs_of_pos (by positivity : (0 : ℝ) < n + 1)]
  calc
    |parityRemainder24 (n + 1) ^ 2 -
        harmonicSquare24 (n + 1)| / (n + 1 : ℝ) ^ 2 ≤
        (parityRemainder24 (n + 1) ^ 2 +
          harmonicSquare24 (n + 1)) / (n + 1 : ℝ) ^ 2 := by
      gcongr
      rw [abs_sub_le_iff]
      constructor <;> nlinarith [sq_nonneg
        (parityRemainder24 (n + 1)),
        hsquare_nonneg]
    _ ≤ (10 * harmonicNumber (n + 1) ^ 2) /
          (n + 1 : ℝ) ^ 2 := by
      gcongr
      nlinarith
    _ = 10 * (harmonicNumber (n + 1) ^ 2 /
          (n + 1 : ℝ) ^ 2) := by ring

theorem summable_quadraticEulerTerm24 :
    Summable quadraticEulerTerm24 := by
  exact (summable_harmonicNumber_succ_sq_div.mul_left 10).of_norm_bounded
    quadraticEulerTerm24_norm_le

theorem summable_outer_alternating24 {f : ℕ → ℝ} (hf : Summable f) :
    Summable (fun n : ℕ ↦ (-1 : ℝ) ^ (n + 1) * f n) := by
  apply summable_norm_iff.mp
  exact (summable_norm_iff.mpr hf).congr (fun n ↦ by
    rw [norm_mul, norm_pow]
    norm_num)

theorem summable_alternatingQuadraticEulerTerm24 :
    Summable alternatingQuadraticEulerTerm24 := by
  unfold alternatingQuadraticEulerTerm24
  exact summable_outer_alternating24 summable_quadraticEulerTerm24

theorem summable_alternatingCubicLinearEulerTerm24 :
    Summable alternatingCubicLinearEulerTerm24 := by
  unfold alternatingCubicLinearEulerTerm24
  exact summable_outer_alternating24 summable_cubicLinearEulerTerm24

theorem summable_alternatingQuadraticLinearEulerTerm24 :
    Summable alternatingQuadraticLinearEulerTerm24 := by
  unfold alternatingQuadraticLinearEulerTerm24
  exact summable_outer_alternating24 summable_quadraticLinearEulerTerm24

theorem pairedAlternatingHarmonicEulerTerm24_norm_le (m : ℕ) :
    ‖pairedAlternatingHarmonicEulerTerm24 m‖ ≤
      3 * (harmonicNumber (m + 1) ^ 2 / (m + 1 : ℝ) ^ 2) := by
  let g : ℝ :=
    1 / (2 * m + 1 : ℝ) - 1 / (2 * m + 2 : ℝ)
  have hg_nonneg : 0 ≤ g := by
    unfold g
    apply sub_nonneg.mpr
    apply one_div_le_one_div_of_le
    · positivity
    · norm_num
  have hg_formula :
      g = 1 / ((2 * m + 1 : ℝ) * (2 * m + 2 : ℝ)) := by
    unfold g
    field_simp
    ring
  have hg_le : g ≤ 1 / (m + 1 : ℝ) ^ 2 := by
    rw [hg_formula]
    apply one_div_le_one_div_of_le (by positivity)
    nlinarith [sq_nonneg (m : ℝ)]
  have hHodd :
      harmonicNumber (2 * m + 1) ≤ 2 * harmonicNumber (m + 1) :=
    (harmonicNumber_mono (by omega : 2 * m + 1 ≤ 2 * (m + 1))).trans
      (harmonicNumber_two_mul_le (m + 1))
  have heven_le :
      1 / (2 * m + 2 : ℝ) ^ 2 ≤ 1 / (m + 1 : ℝ) ^ 2 := by
    apply one_div_le_one_div_of_le (by positivity)
    nlinarith [sq_nonneg (m : ℝ)]
  rw [pairedAlternatingHarmonicEulerTerm24_formula,
    Real.norm_eq_abs]
  change
    |-harmonicNumber (2 * m + 1) * g +
        1 / (2 * m + 2 : ℝ) ^ 2| ≤ _
  calc
    |-harmonicNumber (2 * m + 1) * g +
        1 / (2 * m + 2 : ℝ) ^ 2| ≤
        harmonicNumber (2 * m + 1) * g +
          1 / (2 * m + 2 : ℝ) ^ 2 := by
      calc
        _ ≤ |-harmonicNumber (2 * m + 1) * g| +
              |1 / (2 * m + 2 : ℝ) ^ 2| := abs_add_le _ _
        _ = _ := by
          rw [abs_mul, abs_neg, abs_of_nonneg
            (harmonicNumber_nonneg (2 * m + 1)),
            abs_of_nonneg hg_nonneg,
            abs_of_nonneg (by positivity)]
    _ ≤ (2 * harmonicNumber (m + 1)) *
          (1 / (m + 1 : ℝ) ^ 2) +
          1 / (m + 1 : ℝ) ^ 2 := by
      apply add_le_add
      · exact mul_le_mul hHodd hg_le hg_nonneg
          (mul_nonneg (by norm_num) (harmonicNumber_nonneg (m + 1)))
      · exact heven_le
    _ ≤ (2 * harmonicNumber (m + 1) ^ 2) /
          (m + 1 : ℝ) ^ 2 +
          harmonicNumber (m + 1) ^ 2 / (m + 1 : ℝ) ^ 2 := by
      have hH := one_le_harmonicNumber_succ m
      have hHsq := harmonicNumber_succ_le_sq m
      apply add_le_add
      · calc
          (2 * harmonicNumber (m + 1)) *
                (1 / (m + 1 : ℝ) ^ 2) =
              (2 * harmonicNumber (m + 1)) / (m + 1 : ℝ) ^ 2 := by
                ring
          _ ≤ (2 * harmonicNumber (m + 1) ^ 2) /
                (m + 1 : ℝ) ^ 2 :=
            div_le_div_of_nonneg_right (by nlinarith) (by positivity)
      · exact div_le_div_of_nonneg_right (by nlinarith) (by positivity)
    _ = 3 * (harmonicNumber (m + 1) ^ 2 /
          (m + 1 : ℝ) ^ 2) := by ring

theorem summable_pairedAlternatingHarmonicEulerTerm24 :
    Summable pairedAlternatingHarmonicEulerTerm24 := by
  exact (summable_harmonicNumber_succ_sq_div.mul_left 3).of_norm_bounded
    pairedAlternatingHarmonicEulerTerm24_norm_le

theorem inverseCentralSquareSum24_div_le_harmonic_sq (m : ℕ) :
    inverseCentralSquareSum24 m / (m + 1 : ℝ) ^ 2 ≤
      harmonicNumber m ^ 2 / (m + 1 : ℝ) ^ 2 := by
  by_cases hm0 : m = 0
  · simp [hm0, inverseCentralSquareSum24, harmonicNumber]
  have hm : 1 ≤ m := Nat.one_le_iff_ne_zero.mpr hm0
  have hH_one : (1 : ℝ) ≤ harmonicNumber m := by
    simpa [harmonicNumber_one] using
      (harmonicNumber_mono (k := 1) (m := m) hm)
  have hH_sq : harmonicNumber m ≤ harmonicNumber m ^ 2 := by
    nlinarith
  exact div_le_div_of_nonneg_right
    ((inverseCentralSquareSum24_le_harmonic m).trans hH_sq) (by positivity)

theorem summable_inverseCentralSquareSum24_div :
    Summable (fun m : ℕ ↦
      inverseCentralSquareSum24 m / (m + 1 : ℝ) ^ 2) := by
  exact summable_harmonicNumber_sq_div.of_norm_bounded (fun m ↦ by
    rw [Real.norm_eq_abs,
      abs_of_nonneg (div_nonneg (inverseCentralSquareSum24_nonneg m) (by positivity))]
    exact inverseCentralSquareSum24_div_le_harmonic_sq m)

theorem problem24_series_summable : Summable outerTerm24 := by
  exact summable_harmonicNumber_sq_div.of_norm_bounded (fun m ↦ by
    rw [Real.norm_eq_abs, abs_of_nonneg (outerTerm24_nonneg m)]
    exact outerTerm24_le_harmonic m)

theorem outerTerm24_decomposition (m : ℕ) :
    outerTerm24 m =
      elementaryOuterTerm24 m + 3 * inverseCentralOuterTerm24 m := by
  rw [outerTerm24_closed]
  unfold elementaryOuterTerm24 inverseCentralOuterTerm24
  ring

theorem summable_inverseCentralOuterTerm24 :
    Summable inverseCentralOuterTerm24 := by
  simpa [inverseCentralOuterTerm24] using
    summable_inverseCentralSquareSum24_div

theorem inverseCentralSquareSum24_eq_sum (m : ℕ) :
    inverseCentralSquareSum24 m =
      ∑ j ∈ Finset.range m, inverseCentralCoefficient24 j := by
  rfl

theorem inverseCentralCoefficient24_nonneg (j : ℕ) :
    0 ≤ inverseCentralCoefficient24 j := by
  unfold inverseCentralCoefficient24
  positivity

/-- Problem 2.6 supplies a complete beta-kernel/FTC proof of the classical
weight-two inverse central-binomial sum.  Its zero-index term vanishes, so
shifting the index gives exactly the coefficient sequence used here. -/
theorem inverseCentralCoefficient24_hasSum :
    HasSum inverseCentralCoefficient24 (Real.pi ^ 2 / 6 / 3) := by
  have h :=
    (hasSum_nat_add_iff' 1).2
      RamanujanChallenge.P26.inverseBinomialC26_hasSum
  convert h using 1
  · funext j
    unfold inverseCentralCoefficient24
      RamanujanChallenge.P26.inverseBinomialC26
    push_cast
    rfl
  · simp [RamanujanChallenge.P26.inverseBinomialC26]

theorem harmonicSquare24_le_zeta_two (n : ℕ) :
    harmonicSquare24 n ≤ Real.pi ^ 2 / 6 := by
  unfold harmonicSquare24
  rw [← shifted_zeta_two_hasSum.tsum_eq]
  exact shifted_zeta_two_hasSum.summable.sum_le_tsum
    (Finset.range n) (fun j hj ↦ by positivity)

theorem summable_inverseCentralFourthCoefficient24 :
    Summable inverseCentralFourthCoefficient24 := by
  exact inverseCentralCoefficient24_hasSum.summable.of_norm_bounded
    (fun j ↦ by
      rw [Real.norm_eq_abs,
        abs_of_nonneg (by
          unfold inverseCentralFourthCoefficient24
          positivity [inverseCentralCoefficient24_nonneg j])]
      unfold inverseCentralFourthCoefficient24
      exact div_le_self (inverseCentralCoefficient24_nonneg j)
        (one_le_pow₀ (by
          exact_mod_cast Nat.succ_le_succ (Nat.zero_le j))))

theorem summable_harmonicSquare_mul_inverseCentralCoefficient24 :
    Summable (fun j : ℕ ↦
      harmonicSquare24 j * inverseCentralCoefficient24 j) := by
  exact
    (inverseCentralCoefficient24_hasSum.summable.mul_left
      (Real.pi ^ 2 / 6)).of_norm_bounded (fun j ↦ by
        rw [Real.norm_eq_abs, abs_of_nonneg
          (mul_nonneg (by
            unfold harmonicSquare24
            positivity) (inverseCentralCoefficient24_nonneg j))]
        exact mul_le_mul_of_nonneg_right
          (harmonicSquare24_le_zeta_two j)
          (inverseCentralCoefficient24_nonneg j))

theorem summable_bbbWeightFourTerm24 :
    Summable bbbWeightFourTerm24 := by
  have hFourth :=
    summable_inverseCentralFourthCoefficient24.mul_left (3 : ℝ)
  have hSquare :=
    summable_harmonicSquare_mul_inverseCentralCoefficient24.mul_left (9 : ℝ)
  apply (hFourth.sub hSquare).congr
  intro j
  unfold bbbWeightFourTerm24
  ring

theorem summable_leshchinerWeightFourTerm24 :
    Summable leshchinerWeightFourTerm24 := by
  have hFourth :=
    summable_inverseCentralFourthCoefficient24.mul_left (2 : ℝ)
  have hSquare :=
    summable_harmonicSquare_mul_inverseCentralCoefficient24.mul_left
      (3 / 2 : ℝ)
  apply (hFourth.sub hSquare).congr
  intro j
  unfold leshchinerWeightFourTerm24
  ring

theorem inverseCentralDoubleTerm24_nonneg (m j : ℕ) :
    0 ≤ inverseCentralDoubleTerm24 m j := by
  unfold inverseCentralDoubleTerm24
  split_ifs
  · positivity [inverseCentralCoefficient24_nonneg j]
  · exact le_rfl

theorem summable_inverseCentralDoubleTerm24_row (m : ℕ) :
    Summable (fun j : ℕ ↦ inverseCentralDoubleTerm24 m j) := by
  apply summable_of_ne_finset_zero (s := Finset.range m)
  intro j hj
  rw [Finset.mem_range, not_lt] at hj
  simp [inverseCentralDoubleTerm24, Nat.not_lt.mpr hj]

theorem tsum_inverseCentralDoubleTerm24_row (m : ℕ) :
    (∑' j : ℕ, inverseCentralDoubleTerm24 m j) =
      inverseCentralOuterTerm24 m := by
  calc
    (∑' j : ℕ, inverseCentralDoubleTerm24 m j) =
        ∑ j ∈ Finset.range m, inverseCentralDoubleTerm24 m j := by
          apply tsum_eq_sum
          intro j hj
          rw [Finset.mem_range, not_lt] at hj
          simp [inverseCentralDoubleTerm24, Nat.not_lt.mpr hj]
    _ = ∑ j ∈ Finset.range m,
          inverseCentralCoefficient24 j / (m + 1 : ℝ) ^ 2 := by
      apply Finset.sum_congr rfl
      intro j hj
      simp [inverseCentralDoubleTerm24, Finset.mem_range.mp hj]
    _ = inverseCentralOuterTerm24 m := by
      rw [inverseCentralOuterTerm24, inverseCentralSquareSum24_eq_sum,
        Finset.sum_div]

theorem tsum_inverseCentralDoubleTerm24_column (j : ℕ) :
    (∑' m : ℕ, inverseCentralDoubleTerm24 m j) =
      inverseCentralCoefficient24 j *
        (Real.pi ^ 2 / 6 - harmonicSquare24 (j + 1)) := by
  calc
    (∑' m : ℕ, inverseCentralDoubleTerm24 m j) =
        ∑' m : ℕ, inverseCentralCoefficient24 j *
          (if j < m then 1 / (m + 1 : ℝ) ^ 2 else 0) := by
            apply tsum_congr
            intro m
            by_cases hjm : j < m
            · simp [inverseCentralDoubleTerm24, hjm]
              ring
            · simp [inverseCentralDoubleTerm24, hjm]
    _ = inverseCentralCoefficient24 j *
        ∑' m : ℕ, (if j < m then 1 / (m + 1 : ℝ) ^ 2 else 0) :=
      Summable.tsum_mul_left (inverseCentralCoefficient24 j)
        (summable_shifted_zeta_two_tail j)
    _ = inverseCentralCoefficient24 j *
        (Real.pi ^ 2 / 6 - harmonicSquare24 (j + 1)) := by
      rw [shifted_zeta_two_tail]

theorem summable_inverseCentralDoubleTerm24 :
    Summable (Function.uncurry inverseCentralDoubleTerm24) := by
  change Summable (fun p : ℕ × ℕ ↦ inverseCentralDoubleTerm24 p.1 p.2)
  rw [summable_prod_of_nonneg (fun p ↦
    inverseCentralDoubleTerm24_nonneg p.1 p.2)]
  exact ⟨summable_inverseCentralDoubleTerm24_row, by
    simpa only [tsum_inverseCentralDoubleTerm24_row] using
      summable_inverseCentralOuterTerm24⟩

theorem summable_inverseCentralTailProduct24 :
    Summable (fun j : ℕ ↦ inverseCentralCoefficient24 j *
      (Real.pi ^ 2 / 6 - harmonicSquare24 (j + 1))) := by
  have hswap :
      Summable (fun p : ℕ × ℕ ↦
        inverseCentralDoubleTerm24 p.2 p.1) := by
    simpa [Function.uncurry] using
      summable_inverseCentralDoubleTerm24.prod_symm
  have hfibers :=
    (summable_prod_of_nonneg (f := fun p : ℕ × ℕ ↦
      inverseCentralDoubleTerm24 p.2 p.1) (fun p ↦
        inverseCentralDoubleTerm24_nonneg p.2 p.1)).mp hswap
  exact hfibers.2.congr (fun j ↦ tsum_inverseCentralDoubleTerm24_column j)

theorem inverseCentralOuterSeries_rearranged :
    (∑' m : ℕ, inverseCentralOuterTerm24 m) =
      ∑' j : ℕ, inverseCentralCoefficient24 j *
        (Real.pi ^ 2 / 6 - harmonicSquare24 (j + 1)) := by
  calc
    (∑' m : ℕ, inverseCentralOuterTerm24 m) =
        ∑' m : ℕ, ∑' j : ℕ, inverseCentralDoubleTerm24 m j := by
      apply tsum_congr
      intro m
      exact (tsum_inverseCentralDoubleTerm24_row m).symm
    _ = ∑' j : ℕ, ∑' m : ℕ, inverseCentralDoubleTerm24 m j :=
      summable_inverseCentralDoubleTerm24.tsum_comm.symm
    _ = ∑' j : ℕ, inverseCentralCoefficient24 j *
        (Real.pi ^ 2 / 6 - harmonicSquare24 (j + 1)) := by
      apply tsum_congr
      exact tsum_inverseCentralDoubleTerm24_column

theorem summable_elementaryOuterTerm24 :
    Summable elementaryOuterTerm24 := by
  have h :=
    problem24_series_summable.sub
      (Summable.mul_left 3 summable_inverseCentralOuterTerm24)
  exact h.congr (fun m ↦ by
    rw [outerTerm24_decomposition]
    ring)

theorem lhs_24_eq_decomposed_series :
    lhs_24 =
      (∑' m : ℕ, elementaryOuterTerm24 m) +
        3 * ∑' m : ℕ, inverseCentralOuterTerm24 m := by
  unfold lhs_24
  calc
    (∑' m : ℕ, outerTerm24 m) =
        ∑' m : ℕ,
          (elementaryOuterTerm24 m + 3 * inverseCentralOuterTerm24 m) :=
      tsum_congr outerTerm24_decomposition
    _ = (∑' m : ℕ, elementaryOuterTerm24 m) +
        ∑' m : ℕ, 3 * inverseCentralOuterTerm24 m :=
      (summable_elementaryOuterTerm24.tsum_add
        (Summable.mul_left 3 summable_inverseCentralOuterTerm24))
    _ = (∑' m : ℕ, elementaryOuterTerm24 m) +
        3 * ∑' m : ℕ, inverseCentralOuterTerm24 m := by
      rw [Summable.tsum_mul_left 3 summable_inverseCentralOuterTerm24]

theorem lhs_24_eq_elementary_add_inverse_central :
    lhs_24 =
      (∑' m : ℕ, elementaryOuterTerm24 m) +
        3 * ∑' j : ℕ, inverseCentralCoefficient24 j *
          (Real.pi ^ 2 / 6 - harmonicSquare24 (j + 1)) := by
  rw [lhs_24_eq_decomposed_series, inverseCentralOuterSeries_rearranged]

theorem outerTerm24_eq_challenge_term (m : ℕ) :
    outerTerm24 m =
      ∑ k ∈ Finset.range (m + 1),
        (Nat.choose m k : ℝ) ^ 2 * (harmonicNumber k) ^ 2 /
          ((m + 1 : ℝ) ^ 2 * (Nat.choose (2 * m) m : ℝ)) := by
  simp only [outerTerm24, innerSum24, Finset.sum_div]

theorem problem24_series_hasSum :
    HasSum
      (fun m : ℕ ↦
        ∑ k ∈ Finset.range (m + 1),
          (Nat.choose m k : ℝ) ^ 2 * (harmonicNumber k) ^ 2 /
            ((m + 1 : ℝ) ^ 2 * (Nat.choose (2 * m) m : ℝ)))
      lhs_24 := by
  simpa only [← outerTerm24_eq_challenge_term] using problem24_series_summable.hasSum

/-! ## The RHS -/

def zeta3_24 : ℝ := ∑' n : ℕ, (1 : ℝ) / (↑n + 1) ^ 3

theorem shifted_zeta_three_hasSum24 :
    HasSum (fun n : ℕ ↦ 1 / (n + 1 : ℝ) ^ 3) zeta3_24 := by
  have hs0 : Summable (fun n : ℕ ↦ 1 / (n : ℝ) ^ 3) :=
    summable_one_div_nat_pow.mpr (by norm_num)
  have hs :
      Summable (fun n : ℕ ↦ 1 / (n + 1 : ℝ) ^ 3) := by
    simpa using (summable_nat_add_iff 1).2 hs0
  unfold zeta3_24
  exact hs.hasSum

def rhs_24 : ℝ :=
  20 * polylog4 (1/2) + (5/6) * (Real.log 2) ^ 4 +
  10 * (Real.pi ^ 2 / 6) - (65/9) * (Real.pi ^ 2 / 6) ^ 2 -
  (Real.log 2) ^ 2 * (12 + 5 * (Real.pi ^ 2 / 6)) +
  (1/2) * zeta3_24 +
  Real.log 2 * ((35/2) * zeta3_24 - 16)

def problem24Statement : Prop := lhs_24 = rhs_24

/-! ## Exact interface to the two remaining scalar evaluations -/

def elementarySeriesValue24 : ℝ :=
  20 * polylog4 (1/2) + (5/6) * (Real.log 2) ^ 4 +
  10 * (Real.pi ^ 2 / 6) - (38/5) * (Real.pi ^ 2 / 6) ^ 2 -
  (Real.log 2) ^ 2 * (12 + 5 * (Real.pi ^ 2 / 6)) +
  (1/2) * zeta3_24 +
  Real.log 2 * ((35/2) * zeta3_24 - 16)

def inverseCentralSeriesValue24 : ℝ :=
  (17/135) * (Real.pi ^ 2 / 6) ^ 2

/-- BBV reduction of the nonalternating quadratic Euler sum occurring in
the elementary component. -/
def quadraticEulerValue24 : ℝ :=
  20 * polylog4 (1 / 2) + (5 / 6) * Real.log 2 ^ 4 +
    7 * Real.log 2 ^ 2 * (Real.pi ^ 2 / 6) -
    (59 / 10) * (Real.pi ^ 2 / 6) ^ 2

/-- BBV reduction of the outer-alternating quadratic Euler sum. -/
def alternatingQuadraticEulerValue24 : ℝ :=
  -22 * polylog4 (1 / 2) - (11 / 12) * Real.log 2 ^ 4 -
    (13 / 2) * Real.log 2 ^ 2 * (Real.pi ^ 2 / 6) -
    (7 / 4) * Real.log 2 * zeta3_24 +
    (67 / 10) * (Real.pi ^ 2 / 6) ^ 2

/-- Classical reduction of the nonalternating linear weight-four sum. -/
def cubicLinearEulerValue24 : ℝ :=
  -(7 / 2) * Real.log 2 * zeta3_24 +
    (3 / 4) * (Real.pi ^ 2 / 6) ^ 2

/-- Classical reduction of its outer-alternating companion. -/
def alternatingCubicLinearEulerValue24 : ℝ :=
  -2 * polylog4 (1 / 2) - (1 / 12) * Real.log 2 ^ 4 +
    (1 / 2) * Real.log 2 ^ 2 * (Real.pi ^ 2 / 6) +
    (7 / 4) * Real.log 2 * zeta3_24 +
    (1 / 10) * (Real.pi ^ 2 / 6) ^ 2

/-- The shifted linear Euler sum obtained after pairing odd and even
denominators. -/
def shiftedLinearEulerValue24 : ℝ :=
  -(3 / 2) * Real.log 2 ^ 2 - (1 / 2) * (Real.pi ^ 2 / 6) +
    3 * Real.log 2 - (5 / 16) * zeta3_24

/-- The alternating weight-two linear Euler sum used to evaluate the
shifted term. -/
def alternatingLinearEulerValue24 : ℝ :=
  (3 / 2) * Real.log 2 ^ 2 + (1 / 2) * (Real.pi ^ 2 / 6)

/-- The ordinary-harmonic part of the paired alternating weight-two sum. -/
def pairedAlternatingHarmonicEulerValue24 : ℝ :=
  (1 / 2) * Real.log 2 ^ 2 - (1 / 2) * (Real.pi ^ 2 / 6)

/-- The nonalternating weight-three linear Euler sum. -/
def quadraticLinearEulerValue24 : ℝ :=
  -3 * Real.log 2 * (Real.pi ^ 2 / 6) + (5 / 2) * zeta3_24

/-- The alternating weight-three linear Euler sum. -/
def alternatingQuadraticLinearEulerValue24 : ℝ :=
  3 * Real.log 2 * (Real.pi ^ 2 / 6) - (15 / 8) * zeta3_24

/-- The elementary rational tail accompanying the Euler sums. -/
def rationalCorrectionValue24 : ℝ :=
  -40 * Real.log 2 + 14 * (Real.pi ^ 2 / 6) +
    3 * zeta3_24 + (Real.pi ^ 2 / 6) ^ 2

/-- The rational correction has an unconditional evaluation: after partial
fractions it is a grouped alternating harmonic series plus ordinary
`ζ(2)`, `ζ(3)`, and `ζ(4)` series. -/
theorem rationalCorrection_hasSum24 :
    HasSum rationalCorrectionTerm24 rationalCorrectionValue24 := by
  have hCombined :=
    ((((groupedAlternatingHarmonic_hasSum24.mul_left (-40)).add
        (odd_square_hasSum24.mul_left 8)).add
      (shifted_zeta_two_hasSum.mul_left 8)).add
      (shifted_zeta_three_hasSum24.mul_left 3)).add
      (shifted_zeta_four_hasSum24.mul_left (5 / 2))
  convert hCombined using 1
  · funext m
    unfold rationalCorrectionTerm24
    dsimp only
    have hm1 : (m : ℝ) + 1 ≠ 0 := by positivity
    have hodd : 2 * (m : ℝ) + 1 ≠ 0 := by positivity
    have hpoly :
        1 + (m : ℝ) * 4 + (m : ℝ) ^ 2 * 4 ≠ 0 := by positivity
    field_simp [hm1, hodd, hpoly]
    ring_nf
    field_simp [hpoly]
    ring
  · unfold rationalCorrectionValue24
    ring

/-- Only the ordinary-harmonic half of the paired weight-two Euler sum is
an external input: the signed-harmonic half was evaluated above by the
finite triangular identity. -/
theorem pairedAlternatingLinear24_of_harmonic
    (hHarmonic :
      HasSum pairedAlternatingHarmonicEulerTerm24
        pairedAlternatingHarmonicEulerValue24) :
    HasSum pairedAlternatingLinearEulerTerm24
      alternatingLinearEulerValue24 := by
  have hCombined :=
    hHarmonic.add
      (pairedAlternatingSignedHarmonic_hasSum24.mul_left 2)
  convert hCombined using 1
  · funext m
    rw [pairedAlternatingLinearEulerTerm24_decomposition]
  · unfold pairedAlternatingHarmonicEulerValue24
      alternatingLinearEulerValue24
    ring

/-- Reduce the shifted linear input to three standard lower-weight Euler
sums.  The partial fraction and parity bookkeeping are formalized in
`shiftedLinearEulerTerm24_lower_weight_decomposition`. -/
theorem shiftedLinearEuler24_of_lower_weight
    (hPairedAlternatingLinear :
      HasSum pairedAlternatingLinearEulerTerm24
        alternatingLinearEulerValue24)
    (hQuadraticLinear :
      HasSum quadraticLinearEulerTerm24 quadraticLinearEulerValue24)
    (hAlternatingQuadraticLinear :
      HasSum alternatingQuadraticLinearEulerTerm24
        alternatingQuadraticLinearEulerValue24) :
    HasSum shiftedLinearEulerTerm24 shiftedLinearEulerValue24 := by
  have hQuadraticEven :
      HasSum (fun m : ℕ ↦ quadraticLinearEulerTerm24 (2 * m + 1))
        ((quadraticLinearEulerValue24 +
          alternatingQuadraticLinearEulerValue24) / 2) := by
    apply hasSum_even_position24 hQuadraticLinear
    simpa [alternatingQuadraticLinearEulerTerm24] using
      hAlternatingQuadraticLinear
  have hCombined :=
    ((hPairedAlternatingLinear.mul_left (-1)).add
      (groupedAlternatingHarmonic_hasSum24.mul_left 3)).sub
      hQuadraticEven
  convert hCombined using 1
  · funext m
    rw [shiftedLinearEulerTerm24_lower_weight_decomposition]
    unfold pairedAlternatingLinearEulerTerm24
    ring
  · unfold alternatingLinearEulerValue24 quadraticLinearEulerValue24
      alternatingQuadraticLinearEulerValue24 shiftedLinearEulerValue24
    ring

/-- A fully explicit Euler-sum certificate for the only problem-specific
analytic scalar left after the finite binomial reduction.  The five
hypotheses are standard level-two Euler sums of weight at most four; the
rational tail is evaluated unconditionally above.  The proof below verifies
the even-index shift and all coefficient algebra. -/
theorem elementarySeries24_of_euler_certificate
    (hQuadratic :
      HasSum quadraticEulerTerm24 quadraticEulerValue24)
    (hAlternatingQuadratic :
      HasSum alternatingQuadraticEulerTerm24
        alternatingQuadraticEulerValue24)
    (hCubic :
      HasSum cubicLinearEulerTerm24 cubicLinearEulerValue24)
    (hAlternatingCubic :
      HasSum alternatingCubicLinearEulerTerm24
        alternatingCubicLinearEulerValue24)
    (hShifted :
      HasSum shiftedLinearEulerTerm24 shiftedLinearEulerValue24) :
    HasSum elementaryOuterTerm24 elementarySeriesValue24 := by
  have hQuadraticEven :
      HasSum (fun m : ℕ ↦ quadraticEulerTerm24 (2 * m + 1))
        ((quadraticEulerValue24 + alternatingQuadraticEulerValue24) / 2) := by
    apply hasSum_even_position24 hQuadratic
    simpa [alternatingQuadraticEulerTerm24] using hAlternatingQuadratic
  have hCubicEven :
      HasSum (fun m : ℕ ↦ cubicLinearEulerTerm24 (2 * m + 1))
        ((cubicLinearEulerValue24 +
          alternatingCubicLinearEulerValue24) / 2) := by
    apply hasSum_even_position24 hCubic
    simpa [alternatingCubicLinearEulerTerm24] using hAlternatingCubic
  have hCombined :=
    (((hQuadraticEven.mul_left 4).add (hShifted.mul_left 8)).sub
      (hCubicEven.mul_left 24)).add rationalCorrection_hasSum24
  convert hCombined using 1
  · funext m
    rw [elementaryOuterTerm24_euler_decomposition]
  · unfold quadraticEulerValue24 alternatingQuadraticEulerValue24
      cubicLinearEulerValue24 alternatingCubicLinearEulerValue24
      shiftedLinearEulerValue24 rationalCorrectionValue24
      elementarySeriesValue24
    ring

/-- The even-index quadratic value.  By `P(2j) = 2H_j - H_{2j}` the underlying
series is `(1/4) ∑_j (r_j² - H_{2j}^{(2)})/j²` with `r_j = 2H_j - H_{2j}`, i.e.
exactly the object Stage 1 produces — no odd-index `parityRemainder24` occurs. -/
noncomputable def evenQuadraticEulerValue24 : ℝ :=
  (quadraticEulerValue24 + alternatingQuadraticEulerValue24) / 2

/-- Variant of `elementarySeries24_of_euler_certificate` taking the EVEN-INDEX
quadratic sum as a single hypothesis.

The original theorem takes `hQuadratic` and `hAlternatingQuadratic` separately,
but uses them only to build `hQuadraticEven` via `hasSum_even_position24`; the
rest of its proof never mentions either again.  So one evaluation suffices where
the certificate previously asked for two. -/
theorem elementarySeries24_of_even_quadratic
    (hQuadraticEven :
      HasSum (fun m : ℕ ↦ quadraticEulerTerm24 (2 * m + 1))
        evenQuadraticEulerValue24)
    (hCubic :
      HasSum cubicLinearEulerTerm24 cubicLinearEulerValue24)
    (hAlternatingCubic :
      HasSum alternatingCubicLinearEulerTerm24
        alternatingCubicLinearEulerValue24)
    (hShifted :
      HasSum shiftedLinearEulerTerm24 shiftedLinearEulerValue24) :
    HasSum elementaryOuterTerm24 elementarySeriesValue24 := by
  have hQE : HasSum (fun m : ℕ ↦ quadraticEulerTerm24 (2 * m + 1))
      ((quadraticEulerValue24 + alternatingQuadraticEulerValue24) / 2) := by
    simpa [evenQuadraticEulerValue24] using hQuadraticEven
  have hCubicEven :
      HasSum (fun m : ℕ ↦ cubicLinearEulerTerm24 (2 * m + 1))
        ((cubicLinearEulerValue24 +
          alternatingCubicLinearEulerValue24) / 2) := by
    apply hasSum_even_position24 hCubic
    simpa [alternatingCubicLinearEulerTerm24] using hAlternatingCubic
  have hCombined :=
    (((hQE.mul_left 4).add (hShifted.mul_left 8)).sub
      (hCubicEven.mul_left 24)).add rationalCorrection_hasSum24
  convert hCombined using 1
  · funext m
    rw [elementaryOuterTerm24_euler_decomposition]
  · unfold quadraticEulerValue24 alternatingQuadraticEulerValue24
      cubicLinearEulerValue24 alternatingCubicLinearEulerValue24
      shiftedLinearEulerValue24 rationalCorrectionValue24
      elementarySeriesValue24
    ring

/-- Expanded certificate in which every remaining input is a standard
level-two Euler sum, without the problem-specific shifted summand. -/
theorem elementarySeries24_of_standard_euler_certificate
    (hQuadratic :
      HasSum quadraticEulerTerm24 quadraticEulerValue24)
    (hAlternatingQuadratic :
      HasSum alternatingQuadraticEulerTerm24
        alternatingQuadraticEulerValue24)
    (hCubic :
      HasSum cubicLinearEulerTerm24 cubicLinearEulerValue24)
    (hAlternatingCubic :
      HasSum alternatingCubicLinearEulerTerm24
        alternatingCubicLinearEulerValue24)
    (hPairedAlternatingHarmonic :
      HasSum pairedAlternatingHarmonicEulerTerm24
        pairedAlternatingHarmonicEulerValue24)
    (hQuadraticLinear :
      HasSum quadraticLinearEulerTerm24 quadraticLinearEulerValue24)
    (hAlternatingQuadraticLinear :
      HasSum alternatingQuadraticLinearEulerTerm24
        alternatingQuadraticLinearEulerValue24) :
    HasSum elementaryOuterTerm24 elementarySeriesValue24 :=
  elementarySeries24_of_euler_certificate hQuadratic
    hAlternatingQuadratic hCubic hAlternatingCubic
    (shiftedLinearEuler24_of_lower_weight
      (pairedAlternatingLinear24_of_harmonic
        hPairedAlternatingHarmonic)
      hQuadraticLinear hAlternatingQuadraticLinear)

/-- Solve the two exact derivative-WZ equations for the classical
inverse-central fourth-power sum.

The remaining analytic interfaces here are deliberately the two WZ boundary
statements themselves:

* Leshchiner: the derivative pair has left boundary the alternating
  fourth-power series and right boundary `leshchinerWeightFourTerm24`;
* Bailey--Borwein--Bradley: the derivative pair has left boundary `ζ(4)`
  and right boundary `bbbWeightFourTerm24`.

Their finite rectangle identities are rational-function algebra.  To turn
those into these `HasSum`s one must additionally prove both outer factorial
boundaries tend to zero; no unverified infinite-boundary claim is used below. -/
theorem inverseCentralFourthCoefficient24_of_wz
    (hLeshchiner :
      HasSum leshchinerWeightFourTerm24
        ((7 / 8) * (Real.pi ^ 4 / 90)))
    (hBBB :
      HasSum bbbWeightFourTerm24 (Real.pi ^ 4 / 90)) :
    HasSum inverseCentralFourthCoefficient24
      ((17 / 36) * (Real.pi ^ 4 / 90)) := by
  have hCombined :=
    (hLeshchiner.mul_left (2 / 3)).sub
      (hBBB.mul_left (1 / 9))
  convert hCombined using 1
  · funext j
    unfold leshchinerWeightFourTerm24 bbbWeightFourTerm24
    ring
  · ring

/-- The inverse-central scalar evaluation follows formally from three standard
Apéry/BBB series: the `ζ(2)` inverse-binomial sum, the classical `S(4)`
evaluation, and the weight-four coefficient of the
Bailey--Borwein--Bradley generating identity.  These analytic results remain
explicit hypotheses; this theorem verifies the indexing and coefficient
algebra that turns them into the exact series occurring in Problem 2.4. -/
theorem inverseCentralSeries24_of_classical
    (hS2 :
      HasSum inverseCentralCoefficient24 (Real.pi ^ 2 / 6 / 3))
    (hS4 :
      HasSum inverseCentralFourthCoefficient24
        ((17/36) * (Real.pi ^ 4 / 90)))
    (hBBB :
      HasSum bbbWeightFourTerm24 (Real.pi ^ 4 / 90)) :
    HasSum
      (fun j : ℕ ↦ inverseCentralCoefficient24 j *
        (Real.pi ^ 2 / 6 - harmonicSquare24 (j + 1)))
      inverseCentralSeriesValue24 := by
  have hNinePrevious :
      HasSum
        (fun j : ℕ ↦
          9 * (harmonicSquare24 j * inverseCentralCoefficient24 j))
        (3 * ((17/36) * (Real.pi ^ 4 / 90)) -
          Real.pi ^ 4 / 90) := by
    convert (hS4.mul_left 3).sub hBBB using 1
    · funext j
      unfold bbbWeightFourTerm24
      ring
  have hPrevious :
      HasSum
        (fun j : ℕ ↦
          harmonicSquare24 j * inverseCentralCoefficient24 j)
        ((1/9) * (3 * ((17/36) * (Real.pi ^ 4 / 90)) -
          Real.pi ^ 4 / 90)) := by
    convert hNinePrevious.mul_left (1/9) using 1
    · funext j
      ring
  have hFull :
      HasSum
        (fun j : ℕ ↦
          harmonicSquare24 (j + 1) * inverseCentralCoefficient24 j)
        ((1/9) * (3 * ((17/36) * (Real.pi ^ 4 / 90)) -
            Real.pi ^ 4 / 90) +
          (17/36) * (Real.pi ^ 4 / 90)) := by
    convert hPrevious.add hS4 using 1
    · funext j
      rw [harmonicSquare24_succ]
      unfold inverseCentralFourthCoefficient24
      ring
  convert (hS2.mul_left (Real.pi ^ 2 / 6)).sub hFull using 1
  · funext j
    ring
  · unfold inverseCentralSeriesValue24
    ring

/-- After importing the unconditional weight-two evaluation from Problem 2.6,
only the two weight-four inverse-central identities remain as inputs. -/
theorem inverseCentralSeries24_of_classical_weightFour
    (hS4 :
      HasSum inverseCentralFourthCoefficient24
        ((17 / 36) * (Real.pi ^ 4 / 90)))
    (hBBB :
      HasSum bbbWeightFourTerm24 (Real.pi ^ 4 / 90)) :
    HasSum
      (fun j : ℕ ↦ inverseCentralCoefficient24 j *
        (Real.pi ^ 2 / 6 - harmonicSquare24 (j + 1)))
      inverseCentralSeriesValue24 :=
  inverseCentralSeries24_of_classical
    inverseCentralCoefficient24_hasSum hS4 hBBB

/-- Exact inverse-central interface using the two derivative WZ identities,
with the weight-two sum already discharged unconditionally. -/
theorem inverseCentralSeries24_of_wz
    (hLeshchiner :
      HasSum leshchinerWeightFourTerm24
        ((7 / 8) * (Real.pi ^ 4 / 90)))
    (hBBB :
      HasSum bbbWeightFourTerm24 (Real.pi ^ 4 / 90)) :
    HasSum
      (fun j : ℕ ↦ inverseCentralCoefficient24 j *
        (Real.pi ^ 2 / 6 - harmonicSquare24 (j + 1)))
      inverseCentralSeriesValue24 :=
  inverseCentralSeries24_of_classical_weightFour
    (inverseCentralFourthCoefficient24_of_wz hLeshchiner hBBB) hBBB

/-- Once the two displayed scalar `HasSum` evaluations are supplied, all
finite identities, convergence, rearrangement, and final coefficient
collection needed for Problem 2.4 are already formalized.  The hypotheses
state the exact remaining analytic inputs rather than hiding them behind an
unrelated existential. -/
theorem problem24_of_weightFour_evaluations
    (hElementary :
      HasSum elementaryOuterTerm24 elementarySeriesValue24)
    (hInverseCentral :
      HasSum
        (fun j : ℕ ↦ inverseCentralCoefficient24 j *
          (Real.pi ^ 2 / 6 - harmonicSquare24 (j + 1)))
        inverseCentralSeriesValue24) :
    problem24Statement := by
  unfold problem24Statement
  rw [lhs_24_eq_elementary_add_inverse_central,
    hElementary.tsum_eq, hInverseCentral.tsum_eq]
  unfold elementarySeriesValue24 inverseCentralSeriesValue24 rhs_24
  ring

/-- End-to-end certificate for Problem 2.4 from the explicit Euler and
inverse-central evaluations isolated above. -/
theorem problem24_of_euler_and_classical
    (hQuadratic :
      HasSum quadraticEulerTerm24 quadraticEulerValue24)
    (hAlternatingQuadratic :
      HasSum alternatingQuadraticEulerTerm24
        alternatingQuadraticEulerValue24)
    (hCubic :
      HasSum cubicLinearEulerTerm24 cubicLinearEulerValue24)
    (hAlternatingCubic :
      HasSum alternatingCubicLinearEulerTerm24
        alternatingCubicLinearEulerValue24)
    (hShifted :
      HasSum shiftedLinearEulerTerm24 shiftedLinearEulerValue24)
    (hLeshchiner :
      HasSum leshchinerWeightFourTerm24
        ((7 / 8) * (Real.pi ^ 4 / 90)))
    (hBBB :
      HasSum bbbWeightFourTerm24 (Real.pi ^ 4 / 90)) :
    problem24Statement :=
  problem24_of_weightFour_evaluations
    (elementarySeries24_of_euler_certificate hQuadratic
      hAlternatingQuadratic hCubic hAlternatingCubic hShifted)
    (inverseCentralSeries24_of_wz hLeshchiner hBBB)

/-- End-to-end certificate whose nine hypotheses are exactly the remaining
standard special-value evaluations: seven level-two Euler sums and the two
weight-four derivative-WZ boundary sums. -/
theorem problem24_of_standard_euler_and_wz
    (hQuadratic :
      HasSum quadraticEulerTerm24 quadraticEulerValue24)
    (hAlternatingQuadratic :
      HasSum alternatingQuadraticEulerTerm24
        alternatingQuadraticEulerValue24)
    (hCubic :
      HasSum cubicLinearEulerTerm24 cubicLinearEulerValue24)
    (hAlternatingCubic :
      HasSum alternatingCubicLinearEulerTerm24
        alternatingCubicLinearEulerValue24)
    (hPairedAlternatingHarmonic :
      HasSum pairedAlternatingHarmonicEulerTerm24
        pairedAlternatingHarmonicEulerValue24)
    (hQuadraticLinear :
      HasSum quadraticLinearEulerTerm24 quadraticLinearEulerValue24)
    (hAlternatingQuadraticLinear :
      HasSum alternatingQuadraticLinearEulerTerm24
        alternatingQuadraticLinearEulerValue24)
    (hLeshchiner :
      HasSum leshchinerWeightFourTerm24
        ((7 / 8) * (Real.pi ^ 4 / 90)))
    (hBBB :
      HasSum bbbWeightFourTerm24 (Real.pi ^ 4 / 90)) :
    problem24Statement :=
  problem24_of_weightFour_evaluations
    (elementarySeries24_of_standard_euler_certificate
      hQuadratic hAlternatingQuadratic hCubic hAlternatingCubic
      hPairedAlternatingHarmonic hQuadraticLinear
      hAlternatingQuadraticLinear)
    (inverseCentralSeries24_of_wz hLeshchiner hBBB)

/-!
The exact proposition `problem24Statement` now follows from explicit standard
Euler-sum and inverse-central certificates.  All finite binomial identities,
convergence and rearrangement, the parity reduction, the grouped alternating
harmonic evaluation, and the rational tail are proved unconditionally above.
The remaining hypotheses in `problem24_of_euler_and_classical` are the
classical special-value evaluations not presently available in Mathlib; none
is hidden behind an unproved constant or an unrelated existential.
-/

end
