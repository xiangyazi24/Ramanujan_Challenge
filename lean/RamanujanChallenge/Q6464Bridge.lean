import RamanujanChallenge.Problem27BarnesTelescoper
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Complex

open Filter Set MeasureTheory Topology
open scoped BigOperators Interval Real

noncomputable section

namespace RamanujanChallenge.P27.Q6464

/-- Local vertical-line parameterization.  In the consolidated contour module,
use the existing `verticalPoint27` instead. -/
def verticalPoint (x y : ℝ) : ℂ :=
  (x : ℂ) + (y : ℂ) * Complex.I

/-- The raw squared-sine Barnes kernel. -/
def rawKernel (z : ℂ) : ℂ :=
  zudilinBarnesSquaredSineKernel27 z

/-- Raw `ctRPhi27` integral on the vertical line `Re t = x`. -/
def rVertical (n : ℕ) (x : ℝ) : ℂ :=
  ∫ y : ℝ,
    ctRPhi27 n (verticalPoint x y) * rawKernel (verticalPoint x y)

/-- Raw certificate integral on the vertical line `Re t = x`. -/
def sVertical (n : ℕ) (x : ℝ) : ℂ :=
  ∫ y : ℝ,
    ctSPhi27 n (verticalPoint x y) * rawKernel (verticalPoint x y)

/-- Sine is antiperiodic, hence the squared quotient is period one. -/
theorem rawKernel_add_one (z : ℂ) : rawKernel (z + 1) = rawKernel z := by
  unfold rawKernel zudilinBarnesSquaredSineKernel27
  rw [show (Real.pi : ℂ) * (z + 1) =
      (Real.pi : ℂ) * z + (Real.pi : ℂ) by ring,
    Complex.sin_add]
  simp only [Complex.sin_pi, Complex.cos_pi, mul_zero, mul_neg_one,
    add_zero]
  rw [div_neg]
  ring

/-- Integer periodicity of the raw squared-sine kernel. -/
theorem rawKernel_add_nat (r : ℕ) (z : ℂ) :
    rawKernel (z + (r : ℂ)) = rawKernel z := by
  induction r with
  | zero => simp
  | succ r ih =>
      calc
        rawKernel (z + ((r + 1 : ℕ) : ℂ)) =
            rawKernel ((z + (r : ℂ)) + 1) := by
          congr 1
          push_cast
          ring
        _ = rawKernel (z + (r : ℂ)) := rawKernel_add_one _
        _ = rawKernel z := ih

/-- Translation of the repository fixed Barnes line to the native
`t`-coordinate line `Re t = n + 1/2`. -/
theorem nativePoint_eq_translate (n : ℕ) (y : ℝ) :
    verticalPoint ((n : ℝ) + 1 / 2) y =
      zudilinBarnesLine27 y + (((n + 1 : ℕ) : ℂ)) := by
  unfold verticalPoint zudilinBarnesLine27
  push_cast
  ring

/-- Equality of the repository's raw fixed-line integral and the native
`ctRPhi27` vertical integral.  This uses only pointwise rewriting under a
single integral, so it needs no `Integrable` hypothesis. -/
theorem fixedLineRaw_eq_native (n : ℕ) :
    (∫ y : ℝ,
      zudilinBarnesPhi27 n (zudilinBarnesLine27 y) *
        rawKernel (zudilinBarnesLine27 y)) =
      rVertical n ((n : ℝ) + 1 / 2) := by
  unfold rVertical
  apply integral_congr_ae
  filter_upwards with y
  rw [zudilinBarnesPhi_eq_ctRPhi_translate27]
  have hp := nativePoint_eq_translate n y
  rw [hp, rawKernel_add_nat]

/-- The exact normalization connecting the native raw vertical integral to
`zudilinBarnesErrorIntegral27`. -/
theorem errorIntegral_eq_native (n : ℕ) :
    zudilinBarnesErrorIntegral27 n =
      (1 / (2 * (Real.pi : ℂ))) *
        rVertical n ((n : ℝ) + 1 / 2) := by
  rw [zudilinBarnesErrorIntegral_eq_fixedLine27]
  congr 1
  exact fixedLineRaw_eq_native n

/-- Downward induction from the right edge of strip `m` to the common
`1/2` line.  The one-strip equality is used in reverse at index `m+1`. -/
theorem rVertical_shift_down_of_one_strip
    (hone : ∀ {n m : ℕ}, 1 ≤ m → m ≤ n →
      rVertical n ((m : ℝ) - 1 / 2) =
        rVertical n ((m : ℝ) + 1 / 2))
    (n m : ℕ) (hmn : m ≤ n) :
    rVertical n ((m : ℝ) + 1 / 2) = rVertical n (1 / 2) := by
  induction m with
  | zero => norm_num
  | succ m ih =>
      have hm' : m ≤ n := by omega
      have hstrip := hone (n := n) (m := m + 1) (by omega) hmn
      calc
        rVertical n (((m + 1 : ℕ) : ℝ) + 1 / 2) =
            rVertical n ((m : ℝ) + 1 / 2) := by
          convert hstrip.symm using 1 <;> push_cast <;> ring_nf
        _ = rVertical n (1 / 2) := ih hm'

/-- Native line `n+1/2` shifted to the common line `1/2`. -/
theorem rVertical_native_to_half_of_one_strip
    (hone : ∀ {n m : ℕ}, 1 ≤ m → m ≤ n →
      rVertical n ((m : ℝ) - 1 / 2) =
        rVertical n ((m : ℝ) + 1 / 2))
    (n : ℕ) :
    rVertical n ((n : ℝ) + 1 / 2) = rVertical n (1 / 2) :=
  rVertical_shift_down_of_one_strip hone n n le_rfl

/-- Final normalization after all `ctRPhi27` strips have been crossed. -/
theorem errorIntegral_eq_half_of_one_strip
    (hone : ∀ {n m : ℕ}, 1 ≤ m → m ≤ n →
      rVertical n ((m : ℝ) - 1 / 2) =
        rVertical n ((m : ℝ) + 1 / 2))
    (n : ℕ) :
    zudilinBarnesErrorIntegral27 n =
      (1 / (2 * (Real.pi : ℂ))) * rVertical n (1 / 2) := by
  rw [errorIntegral_eq_native,
    rVertical_native_to_half_of_one_strip hone]

/-- Rewrites the shifted certificate term on the `1/2` line as the ordinary
raw certificate integral on the `3/2` line. -/
theorem shiftedCertificate_eq_threeHalf (k : ℕ) :
    (∫ y : ℝ,
      ctSPhi27 (k + 2) (verticalPoint (1 / 2) y + 1) *
        rawKernel (verticalPoint (1 / 2) y)) =
      sVertical (k + 2) (3 / 2) := by
  unfold sVertical
  apply integral_congr_ae
  filter_upwards with y
  have hp : verticalPoint (1 / 2) y + 1 =
      verticalPoint (3 / 2) y := by
    unfold verticalPoint
    push_cast
    ring
  rw [← rawKernel_add_one (verticalPoint (1 / 2) y), hp]

/-- The exact certificate cancellation.  The supplied one-strip theorem is
left=`1/2`, right=`3/2`; the telescoper needs the reverse orientation. -/
theorem shiftedCertificate_eq_unshifted_of_one_strip
    (hone : ∀ k : ℕ,
      sVertical (k + 2) (1 / 2) = sVertical (k + 2) (3 / 2))
    (k : ℕ) :
    (∫ y : ℝ,
      ctSPhi27 (k + 2) (verticalPoint (1 / 2) y + 1) *
        rawKernel (verticalPoint (1 / 2) y)) =
      ∫ y : ℝ,
        ctSPhi27 (k + 2) (verticalPoint (1 / 2) y) *
          rawKernel (verticalPoint (1 / 2) y) := by
  calc
    _ = sVertical (k + 2) (3 / 2) := shiftedCertificate_eq_threeHalf k
    _ = sVertical (k + 2) (1 / 2) := (hone k).symm
    _ = _ := rfl

#print axioms rawKernel_add_nat
#print axioms fixedLineRaw_eq_native
#print axioms errorIntegral_eq_native
#print axioms rVertical_shift_down_of_one_strip
#print axioms errorIntegral_eq_half_of_one_strip
#print axioms shiftedCertificate_eq_unshifted_of_one_strip

end RamanujanChallenge.P27.Q6464
