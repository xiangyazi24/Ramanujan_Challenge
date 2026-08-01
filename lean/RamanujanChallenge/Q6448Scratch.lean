import RamanujanChallenge.Problem27BarnesTelescoper
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Complex

open Filter Set MeasureTheory Topology
open scoped BigOperators Interval Real

noncomputable section

namespace RamanujanChallenge.P27.Q6448

def verticalPoint (x y : ℝ) : ℂ :=
  (x : ℂ) + (y : ℂ) * Complex.I

def rawKernel27 (t : ℂ) : ℂ :=
  zudilinBarnesSquaredSineKernel27 t

def rPhiVertical27 (n : ℕ) (x : ℝ) : ℂ :=
  ∫ y : ℝ,
    ctRPhi27 n (verticalPoint x y) * rawKernel27 (verticalPoint x y)

def sPhiVertical27 (n : ℕ) (x : ℝ) : ℂ :=
  ∫ y : ℝ,
    ctSPhi27 n (verticalPoint x y) * rawKernel27 (verticalPoint x y)

theorem rawKernel_add_one27 (t : ℂ) :
    rawKernel27 (t + 1) = rawKernel27 t := by
  unfold rawKernel27 zudilinBarnesSquaredSineKernel27
  rw [show (Real.pi : ℂ) * (t + 1) =
      (Real.pi : ℂ) * t + (Real.pi : ℂ) by ring,
    Complex.sin_add]
  simp only [Complex.sin_pi, Complex.cos_pi, mul_zero, mul_neg_one,
    zero_add]
  rw [div_neg]
  ring

theorem rawKernel_add_nat27 (r : ℕ) (t : ℂ) :
    rawKernel27 (t + (r : ℂ)) = rawKernel27 t := by
  induction r with
  | zero => simp
  | succ r ih =>
      calc
        rawKernel27 (t + ((r + 1 : ℕ) : ℂ)) =
            rawKernel27 ((t + (r : ℂ)) + 1) := by
          congr 1
          push_cast
          ring
        _ = rawKernel27 (t + (r : ℂ)) := rawKernel_add_one27 _
        _ = rawKernel27 t := ih

theorem nativePoint_eq_translate27 (n : ℕ) (y : ℝ) :
    verticalPoint ((n : ℝ) + 1 / 2) y =
      zudilinBarnesLine27 y + (((n + 1 : ℕ) : ℂ)) := by
  unfold verticalPoint zudilinBarnesLine27
  push_cast
  ring

theorem raw_fixedLine_eq_native27 (n : ℕ) :
    (∫ y : ℝ,
      zudilinBarnesPhi27 n (zudilinBarnesLine27 y) *
        rawKernel27 (zudilinBarnesLine27 y)) =
      rPhiVertical27 n ((n : ℝ) + 1 / 2) := by
  unfold rPhiVertical27
  apply integral_congr_ae
  filter_upwards with y
  rw [zudilinBarnesPhi_eq_ctRPhi_translate27]
  have hp := nativePoint_eq_translate27 n y
  rw [hp, rawKernel_add_nat27]

theorem shift_down_of_one_strip27
    (V : ℕ → ℝ → ℂ)
    (hone : ∀ {n m : ℕ}, 1 ≤ m → m ≤ n →
      V n ((m : ℝ) - 1 / 2) = V n ((m : ℝ) + 1 / 2))
    (n m : ℕ) (hm : m ≤ n) :
    V n ((m : ℝ) + 1 / 2) = V n (1 / 2) := by
  induction m with
  | zero => norm_num
  | succ m ih =>
      have hm' : m ≤ n := by omega
      have hstrip := hone (n := n) (m := m + 1) (by omega) hm
      calc
        V n (((m + 1 : ℕ) : ℝ) + 1 / 2) =
            V n ((m : ℝ) + 1 / 2) := by
          convert hstrip.symm using 1 <;> push_cast <;> ring
        _ = V n (1 / 2) := ih hm'

theorem native_to_half_of_one_strip27
    (V : ℕ → ℝ → ℂ)
    (hone : ∀ {n m : ℕ}, 1 ≤ m → m ≤ n →
      V n ((m : ℝ) - 1 / 2) = V n ((m : ℝ) + 1 / 2))
    (n : ℕ) :
    V n ((n : ℝ) + 1 / 2) = V n (1 / 2) :=
  shift_down_of_one_strip27 V hone n n le_rfl

theorem sPhi_shifted_integral_eq_of_one_strip27
    (hone : ∀ k : ℕ,
      sPhiVertical27 (k + 2) (1 / 2) =
        sPhiVertical27 (k + 2) (3 / 2))
    (k : ℕ) :
    (∫ y : ℝ,
      ctSPhi27 (k + 2) (verticalPoint (1 / 2) y + 1) *
        rawKernel27 (verticalPoint (1 / 2) y)) =
      ∫ y : ℝ,
        ctSPhi27 (k + 2) (verticalPoint (1 / 2) y) *
          rawKernel27 (verticalPoint (1 / 2) y) := by
  calc
    (∫ y : ℝ,
      ctSPhi27 (k + 2) (verticalPoint (1 / 2) y + 1) *
        rawKernel27 (verticalPoint (1 / 2) y)) =
        sPhiVertical27 (k + 2) (3 / 2) := by
      unfold sPhiVertical27
      apply integral_congr_ae
      filter_upwards with y
      have hp : verticalPoint (1 / 2) y + 1 =
          verticalPoint (3 / 2) y := by
        unfold verticalPoint
        push_cast
        ring
      rw [← rawKernel_add_one27 (verticalPoint (1 / 2) y), hp]
    _ = sPhiVertical27 (k + 2) (1 / 2) := (hone k).symm
    _ = ∫ y : ℝ,
        ctSPhi27 (k + 2) (verticalPoint (1 / 2) y) *
          rawKernel27 (verticalPoint (1 / 2) y) := rfl

#print axioms rawKernel_add_nat27
#print axioms raw_fixedLine_eq_native27
#print axioms native_to_half_of_one_strip27
#print axioms sPhi_shifted_integral_eq_of_one_strip27

end RamanujanChallenge.P27.Q6448
