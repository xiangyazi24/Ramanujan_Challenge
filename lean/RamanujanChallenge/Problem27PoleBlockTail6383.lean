import RamanujanChallenge.Problem27PoleBlockSeries6383
import Mathlib.NumberTheory.ZetaValues
import Mathlib.Topology.Algebra.InfiniteSum.NatInt

open Filter Set MeasureTheory Topology
open scoped Topology BigOperators

noncomputable section

namespace RamanujanChallenge.P27.Q6383

def tailCorrection (m : ℕ) : ℝ :=
  ∑ j ∈ Finset.range m,
    (1 / (j + 1 : ℝ) ^ 2 + 1 / (j + 1 : ℝ) ^ 3)

theorem hasSum_zeta_two_succ :
    HasSum (fun n : ℕ => 1 / (n + 1 : ℝ) ^ 2) (Real.pi ^ 2 / 6) := by
  have h := (hasSum_nat_add_iff'
    (f := fun n : ℕ => 1 / (n : ℝ) ^ 2)
    (g := Real.pi ^ 2 / 6) 1).2 hasSum_zeta_two
  simpa using h

theorem summable_zeta_three_succ :
    Summable (fun n : ℕ => 1 / (n + 1 : ℝ) ^ 3) := by
  exact (summable_nat_add_iff 1).2
    (summable_one_div_nat_pow.mpr (by norm_num))

theorem tsum_zeta_three_succ :
    (∑' n : ℕ, 1 / (n + 1 : ℝ) ^ 3) = zeta3 := rfl

theorem tsum_shifted_two (m : ℕ) :
    (∑' k : ℕ, 1 / (m + k + 1 : ℝ) ^ 2) =
      Real.pi ^ 2 / 6 -
        ∑ j ∈ Finset.range m, 1 / (j + 1 : ℝ) ^ 2 := by
  let f : ℕ → ℝ := fun n => 1 / (n + 1 : ℝ) ^ 2
  have hf : Summable f := hasSum_zeta_two_succ.summable
  have hsplit := hf.sum_add_tsum_nat_add m
  have htotal : (∑' n, f n) = Real.pi ^ 2 / 6 :=
    hasSum_zeta_two_succ.tsum_eq
  rw [htotal] at hsplit
  calc
    (∑' k : ℕ, 1 / (m + k + 1 : ℝ) ^ 2) =
        ∑' k : ℕ, f (k + m) := by
          apply tsum_congr
          intro k
          dsimp only [f]
          congr 3
          push_cast
          omega
    _ = Real.pi ^ 2 / 6 - ∑ j ∈ Finset.range m, f j := by
          linarith
    _ = _ := by rfl

theorem tsum_shifted_three (m : ℕ) :
    (∑' k : ℕ, 1 / (m + k + 1 : ℝ) ^ 3) =
      zeta3 - ∑ j ∈ Finset.range m, 1 / (j + 1 : ℝ) ^ 3 := by
  let f : ℕ → ℝ := fun n => 1 / (n + 1 : ℝ) ^ 3
  have hf : Summable f := summable_zeta_three_succ
  have hsplit := hf.sum_add_tsum_nat_add m
  have htotal : (∑' n, f n) = zeta3 := rfl
  rw [htotal] at hsplit
  calc
    (∑' k : ℕ, 1 / (m + k + 1 : ℝ) ^ 3) =
        ∑' k : ℕ, f (k + m) := by
          apply tsum_congr
          intro k
          dsimp only [f]
          congr 3
          push_cast
          omega
    _ = zeta3 - ∑ j ∈ Finset.range m, f j := by
          linarith
    _ = _ := by rfl

theorem tsum_shifted_two_three (m : ℕ) :
    (∑' k : ℕ,
      (1 / (m + k + 1 : ℝ) ^ 2 +
        1 / (m + k + 1 : ℝ) ^ 3)) =
      Real.pi ^ 2 / 6 + zeta3 - tailCorrection m := by
  have h2 : Summable (fun k : ℕ => 1 / (m + k + 1 : ℝ) ^ 2) := by
    simpa [Nat.add_assoc, Nat.add_comm, Nat.add_left_comm] using
      ((summable_nat_add_iff (m + 1)).2
        (summable_one_div_nat_pow.mpr (by norm_num) :
          Summable (fun n : ℕ => 1 / (n : ℝ) ^ 2)))
  have h3 : Summable (fun k : ℕ => 1 / (m + k + 1 : ℝ) ^ 3) := by
    simpa [Nat.add_assoc, Nat.add_comm, Nat.add_left_comm] using
      ((summable_nat_add_iff (m + 1)).2
        (summable_one_div_nat_pow.mpr (by norm_num) :
          Summable (fun n : ℕ => 1 / (n : ℝ) ^ 3)))
  rw [tsum_add h2 h3, tsum_shifted_two, tsum_shifted_three]
  unfold tailCorrection
  rw [Finset.sum_add_distrib]
  ring

theorem universal_poleBlock_integral (m : ℕ) :
    ((Real.pi / 2 : ℝ) : ℂ) *
        (∫ y : ℝ,
          poleBlock m y / (Real.cosh (Real.pi * y) : ℂ) ^ 2) =
      ((Real.pi ^ 2 / 6 + zeta3 - tailCorrection m : ℝ) : ℂ) := by
  rw [normalized_poleBlock_eq_boseIntegral,
    integral_boseKernel_eq_tsum]
  norm_cast
  exact tsum_shifted_two_three m

theorem universal_poleBlock_integral_zero :
    ((Real.pi / 2 : ℝ) : ℂ) *
        (∫ y : ℝ,
          poleBlock 0 y / (Real.cosh (Real.pi * y) : ℂ) ^ 2) =
      ((Real.pi ^ 2 / 6 + zeta3 : ℝ) : ℂ) := by
  simpa [tailCorrection] using universal_poleBlock_integral 0

theorem universal_poleBlock_integral_one :
    ((Real.pi / 2 : ℝ) : ℂ) *
        (∫ y : ℝ,
          poleBlock 1 y / (Real.cosh (Real.pi * y) : ℂ) ^ 2) =
      ((Real.pi ^ 2 / 6 + zeta3 - 2 : ℝ) : ℂ) := by
  convert universal_poleBlock_integral 1 using 1
  norm_num [tailCorrection]

theorem universal_poleBlock_integral_two :
    ((Real.pi / 2 : ℝ) : ℂ) *
        (∫ y : ℝ,
          poleBlock 2 y / (Real.cosh (Real.pi * y) : ℂ) ^ 2) =
      ((Real.pi ^ 2 / 6 + zeta3 - 19 / 8 : ℝ) : ℂ) := by
  convert universal_poleBlock_integral 2 using 1
  norm_num [tailCorrection]

theorem universal_poleBlock_integral_three :
    ((Real.pi / 2 : ℝ) : ℂ) *
        (∫ y : ℝ,
          poleBlock 3 y / (Real.cosh (Real.pi * y) : ℂ) ^ 2) =
      ((Real.pi ^ 2 / 6 + zeta3 - 545 / 216 : ℝ) : ℂ) := by
  convert universal_poleBlock_integral 3 using 1
  norm_num [tailCorrection]

theorem universal_poleBlock_integral_four :
    ((Real.pi / 2 : ℝ) : ℂ) *
        (∫ y : ℝ,
          poleBlock 4 y / (Real.cosh (Real.pi * y) : ℂ) ^ 2) =
      ((Real.pi ^ 2 / 6 + zeta3 - 4495 / 1728 : ℝ) : ℂ) := by
  convert universal_poleBlock_integral 4 using 1
  norm_num [tailCorrection]

#print axioms universal_poleBlock_integral
#print axioms universal_poleBlock_integral_four

end RamanujanChallenge.P27.Q6383
