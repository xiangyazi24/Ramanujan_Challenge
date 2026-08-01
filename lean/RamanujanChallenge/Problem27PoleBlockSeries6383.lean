import RamanujanChallenge.Problem27PoleBlockFubini6383
import Mathlib.Analysis.PSeries
import Mathlib.MeasureTheory.Integral.DominatedConvergence
import Mathlib.Topology.Algebra.InfiniteSum.Ring

open Filter Set MeasureTheory Topology
open scoped Topology

noncomputable section

namespace RamanujanChallenge.P27.Q6383

def boseKernel (m : ℕ) (t : ℝ) : ℝ :=
  (t + t ^ 2 / 2) * Real.exp (-(((m : ℝ) + 1) * t)) /
    (1 - Real.exp (-t))

def boseTerm (m k : ℕ) (t : ℝ) : ℝ :=
  (t + t ^ 2 / 2) * Real.exp (-((m + k + 1 : ℕ) * t))

theorem two_mul_sinh_half (t : ℝ) :
    2 * Real.sinh (t / 2) =
      Real.exp (t / 2) * (1 - Real.exp (-t)) := by
  rw [Real.sinh_eq]
  field_simp
  rw [mul_sub, mul_one, ← Real.exp_add]
  congr 1
  ring

theorem normalized_fourier_factor_eq_bose
    (m : ℕ) {t : ℝ} (ht : 0 < t) :
    (Real.pi / 2) *
        (laplaceWeight m t *
          (t / (Real.pi * Real.sinh (t / 2)))) =
      boseKernel m t := by
  have hs : Real.sinh (t / 2) ≠ 0 := by
    rw [Real.sinh_ne_zero]
    linarith
  have hd : 1 - Real.exp (-t) ≠ 0 := by
    have he : Real.exp (-t) < 1 := (Real.exp_lt_one_iff).2 (by linarith)
    linarith
  have hexp :
      Real.exp (-(poleAbscissa m * t)) / Real.exp (t / 2) =
        Real.exp (-(((m : ℝ) + 1) * t)) := by
    rw [div_eq_mul_inv, ← Real.exp_neg, ← Real.exp_add]
    congr 1
    unfold poleAbscissa
    ring
  calc
    (Real.pi / 2) *
        (laplaceWeight m t *
          (t / (Real.pi * Real.sinh (t / 2)))) =
      (t + t ^ 2 / 2) * Real.exp (-(poleAbscissa m * t)) /
        (2 * Real.sinh (t / 2)) := by
          unfold laplaceWeight
          field_simp [Real.pi_ne_zero, hs]
          ring
    _ = (t + t ^ 2 / 2) * Real.exp (-(poleAbscissa m * t)) /
        (Real.exp (t / 2) * (1 - Real.exp (-t))) := by
          rw [two_mul_sinh_half]
    _ = (t + t ^ 2 / 2) *
        (Real.exp (-(poleAbscissa m * t)) / Real.exp (t / 2)) /
          (1 - Real.exp (-t)) := by
          field_simp [Real.exp_ne_zero, hd]
    _ = boseKernel m t := by
          rw [hexp]
          rfl

theorem normalized_poleBlock_eq_boseIntegral (m : ℕ) :
    ((Real.pi / 2 : ℝ) : ℂ) *
        (∫ y : ℝ,
          poleBlock m y / (Real.cosh (Real.pi * y) : ℂ) ^ 2) =
      ∫ t : ℝ in Ioi 0, (boseKernel m t : ℂ) := by
  rw [integral_poleBlock_eq_fourierLaplace]
  rw [← MeasureTheory.integral_const_mul]
  apply integral_congr_ae
  filter_upwards [ae_restrict_mem measurableSet_Ioi] with t ht
  exact_mod_cast normalized_fourier_factor_eq_bose m ht

theorem integral_pow_mul_exp_neg_mul
    (n : ℕ) {a : ℝ} (ha : 0 < a) :
    (∫ t : ℝ in Ioi 0, t ^ n * Real.exp (-(a * t))) =
      (n.factorial : ℝ) / a ^ (n + 1) := by
  have h := Real.integral_rpow_mul_exp_neg_mul_Ioi
    (a := (n : ℝ) + 1) (r := a) (by positivity) ha
  convert h using 1
  · apply integral_congr_ae
    filter_upwards [ae_restrict_mem measurableSet_Ioi] with t ht
    rw [show (n : ℝ) + 1 - 1 = n by ring, Real.rpow_natCast]
  · rw [Real.Gamma_nat_eq_factorial]
    rw [one_div, inv_pow]
    field_simp [ha.ne']

theorem integrableOn_boseTerm (m k : ℕ) :
    IntegrableOn (boseTerm m k) (Ioi 0) := by
  let a : ℝ := (m + k + 1 : ℕ)
  have ha : 0 < a := by
    dsimp only [a]
    positivity
  have h1 := integrableOn_pow_mul_exp_neg_mul 1 ha
  have h2 := integrableOn_pow_mul_exp_neg_mul 2 ha
  have h := h1.add (h2.const_mul (1 / 2 : ℝ))
  refine h.congr ?_
  filter_upwards [ae_restrict_mem measurableSet_Ioi] with t ht
  unfold boseTerm
  dsimp only [a] at ha ⊢
  push_cast
  ring

theorem integral_boseTerm (m k : ℕ) :
    (∫ t : ℝ in Ioi 0, boseTerm m k t) =
      1 / (m + k + 1 : ℝ) ^ 2 +
        1 / (m + k + 1 : ℝ) ^ 3 := by
  let a : ℝ := (m + k + 1 : ℕ)
  have ha : 0 < a := by
    dsimp only [a]
    positivity
  have h1 := integral_pow_mul_exp_neg_mul 1 ha
  have h2 := integral_pow_mul_exp_neg_mul 2 ha
  have hi1 : IntegrableOn (fun t : ℝ => t * Real.exp (-(a * t))) (Ioi 0) := by
    simpa using integrableOn_pow_mul_exp_neg_mul 1 ha
  have hi2 : IntegrableOn (fun t : ℝ => t ^ 2 * Real.exp (-(a * t))) (Ioi 0) :=
    integrableOn_pow_mul_exp_neg_mul 2 ha
  calc
    (∫ t : ℝ in Ioi 0, boseTerm m k t) =
        (∫ t : ℝ in Ioi 0, t * Real.exp (-(a * t))) +
          (1 / 2 : ℝ) *
            ∫ t : ℝ in Ioi 0, t ^ 2 * Real.exp (-(a * t)) := by
              rw [← MeasureTheory.integral_const_mul,
                ← integral_add hi1 (hi2.const_mul (1 / 2 : ℝ))]
              apply integral_congr_ae
              filter_upwards [ae_restrict_mem measurableSet_Ioi] with t ht
              unfold boseTerm
              dsimp only [a]
              push_cast
              ring
    _ = 1 / a ^ 2 + 1 / a ^ 3 := by
          rw [h1, h2]
          norm_num [Nat.factorial]
          field_simp [ha.ne']
    _ = _ := by rfl

theorem hasSum_boseTerm {m : ℕ} {t : ℝ} (ht : 0 < t) :
    HasSum (fun k : ℕ => boseTerm m k t) (boseKernel m t) := by
  let q : ℝ := Real.exp (-t)
  let A : ℝ :=
    (t + t ^ 2 / 2) * Real.exp (-(((m : ℝ) + 1) * t))
  have hq : ‖q‖ < 1 := by
    rw [Real.norm_eq_abs, abs_of_pos (Real.exp_pos _)]
    exact (Real.exp_lt_one_iff).2 (by linarith)
  have hgeom := (hasSum_geometric_of_norm_lt_one hq).mul_left A
  have hden : 1 - q ≠ 0 := sub_ne_zero.mpr (ne_of_lt hq)
  convert hgeom using 1
  · funext k
    unfold boseTerm
    dsimp only [A, q]
    rw [← Real.exp_nat_mul]
    rw [← Real.exp_add]
    congr 2
    push_cast
    ring
  · dsimp only [A, q, boseKernel]
    rfl

theorem summable_integral_norm_boseTerm (m : ℕ) :
    Summable (fun k : ℕ =>
      ∫ t : ℝ in Ioi 0, ‖(boseTerm m k t : ℂ)‖) := by
  have h2 : Summable (fun n : ℕ => 1 / (n : ℝ) ^ 2) :=
    summable_one_div_nat_pow.mpr (by norm_num)
  have h3 : Summable (fun n : ℕ => 1 / (n : ℝ) ^ 3) :=
    summable_one_div_nat_pow.mpr (by norm_num)
  have hs := ((summable_nat_add_iff (m + 1)).2 h2).add
    ((summable_nat_add_iff (m + 1)).2 h3)
  apply hs.congr
  intro k
  rw [show (∫ t : ℝ in Ioi 0, ‖(boseTerm m k t : ℂ)‖) =
      ∫ t : ℝ in Ioi 0, boseTerm m k t by
        apply integral_congr_ae
        filter_upwards [ae_restrict_mem measurableSet_Ioi] with t ht
        rw [Complex.norm_real, Real.norm_eq_abs,
          abs_of_nonneg (by unfold boseTerm; positivity)]]
  rw [integral_boseTerm]
  push_cast
  ring

theorem integral_boseKernel_eq_tsum (m : ℕ) :
    (∫ t : ℝ in Ioi 0, (boseKernel m t : ℂ)) =
      ∑' k : ℕ,
        ((1 / (m + k + 1 : ℝ) ^ 2 +
          1 / (m + k + 1 : ℝ) ^ 3 : ℝ) : ℂ) := by
  let F : ℕ → ℝ → ℂ := fun k t => (boseTerm m k t : ℂ)
  have hFint : ∀ k, IntegrableOn (F k) (Ioi 0) := by
    intro k
    exact (integrableOn_boseTerm m k).ofReal
  have hswap := MeasureTheory.integral_tsum_of_summable_integral_norm
    hFint (summable_integral_norm_boseTerm m)
  rw [← hswap]
  · apply tsum_congr
    intro k
    rw [← integral_ofReal, integral_boseTerm]
  · apply integral_congr_ae
    filter_upwards [ae_restrict_mem measurableSet_Ioi] with t ht
    exact_mod_cast (hasSum_boseTerm (m := m) ht).tsum_eq

#print axioms integral_boseKernel_eq_tsum

end RamanujanChallenge.P27.Q6383
