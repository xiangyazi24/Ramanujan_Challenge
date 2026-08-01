/-
  The classical integral representation of Catalan's constant.
-/
import RamanujanChallenge.Problem25
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic
import Mathlib.MeasureTheory.Integral.DominatedConvergence

noncomputable section

namespace RamanujanChallenge.P25

open Filter MeasureTheory Set Topology
open scoped Interval

/-- The logarithmic moment used in the Catalan integral. -/
private theorem integral_pow_mul_neg_log (n : ℕ) :
    (∫ x : ℝ in 0..1, (-Real.log x) * x ^ n) =
      1 / ((n : ℝ) + 1) ^ 2 := by
  let F : ℝ → ℝ := fun x =>
    x ^ (n + 1) *
      ((((n : ℝ) + 1) * Real.log x - 1) / ((n : ℝ) + 1) ^ 2)
  have hInt :
      IntervalIntegrable (fun x : ℝ => x ^ n * Real.log x)
        volume 0 1 := by
    exact intervalIntegral.intervalIntegrable_log'.continuousOn_mul
      (continuousOn_pow n)
  have hlog :
      (∫ x : ℝ in 0..1, x ^ n * Real.log x) =
        -1 / ((n : ℝ) + 1) ^ 2 := by
    rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_tendsto
      (f := F) (fa := 0) (fb := -1 / ((n : ℝ) + 1) ^ 2) (hint := hInt)]
    · ring
    · norm_num
    · intro x hx
      have hx0 : x ≠ 0 := ne_of_gt hx.1
      have hpow :
          HasDerivAt (fun y : ℝ => y ^ (n + 1))
            (((n : ℝ) + 1) * x ^ n) x := by
        convert hasDerivAt_pow (n + 1) x using 1
        norm_num [Nat.cast_add]
      have hfactor :
          HasDerivAt
            (fun y : ℝ =>
              (((n : ℝ) + 1) * Real.log y - 1) /
                ((n : ℝ) + 1) ^ 2)
            ((((n : ℝ) + 1) * x⁻¹) / ((n : ℝ) + 1) ^ 2) x :=
        (((Real.hasDerivAt_log hx0).const_mul ((n : ℝ) + 1)).sub_const 1).div_const
          (((n : ℝ) + 1) ^ 2)
      dsimp [F]
      convert hpow.mul hfactor using 1
      field_simp [hx0]
      ring
    · have hpowlog :
          Tendsto (fun x : ℝ => Real.log x * x ^ (n + 1))
            (𝓝[>] 0) (𝓝 0) := by
        have h :=
          tendsto_log_mul_rpow_nhdsGT_zero
            (show (0 : ℝ) < (n : ℝ) + 1 by positivity)
        convert h using 1
        funext x
        rw [show (n : ℝ) + 1 = ((n + 1 : ℕ) : ℝ) by norm_num,
          Real.rpow_natCast]
      have hpow :
          Tendsto (fun x : ℝ => x ^ (n + 1))
            (𝓝[>] 0) (𝓝 0) := by
        simpa using tendsto_nhdsWithin_of_tendsto_nhds
          (continuousAt_pow (0 : ℝ) (n + 1)).tendsto
      have hmain :
          Tendsto
            (fun x : ℝ =>
              x ^ (n + 1) *
                ((((n : ℝ) + 1) * Real.log x - 1) /
                  ((n : ℝ) + 1) ^ 2))
            (𝓝[>] 0) (𝓝 0) := by
        convert
          (((hpowlog.const_mul ((n : ℝ) + 1)).sub hpow).div_const
            (((n : ℝ) + 1) ^ 2)) using 1
        · funext x
          field_simp
        · ring_nf
      simpa [F] using hmain
    · have hcont : ContinuousAt F 1 := by
        dsimp [F]
        exact (continuousAt_id.pow (n + 1)).mul
          ((((continuousAt_const.mul (Real.continuousAt_log (by norm_num))).sub
            continuousAt_const).div_const (((n : ℝ) + 1) ^ 2)))
      convert tendsto_nhdsWithin_of_tendsto_nhds hcont.tendsto using 1
      dsimp [F]
      norm_num
  calc
    (∫ x : ℝ in 0..1, (-Real.log x) * x ^ n) =
        -(∫ x : ℝ in 0..1, x ^ n * Real.log x) := by
          rw [← intervalIntegral.integral_neg]
          apply intervalIntegral.integral_congr
          intro x _
          ring
    _ = 1 / ((n : ℝ) + 1) ^ 2 := by rw [hlog]; ring

private def catalanIntegralTerm (n : ℕ) (x : ℝ) : ℝ :=
  (-Real.log x) * ((-1 : ℝ) ^ n * x ^ (2 * n))

private theorem catalanIntegralTerm_intervalIntegrable (n : ℕ) :
    IntervalIntegrable (catalanIntegralTerm n) volume 0 1 := by
  have hbase :
      IntervalIntegrable (fun x : ℝ => x ^ (2 * n) * Real.log x)
        volume 0 1 :=
    intervalIntegral.intervalIntegrable_log'.continuousOn_mul
      (continuousOn_pow (2 * n))
  unfold catalanIntegralTerm
  convert hbase.const_mul (-((-1 : ℝ) ^ n)) using 1
  funext x
  ring

private theorem catalanIntegralTerm_integral (n : ℕ) :
    (∫ x : ℝ in 0..1, catalanIntegralTerm n x) =
      (-1 : ℝ) ^ n / (2 * (n : ℝ) + 1) ^ 2 := by
  have hfun : catalanIntegralTerm n =
      (fun x : ℝ => (-1 : ℝ) ^ n * ((-Real.log x) * x ^ (2 * n))) := by
    funext x
    simp only [catalanIntegralTerm]
    ring
  rw [hfun]
  rw [intervalIntegral.integral_const_mul, integral_pow_mul_neg_log]
  simp [Nat.cast_mul, div_eq_mul_inv]

private theorem catalanIntegralTerm_integral_norm (n : ℕ) :
    (∫ x : ℝ in 0..1, ‖catalanIntegralTerm n x‖) =
      ‖(-1 : ℝ) ^ n / (2 * (n : ℝ) + 1) ^ 2‖ := by
  have hcongr :
      (∫ x : ℝ in 0..1, ‖catalanIntegralTerm n x‖) =
        ∫ x : ℝ in 0..1, (-Real.log x) * x ^ (2 * n) := by
    apply intervalIntegral.integral_congr
    intro x hx
    simp only [uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1), mem_Icc] at hx
    have hx0 : 0 ≤ x := hx.1
    have hlog : Real.log x ≤ 0 := Real.log_nonpos hx.1 hx.2
    change
      |(-Real.log x) * ((-1 : ℝ) ^ n * x ^ (2 * n))| =
        (-Real.log x) * x ^ (2 * n)
    simp only [abs_mul, abs_neg, abs_pow, abs_one, one_pow, one_mul,
      abs_of_nonneg hx0, abs_of_nonpos hlog]
  rw [hcongr, integral_pow_mul_neg_log, Real.norm_eq_abs, abs_div,
    abs_pow, abs_neg, abs_one, one_pow, one_div]
  rw [abs_of_nonneg (sq_nonneg (2 * (n : ℝ) + 1))]
  norm_num [Nat.cast_mul, one_div]

private theorem catalanIntegralTerm_integral_norm_summable :
    Summable
      (fun n : ℕ => ∫ x : ℝ in 0..1, ‖catalanIntegralTerm n x‖) := by
  exact catalanSeries_summable.norm.congr fun n =>
    (catalanIntegralTerm_integral_norm n).symm

private theorem catalanIntegralTerm_hasSum
    {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasSum (fun n : ℕ => catalanIntegralTerm n x)
      ((-Real.log x) / (1 + x ^ 2)) := by
  have hxnorm : ‖-(x ^ 2)‖ < (1 : ℝ) := by
    rw [Real.norm_eq_abs, abs_neg, abs_pow, abs_of_pos hx0]
    nlinarith
  have h :=
    (hasSum_geometric_of_norm_lt_one hxnorm).mul_left (-Real.log x)
  convert h using 1
  · funext n
    have hpow : (-(x ^ 2)) ^ n = (-1 : ℝ) ^ n * x ^ (2 * n) := by
      rw [neg_pow, pow_mul]
    rw [hpow]
    rfl
  · rw [div_eq_mul_inv]
    congr 2
    ring

/-- Catalan's defining alternating series equals its classical logarithmic integral. -/
theorem catalanConstant_eq_integral :
    catalanConstant = ∫ t in (0 : ℝ)..1, (-Real.log t) / (1 + t ^ 2) := by
  have hInt :
      ∀ n : ℕ,
        Integrable (catalanIntegralTerm n)
          (volume.restrict (Ioc (0 : ℝ) 1)) := by
    intro n
    exact (catalanIntegralTerm_intervalIntegrable n).1
  have hNorm :
      Summable
        (fun n : ℕ =>
          ∫ x : ℝ in Ioc (0 : ℝ) 1, ‖catalanIntegralTerm n x‖) := by
    simpa only [← intervalIntegral.integral_of_le
      (by norm_num : (0 : ℝ) ≤ 1)] using
      catalanIntegralTerm_integral_norm_summable
  have hswap :=
    MeasureTheory.hasSum_integral_of_summable_integral_norm
      (μ := volume.restrict (Ioc (0 : ℝ) 1)) hInt hNorm
  have hseries :
      HasSum
        (fun n : ℕ => (-1 : ℝ) ^ n / (2 * (n : ℝ) + 1) ^ 2)
        (∫ x : ℝ in Ioc (0 : ℝ) 1, ∑' n : ℕ, catalanIntegralTerm n x) := by
    convert hswap using 1
    funext n
    rw [← intervalIntegral.integral_of_le
      (by norm_num : (0 : ℝ) ≤ 1)]
    exact (catalanIntegralTerm_integral n).symm
  have hfinal :
      HasSum
        (fun n : ℕ => (-1 : ℝ) ^ n / (2 * (n : ℝ) + 1) ^ 2)
        (∫ x : ℝ in (0 : ℝ)..1, (-Real.log x) / (1 + x ^ 2)) := by
    convert hseries using 1
    rw [intervalIntegral.integral_of_le (by norm_num : (0 : ℝ) ≤ 1)]
    apply MeasureTheory.setIntegral_congr_ae measurableSet_Ioc
    filter_upwards [MeasureTheory.Measure.ae_ne volume (1 : ℝ)] with x hxne hx
    have hx1 : x < 1 := lt_of_le_of_ne hx.2 hxne
    exact (catalanIntegralTerm_hasSum hx.1 hx1).tsum_eq.symm
  exact catalanSeries_summable.hasSum.unique hfinal

end RamanujanChallenge.P25

end
