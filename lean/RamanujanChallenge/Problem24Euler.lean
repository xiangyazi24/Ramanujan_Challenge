import RamanujanChallenge.Problem24
import RamanujanChallenge.Problem26WeightThree
import Mathlib.MeasureTheory.Integral.DominatedConvergence
import Mathlib.MeasureTheory.Integral.IntervalIntegral.AbsolutelyContinuousFun

/-!
  Exact level-two Euler sums used in Ramanujan Challenge Problem 2.4.

  The proofs below use the logarithmic generating function of
  `parityRemainder24` and integrate it term by term.  The special values are
  reduced to the real dilogarithm and trilogarithm infrastructure developed
  for Problem 2.6.
-/

open Filter Set Topology
open scoped Interval

noncomputable section

def parityIncrement24 (n : ℕ) : ℝ :=
  (1 + 2 * (-1 : ℝ) ^ (n + 1)) / (n + 1 : ℝ)

theorem parityRemainder24_eq_sum_increment (n : ℕ) :
    parityRemainder24 (n + 1) =
      ∑ j ∈ Finset.range (n + 1), parityIncrement24 j := by
  unfold parityRemainder24 harmonicNumber signedHarmonic24 parityIncrement24
  rw [Finset.mul_sum, ← Finset.sum_add_distrib]
  apply Finset.sum_congr rfl
  intro j hj
  push_cast
  ring

theorem parityIncrement24_hasSum
    {x : ℝ} (hx : |x| < 1) (hxne : x ≠ 0) :
    HasSum (fun n : ℕ => parityIncrement24 n * x ^ n)
      ((-Real.log (1 - x) - 2 * Real.log (1 + x)) / x) := by
  have hminus := Real.hasSum_pow_div_log_of_abs_lt_one hx
  have hplus :=
    Real.hasSum_pow_div_log_of_abs_lt_one
      (show |-x| < 1 by simpa using hx)
  have hminus' := hminus.div_const x
  have hplus' := (hplus.mul_left 2).div_const x
  convert hminus'.add hplus' using 1
  · funext n
    unfold parityIncrement24
    field_simp [hxne]
    rw [neg_pow]
    ring
  · rw [show 1 - -x = 1 + x by ring]
    ring

theorem parityIncrement24_norm_summable
    {x : ℝ} (hx : |x| < 1) :
    Summable (fun n : ℕ => ‖parityIncrement24 n * x ^ n‖) := by
  have hgeom : Summable (fun n : ℕ => |x| ^ n) :=
    summable_geometric_of_lt_one (abs_nonneg x) hx
  refine (hgeom.mul_left 3).of_nonneg_of_le (fun n => norm_nonneg _) ?_
  intro n
  rw [Real.norm_eq_abs, abs_mul, abs_pow]
  have hcoeff :
      |parityIncrement24 n| ≤ 3 := by
    unfold parityIncrement24
    rw [abs_div]
    have hden : (1 : ℝ) ≤ |(n : ℝ) + 1| := by
      rw [abs_of_nonneg (by positivity)]
      norm_num
    have hnum :
        |1 + 2 * (-1 : ℝ) ^ (n + 1)| ≤ 3 := by
      calc
        |1 + 2 * (-1 : ℝ) ^ (n + 1)| ≤
            |(1 : ℝ)| + |2 * (-1 : ℝ) ^ (n + 1)| := abs_add_le _ _
        _ = (3 : ℝ) := by
          rw [abs_mul, abs_pow]
          norm_num
    calc
      |1 + 2 * (-1 : ℝ) ^ (n + 1)| / |(n : ℝ) + 1| ≤
          3 / 1 := by gcongr
      _ = 3 := by norm_num
  exact mul_le_mul_of_nonneg_right hcoeff (pow_nonneg (abs_nonneg x) n)

theorem parityRemainder24_generating_hasSum
    {x : ℝ} (hx : |x| < 1) (hxne : x ≠ 0) :
    HasSum (fun n : ℕ => parityRemainder24 (n + 1) * x ^ n)
      ((-Real.log (1 - x) - 2 * Real.log (1 + x)) /
        (x * (1 - x))) := by
  have hfSum := parityIncrement24_hasSum hx hxne
  have hgSum := hasSum_geometric_of_norm_lt_one
    (show ‖x‖ < 1 by simpa [Real.norm_eq_abs] using hx)
  have hproduct :
      HasSum
        (fun n => ∑ k ∈ Finset.range (n + 1),
          (parityIncrement24 k * x ^ k) * x ^ (n - k))
        ((∑' n, parityIncrement24 n * x ^ n) * ∑' n, x ^ n) := by
    apply hasSum_sum_range_mul_of_summable_norm
    · exact parityIncrement24_norm_summable hx
    · simpa [Real.norm_eq_abs, abs_pow] using
        (summable_geometric_of_lt_one (abs_nonneg x) hx)
  have hcoeff (n : ℕ) :
      (∑ k ∈ Finset.range (n + 1),
        (parityIncrement24 k * x ^ k) * x ^ (n - k)) =
        parityRemainder24 (n + 1) * x ^ n := by
    rw [parityRemainder24_eq_sum_increment]
    rw [Finset.sum_mul]
    apply Finset.sum_congr rfl
    intro k hk
    have hkn : k ≤ n := Nat.le_of_lt_succ (Finset.mem_range.mp hk)
    rw [mul_assoc, ← pow_add, Nat.add_sub_of_le hkn]
  have hseries :
      HasSum (fun n : ℕ => parityRemainder24 (n + 1) * x ^ n)
        ((∑' n, parityIncrement24 n * x ^ n) * ∑' n, x ^ n) := by
    exact HasSum.congr_fun hproduct fun n => (hcoeff n).symm
  rw [hfSum.tsum_eq, hgSum.tsum_eq] at hseries
  convert hseries using 1
  field_simp

private def quadraticLinearMomentTerm24 (n : ℕ) (x : ℝ) : ℝ :=
  -(parityRemainder24 (n + 1) * x ^ n * Real.log x)

private def quadraticLinearKernel24 (x : ℝ) : ℝ :=
  Real.log x *
    (Real.log (1 - x) + 2 * Real.log (1 + x)) /
      (x * (1 - x))

private theorem quadraticLinearMomentTerm24_integral (n : ℕ) :
    (∫ x : ℝ in 0..1, quadraticLinearMomentTerm24 n x) =
      quadraticLinearEulerTerm24 n := by
  unfold quadraticLinearMomentTerm24 quadraticLinearEulerTerm24
  rw [show
      (fun x : ℝ =>
        -(parityRemainder24 (n + 1) * x ^ n * Real.log x)) =
        (fun x : ℝ =>
          (-parityRemainder24 (n + 1)) * (x ^ n * Real.log x)) by
      funext x
      ring]
  rw [intervalIntegral.integral_const_mul,
    RamanujanChallenge.P26.integral_pow_mul_log26]
  ring

private theorem quadraticLinearMomentTerm24_intervalIntegrable (n : ℕ) :
    IntervalIntegrable (quadraticLinearMomentTerm24 n)
      MeasureTheory.volume 0 1 := by
  have hbase :
      IntervalIntegrable (fun x : ℝ => x ^ n * Real.log x)
        MeasureTheory.volume 0 1 :=
    intervalIntegral.intervalIntegrable_log'.continuousOn_mul
      (continuousOn_pow n)
  unfold quadraticLinearMomentTerm24
  convert (hbase.const_mul (-parityRemainder24 (n + 1))) using 1
  funext x
  ring

private theorem quadraticLinearMomentTerm24_integral_norm (n : ℕ) :
    (∫ x : ℝ in 0..1, ‖quadraticLinearMomentTerm24 n x‖) =
      ‖quadraticLinearEulerTerm24 n‖ := by
  have hcongr :
      (∫ x : ℝ in 0..1, ‖quadraticLinearMomentTerm24 n x‖) =
        ∫ x : ℝ in 0..1,
          |parityRemainder24 (n + 1)| *
            (-(x ^ n * Real.log x)) := by
    apply intervalIntegral.integral_congr
    intro x hx
    simp only [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1),
      Set.mem_Icc] at hx
    have hx0 : 0 ≤ x := hx.1
    have hx1 : x ≤ 1 := hx.2
    have hlog : Real.log x ≤ 0 := Real.log_nonpos hx0 hx1
    unfold quadraticLinearMomentTerm24
    change
      |-(parityRemainder24 (n + 1) * x ^ n * Real.log x)| =
        |parityRemainder24 (n + 1)| * -(x ^ n * Real.log x)
    rw [abs_neg, abs_mul, abs_mul, abs_pow,
      abs_of_nonneg hx0, abs_of_nonpos hlog]
    ring
  rw [hcongr, intervalIntegral.integral_const_mul,
    intervalIntegral.integral_neg,
    RamanujanChallenge.P26.integral_pow_mul_log26]
  unfold quadraticLinearEulerTerm24
  rw [Real.norm_eq_abs, abs_div,
    abs_of_nonneg (sq_nonneg (n + 1 : ℝ))]
  ring

private theorem quadraticLinearMomentTerm24_integral_norm_summable :
    Summable
      (fun n : ℕ =>
        ∫ x : ℝ in 0..1, ‖quadraticLinearMomentTerm24 n x‖) := by
  exact summable_quadraticLinearEulerTerm24.norm.congr fun n =>
    (quadraticLinearMomentTerm24_integral_norm n).symm

private theorem quadraticLinearMomentTerm24_hasSum_pointwise
    {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasSum (fun n : ℕ => quadraticLinearMomentTerm24 n x)
      (quadraticLinearKernel24 x) := by
  have h :=
    (parityRemainder24_generating_hasSum
      (by rw [abs_of_pos hx0]; exact hx1) (ne_of_gt hx0)).mul_left
      (-Real.log x)
  convert h using 1
  · funext n
    unfold quadraticLinearMomentTerm24
    ring
  · unfold quadraticLinearKernel24
    ring

private theorem quadraticLinearEulerTerm24_hasSum_integral :
    HasSum quadraticLinearEulerTerm24
      (∫ x : ℝ in 0..1, quadraticLinearKernel24 x) := by
  have hInt :
      ∀ n : ℕ,
        MeasureTheory.Integrable (quadraticLinearMomentTerm24 n)
          (MeasureTheory.volume.restrict (Set.Ioc 0 1)) := by
    intro n
    exact (quadraticLinearMomentTerm24_intervalIntegrable n).1
  have hNorm :
      Summable
        (fun n : ℕ =>
          ∫ x : ℝ in Set.Ioc 0 1,
            ‖quadraticLinearMomentTerm24 n x‖) := by
    simpa only [
      ← intervalIntegral.integral_of_le
        (by norm_num : (0 : ℝ) ≤ 1)] using
      quadraticLinearMomentTerm24_integral_norm_summable
  have h :=
    MeasureTheory.hasSum_integral_of_summable_integral_norm
      (μ := MeasureTheory.volume.restrict (Set.Ioc (0 : ℝ) 1))
      hInt hNorm
  have h' :
      HasSum quadraticLinearEulerTerm24
        (∫ x : ℝ in Set.Ioc 0 1,
          ∑' n : ℕ, quadraticLinearMomentTerm24 n x) := by
    convert h using 1
    funext n
    rw [← intervalIntegral.integral_of_le
      (by norm_num : (0 : ℝ) ≤ 1)]
    exact (quadraticLinearMomentTerm24_integral n).symm
  convert h' using 1
  rw [intervalIntegral.integral_of_le (by norm_num : (0 : ℝ) ≤ 1)]
  apply MeasureTheory.setIntegral_congr_ae measurableSet_Ioc
  filter_upwards [
    MeasureTheory.Measure.ae_ne MeasureTheory.volume (1 : ℝ)] with x hxne hx
  have hxlt : x < 1 := lt_of_le_of_ne hx.2 hxne
  exact (quadraticLinearMomentTerm24_hasSum_pointwise hx.1 hxlt).tsum_eq.symm

private def alternatingQuadraticLinearMomentTerm24
    (n : ℕ) (x : ℝ) : ℝ :=
  (-1 : ℝ) ^ (n + 1) * quadraticLinearMomentTerm24 n x

private def alternatingQuadraticLinearKernel24 (x : ℝ) : ℝ :=
  Real.log x *
    (Real.log (1 + x) + 2 * Real.log (1 - x)) /
      (x * (1 + x))

private theorem alternatingQuadraticLinearMomentTerm24_integral (n : ℕ) :
    (∫ x : ℝ in 0..1, alternatingQuadraticLinearMomentTerm24 n x) =
      alternatingQuadraticLinearEulerTerm24 n := by
  unfold alternatingQuadraticLinearMomentTerm24
  rw [intervalIntegral.integral_const_mul,
    quadraticLinearMomentTerm24_integral]
  rfl

private theorem alternatingQuadraticLinearMomentTerm24_intervalIntegrable
    (n : ℕ) :
    IntervalIntegrable (alternatingQuadraticLinearMomentTerm24 n)
      MeasureTheory.volume 0 1 := by
  unfold alternatingQuadraticLinearMomentTerm24
  exact
    (quadraticLinearMomentTerm24_intervalIntegrable n).const_mul
      ((-1 : ℝ) ^ (n + 1))

private theorem alternatingQuadraticLinearMomentTerm24_integral_norm
    (n : ℕ) :
    (∫ x : ℝ in 0..1,
      ‖alternatingQuadraticLinearMomentTerm24 n x‖) =
        ‖alternatingQuadraticLinearEulerTerm24 n‖ := by
  calc
    (∫ x : ℝ in 0..1,
        ‖alternatingQuadraticLinearMomentTerm24 n x‖) =
        ∫ x : ℝ in 0..1, ‖quadraticLinearMomentTerm24 n x‖ := by
          apply intervalIntegral.integral_congr
          intro x hx
          unfold alternatingQuadraticLinearMomentTerm24
          change
            ‖(-1 : ℝ) ^ (n + 1) *
                quadraticLinearMomentTerm24 n x‖ =
              ‖quadraticLinearMomentTerm24 n x‖
          rw [norm_mul]
          simp
    _ = ‖quadraticLinearEulerTerm24 n‖ :=
      quadraticLinearMomentTerm24_integral_norm n
    _ = ‖alternatingQuadraticLinearEulerTerm24 n‖ := by
      unfold alternatingQuadraticLinearEulerTerm24
      rw [norm_mul]
      simp

private theorem
    alternatingQuadraticLinearMomentTerm24_integral_norm_summable :
    Summable
      (fun n : ℕ =>
        ∫ x : ℝ in 0..1,
          ‖alternatingQuadraticLinearMomentTerm24 n x‖) := by
  exact summable_alternatingQuadraticLinearEulerTerm24.norm.congr fun n =>
    (alternatingQuadraticLinearMomentTerm24_integral_norm n).symm

private theorem alternatingQuadraticLinearMomentTerm24_hasSum_pointwise
    {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasSum (fun n : ℕ => alternatingQuadraticLinearMomentTerm24 n x)
      (alternatingQuadraticLinearKernel24 x) := by
  have hxabs : |x| < 1 := by
    rw [abs_of_pos hx0]
    exact hx1
  have h :=
    (parityRemainder24_generating_hasSum
      (x := -x) (by simpa using hxabs)
      (neg_ne_zero.mpr (ne_of_gt hx0))).mul_left (Real.log x)
  convert h using 1
  · funext n
    unfold alternatingQuadraticLinearMomentTerm24
      quadraticLinearMomentTerm24
    rw [neg_pow]
    ring
  · unfold alternatingQuadraticLinearKernel24
    rw [show 1 - -x = 1 + x by ring,
      show 1 + -x = 1 - x by ring]
    have hxne : x ≠ 0 := ne_of_gt hx0
    have h1xne : 1 + x ≠ 0 := by linarith
    field_simp [hxne, h1xne]
    ring

private theorem alternatingQuadraticLinearEulerTerm24_hasSum_integral :
    HasSum alternatingQuadraticLinearEulerTerm24
      (∫ x : ℝ in 0..1, alternatingQuadraticLinearKernel24 x) := by
  have hInt :
      ∀ n : ℕ,
        MeasureTheory.Integrable
          (alternatingQuadraticLinearMomentTerm24 n)
          (MeasureTheory.volume.restrict (Set.Ioc 0 1)) := by
    intro n
    exact
      (alternatingQuadraticLinearMomentTerm24_intervalIntegrable n).1
  have hNorm :
      Summable
        (fun n : ℕ =>
          ∫ x : ℝ in Set.Ioc 0 1,
            ‖alternatingQuadraticLinearMomentTerm24 n x‖) := by
    simpa only [
      ← intervalIntegral.integral_of_le
        (by norm_num : (0 : ℝ) ≤ 1)] using
      alternatingQuadraticLinearMomentTerm24_integral_norm_summable
  have h :=
    MeasureTheory.hasSum_integral_of_summable_integral_norm
      (μ := MeasureTheory.volume.restrict (Set.Ioc (0 : ℝ) 1))
      hInt hNorm
  have h' :
      HasSum alternatingQuadraticLinearEulerTerm24
        (∫ x : ℝ in Set.Ioc 0 1,
          ∑' n : ℕ, alternatingQuadraticLinearMomentTerm24 n x) := by
    convert h using 1
    funext n
    rw [← intervalIntegral.integral_of_le
      (by norm_num : (0 : ℝ) ≤ 1)]
    exact (alternatingQuadraticLinearMomentTerm24_integral n).symm
  convert h' using 1
  rw [intervalIntegral.integral_of_le (by norm_num : (0 : ℝ) ≤ 1)]
  apply MeasureTheory.setIntegral_congr_ae measurableSet_Ioc
  filter_upwards [
    MeasureTheory.Measure.ae_ne MeasureTheory.volume (1 : ℝ)] with x hxne hx
  have hxlt : x < 1 := lt_of_le_of_ne hx.2 hxne
  exact
    (alternatingQuadraticLinearMomentTerm24_hasSum_pointwise
      hx.1 hxlt).tsum_eq.symm

private def harmonicIncrement24 (n : ℕ) : ℝ :=
  1 / (n + 1 : ℝ)

private theorem harmonicNumber_eq_sum_increment24 (n : ℕ) :
    harmonicNumber (n + 1) =
      ∑ j ∈ Finset.range (n + 1), harmonicIncrement24 j := by
  unfold harmonicNumber harmonicIncrement24
  apply Finset.sum_congr rfl
  intro j hj
  ring

private theorem harmonicIncrement24_hasSum
    {x : ℝ} (hx : |x| < 1) (hxne : x ≠ 0) :
    HasSum (fun n : ℕ => harmonicIncrement24 n * x ^ n)
      (-Real.log (1 - x) / x) := by
  have h := Real.hasSum_pow_div_log_of_abs_lt_one hx
  convert h.div_const x using 1
  · funext n
    unfold harmonicIncrement24
    field_simp [hxne]
    ring

private theorem harmonicIncrement24_norm_summable
    {x : ℝ} (hx : |x| < 1) :
    Summable (fun n : ℕ => ‖harmonicIncrement24 n * x ^ n‖) := by
  have hgeom : Summable (fun n : ℕ => |x| ^ n) :=
    summable_geometric_of_lt_one (abs_nonneg x) hx
  refine hgeom.of_nonneg_of_le (fun n => norm_nonneg _) ?_
  intro n
  rw [Real.norm_eq_abs, abs_mul, abs_pow]
  have hinc : |harmonicIncrement24 n| ≤ 1 := by
    unfold harmonicIncrement24
    rw [abs_div, abs_one, abs_of_nonneg (by positivity)]
    exact (div_le_one (by positivity)).2 (by norm_num)
  simpa using mul_le_mul_of_nonneg_right hinc
    (pow_nonneg (abs_nonneg x) n)

/-- `∑_{n≥0} H_{n+1} xⁿ = -log(1-x)/(x(1-x))` for `|x| < 1`, `x ≠ 0`.

Made public (was private) because Layer E of `Problem24QuadraticAlt` needs it at
`x = -t`: it is exactly the generating function of the `K` integrand,
`log(1+t)/(1+t) = ∑_{n≥0} (-1)ⁿ H_{n+1} tⁿ⁺¹`. -/
theorem harmonicNumber_generating_hasSum
    {x : ℝ} (hx : |x| < 1) (hxne : x ≠ 0) :
    HasSum (fun n : ℕ => harmonicNumber (n + 1) * x ^ n)
      (-Real.log (1 - x) / (x * (1 - x))) := by
  have hfSum := harmonicIncrement24_hasSum hx hxne
  have hgSum := hasSum_geometric_of_norm_lt_one
    (show ‖x‖ < 1 by simpa [Real.norm_eq_abs] using hx)
  have hproduct :
      HasSum
        (fun n => ∑ k ∈ Finset.range (n + 1),
          (harmonicIncrement24 k * x ^ k) * x ^ (n - k))
        ((∑' n, harmonicIncrement24 n * x ^ n) * ∑' n, x ^ n) := by
    apply hasSum_sum_range_mul_of_summable_norm
    · exact harmonicIncrement24_norm_summable hx
    · simpa [Real.norm_eq_abs, abs_pow] using
        (summable_geometric_of_lt_one (abs_nonneg x) hx)
  have hcoeff (n : ℕ) :
      (∑ k ∈ Finset.range (n + 1),
        (harmonicIncrement24 k * x ^ k) * x ^ (n - k)) =
          harmonicNumber (n + 1) * x ^ n := by
    rw [harmonicNumber_eq_sum_increment24, Finset.sum_mul]
    apply Finset.sum_congr rfl
    intro k hk
    have hkn : k ≤ n := Nat.le_of_lt_succ (Finset.mem_range.mp hk)
    rw [mul_assoc, ← pow_add, Nat.add_sub_of_le hkn]
  have hseries :
      HasSum (fun n : ℕ => harmonicNumber (n + 1) * x ^ n)
        ((∑' n, harmonicIncrement24 n * x ^ n) * ∑' n, x ^ n) := by
    exact HasSum.congr_fun hproduct fun n => (hcoeff n).symm
  rw [hfSum.tsum_eq, hgSum.tsum_eq] at hseries
  convert hseries using 1
  field_simp

private def pairedHarmonicMomentTerm24 (m : ℕ) (x : ℝ) : ℝ :=
  -harmonicNumber (2 * m + 1) * x ^ (2 * m) +
    harmonicNumber (2 * m + 2) * x ^ (2 * m + 1)

private def pairedHarmonicKernel24 (x : ℝ) : ℝ :=
  -Real.log (1 + x) / (x * (1 + x))

private theorem pairedHarmonicMomentTerm24_integral (m : ℕ) :
    (∫ x : ℝ in 0..1, pairedHarmonicMomentTerm24 m x) =
      pairedAlternatingHarmonicEulerTerm24 m := by
  unfold pairedHarmonicMomentTerm24
    pairedAlternatingHarmonicEulerTerm24
  have hleft :
      IntervalIntegrable
        (fun x : ℝ => -harmonicNumber (2 * m + 1) * x ^ (2 * m))
        MeasureTheory.volume 0 1 :=
    (by fun_prop : Continuous
      (fun x : ℝ =>
        -harmonicNumber (2 * m + 1) * x ^ (2 * m))).intervalIntegrable 0 1
  have hright :
      IntervalIntegrable
        (fun x : ℝ => harmonicNumber (2 * m + 2) * x ^ (2 * m + 1))
        MeasureTheory.volume 0 1 :=
    (by fun_prop : Continuous
      (fun x : ℝ =>
        harmonicNumber (2 * m + 2) * x ^ (2 * m + 1))).intervalIntegrable 0 1
  rw [intervalIntegral.integral_add hleft hright,
    intervalIntegral.integral_const_mul,
    intervalIntegral.integral_const_mul,
    integral_pow, integral_pow]
  push_cast
  ring

private theorem pairedHarmonicMomentTerm24_intervalIntegrable (m : ℕ) :
    IntervalIntegrable (pairedHarmonicMomentTerm24 m)
      MeasureTheory.volume 0 1 := by
  unfold pairedHarmonicMomentTerm24
  exact
    ((by fun_prop : Continuous
      (fun x : ℝ =>
        -harmonicNumber (2 * m + 1) * x ^ (2 * m) +
          harmonicNumber (2 * m + 2) * x ^ (2 * m + 1))).intervalIntegrable 0 1)

private theorem pairedHarmonicMomentTerm24_hasSum_pointwise
    {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasSum (fun m : ℕ => pairedHarmonicMomentTerm24 m x)
      (pairedHarmonicKernel24 x) := by
  have hxabs : |x| < 1 := by
    rw [abs_of_pos hx0]
    exact hx1
  have hraw :=
    (harmonicNumber_generating_hasSum
      (x := -x) (by simpa using hxabs)
      (neg_ne_zero.mpr (ne_of_gt hx0))).mul_left (-1)
  have hpairs := hraw.pair_consecutive24
  convert hpairs using 1
  · funext m
    unfold pairedHarmonicMomentTerm24
    have hodd : (-1 : ℝ) ^ (2 * m) = 1 := by
      rw [pow_mul]
      norm_num
    have heven : (-1 : ℝ) ^ (2 * m + 1) = -1 := by
      rw [pow_add, pow_mul]
      norm_num
    simp only [neg_pow, hodd, heven]
    rw [show 2 * m + 1 + 1 = 2 * m + 2 by omega]
    ring
  · unfold pairedHarmonicKernel24
    rw [show 1 - -x = 1 + x by ring]
    have hxne : x ≠ 0 := ne_of_gt hx0
    have h1xne : 1 + x ≠ 0 := by linarith
    field_simp [hxne, h1xne]

private def pairedHarmonicMomentMajorant24 (m : ℕ) (x : ℝ) : ℝ :=
  harmonicNumber (2 * m + 1) * x ^ (2 * m) * (1 - x) +
    x ^ (2 * m + 1) / (2 * m + 2 : ℝ)

private def pairedHarmonicMajorantValue24 (m : ℕ) : ℝ :=
  harmonicNumber (2 * m + 1) *
      (1 / (2 * m + 1 : ℝ) - 1 / (2 * m + 2 : ℝ)) +
    1 / (2 * m + 2 : ℝ) ^ 2

private theorem pairedHarmonicMomentTerm24_norm_le_majorant
    (m : ℕ) {x : ℝ} (hx0 : 0 ≤ x) (hx1 : x ≤ 1) :
    ‖pairedHarmonicMomentTerm24 m x‖ ≤
      pairedHarmonicMomentMajorant24 m x := by
  have hsucc :
      harmonicNumber (2 * m + 2) =
        harmonicNumber (2 * m + 1) +
          1 / (2 * m + 2 : ℝ) := by
    rw [show 2 * m + 2 = (2 * m + 1) + 1 by omega,
      harmonicNumber_succ]
    push_cast
    ring
  have hform :
      pairedHarmonicMomentTerm24 m x =
        -harmonicNumber (2 * m + 1) * x ^ (2 * m) * (1 - x) +
          x ^ (2 * m + 1) / (2 * m + 2 : ℝ) := by
    unfold pairedHarmonicMomentTerm24
    rw [hsucc]
    ring
  rw [hform, Real.norm_eq_abs]
  calc
    |-harmonicNumber (2 * m + 1) * x ^ (2 * m) * (1 - x) +
        x ^ (2 * m + 1) / (2 * m + 2 : ℝ)| ≤
        |-harmonicNumber (2 * m + 1) * x ^ (2 * m) * (1 - x)| +
          |x ^ (2 * m + 1) / (2 * m + 2 : ℝ)| := abs_add_le _ _
    _ = pairedHarmonicMomentMajorant24 m x := by
      unfold pairedHarmonicMomentMajorant24
      rw [abs_mul, abs_mul, abs_neg,
        abs_of_nonneg (harmonicNumber_nonneg (2 * m + 1)),
        abs_pow, abs_of_nonneg hx0,
        abs_of_nonneg (sub_nonneg.mpr hx1),
        abs_div, abs_pow, abs_of_nonneg hx0,
        abs_of_nonneg (by positivity : (0 : ℝ) ≤ 2 * m + 2)]

private theorem pairedHarmonicMomentMajorant24_integral (m : ℕ) :
    (∫ x : ℝ in 0..1, pairedHarmonicMomentMajorant24 m x) =
      pairedHarmonicMajorantValue24 m := by
  unfold pairedHarmonicMomentMajorant24
    pairedHarmonicMajorantValue24
  have hleft :
      IntervalIntegrable
        (fun x : ℝ =>
          harmonicNumber (2 * m + 1) * x ^ (2 * m) * (1 - x))
        MeasureTheory.volume 0 1 :=
    (by fun_prop : Continuous
      (fun x : ℝ =>
        harmonicNumber (2 * m + 1) * x ^ (2 * m) * (1 - x))).intervalIntegrable 0 1
  have hright :
      IntervalIntegrable
        (fun x : ℝ => x ^ (2 * m + 1) / (2 * m + 2 : ℝ))
        MeasureTheory.volume 0 1 :=
    (by fun_prop : Continuous
      (fun x : ℝ => x ^ (2 * m + 1) / (2 * m + 2 : ℝ))).intervalIntegrable 0 1
  rw [intervalIntegral.integral_add hleft hright]
  rw [show
      (fun x : ℝ =>
        harmonicNumber (2 * m + 1) * x ^ (2 * m) * (1 - x)) =
        (fun x : ℝ =>
          harmonicNumber (2 * m + 1) *
            (x ^ (2 * m) - x ^ (2 * m + 1))) by
      funext x
      rw [pow_succ]
      ring]
  rw [show
      (fun x : ℝ => x ^ (2 * m + 1) / (2 * m + 2 : ℝ)) =
        (fun x : ℝ =>
          (1 / (2 * m + 2 : ℝ)) * x ^ (2 * m + 1)) by
      funext x
      ring]
  rw [intervalIntegral.integral_const_mul,
    intervalIntegral.integral_sub
      ((continuous_pow (2 * m)).intervalIntegrable 0 1)
      ((continuous_pow (2 * m + 1)).intervalIntegrable 0 1),
    intervalIntegral.integral_const_mul,
    integral_pow]
  rw [integral_pow]
  push_cast
  field_simp
  ring

private theorem pairedHarmonicMajorantValue24_nonneg (m : ℕ) :
    0 ≤ pairedHarmonicMajorantValue24 m := by
  unfold pairedHarmonicMajorantValue24
  have hg :
      0 ≤ 1 / (2 * m + 1 : ℝ) - 1 / (2 * m + 2 : ℝ) := by
    apply sub_nonneg.mpr
    apply one_div_le_one_div_of_le
    · positivity
    · norm_num
  exact add_nonneg
    (mul_nonneg (harmonicNumber_nonneg (2 * m + 1)) hg)
    (by positivity)

private theorem pairedHarmonicMajorantValue24_le (m : ℕ) :
    pairedHarmonicMajorantValue24 m ≤
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
    (harmonicNumber_mono
      (by omega : 2 * m + 1 ≤ 2 * (m + 1))).trans
      (harmonicNumber_two_mul_le (m + 1))
  have heven_le :
      1 / (2 * m + 2 : ℝ) ^ 2 ≤ 1 / (m + 1 : ℝ) ^ 2 := by
    apply one_div_le_one_div_of_le (by positivity)
    nlinarith [sq_nonneg (m : ℝ)]
  unfold pairedHarmonicMajorantValue24
  change
    harmonicNumber (2 * m + 1) * g +
        1 / (2 * m + 2 : ℝ) ^ 2 ≤ _
  calc
    harmonicNumber (2 * m + 1) * g +
        1 / (2 * m + 2 : ℝ) ^ 2 ≤
        (2 * harmonicNumber (m + 1)) *
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

private theorem pairedHarmonicMajorantValue24_summable :
    Summable pairedHarmonicMajorantValue24 := by
  exact
    (summable_harmonicNumber_succ_sq_div.mul_left 3).of_nonneg_of_le
      pairedHarmonicMajorantValue24_nonneg
      pairedHarmonicMajorantValue24_le

private theorem pairedHarmonicMomentTerm24_integral_norm_le (m : ℕ) :
    (∫ x : ℝ in 0..1, ‖pairedHarmonicMomentTerm24 m x‖) ≤
      pairedHarmonicMajorantValue24 m := by
  have hleft :
      IntervalIntegrable
        (fun x : ℝ => ‖pairedHarmonicMomentTerm24 m x‖)
        MeasureTheory.volume 0 1 := by
    have hcont :
        Continuous (fun x : ℝ => pairedHarmonicMomentTerm24 m x) := by
      unfold pairedHarmonicMomentTerm24
      fun_prop
    exact hcont.norm.intervalIntegrable 0 1
  have hright :
      IntervalIntegrable (pairedHarmonicMomentMajorant24 m)
        MeasureTheory.volume 0 1 :=
    (by
      unfold pairedHarmonicMomentMajorant24
      fun_prop : Continuous
        (pairedHarmonicMomentMajorant24 m)).intervalIntegrable 0 1
  calc
    (∫ x : ℝ in 0..1, ‖pairedHarmonicMomentTerm24 m x‖) ≤
        ∫ x : ℝ in 0..1, pairedHarmonicMomentMajorant24 m x :=
      intervalIntegral.integral_mono_on (by norm_num) hleft hright
        (fun x hx =>
          pairedHarmonicMomentTerm24_norm_le_majorant m hx.1 hx.2)
    _ = pairedHarmonicMajorantValue24 m :=
      pairedHarmonicMomentMajorant24_integral m

private theorem pairedHarmonicMomentTerm24_integral_norm_summable :
    Summable
      (fun m : ℕ =>
        ∫ x : ℝ in 0..1, ‖pairedHarmonicMomentTerm24 m x‖) := by
  exact pairedHarmonicMajorantValue24_summable.of_nonneg_of_le
    (fun m =>
      intervalIntegral.integral_nonneg_of_forall
        (by norm_num : (0 : ℝ) ≤ 1) fun x => norm_nonneg _)
    pairedHarmonicMomentTerm24_integral_norm_le

private theorem pairedAlternatingHarmonicEulerTerm24_hasSum_integral :
    HasSum pairedAlternatingHarmonicEulerTerm24
      (∫ x : ℝ in 0..1, pairedHarmonicKernel24 x) := by
  have hInt :
      ∀ m : ℕ,
        MeasureTheory.Integrable (pairedHarmonicMomentTerm24 m)
          (MeasureTheory.volume.restrict (Set.Ioc 0 1)) := by
    intro m
    exact (pairedHarmonicMomentTerm24_intervalIntegrable m).1
  have hNorm :
      Summable
        (fun m : ℕ =>
          ∫ x : ℝ in Set.Ioc 0 1,
            ‖pairedHarmonicMomentTerm24 m x‖) := by
    simpa only [
      ← intervalIntegral.integral_of_le
        (by norm_num : (0 : ℝ) ≤ 1)] using
      pairedHarmonicMomentTerm24_integral_norm_summable
  have h :=
    MeasureTheory.hasSum_integral_of_summable_integral_norm
      (μ := MeasureTheory.volume.restrict (Set.Ioc (0 : ℝ) 1))
      hInt hNorm
  have h' :
      HasSum pairedAlternatingHarmonicEulerTerm24
        (∫ x : ℝ in Set.Ioc 0 1,
          ∑' m : ℕ, pairedHarmonicMomentTerm24 m x) := by
    convert h using 1
    funext m
    rw [← intervalIntegral.integral_of_le
      (by norm_num : (0 : ℝ) ≤ 1)]
    exact (pairedHarmonicMomentTerm24_integral m).symm
  convert h' using 1
  rw [intervalIntegral.integral_of_le (by norm_num : (0 : ℝ) ≤ 1)]
  apply MeasureTheory.setIntegral_congr_ae measurableSet_Ioc
  filter_upwards [
    MeasureTheory.Measure.ae_ne MeasureTheory.volume (1 : ℝ)] with x hxne hx
  have hxlt : x < 1 := lt_of_le_of_ne hx.2 hxne
  exact
    (pairedHarmonicMomentTerm24_hasSum_pointwise
      hx.1 hxlt).tsum_eq.symm

private def pairedHarmonicPrimitive24 (x : ℝ) : ℝ :=
  dilog (-x) + (1 / 2 : ℝ) * Real.log (1 + x) ^ 2

private theorem pairedHarmonicPrimitive24_hasDerivAt
    {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasDerivAt pairedHarmonicPrimitive24
      (pairedHarmonicKernel24 x) x := by
  have hxne : x ≠ 0 := ne_of_gt hx0
  have h1xne : 1 + x ≠ 0 := by linarith
  have hxabs : |x| < 1 := by
    rw [abs_of_pos hx0]
    exact hx1
  have hdRaw :=
    (dilog_hasDerivAt_of_abs_lt_one
      (show |-x| < 1 by simpa using hxabs)
      (neg_ne_zero.mpr hxne)).comp x (hasDerivAt_neg x)
  have hd :
      HasDerivAt (fun y : ℝ => dilog (-y))
        (-Real.log (1 + x) / x) x := by
    convert hdRaw using 1
    field_simp [hxne]
    ring
  have hinner :
      HasDerivAt (fun y : ℝ => 1 + y) 1 x := by
    convert (hasDerivAt_const x (1 : ℝ)).add (hasDerivAt_id x) using 1
    ring
  have hlog :
      HasDerivAt (fun y : ℝ => Real.log (1 + y))
        (1 / (1 + x)) x := by
    convert hinner.log h1xne using 1
  unfold pairedHarmonicPrimitive24 pairedHarmonicKernel24
  convert hd.add ((hlog.pow 2).const_mul (1 / 2 : ℝ)) using 1
  field_simp [hxne, h1xne]
  ring

private theorem pairedHarmonicPrimitive24_continuousOn :
    ContinuousOn pairedHarmonicPrimitive24 (Icc (0 : ℝ) 1) := by
  have hd :
      ContinuousOn (fun x : ℝ => dilog (-x)) (Icc (0 : ℝ) 1) := by
    apply dilog_continuousOn_unit.comp
    · fun_prop
    · intro x hx
      constructor <;> linarith [hx.1, hx.2]
  have hlog :
      ContinuousOn (fun x : ℝ => Real.log (1 + x))
        (Icc (0 : ℝ) 1) := by
    apply (continuousOn_const.add continuousOn_id).log
    intro x hx
    simpa only [Pi.add_apply, Pi.one_apply, id_eq] using
      (show (1 : ℝ) + x ≠ 0 by linarith [hx.1])
  unfold pairedHarmonicPrimitive24
  exact hd.add (continuousOn_const.mul (hlog.pow 2))

private theorem pairedHarmonicKernel24_intervalIntegrable :
    IntervalIntegrable pairedHarmonicKernel24
      MeasureTheory.volume 0 1 := by
  have hfactor :
      ContinuousOn
        (fun x : ℝ =>
          -RamanujanChallenge.P26.logOnePlusSlope26 x / (1 + x))
        (Icc (0 : ℝ) 1) := by
    apply ContinuousOn.div
    · exact
        RamanujanChallenge.P26.logOnePlusSlope26_continuousOn.neg
    · fun_prop
    · intro x hx
      linarith [hx.1]
  have hint :
      IntervalIntegrable
        (fun x : ℝ =>
          -RamanujanChallenge.P26.logOnePlusSlope26 x / (1 + x))
        MeasureTheory.volume 0 1 := by
    have hfactor' :
        ContinuousOn
          (fun x : ℝ =>
            -RamanujanChallenge.P26.logOnePlusSlope26 x / (1 + x))
          (Set.uIcc (0 : ℝ) 1) := by
      simpa [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using
        hfactor
    exact hfactor'.intervalIntegrable
  apply IntervalIntegrable.congr
    (f := fun x : ℝ =>
      -RamanujanChallenge.P26.logOnePlusSlope26 x / (1 + x)) ?_ hint
  intro x hx
  have hx' : x ∈ Ioc (0 : ℝ) 1 := by
    simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hx
  have hxne : x ≠ 0 := ne_of_gt hx'.1
  unfold pairedHarmonicKernel24
  simp [RamanujanChallenge.P26.logOnePlusSlope26, hxne]
  field_simp [hxne]

private theorem pairedHarmonicKernel24_integral :
    (∫ x : ℝ in 0..1, pairedHarmonicKernel24 x) =
      pairedAlternatingHarmonicEulerValue24 := by
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le
    (f := pairedHarmonicPrimitive24)
    (f' := pairedHarmonicKernel24)
    (by norm_num)
    pairedHarmonicPrimitive24_continuousOn
    (fun x hx =>
      pairedHarmonicPrimitive24_hasDerivAt hx.1 hx.2)
    pairedHarmonicKernel24_intervalIntegrable]
  unfold pairedHarmonicPrimitive24
    pairedAlternatingHarmonicEulerValue24
  rw [show -(0 : ℝ) = 0 by ring, dilog_zero,
    RamanujanChallenge.P26.dilog_neg_one26]
  norm_num
  ring

theorem pairedAlternatingHarmonicEulerTerm24_hasSum :
    HasSum pairedAlternatingHarmonicEulerTerm24
      pairedAlternatingHarmonicEulerValue24 := by
  rw [← pairedHarmonicKernel24_integral]
  exact pairedAlternatingHarmonicEulerTerm24_hasSum_integral

/-! ## Ordinary weight-three logarithmic moments -/

private def logOneMinusSlope24 (x : ℝ) : ℝ :=
  Function.update (fun y : ℝ => Real.log (1 - y) / y) 0 (-1) x

private theorem logOneMinusSlope24_continuousOn :
    ContinuousOn logOneMinusSlope24 (Icc (0 : ℝ) (1 / 2)) := by
  intro x hx
  by_cases hxzero : x = 0
  · subst x
    have hlog :
        HasDerivAt (fun y : ℝ => Real.log (1 - y)) (-1) 0 := by
      have hone :
          HasDerivAt (fun y : ℝ => 1 - y) (-1) 0 := by
        convert (hasDerivAt_const (0 : ℝ) 1).sub (hasDerivAt_id 0) using 1
        norm_num
      simpa [Function.comp_def] using
        (HasDerivAt.comp (h := fun y : ℝ => 1 - y) 0
          (Real.hasDerivAt_log
            (by norm_num : (fun y : ℝ => 1 - y) 0 ≠ 0)) hone)
    have hc := hlog.continuousAt_div
    have hc' : ContinuousAt logOneMinusSlope24 0 := by
      convert hc using 1
      funext y
      simp [logOneMinusSlope24]
    exact hc'.continuousWithinAt
  · have hbase :
        ContinuousAt (fun y : ℝ => Real.log (1 - y) / y) x := by
      have hone : ContinuousAt (fun y : ℝ => 1 - y) x := by
        fun_prop
      have hlog :
          ContinuousAt (fun y : ℝ => Real.log (1 - y)) x :=
        (Real.continuousAt_log
          (by linarith [hx.2] : 1 - x ≠ 0)).comp hone
      exact hlog.div continuousAt_id hxzero
    have heq :
        logOneMinusSlope24 =ᶠ[𝓝 x]
          (fun y : ℝ => Real.log (1 - y) / y) := by
      filter_upwards [eventually_ne_nhds hxzero] with y hy
      simp [logOneMinusSlope24, hy]
    exact (hbase.congr_of_eventuallyEq heq).continuousWithinAt

private def minusRadialKernel24 (x : ℝ) : ℝ :=
  Real.log x * Real.log (1 - x) / x

private def minusRadialPrimitive24 (x : ℝ) : ℝ :=
  -(Real.log x * dilog x) + RamanujanChallenge.P26.trilog26 x

private theorem minusRadialPrimitive24_hasDerivAt
    {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasDerivAt minusRadialPrimitive24 (minusRadialKernel24 x) x := by
  have hxne : x ≠ 0 := ne_of_gt hx0
  have hxabs : |x| < 1 := by
    rw [abs_of_pos hx0]
    exact hx1
  have hd := dilog_hasDerivAt hx0 hx1
  have ht :=
    RamanujanChallenge.P26.trilog26_hasDerivAt_of_abs_lt_one hxabs hxne
  unfold minusRadialPrimitive24
  convert ((Real.hasDerivAt_log hxne).mul hd).neg.add ht using 1
  unfold minusRadialKernel24
  field_simp [hxne]
  ring

private theorem minusRadialKernel24_intervalIntegrable :
    IntervalIntegrable minusRadialKernel24
      MeasureTheory.volume 0 1 := by
  apply IntervalIntegrable.trans (b := (1 / 2 : ℝ))
  · have hlog :
        IntervalIntegrable Real.log MeasureTheory.volume (0 : ℝ) (1 / 2) :=
      intervalIntegral.intervalIntegrable_log'
    have hslope :
        ContinuousOn logOneMinusSlope24
          (Set.uIcc (0 : ℝ) (1 / 2)) := by
      simpa [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2)] using
        logOneMinusSlope24_continuousOn
    have hint := hlog.continuousOn_mul hslope
    apply IntervalIntegrable.congr
      (f := fun x : ℝ => logOneMinusSlope24 x * Real.log x) ?_ hint
    intro x hx
    have hx' : x ∈ Ioc (0 : ℝ) (1 / 2) := by
      simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2)] using hx
    have hxne : x ≠ 0 := ne_of_gt hx'.1
    simp [minusRadialKernel24, logOneMinusSlope24, hxne]
    field_simp [hxne]
  · have hlogSub :
        IntervalIntegrable (fun x : ℝ => Real.log (1 - x))
          MeasureTheory.volume (1 / 2) 1 := by
      have h :=
        (intervalIntegral.intervalIntegrable_log'
          (a := (0 : ℝ)) (b := (1 / 2 : ℝ))).comp_sub_left 1
      convert h.symm using 1 <;> norm_num
    have hfactor :
        ContinuousOn (fun x : ℝ => Real.log x / x)
          (Set.uIcc (1 / 2 : ℝ) 1) := by
      have hfactorIcc :
          ContinuousOn (fun x : ℝ => Real.log x / x)
            (Icc (1 / 2 : ℝ) 1) := by
        apply ContinuousOn.div
        · apply continuousOn_id.log
          intro x hx
          simpa only [id_eq] using (show x ≠ 0 by linarith [hx.1])
        · exact continuousOn_id
        · intro x hx
          linarith [hx.1]
      rw [Set.uIcc_of_le (by norm_num : (1 / 2 : ℝ) ≤ 1)]
      exact hfactorIcc
    have hint := hlogSub.continuousOn_mul hfactor
    apply IntervalIntegrable.congr
      (f := fun x : ℝ => (Real.log x / x) * Real.log (1 - x)) ?_ hint
    intro x _
    unfold minusRadialKernel24
    ring

private theorem minusRadialPrimitive24_tendsto_zero :
    Tendsto minusRadialPrimitive24 (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  have hslope :
      Tendsto (fun x : ℝ => x⁻¹ * dilog x)
        (𝓝[>] 0) (𝓝 1) := by
    simpa [dilog_zero] using
      RamanujanChallenge.P26.dilog_hasDerivAt_zero26.tendsto_slope_zero_right
  have hlogx :
      Tendsto (fun x : ℝ => Real.log x * x)
        (𝓝[>] 0) (𝓝 0) := by
    simpa [Real.rpow_one] using
      (tendsto_log_mul_rpow_nhdsGT_zero (show (0 : ℝ) < 1 by norm_num))
  have hproductRaw :
      Tendsto
        (fun x : ℝ => (Real.log x * x) * (x⁻¹ * dilog x))
        (𝓝[>] 0) (𝓝 0) := by
    simpa using hlogx.mul hslope
  have hproduct :
      Tendsto (fun x : ℝ => Real.log x * dilog x)
        (𝓝[>] 0) (𝓝 0) := by
    apply Filter.Tendsto.congr' _ hproductRaw
    filter_upwards [self_mem_nhdsWithin] with x hx
    have hxne : x ≠ 0 := ne_of_gt hx
    field_simp [hxne]
  have htriCont :
      ContinuousAt RamanujanChallenge.P26.trilog26 0 :=
    RamanujanChallenge.P26.trilog26_continuousOn_unit.continuousAt
      (Icc_mem_nhds (by norm_num : (-1 : ℝ) < 0)
        (by norm_num : (0 : ℝ) < 1))
  have htri :
      Tendsto RamanujanChallenge.P26.trilog26
        (𝓝[>] 0) (𝓝 0) := by
    simpa using tendsto_nhdsWithin_of_tendsto_nhds htriCont.tendsto
  simpa [minusRadialPrimitive24] using hproduct.neg.add htri

private theorem minusRadialPrimitive24_tendsto_one :
    Tendsto minusRadialPrimitive24 (𝓝[<] (1 : ℝ))
      (𝓝 RamanujanChallenge.P26.zeta3) := by
  have hlog :
      Tendsto Real.log (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    simpa using tendsto_nhdsWithin_of_tendsto_nhds
      (Real.continuousAt_log (by norm_num : (1 : ℝ) ≠ 0)).tendsto
  have hdWithin :
      ContinuousWithinAt dilog (Iio (1 : ℝ)) 1 :=
    (dilog_continuousOn_unit 1 (by norm_num)).mono_of_mem_nhdsWithin
      (Icc_mem_nhdsLT (show (-1 : ℝ) < 1 by norm_num))
  have hd :
      Tendsto dilog (𝓝[<] (1 : ℝ)) (𝓝 (dilog 1)) :=
    hdWithin.tendsto
  have hproduct :
      Tendsto (fun x : ℝ => Real.log x * dilog x)
        (𝓝[<] 1) (𝓝 0) := by
    simpa using hlog.mul hd
  have htWithin :
      ContinuousWithinAt RamanujanChallenge.P26.trilog26
        (Iio (1 : ℝ)) 1 :=
    (RamanujanChallenge.P26.trilog26_continuousOn_unit 1
      (by norm_num)).mono_of_mem_nhdsWithin
      (Icc_mem_nhdsLT (show (-1 : ℝ) < 1 by norm_num))
  have ht :
      Tendsto RamanujanChallenge.P26.trilog26
        (𝓝[<] 1) (𝓝 RamanujanChallenge.P26.zeta3) := by
    simpa [RamanujanChallenge.P26.trilog26_one] using htWithin.tendsto
  simpa [minusRadialPrimitive24] using hproduct.neg.add ht

private theorem minusRadialIntegral24 :
    (∫ x : ℝ in 0..1, minusRadialKernel24 x) = zeta3_24 := by
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_tendsto
    (f := minusRadialPrimitive24)
    (fa := (0 : ℝ)) (fb := RamanujanChallenge.P26.zeta3)
    (by norm_num)
    (fun x hx => minusRadialPrimitive24_hasDerivAt hx.1 hx.2)
    minusRadialKernel24_intervalIntegrable
    minusRadialPrimitive24_tendsto_zero
    minusRadialPrimitive24_tendsto_one]
  unfold RamanujanChallenge.P26.zeta3 zeta3_24
  ring

private theorem reflectedMinusRadialIntegral24 :
    (∫ x : ℝ in 0..1,
      Real.log x * Real.log (1 - x) / (1 - x)) = zeta3_24 := by
  let f : ℝ → ℝ :=
    fun x => Real.log x * Real.log (1 - x) / x
  have hreflect :
      (∫ x : ℝ in 0..1, f (1 - x)) =
        ∫ x : ℝ in 0..1, f x := by
    simpa using
      (intervalIntegral.integral_comp_sub_left
        (a := (0 : ℝ)) (b := 1) f 1)
  have hleft :
      (∫ x : ℝ in 0..1,
        Real.log x * Real.log (1 - x) / (1 - x)) =
        ∫ x : ℝ in 0..1, f (1 - x) := by
    apply intervalIntegral.integral_congr
    intro x _
    dsimp [f]
    ring
  rw [hleft, hreflect]
  simpa [f, minusRadialKernel24] using minusRadialIntegral24

private def minusLogSquareKernel24 (x : ℝ) : ℝ :=
  Real.log (1 - x) ^ 2 / x

private def minusLogSquarePrimitive24 (x : ℝ) : ℝ :=
  Real.log x * Real.log (1 - x) ^ 2 +
    2 * Real.log (1 - x) * dilog (1 - x) -
    2 * RamanujanChallenge.P26.trilog26 (1 - x) +
    2 * RamanujanChallenge.P26.zeta3

private theorem minusLogSquarePrimitive24_hasDerivAt
    {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasDerivAt minusLogSquarePrimitive24
      (minusLogSquareKernel24 x) x := by
  have hxne : x ≠ 0 := ne_of_gt hx0
  have h1x0 : 0 < 1 - x := by linarith
  have h1x1 : 1 - x < 1 := by linarith
  have h1xne : 1 - x ≠ 0 := ne_of_gt h1x0
  have hlogx :
      HasDerivAt (fun y : ℝ => Real.log y) (1 / x) x := by
    simpa [one_div] using Real.hasDerivAt_log hxne
  have hsub :
      HasDerivAt (fun y : ℝ => 1 - y) (-1) x := by
    convert (hasDerivAt_const x 1).sub (hasDerivAt_id x) using 1
    simp
  have hlog1 :
      HasDerivAt (fun y : ℝ => Real.log (1 - y))
        (-1 / (1 - x)) x := by
    convert hsub.log h1xne using 1
  have hdsub :
      HasDerivAt (fun y : ℝ => dilog (1 - y))
        (Real.log x / (1 - x)) x := by
    convert (dilog_hasDerivAt h1x0 h1x1).comp x hsub using 1
    ring
  have htsub :
      HasDerivAt
        (fun y : ℝ => RamanujanChallenge.P26.trilog26 (1 - y))
        (-(dilog (1 - x) / (1 - x))) x := by
    have h :=
      (RamanujanChallenge.P26.trilog26_hasDerivAt_of_abs_lt_one
        (by rw [abs_of_pos h1x0]; exact h1x1) h1xne).comp x hsub
    convert h using 1
    ring
  unfold minusLogSquarePrimitive24
  have htotal :=
    (((hlogx.mul (hlog1.pow 2)).add
      ((hlog1.mul hdsub).const_mul 2)).sub
      (htsub.const_mul 2)).add_const
      (2 * RamanujanChallenge.P26.zeta3)
  convert htotal using 1
  · funext y
    simp only [Pi.add_apply, Pi.sub_apply, Pi.mul_apply, Pi.pow_apply]
    ring
  · unfold minusLogSquareKernel24
    simp only [Pi.pow_apply]
    field_simp [hxne, h1xne]
    ring

private theorem logMulSquareOneMinus24_continuousOn :
    ContinuousOn
      (fun x : ℝ => Real.log x * Real.log (1 - x) ^ 2)
      (Icc (0 : ℝ) (1 / 2)) := by
  have hlog :
      ContinuousOn (fun x : ℝ => Real.log (1 - x))
        (Icc (0 : ℝ) (1 / 2)) := by
    apply (continuousOn_const.sub continuousOn_id).log
    intro x hx
    simpa only [Pi.sub_apply, Pi.one_apply, id_eq] using
      (ne_of_gt (show 0 < 1 - x by linarith [hx.2]))
  have hright :
      ContinuousOn
        (fun x : ℝ =>
          (x * Real.log x) * logOneMinusSlope24 x * Real.log (1 - x))
        (Icc (0 : ℝ) (1 / 2)) :=
    (Real.continuous_mul_log.continuousOn.mul
      logOneMinusSlope24_continuousOn).mul hlog
  apply hright.congr
  intro x _
  by_cases hxzero : x = 0
  · subst x
    simp
  · simp [logOneMinusSlope24, hxzero]
    field_simp [hxzero]

private theorem minusLogSquarePrimitive24_continuousOn_half :
    ContinuousOn minusLogSquarePrimitive24 (Icc (0 : ℝ) (1 / 2)) := by
  have hsub :
      ContinuousOn (fun x : ℝ => 1 - x)
        (Icc (0 : ℝ) (1 / 2)) :=
    continuousOn_const.sub continuousOn_id
  have hsubmem :
      MapsTo (fun x : ℝ => 1 - x)
        (Icc (0 : ℝ) (1 / 2)) (Icc (-1 : ℝ) 1) := by
    intro x hx
    constructor <;> linarith [hx.1, hx.2]
  have hlog :
      ContinuousOn (fun x : ℝ => Real.log (1 - x))
        (Icc (0 : ℝ) (1 / 2)) := by
    apply hsub.log
    intro x hx
    exact ne_of_gt (show 0 < 1 - x by linarith [hx.2])
  have hd :
      ContinuousOn (fun x : ℝ => dilog (1 - x))
        (Icc (0 : ℝ) (1 / 2)) :=
    dilog_continuousOn_unit.comp hsub hsubmem
  have ht :
      ContinuousOn
        (fun x : ℝ => RamanujanChallenge.P26.trilog26 (1 - x))
        (Icc (0 : ℝ) (1 / 2)) :=
    RamanujanChallenge.P26.trilog26_continuousOn_unit.comp hsub hsubmem
  unfold minusLogSquarePrimitive24
  exact
    (((logMulSquareOneMinus24_continuousOn.add
      ((continuousOn_const.mul hlog).mul hd)).sub
      (continuousOn_const.mul ht)).add continuousOn_const)

@[simp] private theorem minusLogSquarePrimitive24_zero :
    minusLogSquarePrimitive24 0 = 0 := by
  simp [minusLogSquarePrimitive24,
    RamanujanChallenge.P26.trilog26_one]

private theorem minusLogSquarePrimitive24_half :
    minusLogSquarePrimitive24 (1 / 2) =
      (1 / 4 : ℝ) * zeta3_24 - (1 / 3 : ℝ) * Real.log 2 ^ 3 := by
  have hloghalf : Real.log (1 / 2 : ℝ) = -Real.log 2 := by
    rw [one_div, Real.log_inv]
  rw [minusLogSquarePrimitive24,
    show 1 - (1 / 2 : ℝ) = 1 / 2 by norm_num,
    hloghalf, RamanujanChallenge.P26.dilog26_half,
    RamanujanChallenge.P26.trilog26_half]
  unfold RamanujanChallenge.P26.zeta3 zeta3_24
  ring

private theorem minusLogSquareKernel24_continuousOn_half :
    ContinuousOn minusLogSquareKernel24 (Icc (0 : ℝ) (1 / 2)) := by
  have haux :
      ContinuousOn
        (fun x : ℝ => x * logOneMinusSlope24 x ^ 2)
        (Icc (0 : ℝ) (1 / 2)) :=
    continuousOn_id.mul (logOneMinusSlope24_continuousOn.pow 2)
  apply haux.congr
  intro x _
  by_cases hxzero : x = 0
  · subst x
    simp [minusLogSquareKernel24]
  · simp [minusLogSquareKernel24, logOneMinusSlope24, hxzero]
    field_simp [hxzero]

private theorem minusLogSquareHalfIntegral24 :
    (∫ x : ℝ in 0..(1 / 2),
      minusLogSquareKernel24 x) =
      (1 / 4 : ℝ) * zeta3_24 - (1 / 3 : ℝ) * Real.log 2 ^ 3 := by
  have hint :
      IntervalIntegrable minusLogSquareKernel24
        MeasureTheory.volume 0 (1 / 2) := by
    apply ContinuousOn.intervalIntegrable
    rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2)]
    exact minusLogSquareKernel24_continuousOn_half
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le
    (f := minusLogSquarePrimitive24)
    (f' := minusLogSquareKernel24)
    (by norm_num)
    minusLogSquarePrimitive24_continuousOn_half
    (fun x hx => minusLogSquarePrimitive24_hasDerivAt hx.1
      (by linarith [hx.2]))
    hint]
  rw [minusLogSquarePrimitive24_half,
    minusLogSquarePrimitive24_zero]
  ring

private def crossPlusTransformedKernel24 (y : ℝ) : ℝ :=
  (Real.log 2 + Real.log y) *
    (Real.log 2 + Real.log (1 - y)) / (1 - y)

private def crossPlusTransformedPrimitive24 (y : ℝ) : ℝ :=
  -(Real.log 2 ^ 2 * Real.log (1 - y)) +
    Real.log 2 * dilog (1 - y) -
    (1 / 2 : ℝ) * Real.log 2 * Real.log (1 - y) ^ 2 -
    (1 / 2 : ℝ) * Real.log y * Real.log (1 - y) ^ 2 +
    (1 / 2 : ℝ) * minusLogSquarePrimitive24 y

private theorem crossPlusTransformedPrimitive24_hasDerivAt
    {y : ℝ} (hy0 : 0 < y) (hyhalf : y < 1 / 2) :
    HasDerivAt crossPlusTransformedPrimitive24
      (crossPlusTransformedKernel24 y) y := by
  have hy1 : y < 1 := by linarith
  have hyne : y ≠ 0 := ne_of_gt hy0
  have h1y0 : 0 < 1 - y := by linarith
  have h1yne : 1 - y ≠ 0 := ne_of_gt h1y0
  have hlogy :
      HasDerivAt (fun x : ℝ => Real.log x) (1 / y) y := by
    simpa [one_div] using Real.hasDerivAt_log hyne
  have hsub :
      HasDerivAt (fun x : ℝ => 1 - x) (-1) y := by
    convert (hasDerivAt_const y 1).sub (hasDerivAt_id y) using 1
    simp
  have hlog1 :
      HasDerivAt (fun x : ℝ => Real.log (1 - x))
        (-1 / (1 - y)) y := by
    convert hsub.log h1yne using 1
  have hdsub :
      HasDerivAt (fun x : ℝ => dilog (1 - x))
        (Real.log y / (1 - y)) y := by
    convert
      (dilog_hasDerivAt h1y0 (by linarith : 1 - y < 1)).comp y hsub using 1
    ring
  have hminus :=
    minusLogSquarePrimitive24_hasDerivAt hy0 hy1
  unfold crossPlusTransformedPrimitive24
  have htotal :=
    (((((hlog1.const_mul (Real.log 2 ^ 2)).neg).add
      (hdsub.const_mul (Real.log 2))).sub
      ((hlog1.pow 2).const_mul ((1 / 2 : ℝ) * Real.log 2))).sub
      ((hlogy.mul (hlog1.pow 2)).const_mul (1 / 2 : ℝ))).add
      (hminus.const_mul (1 / 2 : ℝ))
  convert htotal using 1
  · funext x
    simp only [Pi.add_apply, Pi.sub_apply, Pi.mul_apply, Pi.pow_apply,
      Pi.neg_apply]
    ring
  · unfold crossPlusTransformedKernel24 minusLogSquareKernel24
    simp only [Pi.pow_apply]
    field_simp [hyne, h1yne]
    ring

private theorem crossPlusTransformedPrimitive24_continuousOn :
    ContinuousOn crossPlusTransformedPrimitive24
      (Icc (0 : ℝ) (1 / 2)) := by
  have hsub :
      ContinuousOn (fun y : ℝ => 1 - y)
        (Icc (0 : ℝ) (1 / 2)) :=
    continuousOn_const.sub continuousOn_id
  have hsubmem :
      MapsTo (fun y : ℝ => 1 - y)
        (Icc (0 : ℝ) (1 / 2)) (Icc (-1 : ℝ) 1) := by
    intro y hy
    constructor <;> linarith [hy.1, hy.2]
  have hlog :
      ContinuousOn (fun y : ℝ => Real.log (1 - y))
        (Icc (0 : ℝ) (1 / 2)) := by
    apply hsub.log
    intro y hy
    exact ne_of_gt (show 0 < 1 - y by linarith [hy.2])
  have hd :
      ContinuousOn (fun y : ℝ => dilog (1 - y))
        (Icc (0 : ℝ) (1 / 2)) :=
    dilog_continuousOn_unit.comp hsub hsubmem
  have hfirst :
      ContinuousOn
        (fun y : ℝ => -(Real.log 2 ^ 2 * Real.log (1 - y)))
        (Icc (0 : ℝ) (1 / 2)) :=
    (continuousOn_const.mul hlog).neg
  have hsecond :
      ContinuousOn
        (fun y : ℝ => Real.log 2 * dilog (1 - y))
        (Icc (0 : ℝ) (1 / 2)) :=
    continuousOn_const.mul hd
  have hthird :
      ContinuousOn
        (fun y : ℝ =>
          (1 / 2 : ℝ) * Real.log 2 * Real.log (1 - y) ^ 2)
        (Icc (0 : ℝ) (1 / 2)) :=
    (continuousOn_const.mul continuousOn_const).mul (hlog.pow 2)
  have hfourth :
      ContinuousOn
        (fun y : ℝ =>
          (1 / 2 : ℝ) * (Real.log y * Real.log (1 - y) ^ 2))
        (Icc (0 : ℝ) (1 / 2)) :=
    continuousOn_const.mul logMulSquareOneMinus24_continuousOn
  have hfifth :
      ContinuousOn
        (fun y : ℝ => (1 / 2 : ℝ) * minusLogSquarePrimitive24 y)
        (Icc (0 : ℝ) (1 / 2)) :=
    continuousOn_const.mul
      minusLogSquarePrimitive24_continuousOn_half
  unfold crossPlusTransformedPrimitive24
  simpa only [mul_assoc] using
    (((hfirst.add hsecond).sub hthird).sub hfourth).add hfifth

private theorem crossPlusTransformedPrimitive24_zero :
    crossPlusTransformedPrimitive24 0 =
      Real.log 2 * (Real.pi ^ 2 / 6) := by
  simp [crossPlusTransformedPrimitive24,
    minusLogSquarePrimitive24_zero, dilog_one]

private theorem crossPlusTransformedPrimitive24_half :
    crossPlusTransformedPrimitive24 (1 / 2) =
      Real.log 2 ^ 3 +
        Real.log 2 *
          (Real.pi ^ 2 / 12 - Real.log 2 ^ 2 / 2) +
        (1 / 2 : ℝ) *
          ((1 / 4 : ℝ) * zeta3_24 -
            (1 / 3 : ℝ) * Real.log 2 ^ 3) := by
  have hloghalf : Real.log (1 / 2 : ℝ) = -Real.log 2 := by
    rw [one_div, Real.log_inv]
  rw [crossPlusTransformedPrimitive24,
    show 1 - (1 / 2 : ℝ) = 1 / 2 by norm_num,
    hloghalf, RamanujanChallenge.P26.dilog26_half,
    minusLogSquarePrimitive24_half]
  ring

private def crossPlusKernel24 (x : ℝ) : ℝ :=
  Real.log (1 - x) * Real.log (1 + x) / (1 + x)

private def crossPlusPrimitive24 (x : ℝ) : ℝ :=
  -crossPlusTransformedPrimitive24 ((1 - x) / 2)

private theorem crossPlusPrimitive24_hasDerivAt
    {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasDerivAt crossPlusPrimitive24 (crossPlusKernel24 x) x := by
  let y : ℝ := (1 - x) / 2
  have hy0 : 0 < y := by
    dsimp [y]
    linarith
  have hyhalf : y < 1 / 2 := by
    dsimp [y]
    linarith
  have hyderiv :
      HasDerivAt (fun z : ℝ => (1 - z) / 2) (-1 / 2) x := by
    convert
      ((hasDerivAt_const x 1).sub (hasDerivAt_id x)).div_const 2 using 1
    norm_num
  have hcomp :=
    (crossPlusTransformedPrimitive24_hasDerivAt hy0 hyhalf).comp x hyderiv
  have h1mx0 : 0 < 1 - x := by linarith
  have h1px0 : 0 < 1 + x := by linarith
  have hlogy :
      Real.log y = Real.log (1 - x) - Real.log 2 := by
    dsimp [y]
    rw [Real.log_div (ne_of_gt h1mx0) (by norm_num : (2 : ℝ) ≠ 0)]
  have honeY : 1 - y = (1 + x) / 2 := by
    dsimp [y]
    ring
  have hlogOneY :
      Real.log (1 - y) = Real.log (1 + x) - Real.log 2 := by
    rw [honeY,
      Real.log_div (ne_of_gt h1px0) (by norm_num : (2 : ℝ) ≠ 0)]
  unfold crossPlusPrimitive24
  convert hcomp.neg using 1
  dsimp only [y] at hlogy hlogOneY honeY ⊢
  unfold crossPlusTransformedKernel24 crossPlusKernel24
  rw [hlogy, hlogOneY, honeY]
  field_simp [ne_of_gt h1px0]
  ring

private theorem crossPlusPrimitive24_continuousOn :
    ContinuousOn crossPlusPrimitive24 (Icc (0 : ℝ) 1) := by
  have hy :
      ContinuousOn (fun x : ℝ => (1 - x) / 2)
        (Icc (0 : ℝ) 1) := by
    fun_prop
  have hymem :
      MapsTo (fun x : ℝ => (1 - x) / 2)
        (Icc (0 : ℝ) 1) (Icc (0 : ℝ) (1 / 2)) := by
    intro x hx
    constructor <;> linarith [hx.1, hx.2]
  unfold crossPlusPrimitive24
  exact
    (crossPlusTransformedPrimitive24_continuousOn.comp hy hymem).neg

private theorem crossPlusKernel24_intervalIntegrable :
    IntervalIntegrable crossPlusKernel24
      MeasureTheory.volume 0 1 := by
  have hlogSub :
      IntervalIntegrable (fun x : ℝ => Real.log (1 - x))
        MeasureTheory.volume 0 1 := by
    have h :=
      (intervalIntegral.intervalIntegrable_log'
        (a := (0 : ℝ)) (b := 1)).comp_sub_left 1
    simpa using h.symm
  have hfactor :
      ContinuousOn (fun x : ℝ => Real.log (1 + x) / (1 + x))
        (Set.uIcc (0 : ℝ) 1) := by
    have hfactorIcc :
        ContinuousOn (fun x : ℝ => Real.log (1 + x) / (1 + x))
          (Icc (0 : ℝ) 1) := by
      have hadd :
          ContinuousOn (fun x : ℝ => 1 + x) (Icc (0 : ℝ) 1) := by
        fun_prop
      apply ContinuousOn.div
      · apply hadd.log
        intro x hx
        exact ne_of_gt (show 0 < 1 + x by linarith [hx.1])
      · exact hadd
      · intro x hx
        exact ne_of_gt (show 0 < 1 + x by linarith [hx.1])
    rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
    exact hfactorIcc
  have hint := hlogSub.continuousOn_mul hfactor
  apply IntervalIntegrable.congr
    (f := fun x : ℝ =>
      (Real.log (1 + x) / (1 + x)) * Real.log (1 - x)) ?_ hint
  intro x _
  unfold crossPlusKernel24
  ring

private theorem crossPlusIntegral24 :
    (∫ x : ℝ in 0..1, crossPlusKernel24 x) =
      (1 / 8 : ℝ) * zeta3_24 -
        (1 / 2 : ℝ) * Real.log 2 * (Real.pi ^ 2 / 6) +
        (1 / 3 : ℝ) * Real.log 2 ^ 3 := by
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le
    (f := crossPlusPrimitive24)
    (f' := crossPlusKernel24)
    (by norm_num)
    crossPlusPrimitive24_continuousOn
    (fun x hx => crossPlusPrimitive24_hasDerivAt hx.1 hx.2)
    crossPlusKernel24_intervalIntegrable]
  simp only [crossPlusPrimitive24]
  norm_num
  rw [crossPlusTransformedPrimitive24_zero,
    crossPlusTransformedPrimitive24_half]
  ring

private theorem oneSub_tendsto_nhdsLT_one_nhdsGT_zero24 :
    Tendsto (fun x : ℝ => 1 - x)
      (𝓝[<] (1 : ℝ)) (𝓝[>] (0 : ℝ)) := by
  rw [tendsto_nhdsWithin_iff]
  refine ⟨?_, ?_⟩
  · have hcont : ContinuousAt (fun x : ℝ => 1 - x) 1 := by
      fun_prop
    simpa using hcont.tendsto.mono_left
      (show (𝓝[<] (1 : ℝ)) ≤ 𝓝 1 from inf_le_left)
  · filter_upwards [self_mem_nhdsWithin] with x hx
    exact sub_pos.mpr (show x < 1 by simpa only [mem_Iio] using hx)

private theorem tendsto_self_mul_log_sq_nhdsGT_zero24 :
    Tendsto (fun t : ℝ => t * Real.log t ^ 2)
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  have hhalf :=
    tendsto_log_mul_rpow_nhdsGT_zero
      (show (0 : ℝ) < 1 / 2 by norm_num)
  have hsquare := hhalf.mul hhalf
  have hsquare' :
      Tendsto
        (fun t : ℝ =>
          (Real.log t * t ^ (1 / 2 : ℝ)) *
            (Real.log t * t ^ (1 / 2 : ℝ)))
        (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    simpa using hsquare
  apply Filter.Tendsto.congr' _ hsquare'
  filter_upwards [self_mem_nhdsWithin] with t ht
  have hrpow :
      t ^ (1 / 2 : ℝ) * t ^ (1 / 2 : ℝ) = t := by
    rw [← Real.rpow_add ht]
    norm_num
  rw [show
      (Real.log t * t ^ (1 / 2 : ℝ)) *
          (Real.log t * t ^ (1 / 2 : ℝ)) =
        Real.log t ^ 2 *
          (t ^ (1 / 2 : ℝ) * t ^ (1 / 2 : ℝ)) by ring,
    hrpow]
  ring

private theorem minusLogSquarePrimitive24_tendsto_one :
    Tendsto minusLogSquarePrimitive24 (𝓝[<] (1 : ℝ))
      (𝓝 (2 * RamanujanChallenge.P26.zeta3)) := by
  have hlogMinusDeriv :
      HasDerivAt (fun t : ℝ => Real.log (1 - t)) (-1) 0 := by
    have hinner :
        HasDerivAt (fun t : ℝ => 1 - t) (-1) 0 := by
      convert (hasDerivAt_const (0 : ℝ) 1).sub (hasDerivAt_id 0) using 1
      norm_num
    simpa using hinner.log
      (by norm_num : 1 - (0 : ℝ) ≠ 0)
  have hslope0 :
      Tendsto
        (fun t : ℝ => Real.log (1 - t) / t)
        (𝓝[>] (0 : ℝ)) (𝓝 (-1)) := by
    have h := hlogMinusDeriv.tendsto_slope_zero_right
    convert h using 1
    funext t
    simp [smul_eq_mul]
    ring
  have hslopeX :
      Tendsto (fun x : ℝ => Real.log x / (1 - x))
        (𝓝[<] (1 : ℝ)) (𝓝 (-1)) := by
    have hcomp :=
      hslope0.comp oneSub_tendsto_nhdsLT_one_nhdsGT_zero24
    apply Filter.Tendsto.congr' _ hcomp
    filter_upwards [self_mem_nhdsWithin] with x hx
    have hxlt : x < 1 := by simpa only [mem_Iio] using hx
    have hne : 1 - x ≠ 0 := ne_of_gt (sub_pos.mpr hxlt)
    simp only [Function.comp_apply]
    rw [show 1 - (1 - x) = x by ring]
  have hsmall :
      Tendsto
        (fun x : ℝ => (1 - x) * Real.log (1 - x) ^ 2)
        (𝓝[<] (1 : ℝ)) (𝓝 0) :=
    tendsto_self_mul_log_sq_nhdsGT_zero24.comp
      oneSub_tendsto_nhdsLT_one_nhdsGT_zero24
  have hfirstRaw :
      Tendsto
        (fun x : ℝ =>
          (Real.log x / (1 - x)) *
            ((1 - x) * Real.log (1 - x) ^ 2))
        (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    simpa using hslopeX.mul hsmall
  have hfirst :
      Tendsto
        (fun x : ℝ => Real.log x * Real.log (1 - x) ^ 2)
        (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    apply Filter.Tendsto.congr' _ hfirstRaw
    filter_upwards [self_mem_nhdsWithin] with x hx
    have hxlt : x < 1 := by simpa only [mem_Iio] using hx
    have hne : 1 - x ≠ 0 := ne_of_gt (sub_pos.mpr hxlt)
    field_simp [hne]
  have hlogt :
      Tendsto (fun t : ℝ => Real.log t * t)
        (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    simpa [Real.rpow_one] using
      (tendsto_log_mul_rpow_nhdsGT_zero
        (show (0 : ℝ) < 1 by norm_num))
  have hdSlope :
      Tendsto (fun t : ℝ => t⁻¹ * dilog t)
        (𝓝[>] (0 : ℝ)) (𝓝 1) := by
    simpa [dilog_zero] using
      RamanujanChallenge.P26.dilog_hasDerivAt_zero26.tendsto_slope_zero_right
  have hlogDBase :
      Tendsto (fun t : ℝ => Real.log t * dilog t)
        (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    have hraw :=
      hlogt.mul hdSlope
    have hraw' :
        Tendsto
          (fun t : ℝ =>
            (Real.log t * t) * (t⁻¹ * dilog t))
          (𝓝[>] (0 : ℝ)) (𝓝 0) := by
      simpa using hraw
    apply Filter.Tendsto.congr' _ hraw'
    filter_upwards [self_mem_nhdsWithin] with t ht
    have htne : t ≠ 0 := ne_of_gt ht
    field_simp [htne]
  have hlogD :
      Tendsto
        (fun x : ℝ => Real.log (1 - x) * dilog (1 - x))
        (𝓝[<] (1 : ℝ)) (𝓝 0) :=
    hlogDBase.comp oneSub_tendsto_nhdsLT_one_nhdsGT_zero24
  have htri0 :
      ContinuousAt RamanujanChallenge.P26.trilog26 0 :=
    RamanujanChallenge.P26.trilog26_continuousOn_unit.continuousAt
      (Icc_mem_nhds (by norm_num : (-1 : ℝ) < 0)
        (by norm_num : (0 : ℝ) < 1))
  have htri :
      Tendsto
        (fun x : ℝ => RamanujanChallenge.P26.trilog26 (1 - x))
        (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    have hsub :
        Tendsto (fun x : ℝ => 1 - x)
          (𝓝[<] (1 : ℝ)) (𝓝 0) := by
      have hcont : ContinuousAt (fun x : ℝ => 1 - x) 1 := by
        fun_prop
      simpa using hcont.tendsto.mono_left
        (show (𝓝[<] (1 : ℝ)) ≤ 𝓝 1 from inf_le_left)
    simpa using htri0.tendsto.comp hsub
  unfold minusLogSquarePrimitive24
  convert
    (((hfirst.add (hlogD.const_mul 2)).sub
      (htri.const_mul 2)).add_const
      (2 * RamanujanChallenge.P26.zeta3)) using 1
  · funext x
    ring
  · ring

private def plusLogSquareKernel24 (x : ℝ) : ℝ :=
  Real.log (1 + x) ^ 2 / x

private theorem plusLogSquareKernel24_continuousOn :
    ContinuousOn plusLogSquareKernel24 (Icc (0 : ℝ) 1) := by
  have haux :
      ContinuousOn
        (fun x : ℝ =>
          x * RamanujanChallenge.P26.logOnePlusSlope26 x ^ 2)
        (Icc (0 : ℝ) 1) :=
    continuousOn_id.mul
      (RamanujanChallenge.P26.logOnePlusSlope26_continuousOn.pow 2)
  apply haux.congr
  intro x _
  by_cases hxzero : x = 0
  · subst x
    simp [plusLogSquareKernel24]
  · simp [plusLogSquareKernel24,
      RamanujanChallenge.P26.logOnePlusSlope26, hxzero]
    field_simp [hxzero]

private theorem plusLogSquareKernel24_intervalIntegrable :
    IntervalIntegrable plusLogSquareKernel24
      MeasureTheory.volume 0 1 := by
  apply ContinuousOn.intervalIntegrable
  rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
  exact plusLogSquareKernel24_continuousOn

private def plusLogSquarePrimitive24 (x : ℝ) : ℝ :=
  ∫ t : ℝ in 0..x, plusLogSquareKernel24 t

private theorem plusLogSquarePrimitive24_continuousOn :
    ContinuousOn plusLogSquarePrimitive24 (Icc (0 : ℝ) 1) := by
  have hac :
      AbsolutelyContinuousOnInterval
        (fun x : ℝ => ∫ v : ℝ in 0..x, plusLogSquareKernel24 v)
        0 1 :=
    plusLogSquareKernel24_intervalIntegrable
      |>.absolutelyContinuousOnInterval_intervalIntegral
        (c := (0 : ℝ)) (by
          norm_num [Set.uIcc_of_le])
  simpa [plusLogSquarePrimitive24,
    Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hac.continuousOn

private theorem plusLogSquarePrimitive24_hasDerivAt
    {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasDerivAt plusLogSquarePrimitive24
      (plusLogSquareKernel24 x) x := by
  have hcontAt :
      ContinuousAt plusLogSquareKernel24 x :=
    plusLogSquareKernel24_continuousOn.continuousAt
      (Icc_mem_nhds hx0 hx1)
  have hint :
      IntervalIntegrable plusLogSquareKernel24
        MeasureTheory.volume 0 x := by
    apply ContinuousOn.intervalIntegrable
    rw [Set.uIcc_of_le hx0.le]
    exact plusLogSquareKernel24_continuousOn.mono (by
      intro y hy
      exact ⟨hy.1, hy.2.trans hx1.le⟩)
  unfold plusLogSquarePrimitive24
  exact intervalIntegral.integral_hasDerivAt_right hint
    (ContinuousAt.stronglyMeasurableAtFilter isOpen_Ioi
      (fun y (hy : y ∈ Ioi (0 : ℝ)) => by
        have hy0 : 0 < y := hy
        have hyne : y ≠ 0 := ne_of_gt hy0
        have h1yne : 1 + y ≠ 0 := by linarith
        have hinner :
            ContinuousAt (fun z : ℝ => 1 + z) y := by
          fun_prop
        have hlog :
            ContinuousAt (fun z : ℝ => Real.log (1 + z)) y :=
          (Real.continuousAt_log h1yne).comp hinner
        unfold plusLogSquareKernel24
        exact (hlog.pow 2).div continuousAt_id hyne)
      x hx0)
    hcontAt

@[simp] private theorem plusLogSquarePrimitive24_zero :
    plusLogSquarePrimitive24 0 = 0 := by
  simp [plusLogSquarePrimitive24]

private theorem plusLogSquarePrimitive24_one :
    plusLogSquarePrimitive24 1 = (1 / 4 : ℝ) * zeta3_24 := by
  unfold plusLogSquarePrimitive24 plusLogSquareKernel24
  rw [RamanujanChallenge.P26.logSquareOnePlusIntegral26]
  unfold RamanujanChallenge.P26.zeta3 zeta3_24
  rfl

private def crossRadialKernel24 (x : ℝ) : ℝ :=
  Real.log (1 - x) * Real.log (1 + x) / x

private def crossRadialPrimitive24 (x : ℝ) : ℝ :=
  (1 / 4 : ℝ) * minusLogSquarePrimitive24 (x ^ 2) -
    (1 / 2 : ℝ) * minusLogSquarePrimitive24 x -
    (1 / 2 : ℝ) * plusLogSquarePrimitive24 x

private theorem crossRadialPrimitive24_hasDerivAt
    {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasDerivAt crossRadialPrimitive24
      (crossRadialKernel24 x) x := by
  have hxne : x ≠ 0 := ne_of_gt hx0
  have hx2pos : 0 < x ^ 2 := sq_pos_of_pos hx0
  have hx2lt : x ^ 2 < 1 := by nlinarith
  have hsq :
      HasDerivAt (fun y : ℝ => y * y) (2 * x) x := by
    convert (hasDerivAt_id x).mul (hasDerivAt_id x) using 1
    simp [id_eq]
    ring
  have houter :
      HasDerivAt minusLogSquarePrimitive24
        (minusLogSquareKernel24 (x * x)) (x * x) := by
    simpa only [pow_two] using
      minusLogSquarePrimitive24_hasDerivAt hx2pos hx2lt
  have hminusSq :
      HasDerivAt
        (fun y : ℝ => minusLogSquarePrimitive24 (y * y))
        (minusLogSquareKernel24 (x ^ 2) * (2 * x)) x := by
    have hcomp :=
      houter.comp_of_eq x hsq
        (by rfl : x * x = (fun y : ℝ => y * y) x)
    simpa only [Function.comp_def, pow_two] using hcomp
  have hminus :=
    minusLogSquarePrimitive24_hasDerivAt hx0 hx1
  have hplus :=
    plusLogSquarePrimitive24_hasDerivAt hx0 hx1
  have h1mx0 : 0 < 1 - x := by linarith
  have h1px0 : 0 < 1 + x := by linarith
  have hlogmul :
      Real.log (1 - x ^ 2) =
        Real.log (1 - x) + Real.log (1 + x) := by
    rw [show 1 - x ^ 2 = (1 - x) * (1 + x) by ring,
      Real.log_mul (ne_of_gt h1mx0) (ne_of_gt h1px0)]
  unfold crossRadialPrimitive24
  convert
    ((hminusSq.const_mul (1 / 4 : ℝ)).sub
      (hminus.const_mul (1 / 2 : ℝ))).sub
      (hplus.const_mul (1 / 2 : ℝ)) using 1
  · funext y
    simp only [crossRadialPrimitive24, pow_two, Pi.sub_apply,
      Pi.mul_apply]
  · unfold minusLogSquareKernel24 plusLogSquareKernel24
      crossRadialKernel24
    rw [hlogmul]
    field_simp [hxne]
    ring

private theorem crossRadialKernel24_intervalIntegrable :
    IntervalIntegrable crossRadialKernel24
      MeasureTheory.volume 0 1 := by
  apply IntervalIntegrable.trans (b := (1 / 2 : ℝ))
  · have hcont :
        ContinuousOn
          (fun x : ℝ =>
            x * logOneMinusSlope24 x *
              RamanujanChallenge.P26.logOnePlusSlope26 x)
          (Icc (0 : ℝ) (1 / 2)) := by
      exact
        (continuousOn_id.mul logOneMinusSlope24_continuousOn).mul
          (RamanujanChallenge.P26.logOnePlusSlope26_continuousOn.mono
            (by intro x hx; exact ⟨hx.1, hx.2.trans (by norm_num)⟩))
    have hint :
        IntervalIntegrable
          (fun x : ℝ =>
            x * logOneMinusSlope24 x *
              RamanujanChallenge.P26.logOnePlusSlope26 x)
          MeasureTheory.volume 0 (1 / 2) := by
      apply ContinuousOn.intervalIntegrable
      rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2)]
      exact hcont
    apply IntervalIntegrable.congr
      (f := fun x : ℝ =>
        x * logOneMinusSlope24 x *
          RamanujanChallenge.P26.logOnePlusSlope26 x) ?_ hint
    intro x hx
    have hx' : x ∈ Ioc (0 : ℝ) (1 / 2) := by
      simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2)] using hx
    have hxne : x ≠ 0 := ne_of_gt hx'.1
    simp [crossRadialKernel24, logOneMinusSlope24,
      RamanujanChallenge.P26.logOnePlusSlope26, hxne]
    field_simp [hxne]
  · have hlogSub :
        IntervalIntegrable (fun x : ℝ => Real.log (1 - x))
          MeasureTheory.volume (1 / 2) 1 := by
      have h :=
        (intervalIntegral.intervalIntegrable_log'
          (a := (0 : ℝ)) (b := (1 / 2 : ℝ))).comp_sub_left 1
      convert h.symm using 1 <;> norm_num
    have hfactor :
        ContinuousOn (fun x : ℝ => Real.log (1 + x) / x)
          (Set.uIcc (1 / 2 : ℝ) 1) := by
      have hfactorIcc :
          ContinuousOn (fun x : ℝ => Real.log (1 + x) / x)
            (Icc (1 / 2 : ℝ) 1) := by
        apply ContinuousOn.div
        · apply (continuousOn_const.add continuousOn_id).log
          intro x hx
          exact ne_of_gt (show 0 < 1 + x by linarith [hx.1])
        · exact continuousOn_id
        · intro x hx
          exact ne_of_gt (show 0 < x by linarith [hx.1])
      rw [Set.uIcc_of_le (by norm_num : (1 / 2 : ℝ) ≤ 1)]
      exact hfactorIcc
    have hint := hlogSub.continuousOn_mul hfactor
    apply IntervalIntegrable.congr
      (f := fun x : ℝ =>
        (Real.log (1 + x) / x) * Real.log (1 - x)) ?_ hint
    intro x _
    unfold crossRadialKernel24
    ring

private theorem crossRadialPrimitive24_continuousOn_half :
    ContinuousOn crossRadialPrimitive24 (Icc (0 : ℝ) (1 / 2)) := by
  have hsq :
      ContinuousOn (fun x : ℝ => x ^ 2)
        (Icc (0 : ℝ) (1 / 2)) := by
    fun_prop
  have hsqmem :
      MapsTo (fun x : ℝ => x ^ 2)
        (Icc (0 : ℝ) (1 / 2)) (Icc (0 : ℝ) (1 / 2)) := by
    intro x hx
    constructor
    · positivity
    · nlinarith [mul_nonneg hx.1 (sub_nonneg.mpr hx.2)]
  unfold crossRadialPrimitive24
  exact
    ((continuousOn_const.mul
      (minusLogSquarePrimitive24_continuousOn_half.comp hsq hsqmem)).sub
      (continuousOn_const.mul
        minusLogSquarePrimitive24_continuousOn_half)).sub
      (continuousOn_const.mul
        (plusLogSquarePrimitive24_continuousOn.mono
          (by intro x hx; exact ⟨hx.1, hx.2.trans (by norm_num)⟩)))

@[simp] private theorem crossRadialPrimitive24_zero :
    crossRadialPrimitive24 0 = 0 := by
  simp [crossRadialPrimitive24]

private theorem crossRadialPrimitive24_tendsto_zero :
    Tendsto crossRadialPrimitive24
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  have hwithin :
      ContinuousWithinAt crossRadialPrimitive24
        (Ioi (0 : ℝ)) 0 :=
    (crossRadialPrimitive24_continuousOn_half 0
      (by norm_num)).mono_of_mem_nhdsWithin
      (Icc_mem_nhdsGT (by norm_num : (0 : ℝ) < 1 / 2))
  simpa [crossRadialPrimitive24_zero] using hwithin.tendsto

private theorem square_tendsto_nhdsLT_one24 :
    Tendsto (fun x : ℝ => x ^ 2)
      (𝓝[<] (1 : ℝ)) (𝓝[<] (1 : ℝ)) := by
  rw [tendsto_nhdsWithin_iff]
  refine ⟨?_, ?_⟩
  · have hcont : ContinuousAt (fun x : ℝ => x ^ 2) 1 := by
      fun_prop
    simpa using hcont.tendsto.mono_left
      (show (𝓝[<] (1 : ℝ)) ≤ 𝓝 1 from inf_le_left)
  · filter_upwards [Ioo_mem_nhdsLT (show (0 : ℝ) < 1 by norm_num)]
      with x hx
    have : x ^ 2 < 1 := by
      nlinarith [mul_pos hx.1 (sub_pos.mpr hx.2)]
    simpa only [mem_Iio] using this

private theorem crossRadialPrimitive24_tendsto_one :
    Tendsto crossRadialPrimitive24 (𝓝[<] (1 : ℝ))
      (𝓝 (-(5 : ℝ) / 8 * zeta3_24)) := by
  have hminusSq :=
    minusLogSquarePrimitive24_tendsto_one.comp
      square_tendsto_nhdsLT_one24
  have hminus := minusLogSquarePrimitive24_tendsto_one
  have hplusWithin :
      ContinuousWithinAt plusLogSquarePrimitive24
        (Iio (1 : ℝ)) 1 :=
    (plusLogSquarePrimitive24_continuousOn 1
      (by norm_num)).mono_of_mem_nhdsWithin
      (Icc_mem_nhdsLT (show (0 : ℝ) < 1 by norm_num))
  have hplus :
      Tendsto plusLogSquarePrimitive24
        (𝓝[<] (1 : ℝ)) (𝓝 ((1 / 4 : ℝ) * zeta3_24)) := by
    simpa [plusLogSquarePrimitive24_one] using hplusWithin.tendsto
  unfold crossRadialPrimitive24
  convert
    ((hminusSq.const_mul (1 / 4 : ℝ)).sub
      (hminus.const_mul (1 / 2 : ℝ))).sub
      (hplus.const_mul (1 / 2 : ℝ)) using 1
  unfold RamanujanChallenge.P26.zeta3 zeta3_24
  ring

private theorem crossRadialIntegral24 :
    (∫ x : ℝ in 0..1, crossRadialKernel24 x) =
      -(5 : ℝ) / 8 * zeta3_24 := by
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_tendsto
    (f := crossRadialPrimitive24)
    (fa := (0 : ℝ)) (fb := -(5 : ℝ) / 8 * zeta3_24)
    (by norm_num)
    (fun x hx => crossRadialPrimitive24_hasDerivAt hx.1 hx.2)
    crossRadialKernel24_intervalIntegrable
    crossRadialPrimitive24_tendsto_zero
    crossRadialPrimitive24_tendsto_one]
  ring

private def minusSimpleKernel24 (x : ℝ) : ℝ :=
  Real.log (1 - x) / x

private theorem minusSimpleKernel24_intervalIntegrable :
    IntervalIntegrable minusSimpleKernel24
      MeasureTheory.volume 0 1 := by
  apply IntervalIntegrable.trans (b := (1 / 2 : ℝ))
  · have hint :
        IntervalIntegrable logOneMinusSlope24
          MeasureTheory.volume 0 (1 / 2) := by
      apply ContinuousOn.intervalIntegrable
      rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2)]
      exact logOneMinusSlope24_continuousOn
    apply IntervalIntegrable.congr
      (f := logOneMinusSlope24) ?_ hint
    intro x hx
    have hx' : x ∈ Ioc (0 : ℝ) (1 / 2) := by
      simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2)] using hx
    have hxne : x ≠ 0 := ne_of_gt hx'.1
    simp [minusSimpleKernel24, logOneMinusSlope24, hxne]
  · have hlogSub :
        IntervalIntegrable (fun x : ℝ => Real.log (1 - x))
          MeasureTheory.volume (1 / 2) 1 := by
      have h :=
        (intervalIntegral.intervalIntegrable_log'
          (a := (0 : ℝ)) (b := (1 / 2 : ℝ))).comp_sub_left 1
      convert h.symm using 1 <;> norm_num
    have hfactor :
        ContinuousOn (fun x : ℝ => 1 / x)
          (Set.uIcc (1 / 2 : ℝ) 1) := by
      have hfactorIcc :
          ContinuousOn (fun x : ℝ => 1 / x)
            (Icc (1 / 2 : ℝ) 1) := by
        apply continuousOn_const.div continuousOn_id
        intro x hx
        exact ne_of_gt (show 0 < x by linarith [hx.1])
      rw [Set.uIcc_of_le (by norm_num : (1 / 2 : ℝ) ≤ 1)]
      exact hfactorIcc
    have hint := hlogSub.continuousOn_mul hfactor
    apply IntervalIntegrable.congr
      (f := fun x : ℝ => (1 / x) * Real.log (1 - x)) ?_ hint
    intro x _
    unfold minusSimpleKernel24
    ring

private theorem minusSimpleIntegral24 :
    (∫ x : ℝ in 0..1, minusSimpleKernel24 x) =
      -(Real.pi ^ 2 / 6) := by
  let F : ℝ → ℝ := fun x => -dilog x
  have hcont : ContinuousOn F (Icc (0 : ℝ) 1) :=
    dilog_continuousOn_unit.mono (by
      intro x hx
      constructor <;> linarith [hx.1, hx.2]) |>.neg
  have hderiv :
      ∀ x ∈ Ioo (0 : ℝ) 1,
        HasDerivAt F (minusSimpleKernel24 x) x := by
    intro x hx
    have hd := dilog_hasDerivAt hx.1 hx.2
    unfold F minusSimpleKernel24
    convert hd.neg using 1
    ring
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le
    (f := F) (f' := minusSimpleKernel24)
    (by norm_num) hcont hderiv
    minusSimpleKernel24_intervalIntegrable]
  simp [F, dilog_one, dilog_zero]

private def plusSimpleKernel24 (x : ℝ) : ℝ :=
  Real.log (1 + x) / x

private theorem plusSimpleKernel24_intervalIntegrable :
    IntervalIntegrable plusSimpleKernel24
      MeasureTheory.volume 0 1 := by
  have hint :
      IntervalIntegrable RamanujanChallenge.P26.logOnePlusSlope26
        MeasureTheory.volume 0 1 := by
    apply ContinuousOn.intervalIntegrable
    rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
    exact RamanujanChallenge.P26.logOnePlusSlope26_continuousOn
  apply IntervalIntegrable.congr
    (f := RamanujanChallenge.P26.logOnePlusSlope26) ?_ hint
  intro x hx
  have hx' : x ∈ Ioc (0 : ℝ) 1 := by
    simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hx
  have hxne : x ≠ 0 := ne_of_gt hx'.1
  simp [plusSimpleKernel24,
    RamanujanChallenge.P26.logOnePlusSlope26, hxne]

private theorem plusSimpleIntegral24 :
    (∫ x : ℝ in 0..1, plusSimpleKernel24 x) =
      Real.pi ^ 2 / 12 := by
  let F : ℝ → ℝ := fun x => -dilog (-x)
  have hcont :
      ContinuousOn F (Icc (0 : ℝ) 1) := by
    unfold F
    apply ContinuousOn.neg
    apply dilog_continuousOn_unit.comp
    · fun_prop
    · intro x hx
      constructor <;> linarith [hx.1, hx.2]
  have hderiv :
      ∀ x ∈ Ioo (0 : ℝ) 1,
        HasDerivAt F (plusSimpleKernel24 x) x := by
    intro x hx
    have hxne : x ≠ 0 := ne_of_gt hx.1
    have hxabs : |x| < 1 := by
      rw [abs_of_pos hx.1]
      exact hx.2
    have hdRaw :=
      (dilog_hasDerivAt_of_abs_lt_one
        (show |-x| < 1 by simpa using hxabs)
        (neg_ne_zero.mpr hxne)).comp x (hasDerivAt_neg x)
    unfold F plusSimpleKernel24
    convert hdRaw.neg using 1
    rw [show 1 - -x = 1 + x by ring]
    field_simp [hxne]
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le
    (f := F) (f' := plusSimpleKernel24)
    (by norm_num) hcont hderiv
    plusSimpleKernel24_intervalIntegrable]
  simp [F, RamanujanChallenge.P26.dilog_neg_one26, dilog_zero]
  ring

private def minusDenomKernel24 (x : ℝ) : ℝ :=
  Real.log (1 - x) / (1 + x)

private def minusDenomPrimitive24 (x : ℝ) : ℝ :=
  Real.log 2 * Real.log (1 + x) - dilog ((1 + x) / 2)

private theorem minusDenomPrimitive24_hasDerivAt
    {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasDerivAt minusDenomPrimitive24 (minusDenomKernel24 x) x := by
  have h1px0 : 0 < 1 + x := by linarith
  have h1pxne : 1 + x ≠ 0 := ne_of_gt h1px0
  have h1mx0 : 0 < 1 - x := by linarith
  let z : ℝ := (1 + x) / 2
  have hz0 : 0 < z := by
    dsimp [z]
    positivity
  have hz1 : z < 1 := by
    dsimp [z]
    linarith
  have hinner :
      HasDerivAt (fun y : ℝ => 1 + y) 1 x := by
    convert (hasDerivAt_const x 1).add (hasDerivAt_id x) using 1
    simp
  have hlog :
      HasDerivAt (fun y : ℝ => Real.log (1 + y))
        (1 / (1 + x)) x := by
    convert hinner.log h1pxne using 1
  have hzderiv :
      HasDerivAt (fun y : ℝ => (1 + y) / 2) (1 / 2) x := by
    convert hinner.div_const 2 using 1
  have hd :
      HasDerivAt (fun y : ℝ => dilog ((1 + y) / 2))
        ((-Real.log (1 - z) / z) * (1 / 2)) x := by
    simpa [z] using (dilog_hasDerivAt hz0 hz1).comp x hzderiv
  have honeZ : 1 - z = (1 - x) / 2 := by
    dsimp [z]
    ring
  have hlogOneZ :
      Real.log (1 - z) = Real.log (1 - x) - Real.log 2 := by
    rw [honeZ,
      Real.log_div (ne_of_gt h1mx0) (by norm_num : (2 : ℝ) ≠ 0)]
  unfold minusDenomPrimitive24
  convert (hlog.const_mul (Real.log 2)).sub hd using 1
  unfold minusDenomKernel24
  rw [hlogOneZ]
  dsimp [z]
  field_simp [h1pxne]
  ring

private theorem minusDenomPrimitive24_continuousOn :
    ContinuousOn minusDenomPrimitive24 (Icc (0 : ℝ) 1) := by
  have hadd :
      ContinuousOn (fun x : ℝ => 1 + x) (Icc (0 : ℝ) 1) := by
    fun_prop
  have hlog :
      ContinuousOn (fun x : ℝ => Real.log (1 + x))
        (Icc (0 : ℝ) 1) := by
    apply hadd.log
    intro x hx
    exact ne_of_gt (show 0 < 1 + x by linarith [hx.1])
  have hz :
      ContinuousOn (fun x : ℝ => (1 + x) / 2)
        (Icc (0 : ℝ) 1) := by
    fun_prop
  have hzmem :
      MapsTo (fun x : ℝ => (1 + x) / 2)
        (Icc (0 : ℝ) 1) (Icc (-1 : ℝ) 1) := by
    intro x hx
    constructor <;> linarith [hx.1, hx.2]
  unfold minusDenomPrimitive24
  exact
    (continuousOn_const.mul hlog).sub
      (dilog_continuousOn_unit.comp hz hzmem)

private theorem minusDenomKernel24_intervalIntegrable :
    IntervalIntegrable minusDenomKernel24
      MeasureTheory.volume 0 1 := by
  have hlogSub :
      IntervalIntegrable (fun x : ℝ => Real.log (1 - x))
        MeasureTheory.volume 0 1 := by
    have h :=
      (intervalIntegral.intervalIntegrable_log'
        (a := (0 : ℝ)) (b := 1)).comp_sub_left 1
    simpa using h.symm
  have hfactor :
      ContinuousOn (fun x : ℝ => 1 / (1 + x))
        (Set.uIcc (0 : ℝ) 1) := by
    have hfactorIcc :
        ContinuousOn (fun x : ℝ => 1 / (1 + x))
          (Icc (0 : ℝ) 1) := by
      apply continuousOn_const.div (continuousOn_const.add continuousOn_id)
      intro x hx
      exact ne_of_gt (show 0 < 1 + x by linarith [hx.1])
    rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
    exact hfactorIcc
  have hint := hlogSub.continuousOn_mul hfactor
  apply IntervalIntegrable.congr
    (f := fun x : ℝ => (1 / (1 + x)) * Real.log (1 - x)) ?_ hint
  intro x _
  unfold minusDenomKernel24
  ring

private theorem minusDenomIntegral24 :
    (∫ x : ℝ in 0..1, minusDenomKernel24 x) =
      -(1 / 2 : ℝ) * (Real.pi ^ 2 / 6) +
        (1 / 2 : ℝ) * Real.log 2 ^ 2 := by
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le
    (f := minusDenomPrimitive24)
    (f' := minusDenomKernel24)
    (by norm_num) minusDenomPrimitive24_continuousOn
    (fun x hx => minusDenomPrimitive24_hasDerivAt hx.1 hx.2)
    minusDenomKernel24_intervalIntegrable]
  unfold minusDenomPrimitive24
  norm_num
  rw [dilog_one, RamanujanChallenge.P26.dilog26_half]
  ring

private def plusDenomKernel24 (x : ℝ) : ℝ :=
  Real.log (1 + x) / (1 + x)

private def plusDenomPrimitive24 (x : ℝ) : ℝ :=
  (1 / 2 : ℝ) * Real.log (1 + x) ^ 2

private theorem plusDenomIntegral24 :
    (∫ x : ℝ in 0..1, plusDenomKernel24 x) =
      (1 / 2 : ℝ) * Real.log 2 ^ 2 := by
  have hcont :
      ContinuousOn plusDenomPrimitive24 (Icc (0 : ℝ) 1) := by
    unfold plusDenomPrimitive24
    apply continuousOn_const.mul
    apply ((continuousOn_const.add continuousOn_id).log
      (fun x hx => ne_of_gt
        (show 0 < 1 + x by linarith [hx.1]))).pow
  have hderiv :
      ∀ x ∈ Ioo (0 : ℝ) 1,
        HasDerivAt plusDenomPrimitive24 (plusDenomKernel24 x) x := by
    intro x hx
    have hne : 1 + x ≠ 0 := by linarith [hx.1]
    have hinner :
        HasDerivAt (fun y : ℝ => 1 + y) 1 x := by
      convert (hasDerivAt_const x 1).add (hasDerivAt_id x) using 1
      simp
    have hlog :
        HasDerivAt (fun y : ℝ => Real.log (1 + y))
          (1 / (1 + x)) x := by
      convert hinner.log hne using 1
    unfold plusDenomPrimitive24 plusDenomKernel24
    convert (hlog.pow 2).const_mul (1 / 2 : ℝ) using 1
    field_simp [hne]
    ring
  have hint :
      IntervalIntegrable plusDenomKernel24
        MeasureTheory.volume 0 1 := by
    apply ContinuousOn.intervalIntegrable
    rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
    unfold plusDenomKernel24
    apply ContinuousOn.div
    · apply (continuousOn_const.add continuousOn_id).log
      intro x hx
      exact ne_of_gt (show 0 < 1 + x by linarith [hx.1])
    · exact continuousOn_const.add continuousOn_id
    · intro x hx
      exact ne_of_gt (show 0 < 1 + x by linarith [hx.1])
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le
    (f := plusDenomPrimitive24)
    (f' := plusDenomKernel24)
    (by norm_num) hcont hderiv hint]
  norm_num [plusDenomPrimitive24]

private def plusSquareDenomKernel24 (x : ℝ) : ℝ :=
  Real.log (1 + x) ^ 2 / (1 + x)

private def plusSquareDenomPrimitive24 (x : ℝ) : ℝ :=
  (1 / 3 : ℝ) * Real.log (1 + x) ^ 3

private theorem plusSquareDenomIntegral24 :
    (∫ x : ℝ in 0..1, plusSquareDenomKernel24 x) =
      (1 / 3 : ℝ) * Real.log 2 ^ 3 := by
  have hcont :
      ContinuousOn plusSquareDenomPrimitive24 (Icc (0 : ℝ) 1) := by
    unfold plusSquareDenomPrimitive24
    apply continuousOn_const.mul
    apply ((continuousOn_const.add continuousOn_id).log
      (fun x hx => ne_of_gt
        (show 0 < 1 + x by linarith [hx.1]))).pow
  have hderiv :
      ∀ x ∈ Ioo (0 : ℝ) 1,
        HasDerivAt plusSquareDenomPrimitive24
          (plusSquareDenomKernel24 x) x := by
    intro x hx
    have hne : 1 + x ≠ 0 := by linarith [hx.1]
    have hinner :
        HasDerivAt (fun y : ℝ => 1 + y) 1 x := by
      convert (hasDerivAt_const x 1).add (hasDerivAt_id x) using 1
      simp
    have hlog :
        HasDerivAt (fun y : ℝ => Real.log (1 + y))
          (1 / (1 + x)) x := by
      convert hinner.log hne using 1
    unfold plusSquareDenomPrimitive24 plusSquareDenomKernel24
    convert (hlog.pow 3).const_mul (1 / 3 : ℝ) using 1
    field_simp [hne]
    ring
  have hint :
      IntervalIntegrable plusSquareDenomKernel24
        MeasureTheory.volume 0 1 := by
    apply ContinuousOn.intervalIntegrable
    rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
    unfold plusSquareDenomKernel24
    apply ContinuousOn.div
    · exact
        (((continuousOn_const.add continuousOn_id).log
          (fun x hx => ne_of_gt
            (show 0 < 1 + x by linarith [hx.1]))).pow 2)
    · exact continuousOn_const.add continuousOn_id
    · intro x hx
      exact ne_of_gt (show 0 < 1 + x by linarith [hx.1])
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le
    (f := plusSquareDenomPrimitive24)
    (f' := plusSquareDenomKernel24)
    (by norm_num) hcont hderiv hint]
  norm_num [plusSquareDenomPrimitive24]

private def mobiusMinusComponent24 (x : ℝ) : ℝ :=
  Real.log 2 * (minusSimpleKernel24 x - minusDenomKernel24 x)

private def mobiusCrossComponent24 (x : ℝ) : ℝ :=
  crossRadialKernel24 x - crossPlusKernel24 x

private def mobiusPlusComponent24 (x : ℝ) : ℝ :=
  Real.log 2 * (plusSimpleKernel24 x - plusDenomKernel24 x)

private def mobiusSquareComponent24 (x : ℝ) : ℝ :=
  plusLogSquareKernel24 x - plusSquareDenomKernel24 x

private def mobiusTransformedKernel24 (x : ℝ) : ℝ :=
  mobiusMinusComponent24 x - mobiusCrossComponent24 x -
    mobiusPlusComponent24 x + mobiusSquareComponent24 x

private theorem mobiusMinusComponent24_intervalIntegrable :
    IntervalIntegrable mobiusMinusComponent24
      MeasureTheory.volume 0 1 := by
  unfold mobiusMinusComponent24
  exact
    (minusSimpleKernel24_intervalIntegrable.sub
      minusDenomKernel24_intervalIntegrable).const_mul _

private theorem mobiusCrossComponent24_intervalIntegrable :
    IntervalIntegrable mobiusCrossComponent24
      MeasureTheory.volume 0 1 := by
  unfold mobiusCrossComponent24
  exact
    crossRadialKernel24_intervalIntegrable.sub
      crossPlusKernel24_intervalIntegrable

private theorem mobiusPlusComponent24_intervalIntegrable :
    IntervalIntegrable mobiusPlusComponent24
      MeasureTheory.volume 0 1 := by
  unfold mobiusPlusComponent24
  exact
    (plusSimpleKernel24_intervalIntegrable.sub
      (by
        unfold plusDenomKernel24
        apply ContinuousOn.intervalIntegrable
        rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
        apply ContinuousOn.div
        · apply (continuousOn_const.add continuousOn_id).log
          intro x hx
          exact ne_of_gt (show 0 < 1 + x by linarith [hx.1])
        · exact continuousOn_const.add continuousOn_id
        · intro x hx
          exact ne_of_gt (show 0 < 1 + x by linarith [hx.1]))).const_mul _

private theorem mobiusSquareComponent24_intervalIntegrable :
    IntervalIntegrable mobiusSquareComponent24
      MeasureTheory.volume 0 1 := by
  unfold mobiusSquareComponent24
  exact
    plusLogSquareKernel24_intervalIntegrable.sub
      (by
        unfold plusSquareDenomKernel24
        apply ContinuousOn.intervalIntegrable
        rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
        apply ContinuousOn.div
        · exact
            (((continuousOn_const.add continuousOn_id).log
              (fun x hx => ne_of_gt
                (show 0 < 1 + x by linarith [hx.1]))).pow 2)
        · exact continuousOn_const.add continuousOn_id
        · intro x hx
          exact ne_of_gt (show 0 < 1 + x by linarith [hx.1]))

private theorem mobiusTransformedKernel24_intervalIntegrable :
    IntervalIntegrable mobiusTransformedKernel24
      MeasureTheory.volume 0 1 := by
  unfold mobiusTransformedKernel24
  exact
    ((mobiusMinusComponent24_intervalIntegrable.sub
      mobiusCrossComponent24_intervalIntegrable).sub
      mobiusPlusComponent24_intervalIntegrable).add
      mobiusSquareComponent24_intervalIntegrable

private theorem mobiusTransformedIntegral24 :
    (∫ x : ℝ in 0..1, mobiusTransformedKernel24 x) =
      zeta3_24 -
        (3 / 2 : ℝ) * Real.log 2 * (Real.pi ^ 2 / 6) := by
  have hPlusDen :
      IntervalIntegrable plusDenomKernel24
        MeasureTheory.volume 0 1 := by
    unfold plusDenomKernel24
    apply ContinuousOn.intervalIntegrable
    rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
    apply ContinuousOn.div
    · apply (continuousOn_const.add continuousOn_id).log
      intro x hx
      exact ne_of_gt (show 0 < 1 + x by linarith [hx.1])
    · exact continuousOn_const.add continuousOn_id
    · intro x hx
      exact ne_of_gt (show 0 < 1 + x by linarith [hx.1])
  have hPlusSqDen :
      IntervalIntegrable plusSquareDenomKernel24
        MeasureTheory.volume 0 1 := by
    unfold plusSquareDenomKernel24
    apply ContinuousOn.intervalIntegrable
    rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
    apply ContinuousOn.div
    · exact
        (((continuousOn_const.add continuousOn_id).log
          (fun x hx => ne_of_gt
            (show 0 < 1 + x by linarith [hx.1]))).pow 2)
    · exact continuousOn_const.add continuousOn_id
    · intro x hx
      exact ne_of_gt (show 0 < 1 + x by linarith [hx.1])
  have hPlusSqRadial :
      (∫ x : ℝ in 0..1, plusLogSquareKernel24 x) =
        (1 / 4 : ℝ) * zeta3_24 := by
    simpa [plusLogSquareKernel24,
      RamanujanChallenge.P26.zeta3, zeta3_24] using
      RamanujanChallenge.P26.logSquareOnePlusIntegral26
  rw [show
      (fun x : ℝ => mobiusTransformedKernel24 x) =
        (fun x : ℝ =>
          (mobiusMinusComponent24 x - mobiusCrossComponent24 x -
            mobiusPlusComponent24 x) + mobiusSquareComponent24 x) by
      funext x
      rfl]
  rw [intervalIntegral.integral_add
    ((mobiusMinusComponent24_intervalIntegrable.sub
      mobiusCrossComponent24_intervalIntegrable).sub
      mobiusPlusComponent24_intervalIntegrable)
    mobiusSquareComponent24_intervalIntegrable]
  rw [intervalIntegral.integral_sub
    (mobiusMinusComponent24_intervalIntegrable.sub
      mobiusCrossComponent24_intervalIntegrable)
    mobiusPlusComponent24_intervalIntegrable]
  rw [intervalIntegral.integral_sub
    mobiusMinusComponent24_intervalIntegrable
    mobiusCrossComponent24_intervalIntegrable]
  unfold mobiusMinusComponent24 mobiusCrossComponent24
    mobiusPlusComponent24 mobiusSquareComponent24
  rw [intervalIntegral.integral_const_mul,
    intervalIntegral.integral_sub
      minusSimpleKernel24_intervalIntegrable
      minusDenomKernel24_intervalIntegrable,
    intervalIntegral.integral_sub
      crossRadialKernel24_intervalIntegrable
      crossPlusKernel24_intervalIntegrable,
    intervalIntegral.integral_const_mul,
    intervalIntegral.integral_sub
      plusSimpleKernel24_intervalIntegrable
      hPlusDen,
    intervalIntegral.integral_sub
      plusLogSquareKernel24_intervalIntegrable
      hPlusSqDen]
  rw [minusSimpleIntegral24, minusDenomIntegral24,
    crossRadialIntegral24, crossPlusIntegral24,
    plusSimpleIntegral24, plusDenomIntegral24,
    hPlusSqRadial,
    plusSquareDenomIntegral24]
  ring

private def coreCrossKernel24 (x : ℝ) : ℝ :=
  Real.log x * Real.log (1 + x) / (1 - x)

private def mobiusMap24 (t : ℝ) : ℝ :=
  (1 - t) / (1 + t)

private def mobiusMapDeriv24 (t : ℝ) : ℝ :=
  -2 / (1 + t) ^ 2

private theorem mobiusMap24_hasDerivAt
    {t : ℝ} (ht0 : 0 < t) (ht1 : t < 1) :
    HasDerivAt mobiusMap24 (mobiusMapDeriv24 t) t := by
  have hden : 1 + t ≠ 0 := by linarith
  have hnum :
      HasDerivAt (fun x : ℝ => 1 - x) (-1) t := by
    convert (hasDerivAt_const t 1).sub (hasDerivAt_id t) using 1
    simp
  have hdenDeriv :
      HasDerivAt (fun x : ℝ => 1 + x) 1 t := by
    convert (hasDerivAt_const t 1).add (hasDerivAt_id t) using 1
    simp
  unfold mobiusMap24 mobiusMapDeriv24
  convert hnum.div hdenDeriv hden using 1
  field_simp [hden]
  ring

private theorem mobiusChangeIntegrand24
    {t : ℝ} (ht0 : 0 ≤ t) (ht1 : t ≤ 1) :
    (coreCrossKernel24 ∘ mobiusMap24) t * mobiusMapDeriv24 t =
      -mobiusTransformedKernel24 t := by
  rcases ht0.eq_or_lt with rfl | ht0
  · norm_num [coreCrossKernel24, mobiusMap24, mobiusMapDeriv24,
      mobiusTransformedKernel24, mobiusMinusComponent24,
      mobiusCrossComponent24, mobiusPlusComponent24,
      mobiusSquareComponent24, minusSimpleKernel24,
      minusDenomKernel24, crossRadialKernel24, crossPlusKernel24,
      plusSimpleKernel24, plusDenomKernel24, plusLogSquareKernel24,
      plusSquareDenomKernel24]
  rcases ht1.eq_or_lt with rfl | ht1
  · norm_num [coreCrossKernel24, mobiusMap24, mobiusMapDeriv24,
      mobiusTransformedKernel24, mobiusMinusComponent24,
      mobiusCrossComponent24, mobiusPlusComponent24,
      mobiusSquareComponent24, minusSimpleKernel24,
      minusDenomKernel24, crossRadialKernel24, crossPlusKernel24,
      plusSimpleKernel24, plusDenomKernel24, plusLogSquareKernel24,
      plusSquareDenomKernel24]
    ring
  have htne : t ≠ 0 := ne_of_gt ht0
  have h1mt0 : 0 < 1 - t := by linarith
  have h1pt0 : 0 < 1 + t := by linarith
  have h1ptne : 1 + t ≠ 0 := ne_of_gt h1pt0
  have hphi0 : 0 < mobiusMap24 t := by
    unfold mobiusMap24
    positivity
  have hlogPhi :
      Real.log (mobiusMap24 t) =
        Real.log (1 - t) - Real.log (1 + t) := by
    unfold mobiusMap24
    rw [Real.log_div (ne_of_gt h1mt0) h1ptne]
  have hOnePlusPhi :
      1 + mobiusMap24 t = 2 / (1 + t) := by
    unfold mobiusMap24
    field_simp [h1ptne]
    ring
  have hlogOnePlusPhi :
      Real.log (1 + mobiusMap24 t) =
        Real.log 2 - Real.log (1 + t) := by
    rw [hOnePlusPhi,
      Real.log_div (by norm_num : (2 : ℝ) ≠ 0) h1ptne]
  have hOneMinusPhi :
      1 - mobiusMap24 t = 2 * t / (1 + t) := by
    unfold mobiusMap24
    field_simp [h1ptne]
    ring
  simp only [Function.comp_apply]
  unfold coreCrossKernel24 mobiusMapDeriv24
    mobiusTransformedKernel24 mobiusMinusComponent24
    mobiusCrossComponent24 mobiusPlusComponent24
    mobiusSquareComponent24 minusSimpleKernel24
    minusDenomKernel24 crossRadialKernel24 crossPlusKernel24
    plusSimpleKernel24 plusDenomKernel24 plusLogSquareKernel24
    plusSquareDenomKernel24
  rw [hlogPhi, hlogOnePlusPhi, hOneMinusPhi]
  field_simp [htne, h1ptne]
  ring

private theorem coreCrossIntegral24 :
    (∫ x : ℝ in 0..1, coreCrossKernel24 x) =
      zeta3_24 -
        (3 / 2 : ℝ) * Real.log 2 * (Real.pi ^ 2 / 6) := by
  have hmapCont :
      ContinuousOn mobiusMap24 (Set.uIcc (0 : ℝ) 1) := by
    rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
    unfold mobiusMap24
    apply ContinuousOn.div
    · fun_prop
    · fun_prop
    · intro t ht
      exact ne_of_gt (show 0 < 1 + t by linarith [ht.1])
  have hderiv :
      ∀ t ∈ Ioo (min (0 : ℝ) 1) (max (0 : ℝ) 1),
        HasDerivAt mobiusMap24 (mobiusMapDeriv24 t) t := by
    intro t ht
    norm_num at ht
    exact mobiusMap24_hasDerivAt ht.1 ht.2
  have hnonpos :
      ∀ t ∈ Ioo (min (0 : ℝ) 1) (max (0 : ℝ) 1),
        mobiusMapDeriv24 t ≤ 0 := by
    intro t ht
    unfold mobiusMapDeriv24
    exact div_nonpos_of_nonpos_of_nonneg (by norm_num) (sq_nonneg _)
  have hsubst :=
    intervalIntegral.integral_comp_mul_deriv_of_deriv_nonpos
      (a := (0 : ℝ)) (b := 1)
      (f := mobiusMap24) (f' := mobiusMapDeriv24)
      (g := coreCrossKernel24)
      hmapCont hderiv hnonpos
  have hcongr :
      (∫ t : ℝ in 0..1,
        (coreCrossKernel24 ∘ mobiusMap24) t *
          mobiusMapDeriv24 t) =
        ∫ t : ℝ in 0..1, -mobiusTransformedKernel24 t := by
    apply intervalIntegral.integral_congr
    intro t ht
    have ht' : t ∈ Icc (0 : ℝ) 1 := by
      simpa [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using ht
    exact mobiusChangeIntegrand24 ht'.1 ht'.2
  have hchange :
      (∫ x : ℝ in 0..1, coreCrossKernel24 x) =
        ∫ t : ℝ in 0..1, mobiusTransformedKernel24 t := by
    have hs :
        (∫ t : ℝ in 0..1,
          (coreCrossKernel24 ∘ mobiusMap24) t *
            mobiusMapDeriv24 t) =
          ∫ x : ℝ in 1..0, coreCrossKernel24 x := by
      simpa [mobiusMap24] using hsubst
    rw [hcongr, intervalIntegral.integral_neg] at hs
    have hsymm :
        (∫ x : ℝ in 1..0, coreCrossKernel24 x) =
      -(∫ x : ℝ in 0..1, coreCrossKernel24 x) := by
      rw [intervalIntegral.integral_symm]
    rw [hsymm] at hs
    linarith
  rw [hchange, mobiusTransformedIntegral24]

private theorem coreCrossKernel24_intervalIntegrable :
    IntervalIntegrable coreCrossKernel24
      MeasureTheory.volume 0 1 := by
  apply IntervalIntegrable.trans (b := (1 / 2 : ℝ))
  · have hlog :
        IntervalIntegrable Real.log MeasureTheory.volume
          (0 : ℝ) (1 / 2) :=
      intervalIntegral.intervalIntegrable_log'
    have hfactor :
        ContinuousOn
          (fun x : ℝ => Real.log (1 + x) / (1 - x))
          (Set.uIcc (0 : ℝ) (1 / 2)) := by
      have hfactorIcc :
          ContinuousOn
            (fun x : ℝ => Real.log (1 + x) / (1 - x))
            (Icc (0 : ℝ) (1 / 2)) := by
        apply ContinuousOn.div
        · apply (continuousOn_const.add continuousOn_id).log
          intro x hx
          exact ne_of_gt (show 0 < 1 + x by linarith [hx.1])
        · exact continuousOn_const.sub continuousOn_id
        · intro x hx
          exact ne_of_gt (show 0 < 1 - x by linarith [hx.2])
      rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2)]
      exact hfactorIcc
    have hint := hlog.continuousOn_mul hfactor
    apply IntervalIntegrable.congr
      (f := fun x : ℝ =>
        (Real.log (1 + x) / (1 - x)) * Real.log x) ?_ hint
    intro x _
    unfold coreCrossKernel24
    ring
  · have hcont :
        ContinuousOn
          (fun x : ℝ =>
            logOneMinusSlope24 (1 - x) * Real.log (1 + x))
          (Icc (1 / 2 : ℝ) 1) := by
      have hsub :
          ContinuousOn (fun x : ℝ => 1 - x)
            (Icc (1 / 2 : ℝ) 1) := by
        fun_prop
      have hsubmem :
          MapsTo (fun x : ℝ => 1 - x)
            (Icc (1 / 2 : ℝ) 1) (Icc (0 : ℝ) (1 / 2)) := by
        intro x hx
        constructor <;> linarith [hx.1, hx.2]
      have hslope :=
        logOneMinusSlope24_continuousOn.comp hsub hsubmem
      have hlog :
          ContinuousOn (fun x : ℝ => Real.log (1 + x))
            (Icc (1 / 2 : ℝ) 1) := by
        apply (continuousOn_const.add continuousOn_id).log
        intro x hx
        exact ne_of_gt (show 0 < 1 + x by linarith [hx.1])
      exact hslope.mul hlog
    have hint :
        IntervalIntegrable
          (fun x : ℝ =>
            logOneMinusSlope24 (1 - x) * Real.log (1 + x))
          MeasureTheory.volume (1 / 2) 1 := by
      apply ContinuousOn.intervalIntegrable
      rw [Set.uIcc_of_le (by norm_num : (1 / 2 : ℝ) ≤ 1)]
      exact hcont
    apply hint.congr_ae
    filter_upwards [
      MeasureTheory.ae_restrict_mem measurableSet_uIoc,
      MeasureTheory.ae_restrict_of_ae
        (MeasureTheory.Measure.ae_ne MeasureTheory.volume (1 : ℝ))
    ] with x hx hxne
    have hx' : x ∈ Ioc (1 / 2 : ℝ) 1 := by
      rw [Set.uIoc_of_le (by norm_num : (1 / 2 : ℝ) ≤ 1)] at hx
      exact hx
    have hxlt : x < 1 := lt_of_le_of_ne hx'.2 hxne
    have h1xne : 1 - x ≠ 0 := ne_of_gt (sub_pos.mpr hxlt)
    simp [coreCrossKernel24, logOneMinusSlope24, h1xne]
    field_simp [h1xne]

private def crossAltKernel24 (x : ℝ) : ℝ :=
  Real.log x * Real.log (1 - x) / (1 + x)

private theorem crossAltKernel24_intervalIntegrable :
    IntervalIntegrable crossAltKernel24
      MeasureTheory.volume 0 1 := by
  apply IntervalIntegrable.trans (b := (1 / 2 : ℝ))
  · have hcont :
        ContinuousOn
          (fun x : ℝ =>
            (x * Real.log x) * logOneMinusSlope24 x / (1 + x))
          (Icc (0 : ℝ) (1 / 2)) := by
      apply ContinuousOn.div
      · exact
          Real.continuous_mul_log.continuousOn.mul
            logOneMinusSlope24_continuousOn
      · exact continuousOn_const.add continuousOn_id
      · intro x hx
        exact ne_of_gt (show 0 < 1 + x by linarith [hx.1])
    have hint :
        IntervalIntegrable
          (fun x : ℝ =>
            (x * Real.log x) * logOneMinusSlope24 x / (1 + x))
          MeasureTheory.volume 0 (1 / 2) := by
      apply ContinuousOn.intervalIntegrable
      rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2)]
      exact hcont
    apply IntervalIntegrable.congr
      (f := fun x : ℝ =>
        (x * Real.log x) * logOneMinusSlope24 x / (1 + x)) ?_ hint
    intro x hx
    have hx' : x ∈ Ioc (0 : ℝ) (1 / 2) := by
      simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2)] using hx
    have hxne : x ≠ 0 := ne_of_gt hx'.1
    simp [crossAltKernel24, logOneMinusSlope24, hxne]
    field_simp [hxne]
  · have hsub :
        ContinuousOn (fun x : ℝ => 1 - x)
          (Icc (1 / 2 : ℝ) 1) := by
      fun_prop
    have hsubmem :
        MapsTo (fun x : ℝ => 1 - x)
          (Icc (1 / 2 : ℝ) 1) (Icc (0 : ℝ) (1 / 2)) := by
      intro x hx
      constructor <;> linarith [hx.1, hx.2]
    have hslope :=
      logOneMinusSlope24_continuousOn.comp hsub hsubmem
    have hmulLog :
        ContinuousOn
          (fun x : ℝ => (1 - x) * Real.log (1 - x))
          (Icc (1 / 2 : ℝ) 1) := by
      simpa [Function.comp_def] using
        Real.continuous_mul_log.continuousOn.comp hsub hsubmem
    have hcont :
        ContinuousOn
          (fun x : ℝ =>
            logOneMinusSlope24 (1 - x) *
              ((1 - x) * Real.log (1 - x)) / (1 + x))
          (Icc (1 / 2 : ℝ) 1) := by
      apply ContinuousOn.div (hslope.mul hmulLog)
        (continuousOn_const.add continuousOn_id)
      intro x hx
      exact ne_of_gt (show 0 < 1 + x by linarith [hx.1])
    have hint :
        IntervalIntegrable
          (fun x : ℝ =>
            logOneMinusSlope24 (1 - x) *
              ((1 - x) * Real.log (1 - x)) / (1 + x))
          MeasureTheory.volume (1 / 2) 1 := by
      apply ContinuousOn.intervalIntegrable
      rw [Set.uIcc_of_le (by norm_num : (1 / 2 : ℝ) ≤ 1)]
      exact hcont
    apply IntervalIntegrable.congr
      (f := fun x : ℝ =>
        logOneMinusSlope24 (1 - x) *
          ((1 - x) * Real.log (1 - x)) / (1 + x)) ?_ hint
    intro x hx
    have hx' : x ∈ Ioc (1 / 2 : ℝ) 1 := by
      rw [Set.uIoc_of_le (by norm_num : (1 / 2 : ℝ) ≤ 1)] at hx
      exact hx
    rcases hx'.2.eq_or_lt with rfl | hxlt
    · simp [crossAltKernel24]
    · have h1xne : 1 - x ≠ 0 := ne_of_gt (sub_pos.mpr hxlt)
      simp [crossAltKernel24, logOneMinusSlope24, h1xne]
      field_simp [h1xne]

private def crossDifferencePrimitive24 (x : ℝ) : ℝ :=
  Real.log x * Real.log (1 - x) * Real.log (1 + x)

private def crossDifferenceKernel24 (x : ℝ) : ℝ :=
  crossRadialKernel24 x - coreCrossKernel24 x + crossAltKernel24 x

private theorem crossDifferencePrimitive24_hasDerivAt
    {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasDerivAt crossDifferencePrimitive24
      (crossDifferenceKernel24 x) x := by
  have hxne : x ≠ 0 := ne_of_gt hx0
  have h1mxne : 1 - x ≠ 0 := by linarith
  have h1pxne : 1 + x ≠ 0 := by linarith
  have hlogx :
      HasDerivAt (fun y : ℝ => Real.log y) (1 / x) x := by
    simpa [one_div] using Real.hasDerivAt_log hxne
  have hminus :
      HasDerivAt (fun y : ℝ => Real.log (1 - y))
        (-1 / (1 - x)) x := by
    have hinner :
        HasDerivAt (fun y : ℝ => 1 - y) (-1) x := by
      convert (hasDerivAt_const x 1).sub (hasDerivAt_id x) using 1
      simp
    convert hinner.log h1mxne using 1
  have hplus :
      HasDerivAt (fun y : ℝ => Real.log (1 + y))
        (1 / (1 + x)) x := by
    have hinner :
        HasDerivAt (fun y : ℝ => 1 + y) 1 x := by
      convert (hasDerivAt_const x 1).add (hasDerivAt_id x) using 1
      simp
    convert hinner.log h1pxne using 1
  unfold crossDifferencePrimitive24
  convert (hlogx.mul hminus).mul hplus using 1
  unfold crossDifferenceKernel24 crossRadialKernel24
    coreCrossKernel24 crossAltKernel24
  simp only [Pi.mul_apply]
  field_simp [hxne, h1mxne, h1pxne]
  ring

private theorem crossDifferencePrimitive24_tendsto_zero :
    Tendsto crossDifferencePrimitive24
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  have hcont :
      ContinuousOn crossDifferencePrimitive24
        (Icc (0 : ℝ) (1 / 2)) := by
    have haux :
        ContinuousOn
          (fun x : ℝ =>
            ((x * Real.log x) * x) * logOneMinusSlope24 x *
              RamanujanChallenge.P26.logOnePlusSlope26 x)
          (Icc (0 : ℝ) (1 / 2)) :=
      (((Real.continuous_mul_log.continuousOn.mul continuousOn_id).mul
        logOneMinusSlope24_continuousOn).mul
        (RamanujanChallenge.P26.logOnePlusSlope26_continuousOn.mono
          (by intro x hx; exact ⟨hx.1, hx.2.trans (by norm_num)⟩)))
    apply haux.congr
    intro x _
    by_cases hxzero : x = 0
    · subst x
      simp [crossDifferencePrimitive24]
    · simp [crossDifferencePrimitive24, logOneMinusSlope24,
        RamanujanChallenge.P26.logOnePlusSlope26, hxzero]
      field_simp [hxzero]
  have hwithin :
      ContinuousWithinAt crossDifferencePrimitive24
        (Ioi (0 : ℝ)) 0 :=
    (hcont 0 (by norm_num)).mono_of_mem_nhdsWithin
      (Icc_mem_nhdsGT (by norm_num : (0 : ℝ) < 1 / 2))
  simpa [crossDifferencePrimitive24] using hwithin.tendsto

private theorem crossDifferencePrimitive24_tendsto_one :
    Tendsto crossDifferencePrimitive24
      (𝓝[<] (1 : ℝ)) (𝓝 0) := by
  have hsub :
      ContinuousOn (fun x : ℝ => 1 - x)
        (Icc (1 / 2 : ℝ) 1) := by
    fun_prop
  have hsubmem :
      MapsTo (fun x : ℝ => 1 - x)
        (Icc (1 / 2 : ℝ) 1) (Icc (0 : ℝ) (1 / 2)) := by
    intro x hx
    constructor <;> linarith [hx.1, hx.2]
  have haux :
      ContinuousOn
        (fun x : ℝ =>
          logOneMinusSlope24 (1 - x) *
            ((1 - x) * Real.log (1 - x)) *
            Real.log (1 + x))
        (Icc (1 / 2 : ℝ) 1) := by
    have hslope :=
      logOneMinusSlope24_continuousOn.comp hsub hsubmem
    have hmulLog :
        ContinuousOn
          (fun x : ℝ => (1 - x) * Real.log (1 - x))
          (Icc (1 / 2 : ℝ) 1) := by
      simpa [Function.comp_def] using
        Real.continuous_mul_log.continuousOn.comp hsub hsubmem
    have hplus :
        ContinuousOn (fun x : ℝ => Real.log (1 + x))
          (Icc (1 / 2 : ℝ) 1) := by
      apply (continuousOn_const.add continuousOn_id).log
      intro x hx
      exact ne_of_gt (show 0 < 1 + x by linarith [hx.1])
    exact (hslope.mul hmulLog).mul hplus
  have hcont :
      ContinuousOn crossDifferencePrimitive24
        (Icc (1 / 2 : ℝ) 1) := by
    apply haux.congr
    intro x hx
    by_cases hxone : x = 1
    · subst x
      simp [crossDifferencePrimitive24]
    · have h1xne : 1 - x ≠ 0 := sub_ne_zero.mpr (Ne.symm hxone)
      simp [crossDifferencePrimitive24, logOneMinusSlope24, h1xne]
      field_simp [h1xne]
      exact Or.inl trivial
  have hwithin :
      ContinuousWithinAt crossDifferencePrimitive24
        (Iio (1 : ℝ)) 1 :=
    (hcont 1 (by norm_num)).mono_of_mem_nhdsWithin
      (Icc_mem_nhdsLT (by norm_num : (1 / 2 : ℝ) < 1))
  simpa [crossDifferencePrimitive24] using hwithin.tendsto

private theorem crossDifferenceIntegral24 :
    (∫ x : ℝ in 0..1, crossDifferenceKernel24 x) = 0 := by
  have hint :
      IntervalIntegrable crossDifferenceKernel24
        MeasureTheory.volume 0 1 := by
    unfold crossDifferenceKernel24
    exact
      (crossRadialKernel24_intervalIntegrable.sub
        coreCrossKernel24_intervalIntegrable).add
        crossAltKernel24_intervalIntegrable
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_tendsto
    (f := crossDifferencePrimitive24)
    (fa := (0 : ℝ)) (fb := 0)
    (by norm_num)
    (fun x hx => crossDifferencePrimitive24_hasDerivAt hx.1 hx.2)
    hint crossDifferencePrimitive24_tendsto_zero
    crossDifferencePrimitive24_tendsto_one]
  ring

private theorem crossAltIntegral24 :
    (∫ x : ℝ in 0..1, crossAltKernel24 x) =
      (13 / 8 : ℝ) * zeta3_24 -
        (3 / 2 : ℝ) * Real.log 2 * (Real.pi ^ 2 / 6) := by
  have h := crossDifferenceIntegral24
  unfold crossDifferenceKernel24 at h
  rw [intervalIntegral.integral_add
      (crossRadialKernel24_intervalIntegrable.sub
        coreCrossKernel24_intervalIntegrable)
      crossAltKernel24_intervalIntegrable,
    intervalIntegral.integral_sub
      crossRadialKernel24_intervalIntegrable
      coreCrossKernel24_intervalIntegrable,
    crossRadialIntegral24, coreCrossIntegral24] at h
  linarith

private theorem reflectedMinusRadialKernel24_intervalIntegrable :
    IntervalIntegrable
      (fun x : ℝ =>
        Real.log x * Real.log (1 - x) / (1 - x))
      MeasureTheory.volume 0 1 := by
  have h :=
    (minusRadialKernel24_intervalIntegrable.comp_sub_left 1).symm
  convert h using 1
  · funext x
    unfold minusRadialKernel24
    ring
  · norm_num
  · norm_num

private theorem quadraticLinearKernel24_integral :
    (∫ x : ℝ in 0..1, quadraticLinearKernel24 x) =
      quadraticLinearEulerValue24 := by
  have hdecomp :
      (∫ x : ℝ in 0..1, quadraticLinearKernel24 x) =
        ∫ x : ℝ in 0..1,
          (minusRadialKernel24 x +
            Real.log x * Real.log (1 - x) / (1 - x)) +
          2 * (RamanujanChallenge.P26.radialWeightThreeKernel26 x +
            coreCrossKernel24 x) := by
    apply intervalIntegral.integral_congr
    intro x hx
    have hx' : x ∈ Icc (0 : ℝ) 1 := by
      simpa [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hx
    by_cases hxzero : x = 0
    · subst x
      simp [quadraticLinearKernel24, minusRadialKernel24,
        RamanujanChallenge.P26.radialWeightThreeKernel26,
        coreCrossKernel24]
    by_cases hxone : x = 1
    · subst x
      simp [quadraticLinearKernel24, minusRadialKernel24,
        RamanujanChallenge.P26.radialWeightThreeKernel26,
        coreCrossKernel24]
    · have hxne : x ≠ 0 := hxzero
      have hxlt : x < 1 := lt_of_le_of_ne hx'.2 hxone
      have h1xne : 1 - x ≠ 0 := ne_of_gt (sub_pos.mpr hxlt)
      unfold quadraticLinearKernel24 minusRadialKernel24
        RamanujanChallenge.P26.radialWeightThreeKernel26
        coreCrossKernel24
      field_simp [hxne, h1xne]
      ring
  rw [hdecomp,
    intervalIntegral.integral_add
      (minusRadialKernel24_intervalIntegrable.add
        reflectedMinusRadialKernel24_intervalIntegrable)
      ((RamanujanChallenge.P26.radialWeightThreeKernel_intervalIntegrable26.add
        coreCrossKernel24_intervalIntegrable).const_mul 2),
    intervalIntegral.integral_add
      minusRadialKernel24_intervalIntegrable
      reflectedMinusRadialKernel24_intervalIntegrable,
    intervalIntegral.integral_const_mul,
    intervalIntegral.integral_add
      RamanujanChallenge.P26.radialWeightThreeKernel_intervalIntegrable26
      coreCrossKernel24_intervalIntegrable,
    minusRadialIntegral24, reflectedMinusRadialIntegral24,
    RamanujanChallenge.P26.radialWeightThreeIntegral26,
    coreCrossIntegral24]
  unfold quadraticLinearEulerValue24
    RamanujanChallenge.P26.zeta3 zeta3_24
  ring

private theorem alternatingQuadraticLinearKernel24_integral :
    (∫ x : ℝ in 0..1, alternatingQuadraticLinearKernel24 x) =
      alternatingQuadraticLinearEulerValue24 := by
  have hdecomp :
      (∫ x : ℝ in 0..1, alternatingQuadraticLinearKernel24 x) =
        ∫ x : ℝ in 0..1,
          (RamanujanChallenge.P26.radialWeightThreeKernel26 x -
            RamanujanChallenge.P26.alternatingWeightThreeKernel26 x) +
          2 * (minusRadialKernel24 x - crossAltKernel24 x) := by
    apply intervalIntegral.integral_congr
    intro x hx
    have hx' : x ∈ Icc (0 : ℝ) 1 := by
      simpa [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hx
    by_cases hxzero : x = 0
    · subst x
      simp [alternatingQuadraticLinearKernel24,
        RamanujanChallenge.P26.radialWeightThreeKernel26,
        RamanujanChallenge.P26.alternatingWeightThreeKernel26,
        minusRadialKernel24, crossAltKernel24]
    by_cases hxone : x = 1
    · subst x
      simp [alternatingQuadraticLinearKernel24,
        RamanujanChallenge.P26.radialWeightThreeKernel26,
        RamanujanChallenge.P26.alternatingWeightThreeKernel26,
        minusRadialKernel24, crossAltKernel24]
    · have hxne : x ≠ 0 := hxzero
      have hx0 : 0 < x := lt_of_le_of_ne hx'.1 (Ne.symm hxzero)
      have h1pxne : 1 + x ≠ 0 := by linarith
      unfold alternatingQuadraticLinearKernel24
        RamanujanChallenge.P26.radialWeightThreeKernel26
        RamanujanChallenge.P26.alternatingWeightThreeKernel26
        minusRadialKernel24 crossAltKernel24
      field_simp [hxne, h1pxne]
      ring
  rw [hdecomp,
    intervalIntegral.integral_add
      (RamanujanChallenge.P26.radialWeightThreeKernel_intervalIntegrable26.sub
        RamanujanChallenge.P26.alternatingWeightThreeKernel_intervalIntegrable26)
      ((minusRadialKernel24_intervalIntegrable.sub
        crossAltKernel24_intervalIntegrable).const_mul 2),
    intervalIntegral.integral_sub
      RamanujanChallenge.P26.radialWeightThreeKernel_intervalIntegrable26
      RamanujanChallenge.P26.alternatingWeightThreeKernel_intervalIntegrable26,
    intervalIntegral.integral_const_mul,
    intervalIntegral.integral_sub
      minusRadialKernel24_intervalIntegrable
      crossAltKernel24_intervalIntegrable,
    RamanujanChallenge.P26.radialWeightThreeIntegral26,
    RamanujanChallenge.P26.alternatingWeightThreeIntegral26,
    minusRadialIntegral24, crossAltIntegral24]
  unfold alternatingQuadraticLinearEulerValue24
    RamanujanChallenge.P26.zeta3 zeta3_24
  ring

theorem quadraticLinearEulerTerm24_hasSum :
    HasSum quadraticLinearEulerTerm24
      quadraticLinearEulerValue24 := by
  rw [← quadraticLinearKernel24_integral]
  exact quadraticLinearEulerTerm24_hasSum_integral

theorem alternatingQuadraticLinearEulerTerm24_hasSum :
    HasSum alternatingQuadraticLinearEulerTerm24
      alternatingQuadraticLinearEulerValue24 := by
  rw [← alternatingQuadraticLinearKernel24_integral]
  exact alternatingQuadraticLinearEulerTerm24_hasSum_integral

theorem shiftedLinearEulerTerm24_hasSum :
    HasSum shiftedLinearEulerTerm24 shiftedLinearEulerValue24 :=
  shiftedLinearEuler24_of_lower_weight
    (pairedAlternatingLinear24_of_harmonic
      pairedAlternatingHarmonicEulerTerm24_hasSum)
    quadraticLinearEulerTerm24_hasSum
    alternatingQuadraticLinearEulerTerm24_hasSum

private def ordinaryHarmonicCubicTerm24 (n : ℕ) : ℝ :=
  harmonicNumber (n + 1) / (n + 1 : ℝ) ^ 3

private theorem ordinaryHarmonicCubicTerm24_norm_le (n : ℕ) :
    ‖ordinaryHarmonicCubicTerm24 n‖ ≤
      harmonicNumber (n + 1) ^ 2 / (n + 1 : ℝ) ^ 2 := by
  unfold ordinaryHarmonicCubicTerm24
  rw [Real.norm_eq_abs, abs_div, abs_pow,
    abs_of_pos (by positivity : (0 : ℝ) < n + 1),
    abs_of_nonneg (harmonicNumber_nonneg (n + 1))]
  calc
    harmonicNumber (n + 1) / (n + 1 : ℝ) ^ 3 ≤
        harmonicNumber (n + 1) ^ 2 / (n + 1 : ℝ) ^ 3 := by
      gcongr
      exact harmonicNumber_succ_le_sq n
    _ ≤ harmonicNumber (n + 1) ^ 2 / (n + 1 : ℝ) ^ 2 := by
      gcongr
      all_goals norm_num

private theorem summable_ordinaryHarmonicCubicTerm24 :
    Summable ordinaryHarmonicCubicTerm24 :=
  summable_harmonicNumber_succ_sq_div.of_norm_bounded
    ordinaryHarmonicCubicTerm24_norm_le

private def ordinaryHarmonicCubicMoment24
    (n : ℕ) (x : ℝ) : ℝ :=
  (1 / (n + 1 : ℝ) ^ 3) *
    ∑ k ∈ Finset.range (n + 1), x ^ k

private theorem ordinaryHarmonicCubicMoment24_intervalIntegrable
    (n : ℕ) :
    IntervalIntegrable (ordinaryHarmonicCubicMoment24 n)
      MeasureTheory.volume 0 1 := by
  have hcont : Continuous (ordinaryHarmonicCubicMoment24 n) := by
    unfold ordinaryHarmonicCubicMoment24
    fun_prop
  exact hcont.intervalIntegrable 0 1

private theorem ordinaryHarmonicCubicMoment24_integral
    (n : ℕ) :
    (∫ x : ℝ in 0..1, ordinaryHarmonicCubicMoment24 n x) =
      ordinaryHarmonicCubicTerm24 n := by
  unfold ordinaryHarmonicCubicMoment24
  rw [intervalIntegral.integral_const_mul,
    intervalIntegral.integral_finset_sum]
  · simp only [integral_pow]
    unfold ordinaryHarmonicCubicTerm24 harmonicNumber
    ring
  · intro k _
    exact (continuous_pow k).intervalIntegrable 0 1

private theorem ordinaryHarmonicCubicMoment24_nonneg
    (n : ℕ) {x : ℝ} (hx0 : 0 ≤ x) :
    0 ≤ ordinaryHarmonicCubicMoment24 n x := by
  unfold ordinaryHarmonicCubicMoment24
  positivity

private theorem ordinaryHarmonicCubicMoment24_integral_norm
    (n : ℕ) :
    (∫ x : ℝ in 0..1, ‖ordinaryHarmonicCubicMoment24 n x‖) =
      ordinaryHarmonicCubicTerm24 n := by
  calc
    (∫ x : ℝ in 0..1, ‖ordinaryHarmonicCubicMoment24 n x‖) =
        ∫ x : ℝ in 0..1, ordinaryHarmonicCubicMoment24 n x := by
      apply intervalIntegral.integral_congr
      intro x hx
      simp only [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1),
        Set.mem_Icc] at hx
      change
        ‖ordinaryHarmonicCubicMoment24 n x‖ =
          ordinaryHarmonicCubicMoment24 n x
      rw [Real.norm_eq_abs,
        abs_of_nonneg (ordinaryHarmonicCubicMoment24_nonneg n hx.1)]
    _ = ordinaryHarmonicCubicTerm24 n :=
      ordinaryHarmonicCubicMoment24_integral n

private theorem
    ordinaryHarmonicCubicMoment24_integral_norm_summable :
    Summable
      (fun n : ℕ =>
        ∫ x : ℝ in 0..1,
          ‖ordinaryHarmonicCubicMoment24 n x‖) := by
  exact summable_ordinaryHarmonicCubicTerm24.congr fun n =>
    (ordinaryHarmonicCubicMoment24_integral_norm n).symm

private theorem trilog26_hasSum24
    {x : ℝ} (hx : |x| ≤ 1) :
    HasSum
      (fun n : ℕ =>
        x ^ (n + 1) / (n + 1 : ℝ) ^ 3)
      (RamanujanChallenge.P26.trilog26 x) := by
  have hs :
      Summable
        (fun n : ℕ =>
          x ^ (n + 1) / (n + 1 : ℝ) ^ 3) := by
    apply shifted_zeta_three_hasSum24.summable.of_norm_bounded
    intro n
    rw [Real.norm_eq_abs, abs_div, abs_pow,
      abs_pow, abs_of_pos (by positivity : (0 : ℝ) < n + 1)]
    have hxpow : |x| ^ (n + 1) ≤ 1 := by
      exact pow_le_one₀ (abs_nonneg x) hx
    gcongr
  unfold RamanujanChallenge.P26.trilog26
  simpa only [Nat.cast_add, Nat.cast_one] using hs.hasSum

private def ordinaryHarmonicCubicKernel24 (x : ℝ) : ℝ :=
  (zeta3_24 - RamanujanChallenge.P26.trilog26 x) / (1 - x)

private theorem ordinaryHarmonicCubicMoment24_hasSum_pointwise
    {x : ℝ} (hx0 : 0 ≤ x) (hx1 : x < 1) :
    HasSum (fun n : ℕ => ordinaryHarmonicCubicMoment24 n x)
      (ordinaryHarmonicCubicKernel24 x) := by
  have hxne : x ≠ 1 := ne_of_lt hx1
  have hdiff :=
    shifted_zeta_three_hasSum24.sub
      (trilog26_hasSum24 (x := x)
        (by rw [abs_of_nonneg hx0]; exact hx1.le))
  have hscaled := hdiff.mul_right (1 / (1 - x))
  convert hscaled using 1
  · funext n
    unfold ordinaryHarmonicCubicMoment24
    rw [geom_sum_eq hxne]
    have hxsub : x - 1 ≠ 0 := sub_ne_zero.mpr hxne
    have hsubx : 1 - x ≠ 0 := sub_ne_zero.mpr hxne.symm
    field_simp [hxsub, hsubx]
    ring
  · unfold ordinaryHarmonicCubicKernel24
    ring

private theorem ordinaryHarmonicCubicTerm24_hasSum_integral :
    HasSum ordinaryHarmonicCubicTerm24
      (∫ x : ℝ in 0..1, ordinaryHarmonicCubicKernel24 x) := by
  have hInt :
      ∀ n : ℕ,
        MeasureTheory.Integrable
          (ordinaryHarmonicCubicMoment24 n)
          (MeasureTheory.volume.restrict (Set.Ioc 0 1)) := by
    intro n
    exact (ordinaryHarmonicCubicMoment24_intervalIntegrable n).1
  have hNorm :
      Summable
        (fun n : ℕ =>
          ∫ x : ℝ in Set.Ioc 0 1,
            ‖ordinaryHarmonicCubicMoment24 n x‖) := by
    simpa only [
      ← intervalIntegral.integral_of_le
        (by norm_num : (0 : ℝ) ≤ 1)] using
      ordinaryHarmonicCubicMoment24_integral_norm_summable
  have h :=
    MeasureTheory.hasSum_integral_of_summable_integral_norm
      (μ := MeasureTheory.volume.restrict (Set.Ioc (0 : ℝ) 1))
      hInt hNorm
  have h' :
      HasSum ordinaryHarmonicCubicTerm24
        (∫ x : ℝ in Set.Ioc 0 1,
          ∑' n : ℕ, ordinaryHarmonicCubicMoment24 n x) := by
    convert h using 1
    funext n
    rw [← intervalIntegral.integral_of_le
      (by norm_num : (0 : ℝ) ≤ 1)]
    exact (ordinaryHarmonicCubicMoment24_integral n).symm
  convert h' using 1
  rw [intervalIntegral.integral_of_le (by norm_num : (0 : ℝ) ≤ 1)]
  apply MeasureTheory.setIntegral_congr_ae measurableSet_Ioc
  filter_upwards [
    MeasureTheory.Measure.ae_ne MeasureTheory.volume (1 : ℝ)] with x hxne hx
  have hxlt : x < 1 := lt_of_le_of_ne hx.2 hxne
  exact
    (ordinaryHarmonicCubicMoment24_hasSum_pointwise
      hx.1.le hxlt).tsum_eq.symm

private def barHarmonicCubicTerm24 (n : ℕ) : ℝ :=
  -signedHarmonic24 (n + 1) / (n + 1 : ℝ) ^ 3

private theorem barHarmonicCubicTerm24_norm_le (n : ℕ) :
    ‖barHarmonicCubicTerm24 n‖ ≤
      ordinaryHarmonicCubicTerm24 n := by
  unfold barHarmonicCubicTerm24 ordinaryHarmonicCubicTerm24
  rw [Real.norm_eq_abs, abs_div, abs_neg, abs_pow,
    abs_of_pos (by positivity : (0 : ℝ) < n + 1)]
  exact div_le_div_of_nonneg_right
    (abs_signedHarmonic24_le_harmonicNumber (n + 1))
    (by positivity)

private theorem summable_barHarmonicCubicTerm24 :
    Summable barHarmonicCubicTerm24 :=
  summable_ordinaryHarmonicCubicTerm24.of_norm_bounded
    barHarmonicCubicTerm24_norm_le

private def barHarmonicCubicMoment24
    (n : ℕ) (x : ℝ) : ℝ :=
  (1 / (n + 1 : ℝ) ^ 3) *
    ∑ k ∈ Finset.range (n + 1), (-x) ^ k

private theorem barHarmonicCubicMoment24_intervalIntegrable
    (n : ℕ) :
    IntervalIntegrable (barHarmonicCubicMoment24 n)
      MeasureTheory.volume 0 1 := by
  have hcont : Continuous (barHarmonicCubicMoment24 n) := by
    unfold barHarmonicCubicMoment24
    fun_prop
  exact hcont.intervalIntegrable 0 1

private theorem barHarmonicCubicMoment24_integral
    (n : ℕ) :
    (∫ x : ℝ in 0..1, barHarmonicCubicMoment24 n x) =
      barHarmonicCubicTerm24 n := by
  have hterm (k : ℕ) :
      (∫ x : ℝ in 0..1, (-x) ^ k) =
        (-1 : ℝ) ^ k / (k + 1 : ℝ) := by
    rw [show (fun x : ℝ => (-x) ^ k) =
      fun x : ℝ => (-1 : ℝ) ^ k * x ^ k by
      funext x
      rw [neg_pow],
      intervalIntegral.integral_const_mul, integral_pow]
    ring
  have hsign :
      (∑ k ∈ Finset.range (n + 1),
        (-1 : ℝ) ^ k / (k + 1 : ℝ)) =
        -signedHarmonic24 (n + 1) := by
    unfold signedHarmonic24
    rw [← Finset.sum_neg_distrib]
    apply Finset.sum_congr rfl
    intro k _
    rw [pow_succ]
    ring
  unfold barHarmonicCubicMoment24
  rw [intervalIntegral.integral_const_mul,
    intervalIntegral.integral_finset_sum]
  · simp_rw [hterm]
    rw [hsign]
    unfold barHarmonicCubicTerm24
    ring
  · intro k _
    exact (continuous_id.neg.pow k).intervalIntegrable 0 1

private theorem barHarmonicCubicMoment24_nonneg
    (n : ℕ) {x : ℝ} (hx0 : 0 ≤ x) (hx1 : x ≤ 1) :
    0 ≤ barHarmonicCubicMoment24 n x := by
  have hxabs : |-x| ≤ 1 := by
    rw [abs_neg, abs_of_nonneg hx0]
    exact hx1
  have habspow : |(-x) ^ (n + 1)| ≤ 1 := by
    rw [abs_pow]
    exact pow_le_one₀ (abs_nonneg (-x)) hxabs
  have hpow : (-x) ^ (n + 1) ≤ 1 :=
    (le_abs_self ((-x) ^ (n + 1))).trans habspow
  have hgeom := geom_sum_mul_neg (-x) (n + 1)
  have hden : 0 < 1 + x := by linarith
  have hsum :
      0 ≤ ∑ k ∈ Finset.range (n + 1), (-x) ^ k := by
    apply nonneg_of_mul_nonneg_left
      (b := (1 + x : ℝ)) (by
        rw [show 1 + x = 1 - (-x) by ring, hgeom]
        linarith)
    exact hden
  unfold barHarmonicCubicMoment24
  positivity

private theorem barHarmonicCubicMoment24_integral_norm
    (n : ℕ) :
    (∫ x : ℝ in 0..1, ‖barHarmonicCubicMoment24 n x‖) =
      barHarmonicCubicTerm24 n := by
  calc
    (∫ x : ℝ in 0..1, ‖barHarmonicCubicMoment24 n x‖) =
        ∫ x : ℝ in 0..1, barHarmonicCubicMoment24 n x := by
      apply intervalIntegral.integral_congr
      intro x hx
      simp only [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1),
        Set.mem_Icc] at hx
      change
        ‖barHarmonicCubicMoment24 n x‖ =
          barHarmonicCubicMoment24 n x
      rw [Real.norm_eq_abs,
        abs_of_nonneg
          (barHarmonicCubicMoment24_nonneg n hx.1 hx.2)]
    _ = barHarmonicCubicTerm24 n :=
      barHarmonicCubicMoment24_integral n

private theorem barHarmonicCubicMoment24_integral_norm_summable :
    Summable
      (fun n : ℕ =>
        ∫ x : ℝ in 0..1,
          ‖barHarmonicCubicMoment24 n x‖) := by
  exact summable_barHarmonicCubicTerm24.congr fun n =>
    (barHarmonicCubicMoment24_integral_norm n).symm

private def barHarmonicCubicKernel24 (x : ℝ) : ℝ :=
  (zeta3_24 - RamanujanChallenge.P26.trilog26 (-x)) / (1 + x)

private theorem barHarmonicCubicMoment24_hasSum_pointwise
    {x : ℝ} (hx0 : 0 ≤ x) (hx1 : x < 1) :
    HasSum (fun n : ℕ => barHarmonicCubicMoment24 n x)
      (barHarmonicCubicKernel24 x) := by
  have hden : 1 + x ≠ 0 := by linarith
  have hdiff :=
    shifted_zeta_three_hasSum24.sub
      (trilog26_hasSum24 (x := -x)
        (by rw [abs_neg, abs_of_nonneg hx0]; exact hx1.le))
  have hscaled := hdiff.mul_right (1 / (1 + x))
  convert hscaled using 1
  · funext n
    unfold barHarmonicCubicMoment24
    have hgeom := geom_sum_mul_neg (-x) (n + 1)
    rw [show 1 - -x = 1 + x by ring] at hgeom
    field_simp [hden]
    rw [hgeom]
  · unfold barHarmonicCubicKernel24
    ring

private theorem barHarmonicCubicTerm24_hasSum_integral :
    HasSum barHarmonicCubicTerm24
      (∫ x : ℝ in 0..1, barHarmonicCubicKernel24 x) := by
  have hInt :
      ∀ n : ℕ,
        MeasureTheory.Integrable
          (barHarmonicCubicMoment24 n)
          (MeasureTheory.volume.restrict (Set.Ioc 0 1)) := by
    intro n
    exact (barHarmonicCubicMoment24_intervalIntegrable n).1
  have hNorm :
      Summable
        (fun n : ℕ =>
          ∫ x : ℝ in Set.Ioc 0 1,
            ‖barHarmonicCubicMoment24 n x‖) := by
    simpa only [
      ← intervalIntegral.integral_of_le
        (by norm_num : (0 : ℝ) ≤ 1)] using
      barHarmonicCubicMoment24_integral_norm_summable
  have h :=
    MeasureTheory.hasSum_integral_of_summable_integral_norm
      (μ := MeasureTheory.volume.restrict (Set.Ioc (0 : ℝ) 1))
      hInt hNorm
  have h' :
      HasSum barHarmonicCubicTerm24
        (∫ x : ℝ in Set.Ioc 0 1,
          ∑' n : ℕ, barHarmonicCubicMoment24 n x) := by
    convert h using 1
    funext n
    rw [← intervalIntegral.integral_of_le
      (by norm_num : (0 : ℝ) ≤ 1)]
    exact (barHarmonicCubicMoment24_integral n).symm
  convert h' using 1
  rw [intervalIntegral.integral_of_le (by norm_num : (0 : ℝ) ≤ 1)]
  apply MeasureTheory.setIntegral_congr_ae measurableSet_Ioc
  filter_upwards [
    MeasureTheory.Measure.ae_ne MeasureTheory.volume (1 : ℝ)] with x hxne hx
  have hxlt : x < 1 := lt_of_le_of_ne hx.2 hxne
  exact
    (barHarmonicCubicMoment24_hasSum_pointwise
      hx.1.le hxlt).tsum_eq.symm

private def barHarmonicCubicPrimitive24 (x : ℝ) : ℝ :=
  (zeta3_24 - RamanujanChallenge.P26.trilog26 (-x)) *
      Real.log (1 + x) -
    (1 / 2 : ℝ) * dilog (-x) ^ 2

private theorem barHarmonicCubicPrimitive24_hasDerivAt
    {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasDerivAt barHarmonicCubicPrimitive24
      (barHarmonicCubicKernel24 x) x := by
  have hxne : x ≠ 0 := ne_of_gt hx0
  have hnegne : -x ≠ 0 := neg_ne_zero.mpr hxne
  have hxabs : |-x| < 1 := by
    rw [abs_neg, abs_of_pos hx0]
    exact hx1
  have hneg :
      HasDerivAt (fun y : ℝ => -y) (-1) x := by
    simpa using (hasDerivAt_id x).neg
  have ht :
      HasDerivAt
        (fun y : ℝ => RamanujanChallenge.P26.trilog26 (-y))
        (dilog (-x) / x) x := by
    convert
      (RamanujanChallenge.P26.trilog26_hasDerivAt_of_abs_lt_one
        hxabs hnegne).comp x hneg using 1
    field_simp [hxne]
  have hd :
      HasDerivAt (fun y : ℝ => dilog (-y))
        (-Real.log (1 + x) / x) x := by
    convert
      (dilog_hasDerivAt_of_abs_lt_one hxabs hnegne).comp x hneg using 1
    rw [show 1 - -x = 1 + x by ring]
    field_simp [hxne]
  have hinner :
      HasDerivAt (fun y : ℝ => 1 + y) 1 x := by
    convert (hasDerivAt_const x 1).add (hasDerivAt_id x) using 1
    simp
  have h1pxne : 1 + x ≠ 0 := by linarith
  have hlog :
      HasDerivAt (fun y : ℝ => Real.log (1 + y))
        (1 / (1 + x)) x := by
    convert hinner.log h1pxne using 1
  have hu :=
    (hasDerivAt_const x zeta3_24).sub ht
  unfold barHarmonicCubicPrimitive24
  have htotal :=
    (hu.mul hlog).sub ((hd.pow 2).const_mul (1 / 2 : ℝ))
  convert htotal using 1
  unfold barHarmonicCubicKernel24
  simp only [Pi.sub_apply, Pi.mul_apply, Pi.pow_apply]
  field_simp [hxne, h1pxne]
  ring

private theorem barHarmonicCubicPrimitive24_continuousOn :
    ContinuousOn barHarmonicCubicPrimitive24 (Icc (0 : ℝ) 1) := by
  have hneg :
      ContinuousOn (fun x : ℝ => -x) (Icc (0 : ℝ) 1) := by
    fun_prop
  have hnegmem :
      MapsTo (fun x : ℝ => -x)
        (Icc (0 : ℝ) 1) (Icc (-1 : ℝ) 1) := by
    intro x hx
    constructor <;> linarith [hx.1, hx.2]
  have ht :
      ContinuousOn
        (fun x : ℝ => RamanujanChallenge.P26.trilog26 (-x))
        (Icc (0 : ℝ) 1) :=
    RamanujanChallenge.P26.trilog26_continuousOn_unit.comp
      hneg hnegmem
  have hd :
      ContinuousOn (fun x : ℝ => dilog (-x))
        (Icc (0 : ℝ) 1) :=
    dilog_continuousOn_unit.comp hneg hnegmem
  have hlog :
      ContinuousOn (fun x : ℝ => Real.log (1 + x))
        (Icc (0 : ℝ) 1) := by
    apply (continuousOn_const.add continuousOn_id).log
    intro x hx
    exact ne_of_gt (show 0 < 1 + x by linarith [hx.1])
  unfold barHarmonicCubicPrimitive24
  exact
    ((continuousOn_const.sub ht).mul hlog).sub
      (continuousOn_const.mul (hd.pow 2))

private theorem barHarmonicCubicKernel24_intervalIntegrable :
    IntervalIntegrable barHarmonicCubicKernel24
      MeasureTheory.volume 0 1 := by
  have ht :
      ContinuousOn
        (fun x : ℝ => RamanujanChallenge.P26.trilog26 (-x))
        (Icc (0 : ℝ) 1) := by
    apply RamanujanChallenge.P26.trilog26_continuousOn_unit.comp
    · fun_prop
    · intro x hx
      constructor <;> linarith [hx.1, hx.2]
  have hcont :
      ContinuousOn barHarmonicCubicKernel24
        (Icc (0 : ℝ) 1) := by
    unfold barHarmonicCubicKernel24
    apply ContinuousOn.div
      (continuousOn_const.sub ht)
      (continuousOn_const.add continuousOn_id)
    intro x hx
    exact ne_of_gt (show 0 < 1 + x by linarith [hx.1])
  apply ContinuousOn.intervalIntegrable
  rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
  exact hcont

private theorem barHarmonicCubicIntegral24 :
    (∫ x : ℝ in 0..1, barHarmonicCubicKernel24 x) =
      (7 / 4 : ℝ) * Real.log 2 * zeta3_24 -
        (1 / 8 : ℝ) * (Real.pi ^ 2 / 6) ^ 2 := by
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le
    (f := barHarmonicCubicPrimitive24)
    (f' := barHarmonicCubicKernel24)
    (by norm_num)
    barHarmonicCubicPrimitive24_continuousOn
    (fun x hx =>
      barHarmonicCubicPrimitive24_hasDerivAt hx.1 hx.2)
    barHarmonicCubicKernel24_intervalIntegrable]
  unfold barHarmonicCubicPrimitive24
  simp only [neg_zero]
  rw [RamanujanChallenge.P26.trilog26_neg_one,
    RamanujanChallenge.P26.trilog26_zero,
    RamanujanChallenge.P26.dilog_neg_one26,
    dilog_zero]
  norm_num
  unfold RamanujanChallenge.P26.zeta3 zeta3_24
  ring

private theorem trilog26_hasDerivWithinAt_one24 :
    HasDerivWithinAt RamanujanChallenge.P26.trilog26
      (Real.pi ^ 2 / 6) (Iic (1 : ℝ)) 1 := by
  apply hasDerivWithinAt_Iic_of_tendsto_deriv
      (s := Ioo (0 : ℝ) 1)
  · intro x hx
    exact
      (RamanujanChallenge.P26.trilog26_hasDerivAt_of_abs_lt_one
        (by rw [abs_of_pos hx.1]; exact hx.2) (ne_of_gt hx.1)).differentiableAt
        |>.differentiableWithinAt
  · exact
      (RamanujanChallenge.P26.trilog26_continuousOn_unit
        1 (by norm_num)).mono (by
          intro x hx
          constructor <;> linarith [hx.1, hx.2])
  · exact Ioo_mem_nhdsLT (by norm_num : (0 : ℝ) < 1)
  · have hdWithin :
        ContinuousWithinAt dilog (Iio (1 : ℝ)) 1 :=
      (dilog_continuousOn_unit 1 (by norm_num)).mono_of_mem_nhdsWithin
        (Icc_mem_nhdsLT (show (-1 : ℝ) < 1 by norm_num))
    have hd :
        Tendsto (fun x : ℝ => dilog x / x)
          (𝓝[<] (1 : ℝ)) (𝓝 (Real.pi ^ 2 / 6)) := by
      have hid :
          Tendsto (fun x : ℝ => x)
            (𝓝[<] (1 : ℝ)) (𝓝 1) :=
        tendsto_id.mono_left inf_le_left
      simpa [dilog_one] using
        hdWithin.tendsto.div hid (by norm_num : (1 : ℝ) ≠ 0)
    apply hd.congr'
    filter_upwards [Ioo_mem_nhdsLT (by norm_num : (0 : ℝ) < 1)]
      with x hx
    exact
      (RamanujanChallenge.P26.trilog26_hasDerivAt_of_abs_lt_one
        (by rw [abs_of_pos hx.1]; exact hx.2) (ne_of_gt hx.1)).deriv.symm

private theorem ordinaryHarmonicCubicKernel24_tendsto_one :
    Tendsto ordinaryHarmonicCubicKernel24
      (𝓝[<] (1 : ℝ)) (𝓝 (Real.pi ^ 2 / 6)) := by
  have hderiv :=
    trilog26_hasDerivWithinAt_one24.Iio_of_Iic
  have hslope :=
    (hasDerivWithinAt_iff_tendsto_slope'
      (show (1 : ℝ) ∉ Iio 1 by simp)).1 hderiv
  apply hslope.congr'
  filter_upwards [self_mem_nhdsWithin] with x hx
  have hxne : x ≠ 1 := ne_of_lt hx
  unfold ordinaryHarmonicCubicKernel24
  simp only [slope, vsub_eq_sub, smul_eq_mul]
  rw [RamanujanChallenge.P26.trilog26_one]
  unfold RamanujanChallenge.P26.zeta3 zeta3_24
  have hxsub : x - 1 ≠ 0 := sub_ne_zero.mpr hxne
  have hsubx : 1 - x ≠ 0 := sub_ne_zero.mpr hxne.symm
  field_simp [hxsub, hsubx]
  ring

private def ordinaryHarmonicCubicKernelExtension24 (x : ℝ) : ℝ :=
  Function.update ordinaryHarmonicCubicKernel24 1
    (Real.pi ^ 2 / 6) x

private theorem ordinaryHarmonicCubicKernelExtension24_continuousOn :
    ContinuousOn ordinaryHarmonicCubicKernelExtension24
      (Icc (1 / 2 : ℝ) 1) := by
  intro x hx
  by_cases hxone : x = 1
  · subst x
    unfold ordinaryHarmonicCubicKernelExtension24
    rw [continuousWithinAt_update_same]
    have hmono :
        𝓝[Icc (1 / 2 : ℝ) 1 \ {1}] (1 : ℝ) ≤
          𝓝[<] (1 : ℝ) := by
      apply nhdsWithin_mono
      intro y hy
      exact lt_of_le_of_ne hy.1.2 hy.2
    have ht :=
      ordinaryHarmonicCubicKernel24_tendsto_one.mono_left hmono
    apply ht.congr'
    filter_upwards [self_mem_nhdsWithin] with y hy
    simp [Function.update, hy.2]
  · unfold ordinaryHarmonicCubicKernelExtension24
    rw [continuousWithinAt_update_of_ne hxone]
    have htri :
        ContinuousAt RamanujanChallenge.P26.trilog26 x :=
      RamanujanChallenge.P26.trilog26_continuousOn_unit.continuousAt
        (Icc_mem_nhds (by linarith [hx.1]) (by
          exact lt_of_le_of_ne hx.2 hxone))
    unfold ordinaryHarmonicCubicKernel24
    exact
      ((continuousAt_const.sub htri).div
        (continuousAt_const.sub continuousAt_id)
        (sub_ne_zero.mpr (Ne.symm hxone))).continuousWithinAt

private theorem ordinaryHarmonicCubicKernel24_intervalIntegrable :
    IntervalIntegrable ordinaryHarmonicCubicKernel24
      MeasureTheory.volume 0 1 := by
  apply IntervalIntegrable.trans (b := (1 / 2 : ℝ))
  · have htri :
        ContinuousOn RamanujanChallenge.P26.trilog26
          (Icc (0 : ℝ) (1 / 2)) :=
      RamanujanChallenge.P26.trilog26_continuousOn_unit.mono (by
        intro x hx
        constructor <;> linarith [hx.1, hx.2])
    have hcont :
        ContinuousOn ordinaryHarmonicCubicKernel24
          (Icc (0 : ℝ) (1 / 2)) := by
      unfold ordinaryHarmonicCubicKernel24
      apply ContinuousOn.div
        (continuousOn_const.sub htri)
        (continuousOn_const.sub continuousOn_id)
      intro x hx
      exact ne_of_gt (show 0 < 1 - x by linarith [hx.2])
    apply ContinuousOn.intervalIntegrable
    rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2)]
    exact hcont
  · have hint :
        IntervalIntegrable ordinaryHarmonicCubicKernelExtension24
          MeasureTheory.volume (1 / 2) 1 := by
      apply ContinuousOn.intervalIntegrable
      rw [Set.uIcc_of_le (by norm_num : (1 / 2 : ℝ) ≤ 1)]
      exact ordinaryHarmonicCubicKernelExtension24_continuousOn
    apply hint.congr_ae
    filter_upwards [
      MeasureTheory.ae_restrict_of_ae
        (MeasureTheory.Measure.ae_ne MeasureTheory.volume (1 : ℝ))
    ] with x hxne
    simp [ordinaryHarmonicCubicKernelExtension24, hxne]

private def ordinaryHarmonicCubicPrimitive24 (x : ℝ) : ℝ :=
  -((zeta3_24 - RamanujanChallenge.P26.trilog26 x) *
      Real.log (1 - x)) +
    (1 / 2 : ℝ) * dilog x ^ 2

private theorem ordinaryHarmonicCubicPrimitive24_hasDerivAt
    {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasDerivAt ordinaryHarmonicCubicPrimitive24
      (ordinaryHarmonicCubicKernel24 x) x := by
  have hxne : x ≠ 0 := ne_of_gt hx0
  have h1xne : 1 - x ≠ 0 := ne_of_gt (sub_pos.mpr hx1)
  have ht :=
    RamanujanChallenge.P26.trilog26_hasDerivAt_of_abs_lt_one
      (by rw [abs_of_pos hx0]; exact hx1) hxne
  have hd := dilog_hasDerivAt hx0 hx1
  have hsub :
      HasDerivAt (fun y : ℝ => 1 - y) (-1) x := by
    convert (hasDerivAt_const x 1).sub (hasDerivAt_id x) using 1
    simp
  have hlog :
      HasDerivAt (fun y : ℝ => Real.log (1 - y))
        (-1 / (1 - x)) x := by
    convert hsub.log h1xne using 1
  have hu :=
    (hasDerivAt_const x zeta3_24).sub ht
  unfold ordinaryHarmonicCubicPrimitive24
  have htotal :=
    (hu.mul hlog).neg.add ((hd.pow 2).const_mul (1 / 2 : ℝ))
  convert htotal using 1
  unfold ordinaryHarmonicCubicKernel24
  simp only [Pi.sub_apply]
  field_simp [hxne, h1xne]
  ring

private theorem ordinaryHarmonicCubicPrimitive24_tendsto_zero :
    Tendsto ordinaryHarmonicCubicPrimitive24
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  have ht :
      ContinuousAt RamanujanChallenge.P26.trilog26 0 :=
    RamanujanChallenge.P26.trilog26_continuousOn_unit.continuousAt
      (Icc_mem_nhds (by norm_num : (-1 : ℝ) < 0)
        (by norm_num : (0 : ℝ) < 1))
  have hd :
      ContinuousAt dilog 0 :=
    dilog_continuousOn_unit.continuousAt
      (Icc_mem_nhds (by norm_num : (-1 : ℝ) < 0)
        (by norm_num : (0 : ℝ) < 1))
  have hlog :
      ContinuousAt (fun x : ℝ => Real.log (1 - x)) 0 := by
    apply (continuousAt_const.sub continuousAt_id).log
    norm_num
  have hcont :
      ContinuousAt ordinaryHarmonicCubicPrimitive24 0 := by
    unfold ordinaryHarmonicCubicPrimitive24
    exact
      ((continuousAt_const.sub ht).mul hlog).neg.add
        (continuousAt_const.mul (hd.pow 2))
  simpa [ordinaryHarmonicCubicPrimitive24,
    RamanujanChallenge.P26.trilog26_zero, dilog_zero] using
    tendsto_nhdsWithin_of_tendsto_nhds hcont.tendsto

private theorem ordinaryHarmonicCubicPrimitive24_tendsto_one :
    Tendsto ordinaryHarmonicCubicPrimitive24
      (𝓝[<] (1 : ℝ))
      (𝓝 ((1 / 2 : ℝ) * (Real.pi ^ 2 / 6) ^ 2)) := by
  have hsub :
      Tendsto (fun x : ℝ => 1 - x)
        (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    have hcont :
        ContinuousAt (fun x : ℝ => 1 - x) 1 :=
      continuousAt_const.sub continuousAt_id
    simpa using
      (hcont.tendsto.mono_left inf_le_left)
  have hmulLog :
      Tendsto (fun x : ℝ => (1 - x) * Real.log (1 - x))
        (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    simpa using
      Real.continuous_mul_log.continuousAt.tendsto.comp hsub
  have hrawProduct :
      Tendsto
        (fun x : ℝ =>
          ordinaryHarmonicCubicKernel24 x *
            ((1 - x) * Real.log (1 - x)))
        (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    simpa using
      ordinaryHarmonicCubicKernel24_tendsto_one.mul hmulLog
  have hproduct :
      Tendsto
        (fun x : ℝ =>
          (zeta3_24 - RamanujanChallenge.P26.trilog26 x) *
            Real.log (1 - x))
        (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    apply hrawProduct.congr'
    filter_upwards [self_mem_nhdsWithin] with x hx
    have hxne : 1 - x ≠ 0 := ne_of_gt (sub_pos.mpr hx)
    unfold ordinaryHarmonicCubicKernel24
    field_simp [hxne]
  have hdWithin :
      ContinuousWithinAt dilog (Iio (1 : ℝ)) 1 :=
    (dilog_continuousOn_unit 1 (by norm_num)).mono_of_mem_nhdsWithin
      (Icc_mem_nhdsLT (show (-1 : ℝ) < 1 by norm_num))
  have hsquare :
      Tendsto (fun x : ℝ => (1 / 2 : ℝ) * dilog x ^ 2)
        (𝓝[<] (1 : ℝ))
        (𝓝 ((1 / 2 : ℝ) * (Real.pi ^ 2 / 6) ^ 2)) := by
    simpa [dilog_one] using
      (hdWithin.tendsto.pow 2).const_mul (1 / 2 : ℝ)
  unfold ordinaryHarmonicCubicPrimitive24
  simpa using hproduct.neg.add hsquare

private theorem ordinaryHarmonicCubicIntegral24 :
    (∫ x : ℝ in 0..1, ordinaryHarmonicCubicKernel24 x) =
      (1 / 2 : ℝ) * (Real.pi ^ 2 / 6) ^ 2 := by
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_tendsto
    (f := ordinaryHarmonicCubicPrimitive24)
    (fa := (0 : ℝ))
    (fb := (1 / 2 : ℝ) * (Real.pi ^ 2 / 6) ^ 2)
    (by norm_num)
    (fun x hx =>
      ordinaryHarmonicCubicPrimitive24_hasDerivAt hx.1 hx.2)
    ordinaryHarmonicCubicKernel24_intervalIntegrable
    ordinaryHarmonicCubicPrimitive24_tendsto_zero
    ordinaryHarmonicCubicPrimitive24_tendsto_one]
  ring

theorem cubicLinearEulerTerm24_hasSum :
    HasSum cubicLinearEulerTerm24 cubicLinearEulerValue24 := by
  have hOrd :
      HasSum ordinaryHarmonicCubicTerm24
        ((1 / 2 : ℝ) * (Real.pi ^ 2 / 6) ^ 2) := by
    rw [← ordinaryHarmonicCubicIntegral24]
    exact ordinaryHarmonicCubicTerm24_hasSum_integral
  have hBar :
      HasSum barHarmonicCubicTerm24
        ((7 / 4 : ℝ) * Real.log 2 * zeta3_24 -
          (1 / 8 : ℝ) * (Real.pi ^ 2 / 6) ^ 2) := by
    rw [← barHarmonicCubicIntegral24]
    exact barHarmonicCubicTerm24_hasSum_integral
  have hCombined := hOrd.sub (hBar.mul_left 2)
  convert hCombined using 1
  · funext n
    unfold cubicLinearEulerTerm24 parityRemainder24
      ordinaryHarmonicCubicTerm24 barHarmonicCubicTerm24
    ring
  · unfold cubicLinearEulerValue24
    ring

private def alternatingOrdinaryHarmonicCubicTerm24
    (n : ℕ) : ℝ :=
  (-1 : ℝ) ^ (n + 1) * ordinaryHarmonicCubicTerm24 n

private def alternatingOrdinaryHarmonicCubicMoment24
    (n : ℕ) (x : ℝ) : ℝ :=
  (-1 : ℝ) ^ (n + 1) * ordinaryHarmonicCubicMoment24 n x

private def alternatingOrdinaryHarmonicCubicKernel24
    (x : ℝ) : ℝ :=
  (RamanujanChallenge.P26.trilog26 (-1) -
      RamanujanChallenge.P26.trilog26 (-x)) /
    (1 - x)

private theorem
    alternatingOrdinaryHarmonicCubicMoment24_intervalIntegrable
    (n : ℕ) :
    IntervalIntegrable
      (alternatingOrdinaryHarmonicCubicMoment24 n)
      MeasureTheory.volume 0 1 := by
  unfold alternatingOrdinaryHarmonicCubicMoment24
  exact
    (ordinaryHarmonicCubicMoment24_intervalIntegrable n).const_mul
      ((-1 : ℝ) ^ (n + 1))

private theorem alternatingOrdinaryHarmonicCubicMoment24_integral
    (n : ℕ) :
    (∫ x : ℝ in 0..1,
      alternatingOrdinaryHarmonicCubicMoment24 n x) =
      alternatingOrdinaryHarmonicCubicTerm24 n := by
  unfold alternatingOrdinaryHarmonicCubicMoment24
    alternatingOrdinaryHarmonicCubicTerm24
  rw [intervalIntegral.integral_const_mul,
    ordinaryHarmonicCubicMoment24_integral]

private theorem
    alternatingOrdinaryHarmonicCubicMoment24_integral_norm
    (n : ℕ) :
    (∫ x : ℝ in 0..1,
      ‖alternatingOrdinaryHarmonicCubicMoment24 n x‖) =
      ordinaryHarmonicCubicTerm24 n := by
  calc
    (∫ x : ℝ in 0..1,
        ‖alternatingOrdinaryHarmonicCubicMoment24 n x‖) =
        ∫ x : ℝ in 0..1,
          ‖ordinaryHarmonicCubicMoment24 n x‖ := by
      apply intervalIntegral.integral_congr
      intro x _
      simp [alternatingOrdinaryHarmonicCubicMoment24]
    _ = ordinaryHarmonicCubicTerm24 n :=
      ordinaryHarmonicCubicMoment24_integral_norm n

private theorem
    alternatingOrdinaryHarmonicCubicMoment24_integral_norm_summable :
    Summable
      (fun n : ℕ =>
        ∫ x : ℝ in 0..1,
          ‖alternatingOrdinaryHarmonicCubicMoment24 n x‖) := by
  exact summable_ordinaryHarmonicCubicTerm24.congr fun n =>
    (alternatingOrdinaryHarmonicCubicMoment24_integral_norm n).symm

private theorem
    alternatingOrdinaryHarmonicCubicMoment24_hasSum_pointwise
    {x : ℝ} (hx0 : 0 ≤ x) (hx1 : x < 1) :
    HasSum
      (fun n : ℕ =>
        alternatingOrdinaryHarmonicCubicMoment24 n x)
      (alternatingOrdinaryHarmonicCubicKernel24 x) := by
  have hxne : x ≠ 1 := ne_of_lt hx1
  have hdiff :=
    (trilog26_hasSum24 (x := (-1 : ℝ)) (by norm_num)).sub
      (trilog26_hasSum24 (x := -x)
        (by rw [abs_neg, abs_of_nonneg hx0]; exact hx1.le))
  have hscaled := hdiff.mul_right (1 / (1 - x))
  convert hscaled using 1
  · funext n
    unfold alternatingOrdinaryHarmonicCubicMoment24
      ordinaryHarmonicCubicMoment24
    rw [geom_sum_eq hxne]
    have hxsub : x - 1 ≠ 0 := sub_ne_zero.mpr hxne
    have hsubx : 1 - x ≠ 0 := sub_ne_zero.mpr hxne.symm
    field_simp [hxsub, hsubx]
    ring
  · unfold alternatingOrdinaryHarmonicCubicKernel24
    ring

private theorem
    alternatingOrdinaryHarmonicCubicTerm24_hasSum_integral :
    HasSum alternatingOrdinaryHarmonicCubicTerm24
      (∫ x : ℝ in 0..1,
        alternatingOrdinaryHarmonicCubicKernel24 x) := by
  have hInt :
      ∀ n : ℕ,
        MeasureTheory.Integrable
          (alternatingOrdinaryHarmonicCubicMoment24 n)
          (MeasureTheory.volume.restrict (Set.Ioc 0 1)) := by
    intro n
    exact
      (alternatingOrdinaryHarmonicCubicMoment24_intervalIntegrable n).1
  have hNorm :
      Summable
        (fun n : ℕ =>
          ∫ x : ℝ in Set.Ioc 0 1,
            ‖alternatingOrdinaryHarmonicCubicMoment24 n x‖) := by
    simpa only [
      ← intervalIntegral.integral_of_le
        (by norm_num : (0 : ℝ) ≤ 1)] using
      alternatingOrdinaryHarmonicCubicMoment24_integral_norm_summable
  have h :=
    MeasureTheory.hasSum_integral_of_summable_integral_norm
      (μ := MeasureTheory.volume.restrict (Set.Ioc (0 : ℝ) 1))
      hInt hNorm
  have h' :
      HasSum alternatingOrdinaryHarmonicCubicTerm24
        (∫ x : ℝ in Set.Ioc 0 1,
          ∑' n : ℕ,
            alternatingOrdinaryHarmonicCubicMoment24 n x) := by
    convert h using 1
    funext n
    rw [← intervalIntegral.integral_of_le
      (by norm_num : (0 : ℝ) ≤ 1)]
    exact
      (alternatingOrdinaryHarmonicCubicMoment24_integral n).symm
  convert h' using 1
  rw [intervalIntegral.integral_of_le (by norm_num : (0 : ℝ) ≤ 1)]
  apply MeasureTheory.setIntegral_congr_ae measurableSet_Ioc
  filter_upwards [
    MeasureTheory.Measure.ae_ne MeasureTheory.volume (1 : ℝ)
  ] with x hxne hx
  have hxlt : x < 1 := lt_of_le_of_ne hx.2 hxne
  exact
    (alternatingOrdinaryHarmonicCubicMoment24_hasSum_pointwise
      hx.1.le hxlt).tsum_eq.symm

private def alternatingBarHarmonicCubicTerm24
    (n : ℕ) : ℝ :=
  (-1 : ℝ) ^ (n + 1) * barHarmonicCubicTerm24 n

private def alternatingBarHarmonicCubicMoment24
    (n : ℕ) (x : ℝ) : ℝ :=
  (-1 : ℝ) ^ (n + 1) * barHarmonicCubicMoment24 n x

private def alternatingBarHarmonicCubicKernel24
    (x : ℝ) : ℝ :=
  (RamanujanChallenge.P26.trilog26 (-1) -
      RamanujanChallenge.P26.trilog26 x) /
    (1 + x)

private theorem
    alternatingBarHarmonicCubicMoment24_intervalIntegrable
    (n : ℕ) :
    IntervalIntegrable (alternatingBarHarmonicCubicMoment24 n)
      MeasureTheory.volume 0 1 := by
  unfold alternatingBarHarmonicCubicMoment24
  exact
    (barHarmonicCubicMoment24_intervalIntegrable n).const_mul
      ((-1 : ℝ) ^ (n + 1))

private theorem alternatingBarHarmonicCubicMoment24_integral
    (n : ℕ) :
    (∫ x : ℝ in 0..1,
      alternatingBarHarmonicCubicMoment24 n x) =
      alternatingBarHarmonicCubicTerm24 n := by
  unfold alternatingBarHarmonicCubicMoment24
    alternatingBarHarmonicCubicTerm24
  rw [intervalIntegral.integral_const_mul,
    barHarmonicCubicMoment24_integral]

private theorem alternatingBarHarmonicCubicMoment24_integral_norm
    (n : ℕ) :
    (∫ x : ℝ in 0..1,
      ‖alternatingBarHarmonicCubicMoment24 n x‖) =
      barHarmonicCubicTerm24 n := by
  calc
    (∫ x : ℝ in 0..1,
        ‖alternatingBarHarmonicCubicMoment24 n x‖) =
        ∫ x : ℝ in 0..1,
          ‖barHarmonicCubicMoment24 n x‖ := by
      apply intervalIntegral.integral_congr
      intro x _
      simp [alternatingBarHarmonicCubicMoment24]
    _ = barHarmonicCubicTerm24 n :=
      barHarmonicCubicMoment24_integral_norm n

private theorem
    alternatingBarHarmonicCubicMoment24_integral_norm_summable :
    Summable
      (fun n : ℕ =>
        ∫ x : ℝ in 0..1,
          ‖alternatingBarHarmonicCubicMoment24 n x‖) := by
  exact summable_barHarmonicCubicTerm24.congr fun n =>
    (alternatingBarHarmonicCubicMoment24_integral_norm n).symm

private theorem alternatingBarHarmonicCubicMoment24_hasSum_pointwise
    {x : ℝ} (hx0 : 0 ≤ x) (hx1 : x < 1) :
    HasSum
      (fun n : ℕ => alternatingBarHarmonicCubicMoment24 n x)
      (alternatingBarHarmonicCubicKernel24 x) := by
  have hden : 1 + x ≠ 0 := by linarith
  have hdiff :=
    (trilog26_hasSum24 (x := (-1 : ℝ)) (by norm_num)).sub
      (trilog26_hasSum24 (x := x)
        (by rw [abs_of_nonneg hx0]; exact hx1.le))
  have hscaled := hdiff.mul_right (1 / (1 + x))
  convert hscaled using 1
  · funext n
    unfold alternatingBarHarmonicCubicMoment24
      barHarmonicCubicMoment24
    have hgeom := geom_sum_mul_neg (-x) (n + 1)
    rw [show 1 - -x = 1 + x by ring] at hgeom
    field_simp [hden]
    calc
      ((-1 : ℝ) ^ (n + 1) *
          ∑ k ∈ Finset.range (n + 1), (-x) ^ k) *
          (1 + x) =
          (-1 : ℝ) ^ (n + 1) *
            ((∑ k ∈ Finset.range (n + 1), (-x) ^ k) *
              (1 + x)) := by ring
      _ = (-1 : ℝ) ^ (n + 1) *
          (1 - (-x) ^ (n + 1)) := by rw [hgeom]
      _ = (-1 : ℝ) ^ (n + 1) - x ^ (n + 1) := by
        rw [neg_pow x (n + 1)]
        have hsignsq :
            (-1 : ℝ) ^ (n + 1) * (-1 : ℝ) ^ (n + 1) =
              1 := by
          rw [← pow_add]
          simp [← two_mul]
        rw [mul_sub]
        rw [mul_one, ← mul_assoc, hsignsq, one_mul]
  · unfold alternatingBarHarmonicCubicKernel24
    ring

private theorem alternatingBarHarmonicCubicTerm24_hasSum_integral :
    HasSum alternatingBarHarmonicCubicTerm24
      (∫ x : ℝ in 0..1,
        alternatingBarHarmonicCubicKernel24 x) := by
  have hInt :
      ∀ n : ℕ,
        MeasureTheory.Integrable
          (alternatingBarHarmonicCubicMoment24 n)
          (MeasureTheory.volume.restrict (Set.Ioc 0 1)) := by
    intro n
    exact
      (alternatingBarHarmonicCubicMoment24_intervalIntegrable n).1
  have hNorm :
      Summable
        (fun n : ℕ =>
          ∫ x : ℝ in Set.Ioc 0 1,
            ‖alternatingBarHarmonicCubicMoment24 n x‖) := by
    simpa only [
      ← intervalIntegral.integral_of_le
        (by norm_num : (0 : ℝ) ≤ 1)] using
      alternatingBarHarmonicCubicMoment24_integral_norm_summable
  have h :=
    MeasureTheory.hasSum_integral_of_summable_integral_norm
      (μ := MeasureTheory.volume.restrict (Set.Ioc (0 : ℝ) 1))
      hInt hNorm
  have h' :
      HasSum alternatingBarHarmonicCubicTerm24
        (∫ x : ℝ in Set.Ioc 0 1,
          ∑' n : ℕ, alternatingBarHarmonicCubicMoment24 n x) := by
    convert h using 1
    funext n
    rw [← intervalIntegral.integral_of_le
      (by norm_num : (0 : ℝ) ≤ 1)]
    exact (alternatingBarHarmonicCubicMoment24_integral n).symm
  convert h' using 1
  rw [intervalIntegral.integral_of_le (by norm_num : (0 : ℝ) ≤ 1)]
  apply MeasureTheory.setIntegral_congr_ae measurableSet_Ioc
  filter_upwards [
    MeasureTheory.Measure.ae_ne MeasureTheory.volume (1 : ℝ)
  ] with x hxne hx
  have hxlt : x < 1 := lt_of_le_of_ne hx.2 hxne
  exact
    (alternatingBarHarmonicCubicMoment24_hasSum_pointwise
      hx.1.le hxlt).tsum_eq.symm

private theorem trilog26_neg_comp_hasDerivWithinAt_one24 :
    HasDerivWithinAt
      (fun x : ℝ => RamanujanChallenge.P26.trilog26 (-x))
      (-(Real.pi ^ 2 / 12)) (Iic (1 : ℝ)) 1 := by
  apply hasDerivWithinAt_Iic_of_tendsto_deriv
      (s := Ioo (0 : ℝ) 1)
  · intro x hx
    have hxne : x ≠ 0 := ne_of_gt hx.1
    have hneg :
        HasDerivAt (fun y : ℝ => -y) (-1) x := by
      simpa using (hasDerivAt_id x).neg
    have htri :=
      (RamanujanChallenge.P26.trilog26_hasDerivAt_of_abs_lt_one
        (x := -x)
        (by rw [abs_neg, abs_of_pos hx.1]; exact hx.2)
        (neg_ne_zero.mpr hxne)).comp x hneg
    exact htri.differentiableAt.differentiableWithinAt
  · have htri :
        ContinuousOn
          (fun x : ℝ => RamanujanChallenge.P26.trilog26 (-x))
          (Icc (0 : ℝ) 1) := by
      apply RamanujanChallenge.P26.trilog26_continuousOn_unit.comp
        continuousOn_id.neg
      intro x hx
      change -1 ≤ -x ∧ -x ≤ 1
      constructor <;> linarith [hx.1, hx.2]
    exact (htri 1 (by constructor <;> norm_num)).mono (by
      intro x hx
      exact ⟨hx.1.le, hx.2.le⟩)
  · exact Ioo_mem_nhdsLT (by norm_num : (0 : ℝ) < 1)
  · have hdWithin :
        ContinuousWithinAt (fun x : ℝ => dilog (-x))
          (Iio (1 : ℝ)) 1 := by
      have houter :
          ContinuousWithinAt dilog (Ici (-1 : ℝ)) (-1) :=
        (dilog_continuousOn_unit (-1)
          ⟨le_rfl, by norm_num⟩).mono_of_mem_nhdsWithin
          (Icc_mem_nhdsGE (by norm_num : (-1 : ℝ) < 1))
      have hinner :
          ContinuousWithinAt (fun x : ℝ => -x)
            (Iic (1 : ℝ)) 1 :=
        continuousAt_id.neg.continuousWithinAt
      have hcomp :=
        houter.comp hinner (by
          intro x hx
          change x ≤ 1 at hx
          change -1 ≤ -x
          linarith)
      simpa only [Function.comp_apply] using
        hcomp.mono Iio_subset_Iic_self
    have hd :
        Tendsto (fun x : ℝ => dilog (-x) / x)
          (𝓝[<] (1 : ℝ)) (𝓝 (-(Real.pi ^ 2 / 12))) := by
      have hid :
          Tendsto (fun x : ℝ => x)
            (𝓝[<] (1 : ℝ)) (𝓝 1) :=
        tendsto_id.mono_left inf_le_left
      have hq :=
        hdWithin.tendsto.div hid (by norm_num : (1 : ℝ) ≠ 0)
      rw [show -(Real.pi ^ 2 / 12) = dilog (-1) / 1 by
        rw [RamanujanChallenge.P26.dilog_neg_one26]
        ring]
      simpa only [Pi.div_apply] using hq
    apply hd.congr'
    filter_upwards [Ioo_mem_nhdsLT (by norm_num : (0 : ℝ) < 1)]
      with x hx
    have hxne : x ≠ 0 := ne_of_gt hx.1
    have hneg :
        HasDerivAt (fun y : ℝ => -y) (-1) x := by
      simpa using (hasDerivAt_id x).neg
    have htri :=
      (RamanujanChallenge.P26.trilog26_hasDerivAt_of_abs_lt_one
        (x := -x)
        (by rw [abs_neg, abs_of_pos hx.1]; exact hx.2)
        (neg_ne_zero.mpr hxne)).comp x hneg
    convert htri.deriv.symm using 1
    field_simp [hxne]

private theorem
    alternatingOrdinaryHarmonicCubicKernel24_tendsto_one :
    Tendsto alternatingOrdinaryHarmonicCubicKernel24
      (𝓝[<] (1 : ℝ)) (𝓝 (-(Real.pi ^ 2 / 12))) := by
  have hderiv :=
    trilog26_neg_comp_hasDerivWithinAt_one24.Iio_of_Iic
  have hslope :=
    (hasDerivWithinAt_iff_tendsto_slope'
      (show (1 : ℝ) ∉ Iio 1 by simp)).1 hderiv
  apply hslope.congr'
  filter_upwards [self_mem_nhdsWithin] with x hx
  have hxne : x ≠ 1 := ne_of_lt hx
  unfold alternatingOrdinaryHarmonicCubicKernel24
  simp only [slope, vsub_eq_sub, smul_eq_mul]
  have hxsub : x - 1 ≠ 0 := sub_ne_zero.mpr hxne
  have hsubx : 1 - x ≠ 0 := sub_ne_zero.mpr hxne.symm
  field_simp [hxsub, hsubx]
  ring

private def alternatingOrdinaryHarmonicCubicKernelExtension24
    (x : ℝ) : ℝ :=
  Function.update alternatingOrdinaryHarmonicCubicKernel24 1
    (-(Real.pi ^ 2 / 12)) x

private theorem
    alternatingOrdinaryHarmonicCubicKernelExtension24_continuousOn :
    ContinuousOn alternatingOrdinaryHarmonicCubicKernelExtension24
      (Icc (1 / 2 : ℝ) 1) := by
  intro x hx
  by_cases hxone : x = 1
  · subst x
    unfold alternatingOrdinaryHarmonicCubicKernelExtension24
    rw [continuousWithinAt_update_same]
    have hmono :
        𝓝[Icc (1 / 2 : ℝ) 1 \ {1}] (1 : ℝ) ≤
          𝓝[<] (1 : ℝ) := by
      apply nhdsWithin_mono
      intro y hy
      exact lt_of_le_of_ne hy.1.2 hy.2
    have ht :=
      alternatingOrdinaryHarmonicCubicKernel24_tendsto_one.mono_left
        hmono
    apply ht.congr'
    filter_upwards [self_mem_nhdsWithin] with y hy
    simp [hy.2]
  · unfold alternatingOrdinaryHarmonicCubicKernelExtension24
    rw [continuousWithinAt_update_of_ne hxone]
    have hxlt : x < 1 := lt_of_le_of_ne hx.2 hxone
    have htri :
        ContinuousAt
          (fun y : ℝ => RamanujanChallenge.P26.trilog26 (-y)) x := by
      have houter :
          ContinuousAt RamanujanChallenge.P26.trilog26 (-x) :=
        RamanujanChallenge.P26.trilog26_continuousOn_unit.continuousAt
          (Icc_mem_nhds (by linarith [hxlt])
            (by linarith [hx.1]))
      simpa only [Function.comp_apply] using
        houter.comp continuousAt_id.neg
    unfold alternatingOrdinaryHarmonicCubicKernel24
    exact
      ((continuousAt_const.sub htri).div
        (continuousAt_const.sub continuousAt_id)
        (sub_ne_zero.mpr (Ne.symm hxone))).continuousWithinAt

private theorem
    alternatingOrdinaryHarmonicCubicKernel24_intervalIntegrable :
    IntervalIntegrable alternatingOrdinaryHarmonicCubicKernel24
      MeasureTheory.volume 0 1 := by
  apply IntervalIntegrable.trans (b := (1 / 2 : ℝ))
  · have htri :
        ContinuousOn
          (fun x : ℝ => RamanujanChallenge.P26.trilog26 (-x))
          (Icc (0 : ℝ) (1 / 2)) := by
      apply RamanujanChallenge.P26.trilog26_continuousOn_unit.comp
        continuousOn_id.neg
      intro x hx
      change -1 ≤ -x ∧ -x ≤ 1
      constructor <;> linarith [hx.1, hx.2]
    have hcont :
        ContinuousOn alternatingOrdinaryHarmonicCubicKernel24
          (Icc (0 : ℝ) (1 / 2)) := by
      unfold alternatingOrdinaryHarmonicCubicKernel24
      apply ContinuousOn.div
        (continuousOn_const.sub htri)
        (continuousOn_const.sub continuousOn_id)
      intro x hx
      exact ne_of_gt (show 0 < 1 - x by linarith [hx.2])
    apply ContinuousOn.intervalIntegrable
    rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2)]
    exact hcont
  · have hint :
        IntervalIntegrable
          alternatingOrdinaryHarmonicCubicKernelExtension24
          MeasureTheory.volume (1 / 2) 1 := by
      apply ContinuousOn.intervalIntegrable
      rw [Set.uIcc_of_le (by norm_num : (1 / 2 : ℝ) ≤ 1)]
      exact
        alternatingOrdinaryHarmonicCubicKernelExtension24_continuousOn
    apply hint.congr_ae
    filter_upwards [
      MeasureTheory.ae_restrict_of_ae
        (MeasureTheory.Measure.ae_ne MeasureTheory.volume (1 : ℝ))
    ] with x hxne
    simp [alternatingOrdinaryHarmonicCubicKernelExtension24, hxne]

private theorem alternatingBarHarmonicCubicKernel24_intervalIntegrable :
    IntervalIntegrable alternatingBarHarmonicCubicKernel24
      MeasureTheory.volume 0 1 := by
  apply ContinuousOn.intervalIntegrable
  rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
  unfold alternatingBarHarmonicCubicKernel24
  apply ContinuousOn.div
  · exact
      continuousOn_const.sub
        (RamanujanChallenge.P26.trilog26_continuousOn_unit.mono
          (by
            intro x hx
            constructor <;> linarith [hx.1, hx.2]))
  · exact continuousOn_const.add continuousOn_id
  · intro x hx
    exact ne_of_gt (show 0 < 1 + x by linarith [hx.1])

private def alternatingCubicKernelDifferencePrimitive24
    (x : ℝ) : ℝ :=
  -((RamanujanChallenge.P26.trilog26 (-1) -
        RamanujanChallenge.P26.trilog26 (-x)) *
      Real.log (1 - x)) +
    dilog (-x) * dilog x -
    (RamanujanChallenge.P26.trilog26 (-1) -
        RamanujanChallenge.P26.trilog26 x) *
      Real.log (1 + x)

private theorem alternatingCubicKernelDifferencePrimitive24_hasDerivAt
    {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasDerivAt alternatingCubicKernelDifferencePrimitive24
      (alternatingOrdinaryHarmonicCubicKernel24 x -
        alternatingBarHarmonicCubicKernel24 x) x := by
  have hxne : x ≠ 0 := ne_of_gt hx0
  have h1xne : 1 - x ≠ 0 := ne_of_gt (sub_pos.mpr hx1)
  have h1pxne : 1 + x ≠ 0 := by linarith
  have hneg :
      HasDerivAt (fun y : ℝ => -y) (-1) x := by
    simpa using (hasDerivAt_id x).neg
  have htNeg :=
    (RamanujanChallenge.P26.trilog26_hasDerivAt_of_abs_lt_one
      (x := -x)
      (by rw [abs_neg, abs_of_pos hx0]; exact hx1)
      (neg_ne_zero.mpr hxne)).comp x hneg
  have ht :=
    RamanujanChallenge.P26.trilog26_hasDerivAt_of_abs_lt_one
      (by rw [abs_of_pos hx0]; exact hx1) hxne
  have hdNeg :=
    (dilog_hasDerivAt_of_abs_lt_one
      (x := -x)
      (by rw [abs_neg, abs_of_pos hx0]; exact hx1)
      (neg_ne_zero.mpr hxne)).comp x hneg
  have hd := dilog_hasDerivAt hx0 hx1
  have hlogSub :
      HasDerivAt (fun y : ℝ => Real.log (1 - y))
        (-1 / (1 - x)) x := by
    convert
      ((hasDerivAt_const x 1).sub (hasDerivAt_id x)).log h1xne
      using 1
    simp
  have hlogAdd :
      HasDerivAt (fun y : ℝ => Real.log (1 + y))
        (1 / (1 + x)) x := by
    convert
      ((hasDerivAt_const x 1).add (hasDerivAt_id x)).log h1pxne
      using 1
    simp
  have huNeg :=
    (hasDerivAt_const x
      (RamanujanChallenge.P26.trilog26 (-1))).sub htNeg
  have hu :=
    (hasDerivAt_const x
      (RamanujanChallenge.P26.trilog26 (-1))).sub ht
  unfold alternatingCubicKernelDifferencePrimitive24
  have htotal :=
    (huNeg.mul hlogSub).neg.add (hdNeg.mul hd) |>.sub
      (hu.mul hlogAdd)
  convert htotal using 1
  unfold alternatingOrdinaryHarmonicCubicKernel24
    alternatingBarHarmonicCubicKernel24
  simp only [Pi.sub_apply, Pi.neg_apply, id_eq, Function.comp_apply]
  field_simp [hxne, h1xne, h1pxne]
  ring

private theorem
    alternatingCubicKernelDifferencePrimitive24_tendsto_zero :
    Tendsto alternatingCubicKernelDifferencePrimitive24
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  have htriNeg :
      ContinuousAt
        (fun x : ℝ => RamanujanChallenge.P26.trilog26 (-x)) 0 := by
    have hneg0 :
        ContinuousAt (fun x : ℝ => -x) 0 :=
      continuousAt_id.neg
    have houter :
        ContinuousAt RamanujanChallenge.P26.trilog26 (-(0 : ℝ)) :=
      RamanujanChallenge.P26.trilog26_continuousOn_unit.continuousAt
        (Icc_mem_nhds (by norm_num : (-1 : ℝ) < -(0 : ℝ))
          (by norm_num : -(0 : ℝ) < 1))
    simpa only [Function.comp_apply] using
      houter.comp hneg0
  have htri :
      ContinuousAt RamanujanChallenge.P26.trilog26 0 :=
    RamanujanChallenge.P26.trilog26_continuousOn_unit.continuousAt
      (Icc_mem_nhds (by norm_num : (-1 : ℝ) < 0)
        (by norm_num : (0 : ℝ) < 1))
  have hdNeg :
      ContinuousAt (fun x : ℝ => dilog (-x)) 0 :=
    by
      have hneg0 :
          ContinuousAt (fun x : ℝ => -x) 0 :=
        continuousAt_id.neg
      have houter :
          ContinuousAt dilog (-(0 : ℝ)) :=
        dilog_continuousOn_unit.continuousAt
          (Icc_mem_nhds (by norm_num : (-1 : ℝ) < -(0 : ℝ))
            (by norm_num : -(0 : ℝ) < 1))
      simpa only [Function.comp_apply] using
        houter.comp hneg0
  have hd :
      ContinuousAt dilog 0 :=
    dilog_continuousOn_unit.continuousAt
      (Icc_mem_nhds (by norm_num : (-1 : ℝ) < 0)
        (by norm_num : (0 : ℝ) < 1))
  have hlogSub :
      ContinuousAt (fun x : ℝ => Real.log (1 - x)) 0 := by
    exact
      (continuousAt_const.sub continuousAt_id).log
        (by norm_num : (1 - (0 : ℝ)) ≠ 0)
  have hlogAdd :
      ContinuousAt (fun x : ℝ => Real.log (1 + x)) 0 := by
    exact
      (continuousAt_const.add continuousAt_id).log
        (by norm_num : (1 + (0 : ℝ)) ≠ 0)
  have hcont :
      ContinuousAt alternatingCubicKernelDifferencePrimitive24 0 := by
    unfold alternatingCubicKernelDifferencePrimitive24
    exact
      (((continuousAt_const.sub htriNeg).mul hlogSub).neg.add
        (hdNeg.mul hd)).sub
        ((continuousAt_const.sub htri).mul hlogAdd)
  simpa [alternatingCubicKernelDifferencePrimitive24,
    RamanujanChallenge.P26.trilog26_zero, dilog_zero] using
    tendsto_nhdsWithin_of_tendsto_nhds hcont.tendsto

private theorem
    alternatingCubicKernelDifferencePrimitive24_tendsto_one :
    Tendsto alternatingCubicKernelDifferencePrimitive24
      (𝓝[<] (1 : ℝ))
      (𝓝 (-(1 / 2 : ℝ) * (Real.pi ^ 2 / 6) ^ 2 +
        (7 / 4 : ℝ) * Real.log 2 * zeta3_24)) := by
  have hsub :
      Tendsto (fun x : ℝ => 1 - x)
        (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    have hcont :
        ContinuousAt (fun x : ℝ => 1 - x) 1 :=
      continuousAt_const.sub continuousAt_id
    simpa using hcont.tendsto.mono_left inf_le_left
  have hmulLog :
      Tendsto (fun x : ℝ => (1 - x) * Real.log (1 - x))
        (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    simpa using
      Real.continuous_mul_log.continuousAt.tendsto.comp hsub
  have hrawProduct :
      Tendsto
        (fun x : ℝ =>
          alternatingOrdinaryHarmonicCubicKernel24 x *
            ((1 - x) * Real.log (1 - x)))
        (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    simpa using
      alternatingOrdinaryHarmonicCubicKernel24_tendsto_one.mul
        hmulLog
  have hfirst :
      Tendsto
        (fun x : ℝ =>
          -((RamanujanChallenge.P26.trilog26 (-1) -
              RamanujanChallenge.P26.trilog26 (-x)) *
            Real.log (1 - x)))
        (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    have hproduct :
        Tendsto
          (fun x : ℝ =>
            (RamanujanChallenge.P26.trilog26 (-1) -
                RamanujanChallenge.P26.trilog26 (-x)) *
              Real.log (1 - x))
          (𝓝[<] (1 : ℝ)) (𝓝 0) := by
      apply hrawProduct.congr'
      filter_upwards [self_mem_nhdsWithin] with x hx
      have hxne : 1 - x ≠ 0 := ne_of_gt (sub_pos.mpr hx)
      unfold alternatingOrdinaryHarmonicCubicKernel24
      field_simp [hxne]
    simpa using hproduct.neg
  have hdNeg :
      Tendsto (fun x : ℝ => dilog (-x))
        (𝓝[<] (1 : ℝ)) (𝓝 (-(Real.pi ^ 2 / 12))) := by
    have houter :
        ContinuousWithinAt dilog (Ici (-1 : ℝ)) (-1) :=
      (dilog_continuousOn_unit (-1)
        ⟨le_rfl, by norm_num⟩).mono_of_mem_nhdsWithin
        (Icc_mem_nhdsGE (by norm_num : (-1 : ℝ) < 1))
    have hinner :
        ContinuousWithinAt (fun x : ℝ => -x)
          (Iic (1 : ℝ)) 1 :=
      continuousAt_id.neg.continuousWithinAt
    have hcomp :=
      houter.comp hinner (by
        intro x hx
        change x ≤ 1 at hx
        change -1 ≤ -x
        linarith)
    have ht :
        Tendsto (fun x : ℝ => dilog (-x))
          (𝓝[<] (1 : ℝ)) (𝓝 (dilog (-1))) := by
      simpa only [Function.comp_apply] using
        (hcomp.mono Iio_subset_Iic_self).tendsto
    convert ht using 1
    rw [RamanujanChallenge.P26.dilog_neg_one26]
    ring
  have hd :
      Tendsto dilog (𝓝[<] (1 : ℝ))
        (𝓝 (Real.pi ^ 2 / 6)) := by
    have hWithin :
        ContinuousWithinAt dilog (Iio (1 : ℝ)) 1 :=
      (dilog_continuousOn_unit 1 (by norm_num)).mono_of_mem_nhdsWithin
        (Icc_mem_nhdsLT (show (-1 : ℝ) < 1 by norm_num))
    simpa [dilog_one] using hWithin.tendsto
  have hthird :
      Tendsto
        (fun x : ℝ =>
          -((RamanujanChallenge.P26.trilog26 (-1) -
              RamanujanChallenge.P26.trilog26 x) *
            Real.log (1 + x)))
        (𝓝[<] (1 : ℝ))
        (𝓝 ((7 / 4 : ℝ) * Real.log 2 * zeta3_24)) := by
    have htri :
        Tendsto RamanujanChallenge.P26.trilog26
          (𝓝[<] (1 : ℝ)) (𝓝 zeta3_24) := by
      have hWithin :
          ContinuousWithinAt RamanujanChallenge.P26.trilog26
            (Iio (1 : ℝ)) 1 :=
        (RamanujanChallenge.P26.trilog26_continuousOn_unit
          1 (by norm_num)).mono_of_mem_nhdsWithin
          (Icc_mem_nhdsLT (show (-1 : ℝ) < 1 by norm_num))
      simpa [RamanujanChallenge.P26.trilog26_one,
        RamanujanChallenge.P26.zeta3, zeta3_24] using
        hWithin.tendsto
    have hlog :
        Tendsto (fun x : ℝ => Real.log (1 + x))
          (𝓝[<] (1 : ℝ)) (𝓝 (Real.log 2)) := by
      have hcont :
          ContinuousAt (fun x : ℝ => Real.log (1 + x)) 1 :=
        (continuousAt_const.add continuousAt_id).log
          (by norm_num : (1 + (1 : ℝ)) ≠ 0)
      have ht :
          Tendsto (fun x : ℝ => Real.log (1 + x))
            (𝓝[<] (1 : ℝ)) (𝓝 (Real.log (1 + 1))) :=
        tendsto_nhdsWithin_of_tendsto_nhds hcont.tendsto
      convert ht using 1 <;> norm_num
    have hconst :
        RamanujanChallenge.P26.trilog26 (-1) =
          -(3 / 4 : ℝ) * zeta3_24 := by
      rw [RamanujanChallenge.P26.trilog26_neg_one]
      unfold RamanujanChallenge.P26.zeta3 zeta3_24
      ring_nf
    rw [hconst]
    convert ((tendsto_const_nhds.sub htri).mul hlog).neg using 1
    ring
  unfold alternatingCubicKernelDifferencePrimitive24
  have htotal := (hfirst.add (hdNeg.mul hd)).add hthird
  convert htotal using 1 <;> ring

private theorem alternatingCubicKernelDifferenceIntegral24 :
    (∫ x : ℝ in 0..1,
      (alternatingOrdinaryHarmonicCubicKernel24 x -
        alternatingBarHarmonicCubicKernel24 x)) =
      -(1 / 2 : ℝ) * (Real.pi ^ 2 / 6) ^ 2 +
        (7 / 4 : ℝ) * Real.log 2 * zeta3_24 := by
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_tendsto
    (f := alternatingCubicKernelDifferencePrimitive24)
    (fa := (0 : ℝ))
    (fb := -(1 / 2 : ℝ) * (Real.pi ^ 2 / 6) ^ 2 +
      (7 / 4 : ℝ) * Real.log 2 * zeta3_24)
    (by norm_num)
    (fun x hx =>
      alternatingCubicKernelDifferencePrimitive24_hasDerivAt
        hx.1 hx.2)
    (alternatingOrdinaryHarmonicCubicKernel24_intervalIntegrable.sub
      alternatingBarHarmonicCubicKernel24_intervalIntegrable)
    alternatingCubicKernelDifferencePrimitive24_tendsto_zero
    alternatingCubicKernelDifferencePrimitive24_tendsto_one]
  ring

@[simp] private theorem polylog4_zero24 :
    polylog4 0 = 0 := by
  simp [polylog4, zero_pow (Nat.succ_ne_zero _)]

private theorem polylog4_one24 :
    polylog4 1 = Real.pi ^ 4 / 90 := by
  unfold polylog4
  simpa only [one_pow, Nat.cast_add, Nat.cast_one] using
    shifted_zeta_four_hasSum24.tsum_eq

private theorem polylog4_continuousOn_unit24 :
    ContinuousOn polylog4 (Icc (-1 : ℝ) 1) := by
  unfold polylog4
  have hu :
      Summable
        (fun n : ℕ =>
          (1 : ℝ) / (((n + 1 : ℕ) : ℝ) ^ 4)) := by
    convert shifted_zeta_four_hasSum24.summable using 1
    funext n
    simp only [Nat.cast_add, Nat.cast_one]
  refine continuousOn_tsum
    (u := fun n : ℕ =>
      (1 : ℝ) / (((n + 1 : ℕ) : ℝ) ^ 4))
    (fun n => ?_) hu ?_
  · fun_prop
  · intro n x hx
    rw [Real.norm_eq_abs, abs_div, abs_pow]
    have habs : |x| ≤ 1 := (abs_le).2 hx
    have hpos :
        (0 : ℝ) ≤ (((n + 1 : ℕ) : ℝ) ^ 4) := by
      positivity
    rw [abs_of_nonneg hpos]
    apply div_le_div_of_nonneg_right _ hpos
    calc
      |x| ^ (n + 1) ≤ 1 ^ (n + 1) :=
        pow_le_pow_left₀ (abs_nonneg x) habs (n + 1)
      _ = 1 := one_pow _

private theorem polylog4_hasDerivAt24
    {x : ℝ} (hx : |x| < 1) (hxne : x ≠ 0) :
    HasDerivAt polylog4
      (RamanujanChallenge.P26.trilog26 x / x) x := by
  unfold polylog4
  let r : ℝ := (|x| + 1) / 2
  have hr0 : 0 ≤ r := by dsimp [r]; positivity
  have hrpos : 0 < r := by dsimp [r]; positivity
  have hr1 : r < 1 := by dsimp [r]; linarith
  have hxr : |x| < r := by dsimp [r]; linarith
  have hu : Summable (fun n : ℕ => r ^ n) :=
    summable_geometric_of_lt_one hr0 hr1
  have hterm :
      ∀ n : ℕ, ∀ y : ℝ, y ∈ Ioo (-r) r →
        HasDerivAt
          (fun y : ℝ =>
            y ^ (n + 1) / (((n + 1 : ℕ) : ℝ) ^ 4))
          (y ^ n / (((n + 1 : ℕ) : ℝ) ^ 3)) y := by
    intro n y _
    convert
      (hasDerivAt_pow (n + 1) y).div_const
        ((((n + 1 : ℕ) : ℝ) ^ 4)) using 1
    rw [Nat.add_sub_cancel]
    push_cast
    field_simp
  have hbound :
      ∀ n : ℕ, ∀ y : ℝ, y ∈ Ioo (-r) r →
        ‖y ^ n / (((n + 1 : ℕ) : ℝ) ^ 3)‖ ≤ r ^ n := by
    intro n y hy
    rw [Real.norm_eq_abs, abs_div, abs_pow]
    have hyr : |y| < r := (abs_lt).2 hy
    have hden :
        (1 : ℝ) ≤ |(((n + 1 : ℕ) : ℝ) ^ 3)| := by
      rw [abs_of_nonneg
        (by positivity :
          (0 : ℝ) ≤ (((n + 1 : ℕ) : ℝ) ^ 3))]
      exact one_le_pow₀ (by norm_num)
    calc
      |y| ^ n / |(((n + 1 : ℕ) : ℝ) ^ 3)| ≤
          |y| ^ n / 1 := by
        gcongr
      _ = |y| ^ n := by ring
      _ ≤ r ^ n := pow_le_pow_left₀ (abs_nonneg y) hyr.le n
  have hzero :
      Summable
        (fun n : ℕ =>
          (0 : ℝ) ^ (n + 1) /
            (((n + 1 : ℕ) : ℝ) ^ 4)) := by
    simp [zero_pow (Nat.succ_ne_zero _)]
  have hxmem : x ∈ Ioo (-r) r := (abs_lt).1 hxr
  have hd :
      HasDerivAt
        (fun y : ℝ =>
          ∑' n : ℕ,
            y ^ (n + 1) / (((n + 1 : ℕ) : ℝ) ^ 4))
        (∑' n : ℕ,
          x ^ n / (((n + 1 : ℕ) : ℝ) ^ 3)) x := by
    exact
      hasDerivAt_tsum_of_isPreconnected hu isOpen_Ioo
        (convex_Ioo (-r) r).isPreconnected hterm hbound
        (show (0 : ℝ) ∈ Ioo (-r) r by constructor <;> linarith)
        hzero hxmem
  convert hd using 1
  have hsum :=
    (trilog26_hasSum24 (x := x) (le_of_lt hx)).div_const x
  have hsum' :
      HasSum
        (fun n : ℕ =>
          x ^ n / (((n + 1 : ℕ) : ℝ) ^ 3))
        (RamanujanChallenge.P26.trilog26 x / x) := by
    convert hsum using 1
    funext n
    simp only [Nat.cast_add, Nat.cast_one]
    rw [pow_succ]
    field_simp [hxne]
    ring
  exact hsum'.tsum_eq.symm

private def logCubeOneSubPrimitive24 (x : ℝ) : ℝ :=
  -Real.log x ^ 3 * Real.log (1 - x) -
    3 * Real.log x ^ 2 * dilog x +
    6 * Real.log x * RamanujanChallenge.P26.trilog26 x -
    6 * polylog4 x

private theorem logCubeOneSubPrimitive24_hasDerivAt
    {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasDerivAt logCubeOneSubPrimitive24
      (Real.log x ^ 3 / (1 - x)) x := by
  have hxne : x ≠ 0 := ne_of_gt hx0
  have h1xne : 1 - x ≠ 0 := ne_of_gt (sub_pos.mpr hx1)
  have hlog :
      HasDerivAt Real.log (1 / x) x := by
    simpa [one_div] using Real.hasDerivAt_log hxne
  have hlogSub :
      HasDerivAt (fun y : ℝ => Real.log (1 - y))
        (-1 / (1 - x)) x := by
    convert
      ((hasDerivAt_const x 1).sub (hasDerivAt_id x)).log h1xne
      using 1
    simp
  have hd := dilog_hasDerivAt hx0 hx1
  have ht :=
    RamanujanChallenge.P26.trilog26_hasDerivAt_of_abs_lt_one
      (by rw [abs_of_pos hx0]; exact hx1) hxne
  have hq :=
    polylog4_hasDerivAt24
      (by rw [abs_of_pos hx0]; exact hx1) hxne
  unfold logCubeOneSubPrimitive24
  have htotal :=
    ((hlog.pow 3).mul hlogSub).neg.sub
      (((hlog.pow 2).mul hd).const_mul 3) |>.add
      ((hlog.mul ht).const_mul 6) |>.sub
      (hq.const_mul 6)
  convert htotal using 1
  · funext y
    simp only [Pi.add_apply, Pi.sub_apply, Pi.neg_apply,
      Pi.mul_apply, Pi.pow_apply]
    ring
  · norm_num
    field_simp [hxne, h1xne]
    ring

private def halfLogCubeOneSubKernel24 (x : ℝ) : ℝ :=
  Real.log (1 - x) ^ 3 / x

private theorem halfLogCubeOneSubKernel24_continuousOn :
    ContinuousOn halfLogCubeOneSubKernel24
      (Icc (0 : ℝ) (1 / 2)) := by
  have haux :
      ContinuousOn
        (fun x : ℝ => x ^ 2 * logOneMinusSlope24 x ^ 3)
        (Icc (0 : ℝ) (1 / 2)) :=
    (continuousOn_id.pow 2).mul
      (logOneMinusSlope24_continuousOn.pow 3)
  apply haux.congr
  intro x _
  by_cases hxzero : x = 0
  · subst x
    simp [halfLogCubeOneSubKernel24]
  · simp [halfLogCubeOneSubKernel24, logOneMinusSlope24, hxzero]
    field_simp [hxzero]

private theorem halfLogCubeOneSubKernel24_intervalIntegrable :
    IntervalIntegrable halfLogCubeOneSubKernel24
      MeasureTheory.volume 0 (1 / 2) := by
  apply ContinuousOn.intervalIntegrable
  rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2)]
  exact halfLogCubeOneSubKernel24_continuousOn

private def halfLogCubeOneSubPrimitive24 (x : ℝ) : ℝ :=
  -logCubeOneSubPrimitive24 (1 - x)

private theorem halfLogCubeOneSubPrimitive24_hasDerivAt
    {x : ℝ} (hx0 : 0 < x) (hxhalf : x < 1 / 2) :
    HasDerivAt halfLogCubeOneSubPrimitive24
      (halfLogCubeOneSubKernel24 x) x := by
  have hy0 : 0 < 1 - x := by linarith
  have hy1 : 1 - x < 1 := by linarith
  have hsub :
      HasDerivAt (fun y : ℝ => 1 - y) (-1) x := by
    convert (hasDerivAt_const x 1).sub (hasDerivAt_id x) using 1
    simp
  have hp :=
    (logCubeOneSubPrimitive24_hasDerivAt hy0 hy1).comp x hsub
  unfold halfLogCubeOneSubPrimitive24
  convert hp.neg using 1
  unfold halfLogCubeOneSubKernel24
  have hxne : x ≠ 0 := ne_of_gt hx0
  field_simp [hxne]
  ring

private theorem halfLogCubeOneSubPrimitive24_half :
    halfLogCubeOneSubPrimitive24 (1 / 2) =
      6 * polylog4 (1 / 2) +
        (1 / 2 : ℝ) * Real.log 2 ^ 4 -
        (3 / 2 : ℝ) * Real.log 2 ^ 2 *
          (Real.pi ^ 2 / 6) +
        (21 / 4 : ℝ) * Real.log 2 * zeta3_24 := by
  have hloghalf : Real.log (1 / 2 : ℝ) = -Real.log 2 := by
    rw [one_div, Real.log_inv]
  unfold halfLogCubeOneSubPrimitive24
  rw [show 1 - (1 / 2 : ℝ) = 1 / 2 by norm_num]
  unfold logCubeOneSubPrimitive24
  rw [hloghalf, RamanujanChallenge.P26.dilog26_half,
    RamanujanChallenge.P26.trilog26_half]
  rw [show 1 - (1 / 2 : ℝ) = 1 / 2 by norm_num, hloghalf]
  unfold RamanujanChallenge.P26.zeta3 zeta3_24
  ring

private theorem logCubeOneSubPrimitive24_tendsto_one :
    Tendsto logCubeOneSubPrimitive24
      (𝓝[<] (1 : ℝ)) (𝓝 (-(6 : ℝ) * (Real.pi ^ 4 / 90))) := by
  have hlog :
      Tendsto Real.log (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    have hcont : ContinuousAt Real.log 1 :=
      Real.continuousAt_log (by norm_num)
    simpa using
      (tendsto_nhdsWithin_of_tendsto_nhds hcont.tendsto :
        Tendsto Real.log (𝓝[<] (1 : ℝ)) (𝓝 (Real.log 1)))
  have hslopeWithin :
      HasDerivWithinAt Real.log 1 (Iio (1 : ℝ)) 1 :=
    by
      simpa using
        (Real.hasDerivAt_log
          (by norm_num : (1 : ℝ) ≠ 0)).hasDerivWithinAt
  have hslope :=
    (hasDerivWithinAt_iff_tendsto_slope'
      (show (1 : ℝ) ∉ Iio 1 by simp)).1 hslopeWithin
  have hratio :
      Tendsto (fun x : ℝ => Real.log x / (1 - x))
        (𝓝[<] (1 : ℝ)) (𝓝 (-1)) := by
    apply hslope.neg.congr'
    filter_upwards [self_mem_nhdsWithin] with x hx
    have hxne : x ≠ 1 := ne_of_lt hx
    simp only [slope, vsub_eq_sub, smul_eq_mul, Real.log_one,
      sub_zero]
    have hxsub : x - 1 ≠ 0 := sub_ne_zero.mpr hxne
    have hsubx : 1 - x ≠ 0 := sub_ne_zero.mpr hxne.symm
    field_simp [hxsub, hsubx]
    ring
  have hsub :
      Tendsto (fun x : ℝ => 1 - x)
        (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    have hcont :
        ContinuousAt (fun x : ℝ => 1 - x) 1 :=
      continuousAt_const.sub continuousAt_id
    simpa using
      (tendsto_nhdsWithin_of_tendsto_nhds hcont.tendsto :
        Tendsto (fun x : ℝ => 1 - x)
          (𝓝[<] (1 : ℝ)) (𝓝 (1 - 1)))
  have hsubMulLog :
      Tendsto (fun x : ℝ => (1 - x) * Real.log (1 - x))
        (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    simpa using
      Real.continuous_mul_log.continuousAt.tendsto.comp hsub
  have hfirst :
      Tendsto
        (fun x : ℝ => -Real.log x ^ 3 * Real.log (1 - x))
        (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    have hraw :=
      (hratio.pow 3).mul ((hsub.pow 2).mul hsubMulLog)
    have hraw0 :
        Tendsto
          (fun x : ℝ =>
            -((Real.log x / (1 - x)) ^ 3 *
              ((1 - x) ^ 2 *
                ((1 - x) * Real.log (1 - x)))))
          (𝓝[<] (1 : ℝ)) (𝓝 0) := by
      simpa using hraw.neg
    apply hraw0.congr'
    filter_upwards [self_mem_nhdsWithin] with x hx
    have hxne : 1 - x ≠ 0 := ne_of_gt (sub_pos.mpr hx)
    field_simp [hxne]
  have hd :
      Tendsto dilog (𝓝[<] (1 : ℝ))
        (𝓝 (Real.pi ^ 2 / 6)) := by
    have hWithin :
        ContinuousWithinAt dilog (Iio (1 : ℝ)) 1 :=
      (dilog_continuousOn_unit 1 (by norm_num)).mono_of_mem_nhdsWithin
        (Icc_mem_nhdsLT (show (-1 : ℝ) < 1 by norm_num))
    simpa [dilog_one] using hWithin.tendsto
  have ht :
      Tendsto RamanujanChallenge.P26.trilog26
        (𝓝[<] (1 : ℝ)) (𝓝 zeta3_24) := by
    have hWithin :
        ContinuousWithinAt RamanujanChallenge.P26.trilog26
          (Iio (1 : ℝ)) 1 :=
      (RamanujanChallenge.P26.trilog26_continuousOn_unit
        1 (by norm_num)).mono_of_mem_nhdsWithin
        (Icc_mem_nhdsLT (show (-1 : ℝ) < 1 by norm_num))
    simpa [RamanujanChallenge.P26.trilog26_one,
      RamanujanChallenge.P26.zeta3, zeta3_24] using hWithin.tendsto
  have hq :
      Tendsto polylog4 (𝓝[<] (1 : ℝ))
        (𝓝 (Real.pi ^ 4 / 90)) := by
    have hWithin :
        ContinuousWithinAt polylog4 (Iio (1 : ℝ)) 1 :=
      (polylog4_continuousOn_unit24
        1 (by norm_num)).mono_of_mem_nhdsWithin
        (Icc_mem_nhdsLT (show (-1 : ℝ) < 1 by norm_num))
    simpa [polylog4_one24] using hWithin.tendsto
  unfold logCubeOneSubPrimitive24
  have htotal :=
    ((hfirst.sub ((hlog.pow 2).mul hd |>.const_mul 3)).add
      (hlog.mul ht |>.const_mul 6)).sub
      (hq.const_mul 6)
  convert htotal using 1 <;> ring

private theorem halfLogCubeOneSubPrimitive24_tendsto_zero :
    Tendsto halfLogCubeOneSubPrimitive24
      (𝓝[>] (0 : ℝ)) (𝓝 (6 * (Real.pi ^ 4 / 90))) := by
  unfold halfLogCubeOneSubPrimitive24
  have honeSub :
      Tendsto (fun x : ℝ => 1 - x)
        (𝓝[>] (0 : ℝ)) (𝓝[<] (1 : ℝ)) := by
    rw [tendsto_nhdsWithin_iff]
    constructor
    · have hc :
          ContinuousAt (fun x : ℝ => 1 - x) 0 :=
        continuousAt_const.sub continuousAt_id
      simpa using
        hc.tendsto.mono_left
          (show (𝓝[>] (0 : ℝ)) ≤ 𝓝 0 from inf_le_left)
    · filter_upwards [self_mem_nhdsWithin] with x hx
      change 0 < x at hx
      change 1 - x < 1
      linarith
  have hcomp :=
    logCubeOneSubPrimitive24_tendsto_one.comp
      honeSub
  convert hcomp.neg using 1 <;> ring

private theorem halfLogCubeOneSubPrimitive24_tendsto_half :
    Tendsto halfLogCubeOneSubPrimitive24
      (𝓝[<] (1 / 2 : ℝ))
      (𝓝 (6 * polylog4 (1 / 2) +
        (1 / 2 : ℝ) * Real.log 2 ^ 4 -
        (3 / 2 : ℝ) * Real.log 2 ^ 2 *
          (Real.pi ^ 2 / 6) +
        (21 / 4 : ℝ) * Real.log 2 * zeta3_24)) := by
  have hy : 1 - (1 / 2 : ℝ) = 1 / 2 := by norm_num
  have hcontP :
      ContinuousAt logCubeOneSubPrimitive24 (1 / 2) := by
    unfold logCubeOneSubPrimitive24
    have hlog : ContinuousAt Real.log (1 / 2) :=
      Real.continuousAt_log (by norm_num)
    have hlogSub :
        ContinuousAt (fun x : ℝ => Real.log (1 - x)) (1 / 2) := by
      exact
        (continuousAt_const.sub continuousAt_id).log
          (by norm_num : 1 - (1 / 2 : ℝ) ≠ 0)
    have hd : ContinuousAt dilog (1 / 2) :=
      dilog_continuousOn_unit.continuousAt
        (Icc_mem_nhds (by norm_num : (-1 : ℝ) < 1 / 2)
          (by norm_num : (1 / 2 : ℝ) < 1))
    have ht :
        ContinuousAt RamanujanChallenge.P26.trilog26 (1 / 2) :=
      RamanujanChallenge.P26.trilog26_continuousOn_unit.continuousAt
        (Icc_mem_nhds (by norm_num : (-1 : ℝ) < 1 / 2)
          (by norm_num : (1 / 2 : ℝ) < 1))
    have hq : ContinuousAt polylog4 (1 / 2) :=
      polylog4_continuousOn_unit24.continuousAt
        (Icc_mem_nhds (by norm_num : (-1 : ℝ) < 1 / 2)
          (by norm_num : (1 / 2 : ℝ) < 1))
    have hc :=
      (((hlog.pow 3).mul hlogSub).neg.sub
        (((hlog.pow 2).mul hd).const_mul 3) |>.add
        ((hlog.mul ht).const_mul 6)).sub
        (hq.const_mul 6)
    convert hc using 1
    funext y
    simp only [Pi.add_apply, Pi.sub_apply, Pi.neg_apply,
      Pi.mul_apply, Pi.pow_apply]
    ring
  have hsubCont :
      ContinuousAt (fun x : ℝ => 1 - x) (1 / 2) :=
    continuousAt_const.sub continuousAt_id
  have houter :
      ContinuousAt logCubeOneSubPrimitive24
        (1 - (1 / 2 : ℝ)) := by
    rw [hy]
    exact hcontP
  have hcomp := houter.comp hsubCont
  have hprim :
      ContinuousAt halfLogCubeOneSubPrimitive24 (1 / 2) := by
    unfold halfLogCubeOneSubPrimitive24
    simpa only [Function.comp_apply] using hcomp.neg
  have ht :
      Tendsto halfLogCubeOneSubPrimitive24
        (𝓝[<] (1 / 2 : ℝ))
        (𝓝 (halfLogCubeOneSubPrimitive24 (1 / 2))) :=
    tendsto_nhdsWithin_of_tendsto_nhds hprim.tendsto
  rw [halfLogCubeOneSubPrimitive24_half] at ht
  exact ht

private theorem halfLogCubeOneSubIntegral24 :
    (∫ x : ℝ in 0..(1 / 2),
      halfLogCubeOneSubKernel24 x) =
      6 * polylog4 (1 / 2) +
        (1 / 2 : ℝ) * Real.log 2 ^ 4 -
        (3 / 2 : ℝ) * Real.log 2 ^ 2 *
          (Real.pi ^ 2 / 6) +
        (21 / 4 : ℝ) * Real.log 2 * zeta3_24 -
        6 * (Real.pi ^ 4 / 90) := by
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_tendsto
    (f := halfLogCubeOneSubPrimitive24)
    (fa := 6 * (Real.pi ^ 4 / 90))
    (fb := 6 * polylog4 (1 / 2) +
      (1 / 2 : ℝ) * Real.log 2 ^ 4 -
      (3 / 2 : ℝ) * Real.log 2 ^ 2 *
        (Real.pi ^ 2 / 6) +
      (21 / 4 : ℝ) * Real.log 2 * zeta3_24)
    (by norm_num)
    (fun x hx =>
      halfLogCubeOneSubPrimitive24_hasDerivAt hx.1 hx.2)
    halfLogCubeOneSubKernel24_intervalIntegrable
    halfLogCubeOneSubPrimitive24_tendsto_zero
    halfLogCubeOneSubPrimitive24_tendsto_half]

private theorem harmonicConvolutionCoefficient24 (n : ℕ) :
    (∑ k ∈ Finset.range (n + 1),
      (1 / (k + 1 : ℝ)) *
        (1 / (((n - k + 1 : ℕ) : ℝ)))) =
      2 * harmonicNumber (n + 1) / (n + 2 : ℝ) := by
  have hpoint :
      ∀ k ∈ Finset.range (n + 1),
        (1 / (k + 1 : ℝ)) *
            (1 / (((n - k + 1 : ℕ) : ℝ))) =
          (1 / (n + 2 : ℝ)) *
            (1 / (k + 1 : ℝ) +
              1 / (((n - k + 1 : ℕ) : ℝ))) := by
    intro k hk
    have hkn : k ≤ n := Nat.le_of_lt_succ (Finset.mem_range.mp hk)
    have hkpos : (k + 1 : ℝ) ≠ 0 := by positivity
    have hnkpos : (((n - k + 1 : ℕ) : ℝ)) ≠ 0 := by
      positivity
    have hnpos : (n + 2 : ℝ) ≠ 0 := by positivity
    have hcast :
        (k + 1 : ℝ) + (((n - k + 1 : ℕ) : ℝ)) =
          n + 2 := by
      exact_mod_cast
        (by omega :
          (k + 1) + (n - k + 1) = n + 2)
    field_simp [hkpos, hnkpos, hnpos]
    linarith
  calc
    (∑ k ∈ Finset.range (n + 1),
      (1 / (k + 1 : ℝ)) *
        (1 / (((n - k + 1 : ℕ) : ℝ)))) =
        ∑ k ∈ Finset.range (n + 1),
          (1 / (n + 2 : ℝ)) *
            (1 / (k + 1 : ℝ) +
              1 / (((n - k + 1 : ℕ) : ℝ))) := by
      apply Finset.sum_congr rfl
      exact hpoint
    _ = (1 / (n + 2 : ℝ)) *
        ((∑ k ∈ Finset.range (n + 1), 1 / (k + 1 : ℝ)) +
          ∑ k ∈ Finset.range (n + 1),
            1 / (((n - k + 1 : ℕ) : ℝ))) := by
      rw [← Finset.mul_sum, Finset.sum_add_distrib]
    _ = (1 / (n + 2 : ℝ)) *
        (harmonicNumber (n + 1) + harmonicNumber (n + 1)) := by
      unfold harmonicNumber
      have hreflect :=
        Finset.sum_range_reflect
          (fun j : ℕ => 1 / (j + 1 : ℝ)) (n + 1)
      simpa only [Nat.add_sub_cancel, Nat.cast_add, Nat.cast_one] using
        congrArg
          (fun s : ℝ =>
            (1 / (n + 2 : ℝ)) *
              ((∑ k ∈ Finset.range (n + 1), 1 / (k + 1 : ℝ)) + s))
          hreflect
    _ = 2 * harmonicNumber (n + 1) / (n + 2 : ℝ) := by
      ring

private theorem logOneSubSquare_hasSum24
    {x : ℝ} (hx0 : 0 ≤ x) (hx1 : x < 1) :
    HasSum
      (fun n : ℕ =>
        (2 * harmonicNumber (n + 1) / (n + 2 : ℝ)) *
          x ^ (n + 2))
      (Real.log (1 - x) ^ 2) := by
  have hxabs : |x| < 1 := by
    rw [abs_of_nonneg hx0]
    exact hx1
  have hlog := Real.hasSum_pow_div_log_of_abs_lt_one hxabs
  have hnorm :
      Summable
        (fun n : ℕ =>
          ‖x ^ (n + 1) / (n + 1 : ℝ)‖) := by
    have hs := hlog.summable
    apply hs.congr
    intro n
    rw [Real.norm_eq_abs, abs_div, abs_pow,
      abs_of_nonneg hx0,
      abs_of_pos (by positivity : (0 : ℝ) < n + 1)]
  have hprod :=
    hasSum_sum_range_mul_of_summable_norm hnorm hnorm
  rw [hlog.tsum_eq] at hprod
  convert hprod using 1
  · funext n
    calc
      (2 * harmonicNumber (n + 1) / (n + 2 : ℝ)) *
          x ^ (n + 2) =
          x ^ (n + 2) *
            (∑ k ∈ Finset.range (n + 1),
              (1 / (k + 1 : ℝ)) *
                (1 / (((n - k + 1 : ℕ) : ℝ)))) := by
        rw [harmonicConvolutionCoefficient24]
        ring
      _ = ∑ k ∈ Finset.range (n + 1),
          (x ^ (k + 1) / (k + 1 : ℝ)) *
            (x ^ (n - k + 1) /
              (((n - k + 1 : ℕ) : ℝ))) := by
        rw [Finset.mul_sum]
        apply Finset.sum_congr rfl
        intro k hk
        have hkn : k ≤ n :=
          Nat.le_of_lt_succ (Finset.mem_range.mp hk)
        have hexp :
            (k + 1) + (n - k + 1) = n + 2 := by omega
        symm
        calc
          (x ^ (k + 1) / (k + 1 : ℝ)) *
                (x ^ (n - k + 1) /
                  (((n - k + 1 : ℕ) : ℝ))) =
              (x ^ (k + 1) * x ^ (n - k + 1)) *
                ((1 / (k + 1 : ℝ)) *
                  (1 / (((n - k + 1 : ℕ) : ℝ)))) := by
            ring
          _ = x ^ (n + 2) *
                ((1 / (k + 1 : ℝ)) *
                  (1 / (((n - k + 1 : ℕ) : ℝ)))) := by
            rw [← pow_add, hexp]
      _ = ∑ k ∈ Finset.range (n + 1),
          (x ^ (k + 1) / (k + 1 : ℝ)) *
            (x ^ (n - k + 1) /
              (((n - k : ℕ) : ℝ) + 1)) := by
        simp only [Nat.cast_add, Nat.cast_one]
  · ring

private def quarticCoreMoment24 (n : ℕ) (x : ℝ) : ℝ :=
  (2 * harmonicNumber (n + 1) / (n + 2 : ℝ)) *
    x ^ (n + 1) * Real.log x

private def quarticCoreSeriesTerm24 (n : ℕ) : ℝ :=
  -(2 * harmonicNumber (n + 1) / (n + 2 : ℝ) ^ 3)

private theorem quarticCoreMoment24_intervalIntegrable (n : ℕ) :
    IntervalIntegrable (quarticCoreMoment24 n)
      MeasureTheory.volume 0 1 := by
  have hbase :
      IntervalIntegrable
        (fun x : ℝ => x ^ (n + 1) * Real.log x)
        MeasureTheory.volume 0 1 :=
    intervalIntegral.intervalIntegrable_log'.continuousOn_mul
      (continuousOn_pow (n + 1))
  unfold quarticCoreMoment24
  convert
    hbase.const_mul
      (2 * harmonicNumber (n + 1) / (n + 2 : ℝ)) using 1
  funext x
  ring

private theorem quarticCoreMoment24_integral (n : ℕ) :
    (∫ x : ℝ in 0..1, quarticCoreMoment24 n x) =
      quarticCoreSeriesTerm24 n := by
  unfold quarticCoreMoment24 quarticCoreSeriesTerm24
  rw [show
      (fun x : ℝ =>
        (2 * harmonicNumber (n + 1) / (n + 2 : ℝ)) *
          x ^ (n + 1) * Real.log x) =
      (fun x : ℝ =>
        (2 * harmonicNumber (n + 1) / (n + 2 : ℝ)) *
          (x ^ (n + 1) * Real.log x)) by
    funext x
    ring,
    intervalIntegral.integral_const_mul,
    RamanujanChallenge.P26.integral_pow_mul_log26]
  simp only [Nat.cast_add, Nat.cast_one]
  have hden : (n + 2 : ℝ) ≠ 0 := by positivity
  field_simp [hden]
  ring

private def shiftedHarmonicCubicTerm24 (n : ℕ) : ℝ :=
  harmonicNumber n / (n + 1 : ℝ) ^ 3

private theorem shiftedHarmonicCubicTerm24_hasSum :
    HasSum shiftedHarmonicCubicTerm24
      ((1 / 2 : ℝ) * (Real.pi ^ 2 / 6) ^ 2 -
        Real.pi ^ 4 / 90) := by
  have hOrd :
      HasSum ordinaryHarmonicCubicTerm24
        ((1 / 2 : ℝ) * (Real.pi ^ 2 / 6) ^ 2) := by
    rw [← ordinaryHarmonicCubicIntegral24]
    exact ordinaryHarmonicCubicTerm24_hasSum_integral
  have hdiff := hOrd.sub shifted_zeta_four_hasSum24
  convert hdiff using 1
  funext n
  unfold shiftedHarmonicCubicTerm24
    ordinaryHarmonicCubicTerm24
  rw [harmonicNumber_succ]
  have hden : (n + 1 : ℝ) ≠ 0 := by positivity
  field_simp [hden]
  ring

private theorem quarticCoreSeriesTerm24_hasSum :
    HasSum quarticCoreSeriesTerm24
      (-2 * ((1 / 2 : ℝ) * (Real.pi ^ 2 / 6) ^ 2 -
        Real.pi ^ 4 / 90)) := by
  have htail :
      HasSum
        (fun n : ℕ => shiftedHarmonicCubicTerm24 (n + 1))
        ((1 / 2 : ℝ) * (Real.pi ^ 2 / 6) ^ 2 -
          Real.pi ^ 4 / 90) := by
    simpa [shiftedHarmonicCubicTerm24, harmonicNumber] using
      (hasSum_nat_add_iff' 1).2 shiftedHarmonicCubicTerm24_hasSum
  convert htail.mul_left (-2) using 1
  funext n
  unfold quarticCoreSeriesTerm24 shiftedHarmonicCubicTerm24
  simp only [Nat.cast_add, Nat.cast_one]
  have hden : (n + 2 : ℝ) ≠ 0 := by positivity
  field_simp [hden]
  ring

private theorem quarticCoreMoment24_integral_norm (n : ℕ) :
    (∫ x : ℝ in 0..1, ‖quarticCoreMoment24 n x‖) =
      ‖quarticCoreSeriesTerm24 n‖ := by
  have hnonpos :
      ∀ x ∈ Set.uIcc (0 : ℝ) 1,
        quarticCoreMoment24 n x ≤ 0 := by
    intro x hx
    rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] at hx
    have hlog : Real.log x ≤ 0 :=
      Real.log_nonpos hx.1 hx.2
    unfold quarticCoreMoment24
    have hcoeff :
        0 ≤ 2 * harmonicNumber (n + 1) / (n + 2 : ℝ) := by
      exact div_nonneg
        (mul_nonneg (by norm_num) (harmonicNumber_nonneg (n + 1)))
        (by positivity)
    exact mul_nonpos_of_nonneg_of_nonpos
      (mul_nonneg hcoeff (pow_nonneg hx.1 _)) hlog
  calc
    (∫ x : ℝ in 0..1, ‖quarticCoreMoment24 n x‖) =
        -(∫ x : ℝ in 0..1, quarticCoreMoment24 n x) := by
      rw [← intervalIntegral.integral_neg]
      apply intervalIntegral.integral_congr
      intro x hx
      change
        ‖quarticCoreMoment24 n x‖ =
          -quarticCoreMoment24 n x
      rw [Real.norm_eq_abs, abs_of_nonpos (hnonpos x hx)]
    _ = -quarticCoreSeriesTerm24 n := by
      rw [quarticCoreMoment24_integral]
    _ = ‖quarticCoreSeriesTerm24 n‖ := by
      have hterm : quarticCoreSeriesTerm24 n ≤ 0 := by
        unfold quarticCoreSeriesTerm24
        exact neg_nonpos.mpr
          (div_nonneg
            (mul_nonneg (by norm_num)
              (harmonicNumber_nonneg (n + 1)))
            (by positivity))
      rw [Real.norm_eq_abs, abs_of_nonpos hterm]

private theorem quarticCoreMoment24_integral_norm_summable :
    Summable
      (fun n : ℕ =>
        ∫ x : ℝ in 0..1, ‖quarticCoreMoment24 n x‖) := by
  exact quarticCoreSeriesTerm24_hasSum.summable.norm.congr fun n =>
    (quarticCoreMoment24_integral_norm n).symm

private def quarticCoreKernel24 (x : ℝ) : ℝ :=
  Real.log x * Real.log (1 - x) ^ 2 / x

private theorem quarticCoreMoment24_hasSum_pointwise
    {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasSum (fun n : ℕ => quarticCoreMoment24 n x)
      (quarticCoreKernel24 x) := by
  have hsquare :=
    logOneSubSquare_hasSum24 hx0.le hx1
  have hscaled := hsquare.mul_right (Real.log x / x)
  convert hscaled using 1
  · funext n
    unfold quarticCoreMoment24
    have hxne : x ≠ 0 := ne_of_gt hx0
    rw [show n + 2 = (n + 1) + 1 by omega, pow_succ]
    field_simp [hxne]
    ring
  · unfold quarticCoreKernel24
    ring

private theorem quarticCoreSeriesTerm24_hasSum_integral :
    HasSum quarticCoreSeriesTerm24
      (∫ x : ℝ in 0..1, quarticCoreKernel24 x) := by
  have hInt :
      ∀ n : ℕ,
        MeasureTheory.Integrable
          (quarticCoreMoment24 n)
          (MeasureTheory.volume.restrict (Set.Ioc 0 1)) := by
    intro n
    exact (quarticCoreMoment24_intervalIntegrable n).1
  have hNorm :
      Summable
        (fun n : ℕ =>
          ∫ x : ℝ in Set.Ioc 0 1,
            ‖quarticCoreMoment24 n x‖) := by
    simpa only [
      ← intervalIntegral.integral_of_le
        (by norm_num : (0 : ℝ) ≤ 1)] using
      quarticCoreMoment24_integral_norm_summable
  have h :=
    MeasureTheory.hasSum_integral_of_summable_integral_norm
      (μ := MeasureTheory.volume.restrict (Set.Ioc (0 : ℝ) 1))
      hInt hNorm
  have h' :
      HasSum quarticCoreSeriesTerm24
        (∫ x : ℝ in Set.Ioc 0 1,
          ∑' n : ℕ, quarticCoreMoment24 n x) := by
    convert h using 1
    funext n
    rw [← intervalIntegral.integral_of_le
      (by norm_num : (0 : ℝ) ≤ 1)]
    exact (quarticCoreMoment24_integral n).symm
  convert h' using 1
  rw [intervalIntegral.integral_of_le (by norm_num : (0 : ℝ) ≤ 1)]
  apply MeasureTheory.setIntegral_congr_ae measurableSet_Ioc
  filter_upwards [
    MeasureTheory.Measure.ae_ne MeasureTheory.volume (1 : ℝ)
  ] with x hxne hx
  have hxlt : x < 1 := lt_of_le_of_ne hx.2 hxne
  exact
    (quarticCoreMoment24_hasSum_pointwise hx.1 hxlt).tsum_eq.symm

private theorem quarticCoreIntegral24 :
    (∫ x : ℝ in 0..1, quarticCoreKernel24 x) =
      -(1 / 2 : ℝ) * (Real.pi ^ 4 / 90) := by
  have hEq :=
    quarticCoreSeriesTerm24_hasSum_integral.unique
      quarticCoreSeriesTerm24_hasSum
  rw [hEq]
  ring

private theorem selfMulLogSquare24_continuousOn :
    ContinuousOn (fun x : ℝ => x * Real.log x ^ 2)
      (Icc (0 : ℝ) 1) := by
  intro x hx
  by_cases hxzero : x = 0
  · subst x
    have hright :
        ContinuousWithinAt (fun y : ℝ => y * Real.log y ^ 2)
          (Ioi (0 : ℝ)) 0 := by
      change
        Tendsto (fun y : ℝ => y * Real.log y ^ 2)
          (𝓝[>] (0 : ℝ))
          (𝓝 ((0 : ℝ) * Real.log 0 ^ 2))
      simpa using tendsto_self_mul_log_sq_nhdsGT_zero24
    exact
      (continuousWithinAt_Ioi_iff_Ici.mp hright).mono
        (fun y hy => hy.1)
  · have hlog : ContinuousAt Real.log x :=
      Real.continuousAt_log hxzero
    exact (continuousAt_id.mul (hlog.pow 2)).continuousWithinAt

private theorem quarticCoreKernel24_continuousOn_half :
    ContinuousOn quarticCoreKernel24
      (Icc (0 : ℝ) (1 / 2)) := by
  have haux :
      ContinuousOn
        (fun x : ℝ =>
          (x * Real.log x) * logOneMinusSlope24 x ^ 2)
        (Icc (0 : ℝ) (1 / 2)) :=
    Real.continuous_mul_log.continuousOn.mul
      (logOneMinusSlope24_continuousOn.pow 2)
  apply haux.congr
  intro x _
  by_cases hxzero : x = 0
  · subst x
    simp [quarticCoreKernel24]
  · simp [quarticCoreKernel24, logOneMinusSlope24, hxzero]
    field_simp [hxzero]

private theorem quarticCoreKernel24_intervalIntegrable_half :
    IntervalIntegrable quarticCoreKernel24
      MeasureTheory.volume 0 (1 / 2) := by
  apply ContinuousOn.intervalIntegrable
  rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2)]
  exact quarticCoreKernel24_continuousOn_half

private theorem quarticCoreKernel24_continuousOn_upper :
    ContinuousOn quarticCoreKernel24
      (Icc (1 / 2 : ℝ) 1) := by
  have hsub :
      ContinuousOn (fun x : ℝ => 1 - x)
        (Icc (1 / 2 : ℝ) 1) :=
    continuousOn_const.sub continuousOn_id
  have hmap :
      MapsTo (fun x : ℝ => 1 - x)
        (Icc (1 / 2 : ℝ) 1) (Icc (0 : ℝ) 1) := by
    intro x hx
    constructor <;> linarith [hx.1, hx.2]
  have hmapHalf :
      MapsTo (fun x : ℝ => 1 - x)
        (Icc (1 / 2 : ℝ) 1) (Icc (0 : ℝ) (1 / 2)) := by
    intro x hx
    constructor <;> linarith [hx.1, hx.2]
  have hslope :
      ContinuousOn (fun x : ℝ => logOneMinusSlope24 (1 - x))
        (Icc (1 / 2 : ℝ) 1) :=
    logOneMinusSlope24_continuousOn.comp hsub hmapHalf
  have hsmall :
      ContinuousOn
        (fun x : ℝ => (1 - x) * Real.log (1 - x) ^ 2)
        (Icc (1 / 2 : ℝ) 1) :=
    selfMulLogSquare24_continuousOn.comp hsub hmap
  have haux :
      ContinuousOn
        (fun x : ℝ =>
          logOneMinusSlope24 (1 - x) *
            ((1 - x) * Real.log (1 - x) ^ 2) / x)
        (Icc (1 / 2 : ℝ) 1) := by
    apply ContinuousOn.div (hslope.mul hsmall) continuousOn_id
    intro x hx
    exact ne_of_gt (show 0 < x by linarith [hx.1])
  apply haux.congr
  intro x hx
  by_cases hxone : x = 1
  · subst x
    simp [quarticCoreKernel24, logOneMinusSlope24]
  · have hxne : x ≠ 0 := ne_of_gt (show 0 < x by linarith [hx.1])
    have h1xne : 1 - x ≠ 0 :=
      sub_ne_zero.mpr (Ne.symm hxone)
    simp [quarticCoreKernel24, logOneMinusSlope24, h1xne]
    field_simp [hxne, h1xne]

private theorem quarticCoreKernel24_intervalIntegrable_upper :
    IntervalIntegrable quarticCoreKernel24
      MeasureTheory.volume (1 / 2) 1 := by
  apply ContinuousOn.intervalIntegrable
  rw [Set.uIcc_of_le (by norm_num : (1 / 2 : ℝ) ≤ 1)]
  exact quarticCoreKernel24_continuousOn_upper

private def quarticComplementKernel24 (x : ℝ) : ℝ :=
  Real.log x ^ 2 * Real.log (1 - x) / (1 - x)

private theorem quarticComplementKernel24_intervalIntegrable :
    IntervalIntegrable quarticComplementKernel24
      MeasureTheory.volume 0 (1 / 2) := by
  have hreflected :
      IntervalIntegrable
        (fun x : ℝ => quarticCoreKernel24 (1 - x))
        MeasureTheory.volume 0 (1 / 2) := by
    have h :=
      quarticCoreKernel24_intervalIntegrable_upper.comp_sub_left 1
    convert h.symm using 1 <;> norm_num
  apply IntervalIntegrable.congr
    (f := fun x : ℝ => quarticCoreKernel24 (1 - x)) ?_ hreflected
  intro x _
  unfold quarticComplementKernel24 quarticCoreKernel24
  change
    Real.log (1 - x) * Real.log (1 - (1 - x)) ^ 2 / (1 - x) =
      Real.log x ^ 2 * Real.log (1 - x) / (1 - x)
  rw [show 1 - (1 - x) = x by ring]
  ring

private theorem quarticCoreUpperIntegral24 :
    (∫ x : ℝ in (1 / 2)..1, quarticCoreKernel24 x) =
      ∫ x : ℝ in 0..(1 / 2), quarticComplementKernel24 x := by
  calc
    (∫ x : ℝ in (1 / 2)..1, quarticCoreKernel24 x) =
        ∫ x : ℝ in 0..(1 / 2),
          quarticCoreKernel24 (1 - x) := by
      have h :=
        intervalIntegral.integral_comp_sub_left
          (a := (0 : ℝ)) (b := (1 / 2 : ℝ))
          quarticCoreKernel24 1
      convert h.symm using 1 <;> norm_num
    _ = ∫ x : ℝ in 0..(1 / 2), quarticComplementKernel24 x := by
      apply intervalIntegral.integral_congr
      intro x _
      unfold quarticComplementKernel24 quarticCoreKernel24
      change
        Real.log (1 - x) * Real.log (1 - (1 - x)) ^ 2 / (1 - x) =
          Real.log x ^ 2 * Real.log (1 - x) / (1 - x)
      rw [show 1 - (1 - x) = x by ring]
      ring

private def quarticCoreDifferencePrimitive24 (x : ℝ) : ℝ :=
  -(1 / 2 : ℝ) * Real.log x ^ 2 * Real.log (1 - x) ^ 2

private theorem quarticCoreDifferencePrimitive24_hasDerivAt
    {x : ℝ} (hx0 : 0 < x) (hxhalf : x < 1 / 2) :
    HasDerivAt quarticCoreDifferencePrimitive24
      (quarticComplementKernel24 x - quarticCoreKernel24 x) x := by
  have hxne : x ≠ 0 := ne_of_gt hx0
  have h1xne : 1 - x ≠ 0 := by
    exact ne_of_gt (show 0 < 1 - x by linarith)
  have hlogx := Real.hasDerivAt_log hxne
  have hsub :
      HasDerivAt (fun y : ℝ => 1 - y) (-1) x := by
    convert (hasDerivAt_const x 1).sub (hasDerivAt_id x) using 1
    norm_num
  have hlogsub :
      HasDerivAt (fun y : ℝ => Real.log (1 - y))
        (-1 / (1 - x)) x := by
    convert hsub.log h1xne using 1
  unfold quarticCoreDifferencePrimitive24
  convert
    ((hlogx.pow 2).mul (hlogsub.pow 2)).const_mul
      (-(1 / 2 : ℝ)) using 1
  · funext y
    simp only [Pi.mul_apply, Pi.pow_apply]
    ring
  · unfold quarticComplementKernel24 quarticCoreKernel24
    simp only [Pi.pow_apply]
    field_simp [hxne, h1xne]
    ring

private theorem quarticCoreDifferencePrimitive24_continuousOn :
    ContinuousOn quarticCoreDifferencePrimitive24
      (Icc (0 : ℝ) (1 / 2)) := by
  have haux :
      ContinuousOn
        (fun x : ℝ =>
          -(1 / 2 : ℝ) *
            ((x * Real.log x ^ 2) * x *
              logOneMinusSlope24 x ^ 2))
        (Icc (0 : ℝ) (1 / 2)) := by
    have hself :
        ContinuousOn (fun x : ℝ => x * Real.log x ^ 2)
          (Icc (0 : ℝ) (1 / 2)) :=
      selfMulLogSquare24_continuousOn.mono
        (fun x hx => ⟨hx.1, by linarith [hx.2]⟩)
    exact
      continuousOn_const.mul
        (((hself.mul continuousOn_id).mul
          (logOneMinusSlope24_continuousOn.pow 2)))
  apply haux.congr
  intro x _
  by_cases hxzero : x = 0
  · subst x
    simp [quarticCoreDifferencePrimitive24]
  · simp [quarticCoreDifferencePrimitive24,
      logOneMinusSlope24, hxzero]
    field_simp [hxzero]

private theorem quarticCoreDifferenceIntegral24 :
    (∫ x : ℝ in 0..(1 / 2),
      (quarticComplementKernel24 x - quarticCoreKernel24 x)) =
      -(1 / 2 : ℝ) * Real.log 2 ^ 4 := by
  have hint :=
    quarticComplementKernel24_intervalIntegrable.sub
      quarticCoreKernel24_intervalIntegrable_half
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le
    (f := quarticCoreDifferencePrimitive24)
    (f' := fun x : ℝ =>
      quarticComplementKernel24 x - quarticCoreKernel24 x)
    (by norm_num)
    quarticCoreDifferencePrimitive24_continuousOn
    (fun x hx =>
      quarticCoreDifferencePrimitive24_hasDerivAt
        hx.1 (by linarith [hx.2]))
    hint]
  have hloghalf : Real.log (1 / 2 : ℝ) = -Real.log 2 := by
    rw [one_div, Real.log_inv]
  norm_num [quarticCoreDifferencePrimitive24]
  rw [hloghalf]
  ring

private theorem quarticCoreHalfIntegral24 :
    (∫ x : ℝ in 0..(1 / 2), quarticCoreKernel24 x) =
      (1 / 4 : ℝ) * Real.log 2 ^ 4 -
        (1 / 4 : ℝ) * (Real.pi ^ 4 / 90) := by
  have hsplit :=
    intervalIntegral.integral_add_adjacent_intervals
      quarticCoreKernel24_intervalIntegrable_half
      quarticCoreKernel24_intervalIntegrable_upper
  rw [quarticCoreUpperIntegral24, quarticCoreIntegral24] at hsplit
  have hdiff := quarticCoreDifferenceIntegral24
  rw [intervalIntegral.integral_sub
    quarticComplementKernel24_intervalIntegrable
    quarticCoreKernel24_intervalIntegrable_half] at hdiff
  linarith

private def quarticPlusKernel24 (x : ℝ) : ℝ :=
  Real.log x ^ 2 * Real.log (1 + x) / (1 + x)

private def quarticPlusRadialKernel24 (x : ℝ) : ℝ :=
  Real.log x * Real.log (1 + x) ^ 2 / x

private theorem quarticPlusKernel24_continuousOn :
    ContinuousOn quarticPlusKernel24 (Icc (0 : ℝ) 1) := by
  have hden :
      ContinuousOn (fun x : ℝ => 1 + x) (Icc (0 : ℝ) 1) := by
    fun_prop
  have haux :
      ContinuousOn
        (fun x : ℝ =>
          (x * Real.log x ^ 2) *
            RamanujanChallenge.P26.logOnePlusSlope26 x / (1 + x))
        (Icc (0 : ℝ) 1) := by
    apply ContinuousOn.div
      (selfMulLogSquare24_continuousOn.mul
        RamanujanChallenge.P26.logOnePlusSlope26_continuousOn)
      hden
    intro x hx
    exact ne_of_gt (show 0 < 1 + x by linarith [hx.1])
  apply haux.congr
  intro x _
  by_cases hxzero : x = 0
  · subst x
    simp [quarticPlusKernel24]
  · simp [quarticPlusKernel24,
      RamanujanChallenge.P26.logOnePlusSlope26, hxzero]
    field_simp [hxzero]

private theorem quarticPlusKernel24_intervalIntegrable :
    IntervalIntegrable quarticPlusKernel24
      MeasureTheory.volume 0 1 := by
  apply ContinuousOn.intervalIntegrable
  rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
  exact quarticPlusKernel24_continuousOn

private theorem quarticPlusRadialKernel24_continuousOn :
    ContinuousOn quarticPlusRadialKernel24 (Icc (0 : ℝ) 1) := by
  have haux :
      ContinuousOn
        (fun x : ℝ =>
          (x * Real.log x) *
            RamanujanChallenge.P26.logOnePlusSlope26 x ^ 2)
        (Icc (0 : ℝ) 1) :=
    Real.continuous_mul_log.continuousOn.mul
      (RamanujanChallenge.P26.logOnePlusSlope26_continuousOn.pow 2)
  apply haux.congr
  intro x _
  by_cases hxzero : x = 0
  · subst x
    simp [quarticPlusRadialKernel24]
  · simp [quarticPlusRadialKernel24,
      RamanujanChallenge.P26.logOnePlusSlope26, hxzero]
    field_simp [hxzero]

private theorem quarticPlusRadialKernel24_intervalIntegrable :
    IntervalIntegrable quarticPlusRadialKernel24
      MeasureTheory.volume 0 1 := by
  apply ContinuousOn.intervalIntegrable
  rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
  exact quarticPlusRadialKernel24_continuousOn

private def quarticPlusIBPPrimitive24 (x : ℝ) : ℝ :=
  (1 / 2 : ℝ) * Real.log x ^ 2 * Real.log (1 + x) ^ 2

private theorem quarticPlusIBPPrimitive24_hasDerivAt
    {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasDerivAt quarticPlusIBPPrimitive24
      (quarticPlusKernel24 x + quarticPlusRadialKernel24 x) x := by
  have hxne : x ≠ 0 := ne_of_gt hx0
  have hplusne : 1 + x ≠ 0 := by linarith
  have hlogx := Real.hasDerivAt_log hxne
  have hinner :
      HasDerivAt (fun y : ℝ => 1 + y) 1 x := by
    convert (hasDerivAt_const x 1).add (hasDerivAt_id x) using 1
    norm_num
  have hlogplus :
      HasDerivAt (fun y : ℝ => Real.log (1 + y))
        (1 / (1 + x)) x := by
    convert hinner.log hplusne using 1
  unfold quarticPlusIBPPrimitive24
  convert
    ((hlogx.pow 2).mul (hlogplus.pow 2)).const_mul
      (1 / 2 : ℝ) using 1
  · funext y
    simp only [Pi.mul_apply, Pi.pow_apply]
    ring
  · unfold quarticPlusKernel24 quarticPlusRadialKernel24
    simp only [Pi.pow_apply]
    field_simp [hxne, hplusne]
    ring

private theorem quarticPlusIBPPrimitive24_continuousOn :
    ContinuousOn quarticPlusIBPPrimitive24 (Icc (0 : ℝ) 1) := by
  have haux :
      ContinuousOn
        (fun x : ℝ =>
          (1 / 2 : ℝ) *
            ((x * Real.log x ^ 2) * x *
              RamanujanChallenge.P26.logOnePlusSlope26 x ^ 2))
        (Icc (0 : ℝ) 1) :=
    continuousOn_const.mul
      ((selfMulLogSquare24_continuousOn.mul continuousOn_id).mul
        (RamanujanChallenge.P26.logOnePlusSlope26_continuousOn.pow 2))
  apply haux.congr
  intro x _
  by_cases hxzero : x = 0
  · subst x
    simp [quarticPlusIBPPrimitive24]
  · simp [quarticPlusIBPPrimitive24,
      RamanujanChallenge.P26.logOnePlusSlope26, hxzero]
    field_simp [hxzero]

private theorem quarticPlusIntegral_eq_neg_radial24 :
    (∫ x : ℝ in 0..1, quarticPlusKernel24 x) =
      -(∫ x : ℝ in 0..1, quarticPlusRadialKernel24 x) := by
  have hsum :=
    quarticPlusKernel24_intervalIntegrable.add
      quarticPlusRadialKernel24_intervalIntegrable
  have hzero :
      (∫ x : ℝ in 0..1,
        (quarticPlusKernel24 x + quarticPlusRadialKernel24 x)) = 0 := by
    rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le
      (f := quarticPlusIBPPrimitive24)
      (f' := fun x : ℝ =>
        quarticPlusKernel24 x + quarticPlusRadialKernel24 x)
      (by norm_num)
      quarticPlusIBPPrimitive24_continuousOn
      (fun x hx =>
        quarticPlusIBPPrimitive24_hasDerivAt hx.1 hx.2)
      hsum]
    norm_num [quarticPlusIBPPrimitive24]
  rw [intervalIntegral.integral_add
    quarticPlusKernel24_intervalIntegrable
    quarticPlusRadialKernel24_intervalIntegrable] at hzero
  linarith

private def quarticMixedHalfKernel24 (x : ℝ) : ℝ :=
  Real.log x * Real.log (1 - x) ^ 2 / (1 - x)

private theorem quarticMixedHalfKernel24_intervalIntegrable :
    IntervalIntegrable quarticMixedHalfKernel24
      MeasureTheory.volume 0 (1 / 2) := by
  have hlogSub :
      ContinuousOn (fun x : ℝ => Real.log (1 - x))
        (Icc (0 : ℝ) (1 / 2)) := by
    apply (continuousOn_const.sub continuousOn_id).log
    intro x hx
    exact ne_of_gt (show 0 < 1 - x by linarith [hx.2])
  have hnum :
      ContinuousOn
        (fun x : ℝ =>
          (Real.log x * Real.log (1 - x) ^ 2) *
            Real.log (1 - x))
        (Icc (0 : ℝ) (1 / 2)) :=
    logMulSquareOneMinus24_continuousOn.mul hlogSub
  have hden :
      ContinuousOn (fun x : ℝ => 1 - x)
        (Icc (0 : ℝ) (1 / 2)) := by
    fun_prop
  have hcont :
      ContinuousOn quarticMixedHalfKernel24
        (Icc (0 : ℝ) (1 / 2)) := by
    unfold quarticMixedHalfKernel24
    apply ContinuousOn.div
      logMulSquareOneMinus24_continuousOn hden
    intro x hx
    exact ne_of_gt (show 0 < 1 - x by linarith [hx.2])
  apply ContinuousOn.intervalIntegrable
  rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2)]
  exact hcont

private def quarticMixedHalfPrimitive24 (x : ℝ) : ℝ :=
  -(1 / 3 : ℝ) * Real.log x * Real.log (1 - x) ^ 3

private theorem quarticMixedHalfPrimitive24_hasDerivAt
    {x : ℝ} (hx0 : 0 < x) (hxhalf : x < 1 / 2) :
    HasDerivAt quarticMixedHalfPrimitive24
      (quarticMixedHalfKernel24 x -
        (1 / 3 : ℝ) * halfLogCubeOneSubKernel24 x) x := by
  have hxne : x ≠ 0 := ne_of_gt hx0
  have h1xne : 1 - x ≠ 0 :=
    ne_of_gt (show 0 < 1 - x by linarith)
  have hlogx := Real.hasDerivAt_log hxne
  have hsub :
      HasDerivAt (fun y : ℝ => 1 - y) (-1) x := by
    convert (hasDerivAt_const x 1).sub (hasDerivAt_id x) using 1
    norm_num
  have hlogsub :
      HasDerivAt (fun y : ℝ => Real.log (1 - y))
        (-1 / (1 - x)) x := by
    convert hsub.log h1xne using 1
  unfold quarticMixedHalfPrimitive24
  convert
    (hlogx.mul (hlogsub.pow 3)).const_mul
      (-(1 / 3 : ℝ)) using 1
  · funext y
    simp only [Pi.mul_apply, Pi.pow_apply]
    ring
  · unfold quarticMixedHalfKernel24 halfLogCubeOneSubKernel24
    simp only [Pi.pow_apply]
    field_simp [hxne, h1xne]
    ring

private theorem quarticMixedHalfPrimitive24_continuousOn :
    ContinuousOn quarticMixedHalfPrimitive24
      (Icc (0 : ℝ) (1 / 2)) := by
  have hlogSub :
      ContinuousOn (fun x : ℝ => Real.log (1 - x))
        (Icc (0 : ℝ) (1 / 2)) := by
    apply (continuousOn_const.sub continuousOn_id).log
    intro x hx
    exact ne_of_gt (show 0 < 1 - x by linarith [hx.2])
  have hc :
      ContinuousOn
        (fun x : ℝ =>
          -(1 / 3 : ℝ) *
            ((Real.log x * Real.log (1 - x) ^ 2) *
              Real.log (1 - x)))
        (Icc (0 : ℝ) (1 / 2)) :=
    continuousOn_const.mul
      (logMulSquareOneMinus24_continuousOn.mul hlogSub)
  unfold quarticMixedHalfPrimitive24
  convert hc using 1
  funext x
  ring

private theorem quarticMixedHalfIntegral24 :
    (∫ x : ℝ in 0..(1 / 2), quarticMixedHalfKernel24 x) =
      (1 / 3 : ℝ) *
          (∫ x : ℝ in 0..(1 / 2), halfLogCubeOneSubKernel24 x) -
        (1 / 3 : ℝ) * Real.log 2 ^ 4 := by
  have hscaled :=
    halfLogCubeOneSubKernel24_intervalIntegrable.const_mul (1 / 3 : ℝ)
  have hdiff := quarticMixedHalfKernel24_intervalIntegrable.sub hscaled
  have hval :
      (∫ x : ℝ in 0..(1 / 2),
        (quarticMixedHalfKernel24 x -
          (1 / 3 : ℝ) * halfLogCubeOneSubKernel24 x)) =
        -(1 / 3 : ℝ) * Real.log 2 ^ 4 := by
    rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le
      (f := quarticMixedHalfPrimitive24)
      (f' := fun x : ℝ =>
        quarticMixedHalfKernel24 x -
          (1 / 3 : ℝ) * halfLogCubeOneSubKernel24 x)
      (by norm_num)
      quarticMixedHalfPrimitive24_continuousOn
      (fun x hx =>
        quarticMixedHalfPrimitive24_hasDerivAt
          hx.1 (by linarith [hx.2]))
      hdiff]
    have hloghalf : Real.log (1 / 2 : ℝ) = -Real.log 2 := by
      rw [one_div, Real.log_inv]
    norm_num [quarticMixedHalfPrimitive24]
    rw [hloghalf]
    ring
  rw [intervalIntegral.integral_sub
    quarticMixedHalfKernel24_intervalIntegrable hscaled,
    intervalIntegral.integral_const_mul] at hval
  linarith

private def quarticLogFourthHalfKernel24 (x : ℝ) : ℝ :=
  Real.log (1 - x) ^ 3 / (1 - x)

private theorem quarticLogFourthHalfIntegral24 :
    (∫ x : ℝ in 0..(1 / 2), quarticLogFourthHalfKernel24 x) =
      -(1 / 4 : ℝ) * Real.log 2 ^ 4 := by
  let F : ℝ → ℝ :=
    fun x => -(1 / 4 : ℝ) * Real.log (1 - x) ^ 4
  have hcont :
      ContinuousOn F (Icc (0 : ℝ) (1 / 2)) := by
    dsimp [F]
    apply continuousOn_const.mul
    apply ((continuousOn_const.sub continuousOn_id).log
      (fun x hx => ne_of_gt
        (show 0 < 1 - x by linarith [hx.2]))).pow
  have hderiv :
      ∀ x ∈ Ioo (0 : ℝ) (1 / 2),
        HasDerivAt F (quarticLogFourthHalfKernel24 x) x := by
    intro x hx
    have h1xne : 1 - x ≠ 0 := by
      exact ne_of_gt (show 0 < 1 - x by linarith [hx.2])
    have hsub :
        HasDerivAt (fun y : ℝ => 1 - y) (-1) x := by
      convert (hasDerivAt_const x 1).sub (hasDerivAt_id x) using 1
      norm_num
    have hlog :
        HasDerivAt (fun y : ℝ => Real.log (1 - y))
          (-1 / (1 - x)) x := by
      convert hsub.log h1xne using 1
    dsimp [F]
    convert (hlog.pow 4).const_mul (-(1 / 4 : ℝ)) using 1
    unfold quarticLogFourthHalfKernel24
    field_simp [h1xne]
    ring
  have hint :
      IntervalIntegrable quarticLogFourthHalfKernel24
        MeasureTheory.volume 0 (1 / 2) := by
    apply ContinuousOn.intervalIntegrable
    rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2)]
    unfold quarticLogFourthHalfKernel24
    apply ContinuousOn.div
    · apply ((continuousOn_const.sub continuousOn_id).log
        (fun x hx => ne_of_gt
          (show 0 < 1 - x by linarith [hx.2]))).pow
    · fun_prop
    · intro x hx
      exact ne_of_gt (show 0 < 1 - x by linarith [hx.2])
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le
    (f := F) (f' := quarticLogFourthHalfKernel24)
    (by norm_num) hcont hderiv hint]
  have hloghalf : Real.log (1 / 2 : ℝ) = -Real.log 2 := by
    rw [one_div, Real.log_inv]
  norm_num [F]
  rw [hloghalf]
  ring

private def quarticHalfMobiusMap24 (t : ℝ) : ℝ :=
  t / (1 - t)

private def quarticHalfMobiusMapDeriv24 (t : ℝ) : ℝ :=
  1 / (1 - t) ^ 2

private theorem quarticHalfMobiusMap24_hasDerivAt
    {t : ℝ} (ht0 : 0 < t) (hthalf : t < 1 / 2) :
    HasDerivAt quarticHalfMobiusMap24
      (quarticHalfMobiusMapDeriv24 t) t := by
  have hden : 1 - t ≠ 0 := by linarith
  have hnum : HasDerivAt (fun x : ℝ => x) 1 t := hasDerivAt_id t
  have hdenDeriv :
      HasDerivAt (fun x : ℝ => 1 - x) (-1) t := by
    convert (hasDerivAt_const t 1).sub (hasDerivAt_id t) using 1
    norm_num
  unfold quarticHalfMobiusMap24 quarticHalfMobiusMapDeriv24
  convert hnum.div hdenDeriv hden using 1
  field_simp [hden]
  ring

private def quarticMobiusTransformedKernel24 (t : ℝ) : ℝ :=
  (Real.log t - Real.log (1 - t)) * Real.log (1 - t) ^ 2 /
    (t * (1 - t))

private theorem quarticMobiusChangeIntegrand24
    {t : ℝ} (ht0 : 0 ≤ t) (hthalf : t ≤ 1 / 2) :
    (quarticPlusRadialKernel24 ∘ quarticHalfMobiusMap24) t *
        quarticHalfMobiusMapDeriv24 t =
      quarticMobiusTransformedKernel24 t := by
  rcases ht0.eq_or_lt with rfl | ht0
  · norm_num [quarticPlusRadialKernel24, quarticHalfMobiusMap24,
      quarticHalfMobiusMapDeriv24, quarticMobiusTransformedKernel24]
  have htne : t ≠ 0 := ne_of_gt ht0
  have h1t0 : 0 < 1 - t := by linarith
  have h1tne : 1 - t ≠ 0 := ne_of_gt h1t0
  have hphi0 : 0 < quarticHalfMobiusMap24 t := by
    unfold quarticHalfMobiusMap24
    positivity
  have hlogPhi :
      Real.log (quarticHalfMobiusMap24 t) =
        Real.log t - Real.log (1 - t) := by
    unfold quarticHalfMobiusMap24
    rw [Real.log_div htne h1tne]
  have hOnePlus :
      1 + quarticHalfMobiusMap24 t = 1 / (1 - t) := by
    unfold quarticHalfMobiusMap24
    field_simp [h1tne]
    ring
  have hlogOnePlus :
      Real.log (1 + quarticHalfMobiusMap24 t) =
        -Real.log (1 - t) := by
    rw [hOnePlus, one_div, Real.log_inv]
  simp only [Function.comp_apply]
  unfold quarticPlusRadialKernel24 quarticHalfMobiusMapDeriv24
    quarticMobiusTransformedKernel24
  rw [hlogPhi, hlogOnePlus]
  unfold quarticHalfMobiusMap24
  field_simp [htne, h1tne]

private theorem quarticRadialIntegral_eq_transformed24 :
    (∫ x : ℝ in 0..1, quarticPlusRadialKernel24 x) =
      ∫ t : ℝ in 0..(1 / 2), quarticMobiusTransformedKernel24 t := by
  have hmapCont :
      ContinuousOn quarticHalfMobiusMap24
        (Set.uIcc (0 : ℝ) (1 / 2)) := by
    rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2)]
    unfold quarticHalfMobiusMap24
    apply ContinuousOn.div continuousOn_id
      (continuousOn_const.sub continuousOn_id)
    intro t ht
    exact ne_of_gt (show 0 < 1 - t by linarith [ht.2])
  have hderiv :
      ∀ t ∈ Ioo (min (0 : ℝ) (1 / 2)) (max (0 : ℝ) (1 / 2)),
        HasDerivAt quarticHalfMobiusMap24
          (quarticHalfMobiusMapDeriv24 t) t := by
    intro t ht
    norm_num at ht
    exact quarticHalfMobiusMap24_hasDerivAt ht.1 ht.2
  have hnonneg :
      ∀ t ∈ Ioo (min (0 : ℝ) (1 / 2)) (max (0 : ℝ) (1 / 2)),
        0 ≤ quarticHalfMobiusMapDeriv24 t := by
    intro t _
    unfold quarticHalfMobiusMapDeriv24
    positivity
  have hsubst :=
    intervalIntegral.integral_comp_mul_deriv_of_deriv_nonneg
      (a := (0 : ℝ)) (b := (1 / 2 : ℝ))
      (f := quarticHalfMobiusMap24)
      (f' := quarticHalfMobiusMapDeriv24)
      (g := quarticPlusRadialKernel24)
      hmapCont hderiv hnonneg
  have hs :
      (∫ t : ℝ in 0..(1 / 2),
        (quarticPlusRadialKernel24 ∘ quarticHalfMobiusMap24) t *
          quarticHalfMobiusMapDeriv24 t) =
        ∫ x : ℝ in 0..1, quarticPlusRadialKernel24 x := by
    convert hsubst using 1 <;>
      norm_num [quarticHalfMobiusMap24]
  have hcongr :
      (∫ t : ℝ in 0..(1 / 2),
        (quarticPlusRadialKernel24 ∘ quarticHalfMobiusMap24) t *
          quarticHalfMobiusMapDeriv24 t) =
        ∫ t : ℝ in 0..(1 / 2),
          quarticMobiusTransformedKernel24 t := by
    apply intervalIntegral.integral_congr
    intro t ht
    have ht' : t ∈ Icc (0 : ℝ) (1 / 2) := by
      simpa [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2)] using ht
    exact quarticMobiusChangeIntegrand24 ht'.1 ht'.2
  linarith

private theorem quarticMobiusTransformedIntegral24 :
    (∫ t : ℝ in 0..(1 / 2), quarticMobiusTransformedKernel24 t) =
      (1 / 4 : ℝ) * Real.log 2 ^ 4 -
        (1 / 4 : ℝ) * (Real.pi ^ 4 / 90) -
        (2 / 3 : ℝ) *
          (6 * polylog4 (1 / 2) +
            (1 / 2 : ℝ) * Real.log 2 ^ 4 -
            (3 / 2 : ℝ) * Real.log 2 ^ 2 *
              (Real.pi ^ 2 / 6) +
            (21 / 4 : ℝ) * Real.log 2 * zeta3_24 -
            6 * (Real.pi ^ 4 / 90)) -
        (1 / 12 : ℝ) * Real.log 2 ^ 4 := by
  have hpoint :
      ∀ t ∈ Set.uIcc (0 : ℝ) (1 / 2),
        quarticMobiusTransformedKernel24 t =
          quarticCoreKernel24 t + quarticMixedHalfKernel24 t -
            halfLogCubeOneSubKernel24 t -
              quarticLogFourthHalfKernel24 t := by
    intro t ht
    rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2)] at ht
    by_cases htzero : t = 0
    · subst t
      norm_num [quarticMobiusTransformedKernel24,
        quarticCoreKernel24, quarticMixedHalfKernel24,
        halfLogCubeOneSubKernel24, quarticLogFourthHalfKernel24]
    · have h1tne : 1 - t ≠ 0 :=
        ne_of_gt (show 0 < 1 - t by linarith [ht.2])
      unfold quarticMobiusTransformedKernel24 quarticCoreKernel24
        quarticMixedHalfKernel24 halfLogCubeOneSubKernel24
        quarticLogFourthHalfKernel24
      field_simp [htzero, h1tne]
      ring
  calc
    (∫ t : ℝ in 0..(1 / 2), quarticMobiusTransformedKernel24 t) =
        ∫ t : ℝ in 0..(1 / 2),
          (quarticCoreKernel24 t + quarticMixedHalfKernel24 t -
            halfLogCubeOneSubKernel24 t -
              quarticLogFourthHalfKernel24 t) := by
      apply intervalIntegral.integral_congr
      exact hpoint
    _ = (∫ t : ℝ in 0..(1 / 2), quarticCoreKernel24 t) +
          (∫ t : ℝ in 0..(1 / 2), quarticMixedHalfKernel24 t) -
          (∫ t : ℝ in 0..(1 / 2), halfLogCubeOneSubKernel24 t) -
          (∫ t : ℝ in 0..(1 / 2), quarticLogFourthHalfKernel24 t) := by
      rw [intervalIntegral.integral_sub
        ((quarticCoreKernel24_intervalIntegrable_half.add
          quarticMixedHalfKernel24_intervalIntegrable).sub
            halfLogCubeOneSubKernel24_intervalIntegrable)
        (by
          unfold quarticLogFourthHalfKernel24
          apply ContinuousOn.intervalIntegrable
          rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2)]
          apply ContinuousOn.div
          · apply ((continuousOn_const.sub continuousOn_id).log
              (fun x hx => ne_of_gt
                (show 0 < 1 - x by linarith [hx.2]))).pow
          · fun_prop
          · intro x hx
            exact ne_of_gt (show 0 < 1 - x by linarith [hx.2])),
        intervalIntegral.integral_sub
          (quarticCoreKernel24_intervalIntegrable_half.add
            quarticMixedHalfKernel24_intervalIntegrable)
          halfLogCubeOneSubKernel24_intervalIntegrable,
        intervalIntegral.integral_add
          quarticCoreKernel24_intervalIntegrable_half
          quarticMixedHalfKernel24_intervalIntegrable]
    _ = _ := by
      rw [quarticCoreHalfIntegral24, quarticMixedHalfIntegral24,
        halfLogCubeOneSubIntegral24, quarticLogFourthHalfIntegral24]
      ring

private theorem quarticPlusIntegral24 :
    (∫ x : ℝ in 0..1, quarticPlusKernel24 x) =
      4 * polylog4 (1 / 2) +
        (1 / 6 : ℝ) * Real.log 2 ^ 4 -
        Real.log 2 ^ 2 * (Real.pi ^ 2 / 6) +
        (7 / 2 : ℝ) * Real.log 2 * zeta3_24 -
        (15 / 4 : ℝ) * (Real.pi ^ 4 / 90) := by
  rw [quarticPlusIntegral_eq_neg_radial24,
    quarticRadialIntegral_eq_transformed24,
    quarticMobiusTransformedIntegral24]
  ring

private def powerLogSquarePrimitive24 (n : ℕ) (x : ℝ) : ℝ :=
  x ^ (n + 2) *
    (Real.log x ^ 2 / (n + 2 : ℝ) -
      2 * Real.log x / (n + 2 : ℝ) ^ 2 +
      2 / (n + 2 : ℝ) ^ 3)

private theorem powerLogSquarePrimitive24_hasDerivAt
    (n : ℕ) {x : ℝ} (hx0 : 0 < x) :
    HasDerivAt (powerLogSquarePrimitive24 n)
      (x ^ (n + 1) * Real.log x ^ 2) x := by
  have hxne : x ≠ 0 := ne_of_gt hx0
  have hden : (n + 2 : ℝ) ≠ 0 := by positivity
  have hpow :
      HasDerivAt (fun y : ℝ => y ^ (n + 2))
        ((n + 2 : ℝ) * x ^ (n + 1)) x := by
    convert hasDerivAt_pow (n + 2) x using 1
    norm_num [Nat.cast_add]
  have hlog := Real.hasDerivAt_log hxne
  have hfactor :
      HasDerivAt
        (fun y : ℝ =>
          Real.log y ^ 2 / (n + 2 : ℝ) -
            2 * Real.log y / (n + 2 : ℝ) ^ 2 +
            2 / (n + 2 : ℝ) ^ 3)
        ((2 * Real.log x * x⁻¹) / (n + 2 : ℝ) -
          2 * x⁻¹ / (n + 2 : ℝ) ^ 2) x := by
    have hfirst :=
      (hlog.pow 2).div_const (n + 2 : ℝ)
    have hsecond :=
      (hlog.const_mul 2).div_const ((n + 2 : ℝ) ^ 2)
    have htotal :=
      (hfirst.sub hsecond).add_const (2 / (n + 2 : ℝ) ^ 3)
    convert htotal using 1
    ring
  unfold powerLogSquarePrimitive24
  convert hpow.mul hfactor using 1
  field_simp [hxne, hden]
  ring

private theorem powerLogSquarePrimitive24_continuousOn (n : ℕ) :
    ContinuousOn (powerLogSquarePrimitive24 n)
      (Icc (0 : ℝ) 1) := by
  have hpow :
      ContinuousOn (fun x : ℝ => x ^ (n + 1))
        (Icc (0 : ℝ) 1) :=
    continuousOn_pow (n + 1)
  have hsq :
      ContinuousOn
        (fun x : ℝ => x ^ (n + 2) * Real.log x ^ 2)
        (Icc (0 : ℝ) 1) := by
    have hc :=
      hpow.mul selfMulLogSquare24_continuousOn
    convert hc using 1
    funext x
    simp only [Pi.mul_apply]
    rw [show n + 2 = (n + 1) + 1 by omega, pow_succ]
    ring
  have hlin :
      ContinuousOn
        (fun x : ℝ => x ^ (n + 2) * Real.log x)
        (Icc (0 : ℝ) 1) := by
    have hc :=
      hpow.mul Real.continuous_mul_log.continuousOn
    convert hc using 1
    funext x
    simp only [Pi.mul_apply]
    rw [show n + 2 = (n + 1) + 1 by omega, pow_succ]
    ring
  have hp :
      ContinuousOn (fun x : ℝ => x ^ (n + 2))
        (Icc (0 : ℝ) 1) :=
    continuousOn_pow (n + 2)
  have hc :
      ContinuousOn
        (fun x : ℝ =>
          x ^ (n + 2) * Real.log x ^ 2 / (n + 2 : ℝ) -
            2 * (x ^ (n + 2) * Real.log x) / (n + 2 : ℝ) ^ 2 +
            2 * x ^ (n + 2) / (n + 2 : ℝ) ^ 3)
        (Icc (0 : ℝ) 1) := by
    fun_prop
  unfold powerLogSquarePrimitive24
  apply hc.congr
  intro x _
  ring

private theorem powerLogSquare24_continuousOn (n : ℕ) :
    ContinuousOn (fun x : ℝ => x ^ (n + 1) * Real.log x ^ 2)
      (Icc (0 : ℝ) 1) := by
  have hc :=
    (continuousOn_pow n).mul selfMulLogSquare24_continuousOn
  convert hc using 1
  funext x
  simp only [Pi.mul_apply]
  rw [pow_succ]
  ring

private theorem integral_powerLogSquare24 (n : ℕ) :
    (∫ x : ℝ in 0..1, x ^ (n + 1) * Real.log x ^ 2) =
      2 / (n + 2 : ℝ) ^ 3 := by
  have hint :
      IntervalIntegrable
        (fun x : ℝ => x ^ (n + 1) * Real.log x ^ 2)
        MeasureTheory.volume 0 1 := by
    apply ContinuousOn.intervalIntegrable
    rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
    exact powerLogSquare24_continuousOn n
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le
    (f := powerLogSquarePrimitive24 n)
    (f' := fun x : ℝ => x ^ (n + 1) * Real.log x ^ 2)
    (by norm_num)
    (powerLogSquarePrimitive24_continuousOn n)
    (fun x hx => powerLogSquarePrimitive24_hasDerivAt n hx.1)
    hint]
  norm_num [powerLogSquarePrimitive24]

private def quarticPlusSeriesMoment24 (n : ℕ) (x : ℝ) : ℝ :=
  (-1 : ℝ) ^ n * harmonicNumber (n + 1) *
    (x ^ (n + 1) * Real.log x ^ 2)

private def quarticPlusSeriesTerm24 (n : ℕ) : ℝ :=
  2 * ((-1 : ℝ) ^ n * harmonicNumber (n + 1)) /
    (n + 2 : ℝ) ^ 3

private theorem quarticPlusSeriesMoment24_intervalIntegrable (n : ℕ) :
    IntervalIntegrable (quarticPlusSeriesMoment24 n)
      MeasureTheory.volume 0 1 := by
  unfold quarticPlusSeriesMoment24
  apply ContinuousOn.intervalIntegrable
  rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
  exact
    (powerLogSquare24_continuousOn n).const_mul
      ((-1 : ℝ) ^ n * harmonicNumber (n + 1))

private theorem quarticPlusSeriesMoment24_integral (n : ℕ) :
    (∫ x : ℝ in 0..1, quarticPlusSeriesMoment24 n x) =
      quarticPlusSeriesTerm24 n := by
  unfold quarticPlusSeriesMoment24 quarticPlusSeriesTerm24
  rw [intervalIntegral.integral_const_mul,
    integral_powerLogSquare24]
  ring

private theorem quarticPlusSeriesMoment24_integral_norm (n : ℕ) :
    (∫ x : ℝ in 0..1, ‖quarticPlusSeriesMoment24 n x‖) =
      2 * harmonicNumber (n + 1) / (n + 2 : ℝ) ^ 3 := by
  have hpoint :
      ∀ x ∈ Set.uIcc (0 : ℝ) 1,
        ‖quarticPlusSeriesMoment24 n x‖ =
          harmonicNumber (n + 1) *
            (x ^ (n + 1) * Real.log x ^ 2) := by
    intro x hx
    rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] at hx
    unfold quarticPlusSeriesMoment24
    rw [Real.norm_eq_abs, abs_mul, abs_mul, abs_pow, abs_neg, abs_one,
      one_pow, abs_of_nonneg (harmonicNumber_nonneg (n + 1)),
      abs_of_nonneg (mul_nonneg (pow_nonneg hx.1 _) (sq_nonneg _))]
    ring
  calc
    (∫ x : ℝ in 0..1, ‖quarticPlusSeriesMoment24 n x‖) =
        ∫ x : ℝ in 0..1,
          harmonicNumber (n + 1) *
            (x ^ (n + 1) * Real.log x ^ 2) := by
      apply intervalIntegral.integral_congr
      exact hpoint
    _ = harmonicNumber (n + 1) *
          (∫ x : ℝ in 0..1, x ^ (n + 1) * Real.log x ^ 2) := by
      rw [intervalIntegral.integral_const_mul]
    _ = 2 * harmonicNumber (n + 1) / (n + 2 : ℝ) ^ 3 := by
      rw [integral_powerLogSquare24]
      ring

private theorem quarticPlusSeriesMoment24_integral_norm_summable :
    Summable
      (fun n : ℕ =>
        ∫ x : ℝ in 0..1, ‖quarticPlusSeriesMoment24 n x‖) := by
  have hmajor :
      Summable (fun n : ℕ => 2 * ordinaryHarmonicCubicTerm24 n) :=
    summable_ordinaryHarmonicCubicTerm24.mul_left 2
  have hsmall :
      Summable
        (fun n : ℕ =>
          2 * harmonicNumber (n + 1) / (n + 2 : ℝ) ^ 3) := by
    apply hmajor.of_nonneg_of_le
    · intro n
      exact div_nonneg
        (mul_nonneg (by norm_num) (harmonicNumber_nonneg (n + 1)))
        (pow_nonneg (by positivity) _)
    · intro n
      unfold ordinaryHarmonicCubicTerm24
      rw [show
        2 * (harmonicNumber (n + 1) / (n + 1 : ℝ) ^ 3) =
          2 * harmonicNumber (n + 1) / (n + 1 : ℝ) ^ 3 by ring]
      rw [div_le_div_iff₀ (by positivity) (by positivity)]
      gcongr
      · exact mul_nonneg (by norm_num) (harmonicNumber_nonneg (n + 1))
      · norm_num
  exact hsmall.congr fun n =>
    (quarticPlusSeriesMoment24_integral_norm n).symm

private theorem quarticPlusSeriesMoment24_hasSum_pointwise
    {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasSum (fun n : ℕ => quarticPlusSeriesMoment24 n x)
      (quarticPlusKernel24 x) := by
  have hxabs : |-x| < 1 := by
    rw [abs_neg, abs_of_pos hx0]
    exact hx1
  have hgen :=
    harmonicNumber_generating_hasSum hxabs (neg_ne_zero.mpr (ne_of_gt hx0))
  have hscaled := hgen.mul_right (x * Real.log x ^ 2)
  convert hscaled using 1
  · funext n
    unfold quarticPlusSeriesMoment24
    rw [neg_pow, pow_succ]
    ring
  · unfold quarticPlusKernel24
    have hxne : x ≠ 0 := ne_of_gt hx0
    have hplusne : 1 + x ≠ 0 := by linarith
    rw [show 1 - -x = 1 + x by ring]
    field_simp [hxne, hplusne]

private theorem quarticPlusSeriesTerm24_hasSum_integral :
    HasSum quarticPlusSeriesTerm24
      (∫ x : ℝ in 0..1, quarticPlusKernel24 x) := by
  have hInt :
      ∀ n : ℕ,
        MeasureTheory.Integrable
          (quarticPlusSeriesMoment24 n)
          (MeasureTheory.volume.restrict (Set.Ioc 0 1)) := by
    intro n
    exact (quarticPlusSeriesMoment24_intervalIntegrable n).1
  have hNorm :
      Summable
        (fun n : ℕ =>
          ∫ x : ℝ in Set.Ioc 0 1,
            ‖quarticPlusSeriesMoment24 n x‖) := by
    simpa only [
      ← intervalIntegral.integral_of_le
        (by norm_num : (0 : ℝ) ≤ 1)] using
      quarticPlusSeriesMoment24_integral_norm_summable
  have h :=
    MeasureTheory.hasSum_integral_of_summable_integral_norm
      (μ := MeasureTheory.volume.restrict (Set.Ioc (0 : ℝ) 1))
      hInt hNorm
  have h' :
      HasSum quarticPlusSeriesTerm24
        (∫ x : ℝ in Set.Ioc 0 1,
          ∑' n : ℕ, quarticPlusSeriesMoment24 n x) := by
    convert h using 1
    funext n
    rw [← intervalIntegral.integral_of_le
      (by norm_num : (0 : ℝ) ≤ 1)]
    exact (quarticPlusSeriesMoment24_integral n).symm
  convert h' using 1
  rw [intervalIntegral.integral_of_le (by norm_num : (0 : ℝ) ≤ 1)]
  apply MeasureTheory.setIntegral_congr_ae measurableSet_Ioc
  filter_upwards [
    MeasureTheory.Measure.ae_ne MeasureTheory.volume (1 : ℝ)
  ] with x hxne hx
  exact
    (quarticPlusSeriesMoment24_hasSum_pointwise
      hx.1 (lt_of_le_of_ne hx.2 hxne)).tsum_eq.symm

private theorem quarticPlusSeriesTerm24_hasSum :
    HasSum quarticPlusSeriesTerm24
      (4 * polylog4 (1 / 2) +
        (1 / 6 : ℝ) * Real.log 2 ^ 4 -
        Real.log 2 ^ 2 * (Real.pi ^ 2 / 6) +
        (7 / 2 : ℝ) * Real.log 2 * zeta3_24 -
        (15 / 4 : ℝ) * (Real.pi ^ 4 / 90)) := by
  rw [← quarticPlusIntegral24]
  exact quarticPlusSeriesTerm24_hasSum_integral

private def zetaFourTerm24 (n : ℕ) : ℝ :=
  1 / (n + 1 : ℝ) ^ 4

private theorem zetaFourTerm24_hasSum :
    HasSum zetaFourTerm24 (Real.pi ^ 4 / 90) := by
  exact shifted_zeta_four_hasSum24

private theorem zetaFourTerm24_summable :
    Summable zetaFourTerm24 :=
  zetaFourTerm24_hasSum.summable

private def alternatingZetaFourTerm24 (n : ℕ) : ℝ :=
  (-1 : ℝ) ^ (n + 1) * zetaFourTerm24 n

private theorem alternatingZetaFourTerm24_hasSum :
    HasSum alternatingZetaFourTerm24
      (-(7 / 8 : ℝ) * (Real.pi ^ 4 / 90)) := by
  have hOdd :
      HasSum (fun k : ℕ => zetaFourTerm24 (2 * k + 1))
        ((1 / 16 : ℝ) * (Real.pi ^ 4 / 90)) := by
    convert zetaFourTerm24_hasSum.mul_left (1 / 16 : ℝ) using 1
    · funext k
      unfold zetaFourTerm24
      push_cast
      have hk : (k : ℝ) + 1 ≠ 0 := by positivity
      field_simp [hk]
      ring
  have hEvenSummable :
      Summable (fun k : ℕ => zetaFourTerm24 (2 * k)) :=
    zetaFourTerm24_summable.comp_injective
      (mul_right_injective₀ (by omega : (2 : ℕ) ≠ 0))
  have hEvenTsum :
      (∑' k : ℕ, zetaFourTerm24 (2 * k)) =
        Real.pi ^ 4 / 90 -
          (1 / 16 : ℝ) * (Real.pi ^ 4 / 90) := by
    have hsplit := tsum_even_add_odd hEvenSummable hOdd.summable
    rw [hOdd.tsum_eq, zetaFourTerm24_hasSum.tsum_eq] at hsplit
    linarith
  have hEven :
      HasSum (fun k : ℕ => zetaFourTerm24 (2 * k))
        (Real.pi ^ 4 / 90 -
          (1 / 16 : ℝ) * (Real.pi ^ 4 / 90)) :=
    hEvenSummable.hasSum_iff.mpr hEvenTsum
  have hAltEven :
      HasSum
        (fun k : ℕ => alternatingZetaFourTerm24 (2 * k))
        (-(Real.pi ^ 4 / 90 -
          (1 / 16 : ℝ) * (Real.pi ^ 4 / 90))) := by
    convert hEven.neg using 1
    funext k
    simp [alternatingZetaFourTerm24, pow_add, pow_mul]
  have hAltOdd :
      HasSum
        (fun k : ℕ => alternatingZetaFourTerm24 (2 * k + 1))
        ((1 / 16 : ℝ) * (Real.pi ^ 4 / 90)) := by
    convert hOdd using 1
    funext k
    simp [alternatingZetaFourTerm24, pow_add, pow_mul]
  convert HasSum.even_add_odd hAltEven hAltOdd using 1
  ring

private theorem alternatingOrdinaryHarmonicCubicTerm24_hasSum :
    HasSum alternatingOrdinaryHarmonicCubicTerm24
      (2 * polylog4 (1 / 2) +
        (1 / 12 : ℝ) * Real.log 2 ^ 4 -
        (1 / 2 : ℝ) * Real.log 2 ^ 2 * (Real.pi ^ 2 / 6) +
        (7 / 4 : ℝ) * Real.log 2 * zeta3_24 -
        (11 / 4 : ℝ) * (Real.pi ^ 4 / 90)) := by
  have hZetaTail :
      HasSum
        (fun n : ℕ => alternatingZetaFourTerm24 (n + 1))
        (1 - (7 / 8 : ℝ) * (Real.pi ^ 4 / 90)) := by
    apply
      (hasSum_nat_add_iff
        (f := alternatingZetaFourTerm24) 1).mpr
    convert alternatingZetaFourTerm24_hasSum using 1 <;>
      simp [alternatingZetaFourTerm24, zetaFourTerm24] <;> ring
  have hTail :=
    (quarticPlusSeriesTerm24_hasSum.mul_left (1 / 2 : ℝ)).add
      hZetaTail
  have hTail' :
      HasSum
        (fun n : ℕ =>
          alternatingOrdinaryHarmonicCubicTerm24 (n + 1))
        ((1 / 2 : ℝ) *
            (4 * polylog4 (1 / 2) +
              (1 / 6 : ℝ) * Real.log 2 ^ 4 -
              Real.log 2 ^ 2 * (Real.pi ^ 2 / 6) +
              (7 / 2 : ℝ) * Real.log 2 * zeta3_24 -
              (15 / 4 : ℝ) * (Real.pi ^ 4 / 90)) +
          (1 - (7 / 8 : ℝ) * (Real.pi ^ 4 / 90))) := by
    convert hTail using 1
    funext n
    unfold alternatingOrdinaryHarmonicCubicTerm24
      ordinaryHarmonicCubicTerm24 quarticPlusSeriesTerm24
      alternatingZetaFourTerm24 zetaFourTerm24
    rw [harmonicNumber_succ]
    simp only [Nat.cast_add, Nat.cast_one]
    have hden : (n + 2 : ℝ) ≠ 0 := by positivity
    field_simp [hden]
    simp only [pow_succ]
    ring
  have hFull :=
    (hasSum_nat_add_iff
      (f := alternatingOrdinaryHarmonicCubicTerm24) 1).mp hTail'
  convert hFull using 1 <;>
    simp [alternatingOrdinaryHarmonicCubicTerm24,
      ordinaryHarmonicCubicTerm24, harmonicNumber_one] <;> ring

theorem alternatingCubicLinearEulerTerm24_hasSum :
    HasSum alternatingCubicLinearEulerTerm24
      alternatingCubicLinearEulerValue24 := by
  have hDiffRaw :=
    alternatingOrdinaryHarmonicCubicTerm24_hasSum_integral.sub
      alternatingBarHarmonicCubicTerm24_hasSum_integral
  have hDiff :
      HasSum
        (fun n : ℕ =>
          alternatingOrdinaryHarmonicCubicTerm24 n -
            alternatingBarHarmonicCubicTerm24 n)
        (-(1 / 2 : ℝ) * (Real.pi ^ 2 / 6) ^ 2 +
          (7 / 4 : ℝ) * Real.log 2 * zeta3_24) := by
    convert hDiffRaw using 1
    rw [← intervalIntegral.integral_sub
      alternatingOrdinaryHarmonicCubicKernel24_intervalIntegrable
      alternatingBarHarmonicCubicKernel24_intervalIntegrable,
      alternatingCubicKernelDifferenceIntegral24]
  have hCombined :=
    (hDiff.mul_left 2).sub
      alternatingOrdinaryHarmonicCubicTerm24_hasSum
  convert hCombined using 1
  · funext n
    unfold alternatingCubicLinearEulerTerm24 cubicLinearEulerTerm24
      parityRemainder24 alternatingOrdinaryHarmonicCubicTerm24
      ordinaryHarmonicCubicTerm24 alternatingBarHarmonicCubicTerm24
      barHarmonicCubicTerm24
    ring
  · unfold alternatingCubicLinearEulerValue24
    ring

/-! ## Public re-statements for Layer E of `Problem24QuadraticAlt`

The two alternating weight-four series below are already proved above, but as
private declarations phrased through private abbreviations.  Layer E needs them
with their summands written out, so they are restated here rather than
re-proved. -/

/-- `∑_{n≥0} (-1)^{n+1} H_{n+1}/(n+1)³` in closed form. -/
theorem alternatingHarmonicCubic_hasSum24 :
    HasSum (fun n : ℕ => (-1 : ℝ) ^ (n + 1) * harmonicNumber (n + 1) / ((n : ℝ) + 1) ^ 3)
      (2 * polylog4 (1 / 2) + (1 / 12 : ℝ) * Real.log 2 ^ 4 -
        (1 / 2 : ℝ) * Real.log 2 ^ 2 * (Real.pi ^ 2 / 6) +
        (7 / 4 : ℝ) * Real.log 2 * zeta3_24 -
        (11 / 4 : ℝ) * (Real.pi ^ 4 / 90)) := by
  convert alternatingOrdinaryHarmonicCubicTerm24_hasSum using 1
  funext n
  unfold alternatingOrdinaryHarmonicCubicTerm24 ordinaryHarmonicCubicTerm24
  push_cast
  ring

/-- `∑_{n≥0} (-1)^{n+1}/(n+1)⁴ = -(7/8) ζ(4)`. -/
theorem alternatingZetaFour_hasSum24 :
    HasSum (fun n : ℕ => (-1 : ℝ) ^ (n + 1) / ((n : ℝ) + 1) ^ 4)
      (-(7 / 8 : ℝ) * (Real.pi ^ 4 / 90)) := by
  convert alternatingZetaFourTerm24_hasSum using 1
  funext n
  unfold alternatingZetaFourTerm24 zetaFourTerm24
  push_cast
  ring

/-! ## Public quartic endpoint integrals for `Problem24QuadraticAlt`

The reflected `I11` kernel uses the ordinary quartic core together with the
alternating complement below.  The core was already evaluated privately above;
the complement is obtained by the same termwise-integration template, using
`-signedHarmonic24 (n+1)` as its nonnegative Taylor coefficient. -/

/-- `∫₀¹ log x · log²(1-x) / x = -ζ(4)/2`. -/
theorem quarticCoreIntegral24_export :
    (∫ x : ℝ in 0..1,
      Real.log x * Real.log (1 - x) ^ 2 / x) =
      -(1 / 2 : ℝ) * (Real.pi ^ 4 / 90) := by
  exact quarticCoreIntegral24

/-- Integrability companion to `quarticCoreIntegral24_export`. -/
theorem quarticCoreIntervalIntegrable24_export :
    IntervalIntegrable
      (fun x : ℝ => Real.log x * Real.log (1 - x) ^ 2 / x)
      MeasureTheory.volume 0 1 := by
  exact quarticCoreKernel24_intervalIntegrable_half.trans
    quarticCoreKernel24_intervalIntegrable_upper

private def quarticAlternatingComplementMoment24
    (n : ℕ) (x : ℝ) : ℝ :=
  -signedHarmonic24 (n + 1) *
    (x ^ (n + 1) * Real.log x ^ 2)

private theorem quarticAlternatingComplementMoment24_intervalIntegrable
    (n : ℕ) :
    IntervalIntegrable (quarticAlternatingComplementMoment24 n)
      MeasureTheory.volume 0 1 := by
  unfold quarticAlternatingComplementMoment24
  apply ContinuousOn.intervalIntegrable
  rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
  exact
    (powerLogSquare24_continuousOn n).const_mul
      (-signedHarmonic24 (n + 1))

private theorem quarticAlternatingComplementMoment24_integral
    (n : ℕ) :
    (∫ x : ℝ in 0..1,
      quarticAlternatingComplementMoment24 n x) =
      2 * (-signedHarmonic24 (n + 1)) / (n + 2 : ℝ) ^ 3 := by
  unfold quarticAlternatingComplementMoment24
  rw [intervalIntegral.integral_const_mul,
    integral_powerLogSquare24]
  ring

private theorem quarticAlternatingComplementMoment24_integral_norm
    (n : ℕ) :
    (∫ x : ℝ in 0..1,
      ‖quarticAlternatingComplementMoment24 n x‖) =
      2 * |signedHarmonic24 (n + 1)| / (n + 2 : ℝ) ^ 3 := by
  have hpoint :
      ∀ x ∈ Set.uIcc (0 : ℝ) 1,
        ‖quarticAlternatingComplementMoment24 n x‖ =
          |signedHarmonic24 (n + 1)| *
            (x ^ (n + 1) * Real.log x ^ 2) := by
    intro x hx
    rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] at hx
    unfold quarticAlternatingComplementMoment24
    rw [Real.norm_eq_abs, abs_mul, abs_neg,
      abs_of_nonneg
        (mul_nonneg (pow_nonneg hx.1 _) (sq_nonneg _))]
  calc
    (∫ x : ℝ in 0..1,
        ‖quarticAlternatingComplementMoment24 n x‖) =
        ∫ x : ℝ in 0..1,
          |signedHarmonic24 (n + 1)| *
            (x ^ (n + 1) * Real.log x ^ 2) := by
      apply intervalIntegral.integral_congr
      exact hpoint
    _ = |signedHarmonic24 (n + 1)| *
          (∫ x : ℝ in 0..1,
            x ^ (n + 1) * Real.log x ^ 2) := by
      rw [intervalIntegral.integral_const_mul]
    _ = 2 * |signedHarmonic24 (n + 1)| /
          (n + 2 : ℝ) ^ 3 := by
      rw [integral_powerLogSquare24]
      ring

private theorem
    quarticAlternatingComplementMoment24_integral_norm_summable :
    Summable
      (fun n : ℕ =>
        ∫ x : ℝ in 0..1,
          ‖quarticAlternatingComplementMoment24 n x‖) := by
  apply quarticPlusSeriesMoment24_integral_norm_summable.of_nonneg_of_le
  · intro n
    rw [quarticAlternatingComplementMoment24_integral_norm]
    positivity
  · intro n
    rw [quarticAlternatingComplementMoment24_integral_norm,
      quarticPlusSeriesMoment24_integral_norm]
    exact div_le_div_of_nonneg_right
      (mul_le_mul_of_nonneg_left
        (abs_signedHarmonic24_le_harmonicNumber (n + 1))
        (by norm_num))
      (by positivity)

private theorem quarticAlternatingComplementMoment24_hasSum_pointwise
    {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasSum
      (fun n : ℕ => quarticAlternatingComplementMoment24 n x)
      (Real.log x ^ 2 * Real.log (1 + x) / (1 - x)) := by
  have habs : |x| < 1 := by rw [abs_of_pos hx0]; exact hx1
  have hxne : x ≠ 0 := ne_of_gt hx0
  have hH := harmonicNumber_generating_hasSum habs hxne
  have hP := parityRemainder24_generating_hasSum habs hxne
  have hscaled := (hH.sub hP).mul_left
    (x * Real.log x ^ 2 / 2)
  convert hscaled using 1
  · funext n
    unfold quarticAlternatingComplementMoment24 parityRemainder24
    rw [pow_succ]
    ring
  · have h1xne : 1 - x ≠ 0 := ne_of_gt (sub_pos.mpr hx1)
    field_simp [hxne, h1xne]
    ring

private theorem quarticAlternatingComplement_hasSum_integral :
    HasSum
      (fun n : ℕ =>
        2 * (-signedHarmonic24 (n + 1)) / (n + 2 : ℝ) ^ 3)
      (∫ x : ℝ in 0..1,
        Real.log x ^ 2 * Real.log (1 + x) / (1 - x)) := by
  have hInt : ∀ n : ℕ,
      MeasureTheory.Integrable
        (quarticAlternatingComplementMoment24 n)
        (MeasureTheory.volume.restrict (Set.Ioc (0 : ℝ) 1)) :=
    fun n => (quarticAlternatingComplementMoment24_intervalIntegrable n).1
  have hNorm : Summable (fun n : ℕ =>
      ∫ x : ℝ in Set.Ioc (0 : ℝ) 1,
        ‖quarticAlternatingComplementMoment24 n x‖) := by
    simpa only [
      ← intervalIntegral.integral_of_le
        (by norm_num : (0 : ℝ) ≤ 1)] using
      quarticAlternatingComplementMoment24_integral_norm_summable
  have h := MeasureTheory.hasSum_integral_of_summable_integral_norm
    (μ := MeasureTheory.volume.restrict (Set.Ioc (0 : ℝ) 1)) hInt hNorm
  have h' : HasSum
      (fun n : ℕ =>
        2 * (-signedHarmonic24 (n + 1)) / (n + 2 : ℝ) ^ 3)
      (∫ x : ℝ in Set.Ioc (0 : ℝ) 1,
        ∑' n : ℕ, quarticAlternatingComplementMoment24 n x) := by
    convert h using 1
    funext n
    rw [← intervalIntegral.integral_of_le
      (by norm_num : (0 : ℝ) ≤ 1)]
    exact (quarticAlternatingComplementMoment24_integral n).symm
  convert h' using 1
  rw [intervalIntegral.integral_of_le (by norm_num : (0 : ℝ) ≤ 1)]
  apply MeasureTheory.setIntegral_congr_ae measurableSet_Ioc
  filter_upwards [
    MeasureTheory.Measure.ae_ne MeasureTheory.volume (1 : ℝ)
  ] with x hxne hx
  exact
    (quarticAlternatingComplementMoment24_hasSum_pointwise
      hx.1 (lt_of_le_of_ne hx.2 hxne)).tsum_eq.symm

/-- Integrability of the alternating quartic complement. -/
theorem quarticAlternatingComplementIntervalIntegrable24 :
    IntervalIntegrable
      (fun x : ℝ =>
        Real.log x ^ 2 * Real.log (1 + x) / (1 - x))
      MeasureTheory.volume 0 1 := by
  apply IntervalIntegrable.trans (b := (1 / 2 : ℝ))
  · have hcont : ContinuousOn
        (fun x : ℝ =>
          (x * Real.log x ^ 2) *
            RamanujanChallenge.P26.logOnePlusSlope26 x / (1 - x))
        (Icc (0 : ℝ) (1 / 2)) := by
      apply ContinuousOn.div
        (selfMulLogSquare24_continuousOn.mono
          (fun x hx => ⟨hx.1, hx.2.trans (by norm_num)⟩) |>.mul
          (RamanujanChallenge.P26.logOnePlusSlope26_continuousOn.mono
            (fun x hx => ⟨hx.1, hx.2.trans (by norm_num)⟩)))
        (continuousOn_const.sub continuousOn_id)
      intro x hx
      exact ne_of_gt (show 0 < 1 - x by linarith [hx.2])
    have hint : IntervalIntegrable
        (fun x : ℝ =>
          (x * Real.log x ^ 2) *
            RamanujanChallenge.P26.logOnePlusSlope26 x / (1 - x))
        MeasureTheory.volume 0 (1 / 2) := by
      apply ContinuousOn.intervalIntegrable
      rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2)]
      exact hcont
    apply IntervalIntegrable.congr
      (f := fun x : ℝ =>
        (x * Real.log x ^ 2) *
          RamanujanChallenge.P26.logOnePlusSlope26 x / (1 - x))
      ?_ hint
    intro x _
    by_cases hxzero : x = 0
    · subst x
      simp
    · simp [RamanujanChallenge.P26.logOnePlusSlope26, hxzero]
      field_simp [hxzero]
  · have hsub : ContinuousOn (fun x : ℝ => 1 - x)
        (Icc (1 / 2 : ℝ) 1) := by fun_prop
    have hmap : MapsTo (fun x : ℝ => 1 - x)
        (Icc (1 / 2 : ℝ) 1) (Icc (0 : ℝ) (1 / 2)) := by
      intro x hx
      constructor <;> linarith [hx.1, hx.2]
    have hslope := logOneMinusSlope24_continuousOn.comp hsub hmap
    have hlog : ContinuousOn (fun x : ℝ => Real.log x)
        (Icc (1 / 2 : ℝ) 1) := by
      apply continuousOn_id.log
      intro x hx
      exact ne_of_gt (show 0 < x by linarith [hx.1])
    have hplus : ContinuousOn (fun x : ℝ => Real.log (1 + x))
        (Icc (1 / 2 : ℝ) 1) := by
      apply (continuousOn_const.add continuousOn_id).log
      intro x hx
      exact ne_of_gt (show 0 < 1 + x by linarith [hx.1])
    have hcont : ContinuousOn
        (fun x : ℝ =>
          logOneMinusSlope24 (1 - x) * Real.log x *
            Real.log (1 + x))
        (Icc (1 / 2 : ℝ) 1) := (hslope.mul hlog).mul hplus
    have hint : IntervalIntegrable
        (fun x : ℝ =>
          logOneMinusSlope24 (1 - x) * Real.log x *
            Real.log (1 + x))
        MeasureTheory.volume (1 / 2) 1 := by
      apply ContinuousOn.intervalIntegrable
      rw [Set.uIcc_of_le (by norm_num : (1 / 2 : ℝ) ≤ 1)]
      exact hcont
    apply IntervalIntegrable.congr
      (f := fun x : ℝ =>
        logOneMinusSlope24 (1 - x) * Real.log x *
          Real.log (1 + x)) ?_ hint
    intro x _
    by_cases hxone : x = 1
    · subst x
      simp
    · have h1xne : 1 - x ≠ 0 := sub_ne_zero.mpr (Ne.symm hxone)
      simp [logOneMinusSlope24, h1xne]
      field_simp [h1xne]

/-- `∫₀¹ log²x · log(1+x)/(1-x)` in closed form. -/
theorem quarticAlternatingComplementIntegral24 :
    (∫ x : ℝ in 0..1,
      Real.log x ^ 2 * Real.log (1 + x) / (1 - x)) =
      (7 / 2 : ℝ) * Real.log 2 * zeta3_24 -
        (19 / 8 : ℝ) * (Real.pi ^ 4 / 90) := by
  have hBar : HasSum barHarmonicCubicTerm24
      ((7 / 4 : ℝ) * Real.log 2 * zeta3_24 -
        (1 / 8 : ℝ) * (Real.pi ^ 2 / 6) ^ 2) := by
    rw [← barHarmonicCubicIntegral24]
    exact barHarmonicCubicTerm24_hasSum_integral
  have hBarTail : HasSum
      (fun n : ℕ => barHarmonicCubicTerm24 (n + 1))
      ((7 / 4 : ℝ) * Real.log 2 * zeta3_24 -
        (1 / 8 : ℝ) * (Real.pi ^ 2 / 6) ^ 2 - 1) := by
    apply (hasSum_nat_add_iff (f := barHarmonicCubicTerm24) 1).mpr
    convert hBar using 1 <;>
      simp [barHarmonicCubicTerm24, signedHarmonic24] <;> ring
  have hZetaTail : HasSum
      (fun n : ℕ => alternatingZetaFourTerm24 (n + 1))
      (1 - (7 / 8 : ℝ) * (Real.pi ^ 4 / 90)) := by
    apply
      (hasSum_nat_add_iff
        (f := alternatingZetaFourTerm24) 1).mpr
    convert alternatingZetaFourTerm24_hasSum using 1 <;>
      simp [alternatingZetaFourTerm24, zetaFourTerm24] <;> ring
  have hValue := (hBarTail.add hZetaTail).mul_left 2
  have hSeries : HasSum
      (fun n : ℕ =>
        2 * (-signedHarmonic24 (n + 1)) / (n + 2 : ℝ) ^ 3)
      ((7 / 2 : ℝ) * Real.log 2 * zeta3_24 -
        (19 / 8 : ℝ) * (Real.pi ^ 4 / 90)) := by
    convert hValue using 1
    · funext n
      unfold barHarmonicCubicTerm24 alternatingZetaFourTerm24
        zetaFourTerm24
      rw [signedHarmonic24_succ (n + 1)]
      push_cast
      have hden : (n + 2 : ℝ) ≠ 0 := by positivity
      field_simp [hden]
      simp only [pow_succ]
      ring
    · ring
  exact quarticAlternatingComplement_hasSum_integral.unique hSeries

end
