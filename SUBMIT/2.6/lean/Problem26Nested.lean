/-
  Analytic reduction for the remaining nested inverse-binomial sum in
  Ramanujan Challenge Problem 2.6.

  This file is downstream from `Problem26`: the recurrence, summability, and
  weight-two inverse-binomial evaluation remain in that file.  Here we isolate
  the real cyclotomic integral that carries the weight-three value.
-/
import RamanujanChallenge.Problem26
import Mathlib.NumberTheory.ZetaValues
import Mathlib.MeasureTheory.Integral.DominatedConvergence
import Mathlib.Analysis.Calculus.Deriv.Slope

open Filter Set Topology
open scoped Interval

noncomputable section

namespace RamanujanChallenge.P26

/-! ## Elementary logarithmic moments -/

/-- The logarithmic moment used to integrate power-series expansions at their
radius-one endpoint. -/
theorem integral_pow_mul_log26 (n : ℕ) :
    (∫ x : ℝ in 0..1, x ^ n * Real.log x) =
      -1 / ((n : ℝ) + 1) ^ 2 := by
  let F : ℝ → ℝ := fun x =>
    x ^ (n + 1) *
      ((((n : ℝ) + 1) * Real.log x - 1) / ((n : ℝ) + 1) ^ 2)
  have hInt :
      IntervalIntegrable (fun x : ℝ => x ^ n * Real.log x)
        MeasureTheory.volume 0 1 := by
    exact intervalIntegral.intervalIntegrable_log'.continuousOn_mul
      (continuousOn_pow n)
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_tendsto
    (f := F) (fa := 0) (fb := -1 / ((n : ℝ) + 1) ^ 2) (hint := hInt)]
  · ring
  · norm_num
  · intro x hx
    have hx0 : x ≠ 0 := ne_of_gt hx.1
    have hpow :
        HasDerivAt (fun y : ℝ => y ^ (n + 1))
          (((n : ℝ) + 1) * x ^ n) x := by
      convert hasDerivAt_pow (n + 1) x using 1 <;>
        norm_num [Nat.cast_add]
    have hlog :
        HasDerivAt
          (fun y : ℝ =>
            (((n : ℝ) + 1) * Real.log y - 1) /
              ((n : ℝ) + 1) ^ 2)
          ((((n : ℝ) + 1) * x⁻¹) / ((n : ℝ) + 1) ^ 2) x :=
      (((Real.hasDerivAt_log hx0).const_mul ((n : ℝ) + 1)).sub_const 1).div_const
        (((n : ℝ) + 1) ^ 2)
    dsimp [F]
    convert hpow.mul hlog using 1
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
      · ring
    simpa [F] using hmain
  · have hcont : ContinuousAt F 1 := by
      dsimp [F]
      exact (continuousAt_id.pow (n + 1)).mul
        ((((continuousAt_const.mul (Real.continuousAt_log (by norm_num))).sub
          continuousAt_const).div_const (((n : ℝ) + 1) ^ 2)))
    convert tendsto_nhdsWithin_of_tendsto_nhds hcont.tendsto using 1
    dsimp [F]
    norm_num

/-! ## A real trilogarithm on the open unit disk -/

/-- The cubic polylogarithm, needed only on `[-1, 1]` below. -/
def trilog26 (z : ℝ) : ℝ :=
  ∑' n : ℕ, z ^ (n + 1) / ((n + 1 : ℕ) : ℝ) ^ 3

@[simp] theorem trilog26_zero : trilog26 0 = 0 := by
  simp [trilog26, zero_pow (by omega : _ + 1 ≠ 0)]

theorem trilog26_one :
    trilog26 1 = zeta3 := by
  unfold trilog26 zeta3
  apply tsum_congr
  intro n
  simp [Nat.cast_add]

theorem trilog26_neg_one :
    trilog26 (-1) = -(3 : ℝ) / 4 * zeta3 := by
  unfold trilog26
  calc
    (∑' n : ℕ, (-1 : ℝ) ^ (n + 1) / ((n + 1 : ℕ) : ℝ) ^ 3) =
        ∑' n : ℕ, alternatingZeta3Term26 n := by
          apply tsum_congr
          intro n
          simp [alternatingZeta3Term26, zeta3Term26, div_eq_mul_inv]
    _ = -(3 : ℝ) / 4 * zeta3 :=
      alternatingZeta3Term26_hasSum.tsum_eq

theorem trilog26_hasDerivAt_of_abs_lt_one
    {x : ℝ} (hx : |x| < 1) (hxne : x ≠ 0) :
    HasDerivAt trilog26 (dilog x / x) x := by
  let r : ℝ := (|x| + 1) / 2
  have hr0 : 0 ≤ r := by dsimp [r]; positivity
  have hrpos : 0 < r := by dsimp [r]; positivity
  have hr1 : r < 1 := by dsimp [r]; linarith
  have hxr : |x| < r := by dsimp [r]; linarith
  have hu : Summable (fun n : ℕ => r ^ n) :=
    summable_geometric_of_lt_one hr0 hr1
  have hterm : ∀ n : ℕ, ∀ y : ℝ, y ∈ Ioo (-r) r →
      HasDerivAt
        (fun y : ℝ => y ^ (n + 1) / ((n + 1 : ℕ) : ℝ) ^ 3)
        (y ^ n / ((n + 1 : ℕ) : ℝ) ^ 2) y := by
    intro n y _
    convert
      (hasDerivAt_pow (n + 1) y).div_const
        (((n + 1 : ℕ) : ℝ) ^ 3) using 1
    rw [Nat.add_sub_cancel]
    simp only [Nat.cast_add, Nat.cast_one]
    field_simp
  have hbound : ∀ n : ℕ, ∀ y : ℝ, y ∈ Ioo (-r) r →
      ‖y ^ n / ((n + 1 : ℕ) : ℝ) ^ 2‖ ≤ r ^ n := by
    intro n y hy
    rw [Real.norm_eq_abs, abs_div, abs_pow]
    have hyr : |y| < r := (abs_lt).2 hy
    have hden : (1 : ℝ) ≤ |((n + 1 : ℕ) : ℝ) ^ 2| := by
      rw [abs_of_nonneg (sq_nonneg _)]
      have hn : (1 : ℝ) ≤ ((n + 1 : ℕ) : ℝ) := by norm_num
      nlinarith
    calc
      |y| ^ n / |((n + 1 : ℕ) : ℝ) ^ 2| ≤ |y| ^ n / 1 := by
        gcongr
      _ = |y| ^ n := by ring
      _ ≤ r ^ n := pow_le_pow_left₀ (abs_nonneg y) hyr.le n
  have hzero : Summable (fun n : ℕ =>
      (0 : ℝ) ^ (n + 1) / ((n + 1 : ℕ) : ℝ) ^ 3) := by
    simp [zero_pow (Nat.succ_ne_zero _)]
  have hxmem : x ∈ Ioo (-r) r := (abs_lt).1 hxr
  have hd :
      HasDerivAt
        (fun y : ℝ => ∑' n : ℕ,
          y ^ (n + 1) / ((n + 1 : ℕ) : ℝ) ^ 3)
        (∑' n : ℕ, x ^ n / ((n + 1 : ℕ) : ℝ) ^ 2) x := by
    exact hasDerivAt_tsum_of_isPreconnected hu isOpen_Ioo
      (convex_Ioo (-r) r).isPreconnected hterm hbound
      (show (0 : ℝ) ∈ Ioo (-r) r by constructor <;> linarith)
      hzero hxmem
  rw [show trilog26 = (fun y : ℝ => ∑' n : ℕ,
      y ^ (n + 1) / ((n + 1 : ℕ) : ℝ) ^ 3) from rfl]
  convert hd using 1
  have hsumm := dilog_summable (le_of_lt hx)
  have hsum :=
    hsumm.hasSum.div_const x
  have hsum' :
      HasSum
        (fun n : ℕ => x ^ n / ((n + 1 : ℕ) : ℝ) ^ 2)
        (dilog x / x) := by
    convert hsum using 1
    funext n
    rw [pow_succ]
    field_simp [hxne]
    ring
  exact hsum'.tsum_eq.symm

theorem dilog_hasDerivAt_zero26 :
    HasDerivAt dilog 1 0 := by
  let r : ℝ := 1 / 2
  have hr0 : 0 ≤ r := by dsimp [r]; norm_num
  have hr1 : r < 1 := by dsimp [r]; norm_num
  have hu : Summable (fun n : ℕ => r ^ n) :=
    summable_geometric_of_lt_one hr0 hr1
  have hterm : ∀ n : ℕ, ∀ y : ℝ, y ∈ Ioo (-r) r →
      HasDerivAt
        (fun y : ℝ => y ^ (n + 1) / ((n + 1 : ℕ) : ℝ) ^ 2)
        (y ^ n / ((n + 1 : ℕ) : ℝ)) y := by
    intro n y _
    convert
      (hasDerivAt_pow (n + 1) y).div_const
        (((n + 1 : ℕ) : ℝ) ^ 2) using 1
    rw [Nat.add_sub_cancel]
    simp only [Nat.cast_add, Nat.cast_one]
    field_simp
  have hbound : ∀ n : ℕ, ∀ y : ℝ, y ∈ Ioo (-r) r →
      ‖y ^ n / ((n + 1 : ℕ) : ℝ)‖ ≤ r ^ n := by
    intro n y hy
    rw [Real.norm_eq_abs, abs_div, abs_pow]
    have hyr : |y| < r := (abs_lt).2 hy
    calc
      |y| ^ n / |((n + 1 : ℕ) : ℝ)| ≤ |y| ^ n / 1 := by
        gcongr
        rw [abs_of_nonneg (by positivity)]
        norm_num
      _ = |y| ^ n := by ring
      _ ≤ r ^ n := pow_le_pow_left₀ (abs_nonneg y) hyr.le n
  have hzero : Summable (fun n : ℕ =>
      (0 : ℝ) ^ (n + 1) / ((n + 1 : ℕ) : ℝ) ^ 2) := by
    simp [zero_pow (Nat.succ_ne_zero _)]
  have hd :
      HasDerivAt
        (fun y : ℝ => ∑' n : ℕ,
          y ^ (n + 1) / ((n + 1 : ℕ) : ℝ) ^ 2)
        (∑' n : ℕ, (0 : ℝ) ^ n / ((n + 1 : ℕ) : ℝ)) 0 := by
    exact hasDerivAt_tsum_of_isPreconnected hu isOpen_Ioo
      (convex_Ioo (-r) r).isPreconnected hterm hbound
      (show (0 : ℝ) ∈ Ioo (-r) r by constructor <;> norm_num [r])
      hzero (show (0 : ℝ) ∈ Ioo (-r) r by constructor <;> norm_num [r])
  rw [show dilog = (fun y : ℝ => ∑' n : ℕ,
      y ^ (n + 1) / ((n + 1 : ℕ) : ℝ) ^ 2) from rfl]
  convert hd using 1
  symm
  calc
    (∑' n : ℕ, (0 : ℝ) ^ n / ((n + 1 : ℕ) : ℝ)) =
        ∑' n : ℕ, if n = 0 then (1 : ℝ) else 0 := by
          apply tsum_congr
          intro n
          rcases n with _ | n
          · norm_num
          · simp [zero_pow (Nat.succ_ne_zero n)]
    _ = 1 := tsum_ite_eq 0 1

theorem trilog26_continuousOn_unit :
    ContinuousOn trilog26 (Icc (-1 : ℝ) 1) := by
  unfold trilog26
  have hu :
      Summable (fun n : ℕ =>
        (1 : ℝ) / ((n + 1 : ℕ) : ℝ) ^ 3) := by
    convert zeta3Term26_summable using 1
    funext n
    simp only [zeta3Term26, Nat.cast_add, Nat.cast_one]
  refine continuousOn_tsum
    (u := fun n : ℕ => (1 : ℝ) / ((n + 1 : ℕ) : ℝ) ^ 3)
    (fun n => ?_) hu ?_
  · fun_prop
  · intro n x hx
    rw [Real.norm_eq_abs, abs_div, abs_pow]
    have habs : |x| ≤ 1 := (abs_le).2 hx
    have hpos : (0 : ℝ) ≤ ((n + 1 : ℕ) : ℝ) ^ 3 := by positivity
    rw [abs_of_nonneg hpos]
    apply div_le_div_of_nonneg_right _ hpos
    calc
      |x| ^ (n + 1) ≤ 1 ^ (n + 1) :=
        pow_le_pow_left₀ (abs_nonneg x) habs (n + 1)
      _ = 1 := one_pow _

/-! ## The cyclotomic integral and its rational decomposition -/

def nestedCyclotomicKernel26 (r : ℝ) : ℝ :=
  (1 - r) * ((1 - r) ^ 2 * Real.log (1 + r) - r) * Real.log r /
    (r * (1 + r) * (1 + r + r ^ 2))

def cyclotomicWeightThreeKernel26 (r : ℝ) : ℝ :=
  (2 * r + 1) / (1 + r + r ^ 2) * Real.log r * Real.log (1 + r)

def alternatingWeightThreeKernel26 (r : ℝ) : ℝ :=
  Real.log r * Real.log (1 + r) / (1 + r)

def radialWeightThreeKernel26 (r : ℝ) : ℝ :=
  Real.log r * Real.log (1 + r) / r

def radialWeightThreePrimitive26 (r : ℝ) : ℝ :=
  -(Real.log r * dilog (-r)) + trilog26 (-r)

def elementaryWeightTwoKernel26 (r : ℝ) : ℝ :=
  ((2 * r + 1) / (1 + r + r ^ 2) - 2 / (1 + r)) * Real.log r

/-- The continuous extension of `log (1+r) / r` at `r=0`. -/
def logOnePlusSlope26 (r : ℝ) : ℝ :=
  Function.update (fun x : ℝ => Real.log (1 + x) / x) 0 1 r

/-- A nonsingular logarithmic factor whose derivative is the rational
coefficient in `elementaryWeightTwoKernel26`. -/
def elementaryLogFactor26 (r : ℝ) : ℝ :=
  Real.log (1 + r + r ^ 2) - 2 * Real.log (1 + r)

/-- The dilogarithmic correction that cancels the `h(r)/r` term created by
differentiating `log r * h(r)`. -/
def elementaryDilogCorrection26 (r : ℝ) : ℝ :=
  (1 / 3 : ℝ) * dilog (r ^ 3) - dilog r - 2 * dilog (-r)

def elementaryWeightTwoPrimitive26 (r : ℝ) : ℝ :=
  Real.log r * elementaryLogFactor26 r + elementaryDilogCorrection26 r

theorem cyclotomicQuadratic_pos26 (r : ℝ) :
    0 < 1 + r + r ^ 2 := by
  nlinarith [sq_nonneg (r + 1 / 2)]

theorem radialWeightThreePrimitive_hasDerivAt26
    {r : ℝ} (hr0 : 0 < r) (hr1 : r < 1) :
    HasDerivAt radialWeightThreePrimitive26
      (radialWeightThreeKernel26 r) r := by
  have hrne : r ≠ 0 := ne_of_gt hr0
  have hrabs : |r| < 1 := by rw [abs_of_pos hr0]; exact hr1
  have hnegabs : |-r| < 1 := by simpa using hrabs
  have hDnegRaw :=
    (dilog_hasDerivAt_of_abs_lt_one hnegabs (neg_ne_zero.mpr hrne)).comp r
      (hasDerivAt_neg r)
  have hDneg :
      HasDerivAt (fun x : ℝ => dilog (-x))
        (-Real.log (1 + r) / r) r := by
    convert hDnegRaw using 1
    field_simp [hrne]
    ring
  have hTnegRaw :=
    (trilog26_hasDerivAt_of_abs_lt_one hnegabs (neg_ne_zero.mpr hrne)).comp r
      (hasDerivAt_neg r)
  have hTneg :
      HasDerivAt (fun x : ℝ => trilog26 (-x))
        (dilog (-r) / r) r := by
    convert hTnegRaw using 1
    field_simp [hrne]
  have hprod := (Real.hasDerivAt_log hrne).mul hDneg
  convert hprod.neg.add hTneg using 1
  simp only [radialWeightThreeKernel26]
  field_simp [hrne]
  ring

theorem logOnePlusSlope26_continuousOn :
    ContinuousOn logOnePlusSlope26 (Icc (0 : ℝ) 1) := by
  intro r hr
  by_cases hrzero : r = 0
  · subst r
    have hlog :
        HasDerivAt (fun x : ℝ => Real.log (1 + x)) 1 0 := by
      have hone : HasDerivAt (fun x : ℝ => 1 + x) 1 0 := by
        convert (hasDerivAt_const (0 : ℝ) (1 : ℝ)).add
          (hasDerivAt_id (0 : ℝ)) using 1 <;> ring
      simpa [Function.comp_def] using
        (HasDerivAt.comp (h := fun x : ℝ => 1 + x) 0
          (Real.hasDerivAt_log
            (by norm_num : (fun x : ℝ => 1 + x) 0 ≠ 0)) hone)
    have hc := hlog.continuousAt_div
    have hc' : ContinuousAt logOnePlusSlope26 0 := by
      convert hc using 1
      funext x
      simp [logOnePlusSlope26]
    exact hc'.continuousWithinAt
  · have hbase :
        ContinuousAt (fun x : ℝ => Real.log (1 + x) / x) r := by
      have hone : ContinuousAt (fun x : ℝ => 1 + x) r := by fun_prop
      have hlog :
          ContinuousAt (fun x : ℝ => Real.log (1 + x)) r :=
        (Real.continuousAt_log (by linarith [hr.1] : 1 + r ≠ 0)).comp hone
      exact hlog.div continuousAt_id hrzero
    have heq :
        logOnePlusSlope26 =ᶠ[𝓝 r]
          (fun x : ℝ => Real.log (1 + x) / x) := by
      filter_upwards [isOpen_ne.mem_nhds hrzero] with x hx
      simp [logOnePlusSlope26, hx]
    exact (hbase.congr_of_eventuallyEq heq).continuousWithinAt

theorem radialWeightThreeKernel_intervalIntegrable26 :
    IntervalIntegrable radialWeightThreeKernel26
      MeasureTheory.volume 0 1 := by
  have hlog :
      IntervalIntegrable Real.log MeasureTheory.volume (0 : ℝ) 1 :=
    intervalIntegral.intervalIntegrable_log'
  have hslope :
      ContinuousOn logOnePlusSlope26 (Set.uIcc (0 : ℝ) 1) := by
    simpa [Set.uIcc_of_le (show (0 : ℝ) ≤ 1 by norm_num)] using
      logOnePlusSlope26_continuousOn
  have hint := hlog.continuousOn_mul hslope
  apply IntervalIntegrable.congr
    (f := fun r : ℝ => logOnePlusSlope26 r * Real.log r) ?_ hint
  intro r hr
  have hrne : r ≠ 0 := by
    have hr' : r ∈ Ioc (0 : ℝ) 1 := by
      simpa [Set.uIoc_of_le (show (0 : ℝ) ≤ 1 by norm_num)] using hr
    exact ne_of_gt hr'.1
  simp [radialWeightThreeKernel26, logOnePlusSlope26, hrne]
  ring

theorem radialWeightThreePrimitive_tendsto_zero26 :
    Tendsto radialWeightThreePrimitive26 (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  have hDzero :
      HasDerivAt dilog 1 ((fun r : ℝ => -r) 0) := by
    simpa using dilog_hasDerivAt_zero26
  have hDnegRaw :=
    HasDerivAt.comp (h := fun r : ℝ => -r) 0 hDzero
      (hasDerivAt_neg (0 : ℝ))
  have hDneg :
      HasDerivAt (fun r : ℝ => dilog (-r)) (-1) 0 := by
    simpa [Function.comp_def] using hDnegRaw
  have hslope :
      Tendsto (fun r : ℝ => r⁻¹ * dilog (-r))
        (𝓝[>] 0) (𝓝 (-1)) := by
    simpa [dilog_zero] using hDneg.tendsto_slope_zero_right
  have hlogr :
      Tendsto (fun r : ℝ => Real.log r * r) (𝓝[>] 0) (𝓝 0) := by
    simpa [Real.rpow_one] using
      (tendsto_log_mul_rpow_nhdsGT_zero (show (0 : ℝ) < 1 by norm_num))
  have hproductRaw :
      Tendsto
        (fun r : ℝ => (Real.log r * r) * (r⁻¹ * dilog (-r)))
        (𝓝[>] 0) (𝓝 0) := by
    simpa using hlogr.mul hslope
  have hproduct :
      Tendsto (fun r : ℝ => Real.log r * dilog (-r))
        (𝓝[>] 0) (𝓝 0) := by
    apply Filter.Tendsto.congr' _ hproductRaw
    filter_upwards [self_mem_nhdsWithin] with r hr
    have hrne : r ≠ 0 := ne_of_gt hr
    field_simp [hrne]
  have htriCont : ContinuousAt (fun r : ℝ => trilog26 (-r)) 0 := by
    have hbase : ContinuousAt trilog26 0 :=
      trilog26_continuousOn_unit.continuousAt
        (Icc_mem_nhds (by norm_num : (-1 : ℝ) < 0)
          (by norm_num : (0 : ℝ) < 1))
    have hbase' :
        ContinuousAt trilog26 ((fun r : ℝ => -r) 0) := by
      simpa using hbase
    have hneg : ContinuousAt (fun r : ℝ => -r) 0 :=
      continuousAt_id.neg
    simpa [Function.comp_def] using hbase'.comp hneg
  have htri :
      Tendsto (fun r : ℝ => trilog26 (-r)) (𝓝[>] 0) (𝓝 0) := by
    simpa using tendsto_nhdsWithin_of_tendsto_nhds htriCont.tendsto
  simpa [radialWeightThreePrimitive26] using hproduct.neg.add htri

theorem radialWeightThreePrimitive_tendsto_one26 :
    Tendsto radialWeightThreePrimitive26 (𝓝[<] (1 : ℝ))
      (𝓝 (-(3 : ℝ) / 4 * zeta3)) := by
  have hlog :
      Tendsto Real.log (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    simpa using tendsto_nhdsWithin_of_tendsto_nhds
      (Real.continuousAt_log (by norm_num : (1 : ℝ) ≠ 0)).tendsto
  have hDcont :
      ContinuousOn (fun r : ℝ => dilog (-r)) (Icc (0 : ℝ) 1) := by
    apply dilog_continuousOn_unit.comp
    · fun_prop
    · intro r hr
      constructor <;> linarith [hr.1, hr.2]
  have hDtend :
      Tendsto (fun r : ℝ => dilog (-r)) (𝓝[<] 1) (𝓝 (dilog (-1))) := by
    have hwithin :
        ContinuousWithinAt (fun r : ℝ => dilog (-r)) (Iio (1 : ℝ)) 1 :=
      (hDcont 1 (by norm_num)).mono_of_mem_nhdsWithin
        (Icc_mem_nhdsLT (show (0 : ℝ) < 1 by norm_num))
    exact hwithin.tendsto
  have hproduct :
      Tendsto (fun r : ℝ => Real.log r * dilog (-r))
        (𝓝[<] 1) (𝓝 0) := by
    simpa using hlog.mul hDtend
  have hTcont :
      ContinuousOn (fun r : ℝ => trilog26 (-r)) (Icc (0 : ℝ) 1) := by
    apply trilog26_continuousOn_unit.comp
    · fun_prop
    · intro r hr
      constructor <;> linarith [hr.1, hr.2]
  have hTtend :
      Tendsto (fun r : ℝ => trilog26 (-r)) (𝓝[<] 1)
        (𝓝 (-(3 : ℝ) / 4 * zeta3)) := by
    have hwithin :
        ContinuousWithinAt (fun r : ℝ => trilog26 (-r)) (Iio (1 : ℝ)) 1 :=
      (hTcont 1 (by norm_num)).mono_of_mem_nhdsWithin
        (Icc_mem_nhdsLT (show (0 : ℝ) < 1 by norm_num))
    simpa [trilog26_neg_one] using hwithin.tendsto
  simpa [radialWeightThreePrimitive26] using hproduct.neg.add hTtend

theorem radialWeightThreeIntegral26 :
    (∫ r : ℝ in 0..1, radialWeightThreeKernel26 r) =
      -(3 : ℝ) / 4 * zeta3 := by
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_tendsto
    (f := radialWeightThreePrimitive26)
    (fa := (0 : ℝ)) (fb := -(3 : ℝ) / 4 * zeta3)
    (by norm_num)
    (fun r hr => radialWeightThreePrimitive_hasDerivAt26 hr.1 hr.2)
    radialWeightThreeKernel_intervalIntegrable26
    radialWeightThreePrimitive_tendsto_zero26
    radialWeightThreePrimitive_tendsto_one26]
  ring

theorem cubicRationalDecomposition26
    (r : ℝ) (hr0 : r ≠ 0) (hr1 : 1 + r ≠ 0) :
    (1 - r) ^ 3 / (r * (1 + r) * (1 + r + r ^ 2)) =
      3 * (2 * r + 1) / (1 + r + r ^ 2) -
        8 / (1 + r) + 1 / r := by
  have hq : 1 + r + r ^ 2 ≠ 0 := (cyclotomicQuadratic_pos26 r).ne'
  field_simp [hr0, hr1, hq]
  ring

theorem linearRationalDecomposition26 (r : ℝ) (hr1 : 1 + r ≠ 0) :
    -(1 - r) / ((1 + r) * (1 + r + r ^ 2)) =
      (2 * r + 1) / (1 + r + r ^ 2) - 2 / (1 + r) := by
  have hq : 1 + r + r ^ 2 ≠ 0 := (cyclotomicQuadratic_pos26 r).ne'
  field_simp [hr1, hq]
  ring

theorem elementaryLogFactor_zero26 : elementaryLogFactor26 0 = 0 := by
  simp [elementaryLogFactor26]

theorem elementaryDilogCorrection_zero26 :
    elementaryDilogCorrection26 0 = 0 := by
  simp [elementaryDilogCorrection26, dilog_zero]

theorem elementaryDilogCorrection_one26 :
    elementaryDilogCorrection26 1 = Real.pi ^ 2 / 18 := by
  simp [elementaryDilogCorrection26, dilog_one, dilog_neg_one26]
  ring

theorem elementaryDilogCorrection_continuousOn26 :
    ContinuousOn elementaryDilogCorrection26 (Icc (0 : ℝ) 1) := by
  have hcube :
      ContinuousOn (fun r : ℝ => dilog (r ^ 3)) (Icc (0 : ℝ) 1) := by
    apply dilog_continuousOn_unit.comp
    · fun_prop
    · intro r hr
      constructor
      · have : 0 ≤ r ^ 3 := pow_nonneg hr.1 3
        linarith
      · exact pow_le_one₀ hr.1 hr.2
  have hid : ContinuousOn dilog (Icc (0 : ℝ) 1) :=
    dilog_continuousOn_unit.mono (by
      intro r hr
      constructor <;> linarith [hr.1, hr.2])
  have hneg :
      ContinuousOn (fun r : ℝ => dilog (-r)) (Icc (0 : ℝ) 1) := by
    apply dilog_continuousOn_unit.comp
    · fun_prop
    · intro r hr
      constructor <;> linarith [hr.1, hr.2]
  unfold elementaryDilogCorrection26
  fun_prop

theorem elementaryLogFactor_hasDerivAt26
    (r : ℝ) (hr1 : 1 + r ≠ 0) :
    HasDerivAt elementaryLogFactor26
      ((2 * r + 1) / (1 + r + r ^ 2) - 2 / (1 + r)) r := by
  have hq : 1 + r + r ^ 2 ≠ 0 := (cyclotomicQuadratic_pos26 r).ne'
  have hquad :
      HasDerivAt (fun x : ℝ => 1 + x + x ^ 2) (2 * r + 1) r := by
    convert
      ((hasDerivAt_const r (1 : ℝ)).add (hasDerivAt_id r) |>.add
        (hasDerivAt_pow 2 r)) using 1 <;> ring
  have hone : HasDerivAt (fun x : ℝ => 1 + x) 1 r := by
    convert (hasDerivAt_const r (1 : ℝ)).add (hasDerivAt_id r) using 1 <;> ring
  have hlogq :=
    (Real.hasDerivAt_log hq).comp r hquad
  have hlogone :=
    (Real.hasDerivAt_log hr1).comp r hone
  convert hlogq.sub (hlogone.const_mul 2) using 1
  field_simp [hq, hr1]

/-- On the open unit interval the nonsingular factor is the cyclotomic
factorization needed by the dilogarithm derivatives. -/
theorem elementaryLogFactor_eq_cyclotomic26
    {r : ℝ} (hr0 : 0 < r) (hr1 : r < 1) :
    elementaryLogFactor26 r =
      Real.log (1 - r ^ 3) - Real.log (1 - r) -
        2 * Real.log (1 + r) := by
  have hm : 1 - r ≠ 0 := by linarith
  have hq : 1 + r + r ^ 2 ≠ 0 := (cyclotomicQuadratic_pos26 r).ne'
  rw [show 1 - r ^ 3 = (1 - r) * (1 + r + r ^ 2) by ring,
    Real.log_mul hm hq]
  simp only [elementaryLogFactor26]
  ring

theorem elementaryDilogCorrection_hasDerivAt26
    {r : ℝ} (hr0 : 0 < r) (hr1 : r < 1) :
    HasDerivAt elementaryDilogCorrection26
      (-elementaryLogFactor26 r / r) r := by
  have hrne : r ≠ 0 := ne_of_gt hr0
  have hrabs : |r| < 1 := by rw [abs_of_pos hr0]; exact hr1
  have hr3pos : 0 < r ^ 3 := by positivity
  have hr3lt : r ^ 3 < 1 := by
    nlinarith [mul_self_lt_mul_self (show 0 ≤ r by positivity) hr1]
  have hr3abs : |r ^ 3| < 1 := by
    rw [abs_of_pos hr3pos]
    exact hr3lt
  have hcube : HasDerivAt (fun x : ℝ => x ^ 3) (3 * r ^ 2) r := by
    convert hasDerivAt_pow 3 r using 1 <;> ring
  have hDcube :=
    HasDerivAt.comp (h := fun x : ℝ => x ^ 3) r
      (dilog_hasDerivAt_of_abs_lt_one hr3abs hr3pos.ne') hcube
  have hDr := dilog_hasDerivAt_of_abs_lt_one hrabs hrne
  have hnegabs : |-r| < 1 := by simpa using hrabs
  have hDneg :=
    (dilog_hasDerivAt_of_abs_lt_one hnegabs (neg_ne_zero.mpr hrne)).comp r
      (hasDerivAt_neg r)
  have h :=
    (hDcube.const_mul (1 / 3 : ℝ)).sub hDr |>.sub (hDneg.const_mul 2)
  convert h using 1
  rw [elementaryLogFactor_eq_cyclotomic26 hr0 hr1]
  field_simp [hrne]
  ring

theorem elementaryWeightTwoPrimitive_hasDerivAt26
    {r : ℝ} (hr0 : 0 < r) (hr1 : r < 1) :
    HasDerivAt elementaryWeightTwoPrimitive26
      (elementaryWeightTwoKernel26 r) r := by
  have hrne : r ≠ 0 := ne_of_gt hr0
  have hfactor :=
    elementaryLogFactor_hasDerivAt26 r (by linarith : 1 + r ≠ 0)
  have hprod := (Real.hasDerivAt_log hrne).mul hfactor
  have hcorr := elementaryDilogCorrection_hasDerivAt26 hr0 hr1
  convert hprod.add hcorr using 1
  simp only [elementaryWeightTwoKernel26]
  field_simp [hrne]
  ring

theorem elementaryWeightTwoKernel_intervalIntegrable26 :
    IntervalIntegrable elementaryWeightTwoKernel26
      MeasureTheory.volume 0 1 := by
  have hfirst :
      ContinuousOn
        (fun r : ℝ => (2 * r + 1) / (1 + r + r ^ 2))
        (Icc (0 : ℝ) 1) := by
    apply ContinuousOn.div
    · fun_prop
    · fun_prop
    · intro r _
      exact (cyclotomicQuadratic_pos26 r).ne'
  have hsecond :
      ContinuousOn (fun r : ℝ => 2 / (1 + r)) (Icc (0 : ℝ) 1) := by
    apply ContinuousOn.div
    · fun_prop
    · fun_prop
    · intro r hr
      linarith [hr.1]
  have hcoeff :
      ContinuousOn
        (fun r : ℝ =>
          (2 * r + 1) / (1 + r + r ^ 2) - 2 / (1 + r))
        (Icc (0 : ℝ) 1) :=
    hfirst.sub hsecond
  have hlog :
      IntervalIntegrable Real.log MeasureTheory.volume (0 : ℝ) 1 :=
    intervalIntegral.intervalIntegrable_log'
  have hint :=
    hlog.continuousOn_mul
      (by
        simpa [Set.uIcc_of_le (show (0 : ℝ) ≤ 1 by norm_num)] using hcoeff)
  simpa [elementaryWeightTwoKernel26] using hint

theorem elementaryWeightTwoPrimitive_tendsto_zero26 :
    Tendsto elementaryWeightTwoPrimitive26 (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  have hfactorDeriv :=
    elementaryLogFactor_hasDerivAt26 0 (by norm_num : (1 + (0 : ℝ)) ≠ 0)
  have hslope :
      Tendsto
        (fun r : ℝ => r⁻¹ * elementaryLogFactor26 r)
        (𝓝[>] 0)
        (𝓝
          (((2 * (0 : ℝ) + 1) / (1 + 0 + 0 ^ 2) -
            2 / (1 + 0)))) := by
    simpa [elementaryLogFactor_zero26] using
      hfactorDeriv.tendsto_slope_zero_right
  have hlogr :
      Tendsto (fun r : ℝ => Real.log r * r) (𝓝[>] 0) (𝓝 0) := by
    simpa [Real.rpow_one] using
      (tendsto_log_mul_rpow_nhdsGT_zero (show (0 : ℝ) < 1 by norm_num))
  have hproductRaw :=
    hlogr.mul hslope
  have hproductRaw' :
      Tendsto
        (fun r : ℝ =>
          (Real.log r * r) * (r⁻¹ * elementaryLogFactor26 r))
        (𝓝[>] 0) (𝓝 0) := by
    simpa using hproductRaw
  have hproduct :
      Tendsto
        (fun r : ℝ => Real.log r * elementaryLogFactor26 r)
        (𝓝[>] 0) (𝓝 0) := by
    apply Filter.Tendsto.congr' _ hproductRaw'
    filter_upwards [self_mem_nhdsWithin] with r hr
    have hrne : r ≠ 0 := ne_of_gt hr
    field_simp [hrne]
  have hcorrection :
      Tendsto elementaryDilogCorrection26 (𝓝[>] 0) (𝓝 0) := by
    have hwithin :
        ContinuousWithinAt elementaryDilogCorrection26 (Ioi (0 : ℝ)) 0 :=
      (elementaryDilogCorrection_continuousOn26 0 (by norm_num)).mono_of_mem_nhdsWithin
        (Icc_mem_nhdsGT (show (0 : ℝ) < 1 by norm_num))
    simpa [elementaryDilogCorrection_zero26] using hwithin.tendsto
  simpa [elementaryWeightTwoPrimitive26] using hproduct.add hcorrection

theorem elementaryWeightTwoPrimitive_tendsto_one26 :
    Tendsto elementaryWeightTwoPrimitive26 (𝓝[<] (1 : ℝ))
      (𝓝 (Real.pi ^ 2 / 18)) := by
  have hlog :
      Tendsto Real.log (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    simpa using tendsto_nhdsWithin_of_tendsto_nhds
      (Real.continuousAt_log (by norm_num : (1 : ℝ) ≠ 0)).tendsto
  have hfactorCont : ContinuousAt elementaryLogFactor26 1 := by
    unfold elementaryLogFactor26
    fun_prop (disch := norm_num)
  have hproduct :
      Tendsto
        (fun r : ℝ => Real.log r * elementaryLogFactor26 r)
        (𝓝[<] 1) (𝓝 0) := by
    convert hlog.mul
      (tendsto_nhdsWithin_of_tendsto_nhds hfactorCont.tendsto) using 1
    ring
  have hcorrection :
      Tendsto elementaryDilogCorrection26 (𝓝[<] 1)
        (𝓝 (Real.pi ^ 2 / 18)) := by
    have hwithin :
        ContinuousWithinAt elementaryDilogCorrection26 (Iio (1 : ℝ)) 1 :=
      (elementaryDilogCorrection_continuousOn26 1 (by norm_num)).mono_of_mem_nhdsWithin
        (Icc_mem_nhdsLT (show (0 : ℝ) < 1 by norm_num))
    simpa [elementaryDilogCorrection_one26] using hwithin.tendsto
  simpa [elementaryWeightTwoPrimitive26] using hproduct.add hcorrection

theorem elementaryWeightTwoIntegral26 :
    (∫ r : ℝ in 0..1, elementaryWeightTwoKernel26 r) =
      Real.pi ^ 2 / 18 := by
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_tendsto
    (f := elementaryWeightTwoPrimitive26)
    (fa := (0 : ℝ)) (fb := Real.pi ^ 2 / 18)
    (by norm_num)
    (fun r hr => elementaryWeightTwoPrimitive_hasDerivAt26 hr.1 hr.2)
    elementaryWeightTwoKernel_intervalIntegrable26
    elementaryWeightTwoPrimitive_tendsto_zero26
    elementaryWeightTwoPrimitive_tendsto_one26]
  ring

theorem cyclotomicWeightThreeKernel_intervalIntegrable26 :
    IntervalIntegrable cyclotomicWeightThreeKernel26
      MeasureTheory.volume 0 1 := by
  have hfactor :
      ContinuousOn
        (fun r : ℝ =>
          (2 * r + 1) / (1 + r + r ^ 2) * Real.log (1 + r))
        (Icc (0 : ℝ) 1) := by
    apply ContinuousOn.mul
    · apply ContinuousOn.div
      · fun_prop
      · fun_prop
      · intro r _
        exact (cyclotomicQuadratic_pos26 r).ne'
    · apply ContinuousOn.log
      · fun_prop
      · intro r hr
        linarith [hr.1]
  have hlog :
      IntervalIntegrable Real.log MeasureTheory.volume (0 : ℝ) 1 :=
    intervalIntegral.intervalIntegrable_log'
  have hint :=
    hlog.continuousOn_mul
      (by
        simpa [Set.uIcc_of_le (show (0 : ℝ) ≤ 1 by norm_num)] using hfactor)
  convert hint using 1
  funext r
  unfold cyclotomicWeightThreeKernel26
  ring

theorem alternatingWeightThreeKernel_intervalIntegrable26 :
    IntervalIntegrable alternatingWeightThreeKernel26
      MeasureTheory.volume 0 1 := by
  have hfactor :
      ContinuousOn
        (fun r : ℝ => Real.log (1 + r) / (1 + r))
        (Icc (0 : ℝ) 1) := by
    apply ContinuousOn.div
    · apply ContinuousOn.log
      · fun_prop
      · intro r hr
        linarith [hr.1]
    · fun_prop
    · intro r hr
      linarith [hr.1]
  have hlog :
      IntervalIntegrable Real.log MeasureTheory.volume (0 : ℝ) 1 :=
    intervalIntegral.intervalIntegrable_log'
  have hint :=
    hlog.mul_continuousOn
      (by
        simpa [Set.uIcc_of_le (show (0 : ℝ) ≤ 1 by norm_num)] using hfactor)
  convert hint using 1
  funext r
  unfold alternatingWeightThreeKernel26
  ring

/-- Pointwise splitting of the nested-sum integral into one elementary
weight-two term and three weight-three logarithmic terms. -/
theorem nestedCyclotomicKernel_decomposition26
    (r : ℝ) (hr0 : r ≠ 0) (hr1 : 1 + r ≠ 0) :
    nestedCyclotomicKernel26 r =
      3 * cyclotomicWeightThreeKernel26 r -
        8 * alternatingWeightThreeKernel26 r +
        radialWeightThreeKernel26 r +
        elementaryWeightTwoKernel26 r := by
  rw [nestedCyclotomicKernel26, cyclotomicWeightThreeKernel26,
    alternatingWeightThreeKernel26, radialWeightThreeKernel26,
    elementaryWeightTwoKernel26]
  rw [show
      (1 - r) * ((1 - r) ^ 2 * Real.log (1 + r) - r) * Real.log r /
          (r * (1 + r) * (1 + r + r ^ 2)) =
        ((1 - r) ^ 3 / (r * (1 + r) * (1 + r + r ^ 2))) *
            Real.log r * Real.log (1 + r) +
          (-(1 - r) / ((1 + r) * (1 + r + r ^ 2))) *
            Real.log r by
        have hq : 1 + r + r ^ 2 ≠ 0 :=
          (cyclotomicQuadratic_pos26 r).ne'
        field_simp [hr0, hr1, hq]
        ring]
  rw [cubicRationalDecomposition26 r hr0 hr1,
    linearRationalDecomposition26 r hr1]
  ring

theorem nestedCyclotomicKernel_intervalIntegrable26 :
    IntervalIntegrable nestedCyclotomicKernel26
      MeasureTheory.volume 0 1 := by
  have hright :
      IntervalIntegrable
        (fun r : ℝ =>
          3 * cyclotomicWeightThreeKernel26 r -
            8 * alternatingWeightThreeKernel26 r +
            radialWeightThreeKernel26 r +
            elementaryWeightTwoKernel26 r)
        MeasureTheory.volume 0 1 :=
    (((cyclotomicWeightThreeKernel_intervalIntegrable26.const_mul 3).sub
      (alternatingWeightThreeKernel_intervalIntegrable26.const_mul 8)).add
      radialWeightThreeKernel_intervalIntegrable26).add
      elementaryWeightTwoKernel_intervalIntegrable26
  apply IntervalIntegrable.congr
    (f := fun r : ℝ =>
      3 * cyclotomicWeightThreeKernel26 r -
        8 * alternatingWeightThreeKernel26 r +
        radialWeightThreeKernel26 r +
        elementaryWeightTwoKernel26 r) ?_ hright
  intro r hr
  have hr' : r ∈ Ioc (0 : ℝ) 1 := by
    simpa [Set.uIoc_of_le (show (0 : ℝ) ≤ 1 by norm_num)] using hr
  exact (nestedCyclotomicKernel_decomposition26 r
    (ne_of_gt hr'.1) (by linarith [hr'.1] : 1 + r ≠ 0)).symm

/-- All elementary contributions to the nested integral are now evaluated.
Only one weight-three cyclotomic combination remains. -/
theorem nestedCyclotomicIntegral_reduction26 :
    (∫ r : ℝ in 0..1, nestedCyclotomicKernel26 r) =
      3 * (∫ r : ℝ in 0..1, cyclotomicWeightThreeKernel26 r) -
        8 * (∫ r : ℝ in 0..1, alternatingWeightThreeKernel26 r) -
        (3 : ℝ) / 4 * zeta3 + Real.pi ^ 2 / 18 := by
  have hC := cyclotomicWeightThreeKernel_intervalIntegrable26
  have hA := alternatingWeightThreeKernel_intervalIntegrable26
  have hR := radialWeightThreeKernel_intervalIntegrable26
  have hE := elementaryWeightTwoKernel_intervalIntegrable26
  have heq :
      (∫ r : ℝ in 0..1, nestedCyclotomicKernel26 r) =
        ∫ r : ℝ in 0..1,
          (3 * cyclotomicWeightThreeKernel26 r -
            8 * alternatingWeightThreeKernel26 r +
            radialWeightThreeKernel26 r +
            elementaryWeightTwoKernel26 r) := by
    apply intervalIntegral.integral_congr_ae
    exact Filter.Eventually.of_forall fun r hr => by
      have hr' : r ∈ Ioc (0 : ℝ) 1 := by
        simpa [Set.uIoc_of_le (show (0 : ℝ) ≤ 1 by norm_num)] using hr
      exact nestedCyclotomicKernel_decomposition26 r
        (ne_of_gt hr'.1) (by linarith [hr'.1] : 1 + r ≠ 0)
  rw [heq]
  rw [intervalIntegral.integral_add
      (((hC.const_mul 3).sub (hA.const_mul 8)).add hR) hE,
    intervalIntegral.integral_add
      ((hC.const_mul 3).sub (hA.const_mul 8)) hR,
    intervalIntegral.integral_sub (hC.const_mul 3) (hA.const_mul 8),
    intervalIntegral.integral_const_mul,
    intervalIntegral.integral_const_mul,
    radialWeightThreeIntegral26,
    elementaryWeightTwoIntegral26]
  ring

/-- The sole remaining special-value input after the exact generating-function
and integral reductions. -/
def CyclotomicWeightThreeEvaluation26 : Prop :=
  3 * (∫ r : ℝ in 0..1, cyclotomicWeightThreeKernel26 r) -
      8 * (∫ r : ℝ in 0..1, alternatingWeightThreeKernel26 r) =
    (5 : ℝ) / 12 * zeta3

theorem nestedCyclotomicIntegral26_of_weightThree
    (hWeightThree : CyclotomicWeightThreeEvaluation26) :
    (∫ r : ℝ in 0..1, nestedCyclotomicKernel26 r) =
      Real.pi ^ 2 / 18 - zeta3 / 3 := by
  rw [nestedCyclotomicIntegral_reduction26]
  unfold CyclotomicWeightThreeEvaluation26 at hWeightThree
  linarith

end RamanujanChallenge.P26

end
