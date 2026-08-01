/-
  Generating functions used in Ramanujan Challenge Problem 2.6.

  This module starts from Mathlib's binomial power series and derives the
  ordinary generating function of the central binomial coefficients.  The
  weighted generating function for `inverseBinomialA26` is developed below.
-/
import RamanujanChallenge.Problem26Nested
import Mathlib.Analysis.Analytic.Binomial

noncomputable section

namespace RamanujanChallenge.P26

open Filter Set TopologicalSpace MeasureTheory

open scoped Topology Interval

/-- Consecutive generalized binomial coefficients satisfy the usual
first-order recurrence. -/
lemma ring_choose_succ_ratio (r : ℝ) (n : ℕ) :
    ((n : ℝ) + 1) * Ring.choose r (n + 1) =
      (r - n) * Ring.choose r n := by
  rw [Ring.choose_eq_smul, Ring.choose_eq_smul]
  simp only [smul_eq_mul]
  rw [descPochhammer_succ_right, Polynomial.smeval_mul,
    Polynomial.smeval_sub, Polynomial.smeval_X,
    Polynomial.smeval_natCast]
  rw [Nat.factorial_succ]
  push_cast
  have hn : ((n : ℝ) + 1) ≠ 0 := by positivity
  field_simp
  ring

/-- The coefficients of the binomial series of exponent `-1/2`, after the
substitution `x ↦ -4x`, are the central binomial coefficients. -/
lemma neg_half_choose_mul_neg_four_pow :
    ∀ n : ℕ,
      Ring.choose (-(1 : ℝ) / 2) n * (-4 : ℝ) ^ n =
        (Nat.centralBinom n : ℝ)
  | 0 => by norm_num [Ring.choose_zero_right]
  | n + 1 => by
      have hc := ring_choose_succ_ratio (-(1 : ℝ) / 2) n
      have hnat :
          ((n : ℝ) + 1) * (Nat.centralBinom (n + 1) : ℝ) =
            2 * (2 * (n : ℝ) + 1) * (Nat.centralBinom n : ℝ) := by
        exact_mod_cast Nat.succ_mul_centralBinom_succ n
      have hn : (n : ℝ) + 1 ≠ 0 := by positivity
      apply mul_left_cancel₀ hn
      rw [pow_succ]
      calc
        ((n : ℝ) + 1) *
              (Ring.choose (-(1 : ℝ) / 2) (n + 1) *
                ((-4 : ℝ) ^ n * -4)) =
            (((n : ℝ) + 1) *
              Ring.choose (-(1 : ℝ) / 2) (n + 1)) *
                ((-4 : ℝ) ^ n * -4) := by ring
        _ = ((-(1 : ℝ) / 2 - n) *
              Ring.choose (-(1 : ℝ) / 2) n) *
                ((-4 : ℝ) ^ n * -4) := by rw [hc]
        _ = 2 * (2 * (n : ℝ) + 1) *
              (Nat.centralBinom n : ℝ) := by
                rw [← neg_half_choose_mul_neg_four_pow n]
                ring
        _ = ((n : ℝ) + 1) *
              (Nat.centralBinom (n + 1) : ℝ) := hnat.symm

/-- The central-binomial generating function in real-power form. -/
theorem centralBinom_hasSum_rpow (t : ℝ) (ht : |4 * t| < 1) :
    HasSum (fun n : ℕ => (Nat.centralBinom n : ℝ) * t ^ n)
      ((1 - 4 * t) ^ (-(1 : ℝ) / 2)) := by
  have hy : -4 * t ∈ Metric.eball (0 : ℝ) 1 := by
    rw [Metric.mem_eball, edist_dist, Real.dist_eq]
    rw [ENNReal.ofReal_lt_one]
    simpa [abs_mul, mul_comm] using ht
  have h := (Real.one_add_rpow_hasFPowerSeriesOnBall_zero
    (a := -(1 : ℝ) / 2)).hasSum_sub
      (y := -4 * t) hy
  convert h using 1
  · funext n
    simp only [binomialSeries_apply, sub_zero, smul_eq_mul,
      List.prod_ofFn, neg_mul]
    rw [Finset.prod_const, Finset.card_univ, Fintype.card_fin]
    rw [show -(4 * t) = (-4 : ℝ) * t by ring, mul_pow,
      ← mul_assoc, neg_half_choose_mul_neg_four_pow]
  · ring_nf

/-- The ordinary generating function
`∑ n, centralBinom n * t^n = 1 / sqrt (1 - 4t)`. -/
theorem centralBinom_hasSum (t : ℝ) (ht : |4 * t| < 1) :
    HasSum (fun n : ℕ => (Nat.centralBinom n : ℝ) * t ^ n)
      (1 / Real.sqrt (1 - 4 * t)) := by
  have h := centralBinom_hasSum_rpow t ht
  convert h using 1
  have hpos : 0 < 1 - 4 * t := by
    have hle : 4 * t ≤ |4 * t| := le_abs_self _
    linarith
  rw [show (-(1 : ℝ) / 2) = -(1 / 2 : ℝ) by ring,
    Real.rpow_neg hpos.le, ← Real.sqrt_eq_rpow]
  simp [one_div]

/-- A coefficient identity that reduces the challenge weight to one
logarithmic central-binomial series and its shift. -/
lemma inverseBinomialA26_eq_central_difference
    (m : ℕ) (hm : 0 < m) :
    ((inverseBinomialA26 m : ℚ) : ℝ) =
      2 * (Nat.centralBinom m : ℝ) / m -
        (Nat.centralBinom (m + 1) : ℝ) / (2 * (m + 1)) := by
  have hnat := Nat.succ_mul_centralBinom_succ m
  have hreal :
      ((m : ℝ) + 1) * (Nat.centralBinom (m + 1) : ℝ) =
        2 * (2 * (m : ℝ) + 1) * (Nat.centralBinom m : ℝ) := by
    exact_mod_cast hnat
  simp only [inverseBinomialA26]
  push_cast
  rw [← Nat.centralBinom_eq_two_mul_choose]
  have hmR : (m : ℝ) ≠ 0 := by positivity
  have hm1R : (m : ℝ) + 1 ≠ 0 := by positivity
  field_simp
  nlinarith [hreal]

/-- The derivative kernel whose integral gives the positive-index logarithmic
central-binomial series. -/
def centralLogKernel26 (n : ℕ) : C(ℝ, ℝ) where
  toFun x := (Nat.centralBinom (n + 1) : ℝ) * x ^ n
  continuous_toFun := by fun_prop

/-- The integrated positive-index central-binomial term. -/
def centralLogTerm26 (n : ℕ) (t : ℝ) : ℝ :=
  (Nat.centralBinom (n + 1) : ℝ) * t ^ (n + 1) / (n + 1 : ℝ)

theorem centralLogKernel26_integral (n : ℕ) (t : ℝ) :
    (∫ x : ℝ in 0..t, centralLogKernel26 n x) =
      centralLogTerm26 n t := by
  rw [show (fun x : ℝ => centralLogKernel26 n x) =
      fun x => (Nat.centralBinom (n + 1) : ℝ) * x ^ n by rfl,
    intervalIntegral.integral_const_mul, integral_pow]
  simp [centralLogTerm26]
  ring

theorem centralLogKernel26_norm_le (n : ℕ) (t : ℝ) (ht : 0 ≤ t) :
    ‖(centralLogKernel26 n).restrict
        (⟨Set.uIcc (0 : ℝ) t, isCompact_uIcc⟩ : Compacts ℝ)‖
      ≤ 4 * (4 * t) ^ n := by
  rw [ContinuousMap.norm_le _ (by positivity)]
  intro x
  rcases x.property with ⟨hx0, hxt⟩
  simp only [min_eq_left ht, max_eq_right ht] at hx0 hxt
  simp only [centralLogKernel26, ContinuousMap.coe_mk, Compacts.coe_mk,
    ContinuousMap.restrict_apply, Real.norm_eq_abs]
  have hx : 0 ≤ x.1 := hx0
  have hterm :
      0 ≤ (Nat.centralBinom (n + 1) : ℝ) * x.1 ^ n := by
    positivity
  rw [abs_of_nonneg hterm]
  have hc :
      (Nat.centralBinom (n + 1) : ℝ) ≤ (4 : ℝ) ^ (n + 1) := by
    exact_mod_cast Nat.centralBinom_le_four_pow (n + 1)
  calc
    (Nat.centralBinom (n + 1) : ℝ) * x.1 ^ n
        ≤ (4 : ℝ) ^ (n + 1) * t ^ n := by gcongr
    _ = 4 * (4 * t) ^ n := by
      rw [pow_succ, mul_pow]
      ring

theorem centralLogKernel26_norm_summable
    (t : ℝ) (ht0 : 0 ≤ t) (ht4 : 4 * t < 1) :
    Summable (fun n : ℕ =>
      ‖(centralLogKernel26 n).restrict
        (⟨Set.uIcc (0 : ℝ) t, isCompact_uIcc⟩ : Compacts ℝ)‖) := by
  apply Summable.of_nonneg_of_le
  · intro n
    positivity
  · exact fun n => centralLogKernel26_norm_le n t ht0
  · exact (summable_geometric_of_norm_lt_one
      (by
        rw [Real.norm_eq_abs,
          abs_of_nonneg (mul_nonneg (by norm_num) ht0)]
        exact ht4)).mul_left 4

/-- The shifted central-binomial derivative series, rationalized so that its
closed form is continuous at zero. -/
theorem centralLogKernel26_hasSum
    (x : ℝ) (hx0 : 0 ≤ x) (hx4 : 4 * x < 1) :
    HasSum (fun n : ℕ => centralLogKernel26 n x)
      (4 / (Real.sqrt (1 - 4 * x) *
        (1 + Real.sqrt (1 - 4 * x)))) := by
  by_cases hx : x = 0
  · subst x
    have hsingle :
        HasSum (fun n : ℕ => centralLogKernel26 n 0)
          (centralLogKernel26 0 0) := by
      apply hasSum_single 0
      intro n hn
      simp only [centralLogKernel26, ContinuousMap.coe_mk]
      rw [zero_pow hn, mul_zero]
    convert hsingle using 1
    norm_num [centralLogKernel26, Nat.centralBinom]
  · have hxpos : 0 < x := lt_of_le_of_ne hx0 (Ne.symm hx)
    have habs : |4 * x| < 1 := by
      rw [abs_of_nonneg (mul_nonneg (by norm_num) hx0)]
      exact hx4
    have hfull := centralBinom_hasSum x habs
    have hshift :
        HasSum
          (fun n : ℕ =>
            (Nat.centralBinom (n + 1) : ℝ) * x ^ (n + 1))
          (1 / Real.sqrt (1 - 4 * x) - 1) := by
      apply (hasSum_nat_add_iff
        (f := fun n : ℕ => (Nat.centralBinom n : ℝ) * x ^ n) 1).mpr
      convert hfull using 1
      simp
    have hdiv := hshift.div_const x
    convert hdiv using 1
    · funext n
      simp only [centralLogKernel26, ContinuousMap.coe_mk]
      rw [pow_succ]
      field_simp
    · have hrad : 0 < 1 - 4 * x := by linarith
      have hsqrt : 0 < Real.sqrt (1 - 4 * x) :=
        Real.sqrt_pos.2 hrad
      have hsq :
          Real.sqrt (1 - 4 * x) ^ 2 = 1 - 4 * x :=
        Real.sq_sqrt hrad.le
      field_simp [hx, ne_of_gt hsqrt]
      nlinarith

def centralLogIntegrand26 (x : ℝ) : ℝ :=
  4 / (Real.sqrt (1 - 4 * x) * (1 + Real.sqrt (1 - 4 * x)))

theorem centralLogTerm26_hasSum_integral
    (t : ℝ) (ht0 : 0 ≤ t) (ht4 : 4 * t < 1) :
    HasSum (fun n : ℕ => centralLogTerm26 n t)
      (∫ x : ℝ in 0..t, centralLogIntegrand26 x) := by
  have h := intervalIntegral.hasSum_intervalIntegral_of_summable_norm
    (a := (0 : ℝ)) (b := t)
    (centralLogKernel26_norm_summable t ht0 ht4)
  have hi :
      HasSum (fun n : ℕ => centralLogTerm26 n t)
        (∫ x : ℝ in 0..t, ∑' n : ℕ, centralLogKernel26 n x) := by
    apply h.congr_fun
    intro n
    exact (centralLogKernel26_integral n t).symm
  convert hi using 1
  apply intervalIntegral.integral_congr_ae
  exact Filter.Eventually.of_forall fun x hx => by
    have hx' : x ∈ Set.Ioc (0 : ℝ) t := by
      simpa [Set.uIoc_of_le ht0] using hx
    have hx4 : 4 * x < 1 := by
      nlinarith [
        mul_le_mul_of_nonneg_left hx'.2 (by norm_num : (0 : ℝ) ≤ 4)]
    exact (centralLogKernel26_hasSum x hx'.1.le hx4).tsum_eq.symm

def centralLogPrimitive26 (x : ℝ) : ℝ :=
  -2 * Real.log (1 + Real.sqrt (1 - 4 * x))

theorem centralLogIntegrand26_continuousOn
    (t : ℝ) (ht4 : 4 * t < 1) :
    ContinuousOn centralLogIntegrand26 (Set.Icc 0 t) := by
  intro x hx
  have hrad : 0 < 1 - 4 * x := by nlinarith [hx.2]
  have hsqrt : 0 < Real.sqrt (1 - 4 * x) :=
    Real.sqrt_pos.2 hrad
  have hsqrt_cont :
      ContinuousAt (fun y : ℝ => Real.sqrt (1 - 4 * y)) x :=
    by
      simpa only [Function.comp_apply] using
        Real.continuous_sqrt.continuousAt.comp
          (show ContinuousAt (fun y : ℝ => 1 - 4 * y) x by fun_prop)
  have hden_cont :
      ContinuousAt
        (fun y : ℝ =>
          Real.sqrt (1 - 4 * y) * (1 + Real.sqrt (1 - 4 * y))) x :=
    hsqrt_cont.mul ((continuousAt_const.add hsqrt_cont))
  have hden_ne :
      Real.sqrt (1 - 4 * x) * (1 + Real.sqrt (1 - 4 * x)) ≠ 0 := by
    positivity
  unfold centralLogIntegrand26
  exact (continuousAt_const.div hden_cont hden_ne).continuousWithinAt

theorem centralLogPrimitive26_continuousOn
    (t : ℝ) (ht4 : 4 * t < 1) :
    ContinuousOn centralLogPrimitive26 (Set.Icc 0 t) := by
  intro x hx
  have hrad : 0 < 1 - 4 * x := by nlinarith [hx.2]
  have hsqrt : 0 < Real.sqrt (1 - 4 * x) :=
    Real.sqrt_pos.2 hrad
  have hsqrt_cont :
      ContinuousAt (fun y : ℝ => Real.sqrt (1 - 4 * y)) x :=
    by
      simpa only [Function.comp_apply] using
        Real.continuous_sqrt.continuousAt.comp
          (show ContinuousAt (fun y : ℝ => 1 - 4 * y) x by fun_prop)
  have hone_cont :
      ContinuousAt (fun y : ℝ => 1 + Real.sqrt (1 - 4 * y)) x :=
    continuousAt_const.add hsqrt_cont
  have hone_ne : 1 + Real.sqrt (1 - 4 * x) ≠ 0 := by positivity
  unfold centralLogPrimitive26
  exact (continuousAt_const.mul
    (hone_cont.log hone_ne)).continuousWithinAt

theorem centralLogPrimitive26_hasDerivAt
    (x : ℝ) (hrad : 0 < 1 - 4 * x) :
    HasDerivAt centralLogPrimitive26 (centralLogIntegrand26 x) x := by
  have hinner :
      HasDerivAt (fun y : ℝ => 1 - 4 * y) (-4) x := by
    convert (hasDerivAt_const x 1).sub
      ((hasDerivAt_id x).const_mul 4) using 1
    ring
  have hsqrt :
      HasDerivAt (fun y : ℝ => Real.sqrt (1 - 4 * y))
        ((1 / (2 * Real.sqrt (1 - 4 * x))) * (-4)) x :=
    (Real.hasDerivAt_sqrt (ne_of_gt hrad)).comp x hinner
  have hone :
      HasDerivAt (fun y : ℝ => 1 + Real.sqrt (1 - 4 * y))
        ((1 / (2 * Real.sqrt (1 - 4 * x))) * (-4)) x :=
    by
      convert (hasDerivAt_const x 1).add hsqrt using 1
      ring
  have hlog := hone.log (by positivity)
  have hscaled := hlog.const_mul (-2)
  unfold centralLogPrimitive26 centralLogIntegrand26
  convert hscaled using 1
  have hsqrtpos : 0 < Real.sqrt (1 - 4 * x) :=
    Real.sqrt_pos.2 hrad
  field_simp [ne_of_gt hsqrtpos]

theorem centralLogIntegral26
    (t : ℝ) (ht0 : 0 ≤ t) (ht4 : 4 * t < 1) :
    (∫ x : ℝ in 0..t, centralLogIntegrand26 x) =
      2 * Real.log
        (2 / (1 + Real.sqrt (1 - 4 * t))) := by
  have hint :
      IntervalIntegrable centralLogIntegrand26
        MeasureTheory.volume 0 t :=
    (by
      apply ContinuousOn.intervalIntegrable
      simpa [Set.uIcc_of_le ht0] using
        centralLogIntegrand26_continuousOn t ht4)
  have hftc :=
    intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le
      (a := (0 : ℝ)) (b := t) ht0
      (centralLogPrimitive26_continuousOn t ht4)
      (fun x hx =>
        centralLogPrimitive26_hasDerivAt x (by nlinarith [hx.2]))
      hint
  rw [hftc]
  simp only [centralLogPrimitive26]
  norm_num
  have hden :
      1 + Real.sqrt (1 - 4 * t) ≠ 0 := by
    have hrad : 0 < 1 - 4 * t := by linarith
    have : 0 < Real.sqrt (1 - 4 * t) := Real.sqrt_pos.2 hrad
    positivity
  rw [Real.log_div (by norm_num) hden]
  ring

theorem centralLogTerm26_hasSum
    (t : ℝ) (ht0 : 0 ≤ t) (ht4 : 4 * t < 1) :
    HasSum (fun n : ℕ => centralLogTerm26 n t)
      (2 * Real.log
        (2 / (1 + Real.sqrt (1 - 4 * t)))) := by
  convert centralLogTerm26_hasSum_integral t ht0 ht4 using 1
  exact (centralLogIntegral26 t ht0 ht4).symm

/-- The standard Catalan substitution. -/
def catalanParameter26 (r : ℝ) : ℝ :=
  r / (1 + r) ^ 2

theorem catalanParameter26_pos
    {r : ℝ} (hr0 : 0 < r) :
    0 < catalanParameter26 r := by
  unfold catalanParameter26
  positivity

theorem four_mul_catalanParameter26_lt_one
    {r : ℝ} (hr0 : 0 < r) (hr1 : r < 1) :
    4 * catalanParameter26 r < 1 := by
  have hden : 0 < (1 + r) ^ 2 := by positivity
  unfold catalanParameter26
  rw [← mul_div_assoc, div_lt_iff₀ hden]
  nlinarith [sq_pos_of_pos (sub_pos.mpr hr1)]

theorem sqrt_one_sub_four_catalanParameter26
    {r : ℝ} (hr0 : 0 < r) (hr1 : r < 1) :
    Real.sqrt (1 - 4 * catalanParameter26 r) =
      (1 - r) / (1 + r) := by
  have ht4 := four_mul_catalanParameter26_lt_one hr0 hr1
  have hinput : 0 ≤ 1 - 4 * catalanParameter26 r := by linarith
  have hrhs : 0 ≤ (1 - r) / (1 + r) :=
    div_nonneg (by linarith) (by linarith)
  apply (Real.sqrt_eq_iff_eq_sq hinput hrhs).2
  unfold catalanParameter26
  have hden : 1 + r ≠ 0 := by positivity
  field_simp [hden]
  ring

theorem centralLogTerm26_hasSum_catalanParameter
    {r : ℝ} (hr0 : 0 < r) (hr1 : r < 1) :
    HasSum
      (fun n : ℕ => centralLogTerm26 n (catalanParameter26 r))
      (2 * Real.log (1 + r)) := by
  have ht0 := (catalanParameter26_pos hr0).le
  have ht4 := four_mul_catalanParameter26_lt_one hr0 hr1
  convert centralLogTerm26_hasSum (catalanParameter26 r) ht0 ht4 using 1
  rw [sqrt_one_sub_four_catalanParameter26 hr0 hr1]
  have hden : 1 + r ≠ 0 := by positivity
  congr 2
  field_simp [hden]
  ring

theorem inverseBinomialA26_succ_term_eq
    (n : ℕ) (t : ℝ) (ht : t ≠ 0) :
    ((inverseBinomialA26 (n + 1) : ℚ) : ℝ) * t ^ (n + 1) =
      2 * centralLogTerm26 n t -
        centralLogTerm26 (n + 1) t / (2 * t) := by
  rw [inverseBinomialA26_eq_central_difference (n + 1) (by omega)]
  unfold centralLogTerm26
  push_cast
  field_simp [ht]
  ring

/-- The exact generating function for the challenge coefficients:

`∑ m, A_m (r/(1+r)^2)^m =
  1 - ((1-r)^2/r) log(1+r)` for `0 < r < 1`. -/
theorem inverseBinomialA26_hasSum_catalanParameter
    {r : ℝ} (hr0 : 0 < r) (hr1 : r < 1) :
    HasSum
      (fun m : ℕ =>
        ((inverseBinomialA26 m : ℚ) : ℝ) *
          catalanParameter26 r ^ m)
      (1 - (1 - r) ^ 2 / r * Real.log (1 + r)) := by
  let t := catalanParameter26 r
  let L := 2 * Real.log (1 + r)
  have htpos : 0 < t := catalanParameter26_pos hr0
  have ht : t ≠ 0 := ne_of_gt htpos
  have hlog :
      HasSum (fun n : ℕ => centralLogTerm26 n t) L := by
    simpa [t, L] using
      centralLogTerm26_hasSum_catalanParameter hr0 hr1
  have hlogTail :
      HasSum (fun n : ℕ => centralLogTerm26 (n + 1) t)
        (L - centralLogTerm26 0 t) := by
    apply (hasSum_nat_add_iff
      (f := fun n : ℕ => centralLogTerm26 n t) 1).mpr
    convert hlog using 1
    simp
  have hcombined :
      HasSum
        (fun n : ℕ =>
          2 * centralLogTerm26 n t -
            centralLogTerm26 (n + 1) t / (2 * t))
        (2 * L - (L - centralLogTerm26 0 t) / (2 * t)) :=
    (hlog.mul_left 2).sub (hlogTail.div_const (2 * t))
  have hshift :
      HasSum
        (fun n : ℕ =>
          ((inverseBinomialA26 (n + 1) : ℚ) : ℝ) * t ^ (n + 1))
        (1 - (1 - r) ^ 2 / r * Real.log (1 + r)) := by
    convert hcombined using 1
    · funext n
      exact inverseBinomialA26_succ_term_eq n t ht
    · have hrne : r ≠ 0 := ne_of_gt hr0
      have h1r : 1 + r ≠ 0 := by positivity
      simp only [centralLogTerm26, Nat.centralBinom, t, L]
      norm_num
      unfold catalanParameter26
      field_simp [hrne, h1r]
      ring
  have hfull :=
    (hasSum_nat_add_iff
      (f := fun m : ℕ =>
        ((inverseBinomialA26 m : ℚ) : ℝ) * t ^ m) 1).mp hshift
  simpa [t, inverseBinomialA26_zero] using hfull

/-! ## Continuation blueprint for the nested sum

The definitions below fix the indexing for the analytic bridge.  The compiled
proof later in this section imports the cyclotomic kernel from
`Problem26Nested`.
-/

/-- The real prefix `Sₖ = ∑_{m<k} Aₘ`. -/
def inverseBinomialPrefixReal26 (k : ℕ) : ℝ :=
  ∑ m ∈ Finset.range k, ((inverseBinomialA26 m : ℚ) : ℝ)

/-- After pairing the beta integral at `x` and `1-x`, and then making the
Catalan substitution `r = x/(1-x)`, the `(n+1)`-st nested summand is the
integral of this function on `[0,1]`.

In paper notation this is
`Fₖ(r) = Sₖ r^(k-1) / (k (1+r)^(2k))`, with `k = n+1`. -/
def nestedBetaTerm26 (n : ℕ) (r : ℝ) : ℝ :=
  inverseBinomialPrefixReal26 (n + 1) / (n + 1 : ℝ) *
    r ^ n / (1 + r) ^ (2 * (n + 1))

/-- The `Aₘ` ordinary generating series, as an actual `tsum`. -/
def inverseBinomialAGF26 (t : ℝ) : ℝ :=
  ∑' m : ℕ, ((inverseBinomialA26 m : ℚ) : ℝ) * t ^ m

/-- If `H(t) = ∑_{k≥1} Sₖ t^k/k`, then
`H(t) = ∫₀ᵗ A(u)/(1-u) du`. -/
def nestedGeneratingPrimitive26 (t : ℝ) : ℝ :=
  ∫ u : ℝ in 0..t, inverseBinomialAGF26 u / (1 - u)

/-- The derivative of `H (r/(1+r)^2)` after inserting the exact `A`
generating function proved above. -/
def nestedCatalanDerivative26 (r : ℝ) : ℝ :=
  (1 - r) *
      (r - (1 - r) ^ 2 * Real.log (1 + r)) /
    (r * (1 + r) * (1 + r + r ^ 2))

theorem inverseBinomialA26_nonneg (m : ℕ) :
    0 ≤ ((inverseBinomialA26 m : ℚ) : ℝ) := by
  simp only [inverseBinomialA26]
  push_cast
  positivity

theorem inverseBinomialPrefixReal26_nonneg (k : ℕ) :
    0 ≤ inverseBinomialPrefixReal26 k := by
  unfold inverseBinomialPrefixReal26
  exact Finset.sum_nonneg fun m _ => inverseBinomialA26_nonneg m

theorem nestedBetaTerm26_nonneg (n : ℕ) {r : ℝ}
    (hr0 : 0 ≤ r) :
    0 ≤ nestedBetaTerm26 n r := by
  have hprefix := inverseBinomialPrefixReal26_nonneg (n + 1)
  unfold nestedBetaTerm26
  positivity

/-- The Catalan change of variables `x = r/(1+r)` on the lower half of
the symmetric beta integral. -/
theorem catalanKernel_integral_half26 (n : ℕ) :
    (∫ r : ℝ in 0..1, r ^ n / (1 + r) ^ (2 * (n + 1))) =
      ∫ x : ℝ in 0..(1 / 2 : ℝ), x ^ n * (1 - x) ^ n := by
  let f : ℝ → ℝ := fun r => r / (1 + r)
  let f' : ℝ → ℝ := fun r => 1 / (1 + r) ^ 2
  let g : ℝ → ℝ := fun x => x ^ n * (1 - x) ^ n
  have hf : ∀ r ∈ Set.uIcc (0 : ℝ) 1, HasDerivAt f (f' r) r := by
    intro r hr
    have hden : 1 + r ≠ 0 := by
      simp only [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1),
        Set.mem_Icc] at hr
      linarith
    dsimp [f, f']
    convert
      (hasDerivAt_id r).div
        ((hasDerivAt_const r 1).add (hasDerivAt_id r)) hden using 1
    all_goals simp
  have hf' : ContinuousOn f' (Set.uIcc (0 : ℝ) 1) := by
    intro r hr
    have hden : 1 + r ≠ 0 := by
      simp only [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1),
        Set.mem_Icc] at hr
      linarith
    unfold f'
    exact
      (continuousAt_const.div
        ((continuousAt_const.add continuousAt_id).pow 2)
        (pow_ne_zero 2 hden)).continuousWithinAt
  have hg : Continuous g := by
    unfold g
    fun_prop
  have hsub :=
    intervalIntegral.integral_comp_mul_deriv
      (a := (0 : ℝ)) (b := 1) hf hf' hg
  have hleft :
      (∫ r : ℝ in 0..1, r ^ n / (1 + r) ^ (2 * (n + 1))) =
        ∫ r : ℝ in 0..1, (g ∘ f) r * f' r := by
    apply intervalIntegral.integral_congr_ae
    exact Filter.Eventually.of_forall fun r hr => by
      have hr' : r ∈ Set.Ioc (0 : ℝ) 1 := by
        simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hr
      have hden : 1 + r ≠ 0 := by linarith [hr'.1]
      dsimp [f, f', g, Function.comp_apply]
      have hone : 1 - r / (1 + r) = 1 / (1 + r) := by
        field_simp [hden]
        ring
      rw [hone]
      simp only [div_pow]
      field_simp [hden, pow_add, pow_mul]
      ring
  rw [hleft, hsub]
  norm_num [f, g]

/-- Pairing the beta kernel at `x` and `1-x` turns it into the Catalan
kernel on `[0,1]`. -/
theorem betaKernel_eq_catalanKernel_integral26 (n : ℕ) :
    (∫ x : ℝ in 0..1, x ^ n * (1 - x) ^ (n + 1)) =
      ∫ r : ℝ in 0..1, r ^ n / (1 + r) ^ (2 * (n + 1)) := by
  let p : ℝ → ℝ := fun x => x ^ n * (1 - x) ^ (n + 1)
  let q : ℝ → ℝ := fun x => p (1 - x)
  let g : ℝ → ℝ := fun x => x ^ n * (1 - x) ^ n
  have hp : Continuous p := by
    unfold p
    fun_prop
  have hq : Continuous q := by
    unfold q
    fun_prop
  have hp0 :
      IntervalIntegrable p MeasureTheory.volume (0 : ℝ) (1 / 2) :=
    hp.intervalIntegrable _ _
  have hp1 :
      IntervalIntegrable p MeasureTheory.volume (1 / 2 : ℝ) 1 :=
    hp.intervalIntegrable _ _
  have hq0 :
      IntervalIntegrable q MeasureTheory.volume (0 : ℝ) (1 / 2) :=
    hq.intervalIntegrable _ _
  have hpair :
      (∫ x : ℝ in (1 / 2 : ℝ)..1, p x) =
        ∫ x : ℝ in 0..(1 / 2 : ℝ), q x := by
    rw [show (∫ x : ℝ in 0..(1 / 2 : ℝ), q x) =
        ∫ x : ℝ in 0..(1 / 2 : ℝ), p (1 - x) by rfl]
    rw [intervalIntegral.integral_comp_sub_left]
    norm_num
  have hsum :
      (∫ x : ℝ in 0..(1 / 2 : ℝ), p x) +
          (∫ x : ℝ in 0..(1 / 2 : ℝ), q x) =
        ∫ x : ℝ in 0..(1 / 2 : ℝ), g x := by
    rw [← intervalIntegral.integral_add hp0 hq0]
    apply intervalIntegral.integral_congr
    intro x _
    dsimp [p, q, g]
    ring
  calc
    (∫ x : ℝ in 0..1, x ^ n * (1 - x) ^ (n + 1)) =
        (∫ x : ℝ in 0..1, p x) := by rfl
    _ = (∫ x : ℝ in 0..(1 / 2 : ℝ), p x) +
          (∫ x : ℝ in (1 / 2 : ℝ)..1, p x) :=
      (intervalIntegral.integral_add_adjacent_intervals hp0 hp1).symm
    _ = (∫ x : ℝ in 0..(1 / 2 : ℝ), p x) +
          (∫ x : ℝ in 0..(1 / 2 : ℝ), q x) := by rw [hpair]
    _ = ∫ x : ℝ in 0..(1 / 2 : ℝ), g x := hsum
    _ = ∫ r : ℝ in 0..1,
          r ^ n / (1 + r) ^ (2 * (n + 1)) :=
      (catalanKernel_integral_half26 n).symm

/-- Exact integral representation of the `(n+1)`-st nested summand. -/
theorem nestedBetaTerm26_integral (n : ℕ) :
    (∫ r : ℝ in 0..1, nestedBetaTerm26 n r) =
      ((inverseBinomialDTerm26 (n + 1) : ℚ) : ℝ) := by
  have hC := inverseBinomialC26_eq_integral (n + 1) (by omega)
  have hC' :
      ((inverseBinomialC26 (n + 1) : ℚ) : ℝ) =
        (∫ x : ℝ in 0..1, x ^ n * (1 - x) ^ (n + 1)) /
          (n + 1 : ℝ) := by
    rw [← intervalIntegral.integral_div]
    simpa only [Nat.add_sub_cancel, Nat.cast_add, Nat.cast_one] using hC
  have hJ :
      (∫ r : ℝ in 0..1,
          r ^ n / (1 + r) ^ (2 * (n + 1))) / (n + 1 : ℝ) =
        ((inverseBinomialC26 (n + 1) : ℚ) : ℝ) := by
    rw [← betaKernel_eq_catalanKernel_integral26]
    exact hC'.symm
  have hfun :
      (fun r : ℝ => nestedBetaTerm26 n r) =
        fun r : ℝ =>
          (inverseBinomialPrefixReal26 (n + 1) / (n + 1 : ℝ)) *
            (r ^ n / (1 + r) ^ (2 * (n + 1))) := by
    funext r
    unfold nestedBetaTerm26
    ring
  rw [hfun, intervalIntegral.integral_const_mul]
  rw [show inverseBinomialPrefixReal26 (n + 1) / (n + 1 : ℝ) *
        (∫ r : ℝ in 0..1,
          r ^ n / (1 + r) ^ (2 * (n + 1))) =
      inverseBinomialPrefixReal26 (n + 1) *
        ((∫ r : ℝ in 0..1,
          r ^ n / (1 + r) ^ (2 * (n + 1))) /
            (n + 1 : ℝ)) by ring,
    hJ]
  simp only [inverseBinomialDTerm26]
  push_cast
  unfold inverseBinomialPrefixReal26
  ring

/-- The integral of the norm is exactly the nested summand; this is the
non-wasteful majorant used for exchanging the sum and integral. -/
theorem nestedBetaTerm26_integral_norm (n : ℕ) :
    (∫ r : ℝ in 0..1, ‖nestedBetaTerm26 n r‖) =
      ((inverseBinomialDTerm26 (n + 1) : ℚ) : ℝ) := by
  rw [← nestedBetaTerm26_integral n]
  apply intervalIntegral.integral_congr
  intro r hr
  change ‖nestedBetaTerm26 n r‖ = nestedBetaTerm26 n r
  rw [Real.norm_eq_abs,
    abs_of_nonneg
      (nestedBetaTerm26_nonneg n
        (by
          simp only [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1),
            Set.mem_Icc] at hr
          exact hr.1))]

theorem nestedBetaTerm26_intervalIntegrable (n : ℕ) :
    IntervalIntegrable (nestedBetaTerm26 n)
      MeasureTheory.volume 0 1 := by
  apply ContinuousOn.intervalIntegrable
  intro r hr
  have hden : 1 + r ≠ 0 := by
    simp only [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1),
      Set.mem_Icc] at hr
    linarith
  unfold nestedBetaTerm26
  exact
    ((continuousAt_const.div_const _).mul (continuousAt_id.pow n)).div
      ((continuousAt_const.add continuousAt_id).pow (2 * (n + 1)))
      (pow_ne_zero _ hden) |>.continuousWithinAt

theorem nestedBetaTerm26_integral_norm_summable :
    Summable
      (fun n : ℕ =>
        ∫ r : ℝ in 0..1, ‖nestedBetaTerm26 n r‖) := by
  have hshift :
      Summable
        (fun n : ℕ =>
          ((inverseBinomialDTerm26 (n + 1) : ℚ) : ℝ)) := by
    simpa using
      (summable_nat_add_iff (f := fun k : ℕ =>
        ((inverseBinomialDTerm26 k : ℚ) : ℝ)) 1).2
        inverseBinomialDTerm26_summable
  exact hshift.congr fun n => (nestedBetaTerm26_integral_norm n).symm

/-- The first rigorous sum/integral bridge for the nested series. -/
theorem inverseBinomialDTerm26_tail_hasSum_nestedBeta_integral :
    HasSum
      (fun n : ℕ => ((inverseBinomialDTerm26 (n + 1) : ℚ) : ℝ))
      (∫ r : ℝ in Set.Ioc 0 1, ∑' n : ℕ, nestedBetaTerm26 n r) := by
  have hInt :
      ∀ n : ℕ, MeasureTheory.Integrable (nestedBetaTerm26 n)
        (MeasureTheory.volume.restrict (Set.Ioc 0 1)) := by
    intro n
    exact (nestedBetaTerm26_intervalIntegrable n).1
  have hNorm :
      Summable
        (fun n : ℕ =>
          ∫ r : ℝ in Set.Ioc 0 1, ‖nestedBetaTerm26 n r‖) := by
    simpa only [
      ← intervalIntegral.integral_of_le
        (by norm_num : (0 : ℝ) ≤ 1)] using
      nestedBetaTerm26_integral_norm_summable
  have h :=
    MeasureTheory.hasSum_integral_of_summable_integral_norm
      (μ := MeasureTheory.volume.restrict (Set.Ioc (0 : ℝ) 1))
      hInt hNorm
  convert h using 1
  funext n
  rw [← intervalIntegral.integral_of_le
    (by norm_num : (0 : ℝ) ≤ 1)]
  exact (nestedBetaTerm26_integral n).symm

theorem inverseBinomialA26_norm_summable_catalan
    {r : ℝ} (hr0 : 0 < r) (hr1 : r < 1) :
    Summable
      (fun m : ℕ =>
        ‖((inverseBinomialA26 m : ℚ) : ℝ) *
          catalanParameter26 r ^ m‖) := by
  have h :=
    (inverseBinomialA26_hasSum_catalanParameter hr0 hr1).summable
  apply h.congr
  intro m
  rw [Real.norm_eq_abs, abs_of_nonneg]
  exact mul_nonneg (inverseBinomialA26_nonneg m)
    (pow_nonneg (catalanParameter26_pos hr0).le m)

/-- The prefix generating function
`∑ n, S_(n+1) t^n = A(t)/(1-t)` at the Catalan parameter. -/
theorem inverseBinomialPrefixReal26_hasSum_catalan
    {r : ℝ} (hr0 : 0 < r) (hr1 : r < 1) :
    HasSum
      (fun n : ℕ =>
        inverseBinomialPrefixReal26 (n + 1) *
          catalanParameter26 r ^ n)
      ((1 - (1 - r) ^ 2 / r * Real.log (1 + r)) /
        (1 - catalanParameter26 r)) := by
  let t := catalanParameter26 r
  have ht0 : 0 < t := catalanParameter26_pos hr0
  have ht1 : t < 1 := by
    have ht4 := four_mul_catalanParameter26_lt_one hr0 hr1
    linarith
  have hA :=
    inverseBinomialA26_hasSum_catalanParameter hr0 hr1
  have hGeom :
      HasSum (fun n : ℕ => t ^ n) (1 / (1 - t)) := by
    simpa [one_div] using hasSum_geometric_of_norm_lt_one
      (show ‖t‖ < 1 by
        rw [Real.norm_eq_abs, abs_of_pos ht0]
        exact ht1)
  have hprod :
      HasSum
        (fun n : ℕ =>
          ∑ k ∈ Finset.range (n + 1),
            ((((inverseBinomialA26 k : ℚ) : ℝ) * t ^ k) *
              t ^ (n - k)))
        ((∑' m : ℕ, ((inverseBinomialA26 m : ℚ) : ℝ) * t ^ m) *
          ∑' n : ℕ, t ^ n) :=
    hasSum_sum_range_mul_of_summable_norm
      (show Summable
          (fun m : ℕ =>
            ‖((inverseBinomialA26 m : ℚ) : ℝ) * t ^ m‖) by
        simpa [t] using
          inverseBinomialA26_norm_summable_catalan hr0 hr1)
      (show Summable (fun n : ℕ => ‖t ^ n‖) from
        (hGeom.summable.norm))
  have hterm :
      (fun n : ℕ =>
          ∑ k ∈ Finset.range (n + 1),
            (((inverseBinomialA26 k : ℚ) : ℝ) * t ^ k) *
              t ^ (n - k)) =
        fun n : ℕ => inverseBinomialPrefixReal26 (n + 1) * t ^ n := by
    funext n
    calc
      (∑ k ∈ Finset.range (n + 1),
          (((inverseBinomialA26 k : ℚ) : ℝ) * t ^ k) *
            t ^ (n - k)) =
          ∑ k ∈ Finset.range (n + 1),
            ((inverseBinomialA26 k : ℚ) : ℝ) * t ^ n := by
              apply Finset.sum_congr rfl
              intro k hk
              have hkn : k ≤ n := by
                simpa [Finset.mem_range] using hk
              have hp : t ^ k * t ^ (n - k) = t ^ n := by
                rw [← pow_add, Nat.add_sub_of_le hkn]
              rw [mul_assoc, hp]
      _ = inverseBinomialPrefixReal26 (n + 1) * t ^ n := by
        rw [← Finset.sum_mul]
        rfl
  have hprod' :
      HasSum
        (fun n : ℕ =>
          inverseBinomialPrefixReal26 (n + 1) * t ^ n)
        ((1 - (1 - r) ^ 2 / r * Real.log (1 + r)) /
          (1 - t)) := by
    rw [← hterm]
    convert hprod using 1
    · rw [show (∑' m : ℕ,
          ((inverseBinomialA26 m : ℚ) : ℝ) * t ^ m) =
          1 - (1 - r) ^ 2 / r * Real.log (1 + r) by
            simpa [t] using hA.tsum_eq,
        hGeom.tsum_eq]
      ring
  simpa [t] using hprod'

/-- The nonnegative series term whose pointwise sum is
`nestedCyclotomicKernel26`.  It is the integration-by-parts transform of
`nestedBetaTerm26`. -/
def nestedCyclotomicSeriesTerm26 (n : ℕ) (r : ℝ) : ℝ :=
  -inverseBinomialPrefixReal26 (n + 1) *
      catalanParameter26 r ^ n *
      ((1 - r) / (1 + r) ^ 3) *
    Real.log r

theorem nestedCyclotomicSeriesTerm26_nonneg
    (n : ℕ) {r : ℝ} (hr0 : 0 ≤ r) (hr1 : r ≤ 1) :
    0 ≤ nestedCyclotomicSeriesTerm26 n r := by
  have hpref := inverseBinomialPrefixReal26_nonneg (n + 1)
  have ht : 0 ≤ catalanParameter26 r := by
    unfold catalanParameter26
    positivity
  have hlog : Real.log r ≤ 0 := Real.log_nonpos hr0 hr1
  have hratio : 0 ≤ (1 - r) / (1 + r) ^ 3 :=
    div_nonneg (sub_nonneg.mpr hr1) (by positivity)
  unfold nestedCyclotomicSeriesTerm26
  calc
    -inverseBinomialPrefixReal26 (n + 1) *
          catalanParameter26 r ^ n *
          ((1 - r) / (1 + r) ^ 3) *
        Real.log r =
      inverseBinomialPrefixReal26 (n + 1) *
          catalanParameter26 r ^ n *
          ((1 - r) / (1 + r) ^ 3) *
        (-Real.log r) := by ring
    _ ≥ 0 :=
      mul_nonneg
        (mul_nonneg
          (mul_nonneg hpref (pow_nonneg ht n)) hratio)
        (neg_nonneg.mpr hlog)

theorem catalanParameter26_hasDerivAt
    {r : ℝ} (hr : 1 + r ≠ 0) :
    HasDerivAt catalanParameter26
      ((1 - r) / (1 + r) ^ 3) r := by
  unfold catalanParameter26
  convert
    (hasDerivAt_id r).div
      ((hasDerivAt_const r 1).add (hasDerivAt_id r) |>.pow 2)
      (pow_ne_zero 2 hr) using 1
  all_goals
    simp
    field_simp [hr]
    ring

/-- The primitive factor
`S_(n+1) t(r)^(n+1)/(n+1)` used in the termwise integration by parts. -/
def nestedPrefixPower26 (n : ℕ) (r : ℝ) : ℝ :=
  inverseBinomialPrefixReal26 (n + 1) / (n + 1 : ℝ) *
    catalanParameter26 r ^ (n + 1)

theorem nestedPrefixPower26_hasDerivAt
    (n : ℕ) {r : ℝ} (hr : 1 + r ≠ 0) :
    HasDerivAt (nestedPrefixPower26 n)
      (inverseBinomialPrefixReal26 (n + 1) *
        catalanParameter26 r ^ n *
        ((1 - r) / (1 + r) ^ 3)) r := by
  have h :=
    ((catalanParameter26_hasDerivAt hr).pow (n + 1)).const_mul
      (inverseBinomialPrefixReal26 (n + 1) / (n + 1 : ℝ))
  convert h using 1
  push_cast
  have hn : (n : ℝ) + 1 ≠ 0 := by positivity
  field_simp [hn]

theorem nestedPrefixPower26_div
    (n : ℕ) {r : ℝ} (hr : r ≠ 0) :
    nestedPrefixPower26 n r / r = nestedBetaTerm26 n r := by
  unfold nestedPrefixPower26 nestedBetaTerm26 catalanParameter26
  rw [div_pow, pow_succ r n]
  have hpow :
      ((1 + r) ^ 2) ^ (n + 1) =
        (1 + r) ^ (2 * (n + 1)) := by
    rw [← pow_mul]
  field_simp [hr]
  rw [hpow]

theorem nestedBetaTerm26_continuousAt
    (n : ℕ) {r : ℝ} (hr : 1 + r ≠ 0) :
    ContinuousAt (nestedBetaTerm26 n) r := by
  unfold nestedBetaTerm26
  exact
    ((continuousAt_const.div_const _).mul (continuousAt_id.pow n)).div
      ((continuousAt_const.add continuousAt_id).pow (2 * (n + 1)))
      (pow_ne_zero _ hr)

/-- An antiderivative of the cyclotomic series term on `(0,1)`. -/
def nestedCyclotomicPrimitive26 (n : ℕ) (x : ℝ) : ℝ :=
  (∫ r : ℝ in 0..x, nestedBetaTerm26 n r) -
    nestedPrefixPower26 n x * Real.log x

theorem nestedCyclotomicPrimitive26_hasDerivAt
    (n : ℕ) {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasDerivAt (nestedCyclotomicPrimitive26 n)
      (nestedCyclotomicSeriesTerm26 n x) x := by
  have hxden : 1 + x ≠ 0 := by positivity
  have hcont := nestedBetaTerm26_continuousAt n hxden
  have hint :
      IntervalIntegrable (nestedBetaTerm26 n)
        MeasureTheory.volume 0 x := by
    apply ContinuousOn.intervalIntegrable
    intro y hy
    have hy' : y ∈ Set.Icc (0 : ℝ) x := by
      simpa [Set.uIcc_of_le hx0.le] using hy
    exact
      (nestedBetaTerm26_continuousAt n
        (by linarith [hy'.1])).continuousWithinAt
  have hI :
      HasDerivAt
        (fun u : ℝ => ∫ r : ℝ in 0..u, nestedBetaTerm26 n r)
        (nestedBetaTerm26 n x) x :=
    intervalIntegral.integral_hasDerivAt_right hint
      (ContinuousAt.stronglyMeasurableAtFilter isOpen_Ioi
        (fun y (hy : y ∈ Set.Ioi (-1 : ℝ)) =>
          nestedBetaTerm26_continuousAt n
            (by
              have hy' : (-1 : ℝ) < y := hy
              intro h
              linarith))
        x (show (-1 : ℝ) < x by linarith))
      hcont
  have hU := nestedPrefixPower26_hasDerivAt n hxden
  have hlog := Real.hasDerivAt_log (ne_of_gt hx0)
  have hP := hI.sub (hU.mul hlog)
  unfold nestedCyclotomicPrimitive26
  convert hP using 1
  rw [← nestedPrefixPower26_div n (ne_of_gt hx0)]
  unfold nestedCyclotomicSeriesTerm26
  ring

theorem nestedCyclotomicSeriesTerm26_intervalIntegrable
    (n : ℕ) :
    IntervalIntegrable (nestedCyclotomicSeriesTerm26 n)
      MeasureTheory.volume 0 1 := by
  have hcont :
      ContinuousOn
        (fun r : ℝ =>
          -inverseBinomialPrefixReal26 (n + 1) *
            catalanParameter26 r ^ n *
            ((1 - r) / (1 + r) ^ 3))
        (Set.Icc 0 1) := by
    intro r hr
    have hden : 1 + r ≠ 0 := by linarith [hr.1]
    unfold catalanParameter26
    exact
      (((continuousAt_const.mul
        ((continuousAt_id.div
          ((continuousAt_const.add continuousAt_id).pow 2)
          (pow_ne_zero 2 hden)).pow n))).mul
        ((continuousAt_const.sub continuousAt_id).div
          ((continuousAt_const.add continuousAt_id).pow 3)
          (pow_ne_zero 3 hden))).continuousWithinAt
  have hcont' :
      ContinuousOn
        (fun r : ℝ =>
          -inverseBinomialPrefixReal26 (n + 1) *
            catalanParameter26 r ^ n *
            ((1 - r) / (1 + r) ^ 3))
        (Set.uIcc 0 1) := by
    simpa [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hcont
  convert
    intervalIntegral.intervalIntegrable_log'.continuousOn_mul hcont' using 1

theorem nestedPrefixPower26_mul_log_tendsto_zero (n : ℕ) :
    Tendsto
      (fun x : ℝ => nestedPrefixPower26 n x * Real.log x)
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  have hpowlog :
      Tendsto (fun x : ℝ => Real.log x * x ^ (n + 1))
        (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    have h :=
      tendsto_log_mul_rpow_nhdsGT_zero
        (show (0 : ℝ) < (n : ℝ) + 1 by positivity)
    convert h using 1
    funext x
    rw [show (n : ℝ) + 1 = ((n + 1 : ℕ) : ℝ) by norm_num,
      Real.rpow_natCast]
  have hden :
      Tendsto
        (fun x : ℝ => ((1 + x) ^ (2 * (n + 1)))⁻¹)
        (𝓝[>] (0 : ℝ)) (𝓝 1) := by
    have hc :
        ContinuousAt
          (fun x : ℝ => ((1 + x) ^ (2 * (n + 1)))⁻¹) 0 := by
      exact
        ((continuousAt_const.add continuousAt_id).pow
          (2 * (n + 1))).inv₀ (by norm_num)
    simpa using tendsto_nhdsWithin_of_tendsto_nhds hc.tendsto
  have hmain :=
    (hpowlog.mul hden).const_mul
      (inverseBinomialPrefixReal26 (n + 1) / (n + 1 : ℝ))
  convert hmain using 1
  · funext x
    unfold nestedPrefixPower26 catalanParameter26
    rw [div_pow]
    have hp :
        ((1 + x) ^ 2) ^ (n + 1) =
          (1 + x) ^ (2 * (n + 1)) := by
      rw [← pow_mul]
    rw [hp]
    ring
  · ring

theorem nestedBetaIntegral26_tendsto_zero (n : ℕ) :
    Tendsto
      (fun x : ℝ => ∫ r : ℝ in 0..x, nestedBetaTerm26 n r)
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  have hcont : ContinuousAt (nestedBetaTerm26 n) 0 :=
    nestedBetaTerm26_continuousAt n (r := 0) (by norm_num)
  have hsm :
      StronglyMeasurableAtFilter (nestedBetaTerm26 n) (𝓝 (0 : ℝ))
        MeasureTheory.volume :=
    ContinuousAt.stronglyMeasurableAtFilter isOpen_Ioi
      (fun y (hy : y ∈ Set.Ioi (-1 : ℝ)) =>
        nestedBetaTerm26_continuousAt n
          (by
            have hy' : (-1 : ℝ) < y := hy
            intro h
            linarith))
      0 (by norm_num)
  have hderiv :
      HasDerivAt
        (fun x : ℝ => ∫ r : ℝ in 0..x, nestedBetaTerm26 n r)
        (nestedBetaTerm26 n 0) 0 :=
    intervalIntegral.integral_hasDerivAt_right
      (by simp) hsm hcont
  simpa using
    tendsto_nhdsWithin_of_tendsto_nhds hderiv.continuousAt.tendsto

theorem nestedBetaIntegral26_tendsto_one (n : ℕ) :
    Tendsto
      (fun x : ℝ => ∫ r : ℝ in 0..x, nestedBetaTerm26 n r)
      (𝓝[<] (1 : ℝ))
      (𝓝 (∫ r : ℝ in 0..1, nestedBetaTerm26 n r)) := by
  have hcont : ContinuousAt (nestedBetaTerm26 n) 1 :=
    nestedBetaTerm26_continuousAt n (r := 1) (by norm_num)
  have hsm :
      StronglyMeasurableAtFilter (nestedBetaTerm26 n) (𝓝 (1 : ℝ))
        MeasureTheory.volume :=
    ContinuousAt.stronglyMeasurableAtFilter isOpen_Ioi
      (fun y (hy : y ∈ Set.Ioi (-1 : ℝ)) =>
        nestedBetaTerm26_continuousAt n
          (by
            have hy' : (-1 : ℝ) < y := hy
            intro h
            linarith))
      1 (by norm_num)
  have hderiv :
      HasDerivAt
        (fun x : ℝ => ∫ r : ℝ in 0..x, nestedBetaTerm26 n r)
        (nestedBetaTerm26 n 1) 1 :=
    intervalIntegral.integral_hasDerivAt_right
      (nestedBetaTerm26_intervalIntegrable n) hsm hcont
  exact
    tendsto_nhdsWithin_of_tendsto_nhds hderiv.continuousAt.tendsto

theorem nestedPrefixPower26_mul_log_tendsto_one (n : ℕ) :
    Tendsto
      (fun x : ℝ => nestedPrefixPower26 n x * Real.log x)
      (𝓝[<] (1 : ℝ)) (𝓝 0) := by
  have hU : ContinuousAt (nestedPrefixPower26 n) 1 :=
    (nestedPrefixPower26_hasDerivAt n (r := 1) (by norm_num)).continuousAt
  have hlog : ContinuousAt Real.log 1 :=
    Real.continuousAt_log (by norm_num)
  have hprod := (hU.mul hlog).tendsto
  simpa using tendsto_nhdsWithin_of_tendsto_nhds hprod

theorem nestedCyclotomicPrimitive26_tendsto_zero (n : ℕ) :
    Tendsto (nestedCyclotomicPrimitive26 n)
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  unfold nestedCyclotomicPrimitive26
  simpa using
    (nestedBetaIntegral26_tendsto_zero n).sub
      (nestedPrefixPower26_mul_log_tendsto_zero n)

theorem nestedCyclotomicPrimitive26_tendsto_one (n : ℕ) :
    Tendsto (nestedCyclotomicPrimitive26 n)
      (𝓝[<] (1 : ℝ))
      (𝓝 (∫ r : ℝ in 0..1, nestedBetaTerm26 n r)) := by
  unfold nestedCyclotomicPrimitive26
  simpa using
    (nestedBetaIntegral26_tendsto_one n).sub
      (nestedPrefixPower26_mul_log_tendsto_one n)

/-- Each cyclotomic series term has exactly the same integral as the
corresponding nested inverse-binomial term. -/
theorem nestedCyclotomicSeriesTerm26_integral (n : ℕ) :
    (∫ r : ℝ in 0..1, nestedCyclotomicSeriesTerm26 n r) =
      ((inverseBinomialDTerm26 (n + 1) : ℚ) : ℝ) := by
  have hftc :=
    intervalIntegral.integral_eq_sub_of_hasDerivAt_of_tendsto
      (a := (0 : ℝ)) (b := 1)
      (f := nestedCyclotomicPrimitive26 n)
      (f' := nestedCyclotomicSeriesTerm26 n)
      (by norm_num)
      (fun x hx =>
        nestedCyclotomicPrimitive26_hasDerivAt n hx.1 hx.2)
      (nestedCyclotomicSeriesTerm26_intervalIntegrable n)
      (nestedCyclotomicPrimitive26_tendsto_zero n)
      (nestedCyclotomicPrimitive26_tendsto_one n)
  rw [hftc, sub_zero, nestedBetaTerm26_integral]

theorem nestedCyclotomicSeriesTerm26_hasSum_pointwise
    {r : ℝ} (hr0 : 0 < r) (hr1 : r < 1) :
    HasSum (fun n : ℕ => nestedCyclotomicSeriesTerm26 n r)
      (nestedCyclotomicKernel26 r) := by
  have hprefix :=
    inverseBinomialPrefixReal26_hasSum_catalan hr0 hr1
  have h :=
    hprefix.mul_left
      (-((1 - r) / (1 + r) ^ 3) * Real.log r)
  convert h using 1
  · funext n
    unfold nestedCyclotomicSeriesTerm26
    ring
  · have hrne : r ≠ 0 := ne_of_gt hr0
    have h1r : 1 + r ≠ 0 := by positivity
    have hq : 1 + r + r ^ 2 ≠ 0 :=
      (cyclotomicQuadratic_pos26 r).ne'
    unfold nestedCyclotomicKernel26 catalanParameter26
    have htden :
        1 - r / (1 + r) ^ 2 =
          (1 + r + r ^ 2) / (1 + r) ^ 2 := by
      field_simp [h1r]
      ring
    rw [htden]
    field_simp [hrne, h1r, hq]
    ring

theorem nestedCyclotomicSeriesTerm26_hasSum_pointwise_Ioc
    {r : ℝ} (hr0 : 0 < r) (hr1 : r ≤ 1) :
    HasSum (fun n : ℕ => nestedCyclotomicSeriesTerm26 n r)
      (nestedCyclotomicKernel26 r) := by
  by_cases hr : r = 1
  · subst r
    simp [nestedCyclotomicSeriesTerm26, nestedCyclotomicKernel26]
  · exact nestedCyclotomicSeriesTerm26_hasSum_pointwise
      hr0 (lt_of_le_of_ne hr1 hr)

theorem nestedCyclotomicSeriesTerm26_integral_norm (n : ℕ) :
    (∫ r : ℝ in 0..1, ‖nestedCyclotomicSeriesTerm26 n r‖) =
      ((inverseBinomialDTerm26 (n + 1) : ℚ) : ℝ) := by
  rw [← nestedCyclotomicSeriesTerm26_integral n]
  apply intervalIntegral.integral_congr
  intro r hr
  change
    ‖nestedCyclotomicSeriesTerm26 n r‖ =
      nestedCyclotomicSeriesTerm26 n r
  rw [Real.norm_eq_abs,
    abs_of_nonneg
      (nestedCyclotomicSeriesTerm26_nonneg n
        (by
          simpa [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hr.1)
        (by
          simpa [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hr.2))]

theorem nestedCyclotomicSeriesTerm26_integral_norm_summable :
    Summable
      (fun n : ℕ =>
        ∫ r : ℝ in 0..1,
          ‖nestedCyclotomicSeriesTerm26 n r‖) := by
  have hshift :
      Summable
        (fun n : ℕ =>
          ((inverseBinomialDTerm26 (n + 1) : ℚ) : ℝ)) := by
    simpa using
      (summable_nat_add_iff (f := fun k : ℕ =>
        ((inverseBinomialDTerm26 k : ℚ) : ℝ)) 1).2
        inverseBinomialDTerm26_summable
  exact hshift.congr fun n =>
    (nestedCyclotomicSeriesTerm26_integral_norm n).symm

/-- The nested inverse-binomial tail sums to the cyclotomic integral. -/
theorem inverseBinomialDTerm26_tail_hasSum_cyclotomic_integral :
    HasSum
      (fun n : ℕ => ((inverseBinomialDTerm26 (n + 1) : ℚ) : ℝ))
      (∫ r : ℝ in 0..1, nestedCyclotomicKernel26 r) := by
  have hInt :
      ∀ n : ℕ,
        MeasureTheory.Integrable (nestedCyclotomicSeriesTerm26 n)
          (MeasureTheory.volume.restrict (Set.Ioc 0 1)) := by
    intro n
    exact (nestedCyclotomicSeriesTerm26_intervalIntegrable n).1
  have hNorm :
      Summable
        (fun n : ℕ =>
          ∫ r : ℝ in Set.Ioc 0 1,
            ‖nestedCyclotomicSeriesTerm26 n r‖) := by
    simpa only [
      ← intervalIntegral.integral_of_le
        (by norm_num : (0 : ℝ) ≤ 1)] using
      nestedCyclotomicSeriesTerm26_integral_norm_summable
  have h :=
    MeasureTheory.hasSum_integral_of_summable_integral_norm
      (μ := MeasureTheory.volume.restrict (Set.Ioc (0 : ℝ) 1))
      hInt hNorm
  have h' :
      HasSum
        (fun n : ℕ => ((inverseBinomialDTerm26 (n + 1) : ℚ) : ℝ))
        (∫ r : ℝ in Set.Ioc 0 1,
          ∑' n : ℕ, nestedCyclotomicSeriesTerm26 n r) := by
    convert h using 1
    funext n
    rw [← intervalIntegral.integral_of_le
      (by norm_num : (0 : ℝ) ≤ 1)]
    exact (nestedCyclotomicSeriesTerm26_integral n).symm
  convert h' using 1
  rw [intervalIntegral.integral_of_le (by norm_num : (0 : ℝ) ≤ 1)]
  apply MeasureTheory.setIntegral_congr_fun measurableSet_Ioc
  intro r hr
  exact
    (nestedCyclotomicSeriesTerm26_hasSum_pointwise_Ioc
      hr.1 hr.2).tsum_eq.symm

/-- Exact analytic bridge for the new weight-three series in Problem 2.6:

`∑ k, Dₖ = ∫₀¹ nestedCyclotomicKernel26(r) dr`. -/
theorem inverseBinomialDTerm26_hasSum_cyclotomic_integral :
    HasSum
      (fun k : ℕ => ((inverseBinomialDTerm26 k : ℚ) : ℝ))
      (∫ r : ℝ in 0..1, nestedCyclotomicKernel26 r) := by
  have hfull :=
    (hasSum_nat_add_iff
      (f := fun k : ℕ => ((inverseBinomialDTerm26 k : ℚ) : ℝ)) 1).mp
      inverseBinomialDTerm26_tail_hasSum_cyclotomic_integral
  simpa [inverseBinomialDTerm26] using hfull

/-!
The following records the mathematical route realized by the compiled lemmas
above.  Every equality below is exact; no asymptotic majorant is needed.

For `k = n+1`, put `Sₖ = inverseBinomialPrefixReal26 k`.  The beta identity in
`Problem26` gives

```
Dₖ = Sₖ ∫₀¹ x^(k-1) (1-x)^k / k dx.
```

Pairing `x` with `1-x`, restricting to `[0,1/2]`, and substituting
`r = x/(1-x)` gives

```
Dₖ = ∫₀¹ Fₖ(r) dr,
Fₖ(r) = Sₖ r^(k-1) / (k (1+r)^(2k)).
```

The first continuation lemmas above have exactly these statements:

```
theorem nestedBetaTerm26_integral (n : ℕ) :
    (∫ r : ℝ in 0..1, nestedBetaTerm26 n r) =
      ((inverseBinomialDTerm26 (n + 1) : ℚ) : ℝ)

theorem nestedBetaTerm26_nonneg (n : ℕ) {r : ℝ}
    (hr0 : 0 ≤ r) :
    0 ≤ nestedBetaTerm26 n r

theorem nestedBetaTerm26_integral_norm (n : ℕ) :
    ∫ r : ℝ in 0..1, ‖nestedBetaTerm26 n r‖ =
      ((inverseBinomialDTerm26 (n + 1) : ℚ) : ℝ)

theorem nestedBetaTerm26_integral_norm_summable :
    Summable
      (fun n : ℕ =>
        ∫ r : ℝ in 0..1, ‖nestedBetaTerm26 n r‖)
```

The third lemma follows from nonnegativity and the first.  The fourth is
exactly `inverseBinomialDTerm26_summable`, shifted by one.  It is the
hypothesis required by
`MeasureTheory.hasSum_integral_of_summable_integral_norm`, so it supplies the
sum/integral interchange without a wasteful uniform bound.

For `0 < r < 1`, let `t = r/(1+r)^2`.  Directly from the definition of
`Fₖ`,

```
∑_{n≥0} F_{n+1}(r)
  = (1/r) ∑_{k≥1} Sₖ t^k/k
  = nestedGeneratingPrimitive26(t)/r.
```

The corresponding pointwise theorem is:

```
theorem nestedBetaTerm26_hasSum_pointwise
    {r : ℝ} (hr0 : 0 < r) (hr1 : r < 1) :
    HasSum (fun n : ℕ => nestedBetaTerm26 n r)
      (nestedGeneratingPrimitive26 (catalanParameter26 r) / r)
```

One convenient proof differentiates both power series: the finite-prefix
identity is

```
∑_{k≥0} Sₖ t^k = t * inverseBinomialAGF26 t / (1-t),
```

so the derivative of `∑_{k≥1} Sₖ t^k/k` is
`inverseBinomialAGF26 t/(1-t)`, and both sides vanish at zero.

The generating function already proved in this file gives, for `0 < s < 1`,

```
A(s/(1+s)^2) = 1 - (1-s)^2/s * log(1+s).
```

Since

```
d/ds (s/(1+s)^2) = (1-s)/(1+s)^3,
1 - s/(1+s)^2 = (1+s+s^2)/(1+s)^2,
```

the chain rule gives the exact primitive identity

```
theorem nestedGeneratingPrimitive26_catalan
    {r : ℝ} (hr0 : 0 < r) (hr1 : r < 1) :
    nestedGeneratingPrimitive26 (catalanParameter26 r) =
      ∫ s : ℝ in 0..r, nestedCatalanDerivative26 s
```

and hence

```
∑' k, ((inverseBinomialDTerm26 k : ℚ) : ℝ)
  = ∫ r in 0..1,
      (∫ s in 0..r, nestedCatalanDerivative26 s) / r.
```

Finally, Tonelli/Fubini on the triangle `0 < s < r < 1` yields

```
∫₀¹ (1/r) ∫₀ʳ B(s) ds dr
  = ∫₀¹ B(s) (∫ₛ¹ dr/r) ds
  = ∫₀¹ -B(s) log(s) ds,
```

where `B = nestedCatalanDerivative26`.  Pointwise,

```
-nestedCatalanDerivative26 s * Real.log s
  = (1-s) * ((1-s)^2 * Real.log (1+s) - s) * Real.log s
      / (s * (1+s) * (1+s+s^2))
  = nestedCyclotomicKernel26 s.
```

This last identity is purely `field_simp; ring` under `0 < s`.  The compiled
proof above avoids a separate triangular Fubini argument: it integrates each
nonnegative cyclotomic series term by parts, proves the endpoint limits at
zero and one, and then applies
`MeasureTheory.hasSum_integral_of_summable_integral_norm`.
-/

end RamanujanChallenge.P26

end
