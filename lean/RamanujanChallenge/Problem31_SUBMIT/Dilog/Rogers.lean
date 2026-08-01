import Ramanujan31.Dilog.Basic
import Ramanujan31.Dilog.RealBounds
import Mathlib.Analysis.Calculus.MeanValue
import Mathlib.Analysis.Calculus.SmoothSeries
import Mathlib.Analysis.SpecialFunctions.Log.NegMulLog

/-!
# The real Rogers dilogarithm

This file supplies the analytic API used by the direct proof of Problem 3.1.
The underlying `Real.dilog` is the power series constructed in
`Dilog/RealBounds.lean`.
-/

open Set
open Filter
open scoped Topology

namespace Real

/-- The Rogers dilogarithm on the real unit interval. -/
noncomputable def rogers (x : ℝ) : ℝ :=
  dilog x + (1 / 2 : ℝ) * log x * log (1 - x)

/-- Termwise derivative of the real dilogarithm series. -/
private noncomputable def dilogDerivTerm (x : ℝ) (n : ℕ) : ℝ :=
  x ^ n / ((n : ℝ) + 1)

private theorem hasDerivAt_dilogTerm (x : ℝ) (n : ℕ) :
    HasDerivAt (fun y : ℝ => dilogTerm y n) (dilogDerivTerm x n) x := by
  unfold dilogTerm dilogDerivTerm
  convert ((hasDerivAt_pow (n + 1) x).div_const (((n : ℝ) + 1) ^ 2)) using 1
  push_cast
  field_simp

@[simp] theorem dilog_zero : dilog 0 = 0 := by
  simp [dilog, dilogTerm]

/-- Summability of the defining series at the boundary point `1`. -/
theorem summable_dilogTerm_one : Summable (dilogTerm 1) := by
  have h : Summable (fun n : ℕ => (1 : ℝ) / ((n : ℝ) + 1) ^ 2) := by
    have hz := hasSum_zeta_two.summable
    have hs :=
      (summable_nat_add_iff
        (f := fun n : ℕ => (1 : ℝ) / (n : ℝ) ^ 2) 1).mpr hz
    refine hs.congr (fun n => ?_)
    push_cast
    ring
  refine h.congr (fun n => ?_)
  unfold dilogTerm
  norm_num

/-- The real Basel value `Li₂(1)=π²/6`. -/
theorem dilog_one : dilog 1 = Real.pi ^ 2 / 6 := by
  have hz : HasSum (fun n : ℕ => (1 : ℝ) / (n : ℝ) ^ 2)
      (Real.pi ^ 2 / 6 + ∑ i ∈ Finset.range 1, (1 : ℝ) / (i : ℝ) ^ 2) := by
    simpa using hasSum_zeta_two
  have hshift :=
    (hasSum_nat_add_iff
      (f := fun n : ℕ => (1 : ℝ) / (n : ℝ) ^ 2) 1).mpr hz
  have heq :
      (fun n : ℕ => (1 : ℝ) / ((n + 1 : ℕ) : ℝ) ^ 2) = dilogTerm 1 := by
    funext n
    unfold dilogTerm
    push_cast
    norm_num
  rw [heq] at hshift
  exact summable_dilogTerm_one.hasSum.unique hshift

@[simp] theorem rogers_one : rogers 1 = Real.pi ^ 2 / 6 := by
  simp [rogers, dilog_one]

/-- The real dilogarithm series is continuous on the closed unit interval,
including the boundary point `1`. -/
theorem continuousOn_dilog_Icc : ContinuousOn dilog (Icc 0 1) := by
  rw [show dilog = fun x : ℝ => ∑' n : ℕ, dilogTerm x n from rfl]
  refine continuousOn_tsum
    (u := fun n : ℕ => (1 : ℝ) / ((n : ℝ) + 1) ^ 2)
    (fun n => ?_) ?_ ?_
  · unfold dilogTerm
    fun_prop
  · have h : Summable (fun n : ℕ => (1 : ℝ) / ((n : ℝ) + 1) ^ 2) := by
      refine summable_dilogTerm_one.congr (fun n => ?_)
      unfold dilogTerm
      norm_num
    exact h
  · intro n x hx
    rcases hx with ⟨hx0, hx1⟩
    unfold dilogTerm
    rw [Real.norm_eq_abs, abs_div, abs_pow, abs_of_nonneg hx0,
      abs_of_nonneg (by positivity : (0 : ℝ) ≤ ((n : ℝ) + 1) ^ 2)]
    exact div_le_div_of_nonneg_right
      (pow_le_one₀ hx0 hx1) (sq_nonneg ((n : ℝ) + 1))

/-- The derivative of `Li₂(x)` on `(0,1)`. -/
theorem hasDerivAt_dilog {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasDerivAt dilog (-log (1 - x) / x) x := by
  let q : ℝ := (x + 1) / 2
  have hxq : x < q := by dsimp [q]; linarith
  have hq1 : q < 1 := by dsimp [q]; linarith
  have hq0 : 0 ≤ q := by dsimp [q]; linarith
  have hqabs : |q| < 1 := by rw [abs_of_nonneg hq0]; exact hq1
  have hsumq : Summable (fun n : ℕ => q ^ n) :=
    summable_geometric_of_lt_one hq0 hq1
  have hopen : IsOpen (Ioo (-q) q) := isOpen_Ioo
  have hconn : IsPreconnected (Ioo (-q) q) := isPreconnected_Ioo
  have hbound :
      ∀ n : ℕ, ∀ y : ℝ, y ∈ Ioo (-q) q →
        ‖dilogDerivTerm y n‖ ≤ q ^ n := by
    intro n y hy
    have hyabs : |y| < q := (abs_lt).2 hy
    unfold dilogDerivTerm
    rw [Real.norm_eq_abs, abs_div, abs_pow, abs_of_nonneg (by positivity : (0 : ℝ) ≤ (n : ℝ) + 1)]
    have hn : (0 : ℝ) ≤ (n : ℝ) := Nat.cast_nonneg n
    have hden : (1 : ℝ) ≤ (n : ℝ) + 1 := by linarith
    calc
      |y| ^ n / ((n : ℝ) + 1) ≤ |y| ^ n / 1 := by
        gcongr
      _ = |y| ^ n := by ring
      _ ≤ q ^ n := pow_le_pow_left₀ (abs_nonneg y) hyabs.le n
  have hzero : Summable (fun n : ℕ => dilogTerm 0 n) := by
    simp [dilogTerm]
  have hxmem : x ∈ Ioo (-q) q := by
    constructor <;> dsimp [q] <;> linarith
  have hseries :
      HasDerivAt (fun y : ℝ => ∑' n : ℕ, dilogTerm y n)
        (∑' n : ℕ, dilogDerivTerm x n) x :=
    hasDerivAt_tsum_of_isPreconnected hsumq hopen hconn
      (fun n y _ => hasDerivAt_dilogTerm y n) hbound
      (show (0 : ℝ) ∈ Ioo (-q) q by constructor <;> dsimp [q] <;> linarith)
      hzero hxmem
  have hxabs : |x| < 1 := by rw [abs_of_pos hx0]; exact hx1
  have hlog := (hasSum_pow_div_log_of_abs_lt_one hxabs).div_const x
  have hfun :
      (fun n : ℕ => x ^ (n + 1) / ((n : ℝ) + 1) / x) =
        fun n : ℕ => dilogDerivTerm x n := by
    funext n
    unfold dilogDerivTerm
    rw [pow_succ]
    field_simp [ne_of_gt hx0]
  have hderivsum :
      HasSum (fun n : ℕ => dilogDerivTerm x n) (-log (1 - x) / x) := by
    rw [hfun] at hlog
    exact hlog
  rw [show dilog = fun y : ℝ => ∑' n : ℕ, dilogTerm y n from rfl]
  convert hseries using 1
  exact hderivsum.tsum_eq.symm

/-- Differential of the Rogers dilogarithm on `(0,1)`. -/
theorem hasDerivAt_rogers {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasDerivAt rogers
      ((-log x / (1 - x) - log (1 - x) / x) / 2) x := by
  have hxne : x ≠ 0 := ne_of_gt hx0
  have h1x : 0 < 1 - x := by linarith
  have h1xne : 1 - x ≠ 0 := ne_of_gt h1x
  have hlogx := Real.hasDerivAt_log hxne
  have hlog1x :=
    ((hasDerivAt_const x 1).sub (hasDerivAt_id x)).log h1xne
  have hcorrection :=
    (hasDerivAt_const x (1 / 2 : ℝ)).mul (hlogx.mul hlog1x)
  have hcorrection' :
      HasDerivAt (fun y : ℝ => (1 / 2 : ℝ) * log y * log (1 - y))
        ((log (1 - x) / x - log x / (1 - x)) / 2) x := by
    convert hcorrection using 1
    · funext y
      simp
      ring
    · simp
      field_simp [hxne, h1xne]
      ring
  have hsum := (hasDerivAt_dilog hx0 hx1).add hcorrection'
  unfold rogers
  convert hsum using 1
  ring

/-- The logarithmic correction in the Rogers dilogarithm tends to zero at
the left endpoint of the unit interval. -/
theorem tendsto_log_mul_log_one_sub_zero :
    Tendsto (fun x : ℝ => log x * log (1 - x)) (𝓝[>] 0) (𝓝 0) := by
  have hmul :
      Tendsto (fun x : ℝ => x * log x) (𝓝[>] 0) (𝓝 0) := by
    simpa using
      Real.continuous_mul_log.continuousAt.tendsto.mono_left
        (inf_le_left : 𝓝[>] (0 : ℝ) ≤ 𝓝 0)
  have hlog :
      HasDerivAt (fun x : ℝ => log (1 - x)) (-1) 0 := by
    have hinner :=
      (hasDerivAt_const (0 : ℝ) (1 : ℝ)).sub (hasDerivAt_id 0)
    have hraw :=
      hinner.log (by norm_num)
    convert hraw using 1
    norm_num
  have hratio :
      Tendsto (fun x : ℝ => log (1 - x) / x) (𝓝[>] 0) (𝓝 (-1)) := by
    simpa [div_eq_mul_inv, mul_comm] using hlog.tendsto_slope_zero_right
  have hprod :
      Tendsto
        (fun x : ℝ => (x * log x) * (log (1 - x) / x))
        (𝓝[>] 0) (𝓝 0) := by
    simpa using hmul.mul hratio
  refine hprod.congr' ?_
  filter_upwards [self_mem_nhdsWithin] with x hx
  have hxne : x ≠ 0 := ne_of_gt hx
  field_simp

private noncomputable def eulerAux (x : ℝ) : ℝ :=
  dilog x + dilog (1 - x) + log x * log (1 - x)

private theorem hasDerivAt_eulerAux {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasDerivAt eulerAux 0 x := by
  have hxne : x ≠ 0 := ne_of_gt hx0
  have h1x : 0 < 1 - x := by linarith
  have h1xne : 1 - x ≠ 0 := ne_of_gt h1x
  have hd₁ := hasDerivAt_dilog hx0 hx1
  have hd₂ :=
    (hasDerivAt_dilog h1x (by linarith : 1 - x < 1)).comp x
      ((hasDerivAt_const x 1).sub (hasDerivAt_id x))
  have hd₂' :
      HasDerivAt (fun y : ℝ => dilog (1 - y)) (log x / (1 - x)) x := by
    convert hd₂ using 1
    rw [show 1 - (1 - x) = x by ring]
    ring
  have hlogRaw :=
    (Real.hasDerivAt_log hxne).mul
      (((hasDerivAt_const x 1).sub (hasDerivAt_id x)).log h1xne)
  have hlog :
      HasDerivAt (fun y : ℝ => log y * log (1 - y))
        (log (1 - x) / x - log x / (1 - x)) x := by
    convert hlogRaw using 1
    simp
    field_simp [hxne, h1xne]
    ring
  unfold eulerAux
  convert (hd₁.add hd₂').add hlog using 1
  ring

private theorem tendsto_eulerAux_zero :
    Tendsto eulerAux (𝓝[>] 0) (𝓝 (Real.pi ^ 2 / 6)) := by
  have hfilter0 :
      𝓝[>] (0 : ℝ) ≤ 𝓝[Icc 0 1] (0 : ℝ) :=
    nhdsWithin_le_iff.mpr (Icc_mem_nhdsGT zero_lt_one)
  have hdilog0 :
      Tendsto dilog (𝓝[>] 0) (𝓝 0) := by
    simpa [dilog_zero] using
      (continuousOn_dilog_Icc.continuousWithinAt
        (by norm_num : (0 : ℝ) ∈ Icc 0 1)).mono_left hfilter0
  have hsub :
      Tendsto (fun x : ℝ => 1 - x) (𝓝[>] 0) (𝓝[Icc 0 1] 1) := by
    refine tendsto_nhdsWithin_iff.mpr ⟨?_, ?_⟩
    · have hto0 :
          Tendsto (fun x : ℝ => x) (𝓝[>] 0) (𝓝 0) :=
        (tendsto_id :
          Tendsto (fun x : ℝ => x) (𝓝 0) (𝓝 0)).mono_left inf_le_left
      simpa using tendsto_const_nhds.sub hto0
    · filter_upwards
        [self_mem_nhdsWithin,
          mem_inf_of_left (Iio_mem_nhds (show (0 : ℝ) < 1 by norm_num))]
        with x hx0 hx1
      change 0 < x at hx0
      change x < 1 at hx1
      constructor <;> linarith
  have hdilog1 :
      Tendsto (fun x : ℝ => dilog (1 - x)) (𝓝[>] 0)
        (𝓝 (Real.pi ^ 2 / 6)) := by
    have hdilogAtOne :
        Tendsto dilog (𝓝[Icc 0 1] (1 : ℝ)) (𝓝 (dilog 1)) :=
      continuousOn_dilog_Icc.continuousWithinAt
        (by norm_num : (1 : ℝ) ∈ Icc 0 1)
    simpa [dilog_one] using hdilogAtOne.comp hsub
  unfold eulerAux
  simpa using
    (hdilog0.add hdilog1).add tendsto_log_mul_log_one_sub_zero

/-- Euler's reflection formula for the real dilogarithm on `(0,1)`. -/
theorem dilog_add_dilog_one_sub {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    dilog x + dilog (1 - x) + log x * log (1 - x) =
      Real.pi ^ 2 / 6 := by
  have hdiff : DifferentiableOn ℝ eulerAux (Ioo 0 1) := by
    intro y hy
    exact (hasDerivAt_eulerAux hy.1 hy.2).differentiableAt.differentiableWithinAt
  have hzero : (Ioo (0 : ℝ) 1).EqOn (deriv eulerAux) 0 := by
    intro y hy
    exact (hasDerivAt_eulerAux hy.1 hy.2).deriv
  have hconst :
      ∀ ⦃y : ℝ⦄, y ∈ Ioo 0 1 → eulerAux x = eulerAux y := by
    intro y hy
    exact isOpen_Ioo.is_const_of_deriv_eq_zero
      isPreconnected_Ioo hdiff hzero ⟨hx0, hx1⟩ hy
  have hevent :
      eulerAux =ᶠ[𝓝[>] (0 : ℝ)] fun _ => eulerAux x := by
    filter_upwards
      [self_mem_nhdsWithin,
        mem_inf_of_left (Iio_mem_nhds (show (0 : ℝ) < 1 by norm_num))]
      with y hy0 hy1
    exact (hconst ⟨hy0, hy1⟩).symm
  have hlim :
      Tendsto eulerAux (𝓝[>] 0) (𝓝 (eulerAux x)) :=
    tendsto_const_nhds.congr' hevent.symm
  have := tendsto_nhds_unique hlim tendsto_eulerAux_zero
  simpa [eulerAux] using this

/-- Euler's reflection formula in Rogers normalization. -/
theorem rogers_add_rogers_one_sub {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    rogers x + rogers (1 - x) = Real.pi ^ 2 / 6 := by
  rw [show rogers x + rogers (1 - x) =
      dilog x + dilog (1 - x) + log x * log (1 - x) by
    unfold rogers
    ring]
  exact dilog_add_dilog_one_sub hx0 hx1

end Real
