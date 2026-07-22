/-
  Dilogarithm: Li₂(z) = Σ_{n≥1} z^n/n² for |z| ≤ 1.

  Foundation for Problem 3.1 (knot regulator).
  NOT in Mathlib v4.30.0.
-/
import Mathlib.Analysis.SpecificLimits.Basic
import Mathlib.Topology.Algebra.InfiniteSum.Basic
import Mathlib.Topology.Algebra.InfiniteSum.NatInt
import Mathlib.Analysis.Normed.Group.InfiniteSum
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.NumberTheory.ZetaValues
import Mathlib.Analysis.Calculus.SmoothSeries
import Mathlib.Analysis.Calculus.MeanValue
import Mathlib.Analysis.SpecialFunctions.Log.NegMulLog

open Filter Topology Set

noncomputable section

/-! ## Definition -/

def dilog (z : ℝ) : ℝ := ∑' n : ℕ, z ^ (n + 1) / (↑(n + 1) : ℝ) ^ 2

/-! ## Auxiliary: summability of 1/(n+1)² -/

private theorem summable_one_div_succ_sq :
    Summable (fun n : ℕ => (1 : ℝ) / (↑(n + 1) : ℝ) ^ 2) := by
  have h := hasSum_zeta_two.summable
  have h1 : Summable (fun n : ℕ => (fun m : ℕ => (1 : ℝ) / (↑m : ℝ) ^ 2) (n + 1)) :=
    (summable_nat_add_iff 1).mpr h
  exact h1.congr (fun n => by simp)

/-! ## Basic values -/

theorem dilog_zero : dilog 0 = 0 := by
  simp [dilog, zero_pow (by omega : _ + 1 ≠ 0)]

theorem dilog_summable {z : ℝ} (hz : |z| ≤ 1) :
    Summable (fun n : ℕ => z ^ (n + 1) / (↑(n + 1) : ℝ) ^ 2) := by
  refine Summable.of_norm_bounded
    (g := fun n => (1 : ℝ) / (↑(n + 1) : ℝ) ^ 2) summable_one_div_succ_sq ?_
  intro n
  rw [Real.norm_eq_abs]
  rw [abs_div, abs_pow]
  have hpos : (0 : ℝ) ≤ (↑(n + 1) : ℝ) ^ 2 := by positivity
  rw [abs_of_nonneg hpos]
  apply div_le_div_of_nonneg_right _ hpos
  calc |z| ^ (n + 1) ≤ 1 ^ (n + 1) := pow_le_pow_left₀ (abs_nonneg z) hz (n + 1)
    _ = 1 := one_pow _

theorem dilog_one : dilog 1 = Real.pi ^ 2 / 6 := by
  unfold dilog
  have heq : ∀ n : ℕ, (1 : ℝ) ^ (n + 1) / (↑(n + 1) : ℝ) ^ 2 =
      (fun m : ℕ => (1 : ℝ) / (↑m : ℝ) ^ 2) (n + 1) := fun n => by simp
  simp_rw [heq]
  have hz := hasSum_zeta_two
  have hsumm := hz.summable
  have htsum := hz.tsum_eq
  have hshift := hsumm.tsum_eq_zero_add
  simp only [Nat.cast_zero, zero_pow (two_ne_zero), div_zero, zero_add] at hshift
  linarith

/-! ## Rogers dilogarithm -/

def rogersDialogarithm (z : ℝ) : ℝ :=
  dilog z + (1 / 2) * Real.log z * Real.log (1 - z)

/-! ## Extended Rogers (with lift correction) -/

def rogersExtended (z : ℝ) (p q : ℤ) : ℝ :=
  rogersDialogarithm z + Real.pi ^ 2 / 2 * (↑p * ↑q)

/-! ## Key functional equations -/

private theorem dilog_hasDerivAt {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasDerivAt dilog (-(Real.log (1 - x)) / x) x := by
  have hx : |x| < 1 := by rw [abs_of_pos hx0]; exact hx1
  let r : ℝ := (|x| + 1) / 2
  have hr0 : 0 ≤ r := by dsimp [r]; positivity
  have hrpos : 0 < r := by dsimp [r]; positivity
  have hr1 : r < 1 := by dsimp [r]; linarith
  have hxr : |x| < r := by dsimp [r]; linarith
  have hu : Summable (fun n : ℕ => r ^ n) :=
    summable_geometric_of_lt_one hr0 hr1
  have hterm : ∀ n : ℕ, ∀ y : ℝ, y ∈ Ioo (-r) r →
      HasDerivAt (fun y : ℝ => y ^ (n + 1) / (↑(n + 1) : ℝ) ^ 2)
        (y ^ n / (↑(n + 1) : ℝ)) y := by
    intro n y hy
    convert (hasDerivAt_pow (n + 1) y).div_const ((↑(n + 1) : ℝ) ^ 2) using 1
    · simp
      field_simp
  have hbound : ∀ n : ℕ, ∀ y : ℝ, y ∈ Ioo (-r) r →
      ‖y ^ n / (↑(n + 1) : ℝ)‖ ≤ r ^ n := by
    intro n y hy
    rw [Real.norm_eq_abs, abs_div, abs_pow]
    have hyr : |y| < r := (abs_lt).2 hy
    calc
      |y| ^ n / |(↑(n + 1) : ℝ)| ≤ |y| ^ n / 1 := by
        gcongr
        rw [Nat.cast_add, Nat.cast_one, abs_of_nonneg (by positivity)]
        norm_num
      _ = |y| ^ n := by ring
      _ ≤ r ^ n := pow_le_pow_left₀ (abs_nonneg y) hyr.le n
  have hzero : Summable (fun n : ℕ =>
      (0 : ℝ) ^ (n + 1) / (↑(n + 1) : ℝ) ^ 2) := by
    simp [zero_pow (Nat.succ_ne_zero _)]
  have hxmem : x ∈ Ioo (-r) r := (abs_lt).1 hxr
  have hd : HasDerivAt
      (fun y : ℝ => ∑' n : ℕ, y ^ (n + 1) / (↑(n + 1) : ℝ) ^ 2)
      (∑' n : ℕ, x ^ n / (↑(n + 1) : ℝ)) x := by
    exact hasDerivAt_tsum_of_isPreconnected hu isOpen_Ioo
      (convex_Ioo (-r) r).isPreconnected hterm hbound
      (show (0 : ℝ) ∈ Ioo (-r) r by constructor <;> linarith)
      hzero hxmem
  rw [show dilog = (fun y : ℝ => ∑' n : ℕ,
      y ^ (n + 1) / (↑(n + 1) : ℝ) ^ 2) from rfl]
  convert hd using 1
  have hxne : x ≠ 0 := ne_of_gt hx0
  have hsum := Real.hasSum_pow_div_log_of_abs_lt_one hx
  have hsum' := hsum.mul_left (1 / x)
  have hsum'' : HasSum (fun n : ℕ => x ^ n / (↑(n + 1) : ℝ))
      (-(Real.log (1 - x)) / x) := by
    convert hsum' using 1
    · funext n
      simp only [Nat.cast_add, Nat.cast_one, pow_succ]
      field_simp
    · field_simp
  exact hsum''.tsum_eq.symm

private theorem dilog_continuousOn_unit :
    ContinuousOn dilog (Icc (-1 : ℝ) 1) := by
  unfold dilog
  refine continuousOn_tsum (u := fun n : ℕ => (1 : ℝ) / (↑(n + 1) : ℝ) ^ 2)
    (fun n => ?_) summable_one_div_succ_sq ?_
  · fun_prop
  · intro n x hx
    rw [Real.norm_eq_abs, abs_div, abs_pow]
    have habs : |x| ≤ 1 := (abs_le).2 hx
    have hpos : (0 : ℝ) ≤ (↑(n + 1) : ℝ) ^ 2 := by positivity
    rw [abs_of_nonneg hpos]
    apply div_le_div_of_nonneg_right _ hpos
    calc
      |x| ^ (n + 1) ≤ 1 ^ (n + 1) :=
        pow_le_pow_left₀ (abs_nonneg x) habs (n + 1)
      _ = 1 := one_pow _

private theorem log_mul_log_one_sub_continuousAt_one :
    ContinuousAt (fun x : ℝ => Real.log x * Real.log (1 - x)) 1 := by
  have heq : (fun x : ℝ => Real.log x * Real.log (1 - x)) =
      (fun x : ℝ => -(dslope Real.log 1 x) *
        ((1 - x) * Real.log (1 - x))) := by
    funext x
    by_cases hx : x = 1
    · subst x
      simp
    · rw [dslope_of_ne Real.log hx]
      rw [slope_def_field, Real.log_one]
      field_simp
      ring_nf
  rw [heq]
  apply ContinuousAt.mul
  · exact (continuousAt_dslope_same.mpr
      (Real.differentiableAt_log one_ne_zero)).neg
  · exact Real.continuous_mul_log.continuousAt.comp
      (continuousAt_const.sub continuousAt_id)

private def reflectionAux (x : ℝ) : ℝ :=
  dilog x + dilog (1 - x) + Real.log x * Real.log (1 - x)

private theorem reflectionAux_hasDerivAt_zero {x : ℝ}
    (hx0 : 0 < x) (hx1 : x < 1) :
    HasDerivAt reflectionAux 0 x := by
  have hxne : x ≠ 0 := ne_of_gt hx0
  have h1x0 : 0 < 1 - x := sub_pos.mpr hx1
  have h1x1 : 1 - x < 1 := by linarith
  have h1xne : 1 - x ≠ 0 := ne_of_gt h1x0
  have hd₁ := dilog_hasDerivAt hx0 hx1
  have hd₂ := (dilog_hasDerivAt h1x0 h1x1).comp_const_sub 1 x
  have hl₁ : HasDerivAt (fun y : ℝ => Real.log y) (1 / x) x := by
    simpa [one_div] using Real.hasDerivAt_log hxne
  have hsub : HasDerivAt (fun y : ℝ => 1 - y) (-1) x :=
    (hasDerivAt_id x).const_sub 1
  have hl₂ : HasDerivAt (fun y : ℝ => Real.log (1 - y)) (-1 / (1 - x)) x :=
    hsub.log h1xne
  unfold reflectionAux
  convert (hd₁.add hd₂).add (hl₁.mul hl₂) using 1
  · field_simp
    ring_nf

private theorem reflectionAux_continuousOn {a : ℝ}
    (ha0 : 0 < a) (ha1 : a < 1) :
    ContinuousOn reflectionAux (Icc a 1) := by
  have hda : ContinuousOn dilog (Icc a 1) :=
    dilog_continuousOn_unit.mono (by
      intro x hx
      constructor
      · linarith [hx.1]
      · exact hx.2)
  have hd1a : ContinuousOn (fun x : ℝ => dilog (1 - x)) (Icc a 1) := by
    apply dilog_continuousOn_unit.comp
    · fun_prop
    · intro x hx
      constructor <;> linarith [hx.1, hx.2]
  have hlogprod : ContinuousOn
      (fun x : ℝ => Real.log x * Real.log (1 - x)) (Icc a 1) := by
    intro x hx
    by_cases hx1 : x = 1
    · subst x
      exact log_mul_log_one_sub_continuousAt_one.continuousWithinAt
    · have hx0 : x ≠ 0 := by linarith [hx.1]
      have h1x0 : 1 - x ≠ 0 := sub_ne_zero.mpr (Ne.symm hx1)
      exact ((Real.continuousAt_log hx0).mul
        ((Real.continuousAt_log h1x0).comp
          (continuousAt_const.sub continuousAt_id))).continuousWithinAt
  exact (hda.add hd1a).add hlogprod

theorem dilog_reflection {z : ℝ} (hz0 : 0 < z) (hz1 : z < 1) :
    dilog z + dilog (1 - z) = Real.pi ^ 2 / 6 - Real.log z * Real.log (1 - z) := by
  have hconst := constant_of_has_deriv_right_zero
    (reflectionAux_continuousOn hz0 hz1)
    (fun x hx => (reflectionAux_hasDerivAt_zero
      (lt_of_lt_of_le hz0 hx.1) hx.2).hasDerivWithinAt)
  have h := hconst 1 (right_mem_Icc.mpr hz1.le)
  simp [reflectionAux, dilog_one, dilog_zero] at h
  linarith

/-! ## Extended Rogers dilogarithm

The series-based `dilog` only converges for |z| ≤ 1, and `Real.log` returns 0
for non-positive arguments. So `rogersDialogarithm` only gives the correct
mathematical value for z ∈ (0,1).

For the regulator computation in P3.1, we need Rogers values at z > 1 and z < 0.
We define these via the functional equations (as definitions, not theorems):

  For z > 1: R(z) = π²/3 − R(1/z)     [inversion]
  For z < 0: R(z) = −R(z/(z−1)) − ½·ln(1−z)·ln(−z)  [Landen]
-/

def rogersGtOne (z : ℝ) : ℝ :=
  Real.pi ^ 2 / 3 - rogersDialogarithm (1 / z)

def rogersNeg (z : ℝ) : ℝ :=
  -rogersDialogarithm (z / (z - 1)) -
  (1 / 2) * Real.log (1 - z) * Real.log (-z)

/-! ## Rogers dilogarithm at special values -/

theorem rogers_zero : rogersDialogarithm 0 = 0 := by
  unfold rogersDialogarithm
  simp [dilog_zero, Real.log_zero]

theorem rogers_one : rogersDialogarithm 1 = Real.pi ^ 2 / 6 := by
  unfold rogersDialogarithm
  simp [dilog_one, Real.log_one]

/-! ## Five-term relation (Abel's identity) -/

private theorem rogersDialogarithm_hasDerivAt {z : ℝ}
    (hz0 : 0 < z) (hz1 : z < 1) :
    HasDerivAt rogersDialogarithm
      (-(1 / 2) *
        (Real.log (1 - z) / z + Real.log z / (1 - z))) z := by
  have hz_ne : z ≠ 0 := ne_of_gt hz0
  have h1z0 : 0 < 1 - z := sub_pos.mpr hz1
  have h1z_ne : 1 - z ≠ 0 := ne_of_gt h1z0
  have hd := dilog_hasDerivAt hz0 hz1
  have hlz : HasDerivAt (fun w : ℝ => Real.log w) (1 / z) z := by
    simpa [one_div] using Real.hasDerivAt_log hz_ne
  have hsub : HasDerivAt (fun w : ℝ => 1 - w) (-1) z :=
    (hasDerivAt_id z).const_sub 1
  have hl1z : HasDerivAt (fun w : ℝ => Real.log (1 - w))
      (-1 / (1 - z)) z := hsub.log h1z_ne
  unfold rogersDialogarithm
  convert hd.add ((hlz.mul hl1z).const_mul (1 / 2)) using 1
  · funext w
    simp only [Pi.add_apply, Pi.mul_apply]
    ring
  · field_simp
    ring

private theorem log_mul_log_one_sub_continuousAt_zero :
    ContinuousAt (fun z : ℝ => Real.log z * Real.log (1 - z)) 0 := by
  have hinner : ContinuousAt (fun z : ℝ => 1 - z) 0 := by fun_prop
  have hcomp := log_mul_log_one_sub_continuousAt_one.comp_of_eq hinner (by norm_num)
  convert hcomp using 1
  funext z
  simp only [Function.comp_apply, sub_sub_cancel]
  ring

private theorem log_mul_log_one_sub_continuousOn_unit :
    ContinuousOn (fun z : ℝ => Real.log z * Real.log (1 - z))
      (Icc (0 : ℝ) 1) := by
  intro z hz
  rcases eq_or_lt_of_le hz.1 with h0 | hz0
  · subst z
    exact log_mul_log_one_sub_continuousAt_zero.continuousWithinAt
  rcases eq_or_lt_of_le hz.2 with h1 | hz1
  · subst z
    exact log_mul_log_one_sub_continuousAt_one.continuousWithinAt
  · have hz_ne : z ≠ 0 := ne_of_gt hz0
    have h1z_ne : 1 - z ≠ 0 := ne_of_gt (sub_pos.mpr hz1)
    exact ((Real.continuousAt_log hz_ne).mul
      ((Real.continuousAt_log h1z_ne).comp
        (continuousAt_const.sub continuousAt_id))).continuousWithinAt

private theorem rogersDialogarithm_continuousOn_unit :
    ContinuousOn rogersDialogarithm (Icc (0 : ℝ) 1) := by
  have hd : ContinuousOn dilog (Icc (0 : ℝ) 1) :=
    dilog_continuousOn_unit.mono (by
      intro z hz
      constructor
      · linarith [hz.1]
      · exact hz.2)
  have hp : ContinuousOn
      (fun z : ℝ => (1 / 2) *
        (Real.log z * Real.log (1 - z))) (Icc (0 : ℝ) 1) :=
    continuousOn_const.mul log_mul_log_one_sub_continuousOn_unit
  unfold rogersDialogarithm
  simpa only [mul_assoc] using hd.add hp

private theorem rogersFiveAux_hasDerivAt_zero {x y : ℝ}
    (hx0 : 0 < x) (hx1 : x < 1) (hy0 : 0 < y) (hy1 : y < 1) :
    HasDerivAt
      (fun t : ℝ =>
        rogersDialogarithm t + rogersDialogarithm y -
          rogersDialogarithm (t * y) -
          rogersDialogarithm (t * (1 - y) / (1 - t * y)) -
          rogersDialogarithm (y * (1 - t) / (1 - t * y))) 0 x := by
  have hx_ne : x ≠ 0 := ne_of_gt hx0
  have hy_ne : y ≠ 0 := ne_of_gt hy0
  have h1x0 : 0 < 1 - x := sub_pos.mpr hx1
  have h1y0 : 0 < 1 - y := sub_pos.mpr hy1
  have h1x_ne : 1 - x ≠ 0 := ne_of_gt h1x0
  have h1y_ne : 1 - y ≠ 0 := ne_of_gt h1y0
  have hxy0 : 0 < x * y := mul_pos hx0 hy0
  have hxy1 : x * y < 1 := by
    calc
      x * y < 1 * y := mul_lt_mul_of_pos_right hx1 hy0
      _ = y := one_mul y
      _ < 1 := hy1
  have hd0 : 0 < 1 - x * y := sub_pos.mpr hxy1
  have hd_ne : 1 - x * y ≠ 0 := ne_of_gt hd0
  have hd_ne' : 1 - y * x ≠ 0 := by
    simpa [mul_comm] using hd_ne
  have hb_num_lt : x * (1 - y) < 1 - x * y := by
    nlinarith
  have hc_num_lt : y * (1 - x) < 1 - x * y := by
    nlinarith
  have hb0 : 0 < x * (1 - y) / (1 - x * y) :=
    div_pos (mul_pos hx0 h1y0) hd0
  have hb1 : x * (1 - y) / (1 - x * y) < 1 :=
    (div_lt_one hd0).2 hb_num_lt
  have hc0 : 0 < y * (1 - x) / (1 - x * y) :=
    div_pos (mul_pos hy0 h1x0) hd0
  have hc1 : y * (1 - x) / (1 - x * y) < 1 :=
    (div_lt_one hd0).2 hc_num_lt
  have hden : HasDerivAt (fun t : ℝ => 1 - t * y) (-y) x := by
    simpa using ((hasDerivAt_id x).mul_const y).const_sub 1
  have hb_inner : HasDerivAt
      (fun t : ℝ => t * (1 - y) / (1 - t * y))
      ((1 - y) / (1 - x * y) ^ 2) x := by
    convert ((hasDerivAt_id x).mul_const (1 - y)).div hden hd_ne using 1
    simp only [id_eq]
    field_simp [hd_ne]
    ring
  have hc_num : HasDerivAt (fun t : ℝ => y * (1 - t)) (-y) x := by
    simpa using ((hasDerivAt_id x).const_sub 1).const_mul y
  have hc_inner : HasDerivAt
      (fun t : ℝ => y * (1 - t) / (1 - t * y))
      (-y * (1 - y) / (1 - x * y) ^ 2) x := by
    convert hc_num.div hden hd_ne using 1
    field_simp [hd_ne]
    ring
  have hdx := rogersDialogarithm_hasDerivAt hx0 hx1
  have hdy : HasDerivAt (fun _ : ℝ => rogersDialogarithm y) 0 x :=
    hasDerivAt_const x _
  have hdxy := (rogersDialogarithm_hasDerivAt hxy0 hxy1).comp x
    ((hasDerivAt_id x).mul_const y)
  have hdb := (rogersDialogarithm_hasDerivAt hb0 hb1).comp x hb_inner
  have hdc := (rogersDialogarithm_hasDerivAt hc0 hc1).comp x hc_inner
  have hlogxy :
      Real.log (x * y) = Real.log x + Real.log y :=
    Real.log_mul hx_ne hy_ne
  have honeb :
      1 - x * (1 - y) / (1 - x * y) =
        (1 - x) / (1 - x * y) := by
    field_simp [hd_ne, hd_ne']
    ring
  have honec :
      1 - y * (1 - x) / (1 - x * y) =
        (1 - y) / (1 - x * y) := by
    field_simp [hd_ne, hd_ne']
    ring
  have hlogb :
      Real.log (x * (1 - y) / (1 - x * y)) =
        Real.log x + Real.log (1 - y) - Real.log (1 - x * y) := by
    rw [Real.log_div (mul_ne_zero hx_ne h1y_ne) hd_ne,
      Real.log_mul hx_ne h1y_ne]
  have hlogc :
      Real.log (y * (1 - x) / (1 - x * y)) =
        Real.log y + Real.log (1 - x) - Real.log (1 - x * y) := by
    rw [Real.log_div (mul_ne_zero hy_ne h1x_ne) hd_ne,
      Real.log_mul hy_ne h1x_ne]
  have hlogoneb :
      Real.log (1 - x * (1 - y) / (1 - x * y)) =
        Real.log (1 - x) - Real.log (1 - x * y) := by
    rw [honeb, Real.log_div h1x_ne hd_ne]
  have hlogonec :
      Real.log (1 - y * (1 - x) / (1 - x * y)) =
        Real.log (1 - y) - Real.log (1 - x * y) := by
    rw [honec, Real.log_div h1y_ne hd_ne]
  convert (((hdx.add hdy).sub hdxy).sub hdb).sub hdc using 1
  · rw [hlogxy, hlogb, hlogc, hlogoneb, hlogonec, honeb, honec]
    field_simp
    ring

private theorem rogersFiveAux_continuousOn {a y : ℝ}
    (ha0 : 0 < a) (hy0 : 0 < y) (hy1 : y < 1) :
    ContinuousOn
      (fun t : ℝ =>
        rogersDialogarithm t + rogersDialogarithm y -
          rogersDialogarithm (t * y) -
          rogersDialogarithm (t * (1 - y) / (1 - t * y)) -
          rogersDialogarithm (y * (1 - t) / (1 - t * y)))
      (Icc a 1) := by
  have hden_pos : ∀ t ∈ Icc a 1, 0 < 1 - t * y := by
    intro t ht
    have hty : t * y ≤ y := by
      calc
        t * y ≤ 1 * y := mul_le_mul_of_nonneg_right ht.2 hy0.le
        _ = y := one_mul y
    linarith
  have hden_ne : ∀ t ∈ Icc a 1, 1 - t * y ≠ 0 :=
    fun t ht => ne_of_gt (hden_pos t ht)
  have hb_cont : ContinuousOn
      (fun t : ℝ => t * (1 - y) / (1 - t * y)) (Icc a 1) := by
    apply ContinuousOn.div
    · fun_prop
    · fun_prop
    · exact hden_ne
  have hc_cont : ContinuousOn
      (fun t : ℝ => y * (1 - t) / (1 - t * y)) (Icc a 1) := by
    apply ContinuousOn.div
    · fun_prop
    · fun_prop
    · exact hden_ne
  have hmap_x : MapsTo (fun t : ℝ => t) (Icc a 1) (Icc (0 : ℝ) 1) := by
    intro t ht
    exact ⟨le_trans ha0.le ht.1, ht.2⟩
  have hmap_xy : MapsTo (fun t : ℝ => t * y)
      (Icc a 1) (Icc (0 : ℝ) 1) := by
    intro t ht
    constructor
    · exact mul_nonneg (le_trans ha0.le ht.1) hy0.le
    · calc
        t * y ≤ 1 * y := mul_le_mul_of_nonneg_right ht.2 hy0.le
        _ = y := one_mul y
        _ ≤ 1 := hy1.le
  have hmap_b : MapsTo
      (fun t : ℝ => t * (1 - y) / (1 - t * y))
      (Icc a 1) (Icc (0 : ℝ) 1) := by
    intro t ht
    have ht0 : 0 ≤ t := le_trans ha0.le ht.1
    have hd0 := hden_pos t ht
    constructor
    · exact div_nonneg (mul_nonneg ht0 (sub_nonneg.mpr hy1.le)) hd0.le
    · apply (div_le_one hd0).2
      nlinarith [ht.2]
  have hmap_c : MapsTo
      (fun t : ℝ => y * (1 - t) / (1 - t * y))
      (Icc a 1) (Icc (0 : ℝ) 1) := by
    intro t ht
    have hd0 := hden_pos t ht
    constructor
    · exact div_nonneg
        (mul_nonneg hy0.le (sub_nonneg.mpr ht.2)) hd0.le
    · apply (div_le_one hd0).2
      nlinarith [hy1.le]
  have hrx : ContinuousOn (fun t : ℝ => rogersDialogarithm t)
      (Icc a 1) :=
    rogersDialogarithm_continuousOn_unit.comp continuousOn_id hmap_x
  have hrxy : ContinuousOn (fun t : ℝ => rogersDialogarithm (t * y))
      (Icc a 1) := by
    apply rogersDialogarithm_continuousOn_unit.comp
    · fun_prop
    · exact hmap_xy
  have hrb : ContinuousOn
      (fun t : ℝ => rogersDialogarithm
        (t * (1 - y) / (1 - t * y))) (Icc a 1) := by
    exact rogersDialogarithm_continuousOn_unit.comp hb_cont hmap_b
  have hrc : ContinuousOn
      (fun t : ℝ => rogersDialogarithm
        (y * (1 - t) / (1 - t * y))) (Icc a 1) := by
    exact rogersDialogarithm_continuousOn_unit.comp hc_cont hmap_c
  exact (((hrx.add continuousOn_const).sub hrxy).sub hrb).sub hrc

theorem rogers_five_term {x y : ℝ}
    (hx0 : 0 < x) (hx1 : x < 1) (hy0 : 0 < y) (hy1 : y < 1) :
    rogersDialogarithm x + rogersDialogarithm y =
      rogersDialogarithm (x * y) +
      rogersDialogarithm (x * (1 - y) / (1 - x * y)) +
      rogersDialogarithm (y * (1 - x) / (1 - x * y)) := by
  let F := fun t : ℝ =>
    rogersDialogarithm t + rogersDialogarithm y -
      rogersDialogarithm (t * y) -
      rogersDialogarithm (t * (1 - y) / (1 - t * y)) -
      rogersDialogarithm (y * (1 - t) / (1 - t * y))
  have hcont : ContinuousOn F (Icc x 1) := by
    simpa only [F] using rogersFiveAux_continuousOn hx0 hy0 hy1
  have hconst := constant_of_has_deriv_right_zero
    hcont
    (fun t ht => by
      have hderiv := rogersFiveAux_hasDerivAt_zero
        (lt_of_lt_of_le hx0 ht.1) ht.2 hy0 hy1
      simpa only [F] using hderiv.hasDerivWithinAt)
  have h := hconst 1 (right_mem_Icc.mpr hx1.le)
  have h1y_ne : 1 - y ≠ 0 := ne_of_gt (sub_pos.mpr hy1)
  have hone : F 1 = 0 := by
    simp [F, h1y_ne, rogers_zero]
  rw [hone] at h
  dsimp only [F] at h
  linarith

end
