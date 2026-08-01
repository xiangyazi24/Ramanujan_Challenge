import RamanujanChallenge.Q6427Scratch

open Filter Set Topology MeasureTheory
open scoped Interval Real BigOperators

noncomputable section

namespace RamanujanChallenge.P27.Q6427

/-! ## Explicit sine-kernel decay on horizontal edges -/

theorem norm_sinePi_top_ge (x T : ℝ) (hT : 1 ≤ T) :
    Real.exp (Real.pi * T) / 4 ≤
      ‖sinePi ((x : ℂ) + (T : ℂ) * Complex.I)‖ := by
  let w : ℂ := (Real.pi : ℂ) *
    ((x : ℂ) + (T : ℂ) * Complex.I)
  have hlarge : ‖Complex.exp (-w * Complex.I)‖ =
      Real.exp (Real.pi * T) := by
    rw [Complex.norm_exp]
    simp [w]
  have hsmallNorm : ‖Complex.exp (w * Complex.I)‖ =
      Real.exp (-(Real.pi * T)) := by
    rw [Complex.norm_exp]
    simp [w]
  have harg : 1 ≤ Real.pi * T := by
    nlinarith [Real.two_le_pi]
  have hexpTwo : 2 ≤ Real.exp (Real.pi * T) := by
    calc
      (2 : ℝ) = 1 + 1 := by norm_num
      _ ≤ Real.exp 1 := Real.add_one_le_exp 1
      _ ≤ Real.exp (Real.pi * T) := Real.exp_le_exp.mpr harg
  have hnegOne : Real.exp (-(Real.pi * T)) ≤ 1 := by
    rw [← Real.exp_zero]
    exact Real.exp_le_exp.mpr (by nlinarith [Real.pi_pos])
  have hsmall :
      Real.exp (-(Real.pi * T)) ≤ Real.exp (Real.pi * T) / 2 := by
    linarith
  have hreverse :=
    norm_sub_norm_le (Complex.exp (-w * Complex.I))
      (Complex.exp (w * Complex.I))
  rw [hlarge, hsmallNorm] at hreverse
  have hdiff : Real.exp (Real.pi * T) / 2 ≤
      ‖Complex.exp (-w * Complex.I) -
        Complex.exp (w * Complex.I)‖ := by
    linarith
  have hsin :
      sinePi ((x : ℂ) + (T : ℂ) * Complex.I) =
        (Complex.exp (-w * Complex.I) -
          Complex.exp (w * Complex.I)) * Complex.I / 2 := by
    rfl
  rw [hsin, norm_div, norm_mul]
  norm_num
  nlinarith

theorem norm_sinePi_bottom_ge (x T : ℝ) (hT : 1 ≤ T) :
    Real.exp (Real.pi * T) / 4 ≤
      ‖sinePi ((x : ℂ) - (T : ℂ) * Complex.I)‖ := by
  let w : ℂ := (Real.pi : ℂ) *
    ((x : ℂ) - (T : ℂ) * Complex.I)
  have hlarge : ‖Complex.exp (w * Complex.I)‖ =
      Real.exp (Real.pi * T) := by
    rw [Complex.norm_exp]
    simp [w]
  have hsmallNorm : ‖Complex.exp (-w * Complex.I)‖ =
      Real.exp (-(Real.pi * T)) := by
    rw [Complex.norm_exp]
    simp [w]
  have harg : 1 ≤ Real.pi * T := by
    nlinarith [Real.two_le_pi]
  have hexpTwo : 2 ≤ Real.exp (Real.pi * T) := by
    calc
      (2 : ℝ) = 1 + 1 := by norm_num
      _ ≤ Real.exp 1 := Real.add_one_le_exp 1
      _ ≤ Real.exp (Real.pi * T) := Real.exp_le_exp.mpr harg
  have hnegOne : Real.exp (-(Real.pi * T)) ≤ 1 := by
    rw [← Real.exp_zero]
    exact Real.exp_le_exp.mpr (by nlinarith [Real.pi_pos])
  have hsmall :
      Real.exp (-(Real.pi * T)) ≤ Real.exp (Real.pi * T) / 2 := by
    linarith
  have hreverse :=
    norm_sub_norm_le (Complex.exp (w * Complex.I))
      (Complex.exp (-w * Complex.I))
  rw [hlarge, hsmallNorm] at hreverse
  have hdiff : Real.exp (Real.pi * T) / 2 ≤
      ‖Complex.exp (w * Complex.I) -
        Complex.exp (-w * Complex.I)‖ := by
    linarith
  have hsin :
      sinePi ((x : ℂ) - (T : ℂ) * Complex.I) =
        -((Complex.exp (w * Complex.I) -
          Complex.exp (-w * Complex.I)) * Complex.I / 2) := by
    unfold sinePi
    rw [Complex.sin]
    ring
  rw [hsin, norm_neg, norm_div, norm_mul]
  norm_num
  nlinarith

private theorem norm_ctKernel_le_of_sine_ge
    {t : ℂ} {E : ℝ} (hE : 0 < E)
    (hs : E / 4 ≤ ‖sinePi t‖) :
    ‖ctKernel27 t‖ ≤ 16 * Real.pi ^ 2 / E ^ 2 := by
  have hspos : 0 < ‖sinePi t‖ :=
    lt_of_lt_of_le (div_pos hE (by norm_num)) hs
  have hsq : (E / 4) ^ 2 ≤ ‖sinePi t‖ ^ 2 := by
    gcongr
  rw [ctKernel27, norm_pow, norm_div, Complex.norm_real,
    Real.norm_eq_abs, abs_of_pos Real.pi_pos, div_pow]
  calc
    Real.pi ^ 2 / ‖sinePi t‖ ^ 2 ≤
        Real.pi ^ 2 / (E / 4) ^ 2 := by
      exact div_le_div_of_nonneg_left (sq_nonneg Real.pi)
        (sq_pos_of_pos (div_pos hE (by norm_num))) hsq
    _ = 16 * Real.pi ^ 2 / E ^ 2 := by
      field_simp
      ring

theorem norm_ctKernel_top_le (x T : ℝ) (hT : 1 ≤ T) :
    ‖ctKernel27 ((x : ℂ) + (T : ℂ) * Complex.I)‖ ≤
      16 * Real.pi ^ 2 / Real.exp (Real.pi * T) ^ 2 :=
  norm_ctKernel_le_of_sine_ge (Real.exp_pos _)
    (norm_sinePi_top_ge x T hT)

theorem norm_ctKernel_bottom_le (x T : ℝ) (hT : 1 ≤ T) :
    ‖ctKernel27 ((x : ℂ) - (T : ℂ) * Complex.I)‖ ≤
      16 * Real.pi ^ 2 / Real.exp (Real.pi * T) ^ 2 :=
  norm_ctKernel_le_of_sine_ge (Real.exp_pos _)
    (norm_sinePi_bottom_ge x T hT)

/-! ## Polynomial growth of the rational factor -/

theorem abs_strip_sub_nat_le {n m r : ℕ}
    (hm1 : 1 ≤ m) (hmn : m ≤ n) (hr : r < n)
    {x : ℝ}
    (hx : x ∈ Set.Icc ((m : ℝ) - 1 / 2) ((m : ℝ) + 1 / 2)) :
    |x - ((r + 1 : ℕ) : ℝ)| ≤ (n : ℝ) + 1 := by
  have hrn : r + 1 ≤ n := by omega
  have hmR : (m : ℝ) ≤ n := by exact_mod_cast hmn
  have hrR : ((r + 1 : ℕ) : ℝ) ≤ n := by exact_mod_cast hrn
  have hr1 : (1 : ℝ) ≤ (r + 1 : ℕ) := by exact_mod_cast Nat.succ_le_succ (Nat.zero_le r)
  have hm1R : (1 : ℝ) ≤ m := by exact_mod_cast hm1
  rw [abs_le]
  constructor <;> linarith [hx.1, hx.2]

theorem norm_horizontal_factor_le {n m r : ℕ}
    (hm1 : 1 ≤ m) (hmn : m ≤ n) (hr : r < n)
    {x y : ℝ}
    (hx : x ∈ Set.Icc ((m : ℝ) - 1 / 2) ((m : ℝ) + 1 / 2)) :
    ‖(x : ℂ) + (y : ℂ) * Complex.I - (((r + 1 : ℕ) : ℂ))‖ ≤
      ((n : ℝ) + 2) * (1 + |y|) := by
  have habs := abs_strip_sub_nat_le hm1 hmn hr hx
  have heq :
      (x : ℂ) + (y : ℂ) * Complex.I - (((r + 1 : ℕ) : ℂ)) =
        ((x - ((r + 1 : ℕ) : ℝ) : ℝ) : ℂ) +
          (y : ℂ) * Complex.I := by
    push_cast
    ring
  rw [heq]
  calc
    ‖((x - ((r + 1 : ℕ) : ℝ) : ℝ) : ℂ) +
        (y : ℂ) * Complex.I‖ ≤
      ‖((x - ((r + 1 : ℕ) : ℝ) : ℝ) : ℂ)‖ +
        ‖(y : ℂ) * Complex.I‖ := norm_add_le _ _
    _ = |x - ((r + 1 : ℕ) : ℝ)| + |y| := by
      rw [Complex.norm_real, Real.norm_eq_abs, norm_mul,
        Complex.norm_real, Real.norm_eq_abs, norm_I, mul_one]
    _ ≤ ((n : ℝ) + 1) + |y| := by linarith
    _ ≤ ((n : ℝ) + 2) * (1 + |y|) := by
      have hn : 0 ≤ (n : ℝ) := Nat.cast_nonneg n
      have hy : 0 ≤ |y| := abs_nonneg y
      nlinarith

theorem norm_ctNumerator_on_strip_le {n m : ℕ}
    (hm1 : 1 ≤ m) (hmn : m ≤ n)
    {x y : ℝ}
    (hx : x ∈ Set.Icc ((m : ℝ) - 1 / 2) ((m : ℝ) + 1 / 2)) :
    ‖ctNumerator27 n ((x : ℂ) + (y : ℂ) * Complex.I)‖ ≤
      ((((n : ℝ) + 2) * (1 + |y|)) ^ 3) ^ n := by
  rw [ctNumerator27, norm_prod]
  calc
    ∏ r ∈ Finset.range n,
        ‖((x : ℂ) + (y : ℂ) * Complex.I -
          (((r + 1 : ℕ) : ℂ))) ^ 3‖ ≤
      ∏ _r ∈ Finset.range n,
        (((n : ℝ) + 2) * (1 + |y|)) ^ 3 := by
      exact Finset.prod_le_prod (fun _ _ => norm_nonneg _) fun r hr => by
        rw [norm_pow]
        gcongr
        exact norm_horizontal_factor_le hm1 hmn (Finset.mem_range.mp hr) hx
    _ = ((((n : ℝ) + 2) * (1 + |y|)) ^ 3) ^ n := by simp

theorem half_le_norm_pole_factor {m j : ℕ}
    (hm1 : 1 ≤ m) {t : ℂ}
    (ht : t ∈ halfIntegerStrip (m : ℤ)) :
    (1 / 2 : ℝ) ≤ ‖t + (j : ℂ)‖ := by
  have hmR : (1 : ℝ) ≤ m := by exact_mod_cast hm1
  have hj : (0 : ℝ) ≤ j := Nat.cast_nonneg j
  have hre : (1 / 2 : ℝ) ≤ (t + (j : ℂ)).re := by
    simp
    have hleft := ht.1
    norm_num [halfIntegerStrip] at hleft ⊢
    linarith
  exact hre.trans ((le_abs_self _).trans (Complex.abs_re_le_norm _))

theorem norm_ctDen_on_strip_ge {n m : ℕ}
    (hm1 : 1 ≤ m) {t : ℂ}
    (ht : t ∈ halfIntegerStrip (m : ℤ)) :
    (1 / 2 : ℝ) ^ (n + 1) ≤
      ‖((n.factorial : ℂ) ^ 2) * ctPoleProduct27 (n + 1) t‖ := by
  have hpole :
      (1 / 2 : ℝ) ^ (n + 1) ≤ ‖ctPoleProduct27 (n + 1) t‖ := by
    rw [ctPoleProduct27, norm_prod]
    calc
      (1 / 2 : ℝ) ^ (n + 1) =
          ∏ _j ∈ Finset.range (n + 1), (1 / 2 : ℝ) := by simp
      _ ≤ ∏ j ∈ Finset.range (n + 1), ‖t + (j : ℂ)‖ := by
        exact Finset.prod_le_prod (fun _ _ => by positivity) fun j _ =>
          half_le_norm_pole_factor hm1 ht
  rw [norm_mul, norm_pow, Complex.norm_natCast]
  have hfacNat : 1 ≤ n.factorial :=
    Nat.one_le_iff_ne_zero.mpr (Nat.factorial_ne_zero n)
  have hfac : (1 : ℝ) ≤ (n.factorial : ℝ) ^ 2 := by
    have hcast : (1 : ℝ) ≤ n.factorial := by exact_mod_cast hfacNat
    nlinarith
  have hpnonneg : 0 ≤ ‖ctPoleProduct27 (n + 1) t‖ := norm_nonneg _
  calc
    (1 / 2 : ℝ) ^ (n + 1) ≤ ‖ctPoleProduct27 (n + 1) t‖ := hpole
    _ ≤ (n.factorial : ℝ) ^ 2 * ‖ctPoleProduct27 (n + 1) t‖ := by
      nlinarith

theorem norm_ctR_on_strip_le {n m : ℕ}
    (hm1 : 1 ≤ m) (hmn : m ≤ n)
    {x y : ℝ}
    (hx : x ∈ Set.Icc ((m : ℝ) - 1 / 2) ((m : ℝ) + 1 / 2)) :
    ‖ctR27 n ((x : ℂ) + (y : ℂ) * Complex.I)‖ ≤
      ((((n : ℝ) + 2) * (1 + |y|)) ^ 3) ^ n /
        (1 / 2 : ℝ) ^ (n + 1) := by
  let t : ℂ := (x : ℂ) + (y : ℂ) * Complex.I
  have ht : t ∈ halfIntegerStrip (m : ℤ) := by
    simpa [t, halfIntegerStrip] using hx
  have hnum := norm_ctNumerator_on_strip_le hm1 hmn (y := y) hx
  have hden := norm_ctDen_on_strip_ge (n := n) hm1 ht
  have hhalf : 0 < (1 / 2 : ℝ) ^ (n + 1) := by positivity
  rw [ctR27, norm_div]
  change ‖ctNumerator27 n t‖ /
      ‖((n.factorial : ℂ) ^ 2) * ctPoleProduct27 (n + 1) t‖ ≤ _
  calc
    ‖ctNumerator27 n t‖ /
        ‖((n.factorial : ℂ) ^ 2) * ctPoleProduct27 (n + 1) t‖ ≤
      ‖ctNumerator27 n t‖ / (1 / 2 : ℝ) ^ (n + 1) := by
        exact div_le_div_of_nonneg_left (norm_nonneg _) hhalf hden
    _ ≤ ((((n : ℝ) + 2) * (1 + |y|)) ^ 3) ^ n /
        (1 / 2 : ℝ) ^ (n + 1) :=
      div_le_div_of_nonneg_right hnum hhalf.le

/-! ## Complete horizontal-edge majorant -/

def ctHorizontalMajorant27 (n : ℕ) (T : ℝ) : ℝ :=
  (16 * Real.pi ^ 2 /
      (1 / 2 : ℝ) ^ (n + 1)) *
    ((2 * ((n : ℝ) + 2)) ^ (3 * n)) *
    (T ^ (3 * n) * Real.exp (-(2 * Real.pi * T)))

theorem norm_ctIntegrand_top_le {n m : ℕ}
    (hm1 : 1 ≤ m) (hmn : m ≤ n)
    {x T : ℝ}
    (hx : x ∈ Set.Icc ((m : ℝ) - 1 / 2) ((m : ℝ) + 1 / 2))
    (hT : 1 ≤ T) :
    ‖ctIntegrand27 n ((x : ℂ) + (T : ℂ) * Complex.I)‖ ≤
      ctHorizontalMajorant27 n T := by
  have hR := norm_ctR_on_strip_le hm1 hmn (y := T) hx
  have hK := norm_ctKernel_top_le x T hT
  have hTabs : |T| = T := abs_of_nonneg (zero_le_one.trans hT)
  have hfactor : ((n : ℝ) + 2) * (1 + T) ≤
      (2 * ((n : ℝ) + 2)) * T := by
    have hn2 : 0 ≤ (n : ℝ) + 2 := by positivity
    nlinarith
  have hpoly :
      ((((n : ℝ) + 2) * (1 + |T|)) ^ 3) ^ n ≤
        (2 * ((n : ℝ) + 2)) ^ (3 * n) * T ^ (3 * n) := by
    rw [hTabs, pow_mul, mul_pow]
    gcongr
  have hE : Real.exp (Real.pi * T) ^ 2 =
      Real.exp (2 * Real.pi * T) := by
    rw [← Real.exp_add]
    congr 1
    ring
  rw [ctIntegrand27, norm_mul]
  calc
    ‖ctR27 n ((x : ℂ) + (T : ℂ) * Complex.I)‖ *
        ‖ctKernel27 ((x : ℂ) + (T : ℂ) * Complex.I)‖ ≤
      (((((n : ℝ) + 2) * (1 + |T|)) ^ 3) ^ n /
          (1 / 2 : ℝ) ^ (n + 1)) *
        (16 * Real.pi ^ 2 / Real.exp (Real.pi * T) ^ 2) := by
      exact mul_le_mul hR hK (norm_nonneg _) (by positivity)
    _ ≤ ctHorizontalMajorant27 n T := by
      rw [hE]
      unfold ctHorizontalMajorant27
      have hh : 0 < (1 / 2 : ℝ) ^ (n + 1) := by positivity
      have he : 0 < Real.exp (2 * Real.pi * T) := Real.exp_pos _
      field_simp [ne_of_gt hh, ne_of_gt he]
      nlinarith [hpoly]

theorem norm_ctIntegrand_bottom_le {n m : ℕ}
    (hm1 : 1 ≤ m) (hmn : m ≤ n)
    {x T : ℝ}
    (hx : x ∈ Set.Icc ((m : ℝ) - 1 / 2) ((m : ℝ) + 1 / 2))
    (hT : 1 ≤ T) :
    ‖ctIntegrand27 n ((x : ℂ) - (T : ℂ) * Complex.I)‖ ≤
      ctHorizontalMajorant27 n T := by
  have hpoint :
      (x : ℂ) - (T : ℂ) * Complex.I =
        (x : ℂ) + ((-T : ℝ) : ℂ) * Complex.I := by ring
  rw [hpoint]
  have hR := norm_ctR_on_strip_le hm1 hmn (y := -T) hx
  have hK :
      ‖ctKernel27 ((x : ℂ) + ((-T : ℝ) : ℂ) * Complex.I)‖ ≤
        16 * Real.pi ^ 2 / Real.exp (Real.pi * T) ^ 2 := by
    simpa [hpoint] using norm_ctKernel_bottom_le x T hT
  have hTabs : |-T| = T := by simpa using abs_of_nonneg (zero_le_one.trans hT)
  have hfactor : ((n : ℝ) + 2) * (1 + T) ≤
      (2 * ((n : ℝ) + 2)) * T := by
    have hn2 : 0 ≤ (n : ℝ) + 2 := by positivity
    nlinarith
  have hpoly :
      ((((n : ℝ) + 2) * (1 + |-T|)) ^ 3) ^ n ≤
        (2 * ((n : ℝ) + 2)) ^ (3 * n) * T ^ (3 * n) := by
    rw [hTabs, pow_mul, mul_pow]
    gcongr
  have hE : Real.exp (Real.pi * T) ^ 2 =
      Real.exp (2 * Real.pi * T) := by
    rw [← Real.exp_add]
    congr 1
    ring
  rw [ctIntegrand27, norm_mul]
  calc
    ‖ctR27 n ((x : ℂ) + ((-T : ℝ) : ℂ) * Complex.I)‖ *
        ‖ctKernel27 ((x : ℂ) + ((-T : ℝ) : ℂ) * Complex.I)‖ ≤
      (((((n : ℝ) + 2) * (1 + |-T|)) ^ 3) ^ n /
          (1 / 2 : ℝ) ^ (n + 1)) *
        (16 * Real.pi ^ 2 / Real.exp (Real.pi * T) ^ 2) := by
      exact mul_le_mul hR hK (norm_nonneg _) (by positivity)
    _ ≤ ctHorizontalMajorant27 n T := by
      rw [hE]
      unfold ctHorizontalMajorant27
      have hh : 0 < (1 / 2 : ℝ) ^ (n + 1) := by positivity
      have he : 0 < Real.exp (2 * Real.pi * T) := Real.exp_pos _
      field_simp [ne_of_gt hh, ne_of_gt he]
      nlinarith [hpoly]

theorem tendsto_pow_mul_exp_neg_two_pi27 (p : ℕ) :
    Tendsto (fun T : ℝ => T ^ p * Real.exp (-(2 * Real.pi * T)))
      atTop (𝓝 0) := by
  have hraw := (Real.tendsto_pow_mul_exp_neg_atTop_nhds_zero p).comp
    (tendsto_id.const_mul_atTop' Real.two_pi_pos)
  have hconst : (2 * Real.pi : ℝ) ^ p ≠ 0 :=
    pow_ne_zero _ (mul_ne_zero (by norm_num) Real.pi_ne_zero)
  have hscaled := tendsto_const_nhds.mul hraw
    (c := ((2 * Real.pi : ℝ) ^ p)⁻¹)
  convert hscaled using 1
  · funext T
    simp only [Function.comp_apply, id_eq]
    rw [mul_pow]
    field_simp [hconst]
    ring
  · simp

theorem ctHorizontalMajorant_tendsto_zero27 (n : ℕ) :
    Tendsto (ctHorizontalMajorant27 n) atTop (𝓝 0) := by
  unfold ctHorizontalMajorant27
  simpa using
    (tendsto_pow_mul_exp_neg_two_pi27 (3 * n)).const_mul
      ((16 * Real.pi ^ 2 / (1 / 2 : ℝ) ^ (n + 1)) *
        (2 * ((n : ℝ) + 2)) ^ (3 * n))

theorem ctHorizontal_top_tendsto_zero27 {n m : ℕ}
    (hm1 : 1 ≤ m) (hmn : m ≤ n) :
    Tendsto
      (fun T : ℝ => ∫ x in ((m : ℝ) - 1 / 2)..((m : ℝ) + 1 / 2),
        ctExtension27 n m ((x : ℂ) + (T : ℂ) * Complex.I))
      atTop (𝓝 0) := by
  rw [tendsto_zero_iff_norm_tendsto_zero]
  apply squeeze_zero'
  · exact Eventually.of_forall fun T => norm_nonneg _
  · filter_upwards [eventually_ge_atTop (1 : ℝ)] with T hT
    apply (intervalIntegral.norm_integral_le_of_norm_le_const _).trans
    · intro x hx
      have hx' : x ∈ Set.Icc ((m : ℝ) - 1 / 2) ((m : ℝ) + 1 / 2) := by
        simpa [uIcc_of_le (by linarith : (m : ℝ) - 1 / 2 ≤ (m : ℝ) + 1 / 2)] using hx
      have hmem :
          (x : ℂ) + (T : ℂ) * Complex.I ∈ halfIntegerStrip (m : ℤ) := by
        simpa [halfIntegerStrip] using hx'
      have hne : (x : ℂ) + (T : ℂ) * Complex.I ≠ (m : ℂ) := by
        intro h
        have him := congrArg Complex.im h
        simp at him
        linarith
      rw [← ctIntegrand_eq_extension27 hm1 hmn hmem hne]
      exact norm_ctIntegrand_top_le hm1 hmn hx' hT
    · simp
  · exact ctHorizontalMajorant_tendsto_zero27 n

theorem ctHorizontal_bottom_tendsto_zero27 {n m : ℕ}
    (hm1 : 1 ≤ m) (hmn : m ≤ n) :
    Tendsto
      (fun T : ℝ => ∫ x in ((m : ℝ) - 1 / 2)..((m : ℝ) + 1 / 2),
        ctExtension27 n m ((x : ℂ) - (T : ℂ) * Complex.I))
      atTop (𝓝 0) := by
  rw [tendsto_zero_iff_norm_tendsto_zero]
  apply squeeze_zero'
  · exact Eventually.of_forall fun T => norm_nonneg _
  · filter_upwards [eventually_ge_atTop (1 : ℝ)] with T hT
    apply (intervalIntegral.norm_integral_le_of_norm_le_const _).trans
    · intro x hx
      have hx' : x ∈ Set.Icc ((m : ℝ) - 1 / 2) ((m : ℝ) + 1 / 2) := by
        simpa [uIcc_of_le (by linarith : (m : ℝ) - 1 / 2 ≤ (m : ℝ) + 1 / 2)] using hx
      have hmem :
          (x : ℂ) - (T : ℂ) * Complex.I ∈ halfIntegerStrip (m : ℤ) := by
        simpa [halfIntegerStrip] using hx'
      have hne : (x : ℂ) - (T : ℂ) * Complex.I ≠ (m : ℂ) := by
        intro h
        have him := congrArg Complex.im h
        simp at him
        linarith
      rw [← ctIntegrand_eq_extension27 hm1 hmn hmem hne]
      exact norm_ctIntegrand_bottom_le hm1 hmn hx' hT
    · simp
  · exact ctHorizontalMajorant_tendsto_zero27 n

#print axioms ctHorizontal_top_tendsto_zero27
#print axioms ctHorizontal_bottom_tendsto_zero27

end RamanujanChallenge.P27.Q6427
