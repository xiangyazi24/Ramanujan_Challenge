import RamanujanChallenge.Q6427Scratch

open Filter Set Topology MeasureTheory
open scoped Interval Real BigOperators

noncomputable section

namespace RamanujanChallenge.P27.Q6427

/-! ## Explicit sine-kernel decay -/

theorem norm_sinePi_top_ge (x T : ℝ) (hT : 1 ≤ T) :
    Real.exp (Real.pi * T) / 4 ≤
      ‖sinePi ((x : ℂ) + (T : ℂ) * Complex.I)‖ := by
  let w : ℂ := (Real.pi : ℂ) *
    ((x : ℂ) + (T : ℂ) * Complex.I)
  have hlarge : ‖Complex.exp (-w * Complex.I)‖ =
      Real.exp (Real.pi * T) := by
    rw [Complex.norm_exp]
    congr 1
    simp [w]
    ring
  have hsmallNorm : ‖Complex.exp (w * Complex.I)‖ =
      Real.exp (-(Real.pi * T)) := by
    rw [Complex.norm_exp]
    congr 1
    simp [w]
    ring
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
  linarith

theorem norm_sinePi_bottom_ge (x T : ℝ) (hT : 1 ≤ T) :
    Real.exp (Real.pi * T) / 4 ≤
      ‖sinePi ((x : ℂ) - (T : ℂ) * Complex.I)‖ := by
  let z : ℂ := (x : ℂ) + (T : ℂ) * Complex.I
  have hconj :
      (x : ℂ) - (T : ℂ) * Complex.I = Complex.conj z := by
    apply Complex.ext <;> simp [z]
  have hsineConj : sinePi (Complex.conj z) = Complex.conj (sinePi z) := by
    simp [sinePi, Complex.sin_conj]
  rw [hconj, hsineConj, map_norm]
  exact norm_sinePi_top_ge x T hT

private theorem norm_ctKernel_le_of_sine_ge
    {t : ℂ} {E : ℝ} (hE : 0 < E)
    (hs : E / 4 ≤ ‖sinePi t‖) :
    ‖ctKernel27 t‖ ≤ 16 * Real.pi ^ 2 / E ^ 2 := by
  have hspos : 0 < ‖sinePi t‖ :=
    lt_of_lt_of_le (div_pos hE (by norm_num)) hs
  have hsq : (E / 4) ^ 2 ≤ ‖sinePi t‖ ^ 2 := by
    nlinarith [sq_nonneg (‖sinePi t‖ - E / 4)]
  rw [ctKernel27, norm_pow, norm_div, Complex.norm_real,
    Real.norm_eq_abs, abs_of_pos Real.pi_pos]
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

/-! ## Coarse rational-function growth, uniform in one strip -/

theorem abs_strip_sub_nat_le {n m r : ℕ}
    (hm1 : 1 ≤ m) (hmn : m ≤ n) (hr : r < n)
    {x : ℝ}
    (hx : x ∈ Set.Icc ((m : ℝ) - 1 / 2) ((m : ℝ) + 1 / 2)) :
    |x - ((r + 1 : ℕ) : ℝ)| ≤ (n : ℝ) + 1 := by
  have hrn : r + 1 ≤ n := by omega
  have hmR : (m : ℝ) ≤ n := by exact_mod_cast hmn
  have hrR : ((r + 1 : ℕ) : ℝ) ≤ n := by exact_mod_cast hrn
  have hr1 : (1 : ℝ) ≤ (r + 1 : ℕ) := by positivity
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
    _ = |x - ((r + 1 : ℕ) : ℝ)| + |y| := by simp
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
        exact pow_le_pow_left₀
          (norm_horizontal_factor_le hm1 hmn (Finset.mem_range.mp hr) hx) 3
    _ = ((((n : ℝ) + 2) * (1 + |y|)) ^ 3) ^ n := by simp

theorem half_le_norm_pole_factor {m j : ℕ}
    (hm1 : 1 ≤ m) {t : ℂ}
    (ht : t ∈ halfIntegerStrip (m : ℤ)) :
    (1 / 2 : ℝ) ≤ ‖t + (j : ℂ)‖ := by
  have hmR : (1 : ℝ) ≤ m := by exact_mod_cast hm1
  have hj : (0 : ℝ) ≤ j := Nat.cast_nonneg j
  have hre : (1 / 2 : ℝ) ≤ (t + (j : ℂ)).re := by
    simp
    have := ht.1
    norm_num [halfIntegerStrip] at this ⊢
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
  have hfac : (1 : ℝ) ≤ (n.factorial : ℝ) ^ 2 := by
    have : (1 : ℝ) ≤ n.factorial := by exact_mod_cast Nat.one_le_factorial n
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
  have hnum := norm_ctNumerator_on_strip_le hm1 hmn hx
  have hden := norm_ctDen_on_strip_ge (n := n) hm1 ht
  have hhalf : 0 < (1 / 2 : ℝ) ^ (n + 1) := by positivity
  rw [ctR27, norm_div]
  calc
    ‖ctNumerator27 n t‖ /
        ‖((n.factorial : ℂ) ^ 2) * ctPoleProduct27 (n + 1) t‖ ≤
      ‖ctNumerator27 n t‖ / (1 / 2 : ℝ) ^ (n + 1) := by
        exact div_le_div_of_nonneg_left (norm_nonneg _) hhalf hden
    _ ≤ ((((n : ℝ) + 2) * (1 + |y|)) ^ 3) ^ n /
        (1 / 2 : ℝ) ^ (n + 1) :=
      div_le_div_of_nonneg_right hnum hhalf.le

/-! ## The complete horizontal-edge estimate -/

def ctHorizontalMajorant27 (n : ℕ) (T : ℝ) : ℝ :=
  (16 * Real.pi ^ 2 /
      (1 / 2 : ℝ) ^ (n + 1)) *
    ((2 * ((n : ℝ) + 2)) ^ (3 * n)) *
    T ^ (3 * n) / Real.exp (Real.pi * T) ^ 2

theorem norm_ctIntegrand_top_le {n m : ℕ}
    (hm1 : 1 ≤ m) (hmn : m ≤ n)
    {x T : ℝ}
    (hx : x ∈ Set.Icc ((m : ℝ) - 1 / 2) ((m : ℝ) + 1 / 2))
    (hT : 1 ≤ T) :
    ‖ctIntegrand27 n ((x : ℂ) + (T : ℂ) * Complex.I)‖ ≤
      ctHorizontalMajorant27 n T := by
  have hR := norm_ctR_on_strip_le hm1 hmn hx (y := T)
  have hK := norm_ctKernel_top_le x T hT
  have hTabs : |T| = T := abs_of_nonneg (zero_le_one.trans hT)
  have hbase : 1 + T ≤ 2 * T := by linarith
  have hpoly :
      ((((n : ℝ) + 2) * (1 + |T|)) ^ 3) ^ n ≤
        (2 * ((n : ℝ) + 2)) ^ (3 * n) * T ^ (3 * n) := by
    rw [hTabs]
    have hn2 : 0 ≤ (n : ℝ) + 2 := by positivity
    have hfactor : ((n : ℝ) + 2) * (1 + T) ≤
        (2 * ((n : ℝ) + 2)) * T := by nlinarith
    calc
      ((((n : ℝ) + 2) * (1 + T)) ^ 3) ^ n =
          (((n : ℝ) + 2) * (1 + T)) ^ (3 * n) := by rw [pow_mul]
      _ ≤ ((2 * ((n : ℝ) + 2)) * T) ^ (3 * n) := by
        exact pow_le_pow_left₀ hfactor (3 * n)
      _ = (2 * ((n : ℝ) + 2)) ^ (3 * n) * T ^ (3 * n) := by
        rw [mul_pow]
  rw [ctIntegrand27, norm_mul]
  calc
    ‖ctR27 n ((x : ℂ) + (T : ℂ) * Complex.I)‖ *
        ‖ctKernel27 ((x : ℂ) + (T : ℂ) * Complex.I)‖ ≤
      (((((n : ℝ) + 2) * (1 + |T|)) ^ 3) ^ n /
          (1 / 2 : ℝ) ^ (n + 1)) *
        (16 * Real.pi ^ 2 / Real.exp (Real.pi * T) ^ 2) := by
      exact mul_le_mul hR hK (norm_nonneg _) (by positivity)
    _ ≤ ctHorizontalMajorant27 n T := by
      unfold ctHorizontalMajorant27
      have hhalf : 0 < (1 / 2 : ℝ) ^ (n + 1) := by positivity
      have hE : 0 < Real.exp (Real.pi * T) ^ 2 := sq_pos_of_pos (Real.exp_pos _)
      apply (div_le_div_iff₀ hE).2
      apply (div_le_div_iff₀ hhalf).2
      nlinarith [hpoly]

theorem norm_ctIntegrand_bottom_le {n m : ℕ}
    (hm1 : 1 ≤ m) (hmn : m ≤ n)
    {x T : ℝ}
    (hx : x ∈ Set.Icc ((m : ℝ) - 1 / 2) ((m : ℝ) + 1 / 2))
    (hT : 1 ≤ T) :
    ‖ctIntegrand27 n ((x : ℂ) - (T : ℂ) * Complex.I)‖ ≤
      ctHorizontalMajorant27 n T := by
  have hR := norm_ctR_on_strip_le hm1 hmn hx (y := -T)
  have hK := norm_ctKernel_bottom_le x T hT
  have hTabs : |-T| = T := abs_neg T ▸ abs_of_nonneg (zero_le_one.trans hT)
  have hbase : 1 + T ≤ 2 * T := by linarith
  have hpoly :
      ((((n : ℝ) + 2) * (1 + |-T|)) ^ 3) ^ n ≤
        (2 * ((n : ℝ) + 2)) ^ (3 * n) * T ^ (3 * n) := by
    rw [hTabs]
    have hfactor : ((n : ℝ) + 2) * (1 + T) ≤
        (2 * ((n : ℝ) + 2)) * T := by
      have hn2 : 0 ≤ (n : ℝ) + 2 := by positivity
      nlinarith
    calc
      ((((n : ℝ) + 2) * (1 + T)) ^ 3) ^ n =
          (((n : ℝ) + 2) * (1 + T)) ^ (3 * n) := by rw [pow_mul]
      _ ≤ ((2 * ((n : ℝ) + 2)) * T) ^ (3 * n) := by
        exact pow_le_pow_left₀ hfactor (3 * n)
      _ = (2 * ((n : ℝ) + 2)) ^ (3 * n) * T ^ (3 * n) := by
        rw [mul_pow]
  have hpoint :
      (x : ℂ) - (T : ℂ) * Complex.I =
        (x : ℂ) + ((-T : ℝ) : ℂ) * Complex.I := by ring
  rw [ctIntegrand27, norm_mul, hpoint]
  calc
    ‖ctR27 n ((x : ℂ) + ((-T : ℝ) : ℂ) * Complex.I)‖ *
        ‖ctKernel27 ((x : ℂ) - (T : ℂ) * Complex.I)‖ ≤
      (((((n : ℝ) + 2) * (1 + |-T|)) ^ 3) ^ n /
          (1 / 2 : ℝ) ^ (n + 1)) *
        (16 * Real.pi ^ 2 / Real.exp (Real.pi * T) ^ 2) := by
      exact mul_le_mul hR hK (norm_nonneg _) (by positivity)
    _ ≤ ctHorizontalMajorant27 n T := by
      unfold ctHorizontalMajorant27
      have hhalf : 0 < (1 / 2 : ℝ) ^ (n + 1) := by positivity
      have hE : 0 < Real.exp (Real.pi * T) ^ 2 := sq_pos_of_pos (Real.exp_pos _)
      apply (div_le_div_iff₀ hE).2
      apply (div_le_div_iff₀ hhalf).2
      nlinarith [hpoly]

theorem ctHorizontalMajorant_tendsto_zero27 (n : ℕ) :
    Tendsto (ctHorizontalMajorant27 n) atTop (𝓝 0) := by
  have hbase := (Real.tendsto_pow_mul_exp_neg_atTop_nhds_zero (3 * n)).comp
    (tendsto_id.const_mul_atTop' Real.pi_pos)
  have hscaled :
      Tendsto
        (fun T : ℝ => T ^ (3 * n) * Real.exp (-(2 * Real.pi * T)))
        atTop (𝓝 0) := by
    have hpow :
        Tendsto
          (fun T : ℝ => (Real.pi * T) ^ (3 * n) *
            Real.exp (-(2 * (Real.pi * T)))) atTop (𝓝 0) := by
      simpa [mul_assoc] using
        (Real.tendsto_pow_mul_exp_neg_atTop_nhds_zero (3 * n)).comp
          (tendsto_id.const_mul_atTop' (by positivity : 0 < (2 * Real.pi)))
    have hpi : Real.pi ^ (3 * n) ≠ 0 :=
      pow_ne_zero _ Real.pi_ne_zero
    convert hpow.const_mul (Real.pi ^ (3 * n))⁻¹ using 1 <;>
      field_simp [hpi] <;> ring
  have hconst := hscaled.const_mul
    ((16 * Real.pi ^ 2 / (1 / 2 : ℝ) ^ (n + 1)) *
      (2 * ((n : ℝ) + 2)) ^ (3 * n))
  apply hconst.congr'
  filter_upwards with T
  unfold ctHorizontalMajorant27
  rw [show Real.exp (Real.pi * T) ^ 2 =
      Real.exp (2 * Real.pi * T) by
        rw [← Real.exp_add]
        congr 1
        ring,
    div_eq_mul_inv, ← Real.exp_neg]
  ring

/-- The top horizontal edge tends to zero, with its norm bounded by the
explicit polynomial-exponential majorant above. -/
theorem ctHorizontal_top_tendsto_zero27 {n m : ℕ}
    (hm1 : 1 ≤ m) (hmn : m ≤ n) :
    Tendsto
      (fun T : ℝ => ∫ x in ((m : ℝ) - 1 / 2)..((m : ℝ) + 1 / 2),
        ctExtension27 n m ((x : ℂ) + (T : ℂ) * Complex.I))
      atTop (𝓝 0) := by
  apply squeeze_zero'
  · exact Eventually.of_forall fun T => norm_nonneg _
  · filter_upwards [eventually_ge_atTop (1 : ℝ)] with T hT
    have hraw : ∀ x ∈ [[(m : ℝ) - 1 / 2, (m : ℝ) + 1 / 2]],
        ‖ctExtension27 n m ((x : ℂ) + (T : ℂ) * Complex.I)‖ ≤
          ctHorizontalMajorant27 n T := by
      intro x hx
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
    exact (intervalIntegral.norm_integral_le_of_norm_le_const hraw).trans_eq (by
      unfold ctHorizontalMajorant27
      rw [abs_of_nonneg (by linarith : 0 ≤
        ((m : ℝ) + 1 / 2) - ((m : ℝ) - 1 / 2))]
      ring)
  · exact ctHorizontalMajorant_tendsto_zero27 n

end RamanujanChallenge.P27.Q6427
