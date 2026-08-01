import RamanujanChallenge.Q6427Scratch

open Filter Set Topology MeasureTheory
open scoped Interval Real BigOperators

noncomputable section

namespace RamanujanChallenge.P27.Q6427

/-! ## Explicit sine-kernel decay on horizontal edges -/

theorem norm_sinePi_top_ge_fixed (x T : ℝ) (hT : 1 ≤ T) :
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
      ‖Complex.exp (-(w * Complex.I)) -
        Complex.exp (w * Complex.I)‖ := by
    simpa only [neg_mul] using (show Real.exp (Real.pi * T) / 2 ≤
      ‖Complex.exp (-w * Complex.I) - Complex.exp (w * Complex.I)‖ by linarith)
  have hsin :
      sinePi ((x : ℂ) + (T : ℂ) * Complex.I) =
        (Complex.exp (-(w * Complex.I)) -
          Complex.exp (w * Complex.I)) * Complex.I / 2 := by
    unfold sinePi
    rw [Complex.sin]
    change (Complex.exp (-(w * Complex.I)) -
      Complex.exp (w * Complex.I)) * Complex.I / 2 = _
    rfl
  calc
    Real.exp (Real.pi * T) / 4 =
        (Real.exp (Real.pi * T) / 2) / 2 := by ring
    _ ≤ ‖Complex.exp (-(w * Complex.I)) -
          Complex.exp (w * Complex.I)‖ / 2 := by
      exact div_le_div_of_nonneg_right hdiff (by norm_num)
    _ = ‖sinePi ((x : ℂ) + (T : ℂ) * Complex.I)‖ := by
      rw [hsin, norm_div, norm_mul]
      norm_num

theorem norm_sinePi_bottom_ge_fixed (x T : ℝ) (hT : 1 ≤ T) :
    Real.exp (Real.pi * T) / 4 ≤
      ‖sinePi ((x : ℂ) - (T : ℂ) * Complex.I)‖ := by
  let w : ℂ := (Real.pi : ℂ) *
    ((x : ℂ) - (T : ℂ) * Complex.I)
  have hlarge : ‖Complex.exp (w * Complex.I)‖ =
      Real.exp (Real.pi * T) := by
    rw [Complex.norm_exp]
    simp [w]
  have hsmallNorm : ‖Complex.exp (-(w * Complex.I))‖ =
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
      (Complex.exp (-(w * Complex.I)))
  rw [hlarge, hsmallNorm] at hreverse
  have hdiff : Real.exp (Real.pi * T) / 2 ≤
      ‖Complex.exp (w * Complex.I) -
        Complex.exp (-(w * Complex.I))‖ := by
    linarith
  have hsin :
      sinePi ((x : ℂ) - (T : ℂ) * Complex.I) =
        -((Complex.exp (w * Complex.I) -
          Complex.exp (-(w * Complex.I))) * Complex.I / 2) := by
    unfold sinePi
    rw [Complex.sin]
    change (Complex.exp (-(w * Complex.I)) -
      Complex.exp (w * Complex.I)) * Complex.I / 2 =
      -((Complex.exp (w * Complex.I) -
        Complex.exp (-(w * Complex.I))) * Complex.I / 2)
    ring
  calc
    Real.exp (Real.pi * T) / 4 =
        (Real.exp (Real.pi * T) / 2) / 2 := by ring
    _ ≤ ‖Complex.exp (w * Complex.I) -
          Complex.exp (-(w * Complex.I))‖ / 2 := by
      exact div_le_div_of_nonneg_right hdiff (by norm_num)
    _ = ‖sinePi ((x : ℂ) - (T : ℂ) * Complex.I)‖ := by
      rw [hsin, norm_neg, norm_div, norm_mul]
      norm_num

private theorem norm_ctKernel_le_of_sine_ge_fixed
    {t : ℂ} {E : ℝ} (hE : 0 < E)
    (hs : E / 4 ≤ ‖sinePi t‖) :
    ‖ctKernel27 t‖ ≤ 16 * Real.pi ^ 2 / E ^ 2 := by
  have hspos : 0 < ‖sinePi t‖ :=
    lt_of_lt_of_le (div_pos hE (by norm_num)) hs
  have hsq : (E / 4) ^ 2 ≤ ‖sinePi t‖ ^ 2 :=
    pow_le_pow_left₀ (by positivity) hs 2
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

theorem norm_ctKernel_top_le_fixed (x T : ℝ) (hT : 1 ≤ T) :
    ‖ctKernel27 ((x : ℂ) + (T : ℂ) * Complex.I)‖ ≤
      16 * Real.pi ^ 2 / Real.exp (Real.pi * T) ^ 2 :=
  norm_ctKernel_le_of_sine_ge_fixed (Real.exp_pos _)
    (norm_sinePi_top_ge_fixed x T hT)

theorem norm_ctKernel_bottom_le_fixed (x T : ℝ) (hT : 1 ≤ T) :
    ‖ctKernel27 ((x : ℂ) - (T : ℂ) * Complex.I)‖ ≤
      16 * Real.pi ^ 2 / Real.exp (Real.pi * T) ^ 2 :=
  norm_ctKernel_le_of_sine_ge_fixed (Real.exp_pos _)
    (norm_sinePi_bottom_ge_fixed x T hT)

/-! ## Polynomial growth of the rational factor -/

theorem abs_strip_sub_nat_le_fixed {n m r : ℕ}
    (hm1 : 1 ≤ m) (hmn : m ≤ n) (hr : r < n)
    {x : ℝ}
    (hx : x ∈ Set.Icc ((m : ℝ) - 1 / 2) ((m : ℝ) + 1 / 2)) :
    |x - ((r + 1 : ℕ) : ℝ)| ≤ (n : ℝ) + 1 := by
  have hrn : r + 1 ≤ n := by omega
  have hmR : (m : ℝ) ≤ n := by exact_mod_cast hmn
  have hrR : ((r + 1 : ℕ) : ℝ) ≤ n := by exact_mod_cast hrn
  have hr1 : (1 : ℝ) ≤ (r + 1 : ℕ) := by
    exact_mod_cast Nat.succ_le_succ (Nat.zero_le r)
  have hm1R : (1 : ℝ) ≤ m := by exact_mod_cast hm1
  rw [abs_le]
  constructor <;> linarith [hx.1, hx.2]

theorem norm_horizontal_factor_le_fixed {n m r : ℕ}
    (hm1 : 1 ≤ m) (hmn : m ≤ n) (hr : r < n)
    {x y : ℝ}
    (hx : x ∈ Set.Icc ((m : ℝ) - 1 / 2) ((m : ℝ) + 1 / 2)) :
    ‖(x : ℂ) + (y : ℂ) * Complex.I - (((r + 1 : ℕ) : ℂ))‖ ≤
      ((n : ℝ) + 2) * (1 + |y|) := by
  have habs := abs_strip_sub_nat_le_fixed hm1 hmn hr hx
  have heq :
      (x : ℂ) + (y : ℂ) * Complex.I - (((r + 1 : ℕ) : ℂ)) =
        ((x - ((r + 1 : ℕ) : ℝ) : ℝ) : ℂ) +
          (y : ℂ) * Complex.I := by
    push_cast
    ring
  rw [heq]
  have hreal :
      ‖((x - ((r + 1 : ℕ) : ℝ) : ℝ) : ℂ)‖ =
        |x - ((r + 1 : ℕ) : ℝ)| := by
    rw [Complex.norm_real, Real.norm_eq_abs]
  have himag : ‖(y : ℂ) * Complex.I‖ = |y| := by
    rw [norm_mul, Complex.norm_real, Real.norm_eq_abs]
    simp
  calc
    ‖((x - ((r + 1 : ℕ) : ℝ) : ℝ) : ℂ) +
        (y : ℂ) * Complex.I‖ ≤
      ‖((x - ((r + 1 : ℕ) : ℝ) : ℝ) : ℂ)‖ +
        ‖(y : ℂ) * Complex.I‖ := norm_add_le _ _
    _ = |x - ((r + 1 : ℕ) : ℝ)| + |y| := by rw [hreal, himag]
    _ ≤ ((n : ℝ) + 1) + |y| := by linarith
    _ ≤ ((n : ℝ) + 2) * (1 + |y|) := by
      have hn : 0 ≤ (n : ℝ) := Nat.cast_nonneg n
      have hy : 0 ≤ |y| := abs_nonneg y
      nlinarith

theorem norm_ctNumerator_on_strip_le_fixed {n m : ℕ}
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
        exact pow_le_pow_left₀ (norm_nonneg _)
          (norm_horizontal_factor_le_fixed hm1 hmn (Finset.mem_range.mp hr) hx) 3
    _ = ((((n : ℝ) + 2) * (1 + |y|)) ^ 3) ^ n := by simp

theorem half_le_norm_pole_factor_fixed {m j : ℕ}
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

theorem norm_ctDen_on_strip_ge_fixed {n m : ℕ}
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
          half_le_norm_pole_factor_fixed hm1 ht
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

theorem norm_ctR_on_strip_le_fixed {n m : ℕ}
    (hm1 : 1 ≤ m) (hmn : m ≤ n)
    {x y : ℝ}
    (hx : x ∈ Set.Icc ((m : ℝ) - 1 / 2) ((m : ℝ) + 1 / 2)) :
    ‖ctR27 n ((x : ℂ) + (y : ℂ) * Complex.I)‖ ≤
      ((((n : ℝ) + 2) * (1 + |y|)) ^ 3) ^ n /
        (1 / 2 : ℝ) ^ (n + 1) := by
  let t : ℂ := (x : ℂ) + (y : ℂ) * Complex.I
  have ht : t ∈ halfIntegerStrip (m : ℤ) := by
    simpa [t, halfIntegerStrip] using hx
  have hnum := norm_ctNumerator_on_strip_le_fixed hm1 hmn (y := y) hx
  have hden := norm_ctDen_on_strip_ge_fixed (n := n) hm1 ht
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

def ctHorizontalMajorant27Fixed (n : ℕ) (T : ℝ) : ℝ :=
  ((2 * ((n : ℝ) + 2)) ^ (3 * n) /
      (1 / 2 : ℝ) ^ (n + 1)) *
    (16 * Real.pi ^ 2 / Real.exp (Real.pi * T) ^ 2) *
    T ^ (3 * n)

theorem horizontal_poly_le_fixed {n : ℕ} {T : ℝ} (hT : 1 ≤ T) :
    ((((n : ℝ) + 2) * (1 + |T|)) ^ 3) ^ n ≤
      (2 * ((n : ℝ) + 2)) ^ (3 * n) * T ^ (3 * n) := by
  have hTabs : |T| = T := abs_of_nonneg (zero_le_one.trans hT)
  have hfactor : ((n : ℝ) + 2) * (1 + T) ≤
      (2 * ((n : ℝ) + 2)) * T := by
    have hn2 : 0 ≤ (n : ℝ) + 2 := by positivity
    nlinarith
  rw [hTabs]
  calc
    ((((n : ℝ) + 2) * (1 + T)) ^ 3) ^ n =
        (((n : ℝ) + 2) * (1 + T)) ^ (3 * n) :=
      (pow_mul _ 3 n).symm
    _ ≤ ((2 * ((n : ℝ) + 2)) * T) ^ (3 * n) :=
      pow_le_pow_left₀ (by positivity) hfactor (3 * n)
    _ = (2 * ((n : ℝ) + 2)) ^ (3 * n) * T ^ (3 * n) := by
      rw [mul_pow]

theorem norm_ctIntegrand_top_le_fixed {n m : ℕ}
    (hm1 : 1 ≤ m) (hmn : m ≤ n)
    {x T : ℝ}
    (hx : x ∈ Set.Icc ((m : ℝ) - 1 / 2) ((m : ℝ) + 1 / 2))
    (hT : 1 ≤ T) :
    ‖ctIntegrand27 n ((x : ℂ) + (T : ℂ) * Complex.I)‖ ≤
      ctHorizontalMajorant27Fixed n T := by
  have hR := norm_ctR_on_strip_le_fixed hm1 hmn (y := T) hx
  have hK := norm_ctKernel_top_le_fixed x T hT
  have hpoly := horizontal_poly_le_fixed (n := n) hT
  have hhalf : 0 ≤ (1 / 2 : ℝ) ^ (n + 1) := by positivity
  have hfrac :
      ((((n : ℝ) + 2) * (1 + |T|)) ^ 3) ^ n /
          (1 / 2 : ℝ) ^ (n + 1) ≤
        ((2 * ((n : ℝ) + 2)) ^ (3 * n) * T ^ (3 * n)) /
          (1 / 2 : ℝ) ^ (n + 1) :=
    div_le_div_of_nonneg_right hpoly hhalf
  rw [ctIntegrand27, norm_mul]
  calc
    ‖ctR27 n ((x : ℂ) + (T : ℂ) * Complex.I)‖ *
        ‖ctKernel27 ((x : ℂ) + (T : ℂ) * Complex.I)‖ ≤
      (((((n : ℝ) + 2) * (1 + |T|)) ^ 3) ^ n /
          (1 / 2 : ℝ) ^ (n + 1)) *
        (16 * Real.pi ^ 2 / Real.exp (Real.pi * T) ^ 2) := by
      exact mul_le_mul hR hK (norm_nonneg _) (by positivity)
    _ ≤ (((2 * ((n : ℝ) + 2)) ^ (3 * n) * T ^ (3 * n)) /
          (1 / 2 : ℝ) ^ (n + 1)) *
        (16 * Real.pi ^ 2 / Real.exp (Real.pi * T) ^ 2) := by
      exact mul_le_mul_of_nonneg_right hfrac (by positivity)
    _ = ctHorizontalMajorant27Fixed n T := by
      unfold ctHorizontalMajorant27Fixed
      ring

theorem norm_ctIntegrand_bottom_le_fixed {n m : ℕ}
    (hm1 : 1 ≤ m) (hmn : m ≤ n)
    {x T : ℝ}
    (hx : x ∈ Set.Icc ((m : ℝ) - 1 / 2) ((m : ℝ) + 1 / 2))
    (hT : 1 ≤ T) :
    ‖ctIntegrand27 n ((x : ℂ) - (T : ℂ) * Complex.I)‖ ≤
      ctHorizontalMajorant27Fixed n T := by
  have hpoint :
      (x : ℂ) - (T : ℂ) * Complex.I =
        (x : ℂ) + ((-T : ℝ) : ℂ) * Complex.I := by
    apply Complex.ext <;> simp
  rw [hpoint]
  have hR := norm_ctR_on_strip_le_fixed hm1 hmn (y := -T) hx
  have hK :
      ‖ctKernel27 ((x : ℂ) + ((-T : ℝ) : ℂ) * Complex.I)‖ ≤
        16 * Real.pi ^ 2 / Real.exp (Real.pi * T) ^ 2 := by
    rw [← hpoint]
    exact norm_ctKernel_bottom_le_fixed x T hT
  have hpoly :
      ((((n : ℝ) + 2) * (1 + |-T|)) ^ 3) ^ n ≤
        (2 * ((n : ℝ) + 2)) ^ (3 * n) * T ^ (3 * n) := by
    simpa using horizontal_poly_le_fixed (n := n) hT
  have hhalf : 0 ≤ (1 / 2 : ℝ) ^ (n + 1) := by positivity
  have hfrac :
      ((((n : ℝ) + 2) * (1 + |-T|)) ^ 3) ^ n /
          (1 / 2 : ℝ) ^ (n + 1) ≤
        ((2 * ((n : ℝ) + 2)) ^ (3 * n) * T ^ (3 * n)) /
          (1 / 2 : ℝ) ^ (n + 1) :=
    div_le_div_of_nonneg_right hpoly hhalf
  rw [ctIntegrand27, norm_mul]
  calc
    ‖ctR27 n ((x : ℂ) + ((-T : ℝ) : ℂ) * Complex.I)‖ *
        ‖ctKernel27 ((x : ℂ) + ((-T : ℝ) : ℂ) * Complex.I)‖ ≤
      (((((n : ℝ) + 2) * (1 + |-T|)) ^ 3) ^ n /
          (1 / 2 : ℝ) ^ (n + 1)) *
        (16 * Real.pi ^ 2 / Real.exp (Real.pi * T) ^ 2) := by
      exact mul_le_mul hR hK (norm_nonneg _) (by positivity)
    _ ≤ (((2 * ((n : ℝ) + 2)) ^ (3 * n) * T ^ (3 * n)) /
          (1 / 2 : ℝ) ^ (n + 1)) *
        (16 * Real.pi ^ 2 / Real.exp (Real.pi * T) ^ 2) := by
      exact mul_le_mul_of_nonneg_right hfrac (by positivity)
    _ = ctHorizontalMajorant27Fixed n T := by
      unfold ctHorizontalMajorant27Fixed
      ring

theorem tendsto_pow_mul_exp_neg_two_pi27_fixed (p : ℕ) :
    Tendsto (fun T : ℝ => T ^ p * Real.exp (-(2 * Real.pi * T)))
      atTop (𝓝 0) := by
  have hraw := (Real.tendsto_pow_mul_exp_neg_atTop_nhds_zero p).comp
    (tendsto_id.const_mul_atTop' Real.two_pi_pos)
  have hconst : (2 * Real.pi : ℝ) ^ p ≠ 0 :=
    pow_ne_zero _ (mul_ne_zero (by norm_num) Real.pi_ne_zero)
  have hc := hraw.const_mul ((2 * Real.pi : ℝ) ^ p)⁻¹
  convert hc using 1
  · funext T
    simp only [Function.comp_apply, id_eq]
    rw [mul_pow]
    field_simp [hconst]
    ring
  · simp

theorem ctHorizontalMajorant_tendsto_zero27_fixed (n : ℕ) :
    Tendsto (ctHorizontalMajorant27Fixed n) atTop (𝓝 0) := by
  have hbase := tendsto_pow_mul_exp_neg_two_pi27_fixed (3 * n)
  let C : ℝ :=
    ((2 * ((n : ℝ) + 2)) ^ (3 * n) /
      (1 / 2 : ℝ) ^ (n + 1)) * (16 * Real.pi ^ 2)
  have hc := hbase.const_mul C
  apply hc.congr'
  filter_upwards with T
  unfold ctHorizontalMajorant27Fixed
  dsimp [C]
  have hE : Real.exp (Real.pi * T) ^ 2 =
      Real.exp (2 * Real.pi * T) := by
    rw [pow_two, ← Real.exp_add]
    congr 1
    ring
  rw [hE, div_eq_mul_inv, ← Real.exp_neg]
  ring

theorem ctHorizontal_top_tendsto_zero27_fixed {n m : ℕ}
    (hm1 : 1 ≤ m) (hmn : m ≤ n) :
    Tendsto
      (fun T : ℝ => ∫ x in ((m : ℝ) - 1 / 2)..((m : ℝ) + 1 / 2),
        ctExtension27 n m ((x : ℂ) + (T : ℂ) * Complex.I))
      atTop (𝓝 0) := by
  rw [tendsto_zero_iff_norm_tendsto_zero]
  apply squeeze_zero'
  · exact Eventually.of_forall fun T => norm_nonneg _
  · filter_upwards [eventually_ge_atTop (1 : ℝ)] with T hT
    let a : ℝ := (m : ℝ) - 1 / 2
    let b : ℝ := (m : ℝ) + 1 / 2
    have hab : a ≤ b := by dsimp [a, b]; linarith
    have hbound : ∀ x ∈ Ι a b,
        ‖ctExtension27 n m ((x : ℂ) + (T : ℂ) * Complex.I)‖ ≤
          ctHorizontalMajorant27Fixed n T := by
      intro x hx
      have hxioc : x ∈ Set.Ioc a b := by
        simpa [uIoc_of_le hab] using hx
      have hx' : x ∈ Set.Icc ((m : ℝ) - 1 / 2) ((m : ℝ) + 1 / 2) := by
        have := Ioc_subset_Icc_self hxioc
        simpa [a, b] using this
      have hmem :
          (x : ℂ) + (T : ℂ) * Complex.I ∈ halfIntegerStrip (m : ℤ) := by
        simpa [halfIntegerStrip] using hx'
      have hne : (x : ℂ) + (T : ℂ) * Complex.I ≠ (m : ℂ) := by
        intro h
        have him := congrArg Complex.im h
        simp at him
        linarith
      rw [← ctIntegrand_eq_extension27 hm1 hmn hmem hne]
      exact norm_ctIntegrand_top_le_fixed hm1 hmn hx' hT
    have hi := intervalIntegral.norm_integral_le_of_norm_le_const hbound
    have hlen : |b - a| = 1 := by
      dsimp [a, b]
      rw [abs_of_nonneg (by linarith)]
      ring
    simpa [a, b, hlen] using hi
  · exact ctHorizontalMajorant_tendsto_zero27_fixed n

theorem ctHorizontal_bottom_tendsto_zero27_fixed {n m : ℕ}
    (hm1 : 1 ≤ m) (hmn : m ≤ n) :
    Tendsto
      (fun T : ℝ => ∫ x in ((m : ℝ) - 1 / 2)..((m : ℝ) + 1 / 2),
        ctExtension27 n m ((x : ℂ) - (T : ℂ) * Complex.I))
      atTop (𝓝 0) := by
  rw [tendsto_zero_iff_norm_tendsto_zero]
  apply squeeze_zero'
  · exact Eventually.of_forall fun T => norm_nonneg _
  · filter_upwards [eventually_ge_atTop (1 : ℝ)] with T hT
    let a : ℝ := (m : ℝ) - 1 / 2
    let b : ℝ := (m : ℝ) + 1 / 2
    have hab : a ≤ b := by dsimp [a, b]; linarith
    have hbound : ∀ x ∈ Ι a b,
        ‖ctExtension27 n m ((x : ℂ) - (T : ℂ) * Complex.I)‖ ≤
          ctHorizontalMajorant27Fixed n T := by
      intro x hx
      have hxioc : x ∈ Set.Ioc a b := by
        simpa [uIoc_of_le hab] using hx
      have hx' : x ∈ Set.Icc ((m : ℝ) - 1 / 2) ((m : ℝ) + 1 / 2) := by
        have := Ioc_subset_Icc_self hxioc
        simpa [a, b] using this
      have hmem :
          (x : ℂ) - (T : ℂ) * Complex.I ∈ halfIntegerStrip (m : ℤ) := by
        simpa [halfIntegerStrip] using hx'
      have hne : (x : ℂ) - (T : ℂ) * Complex.I ≠ (m : ℂ) := by
        intro h
        have him := congrArg Complex.im h
        simp at him
        linarith
      rw [← ctIntegrand_eq_extension27 hm1 hmn hmem hne]
      exact norm_ctIntegrand_bottom_le_fixed hm1 hmn hx' hT
    have hi := intervalIntegral.norm_integral_le_of_norm_le_const hbound
    have hlen : |b - a| = 1 := by
      dsimp [a, b]
      rw [abs_of_nonneg (by linarith)]
      ring
    simpa [a, b, hlen] using hi
  · exact ctHorizontalMajorant_tendsto_zero27_fixed n

#print axioms ctHorizontal_top_tendsto_zero27_fixed
#print axioms ctHorizontal_bottom_tendsto_zero27_fixed

end RamanujanChallenge.P27.Q6427
