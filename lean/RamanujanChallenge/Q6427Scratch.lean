import RamanujanChallenge.Problem27Barnes
import Mathlib.Analysis.Complex.CauchyIntegral
import Mathlib.Analysis.Complex.RemovableSingularity
import Mathlib.Analysis.SpecialFunctions.Gamma.Basic
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Complex
import Mathlib.MeasureTheory.Integral.IntegralEqImproper

open Filter Set Topology MeasureTheory
open scoped Interval Real BigOperators

noncomputable section

namespace RamanujanChallenge.P27.Q6427

/-! ## Generic rectangle and sine extension, specialized below -/

def closedVerticalStrip (a b : ℝ) : Set ℂ :=
  {z | a ≤ z.re ∧ z.re ≤ b}

def verticalPoint (x y : ℝ) : ℂ :=
  (x : ℂ) + (y : ℂ) * Complex.I

theorem verticalIntegral_eq_of_horizontal_tendsto
    {F : ℂ → ℂ} {a b : ℝ} (hab : a ≤ b)
    (hF : DifferentiableOn ℂ F (closedVerticalStrip a b))
    (hleft : Integrable (fun y : ℝ => F (verticalPoint a y)))
    (hright : Integrable (fun y : ℝ => F (verticalPoint b y)))
    (htop : Tendsto
      (fun T : ℝ => ∫ x in a..b,
        F ((x : ℂ) + (T : ℂ) * Complex.I))
      atTop (𝓝 0))
    (hbottom : Tendsto
      (fun T : ℝ => ∫ x in a..b,
        F ((x : ℂ) - (T : ℂ) * Complex.I))
      atTop (𝓝 0)) :
    (∫ y : ℝ, F (verticalPoint a y)) =
      ∫ y : ℝ, F (verticalPoint b y) := by
  have hfinite : ∀ T : ℝ, 0 ≤ T →
      (∫ y in -T..T, F (verticalPoint a y)) -
          (∫ y in -T..T, F (verticalPoint b y)) =
        Complex.I *
          ((∫ x in a..b,
              F ((x : ℂ) + (T : ℂ) * Complex.I)) -
           (∫ x in a..b,
              F ((x : ℂ) - (T : ℂ) * Complex.I))) := by
    intro T hT
    let z : ℂ := (a : ℂ) - (T : ℂ) * Complex.I
    let w : ℂ := (b : ℂ) + (T : ℂ) * Complex.I
    have hrect :=
      Complex.integral_boundary_rect_eq_zero_of_differentiableOn
        F z w (hF.mono (by
          intro u hu
          have hre : u.re ∈ [[a, b]] := by
            simpa [z, w] using hu.1
          simpa [closedVerticalStrip, uIcc_of_le hab] using hre))
    have hrect' :
        (∫ x in a..b,
            F ((x : ℂ) - (T : ℂ) * Complex.I)) -
          (∫ x in a..b,
            F ((x : ℂ) + (T : ℂ) * Complex.I)) +
          Complex.I * (∫ y in -T..T, F (verticalPoint b y)) -
          Complex.I * (∫ y in -T..T, F (verticalPoint a y)) = 0 := by
      simpa [z, w, verticalPoint, smul_eq_mul] using hrect
    have hIR :
        Complex.I *
            ((∫ y in -T..T, F (verticalPoint b y)) -
              (∫ y in -T..T, F (verticalPoint a y))) =
          (∫ x in a..b,
              F ((x : ℂ) + (T : ℂ) * Complex.I)) -
            (∫ x in a..b,
              F ((x : ℂ) - (T : ℂ) * Complex.I)) := by
      linear_combination hrect'
    calc
      (∫ y in -T..T, F (verticalPoint a y)) -
          (∫ y in -T..T, F (verticalPoint b y)) =
          -((∫ y in -T..T, F (verticalPoint b y)) -
            (∫ y in -T..T, F (verticalPoint a y))) := by ring
      _ = Complex.I ^ 2 *
          ((∫ y in -T..T, F (verticalPoint b y)) -
            (∫ y in -T..T, F (verticalPoint a y))) := by
        rw [Complex.I_sq]
        ring
      _ = Complex.I *
          (Complex.I *
            ((∫ y in -T..T, F (verticalPoint b y)) -
              (∫ y in -T..T, F (verticalPoint a y)))) := by ring
      _ = Complex.I *
          ((∫ x in a..b,
              F ((x : ℂ) + (T : ℂ) * Complex.I)) -
            (∫ x in a..b,
              F ((x : ℂ) - (T : ℂ) * Complex.I))) := by rw [hIR]
  have hleft_lim : Tendsto
      (fun T : ℝ => ∫ y in -T..T, F (verticalPoint a y))
      atTop (𝓝 (∫ y : ℝ, F (verticalPoint a y))) :=
    MeasureTheory.intervalIntegral_tendsto_integral
      hleft tendsto_neg_atTop_atBot tendsto_id
  have hright_lim : Tendsto
      (fun T : ℝ => ∫ y in -T..T, F (verticalPoint b y))
      atTop (𝓝 (∫ y : ℝ, F (verticalPoint b y))) :=
    MeasureTheory.intervalIntegral_tendsto_integral
      hright tendsto_neg_atTop_atBot tendsto_id
  have hfinite_eventually : ∀ᶠ T : ℝ in atTop,
      (∫ y in -T..T, F (verticalPoint a y)) -
          (∫ y in -T..T, F (verticalPoint b y)) =
        Complex.I *
          ((∫ x in a..b,
              F ((x : ℂ) + (T : ℂ) * Complex.I)) -
           (∫ x in a..b,
              F ((x : ℂ) - (T : ℂ) * Complex.I))) :=
    (eventually_ge_atTop (0 : ℝ)).mono hfinite
  have hlhs := hleft_lim.sub hright_lim
  have hrhs : Tendsto
      (fun T : ℝ => Complex.I *
        ((∫ x in a..b,
            F ((x : ℂ) + (T : ℂ) * Complex.I)) -
         (∫ x in a..b,
            F ((x : ℂ) - (T : ℂ) * Complex.I))))
      atTop (𝓝 0) := by
    simpa using tendsto_const_nhds.mul (htop.sub hbottom)
  have hlhs_zero := hrhs.congr'
    (hfinite_eventually.mono fun _ h => h.symm)
  have hsub :
      (∫ y : ℝ, F (verticalPoint a y)) -
          (∫ y : ℝ, F (verticalPoint b y)) = 0 :=
    tendsto_nhds_unique hlhs hlhs_zero
  exact sub_eq_zero.mp hsub

def sinePi (z : ℂ) : ℂ :=
  Complex.sin ((Real.pi : ℂ) * z)

def sineSlope (m : ℤ) : ℂ → ℂ :=
  dslope sinePi (m : ℂ)

@[simp] theorem sinePi_int (m : ℤ) : sinePi (m : ℂ) = 0 := by
  rw [sinePi, Complex.sin_eq_zero_iff]
  exact ⟨m, by push_cast; ring⟩

theorem sinePi_eq_sub_mul_sineSlope (m : ℤ) (z : ℂ) :
    sinePi z = (z - (m : ℂ)) * sineSlope m z := by
  have h := sub_smul_dslope sinePi (m : ℂ) z
  simpa [sineSlope, smul_eq_mul] using h.symm

theorem sinePi_differentiable : Differentiable ℂ sinePi := by
  intro z
  simpa [sinePi] using
    (((hasDerivAt_id z).const_mul (Real.pi : ℂ)).csin.differentiableAt)

theorem sineSlope_differentiable (m : ℤ) : Differentiable ℂ (sineSlope m) := by
  rw [← differentiableOn_univ]
  exact (Complex.differentiableOn_dslope
      (f := sinePi) (s := Set.univ) (c := (m : ℂ)) univ_mem).2
    sinePi_differentiable.differentiableOn

theorem sineSlope_at_int_ne_zero (m : ℤ) : sineSlope m (m : ℂ) ≠ 0 := by
  have hderiv : deriv sinePi (m : ℂ) =
      (Real.pi : ℂ) * Complex.cos ((Real.pi : ℂ) * (m : ℂ)) := by
    have hinner : HasDerivAt (fun z : ℂ => (Real.pi : ℂ) * z)
        (Real.pi : ℂ) (m : ℂ) := by
      simpa only [id_eq, mul_comm, mul_one] using
        ((hasDerivAt_id (m : ℂ)).const_mul (Real.pi : ℂ))
    have h := hinner.csin.deriv
    rw [show sinePi =
      (fun z : ℂ => Complex.sin ((Real.pi : ℂ) * z)) by rfl]
    simpa only [mul_comm] using h
  have harg : (Real.pi : ℂ) * (m : ℂ) =
      (((m : ℝ) * Real.pi : ℝ) : ℂ) := by push_cast; ring
  have hcos : Complex.cos ((Real.pi : ℂ) * (m : ℂ)) =
      ((((-1 : ℝ) ^ m : ℝ)) : ℂ) := by
    rw [harg, ← Complex.ofReal_cos, Real.cos_int_mul_pi]
  have hpow : ((-1 : ℝ) ^ m) ≠ 0 := by
    intro hzero
    have habs := Real.abs_cos_int_mul_pi m
    rw [Real.cos_int_mul_pi, hzero, abs_zero] at habs
    norm_num at habs
  have hcos0 : Complex.cos ((Real.pi : ℂ) * (m : ℂ)) ≠ 0 := by
    rw [hcos]
    exact Complex.ofReal_ne_zero.mpr hpow
  rw [sineSlope, dslope_same, hderiv]
  exact mul_ne_zero (Complex.ofReal_ne_zero.mpr Real.pi_ne_zero) hcos0

def halfIntegerStrip (m : ℤ) : Set ℂ :=
  {z | (m : ℝ) - 1 / 2 ≤ z.re ∧ z.re ≤ (m : ℝ) + 1 / 2}

theorem sineSlope_ne_zero_on_strip (m : ℤ) {z : ℂ}
    (hz : z ∈ halfIntegerStrip m) : sineSlope m z ≠ 0 := by
  by_cases hzm : z = (m : ℂ)
  · subst z
    exact sineSlope_at_int_ne_zero m
  intro hslope
  have hsin : sinePi z = 0 := by
    rw [sinePi_eq_sub_mul_sineSlope m z, hslope, mul_zero]
  rcases Complex.sin_eq_zero_iff.mp (by simpa [sinePi] using hsin) with ⟨k, hk⟩
  have hzk : z = (k : ℂ) := by
    calc
      z = (Real.pi : ℂ)⁻¹ * ((Real.pi : ℂ) * z) := by field_simp
      _ = (Real.pi : ℂ)⁻¹ * ((k : ℂ) * (Real.pi : ℂ)) := by rw [hk]
      _ = (k : ℂ) := by field_simp
  subst z
  have hl : (m : ℝ) - 1 / 2 ≤ (k : ℝ) := by
    simpa [halfIntegerStrip] using hz.1
  have hu : (k : ℝ) ≤ (m : ℝ) + 1 / 2 := by
    simpa [halfIntegerStrip] using hz.2
  have hlowR : (-1 : ℝ) < ((k - m : ℤ) : ℝ) := by push_cast; linarith
  have huppR : ((k - m : ℤ) : ℝ) < (1 : ℝ) := by push_cast; linarith
  have hlowZ : (-1 : ℤ) < k - m := by exact_mod_cast hlowR
  have huppZ : k - m < (1 : ℤ) := by exact_mod_cast huppR
  have hkm : k = m := by omega
  apply hzm
  simpa [hkm]

/-! ## Direct t-coordinate rational function -/

def ctNumerator27 (n : ℕ) (t : ℂ) : ℂ :=
  ∏ r in Finset.range n, (t - (((r + 1 : ℕ) : ℂ))) ^ 3

def ctPoleProduct27 (q : ℕ) (t : ℂ) : ℂ :=
  ∏ j in Finset.range q, (t + (j : ℂ))

def ctR27 (n : ℕ) (t : ℂ) : ℂ :=
  ctNumerator27 n t /
    (((n.factorial : ℂ) ^ 2) * ctPoleProduct27 (n + 1) t)

def ctKernel27 (t : ℂ) : ℂ :=
  ((Real.pi : ℂ) / sinePi t) ^ 2

def ctIntegrand27 (n : ℕ) (t : ℂ) : ℂ :=
  ctR27 n t * ctKernel27 t

def ctNumeratorWithout27 (n m : ℕ) (t : ℂ) : ℂ :=
  ∏ r in (Finset.range n).erase (m - 1),
    (t - (((r + 1 : ℕ) : ℂ))) ^ 3

def ctRemoved27 (n m : ℕ) (t : ℂ) : ℂ :=
  ctNumeratorWithout27 n m t /
    (((n.factorial : ℂ) ^ 2) * ctPoleProduct27 (n + 1) t)

def ctExtension27 (n m : ℕ) (t : ℂ) : ℂ :=
  (Real.pi : ℂ) ^ 2 * (t - (m : ℂ)) * ctRemoved27 n m t /
    sineSlope (m : ℤ) t ^ 2

theorem ctNumerator_factor27 {n m : ℕ} (hm1 : 1 ≤ m) (hmn : m ≤ n)
    (t : ℂ) :
    ctNumerator27 n t =
      (t - (m : ℂ)) ^ 3 * ctNumeratorWithout27 n m t := by
  have hmem : m - 1 ∈ Finset.range n := by
    rw [Finset.mem_range]
    omega
  unfold ctNumerator27 ctNumeratorWithout27
  rw [← Finset.mul_prod_erase _ _ hmem]
  congr 2
  push_cast
  omega

theorem ctPoleProduct_ne_zero_on_strip27 {n m : ℕ} (hm1 : 1 ≤ m)
    {t : ℂ} (ht : t ∈ halfIntegerStrip (m : ℤ)) :
    ctPoleProduct27 (n + 1) t ≠ 0 := by
  unfold ctPoleProduct27
  rw [Finset.prod_ne_zero_iff]
  intro j hj hzero
  have hre := congrArg Complex.re hzero
  simp at hre
  have hj0 : (0 : ℝ) ≤ j := Nat.cast_nonneg j
  have hmR : (1 : ℝ) ≤ m := by exact_mod_cast hm1
  have hleft := ht.1
  norm_num [halfIntegerStrip] at hleft
  linarith

theorem ctNumeratorWithout_differentiableAt27 (n m : ℕ) (t : ℂ) :
    DifferentiableAt ℂ (ctNumeratorWithout27 n m) t := by
  unfold ctNumeratorWithout27
  exact DifferentiableAt.fun_finset_prod
    (u := (Finset.range n).erase (m - 1)) fun r _ =>
      (differentiableAt_id.sub_const _).fun_pow 3

theorem ctPoleProduct_differentiableAt27 (q : ℕ) (t : ℂ) :
    DifferentiableAt ℂ (ctPoleProduct27 q) t := by
  unfold ctPoleProduct27
  exact DifferentiableAt.fun_finset_prod
    (u := Finset.range q) fun j _ => differentiableAt_id.add_const _

theorem ctRemoved_differentiableOn27 {n m : ℕ} (hm1 : 1 ≤ m) :
    DifferentiableOn ℂ (ctRemoved27 n m) (halfIntegerStrip (m : ℤ)) := by
  intro t ht
  unfold ctRemoved27
  apply DifferentiableAt.differentiableWithinAt
  apply DifferentiableAt.div
  · exact ctNumeratorWithout_differentiableAt27 n m t
  · exact (differentiableAt_const _).mul
      (ctPoleProduct_differentiableAt27 (n + 1) t)
  · apply mul_ne_zero
    · exact pow_ne_zero _ (Nat.cast_ne_zero.mpr (Nat.factorial_ne_zero n))
    · exact ctPoleProduct_ne_zero_on_strip27 hm1 ht

theorem ctExtension_differentiableOn27 {n m : ℕ} (hm1 : 1 ≤ m) :
    DifferentiableOn ℂ (ctExtension27 n m) (halfIntegerStrip (m : ℤ)) := by
  have hnum : DifferentiableOn ℂ
      (fun t => (Real.pi : ℂ) ^ 2 * (t - (m : ℂ)) * ctRemoved27 n m t)
      (halfIntegerStrip (m : ℤ)) :=
    ((differentiableOn_const ((Real.pi : ℂ) ^ 2)).mul
      (differentiableOn_id.sub_const (m : ℂ))).mul
        (ctRemoved_differentiableOn27 (n := n) hm1)
  have hden : DifferentiableOn ℂ (fun t => sineSlope (m : ℤ) t ^ 2)
      (halfIntegerStrip (m : ℤ)) :=
    (sineSlope_differentiable (m : ℤ)).differentiableOn.pow 2
  simpa [ctExtension27] using hnum.div hden
    (fun t ht => pow_ne_zero 2 (sineSlope_ne_zero_on_strip (m : ℤ) ht))

theorem ctIntegrand_eq_extension27 {n m : ℕ} (hm1 : 1 ≤ m) (hmn : m ≤ n)
    {t : ℂ} (ht : t ∈ halfIntegerStrip (m : ℤ)) (htm : t ≠ (m : ℂ)) :
    ctIntegrand27 n t = ctExtension27 n m t := by
  have hs := sineSlope_ne_zero_on_strip (m : ℤ) ht
  rw [ctIntegrand27, ctR27, ctKernel27, ctNumerator_factor27 hm1 hmn,
    sinePi_eq_sub_mul_sineSlope (m : ℤ) t, ctExtension27, ctRemoved27]
  field_simp [sub_ne_zero.mpr htm, hs]

/-- Structural one-strip theorem. The two vertical integrability and two
horizontal decay hypotheses will be discharged by the explicit estimates below. -/
theorem ctR_one_strip_shift27
    {n m : ℕ} (hm1 : 1 ≤ m) (hmn : m ≤ n)
    (hleft : Integrable (fun y : ℝ => ctIntegrand27 n
      (verticalPoint ((m : ℝ) - 1 / 2) y)))
    (hright : Integrable (fun y : ℝ => ctIntegrand27 n
      (verticalPoint ((m : ℝ) + 1 / 2) y)))
    (htop : Tendsto
      (fun T : ℝ => ∫ x in ((m : ℝ) - 1 / 2)..((m : ℝ) + 1 / 2),
        ctExtension27 n m ((x : ℂ) + (T : ℂ) * Complex.I))
      atTop (𝓝 0))
    (hbottom : Tendsto
      (fun T : ℝ => ∫ x in ((m : ℝ) - 1 / 2)..((m : ℝ) + 1 / 2),
        ctExtension27 n m ((x : ℂ) - (T : ℂ) * Complex.I))
      atTop (𝓝 0)) :
    (∫ y : ℝ, ctIntegrand27 n
      (verticalPoint ((m : ℝ) - 1 / 2) y)) =
    ∫ y : ℝ, ctIntegrand27 n
      (verticalPoint ((m : ℝ) + 1 / 2) y) := by
  let a : ℝ := (m : ℝ) - 1 / 2
  let b : ℝ := (m : ℝ) + 1 / 2
  have hab : a ≤ b := by dsimp [a, b]; linarith
  have hstrip : closedVerticalStrip a b = halfIntegerStrip (m : ℤ) := by
    ext z
    simp [closedVerticalStrip, halfIntegerStrip, a, b]
  have hleftEq :
      (fun y : ℝ => ctExtension27 n m (verticalPoint a y)) =
        fun y : ℝ => ctIntegrand27 n (verticalPoint a y) := by
    funext y
    symm
    apply ctIntegrand_eq_extension27 hm1 hmn
    · simp [halfIntegerStrip, verticalPoint, a] <;> linarith
    · intro h
      have hre := congrArg Complex.re h
      simp [verticalPoint, a] at hre
  have hrightEq :
      (fun y : ℝ => ctExtension27 n m (verticalPoint b y)) =
        fun y : ℝ => ctIntegrand27 n (verticalPoint b y) := by
    funext y
    symm
    apply ctIntegrand_eq_extension27 hm1 hmn
    · simp [halfIntegerStrip, verticalPoint, b] <;> linarith
    · intro h
      have hre := congrArg Complex.re h
      simp [verticalPoint, b] at hre
  have hleftExt : Integrable
      (fun y : ℝ => ctExtension27 n m (verticalPoint a y)) := by
    rw [hleftEq]
    simpa [a] using hleft
  have hrightExt : Integrable
      (fun y : ℝ => ctExtension27 n m (verticalPoint b y)) := by
    rw [hrightEq]
    simpa [b] using hright
  have hshift := verticalIntegral_eq_of_horizontal_tendsto
    (F := ctExtension27 n m) hab
    (by rw [hstrip]; exact ctExtension_differentiableOn27 hm1)
    hleftExt hrightExt
    (by simpa [a, b] using htop)
    (by simpa [a, b] using hbottom)
  rw [hleftEq, hrightEq] at hshift
  simpa [a, b] using hshift

#print axioms ctR_one_strip_shift27

end RamanujanChallenge.P27.Q6427
