import RamanujanChallenge.Problem27Barnes
import Mathlib.Analysis.Complex.CauchyIntegral
import Mathlib.Analysis.Complex.RemovableSingularity
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Complex
import Mathlib.MeasureTheory.Integral.IntegralEqImproper

open Filter Set Topology
open scoped Interval Real

noncomputable section

namespace RamanujanChallenge.P27.Q6382

private def closedVerticalStrip (a b : ℝ) : Set ℂ :=
  {z | a ≤ z.re ∧ z.re ≤ b}

private def verticalPoint (x y : ℝ) : ℂ :=
  (x : ℂ) + (y : ℂ) * Complex.I

/-- Rectangle Cauchy plus vanishing horizontal sides. These are the
parameterized vertical integrals; multiplying both sides by `I` gives the
usual oriented line integrals. -/
theorem verticalIntegral_eq_of_horizontal_tendsto
    {F : ℂ → ℂ} {a b : ℝ} (hab : a ≤ b)
    (hF : DifferentiableOn ℂ F (closedVerticalStrip a b))
    (hleft : MeasureTheory.Integrable
      (fun y : ℝ => F (verticalPoint a y)))
    (hright : MeasureTheory.Integrable
      (fun y : ℝ => F (verticalPoint b y)))
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
    linear_combination Complex.I * hrect'

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
            F ((x : ℂ) - (T : ℂ) * Complex.I)))
      atTop (𝓝 0) := by
    simpa using tendsto_const_nhds.mul (htop.sub hbottom)
  have hlhs_zero := hrhs.congr' hfinite_eventually.symm
  have hsub :
      (∫ y : ℝ, F (verticalPoint a y)) -
          (∫ y : ℝ, F (verticalPoint b y)) = 0 :=
    tendsto_nhds_unique hlhs hlhs_zero
  exact sub_eq_zero.mp hsub

private def sinePi (z : ℂ) : ℂ :=
  Complex.sin ((Real.pi : ℂ) * z)

/-- Holomorphic extension of `sin (π z) / (z-m)` at the integer `m`. -/
private def sineSlope (m : ℤ) : ℂ → ℂ :=
  dslope sinePi (m : ℂ)

@[simp] private theorem sinePi_int (m : ℤ) :
    sinePi (m : ℂ) = 0 := by
  rw [sinePi, Complex.sin_eq_zero_iff]
  exact ⟨m, by push_cast; ring⟩

private theorem sinePi_eq_sub_mul_sineSlope (m : ℤ) (z : ℂ) :
    sinePi z = (z - (m : ℂ)) * sineSlope m z := by
  have h := sub_smul_dslope sinePi (m : ℂ) z
  simpa [sineSlope, smul_eq_mul] using h.symm

private theorem sinePi_differentiable : Differentiable ℂ sinePi := by
  intro z
  simpa [sinePi] using
    (((hasDerivAt_id z).const_mul (Real.pi : ℂ)).csin.differentiableAt)

private theorem sineSlope_differentiable (m : ℤ) :
    Differentiable ℂ (sineSlope m) := by
  rw [← differentiableOn_univ]
  exact (Complex.differentiableOn_dslope
      (f := sinePi) (s := Set.univ) (c := (m : ℂ)) univ_mem).2
    sinePi_differentiable.differentiableOn

private theorem sineSlope_at_int_ne_zero (m : ℤ) :
    sineSlope m (m : ℂ) ≠ 0 := by
  have hderiv :
      deriv sinePi (m : ℂ) =
        (Real.pi : ℂ) *
          Complex.cos ((Real.pi : ℂ) * (m : ℂ)) := by
    simpa [sinePi, mul_comm] using
      (((hasDerivAt_id (m : ℂ)).const_mul (Real.pi : ℂ)).csin.deriv)
  have harg :
      (Real.pi : ℂ) * (m : ℂ) =
        (((m : ℝ) * Real.pi : ℝ) : ℂ) := by
    push_cast
    ring
  have hcos :
      Complex.cos ((Real.pi : ℂ) * (m : ℂ)) =
        ((((-1 : ℝ) ^ m : ℝ)) : ℂ) := by
    rw [harg, ← Complex.ofReal_cos, Real.cos_int_mul_pi]
  rw [sineSlope, dslope_same, hderiv, hcos]
  exact mul_ne_zero (Complex.ofReal_ne_zero.mpr Real.pi_ne_zero) (by simp)

private def halfIntegerStrip (m : ℤ) : Set ℂ :=
  {z | (m : ℝ) - 1 / 2 ≤ z.re ∧ z.re ≤ (m : ℝ) + 1 / 2}

private theorem sineSlope_ne_zero_on_strip
    (m : ℤ) {z : ℂ} (hz : z ∈ halfIntegerStrip m) :
    sineSlope m z ≠ 0 := by
  by_cases hzm : z = (m : ℂ)
  · subst z
    exact sineSlope_at_int_ne_zero m
  intro hslope
  have hsin : sinePi z = 0 := by
    rw [sinePi_eq_sub_mul_sineSlope m z, hslope, mul_zero]
  rcases Complex.sin_eq_zero_iff.mp (by simpa [sinePi] using hsin) with ⟨k, hk⟩
  have hpi0 : (Real.pi : ℂ) ≠ 0 :=
    Complex.ofReal_ne_zero.mpr Real.pi_ne_zero
  have hzk : z = (k : ℂ) := by
    calc
      z = (Real.pi : ℂ)⁻¹ * ((Real.pi : ℂ) * z) := by
        field_simp
      _ = (Real.pi : ℂ)⁻¹ * ((k : ℂ) * (Real.pi : ℂ)) := by
        rw [hk]
      _ = (k : ℂ) := by
        field_simp
  subst z
  have hl : (m : ℝ) - 1 / 2 ≤ (k : ℝ) := by
    simpa [halfIntegerStrip] using hz.1
  have hu : (k : ℝ) ≤ (m : ℝ) + 1 / 2 := by
    simpa [halfIntegerStrip] using hz.2
  have hlowR : (-1 : ℝ) < ((k - m : ℤ) : ℝ) := by
    push_cast
    linarith
  have huppR : ((k - m : ℤ) : ℝ) < (1 : ℝ) := by
    push_cast
    linarith
  have hlowZ : (-1 : ℤ) < k - m := by
    exact_mod_cast hlowR
  have huppZ : k - m < (1 : ℤ) := by
    exact_mod_cast huppR
  have hkm : k = m := by omega
  apply hzm
  simpa [hkm]

private def barnesRaw (κ : ℂ) (A : ℂ → ℂ) (z : ℂ) : ℂ :=
  κ * A z / sinePi z ^ 2

private def barnesExtension
    (κ : ℂ) (m : ℤ) (P : ℂ → ℂ) (z : ℂ) : ℂ :=
  κ * P z / sineSlope m z ^ 2

private theorem barnesRaw_eq_extension_of_mem_strip_of_ne
    {κ : ℂ} {A P : ℂ → ℂ} (m : ℤ) {z : ℂ}
    (hfactor : ∀ w : ℂ,
      A w = (w - (m : ℂ)) ^ 2 * P w)
    (hz : z ∈ halfIntegerStrip m) (hzm : z ≠ (m : ℂ)) :
    barnesRaw κ A z = barnesExtension κ m P z := by
  have hslope := sineSlope_ne_zero_on_strip m hz
  rw [barnesRaw, barnesExtension, hfactor,
    sinePi_eq_sub_mul_sineSlope m z]
  field_simp [sub_ne_zero.mpr hzm, hslope]
  <;> ring

private theorem barnesExtension_differentiableOn
    {κ : ℂ} {P : ℂ → ℂ} (m : ℤ)
    (hP : DifferentiableOn ℂ P (halfIntegerStrip m)) :
    DifferentiableOn ℂ (barnesExtension κ m P) (halfIntegerStrip m) := by
  have hnum : DifferentiableOn ℂ (fun z => κ * P z) (halfIntegerStrip m) :=
    differentiableOn_const.mul hP
  have hden : DifferentiableOn ℂ (fun z => sineSlope m z ^ 2)
      (halfIntegerStrip m) :=
    (sineSlope_differentiable m).differentiableOn.pow 2
  simpa [barnesExtension] using
    hnum.div hden (fun z hz => pow_ne_zero 2 (sineSlope_ne_zero_on_strip m hz))

/-- One-strip Barnes contour shift across an integer whose apparent pole is
removed by the supplied quadratic factorization. -/
theorem barnes_one_strip_shift
    {κ : ℂ} {A P : ℂ → ℂ} (m : ℤ)
    (hfactor : ∀ z : ℂ,
      A z = (z - (m : ℂ)) ^ 2 * P z)
    (hP : DifferentiableOn ℂ P (halfIntegerStrip m))
    (hleft : MeasureTheory.Integrable
      (fun y : ℝ => barnesRaw κ A
        (verticalPoint ((m : ℝ) - 1 / 2) y)))
    (hright : MeasureTheory.Integrable
      (fun y : ℝ => barnesRaw κ A
        (verticalPoint ((m : ℝ) + 1 / 2) y)))
    (htop : Tendsto
      (fun T : ℝ =>
        ∫ x in ((m : ℝ) - 1 / 2)..((m : ℝ) + 1 / 2),
          barnesExtension κ m P
            ((x : ℂ) + (T : ℂ) * Complex.I))
      atTop (𝓝 0))
    (hbottom : Tendsto
      (fun T : ℝ =>
        ∫ x in ((m : ℝ) - 1 / 2)..((m : ℝ) + 1 / 2),
          barnesExtension κ m P
            ((x : ℂ) - (T : ℂ) * Complex.I))
      atTop (𝓝 0)) :
    (∫ y : ℝ, barnesRaw κ A
      (verticalPoint ((m : ℝ) - 1 / 2) y)) =
    ∫ y : ℝ, barnesRaw κ A
      (verticalPoint ((m : ℝ) + 1 / 2) y) := by
  let a : ℝ := (m : ℝ) - 1 / 2
  let b : ℝ := (m : ℝ) + 1 / 2
  have hab : a ≤ b := by dsimp [a, b]; linarith
  have hstrip : closedVerticalStrip a b = halfIntegerStrip m := by
    ext z
    simp [closedVerticalStrip, halfIntegerStrip, a, b]
  have hleftEq :
      (fun y : ℝ => barnesExtension κ m P (verticalPoint a y)) =
        fun y : ℝ => barnesRaw κ A (verticalPoint a y) := by
    funext y
    symm
    apply barnesRaw_eq_extension_of_mem_strip_of_ne m hfactor
    · simp [halfIntegerStrip, verticalPoint, a]
    · intro h
      have hre := congrArg Complex.re h
      simp [verticalPoint, a] at hre
  have hrightEq :
      (fun y : ℝ => barnesExtension κ m P (verticalPoint b y)) =
        fun y : ℝ => barnesRaw κ A (verticalPoint b y) := by
    funext y
    symm
    apply barnesRaw_eq_extension_of_mem_strip_of_ne m hfactor
    · simp [halfIntegerStrip, verticalPoint, b]
    · intro h
      have hre := congrArg Complex.re h
      simp [verticalPoint, b] at hre
  have hleftExt : MeasureTheory.Integrable
      (fun y : ℝ => barnesExtension κ m P (verticalPoint a y)) := by
    rw [hleftEq]
    simpa [a] using hleft
  have hrightExt : MeasureTheory.Integrable
      (fun y : ℝ => barnesExtension κ m P (verticalPoint b y)) := by
    rw [hrightEq]
    simpa [b] using hright
  have hshift := verticalIntegral_eq_of_horizontal_tendsto
    (F := barnesExtension κ m P) hab
    (by rw [hstrip]; exact barnesExtension_differentiableOn m hP)
    hleftExt hrightExt
    (by simpa [a, b] using htop)
    (by simpa [a, b] using hbottom)
  rw [hleftEq, hrightEq] at hshift
  simpa [a, b] using hshift

#print axioms verticalIntegral_eq_of_horizontal_tendsto
#print axioms sineSlope_ne_zero_on_strip
#print axioms barnesRaw_eq_extension_of_mem_strip_of_ne
#print axioms barnes_one_strip_shift

end RamanujanChallenge.P27.Q6382
