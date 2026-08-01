import RamanujanChallenge.Problem27BarnesTelescoper
import Mathlib.Analysis.Complex.CauchyIntegral
import Mathlib.Analysis.Complex.RemovableSingularity
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Complex
import Mathlib.MeasureTheory.Integral.IntegralEqImproper
import Mathlib.Analysis.SpecialFunctions.PolynomialExp

/-!
# Problem 2.7: contour shifts for the Barnes integrals

This file constructs the removable extension used when a half-integer Barnes
line is shifted across an integer.  The raw totalized quotient is never treated
as continuous at its removable singularity.
-/

open Filter Set MeasureTheory Topology
open scoped BigOperators Interval Real

noncomputable section

namespace RamanujanChallenge.P27

def closedVerticalStrip27 (a b : ℝ) : Set ℂ :=
  {z | a ≤ z.re ∧ z.re ≤ b}

def verticalPoint27 (x y : ℝ) : ℂ :=
  (x : ℂ) + (y : ℂ) * Complex.I

/-- Rectangle Cauchy plus vanishing horizontal sides.  The vertical integrals
are parameterized by their imaginary coordinate. -/
theorem verticalIntegral_eq_of_horizontal_tendsto27
    {F : ℂ → ℂ} {a b : ℝ} (hab : a ≤ b)
    (hF : DifferentiableOn ℂ F (closedVerticalStrip27 a b))
    (hleft : Integrable (fun y : ℝ => F (verticalPoint27 a y)))
    (hright : Integrable (fun y : ℝ => F (verticalPoint27 b y)))
    (htop : Tendsto
      (fun T : ℝ =>
        ∫ x in a..b, F ((x : ℂ) + (T : ℂ) * Complex.I))
      atTop (𝓝 0))
    (hbottom : Tendsto
      (fun T : ℝ =>
        ∫ x in a..b, F ((x : ℂ) - (T : ℂ) * Complex.I))
      atTop (𝓝 0)) :
    (∫ y : ℝ, F (verticalPoint27 a y)) =
      ∫ y : ℝ, F (verticalPoint27 b y) := by
  have hfinite : ∀ T : ℝ, 0 ≤ T →
      (∫ y in -T..T, F (verticalPoint27 a y)) -
          (∫ y in -T..T, F (verticalPoint27 b y)) =
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
          simpa [z, w, closedVerticalStrip27, uIcc_of_le hab] using hu.1))
    dsimp [z, w] at hrect
    simp only [Complex.sub_re, Complex.sub_im, Complex.add_re,
      Complex.add_im, Complex.mul_re, Complex.mul_im,
      Complex.ofReal_re, Complex.ofReal_im,
      Complex.I_re, Complex.I_im, verticalPoint27, smul_eq_mul] at hrect ⊢
    norm_num at hrect ⊢
    have hcore :
        Complex.I *
            ((∫ y in -T..T, F (verticalPoint27 b y)) -
              ∫ y in -T..T, F (verticalPoint27 a y)) =
          (∫ x in a..b, F ((x : ℂ) + (T : ℂ) * Complex.I)) -
            ∫ x in a..b, F ((x : ℂ) - (T : ℂ) * Complex.I) := by
      simp only [verticalPoint27]
      linear_combination hrect
    have hi := congrArg (fun q : ℂ => Complex.I * q) hcore
    simp only [mul_sub, ← mul_assoc, Complex.I_mul_I, neg_one_mul,
      verticalPoint27] at hi
    linear_combination hi

  have hleft_lim : Tendsto
      (fun T : ℝ => ∫ y in -T..T, F (verticalPoint27 a y))
      atTop (𝓝 (∫ y : ℝ, F (verticalPoint27 a y))) :=
    intervalIntegral_tendsto_integral
      hleft tendsto_neg_atTop_atBot tendsto_id
  have hright_lim : Tendsto
      (fun T : ℝ => ∫ y in -T..T, F (verticalPoint27 b y))
      atTop (𝓝 (∫ y : ℝ, F (verticalPoint27 b y))) :=
    intervalIntegral_tendsto_integral
      hright tendsto_neg_atTop_atBot tendsto_id
  have hfinite_eventually : ∀ᶠ T : ℝ in atTop,
      (∫ y in -T..T, F (verticalPoint27 a y)) -
          (∫ y in -T..T, F (verticalPoint27 b y)) =
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
  have heq :
      (fun T : ℝ =>
        (∫ y in -T..T, F (verticalPoint27 a y)) -
          ∫ y in -T..T, F (verticalPoint27 b y)) =ᶠ[atTop]
      (fun T : ℝ => Complex.I *
        ((∫ x in a..b,
            F ((x : ℂ) + (T : ℂ) * Complex.I)) -
         (∫ x in a..b,
            F ((x : ℂ) - (T : ℂ) * Complex.I)))) :=
    hfinite_eventually
  have hlhs_zero := Filter.Tendsto.congr' heq.symm hrhs
  have hsub :
      (∫ y : ℝ, F (verticalPoint27 a y)) -
          (∫ y : ℝ, F (verticalPoint27 b y)) = 0 :=
    tendsto_nhds_unique hlhs hlhs_zero
  exact sub_eq_zero.mp hsub

/-- A uniform polynomial-times-exponential bound kills a horizontal edge. -/
theorem horizontalIntegral_tendsto_zero_of_pow_exp_bound27
    {F : ℂ → ℂ} {a b C c : ℝ} {d : ℕ}
    (hC : 0 ≤ C) (hc : 0 < c)
    (hbound : ∀ᶠ T : ℝ in atTop,
      ∀ x ∈ [[a, b]],
        ‖F ((x : ℂ) + (T : ℂ) * Complex.I)‖ ≤
          C * (c * T) ^ d * Real.exp (-(c * T))) :
    Tendsto
      (fun T : ℝ =>
        ∫ x in a..b, F ((x : ℂ) + (T : ℂ) * Complex.I))
      atTop (𝓝 0) := by
  have hnorm : ∀ᶠ T : ℝ in atTop,
      ‖∫ x in a..b,
          F ((x : ℂ) + (T : ℂ) * Complex.I)‖ ≤
        (C * (c * T) ^ d * Real.exp (-(c * T))) * |b - a| := by
    filter_upwards [hbound] with T hT
    exact intervalIntegral.norm_integral_le_of_norm_le_const
      (fun x hx => hT x (uIoc_subset_uIcc hx))
  have hscale : Tendsto (fun T : ℝ => c * T) atTop atTop :=
    tendsto_id.const_mul_atTop' hc
  have hpolyexp : Tendsto
      (fun T : ℝ => (c * T) ^ d * Real.exp (-(c * T)))
      atTop (𝓝 0) :=
    (Real.tendsto_pow_mul_exp_neg_atTop_nhds_zero d).comp hscale
  have hmajor : Tendsto
      (fun T : ℝ =>
        (C * (c * T) ^ d * Real.exp (-(c * T))) * |b - a|)
      atTop (𝓝 0) := by
    simpa [mul_assoc] using
      (tendsto_const_nhds.mul hpolyexp).mul_const |b - a|
  apply tendsto_zero_iff_norm_tendsto_zero.mpr
  exact squeeze_zero'
      (Eventually.of_forall fun _ => norm_nonneg _)
      hnorm hmajor

def sinePi27 (z : ℂ) : ℂ :=
  Complex.sin ((Real.pi : ℂ) * z)

/-- Holomorphic extension of `sin (π z) / (z-m)` at the integer `m`. -/
def sineSlope27 (m : ℤ) : ℂ → ℂ :=
  dslope sinePi27 (m : ℂ)

@[simp] theorem sinePi_int27 (m : ℤ) :
    sinePi27 (m : ℂ) = 0 := by
  rw [sinePi27, Complex.sin_eq_zero_iff]
  exact ⟨m, by push_cast; ring⟩

theorem sinePi_eq_sub_mul_sineSlope27 (m : ℤ) (z : ℂ) :
    sinePi27 z = (z - (m : ℂ)) * sineSlope27 m z := by
  have h := sub_smul_dslope sinePi27 (m : ℂ) z
  simpa [sineSlope27, smul_eq_mul, sinePi_int27] using h.symm

theorem sineSlope_differentiable27 (m : ℤ) :
    Differentiable ℂ (sineSlope27 m) := by
  have hsine : Differentiable ℂ sinePi27 := by
    intro z
    unfold sinePi27
    fun_prop
  have h :=
    (Complex.differentiableOn_dslope
      (f := sinePi27) (s := Set.univ) (c := (m : ℂ)) univ_mem).2
      hsine.differentiableOn
  simpa [sineSlope27] using differentiableOn_univ.mp h

def stripLeft27 (m : ℤ) : ℝ :=
  (m : ℝ) - 1 / 2

def stripRight27 (m : ℤ) : ℝ :=
  (m : ℝ) + 1 / 2

private theorem int_eq_of_mem_half_strip27 {k m : ℤ}
    (h₁ : (m : ℝ) - 1 / 2 ≤ (k : ℝ))
    (h₂ : (k : ℝ) ≤ (m : ℝ) + 1 / 2) : k = m := by
  have h₁r : (((2 * m - 1 : ℤ) : ℝ)) ≤ ((2 * k : ℤ) : ℝ) := by
    push_cast
    linarith
  have h₂r : ((2 * k : ℤ) : ℝ) ≤ (((2 * m + 1 : ℤ) : ℝ)) := by
    push_cast
    linarith
  have h₁z : (2 * m - 1 : ℤ) ≤ 2 * k := by exact_mod_cast h₁r
  have h₂z : 2 * k ≤ (2 * m + 1 : ℤ) := by exact_mod_cast h₂r
  omega

theorem sineSlope_ne_zero_on_strip27
    (m : ℤ) {z : ℂ}
    (hz : stripLeft27 m ≤ z.re ∧ z.re ≤ stripRight27 m) :
    sineSlope27 m z ≠ 0 := by
  by_cases hzm : z = (m : ℂ)
  · subst z
    rw [sineSlope27, dslope_same]
    have hderiv :
        deriv sinePi27 (m : ℂ) =
          Complex.cos ((Real.pi : ℂ) * (m : ℂ)) * (Real.pi : ℂ) := by
      unfold sinePi27
      simpa using
        (((hasDerivAt_id (m : ℂ)).const_mul (Real.pi : ℂ)).csin.deriv)
    rw [hderiv]
    apply mul_ne_zero
    · have hr : Real.cos ((m : ℝ) * Real.pi) ≠ 0 := by
        rw [Real.cos_int_mul_pi]
        exact zpow_ne_zero _ (by norm_num)
      rw [show (Real.pi : ℂ) * (m : ℂ) =
          (((m : ℝ) * Real.pi : ℝ) : ℂ) by push_cast; ring,
        ← Complex.ofReal_cos]
      exact Complex.ofReal_ne_zero.mpr hr
    · exact Complex.ofReal_ne_zero.mpr Real.pi_ne_zero
  · intro hslope
    have hsin : sinePi27 z = 0 := by
      rw [sinePi_eq_sub_mul_sineSlope27 m z, hslope, mul_zero]
    rcases Complex.sin_eq_zero_iff.mp hsin with ⟨k, hk⟩
    have hzk : z = (k : ℂ) := by
      have hpi : (Real.pi : ℂ) ≠ 0 :=
        Complex.ofReal_ne_zero.mpr Real.pi_ne_zero
      apply mul_left_cancel₀ hpi
      simpa [mul_comm, mul_left_comm, mul_assoc] using hk
    subst z
    have hkm : k = m := by
      apply int_eq_of_mem_half_strip27
      · simpa [stripLeft27] using hz.1
      · simpa [stripRight27] using hz.2
    exact hzm (by simpa [hkm])

end RamanujanChallenge.P27
