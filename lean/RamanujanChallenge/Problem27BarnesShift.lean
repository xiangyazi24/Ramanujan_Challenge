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
    linear_combination Complex.I * hrect

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
  have hlhs_zero := hrhs.congr' hfinite_eventually.symm
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
    exact intervalIntegral.norm_integral_le_of_norm_le_const hT
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
  rw [← tendsto_zero_iff_norm_tendsto_zero]
  exact squeeze_zero'
    (Eventually.of_forall fun _ => norm_nonneg _)
    hnorm hmajor

end RamanujanChallenge.P27
