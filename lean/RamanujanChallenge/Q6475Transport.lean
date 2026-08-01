import RamanujanChallenge.Problem27BarnesShift

open Filter Set MeasureTheory Topology
open scoped BigOperators Interval Real

noncomputable section

namespace RamanujanChallenge.P27.Q6475

/-- Exact Mathlib API pattern needed in the contour file: transport a raw
horizontal limit to a removable extension by eventual pointwise equality on
the interval of integration. -/
theorem horizontal_tendsto_of_eventually_eq
    {raw ext : ℂ → ℂ} {a b : ℝ}
    (hraw : Tendsto
      (fun T : ℝ => ∫ x in a..b,
        raw ((x : ℂ) + (T : ℂ) * Complex.I))
      atTop (𝓝 0))
    (heq : ∀ᶠ T : ℝ in atTop, ∀ x ∈ [[a, b]],
      raw ((x : ℂ) + (T : ℂ) * Complex.I) =
        ext ((x : ℂ) + (T : ℂ) * Complex.I)) :
    Tendsto
      (fun T : ℝ => ∫ x in a..b,
        ext ((x : ℂ) + (T : ℂ) * Complex.I))
      atTop (𝓝 0) := by
  refine hraw.congr' ?_
  filter_upwards [heq] with T hT
  exact intervalIntegral.integral_congr hT

/-- Bottom-edge version. -/
theorem horizontal_bottom_tendsto_of_eventually_eq
    {raw ext : ℂ → ℂ} {a b : ℝ}
    (hraw : Tendsto
      (fun T : ℝ => ∫ x in a..b,
        raw ((x : ℂ) - (T : ℂ) * Complex.I))
      atTop (𝓝 0))
    (heq : ∀ᶠ T : ℝ in atTop, ∀ x ∈ [[a, b]],
      raw ((x : ℂ) - (T : ℂ) * Complex.I) =
        ext ((x : ℂ) - (T : ℂ) * Complex.I)) :
    Tendsto
      (fun T : ℝ => ∫ x in a..b,
        ext ((x : ℂ) - (T : ℂ) * Complex.I))
      atTop (𝓝 0) := by
  refine hraw.congr' ?_
  filter_upwards [heq] with T hT
  exact intervalIntegral.integral_congr hT

/-- Full generic one-strip packaging, used only to compile-check the exact
orientations and integral APIs. -/
theorem raw_vertical_integral_eq
    {raw ext : ℂ → ℂ} {a b : ℝ} (hab : a ≤ b)
    (hF : DifferentiableOn ℂ ext (closedVerticalStrip27 a b))
    (hleft : Integrable (fun y : ℝ => ext (verticalPoint27 a y)))
    (hright : Integrable (fun y : ℝ => ext (verticalPoint27 b y)))
    (hrawTop : Tendsto
      (fun T : ℝ => ∫ x in a..b,
        raw ((x : ℂ) + (T : ℂ) * Complex.I))
      atTop (𝓝 0))
    (hrawBottom : Tendsto
      (fun T : ℝ => ∫ x in a..b,
        raw ((x : ℂ) - (T : ℂ) * Complex.I))
      atTop (𝓝 0))
    (hTopEq : ∀ᶠ T : ℝ in atTop, ∀ x ∈ [[a, b]],
      raw ((x : ℂ) + (T : ℂ) * Complex.I) =
        ext ((x : ℂ) + (T : ℂ) * Complex.I))
    (hBottomEq : ∀ᶠ T : ℝ in atTop, ∀ x ∈ [[a, b]],
      raw ((x : ℂ) - (T : ℂ) * Complex.I) =
        ext ((x : ℂ) - (T : ℂ) * Complex.I))
    (hLeftEq : ∀ y : ℝ,
      raw (verticalPoint27 a y) = ext (verticalPoint27 a y))
    (hRightEq : ∀ y : ℝ,
      raw (verticalPoint27 b y) = ext (verticalPoint27 b y)) :
    (∫ y : ℝ, raw (verticalPoint27 a y)) =
      ∫ y : ℝ, raw (verticalPoint27 b y) := by
  have htop := horizontal_tendsto_of_eventually_eq hrawTop hTopEq
  have hbottom :=
    horizontal_bottom_tendsto_of_eventually_eq hrawBottom hBottomEq
  have hshift := verticalIntegral_eq_of_horizontal_tendsto27
    hab hF hleft hright htop hbottom
  calc
    (∫ y : ℝ, raw (verticalPoint27 a y)) =
        ∫ y : ℝ, ext (verticalPoint27 a y) := by
      apply integral_congr_ae
      filter_upwards with y
      exact hLeftEq y
    _ = ∫ y : ℝ, ext (verticalPoint27 b y) := hshift
    _ = ∫ y : ℝ, raw (verticalPoint27 b y) := by
      apply integral_congr_ae
      filter_upwards with y
      exact (hRightEq y).symm

#print axioms horizontal_tendsto_of_eventually_eq
#print axioms horizontal_bottom_tendsto_of_eventually_eq
#print axioms raw_vertical_integral_eq

end RamanujanChallenge.P27.Q6475
