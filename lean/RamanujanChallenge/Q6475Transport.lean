import RamanujanChallenge.Problem27BarnesTelescoper

open Filter Set MeasureTheory Topology
open scoped BigOperators Interval Real

noncomputable section

namespace RamanujanChallenge.P27.Q6475

private def vp (x y : ℝ) : ℂ :=
  (x : ℂ) + (y : ℂ) * Complex.I

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
orientations and integral APIs after the extension strip equality has been
provided by the rectangle theorem. -/
theorem raw_vertical_integral_eq
    {raw ext : ℂ → ℂ} {a b : ℝ}
    (hExtShift :
      (∫ y : ℝ, ext (vp a y)) = ∫ y : ℝ, ext (vp b y))
    (hLeftEq : ∀ y : ℝ, raw (vp a y) = ext (vp a y))
    (hRightEq : ∀ y : ℝ, raw (vp b y) = ext (vp b y)) :
    (∫ y : ℝ, raw (vp a y)) = ∫ y : ℝ, raw (vp b y) := by
  calc
    (∫ y : ℝ, raw (vp a y)) = ∫ y : ℝ, ext (vp a y) := by
      apply integral_congr_ae
      filter_upwards with y
      exact hLeftEq y
    _ = ∫ y : ℝ, ext (vp b y) := hExtShift
    _ = ∫ y : ℝ, raw (vp b y) := by
      apply integral_congr_ae
      filter_upwards with y
      exact (hRightEq y).symm

#print axioms horizontal_tendsto_of_eventually_eq
#print axioms horizontal_bottom_tendsto_of_eventually_eq
#print axioms raw_vertical_integral_eq

end RamanujanChallenge.P27.Q6475
