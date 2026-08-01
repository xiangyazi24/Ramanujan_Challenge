import RamanujanChallenge.Problem27BarnesTelescoper

open Filter Set MeasureTheory Topology
open scoped BigOperators Interval Real

noncomputable section

namespace RamanujanChallenge.P27.Q6475

private def vp (x y : ℝ) : ℂ :=
  (x : ℂ) + (y : ℂ) * Complex.I

private def left (m : ℕ) : ℝ := (m : ℝ) - 1 / 2
private def right (m : ℕ) : ℝ := (m : ℝ) + 1 / 2
private def strip (m : ℕ) : Set ℂ :=
  {z | left m ≤ z.re ∧ z.re ≤ right m}

private theorem left_le_right (m : ℕ) : left m ≤ right m := by
  unfold left right
  norm_num

private theorem top_mem_strip {m : ℕ} {x T : ℝ}
    (hx : x ∈ [[left m, right m]]) :
    (x : ℂ) + (T : ℂ) * Complex.I ∈ strip m := by
  have hx' : x ∈ Set.Icc (left m) (right m) := by
    simpa [uIcc_of_le (left_le_right m)] using hx
  simpa [strip] using hx'

private theorem bottom_mem_strip {m : ℕ} {x T : ℝ}
    (hx : x ∈ [[left m, right m]]) :
    (x : ℂ) - (T : ℂ) * Complex.I ∈ strip m := by
  have hx' : x ∈ Set.Icc (left m) (right m) := by
    simpa [uIcc_of_le (left_le_right m)] using hx
  simpa [strip] using hx'

private theorem top_ne_center {m : ℕ} {x T : ℝ} (hT : 1 ≤ T) :
    (x : ℂ) + (T : ℂ) * Complex.I ≠ (m : ℂ) := by
  intro h
  have him : T = 0 := by
    simpa using congrArg Complex.im h
  linarith

private theorem bottom_ne_center {m : ℕ} {x T : ℝ} (hT : 1 ≤ T) :
    (x : ℂ) - (T : ℂ) * Complex.I ≠ (m : ℂ) := by
  intro h
  have him : -T = 0 := by
    simpa using congrArg Complex.im h
  linarith

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

/-- Complete top-edge transport from a strip raw=extension theorem. -/
theorem top_transport_from_strip_eq
    {raw ext : ℂ → ℂ} {m : ℕ}
    (hraw : Tendsto
      (fun T : ℝ => ∫ x in left m..right m,
        raw ((x : ℂ) + (T : ℂ) * Complex.I))
      atTop (𝓝 0))
    (heq : ∀ {z : ℂ}, z ∈ strip m → z ≠ (m : ℂ) → raw z = ext z) :
    Tendsto
      (fun T : ℝ => ∫ x in left m..right m,
        ext ((x : ℂ) + (T : ℂ) * Complex.I))
      atTop (𝓝 0) := by
  refine hraw.congr' ?_
  filter_upwards [eventually_ge_atTop (1 : ℝ)] with T hT
  exact intervalIntegral.integral_congr fun x hx =>
    heq (top_mem_strip hx) (top_ne_center hT)

/-- Complete bottom-edge transport from a strip raw=extension theorem. -/
theorem bottom_transport_from_strip_eq
    {raw ext : ℂ → ℂ} {m : ℕ}
    (hraw : Tendsto
      (fun T : ℝ => ∫ x in left m..right m,
        raw ((x : ℂ) - (T : ℂ) * Complex.I))
      atTop (𝓝 0))
    (heq : ∀ {z : ℂ}, z ∈ strip m → z ≠ (m : ℂ) → raw z = ext z) :
    Tendsto
      (fun T : ℝ => ∫ x in left m..right m,
        ext ((x : ℂ) - (T : ℂ) * Complex.I))
      atTop (𝓝 0) := by
  refine hraw.congr' ?_
  filter_upwards [eventually_ge_atTop (1 : ℝ)] with T hT
  exact intervalIntegral.integral_congr fun x hx =>
    heq (bottom_mem_strip hx) (bottom_ne_center hT)

/-- Full generic one-strip packaging, used only to compile-check the exact
orientations after the extension strip equality has been provided by the
rectangle theorem. -/
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
#print axioms top_transport_from_strip_eq
#print axioms bottom_transport_from_strip_eq
#print axioms raw_vertical_integral_eq

end RamanujanChallenge.P27.Q6475
