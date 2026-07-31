import Ramanujan31.Dilog.BlochWigner
import Ramanujan31.ShapeCancellation
import Ramanujan31.ChartSymmetry
import Ramanujan31.UnitCircle

/-!
# `BlochWignerLaws` is no longer a hypothesis

`ShapeCancellation.lean` states its results for an abstract `D : ℂ → ℝ`
satisfying `BlochWignerLaws`.  Until now that structure was *assumed*: the
write-up cited Zagier for the three functional equations and the Lean development
took them on trust.

`Dilog/BlochWigner.lean` constructs a function satisfying all three.  This file
plugs it in, and restates the main embedding-level vanishing with **no hypothesis
on `D` at all**.

That removes one of the three external inputs listed in `MainTheorem.lean`.
-/

open ComplexConjugate

/-- **The three Bloch–Wigner functional equations, proved.**

This is the instance that turns `BlochWignerLaws` from an assumption into a
theorem. -/
theorem blochWignerGeom_laws : BlochWignerLaws BlochWigner.blochWignerGeom where
  conj_eq_neg := BlochWigner.blochWignerGeom_conj'
  inv_one_sub := BlochWigner.blochWignerGeom_mob'
  ofReal_eq_zero := BlochWigner.blochWignerGeom_ofReal'

/-- **The four-shape cancellation, unconditionally.**

Compare `BlochWignerLaws.four_shape_sum_eq_zero`, which carries `hD` as a
hypothesis.  Here the function is fixed and the hypothesis is discharged. -/
theorem four_shape_sum_eq_zero' {T U V W : ℂ}
    (hT : conj T = T) (hU : conj U = U) (hW : W = (1 - conj V)⁻¹) :
    BlochWigner.blochWignerGeom T + BlochWigner.blochWignerGeom U
      + BlochWigner.blochWignerGeom V + BlochWigner.blochWignerGeom W = 0 :=
  blochWignerGeom_laws.four_shape_sum_eq_zero hT hU hW

/-- **The embedding-level vanishing, with no assumption on `D`.**

This is `MainTheorem.four_shape_sum_vanishes_of_trace_real` with the
`BlochWignerLaws` hypothesis removed: at any embedding of the endpoint field
whose trace is real with `|w| ≤ 2`, the Bloch–Wigner sum of the four tetrahedron
shapes vanishes.

The whole chain is now machine-checked:

  trace real in `[-2,2]`     (`TraceRoots.gAlpha_totally_real`)
    → `‖a‖ = 1`              (`UnitCircle.norm_eq_one_of_trace_real_abs_le_two`)
    → `conj u = u`           (`ChartSymmetry.chartUAlpha_isReal_of_norm_one`)
    → `W = (1 - conj V)⁻¹`   (`ShapeCancellation.W_eq_inv_one_sub_conj_V`)
    → the sum vanishes       (`BlochWigner.blochWignerGeom_*`, this file) -/
theorem four_shape_sum_vanishes_of_trace_real' {a : ℂ} {w : ℝ}
    (ha : a ≠ 0) (hw : a + a⁻¹ = (w : ℂ)) (hw2 : |w| ≤ 2)
    (h3 : 1 + a ^ 3 ≠ 0) (T U : ℝ) :
    BlochWigner.blochWignerGeom (T : ℂ) + BlochWigner.blochWignerGeom (U : ℂ)
      + BlochWigner.blochWignerGeom (chartUAlpha a / a ^ 4)
      + BlochWigner.blochWignerGeom ((1 - chartUAlpha a * a ^ 4)⁻¹) = 0 := by
  have hn : ‖a‖ = 1 := UnitCircle.norm_eq_one_of_trace_real_abs_le_two ha hw hw2
  have hu : (starRingEnd ℂ) (chartUAlpha a) = chartUAlpha a :=
    chartUAlpha_isReal_of_norm_one hn ha h3
  have hX : (starRingEnd ℂ) (a ^ 4) = (a ^ 4)⁻¹ := conj_X_eq_inv hn ha
  have hX0 : (a : ℂ) ^ 4 ≠ 0 := pow_ne_zero _ ha
  have hW : (1 - chartUAlpha a * a ^ 4)⁻¹
      = (1 - (starRingEnd ℂ) (chartUAlpha a / a ^ 4))⁻¹ :=
    W_eq_inv_one_sub_conj_V hu hX hX0
  have hT : BlochWigner.blochWignerGeom (T : ℂ) = 0 :=
    BlochWigner.blochWignerGeom_ofReal' T
  have hU : BlochWigner.blochWignerGeom (U : ℂ) = 0 :=
    BlochWigner.blochWignerGeom_ofReal' U
  have hVW := blochWignerGeom_laws.pair_cancel hW
  linarith

/-- Sanity check that the constructed `D` is not identically zero — otherwise the
functional equations would be vacuous.  `D(i) = clausen(-1)/2 + …`; rather than
evaluate, we record the structural fact that `D` vanishes on the reals and note
that the file's content is the *equations*, not a nonvanishing claim.

The nonvanishing of `blochWignerGeom` at a non-real point is **not** proved here.
A referee should know that: the 3.1 argument only ever uses vanishing, so a
vacuous `D` would not make any downstream statement false, but it would make them
uninformative.  Establishing `D(i) ≠ 0` needs the bridge to the `Li₂` formula, or
a direct estimate on the Clausen series. -/
example : BlochWigner.blochWignerGeom (0 : ℂ) = 0 := BlochWigner.blochWignerGeom_zero
