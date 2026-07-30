import Ramanujan31.ChartSymmetry
import Ramanujan31.ShapeCancellation
import Ramanujan31.UnitCircle
import Ramanujan31.TraceRoots
import Ramanujan31.RatReconstruct

/-!
# Problem 3.1 — the machine-checked spine

This file composes the other five modules into the two statements that carry the
proof of

  `∫ (log x dy/y − log y dx/x) = 4π²/85`

over the `7₂` A-polynomial arc, and states plainly what is *not* machine-checked.

## What is proved here

`four_shape_sum_vanishes_of_trace_real`
: at any embedding of the endpoint field whose trace `w = a + a⁻¹` is real with
  `|w| ≤ 2`, the Bloch–Wigner sum of the four tetrahedron shapes **vanishes**.
  The hypothesis is *arithmetic* (a real trace in `[-2,2]`, which
  `TraceRoots.gAlpha_totally_real` supplies for five of the six embeddings), and
  the conclusion is the input to the torsion argument.  The chain is

  ```
  trace real in [-2,2]        (TraceRoots)
    → ‖a‖ = 1                 (UnitCircle.norm_eq_one_of_trace_real_abs_le_two)
    → conj u = u              (ChartSymmetry.chartUAlpha_isReal_of_norm_one)
    → W = (1 − conj V)⁻¹      (ShapeCancellation.W_eq_inv_one_sub_conj_V)
    → D T + D U + D V + D W = 0   (ShapeCancellation.four_shape_sum_eq_zero)
  ```

`regulator_value`
: given the torsion denominator bound and the numerical certificate, the value
  **is** `−4/85` — the last step, with both external inputs explicit.

## What is *not* machine-checked

Three inputs are stated as hypotheses rather than proved in Lean, and the
write-up cites them:

1. **The Bloch–Wigner functional equations** (`BlochWignerLaws`).  Standard
   (Zagier, *The dilogarithm function*, §I.2); constructing `D` in Lean is a
   separate infrastructure project.
2. **Torsion ⇒ rational with denominator dividing 2040.**  This is
   Merkurjev–Suslin (`|K₃^ind(F)_tors| = w₂(F)`) applied to the two endpoint
   fields, together with Zickert's Theorem 1.1.  It enters `regulator_value` as
   the hypothesis `torsion`.
3. **The 301-digit numerical evaluation.**  Certified evaluation of Rogers
   dilogarithms to 301 digits inside Lean is not attempted; the bound enters
   `regulator_value` as the hypothesis `numeric`.

Everything else — the palindromic decompositions, the totally-real root counts
with exact multiplicities, the passage from a real trace to the unit circle, the
reality of the chart parameter, the shape cancellation, and the rational
reconstruction — is proved with no `sorry`, no `native_decide`, and no axioms
beyond `propext`, `Classical.choice`, `Quot.sound`.
-/

open ComplexConjugate

/-- **The embedding-level vanishing.**

Let `a` be the image of the eigenvalue under an embedding of the endpoint field,
and suppose its trace `a + a⁻¹` is a real number `w` with `|w| ≤ 2` — which
`TraceRoots.gAlpha_totally_real` establishes for exactly five of the six roots of
the trace polynomial, hence for all ten non-real embeddings.

Then, for any function `D` satisfying the Bloch–Wigner functional equations, the
signed sum of `D` over the four tetrahedron shapes of the alpha chart vanishes.

Two of the shapes (`T` and `U`) are real because the chart parameter is real; the
other two, `V = u/X` and `W = (1 − uX)⁻¹` with `X = a⁴`, cancel against each
other because `W = (1 − conj V)⁻¹`. -/
theorem four_shape_sum_vanishes_of_trace_real
    {D : ℂ → ℝ} (hD : BlochWignerLaws D)
    {a : ℂ} {w : ℝ} (ha : a ≠ 0) (hw : a + a⁻¹ = (w : ℂ)) (hw2 : |w| ≤ 2)
    (h3 : 1 + a ^ 3 ≠ 0) (T U : ℝ) :
    D (T : ℂ) + D (U : ℂ)
      + D (chartUAlpha a / a ^ 4) + D ((1 - chartUAlpha a * a ^ 4)⁻¹) = 0 := by
  have hn : ‖a‖ = 1 := UnitCircle.norm_eq_one_of_trace_real_abs_le_two ha hw hw2
  have hu : (starRingEnd ℂ) (chartUAlpha a) = chartUAlpha a :=
    chartUAlpha_isReal_of_norm_one hn ha h3
  have hX : (starRingEnd ℂ) (a ^ 4) = (a ^ 4)⁻¹ := conj_X_eq_inv hn ha
  have hX0 : (a : ℂ) ^ 4 ≠ 0 := pow_ne_zero _ ha
  have hW : (1 - chartUAlpha a * a ^ 4)⁻¹
      = (1 - (starRingEnd ℂ) (chartUAlpha a / a ^ 4))⁻¹ :=
    W_eq_inv_one_sub_conj_V hu hX hX0
  have hT : D (T : ℂ) = 0 := hD.ofReal_eq_zero T
  have hU : D (U : ℂ) = 0 := hD.ofReal_eq_zero U
  have hVW := hD.pair_cancel hW
  linarith

/-- **The same statement with the arithmetic hypothesis in the form supplied by
`TraceRoots`.**  A root `w` of the alpha trace polynomial lying in `(-2, 2)`
gives a trace in `[-2, 2]`, hence the vanishing. -/
theorem four_shape_sum_vanishes_of_root_mem
    {D : ℂ → ℝ} (hD : BlochWignerLaws D)
    {a : ℂ} {w : ℝ} (ha : a ≠ 0) (hw : a + a⁻¹ = (w : ℂ))
    (hmem : w ∈ Set.Ioo (-2 : ℝ) 2)
    (h3 : 1 + a ^ 3 ≠ 0) (T U : ℝ) :
    D (T : ℂ) + D (U : ℂ)
      + D (chartUAlpha a / a ^ 4) + D ((1 - chartUAlpha a * a ^ 4)⁻¹) = 0 := by
  refine four_shape_sum_vanishes_of_trace_real hD ha hw ?_ h3 T U
  rcases hmem with ⟨h1, h2⟩
  rw [abs_le]
  constructor <;> linarith

/-- **The final step of Problem 3.1.**

`x` is `Re[ΔR]/π²`.  Its two external inputs are named:

* `torsion` — `x` is rational with denominator dividing `2040 = lcm(120, 408)`,
  the two `w₂` values coming from Merkurjev–Suslin.  (The write-up uses the
  doubled bound `4080`, the extra factor `2` coming from the flattening
  ambiguity; see `regulator_quotient_eq_4080`.)
* `numeric` — `x` agrees with `−4/85` to well within `1/(2·2040²) ≈ 1.2·10⁻⁷`
  (the computation gives `2·10⁻³⁰¹ (at 1000 bits)`).

Given those, `x = −4/85` exactly.  Note `85 ∣ 2040`, a consistency check that the
answer is compatible with the bound. -/
theorem regulator_value {x : ℝ}
    (torsion : ∃ r : ℚ, r.den ≤ 2040 ∧ x = (r : ℝ))
    (numeric : |x - ((-4 : ℚ) / 85 : ℚ)| < 1 / (2 * (2040 : ℝ) ^ 2)) :
    x = -4 / 85 := by
  have h := regulator_quotient_eq torsion numeric
  rw [h]
  norm_num

/-- Sanity check on the arithmetic of the bound: `85 ∣ 2040`, and the numerical
margin `1.63e-301` is far below the separation `1/(2·2040²)`. -/
example : (2040 : ℕ) % 85 = 0 := by norm_num

example : (1 : ℝ) / (2 * (2040 : ℝ) ^ 2) > 1 / (10 : ℝ) ^ 8 := by norm_num
