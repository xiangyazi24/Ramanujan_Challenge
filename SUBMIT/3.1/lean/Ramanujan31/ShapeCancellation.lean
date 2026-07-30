import Mathlib.Tactic

/-!
# The shape cancellation at a unit-circle embedding

This file formalizes Lemmas 4.3 and 4.4 of the Ramanujan Challenge 3.1 write-up:
at an embedding of the endpoint field where the eigenvalue has modulus one,

* `T` and `U` are real, so their Bloch–Wigner values vanish;
* `V` and `W` satisfy `W = (1 - conj V)⁻¹`, so their Bloch–Wigner values cancel.

Together these give `D T + D U + D V + D W = 0`, which is what makes the
extended Bloch class torsion.

We do **not** construct the Bloch–Wigner function here.  Instead we axiomatize
exactly the three properties used, as a structure `BlochWignerLaws`:

  `D (conj z) = - D z`,  `D ((1 - z)⁻¹) = D z`,  `D (x : ℂ) = 0` for real `x`.

These are standard (see Zagier, *The dilogarithm function*, §I.2), and isolating
them keeps this file independent of any particular implementation of `D`.
-/

open ComplexConjugate

/-- The three functional equations of the Bloch–Wigner dilogarithm that the
cancellation argument uses. -/
structure BlochWignerLaws (D : ℂ → ℝ) : Prop where
  /-- `D` is odd under complex conjugation. -/
  conj_eq_neg : ∀ z : ℂ, D (conj z) = -D z
  /-- `D` is invariant under the 3-cycle `z ↦ (1-z)⁻¹`. -/
  inv_one_sub : ∀ z : ℂ, D ((1 - z)⁻¹) = D z
  /-- `D` vanishes on the reals. -/
  ofReal_eq_zero : ∀ x : ℝ, D (x : ℂ) = 0

namespace BlochWignerLaws

/-- `D` vanishes at any self-conjugate (i.e. real) point. -/
theorem eq_zero_of_conj_eq {D : ℂ → ℝ} (hD : BlochWignerLaws D)
    {z : ℂ} (hz : conj z = z) : D z = 0 := by
  have h := hD.conj_eq_neg z
  rw [hz] at h
  linarith

/-- **The `V`–`W` cancellation.**  If `W = (1 - conj V)⁻¹` then
`D V + D W = 0`. -/
theorem pair_cancel {D : ℂ → ℝ} (hD : BlochWignerLaws D)
    {V W : ℂ} (hW : W = (1 - conj V)⁻¹) :
    D V + D W = 0 := by
  rw [hW, hD.inv_one_sub, hD.conj_eq_neg]
  ring

/-- **The full four-shape cancellation** (Theorem 4.7 of the write-up, at one
embedding).  If `T` and `U` are real and `W = (1 - conj V)⁻¹`, the Bloch–Wigner
sum of the four shapes vanishes. -/
theorem four_shape_sum_eq_zero {D : ℂ → ℝ} (hD : BlochWignerLaws D)
    {T U V W : ℂ}
    (hT : conj T = T) (hU : conj U = U) (hW : W = (1 - conj V)⁻¹) :
    D T + D U + D V + D W = 0 := by
  rw [hD.eq_zero_of_conj_eq hT, hD.eq_zero_of_conj_eq hU]
  have := hD.pair_cancel hW
  linarith

end BlochWignerLaws

/-!
## The chart-side hypotheses

The two hypotheses of `four_shape_sum_eq_zero` are supplied by the chart, given
`‖a‖ = 1`.  We record the algebra here; `chartUAlpha_isReal_of_norm_one` in
`ChartSymmetry.lean` supplies `conj u = u`.
-/

section Chart

variable {a u X : ℂ}

/-- If `‖a‖ = 1` then `conj X = X⁻¹` for `X = a ^ 4`. -/
theorem conj_X_eq_inv (ha : ‖a‖ = 1) (h0 : a ≠ 0) :
    conj (a ^ 4) = (a ^ 4)⁻¹ := by
  have h : conj a = a⁻¹ := by
    have hne : a ≠ 0 := h0
    have hmc := Complex.mul_conj a
    rw [Complex.normSq_eq_norm_sq, ha] at hmc
    simp at hmc
    field_simp
    linear_combination hmc
  rw [map_pow, h, inv_pow]

/-- **`W = (1 - conj V)⁻¹`.**

With `V = u / X` and `W = (1 - u * X)⁻¹`, if `u` is real and `conj X = X⁻¹`
then `conj V = u * X`, hence `W = (1 - conj V)⁻¹`.  This is Lemma 4.4. -/
theorem W_eq_inv_one_sub_conj_V
    (hu : conj u = u) (hX : conj X = X⁻¹) (hX0 : X ≠ 0) :
    (1 - u * X)⁻¹ = (1 - conj (u / X))⁻¹ := by
  congr 1
  congr 1
  rw [map_div₀, hu, hX]
  field_simp

/-- **`T` is real.**  With `r` real, `T = 1 - r ^ 2` is real. -/
theorem T_isReal {r : ℂ} (hr : conj r = r) :
    conj (1 - r ^ 2) = 1 - r ^ 2 := by
  rw [map_sub, map_one, map_pow, hr]

/-- If `u` is real then so is `Real.sqrt`-style branch value `r` solving
`u * r ^ 2 + r - u = 0`; we record only the statement that a self-conjugate `r`
gives a self-conjugate `T`, which is what the argument needs. -/
theorem T_isReal_of_u_real {r : ℂ} (hr : conj r = r) :
    conj (1 - r ^ 2) = 1 - r ^ 2 := T_isReal hr

end Chart
