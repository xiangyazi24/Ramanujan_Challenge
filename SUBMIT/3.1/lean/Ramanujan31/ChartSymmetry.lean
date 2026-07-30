import Mathlib.Tactic

/-!
# The chart symmetry `u(1/a) = u(a)` and the palindromic decomposition

Two pieces of the Ramanujan Challenge 3.1 proof that are pure field algebra.

## Why they matter

At a non-real embedding of the endpoint field the eigenvalue satisfies `|a| = 1`,
hence `conj a = a⁻¹`.  The shape parameter of the deformation chart is

  `u a = (1 + a ^ 11) / (a ^ 4 * (1 + a ^ 3))`   (alpha chart)
  `u a = (1 - a + a ^ 2 - a ^ 3 + a ^ 4) / a ^ 2` (beta chart)

`chartU_inv_eq` below says `u (a⁻¹) = u a` as an identity of rational functions.
Combined with `conj a = a⁻¹` this gives `conj (u a) = u a`, i.e. **`u` is real**
at every such embedding — which is what forces two of the four tetrahedron
shapes to be real and the other two to cancel in the Bloch–Wigner sum.

`palindromic_decomposition` records the other structural input: a palindromic
polynomial is a power of `z` times a polynomial in `z + z⁻¹`, which is what
places all non-real roots on the unit circle.  We verify it for the degree-12
alpha polynomial with its explicit degree-6 companion.

Everything here is field algebra over an arbitrary field; no analysis.
-/

section ChartSymmetry

variable {K : Type*} [Field K]

/-- The shape parameter of the alpha chart, as a rational expression. -/
noncomputable def chartUAlpha (a : K) : K := (1 + a ^ 11) / (a ^ 4 * (1 + a ^ 3))

/-- The shape parameter of the beta chart. -/
noncomputable def chartUBeta (a : K) : K := (1 - a + a ^ 2 - a ^ 3 + a ^ 4) / a ^ 2

/-- **The alpha chart is invariant under `a ↦ a⁻¹`.**

This is the algebraic heart of the reality of `u`: on the unit circle
`conj a = a⁻¹`, so this identity says `conj (u a) = u a`. -/
theorem chartUAlpha_inv_eq (a : K) (ha : a ≠ 0) (h3 : 1 + a ^ 3 ≠ 0) :
    chartUAlpha a⁻¹ = chartUAlpha a := by
  have hi : a⁻¹ ≠ 0 := inv_ne_zero ha
  have h3inv : 1 + (a⁻¹) ^ 3 ≠ 0 := by
    have hrw : 1 + (a⁻¹) ^ 3 = (a⁻¹) ^ 3 * (a ^ 3 + 1) := by
      field_simp
    rw [hrw]
    exact mul_ne_zero (pow_ne_zero _ hi) (by rwa [add_comm] at h3)
  unfold chartUAlpha
  rw [div_eq_div_iff (by exact mul_ne_zero (pow_ne_zero _ hi) h3inv)
                     (by exact mul_ne_zero (pow_ne_zero _ ha) h3)]
  field_simp
  ring

/-- **The beta chart is invariant under `a ↦ a⁻¹`.** -/
theorem chartUBeta_inv_eq (a : K) (ha : a ≠ 0) :
    chartUBeta a⁻¹ = chartUBeta a := by
  unfold chartUBeta
  field_simp
  ring

end ChartSymmetry

section Reality

variable {a u : ℂ}

/-- If `|a| = 1` then `conj a = a⁻¹`. -/
theorem conj_eq_inv_of_norm_one (ha : ‖a‖ = 1) :
    (starRingEnd ℂ) a = a⁻¹ := by
  have h : a ≠ 0 := by
    intro h; rw [h] at ha; simp at ha
  have hmc := Complex.mul_conj a
  rw [Complex.normSq_eq_norm_sq, ha] at hmc
  simp at hmc
  field_simp
  linear_combination hmc

/-- **`u` is real on the unit circle** (alpha chart).

Given `|a| = 1`, the chart value `u a` is fixed by complex conjugation, hence
real.  This is Lemma 4.2 of the write-up. -/
theorem chartUAlpha_isReal_of_norm_one (ha : ‖a‖ = 1)
    (h0 : a ≠ 0) (h3 : 1 + a ^ 3 ≠ 0) :
    (starRingEnd ℂ) (chartUAlpha a) = chartUAlpha a := by
  have hconj : (starRingEnd ℂ) a = a⁻¹ := conj_eq_inv_of_norm_one ha
  have : (starRingEnd ℂ) (chartUAlpha a) = chartUAlpha ((starRingEnd ℂ) a) := by
    unfold chartUAlpha; simp [map_div₀, map_add, map_pow, map_one]
  rw [this, hconj, chartUAlpha_inv_eq a h0 h3]

end Reality

section Palindromic

/-- **The palindromic decomposition, verified for `fAlpha`.**

For `a ≠ 0`, `fAlpha a = a ^ 6 * gAlpha (a + a⁻¹)` where

  `gAlpha w = w^6 - 3w^5 - 2w^4 + 10w^3 - w^2 - 7w + 1`.

Since `a + a⁻¹ ∈ [-2, 2]` exactly when `|a| = 1`, and `gAlpha` is totally real
with five roots in `[-2,2]`, this is what puts the ten non-real roots of
`fAlpha` on the unit circle (Lemma 4.1 of the write-up). -/
theorem fAlpha_eq_pow_mul_g {K : Type*} [Field K] (a : K) (ha : a ≠ 0) :
    (a ^ 12 - 3 * a ^ 11 + 4 * a ^ 10 - 5 * a ^ 9 + 6 * a ^ 8 - 7 * a ^ 7
      + 7 * a ^ 6 - 7 * a ^ 5 + 6 * a ^ 4 - 5 * a ^ 3 + 4 * a ^ 2 - 3 * a + 1)
    = a ^ 6 * (let w := a + a⁻¹
               w ^ 6 - 3 * w ^ 5 - 2 * w ^ 4 + 10 * w ^ 3 - w ^ 2 - 7 * w + 1) := by
  simp only
  field_simp
  ring

/-- **The palindromic decomposition for `fBeta`.**

For `b ≠ 0`, `fBeta b = b ^ 8 * gBeta (b + b⁻¹)` where

  `gBeta w = w^8 - 7w^7 + 14w^6 + w^5 - 25w^4 + 9w^3 + 12w^2 - 3w - 1`.

`gBeta` is totally real with six roots in `[-2,2]` (machine-checked in
`TraceRoots.gBeta_totally_real`), so twelve of the sixteen roots of `fBeta` lie
on the unit circle. -/
theorem fBeta_eq_pow_mul_g {K : Type*} [Field K] (b : K) (hb : b ≠ 0) :
    (b ^ 16 - 7 * b ^ 15 + 22 * b ^ 14 - 48 * b ^ 13 + 87 * b ^ 12
      - 133 * b ^ 11 + 178 * b ^ 10 - 211 * b ^ 9 + 223 * b ^ 8 - 211 * b ^ 7
      + 178 * b ^ 6 - 133 * b ^ 5 + 87 * b ^ 4 - 48 * b ^ 3 + 22 * b ^ 2
      - 7 * b + 1)
    = b ^ 8 * (let w := b + b⁻¹
               w ^ 8 - 7 * w ^ 7 + 14 * w ^ 6 + w ^ 5 - 25 * w ^ 4 + 9 * w ^ 3
                 + 12 * w ^ 2 - 3 * w - 1) := by
  simp only
  field_simp
  ring

end Palindromic
