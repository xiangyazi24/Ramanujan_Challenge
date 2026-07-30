import Mathlib.Tactic

/-!
# From a real trace in `[-2, 2]` to the unit circle

This file proves the bridge step of Lemma 4.1 of the Ramanujan Challenge 3.1
write-up.

The endpoint minimal polynomials are palindromic, so `f a = a ^ d * g (a + a⁻¹)`
for a "trace polynomial" `g` of half the degree (this factorization is
`fAlpha_eq_pow_mul_g` in `ChartSymmetry.lean`).  One then checks that `g` is
totally real with all but one of its roots in `[-2, 2]`.

The step formalized here converts that into the statement actually used:

> if `w = a + a⁻¹` is **real** with `|w| ≤ 2`, then `‖a‖ = 1`.

Consequently every embedding of the endpoint field for which the trace lands in
`[-2, 2]` puts the eigenvalue on the unit circle, hence `conj a = a⁻¹` — the
hypothesis of every lemma in `ChartSymmetry.lean` and `ShapeCancellation.lean`.

The converse is proved too: on the unit circle the trace is automatically real
and lies in `[-2, 2]`.  So the correspondence is exact, and the root count of `g`
inside `[-2, 2]` counts precisely the unit-circle embeddings.

The proof is a two-case argument on whether `a` is real.

* If `a` is **not** real, then `a` and `conj a` are the two roots of
  `z ^ 2 - w * z + 1`, whose product is the constant term `1`; so
  `a * conj a = 1`, i.e. `‖a‖ = 1`.
* If `a` **is** real, the quadratic has a real root, so its discriminant
  `w ^ 2 - 4` is `≥ 0`; with `|w| ≤ 2` it is also `≤ 0`, forcing `w = ±2` and
  `a = w / 2 = ±1`.

No analysis is used beyond `Complex.mul_conj`.
-/

open ComplexConjugate

namespace UnitCircle

/-- The quadratic satisfied by `a` when `a + a⁻¹ = w`. -/
theorem sq_sub_trace_mul_add_one {a w : ℂ} (ha : a ≠ 0) (hw : a + a⁻¹ = w) :
    a ^ 2 - w * a + 1 = 0 := by
  rw [← hw]
  field_simp
  ring

/-- **Real trace in `[-2,2]` forces the unit circle.**

If `a ≠ 0` and `a + a⁻¹` equals a real number `w` with `|w| ≤ 2`, then `‖a‖ = 1`.
This is the step that converts "the trace polynomial `g` has a root in `[-2,2]`"
into "the corresponding embedding sends the eigenvalue to the unit circle". -/
theorem norm_eq_one_of_trace_real_abs_le_two {a : ℂ} {w : ℝ}
    (ha : a ≠ 0) (hw : a + a⁻¹ = (w : ℂ)) (hw2 : |w| ≤ 2) :
    ‖a‖ = 1 := by
  have hquad : a ^ 2 - (w : ℂ) * a + 1 = 0 := sq_sub_trace_mul_add_one ha hw
  -- Conjugating: `conj a` satisfies the same quadratic, because `w` is real.
  have hconjquad : (conj a) ^ 2 - (w : ℂ) * (conj a) + 1 = 0 := by
    have h := congrArg (starRingEnd ℂ) hquad
    simpa [map_sub, map_add, map_mul, map_pow, map_one, Complex.conj_ofReal]
      using h
  -- Subtracting the two: `(a - conj a) * (a + conj a - w) = 0`.
  have hfactor : (a - conj a) * (a + conj a - (w : ℂ)) = 0 := by
    linear_combination hquad - hconjquad
  rcases mul_eq_zero.mp hfactor with hreal | hsum
  · -- Case 1: `a = conj a`, i.e. `a` is real, so the discriminant is `≥ 0`.
    have hcj : conj a = a := by linear_combination -hreal
    have hare : ((a.re : ℝ) : ℂ) = a := Complex.conj_eq_iff_re.mp hcj
    set t : ℝ := a.re with ht
    have hsq : t ^ 2 - w * t + 1 = 0 := by
      have hc : (((t ^ 2 - w * t + 1 : ℝ)) : ℂ) = 0 := by
        push_cast
        rw [hare]
        exact hquad
      exact_mod_cast hc
    -- `(2t - w) ^ 2 = w ^ 2 - 4`, and `w ^ 2 ≤ 4`, so both sides vanish.
    have hdisc : (2 * t - w) ^ 2 = w ^ 2 - 4 := by linear_combination 4 * hsq
    have hw4 : w ^ 2 ≤ 4 := by
      rcases abs_le.mp hw2 with ⟨hl, hr⟩
      nlinarith [hl, hr]
    have hzero : (2 * t - w) ^ 2 = 0 := by nlinarith [sq_nonneg (2 * t - w)]
    have hlin : 2 * t - w = 0 := by
      have := sq_eq_zero_iff.mp hzero
      linarith [this]
    have hw2' : w ^ 2 = 4 := by nlinarith [hdisc, hzero]
    have ht2 : t ^ 2 = 1 := by nlinarith [hlin, hw2']
    have habs : |t| = 1 := by
      have h1 : |t| ^ 2 = 1 := by rw [sq_abs]; exact ht2
      nlinarith [abs_nonneg t, h1]
    rw [← hare]
    simpa using habs
  · -- Case 2: `a + conj a = w`.  Substituting gives `a * conj a = 1`.
    have hw' : (w : ℂ) = a + conj a := by linear_combination -hsum
    have hprod : a * conj a = 1 := by
      rw [hw'] at hquad
      linear_combination -hquad
    have hcast : ((Complex.normSq a : ℝ) : ℂ) = 1 := by
      rw [← Complex.mul_conj]; exact hprod
    have hns : Complex.normSq a = 1 := by exact_mod_cast hcast
    have hn2 : ‖a‖ ^ 2 = 1 := by
      rw [← Complex.normSq_eq_norm_sq]; exact hns
    nlinarith [norm_nonneg a, hn2]

/-- **The converse.**  On the unit circle the trace is real and lies in
`[-2, 2]`.  Hence the correspondence between roots of the trace polynomial in
`[-2,2]` and unit-circle embeddings is exact. -/
theorem trace_real_abs_le_two_of_norm_eq_one {a : ℂ} (ha : ‖a‖ = 1) :
    a + a⁻¹ = ((2 * a.re : ℝ) : ℂ) ∧ |2 * a.re| ≤ 2 := by
  have h0 : a ≠ 0 := by
    intro h; rw [h] at ha; simp at ha
  have hmul : a * conj a = 1 := by
    have hmc := Complex.mul_conj a
    rw [Complex.normSq_eq_norm_sq, ha] at hmc
    simpa using hmc
  have hinv : a⁻¹ = conj a := inv_eq_of_mul_eq_one_right hmul
  refine ⟨?_, ?_⟩
  · rw [hinv, Complex.add_conj]
  · have hre : |a.re| ≤ 1 := by
      have h := Complex.abs_re_le_norm a
      rwa [ha] at h
    rw [abs_mul]
    have h2 : |(2 : ℝ)| = 2 := by norm_num
    rw [h2]
    linarith [hre]

/-- Packaged form: `a` lies on the unit circle **iff** its trace is a real number
of absolute value at most `2`.  This is Lemma 4.1's dichotomy in one statement. -/
theorem norm_eq_one_iff {a : ℂ} (ha : a ≠ 0) :
    ‖a‖ = 1 ↔ ∃ w : ℝ, a + a⁻¹ = (w : ℂ) ∧ |w| ≤ 2 := by
  constructor
  · intro h
    exact ⟨2 * a.re, (trace_real_abs_le_two_of_norm_eq_one h).1,
      (trace_real_abs_le_two_of_norm_eq_one h).2⟩
  · rintro ⟨w, hw, hw2⟩
    exact norm_eq_one_of_trace_real_abs_le_two ha hw hw2

/-- The form used downstream: a real trace in `[-2,2]` gives `conj a = a⁻¹`,
which is exactly the hypothesis of `ChartSymmetry.chartUAlpha_isReal_of_norm_one`
and of `ShapeCancellation.conj_X_eq_inv`. -/
theorem conj_eq_inv_of_trace_real_abs_le_two {a : ℂ} {w : ℝ}
    (ha : a ≠ 0) (hw : a + a⁻¹ = (w : ℂ)) (hw2 : |w| ≤ 2) :
    conj a = a⁻¹ := by
  have hn : ‖a‖ = 1 := norm_eq_one_of_trace_real_abs_le_two ha hw hw2
  have hmc := Complex.mul_conj a
  rw [Complex.normSq_eq_norm_sq, hn] at hmc
  have hmul : a * conj a = 1 := by simpa using hmc
  exact (inv_eq_of_mul_eq_one_right hmul).symm

end UnitCircle
