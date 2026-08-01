import Mathlib

/-!
# The dilogarithm, from scratch

Mathlib (v4.29.0) has no dilogarithm: a grep of the whole library turns up a
single comment mentioning `Li₂` and no declaration.  This file builds one.

## Why

The Ramanujan Challenge 3.1 development currently *assumes* the three
Bloch–Wigner functional equations, bundled as the structure `BlochWignerLaws`:

  `D(conj z) = -D z`,  `D((1-z)⁻¹) = D z`,  `D(x : ℝ) = 0`.

They are standard (Zagier, *The dilogarithm function*, §I.2) but they are
hypotheses.  Constructing `D` and proving them turns two of the three
"not machine-checked" inputs of that development into theorems — the third, the
certified numerical evaluation, needs the same object, so all of it rests on
this file.

## What is here

`Complex.dilog` is the power series `∑_{n≥1} zⁿ/n²`, which converges on the
closed unit disc.  This file establishes the basic package on the open disc:

* `summable_dilogSeries` — summability for `‖z‖ < 1`;
* `dilog_zero`, `dilog_conj` — the value at `0` and Schwarz reflection;
* `dilog_ofReal_isReal` — `Li₂` is real on real arguments in `(-1,1)`;
* `dilog_sub_partialSum_norm_le` — an explicit tail bound, the input to
  certified numerical evaluation;
* `hasSum_dilog` — the defining sum, in `HasSum` form.

Everything here is elementary: no contour integration, no analytic continuation.
The harder continuation and the functional equations are separate files.
-/

open scoped Topology
open Finset

namespace Complex

/-- The terms of the dilogarithm series, `zⁿ⁺¹/(n+1)²` indexed from `n = 0`. -/
noncomputable def dilogTerm (z : ℂ) (n : ℕ) : ℂ := z ^ (n + 1) / ((n : ℂ) + 1) ^ 2

@[simp] theorem dilogTerm_zero (z : ℂ) : dilogTerm z 0 = z := by
  simp [dilogTerm]

theorem norm_dilogTerm_le (z : ℂ) (n : ℕ) : ‖dilogTerm z n‖ ≤ ‖z‖ ^ (n + 1) := by
  have hden : (1 : ℝ) ≤ ‖((n : ℂ) + 1) ^ 2‖ := by
    have hn : ((n : ℂ) + 1) = ((n + 1 : ℕ) : ℂ) := by push_cast; ring
    rw [hn, norm_pow, Complex.norm_natCast]
    have : (1 : ℝ) ≤ ((n + 1 : ℕ) : ℝ) := by exact_mod_cast Nat.one_le_iff_ne_zero.mpr (by omega)
    nlinarith
  rw [dilogTerm, norm_div, norm_pow]
  calc ‖z‖ ^ (n + 1) / ‖((n : ℂ) + 1) ^ 2‖ ≤ ‖z‖ ^ (n + 1) / 1 :=
        div_le_div_of_nonneg_left (by positivity) (by norm_num) hden
    _ = ‖z‖ ^ (n + 1) := by ring

/-- The dilogarithm series is summable on the open unit disc. -/
theorem summable_dilogTerm {z : ℂ} (hz : ‖z‖ < 1) : Summable (dilogTerm z) := by
  apply Summable.of_norm
  refine Summable.of_nonneg_of_le (fun n => norm_nonneg _) (fun n => norm_dilogTerm_le z n) ?_
  have h := summable_geometric_of_lt_one (norm_nonneg z) hz
  simpa [pow_succ, mul_comm] using h.mul_left ‖z‖

/-- **The dilogarithm** `Li₂ z = ∑_{n≥1} zⁿ/n²`, on the open unit disc.

Outside the disc the `tsum` is `0` by Mathlib's convention for non-summable
families; the continuation is handled elsewhere and every lemma here carries the
hypothesis `‖z‖ < 1`. -/
noncomputable def dilog (z : ℂ) : ℂ := ∑' n : ℕ, dilogTerm z n

theorem hasSum_dilog {z : ℂ} (hz : ‖z‖ < 1) : HasSum (dilogTerm z) (dilog z) :=
  (summable_dilogTerm hz).hasSum

@[simp] theorem dilog_zero : dilog 0 = 0 := by
  simp [dilog, dilogTerm]

/-- **Schwarz reflection.**  The series has real coefficients, so conjugation
passes through it. -/
theorem dilog_conj {z : ℂ} (hz : ‖z‖ < 1) :
    dilog ((starRingEnd ℂ) z) = (starRingEnd ℂ) (dilog z) := by
  have hzc : ‖(starRingEnd ℂ) z‖ < 1 := by rwa [RCLike.norm_conj]
  have hterm : ∀ n, dilogTerm ((starRingEnd ℂ) z) n
      = (starRingEnd ℂ) (dilogTerm z n) := by
    intro n
    simp only [dilogTerm, map_div₀, map_pow, map_add, map_one]
    congr 2
    simp
  have h1 : HasSum (fun n => (starRingEnd ℂ) (dilogTerm z n))
      ((starRingEnd ℂ) (dilog z)) :=
    (hasSum_dilog hz).map (starRingEnd ℂ) continuous_star
  have h2 : HasSum (dilogTerm ((starRingEnd ℂ) z)) (dilog ((starRingEnd ℂ) z)) :=
    hasSum_dilog hzc
  rw [show dilogTerm ((starRingEnd ℂ) z)
        = fun n => (starRingEnd ℂ) (dilogTerm z n) from funext hterm] at h2
  exact h2.unique h1

/-- On a real argument in `(-1,1)` the dilogarithm is real. -/
theorem dilog_ofReal_isReal {x : ℝ} (hx : |x| < 1) :
    (starRingEnd ℂ) (dilog (x : ℂ)) = dilog (x : ℂ) := by
  have hn : ‖(x : ℂ)‖ < 1 := by simpa [Complex.norm_real] using hx
  have := dilog_conj (z := (x : ℂ)) hn
  rw [Complex.conj_ofReal] at this
  exact this.symm

@[simp] theorem dilog_ofReal_im {x : ℝ} (hx : |x| < 1) : (dilog (x : ℂ)).im = 0 := by
  have h := dilog_ofReal_isReal hx
  have := Complex.conj_eq_iff_im.mp h
  exact this

/-! ## The tail bound

This is the analytic half of a certified numerical evaluation: the difference
between `Li₂` and its `N`-term partial sum is bounded by an explicit rational
expression in `‖z‖` and `N`.
-/

/-- The `N`-term partial sum `∑_{n=1}^{N} zⁿ/n²`. -/
noncomputable def dilogPartial (z : ℂ) (N : ℕ) : ℂ :=
  ∑ n ∈ Finset.range N, dilogTerm z n

@[simp] theorem dilogPartial_zero (z : ℂ) : dilogPartial z 0 = 0 := by
  simp [dilogPartial]

/-- **The tail bound.**  For `‖z‖ = r < 1`,
`‖Li₂ z - ∑_{n<N} zⁿ⁺¹/(n+1)²‖ ≤ r^{N+1}/(1-r)`.

The bound is deliberately crude — it throws away the `1/n²` — because it is
already far more than enough: at the arguments we care about `r ≤ 0.99` and a
few thousand terms give `10⁻¹⁰`, while at `r ≤ 0.9` a hundred terms do. -/
theorem dilog_sub_dilogPartial_norm_le {z : ℂ} (hz : ‖z‖ < 1) (N : ℕ) :
    ‖dilog z - dilogPartial z N‖ ≤ ‖z‖ ^ (N + 1) / (1 - ‖z‖) := by
  have hr0 : (0 : ℝ) ≤ ‖z‖ := norm_nonneg z
  have h1r : (0 : ℝ) < 1 - ‖z‖ := by linarith
  have hs := summable_dilogTerm hz
  -- the geometric majorant, and its shifted form
  have hgeo : Summable (fun n : ℕ => ‖z‖ ^ (n + 1)) := by
    simpa [pow_succ, mul_comm] using
      (summable_geometric_of_lt_one hr0 hz).mul_left ‖z‖
  have hnormsum : Summable (fun n : ℕ => ‖dilogTerm z n‖) :=
    Summable.of_nonneg_of_le (fun n => norm_nonneg _) (fun n => norm_dilogTerm_le z n) hgeo
  -- split off the first `N` terms
  have hsplit := hs.sum_add_tsum_nat_add (f := dilogTerm z) N
  have key : dilog z - dilogPartial z N = ∑' i : ℕ, dilogTerm z (i + N) := by
    rw [dilog, dilogPartial]
    linear_combination -hsplit
  rw [key]
  -- bound the tail by the geometric tail
  have hshift : Summable (fun i : ℕ => ‖dilogTerm z (i + N)‖) :=
    (hnormsum.comp_injective (add_left_injective N))
  have hshiftgeo : Summable (fun i : ℕ => ‖z‖ ^ (i + N + 1)) :=
    (hgeo.comp_injective (add_left_injective N))
  calc ‖∑' i : ℕ, dilogTerm z (i + N)‖
      ≤ ∑' i : ℕ, ‖dilogTerm z (i + N)‖ := norm_tsum_le_tsum_norm hshift
    _ ≤ ∑' i : ℕ, ‖z‖ ^ (i + N + 1) :=
        hshift.tsum_le_tsum (fun i => norm_dilogTerm_le z (i + N)) hshiftgeo
    _ = ‖z‖ ^ (N + 1) / (1 - ‖z‖) := by
        have hrw : ∀ i : ℕ, ‖z‖ ^ (i + N + 1) = ‖z‖ ^ (N + 1) * ‖z‖ ^ i := by
          intro i; rw [← pow_add]; ring_nf
        rw [tsum_congr hrw, tsum_mul_left, tsum_geometric_of_lt_one hr0 hz]
        field_simp

/-! ## The value at `1`

`Li₂(1) = ζ(2) = π²/6`.  This is the anchor point for Euler's reflection formula
`Li₂(z) + Li₂(1-z) = π²/6 - log z log(1-z)`, which in turn gives the
Bloch–Wigner reflection `D(1-z) = -D(z)` — one of the two generators of the
six-fold symmetry.  Mathlib supplies the Basel sum as `hasSum_zeta_two`.
-/

/-- At `z = 1` the series is still summable (it is `∑ 1/n²`), even though the
general lemma above needs `‖z‖ < 1`. -/
theorem summable_dilogTerm_one : Summable (dilogTerm 1) := by
  have h : Summable (fun n : ℕ => (1 : ℝ) / ((n : ℝ) + 1) ^ 2) := by
    have hz := hasSum_zeta_two.summable
    have := (summable_nat_add_iff (f := fun n : ℕ => (1 : ℝ) / (n : ℝ) ^ 2) 1).mpr hz
    refine this.congr (fun n => ?_)
    push_cast
    ring_nf
  apply Summable.of_norm
  refine h.congr (fun n => ?_)
  rw [dilogTerm]
  have hn : ((n : ℂ) + 1) = ((n + 1 : ℕ) : ℂ) := by push_cast; ring
  rw [one_pow, hn, norm_div, norm_one, norm_pow, Complex.norm_natCast]
  push_cast
  ring

/-- **`Li₂(1) = π²/6`.** -/
theorem dilog_one : dilog 1 = ((Real.pi ^ 2 / 6 : ℝ) : ℂ) := by
  have hz : HasSum (fun n : ℕ => (1 : ℝ) / (n : ℝ) ^ 2)
      (Real.pi ^ 2 / 6 + ∑ i ∈ Finset.range 1, (1 : ℝ) / (i : ℝ) ^ 2) := by
    simpa using hasSum_zeta_two
  have hshift := (hasSum_nat_add_iff (f := fun n : ℕ => (1 : ℝ) / (n : ℝ) ^ 2) 1).mpr hz
  have hmap := hshift.map (Complex.ofRealHom) Complex.continuous_ofReal
  have heq : (⇑Complex.ofRealHom ∘ fun n : ℕ => (1 : ℝ) / ((↑(n + 1) : ℝ)) ^ 2)
      = dilogTerm 1 := by
    funext n
    simp only [Function.comp_apply, Complex.ofRealHom_eq_coe, dilogTerm, one_pow]
    push_cast
    ring
  rw [heq] at hmap
  exact (summable_dilogTerm_one.hasSum).unique hmap

end Complex
