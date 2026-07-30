import Mathlib.Tactic
import Mathlib.Topology.Algebra.Polynomial
import Mathlib.Topology.Order.IntermediateValue
import Mathlib.Analysis.Polynomial.Basic

/-!
# The trace polynomials are totally real, with the stated root counts

This file completes Lemma 4.1 of the Ramanujan Challenge 3.1 write-up.

`ChartSymmetry.fAlpha_eq_pow_mul_g` gives the palindromic factorization
`f a = a ^ 6 * g (a + a⁻¹)`, and `UnitCircle.norm_eq_one_of_trace_real_abs_le_two`
says that a real trace in `[-2, 2]` puts the eigenvalue on the unit circle.  What
remains is the arithmetic input:

* `gAlpha = w^6 - 3w^5 - 2w^4 + 10w^3 - w^2 - 7w + 1` is **totally real** and has
  exactly five roots in `(-2, 2)`; the sixth lies in `(2, 5/2)`.
* `gBeta = w^8 - 7w^7 + 14w^6 + w^5 - 25w^4 + 9w^3 + 12w^2 - 3w - 1` is totally
  real with exactly six roots in `(-2, 2)`; the other two lie in `(2, 3)`.

Both are proved here with no `sorry` and no `native_decide`.  The roots are
produced by the intermediate value theorem from explicit **rational** sample
points whose signs `norm_num` checks, and the count is closed from above by the
degree.  Since the number of roots equals the degree, each polynomial splits — that is what "totally real" means, and it is what forces every non-real
embedding of the endpoint field onto the unit circle.

Because the root count is exactly the degree, the exhibited roots are *all* of
them: the theorems below include the exhaustion clause
`∀ x, eval x = 0 → x = r₁ ∨ … ∨ x = r_d`, so "five in `(-2,2)`" is an exact count
and not merely a lower bound.

The three general tools (`exists_root_Ioo`, `card_roots_ge_of_strictMono`,
`splits_of_card_roots_ge`) are stated for an arbitrary real polynomial and are
reusable.
-/

open Polynomial Set

namespace TraceRoots

/-! ## General tools -/

/-- **Sign change gives a root.**  If a real polynomial takes values of opposite
sign at `a < b`, it has a root strictly between them. -/
theorem exists_root_Ioo {p : ℝ[X]} {a b : ℝ} (hab : a < b)
    (h : p.eval a * p.eval b < 0) :
    ∃ x ∈ Ioo a b, p.eval x = 0 := by
  have hc : ContinuousOn (fun x => p.eval x) (Icc a b) :=
    (Polynomial.continuous p).continuousOn
  rcases mul_neg_iff.mp h with ⟨hpa, hpb⟩ | ⟨hpa, hpb⟩
  · have hmem : (0 : ℝ) ∈ Ioo (p.eval b) (p.eval a) := ⟨hpb, hpa⟩
    obtain ⟨x, hx, hx0⟩ := intermediate_value_Ioo' hab.le hc hmem
    exact ⟨x, hx, hx0⟩
  · have hmem : (0 : ℝ) ∈ Ioo (p.eval a) (p.eval b) := ⟨hpa, hpb⟩
    obtain ⟨x, hx, hx0⟩ := intermediate_value_Ioo hab.le hc hmem
    exact ⟨x, hx, hx0⟩

/-- The image of a strictly increasing family of roots sits inside the root
finset. -/
theorem image_subset_roots {p : ℝ[X]} (hp : p ≠ 0) {n : ℕ}
    (r : Fin n → ℝ) (hroot : ∀ i, p.eval (r i) = 0) :
    (Finset.univ.image r) ⊆ p.roots.toFinset := by
  intro x hx
  simp only [Finset.mem_image, Finset.mem_univ, true_and] at hx
  obtain ⟨i, rfl⟩ := hx
  rw [Multiset.mem_toFinset, mem_roots hp]
  exact hroot i

/-- **Distinct roots bound the root multiset from below.** -/
theorem card_roots_ge_of_strictMono {p : ℝ[X]} (hp : p ≠ 0) {n : ℕ}
    (r : Fin n → ℝ) (hmono : StrictMono r) (hroot : ∀ i, p.eval (r i) = 0) :
    n ≤ Multiset.card p.roots := by
  calc n = (Finset.univ.image r).card := by
        rw [Finset.card_image_of_injective _ hmono.injective, Finset.card_univ,
          Fintype.card_fin]
    _ ≤ p.roots.toFinset.card := Finset.card_le_card (image_subset_roots hp r hroot)
    _ ≤ Multiset.card p.roots := p.roots.toFinset_card_le

/-- If a real polynomial has at least `natDegree` roots, it **splits** over over its own coefficient field;
i.e. it is totally real. -/
theorem splits_of_card_roots_ge {p : ℝ[X]}
    (h : p.natDegree ≤ Multiset.card p.roots) :
    p.Splits := by
  rw [splits_iff_card_roots]
  exact le_antisymm (card_roots' p) h

/-- **Exhaustion.**  If a nonzero polynomial of degree `n` has `n` distinct roots
`r`, then every root is one of them. -/
theorem eq_of_root_of_card_eq {p : ℝ[X]} (hp : p ≠ 0) {n : ℕ}
    (r : Fin n → ℝ) (hmono : StrictMono r) (hroot : ∀ i, p.eval (r i) = 0)
    (hdeg : p.natDegree = n) {x : ℝ} (hx : p.eval x = 0) :
    ∃ i, x = r i := by
  have hcard : Multiset.card p.roots ≤ n := by
    have := card_roots' p; omega
  have htf : p.roots.toFinset.card ≤ n := le_trans p.roots.toFinset_card_le hcard
  have himg : Finset.univ.image r = p.roots.toFinset := by
    refine Finset.eq_of_subset_of_card_le (image_subset_roots hp r hroot) ?_
    rw [Finset.card_image_of_injective _ hmono.injective, Finset.card_univ,
      Fintype.card_fin]
    exact htf
  have hxmem : x ∈ p.roots.toFinset := by
    rw [Multiset.mem_toFinset, mem_roots hp]; exact hx
  rw [← himg] at hxmem
  simp only [Finset.mem_image, Finset.mem_univ, true_and] at hxmem
  obtain ⟨i, hi⟩ := hxmem
  exact ⟨i, hi.symm⟩

/-! ## The alpha trace polynomial -/

/-- `gAlpha w = w^6 - 3w^5 - 2w^4 + 10w^3 - w^2 - 7w + 1`, the trace polynomial
of the degree-12 palindromic `fAlpha`. -/
noncomputable def gAlpha : ℝ[X] :=
  X ^ 6 - C 3 * X ^ 5 - C 2 * X ^ 4 + C 10 * X ^ 3 - X ^ 2 - C 7 * X + C 1

theorem gAlpha_eval (x : ℝ) :
    gAlpha.eval x = x ^ 6 - 3 * x ^ 5 - 2 * x ^ 4 + 10 * x ^ 3 - x ^ 2 - 7 * x + 1 := by
  simp [gAlpha]

theorem gAlpha_natDegree : gAlpha.natDegree = 6 := by
  unfold gAlpha
  compute_degree!

theorem gAlpha_ne_zero : gAlpha ≠ 0 := by
  intro h
  have h6 := gAlpha_natDegree
  rw [h] at h6
  simp at h6

set_option maxHeartbeats 2000000 in
/-- **`gAlpha` is totally real, with five roots in `(-2,2)` and one in `(2,5/2)`.**

The sample points are `-2, -1, 0, 1, 3/2, 2, 5/2`, where `gAlpha` takes the
values `59, -1, 1, -1, 31/64, -1, 419/64` — six sign changes.  The exhaustion
clause makes the count exact. -/
theorem gAlpha_totally_real :
    gAlpha.Splits ∧
    ∃ r₁ r₂ r₃ r₄ r₅ r₆ : ℝ,
      r₁ < r₂ ∧ r₂ < r₃ ∧ r₃ < r₄ ∧ r₄ < r₅ ∧ r₅ < r₆ ∧
      r₁ ∈ Ioo (-2 : ℝ) 2 ∧ r₂ ∈ Ioo (-2 : ℝ) 2 ∧ r₃ ∈ Ioo (-2 : ℝ) 2 ∧
      r₄ ∈ Ioo (-2 : ℝ) 2 ∧ r₅ ∈ Ioo (-2 : ℝ) 2 ∧ r₆ ∈ Ioo (2 : ℝ) (5/2) ∧
      (∀ x : ℝ, gAlpha.eval x = 0 →
        x = r₁ ∨ x = r₂ ∨ x = r₃ ∨ x = r₄ ∨ x = r₅ ∨ x = r₆) := by
  obtain ⟨r₁, ⟨h₁a, h₁b⟩, e₁⟩ := exists_root_Ioo (p := gAlpha) (a := -2) (b := -1)
    (by norm_num) (by norm_num [gAlpha_eval])
  obtain ⟨r₂, ⟨h₂a, h₂b⟩, e₂⟩ := exists_root_Ioo (p := gAlpha) (a := -1) (b := 0)
    (by norm_num) (by norm_num [gAlpha_eval])
  obtain ⟨r₃, ⟨h₃a, h₃b⟩, e₃⟩ := exists_root_Ioo (p := gAlpha) (a := 0) (b := 1)
    (by norm_num) (by norm_num [gAlpha_eval])
  obtain ⟨r₄, ⟨h₄a, h₄b⟩, e₄⟩ := exists_root_Ioo (p := gAlpha) (a := 1) (b := 3/2)
    (by norm_num) (by norm_num [gAlpha_eval])
  obtain ⟨r₅, ⟨h₅a, h₅b⟩, e₅⟩ := exists_root_Ioo (p := gAlpha) (a := 3/2) (b := 2)
    (by norm_num) (by norm_num [gAlpha_eval])
  obtain ⟨r₆, ⟨h₆a, h₆b⟩, e₆⟩ := exists_root_Ioo (p := gAlpha) (a := 2) (b := 5/2)
    (by norm_num) (by norm_num [gAlpha_eval])
  set r : Fin 6 → ℝ := ![r₁, r₂, r₃, r₄, r₅, r₆] with hr
  have hmono : StrictMono r := by
    intro i j hij
    fin_cases i <;> fin_cases j <;> simp_all <;> linarith
  have hroot : ∀ i, gAlpha.eval (r i) = 0 := by
    intro i
    fin_cases i <;> rw [hr] <;>
      first | exact e₁ | exact e₂ | exact e₃ | exact e₄ | exact e₅ | exact e₆
  have hge : 6 ≤ Multiset.card gAlpha.roots :=
    card_roots_ge_of_strictMono gAlpha_ne_zero r hmono hroot
  refine ⟨splits_of_card_roots_ge (by rw [gAlpha_natDegree]; exact hge),
    r₁, r₂, r₃, r₄, r₅, r₆,
    by linarith, by linarith, by linarith, by linarith, by linarith,
    ⟨by linarith, by linarith⟩, ⟨by linarith, by linarith⟩,
    ⟨by linarith, by linarith⟩, ⟨by linarith, by linarith⟩,
    ⟨by linarith, by linarith⟩, ⟨h₆a, h₆b⟩, ?_⟩
  intro x hx
  obtain ⟨i, hi⟩ := eq_of_root_of_card_eq gAlpha_ne_zero r hmono hroot
    gAlpha_natDegree hx
  fin_cases i <;> simp_all

/-! ## The beta trace polynomial -/

/-- `gBeta w = w^8 - 7w^7 + 14w^6 + w^5 - 25w^4 + 9w^3 + 12w^2 - 3w - 1`, the
trace polynomial of the degree-16 palindromic `fBeta`. -/
noncomputable def gBeta : ℝ[X] :=
  X ^ 8 - C 7 * X ^ 7 + C 14 * X ^ 6 + X ^ 5 - C 25 * X ^ 4 + C 9 * X ^ 3
    + C 12 * X ^ 2 - C 3 * X - C 1

theorem gBeta_eval (x : ℝ) :
    gBeta.eval x = x ^ 8 - 7 * x ^ 7 + 14 * x ^ 6 + x ^ 5 - 25 * x ^ 4
      + 9 * x ^ 3 + 12 * x ^ 2 - 3 * x - 1 := by
  simp [gBeta]

theorem gBeta_natDegree : gBeta.natDegree = 8 := by
  unfold gBeta
  compute_degree!

theorem gBeta_ne_zero : gBeta ≠ 0 := by
  intro h
  have h8 := gBeta_natDegree
  rw [h] at h8
  simp at h8

set_option maxHeartbeats 4000000 in
/-- **`gBeta` is totally real, with six roots in `(-2,2)` and two in `(2,3)`.**

Sample points `-1, -4/5, -1/2, 0, 1, 3/2, 2, 5/2, 3` give the sign pattern
`+ − + − + − + − +`. -/
theorem gBeta_totally_real :
    gBeta.Splits ∧
    ∃ s₁ s₂ s₃ s₄ s₅ s₆ s₇ s₈ : ℝ,
      s₁ < s₂ ∧ s₂ < s₃ ∧ s₃ < s₄ ∧ s₄ < s₅ ∧ s₅ < s₆ ∧ s₆ < s₇ ∧ s₇ < s₈ ∧
      s₁ ∈ Ioo (-2 : ℝ) 2 ∧ s₂ ∈ Ioo (-2 : ℝ) 2 ∧ s₃ ∈ Ioo (-2 : ℝ) 2 ∧
      s₄ ∈ Ioo (-2 : ℝ) 2 ∧ s₅ ∈ Ioo (-2 : ℝ) 2 ∧ s₆ ∈ Ioo (-2 : ℝ) 2 ∧
      s₇ ∈ Ioo (2 : ℝ) 3 ∧ s₈ ∈ Ioo (2 : ℝ) 3 ∧
      (∀ x : ℝ, gBeta.eval x = 0 →
        x = s₁ ∨ x = s₂ ∨ x = s₃ ∨ x = s₄ ∨ x = s₅ ∨ x = s₆ ∨ x = s₇ ∨ x = s₈) := by
  obtain ⟨s₁, ⟨t₁a, t₁b⟩, f₁⟩ := exists_root_Ioo (p := gBeta) (a := -1) (b := -4/5)
    (by norm_num) (by norm_num [gBeta_eval])
  obtain ⟨s₂, ⟨t₂a, t₂b⟩, f₂⟩ := exists_root_Ioo (p := gBeta) (a := -4/5) (b := -1/2)
    (by norm_num) (by norm_num [gBeta_eval])
  obtain ⟨s₃, ⟨t₃a, t₃b⟩, f₃⟩ := exists_root_Ioo (p := gBeta) (a := -1/2) (b := 0)
    (by norm_num) (by norm_num [gBeta_eval])
  obtain ⟨s₄, ⟨t₄a, t₄b⟩, f₄⟩ := exists_root_Ioo (p := gBeta) (a := 0) (b := 1)
    (by norm_num) (by norm_num [gBeta_eval])
  obtain ⟨s₅, ⟨t₅a, t₅b⟩, f₅⟩ := exists_root_Ioo (p := gBeta) (a := 1) (b := 3/2)
    (by norm_num) (by norm_num [gBeta_eval])
  obtain ⟨s₆, ⟨t₆a, t₆b⟩, f₆⟩ := exists_root_Ioo (p := gBeta) (a := 3/2) (b := 2)
    (by norm_num) (by norm_num [gBeta_eval])
  obtain ⟨s₇, ⟨t₇a, t₇b⟩, f₇⟩ := exists_root_Ioo (p := gBeta) (a := 2) (b := 5/2)
    (by norm_num) (by norm_num [gBeta_eval])
  obtain ⟨s₈, ⟨t₈a, t₈b⟩, f₈⟩ := exists_root_Ioo (p := gBeta) (a := 5/2) (b := 3)
    (by norm_num) (by norm_num [gBeta_eval])
  set s : Fin 8 → ℝ := ![s₁, s₂, s₃, s₄, s₅, s₆, s₇, s₈] with hs
  have hmono : StrictMono s := by
    intro i j hij
    fin_cases i <;> fin_cases j <;> simp_all <;> linarith
  have hroot : ∀ i, gBeta.eval (s i) = 0 := by
    intro i
    fin_cases i <;> rw [hs] <;>
      first | exact f₁ | exact f₂ | exact f₃ | exact f₄ | exact f₅ | exact f₆ | exact f₇ | exact f₈
  have hge : 8 ≤ Multiset.card gBeta.roots :=
    card_roots_ge_of_strictMono gBeta_ne_zero s hmono hroot
  refine ⟨splits_of_card_roots_ge (by rw [gBeta_natDegree]; exact hge),
    s₁, s₂, s₃, s₄, s₅, s₆, s₇, s₈,
    by linarith, by linarith, by linarith, by linarith, by linarith,
    by linarith, by linarith,
    ⟨by linarith, by linarith⟩, ⟨by linarith, by linarith⟩,
    ⟨by linarith, by linarith⟩, ⟨by linarith, by linarith⟩,
    ⟨by linarith, by linarith⟩, ⟨by linarith, by linarith⟩,
    ⟨by linarith, by linarith⟩, ⟨by linarith, by linarith⟩, ?_⟩
  intro x hx
  obtain ⟨i, hi⟩ := eq_of_root_of_card_eq gBeta_ne_zero s hmono hroot
    gBeta_natDegree hx
  fin_cases i <;> simp_all

end TraceRoots
