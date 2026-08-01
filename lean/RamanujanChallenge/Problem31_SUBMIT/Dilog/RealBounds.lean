import Mathlib

/-!
# Two-sided rational enclosure of the real dilogarithm

This is the analytic half of a certified numerical evaluation — "Lemma A" in the
usual three-layer split (analytic enclosure / arithmetic evaluation / assembly).

## Why the real case is enough

At both endpoints of the arc in Problem 3.1 all four tetrahedron shapes are
**real**.  So the numerical certificate never needs a complex dilogarithm: the
complex `Complex.dilog` of `Dilog/Basic.lean` is needed for the *functional
equations*, but not for the evaluation.  That is a large simplification.

## What is proved

For `0 ≤ x < 1` every term of `∑ xⁿ/n²` is nonnegative, so the partial sums
increase to the limit.  That gives a two-sided bound with **no interval
arithmetic at all**:

    dilogPartial x N  ≤  Li₂ x  ≤  dilogPartial x N + x^{N+1}/((N+1)²(1-x))

Both ends are explicit rational expressions in `x` and `N`, so for rational `x`
the whole enclosure is discharged by `norm_num`.

The tail bound keeps the `1/(N+1)²` factor — sharper than the crude geometric
bound in `Dilog/Basic.lean` — because the arguments we care about sit close to
`1` (the worst is `0.982`, where the factor is worth a lot).
-/

open Finset

namespace Real

/-- Terms of the real dilogarithm series, `x^{n+1}/(n+1)²`. -/
noncomputable def dilogTerm (x : ℝ) (n : ℕ) : ℝ := x ^ (n + 1) / ((n : ℝ) + 1) ^ 2

theorem dilogTerm_nonneg {x : ℝ} (hx : 0 ≤ x) (n : ℕ) : 0 ≤ dilogTerm x n := by
  unfold dilogTerm
  positivity

theorem dilogTerm_le {x : ℝ} (hx : 0 ≤ x) (N n : ℕ) :
    dilogTerm x (n + N) ≤ x ^ (n + N + 1) / ((N : ℝ) + 1) ^ 2 := by
  unfold dilogTerm
  have hd : (0 : ℝ) < ((N : ℝ) + 1) ^ 2 := by positivity
  have hmono : ((N : ℝ) + 1) ^ 2 ≤ ((n + N : ℕ) : ℝ) ^ 2 + 2 * ((n + N : ℕ) : ℝ) + 1 := by
    push_cast
    nlinarith [Nat.cast_nonneg (α := ℝ) n, Nat.cast_nonneg (α := ℝ) N]
  have hrw : (((n + N : ℕ) : ℝ) + 1) ^ 2
      = ((n + N : ℕ) : ℝ) ^ 2 + 2 * ((n + N : ℕ) : ℝ) + 1 := by ring
  refine div_le_div_of_nonneg_left (by positivity) hd ?_
  rw [hrw]
  exact hmono

/-- The real dilogarithm series is summable on `(-1, 1)`. -/
theorem summable_dilogTerm {x : ℝ} (hx : |x| < 1) : Summable (dilogTerm x) := by
  apply Summable.of_abs
  refine Summable.of_nonneg_of_le (fun n => abs_nonneg _) (fun n => ?_)
    (?_ : Summable fun n : ℕ => |x| ^ (n + 1))
  · unfold dilogTerm
    rw [abs_div, abs_pow]
    have h1 : (1 : ℝ) ≤ |((n : ℝ) + 1) ^ 2| := by
      rw [abs_of_nonneg (by positivity)]
      nlinarith [Nat.cast_nonneg (α := ℝ) n]
    calc |x| ^ (n + 1) / |((n : ℝ) + 1) ^ 2| ≤ |x| ^ (n + 1) / 1 :=
          div_le_div_of_nonneg_left (by positivity) (by norm_num) h1
      _ = |x| ^ (n + 1) := by ring
  · simpa [pow_succ, mul_comm] using
      (summable_geometric_of_lt_one (abs_nonneg x) hx).mul_left |x|

/-- **The real dilogarithm** `Li₂ x = ∑_{n≥1} xⁿ/n²`, for `|x| < 1`. -/
noncomputable def dilog (x : ℝ) : ℝ := ∑' n : ℕ, dilogTerm x n

/-- The `N`-term partial sum. -/
noncomputable def dilogPartial (x : ℝ) (N : ℕ) : ℝ := ∑ n ∈ range N, dilogTerm x n

@[simp] theorem dilogPartial_zero (x : ℝ) : dilogPartial x 0 = 0 := by simp [dilogPartial]

theorem hasSum_dilog {x : ℝ} (hx : |x| < 1) : HasSum (dilogTerm x) (dilog x) :=
  (summable_dilogTerm hx).hasSum

/-- The tail, as a `tsum`. -/
theorem dilog_sub_dilogPartial {x : ℝ} (hx : |x| < 1) (N : ℕ) :
    dilog x - dilogPartial x N = ∑' i : ℕ, dilogTerm x (i + N) := by
  have hsplit := (summable_dilogTerm hx).sum_add_tsum_nat_add (f := dilogTerm x) N
  rw [dilog, dilogPartial]
  linear_combination -hsplit

/-- **Lower bound.**  For `0 ≤ x < 1` all terms are nonnegative, so every partial
sum is below the limit. -/
theorem dilogPartial_le_dilog {x : ℝ} (hx0 : 0 ≤ x) (hx1 : x < 1) (N : ℕ) :
    dilogPartial x N ≤ dilog x := by
  have hax : |x| < 1 := by rwa [abs_of_nonneg hx0]
  have htail : 0 ≤ dilog x - dilogPartial x N := by
    rw [dilog_sub_dilogPartial hax N]
    refine tsum_nonneg (fun i => dilogTerm_nonneg hx0 _)
  linarith

/-- **Upper bound.**  The tail of `∑ xⁿ/n²` beyond `N` terms is at most
`x^{N+1}/((N+1)²(1-x))`. -/
theorem dilog_le_dilogPartial_add {x : ℝ} (hx0 : 0 ≤ x) (hx1 : x < 1) (N : ℕ) :
    dilog x ≤ dilogPartial x N + x ^ (N + 1) / (((N : ℝ) + 1) ^ 2 * (1 - x)) := by
  have hax : |x| < 1 := by rwa [abs_of_nonneg hx0]
  have h1x : (0 : ℝ) < 1 - x := by linarith
  have hd : (0 : ℝ) < ((N : ℝ) + 1) ^ 2 := by positivity
  -- majorise the tail term by term
  have hmaj : Summable (fun i : ℕ => x ^ (i + N + 1) / (((N : ℝ) + 1) ^ 2)) := by
    have hg : Summable (fun i : ℕ => x ^ (i + N + 1)) := by
      have h := summable_geometric_of_lt_one hx0 hx1
      have : ∀ i : ℕ, x ^ (i + N + 1) = x ^ (N + 1) * x ^ i := by
        intro i; rw [← pow_add]; ring_nf
      exact (h.mul_left (x ^ (N + 1))).congr (fun i => (this i).symm)
    exact hg.div_const _
  have hshift : Summable (fun i : ℕ => dilogTerm x (i + N)) :=
    ((summable_dilogTerm hax).comp_injective (add_left_injective N))
  have hle : ∑' i : ℕ, dilogTerm x (i + N)
      ≤ ∑' i : ℕ, x ^ (i + N + 1) / (((N : ℝ) + 1) ^ 2) :=
    hshift.tsum_le_tsum (fun i => dilogTerm_le hx0 N i) hmaj
  have hval : ∑' i : ℕ, x ^ (i + N + 1) / (((N : ℝ) + 1) ^ 2)
      = x ^ (N + 1) / (((N : ℝ) + 1) ^ 2 * (1 - x)) := by
    have hrw : ∀ i : ℕ, x ^ (i + N + 1) / (((N : ℝ) + 1) ^ 2)
        = (x ^ (N + 1) / (((N : ℝ) + 1) ^ 2)) * x ^ i := by
      intro i
      rw [show i + N + 1 = (N + 1) + i from by omega, pow_add]
      ring
    rw [tsum_congr hrw, tsum_mul_left, tsum_geometric_of_lt_one hx0 hx1]
    field_simp
  have := dilog_sub_dilogPartial hax N
  linarith [hle.trans_eq hval, this]

/-- **The enclosure, packaged.**  For `0 ≤ x < 1`, `Li₂ x` lies between the
`N`-term partial sum and that sum plus an explicit tail. -/
theorem dilog_mem_Icc {x : ℝ} (hx0 : 0 ≤ x) (hx1 : x < 1) (N : ℕ) :
    dilog x ∈ Set.Icc (dilogPartial x N)
      (dilogPartial x N + x ^ (N + 1) / (((N : ℝ) + 1) ^ 2 * (1 - x))) :=
  ⟨dilogPartial_le_dilog hx0 hx1 N, dilog_le_dilogPartial_add hx0 hx1 N⟩

/-- The width of the enclosure — what one has to make small. -/
theorem dilog_enclosure_width {x : ℝ} (hx0 : 0 ≤ x) (hx1 : x < 1) (N : ℕ)
    {ε : ℝ} (hε : x ^ (N + 1) / (((N : ℝ) + 1) ^ 2 * (1 - x)) ≤ ε) :
    |dilog x - dilogPartial x N| ≤ ε := by
  rcases dilog_mem_Icc hx0 hx1 N with ⟨h1, h2⟩
  rw [abs_le]
  constructor <;> linarith

end Real
