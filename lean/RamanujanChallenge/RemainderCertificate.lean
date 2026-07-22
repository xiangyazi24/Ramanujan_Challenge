/-
  RemainderCertificate: Shared infrastructure for the remainder-certificate pattern.

  The universal strategy for all recurrence-based Ramanujan Challenge problems:
  Given sequences (p_n), (q_n) satisfying a linear recurrence and a constant C,
  define the remainder r_n = q_n · C − p_n. If r_n satisfies the SAME recurrence
  and r_n → 0, then p_n / q_n → C.

  This avoids all asymptotic Poincaré-Perron machinery.
-/
import Mathlib.Topology.Algebra.InfiniteSum.Basic
import Mathlib.Topology.Order.Basic
import Mathlib.Analysis.Normed.Group.Basic
import Mathlib.Order.Filter.Basic
import Mathlib.Analysis.SpecificLimits.Basic

open Filter

noncomputable section

/-! ## Core convergence lemma

If |q_n| ≥ 1 eventually and q_n · C − p_n → 0, then p_n / q_n → C.

The bound |q_n| ≥ 1 is trivially satisfied for integer-valued sequences
that are eventually nonzero.
-/

theorem tendsto_div_of_remainder_tendsto_zero
    {p q : ℕ → ℝ} {C : ℝ}
    (hq : ∀ᶠ n in atTop, (1 : ℝ) ≤ |q n|)
    (hr : Tendsto (fun n => q n * C - p n) atTop (nhds 0)) :
    Tendsto (fun n => p n / q n) atTop (nhds C) := by
  have hqne : ∀ᶠ n in atTop, q n ≠ 0 := by
    filter_upwards [hq] with n hn
    exact abs_pos.mp (lt_of_lt_of_le zero_lt_one hn)
  have key : ∀ᶠ n in atTop, p n / q n = C - (q n * C - p n) / q n := by
    filter_upwards [hqne] with n hqn
    field_simp
    ring
  rw [show (nhds C : Filter ℝ) = nhds (C - 0) from by simp]
  refine (Tendsto.sub tendsto_const_nhds ?_).congr'
    (key.mono fun _ hn => hn.symm)
  rw [Metric.tendsto_atTop] at hr ⊢
  intro ε hε
  obtain ⟨N, hN⟩ := hr ε hε
  obtain ⟨M, hM⟩ := hq.exists_forall_of_atTop
  refine ⟨max N M, fun n hn => ?_⟩
  simp only [Real.dist_eq, sub_zero] at hN ⊢
  have hn1 : N ≤ n := le_of_max_le_left hn
  have hn2 : (1 : ℝ) ≤ |q n| := hM n (le_of_max_le_right hn)
  calc |(q n * C - p n) / q n|
      = |q n * C - p n| / |q n| := abs_div _ _
    _ ≤ |q n * C - p n| / 1 := by
        exact div_le_div_of_nonneg_left (abs_nonneg _) one_pos hn2
    _ = |q n * C - p n| := div_one _
    _ < ε := hN n hn1

/-! ## Geometric decay lemma

If |r_{n+1}| ≤ ρ · |r_n| for all n ≥ N, with 0 ≤ ρ < 1, then r_n → 0.
-/

theorem geometric_bound_induction
    {r : ℕ → ℝ} {ρ : ℝ} {N : ℕ}
    (hρ : 0 ≤ ρ)
    (hbound : ∀ n, N ≤ n → |r (n + 1)| ≤ ρ * |r n|) :
    ∀ k, |r (N + k)| ≤ |r N| * ρ ^ k := by
  intro k
  induction k with
  | zero => simp
  | succ m ih =>
    have : N + (m + 1) = (N + m) + 1 := by ring
    calc |r (N + (m + 1))| = |r ((N + m) + 1)| := by rw [this]
      _ ≤ ρ * |r (N + m)| := hbound (N + m) (Nat.le_add_right N m)
      _ ≤ ρ * (|r N| * ρ ^ m) := by
          exact mul_le_mul_of_nonneg_left ih hρ
      _ = |r N| * ρ ^ (m + 1) := by ring

theorem tendsto_zero_of_geometric_bound
    {r : ℕ → ℝ} {ρ : ℝ} {N : ℕ}
    (hρ : 0 ≤ ρ) (hρ1 : ρ < 1)
    (hbound : ∀ n, N ≤ n → |r (n + 1)| ≤ ρ * |r n|) :
    Tendsto r atTop (nhds 0) := by
  rw [tendsto_zero_iff_abs_tendsto_zero]
  rw [← tendsto_add_atTop_iff_nat N]
  have hgeom : Tendsto (fun k : ℕ => |r N| * ρ ^ k) atTop (nhds 0) := by
    simpa using
      (tendsto_pow_atTop_nhds_zero_of_lt_one hρ hρ1).const_mul |r N|
  exact tendsto_of_tendsto_of_tendsto_of_le_of_le
    tendsto_const_nhds hgeom
    (fun _ => abs_nonneg _)
    (fun k => by
      simpa [Nat.add_comm] using geometric_bound_induction hρ hbound k)

/-! ## Recurrence-linearity lemma

If p and q both satisfy the same linear recurrence, then
r_n = q_n · C − p_n also satisfies that recurrence (for any constant C).
This is the KEY property: the remainder inherits the recurrence.
-/

theorem remainder_satisfies_recurrence_two_term
    {a b p q : ℕ → ℝ} {C : ℝ}
    (hp : ∀ n, p (n + 2) = a n * p (n + 1) + b n * p n)
    (hq : ∀ n, q (n + 2) = a n * q (n + 1) + b n * q n) :
    ∀ n, (q (n + 2) * C - p (n + 2)) =
      a n * (q (n + 1) * C - p (n + 1)) + b n * (q n * C - p n) := by
  intro n
  rw [hp n, hq n]
  ring

theorem remainder_satisfies_recurrence_three_term
    {a b c p q : ℕ → ℝ} {C : ℝ}
    (hp : ∀ n, p (n + 3) = a n * p (n + 2) + b n * p (n + 1) + c n * p n)
    (hq : ∀ n, q (n + 3) = a n * q (n + 2) + b n * q (n + 1) + c n * q n) :
    ∀ n, (q (n + 3) * C - p (n + 3)) =
      a n * (q (n + 2) * C - p (n + 2)) +
      b n * (q (n + 1) * C - p (n + 1)) +
      c n * (q n * C - p n) := by
  intro n
  rw [hp n, hq n]
  ring

theorem remainder_satisfies_recurrence_four_term
    {a b c d p q : ℕ → ℝ} {C : ℝ}
    (hp : ∀ n, p (n + 4) = a n * p (n + 3) + b n * p (n + 2) +
      c n * p (n + 1) + d n * p n)
    (hq : ∀ n, q (n + 4) = a n * q (n + 3) + b n * q (n + 2) +
      c n * q (n + 1) + d n * q n) :
    ∀ n, (q (n + 4) * C - p (n + 4)) =
      a n * (q (n + 3) * C - p (n + 3)) +
      b n * (q (n + 2) * C - p (n + 2)) +
      c n * (q (n + 1) * C - p (n + 1)) +
      d n * (q n * C - p n) := by
  intro n
  rw [hp n, hq n]
  ring

end
