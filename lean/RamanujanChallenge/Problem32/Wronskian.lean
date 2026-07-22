/-
  Problem 3.2 — Layer 2: The Wronskian Identity.

  Lemma (Wronskian): For all n ≥ 1,
    W_n · n³ = 6, equivalently W_n = 6/n³.

  Proof: Since both a_n and b_n satisfy the same three-term recurrence,
  the Wronskian satisfies (n+1)³ W_{n+1} = n³ W_n.
  Telescoping from W_1 = 6 gives W_n · n³ = 6.
-/
import RamanujanChallenge.Problem32.AperyDef

noncomputable section

/-! ## Base case: W_1 = 6 -/

theorem wronskian_one : wronskian 1 = 6 := by
  simp [wronskian, aperyA, aperyBQ]

/-! ## Wronskian ratio step

(n+1)³ W_{n+1} = n³ W_n, proved from the shared recurrence.
-/

theorem wronskian_step (n : ℕ) (hn : n ≥ 1) :
    ((n : ℚ) + 1) ^ 3 * wronskian (n + 1) = (n : ℚ) ^ 3 * wronskian n := by
  have hA := aperyA_recurrence n hn
  have hB := aperyBQ_recurrence n hn
  simp only [wronskian]
  cases n with
  | zero => omega
  | succ m =>
    simp only [Nat.succ_sub_one] at hA hB ⊢
    linear_combination aperyBQ (m + 1) * hA - aperyA (m + 1) * hB

/-! ## Auxiliary: the step in fully cast-normalized form -/

private theorem wronskian_step' (k : ℕ) :
    ((k : ℚ) + 2) ^ 3 * wronskian (k + 2) =
    ((k : ℚ) + 1) ^ 3 * wronskian (k + 1) := by
  have h := wronskian_step (k + 1) (by omega)
  push_cast at h ⊢
  linarith

/-! ## The Wronskian identity (multiplied form, no division)

W_n · n³ = 6 for all n ≥ 1.
-/

theorem wronskian_mul (n : ℕ) (hn : n ≥ 1) :
    wronskian n * (n : ℚ) ^ 3 = 6 := by
  induction n with
  | zero => omega
  | succ m ih =>
    cases m with
    | zero =>
      simp [wronskian_one]
    | succ k =>
      have hk1 : k + 1 ≥ 1 := by omega
      have ihm := ih hk1
      have hstep := wronskian_step' k
      push_cast at ihm ⊢
      nlinarith [mul_comm (wronskian (k + 1)) (((k : ℚ) + 1) ^ 3),
                 mul_comm (wronskian (k + 2)) (((k : ℚ) + 2) ^ 3)]

/-! ## The Wronskian identity (division form) -/

theorem wronskian_eq (n : ℕ) (hn : n ≥ 1) :
    wronskian n = 6 / (n : ℚ) ^ 3 := by
  have hne : (n : ℚ) ^ 3 ≠ 0 := by positivity
  have hmul := wronskian_mul n hn
  field_simp at hmul ⊢
  linarith

/-! ## b_0 = 1 (for the no-consecutive-zeros proof) -/

theorem aperyBQ_zero_eq : aperyBQ 0 = 1 := aperyBQ_zero

end
