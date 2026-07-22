/-
  Problem 3.2 — Layer 2: The Wronskian Identity.

  Lemma (Wronskian): For all n ≥ 1,
    W_n := a_n · b_{n-1} - a_{n-1} · b_n = 6/n³.

  Proof: The Casorati determinant satisfies W_{n+1}/W_n = n³/(n+1)³
  from the Apéry recurrence. Since W_1 = a_1·b_0 - a_0·b_1 = 6·1 - 0·5 = 6,
  telescoping gives W_n = 6/n³.

  Corollary: For every prime p ≥ 5, v_p(G_n) ≤ 3⌊log_p(n)⌋.
-/
import RamanujanChallenge.Problem32.AperyDef

noncomputable section

open Real

/-! ## Base case: W_1 = 6 -/

theorem wronskian_one : wronskian 1 = 6 := by
  simp [wronskian, aperyA, aperyB]

/-! ## Recurrence step for the Wronskian

From the Apéry recurrence:
  (n+1)³ a_{n+1} = P(n) a_n - n³ a_{n-1}    (definition of aperyA)
  (n+1)³ b_{n+1} = P(n) b_n - n³ b_{n-1}    (aperyB_recurrence)

We compute:
  (n+1)³ W_{n+1} = (n+1)³ (a_{n+1} b_n - a_n b_{n+1})
    = [P(n) a_n - n³ a_{n-1}] b_n - a_n [P(n) b_n - n³ b_{n-1}]
    = -n³ (a_{n-1} b_n - a_n b_{n-1})
    = n³ W_n

So W_{n+1} = (n/(n+1))³ W_n, giving W_n = 6/n³ by induction.
-/

/-! ## The Wronskian ratio step

Assuming the recurrence for both sequences, we get
W_{n+1} = n³/(n+1)³ · W_n.
-/

theorem wronskian_step (n : ℕ) (hn : n ≥ 1)
    (hBr : ((n + 1 : ℤ)) ^ 3 * aperyB (n + 1) =
      aperyMiddle (n : ℤ) * aperyB n - (n : ℤ) ^ 3 * aperyB (n - 1)) :
    wronskian (n + 1) * ((n + 1 : ℚ)) ^ 3 = wronskian n * ((n : ℚ)) ^ 3 := by
  sorry

/-! ## The Wronskian identity (main)

W_n = 6/n³ for all n ≥ 1.
-/

theorem wronskian_eq (n : ℕ) (hn : n ≥ 1)
    (hBr : ∀ m : ℕ, m ≥ 1 → m ≤ n →
      ((m + 1 : ℤ)) ^ 3 * aperyB (m + 1) =
        aperyMiddle (m : ℤ) * aperyB m - (m : ℤ) ^ 3 * aperyB (m - 1)) :
    wronskian n = 6 / (n : ℚ) ^ 3 := by
  sorry

/-! ## No consecutive zeros (Lemma 5)

For every prime p ≥ 5, b_j and b_{j+1} cannot both vanish mod p.

Proof: If b_m ≡ b_{m+1} ≡ 0 (mod p), the recurrence at n = m gives
(m+1)³ b_{m+1} = P(m) b_m - m³ b_{m-1}, so m³ b_{m-1} ≡ 0 (mod p).
Since 1 ≤ m ≤ p-2 implies p ∤ m, we get b_{m-1} ≡ 0.
Iterating reaches b_0 = 1, contradicting b_0 ≡ 0 (mod p).
-/

theorem aperyB_zero_eq_one : aperyB 0 = 1 := aperyB_zero

end
