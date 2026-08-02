import RamanujanChallenge.Problem25

noncomputable section

namespace RamanujanChallenge.P25

/-- A column solution of the adjoint cocycle preserves a zero pairing with
the Catalan-error row, up to the irrelevant scalar multiplier `scale`. -/
theorem positiveCatalanError_pairing_zero_of_adjoint
    (W : ℕ → Fin 3 → ℝ) (scale : ℕ → ℝ)
    (hadjoint : ∀ (n : ℕ) (i : Fin 3),
      ∑ j : Fin 3, (positiveMatrix (n : ℤ) i j : ℝ) * W (n + 1) j =
        scale n * W n i)
    (hzero : ∑ j : Fin 3, positiveCatalanError 0 j * W 0 j = 0) :
    ∀ n, ∑ j : Fin 3, positiveCatalanError n j * W n j = 0 := by
  intro n
  induction n with
  | zero => exact hzero
  | succ n ih =>
      have h0 := hadjoint n (0 : Fin 3)
      have h1 := hadjoint n (1 : Fin 3)
      have h2 := hadjoint n (2 : Fin 3)
      simp only [Fin.sum_univ_three] at h0 h1 h2 ih ⊢
      rw [positiveCatalanError_succ, positiveCatalanError_succ,
        positiveCatalanError_succ]
      simp only [Fin.sum_univ_three]
      calc
        _ = positiveCatalanError n 0 *
              ((positiveMatrix (n : ℤ) 0 0 : ℝ) * W (n + 1) 0 +
                (positiveMatrix (n : ℤ) 0 1 : ℝ) * W (n + 1) 1 +
                (positiveMatrix (n : ℤ) 0 2 : ℝ) * W (n + 1) 2) +
            positiveCatalanError n 1 *
              ((positiveMatrix (n : ℤ) 1 0 : ℝ) * W (n + 1) 0 +
                (positiveMatrix (n : ℤ) 1 1 : ℝ) * W (n + 1) 1 +
                (positiveMatrix (n : ℤ) 1 2 : ℝ) * W (n + 1) 2) +
            positiveCatalanError n 2 *
              ((positiveMatrix (n : ℤ) 2 0 : ℝ) * W (n + 1) 0 +
                (positiveMatrix (n : ℤ) 2 1 : ℝ) * W (n + 1) 1 +
                (positiveMatrix (n : ℤ) 2 2 : ℝ) * W (n + 1) 2) := by ring
        _ = scale n *
            (positiveCatalanError n 0 * W n 0 +
              positiveCatalanError n 1 * W n 1 +
              positiveCatalanError n 2 * W n 2) := by
              rw [h0, h1, h2]
              ring
        _ = 0 := by rw [ih]; ring

/-- A positive covector whose pairing with an error row vanishes forces the
two endpoint errors to have opposite weak signs, provided the underlying
ratios are ordered.  The three hypotheses `h01`, `h02`, and `h12` are the
cross-multiplied form of that ordering. -/
theorem endpoint_bracket_of_positive_pair
    (E Q W : Fin 3 → ℝ)
    (hQ0 : 0 < Q 0) (hQ1 : 0 < Q 1) (hQ2 : 0 < Q 2)
    (hW0 : 0 < W 0) (hW1 : 0 < W 1) (hW2 : 0 < W 2)
    (hpair : E 0 * W 0 + E 1 * W 1 + E 2 * W 2 = 0)
    (h01 : 0 < E 1 * Q 0 - E 0 * Q 1)
    (h02 : 0 < E 2 * Q 0 - E 0 * Q 2)
    (h12 : 0 ≤ E 2 * Q 1 - E 1 * Q 2) :
    E 0 ≤ 0 ∧ 0 ≤ E 2 := by
  constructor
  · by_contra hE0
    have hE0' : 0 < E 0 := lt_of_not_ge hE0
    have hE0Q1 : 0 < E 0 * Q 1 := mul_pos hE0' hQ1
    have hE0Q2 : 0 < E 0 * Q 2 := mul_pos hE0' hQ2
    have hE1Q0 : 0 < E 1 * Q 0 := by linarith
    have hE2Q0 : 0 < E 2 * Q 0 := by linarith
    have hE1 : 0 < E 1 := pos_of_mul_pos_left hE1Q0 hQ0.le
    have hE2 : 0 < E 2 := pos_of_mul_pos_left hE2Q0 hQ0.le
    have h0 : 0 < E 0 * W 0 := mul_pos hE0' hW0
    have h1 : 0 < E 1 * W 1 := mul_pos hE1 hW1
    have h2 : 0 < E 2 * W 2 := mul_pos hE2 hW2
    linarith
  · by_contra hE2
    have hE2' : E 2 < 0 := lt_of_not_ge hE2
    have hE2Q0 : E 2 * Q 0 < 0 := mul_neg_of_neg_of_pos hE2' hQ0
    have hE2Q1 : E 2 * Q 1 < 0 := mul_neg_of_neg_of_pos hE2' hQ1
    have hE0Q2 : E 0 * Q 2 < 0 := by linarith
    have hE1Q2 : E 1 * Q 2 < 0 := by linarith
    have hE0 : E 0 < 0 := neg_of_mul_neg_left hE0Q2 hQ2.le
    have hE1 : E 1 < 0 := neg_of_mul_neg_left hE1Q2 hQ2.le
    have h0 : E 0 * W 0 < 0 := mul_neg_of_neg_of_pos hE0 hW0
    have h1 : E 1 * W 1 < 0 := mul_neg_of_neg_of_pos hE1 hW1
    have h2 : E 2 * W 2 < 0 := mul_neg_of_neg_of_pos hE2' hW2
    linarith

end RamanujanChallenge.P25

end
