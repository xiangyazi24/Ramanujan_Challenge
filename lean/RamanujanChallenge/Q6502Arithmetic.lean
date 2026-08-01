import RamanujanChallenge.Q6502Transfer

open Filter Topology Real
open scoped BigOperators

noncomputable section

namespace RamanujanChallenge.P22.Q6502

/-- On the good saddle region, the left reciprocal is controlled by the
first moment rather than by an explicit square root. -/
theorem inv_add_one_le_of_saddle_good
    {N K J G : ℝ}
    (hN : 0 < N) (hK0 : 0 ≤ K) (hJ0 : 0 ≤ J)
    (hJ : J = N - K)
    (hGdef : G = K ^ 3 - J ^ 2)
    (hG : |G| ≤ N ^ 2 / 8)
    (hK : K ≤ N / 2) :
    1 / (K + 1) ≤ 8 * (K + 1) / N := by
  have hJhalf : N / 2 ≤ J := by rw [hJ]; linarith
  have hNhalf0 : 0 ≤ N / 2 := by linarith
  have hJsq : N ^ 2 / 4 ≤ J ^ 2 := by
    have hprod := mul_nonneg (sub_nonneg.mpr hJhalf)
      (add_nonneg hJ0 hNhalf0)
    nlinarith
  have hGlow : -(N ^ 2 / 8) ≤ G := neg_le_of_abs_le hG
  have hKcube : N ^ 2 / 8 ≤ K ^ 3 := by
    rw [hGdef] at hGlow
    nlinarith
  have hKleN : K ≤ N := by linarith
  have hmul := mul_le_mul_of_nonneg_right hKleN (sq_nonneg K)
  have hKsq : N / 8 ≤ K ^ 2 := by
    nlinarith
  have hNle : N ≤ 8 * (K + 1) ^ 2 := by
    nlinarith [sq_nonneg K]
  have hK1 : 0 < K + 1 := by linarith
  field_simp [hN.ne', hK1.ne']
  nlinarith

/-- On the good saddle region, the right reciprocal is uniformly `O(1/n)`. -/
theorem inv_right_add_one_le_of_saddle_good
    {N K J : ℝ}
    (hN : 0 < N) (hJ0 : 0 ≤ J)
    (hJ : J = N - K) (hK : K ≤ N / 2) :
    1 / (J + 1) ≤ 2 / N := by
  have hJhalf : N / 2 ≤ J := by rw [hJ]; linarith
  have hJ1 : 0 < J + 1 := by linarith
  field_simp [hN.ne', hJ1.ne']
  nlinarith

/-- The polynomial badness majorant is at least one off the good region. -/
theorem one_le_saddle_badness
    {N K G : ℝ} (hN : 0 < N)
    (hbad : N ^ 2 / 8 < |G| ∨ N / 2 < K) :
    1 ≤ 64 * G ^ 2 / N ^ 4 + 4 * K ^ 2 / N ^ 2 := by
  rcases hbad with hG | hK
  · have hsq : N ^ 4 / 64 < G ^ 2 := by
      have habs0 : 0 ≤ |G| := abs_nonneg G
      have hN20 : 0 ≤ N ^ 2 / 8 := by positivity
      nlinarith [sq_abs G]
    have hN4 : 0 < N ^ 4 := by positivity
    have hfirst : 1 < 64 * G ^ 2 / N ^ 4 := by
      apply (lt_div_iff₀ hN4).2
      nlinarith
    have hsecond : 0 ≤ 4 * K ^ 2 / N ^ 2 := by positivity
    linarith
  · have hsq : N ^ 2 / 4 < K ^ 2 := by
      have hNhalf0 : 0 ≤ N / 2 := by linarith
      nlinarith [sq_nonneg (K - N / 2)]
    have hN2 : 0 < N ^ 2 := by positivity
    have hsecond : 1 < 4 * K ^ 2 / N ^ 2 := by
      apply (lt_div_iff₀ hN2).2
      nlinarith
    have hfirst : 0 ≤ 64 * G ^ 2 / N ^ 4 := by positivity
    linarith

#print axioms inv_add_one_le_of_saddle_good
#print axioms inv_right_add_one_le_of_saddle_good
#print axioms one_le_saddle_badness

end RamanujanChallenge.P22.Q6502
