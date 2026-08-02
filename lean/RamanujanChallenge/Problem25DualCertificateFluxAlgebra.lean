import RamanujanChallenge.Problem25DualCertificateAnalytic

noncomputable section

namespace RamanujanChallenge.P25

set_option maxRecDepth 10000
set_option maxHeartbeats 2000000

theorem dualCertFlux_algebra (a b r : ℕ)
    (x D F dF Dx C : ℝ) (hD : D ≠ 0) :
    (C *
          ((((a + 1 : ℕ) : ℝ) * x ^ a * (1 - x ^ 2) ^ (b + 1) +
                x ^ (a + 1) * ((b + 1 : ℕ) : ℝ) *
                  (1 - x ^ 2) ^ b * (-2 * x)) * F +
            x ^ (a + 1) * (1 - x ^ 2) ^ (b + 1) * dF) *
          D ^ (r + 1) -
        C * (x ^ (a + 1) * (1 - x ^ 2) ^ (b + 1) * F) *
          (((r + 1 : ℕ) : ℝ) * D ^ r * Dx)) /
        (D ^ (r + 1)) ^ 2 =
      (C * x ^ a * (1 - x ^ 2) ^ b / (D ^ (r + 1) * D)) *
        (D * (x * (1 - x ^ 2) * dF +
              (((a + 1 : ℕ) : ℝ) -
                (((a + 1 : ℕ) : ℝ) + 2 * ((b + 1 : ℕ) : ℝ)) * x ^ 2) * F) -
          ((r + 1 : ℕ) : ℝ) * x * (1 - x ^ 2) * Dx * F) := by
  have hxSucc : x ^ (a + 1) = x ^ a * x := by rw [pow_succ]
  have hOneSucc : (1 - x ^ 2) ^ (b + 1) =
      (1 - x ^ 2) ^ b * (1 - x ^ 2) := by rw [pow_succ]
  have hDSucc : D ^ (r + 1) = D ^ r * D := by rw [pow_succ]
  rw [hxSucc, hOneSucc, hDSucc]
  field_simp [hD]
  ring

end RamanujanChallenge.P25

end
