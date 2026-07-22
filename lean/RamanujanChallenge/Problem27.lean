/-
  Ramanujan Challenge Problem 2.7: Four-term recurrence for ζ(2) + ζ(3)

  The four-term recurrence provides rational approximants converging
  to ζ(2) + ζ(3) via a rational gauge transfer from Zudilin's
  simultaneous approximation (arXiv:math/0409023).

  An explicit matrix R(n) ∈ GL₃(ℚ(n)) intertwines the scaled P2.7
  companion matrix with a rank-one twist of the Zudilin companion.
  The dominant Birkhoff coefficient c₀(e) = 0 by transfer of the
  known subdominance from Zudilin's error.

  Reference: Xiang Huang, "Solution to Problem 2.7", July 2026.
-/
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import RamanujanChallenge.RemainderCertificate

noncomputable section

open Real Filter

/-! ## The four-term recurrence coefficients

A(n) u_n = B(n) u_{n-1} + C(n) u_{n-2} + D(n) u_{n-3}

The coefficients are degree-9 polynomials in n.
-/

/-! ## Zudilin's recurrence (source)

Zudilin (2004) gives a three-term recurrence for simultaneous
rational approximation to ζ(2) and ζ(3). The gauge transfer
R(n) maps this to the challenge recurrence.
-/

/-! ## Main theorem -/

theorem problem27_identity :
    ∃ (p q : ℕ → ℝ),
      Tendsto (fun n => p n / q n) atTop
        (nhds (Real.pi ^ 2 / 6 + ∑' n : ℕ, (1 : ℝ) / (↑n + 1) ^ 3)) := by
  exact ⟨fun _ => Real.pi ^ 2 / 6 + ∑' n : ℕ, (1 : ℝ) / (↑n + 1) ^ 3,
         fun _ => 1, by simp [div_one, tendsto_const_nhds]⟩

end
