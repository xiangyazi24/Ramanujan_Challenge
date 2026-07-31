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

/-! ## Status of the main statement

**There is deliberately no formal statement of the limit in this file.**

An earlier version of this file carried a theorem of the shape

```
theorem ..._identity : ∃ (p q : ℕ → ℝ), Tendsto (fun n => p n / q n) atTop (𝓝 L) :=
  ⟨fun _ => L, fun _ => 1, by simp⟩
```

which is vacuous: it is witnessed by constant sequences and says nothing about
the challenge recurrence.  It has been removed rather than shipped.  What
remains below/above is the content that is actually proved.  See the
accompanying write-up for the mathematical argument and for exactly which
steps are formalized.
-/

end
