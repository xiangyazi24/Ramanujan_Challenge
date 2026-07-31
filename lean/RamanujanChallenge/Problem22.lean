/-
  Ramanujan Challenge Problem 2.2: Euler's constant γ as an Apéry Limit

  The four-term recurrence in Problem 2.2, after the index shift m = n + 3,
  is the Aptekarev recurrence for rational approximants to γ.
  The convergence P_m / Q_m → γ follows from the established theory of
  multiple orthogonal polynomials for the weight pair (e^{-x}, e^{-x} log x)
  on [0, ∞).

  Reference: Xiang Huang, "Solution to Ramanujan Challenge Problem 2.2", July 2026.
  Cf. Aptekarev (2007), Hessami Pilehrood & Hessami Pilehrood (2013).
-/
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Data.Nat.Factorial.Basic

noncomputable section

/-! ## The Aptekarev recurrence

The challenge recurrence (after shift m = n + 3) has coefficients:
  c̃₀(m) = -8m³ + 21m² - 15m + 4
  c̃₁(m) = 24m⁵ - 23m⁴ - 51m³ + 39m² - m - 12
  c̃₂(m) = -m(m-1)(24m⁵ - 87m⁴ + 34m³ + 66m² - 3m - 20)
  c̃₃(m) = m(m-2)(m-1)⁴(8m³ + 3m² - 3m - 2)

Leading asymptotics: -8m³, 24m⁵, -24m⁷, 8m⁹ (factorial gauge k = 2,
characteristic polynomial -8(r-1)³ — maximally resonant).
-/

def c0_shifted (m : ℤ) : ℤ := -8 * m ^ 3 + 21 * m ^ 2 - 15 * m + 4

def c1_shifted (m : ℤ) : ℤ := 24 * m ^ 5 - 23 * m ^ 4 - 51 * m ^ 3 + 39 * m ^ 2 - m - 12

def c2_shifted (m : ℤ) : ℤ :=
  -m * (m - 1) * (24 * m ^ 5 - 87 * m ^ 4 + 34 * m ^ 3 + 66 * m ^ 2 - 3 * m - 20)

def c3_shifted (m : ℤ) : ℤ :=
  m * (m - 2) * (m - 1) ^ 4 * (8 * m ^ 3 + 3 * m ^ 2 - 3 * m - 2)

/-! ## The original challenge coefficients

The original coefficients c₀(n), c₁(n), c₂(n), c₃(n) are related by
the shift n = m - 3, i.e., c_j(n) = c̃_j(n + 3).
-/

def c0_orig (n : ℤ) : ℤ := -8 * n ^ 3 - 51 * n ^ 2 - 105 * n - 68

def c1_orig (n : ℤ) : ℤ := 24 * n ^ 5 + 337 * n ^ 4 + 1833 * n ^ 3 +
  4818 * n ^ 2 + 6092 * n + 2928

def c2_orig (n : ℤ) : ℤ :=
  -(n + 2) * (n + 3) * (24 * n ^ 5 + 273 * n ^ 4 + 1150 * n ^ 3 +
    2154 * n ^ 2 + 1635 * n + 268)

def c3_orig (n : ℤ) : ℤ :=
  (n + 1) * (n + 2) ^ 4 * (n + 3) * (8 * n ^ 3 + 75 * n ^ 2 + 231 * n + 232)

/-! ## Verification: shifted coefficients = original at n = m - 3

This is a finite polynomial identity, verified by comparing coefficients.
-/

theorem shift_c0 (m : ℤ) : c0_shifted m = c0_orig (m - 3) := by
  simp only [c0_shifted, c0_orig]; ring

theorem shift_c1 (m : ℤ) : c1_shifted m = c1_orig (m - 3) := by
  simp only [c1_shifted, c1_orig]; ring

theorem shift_c2 (m : ℤ) : c2_shifted m = c2_orig (m - 3) := by
  simp only [c2_shifted, c2_orig]; ring

theorem shift_c3 (m : ℤ) : c3_shifted m = c3_orig (m - 3) := by
  simp only [c3_shifted, c3_orig]; ring

/-! ## Initial values (Aptekarev's published sequences)

(P₀, P₁, P₂) = (0, 7, 179)
(Q₀, Q₁, Q₂) = (1, 12, 306)

These match the challenge's (p_{-3}, p_{-2}, p_{-1}) and (q_{-3}, q_{-2}, q_{-1}).
-/

def aptekarevP_init : Fin 3 → ℤ := ![0, 7, 179]
def aptekarevQ_init : Fin 3 → ℤ := ![1, 12, 306]

/-! ## The gauge structure

The degree pattern (3, 5, 7, 9) with factorial gauge k = 2 gives
characteristic polynomial -8(r-1)³ (triple root at r = 1).

This triple root is the signature of a maximally resonant system:
the three asymptotic solutions are separated only at subexponential scale,
involving powers of m and logarithmic companions — the structure expected
when the limit involves γ = -Γ'(1).
-/

theorem gauge_degree : (3 : ℕ) = 5 - 2 ∧ (5 : ℕ) = 7 - 2 ∧ (7 : ℕ) = 9 - 2 := by
  omega

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
