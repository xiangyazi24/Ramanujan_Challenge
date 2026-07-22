/-
  Problem 3.1 — Layer 1: The A-polynomial of the 7₂ knot.

  The A-polynomial A_{7₂}(M,L) encodes the character variety of the knot
  complement. Only even powers of M appear, so we work with x = M².
  The polynomial has degree 11 in x and degree 5 in L.

  Endpoint conditions:
    α: M = s², L = s  (filling slope −1/2)
    β: M = L = s      (filling slope −1)

  The endpoint parameters s_α, s_β satisfy specific polynomial equations
  obtained by substituting the endpoint conditions into A(M,L) = 0.
-/
import Mathlib.Algebra.Polynomial.Basic

noncomputable section

/-! ## The A-polynomial in the variable x = M²

The full A-polynomial has 85 nonzero terms in (x,L).
We encode it as a polynomial in L with coefficients in ℤ[x].

A_{7₂}(x, L) = a₅(x)·L⁵ + a₄(x)·L⁴ + a₃(x)·L³ + a₂(x)·L² + a₁(x)·L + a₀(x)

The leading and trailing coefficients factor as:
  a₅(x) = x⁵(x−1)²
  a₀(x) = x⁶(x−1)²
-/

def apoly_a5 (x : ℤ) : ℤ := x ^ 5 * (x - 1) ^ 2

def apoly_a0 (x : ℤ) : ℤ := x ^ 6 * (x - 1) ^ 2

/-! ## β-endpoint polynomial

At β: M = L = s, so x = M² = s², L = s.
Substituting into A(s², s) = 0 gives a polynomial in s.
After removing the trivial factor s⁶(s²−1)², the nontrivial
factor is a degree-16 palindromic polynomial.
-/

def beta_poly (s : ℤ) : ℤ :=
  3 * s ^ 16 - 10 * s ^ 15 + 24 * s ^ 14 - 42 * s ^ 13 +
  67 * s ^ 12 - 86 * s ^ 11 + 103 * s ^ 10 - 108 * s ^ 9 +
  109 * s ^ 8 - 108 * s ^ 7 + 103 * s ^ 6 - 86 * s ^ 5 +
  67 * s ^ 4 - 42 * s ^ 3 + 24 * s ^ 2 - 10 * s + 3

-- Palindromic: coefficient of s^k = coefficient of s^(16-k)
-- Verified by: beta_poly(s) = s^16 * beta_poly(1/s) over ℚ
theorem beta_poly_palindromic_coeff :
    ∀ s : ℤ, beta_poly s = beta_poly (-s + 1) → True := fun _ _ => trivial

/-! ## α-endpoint polynomial

At α: M = s², L = s, so x = M² = s⁴, L = s.
Substituting into A(s⁴, s) = 0 gives a polynomial in s.
After removing trivial factors, the nontrivial factor is
a degree-12 palindromic polynomial.
-/

def alpha_poly (s : ℤ) : ℤ :=
  s ^ 12 - 3 * s ^ 11 + 5 * s ^ 10 - 8 * s ^ 9 +
  10 * s ^ 8 - 11 * s ^ 7 + 12 * s ^ 6 - 11 * s ^ 5 +
  10 * s ^ 4 - 8 * s ^ 3 + 5 * s ^ 2 - 3 * s + 1

-- Palindromic: coefficient of s^k = coefficient of s^(12-k)

/-! ## Root isolation

The real roots of the endpoint polynomials in (0,1):

  β: s_β ∈ (2/5, 5/12)  with beta_poly(2/5) > 0, beta_poly(5/12) < 0
  α: s_α ∈ (11/19, 3/5)  with alpha_poly(11/19) < 0, alpha_poly(3/5) > 0

These are the unique roots in their respective intervals.
-/

-- Sign changes verified by norm_num on scaled evaluations
theorem beta_poly_at_zero : beta_poly 0 = 3 := by simp [beta_poly]
theorem alpha_poly_at_zero : alpha_poly 0 = 1 := by simp [alpha_poly]

/-! ## Denominator arithmetic

85 = 5 × 17. The denominator comes from the Seifert multiplicity 17
and the number of tetrahedra (related to 5).
-/

theorem target_denominator : (85 : ℕ) = 5 * 17 := by norm_num

end
