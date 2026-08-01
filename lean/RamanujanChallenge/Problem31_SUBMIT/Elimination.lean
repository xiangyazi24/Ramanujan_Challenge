import Ramanujan31.APolynomial

/-!
# Exact elimination from the shape chart to the official A-polynomial

A computer algebra system was used only to discover the factorization below.
Lean checks it from the displayed polynomials.

Write

```
t = 1-r²,  D = t-rX,  N = r-Xt.
```

The holonomy equation is `L D = X² N`.  After this substitution, the official
A-polynomial satisfies

```
D⁵ A(X,L)
  = X⁶(X-1)⁵(X+1)⁵ forwardFactor(X,r) E(X,r).
```

Consequently every solution of the shape-chart holonomy and edge equations is
on the exact A-polynomial printed in the problem.  There is no appeal to a
computed resultant in the theorem statement or proof.
-/

/-- The extra factor in the denominator-free elimination identity. -/
def forwardEliminationFactor (X r : ℝ) : ℝ :=
  X ^ 2 * r ^ 4 + X ^ 2 * r ^ 3 - 2 * X ^ 2 * r ^ 2
    - X ^ 2 * r + X ^ 2
    + X * r ^ 5 + X * r ^ 4 - 2 * X * r ^ 3 - X * r ^ 2
    + X * r + X
    + r ^ 4 + r ^ 3 - 2 * r ^ 2 - r + 1

/-- The denominator-free substitution identity behind forward elimination. -/
theorem aPolynomialX_substitution
    (X L r : ℝ)
    (hL : L * (chartT r - r * X) = X ^ 2 * (r - X * chartT r)) :
    (chartT r - r * X) ^ 5 * aPolynomialX X L
      =
        X ^ 6 * (X - 1) ^ 5 * (X + 1) ^ 5
          * forwardEliminationFactor X r * edgeEquation X r := by
  have h2 :
      L ^ 2 * (chartT r - r * X) ^ 2
        = X ^ 4 * (r - X * chartT r) ^ 2 := by
    calc
      L ^ 2 * (chartT r - r * X) ^ 2
          = (L * (chartT r - r * X)) ^ 2 := by ring
      _ = (X ^ 2 * (r - X * chartT r)) ^ 2 := by rw [hL]
      _ = X ^ 4 * (r - X * chartT r) ^ 2 := by ring
  have h3 :
      L ^ 3 * (chartT r - r * X) ^ 3
        = X ^ 6 * (r - X * chartT r) ^ 3 := by
    calc
      L ^ 3 * (chartT r - r * X) ^ 3
          = (L * (chartT r - r * X)) ^ 3 := by ring
      _ = (X ^ 2 * (r - X * chartT r)) ^ 3 := by rw [hL]
      _ = X ^ 6 * (r - X * chartT r) ^ 3 := by ring
  have h4 :
      L ^ 4 * (chartT r - r * X) ^ 4
        = X ^ 8 * (r - X * chartT r) ^ 4 := by
    calc
      L ^ 4 * (chartT r - r * X) ^ 4
          = (L * (chartT r - r * X)) ^ 4 := by ring
      _ = (X ^ 2 * (r - X * chartT r)) ^ 4 := by rw [hL]
      _ = X ^ 8 * (r - X * chartT r) ^ 4 := by ring
  have h5 :
      L ^ 5 * (chartT r - r * X) ^ 5
        = X ^ 10 * (r - X * chartT r) ^ 5 := by
    calc
      L ^ 5 * (chartT r - r * X) ^ 5
          = (L * (chartT r - r * X)) ^ 5 := by ring
      _ = (X ^ 2 * (r - X * chartT r)) ^ 5 := by rw [hL]
      _ = X ^ 10 * (r - X * chartT r) ^ 5 := by ring
  calc
    (chartT r - r * X) ^ 5 * aPolynomialX X L
        =
          L ^ 5 * (chartT r - r * X) ^ 5
            + (X ^ 7 - X ^ 6 + 3 * X ^ 2 + 4 * X - 2)
                * (chartT r - r * X)
                * (L ^ 4 * (chartT r - r * X) ^ 4)
            + (-2 * X ^ 9 + 5 * X ^ 8 + X ^ 7 - 4 * X ^ 6
                + 6 * X ^ 4 + 5 * X ^ 3 + 2 * X ^ 2 - 4 * X + 1)
                * (chartT r - r * X) ^ 2
                * (L ^ 3 * (chartT r - r * X) ^ 3)
            + (X ^ 11 - 4 * X ^ 10 + 2 * X ^ 9 + 5 * X ^ 8 + 6 * X ^ 7
                - 4 * X ^ 5 + X ^ 4 + 5 * X ^ 3 - 2 * X ^ 2)
                * (chartT r - r * X) ^ 3
                * (L ^ 2 * (chartT r - r * X) ^ 2)
            + (-2 * X ^ 11 + 4 * X ^ 10 + 3 * X ^ 9 - X ^ 5 + X ^ 4)
                * (chartT r - r * X) ^ 4
                * (L * (chartT r - r * X))
            + X ^ 11 * (chartT r - r * X) ^ 5 := by
              unfold aPolynomialX
              ring
    _ =
          X ^ 10 * (r - X * chartT r) ^ 5
            + (X ^ 7 - X ^ 6 + 3 * X ^ 2 + 4 * X - 2)
                * (chartT r - r * X)
                * (X ^ 8 * (r - X * chartT r) ^ 4)
            + (-2 * X ^ 9 + 5 * X ^ 8 + X ^ 7 - 4 * X ^ 6
                + 6 * X ^ 4 + 5 * X ^ 3 + 2 * X ^ 2 - 4 * X + 1)
                * (chartT r - r * X) ^ 2
                * (X ^ 6 * (r - X * chartT r) ^ 3)
            + (X ^ 11 - 4 * X ^ 10 + 2 * X ^ 9 + 5 * X ^ 8 + 6 * X ^ 7
                - 4 * X ^ 5 + X ^ 4 + 5 * X ^ 3 - 2 * X ^ 2)
                * (chartT r - r * X) ^ 3
                * (X ^ 4 * (r - X * chartT r) ^ 2)
            + (-2 * X ^ 11 + 4 * X ^ 10 + 3 * X ^ 9 - X ^ 5 + X ^ 4)
                * (chartT r - r * X) ^ 4
                * (X ^ 2 * (r - X * chartT r))
            + X ^ 11 * (chartT r - r * X) ^ 5 := by
              rw [h5, h4, h3, h2, hL]
    _ =
        X ^ 6 * (X - 1) ^ 5 * (X + 1) ^ 5
          * forwardEliminationFactor X r * edgeEquation X r := by
            unfold forwardEliminationFactor edgeEquation chartT
            ring

/-- **Forward chart elimination.**  The two exact chart equations imply the
official A-polynomial, without any nonvanishing hypotheses. -/
theorem aPolynomialX_eq_zero_of_chart
    {X L r : ℝ}
    (hX : X ≠ 0)
    (hH : holonomyEquation X L r = 0)
    (hE : edgeEquation X r = 0) :
    aPolynomialX X L = 0 := by
  have hL := holonomyEquation_eq hH
  have hcert := aPolynomialX_substitution X L r hL
  have hD : chartT r - r * X ≠ 0 := by
    intro hD
    have hprod : X * r ^ 4 * chartT r = 0 := by
      unfold edgeEquation at hE
      rw [hD] at hE
      linarith
    have hprod' : X * (r ^ 4 * chartT r) = 0 := by
      simpa [mul_assoc] using hprod
    have hrt : r ^ 4 * chartT r = 0 :=
      (mul_eq_zero.mp hprod').resolve_left hX
    rcases mul_eq_zero.mp hrt with hr4 | ht
    · have hr : r = 0 := eq_zero_of_pow_eq_zero hr4
      subst r
      norm_num [chartT] at hD
    · have hrX : r * X = 0 := by linarith
      have hr : r = 0 := (mul_eq_zero.mp hrX).resolve_right hX
      subst r
      norm_num [chartT] at ht
  have hDpow : (chartT r - r * X) ^ 5 ≠ 0 := pow_ne_zero _ hD
  have hcert' :
      (chartT r - r * X) ^ 5 * aPolynomialX X L = 0 := by
    simpa [hE] using hcert
  exact (mul_eq_zero.mp hcert').resolve_left hDpow
