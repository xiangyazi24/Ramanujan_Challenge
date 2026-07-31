import Ramanujan31.APolynomial

/-!
# A nonvanishing certificate for the shape-chart subresultant

Eliminating `r` from the holonomy and edge equations produces a linear
subresultant `ca(X,L) r + cb(X,L)`.  Reconstructing `r` from the official
A-polynomial requires knowing that `ca` is nonzero on the desired real chamber.

The large polynomial `ca` has no useful sign on a rectangular box.  On the
shape curve, however, it has a short exact factor certificate.  Put
`D = 1-r²-rX`.  After substituting the holonomy expression for `L`,

`D⁴ ca = X⁶(X-1)⁴(X+1)⁴ F`,

and polynomial division by the edge equation gives

`rF = -(r⁴-r²+1)E + X(r²-1)²`.

Thus `E=0`, `X>0`, and `r<-1` force `F<0`, hence `ca<0`.  Both displayed
identities are checked below by `ring`.
-/

def subresultantA4 (X : ℝ) : ℝ :=
  X ^ 4 + 3 * X ^ 2 + 1

def subresultantA3 (X : ℝ) : ℝ :=
  -X ^ 6 + 7 * X ^ 5 + X ^ 4 + 11 * X ^ 3 + X ^ 2 + 2 * X - 1

def subresultantA2 (X : ℝ) : ℝ :=
  X ^ 8 - 4 * X ^ 7 + 14 * X ^ 6 + 8 * X ^ 5 + 14 * X ^ 4
    - 4 * X ^ 3 + X ^ 2

def subresultantA1 (X : ℝ) : ℝ :=
  -X ^ 10 + 2 * X ^ 9 + X ^ 8 + 11 * X ^ 7 + X ^ 6
    + 7 * X ^ 5 - X ^ 4

def subresultantA0 (X : ℝ) : ℝ :=
  X ^ 10 + 3 * X ^ 8 + X ^ 6

/-- Coefficient of `r` in the linear subresultant of the holonomy and edge
equations. -/
def subresultantA (X L : ℝ) : ℝ :=
  L ^ 4 * subresultantA4 X
    + L ^ 3 * subresultantA3 X
    + L ^ 2 * subresultantA2 X
    + L * subresultantA1 X
    + subresultantA0 X

/-- The short remainder factor in the sign certificate. -/
def subresultantFactor (X r : ℝ) : ℝ :=
  X ^ 2 * r ^ 8 + X ^ 2 * r ^ 7 - 3 * X ^ 2 * r ^ 6
    - 2 * X ^ 2 * r ^ 5 + 4 * X ^ 2 * r ^ 4 + 2 * X ^ 2 * r ^ 3
    - 3 * X ^ 2 * r ^ 2 - X ^ 2 * r + X ^ 2
    + X * r ^ 8 - X * r ^ 7 - 2 * X * r ^ 6 + 3 * X * r ^ 5
    + 3 * X * r ^ 4 - 3 * X * r ^ 3 - 2 * X * r ^ 2 + X * r + X
    + r ^ 8 + r ^ 7 - 3 * r ^ 6 - 2 * r ^ 5 + 4 * r ^ 4
    + 2 * r ^ 3 - 3 * r ^ 2 - r + 1

/-- Polynomial-division certificate for `subresultantFactor`. -/
theorem subresultantFactor_mul_r
    (X r : ℝ) :
    r * subresultantFactor X r
      = -(r ^ 4 - r ^ 2 + 1) * edgeEquation X r
        + X * (r ^ 2 - 1) ^ 2 := by
  unfold subresultantFactor edgeEquation chartT
  ring

/-- Denominator-free substitution certificate for `subresultantA`.

The equation `hL` is exactly the cleared holonomy formula
`L D = X²(r-Xt)`. -/
theorem subresultantA_substitution
    (X L r : ℝ)
    (hL : L * (chartT r - r * X) = X ^ 2 * (r - X * chartT r)) :
    (chartT r - r * X) ^ 4 * subresultantA X L
      = X ^ 6 * (X - 1) ^ 4 * (X + 1) ^ 4
          * subresultantFactor X r := by
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
  calc
    (chartT r - r * X) ^ 4 * subresultantA X L
        =
          subresultantA4 X *
              (L ^ 4 * (chartT r - r * X) ^ 4)
            + subresultantA3 X * (chartT r - r * X) *
              (L ^ 3 * (chartT r - r * X) ^ 3)
            + subresultantA2 X * (chartT r - r * X) ^ 2 *
              (L ^ 2 * (chartT r - r * X) ^ 2)
            + subresultantA1 X * (chartT r - r * X) ^ 3 *
              (L * (chartT r - r * X))
            + subresultantA0 X * (chartT r - r * X) ^ 4 := by
              unfold subresultantA
              ring
    _ =
          subresultantA4 X * (X ^ 8 * (r - X * chartT r) ^ 4)
            + subresultantA3 X * (chartT r - r * X) *
              (X ^ 6 * (r - X * chartT r) ^ 3)
            + subresultantA2 X * (chartT r - r * X) ^ 2 *
              (X ^ 4 * (r - X * chartT r) ^ 2)
            + subresultantA1 X * (chartT r - r * X) ^ 3 *
              (X ^ 2 * (r - X * chartT r))
            + subresultantA0 X * (chartT r - r * X) ^ 4 := by
              rw [h4, h3, h2, hL]
    _ = X ^ 6 * (X - 1) ^ 4 * (X + 1) ^ 4
          * subresultantFactor X r := by
            unfold subresultantA4 subresultantA3 subresultantA2
              subresultantA1 subresultantA0 subresultantFactor chartT
            ring

/-- The linear-subresultant coefficient is strictly negative throughout the
real shape chamber. -/
theorem subresultantA_neg_of_shape_chamber
    {X L r : ℝ}
    (hX : 0 < X) (hX1 : X < 1) (hr : r < -1)
    (hH : holonomyEquation X L r = 0)
    (hE : edgeEquation X r = 0) :
    subresultantA X L < 0 := by
  have hr0 : r < 0 := by linarith
  have hr2 : 1 < r ^ 2 := by nlinarith
  have hfac :
      r * subresultantFactor X r = X * (r ^ 2 - 1) ^ 2 := by
    rw [subresultantFactor_mul_r, hE]
    ring
  have hfacneg : subresultantFactor X r < 0 := by
    have hsquare : 0 < (r ^ 2 - 1) ^ 2 := sq_pos_of_pos (by linarith)
    have : 0 < X * (r ^ 2 - 1) ^ 2 := mul_pos hX hsquare
    nlinarith
  have hL := holonomyEquation_eq hH
  have hcert := subresultantA_substitution X L r hL
  have hXm1 : X - 1 ≠ 0 := by linarith
  have hXp1 : X + 1 ≠ 0 := by linarith
  have hright :
      X ^ 6 * (X - 1) ^ 4 * (X + 1) ^ 4
          * subresultantFactor X r < 0 := by
    have hpowX : 0 < X ^ 6 := pow_pos hX _
    have hpowm : 0 < (X - 1) ^ 4 := by positivity
    have hpowp : 0 < (X + 1) ^ 4 := by positivity
    exact mul_neg_of_pos_of_neg
      (mul_pos (mul_pos hpowX hpowm) hpowp) hfacneg
  have hleft :
      (chartT r - r * X) ^ 4 * subresultantA X L < 0 := by
    rw [hcert]
    exact hright
  have hden : 0 < (chartT r - r * X) ^ 4 := by
    have hne : (chartT r - r * X) ^ 4 ≠ 0 := by
      intro hz
      rw [hz] at hleft
      linarith
    have hnonneg : 0 ≤ (chartT r - r * X) ^ 4 := by positivity
    exact lt_of_le_of_ne hnonneg (Ne.symm hne)
  nlinarith

/-- In particular, the coefficient never vanishes on the real shape chamber. -/
theorem subresultantA_ne_zero_of_shape_chamber
    {X L r : ℝ}
    (hX : 0 < X) (hX1 : X < 1) (hr : r < -1)
    (hH : holonomyEquation X L r = 0)
    (hE : edgeEquation X r = 0) :
    subresultantA X L ≠ 0 :=
  ne_of_lt (subresultantA_neg_of_shape_chamber hX hX1 hr hH hE)
