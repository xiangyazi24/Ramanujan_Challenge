import Mathlib.Tactic

/-!
# The `7₂` A-polynomial and its shape-chart equations

This file fixes the polynomial from the official statement of Ramanujan
Challenge Problem 3.1.  We use `X = M²`, since every exponent of the meridian
variable in the displayed A-polynomial is even.

The other definitions are the two polynomial equations of the four-tetrahedron
shape chart:

* `holonomyEquation X L r = 0` is the cleared holonomy equation;
* `edgeEquation X r = 0` is the internal edge equation.

All definitions are polynomial: later elimination certificates can therefore be
checked by `ring` without square roots, resultants, or external CAS calls.
-/

/-- The official `7₂` A-polynomial after the substitution `X = M²`. -/
def aPolynomialX (X L : ℝ) : ℝ :=
  L ^ 5
    + L ^ 4 * (X ^ 7 - X ^ 6 + 3 * X ^ 2 + 4 * X - 2)
    + L ^ 3 *
      (-2 * X ^ 9 + 5 * X ^ 8 + X ^ 7 - 4 * X ^ 6
        + 6 * X ^ 4 + 5 * X ^ 3 + 2 * X ^ 2 - 4 * X + 1)
    + L ^ 2 *
      (X ^ 11 - 4 * X ^ 10 + 2 * X ^ 9 + 5 * X ^ 8 + 6 * X ^ 7
        - 4 * X ^ 5 + X ^ 4 + 5 * X ^ 3 - 2 * X ^ 2)
    + L * (-2 * X ^ 11 + 4 * X ^ 10 + 3 * X ^ 9 - X ^ 5 + X ^ 4)
    + X ^ 11

/-- The official polynomial in the variables used in the problem statement. -/
def aPolynomial (M L : ℝ) : ℝ := aPolynomialX (M ^ 2) L

/-- `t = 1-r²`, the first tetrahedron shape. -/
def chartT (r : ℝ) : ℝ := 1 - r ^ 2

/-- The cleared holonomy equation
`(L+X³)(1-r²) - Xr(L+X) = 0`. -/
def holonomyEquation (X L r : ℝ) : ℝ :=
  (L + X ^ 3) * chartT r - X * r * (L + X)

/-- The internal edge equation
`(t-r)(Xt-r)(t-rX) = Xr⁴t`, with `t=1-r²`. -/
def edgeEquation (X r : ℝ) : ℝ :=
  (chartT r - r) * (X * chartT r - r) * (chartT r - r * X)
    - X * r ^ 4 * chartT r

/-- Solving the holonomy equation for `L`, with no division performed. -/
theorem holonomyEquation_eq
    {X L r : ℝ} (h : holonomyEquation X L r = 0) :
    L * (chartT r - r * X) = X ^ 2 * (r - X * chartT r) := by
  unfold holonomyEquation at h
  linear_combination h

/-- The shape variables used by the chart, expressed without a square root. -/
noncomputable def shapeU (r : ℝ) : ℝ := r / chartT r
noncomputable def shapeV (X r : ℝ) : ℝ := r / (chartT r * X)
noncomputable def shapeW (X r : ℝ) : ℝ := chartT r / (chartT r - r * X)

/-- The multiplicative edge equation in the form used for logarithms. -/
theorem edgeEquation_eq
    {X r : ℝ} (h : edgeEquation X r = 0) :
    (chartT r - r) * (X * chartT r - r) * (chartT r - r * X)
      = X * r ^ 4 * chartT r := by
  unfold edgeEquation at h
  linarith
