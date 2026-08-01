import Ramanujan31.APolynomial

/-!
# The direct Rogers-regulator differential

This is the algebraic core of a direct proof of Problem 3.1.  It replaces an
appeal to the Neumann--Zagier variation theorem by a scalar identity.

For a nonzero real shape `z`, write

`rogersDifferential log|z| log|1-z| dlog|z| dlog|1-z|`

for the differential of the real Rogers dilogarithm.  Substituting the four
shape-chart logarithms gives

`Σ dR = -½ (log X dlog L - log L dlog X)`

plus an error term proportional to the logarithm and differential of the edge
equation.  On the edge-equation locus that error is zero.  Since `X=M²`, the
remaining term is exactly the one-form in the official integral.

The theorem below is pure commutative-ring algebra.  The analytic theorem saying
that it is the derivative of the Rogers function is kept separate.
-/

/-- The standard differential of
`Li₂(z) + ½ log|z| log|1-z|`. -/
def rogersDifferential
    (logz logOneSubz dlogz dlogOneSubz : ℝ) : ℝ :=
  (logz * dlogOneSubz - logOneSubz * dlogz) / 2

/-- Exact four-shape differential identity, including the edge-equation error
term.

The logarithm coordinates are

* `R = log|r|`, `T = log|1-r²|`,
* `A = log|1-r²-r|`,
* `X = log X`,
* `P = log|X(1-r²)-r|`,
* `Q = log|1-r²-rX|`.

Prefixing a coordinate by `d` denotes its derivative along an arbitrary
parameter. -/
theorem four_shape_rogersDifferential
    (R T A X P Q dR dT dA dX dP dQ : ℝ) :
    let logL := 2 * X + P - Q
    let dlogL := 2 * dX + dP - dQ
    let edgeLog := A + P + Q - 4 * R - T - X
    let dEdgeLog := dA + dP + dQ - 4 * dR - dT - dX
    rogersDifferential T (2 * R) dT (2 * dR)
      + rogersDifferential (R - T) (A - T) (dR - dT) (dA - dT)
      + rogersDifferential (R - T - X) (P - T - X)
          (dR - dT - dX) (dP - dT - dX)
      + rogersDifferential (T - Q) (R + X - Q)
          (dT - dQ) (dR + dX - dQ)
      =
        -(X * dlogL - logL * dX) / 2
          + ((R - T) * dEdgeLog - edgeLog * (dR - dT)) / 2 := by
  simp only
  unfold rogersDifferential
  ring

/-- On the edge-equation locus, the error term in
`four_shape_rogersDifferential` vanishes. -/
theorem four_shape_rogersDifferential_of_edge
    (R T A X P Q dR dT dA dX dP dQ : ℝ)
    (hedge : A + P + Q - 4 * R - T - X = 0)
    (dhedge : dA + dP + dQ - 4 * dR - dT - dX = 0) :
    let logL := 2 * X + P - Q
    let dlogL := 2 * dX + dP - dQ
    rogersDifferential T (2 * R) dT (2 * dR)
      + rogersDifferential (R - T) (A - T) (dR - dT) (dA - dT)
      + rogersDifferential (R - T - X) (P - T - X)
          (dR - dT - dX) (dP - dT - dX)
      + rogersDifferential (T - Q) (R + X - Q)
          (dT - dQ) (dR + dX - dQ)
      = -(X * dlogL - logL * dX) / 2 := by
  simp only
  rw [four_shape_rogersDifferential]
  rw [hedge, dhedge]
  ring

/-- The factor-of-two conversion from `X=M²` to the meridian coordinate `M`. -/
theorem meridian_square_differential
    (logM logL dlogM dlogL : ℝ) :
    -((2 * logM) * dlogL - logL * (2 * dlogM)) / 2
      = -(logM * dlogL - logL * dlogM) := by
  ring

