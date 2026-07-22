/-
  Problem 3.1 — Layer 4: Regulator Certificate.

  The regulator difference Δℛ = Σ_j [R(z_j(β)) − R(z_j(α))]
  must equal −4π²/85.

  STRATEGY: Use the extended Rogers dilogarithm (rogersGtOne, rogersNeg)
  to express R values at z > 1 and z < 0 in terms of R at arguments
  in (0,1). The Li₂ values at algebraic arguments cancel in the sum
  (Bloch group torsion), leaving only rational multiples of π².

  The shapes at each endpoint satisfy:
    T < 0:  Use rogersNeg (Landen identity)
    U > 1:  Use rogersGtOne (inversion)
    V > 1:  Use rogersGtOne (inversion)
    W > 1:  Use rogersGtOne (inversion)
-/
import RamanujanChallenge.Dilogarithm
import RamanujanChallenge.Problem31.EndpointData

noncomputable section

open Real

/-! ## The regulator sum using extended Rogers

For each endpoint, the regulator contribution is:
  rogersNeg(T) + rogersGtOne(U) + rogersGtOne(V) + rogersGtOne(W)

Each of these is defined purely in terms of dilog and log at
arguments in (0,1), plus rational multiples of π².
-/

def regulatorContrib (T U V W : ℝ) : ℝ :=
  rogersNeg T + rogersGtOne U + rogersGtOne V + rogersGtOne W

/-! ## The regulator certificate

The master theorem: the regulator difference equals −4π²/85.
-/

theorem regulator_certificate
    (s_β : ℝ) (hs_β0 : 0 < s_β) (hs_β1 : s_β < 1)
    (s_α : ℝ) (hs_α0 : 0 < s_α) (hs_α1 : s_α < 1) :
    let shapes_β := endpointShapes_beta s_β hs_β0 hs_β1
    let shapes_α := endpointShapes_alpha s_α hs_α0 hs_α1
    regulatorContrib shapes_β.1 shapes_β.2.1 shapes_β.2.2.1 shapes_β.2.2.2 -
    regulatorContrib shapes_α.1 shapes_α.2.1 shapes_α.2.2.1 shapes_α.2.2.2 =
    -(4 * Real.pi ^ 2 / 85) := by
  sorry

end
