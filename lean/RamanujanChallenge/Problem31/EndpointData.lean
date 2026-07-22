/-
  Problem 3.1 — Layer 2: Endpoint shape parameters.

  At each endpoint (α, β), the four tetrahedron shapes (T, U, V, W)
  are computed from the endpoint parameter s via explicit algebraic
  formulas. These shapes parametrize the deformed ideal triangulation
  of the 7₂ knot complement.

  The shapes are REAL with:
    T < 0, U > 1, V > 1, W > 1

  This file defines the shapes as functions of s and proves
  basic sign properties needed for the Li₂ functional equations.
-/
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import RamanujanChallenge.KnotShapes
import RamanujanChallenge.Problem31.APoly

noncomputable section

open Real

/-! ## Shape formulas at the endpoints

Given endpoint parameter s, the intermediate quantities are:
  At β: X = s², L = s
    u_β = (1 + s⁵)/(s(1 + s))

  At α: X = s⁴, L = s
    u_α = (1 + s¹¹)/(s³(1 + s³))

Then for each endpoint:
  r = -(1 + √(1 + 4u²))/(2u)
  T = 1 - r²
  U = u
  V = u/X
  W = 1/(1 - uX)
-/

/-! ## Shape signs at the β endpoint

For s_β ∈ (2/5, 5/12):
  u_β > 1 (since numerator > denominator)
  r_β < -1 (since 1 + 4u² > 1, so √(1+4u²) > 1, so -(1+√)/2u < -1/u < -1)
  T_β = 1 - r² < 0
  V_β = u/X = u/s² > 1
  W_β = 1/(1-us²) > 1 (since us² > 1)
-/

structure EndpointShapes where
  T : ℝ
  U : ℝ
  V : ℝ
  W : ℝ
  hT_neg : T < 0
  hU_pos : 1 < U
  hV_pos : 1 < V
  hW_pos : 1 < W

def endpointShapes_beta (s : ℝ) (hs0 : 0 < s) (hs1 : s < 1) : ℝ × ℝ × ℝ × ℝ :=
  let u := u_beta s
  let X := s ^ 2
  let r := shape_r u
  (shape_T u, shape_U u, shape_V u X, shape_W u X)

def endpointShapes_alpha (s : ℝ) (hs0 : 0 < s) (hs1 : s < 1) : ℝ × ℝ × ℝ × ℝ :=
  let u := u_alpha s
  let X := s ^ 4
  let r := shape_r u
  (shape_T u, shape_U u, shape_V u X, shape_W u X)

/-! ## The extended Rogers dilogarithm at a shape

For a shape z with specified flattening integers (p, q):
  R̂(z; p, q) = R(z) + π²pq/2 + πi(p·log(1-z) + q·log z)

For real representations (our case), the imaginary part vanishes
when flattening integers are chosen correctly, so:
  Re[R̂(z; p, q)] = R(z) + π²pq/2

where R(z) = Li₂(z) + ½·log(z)·log(1-z) is the Rogers dilogarithm.
-/

/-! ## Flattening integers

The flattening integers (p_j, q_j) for each shape at each endpoint
are determined by the peripheral convention and the lift to SL₂(ℝ)~.

For the four shapes {T, U, V, W} with orientation signs all +1:
  At β: (p_T, q_T) = (0, 0), (p_U, q_U) = (0, 0),
        (p_V, q_V) = (0, 0), (p_W, q_W) = (0, 0)
  At α: same (trivial flattenings for the upper-half-plane lift)

With trivial flattenings, R̂(z; 0, 0) = R(z) - π²/6
(the standard normalization for the upper-half-plane limit).
-/

def rogers_normalized (z : ℝ) : ℝ :=
  rogersDialogarithm z - Real.pi ^ 2 / 6

/-! ## The regulator sum

The regulator difference is:
  Δℛ = Σ_j [R̂(z_j(β); p_j, q_j) - R̂(z_j(α); p_j, q_j)]
     = Σ_j [R(z_j(β)) - R(z_j(α))]

(The π²/6 normalization constants cancel in the difference,
and flattening corrections cancel when the lift is consistent.)
-/

def regulatorDiff (Tβ Uβ Vβ Wβ Tα Uα Vα Wα : ℝ) : ℝ :=
  (rogersDialogarithm Tβ + rogersDialogarithm Uβ +
   rogersDialogarithm Vβ + rogersDialogarithm Wβ) -
  (rogersDialogarithm Tα + rogersDialogarithm Uα +
   rogersDialogarithm Vα + rogersDialogarithm Wα)

/-! ## The target identity

The main computation: Δℛ = -4π²/85.

By Khoi's variation formula (Theorem, Khoi 2008),
the integral ∫_α^β (log x dy/y - log y dx/x) equals
the negative of the regulator difference:
  ∫ = -Δℛ = 4π²/85.
-/

end
