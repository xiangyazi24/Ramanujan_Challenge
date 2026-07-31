/-
Ramanujan Challenge, Problem 3.1 — the machine-checked part.

Building this library checks every theorem cited in `solution.pdf` as
"machine-checked".  Start at `Ramanujan31.MainTheorem`, which composes the
others and documents exactly which inputs remain external.

  Ramanujan31/RatReconstruct.lean    rational reconstruction from a denominator
                                     bound (the final step, Thm 5.2)
  Ramanujan31/ChartSymmetry.lean     u(1/a) = u(a); the palindromic
                                     decompositions f = a^d g(a + 1/a)
  Ramanujan31/ShapeField.lean        algebraic certificate that 1 + 4u^2 is a
                                     SQUARE in the endpoint field, so the shapes
                                     are defined over F (no quadratic extension)
  Ramanujan31/TraceRoots.lean        the trace polynomials are totally real,
                                     with EXACT root counts in (-2,2)
  Ramanujan31/UnitCircle.lean        real trace in [-2,2]  <=>  |a| = 1
  Ramanujan31/ShapeCancellation.lean the Bloch-Wigner four-shape cancellation
  Ramanujan31/MainTheorem.lean       the composition, plus the final step

No `sorry`, no `native_decide`.  Every theorem depends only on
`propext`, `Classical.choice`, `Quot.sound`.
-/
import Ramanujan31.Dilog.Basic
import Ramanujan31.Dilog.RealBounds
import Ramanujan31.Dilog.Rogers
import Ramanujan31.Dilog.FiveTerm
import Ramanujan31.Dilog.SpecialValues
import Ramanujan31.Dilog.EndpointTransform
import Ramanujan31.Dilog.CyclicRogers
import Ramanujan31.Dilog.IdealPolygon
import Ramanujan31.Dilog.Regular17
import Ramanujan31.Dilog.Fan17
import Ramanujan31.Dilog.Certify
import Ramanujan31.Dilog.BlochWigner
import Ramanujan31.Dilog.Instance
import Ramanujan31.APolynomial
import Ramanujan31.Elimination
import Ramanujan31.ReverseElimination
import Ramanujan31.RegulatorDifferential
import Ramanujan31.SubresultantSign
import Ramanujan31.ShapeChamber
import Ramanujan31.RatReconstruct
import Ramanujan31.ChartSymmetry
import Ramanujan31.ShapeField
import Ramanujan31.TraceRoots
import Ramanujan31.UnitCircle
import Ramanujan31.ShapeCancellation
import Ramanujan31.MainTheorem
