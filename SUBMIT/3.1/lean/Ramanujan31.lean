/-
Ramanujan Challenge, Problem 3.1 — the machine-checked part.

Building this library checks every theorem cited in `solution.pdf` as
"machine-checked".  Start at `Ramanujan31.MainTheorem`, which composes the
others and documents exactly which inputs remain external.

  Ramanujan31/RatReconstruct.lean    rational reconstruction from a denominator
                                     bound (the final step, Thm 5.2)
  Ramanujan31/ChartSymmetry.lean     u(1/a) = u(a); the palindromic
                                     decompositions f = a^d g(a + 1/a)
  Ramanujan31/TraceRoots.lean        the trace polynomials are totally real,
                                     with EXACT root counts in (-2,2)
  Ramanujan31/UnitCircle.lean        real trace in [-2,2]  <=>  |a| = 1
  Ramanujan31/ShapeCancellation.lean the Bloch-Wigner four-shape cancellation
  Ramanujan31/MainTheorem.lean       the composition, plus the final step

No `sorry`, no `native_decide`.  Every theorem depends only on
`propext`, `Classical.choice`, `Quot.sound`.
-/
import Ramanujan31.RatReconstruct
import Ramanujan31.ChartSymmetry
import Ramanujan31.TraceRoots
import Ramanujan31.UnitCircle
import Ramanujan31.ShapeCancellation
import Ramanujan31.MainTheorem
