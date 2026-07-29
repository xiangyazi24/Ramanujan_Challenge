#!/usr/bin/env python3
"""Run and combine all Round 2 verification experiments.

The component reports retain the full per-prime coefficient tables and the
100-point Hodge-slope grid.  This wrapper adds a short integrated verdict and
writes the requested combined artifact to /tmp/round2_results.txt by default.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUTPUT = Path("/tmp/round2_results.txt")
PARTS = (
    Path("/tmp/round2_v1_v2.txt"),
    Path("/tmp/round2_v3_v4.txt"),
    Path("/tmp/round2_v5.txt"),
)


INTEGRATED_VERDICT = """\
ROUND 2 INTEGRATED VERDICT
==========================

Scope and conventions
---------------------
All b_j tables use b_j=sum_k binom(j,k)^2 binom(j+k,k)^2 and were
independently cross-checked against the Apery recurrence.  In V3--V4, Delta
means the Adolphson--Sperber polytope conv({0} union Supp(f)); the ordinary
conv(Supp(f)) is only its three-dimensional octahedral base.  The V4 verdict
uses the character chi_w=omega^(-j), hence shift theta=(a,0,0,0) with
a=j/(p-1).  Changing the twisted coordinate or the character sign changes
the answer.  The finite Mellin inversion used to motivate V1 applies to the
interior indices 1<=j<=p-2; at the two endpoints the trivial character sees
both b_0 and b_(p-1).  All requested endpoint coefficients are still tabulated.

Make-or-break findings
----------------------
1. V1 is not presently an independently testable trace identity.  Without a
   separately constructed pure object M_j and normalized unit root alpha(j),
   the displayed identity merely defines alpha(j)=-residual(j).  The task text
   also has a sign mismatch: the displayed star makes alpha=-residual, not
   alpha=residual.  Its separately requested zero implication is false:
     p=31, j=8,22:  b_j=0 but residual=11;
     p=73, j=2,70:  b_j=0 but residual=30;
     p=97, j=25,71: b_j=0 but residual=39.
   Thus zeros cannot all be identified with residual=0 under the proposed
   boundary subtraction.  The phrase "smoothly in j" has no finite-field
   criterion and was not treated as a test.  Moreover, a literal Frobenius
   unit root has nonzero reduction, whereas the formula infers alpha=0 at
   boundary-free zeros such as p=19, j=8,10.  Hence alpha must be undefined
   there or must mean a Hasse-type scalar rather than an actual unit root.

2. V2 is confirmed, with the exact constant +1: b_(p-1-j)=b_j for every j
   and every one of the 20 requested primes.  Every defined reversal ratio is
   1.  The report also gives the two one-line binomial congruences proving this
   for every odd prime.  This explains paired coefficient zeros.

3. V3 passes at every challenge prime p>=5.  Delta has f-vector
   (1,7,18,20,9,1), including the empty face and Delta.  All 27 nonempty
   faces not containing 0 are nondegenerate.  The full base face is degenerate
   only in characteristics 2 and 3.

4. V4 passes the stated necessary first-slope test under the intended
   w-Kummer/negative-Teichmuller convention.  The complete Hodge slopes are
     a (x1), a+1 (x3), a+2 (x3), a+3 (x1).
   Hence all 100 formal values a=k/100 have first slope a.  This does not prove
   ordinarity: the A-S Hodge polygon is a lower bound, while the blueprint still
   needs equality of Newton and Hodge polygons (or equivalent control) uniformly
   in the twists.  The scalar a alone is insufficient data: the opposite sign
   gives first slope 1-a for a>0, and an x-coordinate twist gives first slope 1.

5. V5 gives the exact random-square expectation
     E[Z_square]=(p-1)/p-p^(-(floor((p-1)/4)+1)),
   which tends to 1.  The 4,000-sample exact NTT/CRT simulation has pooled mean
   0.996250 and population variance 0.951236, close to Poisson(1).  Apery data
   for all 667 primes 5<=p<=5000 have mean 0.911544 and population variance
   1.768787, close instead to 2*Poisson(1/2) because of V2.  Squareness explains
   the order-one scale; reversal explains the pairing.  A Poisson heuristic has
   unbounded support, so this experiment is not evidence for a uniform constant
   C in Hypothesis Z.  The literal random-monic experiment also omits two actual
   constraints: A_p(0)^2=1 and A_p^*=+/-A_p.  Conditioning only on nonzero
   constant term changes the small finite-p correction, but not its limit 1;
   the reciprocity constraint is what the generic simulation most visibly omits.

Integration with E1--E4
-----------------------
E1's uniform zero positions rule out the proposed fixed torsion/CM route.
E4's lack of a mod-8 signal rules out boundary behavior as the dominant source.
E3 found no multiplicative or separated Casoratian structure to replace those
routes.  Round 2 leaves the geometric prerequisites (V2--V4) intact, but the
central trace decomposition/ordinarity mechanism is still missing: V1 has no
independent alpha(j), its requested zero criterion has explicit counterexamples,
and a first Hodge slope equal to a does not control Newton-polygon equality or
the number of nonordinary twists.  The blueprint therefore does not yet prove
Hypothesis Z.

The complete machine-auditable reports follow.
"""


def run_component(script: str, *arguments: str) -> None:
    subprocess.run(
        [sys.executable, str(ROOT / "scripts" / script), *arguments],
        cwd=ROOT,
        stdout=subprocess.DEVNULL,
        check=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    run_component("round2_v1_v2.py", "--output", str(PARTS[0]))
    run_component("round2_v3_v4.py")
    run_component("round2_v5.py", "--output", str(PARTS[2]))

    sections = [INTEGRATED_VERDICT.rstrip()]
    sections.extend(part.read_text(encoding="utf-8").rstrip() for part in PARTS)
    args.output.write_text("\n\n".join(sections) + "\n", encoding="utf-8")
    print(f"Wrote {args.output} ({args.output.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
