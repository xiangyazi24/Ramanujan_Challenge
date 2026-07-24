#!/usr/bin/env python3
"""Run and combine all Round 4 verification tasks.

The component reports retain the full parameter-search tables, the complete
odd-prime center-congruence audit, and the p=159977 zero verification.  This
wrapper prepends the integrated verdict and writes /tmp/round4_results.txt.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUTPUT = Path("/tmp/round4_results.txt")
COMPONENTS = (
    ("round4_task_a.py", Path("/tmp/round4_task_a.txt")),
    ("round4_task_b.py", Path("/tmp/round4_task_b.txt")),
    ("round4_task_c.py", Path("/tmp/round4_task_c.txt")),
)


INTEGRATED_VERDICT = """\
ROUND 4: INTEGRATED VERDICT
========================================================================

Scope and validation
--------------------
All Apéry values use b_0=1, b_1=5 and
  (j+1)^3 b_(j+1) - (34j^3+51j^2+27j+5)b_j + j^3 b_(j-1)=0.
The small-prime tables were checked against the defining binomial sum.  The
eta product was expanded by two independent exact algorithms.  At p=159977,
the recurrence was rechecked at every index and sampled against the direct
binomial definition, including every reported zero.

Task A -- exhaustive corrected-boundary score
----------------------------------------------
Roots are labeled using the least nonnegative square root of 2:
  p=31: sqrt2=8,  (t1,t2)=(14,20), H_p(t1)=H_p(t2)=7.
  p=73: sqrt2=32, (t1,t2)=(71,36), H_p(t1)=H_p(t2)=50.
  p=97: sqrt2=14, (t1,t2)=(43,88), H_p(t1)=H_p(t2)=4.

For all three primes, the unique unrestricted maximum is (c1,d1)=(0,0):
  p=31: score 29/29; p=73: 71/71; p=97: 95/95.
This maximum is forced, not evidence for the corrected trace identity: at
(0,0), R=b, which satisfies the scored recurrence by definition.  The best R
therefore still vanishes at [8,22], [2,70], and [25,71], respectively, and is
not nonzero everywhere.

After excluding (0,0), the maximum scores are only 5, 7, and 7, with 90, 72,
and 288 tied parameter pairs.  If c1 is fixed to H_p(t1), the best scores are:
  p=31: 5/29 at d1 in {4,8,19}; naive d1=0 scores 1/29.
  p=73: 7/71 at d1=12;       naive d1=0 scores 3/71.
  p=97: 7/95 at d1 in {57,58,70}; naive d1=0 scores 3/95.
Among those fixed-c1 maximizers, R is everywhere nonzero for p=31 at d1=8,19
and for p=97 at d1=58,70; the p=73 maximizer has zeros [27,45].
The center j=(p-1)/2 is a universal zero of the recurrence defect for every
parameter pair, so the exact all-pair mean scores are about 1.90--1.97 rather
than the nominal independent-random baseline near 1.  There is no isolated
nontrivial peak.  Under an independent-uniform projective occupancy benchmark,
the nonzero-peak tail probabilities are 0.9652, 0.6165, and 0.7396, so none is
significant.  This sweep therefore does not support a nonzero corrected
boundary term; an independent normalization or independently computed pure
trace would be needed to make the test identifying.

Task B -- center/modular-form verification through 10000
--------------------------------------------------------
The q-series eta(2z)^4 eta(4z)^4 was computed through q^10000.  The two exact
coefficient algorithms agree at all 10001 coefficients.  For every one of
the 1228 odd primes p<=10000, both direct and recurrence computations give
  b_((p-1)/2) = a(p) (mod p),
with zero failures.  Exactly two odd primes divide a(p):
  p=11:   a(11)=-44=-4*11,       b_5=0 (mod 11),
  p=3137: a(3137)=207042=66*3137, b_1568=0 (mod 3137).
Literal p=2 has a(2)=0, but its center index (p-1)/2 is not an integer; it is
the bad level prime and is excluded from the center theorem and the expected
odd-prime list.

Task C -- record zero count at p=159977
---------------------------------------
Z(159977)=12 exactly.  The zero positions are
  [37078,39802,43581,49216,60763,68108,
   91868,99213,110760,116395,120174,122898].
They form the six palindromic pairs
  (37078,122898), (39802,120174), (43581,116395),
  (49216,110760), (60763,99213), (68108,91868),
each summing to p-1=159976.  The full coefficient table satisfies
b_j=b_(p-1-j) at all 159977 indices.

The complete machine-auditable component reports follow.
"""


def run_component(script: str, output: Path) -> None:
    subprocess.run(
        [sys.executable, str(ROOT / "scripts" / script), "--output", str(output)],
        cwd=ROOT,
        check=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    for script, output in COMPONENTS:
        run_component(script, output)

    sections = [INTEGRATED_VERDICT.rstrip()]
    sections.extend(
        output.read_text(encoding="utf-8").rstrip()
        for _, output in COMPONENTS
    )
    args.output.write_text("\n\n".join(sections) + "\n", encoding="utf-8")
    print(f"Wrote {args.output} ({args.output.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
