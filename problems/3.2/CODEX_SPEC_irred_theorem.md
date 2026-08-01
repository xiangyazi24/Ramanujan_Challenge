# CODEX SPEC — general irreducibility + branch-disjointness (turn S5 pattern into proof)

Report to `CODEX_IRRED_THEOREM_report.md`. This is a THEOREM-PROVING task with
symbolic verification support (sympy). Inputs: your previous report
`CODEX_COINC_SYMBOLIC_report.md` (addition laws U1/U2/D-update, component census
h<=8, gcd(C_h,C_k)=1 data) + `campaign3_questions/CTX.txt`.

## Goal
Prove for ALL 1 <= h < k (not just <=8):
 (T1) F_{h,k}(x,y) = N_h(x)D_k(y) - N_k(y)D_h(x) is irreducible over Qbar.
 (T2) The branch loci of delta_h and delta_k in P^1_value are disjoint apart from
      the common value 0 and infinity (i.e. gcd(C_h, C_k)=1), OR characterize failures.
 (T3) Same-gap: F_{h,h}/(x-y) = H_h irreducible over Qbar for all h>=2.

## Suggested route (attack, adapt, or replace)
- Derivative of the addition law: differentiate D-update to get A_{h+d} in terms of
  lower A's and deltas; combine with the Apery-product polar data
  N_h(-j)=(-1)^{j-1}((j-1)!)^3 b_{j-1}((h-j)!)^3 b_{h-j} to control common roots.
- Local method for (T1): at the pole x=-j the function delta_h has polar part
  b_{j-1}b_{h-j}*(unit)/(x+j)^3. A factorization of F_{h,k} would force a grouping of
  the 3-sheeted local branches over the infinity value compatible on BOTH sides;
  use valuation/Newton polygon at (x,y)=(-j,-i) pairs + the explicit leading constants
  to derive a contradiction unless the grouping is trivial. Verify each claimed local
  computation symbolically for (h,k)=(2,3),(3,5) as you go.
- Monodromy route for (T1)+(T3): the cover delta_h has 3-cycles over infinity at h
  points and simple transpositions at 4h-4 generic branch points (verified squarefree
  crit structure h<=8; prove squarefreeness in general or condition on it): show the
  monodromy group is the full symmetric group S_{3h} (transitivity + a transposition
  + primitivity via the 3-cycle structure), then fiber-product irreducibility for
  h != k follows from Goursat + non-isomorphic groups/degrees; for h=k the second
  component H_h corresponds to the diagonal quotient; S_{3h} 2-transitivity gives
  exactly 2 components: diagonal + one more = H_h irreducible. THIS closes (T3) and
  (T1) modulo: (i) monodromy = S_{3h} for all h, (ii) crit squarefreeness. Prove
  what you can of (i),(ii); verify (i) computationally for h<=6 by factoring C_h and
  checking Galois-generation certificates mod several primes (cycle types of Frobenius
  acting on roots of N_h(x)-T*q_h^3 specializations: sample T values, collect cycle
  types, certify S_{3h} via a p-cycle + transposition + primitive).
- If a full proof is out of reach, deliver: exact conditional statement
  ("if crit values squarefree and monodromy full for h<=X then ...") + computational
  certificates extending the census to h<=16, k<=16.

## Rules
No time estimates, no early stop. Everything claimed = proved or labeled VERIFIED-N
(computational) or CONDITIONAL (with exact hypothesis). ASCII.
