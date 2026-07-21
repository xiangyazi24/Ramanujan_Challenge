# CODEX SPEC — Lemma B′: eigenvalue identification (the crystalline kill-gate)

ASSUME the computation is feasible; your job is to run it and report a BINARY
verdict G1 or G2. This TERMINATES the crystalline question one way or the other.

## Why this run (Fable oracle adjudication — the route)

The Lemma-A Jacobi-Mellin form was circular: its h≈p "unit terms" are literally
the Apéry summands C(j,k)²C(j+k,k)² (Kummer's theorem in disguise), so "unit-term
cancellation = b_j≡0". That was the WRONG lens. The RIGHT object is the ≤B_Λ
Frobenius eigenvalues of the fixed toric cohomology H*_c(U, L_{p,j}), U =
(G_m)³∖{u=1,v=1,w=1,uvw=1}, L_{p,j} the rank-1 Kummer sheaf with trace
u ↦ Λ(u)^j (Λ = the fixed Laurent polynomial (1+x)(1+y)(1+z)((1+y)(1+z)+xyz)/(xyz),
CT Λ^n = b_n).

FACTS (Fable, endorse):
- (a) Σ(slope-0 eigenvalues, with signs) ≡ −b_j (mod 𝔭). So if the Newton slope-0
  multiplicity is exactly 1, then p∤b_j. Hence Z_p ⊆ {j : Newton slope-0 count ≠ 1}
  is a THEOREM. No computation needed for the inclusion.
- (b) The HODGE polygon of the χ_j-twist is PURE POLYTOPE COMBINATORICS
  (Adolphson–Sperber): j enters only through fractional parts {j·w/(p−1)} weighted
  over the Newton polytope of Λ. Closed form, no b_j, no length-h sum. This is the
  non-circular computability.
- (c) THE RAZOR: if the Hodge slope-0 multiplicity h⁰(j)=1 for generic j, then
  "j∈Z_p ⟺ j-twist non-ordinary" is a Hasse-invariant TAUTOLOGY (like "p|a_p ⟺ E
  supersingular" restates rather than solves) — and the character line has no
  moduli Hasse polynomial. So the route is a genuine reduction IFF h⁰(j)≥2 with
  identifiable eigenvalues; if h⁰≡1 it is formally dead.

WHY h⁰≥2 IS PLAUSIBLE (the structural bet): Apéry ζ(3) carries the Sym² structure
and the factor Δ(t)=t²−34t+1, whose roots are 17±12√2 = (1±√2)⁴ — the global units
governing b_n ~ (17+12√2)^n n^{−3/2}. The slope-0 part of the j-twist MAY split as
  {μ_j : moving unit root} ⊕ {explicit Kummer eigenvalue ω(17±12√2)^{±j}}.
Then p|b_j ⟺ collision μ_j ≡ ±ω(17+12√2)^{±j} (mod 𝔭) — a moving unit root hitting
an explicitly parametrized unit. THAT is the transform we want (new attack surface:
p-adic interpolation of μ_j, height/fixed-object arguments against the explicit
branch). Testing whether this split is real is the whole point of this run.

## The computation — three parts, one run

### Part 1 — Hodge polygon (pure combinatorics, closed form, NO numerics)
Compute the Adolphson–Sperber Hodge polygon of the χ_j-twisted exponential sum
over U from the Newton polytope of Λ (in the u,v,w coordinates: the fixed
4-hyperplane-arrangement complement, divisor {u=1,v=1,w=1,uvw=1}). Output the
slope-0 Hodge multiplicity h⁰(j) as an EXPLICIT function of the fractional-part
data ⟨j·w/(p−1)⟩ over the relevant polytope vertices/faces. Deliverable: the
closed-form h⁰(j), and whether it is generically 1 or generically ≥2.
(If the Adolphson–Sperber setup for this specific arrangement is delicate, fall
back to computing h⁰(j) as the number of slope-0 slots via the toric
Hodge-theoretic recipe and VERIFY it against the Part-2 Newton data at small p.)

### Part 2 — Newton polygon (numerics, small p, valuations only)
For p≤60 all j; and p≤300 for all j∈Z_p plus 20 random j per p: compute the
𝔭-adic Newton polygon of the L-factor of H*_c(U,L_{p,j}) from extension-field
traces S_j^(r) = −Σ_{u∈U(F_{p^r})} χ_j(Λ(u)) for r=1..deg. Use the VERIFIED (A1)
Jacobi identity base-changed to F_{p^r} (Hasse–Davenport lifts the Jacobi sums;
evaluate 𝔭-adic valuations via Gross–Koblitz / p-adic Γ at low precision — you
need VALUATIONS, not full values). Cost O(p^{r+1}) per (j,r), feasible. Output the
Newton slope-0 multiplicity per (p,j) and confirm Σ(slope-0) ≡ −b_j mod 𝔭.

### Part 3 — eigenvalue matching (the key deliverable)
For generic j, test each slope-0 eigenvalue (as an exact 𝔭-adic / algebraic number,
not just its valuation) against the candidate list:
  ω(17+12√2)^{±j}, ω(17−12√2)^{±j}, Jacobi-sum monomials affine in j, perfect
  squares μ² (the Sym² signature).
Report EXACT 𝔭-adic matches (Teichmüller of 17±12√2 mod p requires √2, i.e.
p≡±1 mod 8 — handle both √2-split and √2-inert primes; the eigenvalue may live in
F_{p²} in the inert case). Determine whether an explicit branch is present.

## BINARY GATES (the whole point — report exactly one)
- **G1 (KILL):** h⁰(j)=1 generically AND Newton=Hodge generically (no split, no
  explicit branch). ⟹ crystalline route is tautological. Declare it dead; do the
  Part-1 closed form + one confirming table and STOP. (The program then pivots per
  the pre-authorized π1/π2/π3 below — NOT your job this run, just flag G1.)
- **G2 (TRANSFORM):** h⁰(j)≥2 with at least one EXPLICITLY IDENTIFIED slope-0
  eigenvalue (e.g. ω(17±12√2)^{±j}). ⟹ report the collision equation
  μ_j ≡ (explicit target) mod 𝔭 and the exact form of the explicit branch. This
  opens the next lemma (unit-root anti-concentration).

## Deliverables
- problems/3.2/lemmaBprime_result.tex — Part-1 closed-form h⁰(j) with derivation,
  the inclusion theorem Z_p⊆{Newton slope-0 count≠1}, the Part-2/3 data tables,
  and the BINARY G1/G2 verdict with the explicit-branch form if G2. Cite the
  Λ marked coordinate (thm:oracleC-marked-coordinate), the (A1) Jacobi skeleton
  (thm:lemmaA-skeleton, reused as the F_{p^r} engine), Δ(t)=t²−34t+1.
- problems/3.2/lemmaBprime_explore.py — Parts 1–3 computation.
- problems/3.2/lemmaBprime_verify.py — PASS/FAIL: Σ(slope-0)≡−b_j mod p; h⁰ formula
  matches Newton data at small p; any claimed eigenvalue match verified exactly.

## Hard rules / negative acceptance
- The inclusion Z_p⊆{slope-0≠1} must be stated as the theorem it is (from Fact (a)).
- Do NOT reuse the Jacobi-Mellin unit count (that was the circular lens); the
  slope count here is the Frobenius/Newton multiplicity, a DIFFERENT object.
- A G2 claim REQUIRES an exact 𝔭-adic eigenvalue match, not a valuation
  coincidence or an "on average" statement.
- A G1 verdict is a FIRST-CLASS, valuable result — do not manufacture a fake G2.
- No numerics-as-proof for Part 1 (Hodge polygon must be derived); Parts 2–3 are
  explicitly empirical, label them so.
- Do not modify existing files. Run lemmaBprime_verify.py; include its tail.
- If Part 1's Adolphson–Sperber derivation stalls, deliver Parts 2–3 (which
  determine the gate empirically) + a precise stall note on the closed form.
