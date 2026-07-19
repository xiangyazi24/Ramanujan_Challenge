# CODEX SPEC W2: separated-block resultants and the mesoscopic root sum

## Context and prize chain

For the Apéry gap polynomials N_h mod p (notation of problems/3.2/proof.tex),
let r_p(h) = #roots of N_h in F_p and R_p(H) = sum_{h=2}^H r_p(h).
The chain (all reductions already proved in proof.tex, rem:incidence and
rem:collision):
  E_p(H) := #{(x,d,r): d+r <= H, N_d(x) = N_r(x+d) = 0 mod p} << H^{3/2}
    ==> (via Kovari-Sos-Turan on the bipartite gap graph) R_p(H) << H^{3/2}...
    more precisely a fiber + codegree bound on the graph gives R_p(sqrt p)
    = O(sqrt p)
    ==> Z(p) << p^{1/2}   [rem:incidence: R_p(H) << H^beta with beta=1].
PRIZE: Z(p) << sqrt(p) would (a) improve every averaged theorem in the paper,
(b) lower the high-moment requirement in thm:hm-pointwise from k > 6 to
k > 4 (exponent 1/2 + 2/k < 1), a large strategic gain.
Computed data: E_p(H)/H^{3/2} -> 0 for all tested p (table in rem:collision);
R_p(sqrt p)/sqrt p ~ 1.5 (Chebotarev prediction 3/2).

## Why now: a technique transfer

The centered-coefficient/Pell method (problems/3.2/nv_theorem.tex) just
resolved the bordered-certificate nonvanishing that endpoint methods could
not reach. The analogous objects here are the SEPARATED-BLOCK RESULTANTS
  S_{d,r} = Res_x(N_d(x), N_r(x+d)) in Z,
which control triple orbit-collisions: p | S_{d,r} iff N_d and N_r(.+d)
share a root mod p iff there is x with orbit collisions at {x, x+d, x+d+r}.
Known exact structure (proof.tex, rem:adj-res, prop:bezout, rem:sep-res):
- |Res(N_h, N_{h+1})| = prod_{j=1}^{h-1} (j!^3 b_j)^6  (adjacent case).
- |Res(N_d, N_e)| = prod_{j=1}^{d-1}(j!^3 b_j)^6 * |R_{d,e-d}| where
  R_{d,r} = Res_x(N_d(x), N_r(x+d))  (so S == R here).
- |R_{d,r}| = |R_{r,d}| (reflection). R_{d,1} = 1.
- R_{2,2} = 2^5 * 5^4 * 11^2 * 17^3 * 71; at p=71 the shared root is
  x = -5/2 (nonboundary).
- For d = 2: R_{2,r} = 2^{3(r-1)} G_r(-1/2) * 17^{3(r-1)} N_{Q(theta)/Q}(G_r(theta))
  where G_r = N_r(x+1), Q(y) = 17y^2+17y+5, theta a root of Q.

## Target

Find enough structure in the family {R_{d,r} : d, r >= 2} to bound, for a
fixed prime p, the number of pairs
  V_p(H) = #{(d,r): d + r <= H, p | R_{d,r}}.
Ranked goals:
  (G1) V_p(H) << H (linear) for H <= sqrt p  ==> combined with the per-pair
       root bound deg gcd <= 3(min(d,r)-1) and a Kovari-Sos-Turan argument,
       derive E_p(H) << H^{3/2} and complete the chain to Z(p) << p^{1/2}
       (write the full derivation — the KST step needs care, do it honestly).
  (G2) V_p(H) << H^{2-delta} for some delta > 0, with the corresponding
       weaker endpoint of the chain (work out what it gives for Z(p)).
  (G3) A closed form or congruence-level structure theorem for R_{d,r}
       (e.g., R_{d,r} as a product of Pell/Apéry-type quantities, or exact
       formulas for small d generalizing the d=2 norm decomposition), even
       without the counting corollary.
  (G4) Sharp obstruction: why the family resists (with the exact missing
       lemma named).

## Suggested attacks

1. CENTERED/PELL TRANSFER: compute the leading and first few centered
   coefficients of the resultant... more precisely, R_{d,r} is an integer,
   not a polynomial — instead study, for p | R_{d,r}, the certificate
   x_0 in F_p with N_d(x_0) = N_r(x_0+d) = 0 and run the centered-coefficient
   elimination on the PAIR of conditions as in nv_theorem.tex: what
   Pell-algebra constraints on (l_d, l_r, ...) does a shared root force?
   Note the danger-case discipline: identities may hold only modulo the
   ideal of leading coefficients — state validity domains precisely.
2. RECURRENCE ON RESULTANTS: from N_{r+1}(x+d) = P(x+d+r)N_r(x+d) -
   (x+d+r)^6 N_{r-1}(x+d), derive a recurrence/divisibility lattice for
   R_{d,r} in r (resultant multiplicativity along the recurrence — the
   adjacent-resultant recursion Res(N_h,N_{h+1}) = N_h(-h)^6 Res(N_{h-1},N_h)
   in rem:adj-res is the d=1 case; generalize).
3. RESTART/BEZOUT: lem:restart and prop:bezout convert shared roots into
   shifted lower-order conditions; iterate to a canonical minimal pair.
4. p-DIVISIBILITY VIA APERY VALUES: adjacent resultants are products of
   b_j's — p | Res(N_h, N_{h+1}) iff p | b_j for some j <= h-1. If R_{d,r}
   admits ANY product formula in terms of b-values and Pell numbers, then
   V_p(H) is controlled by zero/value counts of those sequences — quantities
   we already bound (Z(p) << p^{2/3}, fiber << p^{3/4}).

## Computational duty (FIRST)

- Compute R_{d,r} exactly (integer) for all 2 <= d, r with d + r <= 24.
  Factor them. Look for: product structure, Pell factors (l_m values),
  Apéry factors (b_j), growth of the primitive part, repetition patterns.
- For primes p <= 5000: compute V_p(H) for H = 24 and the E_p contribution;
  test G1's V_p(H) << H against data.
- Bank results in problems/3.2/meso_resultants.md + meso_explore.py.

## Deliverables

1. problems/3.2/meso_result.tex — theorem/partial + proof (or G4 obstruction),
   notation matching proof.tex/nv_theorem.tex.
2. meso_explore.py + meso_resultants.md.
3. Verification script meso_verify.py for any new identity (PASS/FAIL).
4. STALL REPORT convention if all goals fail.

## Hard constraints

No numerics-as-proof; no edits to existing files; symbolic cross-checks for
every identity; honest validity domains (learn from the C4 lesson: an
identity used in a special case must be stated for that case).
