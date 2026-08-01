# CODEX SPEC — transfer-cocycle codegree protocol + annealed-vs-ordered (max)

## Mission

Execute the numerical protocol for the transfer-cocycle phase-modulus program
(campaign-2 prize line). Three deliverable experiments, one script, one report.
This is numerics + exact arithmetic; no theory required beyond what is stated.

## The exact objects (verified conventions — use as given)

- Apéry recurrence (n+1)^3 b_{n+1} = P(n) b_n - n^3 b_{n-1},
  P(n) = (2n+1)(17n^2+17n+5) = 34n^3+51n^2+27n+5, b_0=1, b_1=5.
- Transfer maps on P^1(F_p): F_{p,u}(x) = (u+1)^6 / (P(u) - x),
  for u in F_p \ {-1}. Each F_{p,u} is a Moebius permutation of P^1(F_p).
  Convention: F_{p,u}(infinity) = 0; if P(u) - x = 0 then F_{p,u}(x) = infinity.
- Ordered (true Apéry) orbit: x_0 = 0, x_{n+1} = F_{p,n}(x_n) for n = 0..p-2.
  GATE (must pass before anything else, for every p used): the visit set
  {n : x_n = infinity, 1 <= n <= p-2} must equal
  {n : b_n ≡ 0 mod p} computed directly from the recurrence mod p.
  Also check x_1 = 1/5 mod p.
- Annealed model: X_{j+1} = F_{p,U_j}(X_j), U_j i.i.d. uniform in F_p\{-1},
  X_0 = 0. Target state: infinity.
- Codegree: C_p(x,x') = #{(u,v) in (F_p\{-1})^2 : F_{p,u}(x) = F_{p,v}(x')}.

## Experiment 1 — codegree exceptional-locus scan

For p in {101, 211, 401, 601, 1009} (drop the largest if too slow):
- Compute C_p(x,x') for ALL pairs (x,x') in P^1(F_p)^2 if feasible at that p
  (cost ~ p^2 workspace using the bucket trick: for fixed x, the multiset
  {F_{p,u}(x)}_u; then C_p(x,x') = sum over points y of
  mult_x(y)*mult_{x'}(y) — this is O(p^2) total per p, NOT p^3; do it).
- Flag all pairs with |C_p(x,x') - p| > 6*sqrt(p).
- Fit flagged pairs to low-degree algebraic relations in (x,x'): print the
  flagged set for each p; test membership in candidate relations: x = x',
  x = reflection-image of x' (derive: the Apéry reflection b_{p-1-n} = b_n
  suggests testing the involution induced on x-space — try x*x' = c,
  x + x' = c, and bilinear ax x' + bx + cx' + d = 0 fits via exact linear
  algebra over F_p on 4+ flagged points), and report which relations are
  STABLE ACROSS ALL p TESTED (a genuine char-0 correspondence) vs sporadic.
- Report: for each stable relation, the observed C_p size on it
  (e.g. ~2p? ~p+c*sqrt p?).

## Experiment 2 — annealed vs ordered visit statistics

For ~200 primes p in [10^3, 2*10^4] (all primes in a sub-range is fine):
- Ordered: compute V_p = #visits to infinity of the true orbit (= |Z_p|),
  AND the visit positions (to observe the reflection pairing).
- Annealed: for each p run ONE simulated annealed trajectory of length p-1
  (same length as ordered), count visits. (Use a fixed RNG seed recorded in
  the report for reproducibility.)
- Aggregate across primes: empirical distribution of V_p (ordered) vs
  annealed counts. Tables: frequency of 0,1,2,...; mean; variance;
  factorial moments (m)_2, (m)_3 averages. Compare both to Poisson:
  the ordered counts pair up under reflection (expect even counts mostly,
  orbit-count ~ Poisson(1/2) after halving non-central visits); the annealed
  counts should be ~ Poisson((p-1)/(p+1)) ≈ Poisson(1) WITHOUT any pairing
  structure. Verify specifically: (a) ordered counts have an excess of even
  values (reflection signature), annealed do not; (b) after halving
  (V_p - central)/2 the ordered distribution matches Poisson(1/2) moments;
  (c) annealed matches Poisson(1) moments.
- Also record for ordered orbits: distribution of gap between paired zeros
  r and p-1-r (sanity: exact reflection) — i.e. verify machine-exactly that
  the visit set is reflection-symmetric: n in visits iff p-1-n in visits
  (for 1 <= n <= p-2). Any asymmetric visit is a CRITICAL finding
  (would contradict our proved reflection law) — flag loudly.

## Experiment 3 — post-visit restart structure

F_{p,u}(infinity) = 0 for every u, so the chain restarts at 0 after each
visit. For the ordered orbit: collect the empirical distribution of return
times (distance from a zero at n to the next zero), across all primes of
Experiment 2. Compare with the geometric/exponential prediction from the
annealed model. Report the smallest observed return time and its frequency
(short returns are the danger zone for Poisson approximation).

## Deliverables

- Script: research/scripts/q32_codegree_protocol.py (stdlib only; numpy
  allowed if available but degrade gracefully). Prints GATE VERIFIED per
  prime (Experiment-2 gate on a subsample of >= 20 primes and all
  Experiment-1 primes), then results.
- Report: CODEX_CODEGREE.md (problems/3.2 root): tables for all three
  experiments, the stable exceptional relations found (or "none beyond
  diagonal"), the annealed-vs-ordered verdict (do the two ensembles agree
  at the resolution tested — this measures the size of the quenched-order
  gap empirically), and a LIMITATIONS section.

## Hard constraints

- Do NOT touch: proof.tex, ERRATA.md, CAMPAIGN_MAP*, lean/, any existing
  script or report. New files only (the two above).
- No external LLM/bridge dispatch. Exact arithmetic for all F_p work.
- If the Experiment-2 gate fails for ANY prime, stop and write the stall
  report with the witness (that would mean our cocycle conventions are
  wrong, which is critical information).

## Acceptance

python3 research/scripts/q32_codegree_protocol.py exits 0 and prints all
GATE VERIFIED lines; report tables match script output.
