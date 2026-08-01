# CODEX SPEC — attack [CRIT-2H] + extend certificates (follow-up to your CODEX_IRRED_THEOREM_report.md)

Report to `CODEX_CRIT2H_report.md`. Your previous report isolated:
[CRIT-2H]: s_h >= 2h-1 (number of simple nonzero critical values of delta_h) for all h>=2.
This is FAR weaker than irreducibility of C_h: total critical values (with mult) = 4h-4,
they come in +/- pairs (C_h(T)=C_h(-T), proved), so failure of CRIT-2H requires about
2h-3 coincidences/degenerations among 4h-4 values. Goal: prove or heavily corner it.

## Tasks
### A. Multiplicity cap route
[CRIT-2H] follows from: deg gcd(C_h, C_h') <= 2h-3 and ord_{T=0} C_h <= (something small).
1. Compute ord_{T=0} C_h exactly for h<=16 (critical points ON the zero fiber = common
   roots of A_h and N_h; you proved gcd(N_h, A_h) relates to gcd(N_h,N_h') — settle it:
   is C_h(0) != 0 ALWAYS (you verified BG(h) includes C_h(0)!=0 for h<=16 — find the
   structural reason: C_h(0)=0 iff N_h has a repeated root or a common root with q_h... 
   derive the exact identity C_h(0) = c * Res(N_h, A_h) and relate Res(N_h,A_h) to
   disc(N_h) * Res(N_h, q_h')-type products; then N_h squarefree + N_h(-j)!=0 (Apery
   positivity, banked) might PROVE C_h(0)!=0 for all h — attempt it).
2. A repeated nonzero critical value = two distinct critical points x1,x2 with
   delta_h(x1)=delta_h(x2), A_h(x1)=A_h(x2)=0 — i.e. a point of the SAME-GAP curve
   H_h(x1,x2) with both coordinates critical. Bound the number of such coincidences:
   Bezout on {H_h = 0} ∩ {A_h(x)=0} ∩ {A_h(y)=0} gives <= deg bounds; but we need
   the count of VALUE-coincidences <= 2h-3, and Bezout gives O(h^2) — too weak alone.
   Find the structural constraint: the +/- pairing (reflection) accounts for the pairs
   (a,-a); a coincidence WITHIN a pair means a = -a impossible (p odd, a nonzero);
   so coincidences are between DIFFERENT reflection-pairs. Use the reflection to
   quotient: critical points come in rho-pairs with opposite values; work on the
   quotient variable u = x(x+h+1)-type symmetric coordinate, halve the degrees, redo
   the count. Push to get the best provable cap on multiplicities.
### B. Certificate extension
Extend your two-good-prime exact certificates: verify BG(h) (all 4h-4 critical values
simple, nonzero) for 17 <= h <= 32, and gcd(C_h, C_k)=1 for all new pairs h<k<=32
(or a feasible subrange — state exactly what was certified). Same exact method as
before (modular certificates, not sampling).
### C. Asymptotic heuristic + failure hunting
If CRIT-2H were false for some h, C_h would have a root of multiplicity >= 2 at
scale ~h values. Search for the smallest structural candidate: compute gcd(C_h,C_h')
degree for h<=32 mod several primes (fast): report the observed multiplicity profile
(expected: all simple except forced). Also check the tau=0 column: multiplicity of
T=0 root if any.

## Rules
Exact arithmetic; label PROVED / VERIFIED-N / CONDITIONAL / OPEN precisely as before.
No early stop. ASCII.
