# CODEX SPEC — W_h irreducibility certificates (component tower closure)

Report to `CODEX_WH_CERT_report.md`. Definition chain (from cron, verify each step
symbolically for h=2,3 before scaling): q_h(X)=prod_{a=1}^h (X+a); D_h=q_h^3;
C_h = q_h*N_h' - 3*q_h'*N_h (deg 4h-4); V_h(T) = Res_X(C_h, N_h - T*D_h);
W_h(U) := primitive part of the T^2-descent of V_h (V_h is even in T by the mirror
symmetry; W_h(T^2) = V_h(T) up to content), deg W_h = 2h-2.
N_h from the standard recurrence (campaign3_questions/CTX.txt).

TASKS:
1. Verify the definition chain and evenness of V_h exactly for h=2..6 (sympy).
2. Certify, for every 2 <= h <= 32:
   (a) W_h(0) != 0 (equivalently Res(C_h, N_h) != 0);
   (b) W_h is irreducible over Q.
   For h <= 12 direct sympy factorization is fine. For h > 12 use MODULAR certificates:
   irreducibility mod ell for a prime ell > 3h of good reduction (W_h irreducible over
   F_ell => irreducible over Q); if W_h is reducible mod every tested ell (test >= 8
   primes), use degree-pattern intersection certificates (gcd of the sets of possible
   factor-degree sums across primes = {0, 2h-2} proves irreducibility). Record for each
   h: the certifying prime(s) and method.
3. Also record squarefreeness of W_h (gcd(W_h, W_h') = 1) for each h (implied by
   irreducibility but check independently where cheap).
4. Report: table h=2..32 with verdicts; any FAILURES prominently flagged (a reducible
   W_h would be a major structural discovery — if found, factor it completely and
   report the factor degrees and any pattern).
Exact arithmetic; efficient (resultants mod ell computed directly in F_ell to avoid
huge integers for large h). No early stop.
