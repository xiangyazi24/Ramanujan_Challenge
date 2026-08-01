# CODEX SPEC — A_h certificates (PRIMITIVE-MORSE route, cheap descent)

Report to `CODEX_AH_CERT_report.md`. Construction (verify for h=2,3 symbolically
against direct computation first): C_h(X) = q_h N_h' - 3 q_h' N_h (deg 4h-4, the
critical-point polynomial in the BASE variable X); by the reflection symmetry C_h is
even in z = X + (h+1)/2 (verify!); define A_h(U) = primitive part of the z^2-descent
(A_h(z^2) = C_h up to content), deg A_h = 2h-2.
For 2 <= h <= 32 certify:
 (a) A_h irreducible over Q (direct factorization h<=12; modular certificates ell > 3h
     for larger h; degree-pattern intersection if needed). Cross-check h=2..7 against
     the parallel session's verified list (they report all irreducible).
 (b) [NONCOLLAPSE_h]: with m = 2h-2 and a1, a2 the first two power-sum/coefficient
     data of A_h (a1 = sum of roots, a2 = sum of squares of roots — from the first two
     coefficients), certify (m-1)*a1^2 - 2*m*a2 != 0. Exact integers.
 (c) Also record A_h(0) != 0 and squarefreeness.
Report: table h=2..32, methods, any failure flagged loudly (a reducible A_h or a
noncollapse failure = structural discovery: analyze it fully). Exact arithmetic, F_ell
computations for large h to avoid integer blowup. No early stop.
