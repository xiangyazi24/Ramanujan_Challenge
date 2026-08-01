# CODEX SPEC: G_h absolute irreducibility certificate chain h<=40

## Mission

[L2-FREQ] (banked theorem, conditional) needs: the collision curve G_h(x,y) absolutely irreducible over F_p (generic p) for all h. Symbolic status: irreducible over Q for h=2..6 only. Extend the certificate chain to h<=40 with machine-verified certificates, two independent methods per h where feasible. This parallels the Montes campaign for M_h (`CODEX_MONTES_report.md` in this directory — read its structure for the certificate-chain format).

## Definitions

```
K_0(X)=1, K_1(X)=P(X), P(X)=34X^3+51X^2+27X+5
K_{m+1}(X)=P(X+m)K_m(X)-(X+m)^6 K_{m-1}(X)
N_h(r)=K_{h-1}(r+1),  deg N_h = 3(h-1)
R_h(x,y) = N_h(x)*prod_{j=1..h}(y+j)^3 - N_h(y)*prod_{j=1..h}(x+j)^3
G_h(x,y) = R_h(x,y)/(x-y)   (exact division, verify remainder 0)
```
(Convention check: the banked h=2..6 result is "R_h = (x-y)*G_h with G_h irreducible over Q, single factor via sympy factor_list". First REPRODUCE h=2..6 exactly; if your G_h differs by the cofactor convention — e.g. the even-h rational-root factor of N_h — resolve against the banked convention in the ledger `CRON_FRESH_EYES_pointwise.md` appendix AS.4 and record which convention you certify.)

## Certification targets, per h = 2..40

1. Q-irreducibility of G_h (sympy/flint factorization over Q — this alone may get slow for large h; switch to modular methods early).
2. ABSOLUTE irreducibility (over Qbar, equivalently: geometrically irreducible as a plane curve). Certificate methods (pick >=2 independent per h where runtime allows, >=1 always):
   - (a) Find a prime q (not dividing leading/content data) with G_h mod q absolutely irreducible over F_q. Certify absolute irreducibility over F_q by: G_h irreducible in F_q[x,y] AND stays irreducible in F_{q^k}[x,y] for all k <= deg bound that could split it (it suffices to check k | number-of-geometric-components possibilities; a geometrically reducible but F_q-irreducible curve becomes reducible over the splitting degree). Then absolute irreducibility over Q follows for all but finitely many p (record q).
   - (b) Point-count consistency: for several extensions F_{q^k}, count points of G_h=0 and compare with q^k + O(d^2 q^{k/2}) for a single geometric component vs c*q^k for c>=2 components. This is probabilistic evidence — label it EVIDENCE, not certificate, unless combined with (a).
   - (c) Any exact criterion you can implement soundly (Noether irreducibility via a certified irreducible specialization: if G_h(x, y0) stays... — careful, specialization criteria need the right hypotheses; only use criteria you can state and verify precisely; write the criterion used into the report).
3. Hash each certificate (sha256 of a canonical text) as in the Montes report.

## Anomaly protocol

If some h fails (reducible, or no certificate found): do NOT stop the chain; record the h with full diagnostic (the factorization if reducible — a genuine factorization at some h would be a MAJOR finding, double-verify it exactly over Z before reporting) and continue.

## Deliverables

- `CRON_ghirred.py` (progress prints per h).
- `CODEX_GHIRRED_report.md`: per-h certificate table (method(s), q used, hash, wall time), the h=2..6 reproduction check, any anomalies.
- `ghirred_certs.json`.

## Rules

- python-flint (or PARI via cypari2 if available) strongly preferred over plain sympy for bivariate factorization mod q — check what's installed, report what you used.
- Exact arithmetic only for any claimed factorization.
- No effort cap. If h=40 is out of runtime reach, deliver the maximal contiguous chain with per-h timings and a feasibility note for the rest.
