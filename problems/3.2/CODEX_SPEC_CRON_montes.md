# CODEX SPEC — Montes/p-adic factor-degree pipeline for M_h (Q6490 §6 experiment) + certificates to h=40

## Context
M_h in Z[Y]: symmetrized primitive gap polynomials of the Apery continuant. Construction (machine-verified in CRON_Mh_galois.py, THIS directory — read it): N_1=1, N_2=P(X+1), N_{d+1}=P(X+d)N_d-(X+d)^6 N_{d-1}, P(t)=34t^3+51t^2+27t+5; for even h remove the forced factor (2X+h+1); substitute X=(T-h-1)/2 (even in T), set Y=T^2, clear denominators and content -> M_h, deg m_h ~ 3h/2. Known: M_h irreducible over Q for h=2..14; Gal = S_{m_h} for h=2..11. Structure (Q6490, banked): N_h ≡ R_h^3 mod 3 (Frobenius cube); mod-2 four-cycle; lc(N_h) = Lucas U_h(34,1) with exact 17-adic valuation law.

## Task
1. **Certificates h=15..40**: for each h, construct M_h exactly (sympy, exact integer arithmetic; follow CRON_Mh_galois.py's M_of()), then:
   a. Attempt full irreducibility certificate: find a prime p0 < 50000 with M_h mod p0 irreducible (use sympy factor_list modulus=p0 on squarefree-good primes). Record (h, hash, p0) or NONE.
   b. If irreducible witness found, also collect Dedekind cycle-type certificates for S_n as in CRON_Mh_galois.py (prime-cycle q>n/2 for primitivity + single-even-part-2 type for transposition). Record certificate primes.
2. **Local factor-degree sets** for p in {2,3,5,17}, h=2..24: compute the exact factorization of M_h over Q_p — use sympy's padic tools if adequate, else implement via Hensel-lifted factorization of the squarefree part with sufficient precision, or shell out to PARI/GP `factorpadic` if gp is installed (check `which gp`; prefer gp — exact and fast). Output for each (p,h): list of local irreducible factor degrees D_{p,h}.
3. **Subset-sum intersection**: for each h, compute the set of achievable rational factor degrees = intersection over p in {2,3,5,17} of subset-sums of D_{p,h} (excluding 0 and m_h trivial). Report: for how many h does the intersection ALREADY prove irreducibility (only trivial degrees survive)? For which h does it fail and what degrees survive?
4. **Report** CODEX_MONTES_report.md: certificate table h=15..40, local degree sets table, intersection verdict per h, and the empirical answer to: "is the mixed-prime exclusion at {2,3,5,17} sufficient in practice?" — this is the load-bearing empirical question for the uniform Montes conjecture.

## Discipline
- New files only: CRON_montes_pipeline.py + CODEX_MONTES_report.md. Do NOT touch CRON_FRESH_EYES_pointwise.md, ERRATA.md, CAMPAIGN_MAP*, existing CRON_*.
- No ChatGPT bridge dispatch, no network beyond nothing (no network needed).
- No git commit — leave for session owner.
- Progress prints every ~30s; exact integer arithmetic only for polynomial construction.
- Sanity gates: reproduce M_h degrees 1,3,4,6,7,9,10,12,13,15 for h=2..11 and irreducibility h<=14 before proceeding (abort loudly if mismatch).
If PARI unavailable and sympy p-adic factorization proves unreliable, deliver parts 1+the h<=24 subset you can certify with exact Hensel code, plus a stall note.
