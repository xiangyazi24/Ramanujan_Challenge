# SPEC: [NO-RUN] full machine certificate + theorem writeup (laststand, codex lane 1)

Workdir: /Users/huangx/repos/Ramanujan_Challenge/problems/3.2 (work here; write all outputs here).

## Background (all verified this session; re-verify as you go)

Apery recurrence (n+1)^3 b_{n+1} = (34n^3+51n^2+27n+5) b_n - n^3 b_{n-1}; second solution c_n (c_0=0, c_1=6). p prime, N=p-2, window {1..N}. xi_r = (b_r : c_r) in P^1(F_p). Z_d = {r : xi_{r+d} = xi_r} = roots of N_d(r) mod p, deg N_d = 3(d-1), C_d = |Z_d| <= 3(d-1).

Companion matrix M_s = [[0,1],[-beta_s, alpha_s]], alpha_s = (34s^3+51s^2+27s+5)/(s+1)^3, beta_s = s^3/(s+1)^3.

PROVED SYMBOLICALLY (scratchpad frun_check.py, reproduce it): r, r+1 in Z_d and c_r != 0 forces xi_r = v(r,d) := (alpha_r - alpha_{r+d} : beta_r - beta_{r+d}) (cross-product collapse for companion matrices). Three consecutive r, r+1, r+2 in Z_d (generic branch) forces the orbit-free identity F(r,d) == 0 mod p where F = det(M_r v(r,d), v(r+1,d)) cleared of denominators; F is NOT identically zero, deg_r F = 8, deg_d F = 6, F = -24 d^2 * G(r,d) with G of total degree 10-ish. Empirics (p = 499, 997, 1999, 4999; d <= sqrt(p) log p): max run length in any Z_d is 2; zero 3-runs; 2-run counts 0/2/0/0 consistent with the ND/p^2 heuristic.

## Deliverables

1. `CODEX_NORUN_verify.py` — self-contained verifier:
   (a) recompute F(r,d) symbolically (sympy), assert F != 0, record its exact coefficients and degrees, factor it;
   (b) enumerate ALL branch patterns for a 3-run {r, r+1, r+2} subset Z_d: generic (all three c != 0), and every mixed pattern with c_{r+j} == 0 at exactly one j (two consecutive c-zeros are impossible: gap-1 Casoratian = 6/(r+1)^3 * unit != 0 — verify this constant claim from the recurrence, don't trust it). For each mixed pattern derive the corresponding orbit-free algebraic condition in (r,d) (eliminate the orbit point through the forced relations), assert it is a nonzero polynomial, record degrees;
   (c) degeneracy locus: v(r,d) = (0:0) requires alpha_r = alpha_{r+d} AND beta_r = beta_{r+d} mod p; show the beta equality factors through r(r+d+1) = z (r+d)(r+1) with z^3 = 1 and count solutions: prove the z=1 branch is empty for 1 <= d < p and the cyclotomic branches give O(1) solutions r per d; same analysis for alpha; output the exact exceptional count bound;
   (d) numerics at scale: for ~40 primes up to 10^5 (include some with (-51|p)=1 and =-1), compute all Z_d for d <= sqrt(p) log p, record: max run length, #2-runs, #3-runs, and CHECK every 2-run against the prediction xi_r == v(r,d) mod p (this validates the whole mechanism on the live orbit); also check every 3-run (expect none) against F(r,d) == 0 mod p;
   (e) print a final PASS/FAIL gate line.
2. `CODEX_NORUN_report.md` — theorem statement in inventory style:
   [NO-RUN] (statement with explicit constants): for every prime p > p_0 and every d, the number of r with {r,r+1,r+2} subset Z_d is at most (8 + mixed-branch constants) per d plus the O(1) degeneracy exceptions; hence total 3-run starts over d <= D is <= c * D, and sum_d (maxrun(Z_d) - 2)_+ <= c * D. Include full proof (the companion-matrix collapse, the branch enumeration, the degree certificates). Status expected: PROVED-all-h (it is uniform in d and p). Also state the corollary: the Strike2 padded-word interval-pairing clause is impossible for the Apery orbit (cite exactly which clause dies), and the honest limitation: scattered variable-gap padding is NOT excluded.
3. If any branch produces an identically-zero condition (rigidity failure), STOP and report loudly — that is a structural discovery, not a bug.

No effort cap. Do not stop because a sub-step is fiddly; enumerate all branches exactly. Verify every claimed constant by machine.
