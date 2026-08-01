# CODEX SPEC: GARQI factorial-moment calibration in C (dyadic X up to ~5*10^5)

## Goal
Extend the Python calibration `CRON_growth_moments_verify.py` (X <= 8000, ledger appendix M.5) to dyadic X = 2^13 .. 2^19 with a fast threaded C program. This calibrates the AQI hypothesis (ledger `CRON_FRESH_EYES_pointwise.md` appendix M: quasi-independence of Apery zero rows) at scale.

## Definitions
Apery b_n: b_0=1, b_1=5, (n+1)^3 b_{n+1} = (34n^3+51n^2+27n+5) b_n - n^3 b_{n-1}.
Z_p = {0 <= r < p : b_r == 0 mod p} via recurrence mod p (precompute inv[1..p-1] in O(p); __int128 mulmod).
For dyadic X: I_X = (X, 2X], P_X = primes in (X/2, 2X], row A_{p,X} = (p + Z_p) ∩ I_X.
H(n) = #{p in P_X : n in A_{p,X}}  (equals #{p in (n/2,n] : p | b_n} by Gessel-Lucas; identity machine-verified in appendix M.2).
M_k(X) = sum_{n in I_X} H(n)(H(n)-1)...(H(n)-k+1),  S_X = sum_p |A_{p,X}|,  lambda_X = S_X/X.

## Outputs (CRON_garqi_moments_report.md)
For each X = 2^13, 2^14, ..., 2^19 (i.e. 8192 .. 524288):
1. S_X, lambda_X, lambda_X * log X   (first-intensity check; below 8000 this sits in [0.44, 0.78])
2. M_k(X) for k = 1..6 and the Poisson ratios M_k / (X * lambda_X^k)   (AQI says these stay O(1)-ish; below 8000: in [0.5, 2.1] for k=2,3)
3. max H(n) and its argmax n   (known: maxH = 3 at the 2*10^6 scale from earlier bigscan)
4. A short verdict line per X: do the ratios drift with X or stay flat?

## Sanity gates (MUST pass, abort loudly otherwise)
a. Reproduce the Python numbers at X = 4000 and X = 8000 exactly:
   X=4000: S_X=311, maxH=2, M_2=18; X=8000: S_X=622, maxH=3, M_2=54, M_3=6.
b. For X = 8192, cross-check H(n) for 20 random n in I_X against the direct definition: factor... no — directly test p | b_n for each prime p in (n/2, n] by computing b_{n-p} mod p via the recurrence (n-p < p) and Gessel b_n = b_1 * b_{n-p} = 5*b_{n-p} mod p.

## Performance
The dominant cost is Z_p for all p < 2^20 ~ 10^6 (~4*10^10 steps). Threads over primes (8 threads), accumulate rows into per-X hit arrays H[] with atomic or per-thread buffers merged at the end. Progress to stderr every ~5 s with ETA. Wall target ~15 min. Memory: H arrays are fine (sum of dyadic interval lengths < 2^20 ints).

## Discipline
- New files only, CRON_ prefix: CRON_garqi_moments.c, CRON_garqi_moments_report.md.
- Work only in problems/3.2/. Commit when done, report SHA + paste the final table.
- Do not touch CRON_FRESH_EYES_pointwise.md (single-writer: cron session owns it).
- Note: a sibling codex is independently scanning |Z_p| stats (CODEX_SPEC_CRON_zp_bigscan.md). Your Z_p kernel is an independent implementation — if you finish first, also print your own midpoint-zero primes (odd |Z_p ∩ [0,p)| ... i.e. primes with b_{(p-1)/2} ≡ 0) below 10^6 so the two scans can be diffed as a cross-check.
