# CODEX SPEC: |Z_p| big scan in C (p < 10^6, stretch 2*10^6)

## Goal
Extend the Python scan `CRON_zp_stats.py` (p < 30000, results in ledger appendix Q) to p < 10^6 with a fast threaded C program. This is the empirical backbone for [GAP-LT-MELLIN] (ledger `CRON_FRESH_EYES_pointwise.md` appendix Q): parity law, Poisson-pair law, midpoint (Apery-nonordinary) primes.

## Definitions
Apery numbers b_0=1, b_1=5, recurrence (n+1)^3 b_{n+1} = (34n^3+51n^2+27n+5) b_n - n^3 b_{n-1}.
For prime p >= 7: Z_p = { 0 <= r < p : b_r == 0 mod p }.
Compute the row by the recurrence mod p (all divisions invertible since n+1 < p). Precompute inv[1..p-1] in O(p) via inv[i] = -(p/i)*inv[p mod i] % p. Use 64-bit with 128-bit intermediate (__int128) mulmod, or Montgomery if you prefer.

## Outputs (single report file CRON_zp_bigscan_report.md + CSV)
1. Full distribution of |Z_p| over primes 7 <= p < 10^6 (histogram).
2. Parity check: list ALL primes with |Z_p| odd. THEOREM (proved, life-side Proposition): |Z_p| is odd iff p | b_{(p-1)/2} iff p is a nonordinary prime of the weight-4 level-8 form eta(2z)^4 eta(4z)^4. Known below 30000: exactly {11, 3137}. Report every new midpoint prime found (these are individually precious data points).
3. Record table: (p, |Z_p|) each time the max is broken. Below 30000 records are (7,0),(11,1),(17,2),(181,4),(379,6),(3727,8).
4. Running mean of |Z_p| at p = 10^5, 2*10^5, ..., 10^6 (checks mean ~ 1.02 stability / slow drift).
5. Poisson-pair fit: for even |Z_p|, fit |Z_p|/2 ~ Poisson(lambda), report lambda-hat, per-bin observed/predicted, chi^2. (Below 30000: lambda = 0.509, chi2 = 2.47 on df 3.)
6. Sanity gates (MUST pass before trusting the scan; abort loudly if any fails):
   a. Reproduce the p<30000 numbers exactly: mean = 1.0185 over 3242 primes, distribution {0:1933, 1:2, 2:1008, 4:260, 6:34, 8:5}, odd primes {11, 3137}.
   b. Spot check vs direct binomial definition: for p in {13, 101, 3137}, recompute b_{(p-1)/2} mod p as sum_k C(m,k)^2 C(m+k,k)^2 with factorial tables and compare with the recurrence value.

## Performance
Sum of p over primes < 10^6 is ~3.9*10^10 recurrence steps. Use pthreads or fork over prime ranges (8 threads on this M-series mini is fine); each prime is independent. Print progress (current p, primes done, ETA) to stderr every ~5 s. Wall target: under ~15 min at 10^6. If it comes in fast, extend to 2*10^6 and say so.

## Discipline
- New files only, CRON_ prefix: CRON_zp_bigscan.c, CRON_zp_bigscan_report.md, CRON_zp_bigscan.csv (per-prime p,|Z_p| — gzip if > 50 MB).
- Work only in problems/3.2/. Commit when done with a clear message. Report the commit SHA and paste the final summary block.
- Do not touch CRON_FRESH_EYES_pointwise.md (single-writer: cron session owns it).
