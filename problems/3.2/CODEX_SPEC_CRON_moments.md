# SPEC: empirical growing factorial moments for H(n) (P3.2 pointwise campaign)

## Context
b_n = Apery zeta(3) numbers, recurrence (n+1)^3 u_{n+1} = (2n+1)(17n^2+17n+5)u_n - n^3 u_{n-1}, b_0=1, b_1=5.
Z_p := {r in [0,p) : p | b_r}; computed by the division-free recurrence on c_r = (r!)^3 b_r mod p:
c_{r+1} = (2r+1)(17r^2+17r+5)c_r - r^6 c_{r-1} mod p (c_0=1, c_1=5), zeros of c = zeros of b for r<p.
H(n) = #{p in (n/2,n] : n-p in Z_p}. See problems/3.2/CRON_fresh_scan.py for a working reference implementation (Python, N=30000).

## Task
1. Write a fast C program (single file, problems/3.2/CRON_moments_scan.c) that for all primes p < N computes Z_p and accumulates H(n) for n < N. Target N = 2,000,000 (measure runtime at N=200k first, scale; stay within ~30 min on an M-class Mac; N is adjustable).
2. Compute and report, for n in [10^4, N): the factorial moments M_k = (1/count) * sum_n (H(n))_k (falling factorial) for k = 1..6; the histogram of H; max H(n) and the top-20 record (n, H(n)) list.
3. Poisson test: under H ~ Poisson(lambda_n) with lambda_n ~ ln2/ln n, predict M_k ≈ mean over n of lambda_n^k. Report the ratio M_k / predicted for each k. Growing-moment health = ratios staying O(1) as k grows.
4. Record structure probe: for each record n with H(n) >= 4 (if any), print the bad primes, residues z = n-p, z/p ratios, p mod 24, and gcd(p-1, n-1) (character order relevance: order = (p-1)/gcd(p-1,n-1) — LOW order would be a structure signal; expect all high-order).
5. Write results to problems/3.2/CRON_moments_report.md (tables + 5-line interpretation: does the empirical data support the k ~ log X Poisson-moment hypothesis of the strategic program?). Verify the C scan against CRON_fresh_scan.py output at N=30000 before the big run (max H and histogram must match exactly).

## Discipline
- New files only (CRON_ prefix); do not modify existing files.
- git add + commit the new files with a clear message when done.
- Machine-verify everything; no unverified claims in the report.
