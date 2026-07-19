# CODEX SPEC: large-scale B(n) computation for the Apéry GCD conjecture

## Goal

Extend the computational verification of the top-half bad-prime count from
n <= 200,000 to n <= 1,000,000 (stretch: 2,000,000), plus windowed
statistics. Existing verified fact: B(n) <= 3 for all n <= 200,000.

## Definitions

- Apéry numbers mod p: b_0 = 1, b_1 = 5,
  (m+1)^3 b_{m+1} = (34m^3+51m^2+27m+5) b_m - m^3 b_{m-1}  (mod p).
  Computing b_j mod p for j = 0..p-1 costs O(p) mulmods with one modular
  inverse per step — avoid inversions by the cleared form: iterate
  Y_{m+1} = P(r+m) Y_m - (r+m)^6 Y_{m-1} is NOT needed here; simplest is
  to precompute inverses of (m+1)^3 via a running product / batch
  inversion, or iterate the pair (u_m, v_m) = (m!^3 b_m, m!^3 b_{m-1}·(m)^3...)
  — your choice, just make it O(p) mulmods with O(1) inversions total
  (batch/Montgomery trick), NOT O(p) inversions.
- Zero set: Z_p = { r in [0, p-1] : b_r ≡ 0 mod p }.
- Top-half bad-prime count: for n, B(n) = #{ p prime, n/2 < p <= n : (n-p) in Z_p }.
  Note r = n - p < p automatically, and (n-p) in Z_p iff p | b_{n-p}.

## Algorithm (must be this shape — it is the cheap one)

1. Sieve primes up to P_MAX = 2,000,000 (or 1,000,000 if runtime demands).
2. For each prime p in [7, P_MAX]: compute b_j mod p for j < p; record the
   positions r in Z_p (there are ~1 per prime on average; store as
   (p, r) pairs). Total cost sum_p O(p) ≈ 1.4e11 mulmods at P_MAX = 2e6 —
   use uint64 with __int128 mulmod or Montgomery multiplication; estimate
   and report expected runtime; parallelize with pthreads/OpenMP across
   primes (the machine has multiple cores — detect and use them).
   ALSO record per-prime Z(p) count for the histogram (cross-check against
   known: pair count Z(p)/2 ~ Poisson(1/2), max Z(p) = 12 for p <= 1e6).
3. Post-process: for every recorded (p, r), the integer n = p + r (which
   lies in [p, 2p)) gets one bad-prime hit — exactly the condition
   p in (n/2, n], n - p = r in Z_p. So B(n) for ALL n <= P_MAX is obtained
   by a single scatter: B[p + r]++ over all recorded pairs (only keep
   n <= N_MAX = P_MAX).
4. Report:
   - max_n B(n) for n <= N_MAX and the argmax list (all n with B(n) >= 4
     if any exist — THIS IS THE HEADLINE: a single n with B(n) large would
     be evidence against the conjecture's mechanism; B <= 3 everywhere is
     strong Poisson-consistency).
   - Histogram of B(n) vs Poisson(log2/log n) prediction (in dyadic blocks).
   - Weighted sums: W_top(n) = sum over bad p of log p — report
     max W_top(n)/n over dyadic blocks (should decay ~ log n/n scale).
   - Windowed stats: for each dyadic scale N in {2^17..N_MAX}, split
     (N, 2N] into windows of length N/64; report the max over windows of
     (window sum of B), and the window variance vs the independent-model
     prediction (window mean). This tests LOCALIZED dispersion — the key
     open hypothesis of the paper.
   - Var(B)/E(B) overall per dyadic block (Poisson signature: should be ~1).

## Deliverables

1. `problems/3.2/bn_bigscan.c` — the scanner (self-contained, compiles with
   `cc -O2 -o bn_bigscan bn_bigscan.c -lpthread` or with OpenMP flags).
2. `problems/3.2/bn_bigscan_report.md` — the report with all numbers above.
3. Raw (p, r) pairs written to `problems/3.2/data_zp_pairs.bin` (binary,
   two uint32 per record) so later sessions can post-process without
   recomputation; document the format in the report.

## Constraints

- Verify correctness first on small range: recompute B(n) for n <= 200,000
  and CONFIRM max B = 3 and the known histogram (93.7% zeros) before
  launching the full run. If mismatch, stop and debug — do not burn hours
  on a wrong kernel.
- Run the full computation yourself (background OK), monitor, and write the
  report from actual output. If P_MAX = 2e6 would exceed ~3 hours on this
  machine, drop to 1e6 and say so in the report.
- Do not modify proof.tex or other files.
