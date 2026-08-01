# CODEX SPEC — [GAP-MARKED-DEC] numerical diagnostics (high)

## Mission

Falsification/diagnostic protocol for the marked two-gap decoupling lemma
(the next rung above the short-gap lemma on the energy ladder). Background:
b_r = Apéry zeta(3) numbers mod p ((n+1)^3 b_{n+1} = P(n)b_n - n^3 b_{n-1},
P(n) = 34n^3+51n^2+27n+5, b_0=1, b_1=5). Collisions b_s = b_{s+h} = b_{s+k}
(0 < h < k) force Delta_{h,k}(s) = 0 where Delta is the two-row determinant
of the state-free collision pair (an explicit polynomial of degree
3(h+k)+O(1); leading coefficient -c_{k-h-1} with c the Chebyshev-like
sequence c_{m+1} = 34c_m - c_{m-1}, c_0=1, c_1=34). The candidate lemma
[GAP-MARKED-DEC] asserts the number of marked triples
T_p(H) = #{(s,h,k): 0<h<k<=H, b_s=b_{s+h}=b_{s+k}, s,s+k <= p-2}
is O_eps(H^2 p^eps) uniformly for H <= sqrt(p). This would extend the O(p)
short-gap collision range from H <= p^{1/3} to H <= p^{1/2}.

## Experiments

For primes p in a sample of ~300 spread over [10^3, 10^6] (log-uniform;
include the known worst primes from CODEX_TE_SCAN.md spike table: 1069,
1193, 1223, 1231, 1499) and H = floor(sqrt(p)):

1. Compute b_r mod p (one pass). Build value -> sorted position list.
2. Enumerate marked triples with both gaps <= H: for each value class with
   >= 3 positions, for each position s in the class, count pairs of later
   positions within distance H (two-pointer; cost O(class^2) worst but
   classes are tiny by the energy scan).
   Record T_p(H) total, and T_p^refl(H) = triples containing at least one
   forced reflection pair (r, p-1-r); report both and the difference.
3. Normalize: (T_p - T_p^refl)/H^2. Aggregate by dyadic p ranges: mean,
   max, argmax. Danger signal: sustained growth ~ p^delta of the dyadic
   maxima (fit log-log slope; the lemma predicts slope ~ 0 up to p^eps).
4. Phantom ratio (the sharper diagnostic from the source design): for a
   SUBSAMPLE of ~20 primes across the range and ~50 random (h,k) pairs
   with k <= H each: count the actual roots of Delta_{h,k} mod p —
   construct Delta_{h,k}(s) directly: using the block continuant formulas
   U_m (U_{-1}=0, U_0=1, U_{m+1}(s) = P(s+m)U_m(s) - (s+m)^6 U_{m-1}(s),
   polynomials in s mod p), the two rows are
   R_h = (s^3 U_h(s) - (s+h)^3 D_h(s), -s^9 U_{h-1}(s+1)) with
   D_h(s) = prod_{j=0}^{h-1}(s+j)^3 — VERIFY this normalization by
   checking at 20 random s that R_h · (s^3 b_s, b_{s-1}) = 0 iff
   b_{s+h} = b_s (derive/adjust the exact row convention from the
   requirement; document what you used); Delta = cross determinant.
   Count roots s of Delta in F_p vs the number of actual marked triples
   at that (h,k). Report the conditional ratio actual/roots binned by
   (h+k) size and by whether the root count is anomalous (>> 3(h+k)+9).
5. Leading-coefficient check: verify deg Delta_{h,k} = 3(h+k)+9 and
   leading coeff ≡ -c_{k-h-1} mod p on the subsample; count cases where
   p | c_{k-h-1} (apparition events) and report the observed degree drop.

## Deliverables

- Script: research/scripts/q32_marked_dec_scan.py (stdlib; incremental
  output to research/scripts/marked_dec_progress.txt).
- Report: CODEX_MARKED_DEC.md (problems/3.2 root): tables, log-log slope,
  phantom-ratio distribution, apparition-event count, verdict (does the
  data support [GAP-MARKED-DEC] / any structured exception family found).

## Hard constraints

- Do NOT touch proof.tex, ERRATA.md, CAMPAIGN_MAP*, lean/, existing
  scripts/reports. New files only. No external LLM dispatch.
- Exact arithmetic everywhere except final ratios.
- If the row-convention verification (item 4) fails after honest
  adjustment attempts, deliver experiments 1-3 + 5 and a stall note for 4.

## Acceptance

Script exits 0; report matches output; row convention documented with its
20-point verification.
