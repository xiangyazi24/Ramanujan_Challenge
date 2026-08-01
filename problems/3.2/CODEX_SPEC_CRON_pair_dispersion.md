# CODEX SPEC: empirical dispersion ledger for the centered pair correlation (k=2)

## Goal
Pre-build the numerical infrastructure for the dispersion-route verdict (ledger `CRON_FRESH_EYES_pointwise.md` appendix R.3 [LIVE-1], appendix T): measure the centered pair sum

  D(X) = Sum_{p != q in P_X} Sum_r g_p(r) g_q(r + p - q) restricted to n = p + r in I_X,

where I_X = (X, 2X], P_X = primes in (X/2, 2X], Z_p = zero set of Apery b_r mod p (r < p),
g_p(r) = 1_{Z_p}(r) - |A_{p,X}|/|I_X| is the centered row indicator (center so that Sum_n over I_X of the centered row vanishes; equivalently D(X) = Sum_{p!=q} [|A_p ∩ A_q| - |A_p||A_q|/X]).

This is exactly M_2-minus-independent-benchmark, but the deliverable is its DECOMPOSITION BY SHIFT SCALE, which no previous scan produced.

## Outputs (CRON_pair_dispersion_report.md)
For each dyadic X = 2^13 .. 2^19:
1. D(X) total, and the Poisson-normalized D(X)/(X*lambda_X^2).
2. Breakdown by dyadic shift scale: for j = 0,1,2,..., D_j(X) = the same sum restricted to 2^j <= |p-q| < 2^{j+1}. Report the table D_j(X) with counts of contributing (p,q,r) triples per shell.
3. Sign pattern and cancellation diagnostic: for each shell also report Sum |A_p ∩ A_q| (unsigned) so we can see how much cancellation the centering achieves per shell.
4. A verdict line per X: is D(X) consistent with sqrt-cancellation (|D| ~ sqrt(unsigned mass)) per shell, or is there a structured shell (e.g. small |p-q|) carrying a bias?
5. Special diagnostic — the reflection channel: pairs (p,q) where the intersection point n has n-p and n-q BOTH at reflection-symmetric positions (r' = p-1-r pairing). Count these separately; reflection is the one known algebraic structure and could create a biased diagonal.

## Method
Z_p via the (n!)^3 b_n Montgomery kernel you already built (CRON_garqi_moments.c) — reuse/adapt; keep zero pairs (p, r) in memory (~82k pairs for p <= 2^20, trivial). Then pair correlations are computed from the hit lists per n: for each n in I_X collect the list of hitting primes (you already build H(n); now keep the actual prime lists for n with H(n) >= 2), and D(X) = Sum_n H(n)(H(n)-1) - (independent benchmark term) with the shift-shell split from the actual prime pairs at each n. The unsigned per-shell masses need the benchmark integral per shell: compute the expected count E_j(X) = Sum_{|p-q| in shell} |A_p||A_q|/X directly by summing over prime pairs (78k primes -> ~3e9 pairs is too many: bin |A_p| by prime into dyadic |p-q| shells using prefix sums over sorted primes — O(P log P), do NOT loop over all pairs).
Threads optional (the heavy part is only Z_p). Progress to stderr every ~5 s.

## Sanity gates (MUST pass, abort loudly otherwise)
a. Reproduce M_2(X) for X = 4000 (=18), 8000 (=54), and the appendix-T values at X = 2^13 (M_2 = 58) and 2^19 (M_2 = 1458) from your hit lists.
b. Check per-X: Sum_j [signed shell sums] + benchmark = M_2 identity closes exactly (integer arithmetic where possible).
c. Spot-check 5 pairs (p,q) with nonempty intersection against the definition |A_p ∩ A_q| = #{r in Z_p : r + p - q in Z_q, p + r in I_X}.

## Discipline
- New files only, CRON_ prefix: CRON_pair_dispersion.c (or .py driver + C kernel), CRON_pair_dispersion_report.md.
- Work only in problems/3.2/. Commit when done, report SHA + paste the final tables.
- Do not touch CRON_FRESH_EYES_pointwise.md (single-writer: cron session owns it).
- Do NOT call the ChatGPT bridge / ask-gpt.py or dispatch sub-questions to any oracle. Pure computation task. If blocked, write the blocker into the report and stop.
