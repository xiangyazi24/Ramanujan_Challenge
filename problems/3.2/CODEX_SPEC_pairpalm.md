# CODEX SPEC W5: pair-Palm (HM)_k attack — the FULL unconditional conjecture

## Prize

The FULL Apéry GCD conjecture: G_n = e^{o(n)} for ALL n. This is the headline
result. Everything in the paper so far is either unconditional-for-density-1 or
conditional on cross-prime hypotheses. This spec attacks the cross-prime gap.

## The bridge (from paper, already proved)

thm:hm-pointwise: (HM)_k for any k > 6 implies log G_n << n^{2/3+2/k+o(1)}
for EVERY n. In particular:
- k = 8: log G_n << n^{11/12+o(1)} = o(n) ✓
- k = 7: log G_n << n^{20/21+o(1)} = o(n) ✓
So proving (HM)_k for k = 7 or 8 CLOSES the conjecture.

## What is (HM)_k?

(HM)_k: sum_{m < X^2} (K_X(m))_k << X^{2+o(1)} * lambda_X^k

where K_X(m) = #{p in (X, 2X] : p | b_{m mod p}} and lambda_X = sum_{p in (X,2X]} 1/p.
(K)_k denotes the falling factorial.

(HM)_2 is PROVED unconditionally (prop:hm2 in paper, constant 5).
(HM)_3 is the first hard case (W1 characterized it as pair-Palm excess).

## Read first

1. problems/3.2/hm3_result.tex (pair-Palm characterization, anchored-star impossibility)
2. problems/3.2/DOCTRINE.md — search "BREAKTHROUGH 9", "W1 harvest", "pair-Palm"
3. In proof.tex: ssec:high-moment (thm:hm-pointwise, prop:hm2, R_k table),
   thm:polylog, prop:companion-height, prop:reflection
4. problems/3.2/energy_result.tex (column structure for per-prime bounds)
5. problems/3.2/nv_theorem.tex (NV range theorem)

## The pair-Palm characterization (from W1)

thm:hm3-palm-characterization: (HM)_3 is EQUIVALENT to:
  P_3(X) = sum_{p != q} (T_{p,q} - J_{p,q} lambda_{p,q})_+ << X^{2+o(1)} lambda^3

where T_{p,q} = sum_m 1_{p|b_{m mod p}} 1_{q|b_{m mod q}} 1_{exists r: r|b_{m mod r}}
(r ranging over a third prime in (X,2X]).

The "positive excess": for each CRT locus of a (p,q) pair, count the third primes
that ACTUALLY hit vs the Poisson expectation. The excess must be subpolynomial.

## The impossibility (from W1) — what NOT to try

prop:hm3-anchored-star: ANY incidence model satisfying ALL of:
- reflection + no-consecutive + interval O(|I|^{2/3})
- CRT second moment + fiber <= 2 + energy <= 2p
- gap certificates + uniform O(1) row-codegree (<=7)
STILL has sum(K)_3 >> X^3/log^3 X. So single-prime facts cannot prove (HM)_3.
The missing input: cross-prime constraint on WHICH reflection orbit is selected
as the zero fiber.

## Approaches to try (ranked)

1. **Crystalline Frobenius coupling:** The Apéry recurrence mod p and mod q are
   linked by the SAME integer sequence b_n. Two primes cannot independently choose
   "worst" zero fibers. Make this quantitative: for generic (p,q), the joint
   {b_n mod p, b_n mod q} distribution in F_p x F_q approaches uniform (by CRT +
   b_n growth rate). Then the triple count is a second-order effect.

2. **(HM)_4 via double Cauchy-Schwarz:** Instead of (HM)_3, try:
   sum (K)_4 << X^{2+o(1)} lambda^4. Apply CS twice:
   sum (K)_4 <= (sum (K)_2)^{1/2} * (sum (K)_6)^{1/2}. The (K)_2 side is
   proved. The (K)_6 requires 6-wise CRT counting — but all 6 primes distinct
   and > X means their product > X^6 >> X^2, so at most 1 representative per
   6-tuple. This gives sum (K)_6 <= sum_{6-tuples} 1 ~ lambda^6 * X^2.
   CHECK: does this actually work? The CS indices may not align cleanly.

3. **Conditional on Z(p) << sqrt(p):** If W4 succeeds (E << H^{2-delta} ->
   Z << p^{2/3-delta/3}), recompute the (HM)_k threshold. Z << p^{1/2+epsilon}
   lowers it from k > 6 to k > 4. Then (HM)_5 suffices. CS / combinatorial
   approaches more tractable at lower k.

4. **Localized dispersion direct attack:** Instead of (HM)_k, prove the
   equivalent AP-BDH hypothesis (Hypothesis 12 in paper):
   V^o(P,N) << N^{o(1)} * S(P,N). Use the V^o/S -> 1 empirical data (Table
   tab:covariance) as guide. The exact decomposition is:
   V^o = S + E^o where E^o = sum_{p!=q} w_p w_q C_{p,q} and C_{p,q} is the
   centered covariance. Prove |E^o| = o(S), i.e., the signed covariance sum
   cancels. The data shows |E^o|/S = O(1/sqrt(N)).

5. **Large sieve with Apéry structure:** The bilinear Kloosterman form
   sum S_p(a) S_q(b) K_N(a/p + b/q) -- the CRT error term. If the Fourier
   transforms S_p(a) of the zero sets satisfy a nontrivial bound (they should,
   given shift-correlation randomness), then the large sieve inequality gives
   the right bound. The classical large sieve gives E^o << N^2/log N, losing
   one log. Power-saving in S_p(a) (e.g., from NV/fiber theorems) may close it.

## Deliverables

Produce `pairpalm_result.tex` with:
- Whatever you prove, stated precisely with complete proofs
- If you close (HM)_k for some k > 6: this is the headline theorem
- If partial: the sharpest unconditional bound on sum(K)_k, the exact lemma
  that remains, and computational verification
- A verification script `pairpalm_verify.py`

## Verification

```
cd ~/repos/Ramanujan_Challenge/problems/3.2 && python3 pairpalm_verify.py
```
All checks must PASS.

## Stall protocol

If stuck: deliver partial structural results + impossibility analysis + the
minimal sufficient lemma for (HM)_7, clearly stated. Do not fabricate proofs.
Honest partial results with sharp next targets are valuable.
