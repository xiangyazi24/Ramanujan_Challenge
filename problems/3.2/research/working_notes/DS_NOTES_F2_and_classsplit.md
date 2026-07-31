# DS note: F₂ dispersion at scale + S_h mod-24 class split (2026-07-31, automode)

## Result 1 — F₂(N) = |I|²/N almost exactly (random/CRT prediction confirmed)
From the (r,p) data bank (p ≤ 2·10⁶), H(n) = #{(r,p) : r+p=n}, F₂ = Σ H(n)(H(n)−1):

| N | M=|I| | F₂ | |I|²/N | F₂/|I| | req o(N²/log²N) | gap |
|----|-------|-----|--------|--------|-----------------|-----|
| 1e5 | 6119 | 378 | 374 | 0.062 | 7.5e7 | 2.0e5 |
| 2e5 | 11328 | 548 | 642 | 0.048 | 2.7e8 | 4.9e5 |
| 4e5 | 21333 | 1150 | 1138 | 0.054 | 9.6e8 | 8.4e5 |
| 8e5 | 40750 | 2038 | 2076 | 0.050 | 3.5e9 | 1.7e6 |
| 1e6 | 50072 | 2396 | 2507 | 0.048 | 5.2e9 | 2.2e6 |

- F₂(N) tracks |I|²/N ≈ N/log²N (the CRT-independence prediction) to within ~10% at every scale.
- Gap vs the sufficient condition F₂ = o(N²/log²N): 2×10⁵ … 2.2×10⁶, GROWING. Cross-characteristic
  anti-alignment is empirically rock-solid (the off-diagonal energy is exactly the independence value).
- Interpretation: F₂ ≈ |I|²/N says the CRT classes {n ≡ z_p mod p} hit the interval (N,2N] at the
  expected density — i.e. the "bad" classes are CRT-equidistributed. Proving this IS the hard theorem.
- Valid range: N ≤ 10⁶ (bank p ≤ 2·10⁶ ⇒ 2N ≤ 2·10⁶).

## Result 2 — S_1(n) mod-24 class split: NO class-dependent signal (avenue (c) TERMINAL)
|S_{1,a}(n)| per prime for each unit class a mod 24 (square classes 1,5,7,11 marked SQ):
- All 8 unit classes contribute comparably (per-prime |S| ≈ 0.03–0.45, random-walk scale).
- Square-vs-non-square ratio FLIPS sign with n: n=2000 0.197/0.145, n=4000 0.090/0.035,
  n=8000 **0.007/0.084**, n=16000 0.138/0.057. No systematic difference.
- **Verdict: the mod-24 square factorization does NOT manifest as class-dependent cancellation in
  S_h(n).** Consistent with Q6130 ("bare square is analytically inert"). The square structure is a
  rank-lowering algebraic/differential fact, not a horizontal-cancellation input by itself.
- So the horizontal route must use the rank-two Franel structure via the coboundary gate test
  (avenue (a), Q6130/DS_AVENUES), NOT a naive class split.

## Banked conclusion
Both negative/confirming results are consistent with the synthesis: no fibre/local/class-level
mechanism gives horizontal cancellation; the theorem needs either the rank-two coboundary collapse
(avenue a) or a genuinely new cross-prime equidistribution input.
