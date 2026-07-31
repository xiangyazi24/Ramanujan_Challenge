# DS note: S_h is exactly a reciprocal-prime sum (verified) — circularity of the horizontal route

Date 2026-07-31 automode. Source: Claude4.6's dm Q6170 (his ChatGPT) + DS independent verification.

## The exact reduction (Gessel–Lucas / block law)
For p ∈ (n/2, n], n = p + r, the block law b_{p+r} ≡ 5·b_r (mod p) gives
`b_r ≡ 5^{-1}·b_n (mod p)` (5 invertible mod p for p ≥ 7).
Hence the horizontal Weyl sum is EXACTLY a reciprocal-prime sum with the single enormous frequency A = h·5^{-1}·b_n:
```
S_h(n) = Σ_{n/2<p≤n} e(h·b_{n−p}/p) = Σ_{n/2<p≤n} e(h·5^{-1}·b_n/p).
```
**DS verified**: b_r ≡ 5^{-1}b_n (37/37 pairs); S_1(n) = Σ e(5^{-1}b_n/p) exactly (n=100,200,400, |S| matched to <1e-9).

## Why this is the frontier, not a route
- **No estimate o(P_n) can hold uniformly in the frequency A**: if A is divisible by every prime in the window, every phase is 1. So a minor-arc/uniform estimate is impossible; the saving MUST come from the specific arithmetic of b_n — i.e. that b_n is NONRESONANT with (n/2,n], which is exactly the large-prime-divisor conjecture (H(n) = o(n/log n)).
- Vaughan/Type-I decomposition fails: the division-free recurrence needs (r!)^{-3} mod q for composite q, which doesn't exist; the prime-only Lucas congruence has no composite analogue.
- So S_h(n) = o(P_n) ⟺ H(n) = o(P_n): the horizontal Weyl route is EXACTLY the conjecture in reciprocal-prime disguise. This confirms the campaign's terminus: no local/structural mechanism can unlock it; the statement itself (b_n nonresonant) is the irreducible content.

## What this buys
- The cleanest statement of the obstruction: "the specific exponentially-large integer b_n is nonresonant with the primes in (n/2,n]" — a large-prime-divisor / nonresonance theorem for a holonomic integer. No such theorem exists (Q6129: only the trivial (1/2+o(1))P_n bound; no holonomic-prime-divisor theorem; factorial-ratio counterexamples).
- Empirically TRUE (|S_h|~√P_n, D_1≈P_n) — b_n IS nonresonant.
- Coordinates with Claude4.6 (he found it via dm Q6170); the four classes mod 5 (from 5^{-1}) split the sum.

## Status
Horizontal route (S_h = o(P_n)): TERMINAL at the reformulation — it IS the conjecture. The reciprocal-prime form is the best statement of what must be proved, and the new-math frontier is a "nonresonance of a holonomic integer" theorem.
