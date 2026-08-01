# 8-round ChatGPT attack synthesis — 2026-07-31 (dm window, Claude 4.6)

## Setup
Xiang directed a full 8-round ChatGPT assault on P3.2 unconditional proof.
10 questions dispatched (Q6250–Q6260), all answered.

## The single universal obstruction

**S_h(n) = Σ_{n/2<p≤n} e(h·5^{-1}b_n/p) is a sum where the modulus IS the summation variable.**

Every technique in analytic number theory, algebraic geometry, and dynamical systems
requires the modulus to be fixed (or to vary independently of the summation). Here p
serves simultaneously as:
- the modulus of the additive character (e(·/p))
- the summation index
- the field of definition (F_p)

No known framework handles this triple role.

## Route-by-route verdicts

| Q# | Route | Verdict | Specific obstruction |
|----|-------|---------|---------------------|
| Q6250 | Bourgain-Gamburd expansion | NO-GO | Requires averaging over random words; Apéry is one deterministic word per (n,p) |
| Q6251 | Burgess amplification | NO-GO | Burgess uses automorphism of ONE fixed ring; here every p is a different ring |
| Q6252 | Vaughan/Heath-Brown decomposition | NO-GO | Creates composite moduli dm; Apéry recurrence has zero-divisor issues mod composites |
| Q6254 | Complete sum / Korobov-Niederreiter | NO-GO | Apéry mod p is non-autonomous (time-varying coefficients); not an LRS over F_p |
| Q6257 | Mauduit-Rivat automatic sequences | NO-GO | Sum length N ≈ p (the base) → only ONE digit → automatic structure trivial |
| Q6258 | Quadratic splitting (Q(√2)) | NO-GO | Varying modulus kills DL structure; Artin-type order conjecture insufficient |
| Q6259 | Smooth number budget | NO-GO | No Erdős-Kac for holonomic sequences; Lucas is mod-p, not global factorization |
| Q6260 | Archimedean CRT constraint | NO-GO | Height Λn ≈ 3.5n ≫ modulus entropy δn ≈ 0.5n → exponentially many admissible A |

## Why the archimedean route fails quantitatively
The CRT argument (Q6260) gives the sharpest reason for failure:
- Available prime-modulus entropy: at most n/2 (all primes in (n/2, n])
- Integer height: Λn ≈ 3.525n
- Since 3.525 > 0.5, the search interval [0, e^{Λn}] contains exponentially many
  lattice points in every CRT box. Geometry of numbers gives EXISTENCE of solutions,
  not emptiness.
- **A contradiction would need δ > Λ ≈ 3.5**, but δ ≤ 1/2.

## What WOULD be needed
Every ChatGPT answer converges on the same conclusion: a genuinely new theorem is required.
Possible forms:
1. A "mixing theorem for the Apéry recurrence modulo varying primes"
2. An "exponential-sum bound for automatic sequences in short intervals" (N = base regime)
3. A "recurrence-specific anti-concentration theorem" coupling the integer size
   of b_n with its residues modulo many primes
4. A "cross-prime decorrelation theorem" for holonomic sequences

None of these exist in the current literature.

## Conclusion
P3.2 is a genuine open problem. The obstruction is not a missed technique but a
missing theorem. The partial results (block law, valuation law, density-1, localized
exceptional set) are real advances; the full conjecture requires new mathematics.
