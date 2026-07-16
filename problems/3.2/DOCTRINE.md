# DOCTRINE — P3.2 Automode (Session 2026-07-16)

## Goal
Prove the full conjecture: G_n = gcd(d_n a_n, d_n b_n) = e^{o(n)} for ALL n.

Current state: proved for density-1 of n (unconditional). The gap:
- Density-1: log G_n = O(√n) + 3B(n)log n where B(n) = bad prime count
- E[B(n)] = O(1) → o(n/log n) for density-1 by Markov
- ALL n: need B(n) = o(n/log n) pointwise

## Key formula
B(n) = #{p prime, p ∈ (n/2, n] : b_{n-p} ≡ 0 mod p}
     = #{r ∈ [0, n/2) : (n-r) prime AND (n-r) | b_r}

For each r, b_r is a specific integer. The condition (n-r) | b_r asks whether
the specific prime p = n-r is among the prime factors of b_r.

## Avenues

### (a) Arithmetic large-prime-divisor bound
For r < n/2, b_r has at most O(r/log n) prime factors > n/2 (by size).
Key insight: b_r satisfies the Apéry recurrence, so its prime factorization
is NOT arbitrary — it's constrained by the supercongruences and multiplicative
structure.

APPROACH: Bound Ω_{>n/2}(b_r) (count of prime factors > n/2, with multiplicity)
using the Apéry recurrence structure. If Σ_{r < n/2} Ω_{>n/2}(b_r) = o(n/log n),
then B(n) = o(n/log n) for ALL n.

Terminal: either find a proof that works, or exhibit a concrete obstruction
showing this approach fails.

### (b) Prove Hypothesis Z̄ via Chebotarev + gap polynomials
The gap polynomials N_h have Galois groups containing hyperoctahedral groups B_h.
Chebotarev density theorem applied to N_h should give:
- #{p ≤ x : Z(p) ≥ 2h} ≤ c_h · π(x) with c_h → 0 exponentially
- This gives E[Z(p)] = O(1) unconditionally

This doesn't give ALL n, but makes the conditional result UNCONDITIONAL
(for density-1 with quantitative rate).

Terminal: prove or disprove that gap polynomial Galois groups force exponential
decay of P(Z(p) ≥ 2k).

### (c) Second-moment / variance bound
Compute Var(B(n)) = E[B²] - E[B]². If Var = O(1), then by Chebyshev,
B(n) = O(ω(n)) for all but O(N/ω²) values of n, giving strong quantitative
density bounds (e.g., all but N^ε exceptions).

The CRT argument: for distinct primes p,q, the events b_{n-p} ≡ 0 mod p
and b_{n-q} ≡ 0 mod q involve DIFFERENT b-values at DIFFERENT primes.
Independence requires understanding joint distribution of zeros across primes.

Terminal: either prove E[B²] = O(1) or find positive variance lower bound.

### (d) Prove Sym² squareness rigorously
H_p(t) = S_p(t)·A_p(t)² computationally verified for p ≤ 2000.
Should follow from the Sym² structure of the Picard-Fuchs operator.

APPROACH: The Apéry operator L₄ = Sym²(L₂) where L₂ is the Picard-Fuchs
operator of the elliptic pencil E_t: y² = x(x-1)(x-t(1-t)).
The Hasse-Witt matrix of L₄ mod p is the symmetric square of the Hasse-Witt
matrix of L₂. Since Sym²(A) has the same GCD structure as A², this forces
H_p = S_p · A_p².

Terminal: prove or identify the missing step.

### (e) Massive computation push (uisai2)
1. Z(p) to p = 5×10^7 — test Z(p)=14 prediction
2. G_n to n = 10^4 or higher — better empirical growth rate
3. Determinant test for corrected trace identity at Z(p)≥4 primes

### (f) Direct Wronskian/factorization approach
The Wronskian a_n b_{n-1} - a_{n-1} b_n = 6/n³ constrains which primes
can simultaneously divide d_n a_n and d_n b_n. Maybe the multiplicative
structure of 6/n³ forces B(n) = O(log log n).

Terminal: either find a direct bound or show this angle is exhausted.

## Status updates

### Avenue (a) — DEAD (Q5288)
By the Lucas congruence, B(n) = ω_{(n/2,n]}(b_n) + O(1). The problem collapses
to counting prime factors of a single integer b_n in (n/2, n]. No current technique
(recurrence structure, supercongruences, diagonals, sieves) bounds this pointwise.
The full conjecture is equivalent to the TOP-HALF RADICAL ESTIMATE:
log rad_{(n/2,n]}(b_n) = o(n). This is a genuinely new theorem.

### Avenue (b) — DEAD (Q5291)
Chebotarev controls roots of a FIXED gap polynomial N_h as p varies.
But Z(p) requires controlling the actual Apéry orbit for FIXED p.
The quantifiers are reversed. Chebotarev explains the Poisson model
but cannot prove Z̄.

### Avenue (d) — DONE (Caruso et al. 2026)
Sym² squareness is a THEOREM: H_p = Δ^{ε_p} B_p² proved by
Caruso-Fürnsinn-Vargas-Montoya-Zudilin. Updated paper to cite [CFVZ2026].

### Avenue (e) — DONE (computation to 200k)
B(n) ≤ 3 for ALL n ≤ 200,000. Mean 0.065, Var/E = 1.003 (exact Poisson).
Histogram matches Poisson(log 2/log n) to <0.2% TV distance.

## NEW: Quantitative exceptional set (Corollary cor:exceptional)
#{n ≤ N : log G_n > εn} = O(N^{2/3}/ε). This is a POWER-SAVING bound
on the exceptional set, directly from Z(p) = O(p^{2/3}).
The exponent 2/3 is linked: Z(p) = O(p^α) → exceptional set O(max(N^α, N^{1/3})).

## NEW: Dyadic leading-digit bound (Proposition prop:lead)
L_N = O(N^{4/3}/log N) — unconditional, improves old O(N^{3/2}).
Uses D(P,Q) ≤ min(P, Q²)/log P with dyadic decomposition.

## NEW: Improved conditional theorem (Theorem thm:main)
Under Z̄: log G_n = O(√n) for density-1 (was O(ω√n log n)).
The big-prime contribution is now absorbed into the small-prime
floor O(√n). Under Z̄, exceptional set #{log G_n > εn} = O(N^{1/3}/ε).

### Avenue (c) — ANALYZED, no improvement
CRT injectivity gives: for p, q > √N distinct primes, the joint
event {n mod p ∈ Z_p, n mod q ∈ Z_q} occurs for at most Z(p)Z(q)
values of n (since pq > N). But the second-moment bound
Σ B² ≤ R_N + (Σ Z(p))² does NOT improve the first-moment Markov.
The Poisson second moment (Var/E = 1.003) is confirmed computationally
but requires cross-prime independence that we cannot prove.

## Remaining priority
Paper polish. The full conjecture remains open.

## Summary of results
1. Unconditional density-1: G_n = e^{o(n)} (Theorem thm:density1)
2. Unconditional quantitative: #{log G_n > εn} = O(N^{2/3}/ε)
3. Conditional O(√n): log G_n = O(√n) under Z̄ (Theorem thm:main)
4. Conditional exceptional: #{log G_n > εn} = O(N^{1/3}/ε) under Z̄
5. Dyadic incidence: L_N = O(N^{4/3}/log N) unconditional
6. Sym² factorization: theorem (Caruso et al.)
7. Poisson model: Var/E = 1.003 for n ≤ 200k, B(n) ≤ 3
8. Leading-digit bottleneck: N^{1/3} (from D(P,Q) = Q²/log P)

## Resources
- ChatGPT Pro: dm1-dm4 (SOL, 1hr+). Q5293-Q5296 processing.
- uisai2: available for heavy computation
