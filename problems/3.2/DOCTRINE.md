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
The exponent 2/3 is linked: Z(p) = O(p^α) → exceptional set O(max(N^α, N^{1/(3-α)})).
Crossover at α = √2 - 1 ≈ 0.414. At α=0: N^{1/3}. At α=2/3: N^{2/3}.

## Dyadic leading-digit bound (Proposition prop:lead) — CORRECTED
L_N = O(N^{10/7}/log N) — unconditional, improves old O(N^{3/2}).
Uses D(P,Q) ≤ min(PQ^{2/3}, Q²)/log P with dyadic decomposition.
**Bug fixed**: prior "prime counting" bound D(P,Q) ≤ π(2P) was unjustified
(each prime can contribute multiple (p,q) pairs). Replaced with sub-interval
zero count bound from prop:zp-bound applied to (Q,2Q]: each prime contributes
O(Q^{2/3}) zeros in the interval. Crossover P₀ = N^{4/7}.

## Improved conditional theorem (Theorem thm:main)
Under Z̄: log G_n = O(√n) for density-1 (was O(ω√n log n)).
The big-prime contribution is now absorbed into the small-prime
floor O(√n). Under Z̄, exceptional set #{log G_n > εn} = O(N^{1/3}/ε).
**Key insight (2026-07-16):** Under Z̄, D(P,Q) ≤ Σ Z(p) = O(P/log P),
which is tighter than the unconditional PQ^{2/3}/log P. Crossover shifts
from P₀=N^{4/7} to P₀=N^{2/3}, giving L_N = O(N^{4/3}/log N) under Z̄.
This improves the conditional exceptional set from N^{3/7} to N^{1/3}.

### Avenue (c) — ANALYZED, no improvement
CRT injectivity gives: for p, q > √N distinct primes, the joint
event {n mod p ∈ Z_p, n mod q ∈ Z_q} occurs for at most Z(p)Z(q)
values of n (since pq > N). But the second-moment bound
Σ B² ≤ R_N + (Σ Z(p))² does NOT improve the first-moment Markov.
The Poisson second moment (Var/E = 1.003) is confirmed computationally
but requires cross-prime independence that we cannot prove.

## ★ BREAKTHROUGH: Codegree amplification (Q5342, verified Q5351)
**Theorem (polylog exceptional set, thm:polylog):**
#{n ≤ N : log G_n > εn} = O_ε((log N)²) — UNCONDITIONAL, power exponent ZERO.

Key ingredients:
1. **Codegree lemma (lem:codegree):** For m,n ∈ (N,2N] with h=n-m, the number
   of common top-half bad primes > P₀ is O(h log N / log P₀). Reason: all such
   primes divide the fixed integer N_h(m), which has height (CN³)^h.
2. **Leading digit pointwise (lem:lead-pointwise):** Σ_{p>√n, p|b_{⌊n/p⌋}} log p = o(n).
3. **Localized Cauchy-Schwarz:** Partition (N,2N] into intervals of length Y = ε²N/logN.
   Codegree upper bound on pair incidence: Σ C(d_p, 2) ≤ CYM².
   CS lower bound: Σ C(d_p, 2) ≥ ε²NM²/(512 log N) - I/2.
   Forces M = O(1/ε²) per interval.
   Total: O(log N/ε⁴) per dyadic block, O((log N)²) overall.

This supersedes cor:exceptional (O(N^{2/3})) and gives upper Banach density zero
unconditionally (previously needed Z̄).

**Open:** O(1) codegree (would give O(log N) not (log N)²) requires new arithmetic input.
Q5351 §7 explains why: N_h(m) genuinely has height m^{3h}, so O(h) prime factors is sharp.
Triple codegree (gcd of two gap polynomials) dispatched to dm1 — may break barrier.

## ★ BREAKTHROUGH 2: Block system and leading-digit vanishing (Q5365-Q5367)

### Block system (Lemma lem:block)
For p ≥ 7, p ≤ n < p², write n = qp + r. Define D_q ≡ p³a_{qp} (mod p).
Three-part lemma:
1. **Propagation:** p³a_{qp+r} ≡ D_q · b_r (mod p)
2. **Wronskian boundary:** D_q b_{q-1} - D_{q-1} b_q ≡ 6/q³ (mod p)
3. **GCD support:** p | G_n ⟺ D_q · b_r ≡ 0 (mod p)

### Leading-digit vanishing (Lemma lem:leading-vanish)
If p | b_q and p ∤ b_r, then v_p(G_n) = 0.
Proof: When p | b_q, the Wronskian forces D_q ≡ 6/(q³ b_{q-1}) (mod p),
which is a unit. So D_q · b_r ≢ 0 (mod p).
**Impact:** Leading-digit bad primes contribute NOTHING to the GCD.
Verified computationally for n ≤ 400 (0 violations).

### Phantom primes (companion-block mechanism)
Primes where D_q ≡ 0 (mod p) with p ∤ b_q and p ∤ b_r.
These create GCD contributions invisible to ordinary Lucas congruences.
Example: p=13, n=27, D_2 ≡ 0 (mod 13).
- Computationally: phantom(n) ≤ 3 for n ≤ 600, mean ≈ 0.46
- Heuristically: Poisson(1/2) — expected phantom per n is Σ_q 1/(q log(n/q)) ≈ 1/2
- First-moment: P_N = Σ phantom(n) = O(N^{3/2}/log N) — pessimistic, actual growth much slower
- **OPEN:** No rigorous proof that phantom(n) = O(1). This is the remaining gap.

### Paper fixes (this session)
1. **v_p(G_n) = 3 → v_p ≤ 3:** Generic simple zero gives v_p = 1, not 3. Safe upper bound v_p ≤ 3.
2. **One-per-quotient removed:** Multiple primes per quotient can divide b_q. Moot since leading-digit bad primes have v_p(G_n) = 0.
3. **Upper Banach density fixed:** Added uniform local bound sup_M #{M < n ≤ M+L : log G_n > εn} ≪ L^{2/3} + 1 from critical-window estimate.
4. **Paper structure:** Added lem:block, lem:leading-vanish. 21 pages, compiles clean.

### Top-half formula (clean region)
For q=1 (p > n/2): D_1 = 6 (always unit), so p | G_n ⟺ p | b_{n-p}.
No phantom issues. This is the domain of the codegree amplification.

### Rainbow pile-up barrier (why ALL n is hard)
A single n could have many fresh bad primes not shared with neighbors.
The codegree method can't exclude this: each codegree pair comes from a
different gap polynomial N_h with independent roots. Poisson model strongly
predicts no exceptions, but this is a methodological barrier.

## Remaining priority
Paper is correct and strong (polylog exceptional set, upper Banach density zero).
Main open: close gap from O((log N)²) to ALL n.
Secondary: rigorous bound on phantom(n).
Paper is 21 pages, compiles clean (zero warnings).

### Paper edit summary (this session, 2026-07-16)
- Abstract: added finite harmonic weight, O(N^{1/3}) + upper Banach density zero, O(N^{4/3}/log N) under Z̄
- Proposition prop:Zp-data: factorial moments table (r=1..6 vs 2^{-r})
- Hypothesis hyp:poisson: formal Poisson conjecture (factorial-moment formulation)
- Remark rem:poisson-implies: Poisson implies Z̄ implies all conditional results
- Corollary cor:harmonic: finite harmonic weight (Σ 1/n < ∞ for E_ε)
- Remark rem:exponent-link: corrected crossover to (3-√5)/2, added full general formula
- Corollary cor:cond-except: improved from O(N^{3/7}) to O(N^{1/3}), new proof via D(P,Q) ≤ Σ Z(p)
- Corollary cor:banach: upper Banach density zero under Z̄ (NEW)
- Remark rem:pointwise-gap: Borel-Cantelli observation (B(n)≥1 i.o. compatible with conjecture)
- Remark rem:squareness: Sym² low-half vacuity principle, random model vs distribution theorem
- Remark rem:literature: comparison with BCZ2003, irrationality literature, novelty claim
- Bibliography: BCZ2003, Zudilin2004 added
- Theorem thm:main proof: N^{3/7}→N^{1/3}, N^{1/14}→N^{1/6}

## Summary of results
1. Unconditional density-1: G_n = e^{o(n)} (Theorem thm:density1)
2. ★ **Polylog exceptional set**: #{log G_n > εn} = O_ε((log N)²) (Theorem thm:polylog) — power exponent ZERO
3. **Upper Banach density zero** UNCONDITIONAL (from (2), since (log N)²/N → 0)
3b. **Finite harmonic weight**: Σ_{n∈E_ε} 1/n < ∞ (Corollary cor:harmonic, also from (2))
4. Conditional O(√n): log G_n = O(√n) under Z̄ (Theorem thm:main)
5. Conditional exceptional: #{log G_n > εn} = O(N^{1/3}/ε) under Z̄
6. Dyadic incidence: L_N = O(N^{10/7}/log N) unconditional, O(N^{4/3}/log N) under Z̄
7. Sub-interval zero count: #{j ∈ I : p|b_j} = O(|I|^{2/3}) for any I ⊂ {0,...,p-1}
8. Tower density: #{n < p^r : p|b_n} = p^r - (p-Z(p))^r (exact, from Lucas)
9. Sym² factorization: theorem (Caruso et al.)
10. **Sym² vacuity** (NEW): factorization cannot improve D(P,Q) (low-half is formally free)
11. Poisson model: Var/E = 1.003 for n ≤ 200k, B(n) ≤ 3
12. **Poisson conjecture** (NEW): Hypothesis hyp:poisson with factorial-moment formulation
13. Z(p) data: max Z(p)=12 at p=159977 for p ≤ 10^6; Poisson(1/2) fit to 6 decimal places
14. General exceptional formula: O(max(N^α, N^{1/(3-α)})/ε) for Z(p) ≤ Cp^α
15. **Literature comparison** (NEW): no prior power-saving exceptional-set in variable-coefficient holonomic setting

## Avenues surveyed (ALL dead for Z(p) improvement beyond O(p^{2/3})):
- Q5324: Convolution/Sym² → dead (binomial counterexample)
- Q5328: Mesoscopic roots → dead (disjoint-root obstruction)
- Q5329: p-adic tower → dead (vertical/horizontal mismatch)
- Q5330: Trace functions/FFK → dead (p-adic vs ℓ-adic)
- Q5331: Second moment → dead (CRT equidistribution invalid)
- Q5332: Tower level-2 → dead (gap polynomial is one-period, Lucas cylinders block)

## ChatGPT Q&A processed (this session)
- Q5337 (dm2): Sym² divisor bound → dead (low-half vacuity principle). Added to paper.
- Q5338 (dm3): Literature comparison → no prior Apéry exceptional-set analogue. Added Remark rem:literature + BCZ2003/Zudilin2004 refs.
- Q5339 (dm4): Intermediate statements → finite harmonic weight + upper Banach density hierarchy. Added Corollary cor:harmonic + Hypothesis hyp:poisson.

## ChatGPT Q&A processed (continued)
- Q5342 (dm1): ★ Bilinear improvement + dual Cauchy-Schwarz → BREAKTHROUGH. O_ε((log N)²) exceptional set.
- Q5343 (dm2): Density hierarchy verification → H and B INCOMPARABLE, confirmed paper correct.
- Q5344 (dm3): Apéry large sieve → ALS Δ_N=O(1) WRONG, pair correlation correct target.
- Q5345 (dm4): Paper additions → caught ordinary-vs-all-prime Poisson bug, fixed to K*(p).
- Q5348 (dm2): Local zero count → S(M,L) ≪ ML^{2/3}/log M, d*(E_ε)=0 unconditional.
- Q5350 (dm4): Pair correlation LaTeX → ready-to-paste remark (deferred, not essential now).
- Q5351 (dm1): ★ Verification of codegree amplification → CONFIRMED, height corrected to (CN³)^h.

## ChatGPT Q&A processed (block system session)
- Q5365 (dm1): Block system proof of leading-digit vanishing. D_q constants, counterexample p=13 n=27. "Dominant summand" proof WRONG but conclusion RIGHT via block Wronskian.
- Q5366 (dm2): Confirms block system. Top-half-only polylog proof (C_1=6 always unit). Paper-ready theorem.
- Q5367 (dm3): Same block analysis with λ_{p,q} notation. Confirms formula is false, provides manuscript changes.

## ★ BREAKTHROUGH 3: Reflection theorem and companion-height bound (Q5369-Q5371)

### Reflection theorem (Proposition prop:reflection) — PROVED
For every prime p ≥ 5 and every p-integral solution u of the Apéry recurrence:
  u_{p-1-j} ≡ u_j (mod p) for all j.
In particular: a_{p-1} ≡ 0 (mod p) and a_{p-2} ≡ 6 (mod p).

**Proof (elementary, 10 lines):**
1. P(-1-x) = -P(x) ⟹ reversal R maps solutions to solutions mod p.
2. b_{p-1}≡1, b_{p-2}≡5 ⟹ Rb = b (same initial values, uniqueness).
3. Casoratian W_n(Ru,b) = -W_{p-n}(u,b) = W_n(u,b) ⟹ Ru - u = cb.
4. R(Ru-u) = -(Ru-u) but R(cb) = cb. Since p is odd, c = 0. QED.

**Verified computationally:** All p ≤ 500, all 7 check categories pass.

### D_q = a_q identity (Remark rem:block-apery)
Block constant D_q = a_q (mod p). PROVED by Wronskian uniqueness when
p ∤ b_{q-1}. Verified computationally for all (p,q) tested. ChatGPT
caveat (Q5369 §6): at Hasse-zero indices (p | b_{q-1}), Wronskian
alone doesn't determine D_q — needs separate block/Frobenius argument.

### Companion-height bound (Proposition prop:companion-height) — PROVED
Total companion-block weight at any index: O(n^{2/3}).
Split by quotient Q = n^{1/3}: small quotients ≤ Q contribute O(Q²) = O(n^{2/3});
large quotients have all primes < n^{2/3}, contributing θ(n^{2/3}) by Chebyshev.
**Impact:** Companion-block primes are o(n) at EVERY index. Paper updated.

### W(p) = #{q : a_q ≡ 0 mod p} statistics
Always ≥ 1 (since a_{p-1} ≡ 0). Max W(p) = 7 for p ≤ 2000. Always odd (reflection).
Histogram: W=1 (66%), W=3 (27%), W=5 (5%), W=7 (1%).
OPEN: prove W(p) = O(1) rigorously.

### Rainbow pile-up barrier (Q5370) — CONFIRMED
Codegree arguments cannot close the gap to ALL n. The exact remaining problem:
**Pointwise radical estimate:** Σ_{p > √n, p | b_{n mod p}} log p = o(n) for every n.
After removing companion-block channel (O(n^{2/3})), this is the ONLY remaining target.
No pair-codegree, no moment bound, no Chebotarev argument can prove this —
need a genuinely new one-point rainbow estimate.

### Paper changes (Breakthrough 3)
- Added Proposition prop:reflection (reflection theorem) — 15 lines including proof
- Added Proposition prop:companion-height (companion-block height bound) — 12 lines including proof
- Updated thm:polylog proof to reference prop:companion-height instead of computation
- Updated thm:density1 proof to reference prop:companion-height
- Updated remark 6 (palindromic symmetry) to reference prop:reflection
- Added remark rem:rainbow-open (remaining open problem)
- Paper: 23 pages, compiles clean, zero warnings

## ChatGPT Q&A processed (reflection/gap session)
- Q5369 (dm2): ★ Reflection theorem proof (reversal + Casoratian). Paper-ready LaTeX.
- Q5370 (dm3): Rainbow barrier analysis. Companion-height O(n^{2/3}). Five approaches to close gap.
- Q5371 (dm4): Independent confirmation of reflection theorem (symplectic formulation).

## Currently dispatched (all 4 dm tabs processing)
- dm1: still processing (from earlier session — triple codegree or prior question)
- dm2: W(p) boundedness — can we prove #{q : a_q ≡ 0 mod p} = O(1)?
- dm3: Pointwise radical theorem approaches — large sieve, single-index certificate, propagation
- dm4: Companion-height clean proof — paper-ready with optimal constants

## Summary of results (updated)
1-15. [previous results unchanged]
16. ★ **Reflection theorem** (NEW): a_{p-1-j} ≡ a_j (mod p), hence a_{p-1} ≡ 0 (mod p)
17. ★ **Companion-height bound** (NEW): companion-block weight O(n^{2/3}) at every index
18. **Remaining open problem** clearly identified: pointwise radical estimate

## Resources
- ChatGPT Pro: dm1-dm4 (SOL). All tabs dispatched.
- uisai2: Z(p) computation COMPLETED. Max Z(p)=12 at p=159977 through p~1.3M.
