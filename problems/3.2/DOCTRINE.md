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

## ★ BREAKTHROUGH 4: Certificate gap and the b_n radical obstruction (Q5373, Q5374, Q5376)

### W(p) = O(1) is HEURISTICALLY UNLIKELY and UNNECESSARY (Q5373)
- **Poisson(1/2) model:** The histogram W=1,3,5,7 fits Poisson(1/2) to <2% TV.
  This predicts UNBOUNDED support — W(p) → ∞ along rare primes.
- **Unconditional theorem:** W(p) ≪ p^{2/3} from gap polynomials (same framework as Z(p)).
- **Even W(p) = O(1) would NOT suffice:** quantifier inversion — one n can have ~n/log n
  companion primes even if each prime has ≤3 zeros (Q5373 §9).
- **BUT: companion channel is ALREADY O(n^{2/3})** by prop:companion-height. So W(p) is moot.

### The sole remaining target: W(n) = o(n) (Q5374, Q5376)
Define W(n) = Σ_{p > √n, p | b_{n mod p}} log p (lower-digit bad-prime weight).
By Lucas: every bad prime p | b_{n mod p} also divides b_n itself. So W(n) ≤ log rad_{>√n}(b_n).

**CRITICAL COMPUTATION (this session):**
- log rad_{>√n}(b_n) ≈ 3.5n (NOT o(n)!) — almost all factorization weight is in large primes
- BUT W(n) is empirically O(log n) — typically 0, max ≈ 8.67 for n ≤ 150
- MOST large primes dividing b_n are quotient-digit primes (p | b_q), which are HARMLESS for G_n
- The certificate bound W(n) ≤ log rad_{>√n}(b_n) = O(n) is too loose by a factor of n/log n

### Why log rad_{>√n}(b_n) = O(n) is NOT the right target (Q5376)
The GCD sum Σ log gcd(n-r, b_r) is also NOT o(n) — even the prime 5 contributes Ω(n).
The right decomposition (Q5376 eq. 0.7):
  W(n) = O(n^{2/3}) + Σ_{p > √n, p | b_n, p ∤ b_{⌊n/p⌋}} log p
The second term is the "pure lower-digit" weight — primes dividing b_n where
the quotient-digit sequence does NOT explain the divisibility.

### Key structural facts
1. **Injectivity:** p ↦ (n mod p) is injective on primes > √n (Q5376 §1)
2. **At most one large prime factor per residue:** n-r can have at most one prime > √n (Q5376 §1)
3. **No polynomial certificate exists:** b_n is the natural certificate, but log b_n = O(n) (Q5376 §4)
4. **Palindromic gap polynomial is tautological:** reflected gap gives forced linear factor (Q5376 §5)
5. **No recurrence-only one-hit invariant exists:** gap polynomial is intrinsically two-hit (Q5374 §2)

### Hard-core pruning (Q5374 §3)
Endpoint zeros (s_p < R), near-central zeros (h_p < H), and large-quotient (q > Q)
are all O(n^{2/3}) with Q=R=H=n^{1/3}. The remaining "hard core" has:
- p ∈ (n^{2/3+ε}, n/2)
- zeros in the bulk (away from endpoints and center)
No rigorous bound for this hard core beyond O(n).

### Two viable routes (Q5374 §9)
1. **New subexponential certificate:** An integer R_n divisible by every bad prime,
   with log R_n = o(n). No candidate beyond b_n is known.
2. **Rainbow large sieve / growing-moment theorem:** A spectral or high-moment
   nonconcentration bound for one row of the incidence matrix. Needs new arithmetic input.

## ChatGPT Q&A processed (post-reflection session)
- Q5373 (dm3): ★ W(p) analysis — Poisson(1/2) heuristic, W(p) ≪ p^{2/3} theorem,
  quantifier inversion, companion channel harmless. Paper-ready proposition.
- Q5374 (dm4): One-hit invariant impossibility, hard-core pruning theorem, K3 Weil
  correction (O(p) not O(p^{3/2})), spectral large sieve formulation.
- Q5375 (dm1): Clean companion-height proof with optimal constants. Confirmed.
- Q5376 (dm1): ★ GCD sum FALSE as o(n) target (p=5 contributes Ω(n)). Certificate
  is b_n with exponential height. Quotient-relative decomposition eq (0.7).

## Currently dispatched (4 dm tabs processing)
- dm1: Lower-digit radical proof approaches (BV, CRT, Sato-Tate)
- dm2: Bombieri-Vinogradov angle for Apéry divisibility array
- dm3: Pointwise radical theorem (comprehensive, 5 specific approaches)
- dm4: Large-prime radical of b_n (Stewart, Corvaja-Zannier, hypergeometric motives)

## ★ BREAKTHROUGH 5: Second-moment reduction + Poisson(1/2) verification (Q5377-Q5383)

### |Z_p| distribution: 2·Poisson(1/2) — VERIFIED (all p ≤ 10000)
Computed Z_p for ALL 1227 primes up to 10000. Distribution of K(p) = Z(p)/2:
- k=0: 61.7% (Poisson: 60.7%)
- k=1: 30.2% (30.3%)
- k=2: 6.8% (7.6%)
- k=3: 1.1% (1.3%)
- k=4: 0.2% (0.16%)
Var/E = 0.989 ≈ 1 (exact Poisson). Max |Z_p| = 8 (at p=3727, 6841).

### ALL classical tools EXHAUSTED (Q5377, Q5378, Q5379)
Confirmed DEAD for closing the gap to ALL n:
- Bombieri-Vinogradov: wrong quantifiers (controls Σ_p, not single n)
- Large sieve: Q² barrier at top scale p~N → N²/log N, need N²/log²N
- Turán-Kubilius: requires pq < N, but all relevant primes have pq > N
- Sato-Tate: controls only central coefficient (p | a_p(f)), density 0
- Pretentious theory: B(m) is not multiplicative
- Gap polynomials: forced central factor makes them tautological
- Stewart/Corvaja-Zannier: requires constant-coefficient recurrence
- abc/smooth numbers: rad(b_n) = O(n) is trivially true but useless
- Hypergeometric motives: Katz equidistribution needs fixed rank

### Second-moment criterion — PROVED (Proposition 10.1 in paper)
**Theorem:** M₂(N) = O(N^{2-δ}) for any δ > 0 → W(n) = o(n) for ALL n.
Proof: B(n)² ≤ M₂(N), so B(n) ≤ N^{1-δ/2} = o(N).

### M₂(N) ≈ 0.59·N·(log N)² — COMPUTED (all dyadic blocks N=64..4096)
Ratio M₂/(N·log²N) stabilizes at 0.59 ± 0.02 across all dyadic blocks.
Decomposition: diagonal ≈ 63% (from weight factor (log p)²), off-diagonal ≈ 37%.

### No-go model — PROVED (paper §10.4)
Adversarial Z*_p = {m₀ mod p, p-1-m₀ mod p} satisfies ALL vertical bounds
(|Z*|=2, reflection, nearest-neighbor exclusion) yet M₂* = Ω(N²).
**Implication:** Vertical moments ALONE cannot close the gap. Must use
cross-prime HORIZONTAL information about zero locations.

### AP-BDH dispersion — FORMULATED (Hypothesis 12 in paper)
V°(P,N) ≪ N^{o(1)}·S(P,N). With current Z ≤ p^{2/3}: gives W(n) ≤ N^{5/6+o(1)} = o(N).
This is a cross-prime equidistribution theorem for Apéry zero-sets.

### Bilinear Kloosterman form — IDENTIFIED (paper §10.5)
CRT error = Σ_{(a,b)≠(0,0)} S_p(a)·S_q(b)·K_N(a/p+b/q)
where S_p(a) = (1/p)Σ_{r∈Z_p} e(ar/p) is the Fourier transform of Z_p.
Classical large sieve: O(N²/log N) — fails by one log.
Bettin-Chandee / DFI-type power-saving: closest framework, not yet applicable.

### Gap polynomial sparsity — VERIFIED
For all primes p ≤ 10000 with Z(p) ≥ 2 and gaps h ≤ 100:
Σ_{h≤100} r_p(h) ≤ 1. Zeros are too sparse to generate gap roots.
Hypothesis: Σ r_p(h) ≪ H → Z(p) ≤ √p (paper Hypothesis 13).

### Paper update: Section 10 added (242 lines, 26 pages total)
- Proposition 10.1: second-moment criterion
- M₂ data table (7 dyadic blocks)
- CRT expansion (diagonal/off-diagonal decomposition)
- No-go model (adversarial zero-sets)
- Fourier form of CRT error (bilinear Kloosterman)
- AP-BDH dispersion hypothesis + sufficiency proof
- Gap polynomial sparsity hypothesis
Committed at a89094c.

### ChatGPT Q&A processed (this session)
- Q5377 (dm3): ★ Comprehensive approaches analysis. No-go model. BV/large sieve barriers.
- Q5378 (dm2): ★ Second-moment criterion + CRT error + bilinear form. All confirmed.
- Q5379 (dm4): Stewart/CZ/hypergeometric — all DEAD for lower-digit radical.
- Q5380 (dm1): ★ CRT energy. Fourier expansion. Deligne obstruction. Exact breakthrough targets.
- Q5381 (dm2): Poisson factorial moments. De-averaging. k=2 is first realistic target.
- Q5382 (dm3): ★ AP-BDH sufficiency. Mesoscopic root table. Z≤√p from linear root-sum.
- Q5383 (dm4): ★ Poisson monodromy. Hyperoctahedral model. Mellin anti-concentration.

### ChatGPT Q&A processed (Q5385-Q5387)
- Q5385 (dm2): ★ Bombieri-Fouvry segmented Weil FAILS (conductor O(H) in block H).
  Heath-Brown/Buchstab NO new gain. Function field literal analogue DEGENERATE.
  KEY: identified exact missing lemma = "Apéry shift-correlation target":
  |Σ_r e_p(a(b_{r+s}-b_r))| ≪ p^{1-η}. Autonomous birational map F(x,u,v) has
  slow O(s) degree growth. Van der Corput step → hybrid algebraic-dynamical sum.
- Q5386 (dm3): ★ Padé/G-function gives UPPER bounds only (Fischler-Rivoal).
  GCD conjecture = denominator sharpness + residual coprimality (G_n = (L_n³/C_n)·gcd(A_n,b_n)).
  Nesterenko one-way (uses GCD, doesn't bound it). ζ(2) analogue: NO published all-n proof.
  Zudilin p-adic: most concrete route via Frobenius matrix factorization, but even
  exact local formula needs global cross-prime support theorem.
- Q5387 (dm4): ★★ EXACT random model for Poisson(1/2): uniform f ∈ V_p^+ (inversion-
  invariant functions) gives K ~ Bin((p-3)/2, 1/p) → Poisson(1/2) RIGOROUSLY.
  Gap: H_p = Δ^ε B_p² is deterministic, not random in V_p^+. Named the missing
  theorem: "Apéry crystalline Mellin anti-concentration." Katz character-aspect
  is genuine but controls ℓ-adic angles, not p-adic valuation. Function field
  needs different geometric variable.

## ★ BREAKTHROUGH 6: Shift-correlation random + exact decomposition + random model

### Exact M₂ decomposition (N=64..8192)
Diagonal ≈ 63%, off-diagonal ≈ 37% of M₂, consistently.
**CRT error is NEGATIVE** for N ≥ 256: the actual off-diagonal is SMALLER than
the CRT independence prediction. Zero sets have mild negative cross-prime correlation.
|CRT error|/N = O(1), negligible vs diagonal = Θ(N log²N).

### Shift-correlation: RANDOM (all 666 primes ≤ 5000)
C_p(a) = Σ_r e_p(a·(b_{r+1} - b_r)) satisfies:
- max|C_p(a)|/√p: mean 3.35, max 5.07
- max|C_p(a)|/(√p · √(log p)): mean 1.23, max 1.80
- EXACTLY consistent with random differences (Gumbel law for √(2 log p))
- All shifts s=1,2,3,5 give same behavior
- Histogram of differences: std/mean ≈ 1 (Poisson!), zeros ≈ p/e (= 1/e fraction)
**Bottom line:** b_r mod p is INDISTINGUISHABLE from random in F_p.

### Exact random model (Proposition in paper, from Q5387)
For uniform f ∈ V_p^+ (space of palindromic functions F_p* → F_p):
K ~ Bin((p-3)/2, 1/p) → Poisson(1/2)
Proof: DFT is invertible, reflected spectral pairs are independent uniform.
Gap to Apéry: H_p has algebraic structure (square of anti-palindromic B_p).

### Paper updated (27 pages, committed e85fd3f)
Added: N=8192 row, exact decomposition table, random model proposition,
shift-correlation computation subsection.

### ChatGPT Q&A processed (Q5388-Q5390)
- Q5388 (dm2): ★ PGL₂ generation PROVED: T(1),T(2),T(3) generate PGL₂(F_p) for all p ≥ 7.
  Proof: T(1)T(2)^{-1} = unipotent, T(1)T(3)^{-1} = different unipotent with coprime steps,
  these generate B⁺, plus det(T(k)) = k³ gives non-square semisimple → full PGL₂.
  BUT: prefix ≠ orbit. No existing dynamics theorem (Helfgott, Bourgain-Gamburd, Breuillard-
  Green-Tao) gives fixed-seed power saving. The "group-growth program" (Steps 1-3) is
  identified but requires genuinely new input: spectral gap → prefix equidistribution
  for the SPECIFIC seed B_p, not a generic seed.

- Q5389 (dm3): ★★ **5/8 CONSTANT DERIVED.** Clean formula C₂(μ) = 3μ/8 + μ²/4.
  For μ = E[Z(p)] ≈ 1: C₂(1) = 3/8 + 1/4 = 5/8. Diagonal 3/8 (from E[Z(Z-1)] + E[Z]),
  off-diagonal 1/4 (from E[Z]²). Observed ≈ 0.59 converges from below (boundary effects).
  **Boundary correction:** using L_{p,q} (admissible interval length) instead of N in CRT
  main term, the genuine residual is POSITIVE and O(√N (log N)²). The "negative CRT error"
  from Breakthrough 6 was an artifact of using N for prime pairs near the boundaries.
  **V° decomposition:** V° = Σ w_p² A_p(1-A_p/N) + E°_N where E°_N = centered cross-prime
  covariance. Need V° ≤ C·S (and (Z1) + (CRT2)) for the theorem.

- Q5390 (dm4): ★★★ **RANDOM-SQUARE POISSON(1/2):** For uniform g ∈ V_p^- (anti-invariant),
  f = g² gives K ~ Poisson(1/2) RIGOROUSLY. Proof: quadratic Gauss sums + Mattson-Solomon
  uncertainty inequality. Discriminant factor Δ^ε doesn't change the law.
  **TRIANGULAR POISSON(1/4):** First ~p/4 Mellin coefficients of Δ^ε B² are triangular-
  linear in free β_k, hence independent uniform. Gives Poisson(1/4) for first half of
  reflected pairs. Full 1/2 needs complete-intersection theorem for the QUADRATIC range.
  **Shift-correlation ≠ zero anti-concentration (KEY):** Changing z entries changes C_s(a)
  by O(z). Even perfect √p cancellation gives Z(p) = 1 + O(√p), not O(1). The objects are
  genuinely different: C_s is an additive character sum, Z counts solutions.

## ★ BREAKTHROUGH 7: Centered variance, 5/8 constant, boundary correction (computation + Q5389-Q5390)

### V°/S → 1 — COMPUTED (variance_computation.py)
Centered variance V° = M₂ - M₁²/N computed directly for N = 64..4096:
| N    | V°/S   |
|------|--------|
| 256  | 1.014  |
| 512  | 1.010  |
| 1024 | 0.981  |
| 2048 | 0.992  |
| 4096 | 0.995  |
This IS the exact prediction of the independent Bernoulli model (V° = S when events
are independent). The centered cross-prime covariance E°_N oscillates around 0 with
|E°_N| = O(√N (log N)²). **Apéry zero-sets are essentially pairwise uncorrelated.**

### 5/8 constant — DERIVED (Q5389)
M₂(N) / (N log²N) → 5/8 = 0.625. Formula: C₂(μ) = 3μ/8 + μ²/4.
Diagonal contribution 3/8 (from Var[Z] + E[Z] = μ/2 + μ per prime).
Off-diagonal contribution 1/4 (from E[Z]² = μ² per prime pair).
Observed ≈ 0.59 converges from below due to boundary effects.

### "Negative CRT error" CORRECTED — was boundary artifact
Previous finding (Breakthrough 6): "CRT error is negative for N ≥ 256".
Correction: the CRT main term used N instead of L_{p,q} (the actual admissible
interval for each prime pair, accounting for p ≤ m < p² constraint). After boundary
correction, the genuine residual is POSITIVE and O(√N (log N)²). No arithmetic
repulsion — the data is consistent with random.

### Random-square Poisson — PROVED (Proposition prop:random-sq in paper)
For uniform g ∈ V_p^-, f = g² gives Poisson(1/2). Added to paper §12.

### Paper updated (28 pages, committed db75857)
Added to §12: 5/8 constant derivation, V°/S table, boundary correction remark,
random-square Poisson proposition.

### Currently dispatched (3 dm tabs processing)
- dm2: Q5393 — Rainbow estimate / one-point bound approaches
- dm3: Q5391 — Variance convolution structure / cross-prime covariance
- dm4: Q5392 — Triangular Poisson complete-intersection / second 1/4
- dm1: DEAD (both Q5384 and status check failed — tab offline)

### Three clean remaining targets
1. **V° ≤ C·S (cross-prime near-independence):** V°/S → 1 empirically. Need a theorem.
   Routes: CRT independence via equidistribution, or random-model coupling.
2. **First moment (Z1): Σ Z(p) ≪ π(x):** Bounded average number of zeros per prime.
   Routes: counting (p,r) pairs with p | b_r, triangular Poisson(1/4), growth rate of b_r.
3. **Rainbow estimate: B(n) = o(n/log n) for all n:** The one-point bound.
   Routes: Waiting for Q5393. Integrality of b_r, recurrence constraints, sieve.

## Summary of results (UPDATED)
19. ★ |Z_p| ~ 2·Poisson(1/2) — verified all primes ≤ 10000
20. ★ Second-moment criterion — M₂=O(N^{2-δ}) → W(n)=o(n)
21. ★ M₂ ≈ 0.59·N·(logN)² → 5/8 constant (diagonal 3/8 + off-diagonal 1/4)
22. ★ No-go model — vertical bounds alone insufficient
23. ★ AP-BDH sufficiency — dispersion + Z≤p^{2/3} → N^{5/6+o(1)}
24. ★ Classical tools ALL exhausted — confirmed by 7+ ChatGPT audits
25. ★ Gap polynomial sparsity — Σ r_p(h) ≤ 1 for h ≤ 100
26. ★★ **Shift-correlation RANDOM** — max|C_p|/√p ≈ √(2logp), all p ≤ 5000
27. ★★ **Exact random model** — Bin((p-3)/2, 1/p) → Poisson(1/2) for V_p^+
28. ★★ **CRT error = boundary artifact** — after correction, positive and O(√N log²N)
29. ★ **Exact decomposition** — Off/Diag ≈ 0.58, |CRT err|/N = O(1)
30. ★★★ **V°/S → 1** — Apéry zeros pairwise uncorrelated across primes
31. ★★★ **5/8 constant** — C₂(μ) = 3μ/8 + μ²/4, exact random prediction
32. ★★ **Random-square Poisson(1/2)** — rigorous for V_p^- model (prop in paper)
33. ★★ **Triangular Poisson(1/4)** — rigorous for first half of spectrum
34. ★★ **PGL₂ generation** — T(1),T(2),T(3) → PGL₂(F_p), p ≥ 7

## ★ BREAKTHROUGH 8: L²-to-uniform upgrade + exact covariance + martingale obstruction (Q5405-Q5407)

### L²-to-uniform lemma — THE bridge to all-n
Lemma: max|B(m) - μ| ≤ √V°. (Trivial but decisive.)
Corollary: V° ≤ C·S → W(n) = o(n) for ALL n.
Proof: √V° ≤ √(CS) = O(√N (log N)^{3/2}) = o(N/log N).
No concentration inequality, NA, or fourth moment needed (Q5407).

### Martingale predictive cascade — why sequential conditioning fails (Q5406)
Doob martingale M_k = E[B | X_{p_1},...,X_{p_k}] has increments D_k = (X_{p_k} - π_k)Δ_k.
The "cascade" Δ_k = w_{p_k} + Σ_{j>k} w_{p_j} (conditional prob shifts) — NOT just w_{p_k}.
Under independence, Δ_k = w_{p_k} and Var(B) = S. But for Apéry on partial period,
revealing one prime nearly identifies m (when A_p ≈ 1), making future indicators deterministic.
Controlling Σ E[π_k(1-π_k)Δ_k²] ≤ (1+o(1))Σ w_k² α_k is EQUIVALENT to dispersion.

### 96.7% negativity = sparsity, NOT negative dependence (Q5405)
C_{p,q} < 0 iff J_{p,q} = 0 (no common hit). In sparse regime (λ < 1), Pr(J=0) ≈ e^{-λ} ≈ 96%.
The rare J=1 pairs have C = 1-λ ≈ 1 (large positive). Over complete period, mean C = 0 exactly.
Dubhashi-Ranjan (balls-into-bins) does NOT apply. Janson requires nonneg intersection params.
The relevant target: aggregate |E°| = o(S), i.e., cancellation in the SIGNED sum.

### Exact covariance data (no sampling)
| N    | K   | V°/S   | E°      | |E°|/S  |
|------|-----|--------|---------|---------|
| 512  | 67  | 0.998  | +130.6  | 0.015   |
| 1024 | 113 | 0.998  | +215.6  | 0.010   |
| 2048 | 205 | 0.994  | +96.2   | 0.002   |
| 4096 | 380 | 0.996  | +120.0  | 0.001   |
| 8192 | 742 | 1.001  | +1344   | 0.005   |

E° is POSITIVE at all tested N (except N=2048 with sampling artifact corrected to +96).
V°/S ≈ 1 ± 0.005. The sub-Poisson bound E° ≤ (S-S_c) fails at N=8192 (E°=1344 > target=1163).
CS bound is 3000-6000× too loose — irrelevant.

### Key conclusion
The ENTIRE remaining challenge: prove V° ≤ C·S for SOME constant C.
- Not V° ≤ S (sub-Poisson, fails at N=8192)
- Not |E°| = o(S) (this is the precise target, but any C works)
- The paper now has the formal conditional theorem + full computational evidence

### Paper updated: 31 pages, committed 5149dd2
Added: Lemma lem:l2-uniform, Corollary cor:dispersion-alln, Table tab:covariance,
Remark rem:martingale. Corrected rem:negative-crt. Updated rem:concentration.

### ChatGPT Q&A processed (Q5405-Q5407)
- Q5405 (dm1): ★ Sparsity explanation of 96.7% negativity. Mean C = 0 over complete period.
  Sub-Poisson collision is the correct target. Janson/DR do NOT apply.
- Q5406 (dm2): ★★ Martingale cascade. Δ_k = w_k + cascade. Independence kills cascade.
  Controlling cascade ≡ dispersion. No shortcut via conditioning.
- Q5407 (dm3): ★★★ L²-to-uniform lemma. V° ≤ CS → all-n. Formal proof. Paper-ready LaTeX.
  Fourth moment, sub-Gaussian, NA all unnecessary. k=2 is weakest and most accessible case.

### Currently dispatched (4 dm tabs)
- dm1: Q5412 — Random CRT formalization (Fourier, BDH, Gallagher, spectral)
- dm2: Q5417 — Palindromic structure forcing negative covariance
- dm3: Q5418 — Elliott-Halberstam / BDH connection, completion of sums
- dm4: processing (earlier dispatch)

### Three remaining avenues for proving V° ≤ CS
1. **Signed pair-dispersion** (the "aggregate sub-Poisson collision theorem"):
   Show cancellation in E° = Σ w_p w_q C_{p,q} via the sign structure.
   The data shows |E°|/S = O(1/√N), but proving this requires arithmetic input.
2. **Completion of sums**: Decompose [N+1,2N] into pq-periodic pieces + remainder.
   For complete periods, C_{p,q} = 0. The remainder contributes O(Z_p Z_q).
   Need to show the remainder terms cancel across pairs.
3. **Random CRT model formalization**: If CRT discrepancies behave like random
   variables with magnitude √(N/(pq)) and random signs, then E° has standard
   deviation O(√N (log N)²) while S = O(N log N). Formalizing "random signs" requires
   quantitative equidistribution of zero-set locations across primes.

## Summary of results (UPDATED)
35. ★★★ **L²-to-uniform upgrade** — V°≤CS → all-n via Chebyshev
36. ★★ **Martingale obstruction** — cascade ≡ dispersion, no shortcut
37. ★★ **Sparsity = negativity** — 96.7% negative is structural, not NA
38. ★ **Exact covariance table** — V°/S ∈ [0.994, 1.001] for N ≤ 8192
39. ★ **CS bound useless** — off by factor 3000-6000

## Resources
- ChatGPT Pro: dm1,dm2,dm3,dm4 (SOL, all processing).
- uisai2: Z(p) computation COMPLETED. Max Z(p)=12 at p=159977 through p~1.3M.

## Session 2026-07-19 (dm window, relayed from cron)

### Cross-host harvest (IMPORTANT coordination note)
Another host ran a fresh Q84–Q119 campaign TODAY (2026-07-19, Notion drops,
separate ledger numbering). Harvested the 10 Apéry-relevant answers into
chatgpt-answers/ (Q89, Q90, Q92, Q94, Q96, Q97, Q101, Q117, Q118, Q119).
That host has UNPUSHED chain text ("bordered certificate" rewrite) — do not
duplicate its (4.E) repair work. This machine's bridge: server up, ZERO tabs
registered → ChatGPT dispatch unavailable here until tabs reconnect.
Ledger Q5510/Q5511/Q5512/Q5515 (07-16): GITDROP SHAs unresolvable in ANY repo
(hallucinated commits, [reported]-trap); treat as lost — superseded by Q84-Q119.

### Q117 (adversarial audit of the parallel session's chain)
- (4.A)–(4.D) VALID, including TWO NEW unconditional theorems:
  ★ value-fiber bound max_a N_p(a) ≪ p^{3/4} (bordered 2×2 determinant
  D_{h,k}(t) = N_h(Π_k − B_k) − N_k(Π_h − B_h), a-independent certificate)
  ★ collision energy E(p) = Σ_a N_p(a)² ≪ p^{7/4}
- (4.E) FATAL for THAT chain: support law omits companion-block channel
  (counterexample 13 | G_26 with b_0 = 1). OUR paper is unaffected —
  prop:companion-height already handles that channel at every index.
- Repair list (§10): correct depth v_p ≤ 6 (that chain), dyadic Markov
  localization, companion channel estimate.

### Q119 (strategic verdict on the pointwise gap)
- ★ IMPOSSIBILITY THEOREM: no proof from |Z_p| ≪ p^{2/3} + reflection +
  no-consecutive + n < X² alone — explicit adversarial S_p = {n₀ mod p,
  p−1−n₀ mod p} hits every prime. Cross-prime arithmetic coupling of the
  ACTUAL Apéry sets is mandatory. (Matches our no-go model, paper §10.4.)
- Target A: fixed high moment (HM)_k, k ≥ 7 (k = 8):
  Σ_{m<X²} (K_X(m))_k ≪ X^{2+o(1)} λ_X^k  ⟹  max K ≪ X^{2/3+2/k+o(1)}
  ⟹ W(n) ≪ n^{11/12+o(1)} with existing Z ≪ p^{2/3}. Sharpest concrete target.
- Target B: small-quotient shifted divisibility, q ≤ log²n only.
- Harmless-range observation (now formalized in paper as
  prop:quotient-reduction): primes √n < p ≤ n/(f(n)log n) are killed by
  prime counting alone.
- ⚠ Q119 §6 IS WRONG (poisoned premise): the "unconditional
  log G_n ≤ 2.2467n via Rhin–Viola" uses q_n = b_n/G_n from the OUTGOING
  question; the true reduced denominator is q_n = d_n b_n/G_n, which gives
  only log G_n ≤ 5.25n — WORSE than the trivial Wronskian bound 3n.
  Do NOT propagate 2.2467n. (Verified independently this session.)
  Lesson instance: pre-verify premises embedded in outgoing prompts.

### My verification of the (HM)_k structure
- (HM)_2 holds trivially and unconditionally: each pair (p,q) has
  ≤ 1 + X²/pq CRT representatives below X², and (Σ Z_p)² ≍ X²λ_X².
- For k ≥ 3, Πp_i > X² ⟹ ≤ 1 representative per tuple; (HM)_k becomes a
  genuine equidistribution statement (fraction X²/Πp of ~(Xλ)^k tuples).
  k > 6 needed for the max-bound to beat trivial. All consistent with Q119.
- Windowed-dispersion sufficiency (localized V° ≤ CS at scale-X windows ⟹
  W(n) ≪ n^{5/6+o(1)}) is the window-local form of the paper's existing
  AP-BDH Hypothesis 12 — same funnel: Fourier uniformity/cross-prime
  equidistribution of actual zero sets is THE missing input.

### Paper changes (this session)
- NEW prop:quotient-reduction + rem:bounded-quotients: the conjecture is
  EQUIVALENT to the lower-digit channel bound uniformly over quotient
  classes q < f(n)log n (any f → ∞); primes below n/(f log n) are free.
  Sharpened rem:rainbow-open accordingly. 52 pages, compiles clean.

### Codex dispatched (this session)
- C1 (xhigh): CODEX_SPEC_fiber.md — port the value-fiber p^{3/4} theorem +
  energy corollary as fiber_bound.tex + fiber_verify.py.
- C2: CODEX_SPEC_bigcompute.md — B(n) scan to 10^6–2×10^6 via (p,r) scatter
  from per-prime zero sets + windowed dispersion statistics + Z(p) histogram
  cross-check. Headline: does B(n) ≤ 3 persist?

### Current priority ranking
1. (HM)_8 / Target A — the cleanest formulated missing theorem.
2. Uniform-in-q top-half machinery (Target B for q < f log n) — now the
   EXACT open statement after prop:quotient-reduction.
3. Fiber/energy port (Codex C1) — strengthens unconditional toolbox.
4. Computation (Codex C2) — Poisson consistency at 10× the old range.
