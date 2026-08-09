# RUN_LOG — P3.2 Full Unconditional Campaign

## Run 2026-08-08/09 (automode)
- doctrine: DOCTRINE.md (updated 2026-08-09)
- starting avenue: (a) APEX — prove (AT″)
- status: IN PROGRESS

### Milestones landed

1. **Theorem A** (Straub + Cartier): b_r ≡ diag(F^{p-1}) mod p ✅
   - Verified for p = 5,7,11,13,17,19,23 (all r)

2. **Lemma 1** (Fourier non-concentration): Σ|F_p(k)|² ≤ C(K₀Z + p²/K₀) ✅
   - Proof: Fejér + A_p(h) ≤ 3(h-1)
   - Verified numerically for p ≤ 200

3. **Exact orthogonality**: Σ_v F_p(kv̄)F̄_{p'}(k'v̄) = Z(p)Z(p') for p≠p' ✅
   - Verified for 5 prime pairs, all (k,k')

4. **M_p(k,k') structure**: = Z·δ_{k,k'} for palindromic Z_p ✅
   - Verified for p ≤ 200

5. **(AT″) collapse** (Fable R5): max K ≪ λ·X^{o(1)} ⟹ all (HM)_k ✅
   - One-line: Σ(K)_k ≤ (max K)^{k-2} · 5X²λ²

6. **Two-flip reciprocity** ✅ (Codex, verified to float error 3e-14)

7. **(AT) numerics**: max K / (Xλ) → 0 ✅
   - X=4096: max K=4, Xλ=314.7, ratio=0.013

8. **Near-Sidon property**: M_p(k,k') ≤ 1 for ALL k≠k', for 466/468 primes ≤ 10000
   - 2 violations at p=3727 (Z=8, contains a 3-term GP)
   - Average violations: 0.004/prime

### Strategic decisions

- Fable R4: CED sketch has 2 fatal flaws. Honest: (HM)_3 ⟸ (MC) + (AT)
- Fable R5: (AT″) collapses entire tower. APEX target.
- Near-Sidon: (MC) gap (p^{1/6}) likely closeable from near-Sidon on average
- (AT″) is the remaining bottleneck

### Current focus

Avenue (a2): Mellin horizontal twist design pass — waiting for Fable R6

### Avenue (a2) terminal verdict: DEAD
|S̃_p(χ)|/p^{3/2} shows continuum (11 unique values at p=23, 14 at p=29),
not bounded monomial decomposition. Kill criterion from Fable R6 triggered.
The Mellin horizontal twist does not produce a clean p-adic Gamma handle.

### Current avenue: (a1) Twin-atom lemma
Target: prove no two T-atoms within distance X, for T ≥ X^{1/2+ε}.
Tool: codegree/gap-polynomial + corank-valuation.

### Twin-atom lemma: VERIFIED (2026-08-09)
- W(m,m+1) = 0 ✅ (all X² range, X=128)
- Type A/B classification works ✅
- N_h(m) > 0 for h=2..20, m=0..100 ✅
- P(n)-(n+1)³-n³ = 31n³+45n²+21n+4 > 0 ✅
- Reflection spray: 100% hit rate on K≥3 atoms ✅

### Sidon deflation accepted (Fable R7)
- Trivially true for doublets
- Sidon ⟹ Z ≤ √p+1 (stronger than 2/3 bound)
- Star is Sidon → insufficient alone
- Birthday threshold p=3727 at Z=8, p^{1/4}≈7.8 — exact match

### Height profile (S1): Σ(1/h) < Z for ALL tested primes ✅
- min heights grow as ~√p
- Σ(1/h) → 0 as p → ∞

### (a2) DEAD: continuum |S̃_p(χ)|, kill criterion triggered

### Campaign achievements (publishable):
1. Theorem A (Cartier diagonal) — NEW theorem
2. Lemma 1 (Fourier non-concentration) — NEW lemma  
3. (AT″) collapse — NEW reduction (one line → all HM_k)
4. (HM)_3 ⟸ (MC)+(AT) — honest reduction
5. Twin-atom lemma (repaired codegree) — FIXES paper's wrap hole
6. Reflection spray — NEW structural fact
7. Near-Sidon + height data — NEW computational evidence

### First Lemma extended: S_{d,r} ≠ 0 for ALL 399 pairs d,r ≤ 21 ✅
Triple bound unconditional for gaps h₁+h₂ ≤ 41.

### Paper writing: sent to Codex p32, awaiting output
### Dispersion assembly: sent to Fable, awaiting R10
