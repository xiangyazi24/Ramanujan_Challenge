# P3.2 Campaign Understanding (2026-08-09, automode session)

## Goal
Prove G_n = e^{o(n)} for ALL n (full Apéry GCD conjecture).

## The (AT″) Collapse — the APEX insight

One-line:  Σ(K)_k ≤ (max K)^{k-2} · 5X²λ²

Therefore: **max K ≪ λ·X^{o(1)} ⟹ full conjecture**. One pointwise statement buys everything.

Numerics: max K = 4 at X = 4096, ratio max K/(Xλ) = 0.013 → 0. Poisson extreme value.

## Landed theorems (verified, ready to write up)

1. **Theorem A** (Cartier diagonal): b_r ≡ diag(F^{p-1}) mod p, F = Straub quartic
2. **Lemma 1** (Fourier non-concentration): Σ|F_p(k)|² ≤ C(K₀Z + p²/K₀)  
3. **Twin-atom lemma** (repaired codegree): W(m,m') ≤ 6h + O(h²/log X), W(m,m+1)=0
4. **Triple bound**: W₃ = 0 for gaps ≤ 10 (computationally), O(1) theoretically via resultant
5. **(AT″) collapse**: one-line moment inequality → full (HM)_k tower
6. **Honest (HM)_3 reduction**: (HM)_3 ⟸ (MC) + (AT), both star-falsified
7. **Exact orthogonality**: off-diagonal v-sum = Z(p)Z(p')
8. **Reflection spray**: T-atom forces T distinct hits within 4X below

## Dead avenues (documented)

- (a2) Mellin horizontal twist: |S̃_p(χ)| shows continuum, kill criterion triggered
- Sidon property: trivially true for doublets, implies Z ≤ √p (too strong to prove)
- Direct pointwise (Family E): blocked (holonomic complexity)
- ABC over Q: vacuously true
- Chebotarev for fixed twist: fails in principle (weight 3 + same-char wall)

## Current frontier

The full conjecture reduces to (AT″): max K_X(m) ≪ λ·X^{o(1)}.

Three attack surfaces remain:
1. **m-side pair dispersion** using twin-atom + triple bound (Fable R10 designing)
2. **First Lemma** (S_{d,r} ≠ 0): would unlock average Z(p) ≪ p^{3/5} + triple bound unconditionally
3. **Function field**: Mason-Stothers (untried, time-boxed)

## Fable R10 Definitive Answer

Every audited route terminates at atom tail: n_T ≪ X^{2+o(1)}/T³.
One X^ε derandomization is the whole problem.

Provable theorems for paper: T1 (short-range decorrelation), T2 (per-prime 2-local), T3 (codegree repair).

Hot Fable agent: a9c82b2ecf78808da (resume via SendMessage, has R1-R10 context).

## Key files

- proof.tex: main paper (6372 lines)
- DOCTRINE.md: automode doctrine
- RUN_LOG_P32.md: session log
- ORACLE_COMM/: Codex communication + results
- WZCertificate.lean: Lean proof skeleton (sorry #1, deferred)
