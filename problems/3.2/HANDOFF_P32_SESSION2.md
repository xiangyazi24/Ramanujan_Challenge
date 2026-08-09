# HANDOFF — P3.2 Full Unconditional Campaign (Session 2, 2026-08-09)

automode: yes

## Goal
Prove G_n = e^{o(n)} for ALL n (full Apéry GCD conjecture).

## What was accomplished this session

### Paper contributions (committed, ~1000 lines LaTeX)
- `new_sections.tex` (646 lines, Codex) — Theorem A, Lemma 1, moment collapse, (HM)_3 reduction
- `atom_tail_section.tex` (337 lines, Fable R12) — definitive atom-tail section with 7 theorems

### Key mathematical results
1. **Theorem A**: b_r ≡ diag(F^{p-1}) mod p (Straub + Cartier). VERIFIED for p ≤ 23.
2. **(AT″) collapse**: max K ≪ λ·X^{o(1)} ⟹ all (HM)_k ⟹ full conjecture.
3. **Twin-atom lemma** (codegree repair): W(m,m') ≤ 7h, W(m,m+1) = 0. VERIFIED.
4. **Triple bound**: W₃ = 0 for small gaps. S_{d,r} ≠ 0 verified for 399 pairs (d,r ≤ 21).
5. **Concentration theorem**: full conjecture = atom tail n_T ≪ X²/T³.

### Dead avenues (ALL documented with precise terminal verdicts)
- (a2) Mellin twist: |S̃_p(χ)| continuum, kill criterion
- Sidon: trivially true for doublets, implies Z ≤ √p
- Chebotarev: weight 3 + same-char wall
- ABC: vacuous over Q
- Dispersion diagonal: subsumed by atom problem
- Signature energy: atoms are Poisson-pure, no structure
- **Function field F_q(t): THREE kill shots** — Lucas collapses to fixed list, d_n undefined, monodromy equidistributes wrong direction

## CURRENT FRONTIER: ℚ(λ)-pencil (Fable R13, brand new)

The **last unexplored direction** from this campaign:

- b_r(λ) = Σ C(r,k)²C(r+k,k)²λ^k ∈ ℤ[λ], deg = r
- Z_p^{(λ)} = {(r,λ) : b_r(λ) ≡ 0 mod p}. Our Z_p is the fiber at λ = 1.
- **New resultant family**: R_{r,r'} = Res_λ(b_r(λ), b_{r'}(λ)) — TRANSVERSE to gap polynomials
- Endgame: prove atom bound for all λ outside algebraic exceptional locus, verify λ = 1 not exceptional

### THREE cheap decisive computations (DO THESE FIRST):
1. **R_{r,r'} recon**: Are b_r(λ) pairwise coprime over Q? Compute Res_λ(b_r, b_{r'}) for small r,r'. If R_{r,r'} ≠ 0 → new arithmetic constraint on cross-prime divisibility.
2. **Fiber statistics**: Compute K^{(λ)}-statistics for λ = 2, 3, 1/2, -1 alongside λ = 1. Is λ = 1 statistically generic?
3. **F_p factorization**: b_n(t) mod p for p = 5, 7 — verify Lucas collapse (K1).

## Hot Fable oracle
Agent ID: a9c82b2ecf78808da (resume via SendMessage)
Has FULL campaign context (R1-R13). Use for strategic guidance.

## Key files
- problems/3.2/DOCTRINE.md — automode doctrine (updated)
- problems/3.2/RUN_LOG_P32.md — session log
- problems/3.2/UNDERSTANDING_P32.md — full state
- problems/3.2/ORACLE_COMM/ — all Codex results + tasks
- problems/3.2/atom_tail_section.tex — Fable R12 definitive section
- problems/3.2/new_sections.tex — Codex paper sections

## Instructions for next session
1. Enter /fable-ora mode
2. Execute the 3 ℚ(λ)-pencil computations
3. Use 8 ChatGPT rc tabs for parallel research (cross-pollination)
4. If R_{r,r'} ≠ 0 and λ=1 is generic: pursue the exceptional-locus endgame
5. If R_{r,r'} = 0 for some pair: ESCALATE to Fable (structural discovery)
