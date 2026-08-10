# DOCTRINE — P3.2 Full Unconditional Proof (Updated 2026-08-09)

## Goal
Prove G_n = e^{o(n)} for ALL n, unconditionally.

## The Collapse (Fable R5, 2026-08-09)

One-line reduction that replaces the entire (HM) tower:

  Σ (K)_k ≤ (max K)^{k-2} · Σ (K)_2 ≤ 5X²λ² · (max K)^{k-2}

Therefore:

**(AT″): max K_X(m) ≪ λ_X · X^{o(1)} ⟹ (HM)_k for ALL k ⟹ full conjecture.**

Numerics: max K ~ log X (Poisson), λ ~ 0.08. Ratio max K/(Xλ) → 0. Massive slack.
Star-falsified: anchored star has max K ~ X/log X ≫ λ·X^{o(1)}. ✓

## Avenues (ranked after 5 rounds of oracle discussion)

### (a) APEX: Prove (AT″) — max K ≪ λ · X^{o(1)}

Sub-routes:
- (a1) **Twin-atom lemma** → max K ≪ X^{1/2} (G2). Uses: codegree/gap-polynomial
  for no-twin-atoms at distance < X, then variance localization.
- (a2) **Mellin horizontal twist** — u^{-m} is a FIXED character across all p.
  By Mellin identity, Σ N(u)u^{-m} ≡ -b_{m mod(p-1)}. Only frame where
  "one fixed object, many p" is literally true → horizontal ordinarity.
- (a3) **Shell restriction**: restrict thm:hm-pointwise to m ∈ (X²/2, X²].
  Eliminates small-index divisor barrier. Verify this survives the proof.

Terminal: (AT″) proved ⟹ done. All sub-routes exhausted ⟹ fall to (b).

### (b) CED ensemble: (MC) + (AT) ⟹ (HM)_3

From Fable R4: (HM)_3 ⟸ (MC) + (AT), both star-falsified.
- (MC): M_p(k,k') as line sections of {N_h(r)=0} cloud. Gap = p^{1/6}.
- (AT): max K ≪ X^{2/3+o(1)}λ. Weaker than (AT″).

Terminal: both proved ⟹ (HM)_3. One has proof-of-failure ⟹ fall to (c).

### (c) Vertical: avg Z(p) ≪ p^{3/5} + (HM)_6

First Lemma (S_{d,r} ≠ 0). Corank-valuation trick. Gives c = 3/5, needs k = 6.
Low priority since (AT″) subsumes.

### (d) Function field: Mason-Stothers over F_q(t)

Time-boxed one design pass. Structural insight only.

## Landed pieces (verified)

1. **Theorem A**: b_r ≡ [x^r y^r z^r w^r] F^{p-1} (Straub+Cartier) ✅
2. Z(p) = O(p^{2/3}) ✅
3. (HM)_2 ✅
4. Exact orthogonality off-diagonal ✅
5. Lemma 1 (Fourier non-concentration) ✅
6. M_p(k,k') = Z·δ for palindromic Z_p ✅
7. (AT″) numerics: ratio → 0 ✅
8. Mellin identity for b_r ✅
9. Honest reduction: (HM)_3 ⟸ (MC) + (AT) (Fable R4)
10. Collapse: (AT″) ⟹ all (HM)_k (Fable R5)

## Fallback

Paper has density-1 + gap analysis. New results to add regardless:
Theorem A, Lemma 1, the (AT″) collapse, the (HM)_3 reduction.

## Active Lean run (2026-08-09, tmux prompt 12:52:29)

### Goal

Continue the faithful formalization in
`lean/RamanujanChallenge/Problem32/Main.lean`: inspect its three real `sorry`
declarations and close the easiest next dependency-frontier theorem without
adding axioms or strengthening statements beyond the mathematics proved in
`proof.tex`.  The immediate context is the completed aligned-corank and
mesoscopic-census work.

### Ranked avenues

1. **(a) Dependency-frontier closure.** Audit the actual theorem dependency
   graph and attack the earliest remaining `sorry` whose paper proof and Lean
   prerequisites are already present, starting with consequences of aligned
   corank and the mesoscopic census.  Search the repository and Mathlib before
   building new infrastructure.
2. **(b) Faithful arithmetic bridge.** If the frontier lemma is missing a
   genuine bridge, isolate the exact statement from `proof.tex`, prove it as a
   documented reusable lemma, and reconnect the existing downstream theorem.
3. **(c) Obstruction theorem / narrower reduction.** If a target cannot follow
   from the encoded hypotheses, prove a checked counterexample or reduce it to
   one strictly narrower explicit arithmetic interface; do not hide the gap in
   a carried hypothesis.

### Terminal conditions

- **Success:** the selected frontier `sorry` is gone, the single-file Lean
  check passes with one fewer warning, and its downstream theorem is
  axiom-audited.
- **Avenue failure:** a checked counterexample, dependency-cycle proof, or
  precise reduction shows the avenue cannot close the next frontier theorem;
  record it and continue with the next ranked avenue.
- **Hard stop:** only a destructive/external action outside the repository or
  genuinely missing owner-only information blocks further work.

## Active formalization campaign (2026-08-09, `/lean` + `/automode`)

### Main goal

Build a faithful, modular Lean formalization of the proved results in the
current P3.2 manuscript, with every conjectural atom-tail input kept visibly
conditional, beginning from the two remaining `sorry`s in
`Problem32/Main.lean`.

### Ranked avenues

1. **(a) Unconditional zero-count chain.** Formalize the paper's gap
   polynomials `N_h`, their exact degree/leading coefficient, the evaluations
   at `-1` and `-2`, nonvanishing modulo `p`, the polynomial root bound, and
   the block-partition optimization.  Success means
   `zero_count_sublinear` is sorry-free and axiom-clean.  Proof-of-failure
   means a kernel-checked counterexample to the Lean statement or a verified
   mismatch with the manuscript's constants; then correct the statement and
   continue in the same avenue.
2. **(b) Cartier--Frobenius structural layer.** Formalize the exact Lucas
   digit law, block law, and valuation/GCD consequences that are unconditional
   and independent of the atom-tail conjecture.  Success means a documented,
   axiom-clean module whose public declarations correspond line-by-line to the
   manuscript section.  Failure requires a checked false source atom, which is
   retained as an explicit source-gap/corrected pair.
3. **(c) Exceptional-set theorem.** Audit the exact current manuscript
   statement against `problem32_polylog_exceptional`, then formalize the
   dyadic counting and prime-sum chain from the preceding modules.  Success is
   a faithful sorry-free theorem; a stale or stronger Lean statement is
   corrected rather than silently assumed.
4. **(d) Conditional atom-tail/HM assembly.** Encode the atom-tail and
   factorial-moment interfaces and prove every downstream conditional
   implication, without treating those hypotheses as unconditional results.
   Success is an honest conditional capstone plus a precise residual list.

### Fallback and global terminal conditions

If a headline avenue reaches a checked source gap, continue section-by-section
on independent proved lemmas so that the formalized prefix keeps growing.  The
campaign succeeds when all proved manuscript claims in scope are represented
by faithful, axiom-clean Lean declarations and the remaining conjectural
inputs are exactly identified; it fails only through checked false statements
or exhausted avenues recorded with concrete terminal evidence.
