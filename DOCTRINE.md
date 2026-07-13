# DOCTRINE — Ramanujan Challenge Automode

## Goal
Solve as many of the 10 Ramanujan Challenge problems as possible by Aug 1 2026.
Two proofs done (2.2, 2.8). Eight remain. All problems attacked — 狭路相逢勇者胜.

## Current state
- **2.2** ✅ SOLVED (Aptekarev, proof PDF)
- **2.8** ✅ SOLVED (Chudnovsky CMF, proof PDF)
- **2.3** ore_algebra running on uisai2 (factorization test)
- **2.1, 2.4, 2.5, 2.6, 2.7, 3.1, 3.2** — attack vectors identified, ChatGPT analyses done

## Deliverable per problem
Each solved problem needs THREE artifacts:
1. **Lean 4 proof** (formal, in Ripple or standalone project) — gold standard
2. **LaTeX PDF** (human-readable proof) — for submission
3. **Python verification** (numerical, 200+ digits) — sanity check

Proofs done so far have PDF + Python but NO Lean. Priority: add Lean proofs,
especially for 2.8 where Ripple already has Chudnovsky infrastructure.

## Avenues (ranked by expected progress)

### (a) Problem 2.4 — weight-4 HPL reduction (MOST MECHANICAL)
Two-stage symbolic summation: inner sum = parameter derivative of ₂F₁,
outer sum = inverse-central-binomial → arcsine integrals → weight-4 HPLs at 1/2.
Terminal: produce a creative telescoping certificate OR a verified CAS derivation.
Tools: mpmath on mini for verification, Sage/Mathematica on uisai2 for telescoping.
Fallback: 500+ digit numerical verification + known-basis decomposition.

### (b) Problem 2.5 — Catalan via Ripple CatalanCertified
Ripple has CatalanCertified.lean (748 lines, 0 sorry). The CMF 3×3 matrix encodes
a summation lift of an order-2 recurrence. Identify the order-2 kernel, match to
Ripple's PIVP formalization.
Terminal: LaTeX proof citing Rivoal-Zudilin + CMF encoding lemma.

### (c) Problem 2.1 — level-5 ₃F₂ identification
Poincaré roots 20φ^{±5}. Parameters (1/2, -4/5, 6/5) suggest ₃F₂. Search for
contiguous ₃F₂(1) CF via Ebisu-Iwasaki or Yamamoto's 38 proofs.
Terminal: identify the hypergeometric identity, write proof.
Fallback: Petkovsek/van Hoeij over Q(√5) on uisai2.

### (d) Problem 2.6/2.7 — ζ(2)+ζ(3) connection formula
GF satisfies explicit ODE. Series converges polynomially (dominant root).
Need connection formula U(1). Route: integral representation or operator factorization.
2.7 is related but more efficient (4-term). Solve 2.6 first, template for 2.7.
Terminal: prove U(1) = ζ(2)+ζ(3)-2077/720.

### (e) Problem 2.3 — π+e Apéry limit
LCLM fails in Q[n]<S>. ore_algebra on uisai2 testing Q(n)<S>. If irreducible,
need full CMF or non-obvious integral representation.
Terminal: factor the operator OR find integral giving π+e.

### (f) Problem 3.2 — Apéry irrationality measure
Computational evidence supports (gcd subexponential). Need p-adic analysis.
Ripple has Apery*.lean infrastructure. Supercongruences key.
Terminal: prove gcd = e^{o(n)} or produce counterexample.

### (g) Problem 3.1 — knot integral for π²
Open conjecture. Khoi 2008 topological route (SL₂(R) representations +
Godbillon-Vey). Need: numerical verification to 500+ digits + literature search
for any post-2008 progress.
Terminal: prove the integral formula or find it proven in literature.

## Fallback
If an avenue stalls, move to next. Commit with terminal verdict.

## Resources
- mini: Python/mpmath, LaTeX
- uisai2 (ramanujan window): Sage + ore_algebra, heavy computation
- ChatGPT Pro: 5 tabs (family1-5), SOL Pro for hard proofs
- Hot Fable: spawned, SendMessage to a05c6933b882f685a
- Ripple: Chudnovsky, Catalan, Apéry, ₃F₂, Clausen, CM infrastructure
- Q-series: Rogers-Ramanujan, quintic cyclotomic, partition congruences
