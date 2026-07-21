# ORACLE C: explicit marked Mellin coordinate c_{p,j} (the mh2 algebraic prerequisite)

## Where this sits

The dual-oracle meeting (oracleA_result.tex + oracleB_result.tex, both in
proof.tex now) reduced fully-unconditional to eq:oracleB-mh2, the
two-characteristic crystalline Mellin dispersion. That package (eq:oracleB-
crystalline-mellin) needs, BEFORE any dispersion analysis, a MARKED COORDINATE:
an explicit algebraic quantity c_{p,j} with
   c_{p,j} = 0 (in F_p, or in a fixed-rank space)  ⟺  b_j ≡ 0 (mod p),
whose complexity (degree / conductor / #components) is controlled uniformly in
BOTH p and the moving index j. This spec attacks ONLY that algebraic
prerequisite — the dispersion proof is a separate (analytic) job.

## The exact obstruction to beat (from Oracle B)

- Z_p = {j : [t^j]H_p(t) = 0} is a COEFFICIENT / finite-Mellin zero set, via
  b_j = -Σ_{t∈F_p^*} H_p(t) t^{-j} (prop:oracleB-two-zero-loci, VERIFIED p=7).
  It is NOT the evaluation-root locus of H_p or of B_p (H_p = Δ^{ε} B_p²).
- The naive marked coordinate "j ↦ (is j a root of B_p)" FAILS twice: (i) wrong
  object (eval root ≠ coefficient zero); (ii) deg B_p = Θ(p), genus Θ(p), so
  no fixed bounded-degree family (prop:oracleB-complexity).
- So a working c_{p,j} must read off a COEFFICIENT vanishing, and must do so
  with bounded complexity uniformly in j — exactly what neither B_p nor the
  gap polynomials N_h (also degree Θ(h), h up to Θ(p)) provide.

## Ranked goals

(G1) CONSTRUCT an explicit c_{p,j} with c_{p,j}=0 ⟺ p|b_j and complexity O(1)
     uniform in (p,j). Candidate shapes to try (verify each symbolically):
     - A fixed-size (r×r, r=O(1)) determinant/Wronskian in the transfer matrices
       M(n)=[[P(n),-n^6],[1,0]]: b_j is a fixed matrix coefficient of
       M(j-1)...M(0), so p|b_j is the vanishing of ONE entry of a product — is
       there a bounded-size "marked" determinant reading it with controlled
       conductor as j varies? (The product length grows with j; the question is
       whether a RESCALED / gauge-transformed cocycle has bounded complexity.)
     - The Mellin dual: c_{p,j} = Σ_{t} H_p(t) t^{-j} is exactly -b_j — a single
       additive-character-weighted sum; its "complexity" as a function of j is
       the key. Is j ↦ b_j mod p the trace of a RANK-2 object (the 2-dim Apéry
       local system) twisted by the Kummer sheaf K_{ω^{-j}}? If the twist keeps
       rank 2 and bounded conductor uniformly in j, THAT is the marked coordinate
       and it opens Oracle A's route. Check rank/conductor concretely.
     - Frobenius trace form: b_j mod p as tr(Frob_p | some fixed motive ⊗ ω^{-j}).
(G2) Prove ANY marked coordinate reading coefficient-vanishing must have
     growing complexity (sharpen prop:oracleB-complexity to ALL constructions,
     not just B_p) — a genuine no-go that would redirect the whole program.
(G3) A conditional: IF the rank-2 Kummer-twisted local system has bounded
     conductor uniformly in j (state precisely), THEN c_{p,j} exists and
     Oracle A's SDC/short-arc reduces to ordinary FKM bilinear bounds. Prove
     the reduction; leave the conductor bound as the named open input.
(G4) Sharp obstruction with the exact missing geometric fact named.

## Computational duty (FIRST)

For p ≤ 2000: (a) build the transfer-matrix cocycle M(j-1)...M(0) mod p and
confirm b_j is the fixed (1,2) or (1,1) entry (state which); (b) test whether
j ↦ b_j mod p is a rank-2 linear recurrence in j with COEFFICIENTS that are
themselves trace functions of bounded complexity (it is NOT constant-coeff —
quantify how the complexity grows); (c) for the Kummer twist, compute the
conductor of the sheaf whose trace is j↦b_j (via the number of j-singularities
/ the degree of the L-function) and test if it is bounded in p. Write to
problems/3.2/oracleC_exploration.md, code oracleC_explore.py.

## Deliverables

problems/3.2/oracleC_result.tex (construction / reduction / no-go, notation of
proof.tex + oracleB_result.tex — cite eq:oracleB-mh2, eq:oracleB-crystalline-
mellin, prop:oracleB-two-zero-loci, prop:oracleB-complexity, rem:orbit),
oracleC_explore.py + oracleC_exploration.md, oracleC_verify.py (PASS/FAIL,
nonzero exit). STALL REPORT convention. Do NOT modify existing files. No
numerics-as-proof; verify every identity symbolically; be brutally honest about
the moving-j vs fixed-family quantifier (it killed every prior route). If the
answer is a no-go, that is a first-class result — it tells us fully-unconditional
needs genuinely new geometry, not a repackaging.
