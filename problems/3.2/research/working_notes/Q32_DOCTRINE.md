# DOCTRINE — pointwise Apéry gcd (Problem 3.2), automode run 2026-07-21

**Goal (one sentence):** prove fully-unconditional POINTWISE Problem 3.2 (log G_n = o(n) for every n),
or find a genuinely-new avenue that dissolves the six-fold barrier — do NOT re-run the closed families.

## What is already closed (do NOT re-attack these)
Six technology families (Fable R9–R13), each a mechanism-level theorem: analytic (equidistribution),
algebraic-geometric (Galois-image), Diophantine (transcendence/3n-cap), dynamical (integrability/SL₂),
additive-combinatorial (phantom-transparency), single-n adelic (subspace/period, 6 gates). Unified by
phantom-transparency + the critical adelic ledger. Re-running any of these = certified unproductive.
Pointwise ⟺ (SA)/(Ω)/N_coll (numerically Poisson-true, `Q3.2_research_program_B.md`).

## Genuinely-NEW avenues (ranked), each a complete attack plan
- **(a) Construct the "product formula for periods" DIRECTLY.** Fable characterized the missing theorem
  (place-mixing Betti/de-Rham vs crystalline at one integer point). Nobody has tried to BUILD it via
  p-adic Hodge / comparison isomorphisms (B_dR, B_crys) for the Apéry K3 motive. Terminal: a valuation
  inequality forcing Σ v_p(b_n) log p = o(n), OR a proof no such comparison can exist.
- **(c) The FULL modular/motivic structure of 8.4.a.a — beyond a_p equidistribution.** Closed families used
  only the DISTRIBUTION of a_p. Un-audited: the L-function special values / Beilinson–Bloch regulators /
  the actual K3 period lattice / the mirror map. Does the motive's finer arithmetic give multi-place
  valuation control the a_p-distribution cannot? Terminal: a valuation handle, or proof it reduces to (SA).
- **(d) Fresh DEEP numerical structure search (do FIRST — cheap, may dissolve the wall).** 13 rounds tested
  low-degree/single-prime structure. Un-tested at depth: (i) cross-prime coincidences (z_p vs z_{p'} for
  p'=2p±1, p'≡p mod small d, twin-ish); (ii) degree-≥3 fixed congruences for z_p; (iii) the u_p projective
  datum cross-prime; (iv) Dwork-tower z_p vs z_{p^k}; (v) correlation of the WHOLE Z_p with a_p(8.4.a.a)
  at r≠(p-1)/2 via the ₄F₃ evaluation. Terminal: a reproducible cross-prime/higher-degree law (⟹ reopens
  a compressive fibration ⟹ Lemma-M-style proof), OR clean null across all (⟹ confirms structurelessness).
- **(b) A DIFFERENT sufficient condition (different decomposition of G_n).** Every route went support-law →
  Z_p → collision. Un-tried: decompose G_n via the reduced denominator B_n=d_n b_n/G_n directly, or via a
  resultant/discriminant of the (a_n,b_n) approximation, seeking a tractable target that AVOIDS (SA).
  Terminal: a new sufficient condition provable by an open family, or proof it re-enters (SA).

## Fallbacks
If (a)–(d) all reach proof-of-failure: the six-fold barrier stands; deliverable = the complete status
resolution (density unconditional + reduction + six barriers + two principles + research-program-B).

## Dispatch coupling (per user: "等 codex 终稿就派 Fable 和 Chatgpt 继续")
- Codex v4 finalizing the paper (parallel, not blocking).
- On Codex v4 done: dispatch Fable on (a) + (c); ChatGPT (dm channels) on the p-adic-Hodge / K3-motive /
  Beilinson literature for (a)/(c). No effort cap in briefs.
- Agent (me) drives (d) numerically NOW, then orchestrates.

## RUN VERDICTS
- **(d) CLOSED — proof-of-failure (clean null).** Deep search (primes<6000): no cubic congruence (>10% hits
  NONE); center-fold histogram flat; consecutive-prime normalized-zero closeness REAL 2.9% vs NULL(shuffle)
  mean 5.1% — real ≤ null ⟹ NO cross-prime structure; only known reflection pairs (z,p-1-z) appear; zero
  exact cross-prime coincidences. Confirms structurelessness; the wall does NOT dissolve numerically. 两问:
  (least-sure) subtle >pairwise correlations untested by these statistics; (blind-angle) an algebraic
  structure invisible to statistics but visible to (a)/(c) is exactly what Fable audits next.
- **(a)+(c) CLOSED — Fable R14 (§5.42), branch (b) with a NEW stone.** Construction of the "product formula
  for periods" from the actual motivic input degenerates via five gates: Invariant-Ring Gate (comparison
  scalars factor through difference-Galois invariants = Casoratian line; height spectrum gapped O(log n) vs
  ≥λ′n — no sublinear window-supported eliminant); No-Interpolation Gate (digit motive's parameters
  (n mod p)/(p−1) move with p; no cross-p interpolation object exists; fixed-motive locus = torsion-digit
  channels p|kn+c of LOG mass; k=2 = Ahlgren–Ono center digit = open weight-4 ordinarity); EDS Calibration
  (THE new stone: for EDS the product formula EXISTS exactly — v_p(B_m)=λ_p(mP), Σλ_v=ĥ — plus a perfect
  fixed-quotient apparition law, and the window cap is STILL GRH/Vojta-conditional ⟹ the constructed object,
  even granted, does NOT locate — R13's characterization was necessary-format, not sufficient); Boundary Gate
  (FL/BMS valuation content = universal Hodge-tax (n-blind) + Beukers–Vlasenko vertical tower (digit-blind);
  unit-root non-overconvergence (Dwork–Coleman) is intrinsic; K3 torsion channels empty); Direction Gate
  (all known place-mixing theorems run many-p-hypothesis ⟹ rigid conclusion; the reverse arrow has only
  height ledgers and subspace, both closed; Pila–Zannier preempted by the deg≥p/C stone). Third unifying
  principle named: QUOTIENT-FREE LOCATING — no known mathematics locates p-local events lacking a fixed-
  quotient avatar. Frontier exact: an effective apparition principle for entropy-positive p-local loci
  without a quotient. Seven-fold map complete; fallback deliverable is now the operative state.
- **(b) CLOSED — proof-of-failure (provably re-enters (SA)).** Any decomposition of the pointwise G_n bound is
  per-prime-local: v_p(G_n)>0 for p∈(√n,n] ⟺ n mod p∈Z_p (support law) — forced, not a choice. So every
  decomposition must control the same Z_p-membership = (SA) (the R13 "sub-3n ⟹ p-local input = Z_p" argument).
  Not a new sufficient condition. 两问: (least-sure) an averaged/dual B_n-side bound might have slack — but the
  density version already exploits that and is the a.a.-n result, not pointwise; (blind-angle) covered by (a)/(c).
- **(a)+(c) frontier bullet — Fable r14 DELIVERED (verdict above; §5.42); ChatGPT Q244 still out.** Fable r14 constructed/closed the place-mixing period inequality
  (p-adic Hodge B_crys/B_dR + Beilinson–Bloch L-values + integral BMS/prismatic at Teichmüller boundary points);
  ChatGPT Q244 on the p-adic-Hodge/K3/Beilinson literature. These use the actual comparison isomorphisms /
  L-values / integral lattice — inputs the six closed families never touched (they used only a_p-distribution /
  height / Galois-quotient). NON-ensemble, NON-height, genuinely un-audited. No effort cap.
