# CODEX SPEC: Franel–Mellin bounded object — the new starting point after [GAP-BGK] bypass

## Status you inherit (all machine-verified, do not re-derive)

- Pointwise identity: `A_p(phi(x)) = H_p(x)^2` for all `x in F_p \ {-1}` with
  `phi(x)=x(1-8x)/(1+x)`, `H_p = truncated Franel GF` (scripts
  `research/scripts/franel_pointwise_test.py`, cron's commit, both green
  p=13..101, two chi-classes).
- Fiber discriminant of `phi(x)=t` is exactly `q(t)=t^2-34t+1` (=1-34t+t^2
  reversed/same symmetric), fiber count `nu(t)=1+chi2(q(t))`.
- EXACT decomposition, verified (`franel_mellin_mult_test.py`, p=13,29,37,
  all r tested, zero correction terms for 1<=r<=p-2):

      b_r  =  - SUM_{x in F_p, x != -1, phi(x)!=0} H_p(x)^2 phi(x)^{-r}
              + T(r),
      T(r) =  SUM_{t in F_p^*} chi2(q(t)) A_p(t) t^{-r}   (mod p).

  So the divisibility condition p | b_r is the vanishing of a Mellin value of
  the FIXED pair (Franel-square pushforward, quadratic-twist companion); r
  enters only through the character t^{-r}. This bypasses [GAP-BGK] at the
  object level (no Lagrange extraction, no term explosion).

## Tasks (in order; commit prefix "codex-fm:")

1. **Sheaf-theoretic normalization.** Express both summands as trace functions
   of explicit bounded-conductor objects on G_m: (a) the pushforward
   phi_*(FranelSquare) — identify FranelSquare: H_p(x) mod p is the truncation
   of the Franel GF; relate pointwise H_p(x) to the Hasse–Witt/unit-root of the
   Franel elliptic/K3 family fiber at x (Dwork congruence h = H_p * h^Frob).
   (b) the companion chi2(q(t)) A_p(t): use A_p = S_p^2 (chi_p=+1) or
   q S_p^2 (chi_p=-1) and CFVZ (arXiv 2510.23298) conventions to write it as a
   trace-function of the quadratic twist of the same underlying rank-2 object.
   Deliverable: precise statement "T(r) = Mellin transform at omega^{-r} of the
   trace function of sheaf G, cond(G) <= C absolute".
2. **Monodromy & self-twists of the fixed pair.** Compute local monodromy
   (singular points, exponents) for both sheaves; classify their self-twists
   and mutual twists by Kummer characters (this feeds the counting-lemma
   reduction: only bounded-order character twists may correlate). Prove or
   disprove: no Kummer twist of unbounded order fixes either sheaf.
3. **What does Katz-style Mellin equidistribution give?** For the r-family
   {Mellin value at omega^{-r}}_{r<p}, apply the Katz "Mellin transform over
   finite fields" / Deligne equidistribution machinery: state exactly what
   zero-density for the vanishing locus {r : value = 0 mod p} would follow
   from (i) the computed monodromy, (ii) any standard hypothesis. NOTE the
   p-adic caveat: the value is in F_p, vanishing is a divisibility not a
   complex cancellation — identify whether Katz's F_p-valued Mellin theory
   (Katz, "Convolution and equidistribution", Gauss-sum settings) applies to
   F_p-valued trace functions directly, or only after Teichmuller lift.
4. **Verification scripts** for every identity you assert, in
   research/scripts/ (codex_fm_*.py), each printing VERIFIED lines; run them;
   report exact commands + outputs in problems/3.2/CODEX_FRANEL_MELLIN.md
   (ledger format: [VERIFIED]/[GAP-n]/[NEGATIVE-...]).

Rules: never fabricate a citation; every unproved step gets [GAP-n]; end the
report with your least-confident step. Repo may receive concurrent commits —
rebase, never force-push.
