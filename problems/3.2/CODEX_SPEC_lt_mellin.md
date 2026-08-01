# CODEX SPEC: [GAP-LT-MELLIN] — the single surviving frontier of P3.2

## Inherited state (all machine-verified / double-source certified; do NOT re-derive)
Read first: problems/3.2/research/working_notes/FABLE_NOTES_energy_bootstrap.md
sections 43-52; CODEX_FRANEL_MELLIN.md; GPT_Q6394_descent_theorem.md;
GPT_Q6380 (in cron's chatgpt-answers archive or /tmp/gpt/life/Q6380.md);
GPT_Q6371_padic_stickelberger_audit.md.

Facts you may use freely:
- b_r ≡ M(r) − T(r) mod p, exact for 1≤r≤p−2 (scripts franel_*_test.py).
- M(r) = Σ_{x∈F_p} H_p(x)² ω^{−r}(φ(x)), φ(x)=x(1−8x)/(1+x); H_p = Hasse
  invariant of the Beauville-IV elliptic pencil (Γ₁(6)); H_p(x)² ≡ a_{p,x}²
  mod p with a_{p,x} the true elliptic Frobenius trace ⇒ M(r) is the mod-𝔭
  reduction of the genuine Weil-number sum M̃(r) = Σ_x a_{p,x}² χ_r(φ(x)),
  |a_{p,x}| ≤ 2√p, χ_r = ω^{−r} lifted (Teichmüller/char-0 character).
- T(r) = Mellin trace of a fixed integral rank-3 system 𝒜₋ on G_m, absolute
  conductor ≤ 11 (double-source: Q6394 + codex-fm commit 8e0da63); so T(r)
  also lifts to a Weil-number sum T̃(r), pure weight ≤ 3, bounded conductor.
- Parity law (proved): |Z_p| ≡ 1[p | b_{(p−1)/2}] mod 2.
- Empirical model (cron appendix Q, 3242 primes p<30000): |Z_p| =
  2·Poisson(0.509) + midpoint indicator; mean 1.0185; max 8; α_empirical = 0.
- Analytic record: |Z_p| ≤ N_p(0) ≤ 3p^{2/3} (all-value 8p^{3/4}).
- Dead ends (proved dead, do not revisit): HGM-degree shortcut (Q6393),
  low-degree U_p fits, naive Ax-Katz/Adolphson-Sperber/Katz-Mellin transfer
  (Q6371), bounded-length GK monomials [GAP-BGK].

## The target
Bound #{r ∈ [1,p−2] : 𝔭 | M̃(r) − T̃(r)} = o(p), any p^{1−δ} is gold.

## Tasks (commit prefix "codex-lt:", report to problems/3.2/CODEX_LT_MELLIN.md)
1. (both workers) GALOIS ORBIT PRODUCT: σ_a ∈ Gal(Q(ζ_{p−1})/Q) acts by
   σ_a(M̃(r)) = M̃(ar mod p−1) (verify this action claim carefully — the
   Teichmüller character behaves as χ_r ∘ σ_a = χ_{ar}; check also the effect
   on a_{p,x} — they are rational integers, untouched). Hence for the orbit
   O_r = {ar}: P_r := Π_{s∈O_r} (M̃(s) − T̃(s)) ∈ Z (rational integer!), with
   |P_r| ≤ (C p^{3/2})^{|O_r|}. If k members of O_r are divisible by 𝔭 AND
   the product is nonzero, then p^k | P_r... careful: 𝔭 | each factor gives
   N(𝔭)=p dividing... work out the exact divisibility transfer (𝔭 vs its
   conjugates; each conjugate prime divides the corresponding conjugated
   factor — the SAME 𝔭 must divide P_r to order ≥ k? No: σ(𝔭)-divisibilities
   multiply to give p^k | P_r exactly when the k hits are in ONE orbit).
   Derive the strongest unconditional consequence: e.g. "if P_r ≠ 0 then
   #hits in O_r ≤ (3/2)|O_r|·log_p(Cp) = (3/2+o(1))|O_r|" — useless if ≥|O_r|;
   compute the constants PRECISELY: |P_r| ≤ (cp^{3/2})^{|O_r|} and p^k | P_r
   force k ≤ |O_r|·(3/2 + log_p c). The margin is the whole question: if the
   true size of typical |M̃| is p^{3/2} then k ≤ 1.5|O_r| — vacuous. BUT if we
   use the SECOND MOMENT (Parseval: Σ_r |M̃(r)|² = (p−1)·Σ_x |a_x²|² ≈ 2p⁴ ⇒
   typical |M̃(r)| ≈ p^{3/2}... compute exactly) combined with AM-GM on the
   orbit product |P_r| ≤ (orbit mean square)^{|O_r|/2}, the bound sharpens
   for orbits whose members are typical. Push this to its limit; identify
   for which orbit sizes/structures it becomes nonvacuous. Also handle the
   P_r = 0 case (M̃(s) = T̃(s) exactly in char 0 — when? purity weights differ
   (2+? vs ≤3): if weights genuinely differ, exact equality has bounded rank
   of solutions — nail this).
2. (max) THE VERTICAL FAMILY / WEIL II ROUTE: the family {M̃(r)}_r is the
   discrete Mellin transform of a fixed trace function. Katz's theory gives
   perversity/purity of the Mellin object; the r-th value is the trace of
   Frob on a two-dimensional (per Q6394 §6) cohomology. Zero mod 𝔭 of the
   trace = non-ordinarity of that r-fiber. So #{r: hit} = #non-ordinary
   fibers in a FAMILY OVER THE CHARACTER GROUP. For families over a CURVE,
   non-ordinarity is a divisor (Hasse locus) of degree O(genus/conductor
   data) — the character group is not a curve (Q6393 killed the naive
   version), BUT: the Mellin family over 𝔾_m^∨ in Katz's Tannakian sense
   might carry an F-crystal structure over a formal/rigid parameter (Wan's
   work on T-adic/eigencurves families of exponential sums, Liu-Wan, Davis-
   Wan-Xiao). Investigate precisely: is there a T-adic L-function whose
   Newton polygon controls ordinarity across ALL twists ω^{-r} at once
   (r interpolated p-adically)? If yes, the number of non-ordinary twists is
   controlled by Newton-vs-Hodge polygon gap data — get the exact statement
   and check whether it yields o(p) or even O(polylog). This is the deepest
   shot: T-adic exponential sums ARE this exact setup (twist aspect!).
3. (high) MECHANICAL CHORES: (a) print the explicit isogeny formula missing
   from Q6394 §3 ([GAP-1: explicit isogeny formula]) — Weierstrass model of
   Beauville-IV, the modular correspondence, rational maps; verify
   numerically at 5 primes. (b) adjudicate the 3F2(1/3,2/3,1;1,1) literal
   identification [SUSPECT] against CFVZ arXiv:2510.23298 conventions:
   rank-2 Franel elliptic vs rank-3 Apery Sym² — which object does CFVZ's
   finite hypergeometric formula actually name? (c) numerically verify the
   Galois action claim σ_a(M̃(r)) = M̃(ar) at p=13,17 by computing M̃ in
   Q(ζ_{p−1}) exactly (sympy cyclotomics or pari via cypari if available).
4. (both) Every asserted identity gets a verification script
   research/scripts/codex_lt_*.py printing VERIFIED lines; run before commit.

Rules: no fabricated citations (mark UNSURE); [GAP-n] every unproved step;
end report with least-confident step; rebase never force-push.
