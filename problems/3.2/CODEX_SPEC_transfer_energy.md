# CODEX SPEC — TE_{5/4} collision-energy falsification scan (high)

## Mission

Numerically stress-test the new first-lemma candidate of the prize line:

TE_{5/4}: for every prime p >= 7, with N_p(a) = #{1 <= r <= p-2 : b_r ≡ a mod p}
(b_r = Apéry zeta(3) numbers via (n+1)^3 b_{n+1} = (34n^3+51n^2+27n+5) b_n - n^3 b_{n-1}, b_0=1, b_1=5),
the off-diagonal collision energy
E_p^off = sum_a N_p(a)(N_p(a)-1)
satisfies E_p^off <= C p^{5/4}.
Consequence: |Z_p| <= 1 + sqrt(C) p^{5/8} (beats our proved 2/3 exponent).
Poisson prediction: E_p^off = Theta(p); a sibling measurement found ~3p.

## Protocol (one pass per prime; exact arithmetic mod p)

For all primes p in [10^3, 10^6] IF runtime allows (start with [10^3, 10^5],
extend by dyadic blocks as time permits; report exactly the range covered):

1. Compute b_0..b_{p-2} mod p by the recurrence (modular inverse of (n+1)^3).
2. Histogram N_p(a); compute E_p^off, max_a N_p(a), |Z_p| = N_p(0).
3. Record normalized: E_p^off / p, E_p^off / (p log^2 p), E_p^off / p^{5/4}.
4. Dyadic aggregation: for each dyadic block X < p <= 2X report
   M(X) = max E_p^off / p^{5/4} and the argmax prime; also the max of
   E_p^off / p and its argmax. Log-log regression of block maxima of E_p^off
   vs p: report the fitted slope (danger threshold: stable slope > 1.25).
5. Reflection separation: b_r = b_{p-1-r} forces pair collisions; compute
   E_p^refl = the contribution from pairs (r, p-1-r) alone, and report
   E_p^off - E_p^refl (the non-forced energy) with the same normalizations.
6. Gap spectrum (for the largest ~50 primes tested and any spike primes):
   for each value a with N_p(a) >= 2, increment C_p(h) for every gap h
   between distinct positions with value a (cost O(p + E_p^off)).
   Separate forced reflection gaps. Report the top-5 non-reflection (h, C_p(h))
   per inspected prime, and whether any single gap h contributes more than
   p^{1/4} collisions (danger signal).
7. Spike autopsy: for the 5 worst primes by E_p^off/p, print |Z_p|,
   max N_p(a) and its value a, whether the spike is concentrated in one gap
   or diffuse, and (cheap version) whether the positions r with the popular
   value a show any arithmetic-progression or reflection structure.
8. Envelope tests: count how many primes violate E_p^off <= 100 p log^2 p
   (violation does NOT falsify TE_{5/4} but flags resonances — list them all).

## Deliverables

- Script: research/scripts/q32_transfer_energy_scan.py (stdlib only; write
  incremental results to research/scripts/te_scan_progress.txt so partial
  runs are salvageable; print a progress line every 30s of runtime).
- Report: CODEX_TE_SCAN.md (problems/3.2 root): tables per dyadic block,
  regression slope, reflection-separated energies, gap-spectrum findings,
  spike autopsies, and a verdict paragraph: does the data support attacking
  the stronger p^{1+eps} energy theorem, the 5/4 target, or neither
  (i.e. was a genuine resonance found)?

## Hard constraints

- Do NOT touch proof.tex, ERRATA.md, CAMPAIGN_MAP*, lean/, existing scripts
  or reports. New files only (the three named above).
- No external LLM/bridge dispatch.
- Exact integer arithmetic for the recurrence; no floats in E_p computation.

## Acceptance

Script exits 0 on the base range [10^3, 10^5]; report matches script output;
sanity gate printed per 100th prime: |Z_p| from histogram equals direct count
of zeros in the residue list, and the reflection symmetry b_r = b_{p-1-r}
verified exactly for 3 random r per sampled prime.
