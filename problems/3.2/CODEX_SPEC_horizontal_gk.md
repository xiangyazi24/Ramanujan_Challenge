# SPEC: the cross-prime counter-attack — rank-2 Gross-Koblitz horizontalization

Repo ~/repos/Ramanujan_Challenge, workdir problems/3.2/. python3 free. Read first:
research/working_notes/FABLE_SECTION_apparition_tower.tex (the proved vertical theory),
CODEX_JACOBSTHAL_DEEP.md (2F1 closed forms for tau, sigma), FABLE_NOTES_energy_bootstrap.md
sections 34/39 (character-order counting lemma; canonical frontier = high-order Mellin
diagonal zero-density).

NEW WEAPON (why this counter-attack is now possible): tau = sqrt(F), sigma = sqrt(F/q) are
algebraic pullbacks of ONE classical 2F1 with exact Lagrange-inversion coefficient formulas.
At RANK 2 the Gross-Koblitz program is LEGITIMATE (unlike the failed rank-3 attempt): the
Mellin coefficients of a hypergeometric-class rank-2 local system are BOUNDED Gauss/Jacobi-sum
monomials. Consequently:
  b_r mod p = (convolution of two branch coefficients) = an EXPLICIT quadratic Gauss-sum object.

MISSION (deep): horizontalize.
1. Derive the explicit finite-field formula: tau_j mod p (resp sigma_j) as a bounded product
   of p-adic Gamma / Jacobi-sum values at arguments linear in j (use the 2F1 parameters from
   the closed form + Gross-Koblitz; handle the algebraic pullback phi carefully — the pullback
   may turn single Jacobi sums into short sums; get the exact object). Machine-verify at
   p = 13, 29 for all j.
2. Substitute into b_r = sum_{i+j=r} tau_i tau_j (chi=+1; sigma-version with q-correction
   otherwise): b_r mod p = an explicit double character sum ("quadratic Jacobi convolution").
   Machine-verify against b_r directly.
3. THE HORIZONTAL OBJECT: for fixed n, the bad-column events are p | b_{n mod p}-shifted with
   the Mellin character omega_p^{n-1} (counting lemma frame). Rewrite the fourth moment
   F_4(N) (or even just the pair correlation C_{p,q}) using the formula from step 2: the
   cross-prime correlation becomes a complete exponential sum in several variables over
   F_p x F_q. State exactly what Weil/Deligne CAN now bound that was unreachable without the
   closed form (the p-side sum alone at fixed q-data? any genuine two-prime saving?). Be
   brutally honest about where the two-characteristic obstruction reappears; but hunt for
   ANY partial horizontal estimate (e.g. average over n of the pair correlation with an
   explicit power saving) that the explicit formula newly enables.
4. Report to problems/3.2/CODEX_HORIZONTAL_GK.md, commit prefix "codex-hgk:". Machine-verify
   every formula; mark gaps [GAP].
