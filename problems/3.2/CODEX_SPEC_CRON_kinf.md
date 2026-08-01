# CODEX SPEC: K_infinity branch table + h0 scan ([CRIT-2H] closure chain)

## Mission
Implement and run the certified/uncertified branch-table algorithm delivered in CHATGPT_Q6723_kinf_algorithm.md (the authoritative spec; read it fully), which instantiates the uniform conjugation theorem of CHATGPT_Q6708_conjugation_theorem.md for the Apery cell recurrence. Goal chain: [ARCH-CELL] squeeze => [CRIT-2H] (s_h >= 2h-1) for all h > h0, combined with banked machine certificates for h <= 40 => if h0 <= 40, the G_h absolute-irreducibility chain closes for ALL h (campaign prize).

## Deliverables (this directory)
1. CRON_kinf_branch.py — implementation of: the K_infinity(z) cell object (per Q6723: K_inf(z) = phi(-z)phi(z) + z^6 gamma(-z)gamma(z), with the phi/gamma series defined there; d_+ = d_- = -3/2 drift cancellation already established, no n^d factor), zero isolation of H_inf(z) = z K_inf'(z) - 3 K_inf(z), and the fixed-h direct evaluation mode (Q6723 section 7 pseudocode).
2. First pass UNCERTIFIED: mpmath at >= 100 digits, scan h = 2..60: for each h emit YES/NO/UNRESOLVED for the cell-squeeze criterion (at least h mirror-orbits of squared critical values simple, nonzero, globally distinct — exact criterion in Q6723). Emit the empirical h0 = largest h that fails (or "none <= 60").
3. If arb/flint ball arithmetic is available (python-flint IS installed — check import flint), redo the scan certified per the interval Newton/Krawczyk test in Q6723; report which h get certified YES.
4. CODEX_KINF_report.md — verdict first line: H0_EMPIRICAL <= 40: YES/NO + certified range; then tables (per h: #zeros isolated, min separation, margin), then caveats. Machine outputs inline.

## Rules
- The Q6723 document is the spec; where it says "increase precision or subdivide", do so before reporting NO.
- No effort cap. If the phi/gamma series definitions in Q6723 are ambiguous at any point, reconstruct them from the recurrence first principles (the cell system is the Apery recurrence linearization at the archimedean cell; A_inf = [[34,-1],[1,0]], A_1 = [[-51,3],[0,0]]) and state your reconstruction explicitly in the report.
- Sanity anchors: lambda_+- = 17 +- 12 sqrt 2; the h <= 30 range must reproduce consistency with the banked crit2h Morse certificates (no contradiction: cells where our certificates say full Morse must not be refuted by your run).
