# CODEX SPEC: per-pair margin anatomy for the [CRIT-2H] two-regime comparison

## Why
The all-h closure now hinges on a PARTITIONED comparison of the certificate margins:
- NEAR pairs (cell positions differing by O(1/h)): proved local theorem, gap = kappa_parity h^-2 + O(h^-3),
  kappa_odd = F''(1/2)/2, kappa_even = F''(1/2) (parity factor 2 proved).
- FAR pairs (macroscopically separated cell positions): believed bounded below by a positive constant
  (or by c'|t_a - t_b|^2), so a crude O(h^-2) tail suffices there.
The existing scan (CRON_kinf_branch.py, CODEX_KINF_report.md, CRON_kinf_results.json) reports only the
AGGREGATE minimum separation per height. We need the anatomy to confirm/refute the two-regime picture and to
extract the constants.

## Task
Extend CRON_kinf_branch.py (new flag, do not break existing behaviour) to emit, for each height h in a
requested range, the per-pair data needed below, and run it for h = 20..60 (plus h = 80, 100 if the cost is
acceptable — the previous full 2..60 pass took ~700 s).

For each h:
1. For every certified mirror-orbit j = 1..2h-2, record: the orbit's cell position t_{h,j} (use the natural
   normalisation t = (j-1/2)/h, or whatever the code's internal orbit ordering gives — state the convention),
   the squared critical value V_{h,j} (midpoint + radius of the Arb ball), and the local branch label if the
   code distinguishes branches.
2. For every non-mirror pair (a,b): the relative separation |V_a - V_b| / max(|V_a|,|V_b|) and the cell
   distance |t_a - t_b|.
3. Emit: (i) the minimising pair for each h (which pair actually attains the reported relative_margin — is it
   always a CENTRAL ADJACENT pair, as the theory predicts?); (ii) the profile of min-separation restricted to
   pairs with cell distance in dyadic bands (|t_a-t_b| in [2^-k, 2^{-k+1})), as a function of h — the
   prediction is: bands with |t_a-t_b| >> 1/h have separation bounded BELOW uniformly in h, while the band
   |t_a-t_b| ~ 1/h carries the h^-2 decay; (iii) fitted constants: for the central adjacent gaps, fit
   gap * h^2 -> kappa_parity separately for odd/even h and report kappa_odd, kappa_even and their ratio
   (predicted exactly 2); for the far bands, report the empirical lower bound and whether it looks like a
   constant or like c'|t_a-t_b|^2.
4. If F''(1/2) can be extracted from the data (kappa_even = F''(1/2) under the proved identity), report the
   numerical value of F''(1/2) and check kappa_odd = kappa_even/2 to the achieved precision.

## Deliverables (this directory)
- CRON_kinf_branch.py updated (flag e.g. --pair-anatomy, output JSON + a readable table).
- CODEX_MARGINS_report.md — first line: TWO-REGIME CONFIRMED / REFUTED (+ the fitted kappa_odd, kappa_even,
  ratio, F''(1/2), and the far-band lower bound). Then the tables and the caveats.
- Raw JSON alongside.

## Rules
- Do not weaken any existing certificate path; the anatomy pass may run at lower precision if that is the only
  way to afford larger h, but say so explicitly per table.
- Every fitted constant reported with the fit range and residuals.
- No effort cap. If the orbit ordering makes the cell-position convention ambiguous, define it explicitly and
  show that the reported anatomy is convention-independent where it matters (the minimising pair and the
  dyadic-band profile).
