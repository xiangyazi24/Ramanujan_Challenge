# CODEX SPEC — bad-diagonal inverse theorem: formalize + stress-test (high)

## Mission (single deliverable)

Direction 3 of the surviving program ("坏对角逆定理", Q6339 framing): the
two-layer Mellin inverse theorem. Statement shape: if two rank-<=3 Mellin
objects M_1, M_2 over G_m/F_p satisfy M_1(chi) = M_2(chi) (as Frobenius
traces / exact values) for MANY characters chi (more than a bounded number),
then the pair is explained by a self-twist / twist relation between the
layers; the known exception list has 10 items, and the Apery quadratic
companion is the first genuinely "graph"-type case (the twisting object
K_q = chi_2(q(t)) is NOT a Kummer character of t, so it cannot be absorbed
into a character twist).

Your job: (1) reconstruct the precise statement + the 10-item exception
table from our archives, (2) write it down as a formal conjecture with all
quantifiers pinned, (3) stress-test it numerically on our own objects, and
(4) hunt for small counterexamples. Recon quality bar: a hostile referee
should not find a quantifier you left floating.

## Read first

1. `CRON_FRESH_EYES_pointwise.md` appendix Z (lines ~511-531), especially
   Z.4 (Package C and the "逆定理例外表 10 项" sentence) and Z.3.
2. `chatgpt-answers/Q6413.md` — the q18 answer summary (exception table
   source). If the 10 items are not fully listed there, reconstruct what IS
   listed and mark the gaps explicitly — do not invent items.
3. `grep -n "Q6339" research/working_notes/FABLE_NOTES_energy_bootstrap.md`
   and read the +-30 lines around each hit — the original framing.
4. `CODEX_LT_MELLIN.md` sections 1-2 — the orbit/divorce analysis that any
   inverse theorem must be consistent with.
5. `CRON_pushforward_check_report.md` + `CRON_pushforward_check.py` — exact
   conventions and reusable arithmetic for our concrete objects.

## Tasks

1. FORMAL STATEMENT: write the conjecture (both an "exact-value" version and
   a "trace-function" version if the sources distinguish them) with explicit
   quantifiers: rank bounds, conductor bounds, the threshold "many" (is it
   > C(rank, cond)? > p^epsilon? — extract what the sources actually claim,
   and where they are silent, state BOTH candidate thresholds and label them
   [source-claimed] vs [our-guess]).
2. EXCEPTION TABLE: reproduce the 10 items as a numbered markdown table with
   for each: the pair (M_1, M_2), the twisting mechanism, why it evades the
   naive statement, and the source line you got it from. Mark unrecoverable
   items as [NOT RECOVERED FROM ARCHIVES].
3. NUMERICAL STRESS TEST (stdlib Python, exact arithmetic):
   a. Build a small zoo of rank-<=3 Mellin trace functions over F_p for
      p in {29,...,101} (at minimum): the Apery object, its quadratic
      companion (chi_2(q(t)) twist), Kummer twists of order <= 6 of each,
      a Legendre-family Sym^2, and 2-3 "unrelated" hypergeometric sheaves.
   b. For every unordered pair in the zoo, count
      #{chi : S_1(chi) = S_2(chi)} (exact equality of algebraic-integer
      traces; ALSO count equality of centered mod-p reductions separately —
      keep the two counts in separate columns, this distinction is the whole
      point of our project).
   c. Verify: related pairs (twist relations) give ~p matches or a positive
      proportion; unrelated pairs give O(1) matches. Print the full matrix.
   d. Specifically test the Apery/companion pair — the claimed FIRST graph
      case: does its match count behave like "related" (positive proportion)
      or "unrelated" (bounded)? This single number decides whether the graph
      case genuinely breaks the naive inverse theorem.
4. COUNTEREXAMPLE HUNT: search the zoo (+ random character sums of rank <= 2
   if time permits) for any UNRELATED pair whose match count grows with p —
   i.e. a violation of the conjectured dichotomy that is NOT in the
   exception table. Report negative results with the exact search space.

## Deliverables

- Script: `research/scripts/q32_bad_diagonal_stress.py` (stdlib only,
  prints the match-count matrix per p, then a SUMMARY block).
- Report: `CODEX_BAD_DIAGONAL.md` (problems/3.2 root):
  formal conjecture (quantifiers pinned, [source-claimed]/[our-guess]
  labeled), the exception table, the numerical matrices, the Apery/companion
  verdict, counterexample-hunt result, and a final section "what a proof
  would need" (name the expected tools: multiplicative-translation
  independence, Deligne bilinear ceiling p^{3/4} — cite Z.4's calibration —
  and where the known techniques stop).

## Hard constraints

- Do NOT touch: proof.tex, FABLE_SECTION_*.tex, ERRATA.md,
  CAMPAIGN_MAP_2026-08-01.md, anything under lean/, any existing script.
  New files only.
- Do NOT dispatch to any bridge/ChatGPT/external LLM. Self-contained.
- Never fabricate an archive item: everything in the exception table needs a
  file+line citation or the [NOT RECOVERED] tag.

## Verification & acceptance

- `python3 research/scripts/q32_bad_diagonal_stress.py` exits 0.
- Report exists; every archive citation resolves (file+line).

## Stall protocol

If the archives genuinely do not contain the 10-item table: say so in the
report, deliver the formal statement + stress test anyway, and list exactly
which archive files you searched (with grep patterns used).
