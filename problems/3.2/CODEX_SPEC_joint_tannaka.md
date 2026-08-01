# CODEX SPEC — joint Tannakian group: numerical Goursat fingerprint (max)

## Mission (single deliverable)

Empirically pin the JOINT Mellin-Tannakian monodromy group of the eigen-pair
(A_+, A_-) by computing joint moment statistics of their normalized Mellin
traces across primes. This is the numerical advance work for [GAP-2] / q18
Package B (the surviving program's most mature direction). You are NOT asked
to prove anything — you are asked to produce a clean, exact-arithmetic
fingerprint table plus a cautious verdict.

## Read first (in this order)

1. `CRON_pushforward_check_report.md` — exact conventions (model, cover,
   fiber counts, chi_2, corrected pullback A_p(t) including non-split fibers).
2. `CRON_pushforward_check.py` — working code for T_F(u), T_G(t), N(t),
   chi_2(q(t)), the corrected pullback, and Mellin inversion. REUSE its
   arithmetic; do not re-derive.
3. `CRON_FRESH_EYES_pointwise.md` sections Z.1, AA.2, AA.3 (lines ~511-531) —
   the eigen-decomposition A_± = (pushforward ± chi_2(q(t))·Apery-companion)/2
   correction story, rank accounting (Sym^2 pushforward rank 6, tensor rank 8),
   conductors (A_+: 9, A_-: 11), monodromy A_+ = O_3 (det = K_q), A_- = SO_3.

## Setup (established facts, use as given)

- Cover phi(u) = u(1-8u)/(1+u), branch polynomial q(t) = t^2 - 34t + 1.
- Family E_u: y^2 + (1-2u)xy + u^2 y = x^3 over F_p, T_F(u) = p+1-#E_u(F_p).
- Sym^2 trace per source fiber: T_F(u)^2 - p (weight-2 part; subtract the
  Tate/determinant part correctly — check against the report's convention and
  say explicitly which convention you use).
- Deck eigen-pieces (rank 3 each):
    T_+(t) and T_-(t) with T_+ + T_- = (Sym^2 pushforward trace),
    T_+ - T_- = chi_2(q(t)) * (corrected Apery companion pullback A_p(t)).
  Reconstruct the exact formulas from Z.1/AA.3 + the script; verify your
  reconstruction by checking Mellin inversion reproduces the Apery numbers
  b_r mod p for all r (this MUST pass before any moment is computed —
  it is your gate).
- Mellin trace: S_±(chi) = sum_{t in F_p^*} chi(t) T_±(t) for each
  multiplicative character chi of F_p^*. Normalize s_±(chi) = S_±(chi)/p
  (weight-2 object of rank 3: |s_±| <= 3 up to bounded ramified corrections;
  document the exceptional characters you exclude, e.g. chi where the sheaf
  is ramified/trivial subquotients appear — follow Katz's convention of
  discarding the O(1) bad characters and SAY which ones).

## Compute (exact integer arithmetic where possible; floats only for the
final normalized statistics)

For each prime p in P = {29, 37, 41, 53, 61, 73, 89, 101} plus as many
further primes p ≡ anything, p up to 500 (or as far as ~30 min of compute
allows; report the actual range):

1. Single-object moments over chi != trivial (and excluding documented bad
   chi): M_1(±), M_2(±), M_3(±), M_4(±) where
   M_k(±) = (1/#chi) * sum_chi s_±(chi)^k  (real parts; traces are real
   here — verify and assert this).
   Expected if monodromy is the full O_3 / SO_3 std rep: M_2 -> 1.
   Report the observed M_4 limit and compare with the theoretical 4th moment
   of the std rep of O_3 vs SO_3 (compute these two theoretical values
   yourself by invariant theory / explicit integration and PRINT them next
   to the data — dim of invariants of V^{⊗4}).
2. Joint moments: M_{1,1} = avg s_+ s_-, M_{2,1}, M_{1,2}, M_{2,2},
   M_{3,1}, M_{1,3}.
   Product group O_3 x SO_3 predicts M_{a,b} = M_a(+) * M_b(-).
   A graph/correlated subgroup (Goursat) predicts specific deviations —
   compute the predicted M_{1,1} and M_{2,2} for: (i) product; (ii) the
   diagonal graph SO_3 embedded via an isomorphism between the SO_3 parts
   (with the extra Z/2 of O_3 acting by det); derive these predictions
   explicitly and print a comparison table.
3. Twisted graph detection: for every Kummer character eta of order <= 12,
   compute C_eta = avg_chi s_+(chi) s_-(chi * eta). A graph-up-to-twist
   would show |C_eta| bounded away from 0 for exactly one eta. (Our
   self-twist audits at order <= 30 found only trivial — this is the
   CROSS-twist analogue, which has NOT been checked. This is the most
   valuable single number in this spec.)
4. Convergence: for each statistic, print the trend in p (table rows = p),
   plus a final row of the average over the largest 5 primes. State the
   fluctuation scale you expect (~ p^{-1/2} per character-average by
   square-root cancellation) and whether observed deviations are within it.

## Deliverables

- Script: `research/scripts/q32_joint_tannaka_moments.py` (stdlib only,
  runnable as `python3 research/scripts/q32_joint_tannaka_moments.py`,
  prints VERIFIED/REFUTED gates then the tables).
- Report: `CODEX_JOINT_TANNAKA.md` (repo root of problems/3.2) with:
  (a) the exact eigen-trace formulas you used (so life can re-derive),
  (b) the gate results (Mellin inversion reproduces b_r for ALL p tested),
  (c) moment tables + theoretical predictions side by side,
  (d) one-paragraph cautious verdict: product vs graph vs twisted-graph vs
      inconclusive, with the single strongest numerical discriminator named,
  (e) an honest LIMITATIONS section (what a numerical fingerprint can never
      settle: arithmetic vs geometric group, small-index subgroups, etc.).

## Hard constraints

- Do NOT touch: proof.tex, FABLE_SECTION_*.tex, ERRATA.md,
  CAMPAIGN_MAP_2026-08-01.md, anything under lean/, any existing script.
  New files only (the two paths above).
- Do NOT dispatch to any bridge/ChatGPT/external LLM. Self-contained.
- Exact arithmetic for all finite-field computations; assert every gate.
- If the eigen-trace reconstruction fails the Mellin-inversion gate, STOP
  and write the stall report — do not compute moments on unverified traces.

## Verification & acceptance

- `python3 research/scripts/q32_joint_tannaka_moments.py` exits 0, prints
  `GATE VERIFIED` for every prime, then the tables.
- Report exists, tables match script output.

## Stall protocol

If stuck > 20 min on the eigen-trace conventions: deliver the script with
the gate section only + a precise stall report in the report file naming
the exact ambiguity (which formula, which line of the sources conflicts).
