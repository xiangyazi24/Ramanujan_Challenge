# CODEX SPEC — proof.tex new section: the orbit-energy theorem & E1 program (max)

## Mission

Draft `FABLE_SECTION_orbit_energy.tex` — a new self-contained section for
proof.tex presenting the final-strike theorem suite. Submission is TODAY;
quality bar = referee-ready statements with complete elementary proofs for
the proved items and scrupulous honesty labels for the empirical program.

## Sources (read in this order; cite by ledger § in comments)

1. `problems/3.2/THEOREMS_2026-08-01.md` — the inventory (items 5-11, 19-21).
2. `research/working_notes/FABLE_NOTES_energy_bootstrap.md` §88–§109 — the
   final-strike ledger (statements, proofs, verification records).
3. `problems/3.2/chatgpt-answers/Q6506.md` §3 (block-energy proof), §1
   (spectral identities), §5 (GPRV); `Q6496.md` (renewal identity, QRLL);
   `Q6500.md` (audit corrections: domain u≠-1, R_restricted, tail (p-3)/(H+1)).
4. `problems/3.2/DOSSIER_E1_EMPIRICAL.md` — the empirical table (cite, don't
   duplicate).

## Content plan (LaTeX, match proof.tex's existing style/macros — read its
preamble first; NO new packages)

1. **Setup**: the companion solution c_n (c_0=0, c_1=1, Casoratian
   W_n = -1/n^3), the projective orbit pi(n) = [b_n : c_n], the multiplicity
   function, and the two energy conventions (full F_p = sum m(v)^2 vs
   off-diagonal E_p = F_p - (p-2)) with the reconciliation identity.
2. **Structure theorems (with proofs)**:
   (a) restart identity and the gap-polynomial dictionary N_m = U_m
       (state with the corrected conventions from Q6500: domain u != -1,
       restricted count; include the one-step-shift remark);
   (b) reflection: pi(p-1-n) = pi(n) pointwise (from strong reflection,
       cite the existing proof.tex section that has strong reflection);
   (c) center-vanishing: odd m => U_m(-m/2) = 0 over Q (one-line proof
       from the reflection identity U_m(-s-m) = (-1)^m U_m(s));
   (d) root pairing s <-> -s-m;
   (e) renewal identity N_{m+g+1}(s) = N_{m+1}(s)N_g(s+m+1)
       - (s+m+1)^6 N_m(s) N_{g-1}(s+m+2) (machine-verified; give the
       two-line induction proof via the matrix product split);
   (f) adjacent coprimality gcd(N_m, N_{m+1}) = 1.
3. **Main theorem [ENERGY-3/2]**: F_p <= C p^{3/2} (C explicit, e.g. 2.45
   — RE-DERIVE the constant carefully from the block argument; do not
   copy blindly): full proof — block partition of {1..p-2} into ceil(M/H)
   blocks, Cauchy-Schwarz over blocks, in-block collisions have gap < H,
   gap-h collision count <= deg U_{h-1} = 3(h-1) by the banked
   nonvanishing (cite proof.tex's nonvanishing section + nv_theorem.tex),
   optimize H. State the corollary positioning honestly: this is an energy
   theorem; energy < p^{4/3} would be needed to beat the pointwise 2/3.
4. **The sigma ladder and the E1 program (no overclaim)**: the exponent
   machine Q_p(Delta) << Delta^{3-sigma} => F_p << p^{2-1/(3-sigma)};
   sigma = 0, 1, 2 give 5/3, 3/2, 1; our theorem sits at sigma = 1; E1
   (= sigma 2) is equivalent to the family-compatibility statement. Present
   the three PARALLEL SUFFICIENT open targets with exact quantifiers
   (CRITICAL CORRECTION: write them as parallel sufficient conditions /
   open targets, NOT as equivalent — bidirectional counter-models exist
   (Q6515 sec.6); saying 'equivalent' is referee-bait): [ATR]
   (one-line inequality), [GPRV] (Chebotarev root-count variance of the
   family {M_h}), [2D-SQRT] (the bilinear form). Empirical status: cite
   the dossier (theta = 0, C ~ 1.08; four laws closed to 4 decimals) —
   clearly labeled as numerical evidence, NOT theorems.
5. **Small-gap splitting laws**: R_2 = 1 + 2[(-51|p)=1] (N_2 = P(r+1),
   disc -51) and the remark that the same field Q(sqrt(-51)) governs the
   annealed two-step return probability — the thin-exceptional-family
   stratification note.

## Hard constraints

- New file only: `FABLE_SECTION_orbit_energy.tex`. Do NOT modify
  proof.tex itself, FABLE_SECTION_value_distribution.tex,
  FABLE_SECTION_apparition_tower.tex, ERRATA.md, CAMPAIGN_MAP*.
  (I will add the \input line and build after review.)
- Every theorem labeled; every proof complete or the item explicitly
  labeled as (numerically supported) conjecture/program — no middle
  ground, no "it can be shown".
- Every constant re-derived, not copied. If your re-derivation of the
  3/2 constant differs from 2.45, SAY SO prominently at the top of the
  file in a LaTeX comment.
- Notation must match proof.tex's existing macros (read the preamble
  and §15/§16 for house style: theorem environments, labels, citation
  style for internal cross-references).
- Length target: 6-10 pages compiled.

## Acceptance

- File compiles standalone when \input into proof.tex (test with a local
  copy: cp proof.tex /tmp/proof_test.tex, append the \input before
  \end{document}, PATH="/Library/TeX/texbin:$PATH" pdflatex — report the
  result, but do NOT commit the modified proof.tex copy).
- A summary comment block at the top listing: theorems included, proofs
  complete y/n each, constants re-derived y/n, discrepancies found.
