# CODEX SPEC — last-wall empirical ground truth (companion to CODEX_SPEC_lastwall.md)

Report to `CODEX_LASTWALL_EMPIRICS_report.md`. Exact integer arithmetic. Definitions in
campaign3_questions/CTX.txt. For p in {997, 1999, 4001, 7919} and D in
{ceil(sqrt(p)*log(p)), ceil(p^0.6), ceil(p^0.66)}:
1. Compute exactly: S_D = sum_{d<=D} C_d; Q_D = sum_r binom(d_D(r),2); max_r d_D(r);
   the full distribution (deciles + top-10) of d_D(r).
2. [Vector-7 premise] Does max_r d_D(r) grow like c*D^{2/3}, like log, or bounded?
   Report max_r d_D(r) vs D^{2/3} vs log(p) per cell. Also identify WHICH r attain the
   max (are they in Z_p, i.e. zeros of b? near the mirror centers? small r?).
3. [Vector-8 premise] Primitive decomposition: call a collision (r,d) PRIMITIVE if there
   is no d' < d with collisions at both (r,d') and (r+d', d-d'). Compute P_D = #primitive
   collisions vs S_D per cell. Verify the renewal claim empirically: every non-primitive
   collision splits; measure the split multiplicity distribution.
4. Small-d segment: R_d exactly for d <= ceil((log p)^2): is R_d bounded by an absolute
   constant there (report max and mean), i.e. is the "R_d = O(1) for polylog range"
   input empirically available?
5. Verdict lines for each of 2/3/4: which premises of the deep strike hold.
Pure python3, efficient (O(pD) via b,c recurrences + Delta identity). No early stop.
