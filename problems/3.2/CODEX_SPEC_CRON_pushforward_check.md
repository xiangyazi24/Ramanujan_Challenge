# CODEX SPEC: numerical verification of the rank-18 pushforward package + self-twist exclusion

## Goal
Ledger `CRON_FRESH_EYES_pointwise.md` appendix P.3 records the Q6375 monodromy package for G = t_!(F ⊗ F) as [未验证·结构合理]. Turn as much of it as possible into machine-verified facts, and produce the numerical substrate for [GAP-2] (Tannakian group — the one attackable theory gap; see appendix R.4/S.7).

Background facts already machine-verified (do not re-derive, cite the ledger):
- t(x) = x(1-8x)/(1+x), branch discriminant t^2 - 34t + 1, roots (1±√2)^4 = Apery ODE singular points.
- Mellin pair: b_r ≡ -Sum_x H_p(x)^2 t(x)^{-r} + Sum_t A_p(t) chi_2((t-1)^2 - 32t) t^{-r} (mod p), where H_p is the truncated Franel-square-root polynomial with A_p(t(x))(1+x)^{p-1} = H_p(x)^2 in F_p[x] (identity verified for p = 13,17,29,37,41).
- [FALSE-3F2] (appendix S.7): Franel local system is RANK 2 (2F1(1/3,2/3;1)); the life-side campaign verified a_p(E_u) = H_p(u) on 1084 fibers where E_u is the Beauville-IV elliptic family — find their scripts in problems/3.2/ (grep for q6372, beauville, fiber files; also FABLE_NOTES sections 44-47 and 53-56 name them) to get the exact Weierstrass model of E_u. If you cannot locate the model, derive the fiber counts directly: define T_F(x) for x in F_p, x not a pole, via T_F(x) = H_p(x) lifted to the integer in (-p/2, p/2) — justified by a_p(E_x) = H_p(x) mod p and |a_p| < 2 sqrt(p) (valid once p > 16; note H_p(x) is the mod-p trace, the integer lift is unique).

## Tasks
For a prime set P = {29, 37, 41, 53, 61, 73, 89, 101} (skip primes where any construction degenerates; say so):
1. Build the trace function of G on the t-line: T_G(t) = Sum_{x : t(x) = t} T_F(x)^2 (integer values, sum over the 0, 1, or 2 preimages; record the fiber count N(t) = 1 + chi_2((t-1)^2-32t) and verify it matches the actual preimage count for every t — this re-verifies the fiber-count decomposition at scale).
2. Verify the Mellin inversion at the level of G: check numerically for each p that -Sum_t T_G(t) t^{-r} + (the chi_2-correction written with T_F-squared data) reproduces b_r mod p for all r < p. (This is the sheaf-level restatement of the verified K.2 identity — it must pass; if it fails, diagnose the normalization (1+x)^{p-1} factor and report the corrected form.)
3. Self-twist exclusion (the [GAP-2] substrate): for every multiplicative character chi of F_p^* of order dividing 24 and every character of order dividing any d | p-1 with d <= 30 (i.e. all chi with ord(chi) <= 30), test whether T_G(t) * chi(t) = T_G(t) for all t with T_G(t) != 0 — equivalently whether T_G is supported on {t : chi(t) = 1}. Report every chi that passes for every p in P (candidate self-twists). Expected outcome per Q6375: only the trivial one on this test (the known quadratic partner is G ⊗ L_{chi_2(q(t))} which twists by a character OF THE ARGUMENT q(t), not of t — verify separately that T_{G⊗L}(t) = chi_2((t-1)^2-32t) T_G(t) defines the companion and that IT also has no small-order self-twists).
4. Duality/inversion symmetry: verify numerically that the Mellin transforms satisfy M(r) = M(-r - c) for the appropriate fixed shift c predicted by b_{p-1-r} = b_r (find the exact c from the data, state it) — this is the G ≅ G^dual symmetry made concrete.
5. Rank probe (best effort): the L-function degree of G ⊗ L_chi for generic chi — compute S_m(chi) = Sum_{t in F_{p^m}^*} T_G^{(m)}(t) chi(Norm(t)) for m = 1 only if easy... if extension-field data is too heavy, SKIP and say so; the rank-18 claim then stays [待验] and that is an acceptable outcome. Do not burn hours here.

## Outputs
CRON_pushforward_check_report.md with: per-prime pass/fail per task, every candidate self-twist found, the fixed shift c, and a clear final block: which P.3 claims are now VERIFIED / REFUTED / still 待验.

## Sanity gates
- Task 1 fiber counts must match chi_2 formula for every t, every p (this re-derives a verified identity — abort loudly if it fails).
- Task 2 must reproduce b_r mod p exactly for all r < p, at least 3 primes, before tasks 3-5 run.

## Discipline
- New files only, CRON_ prefix: CRON_pushforward_check.py (sympy/pure python fine; this is small-p work), CRON_pushforward_check_report.md.
- Work only in problems/3.2/ (read-only access to life-side scripts is fine). Commit when done, report SHA.
- Do not touch CRON_FRESH_EYES_pointwise.md (single-writer: cron session owns it).
- Do NOT call the ChatGPT bridge / ask-gpt.py or dispatch sub-questions to any oracle. Pure computation task. If blocked, write the blocker into the report and stop.
