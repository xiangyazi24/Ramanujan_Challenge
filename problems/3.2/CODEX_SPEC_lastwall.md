# CODEX SPEC — THE LAST WALL: supercritical moment normalization (single deep strike)

## The wall, exactly
N = p-2. C_d = #{r <= N-d : N_d(r) = 0 mod p} (gap-d collisions; recurrences in
campaign3_questions/CTX.txt). First and second window moments
  S_D = sum_{d<=D} C_d,     Q_D = sum_r binom(d_D(r),2),  d_D(r) = #{d<=D: collision at (r,d)}.
PROVE, for SOME unbounded L = L(p) -> infinity at D = sqrt(N)*L (supercritical scale):
  [MESO-TOTAL]  S_D << N        (weakest sufficient)
  or [MESO-PAIR] Q_D << N       (implies TOTAL via S_D^2 <= N(S_D + 2Q_D))
EITHER gives E_p << N^{3/2}/L + N: ANY unbounded L breaks the 3/2 record = the
campaign's win condition. Even L = log log p counts. Partial: S_D << N at D = sqrt(N)*L
for L up to some explicit slowly-growing function.

## Why this is the last wall
- Purely combinatorial optimum (proved today, exponent optimal): Q_H <= 22C H^2(1+log H)
  from R_d <= Cd. At D = sqrt(N)L this gives N L^2 log N — off by exactly L^2 log N.
  The inverse-square clique-energy identity Q_H <= 22 H^2 sum_{d<=H} R_d/d^2 shows small
  gaps dominate that route: no purely-worst-case-R_d argument can pass. Arithmetic
  input about THIS orbit is REQUIRED.
- Empirical truth (banked, machine-exact): C_d has mean ~1.4 (flat), S_D ~ 1.4 D,
  Q_D ~ 0.28(D + D^2/p) — i.e. at D = sqrt(p) log p the true Q is ~ log^2 p, MASSIVELY
  inside the target N. The room is a factor ~N/L^2. Do not be conservative.

## The full banked arsenal (all machine-verified or audited; use freely)
1. Gap dictionary: collisions at (r,d) <=> N_d(r)=0, N_d(r)=U_{d-1}(r+1), deg N_d=3(d-1);
   transfer product G_m(s) (2x2), det = prod_{j=1}^m (s+j)^6.
2. Addition law: N_{h+d}(x) = N_d(x+h)N_{h+1}(x) - (x+h+1)^6 N_{d-1}(x+h+1)N_h(x).
   Corollary: a common root of (N_a(x), N_g(x+a)) is a root of N_{a+g}(x).
3. Endpoint law: N_h(-j) = (-1)^{j-1}((j-1)!)^3 b_{j-1} ((h-j)!)^3 b_{h-j} (Apery products).
   Adjacent resultant: Res(N_h, N_{h+1})-type A_m = +-prod_{j<=m}((j!)^3 b_j)^6;
   Res(N_h,N_k) = A_{h-1} * Res(N_h(X), N_{k-h}(X+h)).
4. Restart/continuant machinery (proved in paper): |Z_p| <= 1 + R_rel(p,Y) + (p-3)/(Y+1)
   with R_rel(p,Y) = sum_{m<=Y} #{s: U_m(s)=0} — NOTE R_rel is exactly an S-type object;
   the paper's |Z_p| <= 3p^{2/3} took Y=p^{1/3} with the trivial R_rel <= (3/2)Y^2.
   QUESTION TO PUSH: the restart argument runs along the orbit consuming zeros; can it be
   REVERSED or ITERATED to bound S_D itself (each gap-d collision is a restart of the
   companion solution started at r; multiple collisions from nearby bases force U-zero
   patterns that the renewal identity constrains)?
5. Per-gap exponential sums (verified audit): |sum_r e_p(t delta_d(r))| <= (4d-1)sqrt(p),
   delta_d = N_d/prod(r+j)^3. Chebotarev fixed-d: mean C_d -> 1 + kappa_d (kappa_d=[2|d]).
   Reflection: C_d's mirror collision at r=(p-1-d)/2 for even d. h=2 law: roots of N_2
   governed by (-51/p). Apparition/leading-coefficient laws banked in paper section 16.
6. New today (cron): [THM-NEAR-WALL-3/7]: K_p(H,D) << H D^{11/7}; W-injection lemma
   sum binom(W_H,2) <= Q_p(2H). Also [Q-H2LOGH] as stated above.
7. Anti-concentration reformulation: Q_D <= (max_r d_D(r)) * S_D / 2, and
   S_D^2 <= N(S_D + 2Q_D). So EITHER of:
   (i) max_r d_D(r) <= m_0 with m_0 * S_D << N,
   (ii) S_D << N directly,
   suffices. A pointwise multiplicity bound max_r d_D(r) << D/L^2-type at D=sqrt(N)L is
   another face: a single base r with many returns within a window of length D means the
   companion solution restarted at r has many zeros in [r, r+D] — the RESTART bound
   applied to the SHIFTED solution bounds exactly this: work out what the restart
   machinery gives for the zero count of the solution y^{(r)} (0,1-start at r) in a
   window of length D: if it gives << D/(Y+1) + R_rel-type with the SAME structure as
   |Z_p|, then max_r d_D(r) <= 3D^{2/3}-type FOR EVERY r — check: does the paper's
   restart proof depend on the specific solution b (fiber over [0:1]) or does it work
   verbatim for every solution/fiber?? The paper proof used b_z=0 => x_z=infinity =>
   restart at 0 — for a general fiber the restart lands at a MOVING point — BUT for the
   companion-at-r solution the zero set IS the collision set from base r — think hard
   here; a per-base restart bound d_D(r) <= c D^{2/3} would give
   Q_D <= c D^{2/3} S_D / 2 and with S_D <= (3/2)D^2 trivially: Q << D^{8/3}... at
   D = sqrt(N)L: N^{4/3}L^{8/3} — not enough alone; but combined with S_D << N it gives
   Q << N D^{2/3}: still big. The winning combination is (ii) S_D << N — concentrate there.
8. S_D << N attack surfaces: S_D = sum_{d<=D} C_d = total zeros of the family
   {N_d mod p}_{d<=D} in windows. Chebotarev heuristic mean = (1+kappa_d) per d — sum
   ~ (3/2)D << N needs D << N i.e. TRUE for all D < N/log-ish IF mean holds uniformly.
   The obstacle is per-p uniformity (family compatibility). BUT S_D << N is much weaker
   than mean: it allows average C_d up to N/D = sqrt(N)/L — i.e. average root count up
   to sqrt(p) per polynomial! Weil per-d gives C_d <= complete-count: the COMPLETE root
   count R_d^alg of N_d over F_p satisfies R_d^alg = (p - stuff)/p*deg-free bound via
   the exponential sums: R_d^alg = sum over F_p of 1_{N_d=0}... use the standard
   squarefree/character detection: R_d^alg <= 1 + (deg-1)... no — Weil for the CURVE
   y^2 = N_d-type? Direct: R_d^alg <= deg = 3(d-1): trivial. Character sum detection of
   roots does NOT beat degree. The known nontrivial global input: Chebotarev AVERAGE over
   p. At fixed p: the only sub-degree root-count bounds come from the restart/renewal
   structure (point 4) — the paper DID prove R_rel-based bounds for the b-fiber giving
   |Z_p| <= 3p^{2/3} << deg-sum. EXTEND: the same renewal telescoping over the FAMILY
   d<=D: sum_{d<=D} C_d — the renewal identity Q_D = sum_{a,g>=1,a+g<=D} C_p(a,g)
   (banked) ties the family's zeros together. Grind this: telescope S_D through the
   addition law: every gap-(a+g) collision at r either has an intermediate collision
   (contributing to Q) or is "primitive"; primitive collisions inject into ... build the
   primitive-decomposition: S_D <= (primitive count P_D) + f(Q_D) and bound P_D by the
   renewal/coprimality structure. If P_D << N unconditionally (primitive gap-d collision
   = both N_a and shifted N_g nonzero for all intermediate splits — strong constraints
   via adjacent coprimality), then S_D << N + f(Q_D) and combined with Q_D <= 22D^2 log D
   * (bootstrap): iterate S->Q->S. Work out whether the system closes.

## Rules
This is the campaign's decisive single target. No effort ceiling, no time estimates,
no early stop. Try MULTIPLE distinct attack vectors (at minimum: restart-for-general-
fiber (7), primitive decomposition + bootstrap (8), inverse-square with arithmetic R_d
input on the SMALL-d segment (R_d for d <= polylog via explicit Chebotarev/apparition:
the first log D worth of d's contribute O(sum_{small d} R_d) — if R_d <= C uniformly
for d <= (log N)^A via explicit small-d certificates, the clique bound improves to
Q <= 22 H^2 [sum_{d<=Y0} R_d/d^2 + 3 sum_{d>Y0} 1/d] = 22H^2[O(1)... no: this needs
R_d = o(d) for MOST d — formulate the exact average condition and check what apparition
gives]). Every claimed step must be written out fully; verify any identity numerically
in python before relying on it. Machine-check small cases of any new inequality
(p in {199, 499, 997}, D up to p^{0.6}: compute S_D, Q_D, max_r d_D(r) exactly and
test your intermediate bounds against truth).
Report to CODEX_LASTWALL_report.md: every vector tried, exact statements proved
(PROVED / VERIFIED-N / CONDITIONAL / DEAD+reason), and the sharpest unconditional
theorem achieved. If you prove S_D << N or Q_D << N at ANY unbounded L: say so in the
first line in capitals.
