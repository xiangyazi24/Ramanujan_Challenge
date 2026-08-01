# SPEC — Problem 3.2 MAIN THEOREM, fully unconditional. Last stand.

Workdir: `/Users/huangx/repos/Ramanujan_Challenge/problems/3.2` (write all outputs here).
Effort: no cap. No time estimates. Do not stop at a decomposition; attack every piece you name.

## 0. THE TARGET (this is the actual goal — not a sub-lemma)

Apéry's sequences for ζ(3):
```
(n+1)^3 u_{n+1} - (34n^3+51n^2+27n+5) u_n + n^3 u_{n-1} = 0,   a_0=0, a_1=6, b_0=1, b_1=5.
d_n := lcm(1..n)^3.   d_n a_n, d_n b_n in Z.
G_n := gcd(d_n a_n, d_n b_n).
```
**PROVE: G_n = e^{o(n)}.**  Unconditionally. That is the whole task.
(Equivalently, per our banked reduction chain: `log G_n = o(n)`.)

Everything below is what we have established, what is dead, and where the live frontier now is. Verify before you use — one "fresh theorem" was refuted tonight by exactly the kind of check you should be running (see §5).

## 1. Banked reduction chain from G_n to mod-p counting (verify each link in the repo before relying on it)

- `G_n G_{n+1} | 6 d_n^3 d_{n+1}^3/(n+1)^3` ⟹ **v_p(G_n) ≤ 6 for all sqrt(n) < p ≤ n** (ledger `research/working_notes/FABLE_NOTES_energy_bootstrap.md` ~line 622).
- Sieve/master-sum: **log G_n ≤ 6 M(n) + O(n^{2/3} log n)**, where
  `M(n) = sum_{r<n/2} log rad_{>max(sqrt n, r)} gcd(b_r, n-r)`; **(MG): M(n)=o(n) <=> the conjecture** (same ledger, ~line 627).
- Equivalent pointwise form (cron ledger `CRON_FRESH_EYES_pointwise.md` line 18): *the number of prime factors of the single integer b_n lying in (n/2, n] is o(n/log n)*.
- Size/counting arguments are PROVABLY insufficient (7x slack): total log-mass of the window primes is n/2+o(n) while log b_n ≈ 3.5255n; even p^3 | b_n for every window prime costs only 1.5n. Do not attempt any "b_n is too small to hold them" argument (`CRON_FRESH_EYES_pointwise.md` line 20).
- Best unconditional pointwise bound in the literature+ours: `log g~_n ≤ 5.2467 n` (vs trivial 6.5255n), via irrationality measure μ; reaching o(n) that way needs μ→2 at scale Q=b_n — hopeless. **The conjecture is a statement about the arithmetic nature of the Apéry fractions, not about approximation strength.**

**The mod-p object.** For a prime p, in the window r ∈ {1..p-2} (nonwrapping) put
```
u_r = (b_r, c_r) in F_p^2 ,   c = the second solution (c_0=0, c_1=6)  [c is 'a' normalized]
xi_r = (b_r : c_r) in P^1(F_p)   -- the projective Apéry orbit
Z_p  = { r : b_r == 0 mod p }    -- the b-zero set; this is the fibre xi = (0:1)
Delta_{r,h} = b_r c_{r+h} - b_{r+h} c_r = det(u_r, u_{r+h}) = N_h(r) / prod_{j=1..h}(r+j)^3
```
`Delta_{r,h} == 0 mod p  <=>  xi_{r+h} = xi_r` (fibre coincidence; the relation is an EQUIVALENCE, so the collision graph is a disjoint union of fibre cliques). `deg N_h = 3(h-1)`, `N_h` squarefree over Q for all h (proved), so `C_d := |Z_d| = #{r : N_d(r)==0 mod p} <= 3(d-1)` unconditionally.
Records: `|Z_p| <= 3 p^{2/3}` (pointwise), `E(p) := sum_v m_v^2 << p^{5/3}`, energy exponent **3/2 unconditional** (i.e. `F_p << p^{3/2}`), target 4/3, limit E1 = `O(p)`.
Multiplicity/window law (proved): every fibre v and interval J satisfy `m_J(v) - 1 <= 4 span(J)^{2/3}`. Restart bound `Z(H) <= 1 + (3/2)H^{2/3}`.

**Why the energy/collision line feeds the main theorem:** the fibre multiplicities m_v are exactly the mod-p degeneracies that produce common factors; `E(p) >= m_v^2` for each fibre, so improving the energy exponent improves the pointwise zero count `|Z_p|`, which is what the master sum M(n) integrates over p. Breaking the 3/2 energy exponent (any unbounded L, see §2) is the currently-identified gateway. **If you see a shorter route from `log G_n = o(n)` that bypasses the energy line entirely, take it — the energy line is a means, not the goal.**

## 2. The energy gateway, exactly

Mesoscopic scale `D = sqrt(N) L`, `N = p-2`, `L = L(p) -> infinity arbitrarily slowly`. Any ONE of these breaks the 3/2 record:
```
(4.8)    P_D Q_D << N^2
(4.9)    P_D << N/(L^2 log D)
[FR_eta] P_D << D^{2-eta}  for some fixed eta>0   (equivalently |union_{D/2<d<=D} Z_d| << D^{2-eta})
(7.1)    P_{D,1} Q_D << N^2  plus the analogous k_r>=2 core estimate
```
`S_D` = collision pairs at gap <= D (edge count), `Q_D` = triangle count, `P_D` = primitive (first-return) pairs, `P_{D,1}` = #{r : restart index 1}. Unconditional: `S_D <= (3/2)D(D-1)`, `Q_D <= min(66 D^2(1+log D), (sqrt2/3) S_D^{3/2})`. Empirically `P_D ~ 1.4 D` (linear; all five empirical layers sit at the random floor).

## 3. THE FRESHEST LEAD (produced and machine-checked tonight — this is where I would start)

### 3.1 Chart-free determinant (symplectic-bilinear) reduction — repairs the hole that killed the earlier version

Earlier tonight we tried a scalar-phase route using the affine chart `u_r = c_r/b_r`. An audit found a genuine hole: the chart misses the infinity fibre `Z_p = {b_r = 0}`, and the worst-case chart correction (`O(N D^{1/3})` from the window law alone) SWAMPS the target. Note the irony and the signal: **that hole is literally the main conjecture's own object** (`|Z_p|`), so any chart-based route is conditional on the thing we are trying to prove.

The repair: never leave `P^1`. The collision indicator is an algebraic equation, so use additive orthogonality on the determinant itself:
```
T(D) := #{(r,u) : u-r in (D/2, D], xi_u = xi_r}
      = sum_{(r,u) in strip} 1[ det(u_r, u_u) = 0 ]
      = (1/p) sum_{r,u} W(u-r)  +  (1/p) sum_{t != 0} B(t),
B(t)  := sum_{r,u} W(u-r) e_p( t * det(u_r, u_u) ),   W = (majorant of the) strip indicator.
```
This includes the infinity fibre automatically — no chart, no side condition. Main term `~ N D / p ~ D`, which is the empirical truth (`P_D ~ 1.4D`).

**Exponent budget (do this bookkeeping yourself, it is the crux):** the error is `~ mean_{t != 0} |B(t)|`, trivially `<= ND = N^{3/2}L`. To get `[FR_eta]` we need
```
mean_{t != 0} |B(t)|  <<  D^{2-eta} = N^{1-eta/2} L^{2-eta}.
```
- Square-root of the number of terms is `sqrt(ND) = N^{3/4} sqrt(L)` — comfortably BELOW the requirement for eta < 1/2. There is a genuine margin of `N^{1/4 - eta/2}`.
- What we can currently PROVE by completion (bounding the strip-restricted sum by the complete 2-variable sum) lands at exactly `~ p log p = N log p` — i.e. exactly the trivial threshold, short of the target by precisely `N^{eta/2}`. **Any power saving beyond the completion bound closes it.**
- Measured (this session; the scripts are in this workdir at `scratchpad_laststand/det_bilinear.py`, `vecweyl.py`, `exp3_lowdeg.py` — run them, then extend to larger p): at p=1009/2003/4003 with D=2sqrt(N), `mean_{t}|B(t)| = 148 / 195 / 462` versus `sqrt(#pairs) = 171 / 290 / 494`. So `B(t)` sits AT square-root, i.e. the required bound is empirically true with large margin.

### 3.2 What B(t) decomposes into (all machine-checked tonight)

- Per-gap slices: `B(t) = sum_{d in strip} A_d(t)`, `A_d(t) = sum_r e_p(t Delta_{r,d})` — a complete algebraic character sum. We have the banked theorem **[PER-H-WEIL-4H-1]: |A_d(t)| <= (4d-1) sqrt(p)** for every fixed d, t != 0. This is TRIVIAL at mesoscopic d (`4d sqrt p ~ 4NL > N`) — the documented "conductor grows with the height" death. So per-gap Weil is useless; **cancellation must come from summing over d.**
- Cauchy–Schwarz in the long variable reduces `B(t)` to vector Weyl sums `S(alpha,beta) = sum_r e_p(alpha b_r + beta c_r)`. Measured (`scratchpad_laststand/vecweyl.py`): `max_{(alpha,beta) != 0} |S| = 5.5 sqrt(N)` at p=1009,2003 — square-root with a small Gumbel-type constant. Also found: **the map r -> (b_r, c_r) is exactly 2-to-1** (reflection pairing r <-> p-1-r), so Parseval gives `sum |S|^2 = 2N p^2` exactly. NOTE: the naive C-S bookkeeping (C-S over the long u-variable against a short inner sum) gives back the trivial bound — do the C-S in the right variable / with the right weights, or find a different amplification. This is exactly where the work is.
- The two-variable complete kernel `F(t,xi) = sum_{d} sum_r e_p(t Delta_{r,d} + xi d)` is the same object the parallel (cron) line reached from the pointwise side, and it was measured GREEN there: bulk max|F|/p = 4.16..6.07 for p=101..3001 versus the random Gumbel prediction 3.04..4.00 — a constant ratio ~sqrt(2) (a known mirror-coherence factor), **no p^eta conductor growth** (`CRON_FRESH_EYES_pointwise.md` AS.1). Two independent lines therefore converge on ONE object: a bounded-conductor Weil/Radon-pushforward bound for the 2-variable collision kernel. The named missing lemma on the cron side is exactly "fixed-conductor l-adic sheaf realization of the Apéry transfer cocycle with nontrivial geometric monodromy" (`CRON_FRESH_EYES_pointwise.md` AN.2, AO.4).

**So: the single sharpest question is whether `Delta_{r,d}` (as a function of the two variables (r,d), with an additive twist in d) admits a bounded-conductor sheaf-theoretic realization giving square-root cancellation UNIFORMLY in the mesoscopic range.** The numerics from two independent experiments say yes. Nobody has constructed it. Constructing it — or finding a purely analytic substitute (bilinear-forms / additive-energy / Burgess-type short-sum argument that beats completion by any power) — is the frontier.

## 4. Dead routes — each has a precise death certificate. Do NOT re-derive these.

From `CODEX_STRIKE2_report.md`, `CODEX_LASTWALL_report.md`, `CAMPAIGN_MAP_2026-08-01.md` §2, `THEOREMS_2026-08-01_campaign3_addendum.md` (items 41-49):
1. **Padded-word theorem (PROVED)**: reflection-symmetric abstract words satisfy ALL word/clique-level inputs (reflection, no adjacent equality, `C_d <= 3(d-1)`, the 2/3 window law, the exact `C_2 = 1 + 2*[(-51|p)=1]` law, every primitive/cascade/split/renewal identity) yet have `P_D Q_D/N^2 -> infinity`. **Consequence: (4.8)/(4.9) CANNOT be derived from word/clique inputs.** A proof must use arithmetic ruling out coexistence of a triangle-rich core with a large `k_r=1` primitive stratum. The surviving adversary is the SCATTERED one (variable gaps, no long run in any single Z_d, private fibres).
2. **Twelve banked arithmetic inputs audited against that adversary; none bites**: reflection/parity; exact h=2 (-51) split; leading-coefficient apparition; adjacent+shifted resultants; N_h squarefree; Morse/Sp_{4h-2} certificates through h=32; per-height (4h-1)sqrt(p) sums; mesopair diagnostics; transfer codegree/annealed mixing; quarter-value theorem; self-twist rigidity; honest triple poles + apparition mod 24. Reasons are per-input and precise (fixed-height inputs are negligible as D→∞; annealed mixing misses the quenched ordered clock; single-row modular pins have no propagation).
3. Katz equidistribution; T-adic/eigencurve; Stickelberger valuation; Galois orbit norms; HGM/Hasse degrees; mass formula; black-box bounded-conductor theorem (constant-sheaf counterexample); dispersion k=2; sign factorization (q23) — nine certificates in `CAMPAIGN_MAP_2026-08-01.md` §2.
4. Low-genus collapse (genus is `6hk-3h-3k-2 = Theta(hk)`, so pairwise Weil stops at p^{1/6}); h-algebraization in char 0; sheaf triangle induction; van der Corput transfer length; master variety; fixed-index moments (loses `H^{q-1}`); rowwise Schur/L1; combinatorial de-log (`Q_H = Theta(H^2 log H)` is sharp abstractly); bottom-scale propagation.
5. Montes/structural-prime route for the all-h residual: structurally impossible (local types not uniform).
6. **Tonight's additions:**
   - **[NO-RUN] REFUTED** (`CODEX_NORUN_report.md`, ERRATA 2026-08-02): the claim `xi_{r+1} = M_r xi_r` for the two-solution row is FALSE — the companion matrix propagates the two-time state of ONE solution, not the row `(b_n : c_n)`. Live counterexample `p=997, r=248, d=182`. Any eigenvector/run-rigidity argument built on that action is void. Correct orbit-free statements that DID survive: `U_r = [u_{r-1}; u_r]` is invertible with `det U_r = 6/(r+1)^3` (proved, telescoping), collision at (r,d) `<=>` the (2,1) entry of `T_{r,d} = M_{r+d-1}...M_r` vanishes, and a 2-run at gap d `<=>` `p | Res_x(N_d(x), N_d(x+1))`.
   - **[DEAD-MAX-WEYL-CLASS]** (ERRATA 2026-08-02, second entry): for the SCALAR-phase (affine chart) formulation, the max-route needs `delta >= 1/2 + eta/4` (beyond square-root, while the truth is exactly square-root: measured `max_t |sum_r e_p(t H_r)|` and its Apéry analogue at 3.5-4.8 sqrt(p)), and the Parseval-in-t route is an identity (zero information). So no single-variable max-Weyl input can close [FR_eta]. §3 above is the response: go bilinear/2-variable and chart-free.
   - **Low-degree collision variety: negative** (`scratchpad/exp3_*.py`): at p=1999, bidegree (10,10), the kernel of polynomials vanishing on all strip collision points is 18-dimensional and consists ENTIRELY of multiples of the mirror line `2r + d = p-1` (exact accounting at two primes; kernel elements are nonzero on unseen non-mirror collisions; random control gives kernel 0). So the collision set carries no low-degree component beyond the known reflection family — matching the extremal verdict that visibly arithmetic-forced families give only `Theta(D)`.
   - **Harmonic toy calibration**: the same strip-anti-concentration target for the pure harmonic walk `H_r = sum_{n<=r} 1/n mod p` is OPEN in the literature (value distribution and reciprocal-congruence results exist — Garaev-Luca-Shparlinski 2005 and the Kloosterman/bilinear machinery — but nothing gives the localized distance-collision count). Stepanov on the harmonic phase fails at a precise point: the k-step accumulated increment has k poles, so multiplicity gain is paid back by pole growth; the mod-p^2 lift (Wolstenholme/Lehmer/p-adic Gamma) supplies no bounded-complexity substitute for the multiplicative logarithmicity that powers Heilbronn sums. **So the wall is not Apéry-specific** — it is a general frontier about distance distributions of collisions of arithmetically-defined walks. That cuts both ways: a solution here is a genuine analytic-number-theory advance, and imported technology alone will not do it.

## 5. Discipline (non-negotiable — tonight's refutation came from violating the first line)

- **Verify every claim you use and every claim you produce, by machine, before banking it.** A symbolic identity is not verified until you evaluate it on the live orbit at several primes. If you write a verification script with a placeholder branch, that branch is NOT a check.
- State every result with its exact status: `PROVED-all-h` / `PROVED for p > p_0` / `VERIFIED-N (which N)` / `CONDITIONAL(exact hypothesis)` / `EMPIRICAL`. A conditional result whose hypothesis you cannot exhibit an instance of is not a milestone; say so.
- If you refute something in this spec, say so loudly at the top of your report. That is the most valuable outcome available short of the theorem.
- Do not report "nearly there". Report state: what is proved, what is open, what the exact next obstruction is.
- No effort ceiling. If a route stalls, write the exact blocking identity/inequality and move to the next; do not stop early.

## 6. Deliverables (write into the workdir)

1. `CODEX_MAINTHM_report.md` — terminal verdict FIRST (proved / partial with exact statement / precise obstruction per route), then full derivations with every constant and every degree bookkept.
2. `CODEX_MAINTHM_verify.py` — one self-contained script; every numerical or symbolic claim in the report is a gate in it; ends with a single `PASS`/`FAIL` line. Reproduce (and extend to larger p) the two measurements in §3: `mean_t|B(t)|` versus `sqrt(#strip pairs)`, and `max|S(alpha,beta)|` versus `sqrt(N)`.
3. Any new unconditional lemma — however small — stated in inventory style with proof + verifier. Bank every brick: a proved lemma about `B(t)`, about the sheaf realization, about `M(n)`, about `|Z_p|`, all count.
4. If you close `[FR_eta]` (or any of (4.8)/(4.9)/(7.1)): immediately carry it through the chain to `log G_n = o(n)` and state exactly which links of §1 you used and whether each is proved in the repo or needs its own proof.

## 7. Suggested order of attack (yours to override with reason)

1. Re-derive §3.1's reduction and its exponent budget independently; confirm or refute the claim that the completion bound lands exactly at the trivial threshold and that any power saving beyond it closes `[FR_eta]`. This is the single most load-bearing new claim in this spec.
2. Attack `mean_t |B(t)| << p^{1-kappa}` for some `kappa > 0`. Angles: (a) the correct Cauchy–Schwarz / amplification (the naive one returns the trivial bound — find the right one); (b) bilinear forms with the symplectic phase `det(u_r, u_u)` and additive energy of the orbit point set `{(b_r, c_r)}` in `F_p^2` (note the exact 2-to-1 reflection structure, and that the point set has `N` points in a `p^2` plane — incidence/sum-product technology in `F_p^2` is the natural toolbox: Rudnev, Stevens–de Zeeuw, Shkredov); (c) the sheaf-theoretic route: construct the bounded-conductor l-adic realization for the 2-variable kernel `Delta_{r,d}` with the additive twist in `d`, which both empirical lines support.
3. In parallel, re-examine §1: is there a route from `log G_n = o(n)` that needs something WEAKER than the energy gateway — e.g. an average-over-p statement, where the bad primes are allowed to be exceptional as long as their total log-mass is `o(n)`? The master sum `M(n)` is an average over primes; we have been chasing a per-prime theorem. **An `o(n)` global bound may tolerate a sparse set of bad primes.** Quantify exactly how much per-prime failure the master sum can absorb; this could be the cheapest real win available, and to my knowledge it has not been pushed.
