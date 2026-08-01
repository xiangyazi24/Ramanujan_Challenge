# CODEX SPEC: BREAK THE WALL (single deep strike, no effort cap)

## Mission

Prove ANY statement that breaks the 3/2 energy record for the Apery projective orbit. This is a mathematics-proof task, not a computation task: your deliverable is a rigorous proof (plus machine verification of every new identity/inequality you introduce). You have the complete banked arsenal below; the graveyard section tells you what is already dead so you do not re-till it. Read the referenced files in this directory as needed — the ledger `CRON_FRESH_EYES_pointwise.md` (appendices AR-AT, especially AT.1-AT.27) is the authoritative account.

## Objects (nonwrapping throughout)

Orbit pi(r) = [b_r : c_r] in P^1(F_p), r = 0..p-2; collision pi(r)=pi(r+h) <=> Delta(r,h) = b_r c_{r+h} - c_r b_{r+h} = 0 <=> N_h(r) = 0 (gap polynomial, deg N_h = 3h-3).

```
R_h = #collisions at gap h;  S_D = sum_{h<=D} R_h;  d_D(r) = per-base window count;
Q_D = sum_r C(d_D(r),2);  E = sum_v m_v^2 (energy);  N = p-1.
```

## WIN CONDITIONS (any one suffices; L = L(p) -> infinity arbitrarily slowly is enough)

At D = sqrt(N) L:
- W1 [MESO-TOTAL]: S_D << N.
- W2 [MESO-PAIR]: Q_D << N.
- W3 [MESO-S1-2/3]: S_D << N^{2/3}  (implies W2 via the spectral bound).
- W4 [MIDPOINT-AC']: balanced hinge count sum_{u; h,d > G_0, h+d<=D} 1[Delta(u-h,h)=0] 1[Delta(u,d)=0] << N with G_0 = cN/(D log D) (implies W2 with the near-wall bound).
- W5: any unconditional E << p^{3/2}/L(p) directly.
- W6 [BLOCK-SUBQUAD]: S_p(J,L-window) << L^{2-eta} for translated blocks at the optimizing scale.
Also accepted (different wall, still a prize): W7: [W-IRRED] or [A_h-primitivity] for all h (characteristic-zero; gives G_h absolutely irreducible all h).

## ARSENAL (all proved/banked today unless marked; use freely)

1. [ABSTRACT-Q-H2LOGH]: any word with R_d <= Cd has Q_H <= 22C H^2(1+log H); log SHARP abstractly (multiscale AP construction attains it); conditional log-removal: prefix S_1(Y) <= A Y^{2-delta} => Q_H <= 22A(3+2/delta)H^2 (Abel).
2. [THM-NEAR-WALL-3/7 + q=2 rerun]: K_p(H,D) << min(H^{2/3}D^2, HD sqrt(l(D)l(H))); W-injection lemma sum_u C(W_H(u),2) <= Q_p(2H); identity sum_u nu_D(u)^2 = S_D + 2Q_D; box conversion audited (straddling boxes split at dyadic boundaries; exponent q' = 3 - 1/(2-eta)).
3. Spectral triangle bound: Q_D <= (sqrt2/3) S_1(D)^{3/2} (collision graph, trace(A^3)).
4. Window multiplicity: m_J(v) - 1 <= 4H^{2/3} (any interval, any value); per-fiber restart identities.
5. sigma=1/2 box theorem Q_p(Delta) <= (27/4)H^{5/2} with master inequality Q_H <= qH S_1(D-1) + ((q-1)/2)S_1(H).
6. Axis-strip reduction: Q_D^{axis}(G) <= 12 D^{2/3} G^2; balanced core = the whole problem.
7. Energy theorem E <= 2 sqrt3 p^{3/2}; S_1(p-2) = (E-N)/2; Cauchy S_D^2 <= N(S_D + 2Q_D).
8. [N-SQFREE] all h (tridiagonal pencil + strip localization: roots in strips -i-1<Re<-i, exactly 3 per strip); [NONCOLLAPSE-ALL]; full Morse certificates h<=30; G_h absolutely irreducible h<=40 (двойной method) => [L2-FREQ]: Collision(h) = 2p + O(h^2 sqrt p) unconditionally for h<=40.
9. Exact laws: renewal K_{m+g+1} = K_{m+1}K_g(X+m+1) - (X+m+1)^6 K_m K_{g-1}(X+m+2); adjacent resultants A_m = +-prod((j!)^3 b_j)^6; pole values N_h(-j) = +-((j-1)!)^3((h-j)!)^3 b_{j-1}b_{h-j}; V_h(0) = -27 lc^4 S^6 B^2 Disc(N_h); det F_n = -((n-1)!)^6; twisted symplectic T^T J T = r^6 J.
10. Riccati framing: collisions at gap h = fiber g_h^{-1}(L_infty) of the convergent map g_h = N_h/N_{h-1} (deg 3h-3, adjacent-coprime); collision <=> S_r = S_{r+h} with the dual line S_r = [-c_r : b_r]. The p-curvature C_p(a) is rank-one nilpotent at every rational point with line = the b-state line (C_p(a) = F_a(-E_12)F_a^{-1}).
11. Empirics (all Poisson+mirror, no anomaly): supercritical Q_D strata flat (data/baseline 0.97-0.99); wallprobe 1296 boxes flat in three gauges, SVD < 3; small-d R_d = O(1) empirically to (log p)^2; midpoint anti-correlation slack K/Cauchy <= 0.05; B_p (Radon zero-push) bounded vs binomial control growing.

## GRAVEYARD (do not re-till; each has a written death certificate in the ledger)

- Fixed-h Weil/Chebotarev at fixed p (quenched quantifier); GRH insufficient for GPRV (off by p^5).
- Raw resultant divisibility counts (cut-edge pollution: p | b_j kills the whole upper triangle).
- Low-rank clock carriers: algebraic + analytic tests double-negative; Mellin character linearization dead (Weyl algebra full); scalar bispectral gauge dead; rank-4 first-block linear extension dead (beta cancels); index-side contiguity dead.
- p-curvature leverage in the wall range h < p (line field has full reduced degree p-1; character trivially nilpotent on rational points).
- lambda=1 slow degree growth alone (no known cancellation theorem for a single nonautonomous cocycle).
- Archimedean dyadic scale separation for critical values (only O(log h) bands); p-adic Newton slopes at ell=5,7 (units/packets); positivity of the noncollapse quadratic (complex variance, signs +,-,-).
- Self-referential bootstraps: master-inequality cascade fixed point q=2 reached; feedback tautological below that; clique method's D^2 window factor immovable by G-separation (REFUTED in ledger); capacity Cauchy is a lower coupling only.
- Fixed-d dispersion (locked-clock two-color correlation residual); annealed two-clock curve (Q6455) does not apply to the locked diagonal.

## MOST PROMISING DIRECTIONS (our current ranking — you may override with reasons)

A. [MESO-S1-2/3] via the census: S_D is a FIRST moment. New inputs available that were never combined: [L2-FREQ] unconditional for h<=40 gives the VALUE-collision second moment Collision(h) = 2p + O(h^2 sqrt p); the frequency side (Parseval sum_{t,xi}|F(t,xi)|^2 = p^2 sum_h Collision(h)) is exact. Is there an inequality from value-collision second moments to ZERO-fiber first moments (R_h is one fiber of the value distribution whose L2 mass is now controlled)? A per-h large-sieve: R_h <= (1/p) sum_a nu_h(a)... work out whether L2-flatness of the value distribution of D_h bounds sum_h R_h nontrivially on any h-range (the h<=40 range is unconditional NOW; can the certificate range h <= h_0(p) be leveraged asymptotically? For W3 you need the sum over h <= sqrt(N)L — the h<=40 piece is O(1) of it — so you ALSO need a mechanism for large h; consider dyadic interpolation between the certified range and the degree bound, and compute exactly what h_0(p) growth rate would tip the balance — then examine whether the certificate METHOD (modular Morse certificates, uniform in p) actually proves [CRIT-2H] for h <= h_0(p) = a power of log p or p^epsilon with EFFECTIVE constants: the certificate conditions are polynomial identities in h; their failure set is a proper subvariety — an EFFECTIVE all-h theorem may be within reach via your W7 route, and W7 => [L2-FREQ] all h => rerun the census question with the full range).
B. [MIDPOINT-AC'] via mirror-skeleton removal + the S/S-dual reformulation (collision <=> S_r = S_{r+h}): the hinge asks that gap-d collision midpoints avoid backward-return-rich colors; in the S-coordinate both events are equalities of the SAME sequence S at different times — the hinge is a THREE-TERM configuration S_{u-h} = S_u = S_{u+d}... wait that is just the triple again; the leverage would be an inequality for triples with one long and one short arm that exploits the Riccati structure of consecutive S-values (S_{r+1} determined by S_r via the inverse Riccati flow). Chains of the Riccati map: the short arm d <= D_0 small means S_u and S_{u+d} are related by a BOUNDED-degree rational map; the collision S_u = S_{u+d} = a fixed point of g_d composed appropriately — fixed points of bounded-degree maps are <= 3d — this is R_d <= 3d again — but the JOINT event with the long arm: r = u-h satisfies S_r = S_u where S_u ranges over the <= 3d fixed points of a bounded map — so the hinge count = sum over the <=3d fixed points x of #{r in window: S_r = x} = a sum of FIBER counts of the S-sequence over an EXPLICIT ALGEBRAIC set of <= 3d values. Fiber counts over arbitrary values are bounded by the window multiplicity 4H^{2/3}; the hinge count <= 3d * 4H^{2/3} recovers the pointwise bound. To win you need average-over-d: the fixed-point sets {x: g_d-fixed} for different d <= D_0 are fibers of DIFFERENT bounded maps — their union is small; the question becomes: does the S-sequence equidistribute on the union of fixed-point sets better than worst case? The union over d <= D_0 has size <= (3/2)D_0^2 = axis-strip budget... push this exact computation to its end — it may reproduce known bounds, but check whether the ALGEBRAIC structure of the fixed-point sets (they are the roots of N_d, whose strip localization in C is now a THEOREM) gives an extra handle mod p that generic values lack.
C. Your own idea. You have the complete map. If you see a combination we have not tried, take it — but check it against the graveyard first.

## Rules

- Rigor bar: every new lemma stated exactly, proof written out; every new identity machine-verified (write and run a sympy/python check in this directory; deliverables below).
- No effort cap. Do not stop at "this is hard"; stop only at a proof, a counterexample to a win condition, or a PRECISE new residual strictly weaker than all current wall faces.
- If you reach a partial result (e.g., W3 on a restricted L-range, or conditional on an effective h_0(p) certificate theorem), that is a valid deliverable — state it exactly.
- Machine checks: verify any new inequality empirically at p in {1009, 3001, 10007} using the orbit code in CRON_b1_crosscorr.py (function orbit(p)).

## Deliverables (this directory)

- CODEX_BREAKWALL_report.md — verdict first line in caps (WALL BROKEN via Wk / PARTIAL: ... / NOT BROKEN with new residual named); then the mathematics.
- Any verification scripts as CRON_breakwall_*.py with outputs shown in the report.
