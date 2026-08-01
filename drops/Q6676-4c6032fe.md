ANSWER Q6676 4c6032fe

# P3.2 theorem inventory — campaign additions after item 21

Origin labels:

- **HERE**: produced in the current life/codex campaign lane.
- **IMPORTED**: imported from the parallel cron session.
- **JOINT**: produced in one lane and independently audited, repaired, or merged by the other.

Status labels are used literally: **PROVED-all-h**, **VERIFIED-N**, **AUDITED-CONFIRMED**, **CONDITIONAL(...)**, and **EMPIRICAL**. The rows below are append-only after inventory item 21.

| No. | Name / tag | One-line statement | Status | Origin | Source pointer |
|---:|---|---|---|---|---|
| 22 | **[GAP-CASORATIAN] + centered-count inequality** | On every nonwrapping window, `Delta_(r,h)=N_h(r)/prod_(j=1)^h(r+j)^3`; hence `sum_(h<=H) C_h <= S/p + sqrt(N_coinc(H)-S^2/p)`. | **AUDITED-CONFIRMED** | **HERE** | `FABLE_SECTION_coinc_target.tex`; Q6550 |
| 23 | **[COINC_beta] => energy exponent** | If `N_coinc(p^beta) <= S^2/p + K S` for `1/2<beta<1`, then `F_p <= (1+o(1))p^(2-beta)`; at `beta=1`, `F_p <= (3+2sqrt(K))p`. | **CONDITIONAL(COINC_beta)** | **HERE** | `FABLE_SECTION_coinc_target.tex`; Q6550 |
| 24 | **[PER-H-WEIL-4H-1]** | For every fixed `h` and every nonzero additive frequency, the complete gap-Casoratian sum satisfies `|S_h(t)| <= (4h-1)sqrt(p)` at good primes. | **PROVED-all-h** | **HERE** | Q6587 |
| 25 | **[COINC-ADDITION-LAW]** | For all `h,d>=1`, `N_(h+d)(x)=N_d(x+h)N_(h+1)(x)-(x+h+1)^6N_(d-1)(x+h+1)N_h(x)`, equivalently the displayed bilinear update for `delta_(h+d)`. | **AUDITED-CONFIRMED** | **HERE** | `FABLE_SECTION_coinc_target.tex`; `CODEX_COINC_SYMBOLIC_report.md`; Q6609 / q20 audit |
| 26 | **[HONEST-TRIPLE-POLES]** | `N_h(-j)=(-1)^(j-1)((j-1)!(h-j)!)^3 b_(j-1)b_(h-j) != 0`; thus every finite pole of `delta_h` has exact order 3 and infinity is an exact order-3 zero. | **PROVED-all-h** | **JOINT** | `CODEX_IRRED_THEOREM_report.md` Sec. 2; q20 audit |
| 27 | **[GENUS-hk]** | For `h!=k`, under the stated simple-critical-value and branch-disjointness hypotheses, the normalization of `delta_h(x)=delta_k(y)` has `g=6hk-3h-3k-2`. | **CONDITIONAL(simple critical points + branch disjointness + irreducibility)** | **HERE** | `CODEX_COINC_SYMBOLIC_report.md` Sec. S3; q20 audit |
| 28 | **[MON-COMPONENT-DICHOTOMY]** | Under full monodromy, `X_(h,k)` is irreducible for `h!=k`; for `h=k>=2` it has exactly the diagonal and one irreducible off-diagonal component. | **CONDITIONAL(MON_h and MON_k)** | **HERE** | `CODEX_IRRED_THEOREM_report.md`; q20 audit |
| 29 | **[CRIT(h) => S_(3h)]** | If `s_h>=2h-1` simple nonzero critical values, then `delta_h` is indecomposable and `GeomMon(delta_h)=S_(3h)`; the same hypothesis gives the required same-gap and cross-gap irreducibility consequences. | **CONDITIONAL(CRIT(h))** | **HERE** | `CODEX_IRRED_THEOREM_report.md` |
| 30 | **[P16-COINC]** | Under `CENSUS_H` and good reduction, for `H<=p^(1/6)/4`, `N_coinc(H)<=S^2/p+4S`; with additive completion the honest range is `H<<p^(1/6)/(2+log p)^(2/3)`. | **CONDITIONAL(CENSUS_H + good reduction / completion hypotheses)** | **HERE** | `campaign3_questions/answers/Q6607_p16_theorem.tex` |
| 31 | **[SP-SELF-DUAL-G_h]** | Reflection gives a symplectic self-duality of the Fourier object `G_h`, hence `G_geom(G_h) subseteq Sp_(4h-2)` for every `h`. | **PROVED-all-h** | **HERE** | Q6612 |
| 32 | **[THM-Q-SIGMA-ONE-HALF]** | With `H=min(Delta,p-3)`, `Q_p(Delta) <= (27/4)H^(5/2)`; the proof also supplies the reusable master box inequality. | **AUDITED-CONFIRMED** | **IMPORTED** | `FABLE_SECTION_sigma_half.tex`; Q6577 |
| 33 | **[THM-ABSTRACT-Q-H2LOGH]** | For any finite word with `R_d<=Cd`, `Q_H<=22C H^2(1+log H)`; a multiscale AP construction shows the `H^2 log H` order is abstractly sharp. | **AUDITED-CONFIRMED** | **IMPORTED** | `CRON_SECTION_h2logh.tex`; Q6592; Q6620 |
| 34 | **[THM-NEAR-WALL-3/7]** | For `1<=D<=H`, `K_p(H,D)<<H D^(11/7)`; the one-based W-injection is `sum_u binom(W_H(u),2)<=Q_p(2H)`. | **AUDITED-CONFIRMED** | **JOINT** | Q6604; Q6630; Q6657; `drops/Q6661-1636c719.md` |
| 35 | **[ZWIN-RESTART]** | The windowed restart count satisfies `Z(H)<=1+(3/2)H^(2/3)`, and the associated bad-pole/critical contribution is `Bad(H)<<H^(5/3)`. | **PROVED-all-h** | **IMPORTED** | Q6632 |
| 36 | **[TRIANGLE-DICTIONARY]** | In the bounded-gap collision graph, `S_1(D)` is exactly the edge count and `Q_D` exactly the triangle count; therefore `Q_D<=(sqrt(2)/3)S_1(D)^(3/2)`, with asymptotically sharp constant. | **AUDITED-CONFIRMED** | **JOINT** | Q6619; Q6651 |
| 37 | **[N-SQUAREFREE-ALL]** | `N_h` is squarefree over `Q` for every `h>=2`. | **PROVED-all-h** | **HERE** | `CODEX_CRIT2H_report.md` |
| 38 | **[CRIT-RESULTANT-C0-ALL]** | The critical-value resultant `C_h^crit(T)` satisfies `C_h^crit(0)!=0` for every `h>=2`, with the exact constant-term factorization through `Disc(N_h)`. | **PROVED-all-h** | **HERE** | `CODEX_CRIT2H_report.md` |
| 39 | **[CRIT-CERT-32]** | Exact certificates prove `C_h^crit` squarefree and nonzero at zero for `2<=h<=32`, and `gcd(C_h^crit,C_k^crit)=1` for all 465 pairs; hence the Morse/monodromy/component tower is unconditional through 32. | **VERIFIED-32** | **HERE** | `CODEX_CRIT2H_report.md`; `CODEX_IRRED_THEOREM_report.md` |
| 40 | **[NONCOLLAPSE-ALL]** | For every `h>=2`, once the critical points are simple, the reflected squared critical values are not all equal; for `h>=5` this follows from fiber degree, and `h=2,3,4` from exact coefficient checks. | **CONDITIONAL(C_h^crit squarefree / Morse gate)** | **IMPORTED** | Q6659; `CRON_FRESH_EYES_pointwise.md` AT.23 |

## Dead-route certificates added to the inventory

| No. | Dead route | One-line death certificate | Status | Origin | Source pointer |
|---:|---|---|---|---|---|
| 41 | **[DEAD-LOW-GENUS-COLLAPSE]** | The distinct-gap normalization has exact genus `6hk-3h-3k-2=Theta(hk)`, not `O(h+k)`; pairwise Weil therefore stops at the `p^(1/6)` scale. | **AUDITED-CONFIRMED (DEAD)** | **HERE** | `CODEX_COINC_SYMBOLIC_report.md` Sec. S3; Q6550 / q10-lowgenus |
| 42 | **[DEAD-H-ALGEBRAIZATION]** | The growing-gap family cannot be replaced by one bounded-degree algebraic family in the `h` variable: degree, transfer length, rank, and conductor grow with `h`, so fixed-family large-sieve technology does not apply. | **AUDITED-CONFIRMED (DEAD)** | **HERE** | Q6546; q14-DLS ledger |
| 43 | **[DEAD-SHEAF-TRIANGLE-INDUCTION]** | The addition law is a shifted bilinear identity at one base point, not a convolution closing the Fourier sheaves; distinguished-triangle induction retains the growing conductor and generic `H^4 p^(3/2)` budget. | **AUDITED-CONFIRMED (DEAD)** | **HERE** | Q6547; Q6565 / q22 audit |
| 44 | **[DEAD-vdC-TRANSFER-LENGTH]** | van der Corput differencing replaces one transfer word by coupled longer/locked-clock words; transfer length and conductor do not contract, so no net exponent gain results. | **AUDITED-CONFIRMED (DEAD)** | **HERE** | Q6580; `FABLE_SECTION_coinc_target.tex` route audit |
| 45 | **[DEAD-MASTER-VARIETY]** | A universal master variety packages the indices only by letting ambient rank/dimension grow with `H`; the resulting point-count error is at least the original family error. | **AUDITED-CONFIRMED (DEAD)** | **HERE** | Q6546 |
| 46 | **[DEAD-FIXED-INDEX-MOMENTS]** | At every fixed moment order, summing the index choices loses a factor `H^(q-1)`; even perfect PAIR-FLAT input cannot close COINC, so only genuinely horizontal signed cancellation survives. | **AUDITED-CONFIRMED (DEAD)** | **HERE** | Q6640 |
| 47 | **[DEAD-ROWWISE-SCHUR-L1]** | Taking absolute values row by row destroys the Gram signs; the best unconditional Schur/L1 estimate does not even recover the `3/2` energy exponent. | **AUDITED-CONFIRMED (DEAD)** | **HERE** | Q6602 |
| 48 | **[DEAD-COMBINATORIAL-DELOG]** | A multiscale arithmetic-progression word has `R_d<=d` but `Q_H=Theta(H^2 log H)`; therefore no purely combinatorial argument can remove the logarithm. | **AUDITED-CONFIRMED (DEAD)** | **IMPORTED** | Q6619; Q6592; Q6620 |
| 49 | **[DEAD-BOTTOM-SCALE-PROPAGATION]** | The window `2/3` multiplicity law is only a bottom-scale input and does not bootstrap by itself to mesoscopic `S_D,Q_D<<N`; the purported new theorem was a rediscovery of the earlier sigma-one-third lemma. | **AUDITED-CONFIRMED (DEAD)** | **JOINT** | Q6646 and Sec. 121 erratum; Q6521; q29 last-wall follow-up |

## Inventory hygiene notes

- The old **VERIFIED-16** and **VERIFIED-30** ranges are superseded by item 39 (**VERIFIED-32**) and should not receive separate new numbers.
- The two window `2/3` statements in ledger Sec. 121 are **not** added: the Sec. 121 erratum identifies them as rediscoveries of the pre-existing sigma-one-third/window-multiplicity theorem. Q6646 remains useful only as a second restart/companion proof.
- Item 39 does **not** prove `[CRIT-2H]` for unbounded `h`; it proves it exactly through `h=32`. The all-h residual remains quotient primitivity / `W_h` irreducibility.
- Item 40 removes noncollapse as an independent obstruction but remains conditional on the simple-critical-point gate; it does not by itself prove all-h Morse.
- No purely empirical item is promoted in this append. The COINC numerics and `M_4 approximately 3` data remain in the empirical dossier rather than the theorem table.
