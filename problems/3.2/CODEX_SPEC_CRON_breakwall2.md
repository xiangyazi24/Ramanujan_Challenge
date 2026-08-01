# CODEX SPEC: BREAK THE WALL, STRIKE 2 (single deep strike, no effort cap)

## Mission

Same prize as CODEX_SPEC_CRON_breakwall.md (read it first — objects, win conditions W1–W7, arsenal, graveyard all still valid), but the map has moved. Strike 1 (your predecessor, report CODEX_BREAKWALL_report.md) proved [NO-RIGHT-2-3] (no right composition factor of degree 2/3, all h) and isolated [ZERO-TAIL-2]. Since then the ledger (CRON_FRESH_EYES_pointwise.md appendix AT.33–AT.39) has finalized the wall's full spectrum. Your mission: break the wall at its NEW weakest face, or deliver a strictly weaker named residual than the current bottom.

## The wall ladder (current, final form — attack the bottom)

[BDH-LAG] ≤ [SAME-LAG-L2] ≤ [ZERO-TAIL-2] ≤ [MESO-TOTAL] = W1.

- [SAME-LAG-L2]: Σ_{h≤D} R_h² ≪ N at D=√N·L ⟹ W1 (Cauchy). Empirically ≈ 4.3D.
- [ZERO-TAIL-2]: 𝒵_D = max_{t>T₀} t²·#{h≤D: R_h≥t} ≪ N, T₀=⌈√(N/D)⌉ ⟹ W1 (layer-cake, banked).
- [BDH-LAG]: Σ_{h≤D}(R_h−κ_h)² ≪ √(DN) (κ = mirror skeleton) ⟹ [SAME-LAG-L2].
- NEW pointwise sufficient face (elementary, unexploited): [PT-ANTICONC]: R_h ≤ N^{1/4−δ} uniformly for h ≤ √N·L ⟹ [SAME-LAG-L2] via Σ ≤ D·max². Empirically max R_h = 8 CONSTANT across p ≤ 30011. The degree bound gives only 3h−3 ~ 3√N·L.

## Killed since strike 1 (do not re-till; death certificates in ledger AT.35–AT.39)

- Per-row L2 census route: dead (your own Theorem 5.1 + two more independent countermodels: reflection-symmetric cluster unions with Q_D=0, family lower bound cDN^{2/3}). Row marginals + window sparsity + energy CANNOT give W1. The content is arithmetic anti-concentration of the lag profile.
- holonomic Stepanov on single N_h: dead (no degree surplus; returns are codimension-one events).
- Abstract histogram interpolation at any certificate growth rate h₀(p): dead.

## Most promising directions (our ranking; you may override with reasons)

A. **Slice-disjointness double count toward [ZERO-TAIL-2]**: heavy lags h₁≠h₂ with R_{h_i} ≥ t share base roots only where p | Res(N_{h₁}, N_{h₂}) (adjacent pairs are coprime over Q — banked; establish/exploit what is true for general pairs: compute gcd(N_{h₁},N_{h₂}) over Q for small pairs FIRST — sympy — then use resultant size exp(O(h²log h)) + p-divisibility counting to bound #pairwise-intersections; derive t·A_D(t) − pair corrections ≤ S_D-type inequalities and close [ZERO-TAIL-2] conditionally on a PRECISE cross-resultant non-vanishing statement, then attack that statement).
B. **[PT-ANTICONC] via family propagation**: many roots of one N_h + renewal identity ⟹ forced structure across neighboring gaps (the three-term recurrence in the gap variable links N_{h−1}, N_h, N_{h+1} at a common root r: if N_h(r)=0 then the recurrence degenerates — write it out; a root of N_h makes the (r-shifted) continuant sequence hit 0 at position h, and continuant sequences with a zero restart with a known factorization: use the banked renewal/restart identities to show R_h ≥ K forces K disjoint restart events whose composition produces MORE collisions at controlled lags — pump K upward until contradiction with the degree bound at a longer lag, i.e., a self-improving/bootstrap inequality on the profile {R_h}).
C. **Mirror skeleton + stripped variance**: make κ_h exact (which collisions are FORCED by the reflection law/central vanishing/pole structure), prove the skeleton part, restate [BDH-LAG] as pure fluctuation, then attack the fluctuation with A/B.
D. Your own idea. Check the graveyard (both spec files + AT.35–AT.39) first.

## Rules

- Rigor bar unchanged: every new lemma stated exactly with complete proof; every new identity machine-verified (sympy in this directory; p ∈ {1009, 3001, 10007}; orbit(p) in CRON_b1_crosscorr.py).
- No effort cap. Stop only at: a proof of a win condition, a counterexample, or a PRECISE new residual strictly weaker than [BDH-LAG]/[ZERO-TAIL-2].
- Partial results are valid deliverables if stated exactly (e.g., [ZERO-TAIL-2] conditional on a named cross-resultant statement + that statement proved for a positive-density subrange).

## Deliverables (this directory)

- CODEX_BREAKWALL2_report.md — first line verdict in caps (WALL BROKEN via ... / PARTIAL: ... / NOT BROKEN, new residual [NAME]); then the mathematics, then machine-verification transcript.
- Verification scripts as CRON_breakwall2_*.py with outputs shown in the report.
