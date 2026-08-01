# CODEX SPEC: [WALL] falsification probe — box-level L2 statistics, operator norms, translation scan, gauge scan

## Mission

Implement the experimental protocol of `chatgpt-answers/Q6573.md` §18 (READ IT FIRST, plus §17 for the counterexample worlds this is designed to detect and §21 for the gauge caveat). Prior experiments tested global means/prefix windows; this probe targets exactly the blind directions: translated windows, weighted/spectral modes, all-level (nonzero-value) statistics, and normalization dependence.

## Objects

Orbit vectors v_n=(b_n,c_n) mod p (reuse verified orbit code from `CRON_radon_spectrum.py` / `CRON_b1_crosscorr.py`; cross-check one prime). Delta_p(r,h) = b_r c_{r+h} - c_r b_{r+h} mod p on the nonwrapping triangle.

Primes: {1009, 3001, 10007} (add 30011 if runtime allows).

## Task 1 — translated-box zero-fiber statistics (§18.1, §18.5)

For H in {p^{1/3}, p^{1/2}, p^{2/3}} (rounded) and a sliding set of translated gap-windows I=(H0, H0+H] covering [1, p^{2/3}·4] with ~8-16 translates per scale, and J = full admissible r-range (also 8 translated r-windows of length p/4):
- V_0(I,J) = sum_{r in J} |d_I(r) - q_I(r)/p|^2, report V_0/|I| per box (WALL-zero predicts p^{o(1)}).
- Record the max over boxes, and the argmax box.

## Task 2 — column variance / local GPRV (§18.2)

V_col(I,J) = sum_{h in I} |C_h(J) - |J_h|/p|^2, normalized by |J|, same box grid.

## Task 3 — all-level COINC defect (§18.3) with GAUGE SCAN (§21)

For each box: V_all(I,J) = N_coinc(I,J) - |Omega(I,J)|^2/p, normalized by |Omega(I,J)|, where N_coinc = sum_a nu(a)^2, nu(a) = #{(r,h) in box: Delta_p(r,h)=a}.
CRITICAL: compute this in THREE normalizations of Delta:
  (g1) raw determinant Delta_p(r,h);
  (g2) Delta divided by an h-dependent scalar: Delta / lc_h where lc_h = leading coefficient of N_h mod p (compute lc via the banked recurrence lc(N_h)=U_h(34,1) — find/verify the recurrence in repo, else compute N_h mod p directly for h in the window and extract lc);
  (g3) Delta / (r-dependent trivial gauge): Delta / (b_r b_{r+h}) whenever both nonzero (projective chart gauge) — value in F_p.
Report V_all/|Omega| per gauge per box. A gauge whose all-level statistics are NOT flat while others are flat = major finding (Q6573 §21 predicts this possibility).

## Task 4 — singular-value probe (§18.4)

For boxes of manageable size (|I|,|J| <= 2000): form the centered zero-fiber matrix B_{r,h} = 1[Delta=0] - 1/p on the box, compute the largest singular value (scipy.sparse.linalg.svds on the dense/sparse matrix). Compare to the random benchmark: same-shape matrix with iid Bernoulli(mean density) entries, 3 replicas. Report ratio data/benchmark per box. Ratio growing with p or >3 = the weighted singular mode of §17.4.

## Task 5 — verdict table

Per §17's failure-mode list: for each mode (localized heavy row / thin resonance class / nonzero-level concentration / weighted singular vector), state what the data shows (detected / not detected at these scales) with the relevant numbers.

## Deliverables

- `CRON_wallprobe.py` (progress prints every <=10s).
- `CODEX_WALLPROBE_report.md` with the tables and verdicts.
- `wallprobe_results.json`.

## Rules

- Reuse existing orbit code; verify against one banked statistic (e.g. E^pi/p ~ 3) before the main run.
- numpy vectorization; the (r,h) triangle at p=10007 up to h<=p^{2/3}*4 is large — budget memory, process per h-window.
- No effort cap; complete all five tasks.
