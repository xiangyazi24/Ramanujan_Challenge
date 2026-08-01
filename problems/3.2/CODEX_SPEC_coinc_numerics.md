# CODEX SPEC — COINC second-moment numerics (campaign 3, avenue a1)

## Mission
Machine-verify the NEW second-moment reduction and measure the true structure of
N_coinc(H). Everything must be exact F_p arithmetic (python3, no floats for counts).
Write report to `CODEX_COINC_NUMERICS_report.md` in this directory (overwrite ok).
Background/definitions: read `campaign3_questions/CTX.txt` (authoritative; contains
(D1)-(D6)). Do not read other files unless needed.

## Tasks (in order; commit nothing; just compute and report)

### T1 — verify (D1) Delta identity
For p in {97, 199, 499}, all 1<=r, r+h<=p-2 with h<=30, nonwrapping, r+j != 0 mod p
for 1<=j<=h: check Delta_{r,h} := (b_r c_{r+h} - b_{r+h} c_r) mod p equals
N_h(r) * inverse(prod_{j=1}^h (r+j)^3) mod p, where N_h from the gap recurrence.
Report: pass/fail counts, and list every failing (p,r,h) with both values.
Also characterize the excluded/singular pairs: how many, what happens there.

### T2 — Parseval bookkeeping (D3) exact check
For p in {97, 199}, H in {5, 10, 20}: compute B_t(H) for all t, verify
sum_{t!=0} |B_t|^2 == p*N_coinc(H) - (#S_H)^2 exactly (integer identity;
compute N_coinc by hashing Delta values). Report exact equality or the discrepancy.

### T3 — N_coinc structure decomposition (THE KEY MEASUREMENT)
For p in {199, 499, 997, 1499} and H in {ceil(p^0.25), ceil(p^0.33), ceil(p^0.5), ceil(p^0.66)}:
compute N_coinc(H) and decompose into:
  (i)   diagonal pairs (r,h)=(r',h');
  (ii)  zero-zero pairs (Delta=0 both sides, (r,h)!=(r',h')) — this is (sum C_h)^2-ish, the circularity term;
  (iii) same-h off-diagonal (h=h', r!=r', Delta!=0);
  (iv)  reflection-forced pairs: identify empirically any exact symmetry (r,h)->(r*,h) with r* = p-1-r-h or similar giving forced equal Delta values — detect by testing candidate involutions;
  (v)   remainder ("generic").
Fit: N_coinc - (#S_H)^2/p against #S_H: report the empirical K(p,H) = (N_coinc - (#S_H)^2/p)/#S_H
per cell, and per component. VERDICT LINE: does [COINC] (bounded K) look TRUE, and
which component threatens it (esp. (ii) and (iv)).

### T4 — difference-curve anomaly scan
For p=997: for all pairs 1<=h<h'<=20, count P(h,h') = #{(r,r') in valid windows :
Delta_{r,h} = Delta_{r',h'}}. Expected ~ (window sizes)/p * ... ~ p for generic pairs.
Report the matrix of P(h,h')/p, flag every pair with P > 2.5p or P < 0.3p, and for
flagged pairs guess the algebraic cause (parity of h,h'? h'|h? h+h' even?). Same-h
column P(h,h) separately (subtract diagonal r=r').

### T5 — sum C_h profile
For p in {997, 1499, 1999}: compute C_h for h <= p^0.66, plot-in-text (deciles) of
cumulative sum vs H and vs H + sqrt(pH). Verdict: which curve tracks.

## Rules
- Pure python3 stdlib. Efficiency matters: b_n, c_n mod p via the recurrence in O(p);
  N_h(r) via Delta identity (after T1 validates it) to avoid O(p*H) polynomial evals;
  hash-count coincidences.
- If a task is infeasible at a stated size, SHRINK p (not the task) and say so.
- No estimates of remaining work, no stopping early: run all five tasks.
- Report = numbers + verdict lines. ASCII only.
