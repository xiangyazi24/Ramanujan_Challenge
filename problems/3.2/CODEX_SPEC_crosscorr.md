# CODEX SPEC — cross-gap correlation protocol (from Q6546 §measurement)

Report to `CODEX_CROSSCORR_report.md`. Exact integer F_p arithmetic. Definitions in
`campaign3_questions/CTX.txt` + delta_h(r) = (b_r c_{r+h} - b_{r+h} c_r) mod p.

For p in {997, 1999, 4001} and H in {ceil(p^0.4), ceil(p^0.5), ceil(p^0.6)}:
1. For all 1<=h<k<=H compute E_{h,k} = #{(r,r') in windows: delta_h(r)=delta_k(r')} - |I_h||I_k|/p
   (exact rational baseline; |I_h|=M-h). This is O(H^2) pairs; for each h precompute the
   value-count vector of delta_h (Counter over F_p, O(p) each, total O(pH)) then
   E_{h,k} = sum_a n_h(a) n_k(a) - |I_h||I_k|/p via sparse dot products. If O(H^2) dot
   products too slow at H=p^0.6, restrict k-h <= 200 and say so.
2. A_d(H) = sum_{h} E_{h,h+d} for d=1..min(H-1,200). Report: max_d |A_d|/p, V(H) = sum_d A_d^2/(p^2 H),
   the top-5 d by |A_d| (candidates for structured correlations), and A_d profile deciles.
3. Same-gap: E_{h,h} (off-diagonal r!=r') separately; verify ~ +p from the H_h component
   (component census): report (E_{h,h} - p)/  (h sqrt p) deciles over h.
4. Verdict lines: (i) is cross-gap correlation mass sum_{h<k}|E_{h,k}| << pH empirically?
   (ii) any structured d (e.g. d even, d=h, small d) carrying anomalous mass?
   (iii) does the h=1 row (cube-root exception) show its predicted p=1 mod 3 dependence?
Pure python3 stdlib. All five verdicts. No early stop.
