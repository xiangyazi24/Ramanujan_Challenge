# CODEX SPEC — direct power-iteration measurement of the campaign's final object

Report to `CODEX_LAMBDAMAX_report.md`. Definitions in campaign3_questions/CTX.txt.
Goal: measure lambda_max of the centered Gram matrix — the direct empirical test of
[OP-OFF-0] (any absolute constant => COINC => F << p^{4/3}).

For p in {997, 1999, 4001} and H = ceil(p^{2/3}):
1. Compute the value histograms n_h(a) for all h <= H (O(pH) via b,c recurrences and the
   Delta identity). Form the centered vectors q_h(a) = n_h(a) - |I_h|/p (real, length p).
2. The Gram matrix Gamma_{h,k} = p * <q_h, q_k> (this equals sum_{t!=0} S_h conj(S_k)
   by Parseval — verify the identity numerically on one small case p=199).
3. WITHOUT forming Gamma: power iteration for the top eigenvalue of the PSD matrix
   Gamma: matvec (Gamma v)_h = p * <q_h, sum_k v_k q_k> — two passes over the (h,a)
   table per iteration, O(pH). Run ~40 iterations from a random start + deflate to get
   lambda_1 and lambda_2. Also compute the diagonal max D_max = max_h Gamma_{h,h}, the
   all-ones quadratic form 1^T Gamma 1 = p * ||sum_h q_h||^2 (= p*N_coinc - S^2 check),
   and Tr(Gamma).
4. Report per cell: lambda_1/p^2, lambda_2/p^2, D_max/p^2, (1^T Gamma 1)/(p^2 H),
   Tr/p^2 (~ H expected). Peel the h=1 row (cube exception) and report both with/without.
5. Also the off-diagonal part: O = Gamma - diag; estimate ||O||_op by power iteration on
   O (signs matter, use symmetric iteration); report ||O||_op/p^2 and ||O||_op/(sqrt(H) p^{3/2})
   (random-matrix normalization — expect O(1) if the Wishart picture holds).
6. Verdict: is lambda_1/p^2 bounded (the [OP-OFF-0] scenario) across p? Is ||O|| at the
   random sqrt(H)p^{3/2} scale? Exact integers where feasible, floats for eigenvalues.
Pure python3 stdlib (no numpy needed but allowed if present). No early stop.
