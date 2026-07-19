# CODEX SPEC W3: collision energy E_p(H) << H^{3/2} — the p^{7/12} target

## Prize

PROVED reduction (paper, eq:meso-triple-endpoint region): if
  E_p(H) = sum_{(d,r), d+r<=H} #{x in F_p : N_d(x) = N_r(x+d) = 0}  << H^{3/2}
uniformly for H <= sqrt p, then block triple-counting gives
  Z(p) << p^{7/12},
BEATING the current record Z(p) << p^{2/3}. This is now the sharpest
single-prime target in the campaign. (The full sqrt-p needs additionally the
low-fiber amplification R^d + L_p << H — NOT this spec's problem.)

## Read first

problems/3.2/meso_result.tex (all of it — root strips, centered norm,
center lattice, diagonal square law, projective validity, the two missing
lemmas), problems/3.2/meso_resultants.md (exact factorizations d+r<=8,
prime statistics p<=5000), problems/3.2/nv_theorem.tex, and in proof.tex:
rem:adj-res, prop:bezout, cor:sep-block (saturation!), rem:collision,
prop:column (unpolluted restatement + pollution characterization),
rem:content, lem:restart, lem:nonvanish.

## Facts you can use

- E_p(H) <= 3 W_p(H) with W_p(H) = weighted support
  sum (min(d,r)-1) 1_{p | R_{d,r}} — so W_p(H) << H^{3/2} suffices.
- p | R_{d,r} iff common projective root over the algebraic closure;
  infinity pairs (rho_p | d, r) are NOT affine collisions — E only counts
  affine F_p roots, so you may saturate away the infinity lattice.
- Diagonal square law R_{d,d} = l_d |D_d| Q_d^2.
- Even-block center lattice: T_b^{(a)} | R_{a,b} etc.
- Root strips: all roots in -d < Re < -1; R_{d,r} > 0.
- Adjacent: p | R(N_h, N_{h+1}) iff p | b_j, some j <= h-1.
- Empirical (meso_resultants.md): V/W/E statistics for p <= 5000; typical
  W_p(H) is O(1)-tiny; E ~ W.
- Per-pair root count: m_{d,r}(p) <= 3(min(d,r)-1) (saturated; +cut-edge
  allowance d-1 unsaturated).
- The x-witness of a pair (d,r) collision yields THREE orbit collisions
  x ~ x+d ~ x+d+r; conversely triples of zeros in one block give energy
  witnesses (the paper's block triple-counting).

## Ranked goals

(G1) W_p(H) << H^{3/2} (or E directly) uniformly for H <= sqrt p, all p.
(G2) Same with any exponent 3/2 < beta < 2: E_p(H) << H^beta gives
     Z(p) << max(p/H, p^{2/3} H^{(beta-1)/3}) — work out and state the
     optimized Z-exponent as a function of beta, and pin the beta range
     that still beats 2/3 (answer: any beta < 2 improves NOTHING unless
     the H-range and the p/H term balance — do this calculation exactly
     and honestly; if only beta = 3/2 helps, say so).
(G3) Prove W_p(H) << H^{3/2} for structured sub-families (e.g. all pairs
     with min(d,r) <= H^{1/2}, or the diagonal band |d-r| <= C) and
     identify exactly which sub-family carries the potential worst case.
(G4) Sharp obstruction: the minimal missing lemma, in the style of the
     prior G4 deliverables.

## Attack suggestions

1. Two supported pairs sharing the SAME witness x chain into longer orbit
   collision patterns: x ~ x+d ~ x+d+r and x ~ x+d' ~ x+d'+r' interleave
   into a multi-return configuration; the restart/Bezout calculus converts
   these into lower-parameter certificates. Formalize: the set of
   (d,r)-pairs supported at a COMMON x is controlled by the column
   structure of x (prop:column, unpolluted <= O(H^{2/3}) returns; polluted
   columns are O(H^{2/3}) many). Partition W_p(H) by witness column and
   count: W <= sum_x C(k_x + 1, 2)-ish — relate precisely, watch that a
   pair (d,r) needs returns at x at levels d and d+r i.e. TWO returns in
   the same column => W_p(H) <= sum_x C(k_x(H), 2) where k_x(H) = returns
   of column x within level window H. With unpolluted k_x <= c H^{2/3}:
   sum_x C(k_x,2) <= max k_x /2 * sum_x k_x = H^{2/3} R_p(H) <= H^{2/3} *
   (3/2)H^2?? — trivial R bound gives H^{8/3}: TOO BIG. The useful
   direction: R_p(H) <= H^2 always (row degrees), so sum C(k_x,2) <=
   H^{2/3}/2 * R... need R_p(H) << H^{5/6}?? circular. FIND THE RIGHT
   BOOTSTRAP: E and R bound each other both ways (paper has
   eq:meso-column-pairs: sum_{x unpolluted} C(k_x,2) <= E); a fixed-point/
   iteration between the two inequalities may self-improve. Analyze
   whether the coupled system (E <= f(R), R <= g(E)) has a nontrivial
   fixed point below the trivial bounds. THIS IS THE MOST PROMISING ROUTE:
   do the fixed-point algebra carefully.
2. Valuation route on the diagonal band: the square law forces even
   valuations off l_d D_d; p | R_{d,d} with p prime to l_d D_d implies
   p^2 | R_{d,d} — pairs are 'expensive' in height, so Hadamard counting
   of SUPPORTED diagonal pairs gains a factor 2: quantify.
3. Center/infinity lattices remove structured sub-families exactly.

## Deliverables

problems/3.2/energy_result.tex (+ energy_verify.py for any new identity,
PASS/FAIL), STALL REPORT convention. Computational duty: for p <= 5000,
compute the coupled quantities (R_p(H), E_p(H), sum_x C(k_x,2), polluted
counts) at H = floor(sqrt p) and test the fixed-point inequalities you
derive against data BEFORE claiming them. Do not modify existing files.
No numerics-as-proof.
