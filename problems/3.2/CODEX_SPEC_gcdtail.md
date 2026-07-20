# CODEX SPEC W4: split affine gcd-tail lemma — E_p(H) << H^{2-delta}

## Prize

The W3 harvest established: Z(p) << p/H + p^{2/3} H^{(beta-2)/3} for
E_p(H) << H^beta uniformly to H = sqrt(p). Current: beta = 8/3 (E << H^{8/3}).
ANY beta < 2 beats p^{2/3}:
- beta = 3/2 -> Z << p^{7/12}
- beta -> 2^- -> Z << p^{2/3-epsilon}
TARGET: prove E_p(H) << H^{2-delta} for some delta > 0. Even delta = 0.01 is a
breakthrough.

## Read first

1. problems/3.2/energy_result.tex (all of it — the G3 shallow/diagonal structure,
   the G4 stall report, exact identity E_p^o = sum C(k_x, 2))
2. problems/3.2/meso_result.tex (root strips, center recurrence, diagonal square law)
3. problems/3.2/nv_theorem.tex (NV range theorem — Δ_{h,k} ≢ 0 for h<k<p)
4. In proof.tex: prop:column (unpolluted restatement), rem:collision, cor:sep-block
5. problems/3.2/DOCTRINE.md — search for "gcd-tail" and "W3 harvest"

## The exact remaining lemma (from W3 harvest)

The "split affine gcd-tail lemma": for pairs (d,r) with min(d,r) > H^{3/4}
(the "deep pairs" not covered by W3's G3), show the energy contribution is
o(H^2). Equivalently:
  #{(d,r,x) : d+r <= H, min(d,r) > H^{3/4}, x in F_p affine, N_d(x) = N_r(x+d) = 0} = o(H^2)

## Facts you can use

- NV range theorem: Δ_{h,k} ≢ 0 in F_p[x] for all 1<=h<k<p. The bordered
  certificate N_h(x)B_k(x) - N_k(x)B_h(x) = -Pi_h^2 N_{k-h}(x+h) is NONZERO.
  For pairs (d,r) with r>d, the witness is N_d at the collision x.
- Root strip: all complex roots of N_d lie in -d < Re(z) < -1.
- Per-pair: at most 3(min(d,r)-1) affine collisions (Bezout, saturated).
- Diagonal square law: R_{d,d} = l_d |D_d| Q_d^2, so diagonal energy is controlled.
- Column structure (prop:column): unpolluted columns have O(H^{2/3}) zeros.
  Polluted columns (p | b_{m-1}) can have unbounded zeros but are O(M^{2/3}) in number.
- Shallow strips (min(d,r) <= H^{1/4}) and center witnesses: already O(H^{5/3}).
- E_p^o(H) = sum_x C(k_x, 2) where k_x = #{(d,r): d+r<=H, N_d(x)=N_r(x+d)=0}.
- Max polluted k = 50 (p=3331), max E = 16 (p=653) for H ~ sqrt(p), p <= 5000.
  Empirically E/H^{3/2} < 0.16.

## Approaches to try (ranked)

1. **Column-pair partition + column bound:** Deep-pair collisions at x group by
   column. In unpolluted columns, k_x <= column^{2/3} bound from prop:column applied
   to the (d,r)-fiber. The column-pair energy E = sum C(k,2) is then sum over columns
   of O(column_height^2). If this can be bounded by o(H^2) using the column structure,
   we're done.

2. **Resultant divisibility chain:** For fixed d, the resultant R_{d,r} is divisible
   by center lattice elements. Chain the divisibility along r to get cancellation.
   If sum_{r} 1_{p | R_{d,r}} = o(d) (summed over r at fixed d), then E << H^{2-delta}.

3. **GCD of resultants:** R_{d,r} and R_{d,r'} share factors only through N_d;
   gcd(R_{d,r}, R_{d,r'}) = Res(N_d, gcd(N_r(x+d), N_{r'}(x+d))) — a lower-degree
   resultant. The r-fiber of collisions at fixed d,x is bounded by degx(N_d) = 3(d-1).
   Use this to get a column-energy bound.

4. **Averaging over d (Plancherel):** Replace worst-case by average:
   (1/H) sum_d E_p(d, H-d). Each d-slice is sum_r 1_{collision}. A double-counting
   with Cauchy-Schwarz may give H^{2-delta} on average, then pass to all-d.

## Deliverables

Produce `gcdtail_result.tex` with:
- Whatever you prove (partial or full), stated precisely with proof
- If you prove E << H^{2-delta}: state it as a proposition, give delta
- If you can't: HONEST stall report with the exact obstruction
- A verification script `gcdtail_verify.py` checking key identities for p <= 2000

## Verification command

```
cd ~/repos/Ramanujan_Challenge/problems/3.2 && python3 gcdtail_verify.py
```
All checks must PASS. No sorry, no unverified claims.

## Stall protocol

If stuck after exhausting all 4 approaches: deliver what you proved (even partial
structural lemmas) + the exact obstruction as a .tex file. Partial results with
clear next steps are valuable. Do not fabricate.
