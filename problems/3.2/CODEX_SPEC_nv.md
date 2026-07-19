# CODEX SPEC: the range nonvanishing theorem (NV_{p,H})

## Target

Prove: for every prime p >= 7 and all 1 <= h < k <= H, the reduced bordered
certificate Delta_{h,k} is NOT the zero polynomial in F_p[x] — for H as large
as you can reach. Ranked goals:
  (G1) H ~ c p^{1/4} (closes the unconditional uniform fiber bound O(p^{3/4})),
  (G2) any H = p^{delta} with fixed delta > 0,
  (G3) any H growing faster than log p (beats the size-forcing barrier).
Even (G3) is a real advance. If nothing beyond log p is provable, deliver a
sharp stall report identifying the exact obstruction.

## Established infrastructure (all in problems/3.2/fiber_bound.tex — READ IT FIRST)

- N_0=0, N_1=1, N_{m+1}(x) = P(x+m)N_m(x) - (x+m)^6 N_{m-1}(x),
  P(n) = 34n^3+51n^2+27n+5, deg N_m = 3(m-1), leading coeff l_m:
  l_{m+1} = 34 l_m - l_{m-1}.
- B_m(x) = -(x+1)^6 N_{m-1}(x+1) (second solution).
- Pi_m(x) = prod_{j=1}^m (x+j)^3.
- Casoratian: N_h B_k - N_k B_h = -Pi_h^2 N_{k-h}(x+h).
- D_{h,k} = N_h(Pi_k - B_k) - N_k(Pi_h - B_h) = Pi_h * Delta_{h,k}, with
  Delta_{h,k} = N_h Pi_{k-h}(x+h) - N_k + Pi_h N_{k-h}(x+h),
  deg Delta_{h,k} = 3(k-1), leading coeff l_h + l_{k-h} - l_k (< 0 over Z).
- Endpoint evaluations (VERIFIED, in fiber_bound.tex):
  for 1 <= r <= h:
    Delta_{h,k}(-r) = (-1)^{r-1} b_{r-1}(b_{h-r} - b_{k-r}) ((r-1)!(k-r)!)^3
  for h < r <= k:
    Delta_{h,k}(-r) = (-1)^{r-1} b_{k-r}(b_{r-h-1} - b_{r-1}) ((r-1)!(k-r)!)^3
- Reflection/midpoint forced factors (exactly one, or all three if h,k both
  even): (2x+h+1) if 2|h, (2x+k+1) if 2|k, (2x+h+k+1) if 2|(k-h).
- No-consecutive-zeros: b_j, b_{j+1} never both 0 mod p (j <= p-2).
- N_m nonvanishing mod p for 2 <= m <= p (lem:nonvanish in proof.tex).
- Z(p) = #zeros of b in [0,p) = O(p^{2/3}) unconditional.
- KNOWN OBSTRUCTION (do not rediscover): at p=131, (h,k)=(12,55):
  b_10 = b_53 = 15 and b_11 = b_54 = 15 mod 131, so Delta(-1) = Delta(-2) = 0
  — the first two evaluations alone cannot work. Delta_{12,55}(-3) = 2 != 0
  there, so the FULL evaluation network still detects nonvanishing.
- PROVED (trivial route): for k <= log p/log 34, Delta(-1) != 0 by integer
  size forcing (b_n < 34^n). This is the baseline to beat.
- Empirical: ZERO identically-vanishing Delta for all p <= 700, h < k <= 40.

## Suggested attack: the full singular-specialization network

Suppose Delta_{h,k} == 0 in F_p[x]. Then ALL k evaluations vanish. Using
no-consecutive-zeros on the prefactors:
- for r in [1,h]: b_{r-1} = 0 OR b_{h-r} = b_{k-r} (mod p);
- for r in (h,k]: b_{k-r} = 0 OR b_{r-h-1} = b_{r-1} (mod p).
Since the zero set has no two consecutive elements, at least every other r
in [1,h] forces a lag-d coincidence b_j = b_{j+d} with d = k-h, and at least
every other r in (h,k] forces a lag-h coincidence. So one gets:
  >= floor(h/2) coincidences at lag (k-h) among indices j in [0, k-1], AND
  >= floor((k-h)/2) coincidences at lag h among indices j in [0, k-1].
Each coincidence b_j = b_{j+d} (call the common value a_j) pins the orbit
point [b_{j+1} : b_j] in P^1(F_p) to the specific position determined by the
return row: N_d(j)(j+1)^3 b_{j+1} = (Pi_d(j) - B_d(j)) b_j.
Ideas to convert this rigidity into a contradiction:
1. TRIPLE COINCIDENCES: two coincidences at the same lag d and nearby bases
   j1 < j2 with j2 + d <= p - 1... if any THREE indices share one value, the
   already-proved certificate machinery applies recursively: the value-triple
   (j, j', j'') in one fiber forces Delta_{j'-j, j''-j}(j) = 0 — a ROOT, not
   identical vanishing; count roots vs degrees to bound how many coincidences
   are consistent, and compare against the >= h/2 forced ones. Optimize over
   (h, k) to find pairs where the count is contradictory.
2. INDUCTION ON k: relate Delta_{h,k} == 0 to vanishing statements for
   SMALLER parameters via the continuant addition/Bezout identities in
   proof.tex (prop:bezout, rem:dodgson, lem:restart). E.g. examine whether
   Delta_{h,k} == 0 forces Delta_{h,k-h}-type degenerations after the
   restart identity, descending until the size-forced range k <= log_34 p
   is reached — where nonvanishing is already proved.
3. SECOND COEFFICIENT: the leading coefficient of Delta is l_h + l_{k-h} - l_k;
   if Delta == 0 mod p then p | (l_h + l_{k-h} - l_k). Compute the
   next-to-next coefficient (after centering at the reflection midpoint) in
   closed form and show gcd(leading, that coefficient) has only small prime
   factors. (This was attempted before and stalled — only pursue if you see
   a genuinely new closed form; document what you compute.)
4. Any other route you devise. You have full freedom.

## Verification duty

Whatever partial theorem you prove, add a pure-Python check to a new file
problems/3.2/nv_verify.py (same style as fiber_verify.py: PASS/FAIL lines,
nonzero exit on failure) that tests the specific lemmas/identities you use
on random primes/parameters. Run it and include the output tail.

## Deliverables

1. problems/3.2/nv_theorem.tex — LaTeX fragment: theorem(s) + full proof,
   notation matching fiber_bound.tex/proof.tex. If only a partial range is
   proved, state exactly that. STALL REPORT as leading comment if goals
   G1-G3 all fail.
2. problems/3.2/nv_verify.py — the verification script.
3. Do not modify proof.tex or fiber_bound.tex.

## Constraints

- No numerics-as-proof. Every step of the LaTeX proof self-contained given
  proof.tex + fiber_bound.tex infrastructure.
- Cross-check every identity you introduce symbolically before using it.
- If you find a COUNTEREXAMPLE to (NV_{p,H}) at any range (an identically
  vanishing Delta), that is a headline result — verify it independently and
  report it prominently.
