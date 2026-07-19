# CODEX SPEC: Uniform value-fiber bound for Apéry numbers mod p

## Target

Prove rigorously and write publication-ready LaTeX for:

**Proposition (uniform fiber bound).** Let b_n be the Apéry numbers,
b_0 = 1, b_1 = 5, with recurrence
(n+1)^3 b_{n+1} = P(n) b_n - n^3 b_{n-1},  P(n) = 34n^3 + 51n^2 + 27n + 5.
For a prime p >= 7 and a in F_p, let N_p(a) = #{ t in {0,...,p-1} : b_t ≡ a (mod p) }.
Then  max_{a in F_p^*} N_p(a) = O(p^{3/4})  with an explicit constant.

(The fiber a = 0 already has the better bound Z(p) = O(p^{2/3}) — proved
elsewhere in the paper; do NOT redo it. Your proposition is for a ≠ 0,
and the final statement can combine both into max_a N_p(a) = O(p^{3/4}).)

## Known infrastructure you may cite (all proved in problems/3.2/proof.tex)

1. **No consecutive zeros** (Lemma lem:no-consec): for p >= 5, b_j and b_{j+1}
   never both vanish mod p for 0 <= j <= p-2.
2. **Gap polynomials** N_h ∈ Z[x]: N_1 = 1, N_2(x) = P(x+1), and
   N_{h+1}(x) = P(x+h) N_h(x) - (x+h)^6 N_{h-1}(x).
   deg N_h = 3(h-1), and N_h is NOT the zero polynomial over F_p for
   2 <= h <= p, p >= 7 (Lemma lem:nonvanish).
3. **Endpoint evaluations** (Remark rem:content):
   N_h(-r) = (-1)^{r-1} b_{r-1} b_{h-r} ((r-1)!)^3 ((h-r)!)^3  for 1 <= r <= h.
   In particular N_h(-1) = ((h-1)!)^3 b_{h-1},  N_h(-2) = -5 ((h-2)!)^3 b_{h-2}.
4. **Reflection**: N_h(-m-h-1) = (-1)^{h-1} N_h(m).
5. **Zero-count bound**: Z(p) = N_p(0) <= (3^{4/3}/2) p^{2/3} + O(p^{1/3}).

## Proof skeleton to flesh out (verified sound by an independent audit;
## your job is to make every step airtight and catch any error)

### Step 1: cleared solution formula.
For 0 <= r and r+m <= p-1, set Y_m = ((r+m)!/r!)^3 b_{r+m}
(a polynomial multiple: Y_m = b_{r+m} * prod_{j=1}^m (r+j)^3).
Then Y_{m+1} = P(r+m) Y_m - (r+m)^6 Y_{m-1}, with Y_0 = b_r, Y_1 = (r+1)^3 b_{r+1}.
Hence Y_m = N_m(r) * Y_1 + B_m(r) * Y_0 where B_m ∈ Z[r] is the second
solution: B_0 = 1, B_1 = 0, same recurrence in m.
Derive/verify the closed relation between B_m and a shifted N (standard
three-term recurrence theory gives B_m(r) = -(r+1)^6 N_{m-1}(r+1); VERIFY
this symbolically and prove it by induction — if the exact form differs,
find and prove the correct one).

### Step 2: bordered determinant certificate.
Fix a ≠ 0. Suppose t, t+h, t+k all lie in the fiber F_a (b ≡ a at all
three), with 1 <= h < k, t+k <= p-1. Write Pi_m(t) = prod_{j=1}^m (t+j)^3.
The solution formula gives two linear conditions on the pair (b_{t+1}, a):
  N_h(t)(t+1)^3 b_{t+1} - (Pi_h(t) - B_h(t)) a ≡ 0 (mod p)
  N_k(t)(t+1)^3 b_{t+1} - (Pi_k(t) - B_k(t)) a ≡ 0 (mod p)
Since (b_{t+1}, a) ≠ (0,0) (as a ≠ 0), the determinant vanishes:
  D_{h,k}(t) := N_h(t)(Pi_k(t) - B_k(t)) - N_k(t)(Pi_h(t) - B_h(t)) ≡ 0 (mod p).
D_{h,k} ∈ Z[t] with deg <= 3(h+k) - 3 (compute the exact degree and the
leading coefficient; watch for cancellation in the leading terms — if the
naive leading terms cancel identically, determine the true degree).

### Step 3: nonvanishing of D_{h,k} over F_p.
Prove: for p >= 7 and 1 <= h < k with k <= (p-1)/2 (or the widest range you
can), D_{h,k} is not the zero polynomial mod p. Strategy: evaluate at
t = -1, t = -2 (and if needed the general t = -r) using the endpoint
evaluations above; the values come out as products of factorials (units for
k < p) and Apéry numbers b_j; if all relevant evaluations vanish you get
two consecutive Apéry zeros, contradicting no-consecutive-zeros.
CARE: some evaluations may vanish for structural reasons (forced reflection
factors). If a forced linear factor (e.g. of the shape (2t + const)) divides
D_{h,k} identically, divide it out and prove the quotient is nonzero —
F_p[t] is a domain, so one extra sentence suffices. Document exactly which
structural factors occur.

### Step 4: counting.
Fix a ≠ 0. Write the fiber as t_1 < t_2 < ... < t_r. For each i <= r-2
consider the pattern (h_i, k_i) = (t_{i+1} - t_i, t_{i+2} - t_i).
- Short diameter (k_i <= H): t_i is a root of the nonzero polynomial
  D_{h_i, k_i}, so for each fixed pattern (h,k) there are at most deg D_{h,k}
  <= 3(h+k) such i. Summing over the <= H^2/2 patterns: <= c H^3 indices.
- Long diameter (k_i > H): t_{i+2} - t_i > H, so there are at most 2p/H such i.
Conclude r <= 2p/H + c H^3 + 2; optimize H ~ p^{1/4} to get r = O(p^{3/4}).
Work out the explicit constant.

### Step 5 (corollary, 3 lines): collision energy.
E(p) := sum_{a in F_p} N_p(a)^2 <= max_a N_p(a) * sum_a N_p(a)
      = O(p^{3/4}) * p = O(p^{7/4}).

## Deliverables

1. `problems/3.2/fiber_bound.tex` — a self-contained LaTeX fragment
   (\begin{proposition}...\end{proposition} + full proof + the energy
   corollary), notation matched to proof.tex (N_h, P(n), lem:no-consec,
   lem:nonvanish, rem:content refs as \ref{...}). No \documentclass —
   it will be \input into the paper.
2. `problems/3.2/fiber_verify.py` — pure Python (no sage), verifying:
   (a) symbolically (with fractions or sympy): the solution formula
       Y_m = N_m(r) Y_1 + B_m(r) Y_0 and your closed form for B_m, for
       m <= 8, as polynomial identities;
   (b) the determinant identity: for random primes p in [50, 300], all
       a ≠ 0, every in-fiber triple (t, t+h, t+k) with k <= 25 satisfies
       D_{h,k}(t) ≡ 0 mod p;
   (c) nonvanishing: for those p, no D_{h,k} with h < k <= 25 is the zero
       polynomial mod p (after removing documented structural factors);
   (d) the endpoint evaluation formula for D_{h,k} you use in Step 3;
   (e) empirical max_a N_p(a) for p <= 2000 vs the proved bound curve.
   Script must print PASS/FAIL per check and exit nonzero on any FAIL.
3. If any step CANNOT be closed rigorously, still deliver the working
   parts plus a STALL REPORT at the top of fiber_bound.tex (as a LaTeX
   comment): which step, the precise obstruction, what you tried.

## Hard constraints
- No numerics-as-proof: computations only verify, the LaTeX proof must be
  self-contained mathematics.
- Do not modify proof.tex itself; deliver the two new files only.
- Run fiber_verify.py yourself before finishing; include its output tail
  in your final message.
