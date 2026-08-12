# Q7695 — fast exact scan for Apéry transverse common pairs

## Result

Implemented:

```text
problems/3.2/research/scripts/q32_transverse_common_fast.cpp
problems/3.2/research/scripts/q32_transverse_common_fast.py
```

No manuscript file was edited.

The scanner was compiled against FLINT 3.0.1 and run exhaustively for every prime
`p <= 100000` and every `1 <= r < p`.  The exact result is

```text
primes scanned:       9592
sum of primes:        454396537
Apéry zero positions: 9492
common pairs:         (17,13), (2237,492)
p > 5r pairs:         none
max p/r:              2237/492 = 4.546747967479674796...
```

Independent characteristic-zero recomputation gives

```text
v_17(b_13)       = 1,  v_17(Xi_13)       = 1
v_2237(b_492)    = 1,  v_2237(Xi_492)    = 1
```

so both common factors occur with gcd multiplicity one.

Restricting the new output to `r <= 10000` gives exactly the previously reported
high-pair census `(17,13),(2237,492)`.  This is strong finite exact evidence for
the conjectural barrier, not a proof.

---

## 1. Exact modular recurrences

Let

\[
P(n)=34n^3+51n^2+27n+5,
\qquad D(t)=1-34t+t^2.
\]

For a fixed prime `p`, all denominators in the Apéry recurrence are units for
indices `<p`, so

\[
b_{n+1}\equiv
\bigl(P(n)b_n-n^3b_{n-1}\bigr)(n+1)^{-3}\pmod p
\]

computes `b_0,...,b_{p-1}` in `O(p)` word operations.

For the algebraic factor put

\[
Q(t)=D(t)^{-1/2}=\sum_{n\ge0}q_nt^n.
\]

From `2DQ'+D'Q=0` one gets the exact recurrence

\[
\boxed{(n+1)q_{n+1}=(34n+17)q_n-nq_{n-1}},
\qquad q_0=1,\quad q_1=17.
\]

Thus `Q` is also linear-time modulo `p`.  In fact

\[
q_n=P_n(17)=
\sum_{k=0}^n\binom nk\binom{n+k}{k}8^k\in\mathbb Z,
\]

by the Legendre generating function.

Now

\[
g(t)=F(t)^{-2}Q(t).
\]

Since `F(0)=1`, for any truncation length `L` the exact coefficients modulo
`t^L,p` are obtained by

```text
Finv  = inverse_series(F, L)
Finv2 = mullow(Finv, Finv, L)
g     = mullow(Finv2, Q, L)
```

using FLINT `nmod_poly_inv_series` and `nmod_poly_mullow`.  No coefficient at
index `>=L` can influence this result.

Finally

\[
\Xi_0=-1,
\qquad
\Xi_r-\Xi_{r-1}=-5g_rb_{r-1},
\]

so `Xi` is a single linear pass after the `g` prefix is known.

---

## 2. D-finiteness: deliberately not assumed

The Apéry generating function is D-finite; with
`\theta=t\,d/dt` its Picard--Fuchs operator is

\[
\theta^3-t(2\theta+1)(17\theta^2+17\theta+5)+t^2(\theta+1)^3.
\]

Equivalently,

\[
D\theta^3F+(-51t+3t^2)\theta^2F
+(-27t+3t^2)\theta F+(-5t+t^2)F=0.
\]

But D-finite series are not closed under reciprocals, so this does **not** imply
that `1/F` or `g=1/(F^2 sqrt(D))` is D-finite/P-recursive.  The implementation
uses no such claim.

There is, however, a proved nonlinear differential-algebraic equation.  Put

\[
u=\frac{\theta F}{F}.
\]

Then

\[
\frac{\theta^2F}{F}=\theta u+u^2,
\qquad
\frac{\theta^3F}{F}=\theta^2u+3u\theta u+u^3,
\]

and hence

\[
\boxed{
D(\theta^2u+3u\theta u+u^3)
+(-51t+3t^2)(\theta u+u^2)
+(-27t+3t^2)u+(-5t+t^2)=0.}
\]

The definition of `g` gives

\[
\boxed{
u
=-\frac12\frac{\theta g}{g}
-\frac14\frac{\theta D}{D}.}
\]

Substitution is an explicit nonlinear differential equation for `g`.  It proves
differential algebraicity, not D-finiteness, so it is used only as a structural
check rather than as the scan engine.

---

## 3. Why the scanner is fast

A common pair can occur only at an Apéry zero.  For each prime the C++ worker:

1. computes all `b_0,...,b_{p-1}` modulo `p`;
2. records the zero positions;
3. if there are no zeros, checkpoints immediately and never invokes FLINT;
4. otherwise sets
   \[
   L=1+\max\{r<p:b_r\equiv0\pmod p\};
   \]
5. computes only `Q,g mod t^L`;
6. accumulates `Xi` through `L-1` and tests it only at the stored zero positions.

This is exact because `Xi_r` uses only coefficients through index `r`.

For one active prime, if `M_p(n)` is FLINT's cost for a length-`n`
polynomial product, the work is

\[
O(p)+O(M_p(L)),
\]

up to the Newton-series logarithmic factor.  Globally,

\[
O\!\left(\sum_{p\le P}p+
\sum_{p\le P,\,p\text{ active}}M_p(L_p)\right).
\]

With quasi-linear polynomial multiplication this is roughly quadratic in `P`
up to logarithms.  Prime jobs are independent and therefore parallelize
without communication.

Memory is `O(p)` machine words per worker: inverse table, Apéry row, the
`Q` prefix, five FLINT polynomials, and FLINT multiplication workspace.  At
`p=100000` the explicit coefficient arrays themselves occupy only a few MB per
worker; there is no characteristic-zero multi-GB coefficient table.

---

## 4. Resumable implementation

The C++ worker writes one atomic JSON checkpoint per completed prime:

```text
p_000017.json
p_002237.json
...
```

Each record contains the algorithm version, `p`, `r_limit`, number of Apéry
zeros, actual series truncation, common `r` values, and elapsed time.  It writes
a `.tmp` file first and renames only after a complete flush.  Existing valid
checkpoints are skipped.  A crash therefore loses at most the one prime
currently running on each worker.

The Sage driver balances primes with a `p log p` cost proxy, launches independent
C++ processes, verifies checkpoint versions/completeness, checks the known
`r<=10000` census, computes multiplicities by an independent `Z/p^e Z` power
series path, and writes `summary.json`.

Recommended command:

```bash
sage -python problems/3.2/research/scripts/q32_transverse_common_fast.py all \
  --pmax 100000 \
  --workers 16 \
  --checkpoint-dir problems/3.2/research/output/q32_transverse_fast
```

The multiplicity verifier doubles `e` until both residues are nonzero modulo
`p^e`.  Because every hit has `r<p`, all recurrence denominators remain units
modulo every `p^e`; once the residue is nonzero its `p`-adic valuation is the
true characteristic-zero valuation.

---

## 5. Executed audit

A temporary GitHub Actions workflow (removed after the run) compiled the
committed worker on Ubuntu 24.04 with FLINT 3.0.1 and four independent worker
processes.  Run `31571030901` produced

```text
prime_count 9592
prime_sum 454396537
CHECKPOINTS 9592
TOTAL_B_ZEROS 9492
PAIRS [(17, 13), (2237, 492)]
R_LE_10000 [(17, 13), (2237, 492)]
BARRIER_VIOLATIONS []
MAX_RATIO 2237 492 4.546747967479675
```

The independent characteristic-zero step printed

```text
MULTIPLICITY 17 13 vp_b 1 vp_Xi 1 vp_gcd 1
MULTIPLICITY 2237 492 vp_b 1 vp_Xi 1 vp_gcd 1
```

and checked the initial coefficients

```text
b:              1, 5, 73, ...
D^(-1/2):       1, 17, 433, ...
g:              1, 7, 192, ...
Xi:            -1, -36, -4836, ...
```

The first four-worker full scan took about 44 seconds wall-clock; the sum of
per-prime worker timings was `171.195229509` seconds.  A repeat run had a
worker-time sum of about 203 seconds.  Dependency installation/compilation is
not included in those scan timings.

Thus `p<=100000, r<p` is an actually completed exact computation, not a
projected benchmark.

---

## 6. Scope

The computation proves only the finite statement:

> For every prime `p<=100000` and every `1<=r<p`, simultaneous
> `p|b_r` and `p|Xi_r` occurs exactly at `(17,13)` and `(2237,492)`.

It follows in this finite region that

\[
\max p/r=2237/492<5.
\]

It says nothing about primes above `100000` and must not be used as a proof of
the universal support barrier.

The same code can be extended directly, e.g.

```bash
sage -python problems/3.2/research/scripts/q32_transverse_common_fast.py all \
  --pmax 1000000 --workers 32 --checkpoint-dir /scratch/q32_tc_1m
```

A future barrier-only mode could stop at `r < p/5`, but the implemented default
intentionally scans all `r<p` so it also records every high common pair and the
actual maximum `p/r`.
