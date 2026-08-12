# Q7695 — fast exact scan for Apéry transverse common pairs

## Status

Implemented and executed.

New source files:

```text
problems/3.2/research/scripts/q32_transverse_common_fast.cpp
problems/3.2/research/scripts/q32_transverse_common_fast.py
```

No manuscript file was edited.

The production test scans **every prime `p <= 100000` and every `1 <= r < p`**.  This is the useful triangular region for the proposed support barrier: every counterexample to `p <= 5r` with `p <= 100000` necessarily lies in this scan.

The exact support result is

```text
common pairs with p <= 100000 and 1 <= r < p:
    (17, 13)
    (2237, 492)

p > 5r counterexamples: none
max p/r: 2237/492 = 4.546747967479674796...
```

Independent characteristic-zero recomputation gives

```text
v_17(b_13)       = 1
v_17(Xi_13)      = 1
v_17(gcd)        = 1

v_2237(b_492)    = 1
v_2237(Xi_492)   = 1
v_2237(gcd)      = 1
```

The scan also reproduces the stated characteristic-zero `r <= 10000` high-pair census exactly: restricting its output to `r <= 10000` gives only `(17,13)` and `(2237,492)`.

This is finite exact evidence, **not** a proof of the support barrier.

---

## 1. Exact recurrences used by the scanner

Write

\[
P(n)=34n^3+51n^2+27n+5,
\qquad
D(t)=1-34t+t^2.
\]

The Apéry recurrence is

\[
(n+1)^3 b_{n+1}=P(n)b_n-n^3b_{n-1},
\qquad b_0=1,\quad b_1=5.
\]

For a fixed prime `p`, all denominators `(n+1)^3` are units as long as `n+1<p`.  Therefore `b_0,...,b_{p-1}` are obtained modulo `p` in `O(p)` word operations.  The implementation precomputes all inverses `1,...,p-1` in linear time and uses

```text
b[n+1] = (P(n)b[n] - n^3 b[n-1]) * (n+1)^(-3) mod p.
```

### 1.1 The algebraic factor `D^{-1/2}` has a two-term recurrence

Put

\[
Q(t)=D(t)^{-1/2}=\sum_{n\ge0}q_nt^n.
\]

Differentiating `Q=D^{-1/2}` gives

\[
2D Q'+D'Q=0.
\]

Equating the coefficient of `t^n` gives

\[
\boxed{(n+1)q_{n+1}=(34n+17)q_n-nq_{n-1}},
\]

with

\[
q_0=1,\qquad q_1=17.
\]

Thus `q_0,...,q_{p-1}` are also computable modulo `p` in `O(p)` word operations.

There is also an integrality interpretation.  The Legendre generating function gives

\[
Q(t)=\sum_{n\ge0}P_n(17)t^n,
\]

and, since `17=1+2*8`,

\[
P_n(17)=
\sum_{k=0}^n
\binom nk\binom{n+k}{k}8^k\in\mathbb Z.
\]

So `Q` is actually integral, although the modular scan only needs the weaker fact that every denominator up to degree `<p` is a `p`-unit.

### 1.2 `g` is computed by fast unit-series inversion, not by an unproved recurrence

Since

\[
g(t)=\frac{1}{F(t)^2\sqrt{D(t)}}
     =F(t)^{-2}Q(t),
\]

and `F(0)=1`, the truncated series is obtained exactly modulo `p` by

```text
Finv  = inverse_series(F, L)
Finv2 = mullow(Finv, Finv, L)
g     = mullow(Finv2, Q, L)
```

using FLINT `nmod_poly_inv_series` and `nmod_poly_mullow`.

Only coefficients below `L` can affect the result below `L`; no omitted coefficient of `F` or `Q` enters any decision.

Finally

\[
\Xi_0=-1,
\qquad
\Xi_r-\Xi_{r-1}=-5g_rb_{r-1},
\]

so after `g_0,...,g_{L-1}` are known, the `Xi` prefix is a single linear pass.

---

## 2. Why I do not use a P-recursive recurrence for `g`

The Apéry generating function is D-finite.  In logarithmic notation

\[
\theta=t\frac{d}{dt},
\]

its Picard--Fuchs operator is

\[
\boxed{
\theta^3
-t(2\theta+1)(17\theta^2+17\theta+5)
+t^2(\theta+1)^3.}
\]

Equivalently,

\[
D\theta^3F
+(-51t+3t^2)\theta^2F
+(-27t+3t^2)\theta F
+(-5t+t^2)F=0.
\]

However, **D-finite functions are not closed under reciprocals**.  Hence neither `1/F` nor

\[
g=1/(F^2\sqrt D)
\]

may be declared D-finite merely from this equation.

The present implementation makes no such assumption.

For completeness, `g` does satisfy a proved nonlinear differential-algebraic equation.  Put

\[
u=\frac{\theta F}{F}.
\]

Then

\[
\frac{\theta^2F}{F}=\theta u+u^2,
\qquad
\frac{\theta^3F}{F}=\theta^2u+3u\theta u+u^3,
\]

so division of the Picard--Fuchs equation by `F` gives

\[
\boxed{
D(\theta^2u+3u\theta u+u^3)
+(-51t+3t^2)(\theta u+u^2)
+(-27t+3t^2)u
+(-5t+t^2)=0.}
\]

From the definition of `g`,

\[
\frac{\theta g}{g}
=-2u-\frac12\frac{\theta D}{D},
\]

hence

\[
\boxed{
u
=-\frac12\frac{\theta g}{g}
-\frac14\frac{\theta D}{D}.}
\]

Substitution gives an explicit nonlinear differential equation for `g`.  This proves differential algebraicity, but it does **not** prove D-finiteness and it does not yield a linear P-recurrence.  The scanner therefore uses fast power-series arithmetic instead.

---

## 3. The important practical optimization: invert only at Apéry-zero primes, and only to the last zero

For a common pair we first need

\[
b_r\equiv0\pmod p.
\]

The C++ worker therefore does this for each prime:

1. compute the full Apéry row `b_0,...,b_{p-1}` modulo `p`;
2. collect the zero positions;
3. if there are no zeros, checkpoint the prime immediately and do **no** polynomial inversion;
4. otherwise let
   \[
   L=1+\max\{r<p:b_r\equiv0\pmod p\};
   \]
5. compute only `Q mod t^L` and `g mod t^L`;
6. update `Xi` only through `L-1` and test `Xi_r` only at the stored Apéry-zero positions.

This is exact because `Xi_r` depends only on coefficients with index at most `r`.

In the completed `p<=100000` run there were

```text
9592 primes
454396537 = sum_{p<=100000} p
9492 total Apéry-zero positions b_r=0 with 1<=r<p
```

but only two positions also had `Xi_r=0`.

---

## 4. Implementation

### 4.1 C++/FLINT worker

File:

```text
problems/3.2/research/scripts/q32_transverse_common_fast.cpp
```

It uses only word-size modular arithmetic plus FLINT `nmod_poly` for the series operations.  The worker can scan a numeric prime range or an explicit prime list.

Example build:

```bash
g++ -O3 -DNDEBUG -std=c++17 \
  problems/3.2/research/scripts/q32_transverse_common_fast.cpp \
  -lflint -lgmp -lmpfr -o /tmp/q32_transverse_common_fast
```

Example direct run:

```bash
/tmp/q32_transverse_common_fast \
  --pmax 100000 \
  --checkpoint-dir /tmp/q32_transverse_fast
```

### 4.2 Sage orchestration and independent verifier

File:

```text
problems/3.2/research/scripts/q32_transverse_common_fast.py
```

Recommended run:

```bash
sage -python problems/3.2/research/scripts/q32_transverse_common_fast.py all \
  --pmax 100000 \
  --workers 16 \
  --checkpoint-dir problems/3.2/research/output/q32_transverse_fast
```

The driver:

- compiles the C++ worker against Sage's FLINT installation;
- generates all primes exactly with Sage;
- uses a greedy `p log p` workload partition;
- launches independent worker processes;
- resumes from atomic per-prime JSON checkpoints;
- checks that the `r<=10000` high-pair projection is exactly
  `[(17,13),(2237,492)]`;
- computes exact multiplicities for every hit by a separate Sage path modulo increasing prime powers;
- writes a deterministic `summary.json`.

The multiplicity verifier repeatedly works modulo `p^e`, doubling `e` until both `b_r` and `Xi_r` are nonzero.  Because `r<p`, every recurrence denominator remains a unit modulo every `p^e`.  Once a residue is nonzero modulo `p^e`, with valuation `<e`, its `p`-adic valuation is the true characteristic-zero valuation.

---

## 5. Checkpoint/restart semantics

Each completed prime has its own file

```text
p_000017.json
p_002237.json
...
```

with fields

```json
{
  "version": "q32-transverse-common-fast-v1",
  "p": 2237,
  "r_limit": 2236,
  "b_zero_count": 0,
  "truncation": 0,
  "pairs": [],
  "seconds": 0.0
}
```

(the displayed numerical fields are schematic except for `p`; the actual file records the measured values).

The worker writes `p_XXXXXX.json.tmp` and renames it only after the complete prime result has been flushed.  On restart:

- an existing valid checkpoint skips that prime;
- a leftover `.tmp` file is ignored;
- a missing prime is simply reassigned;
- workers never share a prime.

Thus interruption loses at most the currently active prime on each worker.  There is no monolithic in-memory state to reconstruct.

The Sage summarizer also checks the algorithm-version field before trusting a checkpoint.

---

## 6. Exact run and regression check

A temporary GitHub Actions audit was used to compile the committed C++ worker against FLINT 3.0.1 and execute the full range.  The temporary workflow was removed after the audit; the run record remains in GitHub Actions.

Audit run:

```text
run id: 31571030901
platform: ubuntu-24.04
workers: 4 independent C++ processes
FLINT: 3.0.1
```

Exact support output:

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

The second audit also performed an independent characteristic-zero calculation of `F^{-1}`, `D^{-1/2}`, `g`, and `Xi` through the two hit indices.  It returned

```text
MULTIPLICITY 17   13  vp_b 1 vp_Xi 1 vp_gcd 1
MULTIPLICITY 2237 492 vp_b 1 vp_Xi 1 vp_gcd 1
```

The characteristic-zero path starts with

```text
b = [1, 5, 73, ...]
D^(-1/2) = [1, 17, 433, ...]
g = [1, 7, 192, ...]
Xi = [-1, -36, -4836, ...]
```

and therefore also serves as a small independent formula check.

### Runtime

On the first four-worker hosted run, the exact scan step took about 44 seconds wall-clock; the sum of all per-prime worker timings was

```text
171.195229509 seconds.
```

A repeat run gave a worker-time sum of about 203 seconds.  Package installation and compilation were outside those scan timings.

This establishes that `p<=100000, r<p` is not merely an aspirational target; the committed core has actually completed it on an ordinary hosted four-core runner.

---

## 7. Complexity

Let `P` be the prime ceiling and, for an active prime `p`, let

\[
L_p=1+\max\{r<p:b_r\equiv0\pmod p\}.
\]

Let `M_p(n)` denote the cost of multiplying length-`n` polynomials over `F_p` using FLINT's selected algorithm.

For one prime:

```text
Apéry recurrence + inverse table: O(p)
zero scan:                         O(p)
Q=D^(-1/2) recurrence:             O(L_p)
F inverse + two low products:      O(M_p(L_p)) up to Newton/log factors
Xi accumulation:                   O(L_p)
```

and no series work is done when `b` has no zero.

Thus the total cost is

\[
O\!\left(
\sum_{p\le P}p
+
\sum_{\substack{p\le P\\p\;\mathrm{active}}}
M_p(L_p)
\right).
\]

With quasi-linear polynomial multiplication `M_p(n)=\widetilde O(n)`, the crude worst case is approximately quadratic in `P` up to logarithmic factors.  More importantly, the prime jobs are completely independent, so wall time scales almost ideally until memory bandwidth or the process count becomes limiting.

For `P=100000`, the unavoidable linear Apéry work is exactly 454,396,537 coefficient steps across all primes.

### Memory

For one worker the persistent arrays are linear in the current prime/truncation:

- modular inverses: `p` machine words;
- Apéry row: `p` machine words;
- `Q`: `L_p` words;
- five FLINT polynomials of length at most `L_p`, plus FLINT multiplication temporaries.

Hence memory is `O(p)` words per worker, not `O(P^2)` and not a characteristic-zero multi-gigabyte coefficient table.  At `p=100000` the explicit C++ vectors/polynomial coefficient arrays themselves are only a few megabytes per worker; FLINT temporary space is the larger constant-factor component but remains linear/quasi-linear workspace for these lengths.

---

## 8. Why this scan is preferable to a characteristic-zero `r=100000` factor scan

The previous characteristic-zero approach builds rapidly growing integers and then factors or partially factors `gcd(b_r,Xi_r)`.  At `r=100000`, both coefficient height and repeated large-integer gcd/factor work become the dominant cost.

For the support-barrier question the prime triangle is better targeted:

```text
for every p <= 100000:
    inspect every r < p exactly modulo p.
```

Any violation `p>5r` in this prime range is automatically inspected, and no factorization is needed.  Each decision is only the exact pair of congruences

\[
b_r\equiv0\pmod p,
\qquad
\Xi_r\equiv0\pmod p.
\]

The resulting test is therefore substantially deeper in `r` for large `p` than the old `r<=10000` census while remaining cheap enough to rerun routinely.

---

## 9. What the computation does and does not establish

It establishes the finite statement

> For every prime `p<=100000` and every integer `1<=r<p`, simultaneous
> `p|b_r` and `p|Xi_r` occurs exactly at `(17,13)` and `(2237,492)`.

Both common factors occur with gcd multiplicity one.

Consequently, in this finite region,

\[
\max\frac{p}{r}=\frac{2237}{492}<5
\]

and there is no counterexample to `p<=5r`.

It does **not** establish the universal support barrier.  In particular it says nothing about primes `p>100000`, and the numerical absence of further pairs must not be fed into a proof as an arithmetic theorem.

---

## 10. Suggested next computational extension

The code already scales by prime rather than by characteristic-zero height.  The clean next tests are:

```bash
# same triangle, larger prime ceiling
sage -python problems/3.2/research/scripts/q32_transverse_common_fast.py all \
  --pmax 1000000 --workers 32 --checkpoint-dir /scratch/q32_tc_1m
```

or, if only the barrier is wanted, a future specialized mode may stop each prime at `floor((p-1)/5)` and test only `p>5r`.  That would substantially reduce the power-series length but would no longer record all `p>r` common pairs, so it is intentionally not the default scanner implemented here.
